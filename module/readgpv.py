# -*- coding: utf-8 -*-
import numpy as np
import psutil
import setup
import scipy.linalg
from scipy import integrate
import sys
import subprocess

class ReadGPV:
  def __init__(self,dataset,date,ft):
    """使用する各種の設定"""
    ST = setup.Setup(dataset)
    self.nx, self.ny, self.nz, self.mem = ST.set_prm()
    self.surf = 1
    self.press_levels = ST.set_pressure_levels()
    self.date, self.init, self.ft = date,'00',ft
    self.surf_elem, self.elem = self.data_kind()
    self.elem_num = len(self.elem)

  def data_kind(self):
    surf_elem = ( 'UGRD', 'VGRD', 'PRMSL', 'APCP')
    elem      = ( 'UGRD', 'VGRD', 'HGT', 'TMP') 
    return surf_elem, elem

  def init_array(self):
    uwnd_data = np.zeros((self.mem, self.nz, self.ny, self.nx))
    vwnd_data = np.zeros((self.mem, self.nz, self.ny, self.nx))
    hgt_data  = np.zeros((self.mem, self.nz-self.surf, self.ny, self.nx))
    tmp_data  = np.zeros((self.mem, self.nz-self.surf, self.ny, self.nx))
    slp_data  = np.zeros((self.mem, self.surf, self.ny, self.nx))
    rain_data = np.zeros((self.mem, self.surf, self.ny, self.nx))
    return uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data

  def data_read_init_driver(self, data_path:str):
    """
    Args:
      data_path(str): 週間アンサンブルデータのPATH
    Returns:
      data structure: (アンサンブル数, 鉛直層, 緯度, 経度)
      dara (np.ndarray): 各種要素のデータ
    """
    uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = self.init_array()

    for imem in range(self.mem):
      full_mem_data = self.read_gpv(data_path+'/{:03}/{}_{}hr.grd'.format(imem+1,self.date,self.init),self.elem_num)
      uwnd_data[imem,:,:,:] = full_mem_data[0,:,:,:]
      vwnd_data[imem,:,:,:] = full_mem_data[1,:,:,:]
      hgt_data[imem,:,:,:]  = full_mem_data[2,1:,:,:]
      slp_data[imem,:,:,:]  = full_mem_data[2,0,:,:]*0.01
      tmp_data[imem,:,:,:]  = full_mem_data[3,1:,:,:]
      rain_data[imem,:,:,:] = full_mem_data[3,0,:,:]

    return uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data

  def data_read_ft_driver(self, data_path:str):
    """
    Args:
      data_path(str): 週間アンサンブルデータのPATH
    Returns:
      data structure: (アンサンブル数, 鉛直層, 緯度, 経度)
      dara (np.ndarray): 各種要素のデータ
    """
    uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = self.init_array()
    
    for imem in range(self.mem):
      full_mem_data = self.read_gpv(data_path+'/{:03}/{}_{}hr.grd'.format(imem+1,self.date,self.ft),self.elem_num)
      uwnd_data[imem,:,:,:] = full_mem_data[0,:,:,:]
      vwnd_data[imem,:,:,:] = full_mem_data[1,:,:,:]
      hgt_data[imem,:,:,:]  = full_mem_data[2,1:,:,:]
      slp_data[imem,:,:,:]  = full_mem_data[2,0,:,:]*0.01
      tmp_data[imem,:,:,:]  = full_mem_data[3,1:,:,:]
      rain_data[imem,:,:,:] = full_mem_data[3,0,:,:]

    return uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data

  def set_gpv(self, gpv_file, elem):
    return self._open_gpv(gpv_file).reshape(elem, self.nz, self.mem, self.ny, self.nx)

  def read_gpv(self, gpv_file, elem):
    return self._open_gpv(gpv_file).reshape(elem, self.nz, self.ny, self.nx)

  def _open_gpv(self, gpv_file, *, mode='default'):
    print('...... Preparating for {}'.format(gpv_file))
    with open(gpv_file, 'rb') as ifile:
      if mode == 'default':
        data = np.fromfile(ifile, dtype='<f', sep = '')
    return data

  def set_coordinate(self, dx=2.5, dy=2.5):
    lon, lat = [], []
    for ix in range(self.nx):
      lon += [ float('{:.2f}'.format(dx*ix)) ]
    for iy in range(self.ny):
      lat += [ float('{:.2f}'.format(90.0-dy*iy)) ]
    X, Y = np.meshgrid(lon, lat)
    return X, Y

  def weight_latitude(self, lat:np.ndarray) -> np.ndarray:
    return np.sqrt(np.cos(np.deg2rad(lat)))

class Energy_NORM:
  def __init__(self,dataset):
    ST = setup.Setup(dataset)
    self.nx, self.ny, self.nz, self.mem = ST.set_prm()
    self.press_levels = ST.set_pressure_levels()
    self.ctrl, self.surf=1,1
    self.Pr:float=1000.0
    self.Tr:float=270.0
    self.cp:float=1004.0
    self.R:float=287.0

  def init_array(self):
    pertb_uwnd_data = np.zeros((self.mem-self.ctrl, self.nz, self.ny, self.nx))
    pertb_vwnd_data = np.zeros((self.mem-self.ctrl, self.nz, self.ny, self.nx))
    pertb_tmp_data  = np.zeros((self.mem-self.ctrl, self.nz-self.surf, self.ny, self.nx))
    pertb_slp_data  = np.zeros((self.mem-self.ctrl, self.surf, self.ny, self.nx))
    return pertb_uwnd_data, pertb_vwnd_data, pertb_tmp_data, pertb_slp_data

  def data_pertb_driver(self,uwnd,vwnd,tmp,slp):
    """コントロールランからの摂動の作成
    Args:
      ctrl_run(np.ndarray): 摂動を与えていないコントロールランのデータ
      ensm_run(np.ndarray): 各アンサンブルランのデータ
    Returns:
      pretb_elem_data(np.ndarray): アンサンブルランのデータからコントロールランデータを引いた擾乱のデータ
    """
    pertb_uwnd_data, pertb_vwnd_data, pertb_tmp_data, pertb_slp_data = self.init_array()
    for imem in range(self.mem-self.ctrl):
      pertb_uwnd_data[imem,:,:,:] = uwnd[imem+1,:,:,:] - uwnd[self.ctrl-1,:,:,:]
      pertb_vwnd_data[imem,:,:,:] = vwnd[imem+1,:,:,:] - vwnd[self.ctrl-1,:,:,:]
      pertb_tmp_data[imem,:,:,:] = tmp[imem+1,:,:,:] - tmp[self.ctrl-1,:,:,:]
      pertb_slp_data[imem,:,:,:] = slp[imem+1,:,:,:] - slp[self.ctrl-1,:,:,:]

    return pertb_uwnd_data, pertb_vwnd_data, pertb_tmp_data, pertb_slp_data

  def weight_average(self, data:np.ndarray, weight_list:np.ndarray):
    weight_average, sum_of_weight = np.average( 
      a = data, axis = 0, weights = weight_list, returned = True
    )
    #print('..... CALCULATE WEIGHT AVE. SUM OF WEIGHT ', sum_of_weight)
    return weight_average
    
  def calc_dry_EN_NORM(self,
    u_prime:np.ndarray, v_prime:np.ndarray, tmp_prime:np.ndarray, slp_prime:np.ndarray 
  ):
    """乾燥エネルギーノルムの計算
    Args:
      u_prime   (np.ndarray): 東西風のコントロールランからの予測時間における摂動
      v_prime   (np.ndarray): 南北風のコントロールランからの予測時間における摂動
      tmp_prime (np.ndarray): 気温のコントロールランからの予測時間における摂動
      slp_prime (np.ndarray): 海面更生気圧のコントロールランからの予測時間における摂動
    Parameters:
      self.Pr (float) : 経験的に求めた参照気圧. Defaults to 1000 hPa.
      self.Tr (float) : 経験的に求めた参照気温. Defaults to 270 K.
      self.cp (float) : 定圧比熱. Defaults to 1004 J/K*kg.
      self.R  (float) : 気体の状態定数. Defaults to 287.0 J/K*kg.
    Returns:
      dry_energy_norm (np.ndarray): トータル乾燥エネルギーノルム(J/kg)
      constitution -> [緯度, 経度] 
    """

    #Physics
    physical_term = (u_prime)**2+(v_prime)**2
    vint_physical_term = self._vint(physical_term,self.press_levels)/(2*self.Pr)
    
    #Potential
    tmp_term = (self.cp/self.Tr)*((tmp_prime)**2)
    vint_tmp_term = self._vint(tmp_term,self.press_levels[1:])/(2*self.Pr)
    slp_term = (self.R*self.Tr/self.Pr)*(slp_prime**2/self.Pr)*0.5
    vint_potential_term = vint_tmp_term + slp_term

    #SUM OF TERM
    dry_energy_norm = vint_physical_term + vint_potential_term

    return dry_energy_norm, vint_physical_term, vint_potential_term 

  def calc_humid_EN_NORM(self):
    pass

  def verification_region(self, 
    lon, lat, *,
    area_lat_min:float =50.0, area_lat_max:float =20.0,
    area_lon_min:float =120.0, area_lon_max:float =150.0
  ):
    """検証領域のインデックス番号を返す
    Args:
      lon, lat (np.ndarray): 経度, 緯度のnp.ndarray
      area_lat_min (float, optional): 検証領域の緯度最低値. Defaults to 50.
      area_lat_max (float, optional): 検証領域の緯度最高値. Defaults to 20.
      area_lon_min (float, optional): 検証領域の経度最低値. Defaults to 120.
      area_lon_max (float, optional): 検証領域の経度最高値. Defaults to 150.
    Retunrs:
      各値に相当する　listのindex番号.
    Note:
      lambert図法の気象研のoutputでは, おそらく使えないだろう...
    """
    area_lat_min_index, area_lat_max_index = np.where(lat == area_lat_min)[0][0], np.where(lat == area_lat_max)[0][0]
    area_lon_min_index, area_lon_max_index = np.where(lon == area_lon_min)[1][0], np.where(lon == area_lon_max)[1][0]

    return area_lat_min_index, area_lat_max_index, \
           area_lon_min_index, area_lon_max_index

  def _vint(self, x_array:np.ndarray, press_array:np.ndarray) -> np.ndarray:
    """気圧面鉛直積分
    Args:
      x_array (np.ndarray)     : 持っている要素の離散値の値
      press_array (np.ndarray) : 今回使用するp面座標の差の値
    Returns:
      y_array (np.ndarray)     : 鉛直積分した値
    Note:
      今回使用した積分方法はシンプソン則
    """
    y_array = np.empty_like(x_array[0], dtype=np.float32)

    for iy in range(self.ny):
      for ix in range(self.nx):
        y_array[iy,ix] = integrate.simps(x_array[:,iy,ix], press_array[::-1])

    return y_array

  def singular_decomposion(self, array, *, try_num=5):
    """特異値分解
    Args:
      array (np.ndarray)   : 特異値分解したい行列(n,m) 
    Return:
      U_array (np.ndarray)     : ユニタリ行列(n,n)
      sigma_array (np.ndarray) : 特異値の対角成分(r,r) -> sigma_array*sigma_array/m で(array, array.T)の固有値
      V_array (np.ndarray)     : アジョイント行列(n,m)
    """
   
    for _ in range(try_num):
      try:
        #U_array, sigma_array, V_array = numpy.linalg.svd(array,full_matrices=True)
        U_array, sigma_array, V_array = scipy.linalg.svd(array)
        print('..... SUCCESS SINGULAR VECTOR CALCULATION ')
        return U_array, sigma_array, V_array

      except Exception as e:
        print('   >>> ONE MORE TRY <<<   ')
        command = ["sleep", "5.0s"]
        res = subprocess.call(command)
      
      else:
        break
    
    else:
      print('..... NOT SUCCESS SINGULAR VECTOR CALCULATION ')
      mem = psutil.virtual_memory() 
      total, used, available = mem.total, mem.used, mem.available
      percent_total, percent_used, percent_available = (total/total)*100, (used/total)*100, (available/total)*100

      print(
        '..... CHECK MEMORY PERCENT MAX:{:.2f}, USED:{:.2f}, EMPTY:{:.2f}'.format(percent_total,percent_used,percent_available)
        )
      print('')
      sys.exit()



