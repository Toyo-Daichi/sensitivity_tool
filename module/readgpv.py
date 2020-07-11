# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate

class ReadGPV:
  def __init__(self, nx, ny, nz, mem):
    self.nx = nx
    self.ny = ny
    self.nz = nz
    self.surf = 1
    self.ensemble_size = mem

  def data_kind(self, mode='default'):
    if mode is 'default':
      surf_elem = { 'UGRD': 0, 'VGRD': 1, 'PRMSL':2, 'APCP':3 }
      elem      = { 'UGRD': 0, 'VGRD': 1, 'HGT':2, 'TMP':3 } 
    return surf_elem, elem

  def set_coordinate(self, dx=2.5, dy=2.5):
    lon, lat = [], []
    for ix in range(self.nx):
      lon += [ float('{:.2f}'.format(dx*ix)) ]
    for iy in range(self.ny):
      lat += [ float('{:.2f}'.format(90.0-dy*iy)) ]
    X, Y = np.meshgrid(lon, lat)
    return X, Y
    
  def read_gpv(self, gpv_file, elem):
    return self._open_gpv(gpv_file).reshape(elem, self.nz, self.ensemble_size, self.ny, self.nx)

  def _open_gpv(self, gpv_file, *, mode='default'):
    print('...... Preparating for {}'.format(gpv_file))
    with open(gpv_file, 'rb') as ifile:
      if mode == 'default':
        data = np.fromfile(ifile, dtype='<f', sep = '')
    return data

  def calc_prime(self, ctrl_run:np.ndarray, ensm_run:np.ndarray) -> np.ndarray:
    """コントロールランからの摂動の作成
    Args:
      ctrl_run (np.ndarray): 摂動を与えていないコントロールランのデータ
      ensm_run (np.ndarray): 各アンサンブルランのデータ
    Returns:
      (np.ndarray): アンサンブルランのデータからコントロールランデータを引いた擾乱のデータ
    """
    return ensm_run - ctrl_run

  def weight_latitude(self, lat:np.ndarray) -> np.ndarray:
    return np.sqrt(np.cos(np.deg2rad(lat)))

  def weight_average(self, data:np.ndarray, weight_list:np.ndarray):
    weight_average, sum_of_weight = np.average( 
      a = data,
      axis = 0,
      weights = weight_list,
      returned = True
    )
    #print('..... CALCULATE WEIGHT AVE. SUM OF WEIGHT ', sum_of_weight)
    return weight_average

class Energy_norm:
  def __init__(self, nx, ny):
    self.nx = nx
    self.ny = ny
    self.Pr:float=1000.0
    self.Tr:float=270.0
    self.cp:float=1004.0
    self.R:float=287.0 

  def dry_energy_norm(self,
    u_prime:np.ndarray, v_prime:np.ndarray, tmp_prime:np.ndarray, slp_prime:np.ndarray, press_levels:np.ndarray
  ):
    """乾燥エネルギーノルムの計算
    Args:
      u_prime   (np.ndarray): 東西風のコントロールランからの予測時間における摂動
      v_prime   (np.ndarray): 南北風のコントロールランからの予測時間における摂動
      tmp_prime (np.ndarray): 気温のコントロールランからの予測時間における摂動
      slp_prime (np.ndarray): 海面更生気圧のコントロールランからの予測時間における摂動
      press_levels (np.ndarray) : 鉛直積分で用いる気圧面のリスト
    Parameters:
      self.Pr (float) : 経験的に求めた参照気圧. Defaults to 1000 hPa.
      self.Tr (float) : 経験的に求めた参照気温. Defaults to 270 K.
      self.cp (float) : 定圧比熱. Defaults to 1004 J/K*kg.
      self.R  (float) : 気体の状態定数. Defaults to 287.0 J/K*kg.
    Returns:
      dry_energy_norm (np.ndarray): トータル乾燥エネルギーノルム(J/kg)
      constitution -> [緯度, 経度] 
    """

    physical_term = (u_prime)**2 + (v_prime)**2
    potential_term = (self.cp/self.Tr)*((tmp_prime)**2)
    
    first_term = self._vint(physical_term+potential_term, press_levels)
    first_term = first_term/(2*self.Pr)

    sec_term = (self.Tr*self.R/self.Pr**2)*(slp_prime**2)*0.5

    dry_energy_norm = first_term+sec_term 

    return dry_energy_norm

  def humid_energy_norm(self):
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
