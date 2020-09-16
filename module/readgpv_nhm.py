# -*- coding: utf-8 -*-
import numpy as np
import setup

class ReadGPV:
  def __init__(self, dataset, *, read_num:int=0):
    ST = setup.Setup(dataset)
    self.nx, self.ny, self.nz, self.mem = ST.set_prm()
    self.surf = 1
    self.dims = self.nx*self.ny
    self.dataset = dataset
    if 'p_level' in dataset:
      self.elem_list = ( 'UGRD', 'VGRD', 'TMP ', 'SPFH', 'HGT ', 'WWND', 'VOR ', 'SLP ', 'PS  ', 'RAIN')
    elif 'm_level' in dataset:
      self.elem_list = ( 'UGRD', 'VGRD', 'QV  ', 'QC  ', 'QR  ', 'PRS ')

    if not read_num:
      print('..... @ CHECK ELEM LIST @ ')
      print(self.elem_list)

  def init_array(self):
    uwnd_data = np.zeros((self.nz+self.surf, self.ny, self.nx))
    vwnd_data = np.zeros((self.nz+self.surf, self.ny, self.nx))
    tmp_data  = np.zeros((self.nz+self.surf, self.ny, self.nx))
    spfh_data = np.zeros((self.nz, self.ny, self.nx))
    hgt_data  = np.zeros((self.nz, self.ny, self.nx))
    wwnd_data = np.zeros((self.nz, self.ny, self.nx))
    vor_data  = np.zeros((      3, self.ny, self.nx))
    slp_data  = np.zeros((self.ny, self.nx))
    ps_data   = np.zeros((self.ny, self.nx))
    rain_data = np.zeros((self.ny, self.nx))

    return uwnd_data, vwnd_data, tmp_data, spfh_data, hgt_data, wwnd_data, vor_data, slp_data, ps_data, rain_data

  def init_mem_array(self):
    uwnd_data = np.zeros((self.mem, self.nz+self.surf, self.ny, self.nx))
    vwnd_data = np.zeros((self.mem, self.nz+self.surf, self.ny, self.nx))
    tmp_data  = np.zeros((self.mem, self.nz+self.surf, self.ny, self.nx))
    spfh_data = np.zeros((self.mem, self.nz, self.ny, self.nx))
    hgt_data  = np.zeros((self.mem, self.nz, self.ny, self.nx))
    wwnd_data = np.zeros((self.mem, self.nz, self.ny, self.nx))
    vor_data  = np.zeros((self.mem,       3, self.ny, self.nx))
    slp_data  = np.zeros((self.mem, self.ny, self.nx))
    ps_data   = np.zeros((self.mem, self.ny, self.nx))
    rain_data = np.zeros((self.mem, self.ny, self.nx))

    return uwnd_data, vwnd_data, tmp_data, spfh_data, hgt_data, wwnd_data, vor_data, slp_data, ps_data, rain_data 

  def data_read_driver(self, exp_path:str, target_date, *, endian='little'):
    uwnd_data, vwnd_data, tmp_data, spfh_data, _, _, _, _, ps_data, _ = self.init_memarray

    for index, imem in enumerate(self.mem):
      data_path = exp_path + '/{:03}/plevel_grd/{}00.grd'.format(imem+1, target_date)
      _full_data = self._open_gpv(data_path, endian=endian)
      _uwnd_data, _vwnd_data, _tmp_data, _spfh_data, _, _, _, _, _ps_data, _ = self.full_data_cut(_full_data)
      uwnd_data[index,:,:,:] = _uwnd_data[:,:,:]
      vwnd_data[index,:,:,:] = _vwnd_data[:,:,:]
      tmp_data[index,:,:,:]  =  _tmp_data[:,:,:]
      spfh_data[index,:,:,:] = _spfh_data[:,:,:]
      ps_data[index,:,:]     =   _ps_data[:,:]

    return uwnd_data, vwnd_data, tmp_data, spfh_data, ps_data

  def full_data_cut(self, full_data:np.ndarray):
    uwnd_data, vwnd_data, tmp_data, spfh_data, hgt_data, wwnd_data, vor_data, slp_data, ps_data, rain_data = self.init_array()

    # cut
    # (1) full dataset
    uwnd_data = full_data[0:(self.surf+self.nz)*self.dims].reshape(self.surf+self.nz,self.ny,self.nx)
    vwnd_data = full_data[1*(self.surf+self.nz)*self.dims:2*((self.surf+self.nz)*self.dims)].reshape(self.surf+self.nz,self.ny,self.nx)
    tmp_data  = full_data[2*(self.surf+self.nz)*self.dims:3*((self.surf+self.nz)*self.dims)].reshape(self.surf+self.nz,self.ny,self.nx)
    split_full_dims = 3*((self.surf+self.nz)*self.dims)
    # (1) 3d dataset (not include surface)
    spfh_data = full_data[split_full_dims:(split_full_dims+(self.nz)*self.dims)].reshape(self.nz,self.ny,self.nx)
    hgt_data  = full_data[split_full_dims+(1*(self.nz)*self.dims):split_full_dims+(2*(self.nz)*self.dims)].reshape(self.nz,self.ny,self.nx)
    wwnd_data = full_data[split_full_dims+(2*(self.nz)*self.dims):split_full_dims+(3*(self.nz)*self.dims)].reshape(self.nz,self.ny,self.nx)
    split_full_dims = split_full_dims+3*(self.nz)*self.dims
    # (2) 3d dataset (3 layer)
    vor_data  = full_data[split_full_dims:(split_full_dims+3*self.dims)].reshape(3,self.ny,self.nx)
    split_full_dims = split_full_dims+3*self.dims 
    # (3) 2ddataset
    slp_data  = full_data[split_full_dims:(split_full_dims+1*self.dims)].reshape(self.ny,self.nx)
    ps_data   = full_data[split_full_dims+1*self.dims:(split_full_dims+2*self.dims)].reshape(self.ny,self.nx)
    rain_data = full_data[split_full_dims+2*self.dims:(split_full_dims+3*self.dims)].reshape(self.ny,self.nx)
    
    return uwnd_data, vwnd_data, tmp_data, spfh_data, hgt_data, wwnd_data, vor_data, slp_data, ps_data, rain_data 

  def set_coordinate(self, nx:int, ny:int, *, inv:str='off') -> np.ndarray:
    """データの座標系取得
    Args:
        nx (int): x座標の数
        ny (int): y座標の数
    Returns:
        np.ndarray: x座標, y座標
    """
    path = '/work3/daichi/Data'
    grd_dir = path + '/Coordinate/Nhm_grd_{:}-{:}/'.format(nx,ny)
    ifile_lat = grd_dir + 'lat_lambert.grd'
    ifile_lon = grd_dir + 'lon_lambert.grd'

    print('...... Preparating coord. lat, lon data :: {}, {}'.format(ny,nx))
    
    with open(ifile_lon, 'rb') as ifile:
      X = np.fromfile(ifile, dtype='float64', sep = '').reshape(ny, nx)
    
    with open(ifile_lat, 'rb') as ifile:
      if ( inv is 'off' ): 
        Y = np.fromfile(ifile, dtype='float64', sep = '').reshape(ny, nx)
      
      elif ( inv is 'on' ):
        _Y = np.fromfile(ifile, dtype='float64', sep = '').reshape(ny, nx)
        Y = np.flip(_Y, axis=0) 	
    
    return X, Y

  def set_coordinate_txt(self, path:str, nx:int, ny:int, *, inv:str='off') -> np.ndarray:
    """データの座標系取得
    Args:
        path (str): txt形式の緯度経度データ場所のPATH
        nx (int): x座標の数
        ny (int): y座標の数
    Returns:
        np.ndarray: x座標, y座標
    """
    grd_dir = path + '/Coordinate/Nhm_grd_' + str(self.nx) + '-' + str(self.ny) + '/'
    ifile_lat = grd_dir + 'lat_lambertgrd.txt'
    ifile_lon = grd_dir + 'lon_lambertgrd.txt'
    
    print('...... Preparating coord. lat, lon data :: ' + str(self.ny) + ', ' + str(self.nx))

    with open(ifile_lat, 'r') as ifile:
      if ( inv is 'off' ):
        lat = [s.strip() for s in ifile.readlines()]
      elif ( inv is 'on' ):
        lat = [s.strip() for s in ifile.readlines()][::-1]
      lat = [float(s) for s in lat]

    with open(ifile_lon, 'r') as ifile:
      lon = [s.strip() for s in ifile.readlines()]
      lon = [float(s) for s in lon]
      
    X, Y = np.meshgrid(lon, lat)

    return X, Y

  def _open_gpv(self, gpv_file, *, endian='little'):
    print('...... Preparating for {}'.format(gpv_file))
    with open(gpv_file, 'rb') as ifile:
      if endian == 'little':
        data = np.fromfile(ifile, dtype='<f', sep = '')
      elif endian == 'big':
        data = np.fromfile(ifile, dtype='>f', sep = '')
    return data
  
  def weight_latitude(self, lat:np.ndarray) -> np.ndarray:
    return np.sqrt(np.cos(np.deg2rad(np.abs(lat))))

class Energy_NORM:
  def __init__(self,dataset):
    ST = setup.Setup(dataset)
    self.nx, self.ny, self.nz, self.mem = ST.set_prm()
    self.press_levels = ST.set_pressure_levels()
    self.ctrl, self.surf = 1, 1
    self.Pr:float = 800.0
    self.Tr:float = 270.0
    self.cp:float = 1004.0
    self.R:float  = 287.0
    self.Lc:float = 2.5*1000000
    self.wq:float = 1.0

  def init_array(self):
    pertb_uwnd_data = np.zeros((self.mem-self.ctrl, self.nz, self.ny, self.nx))
    pertb_vwnd_data = np.zeros((self.mem-self.ctrl, self.nz, self.ny, self.nx))
    pertb_tmp_data  = np.zeros((self.mem-self.ctrl, self.nz, self.ny, self.nx))
    pertb_spfh_data = np.zeros((self.mem-self.ctrl, self.nz, self.ny, self.nx))
    pertb_ps_data   = np.zeros((self.mem-self.ctrl, self.ny, self.nx))
    return pertb_uwnd_data, pertb_vwnd_data, pertb_tmp_data, pertb_spfh_data, pertb_ps_data

  def data_pertb_driver(self,uwnd,vwnd,tmp,hgt,spfh,ps):
    """コントロールランからの摂動の作成
    Args:
      elem(np.ndarray): 各アンサンブルランのデータ(配列の初期がコントロールラン)
    Returns:
      pretb_elem_data(np.ndarray): アンサンブルランのデータからコントロールランデータを引いた擾乱のデータ
    """
    pertb_uwnd_data, pertb_vwnd_data, pertb_tmp_data, pertb_hgt_data, pertb_spfh_data, pertb_ps_data = self.init_array()
    
    for imem in range(self.mem-self.ctrl):
      pertb_uwnd_data[imem,:,:,:] = uwnd[imem,:,:,:] - uwnd[self.mem,:,:,:]
      pertb_vwnd_data[imem,:,:,:] = vwnd[imem,:,:,:] - vwnd[self.mem,:,:,:]
      pertb_tmp_data[imem,:,:,:]  =  tmp[imem,:,:,:] -  tmp[self.mem,:,:,:]
      pertb_spfh_data[imem,:,:,:] = spfh[imem,:,:,:] - spfh[self.mem,:,:,:]
      pertb_ps_data[imem,:,:]     =   ps[imem,:,:]   -   ps[self.mem,:,:]

    return pertb_uwnd_data, pertb_vwnd_data, pertb_tmp_data, pertb_spfh_data, pertb_ps_data

  def calc_dry_EN_NORM(self,u_prime:np.ndarray, v_prime:np.ndarray, tmp_prime:np.ndarray, ps_prime:np.ndarray):
    """乾燥エネルギーノルムの計算
    Args:
      u_prime   (np.ndarray): 東西風のコントロールランからの予測時間における摂動
      v_prime   (np.ndarray): 南北風のコントロールランからの予測時間における摂動
      tmp_prime (np.ndarray): 気温のコントロールランからの予測時間における摂動
      ps_prime  (np.ndarray): 地表面気圧のコントロールランからの予測時間における摂動
    Parameters:
      self.Pr (float) : 経験的に求めた参照気圧. Defaults to 1000 hPa.
      self.Tr (float) : 経験的に求めた参照気温. Defaults to 270 K.
      self.cp (float) : 定圧比熱. Defaults to 1004 J/K*kg.
      self.R  (float) : 気体の状態定数. Defaults to 287.0 J/K*kg.
    Returns:
      dry_energy_norm (np.ndarray): トータル乾燥エネルギーノルム(J/kg)
      constitution -> [緯度, 経度] 
    """

    #wind term
    physical_term = (u_prime)**2+(v_prime)**2
    #matsueda et al. (2014)
    vint_physical_term = self._vint(physical_term,self.press_levels)*0.5
    #enomoto et al. (2015)
    #vint_physical_term = self._vint(physical_term,self.press_levels)/(2*self.Pr)

    #temperture term
    tmp_term = (self.cp/self.Tr)*((tmp_prime)**2)
    #matsueda et al. (2014)
    vint_tmp_term = self._vint(tmp_term,self.press_levels[:])*0.5
    #enomoto et al. (2015)
    #vint_tmp_term = self._vint(tmp_term,self.press_levels[:])/(2*self.Pr)

    #pressure term
    ps_term = (self.R*self.Tr/self.Pr)*(ps_prime**2/self.Pr)*0.5

    #SUM OF TERM
    vint_potential_term = vint_tmp_term + ps_term
    dry_energy_norm = vint_physical_term + vint_potential_term

    return dry_energy_norm, vint_physical_term, vint_potential_term

  def calc_humid_EN_NORM(self, u_prime, v_prime, tmp_prime, spfh_prime, ps_prime):
    """湿潤エネルギーノルムの計算
    Args:
      u_prime   (np.ndarray) : 東西風のコントロールランからの予測時間における摂動
      v_prime   (np.ndarray) : 南北風のコントロールランからの予測時間における摂動
      tmp_prime (np.ndarray) : 気温のコントロールランからの予測時間における摂動
      spfh_prime (np.ndarray): 比湿のコントロールランからの予測時間における摂動
      ps_prime  (np.ndarray) : 地表面気圧のコントロールランからの予測時間における摂動
    Parameters:
      self.Pr (float) : 経験的に求めた参照気圧. Defaults to 1000 hPa.
      self.Tr (float) : 経験的に求めた参照気温. Defaults to 270 K.
      self.cp (float) : 定圧比熱. Defaults to 1004 J/K*kg.
      self.R  (float) : 気体の状態定数. Defaults to 287.0 J/K*kg.
    Returns:
      dry_energy_norm (np.ndarray): トータル乾燥エネルギーノルム(J/kg)
      constitution -> [緯度, 経度] 
    """
    #wind term
    physical_term = (u_prime)**2+(v_prime)**2
    #matsueda et al. (2014)
    vint_physical_term = self._vint(physical_term,self.press_levels)*0.5
    #enomoto et al. (2015)
    #vint_physical_term = self._vint(physical_term,self.press_levels)/(2*self.Pr)

    #temperture term
    tmp_term = (self.cp/self.Tr)*((tmp_prime)**2)
    #matsueda et al. (2014)
    vint_tmp_term = self._vint(tmp_term,self.press_levels[:])*0.5
    #enomoto et al. (2015)
    #vint_tmp_term = self._vint(tmp_term,self.press_levels[:])/(2*self.Pr)

    #humid term
    spfh_term = self.wq*(self.Lc**2)*(spfh_prime**2)/(self.cp*self.Tr)
    vint_spfh_term = self._vint(spfh_term,self.press_levels[:])*0.5

    #pressure term
    ps_term = (self.R*self.Tr/self.Pr)*(ps_prime**2/self.Pr)*0.5

    #SUM OF TERM
    vint_potential_term = vint_tmp_term + vint_spfh_term + ps_term
    humid_energy_norm = vint_physical_term + vint_potential_term

    return humid_energy_norm, vint_physical_term, vint_potential_term

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
    #print('Pressure shape & LEVEL :: ', x_array.shape, press_array)
    y_array = np.empty_like(x_array[0], dtype=np.float32)

    for iy in range(self.ny):
      for ix in range(self.nx):
        y_array[iy,ix] = integrate.simps(x_array[:,iy,ix], press_array[::-1])

    return y_array

  def singular_decomposion(self, array, *, mode=10, try_num=10):
    """特異値分解
    Args:
      array (np.ndarray) : 特異値分解したい行列(n,m) 
      mode (int)    :　計算したいモード数 
      try_num (int) : メモリ上を考慮してトライする回数
    Return:
      U_array (np.ndarray)     : ユニタリ行列(n,n)
      sigma_array (np.ndarray) : 特異値の対角成分(r,r) -> sigma_array*sigma_array/m で(array, array.T)の固有値
      V_array (np.ndarray)     : アジョイント行列(n,m)
    """

    for _ in range(try_num):
      try:
        #U_array, sigma_array, V_array = numpy.linalg.svd(array,full_matrices=True)
        U_array, sigma_array, V_array = scipy.linalg.svd(array)
        #U_array, sigma_array, V_array = scipy.sparse.linalg.svds(array, k=mode-1)
        return_index = 0
        print('..... SUCCESS SINGULAR VECTOR CALCULATION ')
        return return_index, U_array, sigma_array, V_array

      except Exception as e:
        print('   >>> {:04} cycle failed, ONE MORE TRY <<<   '.format(_+1))
        command = ["sleep", "5.0s"]
        res = subprocess.call(command)
      
      else:
        break
    
    else:
      print('..... NOT SUCCESS SINGULAR VECTOR CALCULATION ')
      return_index, _ = 1, np.zeros((1))
      mem = psutil.virtual_memory() 
      total, used, available = mem.total, mem.used, mem.available
      percent_total, percent_used, percent_available = (total/total)*100, (used/total)*100, (available/total)*100

      print(
        '..... CHECK MEMORY PERCENT MAX:{:.2f}, USED:{:.2f}, EMPTY:{:.2f}'.format(percent_total,percent_used,percent_available)
        )
      print('')

      return return_index, _, _, _

  def eigen_decomposion(self, array):
    """固有値問題
    Args:
        array (np.ndarray): 固有値分解したい行列
    Returns:
        eig_val (np.ndarray]): 固有値
        eig_vec (np.ndarray]): 正規化された固有ベクトル
    """

    eig_val, eig_vec = scipy.linalg.eig(array, left=True, right=False, overwrite_a=True, check_finite=True)
    print('..... SUCCESS EIGEN VECTOR CALCULATION ')
    
    #normalize
    #for i in range(len(eig_vec)): 
    #  eig_vec[i] = eig_vec[i]/np.linalg.norm(eig_vec[i])

    return eig_val, eig_vec

  def eigen_vector_normalization(self, eigen_vector):
    """固有ベクトルの総和を１にずる規格化
    Args:
      eigen_vector(np.ndarray): 固有ベクトル
    Returns:
      normalize_eigen_vector(np.ndarray): 正規化された固有ベクトル
    """

    normalize_eigen_vector = statics_tool.normalize(eigen_vector, axis=None) 
    normalize_eigen_vector = normalize_eigen_vector*np.sqrt(1/normalize_eigen_vector.shape[0])
    print('..... CHECK EIGENVECTOR SUM :: {}'.format((normalize_eigen_vector*normalize_eigen_vector).sum()))

    return normalize_eigen_vector

  def region_normalize_norm(self, normalize_region, lon, lat, array):
    normalize_lat_min_index, normalize_lat_max_index, normalize_lon_min_index, normalize_lon_max_index = \
      self.verification_region(lon,lat,
          area_lat_min=normalize_region[1], area_lat_max=normalize_region[0],
          area_lon_min=normalize_region[2], area_lon_max=normalize_region[3]
      )

    normalize_array = np.zeros_like(array)

    normalize_array[normalize_lat_min_index:normalize_lat_max_index+1,normalize_lon_min_index:normalize_lon_max_index+1] =\
      statics_tool.min_max(
        array[ normalize_lat_min_index:normalize_lat_max_index+1,
               normalize_lon_min_index:normalize_lon_max_index+1 ]
      )

    print('MIN :: ', np.min(normalize_array), 'MAX :: ', np.max(normalize_array))

    return normalize_array


  def region_normalize_3d_norm(self, normalize_region, lon, lat, array):
    normalize_lat_min_index, normalize_lat_max_index, normalize_lon_min_index, normalize_lon_max_index = \
      self.verification_region(lon,lat,
          area_lat_min=normalize_region[1], area_lat_max=normalize_region[0],
          area_lon_min=normalize_region[2], area_lon_max=normalize_region[3]
      )

    normalize_array = np.zeros_like(array)

    normalize_array[:,normalize_lat_min_index:normalize_lat_max_index+1,normalize_lon_min_index:normalize_lon_max_index+1] =\
      statics_tool.min_max(
        array[ :,
               normalize_lat_min_index:normalize_lat_max_index+1,
               normalize_lon_min_index:normalize_lon_max_index+1 ]
      )

    print('MIN :: ', np.min(normalize_array), 'MAX :: ', np.max(normalize_array))

    return normalize_array 


