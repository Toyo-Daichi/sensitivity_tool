# -*- coding: utf-8 -*-
"""
Created from 2020.8.26
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

#my_module
import mapping_draw_NORM
import readgpv_tigge
import statics_tool

class Anl_ENSVSA:
  """Ensemble Singular Vector sensitivity anaysis(ここでは、固有値ベクトルを求める for TIGGE)
    Z.T G Z = (mems, dims) (dims, mems) = (mems, mems) 
    -> 行列サイズから容易に固有値問題を解くことができる。
    詳細は, README.md or Enomoto et al. (2015)に記載されている.
  """

  def __init__(self):
    pass
  
  def init_Z_array(self, dims_xy, *, mode='dry'):
    if mode is 'dry':
      dims = 3*EN.nz*dims_xy + dims_xy #u,v,tmp,ps
      Z_array = np.zeros((dims,EN.mem-EN.ctrl))
    elif mode is 'humid':
      dims = 4*EN.nz*dims_xy + dims_xy #u,v,tmp,spfh,ps
      Z_array = np.zeros((dims,EN.mem-EN.ctrl))
    return Z_array, dims

  def eigen_value_and_vector_driver(self, dims_xy, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps,*,mode='dry'):
    """eigen vector """
    Z_array, _ = self.init_Z_array(dims_xy,mode=mode)

    if mode is 'dry':
      #multi const
      svd_pertb_tmp  = pertb_tmp[:]*np.sqrt(EN.cp/EN.Tr)
      svd_pertb_ps   = pertb_ps[:,EN.surf-1]*np.sqrt((EN.R*EN.Tr)/EN.Pr)

      for imem in range(EN.mem-EN.ctrl):
        Z_array[(0*EN.nz*dims_xy):(1*(EN.nz*dims_xy)),imem] = pertb_uwnd[imem].reshape(-1)
        Z_array[(1*EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
        Z_array[(2*EN.nz*dims_xy):(3*(EN.nz*dims_xy)),imem] = svd_pertb_tmp[imem].reshape(-1)
        Z_array[(3*EN.nz*dims_xy):(3*(EN.nz*dims_xy)+dims_xy),imem] = svd_pertb_ps[imem,EN.surf-1].reshape(-1)

    elif mode is 'humid':
      #multi const
      svd_pertb_tmp  = pertb_tmp[:]*np.sqrt(EN.cp/EN.Tr)
      svd_pertb_spfh = pertb_spfh[:]*np.sqrt(EN.wq*(EN.Lc**2)/(EN.cp*EN.Tr))
      svd_pertb_ps   = pertb_ps[:,EN.surf-1]*np.sqrt((EN.R*EN.Tr)/EN.Pr)

      for imem in range(EN.mem-EN.ctrl):
        Z_array[(0*EN.nz*dims_xy):(1*(EN.nz*dims_xy)),imem] = pertb_uwnd[imem].reshape(-1)
        Z_array[(1*EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
        Z_array[(2*EN.nz*dims_xy):(3*(EN.nz*dims_xy)),imem] = svd_pertb_tmp[imem].reshape(-1)
        Z_array[(3*EN.nz*dims_xy):(4*(EN.nz*dims_xy)),imem] = svd_pertb_spfh[imem].reshape(-1)
        Z_array[(4*EN.nz*dims_xy):(4*(EN.nz*dims_xy)+dims_xy),imem] = svd_pertb_ps[imem,EN.surf-1].reshape(-1)

    array = Z_array.T @ Z_array

    print(np.diag(array))

    eigen_value, eigen_vector = EN.eigen_decomposion(array)
    normalize_eigen_vector = EN.eigen_vector_normalization(eigen_vector)

    return eigen_value, normalize_eigen_vector

  def making_initial_pertb_array(self, 
    dims_xy, 
    pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps, 
    p_array,
    *, mode='dry', 
    start_eigen_mode=0, end_eigen_node=10
    ):
    Z_array, dims = self.init_Z_array(dims_xy,mode=mode)
    array = np.zeros((dims,end_eigen_mode-start_eigen_mode+1))

    if mode is 'dry':
      for imem in range(EN.mem-EN.ctrl):
        Z_array[(0*EN.nz*dims_xy):(1*(EN.nz*dims_xy)),imem] = pertb_uwnd[imem].reshape(-1)
        Z_array[(1*EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
        Z_array[(2*EN.nz*dims_xy):(3*(EN.nz*dims_xy)),imem] = pertb_tmp[imem].reshape(-1)
        Z_array[(3*EN.nz*dims_xy):(3*(EN.nz*dims_xy)+dims_xy),imem] = pertb_ps[imem,EN.surf-1].reshape(-1)

    elif mode is 'humid':
      for imem in range(EN.mem-EN.ctrl):
        Z_array[(0*EN.nz*dims_xy):(1*(EN.nz*dims_xy)),imem] = pertb_uwnd[imem].reshape(-1)
        Z_array[(1*EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
        Z_array[(2*EN.nz*dims_xy):(3*(EN.nz*dims_xy)),imem] = pertb_tmp[imem].reshape(-1)
        Z_array[(3*EN.nz*dims_xy):(4*(EN.nz*dims_xy)),imem] = pertb_spfh[imem].reshape(-1)
        Z_array[(4*EN.nz*dims_xy):(4*(EN.nz*dims_xy)+dims_xy),imem] = pertb_ps[imem,EN.surf-1].reshape(-1)

    for index, _ in enumerate(range(start_eigen_mode,end_eigen_node+1)):
      array[:,index] = Z_array @ p_array[:,_]

    sum_array = np.sum(array, axis=1)
    svd_pertb_uwnd = np.zeros((EN.nz,EN.ny,EN.nx))
    svd_pertb_vwnd = np.zeros((EN.nz,EN.ny,EN.nx))
    svd_pertb_tmp  = np.zeros((EN.nz,EN.ny,EN.nx))
    svd_pertb_spfh = np.zeros((EN.nz,EN.ny,EN.nx))
    svd_pertb_ps   = np.zeros((EN.ny,EN.nx))

    if mode is 'dry':
      svd_pertb_uwnd[:,:,:] = sum_array[(0*EN.nz*dims_xy):(1*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_vwnd[:,:,:] = sum_array[(1*EN.nz*dims_xy):(2*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_tmp[:,:,:]  = sum_array[(2*EN.nz*dims_xy):(3*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_ps[:,:]     = sum_array[(3*EN.nz*dims_xy):(3*EN.nz*dims_xy)+dims_xy].reshape(EN.ny,EN.nx)

    elif mode is 'humid':
      svd_pertb_uwnd[:,:,:] = sum_array[(0*EN.nz*dims_xy):(1*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_vwnd[:,:,:] = sum_array[(1*EN.nz*dims_xy):(2*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_tmp[:,:,:]  = sum_array[(2*EN.nz*dims_xy):(3*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_spfh[:,:,:] = sum_array[(3*EN.nz*dims_xy):(4*EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
      svd_pertb_ps[:,:]     = sum_array[(4*EN.nz*dims_xy):(4*EN.nz*dims_xy)+dims_xy].reshape(EN.ny,EN.nx)

    return svd_pertb_uwnd, svd_pertb_vwnd, svd_pertb_tmp, svd_pertb_spfh, svd_pertb_ps

  def sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps, target_region, *, mode='dry'):
    """Total Energy NORM を計算する """

    energy_norm    = np.zeros((EN.ny,EN.nx))
    physical_term  = np.zeros((EN.ny,EN.nx))
    potential_term = np.zeros((EN.ny,EN.nx))

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    
    lat_grd = lat_max_index-lat_min_index +1
    lon_grd = lon_max_index-lon_min_index +1
    dims = lat_grd*lon_grd

    if mode is 'dry':
      energy_norm, physical_term, potential_term = EN.calc_dry_EN_NORM(
        pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_ps
        )
    elif mode is 'humid':
      energy_norm, physical_term, potential_term = EN.calc_dry_EN_NORM(
        pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps
        )


    print('')
    print('..... Check Vertification area Norm SUM {} {}'.format(
      'SVD MODE', np.sum(energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area physical_term {} {}'.format(
      'SVD MODE', np.sum(physical_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area potential_term {} {}'.format(
      'SVD MODE', np.sum(potential_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))

    print('')

    return energy_norm, physical_term, potential_term
    
if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, init, ft = '2015', '09', '10', '12', '00', '24'
  date = yyyy+mm+dd+hh
  center = 'ECMWF'
  dataset = 'TIGGE_' + center + '_pertb_minus' #'_pertb_plus/minus' or '' 
  mode = 'humid' # 'dry' or 'humid'
  map_prj, set_prj = 'ASIAJPN', 'lcc' # 'CNH', 'lcc' or 'ALL', 'cyl'
  #target_region = ( 30, 37.5, 127.5, 137.5 ) # WJPN
  target_region = ( 32.5, 40, 135, 142.5 )    # EJPN
  #target_region = ( 25, 50, 125, 150 )       # JPN
  start_eigen_mode, end_eigen_mode = 0, 12 #default is 0/9 -> 1-10 mode. 
  normalize_set = 'on' # 'on' or 'on_full' or 'off'
  normalize_region = ( 20, 50, 118.75, 152.5 ) # lat_min/max, lon_min/max

  """Class & parm set """
  DR = Anl_ENSVSA()
  RG = readgpv_tigge.ReadGPV(dataset,date,ft)
  EN = readgpv_tigge.Energy_NORM(dataset)
  MP = mapping_draw_NORM.Mapping_NORM(dataset,date,map_prj)

  """Making pretubation data Vertificate TIME"""
  indir = '/work3/daichi/Data/TIGGE/' + center + '/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, spfh_data, ps_data = RG.data_read_ft_driver(indir+date)
  pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_hgt, pertb_spfh, pertb_ps = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,hgt_data, spfh_data,ps_data)
  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  print('')
  print('..... @ MAKE Pertubation array & REGION Extraction MODE:{} @'.format(mode))
  
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
      pertb_tmp[imem,i_level,:,:]  = pertb_tmp[imem,i_level,:,:]*weight_lat
      pertb_spfh[imem,i_level,:,:] = pertb_spfh[imem,i_level,:,:]*weight_lat
    pertb_ps[imem,EN.surf-1,:,:]   = pertb_ps[imem,EN.surf-1,:,:]*weight_lat

  """Target region setting"""
  lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    
  lat_grd = lat_max_index-lat_min_index +1
  lon_grd = lon_max_index-lon_min_index +1
  dims_xy = lat_grd*lon_grd

  pertb_uwnd = pertb_uwnd[:,:,lat_min_index:lat_max_index+1,lon_min_index:lon_max_index+1]
  pertb_vwnd = pertb_vwnd[:,:,lat_min_index:lat_max_index+1,lon_min_index:lon_max_index+1]
  pertb_tmp  = pertb_tmp[:,:,lat_min_index:lat_max_index+1,lon_min_index:lon_max_index+1]
  pertb_spfh = pertb_spfh[:,:,lat_min_index:lat_max_index+1,lon_min_index:lon_max_index+1]
  pertb_ps   = pertb_ps[:,EN.surf-1,lat_min_index:lat_max_index+1,lon_min_index:lon_max_index+1]

  print('')
  print('..... @ MAKE Eigen VALUE & VECTOR @')
  eigen_value, p_array = DR.eigen_value_and_vector_driver(dims_xy,pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_spfh,pertb_ps,mode=mode)
  print('')

  print('..... @ CHECK Eigen VALUE & NORMALIZE VECTOR @')
  print(eigen_value)
  print('')

  """Calc. Sensitivity Region"""
  uwnd_data, vwnd_data, hgt_data, tmp_data, spfh_data, ps_data = RG.data_read_init_driver(indir+date)
  pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_hgt, pertb_spfh, pertb_ps = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,hgt_data, spfh_data,ps_data)
  dims_xy = EN.ny*EN.nx 
  
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
      pertb_tmp[imem,i_level,:,:]  = pertb_tmp[imem,i_level,:,:]*weight_lat
      pertb_spfh[imem,i_level,:,:] = pertb_spfh[imem,i_level,:,:]*weight_lat
    pertb_ps[imem,EN.surf-1,:,:]   = pertb_ps[imem,EN.surf-1,:,:]*weight_lat
  
  print('')
  print('..... @ MAKE Pertubation array & REGION Extraction @')

  print('')
  print('..... @ MAKE SENSITIVITY REGION @')
  svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_spfh,svd_pertb_ps =\
  DR.making_initial_pertb_array(
    dims_xy,pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_spfh,pertb_ps,p_array,
    mode=mode,start_eigen_mode=start_eigen_mode,end_eigen_node=end_eigen_mode
    )

  energy_norm, _, _ = DR.sensitivity_driver(svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_spfh,svd_pertb_ps,target_region)
  contribute = float((np.sum(eigen_value[start_eigen_mode:end_eigen_mode+1])/np.sum(eigen_value))*100)

  #normalize
  if 'on' in normalize_set:
    print('..... @ MAKE NORMALIZE ENERGY NORM @')
    print('')
    
    if normalize_set is 'on_full':
      #normalize full ver.
      #energy_norm = statics_tool.normalize(energy_norm)
      energy_norm = statics_tool.min_max(energy_norm)
      print('MIN :: ', np.min(energy_norm), 'MAX :: ', np.max(energy_norm))
    
    elif normalize_set is 'on':
      #normalize region ver.
      energy_norm = EN.region_normalize_norm(normalize_region, lon, lat, energy_norm)

  """ Draw function NORM """
  MP.main_norm_driver(
    energy_norm,np.average(hgt_data,axis=0),target_region,ft,date,
    prj=set_prj,label_cfmt='SVD',center=center,TE_mode=mode,
    start_mode=start_eigen_mode+1, end_mode=end_eigen_mode+1, contribute=contribute,
    normalize_set=normalize_set,normalize_region=normalize_region
    )
  
  #MP.each_elem_norm_dry_tigge_driver(
  # svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_ps,
  # target_region,ft,date,
  # center=center,TE_mode=mode
  # )
  
  #MP.each_elem_norm_humid_tigge_driver(
  #  svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_spfh,svd_pertb_ps,
  #  target_region,ft,date,
  #  center=center,TE_mode=mode,
  #  normalize_set=normalize_set,normalize_region=normalize_region
  #  )

  print('Normal END')
