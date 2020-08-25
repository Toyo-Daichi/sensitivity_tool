# -*- coding: utf-8 -*-
"""
Created from 2020.8.16
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import warnings
warnings.filterwarnings('ignore')

#my_module
import mapping_draw_NORM
import readgpv_rish
import statics_tool

class Anl_ENSVSA:
  """
    Basic info.
      Ensemble Singular Vector sensitivity anaysis(ここでは、固有値ベクトルを求める)
      詳細は, README.md or Enomoto et al. (2015)に記載されている.
    Note:
      Z.T G Z = (mems, dims) (dims, mems) = (mems, mems) 
      -> 行列サイズから容易に固有値問題を解くことができる。
  """

  def __init__(self):
    pass
  
  def init_Z_array(self, dims_xy):
    dims = (2*EN.nz*dims_xy)+(EN.nz-EN.surf)*dims_xy+dims_xy # wind+tmp+slp
    Z_array = np.zeros((dims,EN.mem-EN.ctrl))
    return Z_array, dims

  def eigen_value_and_vector_driver(self, dims_xy, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp):
    """ 固有値問題
      Eigen vector
      Args:
        pertb_elem (np.ndarray) : 検証領域の各要素の摂動から, Zの行列を作成
        -> Z.T G Z = (mems, dims) (dims, mems) = (mems, mems) 
        dims_xy (int) : 検証領域の水平グリッド数
      Returns:
        eigen_value (np.ndarray)  : 固有値(各モードの寄与率を計算する際に使用)
        eigen_vector (np.ndarray) : 固有ベクトル(p_vectorに相当)
    """
    Z_array, dims = self.init_Z_array(dims_xy)
    svd_pertb_tmp = pertb_tmp[:]*np.sqrt(EN.cp/EN.Tr)
    svd_pertb_slp = pertb_slp[:,0]*np.sqrt((EN.R*EN.Tr)/EN.Pr)

    for imem in range(EN.mem-EN.ctrl):
      Z_array[(0*dims_xy):(EN.nz*dims_xy),imem] = pertb_uwnd[imem].reshape(-1)
      Z_array[(EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)):(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)),imem] = svd_pertb_tmp[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)):dims,imem] = svd_pertb_slp[imem,0].reshape(-1)

    array = Z_array.T @ Z_array
    eigen_value, eigen_vector = EN.eigen_decomposion(array)

    return eigen_value, eigen_vector

  def making_initial_pertb_array(self, dims_xy, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, p_array, *, mode=10):
    Z_array, dims = self.init_Z_array(dims_xy)
    array = np.zeros((dims,mode))
    for imem in range(EN.mem-EN.ctrl):
      Z_array[(0*dims_xy):(EN.nz*dims_xy),imem] = pertb_uwnd[imem].reshape(-1)
      Z_array[(EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)):(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)),imem] = pertb_tmp[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)):dims,imem] = pertb_slp[imem,0].reshape(-1)

    for _ in range(mode):
      array[:,_] = Z_array @ p_array[:,mode]

    sum_array = np.sum(array, axis=1)
    svd_pertb_uwnd = np.zeros((EN.nz,EN.ny,EN.nx))
    svd_pertb_vwnd = np.zeros((EN.nz,EN.ny,EN.nx))
    svd_pertb_tmp  = np.zeros((EN.nz-EN.surf,EN.ny,EN.nx))
    svd_pertb_slp  = np.zeros((EN.ny,EN.nx))

    svd_pertb_uwnd[:,:,:] = sum_array[(0*dims_xy):(EN.nz*dims_xy)].reshape(EN.nz,EN.ny,EN.nx)
    svd_pertb_vwnd[:,:,:] = sum_array[(EN.nz*dims_xy):(2*(EN.nz*dims_xy))].reshape(EN.nz,EN.ny,EN.nx)
    svd_pertb_tmp[:,:,:] = sum_array[(2*(EN.nz*dims_xy)):(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy))].reshape(EN.nz-EN.surf,EN.ny,EN.nx)
    svd_pertb_slp[:,:] = sum_array[(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)):dims].reshape(EN.surf,EN.ny,EN.nx)

    return svd_pertb_uwnd, svd_pertb_vwnd, svd_pertb_tmp, svd_pertb_slp

  def sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, target_region):
    """Total Energy NORM を計算する """

    dry_energy_norm = np.zeros((EN.ny,EN.nx))
    physical_term   = np.zeros((EN.ny,EN.nx))
    potential_term  = np.zeros((EN.ny,EN.nx))

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    
    lat_grd = lat_max_index-lat_min_index +1
    lon_grd = lon_max_index-lon_min_index +1
    dims = lat_grd*lon_grd

    dry_energy_norm, physical_term, potential_term = EN.calc_dry_EN_NORM(
      pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp
      )

    print('')
    print('..... Check Vertification area Norm SUM {} {}'.format(
      'SVD MODE', np.sum(dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area physical_term {} {}'.format(
      'SVD MODE', np.sum(physical_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area potential_term {} {}'.format(
      'SVD MODE', np.sum(potential_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))

    print('')
    #print(dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])

    return dry_energy_norm, physical_term, potential_term

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, init, ft = '2015', '09', '09', '12', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'EPSW' # 'WFM' or 'EPSW'
  map_prj, set_prj = 'CNH', 'lcc'
  target_region = ( 25, 50, 125, 150 ) # lat_min/max, lon_min/max
  mode = 10

  """Class & parm set """
  DR = Anl_ENSVSA()
  RG = readgpv_rish.ReadGPV(dataset,date,ft)
  EN = readgpv_rish.Energy_NORM(dataset)
  MP = mapping_draw_NORM.Mapping_NORM(dataset, map_prj)

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  """Making pretubation data Vertificate TIME"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_ft_driver(indir+date[0:8])


  #weight on latitude
  print('')
  print('..... @ MAKE Pertubation array & REGION Extraction @')
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    pertb_slp[imem,0,:,:] = pertb_slp[imem,0,:,:]*weight_lat

  """Target region setting """
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
  pertb_slp  = pertb_slp[:,0,lat_min_index:lat_max_index+1,lon_min_index:lon_max_index+1]

  print('')
  print('..... @ MAKE Eigen VALUE & VECTOR @')
  eigen_value, p_array = DR.eigen_value_and_vector_driver(dims_xy,pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp)
  print('')

  print('..... @ CHECK Eigen VALUE @')
  print(eigen_value)
  print('')

  """Calc. Sensitivity Region"""
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_init_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)
  dims_xy = EN.ny*EN.nx 

  #weight on latitude
  print('')
  print('..... @ MAKE Pertubation array & REGION Extraction @')
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    pertb_slp[imem,:,:] = pertb_slp[imem,:,:]*weight_lat

  print('')
  print('..... @ MAKE SENSITIVITY REGION @')
  svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_slp = DR.making_initial_pertb_array(dims_xy,pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp,p_array,mode=mode)
  energy_norm, _, _ = DR.sensitivity_driver(svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_slp,target_region)
  contribute = float((np.sum(eigen_value[:mode+1])/np.sum(eigen_value))*100)

  #normalize
  print('..... @ MAKE NORMALIZE ENERGY NORM @')
  print('')
  #normal_energy_norm = statics_tool.normalize(energy_norm)
  normal_energy_norm = statics_tool.min_max(energy_norm)

  """ Draw function NORM """
  MP.main_norm_driver(normal_energy_norm,np.average(hgt_data,axis=0),target_region,ft,date,label_cfmt='SVD')
  MP.each_elem_norm_dry_rish_driver(svd_pertb_uwnd,svd_pertb_vwnd,svd_pertb_tmp,svd_pertb_slp,EN.press_levels,target_region,ft,date)

  print('Normal END')
