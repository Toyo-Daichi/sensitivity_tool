# -*- coding: utf-8 -*-
"""
Created from 2020.8.16
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import matplotlib.pyplot as plt

#my_module
import readgpv
import setup
import statics_tool

class Anl_ENSVSA:
  """
    Basic info.
      nsemble Singular Vector sensitivity anaysis(ここでは、固有値ベクトルを求める)
      詳細は, README.md or Enomoto et al. (2015)に記載されている.
    Note:
      Z.T G Z = (mems, dims) (dims, mems) = (mems, mems) 
      -> 行列サイズから容易に固有値問題を解くことができる。
  """

  def __init__(self):
    pass
  
  def init_Z_array(self,dims_xy):
    dims = (2*EN.nz*dims_xy)+(EN.nz-EN.surf)*dims_xy+dims_xy # wind+tmp+slp
    Z_array = np.zeros((dims,EN.mem-EN.ctrl))
    return Z_array, dims

  def singular_vector_sensitivity_driver(self, dims_xy, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp,date,ft):
    """singular vector """
    Z_array, dims = self.init_Z_array(dims_xy)
    svd_pertb_tmp = pertb_tmp[:]*np.sqrt(EN.cp/EN.Tr)
    svd_pertb_slp = pertb_slp[:,0]*np.sqrt((EN.R*EN.Tr)/EN.Pr)

    for imem in range(EN.mem-EN.ctrl):
      Z_array[(0*dims_xy):(EN.nz*dims_xy),imem] = pertb_uwnd[imem].reshape(-1)
      Z_array[(EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)):(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)),imem] = svd_pertb_tmp[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)):dims,imem] = svd_pertb_slp[imem,0].reshape(-1)

    array = Z_array.T @ Z_array
    eig_val, eig_vec = EN.eigen_decomposion(array)

    return eig_val, eig_vec
    
if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, init, ft = '2003', '08', '05', '12', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'WFM' # 'WFM' or 'EPSW'
  target_region = ( 25, 50, 125, 150 ) # lat_min/max, lon_min/max

  """Class & parm set """
  DR = Anl_ENSVSA()
  RG = readgpv.ReadGPV(dataset,date,ft)
  EN = readgpv.Energy_NORM(dataset)

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  """Making pretubation data Vertificate TIME"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_ft_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)

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
  eig_val, eig_vec = DR.singular_vector_sensitivity_driver(dims_xy,pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp,date,ft)
  print('')

  print(eig_val)
  print('Normal END')
