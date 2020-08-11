# -*- coding: utf-8 -*-
"""
Created from 2020.8.10
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

class Z_ARRAY:
  """
  Ensemble Singular Vector sensitivity anaysis(特異値ベクトルを求める) 1st STEP.
  メモリー量を減少させるために分けて行う。式の導出などの詳細は, Enomoto et al. (2015)に記載されている.
  """

  def __init__(self):
    pass
  
  def init_Z_array(self):
    dims_xy = EN.ny*EN.nx
    dims = 2*EN.nz*dims_xy+(EN.nz-EN.surf)*dims_xy+dims_xy # wind+tmp+slp
    Z_array = np.zeros((dims,EN.mem-EN.ctrl))
    return Z_array, dims_xy, dims

  def make_z_array(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp,date,ft):
    """singular vector """
    Z_array, dims_xy, dims = self.init_Z_array()
    mode_pertb_uwnd, mode_pertb_vwnd, mode_pertb_tmp, mode_pertb_slp = EN.init_array()

    for imem in range(EN.mem-EN.ctrl):
      Z_array[(0*dims_xy):(EN.nz*dims_xy),imem] = pertb_uwnd[imem].reshape(-1)
      Z_array[(EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)):(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)),imem] = pertb_tmp[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)):dims,imem] = pertb_slp[imem].reshape(-1)

    setup.save_list_ndarray(Z_array,'./work/','Z_array_{}_{}hr'.format(date,ft))

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, init, ft = '2003', '08', '05', '12', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'WFM' # 'WFM' or 'EPSW'

  """Class & parm set """
  DR = Z_ARRAY()
  RG = readgpv.ReadGPV(dataset,date,ft)
  EN = readgpv.Energy_NORM(dataset)

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  """Making pretubation data Vertificate TIME"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_ft_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)   

  #weight on latitude
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    pertb_slp[imem,:,:] = pertb_slp[imem,:,:]*weight_lat

  print('')
  print('..... @ MAKE Z_array @')
  DR.make_z_array(pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp,date,ft)
  print('')

  print('Normal END')
