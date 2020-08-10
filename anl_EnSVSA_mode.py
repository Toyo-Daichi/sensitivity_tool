# -*- coding: utf-8 -*-
"""
Created from 2020.8.10
@author: Toyo_Daichi
"""

import gc
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import matplotlib.pyplot as plt

#my_module
import mapping
import readgpv
import setup
import statics_tool

class Anl_ENSVSA:
  """Ensemble Singular Vector sensitivity anaysis(特異値ベクトルを求める)
    詳細は, README.md or Enomoto et al. (2015)に記載されている.
  """

  def __init__(self):
    pass
  
  def init_Z_array(self):
    dims_xy = EN.ny*EN.nx
    dims = 2*EN.nz*dims_xy+(EN.nz-EN.surf)*dims_xy+dims_xy # wind+tmp+slp
    Z_array = np.zeros((dims,EN.mem-EN.ctrl))
    return Z_array, dims_xy, dims

  def singular_vector_sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp):
    """singular vector """
    Z_array, dims_xy, dims = self.init_Z_array()
    mode_pertb_uwnd, mode_pertb_vwnd, mode_pertb_tmp, mode_pertb_slp = EN.init_array()

    for imem in range(EN.mem-EN.ctrl):
      Z_array[(0*dims_xy):(EN.nz*dims_xy),imem] = pertb_uwnd[imem].reshape(-1)
      Z_array[(EN.nz*dims_xy):(2*(EN.nz*dims_xy)),imem] = pertb_vwnd[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)):(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)),imem] = pertb_tmp[imem].reshape(-1)
      Z_array[(2*(EN.nz*dims_xy)+((EN.nz-EN.surf)*dims_xy)):dims,imem] = pertb_slp[imem].reshape(-1)

    setup.save_list_ndarray(Z_array,'./','Z')
    U_array, sigma_array, V_array = EN.singular_decomposion(Z_array)
    sys.exit()

    return 

  def sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp):
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

    dry_energy_norm, physical_term, potential_term = EN.calc_dry_EN_NORM(pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp)

    print('')
    print('..... Check Vertification area Norm SUM {} {}'.format(
      'Adjoint', np.sum(dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area physical_term {} {}'.format(
      'Adjoint', np.sum(physical_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area potential_term {} {}'.format(
      'Adjoint', np.sum(potential_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))

    print('')
    #print(dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])

    return dry_energy_norm, physical_term, potential_term

  def draw_driver(self, energy_norm, hgt_data, ft, date, mode, mode_percent):
    """Draw sensitivity area @dry enegy norm"""
    fig, ax = plt.subplots()
    mapp = MP.base(projection_mode='lcc')
    lon, lat = RG.set_coordinate() 
    x, y = MP.coord_change(mapp, lon, lat)

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    

    #vertifcation region
    MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)

    #norm draw
    MP.norm_contourf(mapp, x, y, energy_norm, label='scope')
    MP.contour(mapp, x, y, hgt_data[1], elem='500hPa')
    MP.title('NORMALIZE TE [ J/kg ] mode = 1-{}, contributution {f:.3f}%, FT= {}hr, INIT = {}'.format(mode, mode_percent, ft,date))
    plt.show()

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
  MP = mapping.Mapping('NH')

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  """Making pretubation data Vertificate TIME"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_ft_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)   

  #deallocate normal_data
  del uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data
  gc.collect()
 
  #weight on latitude
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    pertb_slp[imem,:,:] = pertb_slp[imem,:,:]*weight_lat


  print('')
  print('..... @ MAKE SINGULAR VECTOR @')
  mode_pertb_uwnd,mode_pertb_vwnd,mode_pertb_tmp,mode_pertb_slp =\
     DR.singular_vector_sensitivity_driver(pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp)
  print('')

  """Calc. Sensitivity Region"""
  #print('')
  #print('..... @ MAKE SENSITIVITY REGION @')
  #energy_norm, _, _ = DR.sensitivity_driver(mode_pertb_uwnd,mode_pertb_vwnd,#mode_pertb_tmp,mode_pertb_slp)
  #
  ##normalize
  #print('..... @ MAKE NORMALIZE ENERGY NORM @')
  #print('')
  #normal_energy_norm = statics_tool.normalize(energy_norm)
  #
  #DR.draw_driver(normal_energy_norm,np.average(hgt_data,axis=0),ft,date)

  print('Normal END')
