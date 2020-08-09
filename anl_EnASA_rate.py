# -*- coding: utf-8 -*-
"""
Created from 2020.8.9
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

#my_module
import mapping
import readgpv

class Anl_ENASA:
  """Adjoint sensitivity anaysis(theta を求める)
    初期時刻の擾乱が線形成長すると仮定して, 検証時刻/検証領域での擾乱を計算して各メンバーの重みを計算する.
    各メンバーの重みで作成した初期大気場でTEを計算すると感度領域とも見なす.
    詳細は, README.md or Enomoto et al. (2015)に記載されている.
  """

  def __init__(self):
    pass

  def adjoint_sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, target_region):
    """Adjoint sensitivity anaysis(theta を求める)
    Args:
      target_region(tuple) : 検証領域の設定
    """
    dry_energy_norm = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
    physical_term   = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
    potential_term  = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
    region_TE       = np.zeros((EN.mem-EN.ctrl))

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    
    lat_grd = lat_max_index-lat_min_index +1
    lon_grd = lon_max_index-lon_min_index +1
    dims = lat_grd*lon_grd

    for imem in range(EN.mem-EN.ctrl):
      dry_energy_norm[imem], physical_term[imem], potential_term[imem] =\
        EN.calc_dry_EN_NORM(pertb_uwnd[imem],pertb_vwnd[imem],pertb_tmp[imem],pertb_slp[imem])

      region_TE[imem] = np.sum(dry_energy_norm[imem,lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
      
      print('')
      print('..... Check Vertification area Norm SUM {:02} {}'.format(
        imem+1, region_TE[imem])
      )
      print('..... Check Vertification area physical_term {:02} {}'.format(
        imem+1, np.sum(physical_term[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
      )
      print('..... Check Vertification area potential_term {:02} {}'.format(
        imem+1, np.sum(potential_term[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
      )

    theta = self._calc_part(region_TE)
    return theta

  def _calc_part(self, region_TE):
    theta = []
    for imem in range(EN.mem-EN.ctrl):
      i_theta = region_TE[imem]/np.sum(region_TE[:])
      theta.append(i_theta)
    return theta

  def sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, theta):
    """Adjoint sensitivity anaysis(theta を求める)
    Args:
      target_region(tuple) : 検証領域の設定
    """
    dry_energy_norm = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
    physical_term   = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
    potential_term  = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
    region_TE       = np.zeros((EN.mem-EN.ctrl))

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    
    lat_grd = lat_max_index-lat_min_index +1
    lon_grd = lon_max_index-lon_min_index +1
    dims = lat_grd*lon_grd

    for imem in range(EN.mem-EN.ctrl):
      dry_energy_norm[imem], physical_term[imem], potential_term[imem] =\
        EN.calc_dry_EN_NORM(pertb_uwnd[imem],pertb_vwnd[imem],pertb_tmp[imem],pertb_slp[imem])

    ave_dry_energy_norm = EN.weight_average(dry_energy_norm,theta)
    ave_physical_term = EN.weight_average(physical_term,theta)
    ave_potential_term = EN.weight_average(potential_term,theta)

    print('')
    print('..... Check Vertification area Norm SUM {} {}'.format(
      'Adjoint', np.sum(ave_dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area physical_term {} {}'.format(
      'Adjoint', np.sum(ave_physical_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))
    print('..... Check Vertification area potential_term {} {}'.format(
      'Adjoint', np.sum(ave_potential_term[lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims
    ))

    print('')
    print(ave_dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index])

    return ave_dry_energy_norm, ave_physical_term, ave_potential_term

  def draw_driver(self, energy_norm, hgt_data):
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
    MP.norm_contourf(mapp, x, y, energy_norm, label='adjoint')
    MP.contour(mapp, x, y, hgt_data[1], elem='500hPa')
    MP.title('TE [ J/kg ] Adjoint sensitivity FT=72hr INIT = 20050902')
    plt.show()

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, init, ft = '2003', '01', '21', '12', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'WFM' # 'WFM' or 'EPSW'
  target_region = ( 20, 50, 120, 150 ) # lat_min/max, lon_min/max

  """Class & parm set """
  DR = Anl_ENASA()
  RG = readgpv.ReadGPV(dataset,date,ft)
  EN = readgpv.Energy_NORM(dataset)
  MP = mapping.Mapping('NH')

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  """Making pretubation data Vertificate TIME"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_ft_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)   
 
  #weight on latitude
  weight_pertb_uwnd, weight_pertb_vwnd, weight_pertb_tmp, weight_pertb_slp = EN.init_array()
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      weight_pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      weight_pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      weight_pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    weight_pertb_slp[imem,:,:] = pertb_slp[imem,:,:]*weight_lat

  print('')
  print('..... @ MAKE EMSEMBLE MEMBER WEIGHT @')
  theta = DR.adjoint_sensitivity_driver(weight_pertb_uwnd,weight_pertb_vwnd,weight_pertb_tmp,weight_pertb_slp,target_region)
  print('')

  """Calc. Sensitivity Region"""
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_init_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)
  
  #weight on latitude
  weight_pertb_uwnd, weight_pertb_vwnd, weight_pertb_tmp, weight_pertb_slp = EN.init_array()
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      weight_pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      weight_pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      weight_pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    weight_pertb_slp[imem,:,:] = pertb_slp[imem,:,:]*weight_lat

  print('')
  print('..... @ MAKE SENSITIVITY REGION @')
  energy_norm, _, _ = DR.sensitivity_driver(weight_pertb_uwnd,weight_pertb_vwnd,weight_pertb_tmp,weight_pertb_slp,theta)

  print('')
  DR.draw_driver(energy_norm,np.average(hgt_data,axis=0))

  print('Normal END')
