# -*- coding: utf-8 -*-
"""
Created from 2020.8.9
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
import readgpv_rish
import statics_tool

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
      pertb_elem (np.ndarray) : 各要素の摂動
      target_region (tuple)   : 検証領域の設定
    Returns:
      theta (list) : 各メンバーの重みの割合
    Note:
      || theta || = 1
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

  def _calc_part_plus_pertb(self, region_TE):
    theta = np.zeros((EN.mem-EN.ctrl))
    for imem in range(0, EN.mem-EN.ctrl, 2):
      theta[imem] = region_TE[imem]/np.sum(region_TE[::2])
    return theta


  def sensitivity_driver(self, pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, theta):
    """Total Energy NORM を計算する """

    dry_energy_norm = np.zeros((EN.ny,EN.nx))
    physical_term   = np.zeros((EN.ny,EN.nx))
    potential_term  = np.zeros((EN.ny,EN.nx))

    ave_pertb_uwnd = np.zeros((EN.nz,EN.ny,EN.nx))
    ave_pertb_vwnd = np.zeros((EN.nz,EN.ny,EN.nx))
    ave_pertb_tmp  = np.zeros((EN.nz-EN.surf,EN.ny,EN.nx))
    ave_pertb_slp  = np.zeros((EN.ny,EN.nx))

    for i_level in range(EN.nz):
      ave_pertb_uwnd[i_level,:,:] = EN.weight_average(pertb_uwnd[:,i_level,:,:],theta)
      ave_pertb_vwnd[i_level,:,:] = EN.weight_average(pertb_vwnd[:,i_level,:,:],theta)
    for i_level in range(EN.nz-EN.surf):
      ave_pertb_tmp[i_level,:,:] = EN.weight_average(pertb_tmp[:,i_level,:,:],theta)
    ave_pertb_slp[:,:] = EN.weight_average(pertb_slp[:,0],theta)

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )
    
    lat_grd = lat_max_index-lat_min_index +1
    lon_grd = lon_max_index-lon_min_index +1
    dims = lat_grd*lon_grd

    dry_energy_norm, physical_term, potential_term = EN.calc_dry_EN_NORM(
      ave_pertb_uwnd, ave_pertb_vwnd, ave_pertb_tmp, ave_pertb_slp
      )

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

    return dry_energy_norm, physical_term, potential_term, \
           ave_pertb_uwnd, ave_pertb_vwnd, ave_pertb_tmp, ave_pertb_slp

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, init, ft = '2018', '07', '04', '12', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'EPSW' # 'WFM' or 'EPSW'
  map_prj, set_prj = 'CNH', 'lcc'
  target_region = ( 20, 50, 120, 150 ) # lat_min/max, lon_min/max

  """Class & parm set """
  DR = Anl_ENASA()
  RG = readgpv_rish.ReadGPV(dataset,date,ft)
  EN = readgpv_rish.Energy_NORM(dataset)
  MP = mapping_draw_NORM.Mapping_NORM(dataset, map_prj)

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
    pertb_slp[imem,0,:,:] = pertb_slp[imem,0,:,:]*weight_lat


  print('')
  print('..... @ MAKE EMSEMBLE MEMBER WEIGHT @')
  theta = DR.adjoint_sensitivity_driver(pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp,target_region)
  print('')

  """Calc. Sensitivity Region"""
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_init_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)

  #weight on latitude
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    pertb_slp[imem,0,:,:] = pertb_slp[imem,0,:,:]*weight_lat

  print('')
  print('..... @ MAKE SENSITIVITY REGION @')
  energy_norm, _, _ , ave_pertb_uwnd,ave_pertb_vwnd,ave_pertb_tmp,ave_pertb_slp =\
     DR.sensitivity_driver(pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp,theta)

  #normalize
  print('..... @ MAKE NORMALIZE ENERGY NORM @')
  print('')
  #normal_energy_norm = statics_tool.normalize(energy_norm)
  normal_energy_norm = statics_tool.min_max(energy_norm)

  print('MIN :: 'np.min(normal_energy_norm), 'MAX :: 'np.max(normal_energy_norm))

  """ Draw function NORM """
  MP.main_norm_driver(energy_norm,np.average(hgt_data,axis=0),target_region,ft,date)
  #MP.each_elem_norm_dry_rish_driver(np.average(pertb_uwnd[::2],axis=0),np.average(pertb_vwnd[::2],axis=0),np.average(pertb_tmp[::2],axis=0),np.average(pertb_slp[::2],axis=0),EN.press_levels,target_region,ft,date)
  #MP.each_elem_norm_dry_rish_driver(pertb_uwnd[0],pertb_vwnd[0],pertb_tmp[0],pertb_slp[0],EN.press_levels,target_region,ft,date)

  print('Normal END')
