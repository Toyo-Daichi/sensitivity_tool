# -*- coding: utf-8 -*-
"""
Created from 2020.8.4
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

"""Ensemble spread anaysis(for TIGGE)
  検証時刻でのトータルエネルギーノルムを作成。
  詳細は, README.md or Enomoto et al. (2015) or Matsueda et al. (2014) に記載されている.
"""

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = '2018', '07', '03', '12', '24'
  date = yyyy+mm+dd+hh
  center = 'JMA'
  dataset = 'TIGGE_' + center
  mode = 'dry' # 'dry' or 'humid'
  map_prj, set_prj = 'CNH', 'lcc' # 'CNH', 'lcc' or 'ALL', 'cyl'
  target_region = ( 25, 50, 125, 150 ) # lat_min/max, lon_min/max

  """Class & parm set """
  RG = readgpv_tigge.ReadGPV(dataset,date,ft)
  EN = readgpv_tigge.Energy_NORM(dataset)
  MP = mapping_draw_NORM.Mapping_NORM(dataset,map_prj)

  """Making pretubation data"""
  indir = '/work3/daichi/Data/TIGGE/' + center + '/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, spfh_data, ps_data = RG.data_read_ft_driver(indir+date)
  pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_hgt, pertb_spfh, pertb_ps = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,hgt_data, spfh_data,ps_data)   
  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
      pertb_tmp[imem,i_level,:,:]  = pertb_tmp[imem,i_level,:,:]*weight_lat
      pertb_spfh[imem,i_level,:,:] = pertb_spfh[imem,i_level,:,:]*weight_lat
      pertb_ps[imem,i_level,:,:]   = pertb_ps[imem,i_level,:,:]*weight_lat

  """ Draw function SPREAD """
  #level_layer = 2
  #MP.spaghetti_diagram_driver(hgt_data,RG.elem[2],target_region,level_layer,ft,date,prj=set_prj,center=center,elem_cfmt='850hPa')
  #for imem in range(EN.mem-EN.ctrl):
  #  MP.pertubation_driver(pertb_hgt[imem],RG.elem[2],target_region,level_layer,ft,date,imem,prj=set_prj,center=center)

  #print('Normal END')
  #sys.exit()

  """Calc. Energy NORM"""
  energy_norm    = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
  physical_term  = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
  potential_term = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))

  lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
    EN.verification_region(lon,lat,
        area_lat_min =50, area_lat_max =20,
        area_lon_min =120, area_lon_max =150
    )
    
  lat_grd = lat_max_index-lat_min_index +1
  lon_grd = lon_max_index-lon_min_index +1
  dims = lat_grd*lon_grd

  for imem in range(0,EN.mem-EN.ctrl,1):
    if mode is 'dry':
      energy_norm[imem], physical_term[imem], potential_term[imem] =\
        EN.calc_dry_EN_NORM(pertb_uwnd[imem],pertb_vwnd[imem],pertb_tmp[imem],pertb_ps[imem,EN.surf-1])
    elif mode is 'humid':
      energy_norm[imem], physical_term[imem], potential_term[imem] =\
        EN.calc_humid_EN_NORM(pertb_uwnd[imem],pertb_vwnd[imem],pertb_tmp[imem],pertb_spfh[imem],pertb_ps[imem,EN.surf-1])

    print('')
    print('..... Check Vertification area Norm SUM {:02} {}'.format(
      imem+1, np.sum(energy_norm[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
    )
    print('..... Check Vertification area physical_term {:02} {}'.format(
      imem+1, np.sum(physical_term[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
    )
    print('..... Check Vertification area potential_term {:02} {}'.format(
      imem+1, np.sum(potential_term[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
    )

  print('')
  print('..... @ MAKE EMSEMBLE MEMBER SPREAD : MODE {} @'.format(mode))
  print('')

  MP.main_norm_driver(np.average(energy_norm,axis=0),np.average(hgt_data,axis=0),target_region, ft, date, center=center, TE_mode=mode)

  print('Normal END')
