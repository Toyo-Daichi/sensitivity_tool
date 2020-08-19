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
import mapping
import readgpv_rish

class Anl_SPREAD:
  """Ensemble spread anaysis(for RISH)
    検証時刻でのトータルエネルギーノルムを作成。
    詳細は, README.md or Enomoto et al. (2015)に記載されている.
  """
  
  def __init__(self):
    pass
     
  def main_driver(self, dry_energy_norm, hgt_data, target_region, ft, date):
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
    
    MP.norm_contourf(mapp, x, y, np.average(dry_energy_norm[1:26:2],axis=0), label='spread_{}hr'.format(ft))
    MP.contour(mapp, x, y, hgt_data[0], elem='850hPa')
    MP.title('TE spread [ J/kg ] FT={}hr INIT = {}'.format(ft,date))
    plt.show()

  def spaghetti_diagram_driver(self, data, elem, target_region, level_layer, ft, date):
    fig, ax = plt.subplots()
    mapp = MP.base(projection_mode='lcc')
    lon, lat = RG.set_coordinate() 
    x, y = MP.coord_change(mapp, lon, lat)
    level = RG.press_levels[level_layer]

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    for _ in range(1,EN.mem-1):
      MP.contour(mapp, x, y, data[_][level_layer], elem='850hPa',colors='blue')

    MP.contour(mapp, x, y, data[0][level_layer], elem='850hPa', linewidths=2.0)

    MP.title('{} SPAGHETTI DIAGRAM level={}hPa FT={}hr INIT = {}'.format(elem,level,ft,date))
    plt.show()

  def pertubation_driver(self, pertb_data, elem, target_region, level_layer, ft, date):
    fig, ax = plt.subplots()
    mapp = MP.base(projection_mode='lcc')
    lon, lat = RG.set_coordinate() 
    x, y = MP.coord_change(mapp, lon, lat)
    level = RG.press_levels[level_layer]

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    MP.diff_contourf(mapp, x, y, pertb_data[_][_], elem='normalize')
    MP.title('{} PERTUBATION level={}hPa, FT={}hr INIT = {}'.format(elem,level,ft,date))
    plt.show()

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = '2018', '07', '04', '12', '00'
  date = yyyy+mm+dd+hh
  dataset = 'EPSW' # 'WFM' or 'EPSW'
  target_region = ( 25, 50, 125, 150 ) # lat_min/max, lon_min/max

  """Class & parm set """
  DR = Anl_SPREAD()
  RG = readgpv_rish.ReadGPV(dataset,date,ft)
  EN = readgpv_rish.Energy_NORM(dataset)
  MP = mapping.Mapping('WNH')

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)

  """Making pretubation data"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_ft_driver(indir+date[0:8])
  pertb_uwnd,pertb_vwnd,pertb_tmp,pertb_slp = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,slp_data)   
 
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
    for i_level in range(EN.nz-EN.surf):
      pertb_tmp[imem,i_level,:,:] = pertb_tmp[imem,i_level,:,:]*weight_lat
    pertb_slp[imem,:,:] = pertb_slp[imem,:,:]*weight_lat

  # Draw spaghetti
  level_layer=0
  DR.spaghetti_diagram_driver(hgt_data,RG.elem[2],target_region,level_layer,ft,date)

  # Draw pertubation
  for _ in range(EN.mem):
    DR.pertubation_driver(pertb_uwnd,RG.elem[0],target_region,level_layer,ft,date)

  print('Normal END')

  """Calc. dry Energy NORM"""
  dry_energy_norm = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
  physical_term   = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))
  potential_term  = np.zeros((EN.mem-EN.ctrl,EN.ny,EN.nx))

  lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
    EN.verification_region(lon,lat,
        area_lat_min =50, area_lat_max =20,
        area_lon_min =120, area_lon_max =150
    )
    
  lat_grd = lat_max_index-lat_min_index +1
  lon_grd = lon_max_index-lon_min_index +1
  dims = lat_grd*lon_grd

  for imem in range(EN.mem-EN.ctrl):
    dry_energy_norm[imem], physical_term[imem], potential_term[imem] =\
    EN.calc_dry_EN_NORM(pertb_uwnd[imem],pertb_vwnd[imem],pertb_tmp[imem],pertb_slp[imem])

    print('')
    print('..... Check Vertification area Norm SUM {:02} {}'.format(
      imem+1, np.sum(dry_energy_norm[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
    )
    print('..... Check Vertification area physical_term {:02} {}'.format(
      imem+1, np.sum(physical_term[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
    )
    print('..... Check Vertification area potential_term {:02} {}'.format(
      imem+1, np.sum(potential_term[imem, lat_min_index:lat_max_index,lon_min_index:lon_max_index])/dims)
    )

  print('')
  print('..... @ MAKE EMSEMBLE MEMBER SPREAD @')
  print('')

  # Draw Energy norm
  DR.main_driver(dry_energy_norm,np.average(hgt_data,axis=0),target_region, ft, date)
  #DR.main_driver(dry_energy_norm,slp_data[0],target_region, ft, date)

  print('Normal END')
