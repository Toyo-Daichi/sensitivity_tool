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
import readgpv_tigge

class Anl_SPREAD:
  """Ensemble spread anaysis(for TIGGE)
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

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = '2003', '01', '21', '12', '00'
  date = yyyy+mm+dd+hh
  center = 'JMA'
  dataset = 'TIGGE_' + center 
  target_region = ( 25, 50, 125, 150 ) # lat_min/max, lon_min/max

  """Class & parm set """
  DR = Anl_SPREAD()
  RG = readgpv_tigge.readgpv_tigge(dataset,date,ft)
  EN = readgpv_tigge.Energy_NORM(dataset)
  MP = mapping.Mapping('CNH')

  lon, lat = RG.set_coordinate()
  weight_lat = RG.weight_latitude(lat)
  
  """Making pretubation data"""
  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, spfh_data, ps_data = RG.data_read_driver(indir+date[0:8])
  pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps = EN.data_pertb_driver(uwnd_data,vwnd_data,tmp_data,spfh_data,ps_data)   
 
  for imem in range(EN.mem-EN.ctrl):
    for i_level in range(EN.nz):
      pertb_uwnd[imem,i_level,:,:] = pertb_uwnd[imem,i_level,:,:]*weight_lat
      pertb_vwnd[imem,i_level,:,:] = pertb_vwnd[imem,i_level,:,:]*weight_lat
      pertb_tmp[imem,i_level,:,:]  = pertb_tmp[imem,i_level,:,:]*weight_lat
      pertb_spfh[imem,i_level,:,:] = pertb_spfh[imem,i_level,:,:]*weight_lat
      pertb_ps[imem,i_level,:,:]   = pertb_ps[imem,i_level,:,:]*weight_lat

  sys.exit()
  
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


  DR.main_driver(dry_energy_norm,np.average(hgt_data,axis=0),target_region, ft, date)
  #DR.main_driver(dry_energy_norm,slp_data[0],target_region, ft, date)

  print('Normal END')
