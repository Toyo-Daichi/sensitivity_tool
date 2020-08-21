# -*- coding: utf-8 -*-

""" If you use local env., please check out. """
#import os, conda
#conda_file_dir = conda.__file__
#conda_dir = conda_file_dir.split('lib')[0]
#proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
#os.environ["PROJ_LIB"] = proj_lib

import numpy as np
import matplotlib.pyplot as plt

#my module
import mapping
import readgpv_rish, readgpv_tigge

class Mapping_NORM:
  def __init__(self, dataset, map_prj):
    self.MP = mapping.Mapping(map_prj)
    if dataset == 'WFM' or 'EPSW':
      self.RG = readgpv_rish.ReadGPV(dataset,'_','_')
      self.EN = readgpv_rish.Energy_NORM(dataset)
    elif 'TIGEE' in dataset:
      self.RG = readgpv_tigge.ReadGPV(dataset,'_','_')
      self.EN = readgpv_tigge.Energy_NORM(dataset)

  def spaghetti_diagram_driver(self, data, elem, target_region, level_layer, ft, date):
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode='lcc')
    lon, lat = self.RG.set_coordinate() 
    x, y = self.MP.coord_change(mapp, lon[0:self.EN.ny,:], lat[0:self.EN.ny,:])
    level = self.RG.level[level_layer]

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    for _ in range(1,self.EN.mem-1):
      self.MP.contour(mapp, x, y, data[_, level_layer, 0:self.EN.ny, :], elem='850hPa',colors='blue')

    self.MP.contour(mapp, x, y, data[0, level_layer, 0:self.EN.ny, :], elem='850hPa', linewidths=2.0)

    self.MP.title('{} SPAGHETTI DIAGRAM level={}hPa FT={}hr INIT = {}'.format(elem,level,ft,date))
    plt.show()

  def pertubation_driver(self, pertb_data, elem, target_region, level_layer, ft, date, prj, imem):
    """
    コントロールランからの差(摂動)を作成
    """
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode=prj)
    lon, lat = self.RG.set_coordinate() 
    lat_size=37
    x, y = self.MP.coord_change(mapp, lon[0:lat_size,:], lat[0:lat_size,:])
    level = self.RG.press_levels[level_layer]

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    self.MP.diff_contourf(mapp, x, y, pertb_data[level_layer, 0:lat_size, :], elem='normalize')
    self.MP.title('{} PERTUBATION level={}hPa, FT={}hr INIT = {}'.format(elem,level,ft,date))
    self.MP.saving('pertubation_{:03}'.format(imem+1),'./work/')
    plt.close("all")

  def main_norm_driver(self, dry_energy_norm, hgt_data, target_region, ft, date):
    """Draw sensitivity area @dry enegy norm"""
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode='lcc')
    lon, lat = self.RG.set_coordinate() 
    x, y = self.MP.coord_change(mapp, lon, lat)
    label_cfmt = 'spread_{}hr'.format(ft)

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    self.MP.norm_contourf(mapp, x, y, np.average(dry_energy_norm[1:26:2],axis=0), label=label_cfmt)
    self.MP.contour(mapp, x, y, hgt_data[0], elem='850hPa')
    self.MP.title('TE spread [ J/kg ] FT={}hr, INIT={}'.format(ft,date))
    plt.show()

  def each_elem_norm_dry_rish_driver(self, 
    pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, # pertbuation data
    press_levels, target_region, ft, date 
    ):

    size_x, size_y = 16, 18
    row, column = 4, 4
    label_cfmt = 'spread_{}hr'.format(ft)

    fig = plt.figure(figsize = (size_x, size_y))
    lon, lat = self.RG.set_coordinate() 
    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    # UWND
    for index, level in enumerate(press_levels):
      ax = fig.add_subplot(row,column,1+4*index)
      mapp = self.MP.base(projection_mode='lcc')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, pertb_uwnd[index]**2 ,label=label_cfmt)
      self.MP.title('UWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH UWND level {}hPa'.format(level))

    # VWND
    for index, level in enumerate(press_levels):
      ax = fig.add_subplot(row,column,1+2*index)
      mapp = self.MP.base(projection_mode='lcc')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, pertb_vwnd[index]**2 ,label=label_cfmt)
      self.MP.title('VWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH VWND level {}hPa'.format(level))

    #TMP
    for index, level in enumerate(press_levels[1:]):
      ax = fig.add_subplot(row,column,3)
      mapp = self.MP.base(projection_mode='lcc')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, (pertb_tmp[index]**2)*(1004/270) ,label=label_cfmt)
      self.MP.title('TMP [ J/kg ] {}hPa '.format(level))
      print('..... FINISH TMP level {}hPa'.format(level))

    #SLP
    ax = fig.add_subplot(row,column,16)
    mapp = self.MP.base(projection_mode='lcc')
    x, y = self.MP.coord_change(mapp, lon, lat)

    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    self.MP.norm_contourf(mapp, x, y, ((pertb_slp[0]**2)/700)*(287*270) ,label=label_cfmt)
    self.MP.title('SLP [ J/kg ] {}hPa '.format(level))
    print('..... FINISH SLP level {}hPa'.format(level))

    #fig.title(' Ensemble based Sensitivity Analysis (JMA) Valid time: {}, Target region: JPN'.format(date))
    plt.show()
