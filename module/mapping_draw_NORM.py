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
import readgpv_rish, readgpv_tigge, readgpv_nhm
import statics_tool

class Mapping_NORM:
  def __init__(self, dataset, date, map_prj):
    self.MP = mapping.Mapping(map_prj)
    if dataset == 'WFM' or dataset == 'EPSW':
      self.RG = readgpv_rish.ReadGPV(dataset,'_','_')
      self.EN = readgpv_rish.Energy_NORM(dataset)
    elif 'TIGGE' in dataset:
      self.RG = readgpv_tigge.ReadGPV(dataset,date,'_')
      self.EN = readgpv_tigge.Energy_NORM(dataset)
    elif 'NHM' in dataset:
      self.RG = readgpv_nhm.ReadGPV(dataset)
      self.EN = readgpv_nhm.Energy_NORM(dataset)

  def spaghetti_diagram_driver(self, 
    data, elem, target_region, level_layer, ft, date, 
    *, prj='lcc', center='JMA', elem_cfmt='500hPa', lat_half_on=1, level_fix=0
    ):
    """(基本は)500or850hPaのスパゲッティ図"""
    
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode=prj)
    lon, lat = self.RG.set_coordinate() 

    if lat_half_on == 0:
      lat_size = self.EN.ny
    elif lat_half_on == 1:
      lat_size = self.EN.ny // 2

    x, y = self.MP.coord_change(mapp, lon[0:lat_size,:], lat[0:lat_size,:])
    
    # if 'WFM' or 'EPSW'
    if level_fix == 1:
      level = self.RG.press_levels[level_layer+1]
    else:
      level = self.RG.press_levels[level_layer]

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    for _ in range(1,self.EN.mem-1):
      self.MP.contour(mapp, x, y, data[_, level_layer, 0:lat_size, :], elem=elem_cfmt,colors='blue')

    self.MP.contour(mapp, x, y, data[0, level_layer, 0:lat_size, :], elem=elem_cfmt, linewidths=2.0, font_on=1)

    self.MP.title('{} SPAGHETTI DIAGRAM level={}hPa FT={}hr INIT={} CENTER={}'.format(elem,level,ft,date,center),fontsize=7.5)
    #plt.show()
    self.MP.saving('{}_spaghetthi_diagram'.format(center),'./work/')

  def pertubation_driver(self,
    pertb_data, elem, target_region, level_layer, ft, date, imem,  
    *, prj='lcc', center='JMA', elem_cfmt='HGT', lat_half_on=1, level_fix=0
    ):
    """(基本は)500or850hPaのコントロールランからの差(摂動)図"""
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode=prj)
    lon, lat = self.RG.set_coordinate() 
    
    if lat_half_on == 0:
      lat_size = self.EN.ny
    elif lat_half_on == 1:
      lat_size = self.EN.ny // 2

    x, y = self.MP.coord_change(mapp, lon[0:lat_size,:], lat[0:lat_size,:])

    # if 'WFM' or 'EPSW'
    if level_fix == 1:
      level = self.RG.press_levels[level_layer+1]
    else:
      level = self.RG.press_levels[level_layer]

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)

    #pertubation
    self.MP.diff_contourf(mapp, x, y, pertb_data[level_layer, 0:lat_size, :], elem='HGT')
    self.MP.hatch_contourf(mapp, x, y, pertb_data[level_layer, 0:lat_size, :], levels=[8.0, 10.0], extend='max')
    self.MP.hatch_contourf(mapp, x, y, pertb_data[level_layer, 0:lat_size, :], levels=[-10.0, -8.0], extend='min')
    self.MP.title('{} PERTUBATION level={}hPa FT={}hr INIT={} CENTER={}'.format(elem,level,ft,date,center),fontsize=7.5)
    self.MP.saving('{}_pertubation_{:03}'.format(center,imem+1),'./work/')
    plt.close("all")

  def main_norm_driver(self,
    energy_norm, hgt_data, target_region, ft, date, 
    *,
    prj='lcc', label_cfmt='spread',  # for title & colorbar
    center='JMA', TE_mode='dry',     # for savefig
    start_mode:int=0, end_mode:int=0, contribute:float=0.0,  # for svd
    normalize_set='off', normalize_region=(0,0,0,0), # for normalize
    ):

    """Draw sensitivity area @dry enegy norm"""
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode=prj)
    lon, lat = self.RG.set_coordinate() 
    x, y = self.MP.coord_change(mapp, lon, lat)
    
    if label_cfmt == 'spread':
      label_cfmt = 'spread_{}hr'.format(ft)
      title_cfmt = 'TE spread [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      save_cfmt = '{}_TE_{}_{}'.format(center,TE_mode,label_cfmt)
    elif label_cfmt == 'adjoint':
      label_cfmt = 'adjoint'
      title_cfmt = ' TE Adjoint sensitivity [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      #title_cfmt = ' NORMALIZE TE Adjoint sensitivity [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      save_cfmt = '{}_TE_{}_{}_{}hr'.format(center,TE_mode,label_cfmt,ft)
    elif label_cfmt == 'SVD':
      label_cfmt = 'svd'
      title_cfmt = ' TE MODE {}-{} contribute:{}% sensitivity [ J/kg ] FT={}hr, INIT={}'.format(start_mode, end_mode,int(contribute),ft,date)
      #title_cfmt = ' NORMALIZE TE MODE sensitivity [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      save_cfmt = '{}_TE_{}_{}_{}hr'.format(center,TE_mode,label_cfmt,ft)

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    if normalize_set == 'on':
      normal_lat_min_index, normal_lat_max_index, normal_lon_min_index, normal_lon_max_index = \
        self.EN.verification_region(lon,lat,
            area_lat_min=normalize_region[1], area_lat_max=normalize_region[0],
            area_lon_min=normalize_region[2], area_lon_max=normalize_region[3]
        )

      self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')

    
    self.MP.norm_contourf(mapp, x, y, energy_norm, label=label_cfmt)
    self.MP.contour(mapp, x, y, hgt_data[1], elem='850hPa', font_on=1)
    self.MP.title(title_cfmt, fontsize=8)
    self.MP.saving(save_cfmt,'./work/')
    plt.show()

  def main_norm_driver_nhm(self,
    energy_norm, hgt_data, target_region, ft, date, 
    *,
    prj='lcc', label_cfmt='spread',  # for title & colorbar
    center='JMA', TE_mode='dry',     # for savefig
    start_mode:int=0, end_mode:int=0, contribute:float=0.0,  # for svd
    ):

    """Draw sensitivity area @dry enegy norm"""
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode=prj)
    lon, lat = self.RG.set_coordinate(self.RG.nx, self.RG.ny) 
    x, y = self.MP.coord_change(mapp, lon, lat)
    
    if label_cfmt == 'spread':
      label_cfmt = 'spread_{}hr'.format(ft)
      title_cfmt = 'TE spread [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      save_cfmt = '{}_TE_{}_{}'.format(center,TE_mode,label_cfmt)
    elif label_cfmt == 'adjoint':
      label_cfmt = 'adjoint'
      title_cfmt = ' TE Adjoint sensitivity [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      #title_cfmt = ' NORMALIZE TE Adjoint sensitivity [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      save_cfmt = '{}_TE_{}_{}_{}hr'.format(center,TE_mode,label_cfmt,ft)
    elif label_cfmt == 'SVD':
      label_cfmt = 'svd'
      title_cfmt = ' TE MODE {}-{} contribute:{}% sensitivity [ J/kg ] FT={}hr, INIT={}'.format(start_mode, end_mode,int(contribute),ft,date)
      #title_cfmt = ' NORMALIZE TE MODE sensitivity [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
      save_cfmt = '{}_TE_{}_{}_{}hr_{:02}-{:02}_{}'.format(center,TE_mode,label_cfmt,ft,start_mode,end_mode,date)

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region_lambert(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    self.MP.point_linear(mapp,x,y,0,self.RG.nx-1,0,self.RG.ny-1, color='black', ls='--')
    self.MP.norm_contourf(mapp, x, y, energy_norm, label=label_cfmt)
    self.MP.contour(mapp, x, y, hgt_data[3], elem='850hPa', font_on=1)
    self.MP.title(title_cfmt, fontsize=8)
    self.MP.saving(save_cfmt,'./work/')
    #plt.show()


  def average_norm_driver(self,
    energy_norm, target_region, date,
    *,
    prj='lcc', label_cfmt='spread',       # for title & colorbar
    center='JMA', TE_mode='dry',                    # for savefig
    normalize_set='off', normalize_region=(0,0,0,0) # for normalize
    ):

    """Draw sensitivity area @dry enegy norm"""
    fig, ax = plt.subplots()
    mapp = self.MP.base(projection_mode=prj)
    lon, lat = self.RG.set_coordinate() 
    x, y = self.MP.coord_change(mapp, lon, lat)
    
    title_cfmt = ' AVERAGE TE MODE sensitivity [ J/kg ] FT=Average, INIT={}'.format(date)
    save_cfmt = '{}_TE_{}_{}_{}hr'.format(center,TE_mode,label_cfmt,'ave')

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    #vertifcation region
    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    
    if normalize_set == 'on':
      normal_lat_min_index, normal_lat_max_index, normal_lon_min_index, normal_lon_max_index = \
        self.EN.verification_region(lon,lat,
            area_lat_min=normalize_region[1], area_lat_max=normalize_region[0],
            area_lon_min=normalize_region[2], area_lon_max=normalize_region[3]
        )

      self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')
    
    self.MP.norm_contourf(mapp, x, y, energy_norm, label=label_cfmt)
    self.MP.title(title_cfmt, fontsize=8)
    self.MP.saving(save_cfmt,'./work/')
    #plt.show()

  def each_elem_norm_dry_rish_driver(self, 
    pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_slp, # pertbuation data
    press_levels, target_region, ft, date,
    *, prj='lcc'
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

    #level setting is bug, needs inverse hPa data.
    # UWND
    for index, level in enumerate(press_levels):
      ax = fig.add_subplot(row,column, (1+3*len(press_levels))-4*index)
      mapp = self.MP.base(projection_mode=prj)
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, 0.5*(pertb_uwnd[index]**2) ,label=label_cfmt)
      self.MP.title('UWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH UWND level {}hPa'.format(level))

    # VWND
    for index, level in enumerate(press_levels):
      ax = fig.add_subplot(row,column,(2+3*len(press_levels))-4*index)
      mapp = self.MP.base(projection_mode=prj)
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, 0.5*(pertb_vwnd[index]**2) ,label=label_cfmt)
      self.MP.title('VWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH VWND level {}hPa'.format(level))

    #TMP
    for index, level in enumerate(press_levels[1:]):
      ax = fig.add_subplot(row,column, (3+3*len(press_levels)-4*index-4))
      mapp = self.MP.base(projection_mode=prj)
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, 0.5*((pertb_tmp[index]**2)*(1004/270)) ,label=label_cfmt)
      self.MP.title('TMP [ J/kg ] {}hPa '.format(level))
      print('..... FINISH TMP  level {}hPa'.format(level))

    #SLP
    ax = fig.add_subplot(row,column,16)
    mapp = self.MP.base(projection_mode=prj)
    x, y = self.MP.coord_change(mapp, lon, lat)

    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    self.MP.norm_contourf(mapp, x, y, 0.5*(((pertb_slp**2)/700)*(287*270)) ,label=label_cfmt)
    self.MP.title('SLP [ J/kg ] ')
    print('..... FINISH SLP ')

    plt.suptitle(' Ensemble based Sensitivity Analysis (JMA) Valid time: {}, Target region: JPN'.format(date))
    plt.show()

  def each_elem_norm_dry_tigge_driver(self, 
    pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_ps, # pertbuation data
    target_region, ft, date, 
    *, prj='lcc', label_cfmt='elem_each',  # for title & colorbar
    center='JMA', TE_mode='dry',           # for savefig
    mode:int=0, contribute:float=0.0,      # for svd
    normalize_opt:str='on'
    ):

    size_x, size_y = 20, 20
    row, column = 6, 4 
    press_levels = (200, 300, 500, 700, 850, 1000)
    press_indexs = (  7,   5,   4,   3,   2,    0)

    title_cfmt = 'TE [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
    save_cfmt = '{}_Each_elem_energy_{}_{}'.format(center,TE_mode,label_cfmt)

    fig = plt.figure(figsize = (size_x, size_y))
    lon, lat = self.RG.set_coordinate() 
    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    uwnd = (pertb_uwnd[:]**2)*0.5
    vwnd = (pertb_vwnd[:]**2)*0.5
    tmp  = (self.EN.cp/self.EN.Tr)*((pertb_tmp[:])**2)*0.5
    ps   = (((self.EN.R*self.EN.Tr)/self.EN.Pr)*(pertb_ps**2/self.EN.Pr)*0.5)
    
    if normalize_opt == 'on':
      for _ in press_indexs:
        uwnd[_] = statics_tool.min_max(uwnd[_])
        vwnd[_] = statics_tool.min_max(vwnd[_])
        tmp[_]  = statics_tool.min_max(tmp[_])
      ps = statics_tool.min_max(ps)
      label_cfmt = 'normalize'

    elif normalize_opt == 'off':
      label_cfmt = 'elem_each'

    # UWND
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,1+4*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, uwnd[index] ,label=label_cfmt)
      self.MP.title('UWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH UWND level {:5} hPa'.format(level))

    # VWND
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,2+4*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, vwnd[index] ,label=label_cfmt)
      self.MP.title('VWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH VWND level {:5} hPa'.format(level))

    #TMP
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,3+4*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.norm_contourf(mapp, x, y, tmp[index], label=label_cfmt)
      self.MP.title('TMP [ J/kg ] {}hPa '.format(level))
      print('..... FINISH TMP  level {:5} hPa'.format(level))

    #PS
    ax = fig.add_subplot(row,column,24)
    mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
    x, y = self.MP.coord_change(mapp, lon, lat)

    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index-1,lat_min_index,lat_max_index-1)
    self.MP.norm_contourf(mapp, x, y, ps,label=label_cfmt)
    self.MP.title('PS [ J/kg ] ')
    print('..... FINISH PS   level surface')

    plt.suptitle(' Ensemble based Sensitivity Analysis ({}) Valid time: {}, Target region: JPN'.format(center,date))
    self.MP.saving(save_cfmt,'./work/')
    plt.show()

  def each_elem_norm_humid_tigge_driver(self, 
    pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps, # pertbuation data
    target_region, ft, date, 
    *, prj='lcc', label_cfmt='elem_each',  # for title & colorbar
    center='JMA', TE_mode='humid',         # for savefig
    mode:int=0, contribute:float=0.0,      # for svd
    normalize_set='off', normalize_region=(0,0,0,0) # for normalize
    ):

    size_x, size_y = 25, 20
    row, column = 6, 5 
    press_levels = (200, 300, 500, 700, 850, 1000)
    press_indexs = (  7,   5,   4,   3,   2,    0)

    title_cfmt = 'TE [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
    save_cfmt = '{}_Each_elem_energy_{}_{}_{}hr'.format(center,TE_mode,date,ft)

    fig = plt.figure(figsize = (size_x, size_y))
    lon, lat = self.RG.set_coordinate() 
    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    uwnd = (pertb_uwnd[:]**2)*0.5
    vwnd = (pertb_vwnd[:]**2)*0.5
    tmp  = (self.EN.cp/self.EN.Tr)*((pertb_tmp[:])**2)*0.5
    spfh = (self.EN.wq*(self.EN.Lc**2)*(pertb_spfh[:]**2)/(self.EN.cp*self.EN.Tr))*0.5
    ps   = (((self.EN.R*self.EN.Tr)/self.EN.Pr)*(pertb_ps**2/self.EN.Pr)*0.5)
    
    if normalize_set == 'on':
      normal_lat_min_index, normal_lat_max_index, normal_lon_min_index, normal_lon_max_index = \
        self.EN.verification_region(lon,lat,
          area_lat_min=normalize_region[1], area_lat_max=normalize_region[0],
          area_lon_min=normalize_region[2], area_lon_max=normalize_region[3]
        )
      uwnd = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,uwnd)
      vwnd = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,vwnd)
      tmp  = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,tmp)
      spfh = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,spfh)
      ps   = self.EN.region_normalize_norm(normalize_region,lon,lat,ps)
      label_cfmt = 'normalize'

    elif normalize_set == 'off':
      label_cfmt = 'elem_each'

    # UWND
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,1+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, uwnd[index],label=label_cfmt)
      self.MP.title('UWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH UWND level {:5} hPa'.format(level))

    # VWND
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,2+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, vwnd[index] ,label=label_cfmt)
      self.MP.title('VWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH VWND level {:5} hPa'.format(level))

    #TMP
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,3+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, tmp[index], label=label_cfmt)
      self.MP.title('TMP [ J/kg ] {}hPa '.format(level))
      print('..... FINISH TMP  level {:5} hPa'.format(level))

    #SPFH
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,4+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)
      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, spfh[index],label=label_cfmt)
      self.MP.title('SPFH [ J/kg ] {}hPa '.format(level))
      print('..... FINISH SPFH level {:5} hPa'.format(level))

    #PS
    ax = fig.add_subplot(row,column,30)
    mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
    x, y = self.MP.coord_change(mapp, lon, lat)

    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    self.MP.point_linear(mapp,x,y,normal_lon_min_index,normal_lon_max_index,normal_lat_min_index,normal_lat_max_index, color='black', ls='--')
    self.MP.norm_contourf(mapp, x, y, ps,label=label_cfmt)
    self.MP.title('PS [ J/kg ] ')
    print('..... FINISH PS   level   surface')

    plt.suptitle(' Ensemble based Sensitivity Analysis ({}) Valid time: {}, Target region: JPN'.format(center,date))
    self.MP.saving(save_cfmt,'./work/')
    plt.show()
  

  def each_elem_norm_humid_Nhm_letkf_driver(self, 
    pertb_uwnd, pertb_vwnd, pertb_tmp, pertb_spfh, pertb_ps, # pertbuation data
    target_region, ft, date, 
    *, prj='lcc', label_cfmt='elem_each',  # for title & colorbar
    TE_mode='humid',                       # for savefig
    mode:int=0, contribute:float=0.0,      # for svd
    normalize_set='off', normalize_region=(0,0,0,0) # for normalize
    ):

    size_x, size_y = 25, 20
    row, column = 6, 5 
    press_levels = (200, 300, 500, 700, 850, 1000)
    press_indexs = ( 13,  11,   9,   7,   5,    0)

    title_cfmt = 'TE [ J/kg ] FT={}hr, INIT={}'.format(ft,date)
    save_cfmt = '{}_Each_elem_energy_{}_{}_{}hr'.format('Nhm-letkf',TE_mode,date,ft)

    fig = plt.figure(figsize = (size_x, size_y))
    lon, lat = self.RG.set_coordinate() 
    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      self.EN.verification_region(lon,lat,
          area_lat_min=target_region[1], area_lat_max=target_region[0],
          area_lon_min=target_region[2], area_lon_max=target_region[3]
      )

    uwnd = (pertb_uwnd[:]**2)*0.5
    vwnd = (pertb_vwnd[:]**2)*0.5
    tmp  = (self.EN.cp/self.EN.Tr)*((pertb_tmp[:])**2)*0.5
    spfh = (self.EN.wq*(self.EN.Lc**2)*(pertb_spfh[:]**2)/(self.EN.cp*self.EN.Tr))*0.5
    ps   = (((self.EN.R*self.EN.Tr)/self.EN.Pr)*(pertb_ps**2/self.EN.Pr)*0.5)
    
    if normalize_set == 'on':
      normal_lat_min_index, normal_lat_max_index, normal_lon_min_index, normal_lon_max_index = \
        self.EN.verification_region(lon,lat,
          area_lat_min=normalize_region[1], area_lat_max=normalize_region[0],
          area_lon_min=normalize_region[2], area_lon_max=normalize_region[3]
        )
      uwnd = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,uwnd)
      vwnd = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,vwnd)
      tmp  = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,tmp)
      spfh = self.EN.region_normalize_3d_norm(normalize_region,lon,lat,spfh)
      ps   = self.EN.region_normalize_norm(normalize_region,lon,lat,ps)
      label_cfmt = 'normalize'

    elif normalize_set == 'off':
      label_cfmt = 'elem_each'

    # UWND
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,1+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,0,self.RG.nx-1,0,self.RG.ny-1, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, uwnd[index],label=label_cfmt)
      self.MP.title('UWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH UWND level {:5} hPa'.format(level))

    # VWND
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,2+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,0,self.RG.nx-1,0,self.RG.ny-1, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, vwnd[index] ,label=label_cfmt)
      self.MP.title('VWND [ J/kg ] {}hPa '.format(level))
      print('..... FINISH VWND level {:5} hPa'.format(level))

    #TMP
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,3+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)

      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,0,self.RG.nx-1,0,self.RG.ny-1, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, tmp[index], label=label_cfmt)
      self.MP.title('TMP [ J/kg ] {}hPa '.format(level))
      print('..... FINISH TMP  level {:5} hPa'.format(level))

    #SPFH
    for place, (index, level) in enumerate(zip(press_indexs, press_levels)):
      ax = fig.add_subplot(row,column,4+5*place)
      mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
      x, y = self.MP.coord_change(mapp, lon, lat)
      self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
      self.MP.point_linear(mapp,x,y,0,self.RG.nx-1,0,self.RG.ny-1, color='black', ls='--')
      self.MP.norm_contourf(mapp, x, y, spfh[index],label=label_cfmt)
      self.MP.title('SPFH [ J/kg ] {}hPa '.format(level))
      print('..... FINISH SPFH level {:5} hPa'.format(level))

    #PS
    ax = fig.add_subplot(row,column,30)
    mapp = self.MP.base(projection_mode=prj,label_cfmt='off')
    x, y = self.MP.coord_change(mapp, lon, lat)

    self.MP.point_linear(mapp,x,y,lon_min_index,lon_max_index,lat_min_index,lat_max_index)
    self.MP.point_linear(mapp,x,y,0,self.RG.nx-1,0,self.RG.ny-1, color='black', ls='--')
    self.MP.norm_contourf(mapp, x, y, ps,label=label_cfmt)
    self.MP.title('PS [ J/kg ] ')
    print('..... FINISH PS   level   surface')

    plt.suptitle(' Ensemble based Sensitivity Analysis ({}) Valid time: {}, Target region: JPN'.format(center,date))
    self.MP.saving(save_cfmt,'./work/')
    plt.show()
  

