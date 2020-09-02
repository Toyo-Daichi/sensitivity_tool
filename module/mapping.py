# -*- coding: utf-8 -*-

""" If you use local env., please check out. """
#import os, conda
#conda_file_dir = conda.__file__
#conda_dir = conda_file_dir.split('lib')[0]
#proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
#os.environ["PROJ_LIB"] = proj_lib

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import seaborn as sns

class Mapping:

  def __init__(self, area):
    if area == 'JPN':
      self.area = area
      self.lon_min, self.lon_max = 120, 155
      self.lat_min, self.lat_max = 17, 50
      self.lat_0, self.lon_0 = 35, 140
    
    elif area == 'EJPN':
      self.area = area
      self.lon_min, self.lon_max = 120, 155
      self.lat_min, self.lat_max = 17, 50
      self.lat_0, self.lon_0 = 35, 140

    elif area == 'WJPN':
      self.area = area
      self.lon_min, self.lon_max = 129, 141
      self.lat_min, self.lat_max = 29, 38
      self.lat_0, self.lon_0 = 35, 140

    elif area == 'WNH':
      self.area = area
      self.lon_min, self.lon_max = 80, 175
      self.lat_min, self.lat_max = 0, 65 
      self.lat_0, self.lon_0 = 35, 140

    elif area == 'CNH':
      self.area = area
      self.lon_min, self.lon_max = 102, 203
      self.lat_min, self.lat_max = 10, 60 
      self.lat_0, self.lon_0 = 35, 135

    elif area == 'ALL':
      self.area = area
      self.lat_0, self.lon_0 = 0, 180

  def base(self, *, projection_mode='lcc', label_cfmt='on'):
    
    if projection_mode is 'lcc':
      mapping = Basemap( 
        projection=projection_mode,
        resolution="i", 
        lat_0=self.lat_0, lon_0=self.lon_0, fix_aspect=(1,1),
        llcrnrlat=self.lat_min, urcrnrlat=self.lat_max, 
        llcrnrlon=self.lon_min, urcrnrlon=self.lon_max
      )

    elif projection_mode is 'cyl':
      mapping = Basemap( 
        projection=projection_mode,
        resolution="l", 
        lat_0=self.lat_0, lon_0=self.lon_0,
       )
    
    mapping.drawcoastlines(color='black', linewidth=0.5)

    if 'NH' in self.area:
      if label_cfmt == 'on':
        mapping.drawmeridians(np.arange(0, 360, 15),  labels=[False, False, False, True], fontsize='small', color='gray', linewidth=0.5)
        mapping.drawparallels(np.arange(-90, 90, 15), labels=[True, False, False, False], fontsize='small', color='gray', linewidth=0.5)
      elif label_cfmt == 'off':
        mapping.drawmeridians(np.arange(0, 360, 15),  labels=[False, False, False, False], fontsize='small', color='gray', linewidth=0.5)
        mapping.drawparallels(np.arange(-90, 90, 15), labels=[False, False, False, False], fontsize='small', color='gray', linewidth=0.5)

    elif 'JPN' in self.area:
      mapping.drawmeridians(np.arange(0, 360, 5),  labels=[False, False, False, True], fontsize='small', color='gray', linewidth=0.5)
      mapping.drawparallels(np.arange(-90, 90, 5), labels=[True, False, False, False], fontsize='small', color='gray', linewidth=0.5)

    elif self.area == 'ALL':
      mapping.drawmeridians(np.arange(0, 360, 40),  labels=[False, False, False, True], fontsize='small', color='gray', linewidth=0.5)
      mapping.drawparallels(np.arange(-90, 90, 25), labels=[False, False, False, True], fontsize='small', color='gray', linewidth=0.5)

    return mapping

  def coord_change(self, basemap, x, y):
    """
    Changing map coordinate: basemap is base.mapping
    """
    return basemap(x, y)

  def rain_contourf(self, basemap, x, y, data, *, hr='default'):
    if hr == 'default':
      rain_levels = [0.4, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    elif hr == '1hr':
      rain_levels = [0.1, 0.5, 1.0, 2.0, 4.0, 10.0, 30.0]
    elif hr == '3hr':
      rain_levels = [1.0, 5.0, 10.0, 20.0, 30.0, 50.0, 70.0]
    elif hr == '24hr':
      rain_levels = [10.0, 80.0, 120.0, 160.0, 240.0, 320.0, 400.0]

    rain_colors = ['#FFFFFF', '#00FFFF', '#000080', '#228B22', '#FFFF00', '#FF8000', '#FF0000', '#FF00FF']

    cmap = plt.contourf(x, y, data, rain_levels, colors=rain_colors, extend='both')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    cbar.set_label('[mm/3hr]', size=6)

  def corr_contourf(self, basemap, x, y, data):
    levels = np.arange(-1.0, 1.1, 0.1) 
    _cmap = sns.diverging_palette(220, 10, as_cmap=True) 
    cmap = plt.contourf(x, y, data, levels, cmap=_cmap, extend='both')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    
  def curl_contourf(self, basemap, x, y, data):
    levels = np.arange(-7.0, 7.0, 0.5) 
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    cmap = plt.contourf(x, y, data, levels, cmap=cm.coolwarm, extend='both')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')

  def qflux_contourf(self, basemap, x, y, data):
    levels = np.arange(10.0, 20.0, 0.5) 
    label = '[ kg/kg x m/s x $ 10^{2} $ ]'
    cmap = plt.contourf(x, y, data, levels, cmap=cm.PuBu, extend='max')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    cbar.set_label(label, size=8)

  def diff_contourf(self, basemap, x, y, data, *, elem='default'):
    if elem == 'default':
      levels = np.arange(-7.0, 7.0, 0.5) 
      label = '[ none ]'
    elif elem == 'normalize':
      levels = np.arange(-1.0, 1.05, 0.05)
      label = '[ none ]'
    elif elem == 'HGT':
      levels = np.arange(-5.0, 5.01, 0.1)
      label = '[ m ]'
    elif elem == 'small':
      levels = np.arange(-0.2, 0.21, 0.01)
      label = '[ none ]'
    elif elem == 'rain':
      levels = np.arange(-15.0, 16.0, 1.0)
      label = '[mm]'
    elif elem == 'curl':
      levels = np.arange(-3.0, 4.0, 0.5)
      label = '[$ s^{-1} $ x $ 10^{5} $]'
    elif elem == 'qflux':
      levels = np.arange(-5.0, 6.0, 0.5)
      label = '[ kg/kg x m/s x $ 10^{2} $ ]'

    _cmap = sns.diverging_palette(220, 10, as_cmap=True) 
    cmap = plt.contourf(x, y, data, levels, cmap=_cmap, extend='both')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    cbar.set_label(label, size=8)

  def norm_contourf(self, basemap, x, y, data, *, label='normalize', cbar_on=0):
    #normalize
    if label == 'normalize':
      #levels = [0.025, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0]
      levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]

    elif label == 'normalize_scope':
      levels = [0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.200]
      #levels = [0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
      #levels = [0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.250]
    
    #simple
    elif label == 'adjoint':
      #levels = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
      #levels = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
      levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]
      #levels = [1.00, 1.25, 1.5, 2.0, 2.5, 3.0, 5.0]
      #levels = [3.0, 4.0, 5.0, 7.5, 10.0, 12.5, 15.0]
      #levels = [7.5, 10.0, 15.0, 30.0, 50.0, 70.0, 100.0]
    
    elif label == 'svd':
      #levels = [0.025, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0]
      levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0]
      #levels = [0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0]
      #levels = [7.5, 10.0, 15.0, 30.0, 50.0, 70.0, 100.0]
      #levels = [20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 300.0]
      #levels = [850.0, 1000.0, 1150.0, 1300.0, 1500.0, 2000.0, 3000.0]
      #levels = [7000.0, 8000.0, 9000.0, 10000.0, 12000.0, 15000.0, 20000.0]
      #levels = [10000.0, 12500.0, 15000.0, 17500.0, 20000.0, 30000.0, 50000.0]

    elif label == 'spread_00hr' or label == 'spread_12hr' or label == 'spread_24hr':
      #levels = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
      #levels = [5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0]
      levels = [250.0, 300.0, 350.0, 400.0, 500.0, 750.0, 1000.0]
      #levels = [850.0, 1000.0, 1150.0, 1300.0, 1500.0, 2000.0, 3000.0]
    
    elif label == 'spread_48hr' or label == 'spread_72hr':
      #levels = [7.5, 10.0, 15.0, 30.0, 50.0, 70.0, 100.0]
      #levels = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0]
      levels = [7000.0, 8000.0, 9000.0, 10000.0, 12000.0, 15000.0, 20000.0]

    elif label == 'elem_each':
      #levels = [7.5, 10.0, 15.0, 30.0, 50.0, 70.0, 100.0]
      levels = [7.5, 10.0, 15.0, 30.0, 50.0, 70.0, 100.0]

    elif label == 'elem_each_spfh_ps':
      levels = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    
    colors = ['#FFFFFF', '#00FFFF', '#000080', '#228B22', '#FFFF00', '#FF8000', '#FF0000', '#FF00FF']
    cmap = plt.contourf(x, y, data, levels, colors=colors, extend='both')
    
    if cbar_on == 0:
      cbar = basemap.colorbar(cmap, 'right', size='2.5%')
      cbar.set_label('[ J/kg ]', size=8)
    elif cbar_on == 1:
      return cmap

  def hatch_contourf(self, basemap, x, y, data, *, levels=[0.0, 5.0], extend='max'):
    basemap.contourf(x, y, data, levels, alpha=0.0, hatches=['xx'], lw=0.5, colors='None', extend=extend)
    
  def contour(self, basemap, x, y, data, *, elem='default', colors="black", linestyles='-', linewidths=0.5, font_on=0):
    if elem == 'default':
      clevs = np.arange(995.0, 1020.0, 5.0) 
    elif elem == 'diff':
      clevs = np.arange(-10.0, 10.0, 1.0) 
    elif elem == 'norm':
      clevs = np.arange(0.0, 5.0, 0.1) 
    elif elem == '500hPa':
      clevs = np.arange(5000.0, 6000.0, 100) 
    elif elem == '850hPa':
      clevs = np.arange(500.0, 2500.0, 50) 
    elif elem == 'slp':
      clevs = np.arange(970.0, 1020.0, 2.5) 

    contour = basemap.contour(x, y, data, clevs, colors=colors, linestyles=linestyles, linewidths=linewidths)
    if font_on == 1:
      contour.clabel(fmt='%1.1f', fontsize=8)
  
  def vector(self, basemap, x, y, u, v, *, skip=10, scale=1.0e-4):
    arr_skip = (slice(None,None,skip),slice(None,None,skip))
    basemap.quiver(x[arr_skip], y[arr_skip], u[arr_skip], v[arr_skip], angles='xy',scale_units='xy',scale=scale, color='gray', alpha=0.7)

  def barb(self, basemap, x, y, u, v):
    pass

  def point_linear(self, basemap, lon, lat, x_min, x_max, y_min, y_max, *, color='red', ls='-', lw=2.0):
    basemap.plot(lon[y_min:y_max+1,x_max], lat[y_min:y_max+1,x_max], color=color, ls=ls, lw=lw)
    basemap.plot(lon[y_min:y_max+1,x_min], lat[y_min:y_max+1,x_min], color=color, ls=ls, lw=lw)
    basemap.plot(lon[y_min,x_min:x_max+1], lat[y_min,x_min:x_max+1], color=color, ls=ls, lw=lw)
    basemap.plot(lon[y_max,x_min:x_max+1], lat[y_max,x_min:x_max+1], color=color, ls=ls, lw=lw)

  def point(self, basemap, lon, lat, x, y, *, marker="*", markersize=30):
    basemap.plot(
       lon[y,x], lat[y,x], marker=marker, markersize=markersize
    )

  def text(self, basemap, x, y, string, *, size=20, color="black"):
    plt.text(x, y, string, size=size, color=color)

  def title(self, txt, *, fontsize=10):
    plt.title(txt, loc='left', fontsize=fontsize)

  def saving(self, title, outpath):
    plt.savefig(outpath + title + '.png')
    print('...... Saving fig :: ', outpath + title + '.png')

