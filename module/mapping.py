# -*- coding: utf-8 -*-

""" If you use local env., please check out. """
#
#import os, conda
#conda_file_dir = conda.__file__
#conda_dir = conda_file_dir.split('lib')[0]
#proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
#os.environ["PROJ_LIB"] = proj_lib

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap

class Mapping:

  def __init__(self, area):
    if area == 'JPN':
      self.lon_min, self.lon_max = 120, 155
      self.lat_min, self.lat_max = 17, 50
      self.lat_0, self.lon_0 = 35, 140
    
    elif area == 'EJPN':
      self.lon_min, self.lon_max = 120, 155
      self.lat_min, self.lat_max = 17, 50
      self.lat_0, self.lon_0 = 35, 140

    elif area == 'WJPN':
      self.lon_min, self.lon_max = 129, 141
      self.lat_min, self.lat_max = 29, 38
      self.lat_0, self.lon_0 = 35, 140

    elif area == 'NH':
      self.lon_min, self.lon_max = 120, 55
      self.lat_min, self.lat_max = 25, 40
      self.lat_0, self.lon_0 = 35, 140

  def base(self, *, projection_mode='lcc'):
    
    mapping = Basemap( 
      projection=projection_mode,
      resolution="i", 
      lat_0=self.lat_0, lon_0=self.lon_0, fix_aspect=(1,1),
      llcrnrlat=self.lat_min, urcrnrlat=self.lat_max, 
      llcrnrlon=self.lon_min, urcrnrlon=self.lon_max
    )
    
    mapping.drawcoastlines(color='black', linewidth=0.5)
    mapping.drawmeridians(np.arange(0, 360, 5),  labels=[False, True, False, True], fontsize='small', color='gray', linewidth=0.5)
    mapping.drawparallels(np.arange(-90, 90, 5), labels=[True, False, False, True], fontsize='small', color='gray', linewidth=0.5)
        
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
    cmap = plt.contourf(x, y, data, levels, cmap=cm.coolwarm)
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    
  def curl_contourf(self, basemap, x, y, data):
    levels = np.arange(-7.0, 7.0, 0.5) 
    cmap = plt.contourf(x, y, data, levels, cmap=cm.coolwarm)
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
      label = '[none]'
    elif elem == 'rain':
      levels = np.arange(-15.0, 16.0, 1.0)
      label = '[mm]'
    elif elem == 'curl':
      levels = np.arange(-3.0, 4.0, 0.5)
      label = '[$ s^{-1} $ x $ 10^{5} $]'
    elif elem == 'qflux':
      levels = np.arange(-5.0, 6.0, 0.5)
      label = '[ kg/kg x m/s x $ 10^{2} $ ]'

    cmap = plt.contourf(x, y, data, levels, cmap=cm.coolwarm, extend='both')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    cbar.set_label(label, size=8)

  def norm_contourf(self, basemap, x, y, data):
    levels = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    colors = ['#FFFFFF', '#00FFFF', '#000080', '#228B22', '#FFFF00', '#FF8000', '#FF0000', '#FF00FF']

    cmap = plt.contourf(x, y, data, levels, colors=colors, extend='both')
    cbar = basemap.colorbar(cmap, 'right', size='2.5%')
    cbar.set_label('[J/kg]', size=6)

  def contour(self, basemap, x, y, data, *, elem='default', colors="black", linestyles='-', linewidths=0.5):
    if elem == 'default':
      clevs = np.arange(995.0, 1020.0, 5.0) 
    elif elem == 'diff':
      clevs = np.arange(-10.0, 10.0, 1.0) 
    elif elem == 'norm':
      clevs = np.arange(0.0, 5.0, 0.1) 
    elif elem == '500hPa':
      clevs = np.arange(500.0, 6000.0, 100) 

    contour = basemap.contour(x, y, data, clevs, colors=colors, linestyles=linestyles, linewidths=linewidths)
    contour.clabel(fmt='%1.1f', fontsize=8)
  
  def vector(self, basemap, x, y, u, v, *, skip=10, scale=1.0e-4):
    arr_skip = (slice(None,None,skip),slice(None,None,skip))
    basemap.quiver(x[arr_skip], y[arr_skip], u[arr_skip], v[arr_skip], angles='xy',scale_units='xy',scale=scale, color='gray', alpha=0.7)

  def barb(self, basemap, x, y, u, v):
    pass

  def point(self, basemap, x, y, *, marker="o", markersize=10):
    map.plot(x, y, marker=marker, markersize=markersize)

  def text(self, basemap, x, y, string, *, size=20, color="black"):
    plt.text(x, y, string, size=size, color=color)

  def title(self, txt):
    plt.title(txt, loc='left', fontsize=10)

  def saving(self, title, outpath):
    plt.savefig(outpath + title + '.png')
    print('...... Saving fig :: ', outpath + title + '.png')

