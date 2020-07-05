# -*- coding: utf-8 -*-
import numpy as np

class ReadGPV:
  def __init__(self, nx, ny, nz, mem):
    self.nx = nx
    self.ny = ny
    self.nz = nz
    self.ensemble_size = mem

  def data_kind(self, mode='default'):
    if mode is 'default':
      surf_elem = { 'UGRD': 0, 'VGRD': 1, 'PRMSL':2, 'APCP':3 }
      elem      = { 'UGRD': 0, 'VGRD': 1, 'HGT':2, 'TMP':3 } 
    return surf_elem, elem

  def set_coordinate(self, dx=2.5, dy=2.5):
    lon, lat = [], []
    for ix in range(self.nx):
      lon += [ float('{:.2f}'.format(dx*ix)) ]
    for iy in range(self.ny):
      lat += [ float('{:.2f}'.format(dy*iy)) ]
    X, Y = np.meshgrid(lon, lat)
    return X, Y
    
  def read_gpv(self, gpv_file, elem):
    return self._open_gpv(gpv_file).reshape(elem, self.nz, self.ensemble_size, self.ny, self.nx)

  def _open_gpv(self, gpv_file, *, mode='default'):
    print('...... Preparating for {}'.format(gpv_file))
    with open(gpv_file, 'rb') as ifile:
      if mode == 'default':
        data = np.fromfile(ifile, dtype='<f', sep = '')
      return data

