# -*- coding: utf-8 -*-
"""
Created from 2020.8.3
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))

import numpy as np
import readgpv_rish
import mapping 
import matplotlib.pyplot as plt


if __name__ == "__main__":
  
  with open('/work3/daichi/Data/GSM_EnData/grib/20180704/ctl.grd', 'rb') as ifile:
    data_ctl = np.fromfile(ifile, dtype='>f', sep = '').reshape(73,144)

  with open('/work3/daichi/Data/GSM_EnData/grib/20180704/10.grd', 'rb') as ifile:
    data_plus = np.fromfile(ifile, dtype='>f', sep = '').reshape(73,144)

  with open('/work3/daichi/Data/GSM_EnData/grib/20180704/10-.grd', 'rb') as ifile:
    data_minus = np.fromfile(ifile, dtype='>f', sep = '').reshape(73,144)

  yyyy, mm, dd, hh, init, ft = '2018', '07', '04', '12', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'EPSW' # 'WFM' or 'EPSW'
  RG = readgpv_rish.ReadGPV(dataset,date,ft)
  lon, lat = RG.set_coordinate()

  MP = mapping.Mapping('CNH')

  fig, ax = plt.subplots()
  mapp = MP.base(projection_mode='lcc')
  lon, lat = RG.set_coordinate() 
  x, y = MP.coord_change(mapp, lon, lat)

  #MP.diff_contourf(mapp, x, y, data_plus-data_ctl, elem='normalize')
  MP.diff_contourf(mapp, x, y, data_minus-data_ctl, elem='normalize')

  plt.show()

