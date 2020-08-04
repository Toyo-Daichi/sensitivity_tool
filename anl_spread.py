# -*- coding: utf-8 -*-
"""
Created from 2020.8.4
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../module'))
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

#import my_module
import mapping
import readgpv

class Anl_SPREAD:
    def __init__(self):
      pass
      
    def main_driver(self, data):
      """Draw sensitivity area @dry enegy norm"""
      fig, ax = plt.subplots()
      mapp = MP.base(projection_mode='lcc')
      x, y = MP.coord_change(mapp, lon, lat)
    
      MP.contour(mapp, x, y,  data[0,2,:,:],elem='500hPa')
      plt.show()

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = '2005', '09', '02', '00', '72'
  date = yyyy+mm+dd+hh
  dataset = 'WFM' # 'WFM' or 'EPSW'

  """Class & parm set """
  DR = Anl_SPREAD()
  RG = readgpv.ReadGPV(dataset,date,ft)
  MP = mapping.Mapping('NH')

  indir = '/work3/daichi/Data/GSM_EnData/bin/'
  uwnd_data, vwnd_data, hgt_data, tmp_data, slp_data, rain_data = RG.data_read_driver(indir)

  print('Normal END')