# -*- coding: utf-8 -*-
"""
Created from 2020.8.3
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import readgpv
import setup

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = 2003, 8, 5, 12, 72 
  date = '{:04}{:02}{:02}{:02}'.format(yyyy,mm,dd,hh)
  dataset = 'WFM' # 'WFM' or 'EPSW'
  
  """Class & data set """
  EN = readgpv.Energy_NORM(dataset)
  data_path = './work/Z_array_{}{:02}{:02}{:02}_{:02}hr.npy'.format(yyyy,mm,dd,hh,ft)
  data_array = np.load(data_path)

  print('..... @ MAKE SVD DECOMPOSITION @')
  u_array,sig_array,v_array = EN.singular_decomposion(data_array, try_num=1000)
  
  setup.save_list_ndarray(u_array,'./work/','U_array_{}_{}hr'.format(date,ft))
  setup.save_list_ndarray(sig_array,'./work/','SIG_array_{}_{}hr'.format(date,ft))
  setup.save_list_ndarray(v_array,'./work/','V_array_{}_{}hr'.format(date,ft))
  
  print('Normal END')
