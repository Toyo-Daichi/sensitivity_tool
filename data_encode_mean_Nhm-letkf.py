# -*- coding: utf-8 -*-
"""
Created from 2020.11.12
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import readgpv_tigge
import setup
import struct
import subprocess

if __name__ == "__main__":
  """Set basic info. """
  yyyy = '2018'; t_yyyy = '2018' ; e_yyyy = '2018'
  mm   = '07'  ; t_mm   = '07'   ; e_mm   = '07'
  dd   = '04'  ; t_dd   = '05'   ; e_dd   = '05'
  hh   = '12'  ; t_hh   = '12'   ; e_hh   = '12'

  dataset = 'NHM_JPN'
  exp_name, exp_type, ft = 'exp005__JPN_05km_mem050_NEST_FCST', 'MF05km_jpn_M050', '24'

  """Class & parm set """
  init_date, target_date, exp_date = yyyy+mm+dd+hh, t_yyyy+t_mm+t_dd+t_hh, e_yyyy+e_mm+e_dd+e_hh+'00'
  RG = readgpv_nhm.ReadGPV(dataset)

  """Making pretubation data Vertificate TIME"""
  indir = '/data1/da01/Toyo/Data/Nhm-letkf/{}/{}/{}/Ges/'.format(exp_name,exp_type,exp_date)
  uwnd_data, vwnd_data, tmp_data, spfh_data, _, ps_data = RG.data_read_driver(indir,target_date)
  
  print('Normal END')
