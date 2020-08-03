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
  yyyy, mm, dd, hh, ft = 2005, 9, 2, 00, 72
  dataset = 'WFM'
  data_dir = '/work3/daichi/Data/GSM_EnData'
  indata = data_dir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)

  """Class & data set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  RG = readgpv.ReadGPV(nx,ny,nz,mem)
  full_data = RG.read_gpv(indata)

  for imem in range(mem):
    outdir = data_dir+'/bin/{:03}/'.format(imem)
    os.makedirs(outdir,exist_ok=True)

    with open(outdir+'{}{:02}{:02}{:02}_{:02}hr.grd'.format(yyyy,mm,dd,hh,ft), 'wb') as ofile:
      ofile.write(full_data[:, :, imem, :, :])

    print('...... Output data on ' + '{}{:02}{:02}{:02}_{:02}hr.grd'.format(yyyy,mm,dd,hh,ft) + ' {}'.format(imem))
  
  print('Normal END')
