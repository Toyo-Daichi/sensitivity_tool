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
import struct
import subprocess

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = 2019, 10, 11, 12, 72 
  date = '{:04}{:02}{:02}{:02}'.format(yyyy,mm,dd,hh)
  dataset = 'EPSW' # 'WFM' or 'EPSW'
  var_list = ('UGRD', 'VGRD', 'HGT', 'TMP') #level=surf, HGT, TMP -> PRMSL, APCP
  make_var = 1 # 0(make each var output) or 1(only full data)
  
  """Class & data set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  data_dir = '/work3/daichi/Data/GSM_EnData'
  indata = data_dir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)

  if (dataset is 'WFM'):
    set_endian = 'default' # little
  elif (dataset is 'EPSW'):
    set_endian = 'big'

  RG = readgpv.ReadGPV(dataset,date,ft)
  full_data = RG.set_gpv(indata,len(var_list), endian=set_endian)
  cfmt = 'f'*(nx*ny*nz)

  for imem in range(mem):
    outdir = data_dir+'/bin/{}{:02}{:02}/{:03}/'.format(yyyy,mm,dd,imem+1)
    os.makedirs(outdir,exist_ok=True)

    for ivar, name in enumerate(var_list):
      with open(outdir+name+'.grd','wb') as ofile:
        _data =np.ravel(full_data[ivar,:,imem,:,:])
        grd = struct.pack(cfmt,*_data)
        ofile.write(grd)
    
    try: 
      command = ["bash","./module/cat.sh",outdir,'{:04}{:02}{:02}{:02}'.format(yyyy,mm,dd,hh),'{:02}'.format(ft)]
      res = subprocess.call(command)
      print('...... Output data on '+'{}{:02}{:02}{:02}_{:02}hr.grd'.format(yyyy,mm,dd,hh,ft)+' MEM::{:03}'.format(imem+1))

    except:
      print('suprocess.check_call() failed')
    
    if (make_var == 1):
      command = ["rm", outdir+'/UGRD.grd', outdir+'/VGRD.grd', outdir+'HGT.grd', outdir+'TMP.grd']
      res = subprocess.call(command)
      
  if (make_var == 1):
    command = ["rm", indata]
    res = subprocess.call(command)
  
  print('Normal END')
