# -*- coding: utf-8 -*-
"""
Created from 2020.8.18
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
  yyyy, mm, dd, hh, ft = 2018, 7, 4, 12, 00
  date = '{:04}{:02}{:02}{:02}'.format(yyyy,mm,dd,hh)
  center = 'JMA'
  dataset = 'TIGGE_' + center 
  set_endian = 'little' 
  xyz_var_list, surf_list = ('UGRD', 'VGRD', 'HGT', 'TMP', 'SPFH'), ('PS')
  xyz_var_num, surf_var_num = 5, 1
  make_var = 0 # 0(make each var output) or 1(only full data)
  
  """Class & data set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  dims = nx*ny*nz
  data_dir = '/work3/daichi/Data/TIGGE/'
  indata = data_dir + center + '/{}{:02}{:02}{:02}/'.format(yyyy,mm,dd,hh) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)

  RG = readgpv_tigge.ReadGPV(dataset,date,ft)
  data = RG.set_gpv(indata,endian=set_endian)
  xyz_data, surf_data = data[:dims*xyz_var_num].reshape(xyz_var_num,nz,mem,ny,nx), data[dims*xyz_var_num:].reshape(mem,ny,nx)

  cfmt = 'f'*(dims)

  for imem in range(mem):
    outdir = data_dir+'/bin/{}{:02}{:02}/{:03}/'.format(yyyy,mm,dd,imem+1)
    os.makedirs(outdir,exist_ok=True)

    for ivar, name in enumerate(xyz_var_list):
      with open(outdir+name+'.grd','wb') as ofile:
        _data =np.ravel(xyz_data[ivar,:,imem,:,:])
      grd = struct.pack(cfmt,*_data)
      ofile.write(grd)
    
    with open(outdir+'PS'+'.grd','wb') as ofile:
      _data =np.ravel(surf_data[imem,:,:])
      grd = struct.pack(cfmt,*_data)

    try: 
      command = ["bash","./module/cat.sh",outdir,'{:04}{:02}{:02}{:02}'.format(yyyy,mm,dd,hh),'{:02}'.format(ft),'TIGGE']
      res = subprocess.call(command)
      print('...... Output data on '+'{}{:02}{:02}{:02}_{:02}hr.grd'.format(yyyy,mm,dd,hh,ft)+' MEM::{:03}'.format(imem+1))

    except:
      print('suprocess.check_call() failed')
    
    if (make_var == 1):
      command = ["rm", outdir+'/UGRD.grd', outdir+'/VGRD.grd', outdir+'/HGT.grd', outdir+'/TMP.grd', outdir+'/SPFH.grd', outdir+'/PS.grd']
      res = subprocess.call(command)
      
  if (make_var == 1):
    command = ["rm", indata]
    res = subprocess.call(command)
  
  print('Normal END')
