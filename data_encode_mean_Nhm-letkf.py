# -*- coding: utf-8 -*-
"""
Created from 2020.11.12
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import readgpv_nhm
import setup
import struct
import subprocess
from time_loop import making_timelist

# I/O
def data_read_driver(exp_path:str, target_date:str, *, endian='big'):
  uwnd, vwnd, tmp, spfh, hgt, wwnd, vor, slp, ps, rain = RG.init_mem_array()
  for index, imem in enumerate(range(RG.mem)):
    data_path = exp_path + '/{:03}/fcst_p.grd/{}.grd'.format(imem+1, target_date)
    full_data = RG._open_gpv(data_path, endian=endian)
    _uwnd, _vwnd, _tmp, _spfh, _hgt, _wwnd, _vor, _slp, _ps, _rain = RG.full_data_cut(full_data)
    del full_data
    uwnd[index,:,:,:] = _uwnd[:,:,:]
    vwnd[index,:,:,:] = _vwnd[:,:,:]
    tmp[index,:,:,:]  = _tmp[:,:,:]
    spfh[index,:,:,:] = _spfh[:,:,:]
    hgt[index,:,:,:]  = _hgt[:,:,:]
    wwnd[index,:,:,:] = _wwnd[:,:,:]
    vor[index,:,:,:]  = _vor[:,:,:]
    slp[index,:,:]    = _slp[:,:]
    ps[index,:,:]     = _ps[:,:]
    rain[index,:,:]   = _rain[:,:]

  return uwnd, vwnd, tmp, spfh, hgt, wwnd, vor, slp, ps, rain

def write_grd(gpv_file,data,cfmt):
  with open(gpv_file, 'wb') as ofile:
    _data = np.ravel(data)
    grd = struct.pack(cfmt,*_data)
    ofile.write(grd)

# Calculate
def mean(array):
  return np.mean(array, axis=0)

if __name__ == "__main__":
  """Set basic info. """
  yyyy = 2018; e_yyyy = 2018
  mm   = 7   ; e_mm   = 7
  dd   = 5   ; e_dd   = 6
  hh   = 12  ; e_hh   = 12

  dt, file_num = 1, 25
  # 5km
  dataset = 'NHM_JPN'
  #exp_name, exp_type = 'exp005__JPN_05km_mem050_NEST_FCST', 'MF05km_jpn_M050'
  #exp_name, exp_type = 'exp007__JPN_05km_mem050_NEST_FCST_LACK', 'MF05km_jpn_M050'
  #exp_name, exp_type = 'exp007__JPN_05km_mem050_NEST_FCST_LACK_2cycle', 'MF05km_jpn_M050'
  #exp_name, exp_type = 'exp007__JPN_05km_mem050_NEST_FCST_LACK_3cycle', 'MF05km_jpn_M050'
  exp_name, exp_type = 'exp007__JPN_05km_mem050_NEST_FCST_LACK_4cycle', 'MF05km_jpn_M050'
  
  # 2km
  #dataset = 'NHM_WJPN'
  #exp_name, exp_type = 'exp008__JPN_02km_mem025_W_NEST_FCST', 'MF2km_nest_M025'
  #exp_name, exp_type = 'exp009__JPN_02km_mem025_W_NEST_FCST_LACK', 'MF2km_nest_M025'

  """Class & parm set """
  RG = readgpv_nhm.ReadGPV(dataset)
  init_date, exp_date = '{}{:02}{:02}{:02}'.format(yyyy,mm,dd,hh), '{}{:02}{:02}{:02}00'.format(e_yyyy,e_mm,e_dd,e_hh)
  dateset = making_timelist(yyyy,mm,dd,hh,dt,file_num)

  """Making pretubation data Vertificate TIME"""
  indir = '/data1/da01/Toyo/Data/Nhm-letkf/{}/{}/{}/Ges/'.format(exp_name,exp_type,exp_date)
  elems = ('UGRD', 'VGRD', 'TMP', 'SPFH', 'HGT', 'WWND', 'VOR', 'SLP', 'PS', 'RAIN')
  outdir = indir + '000' + '/fcst_p.grd/'
  os.makedirs(outdir,exist_ok=True)
  for target_date in dateset:
    print('>>> @ DATE: {} @'.format(target_date))
    uwnd_mems, vwnd_mems, tmp_mems, spfh_mems, hgt_mems, wwnd_mems, vor_mems, slp_mems, ps_mems, rain_mems =\
    data_read_driver(indir,target_date)
    
    # make ensemble mean
    uwnd_mean = mean(uwnd_mems);  vwnd_mean = mean(vwnd_mems) 
    tmp_mean  = mean(tmp_mems) ;  spfh_mean = mean(spfh_mems) 
    hgt_mean  = mean(hgt_mems) ;  wwnd_mean = mean(wwnd_mems) 
    vor_mean  = mean(vor_mems) ;  slp_mean  = mean(slp_mems) 
    ps_mean   = mean(ps_mems)  ;  rain_mean = mean(rain_mems)  
    
    for index, var in enumerate(list((uwnd_mean,vwnd_mean,tmp_mean,spfh_mean,hgt_mean,wwnd_mean,vor_mean,slp_mean,ps_mean,rain_mean))):
      elem = elems[index]
      if index < 7: #3d full set
        nz, ny, nx = np.shape(var)
      else: #2dset
        ny, nx = np.shape(var); nz=1
      dims = nx*ny*nz
      cfmt = '>' + 'f'*dims
      write_grd(outdir+'{}.grd'.format(elem),var,cfmt)

    try: 
      command = ["bash","./module/cat.sh",outdir,target_date,'00','LETKF']
      res = subprocess.call(command)
      print('...... Output data on '+'{}.grd'.format(target_date))
      
    except:
      print('suprocess.check_call() failed')
      quit()
  
    # tidy up
    for elem in elems:
      os.remove(outdir+'/{}.grd'.format(elem))

  print('Normal END')
