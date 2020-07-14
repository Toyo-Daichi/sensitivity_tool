# -*- coding: utf-8 -*-
"""
Created from 2020.7.15
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np

#import my_module
import mapping
import readgpv
import setup

class Anl_EnSVSA:

  def __init__(self):
    pass

  def En_singular_vector_sensitivity_driver(self, data_path:str):
    """
    Args:
      data_path(str) : 週間アンサンブルデータのPATH
    Note:
      full_data(np.ndarray)  : grib形式をバイナリー形式に自作したデータセット
      constitution -> [要素, 高度面, ensemble_mem:0=ctrl_run, 緯度, 経度] 
      *** 気温/高度は, 地表面では積算降水量/海面更生気圧となっているので注意
      pertb_elem(np.ndarray) : コントロールランから各メンバーを引いた摂動データセット
      constitution -> [高度*緯度*経度, ensemble_mem -1(cntl分)]
    """

    """Set parm. """
    surf_elem, elem = RG.data_kind()
    lon, lat = RG.set_coordinate()
    weight_lat = RG.weight_latitude(lat)
    press_levels = ST.set_pressure_levels()
    elem_num = len(elem)
    
    """Set data. & Making pertubation"""
    full_data = RG.read_gpv(data_path, elem_num)
    pertb_uwnd = np.zeros(((RG.nz-RG.surf)*RG.ny*RG.nx, RG.ensemble_size-1))
    pertb_vwnd = np.zeros(((RG.nz-RG.surf)*RG.ny*RG.nx, RG.ensemble_size-1))
    pertb_tmp  = np.zeros(((RG.nz-RG.surf)*RG.ny*RG.nx, RG.ensemble_size-1))
    pertb_slp  = np.zeros((RG.ny*RG.nx, RG.ensemble_size-1))

    for imem in range(1, RG.ensemble_size):
      pertb_uwnd[:, imem] = RG.calc_prime(full_data[elem['UGRD'],1:,0].reshape(-1), full_data[elem['UGRD'],1:,imem].reshape(-1))
    print(eigen_mtx)






      #pertb_vwnd[imem-1, :] = RG.calc_prime(full_data[elem['VGRD'],1:,0].reshape(-1), full_data[elem['VGRD'],1:,imem].reshape(-1))
      #pertb_tmp [imem-1, :] = RG.calc_prime(full_data[elem['TMP'],1:,0].reshape(-1), full_data[elem['TMP'],1:,imem].reshape(-1))
      #pertb_slp [imem-1, :] = RG.calc_prime(full_data[surf_elem['PRMSL'],0,0].reshape(-1)*0.01, full_data[surf_elem['PRMSL'],0,imem].reshape(-1)*0.01)
      ## latitude weight
      #pertb_uwnd[imem-1] = pertb_uwnd[imem-1]*weight_lat
      #pertb_vwnd[imem-1] = pertb_vwnd[imem-1]*weight_lat
      #pertb_tmp[imem-1]  = pertb_tmp[imem-1]*weight_lat*np.sqrt(EN.cp/EN.Tr)
      #pertb_slp[imem-1]  = pertb_slp[imem-1]*weight_lat*np.sqrt(EN.R*EN.Tr)/EN.Pr

    #"""Making singular value of decomposition of Z"""
    #eigen_mtx = np.zeros(((elem_num-1)*(RG.nz-1)*RG.nx*RG.ny+RG.nx*RG.ny, RG.ensemble_size-1))
    
    #for imem in range(RG.ensemble_size):
      #eigen_mtx[ :(RG.nz-1)*RG.nx*RG.ny, imem] = pertb_uwnd[imem, :]
      #eigen_mtx[  (RG.nz-1)*RG.nx*RG.ny:2*(RG.nz-1)*RG.nx*RG.ny, imem] = pertb_vwnd[imem,:]
      #eigen_mtx[2*(RG.nz-1)*RG.nx*RG.ny:3*(RG.nz-1)*RG.nx*RG.ny, imem] = pertb_tmp[imem,:]
      #eigen_mtx[3*(RG.nz-1)*RG.nx*RG.ny:4*RG.nx*RG.ny, imem] = pertb_slp[imem,:]


if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = 2005, 9, 2, 00, 72
  dataset = 'WFM'

  """Class & parm set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  RG = readgpv.ReadGPV(nx,ny,nz,mem)
  EN = readgpv.Energy_norm(nx,ny) 
  MP = mapping.Mapping('JPN')
  DR = Anl_EnSVSA()

  indir = '/work3/daichi/Data/GSM_EnData'
  indata = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)
  DR.En_singular_vector_sensitivity_driver(indata)
  
  print('Normal END')
