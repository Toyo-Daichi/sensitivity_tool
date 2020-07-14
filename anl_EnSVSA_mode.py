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
      *** 気温/高度は, 地表面では積算降水量/海面更生気圧となっているので注意,
          緯度の重み付けをしてしまっているので降水量を書くときは注意(* これで書かない).
      pertb_elem(np.ndarray) : コントロールランから各メンバーを引いた摂動データセット
      constitution -> [高度*緯度*経度, ensemble_mem -1(cntl分)]
    """

    """Set parm. """
    surf_elem, elem = RG.data_kind()
    lon, lat = RG.set_coordinate()
    weight_lat = RG.weight_latitude(lat)
    press_levels = ST.set_pressure_levels()
    elem_num, dims = len(elem), RG.nx*RG.ny

    
    """Set data. & Making pertubation"""
    full_data = RG.read_gpv(data_path, elem_num)
    full_data = full_data[:,:,:]*weight_lat #latitude weight

    eigen_mtx  = np.zeros((10*RG.nx*RG.ny, RG.ensemble_size-1)) 
    pertb_uwnd = np.zeros((3*RG.nx*RG.ny, RG.ensemble_size-1))
    pertb_vwnd = np.zeros((3*RG.nx*RG.ny, RG.ensemble_size-1))
    pertb_tmp  = np.zeros((3*RG.nx*RG.ny, RG.ensemble_size-1))
    pertb_slp  = np.zeros((RG.nx*RG.ny, RG.ensemble_size-1))

    for imem in range(RG.ensemble_size-1):
      pertb_uwnd[0*dims:1*dims, imem] = (RG.calc_prime(full_data[elem['UGRD'],1,0,:,:].reshape(-1), full_data[elem['UGRD'],1,imem,:,:].reshape(-1)))
      pertb_uwnd[1*dims:2*dims, imem] = (RG.calc_prime(full_data[elem['UGRD'],2,0,:,:].reshape(-1), full_data[elem['UGRD'],2,imem,:,:].reshape(-1)))
      pertb_uwnd[2*dims:3*dims, imem] = (RG.calc_prime(full_data[elem['UGRD'],3,0,:,:].reshape(-1), full_data[elem['UGRD'],3,imem,:,:].reshape(-1)))
      pertb_vwnd[0*dims:1*dims, imem] = (RG.calc_prime(full_data[elem['VGRD'],1,0,:,:].reshape(-1), full_data[elem['VGRD'],1,imem,:,:].reshape(-1)))
      pertb_vwnd[1*dims:2*dims, imem] = (RG.calc_prime(full_data[elem['VGRD'],2,0,:,:].reshape(-1), full_data[elem['VGRD'],2,imem,:,:].reshape(-1)))
      pertb_vwnd[2*dims:3*dims, imem] = (RG.calc_prime(full_data[elem['VGRD'],3,0,:,:].reshape(-1), full_data[elem['VGRD'],3,imem,:,:].reshape(-1)))
      pertb_tmp[0*dims:1*dims, imem] = (RG.calc_prime(full_data[elem['TMP'],1,0,:,:].reshape(-1), full_data[elem['TMP'],1,imem,:,:].reshape(-1)))
      pertb_tmp[1*dims:2*dims, imem] = (RG.calc_prime(full_data[elem['TMP'],2,0,:,:].reshape(-1), full_data[elem['TMP'],2,imem,:,:].reshape(-1)))
      pertb_tmp[2*dims:3*dims, imem] = (RG.calc_prime(full_data[elem['TMP'],3,0,:,:].reshape(-1), full_data[elem['TMP'],3,imem,:,:].reshape(-1)))
      pertb_slp[0*dims:1*dims, imem] = (RG.calc_prime(full_data[surf_elem['PRMSL'],1,0,:,:].reshape(-1), full_data[surf_elem['PRMSL'],3,imem,:,:].reshape(-1)))

    """Making singular value of decomposition of Z"""
    eigen_mtx[0*dims:3*dims,:] = pertb_uwnd 
    eigen_mtx[3*dims:6*dims,:] = pertb_vwnd 
    eigen_mtx[6*dims:9*dims,:] = pertb_tmp
    eigen_mtx[9*dims:10*dims,:] = pertb_slp 

    print('...... start SVD CALCULATE')
    U, sigma, V = np.linalg.svd(eigen_mtx, full_matrices=True)
    E_eigen = (U@sigma)/np.sqrt(RG.ensemble_size-1)
    
    np.save('eigen_martix', E_eigen)

      #pertb_uwnd[imem-1] = pertb_uwnd[imem-1]*weight_lat
      #pertb_vwnd[imem-1] = pertb_vwnd[imem-1]*weight_lat
      #pertb_tmp[imem-1]  = pertb_tmp[imem-1]*weight_lat*np.sqrt(EN.cp/EN.Tr)
      #pertb_slp[imem-1]  = pertb_slp[imem-1]*weight_lat*np.sqrt(EN.R*EN.Tr)/EN.Pr

    #eigen_mtx = np.zeros(((elem_num-1)*(RG.nz-1)*RG.nx*RG.ny+RG.nx*RG.ny, RG.ensemble_size-1))
    
    #for imem in range(RG.ensemble_size):
      #eigen_mtx[ :(RG.nz-1)*RG.nx*RG.ny, imem] = pertb_uwnd[imem, :]
      #eigen_mtx[  (RG.nz-1)*RG.nx*RG.ny:2*(RG.nz-1)*RG.nx*RG.ny, imem] = pertb_vwnd[imem,:]
      #eigen_mtx[2*(RG.nz-1)*RG.nx*RG.ny:3*(RG.nz-1)*RG.nx*RG.ny, imem] = pertb_tmp[imem,:]
      #eigen_mtx[3*(RG.nz-1)*RG.nx*RG.ny:4*RG.nx*RG.ny, imem] = pertb_slp[imem,:]


if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = 2005, 9, 2, 00, 'anl' 
  dataset = 'WFM'

  """Class & parm set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  RG = readgpv.ReadGPV(nx,ny,nz,mem)
  EN = readgpv.Energy_norm(nx,ny) 
  MP = mapping.Mapping('JPN')
  DR = Anl_EnSVSA()

  indir = '/work3/daichi/Data/GSM_EnData'
  indata = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)
  DR.En_singular_vector_sensitivity_driver(indata)
  
  print('Normal END')
