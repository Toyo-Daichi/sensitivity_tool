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

class Anl_basem:

  def __init__(self):
    pass

  def En_ajoint_sensitivity_driver(self, data_path:str):
    """
    Args:
      data_path(str) : 週間アンサンブルデータのPATH
    Returns:
      ensemble_rate_list(list) : アンサンブルメンバーに割り振る割合
    Note:
      full_data(np.ndarray)  : grib形式をバイナリー形式に自作したデータセット
      constitution -> [要素, 高度面, ensemble_mem:0=ctrl_run, 緯度, 経度] 
      *** 気温/高度は, 地表面では積算降水量/海面更生気圧となっているので注意
      pertb_elem(np.ndarray) : コントロールランから各メンバーを引いた摂動データセット
      constitution -> [ensemble_mem -1(cntl分), 高度, 緯度, 経度]
      dry_energy_norm(list, np.ndarray) : 各メンバーから求めた乾燥エネルギーノルム
      constitution -> [ensemble_mem -1(cntl分), np.ndarray[緯度, 経度]]
    """

    """Set parm. """
    surf_elem, elem = RG.data_kind()
    lon, lat = RG.set_coordinate()
    weight_lat = RG.weight_latitude(lat)
    press_levels = ST.set_pressure_levels()
    elem_num = len(elem)
    
    """Set data. & Making pertubation"""
    full_data = RG.read_gpv(data_path, elem_num)
    pertb_uwnd = np.zeros((RG.emsemble_size-1, RG.nz, RG.ny, RG.nx))
    pertb_vwnd = np.zeros((RG.emsemble_size-1, RG.nz, RG.ny, RG.nx))
    pertb_tmp  = np.zeros((RG.emsemble_size-1, RG.nz, RG.ny, RG.nx))
    pertb_slp  = np.zeros((RG.emsemble_size-1, RG.ny, RG.nx))

    for imem in range(1, RG.ensemble_size):
      pertb_uwnd[imem-1, :, :, :] = RG.calc_prime(full_data[elem['UGRD'],1:,0], full_data[elem['UGRD'],1:,imem])
      pertb_vwnd[imem-1, :, :, :] = RG.calc_prime(full_data[elem['VGRD'],1:,0], full_data[elem['VGRD'],1:,imem])
      pertb_tmp[imem-1, :, :, :] = RG.calc_prime(full_data[elem['TMP'],1:,0], full_data[elem['TMP'],1:,imem])
      pertb_slp[imem-1, :, :] = RG.calc_prime(np.log(full_data[surf_elem['PRMSL'],0,0]*0.01), np.log(full_data[surf_elem['PRMSL'],0,imem]*0.01))
      # latitude weight
      pertb_uwnd[imem-1] = pertb_uwnd[imem-1]*weight_lat
      pertb_vwnd[imem-1] = pertb_vwnd[imem-1]*weight_lat
      pertb_tmp[imem-1]  = pertb_tmp[imem-1]*weight_lat*np.sqrt(EN.cp/EN.Tr)
      pertb_slp[imem-1]  = pertb_slp[imem-1]*weight_lat*np.sqrt((EN.R*EN.Tr)/EN.Pr)
   
    """Calc. dry enegy norm"""
    dry_energy_norm = [[] for _ in range(RG.ensemble_size)]
    for imem in range(RG.ensemble_size-1):
      dry_energy_norm[imem] = EN.dry_energy_norm(
        pertb_uwnd[imem], pertb_vwnd[imem],
        pertb_tmp[imem],  pertb_slp[imem],
        press_levels 
      )
    
    """Calc. rate of Each ensemble member"""
    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(
        lon, lat,
        area_lat_min =50, area_lat_max =20,
        area_lon_min =120, area_lon_max =150
    )

    vertification_region_norm_list = []
    for imem in range(RG.ensemble_size-1):
      vertification_region_norm_list.append(np.sum(dry_energy_norm[imem][lat_min_index:lat_max_index,lon_min_index:lon_max_index]))

    ensemble_rate_list = []
    for imem in range(RG.ensemble_size-1):
      ensemble_rate_list.append((vertification_region_norm_list[imem]/sum(vertification_region_norm_list))*100)

    return ensemble_rate_list

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh, ft = 2005, 9, 2, 0, 72
  dataset = 'WFM'

  """Class & parm set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  RG = readgpv.ReadGPV(nx,ny,nz,mem)
  EN = readgpv.Energy_norm(nx,ny) 
  MP = mapping.Mapping('JPN')
  DR = Anl_basem()

  indir = '/work3/daichi/Data/GSM_EnData'
  indata = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)

  ensemble_rate_list =  DR.En_ajoint_sensitivity_driver(indata)
  DR.En_ajoint_sensitivity_driver_v1(indata, ensemble_rate_list)

  print('Normal END')
