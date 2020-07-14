# -*- coding: utf-8 -*-
"""
Created from 2020.7.15
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../module'))
import numpy as np
import matplotlib.pyplot as plt

#import my_module
import mapping
import readgpv
import setup

class Anl_basem:

  def __init__(self):
    pass

  def sensitivity_driver(self,  vertifi_data_path:str):
    """
    Args:
      vertifi_data_path(str): 週間アンサンブルデータのPATH(検証時刻)
    Note:
      init_full_data(np.ndarray) : grib形式をバイナリー形式に自作したデータセット(初期時刻)
      vertifi_full_data(np.ndarray) : grib形式をバイナリー形式に自作したデータセット(検証時刻)
      constitution -> [要素, 高度面, ensemble_mem:0=ctrl_run, 緯度, 経度] 
      *** 気温/高度は, 地表面では積算降水量/海面更生気圧となっているので注意
      pertb_elem(np.ndarray) : コントロールランから各メンバーを引いた摂動データセット
      constitution -> [ensemble_mem -1(cntl分), 高度, 緯度, 経度]
      ens_rate_list(np.ndarray) : アンサンブルメンバーに割り振る割合リスト
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
    vertifi_full_data = RG.read_gpv(init_data_path, elem_num)
    #vertifi_full_data = RG.read_gpv(vertifi_data_path, elem_num)
    pertb_uwnd = np.zeros((RG.ensemble_size-1, RG.nz-RG.surf, RG.ny, RG.nx))
    pertb_vwnd = np.zeros((RG.ensemble_size-1, RG.nz-RG.surf, RG.ny, RG.nx))
    pertb_tmp  = np.zeros((RG.ensemble_size-1, RG.nz-RG.surf, RG.ny, RG.nx))
    pertb_slp  = np.zeros((RG.ensemble_size-1, RG.ny, RG.nx))

    for imem in range(1, RG.ensemble_size):
      pertb_uwnd[imem-1, :, :, :] = RG.calc_prime(vertifi_full_data[elem['UGRD'],1:,0], vertifi_full_data[elem['UGRD'],1:,imem])
      pertb_vwnd[imem-1, :, :, :] = RG.calc_prime(vertifi_full_data[elem['VGRD'],1:,0], vertifi_full_data[elem['VGRD'],1:,imem])
      pertb_tmp[imem-1, :, :, :] = RG.calc_prime(vertifi_full_data[elem['TMP'],1:,0], vertifi_full_data[elem['TMP'],1:,imem])
      pertb_slp[imem-1, :, :] = RG.calc_prime(np.log(vertifi_full_data[surf_elem['PRMSL'],0,0]*0.01), np.log(vertifi_full_data[surf_elem['PRMSL'],0,imem]*0.01))
      # latitude weight
      pertb_uwnd[imem-1] = pertb_uwnd[imem-1]*weight_lat
      pertb_vwnd[imem-1] = pertb_vwnd[imem-1]*weight_lat
      pertb_tmp[imem-1]  = pertb_tmp[imem-1]*weight_lat*np.sqrt(EN.cp/EN.Tr)
      pertb_slp[imem-1]  = pertb_slp[imem-1]*weight_lat*np.sqrt((EN.R*EN.Tr))/EN.Pr

    """Calc. dry enegy norm"""
    # check each norm
    imem = 17 s
    dry_energy_norm = EN.dry_energy_norm(
      pertb_uwnd[imem], pertb_vwnd[imem],
      pertb_tmp[imem], pertb_slp[imem],
      press_levels 
      )

    lat_min_index, lat_max_index, lon_min_index, lon_max_index = \
      EN.verification_region(
        lon, lat,
        area_lat_min =50, area_lat_max =20,
        area_lon_min =120, area_lon_max =150
      )

    print('..... Check Vertification area Norm SUM {}'.format(
      np.sum(dry_energy_norm[lat_min_index:lat_max_index,lon_min_index:lon_max_index]))
    )

    """Draw sensitivity area @dry enegy norm"""
    fig, ax = plt.subplots()
    mapp = MP.base(projection_mode='lcc')
    x, y = MP.coord_change(mapp, lon, lat)
    
    MP.norm_contourf(mapp, x, y, dry_energy_norm*10)
    MP.contour(mapp, x, y, init_full_data[elem['HGT'], 2, 0], elem='500hPa')
    MP.title('Test run, ft:72hr')
    plt.show()


if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh = 2005, 9, 2, 00 
  init, ft = 'anl', 72 
  dataset = 'WFM'

  """Class & parm set """
  ST = setup.Setup(dataset)
  nx, ny, nz, mem = ST.set_prm()
  RG = readgpv.ReadGPV(nx,ny,nz,mem)
  EN = readgpv.Energy_norm(nx,ny) 
  MP = mapping.Mapping('NH')
  DR = Anl_basem()

  indir = '/work3/daichi/Data/GSM_EnData'
  vertifi_data = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)
  ensm_rate_list = indir + '/rate/' + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.npy'.format(yyyy,mm,dd,hh,ft,mem)

  DR.sensitivity_driver(vertifi_data, ensm_rate_list)
  print('Normal END')
