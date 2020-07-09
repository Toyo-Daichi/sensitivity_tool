# -*- coding: utf-8 -*-
"""
Created from 2020.7.15
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

#import my_module
import mapping
import readgpv
import setup

class Anl_basem:

  def __init__(self):
    pass

  def main_driver(self, data_path:str):
    """
    Args:
      data_path(str) : 週間アンサンブルデータのPATH 
    Note:
      full_data(np.ndarray)  : grib形式をバイナリー形式に自作したデータセット
      constitution -> [要素, 高度面, ensemble_mem:0=ctrl_run, 緯度, 経度] 
      *** 気温/高度は, 地表面では積算降水量/海面更生気圧となっているので注意
      pertb_elem(np.ndarray) : コントロールランから各メンバーを引いた摂動データセット
      constitution -> [ensemble_mem -1(cntl分), 高度, 緯度, 経度]
    """

    """Set parm. """
    surf_elem, elem = RG.data_kind()
    lon, lat = RG.set_coordinate()
    weight_lat = RG.weight_latitude(lat)
    press_levels = ST.set_pressure_levels()
    elem_num = len(elem)
    
    """Set data. & Making pertubation"""
    full_data = RG.read_gpv(data_path, elem_num)
    pertb_uwnd = [[] for _ in range(RG.ensemble_size-1)]
    pertb_vwnd = [[] for _ in range(RG.ensemble_size-1)]
    pertb_tmp = [[] for _ in range(RG.ensemble_size-1)]
    pertb_slp = [[] for _ in range(RG.ensemble_size-1)]

    for imem in range(1, RG.ensemble_size):
      pertb_uwnd[imem-1] = RG.calc_prime(full_data[elem['UGRD'],1:,0], full_data[elem['UGRD'],1:,imem])
      pertb_vwnd[imem-1] = RG.calc_prime(full_data[elem['VGRD'],1:,0], full_data[elem['VGRD'],1:,imem])
      pertb_tmp[imem-1] = RG.calc_prime(full_data[elem['TMP'],1:,0], full_data[elem['TMP'],1:,imem])
      pertb_slp[imem-1] = RG.calc_prime(np.log(full_data[surf_elem['PRMSL'],0,0]*0.01), np.log(full_data[surf_elem['PRMSL'],0,imem]*0.01))
      # latitude weight
      pertb_uwnd[imem-1] = pertb_uwnd[imem-1][:,:]*weight_lat
      pertb_vwnd[imem-1] = pertb_vwnd[imem-1][:,:]*weight_lat
      pertb_tmp[imem-1]  = pertb_tmp[imem-1][:,:]*weight_lat*np.sqrt(EN.cp/EN.Tr)
      pertb_slp[imem-1]  = pertb_slp[imem-1][:,:]*weight_lat*np.sqrt((EN.R*EN.Tr)/EN.Pr)

    """Pertubation np.array Update"""
    pertb_pertb_uwnd = np.array(pertb_uwnd)
    pertb_pertb_vwnd = np.array(pertb_uwnd)
    pertb_pertb_tmp = np.array(pertb_tmp)
    pertb_pertb_slp = np.array(pertb_slp)

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

    #"""Description func. """
    #fig, ax = plt.subplots()
    #mapp = MP.base(projection_mode='lcc')
    #x, y = MP.coord_change(mapp, lon, lat)

    #MP.rain_contourf(mapp, x, y, full_data[surf_elem['APCP'],0,0], hr='default')
    #MP.contour(mapp, x, y, full_data[elem['HGT'], 2, 0], elem='500hPa')
    #MP.norm_contourf(mapp, x, y, dry_energy_norm[0])
    #MP.contour(mapp, x, y, full_data[surf_elem['PRMSL'],0,1]*0.01)
    #MP.rain_contourf(mapp, x, y, full_data[surf_elem['PRMSL']*0.01,1,24], hr='default')
    #for _ in range(1, RG.ensemble_size):
    #MP.vector(mapp, x, y, surf_data[surf_elem['US']], surf_data[surf_elem['VS']], skip=5)
    #MP.title('Test run, ft:72hr')
    #plt.show()

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

  if ft == 'anl':
    indir = '/work3/daichi/Data/GSM_EnData'
    indata = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_anlhr_{:02}mem.grd'.format(yyyy,mm,dd,hh,mem)
  elif ft != 'anl':
    indir = '/work3/daichi/Data/GSM_EnData'
    indata = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)

  DR = Anl_basem()
  DR.main_driver(indata)
  print('Normal END')
