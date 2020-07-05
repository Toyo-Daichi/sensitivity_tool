# -*- coding: utf-8 -*-
"""
Created from 2020.7.15
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))
import matplotlib.pyplot as plt

#import my_module
import mapping
import readgpv

class Anl_basem:

  def __init__(self):
    pass

  def main_driver(self, data_path):
    """Set parm. """
    surf_elem, elem = RG.data_kind()
    lon, lat = RG.set_coordinate()
    elem_num = len(elem)
    
    """Set data. & Making pertubation"""

    full_data = RG.read_gpv(data_path, elem_num)





    fig, ax = plt.subplots()
    mapp = MP.base(projection_mode='lcc')
    x, y = MP.coord_change(mapp, lon, lat)

    MP.rain_contourf(mapp, x, y, full_data[surf_elem['APCP'],0,24], hr='default')
    MP.contour(mapp, x, y, full_data[surf_elem['PRMSL'],0,24]*0.01)
    #MP.vector(mapp, x, y, surf_data[surf_elem['US']], surf_data[surf_elem['VS']], skip=5)
    #MP.title(self.exp_name + ':: 201807061800, ft:6hr, mem:Mean')		

    plt.show()

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh = 2005, 9, 2, 0
  nx, ny, nz = 144, 37, 4
  ft, mem = 72, 25

  indir = '/work3/daichi/Data/GSM_EnData'
  indata = indir + '/bin/{}{:02}{:02}/'.format(yyyy,mm,dd) + '{}{:02}{:02}{:02}_{:02}hr_{:02}mem.grd'.format(yyyy,mm,dd,hh,ft,mem)

  """Class set"""
  RG = readgpv.ReadGPV(nx,ny,nz,mem)
  MP = mapping.Mapping('JPN')
  #Main_driver
  DR = Anl_basem()

  DR.main_driver(indata)
  print('Normal END')
