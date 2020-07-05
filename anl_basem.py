# -*- coding: utf-8 -*-
"""
Created from 2020.7.15
@author: Toyo_Daichi
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), './module'))

#import my_module
import mapping
import readgpv

class Anl_basem:

  def __init__(self):
    pass

  def main_driver(self, data_path):
    mp = mapping.Mapping('WJPN')
    elem, surf_elem = rg.data_kind()
    data, surf_data = rg.open_grd_letkf(data_path, len(elem), len(surf_elem))
    lon, lat = rg.data_coord(self.indir)

    fig, ax = plt.subplots()
    mapp = mp.base(projection_mode='lcc')
    x, y = mp.coord_change(mapp, lon, lat)

    mp.rain_contourf(mapp, x, y, surf_data[surf_elem['RAIN']], hr='default')
    mp.contour(mapp, x, y, data[elem['P']][5]*0.01)
    mp.vector(mapp, x, y, surf_data[surf_elem['US']], surf_data[surf_elem['VS']], skip=5)

    mp.title(self.exp_name + ':: 201807061800, ft:6hr, mem:Mean')		
    #mp.saving('201807061800ft6hr', self.indir + 'Nhm-letkf/fig/' + self.exp_name + '/')

    plt.show()

if __name__ == "__main__":
  """Set basic info. """
  yyyy, mm, dd, hh = 2005, 9, 2, 0
  nx, ny, nz = 144, 37, 4
  ft, ensemble_size = 72, 25

  indir = '/work3/daichi/Data/GSM_EnData/'
  indata = indir + '/bin/{0}{1}{2}/'.format(yyyy,mm,dd) + '{0}{1}{2}{3}_{4}hr_{5}mem.grd'.format(yyyy,mm,dd,hh,ft,ensemble_size)

  """Class set"""
  RG = readgpv.ReadGPV(nx,ny,nz)
  MP = mapping.Mapping('NH')
  #Main_driver
  DR = Anl_basem()

  data = RG.read_gpv(indata)

  print('Normal END')
