# -*- coding: utf-8 -*-
import numpy as np

class Setup:
  def __init__(self, name):
    self.dataset = name

  def set_prm(self):
    if self.dataset is 'WFM':
      nx, ny, nz, mem = 144, 37, 4, 25
    elif self.dataset is 'EPSW':
      nx, ny, nz, mem = 144, 73, 4, 27
    elif self.dataset is 'TIGGE_JMA':
      nx, ny, nz, mem = 288, 145, 8, 27 
    elif self.dataset is 'TIGGE_NCEP':
      nx, ny, nz, mem = 288, 145, 8, 17 
    elif self.dataset is 'TIGGE_ECMWF':
      nx, ny, nz, mem = 288, 145, 8, 50 
    elif self.dataset is 'TIGGE_CMC':
      nx, ny, nz, mem = 288, 145, 8, 20 
    elif self.dataset is 'TIGGE_UKMO':
      nx, ny, nz, mem = 288, 145, 8, 17 

      print(nx,ny,nz,mem)

    return nx, ny, nz, mem

  def set_pressure_levels(self):
    if self.dataset is 'WFM':
      press_levels = np.array([1000.0, 850.0, 500.0, 300.0])
    elif self.dataset is 'EPSW':
      press_levels = np.array([1000.0, 850.0, 500.0, 300.0])
    elif 'TIGGE' in self.dataset:
      press_levels = np.array([1000.0, 925.0, 850.0, 700.0, 500.0, 300.0, 250.0, 200.0])
    return press_levels

""" simple package"""
def save_list_ndarray(data:list, indir:str, name:str):
  print('..... SAVE NPY FILE ' + indir + name)
  data = np.array(data, dtype='float32')
  np.save(indir+name,data)

# *** memory check
#from memory_profiler import profile
#@profile()
#def _SVD_error_check(self,array):
#  U_array, sigma_array, V_array = svd(array,full_matrices=True)
#  return array

#def print_varsize():
#  import types
#  print('')
#  print("{}{: >15}{}{: >10}{}".format('|','Variable Name','|','  Size','|'))
#  print(" -------------------------- ")
#  for k, v in globals().items():
#    if hasattr(v, 'size') and not k.startswith('_') and not isinstance(v,types.ModuleType):
#      print("{}{: >15}{}{: >10}{}".format('|',k,'|',str(v.size),'|'))
#    elif hasattr(v, '__len__') and not k.startswith('_') and not isinstance(v,types.ModuleType):
#      print("{}{: >15}{}{: >10}{}".format('|',k,'|',str(len(v)),'|'))




