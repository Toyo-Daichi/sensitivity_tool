# -*- coding: utf-8 -*-
import numpy as np

class Setup:
  def __init__(self, name):
    self.dataset = name

  def set_prm(self):
    if self.dataset is 'WFM':
      nx, ny, nz, mem = 144, 37, 4, 25
    elif self.dataset is 'EPSW':
      nx, ny, nz = 144, 73, 4, 27
    return nx, ny, nz, mem

  def set_pressure_levels(self):
    if self.dataset is 'WFM':
      press_levels = np.array([850.0, 500.0, 300.0])
    elif self.dataset is 'EPSW':
      press_levels = np.array([850.0, 500.0, 300.0])
    return press_levels

""" simple package"""

def save_list_ndarray(data:list, indir:str, name:str):
  data = np.array(data, dtype='float32')
  np.save(indir+name,data)
