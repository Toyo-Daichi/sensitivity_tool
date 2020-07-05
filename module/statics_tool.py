#-*- coding: utf-8 -*-

import numpy as np
from sklearn.metrics import *

def mae(true:list, prediction:list) -> list:
  return mean_absolute_error(true, prediction)

def rmse(true:list, prediction:list) -> np.ndarray:
  return np.sqrt(mean_squared_error(true, prediction))

def mse(true:list, prediction:list) -> list:
  return mean_squared_error(true, prediction)

def variance(x) -> float:
  return np.var(x) 

def covariance(x,y) -> float:
  _mtx = _covariance_mtx(x,y)
  return _mtx[0,1]

def _covariance_mtx(x,y):
  return np.cov(x,y)

def correlation_coefficient(sig_xy:float, sig_x:float, sig_y:float) -> float:
  return sig_xy / (sig_x*sig_y)

def gaussian_func(x, *, amp:float=1.0, ave:float=0, std:float=1.0):
  """ガウス関数
  Args:
    amp(float): 振幅
    ave(float): 平均
    std(float): 標準偏差
  Returns:
    gaussian_function(float)
  """
  return amp*np.exp(-1*((x - ave)/2*std)**2)