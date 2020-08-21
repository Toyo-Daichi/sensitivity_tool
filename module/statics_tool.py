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

def normalize(v, axis=-1, order=2):
  """
  平均を0, 標準偏差を1とする正規化
  https://deepage.net/features/numpy-normalize.html
  """
  l2 = np.linalg.norm(v, ord = order, axis=axis, keepdims=True)
  l2[l2==0] = 1
  return v/l2

def min_max(x, axis=None):
  """
  最小値を0, 最大値を1とする正規化
  https://deepage.net/features/numpy-normalize.html
  """
  min = x.min(axis=axis, keepdims=True)
  max = x.max(axis=axis, keepdims=True)
  result = (x-min)/(max-min)
  return result

def weight_average(self, data:np.ndarray, weight_list:np.ndarray):
  weight_average, sum_of_weight = np.average( 
    a = data, axis = 0, weights = weight_list, returned = True
  )
  #print('..... CALCULATE WEIGHT AVE. SUM OF WEIGHT ', sum_of_weight)
  return weight_average
