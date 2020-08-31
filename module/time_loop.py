# -*- coding: utf-8 -*-

def making_timelist(
    yyyy:int, mm:int, dd:int, hh:int,
    during:int, file_num:int
  ) -> list:
  """数時間ごとの時間リストを作成
  Args:
      yyyy (int): 初期時刻の年
      mm (int)  : 初期時刻の月
      dd (int)  : 初期時刻の日
      hh (int)  : 初期時刻の時間
      during (int)  : 次の出力時間までのインターバル(hr)
      file_num (int): 出力するタイムリストの総数
  Returns:
      list: 任意の数のタイムリスト
  Note:
      Created on 2019.12.09 Ver1.0
      Created on 2020.07.19 Ver1.1
  """
  time_list = []
  num_time_list = len(time_list)
  
  month_thirtyone = [ 1, 3, 5, 7, 8, 10, 12 ]
  month_thirty    = [ 4, 6, 9, 11 ]
  month_twntynine = [ 2 ]
  
  while num_time_list < file_num:
    time_list.append(str(yyyy) + str('%02d' % mm) + str('%02d' % dd) + str('%02d' % hh) + '00')
    hh = hh + during

    if hh >= 24:
      hh, dd = 0, dd + 1 
      if mm in month_thirty:
        if dd >= 30:
          mm, dd = mm + 1, 1
      elif mm in month_thirtyone:
        if dd >= 31:
          mm, dd = mm + 1, 1
      elif mm in month_twntynine:
        if yyyy % 4 == 0:
          if dd >= 28:
            mm, dd = mm + 1, 1
        else:
          if dd >=29:
            mm, dd = mm + 1, 1
    num_time_list += 1

  return time_list
