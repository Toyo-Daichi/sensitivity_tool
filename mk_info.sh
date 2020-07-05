#!/bin/csh -f
#-------------------------------------------------------
# sub          :: mk_info.sh
# purpose      :: ifile devide by element 
# main command :: wgrib2
#
#					  	coded by Toyooka
#						 @2018.7.29
# -------------------------------------------
# 
# this script is only for GSM/JMA
# +++ About GSM(Global Spect Model) data
# 		elem = (wnd, upward flow, tmp, hgt)
# 		12 layer data pre 
# 	   >>> 1000, 925, 850, 700, 600, 500, 400 
# 	   		 300, 250, 200, 150, 100 hPa
#
#-------------------------------------------------------

 set s_yy = 2019
 set e_yy = 2019
 while ( ${s_yy} <= ${e_yy} )

 set s_mm = 10
 set e_mm = 10
 while ( ${s_mm} <= ${e_mm} )
 set m0 = `printf %02d ${s_mm}`

 set s_dd = 20 
 set e_dd = 30

 if ( ${m0} == 01 ) set nday = 31
 if ( ${m0} == 02 ) set nday = 28
 if ( ${m0} == 03 ) set nday = 31
 if ( ${m0} == 04 ) set nday = 30
 if ( ${m0} == 05 ) set nday = 31
 if ( ${m0} == 06 ) set nday = 30
 if ( ${m0} == 07 ) set nday = 31
 if ( ${m0} == 08 ) set nday = 31
 if ( ${m0} == 09 ) set nday = 30
 if ( ${m0} == 10 ) set nday = 31
 if ( ${m0} == 11 ) set nday = 30
 if ( ${m0} == 12 ) set nday = 31
 if ( ${s_yy} % 4 == 0 && ${m0} == 2 ) set nday =29

 while ( ${s_dd} <= ${e_dd} )
 #or
 #while ( ${s_dd} <= ${nday} )
 set d0 = `printf %02d ${s_dd}`

 set s_hh = 0
 set e_hh = 18
 while ( ${s_hh} <= ${e_hh} )
 set h0 = `printf %02d ${s_hh}`

 mkdir -p ./bin/${s_yy}${m0}${d0}
 
 set i_dir = /work3/daichi/puff_model/tra_analysis/data/${s_yy}${m0}${d0}
 set i_file = ${i_dir}/Z__C_RJTD_${s_yy}${m0}${d0}${h0}0000_GSM_GPV_Rgl_FD0000_grib2.bin 

 echo ${i_file}

 set o_file = ./bin/${s_yy}${m0}${d0}

 set I = ( 1000 925 850 700 600 500 400 300 250 200 150 100)
 set i = 1

 while ( ${i} <= 12 )
 
 set C = `printf %4d ${I[${i}]}` 
 set level = `eval echo ${C} mb`

 echo ${level}
 
 wgrib2 -v ${i_file} | grep "UGRD" | grep "${level}" |wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/uwnd_${s_yy}${m0}${d0}${h0}.bin
 wgrib2 -v ${i_file} | grep "VGRD" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/vwnd_${s_yy}${m0}${d0}${h0}.bin
 wgrib2 -v ${i_file} | grep "VVEL" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/vvel_${s_yy}${m0}${d0}${h0}.bin
 wgrib2 -v ${i_file} | grep "TMP" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/tmp_${s_yy}${m0}${d0}${h0}.bin
 wgrib2 -v ${i_file} | grep "HGT" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/hgt_${s_yy}${m0}${d0}${h0}.bin

 @ i = ${i} + 1

 end

 cat ${o_file}/uwnd_${s_yy}${m0}${d0}${h0}.bin ${o_file}/vwnd_${s_yy}${m0}${d0}${h0}.bin ${o_file}/vvel_${s_yy}${m0}${d0}${h0}.bin ${o_file}/tmp_${s_yy}${m0}${d0}${h0}.bin ${o_file}/hgt_${s_yy}${m0}${d0}${h0}.bin >! ${o_file}/${s_yy}${m0}${d0}${h0}.bin

 #rm -f ${i_file}
 rm -f ${o_file}/*_*
 
 @ s_hh = ${s_hh} + 6
 end
 set hh = 0
 @ s_dd = ${s_dd} + 1
 end
 @ s_mm = ${s_mm} + 1
 end
 @ s_yy = ${s_yy} + 1
 end

 exit

