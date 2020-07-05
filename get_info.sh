#!/bin/csh -f
#--------------------------------------------
#  sub       :: get_info.sh
#  purpose   :: ypu can get JMA GSM_data
#
#  you decide start date & end date.
#    >(s_yy, s_mm, s_dd  & e_yy, e_mm, e_dd) 
#							
#							coded by Toyooka
#							 @2018.7.29
#--------------------------------------------

 set s_yy = 2019
 set e_yy = 2019
 while ( ${s_yy} <= ${e_yy} ) 

 set s_mm = 10 
 set e_mm = 10
 while ( ${s_mm} <= ${e_mm} )
	 	
 set m0 = ${s_mm}
 if( ${s_mm} <= 9 ) set m0 = 0${s_mm}

 set s_dd = 1 
 set e_dd = 30
 while ( ${s_dd} <= ${e_dd} )
 set d0 = ${s_dd}
 if( ${s_dd} <= 9 ) set d0 = 0${s_dd}

 cd data
 mkdir -p ${s_yy}${m0}${d0}
 set i_dir = /work3/daichi/puff_model/tra_analysis/data/${s_yy}${m0}${d0}
 
 sleep 5

 set o_time = 0000
 
 set s_hh = 0 
 set e_hh = 18
 while ( ${s_hh} <= ${e_hh} )
 set h0 = ${s_hh}
 if( ${s_hh} <= 9 ) set h0 = 0${s_hh}

 wget http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/${s_yy}/${m0}/${d0}/Z__C_RJTD_${s_yy}${m0}${d0}${h0}0000_GSM_GPV_Rgl_FD${o_time}_grib2.bin -P ${i_dir}


 @ s_hh = ${s_hh} + 6
 end
 @ s_dd = ${s_dd} + 1
 cd /work3/daichi/puff_model/tra_analysis/ 
 end
 @ s_mm = ${s_mm} + 1
 end
 @ s_yy = ${s_yy} + 1
 end

exit
