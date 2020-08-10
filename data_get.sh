#!/bin/csh -f
set datapath = '/work3/daichi/Data/GSM_EnData/grib/'

# set date
set s_yy = 2003; set e_yy = 2003
set s_mm = 8   ; set e_mm = 8
set s_dd = 5   ; set e_dd = 5
set s_hh = 12  ; set e_hh = 12 

while ( ${s_yy} <= ${e_yy} ) 
  while ( ${s_mm} <= ${e_mm} )
    set m0 = `printf %02d ${s_mm}`
    while ( ${s_dd} <= ${e_dd} )
      set d0 = `printf %02d ${s_dd}`

      cd ${datapath}
      mkdir -p ./${s_yy}${m0}${d0}
      set i_dir = ${datapath}/${s_yy}${m0}${d0}
 
      sleep 5

      while ( ${s_hh} <= ${e_hh} )
        set h0 = `printf %02d ${s_hh}`
        set accum_day = ${s_yy}${m0}${d0}${h0}
        
        if ( ${accum_day} <= 2006030100 ) then
          wget http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/${s_yy}/${m0}/${d0}/WFM12XSFC -P ${i_dir}
          wget http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/${s_yy}/${m0}/${d0}/WFM12XPLL -P ${i_dir}
          wget http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/${s_yy}/${m0}/${d0}/WFM12XPLM -P ${i_dir}
          wget http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/${s_yy}/${m0}/${d0}/WFM12XPLH -P ${i_dir}

        else if ( ${accum_day} <= 2020032300 ) then
          wget http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/${s_yy}/${m0}/${d0}/Z__C_RJTD_${s_yy}${m0}${d0}${h0}0000_EPSW_GPV_Rgl_FD00-08_grib2.bin -P ${i_dir}

        else 
          echo "CHECK your DATE or wget Z__C_RJTD_yyyyMMddhhmmss_EPSG_GPV_Rgl_Gll1p25deg_FD0000-0100_grib2.bin"

        endif

        @ s_hh = ${s_hh} + 6
      end
    
      @ s_dd = ${s_dd} + 1
      cd ${datapath}
    end

  @ s_mm = ${s_mm} + 1
  end

@ s_yy = ${s_yy} + 1
end

exit
