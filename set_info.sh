#!/bin/csh -f
set datapath = '/work3/daichi/Data/GSM_EnData/'

# set date
set s_yy = 2019; set e_yy = 2019
set s_mm = 10  ; set e_mm = 10
set s_dd = 1   ; set e_dd = 30
set s_hh = 0   ; set e_hh = 18

# set your target info.
set ft = 72 

while ( ${s_yy} <= ${e_yy} )

  while ( ${s_mm} <= ${e_mm} )
  set m0 = `printf %02d ${s_mm}`
  
    while ( ${s_dd} <= ${e_dd} )
    set d0 = `printf %02d ${s_dd}`

      while ( ${s_hh} <= ${e_hh} )
      set h0 = `printf %02d ${s_hh}`

      cd ${datapath}
      mkdir -p ./${s_yy}${m0}${d0}
      set i_dir = ${datapath}/grib/${s_yy}${m0}${d0}
      set o_dir = ${datapath}/bin/

      set accum_day = ${s_yy}${m0}${d0}${h0}
      if ( ${accum_day} <= 2006030100 ) then
        set i_surf = ${i_dir}/WFM12XSFC
        set i_850  = ${i_dir}/WFM12XPLL
        set i_500  = ${i_dir}/WFM12XPLM
        set i_300  = ${i_dir}/WFM12XPLH

        set i_list = ( ${i_surf} ${i_850} ${i_500} ${300} )
        set il     = 1

        while ( ${il} <= 4 )
        set ifile = ${i_list[${il}]}

        wgrib -v ${i_file} | grep "UGRD" | wgrib ${i_file} -i -no_header -append -ieee ${o_file}/uwnd_${s_yy}${m0}${d0}${h0}.grd
        wgrib -v ${i_file} | grep "VGRD" | wgrib ${i_file} -i -no_header -append -ieee ${o_file}/vwnd_${s_yy}${m0}${d0}${h0}.grd
        wgrib -v ${i_file} | grep "TMP"  | wgrib ${i_file} -i -no_header -append -ieee ${o_file}/tmp_${s_yy}${m0}${d0}${h0}.grd

        if ( ${il} == 1 ) wgrib -v ${i_file} | grep "sfc" | wgrib ${i_file} -i -no_header -append -ieee ${o_file}/hgt_${s_yy}${m0}${d0}${h0}.grd

        if ( ${il} != 1 ) wgrib -v ${i_file} | grep "HGT" | wgrib ${i_file} -i -no_header -append -ieee ${o_file}/hgt_${s_yy}${m0}${d0}${h0}.grd

      else
        # Leaving the scalability.
        set i_full = ${i_dir}/grib/Z__C_RJTD_${s_yy}${m0}${d0}${h0}_EPSW_GPV_Rgl_FD00-08_grib2.bin 

        set i_list = ( sfc 850 500 300 )
        set il = 1

        while ( ${i} <= 4 )
 
          set C = `printf %4d ${I[${i}]}` 
          set level = `eval echo ${C} mb`
          wgrib2 -v ${i_file} | grep "UGRD" | grep "${level}" |wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/uwnd_${s_yy}${m0}${d0}${h0}.grd
          wgrib2 -v ${i_file} | grep "VGRD" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/vwnd_${s_yy}${m0}${d0}${h0}.grd
          wgrib2 -v ${i_file} | grep "TMP" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/tmp_${s_yy}${m0}${d0}${h0}.grd
          wgrib2 -v ${i_file} | grep "HGT" | grep "${level}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_file}/hgt_${s_yy}${m0}${d0}${h0}.grd
          @ i = ${i} + 1
        
        end

      endif

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
