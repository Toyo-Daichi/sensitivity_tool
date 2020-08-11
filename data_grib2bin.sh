#!/bin/csh -f

# set command
alias wgrib1 '/usr/bin/wgrib'
alias wgrib2 '/usr/bin/wgrib2'

set datapath = '/work3/daichi/Data/GSM_EnData/'

# set your target info.
set ft_list = ( anl 24 48 72)

foreach ft ( ${ft_list} )
echo ${ft}

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

      while ( ${s_hh} <= ${e_hh} )
      set h0 = `printf %02d ${s_hh}`

      set i_dir = ${datapath}/grib/${s_yy}${m0}${d0}
      set o_dir = ${datapath}/bin/${s_yy}${m0}${d0}
      mkdir -p ${o_dir}

      set accum_day = ${s_yy}${m0}${d0}${h0}
      if ( ${accum_day} <= 2006030100 ) then
        set nx = 144; set ny = 37; set nz = 4; set mem = 25 

        set i_surf = ${i_dir}/WFM12XSFC
        set i_850  = ${i_dir}/WFM12XPLL
        set i_500  = ${i_dir}/WFM12XPLM
        set i_300  = ${i_dir}/WFM12XPLH

        set i_list = ( ${i_surf} ${i_850} ${i_500} ${i_300} )
        set il = 1

        while ( ${il} <= 4 )
        set i_file = ${i_list[${il}]}

        if ( ${ft} == 'anl' ) then
          wgrib1 -v ${i_file} | grep "UGRD"  | grep "${ft}" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_${il}.grd
          wgrib1 -v ${i_file} | grep "VGRD"  | grep "${ft}" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_${il}.grd

          if ( ${il} == 1 ) then
            wgrib1 -v ${i_file} | grep "PRMSL"   | grep "${ft}" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
            wgrib1 -v ${i_file} | grep "APCP" | grep "0-72hr acc" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
            sleep 3s
          else if ( ${il} != 1 ) then
            wgrib1 -v ${i_file} | grep "HGT"  | grep "${ft}" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
            wgrib1 -v ${i_file} | grep "TMP"  | grep "${ft}" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
            sleep 3s
          endif

        else if ( ${ft} != 'anl' ) then
          wgrib1 -v ${i_file} | grep "UGRD"  | grep "${ft}hr fcst" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_${il}.grd
          wgrib1 -v ${i_file} | grep "VGRD"  | grep "${ft}hr fcst" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_${il}.grd

          if ( ${il} == 1 ) then
            wgrib1 -v ${i_file} | grep "PRMSL"   | grep "${ft}hr fcst" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
            wgrib1 -v ${i_file} | grep "APCP" | grep "0-${ft}hr acc" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
            sleep 3s
          else if ( ${il} != 1 ) then
            wgrib1 -v ${i_file} | grep "HGT"  | grep "${ft}hr fcst" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
            wgrib1 -v ${i_file} | grep "TMP"  | grep "${ft}hr fcst" | wgrib1 ${i_file} -i -bin -nh -o ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
            sleep 3s
          endif

        endif
        
        @ il = ${il} + 1
        end

      else if ( ${accum_day} <= 2020032300 ) then
        set nx = 144; set ny = 73; set nz = 4; set mem = 27 

        set i_file = ${i_dir}/Z__C_RJTD_${s_yy}${m0}${d0}${h0}0000_EPSW_GPV_Rgl_FD00-08_grib2.bin 
        set i_list = ( "ground" "850" "500" "300" )
        set il = 1

        while ( ${il} <= 4 )
          set level = ${i_list[${il}]}
          echo ${level}
          if ( ${ft} == 'anl' ) then
            wgrib2 -v ${i_file} | grep "UGRD" | grep "${level} mb" | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_${il}.grd
            wgrib2 -v ${i_file} | grep "VGRD" | grep "${level} mb" | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_${il}.grd
            if ( ${il} == 1 ) then
              wgrib2 -v ${i_file} | grep "PRMSL"   | grep "${ft}"    | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
              wgrib2 -v ${i_file} | grep "APCP"    | grep "0-8 day"  | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
              sleep 3s
            else if ( ${il} != 1 ) then
              wgrib2 -v ${i_file} | grep "HGT"  | grep "${level} " | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
              wgrib2 -v ${i_file} | grep "TMP"  | grep "${level} " | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
              sleep 3s
            endif

          else if ( ${ft} != 'anl' ) then
            @ ft_day = ${ft} / 24
            echo ${ft_day}

            wgrib2 -v ${i_file} | grep "UGRD" | grep "${level} " | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_${il}.grd
            wgrib2 -v ${i_file} | grep "VGRD" | grep "${level} " | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_${il}.grd
            if ( ${il} == 1 ) then
              wgrib2 -v ${i_file} | grep "PRMSL"   | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
              if ( ${ft_day} <= 2 ) then
                wgrib2 -v ${i_file} | grep "APCP"  | grep "0-${ft} hour"  | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
              else if ( ${ft_day} > 2 ) then
                wgrib2 -v ${i_file} | grep "APCP"  | grep "0-${ft_day} day"  | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
              endif
              sleep 3s
            else if ( ${il} != 1 ) then
              wgrib2 -v ${i_file} | grep "HGT"  | grep "${level} mb" | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${il}.grd
              wgrib2 -v ${i_file} | grep "TMP"  | grep "${level} mb" | grep "${ft}" | wgrib2 ${i_file} -i -no_header -append -ieee ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${il}.grd
              sleep 3s
            endif
          endif

          @ il = ${il} + 1
          #echo ${il}
        end

      else if ( ${accum_day} > 2020032300 ) then
        echo "CHECK your DATE or wget (on get_info.sh) Z__C_RJTD_yyyyMMddhhmmss_EPSG_GPV_Rgl_Gll1p25deg_FD0000-0100_grib2.bin"
        exit

      endif
      
      if ( ${ft} == 'anl' ) set ft = 00

      cat ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_?.grd >! ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_${ft}hr.grd
      cat ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_?.grd >! ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_${ft}hr.grd
      cat ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_?.grd  >! ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${ft}hr.grd
      cat ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_?.grd  >! ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${ft}hr.grd
      
      cat ${o_dir}/uwnd_${s_yy}${m0}${d0}${h0}_${ft}hr.grd \
          ${o_dir}/vwnd_${s_yy}${m0}${d0}${h0}_${ft}hr.grd \
          ${o_dir}/hgt_${s_yy}${m0}${d0}${h0}_${ft}hr.grd  \
          ${o_dir}/tmp_${s_yy}${m0}${d0}${h0}_${ft}hr.grd  \
      >!  ${o_dir}/${s_yy}${m0}${d0}${h0}_${ft}hr_${mem}mem.grd

      rm -rf  ${o_dir}/uwnd_*.grd ${o_dir}/vwnd_*.grd ${o_dir}/tmp_*.grd ${o_dir}/hgt_*.grd 
      
      @ s_hh = ${s_hh} + 6
      end
      
    set hh = 0
    @ s_dd = ${s_dd} + 1
    end

  @ s_mm = ${s_mm} + 1
  end

@ s_yy = ${s_yy} + 1
end

end


exit
