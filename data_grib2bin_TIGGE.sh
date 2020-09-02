#!/bin/sh
# Grib2 TIGGE data -> Grads format

# set parm.
center_list=( 'JMA' ) 
#center_list=( 'NCEP' 'CMC' ) #('UKMO' ->  not include Q)
#center_list=( 'JMA' 'ECMWF' 'NCEP' 'CMC' ) #('UKMO' ->  not include Q)
#ft_list=( 'anl' )
#ft_list=( 'anl' '24' '48' '72' )
level_list=( '1000' '925' '850' '700' '500' '300' '250' '200' )

#data info (from 30get_data.sh_FULL2.5deg)
#complist="u/v/t/q/gh";#levelist="1000/925/850/700/500/300/250/200"

for center in ${center_list[@]}; do

indata_path=/work1/mio/tigge_full/data/1.25deg/
outdata_path=/work3/daichi/Data/TIGGE/${center}/

#set date
yyyy=2018; e_yyyy=2018
mm=7     ; e_mm=7
dd=3     ; e_dd=7
hh=12    ; e_hh=12

M=`printf %2.2i ${mm}`
D=`printf %2.2i ${dd}`
H=`printf %2.2i ${hh}`

#ensemble size
case ${center} in
  "NCEP" ) mem=21 ;;
  "ECMWF") mem=51 ;;
  "CMC"  ) mem=21 ;;
  "UKMO" ) mem=18 ;;
esac

if [ "${center}" = 'JMA' ] && [ ${yyyy}${M}${D}${H}00 -lt 2006030100 ]; then
  mem=27
elif [ "${center}" = 'JMA' ] && [ ${yyyy}${M}${D}${H}00 -lt 2014022600 ]; then
  mem=51
elif [ "${center}" = 'JMA' ] && [ ${yyyy}${M}${D}${H}00 -ge 2014022600 ]; then
  mem=27
fi

while [ ${dd} -le ${e_dd} ]; do
  D=`printf %2.2i ${dd}`
  H=`printf %2.2i ${hh}`
  echo '@ READING DAY ' ${yyyy}${M}${D}${H} ${center} '@'
  
  o_dir=${outdata_path}/${yyyy}${M}${D}${H}/
  mkdir -p ${o_dir}
  i_U_file=${indata_path}/${yyyy}${M}/U_${center}_${yyyy}${M}${D}${H}.grib2
  i_V_file=${indata_path}/${yyyy}${M}/V_${center}_${yyyy}${M}${D}${H}.grib2
  i_Z_file=${indata_path}/${yyyy}${M}/Z_${center}_${yyyy}${M}${D}${H}.grib2
  i_T_file=${indata_path}/${yyyy}${M}/T_${center}_${yyyy}${M}${D}${H}.grib2
  i_Q_file=${indata_path}/${yyyy}${M}/Q_${center}_${yyyy}${M}${D}${H}.grib2
  i_PS_file=${indata_path}/${yyyy}${M}/PS_${center}_${yyyy}${M}${D}${H}.grib2

  for ft in ${ft_list[@]}; do
    i_level=1
    
    for level in ${level_list[@]}; do
      L=`printf %2.2i ${i_level}`
      if [ ${ft} = "anl" ]; then
        wgrib2 -v ${i_U_file}  | grep "${level} mb" | grep ":anl:" | wgrib2 ${i_U_file} -i -no_header -append -ieee ${o_dir}/uwnd_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_V_file}  | grep "${level} mb" | grep ":anl:" | wgrib2 ${i_V_file} -i -no_header -append -ieee ${o_dir}/vwnd_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_Z_file}  | grep "${level} mb" | grep ":anl:" | wgrib2 ${i_Z_file} -i -no_header -append -ieee ${o_dir}/hgt_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_T_file}  | grep "${level} mb" | grep ":anl:" | wgrib2 ${i_T_file} -i -no_header -append -ieee ${o_dir}/tmp_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_Q_file}  | grep "${level} mb" | grep ":anl:" | wgrib2 ${i_Q_file} -i -no_header -append -ieee ${o_dir}/spfh_${yyyy}${M}${D}${H}_${L}_level.grd
      elif [ ${ft} != "anl" ]; then
        wgrib2 -v ${i_U_file}  | grep "${level} mb" | grep ":${ft} hour fcst:" | wgrib2 ${i_U_file} -i -no_header -append -ieee ${o_dir}/uwnd_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_V_file}  | grep "${level} mb" | grep ":${ft} hour fcst:" | wgrib2 ${i_V_file} -i -no_header -append -ieee ${o_dir}/vwnd_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_Z_file}  | grep "${level} mb" | grep ":${ft} hour fcst:" | wgrib2 ${i_Z_file} -i -no_header -append -ieee ${o_dir}/hgt_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_T_file}  | grep "${level} mb" | grep ":${ft} hour fcst:" | wgrib2 ${i_T_file} -i -no_header -append -ieee ${o_dir}/tmp_${yyyy}${M}${D}${H}_${L}_level.grd
        wgrib2 -v ${i_Q_file}  | grep "${level} mb" | grep ":${ft} hour fcst:" | wgrib2 ${i_Q_file} -i -no_header -append -ieee ${o_dir}/spfh_${yyyy}${M}${D}${H}_${L}_level.grd
      fi
      
      #surface pressure
      if [ ${ft} = "anl" ]; then
        wgrib2 -v ${i_PS_file} | grep ":surface:" | grep ":anl:" | wgrib2 ${i_PS_file} -i -no_header -append -ieee ${o_dir}/ps_${yyyy}${M}${D}${H}_surf_level.grd
      elif [ ${ft} != "anl" ]; then
        wgrib2 -v ${i_PS_file} | grep ":surface:" | grep ":${ft} hour fcst:" | wgrib2 ${i_PS_file} -i -no_header -append -ieee ${o_dir}/ps_${yyyy}${M}${D}${H}_surf_level.grd
      fi
      i_level=`expr ${i_level} + 1`
    done

    #combine
    if [ ${ft} = "anl" ]; then ft='00' ;fi
    cat ${o_dir}/uwnd_${yyyy}${M}${D}${H}_*_level.grd > ${o_dir}/uwnd_${yyyy}${M}${D}${H}_${ft}hr.grd
    cat ${o_dir}/vwnd_${yyyy}${M}${D}${H}_*_level.grd > ${o_dir}/vwnd_${yyyy}${M}${D}${H}_${ft}hr.grd
    cat ${o_dir}/hgt_${yyyy}${M}${D}${H}_*_level.grd  > ${o_dir}/hgt_${yyyy}${M}${D}${H}_${ft}hr.grd
    cat ${o_dir}/tmp_${yyyy}${M}${D}${H}_*_level.grd  > ${o_dir}/tmp_${yyyy}${M}${D}${H}_${ft}hr.grd
    cat ${o_dir}/spfh_${yyyy}${M}${D}${H}_*_level.grd > ${o_dir}/spfh_${yyyy}${M}${D}${H}_${ft}hr.grd
    cat ${o_dir}/ps_${yyyy}${M}${D}${H}_*_level.grd   > ${o_dir}/ps_${yyyy}${M}${D}${H}_${ft}hr.grd
  
    cat ${o_dir}/uwnd_${yyyy}${M}${D}${H}_${ft}hr.grd \
        ${o_dir}/vwnd_${yyyy}${M}${D}${H}_${ft}hr.grd \
        ${o_dir}/hgt_${yyyy}${M}${D}${H}_${ft}hr.grd  \
        ${o_dir}/tmp_${yyyy}${M}${D}${H}_${ft}hr.grd  \
        ${o_dir}/spfh_${yyyy}${M}${D}${H}_${ft}hr.grd \
        ${o_dir}/ps_${yyyy}${M}${D}${H}_${ft}hr.grd   \
      > ${o_dir}/${yyyy}${M}${D}${H}_${ft}hr_${mem}mem.grd

   rm ${o_dir}/*_${yyyy}${M}${D}${H}_*_level.grd ${o_dir}/*_${yyyy}${M}${D}${H}_${ft}hr.grd
   
   sleep 0.5s
  done
  dd=`expr ${dd} + 1`

done
done

exit
