#!/bin/sh
#set -ex
outdir=$1;date=$2;ft=$3;dataset=$4

if [ ${dataset} = "RISH" ]; then
  cat ${outdir}/UGRD.grd ${outdir}/VGRD.grd ${outdir}/HGT.grd ${outdir}/TMP.grd > ${outdir}/${date}_${ft}hr.grd

elif [ ${dataset} = "TIGGE" ]; then
  cat ${outdir}/UGRD.grd ${outdir}/VGRD.grd ${outdir}/HGT.grd ${outdir}/TMP.grd ${outdir}/SPFH.grd ${outdir}/PS.grd \
  > ${outdir}/${date}_${ft}hr.grd

fi

exit
