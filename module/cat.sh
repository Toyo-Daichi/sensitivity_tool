#!/bin/sh
#set -ex
outdir=$1;date=$2;ft=$3;dataset=$4

if [ ${dataset} = "RISH" ]; then
  cat ${outdir}/UGRD.grd ${outdir}/VGRD.grd ${outdir}/HGT.grd ${outdir}/TMP.grd > ${outdir}/${date}_${ft}hr.grd

elif [ ${dataset} = "TIGGE" ]; then
  cat ${outdir}/UGRD.grd ${outdir}/VGRD.grd ${outdir}/HGT.grd ${outdir}/TMP.grd ${outdir}/SPFH.grd ${outdir}/PS.grd \
  > ${outdir}/${date}_${ft}hr.grd

elif [ ${dataset} = "LETKF" ]; then
  cat ${outdir}/UGRD.grd ${outdir}/VGRD.grd ${outdir}/TMP.grd ${outdir}/SPFH.grd ${outdir}/HGT.grd ${outdir}/WWND.grd \
      ${outdir}/VOR.grd  ${outdir}/SLP.grd  ${outdir}/PS.grd ${outdir}/RAIN.grd \
  > ${outdir}/${date}.grd
fi

exit
