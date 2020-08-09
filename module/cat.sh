#!/bin/sh
#set -ex
outdir=$1;date=$2;ft=$3

cat ${outdir}/UGRD.grd ${outdir}/VGRD.grd ${outdir}/HGT.grd ${outdir}/TMP.grd > ${outdir}/${date}_${ft}hr.grd
exit
