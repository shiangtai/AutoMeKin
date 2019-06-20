#!/bin/bash
if [ -f amk.dat ];then
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
cwd=$PWD
molecule=`awk '{if($1=="molecule") print $2}' $inputfile`
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$molecule'"}' $inputfile`

###view tslist file
if [ -f ${tsdirll}/tslist ]; then
   printf "  ts #  MOPAC file name  w_imag    Energy     w1     w2     w3     w4 traj #   Folder\n"
   printf "  ----  ---------------  ------    ------   ----   ----   ----   ---- ------   ------\n"
   awk '{printf "%6s%17s %6.0fi %9s %6.0f %6.0f %6.0f %6.0f   %4s %8s\n",$2,$3,$4,$5,$6,$7,$8,$9,$11,$13}' ${tsdirll}/tslist
else
   echo "tslist is not in ${tsdirll}"
fi
###

