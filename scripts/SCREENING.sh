#!/bin/bash

source utils.sh
#On exit remove tmp files
tmp_files=(tors.log tmp* sprint.out deg* atsdum2.out ConnMat* labels mingeom* ScalMat) 
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
###
inputfile=$1
if [ $# -eq 0 ]; then
   echo "Please, provide the input file as the argument of this script"
   exit
fi
sampling=` awk '{if($1=="sampling") {if($2=="microcanonical") print "1";if($2=="canonical") print "2";if($2=="association") print "3";if($2=="external") print "4"}}'  $inputfile `

if [ -z "$sampling" ]; then
   echo "No sampling selected in the input file"
   exit
fi
if [ $sampling -eq 3 ]; then
   SPRINT_screening_assoc.sh $inputfile
else
   SPRINT_screening.sh  $inputfile
fi

