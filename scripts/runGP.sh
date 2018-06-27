#!/bin/bash
tsdirll=$3/tsdirLL_$2
tsdirhl=$3/tsdirHL_$2

if [ -d batch$1 ]; then  rm -r batch$1 ; fi
mkdir batch$1
#copy only the neccesary files.
cp tsscds.dat batch$1/tsscds.dat
echo "tsdirll $tsdirll"  >> batch$1/tsscds.dat
echo "tsdirhl $tsdirhl"  >> batch$1/tsscds.dat
cp $2.xyz batch$1
if [ -f thdist  ]; then cp  thdist batch$1 ; fi
(cd batch$1 && tsscds.sh tsscds.dat >tsscds.log)

