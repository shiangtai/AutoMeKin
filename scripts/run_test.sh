#!/bin/bash
# Run this script to test the program for the formic acid example
# You must have loaded amk module before:
# module load amk/2018

# run_test.sh --prefix=path_to_program
# path_to_program: is the path to the installation directory ($HOME/amk-2018 by default)
# ntasks: is number of parallel tasks
# niter: is the number of iterations

if [ $# -ge 1 ]; then
   tpath=${@: -1}
   path_to_program=$(echo $tpath | sed  's/=/ /' | awk '{print $NF"/examples"}' | sed 's/\/\//\//')
else
   path_to_program=$HOME/amk-2018/examples
fi
file=$path_to_program/FA

if [ ! -f ${file}.dat ]; then
   echo "The path to amk-2018 is not correct"
   exit
fi
ntasks=4
niter=2
runningtasks=$ntasks

mkdir FA_test
cd FA_test
cp ${file}.dat ${file}.xyz .

nohup llcalcs.sh FA.dat $ntasks $niter $runningtasks >llcalcs.log 2>&1 &


