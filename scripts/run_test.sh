#!/bin/bash
# Run this script to test the program for the formic acid example
# You must have loaded amk module before:
# module load amk/2018

# Change the values below if needed:
# path_to_program: is the path to the installation directory ($HOME/amk-2018 by default)
# ntasks: is number of parallel tasks
# niter: is the number of iterations

path_to_program=$HOME/amk-2018
ntasks=4
niter=2
runningtasks=$ntasks

mkdir test
cd test
cp $path_to_program/examples/FA.* . 

nohup llcalcs.sh FA.dat $ntasks $niter $runningtasks >llcalcs.log 2>&1 &

