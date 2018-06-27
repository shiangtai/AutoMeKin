#!/bin/bash
#default sbatch FT2 resources
#SBATCH --time=12:00:00
#SBATCH -n 4
#SBATCH --output=hlcalcs-%j.log
#_remove_this_in_ft_SBATCH --partition=cola-corta,thinnodes
#SBATCH --ntasks-per-node=2
#SBATCH -c 12
#
# SBATCH -p shared --qos=shared
# SBATCH --ntasks-per-node=2
# SBATCH -c 10


exe="hlcalcs.sh"
#if no arguments are provided, then a gui pops up 
if [ $# -eq 0 ]; then
   nl=$(whereis zenity | awk '{print NF}')
   if [ $nl -eq 1 ]; then 
      echo "Sorry, but it seems that zenity is not installed in this linux distribution"
      echo "window.sh cannot be employed"
      exit 1
   fi
   FILE="$(zenity --file-selection --filename="$PWD/*.dat" --file-filter="*.dat" --title="Select the input file" 2> /dev/null)"
   inputfile="$(basename $FILE)"
   echo "Selected input file: $inputfile"
   runningtasks="$(zenity --forms --title="hlcalcs.sh GUI" --text="Add input data" \
      --add-entry="Max number of running tasks" 2>/dev/null  )"
elif [ $# -eq 1 ]; then
   inputfile=$1
   if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
      runningtasks=$SLURM_NTASKS
   else
     echo "With one argument it must be run under the SLURM batch system:"
     echo "sbatch $exe inputfile"
     exit 1
   fi
elif [ $# -eq 2 ]; then
   inputfile=$1
   runningtasks=$2
else
   echo You must provide zero or two arguments:
   echo "nohup $exe >hlcalcs.log 2>&1 &"
   echo or
   echo "nohup $exe inputfile runningtasks >hlcalcs.log 2>&1 &"
   exit 1
fi
###Are we in the right folder?
if [ ! -f $inputfile ];then
   echo "$inputfile is not in this folder"
   exit 1
fi
export runningtasks
###
source utils.sh
#tmp_files=(tmp*)
#trap 'err_report $LINENO' ERR
#trap cleanup EXIT INT
# Printing the references of the method
print_ref
##
#set interactive mode to 0
inter=0
export inter
echo $$ > .script.pid
#
system="$(basename $inputfile .dat)"
echo "Start of Gaussian calculations for $system"
echo "*********************"
echo "**TS  calcs**"
echo "*********************"
TS.sh $inputfile > /dev/null
echo "*********************"
echo "**IRC calcs**"
echo "*********************"
IRC.sh  >/dev/null 
echo "*********************"
echo "**MIN calcs**"
echo "*********************"
MIN.sh  >/dev/null 
echo "*********************"
echo "**RXN calcs**"
echo "*********************"
RXN_NETWORK.sh >/dev/null
echo "*********************"
echo "**KMC calcs**"
echo "*********************"
KMC.sh >/dev/null
echo "*********************"
echo "**PROD calcs**"
echo "*********************"
PRODs.sh >/dev/null
echo "*********************"
echo "**FINAL calcs**"
echo "*********************"
FINAL.sh >/dev/null 
echo "Succesfully completed the calculations"
