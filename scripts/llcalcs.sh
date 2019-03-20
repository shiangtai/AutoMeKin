#!/bin/bash
# default sbatch FT2
#SBATCH --output=llcalcs-%j.log
#SBATCH --time=04:00:00
# partition selection

#_remove_this_in_ft_SBATCH -p shared --qos=shared
#SBATCH -c 1 --mem-per-cpu=2048
#SBATCH -n 32

# SBATCH --partition=cola-corta,thinnodes
# SBATCH -c 1
# SBATCH -n 48

# first  arg is inputfile
# second arg is nbatches (200 is a good number)
# third  arg is niter  
#if no arguments are provided, then a gui pops up
exe="llcalcs.sh"
#exe="$(basename $0)"
if [ $# -eq 0 ]; then
#   nl=$(whereis zenity | awk '{print NF}')
#   if [ $nl -eq 1 ]; then 
#      echo "Sorry, but it seems that zenity is not installed in this linux distribution"
#      echo "window.sh cannot be employed"
#      exit 1
#   fi
   FILE="$(zenity --file-selection --filename="$PWD/*.dat" --file-filter="*.dat" --title="Select the input file" 2> /dev/null)"
   inputfile="$(basename $FILE)"
   echo "Selected input file: $inputfile"
   answer="$(zenity --forms --title="llcalcs.sh GUI" --text="Add input data" \
      --add-entry="Number of tasks" \
      --add-entry="Number of iterations" \
      --add-entry="Max number of running tasks" 2>/dev/null | awk 'BEGIN{FS="|"};{print $1,$2,$3}' )"
   nbatch=$(echo "$answer" | awk '{print $1}')
   niter=$(echo "$answer" | awk '{print $2}')
   runningtasks=$(echo "$answer" | awk '{print $3}')
elif [ $# -eq 3 ]; then
   inputfile=$1
   nbatch=$2
   niter=$3
   if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
      runningtasks=$SLURM_NTASKS
   else
     echo "With three arguments it must be run under the SLURM batch system:"
     echo "sbatch $exe inputfile ntasks niter"
     exit 1
   fi
elif [ $# -eq 4 ]; then
   inputfile=$1
   nbatch=$2
   niter=$3
   runningtasks=$4
else
   echo You must provide zero or four arguments:
   echo "nohup $exe >llcalcs.log 2>&1 &"
   echo or
   echo "nohup $exe inputfile ntasks niter runningtasks >llcalcs.log 2>&1 &" 
   exit 1
fi
export runningtasks

###Are we in the right folder?
if [ ! -f $inputfile ];then
   echo "$inputfile is not in this folder"
   exit 1
fi
if [ -z $nbatch ] || [ -z $niter ]; then
   echo "Number of batches and/or number of iterations have not been set"
   exit 1
fi
#EMN. If nbatch=0 do not run dynamics.
if [ $nbatch -eq 0 ]; then niter=1 ; fi
#EMN
###
iter=0
source utils.sh
#tmp_files=(tmp*)
#trap 'err_report $LINENO' ERR
#trap cleanup EXIT INT
# Printing the references of the method
print_ref
##
echo "A total of $nbatch batches of trajectories will be run for each of the $niter iterations"
system="$(basename $inputfile .dat)"
echo "Start of MOPAC search of transition states and intermediates of $system"
iter=1
#set interactive mode to 0
inter=0
export inter
echo $$ > .script.pid
#
while [ $iter -le $niter ]; do
   export iter
   echo "*********************"
   echo "    Iter $iter of $niter"
   echo "*********************"
   echo "$iter/$niter" > iter.txt
   if [ $nbatch -gt 0 ]; then
      tsscds_parallel.sh $inputfile $nbatch >/dev/null
   fi
   echo "**IRC calcs**"
   irc.sh > /dev/null
   echo "*********************"
   track_view.sh
   echo "*********************"
   echo "**MIN calcs**"
   min.sh  > /dev/null
   echo "**RXN_NETWORK calcs**"
   rxn_network.sh >/dev/null
   echo "Succesfully finished $iter iterations"
   echo "*********************"
   echo "**KMC calcs**"
   echo "*********************"
   kmc.sh > /dev/null
   echo "*********************"
   echo "*Making final folder*"
   echo "*********************"
   final.sh > /dev/null
   ((iter=iter+1))
done
echo "Succesfully completed the calculations"
