#!/bin/bash

#Function for usage of tsscds_parallel
function usage {
   echo $*
   echo "Execute this script as in this example:"
   if [ "$exe" == "slurm_script" ]; then exe="sbatch tsscds_parallel.sh";fi
   echo ""
   echo " $exe FA.dat 100"
   echo ""
   echo " where FA.dat is the inputfile and 100 is the total number of tasks"
   echo ""
   exit 1
}

#Function for usage of tsscds
function usages {
   echo $*
   echo "Execute this script as in this example:"
   echo "  $exe tsscds.dat "
   echo "where tsscds.dat is the inputfile "
   exit 1
}

#Function for usage of llcalcs
function usagell {
   echo $*
   echo "Execute this script as in this example:"
   echo "  $exe tsscds.dat nbatches niter"
   echo "where tsscds.dat is the inputfile"
   echo "nbatches is the number of batches"
   echo "and niter the number of interactions"
   exit 1
}

#Function to printout the references
function print_ref {
   build="$(awk '{print $1}' ${TSSCDS}/share/tsscds_build)"
   echo "*************************************************************************"
   echo "                   tsscds version 2018.${build}                           "
   echo "*************************************************************************"
   echo "*                                                                       *"
   echo "*  Cite this program as:                                                *"
   echo "*                                                                       *"
   echo "*  1) E. Martinez-Nunez, Phys. Chem. Chem. Phys., 2015, 17, 14912.      *"
   echo "*                                                                       *"
   echo "*  2) E. Martinez-Nunez, J. Comput. Chem., 2015, 36, 222.               *"
   echo "*                                                                       *"
   echo "*  3) MOPAC2016, James J. P. Stewart, Stewart Computational Chemistry,  *" 
   echo "*     Colorado Springs, CO, USA, HTTP://OpenMOPAC.net (2016).           *" 
   echo "*                                                                       *"
   echo "*************************************************************************"
}

#Function to submit jobs using slurm
function slurm {
#lets the user specify memory
#$SLURM_MEM_PER_CPU defined when option --mem-per-cpu= is used
#$SLURM_MEM_PER_NODE defined when option --mem= is used
   MEMPERCORE=$(sinfo -e -n $SLURM_NODELIST  -N -o "%m %c" -h | awk '{if(NR==1){min=$1/$2}else{new=$1/$2;if(new<min)min=new}}END{print min}')
   corespertask=${SLURM_CPUS_PER_TASK=1}
   if [ ! -z $SLURM_MEM_PER_NODE ]
   then
     #--ntasks-per-node= compulsory
     if [ ! -z $SLURM_NTASKS_PER_NODE ]
     then
       MEMPERCORE=$(( $SLURM_MEM_PER_NODE/ ($SLURM_NTASKS_PER_NODE * $corespertask) ))
     else
       echo "Please specify --ntasks-per-node= at the sbatch invocation"
       echo "or use --mem-per-cpu= option"
       exit 1
     fi
   fi
   if [ ! -z $SLURM_MEM_PER_CPU ]
   then
     MEMPERCORE=$SLURM_MEM_PER_CPU
   fi

#  SRUN="srun --exclusive -N1 -n1 --mem-per-cpu=$MEMPERCORE"
   SRUN="srun -N1 -n1 --mem=$(( $MEMPERCORE*$corespertask )) -c $corespertask --cpu_bind=none"
   runningtasks=$SLURM_NTASKS
}


#Function to submit jobs in parallel
function doparallel {
if [ ! -d tsscds_parallel-logs ];then mkdir tsscds_parallel-logs;fi

#slurm job?
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]
then
  slurm
else
# use as many concurrent tasks as number of cores-1 (if it is not defined or ill-defined)
  if [ -z $runningtasks ] || [ $runningtasks -gt $(nproc --ignore=1) ]; then
     runningtasks=$(nproc --ignore=1)
  fi
fi

# --delay 0.2 prevents overloading the controlling node
# -j is the number of tasks parallel runs
# --joblog makes parallel create a log of tasks that it has already run
# --resume-failed makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off

#progress bar only in interactive mode
if [ -z $inter ] && [ -z "$SRUN" ]; then
   parallel="parallel --bar --delay 0.2 -j $runningtasks --joblog tsscds_parallel-logs/${exe}-${iter}-tasks.log"
else
   parallel="parallel --delay 0.2 -j $runningtasks --joblog tsscds_parallel-logs/${exe}-${iter}-task.log"
fi
# this runs the parallel command we want
# in this case, we are running a script named runGP.sh
# parallel uses ::: to separate options. Here {0..99} is a shell expansion
# so parallel will run the command passing the numbers 0 through 99
# via argument {1}
COMMAND=$1
TASKS=$2
#$parallel "runGP.sh {1} $molecule" ::: $(seq $noj1 $nojf)
if [ -z "$SRUN" ]; then
   if [ -z $inter ]; then
      if [ -z $iter ]; then
         $parallel $COMMAND ::: $TASKS 2> >(zenity --progress --auto-close --no-cancel --title="parallel progress bar $exe" --width=500 --height=100 2> /dev/null) &
      else
         $parallel $COMMAND ::: $TASKS 2> >(zenity --progress --auto-close --no-cancel --title="parallel progress bar $exe iter=$iter" --width=500 --height=100 2> /dev/null)  &
      fi
   else
      nohup $parallel $COMMAND ::: $TASKS  >/dev/null 2>&1 & 
      echo $! > .parallel.pid
      wait
   fi
   echo $! > .parallel.pid
else
   $parallel $SRUN $COMMAND ::: $TASKS
fi
}

#function to printout the line in which the error occurred
function err_report {
    echo "Error on line $1 of $exe"
    rm -rf ${tmp_files[@]} 
    exit 1
}
#function to cleanup on exit
function cleanup {
    rm -rf ${tmp_files[@]} 
    echo ""
    echo "Cleaning up tmp files and exiting $exe"
}

function cleanup2 {
    rm -rf ${tmp_files[@]} 
}


function readhl {
   HLstring0="$(awk '{if($1=="HighLevel") print $2}' $inputfile)"
   HLstring="$(echo "$HLstring0" | sed 's@//@ @')"
   reduce=$(awk 'BEGIN{red=-1};{if($1=="HL_rxn_network") {if($2=="complete") red=0;if($2=="reduced" && NF==3) red=$3}};END{print red}' $inputfile)
   noHLcalc=$(echo $HLstring | awk 'BEGIN{nc=0};{nc=NF};END{print nc}')
   if [ $noHLcalc -eq 0 ]; then echo Please, provide HighLevel keyword ; exit ; fi
   iop=$(awk '{if($1=="iop") print $2}' $inputfile)
   level1=$(echo $HLstring | awk '{print $NF}')
   HLcalc1=$(echo "$level1" | sed 's@/@ @g;s@u@@g' | awk 'BEGIN{IGNORECASE=1};{if($1=="hf") m="HF";else if($1=="mp2") m="MP2"; else if($1=="ccsd(t)") m="CCSDT";else m="DFT"};END{print m}' )
#   echo High level calculations: "$HLstring0"
   if [ $noHLcalc -eq 1 ]; then
      HLcalc=$HLcalc1
   elif [ $noHLcalc -eq 2 ]; then
     level2=$(echo $HLstring | awk '{print $1}')
     HLcalc2=$(echo "$level2" | sed 's@/@ @g;s@u@@g' | awk 'BEGIN{IGNORECASE=1};{if($1=="hf") m="HF";else if($1=="mp2") m="MP2"; else if($1=="ccsd(t)") m="CCSDT";else m="DFT"};END{print m}' )
     HLcalc=$HLcalc2
   fi
}
