#!/bin/bash
#default sbatch FT2 resources
#SBATCH --time=08:00:00
#SBATCH -n 4
#SBATCH --output=MIN-%j.log
#_remove_this_in_ft_SBATCH --partition=cola-corta,thinnodes
#SBATCH --ntasks-per-node=2
#SBATCH -c 12
#
# SBATCH -p shared --qos=shared
# SBATCH --ntasks-per-node=2
# SBATCH -c 10

sharedir=${AMK}/share

exe="MIN.sh"
source utils.sh
#exe=$(basename $0)
#current working dir
cwd=$PWD
#On exit remove tmp files
tmp_files=(atsdum2.out ConnMat deg* labels ScalMat sprint*)
#trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

#check if the inputfile has not been deleted
if [ -f amk.dat ];then
   echo "amk.dat is in the current dir"
   inputfile=amk.dat
else
   echo "amk input file is missing. You sure you are in the right folder?"
   exit
fi
#####

#Make sure g09 is submitted to slurm in ft2
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  t=$(srun -N 1 -n 1 g09<${sharedir}/g09_test.dat | awk 'BEGIN{t=0};/Normal ter/{t=1};END{print t}')
else
  t=$(g09<${sharedir}/g09_test.dat | awk 'BEGIN{t=0};/Normal ter/{t=1};END{print t}')
fi

if [ $t -eq 0 ]; then
   echo "Please check that gaussian09 is installed in your computer and it can be invoked as g09"
   exit 1
fi


molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "tsdirHL_'$molecule'"}' $inputfile`


##############################
##     Now opt the minima   ##
##############################
vdw=`awk 'BEGIN{vdw=0};{if($1=="vdw") vdw=1};END{print vdw}' $inputfile `
if [ $vdw -eq 1 ]; then
   cp thdist thdist_backup
   nvdw=`awk '{if($1=="vdw") print $2}' $inputfile `
   for i in $(seq 1 $nvdw)
   do
      at1[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $1; i++}  }' $inputfile `
      at2[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $2; i++}  }' $inputfile `
      dis[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i') print $3; i++}  }' $inputfile `
      awk '{if($1=="'${at1[$i]}'" && $2=="'${at2[$i]}'") {print $1,$2,"'${dis[$i]}'"}
      else if($2=="'${at1[$i]}'" && $1=="'${at2[$i]}'") {print $1,$2,"'${dis[$i]}'"}
      else {print $0}
      } ' thdist >thdist_vdw
      cp thdist_vdw thdist
   done
fi

echo "now the Minima dir "
n=0
sqlite3 ${tsdirhl}/IRC/inputs.db "drop table if exists gaussian; create table gaussian (id INTEGER PRIMARY KEY,name TEXT, input TEXT, unique (name));"
for i in $(sqlite3 $tsdirhl/TSs/tshl.db "select name from tshl")
do
  if [ -f $tsdirhl/IRC/minf_$i.log ] && [ -f $tsdirhl/IRC/minr_$i.log ]; then
     calc1=$(awk 'BEGIN{calc=1};/Normal termi/{calc=0};/Error termi/{calc=0};END{print calc}' $tsdirhl/IRC/minf_$i.log)
     calc2=$(awk 'BEGIN{calc=1};/Normal termi/{calc=0};/Error termi/{calc=0};END{print calc}' $tsdirhl/IRC/minr_$i.log)
     if [ $calc1 -eq 0 ] && [ $calc2 -eq 0 ]; then
        calc=0
     else
        calc=1
     fi
  else
     calc=1
  fi

#  if [ -f $tsdirhl/IRC/minf_$i.log ] && [ -f $tsdirhl/IRC/minr_$i.log ]; then
  if [ $calc -eq 0 ]; then
    echo "Min opt completed for" $i
  else
    ((n=n+2))
    echo "Submit Mins opt calc for" $i
#gettting the minima structures from the IRC output files
    get_minfminr_g09.sh $i 
#set-up gaussin09 calculation for $i
#    names[$n]=$i
####
  fi 
done
#Perform n parallel calculations
echo Performing a total of $n opt calculations
if [ $n -gt 0 ]; then
#   doparallel "runMIN.sh $PWD/$tsdirhl/IRC" "$(echo ${names[@]})"
   doparallel "runMIN.sh {1} $PWD/$tsdirhl/IRC" "$(seq $n)"
fi

if [ $vdw -eq 1 ]; then
   cp thdist_backup thdist
fi


