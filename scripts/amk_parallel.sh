#!/bin/bash
# default sbatch FT2
#SBATCH --output=amk_parallel-%j.log
#SBATCH --time=04:00:00
# partition selection

#_remove_this_in_ft_SBATCH -p shared --qos=shared
#SBATCH -c 1 --mem-per-cpu=2048
#SBATCH -n 8

# SBATCH --partition=cola-corta,thinnodes
# SBATCH -c 1
# SBATCH -n 24


#exe=$(basename $0)
# under batchs systems the scripts are copied to a generic script (in slurm slurm_script)
exe="amk_parallel.sh"
source utils.sh

#current working dir
cwd=$PWD
# Printing the references of the method
print_ref

#Checking the input of this script and also the presence of some files
re='^[0-9]+$'
if [ $# -eq 2 ]; then
   if [[ ! $2 =~ $re ]]; then usage "Second argument must be a number";fi
   noj1=$(( $(sort -nr <(find . -maxdepth 1 -type d -print | grep 'batch' | sed 's@batch@@;s@./@@') | head -n1) +1 ))
   nojf=$(( $noj1 + $2 -1))
elif [ $# -eq 3 ]; then
   if [[ ! $3 =~ $re ]]; then usage "Third argument must be a number";fi
   noj1=$2
   nojf=$3
else
   usage "At least 2 arguments required"
fi
echo "Making directories batch${noj1}-batch${nojf}"
echo ""

inputfile=$1

#input file exists?
if [ ! -f $inputfile ]
then
   usage "The first argument of this script must be the inputfile"
fi

molecule=`awk '{if($1=="molecule") print $2}'  $inputfile`
if [ -z $molecule ]; then echo "no molecule field, please check inputfile";exit 1;fi
if [ -f ${molecule}.xyz ]; then
   cp ${molecule}.xyz ${molecule}_ref.xyz
else
   echo ${molecule}.xyz is missing
   exit
fi
natom=$(awk 'NR==1{print $1}' ${molecule}.xyz)      

if [ "$inputfile" == "amk.dat" ]; then
   echo "The name of the input file is amk.dat"
   echo ""
else
   echo "The name of the input file is $inputfile"
   echo "amk.dat is a symbolic link of $inputfile"
   echo ""
   ln -sf $inputfile amk.dat
fi

###Adding path to tsdirll and tsdirhl
###Create tsdirll folder
tsdirll=$cwd/tsdirLL_$molecule
if [ ! -d "$tsdirll" ]; then
   echo "tsdirll does not exist. It will be created"
   mkdir $tsdirll 
else
   echo "tsdirll already exists."
fi
###
sqlite3 ${tsdirll}/track.db "create table if not exists track (id INTEGER PRIMARY KEY,nts  INTEGER, noj1 INTEGER, nojf INTEGER, ntraj INTEGER, emin REAL, emax REAL, permin INTEGER, permax INTEGER);"
##sampling must be 1 or 2
sampling=` awk '{if($1=="sampling") {if($2=="microcanonical") print "1";if($2=="canonical") print "2";if($2=="association") print "3";if($2=="external") print "4"}}'  $inputfile `
if [ $sampling -gt 2 ]; then
   echo With this sampling $exe cannot be employed
   exit 1
else
   if [ $sampling -eq 1 ]; then
      et="energies"
      uet="kcal/mol"
      flag="etraj"
   else
      et="temperatures"
      uet="K"
      flag="temp"
   fi
   erange="$(awk '{if($1=="'$flag'") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile | sed 's/-/ /')"

   nf="$(echo "$erange" | awk '{print NF}')"
   if [ $nf -eq 1 ]; then
      emin0="$(echo $erange | awk '{printf "%8.2f",$1}')"
      emax0="$(echo $erange | awk '{printf "%8.2f",$1}')"
   elif [ $nf -eq 2 ]; then
      emin0="$(echo $erange | awk '{printf "%8.2f",$1}')"
      emax0="$(echo $erange | awk '{printf "%8.2f",$2}')"
   elif [ $nf -eq 0 ]; then
      echo Range of $et is not provided in $inputfile and it will automatically be determined
      s=$(echo "3*$natom-6" | bc )
      emin0=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
      emax0=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
      if [ $sampling -eq 2 ]; then
         lstnm=`awk 'BEGIN{lstnm=0};{if($1=="atoms" && NF==3) {lstnm=$3}};END{print lstnm}' $inputfile `
         thmass=`awk 'BEGIN{thmass=0};{if($1=="thmass") {thmass=$2}};END{print thmass}' $inputfile `
         nlms=`awk 'BEGIN{atoms=0};{if($1=="atoms" && $2!="all") atoms=$2;if($1=="atoms" && $2=="all") atoms=0};END{print atoms}' $inputfile `

         awk '{print $0}
         END{
         print "100"
         print '$nlms'
         if('$nlms'>0) print '$lstnm'
         print '$thmass'
         }' ${molecule}.xyz | termo.exe > /dev/null

         natefin=$(awk '/Number of atoms to be excited/{print $NF}' fort.66)
         emin0=$(echo "335.51*$emin0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
         emax0=$(echo "335.51*$emax0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
      fi
   else
      echo Check the value of $flag in $inputfile
      exit 1
   fi
#set value of factor
   id=$(sqlite3 ${tsdirll}/track.db "select max(id) from track" | awk '{print $1+1-1}')
   if [ $id -gt 0 ]; then
      permin=$(sqlite3 ${tsdirll}/track.db "select permin from track where id='$id'")
      permax=$(sqlite3 ${tsdirll}/track.db "select permax from track where id='$id'")
      if [ $permin -gt 60 ]; then
         factormin=0.9
      elif [ $permin -ge 0 ] && [ $permin -le 60 ];then
         factormin=$(echo "10/9" | bc -l)
      else
         factormin=1
      fi
      if [ $permax -gt 60 ]; then
         factormax=$(echo "10/9" | bc -l)
      elif [ $permax -ge 0 ] && [ $permax -le 60 ];then
         factormax=0.9
      else
         factormax=1
      fi
   else
      factormin=1
      factormax=1
   fi
   emin=$(echo "$emin0*$factormin" | bc -l | awk '{printf "%8.2f",$1}')
   emax=$(echo "$emax0*$factormax" | bc -l | awk '{printf "%8.2f",$1}')
#faf is employed to avoid emin>emax situations
   faf=$(echo "$emax-$emin" | bc -l | awk 'BEGIN{faf=1};{if($1<0) faf=0};END{print faf}')
   if [ $faf -eq 0 ]; then
      s=$(echo "3*$natom-6" | bc )
      emin_sug=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
      emax_sug=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
      if [ $sampling -eq 2 ]; then
         emin_sug=$(echo "335.51*$emin_sug/$natom" | bc -l | awk '{printf "%8.2f",$1}')
         emax_sug=$(echo "335.51*$emax_sug/$natom" | bc -l | awk '{printf "%8.2f",$1}')
      fi
      echo You may consider changing the values of $flag in $inputfile
      echo Suggested range of ${et}: ${emin_sug}-${emax_sug}
      emin=$emin0
      emax=$emax0
   fi
   echo New range of ${et}: ${emin}-${emax} ${uet}
   sqlite3 ${tsdirll}/track.db "insert into track (noj1,nojf,emin,emax,permin,permax) values ($noj1,$nojf,$emin,$emax,-1,-1);"
   tmpinp="$(awk '{if($1=="sampling")
      {print $0
      print "'$flag' erange"}
   else if($1!="'$flag'") print $0}' $inputfile)"
   if [[ "$emin" == "$emax" ]]; then
      echo "$tmpinp" | sed 's/erange/'"$emin"'/' >$inputfile
   else
      echo "$tmpinp" | sed 's/erange/'"$emin"'-'"$emax"'/' >$inputfile
   fi
fi


#First we submit a tors.sh calc for every single relevant minimum
tors.sh >tors.log &
#Then we submit the nojf-noj1+1 trajectory jobs
#ft2 slurm
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  if (( $nojf-$noj1+1 < $SLURM_NTASKS )); then 
    echo "WARNING: Number of trajectory jobs ($nojf-$noj1+1) lower than allocated tasks ($SLURM_NTASKS)."
  fi
fi
doparallel "runGP.sh {1} $molecule $cwd" "$(seq $noj1 $nojf)"

