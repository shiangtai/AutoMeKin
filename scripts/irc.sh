#!/bin/bash
# default sbatch FT2
#SBATCH --output=irc-%j.log
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
exe="irc.sh"
source utils.sh
#remove tmp files
tmp_files=(tmp* bbfs.* *.arc *.mop coordir )
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

cwd=$PWD
sharedir=${TSSCDS}/share
#check the arguments of the script
if [ $# -gt 0 ]; then
   ci=$1
else
   ci="proceed"
fi

if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
if [ $ci != "screening" ] && [ $ci != "proceed" ]; then
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "                  Wrong argument                        "
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "To check what screening has done execute this script as:"
   echo "$exe screening"
   echo ""
   echo "To proceed with the irc execute this script as:"
   echo "$exe"
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   exit
fi
###Do screnning before anything else
SCREENING.sh $inputfile
if [ $ci == "screening" ]; then 
   echo ""
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "           Presumably you may want to check what SCREENING.sh has done            "
   echo "   Please check redundant and fragmented structures indicated in screening.log    "
   echo " If there are not what you expected you might change avgerr,bigerr and/or thdiss  "
   echo "Then, you can carry on with the IRC calculations, run this script without argument"
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo ""
   exit
fi
##Once your screening is safe, delete log 
rm -rf screening.log
#####
molecule=`awk '{if($1=="molecule") print $2}' $inputfile`
natom="$(awk 'NR>2{if(NF==4) ++natom};END{print natom}' ${molecule}.xyz)"
minref=$molecule"_ref"
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$molecule'"}' $inputfile`
method=` awk 'BEGIN{llcalc="PM7"};{if($1=="LowLevel") {$1="";llcalc=$0}};END{print llcalc}' $inputfile `
charge=` awk '/charge/{print $2}'  $inputfile `
TKMC=`awk '{if($1 == "TKMC") {print $2;nend=1}};END{if(nend==0) print "298"}' $inputfile`
tslist=${tsdirll}/tslist

if [ ! -d "$tsdirll/MINs" ]; then
   echo "MINs does not exist. It will be created"
   mkdir $tsdirll/MINs
else
   echo "MINs already exists."
fi
##create table for min
sqlite3 ${tsdirll}/MINs/min.db "create table if not exists min (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, sigma INTEGER,unique(name));"

# Optimize minref and calculate freqs

dum1="$(sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template1)"
dumfreq="$(awk 'BEGIN{one="1"};{if(NF==4) print $1,$2,one,$3,one,$4,one}' $minref.xyz)"
dumoldgeofreq="$(sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2)"
echo "$dum1"           > $molecule"_freq.mop"
echo "$dumfreq"       >> $molecule"_freq.mop"
echo "$dumoldgeofreq" >> $molecule"_freq.mop"
mopacl $molecule"_freq.mop" 2>/dev/null

# First we copy min0 in MIN directory 
echo "Moving min0 to its final location"
if [ -f $tsdirll/MINs/min0.out ]; then
   echo "Calcs completed for min0"
else
   geom="$(get_geom_mopac.sh $molecule"_freq.out" | awk '{if(NF==4) print $0}')"
   sed 's/thermo/thermo('$TKMC','$TKMC')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template >  $tsdirll/MINs/min0.mop
   echo "$geom"  >> $tsdirll/MINs/min0.mop
   mopacl $tsdirll/MINs/min0.mop 2>/dev/null
   e0=`awk '/HEAT OF FORMATION =/{e=$5};END{print e}' $tsdirll/MINs/min0.out `
   zpe0=`awk '/          ZERO POINT ENERGY/{zpe=$4};END{print zpe}' $tsdirll/MINs/min0.out `

   g_corr0=`awk 'BEGIN{t='$TKMC'}
      /          ZERO POINT ENERGY/{zpe=$4}
      /CALCULATED THERMODYNAMIC PROPERTIES/{ok=1}
   {if(ok==1 && $1 == '$TKMC') {
   getline;getline;getline;getline;
   h=$3/1000;s=$5/1000;print zpe+h-t*s;exit}
   }' $tsdirll/MINs/min0.out `
#   natom=$(echo "$geom" | awk 'END{print NR}')
   name=min0_0
#   echo $natom > $tsdirll/MINs/${name}.rxyz
#   echo "E=" $e0 "ZPE=" $zpe0 "Gcorr" $g_corr0 >> $tsdirll/MINs/${name}.rxyz
#   echo "$geom" >> $tsdirll/MINs/${name}.rxyz
   freq="$(get_freq_mopac.sh $tsdirll/MINs/min0.out)"
   sigma=$(awk '/SYMMETRY NUMBER/{print $NF;exit}' $tsdirll/MINs/min0.out)
#   echo "$freq" >> $tsdirll/MINs/${name}.rxyz
   sqlite3 ${tsdirll}/MINs/min.db "insert into min (natom,name,energy,zpe,g,geom,freq,sigma) values ($natom,'$name',$e0,$zpe0,$g_corr0,'$geom','$freq',$sigma);"
fi

# Now we do things specific of IRC 
if [ ! -d "$tsdirll/IRC" ]; then
   echo "IRC does not exist. It will be created"
   mkdir $tsdirll/IRC
else
   echo "IRC already exists."
fi
if [ ! -d "$tsdirll/TSs" ]; then
   echo "TSs does not exist. It will be created"
   mkdir $tsdirll/TSs
else
   echo "TSs already exists."
fi
#n=0
m=0
sqlite3 ${tsdirll}/inputs.db "drop table if exists mopac; create table mopac (id INTEGER PRIMARY KEY,name TEXT, unique(name));"
for name in $(awk '{print $3}' $tslist)
do
#  echo $name
#EMN
  if [ -f $tsdirll/TSs/$name"_thermo.out" ] && [ -f $tsdirll/IRC/$name"_ircf.out" ] && [ -f $tsdirll/IRC/$name"_ircr.out" ]; then
     calc1=$(awk 'BEGIN{calc=1};/MOPAC DONE/{calc=0};END{print calc}' $tsdirll/TSs/$name"_thermo.out")
     calc2=$(awk 'BEGIN{calc=1};/MOPAC DONE/{calc=0};END{print calc}' $tsdirll/IRC/$name"_ircf.out")
     calc3=$(awk 'BEGIN{calc=1};/MOPAC DONE/{calc=0};END{print calc}' $tsdirll/IRC/$name"_ircr.out")
     if [ $calc1 -eq 0 ] && [ $calc2 -eq 0 ] && [ $calc3 -eq 0 ]; then
        calc=0
     else
        calc=1
     fi
  else
     calc=1
  fi
#  if [ -f $tsdirll/TSs/$name"_thermo.out" ]; then
  if [ $calc -eq 0 ]; then
#EMN
    echo "Calcs completed for" $name
  else
     ((m=m+1))
     get_geom_mopac.sh $tsdirll/$name".out" | awk '{if(NF==4) print $0}' > tmp_geom
     sed 's/thermo/thermo('$TKMC','$TKMC')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template > $tsdirll/TSs/$name"_thermo.mop"
     cat tmp_geom >> $tsdirll/TSs/$name"_thermo.mop" 
     sed 's/method/'"$method"' charge='$charge' irc=1/g' $sharedir/freq_template1 > $tsdirll/IRC/$name"_ircf.mop"
     sed 's/method/'"$method"' charge='$charge' irc=-1/g' $sharedir/freq_template1 > $tsdirll/IRC/$name"_ircr.mop"
     cat tmp_geom >> $tsdirll/IRC/$name"_ircf.mop"
     cat tmp_geom >> $tsdirll/IRC/$name"_ircr.mop"
     sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2 | sed 's/force/cycles=5000 recalc=1/g'  >> $tsdirll/IRC/$name"_ircf.mop"
     sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2 | sed 's/force/cycles=5000 recalc=1/g'  >> $tsdirll/IRC/$name"_ircr.mop"
#     names[$m]=$name
    echo -e "insert or ignore into mopac values (NULL,'$name');\n.quit" | sqlite3 ${tsdirll}/inputs.db
  fi
done
echo Performing a total of $m irc calculations
#Perform m parallel calculations
if [ $m -gt 0 ]; then
#ft2 slurm
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  if (( $m < $SLURM_NTASKS )); then 
    echo "WARNING: Number of irc calculations ($m) lower than allocated tasks ($SLURM_NTASKS)."
  fi
fi
#   doparallel "runirc.sh $tsdirll" "$(echo ${names[@]})"
   doparallel "runirc.sh {1} $tsdirll" "$(seq $m)"
fi

