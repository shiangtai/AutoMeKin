#!/bin/bash
exe=$(basename $0)
#the number of args is 2 or 3: temp calc and arg-->optional
if [ $# -ne 2 ] && [ $# -ne 3 ]; then
   echo "Execute this script as:"
   echo "$exe temp calc (RXNet_opt)"
   echo "temp is the new temperautre is K"
   echo "calc is either ll (low-level) or hl (high-level)"
   echo "RXNet_opt can be:"
   echo "allstates: the RXNet is constructed including all minima (not grouped by conformational isomers)"
   echo "modify: this allows you to modify the conformational isomers after a first run of $exe"
   echo "Without RXNet_opt, $exe will construct the RXNet with each set of conformational isomers forming a single state"
   exit
##set up the value of ci (RXNet_opt)
elif [ $# -eq 2 ]; then
  ci=1
  echo "Every set of conformational isomers is a single state"
elif [ $# -eq 3 ]; then
  if [ $3 == "allstates" ];then
     ci=2
     echo "Every conformational isomer is a single state"
  elif [ $3 == "modify" ]; then
     ci=0
  else
     echo "Execute this script as:"
     echo "$exe temp calc (RXNet_opt)"
     echo "temp is the new temperautre is K"
     echo "calc is either ll (low-level) or hl (high-level)"
     echo "RXNet_opt can be:"
     echo "allstates: the RXNet is constructed including all minima (not grouped by conformational isomers)"
     echo "modify: this allows you to modify the conformational isomers after a first run of $exe"
     echo "Without RXNet_opt, $exe will construct the RXNet with each set of conformational isomers forming a single state"
     exit
  fi
fi
temperature=$1
calc=$2
calc_up=$(echo $calc | tr '[:lower:]' '[:upper:]' )
cwd=$PWD
source utils.sh
#On exit remove tmp files
tmp_files=(tmp* temp*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
natom=` awk 'NR==1{print $1}'  $molecule.xyz `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$molecule'"}' $inputfile` 
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirHL_'$molecule'"}' $inputfile` 
rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
##Checking that rate and temperature are in the inputfile
if [ $rate -eq 1 ]; then
   echo "$exe recalculates the kinetics for a canonical ensemble for a different temperature"
   echo "And you have not chosen canonical in the kinetics section"
   exit
fi
if [ $calc == "ll" ]; then
   tsdir=$tsdirll
   sort_script=sort.sh
   rxn_script=rxn_network1.sh
   ide_script=identifyprod.sh
   kmc_script=kmc.sh
   fin_script=final.sh
   tablemin=minnr
   tablets=ts
   factor=1
elif [ $calc == "hl" ]; then
   tsdir=$tsdirhl
   sort_script=SORT.sh
   rxn_script=RXN_NETWORK1.sh
   ide_script=IDENTIFYPROD.sh
   kmc_script=KMC.sh
   fin_script=FINAL.sh
   tablemin=minnrhl
   tablets=tshl
   factor=627.51
else
   echo "Wrong choice for the calc type. Please choose between: ll or hl"
   exit 0
fi
###Doing the calcs
#updating ${tablets}.db table
echo Updating $tablets table
for name in $(sqlite3 $tsdir/TSs/${tablets}.db "select name from $tablets")
do
   energy="$(sqlite3 $tsdir/TSs/${tablets}.db "select energy,zpe from $tablets where name='$name'" | sed 's@|@ @' )"
   sigma=$(sqlite3 $tsdir/TSs/${tablets}.db "select sigma from $tablets where name='$name'")
##create temp file tmp_rxyz from sqlite tables
   echo $natom  >tmp_rxyz
   echo E= "$(echo $energy | awk '{print $1}')" zpe= "$(echo $energy | awk '{print $2}')" g_corr= 0 sigma= $sigma >>tmp_rxyz
   sqlite3 $tsdir/TSs/${tablets}.db "select geom from $tablets where name='$name'" >> tmp_rxyz
   sqlite3 $tsdir/TSs/${tablets}.db "select freq from $tablets where name='$name'" | awk '{if($1>0) print $1}' >> tmp_rxyz 
   g_corr=$(thermochem.py tmp_rxyz $temperature $calc | awk '/Thermal correction to Gib/{print $9/'$factor'}' )
#   echo "Updating G value of $name-->$g_corr old g= $g"
   sqlite3 $tsdir/TSs/${tablets}.db "update $tablets set g='$g_corr' where name='$name';"
done
echo Now updating $tablemin table
###updating ${tablemin}.db table
for name in $(sqlite3 $tsdir/MINs/norep/${tablemin}.db "select name from $tablemin")
do
   energy="$(sqlite3 $tsdir/MINs/norep/${tablemin}.db "select energy,zpe from $tablemin where name='$name'" | sed 's@|@ @' )"
   sigma=$(sqlite3 $tsdir/MINs/norep/${tablemin}.db "select sigma from $tablemin where name='$name'")
##create temp file tmp_rxyz from sqlite tables
   echo $natom  >tmp_rxyz
   echo E= "$(echo $energy | awk '{print $1}')" zpe= "$(echo $energy | awk '{print $2}')" g_corr= 0 sigma= $sigma >>tmp_rxyz
   sqlite3 $tsdir/MINs/norep/${tablemin}.db "select geom from $tablemin where name='$name'" >>tmp_rxyz
   sqlite3 $tsdir/MINs/norep/${tablemin}.db "select freq from $tablemin where name='$name'" | awk '{if($1>0) print $1}' >>tmp_rxyz 
   g_corr=$(thermochem.py tmp_rxyz $temperature $calc | awk '/Thermal correction to Gib/{print $9/'$factor'}' )
#   echo "Updating G value of $name-->$g_corr old g= $g"
   sqlite3 $tsdir/MINs/norep/${tablemin}.db "update $tablemin set g='$g_corr' where name='$name';"
done
###
##
echo "Sorting the structures again according the their new G values"
$sort_script
echo "Creating a new RXN network with the new G values"
$rxn_script $ci
echo "Getting the formula of the products"
$ide_script
#change the temp in the original inputfile and save orig in temp file
#since inputfile is a symbolic link, search for the original file
origif=$(readlink -n $inputfile)
orig0=$(basename $origif .dat) 
cp $origif temp_$origif
newif="$(awk '{if($1=="TKMC") print $1,'$temperature';else print $0}' $origif)"
echo "$newif" > $origif
##
echo "Running the kinetics with the new G values and temperature"
$kmc_script
echo "Reconstructing FINAL_${calc_up}_${orig0} folder"
$fin_script
#recover the original inputfile
cp temp_$origif $origif
