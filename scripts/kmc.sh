#!/bin/bash
source utils.sh
#remove tmp files
tmp_files=(tmp*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT
exe=$(basename $0)

if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
molecule=`awk '{if($1=="molecule") print $2}' $inputfile`

#####
##### units employed in the RRKM calc. default is ps-1
#####
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`

units=1e-12
temperature=` awk 'BEGIN{t=0};/TKMC/{t=$2};END{print t}'  $inputfile `
energy=` awk 'BEGIN{e=0};/EKMC/{e=$2};END{print e}'  $inputfile `
nmol=` awk 'BEGIN{nmol=1000};{if($1=="nmol") nmol=$2};END{print nmol}'  $inputfile `
imin=` awk 'BEGIN{imin="min0"};{if($1=="imin") imin=$2};END{print imin}'  $inputfile `
rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
step=` awk 'BEGIN{step=10};{if ($1=="Stepsize") step=$2};END{printf "%8.0f",step}'  $inputfile `
rxnfile=$(awk 'BEGIN{suf="_long.cg_groupedprods"};/PathInfo/{if($2=="Relevant") suf=".relevant"};END{print "RXNet"suf}' $inputfile )
catalysis=`awk 'BEGIN{cat=0};{if($1=="species") cat=1};END{print cat}' $inputfile`


if [ $rate -eq 0 ] && [ $temperature -eq 0 ] ; then
      echo "For a canonical ensemble please provide a temperature (TKMC)"
      exit
fi
if [ $rate -eq 1 ] && [ $energy -eq 0 ] ; then
      echo "For a microcanonical ensemble please provide an energy (EKMC)"
      exit
fi



if [ $imin == "ask" ]; then
   echo -n "Provide the label of the starting minimum: "
   read imin
fi
if [ $imin == "min0" ]; then
   minn=`awk '/min0/{print $2}' $tsdirll/MINs/SORTED/MINlist_sorted`
   imin=`awk 'BEGIN{min='$minn'}
   {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
   j=1
   while(j<=iso[NR]){
      if('$minn'==m[NR,j]) min=m[NR,1]
      j++
      }
   }
   END{print min}' $tsdirll/working/conf_isomer.out `
fi




if [ -f $tsdirll/KMC/$rxnfile ]; then
   echo "RXNfile is:" $tsdirll"/KMC/"$rxnfile
else
   echo "Specify a valid keyword for PathInfo"
   exit
fi
if [ $rate -eq 1 ]; then
   if [ -z "$energy" ]; then
      echo "Energy not given. Please provide an energy (in kcal/mol) using keyword EKMC" 
      exit
   fi 
   if [ ! -d "$tsdirll/KMC/RRKM" ]; then
      echo "KMC/RRKM does not exist. It will be created"
      mkdir $tsdirll/KMC/RRKM
   else
      echo "KMC/RRKM already exists. It will be created again"
      rm -r $tsdirll/KMC/RRKM
      mkdir $tsdirll/KMC/RRKM
   fi

elif [ $rate -eq 0 ]; then
   if [ -z "$temperature" ]; then
      echo "Temperature not given. Please provide a temperature (in K) using keyword TKMC" 
      exit
   fi 
   if [ ! -d "$tsdirll/KMC/TST" ]; then
      echo "KMC/TST does not exist. It will be created"
      mkdir $tsdirll/KMC/TST
   else
      echo "KMC/TST already exists. It will be created again"
      rm -r $tsdirll/KMC/TST
      mkdir $tsdirll/KMC/TST
   fi
fi
#Now, gather all information: freqs, barriers, etc to make TST/KMC analysis
#temperature is the temperature in K
#units is a factor from s-1 to (ps in this case)
#imin is the min from which we start the kinectics
#nmol is the number of molecules in the KMC calculations
#step is the step size for printout in the KMC calc
if [ -f $tsdirll/KMC/processinformation ]; then rm $tsdirll/KMC/processinformation ; fi
zero=0
echo "Starting minimum" $imin
nnpro=0
awk '{if(NR>2) print $0}' $tsdirll/KMC/$rxnfile > tmp_rxn
set `awk '{print NR}' tmp_rxn`
for i
do
  ((nnpro=nnpro+1))
  ts=`awk 'NR=='$i',NR=='$i'{print $3}' tmp_rxn`
  tsn=`awk 'NR=='$i',NR=='$i'{print $2}' tmp_rxn`
  procn=`awk 'NR=='$i',NR=='$i'{print $2}' tmp_rxn`
  min1=`awk 'NR=='$i',NR=='$i'{print $8}'  tmp_rxn`
  nproc=`awk 'BEGIN{nproc=1};NR=='$i',NR=='$i'{if($10~"MIN") nproc=nproc+1};END{print nproc}' tmp_rxn`
  state1=`awk 'NR=='$i',NR=='$i'{print $8}' tmp_rxn`
  state2=`awk 'NR=='$i',NR=='$i'{print $11}' tmp_rxn`
  lmin1="MIN"$min1
  lts="TS"$tsn
  echo "Proc" $nnpro "TS" $tsn $state1 $state2 >> $tsdirll/KMC/processinformation

  if [ $nproc -eq 2 ]; then
     deg1=`awk 'NR=='$i',NR=='$i'{print $15/$14/$12}' tmp_rxn`
  else
     deg1=`awk 'NR=='$i',NR=='$i'{print $14/$13/$12}' tmp_rxn`
  fi

#
  if [ $rate -eq 0 ]; then
     g1="$(awk 'NR=='$i',NR=='$i'{printf "%10.2f\n",$5}' tmp_rxn)"
     g2="$(awk '{if($2=='$min1') printf "%10.2f\n",$4}' $tsdirll/MINs/SORTED/MINlist_sorted)"
     deltag="$(echo "$g1" "$g2" | awk '{printf "%10.2f",$1-$2}')"
     echo $deltag $temperature $deg1 > $tsdirll/KMC/TST/proc1_TS${procn}.dat 
     tst.exe <$tsdirll/KMC/TST/proc1_TS${procn}.dat>$tsdirll/KMC/TST/proc1_TS${procn}.out
     rate1=`awk '{print $0}' $tsdirll/KMC/TST/proc1_TS${procn}.out` 
     echo $rate1 $state1 $state2 >> $tsdirll/KMC/TST/rate$temperature.out
     echo "Running proc1_TS${procn}"
  elif [ $rate -eq 1 ]; then
     ets=`awk 'NR=='$i',NR=='$i'{printf "%10.0f\n",349.75*$5}' tmp_rxn`
     errkm0=` awk '/EKMC/{printf "%10.0f\n",349.75*$2}'  $inputfile `
     e0=`awk '{if(NR==1) printf "%10.0f\n",349.75*$4}' $tsdirll/MINs/SORTED/MINlist_sorted`
#     ((errkm=errkm0-e0+1000))
     errkm=$(echo "$errkm0 - $e0 + 1000" | bc -l)
     egap=`awk '{if($2=='$min1') printf "%10.0f\n",349.75*$4}' $tsdirll/MINs/SORTED/MINlist_sorted`
     egapkcal=`awk '{if($2=='$min1') printf "%15.4f\n",$4}' $tsdirll/MINs/SORTED/MINlist_sorted`
     echo "Direct via TS$procn for process $state1 $state2 $egapkcal $energy" > $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
#     ((ebarrier=ets-egap))
     ebarrier=$(echo "$ets - $egap" | bc -l)
     if [ $ebarrier -lt 0 ]; then ebarrier=0 ; fi
     echo "$errkm,"$ebarrier",100,0,"$deg1",0,0" >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
     echo "0,0" >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
     echo "rrkm" >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
     echo "1.0" >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
###sqlite3
     sqlite3 $tsdirll/MINs/SORTED/mins.db "select freq from mins where name='$lmin1'" | wc -l >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
     sqlite3 $tsdirll/MINs/SORTED/mins.db "select freq from mins where name='$lmin1'" | awk '{printf "%10.0f\n",sqrt($1*$1)}'  >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
     sqlite3 $tsdirll/TSs/SORTED/tss.db "select freq from tss where name='$lts'" | awk 'END{print NR-1}'  >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
     sqlite3 $tsdirll/TSs/SORTED/tss.db "select freq from tss where name='$lts'" | awk '{if($1>0) printf "%10.0f\n",$1}'  >> $tsdirll/KMC/RRKM/proc1_TS${procn}.dat
###
     echo "Running proc1_TS${procn}"
     rrkm.exe <$tsdirll/KMC/RRKM/proc1_TS${procn}.dat>$tsdirll/KMC/RRKM/proc1_TS${procn}.out  
     awk 'BEGIN{ns=100000}
     /via/{p1=$6;p2=$7;e=('$energy'-$8)*349.75;if(e<0) {print "0.0",p1,p2;exit}}
     /k\(E/{ns=NR}
     {if(NR<(ns+2))
        dif=1000
     else {
        dif=e-$1
        if(dif<0 && NR==(ns+2)) {print "0.0",p1,p2;exit}
        if(dif<0 && NR>(ns+2)) {print $2*'$units',p1,p2;exit}} }'  $tsdirll/KMC/RRKM/proc1_TS${procn}.out >> $tsdirll/KMC/RRKM/rate$energy.out
  fi 
  if [ $nproc -eq 2 ]; then
     ((nnpro=nnpro+1))
     echo "Proc" $nnpro "TS" $tsn $state2 $state1 >> $tsdirll/KMC/processinformation
     min2=`awk 'NR=='$i',NR=='$i'{print $11}' tmp_rxn`
     lmin2="MIN"$min2
     deg2=`awk 'NR=='$i',NR=='$i'{print $15/$16/$13}' tmp_rxn`
     if [ $rate -eq 0 ]; then
        g2="$(awk '{if($2=='$min2') printf "%10.3f\n",$4}' $tsdirll/MINs/SORTED/MINlist_sorted)"
        deltag="$(echo "$g1" "$g2" | awk '{printf "%10.2f",$1-$2}')"
        echo $deltag $temperature $deg2 > $tsdirll/KMC/TST/proc2_TS${procn}.dat 
        tst.exe <$tsdirll/KMC/TST/proc2_TS${procn}.dat>$tsdirll/KMC/TST/proc2_TS${procn}.out
        rate2=`awk '{print $0}' $tsdirll/KMC/TST/proc2_TS${procn}.out` 
        echo $rate2 $state2 $state1 >> $tsdirll/KMC/TST/rate$temperature.out
        echo "Running proc2_TS${procn}"
     elif [ $rate -eq 1 ]; then
        egap=`awk '{if($2=='$min2') printf "%10.0f\n",349.75*$4}' $tsdirll/MINs/SORTED/MINlist_sorted`
        egapkcal=`awk '{if($2=='$min2') printf "%15.4f\n",$4}' $tsdirll/MINs/SORTED/MINlist_sorted`
        echo "Reverse via TS$procn for process $state2 $state1 $egapkcal $energy" > $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        deg2=`awk 'NR=='$i',NR=='$i'{print $15/$16/$13}' tmp_rxn`
#        ((ebarrier=ets-egap))
        ebarrier=$(echo "$ets - $egap" | bc -l)
        if [ $ebarrier -lt 0 ]; then ebarrier=0 ; fi
        echo "$errkm,"$ebarrier",100,0,"$deg2",0,0" >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        echo "0,0" >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        echo "rrkm" >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        echo "1.0" >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
###sqlite3
#nfreq_react
        sqlite3 $tsdirll/MINs/SORTED/mins.db "select freq from mins where name='$lmin2'" | wc -l >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        sqlite3 $tsdirll/MINs/SORTED/mins.db "select freq from mins where name='$lmin2'" | awk '{printf "%10.0f\n",sqrt($1*$1)}'  >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
#nfreq_ts
        sqlite3 $tsdirll/TSs/SORTED/tss.db "select freq from tss where name='$lts'" | awk 'END{print NR-1}'  >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        sqlite3 $tsdirll/TSs/SORTED/tss.db "select freq from tss where name='$lts'" | awk '{if($1>0) printf "%10.0f\n",$1}'  >> $tsdirll/KMC/RRKM/proc2_TS${procn}.dat
        echo "Running proc2_TS${procn}"
        rrkm.exe <$tsdirll/KMC/RRKM/proc2_TS${procn}.dat>$tsdirll/KMC/RRKM/proc2_TS${procn}.out  
        awk 'BEGIN{ns=100000}
        /via/{p1=$6;p2=$7;e=('$energy'-$8)*349.75;if(e<0) {print "0.0",p1,p2;exit}}
        /k\(E/{ns=NR}
        {if(NR<(ns+2))
           dif=1000
        else {
           dif=e-$1
           if(dif<0 && NR==(ns+2)) {print "0.0",p1,p2;exit}
           if(dif<0 && NR>(ns+2)) {print $2*'$units',p1,p2;exit}} }'  $tsdirll/KMC/RRKM/proc2_TS${procn}.out >> $tsdirll/KMC/RRKM/rate$energy.out
     fi
  fi
done

if [ $catalysis -eq 1 ]; then
   echo "Rates in KMC/TST"
   echo "Since this is a catalysis job, do not run KMC simulations."
   echo "They will be run at the end for all subsystems."
   exit
fi

if [ $rate -eq 0 ]; then
   echo "kmc calc" > $tsdirll/KMC/kmcT$temperature.dat
   nproc=`awk 'END{print NR}' $tsdirll/KMC/TST/rate$temperature.out`
   nspec=`awk 'BEGIN{max=0};{if($2 >max) max=$2;if($3>max) max=$3};END{print max}' $tsdirll/KMC/TST/rate$temperature.out`
   echo $nproc, $nspec, "1" >> $tsdirll/KMC/kmcT$temperature.dat
   cat $tsdirll/KMC/TST/rate$temperature.out >> $tsdirll/KMC/kmcT$temperature.dat
   for i in $(seq 1 $nspec)
   do
     if [ $i -eq $imin ]; then
        echo $nmol >> $tsdirll/KMC/kmcT$temperature.dat
     else
        echo "0" >> $tsdirll/KMC/kmcT$temperature.dat
     fi
   done
   echo "0" >> $tsdirll/KMC/kmcT$temperature.dat
   echo $step >> $tsdirll/KMC/kmcT$temperature.dat
   echo "Running KMC calc"
   kmc.exe <$tsdirll/KMC/kmcT$temperature.dat>$tsdirll/KMC/kmcT$temperature.out
   kmcfile=$tsdirll/KMC/kmcT$temperature.out
   postb="T$temperature"
elif [ $rate -eq 1 ]; then
   echo "kmc calc" > $tsdirll/KMC/kmcE$energy.dat
   nproc=`awk 'END{print NR}' $tsdirll/KMC/RRKM/rate$energy.out`
   nspec=`awk 'BEGIN{max=0};{if($2 >max) max=$2;if($3>max) max=$3};END{print max}' $tsdirll/KMC/RRKM/rate$energy.out`
   echo $nproc, $nspec, "1" >> $tsdirll/KMC/kmcE$energy.dat
   cat $tsdirll/KMC/RRKM/rate$energy.out >> $tsdirll/KMC/kmcE$energy.dat
   for i in $(seq 1 $nspec)
   do
     if [ $i -eq $imin ]; then
        echo $nmol >> $tsdirll/KMC/kmcE$energy.dat
     else
        echo "0" >> $tsdirll/KMC/kmcE$energy.dat
     fi
   done
   echo "0" >> $tsdirll/KMC/kmcE$energy.dat
   echo $step >> $tsdirll/KMC/kmcE$energy.dat
   kmc.exe <$tsdirll/KMC/kmcE$energy.dat>$tsdirll/KMC/kmcE$energy.out
   kmcfile=$tsdirll/KMC/kmcE$energy.out
   postb="E$energy"
fi

#getting the branching ratios
lastmin=`awk '{lm=$2};END{print lm}' $tsdirll/MINs/SORTED/MINlist_sorted `
awk 'BEGIN{ok=0}
/Population of every species/{point=NR;ok=1}
/counts per process/{ok=0}
{if(ok==1 && NR>point) print $0
}' $kmcfile > $tsdirll/KMC/branching$postb
echo "  %  Products" > $tsdirll/KMC/branching$postb.out
set `awk '{print $1}' $tsdirll/KMC/branching$postb `
for i
do
   popp=`awk 'NR=='$i',NR=='$i'{print $2/'$nmol'*100}' $tsdirll/KMC/branching$postb`
   if [ $i -gt $lastmin ] ; then
      code=`awk '{if($3=='$i') {print $2;exit}}' $tsdirll/PRODs/PRlist_kmc.log`
      name=$(awk '{if($2=='$code') print $3}' $tsdirll/PRODs/PRlist_kmc)
      namen=$(basename $name .rxyz)
      namesql="PR"$code"_"$namen
      prod="$(sqlite3 $tsdirll/PRODs/prod.db "select formula from prod where name='$namesql'")"
#      prod=`FormulaPROD.sh $tsdirll/PRODs/PR$code"_"*`
      printf "%8s %10s\n" $popp "$prod" >> $tsdirll/KMC/branching$postb.out
   fi 
done

