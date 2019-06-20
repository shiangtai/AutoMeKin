#!/bin/bash

source utils.sh
#On exit remove tmp files
tmp_files=(ConnMat tmp* labels mingeom ScalMat *.arc *.mop fort.* fortdir pesdir tswrk *_dyn* *_backup rotate.dat minn black_list*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT


cwd=$PWD
sharedir=${AMK}/share
# Printing the references of the method
print_ref
##
exe=$(basename $0)
echo "Using MOPAC development version for the dynamics"
echo ""
##dt=0.5 fs. ncyles=nfs*dnc
dnc=2
##
if [ $# -eq 0 ]; then usages "One argument is required" ; fi
inputfile=$1
if [ ! -f $inputfile ]; then
   echo "The file $inputfile does not exist"
   exit
fi
if [ "$inputfile" == "amk.dat" ]; then
   echo "Reading amk.dat"
   echo ""
else
   echo "Reading $inputfile"
   echo "amk.dat is a symbolic link of $inputfile"
   echo ""
   ln -sf $inputfile amk.dat
fi

nb=`basename $cwd`
###EMN. seed for srand. 
srandseed=$(echo $nb | awk '/batch/{print $0}' | sed 's@batch@@' | awk '{print $1+100}')
###EMN
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
charge=`awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' $inputfile `
multiple_minima=`awk 'BEGIN{mm=1};{if($1=="multiple_minima" && $2=="yes") mm=1};{if($1=="multiple_minima" && $2=="no") mm=0};END{print mm}' $inputfile `
sampling=` awk '{if($1=="sampling") {if($2=="microcanonical") print "1";if($2=="canonical") print "2";if($2=="association") print "3";if($2=="external") print "4"}}'  $inputfile `

echo "+++General Section+++"
##check that rate is ok
if [ $sampling -ne 3 ]; then
   rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
   tcheck=` awk 'BEGIN{t=0};/TKMC/{t=$2};END{print t}'  $inputfile `
   echeck=` awk 'BEGIN{e=0};/EKMC/{e=$2};END{print e}'  $inputfile `
   if [ -z $rate ]; then echo "Please provide a value for Rate (canonical or microcanonical)"; exit; fi
   if [ $rate -eq 0 ] && [ $tcheck -eq 0 ] ; then
         echo "For a canonical ensemble please provide a temperature (TKMC)"
         exit
   fi
   if [ $rate -eq 1 ] && [ $echeck -eq 0 ] ; then
         echo "For a microcanonical ensemble please provide an energy (EKMC)"
         exit
   fi
   emaxts=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$echeck'};{if($1=="MaxEn") en=$2};END{print 1.5*en}' $inputfile `
#fi
##check that xyz file is present
#if [ $sampling -ne 3 ]; then
   if [ ! -f $molecule.xyz ]; then
      echo "$molecule.xyz file does not exist"
      exit
   else
      natom=` awk 'NR==1{print $1}'  $molecule.xyz `
      awk '{if(NR>2 && NF==4) print $1}' $molecule.xyz > labels
      echo "Number of Atoms of the MS    = " $natom
      echo "Name of the Molecular System = " $molecule
#      min0=$(awk 'BEGIN{min0=0};{if($1=="min0") print $2};END{print min0}' $inputfile)
##sel_mol.sh will choose a minimum from either tsdirll or tsdirhl to perform iamk
#      sel_mol.sh $inputfile $min0
      sel_mol.sh $inputfile $multiple_minima
   fi
##
fi
echo "Charge of the molecule       = " $charge
method=` awk 'BEGIN{llcalc="PM7"};{if($1=="LowLevel") {$1="";llcalc=$0}};END{print llcalc}' $inputfile `
echo "Semiempirical Method         = " $method
echo ""
echo "+++CDS Section+++"
deltat=1
itrajn=` awk '{if($1=="ntraj") print $2}'  $inputfile `
ncycles=` awk 'BEGIN{time=500};{if($1=="fs") time=$2};END{print time*'$dnc'}'  $inputfile `
nfs=` awk 'BEGIN{time=500};{if($1=="fs") time=$2};END{print time}'  $inputfile `
irange=20
irangeo2=`echo "scale=0; $irange/2" | bc` 
dynamics_template="$(cat $sharedir/dynamics_template)"
###This is only in the development version AXD(nbondsfrozen) and a bias potential (nbondsbreak and nbondsform)
factorflipv=` awk '{if($1=="factorflipv") factor=$2};END{print factor}'  $inputfile `
if [ ! -z $factorflipv ]; then
   echo "   A factor of $factorflipv is employed to flip velocities"
   tmp0="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug vflip='$factorflipv'/')"
   dynamics_template="$tmp0"
fi
nbondsfrozen=` awk 'BEGIN{nbf=0};{if($1=="nbondsfrozen") nbf=$2};END{print nbf}'  $inputfile `
if [ $nbondsfrozen -gt 0 ]; then
   echo "   $nbondsfrozen bonds will not break thanks to AXD"
   tmp0="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug nbondsfrozen='$nbondsfrozen'/')"
   dynamics_template="$tmp0"
   for i in $(seq 1 $nbondsfrozen)
   do
      ibf[$i]=`awk '{if($1=="nbondsfrozen") {i=1;while(i<='$i'){getline;if(i=='$i') print $1;++i}}}' $inputfile`
      jbf[$i]=`awk '{if($1=="nbondsfrozen") {i=1;while(i<='$i'){getline;if(i=='$i') print $2;++i}}}' $inputfile`
      echo "   bond $i = ${ibf[$i]}-${jbf[$i]}"
      if [ $i -eq 1 ]; then
         echo ${ibf[$i]} ${jbf[$i]} > fort.67
      else
         echo ${ibf[$i]} ${jbf[$i]} >> fort.67
      fi
   done
fi
nbondsbreak=` awk 'BEGIN{nbb=0};{if($1=="nbondsbreak") nbb=$2};END{print nbb}'  $inputfile `
if [ $nbondsbreak -gt 0 ]; then
   echo "   $nbondsbreak bonds will break using biased dynamics simulations"
   tmp1="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug nbondsbreak='$nbondsbreak'/')"
   dynamics_template="$tmp1"
   kapparep=`awk 'BEGIN{kapparep=200};{if($1=="kapparep") kapparep=$2};END{print kapparep}'  $inputfile`
   kappa=`awk 'BEGIN{kappa=200};{if($1=="kappa") kappa=$2};END{print kappa}'  $inputfile`
   iexp=`awk 'BEGIN{iexp=1};{if($1=="iexp") iexp=$2};END{print iexp}'  $inputfile`
   rmin=`awk 'BEGIN{rmin=0.2};{if($1=="rmin") iexp=$2};END{print rmin}'  $inputfile`
   echo $kapparep $kappa $iexp $rmin > fort.70
   for i in $(seq 1 $nbondsbreak)
   do
      ibb[$i]=`awk '{if($1=="nbondsbreak") {i=1;while(i<='$i'){getline;if(i=='$i') print $1;++i}}}' $inputfile`
      jbb[$i]=`awk '{if($1=="nbondsbreak") {i=1;while(i<='$i'){getline;if(i=='$i') print $2;++i}}}' $inputfile`
      echo "   bond $i = ${ibb[$i]}-${jbb[$i]}"
      if [ $i -eq 1 ]; then
         echo ${ibb[$i]} ${jbb[$i]} > fort.68
      else
         echo ${ibb[$i]} ${jbb[$i]} >> fort.68
      fi
   done
fi
nbondsform=` awk 'BEGIN{nbfo=0};{if($1=="nbondsform") nbfo=$2};END{print nbfo}'  $inputfile `
if [ $nbondsform -gt 0 ]; then
   echo "   $nbondsform bonds will form using biased dynamics simulations"
   tmp2="$(echo "$dynamics_template" | sed 's/itry=100/itry=100 debug nbondsform='$nbondsform'/')"
   dynamics_template="$tmp2"
   kapparep=`awk 'BEGIN{kapparep=200};{if($1=="kapparep") kapparep=$2};END{print kapparep}'  $inputfile`
   kappa=`awk 'BEGIN{kappa=200};{if($1=="kappa") kappa=$2};END{print kappa}'  $inputfile`
   iexp=`awk 'BEGIN{iexp=1};{if($1=="iexp") iexp=$2};END{print iexp}'  $inputfile`
   rmin=`awk 'BEGIN{rmin=0.2};{if($1=="rmin") iexp=$2};END{print rmin}'  $inputfile`
   echo $kapparep $kappa $iexp $rmin > fort.70
   for i in $(seq 1 $nbondsform)
   do
      ibfo[$i]=`awk '{if($1=="nbondsform") {i=1;while(i<='$i'){getline;if(i=='$i') print $1;++i}}}' $inputfile`
      jbfo[$i]=`awk '{if($1=="nbondsform") {i=1;while(i<='$i'){getline;if(i=='$i') print $2;++i}}}' $inputfile`
      echo "   bond $i = ${ibfo[$i]}-${jbfo[$i]}"
      if [ $i -eq 1 ]; then
         echo ${ibfo[$i]} ${jbfo[$i]} > fort.69
      else
         echo ${ibfo[$i]} ${jbfo[$i]} >> fort.69
      fi
   done
fi
if [ $nbondsbreak -gt 0 ] || [ $nbondsform -gt 0 ]; then
   irange=`awk 'BEGIN{irange=20};{if($1=="irange") irange=$2};END{print irange}'  $inputfile`
   irangeo2=`echo "scale=0; $irange/2" | bc` 
   echo "Selected irange= $irange"
   echo "irangeo2= $irangeo2"
fi
###
if [ $sampling -eq 1 ]; then
###Energy can be a single value or a range of values
   erange="$(awk '{if($1=="etraj") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile | sed 's/-/ /')"
   nf="$(echo "$erange" | awk '{print NF}')"
   if [ $nf -eq 1 ]; then
      energy="$(echo $erange | awk '{print $1}')"
   elif [ $nf -eq 2 ]; then
      data="$(echo "$erange" | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')" 
      energy=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=`echo "scale=0; $irange/2" | bc` 
   elif [ $nf -eq 0 ]; then
      echo Range of energies is not provided in $inputfile and it will automatically be determined
      s=$(echo "3*$natom-6" | bc )
      emin0=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
      emax0=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
      data="$(echo $emin0 $emax0 | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')" 
      energy=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=`echo "scale=0; $irange/2" | bc` 
   else
      echo "Check the value of etraj"
      exit
   fi
###
#   energy=` awk '{if($1=="etraj") print $2}'  $inputfile `
   if [ -z "$energy" ]; then
      echo "   Microcanonical ensemble."
      echo "   Please provide an energy for the trajectories using keyword etraj"
      exit
   fi 
   lstnm=`awk 'BEGIN{lstnm=0};{if($1=="modes" && NF==3) {lstnm=$3}};END{print lstnm}' $inputfile `
   nlms=`awk 'BEGIN{modes=0};{if($1=="modes" && $2!="all") modes=$2;if($1=="modes" && $2=="all") modes=0};END{print modes}' $inputfile `
   seed=`awk 'BEGIN{seed=0};/seed/{seed=$2};END{print seed}' $inputfile `
   echo "Microcanonical sampling" 
   if [ $nlms -gt 0 ]; then
      echo "   Number of modes excited =" $nlms 
      echo "   Modes excited:" $lstnm 
   else
      echo "   All normal modes are excited "
   fi
   echo "   Energy (kcal/mol) =" $energy
   echo "   Selected value of irange =" $irange
elif [ $sampling -eq 2 ]; then
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
   trange="$(awk '{if($1=="temp") {for(i=2;i<=NF;i++) printf "%s ",$i}}'  $inputfile | sed 's/-/ /')"
   nf="$(echo "$trange" | awk '{print NF}')"
   if [ $nf -eq 1 ]; then
      energy="$(echo $trange | awk '{print $1}')"
   elif [ $nf -eq 2 ]; then
      data="$(echo "$trange" | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')"
      energy=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=`echo "scale=0; $irange/2" | bc` 
   elif [ $nf -eq 0 ]; then
      echo Range of temperatures is not provided in $inputfile and it will automatically be determined
      s=$(echo "3*$natom-6" | bc )
      emin0=$(echo "16.25*($s-1)" | bc -l | awk '{e=$1;if(e>400) e=400;printf "%8.2f",e}')
      emax0=$(echo "46.25*($s-1)" | bc -l | awk '{e=$1;if(e>1200) e=1200;printf "%8.2f",e}')
      tmin0=$(echo "335.51*$emin0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
      tmax0=$(echo "335.51*$emax0/$natefin" | bc -l | awk '{printf "%8.2f",$1}')
      data="$(echo $tmin0 $tmax0 | awk 'BEGIN{steps=3;srand('$srandseed');n=steps+1;rn=int(n*rand())}
      {le=$1;he=$2;range=he-le}
      END{
      delta=range/steps
      printf "%8.2f %8.0f",le+rn*delta,rn
      }')"
      energy=$(echo "$data" | awk '{printf "%8.2f",$1}' )
      irange=$(echo "$data" | awk '{rn=$2};END{print 20-rn*2}' )
      irangeo2=`echo "scale=0; $irange/2" | bc` 
   else
      echo "Check the value of temp"
      exit
   fi
   if [ -z "$energy" ]; then
      echo "   Canonical ensemble."
      echo "   Please provide a temperature for the trajectories using keyword temp"
      exit
   fi 
   echo "Canonical sampling" 
   if [ $nlms -gt 0 ]; then
      echo "   Number of atoms excited =" $nlms 
      echo "   Atoms excited:" $lstnm 
      echo "   The other atoms are frozen"
      echo "   Only atoms with masses greater than" $thmass "will receive kinetic energy"
   else
      echo "   All atoms are excited "
   fi
   echo "   Temperature (K) =" $energy
   echo "   Selected value of irange =" $irange
elif [ $sampling -eq 3 ]; then
   echo "Looking for A-B complexes"
elif [ $sampling -eq 4 ]; then
   echo "Coupling external dynamics results with BBFS algorithm"
else
   echo "No sampling provided. Please check your inputfile"
   echo "Number of  trajectories  = " $itrajn
   exit
fi
if [ $sampling -lt 3 ]; then
   echo "Number of  trajectories  = " $itrajn
   echo "Time(fs) = " $nfs
fi


freqmin=` awk '/freqmin/{print $2}'  $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$molecule'"}' $inputfile` 
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirHL_'$molecule'"}' $inputfile` 

rcwd=$(echo $tsdirll | sed 's@/tsdir@ @' | awk '{print $1}')
##Create tsll.db database with the low-level transition states
#sqlite3 ${rcwd}/tsll.db "create table if not exists n (id INTEGER PRIMARY KEY,name TEXT,wimag REAL,energy REAL,w1 REAL,w2 REAL,w3 REAL,w4 REAL,traj INTEGER,path TEXT);"

if [ $sampling -eq 3 ]; then
   specA=`awk '/A=/{print $2}' $inputfile`
   specB=`awk '/B=/{print $2}' $inputfile`
   if [ ! -f $specA.xyz ]; then
      echo "$specA.xyz file does not exist"
      exit
   fi
   if [ ! -f $specB.xyz ]; then
      echo "$specB.xyz file does not exist"
      exit
   fi
   assocdir=$cwd"/assoc_"$specA"_"$specB
   atom1rot=`awk '/rotate/{ff=$2;if(ff=="com") ff=-1};END{print ff}' $inputfile`
   atom2rot=`awk '/rotate/{ff=$3;if(ff=="com") ff=-1};END{print ff}' $inputfile`
   dist=`awk '/rotate/{print $4}' $inputfile` 
   distm=`awk '/rotate/{print $5}' $inputfile` 
   if [ ! -d "$assocdir" ]; then
      echo "$assocdir does not exist. It will be created"
      mkdir $assocdir
   else
      echo "$assocdir already exists."
   fi
fi


tslistll=${tsdirll}/tslist
bbfs=`awk 'BEGIN{p=0};/BBFS/{p=1};END{print p}' $inputfile`
#The default is fastmode
wrkmode=`awk 'BEGIN{mode=1};/slowmode/{mode=0};/fastmode/{mode=1};END{print mode}' $inputfile `
if [ $bbfs -eq 1 ]; then
   echo ""
   echo "+++BBFS Section+++"
   echo "TSs with imag freq (cm-1)   >= " $freqmin
   echo "      and Energy (kcal/mol) <= " $emaxts
   echo "TSs at  low level of theory in " $tsdirll
   echo "TSs at high level of theory in " $tsdirhl
   if [ $wrkmode -eq 0 ]; then echo "Slow mode. Standard Hessian update and 3 points per candidate in single paths" ; fi
   if [ $wrkmode -eq 1 ]; then echo "Fast mode. Hessian update every 10 steps and 1 point  per candidate in single paths " ; fi
   echo ""
   if [ ! -d "$tsdirll" ]; then
      echo "tsdirll does not exist. It will be created"
      mkdir $tsdirll 2>tmp_err
      if [ -s tmp_err ]; then
         echo "check the path of tsdirll folder"
         exit
      fi
   else
      echo "tsdirll already exists."
   fi

fi

####create folders only for sampling ne 3

if [ $sampling -ne 3 ]; then
   rm -rf fortdir
   mkdir  fortdir
   rm -rf tswrk
   mkdir  tswrk  
   if [ $sampling -ne 4 ]; then
      rm -rf coordir
      mkdir  coordir
   fi
   rm -rf pesdir
   mkdir  pesdir
fi
#####
frtmplt="$(cat $sharedir/freq_template1 | sed 's/method/'"$method"' charge='$charge'/g')"
if [ $wrkmode -eq 0 ]; then ts_template="$(cat $sharedir/ts_templateslow)" ; fi
if [ $wrkmode -eq 1 ]; then ts_template="$(cat $sharedir/ts_templatefast)" ; fi
dum_ts="$(echo "$ts_template" | sed 's/method/'"$method"' charge='$charge'/g')"
dumoldgeofreq="$(cat $sharedir/freq_template2 | sed 's/method/'"$method"' charge='$charge'/g')"
if [ $sampling -ne 3 ]; then
   echo ""
   echo "Performing frequency calculation"
   dumfreq="$(awk 'BEGIN{one="1"};{if(NR>2 && NF==4) print $1,$2,one,$3,one,$4,one}' $molecule.xyz)"
   echo "$frtmplt"        > $molecule"_freq.mop"
   echo "$dumfreq"       >> $molecule"_freq.mop"
   if [ $sampling -eq 1 ]; then 
      echo "$dumoldgeofreq" >> $molecule"_freq.mop"
   fi
   mopacl $molecule"_freq.mop" 2>/dev/null
fi

if [ $bbfs -eq 1 ]; then
###########min0 is the reference minimum. If it exists, take its energy
   if [ -f $tsdirll/MINs/min.db ]; then
      e0=$(sqlite3 ${tsdirll}/MINs/min.db "select energy from min where name='min0_0'")
   else
      e0=$(awk '/FINAL HEAT OF FORMATION =/{e0=$6};END{print e0}' ${molecule}_freq.out )
   fi
###########3
   edum=`echo "$emaxts+$e0" | bc | awk '{printf "%10.0f",$1}'` 
   emaxts=$edum
   thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
#this is stuff for bbfs only (connectivity matrix and things like that)
   get_geom_mopac.sh $molecule"_freq.out" | awk '{if(NF==4) print $0}' >mingeom
   createthdist.sh $thd
   createMat.sh
fi


if [ $sampling -eq 1 ]; then 
   echo $seed > tmp_qv
   initialqv_mopac_samp1.sh $molecule"_freq.out" $energy $nlms $lstnm >> tmp_qv
elif [ $sampling -eq 2 ]; then
   initialqv_mopac_samp2.sh $molecule"_freq.out" $energy $nlms $lstnm $thmass > tmp_qv
elif [ $sampling -eq 3 ]; then
###First, we need to select the best possible specA 
   select_AandB.sh $inputfile
###
   n1=`awk 'NR==1,NR==1{print $1;exit} ' $specA.xyz`
   n2=`awk 'NR==1,NR==1{print $1;exit} ' $specB.xyz`
   n="$(echo $n1 $n2 | awk '{print $1+$2}')"
   echo $n $n1 $dist $distm > rotate.dat
   echo $atom1rot $atom2rot >> rotate.dat
   awk '{if(NF==4) print $0}' $specA.xyz >>rotate.dat
   awk '{if(NF==4) print $0}' $specB.xyz >>rotate.dat
   for i in $(seq 1 100)
   do
      sed 's/method/'"$method"' charge='$charge' bonds/g' $sharedir/freq_template1 > $assocdir/assoc$i".mop" 
      rotate.exe <rotate.dat>>$assocdir/assoc$i".mop" 
      sed 's/method/'"$method"' charge='$charge'/g' $sharedir/freq_template2 >>$assocdir/assoc$i".mop" 
   done
   inter=0
   doparallel "runAS.sh {1} $assocdir" "$(seq 1 100)" 
   echo "MOPAC Jobs terminated"
   echo "Now the SCREENING"
   SCREENING.sh $inputfile
   exit
fi


echo "Performing accelerated dynamics calculations"
for i in $(seq 1 $itrajn) 
do 
####
#This is only for internal dynamics (MOPAC)
####
if [ $sampling -lt 3 ]; then
  if [ "$(ls -A fortdir)" ]; then rm fortdir/* ; fi
  if [ "$(ls -A pesdir)" ]; then rm pesdir/* ; fi
  if [ "$(ls -A tswrk)" ]; then rm tswrk/* ; fi
  if [ $sampling -eq 1 ]; then 
     nm.exe <tmp_qv>tmp_geomd
     dytem0="$(echo "$dynamics_template" | sed 's/method/'"$method"' charge='$charge'/g')"
  elif [ $sampling -eq 2 ]; then
     termo.exe <tmp_qv | sed 's/ 1.d69/1.d69/g' >tmp_geomd
     dytem0="$(echo "$dynamics_template" | sed 's/method/'"$method"' charge='$charge' large/g')"
  fi
  dytem1="$(echo "$dytem0" | sed 's/ncycles/'$ncycles'/;s/deltat/'$deltat'/')"
  echo "$dytem1"  > $molecule"_dyn"$i".mop"
  cat   tmp_geomd   >> $molecule"_dyn"$i".mop"
  echo ""
  echo ""
  echo "+-+-+-+-+-+-+-+-+- Trajectory ${i} +-+-+-+-+-+-+-+-+-"
  if [ $sampling -eq 1 ]; then 
    mopacl $molecule"_dyn"$i".mop" &> $molecule"_dyn"$i".log"
    if [ ! -f $molecule"_dyn"$i".xyz" ]; then
       echo $molecule"_dyn"$i".xyz does not exist"
       break
    else
       mv $molecule"_dyn"$i".xyz" coordir
    fi
  elif [ $sampling -eq 2 ]; then
    mopacl $molecule"_dyn"$i".mop" &> $molecule"_dyn"$i".log"
    if [ ! -f $molecule"_dyn"$i".xyz" ]; then
       echo $molecule"_dyn"$i".xyz does not exist"
       break
    else
       mv $molecule"_dyn"$i".xyz" coordir
    fi
  fi
  for ik in *_dyn*
  do
   rm $ik
  done
  if [ $bbfs -eq 0 ]; then
     echo "End of traj "$i
     echo "Only trajs. No BBFS  "
     break
  fi
else
  echo "************************************************"
  echo "        External Dynamics. Traj "$i
  echo "************************************************"
fi
###########
#From here everything is common for internal and external dynamics
###########
#  echo "Running BBFS algorithm"
#  echo $i > $filewnotr
  echo $irange > bbfs.dat
  snapshots_mopac.sh coordir/$molecule"_dyn"$i".xyz" $i

  if [ $i -eq 1 ]; then
     bbfs.exe<bbfs.dat>bbfs.out
  else
     bbfs.exe<bbfs.dat>>bbfs.out
  fi
  path=`awk '/Number of paths/{np=$4};END{print np}' bbfs.out `
  if [ $path -eq 0 ]; then
     echo "This traj has no paths "
     continue
  fi

  for j in $(seq 1 $path)
  do
    mv fort.$j?? fortdir
  done
  nc=0
  for j in $(ls fortdir)
  do
    ((nc=nc+1))
    n="$(basename $j | sed 's/fort.//g')"
    paste labels fortdir/$j >tmp_geomp
    echo "$frtmplt"   > pesdir/pes$n
    cat tmp_geomp     >> pesdir/pes$n
#    cat pesdir/pes$n >> pesdir/pes0
#    echo "" >> pesdir/pes0
  done

# poner un flag para hacer todo el pes0
#  echo "Running the whole path"
#  mopacl pesdir/pes0
#  pes.sh pesdir/pes0.out > pesdir/pes0.log
#
  echo "Npaths=" $path
#EMN
#EMN
  chapath[0]=0
  for ip in $(seq $path)
  do
# Find the highest energy point
     ((tspt = $ip*100+$irangeo2))
        ijc=`awk '/Joint path=/{if($2=='$ip') ijc=$5};END{print ijc}' bbfs.out`
##If previous path was multiple, break
        chapath[$ip]=$ijc
        jp=`echo "scale=0; $ip-1" | bc` 
        if [ ${chapath[$jp]} -eq 1 ]; then continue ; fi
##

        if [ $ijc -eq 0 ] && [ $wrkmode -eq 0 ]; then echo "Path" $ip" (Single): 3 attempts to locate the ts" ; fi
        if [ $ijc -eq 0 ] && [ $wrkmode -eq 1 ]; then echo "Path" $ip" (Single): 1 attempt  to locate the ts" ; fi
        if [ $ijc -eq 1 ]; then echo "Path" $ip" (Multiple): several attempts to locate the ts" ; fi
# middle point
        if [ $ijc -eq 1 ]; then
           ll=-1
           ((ul=$irangeo2-1))  
           if [ $irangeo2 -ge 8 ]; then 
              dlt=2
           else
              dlt=1
           fi
        else
           ul=1
           if [ $wrkmode -eq 0 ]; then
              ll=-1
              dlt=1
           else
              ll=0
              dlt=2
           fi
        fi
        npo=0
        for itspt in $(seq $ll $dlt $ul)
        do 
           ((npo=npo+1))
           ((ctspt = tspt + itspt))
           mopacl pesdir/pes$ctspt 2>/dev/null
           geom.sh pesdir/pes$ctspt.out $natom
           if [ -f "tswrk/geom1" ]; then 
              paste labels tswrk/geom1 > tmp_geomts
           else
              echo "Partial opt failed"
              continue
           fi
           echo "$dum_ts"         > tswrk/ts$i"_"$ip"_"$ctspt
           cat tmp_geomts           >> tswrk/ts$i"_"$ip"_"$ctspt
           echo "$dumoldgeofreq" >> tswrk/ts$i"_"$ip"_"$ctspt
           mopacl tswrk/ts$i"_"$ip"_"$ctspt 2>/dev/null
#
           int=`awk 'BEGIN{intv=0};/Too many variables/{intv=1};END{print intv}' tswrk/ts$i"_"$ip"_"$ctspt".out"`
           if [ $int -eq 1 ]; then
              sed -i 's/ts /ts int /g' tswrk/ts$i"_"$ip"_"$ctspt
#              echo "Repeating TS opt with int coordinates!"
              mopacl tswrk/ts$i"_"$ip"_"$ctspt 2>/dev/null
           fi
#           rm tswrk/ts$i"_"$ip"_"$ctspt
           fe="$(mopac_freq_ts.sh tswrk/ts$i"_"$ip"_"$ctspt".out")"
           fi="$(echo "$fe" | awk '{printf "%10.0f",$1}')"
           ei="$(echo "$fe" | awk '{printf "%10.0f",$2}')"
           if [[ ("$fi" -eq -1) ]]; then
              printf "     Pt%2s: failed-->Lowest real freq is negative\n" $npo
              continue
           elif [[ ("$fi" -eq -2) ]]; then
              printf "     Pt%2s: failed-->Sum of 2 lowest real freqs < 10cm-1\n" $npo
              continue
           elif [[ ("$fi" -eq -3) ]]; then
              printf "     Pt%2s: failed-->Stationary point is a minimum\n" $npo
              continue
           elif [[ ("$fi" -eq -4) ]]; then
              printf "     Pt%2s: failed-->EF algorithm was unable to optimize a TS\n" $npo
              continue
           elif [[ ("$fi" -eq -5) ]]; then
              printf "     Pt%2s: failed-->This is a vdw ts (separated products)\n" $npo
              continue
           elif [[ ("$ei" -gt "$emaxts") ]]; then
              printf "     Pt%2s: TS optimized but not added-->E=%4s kcal/mol, which is greater than %4s kcal/mol\n" $npo $ei $emaxts
              continue
           fi
           if [[ ("$fi" -ge "$freqmin") ]]; then
              f="$(echo "$fe" | awk '{printf "%10.2f",$1}')" 
              e="$(echo "$fe" | awk '{printf "%10.2f",$2}')"  
              f1="$(echo "$fe" | awk '{printf "%10.2f",$3}')" 
              f2="$(echo "$fe" | awk '{printf "%10.2f",$4}')" 
              f3="$(echo "$fe" | awk '{printf "%10.2f",$5}')" 
              f4="$(echo "$fe" | awk '{printf "%10.2f",$6}')"
# GLB added lock to tslist so that duplicate numbers are not created
              (
              flock -x 200 || exit 1
              if [ -f "$tslistll" ]; then
                 ok=$(diff.sh $f $e $f1 $f2 $f3 $f4 $tslistll)
                 if [[ ("$ok" -eq "-1") ]]; then
                    nt=`awk '{nt=$2};END{print nt}' $tslistll `
                    ((nt = nt + 1))
                    name=ts${nt}_${nb}
                    printf "ts%5s%18s%9s%9s%9s%9s%9s%9s traj=%4s Path= %10s\n" $nt $name $f $e $f1 $f2 $f3 $f4 $i $nb  >> $tslistll
                    cp tswrk/ts${i}_${ip}_${ctspt}.out  $tsdirll/${name}.out
                    printf "     Pt%2s: TS optimized and added to ts list\n" $npo
                 else
                    printf "     Pt%2s: TS optimized but not added-->redundant with ts %4s\n" $npo $ok
                 fi
              else
                 nt=1
                 name=ts${nt}_${nb}
                 printf "ts%5s%18s%9s%9s%9s%9s%9s%9s traj=%4s Path= %10s\n" $nt $name $f $e $f1 $f2 $f3 $f4 $i $nb  >> $tslistll
                 cp tswrk/ts${i}_${ip}_${ctspt}.out  $tsdirll/${name}.out
                 printf "     Pt%2s: TS optimized and added to ts list\n" $npo
              fi
              ) 200>>${tslistll}.lock
              break
           else
              printf "     Pt%2s: TS optimized but not added-->imag=%4si cm-1, which is lower than %4si cm-1\n" $npo $fi $freqmin
           fi
        done
  done
done



