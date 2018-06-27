#!/bin/bash
sharedir=${TSSCDS}/share
source utils.sh
#On exit remove tmp files
tmp_files=(ConnMat deg* fort.* intern* labels mingeom minn ScatMat ts_tors* ScalMat *_opt.* tors.arc tors.mop tors.out tors.xyz geom*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
cwd=$PWD
if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
mole=` awk '{if($1=="molecule") print $2}'  $inputfile `
molecule=tors
nb=`echo $cwd | awk 'BEGIN{FS="/"};{print $NF}' `
###
wrkmode=`awk 'BEGIN{mode=0};/fastmode/{mode=1};END{print mode}' $inputfile `
charge=` awk '/charge/{print $2}'  $inputfile `
method=` awk 'BEGIN{llcalc="PM7"};{if($1=="LowLevel") {$1="";llcalc=$0}};END{print llcalc}' $inputfile `
if [ $wrkmode -eq 0 ]; then ts_template="$(cat $sharedir/ts_templateslow)" ; fi
if [ $wrkmode -eq 1 ]; then ts_template="$(cat $sharedir/ts_templatefast)" ; fi
dumoldgeofreq="$(cat $sharedir/freq_template2 | sed 's/method/'"$method"' charge='$charge'/g')"
dum_ts="$(echo "$ts_template" | sed 's/method/'"$method"' charge='$charge'/g')"
dum1="$(sed 's/method/'"$method"' charge='$charge' BONDS INT/g' $sharedir/freq_template1)"

###

tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$mole'"}' $inputfile` 
rcwd=$(echo $tsdirll | sed 's@/tsdir@ @' | awk '{print $1}')
rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
energy=` awk 'BEGIN{e=0};/EKMC/{e=$2};END{print e}'  $inputfile `
en=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile `
emaxts=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print 1.5*en}' $inputfile `

tslistll=${tsdirll}/tslist
freqmin=` awk '/freqmin/{print $2}'  $inputfile `
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
#Compute the absolute values of emax and emin 
if [ -f $tsdirll/MINs/min.db ]; then
   e0=$(sqlite3 ${tsdirll}/MINs/min.db "select energy from min where name='min0_0'")
else
   echo "$dum1"                         > ${mole}_opt.mop
   awk 'NR>2{print $0}' ${mole}.xyz    >> ${mole}_opt.mop
   echo "$dumoldgeofreq"               >> ${mole}_opt.mop
   mopacl ${mole}_opt.mop 2>/dev/null
   e0=$(awk '/FINAL HEAT OF FORMATION =/{e0=$6};END{print e0}' ${mole}_opt.out )
fi

echo The value of e0=$e0
echo The relative value of emaxts=$emaxts
edum=`echo "$emaxts+$e0" | bc | awk '{printf "%10.0f",$1}'`
emaxts=$edum
echo The absolute value of emaxts=$emaxts

##Print the absolute values of emax and emin 
#echo "Emaxts (absolute value) is $emaxts"

nmin=0

kmcfilell=$tsdirll/KMC/RXNet_long.cg_groupedprods
minfilell=$tsdirll/MINs/SORTED/MINlist_sorted
confilell=$tsdirll/working/conf_isomer.out
factor=1.5

if [ -f $minfilell ] && [ -f $kmcfilell ]; then
   echo "Select the relevant minima to find conformational isomers"
   minn=`awk '/min0/{print $2}' $minfilell`
   if [ ! -f $confilell ]; then
      echo "File $confilell does not exist and we cannot proceed"
      exit
   fi
   minok=`awk 'BEGIN{min='$minn'}
   {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
   j=1
   while(j<=iso[NR]){
      if('$minn'==m[NR,j]) min=m[NR,1]
      j++
      }
   }
   END{print min}' $confilell `
else
   echo "Files $minfilell and/or $kmcfilell do not exist and we cannot proceed"
   exit
fi

get_minn.sh $kmcfilell $minok $en $factor

set `awk '{print NR}' minn `
for jmin
do
   ((nmin=nmin+1))
   minnb=`awk 'NR=='$nmin',NR=='$nmin'{print $1}' minn`
   echo "Running the tors calculations for  MIN $minnb"
   names="MIN"$minnb
   sqlite3 $tsdirll/MINs/SORTED/mins.db "select geom from mins where name='$names'" >${molecule}.xyz
   natom=` awk 'END{print NR}'  ${molecule}.xyz `
   echo "$dum1"            > $molecule"_opt.mop"
   cat ${molecule}.xyz	  >> $molecule"_opt.mop"
   echo "$dumoldgeofreq"  >> $molecule"_opt.mop"
   mopacl $molecule"_opt.mop" 2>/dev/null
   get_geom_mopac.sh $molecule"_opt.out" | awk '{if(NF==4) print $0}' >mingeom 
   createthdist.sh $thd
   createMat.sh
   deg="$(awk '{si=0;for(i=1;i<=NF;i++) {si+=$i;if(i==NF) print si,$0 }}' ConnMat)" 
   bo="$(awk '/BOND ORDERS AND VALENCIES/{getline
   getline
   getline
   getline
   i=1
   while(i<='$natom'){
       getline
       for(j=3;j<=i+1;j++) {k=j-2;if($j<2.0) print i,k}
       i++
     }
   }' $molecule"_opt.out")" 
   echo "$deg" > deg_bo
   echo "$bo" >> deg_bo
   awk '{
   if(NR<='$natom'){
     deg[NR]=$1
     l=0
     if(deg[NR]>1) {for(i=2;i<=NF;i++) {if($i==1) {++l;jatom[NR,l]=i-1 } }}
     }
   else {
     if(deg[$1]>1 && deg[$2]>1) {
       ok=0
       j=1
       while( j<=deg[$1] ){
         if(jatom[$1,j] != $2) k=jatom[$1,j]
         else ok=1
         ++j
         }
       j=1
       while( j<=deg[$2] ){
         if(jatom[$2,j] != $1) {l=jatom[$2,j];break}
         ++j
         }
       if(ok==1) print k,$0,l
       }
     }
   }' deg_bo  > deg_bo.out
   
   ntor=`wc -l deg_bo.out | awk '{print $1}'`
   echo "Number of torsions $ntor"
   if [ $ntor -eq 0 ]; then exit ; fi
   nol=0
   set `awk '{print NR}' deg_bo.out `
   for itor
   do
      ((nol=nol+1))
      echo "Running the TS search for tors $itor of molecule $molecule"
      lr=`awk 'NR=='$nol',NR=='$nol'{print $1;exit}' deg_bo.out `
      l1=`awk 'NR=='$nol',NR=='$nol'{print $2;exit}' deg_bo.out `
      l2=`awk 'NR=='$nol',NR=='$nol'{print $3;exit}' deg_bo.out `
      l3=`awk 'NR=='$nol',NR=='$nol'{print $4;exit}' deg_bo.out `
      if [ $l1 -gt $lr ]; then ((l1=l1-1)) ; fi
      if [ $l2 -gt $lr ]; then ((l2=l2-1)) ; fi
      if [ $l3 -gt $lr ]; then ((l3=l3-1)) ; fi
      get_geom_mopac.sh $molecule"_opt.out" >intern.dat
      awk 'NR=='$nol',NR=='$nol'{print $0}' deg_bo.out >>intern.dat
      intern.exe <intern.dat>intern.out
      fok=`awk 'BEGIN{fok=1};/Abort/{fok=0};END{print fok}' intern.out`
      if [ $fok -eq 0 ]; then continue ; fi
      internlastatom="$(awk '{print $1,$2,$3,$4,$5,$6,$7,'$l1','$l2','$l3'}' intern.out)"
      sed 's/method/'"$method"' charge='$charge'/g' $sharedir/path_template >tors.mop
      get_geom_mopac.sh $molecule"_opt.out" | awk '{if(NF==4 ) print $0}' | awk '{if(NR!= '$lr') print $0}' >> tors.mop 
      echo "$internlastatom" >> tors.mop
      mopacl tors.mop 2>/dev/null
      awk '/VARIABLE        FUNCTION/{++i;getline
      e[i]=$3
      getline
      getline
      getline
      j=1
      while(j<='$natom'){
        getline
        l[i,j]=$0
        j++
        }
      }
      END{
      i0=2
      while(i0<=i){
        if(e[i0]>e[i0-1] && e[i0]>e[i0+1]) {
           k=1
           proc=1
           while(k<=nmax){
             diff0=e[i0]-emax[k]
             diff=sqrt(diff0*diff0)
             if(diff<=0.01) proc=0
             k++
             } 
           if(proc==0) {i0++;continue}
           ++nmax
           emax[nmax]=e[i0]  
           j=1
           while(j<='$natom'){
             print l[i0,j]
             j++
             }
           }
        i0++
        } 
      }' tors.out>geomts_tors0
      nomax="$(awk 'END{print NR/'$natom'}' geomts_tors0)"
      echo "$nomax possible candidate(s)"
      for inmax in $(seq 1 $nomax)
      do
         echo "$inmax th candidate(s)"
         awk 'BEGIN{nr0='$natom'*('$inmax'-1)+1;nrf='$natom'*'$inmax'};{if(NR>=nr0 && NR<=nrf) print $0}' geomts_tors0 >geomts_tors
         awk '$7=1' geomts_tors >geomts_torsout
         echo "$dum_ts"          > ts_tors$itor".mop"
         paste geomts_torsout   >> ts_tors$itor".mop"
         echo "$dumoldgeofreq"  >> ts_tors$itor".mop"
         mopacl ts_tors$itor".mop" 2>/dev/null

   #check the ts
         fe="$(mopac_freq_ts.sh ts_tors$itor".out")" 
         fi="$(echo "$fe" | awk '{printf "%10.0f",$1}')"
         ei="$(echo "$fe" | awk '{printf "%10.0f",$2}')"
         if [[ ("$fi" -eq -1) ]]; then
            echo "TS opt failed. Lowest real freq < 0 cm-1"
            continue
         elif [[ ("$fi" -eq -2) ]]; then
            echo "TS opt failed. Sum of two Lowest real freqs < 50 cm-1"
            continue
         elif [[ ("$fi" -eq -3) ]]; then
            echo "TS opt failed. Stationary point is a minimum"
            continue
         elif [[ ("$fi" -eq -4) ]]; then
            echo "TS opt failed. Stationary point search crashed"
            continue
         elif [[ ("$fi" -eq -5) ]]; then
            echo "TS opt failed. This is a vdw ts (separated products)"
            continue
         elif [[ ("$ei" -gt "$emaxts") ]]; then
            echo "TS has an energy $ei, which is greater than" $emaxts
            continue
         fi
         if [[ ("$fi" -ge "$freqmin") ]]; then
            f="$(echo "$fe" | awk '{printf "%10.2f",$1}')"
            e="$(echo "$fe" | awk '{printf "%10.2f",$2}')"
            f1="$(echo "$fe" | awk '{printf "%10.2f",$3}')"
            f2="$(echo "$fe" | awk '{printf "%10.2f",$4}')"
            f3="$(echo "$fe" | awk '{printf "%10.2f",$5}')"
            f4="$(echo "$fe" | awk '{printf "%10.2f",$6}')"
            if [ -f "$tslistll" ]; then
               ok=`diff.sh $f $e $f1 $f2 $f3 $f4 $tslistll`
               if [[ ("$ok" -eq "-1") ]]; then
                  nt=`awk '{nt=$2};END{print nt}' $tslistll `
                  ((nt = nt + 1))
                  name=ts${nt}_${nb}
                  printf "ts%5s%18s%9s%9s%9s%9s%9s%9s traj=   0 Path=%10s\n" $nt $name $f $e $f1 $f2 $f3 $f4 $i $nb  >> $tslistll
                  cp ts_tors${itor}.out  $tsdirll/${name}.out 
                  printf "     Pt%2s: TS optimized and added to ts list\n" $itor
                  continue
               else
                  printf "     Pt%2s: TS optimized but not added-->redundant with ts %4s\n" $itor $ok
                  continue
               fi
            else
               nt=1
               name=ts${nt}_${nb}
               cp ts_tors${itor}.out  $tsdirll/${name}.out
               printf "     Pt%2s: TS optimized and added to ts list\n" $itor
               printf "ts%5s%18s%9s%9s%9s%9s%9s%9s traj=    0 Path=%10s\n" $nt $name $f $e $f1 $f2 $f3 $f4 $i $nb  >> $tslistll
               continue
            fi
         else
            printf "     Pt%2s: TS optimized but not added-->imag=%4si cm-1, which is lower than %4si cm-1\n" $itor $fi $freqmin
         fi
      done
##
   done
done

