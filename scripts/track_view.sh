#!/bin/bash
if [ -f tsscds.dat ];then
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
exe=$(basename $0)
cwd=$PWD
molecule=`awk '{if($1=="molecule") print $2}' $inputfile`
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirLL_'$molecule'"}' $inputfile`
sampling=` awk '{if($1=="sampling") {if($2=="microcanonical") print "1";if($2=="canonical") print "2";if($2=="association") print "3";if($2=="external") print "4"}}'  $inputfile `
if [ $sampling -eq 1 ]; then 
   flag="Energy"
   fwor="     Energy"
elif [ $sampling -eq 2 ]; then 
   flag="Temperature" 
   fwor="Temperature" 
fi
nts=$(awk 'END{print NR}' ${tsdirll}/tslist)
###view sqlite dabase
if [ -f ${tsdirll}/track.db ]; then
   id=$(sqlite3 ${tsdirll}/track.db "select max(id) from track")
   noj1=$(sqlite3 ${tsdirll}/track.db "select noj1 from track where id='$id'")
   nojf=$(sqlite3 ${tsdirll}/track.db "select nojf from track where id='$id'")
   batch=$(find . -maxdepth 1 -type d -print | grep 'batch')
   if [ ! -z "$batch" ]; then
      ntraj=$(for i in $(seq $noj1 $nojf); do grep "Trajectory" batch$i/tsscds.log 2>/dev/null ; done | wc -l)
      emin=$(sqlite3 ${tsdirll}/track.db "select emin from track where id='$id'")
      emax=$(sqlite3 ${tsdirll}/track.db "select emax from track where id='$id'")
###
      permin=$(for i in $(seq $noj1 $nojf); do awk '{if($1=="'"$flag"'") e=$4};/Npath/{++npath};/Trajectory/{++ntraj};END{diff=sqrt((e-'$emin')^2);if(diff<=1) print npath,ntraj}' batch$i/tsscds.log 2>/dev/null ; done | awk '{sumi+=$1;sumt+=$2};END{if(sumt>0) print int(sumi/sumt*100);if(sumt==0) print "-1"}')
      permax=$(for i in $(seq $noj1 $nojf); do awk '{if($1=="'"$flag"'") e=$4};/Npath/{++npath};/Trajectory/{++ntraj};END{diff=sqrt((e-'$emax')^2);if(diff<=1) print npath,ntraj}' batch$i/tsscds.log 2>/dev/null ; done | awk '{sumi+=$1;sumt+=$2};END{if(sumt>0) print int(sumi/sumt*100);if(sumt==0) print "-1"}')
###
      sqlite3 ${tsdirll}/track.db "update track set ntraj='$ntraj' where id='$id';update track set nts='$nts' where id='$id';update track set permin='$permin' where id='$id'; update track set permax='$permax' where id='$id'"
   fi
   printf "Iter #          TSs       ntrajs   $fwor range   Percent@Min   Percent@Max\n" 
   sqlite3 ${tsdirll}/track.db "select id,nts,ntraj,emin,emax,permin,permax from track" | sed 's@|@ @g' | awk '{ntraj+=$3;printf "%6.0f       %6.0f       %6.0f  %10.2f-%-10.2f     %6.0f        %6.0f\n",$1,$2,ntraj,$4,$5,$6,$7}' 
else
   echo "${tsdirll}/track.db does not exist"
   echo "This option is only available when running tsscds_parallel"
fi
###
