#!/bin/bash
source utils.sh

inputfile=$1
molecule=` awk '/molecule/{print $2}'  $inputfile `
natom="$(awk 'NR==1,NR==1{print $1}' $molecule.xyz)"
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`
#On exit remove tmp files
tmp_files=(${tsdirll}/tmp* ${tsdirll}/tslist_* tmp* black_list*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT
exe=$(basename $0)
tslistll=${tsdirll}/tslist
tslistlog=${tsdirll}/tslistlog
bu_ts=${tsdirll}/backup
avgerr=`awk '/avgerr/{print $2}' $inputfile`
bigerr=`awk '/bigerr/{print $2}' $inputfile`
thdiss=`awk '/thdiss/{print $2}' $inputfile`
##
if [ ! -d "$bu_ts" ]; then
   echo "Making a backup folder and saving tslist"
   mkdir $bu_ts 
   cp $tslistll $bu_ts
else
   echo "Removing backup folder content and saving tslist"
   rm $bu_ts/*
   cp $tslistll $bu_ts
fi

nrep=0
nerr=0
nfra=0
if [ -f screening.log ]; then
   echo "" >>screening.log
   nu=$(awk '/Summary/{nu=$6};END{print nu+1}' screening.log)
   echo "Summary of screening calculations. Execution $nu of $exe" >>screening.log
else
   nu=1
   echo "Summary of screening calculations. Execution $nu of $exe" >screening.log
fi

if [ -f "black_list.dat" ]; then rm black_list.dat; fi
if [ -f "black_list.out" ]; then rm black_list.out; fi
echo "Screening" > $tsdirll/tslist_screened
echo "List of disconnected (at least two fragments) TS structures" > $tsdirll/tslist_disconnected
for name in $(awk '{print $3}' $tslistll)
do
   file=${tsdirll}/${name}.out
   cp $file $bu_ts
   i=$(echo $name | sed 's@ts@@;s@_@ @' | awk '{print $1}')
   get_geom_mopac.sh $file > tmp_geom
   cherr=`awk 'BEGIN{zero=0.0};/Error/{err=1};END{if(err==0) err=zero;print err}'  tmp_geom`
   if  [[ ("$cherr" -eq "1") ]] 
   then 
     ((nerr=nerr+1))
     echo "Structure ts$i removed-->opt failed"
     echo "Structure ts$i removed-->opt failed" >>screening.log
     echo ts$i"_out Error">> $tsdirll/tslist_screened
### remove ouput files with errors
     rm -rf $file 
     sed -i '/'$name'/d' $tslistll
###
     continue 
   else
     echo ts$i"_out data">> $tsdirll/tslist_screened
   fi 
   echo "Labels" >labels
   awk '
   NR==1{natom=$1
   i=1
   getline
   while(i<=natom){
     getline
     print $0
     print $1 >> "labels"
     i++ 
     }
   }' tmp_geom >mingeom
   createthdist.sh $thd
   createMat2.sh 
   echo "1 $natom" | cat - ConnMat | sprint.exe >sprint.out

   paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out

   deg_form.sh > deg_form.out
   awk '/HEAT OF FORMATION =/{e=$5};END{printf "%9.3f\n",e}' $file > $tsdirll/ts$i"_data"
   format.sh ts$i $tsdirll $thdiss 
   ndis=`awk '{ndis=$1};END{print ndis}' $tsdirll/ts$i"_data" `
### mv TSs where there is 2 or more fragments already formed
   if  [[ ("$ndis" -gt "1") ]] 
   then 
     ((nfra=nfra+1))
     mv $file $tsdirll/DISCNT_ts$i.out
     els="$(cat tmp_ELs)"
     printf "Structure ts%-5s renamed as DISCNT_ts%-5s-->%2s fragments. Values of thdiss: %-40s\n" $i $i $ndis "$els"
     printf "Structure ts%-5s renamed as DISCNT_ts%-5s-->%2s fragments. Values of thdiss: %-40s\n" $i $i $ndis "$els" >>screening.log
     sed -i '/'$name'/d' $tslistll
###remove from sqlite database
   fi
###
   cat $tsdirll/ts$i"_data" >> $tsdirll/tslist_screened 
done
rm $tsdirll/ts*_data
reduce.sh $tsdirll ts
awk '{if($NF==1) print $0}' $tsdirll/tslist_screened.red > $tsdirll/tslist_screened.redconn
awk '{if($NF>1) print $0}' $tsdirll/tslist_screened.red >> $tsdirll/tslist_disconnected
diffGT.sh $tsdirll/tslist_screened.redconn $tsdirll ts $avgerr $bigerr
### mv repeated TSs
if [ -f "black_list.out" ]; then
   set `awk '{print $0}' black_list.out `
   for i
   do 
     file=$tsdirll/ts${i}_*.out
     name=$(basename $file .out)
     mv $file $tsdirll/REPEAT_ts$i.out
###
     ((nrep=nrep+1))
     orig=$(awk '{if($2=='$i') {print $1;exit}}' $tsdirll/tslist_screened.lowdiffs)
     apes="$(awk '{if($2=='$i') {print $3,$4;exit}}' $tsdirll/tslist_screened.lowdiffs)"
     printf "Structure ts%-5s renamed as REPEAT_ts%-5s-->redundant with ts%-5s. Values of avgerr and bigerr: %-40s\n" $i $i $orig "$apes"
     printf "Structure ts%-5s renamed as REPEAT_ts%-5s-->redundant with ts%-5s. Values of avgerr and bigerr: %-40s\n" $i $i $orig "$apes" >>screening.log
     sed -i '/'$name'/d' $tslistll
   done
fi
echo "$nrep repetitions"
echo "$nrep repetitions" >> screening.log
echo "$nfra fragmented"
echo "$nfra fragmented" >> screening.log
echo "$nerr opt failed"
echo "$nerr opt failed" >> screening.log
###
if [ -f ${tsdirll}/track.db ]; then
   track_view.sh 
fi


