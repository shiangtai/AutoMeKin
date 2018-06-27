#!/bin/bash

source utils.sh


n=0
inputfile=$1
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
specA=`awk '/A=/{print $2}' $inputfile`
specB=`awk '/B=/{print $2}' $inputfile`
assocdir=$PWD"/assoc_"$specA"_"$specB
assoclist="$assocdir/assoclist"

#On exit remove tmp files
tmp_files=($assocdir/*.mop $assocdir/*.arc $assocdir/*.den $assocdir/*.res)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
ls $assocdir/assoc*.out > $assoclist
avgerr=`awk '/avgerr/{print $2}' $inputfile`
bigerr=`awk '/bigerr/{print $2}' $inputfile`
thdiss=`awk '/thdiss/{print $2}' $inputfile`
if [ -f "black_list.dat" ]; then rm black_list.dat; fi
if [ -f "black_list.out" ]; then rm black_list.out; fi
echo "Screening" > $assocdir/assoclist_screened
echo "List of errors in the optimizations" > $assocdir/assoclist_opt_err
echo "List of disconnected (at least two fragments) structures" > $assocdir/assoclist_disconnected
set `awk '{print $1}' $assoclist`
for i
do
   file=$i
   ((n=n+1))
   name=`basename $file .out`
   echo $file  $name
   get_geom_mopac.sh $file > tmp_geom
   cherr=`awk 'BEGIN{zero=0.0};/Error/{err=1};END{if(err==0) err=zero;print err}'  tmp_geom`
   if  [[ ("$cherr" -eq "1") ]] 
   then 
     echo "Error in this output file"
     echo $name" Error">> $assocdir/assoclist_screened
     echo $name" Error">> $assocdir/assoclist_opt_err
### mv ouput files with errors
     mv $file $assocdir/OPTERR_$name
###
     continue 
   else
     echo $name" data">> $assocdir/assoclist_screened
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

   natom=`wc -l mingeom| awk '{print $1}' `
   echo "1 $natom" | cat - ConnMat | sprint.exe >sprint.out
 
   paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out

   deg_form.sh > deg_form.out
   awk '/HEAT OF FORMATION =/{e=$5};END{printf "%8.2f\n",e}' $file > $assocdir/$name"_data"
   format.sh $name $assocdir $thdiss 
   ndis=`awk '{ndis=$1};END{print ndis}' $assocdir/$name"_data" `
### mv structures where there is 2 or more fragments already formed
   if  [[ ("$ndis" -gt "1") ]] 
   then 
     mv $file $assocdir/DISCNT_$name
   fi
###
   cat $assocdir/$name"_data" >> $assocdir/assoclist_screened 
done
rm $assocdir/*_data
reduce.sh $assocdir assoc
awk '{if($NF==1) print $0}' $assocdir/assoclist_screened.red > $assocdir/assoclist_screened.redconn
awk '{if($NF>1) print $0}' $assocdir/assoclist_screened.red >> $assocdir/assoclist_disconnected
diffGT.sh $assocdir/assoclist_screened.redconn $assocdir assoc $avgerr $bigerr
### mv repeated structures
if [ -f "black_list.out" ]; then
   set `awk '{print $0}' black_list.out `
   for i
   do 
     mv $assocdir/assoc$i.out $assocdir/REPEAT_assoc$i.out
   done
else
   echo "No repetitions"
fi
###

sort_assoc.sh $inputfile

select_assoc.sh $inputfile
###
