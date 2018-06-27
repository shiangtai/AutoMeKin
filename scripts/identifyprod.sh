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
   echo "tsscds input file is missing"
   exit
fi
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`
rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
energy=` awk 'BEGIN{e=0};/EKMC/{e=$2};END{print e}'  $inputfile `
en=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile `

nlprlist=`wc -l $tsdirll/PRODs/PRlist | awk '{print $1}'`

if [ $nlprlist -eq 1 ];then
   echo "No products for this system"
   minn0=`awk '/min0/{print $2}' $tsdirll/MINs/SORTED/MINlist_sorted`
   minn=`awk 'BEGIN{minn='$minn0'};/ '$minn0' /{minn=$1};END{print minn}' $tsdirll/working/conf_isomer.out`
   linked_paths.sh $tsdirll/KMC/RXNet_long.cg $minn $en  > $tsdirll/KMC/RXNet_long.cg_groupedprods
#   cp tmp_lp $tsdirll/KMC/RXNet_long.cg_groupedprods
   exit
fi

echo "codes of products" > tmp_code
for name in $(sqlite3 $tsdirll/PRODs/prod.db "select name from prod")
#set ` awk '{if(NR>1) print $2}' $tsdirll/PRODs/PRlist `
#for i
do
#   file=$tsdirll/PRODs/PR$i"_"*
#   file=$tsdirll/PRODs/${name}.rxyz
#   sqlite3 $tsdirll/PRODs/prod.db "select natom,'E= XX ZPE= XX Gcorr XX',geom,freq from prod where name='$name'" | sed 's@|@\n@g' >tmp_forminp
#   formula="$(FormulaPROD.sh tmp_forminp)"
   formula="$(sqlite3 $tsdirll/PRODs/prod.db "select natom,'E= XX ZPE= XX Gcorr XX',geom,freq from prod where name='$name'" | sed 's@|@\n@g' | FormulaPROD.sh)"
   echo "$formula" >>tmp_code
#   sqlite3 ${tsdirll}/PRODs/prod.db "insert into prod (formula) values ('$tmp_code');"
   sqlite3 ${tsdirll}/PRODs/prod.db "update prod set formula='$formula' where name='$name';"
   echo "Getting the formula for $name"
done

paste $tsdirll/PRODs/PRlist tmp_code > $tsdirll/PRODs/PRlist_kmc

lastmin=`awk '{lm=$2};END{print lm}' $tsdirll/MINs/SORTED/MINlist_sorted `

sed 's/ + /+/g' $tsdirll/PRODs/PRlist_kmc | awk 'BEGIN{n='$lastmin'} 
{if(NR==1) print $0
if($1=="PROD") {++i
  fl[i]=$NF
  j=1
  p=1
  while(j<=i-1){
   if(fl[i]==fl[j]) {p=0;code[i]=code[j];break}
   j++
   }
  if(p==1) {++n;code[i]=n}
  print $1,$2,code[i]
  }
}' > $tsdirll/PRODs/PRlist_kmc.log


cat $tsdirll/PRODs/PRlist_kmc.log $tsdirll/KMC/RXNet_long.cg | awk '/PROD/{if(NF==3) ncode[$2]=$3} 
/KMC file/{lp=1;print $0}
{if(lp==1) {
  if($10~"MIN" || $1~"number") print $0
  if($10~"PROD") printf "%2s %4.0f %20s %3s %8.3f %5s %4s %3.0f %4s %6s %3.0f %5.0f %15.0f %6.0f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,ncode[$11],$12,$13,$14
  }
}' >tmp

minn0=`awk '/min0/{print $2}' $tsdirll/MINs/SORTED/MINlist_sorted`
minn=`awk 'BEGIN{minn='$minn0'};/ '$minn0' /{minn=$1};END{print minn}' $tsdirll/working/conf_isomer.out`

linked_paths.sh tmp $minn $en > $tsdirll/KMC/RXNet_long.cg_groupedprods

#cp tmp_lp $tsdirll/KMC/RXNet_long.cg_groupedprods
