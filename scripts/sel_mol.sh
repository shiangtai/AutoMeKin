#!/bin/bash
echo ""
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
echo "Selecting starting structure in sel_mol"
###EMN. seed for srand. 
seed=$(basename $PWD | awk '/batch/{print $0}' | sed 's@batch@@' | awk '{print $1+100}')
###EMN
inputfile=$1
mm=$2
rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
energy=` awk 'BEGIN{e=0};/EKMC/{e=$2};END{print e}'  $inputfile `
en=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile `
molecule=` awk '/molecule/{print $2}'  $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile` 
#tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirHL_'$molecule'"}' $inputfile` 
#kmcfilehl=$tsdirhl/KMC/RXNet_long.cg_groupedprods
#minfilehl=$tsdirhl/MINs/SORTED/MINlist_sorted
#mindirhl=$tsdirhl/MINs/SORTED
kmcfilell=$tsdirll/KMC/RXNet_long.cg_groupedprods
minfilell=$tsdirll/MINs/SORTED/MINlist_sorted
mindirll=$tsdirll/MINs/SORTED
confilell=$tsdirll/working/conf_isomer.out
#confilehl=$tsdirhl/working/conf_isomer.out
factor=1.5
#if [ -f $minfilehl ] && [ -f $kmcfilehl ]; then
#   echo "itsscds with minima from tsdirhl"
#   if [ $2 -gt 0 ]; then
#      minn=$2
#   else 
#      minn=`awk '/min0/{print $2}' $minfilehl`
#   fi
###
#   minok=`awk 'BEGIN{min='$minn'}
#   {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
#   j=1
#   while(j<=iso[NR]){
#      if('$minn'==m[NR,j]) min=m[NR,1]
#      j++
#      }
#   }
#   END{print min}' $confilehl `
###
#   get_minn.sh $kmcfilehl $minok $en $factor
#   selm=`awk 'BEGIN{srand('$seed');rn=rand()};{n[NR]=$1};
#   END{
#   for(i=1;i<=NR;i++) den+='$factor'^(NR-i)
#   p[NR]=1/den
#   i=1
#   ptot=0
#   while(i<=NR){
#     p[i]='$factor'^(NR-i)*p[NR]
#     ptot+=p[i]
#     if(rn<ptot) {print n[i];exit}
#     i++ 
#     }
#   }' minn`
#   echo "selected minimum: $selm" 
#   names="MIN"$selm
#   sqlite3 $mindirhl/minshl.db "select natom,geom from minshl where name='$names'" | sed 's@|@\n\n@' > ${molecule}.xyz
#   sqlite3 $mindirhl/minshl.db "select geom from minshl where name='$names'" 
#   awk '{if(NR==1) print $0"\n"};{if(NF==4) print $0}' $mindirhl/MIN$selm"_"*.rxyz > $molecule.xyz
#elif [ -f $minfilell ] && [ -f $kmcfilell ]; then
if [ -f $minfilell ] && [ -f $kmcfilell ] && [ $mm -eq 1 ]; then
   echo "itsscds with minima from tsdirll"
#   if [ $2 -gt 0 ]; then
#      minn=$2
#   else
      minn=`awk '/min0/{print $2}' $minfilell`
#   fi
##
   minok=`awk 'BEGIN{min='$minn'}
   {for(i=1;i<=NF;i++) {m[NR,i]=$i;iso[NR]=NF}
   j=1
   while(j<=iso[NR]){
      if('$minn'==m[NR,j]) min=m[NR,1]
      j++
      }
   }
   END{print min}' $confilell `
##
###
   get_minn.sh $kmcfilell $minok $en $factor
   selm=`awk 'BEGIN{srand('$seed');rn=rand()};{n[NR]=$1};
   END{
   for(i=1;i<=NR;i++) den+='$factor'^(NR-i)
   p[NR]=1/den
   i=1
   ptot=0
   while(i<=NR){
     p[i]='$factor'^(NR-i)*p[NR]
     ptot+=p[i]
     if(rn<ptot) {print n[i];exit}
     i++
     }
   }' minn`
   echo "selected minimum: $selm" 
   names="MIN"$selm
   sqlite3 $mindirll/mins.db "select natom,geom from mins where name='$names'" | sed 's@|@\n\n@g' > ${molecule}.xyz
   sqlite3 $mindirll/mins.db "select geom from mins where name='$names'" 
#   awk '{if(NR==1) print $0"\n"};{if(NF==4) print $0}' $mindirll/MIN$selm"_"*.rxyz > $molecule.xyz
else
   echo "tsscds with only one minimum"
   cp $molecule".xyz" $molecule"_ref.xyz"
fi
echo "Exiting sel_mol"
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
echo ""
