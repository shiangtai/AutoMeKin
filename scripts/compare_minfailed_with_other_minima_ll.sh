#!/bin/bash
#this script serves to identify the failed minima (those that remain even after using MINFAILED.sh)
# it compares the structure of the failed minimum with those of the MIN already optimized
sharedir=${AMK}/share
source utils.sh
#remove tmp files
tmp_files=(tmp* deg* ConnMat atsdum2.out labels mingeom ScalMat sprint.* listtss)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)

inputfile=amk.dat
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "tsdirLL_'$molecule'"}' $inputfile`

thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
vdw=`awk 'BEGIN{vdw=0};{if($1=="vdw") vdw=1};END{print vdw}' $inputfile `
if [ $vdw -eq 1 ]; then
   cp thdist thdist_backup
   nvdw=`awk '{if($1=="vdw") print $2}' $inputfile `
   echo $nvdw "Van der Waals distances to be taken into account"
   for i in $(seq 1 $nvdw)
   do
      at1[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $1; i++}  }' $inputfile `
      at2[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i')print $2; i++}  }' $inputfile `
      dis[$i]=`awk '/vdw/{i=1; while(i<='$i'){getline;if(i=='$i') print $3; i++}  }' $inputfile `
      echo "Distance $i between atoms ${at1[$i]} and ${at2[$i]} is ${dis[$i]}"
      awk '{if($1=="'${at1[$i]}'" && $2=="'${at2[$i]}'") {print $1,$2,"'${dis[$i]}'"}
      else if($2=="'${at1[$i]}'" && $1=="'${at2[$i]}'") {print $1,$2,"'${dis[$i]}'"}
      else {print $0}
      } ' thdist >thdist_vdw
      cp thdist_vdw thdist
   done
fi


awk '/failed/{print $3}' $tsdirll/KMC/RXNet_long.cg > listtss
nfail=`wc -l listtss | awk '{print $1}'`
if [ $nfail -eq 0 ]; then
   echo "You are lucky. No problems with the structures. Go head"
   if [ $vdw -eq 1 ]; then
      cp thdist_backup thdist
   fi
   exit
else
   echo "For" $nfail "paths the minima are problematic" 
   echo "You should edit KMC/RXNet_long.cg according to the suggestions below"
fi

sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt} 
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > atsdum2.out

echo "This file is just for the values of the computed minima" >tmp_values_ref

minpath=$tsdirll/MINs/SORTED

#Do the calcs for the minima only once

totmin=`awk 'END{print $2}' $minpath/MINlist_sorted`
for file in $(sqlite3 ${minpath}/mins.db "select lname from mins")
do
   i=$(sqlite3 ${minpath}/mins.db "select id from mins where lname='$file'")
   echo $i "out of" $totmin "file=" $file
   
   echo "Labels" >labels
   sqlite3 ${minpath}/mins.db "select geom from mins where lname='$file'"  >mingeom
   awk '{print $1}' mingeom  >>labels

   if [ $vdw -eq 0 ]; then
      createthdist.sh $thd
   fi

   natom=`wc -l mingeom| awk '{print $1}' `
   echo "1" $natom > sprint.dat
   createMat.sh
   cat ConnMat >> sprint.dat
   sprint.exe <sprint.dat>sprint.out
   awk '/Natom/{natom=$2}
   /Adjace/{i=1
   while(i<=natom){
     getline
     print $0
     i++
     }
   }' sprint.out >tmp_adja


   awk '{if( NR == FNR) {l[NR]=$1;n[NR]=NR/10+1;tne=NR}}
   {if(FNR > 1 && NR >FNR ) {
      IGNORECASE = 1
      i=1
      while(i<=tne){
         if( $1 == l[i]) print n[i]
         i++
         }
     }
   }' atsdum2.out labels > tmp_atn

   cat tmp_atn tmp_adja >wrk

   awk '{if(NF==1) n[NR]=$1}
   {
   if(NF>1) {++i; for(j=1;j<=NF;j++) a[i,j]=$j}
   }
   END{
   ii=1
   while(ii<=i) {
     sumi=0
     j=1
     while(j<=i) {
        k=1
        while(k<=i) {
           sumi+=a[ii,j]*a[j,k]*n[ii]*n[j]*n[k]
           ++k
           }
        ++j
        }
     sum+=(sumi)^n[ii]
     ++ii
     }
   print '$i',sum}'  wrk >> tmp_values_ref


done


sed 's/.rxyz//g' $tsdirll/MINs/minfailed_list > tmp01


set `awk '{print $0} ' listtss`
for i
do 
  ((n=n+1))
  name=`awk '/'$i'/{print $0} ' $tsdirll/MINs/minfailed_list`
  nmbr=`awk '/'$i'/{print NR} ' $tsdirll/MINs/minfailed_list`
  if [ "$nmbr" == "" ]; then
     echo "the structure that failed is not in the list of failed opts. Check the IRC"
     if [ $vdw -eq 1 ]; then
        cp thdist_backup thdist
     fi
     exit
  fi 
  name2=`awk 'NR=='$nmbr',NR=='$nmbr'{print $0} ' tmp01`
  sqlite3 $tsdirll/MINs/min.db "select geom from min where name='$name2'"  > minfailed_structure

# do the calc for the minfailed_structure
   echo "Labels" >labels
   cp minfailed_structure  mingeom
   awk '{print $1 }' minfailed_structure >> labels
   if [ $vdw -eq 0 ]; then
      createthdist.sh $thd
   fi
   natom=`wc -l mingeom| awk '{print $1}' `
   echo "1" $natom > sprint.dat
   createMat.sh
   cat ConnMat >> sprint.dat
   sprint.exe <sprint.dat>sprint.out
   awk '/Natom/{natom=$2}
   /Adjace/{i=1
   while(i<=natom){
     getline
     print $0
     i++
     }
   }' sprint.out >tmp_adja

   awk '{if( NR == FNR) {l[NR]=$1;n[NR]=NR/10+1;tne=NR}}
   {if(FNR > 1 && NR >FNR ) {
      IGNORECASE = 1
      i=1
      while(i<=tne){
         if( $1 == l[i]) print n[i]
         i++
         }
     }
   }' atsdum2.out labels > tmp_atn


   cat tmp_atn tmp_adja >wrk
   awk '{if(NF==1) n[NR]=$1}
   {
   if(NF>1) {++i; for(j=1;j<=NF;j++) a[i,j]=$j}
   }
   END{
   ii=1
   while(ii<=i) {
     sumi=0
     j=1
     while(j<=i) {
        k=1
        while(k<=i) {
           sumi+=a[ii,j]*a[j,k]*n[ii]*n[j]*n[k]
           ++k
           }
        ++j
        }
     sum+=(sumi)^n[ii]
     ++ii
     }
   print sum}'  wrk > tmp_value_structure

cat tmp_value_structure tmp_values_ref |  awk 'BEGIN{print "Structure '$name2' is similar to MIN"}
{if(NF==1) comp=$1
if(NF==2 && comp==$2) printf "%s ", $1
}
END{
print "\n"
}' 

done

if [ $vdw -eq 1 ]; then
   cp thdist_backup thdist
fi


