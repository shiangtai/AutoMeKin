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
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "tsdirHL_'$molecule'"}' $inputfile`

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


awk '/failed/{print $3}' $tsdirhl/KMC/RXNet_long.cg > listtss
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

minpath=$tsdirhl/MINs/SORTED

#Do the calcs for the minima only once

totmin=`awk 'END{print $2}' $minpath/MINlist_sorted`
#set `awk '{print $2}' $minpath/MINlist_sorted`
#for i
for file in $(sqlite3 ${minpath}/minshl.db "select lname from minshl")
do
   i=$(sqlite3 ${minpath}/minshl.db "select id from minshl where lname='$file'")
#   file=$minpath/MIN$i"_"*.rxyz
   echo $i "out of" $totmin "file=" $file
   
   echo "Labels" >labels
   sqlite3 ${minpath}/minshl.db "select geom from minshl where lname='$file'"  >mingeom
   awk '{print $1}' mingeom  >>labels
#   awk '{if(NF==4) print $0 }' $file > mingeom
#   awk '{if(NF==4) print $1 }' $file >> labels

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

   cat tmp_atn tmp_adja >tmp_wrk

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
   print '$i',sum}'  tmp_wrk >> tmp_values_ref


done


sed 's/.rxyz//g' minfailed_list >tmp01


set `awk '{print $0} ' listtss`
for i
do 
  ((n=n+1))
  name=`awk '/'$i'/{print $0} ' minfailed_list`
  nmbr=`awk '/'$i'/{print NR} ' minfailed_list`
  if [ "$nmbr" == "" ]; then
     echo "the structure that failed is not in the list of failed opts. Check the IRC"
     if [ $vdw -eq 1 ]; then
        cp thdist_backup thdist
     fi
     exit
  fi 
  name2=`awk 'NR=='$nmbr',NR=='$nmbr'{print $0} ' tmp01`
  echo  $i $name $nmbr $name2
  awk '{if( NR == FNR) l[NR]=$1} 
  /Coordinates/{getline
  getline
  i=1
  natom=0
  while(i<=10000){
    getline
    if(NF==6)  ++natom
    if(NF==6)  nl[natom]=$2
    if(NF==6)  x[natom]=$4
    if(NF==6)  y[natom]=$5
    if(NF==6)  z[natom]=$6
    if(NF==1)  break
    ++i
    }
  }
  END{
  n=1
  while(n<=natom){
   print l[nl[n]],x[n],y[n],z[n]
   ++n
   }
  }' atsdum2.out $tsdirhl/IRC/$name2".log"  > minfailed_structure

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


   cat tmp_atn tmp_adja >tmp_wrk
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
   print sum}'  tmp_wrk > tmp_value_structure

cat tmp_value_structure tmp_values_ref | awk 'BEGIN{print "This structure is similar to MIN"} 
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


