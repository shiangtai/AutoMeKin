#!/bin/bash
source utils.sh
#remove tmp files
tmp_files=(tmp*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
inputfile=amk.dat
molecule=` awk '/molecule/{print $2}'  $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`
# toen=0 print DG
# toen=1 print DE
toen=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
if [ $toen -eq 0 ] ; then
   echo "Sort Gibbs Free Energy differences"
elif [ $toen -eq 1 ]; then
   echo "Sort Energy(+ZPE) differences"
else
   echo "Specify a type of rate: canonical or microcanonical"
   exit
fi

if [ $toen -eq 1 ]; then
   dum_min0=$(sqlite3 ${tsdirll}/MINs/norep/minnr.db "select energy,zpe from minnr where name like '%min0%'" | sed 's@|@ @g')
elif [ $toen -eq 0 ]; then
   dum_min0=$(sqlite3 ${tsdirll}/MINs/norep/minnr.db "select energy,g from minnr where name like '%min0%'" | sed 's@|@ @g')
fi
# First of all, we sort the TSs
if [ ! -d "$tsdirll/TSs/SORTED" ]; then
   echo "$tsdirll/TSs/SORTED does not exist. It will be created"
   mkdir $tsdirll/TSs/SORTED
else
   rm -r $tsdirll/TSs/SORTED
   mkdir $tsdirll/TSs/SORTED
   echo "$tsdirll/TSs/SORTED already exists but it was created again"
fi
if [ ! -d "$tsdirll/MINs/SORTED" ]; then
   echo "$tsdirll/MINs/SORTED does not exist. It will be created"
   mkdir $tsdirll/MINs/SORTED
else
   rm -r $tsdirll/MINs/SORTED
   mkdir $tsdirll/MINs/SORTED
   echo "$tsdirll/MINs/SORTED already exists but it was created again"
fi
###create sorted tables
sqlite3 ${tsdirll}/TSs/SORTED/tss.db "drop table if exists tss; create table tss (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,lname TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT);"
sqlite3 ${tsdirll}/MINs/SORTED/mins.db "drop table if exists mins; create table mins (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,lname TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT);"

# We start with the TSs
if [ -f tmp_oe ]; then rm tmp_oe ; fi
if [ -f tmp_noe ]; then rm tmp_noe ; fi


echo "non-ordered energies"
for name in $(sqlite3 ${tsdirll}/TSs/ts.db "select name from ts")
do
  echo $name
  namer=${name}.rxyz
  if [ $toen -eq 1 ]; then
     dum_en=$(sqlite3 ${tsdirll}/TSs/ts.db "select energy,zpe from ts where name='$name'" | sed 's@|@ @g')
  elif [ $toen -eq 0 ]; then
     dum_en=$(sqlite3 ${tsdirll}/TSs/ts.db "select energy,g from ts where name='$name'" | sed 's@|@ @g')
  fi
  echo "$dum_min0"   >tmp_en
  echo "$dum_en"    >>tmp_en
  awk '{e[NR]=$1;ezpe[NR]=$2};END{printf "%10s %10.6f\n","'$namer'",e[2]-e[1]+ezpe[2]-ezpe[1]}' tmp_en >> tmp_noe
done
file=tmp_noe
echo "ordered energies"
awk '{a[++d]=$2}
END{
n = asort(a,b)
for (i=1; i<=n; i++){
 print b[i]
 }
}' $file > tmp_oe
echo "final energies"
paste tmp_oe tmp_noe | awk '{oe[NR]=$1;flag[NR]=$2;e[NR]=$3} 
END{
i=1
while(i<=NR){
  j=1
  while(j<=NR){
    if(oe[i]==e[j]) print "TS",i,flag[j],e[j]
    ++j
    }
  ++i
  }
}' > $tsdirll/TSs/SORTED/TSlist_sorted
set `awk '{print $2}' $tsdirll/TSs/SORTED/TSlist_sorted`
for i
do
    name=`awk 'NR=='$i',NR=='$i'{print $3}' $tsdirll/TSs/SORTED/TSlist_sorted`
    namenr=$(basename $name .rxyz)
    echo $i $name
##keep this for the moment
    name0=TS$i
    names=TS${i}_${name}
#    cp $tsdirll/TSs/${name} $tsdirll/TSs/SORTED/${names}
##insert data into tss 
    sqlite3 "" "attach '${tsdirll}/TSs/ts.db' as ts; attach '${tsdirll}/TSs/SORTED/tss.db' as tss;
    insert into tss (natom,name,lname,energy,zpe,g,geom,freq) select natom,'$name0','$names',energy,zpe,g,geom,freq from ts where name='$namenr';"
done

# We now go on with the MINs
if [ -f tmp_oe ]; then rm tmp_oe ; fi
if [ -f tmp_noe ]; then rm tmp_noe ; fi

echo "non-ordered energies"
for name in $(sqlite3 ${tsdirll}/MINs/norep/minnr.db "select name from minnr")
do
  echo $name
  namer=${name}.rxyz
  if [ $toen -eq 1 ]; then
     dum_en=$(sqlite3 ${tsdirll}/MINs/norep/minnr.db "select energy,zpe from minnr where name='$name'" | sed 's@|@ @g')
  elif [ $toen -eq 0 ]; then
     dum_en=$(sqlite3 ${tsdirll}/MINs/norep/minnr.db "select energy,g from minnr where name='$name'" | sed 's@|@ @g')
  fi
  echo "$dum_min0"   >tmp_en
  echo "$dum_en"    >>tmp_en
  awk '{e[NR]=$1;ezpe[NR]=$2};END{printf "%10s %10.6f\n","'$namer'",e[2]-e[1]+ezpe[2]-ezpe[1]}' tmp_en >> tmp_noe
done
file=tmp_noe
echo "ordered energies"
awk '{a[++d]=$2}
END{
n = asort(a,b)
for (i=1; i<=n; i++){
 print b[i]
 }
}' $file > tmp_oe
paste tmp_oe tmp_noe > tmp_energies
echo "final energies"
awk '{oe[NR]=$1;flag[NR]=$2;e[NR]=$3}
END{
i=1
while(i<=NR){
  j=1
  while(j<=NR){
    if(oe[i]==e[j]) print "MIN",i,flag[j],e[j]
    ++j
    }
  ++i
  }
}' tmp_energies > $tsdirll/MINs/SORTED/MINlist_sorted

sed 's/_min/ min/g' $tsdirll/MINs/SORTED/MINlist_sorted >tmp_mls
set `awk '{print $2}' $tsdirll/MINs/SORTED/MINlist_sorted`
for i
do
    name=`awk 'NR=='$i',NR=='$i'{print $3}' $tsdirll/MINs/SORTED/MINlist_sorted`
    namenr=$(basename $name .rxyz)
    name2=`awk 'NR=='$i',NR=='$i'{print $4}' tmp_mls`
    echo $i $name2
##keep this for the moment
    name0=MIN$i
    names=MIN${i}_${name2}
#    cp $tsdirll/MINs/norep/${name} $tsdirll/MINs/SORTED/${names}
##insert data into mins 
    sqlite3 "" "attach '${tsdirll}/MINs/norep/minnr.db' as minnr; attach '${tsdirll}/MINs/SORTED/mins.db' as mins;
    insert into mins (natom,name,lname,energy,zpe,g,geom,freq) select natom,'$name0','$names',energy,zpe,g,geom,freq from minnr where name='$namenr';"
done

sed 's/_min/ min/g' $tsdirll/MINs/SORTED/MINlist_sorted | awk '{print $1,$2,$4,$5}' > $tsdirll/MINs/SORTED/MINlist_sorted.log

