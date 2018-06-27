#!/bin/bash
source utils.sh
sharedir=${TSSCDS}/share
#remove tmp files
tmp_files=(mingeom labels deg* tmp* black_list* ConnMat sprint.* ScalMat atsdum2.out)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT
exe=$(basename $0)
inputfile=tsscds.dat
molecule=`awk '{if($1=="molecule") print $2}' $inputfile`
natom="$(awk 'NR>2{if(NF==4) ++natom};END{print natom}' ${molecule}.xyz)"
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
avgerr=`awk '/avgerr/{print $2}' $inputfile`
bigerr=`awk '/bigerr/{print $2}' $inputfile`
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`
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

if [ ! -d "$tsdirll/PRODs" ]; then
   echo "PRODs does not exist. It will be created"
   mkdir $tsdirll/PRODs
else
   echo "PRODs already exists."
   rm -r $tsdirll/PRODs
   mkdir $tsdirll/PRODs
fi
if [ ! -d "$tsdirll/MINs/norep" ]; then
   echo "MINs/norep does not exist. It will be created"
   mkdir $tsdirll/MINs/norep
else
   echo "MINs/norep already exists."
   rm -r $tsdirll/MINs/norep
   mkdir $tsdirll/MINs/norep
fi
##create prod and norep databases
sqlite3 ${tsdirll}/PRODs/prod.db "drop table if exists prod; create table prod (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, formula TEXT);"
sqlite3 ${tsdirll}/MINs/norep/minnr.db "drop table if exists minnr; create table minnr (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, sigma INTEGER);"
#
echo "Nomenclature" > $tsdirll/MINs/names_of_minima
echo "Screening" > $tsdirll/MINs/minlist_screened
echo "List of bad MINIMA" > $tsdirll/MINs/minfailed_list
echo "PR list" > $tsdirll/PRODs/PRlist
nmin=0
#n=0
npro=0
nrm=0
for name in $(sqlite3 ${tsdirll}/MINs/min.db "select name from min")
#for i in $(ls $tsdirll/MINs/*.rxyz)
do 
#   name="$(basename $i .rxyz)"
   echo $name
   echo "Labels" >labels
#   awk '{if(NR>=3 && NF==4 ) print $0 }' $i > mingeom 
#   awk '{if(NR>=3 && NF==4 ) print $1 }' $i >> labels 
   sqlite3 ${tsdirll}/MINs/min.db "select geom from min where name='$name'"  >mingeom
   awk '{print $1}' mingeom  >>labels
   if [ $vdw -eq 0 ]; then
      createthdist.sh $thd
   fi
#   natom=`wc -l mingeom| awk '{print $1}' `
   echo "1" $natom > sprint.dat
   createMat2.sh
   cat ConnMat >> sprint.dat
   createMat.sh
   cat ConnMat >> sprint.dat
   sprint2.exe <sprint.dat >sprint.out

   paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out


   deg_form.sh > deg_form.out
###
   sqlite3 ${tsdirll}/MINs/min.db "select energy from min where name='$name'"  >> $tsdirll/MINs/${name}_data
#   awk 'NR==2,NR==2{printf "%9.3f\n",$2}' $i >> $tsdirll/MINs/$name"_data" 
##
   format.sh $name $tsdirll/MINs 0.005
   ndis=`awk '{ndis=$1};END{print ndis}' $tsdirll/MINs/${name}_data `
### mv MINs where there is 2 or more fragments already formed
   if  [[ ("$ndis" -gt "1") ]]
   then 
      ((npro=npro+1)) 
      echo "Products=" $ndis $name 
##conserve rxyz for the moment
      namepr=PR${npro}_${name}
#      mv ${tsdirll}/MINs/${name}.rxyz $tsdirll/PRODs/${namepr}.rxyz
##insert data into prod.db and delete it from min.db
      sqlite3 "" "attach '${tsdirll}/MINs/min.db' as min; attach '${tsdirll}/PRODs/prod.db' as prod;
      insert into prod (natom,name,energy,zpe,g,geom,freq) select natom,'$namepr',energy,zpe,g,geom,freq from min where name='$name';delete from min where name='$name'"
      echo "PROD" $npro $name.rxyz >> $tsdirll/PRODs/PRlist
   else
      freq=$(sqlite3 $tsdirll/MINs/min.db "select freq from min where name='$name'")
      if [ ! -z "$freq" ]; then
         ((nmin=nmin+1))
         echo "min $nmin --> ${name}.rxyz" >> $tsdirll/MINs/names_of_minima
         echo "min$nmin data" >> $tsdirll/MINs/minlist_screened
         cat $tsdirll/MINs/${name}_data >> $tsdirll/MINs/minlist_screened
      else
         echo $name.rxyz >> $tsdirll/MINs/minfailed_list
      fi
   fi
###
  rm $tsdirll/MINs/${name}_data
done 
echo "Total number of minima" $nmin
### if vdw=1 restore thdist
if [ $vdw -eq 1 ]; then
   cp thdist_backup thdist
fi
#reduce output
reduce.sh $tsdirll/MINs min
awk '{if($NF==1) print $0}' $tsdirll/MINs/minlist_screened.red >  $tsdirll/MINs/minlist_screened.redconn
awk '{if($NF> 1) print $0}' $tsdirll/MINs/minlist_screened.red >> $tsdirll/MINs/minlist_disconnected
diffGT.sh $tsdirll/MINs/minlist_screened.redconn $tsdirll/MINs min $avgerr $bigerr
#remove repeated structures in this dir and also in MINs
cp $tsdirll/MINs/names_of_minima $tsdirll/MINs/names_of_minima_norep
if [ -f "black_list.out" ]; then 
   set `awk '{print $0}' black_list.out `
   for i
   do
     echo "Structure min"$i "repeated"
     ((nrm=nrm+1))
     awk '{if($2 != '$i') print $0}' $tsdirll/MINs/names_of_minima_norep > tmp 
     cp tmp $tsdirll/MINs/names_of_minima_norep 
   done
else
   echo "No repetitions"
fi
###
nn=0
for name in $(awk '{if(NR>1) print $4}' ${tsdirll}/MINs/names_of_minima_norep)
do
  ((nn=nn+1))
  namenrxyz=$(basename $name .rxyz)
  number=`awk 'NR=='$nn'+1,NR=='$nn'+1{print $2}'  ${tsdirll}/MINs/names_of_minima_norep`
  namenr=$(basename min${number}_${name} .rxyz)
##for the moment we copy rxyz stuff
#  cp ${tsdirll}/MINs/${name} ${tsdirll}/MINs/norep/min${number}_${name}
##insert data into minnr.db and delete it from min.db
  sqlite3 "" "attach '${tsdirll}/MINs/min.db' as min; attach '${tsdirll}/MINs/norep/minnr.db' as minnr;
  insert into minnr (natom,name,energy,zpe,g,geom,freq,sigma) select natom,'$namenr',energy,zpe,g,geom,freq,sigma from min where name='$namenrxyz';"
done
###
((nfin=nmin-nrm))
echo $nmin "have been optimized, of which" $nrm "removed because of repetitions"
echo "Current number of MINs optimized =" $nfin
#
