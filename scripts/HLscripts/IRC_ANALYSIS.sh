#!/bin/bash
source utils.sh
sharedir=${TSSCDS}/share

if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
#####
echo "Input file" $inputfile
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
natom="$(awk 'NR==1,NR==1{print $1}' $molecule.xyz)"
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "tsdirHL_'$molecule'"}' $inputfile`

#On exit remove tmp files
tmp_files=($tsdirhl/IRC/*.chk black* minfailed_list atsdum2.out labels mingeom sprint.* deg* tmp* ConnMat ScalMat) 
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
if [ ! -d "$tsdirhl/PRODs" ]; then
   echo "$tsdirhl/PRODs does not exist. It will be created"
   mkdir $tsdirhl/PRODs
else
   echo "$tsdirhl/PRODs already exists. Remove PR files"
   rm -f $tsdirhl/PRODs/PR*
fi
if [ ! -d "$tsdirhl/MINs/norep" ]; then
   echo "$tsdirhl/MINs/norep does not exist. It will be created"
   mkdir $tsdirhl/MINs/norep
else
   echo "$tsdirhl/MINs/norep already exists."
   rm -r $tsdirhl/MINs/norep
   mkdir $tsdirhl/MINs/norep
fi
#emilio
##Reading HL stuff
readhl
##

thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `

sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > atsdum2.out

#
avgerr=`awk '/avgerr/{avg=$2;if(avg>0.001) avg=0.001};END{print avg}' $inputfile`
bigerr=`awk '/bigerr/{big=$2;if(big>1) big=1};END{print big}' $inputfile`



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
#remove some stuff at the beginning (in case of repeating the same script)
rm -f $tsdirhl/MINs/*min*ts*

#
echo "Nomenclature" > $tsdirhl/MINs/names_of_minima
echo "Screening" > $tsdirhl/MINs/minlist_screened
dum_min0="$(sqlite3 $tsdirhl/MINs/minhl.db "select energy,zpe from minhl where name='min0'" | sed 's@|@ @g' )"

#dum_min0="$(awk '{if(NR==2) print $2,$4}' $tsdirhl/MINs/min0.rxyz)"
### and prodhl minnrhl table
sqlite3 ${tsdirhl}/PRODs/prodhl.db "drop table if exists prodhl; create table prodhl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, formula TEXT);"
sqlite3 ${tsdirhl}/MINs/norep/minnrhl.db "drop table if exists minnrhl; create table minnrhl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT,energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, sigma INTEGER);"

echo "List of bad MINIMA" > minfailed_list
echo "PR list" > $tsdirhl/PRODs/PRlist
n=0
nmin=0
npro=0
nrm=0
cp $tsdirhl/min0.log $tsdirhl/IRC
for i in $(ls $tsdirhl/IRC/min*.log)
do 
  ((n=n+1))
  name=$(basename $i .log)
# See if it did not crashed and grab geometires 
  energy=`get_energy_g09_$HLcalc.sh $i $noHLcalc`
  awk 'BEGIN{huge=1000000;zero=0;gcorr=0;zpe=0}
  {if( NR == FNR) l[NR]=$1}
  /Frequencies/{++nfreq
  if($3>0 && $4>0 && nfreq==1) ok=1
  }
  /SCF Done/{e=$5}
  /Zero-point vibrational energy/{getline;zpe=$1}
  /Thermal correction to Gibbs Free Energy/{gcorr=$7}
  /orientation:/{ getline
  getline
  getline
  getline
  i=1
  while(i<=huge){
    getline
    if(NF==1) break
    n[i]=$2
    x[i]=$4
    y[i]=$5
    z[i]=$6
    natom=i
    i++
    }
  }
  /Error termi/{ok=0}
  END{
  if(ok==0) ok=zero
  print ok,natom
  printf "%2s %14s %4s %8.2f %5s %14.9f\n","E=","'$energy'","ZPE=",zpe,"Gcorr",gcorr
  i=1
  while(i<=natom){
    print l[n[i]],x[i],y[i],z[i]
    i++
    }
  }' atsdum2.out $i > tmpstr

  ok=`awk 'NR==1,NR==1{print $1}' tmpstr`
  if [ $ok -eq 0 ]; then
     echo "$name check later on-->this is a product or problem in opt"
  fi
  if [ $ok -eq 1 ]; then
     echo "$name optimized correctly"
#     awk '{if(NR==1) 
#             print $2 
#           else
#             print $0}' tmpstr > $tsdirhl/MINs/$name.rxyz
#     get_freq_g09.sh $i >> $tsdirhl/MINs/$name.rxyz
# insert all minima except min0
     if [ $name != "min0" ]; then
        zpe=$(get_ZPE_g09.sh $i)
        g=$(get_G_g09.sh $i)
        geom="$(get_geom_g09.sh $i)"
        freq="$(get_freq_g09.sh $i)"
###EMN
        sqlite3 ${tsdirhl}/MINs/minhl.db "insert or ignore into minhl (natom,name,energy,zpe,g,geom,freq) values ($natom,'$name',$energy,$zpe,$g,'$geom','$freq');"
     fi

#Now we screen the list to rename duplicates
     echo "Labels" >labels
     awk '{if(NR>=3) print $0 }' tmpstr > mingeom 
     awk '{if(NR>=3) print $1 }' tmpstr >> labels 
     if [ $vdw -eq 0 ]; then
        createthdist.sh $thd
     fi
#     natom=`wc -l mingeom| awk '{print $1}' `
     echo "1" $natom > sprint.dat
     createMat2.sh
     cat ConnMat >> sprint.dat
     createMat.sh
     cat ConnMat >> sprint.dat
     sprint2.exe <sprint.dat >sprint.out

     paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out
     deg_form.sh > deg_form.out
###
     dum_en="$(awk 'NR==2,NR==2{print $2,$4}' tmpstr)"
     echo "$dum_min0"   >tmp_en
     echo "$dum_en"    >>tmp_en

##     awk '{e[NR]=$1;ezpe[NR]=$2};END{print 627.509*(e[2]-e[1])+ezpe[2]-ezpe[1]}' tmp_en >> $tsdirhl/MINs/$name"_data"
##use absolute energy instead of relative one
     echo $energy >> $tsdirhl/MINs/${name}_data
##
     format.sh $name $tsdirhl/MINs 0.005
     ndis=`awk '{ndis=$1};END{print ndis}' $tsdirhl/MINs/${name}_data `
### mv MINs where there is 2 or more fragments already formed
     if  [[ ("$ndis" -gt "1") ]]
     then 
        ((npro=npro+1)) 
        echo "Products=" $ndis $name 
        namepr=PR${npro}_${name}
###remove this later on
#        mv $tsdirhl/MINs/$name.rxyz $tsdirhl/PRODs/PR$npro"_"$name.rxyz
##insert data into prodhl.db and delete it from min.db
        sqlite3 "" "attach '${tsdirhl}/MINs/minhl.db' as minhl; attach '${tsdirhl}/PRODs/prodhl.db' as prodhl; insert into prodhl (natom,name,energy,zpe,g,geom,freq) select natom,'$namepr',energy,zpe,g,geom,freq from minhl where name='$name';delete from minhl where name='$name'"
        echo "PROD" $npro $name.rxyz >> $tsdirhl/PRODs/PRlist
     else
        ((nmin=nmin+1))
        echo "min" $nmin "-->" $name.rxyz >> $tsdirhl/MINs/names_of_minima
        echo "min"$nmin "data" >> $tsdirhl/MINs/minlist_screened
        cat $tsdirhl/MINs/${name}_data >> $tsdirhl/MINs/minlist_screened
     fi
###
  else
     echo $name.rxyz >> minfailed_list
     continue
  fi
  rm $tsdirhl/MINs/${name}_data

done 
echo "Total number of minima" $nmin
### if vdw=1 restore thdist
if [ $vdw -eq 1 ]; then
   cp thdist_backup thdist
fi

#reduce output
reduce.sh $tsdirhl/MINs min
awk '{if($NF==1) print $0}' $tsdirhl/MINs/minlist_screened.red >  $tsdirhl/MINs/minlist_screened.redconn
awk '{if($NF> 1) print $0}' $tsdirhl/MINs/minlist_screened.red >> $tsdirhl/MINs/minlist_disconnected
diffGT.sh $tsdirhl/MINs/minlist_screened.redconn $tsdirhl/MINs min $avgerr $bigerr
#remove repeated structures in this dir and also in MINs
cp $tsdirhl/MINs/names_of_minima $tsdirhl/MINs/names_of_minima_norep
if [ -f "black_list.out" ]; then 
   set `awk '{print $0}' black_list.out `
   for i
   do
     echo "Structure min"$i "repeated"
     ((nrm=nrm+1))
     nomnr="$(awk '{if($2 != '$i') print $0}' $tsdirhl/MINs/names_of_minima_norep)"
     echo "$nomnr" > $tsdirhl/MINs/names_of_minima_norep 
   done
else
   echo "No repetitions"
fi
###
nn=0
for name in $(awk '{if(NR>1) print $4}' ${tsdirhl}/MINs/names_of_minima_norep)
do
  ((nn=nn+1))
  namenrxyz=$(basename $name .rxyz)
  number=`awk 'NR=='$nn'+1,NR=='$nn'+1{print $2}'  ${tsdirhl}/MINs/names_of_minima_norep`
  namenr=$(basename min${number}_${name} .rxyz)
##insert data into minnrhl.db from minhl.db
#  cp ${tsdirhl}/MINs/${name} ${tsdirhl}/MINs/norep/min${number}_${name}
##
  sqlite3 "" "attach '${tsdirhl}/MINs/minhl.db' as minhl; attach '${tsdirhl}/MINs/norep/minnrhl.db' as minnrhl; insert into minnrhl (natom,name,energy,zpe,g,geom,freq,sigma) select natom,'$namenr',energy,zpe,g,geom,freq,sigma from minhl where name='$namenrxyz';"
done


###
((nfin=nmin-nrm))
echo $nmin "have been optimized at the HL, of which" $nrm "removed because of repetitions"
echo "Current number of MINs optimized at the HL=" $nfin
#


##################
##################

echo "Now running MINFAILED.sh"

##################
##################

cp $tsdirhl/PRODs/PRlist $tsdirhl/PRODs/PRlist.old
npro=`awk '{npro=$2};END{print npro}' $tsdirhl/PRODs/PRlist `
#file to screen the failed minima and/or products

nfail=`wc -l minfailed_list | awk '{print $1-1}'`
if [ $nfail -eq 0 ]; then
   echo "You are lucky. All minima have been optimized correctly. Exit here"
   exit
else
   echo "number of minima that failed and/or are products" $nfail
fi

for i in $(awk 'NR>1{print $0}' minfailed_list)
do 
  name=$(basename $i .rxyz)
  awk 'BEGIN{huge=1000000;zero=0}
  {if( NR == FNR) l[NR]=$1}
  /Charge =/{ic=9}
  {if(ic==9 && NF==4) ++natom}
  {if(ic==9 && NF==0) ic=0}
  /orientation:/{ getline
  getline
  getline
  getline
  i=1
  print natom
  print "E=",zero,zero
  while(i<=huge){
    getline
    if(NF==1) exit 
    print l[$2],$4,$5,$6
    i++
    }
  }' atsdum2.out $tsdirhl/IRC/$name.log > tmpstr

#Now we screen the list to rename duplicates
  echo "Labels" >labels
  geom=$(awk '{if(NR>=3) print $0 }' tmpstr) 
  awk '{if(NR>=3) print $0 }' tmpstr > mingeom 
  awk '{if(NR>=3) print $1 }' tmpstr >> labels 
  createthdist.sh $thd
#  natom=`wc -l mingeom| awk '{print $1}' `
  echo "1" $natom > sprint.dat
  createMat2.sh
  cat ConnMat >> sprint.dat
  createMat.sh
  cat ConnMat >> sprint.dat
  sprint2.exe <sprint.dat >sprint.out
  paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out
  deg_form.sh > deg_form.out
##
  echo "This is a just to see if there is more than one fragment" > $tsdirhl/MINs/${name}_data
#
  format.sh $name $tsdirhl/MINs 0.005
  ndis=`awk '{ndis=$1};END{print ndis}' $tsdirhl/MINs/${name}_data `
#  echo $name $ndis
### mv MINs where there is 2 or more fragments already formed
  if  [[ ("$ndis" -gt "1") ]]
  then 
     ((npro=npro+1)) 
##remove this later on
#     cp tmpstr $tsdirhl/PRODs/PR$npro"_"$name.rxyz
     namepr="PR"$npro"_"$name
##EMNinsert into prodhl.db
     sqlite3 ${tsdirhl}/PRODs/prodhl.db "insert into prodhl (natom,name,energy,zpe,g,geom,freq) values ($natom,'$namepr',0,0,0,'$geom',0);"
     lp=`awk 'BEGIN{lp=1};/'$name'/{lp=0};END{print lp}' $tsdirhl/PRODs/PRlist `
#     echo $lp 
     if  [[ ("$lp" -eq "1") ]]
     then
         echo "PROD" $npro $name.rxyz >> $tsdirhl/PRODs/PRlist
         echo "Move structure $name to PROD # $npro"
     else
         echo "Structure $name already moved to PROD" 
     fi
  else
     echo "Double check this opt: $name"
  fi
done 
