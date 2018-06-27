#!/bin/bash
#default sbatch FT2 resources
#SBATCH --time=08:00:00
#SBATCH -n 4
#SBATCH --output=IRC-%j.log
#_remove_this_in_ft_SBATCH --partition=cola-corta,thinnodes
#SBATCH --ntasks-per-node=2
#SBATCH -c 12
#
# SBATCH -p shared --qos=shared
# SBATCH --ntasks-per-node=2
# SBATCH -c 10

sharedir=${TSSCDS}/share

exe="IRC.sh"
source utils.sh
#On exit remove tmp files
tmp_files=(atsdum2.out black* ConnMat deg* labels mingeom ScalMat sprint.out tmp* screening.log)
#trap 'err_report $LINENO' ERR
trap cleanup EXIT INT
#current working dir

if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi
#Make sure g09 is submitted to slurm in ft2
if [ ! -z $SLURM_JOB_ID ] && [ ! -z $SLURM_NTASKS ]; then
  t=$(srun -N 1 -n 1 g09<${sharedir}/g09_test.dat | awk 'BEGIN{t=0};/Normal ter/{t=1};END{print t}')
else
  t=$(g09<${sharedir}/g09_test.dat | awk 'BEGIN{t=0};/Normal ter/{t=1};END{print t}')
fi

if [ $t -eq 0 ]; then
   echo "Please check that gaussian09 is installed in your computer and it can be invoked as g09"
   exit 1
fi


#####
#####
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
natom="$(awk 'NR==1,NR==1{print $1}' $molecule.xyz)"
rate=` awk '/Rate/{if ($2=="canonical") print "0";if($2=="microcanonical") print "1" }' $inputfile `
energy=` awk 'BEGIN{e=0};/EKMC/{e=$2};END{print e}'  $inputfile `
maxen=`awk 'BEGIN{if('$rate'==0) en=100;if('$rate'==1) en='$energy'};{if($1=="MaxEn") en=$2};END{print en}' $inputfile `
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
charge=`awk 'BEGIN{ch=0};{if($1=="charge") ch=$2};END{print ch}' $inputfile `
mult=`awk 'BEGIN{mu=1};{if($1=="mult") mu=$2};END{print mu}' $inputfile `

min_name=`awk '/molecule/{print $2}' $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "tsdirLL_'$molecule'"}' $inputfile`
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "tsdirHL_'$molecule'"}' $inputfile`
##
readhl
##

if [ "$HLcalc" == "" ]; then 
  echo "please provide a HL calc type"
  exit
else
  echo "HL calc: " $HLcalc
fi
###create table minhl.db
sqlite3 ${tsdirhl}/MINs/minhl.db "drop table if exists minhl; create table minhl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT, sigma INTEGER, unique(name));"

echo "Molecule name" $min_name
echo "tsdirll is " $tsdirll

sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt} 
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} '  > tmp_as

#
avgerr=`awk '/avgerr/{avg=$2;if(avg>0.001) avg=0.001};END{print avg}' $inputfile`
bigerr=`awk '/bigerr/{big=$2;if(big>1) big=1};END{print big}' $inputfile`
thdiss=`awk '/thdiss/{print $2}' $inputfile`

# First we copy min0 in MIN directory 
echo "Moving min0 to its final location"
energy=`get_energy_g09_$HLcalc.sh $tsdirhl/min0.log $noHLcalc`

###EMN
name=min0
zpe=$(get_ZPE_g09.sh $tsdirhl/min0.log)
g=$(get_G_g09.sh $tsdirhl/min0.log)
geom="$(get_geom_g09.sh $tsdirhl/min0.log)"
freq="$(get_freq_g09.sh $tsdirhl/min0.log)"
sigma=$(awk 'BEGIN{IGNORECASE=1};/SYMMETRY NUMBER/{print $NF;exit}' $tsdirhl/min0.log | sed 's@\.@@' )
sqlite3 ${tsdirhl}/MINs/minhl.db "insert into minhl (natom,name,energy,zpe,g,geom,freq,sigma) values ($natom,'$name',$energy,$zpe,$g,'$geom','$freq',$sigma);"

# Now we do things specific of IRC 
if [ ! -d "$tsdirhl/IRC" ]; then
   echo "$tsdirhl/IRC does not exist. It will be created"
   mkdir $tsdirhl/IRC
else
   echo "$tsdirhl/IRC already exists."
fi
if [ ! -d "$tsdirhl/TSs" ]; then
   echo "$tsdirhl/TSs does not exist. It will be created"
   mkdir $tsdirhl/TSs
else
   echo "$tsdirhl/TSs already exists"
fi
###Move checkpoint files to IRC folder
cp $tsdirhl/ts*.chk $tsdirhl/IRC
###emilio
if [ -f "black_list.dat" ]; then rm black_list.dat; fi
if [ -f "black_list.out" ]; then rm black_list.out; fi
echo "List of disconnected (at least two fragments) TS structures" > $tsdirhl/TSs/tslist_disconnected
echo "Screening" > $tsdirhl/TSs/tslist_screened
nts=0
nrm=0
ndi=0
file=${tsdirll}/tslist

sqlite3 $tsdirhl/MINs/minhl.db "select energy,zpe from minhl where name='min0'" | sed 's@|@ @g' >tmp_min0
en_min0=$(sqlite3 $tsdirhl/MINs/minhl.db "select energy,g from minhl where name='min0'" | sed 's@|@ @g' | awk '{printf "%20.10f\n",$1+$2}')

###tshl and tshldscnt tshlrep tshlhe table
sqlite3 ${tsdirhl}/TSs/tshl.db "drop table if exists tshl; create table tshl (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT,number INTEGER,sigma INTEGER);"
sqlite3 ${tsdirhl}/TSs/tshldscnt.db "drop table if exists tshldscnt; create table tshldscnt (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT,number INTEGER);"
sqlite3 ${tsdirhl}/TSs/tshlrep.db "drop table if exists tshlrep; create table tshlrep (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT,number INTEGER);"
sqlite3 ${tsdirhl}/TSs/tshlhe.db "drop table if exists tshlhe; create table tshlhe (id INTEGER PRIMARY KEY,natom INTEGER, name TEXT, energy REAL,zpe REAL,g REAL,geom TEXT,freq TEXT,number INTEGER);"
for name in $(awk '{print $3}' $file)
do
  echo "Checking $name"
  number=$(echo $name | sed 's@ts@@;s@_@ @' | awk '{print $1}')
  if [ -f $tsdirhl/${name}.log ]; then
    energy=`get_energy_g09_$HLcalc.sh $tsdirhl/${name}.log $noHLcalc`
    awk 'BEGIN{huge=1000000;zero=0;gcorr=0;zpe=0}
    {if( NR == FNR) l[NR]=$1}  
    /Frequencies/{++nfreq
    if($3<0 && $4>0 && nfreq==1) ok=1
    }
    /Zero-point vibrational energy/{getline;zpe=$1}
    /Thermal correction to Gibbs Free Energy/{gcorr=$7}
    /Standard orientation:/{ getline
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
    }' tmp_as $tsdirhl/$name.log > tmp_str
    ok=`awk 'NR==1,NR==1{print $1}' tmp_str`
    if [ $ok -eq 1 ]; then
       ((nts=nts+1))
###insert into tshl table
       zpe=$(get_ZPE_g09.sh $tsdirhl/${name}.log)
       g=$(get_G_g09.sh $tsdirhl/${name}.log)
       geom="$(get_geom_g09.sh $tsdirhl/${name}.log)"
       freq="$(get_freq_g09.sh $tsdirhl/${name}.log)"
       sigma=$(awk 'BEGIN{IGNORECASE=1};/SYMMETRY NUMBER/{print $NF;exit}' $tsdirhl/${name}.log | sed 's@\.@@' )
       sqlite3 ${tsdirhl}/TSs/tshl.db "insert into tshl (natom,name,energy,zpe,g,geom,freq,number,sigma) values ($natom,'$name',$energy,$zpe,$g,'$geom','$freq',$number,$sigma);"
#Now we screen the list to remove duplicates
       echo ts$number"_out data"> $tsdirhl/TSs/$name"_data"
       echo "Labels" >labels
       awk '{if(NR>=3) print $0 }' tmp_str > mingeom 
       awk '{if(NR>=3) print $1 }' tmp_str >> labels 
       createthdist.sh $thd
       createMat2.sh
       echo "1 $natom" | cat - ConnMat | sprint.exe >sprint.out
  
       paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out

       deg_form.sh > deg_form.out
##
       awk 'NR==2,NR==2{print $2,$4}' tmp_str > tmp_ene
       cat tmp_min0 tmp_ene >tmp_en
       awk '{e[NR]=$1;ezpe[NR]=$2};END{print 627.509*(e[2]-e[1])+ezpe[2]-ezpe[1]}' tmp_en  >> $tsdirhl/TSs/$name"_data"
##
       format.sh $name $tsdirhl/TSs $thdiss 
       ndis=`awk '{ndis=$1};END{print ndis}' $tsdirhl/TSs/$name"_data" `
### mv TSs where there is 2 or more fragments already formed
       if  [[ ("$ndis" -gt "1") ]]
       then
         ((ndi=ndi+1))
        echo "Structure $name removed-->fragmented TS"
        sqlite3 "" "attach '${tsdirhl}/TSs/tshl.db' as tshl; attach '${tsdirhl}/TSs/tshldscnt.db' as tshldscnt;
        insert into tshldscnt (natom,name,energy,zpe,g,geom,freq,number) select natom,name,energy,zpe,g,geom,freq,number from tshl where name='$name';delete from tshl where name='$name';"
       fi
###
       cat $tsdirhl/TSs/$name"_data" >> $tsdirhl/TSs/tslist_screened
    else
       echo "failed to optimize $name"
       continue
    fi
  else
    echo $name "has been removed because of previous repetitions or it does not exist"
    continue
  fi
  rm $tsdirhl/TSs/$name"_data"
done 
ntot=$(awk 'END{print NR}' $file)
#reduce output
reduce.sh $tsdirhl/TSs ts
awk '{if($NF==1) print $0}' $tsdirhl/TSs/tslist_screened.red >  $tsdirhl/TSs/tslist_screened.redconn
awk '{if($NF> 1) print $0}' $tsdirhl/TSs/tslist_screened.red >> $tsdirhl/TSs/tslist_disconnected
diffGT.sh $tsdirhl/TSs/tslist_screened.redconn $tsdirhl/TSs ts  $avgerr $bigerr
#remove repeated structures in this dir and also in TSs
if [ -f "black_list.out" ]; then 
   set `awk '{print $0}' black_list.out `
   for i
   do
     orig=$(awk '{if($2=='$i') {print $1;exit}}' $tsdirhl/TSs/tslist_screened.lowdiffs)
     echo "Structure ts"$i "removed-->redundant with ts$orig"
#     rm $tsdirhl/TSs/ts$i"_"*
     sqlite3 "" "attach '${tsdirhl}/TSs/tshl.db' as tshl; attach '${tsdirhl}/TSs/tshlrep.db' as tshlrep;
     insert into tshlrep (natom,name,energy,zpe,g,geom,freq,number) select natom,name,energy,zpe,g,geom,freq,number from tshl where number='$i';delete from tshl where number='$i';"
     ((nrm=nrm+1))
   done
else
   echo "No repetitions"
fi
((nfin=nts-nrm-ndi))
echo "Of the total" $ntot "TSs optimized at the LL, a total of" $nts "have been optimized at the HL"
echo "$nrm removed because of repetitions"
echo "$ndi removed because they are fragmented TSs"
echo "Current number of TSs optimized at the HL=" $nfin

if [ $nfin -eq 0 ]; then
   echo "No TS optimized at the high-level. Check the output files"
   exit
fi


##############################
##       Now run IRC        ##
##############################
echo "now the IRC dir "
sed  's@iop@'$iop'@;s/opt=(ts,noeigentest,calcall,noraman)/guess=read geom=check irc=(forward,maxpoints=100,rcfc,recalc=10,stepsize=10) iop(1\/108=-1)/;s@level1@'$level1'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/hl_input_template_notemp > tmp_hlircf_template0
sed  's@iop@'$iop'@;s/opt=(ts,noeigentest,calcall,noraman)/guess=read geom=check irc=(reverse,maxpoints=100,rcfc,recalc=10,stepsize=10) iop(1\/108=-1)/;s@level1@'$level1'@;s/charge/'$charge'/;s/mult/'$mult'/' $sharedir/hl_input_template_notemp > tmp_hlircr_template0
m=0
sqlite3 ${tsdirhl}/IRC/inputs.db "drop table if exists gaussian; create table gaussian (id INTEGER PRIMARY KEY,name TEXT, input TEXT, unique (name));"
for i in $(sqlite3 $tsdirhl/TSs/tshl.db "select name from tshl")
do
  en_ts=$(sqlite3 $tsdirhl/TSs/tshl.db "select energy,g from tshl where name='$i'" | sed 's@|@ @g' | awk '{printf "%20.10f\n",$1+$2}')
  deltg="$(echo "$en_min0" "$en_ts" | awk '{printf "%20.10f\n",($2-$1)*627.51}')"
  res=$(echo "$deltg < $maxen" | bc )
  if [ -f $tsdirhl/IRC/ircf_$i.log ] && [ -f $tsdirhl/IRC/ircr_$i.log ]; then
    echo "IRC completed for $i"
  elif [ $res -eq 0 ]; then
    echo "The energy of TS $i is: $deltg, which is greater than the threshold: $maxen" 
    sqlite3 "" "attach '${tsdirhl}/TSs/tshl.db' as tshl; attach '${tsdirhl}/TSs/tshlhe.db' as tshlhe;
    insert into tshlhe (natom,name,energy,zpe,g,geom,freq,number) select natom,name,energy,zpe,g,geom,freq,number from tshl where name='$i';delete from tshl where name='$i';"
  else
    ((m=m+1))
    echo "Submit IRC calc for" $i
#set-up gaussin09 calculation for $i
###Rename checkpoint files to be employed in the IRC calculations
    mv $tsdirhl"/IRC/"$i".chk" $tsdirhl"/IRC/ircf_"$i".chk"
    cp $tsdirhl"/IRC/ircf_"$i".chk" $tsdirhl"/IRC/ircr_"$i".chk"
###
    chkf="$(echo "%chk=ircf_"$i)"
    chkr="$(echo "%chk=ircr_"$i)"
# if the imag <100 then stepsize=30
    imag=`awk 'BEGIN{fl=0};/Frequencies/{f0=$3;f=sqrt(f0*f0);if(f<100)fl=1;print fl ;exit}' $tsdirhl/$i.log `

    if [ $imag -eq 1 ] ; then
       echo "ts $i has an imaginary freq lower than 100 cm-1"      
       sed 's/stepsize=10/stepsize=30/g' tmp_hlircf_template0 > tmp_hlircf_template
       sed 's/stepsize=10/stepsize=30/g' tmp_hlircr_template0 > tmp_hlircr_template
    else
       cp tmp_hlircf_template0 tmp_hlircf_template
       cp tmp_hlircr_template0 tmp_hlircr_template
    fi

    calf="$(cat tmp_hlircf_template)"
    calr="$(cat tmp_hlircr_template)"
    ig09f="$(echo -e "$chkf"'\n'"$calf")"
    ig09r="$(echo -e "$chkr"'\n'"$calr")"
#    names[$m]=ircf_"$i"
    echo -e "insert or ignore into gaussian values (NULL,'ircf_$i','$ig09f');\n.quit" | sqlite3 ${tsdirhl}/IRC/inputs.db
    ((m=m+1))
#    names[$m]=ircr_"$i"
    echo -e "insert or ignore into gaussian values (NULL,'ircr_$i','$ig09r');\n.quit" | sqlite3 ${tsdirhl}/IRC/inputs.db
  fi 
done
####
#Perform m parallel calculations
echo Performing a total of $m irc calculations
if [ $m -gt 0 ]; then
#   doparallel "runIRC.sh $tsdirhl" "$(echo ${names[@]})"
   doparallel "runIRC.sh {1} $tsdirhl" "$(seq $m)"
fi


