#!/bin/bash
# script to remove ts $1
source utils.sh
#On exit remove tmp files
tmp_files=(rmts_arg*)
trap 'err_report $LINENO' ERR
trap cleanup EXIT INT

exe=$(basename $0)
cwd=$PWD
if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi

for var in "$@"
do
    echo $var >> rmts_arg
done

molecule=`awk '{if($1=="molecule") print $2}' $inputfile`
tsdirhl=`awk '{if($1 == "tsdirhl") {print $2;nend=1}};END{if(nend==0) print "'$cwd'/tsdirHL_'$molecule'"}' $inputfile`

#backup RXNet files
cp ${tsdirhl}/KMC/RXNet ${tsdirhl}/KMC/RXNet_backup
cp ${tsdirhl}/KMC/RXNet.cg ${tsdirhl}/KMC/RXNet.cg_backup
cp ${tsdirhl}/KMC/RXNet_long.cg ${tsdirhl}/KMC/RXNet_long.cg_backup
cp ${tsdirhl}/KMC/RXNet_long.cg_groupedprods ${tsdirhl}/KMC/RXNet_long.cg_groupedprods_backup
cp ${tsdirhl}/KMC/RXNet.relevant ${tsdirhl}/KMC/RXNet.relevant_backup
#remove $1 ts in rxnet files
awk '{if(NR==FNR){ts[NR]=$1;++nts}};{if(NR>FNR) {lp=0;for(i=1;i<=nts;i++) {if($2==ts[i]) lp=1}; if(lp==0) print $0}}' rmts_arg ${tsdirhl}/KMC/RXNet_backup > ${tsdirhl}/KMC/RXNet
awk '{if(NR==FNR){ts[NR]=$1;++nts}};{if(NR>FNR) {lp=0;for(i=1;i<=nts;i++) {if($2==ts[i]) lp=1}; if(lp==0) print $0}}' rmts_arg ${tsdirhl}/KMC/RXNet.cg_backup > ${tsdirhl}/KMC/RXNet.cg
awk '{if(NR==FNR){ts[NR]=$1;++nts}};{if(NR>FNR) {lp=0;for(i=1;i<=nts;i++) {if($2==ts[i]) lp=1}; if(lp==0) print $0}}' rmts_arg ${tsdirhl}/KMC/RXNet_long.cg_backup > ${tsdirhl}/KMC/RXNet_long.cg
awk '{if(NR==FNR){ts[NR]=$1;++nts}};{if(NR>FNR) {lp=0;for(i=1;i<=nts;i++) {if($2==ts[i]) lp=1}; if(lp==0) print $0}}' rmts_arg ${tsdirhl}/KMC/RXNet_long.cg_groupedprods_backup > ${tsdirhl}/KMC/RXNet_long.cg_groupedprods
awk '{if(NR==FNR){ts[NR]=$1;++nts}};{if(NR>FNR) {lp=0;for(i=1;i<=nts;i++) {if($2==ts[i]) lp=1}; if(lp==0) print $0}}' rmts_arg ${tsdirhl}/KMC/RXNet.relevant_backup > ${tsdirhl}/KMC/RXNet.relevant


#Re-do the kmc simulations and gather results in final directory
KMC.sh
FINAL.sh

#remove tss from TSinfo
awk 'BEGIN{ch=1}
{if(NR==FNR){ts[NR]=$1;++nts}}
/Conformational/{ch=0}
{if(ch==0) print $0}
{if(NR>FNR && ch==1) {lp=0;for(i=1;i<=nts;i++) {if($1==ts[i]) lp=1}; if(lp==0) print $0}}{if(NR==FNR){ts[NR]=$1;++nts}}' rmts_arg FINAL_HL_${molecule}/TSinfo > FINAL_HL_${molecule}/TSinfo_backup 

cp FINAL_HL_${molecule}/TSinfo_backup FINAL_HL_${molecule}/TSinfo

#recover the original rxnet files
cp ${tsdirhl}/KMC/RXNet_backup ${tsdirhl}/KMC/RXNet
cp ${tsdirhl}/KMC/RXNet.cg_backup ${tsdirhl}/KMC/RXNet.cg
cp ${tsdirhl}/KMC/RXNet_long.cg_backup ${tsdirhl}/KMC/RXNet_long.cg
cp ${tsdirhl}/KMC/RXNet_long.cg_groupedprods_backup ${tsdirhl}/KMC/RXNet_long.cg_groupedprods
cp ${tsdirhl}/KMC/RXNet.relevant_backup ${tsdirhl}/KMC/RXNet.relevant
