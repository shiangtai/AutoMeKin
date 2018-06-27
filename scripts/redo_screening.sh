#!/bin/bash
if [ -f tsscds.dat ];then
   echo "tsscds.dat is in the current dir"
   inputfile=tsscds.dat
else
   echo "tsscds input file is missing. You sure you are in the right folder?"
   exit
fi

molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`
bu_ts=${tsdirll}/backup
exe=$(basename $0)

if [ ! -d "$bu_ts" ]; then
   echo "Folder $bu_ts does not exist"
   echo "Have you already run SCREENING.sh?"
   exit
fi

echo "Reverting changes of SCREENING.sh"
rm -rf ${tsdirll}/DIS* ${tsdirll}/REP*
cp ${bu_ts}/* ${tsdirll}
##
##We now redo the screening with the new set of bigerr,avgerr and thdiss
irc.sh screening
echo "Exiting $exe"

