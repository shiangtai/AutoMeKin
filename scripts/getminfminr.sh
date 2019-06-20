#!/bin/bash
sharedir=${AMK}/share
name=$1
inputfile=amk.dat
method=` awk 'BEGIN{llcalc="PM7"};{if($1=="LowLevel") {$1="";llcalc=$0}};END{print llcalc}' $inputfile `
charge=` awk '/charge/{print $2}'  $inputfile `
molecule=` awk '{if($1=="molecule") print $2}'  $inputfile `
TKMC=`awk '{if($1 == "TKMC") {print $2;nend=1}};END{if(nend==0) print "298"}' $inputfile`
tsdirll=`awk '{if($1 == "tsdirll") {print $2;nend=1}};END{if(nend==0) print "'$PWD'/tsdirLL_'$molecule'"}' $inputfile`

geomf="$(get_geom_mopac.sh $tsdirll/IRC/$name"_ircf.out" | awk '{if(NF==4) print $0}')" 
geomr="$(get_geom_mopac.sh $tsdirll/IRC/$name"_ircr.out" | awk '{if(NF==4) print $0}')"
sed 's/thermo/thermo('$TKMC','$TKMC')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template >  $tsdirll/IRC/minf_$name.mop
sed 's/thermo/thermo('$TKMC','$TKMC')/;s/method/'"$method"' charge='$charge'/' $sharedir/thermo_template >  $tsdirll/IRC/minr_$name.mop
echo "$geomf" >> $tsdirll/IRC/minf_$name.mop
echo "$geomr" >> $tsdirll/IRC/minr_$name.mop

