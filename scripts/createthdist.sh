#!/bin/bash
if [ $1 -eq 1 ]; then exit; fi
sharedir=${TSSCDS}/share
awk '
NR==FNR{++i;label[i]=$1;x[i]=$2;y[i]=$3;z[i]=$4}
NR>FNR{++j;r[$1]=$2}
END{
IGNORECASE = 1
i0=1
while(i0<=i-1){
 ii=i0+1
 while(ii<=i){
  print label[i0],label[ii],1.2*( r[label[i0]] + r[label[ii]] )
  ii++
  }  
 i0++
 }
}' mingeom $sharedir/list_of_atomic_radii | awk '{i=1
a1[NR]=$1
a2[NR]=$2
a3[NR]=$3
ok=1
while(i<NR){
 if( a1[NR]==a1[i] && a2[NR]==a2[i] ) {ok=0}
 if( a2[NR]==a1[i] && a1[NR]==a2[i] ) {ok=0}
 i++
 }
if(NR==1) 
 print a1[NR],a2[NR],a3[NR]
else
 if(ok==1) print a1[NR],a2[NR],a3[NR]
}' > thdist
