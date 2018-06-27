#!/bin/bash
sharedir=${TSSCDS}/share
sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > atsdum2.out

awk 'BEGIN{huge=1000000}
{if( NR == FNR) l[NR]=$1}
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
END{
i=1
while(i<=natom){
  print l[n[i]],x[i],y[i],z[i]
  i++
  }
}' atsdum2.out $1
