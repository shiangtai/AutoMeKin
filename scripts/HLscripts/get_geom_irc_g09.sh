#!/bin/bash
sharedir=${TSSCDS}/share
sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > tmp_as

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
  if(NF==6) n[i]=$2
  if(NF==6) x[i]=$4
  if(NF==6) y[i]=$5
  if(NF==6) z[i]=$6
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
print ""
}' tmp_as $1

