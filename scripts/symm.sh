#!/bin/bash
sharedir=${AMK}/share
sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > atsdum2.out

file=$1
#  echo $file
awk '{if (NR == FNR) {l[NR]=$1;tne=NR}}
{if(NR > FNR && FNR==1) print $0}
{IGNORECASE = 1}
{if(NR > FNR && FNR>1) {
    i=1
    while(i<=tne){
      if( $1 == l[i]) print i,$2,$3,$4  
      i++
      }
  }
}' atsdum2.out $file >symm.dat

symm0.exe <symm.dat > tmp_symm
cont=`awk 'BEGIN{cont=1};/No more calc/{cont=0};END{print cont}' tmp_symm `

#echo $cont
if [ $cont -eq 1 ]; then  symm.exe <symm.dat>> tmp_symm ; fi
#  symm.exe <symm.dat>> tmp_symm
#else
#  echo "No more calc"
#fi
