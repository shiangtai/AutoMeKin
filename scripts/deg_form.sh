#!/bin/bash
sharedir=${AMK}/share

sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > atsdum2.out

awk '{if( NR == FNR)  {l[NR]=$1;tne=NR }}
{if(NR>FNR) {
IGNORECASE = 1
i=1
while(i<=tne){
  if( $1 == l[i]) {++n[i];d[i,n[i]]=$2}
  i++
  }
 }
}
END{
i=1
while(i<=105){
  if(n[i]>0) {
     for (j=1; j<=100; j++) a[j]=0
     for (j=1; j<=n[i]; j++) a[j]=d[i,j]
     print "Atom= ",i,n[i] 
     if(n[i]==1) print d[i,1]
     nn = asort(a,b)
     for (ii=1; ii<=nn; ii++){
     if(n[i]>1 && b[ii]>0) print b[ii]
     }
  }
  i++
  }
}' atsdum2.out deg.out
