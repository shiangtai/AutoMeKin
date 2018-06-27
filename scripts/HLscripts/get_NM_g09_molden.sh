#!/bin/bash
sharedir=${TSSCDS}/share
sed "s/'/ /g;s/,/ /g;s/+00/+00 /g" ${sharedir}/atsymb | awk '/ams=/{amlab=1}
/character/{++nt}
{for(i=1;i<=(NF-1);i++) {if($1 != "real" && amlab==1 && nt==0) {++j;m[j]=$i} }}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) {++k;print $i,m[k]}
} ' > tmp_as


echo "[Molden Format]" > $2.freq
echo "[FREQ]"         >> $2.freq
##echo the freqs here
get_freq_g09.sh $1 | awk '{printf "%6.1f\n",$1}' >> $2.freq
echo "       " >> $2.freq
echo "[FR-COORD]" >> $2.freq


awk 'BEGIN{huge=10000;ifreq=0;atobohr=1.889726}
{if( NR == FNR) {l[NR]=$1;m[NR]=$2}}
/Coordinates/{getline
getline
i=1
natom=0
while(i<=huge){
  getline
  if(NF==6)  ++natom
  if(NF==6)  nl[natom]=$2
  if(NF==6)  x[natom]=$4
  if(NF==6)  y[natom]=$5
  if(NF==6)  z[natom]=$6
  if(NF==1)  break
  ++i
  }
}
/Frequencies/{++ifreq;
   getline;getline;getline;getline
   j=1
   while(j<=natom) {
     getline
     iatom=$1
     for(inf=3;inf<=NF;++inf) {i=int(inf/3);k=inf-3*i-2;nfreq=3*ifreq-3+i;nm[3*iatom+k,nfreq]=$inf}
     j++
     }
}
END{
n=1
while(n<=natom){
 print l[nl[n]],atobohr*x[n],atobohr*y[n],atobohr*z[n]  >>  "'$2'.freq"
 ++n
 }
print "" >> "'$2'.freq"
print "[FR-NORM-COORD]" >> "'$2'.freq"
i=1
while(i<=nfreq){
 print "Vibration "i >> "'$2.freq'"
   j=1
   while (j<=natom) {
   print nm[3*j-2,i]*atobohr,nm[3*j-1,i]*atobohr,nm[3*j,i]*atobohr >> "'$2'.freq"
   ++j
   }
 ++i
 }
}' tmp_as $1 




