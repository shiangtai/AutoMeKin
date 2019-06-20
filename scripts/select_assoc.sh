#!/bin/bash
sharedir=${AMK}/share
sed "s/'/ /g;s/,/ /g" $sharedir/atsymb | awk '/character/{++nt}
/end/{exit}
{i0=1;if($1 ~ /character/) i0=3
for(i=i0;i<=(NF-1);i++) if(nt>=1) print $i
} ' > atsdum2.out

inputfile=$1
frA=`awk '/A=/{print $2}' $inputfile`
frB=`awk '/B=/{print $2}' $inputfile`
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
assocdir=$PWD"/assoc_"$frA"_"$frB
rm -f $assocdir/selected_min*
natomA=`awk 'NR==1,NR==1{print $1}' $frA.xyz`
natomB=`awk 'NR==1,NR==1{print $1}' $frB.xyz`
awk '{if(NF==4) print $0}' $frA.xyz >mingeom
createthdist.sh $thd
createMat.sh
cp ConnMat ConnMatA
awk '{if(NF==4) print $0}' $frB.xyz >mingeom
createthdist.sh $thd
createMat.sh
cp ConnMat ConnMatB
n=0
nmin="$(ls $assocdir/assoc*.out | wc -l | awk '{print $1}')"
echo "A total of $nmin minima have been optimized"
for i in $(ls $assocdir/assoc*.out )
do
   ((n=n+1))
   echo "Doing analysis of min $n: $i"
   get_geom_mopac.sh $i | awk '{if(NF==4) print $0}' >mingeom0
   if [ $n -eq 1 ]; then
      met_label=`awk '{if( NR == FNR) {l[NR]=$1;tne=NR}}
      {if(FNR > 1 && NR >FNR ) {
         IGNORECASE = 1
         i=1
         while(i<=tne){
            if( $1 == l[i] && i>=21 && i<=30) print FNR
            if( $1 == l[i] && i>=39 && i<=48) print FNR
            if( $1 == l[i] && i>=72 && i<=80) print FNR
            i++
            }
        }
      }' atsdum2.out mingeom0`
      if [ -z "$met_label" ]; then met_label=0 ; fi
   fi
   awk '{
   if(NR<='$natomA') 
     print $0 > "mingeomA"
   else
     print $0 > "mingeomB"
   }' mingeom0
   cp mingeomA mingeom
   createthdist.sh $thd
   createMat.sh
   cp ConnMat ConnMatAp
   diffA="$(paste ConnMatA ConnMatAp | awk 'BEGIN{diff=0;natomA='$natomA'}
   {for(j=1;j<=natomA;j++){d=$j-$(j+natomA);diff+=sqrt(d*d)}
   }
   END{print diff}')"
   cp mingeomB mingeom
   createthdist.sh $thd
   createMat.sh
   cp ConnMat ConnMatBp
   diffB="$(paste ConnMatB ConnMatBp | awk 'BEGIN{diff=0;natomB='$natomB'}
   {for(j=1;j<=natomB;j++){diff+=$j-$(j+natomB)}
   }
   END{print diff}')"
   diff="$(echo $diffA $diffB | awk '{print $1+$2}')"
   val=`awk 'BEGIN{huge=10^10;val=0;
      if('$met_label'==0) {print "0";exit}
      while(i<=huge){
        low=6*i
        high=6*(i+1)
        if('$met_label'>low && '$met_label'<high) {pos='$met_label'-low;break}
        i++
        }
      }
      /BOND ORDERS AND VALENCIES/{getline
      getline
      getline
      getline
      i=1
      while(i<=huge){
        getline
        if($2=='$met_label' && NF==pos+2) {print $NF;exit}
        i++
        }
      }' $i`
   e=`awk '/FINAL HEAT OF FORMATION/{print $6;exit}' $i`
   echo $i $val $e>> $assocdir/selected_min_$diff
done

for i in $(seq 0 10)
do
   if [ ! -f $assocdir/selected_min_$i ]; then
      echo "No minima found with changes in the geometries of $i atoms"
      continue 
   else
      echo "Minima found with changes in the geometries of $i atoms"
      awk 'BEGIN{min=10^10}
      {name[NR]=$1
      val[NR]=$2
      e[NR]=$3
      if($3<min) min=$3
      }
      END{i=1
      while(i<=NR){
         en=e[i]-min
         point=2^val[i]-0.1*en
         print name[i],point
         i++
         }
      }' $assocdir/selected_min_$i > $assocdir/selected_min.out
      break
   fi
done

smin=`awk 'BEGIN{max=-10^10}
{if($2>max) {max=$2;name=$1 }}
END{print name }' $assocdir/selected_min.out`
echo "Selected minimum=" $smin
echo "That minimum will be now called= "$frA"_"$frB".xyz" 
get_geom_mopac.sh $smin > $frA"_"$frB".xyz"
