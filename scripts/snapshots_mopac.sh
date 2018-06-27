file=$1
awk '{if(NR==1) {natom=$1;print natom} }
{init=(NR-3)/(natom+2)
fin=NR/(natom+2)
if(init == int(init) ) p=1
if(fin  == int(fin)  ) p=0 
if(p==1 || fin  == int(fin)) {++ntot;x[ntot]=$2;y[ntot]=$3;z[ntot]=$4} }
END{
print ntot/natom
i=1
while(i<=ntot){
  print x[i],y[i],z[i] 
  ++i
  }
}' $file >> bbfs.dat




