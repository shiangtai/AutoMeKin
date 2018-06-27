f=$1
e=$2
f1=$3
f2=$4
f3=$5
f4=$6
file=$7
awk 'BEGIN{ok=-1;el='$e'}
{
nfeq=0
df=$4-'$f'   
de=$5-el
df2=$7-'$f2'   
df3=$8-'$f3'   
df4=$9-'$f4'   
df=sqrt(df*df)
de=sqrt(de*de)
df2=sqrt(df2*df2)
df3=sqrt(df3*df3)
df4=sqrt(df4*df4)
if(df< 20) ++nfeq
if(df2<20) ++nfeq
if(df3<20) ++nfeq
if(df4<20) ++nfeq
if(de<=0.2 &&  nfeq==4 ) {ok=$2}}
END{
print ok
}' $file
