file=$1
nc=$2
awk '/sp ca/{++fl};/MP2/{if(fl=='$nc'-1) print $6;exit}' $file >emp2
sed 's/D/ /g' emp2 > emp2s
awk '{printf "%14.9f\n",$1*10^$2} ' emp2s
