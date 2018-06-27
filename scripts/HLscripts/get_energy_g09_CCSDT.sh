file=$1
nc=$2
awk '/sp ca/{++fl};/CCSD\(T/{if(fl=='$nc'-1) print $2;exit}' $file >eccsdt
sed 's/D/ /g' eccsdt > eccsdts
awk '{printf "%14.9f\n",$1*10^$2} ' eccsdts
