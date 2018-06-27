file=$1
echo "Labels" >labels
awk '{if(NF==4) print $0 }' $file > mingeom_$file 
awk '{if(NF==4) print $1 }' $file >> labels 
natom=`wc -l mingeom_$file | awk '{print $1}' `
echo "1" $natom > sprint.dat
createMat3.sh mingeom_$file
cat ConnMat >> sprint.dat
cat ConnMat >> sprint.dat
sprint2.exe <sprint.dat >sprint.out

paste <(awk '{if(NR>1) print $0}' labels) <(deg.sh) >deg.out

deg_form.sh > deg_form.out
##
echo "This is a just to see if there is more than one fragment" > tmp_data
#
format.sh tmp $PWD 0.005

rm -f mingeom_$file

awk '{ndis=$1};END{print ndis}' tmp_data 
