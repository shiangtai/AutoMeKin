natom=$1
natomA=$2
inputfile=amk.dat
thd=`awk 'BEGIN{f=0};/NOcreatethdist/{f=1};END{print f}' $inputfile `
sep_frag.sh $natomA 5 > mingeom
echo "1" $natom >sprint.dat
createthdist.sh $thd
createMat.sh
cat ConnMat >> sprint.dat
sprint.exe<sprint.dat>sprint.out
awk 'BEGIN{res=0}
/Results for the Laplacian/{getline
if($8==0) res=1}
END{print res}' sprint.out
