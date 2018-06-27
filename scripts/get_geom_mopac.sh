file=$1
awk 'BEGIN{lok=0}
/FINAL HE/{lok=1}
/CARTESIAN COORDINATES/{
  if(lok==0) {
  getline
  getline
  getline
  i=1
  while(i<=100000){
    getline
    if(NF==0) break
    i++
    }
  print i-1,"\n"
  }
  if(lok==1) {
   getline
   i=1
   while(i<=100000){
    getline
    if(NF==0) exit
    print $2,$3,$4,$5
    i++
    }
  }
}
END{
if(lok==0) print "Error"
}' $file
