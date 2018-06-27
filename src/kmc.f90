program kmc
implicit none
! KMC computes the population vs time of a number of species involved in several dynamical processes
! rate: rate constant of a given processes  
! p:    population of a given species (p0 is its initial value)
! re (pr): reactant (product) for a given process
! ndisp number of species that dissapear
! ndd species that dissapear
character*80 :: title
integer, parameter :: dp = selected_real_kind(15,307)
integer :: m,nran,nr,nesp,inran,i,j,k,l,mu,ndisp,pd,kk,ijk,tcont,tint
real(dp) :: t,a0,r2a0,suma,ptot
real(dp) :: rnd(2)
integer,dimension(:),allocatable :: p,p0,re,pr,n,ndd,cont
real (dp) ,dimension(:),allocatable :: a,rate
read(*,"(a80)") title
print "(t3,a80)",title
read(*,*) m,nesp,nran
allocate(re(m),pr(m),a(m),rate(m),p0(nesp),p(nesp),n(nesp),cont(m))
n=(/ (l,l=1,nesp) /)
do i=1,m
read(*,*) rate(I),re(i),pr(i)
enddo
ptot=0.d0
print*, "at the beginning"
do ijk=1,m
   cont(ijk)=0
enddo
do i=1,nesp
read(*,*) p0(i)
print*, "P0(",i,")=",p0(i)
ptot=ptot+p0(i)
enddo
read(*,*) ndisp
allocate(ndd(ndisp))
if(ndisp>0) read(*,*) ndd
read(*,*) tint
print "(t3,a,i4,/,t3,a,i4,/,t3,a,i4,/)","# of calcs:",nran,"# of procs:",m,"# of specs:",nesp
print "(t3,a,1p,20(e9.2))","Rates:",rate
print "(t3,a,20(i9))","Reacts:",re
print "(t3,a,20(i9))","Prods :",pr
print "(/,t3,a,1p,i10)","Step size =",tint
big: do inran=1,nran
   print "(/,t3,a,i4,/,t3,a,9999(i7))","Calculation number",inran,"Time(ps)",n
   p=p0
   t=0.d0
   tcont=0
   print "(e10.4,9999(i7))",t,p
   do j=1,m
      a(j)=rate(j)*p(re(j))
   enddo
   a0=sum(a)
   do while(a0>0)
     call random_number(rnd)
     t=t-log(rnd(1))/a0      
     r2a0=rnd(2)*a0
     suma=0.d0
     s1: do mu=1,m
       suma=suma+a(mu)
       if(suma>=r2a0) exit s1
     enddo s1
     p(re(mu))=p(re(mu))-1
     p(pr(mu))=p(pr(mu))+1
     cont(mu)=cont(mu)+1
     tcont=tcont+1
     do j=1,m
        a(j)=rate(j)*p(re(j))
     enddo
     a0=sum(a)
     if(mod(tcont,tint)==0.or.a0==0) print "(e10.4,9999(i7))",t,p
   enddo
enddo big
print*,"Population of every species"
do i=1,nesp
   print "(i6,i7)",i,p(i)
enddo
print*,"counts per process"
do i=1,m
   print "(i6,i20)",i,cont(i)
enddo
end program kmc