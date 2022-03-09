program main
  !an example to create mesh file and friction file with rectangular elements
  implicit None
  integer::i,j,k,nb,ios
  real(8)::r,dep,dipangle,a0,b0,dc0,a_max
  integer,parameter::n=9534
  real(8),parameter::pi=4.d0*atan(1.d0)
  real(8),parameter::ds0=1d0,mu0=0.6d0,vref=1d-6,vs=3.464d0,velinit=1d-9,rigid=32.04d0
  real(8)::xcol(n),ycol(n),zcol(n),ang(n),angd(n),rake(n),xr(n),yr(n),zr(n),ds(n)
  real(8)::a(n),b(n),dc(n),f0(n),vel(n),sigma(n),tau(n),mu(n),taudot(n),sigmadot(n)
  real(8)::xs1(n),xs2(n),xs3(n),ys1(n),ys2(n),ys3(n),zs1(n),zs2(n),zs3(n)
  character(128)::dummy

  !geometry
  open(12,file='bp5t.stl',iostat=ios)
  read(12,*)
  !write(*,*) ios
  k=0
  do while(ios==0)
    k=k+1
    read(12,*,iostat=ios)
    read(12,*,iostat=ios)
    read(12,*,iostat=ios) dummy,xs1(k),ys1(k),zs1(k)
    read(12,*,iostat=ios) dummy,xs2(k),ys2(k),zs2(k)
    read(12,*,iostat=ios) dummy,xs3(k),ys3(k),zs3(k)
    read(12,*,iostat=ios)
    read(12,*,iostat=ios)
    xcol(k)=(xs1(k)+xs2(k)+xs3(k))/3
    ycol(k)=(ys1(k)+ys2(k)+ys3(k))/3
    zcol(k)=(zs1(k)+zs2(k)+zs3(k))/3
    if(ios==0)write(*,*) xcol(k),ycol(k),zcol(k)
  end do
  write(*,*)k-1

   !frictional parameters
   open(111,file='bp5tparam.dat')
   a0=0.004d0
   b0=0.03d0
   dc0=0.14d0
   a_max=0.04d0
   do i=1,N
     f0(i)=mu0
     dep=-zcol(i)

     !for BP5
     if((dep.gt.4d0).and.(dep.lt.16d0).and.(abs(xcol(i)).lt.30d0)) then
       a(i)=a0
     else if((dep.lt.2d0).or.(dep.gt.18d0).or.((abs(xcol(i)).gt.32d0))) then
       a(i)=a_max
     else
       r=max(abs(dep-10d0)-6d0,abs(xcol(i))-30d0)/2d0
       a(i)=a0+r*(a_max-a0)
     end if

     r=max(abs(dep-10d0)-6d0,abs(xcol(i))-30d0)/2d0
     a(i)=min(a0+r*(a_max-a0),a_max)
     a(i)=max(a(i),a0)

     !a(i)=0.02d0
     b(i)=b0
     dc(i)=dc0

     if((abs(xcol(i)+24d0).lt.6d0).and.(abs(dep-10d0).lt.6d0)) dc(i)=0.13d0

     !sigma(i)=min(17.0*zcol(i)+10d0,240d0)
     sigma(i)=25d0
     vel(i)=velinit
     dep=-zcol(i)
     if((abs(xcol(i)+24d0).lt.6d0).and.(abs(dep-10d0).lt.6d0)) vel(i)=3d-2
     !mu(i)=f0(i)+(a(i)-b(i))*log(vel(i)/vref)
     mu(i)=a(i)*asinh(0.5d0*vel(i)/vref*exp((f0(i)+b(i)*log(vref/velinit))/a(i)))+rigid/(2*Vs)*vel(i)
     tau(i)=mu(i)*sigma(i)

     taudot(i)=0d0
     sigmadot(i)=0d0

     write(111,'(10e17.8)') 0d0,a(i),b(i),dc(i),f0(i),tau(i),sigma(i),vel(i),taudot(i),sigmadot(i)
     !write(111,'(8e17.8)') 0d0,a(i),b(i),dc(i),f0(i),xcol(i),ycol(i),zcol(i)
     !write(*,*) i,omega
     !if(my_rank.eq.0)write(*,*)phi(i),sigma(i),vel(i)
   end do

  stop
end program
