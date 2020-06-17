program main
  implicit none
  integer::i,j,imax,jmax,k,kmax
  real(8)::xcol(100000),ycol(100000),xsl(100000),xsr(100000),ysl(100000),ysr(100000),angs(100000),slip(100000)
  real(8)::szx,szy,szx_inc,szy_inc,angle
  real(8)::ds,xloc,yloc
  real(8),parameter::pi=4*atan(1d0)
  character(128)::fname
  ds=0.0025d0
  imax=150
  jmax=100
  kmax=3000

fname='output/slip240.dat'

open(12,file=fname)
do i=1,kmax
 !if(i>10000) ds=0.003d0
  read(12,*) xcol(i),ycol(i),slip(i),angs(i)
  xsl(i)=xcol(i)-0.5d0*ds*cos(angs(i))
  xsr(i)=xcol(i)+0.5d0*ds*cos(angs(i))
  ysl(i)=ycol(i)-0.5d0*ds*sin(angs(i))
  ysr(i)=ycol(i)+0.5d0*ds*sin(angs(i))
  ! xsl(i)=(i-1)*ds
  ! xsr(i)=i*ds
  ! ysl(i)=0d0
  ! ysr(i)=0d0
  ! slip(i)=1d0
end do
close(12)

! kmax=2400
! ds=0.005d0
! angle=0d0/180*pi
! do i=1,2000
!   xsl(i)=cos(angle)*ds*(i-1)
!   xsr(i)=cos(angle)*ds*i
!   ysl(i)=sin(angle)*ds*(i-1)
!   ysr(i)=sin(angle)*ds*i
! end do
! angle=0d0/180*pi
! do i=1,400
!   xsl(2000+i)=cos(angle)*ds*(i-1)+10.0
!   xsr(2000+i)=cos(angle)*ds*i+10.0
!   ysl(2000+i)=sin(angle)*ds*(i-1)-1d0
!   ysr(2000+i)=sin(angle)*ds*i-1d0
! end do

do i=0,imax
  do j=0,jmax
    xloc=-2d0+10d0*(i+0.5)/imax
    yloc=-2d0+10d0*j/jmax
    szx=0d0
    szy=0d0
    do k=1,kmax
      call TY1997B6(xloc,yloc,xsl(k),xsr(k),ysl(k),ysr(k),angs(k),slip(k),szx_inc,szy_inc)
      szx=szx+szx_inc;szy=szy+szy_inc
    end do
    !cfs=
    write(*,'(4e16.5)') xloc,yloc,szx,szy
  end do
  write(*,*)
end do
  stop
contains
subroutine TY1997B6(xloc,yloc,xsl,xsr,ysl,ysr,angs,slip,szx_inc,szy_inc)
  implicit none
  real(8),intent(in)::xloc,yloc,xsl,xsr,ysl,ysr,angs,slip
  real(8),intent(out)::szx_inc,szy_inc
  real(8)::dx(2),dy(2),r(2),gx(2),gy(2),szx(2),szy(2),factor
  real(8),parameter::rigid=40d0,pi=3.140d0,pois=0.25d0

  factor=rigid/(2.d0*pi)
  !dx(1)=xcol(i)-xe(j-1)
  dx(1)=xloc-xsl
  !dy(1)=ycol(i)-ye(j-1)
  dy(1)=yloc-ysl
  r(1)=dsqrt(dx(1)**2+dy(1)**2)
  gx(1)=dx(1)/r(1)
  gy(1)=dy(1)/r(1)

  dx(2)=xloc-xsr
  dy(2)=yloc-ysr
  r(2)=dsqrt(dx(2)**2+dy(2)**2)
  gx(2)=dx(2)/r(2)
  gy(2)=dy(2)/r(2)

  szx(1)=gy(1)/r(1)
  szx(2)=gy(2)/r(2)
  szy(1)=-gx(1)/r(1)
  szy(2)=-gx(2)/r(2)

  szx_inc=(szx(2)-szx(1))*slip*factor
  szy_inc=(szy(2)-szy(1))*slip*factor
  return
end subroutine
end program
