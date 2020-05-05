program main
  implicit none
  integer::i,j,imax,jmax,k
  real(8)::xcol(10000),ycol(10000),xsl(10000),xsr(10000),ysl(10000),ysr(10000),angs(10000),slip(10000)
  real(8)::sxx,syy,sxy,sxx_inc,sxy_inc,syy_inc
  real(8)::ds,xloc,yloc
  character(128)::fname
  ds=0.01d0
  imax=140
  jmax=10

fname='output/slip120.dat'
open(12,file=fname)
do i=1,10000
  read(12,*) xcol(i),ycol(i),slip(i),angs(i)
  xsl(i)=xcol(i)-ds*cos(angs(i))
  xsr(i)=xcol(i)+ds*cos(angs(i))
  ysl(i)=ycol(i)-ds*sin(angs(i))
  ysr(i)=ycol(i)+ds*sin(angs(i))
end do
close(12)

do i=0,imax
  do j=0,jmax
    xloc=-20d0+140d0*i/imax
    yloc=-5d0+10d0*j/jmax
    sxx=0d0
    sxy=0d0
    syy=0d0
    do k=1,10000
      call TY1997B11(xloc,yloc,xsl(k),xsr(k),ysl(k),ysr(k),angs(k),slip(k),sxx_inc,sxy_inc,syy_inc)
      sxx=sxx+sxx_inc;syy=syy+syy_inc;sxy=sxy+sxy_inc
    end do
    !cfs=
    write(*,'(5e16.5)') xloc,yloc,sxx,syy,sxy
  end do
  write(*,*)
end do
  stop
contains
subroutine TY1997B11(xloc,yloc,xsl,xsr,ysl,ysr,angs,slip,sxx_inc,sxy_inc,syy_inc)
  implicit none
  real(8),intent(in)::xloc,yloc,xsl,xsr,ysl,ysr,angs,slip
  real(8),intent(out)::sxx_inc,sxy_inc,syy_inc
  real(8)::dx(2),dy(2),r(2),gx(2),gy(2),nx,ny,sxx(2),sxy(2),syy(2),factor
  real(8),parameter::rigid=40d0,pi=3.140d0,pois=0.25d0

  factor=rigid/(2.d0*pi*(1.d0-pois))
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

  nx=-sin(angs)
  ny=cos(angs)


  sxx(1)=(nx*gx(1)+ny*gy(1))/r(1) - 2*gx(1)*gy(1)*(-ny*gx(1)+nx*gy(1))/r(1)
  sxx(2)=(nx*gx(2)+ny*gy(2))/r(2) - 2*gx(2)*gy(2)*(-ny*gx(2)+nx*gy(2))/r(2)
  sxy(1)=(-gx(1)**2+gy(1)**2)*(-ny*gx(1)+nx*gy(1))/r(1)
  sxy(2)=(-gx(2)**2+gy(2)**2)*(-ny*gx(2)+nx*gy(2))/r(2)
  syy(1)=(nx*gx(1)+ny*gy(1))/r(1) + 2*gx(1)*gy(1)*(-ny*gx(1)+nx*gy(1))/r(1)
  syy(2)=(nx*gx(2)+ny*gy(2))/r(2) + 2*gx(2)*gy(2)*(-ny*gx(2)+nx*gy(2))/r(2)
  sxx_inc=(sxx(2)-sxx(1))*slip*factor
  sxy_inc=(sxy(2)-sxy(1))*slip*factor
  syy_inc=(syy(2)-syy(1))*slip*factor
  return
end subroutine
end program
