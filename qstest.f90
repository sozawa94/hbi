program main
  implicit none
  integer::i,j,imax,jmax
  real(8)::xsl,xsr,ysl,ysr,disp,qsxx,qsxy,qsyy,xloc,yloc
  xsl=-25d0
  xsr=25d0
  ysl=0d0
  ysr=0d0
  disp=1.d0
  imax=100
  jmax=100
  do i=0,imax
    do j=0,jmax
      xloc=23d0+4d0*dble(i)/imax
      yloc=-2d0+4d0*dble(j)/jmax
      call calc_is(xloc,yloc,xsl,xsr,ysl,ysr,disp,qsxx,qsxy,qsyy)
      write(*,'(5e17.5)') xloc,yloc,qsxx,qsxy,qsyy
    end do
    write(*,*)
  end do
  stop
contains
subroutine calc_is(xloc,yloc,xsl,xsr,ysl,ysr,disp,qsxx,qsxy,qsyy)
  real(8),intent(in)::xloc,yloc,xsl,xsr,ysl,ysr,disp
  real(8),intent(out)::qsxx,qsxy,qsyy
  real(8)::dx(2),dy(2),r(2),gx(2),gy(2),angs,nx,ny,sxx(2),sxy(2),syy(2)
!stressing by vpl creep

!stressing by fault slip

 !dx(1)=xcol(i)-xe(j-1)
dx(1)=xloc-xsl
!dy(1)=ycol(i)-ye(j-1)
dy(1)=yloc-ysl
r(1)=sqrt(dx(1)**2+dy(1)**2)
gx(1)=dx(1)/r(1)
gy(1)=dy(1)/r(1)

dx(2)=xloc-xsr
dy(2)=yloc-ysr
r(2)=sqrt(dx(2)**2+dy(2)**2)
gx(2)=dx(2)/r(2)
gy(2)=dy(2)/r(2)

angs=atan((ysr-ysl)/(xsr-xsl))
nx=-sin(angs)
ny=cos(angs)

sxx(1)=(nx*gx(1)+ny*gy(1)+2*gx(1)*gy(1)*(-ny*gx(1)+nx*gy(1)))/r(1)
sxx(2)=(nx*gx(2)+ny*gy(2)+2*gx(2)*gy(2)*(-ny*gx(2)+nx*gy(2)))/r(2)
sxy(1)=(-gx(1)**2+gy(1)**2)*(-ny*gx(1)+nx*gy(1))/r(1)
sxy(2)=(-gx(2)**2+gy(2)**2)*(-ny*gx(2)+nx*gy(2))/r(2)
syy(1)=(nx*gx(1)+ny*gy(1)-2*gx(1)*gy(1)*(-ny*gx(1)+nx*gy(1)))/r(1)
syy(2)=(nx*gx(2)+ny*gy(2)-2*gx(2)*gy(2)*(-ny*gx(2)+nx*gy(2)))/r(2)
qsxx=sxx(2)-sxx(1)
qsxy=sxy(2)-sxy(1)
qsyy=syy(2)-syy(1)
end subroutine
end program
