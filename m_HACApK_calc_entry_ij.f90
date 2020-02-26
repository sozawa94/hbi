module m_HACApK_calc_entry_ij
  !use m_matel_ij
  use mod_constant
  use dtriangular
  type :: st_HACApK_calc_entry
    real*8,pointer :: ao(:)
    character(128)::problem
    integer :: nd,lp61
    real(8),pointer::xcol(:),ycol(:),zcol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    real(8),pointer::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
    real(8),pointer::ys1(:),ys2(:),ys3(:)
    character(128)::v
  end type st_HACApK_calc_entry

  public :: HACApK_entry_ij
contains
  !***HACApK_entry_ij
  real*8 function HACApK_entry_ij(i, j, st_bemv)
    ! zbemv is data needed for calculating the value of the element
    integer::i,j
    type(st_HACApK_calc_entry) :: st_bemv
    select case(st_bemv%problem)
    case('2dp')
      HACApK_entry_ij=matels2dp_ij(i,j,st_bemv%xcol,st_bemv%xel,st_bemv%xer)

    case('2dn_vector')
      HACApK_entry_ij=matel2dn_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang,st_bemv%v)

    case('2dn')
      HACApK_entry_ij=tensor2d_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang,st_bemv%v)
    case('3dp')
      HACApK_entry_ij=matel3dp_ij(i,j,st_bemv%xcol,st_bemv%zcol,&
      & st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%xs4,&
      & st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%zs4)
    case('3dn')
      HACApK_entry_ij=matels1_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%zcol, st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%ys1,st_bemv%ys2,st_bemv%ys3,st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%v)
    end select

  end function HACApK_entry_ij

  real(8) function tensor2d_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang,v)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    character(128),intent(in)::v
    real(8)::xp,xm,yp,ym,kern11,kern12,kern22,sin2,cos2,ds

    xp=cos(ang(j))*(xcol(i)-xel(j))+sin(ang(j))*(ycol(i)-yel(j))
    xm=cos(ang(j))*(xcol(i)-xer(j))+sin(ang(j))*(ycol(i)-yer(j))
    yp=-sin(ang(j))*(xcol(i)-xel(j))+cos(ang(j))*(ycol(i)-yel(j))
    ym=-sin(ang(j))*(xcol(i)-xer(j))+cos(ang(j))*(ycol(i)-yer(j))
    kern11=-0.5d0*rigid/vs*(inte11s(xp,yp)-inte11s(xm,ym))
    kern12=-0.5d0*rigid/vs*(inte12s(xp,yp)-inte12s(xm,ym))
    kern22=-0.5d0*rigid/vs*(inte22s(xp,yp)-inte22s(xm,ym))

    ! xp=dcos(ang(j))*(xcol(i)-xcol(j))+dsin(ang(j))*(ycol(i)-ycol(j))
    ! yp=-dsin(ang(j))*(xcol(i)-xcol(j))+dcos(ang(j))*(ycol(i)-ycol(j))
    ! ds=dsqrt((xer(j)-xel(j))**2+(yer(j)-yel(j))**2)
    ! !write(*,*) xp,yp,ds
    ! kern11=-0.5d0*rigid/vs*(inte11s(xp+0.5d0*ds,yp)-inte11s(xp-0.5d0*ds,yp))
    ! kern12=-0.5d0*rigid/vs*(inte12s(xp+0.5d0*ds,yp)-inte12s(xp-0.5d0*ds,yp))
    ! kern22=-0.5d0*rigid/vs*(inte22s(xp+0.5d0*ds,yp)-inte22s(xp-0.5d0*ds,yp))
    !kern12=-kern12
    !write(*,*) kern11,kern22

    !=>global coordinate system
    sin2=dsin(-2*ang(j))
    cos2=dcos(-2*ang(j))
    select case(v)
    case('xx')
      tensor2d_ij=0.5d0*(kern11+kern22)+0.5d0*(kern11-kern22)*cos2+kern12*sin2
    case('xy')
      tensor2d_ij=-0.5d0*(kern11-kern22)*sin2+kern12*cos2
    case('yy')
      tensor2d_ij=0.5d0*(kern11+kern22)-0.5d0*(kern11-kern22)*cos2-kern12*sin2
    end select
  end function

  real(8) function matel3dp_ij(i,j,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),zcol(:)
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)

    matel3dp_ij=ret3dp(xcol(i)-xs1(j),xcol(i)-xs2(j),xcol(i)-xs3(j),xcol(i)-xs4(j),&
    & zcol(i)-zs1(j),zcol(i)-zs2(j),zcol(i)-zs3(j),zcol(i)-zs4(j))

  end function

  real(8) function ret3dp(dx1,dx2,dx3,dx4,dz1,dz2,dz3,dz4)
    implicit none
    real(8),intent(in)::dx1,dx2,dx3,dx4,dz1,dz2,dz3,dz4
    real(8)::angle,dx,dy,sxx,sxy,syy,gtau(2),gsig(2),factor
    real(8)::dr1,dr2,dr3,dr4,ret1,ret2,ret3
    factor = rigid/(4.d0 * pi)
    dr1=dsqrt(dx1**2+dz1**2)
    dr2=dsqrt(dx2**2+dz2**2)
    dr3=dsqrt(dx3**2+dz3**2)
    dr4=dsqrt(dx4**2+dz4**2)

    ret1=2.d0/3.d0*(dr4/(dx4*dz4)+dr2/(dx2*dz2)-dr3/(dx3*dz3)-dr1/(dx1*dz1))
    ret2=1.d0/3.d0/dx1*(dr4/dz4+dz4/dr4-dr1/dz1-dz1/dr1)
    ret3=1.d0/3.d0/dx2*(dr2/dz2+dz2/dr2-dr3/dz3-dz3/dr3)
    ret3dp = factor*(ret1+ret2+ret3)
  end function

  real(8) function matel2dn_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang,v)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    !real(8),intent(in)::rigid,pois
    real(8)::xc,yc,xsl,ysl,xsr,ysr,angr,angs
    character(128),intent(in)::v

    !write(*,*) i,j

    xc=xcol(i)
    yc=ycol(i)
    !  write(*,*) xc,yc
    xsl=xel(j)
    xsr=xer(j)
    ysl=yel(j)
    ysr=yer(j)

    angr=ang(i)
    angs=ang(j)
    !write(*,*) i,j,xc,yc,angr,angs
    call kern(v,xc,yc,xsl,ysl,xsr,ysr,angr,angs,matel2dn_ij)
    !select case(v)
    !case('s')
    !  matel2dn_ij=rets(xc,yc,xsl,ysl,xsr,ysr,angr,angs)
    !case('n')
    !  matel2dn_ij=retn(xc,yc,xsl,ysl,xsr,ysr,angr,angs)
    !end select
  end function

  subroutine kern(v,xc,yc,xsl,ysl,xsr,ysr,angr,angs,ret)
    implicit none
    real(8),intent(in)::xc,yc,xsl,ysl,xsr,ysr,angr,angs
    character(128),intent(in)::v
    !real(8),parameter::factor=1.d0
    real(8)::dx(4),dy(4),angt(4),r(4),c(4),d(4),gt(4),gn(4)
    real(8)::I1(2),I2(2),nx,ny,gx(2),gy(2),factor
    real(8),intent(out)::ret

    ! rigid=32.04d0
    ! pois=0.25d0
    ! pi=3.14159265d0
    factor=rigid/(2.d0*pi*(1.d0-pois))

    !dx(1)=xcol(i)-xe(j-1)
    dx(1)=xc-xsl
    !dy(1)=ycol(i)-ye(j-1)
    dy(1)=yc-ysl
    r(1)=sqrt(dx(1)**2+dy(1)**2)
    gx(1)=dx(1)/r(1)
    gy(1)=dy(1)/r(1)
    nx=-sin(angr)
    ny=cos(angr)
    !angt(1)=acos(dx(1)/r(1))
    !shear

    !gt(1)=cos(angs-angt(1))
    !gn(1)=-sin(angs-angt(1))
    !c(1)=-cos(2*angr-2*angt(1))
    !d(1)=sin(2*angr-2*angt(1))
    !gtau(1)=factor*c(1)*gt(1)/r(1)
    !gsig(1)=factor*(gn(1)/r(1)+d(1)*gt(1)/r(1))

    dx(2)=xc-xsr
    dy(2)=yc-ysr
    r(2)=sqrt(dx(2)**2+dy(2)**2)
    gx(2)=dx(2)/r(2)
    gy(2)=dy(2)/r(2)

    select case(v)
    case('s')
      I1(1)=-gy(1)/r(1)*(4*nx*ny*gx(1)*gy(1)+(ny**2-nx**2)*(gy(1)**2-gx(1)**2))
      I2(1)=-gx(1)/r(1)*(4*nx*ny*gx(1)*gy(1)+(ny**2-nx**2)*(gy(1)**2-gx(1)**2))
      I1(2)=-gy(2)/r(2)*(4*nx*ny*gx(2)*gy(2)+(ny**2-nx**2)*(gy(2)**2-gx(2)**2))
      I2(2)=-gx(2)/r(2)*(4*nx*ny*gx(2)*gy(2)+(ny**2-nx**2)*(gy(2)**2-gx(2)**2))


      !normal
    case('n')
      I1(1)=gx(1)/r(1)-gy(1)/r(1)*(2*nx*ny*(gy(1)**2-gx(1)**2)-2*(ny**2-nx**2)*gx(1)*gy(1))
      I2(1)=gy(1)/r(1)+gx(1)/r(1)*(2*nx*ny*(gy(1)**2-gx(1)**2)-2*(ny**2-nx**2)*gx(1)*gy(1))
      I1(2)=gx(2)/r(2)-gy(2)/r(2)*(2*nx*ny*(gy(2)**2-gx(2)**2)-2*(ny**2-nx**2)*gx(2)*gy(2))
      I2(2)=gy(2)/r(2)+gx(2)/r(2)*(2*nx*ny*(gy(2)**2-gx(2)**2)-2*(ny**2-nx**2)*gx(2)*gy(2))
    end select

    nx=-sin(angs)
    ny=cos(angs)
    ret=factor*(nx*(-I1(2)+I1(1))+ny*(I2(2)-I2(1)))
  end subroutine

  real(8) function rets(xc,yc,xsl,ysl,xsr,ysr,angr,angs)
    implicit none
    real(8)::xc,yc,xsl,ysl,xsr,ysr,angr,angs,angp
    real(8)::angle,dx,dy,sxx,sxy,syy,gtau(2),gsig(2),factor,dr
    factor=rigid/(2.d0*pi*(1.d0-pois))
    angle=angr-angs
    !dx=xcol(i)-xel(j)
    dx=xc-xsl
    dy=yc-ysl
    dr=dsqrt(dx**2+dy**2)
    !angp=asin(dy/dr)
    !angp=acos(dx/dr)
    !write(*,*)angs,angr,angp
    angp=atan(dy/dx)!*dx/abs(dx)
    if(dy/dx.lt.0d0) angp=atan(abs(dy/dx))
    gtau(1)=-cos(2*angr-2*angp)*cos(angs-angp)/dr


    !sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
    !sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
    !syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
    !write(*,*) sxx,sxy,syy
    !gtau(1)=-0.5d0*(sxx-syy)*dsin(2.d0*angle)+sxy*dcos(2.d0*angle)

    dx=xc-xsr
    dy=yc-ysr
    dr=dsqrt(dx**2+dy**2)
    !angp=asin(dy/dr)
    !angp=acos(dx/dr)
    angp=atan(dy/dx)
    if(dy/dx.lt.0d0) angp=atan(abs(dy/dx))
    gtau(2)=-cos(2*angr-2*angp)*cos(angs-angp)/dr

    !write(*,*) dx,dy
    !sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
    !sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
    !syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
    !write(*,*) sxx,sxy,syy
    !gtau(2)=-0.5d0*(sxx-syy)*dsin(2.d0*angle)+sxy*dcos(2.d0*angle)

    rets = -gtau(2)+gtau(1)
    !write(*,*)rets
  end function rets
  real(8) function retn(xc,yc,xsl,ysl,xsr,ysr,angr,angs)
    implicit none
    real(8)::xc,yc,xsl,ysl,xsr,ysr,angr,angs,angp,dr
    real(8)::angle,dx,dy,sxx,sxy,syy,gtau(2),gsig(2),factor
    factor=rigid/(2.d0*pi*(1.d0-pois))
    angle=angr-angs
    !dx=xcol(i)-xel(j)
    dx=xc-xsl
    dy=yc-ysl
    dr=dsqrt(dx**2+dy**2)
    !angp=asin(dy/dr)
    !if(dx.lt.0d0) angp=acos(-dx/dr)
    angp=atan(dy/dx)
    !write(*,*) angs,angr,angp
    gsig(1)=(-sin(angs-angp)+sin(2*angr-2*angp)*cos(angs-angp))/dr

    !sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
    !sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
    !syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
    !write(*,*) sxx,sxy,syy
    !gsig(1)=0.5d0*(sxx+syy)-0.5d0*(sxx-syy)*dcos(2.d0*angle)-sxy*dsin(2.d0*angle)

    dx=xc-xsr
    dy=yc-ysr
    dr=dsqrt(dx**2+dy**2)
    !angp=acos(dx/dr)
    angp=asin(dy/dr)
    !if(dx.lt.0d0) angp=acos(-dx/dr)
    !angp=atan(dy/dx)
    gsig(2)=(-sin(angs-angp)+sin(2*angr-2*angp)*cos(angs-angp))/dr

    !write(*,*) dx,dy
    !sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
    !sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
    !syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
    !write(*,*) sxx,sxy,syy
    !gsig(2)=0.5d0*(sxx+syy)-0.5d0*(sxx-syy)*dcos(2.d0*angle)-sxy*dsin(2.d0*angle)

    retn = gsig(2)-gsig(1)
  end function retn

  real(8) function matels2dp_ij(i,j,xcol,xel,xer)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),xel(:),xer(:)
    !real(8),intent(in)::rigid,pois
    real(8)::factor

    factor=rigid/(2.d0*pi*(1.d0-pois))
    matels2dp_ij=factor*(1.d0/(xcol(i)-xer(j))-1.d0/(xcol(i)-xel(j)))
  end function matels2dp_ij


  function inte12s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte12s,r
    r=sqrt(x1**2+x2**2)
    pa=vs/vp
    inte12s=2*vs/pi*(1-pa**2)*x1*(x1**2-x2**2)/r**4
    return
  end function

  function inte11s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte11s,r
    r=sqrt(x1**2+x2**2)
    pa=vs/vp
    inte11s=-2*vs/pi*(1-pa**2)*x2*(3*x1**2+x2**2)/r**4
    return
  end function

  function inte22s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte22s,r
    r=sqrt(x1**2+x2**2)
    pa=vs/vp
    inte22s=2*vs/pi*(1-pa**2)*x2*(x1**2-x2**2)/r**4
    return
  end function

  real(8) function matels1_ij(i,j,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,v)
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),zcol(:),xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    character(128),intent(in)::v
    real(8)::P1(3),P2(3),P3(3)
    real(8)::Sxx,Syy,Szz,Sxy,Sxz,Syz

    P1=(/xs1(j),ys1(j),zs1(j)/)
    P2=(/xs2(j),ys2(j),zs2(j)/)
    P3=(/xs3(j),ys3(j),zs3(j)/)

    call TDstressFS(xcol(i),ycol(i),zcol(i),P1,P2,P3,1.d0,0.d0,0.d0,rigid,rigid,&
    & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    select case(v)
    case('xx')
      matels1_ij=Sxx
    case('yy')
      matels1_ij=Syy
    case('zz')
      matels1_ij=Szz
    case('xy')
      matels1_ij=Sxy
    case('xz')
      matels1_ij=Sxz
    case('yz')
      matels1_ij=Syz
    end select
    return
  end function matels1_ij

endmodule m_HACApK_calc_entry_ij
