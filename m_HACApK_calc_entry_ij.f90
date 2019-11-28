module m_HACApK_calc_entry_ij
  !use m_matel_ij
  use mod_constant
  type :: st_HACApK_calc_entry
    real*8,pointer :: ao(:)
    character(128)::problem
    integer :: nd,lp61
    real(8),pointer::xcol(:),ycol(:),zcol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    real(8),pointer::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
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
    case('2dp','2dpv')
      HACApK_entry_ij=matels2dp_ij(i,j,st_bemv%xcol,st_bemv%xel,st_bemv%xer)

    case('2dn','2dnv')
      select case(st_bemv%v)

      case('s')
        HACApK_entry_ij=matels2dn_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
        & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang)

      case('n')
        HACApK_entry_ij=mateln2dn_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
        & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang)

      end select

    case('3dp')
      HACApK_entry_ij=matel3dp_ij(i,j,st_bemv%xcol,st_bemv%zcol,&
      & st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%xs4,&
      & st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%zs4)
    end select

  end function HACApK_entry_ij

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

  real(8) function matels2dn_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    !real(8),intent(in)::rigid,pois
    real(8)::xc,yc,xsl,ysl,xsr,ysr,angr,angs
    character(128)::v

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
    v='s'
    call kern(v,xc,yc,xsl,ysl,xsr,ysr,angr,angs,matels2dn_ij)

    !matels_ij=rets(xc,yc,xsl,ysl,xsr,ysr,angr,angs,rigid,pois)
  end function

  real(8) function mateln2dn_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    !real(8),intent(in)::rigid,pois
    real(8)::xc,yc,xsl,ysl,xsr,ysr,angr,angs
    character(128)::v

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
    v='n'
    call kern(v,xc,yc,xsl,ysl,xsr,ysr,angr,angs,mateln2dn_ij)
    !mateln_ij=retn(v,xc,yc,xsl,ysl,xsr,ysr,angr,angs,rigid,pois)
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

   real(8) function matels2dp_ij(i,j,xcol,xel,xer)
     implicit none
     integer,intent(in)::i,j
     real(8),intent(in)::xcol(:),xel(:),xer(:)
     !real(8),intent(in)::rigid,pois
     real(8)::factor

     factor=rigid/(2.d0*pi*(1.d0-pois))
     matels2dp_ij=factor*(1.d0/(xcol(i)-xer(j))-1.d0/(xcol(i)-xel(j)))
   end function matels2dp_ij

endmodule m_HACApK_calc_entry_ij
