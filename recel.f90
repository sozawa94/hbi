module m_matel_ij
  use mod_constant
  !3D static kernel for rectangular elements
contains
  real(8) function matel_ij(i,j,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4,rigid)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),zcol(:),rigid
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)

    matel_ij=ret(xcol(i)-xs1(j),xcol(i)-xs2(j),xcol(i)-xs3(j),xcol(i)-xs4(j),&
    & zcol(i)-zs1(j),zcol(i)-zs2(j),zcol(i)-zs3(j),zcol(i)-zs4(j),rigid)

  end function matel_ij

  real(8) function ret(dx1,dx2,dx3,dx4,dz1,dz2,dz3,dz4,rigid)
    implicit none
    real(8),intent(in)::dx1,dx2,dx3,dx4,dz1,dz2,dz3,dz4,rigid
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
    ret = factor*(ret1+ret2+ret3)
    ! f=1/(1.d0-pois)
    !
    ! a1=f*dz1*(3*dx1**2+2*dz1**2)/(dx1*(dx1**2+dz1**2)**1.5d0)
    ! b1=f*dx1**3/(dz1*(dx1**2+dz1**2)**1.5d0)
    ! c1=g*dx1/(dz1*(dx1**2+dz1**2)**0.5d0)
    !
    ! a2=f*dz2*(3*dx2**2+2*dz2**2)/(dx2*(dx2**2+dz2**2)**1.5d0)
    ! b2=f*dx2**3/(dz2*(dx2**2+dz2**2)**1.5d0)
    ! c2=g*dx2/(dz2*(dx2**2+dz2**2)**0.5d0)
    !
    ! a3=f*dz3*(3*dx3**2+2*dz3**2)/(dx3*(dx3**2+dz3**2)**1.5d0)
    ! b3=f*dx3**3/(dz3*(dx3**2+dz3**2)**1.5d0)
    ! c3=g*dx3/(dz3*(dx3**2+dz3**2)**0.5d0)
    !
    ! a4=f*dz4*(3*dx4**2+2*dz4**2)/(dx4*(dx4**2+dz4**2)**1.5d0)
    ! b1=f*dx1**3/(dz1*(dx1**2+dz1**2)**1.5d0)
    ! c1=g*dx1/(dz1*(dx1**2+dz1**2)**0.5d0)
  end function ret
end module m_matel_ij
