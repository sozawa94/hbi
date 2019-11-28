module m_matel_ij
  use dtriangular
  contains
  real(8) function matels_ij(i,j,xcol,ycol,zcol,xe,ye,ze)
  integer,intent(in)::i,j
  real(8),intent(in)::xcol(:),ycol(:),zcol(:),xe(:,:),ye(:,:),ze(:,:)
  real(8)::P1(3),P2(3),P3(3)
  real(8)::Sxx,Syy,Szz,Sxy,Sxz,Syz

  P1=(/xe(j,1),ye(j,1),ze(j,1)/)
  P2=(/xe(j,2),ye(j,2),ze(j,2)/)
  P3=(/xe(j,3),ye(j,3),ze(j,3)/)

  call TDstressFS(xcol(i),ycol(i),zcol(i),P1,P2,P3,1.d0,0.d0,0.d0,4.d2,4.d2,&
    & Sxx,Syy,Szz,Sxy,Sxz,Syz)

    matels_ij=Sxy
    return
 end function matels_ij
 ! real(8) function mateln_ij(i,j,xcol,ycol,zcol,xe,ye,ze)
 ! integer,intent(in)::i,j
 ! real(8),intent(in)::xcol(:),ycol(:),zcol(:),xe(:,:),ye(:,:),ze(:,:)
 ! real(8)::P1(3),P2(3),P3(3)
 ! real(8)::Sxx,Syy,Szz,Sxy,Sxz,Syz
 !
 ! P1=(/xe(j,1),ye(j,1),ze(j,1)/)
 ! P2=(/xe(j,2),ye(j,2),ze(j,2)/)
 ! P3=(/xe(j,3),ye(j,3),ze(j,3)/)
 !
 ! call TDstressFS(xcol(i),ycol(i),zcol(i),P1,P2,P3,1.d0,0.d0,0.d0,4.d2,4.d2,&
 !   & Sxx,Syy,Szz,Sxy,Sxz,Syz)
 !
 !   mateln_ij=
 ! end function mateln_ij
 end module m_matel_ij
