module dtriangular
  USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  !stop
contains
  !-------------------------------------------------------------------------------
  subroutine TDstressHS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    implicit none
    real(8),intent(in)::X,Y,Z,Ss,Ds,Ts,mu,lambda
    real(8),intent(in)::P1(3),P2(3),P3(3)
    real(8),intent(out)::Sxx,Syy,Szz,Sxy,Sxz,Syz
    real(8)::xx,yy,zz,Q1(3),Q2(3),Q3(3)
    real(8)::Sxx_m,Syy_m,Szz_m,Sxy_m,Sxz_m,Syz_m
    real(8)::Sxx_h,Syy_h,Szz_h,Sxy_h,Sxz_h,Syz_h
    real(8)::Sxx_i,Syy_i,Szz_i,Sxy_i,Sxz_i,Syz_i

    xx=X
    yy=Y
    zz=Z

    call TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    &Sxx_m,Syy_m,Szz_m,Sxy_m,Sxz_m,Syz_m)

    call TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    &Sxx_h,Syy_h,Szz_h,Sxy_h,Sxz_h,Syz_h)

    Q1=P1
    Q2=P2
    Q3=P3
    Q1(3)=-Q1(3)
    Q2(3)=-Q2(3)
    Q3(3)=-Q3(3)
    call TDstressFS(xx,yy,zz,Q1,Q2,Q3,Ss,Ds,Ts,mu,lambda,&
    &Sxx_i,Syy_i,Szz_i,Sxy_i,Sxz_i,Syz_i)

    Sxx=Sxx_m+Sxx_i+Sxx_h
    Syy=Syy_m+Syy_i+Syy_h
    Szz=Szz_m+Szz_i+Szz_h
    Sxy=Sxy_m+Sxy_i+Sxy_h
    Sxz=Sxz_m+Sxz_i+Sxz_h
    Syz=Syz_m+Syz_i+Syz_h

  end subroutine
  subroutine TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    implicit none
    real(8),intent(in)::X, Y, Z, P1(3), P2(3), P3(3), Ss, Ds, Ts,mu,lambda
    real(8),intent(out):: Sxx,Syy,Szz,Sxy,Sxz,Syz
    real(8) :: nu,exx,eyy,ezz,exy,exz,eyz,Exx_,Eyy_,Ezz_,Exy_,Exz_,Eyz_
    integer :: Trimode
    logical :: casepLog, casenLog, casezLog
    real(8) :: A, B, C, bx, by, bz, Fi, na, nb, nc, &
    & u, u1Tn, u2Tn, u3Tn, u1Tp, u2Tp, u3Tp, v, v1Tn, v2Tn, v3Tn, v1Tp, v2Tp, v3Tp, w, w1Tn, w2Tn, w3Tn, w1Tp, w2Tp, w3Tp, &
    & x_, xn, xp, y_, yn, yp, z_, zn, zp
    REAL*8, DIMENSION(3) :: a_, b_, c_  ! N.B. "_" = "lower-case" (because Fortran is NOT case-sensitive).
    REAL*8, DIMENSION(3) :: e12, e13, e23, eY, eZ, P2mP1, P3mP1, p_1, p_2, p_3, Vnorm, Vstrike, Vdip
    REAL*8, DIMENSION(3, 3) :: At, matrix3x3
    real(8)::Exx1Tp,Eyy1Tp,Ezz1Tp,Exy1Tp,Exz1Tp,Eyz1Tp
    real(8)::Exx2Tp,Eyy2Tp,Ezz2Tp,Exy2Tp,Exz2Tp,Eyz2Tp
    real(8)::Exx3Tp,Eyy3Tp,Ezz3Tp,Exy3Tp,Exz3Tp,Eyz3Tp
    real(8)::Exx1Tn,Eyy1Tn,Ezz1Tn,Exy1Tn,Exz1Tn,Eyz1Tn
    real(8)::Exx2Tn,Eyy2Tn,Ezz2Tn,Exy2Tn,Exz2Tn,Eyz2Tn
    real(8)::Exx3Tn,Eyy3Tn,Ezz3Tn,Exy3Tn,Exz3Tn,Eyz3Tn! TDstressFS


    nu = 1/(1+lambda/mu)/2 ! Poisson's ratio
    !!!write(*,*) 'calc'

    bx = Ts ! Tensile-slip
    by = Ss ! Strike-slip
    bz = Ds ! Dip-slip

    P2mP1 = P2 - P1 ! all 3 components
    P3mP1 = P3 - P1
    CALL DCross(P2mP1, P3mP1, Vnorm) ! but, this still needs to be normalized:
    Vnorm=Vnorm/sqrt(dot_product(Vnorm,Vnorm))

    eY = (/ 0.0D0, 1.0D0, 0.0D0 /)
    eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
    CALL DCross(eZ, Vnorm, Vstrike)
    Vstrike=Vstrike/sqrt(dot_product(Vstrike,Vstrike))

    IF (dot_product(Vstrike,Vstrike) == 0.0D0) THEN
      Vstrike = eY * Vnorm(3)
    END IF
    Vstrike=Vstrike/sqrt(dot_product(Vstrike,Vstrike))
    CALL DCross(Vnorm, Vstrike, Vdip)
    p_1 = 0.0D0
    p_2 = 0.0D0
    p_3 = 0.0D0

    At(1, 1:3) = Vnorm(1:3)
    At(2, 1:3) = Vstrike(1:3)
    At(3, 1:3) = Vdip(1:3)

    CALL CoordTrans(X-P2(1), Y-P2(2), Z-P2(3), At, & ! inputs
    & x_, y_, z_)                      ! outputs
    CALL CoordTrans(P1(1)-P2(1), P1(2)-P2(2), P1(3)-P2(3), At, & ! inputs
    & p_1(1), p_1(2), p_1(3))                      ! outputs
    CALL CoordTrans(P3(1)-P2(1), P3(2)-P2(2), P3(3)-P2(3), At, & ! inputs
    & p_3(1), p_3(2), p_3(3))                      ! outputs

    e12 = p_2 - p_1
    e12=e12/sqrt(dot_product(e12,e12))
    e13 = p_3 - p_1
    e13=e13/sqrt(dot_product(e13,e13))
    e23 = p_3 - p_2
    e23=e23/sqrt(dot_product(e23,e23))



    A = ACOS(Dot_product(e12, e13))
    B = ACOS(Dot_product(-e12, e23))
    C = ACOS(Dot_product(e23, e13))

    !Ts,Ss,Ds-->bx,by,bz
    !!!write(*,*) At
    !CALL CoordTrans(Ts, Ss, Ds, At,&
    !& bx,by,bz)
    !!write(*,*) 'bx,by,bz',bx,by,bz


    CALL trimodefinder(y_, z_, x_, p_1(2:3), p_2(2:3), p_3(2:3), & ! inputs
    & Trimode)
    !write(*,*)Trimode                                 ! output
    casepLog = (Trimode == 1)
    casenLog = (Trimode == -1)
    casezLog = (Trimode == 0)

    IF (casepLog) THEN
      xp = x_
      yp = y_
      zp = z_
    END IF
    IF (casenLog) THEN
      xn = x_
      yn = y_
      zn = z_
    END IF
    IF (casezLog) THEN
      xn = x_
      yn = y_
      zn = z_
    END IF

    IF (casepLog) THEN
      CALL TDSetupD(xp, yp, zp, A, bx, by, bz, nu, p_1, -e13, & ! inputs
      & Exx1Tp,Eyy1Tp,Ezz1Tp,Exy1Tp,Exz1Tp,Eyz1Tp)
      CALL TDSetupD(xp, yp, zp, B, bx, by, bz, nu, p_2, e12, & ! inputs
      & Exx2Tp,Eyy2Tp,Ezz2Tp,Exy2Tp,Exz2Tp,Eyz2Tp)
      CALL TDSetupD(xp, yp, zp, C, bx, by, bz, nu, p_3, e23, & ! inputs
      & Exx3Tp,Eyy3Tp,Ezz3Tp,Exy3Tp,Exz3Tp,Eyz3Tp)
    END IF

    IF (caseNLog) THEN
      CALL TDSetupD(xn, yn, zn, A, bx, by, bz, nu, p_1, e13, & ! inputs
      & Exx1Tn,Eyy1Tn,Ezz1Tn,Exy1Tn,Exz1Tn,Eyz1Tn)
      CALL TDSetupD(xn, yn, zn, B, bx, by, bz, nu, p_2, -e12, & ! inputs
      & Exx2Tn,Eyy2Tn,Ezz2Tn,Exy2Tn,Exz2Tn,Eyz2Tn)
      CALL TDSetupD(xn, yn, zn, C, bx, by, bz, nu, p_3, -e23, & ! inputs
      & Exx3Tn,Eyy3Tn,Ezz3Tn,Exy3Tn,Exz3Tn,Eyz3Tn)
    END IF

     IF (casezLog) THEN
      CALL TDSetupD(xn, yn, zn, A, bx, by, bz, nu, p_1, e13, & ! inputs
      & Exx1Tn,Eyy1Tn,Ezz1Tn,Exy1Tn,Exz1Tn,Eyz1Tn)
      CALL TDSetupD(xn, yn, zn, B, bx, by, bz, nu, p_2, -e12, & ! inputs
      & Exx2Tn,Eyy2Tn,Ezz2Tn,Exy2Tn,Exz2Tn,Eyz2Tn)
      CALL TDSetupD(xn, yn, zn, C, bx, by, bz, nu, p_3, -e23, & ! inputs
      & Exx3Tn,Eyy3Tn,Ezz3Tn,Exy3Tn,Exz3Tn,Eyz3Tn)
    END IF

    ! Calculate the strain tensor components in TDCS
    if (casepLog) then
      exx = Exx1Tp+Exx2Tp+Exx3Tp
      eyy = Eyy1Tp+Eyy2Tp+Eyy3Tp
      ezz = Ezz1Tp+Ezz2Tp+Ezz3Tp
      exy = Exy1Tp+Exy2Tp+Exy3Tp
      exz = Exz1Tp+Exz2Tp+Exz3Tp
      eyz = Eyz1Tp+Eyz2Tp+Eyz3Tp
    end if
    if (casenLog) then
      exx = Exx1Tn+Exx2Tn+Exx3Tn
      eyy = Eyy1Tn+Eyy2Tn+Eyy3Tn
      ezz = Ezz1Tn+Ezz2Tn+Ezz3Tn
      exy = Exy1Tn+Exy2Tn+Exy3Tn
      exz = Exz1Tn+Exz2Tn+Exz3Tn
      eyz = Eyz1Tn+Eyz2Tn+Eyz3Tn
    end if
    if (casezLog) then
      exx = Exx1Tn+Exx2Tn+Exx3Tn
      eyy = Eyy1Tn+Eyy2Tn+Eyy3Tn
      ezz = Ezz1Tn+Ezz2Tn+Ezz3Tn
      exy = Exy1Tn+Exy2Tn+Exy3Tn
      exz = Exz1Tn+Exz2Tn+Exz3Tn
      eyz = Eyz1Tn+Eyz2Tn+Eyz3Tn
    end if
    if (casezLog) then
      exx = IEEE_VALUE(u, IEEE_QUIET_NAN)
      eyy = IEEE_VALUE(u, IEEE_QUIET_NAN)
      ezz = IEEE_VALUE(u, IEEE_QUIET_NAN)
      exy = IEEE_VALUE(u, IEEE_QUIET_NAN)
      exz = IEEE_VALUE(u, IEEE_QUIET_NAN)
      eyz = IEEE_VALUE(u, IEEE_QUIET_NAN)
    end if

    !Sxx = 2*mu*exx+lambda*(exx+eyy+ezz)
    !Syy = 2*mu*eyy+lambda*(exx+eyy+ezz)
    !Szz = 2*mu*ezz+lambda*(exx+eyy+ezz)
    !Sxy = 2*mu*exy
    !Sxz = 2*mu*exz
    !Syz = 2*mu*eyz
    !!write(*,*)'Sxx,Syy,Szz,Sxy,Sxz,Syz in TDCS'
    !write(*,*)Sxx,Syy,Szz,Sxy,Sxz,Syz


    ! Transform the strain tensor components from TDCS into EFCS
    !write(*,*)At
    call TensTrans(exx,eyy,ezz,exy,exz,eyz,At,&
    & Exx_,Eyy_,Ezz_,Exy_,Exz_,Eyz_)

    ! Calculate the stress tensor components in EFCS
    Sxx = 2*mu*Exx_+lambda*(Exx_+Eyy_+Ezz_)
    Syy = 2*mu*Eyy_+lambda*(Exx_+Eyy_+Ezz_)
    Szz = 2*mu*Ezz_+lambda*(Exx_+Eyy_+Ezz_)
    Sxy = 2*mu*Exy_
    Sxz = 2*mu*Exz_
    Syz = 2*mu*Eyz_

    !!!write(*,*)Sxx,Syy,Szz,Sxy,Sxz,Syz

    !!!write(*,*)Sxy

  end subroutine
  subroutine TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    &Sxx,Syy,Szz,Sxy,Sxz,Syz)
    implicit none
    real(8),intent(in)::X,Y,Z,P1(3),P2(3),P3(3),Ss,Ds,Ts,mu,lambda
    real(8),intent(out)::Sxx,Syy,Szz,Sxy,Sxz,Syz
    real(8)::bx,by,bz,P2mP1(3),P3mP1(3),Vstrike(3),Vnorm(3),Vdip(3),eY(3),eZ(3)
    real(8)::At(3,3),bxt,byt,bzt
    real(8)::Sxx1,Syy1,Szz1,Sxy1,Sxz1,Syz1,Sxx2,Syy2,Szz2,Sxy2,Sxz2,Syz2,Sxx3,Syy3,Szz3,Sxy3,Sxz3,Syz3

  bx = Ts!; % Tensile-slip
  by = Ss!; % Strike-slip
  bz = Ds!; % Dip-slip

  P2mP1 = P2 - P1
  P3mP1 = P3 - P1
  !write(*,*) P2mP1,P3mP1
  CALL DCross(P2mP1, P3mP1, Vnorm)
  Vnorm=Vnorm/sqrt(dot_product(Vnorm,Vnorm))
  !write(*,*)Vnorm
  eY = (/ 0.0D0, 1.0D0, 0.0D0 /)
  eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
  CALL DCross(eZ, Vnorm, Vstrike)
  Vstrike=Vstrike/sqrt(dot_product(Vstrike,Vstrike))

  IF (dot_product(Vstrike,Vstrike) == 0.0D0) THEN
    Vstrike = eY * Vnorm(3)
  END IF
  Vstrike=Vstrike/sqrt(dot_product(Vstrike,Vstrike))
  !write(*,*)Vstrike
  CALL DCross(Vnorm, Vstrike, Vdip)
  At(1:3, 1) = Vnorm(1:3)
  At(1:3, 2) = Vstrike(1:3)
  At(1:3, 3) = Vdip(1:3)
  CALL CoordTrans(bx, by, bz, At, & ! inputs
  & bxt, byt, bzt)

  call AngSetupFSC_S(X,Y,Z,bXt,bYt,bZt,P1,P2,mu,lambda,&
  &Sxx1,Syy1,Szz1,Sxy1,Sxz1,Syz1)
  !write(*,*)Sxx1,Syy1,Szz1,Sxy1,Sxz1,Syz1
  ! write(*,*)
  call AngSetupFSC_S(X,Y,Z,bXt,bYt,bZt,P2,P3,mu,lambda,&
  &Sxx2,Syy2,Szz2,Sxy2,Sxz2,Syz2)
  !write(*,*)Sxx2,Syy2,Szz2,Sxy2,Sxz2,Syz2
  ! write(*,*)
  call AngSetupFSC_S(X,Y,Z,bXt,bYt,bZt,P3,P1,mu,lambda,&
  &Sxx3,Syy3,Szz3,Sxy3,Sxz3,Syz3)
  !write(*,*)Sxx3,Syy3,Szz3,Sxy3,Sxz3,Syz3
  !write(*,'(4f15.9)') Sxy1,Sxy2,Sxy3,Sxy1+Sxy2+Sxy3
  !
  Sxx=Sxx1+Sxx2+Sxx3
  Syy=Syy1+Syy2+Syy3
  Szz=Szz1+Szz2+Szz3
  Sxy=Sxy1+Sxy2+Sxy3
  Sxz=Sxz1+Sxz2+Sxz3
  Syz=Syz1+Syz2+Syz3
  !write(*,'(3f15.9)') Sxy,Sxz,Syz
  return
  end subroutine
  !------------------------------------------------------------------------------
  subroutine TensTrans(Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1,A,&
    & Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2)
    implicit none
    REAL*8, INTENT(IN) :: Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1
    REAL*8, DIMENSION(3, 3) :: A
    REAL*8, INTENT(OUT) :: Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2

    Txx2 = A(1,1)**2*Txx1+2*A(1,1)*A(2,1)*Txy1+2*A(1,1)*A(3,1)*Txz1+2*A(2,1)*A(3,1)*Tyz1+&
    &  A(2,1)**2*Tyy1+A(3,1)**2*Tzz1
    Tyy2 = A(1,2)**2*Txx1+2*A(1,2)*A(2,2)*Txy1+2*A(1,2)*A(3,2)*Txz1+2*A(2,2)*A(3,2)*Tyz1+&
    &  A(2,2)**2*Tyy1+A(3,2)**2*Tzz1
    Tzz2 = A(1,3)**2*Txx1+2*A(1,3)*A(2,3)*Txy1+2*A(1,3)*A(3,3)*Txz1+2*A(2,3)*A(3,3)*Tyz1+&
    &  A(2,3)**2*Tyy1+A(3,3)**2*Tzz1
    Txy2 = A(1,1)*A(1,2)*Txx1+(A(1,1)*A(2,2)+A(2,1)*A(1,2))*Txy1+(A(1,1)*A(3,2)+&
    &  A(1,2)*A(3,1))*Txz1+(A(3,2)*A(2,1)+A(3,1)*A(2,2))*Tyz1+A(2,2)*A(2,1)*Tyy1+A(3,1)*A(3,2)*Tzz1
    Txz2 = A(1,1)*A(1,3)*Txx1+(A(1,1)*A(2,3)+A(1,3)*A(2,1))*Txy1+(A(1,1)*A(3,3)+&
    &  A(3,1)*A(1,3))*Txz1+(A(3,3)*A(2,1)+A(3,1)*A(2,3))*Tyz1+A(2,3)*A(2,1)*Tyy1+A(3,1)*A(3,3)*Tzz1
    Tyz2 = A(1,2)*A(1,3)*Txx1+(A(1,3)*A(2,2)+A(1,2)*A(2,3))*Txy1+(A(1,3)*A(3,2)+&
    &  A(1,2)*A(3,3))*Txz1+(A(2,3)*A(3,2)+A(3,3)*A(2,2))*Tyz1+A(2,2)*A(2,3)*Tyy1+A(3,2)*A(3,3)*Tzz1
  end subroutine
  !------------------------------------------------------------------------------
  SUBROUTINE CoordTrans(x1, x2, x3, A, &
    & Xcap1, Xcap2, Xcap3)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, x2, x3
    REAL*8, DIMENSION(3, 3) :: A
    REAL*8, INTENT(OUT) :: Xcap1, Xcap2, Xcap3

    Xcap1 = A(1, 1) * x1 + A(1, 2) * x2 + A(1, 3) * x3
    Xcap2 = A(2, 1) * x1 + A(2, 2) * x2 + A(2, 3) * x3
    Xcap3 = A(3, 1) * x1 + A(3, 2) * x2 + A(3, 3) * x3

  END SUBROUTINE CoordTrans
  !------------------------------------------------------------------------------
  SUBROUTINE trimodefinder(x, y, z, p1, p2, p3, & ! inputs
    & trimode)               ! output
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x, y, z
    REAL*8, DIMENSION(2), INTENT(IN) :: p1, p2, p3
    INTEGER, INTENT(OUT) :: trimode
    REAL*8 :: a, b, c
    a = ((p2(2)-p3(2))*(x-p3(1))+(p3(1)-p2(1))*(y-p3(2)))/  &
    & ((p2(2)-p3(2))*(p1(1)-p3(1))+(p3(1)-p2(1))*(p1(2)-p3(2)))
    b = ((p3(2)-p1(2))*(x-p3(1))+(p1(1)-p3(1))*(y-p3(2)))/  &
    & ((p2(2)-p3(2))*(p1(1)-p3(1))+(p3(1)-p2(1))*(p1(2)-p3(2)))
    c = 1.0D0 - a - b

    trimode=1
    IF ((a <= 0.0D0) .AND. (b > c) .AND. (c > a)) trimode = -1
    IF ((b <= 0.0D0) .AND. (c > a) .AND. (a > b)) trimode = -1
    IF ((c <= 0.0D0) .AND. (a > b) .AND. (b > c)) trimode = -1
    IF ((a == 0.0D0) .AND. (b >= 0.0D0) .AND. (c >= 0.0D0)) trimode = 0
    IF ((a >= 0.0D0) .AND. (b == 0.0D0) .AND. (c >= 0.0D0)) trimode = 0
    IF ((a >= 0.0D0) .AND. (b >= 0.0D0) .AND. (c == 0.0D0)) trimode = 0
    IF ((trimode == 0) .AND. (z /= 0.0D0)) trimode = 1
  end subroutine
  !------------------------------------------------------------------------------
  SUBROUTINE TDSetupD(x, y, z, alpha, bx, by, bz, nu, TriVertex, SideVec, & ! inputs
    & exx,eyy,ezz,exy,exz,eyz)                              ! outputs
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: alpha, bx, by, bz, nu, x, y, z
    REAL*8, DIMENSION(3), INTENT(IN) :: SideVec, TriVertex
    REAL*8, INTENT(OUT) :: exx,eyy,ezz,exy,exz,eyz
    REAL*8 :: exx_,eyy_,ezz_,exy_,exz_,eyz_
    real(8),parameter::pi=4*atan(1.d0)

    REAL*8 :: by1, bz1, v0, w0, y1, z1
    REAL*8, DIMENSION(2) :: r1, r2, r3
    REAL*8, DIMENSION(2, 2) :: A
    REAL*8, DIMENSION(3, 3) :: B

    ! Transformation matrix
    A(1, 1) = SideVec(3)
    A(1, 2) = -SideVec(2)
    A(2, 1) = SideVec(2)
    A(2, 2) = SideVec(3)

    ! Transform coordinates of the calculation points from TDCS into ADCS
    r1(1) = A(1, 1) * (y - TriVertex(2)) + A(1, 2) * (z -TriVertex(3))
    r1(2) = A(2, 1) * (y - TriVertex(2)) + A(2, 2) * (z -TriVertex(3))
    y1 = r1(1)
    z1 = r1(2)

    r2(1) = A(1, 1) * by + A(1, 2) * bz
    r2(2) = A(2, 1) * by + A(2, 2) * bz
    by1 = r2(1)
    bz1 = r2(2)

    ! Calculate strains associated with an angular dislocation in ADCS
    call AngDisStrain(x, y1, z1, -pi+alpha, bx, by1, bz1, nu, & ! inputs
    & Exx,Eyy,Ezz,Exy,Exz,Eyz) ! outputs

    B=0
    B(1,1)=1
    B(2,2)=A(1,1)
    B(2,3)=A(1,2)
    B(3,2)=A(2,1)
    B(3,3)=A(2,2)

    exx_=exx
    eyy_=eyy
    ezz_=ezz
    exy_=exy
    exz_=exz
    eyz_=eyz
    call TensTrans(exx_,eyy_,ezz_,exy_,exz_,eyz_,B,& !inputs
    & Exx,Eyy,Ezz,Exy,Exz,Eyz) !outputs
    return
  end subroutine
  ! function [Stress,Strain]=AngSetupFSC_S(X,Y,Z,bX,bY,bZ,PA,PB,mu,lambda)
  subroutine AngSetupFSC_S(X,Y,Z,bX,bY,bZ,PA,PB,mu,lambda,Sxx,Syy,Szz,Sxy,Sxz,Syz)
    implicit none
    real(8),intent(in)::X,Y,Z,bX,bY,bZ,PA(3),PB(3),mu,lambda
    real(8),intent(out)::Sxx,Syy,Szz,Sxy,Sxz,Syz
    logical::configI
    real(8)::A(3,3),nu,SIdeVec(3),eZ(3),beta,ey1(3),ey2(3),ey3(3)
    real(8)::y1A,y2A,y3A,y1AB,y2AB,y3AB,y1B,y2B,y3B,b1,b2,b3,bxt,byt,bzt
    real(8)::v11A,v22A,v33A,v12A,v13A,v23A,v11B,v22B,v33B,v12B,v13B,v23B
    real(8)::v11,v22,v33,v12,v13,v23,exx,eyy,ezz,exy,exz,eyz
    real(8),parameter::pi=4*atan(1d0)
  nu = lambda/(lambda+mu)/2!; % Poisson's ratio
  !
  SideVec = PB-PA
  !write(*,*)
  !write(*,*)'SideVec',SideVec
  eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
  beta = acos(-SideVec(3)/sqrt(dot_product(SideVec,SideVec)))!;
  !write(*,*)beta
  !
   Sxx=0d0;Syy=0d0;Szz=0d0;Sxy=0d0;Sxz=0d0;Syz=0d0
  if ((abs(beta)<1d-3).or.(abs(pi-beta)<1d-3)) then
  !     Stress = zeros(length(X),6);
  return
  !Sxx=0d0;Syy=0d0;Szz=0d0;Sxy=0d0;Sxz=0d0;Syz=0d0
  !     Strain = zeros(length(X),6);
  else
       ey1 = (/SideVec(1),SideVec(2),0d0/)

       ey1 = ey1/sqrt(dot_product(ey1,ey1))
       !write(*,*) 'ey1',ey1
       ey3 = -eZ
       !write(*,*) 'ey3',ey3
       call Dcross(ey3,ey1,ey2)
       !write(*,*) 'ey2',ey2
       !ey2 = cross(ey3,ey1);
       !A = [ey1,ey2,ey3]!; % Transformation matrix
       A(1,1:3)=ey1(1:3)
       A(2,1:3)=ey2(1:3)
       A(3,1:3)=ey3(1:3)
  !
  !     % Transform coordinates from EFCS to the first ADCS
  !     [y1A,y2A,y3A] = CoordTrans(X-PA(1),Y-PA(2),Z-PA(3),A);
  CALL CoordTrans(X-PA(1),Y-PA(2),Z-PA(3), A, & ! inputs
  & y1A,y2A,y3A)
  !if((abs(beta)<1d-1).or.(abs(pi-beta)<1d-1)) then
  !if(abs(y1A*Y1A+Y2A*Y2A).le.1d-6.and.(abs(beta)<1d-3.or.abs(pi-beta)>1d-3))return
  !end if
  !write(*,*)'y1A,y2A,y3A',y1A,y2A,y3A
  !     % Transform coordinates from EFCS to the second ADCS
  CALL CoordTrans(SideVec(1),SideVec(2),SideVec(3), A, & ! inputs
  & y1AB,y2AB,y3AB)
  !     [y1AB,y2AB,y3AB] = CoordTrans(SideVec(1),SideVec(2),SideVec(3),A);
       y1B = y1A-y1AB
       y2B = y2A-y2AB
       y3B = y3A-y3AB
  !write(*,*) 'y1B,y2B,y3B',y1B,y2B,y3B
  !if(abs(y1B*Y1B+y2B*y2B).le.1d-6.and.(abs(beta)<1d-3.or.abs(pi-beta)<1d-3)) return
  !     % Transform slip vector components from EFCS to ADCS
  !     [b1,b2,b3] = CoordTrans(bX,bY,bZ,A);
  CALL CoordTrans(bx, by, bz, A, & ! inputs
  & b1, b2, b3)
  !write(*,*)'b1,b2,b3',b1,b2,b3
  !
  !     % Determine the best artefact-free configuration for the calculation
  !     % points near the free surface
  !     I = (beta*y1A)>=0;
  configI=.false.
  if (beta*y1A.ge.0d0) configI=.true.
  !write(*,*)beta*y1A,configI
  !     % For singularities at surface
       v11A = 0d0
       v22A = 0d0
       v33A = 0d0
       v12A = 0d0
       v13A = 0d0
       v23A = 0d0
  !
       v11B = 0d0
       v22B = 0d0
       v33B = 0d0
       v12B = 0d0
       v13B = 0d0
       v23B = 0d0
  !
  !     % Configuration I
  !     [v11A(I),v22A(I),v33A(I),v12A(I),v13A(I),v23A(I)] = ...
  !         AngDisStrainFSC(y1A(I),y2A(I),y3A(I),-pi+beta,b1,b2,b3,nu,-PA(3));
  if(configI) then
  call AngDisStrainFSC(y1A,y2A,y3A,-pi+beta,b1,b2,b3,nu,-PA(3),&
  &v11A,v22A,v33A,v12A,v13A,v23A)
  !write(*,*)'v11A,v22A,v33A,v12A,v13A,v23A'
  !write(*,*)v11A,v22A,v33A,v12A,v13A,v23A
  !
  !     [v11B(I),v22B(I),v33B(I),v12B(I),v13B(I),v23B(I)] = ...
  !         AngDisStrainFSC(y1B(I),y2B(I),y3B(I),-pi+beta,b1,b2,b3,nu,-PB(3));
  call AngDisStrainFSC(y1B,y2B,y3B,-pi+beta,b1,b2,b3,nu,-PB(3),&
  &v11B,v22B,v33B,v12B,v13B,v23B)
  !write(*,*)'v11B,v22B,v33B,v12B,v13B,v23B'
  !write(*,*)v11B,v22B,v33B,v12B,v13B,v23B

  else
  !     % Configuration II
  !     [v11A(~I),v22A(~I),v33A(~I),v12A(~I),v13A(~I),v23A(~I)] = ...
  !         AngDisStrainFSC(y1A(~I),y2A(~I),y3A(~I),beta,b1,b2,b3,nu,-PA(3));
  call AngDisStrainFSC(y1A,y2A,y3A,beta,b1,b2,b3,nu,-PA(3),&
  &v11A,v22A,v33A,v12A,v13A,v23A)
  !write(*,*)'v11A,v22A,v33A,v12A,v13A,v23A'
  !write(*,*)v11A,v22A,v33A,v12A,v13A,v23A
  !
  !     [v11B(~I),v22B(~I),v33B(~I),v12B(~I),v13B(~I),v23B(~I)] = ...
  !         AngDisStrainFSC(y1B(~I),y2B(~I),y3B(~I),beta,b1,b2,b3,nu,-PB(3));
  call AngDisStrainFSC(y1B,y2B,y3B,beta,b1,b2,b3,nu,-PB(3),&
  &v11B,v22B,v33B,v12B,v13B,v23B)
  !write(*,*)'v11B,v22B,v33B,v12B,v13B,v23B'
  !write(*,*)v11B,v22B,v33B,v12B,v13B,v23B
  end if
  !     % Calculate total Free Surface Correction to strains in ADCS
      v11 = v11B-v11A
      v22 = v22B-v22A
       v33 = v33B-v33A
       v12 = v12B-v12A
       v13 = v13B-v13A
       v23 = v23B-v23A
       !rite(*,*)v11,v22,v33,v12,v13,v23
  !
  !     % Transform total Free Surface Correction to strains from ADCS to EFCS
  !     [Exx,Eyy,Ezz,Exy,Exz,Eyz] = TensTrans(v11,v22,v33,v12,v13,v23,A');
  !debug
  ! A(1:3,1)=ey1(1:3)
  ! A(1:3,2)=ey2(1:3)
  ! A(1:3,3)=ey3(1:3)
  !debug end

  call TensTrans(v11,v22,v33,v12,v13,v23,A,&
    & Exx,Eyy,Ezz,Exy,Exz,Eyz)
    !write(*,*)Exx,Eyy,Ezz,Exy,Exz,Eyz
  !
  !     % Calculate total Free Surface Correction to stresses in EFCS
       Sxx = 2*mu*Exx+lambda*(Exx+Eyy+Ezz)
       Syy = 2*mu*Eyy+lambda*(Exx+Eyy+Ezz)
       Szz = 2*mu*Ezz+lambda*(Exx+Eyy+Ezz)
       Sxy = 2*mu*Exy
       Sxz = 2*mu*Exz
       Syz = 2*mu*Eyz
  !
  !     Strain = [Exx,Eyy,Ezz,Exy,Exz,Eyz];
  !     Stress = [Sxx,Syy,Szz,Sxy,Sxz,Syz];
  ! end
  end if
  end subroutine
  !------------------------------------------------------------------------------
  SUBROUTINE AngDisStrain(x, y, z, alpha, bx, by, bz, nu, & ! inputs
    & Exx,Eyy,Ezz,Exy,Exz,Eyz) ! outputs

    IMPLICIT NONE
    REAL*8, PARAMETER :: pi = 4*atan(1.d0)
    REAL*8, INTENT(IN) :: x, y, z, alpha, bx, by, bz, nu
    REAL*8, INTENT(OUT) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
    REAL*8 :: cosA, eta, r, sinA, ux, uy, uz, vx, vy, vz, wx, wy, wz, zz, zeta
    Real*8 :: C,r2,r2z2,r3,r3z,rFi_rx,rFi_ry,rFi_rz,rz,S,W,W2,W2r,W2r2,Wr,Wr3,x2,y2,z2

    sinA = dsin(alpha)
    cosA = dcos(alpha)
    eta = y*cosA-z*sinA
    zeta = y*sinA+z*cosA

    x2 = x**2
    y2 = y**2
    z2 = z**2
    r2 = x2+y2+z2
    r = dsqrt(r2)
    r3 = r*r2
    rz = r*(r-z)
    r2z2 = r2*(r-z)**2
    r3z = r3*(r-z)

    W = zeta-r
    W2 = W**2
    Wr = W*r
    W2r = W2*r
    Wr3 = W*r3
    W2r2 = W2*r2

    C = (r*cosA-z)/Wr
    S = (r*sinA-y)/Wr

    ! Partial derivatives of the Burgers' function
    rFi_rx = (eta/r/(r-zeta)-y/r/(r-z))/4/pi
    rFi_ry = (x/r/(r-z)-cosA*x/r/(r-zeta))/4/pi
    rFi_rz = (sinA*x/r/(r-zeta))/4/pi

    Exx = bx*(rFi_rx)+bx/8/pi/(1-nu)*(eta/Wr+eta*x2/W2r2-eta*x2/Wr3+y/rz-&
    & x2*y/r2z2-x2*y/r3z)-by*x/8/pi/(1-nu)*(((2*nu+1)/Wr+x2/W2r2-x2/Wr3)*cosA+&
    & (2*nu+1)/rz-x2/r2z2-x2/r3z)+bz*x*sinA/8/pi/(1-nu)*((2*nu+1)/Wr+x2/W2r2-x2/Wr3)

    Eyy = by*(rFi_ry)+bx/8/pi/(1-nu)*((1/Wr+S**2-y2/Wr3)*eta+(2*nu+1)*y/rz-y**3/r2z2-&
    & y**3/r3z-2*nu*cosA*S)-by*x/8/pi/(1-nu)*(1/rz-y2/r2z2-y2/r3z+&
    & (1/Wr+S**2-y2/Wr3)*cosA)+bz*x*sinA/8/pi/(1-nu)*(1/Wr+S**2-y2/Wr3)

    Ezz = bz*(rFi_rz)+bx/8/pi/(1-nu)*(eta/W/r+eta*C**2-eta*z2/Wr3+y*z/r3+&
    & 2*nu*sinA*C)-by*x/8/pi/(1-nu)*((1/Wr+C**2-z2/Wr3)*cosA+z/r3)+&
    & bz*x*sinA/8/pi/(1-nu)*(1/Wr+C**2-z2/Wr3)

    Exy = bx*(rFi_ry)/2+by*(rFi_rx)/2-bx/8/pi/(1-nu)*(x*y2/r2z2-nu*x/rz+&
    & x*y2/r3z-nu*x*cosA/Wr+eta*x*S/Wr+eta*x*y/Wr3)+ &
    & by/8/pi/(1-nu)*(x2*y/r2z2-nu*y/rz+x2*y/r3z+nu*cosA*S+ &
    & x2*y*cosA/Wr3+x2*cosA*S/Wr)-bz*sinA/8/pi/(1-nu)*(nu*S+x2*S/Wr+x2*y/Wr3)

    Exz = bx*(rFi_rz)/2+bz*(rFi_rx)/2-bx/8/pi/(1-nu)*(-x*y/r3+nu*x*sinA/Wr+eta*x*C/Wr+ &
    & eta*x*z/Wr3)+by/8/pi/(1-nu)*(-x2/r3+nu/r+nu*cosA*C+x2*z*cosA/Wr3+&
    & x2*cosA*C/Wr)-bz*sinA/8/pi/(1-nu)*(nu*C+x2*C/Wr+x2*z/Wr3)

    Eyz = by*(rFi_rz)/2+bz*(rFi_ry)/2+bx/8/pi/(1-nu)*(y2/r3-nu/r-nu*cosA*C+ &
    & nu*sinA*S+eta*sinA*cosA/W2-eta*(y*cosA+z*sinA)/W2r+eta*y*z/W2r2-eta*y*z/Wr3)- &
    & by*x/8/pi/(1-nu)*(y/r3+sinA*cosA**2/W2-cosA*(y*cosA+z*sinA)/&
    & W2r+y*z*cosA/W2r2-y*z*cosA/Wr3)-bz*x*sinA/8/pi/(1-nu)*&
    &(y*z/Wr3-sinA*cosA/W2+(y*cosA+z*sinA)/W2r-y*z/W2r2)

  end subroutine

  subroutine AngDisStrainFSC(y1,y2,y3,beta,b1,b2,b3,nu,a,&
    &v11,v22,v33,v12,v13,v23)
    implicit none
    real(8),intent(in)::y1,y2,y3,beta,b1,b2,b3,nu,a
    real(8),intent(out)::v11,v22,v33,v12,v13,v23
    real(8)::sinB,cosB,cotB,y3b,z1b,z3b,rb2,rb,W1,W2,W3,W4,W5,W6,W7,W8,W9,N1
    real(8)::rFib_ry2,rFib_ry1,rFib_ry3
    REAL*8, PARAMETER :: pi = 4*datan(1.d0)
  sinB = dsin(beta)
  cosB = dcos(beta)
  cotB = 1d0/dtan(beta)
  !write(*,*)'cotB',cotB
  y3b = y3+2*a
  z1b = y1*cosB+y3b*sinB
  z3b = -y1*sinB+y3b*cosB
  rb2 = y1*y1+y2*y2+y3b*y3b
  rb = dsqrt(rb2)
  !write(*,*)rb+z3b,rb+y3b!
  W1 = rb*cosB+y3b
  W2 = cosB+a/rb
  W3 = cosB+y3b/rb
  W4 = nu+a/rb
  W5 = 2*nu+a/rb
  W6 = rb+y3b
  W7 = rb+z3b
  W8 = y3+a
  W9 = 1+a/rb/cosB
  !
  N1 = 1-2*nu
  !write(*,*)
  ! % Partial derivatives of the Burgers' function
  rFib_ry2 = z1b/rb/(rb+z3b)-y1/rb/(rb+y3b)!; % y2 = x in ADCS
  rFib_ry1 = y2/rb/(rb+y3b)-cosB*y2/rb/(rb+z3b)!; % y1 = y in ADCS
  rFib_ry3 = -sinB*y2/rb/(rb+z3b)!; % y3 = z in ADCS
  !write(*,*)rFib_ry2,rFib_ry1,rFib_ry3

  v11 = b1*(0.25d0*((-2+2*nu)*N1*rFib_ry1*cotB**2-N1*y2/W6**2*((1-W5)*cotB-&
      &y1/W6*W4)/rb*y1+N1*y2/W6*(a/rb**3*y1*cotB-1/W6*W4+y1**2/&
      &W6**2*W4/rb+y1**2/W6*a/rb**3)-N1*y2*cosB*cotB/W7**2*W2*(y1/&
      &rb-sinB)-N1*y2*cosB*cotB/W7*a/rb**3*y1-3*a*y2*W8*cotB/rb**5*&
      &y1-y2*W8/rb**3/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1-y2*W8/&
      &rb2/W6**2*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1+y2*W8/rb/W6*&
      &(1/W6*W5-y1**2/W6**2*W5/rb-y1**2/W6*a/rb**3+a/rb2-2*a*y1**&
      &2/rb2**2)-y2*W8/rb**3/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+&
      &(2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*y1-y2*W8/rb/&
      &W7**2*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*&
      &cosB)-a*y3b*cosB*cotB/rb2)*(y1/rb-sinB)+y2*W8/rb/W7*(-cosB/&
      &W7**2*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)*(y1/&
      &rb-sinB)+cosB/W7*(1/rb*cosB*y1*(N1*cosB-a/rb)*cotB+W1*a/rb**&
      &3*y1*cotB+(2-2*nu)*(1/rb*sinB*y1-1)*cosB)+2*a*y3b*cosB*cotB/&
      &rb2**2*y1))/pi/(1-nu))+&
      &b2*(0.25d0*(N1*(((2-2*nu)*cotB**2+nu)/rb*y1/W6-((2-2*nu)*cotB**2+1)*&
      &cosB*(y1/rb-sinB)/W7)-N1/W6**2*(-N1*y1*cotB+nu*y3b-a+a*y1*&
      &cotB/rb+y1**2/W6*W4)/rb*y1+N1/W6*(-N1*cotB+a*cotB/rb-a*&
      &y1**2*cotB/rb**3+2*y1/W6*W4-y1**3/W6**2*W4/rb-y1**3/W6*a/&
      &rb**3)+N1*cotB/W7**2*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*(y1/&
      &rb-sinB)-N1*cotB/W7*(cosB**2-a*(1/rb*sinB*y1-1)/rb/cosB+a*&
      &(rb*sinB-y1)/rb**3/cosB*y1)-a*W8*cotB/rb**3+3*a*y1**2*W8*&
      &cotB/rb**5-W8/W6**2*(2*nu+1/rb*(N1*y1*cotB+a)-y1**2/rb/W6*&
      &W5-a*y1**2/rb**3)/rb*y1+W8/W6*(-1/rb**3*(N1*y1*cotB+a)*y1+&
      &1/rb*N1*cotB-2*y1/rb/W6*W5+y1**3/rb**3/W6*W5+y1**3/rb2/&
      &W6**2*W5+y1**3/rb2**2/W6*a-2*a/rb**3*y1+3*a*y1**3/rb**5)-W8*&
      &cotB/W7**2*(-cosB*sinB+a*y1*y3b/rb**3/cosB+(rb*sinB-y1)/rb*&
      &((2-2*nu)*cosB-W1/W7*W9))*(y1/rb-sinB)+W8*cotB/W7*(a*y3b/&
      &rb**3/cosB-3*a*y1**2*y3b/rb**5/cosB+(1/rb*sinB*y1-1)/rb*&
      &((2-2*nu)*cosB-W1/W7*W9)-(rb*sinB-y1)/rb**3*((2-2*nu)*cosB-W1/&
      &W7*W9)*y1+(rb*sinB-y1)/rb*(-1/rb*cosB*y1/W7*W9+W1/W7**2*&
      &W9*(y1/rb-sinB)+W1/W7*a/rb**3/cosB*y1)))/pi/(1-nu))+&
      &b3*(0.25d0*(N1*(-y2/W6**2*(1+a/rb)/rb*y1-y2/W6*a/rb**3*y1+y2*&
      &cosB/W7**2*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb**3*y1)+y2*W8/&
      &rb**3*(a/rb2+1/W6)*y1-y2*W8/rb*(-2*a/rb2**2*y1-1/W6**2/&
      &rb*y1)-y2*W8*cosB/rb**3/W7*(W1/W7*W2+a*y3b/rb2)*y1-y2*W8*&
      &cosB/rb/W7**2*(W1/W7*W2+a*y3b/rb2)*(y1/rb-sinB)+y2*W8*&
      &cosB/rb/W7*(1/rb*cosB*y1/W7*W2-W1/W7**2*W2*(y1/rb-sinB)-&
      &W1/W7*a/rb**3*y1-2*a*y3b/rb2**2*y1))/pi/(1-nu))

  v22 = b1*(0.25d0*(N1*(((2-2*nu)*cotB**2-nu)/rb*y2/W6-((2-2*nu)*cotB**2+1-&
      &2*nu)*cosB/rb*y2/W7)+N1/W6**2*(y1*cotB*(1-W5)+nu*y3b-a+y2**&
      &2/W6*W4)/rb*y2-N1/W6*(a*y1*cotB/rb**3*y2+2*y2/W6*W4-y2**&
      &3/W6**2*W4/rb-y2**3/W6*a/rb**3)+N1*z1b*cotB/W7**2*W2/rb*&
      &y2+N1*z1b*cotB/W7*a/rb**3*y2+3*a*y2*W8*cotB/rb**5*y1-W8/&
      &W6**2*(-2*nu+1/rb*(N1*y1*cotB-a)+y2**2/rb/W6*W5+a*y2**2/&
      &rb**3)/rb*y2+W8/W6*(-1/rb**3*(N1*y1*cotB-a)*y2+2*y2/rb/&
      &W6*W5-y2**3/rb**3/W6*W5-y2**3/rb2/W6**2*W5-y2**3/rb2**2/W6*&
      &a+2*a/rb**3*y2-3*a*y2**3/rb**5)-W8/W7**2*(cosB**2-1/rb*(N1*&
      &z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb**3-1/rb/W7*(y2**2*cosB**2-&
      &a*z1b*cotB/rb*W1))/rb*y2+W8/W7*(1/rb**3*(N1*z1b*cotB+a*&
      &cosB)*y2-3*a*y3b*z1b*cotB/rb**5*y2+1/rb**3/W7*(y2**2*cosB**2-&
      &a*z1b*cotB/rb*W1)*y2+1/rb2/W7**2*(y2**2*cosB**2-a*z1b*cotB/&
      &rb*W1)*y2-1/rb/W7*(2*y2*cosB**2+a*z1b*cotB/rb**3*W1*y2-a*&
      &z1b*cotB/rb2*cosB*y2)))/pi/(1-nu))+&
      &b2*(0.25d0*((2-2*nu)*N1*rFib_ry2*cotB**2+N1/W6*((W5-1)*cotB+y1/W6*&
      &W4)-N1*y2**2/W6**2*((W5-1)*cotB+y1/W6*W4)/rb+N1*y2/W6*(-a/&
      &rb**3*y2*cotB-y1/W6**2*W4/rb*y2-y2/W6*a/rb**3*y1)-N1*cotB/&
      &W7*W9+N1*y2**2*cotB/W7**2*W9/rb+N1*y2**2*cotB/W7*a/rb**3/&
      &cosB-a*W8*cotB/rb**3+3*a*y2**2*W8*cotB/rb**5+W8/rb/W6*(N1*&
      &cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))-y2**2*W8/rb**3/W6*&
      &(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))-y2**2*W8/rb2/W6**&
      &2*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))+y2*W8/rb/W6*&
      &(2*nu*y1/W6**2/rb*y2+a*y1/rb**3*(1/rb+1/W6)*y2-a*y1/rb*&
      &(-1/rb**3*y2-1/W6**2/rb*y2))+W8*cotB/rb/W7*((-2+2*nu)*cosB+&
      &W1/W7*W9+a*y3b/rb2/cosB)-y2**2*W8*cotB/rb**3/W7*((-2+2*nu)*&
      &cosB+W1/W7*W9+a*y3b/rb2/cosB)-y2**2*W8*cotB/rb2/W7**2*((-2+&
      &2*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)+y2*W8*cotB/rb/W7*(1/&
      &rb*cosB*y2/W7*W9-W1/W7**2*W9/rb*y2-W1/W7*a/rb**3/cosB*y2-&
      &2*a*y3b/rb2**2/cosB*y2))/pi/(1-nu))+&
      &b3*(0.25d0*(N1*(-sinB/rb*y2/W7+y2/W6**2*(1+a/rb)/rb*y1+y2/W6*&
      &a/rb**3*y1-z1b/W7**2*W2/rb*y2-z1b/W7*a/rb**3*y2)-y2*W8/&
      &rb**3*(a/rb2+1/W6)*y1+y1*W8/rb*(-2*a/rb2**2*y2-1/W6**2/&
      &rb*y2)+W8/W7**2*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/&
      &rb/W7*(y2**2*cosB*sinB-a*z1b/rb*W1))/rb*y2-W8/W7*(sinB*a/&
      &rb**3*y2-z1b/rb**3*(1+a*y3b/rb2)*y2-2*z1b/rb**5*a*y3b*y2+&
      &1/rb**3/W7*(y2**2*cosB*sinB-a*z1b/rb*W1)*y2+1/rb2/W7**2*&
      &(y2**2*cosB*sinB-a*z1b/rb*W1)*y2-1/rb/W7*(2*y2*cosB*sinB+a*&
      &z1b/rb**3*W1*y2-a*z1b/rb2*cosB*y2)))/pi/(1-nu))

  v33 = b1*(0.25d0*((2-2*nu)*(N1*rFib_ry3*cotB-y2/W6**2*W5*(y3b/rb+1)-&
      &0.5d0*y2/W6*a/rb**3*2*y3b+y2*cosB/W7**2*W2*W3+0.5d0*y2*cosB/W7*&
      &a/rb**3*2*y3b)+y2/rb*(2*nu/W6+a/rb2)-0.5d0*y2*W8/rb**3*(2*&
      &nu/W6+a/rb2)*2*y3b+y2*W8/rb*(-2*nu/W6**2*(y3b/rb+1)-a/&
      &rb2**2*2*y3b)+y2*cosB/rb/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-&
      &0.5d0*y2*W8*cosB/rb**3/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)*2*&
      &y3b-y2*W8*cosB/rb/W7**2*(1-2*nu-W1/W7*W2-a*y3b/rb2)*W3+y2*&
      &W8*cosB/rb/W7*(-(cosB*y3b/rb+1)/W7*W2+W1/W7**2*W2*W3+0.5d0*&
      &W1/W7*a/rb**3*2*y3b-a/rb2+a*y3b/rb2**2*2*y3b))/pi/(1-nu))+&
      &b2*(0.25d0*((-2+2*nu)*N1*cotB*((y3b/rb+1)/W6-cosB*W3/W7)+(2-2*nu)*&
      &y1/W6**2*W5*(y3b/rb+1)+0.5d0*(2-2*nu)*y1/W6*a/rb**3*2*y3b+(2-&
      &2*nu)*sinB/W7*W2-(2-2*nu)*z1b/W7**2*W2*W3-0.5d0*(2-2*nu)*z1b/&
      &W7*a/rb**3*2*y3b+1/rb*(N1*cotB-2*nu*y1/W6-a*y1/rb2)-0.5d0*&
      &W8/rb**3*(N1*cotB-2*nu*y1/W6-a*y1/rb2)*2*y3b+W8/rb*(2*nu*&
      &y1/W6**2*(y3b/rb+1)+a*y1/rb2**2*2*y3b)-1/W7*(cosB*sinB+W1*&
      &cotB/rb*((2-2*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*&
      &W1/rb/W7))+W8/W7**2*(cosB*sinB+W1*cotB/rb*((2-2*nu)*cosB-W1/&
      &W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*W3-W8/W7*((cosB*&
      &y3b/rb+1)*cotB/rb*((2-2*nu)*cosB-W1/W7)-0.5d0*W1*cotB/rb**3*&
      &((2-2*nu)*cosB-W1/W7)*2*y3b+W1*cotB/rb*(-(cosB*y3b/rb+1)/W7+&
      &W1/W7**2*W3)-0.5d0*a/rb**3*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7)*&
      &2*y3b+a/rb*(-z1b/rb2-y3b*sinB/rb2+y3b*z1b/rb2**2*2*y3b-&
      &sinB*W1/rb/W7-z1b*(cosB*y3b/rb+1)/rb/W7+0.5d0*z1b*W1/rb**3/&
      &W7*2*y3b+z1b*W1/rb/W7**2*W3)))/pi/(1-nu))+&
      &b3*(0.25d0*((2-2*nu)*rFib_ry3-(2-2*nu)*y2*sinB/W7**2*W2*W3-0.5d0*&
      &(2-2*nu)*y2*sinB/W7*a/rb**3*2*y3b+y2*sinB/rb/W7*(1+W1/W7*&
      &W2+a*y3b/rb2)-0.5d0*y2*W8*sinB/rb**3/W7*(1+W1/W7*W2+a*y3b/&
      &rb2)*2*y3b-y2*W8*sinB/rb/W7**2*(1+W1/W7*W2+a*y3b/rb2)*W3+&
      &y2*W8*sinB/rb/W7*((cosB*y3b/rb+1)/W7*W2-W1/W7**2*W2*W3-&
      &0.5d0*W1/W7*a/rb**3*2*y3b+a/rb2-a*y3b/rb2**2*2*y3b))/pi/(1-nu))

  v12 =b1/2*(0.25d0*((-2+2*nu)*N1*rFib_ry2*cotB**2+N1/W6*((1-W5)*cotB-y1/W6*W4)-&
      &N1*y2**2/W6**2*((1-W5)*cotB-y1/W6*W4)/rb+N1*y2/W6*&
      &(a/rb**3*y2*cotB+y1/W6**2*W4/rb*y2+y2/W6*a/rb**3*y1)+N1*&
      &cosB*cotB/W7*W2-N1*y2**2*cosB*cotB/W7**2*W2/rb-N1*y2**2*cosB*&
      &cotB/W7*a/rb**3+a*W8*cotB/rb**3-3*a*y2**2*W8*cotB/rb**5+W8/&
      &rb/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)-y2**2*W8/rb**3/W6*(-N1*&
      &cotB+y1/W6*W5+a*y1/rb2)-y2**2*W8/rb2/W6**2*(-N1*cotB+y1/&
      &W6*W5+a*y1/rb2)+y2*W8/rb/W6*(-y1/W6**2*W5/rb*y2-y2/W6*&
      &a/rb**3*y1-2*a*y1/rb2**2*y2)+W8/rb/W7*(cosB/W7*(W1*(N1*&
      &cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/&
      &rb2)-y2**2*W8/rb**3/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-&
      &2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)-y2**2*W8/rb2/&
      &W7**2*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*&
      &cosB)-a*y3b*cosB*cotB/rb2)+y2*W8/rb/W7*(-cosB/W7**2*(W1*&
      &(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)/rb*y2+cosB/&
      &W7*(1/rb*cosB*y2*(N1*cosB-a/rb)*cotB+W1*a/rb**3*y2*cotB+(2-2*&
      &nu)/rb*sinB*y2*cosB)+2*a*y3b*cosB*cotB/rb2**2*y2))/pi/(1-nu))+&
      &b2/2*(0.25d0*(N1*(((2-2*nu)*cotB**2+nu)/rb*y2/W6-((2-2*nu)*cotB**2+1)*&
      &cosB/rb*y2/W7)-N1/W6**2*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/rb+&
      &y1**2/W6*W4)/rb*y2+N1/W6*(-a*y1*cotB/rb**3*y2-y1**2/W6**&
      &2*W4/rb*y2-y1**2/W6*a/rb**3*y2)+N1*cotB/W7**2*(z1b*cosB-a*&
      &(rb*sinB-y1)/rb/cosB)/rb*y2-N1*cotB/W7*(-a/rb2*sinB*y2/&
      &cosB+a*(rb*sinB-y1)/rb**3/cosB*y2)+3*a*y2*W8*cotB/rb**5*y1-&
      &W8/W6**2*(2*nu+1/rb*(N1*y1*cotB+a)-y1**2/rb/W6*W5-a*y1**2/&
      &rb**3)/rb*y2+W8/W6*(-1/rb**3*(N1*y1*cotB+a)*y2+y1**2/rb**&
      &3/W6*W5*y2+y1**2/rb2/W6**2*W5*y2+y1**2/rb2**2/W6*a*y2+3*&
      &a*y1**2/rb**5*y2)-W8*cotB/W7**2*(-cosB*sinB+a*y1*y3b/rb**3/&
      &cosB+(rb*sinB-y1)/rb*((2-2*nu)*cosB-W1/W7*W9))/rb*y2+W8*cotB/&
      &W7*(-3*a*y1*y3b/rb**5/cosB*y2+1/rb2*sinB*y2*((2-2*nu)*cosB-&
      &W1/W7*W9)-(rb*sinB-y1)/rb**3*((2-2*nu)*cosB-W1/W7*W9)*y2+(rb*&
      &sinB-y1)/rb*(-1/rb*cosB*y2/W7*W9+W1/W7**2*W9/rb*y2+W1/W7*&
      &a/rb**3/cosB*y2)))/pi/(1-nu))+&
      &b3/2*(0.25d0*(N1*(1/W6*(1+a/rb)-y2**2/W6**2*(1+a/rb)/rb-y2**2/&
      &W6*a/rb**3-cosB/W7*W2+y2**2*cosB/W7**2*W2/rb+y2**2*cosB/W7*&
      &a/rb**3)-W8/rb*(a/rb2+1/W6)+y2**2*W8/rb**3*(a/rb2+1/W6)-&
      &y2*W8/rb*(-2*a/rb2**2*y2-1/W6**2/rb*y2)+W8*cosB/rb/W7*&
      &(W1/W7*W2+a*y3b/rb2)-y2**2*W8*cosB/rb**3/W7*(W1/W7*W2+a*&
      &y3b/rb2)-y2**2*W8*cosB/rb2/W7**2*(W1/W7*W2+a*y3b/rb2)+y2*&
      &W8*cosB/rb/W7*(1/rb*cosB*y2/W7*W2-W1/W7**2*W2/rb*y2-W1/&
      &W7*a/rb**3*y2-2*a*y3b/rb2**2*y2))/pi/(1-nu))+&
      &b1/2*(0.25d0*(N1*(((2-2*nu)*cotB**2-nu)/rb*y1/W6-((2-2*nu)*cotB**2+1-&
      &2*nu)*cosB*(y1/rb-sinB)/W7)+N1/W6**2*(y1*cotB*(1-W5)+nu*y3b-&
      &a+y2**2/W6*W4)/rb*y1-N1/W6*((1-W5)*cotB+a*y1**2*cotB/rb**3-&
      &y2**2/W6**2*W4/rb*y1-y2**2/W6*a/rb**3*y1)-N1*cosB*cotB/W7*&
      &W2+N1*z1b*cotB/W7**2*W2*(y1/rb-sinB)+N1*z1b*cotB/W7*a/rb**&
      &3*y1-a*W8*cotB/rb**3+3*a*y1**2*W8*cotB/rb**5-W8/W6**2*(-2*&
      &nu+1/rb*(N1*y1*cotB-a)+y2**2/rb/W6*W5+a*y2**2/rb**3)/rb*&
      &y1+W8/W6*(-1/rb**3*(N1*y1*cotB-a)*y1+1/rb*N1*cotB-y2**2/&
      &rb**3/W6*W5*y1-y2**2/rb2/W6**2*W5*y1-y2**2/rb2**2/W6*a*y1-&
      &3*a*y2**2/rb**5*y1)-W8/W7**2*(cosB**2-1/rb*(N1*z1b*cotB+a*&
      &cosB)+a*y3b*z1b*cotB/rb**3-1/rb/W7*(y2**2*cosB**2-a*z1b*cotB/&
      &rb*W1))*(y1/rb-sinB)+W8/W7*(1/rb**3*(N1*z1b*cotB+a*cosB)*&
      &y1-1/rb*N1*cosB*cotB+a*y3b*cosB*cotB/rb**3-3*a*y3b*z1b*cotB/&
      &rb**5*y1+1/rb**3/W7*(y2**2*cosB**2-a*z1b*cotB/rb*W1)*y1+1/&
      &rb/W7**2*(y2**2*cosB**2-a*z1b*cotB/rb*W1)*(y1/rb-sinB)-1/rb/&
      &W7*(-a*cosB*cotB/rb*W1+a*z1b*cotB/rb**3*W1*y1-a*z1b*cotB/&
      &rb2*cosB*y1)))/pi/(1-nu))+&
      &b2/2*(0.25d0*((2-2*nu)*N1*rFib_ry1*cotB**2-N1*y2/W6**2*((W5-1)*cotB+&
      &y1/W6*W4)/rb*y1+N1*y2/W6*(-a/rb**3*y1*cotB+1/W6*W4-y1**&
      &2/W6**2*W4/rb-y1**2/W6*a/rb**3)+N1*y2*cotB/W7**2*W9*(y1/&
      &rb-sinB)+N1*y2*cotB/W7*a/rb**3/cosB*y1+3*a*y2*W8*cotB/rb**&
      &5*y1-y2*W8/rb**3/W6*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/&
      &W6))*y1-y2*W8/rb2/W6**2*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/&
      &rb+1/W6))*y1+y2*W8/rb/W6*(-2*nu/W6+2*nu*y1**2/W6**2/rb-a/&
      &rb*(1/rb+1/W6)+a*y1**2/rb**3*(1/rb+1/W6)-a*y1/rb*(-1/&
      &rb**3*y1-1/W6**2/rb*y1))-y2*W8*cotB/rb**3/W7*((-2+2*nu)*&
      &cosB+W1/W7*W9+a*y3b/rb2/cosB)*y1-y2*W8*cotB/rb/W7**2*((-2+&
      &2*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)*(y1/rb-sinB)+y2*W8*&
      &cotB/rb/W7*(1/rb*cosB*y1/W7*W9-W1/W7**2*W9*(y1/rb-sinB)-&
      &W1/W7*a/rb**3/cosB*y1-2*a*y3b/rb2**2/cosB*y1))/pi/(1-nu))+&
      &b3/2*(0.25d0*(N1*(-sinB*(y1/rb-sinB)/W7-1/W6*(1+a/rb)+y1**2/W6**&
      &2*(1+a/rb)/rb+y1**2/W6*a/rb**3+cosB/W7*W2-z1b/W7**2*W2*&
      &(y1/rb-sinB)-z1b/W7*a/rb**3*y1)+W8/rb*(a/rb2+1/W6)-y1**2*&
      &W8/rb**3*(a/rb2+1/W6)+y1*W8/rb*(-2*a/rb2**2*y1-1/W6**2/&
      &rb*y1)+W8/W7**2*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/&
      &rb/W7*(y2**2*cosB*sinB-a*z1b/rb*W1))*(y1/rb-sinB)-W8/W7*&
      &(sinB*a/rb**3*y1+cosB/rb*(1+a*y3b/rb2)-z1b/rb**3*(1+a*y3b/&
      &rb2)*y1-2*z1b/rb**5*a*y3b*y1+1/rb**3/W7*(y2**2*cosB*sinB-a*&
      &z1b/rb*W1)*y1+1/rb/W7**2*(y2**2*cosB*sinB-a*z1b/rb*W1)*&
      &(y1/rb-sinB)-1/rb/W7*(-a*cosB/rb*W1+a*z1b/rb**3*W1*y1-a*&
      &z1b/rb2*cosB*y1)))/pi/(1-nu))

  v13 = b1/2*(0.25d0*((-2+2*nu)*N1*rFib_ry3*cotB**2-N1*y2/W6**2*((1-W5)*&
      &cotB-y1/W6*W4)*(y3b/rb+1)+N1*y2/W6*(0.5d0*a/rb**3*2*y3b*cotB+&
      &y1/W6**2*W4*(y3b/rb+1)+0.5d0*y1/W6*a/rb**3*2*y3b)-N1*y2*cosB*&
      &cotB/W7**2*W2*W3-0.5d0*N1*y2*cosB*cotB/W7*a/rb**3*2*y3b+a/&
      &rb**3*y2*cotB-1.5d0*a*y2*W8*cotB/rb**5*2*y3b+y2/rb/W6*(-N1*&
      &cotB+y1/W6*W5+a*y1/rb2)-0.5d0*y2*W8/rb**3/W6*(-N1*cotB+y1/&
      &W6*W5+a*y1/rb2)*2*y3b-y2*W8/rb/W6**2*(-N1*cotB+y1/W6*W5+&
      &a*y1/rb2)*(y3b/rb+1)+y2*W8/rb/W6*(-y1/W6**2*W5*(y3b/rb+&
      &1)-0.5d0*y1/W6*a/rb**3*2*y3b-a*y1/rb2**2*2*y3b)+y2/rb/W7*&
      &(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)-&
      &a*y3b*cosB*cotB/rb2)-0.5d0*y2*W8/rb**3/W7*(cosB/W7*(W1*(N1*&
      &cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/&
      &rb2)*2*y3b-y2*W8/rb/W7**2*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+&
      &(2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*W3+y2*W8/rb/&
      &W7*(-cosB/W7**2*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*&
      &cosB)*W3+cosB/W7*((cosB*y3b/rb+1)*(N1*cosB-a/rb)*cotB+0.5d0*W1*&
      &a/rb**3*2*y3b*cotB+0.5d0*(2-2*nu)/rb*sinB*2*y3b*cosB)-a*cosB*&
      &cotB/rb2+a*y3b*cosB*cotB/rb2**2*2*y3b))/pi/(1-nu))+&
      &b2/2*(0.25d0*(N1*(((2-2*nu)*cotB**2+nu)*(y3b/rb+1)/W6-((2-2*nu)*cotB**&
      &2+1)*cosB*W3/W7)-N1/W6**2*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/&
      &rb+y1**2/W6*W4)*(y3b/rb+1)+N1/W6*(nu-0.5d0*a*y1*cotB/rb**3*2*&
      &y3b-y1**2/W6**2*W4*(y3b/rb+1)-0.5d0*y1**2/W6*a/rb**3*2*y3b)+&
      &N1*cotB/W7**2*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*W3-N1*cotB/&
      &W7*(cosB*sinB-0.5d0*a/rb2*sinB*2*y3b/cosB+0.5d0*a*(rb*sinB-y1)/&
      &rb**3/cosB*2*y3b)-a/rb**3*y1*cotB+1.5d0*a*y1*W8*cotB/rb**5*2*&
      &y3b+1/W6*(2*nu+1/rb*(N1*y1*cotB+a)-y1**2/rb/W6*W5-a*y1**2/&
      &rb**3)-W8/W6**2*(2*nu+1/rb*(N1*y1*cotB+a)-y1**2/rb/W6*W5-a*&
      &y1**2/rb**3)*(y3b/rb+1)+W8/W6*(-0.5d0/rb**3*(N1*y1*cotB+a)*2*&
      &y3b+0.5d0*y1**2/rb**3/W6*W5*2*y3b+y1**2/rb/W6**2*W5*(y3b/rb+&
      &1)+0.5d0*y1**2/rb2**2/W6*a*2*y3b+1.5d0*a*y1**2/rb**5*2*y3b)+&
      &cotB/W7*(-cosB*sinB+a*y1*y3b/rb**3/cosB+(rb*sinB-y1)/rb*((2-&
      &2*nu)*cosB-W1/W7*W9))-W8*cotB/W7**2*(-cosB*sinB+a*y1*y3b/rb**&
      &3/cosB+(rb*sinB-y1)/rb*((2-2*nu)*cosB-W1/W7*W9))*W3+W8*cotB/&
      &W7*(a/rb**3/cosB*y1-1.5d0*a*y1*y3b/rb**5/cosB*2*y3b+0.5d0/&
      &rb2*sinB*2*y3b*((2-2*nu)*cosB-W1/W7*W9)-0.5d0*(rb*sinB-y1)/rb**&
      &3*((2-2*nu)*cosB-W1/W7*W9)*2*y3b+(rb*sinB-y1)/rb*(-(cosB*y3b/&
      &rb+1)/W7*W9+W1/W7**2*W9*W3+0.5d0*W1/W7*a/rb**3/cosB*2*&
      &y3b)))/pi/(1-nu))+&
      &b3/2*(0.25d0*(N1*(-y2/W6**2*(1+a/rb)*(y3b/rb+1)-0.5d0*y2/W6*a/&
      &rb**3*2*y3b+y2*cosB/W7**2*W2*W3+0.5d0*y2*cosB/W7*a/rb**3*2*&
      &y3b)-y2/rb*(a/rb2+1/W6)+0.5d0*y2*W8/rb**3*(a/rb2+1/W6)*2*&
      &y3b-y2*W8/rb*(-a/rb2**2*2*y3b-1/W6**2*(y3b/rb+1))+y2*cosB/&
      &rb/W7*(W1/W7*W2+a*y3b/rb2)-0.5d0*y2*W8*cosB/rb**3/W7*(W1/&
      &W7*W2+a*y3b/rb2)*2*y3b-y2*W8*cosB/rb/W7**2*(W1/W7*W2+a*&
      &y3b/rb2)*W3+y2*W8*cosB/rb/W7*((cosB*y3b/rb+1)/W7*W2-W1/&
      &W7**2*W2*W3-0.5d0*W1/W7*a/rb**3*2*y3b+a/rb2-a*y3b/rb2**2*2*&
      &y3b))/pi/(1-nu))+&
      &b1/2*(0.25d0*((2-2*nu)*(N1*rFib_ry1*cotB-y1/W6**2*W5/rb*y2-y2/W6*&
      &a/rb**3*y1+y2*cosB/W7**2*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb**&
      &3*y1)-y2*W8/rb**3*(2*nu/W6+a/rb2)*y1+y2*W8/rb*(-2*nu/W6**&
      &2/rb*y1-2*a/rb2**2*y1)-y2*W8*cosB/rb**3/W7*(1-2*nu-W1/W7*&
      &W2-a*y3b/rb2)*y1-y2*W8*cosB/rb/W7**2*(1-2*nu-W1/W7*W2-a*&
      &y3b/rb2)*(y1/rb-sinB)+y2*W8*cosB/rb/W7*(-1/rb*cosB*y1/W7*&
      &W2+W1/W7**2*W2*(y1/rb-sinB)+W1/W7*a/rb**3*y1+2*a*y3b/rb2**&
      &2*y1))/pi/(1-nu))+&
      &b2/2*(0.25d0*((-2+2*nu)*N1*cotB*(1/rb*y1/W6-cosB*(y1/rb-sinB)/W7)-&
      &(2-2*nu)/W6*W5+(2-2*nu)*y1**2/W6**2*W5/rb+(2-2*nu)*y1**2/W6*&
      &a/rb**3+(2-2*nu)*cosB/W7*W2-(2-2*nu)*z1b/W7**2*W2*(y1/rb-&
      &sinB)-(2-2*nu)*z1b/W7*a/rb**3*y1-W8/rb**3*(N1*cotB-2*nu*y1/&
      &W6-a*y1/rb2)*y1+W8/rb*(-2*nu/W6+2*nu*y1**2/W6**2/rb-a/rb2+&
      &2*a*y1**2/rb2**2)+W8/W7**2*(cosB*sinB+W1*cotB/rb*((2-2*nu)*&
      &cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*(y1/rb-&
      &sinB)-W8/W7*(1/rb2*cosB*y1*cotB*((2-2*nu)*cosB-W1/W7)-W1*&
      &cotB/rb**3*((2-2*nu)*cosB-W1/W7)*y1+W1*cotB/rb*(-1/rb*cosB*&
      &y1/W7+W1/W7**2*(y1/rb-sinB))-a/rb**3*(sinB-y3b*z1b/rb2-&
      &z1b*W1/rb/W7)*y1+a/rb*(-y3b*cosB/rb2+2*y3b*z1b/rb2**2*y1-&
      &cosB*W1/rb/W7-z1b/rb2*cosB*y1/W7+z1b*W1/rb**3/W7*y1+z1b*&
      &W1/rb/W7**2*(y1/rb-sinB))))/pi/(1-nu))+&
      &b3/2*(0.25d0*((2-2*nu)*rFib_ry1-(2-2*nu)*y2*sinB/W7**2*W2*(y1/rb-&
      &sinB)-(2-2*nu)*y2*sinB/W7*a/rb**3*y1-y2*W8*sinB/rb**3/W7*(1+&
      &W1/W7*W2+a*y3b/rb2)*y1-y2*W8*sinB/rb/W7**2*(1+W1/W7*W2+&
      &a*y3b/rb2)*(y1/rb-sinB)+y2*W8*sinB/rb/W7*(1/rb*cosB*y1/&
      &W7*W2-W1/W7**2*W2*(y1/rb-sinB)-W1/W7*a/rb**3*y1-2*a*y3b/&
      &rb2**2*y1))/pi/(1-nu))

  v23 = b1/2*(0.25d0*(N1*(((2-2*nu)*cotB**2-nu)*(y3b/rb+1)/W6-((2-2*nu)*&
      &cotB**2+1-2*nu)*cosB*W3/W7)+N1/W6**2*(y1*cotB*(1-W5)+nu*y3b-a+&
      &y2**2/W6*W4)*(y3b/rb+1)-N1/W6*(0.5d0*a*y1*cotB/rb**3*2*y3b+&
      &nu-y2**2/W6**2*W4*(y3b/rb+1)-0.5d0*y2**2/W6*a/rb**3*2*y3b)-N1*&
      &sinB*cotB/W7*W2+N1*z1b*cotB/W7**2*W2*W3+0.5d0*N1*z1b*cotB/W7*&
      &a/rb**3*2*y3b-a/rb**3*y1*cotB+1.5d0*a*y1*W8*cotB/rb**5*2*y3b+&
      &1/W6*(-2*nu+1/rb*(N1*y1*cotB-a)+y2**2/rb/W6*W5+a*y2**2/&
      &rb**3)-W8/W6**2*(-2*nu+1/rb*(N1*y1*cotB-a)+y2**2/rb/W6*W5+&
      &a*y2**2/rb**3)*(y3b/rb+1)+W8/W6*(-0.5d0/rb**3*(N1*y1*cotB-a)*&
      &2*y3b-0.5d0*y2**2/rb**3/W6*W5*2*y3b-y2**2/rb/W6**2*W5*(y3b/&
      &rb+1)-0.5d0*y2**2/rb2**2/W6*a*2*y3b-1.5d0*a*y2**2/rb**5*2*y3b)+&
      &1/W7*(cosB**2-1/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb**&
      &3-1/rb/W7*(y2**2*cosB**2-a*z1b*cotB/rb*W1))-W8/W7**2*(cosB**2-&
      &1/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb**3-1/rb/W7*&
      &(y2**2*cosB**2-a*z1b*cotB/rb*W1))*W3+W8/W7*(0.5d0/rb**3*(N1*&
      &z1b*cotB+a*cosB)*2*y3b-1/rb*N1*sinB*cotB+a*z1b*cotB/rb**3+a*&
      &y3b*sinB*cotB/rb**3-1.5d0*a*y3b*z1b*cotB/rb**5*2*y3b+0.5d0/rb**&
      &3/W7*(y2**2*cosB**2-a*z1b*cotB/rb*W1)*2*y3b+1/rb/W7**2*(y2**&
      &2*cosB**2-a*z1b*cotB/rb*W1)*W3-1/rb/W7*(-a*sinB*cotB/rb*W1+&
      &0.5d0*a*z1b*cotB/rb**3*W1*2*y3b-a*z1b*cotB/rb*(cosB*y3b/rb+&
      &1))))/pi/(1-nu))+&
      &b2/2*(0.25d0*((2-2*nu)*N1*rFib_ry3*cotB**2-N1*y2/W6**2*((W5-1)*cotB+&
      &y1/W6*W4)*(y3b/rb+1)+N1*y2/W6*(-0.5d0*a/rb**3*2*y3b*cotB-y1/&
      &W6**2*W4*(y3b/rb+1)-0.5d0*y1/W6*a/rb**3*2*y3b)+N1*y2*cotB/&
      &W7**2*W9*W3+0.5d0*N1*y2*cotB/W7*a/rb**3/cosB*2*y3b-a/rb**3*&
      &y2*cotB+1.5d0*a*y2*W8*cotB/rb**5*2*y3b+y2/rb/W6*(N1*cotB-2*&
      &nu*y1/W6-a*y1/rb*(1/rb+1/W6))-0.5d0*y2*W8/rb**3/W6*(N1*&
      &cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))*2*y3b-y2*W8/rb/W6**&
      &2*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))*(y3b/rb+1)+y2*&
      &W8/rb/W6*(2*nu*y1/W6**2*(y3b/rb+1)+0.5d0*a*y1/rb**3*(1/rb+&
      &1/W6)*2*y3b-a*y1/rb*(-0.5d0/rb**3*2*y3b-1/W6**2*(y3b/rb+&
      &1)))+y2*cotB/rb/W7*((-2+2*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)-&
      &0.5d0*y2*W8*cotB/rb**3/W7*((-2+2*nu)*cosB+W1/W7*W9+a*y3b/&
      &rb2/cosB)*2*y3b-y2*W8*cotB/rb/W7**2*((-2+2*nu)*cosB+W1/W7*&
      &W9+a*y3b/rb2/cosB)*W3+y2*W8*cotB/rb/W7*((cosB*y3b/rb+1)/&
      &W7*W9-W1/W7**2*W9*W3-0.5d0*W1/W7*a/rb**3/cosB*2*y3b+a/rb2/&
      &cosB-a*y3b/rb2**2/cosB*2*y3b))/pi/(1-nu))+&
      &b3/2*(0.25d0*(N1*(-sinB*W3/W7+y1/W6**2*(1+a/rb)*(y3b/rb+1)+&
      &0.5d0*y1/W6*a/rb**3*2*y3b+sinB/W7*W2-z1b/W7**2*W2*W3-0.5d0*&
      &z1b/W7*a/rb**3*2*y3b)+y1/rb*(a/rb2+1/W6)-0.5d0*y1*W8/rb**&
      &3*(a/rb2+1/W6)*2*y3b+y1*W8/rb*(-a/rb2**2*2*y3b-1/W6**2*&
      &(y3b/rb+1))-1/W7*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/&
      &rb/W7*(y2**2*cosB*sinB-a*z1b/rb*W1))+W8/W7**2*(sinB*(cosB-&
      &a/rb)+z1b/rb*(1+a*y3b/rb2)-1/rb/W7*(y2**2*cosB*sinB-a*z1b/&
      &rb*W1))*W3-W8/W7*(0.5d0*sinB*a/rb**3*2*y3b+sinB/rb*(1+a*y3b/&
      &rb2)-0.5d0*z1b/rb**3*(1+a*y3b/rb2)*2*y3b+z1b/rb*(a/rb2-a*&
      &y3b/rb2**2*2*y3b)+0.5d0/rb**3/W7*(y2**2*cosB*sinB-a*z1b/rb*&
      &W1)*2*y3b+1/rb/W7**2*(y2**2*cosB*sinB-a*z1b/rb*W1)*W3-1/&
      &rb/W7*(-a*sinB/rb*W1+0.5d0*a*z1b/rb**3*W1*2*y3b-a*z1b/rb*&
      &(cosB*y3b/rb+1))))/pi/(1-nu))+&
      &b1/2*(0.25d0*((2-2*nu)*(N1*rFib_ry2*cotB+1/W6*W5-y2**2/W6**2*W5/&
      &rb-y2**2/W6*a/rb**3-cosB/W7*W2+y2**2*cosB/W7**2*W2/rb+y2**2*&
      &cosB/W7*a/rb**3)+W8/rb*(2*nu/W6+a/rb2)-y2**2*W8/rb**3*(2*&
      &nu/W6+a/rb2)+y2*W8/rb*(-2*nu/W6**2/rb*y2-2*a/rb2**2*y2)+&
      &W8*cosB/rb/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-y2**2*W8*cosB/&
      &rb**3/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-y2**2*W8*cosB/rb2/W7**&
      &2*(1-2*nu-W1/W7*W2-a*y3b/rb2)+y2*W8*cosB/rb/W7*(-1/rb*&
      &cosB*y2/W7*W2+W1/W7**2*W2/rb*y2+W1/W7*a/rb**3*y2+2*a*&
      &y3b/rb2**2*y2))/pi/(1-nu))+&
      &b2/2*(0.25d0*((-2+2*nu)*N1*cotB*(1/rb*y2/W6-cosB/rb*y2/W7)+(2-&
      &2*nu)*y1/W6**2*W5/rb*y2+(2-2*nu)*y1/W6*a/rb**3*y2-(2-2*&
      &nu)*z1b/W7**2*W2/rb*y2-(2-2*nu)*z1b/W7*a/rb**3*y2-W8/rb**&
      &3*(N1*cotB-2*nu*y1/W6-a*y1/rb2)*y2+W8/rb*(2*nu*y1/W6**2/&
      &rb*y2+2*a*y1/rb2**2*y2)+W8/W7**2*(cosB*sinB+W1*cotB/rb*((2-&
      &2*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))/&
      &rb*y2-W8/W7*(1/rb2*cosB*y2*cotB*((2-2*nu)*cosB-W1/W7)-W1*&
      &cotB/rb**3*((2-2*nu)*cosB-W1/W7)*y2+W1*cotB/rb*(-cosB/rb*&
      &y2/W7+W1/W7**2/rb*y2)-a/rb**3*(sinB-y3b*z1b/rb2-z1b*W1/&
      &rb/W7)*y2+a/rb*(2*y3b*z1b/rb2**2*y2-z1b/rb2*cosB*y2/W7+&
      &z1b*W1/rb**3/W7*y2+z1b*W1/rb2/W7**2*y2)))/pi/(1-nu))+&
      &b3/2*(0.25d0*((2-2*nu)*rFib_ry2+(2-2*nu)*sinB/W7*W2-(2-2*nu)*y2**2*&
      &sinB/W7**2*W2/rb-(2-2*nu)*y2**2*sinB/W7*a/rb**3+W8*sinB/rb/&
      &W7*(1+W1/W7*W2+a*y3b/rb2)-y2**2*W8*sinB/rb**3/W7*(1+W1/&
      &W7*W2+a*y3b/rb2)-y2**2*W8*sinB/rb2/W7**2*(1+W1/W7*W2+a*&
      &y3b/rb2)+y2*W8*sinB/rb/W7*(1/rb*cosB*y2/W7*W2-W1/W7**2*&
      &W2/rb*y2-W1/W7*a/rb**3*y2-2*a*y3b/rb2**2*y2))/pi/(1-nu))
end subroutine
  subroutine DCross(v1,v2,v3)
    implicit none
    real(8),dimension(3)::v1,v2,v3
    v3(1)=v1(2)*v2(3)-v2(2)*v1(3)
    v3(2)=v1(3)*v2(1)-v2(3)*v1(1)
    v3(3)=v1(1)*v2(2)-v2(1)*v1(2)

  end subroutine
end module dtriangular
