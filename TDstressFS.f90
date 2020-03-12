module dtriangular
  USE, INTRINSIC :: IEEE_ARITHMETIC
  implicit none
  !real(8)::Sxx,Syy,Szz,Sxy,Sxz,Syz
  !REAL*8, DIMENSION(3):: P1, P2, P3
  !P1=(/-1.d0,0.d0,0.d0/)
  !P2=(/1.d0,-1.d0,-1.d0/)
  !P3=(/0.d0,15.d0,5.d0/)
  !call  TDstressFS(0.d0,0.d0,0.d0,P1,P2,P3,-1.d0,2.d0,3.d0,1d0,1d0, &
  !& Sxx,Syy,Szz,Sxy,Sxz,Syz)
  !write(*,'(6e12.4)') Sxx,Syy,Szz,Sxy,Sxz,Syz
  !stop
contains
  !-------------------------------------------------------------------------------
  subroutine TDstressHS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    implicit none
    real(8),intent(in)::X,Y,Z,Ss,Ds,Ts,mu,lambda
    real(8),intent(inout)::P1(3),P2(3),P3(3)
    real(8),intent(out)::Sxx,Syy,Szz,Sxy,Sxz,Syz
    real(8)::Sxx_m,Syy_m,Szz_m,Sxy_m,Sxz_m,Syz_m
    real(8)::Sxx_h,Syy_h,Szz_h,Sxy_h,Sxz_h,Syz_h
    real(8)::Sxx_i,Syy_i,Szz_i,Sxy_i,Sxz_i,Syz_i

    ! if any(Z>0 | P1(3)>0 | P2(3)>0 | P3(3)>0)
    !     error('Half-space solution: Z coordinates must be negative!')
    ! elseif P1(3)==0 && P2(3)==0 && P3(3)==0
    !     Stress = zeros(numel(X),6);
    !     Strain = zeros(numel(X),6);
    !     return
    ! end
    ! % Calculate main dislocation contribution to strains and stresses
    ! [StsMS,StrMS] = TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);
    call TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    &Sxx_m,Syy_m,Szz_m,Sxy_m,Sxz_m,Syz_m)
    !write(*,*)Sxx_m,Syy_m,Szz_m,Sxy_m,Sxz_m,Syz_m
    !
    ! % Calculate harmonic function contribution to strains and stresses
    ! [StsFSC,StrFSC] = TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);
    call TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    &Sxx_h,Syy_h,Szz_h,Sxy_h,Sxz_h,Syz_h)
    !write(*,*)'Sxx_h,Syy_h,Szz_h,Sxy_h,Sxz_h,Syz_h'
    !write(*,*)Sxx_h,Syy_h,Szz_h,Sxy_h,Sxz_h,Syz_h
    ! % Calculate image dislocation contribution to strains and stresses
    P1(3) = -P1(3)
    P2(3) = -P2(3)
    P3(3) = -P3(3)
    ! [StsIS,StrIS] = TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);
    call TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    &Sxx_i,Syy_i,Szz_i,Sxy_i,Sxz_i,Syz_i)
    !write(*,*)Sxx_i,Syy_i,Szz_i,Sxy_i,Sxz_i,Syz_i
    ! % Calculate the complete stress and strain tensor components in EFCS
    Sxx=Sxx_m+Sxx_h+Sxx_i
    Syy=Syy_m+Syy_h+Syy_i
    Szz=Szz_m+Szz_h+Szz_i
    Sxy=Sxy_m+Sxy_h+Sxy_i
    Sxz=Sxz_m+Sxz_h+Sxz_i
    Syz=Syz_m+Syz_h+Syz_i
    ! Stress = StsMS+StsIS+StsFSC;
    ! Strain = StrMS+StrIS+StrFSC;

  end subroutine
  subroutine TDstressFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda,&
    & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: X, Y, Z ! <==== In my version, these are scalars,
    !       and (X, Y, Z) is the single observation point.
    !       However, in the MatLab version, these are (potentially very long) column vectors,
    !       with values rearranged from 0-D, 1-D, 2-D, or 3-D arrays!
    REAL*8, DIMENSION(3), INTENT(IN) :: P1, P2, P3 ! The three vertices of the triangular dislocation.
    REAL*8, INTENT(IN) :: Ss, Ds, Ts ! 3 components of Burger's vector: Strike-slip, Dip-slip, Tensile-slip
    REAL*8 :: nu,mu,lambda
    REAL*8, INTENT(OUT) :: Sxx,Syy,Szz,Sxy,Sxz,Syz ! 3 components of displacement at the observation point (East, North, Vertical/Up).
    REAL*8 :: exx,eyy,ezz,exy,exz,eyz
    REAL*8 :: Exx_,Eyy_,Ezz_,Exy_,Exz_,Eyz_

    INTEGER :: Trimode
    LOGICAL :: casepLog, casenLog, casezLog
    REAL*8 :: A, B, C, bx, by, bz, Fi, na, nb, nc, &
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
    ! Calculates stresses and strains associated with a triangular dislocation
    ! in an elastic full-space
    !
    ! TD: Triangular Dislocation
    ! EFCS: Earth-Fixed Coordinate System
    ! TDCS: Triangular Dislocation Coordinate System
    ! ADCS: Angular Dislocation Coordinate System
    !
    ! INPUTS
    ! X, Y and Z:
    ! Coordinates of calculation points in EFCS (East, North, Up) X, Y and Z
    ! must have the same size
    !
    ! P1,P2 and P3:
    ! Coordinates of TD vertices in EFCS
    !
    ! Ss, Ds and Ts:
    ! TD slip vector components (Strike-slip, Dip-slip, Tensile-slip)
    !
    ! mu and lambda:
    ! Lame constants
    !
    ! OUTPUTS
    ! Stress:
    ! Calculated stress tensor components in EFCS The six columns of Stress
    ! are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively The stress components
    ! have the same unit as Lame constants
    !
    ! Strain:
    ! Calculated strain tensor components in EFCS The six columns of Strain
    ! are Exx, Eyy, Ezz, Exy, Exz and Eyz, respectively The strain components
    ! are dimensionless

    nu = 1/(1+lambda/mu)/2 ! Poisson's ratio
    !!!write(*,*) 'calc'

    bx = Ts ! Tensile-slip
    by = Ss ! Strike-slip
    bz = Ds ! Dip-slip
    !!write(*,*) 'bx,by,bz',bx,by,bz


    ! Calculate unit strike, dip and normal to TD vectors: For a horizontal TD
    ! as an exception, if the normal vector points upward, the strike and dip
    ! vectors point Northward and Westward, whereas if the normal vector points
    ! downward, the strike and dip vectors point Southward and Westward,
    ! respectively

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
    CALL DCross(Vnorm, Vstrike, Vdip) ! (apparently, no normalization is needed here)
    !!write(*,*) Vnorm,Vstrike,Vdip
    p_1 = 0.0D0 ! "_" = "lower-case"
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

    ! Determine the best arteact-free configuration for each calculation point
    CALL trimodefinder(y_, z_, x_, p_1(2:3), p_2(2:3), p_3(2:3), & ! inputs
    & Trimode)
    !!!write(*,*)Trimode                                 ! output
    casepLog = (Trimode == 1) ! Note that Fortran TRUE and FALSE are not necessarily represented the same way as MatLab 1 and 0 !
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
    !if (casezLog) then
    !  exx = IEEE_VALUE(u, IEEE_QUIET_NAN)
    !  eyy = IEEE_VALUE(u, IEEE_QUIET_NAN)
    !  ezz = IEEE_VALUE(u, IEEE_QUIET_NAN)
    !  exy = IEEE_VALUE(u, IEEE_QUIET_NAN)
    !  exz = IEEE_VALUE(u, IEEE_QUIET_NAN)
    !  eyz = IEEE_VALUE(u, IEEE_QUIET_NAN)
    !end if

    Sxx = 2*mu*exx+lambda*(exx+eyy+ezz)
    Syy = 2*mu*eyy+lambda*(exx+eyy+ezz)
    Szz = 2*mu*ezz+lambda*(exx+eyy+ezz)
    Sxy = 2*mu*exy
    Sxz = 2*mu*exz
    Syz = 2*mu*eyz
    !!write(*,*)'Sxx,Syy,Szz,Sxy,Sxz,Syz in TDCS'
    !!write(*,*)Sxx,Syy,Szz,Sxy,Sxz,Syz


    ! Transform the strain tensor components from TDCS into EFCS

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
  ! function [Stress,Strain]=TDstress_HarFunc(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda)
  ! % TDstress_HarFunc calculates the harmonic function contribution to the
  ! % strains and stresses associated with a triangular dislocation in a
  ! % half-space. The function cancels the surface normal tractions induced by
  ! % the main and image dislocations.
  !
  bx = Ts!; % Tensile-slip
  by = Ss!; % Strike-slip
  bz = Ds!; % Dip-slip
  !
  ! % Calculate unit strike, dip and normal to TD vectors: For a horizontal TD
  ! % as an exception, if the normal vector points upward, the strike and dip
  ! % vectors point Northward and Westward, whereas if the normal vector points
  ! % downward, the strike and dip vectors point Southward and Westward,
  ! % respectively.
  ! Vnorm = cross(P2-P1,P3-P1);
  ! Vnorm = Vnorm/norm(Vnorm);
  !
  ! eY = [0 1 0]';
  ! eZ = [0 0 1]';
  ! Vstrike = cross(eZ,Vnorm);
  !
  ! if norm(Vstrike)==0
  !     Vstrike = eY*Vnorm(3);
  ! end
  ! Vstrike = Vstrike/norm(Vstrike);
  ! Vdip = cross(Vnorm,Vstrike);
  P2mP1 = P2 - P1 ! all 3 components
  P3mP1 = P3 - P1
  !write(*,*) P2mP1,P3mP1
  CALL DCross(P2mP1, P3mP1, Vnorm) ! but, this still needs to be normalized:
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
  CALL DCross(Vnorm, Vstrike, Vdip) ! (apparently, no normalization is needed here)
  !write(*,*)Vdip
  ! % Transform slip vector components from TDCS into EFCS
  ! A = [Vnorm Vstrike Vdip];
  ! [bX,bY,bZ] = CoordTrans(bx,by,bz,A);
  At(1:3, 1) = Vnorm(1:3)
  At(1:3, 2) = Vstrike(1:3)
  At(1:3, 3) = Vdip(1:3)
  CALL CoordTrans(bx, by, bz, At, & ! inputs
  & bxt, byt, bzt)
  !write(*,*) bx,by,bz,bxt,byt,bzt
  !
  ! % Calculate contribution of angular dislocation pair on each TD side
  ! [Stress1,Strain1] = AngSetupFSC_S(X,Y,Z,bX,bY,bZ,P1,P2,mu,lambda); % P1P2
  ! [Stress2,Strain2] = AngSetupFSC_S(X,Y,Z,bX,bY,bZ,P2,P3,mu,lambda); % P2P3
  ! [Stress3,Strain3] = AngSetupFSC_S(X,Y,Z,bX,bY,bZ,P3,P1,mu,lambda); % P3P1
  call AngSetupFSC_S(X,Y,Z,bXt,bYt,bZt,P1,P2,mu,lambda,&
  &Sxx1,Syy1,Szz1,Sxy1,Sxz1,Syz1)
  !write(*,*)Sxx1,Syy1,Szz1,Sxy1,Sxz1,Syz1
  !write(*,*)
  call AngSetupFSC_S(X,Y,Z,bXt,bYt,bZt,P2,P3,mu,lambda,&
  &Sxx2,Syy2,Szz2,Sxy2,Sxz2,Syz2)
  !write(*,*)Sxx2,Syy2,Szz2,Sxy2,Sxz2,Syz2
  !write(*,*)
  call AngSetupFSC_S(X,Y,Z,bXt,bYt,bZt,P3,P1,mu,lambda,&
  &Sxx3,Syy3,Szz3,Sxy3,Sxz3,Syz3)
  !write(*,*)Sxx3,Syy3,Szz3,Sxy3,Sxz3,Syz3
  !write(*,*)
  !
  ! % Calculate total harmonic function contribution to strains and stresses
  ! Stress = Stress1+Stress2+Stress3;
  ! Strain = Strain1+Strain2+Strain3;
  Sxx=Sxx1+Sxx2+Sxx3
  Syy=Syy1+Syy2+Syy3
  Szz=Szz1+Szz2+Szz3
  Sxy=Sxy1+Sxy2+Sxy3
  Sxz=Sxz1+Sxz2+Sxz3
  Syz=Syz1+Syz2+Syz3
  return
  end subroutine
  !------------------------------------------------------------------------------
  subroutine TensTrans(Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1,A,&
    & Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2)
    implicit none
    REAL*8, INTENT(IN) :: Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1
    REAL*8, DIMENSION(3, 3) :: A
    REAL*8, INTENT(OUT) :: Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2
    ! TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate
    ! system to x2y2z2 coordinate system "A" is the transformation matrix,
    ! whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1 The
    ! coordinates of e1,e2 and e3 in A must be given in x2y2z2 The transpose
    ! of A (ie, A') does the transformation from x2y2z2 into x1y1z1

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
  SUBROUTINE CoordTrans(x1, x2, x3, A, & ! inputs
    & Xcap1, Xcap2, Xcap3)      ! outputs
    ! CoordTrans transforms the coordinates of the vectors, from
    ! x1x2x3 coordinate system to X1X2X3 coordinate system "A" is the
    ! transformation matrix, whose columns e1,e2 and e3 are the unit base
    ! vectors of the x1x2x3 The coordinates of e1,e2 and e3 in A must be given
    ! in X1X2X3 The transpose of A (ie, A') will transform the coordinates
    ! from X1X2X3 into x1x2x3
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1, x2, x3
    REAL*8, DIMENSION(3, 3) :: A
    REAL*8, INTENT(OUT) :: Xcap1, Xcap2, Xcap3
    !  "cap" = "capital" added because Fortran does not distinguish upper- from lower-case

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
    ! trimodefinder calculates the normalized barycentric coordinates of
    ! the points with respect to the TD vertices and specifies the appropriate
    ! artefact-free configuration of the angular dislocations for the
    ! calculations The input matrices x, y and z share the same size and
    ! correspond to the y, z and x coordinates in the TDCS, respectively p1,
    ! p2 and p3 are two-component matrices representing the y and z coordinates
    ! of the TD vertices in the TDCS, respectively
    ! The components of the output (trimode) corresponding to each calculation
    ! points, are 1 for the first configuration, -1 for the second
    ! configuration and 0 for the calculation point that lie on the TD sides


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
    ! TDSetupS transforms coordinates of the calculation points as well as
    ! slip vector components from ADCS into TDCS It then calculates the
    ! strains in ADCS and transforms them into TDCS
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: alpha, bx, by, bz, nu, x, y, z
    REAL*8, DIMENSION(3), INTENT(IN) :: SideVec, TriVertex
    REAL*8, INTENT(OUT) :: exx,eyy,ezz,exy,exz,eyz
    REAL*8 :: exx_,eyy_,ezz_,exy_,exz_,eyz_
    real(8),parameter::pi=3.14159265358979d0

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

  ! % AngSetupFSC_S calculates the Free Surface Correction to strains and
  ! % stresses associated with angular dislocation pair on each TD side.
  !
  nu = lambda/(lambda+mu)/2!; % Poisson's ratio
  !
  ! % Calculate TD side vector and the angle of the angular dislocation pair
  SideVec = PB-PA
  eZ = (/ 0.0D0, 0.0D0, 1.0D0 /)
  beta = acos(-SideVec(3)/sqrt(dot_product(SideVec,SideVec)))!;
  !
  if ((abs(beta)<1d-8).or.(abs(pi-beta)<1d-8)) then
  !     Stress = zeros(length(X),6);
  Sxx=0d0;Syy=0d0;Szz=0d0;Sxy=0d0;Sxz=0d0;Syz=0d0
  !     Strain = zeros(length(X),6);
  else
       ey1 = (/SideVec(1),SideVec(2),0d0/)
       ey1 = ey1/sqrt(dot_product(ey1,ey1))
       ey3 = -eZ
       call Dcross(ey3,ey1,ey2)
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
  !write(*,*)'y1A,y2A,y3A',y1A,y2A,y3A
  !     % Transform coordinates from EFCS to the second ADCS
  CALL CoordTrans(SideVec(1),SideVec(2),SideVec(3), A, & ! inputs
  & y1AB,y2AB,y3AB)
  !     [y1AB,y2AB,y3AB] = CoordTrans(SideVec(1),SideVec(2),SideVec(3),A);
       y1B = y1A-y1AB
       y2B = y2A-y2AB
       y3B = y3A-y3AB
  !
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
  if (beta*y1A.ge.0) configI=.true.
  !
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
       !write(*,*)v11,v22,v33,v12,v13,v23
  !
  !     % Transform total Free Surface Correction to strains from ADCS to EFCS
  !     [Exx,Eyy,Ezz,Exy,Exz,Eyz] = TensTrans(v11,v22,v33,v12,v13,v23,A');
  call TensTrans(v11,v22,v33,v12,v13,v23,A,&
    & Exx,Eyy,Ezz,Exy,Exz,Eyz)
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
    ! AngDisStrain calculates the strains associated with an angular
    ! dislocation in an elastic full-space
    !Note that the orginal MatLab version allows x, y, and z to be arrays of test points
    !(in 0-D, 1-D, 2-D, or 3-D) however, in this Fortran version there is only a single
    !test point at (x, y, z), and each of these is a simple REAL*8 scalar number

    IMPLICIT NONE
    REAL*8, PARAMETER :: pi = 3.14159265358979D0
    REAL*8, INTENT(IN) :: x, y, z, alpha, bx, by, bz, nu
    REAL*8, INTENT(OUT) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
    REAL*8 :: cosA, eta, r, sinA, ux, uy, uz, vx, vy, vz, wx, wy, wz, zz, zeta
    Real*8 :: C,r2,r2z2,r3,r3z,rFi_rx,rFi_ry,rFi_rz,rz,S,W,W2,W2r,W2r2,Wr,Wr3,x2,y2,z2

    sinA = sin(alpha)
    cosA = cos(alpha)
    eta = y*cosA-z*sinA
    zeta = y*sinA+z*cosA

    x2 = x**2
    y2 = y**2
    z2 = z**2
    r2 = x2+y2+z2
    r = sqrt(r2)
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
  !-------------------------------------------------------------------------------
  ! function [v11,v22,v33,v12,v13,v23] = AngDisStrainFSC(y1,y2,y3,beta,...
  !     b1,b2,b3,nu,a)
  subroutine AngDisStrainFSC(y1,y2,y3,beta,b1,b2,b3,nu,a,&
    &v11,v22,v33,v12,v13,v23)
    implicit none
    real(8),intent(in)::y1,y2,y3,beta,b1,b2,b3,nu,a
    real(8),intent(out)::v11,v22,v33,v12,v13,v23
    real(8)::sinB,cosB,cotB,y3b,z1b,z3b,rb2,rb,W1,W2,W3,W4,W5,W6,W7,W8,W9,N1
    real(8)::rFib_ry2,rFib_ry1,rFib_ry3
    REAL*8, PARAMETER :: pi = 3.14159265358979D0
  ! % AngDisStrainFSC calculates the harmonic function contribution to the
  ! % strains associated with an angular dislocation in an elastic half-space.
  !
  sinB = sin(beta)
  cosB = cos(beta)
  cotB = 1d0/tan(beta)
  !write(*,*)'cotB',cotB
  y3b = y3+2*a
  z1b = y1*cosB+y3b*sinB
  z3b = -y1*sinB+y3b*cosB
  rb2 = y1*y1+y2*y2+y3b*y3b
  rb = sqrt(rb2)
  !
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
  !
  ! % Partial derivatives of the Burgers' function
  rFib_ry2 = z1b/rb/(rb+z3b)-y1/rb/(rb+y3b)!; % y2 = x in ADCS
  rFib_ry1 = y2/rb/(rb+y3b)-cosB*y2/rb/(rb+z3b)!; % y1 = y in ADCS
  rFib_ry3 = -sinB*y2/rb/(rb+z3b)!; % y3 = z in ADCS
  !
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
