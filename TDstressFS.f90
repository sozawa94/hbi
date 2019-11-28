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
    !write(*,*) 'calc'

    bx = Ts ! Tensile-slip
    by = Ss ! Strike-slip
    bz = Ds ! Dip-slip


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
    !write(*,*) Vnorm(1),Vstrike(2),Vdip(3)
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
    !write(*,*) At
    CALL CoordTrans(Ts, Ss, Ds, At,&
    & bx,by,bz)
  !write(*,*) bx,by,bz

    ! Determine the best arteact-free configuration for each calculation point
    CALL trimodefinder(y_, z_, x_, p_1(2:3), p_2(2:3), p_3(2:3), & ! inputs
    & Trimode)
    !write(*,*)Trimode                                 ! output
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
    !IF (casezLog) THEN
    !  xn = x_
    !  yn = y_
    !  zn = z_
    !END IF

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

    ! IF (casezLog) THEN
    !  CALL TDSetupD(xn, yn, zn, A, bx, by, bz, nu, p_1, e13, & ! inputs
    !  & Exx1Tn,Eyy1Tn,Ezz1Tn,Exy1Tn,Exz1Tn,Eyz1Tn)
    !  CALL TDSetupD(xn, yn, zn, B, bx, by, bz, nu, p_2, -e12, & ! inputs
    !  & Exx2Tn,Eyy2Tn,Ezz2Tn,Exy2Tn,Exz2Tn,Eyz2Tn)
    !  CALL TDSetupD(xn, yn, zn, C, bx, by, bz, nu, p_3, -e23, & ! inputs
    !  & Exx3Tn,Eyy3Tn,Ezz3Tn,Exy3Tn,Exz3Tn,Eyz3Tn)
    !END IF

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
    !if (casezLog) then
    !  exx = Exx1Tn+Exx2Tn+Exx3Tn
    !  eyy = Eyy1Tn+Eyy2Tn+Eyy3Tn
    !  ezz = Ezz1Tn+Ezz2Tn+Ezz3Tn
    !  exy = Exy1Tn+Exy2Tn+Exy3Tn
    !  exz = Exz1Tn+Exz2Tn+Exz3Tn
    !  eyz = Eyz1Tn+Eyz2Tn+Eyz3Tn
    !end if
    if (casezLog) then
      exx = IEEE_VALUE(u, IEEE_QUIET_NAN)
      eyy = IEEE_VALUE(u, IEEE_QUIET_NAN)
      ezz = IEEE_VALUE(u, IEEE_QUIET_NAN)
      exy = IEEE_VALUE(u, IEEE_QUIET_NAN)
      exz = IEEE_VALUE(u, IEEE_QUIET_NAN)
      eyz = IEEE_VALUE(u, IEEE_QUIET_NAN)
    end if

    Sxx = 2*mu*exx+lambda*(exx+eyy+ezz)
    Syy = 2*mu*eyy+lambda*(exx+eyy+ezz)
    Szz = 2*mu*ezz+lambda*(exx+eyy+ezz)
    Sxy = 2*mu*exy
    Sxz = 2*mu*exz
    Syz = 2*mu*eyz

    ! Transform the strain tensor components from TDCS into EFCS

    !call TensTrans(exx,eyy,ezz,exy,exz,eyz,At,&
    !& Exx_,Eyy_,Ezz_,Exy_,Exz_,Eyz_)

    ! Calculate the stress tensor components in EFCS
    !Sxx = 2*mu*Exx_+lambda*(Exx_+Eyy_+Ezz_)
    !Syy = 2*mu*Eyy_+lambda*(Exx_+Eyy_+Ezz_)
    !Szz = 2*mu*Ezz_+lambda*(Exx_+Eyy_+Ezz_)
    !Sxy = 2*mu*Exy_
    !Sxz = 2*mu*Exz_
    !Syz = 2*mu*Eyz_

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
  subroutine DCross(v1,v2,v3)
    implicit none
    real(8),dimension(3)::v1,v2,v3
    v3(1)=v1(2)*v2(3)-v2(2)*v1(3)
    v3(2)=v1(3)*v2(1)-v2(3)*v1(1)
    v3(3)=v1(1)*v2(2)-v2(1)*v1(2)

  end subroutine
end module dtriangular
