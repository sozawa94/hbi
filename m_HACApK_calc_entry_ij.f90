module m_HACApK_calc_entry_ij
  !use m_matel_ij
  use mod_constant
  use dtriangular
  use mod_okada
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  ! Immutable after initialization: safe for concurrent ACA entry evaluations.
  type :: layered_antiplane_model
    real(8), allocatable :: bottoms(:), moduli(:), wave(:), weight(:)
    real(8), allocatable :: densities(:)
  end type

  type :: st_HACApK_calc_entry
    real*8,pointer :: ao(:)
    character(128)::problem
    integer :: nd,lp61
    real(8),pointer::xcol(:),ycol(:),zcol(:),xel(:),xer(:),yel(:),yer(:),ang(:),angd(:),rake(:)
    real(8),pointer::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
    real(8),pointer::ys1(:),ys2(:),ys3(:),ys4(:)
    real(8),pointer::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:)
    real(8),pointer::ds(:),dsl(:),dsd(:)
    real(8)::w,fwid
    !real(8),pointer::ds
    character(128)::v,md
    type(layered_antiplane_model) :: layered
  end type st_HACApK_calc_entry

  public :: HACApK_entry_ij
contains
  subroutine geological_frame(input_normal,normal,strike,updip)
    ! z is upward. Rake is measured from RHR strike towards updip.
    ! An unoriented plane has a unique upward normal except at vertical dip.
    ! Vertical tie: choose ny<0, then nx>0; horizontal tie: strike is +x.
    implicit none
    real(8),intent(in)::input_normal(3)
    real(8),intent(out)::normal(3),strike(3),updip(3)
    real(8),parameter::tol=1d-12
    real(8)::length,horizontal
    length=sqrt(dot_product(input_normal,input_normal))
    if(length<=0d0) error stop 'Degenerate fault normal'
    normal=input_normal/length
    if(normal(3)<-tol) then
      normal=-normal
    else if(abs(normal(3))<=tol) then
      if(normal(2)>tol.or.(abs(normal(2))<=tol.and.normal(1)<0d0)) normal=-normal
    end if
    horizontal=sqrt(normal(1)**2+normal(2)**2)
    if(horizontal>tol) then
      strike=(/-normal(2),normal(1),0d0/)/horizontal
    else
      strike=(/1d0,0d0,0d0/)
    end if
    call DCross(normal,strike,updip)
    updip=updip/sqrt(dot_product(updip,updip))
  end subroutine geological_frame

  subroutine geological_triangle_slip(p1,p2,p3,rake,ss,ds)
    ! Convert geological slip to the native, vertex-order-dependent TD basis.
    ! The native displacement jump changes side when its normal is reversed.
    implicit none
    real(8),intent(in)::p1(3),p2(3),p3(3),rake
    real(8),intent(out)::ss,ds
    real(8)::raw_normal(3),native_strike(3),native_dip(3),normal(3),strike(3),updip(3),jump(3),length,side
    call DCross(p2-p1,p3-p1,raw_normal)
    length=sqrt(dot_product(raw_normal,raw_normal))
    if(length<=0d0) error stop 'Degenerate triangular fault element'
    raw_normal=raw_normal/length
    call geological_frame(raw_normal,normal,strike,updip)
    call DCross((/0d0,0d0,1d0/),raw_normal,native_strike)
    length=sqrt(dot_product(native_strike,native_strike))
    if(length>0d0) then
      native_strike=native_strike/length
    else
      ! Same horizontal-plane fallback as the TD kernels.
      native_strike=(/0d0,raw_normal(3),0d0/)
    end if
    call DCross(raw_normal,native_strike,native_dip)
    side=sign(1d0,dot_product(raw_normal,normal))
    jump=side*(strike*cos(rake)+updip*sin(rake))
    ss=dot_product(jump,native_strike)
    ds=dot_product(jump,native_dip)
  end subroutine geological_triangle_slip

  real(8) function geological_triangle_entry(i,j,st_bemv,halfspace) result(value)
    implicit none
    integer,intent(in)::i,j
    logical,intent(in)::halfspace
    type(st_HACApK_calc_entry),intent(in)::st_bemv
    real(8)::p1(3),p2(3),p3(3),ss,ds,ts,source_rake,lambda
    real(8)::sxx,syy,szz,sxy,sxz,syz,normal(3),strike(3),updip(3),traction(3),slip_direction(3)
    p1=(/st_bemv%xs1(j),st_bemv%ys1(j),st_bemv%zs1(j)/)
    p2=(/st_bemv%xs2(j),st_bemv%ys2(j),st_bemv%zs2(j)/)
    p3=(/st_bemv%xs3(j),st_bemv%ys3(j),st_bemv%zs3(j)/)
    ts=0d0
    select case(trim(st_bemv%md))
    case('st')
      source_rake=0d0
    case('dp')
      source_rake=pi/2
    case('s','o')
      source_rake=st_bemv%rake(j)
    case('open')
      ts=1d0;source_rake=0d0
    case default
      error stop 'Unsupported triangular slip mode'
    end select
    call geological_triangle_slip(p1,p2,p3,source_rake,ss,ds)
    if(ts/=0d0) then
      ss=0d0;ds=0d0
    end if
    lambda=2*pois*rigid/(1-2*pois)
    if(halfspace) then
      call TDstressHS(st_bemv%xcol(i),st_bemv%ycol(i),st_bemv%zcol(i),p1,p2,p3,ss,ds,ts,rigid,lambda, &
        sxx,syy,szz,sxy,sxz,syz)
    else
      call TDstressFS(st_bemv%xcol(i),st_bemv%ycol(i),st_bemv%zcol(i),p1,p2,p3,ss,ds,ts,rigid,lambda, &
        sxx,syy,szz,sxy,sxz,syz)
    end if
    call geological_frame((/st_bemv%ev31(i),st_bemv%ev32(i),st_bemv%ev33(i)/),normal,strike,updip)
    traction=(/sxx*normal(1)+sxy*normal(2)+sxz*normal(3), &
      sxy*normal(1)+syy*normal(2)+syz*normal(3),sxz*normal(1)+syz*normal(2)+szz*normal(3)/)
    select case(trim(st_bemv%v))
    case('n')
      value=-dot_product(normal,traction)
    case('s','o')
      slip_direction=strike*cos(st_bemv%rake(i))+updip*sin(st_bemv%rake(i))
      value=dot_product(slip_direction,traction)
    case default
      error stop 'Unsupported triangular traction component'
    end select
  end function geological_triangle_entry

  !***HACApK_entry_ij
  real*8 function HACApK_entry_ij(i, j, st_bemv)
    ! zbemv is data needed for calculating the value of the element
    implicit none
    integer::i,j
    real(8)::tmp1,tmp2
    logical::fullspace
    type(st_HACApK_calc_entry) :: st_bemv
    select case(st_bemv%problem)
    case('2dp')
      HACApK_entry_ij=matels2dp_ij(i,j,st_bemv%xcol,st_bemv%xel,st_bemv%xer,st_bemv%v,st_bemv%md)

    case('2dpfz')
      HACApK_entry_ij=matels2dpFZ_ij(i,j,st_bemv%xcol,st_bemv%xel,st_bemv%xer,st_bemv%w)

    case('2dn_vector')
      HACApK_entry_ij=matel2dn_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang,st_bemv%v)

    case('2dn')
      HACApK_entry_ij=tensor2d_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang,st_bemv%v,st_bemv%md,st_bemv%rake)

    case('2dph','2dnh')
      HACApK_entry_ij=matel2dph2_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang,st_bemv%v)
      !HACApK_entry_ij=matel2dph_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%ang,st_bemv%v)

    case('2dna')
      HACApK_entry_ij=tensor2d3_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      & st_bemv%xel,st_bemv%xer,st_bemv%yel,st_bemv%yer,st_bemv%ang)

    case('2dvs')
      if(allocated(st_bemv%layered%moduli)) then
        select case(trim(st_bemv%v))
        case('s')
          if(trim(st_bemv%md)/='s') error stop 'Layered 2dvs supports strike slip only'
          HACApK_entry_ij=matel2dvs_layered_ij(i,j,st_bemv)
        case('n')
          HACApK_entry_ij=0d0 ! Antiplane slip produces no normal traction.
        case default
          error stop 'Unsupported layered 2dvs traction component'
        end select
      else
        HACApK_entry_ij=matels2dpa_ij(i,j,st_bemv%xcol,st_bemv%xel,st_bemv%xer)-&
        & matels2dpa_ij(i,j,st_bemv%xcol,-st_bemv%xel,-st_bemv%xer)
      end if
      !   &tensor2d3_ij(i,j,st_bemv%xcol,st_bemv%ycol,&
      !   & st_bemv%xel,st_bemv%xer,-st_bemv%yel,-st_bemv%yer,st_bemv%ang)

    case('3dp')
      HACApK_entry_ij=matel3dp_ij(i,j,st_bemv%xcol,st_bemv%zcol,&
      & st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%xs4,&
      & st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%zs4,st_bemv%v, st_bemv%md)
    ! case('3dph')
    !   HACApK_entry_ij=matel3dph_ij(i,j,st_bemv%xcol,st_bemv%zcol,&
    !   & st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%xs4,&
    !   & st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%zs4)
    !   !write(*,*)HACApK_entry_ij

    case('3dhr','3dph')
      fullspace=.false.
      HACApK_entry_ij=okada_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%zcol,st_bemv%ang,st_bemv%angd,st_bemv%dsl,st_bemv%dsd,st_bemv%v,st_bemv%rake,st_bemv%w,fullspace)
      !write(*,*)HACApK_entry_ij

    case('3dnr')
      fullspace=.true.
      HACApK_entry_ij=okada_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%zcol,st_bemv%ang,st_bemv%angd,st_bemv%dsl,st_bemv%dsd,st_bemv%v,st_bemv%rake,st_bemv%w,fullspace)
      !write(*,*)HACApK_entry_ij
    case('3dn','3dnt')
      !HACApK_entry_ij=matels1_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%zcol, st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%ys1,st_bemv%ys2,st_bemv%ys3,st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%v,st_bemv%md)
      HACApK_entry_ij=matel3dn_ij(i,j,st_bemv)
    case('3dh','3dht') !half space
      HACApK_entry_ij=matel3dh_ij(i,j,st_bemv)
      !HACApK_entry_ij=matelh1_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%zcol,st_bemv%strike,st_bemv%dip,st_bemv%v,st_bemv%md)
    case('25d')
      !HACApK_entry_ij=matels1_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%zcol, st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%ys1,st_bemv%ys2,st_bemv%ys3,st_bemv%zs1,st_bemv%zs2,st_bemv%zs3,st_bemv%v,st_bemv%md)
      HACApK_entry_ij=matelrec_ij(i,j,st_bemv%xcol,st_bemv%ycol,st_bemv%ds,st_bemv%ang,st_bemv%fwid,st_bemv%v)
    end select
    !stop
    return
  end function HACApK_entry_ij

  subroutine init_layered_2dvs(st_bemv,bottoms,moduli,wave_order,cutoff,densities)
    ! depths in km, mu in GPa: returned GPa/km equals HBI's MPa/m.
    ! xcol/xel/xer are POSITIVE DOWNWARD depths for problem='2dvs'.
    ! Initialize on every MPI rank, before HACApK_generate (never inside ACA).
    implicit none
    type(st_HACApK_calc_entry), intent(inout) :: st_bemv
    real(8), intent(in) :: bottoms(:),moduli(:)
    integer, optional, intent(in) :: wave_order
    real(8), optional, intent(in) :: cutoff
    real(8), optional, intent(in) :: densities(:) ! kg/m^3; optional for static-only API use
    integer :: n,nq,l,p,panels,offset
    real(8) :: d,scale,qmax,lo,hi,cut,tol
    real(8), allocatable :: nodes(:),weights(:)
    n=size(bottoms); nq=32; cut=30d0
    if(present(wave_order)) nq=wave_order
    if(present(cutoff)) cut=cutoff
    if(trim(st_bemv%problem)/='2dvs') error stop 'Layered elasticity requires problem=2dvs'
    if(size(moduli)/=n+1) error stop 'Need n layer bottoms and n+1 moduli'
    if(any(.not.ieee_is_finite(bottoms)).or.any(.not.ieee_is_finite(moduli))) &
      error stop 'Nonfinite layered material'
    if(any(bottoms<=0d0).or.any(moduli<=0d0)) error stop 'Layer depths and moduli must be positive'
    if(present(densities)) then
      if(size(densities)/=n+1) error stop 'Need n+1 layer densities'
      if(any(.not.ieee_is_finite(densities)).or.any(densities<=0d0)) &
        error stop 'Layer densities must be finite and positive'
    end if
    do l=2,n
      if(bottoms(l)<=bottoms(l-1)) error stop 'Layer bottoms must increase'
    end do
    if(nq<2.or.nq>512.or..not.ieee_is_finite(cut).or.cut<=0d0) error stop 'Invalid wave quadrature'
    if(.not.associated(st_bemv%xcol).or..not.associated(st_bemv%xel).or..not.associated(st_bemv%xer)) &
      error stop 'Initialize 2dvs geometry before layered elasticity'
    if(size(st_bemv%xcol)==0.or.size(st_bemv%xel)/=size(st_bemv%xcol).or. &
      size(st_bemv%xer)/=size(st_bemv%xcol)) error stop 'Invalid layered 2dvs geometry sizes'
    if(any(.not.ieee_is_finite(st_bemv%xcol)).or.any(.not.ieee_is_finite(st_bemv%xel)).or. &
      any(.not.ieee_is_finite(st_bemv%xer))) error stop 'Nonfinite 2dvs geometry'
    if(any(st_bemv%xel<0d0).or.any(st_bemv%xer<=st_bemv%xel).or. &
      any(st_bemv%xcol<=st_bemv%xel).or.any(st_bemv%xcol>=st_bemv%xer)) &
      error stop 'Require positive-depth elements and interior collocation points'
    do l=2,size(st_bemv%xcol)
      if(st_bemv%xel(l)<st_bemv%xer(l-1)) error stop 'Layered 2dvs elements must be ordered and nonoverlapping'
    end do
    scale=maxval(st_bemv%xer)
    if(n>0) scale=max(scale,maxval(bottoms))
    tol=64d0*epsilon(1d0)*scale
    ! For ordered nonoverlapping intervals the closest tip is an own-element tip.
    d=min(minval(st_bemv%xcol-st_bemv%xel),minval(st_bemv%xer-st_bemv%xcol))
    do l=1,size(st_bemv%xcol)
      if(n>0) d=min(d,minval(abs(st_bemv%xcol(l)-bottoms)))
    end do
    if(d<=tol) error stop 'Collocation on a source tip or material interface is unsupported'
    if(allocated(st_bemv%layered%moduli)) then
      deallocate(st_bemv%layered%bottoms,st_bemv%layered%moduli,st_bemv%layered%wave,st_bemv%layered%weight)
    end if
    if(allocated(st_bemv%layered%densities)) deallocate(st_bemv%layered%densities)
    ! Fujitsu's default flags do not enable automatic allocatable assignment.
    allocate(st_bemv%layered%bottoms(n),st_bemv%layered%moduli(n+1))
    st_bemv%layered%bottoms(:)=bottoms; st_bemv%layered%moduli(:)=moduli
    if(present(densities)) then
      allocate(st_bemv%layered%densities(n+1))
      st_bemv%layered%densities(:)=densities
    end if
    ! Geometric panels resolve both long wavelengths and small tip/interface distances.
    qmax=cut/d; lo=0d0; hi=min(1d0/scale,qmax); panels=1
    do while(hi<qmax)
      hi=min(2d0*hi,qmax); panels=panels+1
    end do
    allocate(nodes(nq),weights(nq))
    call layered_gauss_legendre(nodes,weights)
    allocate(st_bemv%layered%wave(nq*panels),st_bemv%layered%weight(nq*panels))
    hi=min(1d0/scale,qmax)
    do p=1,panels
      offset=(p-1)*nq
      st_bemv%layered%wave(offset+1:offset+nq)=lo+(nodes+1d0)*(hi-lo)/2d0
      st_bemv%layered%weight(offset+1:offset+nq)=weights*(hi-lo)/(2d0*pi)
      lo=hi; hi=min(2d0*hi,qmax)
    end do
  end subroutine

  subroutine read_layered_2dvs(st_bemv,filename,wave_order,cutoff)
    ! File: n; n rows (bottom[km], mu[GPa], rho[kg/m^3]); half-space (mu,rho).
    implicit none
    type(st_HACApK_calc_entry), intent(inout) :: st_bemv
    character(*), intent(in) :: filename
    integer, optional, intent(in) :: wave_order
    real(8), optional, intent(in) :: cutoff
    integer :: unit,ios,n,l
    real(8), allocatable :: bottoms(:),moduli(:),densities(:)
    character(512) :: line
    open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0) error stop 'Cannot open elastic_layer_file'
    read(unit,*,iostat=ios) n
    if(ios/=0) error stop 'Cannot read elastic layer count'
    if(n<0.or.n>10000) error stop 'Invalid elastic layer count'
    allocate(bottoms(n),moduli(n+1),densities(n+1))
    do l=1,n
      read(unit,'(a)',iostat=ios) line
      if(ios/=0) error stop 'Cannot read elastic layer row'
      read(line,*,iostat=ios) bottoms(l),moduli(l),densities(l)
      if(ios/=0) error stop 'Layer row requires bottom[km], mu[GPa], rho[kg/m^3]'
    end do
    read(unit,'(a)',iostat=ios) line
    if(ios/=0) error stop 'Cannot read half-space row'
    read(line,*,iostat=ios) moduli(n+1),densities(n+1)
    if(ios/=0) error stop 'Half-space row requires mu[GPa], rho[kg/m^3]'
    close(unit)
    call init_layered_2dvs(st_bemv,bottoms,moduli,wave_order,cutoff,densities)
  end subroutine

  subroutine layered_material_at(model,z,mu,rho,vs_local,eta)
    ! mu [GPa], rho [kg/m^3], Vs [km/s], eta [MPa*s/m].
    ! sqrt(mu*1e9/rho)/1000 = sqrt(mu*1000/rho).
    ! GPa/(km/s) has the same numerical value as MPa*s/m.
    implicit none
    type(layered_antiplane_model), intent(in) :: model
    real(8), intent(in) :: z
    real(8), intent(out) :: mu,rho,vs_local,eta
    integer :: l
    if(.not.allocated(model%densities)) error stop 'Radiation damping requires layer densities'
    if(.not.ieee_is_finite(z).or.z<0d0) error stop 'Invalid material evaluation depth'
    l=1
    do while(l<=size(model%bottoms))
      if(z<model%bottoms(l)) exit
      l=l+1
    end do
    mu=model%moduli(l); rho=model%densities(l)
    vs_local=sqrt(mu*1000d0/rho)
    eta=0.5d0*mu/vs_local
  end subroutine

  real(8) function matel2dvs_layered_ij(i,j,st_bemv) result(value)
    ! Midpoint collocation, unit CONSTANT source slip; no extra element width.
    ! Analytic homogeneous finite part + inverse Fourier layered correction.
    ! Self traction includes subtraction of the plastic slip jump, never zeroed.
    implicit none
    integer, intent(in) :: i,j
    type(st_HACApK_calc_entry), intent(in) :: st_bemv
    real(8) :: z,a,b,mu,q,spectrum,homogeneous
    integer :: k
    z=st_bemv%xcol(i); a=st_bemv%xel(j); b=st_bemv%xer(j)
    mu=layered_modulus(st_bemv%layered,z)
    value=mu/(2d0*pi)*(1d0/(z-b)-1d0/(z-a)-1d0/(z+b)+1d0/(z+a))
    if(all(st_bemv%layered%moduli==mu)) return
    do k=1,size(st_bemv%layered%wave)
      q=st_bemv%layered%wave(k)
      spectrum=layered_patch_spectrum(st_bemv%layered,q,a,b,z)
      if(z<a) then
        homogeneous=0.5d0*mu*(exp(-q*(a-z))-exp(-q*(b-z)))
      else if(z>b) then
        homogeneous=0.5d0*mu*(exp(-q*(z-b))-exp(-q*(z-a)))
      else
        homogeneous=mu*(1d0-0.5d0*(exp(-q*(z-a))+exp(-q*(b-z))))
      end if
      homogeneous=homogeneous+0.5d0*mu*(exp(-q*(z+a))-exp(-q*(z+b)))
      ! Both spectra contain the same -mu*slip plastic term: it cancels here.
      value=value+st_bemv%layered%weight(k)*(spectrum-homogeneous)
    end do
  end function

  real(8) function layered_modulus(model,z) result(mu)
    implicit none
    type(layered_antiplane_model), intent(in) :: model
    real(8), intent(in) :: z
    integer :: l
    l=1
    do while(l<=size(model%bottoms))
      if(z<model%bottoms(l)) exit
      l=l+1
    end do
    mu=model%moduli(l) ! At an interface use the deeper material.
  end function

  real(8) function layered_admittance(model,q,z,upward) result(y)
    ! Y=sigma23/(q*u3); free solution downwards, decaying solution upwards.
    implicit none
    type(layered_antiplane_model), intent(in) :: model
    real(8), intent(in) :: q,z
    logical, intent(in) :: upward
    real(8) :: lo,hi,mu,t
    integer :: l,n
    n=size(model%bottoms)
    if(upward) then
      y=-model%moduli(n+1)
      do l=n,1,-1
        lo=0d0
        if(l>1) lo=model%bottoms(l-1)
        lo=max(lo,z); hi=model%bottoms(l)
        if(hi<=lo) cycle
        mu=model%moduli(l); t=tanh(q*(hi-lo))
        y=(y-mu*t)/(1d0-y/mu*t)
      end do
    else
      y=0d0; lo=0d0
      do l=1,n+1
        hi=z
        if(l<=n) hi=min(z,model%bottoms(l))
        if(hi>lo) then
          mu=model%moduli(l); t=tanh(q*(hi-lo))
          y=(y+mu*t)/(1d0+y/mu*t)
        end if
        if(hi>=z) exit
        lo=hi
      end do
    end if
  end function

  real(8) function layered_log_gain(model,q,a,b,y_start,upward) result(gain)
    ! Only propagate in a growing direction; log amplitudes avoid overflow.
    implicit none
    type(layered_antiplane_model), intent(in) :: model
    real(8), intent(in) :: q,a,b,y_start
    logical, intent(in) :: upward
    real(8) :: lo,hi,mu,r,t,y,direction
    integer :: l,first,last,step,n
    n=size(model%bottoms); gain=0d0; y=y_start
    first=1; last=n+1; step=1; direction=1d0
    if(upward) then
      first=n+1; last=1; step=-1; direction=-1d0
    end if
    do l=first,last,step
      lo=0d0
      if(l>1) lo=model%bottoms(l-1)
      lo=max(a,lo); hi=b
      if(l<=n) hi=min(b,model%bottoms(l))
      if(hi<=lo) cycle
      mu=model%moduli(l); r=q*(hi-lo); t=tanh(r)
      gain=gain+r+log(0.5d0*(1d0+exp(-2d0*r)))+log(1d0+direction*y/mu*t)
      y=(y+direction*mu*t)/(1d0+direction*y/mu*t)
    end do
  end function

  real(8) function layered_line_traction(model,q,s,z,mu_source) result(value)
    implicit none
    type(layered_antiplane_model), intent(in) :: model
    real(8), intent(in) :: q,s,z,mu_source
    real(8) :: yf,yd,y,gain
    yf=layered_admittance(model,q,s,.false.)
    yd=layered_admittance(model,q,s,.true.)
    if(z<s) then
      y=layered_admittance(model,q,z,.false.)
      gain=layered_log_gain(model,q,z,s,y,.false.)
    else
      y=layered_admittance(model,q,z,.true.)
      gain=layered_log_gain(model,q,s,z,y,.true.)
    end if
    value=layered_modulus(model,z)*q*(mu_source/(yf-yd))*exp(-gain)
  end function

  real(8) function layered_patch_spectrum(model,q,a,b,z) result(value)
    ! Exact source-depth integral within each layer and each side of receiver.
    ! Integral of f(s)/f(hi) or d(s)/d(lo), using bounded hyperbolic ratios.
    implicit none
    type(layered_antiplane_model), intent(in) :: model
    real(8), intent(in) :: q,a,b,z
    real(8) :: lo,hi,mu,s,y,r,t,e,one_minus_sech,integral,end_hi
    integer :: l,side
    value=0d0
    do l=1,size(model%moduli)
      lo=0d0
      if(l>1) lo=model%bottoms(l-1)
      lo=max(a,lo); end_hi=b
      if(l<=size(model%bottoms)) end_hi=min(b,model%bottoms(l))
      if(end_hi<=lo) cycle
      mu=model%moduli(l)
      do side=1,2
        if(side==1) then
          hi=min(end_hi,z)
          if(hi<=lo) cycle
          s=hi; y=layered_admittance(model,q,lo,.false.)/mu
        else
          hi=end_hi; lo=max(lo,z)
          if(hi<=lo) cycle
          s=lo; y=-layered_admittance(model,q,hi,.true.)/mu
        end if
        r=q*(hi-lo); t=tanh(r); e=exp(-r)
        one_minus_sech=(1d0-e)**2/(1d0+e*e)
        integral=(t+y*one_minus_sech)/(q*(1d0+y*t))
        value=value+layered_line_traction(model,q,s,z,mu)*integral
      end do
    end do
  end function

  subroutine layered_gauss_legendre(nodes,weights)
    implicit none
    real(8), intent(out) :: nodes(:),weights(:)
    real(8) :: x,old,p0,p1,p2,derivative
    integer :: n,i,j,iteration
    n=size(nodes)
    do i=1,(n+1)/2
      x=cos(pi*(real(i,8)-0.25d0)/(real(n,8)+0.5d0))
      do iteration=1,100
        p0=1d0; p1=x
        do j=2,n
          p2=((2d0*j-1d0)*x*p1-(j-1d0)*p0)/j
          p0=p1; p1=p2
        end do
        derivative=n*(x*p1-p0)/(x*x-1d0)
        old=x; x=old-p1/derivative
        if(abs(x-old)<=4d0*epsilon(1d0)) exit
      end do
      if(iteration>100) error stop 'Gauss-Legendre iteration failed'
      nodes(i)=-x; nodes(n+1-i)=x
      weights(i)=2d0/((1d0-x*x)*derivative**2); weights(n+1-i)=weights(i)
    end do
  end subroutine

  real(8) function matel2dph2_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang,v)
  implicit none
  real(8)::angle,dx
  real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
  integer,intent(in)::i,j
  real(8)::sxx,sxy,syy,sxxp,sxyp,syyp,sxxq,sxyq,syyq,dip,s1,s2
  character(128)::v

  dip=ang(j)
  s1=dcos(dip)
  s2=dsin(dip)
  Call D2dip2(xcol(i),ycol(i),xel(j),yel(j),s1,s2,sxxp,syyp,sxyp)
  Call D2dip2(xcol(i),ycol(i),xer(j),yer(j),s1,s2,sxxq,syyq,sxyq)
  
  sxx=sxxp-sxxq
  syy=syyp-syyq
  sxy=sxyp-sxyq

  select case(v)
  case('s')
    matel2dph2_ij=-0.5d0*(sxx-syy)*dsin(2*ang(i))+sxy*dcos(2*ang(i))
  case('n')
    matel2dph2_ij=-0.5d0*(sxx+syy)+0.5d0*(sxx-syy)*dcos(2*ang(i))+sxy*dsin(2*ang(i))
  case('k')
    matel2dph2_ij=sxx+syy+pois*(sxx+syy)
  end select
end function

real(8) function load2dnh(xr,yr,xs,ys,ang,v)
  implicit none
  real(8)::angle,dx
  real(8),intent(in)::xr,yr,xs,ys,ang
  character(128),intent(in)::v
  real(8)::sxx,sxy,syy,sxxp,sxyp,syyp,sxxq,sxyq,syyq,dip,s1,s2

  dip=ang
  s1=dcos(dip)
  s2=dsin(dip)
  Call D2dip2(xr,yr,xs,ys,s1,s2,sxxp,syyp,sxyp)
  !Call D2dip2(xcol(i),ycol(i),xer(j),yer(j),s1,s2,sxxq,syyq,sxyq)
  
  sxx=sxxp!-sxxq
  syy=syyp!-syyq
  sxy=sxyp!-sxyq
  !write(*,*) sxx,sxy,syy


  select case(v)
  case('s')
    load2dnh=-0.5d0*(sxx-syy)*dsin(2*ang)+sxy*dcos(2*ang)
  case('n')
    load2dnh=-0.5d0*(sxx+syy)+0.5d0*(sxx-syy)*dcos(2*ang)+sxy*dsin(2*ang)
  end select
end function

Subroutine D2dip2(X,Y,xs,ys,s1,s2,sxx,syy,sxy)
  Implicit None
  ! input
  Real(8),intent(in):: X,Y,xs,ys,s1,s2
  ! output
  Real(8),intent(out):: sxx,sxy,syy
  real(8)::r1,r2,sxx1,sxx2,syy1,syy2,sxy1,sxy2,coef

  r1=sqrt((X-xs)**2+(Y-ys)**2)
  r2=sqrt((X-xs)**2+(Y+ys)**2)

  coef=rigid/2d0/pi/(1d0-pois)

  sxx1=s1*coef*((Y-ys)*((Y-ys)**2+3*(X-xs)**2)/r1**4-(Y+ys)*((Y+ys)**2+3*(X-xs)**2)/r2**4+2*ys/r2**6*(6*Y*(Y+ys)*(X-xs)**2-(Y-ys)*(Y+ys)**3-(X-xs)**4))
  sxx2=s2*coef*((X-xs)*((Y-ys)**2-(X-xs)**2)/r1**4-(X-xs)*((Y+ys)**2-(X-xs)**2)/r2**4 + 4*ys*(X-xs)/r2**6 *((2*ys-Y)*(Y+ys)**2+(3*Y+2*ys)*(X-xs)**2))
  syy1=s1*coef*((Y-ys)*((Y-ys)**2-(X-xs)**2)/r1**4-(Y+ys)*((Y+ys)**2-(X-xs)**2)/r2**4 - 2*ys/r2**6*(6*Y*(Y+ys)*(X-xs)**2-(3*Y+ys)*(Y+ys)**3+(X-xs)**4))
  syy2=-s2*coef*((X-xs)*(3*(Y-ys)**2+(X-xs)**2)/r1**4-(X-xs)*(3*(Y+ys)**2+(X-xs)**2)/r2**4 - 4*ys*(X-xs)*Y/r2**6 *(3*(Y+ys)**2-(X-xs)**2))
  sxy1=s1*coef*((X-xs)*((Y-ys)**2-(X-xs)**2)/r1**4-(X-xs)*((Y+ys)**2-(X-xs)**2)/r2**4 + 4*ys*(X-xs)*Y/r2**6*(3*(Y+ys)**2-(X-xs)**2))
  sxy2=s2*coef*((Y-ys)*((Y-ys)**2-(X-xs)**2)/r1**4-(Y+ys)*((Y+ys)**2-(X-xs)**2)/r2**4 + 2*ys/r2**6*(6*Y*(Y+ys)*(X-xs)**2-(X-xs)**4+(ys-Y)*(Y+ys)**3))
  sxx=sxx1+sxx2; syy=syy1+syy2; sxy=sxy1+sxy2

End Subroutine D2dip2

  ! real(8) function load2dnh(xs,ys,edge,dipangle,v)
  !   implicit none
  !   real(8)::angle,dx
  !   real(8),intent(in)::xs,ys,edge,dipangle
  !   real(8)::p22,p23,p33,u1,u2
  !   character(128),intent(in)::v
  !   Call D2dip(xs,ys,edge,1d8,dipangle*pi/180d0,1d0,p22,p23,p33)
  !   select case(v)
  !   case('xx')
  !     load2dnh=p22
  !   case('xy')
  !     load2dnh=p23
  !   case('yy')
  !     load2dnh=p33
  !   end select
  ! end function

  real(8) function matel2dph_ij(i,j,xcol,ycol,ang,v)
    implicit none
    real(8)::angle,dx
    real(8),intent(in)::xcol(:),ycol(:),ang(:)
    integer,intent(in)::i,j
    real(8)::p22,p23,p33,u1,u2
    character(128)::v

    angle=ang(1)
    dx=sqrt((xcol(1)-xcol(2))**2+(ycol(1)-ycol(2))**2)
    !write(*,*) angle,dx
    Call D2dip(xcol(i),ycol(i),(j-1)*dx,j*dx,angle,1d0,p22,p23,p33)
    select case(v)
    case('xx')
      matel2dph_ij=p22
    case('xy')
      matel2dph_ij=p23
    case('yy')
      matel2dph_ij=p33
    case('s')
      matel2dph_ij=0.5d0*(p22-p33)*dsin(-2*ang(i))+p23*dcos(-2*ang(i))
    case('n')
      matel2dph_ij=-(0.5d0*(p22+p33)-0.5d0*(p22-p33)*dcos(2*ang(i))-p23*dsin(2*ang(i)))
    end select
  end function

  Subroutine D2dip(x2,x3,s1,s2,angle,disl,p22,p23,p33)
    Implicit None
    ! input
    Real(8),intent(in):: x2,x3,s1,s2,angle,disl
    ! output
    Real(8),intent(out):: p22,p23,p33

    Real*8 p22_1,p22_2,p23_1,p23_2,p33_1,p33_2
    Real*8 u2_1,u2_2,u3_1,u3_2

    Call res(x2,x3,s1,angle,disl,p22_1,p23_1,p33_1)
    Call res(x2,x3,s2,angle,disl,p22_2,p23_2,p33_2)

    p22=p22_2-p22_1
    p23=p23_2-p23_1
    p33=p33_2-p33_1

  End Subroutine D2dip

  Subroutine res(x2,x3,s,angle,disl,p22,p23,p33)
    Implicit None
    Real(8),intent(in):: x2,x3,s,angle,disl
    Real(8),intent(out):: p22,p23,p33

    Real*8 r12,r22,r12i,r22i,r14,r24,r16,r26
    Real*8 cosd,sind,cos2d,sin2d
    Real*8 coef

    Real*8 p221,p222,p223,p224,p225
    Real*8 p231,p232,p233,p234,p235,p236
    Real*8 p331,p332,p333,p334,p335

    Real*8 coef2
    Real*8 u21,u22,u23,u24,u25,u26

    cosd=dcos(angle)
    sind=dsin(angle)
    cos2d=dcos(2d0*angle)
    sin2d=dsin(2d0*angle)
    r12=(x2-s*cosd)**2+(x3-s*sind)**2
    r22=(x2-s*cosd)**2+(x3+s*sind)**2
    r12i=1d0/r12
    r22i=1d0/r22
    r14=r12**2
    r24=r22**2
    r16=r12**3
    r26=r22**3

    coef=rigid*disl/2d0/pi/(1d0-pois)


    p221=(x2*sind-3d0*x3*cosd)*(r12i-r22i)
    p222=sin2d*(s*(r12i+3d0*r22i))
    p223=-2d0*(x2*sind-x3*cosd)*(((x3-s*sind)**2)/r14-((x3+s*sind)**2)/r24)
    p224=-4d0*sind*(s/r24)*(3d0*x2*x3*sind+5d0*x3**2*cosd+2d0*s*sind*(2d0*x3*cosd+x2*sind))
    p225=16d0*(x2*sind+x3*cosd)*x3*sind*(s*(x3+s*sind)**2/r26)

    p22=coef*(p221+p222+p223+p224+p225)


    p231=(x2*cosd-x3*sind)*(r12i-r22i)
    p232=-cos2d*(s*(r12i-r22i))
    p233=2d0*(x2*sind-x3*cosd)*(x2-s*cosd)*((x3-s*sind)/r14-(x3+s*sind)/r24)
    p234=4d0*sind*(x2*sind+2d0*x3*cosd)*(s*(x2-s*cosd)/r24)
    p235=4d0*x3*sind*sind*(s*(x3+s*sind)/r24)
    p236=-16d0*x3*sind*(x2*sind+x3*cosd)*(s*(x3+s*sind)*(x2-s*cosd)/r26)

    p23=coef*(p231+p232+p233+p234+p235+p236)


    p331=(x2*sind+x3*cosd)*(r12i-r22i)
    p332=-sin2d*(s*(r12i-r22i))
    p333=2d0*(x2*sind-x3*cosd)*(((x3-s*sind)**2)/r14-((x3+s*sind)**2)/r24)
    p334=4d0*x3*sind*(s*(x2*sind+3d0*x3*cosd+s*sin2d)/r24)
    p335=-16d0*x3*sind*(x2*sind+x3*cosd)*(s*(x3+s*sind)**2/r26)

    p33=coef*(p331+p332+p333+p334+p335)

    ! coef2=disl/2d0/pi
    !
    ! u21=(1d0-2d0*pois)/2d0/(1d0-pois)*sind*log(dsqrt(r12)/dsqrt(r22))
    ! u22=cosd*datan((x2-s*cosd)/(x3+s*sind))
    ! u23=-cosd*datan((x2-s*cosd)/(x3-s*sind))
    ! u24=-1d0/2d0/(1d0-pois)*(x2*sind-x3*cosd)*(x2-s*cosd)*(r12i-r22i)
    ! u25=2d0*s*sind*((1d0-2d0*pois)/2d0/(1d0-pois)*x3*sind+s-x2*cosd)/r22
    ! u26=2d0/(1d0-pois)
    return
  End Subroutine res

  real(8) function tensor2d_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang,v,md,rake)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:),rake(:)
    character(128),intent(in)::v,md
    real(8)::xp,xm,yp,ym,kern11,kern12,kern22,sin2,cos2,ds,sum_xx,sum_xy,sum_yy

    xp=cos(ang(j))*(xcol(i)-xel(j))+sin(ang(j))*(ycol(i)-yel(j))
    xm=cos(ang(j))*(xcol(i)-xer(j))+sin(ang(j))*(ycol(i)-yer(j))
    yp=-sin(ang(j))*(xcol(i)-xel(j))+cos(ang(j))*(ycol(i)-yel(j))
    ym=-sin(ang(j))*(xcol(i)-xer(j))+cos(ang(j))*(ycol(i)-yer(j))

    kern11=-0.5d0*rigid/vs*(inte11s(xp,yp)-inte11s(xm,ym))
    kern12=-0.5d0*rigid/vs*(inte12s(xp,yp)-inte12s(xm,ym))
    kern22=-0.5d0*rigid/vs*(inte22s(xp,yp)-inte22s(xm,ym))
    if(md.eq.'open')then
      kern11=-0.5d0*rigid/vs*(inte11o(xp,yp)-inte11o(xm,ym))
      kern12=-0.5d0*rigid/vs*(inte12o(xp,yp)-inte12o(xm,ym))
      kern22=-0.5d0*rigid/vs*(inte22o(xp,yp)-inte22o(xm,ym))
    end if

    !=>global coordinate system
    sin2=dsin(-2*ang(j))
    cos2=dcos(-2*ang(j))
    !select case(v)
    !case('xx')
      sum_xx=0.5d0*(kern11+kern22)+0.5d0*(kern11-kern22)*cos2+kern12*sin2
    !case('xy')
      sum_xy=-0.5d0*(kern11-kern22)*sin2+kern12*cos2
    !case('yy')
      sum_yy=0.5d0*(kern11+kern22)-0.5d0*(kern11-kern22)*cos2-kern12*sin2
    !end select

    select case(v)
    case('s')
      tensor2d_ij=(0.5d0*(sum_xx-sum_yy)*dsin(-2*ang(i))+sum_xy*dcos(-2*ang(i)))*rake(i)*rake(j)
    case('n')
      tensor2d_ij=-(0.5d0*(sum_xx+sum_yy)-0.5d0*(sum_xx-sum_yy)*dcos(2*ang(i))-sum_xy*dsin(2*ang(i)))*rake(j)
    end select
  end function

  real(8) function tensor2d_load(x,y,xsl,xsr,ysl,ysr,angs,v)
    implicit none
    real(8),intent(in)::x,y,xsl,xsr,ysl,ysr,angs
    character(128),intent(in)::v
    real(8)::xp,xm,yp,ym,kern11,kern12,kern22,sin2,cos2,ds

    xp=cos(angs)*(x-xsl)+sin(angs)*(y-ysl)
    xm=cos(angs)*(x-xsr)+sin(angs)*(y-ysr)
    yp=-sin(angs)*(x-xsl)+cos(angs)*(y-ysl)
    ym=-sin(angs)*(x-xsr)+cos(angs)*(y-ysr)
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
    sin2=dsin(-2*angs)
    cos2=dcos(-2*angs)
    select case(v)
    case('xx')
      tensor2d_load=0.5d0*(kern11+kern22)+0.5d0*(kern11-kern22)*cos2+kern12*sin2
    case('xy')
      tensor2d_load=-0.5d0*(kern11-kern22)*sin2+kern12*cos2
    case('yy')
      tensor2d_load=0.5d0*(kern11+kern22)-0.5d0*(kern11-kern22)*cos2-kern12*sin2
    end select
  end function

  real(8) function tensor2d3_ij(i,j,xcol,ycol,xel,xer,yel,yer,ang)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),ang(:)
    real(8)::xp,xm,yp,ym,kern31,kern32,sin2,cos2,ds

    xp=cos(ang(j))*(xcol(i)-xel(j))+sin(ang(j))*(ycol(i)-yel(j))
    xm=cos(ang(j))*(xcol(i)-xer(j))+sin(ang(j))*(ycol(i)-yer(j))
    yp=-sin(ang(j))*(xcol(i)-xel(j))+cos(ang(j))*(ycol(i)-yel(j))
    ym=-sin(ang(j))*(xcol(i)-xer(j))+cos(ang(j))*(ycol(i)-yer(j))

    kern31=-0.5d0*rigid/vs*(inte31s(xp,yp)-inte31s(xm,ym))
    kern32=-0.5d0*rigid/vs*(inte32s(xp,yp)-inte32s(xm,ym))
    !=>global coordinate system
    tensor2d3_ij=sin(ang(i))*kern31+cos(ang(i))*kern32
  end function

  real(8) function tensor2d3_load(x,y,xsl,xsr,ysl,ysr,angs)
    implicit none
    real(8),intent(in)::x,y,xsl,xsr,ysl,ysr,angs
    real(8)::xp,xm,yp,ym,kern31,kern32,sin2,cos2,ds

    xp=cos(angs)*(x-xsl)+sin(angs)*(y-ysl)
    xm=cos(angs)*(x-xsr)+sin(angs)*(y-ysr)
    yp=-sin(angs)*(x-xsl)+cos(angs)*(y-ysl)
    ym=-sin(angs)*(x-xsr)+cos(angs)*(y-ysr)

    kern31=-0.5d0*rigid/vs*(inte31s(xp,yp)-inte31s(xm,ym))
    kern32=-0.5d0*rigid/vs*(inte32s(xp,yp)-inte32s(xm,ym))
    !=>global coordinate system
    tensor2d3_load=sin(angs)*kern31+cos(angs)*kern32

  end function

  real(8) function matel3dp_ij(i,j,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4,v,md)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),zcol(:)
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
    character(128),intent(in)::md,v

    matel3dp_ij=ret3dp(xcol(i)-xs1(j),xcol(i)-xs2(j),xcol(i)-xs3(j),xcol(i)-xs4(j),&
    & zcol(i)-zs1(j),zcol(i)-zs2(j),zcol(i)-zs3(j),zcol(i)-zs4(j),v,md)

  end function
  real(8) function matel3dph_ij(i,j,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4,v,md)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),zcol(:)
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
    character(128),intent(in)::md,v
    real(8)::tmp1,tmp2

    tmp1=ret3dp(xcol(i)-xs1(j),xcol(i)-xs2(j),xcol(i)-xs3(j),xcol(i)-xs4(j),&
    & zcol(i)-zs1(j),zcol(i)-zs2(j),zcol(i)-zs3(j),zcol(i)-zs4(j),v,md)
    tmp2=ret3dp(xcol(i)-xs1(j),xcol(i)-xs2(j),xcol(i)-xs3(j),xcol(i)-xs4(j),&
    & zcol(i)+zs1(j),zcol(i)+zs2(j),zcol(i)+zs3(j),zcol(i)+zs4(j),v,md)
    matel3dph_ij=tmp1-tmp2
   ! write(*,*)tmp1,tmp2
    !write(*,*) matel3dph_ij
  end function

  real(8) function ret3dp(dx1,dx2,dx3,dx4,dz1,dz2,dz3,dz4,v,md)
    implicit none
    real(8),intent(in)::dx1,dx2,dx3,dx4,dz1,dz2,dz3,dz4
    character(128),intent(in)::md,v
    real(8)::angle,dx,dy,sxx,sxy,syy,gtau(2),gsig(2),factor
    real(8)::dr1,dr2,dr3,dr4,ret1,ret2,ret3,alpha
    ret3dp = 0.0

    factor = rigid/(4.d0 * pi)
    dr1=dsqrt(dx1**2+dz1**2)
    dr2=dsqrt(dx2**2+dz2**2)
    dr3=dsqrt(dx3**2+dz3**2)
    dr4=dsqrt(dx4**2+dz4**2)

    alpha=1/(2-2*pois)
    !alpha=0.5

    select case(md)
    case("s","o")
      if(v=="s") then
        ret1=2.d0*(1-alpha)*(dr4/(dx4*dz4)+dr2/(dx2*dz2)-dr3/(dx3*dz3)-dr1/(dx1*dz1))
        ret2=-(1-2*alpha)/dx1*(dr4/dz4+dz4/dr4-dr1/dz1-dz1/dr1)
        ret3=-(1-2*alpha)/dx2*(dr2/dz2+dz2/dr2-dr3/dz3-dz3/dr3)
        ret3dp = factor*(ret1+ret2+ret3)
      end if
    case('open')
      if(v=="n") then
        ret1=2.d0*alpha*(dr4/(dx4*dz4)+dr2/(dx2*dz2)-dr3/(dx3*dz3)-dr1/(dx1*dz1))
        ret3dp = factor*ret1
      end if
    end select

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
      !I1(1)=gx(1)/r(1)-gy(1)/r(1)*(2*nx*ny*(gy(1)**2-gx(1)**2)-2*(ny**2-nx**2)*gx(1)*gy(1))
      I1(1)=gx(1)/r(1)+gy(1)/r(1)*(2*nx*ny*(gy(1)**2-gx(1)**2)-2*(ny**2-nx**2)*gx(1)*gy(1))
      I2(1)=gy(1)/r(1)+gx(1)/r(1)*(2*nx*ny*(gy(1)**2-gx(1)**2)-2*(ny**2-nx**2)*gx(1)*gy(1))
      !I1(2)=gx(2)/r(2)-gy(2)/r(2)*(2*nx*ny*(gy(2)**2-gx(2)**2)-2*(ny**2-nx**2)*gx(2)*gy(2))
      I1(2)=gx(2)/r(2)+gy(2)/r(2)*(2*nx*ny*(gy(2)**2-gx(2)**2)-2*(ny**2-nx**2)*gx(2)*gy(2))
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
    angp=atan2(dy,dx)!*dx/abs(dx)
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
    angp=atan2(dy,dx)
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
    angp=atan2(dy,dx)
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

  real(8) function matels2dp_ij(i,j,xcol,xel,xer,v,md)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),xel(:),xer(:)
    character(128),intent(in)::v,md
    !real(8),intent(in)::rigid,pois
    real(8)::factor

    !inplane shear
    factor=0.d0
    if(v.eq.'s'.and.md.eq.'s') then
       factor=rigid/(2.d0*pi*(1.d0-pois))
       else if(v.eq.'n'.and.md.eq.'open') then
       factor=-rigid/(2.d0*pi*(1.d0-pois))
    end if
    !antiplane shear
    !factor=rigid/(2.d0*pi)

    matels2dp_ij=factor*(1.d0/(xcol(i)-xer(j))-1.d0/(xcol(i)-xel(j)))
  end function matels2dp_ij

  real(8) function matels2dpFZ_ij(i,j,xcol,xel,xer,fwid)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),xel(:),xer(:),fwid
    !real(8),intent(in)::rigid,pois
    integer::k
    real(8)::factor,a,dx1,dx2

    !antiplane shear
    factor=rigid/(2.d0*pi)

    matels2dpFZ_ij=1.d0/(xcol(i)-xer(j))-1.d0/(xcol(i)-xel(j))
  do k=1,10
    a=k*fwid
    dx1=xcol(i)-xer(j)
    dx2=xcol(i)-xel(j)
    matels2dpFZ_ij=matels2dpFZ_ij+2*dx1*(dx1**3-3*a**2)/(a**2+dx1**2)**3-2*dx2*(dx2**3-3*a**2)/(a**2+dx2**2)**3
  end do   
  matels2dpFZ_ij=matels2dpFZ_ij*factor

  end function matels2dpFZ_ij

  real(8) function matels2dpa_ij(i,j,xcol,xel,xer)
  implicit none
  integer,intent(in)::i,j
  real(8),intent(in)::xcol(:),xel(:),xer(:)
  !real(8),intent(in)::rigid,pois
  real(8)::factor

  !antiplane shear
  factor=rigid/(2.d0*pi)
  matels2dpa_ij=factor*(1.d0/(xcol(i)-xer(j))-1.d0/(xcol(i)-xel(j)))
end function matels2dpa_ij

  function inte12s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte12s,r
    r=sqrt(x1**2+x2**2)
    pa=sqrt((1-2*pois)/(2-2*pois))
    inte12s=2*vs/pi*(1-pa**2)*x1*(x1**2-x2**2)/r**4
    return
  end function

  function inte11s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte11s,r
    r=sqrt(x1**2+x2**2)
    pa=sqrt((1-2*pois)/(2-2*pois))
    inte11s=-2*vs/pi*(1-pa**2)*x2*(3*x1**2+x2**2)/r**4
    return
  end function

  function inte22s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte22s,r
    r=sqrt(x1**2+x2**2)
    pa=sqrt((1-2*pois)/(2-2*pois))
    inte22s=2*vs/pi*(1-pa**2)*x2*(x1**2-x2**2)/r**4
    return
  end function

  function inte12o(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte12o,r
    r=sqrt(x1**2+x2**2)
    pa=sqrt((1-2*pois)/(2-2*pois))
    inte12o=2*vs/pi*(1-pa**2)*x2*(x1**2-x2**2)/r**4
    return
  end function

  function inte11o(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte11o,r
    r=sqrt(x1**2+x2**2)
    pa=sqrt((1-2*pois)/(2-2*pois))
    inte11o=2*vs/pi*(1-pa**2)*x1*(x1**2-x2**2)/r**4
    return
  end function

  function inte22o(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte,pa,inte22o,r
    r=sqrt(x1**2+x2**2)
    pa=sqrt((1-2*pois)/(2-2*pois))
    inte22o=2*vs/pi*(1-pa**2)*x1*(3*x1**2+x2**2)/r**4
    return
  end function

  function inte31s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,pa,r,inte31s
    r=sqrt(x1**2+x2**2)
    inte31s=-vs*x2/pi/r**2
    return
  end function

  function inte32s(x1,x2)
    implicit none
    real(8)::x1,x2,t,ss,sp,pa,r,inte32s
    r=sqrt(x1**2+x2**2)
    inte32s=vs/pi*x1/r**2
    return
  end function

  real(8) function matels1_ij(i,j,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,v,md)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),zcol(:),xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    character(128),intent(in)::v,md
    real(8)::P1(3),P2(3),P3(3)
    real(8)::Sxx,Syy,Szz,Sxy,Sxz,Syz

    P1=(/xs1(j),ys1(j),zs1(j)/)
    P2=(/xs2(j),ys2(j),zs2(j)/)
    P3=(/xs3(j),ys3(j),zs3(j)/)
    select case(md)
    case('st')
      call TDstressFS(xcol(i),ycol(i),zcol(i),P1,P2,P3,1.d0,0.d0,0.d0,rigid,rigid,&
      & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    case('dp')
      call TDstressFS(xcol(i),ycol(i),zcol(i),P1,P2,P3,0.d0,1.d0,0.d0,rigid,rigid,&
      & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    end select
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
  real(8) function matel3dn_ij(i,j,st_bemv)
    implicit none
    integer,intent(in)::i,j
    type(st_HACApK_calc_entry)::st_bemv
    matel3dn_ij=geological_triangle_entry(i,j,st_bemv,.false.)
  end function matel3dn_ij
  real(8) function matel3dh_ij(i,j,st_bemv)
    implicit none
    integer,intent(in)::i,j
    type(st_HACApK_calc_entry)::st_bemv
    matel3dh_ij=geological_triangle_entry(i,j,st_bemv,.true.)
  end function matel3dh_ij

  real(8) function matelh1_ij(i,j,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,v,md)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),zcol(:),xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    character(128),intent(in)::v,md
    real(8)::P1(3),P2(3),P3(3)
    real(8)::Sxx,Syy,Szz,Sxy,Sxz,Syz

    P1=(/xs1(j),ys1(j),zs1(j)/)
    P2=(/xs2(j),ys2(j),zs2(j)/)
    P3=(/xs3(j),ys3(j),zs3(j)/)
    select case(md)
    case('st')
      call TDstressHS(xcol(i),ycol(i),zcol(i),P1,P2,P3,1.d0,0.d0,0.d0,rigid,rigid,&
      & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    case('dp')
      call TDstressHS(xcol(i),ycol(i),zcol(i),P1,P2,P3,0.d0,1.d0,0.d0,rigid,rigid,&
      & Sxx,Syy,Szz,Sxy,Sxz,Syz)
    end select
    select case(v)
    case('xx')
      matelh1_ij=Sxx
    case('yy')
      matelh1_ij=Syy
    case('zz')
      matelh1_ij=Szz
    case('xy')
      matelh1_ij=Sxy
    case('xz')
      matelh1_ij=Sxz
    case('yz')
      matelh1_ij=Syz
    end select
    return
  end function matelh1_ij

  real(8) function matelrec_ij(i,j,xcol,ycol,ds,strike,w,v)
    implicit none
    integer,intent(in)::i,j
    real(8),intent(in)::xcol(:),ycol(:),ds(:),strike(:),w
    character(128),intent(in)::v
    integer::iret
    real(8)::dx,dy,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,sxx,syy,szz,sxy,sxz,syz,alpha,sxx2,sxy2,syy2
    real(8)::exx,eyy,ezz,exy,eyz,ezx

    alpha=(1d0+(0.5d0/pois-1d0))/(1d0+2d0*(0.5d0/pois-1d0))
    !rotation so that strike is parallel to y axis
    dx=cos(strike(j))*(xcol(i)-xcol(j))+sin(strike(j))*(ycol(i)-ycol(j))
    dy=-sin(strike(j))*(xcol(i)-xcol(j))+cos(strike(j))*(ycol(i)-ycol(j))

    call dc3d25(alpha,dx,dy,ds(j),w,sxx,sxy,syy)

    !rerotation
    select case(v)
    case('s')
      matelrec_ij=-0.5d0*(sxx-syy)*dsin(2*(strike(i)-strike(j)))+sxy*dcos(2*(strike(i)-strike(j)))
    case('n')
      matelrec_ij=-0.5d0*(sxx+syy)+0.5d0*(sxx-syy)*dcos(2*(strike(i)-strike(j)))+sxy*dsin(2*(strike(i)-strike(j)))
      !matelrec_ij=sin(strike(i))**2*sxx2+cos(strike(i))**2*syy2+sxy2*sin(2*strike(i))
    end select


    return
  end function matelrec_ij
  
  subroutine  dc3d25(alpha,x,y,dl,dw,sxx,sxy,syy)
    implicit none
    real(8),intent(in)::alpha,x,y,dl,dw
    real(8),intent(out)::sxx,sxy,syy
    real(8)::ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    integer::i,j,k
    real(8)::  xi(2),et(2),d,p,q
    real(8)::  u(12),du(12),dua(12),exx,exy,eyy,ezz

    xi(1)=x+0.5d0*dl
    xi(2)=x-0.5d0*dl
    q=y
    p=0.d0
    et(1)=p+0.5*dw
    et(2)=p-0.5*dw
 

    u=0d0
    do k=1,2
      do j=1,2
        call ua25(alpha,xi(j),et(k),q,dua)
        !du(1)  =-dua(1)
        !du(2)=dua(3)
        !du(3)=-dua(2)
        du(4)  =-dua(4) !ux,x
        du(5)=dua(6) !uy,x
        !du(6)=-dua(5)
        du(7)  =-dua(7)!ux,y
        du(8)=dua(9)!uy,y
        !du(9)=-dua(8)
        !du(10)  =dua(10)
        !du(11)=dua(12)
        du(12)=dua(11)!uz,z
        do i=1,12
          if(j+k.ne.3) u(i)=u(i)+du(i)
          if(j+k.eq.3) u(i)=u(i)-du(i)
        end do
        !----
      end do
    end do

    ! !====
    ux=u(1)
    uy=u(2)
    uz=u(3)
    uxx=u(4)
    uyx=u(5)
    uzx=u(6)
    uxy=u(7)
    uyy=u(8)
    uzy=u(9)
    uxz=u(10)
    uyz=u(11)
    uzz=u(12)

    exx=-uxx
    eyy=-uyy
    ezz=-uzz
    exy=-0.5d0*(uxy+uyx)
    !eyz=-0.5d0*(uyz+uzy)
    !ezx=-0.5d0*(uzx+uxz)
    sxx=rigid*(exx+eyy+ezz)+2*rigid*exx
    syy=rigid*(exx+eyy+ezz)+2*rigid*eyy
    !szz=rigid*(exx+eyy+ezz)+2*rigid*ezz
    sxy=2*rigid*exy
    return

  end subroutine dc3d25
  subroutine  ua25(alpha,xi,et,q,du)
    implicit none
    real(8),intent(in)::xi,et,q,alpha
    !integer,intent(in)::kxi,ket
    real(8),intent(out)::du(12)
    real(8)::qx,qy,xy,rxi,ret,r2,r,x11,zle,y11,y32,yp,ale
    integer::i

    yp=q
    r2=xi**2+et**2+q**2
    r =dsqrt(r2)
    rxi=r+xi
    x11=1d0/(r*rxi)
    ret=r+et
    y11=1d0/(r*ret)
    y32=(r+ret)*y11*y11/r
    !----

    !du( 1)= datan2(xi*et,q*r)/2d0 +alpha/2*xi*q*y11
    !du( 2)= alpha/2*q/r
    !du( 3)= (1d0-alpha)/2d0*dlog(r+et)-alpha/2*q*q*y11
    du( 4)=-(1d0-alpha)/2d0*q*y11-alpha/2*xi**2*q*y32 !-ux,x
    !du( 5)=-alpha/2*xi*q/r**3
    du( 6)= (1d0-alpha)/2d0*xi*y11+alpha/2*xi*q**2*y32 !uy,x
    du( 7)= (1d0-alpha)/2d0*xi*y11+alpha/2*xi*(et/r**3+xi**2*y32) +et/2d0*x11 !-ux,y
    !du( 8)= alpha/2*(1d0/r-q*q/r**3)
    du( 9)= (1d0-alpha)/2d0*(q*y11) -alpha/2*q*(et/r**3+xi**2*y32) !uy,y
    !du(10)= alpha/2*xi*(yp/r**3)+yp/2d0*x11
    du(11)= alpha/2*(et*q/r**3) !uz,z
    !du(12)=-(1d0-alpha)/2d0*(1d0/r-q*y11) -alpha/2*q*(yp/r**3+xi**2*y32)
    !write(*,*)'du',du
    du=du/(pi*2)

    return
  end subroutine
  real(8) function okada_ij(i,j,xcol,ycol,zcol,ang,angd,dsl,dsd,v,rake,w,fullspace)
    implicit none
    integer,intent(in)::i,j
    !type(st_HACApK_calc_entry) :: st_bemv
    real(8),intent(in)::xcol(:),ycol(:),zcol(:),ang(:),angd(:),dsl(:),dsd(:),rake(:),w
    character(128),intent(in)::v
    logical,intent(in)::fullspace
    integer::iret
    real(8)::dx,dy,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,sxx,syy,szz,sxy,sxz,syz,alpha
    real(8)::exx,eyy,ezz,exy,eyz,ezx,rotang,dpang,Arot(3,3),p(6),rr,z,depth,dip,rs,rd,lambda

    alpha=(1d0+(0.5d0/pois-1d0))/(1d0+2d0*(0.5d0/pois-1d0))
    !rotation so that strike is parallel to y axis
   ! write(*,*) st_bemv%strike(j)
   ! dx=cos(st_bemv%ang(j))*(st_bemv%xcol(i)-st_bemv%xcol(j))+sin(st_bemv%ang(j))*(st_bemv%ycol(i)-st_bemv%ycol(j))
   ! dy=-sin(st_bemv%ang(j))*(st_bemv%xcol(i)-st_bemv%xcol(j))+cos(st_bemv%ang(j))*(st_bemv%ycol(i)-st_bemv%ycol(j))
     dx=cos(ang(j))*(xcol(i)-xcol(j))+sin(ang(j))*(ycol(i)-ycol(j))
    dy=-sin(ang(j))*(xcol(i)-xcol(j))+cos(ang(j))*(ycol(i)-ycol(j))

    !select case(md)
    !case('st')
    !call dc3d(alpha,dx,dy,ds(j),w,&
    !& ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)
    !w=st_bemv%w
    z=zcol(i)
    depth=-zcol(j)
    dip=angd(j)*180/pi
    rs=dsl(j)/w
    rd=dsd(j)/w

    !rake=st_bemv%rake
   !st_bemv%dip(j)=pi/180
   !write(*,*)w
   !

    call okada(alpha,dx/w,dy/w,z/w,depth/w,dip,-0.5d0*rs,0.5d0*rs,-0.5d0*rd,0.5d0*rd,cos(rake(j)),sin(rake(j)),0d0,&
   &ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,fullspace)
  ! write(*,*) st_bemv%zcol(i),-st_bemv%zcol(j),st_bemv%angd(j)*180/pi
   ! write(*,*)uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    ! write(*,*)iret    ! case('dp')
    !   call dc3d(alpha,dx,dy,zcol(i),zcol(j),dip(j),-0.5d0*ds,0.5d0*ds,-0.5d0*ds,0.5d0*ds,0d0,1d-8,0d0, ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    !end select
    exx=uxx
    eyy=uyy
    ezz=uzz
    exy=0.5d0*(uxy+uyx)
    eyz=0.5d0*(uyz+uzy)
    ezx=0.5d0*(uzx+uxz)
    ! Isotropic Hooke law: lambda equals rigid only for pois=0.25.
    lambda=2*rigid*pois/(1-2*pois)
    sxx=lambda*(exx+eyy+ezz)+2*rigid*exx
    syy=lambda*(exx+eyy+ezz)+2*rigid*eyy
    szz=lambda*(exx+eyy+ezz)+2*rigid*ezz
    sxy=2*rigid*exy
    sxz=2*rigid*ezx
    syz=2*rigid*eyz

    !rotang=ang(i)-ang(j) !wrong
    rotang=ang(j)-ang(i) !fixed Nov 24 2024
    dpang=angd(i)

    Arot(:,1)=(/cos(rotang),-sin(rotang),0d0/)
    Arot(:,2)=(/sin(rotang)*cos(dpang),cos(rotang)*cos(dpang),sin(dpang)/)
    Arot(:,3)=(/sin(rotang)*sin(dpang),cos(rotang)*sin(dpang),-cos(dpang)/)
   ! write(*,*)Arot

   ! write(*,*)Arot
   call TensTrans(Sxx,Syy,Szz,Sxy,Sxz,Syz,Arot,&
   &p(1),p(2),p(3),p(4),p(5),p(6))

   !write(*,*)st_bemv%v
    select case(v)
      case('n')
        ! Compression positive, consistent with the 2-D and triangular entries.
        okada_ij=-p(3)
      case('s','o')
        ! The local third axis is downward; geological traction uses its opposite.
        okada_ij=-p(5)*cos(rake(i))-p(6)*sin(rake(i))
    end select
 ! write(*,*) okada_ij
    ! rr=(st_bemv%xcol(i)-st_bemv%xs1(j))**2+(st_bemv%zcol(i)-st_bemv%zs1(j))**2
    ! okada_ij=1/rr
   return

  end function okada_ij

  real(8) function okada_load(x,y,z,angs,angd,rake,xs1,ys1,lengx,depth,cangs,cangd,crake,v)
    implicit none
    real(8),intent(in)::x,y,z,rake,angs,angd,xs1,ys1,lengx,depth,cangs,cangd,crake
    character(1),intent(in)::v
    logical::fullspace=.false.
    integer::iret
    real(8)::dx,dy,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,sxx,syy,szz,sxy,sxz,syz,alpha
    real(8)::exx,eyy,ezz,exy,eyz,ezx,rotang,dpang,Arot(3,3),p(6),rr,dip,lambda,fwid=100d0
    alpha=(1d0+(0.5d0/pois-1d0))/(1d0+2d0*(0.5d0/pois-1d0))
    !rotation so that strike is parallel to y axis
    dx=cos(cangs)*(x-xs1)+sin(cangs)*(y-ys1)
    dy=-sin(cangs)*(x-xs1)+cos(cangs)*(y-ys1)

    dip=cangd*180/pi

    !write(*,*)'okada_load',dx,dy,z,depth,ang,dip,crake,rake
    call okada(alpha,dx,dy,z,depth,dip,-lengx/2,lengx/2,-fwid,0d0,cos(crake),sin(crake),0d0,&
   &ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,fullspace)

   !write(*,*)'okada_load',ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz

    exx=uxx
    eyy=uyy
    ezz=uzz
    exy=0.5d0*(uxy+uyx)
    eyz=0.5d0*(uyz+uzy)
    ezx=0.5d0*(uzx+uxz)
    ! Use the same elastic constants and normal-stress convention as okada_ij.
    lambda=2*rigid*pois/(1-2*pois)
    sxx=lambda*(exx+eyy+ezz)+2*rigid*exx
    syy=lambda*(exx+eyy+ezz)+2*rigid*eyy
    szz=lambda*(exx+eyy+ezz)+2*rigid*ezz
    sxy=2*rigid*exy
    sxz=2*rigid*ezx
    syz=2*rigid*eyz
    
    !write(*,*)sxx,syy,szz,sxy,syz,sxz

    rotang=cangs-angs !fixed Nov 24 2024
    dpang=angd
    !write(*,*)'okada_load',rotang,dpang

    Arot(:,1)=(/cos(rotang),-sin(rotang),0d0/)
    Arot(:,2)=(/sin(rotang)*cos(dpang),cos(rotang)*cos(dpang),sin(dpang)/)
    Arot(:,3)=(/sin(rotang)*sin(dpang),cos(rotang)*sin(dpang),-cos(dpang)/)
   ! write(*,*)Arot

   ! write(*,*)Arot
   call TensTrans(Sxx,Syy,Szz,Sxy,Sxz,Syz,Arot,&
   &p(1),p(2),p(3),p(4),p(5),p(6))
   !write(*,*)'okada_load',p(1),p(2),p(3),p(4),p(5),p(6)

   select case(v)
      case('n')
        okada_load=-p(3)
      case('s','o')
        okada_load=-p(5)*cos(rake)-p(6)*sin(rake)
    end select

   return
 end function okada_load
 subroutine okada_surface(xst,yst,x,y,z,ang,angd,rake,dsl,dsd,ux,uy,uz)
    implicit none
    real(8),intent(in)::x,y,z,xst,yst,ang,angd,rake,dsl,dsd
    real(8),intent(out)::ux,uy,uz
    logical::fullspace=.false.
    integer::iret
    real(8)::dx,dy,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,sxx,syy,szz,sxy,sxz,syz,alpha,ux1,uy1
    real(8)::exx,eyy,ezz,exy,eyz,ezx,rotang,dpang,Arot(3,3),p(6),rr,dip,fwid=100d0
    
    alpha=(1d0+(0.5d0/pois-1d0))/(1d0+2d0*(0.5d0/pois-1d0))

    !rotation so that strike is parallel to y axis
    dx=cos(ang)*(xst-x)+sin(ang)*(yst-y)
    dy=-sin(ang)*(xst-x)+cos(ang)*(yst-y)
    dip=angd*180/pi

    call okada(alpha,dx,dy,0.d0,-z,dip,-0.5d0*dsl,0.5d0*dsl,-0.5d0*dsd,0.5d0*dsd,cos(rake),sin(rake),0d0,&
   &ux1,uy1,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,fullspace)

    !rotation so that strike is parallel to y axis
    ux=ux1*cos(ang)-uy1*sin(ang)
    uy=ux1*sin(ang)+uy1*cos(ang)
    
    ! Arot(:,1)=(/cos(rotang),-sin(rotang),0d0/)
    ! Arot(:,2)=(/sin(rotang)*cos(dpang),cos(rotang)*cos(dpang),sin(dpang)/)
    ! Arot(:,3)=(/sin(rotang)*sin(dpang),cos(rotang)*sin(dpang),-cos(dpang)/)
 
   return
 end subroutine okada_surface
 subroutine surface_disp_matrix(matsurf,st_bemv,xst,yst)
    implicit none
    real(8),intent(in)::xst(:),yst(:)
    real(8),intent(out)::matsurf(:,:,:) 
    integer::nstation,ncellg,i,j
    type(st_HACApK_calc_entry) :: st_bemv
    nstation=size(xst)
    ncellg=size(st_bemv%xs1)
 
    do i=1,ncellg
      do j=1,nstation
        call okada_surface(xst(j),yst(j),st_bemv%xcol(i),st_bemv%ycol(i),st_bemv%zcol(i),st_bemv%ang(i),st_bemv%angd(i),&
        & st_bemv%rake(i),st_bemv%dsl(i),st_bemv%dsd(i),matsurf(j,1,i),matsurf(j,2,i),matsurf(j,3,i))
      end do
    end do

 end subroutine surface_disp_matrix
  end module
