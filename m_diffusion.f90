module mod_diffusion
use mod_constant
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  type :: t_params
  integer::nwell,nn,npoint,n1,i1,i2,nfault,ninj,ncell
  integer::nconnect=0
  integer,pointer::kleng(:),iwell(:),jwell(:),connect(:,:),ns(:)
  real(8)::phi0,beta,eta,sigmastar,kp0,kpmin,kpmax,kL,kT,pinj,pbcl,pbcr,pbct,pbcb,qinj,q0
  real(8)::qbcl,qbcr,qbct,qbcb
  real(8)::w0=1d0 ! Hydraulic fault-zone thickness [m].
  real(8),allocatable::wG(:) ! Original geometry order; fixed in time.
  real(8)::tinj=1d5,xinj=0d0,yinj=0d0,zinj=0d0,sigmaHF=0.0
  real(8),pointer::kp(:),kpG(:),qtimes(:),qvals(:,:),pfhyd(:,:),phi(:),phiG(:),dep(:),qs(:)
  character(128)::bc,bcl,bcr,bct,bcb,setting,injection,injection_file,network_file
  logical::injectionfromfile,switch,permev,permsigma,network,initialize_steadystate
  end type t_params
contains
  subroutine diffusion_width(param_diff,n,width)
    implicit none
    type(t_params),intent(in)::param_diff
    integer,intent(in)::n
    real(8),intent(out)::width(n)
    width=param_diff%w0
    if(allocated(param_diff%wG)) then
      if(size(param_diff%wG)/=n) error stop 'Hydraulic w array has wrong size'
      width=param_diff%wG
    end if
    if(any(.not.ieee_is_finite(width)).or.any(width<=0d0)) &
      error stop 'Hydraulic w must be finite and positive [m]'
  end subroutine

  real(8) function well_rate(param_diff,kwell,time) result(q)
    implicit none
    type(t_params),intent(in)::param_diff
    integer,intent(in)::kwell
    real(8),intent(in)::time
    integer::j,n
    n=param_diff%npoint
    if(n<2) error stop 'Well history needs at least two time points'
    q=0d0
    if(time<param_diff%qtimes(1)) return
    q=param_diff%qvals(kwell,n)
    do j=1,n-1
      if(param_diff%qtimes(j+1)<param_diff%qtimes(j)) error stop 'Well times must not decrease'
      if(time<param_diff%qtimes(j+1)) then
        q=param_diff%qvals(kwell,j)+(param_diff%qvals(kwell,j+1)-param_diff%qvals(kwell,j))* &
          (time-param_diff%qtimes(j))/(param_diff%qtimes(j+1)-param_diff%qtimes(j))
        exit
      end if
    end do
    q=max(0d0,q) ! Injection-only convention retained; raw volume rate, not Q/w.
  end function

  subroutine diffusion_next_knot(param_diff,time,dtnxt)
    implicit none
    type(t_params),intent(in)::param_diff
    real(8),intent(in)::time
    real(8),intent(inout)::dtnxt
    integer::j
    if(.not.param_diff%injectionfromfile) return
    do j=1,param_diff%npoint
      if(param_diff%qtimes(j)>time) then
        dtnxt=min(dtnxt,param_diff%qtimes(j)-time)
        exit
      end if
    end do
  end subroutine

  subroutine conservative_diffusion(nx,ny,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter,dep,surface)
    ! Backward Euler finite volumes: (M+h*L) p_new=M*p_old+source.
    ! M=beta*phi*w, T=k*w/eta_f. Pressure MPa, ds0 km, w m.
    ! Harmonic face mobility and storage-weighted SPD matrix; Jacobi PCG.
    implicit none
    integer,intent(in)::nx,ny
    real(8),intent(in)::ds0,pf(:),cdiff(:),str(:),h,time
    real(8),intent(out)::pfnew(:)
    integer,intent(out)::niter
    type(t_params),intent(in)::param_diff
    real(8),optional,intent(in)::dep(:)
    logical,optional,intent(in)::surface
    logical::two_dimensional
    real(8)::width(size(pf)),mass(size(pf)),mobility(size(pf)),diag(size(pf)),rhs(size(pf))
    real(8)::r(size(pf)),z(size(pf)),p(size(pf)),ap(size(pf)),edge(2*size(pf)+4*param_diff%nconnect)
    integer::left(2*size(pf)+4*param_diff%nconnect),right(2*size(pf)+4*param_diff%nconnect)
    real(8)::dx,area,scale,rz,rznew,alpha,denom,tol,q,g
    integer::n,i,j,k,u,v,ne,first,last,kw,iteration,itermax
    n=size(pf); ne=0; niter=0
    two_dimensional=ny>1
    if(present(surface)) two_dimensional=surface
    if(nx*ny/=n.or.size(cdiff)/=n.or.size(str)/=n) error stop 'Invalid diffusion dimensions'
    if(.not.ieee_is_finite(ds0).or.ds0<=0d0.or..not.ieee_is_finite(h).or.h<0d0) &
      error stop 'Invalid diffusion spacing or step'
    if(any(.not.ieee_is_finite(cdiff)).or.any(cdiff<0d0).or. &
      any(.not.ieee_is_finite(str)).or.any(str<=0d0)) error stop 'Invalid diffusion transport/storage'
    call diffusion_width(param_diff,n,width)
    dx=ds0*1000d0; area=dx
    if(two_dimensional) area=dx*dx
    mass=str*width
    mobility=cdiff*1d6*mass ! cdiff supplied in km^2/s; mobility in m^3/(Pa s).
    scale=maxval(mass)
    diag=mass/scale; rhs=diag*pf
    if(present(dep)) then
      if(ny/=1.or.size(dep)/=n) error stop 'Invalid network dimensions'
      if(param_diff%ns(1)/=0.or.param_diff%ns(param_diff%nfault+1)/=n) &
        error stop 'Network element counts do not match geometry'
      do k=1,param_diff%nfault
        first=param_diff%ns(k)+1; last=param_diff%ns(k+1)
        do i=first,last-1
          call add_edge(i,i+1)
        end do
        if(dep(first)<0.1d0) call boundary(first,'Dirichlet',0d0,0d0)
        if(dep(last)<0.1d0) call boundary(last,'Dirichlet',0d0,0d0)
      end do
      do k=1,param_diff%nconnect
        u=param_diff%connect(k,1); v=param_diff%connect(k,2)
        if(u<1.or.u>=n.or.v<1.or.v>=n) error stop 'Invalid network connection'
        call add_edge(u,v); call add_edge(u+1,v)
        call add_edge(u,v+1); call add_edge(u+1,v+1)
      end do
    else
      do i=1,nx
        do j=1,ny
          u=(i-1)*ny+j
          if(i<nx) call add_edge(u,u+ny)
          if(j<ny) call add_edge(u,u+1)
        end do
      end do
      do j=1,ny
        call boundary(j,param_diff%bcl,param_diff%pbcl,param_diff%qbcl)
        q=param_diff%qbcr
        if(.not.two_dimensional) q=-q ! 1-D qbc is positive in the +coordinate direction.
        call boundary((nx-1)*ny+j,param_diff%bcr,param_diff%pbcr,q)
      end do
      if(two_dimensional) then
        do i=1,nx
          call boundary((i-1)*ny+1,param_diff%bct,param_diff%pbct,param_diff%qbct)
          call boundary(i*ny,param_diff%bcb,param_diff%pbcb,param_diff%qbcb)
        end do
      end if
    end if
    if(param_diff%injectionfromfile) then
      do kw=1,param_diff%nwell
        i=param_diff%iwell(kw); j=1
        if(two_dimensional) j=param_diff%jwell(kw)
        if(i<1.or.i>nx.or.j<1.or.j>ny) error stop 'Well index outside diffusion grid'
        u=(i-1)*ny+j
        rhs(u)=rhs(u)+h*well_rate(param_diff,kw,time)/(area*1d6*scale)
      end do
    else if(time<param_diff%tinj*31536000d0) then
      u=param_diff%ninj
      if(two_dimensional) u=(max(1,nx/2)-1)*ny+max(1,ny/2)
      select case(trim(param_diff%injection))
      case('flowrate','pressure')
        if(u<1.or.u>n) error stop 'Injection element outside diffusion grid'
        if(param_diff%injection=='flowrate') then
          rhs(u)=rhs(u)+h*param_diff%qinj/(area*1d6*scale)
        else
          if(two_dimensional) error stop 'Pressure injection is supported only in 1-D diffusion'
          g=h*mass(u)/(10d0*scale)
          diag(u)=diag(u)+g; rhs(u)=rhs(u)+g*param_diff%pinj
        end if
      case('distributed')
        if(.not.present(dep)) error stop 'Distributed injection requires a fault network'
        rhs=rhs+h*param_diff%qs/(dx*1d6*scale) ! Per-element raw rate per unit strike.
      end select
    end if
    pfnew=pf
    call multiply(pfnew,ap)
    r=rhs-ap
    tol=1d-12*max(sqrt(sum(rhs*rhs)),sqrt(sum(ap*ap)),1d-12)
    if(sqrt(sum(r*r))<=tol) return
    z=r/diag; p=z; rz=sum(r*z)
    itermax=max(200,10*n)
    do iteration=1,itermax
      call multiply(p,ap)
      denom=sum(p*ap)
      if(denom<=0d0.or..not.ieee_is_finite(denom)) error stop 'Invalid diffusion PCG matrix'
      alpha=rz/denom; pfnew=pfnew+alpha*p; r=r-alpha*ap
      niter=iteration
      if(sqrt(sum(r*r))<=tol) exit
      z=r/diag; rznew=sum(r*z); p=z+(rznew/rz)*p; rz=rznew
    end do
    if(iteration>itermax) error stop 'Conservative diffusion PCG failed to converge'
    if(any(.not.ieee_is_finite(pfnew))) error stop 'Nonfinite diffusion pressure'
  contains
    subroutine add_edge(a,b)
      integer,intent(in)::a,b
      real(8)::t,small,large
      if(a==b) return
      small=min(mobility(a),mobility(b)); large=max(mobility(a),mobility(b))
      t=0d0
      if(large>0d0) t=2d0*small/(1d0+small/large)
      ne=ne+1; left(ne)=a; right(ne)=b; edge(ne)=h*t/(dx*dx*scale)
      diag(a)=diag(a)+edge(ne); diag(b)=diag(b)+edge(ne)
    end subroutine
    subroutine boundary(a,kind,pressure,flux)
      integer,intent(in)::a
      character(*),intent(in)::kind
      real(8),intent(in)::pressure,flux
      real(8)::c
      select case(trim(kind))
      case('Dirichlet')
        c=h*2d0*mobility(a)/(dx*dx*scale)
        diag(a)=diag(a)+c; rhs(a)=rhs(a)+c*pressure
      case('Neumann')
        ! Boundary values remain Darcy flux [m/s]; multiply by hydraulic width.
        rhs(a)=rhs(a)+h*flux*width(a)/(dx*1d6*scale)
      case default
        error stop 'Unknown diffusion boundary type'
      end select
    end subroutine
    subroutine multiply(x,y)
      real(8),intent(in)::x(:)
      real(8),intent(out)::y(:)
      integer::e,a,b
      y=diag*x
      do e=1,ne
        a=left(e); b=right(e)
        y(a)=y(a)-edge(e)*x(b); y(b)=y(b)-edge(e)*x(a)
      end do
    end subroutine
  end subroutine

  subroutine Beuler1d(ncell,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
    implicit none
    integer,intent(in)::ncell
    real(8),intent(in)::ds0,pf(:),cdiff(:),str(:),h,time
    real(8),intent(out)::pfnew(:)
    integer,intent(out)::niter
    type(t_params),intent(in)::param_diff
    call conservative_diffusion(ncell,1,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
  end subroutine

  subroutine Beuler1dn(ncell,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter,dep)
    implicit none
    integer,intent(in)::ncell
    real(8),intent(in)::ds0,pf(:),cdiff(:),str(:),h,time,dep(:)
    real(8),intent(out)::pfnew(:)
    integer,intent(out)::niter
    type(t_params),intent(in)::param_diff
    call conservative_diffusion(ncell,1,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter,dep)
  end subroutine

  subroutine Beuler2d(imax,jmax,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
    implicit none
    integer,intent(in)::imax,jmax
    real(8),intent(in)::ds0,pf(:,:),cdiff(:,:),str(:,:),h,time
    real(8),intent(out)::pfnew(:,:)
    integer,intent(out)::niter
    type(t_params),intent(in)::param_diff
    real(8)::p(imax*jmax),d(imax*jmax),s(imax*jmax),out(imax*jmax)
    integer::i,j,k
    do i=1,imax
      do j=1,jmax
        k=(i-1)*jmax+j; p(k)=pf(i,j); d(k)=cdiff(i,j); s(k)=str(i,j)
      end do
    end do
    call conservative_diffusion(imax,jmax,ds0,p,d,s,param_diff,h,out,time,niter,surface=.true.)
    do i=1,imax
      do j=1,jmax
        pfnew(i,j)=out((i-1)*jmax+j)
      end do
    end do
  end subroutine

subroutine setup_network(number,param_diff,my_rank)
  implicit none
  integer,intent(in)::number
  type(t_params):: param_diff
  integer,intent(in)::my_rank
  integer::k,ios
  integer,allocatable::ncell1(:)
  if(param_diff%network_file=='default') then
    write(param_diff%network_file,'("examples/connection",i0,".dat")') number
  end if
  write(*,*) param_diff%network_file
  open(77,file=param_diff%network_file,iostat=ios)
  if(ios /= 0) then
    if(my_rank==0)write(*,*) 'ERROR: Failed to open network file'
    stop
  end if
  read(77,*) param_diff%nfault
  allocate(ncell1(param_diff%nfault),param_diff%ns(param_diff%nfault+1))
  read(77,*) ncell1
  param_diff%ns(1)=0
  do k=1,param_diff%nfault
    param_diff%ns(k+1)=param_diff%ns(k)+ncell1(k)
  end do
  read(77,*) param_diff%nconnect
  write(*,*) param_diff%ns
  param_diff%ncell=param_diff%ns(param_diff%nfault+1)

  allocate(param_diff%connect(param_diff%nconnect,2))
  do k=1,param_diff%nconnect
    read(77,*) param_diff%connect(k,:)
  end do
  write(*,*) param_diff%connect

  close(77)
  return
end subroutine
subroutine input_well(problem,param_diff,my_rank)
  implicit none
  character(*),intent(in)::problem
  type(t_params):: param_diff
  integer,intent(in)::my_rank
  integer::k,kwell,ios
  open(77,file=param_diff%injection_file,iostat=ios)
  if(ios /= 0) then
    if(my_rank==0)write(*,*) 'ERROR: Failed to open well file'
    stop
  end if
  read(77,*) param_diff%nwell
  read(77,*) param_diff%npoint
  if(param_diff%nwell<1.or.param_diff%npoint<2) error stop 'Invalid well history dimensions'
  allocate(param_diff%iwell(param_diff%nwell),param_diff%jwell(param_diff%nwell),param_diff%qvals(param_diff%nwell,param_diff%npoint))
  allocate(param_diff%qtimes(param_diff%npoint))
  read(77,*) param_diff%qtimes(1:param_diff%npoint)
  if(any(.not.ieee_is_finite(param_diff%qtimes))) error stop 'Nonfinite well times'
  if(any(param_diff%qtimes(2:)<param_diff%qtimes(:param_diff%npoint-1))) &
    error stop 'Well times must not decrease'
  select case(problem)
  case('3dp', '3dph')
    do kwell=1,param_diff%nwell
      read(77,*) param_diff%iwell(kwell),param_diff%jwell(kwell)
    end do
  case('2dp','2dn','2dvs','2dph','2dnh')
    do kwell=1,param_diff%nwell
      read(77,*) param_diff%iwell(kwell)
    end do
  end select
  do kwell=1,param_diff%nwell
    read(77,*) param_diff%qvals(kwell,1:param_diff%npoint)
  end do
  if(any(.not.ieee_is_finite(param_diff%qvals))) error stop 'Nonfinite well rates'
  
  close(77)
  param_diff%switch=.false.
  param_diff%nn=2
end subroutine 
     !pressure dependent permeability
  subroutine diffusion2dwp(pf,sigma,h,ds0,time,dtnxt,param_diff)
    implicit none
    integer::kit,errloc(1),i,n,niter
    real(8),intent(inout)::pf(:),sigma(:),dtnxt
    real(8),intent(in)::h,time,ds0
    real(8),dimension(size(pf))::dpf,pftry,pfnew,sigmae,cdiff,str,pfinit
    real(8)::err,err0,cc
    real(8),parameter::dpth=0.2
    type(t_params):: param_diff
    pfinit=pf
    pftry=pf
    n=size(pf)

    err=0d0
    do kit=1,50
      !calculate diffusion coefficient
      cdiff=param_diff%kpG*exp(-sigma/param_diff%sigmastar)*(1+param_diff%sigmaHF/sigma)/(param_diff%eta*param_diff%beta*param_diff%phiG)*1d-6
      str=param_diff%beta*param_diff%phiG
      !write(*,*) cdiff(1),str(1)
      if(param_diff%network) then
        call Beuler1dn(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter,param_diff%dep)
      else
        call Beuler1d(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
      end if
      err=maxval(abs(pftry-pfnew))
      !write(*,*) pftry
      !write(*,*) pfnew
      errloc=maxloc(abs(pftry-pfnew))
      !write(*,*) 'err',err
      if(err<1e-3) exit
      sigma=sigma+pftry-pfnew
      pftry=pfnew
      err0=err
    end do

    do i=1,size(pf)
      dpf(i)=pfnew(i)-pf(i)
      pf(i)=pfnew(i)
      !y(4*i-2)=y(4*i-2)-mu(i)*dpf(i)
    end do

    if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)
    return
  end subroutine

  ! pressure-independent permeability (ks=kp)
  subroutine diffusion2dwop(pf,h,ds0,time,dtnxt,param_diff)
    implicit none
    integer::kit,errloc(1),i,n,niter
    real(8),intent(inout)::pf(:),dtnxt
    real(8),intent(in)::h,time,ds0
    real(8),dimension(size(pf))::dpf,pftry,pfnew,sigmae,sigma,cdiff,str
    real(8)::err,err0,cc
    real(8),parameter::dpth=0.2,tny=1d0
    type(t_params):: param_diff

    n=size(pf)

    !calculate diffusion coefficient
    cdiff=param_diff%kpG/(param_diff%eta*param_diff%beta*param_diff%phiG)*1d-6
    str=param_diff%beta*param_diff%phiG
    !write(*,*) cdiff(1),str(1)
    if(param_diff%network) then
      call Beuler1dn(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter,param_diff%dep)
    else
      call Beuler1d(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
    end if

    !$omp parallel do
    do i=1,size(pf)
      dpf(i)=pfnew(i)-pf(i)
      pf(i)=pfnew(i)
    end do
    !$omp end parallel do

    !write(*,*) 'niter',niter
    !write(*,*) 'dpf',maxval(dpf)
    if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)

     !check if next time step is after the change in the injection rate
    call diffusion_next_knot(param_diff,time,dtnxt)

    return
  end subroutine

  ! pressure-independent permeability (ks=kp)
  subroutine diffusionpseudo2dwop(pf,h,ds0,time,dtnxt,param_diff)
    implicit none
    integer::kit,errloc(1),i,n,niter
    real(8),intent(inout)::pf(:),dtnxt
    real(8),intent(in)::h,time,ds0
    real(8),dimension(size(pf))::dpf,pftry,pfnew,sigmae,sigma,cdiff,str
    real(8)::err,err0,cc
    real(8),parameter::dpth=0.2,tny=1d0
    type(t_params):: param_diff

    n=size(pf)
    !calculate diffusion coefficient
    cdiff=param_diff%kpG/(param_diff%eta*param_diff%beta*param_diff%phiG)*1d-6
    str=param_diff%beta*param_diff%phiG
    !write(*,*) cdiff(1),str(1)
    if(param_diff%network) then
      call Beuler1dn(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter,param_diff%dep)
    else
      call Beuler1d(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
    end if

    !$omp parallel do
    do i=1,size(pf)
      dpf(i)=pfnew(i)-pf(i)
      pf(i)=pfnew(i)
    end do
    !$omp end parallel do

    !write(*,*) 'niter',niter
    !write(*,*) 'dpf',maxval(dpf)
    if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)

     !check if next time step is after the change in the injection rate
    call diffusion_next_knot(param_diff,time,dtnxt)

    return
  end subroutine


  subroutine diffusion3dwop(imax,jmax,pf,h,ds0,time,dtnxt,param_diff)
    implicit none
    integer::kit,errloc(1),l,l_,i,j,niter,nn
    integer,intent(in)::imax,jmax
    real(8),intent(inout)::pf(:),dtnxt
    real(8),intent(in)::h,time,ds0
    real(8)::dpf(imax*jmax),pfd(imax,jmax),pfnew(imax,jmax),sigmae(imax*jmax),cdiff(imax,jmax),err,err0,adpf(imax*jmax),pfhydd(imax,jmax)
    real(8)::cc,str(imax,jmax),x
    real(8),parameter::dpth=0.1,tny=1d0
    type(t_params):: param_diff

    !real(8),parameter::cc=1d-12 !beta(1e-8)*phi(1e-1)*eta(1e-3)
    cdiff=0d0
    
    do l=1,size(pf)
      !sigma(l)=sigmae(l)+pf(l)
      i=(l-1)/jmax+1
      j=l-(i-1)*jmax
      pfd(i,j)=pf(l)
      !pfhydd(i,j)=pfhyd(l)
      cc=param_diff%eta*param_diff%beta*param_diff%phiG(l)
      str(i,j)=param_diff%beta*param_diff%phiG(l)
      cdiff(i,j)=param_diff%kpG(l)/cc*1d-6 !Pa-->MPa
      !cdiff(i,j)=kpmax/cc*1d-6 !Pa-->MPa
      !write(*,*) l_,i,j
    end do

    err=0d0

    call Beuler2d(imax,jmax,ds0,pfd,cdiff,str,param_diff,h,pfnew,time,niter)!,pfhydd)

    do l=1,imax*jmax
      i=(l-1)/jmax+1
      j=l-(i-1)*jmax
      dpf(l)=pfnew(i,j)-pf(l)
      pf(l)=pfnew(i,j)
    end do
      !write(*,*)sum(pf)

    !write(*,*) param_diff%setting
      adpf=abs(dpf)
    !write(*,*) 'dpf',maxval(adpf)
    if(dtnxt/h*maxval(adpf)>dpth)  dtnxt=dpth*h/maxval(adpf)

    !check if next time step is after the change in the injection rate
    call diffusion_next_knot(param_diff,time,dtnxt)
    return
  end subroutine

    subroutine initialize_steadystate(param_diff)
    type(t_params):: param_diff
    end subroutine

end module
