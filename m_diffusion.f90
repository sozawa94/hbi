module mod_diffusion
use mod_constant
  type :: t_params
  integer::nwell,nn
  integer,pointer::kleng(:),iwell(:),jwell(:)
  real(8)::phi0,beta,eta,s0,kp0,kpmin,kpmax,kL,kT,pinj,pbcl,pbcr,pbct,pbcb,qinj,q0
  real(8)::qbcl,qbcr,qbct,qbcb
  real(8)::tinj=1d5
  real(8),pointer::kp(:),kpG(:),qtimes(:,:),qvals(:,:),pfhyd(:,:),phi(:),phiG(:)
  character(128)::bc,bcl,bcr,bct,bcb,setting,injection,injection_file
  logical::injectionfromfile,switch,permev
  end type t_params
contains
subroutine input_well(param_diff)
  implicit none
  type(t_params):: param_diff
  integer::k,kwell
  open(77,file=param_diff%injection_file)
  write(*,*) "Read injection file"
  read(77,*) param_diff%nwell
  allocate(param_diff%iwell(param_diff%nwell),param_diff%jwell(param_diff%nwell),param_diff%qvals(param_diff%nwell,50))
  allocate(param_diff%qtimes(param_diff%nwell,50),param_diff%kleng(param_diff%nwell))
  do kwell=1,param_diff%nwell
    read(77,*) param_diff%iwell(kwell),param_diff%jwell(kwell),param_diff%kleng(kwell)
    !write(*,*) iwell(kwell),jwell(kwell)
    read(77,*) param_diff%qvals(kwell,1:param_diff%kleng(kwell))
    read(77,*) param_diff%qtimes(kwell,1:param_diff%kleng(kwell))
  end do
  close(77)
  param_diff%switch=.false.
  param_diff%nn=2
end subroutine 
     !pressure dependent permeability
  ! subroutine implicitsolver(pf,sigma,ks,phi,beta,eta,s0,h,time,dtnxt)
  !   implicit none
  !   integer::kit,errloc(1),niter,i,n
  !   real(8),intent(inout)::pf(:),dtnxt
  !   real(8),intent(in)::h,time,sigma(:),ks(:),phi,eta,beta,s0
  !   real(8),dimension(size(pf))::dpf,pftry,pfnew,sigmae,cdiff
  !   real(8)::err,err0,cc
  !   real(8),parameter::dpth=0.3
  !   !real(8),parameter::cc=1d-12 !beta(1e-8)*phi(1e-1)*eta(1e-3) !BP6
  !   !real(8),parameter::cc=1d-15 !beta(1e-9)*phi(1e-2)*eta(1e-4) !Zhu et al.
  !   pftry=pf
  !   n=size(pf)

  !   err=0d0
  !   do kit=1,20
  !     !calculate diffusion coefficient
  !     cc=eta*beta*phi
  !     cdiff=ks*exp(-(sigma-pftry)/s0)/cc*1d-6
  !       !cdiff(i)=ks(i)/cc*1d-6

  !     call Beuler(n,pf,cdiff,h,pfnew,time,niter,bc)
  !     err=maxval(abs(pftry-pfnew))
  !     !write(*,*) pftry
  !     !write(*,*) pfnew
  !     errloc=maxloc(abs(pftry-pfnew))
  !     !write(*,*) 'err',err
  !     if(err<1e-3) exit
  !     pftry=pfnew
  !     err0=err
  !   end do

  !   do i=1,size(pf)
  !     dpf(i)=pfnew(i)-pf(i)
  !     pf(i)=pfnew(i)
  !     !y(4*i-2)=y(4*i-2)-mu(i)*dpf(i)
  !   end do

  !   if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)
  !   return
  ! end subroutine

  ! pressure-independent permeability (ks=kp)
  subroutine diffusion2dwop(pf,h,ds0,time,dtnxt,param_diff)
    implicit none
    integer::kit,errloc(1),i,n,niter
    real(8),intent(inout)::pf(:),dtnxt
    real(8),intent(in)::h,time,ds0
    real(8),dimension(size(pf))::dpf,pftry,pfnew,sigmae,sigma,cdiff,str
    real(8)::err,err0,cc
    real(8),parameter::dpth=0.2
    type(t_params):: param_diff

    !real(8),parameter::cc=1d-12 !beta(1e-8)*phi(1e-1)*eta(1e-3)
    n=size(pf)


    !$omp parallel do
    do i=1,size(pf)
      pftry(i)=pf(i)
      !sigmae(i)=y(4*i-1)
      !sigma(i)=sigmae(i)+pf(i)
    end do
    !$omp end parallel do

    err=0d0

      !calculate diffusion coefficient
    cdiff=param_diff%kpG/(param_diff%eta*param_diff%beta*param_diff%phiG)*1d-6
    str=param_diff%beta*param_diff%phiG
    !write(*,*) cdiff(1),str(1)

    call Beuler1d(n,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)

    !$omp parallel do
    do i=1,size(pf)
      dpf(i)=pfnew(i)-pf(i)
      pf(i)=pfnew(i)
      !y(4*i-1)=y(4*i-1)-dpf(i) !update effective normal stress
    end do
    !$omp end parallel do

    !write(*,*) 'niter',niter
    !write(*,*) 'dpf',maxval(dpf)
    if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)
    return
  end subroutine

  subroutine Beuler1d(ncell,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
    integer,parameter::itermax=2000
    integer,intent(in)::ncell
    real(8),intent(in)::pf(:),h,cdiff(:),str(:),time,ds0
    real(8),intent(out)::pfnew(:)
    integer,intent(out)::niter
    real(8)::eta, tmp1,tmp2,rsnew,rsold
    real(8),dimension(ncell)::pf0,m,dpf,p,r,b,x,alpha,sat
    real(8)::Dxx(ncell,3),Am(ncell,3)
    integer::n,iter
    real(8)::p0=0.0,td=1d6,tol=1e-4
    type(t_params):: param_diff
    !real(8),parameter::str=1e-11 !beta(1e-9)*phi(1e-2)
    n=ncell
    niter=0
    Dxx=0d0

    !compute Dxx for Dirichlet BC
    select case(param_diff%bcl)
    case('Dirichlet')
      Dxx(1,1)=-cdiff(1)-cdiff(2)
      Dxx(1,2)=cdiff(2)-cdiff(1)
      Dxx(2,1:3)=(/cdiff(2)/2-cdiff(1)/2, -cdiff(1)/2-cdiff(2)-cdiff(3)/2, cdiff(2)/2+cdiff(3)/2/)
    case('Neumann')
      Dxx(1,1)=-cdiff(1)-cdiff(2)
      Dxx(1,2)=-Dxx(1,1)
      Dxx(2,1:3)=(/cdiff(1)/2+cdiff(2)/2, -cdiff(1)/2-cdiff(2)-cdiff(3)/2, cdiff(2)/2+cdiff(3)/2/)
    end select

    do i=3,n-2
      Dxx(i,1:3)=(/cdiff(i-1)/2+cdiff(i)/2, -cdiff(i-1)/2-cdiff(i)-cdiff(i+1)/2, cdiff(i)/2+cdiff(i+1)/2/)
    end do

    select case(param_diff%bcr)
    case('Dirichlet')
      Dxx(n-1,1:3)=(/cdiff(n-2)/2+cdiff(n-1)/2, -cdiff(n-2)/2-cdiff(n-1)-cdiff(n)/2, -cdiff(n)/2+cdiff(n-1)/2/)
      Dxx(n,3)=-cdiff(n)-cdiff(n-1)
      Dxx(n,2)=cdiff(n-1)-cdiff(n)

    case('Neumann')
      Dxx(n-1,1:3)=(/cdiff(n-2)/2+cdiff(n-1)/2, -cdiff(n-2)/2-cdiff(n-1)-cdiff(n)/2, cdiff(n-1)/2+cdiff(n)/2/)
      Dxx(n,3)=-cdiff(n)-cdiff(n-1)
      Dxx(n,2)=-Dxx(n,3)
    end select

    Dxx=Dxx/ds0/ds0

    !write(*,*)

    Am=0d0
    Am(1,1)=1.0-h*Dxx(1,1)
    Am(1,2)=-h*Dxx(1,2)
    Am(n,3)=1.0-h*Dxx(n,3)
    Am(n,2)=-h*Dxx(n,2)

     !$omp parallel do
    do i=2,n-1
      Am(i,1)=-h*Dxx(i,1)
      Am(i,2)=1.0-h*Dxx(i,2)
      Am(i,3)=-h*Dxx(i,3)
    end do
     !$omp end parallel do

    do i=1,n
      !write(*,'(3e15.6)') Am(i,:)
    end do

    SAT=0d0

    !penalty vector
    select case(param_diff%bcl)
    case('Dirichlet')
      !SAT(1)=-1.0*cdiff(1)/ds0/ds0*(pbcl-pfhyd(1))*h*2
      !SAT(2)=-0.5*cdiff(1)/ds0/ds0*(pbcl-pfhyd(1))*h*2
      SAT(1)=-1.0*cdiff(1)/ds0/ds0*pbcl*h*2
      SAT(2)=-0.5*cdiff(1)/ds0/ds0*pbcl*h*2
    case('Neumann')
      SAT(1)=-param_diff%qbcl/str(1)*1e-9*h/ds0*2
    end select

    select case(param_diff%bcr)
    case('Dirichlet')
      !SAT(n)=-1*cdiff(n)/ds0/ds0*(pbcr-pfhyd(n))*h*2
      SAT(n)=-1*cdiff(n)/ds0/ds0*pbcr*h*2
      SAT(n-1)=SAT(n)/2
    case('Neumann')
      SAT(n)=param_diff%qbcr/str(n)*1e-9*h/ds0*2
      !if(setting=='thrust') SAT(n)=qin(2000)/str(n)*1e-9*h/ds0*2
    end select

    x=pf!-pfhyd !initial guess

    b=pf-SAT!-pfhyd   
    !injection at the center of the fault
    if(param_diff%injection=='pressure' .and. time<param_diff%tinj*365*24*3600) then
      b(N/2)=b(N/2)+h*param_diff%pinj/1e1 !injection pressure
      Am(N/2,2)=Am(N/2,2)+h/1e1
    else if(param_diff%injection=='flowrate' .and. time<param_diff%tinj*365*24*3600) then
      b(N/2)=b(N/2)+h*param_diff%qinj/str(N/2)*1e-9/ds0 !injection rate
    end if
    !write(*,*) h,str(N/2),ds0

    b(1)=b(1)/2
    b(n)=b(n)/2

    Am(1,:)=Am(1,:)/2
    Am(n,:)=Am(n,:)/2

    m=0d0
    m(1)=x(1)*Am(1,1)+x(2)*Am(1,2)

    !$omp parallel do
    do i=2,n-1
      m(i)=x(i-1)*Am(i,1)+x(i)*Am(i,2)+x(i+1)*Am(i,3)
    end do
    !$omp end parallel do

    m(n)=x(n)*Am(n,3)+x(n-1)*Am(n,2)
    !write(*,'(9e15.6)')m

    r=b-m
    p=r
    rsold=sum(r*r)
    ! write(*,*) rsold
    if(rsold<tol**2*n)  then
      go to 100
    end if
    niter=itermax
    do iter=1,itermax
      tmp1=sum(r*r)
      m=0d0
      m(1)=p(1)*Am(1,1)+p(2)*Am(1,2)
      !$omp parallel do
      do i=2,n-1
        m(i)=p(i-1)*Am(i,1)+p(i)*Am(i,2)+p(i+1)*Am(i,3)
      end do
      !$omp end parallel do
      m(n)=p(n)*Am(n,3)+p(n-1)*Am(n,2)
      !write(*,'(9e15.6)')m

      tmp2=sum(m*p)
      alpha=tmp1/tmp2
      x=x+alpha*p
      r=r-alpha*m
      !write(*,'(9e15.6)')r
      rsnew = sum(r*r)
      write(*,*)iter,rsnew
      if(rsnew<tol**2*n) then
          niter=iter
        exit
      end if
      p = r + (rsnew / rsold) * p
      rsold = rsnew
      !write(*,'(9e15.6)')x

    end do

    if(niter==itermax) write(*,*) "Maximum iteration"
    100 pfnew=x!+pfhyd
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
    if(param_diff%switch) then
        dtnxt=2e1
        param_diff%switch=.false.
    end if
    if(param_diff%injectionfromfile) then
      if(time+dtnxt>param_diff%qtimes(1,param_diff%nn)) then
        dtnxt=param_diff%qtimes(1,param_diff%nn)-time-tny
        param_diff%switch=.true.
        param_diff%nn=param_diff%nn+1
        if(param_diff%qtimes(1,param_diff%nn)-param_diff%qtimes(1,param_diff%nn-1)<1e0) param_diff%nn=param_diff%nn+1
        !write(*,*) param_diff%nn,param_diff%qtimes(1,param_diff%nn),time+dtnxt
      end if
    end if
    
    return
  end subroutine

  subroutine Beuler2d(imax,jmax,ds0,pf,cdiff,str,param_diff,h,pfnew,time,niter)
    implicit none
    integer,parameter::itermax=1000
    integer,intent(in)::imax,jmax
    real(8),intent(in)::pf(:,:),h,cdiff(:,:),str(:,:),time,ds0!,pfhyd(:,:)
    real(8),intent(out)::pfnew(:,:)
    integer,intent(out)::niter
    real(8)::Dxx(imax,jmax,3),Dyy(imax,jmax,3),Amx(imax,jmax,3),Amy(imax,jmax,3),mx(imax,jmax),my(imax,jmax)
    real(8)::p(imax,jmax),m(imax,jmax),r(imax,jmax),x(imax,jmax),b(imax,jmax),SAT(imax,jmax)
    integer::n,iter,i,j,k,kwell
    real(8)::p0=0.0,rsold,rsnew,tmp1,tmp2,alpha,v1,v0,t1,t0,qtmp,qdt
    real(8),parameter::tol=1e-6
    type(t_params):: param_diff
    !real(8),parameter::str=1e-11 !beta(1e-9)*phi(1e-2)
    niter=0
    p=0d0;m=0d0;r=0d0;x=0d0;b=0d0

    Dxx=0d0 
    do i=1,imax
        !compute Dxx for Dirichlet BC
        select case(param_diff%bct)
        case('Dirichlet')
        Dxx(i,1,1)=-cdiff(i,1)-cdiff(i,2)
        Dxx(i,1,2)=cdiff(i,2)-cdiff(i,1)
        Dxx(i,2,1:3)=(/cdiff(i,2)/2-cdiff(i,1)/2, -cdiff(i,1)/2-cdiff(i,2)-cdiff(i,3)/2, cdiff(i,2)/2+cdiff(i,3)/2/)
        case('Neumann')
        Dxx(i,1,1)=-cdiff(i,1)-cdiff(i,2)
        Dxx(i,1,2)=-Dxx(i,1,1)
        Dxx(i,2,1:3)=(/cdiff(i,1)/2+cdiff(i,2)/2, -cdiff(i,1)/2-cdiff(i,2)-cdiff(i,3)/2, cdiff(i,2)/2+cdiff(i,3)/2/)
        end select

        do j=3,jmax-2
        Dxx(i,j,1:3)=(/cdiff(i,j-1)/2+cdiff(i,j)/2, -cdiff(i,j-1)/2-cdiff(i,j)-cdiff(i,j+1)/2, cdiff(i,j)/2+cdiff(i,j+1)/2/)
        end do

        select case(param_diff%bcb)
        case('Dirichlet')
        Dxx(i,jmax-1,1:3)=(/cdiff(i,jmax-2)/2+cdiff(i,jmax-1)/2, -cdiff(i,jmax-2)/2-cdiff(i,jmax-1)-cdiff(i,jmax)/2, -cdiff(i,jmax)/2+cdiff(i,jmax-1)/2/)
        Dxx(i,jmax,3)=-cdiff(i,jmax)-cdiff(i,jmax-1)
        Dxx(i,jmax,2)=cdiff(i,jmax-1)-cdiff(i,jmax)
        case('Neumann')
        Dxx(i,jmax-1,1:3)=(/cdiff(i,jmax-2)/2+cdiff(i,jmax-1)/2, -cdiff(i,jmax-2)/2-cdiff(i,jmax-1)-cdiff(i,jmax)/2, cdiff(i,jmax-1)/2+cdiff(i,jmax)/2/)
        Dxx(i,jmax,3)=-cdiff(i,jmax)-cdiff(i,jmax-1)
        Dxx(i,jmax,2)=-Dxx(i,jmax,3)
        end select
    end do

    Dxx=Dxx/ds0/ds0

    Dyy=0d0
    do j=1,jmax
        !compute Dxx for Dirichlet BC
        select case(param_diff%bcl)
        case('Dirichlet')
        Dyy(1,j,1)=-cdiff(1,j)-cdiff(2,j)
        Dyy(1,j,2)=cdiff(2,j)-cdiff(1,j)
        Dyy(2,j,1:3)=(/cdiff(2,j)/2-cdiff(1,j)/2, -cdiff(1,j)/2-cdiff(2,j)-cdiff(3,j)/2, cdiff(2,j)/2+cdiff(3,j)/2/)
        case('Neumann')
        Dyy(1,j,1)=-cdiff(1,j)-cdiff(2,j)
        Dyy(1,j,2)=-Dyy(1,j,1)
        Dyy(2,j,1:3)=(/cdiff(1,j)/2+cdiff(2,j)/2, -cdiff(1,j)/2-cdiff(2,j)-cdiff(3,j)/2, cdiff(2,j)/2+cdiff(3,j)/2/)
        end select

        do i=3,imax-2
        Dyy(i,j,1:3)=(/cdiff(i-1,j)/2+cdiff(i,j)/2, -cdiff(i-1,j)/2-cdiff(i,j)-cdiff(i+1,j)/2, cdiff(i,j)/2+cdiff(i+1,j)/2/)
        end do

        select case(param_diff%bcr)
        case('Dirichlet')
        Dyy(imax-1,j,1:3)=(/cdiff(imax-2,j)/2+cdiff(imax-1,j)/2, -cdiff(imax-2,j)/2-cdiff(imax-1,j)-cdiff(imax,j)/2, -cdiff(imax,j)/2+cdiff(imax-1,j)/2/)
        Dyy(imax,j,3)=-cdiff(imax,j)-cdiff(imax-1,j)
        Dyy(imax,j,2)=cdiff(imax-1,j)-cdiff(imax,j)
        case('Neumann')
        Dyy(imax-1,j,1:3)=(/cdiff(imax-2,j)/2+cdiff(imax-1,j)/2, -cdiff(imax-2,j)/2-cdiff(imax-1,j)-cdiff(imax,j)/2, cdiff(imax-1,j)/2+cdiff(imax,j)/2/)
        Dyy(imax,j,3)=-cdiff(imax,j)-cdiff(imax-1,j)
        Dyy(imax,j,2)=-Dyy(imax,j,3)
        end select
    end do

    Dyy=Dyy/ds0/ds0

    Amx=0d0
    Amx(:,1,1)=0.5-h*Dxx(:,1,1)
    Amx(:,1,2)=-h*Dxx(:,1,2)
    do j=2,jmax-1
      Amx(:,j,1)=-h*Dxx(:,j,1)
      Amx(:,j,2)=0.5-h*Dxx(:,j,2)
      Amx(:,j,3)=-h*Dxx(:,j,3)
    end do
    Amx(:,jmax,3)=0.5-h*Dxx(:,jmax,3)
    Amx(:,jmax,2)=-h*Dxx(:,jmax,2)

    Amy=0d0
    Amy(1,:,1)=0.5-h*Dyy(1,:,1)
    Amy(1,:,2)=-h*Dyy(1,:,2)
    do i=2,imax-1
      Amy(i,:,1)=-h*Dyy(i,:,1)
      Amy(i,:,2)=0.5-h*Dyy(i,:,2)
      Amy(i,:,3)=-h*Dyy(i,:,3)
    end do
    Amy(imax,:,3)=0.5-h*Dyy(imax,:,3)
    Amy(imax,:,2)=-h*Dyy(imax,:,2)

    SAT=0d0
  
    x=pf!-pfhyd !initial guess

    select case(param_diff%bct)
      case('Neumann')
        SAT(:,1)=-param_diff%qbct/str(:,jmax)*1e-9*h/ds0*2
      case('Dirichlet')
        SAT(:,1)=-1.0*cdiff(:,jmax)/ds0/ds0*param_diff%pbct*h*2
        SAT(:,2)=-0.5*cdiff(:,jmax)/ds0/ds0*param_diff%pbct*h*2
    end select

    select case(param_diff%bcb)
      case('Neumann')
        SAT(:,jmax)=-param_diff%qbcb/str(:,jmax)*1e-9*h/ds0*2
      case('Dirichlet')
        SAT(:,jmax)=-1.0*cdiff(:,jmax)/ds0/ds0*param_diff%pbcb*h*2
        SAT(:,jmax-1)=-0.5*cdiff(:,jmax)/ds0/ds0*param_diff%pbcb*h*2
    end select

    select case(param_diff%bcl)
      case('Neumann')
        SAT(1,:)=-param_diff%qbcl/str(1,:)*1e-9*h/ds0*2
      case('Dirichlet')
        SAT(1,:)=-1.0*cdiff(1,:)/ds0/ds0*param_diff%pbcl*h*2
        SAT(2,:)=-0.5*cdiff(1,:)/ds0/ds0*param_diff%pbcl*h*2
    end select

    select case(param_diff%bcr)
      case('Neumann')
        SAT(imax,:)=-param_diff%qbcr/str(imax,:)*1e-9*h/ds0*2
      case('Dirichlet')
        SAT(imax,:)=-1.0*cdiff(imax,:)/ds0/ds0*param_diff%pbcr*h*2
        SAT(imax-1,:)=-0.5*cdiff(imax,:)/ds0/ds0*param_diff%pbcr*h*2
    end select

    b=pf-SAT!-pfhyd
    !write(*,*) param_diff%setting
    if(param_diff%injectionfromfile) then
      qtmp=0d0
      do kwell=1,param_diff%nwell
        do k=1,param_diff%kleng(kwell)-1
          t0=param_diff%qtimes(kwell,k)
          t1=param_diff%qtimes(kwell,k+1)
          v0=param_diff%qvals(kwell,k)
          v1=param_diff%qvals(kwell,k+1)
          if (time >= t0 .and. time <= t1) then
            qtmp=(v1-v0)/(t1-t0)*(time-t0)+v0
          else if (time> t1.and. k == param_diff%kleng(kwell)-1) then
            qtmp=v1
          end if
        end do
        if(qtmp<0) qtmp=0d0
        i=param_diff%iwell(kwell)
        j=param_diff%jwell(kwell)
        !write(*,*)i,j,qtmp
        b(i,j)=b(i,j)+h*qtmp/str(i,j)*1e-12/ds0/ds0
      end do
    else if(param_diff%injection=='flowrate' .and. time<param_diff%tinj*365*24*3600) then
    !write(*,*) param_diff%tinj,"injection"
      b(imax/2,jmax/2)=b(imax/2,jmax/2)+h*param_diff%qinj/str(imax/2,jmax/2)*1e-12/ds0/ds0
    ! end if
    end if

    !b(imax/2-5:imax/2+5,jmax/2-5:jmax/2+5)=b(imax/2-5:imax/2+5,jmax/2-5:jmax/2+5)+h*q0/beta/phi0*1e-12/ds0/ds0
    

    mx=0d0
    do i=1,imax
        mx(i,1)=x(i,1)*Amx(i,1,1)+x(i,2)*Amx(i,1,2)
        do j=2,jmax-1
            mx(i,j)=x(i,j-1)*Amx(i,j,1)+x(i,j)*Amx(i,j,2)+x(i,j+1)*Amx(i,j,3)
        end do
        mx(i,jmax)=x(i,jmax)*Amx(i,jmax,3)+x(i,jmax-1)*Amx(i,jmax,2)
    end do

    my=0d0
    do j=1,jmax
        my(1,j)=x(1,j)*Amy(1,j,1)+x(2,j)*Amy(1,j,2)
        do i=2,imax-1
            my(i,j)=x(i-1,j)*Amy(i,j,1)+x(i,j)*Amy(i,j,2)+x(i+1,j)*Amy(i,j,3)
        end do
        my(imax,j)=x(imax,j)*Amy(imax,j,3)+x(imax-1,j)*Amy(imax,j,2)
    end do

    m=mx+my
    r=b-m
    p=r
    rsold=sum(r*r)
    !write(*,*)rsold
    if(rsold<tol**2*imax*jmax)  then
      go to 100
    end if

    niter=itermax
    do iter=1,itermax
        tmp1=sum(r*r)

        mx=0d0
        do i=1,imax
            mx(i,1)=p(i,1)*Amx(i,1,1)+p(i,2)*Amx(i,1,2)
            do j=2,jmax-1
                mx(i,j)=p(i,j-1)*Amx(i,j,1)+p(i,j)*Amx(i,j,2)+p(i,j+1)*Amx(i,j,3)
            end do
            mx(i,jmax)=p(i,jmax)*Amx(i,jmax,3)+p(i,jmax-1)*Amx(i,jmax,2)
        end do
    
        my=0d0
        do j=1,jmax
            my(1,j)=p(1,j)*Amy(1,j,1)+p(2,j)*Amy(1,j,2)
            do i=2,imax-1
                my(i,j)=p(i-1,j)*Amy(i,j,1)+p(i,j)*Amy(i,j,2)+p(i+1,j)*Amy(i,j,3)
            end do
            my(imax,j)=p(imax,j)*Amy(imax,j,3)+p(imax-1,j)*Amy(imax,j,2)
        end do

        m=mx+my
  
        tmp2=sum(m*p)
        alpha=tmp1/tmp2
        x=x+alpha*p
        r=r-alpha*m
        rsnew = sum(r*r)
        !write(*,*)iter,rsnew
        if(rsnew<tol**2*imax*jmax) then
            niter=iter
          exit
        end if
        p = r + (rsnew / rsold) * p
        rsold = rsnew
        !write(*,'(9e15.6)')x
  
      end do

      if(niter==itermax) write(*,*) "Maximum iteration"
      100 pfnew=x!+pfhyd

    return
    end subroutine

end module