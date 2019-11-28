module mod_derivs
  implicit none
  type::t_deriv
    real(8),pointer::a(:),b(:),dc(:),vc(:),sigma(:)
    real(8)::rigid,vs,pois,mu0,vref,vw,fw,vpl,sr
    integer::load
    character(128)::law
  end type

contains
  subroutine derivs(x, y, dydx,&
    &NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    type(st_HACApK_lcontrol),intent(in) :: st_ctl
    type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    type(st_HACApK_calc_entry) :: st_bemv
    type(t_deriv),intent(in) :: st_deriv
    integer,intent(in) :: NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in) :: x
    real(8),intent(in) ::y(:)
    real(8),intent(out) :: dydx(:)
    real(8) :: veltmp(NCELL),tautmp(NCELL),veltmpg(NCELLg)
    real(8) :: dlnvdt(NCELL),dtaudt(NCELL)
    real(8) :: sum_gs(NCELL),sum_gsg(NCELLg)
    real(8):: c1, c2, c3, arg,c,g,tauss,dphidt
    integer :: i, j, nc,ierr,lrtrn

    !if(my_rank.eq.0) then
    do i = 1, NCELL
      veltmp(i) = exp (y(2*i-1))
      tautmp(i) = y(2*i)
    enddo
    !end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

    !matrix-vector mutiplation
    ! st_bemv%v='s'
    !veltmp=1.d-6
    !write(*,*) st_deriv%load
    select case(st_deriv%load)
      case(0)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,sum_gsG,veltmpG)
      case(1)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,sum_gsG,veltmpG-st_deriv%vpl)
    end select
  !write(*,*) sum_gsg(0),sum_gsg(NCELLG)
  !stop
    !call MPI_SCATTERv(sum_gsG,NCELL,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    ! call MPI_SCATTER(sum_gnG,NCELL,MPI_REAL8,sum_gn,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    !do i=1,NCELL
      ! dsigdt(i)=sum_gn(i)+sigdot(i)
      !dsigdt(i)=sum_gn(i)
      do i=1,NCELL
!        sum_gs(i)=sum_gs(i)+st_deriv%taudot(i)
       sum_gs(i)=sum_gs(i)+st_deriv%sr
      end do

      select case(st_deriv%law)
        ! aging law (Linker & Dieterich, 1992)
      case('a')
        call deriv_a(sum_gs,veltmp,tautmp,st_deriv%a,st_deriv%b,st_deriv%dc,&
        &st_deriv%vc,st_deriv%mu0,st_deriv%sigma,st_deriv%vref,st_deriv%rigid,st_deriv%vs,&
        &dlnvdt,dtaudt)

        ! slip law
      case('s')
        call deriv_s(sum_gs,veltmp,tautmp,st_deriv%a,st_deriv%b,st_deriv%dc,&
        &st_deriv%vc,st_deriv%mu0,st_deriv%sigma,st_deriv%vref,st_deriv%rigid,st_deriv%vs,&
        &dlnvdt,dtaudt)

      !case('d') !regularized slip law + dynamic weakening (Dunham+ 2011)
      !  call deriv_d(sum_gs,veltmp,tautmp,st_deriv%a,st_deriv%b,st_deriv%dc,&
      !  &st_deriv%vc,st_deriv%mu0,st_deriv%sigma,st_deriv%vref,st_deriv%rigid,st_deriv%vs,&
      !  &st_deriv%fw,st_deriv%vw,dlnvdt,dtaudt)

       end select

    do i = 1, NCELL
      dydx(2*i-1) = dlnvdt( i )
      dydx(2*i) = dtaudt( i )
      ! dydx(3*i) = dsigdt( i )
    enddo

    return
  end subroutine

  subroutine deriv_a(sum_gs,veltmp,tautmp,a,b,dc,vc,mu0,sigma,vref,rigid,vs,dlnvdt,dtaudt)
    implicit none
    integer::i
    real(8)::arg
    real(8),intent(in)::sum_gs(:),veltmp(:),tautmp(:),a(:),b(:),dc(:),vc(:),sigma(:)
    real(8),intent(in)::mu0,vref,vs,rigid
    real(8),intent(out)::dlnvdt(:),dtaudt(:)
    do i=1,size(sum_gs)
      arg=dc(i)/vref*(exp((tautmp(i)/sigma(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
      dlnvdt(i)=sum_gs(i)-b(i)*sigma(i)/(arg+dc(i)/vc(i))*(1.d0-veltmp(i)*arg/dc(i))
      dlnvdt(i)=dlnvdt(i)/(a(i)*sigma(i)+0.5d0*rigid*veltmp(i)/vs)
      dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
    end do
  end subroutine

  subroutine deriv_s(sum_gs,veltmp,tautmp,a,b,dc,vc,mu0,sigma,vref,rigid,vs,dlnvdt,dtaudt)
    implicit none
    integer::i
    real(8)::arg
    real(8),intent(in)::sum_gs(:),veltmp(:),tautmp(:),a(:),b(:),dc(:),vc(:),sigma(:)
    real(8),intent(in)::mu0,vref,vs,rigid
    real(8),intent(out)::dlnvdt(:),dtaudt(:)
    do i=1,size(sum_gs)
      arg=dc(i)/vref*(exp((tautmp(i)/sigma(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
      !c1 = veltmp(i) / dc(i) * arg
      dlnvdt(i)=sum_gs(i)+b(i)*sigma(i)*veltmp(i)*arg/(arg*dc(i)+dc(i)**2/vc(i))*dlog(veltmp(i)*arg/dc(i))
      dlnvdt(i)=dlnvdt(i)/(a(i)*sigma(i)+0.5d0*rigid*veltmp(i)/vs)
      dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
    end do
  end subroutine

  ! subroutine deriv_d(sum_gs,veltmp,tautmp,a,b,dc,vc,mu0,sigma,rigid,vref,vs,fw,vw,dlnvdt,dtaudt)
  !   implicit none
  !   integer::i
  !   real(8)::theta,g,tauss,dphidt
  !   real(8),intent(in)::sum_gs(:),veltmp(:),tautmp(:),a(:),b(:),dc(:),vc(:),sigma(:)
  !   real(8),intent(in)::mu0,vref,vs,fw,vw,rigid
  !   real(8),intent(out)::dlnvdt(:),dtaudt(:)
  !   do i=1,size(sum_gs)
  !     theta=tautmp(i)/sigma(i)-a(i)*dlog(veltmp(i)/vref)
  !     g=0.5d0*veltmp(i)/vref*exp(theta/a(i))
  !     tauss=mu0+(a(i)-b(i))*dlog(veltmp(i)/vref)
  !     tauss=fw+(tauss-fw)/(1.d0+(veltmp(i)/vw)**8)**0.125d0 !flash heating
  !     dphidt=-veltmp(i)/dc(i)*(tautmp(i)-tauss)
  !     dlnvdt(i)=sum_gs(i)-g/sqrt(1.d0+g**2)*dphidt
  !     dlnvdt(i)=dlnvdt(i)/(a(i)*g/sqrt(1.d0+g**2)+0.5d0*rigid*veltmp(i)/vs)
  !     dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
  !   end do
  ! end subroutine
  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    integer::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::dydx(:),yscal(:),htry,eps
    real(8),intent(inout)::y(:),x
    real(8),intent(out)::hdid,hnext
    type(st_HACApK_lcontrol),intent(in) :: st_ctl
    type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    type(st_HACApK_calc_entry) :: st_bemv
    type(t_deriv),intent(in) :: st_deriv
    integer :: i,ierr
    real(8) :: errmax,h,xnew,htemp,errmax_gb
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89e-4,hmax=1d5

    h=htry
    do while(.true.)
      call rkck(y,dydx,x,h,ytemp,yerr,NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)!,derivs)
      errmax=maxval(abs(yerr(:)/yscal(:)))/eps

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,                  &
      &     MPI_MAX,MPI_COMM_WORLD,ierr)

      if(errmax_gb.lt.1.d0) then
        exit
      end if

      htemp=SAFETY*h*(errmax_gb**PSHRNK)
      h=sign(max(abs(htemp),0.1*abs(h)),h)
      xnew=x+h
      if(xnew-x<1.d-8) stop
    end do

    hnext=min(3*h,SAFETY*h*(errmax_gb**PGROW),5d6)

    hdid=h
    x=x+h
    y(:)=ytemp(:)
    return
  end subroutine

  !---------------------------------------------------------------------
  subroutine rkck(y,dydx,x,h,yout,yerr,NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    integer,intent(in)::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::y(:),dydx(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    type(st_HACApK_lcontrol),intent(in) :: st_ctl
    type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    type(st_HACApK_calc_entry) :: st_bemv
    type(t_deriv),intent(in) :: st_deriv
    integer ::i
    integer,parameter::nmax=100000
    real(8) :: ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax),ytemp(nmax)
    real(8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51
    real(8) :: B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3./40.)
    parameter (B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5)
    parameter (B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.)
    parameter (B63=575./13824.,B64=44275./110592.,B65=253./4096.)
    parameter (C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.)
    parameter (DC1=C1-2825./27648.,DC3=C3-18575./48384.)
    parameter (DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)

    !     -- 1st step --
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+B21*h*dydx(i)
    end do

    !    -- 2nd step --
    call derivs(x+a2*h, ytemp, ak2,&
    & NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
    end do

    !     -- 3rd step --
    call derivs(x+a3*h, ytemp, ak3,&
    & NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
    end do

    !     -- 4th step --
    call derivs(x+a4*h, ytemp, ak4,&
    & NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
    end do

    !     -- 5th step --
    call derivs(x+a5*h, ytemp, ak5,&
    & NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do

    !     -- 6th step --
    call derivs(x+a6*h, ytemp, ak6,&
    & NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
    end do

    !$omp parallel do
    do i=1,size(y)
      yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
    end do
    return
  end subroutine

  subroutine coordinate(imax,jmax,ds,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
    implicit none
    integer,intent(in)::imax,jmax
    real(8),intent(in)::ds
    real(8),intent(out)::xcol(:),zcol(:)
    real(8),intent(out)::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
    integer::i,j,k
    do i=1,imax
      do j=1,jmax
        k=(i-1)*jmax+j
        xcol(k)=(i-0.5d0)*ds
        zcol(k)=(j-0.5d0)*ds
        xs1(k)=xcol(k)+0.5d0*ds
        xs2(k)=xcol(k)-0.5d0*ds
        xs3(k)=xcol(k)-0.5d0*ds
        xs4(k)=xcol(k)+0.5d0*ds
        zs1(k)=zcol(k)+0.5d0*ds
        zs2(k)=zcol(k)+0.5d0*ds
        zs3(k)=zcol(k)-0.5d0*ds
        zs4(k)=zcol(k)-0.5d0*ds
      end do
    end do

  end subroutine coordinate
  subroutine params(NCELLg,a0,b0,dc0,vc0,ag,bg,dcg,vcg)
    integer,intent(in)::NCELLg
    real(8),intent(in)::a0,b0,dc0,vc0
    real(8),intent(out)::ag(:),bg(:),dcg(:),vcg(:)
    integer::i
    !real(8),allocatable::xc(:),zc(:)
    !yc=0.5d0*ds;zc=0.5d0*ds;dr=0.1d0*ds
    ! allocate(xc(kmax),zc(kmax))
    ! do k=1,kmax
    !   100    call random_number(r)
    !   xc(k)=ds*r(1)*imax
    !   zc(k)=ds*r(2)*jmax
    !   lad=dr*3.d0
    !   do i=1,k-1
    !     if((xc(k)-xc(i))**2+(zc(k)-zc(i))**2.lt.(lad)**2) go to 100
    !   end do
    !   write(*,*) xc(k),zc(k)
    ! end do
    !
    ! do i=1,NCELLg
    !   vmax(i)=0.d0
    !   vcg(i)=vc0
    !   do k=1,kmax
    !     dx=xcol(i)-xc(k);dz=zcol(i)-zc(k)
    !     !if(dx**2+dz**2.lt.(1.2*dr)**2) vcg(i)=vc0*1d3
    !     !if(dx**2+dz**2.lt.dr**2) vcg(i)=vc0*1d6
    !     if(dx**2+dz**2.lt.dr**2) vcg(i)=vc0*1d-6
    !   end do
    ! end do

    do i=1,NCELLg
      ag(i)=a0
      !do k=1,kmax
      !dx=xcoll(i)-xc(k);dz=zcoll(i)-zc(k)
      !if(dx**2+dz**2.lt.(1.2*dr)**2) vcg(i)=vc0*1d3
      !if(dx**2+dz**2.lt.dr**2) vcg(i)=vc0*1d6
      !if(dx**2+dz**2.lt.dr**2) then
      !  a(i)=a0+0.005; vc(i)=vc0*1d6
      !end if
      !end do
      bg(i)=b0
      dcg(i)=dc0
      vcg(i)=vc0
      ! do k=1,kmax
      !   dy=ycoll(i)-yc(k);dz=zcoll(i)-zc(k)
      !   if(dy**2+dz**2.lt.dr**2) vc(i)=vc0*1d4
      ! end do
      !taudot(i)=vpl*rigid/(1.5d0*pi)/xcoll(i)
    end do
  end subroutine
end module
