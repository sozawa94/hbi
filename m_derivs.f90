module mod_derivs
  !module for time integration
  implicit none
  type::t_deriv
    real(8),pointer::a(:),b(:),dc(:),vc(:),sigma(:),taudot(:),sigmdot(:)
    real(8)::rigid,vs,pois,mu0,vref,vw,fw
    integer::NCELL,NCELLg,load
    integer,pointer::rcounts(:),displs(:)
    character(128)::law
  end type
contains
  subroutine derivs(x, y, dydx,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    type(t_deriv),intent(in) :: st_deriv
    !integer,intent(in) :: NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in) :: x
    real(8),intent(in) ::y(:)
    real(8),intent(out) :: dydx(:)
    real(8) :: veltmp(st_deriv%NCELL),tautmp(st_deriv%NCELL),sigmatmp(st_deriv%NCELL)
    real(8) :: dlnvdt(st_deriv%NCELL),dtaudt(st_deriv%NCELL),dsigdt(st_deriv%NCELL)
    real(8) :: sum_gs(st_deriv%NCELL),sum_gn(st_deriv%NCELL)
    real(8)::veltmpg(st_deriv%NCELLg),sum_gsg(st_deriv%NCELLg),sum_gng(st_deriv%NCELLg)
    real(8):: c1, c2, c3, arg,c,g,tauss,dphidt
    integer :: i, j, nc,ierr,lrtrn

    !if(my_rank.eq.0) then
    select case(problem)
    case('2dp','3dp')
      do i = 1, st_deriv%NCELL
        veltmp(i) = exp (y(2*i-1))
        tautmp(i) = y(2*i)
        sigmatmp(i)=sigma(i)
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(veltmp,st_deriv%NCELL,MPI_REAL8,veltmpG,st_deriv%rcounts,st_deriv%displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      select case(st_deriv%load)
      case(0)
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
      case(1)
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG-st_deriv%vpl)
      end select
      !call MPI_SCATTERv(sum_gsG,NCELL,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_gsg,st_deriv%rcounts,st_deriv%displs,MPI_REAL8,sum_gs,st_deriv%NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      do i=1,st_deriv%NCELL
        sum_gn(i)=0.d0
        sum_gs(i)=sum_gs(i)+st_deriv%taudot(i)
      end do

    case('2dn')
      do i = 1, st_deriv%NCELL
        veltmp(i) = exp (y(3*i-2))
        tautmp(i) = y(3*i-1)
        sigmatmp(i) = y(3*i)
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(Veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,                &
      &     MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      st_bemv%v='s'
      !veltmp=1.d-6
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
      !if(my_rank.eq.0) write(*,*) sum_gs
      !stop
      st_bemv%v='n'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxpn,st_bemv,st_ctl,sum_gnG,veltmpG)
      !if(my_rank.eq.10) write(*,*) sum_gn(1)

      call MPI_SCATTERv(sum_gsG,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_gnG,rcounts,displs,MPI_REAL8,sum_gn,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      do i=1,st_deriv%NCELL
        sum_gs(i)=sum_gs(i)+st_deriv%taudot(i)
        sum_gn(i)=sum_gn(i)+st_deriv%sigdot(i)
      end do
    end select
    !end if



    select case(st_deriv%law)
      ! aging law (Linker & Dieterich, 1992)
    case('a')
      call deriv_a(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,st_deriv,dlnvdt,dtaudt,dsigdt)

      ! slip law
    case('s')
      call deriv_s(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,st_deriv,dlnvdt,dtaudt,dsigdt)

      ! RFL in FDMAP (Dunham+ 2011)
    case('d')
      call deriv_d(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,st_deriv,dlnvdt,dtaudt,dsigdt)

    end select

    select case(problem)
    case('2dp','3dp')
      do i = 1, st_deriv%NCELL
        dydx(2*i-1) = dlnvdt( i )
        dydx(2*i) = dtaudt( i )
      enddo
    case('2dn')
      do i = 1, st_deriv%NCELL
        dydx(3*i-2) = dlnvdt( i )
        dydx(3*i-1) = dtaudt( i )
        dydx(3*i) = dsigdt( i )
      enddo
    end select

    return
  end subroutine

  subroutine deriv_a(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,st_deriv,dlnvdt,dtaudt,dsigdt)
    implicit none
    integer::i
    real(8)::arg
    type(t_deriv),intent(in) :: st_deriv
    real(8),intent(in)::sum_gs(:),sum_gn(:),veltmp(:),tautmp(:),sigmatmp(:)
    !real(8),intent(in)::a(:),b(:),dc(:),vc(:)
    !real(8),intent(in)::mu0,vref,vs,rigid,alpha
    real(8),intent(out)::dlnvdt(:),dtaudt(:),dsigdt(:)
    do i=1,size(sum_gs)
      dsigdt(i)=sum_gn(i)
      arg=dc(i)/vref*(exp((tautmp(i)/sigmatmp(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
      dlnvdt(i)=sum_gs(i)-b(i)*sigmatmp(i)/(arg+dc(i)/vc(i))*(1.d0-veltmp(i)*arg/dc(i))+(tautmp(i)/sigmatmp(i)-alpha)*dsigdt(i)
      dlnvdt(i)=dlnvdt(i)/(a(i)*sigmatmp(i)+0.5d0*rigid*veltmp(i)/vs)
      dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
    end do
  end subroutine

  subroutine deriv_s(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,st_deriv,dlnvdt,dtaudt,dsigdt)
    implicit none
    integer::i
    real(8)::arg
    type(t_deriv),intent(in) :: st_deriv
    real(8),intent(in)::sum_gs(:),sum_gn(:),veltmp(:),tautmp(:),sigmatmp(:)
    !real(8),intent(in)::a(:),b(:),dc(:),vc(:)
    !real(8),intent(in)::mu0,vref,vs,rigid,alpha
    real(8),intent(out)::dlnvdt(:),dtaudt(:),dsigdt(:)
    do i=1,size(sum_gs)
      dsigdt(i)=sum_gn(i)
      arg=dc(i)/vref*(exp((tautmp(i)/sigmatmp(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
      dlnvdt(i)=sum_gs(i)+b(i)*sigmatmp(i)*veltmp(i)*dlog(veltmp(i)*arg/dc(i))+(tautmp(i)/sigmatmp(i)-alpha)*dsigdt(i)
      dlnvdt(i)=dlnvdt(i)/(a(i)*sigmatmp(i)+0.5d0*rigid*veltmp(i)/vs)
      dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
    end do
  end subroutine

  subroutine deriv_d(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,st_deriv,dlnvdt,dtaudt,dsigdt)
    implicit none
    integer::i
    real(8)::theta,g,fss,dphidt,psi
    type(t_deriv),intent(in) :: st_deriv
    real(8),intent(in)::sum_gs(:),sum_gn(:),veltmp(:),tautmp(:),sigmatmp(:)
    !real(8),intent(in)::a(:),b(:),dc(:)
    !real(8),intent(in)::mu0,vref,vs,fw,vw,rigid
    real(8),intent(out)::dlnvdt(:),dtaudt(:),dsigdt(:)
    do i=1,size(sum_gs)
      dsigdt(i)=sum_gn(i)
      !theta=tautmp(i)/sigmatmp(i)-a(i)*dlog(veltmp(i)/vref)
      psi=a(i)*dlog(2.d0*vref/veltmp(i)*dsinh(tautmp(i)/sigmatmp(i)/a(i)))
      !if(i.eq.500) write(*,*) 'vel,tau,sigma,psi',veltmp(500),tautmp(500),sigmatmp(500),psi
      g=0.5d0*veltmp(i)/vref*exp(psi/a(i))
      fss=mu0+(a(i)-b(i))*dlog(veltmp(i)/vref)
      fss=fw+(fss-fw)/(1.d0+(veltmp(i)/vw)**8)**0.125d0 !flash heating
      dphidt=-veltmp(i)/dc(i)*(tautmp(i)/sigmatmp(i)-fss)
      dlnvdt(i)=sum_gs(i)/sigmatmp(i)-g/sqrt(1.d0+g**2)*dphidt-tautmp(i)/sigmatmp(i)**2*dsigdt(i)
      dlnvdt(i)=dlnvdt(i)/(a(i)*g/sqrt(1.d0+g**2)+0.5d0*rigid*veltmp(i)/vs/sigmatmp(i))
      dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
    end do

  end subroutine

  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    !integer::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::yscal(:),htry,eps
    real(8),intent(inout)::y(:),x,dydx(:)
    real(8),intent(out)::hdid,hnext
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    type(t_deriv),intent(in) :: st_deriv
    integer :: i,ierr
    real(8) :: errmax,h,xnew,htemp,errmax_gb
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4,hmax=1d5

    h=htry
    do while(.true.)
      call rkck(y,dydx,x,h,ytemp,yerr,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)!,derivs)
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
  subroutine rkck(y,dydx,x,h,yout,yerr,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    !integer,intent(in)::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::y(:),dydx(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
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
    call derivs(x+a2*h, ytemp, ak2,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
    end do

    !     -- 3rd step --
    call derivs(x+a3*h, ytemp, ak3,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
    end do

    !     -- 4th step --
    call derivs(x+a4*h, ytemp, ak4,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
    end do

    !     -- 5th step --
    call derivs(x+a5*h, ytemp, ak5,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do

    !     -- 6th step --
    call derivs(x+a6*h, ytemp, ak6,st_deriv)!,st_deriv,st_leafmtxp,st_bemv,st_ctl)
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

end module
