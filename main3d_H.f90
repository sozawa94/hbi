program main
  !This program solves the EQ nucleation process on a 3D nonplanar fault.
  !The fault is discretized into many triangular elements to calculate the elastic
  !interaction using the triangular dislocation theory (TDstressFS.f90)
  !Rate and State Friction (aging law, slip law, Nagata law)
  !length is normalized by Dc/theta*=10^(-5)m
  !time is normalized by theta*=1(sec)
  !stress is normalized by sigma0=100MPa

  !developed by SO OZAWA (Master Student, Earthquake Research Institute, UTokyo)

  !$ use omp_lib
  use dtriangular
  implicit none
   include 'mpif.h'
  integer::NCELL, nstep1, lp, i,j,k,m,counts,interval,number,rc,array(2)
  integer::clock,cr,counts2,NVER
  integer::r1,r2,r3,s1,s2,s3
  real(8)::p1(3),p2(3),p3(3)
  character*128::filename,dum,law
  real(8)::a0,b0,sr,omega,theta,dtau,tiny,x,time1,time2,moment,aslip,avv
  real(8)::factor,rigid
  real(8)::eps,vpl,ST(6),r
  real(8)::dtime,dtnxt,dttry,dtdid,alpha,ds,amp,mui,strinit,velinit,velmax
  real(8),parameter::pi=4.d0*atan(1.d0),pois=0.25d0
  real(8),parameter::dc0=1.d0,sigma0=1.0d0,mu0=0.60d0,tref=1.d9
  real(8),parameter::dtinit=10d0,vref = 1.d0,minsig=0.05
  integer,parameter::int1=1
  real(8),allocatable::a(:),b(:),dc(:),sigma(:),vel(:),tau(:),disp(:),gshear(:,:)
  real(8),allocatable::xcol(:),ycol(:),xe(:),ye(:),taudot(:),sigdot(:),ang(:),yr(:)
  real(8),allocatable::gnorm(:,:),y(:),yscal(:),s(:),xtime1(:), dtime1(:),dydx(:)
  real(8),allocatable::atau(:),mu(:),xr(:),angs(:),ze(:),zcol(:)
  real(8),allocatable::ai(:,:)
  integer,allocatable::bi(:,:)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )

  !input parameters
  open(33,file='input.dat')
  read(33,*) dum,NCELL !number of cells
  read(33,*) dum,NVER !number of vertex
  read(33,*) dum,nstep1 !maxmimum time step
  read(33,*) dum,number !sample number
  read(33,*) dum,interval !minimum roughness scale
  read(33,*) dum,amp !roughness amplitude beta
  read(33,*) dum,strinit !initial applied stress(normalized)
  read(33,*) dum,velinit !initial slip velocity
  read(33,*) dum,omega !initial omega=V*theta
  read(33,*) dum,ds !mesh interval(normalized by Dc)
  read(33,*) dum,a0 !a in RSF
  read(33,*) dum,b0 !b in RSF
  read(33,*) dum,alpha !alpha in RSF
  read(33,*) dum,sr !loading rate
  read(33,*) dum,velmax !stop when Vmax exceeds this value
  read(33,*) dum,rigid !rigidity
  read(33,*) dum,eps !error allowance
  read(33,*) dum,law ! evolution law
  read(*,*) number
  !read(*,*) omega
  close(33)
  !write(*,*) strinit,velinit,sr,law

  allocate(a(NCELL),b(NCELL),dc(NCELL),sigma(NCELL),vel(NCELL),tau(NCELL))
  allocate(disp(NCELL),xcol(NCELL),ycol(NCELL),zcol(NCELL),xe(Nver),ye(Nver),ze(Nver))
  allocate(gshear(NCELL,NCELL),taudot(NCELL),sigdot(NCELL),ang(NCELL),angs(NCELl))
  allocate(gnorm(NCELL,NCELL),y(3*NCELL),yscal(3*NCELL),s(NCELL))
  allocate(xtime1(NSTEP1), dtime1(NSTEP1),dydx(3 * NCELL),yr(0:NCELL),xr(0:NCELL))
  allocate(atau(NCELL),mu(NCELL))
  allocate(ai(NCELl,3),bi(NCELl,3))

  !fault's geometry
  open(15,file='mesh.dat')
  open(16,file='tri.dat')

  do i=1,NVER
    read(15,*) ye(i),ze(i)
    call random_number(r)
    ye(i)=ye(i)*ds
    ze(i)=ze(i)*ds
    xe(i)=r*1d-15
    !xe(i)=0.d0
  end do
  do k=1,NCELL
    read(16,*) ai(k,1:3)
    bi(k,1:3)=int(ai(k,1:3))
    !write(*,*) a(k,1:3),b(k,1:3)
  end do
  close(15)
  close(16)

  !The calculation point is the center of each triangle
  do i=1,NCELL
    r1=bi(i,1)
    r2=bi(i,2)
    r3=bi(i,3)
    xcol(i)=(xe(r1)+xe(r2)+xe(r3))/3.d0
    ycol(i)=(ye(r1)+ye(r2)+ye(r3))/3.d0
    zcol(i)=(ze(r1)+ze(r2)+ze(r3))/3.d0
  end do
  write(*,*) 'mesh generated'

  !frictional parameters
  a=a0
  b=b0
  dc=dc0

  !calculate green's function
  do i=1,NCELL
    do j=1,NCELL
      s1=bi(j,1)
      s2=bi(j,2)
      s3=bi(j,3)
      p1=(/xe(s1),ye(s1),ze(s1)/)
      p2=(/xe(s2),ye(s2),ze(s2)/)
      p3=(/xe(s3),ye(s3),ze(s3)/)
      call  TDstressFS(xcol(i),ycol(i),zcol(i),P1,P2,P3,1.d0,0.d0,0.d0,rigid,rigid, &
      & St(1),St(2),St(3),St(4),St(5),St(6))
      gshear(i,j)=St(4)
      gnorm(i,j)=St(1)
    end do
  end do
  write(*,*) maxval(gshear)
  write(*,*) 'kernel calculated'

  !stressing rate
  do i=1,NCELL
    taudot(i)=sr
    !taudot(i)=sr*cos(2*(ang(i)-0.5d0*atan(0.63d0)))
    !taudot(i)=sr*cos(2*ang(i))
    !sigdot(i)=sr*sin(2*(ang(i)-0.5d0*atan(0.63d0)))
    !sigdot(i)=sr*sin(2*ang(i))
    sigdot(i)=0
  end do

  !output settings
  !a:spatial information
  !b:temporal information
  !call system_clock(clock)
  write(filename, '("output/",i0,"a.dat")') number
  open(11,file=filename)
  write(filename, '("output/",i0,"b.dat")') number
  open(12,file=filename)

  open(17,file='initial.dat')

  !initial condition
  theta=omega/velinit
  mui=mu0+a0*dlog(velinit/vref)+b0*dlog(theta)

  do i=1,NCELL
    sigma(i)=sigma0
    !sigma(i)=sigma0+mui*sigma0*sin(2*ang(i))
    tau(i)=mui*sigma(i)
    s(i)=theta
    vel(i)=exp((tau(i)/sigma(i)-mu0-b(i)*log(s(i)))/a(i))
    disp(i)=0.d0
    write(17,'(i8,6e15.6)') i,xcol(i),ycol(i),zcol(i),tau(i),sigma(i),vel(i)
    !& taudot(i),sigdot(i),vel(i),s(i),vel(i)*s(i)
  end do
  close(17)

  write(*,*) 'initial velocity',velinit
  write(*,*) 'initial theta',theta
  write(*,*) 'initial omega',velinit*theta
  !stop

  write(*,*) 'start time integration!'
  !call cpu_time(time1)
  !$ time1=omp_get_wtime()
  !write(*,*) 'step,time,velocity'

  x=0.d0
  dtnxt = dtinit

  !start time integration
  do k=1,NSTEP1 !k : time step
    do i=1,NCELL
      y(3*i-2) = dlog(vel(i))
      y(3*i-1) = tau(i)
      y(3*i) = sigma(i)
    end do

    dttry = dtnxt
    dtime = 0.0

    call derivs( NCELL, x, y, dydx,a, b, dc, sigma, gshear, gnorm,sr, vref, mu0,law,rigid)
    !write(*,*) 'call derivs'
    do i = 1, 3 * NCELL
      yscal(i) = abs( y(i) ) + abs( dttry * dydx(i) ) + tiny
    end do

    call rkqs( NCELL, y, dydx, 3 * NCELL, x, dttry, eps,yscal, dtdid, &
    & dtnxt, derivs, a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid)

    dtime = dtime + dtdid

    do i = 1, NCELL
      vel(i) = exp(y(3*i-2))
      tau(i) = y(3*i-1)
      sigma(i) = y(3*i)
      disp(i) = disp(i) + exp( y(3 * i - 2) ) * dtdid
      s(i)=exp((tau(i)/sigma(i)-mu0-a(i)*dlog(vel(i)/vref))/b(i))
      mu(i)=tau(i)/sigma(i)
    end do

    ! If there is negative normal stress, stop
    if(minval(sigma).lt.minsig) stop

    dttry = dtnxt

    !array=maxloc(vel)
    m=NCELL/2
    !nucleation lemgth
    counts=0
    do i=1,NCELL
      if(vel(i).gt.1d0) counts=counts+1
    end do

    !iter,time,maxvel,maxst,maxloc
    !adisp=sum(disp)/NCELL
    !avel=rigid*ds*sum(vel)
    write(12,'(i5,f16.4,e12.4)')k,x,maxval(vel)
    if(mod(k,30).eq.0) then
      do i=1,NCELL
        write(11,'(8e15.6,i7,2f15.4)') ycol(i),zcol(i),log10(vel(i)),tau(i),sigma(i),disp(i),&
        & s(i),tau(i)/sigma(i),k,x,s(i)*vel(i)
        !write(13,'(i5,6e15.6,i7,2f15.4)') i,log10(vel(i)),tau(i),sigma(i),disp(i),&
        !& s(i),tau(i)/sigma(i),k,x,s(i)*vel(i)
      end do
      write(11,*)
      write(11,*)
      !write(13,*)


      write(*,'(i10,f17.5,e12.5)') k,x,maxval(vel)
    end if
    !write(*,*) 'nucleation size',counts
    ! m=maxloc(vel)
    if(maxval(vel).gt.velmax) then
      cr=cr+1
    else
      cr=0
    end if

    if(cr.gt.10) exit
    !write(13,*)
    dtime1(k) = dtime
    xtime1(k) = x
    !     write(12,'(2e12.5)') xtime1(k)*tref, dtime1(k)*tref
  end do                    ! k = 1, NSTEP1 : outside time loop end

  !-----final results---------------
  !open(30,file='result.txt',position='append')

  !nucleation lemgth
  !counts=0
  !do i=1,NCELL
  !  if(vel(i).gt.1d0) counts=counts+1
  !end do

  !counts2=0
  !do i=1,NCELL
  !  if(vel(i).gt.1d-4) counts2=counts2+1
  !end do

  !released moment
  !rc=0
  !aslip=0
  !do i=1,10
  !  if(disp(250*i).lt.5.0) then
  !    rc=rc+1
  !    aslip=aslip+disp(250*i)
  !  end if
  !end do
  !aslip=aslip/rc
  !moment = (sum(disp)/NCELL-aslip)*NCELL

  !write(30,'(f10.4,3i6)') amp,number,counts,counts2
  !write(*,*) 'nucleation size',counts
  !write(*,*) 'nucleation length',counts
  !call cpu_time(time2)
  !$ time2=omp_get_wtime()
  !write(*,*) 'cpu time', time2-time1
  !go to 100

  stop
contains
  !---------------------------------------------------------------------
  subroutine derivs( NCELL, x, y, dydx,a, b, dc, sigma, gshear, gnorm,sr, vref, mu0,law,rigid)
    !---------------------------------------------------------------------
    implicit none
    !$use omp_lib
    integer :: NCELL
    real(8) :: x, y(3 * NCELL), dydx(3 * NCELL)
    real(8) :: a(NCELL), b(NCELL), dc(NCELL), sigma(NCELL)
    real(8) :: gshear(NCELL,NCELL), gnorm(NCELL,NCELL)
    real(8) :: sr, vref, mu0, s(NCELL)
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),vs,vc,v2
    real(8) :: dlnvdt(NCELL),dtaudt(NCELL),dsigdt(NCELL)
    real(8) :: sum_gs(NCELL), sum_gn(NCELL), tau0, c, c1, c2, c3, arg, rigid, alpha
    character(128) :: law
    integer :: i, j, nc

    !     -- Runge-Kutta copled ODE dy(i)/dx definition --
    do i = 1, NCELL
      veltmp(i) = exp(y(3*i-2))
      tautmp(i) = y(3*i-1)
      sigmatmp(i) = y(3*i)
    enddo
    !write(*,*) 'maxval(veltmp)',maxval(veltmp)
    vs=3.d8
    sum_gs = 0.d0
    sum_gn = 0.d0
    !write(*,*) 'call omp'
    !$omp parallel do private(j)
    do i = 1, NCELL
      do j = 1, NCELL
        sum_gn(i) = sum_gn(i) + gnorm(i,j) * veltmp(j)
        sum_gs(i) = sum_gs(i) + gshear(i,j) * veltmp(j)
      end do
    end do
    !$omp end parallel do
    !write(*,*) maxval(sum_gs)
    !$omp parallel do
    do i=1,NCELL
      dsigdt(i)=sum_gn(i)+sigdot(i)
      !dsigdt(i)=sum_gn(i)
      sum_gs(i)=sum_gs(i)+taudot(i)
      !    -- equation dtau/dt --
      !write(*,*) i,veltmp(i)
      !dsigdt(i) = sigdot(i) + sum(gnorm(i,:) * veltmp(:))
      !sum_gn=sigdot(i) + sum(gnorm(i,:) * veltmp(:))
      !dtaudt(i) = taudot(i) + sum(gshear(i,:) * veltmp(:))
      !sum_gs=taudot(i) + sum(gshear(i,:) * veltmp(:))
      !write(*,*) dsigdt(i),dtaudt(i)
      !radiation damping
      !rigid=1.d4
      !     -- equation dtau/dt --
      select case(law)
        ! aging law (Linker & Dieterich, 1992)
      case('a')
        c1 = b(i) / dc(i) * veltmp(i)
        arg = tautmp(i)/sigmatmp(i) - mu0 - a(i) * dlog( veltmp(i) / vref)
        c2 = b(i) / dc(i) * vref * exp( - arg / b(i) ) !b/theta
        c3 = alpha/sigmatmp(i)*dsigdt(i)

        dtaudt(i)=sum_gs(i) + 0.5d0*rigid/vs*veltmp(i)/a(i)*((c2-c1-c3)+ &
        & tautmp(i)/sigmatmp(i)**2*dsigdt(i))
        dtaudt( i ) = dtaudt(i)/(1.d0+rigid*veltmp(i)*0.5d0/vs/a(i)/sigmatmp(i))
        !dtaudt( i ) = (sr+sum_gs)*a(i)/veltmp(i)+rigid*0.5d0/vs*(c2-c1)
        !dlnvdt(i)=2.d0*vs/rigid/veltmp(i)*(dtaudt(i)-sum_gs)
        !dlnvdt(i)=(dtaudt(i)/sigma(i)+c1-c2+c3)/a(i)
        dlnvdt( i ) = ( dtaudt(i) /sigmatmp(i) - dsigdt(i)*tautmp(i)/ &
        & sigmatmp(i)**2+ c1 - c2 +c3) / a(i)

        ! slip law
      case('s')
        !if(vel(i).gt.1d4) dc(i)=
        arg = tautmp(i)/sigmatmp(i) - mu0 - (a(i)-b(i)) * dlog( veltmp(i) / vref)
        c1 = veltmp(i) / dc(i) * arg
        c3 = alpha/sigmatmp(i)*dsigdt(i)
        dtaudt(i)=sum_gs(i) + 0.5d0*rigid/vs*veltmp(i)/a(i)*((-c1-c3)+ &
        & tautmp(i)/sigmatmp(i)**2*dsigdt(i))
        dtaudt( i ) = dtaudt(i)/(1.d0+rigid*veltmp(i)*0.5d0/vs/a(i)/sigmatmp(i))
        !dtaudt( i ) = (sr+sum_gs)*a(i)/veltmp(i)-rigid*0.5d0/vs*c1
        !dtaudt( i ) = dtaudt(i)/(rigid*0.5d0/vs+a(i)/veltmp(i))
        dlnvdt( i ) = ( dtaudt(i) /sigmatmp(i) - dsigdt(i)*tautmp(i)/ &
        & sigmatmp(i)**2+ c1 + c3) / a(i)
        !dlnvdt( i ) = ( dtaudt(i) / sigma(i)- + c1 )/  a(i)

        !composite law
      case('c')
        dtaudt(i)=sr+sum_gs(i)
        vc=1d-2
        arg = tautmp(i) - tau0 - (a(i)-b(i)) * dlog( veltmp(i) / vref)
        s(i)=exp((tautmp(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))
        c1 = veltmp(i) / dc(i) * arg-b(i)/s(i)*exp(-veltmp(i)/vc)
        dlnvdt( i ) = ( dtaudt(i)  + c1 ) / a(i)

        !nagata law
      case('n')
        c=2.0
        c1 = b(i) / dc(i) * veltmp(i)
        arg = tautmp(i)/sigmatmp(i) - mu0 - a(i) * dlog( veltmp(i) / vref)
        c2 = b(i) / dc(i) * vref * exp( - arg / b(i) ) !b/theta
        c3 = alpha/sigmatmp(i)*dsigdt(i)
        dtaudt(i)=sum_gs(i)
        dlnvdt(i)=(dtaudt(i)*sigmatmp(i)-tautmp(i)*dsigdt(i))/sigmatmp(i)**2- &
        & c2+c1+c3+c*dtaudt(i)
        dlnvdt(i)=dlnvdt(i)/a(i)
      end select
      !write(*,*) dlnvdt(1000)
    enddo
    !write(*,*) maxval(dlnvdt),maxval(dtaudt),maxval(dsigdt)

    do i = 1, NCELL
      dydx(3*i-2) = dlnvdt( i )
      dydx(3*i-1) = dtaudt( i )
      dydx(3*i) = dsigdt( i )
    enddo

    return
  end subroutine
  !---------------------------------------------------------------------
  subroutine rkqs( NCELL, y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs,&
    &   a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid)
    !---------------------------------------------------------------------
    implicit none
    integer :: NCELL
    real(8) :: a(NCELL), b(NCELL), dc(NCELL), sigma(NCELL)
    real(8) :: gshear(NCELL,NCELL),gnorm(NCELL,NCELL), sr, vref, mu0,rigid
    character*128::law

    INTEGER :: n,nmax
    REAL(8) :: eps,hdid,hnext,htry,x,htemp
    real(8) :: dydx(3 * NCELL), y(3 * NCELL), yscal(3 * NCELL)
    EXTERNAL derivs
    PARAMETER  (NMAX = 30000)
    INTEGER :: i
    REAL(8) :: errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON,hmax
    PARAMETER (SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89e-4)
    parameter(hmax=1d5)
    h=htry
    !write(*,*) 'call rkqs'
    do while(.true.)
      !write(*,*) 'call rkck'
      call rkck( NCELL,y, dydx, n, x, h, ytemp, yerr, derivs,a, b, dc, sigma,&
      & gshear, sr, vref, mu0,law,rigid)
      errmax=maxval(abs(yerr(:)/yscal(:)))/eps
      !write(*,*) errmax
      if(errmax.lt.1.d0) exit
      htemp=SAFETY*h*(errmax**PSHRNK)
      h=sign(max(abs(htemp),0.1*abs(h)),h)
      xnew=x+h
      !write(*,*) h,errmax
      if(xnew-x<0.000000001) stop
      !    write(*,*) x,k
      !    stop
      !  end if
    end do

    !if(errmax.gt.ERRCON) then
    hnext=min(1.5*h,SAFETY*h*(errmax**PGROW),5e6)

    hdid=h
    x=x+h
    y(:)=ytemp(:)
    return
  end subroutine

  !---------------------------------------------------------------------
  subroutine rkck( NCELL,y, dydx, n, x, h, yout, yerr, derivs,&
    &     a, b, dc, sigma, gshear, sr, vref, mu0,law,rigid)
    !---------------------------------------------------------------------
    implicit none
    integer :: NCELL
    real(8) :: a(NCELL), b(NCELL), dc(NCELL), sigma(NCELL)
    real(8) :: gshear(NCELL,NCELL), sr, vref, mu0,rigid
    character*128::law

    INTEGER :: n, NMAX
    REAL(8) :: h,x
    real(8) :: dydx(3 * NCELL), y(3 * NCELL)
    real(8) :: yerr(n), yout(n)
    EXTERNAL :: derivs
    PARAMETER (NMAX = 30000)
    INTEGER :: i
    REAL(8) :: ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX)
    real(8) :: ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51
    real(8) :: B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.)
    parameter (B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5)
    parameter (B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.)
    parameter (B63=575./13824.,B64=44275./110592.,B65=253./4096.)
    parameter (C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.)
    parameter (DC1=C1-2825./27648.,DC3=C3-18575./48384.)
    parameter (DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)

    !     -- 1st step --
    !$omp parallel do
    do i=1,n
      ytemp(i)=y(i)+B21*h*dydx(i)
      !write(*,*) i,ytemp(i)
    end do
    !stop
    !write(*,*) maxval(ytemp)

    !    -- 2nd step --

    call derivs(NCELL,x+A2*h,ytemp,ak2, a, b, dc, sigma, gshear, gnorm, sr, vref, mu0,law ,rigid)
    !$omp parallel do
    do i=1,n
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      !write(*,*) i,ytemp(i)
    end do

    !write(*,*) maxval(ytemp),maxval(ak2)

    !     -- 3rd step --
    call derivs(NCELL,x+A3*h,ytemp,ak3, a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid )
    !write(*,*) '3rd step',maxval(ak3),maxval(ak2),maxval(dydx),maxval(y)
    !$omp parallel do
    do i=1,n
      ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      !write(*,*) i,y(i),ytemp(i),dydx(i),ak2(i),ak3(i)
    end do
    !write(*,*) ytemp
    !stop
    !     -- 4th step --
    call derivs(NCELL,x+A4*h,ytemp,ak4,a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid )
    !write(*,*) '4th step',maxval(y),maxval(ak4)
    !$omp parallel do
    do i=1,n
      ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
      !write(*,*) i,ytemp(i)
    end do
    !write(*,*) maxval(ytemp)

    !     -- 5th step --
    call derivs(NCELL,x+A5*h,ytemp,ak5,a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law ,rigid)
    !$omp parallel do
    do i=1,n
      ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do
    !     -- 6th step --
    call derivs(NCELL,x+A6*h,ytemp,ak6, a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid )
    !$omp parallel do
    do i=1,n
      yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
    end do

    !$omp parallel do
    do i=1,n
      yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
    end do
    !write(*,*) yerr

    return
  end subroutine
end program
