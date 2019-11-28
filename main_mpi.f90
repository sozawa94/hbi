  program main
    !This program solves the EQ nucleation process on a 2D nonplanar fault with
    !numerous microfaults.
    !MPI calculation.
    !Rate and State Friction (aging law, slip law, Nagata law)
    !length is normalized by Dc/theta*=10^(-5)m
    !time is normalized by theta*=1(sec)
    !stress is normalized by sigma0=100MPa

    !developed by SO OZAWA (Master Student, Earthquake Research Institute, UTokyo)
  !$ use omp_lib
   implicit none
   include 'mpif.h'

   integer::NCELL, nstep1, lp, i,j,k,m,counts,interval,number,rc,NCELLg
   integer::clock,cr,counts2,imax,jmax,NCELLm,ierr,my_rank,np,nfile,nfmax
   integer::istart,iend
   character*128::fname,dum,law
   real(8)::a0,b0,sr,omega,theta,dtau,tiny,x,time1,time2,moment,aslip,avv
   real(8)::factor,dx,dy,angle,sxy,sxx,syy,gtau(2),gsig(2),rigid
   real(8)::c,d,gx,gy,gt,gn,r(2),avel,eps,adisp,vpl,w
   real(8)::dtime,dtnxt,dttry,dtdid,alpha,ds,amp,mui,strinit,velinit,velmax
   real(8),parameter::pi=4.d0*atan(1.d0),pois=0.25d0
   real(8),parameter::dc0=1.d0,sigma0=1.0d0,mu0=0.60d0,tref=1.d9
   real(8),parameter::dtinit=10d0,vref = 1.d0,minsig=0.05
   integer,parameter::int1=1
   real(8),allocatable::a(:),b(:),dc(:),sigma(:),vel(:),tau(:),disp(:),gshear(:,:)
   real(8),allocatable::xcol(:),ycol(:),taudot(:),sigdot(:),ang(:)
   real(8),allocatable::xel(:),xer(:),yel(:),yer(:),angs(:),xcolg(:),ycolg(:)
   real(8),allocatable::gnorm(:,:),y(:),yscal(:),s(:),xtime1(:), dtime1(:),dydx(:)
   real(8),allocatable::atau(:),mu(:)
   real(8),allocatable::sigmag(:),velg(:)

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
   call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )

  !if(my_rank.eq.0) then
    open(33,file='input.dat')
    read(33,*) dum,NCELLm !number of cells (main fault)
    read(33,*) dum,imax
    read(33,*) dum,jmax
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
    read(33,*) dum,rigid !rigidity normalized by 100MPa
    read(33,*) dum,eps !error allowance
    read(33,*) dum,law ! evolution law
    !read(*,*) number
    !read(*,*) omega
    close(33)
    write(*,*) 'read'
    NCELLg=NCELLm+imax*jmax
   !end if

   NCELL=NCELLg/np
   istart=my_rank*NCELL+1
   iend=(my_rank+1)*NCELL
   !scatter NCELL,NCELLg,a,b,alpha,dc,sr,rigid,omega,velinit,eps,law
   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call MPI_SCATTER(a0g,1,MPI_REAL8,a0,1,MPI_REAL8,iroot,MPI_COMM_WORLD,ierr)

   allocate(a(NCELL),b(NCELL),dc(NCELL),sigma(NCELL),vel(NCELL),tau(NCELL))
   allocate(disp(NCELL),xcol(NCELL),ycol(NCELL),ang(NCELL),angs(NCELLg))
   allocate(xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg),xcolg(NCELLg),ycolg(NCELLg))
   allocate(gshear(NCELL,NCELLg),gnorm(NCELL,NCELLg),taudot(NCELL),sigdot(NCELL))
   allocate(y(3*NCELL),yscal(3*NCELL),s(NCELL),mu(NCELL))
   allocate(xtime1(NSTEP1), dtime1(NSTEP1),dydx(3*NCELL))
   allocate(sigmaG(NCELLg),velG(NCELLg))

   !geometry
   if(my_rank.eq.0) then
   Open(11,file='top.dat')
   do i=1,NCELLg
     read(11,*) xcolg(i),ycolg(i),angs(i),xel(i),xer(i),yel(i),yer(i)
   end do
   close(11)
   end if

   call MPI_SCATTER(xcolg,NCELL,MPI_REAL8,xcol,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
   call MPI_SCATTER(ycolg,NCELL,MPI_REAL8,ycol,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
   call MPI_SCATTER(angs,NCELL,MPI_REAL8,ang,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(xel, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(xer, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(yel, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(yer, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

   !green function
   factor = rigid / ( 2.d0 * pi * ( 1.d0 - pois ) )
   do i=1, NCELL
     do j=1, NCELLg
       angle=ang(i)-angs(j) ! check sign
       !write(*,*) angle
        dx=xcol(i)-xel(j)
        dy=ycol(i)-yel(j)
        !write(*,*) dx,dy
        sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
        sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
        syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
        !write(*,*) sxx,sxy,syy
        gtau(1)=-0.5d0*(sxx-syy)*sin(2.d0*angle)+sxy*cos(2.d0*angle)
        gsig(1)=0.5d0*(sxx+syy)-0.5d0*(sxx-syy)*cos(2.d0*angle)-sxy*sin(2.d0*angle)

        dx=xcol(i)-xer(j)
        dy=ycol(i)-yer(j)
        !write(*,*) dx,dy
        sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
        sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
        syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
        !write(*,*) sxx,sxy,syy
        gtau(2)=-0.5d0*(sxx-syy)*sin(2.d0*angle)+sxy*cos(2.d0*angle)
        gsig(2)=0.5d0*(sxx+syy)-0.5d0*(sxx-syy)*cos(2.d0*angle)-sxy*sin(2.d0*angle)

      gshear(i,j)=gtau(2)-gtau(1)
      gnorm(i,j)=gsig(2)-gsig(1)
     end do
   end do
   !call green(gshear,gnorm,NCELLg,NCELL,xcol,ycol,xel,xer,yel,yer,ang,angs,factor)
   write(*,*) 'kern calculated'
   !if(my_rank.eq.0) write(*,*) gshear(1,1),gnorm(1,1)
   !frictional parameters
    a=a0
    b=b0
    dc=dc0

   !stressing rate
   do i=1,NCELL
     !taudot(i)=sr*cos(2*(ang(i)-0.5d0*atan(0.63d0)))
     taudot(i)=sr*cos(2*ang(i))
     !sigdot(i)=sr*sin(2*(ang(i)-0.5d0*atan(0.63d0)))
     sigdot(i)=sr*sin(2*ang(i))
   end do

   !Do nfile=11,nfmax
     nfile=my_rank+11
     Write(fname,'(I2.2,A)') my_rank,'b.dat'
     Open(nfile,file="output/"//fname)
   !EndDo
   !initial condition
   !mu=mu0+strinit
   !theta=exp((strinit-a0*dlog(velinit/vref))/b0)
   theta=omega/velinit
   mui=mu0+a0*dlog(velinit/vref)+b0*dlog(theta)

   do i=1,NCELL
    sigma(i)=sigma0+mui*sigma0*sin(2*ang(i))
    tau(i)=mui*sigma(i)
    s(i)=theta
    !s(i)=exp((tau(i)/sigma(i)-mu0-a(i)*dlog(vel(i)/vref))/b(i))
    vel(i)=exp((tau(i)/sigma(i)-mu0-b(i)*log(s(i)))/a(i))
    !mu=mu0+a(i)*dlog(vel(i)/vref)+b(i)*dlog(s(i))
    disp(i)=0.d0
    !write(nfile,'(i5,8e15.6)') i,log10(vel(i)),tau(i),sigma(i),disp(i),&
    !& s(i),tau(i)/sigma(i),xcol(i),ycol(i)
    !write(36,'(i6,11e15.6)') i,xcol(i),ycol(i),ang(i)*180/pi,tau(i),sigma(i),tau(i)/sigma(i),&
    !& taudot(i),sigdot(i),vel(i),s(i),vel(i)*s(i)
   end do
   !close(36)

   !output setting
   if(my_rank.eq.0) then
     open(52,file='monitor.dat')
   end if

   !nfmax=NCELLg/NCELL+10

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!-------start time integration--------
write(*,*) 'start time integration'
   x=0.d0
   dtnxt = dtinit
   time1= MPI_Wtime()

   !outside time loop
    do k=1,NSTEP1
       do i=1,NCELL
          y(3*i-2) = dlog(vel(i))
          y(3*i-1) = tau(i)
          y(3*i) = sigma(i)
       end do

       dttry = dtnxt
       dtime = 0.0

          call derivs( NCELL,NCELLg ,x, y, dydx,a, b, dc, gshear, gnorm,sr, vref, mu0,law,rigid)
          !write(*,*) 'call derivs'
          do i = 1, 3 * NCELL
             yscal(i) = abs( y(i) ) + abs( dttry * dydx(i) ) + tiny
          end do

          call rkqs( NCELL,NCELLg,y,dydx,3 * NCELL, x, dttry, eps,yscal, dtdid, &
          & dtnxt, derivs, a, b, dc,gshear,gnorm, sr, vref, mu0,law,rigid)
          !write(*,*) 'call rkqs'
          dtime = dtime + dtdid

       do i = 1, NCELL
          vel(i) = exp(y(3*i-2))
          tau(i) = y(3*i-1)
          sigma(i) = y(3*i)
          disp(i) = disp(i) + exp( y(3*i-2) ) * dtdid
          s(i)=exp((tau(i)/sigma(i)-mu0-a(i)*dlog(vel(i)/vref))/b(i))
          mu(i)=tau(i)/sigma(i)
       end do

       Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(sigma,NCELL,MPI_REAL8,sigmaG,NCELL,                &
       &     MPI_REAL8,MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(vel,NCELL,MPI_REAL8,velG,NCELL,                &
       &     MPI_REAL8,MPI_COMM_WORLD,ierr)

       if(minval(sigmag).lt.minsig) exit
       if(maxval(velg).gt.velmax) exit

       dttry = dtnxt

       if(my_rank.eq.0) then
       time2= MPI_Wtime()
       write(52,'(i5,f16.4,e16.4,f16.4)')k,x,maxval(velG),time2-time1
       end if

       if(mod(k,10).eq.0) then
         do i=1,NCELL
          write(nfile,'(i5,6e15.6,i7,2f15.4)') i+my_rank*NCELL,log10(vel(i)),tau(i),sigma(i),disp(i),&
          & s(i),tau(i)/sigma(i),k,x,s(i)*vel(i)
         end do
         write(nfile,*)
       end if

       !write(13,*)
       dtime1(k) = dtime
       xtime1(k) = x
  !     write(12,'(2e12.5)') xtime1(k)*tref, dtime1(k)*tref
    end do                    ! k = 1, NSTEP1 : outside time loop end

  !-----final results---------------
    time2= MPI_Wtime()
    !$ time2=omp_get_wtime()
    !write(*,*) 'cpu time', time2-time1

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    Call MPI_FINALIZE(ierr)
    !go to 100

    stop
  contains
  !---------------------------------------------------------------------
    subroutine derivs( NCELL,NCELLg, x, y, dydx,a, b, dc,gshear, gnorm,sr, vref, mu0,law,rigid)
  !---------------------------------------------------------------------
    implicit none
    include 'mpif.h'
    !$use omp_lib
    integer :: NCELL,NCELLg
    real(8) :: x, y(3 * NCELL), dydx(3 * NCELL)
    real(8) :: a(NCELL), b(NCELL), dc(NCELL)
    real(8) :: gshear(NCELL,NCELLg), gnorm(NCELL,NCELLg)
    real(8) :: sr, vref, mu0, s(NCELL)
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),vs,vc,v2
    real(8) :: dlnvdt(NCELL),dtaudt(NCELL),dsigdt(NCELL),veltmpG(NCELLg)
    real(8) :: sum_gs(NCELL), sum_gn(NCELL), tau0, c, c1, c2, c3, arg, rigid, alpha
    character :: law
    integer :: i, j, nc

  !     -- Runge-Kutta copled ODE dy(i)/dx definition --
  !if(my_rank.eq.0) write(*,*) 'call derivs'
    do i = 1, NCELL
       veltmp(i) = exp (y(3*i-2))
       tautmp(i) = y(3*i-1)
       sigmatmp(i) = y(3*i)
    enddo
   !write(*,*) 'maxval(veltmp)',maxval(veltmp)

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_ALLGATHER(Veltmp,NCELL,MPI_REAL8,veltmpG,NCELL,                &
  &     MPI_REAL8,MPI_COMM_WORLD,ierr)

    vs=3.d8
    sum_gs = 0.d0
    sum_gn = 0.d0
   !write(*,*) 'call omp'

    do i = 1, NCELL
       do j = 1, NCELLg
          sum_gn(i) = sum_gn(i) + gnorm(i,j) * veltmpG(j)
          sum_gs(i) = sum_gs(i) + gshear(i,j) * veltmpG(j)
       end do
    end do
    !if(my_rank.eq.0) write(*,*) sum_gs(1)
  !write(*,*) maxval(sum_gs)
  !$omp parallel do
    do i=1,NCELL
       dsigdt(i)=sum_gn(i)+sigdot(i)
       !dsigdt(i)=sum_gn(i)
       sum_gs(i)=sum_gs(i)+taudot(i)

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
   subroutine rkqs( NCELL,NCELLg,y, dydx, n,x, htry, eps, yscal, hdid, hnext, derivs,&
        &   a, b, dc, gshear,gnorm, sr, vref, mu0,law,rigid)
  !---------------------------------------------------------------------
   implicit none
   include 'mpif.h'
   integer :: NCELL,NCELLg
   real(8) :: a(NCELL), b(NCELL), dc(NCELL)
   real(8) :: gshear(NCELL,NCELLg),gnorm(NCELL,NCELLg), sr, vref, mu0,rigid
   character*128::law

   INTEGER :: n,nmax
   REAL(8) :: eps,hdid,hnext,htry,x,htemp,err
   real(8) :: dydx(3 * NCELL), y(3 * NCELL), yscal(3 * NCELL)
   EXTERNAL derivs
   PARAMETER  (NMAX = 30000)
   INTEGER :: i
   REAL(8) :: errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON,hmax,errmax_gb
   PARAMETER (SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89e-4)
   parameter(hmax=1d5)
   h=htry
   !if(my_rank.eq.0) write(*,*) 'call rkqs'
   do while(.true.)
       call rkck( NCELL,NCELLg,y, dydx, n, x, h, ytemp, yerr, derivs,a, b, dc,&
   & gshear,gnorm, sr, vref, mu0,law,rigid)

   errmax=0d0
   Do i=1,3*NCELL
     err=dabs(yerr(i)/yscal(i))
     errmax=max(err,errmax)
   EndDo
   !write(*,*) errmax

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,                  &
  &     MPI_MAX,MPI_COMM_WORLD,ierr)

   errmax=errmax_gb/eps

   !write(*,*) errmax
   if(errmax.lt.1.d0) exit
   htemp=SAFETY*h*(errmax**PSHRNK)
   h=sign(max(abs(htemp),0.1*abs(h)),h)
   xnew=x+h
   !write(*,*) h,errmax
   if(xnew-x<1.d-8) stop
   end do

   !if(errmax.gt.ERRCON) then
     hnext=min(1.5*h,SAFETY*h*(errmax**PGROW),3d6)

     hdid=h
     x=x+h
     y(:)=ytemp(:)
   return
   end subroutine

  !---------------------------------------------------------------------
   subroutine rkck( NCELL,NCELLg,y, dydx, n, x, h, yout, yerr, derivs,&
   &     a, b, dc, gshear,gnorm, sr, vref, mu0,law,rigid)
  !---------------------------------------------------------------------
   implicit none
   include 'mpif.h'
   integer :: NCELL,NCELLg
   real(8) :: a(NCELL), b(NCELL), dc(NCELL)
   real(8) :: gshear(NCELL,NCELLg),gnorm(NCELL,NCELLg), sr, vref, mu0,rigid
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
  !if(my_rank.eq.0) write(*,*) 'call rkck'
  !$omp parallel do
   do i=1,n
     ytemp(i)=y(i)+B21*h*dydx(i)
     !write(*,*) i,ytemp(i)
   end do
   !stop
   !write(*,*) maxval(ytemp)

  !    -- 2nd step --

   call derivs(NCELL,NCELLg,x+A2*h,ytemp,ak2, a, b, dc, gshear, gnorm, sr, vref, mu0,law ,rigid)
   !$omp parallel do
   do i=1,n
     ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
     !write(*,*) i,ytemp(i)
   end do

   !write(*,*) maxval(ytemp),maxval(ak2)

  !     -- 3rd step --
   call derivs(NCELL,NCELLg,x+A3*h,ytemp,ak3, a, b, dc,gshear,gnorm, sr, vref, mu0,law,rigid )
   !write(*,*) '3rd step',maxval(ak3),maxval(ak2),maxval(dydx),maxval(y)
   !$omp parallel do
   do i=1,n
     ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
     !write(*,*) i,y(i),ytemp(i),dydx(i),ak2(i),ak3(i)
   end do
   !write(*,*) ytemp
  !stop
  !     -- 4th step --
   call derivs(NCELL,NCELLg,x+A4*h,ytemp,ak4,a, b, dc, gshear,gnorm, sr, vref, mu0,law,rigid )
   !write(*,*) '4th step',maxval(y),maxval(ak4)
  !$omp parallel do
   do i=1,n
     ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
     !write(*,*) i,ytemp(i)
   end do
   !write(*,*) maxval(ytemp)

  !     -- 5th step --
   call derivs(NCELL,NCELLg,x+A5*h,ytemp,ak5,a, b, dc,gshear,gnorm, sr, vref, mu0,law ,rigid)
  !$omp parallel do
   do i=1,n
     ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
   end do
  !     -- 6th step --
   call derivs(NCELL,NCELLg,x+A6*h,ytemp,ak6, a, b, dc,gshear,gnorm, sr, vref, mu0,law,rigid )
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

   subroutine green(gshear,gnorm,NCELLg,NCELL,xcol,ycol,xel,xer,yel,yer,ang,angs,factor)
     implicit none
     integer::i,j,NCELL,NCELLg
     real(8)::gshear(NCELL,NCELLG),gnorm(NCELL,NCELLg),ang(NCELL),angs(NCELLg)
     real(8)::xcol(NCELL),ycol(NCELL),xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg)
     real(8)::dx,dy,angle,sxx,sxy,syy,gtau(2),gsig(2),factor

     do i=1, NCELL
       do j=1, NCELLg
         angle=ang(i)-angs(j) ! check sign
         !write(*,*) angle
          dx=xcol(i)-xel(j)
          dy=ycol(i)-yel(j)
          !write(*,*) dx,dy
          sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
          sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
          syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
          !write(*,*) sxx,sxy,syy
          gtau(1)=-0.5d0*(sxx-syy)*sin(2.d0*angle)+sxy*cos(2.d0*angle)
          gsig(1)=0.5d0*(sxx+syy)-0.5d0*(sxx-syy)*cos(2.d0*angle)-sxy*sin(2.d0*angle)

          dx=xcol(i)-xer(j)
          dy=ycol(i)-yer(j)
          !write(*,*) dx,dy
          sxy=factor*dx*(dx**2-dy**2)/(dx**2+dy**2)**2
          sxx=factor*dy*(3.d0*dx**2+dy**2)/(dx**2+dy**2)**2
          syy=factor*dy*(dy**2-dx**2)/(dy**2+dx**2)**2
          !write(*,*) sxx,sxy,syy
          gtau(2)=-0.5d0*(sxx-syy)*sin(2.d0*angle)+sxy*cos(2.d0*angle)
          gsig(2)=0.5d0*(sxx+syy)-0.5d0*(sxx-syy)*cos(2.d0*angle)-sxy*sin(2.d0*angle)

        gshear(i,j)=gtau(2)-gtau(1)
        gnorm(i,j)=gsig(2)-gsig(1)
       end do
     end do
   end subroutine
  end program
