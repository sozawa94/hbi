program main
  !BIEM on a 3D planar fault.
  !RSF with spatial heterogeneous cut-off velocity
  !H-matrix is used for fast computation.
  !Rate and State Friction (aging law, slip law, Nagata law)
  !length is normalized by Dc/theta*=10^(-5)m
  !time is normalized by theta*=1(sec)
  !stress is normalized by sigma0=100MPa

  !developed by SO OZAWA (Master Student, Earthquake Research Institute, UTokyo)
  !$ use omp_lib
  use m_HACApK_solve
  use m_HACApK_base
  use m_HACApK_use
  use dtriangular
  implicit none
  include 'mpif.h'
  integer::NCELL, nstep1, lp, i,j,k,m,counts,interval,number,lrtrn,nl,NCELLg
  integer::clock,cr,counts2,imax,jmax,NCELLm,seedsize,icomm,np,ierr,my_rank
  integer,allocatable::seed(:)
  character*128::fname,dum,law
  real(8)::a0,b0,sr,omega,theta,dtau,tiny,x,time1,time2,moment,aslip,avv
  real(8)::factor,dx,dy,angle,sxy,sxx,syy,gtau(2),gsig(2),rigid,vc0,dz
  real(8)::c,d,dr,gx,gy,gt,gn,r(3),avel,eps,adisp,vpl,w,rtol,av,outv
  real(8)::dtime,dtnxt,dttry,dtdid,alpha,ds,amp,mui,strinit,velinit,velmax
  real(8),parameter::pi=4.d0*atan(1.d0),pois=0.25d0
  real(8),parameter::dc0=1.d0,sigma0=1.0d0,mu0=0.60d0,tref=1.d9
  real(8),parameter::dtinit=1000d0,vref = 1.d0,minsig=0.05,vs=3.d8
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_calc_entry) :: st_bemv
  real(8),allocatable :: coord(:,:)
  real(8),allocatable::a(:),b(:),dc(:),sigma(:),taudot(:),vc(:),vcg(:)
  real(8),allocatable::vel(:),tau(:),disp(:),mu(:),s(:)
  real(8),allocatable::velG(:),dispG(:),muG(:)
  real(8),allocatable::xcol(:),ycol(:),zcol(:),xe(:,:),ye(:,:),ze(:,:)
  real(8),allocatable::y(:),yscal(:),dydx(:),ycoll(:),zcoll(:)
  real(8),allocatable::xr(:),yr(:),zr(:),ai(:,:),yc(:),zc(:)
  integer::r1,r2,r3,NVER,amari,out,kmax
  integer,allocatable::displs(:),rcounts(:),bi(:,:)

  icomm=MPI_COMM_WORLD
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )
  if(my_rank.eq.0) time1=MPI_Wtime()


  !input parameters
  open(33,file='input.dat')
  read(33,*) dum,NL !length
  read(33,*) dum,NVER
  read(33,*) dum,kmax !number of patches
  read(33,*) dum,nstep1 !maxmimum time step
  read(33,*) dum,number !sample number
  read(33,*) dum,interval !output interval
  read(33,*) dum,velinit !initial slip velocity
  read(33,*) dum,omega !initial omega=V*theta
  read(33,*) dum,ds !radius(normalized by Dc)
  read(33,*) dum,a0 !a in RSF
  read(33,*) dum,b0 !b in RSF
  read(33,*) dum,vc0 !Vc in RSF
  read(33,*) dum,alpha !alpha in RSF
  read(33,*) dum,sr !loading rate
  read(33,*) dum,velmax !stop when Vmax exceeds this value
  read(33,*) dum,rigid !rigidity normalized by 100MPa
  read(33,*) dum,eps !error allowance
  read(33,*) dum,law ! evolution law
  !read(*,*) omega
  close(33)
  allocate(xr(NVER),yr(NVER),zr(NVER),ai(NL,3),bi(NL,3))
  allocate(yc(kmax),zc(kmax))

  !write(*,*) strinit,velinit,sr,law

  !MPI setting
  !NCELLg=2*NL*NL
  allocate(rcounts(np),displs(np+1))
  NCELLg=NL
  amari=mod(NCELLg,np)
  do k=1,amari
    rcounts(k)=NCELLg/np+1
  end do
  do k=amari+1,np
    rcounts(k)=NCELLg/np
  end do
  displs(1)=0
  do k=2,np+1
    displs(k)=displs(k-1)+rcounts(k-1)
  end do
  NCELL=rcounts(my_rank+1)
  if(my_rank.eq.0) write(*,*) rcounts,displs

  allocate(a(NCELL),b(NCELL),dc(NCELL),vc(NCELL),sigma(NCELL),taudot(NCELL))
  allocate(vel(NCELL),tau(NCELL),mu(NCELL),s(NCELL),disp(NCELL))
  allocate(velG(NCELLg),dispG(NCELLg),muG(NCELLg),vcg(NCELLg))
  allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg),xe(NCELLg,3),ye(NCELLg,3),ze(NCELLg,3))
  allocate(ycoll(NCELL),zcoll(NCELL))
  allocate(y(2*NCELL),yscal(2*NCELL),dydx(2 * NCELL))

  !mesh generation
  if(my_rank.eq.0) then
    !call mesh(NL,xe,ye,ze,xcol,ycol,zcol)
    !fault's geometry
  open(15,file='mesh005.dat')
  open(16,file='tri005.dat')

  do i=1,NVER
    read(15,*) yr(i),zr(i)
    call random_number(r)
    yr(i)=yr(i)*ds
    zr(i)=zr(i)*ds
    xr(i)=r(1)*1d-10
    !xe(i)=0.d0
  end do
  do k=1,NCELLg
    read(16,*) ai(k,1:3)
    bi(k,1:3)=int(ai(k,1:3))
    !write(*,*) a(k,1:3),b(k,1:3)
  end do
  close(15)
  close(16)

  !The calculation point is the center of each triangle
  do i=1,NCELLg
    r1=bi(i,1)
    r2=bi(i,2)
    r3=bi(i,3)
    xe(i,:)=(/xr(r1),xr(r2),xr(r3)/)
    ye(i,:)=(/yr(r1),yr(r2),yr(r3)/)
    ze(i,:)=(/zr(r1),zr(r2),zr(r3)/)
    xcol(i)=(xr(r1)+xr(r2)+xr(r3))/3.d0
    ycol(i)=(yr(r1)+yr(r2)+yr(r3))/3.d0
    zcol(i)=(zr(r1)+zr(r2)+zr(r3))/3.d0
  end do
  write(*,*) 'mesh generated'

  do k=1,kmax
    call random_number(r)
    yc(k)=ds*(2*r(1)-1.d0)
    zc(k)=ds*(2*r(2)-1.d0)
    write(*,*) yc(k),zc(k)
  end do
  dr=0.1d0*ds

  do i=1,NCELLg
    vcg(i)=vc0
    do k=1,kmax
      dy=ycol(i)-yc(k);dz=zcol(i)-zc(k)
      if(dy**2+dz**2.lt.dr**2) vcg(i)=vc0*1d4
    end do
  end do
  end if
  !call MPI_FINALIZE(ierr)
  !stop
  !stop

  !MPI communication
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xcol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ycol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(zcol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(xe, 3*NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ye, 3*NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ze, 3*NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(yc, kmax, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(zc, kmax, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERv(ycol,rcounts,displs,MPI_REAL8,ycoll,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(zcol,rcounts,displs,MPI_REAL8,zcoll,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(vcg,rcounts,displs,MPI_REAL8,vc,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  st_bemv%v='s'

  !HACApK setting
  lrtrn=HACApK_init(NCELLg,st_ctl,st_bemv,icomm)
  allocate(coord(NCELLg,3))
  allocate(st_bemv%xcol(NCELLg),st_bemv%ycol(NCELLg),st_bemv%zcol(NCELLg))
  allocate(st_bemv%xe(NCELLg,3),st_bemv%ye(NCELLg,3),st_bemv%ze(NCELLg,3))
  st_bemv%xcol=xcol
  st_bemv%ycol=ycol
  st_bemv%zcol=zcol
  st_bemv%xe=xe
  st_bemv%ye=ye
  st_bemv%ze=ze

  do i=1,NCELLg
    coord(i,1)=xcol(i)
    coord(i,2)=ycol(i)
    coord(i,3)=zcol(i)
  end do

  !generate kernel (Hmatrix is used)
  lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,1d-4)
  if(my_rank.eq.0) write(*,*) 'H-matrix generated'
  !stop


  !frictional parameters and stressing rate (uniform in this simulation)
  !yc=0.5d0*ds;zc=0.5d0*ds;dr=0.1d0*ds
  do i=1,NCELL
    a(i)=a0
    b(i)=b0
    dc(i)=dc0
    !vc(i)=vc0
    ! do k=1,kmax
    !   dy=ycoll(i)-yc(k);dz=zcoll(i)-zc(k)
    !   if(dy**2+dz**2.lt.dr**2) vc(i)=vc0*1d4
    ! end do
    taudot(i)=sr
  end do

  !output setting
  if(my_rank.eq.0) then
    write(fname,'("monitor",i0,".dat")') number
    open(52,file=fname)
    write(fname,'("output/",i0,".dat")') number
    open(50,file=fname)
  end if

  !nfmax=NCELL/NCELL+10
  !write(*,*) 'prepared'

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !initial condition
  theta=omega/velinit
  !mui=mu0+a0*dlog(velinit/vref)+b0*dlog(theta)

  do i=1,NCELL
    mui=mu0+a(i)*dlog(velinit/vref)+b(i)*dlog(theta*vref/dc(i)+vref/vc(i))
    sigma(i)=sigma0
    tau(i)=mui*sigma(i)
    !s(i)=theta
    vel(i)=velinit
    !vel(i)=exp((tau(i)/sigma(i)-mu0-b(i)*log(s(i)))/a(i))
    disp(i)=0.d0
    !write(16,'(i6,11e15.6)') i,xcol(i),ycol(i),ang(i)*180/pi,tau(i),sigma(i),tau(i)/sigma(i),&
    !& taudot(i),sigdot(i),vel(i),s(i),vel(i)*s(i)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_ALLGATHERv(vc,NCELL,MPI_REAL8,vcG,rcounts,displs,                &
  &     MPI_REAL8,MPI_COMM_WORLD,ierr)

  if(my_rank.eq.0) then
    open(19,file='initial.dat')
    do i=1,NCELLg
      write(19,*) ycol(i),zcol(i),vcg(i)
    end do
    close(19)
    write(*,*) 'start time integration'
  end if


  x=0.d0
  dtnxt = dtinit
  outv=1d-6

  !time integration
  do k=1,NSTEP1
    dttry = dtnxt
    do i=1,NCELL
      y(2*i-1) = dlog(vel(i))
      y(2*i) = tau(i)
    end do
    !dtime = 0.0
    call derivs(x, y, dydx,&
    & NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !call derivs( NCELL, x, y, dydx,a, b, dc, sigma, gshear, gnorm,sr, vref, mu0,law,rigid)
    !write(*,*) 'call derivs'
    do i = 1, 2*NCELL
      yscal(i) = abs(y(i)) + abs( dttry * dydx(i) ) + tiny
    end do
    !yscal=y

    !call rkqs( NCELL, y, dydx, 3 * NCELL, x, dttry, eps,yscal, dtdid, &
    !& dtnxt, derivs, a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid)
    call rkqs(y,dydx,x,dttry,eps,yscal,dtdid,dtnxt,derivs)
    !dtime = dtime + dtdid

    !end do

    do i = 1, NCELL
      vel(i) = exp(y(2*i-1))
      tau(i) = y(2*i)
      disp(i) = disp(i) + exp( y(2*i-1) ) * dtdid
      !s(i)=exp((tau(i)/sigma(i)-mu0-a(i)*dlog(vel(i)/vref))/b(i))
      mu(i)=tau(i)/sigma(i)
    end do

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,                &
    &     MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(mu,NCELL,MPI_REAL8,muG,rcounts,displs,                &
    &     MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(disp,NCELL,MPI_REAL8,dispG,rcounts,displs,                &
    &     MPI_REAL8,MPI_COMM_WORLD,ierr)
    !if(my_rank.eq.0) write(*,*) 'allgatherv'

    !output
    !$ time2=omp_get_wtime()
    if(my_rank.eq.0) then
      time2= MPI_Wtime()
      write(52,'(i6,f16.4,3e16.4,f16.4)')k,x,maxval(velG),sum(dispG)/NCELLg,sum(muG)/NCELLg,time2-time1
      !if(mod(k,interval).eq.0) then
      !output control
      out=1
      !A : iteration number
      if(mod(k,interval).eq.0) out=0

      !B : slip velocity
      if(maxval(velG).gt.outv) then
       out=0
       outv=outv*(10.d0)**(0.2d0)
      end if


        !$ time2=omp_get_wtime()
      if(out.eq.0) then
        do i=1,NCELLg
          write(50,'(5e15.6,i7)') ycol(i),zcol(i),log10(velG(i)),muG(i),dispG(i),k
        end do
        write(50,*)
        write(50,*)
      end if

      !  write(*,'(i10,f17.5,2e14.5)') k,x,maxval(velG),time2-time1
        outv=outv*(10.d0)**(0.2d0)

      if(abs(maxval(velG)).gt.velmax) then
        if(my_rank .eq. 0) write(*,*) 'nucleation end'
        exit
      end if

    end if
    dttry = dtnxt
    !dtime1(k) = dtime
    !xtime1(k) = x
    !     write(12,'(2e12.5)') xtime1(k)*tref, dtime1(k)*tref
  end do

  !finalize
  !$ time2=omp_get_wtime()
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp)
  lrtrn=HACApK_finalize(st_ctl)
  if(my_rank.eq.0) then
  time2= MPI_Wtime()
  write(*,*) 'time(s)', time2-time1
  end if
  Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  Call MPI_FINALIZE(ierr)
  !go to 100

  stop
contains
  !---------------------------------------------------------------------
  subroutine derivs(x, y, dydx,&
    &NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    use dtriangular
    implicit none
    include 'mpif.h'
    !$use omp_lib
    integer,intent(in) :: NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in) :: x
    real(8),intent(in) ::y(:)
    real(8),intent(out) :: dydx(:)
    real(8),intent(in) :: a(:), b(:), dc(:),taudot(:),sigma(:),vc(:)
    type(st_HACApK_lcontrol) :: st_ctl
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    type(st_HACApK_calc_entry) :: st_bemv
    real(8),intent(in) :: vref, mu0,vs,rigid,alpha
    character(128),intent(in) :: law
    real(8) :: veltmp(NCELL),tautmp(NCELL)
    real(8) :: dlnvdt(NCELL),dtaudt(NCELL)
    real(8) :: sum_gs(NCELL),sum_gsG(NCELLg),veltmpG(NCELLg)
    real(8):: c1, c2, c3, arg,c
    integer :: i, j, nc

    !     -- Runge-Kutta copled ODE dy(i)/dx definition --
    do i = 1, NCELL
      veltmp(i) = exp(y(2*i-1))
      veldft(i)=veltmp(i)-1.d-4
      tautmp(i) = y(2*i)
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !call MPI_ALLGATHERv(Veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,                &
    !&     MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(Veldft,NCELL,MPI_REAL8,veltmpG,rcounts,displs,                &
    &     MPI_REAL8,MPI_COMM_WORLD,ierr)
    !if(my_rank.eq.0) write(*,*) 'allgatherv'

    lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,sum_gsG,veltmpG)
    !if(my_rank.eq.0) write(*,*) sum_gsG

    call MPI_SCATTERv(sum_gsG,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !if(my_rank.eq.0) write(*,*) 'scatterv'

    !$omp parallel do
    do i=1,NCELL
!      sum_gs(i)=sum_gs(i)+taudot(i)
      !     -- equation dtau/dt --
      select case(law)
        ! aging law (Linker & Dieterich, 1992)
      case('a')
        arg=dc(i)/vref*(exp((tautmp(i)/sigma(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
        dlnvdt(i)=sum_gs(i)-b(i)*sigma(i)/(arg+dc(i)/vc(i))*(1.d0-veltmp(i)*arg/dc(i))
        dlnvdt(i)=dlnvdt(i)/(a(i)*sigma(i)+0.5d0*rigid*veltmp(i)/vs)
        dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)

        !c1 = b(i) / dc(i) * veltmp(i)
        !arg = tautmp(i)/sigma(i) - mu0 - a(i) * dlog( veltmp(i) / vref)
        !c2 = b(i) / dc(i) * vref * exp( - arg / b(i) ) !b/theta
        !c3 = alpha/sigmatmp(i)*dsigdt(i)

        !dtaudt(i)=sum_gs(i) + 0.5d0*rigid/vs*veltmp(i)/a(i)*((c2-c1))
        !dtaudt( i ) = dtaudt(i)/(1.d0+rigid*veltmp(i)*0.5d0/vs/a(i)/sigma(i))
        !dtaudt( i ) = (sr+sum_gs)*a(i)/veltmp(i)+rigid*0.5d0/vs*(c2-c1)
        !dlnvdt(i)=2.d0*vs/rigid/veltmp(i)*(dtaudt(i)-sum_gs)
        !dlnvdt(i)=(dtaudt(i)/sigma(i)+c1-c2+c3)/a(i)
        !dlnvdt( i ) = ( dtaudt(i) /sigma(i) + c1 - c2) / a(i)

        ! slip law
      case('s')
        !if(vel(i).gt.1d4) dc(i)=
        !arg = tautmp(i)/sigma(i) - mu0 - (a(i)-b(i)) * dlog( veltmp(i) / vref)
        arg=dc(i)/vref*(exp((tautmp(i)/sigma(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
        !c1 = veltmp(i) / dc(i) * arg
        dlnvdt(i)=sum_gs(i)+b(i)*sigma(i)*veltmp(i)*arg/(arg*dc(i)+dc(i)**2/vc(i))*dlog(veltmp(i)*arg/dc(i))
        dlnvdt(i)=dlnvdt(i)/(a(i)*sigma(i)+0.5d0*rigid*veltmp(i)/vs)
        dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
        !c3 = alpha/sigmatmp(i)*dsigdt(i)
        !dtaudt(i)=sum_gs(i) + 0.5d0*rigid/vs*veltmp(i)/a(i)*(-c1)
        !dtaudt( i ) = dtaudt(i)/(1.d0+rigid*veltmp(i)*0.5d0/vs/a(i)/sigma(i))
        !dtaudt( i ) = (sr+sum_gs)*a(i)/veltmp(i)-rigid*0.5d0/vs*c1
        !dtaudt( i ) = dtaudt(i)/(rigid*0.5d0/vs+a(i)/veltmp(i))
        !dlnvdt( i ) = ( dtaudt(i) /sigma(i) +c1) / a(i)
        !dlnvdt( i ) = ( dtaudt(i) / sigma(i)- + c1 )/  a(i)

        !composite law
        !case('c')
        !  dtaudt(i)=sum_gs(i)
        !  vc=1d-2
        !  arg = tautmp(i) - tau0 - (a(i)-b(i)) * dlog( veltmp(i) / vref)
        !  s(i)=exp((tautmp(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))
        !  c1 = veltmp(i) / dc(i) * arg-b(i)/s(i)*exp(-veltmp(i)/vc)
        !  dlnvdt( i ) = ( dtaudt(i)  + c1 ) / a(i)

        !nagata law
      case('n')
        c=2.0
        c1 = b(i) / dc(i) * veltmp(i)
        arg = tautmp(i)/sigma(i) - mu0 - a(i) * dlog( veltmp(i) / vref)
        c2 = b(i) / dc(i) * vref * exp( - arg / b(i) ) !b/theta
        !c3 = alpha/sigmatmp(i)*dsigdt(i)
        dtaudt(i)=sum_gs(i)
        dlnvdt(i)=(dtaudt(i)*sigma(i))/sigma(i)**2- &
        & c2+c1+c*dtaudt(i)
        dlnvdt(i)=dlnvdt(i)/a(i)
      end select
      !write(*,*) dlnvdt(1000)
    enddo
    !write(*,*) maxval(dlnvdt),maxval(dtaudt),maxval(dsigdt)

    do i = 1, NCELL
      dydx(2*i-1) = dlnvdt( i )
      dydx(2*i) = dtaudt( i )
    enddo

    return
  end subroutine derivs
  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
    !---------------------------------------------------------------------
    implicit none
    include 'mpif.h'
    real(8),intent(in)::dydx(:),yscal(:),htry,eps
    real(8),intent(inout)::y(:),x
    real(8),intent(out)::hdid,hnext

    interface
      !subroutine derivs(x,y,dydx)
      subroutine derivs(x, y, dydx,&
        &NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
        use m_HACApK_solve
        use m_HACApK_base
        use m_HACApK_use
        use dtriangular
        implicit none
        include 'mpif.h'
        real(8),intent(in)::x
        real(8),intent(in)::y(:)
        real(8),intent(in)::a(:),b(:),dc(:),taudot(:),sigma(:),vc(:)
        real(8),intent(in)::rigid,vref,mu0,vs,alpha
        integer,intent(in)::NCELL,NCELLg,rcounts(:),displs(:)
        type(st_HACApK_lcontrol) :: st_ctl
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_calc_entry) :: st_bemv
        real(8),intent(out)::dydx(:)
        character(128),intent(in)::law
      end subroutine derivs
    end interface
    INTEGER :: i,ierr,iter
    REAL(8) :: errmax,h,xnew,htemp,errmax_gb
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89e-4,hmax=1d5

    h=htry
    iter=0
    !write(*,*) 'call rkqs'
    do while(.true.)
      iter=iter+1
      !write(*,*) 'call rkck'
      call rkck(y,dydx,x,h,ytemp,yerr,derivs)
      errmax=maxval(abs(yerr(:)/yscal(:)))/eps

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,                  &
      &     MPI_MAX,MPI_COMM_WORLD,ierr)
      !if(my_rank.eq.0) write(52,*) iter,h,errmax_gb
      !write(*,*) errmax
      if(errmax_gb.lt.1.d0) exit
      htemp=SAFETY*h*(errmax_gb**PSHRNK)
      h=sign(max(abs(htemp),0.1*abs(h)),h)
      xnew=x+h
      !write(*,*) h,errmax
      if(xnew-x<1d-8) stop
      !    write(*,*) x,k
      !    stop
      !  end if
    end do

    !if(errmax.gt.ERRCON) then
    hnext=min(3*h,SAFETY*h*(errmax_gb**PGROW),1e7)

    hdid=h
    x=x+h
    y(:)=ytemp(:)
    return
  end subroutine

  !---------------------------------------------------------------------
  subroutine rkck(y,dydx,x,h,yout,yerr,derivs)
    !---------------------------------------------------------------------
    implicit none
    include 'mpif.h'
    real(8),intent(in)::y(:),dydx(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    interface
      !subroutine derivs(x,y,dydx)
      subroutine derivs(x, y, dydx,&
        &NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
        use m_HACApK_solve
        use m_HACApK_base
        use m_HACApK_use
        use dtriangular
        implicit none
        include 'mpif.h'
        real(8),intent(in)::x
        real(8),intent(in)::y(:)
        real(8),intent(in)::a(:),b(:),dc(:),taudot(:),sigma(:),vc(:)
        real(8),intent(in)::rigid,vref,mu0,vs,alpha
        integer,intent(in)::NCELL,NCELLg,rcounts(:),displs(:)
        type(st_HACApK_lcontrol) :: st_ctl
        type(st_HACApK_leafmtxp) :: st_leafmtxp
        type(st_HACApK_calc_entry) :: st_bemv
        real(8),intent(out)::dydx(:)
        character(128),intent(in)::law
      end subroutine derivs
    end interface
    INTEGER ::i
    integer,parameter::nmax=100000
    REAL(8) :: ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax),ytemp(nmax)
    real(8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51
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
    do i=1,size(y)
      ytemp(i)=y(i)+B21*h*dydx(i)
      !write(*,*) i,ytemp(i)
    end do
    !stop
    !write(*,*) maxval(ytemp)

    !    -- 2nd step --
    call derivs(x+a2*h, ytemp, ak2,&
    & NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !call derivs(NCELL,x+A2*h,ytemp,ak2, a, b, dc, sigma, gshear, gnorm, sr, vref, mu0,law ,rigid)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      !write(*,*) i,ytemp(i)
    end do

    !write(*,*) maxval(ytemp),maxval(ak2)

    !     -- 3rd step --
    call derivs(x+a3*h, ytemp, ak3,&
    & NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !call derivs(NCELL,x+A3*h,ytemp,ak3, a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid )
    !write(*,*) '3rd step',maxval(ak3),maxval(ak2),maxval(dydx),maxval(y)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      !write(*,*) i,y(i),ytemp(i),dydx(i),ak2(i),ak3(i)
    end do
    !write(*,*) ytemp
    !stop
    !     -- 4th step --
    call derivs(x+a4*h, ytemp, ak4,&
    & NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !call derivs(NCELL,x+A4*h,ytemp,ak4,a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid )
    !write(*,*) '4th step',maxval(y),maxval(ak4)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
      !write(*,*) i,ytemp(i)
    end do
    !write(*,*) maxval(ytemp)

    !     -- 5th step --
    call derivs(x+a5*h, ytemp, ak5,&
    & NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !call derivs(NCELL,x+A5*h,ytemp,ak5,a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law ,rigid)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do
    !     -- 6th step --
    call derivs(x+a6*h, ytemp, ak6,&
    & NCELL,NCELLg,rcounts,displs,a,b,dc,vc,st_leafmtxp,st_bemv,st_ctl,taudot,vref,mu0,law,rigid,vs,alpha,sigma)
    !call derivs(NCELL,x+A6*h,ytemp,ak6, a, b, dc, sigma, gshear,gnorm, sr, vref, mu0,law,rigid )
    !$omp parallel do
    do i=1,size(y)
      yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
    end do

    !$omp parallel do
    do i=1,size(y)
      yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
    end do
    !write(*,*) yerr

    return
  end subroutine
  subroutine mesh(NL,xe,ye,ze,xcol,ycol,zcol)
    implicit none
    integer,intent(in)::NL
    real(8),intent(out)::xcol(:),ycol(:),zcol(:),xe(:,:),ye(:,:),ze(:,:)
    integer::i,j,k
    real(8)::r(3),c
    c=1.d0/3.d0
    do i=1,NL/2
      do j=1,NL/2
        call random_number(r)
        k=(i-1)*NL*4+8*j-7
        ye(k,1)=(2*i-2)*ds
        ze(k,1)=(2*j-2)*ds
        ye(k,2)=(2*i-1)*ds
        ze(k,2)=(2*j-2)*ds
        ye(k,3)=(2*i-2)*ds
        ze(k,3)=(2*j-1)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j-6
        ye(k,1)=(2*i-1)*ds
        ze(k,1)=(2*j-1)*ds
        ye(k,2)=(2*i-1)*ds
        ze(k,2)=(2*j-2)*ds
        ye(k,3)=(2*i-2)*ds
        ze(k,3)=(2*j-1)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j-5
        ye(k,1)=(2*i-1)*ds
        ze(k,1)=(2*j-1)*ds
        ye(k,2)=(2*i-1)*ds
        ze(k,2)=(2*j-2)*ds
        ye(k,3)=(2*i)*ds
        ze(k,3)=(2*j-1)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j-4
        ye(k,1)=(2*i)*ds
        ze(k,1)=(2*j-2)*ds
        ye(k,2)=(2*i-1)*ds
        ze(k,2)=(2*j-2)*ds
        ye(k,3)=(2*i)*ds
        ze(k,3)=(2*j-1)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j-3
        ye(k,1)=(2*i-2)*ds
        ze(k,1)=(2*j)*ds
        ye(k,2)=(2*i-2)*ds
        ze(k,2)=(2*j-1)*ds
        ye(k,3)=(2*i-1)*ds
        ze(k,3)=(2*j)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j-2
        ye(k,1)=(2*i-1)*ds
        ze(k,1)=(2*j-1)*ds
        ye(k,2)=(2*i-2)*ds
        ze(k,2)=(2*j-1)*ds
        ye(k,3)=(2*i-1)*ds
        ze(k,3)=(2*j)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j-1
        ye(k,1)=(2*i-1)*ds
        ze(k,1)=(2*j-1)*ds
        ye(k,2)=(2*i)*ds
        ze(k,2)=(2*j-1)*ds
        ye(k,3)=(2*i-1)*ds
        ze(k,3)=(2*j)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

        call random_number(r)
        k=(i-1)*NL*4+8*j
        ye(k,1)=(2*i)*ds
        ze(k,1)=(2*j)*ds
        ye(k,2)=(2*i)*ds
        ze(k,2)=(2*j-1)*ds
        ye(k,3)=(2*i-1)*ds
        ze(k,3)=(2*j)*ds
        xe(k,:)=r*1d-8
        xcol(k)=c*(xe(k,1)+xe(k,2)+xe(k,3))
        ycol(k)=c*(ye(k,1)+ye(k,2)+ye(k,3))
        zcol(k)=c*(ze(k,1)+ze(k,2)+ze(k,3))
        !write(*,*) k,ycol(k),zcol(k)

      end do
    end do
  end subroutine mesh
end program
