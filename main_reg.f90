program main
  !$ use omp_lib
  use m_HACApK_solve
  use m_HACApK_base
  use m_HACApK_use
  !use mod_derivs
  use mod_constant
  use m_HACApK_calc_entry_ij
  implicit none
  include 'mpif.h'
  integer:: date_time(8)
  character(len=10):: sys_time(3)
  integer::NCELL, nstep1, lp, i,j,k,m,counts,interval,number,lrtrn,nl,NCELLg
  integer::clock,cr,counts2,imax,jmax,NCELLm,seedsize,icomm,np,ierr,my_rank,load,eventcount,thec
  logical::slipping,outfield,limitsigma
  integer,allocatable::seed(:)
  character*128::fname,dum,law,input_file,problem,geofile
  real(8)::a0,b0,dc0,sr,omega,theta,dtau,tiny,x,time1,time2,moment,aslip,avv
  real(8)::vc0,mu0,dtinit,onset_time,tr,vw0,fw0,velmin,muinit
  real(8)::r,eps,vpl,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxevent
  real(8)::dtime,dtnxt,dttry,dtdid,dtmin,alpha,ds,amp,mui,strinit,velinit,velmax
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxps,st_leafmtxpn
  type(st_HACApK_calc_entry) :: st_bemv
  !type(t_deriv)::
  real(8),allocatable ::coord(:,:),vmax(:)
  real(8),allocatable::ag(:),bg(:),dcg(:),vcg(:),fwg(:),vwg(:),taudotg(:),sigdotg(:)
  real(8),allocatable::a(:),b(:),dc(:),vc(:),fw(:),vw(:),ac(:),taudot(:),sigdot(:)
  real(8),allocatable::phi(:),vel(:),tau(:),sigma(:),disp(:),mu(:)
  real(8),allocatable::phiG(:),velG(:),tauG(:),sigmaG(:),dispG(:),muG(:),ruptG(:)
  real(8),allocatable::xcol(:),ycol(:),zcol(:)
  real(8),allocatable::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
  real(8),allocatable::xel(:),xer(:),yel(:),yer(:),ang(:)
  real(8),allocatable::y(:),yscal(:),dydx(:),xcoll(:),zcoll(:),yg(:)
  real(8),allocatable::xr(:),yr(:),zr(:),ai(:,:)
  integer::r1,r2,r3,NVER,amari,out,kmax,loci,locj,loc,stat
  integer,allocatable::displs(:),rcounts(:),vmaxin(:),rupsG(:)

  !initialize
  icomm=MPI_COMM_WORLD
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )
  if(my_rank.eq.0) time1=MPI_Wtime()

  !input file must be specified when running
  !example) mpirun -np 16 ./ha.out default.in
  call get_command_argument(1,input_file,status=stat)

  !reading input file
  if(my_rank.eq.0) write(*,*) 'Reading input file'
  open(33,file=input_file)

  !information of fault geometry
  read(33,*) dum,problem

  select case(problem)
  case('2dp')
    read(33,*) dum,NCELLg
  case('2dn')
    read(33,*) dum,NCELLg
    read(33,*) dum,geofile
  case('3dp')
    read(33,*) dum,imax !x length
    read(33,*) dum,jmax !z length
  end select

  !read(33,*) dum,kmax !number of VS patches
  !read(33,*) dum,dr !radius of VS patches (unit:ds)
  read(33,*) dum,ds !mesh interval(normalized by Dc)

  !output control
  read(33,*) dum,number !output filename
  read(33,*) dum,interval !output frequency
  !read(33,*) dum,dlapse !output per time
  !read(33,*) dum,loci !output localdata i
  !read(33,*) dum,locj !output localdata j

  !continue or stop control
  read(33,*) dum,nstep1 !maxmimum time step
  !read(33,*) dum,thec !stop when the eventcount exceeds this value
  read(33,*) dum,velmax !stop when slip rate exceeds this value
  read(33,*) dum,velmin

  !physical parameters in calculation
  read(33,*) dum,law ! evolution law in RSF a: aging s: slip
  read(33,*) dum,a0 !a in RSF
  read(33,*) dum,b0 !b in RSF
  read(33,*) dum,dc0
  read(33,*) dum,vw0
  read(33,*) dum,fw0
  read(33,*) dum,vc0 !cut-off velocity in RSF (should be large enough when assuming conventional RSF)
  read(33,*) dum,mu0 !mu0 in RSF
  read(33,*) dum,load !loading type 0: uniform stress 1: uniform slip deficit
  read(33,*) dum,sr !if load=0: stressing rate (sigma0/sec) load=1: unused
  read(33,*) dum,vpl !if load=1: loading velocity load=0: unused
  read(33,*) dum,tr !viscoelastic relaxation (unused)
  !read(33,*) dum,rigid !shear modulus normalized by sigma0
  !read(33,*) dum,vs !Swave speed (dc/s)
  !read(33,*) dum,pois !poisson ratio

  !initial values
  read(33,*) dum,velinit !initial slip velocity
  read(33,*) dum,muinit !initial omega=V*theta
  read(33,*) dum,dtinit !initial timestep
  read(33,*) dum,limitsigma
  read(33,*) dum,eps !error allowance in time integration in Runge-Kutta
  !read(*,*) amp
  !read(*,*) omega

  !for FDMAP
  !read(33,*) coordinate_file
  close(33)

  !MPI setting
  !NCELLg=2*NL*NL
  if(problem.eq.'3dp') then
    NCELLg=imax*jmax
    loc=loci*(imax-1)+locj
  end if

  allocate(rcounts(np),displs(np+1))
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

  !allocation
  allocate(a(NCELL),b(NCELL),dc(NCELL),vc(NCELL),fw(NCELL),vw(NCELL),taudot(NCELL),sigdot(NCELL))
  allocate(ag(NCELLg),bg(NCELLg),dcg(NCELLg),vcg(NCELLg),fwg(NCELLg),vwg(NCELLg),taudotg(NCELLg),sigdotg(NCELLg))
  allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL))
  allocate(phiG(NCELLg),velG(NCELLg),tauG(NCELLg),sigmaG(NCELLg),dispG(NCELLg),muG(NCELLg),ruptG(NCELLg),rupsG(NCELLg))
  select case(problem)
  case('2dp')
    allocate(xcol(NCELLg),xel(NCELLg),xer(NCELLg))
    xcol=0d0;xel=0d0;xer=0d0
  case('2dn')
    allocate(xcol(NCELLg),ycol(NCELLg),ang(NCELLg))
    allocate(xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg))
    xcol=0d0;ycol=0d0;ang=0d0;xel=0d0;xer=0d0;yel=0d0;yer=0d0
  case('3dp')
    allocate(xcol(NCELLg),zcol(NCELLg))
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
    allocate(xcoll(NCELL),zcoll(NCELL))
    xcol=0d0; zcol=0d0
    xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
    zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
  end select

  select case(problem)
  case('2dp','3dp')
    allocate(y(2*NCELL),yscal(2*NCELL),dydx(2*NCELL))
  case('2dn')
    allocate(y(3*NCELL),yscal(3*NCELL),dydx(3*NCELL))
  !case('2dnv')
  !  allocate(y(4*NCELL),yscal(4*NCELL),dydx(4*NCELL))
  end select

  allocate(vmax(NCELLg),vmaxin(NcELLg))

  !mesh generation (rectangular assumed)
  if(my_rank.eq.0) write(*,*) 'Generating mesh'
  if(my_rank.eq.0) then
    select case(problem)
    case('2dp')
      call coordinate2dp(NCELLg,ds,xel,xer,xcol)
    case('2dn') !geometry file is necessary
      call coordinate2dn(geofile,NCELLg,xel,xer,yel,yer,xcol,ycol,ang)
    case('3dp')
      call coordinate3d(imax,jmax,ds,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
    end select
  end if

  !random number seed
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do i = 1, seedsize
    call system_clock(count=seed(i))
  end do
  call random_seed(put=seed(:))

  write(*,*) 'MPI communication'
  !MPI communication
  select case(problem)
  case('2dp')
    call MPI_BCAST(xcol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xel, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xer, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  case('2dn')
    call MPI_BCAST(xcol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ycol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xel, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xer, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(yel, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(yer, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ang, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  case('3dp')
    call MPI_BCAST(xcol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zcol, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xs1, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xs2, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xs3, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(xs4, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zs1, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zs2, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zs3, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(zs4, NCELLg, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_SCATTERv(xcol,rcounts,displs,MPI_REAL8,xcoll,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_SCATTERv(zcol,rcounts,displs,MPI_REAL8,zcoll,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end select


  !HACApK setting
  lrtrn=HACApK_init(NCELLg,st_ctl,st_bemv,icomm)
  allocate(coord(NCELLg,3))
  select case(problem)
  case('2dp')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%problem=problem

  case('2dn')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg))
    allocate(st_bemv%ycol(NCELLg),st_bemv%yel(NCELLg),st_bemv%yer(NCELLg),st_bemv%ang(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%ycol=ycol;st_bemv%yel=yel;st_bemv%yer=yer
    st_bemv%ang=ang
    st_bemv%problem=problem

  case('3dp')
    allocate(st_bemv%xcol(NCELLg),st_bemv%zcol(NCELLg))
    allocate(st_bemv%xs1(NCELLg),st_bemv%xs2(NCELLg),st_bemv%xs3(NCELLg),st_bemv%xs4(NCELLg))
    allocate(st_bemv%zs1(NCELLg),st_bemv%zs2(NCELLg),st_bemv%zs3(NCELLg),st_bemv%zs4(NCELLg))
    st_bemv%xcol=xcol
    st_bemv%zcol=zcol
    st_bemv%xs1=xs1
    st_bemv%xs2=xs2
    st_bemv%xs3=xs3
    st_bemv%xs4=xs4
    st_bemv%zs1=zs1
    st_bemv%zs2=zs2
    st_bemv%zs3=zs3
    st_bemv%zs4=zs4
    st_bemv%problem=problem

  end select


  !generate kernel (H-matrix aprrox)
  if(my_rank.eq.0) write(*,*) 'Generating kernel'
  select case(problem)
  case('2dp')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=0.d0
      coord(i,3)=0.d0
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,1d-4)

  case('2dn')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=ycol(i)
      coord(i,3)=0.d0
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    st_bemv%v='s'
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,1d-4)
    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxpn,st_bemv,st_ctl,coord,1d-4)

  case('3dp')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=0.d0
      coord(i,3)=zcol(i)
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,1d-4)
  end select

  !setting frictional parameters
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(my_rank.eq.0) write(*,*) 'Setting fault parameters'
  if(my_rank.eq.0) call params(problem,NCELLg,a0,b0,dc0,vc0,ag,bg,dcg,vcg,fwg,vwg)
  if(my_rank.eq.0) call loading(problem,NCELLg,sr,taudotG,sigdotG)

  call MPI_SCATTERv(vcg,rcounts,displs,MPI_REAL8,vc,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(ag,rcounts,displs,MPI_REAL8,a,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(bg,rcounts,displs,MPI_REAL8,b,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(dcg,rcounts,displs,MPI_REAL8,dc,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(fwg,rcounts,displs,MPI_REAL8,fw,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(vwg,rcounts,displs,MPI_REAL8,vw,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(taudotg,rcounts,displs,MPI_REAL8,taudot,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(sigdotg,rcounts,displs,MPI_REAL8,sigdot,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


  !setting initial condition
  call initcond(phiG,sigmaG,tauG)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(tauG,rcounts,displs,MPI_REAL8,tau,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(sigmaG,rcounts,displs,MPI_REAL8,sigma,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(phiG,rcounts,displs,MPI_REAL8,phi,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  !mui=mu0+a0*dlog(velinit/vref)+b0*dlog(theta)
  !xc=0.5d0*imax*ds;zc=0.5d0*jmax*ds;dr=30*ds

  !do i=1,NCELL
    !call random_number(omega)
    !omega=omega*0.5d0+0.5d0
    !theta=omega/velinit
    !s(i)=omega+0.1d0*s(i)
    !mu(i)=mu0+(a(i)-b(i))*dlog(abs(velinit)/vref)
  !  phi(i)=omega
  !  mu(i)=a(i)*asinh(0.5d0*velinit/vref*exp(phi(i)/a(i)))
  !  vel(i)=velinit
  !  sigma(i)=sigma0
  !  tau(i)=mu(i)*sigma(i)*sign(vel(i),1.d0)/vel(i)
    !phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(tau(i)/sigma(i)/a(i)))
    !write(*,*) tau(i),vel(i),phi(i)
  !end do

  !for FDMAP-BIEM simulation
  !call input_from_FDMAP()

  !setting output files
  if(my_rank.eq.0) then
    write(fname,'("output/monitor",i0,".dat")') number
    open(52,file=fname)
    write(fname,'("output/",i0,".dat")') number
    open(50,file=fname)
    write(fname,'("output/rupt",i0,".dat")') number
    open(48,file=fname)
    open(19,file='job.log',position='append')
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


  !record parameters in file
  if(my_rank.eq.0) then
    call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    write(19,*)
    write(19,'(3i5)') date_time(1),date_time(2),date_time(3)
    write(19,*) 'job number',number !output filename
    !add anything you want

    write(*,*) 'start time integration'
  end if

  !setting minimum time step by CFL condition
  !dtmin=0.5d0*ds/(vs*sqrt(3.d0))

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !start time integration
  time1=MPI_Wtime()
  x=0.d0 !x is time
  ruptG=0d0
  rupsG=0
  disp=0.d0
  dtnxt = dtinit
  !outv=1d-6
  slipping=.false.
  eventcount=0
   !call MPI_ALLGATHERv(sigma,NCELL,MPI_REAL8,sigmaG,rcounts,displs, MPI_REAL8,MPI_COMM_WORLD,ierr)
    !call MPI_ALLGATHERv(tau,NCELL,MPI_REAL8,tauG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    !call MPI_ALLGATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    !call MPI_ALLGATHERv(mu,NCELL,MPI_REAL8,muG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    !call MPI_ALLGATHERv(disp,NCELL,MPI_REAL8,dispG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

  !do i=1,NCELLg
  !  write(50,'(8e15.6,i6)') xcol(i),ycol(i),velG(i),tauG(i),sigmaG(i),muG(i),dispG(i),x,k
  !end do
  !write(50,*)

  do k=1,NSTEP1
    dttry = dtnxt

    select case(problem)
    case('2dp','3dp')
      do i=1,NCELL
        y(2*i-1) = phi(i)
        y(2*i) = tau(i)
      end do

    case('2dn')
      do i=1,NCELL
        y(3*i-2) = phi(i)
        y(3*i-1) = tau(i)
        y(3*i)=sigma(i)
      end do
      ! call MPI_SCATTERv(yg,3*rcounts,3*displs,MPI_REAL8,y,3*NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    end select

    call derivs(x, y, dydx)!,,st_leafmtxps,st_leafmtxpn,st_bemv,st_ctl)

    do i = 1, size(yscal)
      yscal(i)=abs(y(i))+abs(dttry*dydx(i))+tiny
    end do

    call rkqs(y,dydx,x,dttry,eps,yscal,dtdid,dtnxt)

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    select case(problem)
    case('2dp','3dp')
      ! call MPI_ALLGATHERv(y,2*NCELL,MPI_REAL8,yG,2*rcounts,2*displs,                &
      ! &     MPI_REAL8,MPI_COMM_WORLD,ierr)
      do i = 1, NCELL
        phi(i) = y(2*i-1)
        tau(i) = y(2*i)
        !disp(i) = disp(i)+exp(y(2*i-1))*dtdid
        vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
        !write(*,*)vel(i),dtdid
        disp(i)=disp(i)+vel(i)*dtdid
        mu(i)=tau(i)/sigma(i)
      end do
    case('2dn')
      ! call MPI_ALLGATHERv(y,3*NCELL,MPI_REAL8,yG,3*rcounts,3*displs,                &
      ! &     MPI_REAL8,MPI_COMM_WORLD,ierr)
      do i = 1, NCELL
        phi(i) = y(3*i-2)
        tau(i) = y(3*i-1)
        sigma(i) = y(3*i)
        !artificial limit of normal stress motivated by plastic simulation
        if(limitsigma) then
          if(sigma(i).lt.30d0) sigma(i)=30d0
          if(sigma(i).gt.170d0) sigma(i)=170d0
        end if
        !disp(i) = disp(i) + exp( y(3*i-2) i)*dtdid
        vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
        disp(i)=disp(i)+vel(i)*dtdid
        !s(i)=a(i)*dlog(2.d0*vref/vel(i)*dsinh(tau(i)/sigma(i)/a(i)))
        !s(i)=exp((tau(i)/sigma(i)-mu0-a(i)*dlog(vel(i)/vref))/b(i))
        mu(i)=tau(i)/sigma(i)
      end do

    end select

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(sigma,NCELL,MPI_REAL8,sigmaG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(tau,NCELL,MPI_REAL8,tauG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(mu,NCELL,MPI_REAL8,muG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(disp,NCELL,MPI_REAL8,dispG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(phi,NCELL,MPI_REAL8,phiG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    !if(my_rank.eq.0) write(*,*) 'allgatherv'

    !output
    if(my_rank.eq.0) then
      time2= MPI_Wtime()
      write(52,'(i6,f18.4,3e16.5,f16.4)')k,x,maxval(log10(abs(velG))),sum(dispG)/NCELLg,sum(muG)/NCELLg,time2-time1
      !do i=1,size(vmax)
      !  vmax(i)=max(vmax(i),velG(i))
      !  vmaxin(i)=k
      !end do
      !if(mod(k,interval).eq.0) then

      !for FDMAP
      !PsiG=ag*dlog(2.d0*vref/velG*dsinh(tauG/sigmaG/ag))
      do i=1,NCELLg
        if(velG(i).gt.0.05d0.and.ruptG(i).le.1d-6) then
          ruptG(i)=x
          rupsG(i)=k
        end if
      end do


      !output distribution control
      outfield=.false.
      !A : iteration number
      if(mod(k,interval).eq.0) outfield=.true.
      !if(k.lt.18000) out=1

      !B : slip velocity
      !if(maxval(velG).gt.outv) then
      !  out=0
      !  outv=outv*(10.d0)**(0.5d0)
      !end if

      if(outfield) then
        write(*,*) 'time step=' ,k
        select case(problem)
        case('3dp')
          do i=1,NCELLg
            write(50,'(5e15.6,i7)') xcol(i),zcol(i),log10(velG(i)),muG(i),dispG(i),k
          end do
          write(50,*)
          write(50,*)
        case('2dp','2dn')
          do i=1,NCELLg
            write(50,'(i0,8e15.6,i6)') i,xcol(i),log10(abs(velG(i))),tauG(i),sigmaG(i),muG(i),dispG(i),phiG(i),x,k
          end do
          write(50,*)
        end select
      end if

    end if


    ! if(.not.slipping) then
    !   vmaxevent=0.d0
    !   if(abs(maxval(velG)).gt.1d-5) then
    !     slipping=.true.
    !     idisp=dispG
    !     onset_time=x
    !     lapse=0.d0
    !     if(my_rank.eq.0) write(fname,'("output/event",i0,".dat")') number
    !     if(my_rank.eq.0) open(53,file=fname)
    !   end if
    ! end if
    !
    ! if(slipping) then
    !   if(abs(maxval(velG)).lt.1d-5) then
    !     slipping=.false.
    !     eventcount=eventcount+1
    !     if(my_rank.eq.0) write(19,*) eventcount,vmaxevent
    !   end if
    !   vmaxevent=max(vmaxevent,maxval(velG))
    !   !write(53,'(i6,4e16.6)') !k,x-onset_time,sum(dispG-idisp),sum(velG),sum(acg**2)
    !   !if(x-onset_time.gt.lapse) then
    !   !  lapse=lapse+dlapse
    !   !end if
    ! end if

    !simulation ends before nstep1 when
    !(1) eventcount exceeds threshold
    ! if(eventcount.eq.thec) then
    !   if(my_rank .eq. 0) write(*,*) 'eventcount 10'
    !   go to 200
    ! end if

    !(2) slip velocity exceeds threshold (for nucleation)
    if(maxval(abs(velG)).gt.velmax) then
      if(my_rank .eq. 0) write(*,*) 'slip rate above threshold'
      exit
    end if

    if(maxval(abs(velG)).lt.velmin) then
      if(my_rank .eq. 0) write(*,*) 'slip rate below threshold'
      exit
    end if

    dttry = dtnxt
  end do


  !output for FDMAP communication
  !call output_to_FDMAP()
  if(my_rank.eq.0) then
    do i=1,NCELLg
      write(48,'(3f16.4,i6)') xcol(i),ycol(i),ruptG(i),rupsG(i)
    end do
  end if

  200  if(my_rank.eq.0) then
  time2= MPI_Wtime()
  write(*,*) 'time(s)', time2-time1
  close(52)
  close(50)
  close(19)
end if
Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
select case(problem)
case('2dp','3dp')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxps)
case('2dn')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxps)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxpn)
end select
lrtrn=HACApK_finalize(st_ctl)
Call MPI_FINALIZE(ierr)
stop
contains
  subroutine initcond(phiG,sigmaG,tauG)
    implicit none
    real(8),intent(out)::phiG(:),sigmaG(:),tauG(:)
    real(8)::sxx0,sxy0,syy0,psi,theta
    syy0=sigma0
    sxy0=syy0*muinit
    psi=45d0
    sxx0=syy0*(1d0+2*sxy0/(syy0*dtan(2*psi/180d0*pi)))
    write(*,*) 'sxx0,sxy0,syy0'
    write(*,*) sxx0,sxy0,syy0
    do i=1,size(velG)
        tauG(i)=sxy0*cos(2*ang(i))+0.5d0*(sxx0-syy0)*sin(2*ang(i))
        sigmaG(i)=sin(ang(i))**2*sxx0+cos(ang(i))**2*syy0+sxy0*sin(2*ang(i))
        phiG(i)=ag(i)*dlog(2*vref/velinit*sinh(tauG(i)/sigmaG(i)/ag(i)))
        omega=exp((PhiG(i)-mu0)/bG(i))*velinit/vref
        !write(*,*) Omega
    end do

  end subroutine initcond

  subroutine coordinate3d(imax,jmax,ds,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
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

  end subroutine coordinate3d

  subroutine coordinate2dn(geofile,NCELLg,xel,xer,yel,yer,xcol,ycol,ang)
    implicit none
    integer,intent(in)::NCELLg
    character(128),intent(in)::geofile
    real(8),intent(out)::xel(:),xer(:),yel(:),yer(:),xcol(:),ycol(:),ang(:)
    integer::i,j,k,file_size,n,Np,Nm,ncellf
    real(8)::dx,xr(0:NCELLg),yr(0:NCELLg),nx(NCELLg),ny(NCELLg),r(NCELLg)
    real(8),allocatable::data(:)

    open(20,file=geofile,access='stream')
    read(20) xel,xer,yel,yer
    !write(*, '(a, i12)') "file size = ", file_size

    do i=1,NCELLg
      ang(i)=datan((yer(i)-yel(i))/(xer(i)-xel(i)))
      xcol(i)=0.5d0*(xel(i)+xer(i))
      ycol(i)=0.5d0*(yel(i)+yer(i))
    end do

    ! open(14,file='top3.dat')
    ! do i=1,NCELLg
    !   write(14,'(7e16.6)') xcol(i),ycol(i),ang(i),xel(i),xer(i),yel(i),yer(i)
    ! end do
    ! close(14)

    return
  end subroutine


  subroutine coordinate2dp(NCELLg,ds,xel,xer,xcol)
    implicit none
    integer,intent(in)::NCELLg
    real(8),intent(in)::ds
    real(8),intent(out)::xel(:),xer(:),xcol(:)
    integer::i,j,k

    !flat fault with element size ds
    !open(14,file='top.dat')
    do i=1,NCELLg
      xel(i)=(i-1)*ds
      xer(i)=i*ds
      xcol(i)=0.5d0*(xel(i)+xer(i))
      !write(14,'(3e16.6)') xcol(i),xel(i),xer(i)
    enddo
    !close(14)
    return
  end subroutine

  subroutine params(problem,NCELLg,a0,b0,dc0,vc0,ag,bg,dcg,vcg,fwg,vwg)
    character(128),intent(in)::problem
    integer,intent(in)::NCELLg
    real(8),intent(in)::a0,b0,dc0,vc0
    real(8),intent(out)::ag(:),bg(:),dcg(:),vcg(:),fwg(:),vwg(:)
    integer::i
    do i=1,NCELLg
      ag(i)=a0
      !if(abs(i-NCELLg/2).gt.NCELLg/4) ag(i)=0.024d0
      bg(i)=b0
      dcg(i)=dc0
      vcg(i)=vc0
      fwg(i)=fw0
      vwg(i)=vw0
      !write(*,*)ag(i),bg(i),dcg(i)
    end do

  end subroutine
  subroutine loading(problem,NCELLg,sr,taudotG,sigdotG)
    character(128),intent(in)::problem
    integer,intent(in)::NCELLg
    real(8),intent(in)::sr
    real(8),intent(out)::taudotg(:),sigdotg(:)
    real(8)::factor,edge,ret1,ret2
    integer::i
    character(128)::v
    select case(problem)
    case('2dp','3dp')
      taudotG=sr
      sigdotG=0d0
    case('2dn')
      open(15,file='sr')
      do i=1,NCELLg
      select case(load)
      case(0)
        taudotG(i)=sr*cos(2*ang(i))
        sigdotG(i)=sr*sin(2*ang(i))
      case(1)
        edge=ds*NCELLg/2

          v='s'
          call kern(v,xcol(i),ycol(i),-99*edge,ycol(1),-edge,ycol(1),ang(i),0d0,ret1)
          call kern(v,xcol(i),ycol(i),edge,ycol(NCELLg),99*edge,ycol(NCELLg),ang(i),0d0,ret2)
          taudotG(i)=vpl*(ret1+ret2)
          v='n'
          call kern(v,xcol(i),ycol(i),-99*edge,ycol(1),-edge,ycol(1),ang(i),0d0,ret1)
          call kern(v,xcol(i),ycol(i),edge,ycol(NCELLg),99*edge,ycol(NCELLg),ang(i),0d0,ret2)
          sigdotG(i)=vpl*(ret1+ret2)
      end select
      write(15,*) taudotG(i),sigdotG(i)
    end do
    close(15)
    !case('2dpv','2dnv')
    !  factor=rigid/(2.d0*pi*(1.d0-pois))
    !  edge=-ds*NCELLg
    !  do i=1,NCELLg
    !    taudotg(i)=vpl*factor*(1.d0/xcol(i)-1.d0/(xcol(i)-edge))
    !    sigdotG(i)=0d0
    !  end do
    end select
  end subroutine

  subroutine derivs(x, y, dydx)!,,st_leafmtxp,st_bemv,st_ctl)
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    !integer,intent(in) :: NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in) :: x
    real(8),intent(in) ::y(:)
    real(8),intent(out) :: dydx(:)
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),efftmp(NCELL),phitmp(NCELL)
    real(8) :: dlnvdt(NCELL),dtaudt(NCELL),dsigdt(NCELL),deffdt(NCELL),dphidt(NCELL)
    real(8) :: sum_gs(NCELL),sum_gn(NCELL)
    real(8)::veltmpg(NCELLg),sum_gsg(NCELLg),sum_gng(NCELLg),efftmpG(NCELLg)
    real(8):: c1, c2, c3, arg,c,g,tauss
    integer :: i, j, nc,ierr,lrtrn

    !if(my_rank.eq.0) then
    select case(problem)
    case('2dp','3dp')
      do i = 1, NCELL
        phitmp(i) = y(2*i-1)
        tautmp(i) = y(2*i)
        sigmatmp(i)=sigma(i)
        veltmp(i) = 2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      select case(load)
      case(0)
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
        call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        do i=1,NCELL
          sum_gn(i)=0.d0
          sum_gs(i)=sum_gs(i)+taudot(i)
        end do
      case(1)
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG-vpl)
        call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      end select
      !call MPI_SCATTERv(sum_gsG,NCELL,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)



    case('2dn')
      do i = 1, NCELL
        phitmp(i) = y(3*i-2)
        tautmp(i) = y(3*i-1)
        sigmatmp(i) = y(3*i)
        veltmp(i) = 2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))
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

      do i=1,NCELL
        sum_gs(i)=sum_gs(i)+taudot(i)
        sum_gn(i)=sum_gn(i)+sigdot(i)
      end do

    end select

    ! select case(problem)
    ! case('2dp','3dp','2dpv')
    !   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !   call MPI_ALLGATHERv(veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
    !
    !   !matrix-vector mutiplation
    !   select case(load)
    !   case(0)
    !     lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
    !   case(1)
    !     lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG-vpl)
    !   end select
    !   !call MPI_SCATTERv(sum_gsG,NCELL,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !   call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !   do i=1,NCELL
    !     sum_gn(i)=0.d0
    !     sum_gs(i)=sum_gs(i)+taudot(i)
    !   end do
    !
    ! case('2dn','2dnv')
    !   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !   call MPI_ALLGATHERv(Veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,                &
    !   &     MPI_REAL8,MPI_COMM_WORLD,ierr)
    !
    !   !matrix-vector mutiplation
    !   st_bemv%v='s'
    !   !veltmp=1.d-6
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
    !   !if(my_rank.eq.0) write(*,*) sum_gs
    !   !stop
    !   st_bemv%v='n'
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxpn,st_bemv,st_ctl,sum_gnG,veltmpG)
    !   !if(my_rank.eq.10) write(*,*) sum_gn(1)
    !
    !   call MPI_SCATTERv(sum_gsG,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !   call MPI_SCATTERv(sum_gnG,rcounts,displs,MPI_REAL8,sum_gn,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !
    !   do i=1,NCELL
    !     sum_gs(i)=sum_gs(i)+taudot(i)
    !     sum_gn(i)=sum_gn(i)+sigdot(i)
    !   end do
    !
    ! end select
    !end if


    select case(law)
      ! aging law (Linker & Dieterich, 1992)
    case('a')
      call deriv_a(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
      ! slip law
    case('s')
      call deriv_s(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
      ! RFL in FDMAP (Dunham+ 2011)
    case('d')
      call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)

    end select

    select case(problem)
    case('2dp','3dp')
      do i = 1, NCELL
        dydx(2*i-1) = dphidt( i )
        dydx(2*i) = dtaudt( i )
      enddo
    case('2dn')
      do i = 1, NCELL
        dydx(3*i-2) = dphidt( i )
        dydx(3*i-1) = dtaudt( i )
        dydx(3*i) = dsigdt( i )
      enddo
    end select

    return
  end subroutine

  subroutine deriv_a(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
    implicit none
    integer::i
    real(8)::arg
    !type(t_deriv),intent(in) ::
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

  subroutine deriv_va(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,efftmp,dlnvdt,dtaudt,dsigdt,deffdt)
    implicit none
    integer::i
    real(8)::arg
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs(:),sum_gn(:),veltmp(:),tautmp(:),sigmatmp(:),efftmp(:)
    !real(8),intent(in)::a(:),b(:),dc(:),vc(:)
    !real(8),intent(in)::mu0,vref,vs,rigid,alpha
    real(8),intent(out)::dlnvdt(:),dtaudt(:),dsigdt(:),deffdt(:)
    do i=1,size(sum_gs)
      dsigdt(i)=sum_gn(i)
      arg=dc(i)/vref*(exp((tautmp(i)/sigmatmp(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
      dlnvdt(i)=sum_gs(i)-b(i)*sigmatmp(i)/(arg+dc(i)/vc(i))*(1.d0-veltmp(i)*arg/dc(i))+(tautmp(i)/sigmatmp(i)-alpha)*dsigdt(i)
      dlnvdt(i)=dlnvdt(i)/(a(i)*sigmatmp(i)+0.5d0*rigid*veltmp(i)/vs)
      dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
      deffdt(i)=veltmp(i)-efftmp(i)/tr
    end do
  end subroutine

  subroutine deriv_s(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
    implicit none
    integer::i
    real(8)::arg
    !type(t_deriv),intent(in) ::
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

  subroutine deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)
    implicit none
    integer::i
    real(8)::fss,dvdtau,dvdsig,dvdphi
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs(:),sum_gn(:),phitmp(:),tautmp(:),sigmatmp(:),veltmp(:)
    real(8),intent(out)::dphidt(:),dtaudt(:),dsigdt(:)
    do i=1,size(sum_gs)
      !write(*,*) 'vel',veltmp(i)
      dsigdt(i)=sum_gn(i)
      !write(*,*) 'dsigdt',dsigdt(i)

      fss=mu0+(a(i)-b(i))*dlog(abs(veltmp(i))/vref)
      !fss=fw(i)+(fss-fw(i))/(1.d0+(veltmp(i)/vw(i))**8)**0.125d0 !flash heating
      !slip law
      !dphidt(i)=-abs(veltmp(i))/dc(i)*(abs(tautmp(i))/sigmatmp(i)-fss)
      !aging law
      dphidt(i)=b(i)*vref/dc(i)*exp((mu0-phitmp(i))/b(i))-abs(veltmp(i))/dc(i)

      dvdtau=2*vref*exp(-phitmp(i)/a(i))*cosh(tautmp(i)/sigmatmp(i)/a(i))/(a(i)*sigmatmp(i))
      dvdsig=-2*vref*exp(-phitmp(i)/a(i))*cosh(tautmp(i)/sigmatmp(i)/a(i))*tautmp(i)/(a(i)*sigmatmp(i)**2)
      !dvdphi=2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
      dvdphi=-veltmp(i)/a(i)
      !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*phitmp(i)*dvdsig*sigmatmp(i))
      dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))
      dtaudt(i)=dtaudt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      !write(*,*) rigid/vs*dvdtau
      if(veltmp(i).le.0d0) then
        dvdtau=2*vref*exp(-phitmp(i)/a(i))*cosh(tautmp(i)/sigmatmp(i)/a(i))/(a(i)*sigmatmp(i))
        dvdsig=-2*vref*exp(-phitmp(i)/a(i))*cosh(tautmp(i)/sigmatmp(i)/a(i))*tautmp(i)/(a(i)*sigmatmp(i)**2)
        !sign ok?
        dvdphi=2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
        !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(-dvdphi*phitmp(i)*dvdsig*sigmatmp(i))
        dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))
        dtaudt(i)=dtaudt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      end if
    end do

  end subroutine

  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
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
    integer :: i,ierr
    real(8) :: errmax,h,xnew,htemp,errmax_gb
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4,hmax=1d5

    h=htry
    do while(.true.)
      call rkck(y,dydx,x,h,ytemp,yerr)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
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

    hnext=min(2*h,SAFETY*h*(errmax_gb**PGROW),1d8)

    hdid=h
    x=x+h
    y(:)=ytemp(:)
    return
  end subroutine

  !---------------------------------------------------------------------
  subroutine rkck(y,dydx,x,h,yout,yerr)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
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
    call derivs(x+a2*h, ytemp, ak2)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
    end do

    !     -- 3rd step --
    call derivs(x+a3*h, ytemp, ak3)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
    end do

    !     -- 4th step --
    call derivs(x+a4*h, ytemp, ak4)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
    end do

    !     -- 5th step --
    call derivs(x+a5*h, ytemp, ak5)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do

    !     -- 6th step --
    call derivs(x+a6*h, ytemp, ak6)!,,st_leafmtxp,st_bemv,st_ctl)
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
end program
