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

 !job ID
  integer::number
  !# of elements and timestep
  integer::NCELL,NCELLg,NSTEP
  integer::imax,jmax !for 3dp,3dhr

  !for HACApK
  real(8),allocatable ::coord(:,:),vmax(:)
  real(8)::eps_h
  type(st_HACApK_lcontrol) :: st_ctl,st_ctl2
  type(st_HACApK_leafmtxp) :: st_leafmtxps,st_leafmtxpn
  type(st_HACApK_leafmtxp) :: st_leafmtxp_s,st_leafmtxp_n,st_leafmtxp_d,st_leafmtxp_c
  type(st_HACApK_leafmtxp) :: st_leafmtxp_s2,st_leafmtxp_n2,st_leafmtxp_d2
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xx,st_leafmtxp_xy,st_leafmtxp_yy
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xz,st_leafmtxp_yz,st_leafmtxp_zz
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xx2,st_leafmtxp_xy2,st_leafmtxp_yy2
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xz2,st_leafmtxp_yz2,st_leafmtxp_zz2
  type(st_HACApK_calc_entry) :: st_bemv

  !for Lattice H matrix
  real(8),allocatable::wws(:)
  type(st_HACApK_latticevec) :: st_vel,st_sum
  type(st_HACApK_LHp) :: st_LHp,st_LHp_s,st_LHp_d,st_LHp_n,st_LHp_xx,st_LHp_xy,st_LHp_yy
  type(st_HACApK_LHp) :: st_LHp_s2,st_LHp_d2,st_LHp_n2

  !for MPI communication and time
  integer::counts2,icomm,np,npd,ierr,my_rank,npgl
  integer,allocatable::displs(:),rcounts(:),vars(:)
  integer:: date_time(8)
  character(len=10):: sys_time(3)
  real(8)::time1,time2,time3,time4,timer,timeH

  !for fault geometry
  real(8),allocatable::xcol(:),ycol(:),zcol(:),ds(:)
  real(8),allocatable::xs1(:),xs2(:),xs3(:),xs4(:) !for 3dp
  real(8),allocatable::zs1(:),zs2(:),zs3(:),zs4(:) !for 3dp
  real(8),allocatable::ys1(:),ys2(:),ys3(:),ys4(:) !for 3dn
  real(8),allocatable::xel(:),xer(:),yel(:),yer(:),ang(:),angd(:)
  real(8),allocatable::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:)

  !parameters for each elements
  real(8),allocatable::a(:),b(:),dc(:),f0(:),fw(:),vw(:),vc(:),taudot(:),tauddot(:),sigdot(:),kLv(:),kTv(:)
  real(8),allocatable::ag(:),bg(:),dcg(:),f0g(:),taug(:),sigmag(:),velg(:),taudotg(:),sigdotg(:),pfg(:),qin(:)

  !variables
  real(8),allocatable::psi(:),vel(:),tau(:),sigma(:),disp(:),mu(:),rupt(:),idisp(:),velp(:),pfhyd(:),cslip(:)
  real(8),allocatable::taus(:),taud(:),vels(:),veld(:),disps(:),dispd(:),rake(:),pf(:),sigmae(:),ks(:),qflow(:),kp(:),phi(:)

  real(8),allocatable::rdata(:),qvals(:,:),qtimes(:,:)
  integer,allocatable::iwell(:),jwell(:),kleng(:)
  integer::lp,i,i_,j,k,kstart,kend,m,counts,interval,lrtrn,nl,ios,nmain,rk,nout(10),file_size,nwell
  integer::hypoloc(1),load,eventcount,thec,inloc,sw,errmaxloc,niter,mvelloc(1),locid(10)

  !controls
  logical::aftershock,buffer,nuclei,slipping,outfield,slipevery,limitsigma,dcscale,slowslip,slipfinal,outpertime
  logical::initcondfromfile,parameterfromfile,injectionfromfile,backslip,sigmaconst,foward,inverse,geofromfile,restart,latticeh,pressuredependent,pfconst
  character*128::fname,dum,law,input_file,problem,geofile,param,pvalue,slipmode,project,parameter_file,outdir,command,bcl,bcr,bc,evlaw,setting
  character(128)::injection_file
  real(8)::a0,a1,b0,dc0,sr,omega,theta,dtau,tiny,moment,wid,normal,ieta,meanmu,meanmuG,meandisp,meandispG,moment0,mvel,mvelG
  real(8)::vc0,mu0,onset_time,tr,vw0,fw0,velmin,tauinit,intau,trelax,maxnorm,maxnormG,minnorm,minnormG,sigmainit,muinit
  real(8)::r,vpl,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxeventi,sparam,tmax,dtmax,tout,dummy(10)
  real(8)::ds0,amp,mui,velinit,phinit,velmax,maxsig,minsig,v1,dipangle,crake,s,sg,q0
  real(8)::kpmax,kpmin,kp0,kT,kL,s0,ksinit,dtout,pfinit,pbcl,pbcr,lf,eta,beta,phi0,str,cc,td,cd

  !random_number
  integer,allocatable::seed(:)
  integer::seedsize

  !for time integration
  real(8)::x !time
  real(8),allocatable::y(:),yscal(:),dydx(:),yg(:)
  real(8)::eps_r,errmax_gb,dtinit,dtnxt,dttry,dtdid,dtmin,tp,fwid

  integer::r1,r2,r3,NVER,amari,out,kmax,loci,locj,loc,stat,nth
  integer,allocatable::rupsG(:)

  !initialize
  icomm=MPI_COMM_WORLD
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )

  if(my_rank.eq.0) then
    write(*,*) '# of MPI cores', np
  end if
  !input file must be specified when running
  !example) mpirun -np 16 ./ha.out default.in
  call get_command_argument(1,input_file,status=stat)

  open(33,file=input_file,iostat=ios)
  !if(my_rank.eq.0) write(*,*) 'input_file',input_file
  if(ios /= 0) then
    write(*,*) 'Failed to open inputfile'
    stop
  end if

  !get filenumber
  number=0
  if(input_file(1:11)=='examples/in') then
    input_file=adjustl(input_file(12:))
    !write(*,*) input_file
    read(input_file,*) number
    !write(*,*) number
  end if

  time1=MPI_Wtime()

  !default parameters
  nmain=1000000
  eps_r=1d-4
  eps_h=1d-4
  velmax=1d7
  velmin=1d-16
  tmax=1d12
  td=tmax
  dipangle=0d0
  sigmaconst=.false.
  foward=.false.
  inverse=.false.
  slipfinal=.false.
  restart=.false.
  maxsig=300d0
  minsig=0d0
  dtinit=1d0
  tp=86400d0
  initcondfromfile=.false.
  parameterfromfile=.false.
  injectionfromfile=.false.
  pressuredependent=.false.
  pfconst=.false.
  bc='Neumann'
  evlaw='aging'
  !number=0


  do while(ios==0)
    read(33,*,iostat=ios) param,pvalue
    !write(*,*) param,pvalue
    select case(param)
    case('problem')
      read (pvalue,*) problem
    case('ncellg')
      read (pvalue,*) ncellg
    case('imax')
      read (pvalue,*) imax
    case('jmax')
      read (pvalue,*) jmax
    case('nstep')
      read (pvalue,*) nstep
    case('loc')
      read (pvalue,*) loc
    case('filenumber')
      read (pvalue,*) number
    case('ds')
      read (pvalue,*) ds0
    case('velmax')
      read (pvalue,*) velmax
    case('velmin')
      read (pvalue,*) velmin
    case('a')
      read (pvalue,*) a0
    case('a1')
      read (pvalue,*) a1
    case('b')
      read (pvalue,*) b0
    case('dc')
      read (pvalue,*) dc0
    case('f0')
      read (pvalue,*) mu0
    case('load')
      read (pvalue,*) load
    case('sr')
      read (pvalue,*) sr
    case('vpl')
      read (pvalue,*) vpl
    case('interval')
      read (pvalue,*) interval
    case('geometry_file')
      read (pvalue,'(a)') geofile
    case('velinit')
      read (pvalue,*) velinit
    case('tauinit')
      read (pvalue,*) tauinit
    case('sigmainit')
      read (pvalue,*) sigmainit
    case('pfinit')
      read (pvalue,*) pfinit
    case('muinit')
      read (pvalue,*) muinit
    case('ksinit')
      read (pvalue,*) ksinit
    case('dtinit')
      read (pvalue,*) dtinit
    case('sparam')
      read (pvalue,*) sparam
    case('q0')
      read (pvalue,*) q0
    case('eta')
      read (pvalue,*) eta
    case('phi0')
      read (pvalue,*) phi0
    case('beta')
      read (pvalue,*) beta
    case('kpmax')
      read (pvalue,*) kpmax
    case('kpmin')
      read (pvalue,*) kpmin
    case('kp0')
      read (pvalue,*) kp0
    case('s0')
      read (pvalue,*) s0
    case('kL')
      read (pvalue,*) kL
    case('kT')
      read (pvalue,*) kT
    case('cd')
      read (pvalue,*) cd
    case('td')
      read (pvalue,*) td
    case('pbcl')
      read (pvalue,*) pbcl
    case('pbcr')
      read (pvalue,*) pbcr
    case('tmax')
      read (pvalue,*) tmax
    case('eps_r')
      read (pvalue,*) eps_r
    case('eps_h')
      read (pvalue,*) eps_h
    case('slipevery')
      read (pvalue,*) slipevery
    case('backslip')
      read (pvalue,*) backslip
    case('pressuredependent')
      read (pvalue,*) pressuredependent
    case('pfconst')
      read (pvalue,*) pfconst
    case('limitsigma')
      read (pvalue,*) limitsigma
    case('sigmaconst')
      read(pvalue,*) sigmaconst
    case('foward')
      read(pvalue,*) foward
    case('inverse')
      read(pvalue,*) inverse
    case('geofromfile')
      read(pvalue,*) geofromfile
    case('maxsig')
      read(pvalue,*) maxsig
    case('minsig')
      read(pvalue,*) minsig
    case('crake')
      read(pvalue,*) crake
    case('dipangle')
      read(pvalue,*) dipangle
    case('outpertime')
      read(pvalue,*) outpertime
    case('restart')
      read(pvalue,*) restart
    case('latticeh')
      read(pvalue,*) latticeh
    case('dtout')
      read(pvalue,*) dtout
    case('parameterfromfile')
      read(pvalue,*) parameterfromfile
    case('parameter_file')
      read(pvalue,'(a)') parameter_file
    case('bc')
      read (pvalue,'(a)') bc
    case('evlaw')
      read (pvalue,'(a)') evlaw
    case('setting')
      read (pvalue,'(a)') setting
    case('injectionfromfile')
      read(pvalue,*) injectionfromfile
    case('injection_file')
      read(pvalue,'(a)') injection_file
    end select
  end do
  close(33)
  tmax=tmax*365*24*3600

  !limitsigma=.true.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !MPI setting
  !NCELLg=2*NL*NL
  !if(problem.eq.'3dp') then
  nmain=ncellg
  select case(problem)
  case('3dp','3dph','3dhfr')
    NCELLg=imax*jmax
    loc=loci*(imax-1)+locj
  ! case('3dnf','3dn')
  !   NCELLg=imax*jmax*2
  end select
  !end if

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
  allocate(vars(NCELL))
  do i=1,NCELL
    vars(i)=displs(my_rank+1)+i
    !write(*,*) displs(my_rank+1),i,vars(i)
  end do

  !stop
  !call varscalc(NCELL,displs,vars)
  if(my_rank.eq.0) then
    write(*,*) 'job number',number
  end if

  !allocation
  !allocate(a(NCELLg),b(NCELLg),dc(NCELLg),f0(NCELLg),fw(NCELLg),vw(NCELLg),vc(NCELLg),taudot(NCELLg),tauddot(NCELLg),sigdot(NCELLg))
  allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg),ds(NCELLg))
  xcol=0d0;ycol=0d0;zcol=0d0

  select case(problem)
  ! case('2dp','2dh')
  !   allocate(xel(NCELLg),xer(NCELLg))
  !   xel=0d0;xer=0d0
  !   allocate(phi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg),velp(NCELLg))
  ! case('2dn','2dnh','2dn3','25d')
  !   allocate(ang(NCELLg),xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg))
  !   ang=0d0;xel=0d0;xer=0d0;yel=0d0;yer=0d0
  !   allocate(phi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg),velp(NCELLg))
  !   allocate(taui(NCELLg),sigmai(NCELLg))
  case('3dp','3dph')
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
    xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
    zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
  allocate(psi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg),cslip(NCELLg),pf(NCELLg),sigmae(NCELLg),pfhyd(NCELLg))
  psi=0d0;vel=0d0;tau=0d0;sigma=0d0;disp=0d0
  allocate(a(NCELLg),b(NCELLg),dc(NCELLg),f0(NCELLg),taudot(NCELLg),sigdot(NCELLg),ks(NCELLg),kp(NCELLg),qflow(NCELLg),kLv(NCELLg),kTv(NCELLg),phi(NCELLg))
  taudot=0d0;sigdot=0d0;ks=0d0

  ! case('3dn','3dh','3dnf','3dhf')
  !   allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg))
  !   allocate(ys1(NCELLg),ys2(NCELLg),ys3(NCELLg))
  !   allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg))
  !   allocate(ev11(NCELLg),ev12(NCELLg),ev13(NCELLg))
  !   allocate(ev21(NCELLg),ev22(NCELLg),ev23(NCELLg))
  !   allocate(ev31(NCELLg),ev32(NCELLg),ev33(NCELLg))
  !   xs1=0d0; xs2=0d0; xs3=0d0
  !   ys1=0d0; ys2=0d0; ys3=0d0
  !   zs1=0d0; zs2=0d0; zs3=0d0
  !   allocate(phi(NCELLg),vels(NCELLg),veld(NCELLG),taus(NCELLg),taud(NCELLg),sigma(NCELLg),disp(NCELLg),disps(NCELLg),dispd(NCELLG),mu(NCELLg),rake(NCELLg),vel(NCELLG),tau(NCELLg),idisp(NCELLg),velp(NCELLg))
  ! case('3dhfr')
  !   allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
  !   allocate(ys1(NCELLg),ys2(NCELLg),ys3(NCELLg),ys4(NCELLg))
  !   allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
  !   allocate(ang(NCELLg),angd(NCELLg))
  !   xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
  !   ys1=0d0; ys2=0d0; ys3=0d0; ys4=0d0
  !   zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
  !   angd=0d0; ang=0d0
  !   allocate(phi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg),velp(NCELLg))
  end select

  allocate(y(4*NCELL),yscal(4*NCELL),dydx(4*NCELL),yg(4*NCELLg))

  !mesh generation
  if(my_rank==0) write(*,*) 'Generating mesh'
  select case(problem)
  case('2dp','2dpah')
    call coordinate2dp(NCELLg,ds0,xel,xer,xcol)
  case('2dph')
    call coordinate2dph()
  case('2dn')
    open(20,file=geofile,status='old',iostat=ios)
    if(ios /= 0) then
      if(my_rank==0)write(*,*) 'error: Failed to open geometry file'
      stop
    end if
    do i=1,NCELLg
      read(20,*) xel(i),xer(i),yel(i),yer(i)
    end do
    close(20)
    call coordinate2dn()
  case('3dp')
    call coordinate3dp(imax,jmax,ds0,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
    ds=ds0*ds0
  case('3dnr','3dhr')
    open(20,file=geofile,status='old',iostat=ios)
    if(ios /= 0) then
     if(my_rank==0)write(*,*) 'error: Failed to open geometry file'
     stop
    end if
    do i=1,NCELLg
     read(20,*) xcol(i),ycol(i),zcol(i),ang(i),angd(i)
     ds(i)=ds0*ds0
    end do
    close(20)
  case('3dph')
    call coordinate3ddip(imax,jmax,ds0,dipangle)
    ds=ds0*ds0
  case('3dnt','3dht')
    !.stl format
    open(12,file=geofile,iostat=ios)
    if(ios /= 0) then
      if(my_rank==0)write(*,*) 'error: Failed to open geometry file'
      stop
    end if
    do while(.true.)
      read(12,*) dum
      if(dum=='facet') exit
    end do
    !write(*,*) ios
    do k=1,ncellg
      !read(12,*)
      read(12,*) !outer loop
      read(12,*) dum,xs1(k),ys1(k),zs1(k) !vertex
      read(12,*) dum,xs2(k),ys2(k),zs2(k) !vertex
      read(12,*) dum,xs3(k),ys3(k),zs3(k) !vertex
      read(12,*) !end loop
      read(12,*) !endfacet
      read(12,*) !facet
      xcol(k)=(xs1(k)+xs2(k)+xs3(k))/3
      ycol(k)=(ys1(k)+ys2(k)+ys3(k))/3
      zcol(k)=(zs1(k)+zs2(k)+zs3(k))/3
    !  write(*,*)ios
      !if(my_rank==0)write(*,'(9e17.8)') xs1(k),ys1(k),zs1(k),xs2(k),ys2(k),zs2(k),xs3(k),ys3(k),zs3(k)
    end do
    call evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33,ds)
  end select

  !call initcond3dn(phi,sigma,taus,taud)
  !stop
  !random number seed
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do i = 1, seedsize
    call system_clock(count=seed(i))
  end do
  call random_seed(put=seed(:))

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !stop

  !HACApK setting
  lrtrn=HACApK_init(NCELLg,st_ctl,st_bemv,icomm)
  allocate(coord(NCELLg,3))

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

  !generate kernel (H-matrix aprrox)
  if(my_rank.eq.0) write(*,*) 'Generating kernel'
  do i=1,NCELLg
    coord(i,1)=xcol(i)
    coord(i,2)=ycol(i)
    coord(i,3)=zcol(i)
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  st_ctl%param(8)=1

  select case(problem)
  case('2dp','2dpah','2dn3','3dp')
    sigmaconst=.true.
  end select
  st_bemv%md='s'
  st_bemv%v='s'
  lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,eps_h)
  if(.not.sigmaconst) then
    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxpn,st_bemv,st_ctl,coord,eps_h)
  end if

  !end if
  !write(*,*) my_rank,st_ctl%lpmd(33),st_ctl%lpmd(37)

  !setting frictional parameters
  a=a0
  b=b0
  dc=dc0
  f0=mu0

  if(parameterfromfile) then
    open(99,file=parameter_file)
    write(*,*) "reading friction parameters from file"
    do i=1,ncellg
      read(99,*) a(i),b(i),dc(i),f0(i)
    end do
    close(99)
  end if

  if(.not.backslip) then
    taudot=sr
    sigdot=0d0
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !setting initial condition
  call initcond_bgflow()
  select case(setting)
  case('injection')
    call initcond_injection()
  end select

  if(pressuredependent) then
    do i=1,ncellg
      kp(i)=ks(i)*exp(-sigmae(i)/s0)
    end do
  else 
    kp=ks
  end if
  qflow=q0

  if(my_rank==0) then
    ! write(fname,'("output/ind",i0,"_",i0,".dat")') number,my_rank
    ! nout=my_rank+100
    ! open(nout(1),file=fname)
    ! !write(nout)st_sum%lodc(1:NCELL)
    ! do i=1, ncell
    !   write(nout(1),'(i0)') st_sum%lodc(i)
    ! end do
    ! close(nout(1))

    ! write(fname,'("output/xyz",i0,"_",i0,".dat")') number,my_rank
    ! nout(1)=my_rank+100
    ! open(nout(1),file=fname)
    ! do i=1,ncell
    !   i_=st_sum%lodc(i)
    !   write(nout(1),'(3e15.6)') xcol(i_),ycol(i_),zcol(i_)
    ! end do
    ! close(nout(1))
    nout(1)=100

    write(fname,'("output/vel",i0,".dat")') number
    open(nout(1),file=fname,form='unformatted',access='stream',status='replace')
    nout(2)=nout(1)+np
    write(fname,'("output/slip",i0,".dat")') number
    open(nout(2),file=fname,form='unformatted',access='stream',status='replace')
    nout(3)=nout(2)+np
    write(fname,'("output/sigma",i0,".dat")') number
    open(nout(3),file=fname,form='unformatted',access='stream',status='replace')
    nout(4)=nout(3)+np
    write(fname,'("output/tau",i0,".dat")') number
    open(nout(4),file=fname,form='unformatted',access='stream',status='replace')
    nout(5)=nout(4)+np
    write(fname,'("output/ks",i0,".dat")') number
    open(nout(5),file=fname,form='unformatted',access='stream',status='replace')
    nout(6)=nout(5)+np
    write(fname,'("output/pf",i0,".dat")') number
    open(nout(6),file=fname,form='unformatted',access='stream',status='replace')
    nout(7)=nout(6)+np
    write(fname,'("output/qflow",i0,".dat")') number
    open(nout(7),file=fname,form='unformatted',access='stream',status='replace')
    nout(8)=nout(7)+np
    write(fname,'("output/EQslip",i0,".dat")') number
    open(nout(8),file=fname,form='unformatted',access='stream',status='replace')

    write(fname,'("output/time",i0,".dat")') number
    open(50,file=fname)
    write(fname,'("output/event",i0,".dat")') number
    open(51,file=fname)
    write(fname,'("output/monitor",i0,".dat")') number
    open(52,file=fname)
    nl=5

    call find_locid(imax/2,jmax/2,locid(1))
    call find_locid(imax/4,jmax/2,locid(2))
    call find_locid(imax*3/4,jmax/2,locid(3))
    call find_locid(imax/2,jmax/4,locid(4))
    call find_locid(imax/2,jmax*3/4,locid(5))

    do i=1,nl
      write(fname,'("output/local",i0,"-",i0,".dat")') number,i
      open(52+i,file=fname)
    end do

    open(19,file='job.log',position='append')
    call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    write(19,'(a20,i0,a6,a12,a6,a12,a4,i0)') 'Starting job number=',number,'date',sys_time(1),'time',sys_time(2),'np',np
    close(19)
    !allocate(locid(9))
    !ocid=(/1850,2000,2050,2100,2150,2250,2350,2500,2750/)
    !call open_BP6()

  end if

  !setting minimum time step by CFL condition
  !dtmin=0.5d0*ds/(vs*sqrt(3.d0))
  
  k=0
  errmax_gb=0d0
  dtdid=0d0
  mvelG=maxval(vel)
  meanmuG=sum(mu)/NCELLg
  minnormG=minval(sigmae)
  maxnormG=maxval(sigmae)
  meandispG=sum(disp)/NCELLg

  if(injectionfromfile) call input_well()
  
  if(my_rank==0) then
    write(50,'(i7,f19.4)') k,x
    call output_monitor()
    write(nout(1)) vel
    write(nout(2)) disp
    write(nout(3)) sigmae
    write(nout(4)) tau
    write(nout(5)) kp
    write(nout(6)) pf
    write(nout(7)) qflow
    !call output_BP6()
  end if
  dtnxt = dtinit

tout=dtout*365*24*60*60

  x=0.d0 !x is time
  timer=0d0
  k=0
  rk=0

  dtnxt = dtinit
  !outv=1d-6
  slipping=.false.
  eventcount=0
  sg=ds0*ds0*NCELLg
   

    !$omp parallel do
    do i=1,NCELL
      i_=vars(i)
      y(4*i-3) = psi(i_)
      y(4*i-2) = tau(i_)
      y(4*i-1)=sigmae(i_)
      y(4*i) = ks(i_)
    end do
      !$omp end parallel do
    !call MPI_SCATTERv(yG,3*rcounts,3*displs,MPI_REAL8,y,3*NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  
  !stop
  time2=MPI_Wtime()
  if(my_rank.eq.0) write(*,*) 'Finished all initial processing, time(s)=',time2-time1
  call MPI_BARRIER(MPI_COMM_WORLD,ierr);time1=MPI_Wtime()
  do k=1,NSTEP
       !parallel computing for Runge-Kutta
    dttry = dtnxt
    !time3=MPI_Wtime()

    !update state and stress with explicit solver
    !write(*,*) vel(1489),y(1489),sigmae(1489)
    call rkqs(y,dydx,x,dttry,eps_r,dtdid,dtnxt,errmax_gb,errmaxloc)
    call MPI_ALLGATHERv(y,4*NCELL,MPI_REAL8,yG,4*rcounts,4*displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

    do i=1,NCELLg
      psi(i) = yg(4*i-3)
      tau(i) = yg(4*i-2)
      sigmae(i)=yg(4*i-1)
      !sigma(i)=sigmainit
      !sigmae(i)=sigma(i)-pf(i)

      mu(i)=tau(i)/sigmae(i)
      vel(i)= 2*vref*exp(-psi(i)/a(i))*sinh(tau(i)/sigmae(i)/a(i))
      !disp(i)=disp(i)+vel(i)*dtdid*0.5d0
      !sigma(i)=yg(4*i-1)+pf(i)
      ks(i)=yg(4*i)
    end do
    !write(*,*) minval(tau)

    !update pf with implicit solver
    !write(*,*)maxval(pf)
    !write(*,*)'call implicit solver'
    ! if(.not.pfconst) then
    !   if(pressuredependent) then
    !     !call implicitsolver(pf,sigma,ks,dtdid,x,dtnxt)
    !   else
        call implicitsolver2(pf,sigma,ks,dtdid,x,dtnxt,niter)
      ! end if
    ! end if
    !write(*,*)maxval(pf)

    !time4=MPI_Wtime()
    !timer=timer+time4-time3

    !limitsigma
    ! if(limitsigma) then
    !   do i=1,NCELLg
    !     if(yg(4*i-1)<minsig) sigma(i)=pf(i)+minsig
    !     if(yg(4*i-1)>maxsig) yg(4*i-1)=maxsig
    !   end do
    ! end if

    !compute physical values for control and output
    !write(*,*) vel(1489),pf(1489),sigmae(1489)
    !$omp parallel do
    do i = 1, NCELLg
      ! psi(i) = yg(4*i-3)
      ! tau(i) = yg(4*i-2)
      sigmae(i)=sigma(i)-pf(i)
      yg(4*i-1)=sigmae(i)
      ! ks(i)=yg(4*i)
      disp(i)=disp(i)+vel(i)*dtdid*0.5d0 !2nd order
      vel(i)= 2*vref*exp(-psi(i)/a(i))*sinh(tau(i)/sigmae(i)/a(i))
      disp(i)=disp(i)+vel(i)*dtdid*0.5d0
      mu(i)=tau(i)/sigmae(i)
      kp(i)=ks(i)
    end do
    !$omp end parallel do
    !write(*,*) vel(1489),pf(1489),sigmae(1489)
    call MPI_SCATTERv(yG,4*rcounts,4*displs,MPI_REAL8,y,4*NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)


    mvelG=maxval(vel)
    meanmuG=sum(mu)/NCELLg
    minnormG=minval(sigmae)
    maxnormG=maxval(sigmae)
    meandispG=sum(disp)/NCELLg

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    !stop controls
    if(mvelG>velmax) then
      if(my_rank == 0) write(*,*) 'slip rate above vmax'
      exit
    end if
    if(mvelG<velmin) then
      if(my_rank == 0) write(*,*) 'slip rate below vmin'
      exit
    end if
    if(x>tmax) then
      if(my_rank == 0) write(*,*) 'time exceeds tmax'
      exit
    end if
    if(minnormG<0.2) then
      if(my_rank == 0) write(*,*) 'normal stress below 0.2MPa'
      exit
    end if

    
    !output distribution control
    outfield=.false.
    if(mod(k,interval)==0) then
      outfield=.true.
    end if
    if(x>tout) then
      outfield=.true.
      tout=x+dtout*365*24*60*60
    end if

    if(outfield) then
      if(my_rank==0) then
        write(*,'(a,i0,f17.8,a)') 'time step=' ,k,x/365/24/60/60, ' yr'
        write(50,'(i7,f19.4)') k,x
        !if(slipping) then
        !  write(53,*) k,x/365/24/60/60,1
        !else
        !  write(53,*) k,x/365/24/60/60,0
        !end if

        close(52)
        write(fname,'("output/monitor",i0,".dat")') number
        open(52,file=fname,position='append')
        do i=1,nl
          close(52+i)
          write(fname,'("output/local",i0,"-",i0,".dat")') number,i
          open(52+i,file=fname,position='append')
        end do
      end if
      !lattice H
      if(my_rank==0) then
        write(nout(1)) vel
        write(nout(2)) disp
        write(nout(3)) sigmae
        write(nout(4)) tau
        write(nout(5)) kp
        write(nout(6)) pf
        write(nout(7)) qflow
      end if
    end if

    if(my_rank==0) then
      call output_monitor()
      do i=1,nl
        call output_local(52+i,locid(i))
      end do
    end if
    !call output_BP6()
    time4=MPI_Wtime()
    timer=timer+time4-time3

    !event list
    if(.not.slipping) then
      if(mvelG>1d-2) then
        slipping=.true.
        eventcount=eventcount+1
        moment0=meandispG
        idisp=disp
        hypoloc=maxloc(abs(vel))
        onset_time=x
        !tout=onset_time

        !onset save
        if(slipevery.and.(my_rank<npd)) then
          ! write(nout) vel
          ! write(nout2) disp
          ! write(nout3) sigmae
          ! write(nout4) tau
          ! write(nout5) ks
        end if

      end if
    end if
    !
    if(slipping) then
      if(mvelG<1d-3) then
        slipping=.false.
        !tout=x
        moment=meandispG-moment0
        !eventcount=eventcount+1
        !end of an event
        if(my_rank==0) then
          write(51,'(i0,i7,f17.2,i7,e15.6,f14.4)') eventcount,k,onset_time,hypoloc,moment,(log10(moment*rigid*sg)+5.9)/1.5
        end if
        cslip=disp-idisp
        if(my_rank==0) write(nout(8)) cslip

      end if
      !   vmaxevent=max(vmaxevent,maxval(vel))
      !   !write(53,'(i6,4e16.6)') !k,x-onset_time,sum(disp-idisp),sum(vel),sum(acg**2)
      !   !if(x-onset_time>lapse) then
      !   !  lapse=lapse+dlapse
      !   !end if
    end if


    dttry = dtnxt
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierr);time2= MPI_Wtime()
  !write(*,'(i0,5f16.4)')my_rank,time2-time1,timeo,st_ctl%time(1),st_ctl%time(2),timer

  !if(SEAS.and.my_rank.eq.0)  call output_rupt_BP()
  200  if(my_rank.eq.0) then
  write(*,*) 'time(s)', time2-time1,timer
  write(*,*) 'time for matvec(s)', sum(st_ctl%time),st_ctl%time(1)

  open(19,file='job.log',position='append')
  write(19,'(a20,i0,f16.2)') 'Finished job number=',number,time2-time1
  !open(19,file='job.log',position='append')
  close(52)
  close(50)
  close(48)
  close(47)
  close(46)
  close(44)
  close(19)
  do i=121,123
    close(i)
  end do
end if
!if(my_rank.eq.0) write(19,'(a20,i0,f16.2)')'Finished job number=',number,time2-time1
Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_s)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_n)

lrtrn=HACApK_finalize(st_ctl)
Call MPI_FINALIZE(ierr)
stop
contains
  !------------output-----------------------------------------------------------!
subroutine output_monitor()
  implicit none
  time2=MPi_Wtime()
  write(52,'(i7,f19.4,7e16.5,f16.4,i6)')k,x,log10(mvelG),meandispG,meanmuG,maxnormG,minnormG,errmax_gb,dtdid,time2-time1,niter
end subroutine
subroutine output_local(nf,loc)
  implicit none
  integer,intent(in)::nf,loc
  write(nf,'(i7,f19.4,8e16.6)')k,x,log10(vel(loc)),disp(loc),sigmae(loc),tau(loc),pf(loc),mu(loc),psi(loc),kp(loc)
end subroutine
! subroutine output_field()
!   do i=1,NCELL
!     i_=st_sum%lodc(i)
!     write(nout(1),'(i0,10e14.5,i10)') i_,xcol(i_),ycol(i_),zcol(i_),log10(vel(i)),tau(i),sigmae(i),mu(i),disp(i),psi(i),x,k
!   end do
!   write(nout(1),*)
!   write(nout(1),*)
! end subroutine
subroutine find_locid(i,j,loc_)
  integer::loc
  integer,intent(in)::i,j
  integer,intent(out)::loc_
  loc_=(i-1)*jmax+j
  return
end subroutine


  !------------initond-----------------------------------------------------------!
subroutine initcond_injection()
  implicit none
  real(8)::rr,rand
  phi=phi0
  kp=kpmin
  ks=kp
  kTv=kT
  pf=pfinit
  sigma=sigmainit
  pfhyd=0.d0
  pf=0d0
  ! !pf(ncellg/2)=2.5
  ! do i=1,ncellg
  !   i_=st_sum%lodc(i)
  !   rr=xcol(i_)**2+zcol(i_)**2
  !   if(rr<0.04**2) pf(i)=1.0
  ! end do

  sigmae=sigma-pf
  vel=velinit
  tau=sigmae*muinit
  ! vel=tau/abs(tau)*velinit
  mu=muinit
  psi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))

  !randomize initial state
  ! do i=1,ncellg
  !   call random_number(rand)
  !   psi(i)=psi(i)+(rand-0.5)*0.1
  !   vel(i)=2*vref*exp(-psi(i)/a(i))*sinh(tau(i)/sigmae(i)/a(i))
  ! end do

  !random initial shear stress
  ! do i=1,ncellg
  !   call random_number(rand)
  !   tau(i)=sigmae(i)*(muinit+(rand-0.5)*0.0)
  ! end do
  ! mu=tau/sigmae
  ! psi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))

  disp=0d0
end subroutine

subroutine initcond_bgflow()
  implicit none
  lf=imax*ds0

  !q0=(pbcr-pbcl)/lf*ks(1)/1d-12*1d-6*2

  !kp0=kp0*exp(sigmainit/s0)
  kpmax=kp0*(1+kL/kT/Vpl)-kpmin*kL/kT/Vpl

  kTv=kT
  kp=kp0
  ks=kp0
  pbcr=0d0
  pbcl=pbcr+q0*lf/ks(1)*eta*1d-3
  phi=phi0

  pf=pbcl+(pbcr-pbcl)*xcol/lf
  sigma=sigmainit+pf
  sigmae=sigma-pf
  pfhyd=0.d0
  vel=velinit
  ! do i=1,ncell
  !   if(abs(xcol(i)-5.0)<1.0) vel(i)=velinit*(1.0+0.001*sin(2*pi/2.0*xcol(i)))
  ! end do
  tau=sigmae*(f0+(a-b)*log(vel/vref))*1.001
  ! vel=tau/abs(tau)*velinit
  mu=tau/sigmae
  psi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))
  disp=0d0
end subroutine

subroutine input_well()
  implicit none
  integer::k,kwell
  open(77,file=injection_file)
  read(77,*) nwell
  write(*,*) 'nwell',nwell
  allocate(iwell(nwell),jwell(nwell),qvals(nwell,50),qtimes(nwell,50),kleng(nwell))
  do kwell=1,nwell
    read(77,*) iwell(kwell),jwell(kwell),kleng(kwell)
    !write(*,*) iwell(kwell),jwell(kwell)
    read(77,*) qvals(kwell,1:kleng(kwell))
    read(77,*) qtimes(kwell,1:kleng(kwell))
    !write(*,*) kleng(kwell)
    !write(*,*) qvals(kwell,:)
    !write(*,*) qtimes(kwell,:)
  end do
  close(77)
end subroutine 


  !------------coordinate-----------------------------------------------------------!
subroutine coordinate2dp(NCELLg,ds0,xel,xer,xcol)
  implicit none
  integer,intent(in)::NCELLg
  real(8),intent(in)::ds0
  real(8),intent(out)::xel(:),xer(:),xcol(:)
  integer::i,j,k

  !flat fault with element size ds
  do i=1,NCELLg
    ds(i)=ds0
    xel(i)=(i-1)*ds0
    xer(i)=i*ds0
    xcol(i)=0.5d0*(xel(i)+xer(i))
    !write(14,'(3e16.6)') xcol(i),xel(i),xer(i)
  enddo
  !close(14)
  return
end subroutine

subroutine coordinate2dph()
  implicit none
  integer::i,j,k
  real(8)::ofs=0.0
  !planar fault with element size ds and dipangle=dipangle
  do i=1,NCELLg
    ds(i)=ds0
    xel(i)=(ofs+(i-1)*ds0)*cos(dipangle*pi/180)
    xer(i)=(ofs+i*ds0)*cos(dipangle*pi/180)
    yel(i)=(ofs+(i-1)*ds0)*sin(dipangle*pi/180)
    yer(i)=(ofs+i*ds0)*sin(dipangle*pi/180)
    xcol(i)=0.5d0*(xel(i)+xer(i))
    ycol(i)=0.5d0*(yel(i)+yer(i))
    ang(i)=datan2(yer(i)-yel(i),xer(i)-xel(i))
    write(*,*) xcol(i),ycol(i)
  enddo
  !close(14)
  return
end subroutine

subroutine coordinate2dn()
  implicit none
  integer::i
  do i=1,NCELLg
    ds(i)=sqrt((xer(i)-xel(i))**2+(yer(i)-yel(i))**2)
    ang(i)=datan2(yer(i)-yel(i),xer(i)-xel(i))
    xcol(i)=0.5d0*(xel(i)+xer(i))
    ycol(i)=0.5d0*(yel(i)+yer(i))
    write(*,*) ds(i),ang(i)
  enddo
  return
end subroutine

subroutine coordinate3dp(imax,jmax,ds0,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
  implicit none
  integer,intent(in)::imax,jmax
  real(8),intent(in)::ds0
  real(8),intent(out)::xcol(:),zcol(:)
  real(8),intent(out)::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
  real(8)::dx,dz
  integer::i,j,k

  dx=ds0
  dz=ds0
  do i=1,imax
    do j=1,jmax
      k=(i-1)*jmax+j
      xcol(k)=(i-imax/2-0.5d0)*dx
      zcol(k)=-(j-jmax/2-0.5d0)*dz
      xs1(k)=xcol(k)+0.5d0*dx
      xs2(k)=xcol(k)-0.5d0*dx
      xs3(k)=xcol(k)-0.5d0*dx
      xs4(k)=xcol(k)+0.5d0*dx
      zs1(k)=zcol(k)+0.5d0*dz
      zs2(k)=zcol(k)+0.5d0*dz
      zs3(k)=zcol(k)-0.5d0*dz
      zs4(k)=zcol(k)-0.5d0*dz
    end do
  end do
  return
end subroutine coordinate3dp

  subroutine coordinate3ddip(imax,jmax,ds0,dipangle)
    implicit none
    integer,intent(in)::imax,jmax
    real(8),intent(in)::ds0,dipangle
    !integer,intent(in)::NCELLg
    !real(8),intent(out)::xcol(:),ycol(:),zcol(:)
    !real(8),intent(out)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    integer::i,j,k
    real(8)::xc,yc,zc,amp,yr(0:imax),zr(0:imax),stangle


      !dipangle=dipangle*pi/180d0
      stangle=0d0*pi/180d0
       k=0
       yr(0)=0d0
       zr(0)=0d0
       !nb=int(50.0*jmax/320.0)
       !if(my_rank==0) write(*,*)nb
       do j=1,jmax
       yr(j)=yr(j-1)-ds0*cos(dipangle*pi/180d0)
       zr(j)=zr(j-1)-ds0*sin(dipangle*pi/180d0)
       end do
       do i=1,imax
         do j=1,jmax
         k=k+1
           xcol(k)=(i-imax/2-0.5d0)*ds0
           ycol(k)=0.5d0*(yr(j-1)+yr(j)) !-(j-0.5d0)*ds0*cos(dipangle)
           zcol(k)=0.5d0*(zr(j-1)+zr(j)) !-(j-0.5d0)*ds0*sin(dipangle)

           angd(k)=datan2(zr(j-1)-zr(j),yr(j-1)-yr(j))
           ang(k)=0d0
           !write(*,*)angd(k)
           if(my_rank==0)write(111,*)xcol(k),ycol(k),zcol(k)
           end do
       end do

    return
  end subroutine coordinate3ddip

  subroutine varscalc(NCELL,displs,vars)
    implicit none
    integer,intent(in)::NCELL,displs(:)
    integer,intent(out)::vars(:)
    do i=1,NCELL
      vars(i)=displs(i-1)+i
      !write(*,*) my_rank,i,vars(i)
    end do
    return
  end subroutine
  subroutine evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33,ds)
    !calculate ev for each element
    implicit none
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    real(8),intent(out)::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:),ds(:)
    real(8)::rr,vba(0:2),vca(0:2),tmp1,tmp2,tmp3

    do k=1,NCELLg
      vba(0) = xs2(k)-xs1(k)
      vba(1) = ys2(k)-ys1(k)
      vba(2) = zs2(k)-zs1(k)
      vca(0) = xs3(k)-xs1(k)
      vca(1) = ys3(k)-ys1(k)
      vca(2) = zs3(k)-zs1(k)

      ev31(k) = vba(1)*vca(2)-vba(2)*vca(1)
      ev32(k) = vba(2)*vca(0)-vba(0)*vca(2)
      ev33(k) = vba(0)*vca(1)-vba(1)*vca(0)
      rr = sqrt(ev31(k)*ev31(k)+ev32(k)*ev32(k)+ev33(k)*ev33(k))
      !// unit vectors for local coordinates of elements
      ev31(k) = ev31(k)/rr ; ev32(k) = ev32(k)/rr ; ev33(k) = ev33(k)/rr
      !if(my_rank==0) write(*,'(i0,3e15.6)') k,ev31(k),ev32(k),ev33(k)

      if( abs(ev33(k)) < 1.0d0 ) then
        ev11(k) = -ev32(k) ; ev12(k) = ev31(k) ; ev13(k) = 0.0d0
        rr = sqrt(ev11(k)*ev11(k) + ev12(k)*ev12(k))
        ev11(k) = ev11(k)/rr ; ev12(k) = ev12(k)/rr;
      else
        ev11(k) = 1.0d0 ; ev12(k) = 0.0d0 ; ev13(k) = 0.0d0
      end if
      !if(my_rank==0) write(*,*) ev11(k),ev12(k),ev13(k)

      ev21(k) = ev32(k)*ev13(k)-ev33(k)*ev12(k)
      ev22(k) = ev33(k)*ev11(k)-ev31(k)*ev13(k)
      ev23(k) = ev31(k)*ev12(k)-ev32(k)*ev11(k)

      tmp1=vba(0)*vba(0)+vba(1)*vba(1)+vba(2)*vba(2)
      tmp2=vca(0)*vca(0)+vca(1)*vca(1)+vca(2)*vca(2)
      tmp3=vba(0)*vca(0)+vba(1)*vca(1)+vba(2)*vca(2)
      ds(k)=0.5d0*sqrt(tmp1*tmp2-tmp3*tmp3)
      !if(my_rank==0) write(*,*)ev21(k),ev22(k),ev23(k)
    end do

  end subroutine
  subroutine implicitsolver2(pf,sigma,ks,h,time,dtnxt,niter)
    implicit none
    integer::kit,errloc(1),l,l_
    integer,intent(out)::niter
    real(8),intent(inout)::pf(:),dtnxt
    real(8),intent(in)::h,time,sigma(:),ks(:)
    real(8)::dpf(ncellg),pfd(imax,jmax),pfnew(imax,jmax),sigmae(ncellg),cdiff(imax,jmax),err,err0,adpf(ncellg)
    real(8),parameter::dpth=0.1
    !real(8),parameter::cc=1d-12 !beta(1e-8)*phi(1e-1)*eta(1e-3)
    cdiff=0d0

    do l=1,ncellg
      !sigma(l)=sigmae(l)+pf(l)
      i=(l-1)/jmax+1
      j=l-(i-1)*jmax
      cc=eta*beta*phi0
      pfd(i,j)=pf(l)
      cdiff(i,j)=ks(l)/cc*1d-6 !Pa-->MPa
      !cdiff(i,j)=kpmax/cc*1d-6 !Pa-->MPa
      !write(*,*) l_,i,j
    end do
    !write(*,*) pfd

    err=0d0
      call Beuler(pfd,cdiff,h,pfnew,time,niter)

      do l=1,ncellg
        i=(l-1)/jmax+1
        j=l-(i-1)*jmax
        dpf(l)=pfnew(i,j)-pf(l)
        pf(l)=pfnew(i,j)
      end do
      !write(*,*)sum(pf)

    !write(*,*) 'niter',niter
      adpf=abs(dpf)
    !write(*,*) 'dpf',maxval(adpf)
    if(dtnxt/h*maxval(adpf)>dpth)  dtnxt=dpth*h/maxval(adpf)
    return
  end subroutine

  subroutine Beuler(pf,cdiff,h,pfnew,time,niter)
    implicit none
    integer,parameter::itermax=1000
    real(8),intent(in)::pf(:,:),h,cdiff(:,:),time
    real(8),intent(out)::pfnew(:,:)
    integer,intent(out)::niter
    real(8)::Dxx(imax,jmax,3),Dyy(imax,jmax,3),Amx(imax,jmax,3),Amy(imax,jmax,3),mx(imax,jmax),my(imax,jmax)
    real(8)::p(imax,jmax),m(imax,jmax),r(imax,jmax),x(imax,jmax),b(imax,jmax),SAT(imax,jmax)
    integer::n,iter,i,j,k,kwell
    real(8)::p0=0.0,rsold,rsnew,tmp1,tmp2,alpha,v1,v0,t1,t0,qtmp
    real(8),parameter::tol=1e-4
    !real(8),parameter::str=1e-11 !beta(1e-9)*phi(1e-2)
    niter=0
    p=0d0;m=0d0;r=0d0;x=0d0;b=0d0

    Dxx=0d0 
    do i=1,imax
        !compute Dxx for Dirichlet BC
        select case(bc)
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

        select case(bc)
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
        select case(bc)
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

        select case(bc)
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

    SAT(1,:)=-q0/beta/phi0*1e-9*h/ds0*2
    SAT(imax,:)=q0/beta/phi0*1e-9*h/ds0*2


    x=pf !initial guess
    b=pf-SAT

    if(injectionfromfile) then
      qtmp=0d0
      do kwell=1,nwell
        do k=1,kleng(kwell)-1
        t0=qtimes(kwell,k)
        t1=qtimes(kwell,k+1)
        v0=qvals(kwell,k)
        v1=qvals(kwell,k+1)
        if (time >= t0 .and. time <= t1) then
          qtmp=(v1-v0)/(t1-t0)*(time-t0)+v0
        else if (time > t1.and. k == kleng(kwell)-1) then
          qtmp=v1
        end if
        end do
        i=iwell(kwell)
        j=jwell(kwell)
        !write(*,*)i,j,qtmp
        b(i,j)=b(i,j)+h*qtmp/beta/phi0*1e-12/ds0/ds0
      end do
    else if(time<td) then
      b(imax/2,jmax/2)=b(imax/2,jmax/2)+h*q0/beta/phi0*1e-12/ds0/ds0
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
      100 pfnew=x

    return
    end subroutine

!computing dydx for time integration
    subroutine derivs(x, y, dydx)
      use m_HACApK_solve
      use m_HACApK_base
      use m_HACApK_use
      implicit none
      real(8),intent(in) :: x
      real(8),intent(in) ::y(:)
      real(8),intent(out) :: dydx(:)
      real(8) :: veltmp(NCELL),tautmp(NCELL),sigmaetmp(NCELL),psitmp(NCELL),kstmp(NCELL)
      real(8) :: sum_gs(NCELL),sum_gn(NCELL),veltmpG(NCELLg),sum_gsg(NCELLg),sum_gng(NCELLg)
      real(8) :: time3,time4,c1, c2, c3, fss,arg,arg2,c,g,tauss,Arot(3,3),p(6),fac
      real(8)::dvdsig,sxx0,sxy0,syy0,dvdpsi,dpsidt,dsigdt,dtaudt,dvdtau,dksdt
      integer :: i, j, nc,ierr,lrtrn,i_
      real(8)::Atide=1e-2,Ttide=24*3600d0
  
      !if(latticeh) then
  
      !$omp parallel do
      do i = 1, NCELL
        i_=vars(i)
        psitmp(i) = y(4*i-3)
        tautmp(i) = y(4*i-2)
        sigmaetmp(i) = y(4*i-1)
        kstmp(i)=y(4*i)
        veltmp(i) = 2*vref*dexp(-psitmp(i)/a(i_))*dsinh(tautmp(i)/sigmaetmp(i)/a(i_))
        !if(abs(i-ncellg/2)<5) write(*,*) psitmp(i)
      enddo
      !$omp end parallel do
     ! write(*,*) y(1489*4-1)
     ! write(*,*) veltmp(1489),tautmp(1489),sigmaetmp(1489)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      sum_gn=0d0

      call MPI_ALLGATHERv(veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
      if(backslip) then
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG-vpl)
      else
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
      end if
      call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     

    !$omp parallel do
      do i=1,NCELL
        i_=vars(i)
        !permeability evolution
        dksdt=-veltmp(i)/kL*(kstmp(i)-kpmax)-(kstmp(i)-kpmin)/kTv(i_)
  
        !aging law
        select case(evlaw)
        case('aging')
          if(b(i_)==0) then
            dpsidt=0d0
          else 
            dpsidt=b(i_)/dc(i_)*vref*dexp((f0(i_)-psitmp(i))/b(i_))-b(i_)*veltmp(i)/dc(i_)
          end if
        !slip law
        case('slip')
          fss=f0(i_)+(a(i_)-b(i_))*dlog(abs(veltmp(i))/vref)
          dpsidt=-abs(veltmp(i))/dc(i_)*(abs(tautmp(i))/sigmaetmp(i)-fss)
        end select
  
        dsigdt=sum_gn(i)
        !dsigdt=dsigdt-cd*dphidt
  
        dvdtau=2*vref*dexp(-psitmp(i)/a(i_))*dcosh(tautmp(i)/sigmaetmp(i)/a(i_))/(a(i_)*sigmaetmp(i))
        dvdsig=-2*vref*dexp(-psitmp(i)/a(i_))*dcosh(tautmp(i)/sigmaetmp(i)/a(i_))*tautmp(i)/(a(i_)*sigmaetmp(i)**2)
        dvdpsi=-veltmp(i)/a(i_)
        dtaudt=sum_gs(i)-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt)
        !if(i==NCELL/2) dtaudt=dtaudt-tautmp(i)/sigmatmp(i)*q0
        dtaudt=dtaudt/(1d0+0.5d0*rigid/vs*dvdtau)
  
        dydx(4*i-3)=dpsidt
        dydx(4*i-2)=dtaudt
        dydx(4*i-1)=dsigdt
        dydx(4*i)=dksdt
        !dydx(4*i)=0d0
       ! if(i==1489) write(*,*) dtaudt,dsigdt
        !call deriv_d(sum_gs(i),sum_gn(i),phitmp(i),tautmp(i),sigmatmp(i),veltmp(i),a(i),b(i),dc(i),f0(i),dydx(3*i-2),dydx(3*i-1),dydx(3*i))
      enddo
        !$omp end parallel do
      return
    end subroutine

 !---------------------------------------------------------------------
    subroutine rkqs(y,dydx,x,htry,eps,hdid,hnext,errmax_gb,errmaxloc)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
      !---------------------------------------------------------------------
      use m_HACApK_solve
      use m_HACApK_base
      use m_HACApK_use
      implicit none
      !include 'mpif.h'
      !integer::NCELL,NCELLg,rcounts(:),displs(:)
      real(8),intent(in)::htry,eps
      real(8),intent(inout)::y(:),x,dydx(:)
      real(8),intent(out)::hdid,hnext,errmax_gb !hdid: resulatant dt hnext: htry for the next
      integer,intent(out)::errmaxloc
      !type(st_HACApK_lcontrol),intent(in) :: st_ctl
      !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
      !type(st_HACApK_calc_entry) :: st_bemv
      integer :: i,ierr,loc
      real(8) :: errmax,h,xnew,htemp,dtmin
      real(8),dimension(size(y))::yerr,ytemp
      real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4
  
      h=htry
      !dtmin=0.5d0*minval(ds)/vs
      !call derivs(x,y,dydx)
      do while(.true.)
  
        call MPI_BARRIER(MPI_COMM_WORLD,ierr);time3=MPI_Wtime()
        call rkck(y,x,h,ytemp,yerr)
        !call rk2(y,x,h,ytemp,yerr)
        time4=MPI_Wtime()
        timeH=timeH+time4-time3
  
        errmax=0d0
        !do i=1,NCELL
        !  if(abs(yerr(3*i-2)/ytemp(3*i-2))/eps>errmax) errmax=abs(yerr(3*i-2)/ytemp(3*i-2))/eps
          !errmax=errmax+yerr(3*i-2)**2
        !end do
        errmaxloc=0
        ! do i=1,ncell
        !   if(abs(yerr(4*i-2)/ytemp(4*i-2))/eps>errmax) then
        !     errmax=abs(yerr(4*i-2)/ytemp(4*i-2))/eps
        !     errmaxloc=i
        !   end if!errmax=errmax+yerr(3*i-2)**2
        ! end do
  
        do i=1,4*ncell
          if(abs(yerr(i)/ytemp(i))/eps>errmax) then
            errmax=abs(yerr(i)/ytemp(i))/eps
            errmaxloc=i
          end if!errmax=errmax+yerr(3*i-2)**2
        end do

        call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  
        !if(my_rank==0)write(*,*) h,errmax_gb
        !if(h<0.25d0*ds0/vs)exit
        if((errmax_gb<1.d0).and.(errmax_gb>1d-15)) then
          exit
        end if
  
        if(errmax_gb>1d-15) then
          h=max(0.5d0*h,SAFETY*h*(errmax_gb**PSHRNK))
        else
          h=0.5*h
        end if
  
  
  
        xnew=x+h
        if(xnew-x<1.d-15) then
          if(my_rank.eq.0)write(*,*)'error: dt is too small'
          stop
        end if
  
      end do
  
      hnext=min(2*h,SAFETY*h*(errmax_gb**PGROW))
      if(outpertime) then
        hnext=min(hnext,dtout*365*24*3600)
      end if
      !if(load==0)hnext=min(hnext,dtmax)
      !hnext=max(0.249d0*ds0/vs,SAFETY*h*(errmax_gb**PGROW))
  
      !hnext=min(,1d9)
  
      hdid=h
      x=x+h
      y(:)=ytemp(:)
  
    end subroutine

 
    !---------------------------------------------------------------------
  subroutine rkck(y,x,h,yout,yerr)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    !include 'mpif.h'
    !integer,intent(in)::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::y(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    !integer,intent(out)::ierr
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    integer ::i
    real(8) :: ak1(4*NCELL),ak2(4*NCELL),ak3(4*NCELL),ak4(4*NCELL),ak5(4*NCELL),ak6(4*NCELL),ytemp(4*NCELL)
    real(8) :: A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51
    real(8) :: B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3./40.)
    parameter (B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5)
    parameter (B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.)
    parameter (B63=575./13824.,B64=44275./110592.,B65=253./4096.)
    parameter (C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.)
    parameter (DC1=C1-2825./27648.,DC3=C3-18575./48384.)
    parameter (DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)
    !ierr=0
    !     -- 1st step --
    call derivs(x, y, ak1)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+B21*h*ak1(i)
    end do
    !$omp end parallel do

    !    -- 2nd step --
    call derivs(x+a2*h, ytemp, ak2)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B31*ak1(i)+B32*ak2(i))
    end do
    !$omp end parallel do

    !     -- 3rd step --
    call derivs(x+a3*h, ytemp, ak3)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B41*ak1(i)+B42*ak2(i)+B43*ak3(i))
    end do
    !$omp end parallel do

    !     -- 4th step --
    call derivs(x+a4*h, ytemp, ak4)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B51*ak1(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
    end do
    !$omp end parallel do

    !     -- 5th step --
    call derivs(x+a5*h, ytemp, ak5)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B61*ak1(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do
    !$omp end parallel do

    !     -- 6th step --
    call derivs(x+a6*h, ytemp, ak6)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      yout(i)=y(i)+h*(C1*ak1(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
    end do
    !$omp end parallel do


    !$omp parallel do
    do i=1,size(y)
      yerr(i)=h*(DC1*ak1(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
      !if(abs(yerr(i))>=1d6)ierr=1
    end do
    !$omp end parallel do

    return
  end subroutine

  subroutine foward_check()
    implicit none
    real(8)::rr,lc
    integer::p

    ! vel=0d0
    ! lc=0.3d0
    ! do i=1,NCELLg
    !   rr=ycol(i)**2+zcol(i)**2
    !   if(rr<lc**2) vel(i)=5d0/rigid*sqrt(lc**2-rr)
    ! end do
    vel=1d0
    !vel(21299:)=1d0
    !vel(1)=1d0! p=532
    ! vel(p)=1d0


    write(fname,'("stress",i0)') number
    open(29,file=fname)

    select case(problem)
    case('2dn','25d')
      !slip from file
      ! open(45,file='../fd2d/rupt2.dat')
      ! do i=1,NCELLg
      !   read(45,*) a(i),vel(i),b(i)
      ! end do

      st_bemv%v='xx'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xx,st_bemv,st_ctl,a,vel)
      st_bemv%v='xy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xy,st_bemv,st_ctl,b,vel)
      st_bemv%v='yy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yy,st_bemv,st_ctl,dc,vel)
      if(my_rank.eq.0) then
        do i=1,NCELLg
          taudot(i)=0.5d0*(a(i)-dc(i))*dsin(-2*ang(i))+b(i)*dcos(-2*ang(i))
          sigdot(i)=-(0.5d0*(a(i)+dc(i))-0.5d0*(a(i)-dc(i))*dcos(2*ang(i))-b(i)*dsin(2*ang(i)))
          write(29,'(4e16.4)') xcol(i),ang(i),taudot(i),sigdot(i)
        end do
      end if
    case('3dp','3dph')
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,a,vel)

      if(my_rank.eq.0) then
        do i=1,NCELLg
          write(29,'(3e16.4)') xcol(i),zcol(i),a(i)
        end do
      end if
    case('3dn','3dh')
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_s2,st_bemv,st_ctl,a,vel)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_d2,st_bemv,st_ctl,b,vel)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n2,st_bemv,st_ctl,dc,vel)
      if(my_rank.eq.0) then
        do i=1,NCELLg
          write(29,'(6e16.4)') xcol(i),ycol(i),zcol(i),a(i),b(i),dc(i)
        end do
      end if
    case('3dnf','3dhf')
      st_bemv%md='st'
      st_bemv%v='s'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_s,st_bemv,st_ctl,a,vel)
      st_bemv%v='n'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n,st_bemv,st_ctl,dc,vel)
      if(my_rank.eq.0) then
        do i=1,NCELLg
          write(29,'(6e16.4)') xcol(i),ycol(i),zcol(i),a(i),b(i),dc(i)
        end do
      end if
    case('2dnh')
      st_bemv%v='xx'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xx,st_bemv,st_ctl,a,vel)
      st_bemv%v='xy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xy,st_bemv,st_ctl,b,vel)
      st_bemv%v='yy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yy,st_bemv,st_ctl,dc,vel)
      if(my_rank.eq.0) then
        do i=1,NCELLg
          taudot(i)=0.5d0*(a(i)-dc(i))*dsin(-2*ang(i))+b(i)*dcos(-2*ang(i))
          sigdot(i)=-(0.5d0*(a(i)+dc(i))-0.5d0*(a(i)-dc(i))*dcos(2*ang(i))-b(i)*dsin(2*ang(i)))
          write(29,'(4e16.4)') xcol(i),ycol(i),taudot(i),sigdot(i)
        end do
      end if
    end select
    Call MPI_FINALIZE(ierr)
    stop
  end subroutine

  subroutine inverse_problem()
    write(*,*) 'slip from stress drop'
    write(fname,'("stress",i0)') number
    open(29,file=fname)

    select case(problem)
    case('2dp')
      taudot=-1d0
      lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,eps_h)
      !lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sigdot,taudot)
      lrtrn=HACApK_gensolv(st_leafmtxp_c,st_bemv,st_ctl,coord,taudot,sigdot,eps_h)
      if(my_rank.eq.0) then
        do i=1,NCELLg
          write(29,'(2e16.4)') xcol(i),sigdot(i)
        end do
      end if
    case('3dhf')
      do i=1,ncellg
        taudot(i)=-1d0
      end do
      st_bemv%v='s'
      st_bemv%md='st'
      lrtrn=HACApK_generate(st_leafmtxp_c,st_bemv,st_ctl,coord,eps_h)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_c,st_bemv,st_ctl,sigdot,taudot)
      !lrtrn=HACApK_gensolv(st_leafmtxp_c,st_bemv,st_ctl,coord,taudot,sigdot,eps_h)
      if(my_rank.eq.0) then
        do i=1,NCELLg
          write(29,'(4e16.4)') xcol(i),ycol(i),zcol(i),sigdot(i)
        end do
      end if
    end select
    Call MPI_FINALIZE(ierr)
    stop
  end subroutine

!   function rtnewt(prev,eps,nst,p,t0,sum)
!     integer::j
!     integer,parameter::jmax=20
!     real(8)::rtnewt,prev,eps
!     real(8)::f,df,dx,sum,nst,p,t0
!     rtnewt=prev
!     !write(*,*) rtnewt
!     do j=1,jmax
!       x=rtnewt
!       f=x+ieta*sigma0*(mu0+(a0-b0)*log(x/vref))-vpl
!       df=1+ieta*sigma0*(a0-b0)/x
!       dx=f/df
!       rtnewt=rtnewt-dx
!       !write(*,*) rtnewt
!       if(abs(dx).lt.eps) return
!     end do
!     write(*,*) 'maximum iteration'
!     stop
!   end function
! end program

! subroutine open_bp(problem)
!   character(128),intent(in)::problem
!   real(8)::xd(81)
!   select case(problem)
!   !SEAS BP5
!   case('3dph')
!     open(101,file="output/fltst_strk-36dp+00")
!     open(102,file="output/fltst_strk-16dp+00")
!     open(103,file="output/fltst_strk+00dp+00")
!     open(104,file="output/fltst_strk+16dp+00")
!     open(105,file="output/fltst_strk+36dp+00")
!     open(106,file="output/fltst_strk-24dp+10")
!     open(107,file="output/fltst_strk-16dp+10")
!     open(108,file="output/fltst_strk+00dp+10")
!     open(109,file="output/fltst_strk+16dp+10")
!     open(110,file="output/fltst_strk+00dp+22")
!     do i=101,110
!       write(i,*)"# This is the header:"
!       write(i,*)"# problem=SEAS Benchmark BP5-QD"
!       write(i,*)"# code=hbi"
!       write(i,*)"# modeler=So Ozawa"
!       write(i,*)"# date=2021/03/19"
!       write(i,*)"# element_size=500m"
!       write(i,*)"# Column #1 = Time (s)"
!       write(i,*)"# Column #2 = Slip_2(m)"
!       write(i,*)"# Column #3 = Slip_3(m)"
!       write(i,*)"# Column #4 = Slip_rate_2(log10 m/s)"
!       write(i,*)"# Column #5 = Slip_rate_3(log10 m/s)"
!       write(i,*)"# Column #6 = Shear_stress_2 (MPa)"
!       write(i,*)"# Column #7 = Shear_stress_3 (MPa)"
!       write(i,*)"# Column #8 = State (log10 s)"
!       write(i,*)"# The line below lists the names of the data fields"
!       write(i,*)"t slip_2 slip_3 slip_rate_2 slip_rate_3 shear_stress_2 shear_stress_3 state"
!       write(i,*)"# Here is the time-series data."
!     end do

!     open(120,file="output/global.dat")
!     i=120
!     write(i,*)"# This is the file header:"
!     write(i,*)"# problem=SEAS Benchmark BP4-QD"
!     write(i,*)"# code=hbi"
!     write(i,*)"# modeler=So Ozawa"
!     write(i,*)"# date=2021/03/19"
!     write(i,*)"# element_size=500m"
!     write(i,*)"# Column #1 = Time (s)"
!     write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
!     write(i,*)"# Column #3 = Moment rate (N-m/s)"
!     write(i,*)"# The line below lists the names of the data fields"
!     write(i,*)"t max_slip_rate moment_rate"
!     write(i,*)"# Here is the time-series data."

!     open(130,file="output/rupture.dat")
!     i=130
!     write(i,*)"# This is the file header:"
!     write(i,*)"# problem=SEAS Benchmark BP4-QD"
!     write(i,*)"# code=hbi"
!     write(i,*)"# modeler=So Ozawa"
!     write(i,*)"# date=2021/03/19"
!     write(i,*)"# element_size=500m"
!     write(i,*)"# Column #1 = x2 (m)"
!     write(i,*)"# Column #2 = x3 (m)"
!     write(i,*)"# Column #3 = time (s)"
!     write(i,*)"# The line below lists the names of the data fields"
!     write(i,*)"x2 x3 t"
!     write(i,*)"# Here is the data."

!   !SEAS BP4
!     case('3dp')
!     open(101,file="output/fltst_strk-360dp+000")
!     open(102,file="output/fltst_strk-225dp-750")
!     open(103,file="output/fltst_strk-165dp-120")
!     open(104,file="output/fltst_strk-165dp+000")
!     open(105,file="output/fltst_strk-165dp+120")
!     open(106,file="output/fltst_strk+000dp-210")
!     open(107,file="output/fltst_strk+000dp-120")
!     open(108,file="output/fltst_strk+000dp+000")
!     open(109,file="output/fltst_strk+000dp+120")
!     open(110,file="output/fltst_strk+000dp+210")
!     open(111,file="output/fltst_strk+165dp-120")
!     open(112,file="output/fltst_strk+165dp+000")
!     open(113,file="output/fltst_strk+165dp+120")
!     open(114,file="output/fltst_strk+360dp+000")
!     do i=101,114
!       write(i,*)"# This is the header:"
!       write(i,*)"# problem=SEAS Benchmark BP4-QD"
!       write(i,*)"# code=hbi"
!       write(i,*)"# modeler=So Ozawa"
!       write(i,*)"# date=2021/03/19"
!       write(i,*)"# element_size=500m"
!       write(i,*)"# Column #1 = Time (s)"
!       write(i,*)"# Column #2 = Slip_2(m)"
!       write(i,*)"# Column #3 = Slip_3(m)"
!       write(i,*)"# Column #4 = Slip_rate_2(log10 m/s)"
!       write(i,*)"# Column #5 = Slip_rate_3(log10 m/s)"
!       write(i,*)"# Column #6 = Shear_stress_2 (MPa)"
!       write(i,*)"# Column #7 = Shear_stress_3 (MPa)"
!       write(i,*)"# Column #8 = State (log10 s)"
!       write(i,*)"# The line below lists the names of the data fields"
!       write(i,*)"t slip_2 slip_3 slip_rate_2 slip_rate_3 shear_stress_2 shear_stress_3 state"
!       write(i,*)"# Here is the time-series data."
!     end do

!     open(120,file="output/global.dat")
!     i=120
!     write(i,*)"# This is the file header:"
!     write(i,*)"# problem=SEAS Benchmark BP4-QD"
!     write(i,*)"# code=hbi"
!     write(i,*)"# modeler=So Ozawa"
!     write(i,*)"# date=2021/03/19"
!     write(i,*)"# element_size=500m"
!     write(i,*)"# Column #1 = Time (s)"
!     write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
!     write(i,*)"# Column #3 = Moment rate (N-m/s)"
!     write(i,*)"# The line below lists the names of the data fields"
!     write(i,*)"t max_slip_rate moment_rate"
!     write(i,*)"# Here is the time-series data."

!     open(130,file="output/rupture.dat")
!     i=130
!     write(i,*)"# This is the file header:"
!     write(i,*)"# problem=SEAS Benchmark BP4-QD"
!     write(i,*)"# code=hbi"
!     write(i,*)"# modeler=So Ozawa"
!     write(i,*)"# date=2021/03/19"
!     write(i,*)"# element_size=500m"
!     write(i,*)"# Column #1 = x2 (m)"
!     write(i,*)"# Column #2 = x3 (m)"
!     write(i,*)"# Column #3 = time (s)"
!     write(i,*)"# The line below lists the names of the data fields"
!     write(i,*)"x2 x3 t"
!     write(i,*)"# Here is the data."

!     !SEAS BP3
!     case('2dnh')
!     open(101,file="output/fltst_dp000")
!     open(102,file="output/fltst_dp025",status='replace')
!     open(103,file="output/fltst_dp050",status='replace')
!     open(104,file="output/fltst_dp075",status='replace')
!     open(105,file="output/fltst_dp100",status='replace')
!     open(106,file="output/fltst_dp125",status='replace')
!     open(107,file="output/fltst_dp150",status='replace')
!     open(108,file="output/fltst_dp175",status='replace')
!     open(109,file="output/fltst_dp200",status='replace')
!     open(110,file="output/fltst_dp250",status='replace')
!     open(111,file="output/fltst_dp300",status='replace')
!     open(112,file="output/fltst_dp350",status='replace')
!     do i=101,112
!       write(i,*)"# This is the header:"
!       write(i,*)"# problem=SEAS Benchmark BP3-QD"
!       write(i,*)"# code=hbi"
!       write(i,*)"# modeler=So Ozawa"
!       write(i,*)"# date=2021/01/22"
!       write(i,*)"# element_size=25m"
!       write(i,*)"# location= on fault, 0km down-dip distance"
!       write(i,*)"# Column #1 = Time (s)"
!       write(i,*)"# Column #2 = Slip (m)"
!       write(i,*)"# Column #3 = Slip rate (log10 m/s)"
!       write(i,*)"# Column #4 = Shear stress (MPa)"
!       write(i,*)"# Column #5 = Normal stress (MPa)"
!       write(i,*)"# Column #6 = State (log10 s)"
!       write(i,*)"# The line below lists the names of the data fields"
!       write(i,*)"t slip slip_rate shear_stress normal_stress state"
!       write(i,*)"# Here is the time-series data."
!     end do
!     open(121,file="output/slip.dat",status='replace')
!     open(122,file="output/shear_stress.dat",status='replace')
!     open(123,file="output/normal_stress.dat",status='replace')

!     do i=121,123
!     write(i,*)"# This is the file header:"
!     write(i,*)"# problem=SEAS Benchmark BP3-QD"
!     write(i,*)"# code=hbi"
!     write(i,*)"# modeler=So Ozawa"
!     write(i,*)"# date=2021/03/16"
!     write(i,*)"# element_size=25m"
!     write(i,*)"# Column #1 = Time (s)"
!     write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
!     end do

!     write(121,*)"# Column #3-83 = Slip (m)"
!     write(122,*)"# Column #3-83 = Shear stress (MPa)"
!     write(123,*)"# Column #3-83 = Normal stress (MPa)"

!     do i=121,123
!     write(i,*)"# The line below lists the names of the data fields"
!     write(i,*)"xd"
!     end do
!     write(121,*)"t max_slip_rate slip"
!     write(122,*)"t max_slip_rate shear_stress"
!     write(123,*)"t max_slip_rate normal_stress"
!     do i=121,123
!     write(i,*)"# Here are the data."
!     end do
!     do i=1,81
!       xd(i)=(i-1)*500d0
!     end do
!     write(121,'(83e22.14)') 0d0,0d0,xd
!     write(122,'(83e22.14)') 0d0,0d0,xd
!     write(123,'(83e22.14)') 0d0,0d0,xd
!   end select
! return
! end subroutine
! subroutine debug()

! end subroutine
end program
