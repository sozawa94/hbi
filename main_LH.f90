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
  integer::NCELL,NCELLg,NSTEP1
  integer::imax,jmax !for 3dp

  !for HACApK
  real(8),allocatable ::coord(:,:),vmax(:)
  real(8)::eps_h
  type(st_HACApK_lcontrol) :: st_ctl
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
  real(8)::time1,time2,time3,time4,timer

  !for fault geometry
  real(8),allocatable::xcol(:),ycol(:),zcol(:),ds(:)
  real(8),allocatable::xs1(:),xs2(:),xs3(:),xs4(:) !for 3dp
  real(8),allocatable::zs1(:),zs2(:),zs3(:),zs4(:) !for 3dp
  real(8),allocatable::ys1(:),ys2(:),ys3(:) !for 3dn
  real(8),allocatable::xel(:),xer(:),yel(:),yer(:),ang(:)
  real(8),allocatable::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:)

  !parameters for each elements
  real(8),allocatable::a(:),b(:),dc(:),f0(:),fw(:),vw(:),vc(:),taudot(:),tauddot(:),sigdot(:)

  !variables
  real(8),allocatable::phi(:),vel(:),tau(:),sigma(:),disp(:),mu(:),rupt(:),idisp(:),velp(:)
  real(8),allocatable::taus(:),taud(:),vels(:),veld(:),disps(:),dispd(:),rake(:)


  integer::lp,i,i_,j,k,m,counts,interval,lrtrn,nl,ios,nmain,rk,nout
  integer,allocatable::locid(:)
  integer::hypoloc(1),load,eventcount,thec,inloc,sw

  !controls
  logical::aftershock,buffer,nuclei,slipping,outfield,slipevery,limitsigma,dcscale,slowslip,slipfinal
  logical::nonuniformstress,backslip,sigmaconst,foward,inverse,geofromfile,melange,creep,SEAS
  character*128::fname,dum,law,input_file,problem,geofile,param,pvalue,slipmode,project
  real(8)::a0,b0,dc0,sr,omega,theta,dtau,tiny,moment,wid,normal,ieta,meanmu,meanmuG,meandisp,meandispG,moment0,mvel,mvelG
  real(8)::psi,vc0,mu0,onset_time,tr,vw0,fw0,velmin,muinit,intau,trelax,maxnorm,maxnormG,minnorm,minnormG
  real(8)::r,vpl,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxeventi,sparam,tmax,dtmax
  real(8)::alpha,ds0,amp,mui,velinit,phinit,velmax,maxsig,minsig,v1,dipangle,crake

  !temporal variable

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
  npd=int(sqrt(dble(np)))
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )

  if(my_rank==0) then
    write(*,*) '# of MPI cores', np
  end if
  !input file must be specified when running
  !example) mpirun -np 16 ./ha.out default.in
  call get_command_argument(1,input_file,status=stat)

  open(33,file=input_file,iostat=ios)
  !if(my_rank==0) write(*,*) 'input_file',input_file
  if(ios /= 0) then
    write(*,*) 'Failed to open inputfile'
    stop
  end if

  !get filenumber
  number=0
  if(input_file(1:2)=='in') then
    input_file=adjustl(input_file(7:))
    write(*,*) input_file
    read(input_file,*) number
    write(*,*) number
  end if
  time1=MPI_Wtime()

  !default parameters
  nmain=1000000
  eps_r=1d-5
  eps_h=1d-5
  velmax=1d7
  velmin=1d-16
  law='d'
  tmax=1d12
  nuclei=.false.
  slipevery=.false.
  foward=.false.
  inverse=.false.
  slipfinal=.false.
  nonuniformstress=.false.
  maxsig=300d0
  minsig=20d0
  amp=0d0
  vc0=1d6
  vw0=1d6
  fw0=0.3d0
  dtinit=1d0
  tp=86400d0
  trelax=1d18
  project="none"
  !number=0


  do while(ios==0)
    read(33,*,iostat=ios) param,pvalue
    !write(*,*) param,pvalue
    select case(param)
    case('problem')
      read (pvalue,*) problem
    case('NCELLg')
      read (pvalue,*) ncellg
    case('imax')
      read (pvalue,*) imax
    case('jmax')
      read (pvalue,*) jmax
    case('NSTEP1')
      read (pvalue,*) nstep1
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
    case('b')
      read (pvalue,*) b0
    case('dc')
      read (pvalue,*) dc0
    case('vw')
      read (pvalue,*) vw0
    case('fw')
      read (pvalue,*) fw0
    case('vc')
      read (pvalue,*) vc0
    case('mu0')
      read (pvalue,*) mu0
    case('load')
      read (pvalue,*) load
    case('sr')
      read (pvalue,*) sr
    case('vpl')
      read (pvalue,*) vpl
    case('interval')
      read (pvalue,*) interval
    case('geometry')
      read (pvalue,*) geofile
    case('velinit')
      read (pvalue,*) velinit
    case('muinit')
      read (pvalue,*) muinit
    case('phinit')
      read (pvalue,*) phinit
    case('psi')
      read (pvalue,*) psi
    case('dtinit')
      read (pvalue,*) dtinit
    case('intau')
      read (pvalue,*) intau
    case('inloc')
      read (pvalue,*) inloc
    case('sparam')
      read (pvalue,*) sparam
    case('tmax')
      read (pvalue,*) tmax
    case('eps_r')
      read (pvalue,*) eps_r
    case('eps_h')
      read (pvalue,*) eps_h
    case('amp')
      read(pvalue,*) amp
    case('wid')
      read(pvalue,*) wid
    case('fwid')
      read(pvalue,*) fwid
    case('nuclei')
      read (pvalue,*) nuclei
    case('slipevery')
      read (pvalue,*) slipevery
    case('slipfinal')
      read (pvalue,*) slipfinal
    case('limitsigma')
      read (pvalue,*) limitsigma
    case('buffer')
      read (pvalue,*) buffer
    case('aftershock')
      read(pvalue,*) aftershock
    case('nmain')
      read(pvalue,*) nmain
    case('slipmode')
      read(pvalue,*) slipmode
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
    case('dipangle')
      read(pvalue,*) dipangle
    case('trelax')
      read(pvalue,*) trelax
    case('nonuniformstress')
      read(pvalue,*) nonuniformstress
    case('npgl')
      read(pvalue,*) npgl
    case('project')
      read(pvalue,*) project
    case('crake')
      read(pvalue,*) crake
    end select
  end do
  close(33)
  !limitsigma=.true.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !MPI setting
  !NCELLg=2*NL*NL
  !if(problem=='3dp') then
  nmain=ncellg
  select case(problem)
  case('3dp','3dph')
    NCELLg=imax*jmax
    loc=loci*(imax-1)+locj
  ! case('3dnf','3dn')
  !   NCELLg=imax*jmax*2
  end select
  !end if
  !stop
  !call varscalc(NCELL,displs,vars)
  if(my_rank==0) then
    write(*,*) 'job number',number
    write(*,*) 'project',project
  end if

  !allocation
  allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg),ds(NCELLg))
  xcol=0d0;ycol=0d0;zcol=0d0

  select case(problem)
  case('2dp','2dh')
    allocate(xel(NCELLg),xer(NCELLg))
    xel=0d0;xer=0d0
  case('2dn','2dnh','2dn3','25d')
    allocate(ang(NCELLg),xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg))
    ang=0d0;xel=0d0;xer=0d0;yel=0d0;yer=0d0
  case('3dp','3dph')
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
    xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
    zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
  case('3dnf','3dhf')
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg))
    allocate(ys1(NCELLg),ys2(NCELLg),ys3(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg))
    allocate(ev11(NCELLg),ev12(NCELLg),ev13(NCELLg))
    allocate(ev21(NCELLg),ev22(NCELLg),ev23(NCELLg))
    allocate(ev31(NCELLg),ev32(NCELLg),ev33(NCELLg))
    xs1=0d0; xs2=0d0; xs3=0d0
    ys1=0d0; ys2=0d0; ys3=0d0
    zs1=0d0; zs2=0d0; zs3=0d0
  case('3dn','3dh')
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg))
    allocate(ys1(NCELLg),ys2(NCELLg),ys3(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg))
    allocate(ev11(NCELLg),ev12(NCELLg),ev13(NCELLg))
    allocate(ev21(NCELLg),ev22(NCELLg),ev23(NCELLg))
    allocate(ev31(NCELLg),ev32(NCELLg),ev33(NCELLg))
    xs1=0d0; xs2=0d0; xs3=0d0
    ys1=0d0; ys2=0d0; ys3=0d0
    zs1=0d0; zs2=0d0; zs3=0d0

  end select

  !allocate(vmax(NCELLg),vmaxin(NcELLg))

  !mesh generation (rectangular assumed)
  if(my_rank==0) write(*,*) 'Generating mesh'
  select case(problem)
  case('2dp','2dh')
    call coordinate2dp(NCELLg,ds0,xel,xer,xcol)
  case('2dnh')
    call coordinate2dnh()
  case('3dp')
    call coordinate3dp(imax,jmax,ds0,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
  case('3dph')
    call coordinate3dph(imax,jmax,ds0,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
  case('3dn','3dh','3dnf','3dhf','fdph_FP11')
    open(20,file=geofile)
    do i=1,NCELLg
      read(20,*) k,xs1(i),ys1(i),zs1(i),xs2(i),ys2(i),zs2(i),xs3(i),ys3(i),zs3(i),xcol(i),ycol(i),zcol(i)
    end do
    !zs1=zs1-0.01d0
    !zs2=zs2-0.01d0
    !zs3=zs3-0.01d0
    !zcol=zcol-0.01d0
    !call coordinate3dn(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    !call coordinate3dns(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    !call coordinate3dns2(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    call evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33)
  end select

  if(project=='2DBEND') call coordinate_2DBEND()
  if(project=='SEAMOUNT')call coordinate3dn_SEAMOUNT(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)

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
  select case(problem)
  case('2dp','2dh')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%problem=problem

  case('2dn','2dn3','2dnh')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg),st_bemv%ds(NCELLg))
    allocate(st_bemv%ycol(NCELLg),st_bemv%yel(NCELLg),st_bemv%yer(NCELLg),st_bemv%ang(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%ycol=ycol;st_bemv%yel=yel;st_bemv%yer=yer
    st_bemv%ang=ang; st_bemv%ds=ds
    st_bemv%problem=problem

  case('25d')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg),st_bemv%ds(NCELLg))
    allocate(st_bemv%ycol(NCELLg),st_bemv%yel(NCELLg),st_bemv%yer(NCELLg),st_bemv%ang(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%ycol=ycol;st_bemv%yel=yel;st_bemv%yer=yer
    st_bemv%ang=ang; st_bemv%ds=ds
    st_bemv%problem=problem
    st_bemv%w=fwid

  case('3dp','3dph')
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

  case('3dn','3dh','3dnf','3dhf','fdph_FP11')
    allocate(st_bemv%xcol(NCELLg),st_bemv%ycol(NCELLg),st_bemv%zcol(NCELLg))
    allocate(st_bemv%xs1(NCELLg),st_bemv%xs2(NCELLg),st_bemv%xs3(NCELLg))
    allocate(st_bemv%ys1(NCELLg),st_bemv%ys2(NCELLg),st_bemv%ys3(NCELLg))
    allocate(st_bemv%zs1(NCELLg),st_bemv%zs2(NCELLg),st_bemv%zs3(NCELLg))
    allocate(st_bemv%ev11(NCELLg),st_bemv%ev12(NCELLg),st_bemv%ev13(NCELLg))
    allocate(st_bemv%ev21(NCELLg),st_bemv%ev22(NCELLg),st_bemv%ev23(NCELLg))
    allocate(st_bemv%ev31(NCELLg),st_bemv%ev32(NCELLg),st_bemv%ev33(NCELLg))
    st_bemv%xcol=xcol
    st_bemv%ycol=ycol
    st_bemv%zcol=zcol
    st_bemv%xs1=xs1
    st_bemv%xs2=xs2
    st_bemv%xs3=xs3
    st_bemv%ys1=ys1
    st_bemv%ys2=ys2
    st_bemv%ys3=ys3
    st_bemv%zs1=zs1
    st_bemv%zs2=zs2
    st_bemv%zs3=zs3
    st_bemv%ev11=ev11; st_bemv%ev12=ev12; st_bemv%ev13=ev13
    st_bemv%ev21=ev21; st_bemv%ev22=ev22; st_bemv%ev23=ev23
    st_bemv%ev31=ev31; st_bemv%ev32=ev32; st_bemv%ev33=ev33
    st_bemv%problem=problem
    st_bemv%rake=crake
  end select

  ! i=4998
  ! j=2109
  ! st_bemv%v='s'
  ! st_bemv%md='st'
  ! write(*,*) j,matel3dh_ij(i,j,st_bemv)
  ! stop
  !open(29,file='tmp')
  !do j=1,NCELLg
  !  write(*,*) j,matel3dh_ij(i,j,st_bemv)
  !end do
  !stop

  !generate kernel (H-matrix aprrox)
  if(my_rank==0) write(*,*) 'Generating kernel'
  do i=1,NCELLg
    coord(i,1)=xcol(i)
    coord(i,2)=ycol(i)
    coord(i,3)=zcol(i)
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !st_ctl%param(21)=1
  select case(problem)
  case('2dp','2dh','2dn3','3dp','3dph')
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp,st_leafmtxps,st_bemv,st_ctl,coord,eps_h)
    allocate(wws(st_leafmtxps%ndlfs))
    lrtrn=HACApK_gen_lattice_vector(st_vel,st_leafmtxps,st_ctl)
    lrtrn=HACApK_gen_lattice_vector(st_sum,st_leafmtxps,st_ctl)

    NCELL=st_vel%ndc
    allocate(y(2*NCELL),yscal(2*NCELL),dydx(2*NCELL))
    allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL),idisp(NCELL))


  case('2dn','2dnh','25d')
    st_bemv%v='xx'
    lrtrn=HACApK_generate(st_leafmtxp_xx,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_xx,st_leafmtxp_xx,st_bemv,st_ctl,coord,eps_h)
    allocate(wws(st_leafmtxp_xx%ndlfs))
    lrtrn=HACApK_gen_lattice_vector(st_vel,st_leafmtxp_xx,st_ctl)
    lrtrn=HACApK_gen_lattice_vector(st_sum,st_leafmtxp_xx,st_ctl)

    st_bemv%v='xy'
    lrtrn=HACApK_generate(st_leafmtxp_xy,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_xy,st_leafmtxp_xy,st_bemv,st_ctl,coord,eps_h)

    st_bemv%v='yy'
    lrtrn=HACApK_generate(st_leafmtxp_yy,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_yy,st_leafmtxp_yy,st_bemv,st_ctl,coord,eps_h)

    NCELL=st_vel%ndc
    allocate(y(3*NCELL),yscal(3*NCELL),dydx(3*NCELL))
    allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL),idisp(NCELL))

  case('3dnf','3dhf')
    st_bemv%md='st'
    if(slipmode.eq.'dip') st_bemv%md='dp'
    st_bemv%v='s'
    if(slipmode=='dip') st_bemv%v='d'
    rake=rake/180*pi
    st_bemv%md='o'
    st_bemv%v='o'

    lrtrn=HACApK_generate(st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_s,st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
    allocate(wws(st_leafmtxp_s%ndlfs))
    lrtrn=HACApK_gen_lattice_vector(st_vel,st_leafmtxp_s,st_ctl)
    lrtrn=HACApK_gen_lattice_vector(st_sum,st_leafmtxp_s,st_ctl)

    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_n,st_leafmtxp_n,st_bemv,st_ctl,coord,eps_h)

    NCELL=st_vel%ndc
    allocate(y(3*NCELL),yscal(3*NCELL),dydx(3*NCELL))
    allocate(phi(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL),vel(NCELL),tau(NCELL),idisp(NCELL))

  case('3dn','3dh')
    !kernel for strike slip
    st_bemv%md='st'
    st_bemv%v='s'
    lrtrn=HACApK_generate(st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_s,st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='d'
    lrtrn=HACApK_generate(st_leafmtxp_d,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_d,st_leafmtxp_d,st_bemv,st_ctl,coord,eps_h)

    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_n,st_leafmtxp_n,st_bemv,st_ctl,coord,eps_h)

    !kernel for dip slip
    st_bemv%md='dp'
    st_bemv%v='s'
    lrtrn=HACApK_generate(st_leafmtxp_s2,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_s2,st_leafmtxp_s2,st_bemv,st_ctl,coord,eps_h)

    st_bemv%v='d'
    lrtrn=HACApK_generate(st_leafmtxp_d2,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_d2,st_leafmtxp_d2,st_bemv,st_ctl,coord,eps_h)

    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n2,st_bemv,st_ctl,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_n2,st_leafmtxp_n2,st_bemv,st_ctl,coord,eps_h)


   allocate(wws(st_leafmtxp_s%ndlfs))
   lrtrn=HACApK_gen_lattice_vector(st_vel,st_leafmtxp_s,st_ctl)
   lrtrn=HACApK_gen_lattice_vector(st_sum,st_leafmtxp_s,st_ctl)

    NCELL=st_vel%ndc
   allocate(y(4*NCELL),yscal(4*NCELL),dydx(4*NCELL))
   allocate(phi(NCELL),vels(NCELL),veld(NCELL),taus(NCELL),taud(NCELL),sigma(NCELL))
   allocate(disp(NCELL),disps(NCELL),dispd(NCELL),mu(NCELL),rake(NCELL),vel(NCELL),tau(NCELL),idisp(NCELL),velp(NCELL))

  end select

  allocate(a(NCELL),b(NCELL),dc(NCELL),f0(NCELL),taudot(NCELL),tauddot(NCELL),sigdot(NCELL))

  !setting frictional parameters
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(my_rank==0) write(*,*) 'Setting fault parameters'
  !uniform
  a=a0
  b=b0
  dc=dc0
  f0=mu0
  !nonuniform parameters
  call params()
  call loading()

  !max time step
  select case(load)
  case(0)
    dtmax=0.02d0*10d0/sr
  end select

  if(foward) call foward_check()
  if(inverse) call inverse_problem()

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !setting initial condition
  !uniform
  sigma=sigma0
  tau=sigma*muinit
  mu=tau/sigma
  vel=tau/abs(tau)*velinit
  phi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
  disp=0d0
  select case(problem)
  case('3dn','3dh')
    taus=tau
    taud=0d0
  end select

  !non-uniform initial stress from subroutine initcond()
  if(project=="REGIONAL") call initcond_REGIONAL()

  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !for FDMAP-BIEM simulation
  !call input_from_FDMAP()

  !setting output files
  if(my_rank.lt.npd) then
    ! write(fname,'("output/",i0,"_",i0,".dat")') number,my_rank
    ! nout=my_rank+100
    ! open(nout,file=fname)
    write(fname,'("output/xyz",i0,"_",i0,".dat")') number,my_rank
    nout=my_rank+100
    open(nout,file=fname)
    do i=1,ncell
      i_=st_sum%lodc(i)
      write(nout,*) xcol(i_),ycol(i_),zcol(i_)
    end do
    close(nout)
    write(fname,'("output/vel",i0,"_",i0,".dat")') number,my_rank
    open(nout,file=fname,form='unformatted',access='stream')
  end if
  if(my_rank.eq.0) then
    write(fname,'("output/monitor",i0,".dat")') number
    open(52,file=fname)
    !write(fname,'("output/slip",i0,".dat")') number
    !open(46,file=fname,form='unformatted',access='stream')
    !write(fname,'("output/vel",i0,".dat")') number
    !open(47,file=fname,form='unformatted',access='stream')
    write(fname,'("output/event",i0,".dat")') number
    open(44,file=fname)
    open(19,file='job.log',position='append')
    call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    write(19,'(a20,i0,a6,a12,a6,a12)') 'Starting job number=',number,'date',sys_time(1),'time',sys_time(2)
    close(19)
    !open(73,file='output/tofd2d',access='stream')

    !if(SEAS) call open_BP(problem)
  end if

  !setting minimum time step by CFL condition
  !dtmin=0.5d0*ds/(vs*sqrt(3.d0))

  x=0.d0 !x is time
  k=0
  rk=0
  dtnxt = dtinit
  !outv=1d-6
  slipping=.false.
  eventcount=0
  sw=0
  timer=0d0
  mvelG=maxval(abs(vel))
  !output intiial condition
  if(my_rank<npd) then
    !call output_field()
    ! do i=1,ncell
    !   i_=st_sum%lodc(i)
    !   write(nout)i_
    ! end do
    write(nout) vel
  end if
  if(my_rank==0) then
    !call output_field_fd2d()
    !write(46) disp
    !write(47) vel
    !write(48) tau
    !write(49) sigma
    call output_monitor()
  end if
  !time2=MPI_Wtime()
  !output initial values


  !do i=1,NCELLg
  !  write(50,'(8e15.6,i6)') xcol(i),ycol(i),vel(i),tau(i),sigma(i),mu(i),disp(i),x,k
  !end do
  !write(50,*)
  select case(problem)
  case('2dp','2dh','2dn3','3dp','3dph')
    do i=1,NCELL
      y(2*i-1) = phi(i)
      y(2*i) = tau(i)
    end do

  case('2dn','2dnh','3dnf','3dhf','25d')
    do i=1,NCELL
      y(3*i-2) = phi(i)
      y(3*i-1) = tau(i)
      y(3*i)=sigma(i)
    end do

  case('3dn','3dh')
    do i=1,NCELL
      y(4*i-3) = phi(i)
      y(4*i-2) = taus(i)
      y(4*i-1) = taud(i)
      y(4*i)=sigma(i)
    end do

  end select
  !stop
  time2=MPI_Wtime()
  if(my_rank==0) write(*,*) 'Finished all initial processing, time(s)=',time2-time1
  time1=MPI_Wtime()
  do k=1,NSTEP1
    !parallel computing for Runge-Kutta
    dttry = dtnxt
    time3=MPI_Wtime()
    call rkqs(y,dydx,x,dttry,eps_r,dtdid,dtnxt,errmax_gb)
    time4=MPI_Wtime()
    timer=timer+time4-time3

    !limitsigma
    if(limitsigma) then
      select case(problem)
      case('2dn','3dnf','3dhf','25d')
        do i=1,NCELL
          if(y(3*i)<minsig) y(3*i)=minsig
          if(y(3*i)>maxsig) y(3*i)=maxsig
        end do
      case('3dn','3dh')
        do i=1,NCELL
          if(y(4*i)<minsig) y(4*i)=minsig
          if(y(4*i)>maxsig) y(4*i)=maxsig
        end do
      end select
    end if

    !compute physical values for control and output
    select case(problem)
    case('2dp','2dh','2dn3','3dp','3dph')
      do i = 1, NCELL
        phi(i) = y(2*i-1)
        tau(i) = y(2*i)
        disp(i)=disp(i)+vel(i)*dtdid*0.5d0 !2nd order
        vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
        disp(i)=disp(i)+vel(i)*dtdid*0.5d0
        mu(i)=tau(i)/sigma(i)
      end do

    case('2dn','2dnh','3dnf','3dhf','25d')
      do i = 1, NCELL
        phi(i) = y(3*i-2)
        tau(i) = y(3*i-1)
        sigma(i)=y(3*i)
        disp(i)=disp(i)+vel(i)*dtdid*0.5d0 !2nd order
        vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
        disp(i)=disp(i)+vel(i)*dtdid*0.5d0
        mu(i)=tau(i)/sigma(i)
      end do
    case('3dn','3dh')

    end select

    !if(outfield.and.(my_rank.lt.npd)) call output_field()

    mvel=maxval(abs(vel))
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(mvel,mvelG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(mvelG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    maxnorm=maxval(sigma)
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(maxnorm,maxnormG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(maxnormG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    minnorm=minval(sigma)
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(minnorm,minnormG,1,MPI_REAL8,MPI_MIN,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(minnormG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    meandisp=sum(disp)
    !call MPI_ALLREDUCE(meandisp,meandispG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meandisp,meandispG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meandispG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meandispG=meandispG/ncellg

    meanmu=sum(mu)
    !call MPI_ALLREDUCE(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meanmuG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meanmuG=meanmuG/ncellg

    !stop
    !output
    if(my_rank==0) call output_monitor()
      !output distribution control
    outfield=.false.
    if(mod(k,interval)==0) outfield=.true.
    if(outfield) then
      if(my_rank==0) write(*,*) 'time step=' ,k,x/365/24/60/60
      if(my_rank<npd) write(nout) vel

    end if


    !event list
    if(.not.slipping) then
      if(mvelG>1d-2) then
        slipping=.true.
        eventcount=eventcount+1
        moment0=meandispG
        !hypoloc=maxloc(abs(vel))
        onset_time=x

        !onset save
        if(slipevery.and.(my_rank==0)) then
          !write(46) disp
          !write(47) vel
          !if(project=='2DBEND')then
          !  write(48) tau-taudot*x
          !  write(49) sigma-sigdot*x
          !else
          !  write(48) tau!-taudot*x
          !  write(49) sigma!-sigdot*x
          !end if
          !call output_field()
        end if

      end if
    end if
    !
    if(slipping) then
      if(mvelG<1d-3) then
        slipping=.false.
        moment=meandispG-moment0
        !eventcount=eventcount+1
        !end of an event
        if(my_rank==0) then
          write(44,'(i0,3e15.6)') eventcount,onset_time,moment,x-onset_time
          if(slipevery) then
            !call output_field()
            !write(46) disp
            !write(47) vel
            !if(project=='2DBEND')then
            !  write(48) tau-taudot*x
            !  write(49) sigma-sigdot*x
            !else
            !  write(48) tau!-taudot*x
            !  write(49) sigma!-sigdot*x
            !end if
          end if
        end if
      end if
      !   vmaxevent=max(vmaxevent,maxval(vel))
      !   !write(53,'(i6,4e16.6)') !k,x-onset_time,sum(disp-idisp),sum(vel),sum(acg**2)
      !   !if(x-onset_time>lapse) then
      !   !  lapse=lapse+dlapse
      !   !end if
    end if

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
    !if(maxval(sigma)>=maxsig) then
    !  if(my_rank == 0) write(*,*) 'sigma exceeds maxsig'
      !exit
    !end if
  end do

  !output for FDMAP communication
  !call output_to_FDMAP()

  time2= MPI_Wtime()
  200  if(my_rank==0) then
  write(*,*) 'time(s)', time2-time1,timer
  write(*,*) 'time for matvec(s)', sum(st_ctl%time)
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
end if
!if(my_rank==0) write(19,'(a20,i0,f16.2)')'Finished job number=',number,time2-time1
Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
select case(problem)
case('2dp','2dh','2dn3','3dp','3dph')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxps)
case('2dn','2dnh','25d')
  !lrtrn=HACApK_free_leafmtxp(st_leafmtxps)
  !lrtrn=HACApK_free_leafmtxp(st_leafmtxpn)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xx)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xy)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_yy)
case('3dnf','3dhf')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_s)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_n)
case('3dn','3dh')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_s)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_d)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_n)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_s2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_d2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_n2)
end select
lrtrn=HACApK_finalize(st_ctl)
Call MPI_FINALIZE(ierr)
stop
contains
  !------------output-----------------------------------------------------------!
  subroutine output_monitor()
    implicit none
    time2=MPi_Wtime()
    write(52,'(i7,f19.4,7e16.5,f16.4)')k,x,log10(mvelG),meandispG,meanmuG,maxnormG,minnormG,errmax_gb,dtdid,time2-time1
  end subroutine
  subroutine output_field()
    implicit none
    select case(problem)
    case('3dp','3dph')
      do i=1,NCELL
        i_=st_sum%lodc(i)
        write(nout,'(6e15.6,i10)') xcol(i_),zcol(i_),log10(vel(i_)),mu(i_),disp(i_),phi(i_),k
      end do
      write(nout,*)
      write(nout,*)
    case('2dp','2dh','2dn','2dnh','2dn3','25d')
      do i=1,NCELL
        i_=st_sum%lodc(i)
        write(nout,'(i0,10e15.6,i10)') i_,xcol(i_),ycol(i_),log10(abs(vel(i))),tau(i),sigma(i),mu(i),disp(i),phi(i),x,k
      end do
      write(nout,*)
    case('3dnf','3dhf')
      do i=1,NCELL
        i_=st_sum%lodc(i)
        write(nout,'(10e14.5,i10)') xcol(i_),ycol(i_),zcol(i_),log10(vel(i)),tau(i),sigma(i),mu(i),disp(i),phi(i),x,k
      end do
      write(nout,*)
      write(nout,*)
    case('3dn','3dh')
      do i=1,NCELL
        i_=st_sum%lodc(i)
        write(nout,'(12e14.5,i10)') xcol(i_),ycol(i_),zcol(i_),log10(vel(i)),taus(i),taud(i),phi(i),mu(i),sigma(i),disps(i),dispd(i),rake(i),k
      end do
      write(nout,*)
      write(nout,*)
    end select
  end subroutine

  !------------initond-----------------------------------------------------------!
  subroutine initcond2dh(phi,sigma,tau,disp,vel)
    implicit none
    real(8),intent(out)::phi(:),sigma(:),tau(:),disp(:),vel(:)
    do i=1,NCELLg
      sigma(i)=min(17.0*xcol(i)+10d0,240d0)
      tau(i)=muinit*sigma(i)
      disp(i)=0d0
      vel(i)=velinit
      !phi(i)=a(i)*dlog(2*vref/velinit*sinh(abs(tau(i))/sigma(i)/a(i)))
      phi(i)=f0(i)+b(i)*log(b(i)*vref/vel(i))-0.001
      omega=exp((phi(i)-f0(i))/b(i))*vel(i)/vref/b(i)
      !if(my_rank==0)write(*,*) omega
    end do
  end subroutine

  subroutine initcond_REGIONALplus()
    implicit none
    real(8)::PS11,PS22,PS33,PS12,tp,tr,svalue
    !uniform tensor in a full-space
    open(97,file='psd.dat')
    do i=1,NCELL
      i_=st_sum%lodc(i)
      PS11=sigma0
      if(problem=="3dh") PS11=-zcol(i)*16.7d0+10d0
      PS22=sigma0
      PS33=sigma0
      PS12=PS22*muinit
      taus(i) = ev11(i_)*ev31(i_)*PS11 + ev12(i_)*ev32(i_)*PS22+ (ev11(i_)*ev32(i_)+ev12(i_)*ev31(i_))*PS12 + ev13(i_)*ev33(i_)*PS33
      taud(i) = ev21(i_)*ev31(i_)*PS11 + ev22(i_)*ev32(i_)*PS22+ (ev21(i_)*ev32(i_)+ev22(i_)*ev31(i_))*PS12 + ev23(i_)*ev33(i_)*PS33
      sigma(i) = ev31(i_)*ev31(i_)*PS11 + ev32(i_)*ev32(i_)*PS22+ (ev31(i_)*ev32(i_)+ev32(i_)*ev31(i_))*PS12 + ev33(i_)*ev33(i_)*PS33
      !vel(i)=velinit
      !phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(sqrt(taus(i)**2+taud(i)**2)/sigma(i)/a(i)))
      tp=sigma(i)*0.5d0
      tr=sigma(i)*0.3d0
      svalue=(tp-taus(i))/(taus(i)-tr)
      write(97,*) ycol(i),zcol(i),svalue
    end do
    close(97)
  end subroutine

  subroutine initcond_REGIONAL()
    implicit none
    real(8)::PS11,PS22,PS33,PS12
    select case(slipmode)
    case('strike')
      do i=1,NCELL
        PS11=sigma0
        if(problem=="3dhf")PS11=-zcol(i)*16.7d0+10d0
        PS22=PS11
        PS33=PS11
        PS12=-PS22*muinit
        tau(i) = ev11(i_)*ev31(i_)*PS11 + ev12(i_)*ev32(i_)*PS22+ (ev11(i_)*ev32(i_)+ev12(i_)*ev31(i_))*PS12 + ev13(i_)*ev33(i_)*PS33
        sigma(i) = ev31(i_)*ev31(i_)*PS11 + ev32(i_)*ev32(i_)*PS22+ (ev31(i_)*ev32(i_)+ev32(i_)*ev31(i_))*PS12 + ev33(i_)*ev33(i_)*PS33
        !vel(i)=velinit
        vel(i)=velinit
        phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(tau(i)/sigma(i)/a(i)))
      end do
    case('dip')
      do i=1,NCELL
        tau(i) = sigma0*muinit
        sigma(i) = sigma0
        vel(i)=velinit
        phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(tau(i)/sigma(i)/a(i)))
      end do
    end select
  end subroutine

  subroutine add_nuclei(tau,intau,inloc)
    implicit none
    real(8),intent(in)::intau
    integer,intent(in)::inloc
    real(8),intent(inout)::tau(:)
    real(8)::ra
    integer::lc
    ra=sqrt((xcol(2)-xcol(1))**2+(ycol(2)-ycol(1))**2)
    lc=int(rigid*(1.d0-pois)/pi*dc0*b0/(b0-a0)**2/sigma0/ra)
    lc=100
    !write(*,*) 'lc=',lc
    do i=1,min(nmain,ncellg)
      tau(i)=tau(i)+exp(-dble(i-inloc)**2/lc**2)*intau*tau(inloc)/abs(tau(inloc))
      !write(*,*) i,tau(i)
    end do
    return
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

  subroutine coordinate_2DBEND()
    implicit none
    character(128)::geofile2,geom
    integer::i,j,k,file_size,n,Np,Nm,ncellf,q
    real(8),allocatable::data(:),yr(:)

    geom='dbend'
    do i=1,Ncellg
      select case(geom)
        !flat fault approx
      case('bump')
        !xel(i)=5.12d0+ds0*(i-1-NCELLg/2)
        !xer(i)=5.12d0+ds0*(i-NCELLg/2)
        xel(i)=ds0*(i-1-NCELLg/2)
        xer(i)=ds0*(i-NCELLg/2)
        yel(i)=amp*exp(-(xel(i)-0.0)**2/wid**2)
        yer(i)=amp*exp(-(xer(i)-0.0)**2/wid**2)
        !write(*,*) xel(i),yel(i)
        !double bend
      case('dbend')
        xel(i)=ds0*(i-1-NCELLg/2)
        xer(i)=ds0*(i-NCELLg/2)
        yel(i)=amp*tanh((xel(i)-0d0)/wid)
        yer(i)=amp*tanh((xer(i)-0d0)/wid)
        !yel(i)=2.5*tanh((xel(i)-25d0)/5.0)-2.5*tanh((xel(i)+25d0)/5.0)
        !yer(i)=2.5*tanh((xer(i)-25d0)/5.0)-2.5*tanh((xer(i)+25d0)/5.0)
      case('sbend')
        xel(i)=ds0*(i-1-NCELLg/2)!/sqrt(1+amp**2)
        xer(i)=ds0*(i-NCELLg/2)!/sqrt(1+amp**2)
        yel(i)=amp*wid*sqrt(1.0+((xel(i)-0d0)/wid)**2)
        yer(i)=amp*wid*sqrt(1.0+((xer(i)-0d0)/wid)**2)
      end select
    end do

    !reading mesh data from mkelm.f90
    if(geofromfile) then
      geofile2='geos/'//geofile
      open(20,file=geofile2,access='stream')
      read(20) xel,xer,yel,yer
    end if


    ! geofile2='alpha0.001Lmin2N5001seed1.curve'
    ! open(32,file=geofile2,access='stream')
    ! inquire(32, size=file_size)
    ! q=file_size/8
    ! write(*,*) 'q=',q
    ! allocate(yr(q/4))
    ! allocate(data(q))
    ! read(32) data
    ! close(32)
    ! yr(1:q/4)=data(q/4+1:q/2)
    ! amp=1d-3
    ! do i=1,NCELLg
    !   !xel(i)=5.12d0+ds0*(i-1-NCELLg/2)
    !   !xer(i)=5.12d0+ds0*(i-NCELLg/2)
    !   yel(i)=yel(i)+yr(i-1)*amp
    !   yer(i)=yer(i)+yr(i)*amp
    !   yel(i)=yr(i-1)*amp
    !   yer(i)=yr(i)*amp
    ! end do

    !computing local angles and collocation points
    do i=1,NCELLg
      ds(i)=sqrt((xer(i)-xel(i))**2+(yer(i)-yel(i))**2)
      ang(i)=datan2(yer(i)-yel(i),xer(i)-xel(i))
      xcol(i)=0.5d0*(xel(i)+xer(i))
      ycol(i)=0.5d0*(yel(i)+yer(i))
      !write(*,*) xcol(i),ycol(i)
    end do

    return
  end subroutine
  subroutine coordinate2dnh()
    implicit none
    integer::i,j,k

    !flat fault with element size ds
    do i=1,NCELLg
      ds(i)=ds0
      xel(i)=(i-1)*ds0*cos(dipangle*pi/180)
      xer(i)=i*ds0*cos(dipangle*pi/180)
      yel(i)=(i-1)*ds0*sin(dipangle*pi/180)
      yer(i)=i*ds0*sin(dipangle*pi/180)
      xcol(i)=0.5d0*(xel(i)+xer(i))
      ycol(i)=0.5d0*(yel(i)+yer(i))
      ang(i)=datan2(yer(i)-yel(i),xer(i)-xel(i))
      !write(14,'(3e16.6)') xcol(i),xel(i),xer(i)
    enddo
    !close(14)
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

  subroutine coordinate3dph(imax,jmax,ds0,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
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
        zcol(k)=-(j-0.5d0)*dz-1d-9
        !xcol(k)=(i-imax/2-0.5d0)*ds0
        !zcol(k)=(j-jmax/2-0.5d0)*ds0
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
  end subroutine coordinate3dph

  subroutine coordinate3dns(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    implicit none
    integer,intent(in)::NCELLg
    real(8),intent(out)::xcol(:),ycol(:),zcol(:)
    real(8),intent(out)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    integer::i,j,k,imax,jmax
    real(8)::dipangle,xc,yc,zc,amp
    !real(4)::xl(0:2048,0:2048)

    imax=50
    jmax=50
    dipangle=30d0*pi/180d0
    do i=1,imax
      do j=1,jmax
        k=(i-1)*jmax+j
        !xcol(k)=(i-imax/2-0.5d0)*ds0
        !zcol(k)=-(j-0.5d0)*ds0-0.001d0
        xc=(i-imax/2-0.5)*ds0
        yc=-(j-0.5d0)*ds0*cos(dipangle)
        zc=-(j-0.5d0)*ds0*sin(dipangle)-1d-3!-100d0

        xs1(2*k-1)=xc-0.5d0*ds0
        xs2(2*k-1)=xc+0.5d0*ds0
        xs3(2*k-1)=xc-0.5d0*ds0
        zs1(2*k-1)=zc+0.5d0*ds0*sin(dipangle)
        zs2(2*k-1)=zc+0.5d0*ds0*sin(dipangle)
        zs3(2*k-1)=zc-0.5d0*ds0*sin(dipangle)
        ys1(2*k-1)=yc+0.5d0*ds0*cos(dipangle)
        ys2(2*k-1)=yc+0.5d0*ds0*cos(dipangle)
        ys3(2*k-1)=yc-0.5d0*ds0*cos(dipangle)

        xs2(2*k)=xc+0.5d0*ds0
        xs1(2*k)=xc+0.5d0*ds0
        xs3(2*k)=xc-0.5d0*ds0
        zs2(2*k)=zc-0.5d0*ds0*sin(dipangle)
        zs1(2*k)=zc+0.5d0*ds0*sin(dipangle)
        zs3(2*k)=zc-0.5d0*ds0*sin(dipangle)
        ys2(2*k)=yc-0.5d0*ds0*cos(dipangle)
        ys1(2*k)=yc+0.5d0*ds0*cos(dipangle)
        ys3(2*k)=yc-0.5d0*ds0*cos(dipangle)

      end do
    end do
    do k=1,ncellg
      xcol(k)=(xs1(k)+xs2(k)+xs3(k))/3.d0
      ycol(k)=(ys1(k)+ys2(k)+ys3(k))/3.d0
      zcol(k)=(zs1(k)+zs2(k)+zs3(k))/3.d0
      write(*,*) xcol(k),ycol(k),zcol(k)
    end do

    ! open(30,file='roughsurf.txt')
    ! do k=0,2048
    !   read(30,*) xl(k,0:2048)
    ! end do
    ! close(30)
    ! amp=0.000d0
    ! if(my_rank==0) open(32,file='tmp')
    ! do i=1,NCELLg
    !   xcol(i)=(xs1(i)+xs2(i)+xs3(i))/3.d0
    !   zcol(i)=(zs1(i)+zs2(i)+zs3(i))/3.d0
    !
    !   j=int((xs1(i)+10)*102.4)
    !   k=int(-102.4*zs1(i))
    !   ys1(i)=xl(j,k)*amp
    !   j=int((xs2(i)+10)*102.4)
    !   k=int(-102.4*zs2(i))
    !   ys2(i)=xl(j,k)*amp
    !   j=int((xs3(i)+10)*102.4)
    !   k=int(-102.4*zs3(i))
    !   ys3(i)=xl(j,k)*amp
    !   ycol(i)=(ys1(i)+ys2(i)+ys3(i))/3.d0
    !   if(my_rank==0) write(32,*) xcol(i),ycol(i),zcol(i)
    ! end do

    return
  end subroutine coordinate3dns

  subroutine coordinate3dns2(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    implicit none
    integer,intent(in)::NCELLg
    real(8),intent(out)::xcol(:),ycol(:),zcol(:)
    real(8),intent(out)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    integer::i,j,k
    real(8)::dipangle,xc,yc,zc

    !imax=150
    !jmax=150
    amp=0.1
    wid=0.1
write(*,*) imax,jmax
    do i=1,imax
      do j=1,jmax
        k=(i-1)*jmax+j
        !xcol(k)=(i-imax/2-0.5d0)*ds0
        !zcol(k)=-(j-0.5d0)*ds0-0.001d0
        xc=(i-imax/2-0.5)*ds0
        zc=-(j-0.5d0)*ds0-1d-9!-100d0
        !yc=0d0

        xs1(2*k-1)=xc-0.5d0*ds0
        xs2(2*k-1)=xc+0.5d0*ds0
        xs3(2*k-1)=xc-0.5d0*ds0
        zs1(2*k-1)=zc+0.5d0*ds0
        zs2(2*k-1)=zc+0.5d0*ds0
        zs3(2*k-1)=zc-0.5d0*ds0
        ys1(2*k-1)=sbend(xs1(2*k-1),amp,wid)
        ys2(2*k-1)=sbend(xs2(2*k-1),amp,wid)
        ys3(2*k-1)=sbend(xs3(2*k-1),amp,wid)

        xs2(2*k)=xc+0.5d0*ds0
        xs1(2*k)=xc+0.5d0*ds0
        xs3(2*k)=xc-0.5d0*ds0
        zs2(2*k)=zc-0.5d0*ds0
        zs1(2*k)=zc+0.5d0*ds0
        zs3(2*k)=zc-0.5d0*ds0
        ys1(2*k)=sbend(xs1(2*k),amp,wid)
        ys2(2*k)=sbend(xs2(2*k),amp,wid)
        ys3(2*k)=sbend(xs3(2*k),amp,wid)

      end do
    end do

    do k=1,ncellg

      xcol(k)=(xs1(k)+xs2(k)+xs3(k))/3.d0
      ycol(k)=(ys1(k)+ys2(k)+ys3(k))/3.d0
      zcol(k)=(zs1(k)+zs2(k)+zs3(k))/3.d0
      write(*,*) xcol(k),ycol(k),zcol(k)
    end do
    return
  end subroutine coordinate3dns2

  ! subroutine coordinate3dn_ROUGH(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
  !   implicit none
  !   integer,intent(in)::NCELLg
  !   real(8),intent(inout)::xcol(:),ycol(:),zcol(:)
  !   real(8),intent(inout)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
  !   real(4)::xl(0:2048,0:2048)
  !   !real(8),parameter::amp=1d-4
  !   real(8)::area
  !   integer::i,j,k
  !
  !     open(30,file='roughsurf.txt')
  !     do k=0,2048
  !       read(30,*) xl(k,0:2048)
  !     end do
  !     close(30)
  !     !if(my_rank==0) open(32,file='tmp')
  !     do i=1,NCELLg
  !       j=int((ys1(i)+10)*102.4)
  !       k=int(-102.4*zs1(i))
  !       xs1(i)=xl(j,k)*amp
  !       j=int((ys2(i)+10)*102.4)
  !       k=int(-102.4*zs2(i))
  !       xs2(i)=xl(j,k)*amp
  !       j=int((ys3(i)+10)*102.4)
  !       k=int(-102.4*zs3(i))
  !       xs3(i)=xl(j,k)*amp
  !       xcol(i)=(xs1(i)+xs2(i)+xs3(i))/3.d0
  !       !if(my_rank==0) write(32,*) xcol(i),ycol(i),zcol(i)
  !     end do
  !   return
  ! end subroutine

  subroutine coordinate3dn_SEAMOUNT(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    implicit none
    integer,intent(in)::NCELLg
    real(8),intent(inout)::xcol(:),ycol(:),zcol(:)
    real(8),intent(inout)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    real(8)::area,amp,wid
    integer::i,j,k
    !amp=0.0d0
    wid=0.5d0
    do i=1,ncellg
      zs1(i)=zs1(i)+bump(xs1(i),ys1(i),amp,wid)
      zs2(i)=zs2(i)+bump(xs2(i),ys2(i),amp,wid)
      zs3(i)=zs3(i)+bump(xs3(i),ys3(i),amp,wid)
      zcol(i)=(zs1(i)+zs2(i)+zs3(i))/3d0
    end do
    return
  end subroutine

  function bump(x,y,amp,wid)
    implicit none
    real(8)::x,y,amp,wid,bump,rr
    rr=(x-25d0)**2+(y+5d0)**2
    bump=amp*exp(-rr/wid**2)
    return
  end function

  function dbend(y,amp,wid)
    implicit none
    real(8)::y,amp,wid,dbend
    dbend=amp*tanh((y-0d0)/wid)
  end function

  function sbend(y,amp,wid)
    implicit none
    real(8)::y,amp,wid,sbend
    !sbend=amp*tanh((y-0d0)/wid)
    sbend=amp*wid*sqrt(1.0+((y-0d0)/wid)**2)
  end function

  subroutine params()
    implicit none
    real(8)::len,cent,dep,a_max,xd,r
    integer::i

    !uniform
    if(project=='SEAMOUNT')then
    a_max=0.030d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      dep=sqrt(ycol(i_)**2+zcol(i_)**2)
      r=max(abs(dep-10d0)-7d0,abs(xcol(i_)-22d0)-24d0)/3d0
      a(i)=a0+r*(a_max-a0)
      a(i)=max(a(i),a0)
      b(i)=b0
      dc(i)=dc0
      f0(i)=mu0
      !if(aftershock.and.i<=nmain) f0(i)=mu0-0.05d0
      !dc is proportional to fault size
      !if(dcscale) dc(i)=ds(i)/0.004d0*0.001d0
      !if(aftershock.and.(i>nmain)) dc(i)=0.001d0
      !if(creep.and.abs(i-ncellg/2)<100) a(i)=0.030d0
    end do
    end if

    !if(project=='SEAMOUNT')then
    !  do i=1,ncellg
    !    if(abs(xcol(i))>40d0.or.abs(zcol(i)+0.8)>0.7d0) a(i)=0.024
    !  end do
    !end if
  return
  end subroutine

  subroutine loading()
    implicit none
    !character(128),intent(in)::problem
    !integer,intent(in)::NCELLg
    !real(8),intent(in)::sr
    !real(8),intent(out)::taudot(:),tauddot(:),sigdot(:)
    real(8)::factor,edge,ret1,ret2,xx1,xx2,xy1,xy2,yy1,yy2,lang
    integer::i
    character(128)::v
    open(15,file='sr')
    select case(problem)
    case('2dp','3dp','2dh','3dhf','3dph')
      taudot=sr
      tauddot=0d0
      sigdot=0d0
    case('3dnf')
      do i=1,NCELL
        !taudot(i) = -(ev11(i)*ev32(i)+ev12(i)*ev31(i))*sr
        !sigdot(i) = -(ev31(i)*ev32(i)+ev32(i)*ev31(i))*sr
        taudot(i)=sr
        sigdot(i)=0d0
      end do

    case('3dn','3dh')
      taudot=0d0
      tauddot=0d0
      sigdot=0d0
      do i=1,ncell
        taudot(i) = -(ev11(i)*ev32(i)+ev12(i)*ev31(i))*sr
        sigdot(i) = -(ev31(i)*ev32(i)+ev32(i)*ev31(i))*sr
      end do
    end select
    close(15)
  end subroutine
  subroutine evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33)
    !calculate ev for each element
    implicit none
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    real(8),intent(out)::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:)
    real(8)::rr,vba(0:2),vca(0:2)

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
      !if(my_rank==0) write(*,*)ev21(k),ev22(k),ev23(k)
    end do

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
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),phitmp(NCELL)
    real(8) :: dtaudt(NCELL),dsigdt(NCELL),dphidt(NCELL)
    real(8) :: taustmp(NCELL),taudtmp(NCELL),velstmp(NCELL),veldtmp(NCELL),dtausdt(NCELL),dtauddt(NCELL)
    real(8) :: sum_gs(NCELL),sum_gn(NCELL),sum_gd(NCELL)!,velstmpG(NCELLg),veldtmpG(NCELLg)
    real(8) :: sum_xx(NCELL),sum_xy(NCELL),sum_yy(NCELL)!,sum_xz(NCELL),sum_yz(NCELL),sum_zz(NCELL)
    !real(8) :: sum_xxG(NCELLg),sum_xyG(NCELLg),sum_yyG(NCELLg)!,sum_xzG(NCELLg),sum_yzG(NCELLg),sum_zzG(NCELLg)
    !real(8) :: sum_xx2G(NCELLg),sum_xy2G(NCELLg),sum_yy2G(NCELLg),sum_xz2G(NCELLg),sum_yz2G(NCELLg),sum_zz2G(NCELLg)
    !real(8) :: veltmpG(NCELLg),sum_gsg(NCELLg),sum_gng(NCELLg),sum_gdg(NCELLg)!,efftmpG(NCELLg)
    !real(8) :: sum_gs2G(NCELLg),sum_gd2G(NCELLg),sum_gn2G(NCELLg)
    real(8) :: time3,time4,c1, c2, c3, arg,arg2,c,g,tauss,Arot(3,3),p(6),fac,sxx0,sxy0,syy0
    integer :: i, j, nc,ierr,lrtrn,i_

    !if(my_rank==0) then
    select case(problem)
    case('2dp','2dn3','3dp','2dh','3dph')
      do i = 1, NCELL
        phitmp(i) = y(2*i-1)
        tautmp(i) = y(2*i)
        sigmatmp(i)=sigma(i) !normal stress is constant for planar fault
        veltmp(i) = 2*vref*dexp(-phitmp(i)/a(i))*dsinh(tautmp(i)/sigmatmp(i)/a(i))
        !write(*,*) veltmp(i)
      enddo

      if(load==1) then
        st_vel%vs=veltmp-vpl
      else
        st_vel%vs=veltmp
      end if
      call HACApK_adot_lattice_hyp(st_sum,st_LHp,st_ctl,wws,st_vel)
      sum_gs(:)=st_sum%vs(:)

      !if(my_rank.eq.0) write(*,*) time4-time3
      do i=1,NCELL
        sum_gn(i)=0.d0
        sum_gs(i)=sum_gs(i)+taudot(i)
      end do

      call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)
      !call deriv_c(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)

      do i = 1, NCELL
        dydx(2*i-1) = dphidt(i)
        dydx(2*i) = dtaudt(i)
      enddo
      !call MPI_SCATTERv(sum_gsG,NCELL,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    case('2dn','2dnh','25d')

    case('2dn_vector','3dnf','3dhf')
      do i = 1, NCELL
        phitmp(i) = y(3*i-2)
        tautmp(i) = y(3*i-1)
        sigmatmp(i) = y(3*i)
        veltmp(i) = 2*vref*dexp(-phitmp(i)/a(i))*dsinh(tautmp(i)/sigmatmp(i)/a(i))
      enddo
      !matrix-vector mutiplation

      if(load==1) then
        st_vel%vs=veltmp-vpl
      else
        st_vel%vs=veltmp
      end if
      call HACApK_adot_lattice_hyp(st_sum,st_LHp_s,st_ctl,wws,st_vel)
      sum_gs(:)=st_sum%vs(:)
      call HACAPK_adot_lattice_hyp(st_sum,st_LHP_n,st_ctl,wws,st_vel)
      sum_gn(:)=st_sum%vs(:)
      !sum_gn=0d0

      do i=1,NCELL
        sum_gs(i)=sum_gs(i)+taudot(i)
        sum_gn(i)=sum_gn(i)+sigdot(i)
      end do
      !write(*,*)sum_gs(1)
      call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)

      do i = 1, NCELL
        dydx(3*i-2) = dphidt(i)
        dydx(3*i-1) = dtaudt(i)
        dydx(3*i) = dsigdt(i)
      enddo

    case('3dn','3dh')

    end select

    return
  end subroutine

  ! subroutine deriv_a(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
  !   implicit none
  !   integer::i
  !   real(8)::arg
  !   !type(t_deriv),intent(in) ::
  !   real(8),intent(in)::sum_gs(:),sum_gn(:),veltmp(:),tautmp(:),sigmatmp(:)
  !   !real(8),intent(in)::a(:),b(:),dc(:),vc(:)
  !   !real(8),intent(in)::mu0,vref,vs,rigid,alpha
  !   real(8),intent(out)::dlnvdt(:),dtaudt(:),dsigdt(:)
  !   do i=1,size(sum_gs)
  !     dsigdt(i)=sum_gn(i)
  !     arg=dc(i)/vref*(exp((tautmp(i)/sigmatmp(i)-mu0-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
  !     dlnvdt(i)=sum_gs(i)-b(i)*sigmatmp(i)/(arg+dc(i)/vc(i))*(1.d0-veltmp(i)*arg/dc(i))+(tautmp(i)/sigmatmp(i)-alpha)*dsigdt(i)
  !     dlnvdt(i)=dlnvdt(i)/(a(i)*sigmatmp(i)+0.5d0*rigid*veltmp(i)/vs)
  !     dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
  !   end do
  ! end subroutine
  ! subroutine deriv_s(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
  !   implicit none
  !   integer::i
  !   real(8)::arg
  !   !type(t_deriv),intent(in) ::
  !   real(8),intent(in)::sum_gs(:),sum_gn(:),veltmp(:),tautmp(:),sigmatmp(:)
  !   !real(8),intent(in)::a(:),b(:),dc(:),vc(:)
  !   !real(8),intent(in)::mu0,vref,vs,rigid,alpha
  !   real(8),intent(out)::dlnvdt(:),dtaudt(:),dsigdt(:)
  !   do i=1,size(sum_gs)
  !     dsigdt(i)=sum_gn(i)
  !     arg=dc(i)/vref*(exp((tautmp(i)/sigmatmp(i)-f0(i)-a(i)*dlog(veltmp(i)/vref))/b(i))-vref/vc(i))
  !     dlnvdt(i)=sum_gs(i)+b(i)*sigmatmp(i)*veltmp(i)*dlog(veltmp(i)*arg/dc(i))+(tautmp(i)/sigmatmp(i)-alpha)*dsigdt(i)
  !     dlnvdt(i)=dlnvdt(i)/(a(i)*sigmatmp(i)+0.5d0*rigid*veltmp(i)/vs)
  !     dtaudt(i)=sum_gs(i)-0.5d0*rigid*veltmp(i)/vs*dlnvdt(i)
  !   end do
  ! end subroutine
  subroutine deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)
    implicit none
    integer::i,i_
    real(8)::fss,dvdtau,dvdsig,dvdphi
    !real(8),parameter::fw=0.2
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs(:),sum_gn(:),phitmp(:),tautmp(:),sigmatmp(:),veltmp(:)
    real(8),intent(out)::dphidt(:),dtaudt(:),dsigdt(:)
    do i=1,ncell
      dsigdt(i)=sum_gn(i)
      !write(*,*) 'dsigdt',dsigdt(i)

      !regularized slip law
      !fss=mu0+(a(i_)-b(i_))*dlog(abs(veltmp(i))/vref)
      !fss=fw(i_)+(fss-fw(i_))/(1.d0+(veltmp(i)/vw(i_))**8)**0.125d0 !flash heating
      !dphidt(i)=-abs(veltmp(i))/dc(i_)*(abs(tautmp(i))/sigmatmp(i)-fss)

      !regularized aing law
      dphidt(i)=b(i)/dc(i)*vref*dexp((f0(i)-phitmp(i))/b(i))-b(i)*abs(veltmp(i))/dc(i)

      !regularized aging law with cutoff velocity for evolution
      !dphidt(i)=b(i_)/dc(i_)*vref*dexp((f0(i_)-phitmp(i))/b(i_))*(1d0-abs(veltmp(i))/vref*(exp((phitmp(i)-f0(i_))/b(i_))-vref/vc(i_)))


      dvdtau=2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmatmp(i)/a(i))/(a(i)*sigmatmp(i))
      dvdsig=-2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmatmp(i)/a(i))*tautmp(i)/(a(i)*sigmatmp(i)**2)
      !dvdphi=2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
      dvdphi=-veltmp(i)/a(i)
      !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*phitmp(i)*dvdsig*sigmatmp(i))
      dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))
      dtaudt(i)=dtaudt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      !write(*,*) rigid/vs*dvdtau
      if(veltmp(i)<=0d0) then
        dvdtau=2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmatmp(i)/a(i))/(a(i)*sigmatmp(i))
        dvdsig=-2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmatmp(i)/a(i))*tautmp(i)/(a(i)*sigmatmp(i)**2)
        !sign ok?
        !dvdphi=2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
        dvdphi=-veltmp(i)/a(i)
        !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(-dvdphi*phitmp(i)*dvdsig*sigmatmp(i))
        dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))
        dtaudt(i)=dtaudt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      end if
    end do
  end subroutine
  subroutine deriv_3dn(sum_gs,sum_gd,sum_gn,phitmp,taustmp,taudtmp,tautmp,sigmatmp,veltmp,dphidt,dtausdt,dtauddt,dsigdt)
    implicit none
    integer::i
    real(8)::fss,dvdtau,dvdsig,dvdphi,absV
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs(:),sum_gd(:),sum_gn(:),phitmp(:),taustmp(:),taudtmp(:),tautmp(:),sigmatmp(:),veltmp(:)
    real(8),intent(out)::dphidt(:),dtausdt(:),dtauddt(:),dsigdt(:)
    do i=1,ncell
      !write(*,*) 'vel',veltmp(i)
      dsigdt(i)=sum_gn(i)
      !fss=mu0+(a(i)-b(i))*dlog(abs(veltmp(i))/vref)
      !fss=fw(i)+(fss-fw(i))/(1.d0+(veltmp(i)/vw(i))**8)**0.125d0 !flash heating
      !slip law
      !dphidt(i)=-abs(veltmp(i))/dc(i)*(abs(tautmp(i))/sigmatmp(i)-fss)
      !aing law
      dphidt(i)=b(i)*vref/dc(i)*exp((f0(i)-phitmp(i))/b(i))-b(i)*veltmp(i)/dc(i)
      dvdtau=2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmatmp(i)/a(i))/(a(i)*sigmatmp(i))
      dvdsig=-2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmatmp(i)/a(i))*tautmp(i)/(a(i)*sigmatmp(i)**2)
      dvdphi=-veltmp(i)/a(i)
      dtausdt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))*(taustmp(i)/tautmp(i))
      dtausdt(i)=dtausdt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      dtauddt(i)=sum_gd(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))*(taudtmp(i)/tautmp(i))
      dtauddt(i)=dtauddt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
    end do
  end subroutine

  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,hdid,hnext,errmax_gb)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    !integer::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::htry,eps
    real(8),intent(inout)::y(:),x,dydx(:)
    real(8),intent(out)::hdid,hnext,errmax_gb !hdid: resulatant dt hnext: htry for the next
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    integer :: i,ierr,loc
    real(8) :: errmax,h,xnew,htemp,dtmin
    real(8),dimension(size(y))::yerr,ytemp,yscal
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4

    h=htry
    !dtmin=0.5d0*minval(ds)/vs
    !call derivs(x,y,dydx)
    do while(.true.)

      call rkck(y,x,h,ytemp,yerr)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
      !if(ierr==1) then

      !end if

      errmax=0d0
      select case(problem)
      case('3dn','3dh')
        do i=1,NCELL
          if(abs(yerr(4*i-3)/ytemp(4*i-3))/eps>errmax) errmax=abs(yerr(4*i-3)/ytemp(4*i-3))/eps
          !errmax=errmax+yerr(4*i-3)**2
        end do
      case('3dnf','3dhf','2dn','2dnh','25d')
        do i=1,NCELL
          if(abs(yerr(3*i-2)/ytemp(3*i-2))/eps>errmax) errmax=abs(yerr(3*i-2)/ytemp(3*i-2))/eps
          !errmax=errmax+yerr(3*i-2)**2
        end do
      case('2dh','2dp','2dn3','3dph','3dp')
        do i=1,NCELL
          if(abs(yerr(2*i-1))/ytemp(2*i-1)/eps>errmax) errmax=abs(yerr(2*i-1))/ytemp(2*i-1)/eps
          !errmax=errmax+yerr(2*i-1)**2
        end do
      end select
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      !call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      !errmax_gb=sqrt(errmax_gb/NCELLg)/eps

      call MPI_reduce(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      call MPI_bcast(errmax_gb,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
      !write(*,*) h,errmax,errmax_gb


      ! call MPI_reduce(errmax,errmax_gb,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      ! errmax_gb=sqrt(errmax_gb/NCELLg)/eps
      ! call MPI_bcast(errmax_gb,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)



      !write(*,*) h,errmax_gb
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
        write(*,*)'dt is too small'
        stop
      end if

    end do

    hnext=min(2*h,SAFETY*h*(errmax_gb**PGROW))
    if(load==0)hnext=min(hnext,dtmax)
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
    include 'mpif.h'
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

    !    -- 2nd step --
    call derivs(x+a2*h, ytemp, ak2)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B31*ak1(i)+B32*ak2(i))
    end do

    !     -- 3rd step --
    call derivs(x+a3*h, ytemp, ak3)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B41*ak1(i)+B42*ak2(i)+B43*ak3(i))
    end do

    !     -- 4th step --
    call derivs(x+a4*h, ytemp, ak4)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B51*ak1(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
    end do

    !     -- 5th step --
    call derivs(x+a5*h, ytemp, ak5)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+h*(B61*ak1(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
    end do

    !     -- 6th step --
    call derivs(x+a6*h, ytemp, ak6)!,,st_leafmtxp,st_bemv,st_ctl)
    !$omp parallel do
    do i=1,size(y)
      yout(i)=y(i)+h*(C1*ak1(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
    end do


    !$omp parallel do
    do i=1,size(y)
      yerr(i)=h*(DC1*ak1(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
      !if(abs(yerr(i))>=1d6)ierr=1
    end do
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
    !vel(15)=1d0
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
      if(my_rank==0) then
        do i=1,NCELLg
          taudot(i)=0.5d0*(a(i)-dc(i))*dsin(-2*ang(i))+b(i)*dcos(-2*ang(i))
          sigdot(i)=-(0.5d0*(a(i)+dc(i))-0.5d0*(a(i)-dc(i))*dcos(2*ang(i))-b(i)*dsin(2*ang(i)))
          write(29,'(4e16.4)') xcol(i),ang(i),taudot(i),sigdot(i)
        end do
      end if
    case('3dp','3dph')
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,a,vel)

      if(my_rank==0) then
        do i=1,NCELLg
          write(29,'(3e16.4)') xcol(i),zcol(i),a(i)
        end do
      end if
    case('3dn','3dh')
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_s2,st_bemv,st_ctl,a,vel)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_d2,st_bemv,st_ctl,b,vel)
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n2,st_bemv,st_ctl,dc,vel)
      if(my_rank==0) then
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
      if(my_rank==0) then
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
      if(my_rank==0) then
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
      if(my_rank==0) then
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
      if(my_rank==0) then
        do i=1,NCELLg
          write(29,'(4e16.4)') xcol(i),ycol(i),zcol(i),sigdot(i)
        end do
      end if
    end select
    Call MPI_FINALIZE(ierr)
    stop
  end subroutine

  function rtnewt(prev,eps,nst,p,t0,sum)
    integer::j
    integer,parameter::jmax=20
    real(8)::rtnewt,prev,eps
    real(8)::f,df,dx,sum,nst,p,t0
    rtnewt=prev
    !write(*,*) rtnewt
    do j=1,jmax
      x=rtnewt
      f=x+ieta*sigma0*(mu0+(a0-b0)*log(x/vref))-vpl
      df=1+ieta*sigma0*(a0-b0)/x
      dx=f/df
      rtnewt=rtnewt-dx
      !write(*,*) rtnewt
      if(abs(dx)<eps) return
    end do
    write(*,*) 'maximum iteration'
    stop
  end function
end program

subroutine open_bp(problem)
  character(128),intent(in)::problem
  real(8)::xd(81)
  select case(problem)
  !SEAS BP5
  case('3dph')
    open(101,file="output/fltst_strk-36dp+00")
    open(102,file="output/fltst_strk-16dp+00")
    open(103,file="output/fltst_strk+00dp+00")
    open(104,file="output/fltst_strk+16dp+00")
    open(105,file="output/fltst_strk+36dp+00")
    open(106,file="output/fltst_strk-24dp+10")
    open(107,file="output/fltst_strk-16dp+10")
    open(108,file="output/fltst_strk+00dp+10")
    open(109,file="output/fltst_strk+16dp+10")
    open(110,file="output/fltst_strk+00dp+22")
    do i=101,110
      write(i,*)"# This is the header:"
      write(i,*)"# problem=SEAS Benchmark BP5-QD"
      write(i,*)"# code=hbi"
      write(i,*)"# modeler=So Ozawa"
      write(i,*)"# date=2021/03/19"
      write(i,*)"# element_size=500m"
      write(i,*)"# Column #1 = Time (s)"
      write(i,*)"# Column #2 = Slip_2(m)"
      write(i,*)"# Column #3 = Slip_3(m)"
      write(i,*)"# Column #4 = Slip_rate_2(log10 m/s)"
      write(i,*)"# Column #5 = Slip_rate_3(log10 m/s)"
      write(i,*)"# Column #6 = Shear_stress_2 (MPa)"
      write(i,*)"# Column #7 = Shear_stress_3 (MPa)"
      write(i,*)"# Column #8 = State (log10 s)"
      write(i,*)"# The line below lists the names of the data fields"
      write(i,*)"t slip_2 slip_3 slip_rate_2 slip_rate_3 shear_stress_2 shear_stress_3 state"
      write(i,*)"# Here is the time-series data."
    end do

    open(120,file="output/global.dat")
    i=120
    write(i,*)"# This is the file header:"
    write(i,*)"# problem=SEAS Benchmark BP4-QD"
    write(i,*)"# code=hbi"
    write(i,*)"# modeler=So Ozawa"
    write(i,*)"# date=2021/03/19"
    write(i,*)"# element_size=500m"
    write(i,*)"# Column #1 = Time (s)"
    write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
    write(i,*)"# Column #3 = Moment rate (N-m/s)"
    write(i,*)"# The line below lists the names of the data fields"
    write(i,*)"t max_slip_rate moment_rate"
    write(i,*)"# Here is the time-series data."

    open(130,file="output/rupture.dat")
    i=130
    write(i,*)"# This is the file header:"
    write(i,*)"# problem=SEAS Benchmark BP4-QD"
    write(i,*)"# code=hbi"
    write(i,*)"# modeler=So Ozawa"
    write(i,*)"# date=2021/03/19"
    write(i,*)"# element_size=500m"
    write(i,*)"# Column #1 = x2 (m)"
    write(i,*)"# Column #2 = x3 (m)"
    write(i,*)"# Column #3 = time (s)"
    write(i,*)"# The line below lists the names of the data fields"
    write(i,*)"x2 x3 t"
    write(i,*)"# Here is the data."

  !SEAS BP4
    case('3dp')
    open(101,file="output/fltst_strk-360dp+000")
    open(102,file="output/fltst_strk-225dp-750")
    open(103,file="output/fltst_strk-165dp-120")
    open(104,file="output/fltst_strk-165dp+000")
    open(105,file="output/fltst_strk-165dp+120")
    open(106,file="output/fltst_strk+000dp-210")
    open(107,file="output/fltst_strk+000dp-120")
    open(108,file="output/fltst_strk+000dp+000")
    open(109,file="output/fltst_strk+000dp+120")
    open(110,file="output/fltst_strk+000dp+210")
    open(111,file="output/fltst_strk+165dp-120")
    open(112,file="output/fltst_strk+165dp+000")
    open(113,file="output/fltst_strk+165dp+120")
    open(114,file="output/fltst_strk+360dp+000")
    do i=101,114
      write(i,*)"# This is the header:"
      write(i,*)"# problem=SEAS Benchmark BP4-QD"
      write(i,*)"# code=hbi"
      write(i,*)"# modeler=So Ozawa"
      write(i,*)"# date=2021/03/19"
      write(i,*)"# element_size=500m"
      write(i,*)"# Column #1 = Time (s)"
      write(i,*)"# Column #2 = Slip_2(m)"
      write(i,*)"# Column #3 = Slip_3(m)"
      write(i,*)"# Column #4 = Slip_rate_2(log10 m/s)"
      write(i,*)"# Column #5 = Slip_rate_3(log10 m/s)"
      write(i,*)"# Column #6 = Shear_stress_2 (MPa)"
      write(i,*)"# Column #7 = Shear_stress_3 (MPa)"
      write(i,*)"# Column #8 = State (log10 s)"
      write(i,*)"# The line below lists the names of the data fields"
      write(i,*)"t slip_2 slip_3 slip_rate_2 slip_rate_3 shear_stress_2 shear_stress_3 state"
      write(i,*)"# Here is the time-series data."
    end do

    open(120,file="output/global.dat")
    i=120
    write(i,*)"# This is the file header:"
    write(i,*)"# problem=SEAS Benchmark BP4-QD"
    write(i,*)"# code=hbi"
    write(i,*)"# modeler=So Ozawa"
    write(i,*)"# date=2021/03/19"
    write(i,*)"# element_size=500m"
    write(i,*)"# Column #1 = Time (s)"
    write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
    write(i,*)"# Column #3 = Moment rate (N-m/s)"
    write(i,*)"# The line below lists the names of the data fields"
    write(i,*)"t max_slip_rate moment_rate"
    write(i,*)"# Here is the time-series data."

    open(130,file="output/rupture.dat")
    i=130
    write(i,*)"# This is the file header:"
    write(i,*)"# problem=SEAS Benchmark BP4-QD"
    write(i,*)"# code=hbi"
    write(i,*)"# modeler=So Ozawa"
    write(i,*)"# date=2021/03/19"
    write(i,*)"# element_size=500m"
    write(i,*)"# Column #1 = x2 (m)"
    write(i,*)"# Column #2 = x3 (m)"
    write(i,*)"# Column #3 = time (s)"
    write(i,*)"# The line below lists the names of the data fields"
    write(i,*)"x2 x3 t"
    write(i,*)"# Here is the data."

    !SEAS BP3
    case('2dnh')
    open(101,file="output/fltst_dp000")
    open(102,file="output/fltst_dp025",status='replace')
    open(103,file="output/fltst_dp050",status='replace')
    open(104,file="output/fltst_dp075",status='replace')
    open(105,file="output/fltst_dp100",status='replace')
    open(106,file="output/fltst_dp125",status='replace')
    open(107,file="output/fltst_dp150",status='replace')
    open(108,file="output/fltst_dp175",status='replace')
    open(109,file="output/fltst_dp200",status='replace')
    open(110,file="output/fltst_dp250",status='replace')
    open(111,file="output/fltst_dp300",status='replace')
    open(112,file="output/fltst_dp350",status='replace')
    do i=101,112
      write(i,*)"# This is the header:"
      write(i,*)"# problem=SEAS Benchmark BP3-QD"
      write(i,*)"# code=hbi"
      write(i,*)"# modeler=So Ozawa"
      write(i,*)"# date=2021/01/22"
      write(i,*)"# element_size=25m"
      write(i,*)"# location= on fault, 0km down-dip distance"
      write(i,*)"# Column #1 = Time (s)"
      write(i,*)"# Column #2 = Slip (m)"
      write(i,*)"# Column #3 = Slip rate (log10 m/s)"
      write(i,*)"# Column #4 = Shear stress (MPa)"
      write(i,*)"# Column #5 = Normal stress (MPa)"
      write(i,*)"# Column #6 = State (log10 s)"
      write(i,*)"# The line below lists the names of the data fields"
      write(i,*)"t slip slip_rate shear_stress normal_stress state"
      write(i,*)"# Here is the time-series data."
    end do
    open(121,file="output/slip.dat",status='replace')
    open(122,file="output/shear_stress.dat",status='replace')
    open(123,file="output/normal_stress.dat",status='replace')

    do i=121,123
    write(i,*)"# This is the file header:"
    write(i,*)"# problem=SEAS Benchmark BP3-QD"
    write(i,*)"# code=hbi"
    write(i,*)"# modeler=So Ozawa"
    write(i,*)"# date=2021/03/16"
    write(i,*)"# element_size=25m"
    write(i,*)"# Column #1 = Time (s)"
    write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
    end do

    write(121,*)"# Column #3-83 = Slip (m)"
    write(122,*)"# Column #3-83 = Shear stress (MPa)"
    write(123,*)"# Column #3-83 = Normal stress (MPa)"

    do i=121,123
    write(i,*)"# The line below lists the names of the data fields"
    write(i,*)"xd"
    end do
    write(121,*)"t max_slip_rate slip"
    write(122,*)"t max_slip_rate shear_stress"
    write(123,*)"t max_slip_rate normal_stress"
    do i=121,123
    write(i,*)"# Here are the data."
    end do
    do i=1,81
      xd(i)=(i-1)*500d0
    end do
    write(121,'(83e22.14)') 0d0,0d0,xd
    write(122,'(83e22.14)') 0d0,0d0,xd
    write(123,'(83e22.14)') 0d0,0d0,xd
  end select
return
end subroutine
subroutine debug()

end subroutine
