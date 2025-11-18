program main
  !$ use omp_lib
  use m_HACApK_solve
  use m_HACApK_base
  use m_HACApK_use
  use mod_constant
  use m_HACApK_calc_entry_ij
  use mod_diffusion
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
  type(st_HACApK_lcontrol) :: st_ctl,st_ctl2,st_ctl3,st_ctl4
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

  !diffusion
  type(t_params) :: param_diff

  !for MPI communication and time
  integer::counts2,icomm,np,npd,ierr,my_rank,npgl,noutloc,locid(10),rankloc(10),loc_(10),rank_mvel,mvel_loc(1),mvel_loc_
  integer,allocatable::displs(:),rcounts(:),vars(:)
  integer:: date_time(8)
  character(len=10):: sys_time(3)
  real(8)::time1,time2,time3,time4,timer,timeH

  !for fault geometry
  ! real(8),allocatable::xcol(:),ycol(:),zcol(:),ds(:),dsl(:),dsd(:)
  ! real(8),allocatable::xs1(:),xs2(:),xs3(:),xs4(:) !for 3dp
  ! real(8),allocatable::zs1(:),zs2(:),zs3(:),zs4(:) !for 3dp
  ! real(8),allocatable::ys1(:),ys2(:),ys3(:),ys4(:) !for 3dn
  ! real(8),allocatable::xel(:),xer(:),yel(:),yer(:),ang(:),angd(:)
  ! real(8),allocatable::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:)

  !parameters for each elements
  real(8),allocatable::a(:),b(:),dc(:),f0(:),fw(:),vw(:),vc(:),taudot(:),tauddot(:),sigdot(:)
  real(8),allocatable::vpl(:),values(:,:),etav(:),etab(:),pre(:),vflow(:),veln(:),slipn(:)
  real(8),allocatable::ag(:),bg(:),dcg(:),f0g(:),vcg(:),vwg(:),fwg(:),taug(:),sigmag(:),velg(:),taudotg(:),sigdotg(:)
  real(8),allocatable::vplg(:),cslipG(:),etavg(:),etabg(:),velnG(:),pfG(:),tmparray(:)
  integer,allocatable::listG(:)

  !variables
  real(8),allocatable::psi(:),vel(:),tau(:),sigma(:),slip(:),mu(:),rupt(:),islip(:),velp(:),cslip(:),sigma0(:),tau0(:),vslip(:)
  real(8),allocatable::taus(:),taud(:),vels(:),veld(:),slips(:),slipd(:),rake(:),lbds(:),pf(:),sigmat(:)

  real(8),allocatable::rdata(:)
  integer::lp,i,i_,j,k,kstart,kend,m,counts,interval,lrtrn,nl,ios,nmain,rk,nout(20),file_size,nrjct,ncol
  integer::hypoloc(1),load,eventcount,thec,inloc,sw,onset_step

  !controls
  logical::dilatancy,buffer,nuclei,slipping,outfield,structured,limitsigma,dcscale,slowslip,slipfinal,deepcreep,rakefromglobal,viscous
  logical::initcondfromfile,parameterfromfile,backslip,sigmaconst,forward,inverse,geofromfile,restart,latticeh,debug,bgstress,relax,analyticalinjection
  logical::opening,sorted,bingham,meshisinmeter,pressurediffusion
  character*128::fname,dum,law,input_file,problem,geofile,param,pvalue,slipmode,project,parameter_file,outdir,command,evlaw,param2(20)
  real(8)::a0,b0,dc0,sr,omega,theta,dtau,tiny,moment,wid,normal,ieta,meanmu,meanmuG,meanslip,meanslipG,moment0,mvel,mvelG,etav0,etab0
  real(8)::vc0,mu0,onset_time,tr,vw0,fw0,velmin,tauinit,intau,trelax,maxnorm,maxnormG,minnorm,minnormG,sigmainit,muinit,pfinit,tmp
  real(8)::r,vpl0,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxeventi,sparam,tmax,dtmax,tout,dtout,dtout_co,dtout_inter,dummy(10),tdil,cdil,nflow,MCNS,vref
  real(8)::cdiff,pf0,ds0,amp,mui,velinit,psinit,velmax,maxsig,minsig,v1,dipangle,crake,s,sg,errold,xhypo,yhypo,zhypo,convangle,velth,ztop

  !temporal variable

  !random_number
  integer,allocatable::seed(:)
  integer::seedsize

  !for time integration
  real(8)::x !time
  real(8),allocatable::y(:),yscal(:),dydx(:),yg(:)
  real(8)::eps_r,errmax_gb,dtinit,dtnxt,dttry,dtdid,dtmin,tp

  integer::r1,r2,r3,NVER,amari,out,kmax,loci,locj,stat,nth
  integer,allocatable::rupsG(:)

  !initialize
  icomm=MPI_COMM_WORLD
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  npd=int(sqrt(dble(np)))
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )
  allocate(rcounts(npd),displs(npd+1))

  if(my_rank==0) then
    write(*,*) 'HBI ver. 2025.10.0'
    write(*,*) '# of MPI', np
  end if
  !input file must be specified when running
  !the number of cores must be a squared number
  !example) mpirun -np 16 ./lhbiem default.in

  call get_command_argument(1,input_file,status=stat)

  open(33,file=input_file,status='old',iostat=ios)
  !if(my_rank==0) write(*,*) 'input_file',input_file
  if(ios /= 0) then
    if(my_rank==0)write(*,*) 'ERROR: Failed to open input file'
    stop
  end if 

  !get filenumber
  number=0
  if(input_file(1:2)=='in') then
    input_file=adjustl(input_file(7:))
    !write(*,*) input_file
    read(input_file,*) number
    !write(*,*) number
  end if
  if(input_file(1:11)=='examples/in') then
    input_file=adjustl(input_file(12:))
    !write(*,*) input_file
    read(input_file,*) number
    !write(*,*) number
  end if

  if(my_rank==0) then
  outdir='output'
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  !write(*, *) trim(command)
  call system(command)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr);time1=MPI_Wtime()

  !default parameters
  nmain=1000000
  eps_r=1d-4; eps_h=1d-4
  vref=1e-6
  velmax=1d7; velmin=1d-20; velth=1e-2
  tmax=1d4
  interval=0

  bgstress=.false.; sigmaconst=.false.; forward=.false.; inverse=.false.
  restart=.false.; deepcreep=.false.; limitsigma=.true.; analyticalinjection=.false.
  opening=.false.;  viscous=.false.; bingham=.false.; relax=.false.
  meshisinmeter=.false.; parameterfromfile=.false.; rakefromglobal=.false.
  debug=.false.; structured=.false.; dilatancy=.false.
  maxsig=300d0; minsig=1d0; muinit=0d0
  dtout=365*24*3600; dtout_co=1000.0; dtinit=1d0;  dtmax=1e10
  tp=86400d0
  ncol=0; noutloc=0; locid=0; cdiff=0.0; pf0=0.0
  pfinit=0d0; nflow=1.0
  geofile="default"; parameter_file="default"
  errold=1.0; ztop=0.0
  !st_bemv%fwid=1e8
  evlaw='aging'
  param_diff%bc='none'
  param_diff%bcl='Neumann';param_diff%bcr='Neumann';param_diff%bct='Neumann';param_diff%bcb='Neumann'
  param_diff%injection='none'
  param_diff%tinj=1e8;param_diff%permev=.false.

  !read input file
  call read_inputfile()

  !check inconsistency in input parameters
  if(pressurediffusion) then
    if(param_diff%injectionfromfile) then
      call input_well(param_diff)
    end if
    select case(problem)
    case('2dnh','3dht','3dnt')
      if(my_rank == 0) write(*,*) 'ERROR: pressurediffusion is impossible for problem ', problem
    stop
    end select
  end if

  tmax=tmax*365*24*3600
  dtout_inter=dtout*365*24*3600
  cdiff=1e-6*param_diff%kp0/param_diff%eta/param_diff%beta/param_diff%phi0
  if(interval==0) interval=Nstep

  if(geofile=='default') then
    write(geofile,'("examples/geo",i0,".dat")') number
  end if

  if(parameter_file=='default') then
    write(parameter_file,'("examples/param",i0,".dat")') number
  end if

  !limitsigma=.true.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !stop
  !call varscalc(NCELL,displs,vars)
  if(my_rank==0) then
    write(*,*) 'job number',number
  end if


  select case(problem)
  case('2dn','2dnh','25d')
      open(12,file=geofile,iostat=ios)
      if(ios /= 0) then
        if(my_rank==0)write(*,*) 'ERROR: Failed to open geometry file'
        stop
      end if
      nl=0
      do
        read(12,'()',end=100)
        nl = nl + 1
      end do
100   close(12)
      ncellg=nl
    if(my_rank==0) write(*,*) 'NCELLg',ncellg
  case('3dp','3dph')
    NCELLg=imax*jmax
  case('3dht','3dnt')
    if(structured) then
      NCELLg=imax*jmax*2
    else
      open(12,file=geofile,iostat=ios)
      if(ios /= 0) then
        if(my_rank==0)write(*,*) 'ERROR: Failed to open geometry file'
        stop
      end if
      nl=0
      do
        read(12,'()',end=101)
        nl = nl + 1
      end do
101   close(12)
      ncellg=nl/7
    end if
    if(my_rank==0) write(*,*) 'NCELLg',ncellg
  end select

  if(ncellg == 0) then
    if(my_rank == 0) write(*,*) 'ERROR: Ncellg is zero'
    stop
  end if

  !HACApK initialize
  lrtrn=HACApK_init(NCELLg,st_ctl,st_bemv,icomm)
  st_ctl%param(8)=20

  !allocation
  allocate(st_bemv%xcol(NCELLg),st_bemv%ycol(NCELLg),st_bemv%zcol(NCELLg))
  allocate(st_bemv%xel(NCELLg),st_bemv%xer(NCELLg),st_bemv%yel(NCELLg),st_bemv%yer(Ncellg))
  allocate(st_bemv%ang(NCELLg),st_bemv%angd(NCELLg),st_bemv%rake(NCELLg))
  allocate(st_bemv%ds(NCELLg),st_bemv%dsl(NCELLg),st_bemv%dsd(NCELLg))
  allocate(st_bemv%xs1(NCELLg),st_bemv%xs2(NCELLg),st_bemv%xs3(NCELLg),st_bemv%xs4(NCELLg))
  allocate(st_bemv%ys1(NCELLg),st_bemv%ys2(NCELLg),st_bemv%ys3(NCELLg),st_bemv%ys4(NCELLg))
  allocate(st_bemv%zs1(NCELLg),st_bemv%zs2(NCELLg),st_bemv%zs3(NCELLg),st_bemv%zs4(NCELLg))
  allocate(st_bemv%ev11(NCELLg),st_bemv%ev12(NCELLg),st_bemv%ev13(NCELLg))
  allocate(st_bemv%ev21(NCELLg),st_bemv%ev22(NCELLg),st_bemv%ev23(NCELLg))
  allocate(st_bemv%ev31(NCELLg),st_bemv%ev32(NCELLg),st_bemv%ev33(NCELLg))
  allocate(listG(NCELLg),lbds(NCELL))

  !allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg),ds(NCELLg),dsl(NCELLg),dsd(NCELLg))
  allocate(ag(NCELLg),bg(NCELLg),dcg(NCELLg),f0g(NCELLg),etavg(NCELLg),etabg(NCELLg),vcg(NCELLg),vwg(NCELLg),fwg(NCELLg))
  allocate(taug(NCELLg),sigmag(NCELLg),pfG(Ncellg),velG(NCELLg),rake(NCELLg),cslipG(NCELLg),velnG(NCELLg),tmparray(NCELLg))
  allocate(taudotg(NCELLg),sigdotg(NCELLg),vplg(NCELLg),param_diff%kpG(NCELLg),param_diff%phiG(NCELLg))
  param_diff%kpG=param_diff%kp0;param_diff%phiG=param_diff%phi0;pfG=pfinit

  st_bemv%xcol=0d0;st_bemv%ycol=0d0;st_bemv%zcol=0d0;st_bemv%ds=0d0
  st_bemv%w=ds0
  st_bemv%problem=problem

  !allocate(vmax(NCELLg),vmaxin(NcELLg))

  !mesh generation (rectangular assumed)
  if(my_rank==0) write(*,*) 'Generating mesh'
  select case(problem)
  case('2dp','2dvs')
    call coordinate2dp(NCELLg, ds0)
    st_bemv%dsl=st_bemv%ds

  case('2dph')
    call coordinate2dph()
    st_bemv%dsl=st_bemv%ds
  
  case('2dn','2dnh','25d')
    open(20,file=geofile,status='old',iostat=ios)
    if(ios /= 0) then
      if(my_rank==0)write(*,*) 'ERROR: Failed to open geometry file'
      stop
    end if
    do i=1,NCELLg
      read(20,*) st_bemv%xel(i),st_bemv%xer(i),st_bemv%yel(i),st_bemv%yer(i)
    end do
    close(20)
    call coordinate2dn()
    st_bemv%dsl=st_bemv%ds

  case('3dp')
    call coordinate3dp(imax,jmax,ds0)
    st_bemv%ds=ds0*ds0
    st_bemv%dsl=ds0
    st_bemv%dsd=ds0

  case('3dnr','3dhr')
    open(20,file=geofile,status='old',iostat=ios)
    if(ios /= 0) then
     if(my_rank==0)write(*,*) 'ERROR: Failed to open geometry file'
     stop
    end if
    do i=1,NCELLg
     read(20,*) st_bemv%xcol(i),st_bemv%ycol(i),st_bemv%zcol(i),st_bemv%ang(i),st_bemv%angd(i),st_bemv%dsl(i),st_bemv%dsd(i)
     st_bemv%ds(i)=st_bemv%dsl(i)*st_bemv%dsd(i)
    end do
    ds0=minval(st_bemv%dsl)
    close(20)

  case('3dph')
    call coordinate3ddip(imax,jmax,ds0,dipangle)
    st_bemv%ds=ds0*ds0
    st_bemv%dsl=ds0
    st_bemv%dsd=ds0

  case('3dnt','3dht')

    if(structured) then
      call coordinate3dns(NCELLg)
    else
      !.stl format
      open(12,file=geofile,iostat=ios)
      if(ios /= 0) then
        if(my_rank==0)write(*,*) 'ERROR: Failed to open geometry file'
        stop
      end if

      do while(.true.)
        read(12,*) dum
        if(dum=='facet') exit
      end do
      !write(*,*) ios
      do k=1,ncellg
        read(12,*) !outer loop
        read(12,*) dum,st_bemv%xs1(k),st_bemv%ys1(k),st_bemv%zs1(k) !vertex
        read(12,*) dum,st_bemv%xs2(k),st_bemv%ys2(k),st_bemv%zs2(k) !vertex
        read(12,*) dum,st_bemv%xs3(k),st_bemv%ys3(k),st_bemv%zs3(k) !vertex
        read(12,*) !end loop
        read(12,*) !endfacet
        read(12,*) !facet
        if(meshisinmeter) then
          st_bemv%xs1(k)=st_bemv%xs1(k)/1e3; st_bemv%ys1(k)=st_bemv%ys1(k)/1e3; st_bemv%zs1(k)=st_bemv%zs1(k)/1e3
          st_bemv%xs2(k)=st_bemv%xs2(k)/1e3; st_bemv%ys2(k)=st_bemv%ys2(k)/1e3; st_bemv%zs2(k)=st_bemv%zs2(k)/1e3
          st_bemv%xs3(k)=st_bemv%xs3(k)/1e3; st_bemv%ys3(k)=st_bemv%ys3(k)/1e3; st_bemv%zs3(k)=st_bemv%zs3(k)/1e3
        end if
        st_bemv%xcol(k)=(st_bemv%xs1(k)+st_bemv%xs2(k)+st_bemv%xs3(k))/3
        st_bemv%ycol(k)=(st_bemv%ys1(k)+st_bemv%ys2(k)+st_bemv%ys3(k))/3
        st_bemv%zcol(k)=(st_bemv%zs1(k)+st_bemv%zs2(k)+st_bemv%zs3(k))/3
      !  write(*,*)ios
        !if(my_rank==0)write(*,'(9e17.8)') xs1(k),ys1(k),zs1(k),xs2(k),ys2(k),zs2(k),xs3(k),ys3(k),zs3(k)
      end do
    end if

    !mesh format created by .msh => mkelm.c
    ! open(20,file=geofile,status='old',iostat=ios)
    ! if(ios /= 0) then
    !   if(my_rank==0)write(*,*) 'error: Failed to open geometry file'
    !   stop
    ! end if
    ! do i=1,NCELLg
    !   read(20,*) k,xs1(i),ys1(i),zs1(i),xs2(i),ys2(i),zs2(i),xs3(i),ys3(i),zs3(i),xcol(i),ycol(i),zcol(i)
    ! end do
    ! close(20)
    call evcalc(st_bemv%xs1,st_bemv%xs2,st_bemv%xs3,st_bemv%ys1,st_bemv%ys2,st_bemv%ys3,st_bemv%zs1,st_bemv%zs2,st_bemv%zs3 &
    & ,st_bemv%ev11,st_bemv%ev12,st_bemv%ev13,st_bemv%ev21,st_bemv%ev22,st_bemv%ev23,st_bemv%ev31,st_bemv%ev32,st_bemv%ev33,st_bemv%ds)
    st_bemv%dsl=sqrt(2*st_bemv%ds)
  end select

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   rake=crake
   st_bemv%rake=rake

  !nonuniform parameters from file
  if(parameterfromfile) then
    open(99,file=parameter_file,iostat=ios)
    if(ios /= 0) then
      if(my_rank==0) write(*,*) 'ERROR: Failed to open parameter file'
      stop
    end if
    read(99, *, iostat=ios) param2(1:ncol)
    allocate(values(NCELLg,ncol))
    do i=1,ncellg
      read(99,*) values(i,1:ncol)
    end do

    do k=1,ncol
      if(my_rank==0) write(*,*) param2(k)
      select case(param2(k))
      case('rake')
        st_bemv%rake(:)=values(:,k)
      case('a')
        ag(:)=values(:,k)
      case('b')
        bg(:)=values(:,k)
      case('dc')
        dcg(:)=values(:,k)
      case('f0')
        f0g(:)=values(:,k)
      case('vc')
        vcg(:)=values(:,k)
      case('vw')
        vwg(:)=values(:,k)
      case('fw')
        fwg(:)=values(:,k)
      case('etav')
        etavg(:)=values(:,k)
      case('etab')
        etabg(:)=values(:,k)
      case('tau')
        taug(:)=values(:,k)
      case('sigma')
        sigmag(:)=values(:,k)
      case('vel')
        velg(:)=values(:,k)
      case('kp')
        param_diff%kpG(:)=values(:,k)
      case('phi')
        param_diff%phiG(:)=values(:,k)
      case('taudot')
        taudotg(:)=values(:,k)
      case('sigmadot')
        sigdotg(:)=values(:,k)
      case('vpl')
        vplg(:)=values(:,k)
      case('pf')
        pfg(:)=values(:,k)
      end select
    end do

    deallocate(values)

    ! do i=1,ncellg
    !   read(99,*) rake(i),ag(i),bg(i),dcg(i),f0g(i),taug(i),sigmag(i),velg(i),taudotg(i),sigdotg(i)
    ! end do

    close(99)
  end if

  if(problem .ne. '2dn') st_bemv%rake=st_bemv%rake/180d0*pi

  if(rakefromglobal) then
    call calcrake(st_bemv%ev11,st_bemv%ev12,convangle,st_bemv%rake)
  end if

  !HACApK setting
  !lrtrn=HACApK_init(NCELLg,st_ctl,st_bemv,icomm)
  !st_ctl%param(8)=20
  !if(latticeh) st_ctl%param(8)=20
  !lrtrn=HACApK_init(NCELLg,st_ctl2,st_bemv,icomm)
  allocate(coord(NCELLg,3))

  !generate kernel (H-matrix aprrox)
  if(my_rank==0) write(*,*) 'Generating kernel'
  do i=1,NCELLg
    coord(i,1)=st_bemv%xcol(i)
    coord(i,2)=st_bemv%ycol(i)
    coord(i,3)=st_bemv%zcol(i)
  end do

  !ycol=0d0
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !if(np==1) st_ctl%param(43)=1
  st_ctl%param(8)=20
  !st_ctl%param(7)=1
  !st_ctl%param(12)=1
  select case(problem)
  case('2dp','2dn3','3dp','2dvs')
    sigmaconst=.true.
  end select
  
  if(opening) sigmaconst=.false.

  st_bemv%md='s'
  st_bemv%v='s'
  lrtrn=HACApK_generate(st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
  !if(latticeh) then
  lrtrn=HACApK_construct_LH(st_LHp_s,st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
  allocate(wws(st_leafmtxp_s%ndlfs))
  lrtrn=HACApK_gen_lattice_vector(st_vel,st_leafmtxp_s,st_ctl)
  lrtrn=HACApK_gen_lattice_vector(st_sum,st_leafmtxp_s,st_ctl)
  NCELL=st_vel%ndc
  if(.not.sigmaconst) then
    lrtrn=HACApK_init(NCELLg,st_ctl2,st_bemv,icomm)
   ! st_ctl2%param(8)=20
    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n,st_bemv,st_ctl2,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_n,st_leafmtxp_n,st_bemv,st_ctl2,coord,eps_h)
  end if
  if(opening) then
    lrtrn=HACApK_init(NCELLg,st_ctl3,st_bemv,icomm)
    st_bemv%md='open'
    st_bemv%v='s'
    lrtrn=HACApK_generate(st_leafmtxp_s2,st_bemv,st_ctl3,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_s2,st_leafmtxp_s2,st_bemv,st_ctl3,coord,eps_h)
    lrtrn=HACApK_init(NCELLg,st_ctl4,st_bemv,icomm)
    st_bemv%md='open'
    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n2,st_bemv,st_ctl4,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_n2,st_leafmtxp_n2,st_bemv,st_ctl4,coord,eps_h)
  end if

  call MPI_GATHER(ncell,1,MPI_INT,rcounts,1,MPI_INT,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  displs(1)=0
  do i=2,npd+1
    displs(i)=displs(i-1)+rcounts(i-1)
  end do
  call MPI_GATHERv(st_sum%lodc,NCELL,MPI_INTEGER,listG,rcounts,displs,MPI_INTEGER,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

  !write(*,*) my_rank,st_ctl%lpmd(33),st_ctl%lpmd(37)

  !allocate local variables
  allocate(y(3*NCELL),yscal(3*NCELL),dydx(3*NCELL))
  allocate(psi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),slip(NCELL),mu(NCELL),veln(NCELL),slipn(ncell),pf(Ncell),sigmat(Ncell))
  allocate(islip(NCELL),cslip(NCELL),sigma0(NCELL),tau0(NCELL),etav(NCELL),etab(NCELL),pre(NCELL),vflow(NCELL),vslip(NCELL))
  psi=0d0;vel=0d0;tau=0d0;sigma=0d0;slip=0d0;etav=0d0;etab=0d0;pre=0d0;vflow=0d0;vslip=0d0;pf=0d0
  allocate(a(NCELL),b(NCELL),dc(NCELL),f0(NCELL),vc(NCELL),vw(NCELL),fw(NCELL),taudot(NCELL),tauddot(NCELL),sigdot(NCELL),vpl(NCELL))
  taudot=0d0;sigdot=0d0
  allocate(param_diff%kp(ncell))

   

  !uniform frictional parameters
  a=a0; b=b0; dc=dc0; f0=mu0; vc=vc0; vpl=vpl0; vw=vw0; fw=fw0
  param_diff%kp=param_diff%kp0
  if(bingham) etab=etab0
  if(viscous) pre=1/etav0

  if(.not.backslip) then
    taudot=sr
    sigdot=0d0
  end if

  if(parameterfromfile) then
   call set_params()
  end if

  if(deepcreep) then
    if(my_rank == 0) write(*,*) 'Loading rate is calculated from deep creep'
    select case(problem)
    case('3dph')
      call taudot_3dph()
    case('2dnh','2dph')
      call taudot_2dnh()
    end select
  end if

  !call taudot_3dph()

  ! dc=2e-5
  ! f0=0.35
  ! taudot=taudot*1e-3
  ! sigdot=sigdot*1e-3

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !stop
  !max time step
  !if(load==0) dtmax=0.02d0*10d0/sr
  !write(*,*) forward
  if(forward) then
    if(my_rank==0)write(*,*) 'Debug mode: calculate stress from uniform slip'
    call forward_check()
  end if
  if(inverse) call inverse_problem()

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !restart
  if(restart) then
    if(my_rank==0)then
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname,status='old',iostat=ios)
      if(ios /= 0) then
        write(*,*) 'ERROR: monitor file is not found'
        stop
      end if

      do while(.true.)
        read(52,*,iostat=ios)k,dummy(1:9)
        kstart=k+1
        x=dummy(1)
        dtnxt=min(2*dummy(8),0.9*dummy(8)*(dummy(7)**(-0.2)))
        if(ios<0) exit
      end do
      close(52)
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname,status='old',position='append')
      write(fname,'("output/time",i0,".dat")') number
      open(50,file=fname,status='old',position='append')

      kend=k
      write(*,*) "Restarting: time step=", kstart
      write(*,*) "Time (yr)=", x/365/24/3600

      write(fname,'("output/event",i0,".dat")') number
      open(44,file=fname)
      do while(.true.)
        read(44,*,iostat=ios) k
        if(ios<0) exit
      end do
      eventcount=k
      close(44)

      write(fname,'("output/event",i0,".dat")') number
      open(44,file=fname,status='old',position='append')


      open(19,file='job.log',position='append')
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
      write(19,'(a20,i0,a6,a12,a6,a12)') 'Restarting job number=',number,'date',sys_time(1),'time',sys_time(2)
      close(19)

    end if
    call MPI_bcast(kstart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_bcast(dtnxt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    nout(1)=my_rank+100
    write(fname,'("output/vel",i0,".dat")') number
    open(nout(1),file=fname,form='unformatted',access='stream')
    inquire(nout(1),size=file_size)
    m=file_size/8
    allocate(rdata(m))
    read(nout(1)) rdata
    close(nout(1))
    do i=1,NCELL
      i_=st_sum%lodc(i)
      vel(i)=rdata(m-NCELLg+i_)
    end do
    

    write(fname,'("output/slip",i0,".dat")') number
    open(nout(1),file=fname,form='unformatted',access='stream')
    read(nout(1)) rdata
    close(nout(1))
    do i=1,NCELL
      i_=st_sum%lodc(i)
      slip(i)=rdata(m-NCELLg+i_)
    end do

    write(fname,'("output/sigma",i0,".dat")') number
    open(nout(1),file=fname,form='unformatted',access='stream')
    read(nout(1)) rdata
    close(nout(1))
    do i=1,NCELL
      i_=st_sum%lodc(i)
      sigma(i)=rdata(m-NCELLg+i_)
    end do

    write(fname,'("output/tau",i0,".dat")') number
    open(nout(1),file=fname,form='unformatted',access='stream')
    read(nout(1)) rdata
    close(nout(1))
    do i=1,NCELL
      i_=st_sum%lodc(i)
      tau(i)=rdata(m-NCELLg+i_)
    end do

    if(viscous) then
      write(fname,'("output/vflow",i0,".dat")') number
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      close(nout(1))
      do i=1,NCELL
        i_=st_sum%lodc(i)
        vflow(i)=rdata(m-NCELLg+i_)
      end do

      write(fname,'("output/vslip",i0,".dat")') number
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      close(nout(1))
      do i=1,NCELL
        i_=st_sum%lodc(i)
        vslip(i)=rdata(m-NCELLg+i_)
      end do
    end if

    if(opening) then
     write(fname,'("output/veln",i0,".dat")') number
     open(nout(1),file=fname,form='unformatted',access='stream')
     read(nout(1)) rdata
     close(nout(1))
     do i=1,NCELL
       i_=st_sum%lodc(i)
       veln(i)=rdata(m-NCELLg+i_)
     end do
     write(fname,'("output/slipn",i0,".dat")') number
     open(nout(1),file=fname,form='unformatted',access='stream')
     read(nout(1)) rdata
     close(nout(1))
     do i=1,NCELL
       i_=st_sum%lodc(i)
       slipn(i)=rdata(m-NCELLg+i_)
     end do
  
    end if

    if(pressurediffusion) then
      write(fname,'("output/pf",i0,".dat")') number
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      close(nout(1))
      do i=1,NCELL
        i_=st_sum%lodc(i)
        pf(i)=rdata(m-NCELLg+i_)
      end do
    end if

    if(param_diff%permev) then
      write(fname,'("output/kp",i0,".dat")') number
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      close(nout(1))
      do i=1,NCELL
        i_=st_sum%lodc(i)
        param_diff%kp(i)=rdata(m-NCELLg+i_)
      end do
    end if
    !write(*,*) my_rank,m

    psi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
    !write(*,*) tau
    write(fname,'("output/vel",i0,".dat")') number
    open(nout(1),file=fname,form='unformatted',access='stream',status='old',position='append')
    nout(2)=nout(1)+np
    write(fname,'("output/slip",i0,".dat")') number
    open(nout(2),file=fname,form='unformatted',access='stream',status='old',position='append')
    nout(3)=nout(2)+np
    write(fname,'("output/sigma",i0,".dat")') number
    open(nout(3),file=fname,form='unformatted',access='stream',status='old',position='append')
    nout(4)=nout(3)+np
    write(fname,'("output/tau",i0,".dat")') number
    open(nout(4),file=fname,form='unformatted',access='stream',status='old',position='append')
    nout(5)=nout(4)+np
    write(fname,'("output/EQslip",i0,".dat")') number
    open(nout(5),file=fname,form='unformatted',access='stream',status='replace')

    if(viscous) then
      nout(6)=nout(5)+1
      write(fname,'("output/vflow",i0,".dat")') number
      open(nout(6),file=fname,form='unformatted',access='stream',status='old',position='append')
      nout(7)=nout(6)+1
      write(fname,'("output/vslip",i0,".dat")') number
      open(nout(7),file=fname,form='unformatted',access='stream',status='old',position='append')
    end if

    if(opening) then
      nout(8)=nout(7)+1
      write(fname,'("output/veln",i0,".dat")') number
      open(nout(8),file=fname,form='unformatted',access='stream',status='old',position='append')
      nout(9)=nout(8)+1
      write(fname,'("output/slipn",i0,".dat")') number
      open(nout(9),file=fname,form='unformatted',access='stream',status='old',position='append')
    end if

    if(pressurediffusion) then
      nout(10)=nout(9)+1
      write(fname,'("output/pf",i0,".dat")') number
      open(nout(10),file=fname,form='unformatted',access='stream',status='old',position='append')
    end if

    if(param_diff%permev) then
      nout(11)=nout(10)+1
      write(fname,'("output/kp",i0,".dat")') number
      open(nout(10),file=fname,form='unformatted',access='stream',status='replace',position='append')
    end if

    rankloc=-1
    if(my_rank<npd) then
      do k=1,noutloc
      do i=1,ncell
        if(st_sum%lodc(i)==locid(k)) then
          rankloc(k)=my_rank
          loc_(k)=i
          !write(*,*) rankloc(k),loc_(k)
          write(*,*)  st_bemv%xcol(locid(k)),st_bemv%ycol(locid(k)),st_bemv%zcol(locid(k))
          exit
        end if
      end do
     end do
    end if

    do k=1,10
      if(my_rank==rankloc(k)) then
        write(fname,'("output/local",i0,"-",i0,".dat")') number,locid(k)
        open(53+k,file=fname,position='append')
        !write(*,*) my_rank,locid(k)
     end if
   end do

    s=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      s=s+st_bemv%ds(i_)
    end do

    call MPI_reduce(s,sG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(mvelG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr); time2=MPI_Wtime()
    if(my_rank==0) write(*,*) 'Finished all initial processing, time (s)=',time2-time1
    time1=MPI_Wtime()
  
  !no restart--------------------------------------------------------------------------------------------
  else 

    !setting initial condition

    !uniform values
    sigma=sigmainit
    tau=tauinit
    vel=velinit
    veln=0
    pf=pfinit
  
    !if(my_rank==0) write(*,*) tau
    if(muinit.ne.0d0) tau=sigma*muinit

    !non-uniform initial stress from file
    if(parameterfromfile) then
      do k=1,ncol
      select case(param2(k))
        case('tau')
          do i=1,NCELL
            i_=st_sum%lodc(i)
            tau(i)=taug(i_)
          end do
        case('sigma')
          do i=1,NCELL
            i_=st_sum%lodc(i)
            sigma(i)=sigmag(i_)
          end do
        case('vel')
          do i=1,NCELL
            i_=st_sum%lodc(i)
            vel(i)=velg(i_)
          end do
        case('kp')
          do i=1,NCELL
            i_=st_sum%lodc(i)
            param_diff%kp(i)=param_diff%kpG(i_)
          end do
        end select
      end do
    end if

    if(viscous) vflow=pre*tau**nflow

    if(bgstress) then
      call initcond()
    end if

    sigma0=sigma
    tau0=tau
    
    if(param_diff%injection=='pressure') then
      pfG(Ncellg/2)=param_diff%pinj
      !sigma(Ncellg/2)=sigmainit-pf(Ncellg/2)
      !tau(Ncellg/2)=muinit*sigma(Ncellg/2)
    end if

    mu=tau/sigma
    psi=a*dlog(2*vref/vel*sinh(tau/sigma/a))

    if(evlaw=='mCNS') psi=tau/sigma-a*dlog(vel/vref)
    slip=0d0
    slipn=0d0

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
    x=0d0
    kstart=1
    kend=0
    !if(sorted) then
      call output_coord_sorted()
    !else
    !  call output_coord()
    !end if
    if(my_rank==0) then
      nout(1)=80
      write(fname,'("output/vel",i0,".dat")') number
      open(nout(1),file=fname,form='unformatted',access='stream',status='replace')
      nout(2)=nout(1)+1
      write(fname,'("output/slip",i0,".dat")') number
      open(nout(2),file=fname,form='unformatted',access='stream',status='replace')
      nout(3)=nout(2)+1
      write(fname,'("output/sigma",i0,".dat")') number
      open(nout(3),file=fname,form='unformatted',access='stream',status='replace')
      nout(4)=nout(3)+1
      write(fname,'("output/tau",i0,".dat")') number
      open(nout(4),file=fname,form='unformatted',access='stream',status='replace')
      nout(5)=nout(4)+1
      write(fname,'("output/EQslip",i0,".dat")') number
      open(nout(5),file=fname,form='unformatted',access='stream',status='replace')
      
      if(viscous) then
        nout(6)=nout(5)+1
        write(fname,'("output/vflow",i0,".dat")') number
        open(nout(6),file=fname,form='unformatted',access='stream',status='replace')
        nout(7)=nout(6)+1
        write(fname,'("output/vslip",i0,".dat")') number
        open(nout(7),file=fname,form='unformatted',access='stream',status='replace')
      end if

      if(opening) then
        nout(8)=nout(7)+1
        write(fname,'("output/veln",i0,".dat")') number
        open(nout(8),file=fname,form='unformatted',access='stream',status='replace')
        nout(9)=nout(8)+1
        write(fname,'("output/slipn",i0,".dat")') number
        open(nout(9),file=fname,form='unformatted',access='stream',status='replace')
      end if

      if(pressurediffusion) then
        nout(10)=nout(9)+1
        write(fname,'("output/pf",i0,".dat")') number
        open(nout(10),file=fname,form='unformatted',access='stream',status='replace')
        if(param_diff%permev) then
          nout(11)=nout(10)+1
          write(fname,'("output/kp",i0,".dat")') number
          open(nout(11),file=fname,form='unformatted',access='stream',status='replace')
        end if
      end if      

      write(fname,'("output/time",i0,".dat")') number
      open(50,file=fname)
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname)
      write(fname,'("output/event",i0,".dat")') number
      open(44,file=fname,status='replace')
      open(19,file='job.log',position='append')
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
      write(19,'(a20,i0,a6,a12,a6,a12,a4,i4)') 'Starting job number=',number,'date',sys_time(1),'time',sys_time(2),'np',np
      close(19)
    end if

    !local data output
    rankloc=-1
    if(my_rank<npd) then
      do k=1,noutloc
      do i=1,ncell
        if(st_sum%lodc(i)==locid(k)) then
          rankloc(k)=my_rank
          loc_(k)=i
          !write(*,*) rankloc(k),loc_(k)
          write(*,*)  st_bemv%xcol(locid(k)),st_bemv%ycol(locid(k)),st_bemv%zcol(locid(k))
          exit
        end if
      end do
     end do
    end if

    do k=1,10
      if(my_rank==rankloc(k)) then
        write(fname,'("output/local",i0,"-",i0,".dat")') number,locid(k)
        open(53+k,file=fname,position='append')
        !write(*,*) my_rank,locid(k)
     end if
   end do

    s=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      s=s+st_bemv%ds(i_)
    end do

    call MPI_reduce(s,sG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(mvelG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)


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

    ! meanslip=sum(slip*ds)/sum(ds)
    meanslip=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meanslip=meanslip+slip(i)*st_bemv%ds(i_)
    end do
    !call MPI_ALLREDUCE(meanslip,meanslipG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meanslip,meanslipG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meanslipG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meanslipG=meanslipG/sg

    ! meanmu=sum(mu*ds)
    meanmu=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meanmu=meanmu+mu(i)*st_bemv%ds(i_)
    end do
    !call MPI_ALLREDUCE(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meanmuG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meanmuG=meanmuG/sg

    call MPI_BARRIER(MPI_COMM_WORLD,ierr); time2=MPI_Wtime()
    if(my_rank==0) write(*,*) 'Finished all initial processing, time(s)=',time2-time1
    time1=MPI_Wtime()

    !output intiial condition
    k=0
    errmax_gb=0d0
    dtdid=0d0
    if(my_rank==0) then
      write(50,'(i7,f19.4)') k,x
      write(*,*) 'initial Vmax (m/s)',mvelG
      write(*,*) 'initial mu_avg',meanmuG
      call output_monitor()
    end if
    call output_field_sorted()
    
    dtnxt = dtinit
  end if

  dtout=dtout_inter
  tout=dtout
  !tout=20*365*24*60*60
  rk=0

  !outv=1d-6
  slipping=.false.
  eventcount=0
  sw=0
  timer=0d0
  timeH=0d0
  !time2=MPI_Wtime()

  !$omp parallel do
  do i=1,NCELL
    y(3*i-2) = psi(i)
    y(3*i-1) = tau(i)
    y(3*i) = sigma(i)
    !if(my_rank==53)write(*,*) psi(i),tau(i),sigma(i)
  end do
  !$omp end parallel do
  !stop

  do k=kstart,NSTEP
    !parallel computing for Runge-Kutta
    !write(*,*) dc(1)/mvelG
    !dttry = min(dtnxt,0.05*dc(1)/mvelG)
    dttry = dtnxt
    !time3=MPI_Wtime()
    call rkqs(y,dydx,x,dttry,eps_r,dtdid,dtnxt,errmax_gb,nrjct)
    !call rkqs2(y,dydx,x,dttry,eps_r,errold,dtdid,dtnxt,errmax_gb,nrjct)
    !time4=MPI_Wtime()
    !timer=timer+time4-time3

    if(analyticalinjection) then
      select case(problem)
      case('2dp','2dn')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          !along-fault pressurediffusion
          sigma0(i)=sigma0(i)+y(3*i)-sigma(i)
          !y(3*i)=sigma0(i)-pf0*erfc(xcol(i_)/(2*sqrt(cdiff*x)))
          if(param_diff%injection=='pressure') then
            y(3*i)=sigma0(i)-param_diff%pinj*erfc(abs(st_bemv%xcol(i_)-st_bemv%xcol(ncellg/2))/(2*sqrt(cdiff*x)))
          else if(param_diff%injection=='flowrate') then
            y(3*i)=sigma0(i)-pf1d(param_diff%pinj,cdiff,x,st_bemv%xcol(i_))
          end if
          !y(3*i)=max(minsig,sigma0(i)-pf0*erfc(xcol(i_)/(2*sqrt(cdiff*x))))
          !homogenous pressurediffusion
          !y(3*i)=sigma0(i)-pf2d(pf0,cdiff,x,st_bemvs%xcol(i_),st_bemv%ycol(i_)-0.2)
        end do

      case('3dp')
        do i=1,NCELL
          i_=st_sum%lodc(i)
           !along-fault pressurediffusion
          sigma0(i)=sigma0(i)+y(3*i)-sigma(i)
          if(param_diff%injection=='flowrate') then
            y(3*i)=sigma0(i)-pf2d(pf0,cdiff,x,st_bemv%xcol(i_),st_bemv%zcol(i_))
          end if
          !y(3*i)=max(minsig,sigma0(i)-pf2d(pf0,cdiff,x,xcol(i_),zcol(i_)))
          !homogenous pressurediffusion
          !r=sqrt(st_bemv%xcol(i_)**2+st_bemv%zcol(i_)**2)
          !y(3*i)=max(minsig,sigma0(i)-pf0*erfc(r/(2*sqrt(cdiff*x))/r))
        end do
      end select 
    end if

    !compute physical values for control and output
    !$omp parallel do
    do i = 1, NCELL
      psi(i) = y(3*i-2)
      tau(i) = y(3*i-1)
      sigma(i)=y(3*i)

      slip(i)=slip(i)+(vel(i)+vflow(i))*dtdid*0.5d0 !2nd order
      vel(i)= 2*vref*exp(-psi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
      if(evlaw=='mCNS') vel(i) = vref*dexp((tau(i)/sigma(i)-psi(i))/a(i))
      if(viscous) then
        vslip(i)=vslip(i)+vflow(i)*dtdid*0.5d0
        vflow(i)=pre(i)*tau(i)**nflow
        vslip(i)=vslip(i)+vflow(i)*dtdid*0.5d0
      end if
      if(opening) then
        veln(i)=vel(i)*(max(0.0,minsig-sigma(i))*1.0+min(0.0,maxsig-sigma(i))*1.0)
        !veln(i)=-(sigma(i)-sigma0(i))*1e-12
        slipn(i)=slipn(i)+veln(i)*dtdid
      end if
      
      slip(i)=slip(i)+(vel(i)+vflow(i))*dtdid*0.5d0
      mu(i)=tau(i)/sigma(i)
      sigmat(i)=sigma(i)+pf(i)
    end do
    !$omp end parallel do

    if(pressurediffusion) then
      if(param_diff%permev) then
        do i=1,ncell
          tmp=-vel(i)/param_diff%kL*(param_diff%kp(i)-param_diff%kpmax)-(param_diff%kp(i)-param_diff%kpmin)/param_diff%kT
          param_diff%kp(i)=param_diff%kp(i)+dtdid*tmp
        end do
        call MPI_GATHERv(param_diff%kp,NCELL,MPI_REAL8,tmparray,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)   
        do i=1, NCELLg
          i_=listG(i)
          param_diff%kpG(i_)=tmparray(i)
        end do
      end if
      call MPI_GATHERv(pf,NCELL,MPI_REAL8,tmparray,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      if(my_rank==0) then
        do i=1, NCELLg
          i_=listG(i)
          pfG(i_)=tmparray(i)
        end do
       !write(*,*) minval(param_diff%kpG),maxval(param_diff%kpG), maxval(param_diff%phiG),maxval(pfG),minval(pfG)
        select case(problem)
        case('2dp','2dn')
          call diffusion2dwop(pfG,dtdid,ds0,x,dtnxt,param_diff)
        case('3dp','3dph')
          call diffusion3dwop(imax,jmax,pfG,dtdid,ds0,x,dtnxt,param_diff)
        end select
      !end if
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pfG,NCELLg,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtnxt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      do i=1,ncell
        i_=st_sum%lodc(i)
        pf(i)=pfG(i_)
        sigma(i)=sigmat(i)-pf(i)
        if(limitsigma) then
          if(sigma(i)<minsig) sigma(i)=minsig
          if(sigma(i)>maxsig) sigma(i)=maxsig
        end if
        y(3*i)=sigma(i)
      end do
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr); time3=MPI_Wtime()

    mvel=maxval(abs(vel))
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(mvel,mvelG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    !call MPI_reduce(mvel,rank_mvel,1,MPI_REAL8_INT,MPI_MAXLOC,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(mvelG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    !call MPI_bcast(rank_mvel,1,MPI_INTEGER,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    !if(my_rank==0) write(*,*) rank_mvel
    !if(my_rank==rank_mvel) then
    !  mvel_loc=maxloc(abs(vel))
    !  mvel_loc_=st_sum%lodc(mvel_loc(1))
    !end if
    !call MPI_BCAST(mvel_loc_,1,MPI_INTEGER,rank_mvel,MPI_COMM_WORLD,ierr)
    !if(my_rank==0) write(*,*) mvel_loc_


    maxnorm=maxval(sigma)
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(maxnorm,maxnormG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(maxnormG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    minnorm=minval(sigma)
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(minnorm,minnormG,1,MPI_REAL8,MPI_MIN,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(minnormG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    ! meanslip=sum(slip*ds)/sum(ds)
    meanslip=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meanslip=meanslip+abs(slip(i))*st_bemv%ds(i_)
    end do
    !call MPI_ALLREDUCE(meanslip,meanslipG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meanslip,meanslipG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meanslipG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meanslipG=meanslipG/sg

    ! meanmu=sum(mu*ds)
    meanmu=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meanmu=meanmu+mu(i)*st_bemv%ds(i_)
    end do
    !call MPI_ALLREDUCE(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meanmuG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meanmuG=meanmuG/sg
    !if(outfield.and.(my_rank.lt.npd)) call output_field()
 
    outfield=.false.
    !event list
    if(.not.slipping) then
      if(mvelG>velth) then
        outfield=.true.
        slipping=.true.
        dtout=dtout_co
        tout=x+dtout
        eventcount=eventcount+1
        moment0=meanslipG
        islip=slip
        call get_mvelloc(mvel_loc)
        hypoloc=mvel_loc(1)
        onset_time=x
        onset_step=k
        !tout=onset_time

      end if
    end if
    !
    if(slipping) then
      if(mvelG<0.5*velth) then
        outfield=.true.
        slipping=.false.
        dtout=dtout_inter
        tout=x+dtout
        moment=meanslipG-moment0
        !eventcount=eventcount+1
        !end of an event
        if(my_rank==0) then
          write(44,'(i0,i8,f17.2,f14.4,i8)') eventcount,onset_step,onset_time,(log10(moment*rigid*sg)+5.9)/1.5,hypoloc
        end if
        cslip=abs(slip-islip)
        call output_EQslip_sorted()

      end if
      !   vmaxevent=max(vmaxevent,maxval(vel))
      !   !write(53,'(i6,4e16.6)') !k,x-onset_time,sum(slip-islip),sum(vel),sum(acg**2)
      !   !if(x-onset_time>lapse) then
      !   !  lapse=lapse+dlapse
      !   !end if
    end if

    if(k==1) outfield=.true.

    if(mod(k,interval)==0) then
      outfield=.true.
    end if
    if(x>tmax) then
      outfield=.true.
    end if
    !if(k>3800.and.k<4200.and.mod(k,10)==0) outfield=.true.
    if(x>tout) then
      outfield=.true.
      tout=x+dtout
      !if(x<300*356*23*3600) outfield=.false.
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
      end if
      
      call output_field_sorted()

    end if

    if(my_rank==0 .and. k>kend) then
      call output_monitor()
    end if
    do j=1,10
      if(my_rank==rankloc(j) .and. k>kend) then
        call output_local(53+j,loc_(j))
      end if
    end do
    time4=MPI_Wtime()
    timer=timer+time4-time3

    !stop controls
    if(mvelG>velmax) then
      if(my_rank == 0) write(*,*) 'Slip rate above vmax at time step=', k
      exit
    end if
    if(mvelG<velmin) then
      if(my_rank == 0) write(*,*) 'Slip rate below vmin at time step=', k
      exit
    end if
    if(x>tmax) then
      if(my_rank == 0) write(*,*) 'Time exceeds tmax at time step=', k
      exit
    end if
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr); time2= MPI_Wtime()
  200  if(my_rank==0) then
  write(*,*) 'time(s)', time2-time1,timer,timeH
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
end if
!if(my_rank==0) write(19,'(a20,i0,f16.2)')'Finished job number=',number,time2-time1
Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
select case(problem)
case('2dp','2dn3','3dp','2dvs')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_s)
case('2dn','2dph','3dnt','3dht','3dnr','3dhr','3dph')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_s)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_n)
end select
lrtrn=HACApK_finalize(st_ctl)
Call MPI_FINALIZE(ierr)
stop
contains
  !------------output-----------------------------------------------------------!
  subroutine output_monitor()
    implicit none
    time2=MPi_Wtime()
    write(52,'(i7,f19.4,7e16.5,i4,f16.4)')k,x,log10(mvelG),meanslipG,meanmuG,maxnormG,minnormG,errmax_gb,dtdid,nrjct,time2-time1
  end subroutine

  subroutine output_local(nf,loc_)
    implicit none
    integer,intent(in)::nf,loc_
    write(nf,'(i7,f19.4,9e16.5)')k,x,log10(abs(vel(loc_))),slip(loc_),tau(loc_),sigma(loc_),mu(loc_),psi(loc_),vflow(loc_),veln(loc_),slipn(loc_)
  end subroutine

  subroutine get_mvelloc(mvel_loc)
    implicit none
    integer::i,i_,nn!,rcounts(npd),displs(npd+1)
    integer::mvel_loc(1)
    real(8)::velG(NCELLg),velG2(NCELLg)
    call MPI_GATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    if(my_rank==0) then
    do i=1, NCELLg
      i_=listG(i)
      velG2(i_)=velG(i)
    end do

    mvel_loc=maxloc(abs(velG2))
    end if
    return
  end subroutine

  ! subroutine output_field()
  !   implicit none
  !   integer::nn,rcounts(npd),displs(npd+1)
  !   integer::mvel_loc(1)
  !   real(8)::velG(NCELLg),tauG(NCELLg),sigmaG(Ncellg),slipG(ncellg),vflowg(ncellg),vslipG(ncellg),velnG(ncellg),slipnG(Ncellg),pfG(Ncellg)
  !   call MPI_GATHER(ncell,1,MPI_INT,rcounts,1,MPI_INT,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

  !   displs(1)=0
  !   do nn=2,npd+1
  !    displs(nn)=displs(nn-1)+rcounts(nn-1)
  !   end do

  !   call MPI_GATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   call MPI_GATHERv(tau,NCELL,MPI_REAL8,tauG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   call MPI_GATHERv(sigma,NCELL,MPI_REAL8,sigmaG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   call MPI_GATHERv(slip,NCELL,MPI_REAL8,slipG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   if(viscous) then
  !     call MPI_GATHERv(vflow,NCELL,MPI_REAL8,vflowG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !     call MPI_GATHERv(vslip,NCELL,MPI_REAL8,vslipG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   end if
  !   if(opening) then
  !     call MPI_GATHERv(veln,NCELL,MPI_REAL8,velnG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !     call MPI_GATHERv(slipn,NCELL,MPI_REAL8,slipnG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   end if
  !   if(pressurediffusion) then
  !     call MPI_GATHERv(pf,NCELL,MPI_REAL8,pfG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   end if
  !   !write(*,*) velG

  !   if(my_rank==0) then
  !     mvel_loc=maxloc(abs(velG))
  !     write(nout(1)) velG
  !     write(nout(2)) slipG
  !     write(nout(3)) sigmaG
  !     write(nout(4)) tauG
  !     if(viscous) then
  !       write(nout(6)) vflowG
  !       write(nout(7)) vslipG
  !     end if
  !     if(opening) then
  !       write(nout(8)) velnG
  !       write(nout(9)) slipnG
  !     end if
  !     if(pressurediffusion) then
  !       write(nout(10)) pfG     
  !     end if
  !   end if
  ! end subroutine

  subroutine output_field_sorted()
    implicit none
    !integer::nn,rcounts(npd),displs(npd+1)
    integer::mvel_loc(1),i_!,listG(NCELLg),i_
    real(8)::velG(NCELLg),tauG(NCELLg),sigmaG(Ncellg),slipG(ncellg),vflowg(ncellg),vslipG(ncellg),velnG(ncellg),slipnG(Ncellg)
    real(8)::velG2(NCELLg),tauG2(NCELLg),sigmaG2(Ncellg),slipG2(ncellg),vflowg2(ncellg),vslipG2(ncellg),velnG2(ncellg),slipnG2(Ncellg)

    call MPI_GATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_GATHERv(tau,NCELL,MPI_REAL8,tauG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_GATHERv(sigma,NCELL,MPI_REAL8,sigmaG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_GATHERv(slip,NCELL,MPI_REAL8,slipG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    if(viscous) then
      call MPI_GATHERv(vflow,NCELL,MPI_REAL8,vflowG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      call MPI_GATHERv(vslip,NCELL,MPI_REAL8,vslipG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    end if
    if(opening) then
      call MPI_GATHERv(veln,NCELL,MPI_REAL8,velnG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      call MPI_GATHERv(slipn,NCELL,MPI_REAL8,slipnG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    end if

    if(my_rank==0) then
      do i=1, NCELLg
        i_=listG(i)
        velG2(i_)=velG(i)
        slipG2(i_)=slipG(i)
        sigmaG2(i_)=sigmaG(i)
        tauG2(i_)=tauG(i)
      end do
      write(nout(1)) velG2
      write(nout(2)) slipG2
      write(nout(3)) sigmaG2
      write(nout(4)) tauG2

      if(viscous) then
        do i=1, NCELLg
          i_=listG(i)
          vflowG2(i_)=vflowG(i)
          vslipG2(i_)=vslipG(i)
        end do
        write(nout(6)) vflowG2
        write(nout(7)) vslipG2
      end if

      if(opening) then
        do i=1, NCELLg
          i_=listG(i)
          velnG2(i_)=velnG(i)
          slipnG2(i_)=slipnG(i)
        end do
        write(nout(8)) velnG2
        write(nout(9)) slipnG2
      end if
      if(pressurediffusion) then
        write(nout(10)) pfG
        if(param_diff%permev) then
          write(nout(11)) param_diff%kpG
        end if
      end if
    end if

  end subroutine

  ! subroutine output_EQslip()
  !   implicit none
  !   !integer::nn,rcounts(npd),displs(npd+1)
  !   !integer::mvel_loc(1),listG(NCELLg),i_
  !   real(8)::cslipG(ncellg),cslipG2(ncellg)

  !   call MPI_GATHERv(cslip,NCELL,MPI_REAL8,cslipG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
  !   if(my_rank==0) then
  !     write(nout(5)) cslipG
  !   end if
  !   return
  ! end subroutine
  subroutine output_EQslip_sorted()
    implicit none
    !integer::nn,rcounts(npd),displs(npd+1)
    !integer::mvel_loc(1),listG(NCELLg),i_
    real(8)::cslipG(ncellg),cslipG2(ncellg)

    call MPI_GATHERv(cslip,NCELL,MPI_REAL8,cslipG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

    if(my_rank==0) then
      do i=1, NCELLg
        i_=listG(i)
        cslipG2(i_)=cslipG(i)
      end do
      write(nout(5)) cslipG2
    end if

    return
  end subroutine
!   subroutine output_coord()
!     implicit none
!     integer::nn
!     nout(1)=100
!     nout(2)=101

!     if(my_rank==0) then
!       write(fname,'("output/xyz",i0,".dat")') number
!       open(nout(1),file=fname,status='replace')
!       write(fname,'("output/ind",i0,".dat")') number
!       open(nout(2),file=fname,status='replace')
!       do i=1,ncell
!         i_=st_sum%lodc(i)
!         write(nout(2),*) i_
!         write(nout(1),'(3e15.6)') st_bemv%xcol(i_),st_bemv%ycol(i_),st_bemv%zcol(i_)
!       end do
!       close(nout(1))
!       close(nout(2))
!     end if
 
!     Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     do nn=1,npd-1
!       if(my_rank==nn) then
!       write(fname,'("output/xyz",i0,".dat")') number
!       open(nout(1),file=fname,position='append')
!       write(fname,'("output/ind",i0,".dat")') number
!       open(nout(2),file=fname,position='append')
!       do i=1,ncell
!         i_=st_sum%lodc(i)
!         write(nout(2),*) i_
!         write(nout(1),'(3e15.6)') st_bemv%xcol(i_),st_bemv%ycol(i_),st_bemv%zcol(i_)
!       end do
!       close(nout(1))
!       close(nout(2))
!       end if
!       Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!     end do
!     !Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!     return
! end subroutine
subroutine output_coord_sorted()
    implicit none
    nout(1)=100
    nout(2)=101

      if(my_rank==0) then
        write(fname,'("output/xyz",i0,".dat")') number
        open(nout(1),file=fname,status='replace')
        write(fname,'("output/ind",i0,".dat")') number
        open(nout(2),file=fname,status='replace')
        do i=1,ncellg
          write(nout(1),'(3e15.6)') st_bemv%xcol(i),st_bemv%ycol(i),st_bemv%zcol(i)
        end do
        close(nout(1))
        close(nout(2))
      end if
    return
end subroutine

  !------------coordinate-----------------------------------------------------------!
  subroutine coordinate2dp(NCELLg,ds0)
    implicit none
    integer,intent(in)::NCELLg
    real(8),intent(in)::ds0
    integer::i,j,k

    !flat fault with element size ds
    do i=1,NCELLg
      st_bemv%ds(i)=ds0
      st_bemv%xel(i)=(i-1)*ds0
      st_bemv%xer(i)=i*ds0
      st_bemv%xcol(i)=0.5d0*(st_bemv%xel(i)+st_bemv%xer(i))
      !write(14,'(3e16.6)') xcol(i),xel(i),xer(i)
    enddo
    !close(14)
    return
  end subroutine

  subroutine coordinate2dph()
    implicit none
    integer::i,j,k
    !planar fault with element size ds and dipangle=dipangle
    do i=1,NCELLg
      st_bemv%ds(i)=ds0
      st_bemv%xel(i)=(i-1)*ds0*cos(dipangle*pi/180)
      st_bemv%xer(i)=i*ds0*cos(dipangle*pi/180)
      st_bemv%yel(i)=(i-1)*ds0*sin(dipangle*pi/180)
      st_bemv%yer(i)=i*ds0*sin(dipangle*pi/180)
      st_bemv%xcol(i)=0.5d0*(st_bemv%xel(i)+st_bemv%xer(i))
      st_bemv%ycol(i)=0.5d0*(st_bemv%yel(i)+st_bemv%yer(i))
      st_bemv%ang(i)=datan2(st_bemv%yer(i)-st_bemv%yel(i),st_bemv%xer(i)-st_bemv%xel(i))
      !write(14,'(3e16.6)') xcol(i),xel(i),xer(i)
    enddo
    !close(14)
    return
  end subroutine

  subroutine coordinate2dn()
    implicit none
    integer::i
    do i=1,NCELLg
      st_bemv%ds(i)=sqrt((st_bemv%xer(i)-st_bemv%xel(i))**2+(st_bemv%yer(i)-st_bemv%yel(i))**2)
      st_bemv%ang(i)=datan2(st_bemv%yer(i)-st_bemv%yel(i),st_bemv%xer(i)-st_bemv%xel(i))
      st_bemv%xcol(i)=0.5d0*(st_bemv%xel(i)+st_bemv%xer(i))
      st_bemv%ycol(i)=0.5d0*(st_bemv%yel(i)+st_bemv%yer(i))
      !write(*,*) ds(i),ang(i)
    enddo
    return
  end subroutine

  subroutine coordinate3dp(imax,jmax,ds0)
    implicit none
    integer,intent(in)::imax,jmax
    real(8),intent(in)::ds0
    real(8)::dx,dz
    integer::i,j,k

    dx=ds0
    dz=ds0
    do i=1,imax
      do j=1,jmax
        k=(i-1)*jmax+j
        st_bemv%xcol(k)=(i-imax/2-0.5d0)*dx
        st_bemv%zcol(k)=-(j-jmax/2-0.5d0)*dz
        st_bemv%xs1(k)=st_bemv%xcol(k)+0.5d0*dx
        st_bemv%xs2(k)=st_bemv%xcol(k)-0.5d0*dx
        st_bemv%xs3(k)=st_bemv%xcol(k)-0.5d0*dx
        st_bemv%xs4(k)=st_bemv%xcol(k)+0.5d0*dx
        st_bemv%zs1(k)=st_bemv%zcol(k)+0.5d0*dz
        st_bemv%zs2(k)=st_bemv%zcol(k)+0.5d0*dz
        st_bemv%zs3(k)=st_bemv%zcol(k)-0.5d0*dz
        st_bemv%zs4(k)=st_bemv%zcol(k)-0.5d0*dz
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
      real(8)::xc,yc,zc,amp,yr(0:jmax),zr(0:jmax),stangle


        !dipangle=dipangle*pi/180d0
        stangle=0d0*pi/180d0
         k=0
         yr(0)=0d0
         zr(0)=ztop
         !nb=int(50.0*jmax/320.0)
         !if(my_rank==0) write(*,*)nb
         do j=1,jmax
         yr(j)=yr(j-1)-ds0*cos(dipangle*pi/180d0)
         zr(j)=zr(j-1)-ds0*sin(dipangle*pi/180d0)
         end do
         do i=1,imax
           do j=1,jmax
           k=k+1
             st_bemv%xcol(k)=(i-imax/2-0.5d0)*ds0
             !ys1(k)=yr(j-1)
             !ys2(k)=yr(j-1)
             !ys3(k)=yr(j)
             !ys4(k)=yr(j)
             st_bemv%ycol(k)=0.5d0*(yr(j-1)+yr(j)) !-(j-0.5d0)*ds0*cos(dipangle)
             st_bemv%zcol(k)=0.5d0*(zr(j-1)+zr(j)) !-(j-0.5d0)*ds0*sin(dipangle)
             !xcol(k)=(i-imax/2-0.5d0)*ds0
             !zcol(k)=(j-jmax/2-0.5d0)*ds0
             !xs1(k)=xcol(k)+0.5d0*ds0
             !xs2(k)=xcol(k)-0.5d0*ds0
             !xs3(k)=xcol(k)-0.5d0*ds0
             !xs4(k)=xcol(k)+0.5d0*ds0
             !zs1(k)=zcol(k)+0.5d0*ds0*sin(dipangle)
             !zs2(k)=zcol(k)+0.5d0*ds0*sin(dipangle)
             !zs3(k)=zcol(k)-0.5d0*ds0*sin(dipangle)
             !zs4(k)=zcol(k)-0.5d0*ds0*sin(dipangle)
             !ys1(k)=ycol(k)+0.5d0*ds0*cos(dipangle)*sin(stangle)
             !ys2(k)=ycol(k)-0.5d0*ds0*cos(dipangle)*sin(stangle)
             !ys3(k)=ycol(k)-0.5d0*ds0*cos(dipangle)*sin(stangle)
             !ys4(k)=ycol(k)+0.5d0*ds0*cos(dipangle)*sin(stangle)
             st_bemv%angd(k)=datan2(zr(j-1)-zr(j),yr(j-1)-yr(j))
             st_bemv%ang(k)=0d0
             !write(*,*)angd(k)
             if(my_rank==0)write(111,*)st_bemv%xcol(k),st_bemv%ycol(k),st_bemv%zcol(k)
             end do
         end do

      return
    end subroutine coordinate3ddip

    subroutine coordinate3dns(NCELLg)
      implicit none
      integer,intent(in)::NCELLg
      integer::i,j,k
      integer,parameter::ndata=256,nskip=22
      real(8)::xc,yc,zc,amp
      real(4)::xl(ndata+1,ndata+1)

      open(30,file='examples/roughsurf.txt')
      do k=1,ndata+1
        read(30,*) xl(k,1:ndata+1)
      end do
      close(30)

      do i=1,imax
        do j=1,jmax
          k=(i-1)*jmax+j
          !xcol(k)=(i-imax/2-0.5d0)*ds0
          !zcol(k)=-(j-0.5d0)*ds0-0.001d0
          xc=(i-imax/2-0.5)*ds0
          yc=-(j-0.5)*ds0*cos(dipangle*pi/180d0)
          zc=-(j-0.5)*ds0*sin(dipangle*pi/180d0)
  
          st_bemv%xs1(2*k-1)=xc-0.5d0*ds0; st_bemv%xs2(2*k-1)=xc+0.5d0*ds0; st_bemv%xs3(2*k-1)=xc-0.5d0*ds0
          st_bemv%zs1(2*k-1)=zc+0.5d0*ds0; st_bemv%zs2(2*k-1)=zc+0.5d0*ds0; st_bemv%zs3(2*k-1)=zc-0.5d0*ds0
          st_bemv%ys1(2*k-1)=yc-0.5d0*ds0; st_bemv%ys2(2*k-1)=yc-0.5d0*ds0; st_bemv%ys3(2*k-1)=yc+0.5d0*ds0

          ! ys2(2*k-1)=xl(i+1,j)
          ! ys3(2*k-1)=xl(i,j+1)
          st_bemv%ys1(2*k-1)=xl(i+nskip,j+nskip); st_bemv%ys2(2*k-1)=xl(i+1+nskip,j+nskip); st_bemv%ys3(2*k-1)=xl(i+nskip,j+1+nskip)
  
          st_bemv%xs2(2*k)=xc+0.5d0*ds0; st_bemv%xs1(2*k)=xc+0.5d0*ds0; st_bemv%xs3(2*k)=xc-0.5d0*ds0
          st_bemv%zs2(2*k)=zc-0.5d0*ds0; st_bemv%zs1(2*k)=zc+0.5d0*ds0; st_bemv%zs3(2*k)=zc-0.5d0*ds0
          st_bemv%ys2(2*k)=yc+0.5d0*ds0; st_bemv%ys1(2*k)=yc-0.5d0*ds0; st_bemv%ys3(2*k)=yc+0.5d0*ds0
          st_bemv%ys2(2*k)=xl(i+1+nskip,j+1+nskip); st_bemv%ys1(2*k)=xl(i+1+nskip,j+nskip); st_bemv%ys3(2*k)=xl(i+nskip,j+1+nskip)
  
        end do
      end do
      do k=1,ncellg
        st_bemv%xcol(k)=(st_bemv%xs1(k)+st_bemv%xs2(k)+st_bemv%xs3(k))/3.d0
        st_bemv%ycol(k)=(st_bemv%ys1(k)+st_bemv%ys2(k)+st_bemv%ys3(k))/3.d0
        st_bemv%zcol(k)=(st_bemv%zs1(k)+st_bemv%zs2(k)+st_bemv%zs3(k))/3.d0
        !write(*,*) xcol(k),ycol(k),zcol(k)
      end do
  
      ! open(30,file='roughsurf.txt')
      ! do k=0,255
      !   read(30,*) xl(k,0:255)
      ! end do
      ! close(30)
      ! amp=0.000d0
      ! if(my_rank.eq.0) open(32,file='tmp')
      ! do i=1,NCELLg
      !   ! xcol(i)=(xs1(i)+xs2(i)+xs3(i))/3.d0
      !   ! zcol(i)=(zs1(i)+zs2(i)+zs3(i))/3.d0
      
      !   ! j=int((xs1(i)+10)*102.4)
      !   ! k=int(-102.4*zs1(i))
      !   ! ys1(i)=xl(j,k)*amp
      !   ! j=int((xs2(i)+10)*102.4)
      !   ! k=int(-102.4*zs2(i))
      !   ! ys2(i)=xl(j,k)*amp
      !   ! j=int((xs3(i)+10)*102.4)
      !   ! k=int(-102.4*zs3(i))
      !   ! ys3(i)=xl(j,k)*amp
      !   ! ycol(i)=(ys1(i)+ys2(i)+ys3(i))/3.d0

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


      !   if(my_rank.eq.0) write(32,*) xcol(i),ycol(i),zcol(i)
      ! end do
  
      return
    end subroutine coordinate3dns

  subroutine evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33,ds)
    !calculate ev for each element
    implicit none
    real(8),intent(in)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    real(8),intent(out)::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:),ds(:)
    real(8)::rr,vba(0:2),vca(0:2),tmp1,tmp2,tmp3,theta
    theta=0d0 
    !open(50,file='cas3rake')
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

      ev21(k) = ev32(k)*ev13(k)-ev33(k)*ev12(k)
      ev22(k) = ev33(k)*ev11(k)-ev31(k)*ev13(k)
      ev23(k) = ev31(k)*ev12(k)-ev32(k)*ev11(k)

      !tmp1=(ev12(k)-atan(theta)*ev11(k))/(atan(theta)*ev21(k)-ev22(k))
      !tmp2=atan(tmp1)
      !tmp1=atan2(-ev11(k),ev12(k))
      !if(my_rank==0) write(50,*) ev11(k),ev12(k),90-tmp1/pi*180
      !rake(k)=pi/2-tmp1+convangle*pi/180


      tmp1=vba(0)*vba(0)+vba(1)*vba(1)+vba(2)*vba(2)
      tmp2=vca(0)*vca(0)+vca(1)*vca(1)+vca(2)*vca(2)
      tmp3=vba(0)*vca(0)+vba(1)*vca(1)+vba(2)*vca(2)
      ds(k)=0.5d0*sqrt(tmp1*tmp2-tmp3*tmp3)
      !if(my_rank==0) write(*,*)ev21(k),ev22(k),ev23(k)
    end do

  end subroutine

  subroutine calcrake(ev11,ev12,convangle,rake)
    real(8),intent(in)::ev11(:),ev12(:),convangle
    real(8),intent(out)::rake(:)
    real(8)::tmp1

    do k=1,ncellg
      tmp1=atan2(-ev11(k),ev12(k))
      !if(my_rank==0) write(50,*) ev11(k),ev12(k),90-tmp1/pi*180
      rake(k)=pi/2-tmp1+convangle*pi/180
    end do
    return
  end subroutine

  subroutine taudot_3dph()
    real(8)::rate
    !open(111,file='tmp')
    do i=1,NCELL
      i_=st_sum%lodc(i)
      taudot(i)=okada_load(st_bemv%xcol(i_),st_bemv%ycol(i_),st_bemv%zcol(i_),-ds0*imax/2,ds0*imax/2,ds0*jmax,0d0,pi/2.0,0d0)*vpl0
    end do
    !close(111)
  end subroutine

  subroutine taudot_2dnh()
    real(8)::rate
    integer::nload=1000
    character(128)::vv
    do i=1,NCELL
      i_=st_sum%lodc(i)
      vv="s"
      taudot(i)=load2dnh(st_bemv%xcol(i_),st_bemv%ycol(i_),st_bemv%xer(Nload),st_bemv%yer(Nload),pi*dipangle/180,vv)*vpl0
      vv="n"
      sigdot(i)=load2dnh(st_bemv%xcol(i_),st_bemv%ycol(i_),st_bemv%xer(Nload),st_bemv%yer(Nload),pi*dipangle/180,vv)*vpl0

      taudot(i)=taudot(i)+sr*0.5*sin(2*st_bemv%ang(i_))
      sigdot(i)=sigdot(i)+sr*sin(st_bemv%ang(i_))**2

      !write(*,*) xcol(i_),taudot(i),sigdot(i)
    end do
  end subroutine

  subroutine initcond()
    implicit none
    real(8),parameter::rho=1.6e3,g=9.80
    real(8)::sxx0,sxy0,syy0
  do i=1,ncell
    i_=st_sum%lodc(i)
    syy0=min(50.0,rho*g*st_bemv%ycol(i_)/1e3)
    sxy0=0.0
    sxx0=syy0*3.0
    tau(i)=sxy0*cos(2*st_bemv%ang(i_))+0.5*(sxx0-syy0)*sin(2*st_bemv%ang(i_))
    sigma(i)=sin(st_bemv%ang(i_))**2*sxx0+cos(st_bemv%ang(i_))**2*syy0+sxy0*sin(2*st_bemv%ang(i_))
  end do
    mu=tau/sigma
    vel=velinit
    psi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
    slip=0d0

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
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),psitmp(NCELL),vflowtmp(NCELL),velntmp(NCELL)
    real(8) :: sum_gs(NCELL),sum_gn(NCELL)
    real(8) :: time3,time4,c1, c2, c3, arg,arg2,c,g,tauss,Arot(3,3),p(6),fac,sxx0,sxy0,syy0
    integer :: i, j, nc,ierr,lrtrn,i_

    !if(latticeh) then
    vflowtmp=0d0

    !$omp parallel do
    do i = 1, NCELL
      psitmp(i) = y(3*i-2)
      tautmp(i) = y(3*i-1)
      sigmatmp(i) = y(3*i)
      veltmp(i) = 2*vref*dexp(-psitmp(i)/a(i))*dsinh(tautmp(i)/sigmatmp(i)/a(i))
      if(evlaw=='mCNS') veltmp(i) = vref*dexp((tautmp(i)/sigmatmp(i)-psitmp(i))/a(i))
      !if(my_rank==0)write(*,*) veltmp(i)
    enddo
    !$omp end parallel do

    if(viscous) then
      do i=1,NCELL
        vflowtmp(i)=pre(i)*tautmp(i)**nflow
      end do
    end if

    !opening/wear
    if(opening) then
      do i=1,NCELL
        velntmp(i)=max(0.0,minsig-sigmatmp(i))*abs(veltmp(i))
        !velntmp(i)=-(sigmatmp(i)-sigma0(i))*1e-12
      end do
    end if

    !matrix-vector mutiplation
    if(backslip) then
      st_vel%vs=veltmp+vflowtmp-vpl
    else
      st_vel%vs=veltmp+vflowtmp
    end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr);time3=MPI_Wtime()
    call HACApK_adot_lattice_hyp(st_sum,st_LHp_s,st_ctl,wws,st_vel)
    if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') then
      sum_gs(:)=st_sum%vs(:)/ds0
    else
      sum_gs(:)=st_sum%vs(:)
    end if
   

    if(sigmaconst) then
      sum_gn=0d0
    else
      call HACAPK_adot_lattice_hyp(st_sum,st_LHP_n,st_ctl2,wws,st_vel)
      if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') then
        sum_gn(:)=st_sum%vs(:)/ds0
      else
        sum_gn(:)=st_sum%vs(:)
      end if

      if(opening) then
        st_vel%vs=velntmp
        call HACAPK_adot_lattice_hyp(st_sum,st_LHP_s2,st_ctl3,wws,st_vel)
        if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') then
          sum_gs(:)=sum_gs(:)+st_sum%vs(:)/ds0
        else
          sum_gs(:)=sum_gs(:)+st_sum%vs(:)
        end if

        call HACAPK_adot_lattice_hyp(st_sum,st_LHP_n2,st_ctl4,wws,st_vel)
        if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') then
          sum_gn(:)=sum_gn(:)+st_sum%vs(:)/ds0
        else
          sum_gn(:)=sum_gn(:)+st_sum%vs(:)
        end if
      end if
    end if
    !time4=MPI_Wtime()
    !timeH=timeH+time4-time3

    !$omp parallel do
    do i=1,NCELL
      sum_gn(i)=sum_gn(i)+sigdot(i)
      sum_gs(i)=sum_gs(i)+taudot(i)
      if(relax) then
        sum_gn(i)=sum_gn(i)-(sigmatmp(i)-sigma0(i))/trelax
        sum_gs(i)=sum_gs(i)-(tautmp(i)-tau0(i))/trelax
      end if
  
      call deriv(sum_gs(i),sum_gn(i),psitmp(i),tautmp(i),sigmatmp(i),veltmp(i),a(i),b(i),dc(i),f0(i),vc(i),vw(i),fw(i),etab(i),dydx(3*i-2),dydx(3*i-1),dydx(3*i))
    enddo
    !$omp end parallel do
    !write(*,*) dydx
    return
  end subroutine

  subroutine deriv(sum_gs,sum_gn,psitmp,tautmp,sigmatmp,veltmp,a,b,dc,f0,vc,vw,fw,etab,dpsidt,dtaudt,dsigdt)
    implicit none
    real(8)::fss,dvdtau,dvdsig,dvdpsi,mu,psiss,dcv,f
    real(8),parameter::V0=0.5,n=4,Cr=1.0
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs,sum_gn,psitmp,tautmp,sigmatmp,veltmp,a,b,dc,f0,vc,vw,fw,etab
    real(8),intent(out)::dpsidt,dtaudt,dsigdt
    dsigdt=sum_gn
    if(limitsigma)then
      if(sigmatmp<minsig) dsigdt=0
      if(sigmatmp>maxsig) dsigdt=0
    end if
   
    !regularized slip law
    !fss=mu0+(a(i_)-b(i_))*dlog(abs(veltmp(i))/vref)
    !fss=fw(i_)+(fss-fw(i_))/(1.d0+(veltmp(i)/vw(i_))**8)**0.125d0 !flash heating
    !dpsidt(i)=-abs(veltmp(i))/dc(i_)*(abs(tautmp(i))/sigmatmp(i)-fss)

    select case(evlaw)
    case('aging')
      if(b==0) then
        dpsidt=0d0
      else 
        dpsidt=b/dc*vref*dexp((f0-psitmp)/b)-b*abs(veltmp)/dc
      end if
    case('agingfh')
      if(b==0) then
        dpsidt=0d0
      else 
        dpsidt=b/dc*vref*dexp((f0-psitmp)/b)-b*abs(veltmp)/dc
      end if
    case('slip')
      fss=f0+(a-b)*dlog(abs(veltmp)/vref)
      dpsidt=-abs(veltmp)/dc*(abs(tautmp)/sigmatmp-fss)
    case('flashheating')
      fss=f0+(a-b)*dlog(abs(veltmp)/vref)
      fss=fw+(fss-fw)/(1.d0+(veltmp/vw)**8)**0.125d0 !flash heating
      psiss=a*dlog(2*vref/veltmp*sinh(fss/a))
      !dcv=dc*(1+(veltmp/vref)**2)**(0.5)
      !dcv=dc*(1+log((veltmp/vref)+1))
      dpsidt=-veltmp/dcv*(psitmp-psiss)
    case('cutoff')
      if(b==0) then
        dpsidt=0d0
      else 
        dpsidt=b/dc*vref*(1.0+abs(veltmp)/vc)*dexp((f0-psitmp)/b)-b*abs(veltmp)/dc
      end if
    case('mCNS')
      if(b==0) then
        dpsidt=0d0
      else
        dpsidt=b*Vref/dc*((psitmp-f0)/b/(MCNS+1)+1)*(((psitmp-f0)/b/(MCNS+1)+1)**(-MCNS-1)-veltmp/vref)
      end if
    end select
    
    dvdtau=2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)/(a*sigmatmp)
    dvdsig=-2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)*tautmp/(a*sigmatmp**2)
    dvdpsi=-veltmp/a
    dtaudt=(sum_gs-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt))/(1d0+0.5d0*rigid/vs*dvdtau)
    if(bingham) dtaudt=(sum_gs-(etab+0.5d0*rigid/vs)*(dvdpsi*dpsidt+dvdsig*dsigdt))/(1d0+(etab+0.5d0*rigid/vs)*dvdtau)
    !write(*,*) rigid/vs*dvdtau
    ! if(veltmp<=0d0) then
    !   dvdtau=2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)/(a*sigmatmp)
    !   dvdsig=-2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)*tautmp/(a*sigmatmp**2)
    !   !sign ok?
    !   !dvdpsi=2*vref*exp(-psitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
    !   dvdpsi=-veltmp/a
    !   dtaudt=sum_gs-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt)
    !   dtaudt=dtaudt/(1d0+0.5d0*rigid/vs*dvdtau)
    ! end if
  end subroutine
  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,hdid,hnext,errmax_gb,nrjct)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
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
    integer,intent(out)::nrjct
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    integer :: i,ierr,loc
    real(8) :: errmax,h,xnew,htemp,dtmin
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4,kp=0.08,tiny=1e-20

    nrjct=0
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

      !check nan
      do i=1,size(ytemp)
        if(abs(yerr(i)/ytemp(i))/eps>errmax) errmax=abs(yerr(i)/ytemp(i))/eps
        if(ytemp(i)-1 ==ytemp(i)) errmax=10.0
        !errmax=errmax+yerr(3*i-2)**2
      end do
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      !call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      !errmax_gb=sqrt(errmax_gb/NCELLg)/eps
      call MPI_BARRIER(MPI_COMM_WORLD,ierr); time3=MPI_Wtime()
      call MPI_reduce(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      call MPI_bcast(errmax_gb,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
      time4=MPI_Wtime()
      timer=timer+time4-time3
      !if(my_rank==0)write(*,*) h,errmax,errmax_gb



      ! call MPI_reduce(errmax,errmax_gb,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      ! errmax_gb=sqrt(errmax_gb/NCELLg)/eps
      ! call MPI_bcast(errmax_gb,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)



      !if(my_rank==0)write(*,*) h,errmax_gb
      !if(h<0.25d0*ds0/vs)exit
      if((errmax_gb<1.d0).and.(errmax_gb>tiny)) then
        exit
      end if

      nrjct=nrjct+1
      if(errmax_gb>tiny) then
        h=max(0.5d0*h,SAFETY*h*(errmax_gb**PSHRNK))
      else
        h=0.5*h
      end if

      xnew=x+h
      if(xnew-x<1.d-15) then
        if(my_rank.eq.0) write(*,*) 'ERROR: Runge-Kutta method did not converge'
        stop
      end if

    end do

    hdid=h
    x=x+h
    y(:)=ytemp(:)

    hnext=min(2*h,SAFETY*h*(errmax_gb**PGROW))
    hnext=min(hnext,dtmax)
    if(x+hnext>tmax) hnext=tmax-x+1e-5

  end subroutine
  !---------------------------------------------------------------------
  subroutine rkqs2(y,dydx,x,htry,eps,errold,hdid,hnext,errmax_gb,nrjct)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    real(8),intent(in)::htry,eps,errold
    real(8),intent(inout)::y(:),x,dydx(:)
    real(8),intent(out)::hdid,hnext,errmax_gb !hdid: resulatant dt hnext: htry for the next
    integer,intent(out)::nrjct
    integer :: i,ierr,loc
    real(8) :: errmax,h,xnew,htemp,dtmin,tmp
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4,kp=-0.08

    nrjct=0
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

      do i=1,size(ytemp)
        if(abs(yerr(i)/ytemp(i))/eps>errmax) errmax=abs(yerr(i)/ytemp(i))/eps
        !errmax=errmax+yerr(3*i-2)**2
      end do
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      !call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      !errmax_gb=sqrt(errmax_gb/NCELLg)/eps
      call MPI_BARRIER(MPI_COMM_WORLD,ierr); time3=MPI_Wtime()
      call MPI_reduce(errmax,errmax_gb,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      call MPI_bcast(errmax_gb,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
      time4=MPI_Wtime()
      timer=timer+time4-time3
      !if(my_rank==0)write(*,*) h,errmax,errmax_gb



      ! call MPI_reduce(errmax,errmax_gb,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
      ! errmax_gb=sqrt(errmax_gb/NCELLg)/eps
      ! call MPI_bcast(errmax_gb,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)



      !if(my_rank==0)write(*,*) h,errmax_gb
      !if(h<0.25d0*ds0/vs)exit
      if((errmax_gb<1.d0).and.(errmax_gb>1d-15)) then
        exit
      end if

      nrjct=nrjct+1
      if(errmax_gb>1d-15) then
        h=max(0.5d0*h,SAFETY*h*(errmax_gb**PSHRNK))
      else
        h=0.5*h
      end if



      xnew=x+h
      if(xnew-x<1.d-15) then
        if(my_rank.eq.0)write(*,*)'ERROR: Runge-Kutta method did not converge'
        stop
      end if

    end do
    if(nrjct>0) then
      tmp=h**2/htry
    end if
    tmp=h
    hnext=min(2*h,tmp*(errmax_gb**PGROW)*(errold**kp))
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
    real(8),intent(in)::y(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    integer ::i
    real(8) :: ak1(3*NCELL),ak2(3*NCELL),ak3(3*NCELL),ak4(3*NCELL),ak5(3*NCELL),ak6(3*NCELL),ytemp(3*NCELL)
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

  subroutine rk2(y,x,h,yout,yerr)
    implicit none
    real(8),intent(in)::y(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    real(8)::ytemp(3*NCELL),y1(3*NCELL),ak1(3*NCELL),ak2(3*NCELL)
    integer::i

    !1st-order solution
    call derivs(x, y, ak1)
    !$omp parallel do
    do i=1,size(y)
      y1(i)=y(i)+h*ak1(i)
    end do

    !2nd-order solution
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+0.5d0*h*ak1(i)
    end do
    call derivs(x+0.5*h, ytemp, ak2)

    !$omp parallel do
    do i=1,size(y)
      yout(i)=y(i)+h*ak2(i)
      yerr(i)=yout(i)-y1(i)
    end do

    return
  end subroutine

  subroutine forward_check()
    implicit none
    real(8)::rr,lc,ret1(NCELLg),ret2(NCELLg),vec(NCELLg)
    integer::p
    ret1=0d0
    ret2=0d0

    vec=1d0
    !vec(2)=-1d0
    !do i=1,NCELL
    !  i_=st_sum%lodc(i)
    !  if(ycol(i_)>5d0) vec(i)=1d0
    !end do
    write(fname,'("output/stress",i0)') number
    open(29,file=fname)

    ! select case(problem)
    ! case('2dn')
    !   !slip from file
    !   ! open(45,file='../fd2d/rupt2.dat')
    !   ! do i=1,NCELLg
    !   !   read(45,*) a(i),vel(i),b(i)
    !   ! end do
    !
    !   st_bemv%v='xx'
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xx,st_bemv,st_ctl,a,vel)
    !   st_bemv%v='xy'
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xy,st_bemv,st_ctl,b,vel)
    !   st_bemv%v='yy'
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yy,st_bemv,st_ctl,dc,vel)
    !   if(my_rank==0) then
    !     do i=1,NCELLg
    !       taudot(i)=0.5d0*(a(i)-dc(i))*dsin(-2*ang(i))+b(i)*dcos(-2*ang(i))
    !       sigdot(i)=-(0.5d0*(a(i)+dc(i))-0.5d0*(a(i)-dc(i))*dcos(2*ang(i))-b(i)*dsin(2*ang(i)))
    !       write(29,'(4e16.4)') xcol(i),ang(i),taudot(i),sigdot(i)
    !     end do
    !   end if
    ! case('3dp')
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,a,vel)
    !
    !   if(my_rank==0) then
    !     do i=1,NCELLg
    !       write(29,'(3e16.4)') xcol(i),zcol(i),a(i)
    !     end do
    !   end if
    ! case('3dn','3dh')
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_s2,st_bemv,st_ctl,a,vel)
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_d2,st_bemv,st_ctl,b,vel)
    !   lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n2,st_bemv,st_ctl,dc,vel)
    !   if(my_rank==0) then
    !     do i=1,NCELLg
    !       write(29,'(6e16.4)') xcol(i),ycol(i),zcol(i),a(i),b(i),dc(i)
    !     end do
    !   end if
    ! case('3dnt','3dht','3dnr','3dhr')
      !lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_s,st_bemv,st_ctl,ret1,vec)
      st_vel%vs=vec
      call HACApK_adot_lattice_hyp(st_sum,st_LHp_s,st_ctl,wws,st_vel)
      ret1(:)=st_sum%vs(:)
      if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') ret1(:)=ret1(:)/ds0

      if(.not.sigmaconst) then
        !lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n,st_bemv,st_ctl,ret2,vec)
        st_vel%vs=vec
        call HACApK_adot_lattice_hyp(st_sum,st_LHp_n,st_ctl2,wws,st_vel)
        ret2(:)=st_sum%vs(:)
        if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') ret2(:)=ret2(:)/ds0
      end if


      !lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n,st_bemv,st_ctl,ret2,vec)
      if(my_rank==0) then
        do i=1,NCELL
          i_=st_sum%lodc(i)
          write(29,'(6e16.4)') st_bemv%xcol(i_),st_bemv%ycol(i_),st_bemv%zcol(i_),vec(i_),ret1(i),ret2(i)
          ! write(29,'(2i4,3e16.4)')i,j,a(i),b(i),dc(i)
        end do
      end if

    !end select
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
          write(29,'(2e16.4)') st_bemv%xcol(i),sigdot(i)
        end do
      end if
    case('3dht')
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
          write(29,'(4e16.4)') st_bemv%xcol(i),st_bemv%ycol(i),st_bemv%zcol(i),sigdot(i)
        end do
      end if
    end select
    Call MPI_FINALIZE(ierr)
    stop
  end subroutine

  FUNCTION pf2d(pf0,alpha,time,x,z)
    implicit none
    real(8)::pf2d,time,x,z,alpha,pf0
    pf2d=pf0*expint(1,((x-0.0)**2+(z-0.0)**2)/4/alpha/time)
    return
  end function

  FUNCTION pf1d(pf0,alpha,time,x)
    implicit none
    real(8)::pf1d,time,x,z,alpha,pf0
    pf1d=pf0*sqrt(time)*(exp(-x**2/(4*alpha*time))/sqrt(pi)-abs(x)/sqrt(4*alpha*time)*erfc(abs(x)/sqrt(4*alpha*time)))
    return
  end function

   FUNCTION pf1d2(pf0,alpha,time,x)
    implicit none
    real(8)::pf1d2,time,x,z,alpha,pf0
    pf1d2=pf0*erfc(abs(x)/sqrt(4*alpha*time))
    return
  end function

  FUNCTION expint(n,x)
    INTEGER ::n
    REAL(8) ::expint,x
    integer,parameter:: MAXIT=100
    real(8),parameter::EPS=1d-7,FPMIN=1d-30,EULER=.5772156649
    INTEGER ::i,ii,nm1
    REAL(8) ::a,b,c,d,del,fact,h,psi
    nm1=n-1 
    if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1))) then
        stop
    else if(n.eq.0)then
    expint=exp(-x)/x 
    else if(x.eq.0.) then
    expint=1./nm1 
    else if(x.gt.1.)then
        b=x+n   
        c=1./FPMIN
        d=1./b
        h=d
    
    
    do i=1,MAXIT
        a=-i*(nm1+i)
        b=b+2.
        d=1./(a*d+b)
        c=b+a/c
        del=c*d
        h=h*del 
    if(abs(del-1.).lt.EPS)then
        expint=h*exp(-x)
    return 
    endif
    enddo
    
    else 
        if(nm1.ne.0)then
        expint=1./nm1
    else
        expint=-log(x)-EULER 
    endif
    fact=1.
    
    do i=1,MAXIT 
        fact=-fact*x/i 
        if(i.ne.nm1)then
    del=-fact/(i-nm1) 
    else
    psi=-EULER
    do ii=1,nm1
    psi=psi+1./ii 
    enddo
    del=fact*(-log(x)+psi)
    endif
      expint=expint+del
      if(abs(del).lt.abs(expint)*EPS) return 
    enddo 
    endif
    return 
  end 
  !parameters
  subroutine set_params()
    implicit none
    do k=1,ncol
      !write(*,*) param2(k)
      select case(param2(k))
      case('a')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          a(i)=ag(i_)
        end do
      case('b')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          b(i)=bg(i_)
        end do      
      case('dc')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          dc(i)=dcg(i_)
        end do
      case('f0')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          f0(i)=f0g(i_)
        end do
      case('vc')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          vc(i)=vcg(i_)
        end do
      case('fw')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          fw(i)=fwg(i_)
        end do
      case('vw')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          vw(i)=vwg(i_)
        end do
      case('etav')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          etav(i)=etavg(i_)
        end do
        pre=1.0/etav
      case('taudot')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          taudot(i)=taudotg(i_)
        end do
      case('sigmadot')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          sigdot(i)=sigdotg(i_)
        end do
      case('vpl')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          vpl(i)=vplg(i_)
        end do
      case('pf')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          pf(i)=pfg(i_)
        end do
      case('etab')
        do i=1,NCELL
          i_=st_sum%lodc(i)
          etab(i)=etabg(i_)
        end do
      end select
    end do
  end subroutine

  subroutine read_inputfile()
  implicit none
  i=0
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
    case('f0')
      read (pvalue,*) mu0
    case('vref')
      read (pvalue,*) vref
    case('vc')
      read (pvalue,*) vc0
    case('vw')
      read (pvalue,*) vw0
    case('fw')
      read (pvalue,*) fw0
    case('MCNS')
      read (pvalue,*) MCNS
    case('etav')
      read (pvalue,*) etav0
    case('wid')
      read (pvalue,*) wid
    case('nflow')
      read (pvalue,*) nflow
    case('etab')
      read (pvalue,*) etab0
    case('sr')
      read (pvalue,*) sr
    case('vpl')
      read (pvalue,*) vpl0
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
    case('muinit')
      read (pvalue,*) muinit
    case('pfinit')
      read (pvalue,*) pfinit
    case('dtinit')
      read (pvalue,*) dtinit
    case('tmax')
      read (pvalue,*) tmax
    case('dtout')
      read (pvalue,*) dtout
    case('dtout_co')
      read (pvalue,*) dtout_co
    case('velth')
      read (pvalue,*) velth
    case('dtmax')
      read (pvalue,*) dtmax
    case('eps_r')
      read (pvalue,*) eps_r
    case('eps_h')
      read (pvalue,*) eps_h
    case('backslip')
      read (pvalue,*) backslip
    case('limitsigma')
      read (pvalue,*) limitsigma
    case('sigmaconst')
      read(pvalue,*) sigmaconst
    case('opening')
      read(pvalue,*) opening
    case('forward')
      read(pvalue,*) forward
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
    case('convangle')
      read(pvalue,*) convangle
    case('dipangle')
      read(pvalue,*) dipangle
    case('ztop')
      read(pvalue,*) ztop
    case('fwid')
      read(pvalue,*) st_bemv%fwid
    case('tdil')
      read(pvalue,*) tdil
    case('cdil')
      read(pvalue,*) cdil
    case('cdiff')
      read(pvalue,*) cdiff
    case('injection')
      read(pvalue,*) param_diff%injection
    case('permev')
      read(pvalue,*) param_diff%permev
    case('beta')
      read(pvalue,*) param_diff%beta
    case('phi')
      read(pvalue,*) param_diff%phi0
    case('eta')
      read(pvalue,*) param_diff%eta
    case('kp')
      read(pvalue,*) param_diff%kp0
    case('kpmin')
      read(pvalue,*) param_diff%kpmin
    case('kpmax')
      read(pvalue,*) param_diff%kpmax
    case('kL')
      read(pvalue,*) param_diff%kL
    case('kT')
      read(pvalue,*) param_diff%kT
    case('qinj')
      read(pvalue,*) param_diff%qinj
    case('pinj')
      read(pvalue,*) param_diff%pinj
    case('bc')
      read(pvalue,*) param_diff%bc
    case('bcl')
      read(pvalue,*) param_diff%bcl
    case('bcr')
      read(pvalue,*) param_diff%bcr
    case('bct')
      read(pvalue,*) param_diff%bct
    case('bcb')
      read(pvalue,*) param_diff%bcb
    case('pbcl')
      read(pvalue,*) param_diff%pbcl
    case('pbcr')
      read(pvalue,*) param_diff%pbcr
    case('pbct')
      read(pvalue,*) param_diff%pbct
    case('pbcb')
      read(pvalue,*) param_diff%pbcb
    case('qbcl')
      read(pvalue,*) param_diff%qbcl
    case('qbcr')
      read(pvalue,*) param_diff%qbcr
    case('qbct')
      read(pvalue,*) param_diff%qbct
    case('qbcb')
      read(pvalue,*) param_diff%qbcb
    case('tinj')
      read(pvalue,*) param_diff%tinj
    case('injectionfromfile')
      read(pvalue,*) param_diff%injectionfromfile
    case('injection_file')
      read(pvalue,'(a)') param_diff%injection_file
    case('restart')
      read(pvalue,*) restart
    case('parameterfromfile')
      read(pvalue,*) parameterfromfile
    case('parameter_file')
      read(pvalue,'(a)') parameter_file
    case('evlaw')
      read(pvalue,'(a)') evlaw
    case('trelax')
      read(pvalue,*) trelax
    case('relax')
      read(pvalue,*) relax
    case('deepcreep')
      read(pvalue,*) deepcreep
    case('bgstress')
      read(pvalue,*) bgstress
    case('rakefromglobal')
      read(pvalue,*) rakefromglobal
    case('structured')
      read(pvalue,*) structured
    case('viscous')
      read(pvalue,*) viscous
    case('analyticalinjection')
      read(pvalue,*) analyticalinjection
    case('pressurediffusion')
      read(pvalue,*) pressurediffusion
    case('bingham')
      read(pvalue,*) bingham
    case('meshisinmeter')
      read(pvalue,*) meshisinmeter
    case('parameter_file_ncol')
      read(pvalue,*) ncol
    case('outloc')
      noutloc=noutloc+1
      read(pvalue,*) locid(noutloc)
    case default
      if(my_rank==0) write(*,*) 'WARNING: ', trim(param), ' is an unknown parameter'
    end select
  end do
  close(33)
  end subroutine
end program
