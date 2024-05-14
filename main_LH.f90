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
  integer::counts2,icomm,np,npd,ierr,my_rank,npgl,loc,rankloc,loc_,rank_mvel,mvel_loc(1),mvel_loc_
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
  real(8),allocatable::a(:),b(:),dc(:),f0(:),fw(:),vw(:),vc(:),taudot(:),tauddot(:),sigdot(:)
  real(8),allocatable::ag(:),bg(:),dcg(:),f0g(:),taug(:),sigmag(:),velg(:),taudotg(:),sigdotg(:)

  !variables
  real(8),allocatable::psi(:),vel(:),tau(:),sigma(:),disp(:),mu(:),rupt(:),idisp(:),velp(:),cslip(:)
  real(8),allocatable::taus(:),taud(:),vels(:),veld(:),disps(:),dispd(:),rake(:)

  real(8),allocatable::rdata(:)
  integer::lp,i,i_,j,k,kstart,kend,m,counts,interval,lrtrn,nl,ios,nmain,rk,nout,nout2,nout3,nout4,nout5,file_size,nrjct
  integer,allocatable::locid(:)
  integer::hypoloc(1),load,eventcount,thec,inloc,sw

  !controls
  logical::aftershock,buffer,nuclei,slipping,outfield,slipevery,limitsigma,dcscale,slowslip,slipfinal,deepcreep
  logical::initcondfromfile,parameterfromfile,backslip,sigmaconst,forward,inverse,geofromfile,restart,latticeh,debug,bgstress
  character*128::fname,dum,law,input_file,problem,geofile,param,pvalue,slipmode,project,parameter_file,outdir,command,evlaw
  real(8)::a0,b0,dc0,sr,omega,theta,dtau,tiny,moment,wid,normal,ieta,meanmu,meanmuG,meandisp,meandispG,moment0,mvel,mvelG
  real(8)::vc0,mu0,onset_time,tr,vw0,fw0,velmin,tauinit,intau,trelax,maxnorm,maxnormG,minnorm,minnormG,sigmainit,muinit
  real(8)::r,vpl,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxeventi,sparam,tmax,dtmax,tout,dtout,dummy(10)
  real(8)::alpha,ds0,amp,mui,velinit,psinit,velmax,maxsig,minsig,v1,dipangle,crake,s,sg,errold,xhypo,yhypo,zhypo

  !temporal variable

  !random_number
  integer,allocatable::seed(:)
  integer::seedsize

  !for time integration
  real(8)::x !time
  real(8),allocatable::y(:),yscal(:),dydx(:),yg(:)
  real(8)::eps_r,errmax_gb,dtinit,dtnxt,dttry,dtdid,dtmin,tp,fwid

  integer::r1,r2,r3,NVER,amari,out,kmax,loci,locj,stat,nth
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

  open(33,file=input_file,status='old',iostat=ios)
  !if(my_rank==0) write(*,*) 'input_file',input_file
  if(ios /= 0) then
    if(my_rank==0)write(*,*) 'error: Failed to open input file'
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
  write(*, *) trim(command)
  call system(command)
  end if

  call MPI_BARRIER(MPI_COMM_WORLD,ierr);time1=MPI_Wtime()

  !default parameters
  nmain=1000000
  eps_r=1d-4
  eps_h=1d-4
  velmax=1d7
  velmin=1d-16
  law='d'
  tmax=1d4
  interval=0
  bgstress=.false.
  nuclei=.false.
  slipevery=.false.
  sigmaconst=.false.
  forward=.false.
  inverse=.false.
  slipfinal=.false.
  restart=.false.
  deepcreep=.false.
  limitsigma=.false.
  maxsig=300d0
  minsig=0.2d0
  amp=0d0
  muinit=0d0
  dtout=365*24*3600
  dtinit=1d0
  tp=86400d0
  trelax=1d18
  project="none"
  initcondfromfile=.false.
  parameterfromfile=.false.
  errold=1.0
  fwid=1e8
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
    case('muinit')
      read (pvalue,*) muinit
    case('dtinit')
      read (pvalue,*) dtinit
    case('tmax')
      read (pvalue,*) tmax
    case('dtout')
      read (pvalue,*) dtout
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
    case('dipangle')
      read(pvalue,*) dipangle
    case('fwid')
      read(pvalue,*) fwid
    case('restart')
      read(pvalue,*) restart
    case('latticeh')
      read(pvalue,*) latticeh
    case('parameterfromfile')
      read(pvalue,*) parameterfromfile
    case('parameter_file')
      read(pvalue,'(a)') parameter_file
    case('evlaw')
      read(pvalue,'(a)') evlaw
    case('loc')
      read(pvalue,*) loc
    case('debug')
      read(pvalue,*) debug
    case('deepcreep')
      read(pvalue,*) deepcreep
    case('bgstress')
      read(pvalue,*) bgstress
    end select
  end do
  close(33)
  tmax=tmax*365*24*3600
  if(interval==0) interval=Nstep

  !limitsigma=.true.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !stop
  !call varscalc(NCELL,displs,vars)
  if(my_rank==0) then
    write(*,*) 'job number',number
  end if

  select case(problem)
  case('3dp','3dph')
    NCELLg=imax*jmax
  case('3dht')
    open(12,file=geofile,iostat=ios)
    if(ios /= 0) then
      if(my_rank==0)write(*,*) 'error: Failed to open geometry file'
      stop
    end if
    nl=0
    do
     read(12,'()',end=100)
     nl = nl + 1
    end do
100 close(12)
    ncellg=nl/7
    if(my_rank==0) write(*,*) 'NCELLg',ncellg
  end select

  if(ncellg==0.and.my_rank==0) then
    write(*,*) 'error: Ncellg is zero'
    stop
  end if
  !allocation
  allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg),ds(NCELLg))
  allocate(ag(NCELLg),bg(NCELLg),dcg(NCELLg),f0g(NCELLg))
  allocate(taug(NCELLg),sigmag(NCELLg),velG(NCELLg),rake(NCELLg))
  allocate(taudotg(NCELLg),sigdotg(NCELLg))

  xcol=0d0;ycol=0d0;zcol=0d0;ds=0d0

  select case(problem)
  case('2dp','2dvs')
    allocate(xel(NCELLg),xer(NCELLg))
    xel=0d0;xer=0d0
  case('2dn','2dph','2dn3','2dnh','25d')
    allocate(ang(NCELLg),xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg))
    ang=0d0;xel=0d0;xer=0d0;yel=0d0;yer=0d0
  case('3dp')
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
    xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
    zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
  case('3dnt','3dht')
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg))
    allocate(ys1(NCELLg),ys2(NCELLg),ys3(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg))
    allocate(ev11(NCELLg),ev12(NCELLg),ev13(NCELLg))
    allocate(ev21(NCELLg),ev22(NCELLg),ev23(NCELLg))
    allocate(ev31(NCELLg),ev32(NCELLg),ev33(NCELLg))
    xs1=0d0; xs2=0d0; xs3=0d0
    ys1=0d0; ys2=0d0; ys3=0d0
    zs1=0d0; zs2=0d0; zs3=0d0
    rake=0d0
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
    rake=0d0
  case('3dnr','3dhr','3dph')
    allocate(ang(NCELLg),angd(NCELLg))
    !xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
    !ys1=0d0; ys2=0d0; ys3=0d0; ys4=0d0
    !zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
    angd=0d0; ang=0d0; rake=0d0

  end select
  !allocate(vmax(NCELLg),vmaxin(NcELLg))

  !mesh generation (rectangular assumed)
  if(my_rank==0) write(*,*) 'Generating mesh'
  select case(problem)
  case('2dp','2dvs')
    call coordinate2dp(NCELLg,ds0,xel,xer,xcol)

  case('2dph')
    call coordinate2dph()
  
  case('2dn','2dnh','25d')
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
    !coordinate3dns(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)

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
    call evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33,ds)
  end select

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !stop

  rake=crake
  !nonuniform parameters from file
  if(parameterfromfile) then
    open(99,file=parameter_file)
    do i=1,ncellg
      read(99,*) rake(i),ag(i),bg(i),dcg(i),f0g(i),taug(i),sigmag(i),velg(i),taudotg(i),sigdotg(i)
    end do
    close(99)
  end if
  rake=rake/180d0*pi

  !HACApK setting
  lrtrn=HACApK_init(NCELLg,st_ctl,st_bemv,icomm)
  st_ctl%param(8)=20
  !if(latticeh) st_ctl%param(8)=20
  !lrtrn=HACApK_init(NCELLg,st_ctl2,st_bemv,icomm)
  allocate(coord(NCELLg,3))
  select case(problem)
  case('2dp','2dvs')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%problem=problem

  case('2dn','2dn3','2dph','2dnh','25d')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg),st_bemv%ds(NCELLg))
    allocate(st_bemv%ycol(NCELLg),st_bemv%yel(NCELLg),st_bemv%yer(NCELLg),st_bemv%ang(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%ycol=ycol;st_bemv%yel=yel;st_bemv%yer=yer
    st_bemv%ang=ang; st_bemv%ds=ds; st_bemv%w=fwid
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

  case('3dhr','3dnr','3dph')
    allocate(st_bemv%xcol(NCELLg),st_bemv%ycol(NCELLg),st_bemv%zcol(NCELLg))
    allocate(st_bemv%ang(NCELLg),st_bemv%angd(NCELLg),st_bemv%rake(NCELLg))
    st_bemv%xcol=xcol
    st_bemv%ycol=ycol
    st_bemv%zcol=zcol
    st_bemv%angd=angd
    st_bemv%ang=ang
    st_bemv%problem=problem
    st_bemv%rake=rake
    st_bemv%w=ds0
  case('3dht','3dnt')
    allocate(st_bemv%xcol(NCELLg),st_bemv%ycol(NCELLg),st_bemv%zcol(NCELLg))
    allocate(st_bemv%xs1(NCELLg),st_bemv%xs2(NCELLg),st_bemv%xs3(NCELLg))
    allocate(st_bemv%ys1(NCELLg),st_bemv%ys2(NCELLg),st_bemv%ys3(NCELLg))
    allocate(st_bemv%zs1(NCELLg),st_bemv%zs2(NCELLg),st_bemv%zs3(NCELLg))
    allocate(st_bemv%ev11(NCELLg),st_bemv%ev12(NCELLg),st_bemv%ev13(NCELLg))
    allocate(st_bemv%ev21(NCELLg),st_bemv%ev22(NCELLg),st_bemv%ev23(NCELLg))
    allocate(st_bemv%ev31(NCELLg),st_bemv%ev32(NCELLg),st_bemv%ev33(NCELLg))
    allocate(st_bemv%rake(NCELLg))
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
    st_bemv%rake=rake
  end select

  !generate kernel (H-matrix aprrox)
  if(my_rank==0) write(*,*) 'Generating kernel'
  do i=1,NCELLg
    coord(i,1)=xcol(i)
    coord(i,2)=ycol(i)
    coord(i,3)=zcol(i)
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
    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n,st_bemv,st_ctl2,coord,eps_h)
    lrtrn=HACApK_construct_LH(st_LHp_n,st_leafmtxp_n,st_bemv,st_ctl2,coord,eps_h)
  end if
  !write(*,*) my_rank,st_ctl%lpmd(33),st_ctl%lpmd(37)
  allocate(y(3*NCELL),yscal(3*NCELL),dydx(3*NCELL))
  allocate(psi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL),idisp(NCELL),cslip(NCELL))
  psi=0d0;vel=0d0;tau=0d0;sigma=0d0;disp=0d0
  allocate(a(NCELL),b(NCELL),dc(NCELL),f0(NCELL),taudot(NCELL),tauddot(NCELL),sigdot(NCELL))
  taudot=0d0;sigdot=0d0

  !uniform frictional parameters
  a=a0
  b=b0
  dc=dc0
  f0=mu0

  if(.not.backslip) then
    taudot=sr
    sigdot=0d0
  end if


  if(parameterfromfile) then
    do i=1,NCELL
      i_=st_sum%lodc(i)
      a(i)=ag(i_)
      b(i)=bg(i_)
      dc(i)=dcg(i_)
      f0(i)=f0g(i_)
      taudot(i)=taudotg(i_)
      sigdot(i)=sigdotg(i_)
    end do
  end if

  if(deepcreep) then
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
  if(forward) call forward_check()
  if(inverse) call inverse_problem()

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !restart
  if(restart) then
    if(my_rank==0)then
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname,status='old')
      do while(.true.)
        read(52,*,iostat=ios)k,dummy(1:9)
        if(mod(k,interval)==0) then
          kstart=k+1
          x=dummy(1)
          dtnxt=min(2*dummy(8),0.9*dummy(8)*(dummy(7)**(-0.2)))
        end if
        if(ios<0) exit
      end do
      close(52)
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname,status='old',position='append')

      kend=k
      write(*,*) kstart,x,dtnxt

      write(fname,'("output/event",i0,".dat")') number
      open(44,file=fname,position='append')
      open(19,file='job.log',position='append')
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
      write(19,'(a20,i0,a6,a12,a6,a12)') 'Starting job number=',number,'date',sys_time(1),'time',sys_time(2)
      close(19)

    end if
    call MPI_bcast(kstart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_bcast(dtnxt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    if(my_rank<npd) then
      nout=my_rank+100
      write(fname,'("output/vel",i0,"_",i0,".dat")') number,my_rank
      open(nout,file=fname,form='unformatted',access='stream')
      inquire(nout,size=file_size)
      m=file_size/8
      allocate(rdata(m))
      read(nout) rdata
      vel(1:NCELL)=rdata(m-NCELL+1:m)
      !write(*,*) my_rank,vel
      close(nout)

      write(fname,'("output/slip",i0,"_",i0,".dat")') number,my_rank
      open(nout,file=fname,form='unformatted',access='stream')
      read(nout) rdata
      disp(1:NCELL)=rdata(m-NCELL+1:m)
      close(nout)

      write(fname,'("output/sigma",i0,"_",i0,".dat")') number,my_rank
      open(nout,file=fname,form='unformatted',access='stream')
      read(nout) rdata
      sigma(1:NCELL)=rdata(m-NCELL+1:m)
      close(nout)

      write(fname,'("output/tau",i0,"_",i0,".dat")') number,my_rank
      open(nout,file=fname,form='unformatted',access='stream')
      read(nout) rdata
      tau(1:NCELL)=rdata(m-NCELL+1:m)
      close(nout)
      !write(*,*) my_rank,m

      psi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
      !write(*,*) tau
      write(fname,'("output/vel",i0,"_",i0,".dat")') number,my_rank
      open(nout,file=fname,form='unformatted',access='stream',status='old',position='append')
      nout2=nout+np
      write(fname,'("output/slip",i0,"_",i0,".dat")') number,my_rank
      open(nout2,file=fname,form='unformatted',access='stream',status='old',position='append')
      nout3=nout2+np
      write(fname,'("output/sigma",i0,"_",i0,".dat")') number,my_rank
      open(nout3,file=fname,form='unformatted',access='stream',status='old',position='append')
      nout4=nout3+np
      write(fname,'("output/tau",i0,"_",i0,".dat")') number,my_rank
      open(nout4,file=fname,form='unformatted',access='stream',status='old',position='append')

      nout5=nout4+np
      write(fname,'("output/EQslip",i0,"_",i0,".dat")') number,my_rank
      open(nout5,file=fname,form='unformatted',access='stream',status='replace')

      rankloc=-1
    if(my_rank<npd) then
      do i=1,ncell
        if(st_sum%lodc(i)==loc) then
          rankloc=my_rank
          loc_=i
          write(*,*) rankloc,loc_
          write(*,*)  xcol(loc),ycol(loc),zcol(loc)
          exit
        end if
      end do
    end if

    if(my_rank==rankloc) then
      write(fname,'("output/local",i0,"-",i0,".dat")') number,loc
      open(53,file=fname)
      write(*,*) my_rank,loc
    end if

    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_bcast(psi,NCELL,MPI_REAL8,0,st_ctl%lpmd(35),ierr)
    call MPI_bcast(vel,NCELL,MPI_REAL8,0,st_ctl%lpmd(35),ierr)
    call MPI_bcast(tau,NCELL,MPI_REAL8,0,st_ctl%lpmd(35),ierr)
    call MPI_bcast(sigma,NCELL,MPI_REAL8,0,st_ctl%lpmd(35),ierr)
    call MPI_bcast(disp,NCELL,MPI_REAL8,0,st_ctl%lpmd(35),ierr)

    s=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      s=s+ds(i_)
    end do

    call MPI_reduce(s,sG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(mvelG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr); time2=MPI_Wtime()
    if(my_rank==0) write(*,*) 'Finished all initial processing, time(s)=',time2-time1
    time1=MPI_Wtime()
  
  !no restart
  else 

    !setting initial condition

    !uniform values
    sigma=sigmainit
    tau=tauinit
    !if(my_rank==0) write(*,*) tau

    if(muinit.ne.0d0) tau=sigma*muinit
    vel=tau/abs(tau)*velinit
    mu=tau/sigma
    psi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
    disp=0d0

    !non-uniform initial stress from file
    if(parameterfromfile) then
      do i=1,NCELL
        i_=st_sum%lodc(i)
        tau(i)=taug(i_)
        sigma(i)=sigmag(i_)
        vel(i)=velg(i_)
        mu(i)=tau(i)/sigma(i)
        psi(i)=a(i)*dlog(2*vref/vel(i)*sinh(tau(i)/sigma(i)/a(i)))
        !write(*,*) tau(i),sigma(i),vel(i)
      end do
    end if

    if(bgstress) then
      call initcond()
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !calculate Lb/ds
    if(my_rank==0) write(*,*) 'Lb/ds~',rigid*dc(1)/b(1)/sigma(1) 

    x=0d0
    kstart=1
    kend=0
    call output_coord()
    if(my_rank==0) then
      ! write(fname,'("output/ind",i0,"_",i0,".dat")') number,my_rank
      ! nout=my_rank+100
      ! open(nout,file=fname)
      ! !write(nout)st_sum%lodc(1:NCELL)
      ! do i=1,ncell
      !   write(nout,'(i0)') st_sum%lodc(i)
      ! end do
      ! close(nout)

      ! write(fname,'("output/xyz",i0,".dat")') number
      ! nout=my_rank+100
      ! open(nout,file=fname)
      ! do i=1,ncellg
      !   write(nout,'(3e15.6)') xcol(i),ycol(i),zcol(i)
      ! end do
      ! close(nout)
      ! ! write(fname,'("output/prm",i0,"_",i0,".dat")') number,my_rank
      ! ! nout=my_rank+100
      ! ! open(nout,file=fname)
      ! ! do i=1,ncell
      ! !   write(nout,*) a(i),taudot(i),sigdot(i_)
      ! ! end do
      ! ! close(nout)

      nout=80
      write(fname,'("output/vel",i0,".dat")') number
      open(nout,file=fname,form='unformatted',access='stream',status='replace')
      nout2=nout+1
      write(fname,'("output/slip",i0,".dat")') number
      open(nout2,file=fname,form='unformatted',access='stream',status='replace')
      nout3=nout2+1
      write(fname,'("output/sigma",i0,".dat")') number
      open(nout3,file=fname,form='unformatted',access='stream',status='replace')
      nout4=nout3+1
      write(fname,'("output/tau",i0,".dat")') number
      open(nout4,file=fname,form='unformatted',access='stream',status='replace')
      
      nout5=nout4+1
      write(fname,'("output/EQslip",i0,".dat")') number
      open(nout5,file=fname,form='unformatted',access='stream',status='replace')

      write(fname,'("output/time",i0,".dat")') number
      open(50,file=fname)
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname)
      write(fname,'("output/event",i0,".dat")') number
      open(44,file=fname)
      open(19,file='job.log',position='append')
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
      write(19,'(a20,i0,a6,a12,a6,a12,a4,i0)') 'Starting job number=',number,'date',sys_time(1),'time',sys_time(2),'np',np
      close(19)
    end if


    rankloc=-1
    if(my_rank<npd) then
      do i=1,ncell
        if(st_sum%lodc(i)==loc) then
          rankloc=my_rank
          loc_=i
          write(*,*) rankloc,loc_
          write(*,*)  xcol(loc),ycol(loc),zcol(loc)
          exit
        end if
      end do
    end if

    if(my_rank==rankloc) then
      write(fname,'("output/local",i0,"-",i0,".dat")') number,loc
      open(53,file=fname)
      write(*,*) my_rank,loc
    end if

    s=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      s=s+ds(i_)
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

    ! meandisp=sum(disp*ds)/sum(ds)
    meandisp=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meandisp=meandisp+disp(i)*ds(i_)
    end do
    !call MPI_ALLREDUCE(meandisp,meandispG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meandisp,meandispG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meandispG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meandispG=meandispG/sg

    ! meanmu=sum(mu*ds)
    meanmu=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meanmu=meanmu+mu(i)*ds(i_)
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
      write(*,*) 'initial Vmax',mvelG
      write(*,*) 'initial mu_avg',meanmuG
      call output_monitor()
    end if
    call output_field(mvel_loc)
    dtnxt = dtinit
  end if
  tout=dtout*365*24*60*60
  !tout=20*365*24*60*60
  rk=0

  !outv=1d-6
  slipping=.false.
  eventcount=0
  sw=0
  timer=0d0
  timeH=0d0
  !time2=MPI_Wtime()
  !output initial values


  !do i=1,NCELLg
  !  write(50,'(8e15.6,i6)') xcol(i),ycol(i),vel(i),tau(i),sigma(i),mu(i),disp(i),x,k
  !end do
  !write(50,*)

  !$omp parallel do
  do i=1,NCELL
    y(3*i-2) = psi(i)
    y(3*i-1) = tau(i)
    y(3*i)=sigma(i)
    !if(my_rank==53)write(*,*) psi(i),tau(i),sigma(i)
  end do
  !$omp end parallel do
  !stop


  do k=kstart,NSTEP
    !parallel computing for Runge-Kutta
    dttry = dtnxt
    !time3=MPI_Wtime()
    call rkqs(y,dydx,x,dttry,eps_r,dtdid,dtnxt,errmax_gb,nrjct)
    !call rkqs2(y,dydx,x,dttry,eps_r,errold,dtdid,dtnxt,errmax_gb,nrjct)
    !time4=MPI_Wtime()
    !timer=timer+time4-time3

    !if normal stress is bounded
      do i=1,NCELL
        if(y(3*i)<minsig) y(3*i)=minsig
        if(y(3*i)>maxsig) y(3*i)=maxsig
      end do

    !compute physical values for control and output
    !$omp parallel do
    do i = 1, NCELL
      psi(i) = y(3*i-2)
      tau(i) = y(3*i-1)
      sigma(i)=y(3*i)
      disp(i)=disp(i)+vel(i)*dtdid*0.5d0 !2nd order
      vel(i)= 2*vref*exp(-psi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
      disp(i)=disp(i)+vel(i)*dtdid*0.5d0
      mu(i)=tau(i)/sigma(i)
    end do
    !$omp end parallel do

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

    ! meandisp=sum(disp*ds)/sum(ds)
    meandisp=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meandisp=meandisp+disp(i)*ds(i_)
    end do
    !call MPI_ALLREDUCE(meandisp,meandispG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meandisp,meandispG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meandispG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meandispG=meandispG/sg

    ! meanmu=sum(mu*ds)
    meanmu=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      meanmu=meanmu+mu(i)*ds(i_)
    end do
    !call MPI_ALLREDUCE(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_reduce(meanmu,meanmuG,1,MPI_REAL8,MPI_SUM,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(meanmuG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    meanmuG=meanmuG/sg
    !if(outfield.and.(my_rank.lt.npd)) call output_field()
 
    !output distribution control
    outfield=.false.
    if(mod(k,interval)==0) then
      outfield=.true.
    end if
    if(x>tout) then
      outfield=.true.
      tout=x+dtout*365*24*60*60
      if(x<300*356*23*3600) outfield=.false.
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
      !lattice H
       call output_field(mvel_loc)
    end if

    if(my_rank==0.and.k>kend) then
      call output_monitor()
    end if
    if(my_rank==rankloc.and.k>kend) then
      call output_local(loc_)
    end if
    time4=MPI_Wtime()
    timer=timer+time4-time3

    !event list
    if(.not.slipping) then
      if(mvelG>1d-1) then
        slipping=.true.
        eventcount=eventcount+1
        moment0=meandispG
        idisp=disp
        hypoloc=mvel_loc(1)
        onset_time=x
        !tout=onset_time

        !onset save
        ! if(slipevery.and.(my_rank<npd)) then
        !   write(nout) vel
        !   write(nout2) disp
        !   write(nout3) sigma
        !   write(nout4) tau
        ! end if

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
          write(44,'(i0,i7,f17.2,f14.4,i8)') eventcount,k,onset_time,(log10(moment*rigid*sg)+5.9)/1.5,hypoloc
        end if
        cslip=disp-idisp
        if(my_rank<npd) then
          write(nout5) cslip
        end if
        ! if(slipevery.and.(my_rank<npd)) then
        !   write(nout) vel
        !   write(nout2) disp
        !   write(nout3) sigma
        !   write(nout4) tau
        ! end if

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
    write(52,'(i7,f19.4,7e16.5,i4,f16.4)')k,x,log10(mvelG),meandispG,meanmuG,maxnormG,minnormG,errmax_gb,dtdid,nrjct,time2-time1
  end subroutine

  subroutine output_local(loc_)
    implicit none
    integer,intent(in)::loc_
    write(53,'(i7,f19.4,6e16.5)')k,x,log10(vel(loc_)),disp(loc_),tau(loc_),sigma(loc_),mu(loc_),psi(loc_)
  end subroutine

  subroutine output_field(mvel_loc)
    implicit none
    integer::nn,rcounts(npd),displs(npd+1)
    integer::mvel_loc(1)
    real(8)::velG(NCELLg),tauG(NCELLg),sigmaG(Ncellg),dispG(ncellg)
    call MPI_GATHER(ncell,1,MPI_INT,rcounts,1,MPI_INT,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

    displs(1)=0
    do nn=2,npd+1
      displs(nn)=displs(nn-1)+rcounts(nn-1)
    end do

    call MPI_GATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_GATHERv(tau,NCELL,MPI_REAL8,tauG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_GATHERv(sigma,NCELL,MPI_REAL8,sigmaG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_GATHERv(disp,NCELL,MPI_REAL8,dispG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

    if(my_rank==0) then
      mvel_loc=maxloc(abs(velG))
      write(nout) velG
      write(nout2) dispG
      write(nout3) sigmaG
      write(nout4) tauG
    end if

  end subroutine

  ! subroutine output_EQslip()
  !   implicit none
  !   integer::nn,rcounts(npd),displs(npd+1)
  !   real(8)::EQslipG(ncellg)
  !   nout=100

  !   call MPI_GATHER(ncell,1,MPI_INT,rcounts,1,MPI_INT,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

  !   displs(1)=0
  !   do nn=2,npd+1
  !     displs(nn)=displs(nn-1)+rcounts(nn-1)
  !   end do

  !   call MPI_GATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

  !   if(my_rank==0) then
  !     write(nout) velG
  !   end if

  ! end subroutine
  ! subroutine output_field_init()
  !   implicit none
  !   integer::nn,rcounts(npd),displs(npd+1)
  !   !real(8),allocatable::vel(:)
  !   real(8)::velG(NCELLg)
  !   nout=100
  !   if(my_rank==0) then
  !      write(fname,'("output/vel",i0,".dat")') number
  !      open(nout,file=fname,form='unformatted',access='stream',status='replace')
  !   end if

  !   write(*,*) my_rank,st_ctl%lpmd(37),st_ctl%lpmd(31)

  !   call MPI_GATHER(ncell,1,MPI_INT,rcounts,1,MPI_INT,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

  !   displs(1)=0
  !   do nn=2,npd+1
  !     displs(nn)=displs(nn-1)+rcounts(nn-1)
  !   end do

  !   call MPI_GATHERv(vel,NCELL,MPI_REAL8,velG,rcounts,displs,MPI_REAL8,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)

  !   if(my_rank==0) write(nout) velG

  ! end subroutine


  subroutine output_coord()
    implicit none
    integer::nn
    nout=100
    nout2=101

   if(my_rank==0) then
    write(fname,'("output/xyz",i0,".dat")') number
    open(nout,file=fname,status='replace')
    write(fname,'("output/ind",i0,".dat")') number
    open(nout2,file=fname,status='replace')
      do i=1,ncell
        i_=st_sum%lodc(i)
        write(nout2,*) i_
        write(nout,'(3e15.6)') xcol(i_),ycol(i_),zcol(i_)
      end do
    close(nout)
    close(nout2)

    end if
    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    do nn=1,npd-1
      if(my_rank==nn) then
      write(fname,'("output/xyz",i0,".dat")') number
      open(nout,file=fname,position='append')
      write(fname,'("output/ind",i0,".dat")') number
      open(nout2,file=fname,position='append')
      do i=1,ncell
        i_=st_sum%lodc(i)
        write(nout2,*) i_
        write(nout,'(3e15.6)') xcol(i_),ycol(i_),zcol(i_)
      end do
      close(nout)
      close(nout2)
      end if
      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end do

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
    !planar fault with element size ds and dipangle=dipangle
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
      real(8)::xc,yc,zc,amp,yr(0:jmax),zr(0:jmax),stangle


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
             !ys1(k)=yr(j-1)
             !ys2(k)=yr(j-1)
             !ys3(k)=yr(j)
             !ys4(k)=yr(j)
             ycol(k)=0.5d0*(yr(j-1)+yr(j)) !-(j-0.5d0)*ds0*cos(dipangle)
             zcol(k)=0.5d0*(zr(j-1)+zr(j)) !-(j-0.5d0)*ds0*sin(dipangle)
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
             angd(k)=datan2(zr(j-1)-zr(j),yr(j-1)-yr(j))
             ang(k)=0d0
             !write(*,*)angd(k)
             if(my_rank==0)write(111,*)xcol(k),ycol(k),zcol(k)
             end do
         end do

      return
    end subroutine coordinate3ddip

    subroutine coordinate3dns(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
      implicit none
      integer,intent(in)::NCELLg
      real(8),intent(out)::xcol(:),ycol(:),zcol(:)
      real(8),intent(out)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
      integer::i,j,k,imax,jmax
      real(8)::dipangle,xc,yc,zc,amp
      real(4)::xl(0:2048,0:2048)

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
      ! if(my_rank.eq.0) open(32,file='tmp')
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
      !   if(my_rank.eq.0) write(32,*) xcol(i),ycol(i),zcol(i)
      ! end do
  
      return
    end subroutine coordinate3dns

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

  subroutine taudot_3dph()
    real(8)::rate
    !open(111,file='tmp')
    do i=1,NCELL
      i_=st_sum%lodc(i)
      taudot(i)=okada_load(xcol(i_),ycol(i_),zcol(i_),-ds0*imax/2,ds0*imax/2,ds0*jmax,0d0,pi/2.0,0d0)*vpl
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
      taudot(i)=load2dnh(xcol(i_),ycol(i_),xer(Nload),yer(Nload),pi*dipangle/180,vv)*vpl
      vv="n"
      sigdot(i)=load2dnh(xcol(i_),ycol(i_),xer(Nload),yer(Nload),pi*dipangle/180,vv)*vpl

      taudot(i)=taudot(i)+sr*0.5*sin(2*ang(i_))
      sigdot(i)=sigdot(i)+sr*sin(ang(i_))**2

      write(*,*) xcol(i_),taudot(i),sigdot(i)
    end do
  end subroutine

  subroutine initcond()
    implicit none
    real(8),parameter::rho=1.6e3,g=9.80
    real(8)::sxx0,sxy0,syy0
  do i=1,ncell
    i_=st_sum%lodc(i)
    syy0=min(50.0,rho*g*ycol(i_)/1e3)
    sxy0=0.0
    sxx0=syy0*3.0
    tau(i)=sxy0*cos(2*ang(i_))+0.5*(sxx0-syy0)*sin(2*ang(i_))
    sigma(i)=sin(ang(i_))**2*sxx0+cos(ang(i_))**2*syy0+sxy0*sin(2*ang(i_))
  end do
    mu=tau/sigma
    vel=velinit
    psi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
    disp=0d0

  end subroutine

  !computing dydx for time integration
  subroutine derivs(x, y, dydx)
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    !include 'mpif.h'
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    !integer,intent(in) :: NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in) :: x
    real(8),intent(in) ::y(:)
    real(8),intent(out) :: dydx(:)
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),psitmp(NCELL)
    real(8) :: sum_gs(NCELL),sum_gn(NCELL)!,velstmpG(NCELLg),veldtmpG(NCELLg)
    !real(8) :: sum_xx(NCELL),sum_xy(NCELL),sum_yy(NCELL)!,sum_xz(NCELL),sum_yz(NCELL),sum_zz(NCELL)
    !real(8) :: sum_xxG(NCELLg),sum_xyG(NCELLg),sum_yyG(NCELLg)!,sum_xzG(NCELLg),sum_yzG(NCELLg),sum_zzG(NCELLg)
    !real(8) :: sum_xx2G(NCELLg),sum_xy2G(NCELLg),sum_yy2G(NCELLg),sum_xz2G(NCELLg),sum_yz2G(NCELLg),sum_zz2G(NCELLg)
    !real(8) :: veltmpG(NCELLg),sum_gsg(NCELLg),sum_gng(NCELLg),sum_gdg(NCELLg)!,efftmpG(NCELLg)
    !real(8) :: sum_gs2G(NCELLg),sum_gd2G(NCELLg),sum_gn2G(NCELLg)
    real(8) :: time3,time4,c1, c2, c3, arg,arg2,c,g,tauss,Arot(3,3),p(6),fac,sxx0,sxy0,syy0
    integer :: i, j, nc,ierr,lrtrn,i_

    !if(latticeh) then

    !$omp parallel do
    do i = 1, NCELL
      psitmp(i) = y(3*i-2)
      tautmp(i) = y(3*i-1)
      sigmatmp(i) = y(3*i)
      veltmp(i) = 2*vref*dexp(-psitmp(i)/a(i))*dsinh(tautmp(i)/sigmatmp(i)/a(i))
      !if(my_rank==0)write(*,*) veltmp(i)
    enddo
    !$omp end parallel do

    !matrix-vector mutiplation
    if(backslip) then
      st_vel%vs=veltmp-vpl
    else
      st_vel%vs=veltmp
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
      call HACAPK_adot_lattice_hyp(st_sum,st_LHP_n,st_ctl,wws,st_vel)
      if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') then
        sum_gn(:)=st_sum%vs(:)/ds0
      else
        sum_gn(:)=st_sum%vs(:)
      end if
    end if
    !time4=MPI_Wtime()
    !timeH=timeH+time4-time3

    !$omp parallel do
    do i=1,NCELL
      sum_gs(i)=sum_gs(i)+taudot(i)
      sum_gn(i)=sum_gn(i)+sigdot(i)
      !write(*,*) sum_gs(i)
      call deriv_d(sum_gs(i),sum_gn(i),psitmp(i),tautmp(i),sigmatmp(i),veltmp(i),a(i),b(i),dc(i),f0(i),dydx(3*i-2),dydx(3*i-1),dydx(3*i))
    enddo
    !$omp end parallel do
    !write(*,*) dydx
    return
  end subroutine

  subroutine deriv_d(sum_gs,sum_gn,psitmp,tautmp,sigmatmp,veltmp,a,b,dc,f0,dpsidt,dtaudt,dsigdt)
    implicit none
    real(8)::fss,dvdtau,dvdsig,dvdpsi,mu,psiss,dcv
    real(8),parameter::fw=0.2,vw=0.1
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs,sum_gn,psitmp,tautmp,sigmatmp,veltmp,a,b,dc,f0
    real(8),intent(out)::dpsidt,dtaudt,dsigdt
    dsigdt=sum_gn
    !regularized slip law
    !fss=mu0+(a(i_)-b(i_))*dlog(abs(veltmp(i))/vref)
    !fss=fw(i_)+(fss-fw(i_))/(1.d0+(veltmp(i)/vw(i_))**8)**0.125d0 !flash heating
    !dpsidt(i)=-abs(veltmp(i))/dc(i_)*(abs(tautmp(i))/sigmatmp(i)-fss)

    select case(evlaw)
    case('aging')
      dpsidt=b/dc*vref*dexp((f0-psitmp)/b)-b*abs(veltmp)/dc    !slip law
    case('slip')
      fss=f0+(a-b)*dlog(abs(veltmp)/vref)
      dpsidt=-abs(veltmp)/dc*(abs(tautmp)/sigmatmp-fss)
    case('dynamic')
      fss=f0+(a-b)*dlog(abs(veltmp)/vref)
      fss=fw+(fss-fw)/(1.d0+(veltmp/vw)**8)**0.125d0 !flash heating
      psiss=a*dlog(2*vref/veltmp*sinh(fss/a))
      !dcv=dc*(1+(veltmp/vref)**2)**(0.5)
      dcv=dc*(1+log((veltmp/vref)+1))
      dpsidt=-veltmp/dcv*(psitmp-psiss)
    end select

    !regularized aing law
    !dpsidt=b/dc*vref*dexp((f0-psitmp)/b)-b*abs(veltmp)/dc

    !regularized aging law with cutoff velocity for evolution
    !dpsidt(i)=b(i_)/dc(i_)*vref*dexp((f0(i_)-psitmp(i))/b(i_))*(1d0-abs(veltmp(i))/vref*(exp((psitmp(i)-f0(i_))/b(i_))-vref/vc(i_)))

    dvdtau=2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)/(a*sigmatmp)
    !mu=tautmp(i)/sigmatmp(i)
    !dvdtau=vref*(exp((mu-psitmp(i))/a(i))+exp((-mu-psitmp(i))/a(i)))/(a(i)*sigmatmp(i))
    dvdsig=-2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)*tautmp/(a*sigmatmp**2)
    !dvdpsi=2*vref*exp(-psitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
    dvdpsi=-veltmp/a
    !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdpsi*psitmp(i)*dvdsig*sigmatmp(i))
    dtaudt=sum_gs-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt)
    dtaudt=dtaudt/(1d0+0.5d0*rigid/vs*dvdtau)
    !write(*,*) rigid/vs*dvdtau
    if(veltmp<=0d0) then
      dvdtau=2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)/(a*sigmatmp)
      dvdsig=-2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)*tautmp/(a*sigmatmp**2)
      !sign ok?
      !dvdpsi=2*vref*exp(-psitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
      dvdpsi=-veltmp/a
      !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(-dvdpsi*psitmp(i)*dvdsig*sigmatmp(i))
      dtaudt=sum_gs-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt)
      dtaudt=dtaudt/(1d0+0.5d0*rigid/vs*dvdtau)
    end if
  end subroutine
  subroutine deriv_3dn(sum_gs,sum_gd,sum_gn,psitmp,taustmp,taudtmp,tautmp,sigmatmp,veltmp,a,b,dc,f0,dpsidt,dtausdt,dtauddt,dsigdt)
    implicit none
    integer::i
    real(8)::fss,dvdtau,dvdsig,dvdpsi,absV
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs,sum_gd,sum_gn,psitmp,taustmp,taudtmp,tautmp,sigmatmp,veltmp,a,b,dc,f0
    real(8),intent(out)::dpsidt,dtausdt,dtauddt,dsigdt
    !write(*,*) 'vel',veltmp(i)
    dsigdt=sum_gn
    !fss=mu0+(a(i)-b(i))*dlog(abs(veltmp(i))/vref)
    !fss=fw(i)+(fss-fw(i))/(1.d0+(veltmp(i)/vw(i))**8)**0.125d0 !flash heating
    !slip law
    !dpsidt(i)=-abs(veltmp(i))/dc(i)*(abs(tautmp(i))/sigmatmp(i)-fss)
    !aing law
    dpsidt=b*vref/dc*exp((f0-psitmp)/b)-b*veltmp/dc
    dvdtau=2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)/(a*sigmatmp)
    dvdsig=-2*vref*dexp(-psitmp/a)*dcosh(tautmp/sigmatmp/a)*tautmp/(a*sigmatmp**2)
    dvdpsi=-veltmp/a
    dtausdt=sum_gs-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt)*(taustmp/tautmp)
    dtausdt=dtausdt/(1d0+0.5d0*rigid/vs*dvdtau)
    dtauddt=sum_gd-0.5d0*rigid/vs*(dvdpsi*dpsidt+dvdsig*dsigdt)*(taudtmp/tautmp)
    dtauddt=dtauddt/(1d0+0.5d0*rigid/vs*dvdtau)
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
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4,kp=0.08

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

      do i=1,3*NCELL
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
        if(my_rank.eq.0)write(*,*)'error: dt is too small'
        stop
      end if

    end do

    hnext=min(2*h,SAFETY*h*(errmax_gb**PGROW))
    ! if(outpertime) then
    !   hnext=min(hnext,dtout*365*24*3600)
    ! end if
    !if(load==0)hnext=min(hnext,dtmax)
    !hnext=max(0.249d0*ds0/vs,SAFETY*h*(errmax_gb**PGROW))

    !hnext=min(,1d9)

    hdid=h
    x=x+h
    y(:)=ytemp(:)

  end subroutine
  !---------------------------------------------------------------------
  subroutine rkqs2(y,dydx,x,htry,eps,errold,hdid,hnext,errmax_gb,nrjct)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    !include 'mpif.h'
    !integer::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::htry,eps,errold
    real(8),intent(inout)::y(:),x,dydx(:)
    real(8),intent(out)::hdid,hnext,errmax_gb !hdid: resulatant dt hnext: htry for the next
    integer,intent(out)::nrjct
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
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

      do i=1,3*NCELL
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
        if(my_rank.eq.0)write(*,*)'error: dt is too small'
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
    !include 'mpif.h'
    !integer,intent(in)::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::y(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    !integer,intent(out)::ierr
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
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
    if(my_rank==0) write(*,*) 'Calculating stress field from uniform slip'
    ret1=0d0
    ret2=0d0

    vec=1d0
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
        call HACApK_adot_lattice_hyp(st_sum,st_LHp_n,st_ctl,wws,st_vel)
        ret2(:)=st_sum%vs(:)
        if(problem=='3dnr'.or.problem=='3dhr'.or.problem=='3dph') ret2(:)=ret2(:)/ds0
      end if


      !lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n,st_bemv,st_ctl,ret2,vec)
      if(my_rank==0) then
        do i=1,NCELL
          i_=st_sum%lodc(i)
          write(29,'(6e16.4)') xcol(i_),ycol(i_),zcol(i_),vec(i_),ret1(i),ret2(i)
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
          write(29,'(2e16.4)') xcol(i),sigdot(i)
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
          write(29,'(4e16.4)') xcol(i),ycol(i),zcol(i),sigdot(i)
        end do
      end if
    end select
    Call MPI_FINALIZE(ierr)
    stop
  end subroutine
end program
