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
  real(8),allocatable::phi(:),vel(:),tau(:),sigma(:),disp(:),mu(:),rupt(:),idisp(:),velp(:),pfhyd(:)
  real(8),allocatable::taus(:),taud(:),vels(:),veld(:),disps(:),dispd(:),rake(:),pf(:),sigmae(:),ks(:),qflow(:)

  real(8),allocatable::rdata(:)
  integer::lp,i,i_,j,k,kstart,kend,m,counts,interval,lrtrn,nl,ios,nmain,rk,nout(10),file_size
  integer,allocatable::locid(:)
  integer::hypoloc(1),load,eventcount,thec,inloc,sw,errmaxloc,niter,mvelloc(1)

  !controls
  logical::aftershock,buffer,nuclei,slipping,outfield,slipevery,limitsigma,dcscale,slowslip,slipfinal,outpertime
  logical::initcondfromfile,parameterfromfile,backslip,sigmaconst,foward,inverse,geofromfile,restart,latticeh,pressuredependent
  character*128::fname,dum,law,input_file,problem,geofile,param,pvalue,slipmode,project,parameter_file,outdir,command,bcl,bcr,evlaw
  real(8)::a0,a1,b0,dc0,sr,omega,theta,dtau,tiny,moment,wid,normal,ieta,meanmu,meanmuG,meandisp,meandispG,moment0,mvel,mvelG
  real(8)::psi,vc0,mu0,onset_time,tr,vw0,fw0,velmin,tauinit,intau,trelax,maxnorm,maxnormG,minnorm,minnormG,sigmainit,muinit
  real(8)::r,vpl,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxeventi,sparam,tmax,dtmax,tout,dummy(10)
  real(8)::alpha,ds0,amp,mui,velinit,phinit,velmax,maxsig,minsig,v1,dipangle,crake,s,sg,cdiff,q0,qin
  real(8)::kpmax,kpmin,kp0,kT,kL,s0,ksinit,dtout,pfinit,pbcl,pbcr,lf,eta,beta,nr,str,cc
  real(8)::b0list(9)=(/0.012,0.0115,0.011,0.0105,0.010,0.0095,0.009,0.0085,0.008/)

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

  open(33,file=input_file,status='old',iostat=ios)
  !if(my_rank==0) write(*,*) 'input_file',input_file
  if(ios /= 0) then
    if(my_rank==0)write(*,*) 'error: Failed to open input file'
    stop
  end if

  !get filenumber
  number=0
  if(input_file(1:11)=='examples/in') then
    input_file=adjustl(input_file(12:))
    write(*,*) input_file
    read(input_file,*) number
    write(*,*) number
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
  tmax=1d12
  dipangle=0d0
  nuclei=.false.
  slipevery=.false.
  sigmaconst=.false.
  foward=.false.
  inverse=.false.
  slipfinal=.false.
  restart=.false.
  latticeh=.false.
  outpertime=.false.
  maxsig=300d0
  minsig=10d0
  amp=0d0
  dtinit=1d0
  tp=86400d0
  trelax=1d18
  project="none"
  initcondfromfile=.false.
  parameterfromfile=.false.
  pressuredependent=.true.
  bcl='Neumann'
  bcr='Neumann'
  evlaw='aging'
  alpha=1.0
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
    case('cdiff')
      read (pvalue,*) cdiff
    case('q0')
      read (pvalue,*) q0
    case('eta')
      read (pvalue,*) eta
    case('nr')
      read (pvalue,*) nr
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
    case('alpha')
      read (pvalue,*) alpha
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
    case('bcl')
      read (pvalue,'(a)') bcl
    case('bcr')
      read (pvalue,'(a)') bcr
    case('evlaw')
      read (pvalue,'(a)') evlaw
    end select
  end do
  close(33)
  str=nr*beta
  cc=nr*beta*eta

  !limitsigma=.true.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  select case(problem)
  case('3dp','3dph')
    NCELLg=imax*jmax
  end select
  !stop
  !call varscalc(NCELL,displs,vars)
  if(my_rank==0) then
    write(*,*) 'job number',number
    write(*,*) 'project',project
  end if

  if(ncellg==0.and.my_rank==0) then
    write(*,*) 'error: Ncell is zero'
    stop
  end if
  !allocation
  allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg),ds(NCELLg))
  allocate(ag(NCELLg),bg(NCELLg),dcg(NCELLg),f0g(NCELLg),qin(ncellg))
  allocate(taug(NCELLg),sigmag(NCELLg),velG(NCELLg),rake(NCELLg),pfg(NCELLg))
  allocate(taudotg(NCELLg),sigdotg(NCELLg))
  qin=0d0

  nl=5
  allocate(locid(nl))

  xcol=0d0;ycol=0d0;zcol=0d0;ds=0d0

  select case(problem)
  case('2dp','2dpah')
    allocate(xel(NCELLg),xer(NCELLg))
    xel=0d0;xer=0d0
  case('2dn','2dph','2dn3')
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
      read(99,*) rake(i),ag(i),bg(i),dcg(i),f0g(i),taug(i),sigmag(i),pfg(i),velg(i),taudotg(i),sigdotg(i)
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
  case('2dp','2dpah')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%problem=problem

  case('2dn','2dn3','2dph')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg),st_bemv%ds(NCELLg))
    allocate(st_bemv%ycol(NCELLg),st_bemv%yel(NCELLg),st_bemv%yer(NCELLg),st_bemv%ang(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%ycol=ycol;st_bemv%yel=yel;st_bemv%yer=yer
    st_bemv%ang=ang; st_bemv%ds=ds
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
  case('2dp','2dpah','2dn3','3dp')
    sigmaconst=.true.
  end select
  st_bemv%md='s'
  st_bemv%v='s'
  lrtrn=HACApK_generate(st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
  if(.not.sigmaconst) then
    st_bemv%v='n'
    lrtrn=HACApK_generate(st_leafmtxp_n,st_bemv,st_ctl,coord,eps_h)
  end if
  !if(latticeh) then
  lrtrn=HACApK_construct_LH(st_LHp_s,st_leafmtxp_s,st_bemv,st_ctl,coord,eps_h)
  allocate(wws(st_leafmtxp_s%ndlfs))
  lrtrn=HACApK_gen_lattice_vector(st_vel,st_leafmtxp_s,st_ctl)
  lrtrn=HACApK_gen_lattice_vector(st_sum,st_leafmtxp_s,st_ctl)
  NCELL=st_vel%ndc
  if(.not.sigmaconst) lrtrn=HACApK_construct_LH(st_LHp_n,st_leafmtxp_n,st_bemv,st_ctl,coord,eps_h)
  !end if
  !write(*,*) my_rank,st_ctl%lpmd(33),st_ctl%lpmd(37)
  allocate(y(4*NCELL),yscal(4*NCELL),dydx(4*NCELL))
  allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL),idisp(NCELL),pf(NCELL),sigmae(NCELL),pfhyd(NCELL))
  phi=0d0;vel=0d0;tau=0d0;sigma=0d0;disp=0d0
  allocate(a(NCELL),b(NCELL),dc(NCELL),f0(NCELL),taudot(NCELL),tauddot(NCELL),sigdot(NCELL),ks(NCELL),qflow(NCELL),kLv(NCELL),kTv(NCELL))
  taudot=0d0;sigdot=0d0;ks=0d0

  !uniform
  lf=ds0*ncellg*cos(dipangle*pi/180)
  loc=ncellg/2

  ! number=number*10
  ! counts=0
  ! 777 number=number+1
  ! counts=counts+1
  ! b0=b0list(counts)

  a=a0
  b=b0
  dc=dc0
  f0=mu0

  !adding velocity-strengthening buffer at the edge
  ! do i=1,ncellg
  !   a(i)=max(a0,a0+(abs(i-800)-600)*2d-5)
  !   a(i)=min(a(i),0.015)
  !   !if(abs(xcol(i)-lf/2)>lf/2-lf/20) a(i)=a1
  !   !if(i<ncellg/10.or.i>ncellg*5/8) a(i)=a1
  !   !a(i)=a1
  ! end do

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

  !setting frictional parameters
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !stop
  !max time step
  !if(load==0) dtmax=0.02d0*10d0/sr

  if(foward) call foward_check()
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
      nout(1)=my_rank+100
      write(fname,'("output/vel",i0,"_",i0,".dat")') number,my_rank
      open(nout(1),file=fname,form='unformatted',access='stream')
      inquire(nout(1),size=file_size)
      m=file_size/8
      allocate(rdata(m))
      read(nout(1)) rdata
      vel(1:NCELL)=rdata(m-NCELL+1:m)
      !write(*,*) my_rank,vel
      close(nout(1))

      write(fname,'("output/slip",i0,"_",i0,".dat")') number,my_rank
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      disp(1:NCELL)=rdata(m-NCELL+1:m)
      close(nout(1))

      write(fname,'("output/sigma",i0,"_",i0,".dat")') number,my_rank
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      sigma(1:NCELL)=rdata(m-NCELL+1:m)
      close(nout(1))

      write(fname,'("output/tau",i0,"_",i0,".dat")') number,my_rank
      open(nout(1),file=fname,form='unformatted',access='stream')
      read(nout(1)) rdata
      tau(1:NCELL)=rdata(m-NCELL+1:m)
      close(nout(1))
      !write(*,*) my_rank,m

      phi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
      !write(*,*) tau

    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_bcast(phi,NCELL,MPI_REAL8,0,st_ctl%lpmd(35),ierr)
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
    !ks=(kpmax*vpl/kL+kpmin/kT)/(Vpl/kL+1/kT)
    !ks=(kpmax*vref/kL*(vpl/vref)**alpha+kpmin/kT)/(Vref/kL*(vpl/vref)**alpha+1/kT)
    !ks=kp0
    !kpmax=kp0*(1+kL/kT/Vpl)-kpmin*kL/kT/Vpl



    !depth dependent initial condition
    call initcond_thrust()
    !steady state background flow and
    !call initcond_ss()
    !injection-induced seismicity (BP&)
    !call initcond_injection()
    !tectonic swarm
    !call initcond_swarm()

    !pfhyd=sigma*0.4d0
    !pf(NCELL/2)=3d0
    !tau=tauinit
    qflow=q0


    !ks=1d-13*exp(50d0/30d0)
    !ks=ksinit

    !non-uniform initial stress from file
    if(parameterfromfile) then
      do i=1,NCELL
        i_=st_sum%lodc(i)
        tau(i)=taug(i_)
        sigma(i)=sigmag(i_)
        pfhyd(i)=(sigma(i)-0d0)*0.4d0+0d0
        pf(i)=pfg(i_)
        sigmae(i)=sigma(i)-pf(i)
        vel(i)=velg(i_)
        mu(i)=tau(i)/sigmae(i)
        phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(tau(i)/sigmae(i)/a(i)))
        !write(*,*) tau(i),sigma(i),vel(i)
      end do
    end if

    x=0d0
    kstart=1
    kend=0
    if(my_rank<npd) then
      write(fname,'("output/ind",i0,"_",i0,".dat")') number,my_rank
      nout(1)=my_rank+100
      open(nout(1),file=fname,form='unformatted',access='stream')
      write(nout(1))st_sum%lodc(1:NCELL)
      close(nout(1))

      write(fname,'("output/xyz",i0,"_",i0,".dat")') number,my_rank
      nout(1)=my_rank+100
      open(nout(1),file=fname)
      do i=1,ncell
        i_=st_sum%lodc(i)
        write(nout(1),'(3e15.6)') xcol(i_),ycol(i_),zcol(i_)
      end do
      close(nout(1))
      ! write(fname,'("output/prm",i0,"_",i0,".dat")') number,my_rank
      ! nout=my_rank+100
      ! open(nout,file=fname)
      ! do i=1,ncell
      !   write(nout,*) a(i),taudot(i),sigdot(i_)
      ! end do
      ! close(nout)
      write(fname,'("output/vel",i0,"_",i0,".dat")') number,my_rank
      open(nout(1),file=fname,form='unformatted',access='stream',status='replace')
      nout(2)=nout(1)+np
      write(fname,'("output/slip",i0,"_",i0,".dat")') number,my_rank
      open(nout(2),file=fname,form='unformatted',access='stream',status='replace')
      nout(3)=nout(2)+np
      write(fname,'("output/sigma",i0,"_",i0,".dat")') number,my_rank
      open(nout(3),file=fname,form='unformatted',access='stream',status='replace')
      nout(4)=nout(3)+np
      write(fname,'("output/tau",i0,"_",i0,".dat")') number,my_rank
      open(nout(4),file=fname,form='unformatted',access='stream',status='replace')
      nout(5)=nout(4)+np
      write(fname,'("output/ks",i0,"_",i0,".dat")') number,my_rank
      open(nout(5),file=fname,form='unformatted',access='stream',status='replace')
      nout(6)=nout(5)+np
      write(fname,'("output/pf",i0,"_",i0,".dat")') number,my_rank
      open(nout(6),file=fname,form='unformatted',access='stream',status='replace')
      nout(7)=nout(6)+np
      write(fname,'("output/qflow",i0,"_",i0,".dat")') number,my_rank
      open(nout(7),file=fname,form='unformatted',access='stream',status='replace')
    end if
    if(my_rank.eq.0) then
      write(fname,'("output/event",i0,".dat")') number
      open(51,file=fname)
      write(fname,'("output/monitor",i0,".dat")') number
      open(52,file=fname)
      do i=1,nl
        locid(i)=int(dble(i-0.5)/nl*ncellg)
      end do
      do i=1,nl
        write(fname,'("output/local",i0,"-",i0,".dat")') number,locid(i)
        open(52+i,file=fname)
      end do
      !write(fname,'("output/out",i0,".dat")') number
      !open(53,file=fname)
      !write(fname,'("output/vel",i0,".dat")') number
      !open(47,file=fname,form='unformatted',access='stream')

      open(19,file='job.log',position='append')
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
      write(19,'(a20,i0,a6,a12,a6,a12,a4,i0)') 'Starting job number=',number,'date',sys_time(1),'time',sys_time(2),'np',np
      close(19)
      !allocate(locid(9))
      !ocid=(/1850,2000,2050,2100,2150,2250,2350,2500,2750/)
      !call open_BP6()

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

    maxnorm=maxval(sigmae)
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(maxnorm,maxnormG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(maxnormG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    minnorm=minval(sigmae)
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
    if(my_rank<npd) then
      !call output_field()
      ! do i=1,ncell
      !   i_=st_sum%lodc(i)
      !   write(nout)i_
      ! end do
      write(nout(1)) vel
      write(nout(2)) disp
      write(nout(3)) sigmae
      write(nout(4)) tau
      write(nout(5)) ks
      write(nout(6)) pf
      write(nout(7)) qflow
    end if
    k=0
    errmax_gb=0d0
    dtdid=0d0
    if(my_rank==0) then
      call output_monitor()
      !call output_BP6()
    end if
    dtnxt = dtinit
  end if
  tout=dtout*365*24*60*60
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
    y(4*i-3) = phi(i)
    y(4*i-2) = tau(i)
    y(4*i-1) = sigmae(i)
    y(4*i) = ks(i)
    !if(my_rank==53)write(*,*) phi(i),tau(i),sigma(i)
  end do
  !$omp end parallel do
  !stop


  do k=kstart,NSTEP
    !parallel computing for Runge-Kutta
    dttry = dtnxt
    !time3=MPI_Wtime()

    !update state and stress with explicit solver
    call rkqs(y,dydx,x,dttry,eps_r,dtdid,dtnxt,errmax_gb,errmaxloc)
    do i=1,NCELL
      sigmae(i)=y(4*i-1)
      sigma(i)=y(4*i-1)+pf(i)
      ks(i)=y(4*i)
    end do
    write(*,*) sigma(ncell/2),sigmae(ncell/2),pf(ncell/2)

    !update pf with implicit solver
    !write(*,*)maxval(pf)
    !write(*,*)'call implicit solver'
    if(pressuredependent) then
      call implicitsolver(pf,sigma,ks,dtdid,x,dtnxt)
    else
      call implicitsolver2(pf,y,dtdid,x,dtnxt,niter)
    end if
    !write(*,*)maxval(pf)

    !time4=MPI_Wtime()
    !timer=timer+time4-time3

    !limitsigma
    ! if(limitsigma) then
    !   do i=1,NCELL
    !     if(y(3*i)<minsig) y(3*i)=minsig
    !     if(y(3*i)>maxsig) y(3*i)=maxsig
    !   end do
    ! end if

    !compute physical values for control and output
    !$omp parallel do
    do i = 1, NCELL
      !y(3*i)=sigma(i)-pf(i)
      phi(i) = y(4*i-3)
      tau(i) = y(4*i-2)
      sigmae(i)=sigma(i)-pf(i)
      y(4*i-1)=sigmae(i)
      ks(i)=y(4*i)
      disp(i)=disp(i)+vel(i)*dtdid*0.5d0 !2nd order
      vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigmae(i)/a(i))
      !if(abs(i-2000)<5)write(*,*) vel(i),tau(i)
      disp(i)=disp(i)+vel(i)*dtdid*0.5d0
      mu(i)=tau(i)/sigmae(i)
    end do
    !$omp end parallel do
    write(*,*)minval(sigmae)

    !computing flow rate
    qflow(1)=-ks(1)/eta*1e3*(pf(2)-pf(1))/ds0
    do i=2,ncell-1
      qflow(i)=-ks(i)/eta*1e3*(pf(i+1)-pf(i-1))/2/ds0
    end do
    qflow(ncell)=-ks(ncell)/eta*1e3*(pf(ncell)-pf(ncell-1))/ds0

    !correcting permeability by pressure dependence
    if(pressuredependent) then
      do i=1,ncell
        ks(i)=ks(i)*exp(-sigmae(i)/s0)
      end do
    end if

    write(*,*) sigma(ncellg),sigmae(ncellg),pf(ncellg)

    !write(*,*)minval(sigmae),minval(sigma)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr); time3=MPI_Wtime()

    mvel=maxval(abs(vel))
    mvelloc=maxloc(abs(vel))

    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(mvel,mvelG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(mvelG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)
    maxnorm=maxval(sigmae)
    !call MPI_ALLREDUCE(mvel,mvelG,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_reduce(maxnorm,maxnormG,1,MPI_REAL8,MPI_MAX,st_ctl%lpmd(37),st_ctl%lpmd(31),ierr)
    call MPI_bcast(maxnormG,1,MPI_REAL8,st_ctl%lpmd(33),st_ctl%lpmd(35),ierr)

    minnorm=minval(sigmae)
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
    if(mod(k,interval)==0) outfield=.true.
    if(outpertime.and.x>tout) then
      outfield=.true.
      !if(slipping) then
      !  tout=tout+5d0
      !else
        tout=tout+dtout*365*24*60*60
      !end if
    end if

    if(outfield) then
      if(my_rank==0) then
        write(*,'(a,i0,f17.8,a)') 'time step=' ,k,x/365/24/60/60, ' yr'
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
          write(fname,'("output/local",i0,"-",i0,".dat")') number,locid(i)
          open(52+i,file=fname,position='append')
        end do
      end if
      !lattice H
      if(my_rank<npd) then
        write(nout(1)) vel
        write(nout(2)) disp
        write(nout(3)) sigmae
        write(nout(4)) tau
        write(nout(5)) ks
        write(nout(6)) pf
        write(nout(7)) qflow
      end if
    end if

    if(my_rank==0.and.k>kend) then
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
        hypoloc=maxloc(abs(vel))
        onset_time=x
        tout=onset_time

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
        tout=x
        moment=meandispG-moment0
        !eventcount=eventcount+1
        !end of an event
        if(my_rank==0) then
          write(51,'(i0,i7,f17.2,i5,e15.6,f14.4)') eventcount,k,onset_time,hypoloc,moment,(log10(moment*rigid*sg)+5.9)/1.5
        end if
        if(slipevery.and.(my_rank<npd)) then
          ! write(nout) vel
          ! write(nout2) disp
          ! write(nout3) sigmae
          ! write(nout4) tau
          ! write(nout5) ks
        end if

      end if
      !   vmaxevent=max(vmaxevent,maxval(vel))
      !   !write(53,'(i6,4e16.6)') !k,x-onset_time,sum(disp-idisp),sum(vel),sum(acg**2)
      !   !if(x-onset_time>lapse) then
      !   !  lapse=lapse+dlapse
      !   !end if
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
  do i=1,10
    close(nout(i))
  end do
end if
!if(my_rank==0) write(19,'(a20,i0,f16.2)')'Finished job number=',number,time2-time1
Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! if(counts>=9) then
!   go to 888
! else
!   go to 777
! end if

888 select case(problem)
case('2dp','2dpah','2dn3','3dp')
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
    write(52,'(i7,f19.4,8e16.5,f16.4,i6)')k,x,log10(mvelG),meandispG,meanmuG,maxnormG,minnormG,xcol(mvelloc(1)),errmax_gb,dtdid,time2-time1,niter
  end subroutine
  subroutine output_local(nf,loc)
    integer,intent(in)::nf,loc
    implicit none
    write(nf,'(i7,f19.4,8e16.6)')k,x,log10(vel(loc)),disp(loc),sigmae(loc),tau(loc),pf(loc),mu(loc),phi(loc),ks(loc)
  end subroutine
  subroutine output_field()
    implicit none
    do i=1,NCELL
      i_=st_sum%lodc(i)
      write(nout(1),'(i0,10e14.5,i10)') i_,xcol(i_),ycol(i_),zcol(i_),log10(vel(i)),tau(i),sigmae(i),mu(i),disp(i),phi(i),x,k
    end do
    write(nout(1),*)
    write(nout(1),*)
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
  subroutine initcond_thrust()
    implicit none
    real(8),parameter::g=9.80,rhor=2.6d3,rhof=1.0d3,ofs=5.0,Qe=1d4,T0=1e0
    real(8)::q(ncellg),Temp,qtmp

    do i=1,ncellg
      a(i)=max(a0,a0+(abs(i-600)-400)*2d-5)
      a(i)=min(a(i),b0+0.002)
      !a(i)=0.02
      !if(abs(xcol(i)-lf/2)>lf/2-lf/20) a(i)=a1
      !if(i<ncellg/10.or.i>ncellg*5/8) a(i)=a1
      !a(i)=a1
      if(i<ncellg) qin(i)=q0/ncellg
    end do
    qtmp=sum(qin)
    q0=qtmp

    do i=1,ncellg
        q(i)=qtmp-qin(i)
        qtmp=q(i)
        Temp=300d0+12d0*ycol(i)
        kTv(i)=min(T0*exp(Qe/Temp),1d12)
        ks(i)=(kpmax*Vpl/kL+kpmin/kTv(i))/(Vpl/kL+1/kTv(i))
        write(*,*) ks(i),q(i),KTv(i)
    end do


    sigma=g*rhor*(ycol+0.25)*1e-3
    pfhyd=g*rhof*(ycol+0.25)*1e-3
    pf(1)=pfhyd(1)
    sigmae(1)=sigma(1)-pf(1)
    pbcl=pfhyd(1)

    do i=2,ncellg
      sigmae(i)=sigmae(i-1)+ds0*sin(dipangle*pi/180d0)*((rhor-rhof)*g*1e-3-q(i)*eta/ks(i)*1e-3*exp(sigmae(i-1)/s0))
      pf(i)=sigma(i)-sigmae(i)
      write(*,*)sigma(i),sigmae(i),pf(i)
    end do
    pbcr=pf(ncellg)
    vel=velinit
    tau=sigmae*(f0+(a-b)*log(vel/vref))*1.001
    !tau=sigmae*0.5
    ! vel=tau/abs(tau)*velinit
    mu=tau/sigmae
    phi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))
    disp=0d0
    !stop
  end subroutine

  subroutine initcond_ss()
    implicit none
    !pbcr is given
    !q0=(pbcr-pbcl)/lf*ks(1)/1d-12*1d-6*2
    kpmax=kp0*(1+kL/kT/Vpl)-kpmin*kL/kT/Vpl
    kTv=kT
    ks=kp0

    !q0 is given
    !q=q0
    pbcr=0d0
    lf=ncellg*ds0
    pbcl=pbcr+q0*lf/ks(1)*eta*1d-3


    pf=pbcl+(pbcr-pbcl)*xcol/lf
    sigma=sigmainit+pf
    sigmae=sigma-pf
    pfhyd=0.d0
    vel=velinit
    ! do i=1,ncell
    !   if(abs(xcol(i)-5.0)<1.0) vel(i)=velinit*(1.0+0.001*sin(2*pi/2.0*xcol(i)))
    ! end do
    tau=sigmae*(f0+(a-b)*log(vel/vref))
    ! vel=tau/abs(tau)*velinit
    mu=tau/sigmae
    phi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))
    disp=0d0
  end subroutine

  subroutine initcond_injection()
    !BP6
    implicit none
    qin(ncell/2)=1.25e-6
    ks=kpmin
    kTv=kT
    pf=pfinit
    sigma=sigmainit
    sigmae=sigma-pf
    pfhyd=0.d0
    vel=velinit
    tau=sigmae*muinit
    ! vel=tau/abs(tau)*velinit
    mu=muinit
    phi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))
    disp=0d0
  end subroutine

  subroutine initcond_swarm()
    implicit none
    real(8),parameter::g=9.80,rhor=2.6d3,rhof=1.0d3,ofs=5.0

    sigma=g*rhor*(ycol+0.25)*1e-3
    pfhyd=g*rhof*(ycol+0.25)*1e-3

    do i=1,int(ncell*4/5)
      pf(i)=pfhyd(i)
    end do
    do i=int(ncell*4/5)+1,ncell
      pf(i)=sigma(i)-sigmainit
    end do
    pbcl=pfhyd(1)
    kLv=kL
    kTv=kT
    ! do i=900,1000
    !   kLv(i)=10.0
    ! end do
    ks=kpmin
    sigmae=sigma-pf
    vel=velinit
    tau=sigmae*muinit
    ! vel=tau/abs(tau)*velinit
    mu=muinit
    phi=a*dlog(2*vref/vel*sinh(tau/sigmae/a))
    disp=0d0
  end subroutine

  subroutine implicitsolver(pf,sigma,ks,h,time,dtnxt)
    implicit none
    integer::kit,errloc(1),niter
    real(8),intent(inout)::pf(:),dtnxt
    real(8),intent(in)::h,time,sigma(:),ks(:)
    real(8)::dpf(ncell),pftry(ncell),pfnew(ncell),sigmae(ncell),cdiff(ncell),err,err0
    real(8),parameter::dpth=0.1
    !real(8),parameter::cc=1d-12 !beta(1e-8)*phi(1e-1)*eta(1e-3) !BP6
    !real(8),parameter::cc=1d-15 !beta(1e-9)*phi(1e-2)*eta(1e-4) !Zhu et al.


    do i=1,ncell
      pftry(i)=pf(i)
    end do

    err=0d0
    do kit=1,20
      !call CrankNicolson(pf,ks,sigma,pftry,h,pfnew)
      !calculate diffusion coefficient
      do i=1,ncell
        !cdiff(i)=((ks(i)-kpmin)*exp(-(sigma(i)-pftry(i))/s0)+kpmin)/cc*1d-6
        cdiff(i)=ks(i)*exp(-(sigma(i)-pftry(i))/s0)/cc*1d-6
        !cdiff(i)=ks(i)/cc*1d-6
      end do

      call Beuler(pf,cdiff,h,pfnew,time,niter)
      err=maxval(abs(pftry-pfnew))
      !write(*,*) pftry
      !write(*,*) pfnew
      errloc=maxloc(abs(pftry-pfnew))
      write(*,*) 'err',err
      if(err<1e-3) exit
      pftry=pfnew
      err0=err
    end do
    do i=1,ncell
      dpf(i)=pfnew(i)-pf(i)
      pf(i)=pfnew(i)
      !y(4*i-2)=y(4*i-2)-0.6*dpf(i)
    end do

    if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)
    return
  end subroutine

  subroutine implicitsolver2(pf,y,h,time,dtnxt,niter)
    implicit none
    integer::kit,errloc(1)
    integer,intent(out)::niter
    real(8),intent(inout)::pf(:),y(:),dtnxt
    real(8),intent(in)::h,time
    real(8)::dpf(ncell),pftry(ncell),pfnew(ncell),sigmae(ncell),ks(ncell),cdiff(ncell),err,err0
    real(8),parameter::dpth=0.1
    !real(8),parameter::cc=1d-12 !beta(1e-8)*phi(1e-1)*eta(1e-3)

    do i=1,ncell
      pftry(i)=pf(i)
      sigmae(i)=y(4*i-1)
      ks(i)=y(4*i)
      sigma(i)=sigmae(i)+pf(i)
    end do

    err=0d0

      !calculate diffusion coefficient
      do i=1,ncell
        !cdiff(i)=((ks(i)-kpmin)*exp(-(sigma(i)-pftry(i))/s0)+kpmin)/cc*1d-6
        !cdiff(i)=ks(i)*exp(-(sigma(i)-pftry(i))/s0)/cc*1d-6
        cdiff(i)=ks(i)/cc*1d-6 !Pa-->MPa
        !write(*,*) cdiff(i)
      end do

      call Beuler(pf,cdiff,h,pfnew,time,niter)

    do i=1,ncell
      dpf(i)=pfnew(i)-pf(i)
      pf(i)=pfnew(i)
      y(4*i-1)=y(4*i-1)-dpf(i) !update effective normal stress
      !y(4*i-2)=y(4*i-2)-0.6*dpf(i)
    end do
    write(*,*) 'niter',niter
    write(*,*) 'dpf',maxval(dpf)
    if(dtnxt/h*maxval(dpf)>dpth)  dtnxt=dpth*h/maxval(dpf)
    return
  end subroutine

  subroutine Beuler(pf,cdiff,h,pfnew,time,niter)
    integer,parameter::itermax=500
    real(8),intent(in)::pf(:),h,cdiff(:),time
    real(8),intent(out)::pfnew(:)
    integer,intent(out)::niter
    real(8)::eta,pf0(NCELL),m(ncell),dpf(ncell),p(ncell),r(ncell),b(ncell),x(ncell),tmp1,tmp2,rsnew,rsold,alpha(ncell)
    real(8)::Dxx(ncell,ncell),Am(ncell,3),sat(ncell)
    integer::n,iter
    real(8)::p0=0.0,td=1d6
    !real(8),parameter::str=1e-11 !beta(1e-9)*phi(1e-2)
    n=ncell
    niter=0
    ! do i=1,n
    !   !cdiff(i)=((ks(i)-kpmin)*exp(-(sigma(i)-pftry(i))/s0)+kpmin)/cc*1d-6
    !   cdiff(i)=ks(i)*exp(-(sigma(i)-pftry(i))/s0)/cc*1d-6
    !   !cdiff(i)=c0*exp(-pftry(i)/s0)
    !   !cdiff(i)=(1.d0+y(i))*1d-2
    ! end do
    Dxx=0d0

    !compute Dxx for Dirichlet BC
    select case(bcl)
    case('Dirichlet')
      Dxx(1,1)=-cdiff(1)-cdiff(2)
      Dxx(1,2)=cdiff(2)-cdiff(1)
      Dxx(2,1:3)=(/cdiff(2)/2-cdiff(1)/2, -cdiff(1)/2-cdiff(2)-cdiff(3)/2, cdiff(2)/2+cdiff(3)/2/)
    case('Neumann')
      Dxx(1,1)=-cdiff(1)-cdiff(2)
      Dxx(1,2)=-Dxx(1,1)
      Dxx(2,1:3)=(/cdiff(1)/2+cdiff(2)/2, -cdiff(1)/2-cdiff(2)-cdiff(3)/2, cdiff(2)/2+cdiff(3)/2/)
    end select

    do i=3,n-2
      Dxx(i,1:3)=(/cdiff(i-1)/2+cdiff(i)/2, -cdiff(i-1)/2-cdiff(i)-cdiff(i+1)/2, cdiff(i)/2+cdiff(i+1)/2/)
    end do

    select case(bcr)
    case('Dirichlet')
      Dxx(n-1,1:3)=(/cdiff(n-2)/2+cdiff(n-1)/2, -cdiff(n-2)/2-cdiff(n-1)-cdiff(n)/2, -cdiff(n)/2+cdiff(n-1)/2/)
      Dxx(n,3)=-cdiff(n)-cdiff(n-1)
      Dxx(n,2)=cdiff(n-1)-cdiff(n)
    case('Neumann')
      Dxx(n-1,1:3)=(/cdiff(n-2)/2+cdiff(n-1)/2, -cdiff(n-2)/2-cdiff(n-1)-cdiff(n)/2, cdiff(n-1)/2+cdiff(n)/2/)
      Dxx(n,3)=-cdiff(n)-cdiff(n-1)
      Dxx(n,2)=-Dxx(n,3)
    end select

    Dxx=Dxx/ds0/ds0

    !write(*,*)

    Am=0d0
    Am(1,1)=1.0-h*Dxx(1,1)
    Am(1,2)=-h*Dxx(1,2)
    Am(n,3)=1.0-h*Dxx(n,3)
    Am(n,2)=-h*Dxx(n,2)
    do i=2,n-1
      Am(i,1)=-h*Dxx(i,1)
      Am(i,2)=1.0-h*Dxx(i,2)
      Am(i,3)=-h*Dxx(i,3)
    end do
    !do i=0,n
    !  write(*,'(3e15.6)') A(i,:)
    !end do

    !membrane diffusion
    Am(n/2,2)=Am(n/2,2)!+h/td

    SAT=0d0

    !penalty vector for Neumann BC
    !SAT(1)=q0*h/ds0
    !SAT(n)=-q0*h/ds0

    !penalty vector for Dirichlet BC
    select case(bcl)
    case('Dirichlet')
      SAT(1)=-cdiff(1)/ds0/ds0*(pbcl-pfhyd(1))*h*2
      SAT(2)=SAT(1)/2
    case('Neumann')
      SAT(1)=-q0/str*1e-9*h/ds0
    end select

    select case(bcr)
    case('Dirichlet')
      SAT(n)=-cdiff(n)/ds0/ds0*(pbcr-pfhyd(n))*h*2
      SAT(n-1)=SAT(n)/2
    case('Neumann')
      SAT(n)=-q0/str*1e-9*h/ds0
    end select

    ! !source term
    !SAT(1500)=-0.5*q0/str*1e-9*h/ds0

    x=pf-SAT-pfhyd !initial guess
    b=pf-SAT-pfhyd

    !injection at the center of the fault (if nonzero qin is used)
    !if(time<30*24*3600)
    !constant fluid source
    ! do i=2,n-1
    !   b(i)=b(i)+h*0.5*qin(i)/str*1e-9/ds0
    ! end do

    !b(n/2)=b(n/2)+h*q0
    !write(*,*) b(n/2)
    if(time<100*24*3600) b(N/2)=b(N/2)+h*qin(n/2)/str*1e-9/ds0

    m=0d0
    m(1)=x(1)*Am(1,1)+x(2)*Am(1,2)
    do i=2,n-1
      m(i)=x(i-1)*Am(i,1)+x(i)*Am(i,2)+x(i+1)*Am(i,3)
    end do
    m(n)=x(n)*Am(n,3)+x(n-1)*Am(n,2)
    !write(*,'(9e15.6)')m

    r=b-m
    !write(*,'(9e15.6)')r
    p=r
    rsold=sum(r*r)
    if(rsold<1e-12*n)  then
      go to 100
    end if

    do iter=1,itermax
      tmp1=sum(r*r)
      m=0d0
      m(1)=p(1)*Am(1,1)+p(2)*Am(1,2)
      do i=2,n-1
        m(i)=p(i-1)*Am(i,1)+p(i)*Am(i,2)+p(i+1)*Am(i,3)
      end do
      m(n)=p(n)*Am(n,3)+p(n-1)*Am(n,2)
      !write(*,'(9e15.6)')m

      tmp2=sum(m*p)
      alpha=tmp1/tmp2
      x=x+alpha*p
      r=r-alpha*m
      !write(*,'(9e15.6)')r
      rsnew = sum(r*r)
      !write(*,*)iter,rsnew
      if(rsnew<1e-12*n) then
          niter=iter
        exit
      end if
      p = r + (rsnew / rsold) * p
      rsold = rsnew
      !write(*,'(9e15.6)')x

    end do
    100 pfnew=x+pfhyd
    return
  end subroutine

    subroutine CrankNicolson(pf,ks,sigma,pftry,h,pfnew)
      implicit none
      integer::i,j,ite,iter
      integer,parameter::itermax=300
      real(8),intent(in)::pf(:),ks(:),pftry(:),sigma(:),h
      real(8),intent(out)::pfnew(:)
      real(8)::eta,pf0(NCELL),Am(ncell,ncell),m(ncell),dpf(ncell),p(ncell),r(ncell),b(ncell),x(ncell),tmp1,tmp2,rsnew,rsold,alpha(ncell)
      real(8)::dvdtau,dvdsig,diag(ncell),left(ncell),right(ncell),q(ncell),gamma
      real(8),parameter::cc=1d-15 !beta(1e-9)*phi(1e-2)*eta(1e-4)

      do i=1,ncell
        cdiff=((ks(i)-kpmin)*exp(-(sigma(i)-pftry(i))/s0)+kpmin)/cc*1d-6
        alpha(i)=cdiff*h/(ds0**2)
        !write(*,*) cdiff
      end do

      !construction of matrix
      diag=0d0;left=0d0;right=0d0

      do i=1,ncell-1
        diag(i)=1d0+alpha(i)
      end do
      diag(ncell)=1d0+alpha(ncell)/2


      do i=2,ncell
        left(i)=-0.5d0*alpha(i)
        right(i-1)=-0.5d0*alpha(i)
      end do

      do i=1,ncell
        pf0(i)=pf(i)-pfhyd(i)
      end do

      q(1)=(1-alpha(1))*pf0(1)+alpha(2)/2*pf0(2)
      q(ncell)=alpha(ncell-1)/2*pf0(ncell-1)+(1-alpha(ncell)/2)*pf0(ncell)
      do i=2,ncell-1
        q(i)=0.5d0*(alpha(i+1)*pf0(i+1)+alpha(i-1)*pf0(i-1))+(1-alpha(i))*pf0(i)
      end do

      !r=pf
      !write(*,*)maxval(pf)
      x=q
      b=q
      b(NCELL)=b(NCELL)!+h*q0*1e4 !beta*phi/ds
      !if(time<100*24*3600) pf(NCELL/2)=pf0(NCELL/2)+h*q0

      m=0d0
      m(1)=diag(1)*x(1)+right(1)*x(2)
      m(ncell)=left(ncell)*x(ncell-1)+diag(ncell)*x(ncell)
      do i=2,ncell-1
        m(i)=left(i)*x(i-1)+diag(i)*x(i)+right(i)*x(i+1)
      end do

      r=b-m
      p=r
      rsold=sum(r*r)

      do iter=1,itermax
        tmp1=sum(r*r)

        m=0d0
        m(1)=diag(1)*p(1)+right(1)*p(2)
        m(ncell)=left(ncell)*p(ncell-1)+diag(ncell)*p(ncell)
        do i=2,ncell-1
          m(i)=left(i)*p(i-1)+diag(i)*p(i)+right(i)*p(i+1)
        end do

        tmp2=sum(m*p)
        gamma=tmp1/tmp2
        x=x+gamma*p
        r=r-gamma*m
        rsnew = sum(r*r)
        !write(*,*)iter,rsnew
        if(rsnew<1e-10*NCELL) then
          ite=iter
          exit
        end if
        p = r + (rsnew / rsold) * p
        rsold = rsnew
      end do
      eta=0.5d0*rigid/vs

      do i=1,NCELL
        pfnew(i)=x(i)+pfhyd(i)
        !m(i)=-tau(i)/sigmae(i)*dpf(i)
        !dvdtau=2*vref*dexp(-y(3*i-2)/a(i))*dcosh(y(3*i-1)/y(3*i)/a(i))/(a(i)*y(3*i))
        !dvdsig=-2*vref*dexp(-y(3*i-2)/a(i))*dcosh(y(3*i-1)/y(3*i)/a(i))*y(3*i-1)/(a(i)*y(3*i)**2)
        !m(i)=eta*dvdsig/(1d0+eta*dvdtau)*dpf(i)
        !y(3*i-1)=y(3*i-1)!-y(3*i-1)/y(3*i)*dpf(i)
      end do
      !write(*,*)maxval(m),maxval(dpf)

      return
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
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmaetmp(NCELL),phitmp(NCELL),kstmp(NCELL)
    real(8) :: sum_gs(NCELL),sum_gn(NCELL)!,velstmpG(NCELLg),veldtmpG(NCELLg)
    !real(8) :: sum_xx(NCELL),sum_xy(NCELL),sum_yy(NCELL)!,sum_xz(NCELL),sum_yz(NCELL),sum_zz(NCELL)
    !real(8) :: sum_xxG(NCELLg),sum_xyG(NCELLg),sum_yyG(NCELLg)!,sum_xzG(NCELLg),sum_yzG(NCELLg),sum_zzG(NCELLg)
    !real(8) :: sum_xx2G(NCELLg),sum_xy2G(NCELLg),sum_yy2G(NCELLg),sum_xz2G(NCELLg),sum_yz2G(NCELLg),sum_zz2G(NCELLg)
    !real(8) :: veltmpG(NCELLg),sum_gsg(NCELLg),sum_gng(NCELLg),sum_gdg(NCELLg)!,efftmpG(NCELLg)
    !real(8) :: sum_gs2G(NCELLg),sum_gd2G(NCELLg),sum_gn2G(NCELLg)
    real(8) :: time3,time4,c1, c2, c3, fss,arg,arg2,c,g,tauss,Arot(3,3),p(6),fac
    real(8)::dvdsig,sxx0,sxy0,syy0,dvdphi,dphidt,dsigdt,dtaudt,dvdtau,dksdt
    integer :: i, j, nc,ierr,lrtrn,i_

    !if(latticeh) then

    !$omp parallel do
    do i = 1, NCELL
      phitmp(i) = y(4*i-3)
      tautmp(i) = y(4*i-2)
      sigmaetmp(i) = y(4*i-1)
      kstmp(i)=y(4*i)
      veltmp(i) = 2*vref*dexp(-phitmp(i)/a(i))*dsinh(tautmp(i)/sigmaetmp(i)/a(i))
      !if(abs(i-2000)<5) write(*,*) veltmp(i)
    enddo
    !$omp end parallel do
    !write(*,*) maxval(kstmp)

    !matrix-vector mutiplation
    if(backslip) then
      st_vel%vs=veltmp-vpl
    else
      st_vel%vs=veltmp
    end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr);time3=MPI_Wtime()
    call HACApK_adot_lattice_hyp(st_sum,st_LHp_s,st_ctl,wws,st_vel)
    if(problem=='3dhr'.or.problem=='3dph') then
      sum_gs(:)=st_sum%vs(:)/ds0
    else
      sum_gs(:)=st_sum%vs(:)
    end if

    if(sigmaconst) then
      sum_gn=0d0
    else
      call HACAPK_adot_lattice_hyp(st_sum,st_LHP_n,st_ctl,wws,st_vel)
      if(problem=='3dhr'.or.problem=='3dph') then
        sum_gn(:)=st_sum%vs(:)/ds0
      else
        sum_gn(:)=st_sum%vs(:)
      end if
    end if
    !time4=MPI_Wtime()
    !timeH=timeH+time4-time3


    do i=1,NCELL
      !permeability evolution
      !dksdt=-veltmp(i)/kL*(kstmp(i)-kpmax)-(kstmp(i)-kpmin)/kT
      dksdt=-(veltmp(i)/vref)**alpha*vref/kL*(kstmp(i)-kpmax)-(kstmp(i)-kpmin)/kTv(i)

      sum_gs(i)=sum_gs(i)+taudot(i)
      sum_gn(i)=sum_gn(i)+sigdot(i)
      dsigdt=sum_gn(i)
      !dsigdt=sum_gn(i)-cdiff/(ds0**2)*(sigmatmp(i+1)-2*sigmatmp(i)+sigmatmp(i-1))
      !if(i==1) dsigdt=sum_gn(i)-cdiff/(ds0**2)*(sigmatmp(i+1)-sigmatmp(i))
      !if(i==NCELLg) dsigdt=sum_gn(i)-cdiff/(ds0**2)*(-sigmatmp(i)+sigmatmp(i-1))
      !if(i==NCELLg/2.and.x<100*24*3600) dsigdt=-q0

      !aging law
      select case(evlaw)
      case('aging')
        dphidt=b(i)/dc(i)*vref*dexp((f0(i)-phitmp(i))/b(i))-b(i)*veltmp(i)/dc(i)
      !slip law
      case('slip')
        fss=f0(i)+(a(i)-b(i))*dlog(abs(veltmp(i))/vref)
        dphidt=-abs(veltmp(i))/dc(i)*(abs(tautmp(i))/sigmaetmp(i)-fss)
      end select

      dvdtau=2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmaetmp(i)/a(i))/(a(i)*sigmaetmp(i))
      dvdsig=-2*vref*dexp(-phitmp(i)/a(i))*dcosh(tautmp(i)/sigmaetmp(i)/a(i))*tautmp(i)/(a(i)*sigmaetmp(i)**2)
      dvdphi=-veltmp(i)/a(i)
      dtaudt=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt+dvdsig*dsigdt)
      !if(i==NCELL/2) dtaudt=dtaudt-tautmp(i)/sigmatmp(i)*q0
      dtaudt=dtaudt/(1d0+0.5d0*rigid/vs*dvdtau)

      dydx(4*i-3)=dphidt
      dydx(4*i-2)=dtaudt
      dydx(4*i-1)=dsigdt
      dydx(4*i)=dksdt
      !dydx(4*i)=0d0
      !if(i==1) write(*,*) dksdt
      !call deriv_d(sum_gs(i),sum_gn(i),phitmp(i),tautmp(i),sigmatmp(i),veltmp(i),a(i),b(i),dc(i),f0(i),dydx(3*i-2),dydx(3*i-1),dydx(3*i))
    enddo
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
      do i=1,4*ncell
        if(abs(yerr(i)/ytemp(i))/eps>errmax) then
          errmax=abs(yerr(i)/ytemp(i))/eps
          errmaxloc=i
        end if!errmax=errmax+yerr(3*i-2)**2
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

    hnext=min(2*h,SAFETY*h*(errmax_gb**PGROW),5d7)
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

  subroutine rk2(y,x,h,yout,yerr)
    implicit none
    real(8),intent(in)::y(:),x,h
    real(8),intent(out)::yout(:),yerr(:)
    real(8)::ytemp(4*NCELL),y1(4*NCELL),ak1(4*NCELL),ak2(4*NCELL)
    integer::i

    !1st-order solution
    call derivs(x, y, ak1)
    !$omp parallel do
    do i=1,size(y)
      y1(i)=y(i)+h*ak1(i)
    end do
    !$omp end parallel do

    !2nd-order solution
    !$omp parallel do
    do i=1,size(y)
      ytemp(i)=y(i)+0.5d0*h*ak1(i)
    end do
    !$omp end parallel do
    call derivs(x+0.5*h, ytemp, ak2)

    !call beuler(pf,h,iter)

    !$omp parallel do
    do i=1,size(y)
      yout(i)=y(i)+h*ak2(i)
      yerr(i)=yout(i)-y1(i)
    end do
    !$omp end parallel do

    return
  end subroutine

  subroutine foward_check()
    implicit none
    real(8)::rr,lc,ret1(NCELLg),ret2(NCELLg),vec(NCELLg)
    integer::p

    vec=0d0
    do i=1,NCELL
      i_=st_sum%lodc(i)
      if(ycol(i_)>5d0) vec(i)=1d0
    end do
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
      if(problem=='3dhr'.or.problem=='3dph') ret1(:)=ret1(:)/ds0

      if(.not.sigmaconst) then
        !lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_n,st_bemv,st_ctl,ret2,vec)
        st_vel%vs=vec
        call HACApK_adot_lattice_hyp(st_sum,st_LHp_n,st_ctl,wws,st_vel)
        ret2(:)=st_sum%vs(:)
        if(problem=='3dhr'.or.problem=='3dph') ret2(:)=ret2(:)/ds0
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

  subroutine open_BP6()
    real(8)::xd(201)

  open(31,file="output/fltst_strk-15",status='replace')
  open(32,file="output/fltst_strk+00",status='replace')
  open(33,file="output/fltst_strk+05",status='replace')
  open(34,file="output/fltst_strk+10",status='replace')
  open(35,file="output/fltst_strk+15",status='replace')
  open(36,file="output/fltst_strk+25",status='replace')
  open(37,file="output/fltst_strk+35",status='replace')
  open(38,file="output/fltst_strk+50",status='replace')
  open(39,file="output/fltst_strk+75",status='replace')

  open(40,file="output/global.dat",status='replace')


  do i=31,39
    write(i,*)"# This is the header:"
    write(i,*)"# problem=SEAS Benchmark BP6-A"
    write(i,*)"# code=HBI"
    write(i,*)"# modeler=So Ozawa"
    write(i,*)"# date=2022/07/11"
    write(i,*)"# element_size=10m"
    write(i,*)"# location= on fault, 0km down-dip distance"
    write(i,*)"# Column #1 = Time (s)"
    write(i,*)"# Column #2 = Slip (m)"
    write(i,*)"# Column #3 = Slip rate (log10 m/s)"
    write(i,*)"# Column #4 = Shear stress (MPa)"
    write(i,*)"# Column #5 = Pore pressure (MPa)"
    write(i,*)"# Column #6 = Darcy velocity (m/s)"
    write(i,*)"# Column #7 = State (log10 s)"
    write(i,*)"# The line below lists the names of the data fields"
    write(i,*)"t slip slip_rate shear_stress pore_pressure darcy_vel state"
    write(i,*)"# Here is the time-series data."
  end do

  write(40,*)"# This is the file header:"
  write(40,*)"# problem=SEAS Benchmark BP6-A"
  write(40,*)"# code=HBI"
  write(40,*)"# modeler=So Ozawa"
  write(40,*)"# date=2022/07/11"
  write(40,*)"# element_size=10m"
  write(40,*)"# Column #1 = Time (s)"
  write(40,*)"# Column #2 = Max Slip rate (log10 m/s)"
  write(40,*)"# Column #3 = Moment rate (N-m/s)"
  write(40,*)"# The line below lists the names of the data fields"
  write(40,*)"t max_slip_rate moment_rate"
  write(40,*)"# Here is the time-series data."

  open(41,file="output/slip.dat",status='replace')
  open(42,file="output/shear_stress.dat",status='replace')
  open(43,file="output/normal_stress.dat",status='replace')
  open(44,file="output/darcy_vel.dat",status='replace')

  do i=41,44
  write(i,*)"# This is the file header:"
  write(i,*)"# problem=SEAS Benchmark BP6-A"
  write(i,*)"# code=HBI"
  write(i,*)"# modeler=So Ozawa"
  write(i,*)"# date=2022/07/11"
  write(i,*)"# element_size=10m"
  write(i,*)"# Column #1 = Time (s)"
  write(i,*)"# Column #2 = Max Slip rate (log10 m/s)"
  end do

  write(41,*)"# Column #3-203 = Slip (m)"
  write(42,*)"# Column #3-203 = Shear stress (MPa)"
  write(43,*)"# Column #3-203 = Normal stress (MPa)"
  write(44,*)"# Column #3-203 = Darcy velocity (MPa)"

  do i=41,44
  write(i,*)"Computational domain size: -20 km < z < 20 km"

  write(i,*)"# The line below lists the names of the data fields"
  write(i,*)"z"
  write(i,*)"t"
  write(i,*)"max_slip_rate"

  end do
  write(41,*)"slip"
  write(42,*)"shear_stress"
  write(43,*)"normal_stress"
  write(44,*)"darcy_vel"

  do i=41,44
  write(i,*)"# Here are the data."
  end do
  do i=1,201
    xd(i)=-10000d0+(i-1)*100d0
  end do
  write(41,'(203e22.14)') 0d0,0d0,xd
  write(42,'(203e22.14)') 0d0,0d0,xd
  write(43,'(203e22.14)') 0d0,0d0,xd
  write(44,'(203e22.14)') 0d0,0d0,xd
return
end subroutine
subroutine output_BP6()
  implicit none
  do i=1,9
    write(30+i,'(e22.14,6e15.7)') x,disp(locid(i)),log10(abs(vel(locid(i)))),tau(locid(i)),pf(locid(i)),&
    &qflow(locid(i)),log10(dc(locid(i))/vref*exp((phi(locid(i))-f0(locid(i)))/b(locid(i))))
  end do
  write(40,'(e22.14,2e15.7)') x,log10(maxval(vel)),sum(vel)*rigid*1e9*ds0*1e3
  if(mod(k,interval)==3) then
    write(41,'(203e22.14)') x,maxval(log10(abs(vel))),disp(1000:3000:10)
    write(42,'(203e22.14)') x,maxval(log10(abs(vel))),tau(1000:3000:10)
    write(43,'(203e22.14)') x,maxval(log10(abs(vel))),sigmae(1000:3000:10)
    write(44,'(203e22.14)') x,maxval(log10(abs(vel))),qflow(1000:3000:10)
  end if
end subroutine

end program
