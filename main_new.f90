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
  integer::NCELL, nstep1, lp, i,i_,j,k,m,counts,interval,number,lrtrn,nl,NCELLg,ios
  integer::clock,cr,counts2,imax,jmax,NCELLm,seedsize,icomm,np,ierr,my_rank
  integer::hypoloc(1),load,eventcount,thec,inloc
  logical::slipping,outfield,limitsigma
  integer,allocatable::seed(:)
  character*128::fname,dum,law,input_file,problem,geofile,param,pvalue
  real(8)::a0,b0,dc0,sr,omega,theta,dtau,tiny,x,time1,time2,moment,aslip,avv
  real(8)::psi,vc0,mu0,dtinit,onset_time,tr,vw0,fw0,velmin,muinit,intau,errmax_gb
  real(8)::r,eps,vpl,outv,xc,zc,dr,dx,dz,lapse,dlapse,vmaxeventi,sparam,tmax,eps_r,eps_h
  real(8)::dtime,dtnxt,dttry,dtdid,dtmin,alpha,ds,amp,mui,strinit,velinit,phinit,velmax
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxps,st_leafmtxpn
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xx,st_leafmtxp_xy,st_leafmtxp_yy
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xz,st_leafmtxp_yz,st_leafmtxp_zz
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xx2,st_leafmtxp_xy2,st_leafmtxp_yy2
  type(st_HACApK_leafmtxp) :: st_leafmtxp_xz2,st_leafmtxp_yz2,st_leafmtxp_zz2
  type(st_HACApK_calc_entry) :: st_bemv
  !type(t_deriv)::
  real(8),allocatable ::coord(:,:),vmax(:)
  real(8),allocatable::a(:),b(:),dc(:),vc(:),fw(:),vw(:),ac(:),taudot(:),tauddot(:),sigdot(:)
  !real(8),allocatable::phi(:),vel(:),tau(:),sigma(:),disp(:),mu(:)
  real(8),allocatable::phi(:),vel(:),tau(:),sigma(:),disp(:),mu(:),rupt(:),idisp(:)
  real(8),allocatable::taus(:),taud(:),vels(:),veld(:),disps(:),dispd(:),rake(:)
  real(8),allocatable::xcol(:),ycol(:),zcol(:)
  real(8),allocatable::xs1(:),xs2(:),xs3(:),xs4(:) !for 3dp
  real(8),allocatable::zs1(:),zs2(:),zs3(:),zs4(:) !for 3dp
  real(8),allocatable::ys1(:),ys2(:),ys3(:) !for 3dn
  real(8),allocatable::xel(:),xer(:),yel(:),yer(:),ang(:)
  real(8),allocatable::ev11(:),ev12(:),ev13(:),ev21(:),ev22(:),ev23(:),ev31(:),ev32(:),ev33(:)
  real(8),allocatable::y(:),yscal(:),dydx(:),xcoll(:),zcoll(:),yg(:)
  !real(8),allocatable::xr(:),yr(:),zr(:),ai(:,:)
  integer::r1,r2,r3,NVER,amari,out,kmax,loci,locj,loc,stat
  integer,allocatable::displs(:),rcounts(:),vmaxin(:),rupsG(:),vars(:)

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
  !if(my_rank.eq.0) write(*,*) 'Reading input file'
  open(33,file=input_file,iostat=ios)
  if (ios /= 0) then
     write(*,*) 'Failed to open inputfile'
     stop
  end if

  !information of fault geometry
  ! read(33,*) dum,problem
  !
  ! select case(problem)
  ! case('2dp')
  !   read(33,*) dum,NCELLg
  ! case('2dn','3dn','3dh')
  !   read(33,*) dum,NCELLg
  !   read(33,*) dum,geofile
  ! case('3dp')
  !   read(33,*) dum,imax !x length
  !   read(33,*) dum,jmax !z length
  ! end select

  !read(33,*) dum,kmax !number of VS patches
  !read(33,*) dum,dr !radius of VS patches (unit:ds)
  ! read(33,*) dum,ds !mesh interval(normalized by Dc)

  !output control
  ! read(33,*) dum,number !output filename
  ! read(33,*) dum,interval !output frequency
  !read(33,*) dum,dlapse !output per time
  !read(33,*) dum,loci !output localdata i
  !read(33,*) dum,locj !output localdata j

  !continue or stop control
  ! read(33,*) dum,nstep1 !maxmimum time step
  ! !read(33,*) dum,thec !stop when the eventcount exceeds this value
  ! read(33,*) dum,velmax !stop when slip rate exceeds this value
  ! read(33,*) dum,velmin

  !physical parameters in calculation
  ! read(33,*) dum,law ! evolution law in RSF
  ! read(33,*) dum,a0 !a in RSF
  ! read(33,*) dum,b0 !b in RSF
  ! read(33,*) dum,dc0
  ! read(33,*) dum,vw0
  ! read(33,*) dum,fw0
  ! read(33,*) dum,vc0 !cut-off velocity in RSF (should be large enough when assuming conventional RSF)
  ! read(33,*) dum,mu0 !mu0 in RSF
  ! read(33,*) dum,load !loading type 0: uniform stress 1: uniform slip deficit
  ! read(33,*) dum,sr !if load=0: stressing rate (sigma0/sec) load=1: unused
  ! read(33,*) dum,vpl !if load=1: loading velocity load=0: unused
  !read(33,*) dum,tr !viscoelastic relaxation (unused)
  !read(33,*) dum,rigid !shear modulus normalized by sigma0
  !read(33,*) dum,vs !Swave speed (dc/s)
  !read(33,*) dum,pois !poisson ratio

  !initial values & nucleation
  ! read(33,*) dum,velinit !initial slip velocity
  ! read(33,*) dum,muinit !initial omega=V*theta
  ! read(33,*) dum,psi !stress angle
  ! read(33,*) dum,dtinit !initial timestep
  ! read(33,*) dum,intau !initial timestep
  ! read(33,*) dum,inloc !initial timestep
  !
  ! read(33,*) dum,limitsigma
  ! read(33,*) dum,sparam !for aftershock difference of main_sub fault
  ! read(33,*) dum,eps !error allowance in time integration in Runge-Kutta
  !read(*,*) amp
  !read(*,*) omega

  !for FDMAP
  !read(33,*) coordinate_file

  !new input reading system(under construction)
  eps_r=1d-5
  eps_h=1d-5
  tmax=1d12
  limitsigma=.false.

  do while(ios==0)
    read(33,*,iostat=ios) param,pvalue
    !write(*,*) param,pvalue
    select case(param)
    case('problem')
      read (pvalue,*) problem
    case('NCELLg')
      read (pvalue,*) ncellg
    case('NSTEP1')
      read (pvalue,*) nstep1
    case('filenumber')
      read (pvalue,*) number
    case('ds')
      read (pvalue,*) ds
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
    case('psi')
      read (pvalue,*) psi
    case('dtinit')
      read (pvalue,*) dtinit
    case('intau')
      read (pvalue,*) intau
    case('inloc')
      read (pvalue,*) inloc
    case('limitsigma')
      read (pvalue,*) limitsigma
    case('sparam')
      read (pvalue,*) sparam
    case('tmax')
      read (pvalue,*) tmax
    case('eps_r')
      read (pvalue,*) eps_r
    case('eps_h')
      read (pvalue,*) eps_h
    end select
  end do
  close(33)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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
  allocate(vars(NCELL))
  do i=1,NCELL
    vars(i)=displs(my_rank+1)+i
    !write(*,*) displs(my_rank+1),i,vars(i)
  end do

  !stop
  !call varscalc(NCELL,displs,vars)
  if(my_rank.eq.0) write(*,*) rcounts,displs

  !allocation
  allocate(a(NCELLg),b(NCELLg),dc(NCELLg),vc(NCELLg),fw(NCELLg),vw(NCELLg),taudot(NCELLg),tauddot(NCELLg),sigdot(NCELLg))
  allocate(rupt(NCELLg),rupsG(NCELLg))

  select case(problem)
  case('2dp','2dpv2','2dpv3')
    allocate(xcol(NCELLg),xel(NCELLg),xer(NCELLg))
    xcol=0d0;xel=0d0;xer=0d0
    !allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL))
    allocate(phi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg))
  case('2dn','2dn3')
    allocate(xcol(NCELLg),ycol(NCELLg),ang(NCELLg))
    allocate(xel(NCELLg),xer(NCELLg),yel(NCELLg),yer(NCELLg))
    xcol=0d0;ycol=0d0;ang=0d0;xel=0d0;xer=0d0;yel=0d0;yer=0d0
    !allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL))
    allocate(phi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg))
  case('3dp')
    allocate(xcol(NCELLg),zcol(NCELLg))
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
    xcol=0d0; zcol=0d0
    xs1=0d0; xs2=0d0; xs3=0d0; xs4=0d0
    zs1=0d0; zs2=0d0; zs3=0d0; zs4=0d0
    !allocate(phi(NCELL),vel(NCELL),tau(NCELL),sigma(NCELL),disp(NCELL),mu(NCELL))
    allocate(phi(NCELLg),vel(NCELLg),tau(NCELLg),sigma(NCELLg),disp(NCELLg),mu(NCELLg),idisp(NCELLg))
  case('3dn','3dh')
    allocate(xcol(NCELLg),ycol(NCELLg),zcol(NCELLg))
    allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg))
    allocate(ys1(NCELLg),ys2(NCELLg),ys3(NCELLg))
    allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg))
    allocate(ev11(NCELLg),ev12(NCELLg),ev13(NCELLg))
    allocate(ev21(NCELLg),ev22(NCELLg),ev23(NCELLg))
    allocate(ev31(NCELLg),ev32(NCELLg),ev33(NCELLg))
    xcol=0d0; ycol=0d0; zcol=0d0
    xs1=0d0; xs2=0d0; xs3=0d0
    ys1=0d0; ys2=0d0; ys3=0d0
    zs1=0d0; zs2=0d0; zs3=0d0
    !allocate(phi(NCELL),vels(NCELL),veld(NCELL),taus(NCELL),taud(NCELL),sigma(NCELL),disps(NCELL),dispd(NCELL),mu(NCELL))
    allocate(phi(NCELLg),vels(NCELLg),veld(NCELLG),taus(NCELLg),taud(NCELLg),sigma(NCELLg),disps(NCELLg),dispd(NCELLG),mu(NCELLg),rake(NCELLg),vel(NCELLG),tau(NCELLg),idisp(NCELLg))
  end select

  select case(problem) !for Runge-Kutta
  case('2dp','2dn3','3dp')
    allocate(y(2*NCELL),yscal(2*NCELL),dydx(2*NCELL),yg(2*NCELLg))
  case('2dn')
    allocate(y(3*NCELL),yscal(3*NCELL),dydx(3*NCELL),yg(3*NCELLg))
  case('3dn','3dh')
    allocate(y(4*NCELL),yscal(4*NCELL),dydx(4*NCELL),yg(4*NCELLg))
  end select

  !allocate(vmax(NCELLg),vmaxin(NcELLg))

  !mesh generation (rectangular assumed)
  if(my_rank.eq.0) write(*,*) 'Generating mesh'
    select case(problem)
    case('2dp')
      call coordinate2dp(NCELLg,ds,xel,xer,xcol)
    case('2dn','2dn3') !geometry file is necessary
      call coordinate2dn(geofile,NCELLg,xel,xer,yel,yer,xcol,ycol,ang)
    case('3dp')
      call coordinate3dp(imax,jmax,ds,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
    case('3dn','3dh')
      call coordinate3dn(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
      call evcalc(xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3,ev11,ev12,ev13,ev21,ev22,ev23,ev31,ev32,ev33)
    end select

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
  case('2dp')
    allocate(st_bemv%xcol(NCELLg),st_bemv%xel(NCELLg),st_bemv%xer(NCELLg))
    st_bemv%xcol=xcol;st_bemv%xel=xel;st_bemv%xer=xer
    st_bemv%problem=problem

  case('2dn','2dn3')
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

  case('3dn','3dh')
    allocate(st_bemv%xcol(NCELLg),st_bemv%ycol(NCELLg),st_bemv%zcol(NCELLg))
    allocate(st_bemv%xs1(NCELLg),st_bemv%xs2(NCELLg),st_bemv%xs3(NCELLg))
    allocate(st_bemv%ys1(NCELLg),st_bemv%ys2(NCELLg),st_bemv%ys3(NCELLg))
    allocate(st_bemv%zs1(NCELLg),st_bemv%zs2(NCELLg),st_bemv%zs3(NCELLg))
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
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,eps_h)

  case('2dn3')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=ycol(i)
      coord(i,3)=0.d0
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,eps_h)

  case('2dn')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=ycol(i)
      coord(i,3)=0.d0
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    ! st_bemv%v='s'
    ! lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,1d-4)
    ! st_bemv%v='n'
    ! lrtrn=HACApK_generate(st_leafmtxpn,st_bemv,st_ctl,coord,1d-4)
    st_bemv%v='xx'
    lrtrn=HACApK_generate(st_leafmtxp_xx,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='xy'
    lrtrn=HACApK_generate(st_leafmtxp_xy,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='yy'
    lrtrn=HACApK_generate(st_leafmtxp_yy,st_bemv,st_ctl,coord,eps_h)

  case('3dp')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=0.d0
      coord(i,3)=zcol(i)
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    lrtrn=HACApK_generate(st_leafmtxps,st_bemv,st_ctl,coord,eps_h)

  case('3dn','3dh')
    do i=1,NCELLg
      coord(i,1)=xcol(i)
      coord(i,2)=ycol(i)
      coord(i,3)=zcol(i)
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !kernel for strike slip
    st_bemv%md='st'
    st_bemv%v='xx'
    lrtrn=HACApK_generate(st_leafmtxp_xx,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='xy'
    lrtrn=HACApK_generate(st_leafmtxp_xy,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='yy'
    lrtrn=HACApK_generate(st_leafmtxp_yy,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='xz'
    lrtrn=HACApK_generate(st_leafmtxp_xz,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='yz'
    lrtrn=HACApK_generate(st_leafmtxp_yz,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='zz'
    lrtrn=HACApK_generate(st_leafmtxp_zz,st_bemv,st_ctl,coord,eps_h)

    !kernel for dip slip
    st_bemv%md='dp'
    st_bemv%v='xx'
    lrtrn=HACApK_generate(st_leafmtxp_xx2,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='xy'
    lrtrn=HACApK_generate(st_leafmtxp_xy2,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='yy'
    lrtrn=HACApK_generate(st_leafmtxp_yy2,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='xz'
    lrtrn=HACApK_generate(st_leafmtxp_xz2,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='yz'
    lrtrn=HACApK_generate(st_leafmtxp_yz2,st_bemv,st_ctl,coord,eps_h)
    st_bemv%v='zz'
    lrtrn=HACApK_generate(st_leafmtxp_zz2,st_bemv,st_ctl,coord,eps_h)
  end select

  !setting frictional parameters
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if(my_rank.eq.0) write(*,*) 'Setting fault parameters'
  call params(problem,NCELLg,a0,b0,dc0,vc0,a,b,dc,vc,fw,vw)
  call loading(problem,NCELLg,sr,taudot,tauddot,sigdot)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


  !setting initial condition
  select case(problem)
  case('2dp','3dp','2dn3')
    sigma=sigma0
    tau=sigma*muinit
    mu=tau/sigma
    vel=tau/abs(tau)*velinit
    phi=a*dlog(2*vref/vel*sinh(tau/sigma/a))
    !omega=exp((phi(1)-mu0)/b(1))*abs(vel(1))/vref/b(1)
    !write(*,*) 'Omega',Omega
  case('2dn')
    call initcond2d(psi,muinit,phi,sigma,tau,disp)
    call add_nuclei(tau,intau,inloc)
  case('3dn','3dh')
    call initcond3d(phi,sigma,taus,taud)
  end select
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
    write(fname,'("output/slip",i0,".dat")') number
    open(46,file=fname)
    write(fname,'("output/event",i0,".dat")') number
    open(44,file=fname)
    write(fname,'("output/local",i0,".dat")') number
    open(42,file=fname)
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
  rupt=0d0
  rupsG=0
  dtnxt = dtinit
  !outv=1d-6
  slipping=.false.
  eventcount=0

  time2=MPI_Wtime()
  select case(problem)
  case('2dp','2dn','2dn3','3dp')
  write(52,'(i7,f18.5,3e16.5,i7,e16.5,f16.4)')0,x,maxval(log10(abs(vel))),sum(disp)/NCELLg,sum(mu)/NCELLg,maxloc(abs(vel)),log10(maxval(vel(1:10000))),time2-time1
  case('3dn','3dh')
    write(52,'(i7,f18.5,3e16.5,i7,e16.5,f16.4)')0,x,maxval(log10(vel)),sum(disps)/NCELLg,sum(mu)/NCELLg,maxloc(vel),dtdid*maxval(vel),time2-time1
  end select

  select case(problem)
  case('3dp')
    do i=1,NCELLg
      write(50,'(6e15.6,i10)') xcol(i),zcol(i),log10(vel(i)),mu(i),disp(i),k
    end do
    write(50,*)
    write(50,*)
  case('2dn','2dn3')
    do i=1,NCELLg
      write(50,'(i0,9e15.6,i10)') i,xcol(i),ycol(i),log10(abs(vel(i))),tau(i),sigma(i),mu(i),disp(i),phi(i),x,k
    end do
    write(50,*)
  case('2dp')
    do i=1,NCELLg
      write(50,'(i0,8e15.6,i10)') i,xcol(i),log10(abs(vel(i))),tau(i),sigma(i),mu(i),disp(i),phi(i),x,k
    end do
    write(50,*)
  case('3dn','3dh')
    do i=1,NCELLg
      write(50,'(12e14.5,i10)') xcol(i),ycol(i),zcol(i),log10(vel(i)),taus(i),taud(i),phi(i),mu(i),sigma(i),disps(i),dispd(i),rake(i),k
    end do
    write(50,*)
    write(50,*)
  end select

  !do i=1,NCELLg
  !  write(50,'(8e15.6,i6)') xcol(i),ycol(i),vel(i),tau(i),sigma(i),mu(i),disp(i),x,k
  !end do
  !write(50,*)
  select case(problem)
  case('2dp','2dn3','3dp')
    do i=1,NCELL
      i_=vars(i)
      y(2*i-1) = phi(i_)
      y(2*i) = tau(i_)
    end do
    !call MPI_SCATTERv(yG,2*rcounts,2*displs,MPI_REAL8,y,2*NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  case('2dn')
    do i=1,NCELL
      i_=vars(i)
      !write(*,*) my_rank,i_
      y(3*i-2) = phi(i_)
      y(3*i-1) = tau(i_)
      y(3*i)=sigma(i_)
    end do
    !call MPI_SCATTERv(yG,3*rcounts,3*displs,MPI_REAL8,y,3*NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  case('3dn','3dh')
    do i=1,NCELL
      i_=vars(i)
      !write(*,*) my_rank,i_
      y(4*i-3) = phi(i_)
      y(4*i-2) = taus(i_)
      y(4*i-1) = taud(i_)
      y(4*i)=sigma(i_)
    end do
    !call MPI_SCATTERv(yG,4*rcounts,4*displs,MPI_REAL8,y,4*NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  end select
!stop
  do k=1,NSTEP1
    dttry = dtnxt

    call derivs(x, y, dydx)!,,st_leafmtxps,st_leafmtxpn,st_bemv,st_ctl)
    do i = 1, size(yscal)
      yscal(i)=abs(y(i))+abs(dttry*dydx(i))!+tiny
    end do

    !parallel computing for Runge-Kutta
    call rkqs(y,dydx,x,dttry,eps_r,yscal,dtdid,dtnxt,errmax_gb)

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    select case(problem)
    case('2dp','2dn3','3dp')
       call MPI_ALLGATHERv(y,2*NCELL,MPI_REAL8,yG,2*rcounts,2*displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
      do i = 1, NCELLg
        !i=vars(i_)
        !write(*,*) my_rank,i
        phi(i) = yg(2*i-1)
        tau(i) = yg(2*i)
        !write(*,*) my_rank,i,phi(i)
        !disp(i) = disp(i)+exp(y(2*i-1))*dtdid
        vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
        !write(*,*)vel(i),dtdid
        disp(i)=disp(i)+vel(i)*dtdid
        mu(i)=tau(i)/sigma(i)
      end do
    case('2dn')
       call MPI_ALLGATHERv(y,3*NCELL,MPI_REAL8,yG,3*rcounts,3*displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
      do i = 1, NCELLg
        phi(i) = yG(3*i-2)
        tau(i) = yG(3*i-1)
        sigma(i) = yG(3*i)
        vel(i)= 2*vref*dexp(-phi(i)/a(i))*dsinh(tau(i)/sigma(i)/a(i))
        disp(i)=disp(i)+vel(i)*dtdid
        !s(i)=a(i)*dlog(2.d0*vref/vel(i)*dsinh(tau(i)/sigma(i)/a(i)))
        !s(i)=exp((tau(i)/sigma(i)-mu0-a(i)*dlog(vel(i)/vref))/b(i))
        mu(i)=tau(i)/sigma(i)
      end do
    case('3dn','3dh')
      call MPI_ALLGATHERv(y,4*NCELL,MPI_REAL8,yG,4*rcounts,4*displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
      do i = 1, NCELLg
        phi(i) = yG(4*i-3)
        taus(i) = yG(4*i-2)
        taud(i) = yG(4*i-1)
        sigma(i) = yG(4*i)
        !disp(i) = disp(i) + exp( y(3*i-2) i)*dtdid
        tau(i)=sqrt(taus(i)**2+taud(i)**2)
        vel(i)= 2*vref*dexp(-phi(i)/a(i))*dsinh(tau(i)/sigma(i)/a(i))
        vels(i)= vel(i)*taus(i)/tau(i)
        veld(i)= vel(i)*taud(i)/tau(i)
        disps(i)=disps(i)+vels(i)*dtdid
        dispd(i)=dispd(i)+veld(i)*dtdid
        rake(i)=atan2(veld(i),vels(i))/pi*180d0
        mu(i)=sqrt(taus(i)**2+taud(i)**2)/sigma(i)
      end do

    end select

    Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !stop

    !output
    if(my_rank.eq.0) then
      time2= MPI_Wtime()
      select case(problem)
      case('2dp','2dn','2dn3','3dp')
      write(52,'(i7,f19.4,3e16.5,i7,e16.5,f16.4)')k,x,maxval(log10(abs(vel))),sum(disp)/NCELLg,sum(mu)/NCELLg,maxloc(abs(vel)),log10(maxval(vel(1:10000))),time2-time1
      case('3dn','3dh')
        write(52,'(i7,f18.5,3e16.5,i7,e16.5,f16.4)')k,x,maxval(log10(vel)),sum(disps)/NCELLg,sum(mu)/NCELLg,maxloc(vel),dtdid*maxval(vel),time2-time1
      end select
      !do i=1,size(vmax)
      !  vmax(i)=max(vmax(i),vel(i))
      !  vmaxin(i)=k
      !end do
      !if(mod(k,interval).eq.0) then

      !for FDMAP
      !PsiG=a*dlog(2.d0*vref/vel*dsinh(tau/sigma/a))
      do i=1,NCELLg
        if(abs(vel(i)).gt.1d-2.and.rupt(i).le.1d-6) then
          rupt(i)=x
          !rupsG(i)=k
        end if
      end do


      !output distribution control
      outfield=.false.
      !A : iteration number
      if(mod(k,interval).eq.0) outfield=.true.
      !if(k.lt.18000) out=1

      !B : slip velocity
      !if(maxval(vel).gt.outv) then
      !  out=0
      !  outv=outv*(10.d0)**(0.5d0)
      !end if

      if(outfield) then
        write(*,*) 'time step=' ,k
        select case(problem)
        case('3dp')
          do i=1,NCELLg
            write(50,'(6e15.6,i10)') xcol(i),zcol(i),log10(vel(i)),mu(i),disp(i),k
          end do
          write(50,*)
          write(50,*)
        case('2dn','2dn3')
          do i=1,NCELLg
            write(50,'(i0,9e15.6,i10)') i,xcol(i),ycol(i),log10(abs(vel(i))),tau(i),sigma(i),mu(i),disp(i),phi(i),x,k
          end do
          write(50,*)
        case('2dp')
          do i=1,NCELLg
            write(50,'(i0,8e15.6,i10)') i,xcol(i),log10(abs(vel(i))),tau(i),sigma(i),mu(i),disp(i),phi(i),x,k
          end do
          write(50,*)
        case('3dn','3dh')
          do i=1,NCELLg
            write(50,'(12e14.5,i10)') xcol(i),ycol(i),zcol(i),log10(vel(i)),taus(i),taud(i),phi(i),mu(i),sigma(i),disps(i),dispd(i),rake(i),k
          end do
          write(50,*)
          write(50,*)
        end select
      end if

    end if


    !event list
    if(.not.slipping) then
     if(maxval(abs(vel)).gt.1d-2) then
         slipping=.true.
         eventcount=eventcount+1
         idisp=disp
         hypoloc=maxloc(abs(vel))
         onset_time=x
         if(my_rank.eq.0) then
           do i=1,NCELLg
             write(46,*) i,disp(i)
           end do
           write(46,*)
         end if
    !     lapse=0.d0
    !     if(my_rank.eq.0) write(44,*) eventcount,x,maxloc(abs(vel))
    !     if(my_rank.eq.0) write(fname,'("output/event",i0,".dat")') number
    !     if(my_rank.eq.0) open(53,file=fname)
       end if
     end if
    !
    if(slipping) then
      if(maxval(abs(vel)).lt.5d-3) then
         slipping=.false.
         moment=sum(disp-idisp)
         !eventcount=eventcount+1
         if(my_rank.eq.0) then
           write(44,'(i0,f19.4,i7,e15.6)') eventcount,onset_time,hypoloc,moment
           do i=1,NCELLg
             write(46,*) i,disp(i)
           end do
           write(46,*)
        end if
      end if
    !   vmaxevent=max(vmaxevent,maxval(vel))
    !   !write(53,'(i6,4e16.6)') !k,x-onset_time,sum(disp-idisp),sum(vel),sum(acg**2)
    !   !if(x-onset_time.gt.lapse) then
    !   !  lapse=lapse+dlapse
    !   !end if
    end if

    !simulation ends before nstep1 when
    !(1) eventcount exceeds threshold
    ! if(eventcount.eq.thec) then
    !   if(my_rank .eq. 0) write(*,*) 'eventcount 10'
    !   go to 200
    ! end if

    !(2) slip velocity exceeds threshold (for nucleation)
    if(maxval(abs(vel)).gt.velmax) then
      if(my_rank .eq. 0) write(*,*) 'slip rate above vmax'
      exit
    end if

    if(maxval(abs(vel)).lt.velmin) then
      if(my_rank .eq. 0) write(*,*) 'slip rate below vmin'
      exit
    end if
    if(x.gt.tmax) then
      if(my_rank .eq. 0) write(*,*) 'time exceeds tmax'
      exit
    end if

    dttry = dtnxt
  end do


  !output for FDMAP communication
  !call output_to_FDMAP()

   if(my_rank.eq.0) then
     do i=1,NCELLg
       if(problem.eq.'2dp') write(46,*) i,disp(i)
       if(problem.eq.'2dn') write(46,'(4f16.4)') xcol(i),ycol(i),disp(i),ang(i)
     end do
     do i=10076,NCELLg,150
       if(problem.eq.'2dp') write(48,*) i,rupt(i)
       if(problem.eq.'2dn') write(48,'(4f16.4)') xcol(i),ycol(i),rupt(i),ang(i)
     end do
   end if

  200  if(my_rank.eq.0) then
  time2= MPI_Wtime()
  write(*,*) 'time(s)', time2-time1
  close(52)
  close(50)
  close(48)
  close(46)
  close(44)
  close(19)
end if
Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
select case(problem)
case('2dp','2dn3','3dp')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxps)
case('2dn')
  !lrtrn=HACApK_free_leafmtxp(st_leafmtxps)
  !lrtrn=HACApK_free_leafmtxp(st_leafmtxpn)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xx)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xy)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_yy)
case('3dn','3dh')
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xx)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xy)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_yy)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xz)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_yz)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_zz)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xx2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xy2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_yy2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_xz2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_yz2)
  lrtrn=HACApK_free_leafmtxp(st_leafmtxp_zz2)
end select
lrtrn=HACApK_finalize(st_ctl)
Call MPI_FINALIZE(ierr)
stop
contains
  subroutine initcond2d(psi,muinit,phi,sigma,tau,disp)
    implicit none
    real(8),intent(in)::psi,muinit
    real(8),intent(out)::phi(:),sigma(:),tau(:),disp(:)
    real(8)::sxx0,sxy0,syy0,theta
    disp=0d0
    !initial tractions from uniform stress tensor
    syy0=sigma0
    sxy0=syy0*muinit
    !psi=37d0
    !psi=30d0
    !psi=42d0
    sxx0=syy0*(1d0+2*sxy0/(syy0*dtan(2*psi/180d0*pi)))
    write(*,*) 'sxx0,sxy0,syy0'
    write(*,*) sxx0,sxy0,syy0
    if(my_rank.eq.0) open(16,file='initomega')
    do i=1,size(vel)
      !i_=vars(i)
        tau(i)=sxy0*cos(2*ang(i))+0.5d0*(sxx0-syy0)*sin(2*ang(i))
        sigma(i)=sin(ang(i))**2*sxx0+cos(ang(i))**2*syy0+sxy0*sin(2*ang(i))
        !constant velocity
        !vel(i)=velinit*tau(i)/abs(tau(i))
        !phi(i)=a(i)*dlog(2*vref/velinit*sinh(abs(tau(i))/sigma(i)/a(i)))

        !constant Phi
        phi(i)=phinit
        vel(i)= 2*vref*exp(-phi(i)/a(i))*sinh(tau(i)/sigma(i)/a(i))
        omega=exp((phi(i)-mu0)/b(i))*vel(i)/vref/b(i)
        if(my_rank.eq.0) write(16,'(4e16.4)') ang(i)*180/pi,omega,log10(abs(vel(i))),tau(i)/sigma(i)
    end do
    if(my_rank.eq.0) close(16)

    !predefined sigma and tau(debug)
    !sigma=sigma0
    !tau=sigma*muinit
    !vel=velinit
    !disp=0d0

  end subroutine
  subroutine initcond3d(phi,sigma,taus,taud)
    implicit none
    real(8),intent(out)::phi(:),sigma(:),taus(:),taud(:)
    real(8)::vel(NCELLg),PS11,PS22,PS33,PS12

    !uniform
    sigma=sigma0
    taus=sigma*muinit
    !taus=28d0
    taud=0d0
    vel=taus/abs(taus)*velinit
    phi=a*dlog(2*vref/vel*sinh(sqrt(taus**2+taud**2)/sigma/a))
    ! if(my_rank.eq.0) then
    ! do i=1,NCELLg
    !   write(*,*) i,taus(i),phi(i)
    ! end do!omega=exp((phi(1)-mu0)/b(1))*abs(vel(1))/vref/b(1)
    ! end if

    !uniform tensor in a full-space
    PS11=sigma0
    PS22=sigma0
    PS33=sigma0
    PS12=PS22*muinit
    do i=1,NCELLg
      taus(i) = ev11(i)*ev31(i)*PS11 + ev12(i)*ev32(i)*PS22+ (ev11(i)*ev32(i)+ev12(i)*ev31(i))*PS12 + ev13(i)*ev33(i)*PS33
      taud(i) = ev21(i)*ev31(i)*PS11 + ev22(i)*ev32(i)*PS22+ (ev21(i)*ev32(i)+ev22(i)*ev31(i))*PS12 + ev23(i)*ev33(i)*PS33
      sigma(i) = ev31(i)*ev31(i)*PS11 + ev32(i)*ev32(i)*PS22+ (ev31(i)*ev32(i)+ev32(i)*ev31(i))*PS12 + ev33(i)*ev33(i)*PS33
      vel(i)=velinit
      phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(sqrt(taus(i)**2+taud(i)**2)/sigma(i)/a(i)))
    end do

    !depth dependent stress in a half-space
    ! do i=1,NCELLg
    !   PS11=-zcol(i)*17d0
    !   PS22=PS11
    !   PS33=PS11
    !   PS12=PS22*muinit
    !   taus(i) = ev11(i)*ev31(i)*PS11 + ev12(i)*ev32(i)*PS22+ (ev11(i)*ev32(i)+ev12(i)*ev31(i))*PS12 + ev13(i)*ev33(i)*PS33
    !   taud(i) = ev21(i)*ev31(i)*PS11 + ev22(i)*ev32(i)*PS22+ (ev21(i)*ev32(i)+ev22(i)*ev31(i))*PS12 + ev23(i)*ev33(i)*PS33
    !   sigma(i) = ev31(i)*ev31(i)*PS11 + ev32(i)*ev32(i)*PS22+ (ev31(i)*ev32(i)+ev32(i)*ev31(i))*PS12 + ev33(i)*ev33(i)*PS33
    !   vel(i)=velinit
    !   phi(i)=a(i)*dlog(2*vref/vel(i)*sinh(sqrt(taus(i)**2+taud(i)**2)/sigma(i)/a(i)))
    ! end do

  end subroutine

  subroutine add_nuclei(tau,intau,inloc)
    implicit none
    real(8),intent(in)::intau
    integer,intent(in)::inloc
    real(8),intent(inout)::tau(:)
    real(8)::ra
    integer::lc
    ra=sqrt((xcol(2)-xcol(1))**2+(ycol(2)-ycol(1))**2)
    lc=int(0.15d0*rigid*(1.d0-pois)*dc0/(b0-a0)/sigma0/ra)
    write(*,*) 'lc=',lc
    tau(inloc-lc:inloc+lc)=tau(inloc-lc:inloc+lc)+intau*tau(inloc)/abs(tau(inloc))

  end subroutine

  subroutine coordinate2dp(NCELLg,ds,xel,xer,xcol)
      implicit none
      integer,intent(in)::NCELLg
      real(8),intent(in)::ds
      real(8),intent(out)::xel(:),xer(:),xcol(:)
      integer::i,j,k

      !flat fault with element size ds
      do i=1,NCELLg
        xel(i)=(i-1)*ds
        xer(i)=i*ds
        xcol(i)=0.5d0*(xel(i)+xer(i))
        !write(14,'(3e16.6)') xcol(i),xel(i),xer(i)
      enddo
      !close(14)
      return
  end subroutine

  subroutine coordinate2dn(geofile,NCELLg,xel,xer,yel,yer,xcol,ycol,ang)
    implicit none
    integer,intent(in)::NCELLg
    character(128),intent(in)::geofile
    real(8),intent(out)::xel(:),xer(:),yel(:),yer(:),xcol(:),ycol(:),ang(:)
    integer::i,j,k,file_size,n,Np,Nm,ncellf
    real(8)::dx,xr(0:NCELLg),yr(0:NCELLg),nx(NCELLg),ny(NCELLg),r(NCELLg)
    real(8),allocatable::data(:)

    !reading mesh data from mkelm.f90
    open(20,file=geofile,access='stream')
    read(20) xel,xer,yel,yer

    !computing local angles and collocation points
    do i=1,NCELLg
      ang(i)=datan2(yer(i)-yel(i),xer(i)-xel(i))
      xcol(i)=0.5d0*(xel(i)+xer(i))
      ycol(i)=0.5d0*(yel(i)+yer(i))
    end do

    !output to file
    ! open(14,file='top3.dat')
    ! do i=1,NCELLg
    !   write(14,'(7e16.6)') xcol(i),ycol(i),ang(i),xel(i),xer(i),yel(i),yer(i)
    ! end do
    ! close(14)

    return
  end subroutine

  subroutine coordinate3dp(imax,jmax,ds,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
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
    return
  end subroutine coordinate3dp

  subroutine coordinate3dn(NCELLg,xcol,ycol,zcol,xs1,xs2,xs3,ys1,ys2,ys3,zs1,zs2,zs3)
    implicit none
    integer,intent(in)::NCELLg
    real(8),intent(out)::xcol(:),ycol(:),zcol(:)
    real(8),intent(out)::xs1(:),xs2(:),xs3(:),ys1(:),ys2(:),ys3(:),zs1(:),zs2(:),zs3(:)
    real(4)::xl(0:2048,0:2048)
    real(8),parameter::amp=1d-4
    integer::i,j,k
    logical::rough

    !reading mesh data from in_fgeom.dat of mkelm.c of Ando's code
    open(20,file=geofile)
    do i=1,NCELLg
    read(20,*) k,xs1(i),ys1(i),zs1(i),xs2(i),ys2(i),zs2(i),xs3(i),ys3(i),zs3(i),xcol(i),ycol(i),zcol(i)
    !bump
    ! xs1(i)=1d0*bump(ys1(i),zs1(i))
    ! xs2(i)=1d0*bump(ys2(i),zs2(i))
    ! xs3(i)=1d0*bump(ys3(i),zs3(i))
    ! xcol(i)=(xs1(i)+xs2(i)+xs3(i))/3.d0
    end do

rough=.true.
    !rough fault
    if(rough) then
    open(30,file='roughsurf.txt')
    do k=0,2048
      read(30,*) xl(k,0:2048)
    end do
    close(30)
    if(my_rank.eq.0) open(32,file='tmp')
    do i=1,NCELLg
      j=int((ys1(i)+10)*102.4)
      k=int(-102.4*zs1(i))
      xs1(i)=xl(j,k)*amp
      j=int((ys2(i)+10)*102.4)
      k=int(-102.4*zs2(i))
      xs2(i)=xl(j,k)*amp
      j=int((ys3(i)+10)*102.4)
      k=int(-102.4*zs3(i))
      xs3(i)=xl(j,k)*amp
      xcol(i)=(xs1(i)+xs2(i)+xs3(i))/3.d0
      if(my_rank.eq.0) write(32,*) xcol(i),ycol(i),zcol(i)
    end do
    end if
    return
  end subroutine coordinate3dn

  function bump(y,z)
    implicit none
    real(8)::y,z,bump,rr
    rr=(y-5d0)**2+(z+10d0)**2
    bump=max(exp(-rr/4d0)-1d-4,0d0)
    return
  end function

  subroutine params(problem,NCELLg,a0,b0,dc0,vc0,a,b,dc,vc,fw,vw)
    character(128),intent(in)::problem
    integer,intent(in)::NCELLg
    real(8),intent(in)::a0,b0,dc0,vc0
    real(8),intent(out)::a(:),b(:),dc(:),vc(:),fw(:),vw(:)
    integer::i

    !uniform
    do i=1,NCELLg
      a(i)=a0
      !if(abs(xcol(i)-50d0).gt.30d0) a(i)=0.024d0 !for cycle
      b(i)=b0
      dc(i)=dc0
      !if((problem.eq.'2dn').and.i.gt.10000) dc(i)=sparam*dc0
      vc(i)=vc0
      fw(i)=fw0
      vw(i)=vw0
      !write(*,*)a(i),b(i),dcg(i)
    end do

    !depth-dependent
    ! do i=1,NCELLg
    !   if(zcol(i).gt.-4d0) then
    !     a(i)=0.015d0+0.0025d0*(zcol(i)+4d0)
    !   else if(zcol(i).gt.-15d0) then
    !     a(i)=0.015d0
    !   else
    !     a(i)=0.015d0-0.0025d0*(zcol(i)+15d0)
    !   end if
    !   b(i)=0.020d0
    !   dc(i)=dc0
    !   vc(i)=vc0
    !   fw(i)=fw0
    !   vw(i)=vw0
    ! end do

  end subroutine

  subroutine loading(problem,NCELLg,sr,taudot,tauddot,sigdot)
    character(128),intent(in)::problem
    integer,intent(in)::NCELLg
    real(8),intent(in)::sr
    real(8),intent(out)::taudot(:),tauddot(:),sigdot(:)
    real(8)::factor,edge,ret1,ret2
    integer::i
    character(128)::v
    select case(problem)
    case('2dp','2dn3','3dp')
      taudot=sr
      tauddot=0d0
      sigdot=0d0
    case('2dn')
      !open(15,file='sr')
      do i=1,NCELLg
      select case(load)
      case(0)
        taudot(i)=sr*cos(2*ang(i))
        sigdot(i)=sr*sin(2*ang(i))
      case(1)
        !edge=ds*NCELLg/2

          v='s'
          call kern(v,xcol(i),ycol(i),-500d0,ycol(1),0d0,ycol(1),ang(i),0d0,ret1)
          call kern(v,xcol(i),ycol(i),100d0,ycol(10000),600d0,ycol(10000),ang(i),0d0,ret2)
          taudot(i)=vpl*(ret1+ret2)

          v='n'
          call kern(v,xcol(i),ycol(i),-500d0,ycol(1),0d0,ycol(1),ang(i),0d0,ret1)
          call kern(v,xcol(i),ycol(i),100d0,ycol(10000),600d0,ycol(10000),ang(i),0d0,ret2)
          sigdot(i)=vpl*(ret1+ret2)
      end select
      !write(15,*) taudot(i),sigdot(i)
    end do
    !close(15)
    tauddot=0d0
  case('3dn','3dh')
    !taudot=sr
    taudot=0d0
    tauddot=0d0
    sigdot=0d0
    !case('2dpv','2dnv')
    !  factor=rigid/(2.d0*pi*(1.d0-pois))
    !  edge=-ds*NCELLg
    !  do i=1,NCELLg
    !    taudotg(i)=vpl*factor*(1.d0/xcol(i)-1.d0/(xcol(i)-edge))
    !    sigdotG(i)=0d0
    !  end do
    end select
  end subroutine
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
        !if(my_rank.eq.0) write(*,*) ev31(k),ev32(k),ev33(k)
        rr = sqrt(ev31(k)*ev31(k)+ev32(k)*ev32(k)+ev33(k)*ev33(k))
        !// unit vectors for local coordinates of elements
        ev31(k) = ev31(k)/rr ; ev32(k) = ev32(k)/rr ; ev33(k) = ev33(k)/rr
        !write(*,*) ev31(k),ev32(k),ev33(k)

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
        !write(*,*)ev21(k),ev22(k),ev23(k)
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
    real(8) :: veltmp(NCELL),tautmp(NCELL),sigmatmp(NCELL),efftmp(NCELL),phitmp(NCELL)
    real(8) :: dlnvdt(NCELL),dtaudt(NCELL),dsigdt(NCELL),deffdt(NCELL),dphidt(NCELL)
    real(8) :: taustmp(NCELL),taudtmp(NCELL),velstmp(NCELL),veldtmp(NCELL),dtausdt(NCELL),dtauddt(NCELL)
    real(8) :: sum_gs(NCELL),sum_gn(NCELL),sum_gd(NCELL),velstmpG(NCELLg),veldtmpG(NCELLg)
    real(8) :: sum_xx(NCELL),sum_xy(NCELL),sum_yy(NCELL),sum_xz(NCELL),sum_yz(NCELL),sum_zz(NCELL)
    real(8) :: sum_xxG(NCELLg),sum_xyG(NCELLg),sum_yyG(NCELLg),sum_xzG(NCELLg),sum_yzG(NCELLg),sum_zzG(NCELLg)
    real(8) :: sum_xx2G(NCELLg),sum_xy2G(NCELLg),sum_yy2G(NCELLg),sum_xz2G(NCELLg),sum_yz2G(NCELLg),sum_zz2G(NCELLg)
    real(8)::veltmpG(NCELLg),sum_gsg(NCELLg),sum_gng(NCELLg),efftmpG(NCELLg)
    real(8):: c1, c2, c3, arg,c,g,tauss,Arot(3,3),p(6)
    integer :: i, j, nc,ierr,lrtrn,i_

    !if(my_rank.eq.0) then
    select case(problem)
    case('2dp','2dn3','3dp')
      do i = 1, NCELL
        i_=vars(i)
        phitmp(i) = y(2*i-1)
        tautmp(i) = y(2*i)
        sigmatmp(i)=sigma0 !normal stress is constant for planar fault
        veltmp(i) = 2*vref*dexp(-phitmp(i)/a(i_))*dsinh(tautmp(i)/sigmatmp(i)/a(i_))
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      select case(load)
      case(0)
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
        call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        do i=1,NCELL
          i_=vars(i)
          sum_gn(i)=0.d0
          sum_gs(i)=sum_gs(i)+taudot(i_)
        end do
      case(1)
        lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG-vpl)
        call MPI_SCATTERv(sum_gsg,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      end select

      call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)

      do i = 1, NCELL
        dydx(2*i-1) = dphidt(i)
        dydx(2*i) = dtaudt(i)
      enddo
      !call MPI_SCATTERv(sum_gsG,NCELL,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    case('2dn')
      do i = 1, NCELL
        i_=vars(i)
        phitmp(i) = y(3*i-2)
        tautmp(i) = y(3*i-1)
        sigmatmp(i) = y(3*i)
        veltmp(i) = 2*vref*dexp(-phitmp(i)/a(i_))*dsinh(tautmp(i)/sigmatmp(i)/a(i_))
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(Veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,                &
      &     MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      st_bemv%v='xx'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xx,st_bemv,st_ctl,sum_xxG,veltmpG)
      st_bemv%v='xy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xy,st_bemv,st_ctl,sum_xyG,veltmpG)
      st_bemv%v='yy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yy,st_bemv,st_ctl,sum_yyG,veltmpG)

      call MPI_SCATTERv(sum_xxG,rcounts,displs,MPI_REAL8,sum_xx,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_xyG,rcounts,displs,MPI_REAL8,sum_xy,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_yyG,rcounts,displs,MPI_REAL8,sum_yy,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      do i=1,NCELL
        i_=vars(i)
        sum_gs(i)=0.5d0*(sum_xx(i)-sum_yy(i))*dsin(-2*ang(i_))+sum_xy(i)*dcos(-2*ang(i_))
        !sum_gn(i)=0.5d0*(sum_xx(i)+sum_yy(i))-0.5d0*(sum_xx(i)-sum_yy(i))*dcos(2*ang(i))-sum_xy(i)*dsin(2*ang(i))
        sum_gn(i)=-(0.5d0*(sum_xx(i)+sum_yy(i))-0.5d0*(sum_xx(i)-sum_yy(i))*dcos(2*ang(i_))-sum_xy(i)*dsin(2*ang(i_)))
      end do
      do i=1,NCELL
        i_=vars(i)
        sum_gs(i)=sum_gs(i)+taudot(i_)
        sum_gn(i)=sum_gn(i)+sigdot(i_)
      end do
      call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)
      do i = 1, NCELL
        dydx(3*i-2) = dphidt(i)
        dydx(3*i-1) = dtaudt(i)
        dydx(3*i) = dsigdt(i)
        if(limitsigma.and.(sigmatmp(i).lt.30d0)) dsigdt=0d0
        if(limitsigma.and.(sigmatmp(i).gt.170d0)) dsigdt=0d0
      enddo

    case('2dn_vector')
      do i = 1, NCELL
        i_=vars(i)
        phitmp(i) = y(3*i-2)
        tautmp(i) = y(3*i-1)
        sigmatmp(i) = y(3*i)
        veltmp(i) = 2*vref*dexp(-phitmp(i)/a(i_))*dsinh(tautmp(i)/sigmatmp(i)/a(i_))
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(Veltmp,NCELL,MPI_REAL8,veltmpG,rcounts,displs,                &
      &     MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      st_bemv%v='s'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxps,st_bemv,st_ctl,sum_gsG,veltmpG)
      st_bemv%v='n'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxpn,st_bemv,st_ctl,sum_gnG,veltmpG)


      call MPI_SCATTERv(sum_gsG,rcounts,displs,MPI_REAL8,sum_gs,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_gnG,rcounts,displs,MPI_REAL8,sum_gn,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

      do i=1,NCELL
        sum_gs(i)=sum_gs(i)+taudot(i)
        sum_gn(i)=sum_gn(i)+sigdot(i)
      end do
      call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)

      do i = 1, NCELL
        dydx(3*i-2) = dphidt(i)
        dydx(3*i-1) = dtaudt(i)
        dydx(3*i) = dsigdt(i)
      enddo

    case('3dn','3dh')
      do i = 1, NCELL
        i_=vars(i)
        phitmp(i) = y(4*i-3)
        taustmp(i) = y(4*i-2)
        taudtmp(i) = y(4*i-1)
        sigmatmp(i) = y(4*i)
        tautmp(i)=sqrt(taustmp(i)**2+taudtmp(i)**2)
        veltmp(i)=2*vref*dexp(-phitmp(i)/a(i_))*dsinh(tautmp(i)/sigmatmp(i)/a(i_))
        velstmp(i)=veltmp(i)*taustmp(i)/tautmp(i)
        veldtmp(i)=veltmp(i)*taudtmp(i)/tautmp(i)
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(Velstmp,NCELL,MPI_REAL8,velstmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERv(Veldtmp,NCELL,MPI_REAL8,veldtmpG,rcounts,displs,MPI_REAL8,MPI_COMM_WORLD,ierr)

      !matrix-vector mutiplation
      st_bemv%md='st'
      st_bemv%v='xx'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xx,st_bemv,st_ctl,sum_xxG,velstmpG)
      st_bemv%v='xy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xy,st_bemv,st_ctl,sum_xyG,velstmpG-vpl)
      st_bemv%v='yy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yy,st_bemv,st_ctl,sum_yyG,velstmpG)
      st_bemv%v='xz'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xz,st_bemv,st_ctl,sum_xzG,velstmpG)
      st_bemv%v='yz'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yz,st_bemv,st_ctl,sum_yzG,velstmpG)
      st_bemv%v='zz'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_zz,st_bemv,st_ctl,sum_zzG,velstmpG)
      !write(*,*) 'max_sum',maxval(sum_xyG)

      st_bemv%md='dp'
      st_bemv%v='xx'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xx2,st_bemv,st_ctl,sum_xx2G,veldtmpG)
      st_bemv%v='xy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xy2,st_bemv,st_ctl,sum_xy2G,veldtmpG)
      st_bemv%v='yy'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yy2,st_bemv,st_ctl,sum_yy2G,veldtmpG)
      st_bemv%v='xz'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_xz2,st_bemv,st_ctl,sum_xz2G,veldtmpG)
      st_bemv%v='yz'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_yz2,st_bemv,st_ctl,sum_yz2G,veldtmpG)
      st_bemv%v='zz'
      lrtrn=HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp_zz2,st_bemv,st_ctl,sum_zz2G,veldtmpG)
      sum_xxG=sum_xxG+sum_xx2G
      sum_xyG=sum_xyG+sum_xy2G
      sum_yyG=sum_yyG+sum_yy2G
      sum_xzG=sum_xzG+sum_xz2G
      sum_yzG=sum_yzG+sum_yz2G
      sum_zzG=sum_zzG+sum_zz2G
      ! if(my_rank.eq.0) then
      !   open(32,file='tmp')
      ! do i=1,NCELLg
      ! write(32,*) xcol(i),ycol(i),sum_xxG(i)
      ! end do
      ! end if

      call MPI_SCATTERv(sum_xxG,rcounts,displs,MPI_REAL8,sum_xx,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_xyG,rcounts,displs,MPI_REAL8,sum_xy,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_yyG,rcounts,displs,MPI_REAL8,sum_yy,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_xzG,rcounts,displs,MPI_REAL8,sum_xz,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_yzG,rcounts,displs,MPI_REAL8,sum_yz,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTERv(sum_zzG,rcounts,displs,MPI_REAL8,sum_zz,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      !stop
      !stress --> traction
      do i=1,NCELL
        i_=vars(i)
        Arot(1,:)=(/ev11(i_),ev21(i_),ev31(i_)/)
        Arot(2,:)=(/ev12(i_),ev22(i_),ev32(i_)/)
        Arot(3,:)=(/ev13(i_),ev23(i_),ev33(i_)/)
        call TensTrans(sum_xx(i),sum_yy(i),sum_zz(i),sum_xy(i),sum_xz(i),sum_yz(i),Arot,&
                &p(1),p(2),p(3),p(4),p(5),p(6))
        sum_gn(i)=p(3)+sigdot(i_)
        sum_gs(i)=p(5)+taudot(i_)
        sum_gd(i)=p(6)+tauddot(i_)
      end do

      !no dip slip allowed
      dtauddt=0d0
      !call deriv_d(sum_gs,sum_gn,phitmp,taustmp,sigmatmp,veltmp,dphidt,dtausdt,dsigdt)
      !dsigdt=0d0
      !slip rate is parallel to shear traction
      call deriv_3dn(sum_gs,sum_gd,sum_gn,phitmp,taustmp,taudtmp,tautmp,sigmatmp,veltmp,dphidt,dtausdt,dtauddt,dsigdt)

      do i = 1, NCELL
        dydx(4*i-3) = dphidt(i)
        dydx(4*i-2) = dtausdt(i)
        dydx(4*i-1) = dtauddt(i)
        dydx(4*i) = dsigdt(i)

        if(limitsigma.and.(sigmatmp(i).lt.30d0)) dsigdt=0d0
        if(limitsigma.and.(sigmatmp(i).gt.170d0)) dsigdt=0d0
      enddo
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


    ! select case(law)
    !   ! aing law (Linker & Dieterich, 1992)
    ! case('a')
    !   call deriv_a(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
    !   ! slip law
    ! case('s')
    !   call deriv_s(sum_gs,sum_gn,veltmp,tautmp,sigmatmp,dlnvdt,dtaudt,dsigdt)
    !   ! RFL in FDMAP (Dunham+ 2011)
    ! case('d')
    !   call deriv_d(sum_gs,sum_gn,phitmp,tautmp,sigmatmp,veltmp,dphidt,dtaudt,dsigdt)
    !
    ! end select

    ! select case(problem)
    ! case('2dp','3dp')
    !   do i = 1, NCELL
    !     dydx(2*i-1) = dphidt(i)
    !     dydx(2*i) = dtaudt(i)
    !   enddo
    ! case('2dn')
    !   do i = 1, NCELL
    !     dydx(3*i-2) = dphidt(i)
    !     dydx(3*i-1) = dtaudt(i)
    !     dydx(3*i) = dsigdt(i)
    !   enddo
    ! case('3dn','3dh')
    !   do i = 1, NCELL
    !     dydx(4*i-3) = dphidt(i)
    !     dydx(4*i-2) = dtausdt(i)
    !     dydx(4*i-1) = dtauddt(i)
    !     dydx(4*i) = dsigdt(i)
    !   enddo
    ! end select

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
    integer::i,i_
    real(8)::fss,dvdtau,dvdsig,dvdphi
    !type(t_deriv),intent(in) ::
    real(8),intent(in)::sum_gs(:),sum_gn(:),phitmp(:),tautmp(:),sigmatmp(:),veltmp(:)
    real(8),intent(out)::dphidt(:),dtaudt(:),dsigdt(:)
    do i=1,size(sum_gs)
      i_=vars(i)
      !write(*,*) 'vel',veltmp(i)
      dsigdt(i)=sum_gn(i)
      !write(*,*) 'dsigdt',dsigdt(i)

      !fss=mu0+(a(i)-b(i))*dlog(abs(veltmp(i))/vref)
      !fss=fw(i)+(fss-fw(i))/(1.d0+(veltmp(i)/vw(i))**8)**0.125d0 !flash heating
      !regularized slip law
      !dphidt(i)=-abs(veltmp(i))/dc(i)*(abs(tautmp(i))/sigmatmp(i)-fss)
      !regularized aing law
      dphidt(i)=b(i_)*vref/dc(i_)*dexp((mu0-phitmp(i))/b(i_))-abs(veltmp(i))/dc(i_)

      dvdtau=2*vref*dexp(-phitmp(i)/a(i_))*dcosh(tautmp(i)/sigmatmp(i)/a(i_))/(a(i_)*sigmatmp(i))
      dvdsig=-2*vref*dexp(-phitmp(i)/a(i_))*dcosh(tautmp(i)/sigmatmp(i)/a(i_))*tautmp(i)/(a(i_)*sigmatmp(i)**2)
      !dvdphi=2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
      dvdphi=-veltmp(i)/a(i_)
      !dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*phitmp(i)*dvdsig*sigmatmp(i))
      dtaudt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))
      dtaudt(i)=dtaudt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      !write(*,*) rigid/vs*dvdtau
      if(veltmp(i).le.0d0) then
        dvdtau=2*vref*dexp(-phitmp(i)/a(i_))*dcosh(tautmp(i)/sigmatmp(i)/a(i_))/(a(i_)*sigmatmp(i))
        dvdsig=-2*vref*dexp(-phitmp(i)/a(i_))*dcosh(tautmp(i)/sigmatmp(i)/a(i_))*tautmp(i)/(a(i_)*sigmatmp(i)**2)
        !sign ok?
        !dvdphi=2*vref*exp(-phitmp(i)/a(i))*sinh(tautmp(i)/sigmatmp(i)/a(i))/a(i)
        dvdphi=-veltmp(i)/a(i_)
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
      do i=1,size(phitmp)
        !write(*,*) 'vel',veltmp(i)
        dsigdt(i)=sum_gn(i)
        !fss=mu0+(a(i)-b(i))*dlog(abs(veltmp(i))/vref)
        !fss=fw(i)+(fss-fw(i))/(1.d0+(veltmp(i)/vw(i))**8)**0.125d0 !flash heating
        !slip law
        !dphidt(i)=-abs(veltmp(i))/dc(i)*(abs(tautmp(i))/sigmatmp(i)-fss)
        !aing law
        dphidt(i)=b(i_)*vref/dc(i_)*exp((mu0-phitmp(i))/b(i_))-veltmp(i)/dc(i_)
        dvdtau=2*vref*dexp(-phitmp(i)/a(i_))*dcosh(tautmp(i)/sigmatmp(i)/a(i_))/(a(i_)*sigmatmp(i))
        dvdsig=-2*vref*dexp(-phitmp(i)/a(i_))*dcosh(tautmp(i)/sigmatmp(i)/a(i_))*tautmp(i)/(a(i_)*sigmatmp(i)**2)
        dvdphi=-veltmp(i)/a(i_)
        dtausdt(i)=sum_gs(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))*(taustmp(i)/tautmp(i))
        dtausdt(i)=dtausdt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
        dtauddt(i)=sum_gd(i)-0.5d0*rigid/vs*(dvdphi*dphidt(i)+dvdsig*dsigdt(i))*(taudtmp(i)/tautmp(i))
        dtauddt(i)=dtauddt(i)/(1d0+0.5d0*rigid/vs*dvdtau)
      end do
  end subroutine

  !---------------------------------------------------------------------
  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,errmax_gb)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
    !---------------------------------------------------------------------
    use m_HACApK_solve
    use m_HACApK_base
    use m_HACApK_use
    implicit none
    include 'mpif.h'
    !integer::NCELL,NCELLg,rcounts(:),displs(:)
    real(8),intent(in)::yscal(:),htry,eps
    real(8),intent(inout)::y(:),x,dydx(:)
    real(8),intent(out)::hdid,hnext,errmax_gb !hdid: resulatant dt hnext: htry for the next
    !type(st_HACApK_lcontrol),intent(in) :: st_ctl
    !type(st_HACApK_leafmtxp),intent(in) :: st_leafmtxp
    !type(st_HACApK_calc_entry) :: st_bemv
    integer :: i,ierr
    real(8) :: errmax,h,xnew,htemp
    real(8),dimension(size(y))::yerr,ytemp
    real(8),parameter::SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,ERRCON=1.89d-4

    h=htry
    do while(.true.)
      call rkck(y,dydx,x,h,ytemp,yerr)!,,st_leafmtxp,st_bemv,st_ctl)!,derivs)
      errmax=0d0
      select case(problem)
      case('3dn','3dh')
      do i=1,NCELL
        if(abs(yerr(4*i-3)/yscal(4*i-3))/eps.gt.errmax) errmax=abs(yerr(4*i-3)/yscal(4*i-3))/eps
      end do
    case('2dp','3dp','2dn','2dn3')
      errmax=maxval(abs(yerr(:)/yscal(:)))/eps
      end select
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(errmax,errmax_gb,1,MPI_REAL8,                  &
      &     MPI_MAX,MPI_COMM_WORLD,ierr)

      if((errmax_gb.lt.1.d0).and.(errmax_gb.gt.1d-15)) then
        exit
      end if

      !htemp=SAFETY*h*(errmax_gb**PSHRNK)
      h=0.33d0*h
      !h=sign(max(abs(htemp),0.1*abs(h)),h)
      xnew=x+h
      !if(xnew-x<1.d-8) stop
    end do

    hnext=min(1.5*h,SAFETY*h*(errmax_gb**PGROW),1d8)

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
