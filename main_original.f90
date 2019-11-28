program main
  !quasi-dynamic BIEM on a 3D planar fault in homogeneous & infinite medium
  !each cell is rectangular
  !MPI calculation
  !H-matrix approximation (using HACApK)
  !Rate and State Friction
  !$ use omp_lib
  use m_HACApK_solve
  use m_HACApK_base
  use m_HACApK_use
  use mod_derivs
  use dtriangular
  implicit none
  include 'mpif.h'
  integer::NCELL, nstep1, lp, i,j,k,m,counts,interval,number,lrtrn,nl,NCELLg,loc
  integer::clock,cr,counts2,imax,jmax,NCELLm,seedsize,icomm,np,ierr,my_rank,load
  integer,allocatable::seed(:)
  character*128::fname,dum,law
  real(8)::a0,b0,sr,omega,theta,dtau,tiny,x,time1,time2,moment,aslip,avv
  real(8)::rigid,vc0,dz,mu0,pois,dtinit
  real(8)::r(3),eps,vpl,outv,vs
  real(8)::dtime,dtnxt,dttry,dtdid,alpha,ds,amp,mui,strinit,velinit,velmax
  real(8),parameter::pi=4.d0*atan(1.d0),dc0=1.d0,sigma0=1.0d0,vref=1.d0
  type(st_HACApK_lcontrol) :: st_ctl
  type(st_HACApK_leafmtxp) :: st_leafmtxp
  type(st_HACApK_calc_entry) :: st_bemv
  type(t_deriv)::st_deriv
  real(8),allocatable ::coord(:,:),vmax(:)
  real(8),allocatable::ag(:),bg(:),dcg(:),vcg(:),taudot(:)
  real(8),allocatable::a(:),b(:),dc(:),sigma(:),vc(:)
  real(8),allocatable::vel(:),tau(:),disp(:),mu(:),s(:)
  real(8),allocatable::velG(:),dispG(:),muG(:)
  real(8),allocatable::xcol(:),zcol(:)
  real(8),allocatable::xs1(:),xs2(:),xs3(:),xs4(:),zs1(:),zs2(:),zs3(:),zs4(:)
  real(8),allocatable::y(:),yscal(:),dydx(:),xcoll(:),zcoll(:)
  real(8),allocatable::xr(:),yr(:),zr(:),ai(:,:),xc(:),zc(:)
  integer::r1,r2,r3,NVER,amari,out,kmax
  integer,allocatable::displs(:),rcounts(:),vmaxin(:)

  icomm=MPI_COMM_WORLD
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )
  if(my_rank.eq.0) time1=MPI_Wtime()

  !reading input file
  open(33,file='input.dat')
  !DOF, fault geometry
  read(33,*) dum,imax !x length
  read(33,*) dum,jmax !z length
  read(33,*) dum,ds !mesh interval(normalized by Dc) ><50: slip law 2500
  !read(33,*) dum,bc !boundary condition. p: periodic, r:rigid boundary
  !read(33,*) dum,kmax

  !output control
  read(33,*) dum,number !output filename
  read(33,*) dum,interval !output frequency
  read(33,*) dum,loc !output localdata

  !continue or stop control
  read(33,*) dum,nstep1 !maxmimum time step
  read(33,*) dum,velmax !stop when slip rate exceeds this value

  !physical parameters in calculation
  read(33,*) dum,law ! evolution law in RSF a: aging s: slip d:dynamic weakening
  read(33,*) dum,a0 !a in RSF
  read(33,*) dum,b0 !b in RSF
  read(33,*) dum,vc0 !cut-off velocity in RSF
  read(33,*) dum,mu0 !mu0 in RSF
  read(33,*) dum,load !loading type 0: uniform stress 1: uniform slip deficit
  read(33,*) dum,sr !if load=0: stressing rate load=1: unused
  read(33,*) dum,vpl !if load=1: loading velocity load=0: unused
  read(33,*) dum,rigid !shear modulus normalized by 100MPa
  read(33,*) dum,vs !Swave speed (dc/s)
  read(33,*) dum,pois !poisson ratio

  !initial values
  !read(33,*) dum,strinit !initial applied stress(normalized)
  read(33,*) dum,velinit !initial slip velocity
  read(33,*) dum,omega !initial omega=V*theta
  read(33,*) dum,dtinit !initial time step

  read(33,*) dum,eps !error allowance
  !read(*,*) amp
  !read(*,*) omega
  close(33)

  !MPI setting
  !NCELLg=2*NL*NL
  NCELLg=imax*jmax

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

  allocate(a(NCELL),b(NCELL),dc(NCELL),vc(NCELL),sigma(NCELL),taudot(NCELL))
  allocate(ag(NCELLg),bg(NCELLg),dcg(NCELLg),vcg(NCELLg))
  allocate(vel(NCELL),tau(NCELL),mu(NCELL),s(NCELL),disp(NCELL))
  allocate(velG(NCELLg),dispG(NCELLg),muG(NCELLg))
  allocate(xcol(NCELLg),zcol(NCELLg))
  allocate(xs1(NCELLg),xs2(NCELLg),xs3(NCELLg),xs4(NCELLg))
  allocate(zs1(NCELLg),zs2(NCELLg),zs3(NCELLg),zs4(NCELLg))
  allocate(xcoll(NCELL),zcoll(NCELL))
  allocate(y(2*NCELL),yscal(2*NCELL),dydx(2 * NCELL))
  allocate(vmax(NCELLg),vmaxin(NcELLg))
  !dr=dr*ds
  !mesh generation
  if(my_rank.eq.0) then
    call coordinate(imax,jmax,ds,xcol,zcol,xs1,xs2,xs3,xs4,zs1,zs2,zs3,zs4)
    write(*,*) 'mesh generated'

    !random number seed
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do i = 1, seedsize
      call system_clock(count=seed(i))
    end do
    call random_seed(put=seed(:))

    !frictional parameters and stressing rates
    call params(NCELLg,a0,b0,dc0,vc0,ag,bg,dcg,vcg)

  end if

  !MPI communication
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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
  !call MPI_BCAST(xc, kmax, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  !call MPI_BCAST(zc, kmax, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERv(xcol,rcounts,displs,MPI_REAL8,xcoll,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(zcol,rcounts,displs,MPI_REAL8,zcoll,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(vcg,rcounts,displs,MPI_REAL8,vc,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(ag,rcounts,displs,MPI_REAL8,a,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(bg,rcounts,displs,MPI_REAL8,b,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTERv(dcg,rcounts,displs,MPI_REAL8,dc,NCELL,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

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

  do i=1,NCELLg
    coord(i,1)=xcol(i)
    coord(i,2)=0.d0
    coord(i,3)=zcol(i)
  end do

  !generate kernel (H-matrix aprrox)
  lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,1d-4)
  if(my_rank.eq.0) write(*,*) 'H-matrix generated'
  !stop

  taudot=1d-10

  !nfmax=NCELL/NCELL+10
  !write(*,*) 'prepared'

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !initial condition
  theta=omega/velinit
  !mui=mu0+a0*dlog(velinit/vref)+b0*dlog(theta)

  do i=1,NCELL
    !call random_number(omega)
    !omega=omega*0.5d0+0.5d0
    theta=omega/velinit
    mui=mu0+a(i)*dlog(velinit/vref)+b(i)*dlog(theta*vref/dc(i)+vref/vc(i))
    sigma(i)=sigma0
    tau(i)=mui*sigma(i)
    vel(i)=velinit
  end do

  !output setting
  if(my_rank.eq.0) then
    write(fname,'("output/monitor",i0,".dat")') number
    open(52,file=fname)
    write(fname,'("output/",i0,".dat")') number
    open(50,file=fname)
  end if

  allocate(st_deriv%a(NCELL),st_deriv%b(NCELL),st_deriv%dc(NCELL))
  allocate(st_deriv%vc(NCELL),st_deriv%sigma(NCELL))

  st_deriv%a=a; st_deriv%b=b; st_deriv%dc=dc;st_deriv%vc=vc
  st_deriv%sigma=sigma
  st_deriv%vref=vref;st_deriv%vpl=vpl;st_deriv%sr=sr
  st_deriv%rigid=rigid;st_deriv%vs=vs;st_deriv%mu0=mu0;st_deriv%pois=pois
  st_deriv%law=law
  st_deriv%load=load

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !call MPI_ALLGATHERv(vc,NCELL,MPI_REAL8,vcG,rcounts,displs,                &
  !&     MPI_REAL8,MPI_COMM_WORLD,ierr)

  !if(my_rank.eq.0) then
  !  write(fname,'("output/initial",i0,".dat")') number
  !  open(19,file=fname)
  !  do i=1,NCELLg
  !    write(19,*) xcol(i),zcol(i),vcg(i)
  !  end do
  !  close(19)
  !  write(*,*) 'start time integration'
  !end if

  !start time integration
  x=0.d0
  disp=0.d0
  dtnxt = dtinit
  time1=MPI_Wtime()
  outv=1d-6

  do k=1,NSTEP1
    dttry = dtnxt
    do i=1,NCELL
      y(2*i-1) = dlog(vel(i))
      y(2*i) = tau(i)
    end do

    call derivs(x, y, dydx,NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)

    do i = 1, 2*NCELL
      yscal(i)=abs(y(i))+abs(dttry*dydx(i))+tiny
    end do

    call rkqs(y,dydx,x,dttry,eps,yscal,dtdid,dtnxt,NCELL,NCELLg,rcounts,displs,st_deriv,st_leafmtxp,st_bemv,st_ctl)

    do i = 1, NCELL
      vel(i) = exp(y(2*i-1))
      tau(i) = y(2*i)
      disp(i) = disp(i)+exp(y(2*i-1))*dtdid
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
      write(52,'(i6,f16.4,4e16.4,f16.4)')k,x,maxval(velG),sum(dispG)/NCELLg,sum(muG)/NCELLg,sum(VelG)/NCELLg,time2-time1
      do i=1,NCELLg
        vmax(i)=max(vmax(i),velG(i))
        vmaxin(i)=k
      end do
      !if(mod(k,interval).eq.0) then

      !output control
      out=1
      !A : iteration number
      if(mod(k,interval).eq.0) out=0

      !B : slip velocity
      !if(maxval(velG).gt.outv) then
      !  out=0
      !  outv=outv*(10.d0)**(0.5d0)
      !end if


      !$ time2=omp_get_wtime()
      if(out.eq.0) then
        write(*,*) 'time step=' ,k
        do i=1,NCELLg
          write(50,'(5e15.6,i7)') xcol(i),zcol(i),log10(velG(i)),muG(i),dispG(i),k
        end do
        write(50,*)
        write(50,*)
      end if
    end if

    if(abs(maxval(velG)).gt.velmax) then
      if(my_rank .eq. 0) write(*,*) 'dlip rate exceeds threshold'
      exit
    end if

    dttry = dtnxt
  end do

  !finalize
  !$ time2=omp_get_wtime()

  ! if(my_rank.eq.0) then
  !   write(fname,'("output/vmax",i0,".dat")') number
  !   open(71,file=fname)
  !   do i=1,NCELLg
  !     write(71,'(3e15.6,i6)') xcol(i),zcol(i),log10(vmax(i)),vmaxin(i)
  !   end do
  !   close(71)
  ! end if
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
end program
