program main
  implicit none
  integer::sp,job,nrank,rank,ios,ncell,i,ncellg,nstep,n,j,file_size,stat
  integer,parameter::nmax=500000
  real(8)::x(nmax),y(nmax),z(nmax)
  real(8),allocatable::data1(:),vel(:,:),sigma(:,:),tau(:,:),slip(:,:)
  character(128)::fname,tmp1,tmp2
  !read(*,*)job
  !read(*,*)nrank
  call get_command_argument(1,tmp1,status=stat)
  !write(*,*)tmp1
  read(tmp1,*) job
  call get_command_argument(2,tmp2,status=stat)
  !write(*,*)tmp2
  read(tmp2,*) nrank
  !write(*,*) job,nrank

  sp=1
  !nrank=64
  ncellg=0
  do rank=0,nrank-1
    write(fname,'("output/xyz",i0,"_",i0,".dat")') job,rank
    open(25,file=fname,iostat=ios)
    ncell=0
    do while(ios==0)
      read(25,*,iostat=ios) x(ncellg+ncell+1),y(ncellg+ncell+1),z(ncellg+ncell+1)
      ncell=ncell+1
    end do
    ncell=ncell-1
    close(25)

    write(fname,'("output/vel",i0,"_",i0,".dat")') job,rank
    open(25,file=fname,access='stream')
 !open(25,file="output/tau44.dat",access='stream')
    inquire(25, size=file_size)
    j=file_size/8
    allocate(data1(j))
    read(25) data1
    close(25)
    nstep=j/ncell
    write(*,*)ncell,nstep,ncellg
    if(rank==0) then
      allocate(vel(nstep,nmax))
      allocate(sigma(nstep,nmax))
      allocate(tau(nstep,nmax))
      allocate(slip(nstep,nmax))
    end if

    do n=1,nstep
      do i=1,ncell
        vel(n,ncellg+i)=data1(i+ncell*(n-1))
      end do
    end do

    write(fname,'("output/sigma",i0,"_",i0,".dat")') job,rank
    open(25,file=fname,access='stream')
    read(25) data1
    close(25)
    do n=1,nstep
      do i=1,ncell
        sigma(n,ncellg+i)=data1(i+ncell*(n-1))
      end do
    end do

    write(fname,'("output/tau",i0,"_",i0,".dat")') job,rank
    open(25,file=fname,access='stream')
    read(25) data1
    close(25)
    do n=1,nstep
      do i=1,ncell
        tau(n,ncellg+i)=data1(i+ncell*(n-1))
      end do
    end do

    write(fname,'("output/slip",i0,"_",i0,".dat")') job,rank
    open(25,file=fname,access='stream')
    read(25) data1
    close(25)
    nstep=j/ncell
    do n=1,nstep
      do i=1,ncell
        slip(n,ncellg+i)=data1(i+ncell*(n-1))
      end do
    end do
    deallocate(data1)

    ncellg=ncellg+ncell

  end do

  write(*,*)ncellg
  write(fname,'("output/field",i0,".dat")') job
  open(32,file=fname,status='replace')
  do n=1,nstep
    do i=1,ncellg
      if(mod(i,sp)==0) write(32,'(7e15.6)') x(i),y(i),z(i),vel(n,i),sigma(n,i),tau(n,i),slip(n,i)
    end do
    write(32,*)
    write(32,*)
  end do
  close(32)

  stop
end program
