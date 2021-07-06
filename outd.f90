program main
  implicit none
  integer::job,nrank,rank,ios,ncell,i,ncellg,nstep,n,j,file_size
  real(8)::x(200000),y(200000),z(200000)
  real(8),allocatable::data1(:),vel(:,:)
  character(128)::fname
  job=2
  nrank=60
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
    if(rank==0) allocate(vel(nstep,200000))
    do n=1,nstep
    do i=1,ncell
      vel(n,ncellg+i)=data1(i+ncell*(n-1))
      !write(*,*)n,ncellg+i,vel(n,ncellg+1)
    end do
    end do
    ncellg=ncellg+ncell
    deallocate(data1)

  end do
  write(*,*)ncellg
  write(fname,'("output/vel",i0,".dat")') job
  open(32,file=fname,status='replace')
  do n=1,nstep
    do i=1,ncellg
      write(32,'(4e15.6)') x(i),y(i),z(i),vel(n,i)
    end do
    write(32,*)
    write(32,*)
  end do
  stop
end program
