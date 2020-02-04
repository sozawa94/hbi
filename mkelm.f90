program main
  !make meshes for "2dn" of hbi
  implicit none
  integer::i,NCELLg,nm,ns(5),file_size,n,count,j,jmax
  real(8)::xel(100000),xer(100000),yel(100000),yer(100000)
  real(8)::xr(10000),yr(10000)
  real(8),allocatable::data(:)
  real(8)::ds,xc(5),yc(5),r,amp
  character(128)::geofile,type

  xel=0d0;xer=0d0;yel=0d0;yer=0d0
  type='multi'
  nm=3000
  jmax=5
  ds=0.025d0
  amp=10.0
  ns=(/500,300,200,200,200/)
  xc=(/0.5,0.2,0.8,0.3,0.0/)
  yc=(/0.05,-0.02,0.04,-0.03,-0.04/)
  xc=xc*nm*ds
  yc=yc*nm*ds
  !read rough dataset
  geofile='alpha0.001Lmin2N5001seed1.curve'
  open(20,file=geofile,access='stream')
  inquire(20, size=file_size)
  n=file_size/8
  allocate(data(n))
  read(20) data
  close(20)
  xr(1:n/4)=data(1:n/4)
  yr(1:n/4)=data(n/4+1:n/2)
  !write(*,*) xr(1:n/4),yr(1:n/4)


  select case(type)
    case('multi')
      do i=1,nm
        xel(i)=ds*(i-1)
        xer(i)=ds*i
        yel(i)=amp*yr(i)!+r*0.0001d0
        !call random_number(r)
        yer(i)=amp*yr(i+1)!+r*0.0001d0
        write(*,*) xel(i),yel(i)
      end do
      write(*,*)
      count=nm
      do j=1,jmax
      !xc(j)=0.5*nm*ds
      !yc(j)=0.05*nm*ds
      !ns(j)=500
      do i=1,ns(j)
        xel(count+i)=xc(j)+ds*(i-1)
        xer(count+i)=xc(j)+ds*i
        yel(count+i)=yc(j)+amp*yr(nm+i)!+r*0.0001d0
        !call random_number(r)
        yer(count+i)=yc(j)+amp*yr(nm+i+1)!+r*0.0001d0
        write(*,*) xel(count+i),yel(count+i)
      end do
      write(*,*)
      count=count+ns(j)
      end do
  end select
  NCELLg=count
  !output
  geofile='geotmp'
  !write(*,*) "No. of elements=", NCELLg
  open(20,file=geofile,access='stream')
  write(20) xel(1:NCELLg),xer(1:NCELLg),yel(1:NCELLg),yer(1:NCELLg)
  close(20)
  stop
end program
