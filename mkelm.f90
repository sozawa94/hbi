program main
  !make meshes for "2dn" of hbi
  implicit none
  integer::i,NCELLg,nm,ns(15),file_size,n,count,j,jmax,seedsize
  real(8)::xel(100000),xer(100000),yel(100000),yer(100000)
  real(8)::xr(40001),yr(40001)
  real(8),parameter::pi=4.d0*atan(1.d0)
  real(8),allocatable::data(:)
  real(8)::ds,xc(15),yc(15),ang(15),r,amp,wid
  character(128)::geofile,type
  integer,allocatable::seed(:)

  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do i = 1, seedsize
    call system_clock(count=seed(i))
  end do
  call random_seed(put=seed(:))

  xel=0d0;xer=0d0;yel=0d0;yer=0d0
  type='multi'

  !parameters
  nm=6000
  jmax=size(ns)
  ds=0.025d0
  amp=0.5
  wid=5.d-2

  !length
  !ns=(/300,300,300,300,300,200,200,200,200,200/)
  ns(1:5)=300
  ns(6:15)=200
  call random_number(xc)
  call random_number(yc)
  yc=(yc-0.5d0)*wid
  !write(*,*) yc
  !xc=(/0.5,0.2,0.8,0.3,0.7/)
  !yc=(/-0.06,0.0,0.03,-0.1,-0.02/)
  !ang=(/0,10,-10,10,-10/)
  call random_number(ang)
  ang=(ang-0.5d0)*10
  !ang=0d0
  ang=ang/180*pi
  xc=(xc-0.1d0)*nm*ds*1.2d0
  yc=yc*nm*ds
  !read rough dataset
  geofile='40001seed17.curve'
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
      yel(i)=amp*yr(i)-i*amp*yr(nm)/nm!+r*0.0001d0
      !call random_number(r)
      yer(i)=amp*yr(i+1)-(i+1)*amp*yr(nm)/nm!+r*0.0001d0
      write(*,*) xel(i),yel(i)
    end do
    write(*,*)
    count=nm
    do j=1,jmax
      do i=1,ns(j)
        xel(count+i)=xc(j)+cos(ang(j))*ds*(i-1)
        xer(count+i)=xc(j)+cos(ang(j))*ds*i
        yel(count+i)=yc(j)+sin(ang(j))*ds*(i-1)+amp*(yr(count+i)-yr(count))!+r*0.0001d0
        !call random_number(r)
        yer(count+i)=yc(j)+sin(ang(j))*ds*i+amp*(yr(count+i+1)-yr(count))!+r*0.0001d0
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
