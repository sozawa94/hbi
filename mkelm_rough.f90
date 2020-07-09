program main
  !make meshes for "2dn" of hbi
  implicit none
  integer::i,NCELLg,nm,file_size,n,count,j,jmax,seedsize,k
  real(8)::xel(10000),xer(10000),yel(10000),yer(10000)
  real(8)::xr(160001),yr(160001)
  real(8),parameter::pi=4.d0*atan(1.d0)
  real(8),allocatable::data(:),ang(:),xg(:),yg(:),leng(:)
  integer,allocatable::ns(:)
  real(8)::ds,r,amp,wid,xc,yc,rr,yf,angle
  real(8)::x1,x2,y1,y2,dx1,dx2,dxa,t1,dy1,dy2,dya,t2
  character(128)::geofile,type
  integer,allocatable::seed(:)

  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  do i = 1, seedsize
    call system_clock(count=seed(i))
  end do
  call random_seed(put=seed(:))

  xel=0d0;xer=0d0;yel=0d0;yer=0d0
  ncellg=10000
  nm=ncellg
  !xc=(xc-0.1d0)*nm*ds*1.2d0
  !yc=yc*nm*ds
  !read rough dataset
  geofile='alpha0.005Lmin0.5N40001seed1.curve'
  open(20,file=geofile,access='stream')
  inquire(20, size=file_size)
  n=file_size/8
  allocate(data(n))
  read(20) data
  close(20)
  xr(1:n/4)=data(1:n/4)
  yr(1:n/4)=data(n/4+1:n/2)
  !write(*,*) xr(1:n/4),yr(1:n/4)
    do i=1,nm
      xel(i)=ds*(i-1)
      xer(i)=ds*i
      yel(i)=amp*yr(4*i)-i*amp*yr(4*nm)/4/nm!-!+r*0.0001d0
      !call random_number(r)
      yer(i)=amp*yr(4*i+4)-(i+1)*amp*yr(4*nm)/4/nm!+r*0.0001d0
    end do

  !output
  geofile='geos/geotmp'
  !write(*,*) "No. of elements=", NCELLg
  open(20,file=geofile,access='stream',status='replace')
  write(20) xel(1:NCELLg),xer(1:NCELLg),yel(1:NCELLg),yer(1:NCELLg)
  close(20)
  stop
end program
