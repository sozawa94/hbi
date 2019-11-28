module mod_iofdmap
  use mod_constant

contains
  subroutine input_from_FDMAP(number)
    real(8),intent(in)::
    real(8),intent(out)::tauG(:),sigmaG(:),velG(:)

    k=number/10
    iter=number-k*10
    if(iter.eq.1) then
      !Uniform initial stress field
      syy0=-sigma0
      sxy0=-syy0*mu0
      sxx0=syy0*(1d0+2*sxy0/(syy0*dtan(2*phi/180d0*pi)))
      szz0=0.5d0*(sxx0+syy0)
      write(*,*) 'sxx0,sxy0,syy0,szz0',sxx0,sxy0,syy0,szz0
      do i=1,NCELLg
        tauG(i)=sxy0*cos(2*ang(i))+0.5d0*(sxx0-syy0)*sin(2*ang(i))
        sigmaG(i)=sin(ang(i))**2*sxx0+cos(ang(i))**2*syy0+sxy0*sin(2*ang(i))
        velG(i)=vpl!+velinit*1d4*exp(-(i-NCELLg/2)**2/(25d0)**2)
      end do
    else
      !input data from FDMAP's output
      !number=number-1
      write(fname,'("../fdmap/data/",i0,"_I_V.dat")') number-1
      open(25,file=fname,access='stream')
      inquire(25, size=file_size)
      write(*, '(a, i12)') "file size = ", file_size
      n=file_size/8
      allocate(fdmapout(n))
      read(25) fdmapout
      velG=fdmapout(n-NCELLg+1:n)

      !write(*,*) velG
      !deallocate(fdmapout)
      close(25)

      write(fname,'("../fdmap/data/",i0,"_I_S.dat")') number-1!,iter
      open(25,file=fname,access='stream')
      read(25) fdmapout
      tauG(1:NCELLg)=fdmapout(n-NCELLg+1:n)
      close(25)

      write(fname,'("../fdmap/data/",i0,"_I_N.dat")') number-1
      open(25,file=fname,access='stream')
      read(25) fdmapout
      sigmaG(1:NCELLg)=fdmapout(n-NCELLg+1:n)
      close(25)
    end if

  end subroutine

  subroutine output_to_FDMAP(number,coordinate_file,dispG,tauG,sigmaG,psiG,stime,sxx0,sxy0,syy0,szz0,ag,bg,dcg,)
    real(8),intent(in)::dispG(:),tauG(:),sigmaG(:),psiG(:),stime
    integer,intent(in)::number
    character(128),intent(in)::coordinate_file
    integer::
    real(8)::Bx(),By()

    k=number/10
    iter=number-k*10

    !initial stress & state files for fdmap
    open(99,file='../fdmap/problems/prestress',access='stream')
    write(99) tauG,sigmaG
    close(99)

    open(77,file='../fdmap/problems/prestate',access='stream')
    write(77) psiG
    close(77)
    !write(*,*) Psi


    !read initial body stress
    if(iter.eq.1) then
      sxx=-sxx0
      sxy=sxy0
      syy=-syy0
      szz=-szz0

    else
      !getting time step information
      write(fname,'("../fdmap/data/",i0,"_I_V.dat")') number-1
      open(25,file=fname,access='stream')
      inquire(25, size=file_size)
      write(*, '(a, i12)') "file size = ", file_size
      n=file_size/8
      ts=nbody/NCELLg

      !stress change during dynamic rupture
      write(fname,'("data/",i0,"F.ckpt",i0)') number-1,ts
      open(41,file=fname,access='stream')
      inquire(41, size=file_size)
      n=file_size/8
      write(*,*) n
      allocate(fdmap(n))
      read(41) fdmap
      sxx(1:k)=fdmap(2*k+1:3*k)
      sxy(1:k)=fdmap(3*k+1:4*k)
      syy(1:k)=fdmap(4*k+1:5*k)
      szz(1:k)=fdmap(5*k+1:6*k)
      close(41)
      deallocate(fdmap)

      !offset added to above
      write(fname,'("../fdmap/data/",i0,"_B_sxx0.dat")') number-1
      open(41,file=fname,access='stream')
      !open(41,file='data/try_B_sxx0.dat',access='stream')
      inquire(41, size=file_size)
      n=file_size/8
      allocate(fdmap(n))
      read(41) dsxx
      !  sxx(1:k)=fdmap(n-k+1:n)
      close(41)

      write(fname,'("../fdmap/data/",i0,"_B_sxy0.dat")') number-1
      open(41,file=fname,access='stream')
      !open(41,file='data/try_B_sxy0.dat',access='stream')
      read(41) dsxy
      !sxy(1:k)=fdmap(n-k+1:n)
      close(41)

      write(fname,'("../fdmap/data/",i0,"_B_syy0.dat")') number-1
      open(41,file=fname,access='stream')
      !open(41,file='data/try_B_syy0.dat',access='stream')
      read(41) dsyy
      !syy(1:k)=fdmap(n-k+1:n)
      close(41)

      write(fname,'("../fdmap/data/",i0,"_B_szz0.dat")') number-1
      open(41,file=fname,access='stream')
      !open(41,file='data/try_B_szz0.dat',access='stream')
      read(41) dszz
      ! szz(1:k)=fdmap(n-k+1:n)
      close(41)

      deallocate(fdmap)
    end if

    !stress change during interseismic period
    !read B_x,B_y coordinates
    open(41,file=coordinate_file,access='stream')
    read(41) Bx
    read(41) By
    do i=1,nbody
      xloc=Bx(i)
      yloc=By(i)
      call calc_is(xloc,yloc,xel,xer,yel,yer,ang,stime,disp,qsxx,qsxy,qsyy)
      sxx(i)=sxx(i)+dsxx(i)+qsxx
      sxy(i)=sxy(i)+dsxy(i)+qsxy
      syy(i)=syy(i)+dsyy(i)+qsyy
      szz(i)=szz(i)+dszz(i)
    end do

    !write to FDMAP input
    open(43,file='../fdmap/problems/prestress_body',access='stream',status='replace')
    write(43) sxx
    write(43) sxy
    write(43) syy
    write(43) szz
    close(43)

    !save frictinal aprameters
    fv0=vref
    ff0=mu0
    open(32,file='../fdmap/problems/friction',access='stream')
    write(32) ag,bg,fv0,ff0,dcg,fwg,vwg
    close(32)

  end subroutine

  subroutine calc_is(xloc,yloc,xel,xer,yel,yer,ang,stime,disp,qsxx,qsxy,qsyy)
    real(8),intent(in)::xloc,yloc,xel(:),xer(:),yel(:),yer(:),ang(:),disp(:),stime
    real(8),intent(out)::qsxx,qsxy,qsyy
    real(8)::sxx_inc,sxy_inc,syy_inc

    qsxx=0d0;qsxy=0d0;qsyy=0d0
    xsl=-11*edge
    xsr=-edge
    ysl=ycol(1)
    ysr=ycol(1)
    angs=0d0
    slip=stime*vpl
    call TY1997B11(xsl,xsr,ysl,ysr,angs,slip,sxx_inc,sxy_inc,syy_inc)
    qsxx=qsxx+sxx_inc;qsxy=qsxy+sxy_inc;qsyy=qsyy+sxy_inc

    xsl=edge
    xsr=11*edge
    ysl=ycol(NCELLg)
    ysr=ycol(NCELLg)
    call TY1997B11(xsl,xsr,ysl,ysr,angs,slip,sxx_inc,sxy_inc,syy_inc)
    qsxx=qsxx+sxx_inc;qsxy=qsxy+sxy_inc;qsyy=qsyy+sxy_inc

    do k=1,NCELLg
      xsl=xel(k)
      xsr=xer(k)
      ysl=yel(k)
      ysr=yer(k)
      angs=ang(k)
      slip=disp(k)
      call TY1997B11(xsl,xsr,ysl,ysr,angs,slip,sxx_inc,sxy_inc,syy_inc)
      qsxx=qsxx+sxx_inc;qsxy=qsxy+sxy_inc;qsyy=qsyy+sxy_inc
    end do

  end subroutine
  subroutine TY1997B11(xsl,xsr,ysl,ysr,angs,slip,sxx_inc,sxy_inc,syy_inc)
    real(8),intent(in)::
    real(8),intent(out)::
    real(8)::dx(2),dy(2),r(2),gx(2),gy(2),angs,nx,ny,sxx(2),sxy(2),syy(2)

    factor=rigid/(2.d0*pi*(1.d0-pois))
    !dx(1)=xcol(i)-xe(j-1)
    dx(1)=xloc-xsl
    !dy(1)=ycol(i)-ye(j-1)
    dy(1)=yloc-ysl
    r(1)=sqrt(dx(1)**2+dy(1)**2)
    gx(1)=dx(1)/r(1)
    gy(1)=dy(1)/r(1)

    dx(2)=xc-xsr(k)
    dy(2)=yc-ysr(k)
    r(2)=sqrt(dx(2)**2+dy(2)**2)
    gx(2)=dx(2)/r(2)
    gy(2)=dy(2)/r(2)

    nx=-sin(angs)
    ny=cos(angs)

    sxx(1)=(nx*gx(1)+ny*gy(1)+2*gx(1)*gy(1)*(-ny*gx(1)+nx*gy(1)))/r(1)
    sxx(2)=(nx*gx(2)+ny*gy(2)+2*gx(2)*gy(2)*(-ny*gx(2)+nx*gy(2)))/r(2)
    sxy(1)=(-gx(1)**2+gy(1)**2)*(-ny*gx(1)+nx*gy(1))/r(1)
    sxy(2)=(-gx(2)**2+gy(2)**2)*(-ny*gx(2)+nx*gy(2))/r(2)
    syy(1)=(nx*gx(1)+ny*gy(1)-2*gx(1)*gy(1)*(-ny*gx(1)+nx*gy(1)))/r(1)
    syy(2)=(nx*gx(2)+ny*gy(2)-2*gx(2)*gy(2)*(-ny*gx(2)+nx*gy(2)))/r(2)
    sxx_inc=(sxx(2)-sxx(1))*slip*factor
    sxy_inc=(sxy(2)-sxy(1))*slip*factor
    syy_inc=(syy(2)-syy(1))*slip*factor
  end subroutine
end module
