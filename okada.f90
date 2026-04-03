module mod_okada
  implicit none
  !real(8),parameter::pi2=8.d0*datan(1.d0)
  !real(8),parameter::epso=1d-6
  !real(8),private::xi2,et2,q2,r,r2,r3,r5,d,yp,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
  !real(8),private::sd,cd,sdsd,cdcd,sdcd,s2d,c2d
  !real(8)::p,q,s,t,xy,yp,x2,y2,d2,qr,qrx,a3,a5,b3,c3,fuy,vy,wy,fuz,vz,wz

contains
  subroutine okada(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,disl1,disl2,disl3,&
    &ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,fullspace)
    implicit none
    real(8),intent(in)::alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,disl1,disl2,disl3
    logical,intent(in)::fullspace
    real(8),intent(out)::ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    integer::i,j,k
    real(8),parameter::pi2=8.d0*datan(1.d0),epso=1d-6
    real(8)::p18,ddip,aalpha,zz,d,dd1,dd2,dd3,r12,r21,r22,p,q
    real(8):: dummy(5),sd,cd
    real(8)::  xi(2),et(2)
    integer::kxi(2),ket(2)
    real(8)::  u(12),du(12),dua(12),dub(12),duc(12)
    !real(8)::xi2,et2,q2,r,r2,r3,r5,d,yp,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    !real(8):alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
        ! alp1=(1d0-alpha)/2d0
        ! alp2= alpha/2d0
        ! alp3=(1d0-alpha)/alpha
        ! alp4= 1d0-alpha
        ! alp5= alpha
        !----
        p18=pi2/360.d0
        sd=dsin(dip*p18)
        cd=dcos(dip*p18)
        if(dabs(cd).lt.epso) then
          cd=0d0
          if(sd.gt.0d0) sd= 1d0
          if(sd.lt.0d0) sd=-1d0
        endif
        ! sdsd=sd*sd
        ! cdcd=cd*cd
        ! sdcd=sd*cd
        ! s2d=2d0*sdcd
        ! c2d=cdcd-sdsd
       ! write(29,*)sd,cd

   ! zz=z
   ! dd1=disl1
   ! dd2=disl2
   ! dd3=disl3
    xi(1)=x-al1
    xi(2)=x-al2
    !write(*,*) 'xi',xi
    !if(dabs(xi(1)).lt.epso) xi(1)=0d0
    !if(dabs(xi(2)).lt.epso) xi(2)=0d0
    !=====================================
    !====  real-source contribution  =====
    !=====================================
    d=depth+z
    !write(*,*)sd,cd
    !cd=cos(dip)
    !sd=sin(dip)
    p=y*cd+d*sd
    q=y*sd-d*cd
    et(1)=p-aw1
    et(2)=p-aw2
   ! write(*,*) 'q,et',q,et
   ! if(dabs(q).lt.epso)  q=0d0
   ! if(dabs(et(1)).lt.epso) et(1)=0d0
   ! if(dabs(et(2)).lt.epso) et(2)=0d0
    !-------------------------------
    !---- reject singular case -----
    !-------------------------------
    !---- on fault edge
   ! if(q.eq.0d0 .and.((xi(1)*xi(2).le.0d0 .and. et(1)*et(2).eq.0d0).or.(et(1)*et(2).le.0d0 .and. xi(1)*xi(2).eq.0d0)))  then
    !  iret=1
    !  go to 99
    !endif

    !---- on negative extension of fault edge
    ! kxi(1)=0
    ! kxi(2)=0
    ! ket(1)=0
    ! ket(2)=0
    ! r12=dsqrt(xi(1)*xi(1)+et(2)*et(2)+q*q)
    ! r21=dsqrt(xi(2)*xi(2)+et(1)*et(1)+q*q)
    ! r22=dsqrt(xi(2)*xi(2)+et(2)*et(2)+q*q)
    ! !write(*,*) 'r12,r21,r22',r12,r21,r22
    ! if(xi(1).lt.0d0 .and. r21+xi(2).lt.epso) kxi(1)=1
    ! if(xi(1).lt.0d0 .and. r22+xi(2).lt.epso) kxi(2)=1
    ! if(et(1).lt.0d0 .and. r12+et(2).lt.epso) ket(1)=1
    ! if(et(1).lt.0d0 .and. r22+et(2).lt.epso) ket(2)=1


    !====
    u=0d0
    do k=1,2
      do j=1,2
        !call dccon2(xi(j),et(k),q,sd,cd,kxi(k),ket(j))
       ! call ua(xi(j),et(k),q,dd1,dd2,dd3,dua)
         call ua(xi(j),et(k),q,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,dua)

        !----
        !!write(*,*)'sd,cd',sd,cd
        !do 220 i=1,10,3
        du(1)  =-dua(1)
        du(2)=-dua(2)*cd+dua(3)*sd
        du(3)=-dua(2)*sd-dua(3)*cd
        du(4)  =-dua(4)
        du(5)=-dua(5)*cd+dua(6)*sd
        du(6)=-dua(5)*sd-dua(6)*cd
        du(7)  =-dua(7)
        du(8)=-dua(8)*cd+dua(9)*sd
        du(9)=-dua(8)*sd-dua(9)*cd
        du(10)  =dua(10)
        du(11)=dua(11)*cd-dua(12)*sd
        du(12)=dua(11)*sd+dua(12)*cd
        ! du(10)  =dua(10)
        ! du(11)=dua(11)*cd-dua(12)*sd
        ! du(12)=dua(11)*sd+dua(12)*cd
        ! if(i.lt.10) go to 220
        ! du(i)  =-du(i)
        ! du(i+1)=-du(i+1)
        ! du(i+2)=-du(i+2)
        !end do
        !write(*,*)'du',du
        do i=1,12
          if(j+k.ne.3) u(i)=u(i)+du(i)
          if(j+k.eq.3) u(i)=u(i)-du(i)
        end do
        !----
      end do
    end do
   ux=u(1)
    uy=u(2)
    uz=u(3)
    uxx=u(4)
    uyx=u(5)
    uzx=u(6)
    uxy=u(7)
    uyy=u(8)
    uzy=u(9)
    uxz=u(10)
    uyz=u(11)
    uzz=u(12)
    if(fullspace) return

    !return

    !======================================
    !====  image-source contribution  =====
    !======================================
    !write(*,*)'image-source contribution'
    d=depth-z
    p=y*cd+d*sd
    q=y*sd-d*cd
    et(1)=p-aw1
    et(2)=p-aw2
 !   if(dabs(q).lt.epso)  q=0d0
  !  if(dabs(et(1)).lt.epso) et(1)=0d0
  !  if(dabs(et(2)).lt.epso) et(2)=0d0
    !-------------------------------
    !---- reject singular case -----
    !-------------------------------
    !---- on fault edge
    !if(q.eq.0d0 .and.((xi(1)*xi(2).le.0d0 .and. et(1)*et(2).eq.0d0).or.(et(1)*et(2).le.0d0 .and. xi(1)*xi(2).eq.0d0)))  then
     ! iret=1
     ! go to 99
    !endif
    !---- on negative extension of fault edge
    kxi(1)=0
    kxi(2)=0
    ket(1)=0
    ket(2)=0
    r12=dsqrt(xi(1)*xi(1)+et(2)*et(2)+q*q)
    r21=dsqrt(xi(2)*xi(2)+et(1)*et(1)+q*q)
    r22=dsqrt(xi(2)*xi(2)+et(2)*et(2)+q*q)
 !   if(xi(1).lt.0d0 .and. r21+xi(2).lt.epso) kxi(1)=1
 !   if(xi(1).lt.0d0 .and. r22+xi(2).lt.epso) kxi(2)=1
 !   if(et(1).lt.0d0 .and. r12+et(2).lt.epso) ket(1)=1
 !   if(et(1).lt.0d0 .and. r22+et(2).lt.epso) ket(2)=1
    !====
    do k=1,2
      do j=1,2
        !call dccon2(xi(j),et(k),q,sd,cd,kxi(k),ket(j))
        call ua(xi(j),et(k),q,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,dua)
        call ub(xi(j),et(k),q,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,dub)
        call uc(xi(j),et(k),q,z,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,duc)
        !----
        du(1)=dua(1)+dub(1)+z*duc(1)
        du(2)=(dua(2)+dub(2)+z*duc(2))*cd-(dua(3)+dub(3)+z*duc(3))*sd
        du(3)=(dua(2)+dub(2)-z*duc(2))*sd+(dua(3)+dub(3)-z*duc(3))*cd
        du(4)=dua(4)+dub(4)+z*duc(4)
        du(5)=(dua(5)+dub(5)+z*duc(5))*cd-(dua(6)+dub(6)+z*duc(6))*sd
        du(6)=(dua(5)+dub(5)-z*duc(5))*sd+(dua(6)+dub(6)-z*duc(6))*cd
        du(7)=dua(7)+dub(7)+z*duc(7)
        du(8)=(dua(8)+dub(8)+z*duc(8))*cd-(dua(9)+dub(9)+z*duc(9))*sd
        du(9)=(dua(8)+dub(8)-z*duc(8))*sd+(dua(9)+dub(9)-z*duc(9))*cd
        du(10)=dua(10)+dub(10)+z*duc(10)
        du(11)=(dua(11)+dub(11)+z*duc(11))*cd-(dua(12)+dub(12)+z*duc(12))*sd
        du(12)=(dua(11)+dub(11)-z*duc(11))*sd+(dua(12)+dub(12)-z*duc(12))*cd
        du(10)=du(10)+duc(1)
        du(11)=du(11)+duc(2)*cd-duc(3)*sd
        du(12)=du(12)-duc(2)*sd-duc(3)*cd
        !write(*,*)'du_main',du
        do i=1,12
          if(j+k.ne.3) u(i)=u(i)+du(i)
          if(j+k.eq.3) u(i)=u(i)-du(i)
        end do
        !----
      end do
    end do
    !====
       ux=u(1)
    uy=u(2)
    uz=u(3)
    uxx=u(4)
    uyx=u(5)
    uzx=u(6)
    uxy=u(7)
    uyy=u(8)
    uzy=u(9)
    uxz=u(10)
    uyz=u(11)
    uzz=u(12)
    return
  end subroutine okada
  subroutine  ua(xi,et,q,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,dua)
    implicit none
    real(8),intent(in)::xi,et,q,alpha,disl1,disl2,disl3,sd,cd
    integer,intent(in)::kxi(2),ket(2)
    real(8),intent(out)::dua(12)
    real(8)::du(12),qx,qy,xy,alp1,alp2,alp3,alp4,alp5,rxi,ret
    real(8)::xi2,et2,q2,r,r2,r3,r5,d,yp,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    integer::i
     real(8),parameter::pi2=8.d0*datan(1.d0),epso=1d-6

    !c
    !*******************************************************************
    !****    displacement and strain at depth (part-a)             *****
    !****    due to buried finite fault in a semiinfinite medium   *****
    !*******************************************************************
    !c
    !**** input
    !****   xi,et,q : station coordinates in fault system
    !****   disl1-disl3 : strike-, dip-, tensile-dislocations
    !**** output
    !****   u(12) : displacement and their derivatives
    !c
    !real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
    !real(8)::xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    !----
    alp1=(1d0-alpha)/2d0
    alp2= alpha/2d0
    alp3=(1d0-alpha)/alpha
    alp4= 1d0-alpha
    alp5= alpha

    xi2=xi*xi
    et2=et*et
    q2=q*q
    r2=xi2+et2+q2
    r =dsqrt(r2)
    r3=r *r2
    r5=r3*r2
    yp =et*cd+q*sd
    d =et*sd-q*cd
    tt=datan(xi*et/(q*r))
    rxi=r+xi
    alx=dlog(rxi)
    x11=1d0/(r*rxi)
    x32=(r+rxi)*x11*x11/r
    ret=r+et
    ale=dlog(ret)
    y11=1d0/(r*ret)
    y32=(r+ret)*y11*y11/r
    ey=sd/r-yp*q/r3
    ez=cd/r+d*q/r3
    fy=d/r3+xi2*y32*sd
    fz=yp/r3+xi2*y32*cd
    gy=2d0*x11*sd-yp*q*x32
    gz=2d0*x11*cd+d*q*x32
    hy=d*q*x32+xi*q*y32*sd
    hz=yp*q*x32+xi*q*y32*cd

    do i=1,12
      dua(i)=0d0
    end do
    xy=xi*y11
    qx=q *x11
    qy=q *y11
    ! !write(*,*)'xy,qx,qy'
    ! !write(*,*)xy,qx,qy
    !=====================================
    !====  strike-slip contribution  =====
    !=====================================
      du( 1)=    tt/2d0 +alp2*xi*qy
      du( 2)=           alp2*q/r
      du( 3)= alp1*ale -alp2*q*qy
      du( 4)=-alp1*qy  -alp2*xi2*q*y32
      du( 5)=          -alp2*xi*q/r3
      du( 6)= alp1*xy  +alp2*xi*q2*y32
      du( 7)= alp1*xy*sd        +alp2*xi*fy+d/2d0*x11
      du( 8)=                    alp2*ey
      du( 9)= alp1*(cd/r+qy*sd) -alp2*q*fy
      du(10)= alp1*xy*cd        +alp2*xi*fz+yp/2d0*x11
      du(11)=                    alp2*ez
      du(12)=-alp1*(sd/r-qy*cd) -alp2*q*fz
      !write(*,*)'du',du
      do i=1,12
        dua(i)=dua(i)+disl1/pi2*du(i)
      end do

      du( 1)=           alp2*q/r
      du( 2)=    tt/2d0 +alp2*et*qx
      du( 3)= alp1*alx -alp2*q*qx
      du( 4)=        -alp2*xi*q/r3
      du( 5)= -qy/2d0 -alp2*et*q/r3
      du( 6)= alp1/r +alp2*q2/r3
      du( 7)=                      alp2*ey
      du( 8)= alp1*d*x11+xy/2d0*sd +alp2*et*gy
      du( 9)= alp1*yp*x11          -alp2*q*gy
      du(10)=                      alp2*ez
      du(11)= alp1*yp*x11+xy/2d0*cd +alp2*et*gz
      du(12)=-alp1*d*x11          -alp2*q*gz
      !write(*,*)'du',du
      do i=1,12
        dua(i)=dua(i)+disl2/pi2*du(i)
      end do

    !write(*,*)'u',u
    return
  end subroutine
  subroutine  ub(xi,et,q,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,dub)
    implicit none
    real(8),intent(in)::xi,et,q,disl1,disl2,disl3,alpha,sd,cd
    integer,intent(in)::kxi(2),ket(2)
    real(8),intent(out)::dub(12)
    real(8)::du(12),d11,rd,ai4,aj2,aj5,ai3,ak1,aj3,aj6,rd2,ai1,ai2,ak2,ak4,aj1,aj4,qx,qy
    real(8)::xs,ak3,xy,alp1,alp2,alp3,alp4,alp5,rxi,ret
    real(8)::xi2,et2,q2,r,r2,r3,r5,d,yp,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    integer::i
 real(8),parameter::pi2=8.d0*datan(1.d0),epso=1d-6

    !*******************************************************************
    !****    displacement and strain at depth (part-b)             *****
    !****    due to buried finite fault in a semiinfinite medium   *****
    !*******************************************************************
    !c
    !**** input
    !****   xi,et,q : station coordinates in fault system
    !****   disl1-disl3 : strike-, dip-, tensile-dislocations
    !**** output
    !****   u(12) : displacement and their derivatives
    !c
    !real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
    !real(8)::xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    !----
    ! !write(*,*)'sd,cd',sd,cd
    alp1=(1d0-alpha)/2d0
    alp2= alpha/2d0
    alp3=(1d0-alpha)/alpha
    alp4= 1d0-alpha
    alp5= alpha

    xi2=xi*xi
    et2=et*et
    q2=q*q
    r2=xi2+et2+q2
    r =dsqrt(r2)
    r3=r *r2
    r5=r3*r2
    yp =et*cd+q*sd
    d =et*sd-q*cd
    tt=datan(xi*et/(q*r))
    rxi=r+xi
    alx=dlog(rxi)
    x11=1d0/(r*rxi)
    x32=(r+rxi)*x11*x11/r
    ret=r+et
    ale=dlog(ret)
    y11=1d0/(r*ret)
    y32=(r+ret)*y11*y11/r
    ey=sd/r-yp*q/r3
    ez=cd/r+d*q/r3
    fy=d/r3+xi2*y32*sd
    fz=yp/r3+xi2*y32*cd
    gy=2d0*x11*sd-yp*q*x32
    gz=2d0*x11*cd+d*q*x32
    hy=d*q*x32+xi*q*y32*sd
    hz=yp*q*x32+xi*q*y32*cd

    rd=r+d
    ! !write(*,*) 'rd',rd
    d11=1d0/(r*rd)
    aj2=xi*yp/rd*d11
    aj5=-(d+yp*yp/rd)*d11
    if(cd.ne.0d0) then
      if(xi.eq.0d0) then
        ai4=0d0
      else
        xs=dsqrt(xi2+q2)
        ai4=1d0/cd/cd*( xi/rd*sd*cd+2d0*datan((et*(xs+q*cd)+xs*(r+xs)*sd)/(xi*(r+xs)*cd)) )
      endif
      ai3=(yp*cd/rd-ale+sd*dlog(rd))/cd/cd
      ak1=xi*(d11-y11*sd)/cd
      ak3=(q*y11-yp*d11)/cd
      aj3=(ak1-aj2*sd)/cd
      aj6=(ak3-aj5*sd)/cd
    else
      rd2=rd*rd
      !!write(*,*)'rd2',rd2
      ai3=(et/rd+yp*q/rd2-ale)/2d0
      ai4=xi*yp/rd2/2d0
      ak1=xi*q/rd*d11
      ak3=sd/rd*(xi2*d11-1d0)
      aj3=-xi/rd2*(q2*d11-1d0/2d0)
      aj6=-yp/rd2*(xi2*d11-1d0/2d0)
    endif
    !----
    xy=xi*y11
    ai1=-xi/rd*cd-ai4*sd
    ai2= dlog(rd)+ai3*sd
    ak2= 1d0/r+ak3*sd
    ak4= xy*cd-ak1*sd
    aj1= aj5*cd-aj6*sd
    aj4=-xy-aj2*cd+aj3*sd
    ! !write(*,*)'d11,rd,r,d,ai4,aj2,aj5,ai3,ak1,aj3,aj6,rd2,ai1,ai2,ak2,ak4,aj1,aj4,qx,qy'
    ! !write(*,*)d11,rd,r,d,ai4,aj2,aj5,ai3,ak1,aj3,aj6,rd2,ai1,ai2,ak2,ak4,aj1,aj4,qx,qy

    !====
    do i=1,12
      dub(i)=0d0
    end do
    qx=q*x11
    qy=q*y11
    !=====================================
    !====  strike-slip contribution  =====
    !=====================================
      du( 1)=-xi*qy-tt -alp3*ai1*sd
      du( 2)=-q/r      +alp3*yp/rd*sd
      du( 3)= q*qy     -alp3*ai2*sd
      du( 4)= xi2*q*y32 -alp3*aj1*sd
      du( 5)= xi*q/r3   -alp3*aj2*sd
      du( 6)=-xi*q2*y32 -alp3*aj3*sd
      du( 7)=-xi*fy-d*x11 +alp3*(xy+aj4)*sd
      du( 8)=-ey          +alp3*(1d0/r+aj5)*sd
      du( 9)= q*fy        -alp3*(qy-aj6)*sd
      du(10)=-xi*fz-yp*x11 +alp3*ak1*sd
      du(11)=-ez          +alp3*yp*d11*sd
      du(12)= q*fz        +alp3*ak2*sd
      !write(*,*)'ub,du',du
      do i=1,12
        dub(i)=dub(i)+disl1/pi2*du(i)
      end do
    !=====================================
    !====    dip-slip contribution   =====
    !=====================================
      du( 1)=-q/r      +alp3*ai3*sd*cd
      du( 2)=-et*qx-tt -alp3*xi/rd*sd*cd
      du( 3)= q*qx     +alp3*ai4*sd*cd
      du( 4)= xi*q/r3     +alp3*aj4*sd*cd
      du( 5)= et*q/r3+qy  +alp3*aj5*sd*cd
      du( 6)=-q2/r3       +alp3*aj6*sd*cd
      du( 7)=-ey          +alp3*aj1*sd*cd
      du( 8)=-et*gy-xy*sd +alp3*aj2*sd*cd
      du( 9)= q*gy        +alp3*aj3*sd*cd
      du(10)=-ez          -alp3*ak3*sd*cd
      du(11)=-et*gz-xy*cd -alp3*xi*d11*sd*cd
      du(12)= q*gz        -alp3*ak4*sd*cd
      !write(*,*)'ub,du',du
      do i=1,12
        dub(i)=dub(i)+disl2/pi2*du(i)
      end do
    return
  end
  subroutine  uc(xi,et,q,z,sd,cd,kxi,ket,disl1,disl2,disl3,alpha,duc)
    implicit none
    real(8),intent(in)::xi,et,q,z,disl1,disl2,disl3,alpha,sd,cd
    integer,intent(in)::kxi(2),ket(2)
    real(8),intent(out)::duc(12)
    real(8)::du(12),c,x53,y53,h,z32,z53,y0,z0,ppy,ppz,qq,qqy,qqz,xy,qx,qy,qr,cqx,cdr,yy0
    real(8)::alp1,alp2,alp3,alp4,alp5,rxi,ret
    real(8)::xi2,et2,q2,r,r2,r3,r5,d,yp,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    integer::i
     real(8),parameter::pi2=8.d0*datan(1.d0),epso=1d-6
!c
    !*******************************************************************
    !****    displacement and strain at depth (part-c)             *****
    !****    due to buried finite fault in a semiinfinite medium   *****
    !*******************************************************************
    !c
    !**** input
    !****   xi,et,q,z   : station coordinates in fault system
    !****   disl1-disl3 : strike-, dip-, tensile-dislocations
    !**** output
    !****   u(12) : displacement and their derivatives
    !c
    !!real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
    !real(8)::xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
    !----
    alp1=(1d0-alpha)/2d0
    alp2= alpha/2d0
    alp3=(1d0-alpha)/alpha
    alp4= 1d0-alpha
    alp5= alpha

    xi2=xi*xi
    et2=et*et
    q2=q*q
    r2=xi2+et2+q2
    r =dsqrt(r2)
    r3=r *r2
    r5=r3*r2
    yp =et*cd+q*sd
    d =et*sd-q*cd
    tt=datan(xi*et/(q*r))
    rxi=r+xi
    alx=dlog(rxi)
    x11=1d0/(r*rxi)
    x32=(r+rxi)*x11*x11/r
    ret=r+et
    ale=dlog(ret)
    y11=1d0/(r*ret)
    y32=(r+ret)*y11*y11/r
    ey=sd/r-yp*q/r3
    ez=cd/r+d*q/r3
    fy=d/r3+xi2*y32*sd
    fz=yp/r3+xi2*y32*cd
    gy=2d0*x11*sd-yp*q*x32
    gz=2d0*x11*cd+d*q*x32
    hy=d*q*x32+xi*q*y32*sd
    hz=yp*q*x32+xi*q*y32*cd

    c=d+z
    x53=(8.d0*r2+9.d0*r*xi+3d0*xi2)*x11*x11*x11/r2
    y53=(8.d0*r2+9.d0*r*et+3d0*et2)*y11*y11*y11/r2
    h=q*cd-z
    z32=sd/r3-h*y32
    z53=3d0*sd/r5-h*y53
    y0=y11-xi2*y32
    z0=z32-xi2*z53
    ppy=cd/r3+q*y32*sd
    ppz=sd/r3-q*y32*cd
    qq=z*y32+z32+z0
    qqy=3d0*c*d/r5-qq*sd
    qqz=3d0*c*yp/r5-qq*cd+q*y32
    xy=xi*y11
    qx=q*x11
    qy=q*y11
    qr=3d0*q/r5
    cqx=c*q*x53
    cdr=(c+d)/r3
    yy0=yp/r3-y0*cd
    !====
    do  i=1,12
      duc(i)=0d0
    end do
    !=====================================
    !====  strike-slip contribution  =====
    !=====================================
      du( 1)= alp4*xy*cd           -alp5*xi*q*z32
      du( 2)= alp4*(cd/r+2d0*qy*sd) -alp5*c*q/r3
      du( 3)= alp4*qy*cd           -alp5*(c*et/r3-z*y11+xi2*z32)
      du( 4)= alp4*y0*cd                  -alp5*q*z0
      du( 5)=-alp4*xi*(cd/r3+2d0*q*y32*sd) +alp5*c*xi*qr
      du( 6)=-alp4*xi*q*y32*cd            +alp5*xi*(3d0*c*et/r5-qq)
      du( 7)=-alp4*xi*ppy*cd    -alp5*xi*qqy
      du( 8)= alp4*2d0*(d/r3-y0*sd)*sd-yp/r3*cd-alp5*(cdr*sd-et/r3-c*yp*qr)
      du( 9)=-alp4*q/r3+yy0*sd  +alp5*(cdr*cd+c*d*qr-(y0*cd+q*z0)*sd)
      du(10)= alp4*xi*ppz*cd    -alp5*xi*qqz
      du(11)= alp4*2d0*(yp/r3-y0*cd)*sd+d/r3*cd -alp5*(cdr*cd+c*d*qr)
      du(12)=         yy0*cd    -alp5*(cdr*sd-c*yp*qr-y0*sd*sd+q*z0*cd)
      do i=1,12
        duc(i)=duc(i)+disl1/pi2*du(i)
      end do
    !=====================================
    !====    dip-slip contribution   =====
    !=====================================
      du( 1)= alp4*cd/r -qy*sd -alp5*c*q/r3
      du( 2)= alp4*yp*x11       -alp5*c*et*q*x32
      du( 3)=     -d*x11-xy*sd -alp5*c*(x11-q2*x32)
      du( 4)=-alp4*xi/r3*cd +alp5*c*xi*qr +xi*q*y32*sd
      du( 5)=-alp4*yp/r3     +alp5*c*et*qr
      du( 6)=    d/r3-y0*sd +alp5*c/r3*(1d0-3d0*q2/r2)
      du( 7)=-alp4*et/r3+y0*sd*sd -alp5*(cdr*sd-c*yp*qr)
      du( 8)= alp4*(x11-yp*yp*x32) -alp5*c*((d+2d0*q*cd)*x32-yp*et*q*x53)
      du( 9)=  xi*ppy*sd+yp*d*x32 +alp5*c*((yp+2d0*q*sd)*x32-yp*q2*x53)
      du(10)=      -q/r3+y0*sd*cd -alp5*(cdr*cd+c*d*qr)
      du(11)= alp4*yp*d*x32       -alp5*c*((yp-2d0*q*sd)*x32+d*et*q*x53)
      du(12)=-xi*ppz*sd+x11-d*d*x32-alp5*c*((d-2d0*q*cd)*x32-d*q2*x53)
      do i=1,12
        duc(i)=duc(i)+disl2/pi2*du(i)
      end do

    !write(*,*)u
    return
  end subroutine
  ! subroutine  dccon2(xi,et,q,sd,cd,kxi,ket)
  !   implicit none
  !   real(8),intent(in)::xi,et,q,sd,cd
  !   integer,intent(in)::kxi,ket
  !   real(8)::rxi,ret
  !   !c
  !   !*********************************************************************
  !   !****   calculate station geometry constants for finite source   *****
  !   !*********************************************************************
  !   !c
  !   !**** input
  !   !****   xi,et,q : station coordinates in fault system
  !   !****   sd,cd   : sin, cos of dip-angle
  !   !****   kxi,ket : kxi=1, ket=1 means r+xi<epso, r+et<epso, respectively
  !   !c
  !   !c### caution ### if xi,et,q are sufficiently small, they are set to zer0
  !   !c
  !   !real(8)::sd,cd
  !   !----
  !  ! if(dabs(xi).lt.epso) xi=0d0
  !  ! if(dabs(et).lt.epso) et=0d0
  !  ! if(dabs( q).lt.epso)  q=0d0
  !   xi2=xi*xi
  !   et2=et*et
  !   q2=q*q
  !   r2=xi2+et2+q2
  !   r =dsqrt(r2)
  !  ! if(r.eq.0d0) return
  !   r3=r *r2
  !   r5=r3*r2
  !   yp =et*cd+q*sd
  !   d =et*sd-q*cd
  !   !----
  !  ! if(q.eq.0d0) then
  !  !   tt=0d0
  !  ! else
  !     tt=datan(xi*et/(q*r))
  !  ! endif
  !   !----
  !  ! if(kxi.eq.1) then
  !  !   alx=-dlog(r-xi)
  !  !   x11=0d0
  !  !   x32=0d0
  !  ! else
  !     rxi=r+xi
  !     alx=dlog(rxi)
  !     x11=1d0/(r*rxi)
  !     x32=(r+rxi)*x11*x11/r
  !  ! endif
  !   !----
  !  ! if(ket.eq.1) then
  !  !   ale=-dlog(r-et)
  !  !   y11=0d0
  !  !   y32=0d0
  !  ! else
  !     ret=r+et
  !     ale=dlog(ret)
  !     y11=1d0/(r*ret)
  !     y32=(r+ret)*y11*y11/r
  !  ! endif
  !   !----
  !   ey=sd/r-yp*q/r3
  !   ez=cd/r+d*q/r3
  !   fy=d/r3+xi2*y32*sd
  !   fz=yp/r3+xi2*y32*cd
  !   gy=2d0*x11*sd-yp*q*x32
  !   gz=2d0*x11*cd+d*q*x32
  !   hy=d*q*x32+xi*q*y32*sd
  !   hz=yp*q*x32+xi*q*y32*cd
  !   ! !write(*,*) 'dccon2'
  !   ! !write(*,*)xi2,et2,q2,r,r2,r3,r5,d,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
  !
  !   return
  ! end subroutine
end module mod_okada
