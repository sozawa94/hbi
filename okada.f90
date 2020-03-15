module mod_okada
  real(8),parameter::pi2=8.d0*datan(1.d0)
  real(8),parameter::epso=1d-6
  real(8),private::xi2,et2,q2,r,r2,r3,r5,d,yp,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz
  real(8),private::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
  !real(8)::p,q,s,t,xy,yp,x2,y2,d2,qr,qrx,a3,a5,b3,c3,fuy,vy,wy,fuz,vz,wz

contains
  ! subroutine  dc3d0(alpha,x,y,z,depth,dip,pot1,pot2,pot3,pot4,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
  !   implicit none
  !   real(8),intent(in)::alpha,x,y,z,depth,dip,pot1,pot2,pot3,pot4
  !   real(8),intent(out)::ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
  !   integer,intent(out)::iret
  !   ! c
  !   ! !*******************************************************************
  !   ! !****                                                          *****
  !   ! !****    displacement and strain at depth                      *****
  !   ! !****    due to buried point source in a semiinfinite medium   *****
  !   ! !****                         coded by  y.okada ... sep.1991   *****
  !   ! !****                         revised     nov.1991, may.2002   *****
  !   ! !****                                                          *****
  !   ! !*******************************************************************
  !   ! c
  !   ! !**** input
  !   ! !****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)
  !   ! !****   x,y,z : coordinate of observing point
  !   ! !****   depth : source depth
  !   ! !****   dip   : dip-angle (degree)
  !   ! !****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
  !   ! !****       potency=(  moment of double-couple  )/myu     for pot1,2
  !   ! !****       potency=(intensity of isotropic part)/lambda  for pot3
  !   ! !****       potency=(intensity of linear dipole )/myu     for pot4
  !   ! c
  !   ! !**** output
  !   ! !****   ux, uy, uz  : displacement ( unit=(unit of potency) /
  !   ! !****               :                     (unit of x,y,z,depth)**2  )
  !   ! !****   uxx,uyx,uzx : x-derivative ( unit= unit of potency) /
  !   ! !****   uxy,uyy,uzy : y-derivative        (unit of x,y,z,depth)**3  )
  !   ! !****   uxz,uyz,uzz : z-derivative
  !   ! !****   iret        : return code
  !   ! !****               :   =0....normal
  !   ! !****               :   =1....singular
  !   ! !****               :   =2....positive z was given
  !   ! c
  !   real(8)::dummy(8),r,aalpha,ddip,xx,yy,zz,pp1,pp2,pp3,pp4,dd,du
  !   real(8)::  u(12),dua(12),dub(12),duc(12)
  !   integer::i
  !
  !   iret=0
  !   if(z.gt.0.) then
  !     iret=2
  !     go to 99
  !   endif
  !
  !   do i=1,12
  !     u(i)=0d0
  !     dua(i)=0d0
  !     dub(i)=0d0
  !     duc(i)=0d0
  !   end do
  !   aalpha=alpha
  !   ddip=dip
  !   call dccon0(aalpha,ddip)
  !
  !   !=====  real-source contribution  =====
  !   xx=x
  !   yy=y
  !   zz=z
  !   dd=depth+z
  !   call dccon1(xx,yy,dd)
  !   if(r.eq.0d0) then
  !     iret=1
  !     go to 99
  !   endif
  !   !!----
  !   pp1=pot1
  !   pp2=pot2
  !   pp3=pot3
  !   pp4=pot4
  !   call ua0(xx,yy,dd,pp1,pp2,pp3,pp4,dua)
  !   !!----
  !   do i=1,12
  !     if(i.lt.10) u(i)=u(i)-dua(i)
  !     if(i.ge.10) u(i)=u(i)+dua(i)
  !   end do
  !   !!====  image-source contribution  =====
  !   dd=depth-z
  !   call dccon1(xx,yy,dd)
  !   call ua0(xx,yy,dd,pp1,pp2,pp3,pp4,dua)
  !   call ub0(xx,yy,dd,zz,pp1,pp2,pp3,pp4,dub)
  !   call uc0(xx,yy,dd,zz,pp1,pp2,pp3,pp4,duc)
  !   !----
  !   do i=1,12
  !     du=dua(i)+dub(i)+zz*duc(i)
  !     if(i.ge.10) du=du+duc(i-9)
  !     u(i)=u(i)+du
  !   end do
  !   !====
  !   ux=u(1)
  !   uy=u(2)
  !   uz=u(3)
  !   uxx=u(4)
  !   uyx=u(5)
  !   uzx=u(6)
  !   uxy=u(7)
  !   uyy=u(8)
  !   uzy=u(9)
  !   uxz=u(10)
  !   uyz=u(11)
  !   uzz=u(12)
  !   return
  !   !======================================
  !   !====  in case of singular (r=0)  =====
  !   !======================================
  !   99 ux=0d0
  !   uy=0d0
  !   uz=0d0
  !   uxx=0d0
  !   uyx=0d0
  !   uzx=0d0
  !   uxy=0d0
  !   uyy=0d0
  !   uzy=0d0
  !   uxz=0d0
  !   uyz=0d0
  !   uzz=0d0
  !   return
  ! end subroutine dc3d0
  ! subroutine  ua0(x,y,d,pot1,pot2,pot3,pot4,u)
  !   implicit none
  !   real(8),intent(in)::x,y,d,pot1,pot2,pot3,pot4
  !   real(8),intent(out)::u(12)
  !   real(8)::du(12)
  !   integer::i
  !   !real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
  !   !real(8)::p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,fuy,vy,wy,fuz,vz,wz
  !   !c
  !   !*******************************************************************
  !   !****    displacement and strain at depth (part-a)             *****
  !   !****    due to buried point source in a semiinfinite medium   *****
  !   !*******************************************************************
  !   !c
  !   !**** input
  !   !****   x,y,d : station coordinates in fault system
  !   !****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
  !   !**** output
  !   !****   u(12) : displacement and their derivatives
  !   !c
  !
  !   !----
  !   do i=1,12
  !     u(i)=0d0
  !   end do
  !   !=====================================
  !   !====  strike-slip contribution  =====
  !   !=====================================
  !   if(pot1.ne.0d0) then
  !     du( 1)= alp1*q/r3    +alp2*x2*qr
  !     du( 2)= alp1*x/r3*sd +alp2*xy*qr
  !     du( 3)=-alp1*x/r3*cd +alp2*x*d*qr
  !     du( 4)= x*qr*(-alp1 +alp2*(1d0+a5) )
  !     du( 5)= alp1*a3/r3*sd +alp2*y*qr*a5
  !     du( 6)=-alp1*a3/r3*cd +alp2*d*qr*a5
  !     du( 7)= alp1*(sd/r3-y*qr) +alp2*3d0*x2/r5*fuy
  !     du( 8)= 3d0*x/r5*(-alp1*y*sd +alp2*(y*fuy+q) )
  !     du( 9)= 3d0*x/r5*( alp1*y*cd +alp2*d*fuy )
  !     du(10)= alp1*(cd/r3+d*qr) +alp2*3d0*x2/r5*fuz
  !     du(11)= 3d0*x/r5*( alp1*d*sd +alp2*y*fuz )
  !     du(12)= 3d0*x/r5*(-alp1*d*cd +alp2*(d*fuz-q) )
  !     do  i=1,12
  !       u(i)=u(i)+pot1/pi2*du(i)
  !     end do
  !   endif
  !   !==================================
  !   !====  dip-slip contribution  =====
  !   !==================================
  !   if(pot2.ne.0d0) then
  !     du( 1)=            alp2*x*p*qr
  !     du( 2)= alp1*s/r3 +alp2*y*p*qr
  !     du( 3)=-alp1*t/r3 +alp2*d*p*qr
  !     du( 4)=                 alp2*p*qr*a5
  !     du( 5)=-alp1*3d0*x*s/r5 -alp2*y*p*qrx
  !     du( 6)= alp1*3d0*x*t/r5 -alp2*d*p*qrx
  !     du( 7)=                          alp2*3d0*x/r5*vy
  !     du( 8)= alp1*(s2d/r3-3d0*y*s/r5) +alp2*(3d0*y/r5*vy+p*qr)
  !     du( 9)=-alp1*(c2d/r3-3d0*y*t/r5) +alp2*3d0*d/r5*vy
  !     du(10)=                          alp2*3d0*x/r5*vz
  !     du(11)= alp1*(c2d/r3+3d0*d*s/r5) +alp2*3d0*y/r5*vz
  !     du(12)= alp1*(s2d/r3-3d0*d*t/r5) +alp2*(3d0*d/r5*vz-p*qr)
  !     do i=1,12
  !       u(i)=u(i)+pot2/pi2*du(i)
  !     end do
  !   endif
  !   !=======================================
  !   !====  tensile-fault contribution  =====
  !   !=======================================
  !   if(pot3.ne.0d0) then
  !     du( 1)= alp1*x/r3 -alp2*x*q*qr
  !     du( 2)= alp1*t/r3 -alp2*y*q*qr
  !     du( 3)= alp1*s/r3 -alp2*d*q*qr
  !     du( 4)= alp1*a3/r3     -alp2*q*qr*a5
  !     du( 5)=-alp1*3d0*x*t/r5 +alp2*y*q*qrx
  !     du( 6)=-alp1*3d0*x*s/r5 +alp2*d*q*qrx
  !     du( 7)=-alp1*3d0*xy/r5           -alp2*x*qr*wy
  !     du( 8)= alp1*(c2d/r3-3d0*y*t/r5) -alp2*(y*wy+q)*qr
  !     du( 9)= alp1*(s2d/r3-3d0*y*s/r5) -alp2*d*qr*wy
  !     du(10)= alp1*3d0*x*d/r5          -alp2*x*qr*wz
  !     du(11)=-alp1*(s2d/r3-3d0*d*t/r5) -alp2*y*qr*wz
  !     du(12)= alp1*(c2d/r3+3d0*d*s/r5) -alp2*(d*wz-q)*qr
  !     do i=1,12
  !       u(i)=u(i)+pot3/pi2*du(i)
  !     end do
  !   endif
  !   !========================================
  !   !====  inflate source contribution  =====
  !   !========================================
  !   if(pot4.ne.0d0) then
  !     du( 1)=-alp1*x/r3
  !     du( 2)=-alp1*y/r3
  !     du( 3)=-alp1*d/r3
  !     du( 4)=-alp1*a3/r3
  !     du( 5)= alp1*3d0*xy/r5
  !     du( 6)= alp1*3d0*x*d/r5
  !     du( 7)= du(5)
  !     du( 8)=-alp1*b3/r3
  !     du( 9)= alp1*3d0*y*d/r5
  !     du(10)=-du(6)
  !     du(11)=-du(9)
  !     du(12)= alp1*c3/r3
  !     do i=1,12
  !       u(i)=u(i)+pot4/pi2*du(i)
  !     end do
  !   endif
  !   return
  ! end subroutine
  ! subroutine  ub0(x,y,d,z,pot1,pot2,pot3,pot4,u)
  !   implicit none
  !   real(8),intent(in)::x,y,d,z,pot1,pot2,pot3,pot4
  !   real(8),intent(out)::u(12)
  !   real(8)::du(12),rd,d12,d32,d33,d53,d54,fi1,fi2,fi3,fi4,fi5,fj1,fj2,fj3,fj4
  !   real(8)::fk1,fk2,fk3,c
  !   integer::i
  !   !c
  !   !*******************************************************************
  !   !****    displacement and strain at depth (part-b)             *****
  !   !****    due to buried point source in a semiinfinite medium   *****
  !   !*******************************************************************
  !   !c
  !   !**** input
  !   !****   x,y,d,z : station coordinates in fault system
  !   !****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
  !   !**** output
  !   !****   u(12) : displacement and their derivatives
  !   !c
  !   !real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
  !   !real(8)::p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,fuy,vy,wy,fuz,vz,wz
  !   !----
  !   !d+z
  !   rd=r+d
  !   d12=1d0/(r*rd*rd)
  !   d32=d12*(2d0*r+d)/r2
  !   d33=d12*(3d0*r+d)/(r2*rd)
  !   d53=d12*(8d0*r2+9d0*r*d+3d0*d2)/(r2*r2*rd)
  !   d54=d12*(5d0*r2+4d0*r*d+d2)/r3*d12
  !   !----
  !   fi1= y*(d12-x2*d33)
  !   fi2= x*(d12-y2*d33)
  !   fi3= x/r3-fi2
  !   fi4=-xy*d32
  !   fi5= 1d0/(r*rd)-x2*d32
  !   fj1=-3d0*xy*(d33-x2*d54)
  !   fj2= 1d0/r3-3d0*d12+3d0*x2*y2*d54
  !   fj3= a3/r3-fj2
  !   fj4=-3d0*xy/r5-fj1
  !   fk1=-y*(d32-x2*d53)
  !   fk2=-x*(d32-y2*d53)
  !   fk3=-3d0*x*d/r5-fk2
  !   !----
  !   do i=1,12
  !     u(i)=0d0
  !   end do
  !   !=====================================
  !   !====  strike-slip contribution  =====
  !   !=====================================
  !   if(pot1.ne.0d0) then
  !     du( 1)=-x2*qr  -alp3*fi1*sd
  !     du( 2)=-xy*qr  -alp3*fi2*sd
  !     du( 3)=-c*x*qr -alp3*fi4*sd
  !     du( 4)=-x*qr*(1d0+a5) -alp3*fj1*sd
  !     du( 5)=-y*qr*a5      -alp3*fj2*sd
  !     du( 6)=-c*qr*a5      -alp3*fk1*sd
  !     du( 7)=-3d0*x2/r5*fuy      -alp3*fj2*sd
  !     du( 8)=-3d0*xy/r5*fuy-x*qr -alp3*fj4*sd
  !     du( 9)=-3d0*c*x/r5*fuy     -alp3*fk2*sd
  !     du(10)=-3d0*x2/r5*fuz  +alp3*fk1*sd
  !     du(11)=-3d0*xy/r5*fuz  +alp3*fk2*sd
  !     du(12)= 3d0*x/r5*(-c*fuz +alp3*y*sd)
  !     do i=1,12
  !       u(i)=u(i)+pot1/pi2*du(i)
  !     end do
  !   endif
  !   !==================================
  !   !====  dip-slip contribution  =====
  !   !==================================
  !   if(pot2.ne.0d0) then
  !     du( 1)=-x*p*qr +alp3*fi3*sdcd
  !     du( 2)=-y*p*qr +alp3*fi1*sdcd
  !     du( 3)=-c*p*qr +alp3*fi5*sdcd
  !     du( 4)=-p*qr*a5 +alp3*fj3*sdcd
  !     du( 5)= y*p*qrx +alp3*fj1*sdcd
  !     du( 6)= c*p*qrx +alp3*fk3*sdcd
  !     du( 7)=-3d0*x/r5*vy      +alp3*fj1*sdcd
  !     du( 8)=-3d0*y/r5*vy-p*qr +alp3*fj2*sdcd
  !     du( 9)=-3d0*c/r5*vy      +alp3*fk1*sdcd
  !     du(10)=-3d0*x/r5*vz -alp3*fk3*sdcd
  !     du(11)=-3d0*y/r5*vz -alp3*fk1*sdcd
  !     du(12)=-3d0*c/r5*vz +alp3*a3/r3*sdcd
  !     do i=1,12
  !       u(i)=u(i)+pot2/pi2*du(i)
  !     end do
  !   endif
  !   !=======================================
  !   !====  tensile-fault contribution  =====
  !   !=======================================
  !   if(pot3.ne.0d0) then
  !     du( 1)= x*q*qr -alp3*fi3*sdsd
  !     du( 2)= y*q*qr -alp3*fi1*sdsd
  !     du( 3)= c*q*qr -alp3*fi5*sdsd
  !     du( 4)= q*qr*a5 -alp3*fj3*sdsd
  !     du( 5)=-y*q*qrx -alp3*fj1*sdsd
  !     du( 6)=-c*q*qrx -alp3*fk3*sdsd
  !     du( 7)= x*qr*wy     -alp3*fj1*sdsd
  !     du( 8)= qr*(y*wy+q) -alp3*fj2*sdsd
  !     du( 9)= c*qr*wy     -alp3*fk1*sdsd
  !     du(10)= x*qr*wz +alp3*fk3*sdsd
  !     du(11)= y*qr*wz +alp3*fk1*sdsd
  !     du(12)= c*qr*wz -alp3*a3/r3*sdsd
  !     do i=1,12
  !       u(i)=u(i)+pot3/pi2*du(i)
  !     end do
  !   endif
  !   !========================================
  !   !====  inflate source contribution  =====
  !   !========================================
  !   if(pot4.ne.0d0) then
  !     du( 1)= alp3*x/r3
  !     du( 2)= alp3*y/r3
  !     du( 3)= alp3*d/r3
  !     du( 4)= alp3*a3/r3
  !     du( 5)=-alp3*3d0*xy/r5
  !     du( 6)=-alp3*3d0*x*d/r5
  !     du( 7)= du(5)
  !     du( 8)= alp3*b3/r3
  !     du( 9)=-alp3*3d0*y*d/r5
  !     du(10)=-du(6)
  !     du(11)=-du(9)
  !     du(12)=-alp3*c3/r3
  !     do i=1,12
  !       u(i)=u(i)+pot4/pi2*du(i)
  !     end do
  !   endif
  !   return
  ! end subroutine
  ! subroutine  uc0(x,y,d,z,pot1,pot2,pot3,pot4,u)
  !   implicit none
  !   real(8),intent(in)::x,y,d,z,pot1,pot2,pot3,pot4
  !   real(8),intent(out)::u(12)
  !   real(8)::du(12),q2,r7,a7,b5,b7,c5,c7,d7,qr5,qr7,dr5,c
  !   integer::i    !c
  !   !*******************************************************************
  !   !****    displacement and strain at depth (part-b)             *****
  !   !****    due to buried point source in a semiinfinite medium   *****
  !   !*******************************************************************
  !   !c
  !   !**** input
  !   !****   x,y,d,z : station coordinates in fault system
  !   !****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
  !   !**** output
  !   !****   u(12) : displacement and their derivatives
  !   !c
  !   !real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
  !   !real(8)::p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3
  !   !----
  !   !d+z
  !   q2=q*q
  !   r7=r5*r2
  !   a7=1d0-7d0*x2/r2
  !   b5=1d0-5d0*y2/r2
  !   b7=1d0-7d0*y2/r2
  !   c5=1d0-5d0*d2/r2
  !   c7=1d0-7d0*d2/r2
  !   d7=2d0-7d0*q2/r2
  !   qr5=5d0*q/r2
  !   qr7=7d0*q/r2
  !   dr5=5d0*d/r2
  !   !----
  !   do i=1,12
  !     u(i)=0d0
  !   end do
  !   !=====================================
  !   !====  strike-slip contribution  =====
  !   !=====================================
  !   if(pot1.ne.0d0) then
  !     du( 1)=-alp4*a3/r3*cd  +alp5*c*qr*a5
  !     du( 2)= 3d0*x/r5*( alp4*y*cd +alp5*c*(sd-y*qr5) )
  !     du( 3)= 3d0*x/r5*(-alp4*y*sd +alp5*c*(cd+d*qr5) )
  !     du( 4)= alp4*3d0*x/r5*(2d0+a5)*cd   -alp5*c*qrx*(2d0+a7)
  !     du( 5)= 3d0/r5*( alp4*y*a5*cd +alp5*c*(a5*sd-y*qr5*a7) )
  !     du( 6)= 3d0/r5*(-alp4*y*a5*sd +alp5*c*(a5*cd+d*qr5*a7) )
  !     du( 7)= du(5)
  !     du( 8)= 3d0*x/r5*( alp4*b5*cd -alp5*5d0*c/r2*(2d0*y*sd+q*b7) )
  !     du( 9)= 3d0*x/r5*(-alp4*b5*sd +alp5*5d0*c/r2*(d*b7*sd-y*c7*cd) )
  !     du(10)= 3d0/r5*   (-alp4*d*a5*cd +alp5*c*(a5*cd+d*qr5*a7) )
  !     du(11)= 15*x/r7*( alp4*y*d*cd  +alp5*c*(d*b7*sd-y*c7*cd) )
  !     du(12)= 15*x/r7*(-alp4*y*d*sd  +alp5*c*(2d0*d*cd-q*c7) )
  !     do i=1,12
  !       u(i)=u(i)+pot1/pi2*du(i)
  !     end do
  !   endif
  !   !==================================
  !   !====  dip-slip contribution  =====
  !   !==================================
  !   if(pot2.ne.0d0) then
  !     du( 1)= alp4*3d0*x*t/r5          -alp5*c*p*qrx
  !     du( 2)=-alp4/r3*(c2d-3d0*y*t/r2) +alp5*3d0*c/r5*(s-y*p*qr5)
  !     du( 3)=-alp4*a3/r3*sdcd         +alp5*3d0*c/r5*(t+d*p*qr5)
  !     du( 4)= alp4*3d0*t/r5*a5              -alp5*5d0*c*p*qr/r2*a7
  !     du( 5)= 3d0*x/r5*(alp4*(c2d-5d0*y*t/r2)-alp5*5d0*c/r2*(s-y*p*qr7))
  !     du( 6)= 3d0*x/r5*(alp4*(2d0+a5)*sdcd   -alp5*5d0*c/r2*(t+d*p*qr7))
  !     du( 7)= du(5)
  !     du( 8)= 3d0/r5*(alp4*(2d0*y*c2d+t*b5)+alp5*c*(s2d-10*y*s/r2-p*qr5*b7))
  !     du( 9)= 3d0/r5*(alp4*y*a5*sdcd-alp5*c*((3d0+a5)*c2d+y*p*dr5*qr7))
  !     du(10)= 3d0*x/r5*(-alp4*(s2d-t*dr5) -alp5*5d0*c/r2*(t+d*p*qr7))
  !     du(11)= 3d0/r5*(-alp4*(d*b5*c2d+y*c5*s2d)-alp5*c*((3d0+a5)*c2d+y*p*dr5*qr7))
  !     du(12)= 3d0/r5*(-alp4*d*a5*sdcd-alp5*c*(s2d-10*d*t/r2+p*qr5*c7))
  !     do i=1,12
  !       u(i)=u(i)+pot2/pi2*du(i)
  !     end do
  !   endif
  !   !=======================================
  !   !====  tensile-fault contribution  =====
  !   !=======================================
  !   if(pot3.ne.0d0) then
  !     du( 1)= 3d0*x/r5*(-alp4*s +alp5*(c*q*qr5-z))
  !     du( 2)= alp4/r3*(s2d-3d0*y*s/r2)+alp5*3d0/r5*(c*(t-y+y*q*qr5)-y*z)
  !     du( 3)=-alp4/r3*(1d0-a3*sdsd)   -alp5*3d0/r5*(c*(s-d+d*q*qr5)-d*z)
  !     du( 4)=-alp4*3d0*s/r5*a5 +alp5*(c*qr*qr5*a7-3d0*z/r5*a5)
  !     du( 5)= 3d0*x/r5*(-alp4*(s2d-5d0*y*s/r2)-alp5*5d0/r2*(c*(t-y+y*q*qr7)-y*z))
  !     du( 6)= 3d0*x/r5*( alp4*(1d0-(2d0+a5)*sdsd)+alp5*5d0/r2*(c*(s-d+d*q*qr7)-d*z))
  !     du( 7)= du(5)
  !     du( 8)= 3d0/r5*(-alp4*(2d0*y*s2d+s*b5)-alp5*(c*(2d0*sdsd+10*y*(t-y)/r2-q*qr5*b7)+z*b5))
  !     du( 9)= 3d0/r5*( alp4*y*(1d0-a5*sdsd)+alp5*(c*(3d0+a5)*s2d-y*dr5*(c*d7+z)))
  !     du(10)= 3d0*x/r5*(-alp4*(c2d+s*dr5)+alp5*(5d0*c/r2*(s-d+d*q*qr7)-1d0-z*dr5))
  !     du(11)= 3d0/r5*( alp4*(d*b5*s2d-y*c5*c2d)+alp5*(c*((3d0+a5)*s2d-y*dr5*d7)-y*(1d0+z*dr5)))
  !     du(12)= 3d0/r5*(-alp4*d*(1d0-a5*sdsd)-alp5*(c*(c2d+10*d*(s-d)/r2-q*qr5*c7)+z*(1d0+c5)))
  !     do i=1,12
  !       u(i)=u(i)+pot3/pi2*du(i)
  !     end do
  !   endif
  !   !========================================
  !   !====  inflate source contribution  =====
  !   !========================================
  !   if(pot4.ne.0d0) then
  !     du( 1)= alp4*3d0*x*d/r5
  !     du( 2)= alp4*3d0*y*d/r5
  !     du( 3)= alp4*c3/r3
  !     du( 4)= alp4*3d0*d/r5*a5
  !     du( 5)=-alp4*15*xy*d/r7
  !     du( 6)=-alp4*3d0*x/r5*c5
  !     du( 7)= du(5)
  !     du( 8)= alp4*3d0*d/r5*b5
  !     du( 9)=-alp4*3d0*y/r5*c5
  !     du(10)= du(6)
  !     du(11)= du(9)
  !     du(12)= alp4*3d0*d/r5*(2d0+c5)
  !     do i=1,12
  !       u(i)=u(i)+pot4/pi2*du(i)
  !     end do
  !   endif
  !   return
  ! end subroutine
  subroutine  dc3d(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,disl1,disl2,disl3,&
    &ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    implicit none
    real(8),intent(in)::alpha,x,y,z,depth,dip,al1,al2,aw1,aw2,disl1,disl2,disl3
    real(8),intent(out)::ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    integer,intent(out)::iret
    integer::i,j,k
    real(8)::ddip,aalpha,zz,dd1,dd2,dd3,r12,r21,r22,p,q
    !c
    !*******************************************************************
    !****                                                          *****
    !****    displacement and strain at depth                      *****
    !****    due to buried finite fault in a semiinfinite medium   *****
    !****              coded by  y.okada ... sep.1991              *****
    !****              revised ... nov.1991, apr.1992, may.1993,   *****
    !****                          jul.1993, may.2002              *****
    !*******************************************************************
    !c
    !**** input
    !****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)
    !****   x,y,z : coordinate of observing point
    !****   depth : depth of reference point
    !****   dip   : dip-angle (degree)
    !****   al1,al2   : fault length range
    !****   aw1,aw2   : fault width range
    !****   disl1-disl3 : strike-, dip-, tensile-dislocations
    !c
    !**** output
    !****   ux, uy, uz  : displacement ( unit=(unit of disl)
    !****   uxx,uyx,uzx : x-derivative ( unit=(unit of disl) /
    !****   uxy,uyy,uzy : y-derivative        (unit of x,y,z,depth,al,aw) )
    !****   uxz,uyz,uzz : z-derivative
    !****   iret        : return code
    !****               :   =0....normal
    !****               :   =1....singular
    !****               :   =2....positive z was given
    !c
    real(8):: dummy(5)
    real(8)::  xi(2),et(2)
    integer::kxi(2),ket(2)
    real(8)::  u(12),du(12),dua(12),dub(12),duc(12)
    !real(8)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
    !----
    iret=0
    if(z.gt.0.) then
      iret=2
      go to 99
    endif
    !----
    do i=1,12
      u(i)=0d0
      dua(i)=0d0
      dub(i)=0d0
      duc(i)=0d0
    end do
    aalpha=alpha
    ddip=dip
    call dccon0(aalpha,ddip)
    !----
    !write(*,*)'sd,cd',sd,cd

    zz=z
    dd1=disl1
    dd2=disl2
    dd3=disl3
    xi(1)=x-al1
    xi(2)=x-al2
    !write(*,*) 'xi',xi
    if(dabs(xi(1)).lt.epso) xi(1)=0d0
    if(dabs(xi(2)).lt.epso) xi(2)=0d0
    !=====================================
    !====  real-source contribution  =====
    !=====================================
    d=depth+z
    !write(*,*)sd,cd
    p=y*cd+d*sd
    q=y*sd-d*cd
    et(1)=p-aw1
    et(2)=p-aw2
    !write(*,*) 'q,et',q,et
    if(dabs(q).lt.epso)  q=0d0
    if(dabs(et(1)).lt.epso) et(1)=0d0
    if(dabs(et(2)).lt.epso) et(2)=0d0
    !-------------------------------
    !---- reject singular case -----
    !-------------------------------
    !---- on fault edge
    if(q.eq.0d0 .and.((xi(1)*xi(2).le.0d0 .and. et(1)*et(2).eq.0d0).or.(et(1)*et(2).le.0d0 .and. xi(1)*xi(2).eq.0d0)))  then
      iret=1
      go to 99
    endif
    !---- on negative extension of fault edge
    kxi(1)=0
    kxi(2)=0
    ket(1)=0
    ket(2)=0
    r12=dsqrt(xi(1)*xi(1)+et(2)*et(2)+q*q)
    r21=dsqrt(xi(2)*xi(2)+et(1)*et(1)+q*q)
    r22=dsqrt(xi(2)*xi(2)+et(2)*et(2)+q*q)
    !write(*,*) 'r12,r21,r22',r12,r21,r22
    if(xi(1).lt.0d0 .and. r21+xi(2).lt.epso) kxi(1)=1
    if(xi(1).lt.0d0 .and. r22+xi(2).lt.epso) kxi(2)=1
    if(et(1).lt.0d0 .and. r12+et(2).lt.epso) ket(1)=1
    if(et(1).lt.0d0 .and. r22+et(2).lt.epso) ket(2)=1
    !====
    do k=1,2
      do j=1,2
        call dccon2(xi(j),et(k),q,sd,cd,kxi(k),ket(j))
        call ua(xi(j),et(k),q,dd1,dd2,dd3,dua)
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
    !======================================
    !====  image-source contribution  =====
    !======================================
    !write(*,*)'image-source contribution'
    d=depth-z
    p=y*cd+d*sd
    q=y*sd-d*cd
    et(1)=p-aw1
    et(2)=p-aw2
    if(dabs(q).lt.epso)  q=0d0
    if(dabs(et(1)).lt.epso) et(1)=0d0
    if(dabs(et(2)).lt.epso) et(2)=0d0
    !-------------------------------
    !---- reject singular case -----
    !-------------------------------
    !---- on fault edge
    if(q.eq.0d0 .and.((xi(1)*xi(2).le.0d0 .and. et(1)*et(2).eq.0d0).or.(et(1)*et(2).le.0d0 .and. xi(1)*xi(2).eq.0d0)))  then
      iret=1
      go to 99
    endif
    !---- on negative extension of fault edge
    kxi(1)=0
    kxi(2)=0
    ket(1)=0
    ket(2)=0
    r12=dsqrt(xi(1)*xi(1)+et(2)*et(2)+q*q)
    r21=dsqrt(xi(2)*xi(2)+et(1)*et(1)+q*q)
    r22=dsqrt(xi(2)*xi(2)+et(2)*et(2)+q*q)
    if(xi(1).lt.0d0 .and. r21+xi(2).lt.epso) kxi(1)=1
    if(xi(1).lt.0d0 .and. r22+xi(2).lt.epso) kxi(2)=1
    if(et(1).lt.0d0 .and. r12+et(2).lt.epso) ket(1)=1
    if(et(1).lt.0d0 .and. r22+et(2).lt.epso) ket(2)=1
    !====
    do k=1,2
      do j=1,2
        call dccon2(xi(j),et(k),q,sd,cd,kxi(k),ket(j))
        call ua(xi(j),et(k),q,dd1,dd2,dd3,dua)
        call ub(xi(j),et(k),q,dd1,dd2,dd3,dub)
        call uc(xi(j),et(k),q,zz,dd1,dd2,dd3,duc)
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
    !==========================================
    !====  in case of singular (on edge)  =====
    !==========================================
    99 ux=0d0
    uy=0d0
    uz=0d0
    uxx=0d0
    uyx=0d0
    uzx=0d0
    uxy=0d0
    uyy=0d0
    uzy=0d0
    uxz=0d0
    uyz=0d0
    uzz=0d0
    return
  end subroutine dc3d
  subroutine  ua(xi,et,q,disl1,disl2,disl3,u)
    implicit none
    real(8),intent(in)::xi,et,q,disl1,disl2,disl3
    real(8),intent(out)::u(12)
    real(8)::du(12),qx,qy,xy
    integer::i
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
    do i=1,12
      u(i)=0d0
    end do
    xy=xi*y11
    qx=q *x11
    qy=q *y11
    ! !write(*,*)'xy,qx,qy'
    ! !write(*,*)xy,qx,qy
    !=====================================
    !====  strike-slip contribution  =====
    !=====================================
    if(disl1.ne.0d0) then
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
        u(i)=u(i)+disl1/pi2*du(i)
      end do
      !!write(*,*)'u',u
    endif
    !=====================================
    !====    dip-slip contribution   =====
    !=====================================
    if(disl2.ne.0d0) then
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
        u(i)=u(i)+disl2/pi2*du(i)
      end do
    endif
    !=======================================
    !====  tensile-fault contribution  =====
    !=======================================
    if(disl3.ne.0d0) then
      du( 1)=-alp1*ale -alp2*q*qy
      du( 2)=-alp1*alx -alp2*q*qx
      du( 3)=    tt/2d0 -alp2*(et*qx+xi*qy)
      du( 4)=-alp1*xy  +alp2*xi*q2*y32
      du( 5)=-alp1/r   +alp2*q2/r3
      du( 6)=-alp1*qy  -alp2*q*q2*y32
      du( 7)=-alp1*(cd/r+qy*sd)  -alp2*q*fy
      du( 8)=-alp1*yp*x11         -alp2*q*gy
      du( 9)= alp1*(d*x11+xy*sd) +alp2*q*hy
      du(10)= alp1*(sd/r-qy*cd)  -alp2*q*fz
      du(11)= alp1*d*x11         -alp2*q*gz
      du(12)= alp1*(yp*x11+xy*cd) +alp2*q*hz
      !write(*,*)'du',du
      do i=1,12
        u(i)=u(i)+disl3/pi2*du(i)
      end do
    endif
    !write(*,*)'u',u
    return
  end subroutine
  subroutine  ub(xi,et,q,disl1,disl2,disl3,u)
    implicit none
    real(8),intent(in)::xi,et,q,disl1,disl2,disl3
    real(8),intent(out)::u(12)
    real(8)::du(12),d11,rd,ai4,aj2,aj5,ai3,ak1,aj3,aj6,rd2,ai1,ai2,ak2,ak4,aj1,aj4,qx,qy
    real(8)::xs,ak3,xy
    integer::i

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
        ai4=1d0/cdcd*( xi/rd*sdcd+2d0*datan((et*(xs+q*cd)+xs*(r+xs)*sd)/(xi*(r+xs)*cd)) )
      endif
      ai3=(yp*cd/rd-ale+sd*dlog(rd))/cdcd
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
      u(i)=0d0
    end do
    qx=q*x11
    qy=q*y11
    !=====================================
    !====  strike-slip contribution  =====
    !=====================================
    if(disl1.ne.0d0) then
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
        u(i)=u(i)+disl1/pi2*du(i)
      end do
    endif
    !=====================================
    !====    dip-slip contribution   =====
    !=====================================
    if(disl2.ne.0d0) then
      du( 1)=-q/r      +alp3*ai3*sdcd
      du( 2)=-et*qx-tt -alp3*xi/rd*sdcd
      du( 3)= q*qx     +alp3*ai4*sdcd
      du( 4)= xi*q/r3     +alp3*aj4*sdcd
      du( 5)= et*q/r3+qy  +alp3*aj5*sdcd
      du( 6)=-q2/r3       +alp3*aj6*sdcd
      du( 7)=-ey          +alp3*aj1*sdcd
      du( 8)=-et*gy-xy*sd +alp3*aj2*sdcd
      du( 9)= q*gy        +alp3*aj3*sdcd
      du(10)=-ez          -alp3*ak3*sdcd
      du(11)=-et*gz-xy*cd -alp3*xi*d11*sdcd
      du(12)= q*gz        -alp3*ak4*sdcd
      !write(*,*)'ub,du',du
      do i=1,12
        u(i)=u(i)+disl2/pi2*du(i)
      end do
    endif
    !=======================================
    !====  tensile-fault contribution  =====
    !=======================================
    if(disl3.ne.0d0) then
      du( 1)= q*qy           -alp3*ai3*sdsd
      du( 2)= q*qx           +alp3*xi/rd*sdsd
      du( 3)= et*qx+xi*qy-tt -alp3*ai4*sdsd
      du( 4)=-xi*q2*y32 -alp3*aj4*sdsd
      du( 5)=-q2/r3     -alp3*aj5*sdsd
      du( 6)= q*q2*y32  -alp3*aj6*sdsd
      du( 7)= q*fy -alp3*aj1*sdsd
      du( 8)= q*gy -alp3*aj2*sdsd
      du( 9)=-q*hy -alp3*aj3*sdsd
      du(10)= q*fz +alp3*ak3*sdsd
      du(11)= q*gz +alp3*xi*d11*sdsd
      du(12)=-q*hz +alp3*ak4*sdsd
      !write(*,*)'ub,du',du
      do i=1,12
        u(i)=u(i)+disl3/pi2*du(i)
      end do
    endif
    !write(*,*)'ub,u',u
    return
  end
  subroutine  uc(xi,et,q,z,disl1,disl2,disl3,u)
    implicit none
    real(8),intent(in)::xi,et,q,z,disl1,disl2,disl3
    real(8),intent(out)::u(12)
    real(8)::du(12),c,x53,y53,h,z32,z53,y0,z0,ppy,ppz,qq,qqy,qqz,xy,qx,qy,qr,cqx,cdr,yy0
    integer::i    !c
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
      u(i)=0d0
    end do
    !=====================================
    !====  strike-slip contribution  =====
    !=====================================
    if(disl1.ne.0d0) then
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
      du(12)=         yy0*cd    -alp5*(cdr*sd-c*yp*qr-y0*sdsd+q*z0*cd)
      do i=1,12
        u(i)=u(i)+disl1/pi2*du(i)
      end do
    end if
    !=====================================
    !====    dip-slip contribution   =====
    !=====================================
    if(disl2.ne.0d0) then
      du( 1)= alp4*cd/r -qy*sd -alp5*c*q/r3
      du( 2)= alp4*yp*x11       -alp5*c*et*q*x32
      du( 3)=     -d*x11-xy*sd -alp5*c*(x11-q2*x32)
      du( 4)=-alp4*xi/r3*cd +alp5*c*xi*qr +xi*q*y32*sd
      du( 5)=-alp4*yp/r3     +alp5*c*et*qr
      du( 6)=    d/r3-y0*sd +alp5*c/r3*(1d0-3d0*q2/r2)
      du( 7)=-alp4*et/r3+y0*sdsd -alp5*(cdr*sd-c*yp*qr)
      du( 8)= alp4*(x11-yp*yp*x32) -alp5*c*((d+2d0*q*cd)*x32-yp*et*q*x53)
      du( 9)=  xi*ppy*sd+yp*d*x32 +alp5*c*((yp+2d0*q*sd)*x32-yp*q2*x53)
      du(10)=      -q/r3+y0*sdcd -alp5*(cdr*cd+c*d*qr)
      du(11)= alp4*yp*d*x32       -alp5*c*((yp-2d0*q*sd)*x32+d*et*q*x53)
      du(12)=-xi*ppz*sd+x11-d*d*x32-alp5*c*((d-2d0*q*cd)*x32-d*q2*x53)
      do i=1,12
        u(i)=u(i)+disl2/pi2*du(i)
      end do
    endif
    !=======================================
    !====  tensile-fault contribution  =====
    !=======================================
    if(disl3.ne.0d0) then
      du( 1)=-alp4*(sd/r+qy*cd)   -alp5*(z*y11-q2*z32)
      du( 2)= alp4*2d0*xy*sd+d*x11 -alp5*c*(x11-q2*x32)
      du( 3)= alp4*(yp*x11+xy*cd)  +alp5*q*(c*et*x32+xi*z32)
      du( 4)= alp4*xi/r3*sd+xi*q*y32*cd+alp5*xi*(3d0*c*et/r5-2d0*z32-z0)
      du( 5)= alp4*2d0*y0*sd-d/r3 +alp5*c/r3*(1d0-3d0*q2/r2)
      du( 6)=-alp4*yy0           -alp5*(c*et*qr-q*z0)
      du( 7)= alp4*(q/r3+y0*sdcd)   +alp5*(z/r3*cd+c*d*qr-q*z0*sd)
      du( 8)=-alp4*2d0*xi*ppy*sd-yp*d*x32 +alp5*c*((yp+2d0*q*sd)*x32-yp*q2*x53)
      du( 9)=-alp4*(xi*ppy*cd-x11+yp*yp*x32)    +alp5*(c*((d+2d0*q*cd)*x32-yp*et*q*x53)+xi*qqy)
      du(10)=  -et/r3+y0*cdcd -alp5*(z/r3*sd-c*yp*qr-y0*sdsd+q*z0*cd)
      du(11)= alp4*2d0*xi*ppz*sd-x11+d*d*x32 -alp5*c*((d-2d0*q*cd)*x32-d*q2*x53)
      du(12)= alp4*(xi*ppz*cd+yp*d*x32)     +alp5*(c*((yp-2d0*q*sd)*x32+d*et*q*x53)+xi*qqz)
      do i=1,12
        u(i)=u(i)+disl3/pi2*du(i)
      end do
    endif
    !write(*,*)u
    return
  end subroutine
  subroutine  dccon0(alpha,dip)
    implicit none
    !c
    !******************************************************************
    !****   calculate medium constants and fault-dip constants    *****
    !******************************************************************
    !c
    !**** input
    !****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)
    !****   dip   : dip-angle (degree)
    !c### caution ### if cos(dip) is sufficiently small, it is set to zero
    !c
    real(8),intent(in)::dip,alpha
    real(8)::p18
    !real(8),intent(out)::alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
    !----
    alp1=(1d0-alpha)/2d0
    alp2= alpha/2d0
    alp3=(1d0-alpha)/alpha
    alp4= 1d0-alpha
    alp5= alpha
    !----
    p18=pi2/360.d0
    sd=dsin(dip*p18)
    cd=dcos(dip*p18)
    if(dabs(cd).lt.epso) then
      cd=0d0
      if(sd.gt.0d0) sd= 1d0
      if(sd.lt.0d0) sd=-1d0
    endif
    sdsd=sd*sd
    cdcd=cd*cd
    sdcd=sd*cd
    s2d=2d0*sdcd
    c2d=cdcd-sdsd
    !write(*,*)'dccon0'
    !write(*,*)alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d
    return
  end subroutine
  ! subroutine  dccon1(x,y,d)
  !   implicit none
  !   real(8),intent(inout)::x,y,d
  !   real(8)::r7
  !   !c
  !   !*********************************************************************
  !   !****   calculate station geometry constants for point source    *****
  !   !*********************************************************************
  !   !c
  !   !**** input
  !   !****   x,y,d : station coordinates in fault system
  !   !c### caution ### if x,y,d are sufficiently small, they are set to zero
  !   !c
  !   !real(8)::dummy(5),sd,cd
  !   !real(8),intent(out)::p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,fuy,vy,wy,fuz,vz,wz
  !   !----
  !   if(dabs(x).lt.epso) x=0d0
  !   if(dabs(y).lt.epso) y=0d0
  !   if(dabs(d).lt.epso) d=0d0
  !   p=y*cd+d*sd
  !   q=y*sd-d*cd
  !   s=p*sd+q*cd
  !   t=p*cd-q*sd
  !   xy=x*y
  !   x2=x*x
  !   y2=y*y
  !   d2=d*d
  !   r2=x2+y2+d2
  !   r =dsqrt(r2)
  !   if(r.eq.0d0) return
  !   r3=r *r2
  !   r5=r3*r2
  !   r7=r5*r2
  !   !----
  !   a3=1d0-3d0*x2/r2
  !   a5=1d0-5d0*x2/r2
  !   b3=1d0-3d0*y2/r2
  !   c3=1d0-3d0*d2/r2
  !   !----
  !   qr=3d0*q/r5
  !   qrx=5d0*qr*x/r2
  !   !----
  !   fuy=sd-5d0*y*q/r2
  !   fuz=cd+5d0*d*q/r2
  !   vy=s-5d0*y*p*q/r2
  !   vz=t+5d0*d*p*q/r2
  !   wy=fuy+sd
  !   wz=fuz+cd
  !   !write(*,*) 'dccon1'
  !   !write(*,*)p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,fuy,vy,wy,fuz,vz,wz
  !   return
  ! end subroutine
  subroutine  dccon2(xi,et,q,sd,cd,kxi,ket)
    implicit none
    real(8),intent(inout)::xi,et,q,sd,cd
    integer,intent(in)::kxi,ket
    real(8)::rxi,ret
    !c
    !*********************************************************************
    !****   calculate station geometry constants for finite source   *****
    !*********************************************************************
    !c
    !**** input
    !****   xi,et,q : station coordinates in fault system
    !****   sd,cd   : sin, cos of dip-angle
    !****   kxi,ket : kxi=1, ket=1 means r+xi<epso, r+et<epso, respectively
    !c
    !c### caution ### if xi,et,q are sufficiently small, they are set to zer0
    !c
    !real(8)::sd,cd
    !----
    if(dabs(xi).lt.epso) xi=0d0
    if(dabs(et).lt.epso) et=0d0
    if(dabs( q).lt.epso)  q=0d0
    xi2=xi*xi
    et2=et*et
    q2=q*q
    r2=xi2+et2+q2
    r =dsqrt(r2)
    if(r.eq.0d0) return
    r3=r *r2
    r5=r3*r2
    yp =et*cd+q*sd
    d =et*sd-q*cd
    !----
    if(q.eq.0d0) then
      tt=0d0
    else
      tt=datan(xi*et/(q*r))
    endif
    !----
    if(kxi.eq.1) then
      alx=-dlog(r-xi)
      x11=0d0
      x32=0d0
    else
      rxi=r+xi
      alx=dlog(rxi)
      x11=1d0/(r*rxi)
      x32=(r+rxi)*x11*x11/r
    endif
    !----
    if(ket.eq.1) then
      ale=-dlog(r-et)
      y11=0d0
      y32=0d0
    else
      ret=r+et
      ale=dlog(ret)
      y11=1d0/(r*ret)
      y32=(r+ret)*y11*y11/r
    endif
    !----
    ey=sd/r-yp*q/r3
    ez=cd/r+d*q/r3
    fy=d/r3+xi2*y32*sd
    fz=yp/r3+xi2*y32*cd
    gy=2d0*x11*sd-yp*q*x32
    gz=2d0*x11*cd+d*q*x32
    hy=d*q*x32+xi*q*y32*sd
    hz=yp*q*x32+xi*q*y32*cd
    ! !write(*,*) 'dccon2'
    ! !write(*,*)xi2,et2,q2,r,r2,r3,r5,d,tt,alx,ale,x11,y11,x32,y32,ey,ez,fy,fz,gy,gz,hy,hz

    return
  end subroutine
end module
