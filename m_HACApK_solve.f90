!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.0.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
!C***********************************************************************
!C  This file includes routines for utilizing H-matrices, such as solving
!C  linear system with an H-matrix as the coefficient matrix and 
!C  multiplying an H-matrix and a vector,
!C  created by Akihiro Ida at Kyoto University on May 2012,
!C  last modified by Akihiro Ida on Sep 2014,
!C***********************************************************************
module m_HACApK_solve
 use m_HACApK_base
 implicit real*8(a-h,o-z)
 implicit integer*4(i-n)
contains

!***HACApK_adot_lfmtx_p
 subroutine HACApK_adot_lfmtx_p(zau,st_leafmtxp,st_ctl,zu,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(nd),zu(nd)
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4 :: ISTATUS(MPI_STATUS_SIZE),isct(2),irct(2)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
 allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 zau(:)=0.0d0
 call HACApK_adot_body_lfmtx(zau,st_leafmtxp,st_ctl,zu,nd)
 if(nrank==1) return
 wws(1:lnp(mpinr))=zau(lsp(mpinr):lsp(mpinr)+lnp(mpinr)-1)
 ncdp=mod(mpinr+1,nrank)
 ncsp=mod(mpinr+nrank-1,nrank)
! write(mpilog,1000) 'destination process=',ncdp,'; source process=',ncsp
 isct(1)=lnp(mpinr);isct(2)=lsp(mpinr); 
! irct=lnp(ncsp)
 do ic=1,nrank-1
!   idp=mod(mpinr+ic,nrank) ! rank of destination process
!   isp=mod(mpinr+nrank+ic-2,nrank) ! rank of source process
   call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
                     irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
!   write(mpilog,1000) 'ISTATUS=',ISTATUS,'; ierr=',ierr
!   write(mpilog,1000) 'ic=',ic,'; isct=',isct(1),'; irct=',irct(1),'; ivsps=',isct(2),'; ivspr=',irct(2)

   call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
                     wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
!   write(mpilog,1000) 'ISTATUS=',ISTATUS,'; ierr=',ierr
   
   zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
   wws(:irct(1))=wwr(:irct(1))
   isct=irct
!   write(mpilog,1000) 'ic=',ic,'; isct=',isct
 enddo
 deallocate(wws,wwr)
 end subroutine HACApK_adot_lfmtx_p
 
!***HACApK_adot_lfmtx_hyp
 subroutine HACApK_adot_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*),wws(*),wwr(*)
 integer*4 :: isct(*),irct(*)
 integer*4 :: ISTATUS(MPI_STATUS_SIZE)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
! integer*4,dimension(:),allocatable :: ISTATUS
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
! allocate(ISTATUS(MPI_STATUS_SIZE))
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 if(nrank>1)then
   wws(1:lnp(mpinr))=zau(lsp(mpinr):lsp(mpinr)+lnp(mpinr)-1)
   ncdp=mod(mpinr+1,nrank)
   ncsp=mod(mpinr+nrank-1,nrank)
   isct(1)=lnp(mpinr);isct(2)=lsp(mpinr); 
   do ic=1,nrank-1
     call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
                       irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
     call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
                       wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
     zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
     wws(:irct(1))=wwr(:irct(1))
     isct(:2)=irct(:2)
   enddo
 endif
!$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp

!***HACApK_adot_body_lfmtx
 RECURSIVE subroutine HACApK_adot_body_lfmtx(zau,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(nd),zu(nd)
 real*8,dimension(:),allocatable :: zbu
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 nlf=st_leafmtxp%nlf
 do ip=1,nlf
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     allocate(zbu(kt)); zbu(:)=0.0d0
     do il=1,kt
       do it=1,ndt; itt=it+nstrtt-1
         zbu(il)=zbu(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
     do il=1,kt
       do it=1,ndl; ill=it+nstrtl-1
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbu(il)
       enddo
     enddo
     deallocate(zbu)
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         zau(ill)=zau(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
   endif
 enddo
 end subroutine HACApK_adot_body_lfmtx

!***HACApK_adot_body_lfmtx_hyp
 subroutine HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*)
 real*8,dimension(:),allocatable :: zbut
 real*8,dimension(:),allocatable :: zaut
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),ltmp(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;ltmp(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 nlf=st_leafmtxp%nlf; ktmax=st_leafmtxp%ktmax
 ith = omp_get_thread_num()
 ith1 = ith+1
 nths=ltmp(ith); nthe=ltmp(ith1)-1
 allocate(zaut(nd)); zaut(:)=0.0d0
 allocate(zbut(ktmax)) 
 ls=nd; le=1
 do ip=nths,nthe
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(nstrtl<ls) ls=nstrtl; if(nstrtl+ndl-1>le) le=nstrtl+ndl-1
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     zbut(1:kt)=0.0d0
     do il=1,kt
       do it=1,ndt; itt=it+nstrtt-1
         zbut(il)=zbut(il)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
     do il=1,kt
       do it=1,ndl; ill=it+nstrtl-1
         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a2(it,il)*zbut(il)
       enddo
     enddo
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         zaut(ill)=zaut(ill)+st_leafmtxp%st_lf(ip)%a1(it,il)*zu(itt)
       enddo
     enddo
   endif
 enddo
 deallocate(zbut)
 
 do il=ls,le
!$omp atomic
   zau(il)=zau(il)+zaut(il)
 enddo
 end subroutine HACApK_adot_body_lfmtx_hyp
 
!***HACApK_adotsub_lfmtx_p
 subroutine HACApK_adotsub_lfmtx_p(zr,st_leafmtxp,st_ctl,zu,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zu(nd),zr(nd)
 real*8,dimension(:),allocatable :: zau
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;
 allocate(zau(nd))
 call HACApK_adot_lfmtx_p(zau,st_leafmtxp,st_ctl,zu,nd)
 zr(1:nd)=zr(1:nd)-zau(1:nd)
 deallocate(zau)
 end subroutine HACApK_adotsub_lfmtx_p
 
!***HACApK_adotsub_lfmtx_hyp
 subroutine HACApK_adotsub_lfmtx_hyp(zr,zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zr(*),zau(*),zu(*),wws(*),wwr(*)
 integer*4 :: isct(*),irct(*)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 call HACApK_adot_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp workshare
 zr(1:nd)=zr(1:nd)-zau(1:nd)
!$omp end workshare
 end subroutine HACApK_adotsub_lfmtx_hyp
 
!***HACApK_bicgstab_lfmtx
 subroutine HACApK_bicgstab_lfmtx(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: u(nd),b(nd)
 real*8 :: param(*)
 real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   call MPI_Barrier( icomm, ierr )
   st_measure_time=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx start'
 mstep=param(83)
 eps=param(91)
 allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
 zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
 alpha = 0.0;  beta = 0.0;  zeta = 0.0;
 zz=HACApK_dotp_d(nd, b, b); bnorm=dsqrt(zz);
 zr(:nd)=b(:nd)
 call HACApK_adotsub_lfmtx_p(zr,st_leafmtxp,st_ctl,u,nd)
 zshdw(:nd)=zr(:nd)
 zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm
 if(zrnorm/bnorm<eps) return
! mstep=1
 do in=1,mstep
   zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
   zkp(:nd)=zp(:nd)
   call HACApK_adot_lfmtx_p(zakp,st_leafmtxp,st_ctl,zkp,nd)
! exit
   znom=HACApK_dotp_d(nd,zshdw,zr); zden=HACApK_dotp_d(nd,zshdw,zakp);
   alpha=znom/zden; znomold=znom;
   zt(:nd)=zr(:nd)-alpha*zakp(:nd)
   zkt(:nd)=zt(:nd)
   call HACApK_adot_lfmtx_p(zakt,st_leafmtxp,st_ctl,zkt,nd)
   znom=HACApK_dotp_d(nd,zakt,zt); zden=HACApK_dotp_d(nd,zakt,zakt);
   zeta=znom/zden;
   u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
   zr(:nd)=zt(:nd)-zeta*zakt(:nd)
   beta=alpha/zeta*HACApK_dotp_d(nd,zshdw,zr)/znomold;
   zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
   if(zrnorm/bnorm<eps) exit
 enddo
end subroutine HACApK_bicgstab_lfmtx

!***HACApK_bicgstab_lfmtx_hyp
 subroutine HACApK_bicgstab_lfmtx_hyp(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: u(nd),b(nd)
 real*8 :: param(*)
 real*8,dimension(:),allocatable :: zr,zshdw,zp,zt,zkp,zakp,zkt,zakt
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 integer*4 :: isct(2),irct(2)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   call MPI_Barrier( icomm, ierr )
   st_measure_time=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_bicgstab_lfmtx_hyp start'
 mstep=param(83)
 eps=param(91)
 allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 allocate(zr(nd),zshdw(nd),zp(nd),zt(nd),zkp(nd),zakp(nd),zkt(nd),zakt(nd))
 alpha = 0.0;  beta = 0.0;  zeta = 0.0;
 zz=HACApK_dotp_d(nd, b, b); bnorm=dsqrt(zz);
!$omp parallel
!$omp workshare
 zp(1:nd)=0.0d0; zakp(1:nd)=0.0d0
 zr(:nd)=b(:nd)
!$omp end workshare
 call HACApK_adotsub_lfmtx_hyp(zr,zshdw,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp workshare
 zshdw(:nd)=zr(:nd)
!$omp end workshare
!$omp single
 zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
 if(mpinr==0) print*,'Original relative residual norm =',zrnorm/bnorm
!$omp end single
 do in=1,mstep
   if(zrnorm/bnorm<eps) exit
!$omp workshare
   zp(:nd) =zr(:nd)+beta*(zp(:nd)-zeta*zakp(:nd))
   zkp(:nd)=zp(:nd)
!$omp end workshare
   call HACApK_adot_lfmtx_hyp(zakp,st_leafmtxp,st_ctl,zkp,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp single
   znom=HACApK_dotp_d(nd,zshdw,zr); zden=HACApK_dotp_d(nd,zshdw,zakp);
   alpha=znom/zden; znomold=znom;
!$omp end single
!$omp workshare
   zt(:nd)=zr(:nd)-alpha*zakp(:nd)
   zkt(:nd)=zt(:nd)
!$omp end workshare
   call HACApK_adot_lfmtx_hyp(zakt,st_leafmtxp,st_ctl,zkt,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp single
   znom=HACApK_dotp_d(nd,zakt,zt); zden=HACApK_dotp_d(nd,zakt,zakt);
   zeta=znom/zden;
!$omp end single
!$omp workshare
   u(:nd)=u(:nd)+alpha*zkp(:nd)+zeta*zkt(:nd)
   zr(:nd)=zt(:nd)-zeta*zakt(:nd)
!$omp end workshare
!$omp single
   beta=alpha/zeta*HACApK_dotp_d(nd,zshdw,zr)/znomold;
   zrnorm=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm)
   nstp=in
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
!$omp end single
 enddo
!$omp end parallel
end subroutine HACApK_bicgstab_lfmtx_hyp

!***HACApK_gcrm_lfmtx
 subroutine HACApK_gcrm_lfmtx(st_leafmtxp,st_ctl,st_bemv,u,b,param,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: u(nd),b(nd)
 real*8 :: param(*)
 real*8,dimension(:),allocatable :: zr,zar,capap
 real*8,dimension(:,:),allocatable,target :: zp,zap
 real*8,pointer :: zq(:)
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4 :: isct(2),irct(2)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
   call MPI_Barrier( icomm, ierr )
   st_measure_time=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'gcr_lfmtx_hyp start'
 mstep=param(83)
 mreset=param(87)
 eps=param(91)
 allocate(wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
 allocate(zr(nd),zar(nd),zp(nd,mreset),zap(nd,mreset),capap(mreset))
 alpha = 0.0
 zz=HACApK_dotp_d(nd, b, b); bnorm=dsqrt(zz);
 call HACApK_adot_lfmtx_hyp(zar,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
 zr(:nd)=b(:nd)-zar(:nd)
 zp(:nd,1)=zr(:nd)
 zrnorm2=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm2)
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,0,time,log10(zrnorm/bnorm)
 if(zrnorm/bnorm<eps) return
 call HACApK_adot_lfmtx_hyp(zap(:nd,1),st_leafmtxp,st_ctl,zp(:nd,1),wws,wwr,isct,irct,nd)
 do in=1,mstep
   ik=mod(in-1,mreset)+1
   zq=>zap(:nd,ik)
   znom=HACApK_dotp_d(nd,zq,zr); capap(ik)=HACApK_dotp_d(nd,zq,zq)
   alpha=znom/capap(ik)
   u(:nd)=u(:nd)+alpha*zp(:nd,ik)
   zr(:nd)=zr(:nd)-alpha*zq(:nd)
   zrnomold=zrnorm2
   zrnorm2=HACApK_dotp_d(nd,zr,zr); zrnorm=dsqrt(zrnorm2)
   call MPI_Barrier( icomm, ierr )
   en_measure_time=MPI_Wtime()
   time = en_measure_time - st_measure_time
   if(st_ctl%param(1)>0 .and. mpinr==0) print*,in,time,log10(zrnorm/bnorm)
   if(zrnorm/bnorm<eps .or. in==mstep) exit
   call HACApK_adot_lfmtx_hyp(zar,st_leafmtxp,st_ctl,zr,wws,wwr,isct,irct,nd)
   ikn=mod(in,mreset)+1
   zp(:nd,ikn)=zr(:nd)
   zap(:nd,ikn)=zar(:nd)
   do il=1,ik
     zq=>zap(:nd,il)
     znom=HACApK_dotp_d(nd,zq,zar)
     beta=-znom/capap(il)
     zp(:nd,ikn) =zp(:nd,ikn)+beta*zp(:nd,il)
     zap(:nd,ikn)=zap(:nd,ikn)+beta*zq(:nd)
   enddo
 enddo
 nstp=in
end subroutine

!***HACApK_measurez_time_ax_lfmtx
 subroutine HACApK_measurez_time_ax_lfmtx(st_leafmtxp,st_ctl,nd,nstp,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 integer*4 :: isct(2),irct(2)
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr; param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 allocate(u(nd),b(nd),wws(maxval(lnp(0:nrank-1))),wwr(maxval(lnp(0:nrank-1))))
!$omp parallel private(il)
 do il=1,mstep
   u(:)=1.0; b(:)=1.0
   call HACApK_adot_lfmtx_hyp(u,st_leafmtxp,st_ctl,b,wws,wwr,isct,irct,nd)
 enddo
!$omp end parallel
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_lfmtx

!***HACApK_adot_pmt_lfmtx_p
 integer function HACApK_adot_pmt_lfmtx_p(st_leafmtxp,st_bemv,st_ctl,aww,ww)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: ww(st_bemv%nd),aww(st_bemv%nd)
 real*8,dimension(:),allocatable :: u,au
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:)
 mpinr=st_ctl%lpmd(3); icomm=st_ctl%lpmd(1); nd=st_bemv%nd
 allocate(u(nd),au(nd)); u(:nd)=ww(st_ctl%lod(:nd))
 call MPI_Barrier( icomm, ierr )
 call HACApK_adot_lfmtx_p(au,st_leafmtxp,st_ctl,u,nd)
 aww(st_ctl%lod(:nd))=au(:nd)
 HACApK_adot_pmt_lfmtx_p=lrtrn
end function HACApK_adot_pmt_lfmtx_p

!***HACApK_adot_pmt_lfmtx_hyp
 integer function HACApK_adot_pmt_lfmtx_hyp(st_leafmtxp,st_bemv,st_ctl,aww,ww)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: ww(*),aww(*)
 real*8,dimension(:),allocatable :: u,au,wws,wwr
 integer*4,dimension(:),allocatable :: isct,irct
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:)
 mpinr=st_ctl%lpmd(3); icomm=st_ctl%lpmd(1); nd=st_bemv%nd; nrank=st_ctl%lpmd(2)
 allocate(u(nd),au(nd),isct(2),irct(2)); u(:nd)=ww(st_ctl%lod(:nd))
 allocate(wws(maxval(st_ctl%lnp(:nrank))),wwr(maxval(st_ctl%lnp(:nrank))))
 call MPI_Barrier( icomm, ierr )
!$omp parallel
!$omp barrier
 call HACApK_adot_lfmtx_hyp(au,st_leafmtxp,st_ctl,u,wws,wwr,isct,irct,nd)
!$omp barrier
!$omp end parallel
 call MPI_Barrier( icomm, ierr )
 aww(st_ctl%lod(:nd))=au(:nd)
 HACApK_adot_pmt_lfmtx_hyp=lrtrn
end function HACApK_adot_pmt_lfmtx_hyp

endmodule m_HACApK_solve

