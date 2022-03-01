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
 
!***HACApK_adot_lfmtx_hyp2
 subroutine HACApK_adot_lfmtx_hyp2(zau,st_leafmtxp,st_ctl,zu,wwr,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*),wwr(*)
 integer*4 :: ISTATUS(MPI_STATUS_SIZE)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 real*8,pointer :: time(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 time => st_ctl%time(:)
! allocate(ISTATUS(MPI_STATUS_SIZE))
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 wwr(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(wwr,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 st_time=en_time
 call MPI_Allreduce(wwr, zau, nd, MPI_REAL8, MPI_SUM, icomm, ierr)
 en_time=MPI_Wtime()
 time(2)=time(2)+en_time-st_time
!$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp2

!***HACApK_adot_lfmtx_hyp2_
 subroutine HACApK_adot_lfmtx_hyp2_(zau,st_leafmtxp,st_ctl,zu,wwr,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*),wwr(*)
 integer*4 :: ISTATUS(MPI_STATUS_SIZE)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 real*8,pointer :: time(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 time => st_ctl%time(:)
! allocate(ISTATUS(MPI_STATUS_SIZE))
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 st_time=en_time
 call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icomm, ierr)
 zau(1:nd)=wwr(1:nd)
 en_time=MPI_Wtime()
 time(2)=time(2)+en_time-st_time
!$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp2_
 
!***HACApK_adot_lfmtx_hyp
 subroutine HACApK_adot_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,wws,wwr,isct,irct,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(*),zu(*),wws(*),wwr(*)
 integer*4 :: isct(*),irct(*)
 integer*4 :: ISTATUS(MPI_STATUS_SIZE)
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:)
 real*8,pointer :: time(:)
! integer*4,dimension(:),allocatable :: ISTATUS
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr
 time => st_ctl%time(:)
! allocate(ISTATUS(MPI_STATUS_SIZE))
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master
 
 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
   st_time=en_time
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
   call MPI_Barrier( icomm, ierr )
   en_time=MPI_Wtime()
   time(2)=time(2)+en_time-st_time
 endif
!$omp end master
! stop
 end subroutine HACApK_adot_lfmtx_hyp

!***HACApK_adot_lattice_hyp
 subroutine HACApK_adot_lattice_hyp(st_au,st_LHp,st_ctl,wws,st_u)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_LHp) :: st_LHp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_latticevec) :: st_au,st_u
 real*8 wws(:)
 real*8,allocatable :: zww
 integer*4,pointer :: lpmd(:)
 real*8,pointer :: time(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

!!! print*,'HACApK_adot_lattice_hyp is called'

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpilog=lpmd(4); nrank=lpmd(2); mpinrl=lpmd(37); mpinrt=lpmd(33)
 icommt=lpmd(31); icomml=lpmd(35); icomm=lpmd(1)
 ndlfs=st_LHp%ndlfs; ndtfs=st_LHp%ndtfs

!!! write(mpilog,1000) 'adot_lattice_hyp; ndlfs=',ndlfs,'; ndtfs=',ndtfs
! write(mpilog,1000) '  mpinrl=',mpinrl,'; mpinrt=',mpinrt

!$omp master
 st_time=MPI_Wtime()
!$omp end master

 if(nrank==1)then
 st_au%vs(:)=0.0d0
!$omp parallel
   call HACApK_adot_body_lattice_hyp_tmp(st_au%vs,st_LHp,st_ctl,st_u)
!$omp end parallel
 else
!$omp barrier
   wws(:)=0.0d0
!$omp parallel
   call HACApK_adot_body_lattice_hyp_tmp(wws,st_LHp,st_ctl,st_u)
!$omp end parallel
!$omp barrier

!!! call MPI_Barrier( icomm, ierr )

!$omp master
!   write(mpilog,*) 'wws=',wws(1:3)
   en_time=MPI_Wtime()
   time(1)=time(1)+en_time-st_time
   st_time=MPI_Wtime()
   call MPI_reduce(wws,st_au%vs,ndlfs,MPI_REAL8,MPI_SUM,mpinrl,icommt,ierr) ! The mpinrl coincides with the rank of diagonal in the row of process grid
   en_time=MPI_Wtime()
   time(2)=time(2)+en_time-st_time
   call MPI_bcast(st_au%vs,ndtfs,MPI_REAL8,mpinrt,icomml,ierr)
   en2_time=MPI_Wtime()
   time(3)=time(3)+en2_time-en_time
!   write(mpilog,*) 'st_au%vs=',st_au%vs(1:3)
!$omp end master
 endif

 return

!$omp master
 do itb=1,st_LHp%nlft
   nbgt=st_LHp%lbstrtt(itb)-1
   nbndst=st_LHp%lbndcst(itb)-1
   nditb=st_LHp%lbndcst(itb+1)-st_LHp%lbndcst(itb)
   do i=1, nditb
     ic=i+nbndst
     index = st_u%lodc(ic)
     igt=i+nbgt
     write(mpilog,*) igt,index,st_au%vs(ic)
   enddo
 enddo
!$omp end master

 call MPI_Barrier( icomm, ierr )

! print*,'HACApK_adot_blrmtx_hyp4 end'
! stop
 end subroutine

!***HACApK_adot_blrmtx_hyp4
 subroutine HACApK_adot_blrmtx_hyp4(zau,st_leafmtxp,st_ctl,zu,wws,wwr,nd)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wws(:),wwr(:)
 integer*4,pointer :: lpmd(:),lbl2t(:)
 integer*4,allocatable :: ireqs(:),ireqr(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:),zw2(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

! print*,'HACApK_adot_blrmtx_hyp4 is called'

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:); lbl2t(0:) => st_leafmtxp%lbl2t(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); nrank_l=lpmd(36); mpinrl=lpmd(37); mpinrt=lpmd(33)
 icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft; nlfalt=st_leafmtxp%nlfalt
 ndlfs=st_leafmtxp%ndlfs
 allocate(ireqs(0:nrank_l-1),ireqr(0:nrank_l-1))

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
    st_time=MPI_Wtime()
       ilp=1
       do irl=1,nlfl
         ibl=st_leafmtxp%lnlfl2g(1,irl); ibll=(ibl-1)/nlfalt+1
         irlstrtl=st_leafmtxp%lbstrtl(ibll) ; irlnd=st_leafmtxp%lbndl(ibll)
         wws(ilp:ilp+irlnd-1)=zau(irlstrtl:irlstrtl+irlnd-1)
         ilp=ilp+irlnd
       enddo
!!!     call MPI_Allreduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, icommt, ierr)
     call MPI_reduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, mpinrl, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time

     allocate(zw(nd)); 
     
     is=1
     do irl=0,nrank_l-1
       itag=1; irlnd=st_leafmtxp%lbndlfs(irl)
       if(mpinrl==irl .and. mpinrl==mpinrt) then
         zw(is:is+irlnd-1)=wwr(1:irlnd)
       else
!         print*,'mpinrl=',mpinrl,' ;irl=',irl,' ;is=',is, ' ;irlnd=',irlnd
         if(lbl2t(irl)==1) call MPI_IRECV(zw(is), irlnd, MPI_REAL8, irl, itag, icomml, ireqr(irl),ierr)
       endif
       is=is+irlnd
     enddo
     if(lbl2t(mpinrl)==1)then
       do irl=0,nrank_l-1
         if(mpinrl==irl) cycle
         call MPI_ISEND(wwr, ndlfs, MPI_REAL8, irl, itag, icomml, ireqs(irl),ierr)
       enddo
     endif
     
       ilp=1
       do irl=1,nlfl
         ibl=st_leafmtxp%lnlfl2g(1,irl); ibll=(ibl-1)/nlfalt+1
         irlstrtl=st_leafmtxp%lbstrtl(ibll) ; irlnd=st_leafmtxp%lbndl(ibll)
         zau(irlstrtl:irlstrtl+irlnd-1)=wwr(ilp:ilp+irlnd-1)
         ilp=ilp+irlnd
       enddo
     
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
       if(lbl2t(irl)==0) cycle
       call MPI_WAIT(ireqr(irl),ISTATUS,ierr)
     enddo
       
     is=1
     do icl=1,nrank_l; irl=icl-1
       do ibl=0,nlfl
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         if(mpinrl/=irl .and. lbl2t(irl)==1) zau(ips:ipe)=zw(is:ie)
         is=ie+1
       enddo
     enddo

     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
 endif
!$omp end master
! print*,'HACApK_adot_blrmtx_hyp4 end'
! stop
 end subroutine

!***HACApK_adot_blrmtx_hyp31
 subroutine HACApK_adot_blrmtx_hyp31(zau,st_leafmtxp,st_ctl,zu,wws,wwr,nd)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wws(:),wwr(:)
 integer*4,pointer :: lpmd(:),lbl2t(:)
 integer*4,allocatable :: ireqs(:),ireqr(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:),zw2(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

! print*,'HACApK_adot_blrmtx_hyp31 is called'

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:); lbl2t(0:) => st_leafmtxp%lbl2t(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); nrank_l=lpmd(36); mpinrl=lpmd(37)
 icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft; nlfalt=st_leafmtxp%nlfalt
 ndlfs=st_leafmtxp%ndlfs
 allocate(ireqs(0:nrank_l-1),ireqr(0:nrank_l-1))

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
    st_time=MPI_Wtime()
   if(st_ctl%param(41)==1)then
     call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time
     zau(1:nd)=wwr(1:nd)
     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
   else
       ilp=1
       do irl=1,nlfl
         ibl=st_leafmtxp%lnlfl2g(1,irl); ibll=(ibl-1)/nlfalt+1
         irlstrtl=st_leafmtxp%lbstrtl(ibll) ; irlnd=st_leafmtxp%lbndl(ibll)
         wws(ilp:ilp+irlnd-1)=zau(irlstrtl:irlstrtl+irlnd-1)
         ilp=ilp+irlnd
       enddo
     call MPI_Allreduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time

     allocate(zw(nd)); 
     is=1
     do irl=0,nrank_l-1
       itag=1; irlnd=st_leafmtxp%lbndlfs(irl)
       if(mpinrl==irl) then
         zw(is:is+irlnd-1)=wwr(1:irlnd)
       else
!         print*,'mpinrl=',mpinrl,' ;irl=',irl,' ;is=',is, ' ;irlnd=',irlnd
         if(lbl2t(irl)==1) call MPI_IRECV(zw(is), irlnd, MPI_REAL8, irl, itag, icomml, ireqr(irl),ierr)
       endif
       is=is+irlnd
     enddo
     if(lbl2t(mpinrl)==1)then
       do irl=0,nrank_l-1
         if(mpinrl==irl) cycle
         call MPI_ISEND(wwr, ndlfs, MPI_REAL8, irl, itag, icomml, ireqs(irl),ierr)
       enddo
     endif
     
       ilp=1
       do irl=1,nlfl
         ibl=st_leafmtxp%lnlfl2g(1,irl); ibll=(ibl-1)/nlfalt+1
         irlstrtl=st_leafmtxp%lbstrtl(ibll) ; irlnd=st_leafmtxp%lbndl(ibll)
         zau(irlstrtl:irlstrtl+irlnd-1)=wwr(ilp:ilp+irlnd-1)
         ilp=ilp+irlnd
       enddo
     
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
       if(lbl2t(irl)==0) cycle
       call MPI_WAIT(ireqr(irl),ISTATUS,ierr)
     enddo
       
     is=1
     do icl=1,nrank_l; irl=icl-1
       do ibl=0,nlfl
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         if(mpinrl/=irl .and. lbl2t(irl)==1) zau(ips:ipe)=zw(is:ie)
         is=ie+1
       enddo
     enddo

     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
  endif
 endif
!$omp end master
! stop
 end subroutine

!***HACApK_adot_blrmtx_hyp3
 subroutine HACApK_adot_blrmtx_hyp3(zau,st_leafmtxp,st_ctl,zu,wws,wwr,nd)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wws(:),wwr(:)
 integer*4,pointer :: lpmd(:),lbl2t(:)
 integer*4,allocatable :: ireqs(:),ireqr(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:),zw2(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 print*,'HACApK_adot_blrmtx_hyp3 is called'

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:); lbl2t(0:) => st_leafmtxp%lbl2t(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); nrank_l=lpmd(36); mpinrl=lpmd(37)
 icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft; nlfalt=st_leafmtxp%nlfalt
 ndlfs=st_leafmtxp%ndlfs
 allocate(ireqs(0:nrank_l-1),ireqr(0:nrank_l-1))

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
    st_time=MPI_Wtime()
   if(st_ctl%param(41)==1)then
     call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time
     zau(1:nd)=wwr(1:nd)
     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
   else
     ilp=1
     do ilf=1,nlfl
       ip=(ilf-1)*nlft+1
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       wws(ilp:ilp+ndl-1)=zau(nstrtl:nstrtl+ndl-1)
       ilp=ilp+ndl
     enddo
     call MPI_Allreduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time

     allocate(zw(nd)); 
     is=1
     do irl=0,nrank_l-1
       itag=1; irlnd=st_leafmtxp%lbndlfs(irl)
       if(mpinrl==irl) then
         zw(is:is+irlnd-1)=wwr(1:irlnd)
       else
!         print*,'mpinrl=',mpinrl,' ;irl=',irl,' ;is=',is, ' ;irlnd=',irlnd
         if(lbl2t(irl)==1) call MPI_IRECV(zw(is), irlnd, MPI_REAL8, irl, itag, icomml, ireqr(irl),ierr)
       endif
       is=is+irlnd
     enddo
     if(lbl2t(mpinrl)==1)then
       do irl=0,nrank_l-1
         if(mpinrl==irl) cycle
         call MPI_ISEND(wwr, ndlfs, MPI_REAL8, irl, itag, icomml, ireqs(irl),ierr)
       enddo
     endif
     nlfla=nlfalt/nrank_l
     is=1; icl=mpinrl+1
       do ibl=0,nlfla
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         zau(ips:ipe)=wwr(is:ie)
         is=ie+1
       enddo
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
       if(lbl2t(irl)==0) cycle
       call MPI_WAIT(ireqr(irl),ISTATUS,ierr)
     enddo
       
     nlfla=nlfalt/nrank_l
     is=1
     do icl=1,nrank_l; irl=icl-1
       do ibl=0,nlfla
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         if(mpinrl/=irl .and. lbl2t(irl)==1) zau(ips:ipe)=zw(is:ie)
         is=ie+1
       enddo
     enddo

     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
  endif
 endif
!$omp end master
! stop
 end subroutine

!***HACApK_adot_blrmtx_hyp22
 subroutine HACApK_adot_blrmtx_hyp22(zau,st_leafmtxp,st_ctl,zu,wwr,nd)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wwr(:)
 integer*4,pointer :: lpmd(:)
 integer*4,allocatable :: ireqs(:),ireqr(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:),zw2(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); nrank_l=lpmd(36); mpinrl=lpmd(37)
 icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft; nlfalt=st_leafmtxp%nlfalt
 ndlfs=st_leafmtxp%ndlfs
 allocate(ireqs(0:nrank_l-1),ireqr(0:nrank_l-1))

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
!     st_time=en_time
     call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icommt, ierr)
!     en_time=MPI_Wtime()
! call MPI_Barrier( icomm, ierr )
!     time(2)=time(2)+en_time-st_time
     call MPI_Allreduce(wwr, zau, nd, MPI_REAL8, MPI_SUM, icomml, ierr)
!     en2_time=MPI_Wtime()
!     time(3)=time(3)+en2_time-en_time
!$omp end master
! stop
 end subroutine

!***HACApK_adot_blrmtx_hyp21
 subroutine HACApK_adot_blrmtx_hyp21(zau,st_leafmtxp,st_ctl,zu,wws,wwr,nd)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wws(:),wwr(:)
 integer*4,pointer :: lpmd(:)
 integer*4,allocatable :: ireqs(:),ireqr(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:),zw2(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); nrank_l=lpmd(36); mpinrl=lpmd(37)
 icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft; nlfalt=st_leafmtxp%nlfalt
 ndlfs=st_leafmtxp%ndlfs
 allocate(ireqs(0:nrank_l-1),ireqr(0:nrank_l-1))

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
    st_time=MPI_Wtime()
   if(st_ctl%param(41)==1)then
     call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time
     zau(1:nd)=wwr(1:nd)
     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
   else
       ilp=1
       do irl=1,nlfl
         ibl=st_leafmtxp%lnlfl2g(1,irl); ibll=(ibl-1)/nlfalt+1
         irlstrtl=st_leafmtxp%lbstrtl(ibll) ; irlnd=st_leafmtxp%lbndl(ibll)
         wws(ilp:ilp+irlnd-1)=zau(irlstrtl:irlstrtl+irlnd-1)
         ilp=ilp+irlnd
       enddo
     call MPI_Allreduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time

     allocate(zw(nd)); 
     is=1
     do irl=0,nrank_l-1
       itag=1; irlnd=st_leafmtxp%lbndlfs(irl)
       if(mpinrl==irl) then
         zw(is:is+irlnd-1)=wwr(1:irlnd)
       else
!         print*,'mpinrl=',mpinrl,' ;irl=',irl,' ;is=',is, ' ;irlnd=',irlnd
         call MPI_IRECV(zw(is), irlnd, MPI_REAL8, irl, itag, icomml, ireqr(irl),ierr)
       endif
       is=is+irlnd
     enddo
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
       call MPI_ISEND(wwr, ndlfs, MPI_REAL8, irl, itag, icomml, ireqs(irl),ierr)
     enddo
     
     !! my color
       ilp=1
       do irl=1,nlfl
         ibl=st_leafmtxp%lnlfl2g(1,irl); ibll=(ibl-1)/nlfalt+1
         irlstrtl=st_leafmtxp%lbstrtl(ibll) ; irlnd=st_leafmtxp%lbndl(ibll)
         zau(irlstrtl:irlstrtl+irlnd-1)=wwr(ilp:ilp+irlnd-1)
         ilp=ilp+irlnd
       enddo
       
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
       call MPI_WAIT(ireqr(irl),ISTATUS,ierr)
     enddo
       
     is=1
     do icl=1,nrank_l
       do ibl=0,nlfl
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         if(mpinrl/=icl-1) zau(ips:ipe)=zw(is:ie)
         is=ie+1
       enddo
     enddo

     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
  endif
 endif
!$omp end master
! stop
 end subroutine

!***HACApK_adot_blrmtx_hyp2
 subroutine HACApK_adot_blrmtx_hyp2(zau,st_leafmtxp,st_ctl,zu,wws,wwr,nd)
 include 'mpif.h'
 integer ISTATUS(MPI_STATUS_SIZE)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wws(:),wwr(:)
 integer*4,pointer :: lpmd(:)
 integer*4,allocatable :: ireqs(:),ireqr(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:),zw2(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); nrank_l=lpmd(36); mpinrl=lpmd(37)
 icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft; nlfalt=st_leafmtxp%nlfalt
 ndlfs=st_leafmtxp%ndlfs
 allocate(ireqs(0:nrank_l-1),ireqr(0:nrank_l-1))

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
    st_time=MPI_Wtime()
   if(st_ctl%param(41)==1)then
     call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time
     zau(1:nd)=wwr(1:nd)
     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
   else
     ilp=1
     do ilf=1,nlfl
       ip=(ilf-1)*nlft+1
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       wws(ilp:ilp+ndl-1)=zau(nstrtl:nstrtl+ndl-1)
       ilp=ilp+ndl
     enddo
     call MPI_Allreduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time

     allocate(zw(nd)); 
     is=1
     do irl=0,nrank_l-1
       itag=1; irlnd=st_leafmtxp%lbndlfs(irl)
       if(mpinrl==irl) then
         zw(is:is+irlnd-1)=wwr(1:irlnd)
       else
!         print*,'mpinrl=',mpinrl,' ;irl=',irl,' ;is=',is, ' ;irlnd=',irlnd
         call MPI_IRECV(zw(is), irlnd, MPI_REAL8, irl, itag, icomml, ireqr(irl),ierr)
       endif
       is=is+irlnd
     enddo
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
       call MPI_ISEND(wwr, ndlfs, MPI_REAL8, irl, itag, icomml, ireqs(irl),ierr)
     enddo
     nlfla=nlfalt/nrank_l
     is=1; icl=mpinrl+1
       do ibl=0,nlfla
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         zau(ips:ipe)=wwr(is:ie)
         is=ie+1
       enddo
     do irl=0,nrank_l-1
       if(mpinrl==irl) cycle
!       call MPI_WAIT(ireqs(irl),ISTATUS,ierr)
       call MPI_WAIT(ireqr(irl),ISTATUS,ierr)
     enddo
     
!     write(mpilog,*) 'zw='
!     write(mpilog,*) zw
!     write(mpilog,*) 'wwr='
!     write(mpilog,*) wwr(1:ndlfs)
     
     nlfla=nlfalt/nrank_l
     is=1
     do icl=1,nrank_l
       do ibl=0,nlfla
         ibpl=ibl*nrank_l+icl; if(ibpl>nlfalt) exit
         ips=st_leafmtxp%lbstrtl(ibpl)
         ipe=st_leafmtxp%lbstrtl(ibpl+1)-1
         ie=is+ipe-ips
         if(mpinrl/=icl-1) zau(ips:ipe)=zw(is:ie)
! if(mpinr==0) write(*,1000) 'ibpl=',ibpl,' ;ips=',ips,' ;ipe=',ipe,' ;is=',is,' ;ie=',ie
         is=ie+1
       enddo
     enddo

     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
  endif
 endif
!$omp end master
! stop
 end subroutine

!***HACApK_adot_body_lattice
 subroutine HACApK_adot_body_lattice(zau,st_leafmtxp,st_ctl,zu,nd)
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
 end subroutine

!***HACApK_adot_blrmtx_hyp
 subroutine HACApK_adot_blrmtx_hyp(zau,st_leafmtxp,st_ctl,zu,wws,wwr,nd)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: zau(:),zu(:),wws(:),wwr(:)
 integer*4,pointer :: lpmd(:)
 real*8,pointer :: time(:)
 real*8, allocatable :: zw(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); icommt=lpmd(31); icomml=lpmd(35); 
 nlf=st_leafmtxp%nlf; nlfl=st_leafmtxp%nlfl; nlft=st_leafmtxp%nlft
 ndlfs=st_leafmtxp%ndlfs

!$omp master
 call MPI_Barrier( icomm, ierr )
 st_time=MPI_Wtime()
!$omp end master

 zau(:nd)=0.0d0
!$omp barrier
 call HACApK_adot_body_lfmtx_hyp(zau,st_leafmtxp,st_ctl,zu,nd)
!$omp barrier
!$omp master
 call MPI_Barrier( icomm, ierr )
 en_time=MPI_Wtime()
 time(1)=time(1)+en_time-st_time
 if(nrank>1)then
    st_time=MPI_Wtime()
   if(st_ctl%param(41)==1)then
     call MPI_Allreduce(zau, wwr, nd, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time
     zau(1:nd)=wwr(1:nd)
     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
   else
     ilp=1
     do ilf=1,nlfl
       ip=(ilf-1)*nlft+1
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       wws(ilp:ilp+ndl-1)=zau(nstrtl:nstrtl+ndl-1)
       ilp=ilp+ndl
     enddo
     call MPI_Allreduce(wws, wwr, ndlfs, MPI_REAL8, MPI_SUM, icommt, ierr)
     en_time=MPI_Wtime()
     time(2)=time(2)+en_time-st_time

     allocate(zw(nd)); zw(:nd)=0.0d0
     ilp=1
     do ilf=1,nlfl
       ip=(ilf-1)*nlft+1
       ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
       nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
       zw(nstrtl:nstrtl+ndl-1)=wwr(ilp:ilp+ndl-1)
       ilp=ilp+ndl
     enddo
     call MPI_Allreduce(zw, zau, nd, MPI_REAL8, MPI_SUM, icomml, ierr)
     en2_time=MPI_Wtime()
     time(3)=time(3)+en2_time-en_time
  endif
 endif
!$omp end master
! stop
 end subroutine

!***HACApK_adot_body_lattice_hyp_tmp
 subroutine HACApK_adot_body_lattice_hyp_tmp(zau,st_LHp,st_ctl,st_u)
 type(st_HACApK_LHp) :: st_LHp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_latticevec) :: st_u
 type(st_HACApK_leafmtx),pointer :: st_lf
 real*8 :: zau(:)
 real*8,dimension(:),allocatable :: zbu
 real*8,pointer :: zup(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); nthr=lpmd(20)
 nlfl=st_LHp%nlfl; nlft=st_LHp%nlft
 nthr = omp_get_num_threads(); nlflc=nlfl/nthr
 ith = omp_get_thread_num(); ith1 = ith+1
 if(nthr==ith1)then
   nlthe=nlfl
 else
   nlthe=nlflc*ith1
 endif
!#!#!$omp critical
!#!#   print*, 'sub HACApK_adot_body_lattice_hyp_tmp; ith=',ith,'; nlthe=',nlthe,'; nlfl=',nlfl
!#!#!$omp end critical
 
 zup=>st_u%vs
 do ilb=nlflc*ith+1,nlthe; nbndsl=st_LHp%lbndcsl(ilb)-1
   do itb=1,nlft; nbndst=st_LHp%lbndcst(itb)-1
     do ip=1,st_LHp%st_lfp(ilb,itb)%nlf
       st_lf=>st_LHp%st_lfp(ilb,itb)%st_lf(ip)
       ndl   =st_lf%ndl   ; ndt   =st_lf%ndt
       nstrtl=st_lf%nstrtl-1; nstrtt=st_lf%nstrtt-1
!!! write(mpilog,1000) '      nstrtl=',nstrtl,'; nstrtt=',nstrtt,'; ltmtx=',st_lf%ltmtx
       if(st_lf%ltmtx==1)then
         kt=st_lf%kt
         allocate(zbu(kt)); zbu(:)=0.0d0
         do il=1,kt
           do it=1,ndt; itt=it+nstrtt+nbndst
             zbu(il)=zbu(il)+st_lf%a1(it,il)*zup(itt)
           enddo
         enddo
         do il=1,kt
           do it=1,ndl; ill=it+nstrtl+nbndsl
             zau(ill)=zau(ill)+st_lf%a2(it,il)*zbu(il)
           enddo
         enddo
         deallocate(zbu)
       elseif(st_lf%ltmtx==2)then
         do il=1,ndl; ill=il+nstrtl+nbndsl
           do it=1,ndt; itt=it+nstrtt+nbndst
             zau(ill)=zau(ill)+st_lf%a1(it,il)*zup(itt)
           enddo
         enddo
       endif
     enddo
   enddo
 enddo
 end subroutine

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
 ith = omp_get_thread_num(); ith1 = ith+1
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

!***HACApK_measurez_time_ax_lfmtx_o
 subroutine HACApK_measurez_time_ax_lfmtx_o(st_leafmtxp,st_ctl,nd,lrtrn)
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
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
 do il=1,mstep
!$omp workshare
   u(:)=1.0; b(:)=1.0
!$omp end workshare
   call HACApK_adot_lfmtx_hyp(u,st_leafmtxp,st_ctl,b,wws,wwr,isct,irct,nd)
 enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'lfmtx; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'lfmtx;         body  =',ztime_body/mstep
     write(6,2000) 'lfmtx; mpi0_t_reduce  =',st_ctl%time(2)/(mstep*nit)
   endif
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_lfmtx_o

!***HACApK_measurez_time_ax_lfmtx
 subroutine HACApK_measurez_time_ax_lfmtx(st_leafmtxp,st_ctl,nd,lrtrn)
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
 allocate(u(nd),b(nd),wwr(nd))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
 do il=1,mstep
!$omp workshare
   u(:)=1.0; b(:)=1.0
!$omp end workshare
   call HACApK_adot_lfmtx_hyp2(u,st_leafmtxp,st_ctl,b,wwr,nd)
 enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'lfmtx; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'lfmtx;         body  =',ztime_body/mstep
     write(6,2000) 'lfmtx;  t_allreduce  =',st_ctl%time(2)/(mstep*nit)
   endif
 deallocate(wwr)
end subroutine

!***HACApK_measurez_time_ax_blrmtx
 subroutine HACApK_measurez_time_ax_blrmtx(st_leafmtxp,st_ctl,nd,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 ndlfs=st_leafmtxp%ndlfs
 allocate(u(nd),b(nd),wws(ndlfs),wwr(ndlfs))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
   do il=1,mstep
!$omp workshare
     u(:)=1.0; b(:)=1.0
!$omp end workshare
     call HACApK_adot_blrmtx_hyp(u,st_leafmtxp,st_ctl,b,wws,wwr,nd)
   enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'blrmtx; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'blrmtx;         body  =',ztime_body/mstep
     write(6,2000) 'blrmtx; mpi0_t_reduce  =',st_ctl%time(2)/(mstep*nit)
     write(6,2000) 'blrmtx; mpi0_l_reduce  =',st_ctl%time(3)/(mstep*nit)
   endif
 deallocate(wws,wwr)
end subroutine

!***HACApK_measurez_time_ax_blrmtx2
 subroutine HACApK_measurez_time_ax_blrmtx2(st_leafmtxp,st_ctl,nd,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 ndlfs=st_leafmtxp%ndlfs
 allocate(u(nd),b(nd),wws(ndlfs),wwr(ndlfs))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
   do il=1,mstep
!$omp workshare
     u(:)=1.0; b(:)=1.0
!$omp end workshare
     call HACApK_adot_blrmtx_hyp21(u,st_leafmtxp,st_ctl,b,wws,wwr,nd)
   enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'blrmtx2; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'blrmtx2;         body  =',ztime_body/mstep
     write(6,2000) 'blrmtx2; mpi0_t_reduce  =',st_ctl%time(2)/(mstep*nit)
     write(6,2000) 'blrmtx2; mpi0_l_reduce  =',st_ctl%time(3)/(mstep*nit)
   endif
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_blrmtx2

!***HACApK_measurez_time_ax_blrmtx22
 subroutine HACApK_measurez_time_ax_blrmtx22(st_leafmtxp,st_ctl,nd,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 ndlfs=st_leafmtxp%ndlfs
 allocate(u(nd),b(nd),wws(nd),wwr(nd))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
   do il=1,mstep
!$omp workshare
     u(:)=1.0; b(:)=1.0
!$omp end workshare
     call HACApK_adot_blrmtx_hyp22(u,st_leafmtxp,st_ctl,b,wwr,nd)
   enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'blrmtx22; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'blrmtx22;         body  =',ztime_body/mstep
     write(6,2000) 'blrmtx22; mpi0_t_reduce  =',st_ctl%time(2)/(mstep*nit)
     write(6,2000) 'blrmtx22; mpi0_l_reduce  =',st_ctl%time(3)/(mstep*nit)
   endif
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_blrmtx22

!***HACApK_measurez_time_ax_blrmtx3
 subroutine HACApK_measurez_time_ax_blrmtx3(st_leafmtxp,st_ctl,nd,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 ndlfs=st_leafmtxp%ndlfs
 allocate(u(nd),b(nd),wws(ndlfs),wwr(ndlfs))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
   do il=1,mstep
!$omp workshare
     u(:)=1.0; b(:)=1.0
!$omp end workshare
     call HACApK_adot_blrmtx_hyp31(u,st_leafmtxp,st_ctl,b,wws,wwr,nd)
   enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'blrmtx3; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'blrmtx3;         body  =',ztime_body/mstep
     write(6,2000) 'blrmtx3; mpi0_t_reduce  =',st_ctl%time(2)/(mstep*nit)
     write(6,2000) 'blrmtx3; mpi0_l_reduce  =',st_ctl%time(3)/(mstep*nit)
   endif
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_blrmtx3

!***HACApK_measurez_time_ax_blrmtx4
 subroutine HACApK_measurez_time_ax_blrmtx4(st_leafmtxp,st_ctl,nd,lrtrn)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 real*8,dimension(:),allocatable :: wws,wwr,u,b
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lpmd => st_ctl%lpmd(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 mstep=param(99)
 ndlfs=st_leafmtxp%ndlfs
 allocate(u(nd),b(nd),wws(ndlfs),wwr(ndlfs))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 nit=param(98)
 do it=1,nit
   call MPI_Barrier( icomm, ierr )
   st_measure_time_ax=MPI_Wtime()
!$omp parallel private(il)
   do il=1,mstep
!$omp workshare
     u(:)=1.0; b(:)=1.0
!$omp end workshare
     call HACApK_adot_blrmtx_hyp4(u,st_leafmtxp,st_ctl,b,wws,wwr,nd)
   enddo
!$omp end parallel
   call MPI_Barrier( icomm, ierr )
   en_measure_time_ax=MPI_Wtime()
   ztime=en_measure_time_ax - st_measure_time_ax
   ztime_ax=min(ztime,ztime_ax)
   ztime_body=min(st_ctl%time(1),ztime_body)
 enddo
   if(st_ctl%param(1)>0 .and. mpinr==0) then
     write(6,2000) 'blrmtx4; time_AX_once  =',ztime_ax/mstep
     write(6,2000) 'blrmtx4;         body  =',ztime_body/mstep
     write(6,2000) 'blrmtx4; mpi0_t_reduce  =',st_ctl%time(2)/(mstep*nit)
     write(6,2000) 'blrmtx4; mpi0_l_reduce  =',st_ctl%time(3)/(mstep*nit)
   endif
 deallocate(wws,wwr)
end subroutine HACApK_measurez_time_ax_blrmtx4

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
 real*8,dimension(:),allocatable :: u,au,wwr,wws
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:)
 mpinr=st_ctl%lpmd(3); icomm=st_ctl%lpmd(1); nd=st_bemv%nd; nrank=st_ctl%lpmd(2)
 allocate(u(nd),au(nd)); u(:nd)=ww(st_ctl%lod(:nd))
 allocate(wwr(nd),wws(nd))
 call MPI_Barrier( icomm, ierr )
!$omp parallel
!$omp barrier
 if(st_ctl%param(8)==10 .or. st_ctl%param(8)==20)then
     call HACApK_adot_lfmtx_hyp2(au,st_leafmtxp,st_ctl,u,wwr,nd)
!    call HACApK_adot_blrmtx_hyp2(au,st_leafmtxp,st_ctl,u,wws,wwr,nd)
!   call HACApK_adot_blrmtx_hyp22(au,st_leafmtxp,st_ctl,u,wwr,nd)
 else
     call HACApK_adot_lfmtx_hyp2(au,st_leafmtxp,st_ctl,u,wwr,nd)
 endif
!$omp barrier
!$omp end parallel
 call MPI_Barrier( icomm, ierr )
 aww(st_ctl%lod(:nd))=au(:nd)
 HACApK_adot_pmt_lfmtx_hyp=lrtrn
end function HACApK_adot_pmt_lfmtx_hyp

!***HACApK_adot_pmt_lfmtx_hyp_o
 integer function HACApK_adot_pmt_lfmtx_hyp_o(st_leafmtxp,st_bemv,st_ctl,aww,ww)
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
 HACApK_adot_pmt_lfmtx_hyp_o=lrtrn
end function HACApK_adot_pmt_lfmtx_hyp_o

!***HACApK_adot_pmt_blr_hyp
 integer function HACApK_adot_pmt_blr_hyp(st_leafmtxp,st_bemv,st_ctl,aww,ww)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: ww(*),aww(*)
 real*8,dimension(:),allocatable :: u,au,wws,wwr
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 lpmd => st_ctl%lpmd(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 nd=st_bemv%nd
 ndlfs=st_leafmtxp%ndlfs
 allocate(u(nd),au(nd)); u(:nd)=ww(st_ctl%lod(:nd))
 allocate(wws(ndlfs),wwr(ndlfs))
 st_ctl%time(:)=0.0d0; ztime_ax=1.0d10; ztime_body=1.0d10
 call MPI_Barrier( icomm, ierr )
!$omp parallel
!$omp barrier
 call HACApK_adot_blrmtx_hyp2(au,st_leafmtxp,st_ctl,u,wws,wwr,nd)
!$omp barrier
!$omp end parallel
 call MPI_Barrier( icomm, ierr )
 aww(st_ctl%lod(:nd))=au(:nd)
 HACApK_adot_pmt_blr_hyp=lrtrn
end function HACApK_adot_pmt_blr_hyp

endmodule m_HACApK_solve

