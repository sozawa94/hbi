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
!C**************************************************************************
!C  This file includes examples of integrating routines for H-matrices
!C  created by Akihiro Ida at Kyoto University on May 2012
!C  last modified by Akihiro Ida on Sep. 2014
!C**************************************************************************
module m_HACApK_use
 use m_HACApK_solve
 use m_HACApK_base
 implicit real*8(a-h,o-z)
contains

!*** HACApK_gensolv
 integer function HACApK_gensolv(st_leafmtxp,st_bemv,st_ctl,gmid,rhs,sol,ztol)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: gmid(st_bemv%nd,3),rhs(st_bemv%nd),sol(st_bemv%nd),ztol
 real*8 :: ww(st_bemv%nd),aww(st_bemv%nd),aww1(st_bemv%nd)
 1000 format(5(a,i10)/)
 2000 format(5(a,1pe15.8)/)
 
 mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
 icomm=st_ctl%lpmd(1)
 lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,gmid,ztol)
 call MPI_Barrier( icomm, ierr )
 lrtrn=HACApK_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
 call MPI_Barrier( icomm, ierr )
 st_measure_time_ax=MPI_Wtime()
 call HACApK_measurez_time_ax_lfmtx(st_leafmtxp,st_ctl,st_bemv%nd,nstp,lrtrn)
 en_measure_time_ax=MPI_Wtime()
 if(st_ctl%param(1)>0 .and. mpinr==0)  write(6,2000) 'lfmtx; time_AX_once  =',(en_measure_time_ax - st_measure_time_ax)/st_ctl%param(99)
9999 continue
 HACApK_gensolv=lrtrn
 endfunction

!*** HACApK_generate
 integer function HACApK_generate(st_leafmtxp,st_bemv,st_ctl,coord,ztol)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: coord(st_bemv%nd,*)
 integer*8 :: mem8,nth1_mem,imem
 integer*4 :: ierr
 integer*4,dimension(:),allocatable :: lnmtx(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,1pe15.8)/)
 
 lrtrn=0
 nofc=st_bemv%nd; nffc=1; ndim=3
 mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
 st_ctl%param(71)=ztol
 
 call HACApK_chk_st_ctl(st_ctl)
 
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'***************** HACApK start ********************'
 if(st_ctl%param(1)>0)  write(mpilog,1000) 'irank=',mpinr,', nrank=',nrank
 nd=nofc*nffc
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nd=',nd,' nofc=',nofc,' nffc=',nffc
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nrank=',nrank,' nth=',nthr
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'param:'
 if(st_ctl%param(1)>1 .and. mpinr==0) write(*,10000) st_ctl%param(1:100)
 10000 format(10(1pe10.3)/)
 allocate(lnmtx(3))
 call MPI_Barrier( icomm, ierr )
 st_s=MPI_Wtime()

 call HACApK_generate_frame_leafmtx(st_leafmtxp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)

 call MPI_Barrier( icomm, ierr )
 st_create_hmtx=MPI_Wtime()
 st_bemv%lp61=0
 if(st_ctl%param(61)==2)then
   call HACApK_cal_matnorm(znrm2,st_bemv,st_ctl%lpmd,nd)
   call MPI_Barrier( icomm, ierr )
   call MPI_Allreduce( znrm2, znrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr );
   znrm=dsqrt(znrm)/nd
!   print*,'irank=',mpinr,'znrm2=',znrm2,' znrm=',znrm
 elseif(st_ctl%param(61)==3)then
   ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)
   allocate(st_bemv%ao(nd)); st_bemv%ao(:)=0.0d0; zsqnd=sqrt(real(nd))
   do il=ndnr_s,ndnr_e
     zad=HACApK_entry_ij(il,il,st_bemv)
     st_bemv%ao(il)=1.0d0/dsqrt(zad/zsqnd)
   enddo
   call MPI_Barrier( icomm, ierr )
   call HACApK_impi_allgv(st_bemv%ao,st_ctl%lpmd,nd)
!   call MPI_Barrier( icomm, ierr )
   znrm=1.0/nd
   st_bemv%lp61=3
 else
   znrm=0.0d0
 endif
 call MPI_Barrier( icomm, ierr )
 st_cal_matnorm=MPI_Wtime()
 if(st_ctl%param(1)>0)  write(mpilog,1000) 'ndnr_s=',st_ctl%lpmd(6),', ndnr_e=',st_ctl%lpmd(7),', ndnr=',st_ctl%lpmd(5)
 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,' ndlf_s=',st_ctl%lpmd(11),', ndlf_e=',st_ctl%lpmd(12),', ndlf=',st_leafmtxp%nlf
 lnps=nd+1; lnpe=0
 if(st_leafmtxp%nlf<1)then
   print*,'ERROR!; sub HACApK_generate; irank=',mpinr,' nlf=',st_leafmtxp%nlf
 endif
 if(st_ctl%param(10)==0) return
 call HACApK_fill_leafmtx_hyp(st_leafmtxp%st_lf,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe,st_ctl%lthr)
! call HACApK_fill_leafmtx(st_leafmtxp%st_lf,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe)
 call MPI_Barrier( icomm, ierr )
 if(st_ctl%param(1)>0)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:3)
 ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)

 st_fill_hmtx=MPI_Wtime()
 if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_supermatrix             =',st_create_hmtx- st_s
 if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_fill_hmtx               =',st_fill_hmtx-st_cal_matnorm
 if(st_ctl%param(1)>0)  write(mpilog,2000)  'time_construction_Hmatrix    =',st_fill_hmtx-st_s

 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_supermatrix             =',st_create_hmtx - st_s
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_fill_hmtx               =',st_fill_hmtx - st_cal_matnorm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_construction_Hmatrix    =',st_fill_hmtx - st_s

 call MPI_Barrier( icomm, ierr )

 call HACApK_chk_leafmtx(st_leafmtxp,st_ctl,lnmtx,nd,mem8)

 ktp=0
 call HACApK_setcutthread(st_ctl%lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
      
 call MPI_Barrier( icomm, ierr )
! print*,'mpinr=',mpinr,lnps,lnpe
 st_ctl%lnp(mpinr+1)=lnpe-lnps
 call MPI_Barrier( icomm, ierr )
 call MPI_Allgather(lnpe-lnps,1,MPI_INTEGER,st_ctl%lnp,1, MPI_INTEGER, icomm, ierr )
 st_ctl%lsp(mpinr+1)=lnps
 call MPI_Allgather(lnps,1,MPI_INTEGER,st_ctl%lsp,1, MPI_INTEGER, icomm, ierr )
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'lnp=',st_ctl%lnp(:nrank)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'lsp=',st_ctl%lsp(:nrank)
 
 if(st_ctl%param(11)/=0) then
   call MPI_Barrier( icomm, ierr )
   call HACApK_accuracy_leafmtx(st_leafmtxp,st_bemv,st_ctl,st_ctl%lod,st_ctl%lod,st_ctl%lpmd,nofc,nffc)
 endif
9999 continue
 HACApK_generate=lrtrn
 endfunction

!*** HACApK_solve
 integer function HACApK_solve(st_leafmtxp,st_bemv,st_ctl,rhs,sol,ztol)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: rhs(st_bemv%nd),sol(st_bemv%nd),ztol
 real*8,pointer :: param(:)
 real*8,dimension(:),allocatable :: u,b,www,ao
 integer*4,pointer :: lpmd(:),lnp(:),lsp(:),lthr(:),lod(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,1pe15.8)/)

 lpmd => st_ctl%lpmd(:); lnp(0:) => st_ctl%lnp; lsp(0:) => st_ctl%lsp;lthr(0:) => st_ctl%lthr;lod => st_ctl%lod(:); param=>st_ctl%param(:)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); nthr=lpmd(20)
 param(91)=ztol
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'HACApK_solve start'
 nofc=st_bemv%nd;nffc=1;ndim=3
 nd=nofc*nffc
 if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,' lthr=',lthr(0:nthr-1)
 allocate(u(nd),b(nd)); u(:nd)=sol(lod(:nd)); b(:nd)=rhs(lod(:nd))
 if(param(61)==3)then
!   do il=ndnr_s,ndnr_e
   do il=1,nd
     u(il)=u(il)/st_bemv%ao(lod(il))
     b(il)=b(il)*st_bemv%ao(lod(il))
   enddo
 endif
 if(param(83)>0) then
   allocate(ao(nd))
   do il=1,nd
     zzz=HACApK_entry_ij(il,il,st_bemv)
     ao(il)=1.0d0/zzz
   enddo
   
   call MPI_Barrier( icomm, ierr )
   st_measure_time_bicgstab=MPI_Wtime()
   if(param(85)==1)then
!     call HACApK_bicgstab_lfmtx(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
     call HACApK_bicgstab_lfmtx_hyp(st_leafmtxp,st_ctl,u,b,param,nd,nstp,lrtrn)
   elseif(param(85)==2)then
     call HACApK_gcrm_lfmtx(st_leafmtxp,st_ctl,st_bemv,u,b,param,nd,nstp,lrtrn)
   else
   endif
   call MPI_Barrier( icomm, ierr )
   en_measure_time_bicgstab=MPI_Wtime()
   time_bicgstab = en_measure_time_bicgstab - st_measure_time_bicgstab
   if(st_ctl%param(1)>0 .and. mpinr==0)  write(6,2000)              'time_HACApK_solve  =',time_bicgstab
   if(st_ctl%param(1)>0 .and. mpinr==0 .and. nstp>1)  write(6,2000) '       time_1step  =',time_bicgstab/nstp
   allocate(www(nd))
   sol(:nd)=0.0d0; www(lod(:nd))=u(:nd); sol(:nd)=www(:nd)
   deallocate(www)
   if(param(61)==3)then
     do il=1,nd
       sol(il)=sol(il)*st_bemv%ao(il)
     enddo
   endif
 endif
9999 continue
 HACApK_solve=lrtrn
 endfunction

endmodule m_HACApK_use
