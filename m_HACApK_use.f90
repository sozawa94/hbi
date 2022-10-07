!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 4.0.0                                             *
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
!C  This file includes example routines for Lattice H-matrices
!C  created by Akihiro Ida on March 2020
!C  added functions related to thread tarallelized HMVM to HACApK3.0.0 on May 2022
!C  last modified by Akihiro Ida on May 2022
!C**************************************************************************
module m_HACApK_use
 use m_HACApK_solve
 use m_HACApK_base
 implicit real*8(a-h,o-z)
contains

!*** HACApK_gensolv
 integer function     HACApK_gensolv(st_leafmtxp,st_bemv,st_ctl,gmid,rhs,sol,ztol)
   implicit none
   include 'mpif.h'
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_leafmtxp),pointer :: st_lfpp(:,:)
   type(st_HACApK_LHp) :: st_LHp
   type(st_HACApK_calc_entry) :: st_bemv
   type(st_HACApK_lcontrol) :: st_ctl
   type(st_HACApK_lcontrol),pointer :: st_ctlp(:,:)
   type(st_HACApK_latticevec),target :: st_au,st_u,wws
   real*8 :: gmid(st_bemv%nd,3),rhs(st_bemv%nd),sol(st_bemv%nd),ztol
   real*8,dimension(:),allocatable :: www
   integer*4 :: mpinr,mpilog,nrank,icomm,nthr,mpinrt,mpinrl,icommt,icomml,lrtrn,i,ierr
   integer*4 :: nlfl,nlft,itb,ilb,nbndst,nbndsl,ndlfs,ndtfs
   real*8 :: ans,ans0,sign,etime,etime_all,etimec,etimer,etimeb
!   real*8,dimension(:),pointer :: u,au,tmp
   type(st_HACApK_latticevec),pointer :: u,au,tmp
 integer*4,pointer :: lpmd(:)
 real*8,pointer :: time(:)
1000 format(5(a,i10)/)
2000 format(5(a,1pe15.8)/)
3000 format(5(a,f13.6)/)

 
   lrtrn=HACApK_generate(st_leafmtxp,st_bemv,st_ctl,gmid,ztol)
   lrtrn=HACApK_construct_LH(st_LHp,st_leafmtxp,st_bemv,st_ctl,gmid,ztol)
 
 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpilog=lpmd(4); nrank=lpmd(2); mpinrl=lpmd(37); mpinrt=lpmd(33)
 icommt=lpmd(31); icomml=lpmd(35); icomm=lpmd(1)

!   print *, "finish construct"
!   allocate(wws(st_leafmtxp%ndlfs+8))
   lrtrn=HACApK_gen_lattice_vectort(st_u,st_leafmtxp,st_ctl)
   lrtrn=HACApK_gen_lattice_vectort(st_au,st_leafmtxp,st_ctl)
   lrtrn=HACApK_gen_lattice_vectorl(wws,st_leafmtxp,st_ctl)
  allocate(www(st_leafmtxp%ndlfs))

!   print *, "finish gen_lattice_vec"
   mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
   mpinrt=st_ctl%lpmd(33); mpinrl=st_ctl%lpmd(37) 
   icommt=st_ctl%lpmd(31); icomml=st_ctl%lpmd(35)
 ndlfs=st_LHp%ndlfs; ndtfs=st_LHp%ndtfs

   ! sign = 1.0d0
   ! do i = 1, size(st_u%vs)
   !    st_u%vs(i) = sign * 1.0d0/i/10000
   !    st_au%vs(i) = sign * 1.0d0/i/10000
   !    ! st_u%vs(i) = sign * 1.0d0
   !    ! st_au%vs(i) = sign * 1.0d0
   !    sign = sign * -1.0d0
   ! end do
   st_u%vs(:) = 1.0
   st_au%vs(:) = 1.0
   st_ctl%time(:) = 0.0
   u => st_u
   au => st_au
   do i = 1, 1
      ans = 0.0d0
      ans0 = 0.0d0
      call HACApK_adot_lattice_hyp(au,st_LHp,st_ctl,www,u)
      ans = sum(au%vs(1:au%ndc))
!!!      print *, mpinr,mpinrt,ans,au%ndc
      if(0 .ne. mpinrl) then 
         ans = 0.0d0
      end if
      call MPI_reduce(ans,ans0,1,MPI_REAL8,MPI_SUM,0,icomm,ierr)
      if(mpinr == 0) print *, "reduce ans", i, ans0
      tmp => u
      u => au
      au => tmp
   end do

!   write(mpilog,*) '  '
!   write(mpilog,*) u%vs(1:10)
!   write(mpilog,*) '  '
!   write(mpilog,*) u%vs(511:521)

   print *, "etime_process", mpinr, st_ctl%time(1), st_ctl%time(2), st_ctl%time(3)
   etime = st_ctl%time(1)
   call MPI_reduce(etime,etimec,1,MPI_REAL8,MPI_SUM,0,icomm,ierr)
   etime = st_ctl%time(2)
   call MPI_reduce(etime,etimer,1,MPI_REAL8,MPI_SUM,0,icomm,ierr)
   etime = st_ctl%time(3)
   call MPI_reduce(etime,etimeb,1,MPI_REAL8,MPI_SUM,0,icomm,ierr)
   etime = st_ctl%time(1) + st_ctl%time(2) + st_ctl%time(3)
   call MPI_reduce(etime,etime_all,1,MPI_REAL8,MPI_MAX,0,icomm,ierr)
   if(mpinr == 0) write (*,'(A,7F12.6)') "time(all/itr,comp,reduce,bcast),hmat[GB],throughput[GB/s],GFlops", etime_all/30.0, etimec/nrank,etimer/nrank,etimeb/nrank,st_ctl%param(97), st_ctl%param(97)*30/etime_all, st_ctl%param(97)*30/4.0/etime_all

   HACApK_gensolv=lrtrn
   call MPI_Barrier( icomm, ierr )
   stop
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
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'irank=',mpinr,', nrank=',nrank
 nd=nofc*nffc
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nd=',nd,' nofc=',nofc,' nffc=',nffc
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'nrank=',nrank,' nth=',nthr
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'param:'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,10000) st_ctl%param(1:100)
 10000 format(10(1pe10.3)/)
 allocate(lnmtx(3))
 call MPI_Barrier( icomm, ierr )
 st_s=MPI_Wtime()

 !do il=1,nd
 !  write(42,*) coord(il,1:3)
 !enddo

 if(st_ctl%param(8)==10)then
   call HACApK_generate_frame_blrmtx(st_leafmtxp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)
 elseif(st_ctl%param(8)==20)then
   call HACApK_generate_frame_blrleaf(st_leafmtxp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)
 else
   call HACApK_generate_frame_leafmtx(st_leafmtxp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)
 endif
 
!!! stop
 
 if(st_leafmtxp%nlf<1)then
   print*,'ERROR!; sub HACApK_generate; irank=',mpinr,' nlf=',st_leafmtxp%nlf; goto 9999
 endif

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
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'ndnr_s=',st_ctl%lpmd(6),', ndnr_e=',st_ctl%lpmd(7),', ndnr=',st_ctl%lpmd(5)
 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,' ndlf_s=',st_ctl%lpmd(11),', ndlf_e=',st_ctl%lpmd(12),', ndlf=',st_leafmtxp%nlf
 lnps=nd+1; lnpe=0
 if(st_leafmtxp%nlf<1)then
   print*,'ERROR!; sub HACApK_generate; irank=',mpinr,' nlf=',st_leafmtxp%nlf
 endif
 
 if(st_ctl%param(7)==1)call HACApK_gen_mat_plot(st_leafmtxp,st_ctl%lpmd,st_ctl%lthr)
 
 if(st_ctl%param(10)==0) return
 call HACApK_fill_leafmtx_hyp(st_leafmtxp%st_lf,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe,st_ctl%lthr)
! call HACApK_fill_leafmtx(st_leafmtxp%st_lf,st_bemv,st_ctl%param,znrm,st_ctl%lpmd,lnmtx,st_ctl%lod,st_ctl%lod,nd,st_leafmtxp%nlf,lnps,lnpe)
 call MPI_Barrier( icomm, ierr )
 if(st_ctl%param(1)>1)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:3)
 ndnr_s=st_ctl%lpmd(6); ndnr_e=st_ctl%lpmd(7); ndnr=st_ctl%lpmd(5)

 st_fill_hmtx=MPI_Wtime()
 if(st_ctl%param(1)>1)  write(mpilog,2000)  'time_supermatrix             =',st_create_hmtx- st_s
 if(st_ctl%param(1)>1)  write(mpilog,2000)  'time_fill_hmtx               =',st_fill_hmtx-st_cal_matnorm
 if(st_ctl%param(1)>1)  write(mpilog,2000)  'time_construction_Hmatrix    =',st_fill_hmtx-st_s

 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_supermatrix             =',st_create_hmtx - st_s
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_fill_hmtx               =',st_fill_hmtx - st_cal_matnorm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'time_construction_Hmatrix    =',st_fill_hmtx - st_s

 call MPI_Barrier( icomm, ierr )

 call HACApK_chk_leafmtx(st_leafmtxp,st_ctl,lnmtx,nd,mem8)

 ktp=0
 call HACApK_setcutthread(st_ctl%lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
      
 call MPI_Barrier( icomm, ierr )
 if(st_ctl%param(8)==10 .or. st_ctl%param(8)==20)then
 else
!   print*,'mpinr=',mpinr,lnps,lnpe
   st_ctl%lnp(mpinr+1)=lnpe-lnps
   call MPI_Barrier( icomm, ierr )
   call MPI_Allgather(lnpe-lnps,1,MPI_INTEGER,st_ctl%lnp,1, MPI_INTEGER, icomm, ierr )
   st_ctl%lsp(mpinr+1)=lnps
   call MPI_Allgather(lnps,1,MPI_INTEGER,st_ctl%lsp,1, MPI_INTEGER, icomm, ierr )

   if(st_ctl%param(1)>1 .and. mpinr==0) write(*,*) 'lnp=',st_ctl%lnp(:)
   if(st_ctl%param(1)>1 .and. mpinr==0) write(*,*) 'lsp=',st_ctl%lsp(:)
 endif
 
 if(st_ctl%param(11)/=0) then
   call MPI_Barrier( icomm, ierr )
   call HACApK_accuracy_leafmtx(st_leafmtxp,st_bemv,st_ctl,st_ctl%lod,st_ctl%lod,st_ctl%lpmd,nofc,nffc)
 endif
9999 continue
 HACApK_generate=lrtrn
 endfunction
 
!*** HACApK_construct_LH
 integer function HACApK_construct_LH(st_LHp,st_leafmtxp,st_bemv,st_ctl,coord,ztol)
 include 'mpif.h'
 type(st_HACApK_LHp) :: st_LHp
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_lf9lttc),pointer :: st_lfp
 real*8 :: coord(st_bemv%nd,*)
 integer*8 :: mem8,nth1_mem,imem
 integer*4 :: ierr
 integer*4,dimension(:),allocatable :: lnmtx(:)
 integer*4,pointer :: lpmd(:)
 real*8,pointer :: time(:)
 1000 format(5(a,i10)/)
 2000 format(5(a,1pe15.8)/)
 
 lrtrn=0
 nofc=st_bemv%nd; nffc=1; ndim=3
 mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); nthr=st_ctl%lpmd(20)
 lpmd => st_ctl%lpmd(:); time => st_ctl%time(:)
 mpilog=lpmd(4); nrank=lpmd(2); mpinrl=lpmd(37); mpinrt=lpmd(33)
 icommt=lpmd(31); icomml=lpmd(35); icomm=lpmd(1)
 st_ctl%param(71)=ztol
 allocate(lnmtx(3))
 
!#!$omp master
! print*,'HACApK_construct_LH; chck point 01'
!#!$omp end master
 
 nd=nofc*nffc
 st_s=MPI_Wtime()
 if(st_ctl%param(8)==20)then
   call HACApK_generate_frame_LH(st_LHp,st_bemv,st_ctl,coord,lnmtx,nofc,nffc,ndim)
 else
   print*,'ERROR!; sub HACApK_generate_LH; irank=',mpinr,' st_ctl%param(8)=',st_ctl%param(8); goto 9999
 endif
 
 nlf=st_leafmtxp%nlf
 do is=1,nlf
   ilf=st_leafmtxp%st_lf(is)%lttcl
   itf=st_leafmtxp%st_lf(is)%lttct
   st_lfp=>st_LHp%st_lfp(ilf,itf)
   ip=st_lfp%nlf+1
   st_lfp%nlf=ip
   st_lfp%st_lf(ip)=st_leafmtxp%st_lf(is)
   st_lfp%st_lf(ip)%nstrtl=st_lfp%st_lf(ip)%nstrtl-st_LHp%lbstrtl(ilf)+1
   st_lfp%st_lf(ip)%nstrtt=st_lfp%st_lf(ip)%nstrtt-st_LHp%lbstrtt(itf)+1
 enddo
 
 if(st_ctl%param(1)>1) then
   write(mpilog,1000) 'HACApK_construct_LH; nlf=',st_LHp%nlf
   write(mpilog,1000) '       nlfl=',st_LHp%nlfl,'; nlft=',st_LHp%nlft
   write(mpilog,1000) '       ndlfs=',st_LHp%ndlfs,'; ndtfs=',st_LHp%ndtfs
   write(mpilog,*) 'st_LHp%lbstrtl='
     write(mpilog,'(10(i9))') st_LHp%lbstrtl(:)
   write(mpilog,*) 'st_LHp%lbstrtt='
     write(mpilog,'(10(i9))') st_LHp%lbstrtt(:)
   write(mpilog,*) 'st_LHp%lbndlfs='
     write(mpilog,'(10(i9))') st_LHp%lbndlfs(:)
   write(mpilog,*) 'st_LHp%lbndtfs='
     write(mpilog,'(10(i9))') st_LHp%lbndtfs(:)
   write(mpilog,*) 'st_LHp%latticel='
     write(mpilog,'(10(i9))') st_LHp%latticel(:)
   write(mpilog,*) 'st_LHp%latticet='
     write(mpilog,'(10(i9))') st_LHp%latticet(:)
   write(mpilog,*) 'st_LHp%lbndcsl='
     write(mpilog,'(10(i9))') st_LHp%lbndcsl(:)
   write(mpilog,*) 'st_LHp%lbndcst='
     write(mpilog,'(10(i9))') st_LHp%lbndcst(:)
!   do il=1,st_LHp%nlfl
!     write(mpilog,*) 'st_lfp_ndl=',st_LHp%st_lfp(il,:)%ndl
!     write(mpilog,*) 'st_lfp_ndt=',st_LHp%st_lfp(il,:)%ndt
!   enddo
 endif
 
 
 call MPI_Barrier( icomm, ierr )
        
!!!!$omp master
 nlfl=st_LHp%nlfl; nlft=st_LHp%nlft
 allocate(st_LHp%st_lfpp(nlfl,nlft),st_LHp%st_ctlp(nlfl,nlft))
 call HACApK_set_st_ctlp(st_LHp%st_ctlp,st_ctl,nlft,nlfl)
 do itb=1,nlft
   do ilb=1,nlfl
     st_LHp%st_lfpp(ilb,itb)%ktmax=st_leafmtxp%ktmax
     st_LHp%st_lfpp(ilb,itb)%nlf=st_LHp%st_lfp(ilb,itb)%nlf
     st_LHp%st_lfpp(ilb,itb)%st_lf => st_LHp%st_lfp(ilb,itb)%st_lf
     call HACApK_setcutthread_lattice(st_LHp%st_ctlp(ilb,itb)%lthr,st_LHp%st_lfpp(ilb,itb))
!!!     call HACApK_first_touch(st_LHp%st_lfpp(ilb,itb),st_LHp%st_ctlp(ilb,itb))
     call HACApK_contiguous(st_LHp%st_lfpp(ilb,itb),st_LHp%st_ctlp(ilb,itb))
!!!     call HACApK_transform(st_LHp%st_lfpp(ilb,itb),st_LHp%st_ctlp(ilb,itb))
!     write(mpilog,*) 'HACApK_gensolv; mpinrl=',mpinrl,'; mpinrt=',mpinrt,'; ilb=',ilb,'; itb=',itb,'; nlf=',st_LHp%st_lfpp(ilb,itb)%nlf
!     write(mpilog,*) st_LHp%st_ctlp(ilb,itb)%lthr
   enddo
 enddo 
!!!!$omp end master

!#!$omp master
!# print*,'HACApK_construct_LH; chck point 20'
!#!$omp end master

9999 continue
 HACApK_construct_LH=lrtrn
 endfunction

!***HACApK_set_st_ctlp
 subroutine HACApK_set_st_ctlp(st_ctlp,st_ctl,nlft,nlfl)
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_lcontrol) :: st_ctlp(:,:)
 
!$omp parallel
  nthr = omp_get_num_threads()
!$omp end parallel
! print*,'HACApK_set_st_ctlp',nlft,nlfl
 do itb=1,nlft
   do ilb=1,nlfl
     allocate(st_ctlp(ilb,itb)%lthr(nthr+1),st_ctlp(ilb,itb)%lpmd(50),stat = ierr)
     st_ctlp(ilb,itb)%lthr(:)=0
     st_ctlp(ilb,itb)%lpmd(:)=st_ctl%lpmd(:)
   enddo
 enddo
 end subroutine

!***HACApK_setcutthread_lattice
 subroutine HACApK_setcutthread_lattice(lthr,st_leafmtxp)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 integer :: lthr(0:)
 integer*8 :: mem8,nth1_mem,imem

 nthr=1
!$omp parallel
  nthr = omp_get_num_threads()
!$omp end parallel

 nlf=st_leafmtxp%nlf
 mem8=0
 do ip=1,nlf
   ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     mem8=mem8+(ndt+ndl)*kt
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     mem8=mem8+ns
   endif
 enddo

 nth1_mem=mem8/nthr
 lthr(0)=1; lthr(1:)=nlf+1
 imem=0; ith=1
 do il=1,nlf
   ltmtx=st_leafmtxp%st_lf(il)%ltmtx
   ndl=st_leafmtxp%st_lf(il)%ndl; ndt=st_leafmtxp%st_lf(il)%ndt; ns=ndl*ndt
   if(ltmtx==1)then
     kt=st_leafmtxp%st_lf(il)%kt
     imem=imem+(ndl+ndt)*kt
   elseif(ltmtx==2)then
     imem=imem+ns
   endif
   if(imem>nth1_mem*ith)then
     lthr(ith)=il
     ith=ith+1
     if(ith==nthr) exit
   endif
 enddo
 
! print*,'HACApK_setcutthread_lattice; mem8=',mem8,nth1_mem*ith
! print*,'HACApK_setcutthread_lattice; nlf=',nlf,lthr(:)
 
 end subroutine HACApK_setcutthread_lattice

endmodule m_HACApK_use
