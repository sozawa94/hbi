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
!C  This file includes basic routines for H-matrices
!C  created by Akihiro Ida at Kyoto University on May 2012
!C  last modified by Akihiro Ida on Sep. 2014
!C**************************************************************************
module m_HACApK_base
 use m_HACApK_calc_entry_ij
 implicit real*8(a-h,o-z)
 implicit integer*4(i-n)
 integer omp_get_thread_num, omp_get_num_threads
 external omp_get_thread_num, omp_get_num_threads

!*** type :: st_HACApK_cluster
  type :: st_HACApK_cluster
    integer*4 ndim
    integer*4 nstrt,nsize,ndpth,nnson,nmbr
    integer*4 ndscd ! number of descendants
    real*8,pointer :: bmin(:)=>null() ! Bounding box
    real*8,pointer :: bmax(:)=>null()
    real*8 :: zwdth
    type(st_HACApK_cluster),pointer :: pc_sons(:)=>null()
  end type st_HACApK_cluster

!*** type :: st_HACApK_leafmtx
  type :: st_HACApK_leafmtx
    integer*4 ltmtx  ! kind of the matrix; 1:rk 2:full
    integer*4 kt
    integer*4 nstrtl,ndl;
    integer*4 nstrtt,ndt;
    real*8,pointer :: a1(:,:)=>null(),a2(:,:)=>null()
  end type st_HACApK_leafmtx

!*** type :: st_HACApK_leafmtxp
  type :: st_HACApK_leafmtxp
    integer*4 nd ! number of unknowns
    integer*4 nlf ! number of sub-matrices
    integer*4 nlfkt ! number of low-rank sub matrices
    integer*4 ktmax
    type(st_HACApK_leafmtx),pointer :: st_lf(:)=>null()
  end type st_HACApK_leafmtxp

!*** type :: st_HACApK_lcontrol
  type :: st_HACApK_lcontrol
    integer*4,pointer :: lod(:)=>null(),lsp(:)=>null(),lnp(:)=>null(),lthr(:)=>null(),lpmd(:)=>null()
    real*8,   pointer :: param(:)=>null()
    integer :: lf_umpi
  end type st_HACApK_lcontrol

 interface
    real*8 function HACApK_unrm_d(nd,za)
       implicit real*8(a-h,o-z)
       real*8 :: za(:)
    end function

    subroutine HACApK_adotsub_dsm(zr,zaa,zu,ndl,ndt,mdl)
      implicit real*8(a-h,o-z)
      real*8 :: zaa(:,:)
      real*8 :: zu(:),zr(:)
    end subroutine

    subroutine HACApK_maxabsvalloc_d(za,zz,il,nd)
      implicit real*8(a-h,o-z)
      real*8 :: za(:)
    endsubroutine

    subroutine HACApK_maxabsvallocm_d(za,zz,il,nd,lmask)
      implicit real*8(a-h,o-z)
      real*8 :: za(:)
      integer lmask(:)
    endsubroutine

    subroutine HACApK_minabsvalloc_d(za,zz,il,nd)
      implicit real*8(a-h,o-z)
      real*8 za(:)
    end subroutine
   integer function HACApK_med3(nl,nr,nlr2)
   endfunction
 end interface

contains
!***HACApK_init
integer function HACApK_init(nd,st_ctl,st_bemv,icomma)
 implicit real*8(a-h,o-z)
 include 'mpif.h'
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_lcontrol) :: st_ctl
 integer,optional :: icomma
 character*32 logfile
 allocate(st_ctl%param(100))
 st_ctl%param(1:100)=0.0
 st_ctl%param(1) =1;        ! Print : 0:Only Error 1:STD 2:Dubug
 st_ctl%param(9) =1;        ! 1:load balancer
 st_ctl%param(10)=1;        ! 1:fulfill the matrix  0: not fulfill
 st_ctl%param(11)=0;        ! 1:check accuracy of H-matrix 0: not check
 st_ctl%param(12)=0;        ! 1:write down the H-matrix to file 0: not write
 st_ctl%param(21)=15;       ! cluster : leaf size 15
 st_ctl%param(22)=1.0;      ! cluster : max leaf size 1.0*nffc
 st_ctl%param(51)=2.0;      ! H-matrix : dicision param of distance 2.0
 st_ctl%param(60)=1         ! 1:ACA  2:ACA+
 st_ctl%param(61)=1         ! ACA norm 1:MREM  2:test 3:norm
 st_ctl%param(62)=7         ! ACA : predictive average of k
 st_ctl%param(63)=1000;     ! ACA : k-max of R_k-matrix 30
 st_ctl%param(64)=1;        ! ACA : minimun kt
 st_ctl%param(72)=1.0e-15;   ! ACA_EPS
 st_ctl%param(83)=500;      ! solver : maximum iterative number
 st_ctl%param(85)=1;        ! solver : 1:BiCGSTAB 2:GCR(m)
 st_ctl%param(87)=8;        ! solver : number of iteration for reset
 st_ctl%param(99)=100       ! Measure the time of Ax; iterative number
 ierr=0; lrtrn=0
 if(present(icomma))then
   icomm=icomma; st_ctl%lf_umpi=1
 else
   icomm=MPI_COMM_WORLD; lf_umpi=0
   call MPI_Init ( ierr )
   if( ierr .ne. 0 )  print*, 'HACApK_init; Error: MPI_Init failed !!!'
 endif
 call MPI_Comm_size ( icomm, nrank, ierr )
 if(ierr.ne.0) then
    print*, 'Error: MPI_Comm_size failed !!!'
 endif
 call MPI_Comm_rank ( icomm, irank, ierr )
 if(ierr.ne.0) then
    print*, 'Error: MPI_Comm_rank failed !!!'
 endif
 allocate(st_ctl%lpmd(30)); st_ctl%lpmd(:)=0
 st_ctl%lpmd(1)=icomm; st_ctl%lpmd(2)=nrank; st_ctl%lpmd(3)=irank; st_ctl%lpmd(4)=20
 nthr=1
!$omp parallel
  nthr = omp_get_num_threads()
  if(nthr>0) st_ctl%lpmd(20)=nthr
!$omp end parallel
 allocate(st_ctl%lod(nd),st_ctl%lthr(nthr),st_ctl%lnp(nrank),st_ctl%lsp(nrank),stat = ierr)

 call MPI_Barrier( icomm, ierr )

 if(st_ctl%param(1)>0)then
   write(logfile,'(a,i4.4,a)') 'log',irank,'.txt'
   open(st_ctl%lpmd(4),file=logfile)
 endif

 st_bemv%nd=nd
 if(st_bemv%lp61==1) then; endif
 if(ierr/=0)then; goto 9999; endif

9999 continue
 HACApK_init=lrtrn
 endfunction

!***HACApK_finalize
integer function HACApK_finalize(st_ctl)
 implicit real*8(a-h,o-z)
 include 'mpif.h'
 type(st_HACApK_lcontrol) :: st_ctl
 if(st_ctl%param(1)>0) close(st_ctl%lpmd(4))
 if(st_ctl%lf_umpi==0) call MPI_Finalize (ierr)
 lrtrn=HACApK_free_lcontrol(st_ctl)
 HACApK_finalize=lrtrn
endfunction

!***HACApK_chk_st_ctl
 subroutine HACApK_chk_st_ctl(st_ctl)
 type(st_HACApK_lcontrol) :: st_ctl
 if(st_ctl%param(21)<2) then
   st_ctl%param(21)=3
   print*,'Warning: Invalid param(21)=',st_ctl%param(21)
   print*,'It is changed to param(21)=3.'
 endif

 end subroutine HACApK_chk_st_ctl

!***HACApK_generate_frame_leafmtx
 subroutine HACApK_generate_frame_leafmtx(st_leafmtxp,st_bemv,st_ctl,gmid,lnmtx,nofc,nffc,ndim)
 include 'mpif.h'
 type(st_HACApK_cluster) :: st_clt
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_leafmtx),dimension(:), allocatable :: st_leafmtx
 real*8 :: gmid(nofc,ndim)
 integer*8 :: mem8
 integer*4 :: lnmtx(3)
 integer, dimension(:), allocatable :: lhp, lnp
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:),lod(:),lthr(:),lodfc(:)
 1000 format(5(a,i12)/)
 2000 format(5(a,e15.7)/)

 param => st_ctl%param(:)
 lpmd => st_ctl%lpmd(:); lod => st_ctl%lod(:); lthr(0:) => st_ctl%lthr
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1); nthr=lpmd(20)
 nd=nofc*nffc
 allocate(lodfc(nofc),stat=ierr)
 if(ierr/=0)then
   print*,'HACApK_generate_frame_leafmtx: ierr=',ierr
 endif
 do il=1,nofc
   lodfc(il)=il
 enddo
!!!!!!!!!!!!!!!!!! start clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nsrt=1; ndf=nofc; nclst=0; ndpth=0; ndscd=0
 call HACApK_generate_cbitree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst)
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'No. of cluster=',nclst
 if(st_ctl%param(1)>0)  write(mpilog,1000) 'No. of cluster=',nclst

 call HACApK_bndbox(st_clt,gmid,lodfc,nofc)
!!!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ndpth=0; lnmtx(1:3)=0
 call HACApK_count_lntmx(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'No. of nsmtx',lnmtx(1:3)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'   1:Rk-matrix 2: dense-mat 3:H-matrix'
 st_leafmtxp%nlfkt=lnmtx(1)
 nlf=lnmtx(1)+lnmtx(2)
 allocate(st_leafmtx(nlf))
 st_leafmtxp%nlf=nlf; nlf=0

 call HACApK_generate_leafmtx(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_generate_leafmtx end'
 call HACApK_qsort_row_leafmtx(st_leafmtx,1,nlf)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_qsort_row_leafmtx end'
 ilp=1; ips=1
 do ip=1,nlf
   il=st_leafmtx(ip)%nstrtl
   if(il<ilp)then    ; print *,'Error!; HACApK_generate_frame_leafmtx row_sort';
   elseif(il>ilp)then;
     call HACApK_qsort_col_leafmtx(st_leafmtx,ips,ip-1)
     ilp=il; ips=ip
   endif
 enddo

 call HACApK_free_st_clt(st_clt)

 do il=1,nofc
   do ig=1,nffc
     is=ig+(il-1)*nffc
     lod(is)=lodfc(il)
   enddo
 enddo

!!!!!!!!!!!!!!!!!! start MPI load balance  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(lhp(nrank+1), lnp(nrank), stat = ierr)
  if(ierr.ne.0) then
    print*, 'Memory allocation #11 failed !!!'
    goto 9999
  endif

  mdim  = nd / nrank
  mdim1 = mdim + 1
  mrank = nd - mdim * nrank

  lbrns=param(9)
  if(nrank==2) lbrns=0
  if(lbrns==0)then
    do jrank=1,mrank; lnp(jrank) = mdim1; enddo
    do jrank=mrank+1,nrank; lnp(jrank) = mdim; enddo

    lhp(1)=0; lhp(nrank+1)=nd
    do jrank=2,nrank
      lhp(jrank)=lhp(jrank-1)+lnp(jrank-1)
    enddo
  else
    ktp=param(62)
    call HACApK_setcutrow(lhp,mem8,lpmd,st_leafmtx,st_ctl,nd,st_leafmtxp%nlf,ktp)
    do jrank=1,nrank
      lnp(jrank)=lhp(jrank+1)-lhp(jrank)
    enddo
  endif

  lpmd(6)=lhp(mpinr+1)+1; lpmd(7)=lhp(mpinr+2); lpmd(5)=lpmd(7)-lpmd(6)+1
   call MPI_Barrier( icomm, ierr )

 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'nlf=',st_leafmtxp%nlf
! do il=1,mmm
!   if(mpinr==0) write(*,1000) 'il=',il,' nstrtl=',st_leafmtx(il)%nstrtl,' nstrtt=',st_leafmtx(il)%nstrtt
! enddo
! stop

  lpmd(11)=0; lpmd(12)=st_leafmtxp%nlf
  do ip=1,st_leafmtxp%nlf
    if(st_leafmtx(ip)%nstrtl>=lpmd(6) .and. st_leafmtx(ip)%nstrtl<=lpmd(7) .and. lpmd(11)==0) lpmd(11)=ip
    if(st_leafmtx(ip)%nstrtl>lpmd(7) .and. lpmd(11)/=0) then
      lpmd(12)=ip-1; exit
    endif
  enddo
!  write(*,1000) 'irank=',mpinr,' ndlf_s=',lpmd(11),', ndlf_e=',lpmd(12)
  st_leafmtxp%nlf=lpmd(12)-lpmd(11)+1
  allocate(st_leafmtxp%st_lf(st_leafmtxp%nlf)); st_leafmtxp%st_lf(1:st_leafmtxp%nlf)=st_leafmtx(lpmd(11):lpmd(12))
  lnmtx(:)=0; mem8=0; ktp=param(62)
  do ip=lpmd(11),lpmd(12)
    ltmtx=st_leafmtx(ip)%ltmtx; ndl=st_leafmtx(ip)%ndl; ndt=st_leafmtx(ip)%ndt; ns=ndl*ndt
    if(ltmtx==1)then
      lnmtx(1)=lnmtx(1)+1; mem8=mem8+(ndt+ndl)*ktp
    else
      lnmtx(2)=lnmtx(2)+1; mem8=mem8+ns
    endif
  enddo
  deallocate(st_leafmtx)
  call HACApK_setcutthread(lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
 9999 continue
! stop
 end subroutine HACApK_generate_frame_leafmtx

!***HACApK_setcutrow
 subroutine HACApK_setcutrow(lhp,ncpc1,lpmd,st_leafmtx,st_ctl,nd,nlf,kt)
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 type(st_HACApK_lcontrol) :: st_ctl
 integer :: lhp(*),lpmd(*)
 integer*8 :: ncpc,ncpc1
 1000 format(5(a,i12)/)
 2000 format(5(a,e15.7)/)

 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)

 call MPI_Barrier( icomm, ierr )
 ncpc=0
 do il=1,nlf
   ltmtx=st_leafmtx(il)%ltmtx
   ndl=st_leafmtx(il)%ndl; ndt=st_leafmtx(il)%ndt
   if(ltmtx==1)then
     ncpc=ncpc+(ndl+ndt)*kt
   else
     ncpc=ncpc+ndl*ndt
   endif
 enddo
 ncpc1=ncpc/nrank

 if(st_ctl%param(1)>0 .and. mpinr==0) write(6,1000) 'ncpc=',ncpc,'  ncpc/nrank=',ncpc1

 lhp(1)=0; lhp(nrank+1)=nd
 ncpc=0; ip=1
 do il=1,nlf
   ltmtx=st_leafmtx(il)%ltmtx; nstrtl=st_leafmtx(il)%nstrtl
   ndl=st_leafmtx(il)%ndl; ndt=st_leafmtx(il)%ndt
   if(ltmtx==1)then
     ncpc=ncpc+(ndl+ndt)*kt
   else
     ncpc=ncpc+ndl*ndt
   endif
   if(ncpc>ncpc1*ip .and. nstrtl-1>lhp(ip))then
     lhp(ip+1)=nstrtl-1
!     if(mpinr==0) write(6,1000) 'ip=',ip,' lhp(ip+1)=',lhp(ip+1)
     ip=ip+1
     if(ip==nrank) exit
   endif
 enddo
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_setcutrow end'
 end subroutine HACApK_setcutrow

!***HACApK_setcutthread
 subroutine HACApK_setcutthread(lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 integer :: lthr(0:*)
 integer*8 :: mem8,nth1_mem,imem

 nlf=st_leafmtxp%nlf
 nth1_mem=mem8/nthr
 if(st_ctl%param(1)>1) print*,'HACApK_setcutthread; nlf=',nlf,' mem8=',mem8,' nthr=',nthr
 lthr(0)=1; lthr(nthr)=nlf+1
 imem=0; ith=1; kt=ktp
 do il=1,nlf
   ltmtx=st_leafmtxp%st_lf(il)%ltmtx
   ndl=st_leafmtxp%st_lf(il)%ndl; ndt=st_leafmtxp%st_lf(il)%ndt
   !print*,'HACApK_setcutthread; ndl=',ndl,' ndt=',ndt
   if(ltmtx==1)then
     if(ktp==0) kt=st_leafmtxp%st_lf(il)%kt
     imem=imem+(ndl+ndt)*kt
   else
     imem=imem+ndl*ndt
   endif
   if(imem>nth1_mem*ith)then
     lthr(ith)=il
     ith=ith+1
     if(ith==nthr) exit
   endif
 enddo
 if(st_ctl%param(1)>1) print*,'HACApK_setcutthread; lthr=',lthr(0:nthr)
 end subroutine HACApK_setcutthread

!***HACApK_accuracy_leafmtx_body
 subroutine HACApK_accuracy_leafmtx_body(zhnrm,zanrm,st_leafmtxp,st_bemv,lodl,lodt,lpmd,nofc,nffc)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_calc_entry) :: st_bemv
 integer*4 :: lodl(nofc*nffc),lodt(nofc*nffc),lpmd(*)
 1000 format(5(a,i12)/)
 2000 format(5(a,1pe15.7)/)
 do ip=1,st_leafmtxp%nlf
   ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         zz=HACApK_entry_ij(lodl(ill),lodt(itt),st_bemv)
         zanrm=zanrm+zz*zz
         do ik=1,kt
           zz=zz-st_leafmtxp%st_lf(ip)%a2(il,ik)*st_leafmtxp%st_lf(ip)%a1(it,ik)
         enddo
         zhnrm=zhnrm+zz*zz
       enddo
     enddo
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         zanrm=zanrm+st_leafmtxp%st_lf(ip)%a1(it,il)*st_leafmtxp%st_lf(ip)%a1(it,il)
       enddo
     enddo
   endif
 enddo
 end subroutine HACApK_accuracy_leafmtx_body

!***HACApK_accuracy_leafmtx
 subroutine HACApK_accuracy_leafmtx(st_leafmtxp,st_bemv,st_ctl,lodl,lodt,lpmd,nofc,nffc)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_lcontrol) :: st_ctl
 integer*4 :: lodl(nofc*nffc),lodt(nofc*nffc),lpmd(*)
 1000 format(5(a,i12)/)
 2000 format(5(a,1pe15.8)/)

 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx start'
 zhnrm=0.0d0; zanrm=0.0d0; ktmax=0
 call HACApK_accuracy_leafmtx_body(zhnrm1,zanrm1,st_leafmtxp,st_bemv,lodl,lodt,lpmd,nofc,nffc)
 call MPI_reduce( zanrm1, zanrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 call MPI_reduce( zhnrm1, zhnrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 zanrm=dsqrt(zanrm); zhnrm=dsqrt(zhnrm)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx; zanrm=',zanrm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx; zhnrm=',zhnrm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx; zhnrm/zanrm=',zhnrm/zanrm
 if(st_ctl%param(1)>0)  write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zanrm1=',dsqrt(zanrm1)
 if(st_ctl%param(1)>0)  write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm1=',dsqrt(zhnrm1)
 if(st_ctl%param(1)>0)  write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm1/zanrm1=',dsqrt(zhnrm1)/dsqrt(zanrm1)
  if(mpinr==0 .and. st_ctl%param(1)>0) write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zanrm=',zanrm
  if(mpinr==0 .and. st_ctl%param(1)>0) write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm=',zhnrm
  if(mpinr==0 .and. st_ctl%param(1)>0) write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm/zanrm=',zhnrm/zanrm

 end subroutine HACApK_accuracy_leafmtx

!***HACApK_calc_vec
! ld==0: row direction, ld==1: column direction
 subroutine HACApK_calc_vec(zaa, zab, ndp, k, ip, vec, nstrtl, nstrtt,lod, st_bemv, lmsk, ld)
 type(st_HACApK_calc_entry) :: st_bemv
 real*8,target :: zab(:,:)
 real*8 :: vec(:),zaa(:,:)
 real*8,pointer :: zz(:)
 integer*4 :: lmsk(:),lod(:)

 do ii=1,ndp
   if(lmsk(ii)==0) then
     if(ld==0)then
       ill=ip+nstrtl-1; itt=ii+nstrtt-1
     else
       ill=ii+nstrtl-1; itt=ip+nstrtt-1
     endif
     vec(ii)=HACApK_entry_ij(lod(ill),lod(itt),st_bemv)
   endif
 enddo
 if(k==0) return
 zz => zab(ip,1:k)
 call HACApK_adotsub_dsm(vec,zaa,zz,ndp,k,ndp)
 where(lmsk==1) vec=0.0d0
 endsubroutine

!***HACApK_aca
 integer function HACApK_aca(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lod,st_bemv,kmax,eps,znrmmat,pACA_EPS)
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: param(:)
 real*8,target :: zaa(ndl,kmax),zab(ndt,kmax)
 real*8,pointer :: prow(:),pcol(:)
 integer*4 :: lod(:)
 integer*4,dimension(:), allocatable :: lrow_msk,lcol_msk
 1000 format(5(a,i12)/)
 2000 format(5(a,1pe15.8)/)

! write(6,1000) 'nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt, 'kmax=',kmax
 krank=min(ndl,ndt)
 znrm=znrmmat*sqrt(real(ndl)*real(ndt))
 if(param(61)==1) ACA_EPS=pACA_EPS
 if(param(61)==2 .or. param(61)==3) ACA_EPS=pACA_EPS*znrm
 allocate(lrow_msk(ndl),lcol_msk(ndt)); lrow_msk(:)=0; lcol_msk(:)=0; nrow_done=0; ncol_done=0
 HACApK_aca=0; k=1; lstop_aca=0
 kstop=min(kmax,krank)
 if(nstrtl>nstrtt)then; ist=1
 else; ist=ndl
 endif
 do
  !print*,'k=',k,' zeps=',zeps,' ist=',ist
   pcol => zaa(:,k); prow => zab(:,k)
    call HACApK_calc_vec(zab, zaa, ndt, k-1, ist, prow, nstrtl, nstrtt,lod, st_bemv, lcol_msk,0)
   !print*, 'prow=',prow
   call HACApK_maxabsvallocm_d(prow,row_maxval,jst,ndt,lcol_msk)
   !print*,'jst=',jst,' row_maxval=',row_maxval
 !write(6,1000) 'ist=',ist,'  jst=',jst
     zdltinv=1.0d0/prow(jst); prow(:)=prow(:)*zdltinv
    call HACApK_calc_vec(zaa, zab, ndl, k-1, jst, pcol, nstrtl, nstrtt,lod, st_bemv, lrow_msk,1)
   lrow_msk(ist)=1
   call HACApK_maxabsvallocm_d(pcol,col_maxval,istn,ndl,lrow_msk)
   !print*, 'pcol=',pcol
   !print *,'zeps=',zeps
   !print*,'istn=',istn,' col_maxval=',col_maxval
   lcol_msk(jst)=1
   ist=istn; nrow_done=nrow_done+1; ncol_done=ncol_done+1
   !print *, 'abs(row_maxval)=',abs(row_maxval)
   !print *, 'abs(col_maxval)=',abs(col_maxval)
   if(abs(row_maxval)<ACA_EPS .and. abs(col_maxval)<ACA_EPS .and. k>=param(64)) then
!$omp critical
     print *, 'ACA_EPS=',ACA_EPS
     print *, 'abs(row_maxval)=',abs(row_maxval)
     print *, 'abs(col_maxval)=',abs(col_maxval)
     print *, 'stop HACApK_aca 3';
!!!     stop
!$omp end critical
     goto 9999
   endif
   zeps=HACApK_unrm_d(ndl,pcol)*HACApK_unrm_d(ndt,prow)
!   zcolm=HACApK_unrm_d(ndl,pcol); zrowm=HACApK_unrm_d(ndt,prow)
!   zeps=max(zcolm,zrowm,zcolm*zrowm)
   if(k==1 .and. param(61)==1) znrm=zeps
   zeps=zeps/znrm
!  print*,'pcol'; print*,pcol
!  print*,'prow'; print*,prow
!  print*,'lcol_msk',lcol_msk
!  print*,'lrow_msk',lrow_msk
!  write(6,2000) 'zeps=',zeps
  if(zeps<eps .or. k==kstop) lstop_aca = 1
  if(lstop_aca==1 .and. k>=param(64)) then
!    print *,'k=',k, 'param(64)=',param(64)
!    print *,'zeps=',zeps
!    print *,'eps=',eps
!    print*,'???????????'
    exit
  endif
  k=k+1
 enddo
 9999 continue
 deallocate(lrow_msk,lcol_msk)
 HACApK_aca=k
! print*,'HACApK_aca=',aca
  if(zeps>eps .and. k<krank)then
!$omp critical
    print *,'k=',k
    print *,'zeps=',zeps
    print *,'eps=',eps
    write(6,1000) 'nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt
    print*,'znrm=',znrm
    stop
!$omp end critical
  endif
! stop
 endfunction

!***HACApK_fill_leafmtx_hyp
 subroutine HACApK_fill_leafmtx_hyp(st_lf,st_bemv,param,znrmmat,lpmd,lnmtx,lodl,lodt,nd,nlf,lnps,lnpe,lthr)
! type(st_HACApK_leafmtxp) ::  st_leafmtxp
 type(st_HACApK_leafmtx) :: st_lf(:)
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 ::param(:)
 integer*4 :: lodl(nd),lodt(nd),lpmd(:),lnmtx(:),lthr(0:)
 real*8, allocatable :: zab(:,:),zaa(:,:)
 1000 format(5(a,i12)/)
 eps=param(71); ACA_EPS=param(72)*eps; kparam=param(63)
!$OMP parallel default(none) &
!$OMP          shared(st_lf,lodl,st_bemv,znrmmat,lodt,lthr,param) &
!$OMP          private(zab,zaa,kt,ith,ith1,nths,nthr,nthe,ltmtx,ierr,ndl,ndt,ns,nstrtl,nstrtt,ip,il,it,ill,itt) &
!$OMP          firstprivate(eps, ACA_EPS, kparam)
 ith = omp_get_thread_num()
 nthr = omp_get_num_threads()
 if(nthr == 0) write(* ,*) nthr ,ith
 !$OMP barrier
 ith1 = ith+1
 nths=lthr(ith); nthe=lthr(ith1)-1
 ierr=0
 do ip=nths,nthe
   ndl   =st_lf(ip)%ndl   ; ndt   =st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_lf(ip)%nstrtl; nstrtt=st_lf(ip)%nstrtt; ltmtx=st_lf(ip)%ltmtx
   if(ltmtx==1)then
     allocate(zab(ndt,kparam),zaa(ndl,kparam),stat=ierr)
     if(ierr.ne.0) then
!$omp critical
        write(*,*) 'sub HACApK_fill_leafmtx_p; zab,zaa Memory allocation failed !'
        write(*,1000) 'ip=ip',ip,' ierr=',ierr
!$omp end critical
        stop 10
     endif
     if(param(60)==1)then
       kt=HACApK_aca(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
     else
       print*,'Only ACA is avairable! Set param(60)=1.'
       stop
     endif
     if(kt>kparam-1) then
!!!$omp critical
!        write(*,1000) 'WARNING: Insufficient k: kt=',kt,', kparam=',kparam, &
!                      ' nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt
!!!$omp end critical
     endif
     st_lf(ip)%kt=kt
     allocate(st_lf(ip)%a1(ndt,kt),st_lf(ip)%a2(ndl,kt),stat=ierr)
     if(ierr.ne.0) then
!$omp critical
        write(*,*) 'sub HACApK_fill_leafmtx_p; a1,a2 Memory allocation failed !'
        write(*,*) 'ip=ip',ip,' ierr=',ierr
!$omp end critical
        stop 20
     endif
     st_lf(ip)%a1(1:ndt,1:kt)=zab(1:ndt,1:kt)
     st_lf(ip)%a2(1:ndl,1:kt)=zaa(1:ndl,1:kt)
     deallocate(zab,zaa)
   elseif(ltmtx==2)then
     allocate(st_lf(ip)%a1(ndt,ndl),stat=ierr)
     if(ierr.ne.0) then
!$omp critical
        write(*,*) 'sub HACApK_fill_leafmtx_p; a2 Memory allocation failed !'
        write(*,*) 'ip=ip',ip,' ierr=',ierr
!$omp end critical
        stop 30
     endif
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         st_lf(ip)%a1(it,il)=HACApK_entry_ij(lodl(ill),lodt(itt),st_bemv)
       enddo
     enddo
   else
!$omp critical
      write(*,1000) 'HACApK_fill_leafmtx_hyp; ip=',ip,' ltmtx=',ltmtx
!$omp end critical
   endif
 enddo
!$omp end parallel
 do ip=1,nlf
   ndl=st_lf(ip)%ndl; nstrtl=st_lf(ip)%nstrtl
   if(nstrtl<lnps) lnps=nstrtl
   if(nstrtl+ndl>lnpe) lnpe=nstrtl+ndl
 enddo
 end subroutine HACApK_fill_leafmtx_hyp

!***HACApK_fill_leafmtx
 RECURSIVE subroutine HACApK_fill_leafmtx(st_lf,st_bemv,param,znrmmat,lpmd,lnmtx,lodl,lodt,nd,nlf,lnps,lnpe)
 type(st_HACApK_leafmtx) ::  st_lf(:)
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 ::param(:)
 integer*4 :: lodl(nd),lodt(nd),lpmd(:),lnmtx(:)
 real*8,pointer :: zab(:,:),zaa(:,:)
 integer*4,dimension(:),allocatable :: lodc
 1000 format(5(a,i12)/)
 eps=param(71); ACA_EPS=param(72)*eps; kparam=param(63)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5); mpinr=lpmd(3)

 do ip=1,nlf
   ndl   =st_lf(ip)%ndl   ; ndt   =st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_lf(ip)%nstrtl; nstrtt=st_lf(ip)%nstrtt
   if(nstrtl<lnps) lnps=nstrtl
   if(nstrtl+ndl>lnpe) lnpe=nstrtl+ndl
   if(st_lf(ip)%ltmtx==1)then
     allocate(zab(ndt,kparam),zaa(ndl,kparam),lodc(nd))
     lodc(:)=lodl(:)-1
     nstrtt_c=nstrtt-1; nstrtl_c=nstrtl-1
     if(param(60)==1)then
       kt=HACApK_aca(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
     else
       print*,'Only ACA is avairable! Set param(60)=1.'
       stop
     endif
     if(kt>kparam-1) then
        if(param(1)>1 .and. mpinr==0) write(*,1000) 'WARNING: HACApK_fill_leafmtx; Insufficient k: kt=',kt,', kparam=',kparam, &
                      ' nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt
!      stop
     endif
     st_lf(ip)%kt=kt
     allocate(st_lf(ip)%a1(ndt,kt),st_lf(ip)%a2(ndl,kt),stat=ierr)
     if(ierr.ne.0) then
        print*, 'sub HACApK_fill_leafmtx_p; Memory allocation failed !'
        stop
     endif
     st_lf(ip)%a1(1:ndt,1:kt)=zab(1:ndt,1:kt)
     st_lf(ip)%a2(1:ndl,1:kt)=zaa(1:ndl,1:kt)
     deallocate(zab,zaa,lodc)
   elseif(st_lf(ip)%ltmtx==2)then
     allocate(st_lf(ip)%a1(ndt,ndl),stat=ierr)
     do il=1,ndl; ill=il+nstrtl-1
       do it=1,ndt; itt=it+nstrtt-1
         st_lf(ip)%a1(it,il)=HACApK_entry_ij(lodl(ill),lodt(itt),st_bemv)
       enddo
     enddo
   endif
 enddo
 end subroutine HACApK_fill_leafmtx

!***HACApK_chk_leafmtx
 subroutine HACApK_chk_leafmtx(st_leafmtxp,st_ctl,lnmtx,nd,mem8)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 integer :: lnmtx(3)
 integer*8 :: mem8,memh
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10))
 2000 format(5(a,f12.2))

 param => st_ctl%param(:)
 lpmd => st_ctl%lpmd(:)
 mpinr=lpmd(3); mpilog=lpmd(4); icomm=lpmd(1)
 kparam=param(63)
 ktofl=0 ;ktoft=0 ;ktmax=0; mem8=0; ktav=0
 do ip=1,st_leafmtxp%nlf
   ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
     kt=st_leafmtxp%st_lf(ip)%kt
     if(kt>ktmax) ktmax=kt
     if((nstrtl+ndl)>nd .and. nstrtt==1) then
        ktofl=kt
        if(st_ctl%param(1)>1 .and. mpinr==0) write(*,1000) 'leftlower  value of kt=',ktofl,' ndl=',ndl,' ndt=',ndt
     endif
     if((nstrtt+ndt)>nd .and. nstrtl==1) then
        ktoft=kt
        if(st_ctl%param(1)>1 .and. mpinr==0) write(*,1000) 'rightupper value of kt=',ktoft,' ndl=',ndl,' ndt=',ndt
     endif
     mem8=mem8+(ndt+ndl)*kt
     ktav=ktav+kt
   elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
     mem8=mem8+ns
   endif
 enddo
 st_leafmtxp%ktmax=ktmax
 if(st_ctl%param(1)>0)  write(mpilog,1000) ' ktmax=',ktmax
 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( ktmax, ktmaxmax, 1, MPI_INTEGER, MPI_MAX,0, icomm, ierr );
 if(st_ctl%param(1)>1 .and. mpinr==0) write(*,1000) 'Maximun of kt=',ktmaxmax
 call MPI_reduce( ktav, ktavav, 1, MPI_INTEGER, MPI_SUM,0, icomm, ierr );
 if(st_ctl%param(1)>1 .and. mpinr==0) write(*,1000) 'Average of kt=',ktavav/st_leafmtxp%nlfkt

 znndense=real(nd)*real(nd)*8.0/1024.0/1024.0
! if(st_ctl%param(1)>0)  write(mpilog,2000) 'dense matrix memory=',znndense,'(Mbyte)'
 memh=(mem8+sum(lnmtx(:)))
 znn=real(memh)*8.0/1024.0/1024.0
! write(*,*) 'mpinr=',mpinr,'H matrix memory=',znn,'(Mbyte)'
 if(st_ctl%param(1)>0)  write(mpilog,2000) '    H matrix memory=',znn,'(Mbyte)'

 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( znn, znnmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX,0, icomm, ierr );
 call MPI_reduce( znn, znnmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN,0, icomm, ierr );
 call MPI_reduce( znn, znnall, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Memory of the H-matrix=',znnall,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Memory compression v.s. dense matrix=',znnall/znndense*100,'(%)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Maximun memory of sub-matrices for a MPI=',znnmax,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Minimun memory of sub-matrices for a MPI=',znnmin,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Minimun memory/Maximun memory=',znnmin/znnmax

 end subroutine HACApK_chk_leafmtx

!***HACApK_count_lntmx
 RECURSIVE subroutine HACApK_count_lntmx(st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 integer*4 :: lnmtx(3),lpmd(*)
 real*8 :: param(*)
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson

 nleaf=param(21)+1; nlmax=param(22)*nofc
 zs=0.0d0
 do id=1,st_cltl%ndim
   if(st_cltl%bmax(id)<st_cltt%bmin(id))then
     zs=zs+(st_cltt%bmin(id)-st_cltl%bmax(id))*(st_cltt%bmin(id)-st_cltl%bmax(id))
   elseif(st_cltt%bmax(id)<st_cltl%bmin(id))then
     zs=zs+(st_cltl%bmin(id)-st_cltt%bmax(id))*(st_cltl%bmin(id)-st_cltt%bmax(id))
   else
   endif
 enddo
 zdistlt=dsqrt(zs)
 zeta=param(51)
 if((st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt).and.(ndl>=nleaf .and. ndt>=nleaf .and. ndl<=nlmax .and. ndt<=nlmax))then
   lnmtx(1)=lnmtx(1)+1
 else
   if(nnsonl==0 .or. nnsont==0 .or. ndl<=nleaf .or. ndt<=nleaf)then
     lnmtx(2)=lnmtx(2)+1
   else
     lnmtx(3)=lnmtx(3)+1
     do il=1,nnsonl
       do it=1,nnsont
         call HACApK_count_lntmx(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc)
       enddo
     enddo
   endif
 endif
 end subroutine HACApK_count_lntmx

!***HACApK_free_leafmtxp
 integer function HACApK_free_leafmtxp(st_leafmtxp)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
   do il=1,st_leafmtxp%nlf
     if(st_leafmtxp%st_lf(il)%ltmtx==1)then
       deallocate(st_leafmtxp%st_lf(il)%a1)
       deallocate(st_leafmtxp%st_lf(il)%a2)
     elseif(st_leafmtxp%st_lf(il)%ltmtx==2)then
       deallocate(st_leafmtxp%st_lf(il)%a1)
     endif
   enddo
   if(associated(st_leafmtxp%st_lf)) deallocate(st_leafmtxp%st_lf)
   HACApK_free_leafmtxp=0
 end function HACApK_free_leafmtxp

!***HACApK_free_lcontrol
 integer function HACApK_free_lcontrol(st_ctl)
 type(st_HACApK_lcontrol) :: st_ctl
   if(associated(st_ctl%lod)) deallocate(st_ctl%lod)
   if(associated(st_ctl%lsp)) deallocate(st_ctl%lsp)
   if(associated(st_ctl%lnp)) deallocate(st_ctl%lnp)
   if(associated(st_ctl%lthr)) deallocate(st_ctl%lthr)
   if(associated(st_ctl%lpmd)) deallocate(st_ctl%lpmd)
   deallocate(st_ctl%param)
   HACApK_free_lcontrol=0
 end function HACApK_free_lcontrol

!***HACApK_generate_leafmtx
 RECURSIVE subroutine HACApK_generate_leafmtx(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(3),lpmd(*)
 real*8 :: param(*)

 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt

 nleaf=param(21)+1; nlmax=param(22)*nofc
! print*,'nlmax=',nlmax; stop
 zs=0.0d0
 do id=1,st_cltl%ndim
   if(st_cltl%bmax(id)<st_cltt%bmin(id))then
     zs=zs+(st_cltt%bmin(id)-st_cltl%bmax(id))*(st_cltt%bmin(id)-st_cltl%bmax(id))
   elseif(st_cltt%bmax(id)<st_cltl%bmin(id))then
     zs=zs+(st_cltl%bmin(id)-st_cltt%bmax(id))*(st_cltl%bmin(id)-st_cltt%bmax(id))
   else
   endif
 enddo
! zdistlt=max(dsqrt(zs)-st_cltl%zwdth/ndl-st_cltt%zwdth/ndt,0.0)
 zdistlt=dsqrt(zs)
 zeta=param(51)
 if((st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt).and.(ndl>=nleaf .and. ndt>=nleaf .and. ndl<=nlmax .and. ndt<=nlmax))then
   nlf=nlf+1
   st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
   st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
   st_leafmtx(nlf)%kt=0
   st_leafmtx(nlf)%ltmtx=1
 else
   if(nnsonl==0 .or. nnsont==0 .or. ndl<=nleaf .or. ndt<=nleaf)then
     nlf=nlf+1
     st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
     st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
     st_leafmtx(nlf)%ltmtx=2
!     allocate(st_leafmtx(nlf)%a1(ndt,ndl))
   else
     do il=1,nnsonl
       do it=1,nnsont
         call HACApK_generate_leafmtx(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf)
       enddo
     enddo
   endif
 endif
 end subroutine HACApK_generate_leafmtx

!***HACApK_qsort_col_leafmtx
 RECURSIVE subroutine HACApK_qsort_col_leafmtx(st_leafmtx,nlf_s,nlf_e)
 type(st_HACApK_leafmtx) :: st_leafmtx(*),st_www
    if(nlf_s>=nlf_e) return
    nl = nlf_s; nr = nlf_e; nlr2=nl+(nr-nl)/2
    nlt=st_leafmtx(nl)%nstrtt; nrt=st_leafmtx(nr)%nstrtt; nlrt=st_leafmtx(nlr2)%nstrtt
    nmid=HACApK_med3(nlt,nrt,nlrt)
!    print*,'nlf_s=',nlf_s,'nlf_e=',nlf_e,'nlr2=',nlr2,'nmid=',nmid
    do
      do while(st_leafmtx(nl)%nstrtt < nmid); nl=nl+1; enddo
      do while(st_leafmtx(nr)%nstrtt > nmid); nr=nr-1; enddo
      if(nl >= nr) exit
      st_www = st_leafmtx(nl); st_leafmtx(nl) = st_leafmtx(nr); st_leafmtx(nr) = st_www
      nl=nl+1; nr=nr-1
    enddo
    call HACApK_qsort_col_leafmtx(st_leafmtx,nlf_s,nl-1)
    call HACApK_qsort_col_leafmtx(st_leafmtx,nr+1 ,nlf_e)
 end subroutine HACApK_qsort_col_leafmtx

!***HACApK_qsort_row_leafmtx
 RECURSIVE subroutine HACApK_qsort_row_leafmtx(st_leafmtx,nlf_s,nlf_e)
 type(st_HACApK_leafmtx) :: st_leafmtx(*),st_www
    if(nlf_s>=nlf_e) return
    nl = nlf_s; nr = nlf_e; nlr2=nl+(nr-nl)/2
    nmid=HACApK_med3(st_leafmtx(nl)%nstrtl,st_leafmtx(nr)%nstrtl,st_leafmtx(nlr2)%nstrtl)
!    nmid=st_leafmtx(nlr)%nstrtl
!    print*,'nlf_s=',nlf_s,'nlf_e=',nlf_e,'nlr2=',nlr2,'nmid=',nmid
    do
      do while(st_leafmtx(nl)%nstrtl < nmid); nl=nl+1; enddo
      do while(st_leafmtx(nr)%nstrtl > nmid); nr=nr-1; enddo
      if(nl >= nr) exit
      st_www = st_leafmtx(nl); st_leafmtx(nl) = st_leafmtx(nr); st_leafmtx(nr) = st_www
      nl=nl+1; nr=nr-1
    enddo
    call HACApK_qsort_row_leafmtx(st_leafmtx,nlf_s,nl-1)
    call HACApK_qsort_row_leafmtx(st_leafmtx,nr+1 ,nlf_e)
 end subroutine HACApK_qsort_row_leafmtx

!***HACApK_free_st_clt
  RECURSIVE subroutine HACApK_free_st_clt(st_clt)
   type(st_HACApK_cluster) :: st_clt
   nnson=st_clt%nnson
   do ic=1,nnson
     call HACApK_free_st_clt(st_clt%pc_sons(ic))
   enddo
   deallocate(st_clt%bmin,st_clt%bmax)
   deallocate(st_clt%pc_sons)
 end subroutine HACApK_free_st_clt

!***HACApK_generate_cluster
  type(st_HACApK_cluster) function HACApK_generate_cluster(nmbr,ndpth,nstrt,nsize,ndim,nson)
   type(st_HACApK_cluster) :: st_clt
   nmbr=nmbr+1
   st_clt%nstrt=nstrt; st_clt%nsize=nsize; st_clt%ndim=ndim; st_clt%nnson=nson
   st_clt%nmbr=nmbr; st_clt%ndpth=ndpth
   allocate(st_clt%pc_sons(nson))
   HACApK_generate_cluster=st_clt
 endfunction

!***HACApK_bndbox
 RECURSIVE subroutine HACApK_bndbox(st_clt,zgmid,lod,nofc)
 type(st_HACApK_cluster) :: st_clt
 real*8 :: zgmid(nofc,st_clt%ndim)
 integer*4 :: lod(nofc)
 do ic=1,st_clt%nnson
   if(ic==1)then; l=1;
   else; l=l+st_clt%pc_sons(ic-1)%nsize;
   endif
   call HACApK_bndbox(st_clt%pc_sons(ic),zgmid,lod(l),nofc);
 enddo
 ndim=st_clt%ndim
 allocate(st_clt%bmin(ndim),st_clt%bmax(ndim))
 if(st_clt%nnson == 0)then
   do id=1,ndim
      st_clt%bmin(id)=zgmid(lod(1),id); st_clt%bmax(id)=zgmid(lod(1),id)
   enddo
   do id=1,ndim
     do il=2,st_clt%nsize
       if(zgmid(lod(il),id) < st_clt%bmin(id)) st_clt%bmin(id) = zgmid(lod(il),id)
       if(st_clt%bmax(id) < zgmid(lod(il),id)) st_clt%bmax(id) = zgmid(lod(il),id)
     enddo
   enddo
 else
   do id=1,ndim
      st_clt%bmin(id)=st_clt%pc_sons(1)%bmin(id); st_clt%bmax(id)=st_clt%pc_sons(1)%bmax(id)
   enddo
   do il=2,st_clt%nnson
     do id=1,ndim
       if(st_clt%pc_sons(il)%bmin(id) < st_clt%bmin(id)) st_clt%bmin(id)=st_clt%pc_sons(il)%bmin(id)
       if(st_clt%bmax(id) < st_clt%pc_sons(il)%bmax(id)) st_clt%bmax(id)=st_clt%pc_sons(il)%bmax(id)
     enddo
   enddo
 endif
 zwdth=(st_clt%bmax(1)-st_clt%bmin(1))*(st_clt%bmax(1)-st_clt%bmin(1))
 do id=2,ndim
   zwdth=zwdth+(st_clt%bmax(id)-st_clt%bmin(id))*(st_clt%bmax(id)-st_clt%bmin(id))
 enddo
 st_clt%zwdth=dsqrt(zwdth)
 end subroutine HACApK_bndbox

!***HACApK_generate_cbitree
 RECURSIVE subroutine HACApK_generate_cbitree(st_clt,zgmid,param,lpmd,lod,ndpth,ndscd,nsrt,nd,md,ndim,nclst)
   type(st_HACApK_cluster) :: st_clt
   real*8 :: zgmid(md,ndim)
   real*8,dimension(:),allocatable :: zlmin,zlmax
   integer*4 :: lod(md),lpmd(*)
   real*8 :: param(*)
   minsz=param(21)
   ndpth=ndpth+1
!   ndscd=ndscd+1
!  if(i>26) stop
!   print*,''
!   print*,'nsrt=',nsrt,' nd=',nd
   if(nd <= minsz)then
     nson=0
!     nclst=nclst+1
     st_clt=HACApK_generate_cluster(nclst,ndpth,nsrt,nd,ndim,nson)
   else
     allocate(zlmin((ndim)),zlmax((ndim)))
     do id=1,ndim
        zlmin(id)=zgmid(lod(1),id); zlmax(id)=zlmin(id)
       do il=2,nd; zg=zgmid(lod(il),id)
         if    (zg<zlmin(id))then; zlmin(id)=zg;
         elseif(zlmax(id)<zg)then; zlmax(id)=zg;
         endif
       enddo
     enddo
!     print*,'zlmin=',zlmin
!     print*,'zlmax=',zlmax

     zdiff=zlmax(1)-zlmin(1); ncut = 1
     do id=1,ndim
       zidiff=zlmax(id)-zlmin(id)
       if(zidiff>zdiff)then
         zdiff =zidiff; ncut=id
       endif
     enddo
     zlmid= (zlmax(ncut)+zlmin(ncut))/2
!     print*,'ncut=',ncut, '; zlmid=',zlmid

    nl = 1; nr = nd;
    do while(nl < nr)
      do while(nl < nd .and. zgmid(lod(nl),ncut) <= zlmid); nl=nl+1; enddo
      do while(nr >= 0 .and. zgmid(lod(nr),ncut) > zlmid) ; nr=nr-1; enddo
      if(nl < nr) then; nh = lod(nl); lod(nl) = lod(nr); lod(nr) = nh; endif
    enddo

! print*,'nd=',nd, ';ncut=',ncut, '; nsrt=',nsrt,'; nl=',nl

    nson=2
    st_clt=HACApK_generate_cluster(nclst,ndpth,nsrt,nd,ndim,nson)
    nsrt1=nsrt; nd1=nl-1
    call HACApK_generate_cbitree(st_clt%pc_sons(1),zgmid,param,lpmd,lod,ndpth,ndscd,nsrt1,nd1,md,ndim,nclst)
    ndpth=ndpth-1
!    ndscd=ndscd+st_clt%pc_sons(1)%ndscd
    nsrt1=nsrt+nl-1; nd1=nd-nl+1
    call HACApK_generate_cbitree(st_clt%pc_sons(2),zgmid,param,lpmd,lod(nl),ndpth,ndscd,nsrt1,nd1,md,ndim,nclst)
    ndpth=ndpth-1
!    ndscd=ndscd+st_clt%pc_sons(2)%ndscd
  endif
  st_clt%ndscd=nd
 end subroutine HACApK_generate_cbitree

!***HACApK_cal_matnorm
 subroutine HACApK_cal_matnorm(znrm,st_bemv,lpmd,nd)
 type(st_HACApK_calc_entry) :: st_bemv
 integer*4 :: lpmd(:)
 znrm=0.0d0
 do il=lpmd(6),lpmd(7)
   zz=HACApK_entry_ij(il,il,st_bemv)
   znrm=znrm+zz*zz
 enddo
 endsubroutine

!***HACApK_impi_allgv
 subroutine HACApK_impi_allgv(zau,lpmd,nd)
 include 'mpif.h'
 real*8 :: zau(*)
 real*8,dimension(:),allocatable :: wws,wwr
 integer*4 :: lpmd(*),ISTATUS(MPI_STATUS_SIZE),isct(2),irct(2)
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 ndnr_s=lpmd(6); ndnr_e=lpmd(7); ndnr=lpmd(5)
 allocate(wws(nd),wwr(nd))
 wws(1:ndnr)=zau(ndnr_s:ndnr_e)
 ncdp=mod(mpinr+1,nrank)
 ncsp=mod(mpinr+nrank-1,nrank)

! write(mpilog,1000) 'ndnr_s=',ndnr_s,' ndnr_e=',ndnr_e,' ndnr=',ndnr,' ncdp=',ncdp,' ncsp=',ncsp

 isct(1)=ndnr; isct(2)=ndnr_s
 do ic=1,nrank-1
   call MPI_SENDRECV(isct,2,MPI_INTEGER,ncdp,1, &
                     irct,2,MPI_INTEGER,ncsp,1,icomm,ISTATUS,ierr)
   call MPI_SENDRECV(wws,isct,MPI_DOUBLE_PRECISION,ncdp,1, &
                     wwr,irct,MPI_DOUBLE_PRECISION,ncsp,1,icomm,ISTATUS,ierr)
   zau(irct(2):irct(2)+irct(1)-1)=zau(irct(2):irct(2)+irct(1)-1)+wwr(:irct(1))
   wws(:irct(1))=wwr(:irct(1))
   isct=irct
 enddo
 end subroutine

endmodule m_HACApK_base
