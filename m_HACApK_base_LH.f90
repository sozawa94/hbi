!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 3.0.0                                             *
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
!C  added functions related to ACA+ to HACApK1.0.0 on Nov. 2016
!C  corrected the allocation for st_ctl%lthr on Nov. 2016
!C  added a function related to HACApK_view to HACApK1.1.0 on May 2017
!C  added a function related to writting H-matrix to HACApK1.2.0 on May 2017
!C  added functions related to Block clustering to HACApK1.2.0 on May 2017
!C  added functions related to Lattice H-matrix to HACApK1.3.0 on Aug 2017
!C  added functions related to Uniform Lattice size to HACApK2.3.0 on Aug 2018
!C  last modified by Akihiro Ida on Aug 2018
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
    integer*4 lttcl,lttct
    
    integer*4 nlf ! number of leaves(sub-matrices) in the MPI process
    !!!
    integer*8 a1size !!!
    !!!
    real*8,pointer :: a1(:,:)=>null()
    real*8,pointer :: a2(:,:)=>null()
    !!!
    type(st_HACApK_leafmtx),pointer :: st_lf(:)=>null()
  end type st_HACApK_leafmtx

!*** type :: st_HACApK_leafmtxp
  type st_HACApK_leafmtxp
    integer*4 nd    ! number of unknowns of whole matrix
    integer*4 nlf   ! number of leaves(sub-matrices) in the MPI process
    integer*4 nlfkt ! number of low-rank sub matrices in the MPI process
    integer*4 ktmax
    integer*4 nbl           ! number of blocks for MPI assignment
    integer*4 nlfalt        ! number of leaves in row(column) of whole matrix
    integer*4 nlfl, nlft    ! number of leaves in row and column in the MPI process
    integer*4 ndlfs, ndtfs  ! vector sizes in the MPI process
    !!!
    integer*4 st_lf_stride 
    !!!
    type(st_HACApK_leafmtx),pointer :: st_lf(:)=>null()
    !!!
    integer*8,pointer :: lnlfl2g(:,:)=>null()
    !
    integer*4,pointer :: lbstrtl(:)=>null() ! Start points of each block in row
    integer*4,pointer :: lbstrtt(:)=>null() ! Start points of each block in column
    !
    integer*4,pointer :: lbndl(:)=>null() ! vector sizes of each block in row
    integer*4,pointer :: lbndt(:)=>null() ! vector sizes of each block in column
    !
    integer*4,pointer :: lbndlfs(:)=>null() ! vector sizes of each MPI process in row
    integer*4,pointer :: lbndtfs(:)=>null() ! vector sizes of each MPI process in column
    !
    integer*4,pointer :: lbl2t(:)=>null() ! bit vector for recieving data on each MPI process
  end type st_HACApK_leafmtxp

!*** type :: st_HACApK_lf9lttc
  type st_HACApK_lf9lttc
    integer*4 nlf   ! number of leaves(sub-matrices) in the MPI process
    integer*4 ndl,ndt  ! lattice sizes
    type(st_HACApK_leafmtx),pointer :: st_lf(:)=>null()
  end type st_HACApK_lf9lttc

!*** type :: st_HACApK_lcontrol
  type :: st_HACApK_lcontrol
    integer*4 :: lf_umpi
    integer*4 time_offset  ! for C-interface
    integer*4 lpmd_offset  ! for C-interface
    integer*4 lod_offset   ! for C-interface
    integer*4 lsp_offset   ! for C-interface
    integer*4 lnp_offset   ! for C-interface
    integer*4 lthr_offset  ! for C-interface
    real*8,   pointer :: param(:) =>null()
    real*8,   pointer :: time(:)  =>null()
    integer*4,pointer :: lpmd(:)  =>null()
    integer*4,pointer :: lod(:)   =>null()
    integer*4,pointer :: lsp(:)   =>null()
    integer*4,pointer :: lnp(:)   =>null()
    integer*4,pointer :: lthr(:)  =>null()
  end type st_HACApK_lcontrol

!*** type :: st_HACApK_LHp
  type st_HACApK_LHp
    integer*4 nd    ! number of unknowns of whole matrix
    integer*4 nlf   ! number of leaves(sub-matrices) in the MPI process
    integer*4 nlfalt  ! number of lattice blocks in row(column) of whole matrix
    integer*4 nlfl, nlft    ! number of lattice blocks in row and column in the MPI process
    integer*4 ndlfs, ndtfs  ! vector sizes in the MPI process
    type(st_HACApK_lf9lttc),pointer :: st_lfp(:,:)=>null()
    integer*8,pointer :: lnlfl2g(:,:)=>null()
    integer*4,pointer :: lbstrtl(:)=>null() ! Start points of lattice block in row in the whole matrix
    integer*4,pointer :: lbstrtt(:)=>null() ! Start points of lattice block in column in the whole matrix
    integer*4,pointer :: lbndlfs(:)=>null() ! row vector sizes of each MPI process
    integer*4,pointer :: lbndtfs(:)=>null() ! column vector sizes of each MPI process in column
    integer*4,pointer :: latticel(:)=>null() ! row index in whole matrix
    integer*4,pointer :: latticet(:)=>null() ! column index in whole matrix
    integer*4,pointer :: lbndcsl(:)=>null() ! start point of lattice block in row in the MPI process
    integer*4,pointer :: lbndcst(:)=>null() ! start point of lattice block in column in the MPI process
  end type st_HACApK_LHp

!*** type :: st_HACApK_latticevec
  type st_HACApK_latticevec
    integer*4 ndc    ! number of unknowns in the MPI process
    integer*4 nlfc    ! number of blocks in the MPI process
    integer*4,pointer :: lbstrtc(:)=>null() ! Start points of each block in the MPI process in whole vector
    integer*4,pointer :: lbndc(:)=>null() ! Vector sizes of each block in the MPI process
    integer*4,pointer :: lodc(:)=>null() ! Permutation info.
    real*8,pointer :: vs(:)=>null() ! Entry values of the vector in the MPI process
    integer*4 ncomm    ! number of neighbor Umag compornents in the MPI process
    integer*4 ncoup    ! number of neighbor Umag array in the MPI process
    integer*4 nel      ! number of Umag array in the MPI process (==ndc/3)
    integer*4 nnb      ! number of neighbor Umag index array in the MPI process
    integer*4 nlocal   ! number of local Umag compornents in the MPI process
    integer*4,pointer :: iord(:)=>null() ! Index of 'vs' vector with ordering.
    integer*4,pointer :: commnum(:)=>null() ! Neighbor Umag elements for all MPI process.
    integer*4,pointer :: commdisp(:)=>null() ! for MPI_Gatherv.
    integer*4,pointer :: commlist(:)=>null() ! Send list for other process.
    integer*4,pointer :: nblist(:)=>null() ! Neighbor Umag index in the MPI process.
    integer*4,pointer :: nbptr(:)=>null() ! Neighbor Umag index pointer in the MPI process.
    integer*4,pointer :: nbindex(:)=>null() ! Neighbor Umag index in the MPI process.
    integer*4,pointer :: localnum(:)=>null() ! Local UMAG elements in the MPI process.
    integer*4,pointer :: localdisp(:)=>null() ! for MPI_Gatherv.
    real*8,pointer :: vsdf(:)=>null() ! diff of 'vs'
    real*8,pointer :: vssd(:)=>null() ! Send values of the vector in the MPI process
    real*8,pointer :: vsrv(:)=>null() ! Recv values of the vector in the MPI process
    real*8,pointer :: vslc(:)=>null() ! Local value
    real*8,pointer :: nbcoef(:)=>null() ! Exchange coupling coefficient
  end type st_HACApK_latticevec

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

 allocate(st_ctl%param(100),st_ctl%time(10))
 st_ctl%param(1:100)=0.0
 st_ctl%param(1) =1;        ! Print : 0:Only Error 1:STD 2:Dubug
 st_ctl%param(7) =0;        ! 1;Make files for HACApK_view, 0; not make
 st_ctl%param(8) =20;       ! > 10:BLR, 20:LH, other:Leafmtx
 st_ctl%param(9) =1;        ! 1:load balancer
 st_ctl%param(10)=1;        ! 1:fulfill the matrix,  0: not fulfill
 st_ctl%param(11)=0;        ! 1:check accuracy of H-matrix 0: not check
 st_ctl%param(12)=0;        ! 1:write down the H-matrix to file 0: not write
 st_ctl%param(21)=15;       ! > cluster : leaf size 15
 st_ctl%param(22)=1.0;      ! cluster : max leaf size 1.0*nffc
 st_ctl%param(23)=0;        ! 0:traditional clustering 1:Uniform size of block or lattice
 st_ctl%param(41)=0;        ! > p of p-by-q grid. if 0, then p=floor(sqrt(#MPI processes)). if -1 then p=#MPI processes
 st_ctl%param(42)=0;        ! > BLR:block size  LH:latice size. if 0, it is set to N/param(43)/param(41)
 st_ctl%param(43)=6;       ! > BLR, LH : Number of rows in each lattice, used only if param(42)=0 
 if(nrank==1) st_ctl%param(43)=1
 st_ctl%param(51)=2.0;      ! H-matrix : decision param of distance 2.0
 st_ctl%param(52)=0;        ! H-matrix : 0:weak admissibility 1:strong
 st_ctl%param(53)=100;      ! H-matrix : maximun depth of block cluster tree 100
 st_ctl%param(54)=0;        ! H-matrix : 0: quad(Ver1.2), 1:BLR, 2: quad or bi
 st_ctl%param(60)=2         ! 1:ACA,  2:ACA+
 st_ctl%param(61)=1         ! ACA norm 1:MREM,  2:test, 3:norm
 st_ctl%param(62)=7         ! ACA : predictive average of k
 st_ctl%param(63)=1000;     ! ACA : k-max of R_k-matrix 30
 st_ctl%param(64)=1;        ! ACA : minimun kt
 st_ctl%param(72)=1.0e-9;   ! ACA_EPS
 st_ctl%param(83)=10;       ! solver : maximum iterative number
 st_ctl%param(85)=1;        ! solver : 1:BiCGSTAB, 2:GCR(m)
 st_ctl%param(87)=8;        ! solver : number of iteration for reset
 st_ctl%param(91)=1.0e-6    ! Required accuracy iterative solver
 st_ctl%param(98)=50        ! Measure the time of Ax; number of trial of Ax
 st_ctl%param(99)=10        ! Measure the time of Ax; number of iteration

 allocate(st_ctl%lpmd(50)); st_ctl%lpmd(:)=0
 st_ctl%lpmd(1)=icomm; ! MPI communicator for all MPI processors
 st_ctl%lpmd(2)=nrank; ! # of MPI processors
 st_ctl%lpmd(3)=irank; ! my rank
 st_ctl%lpmd(4)=20;    ! Unit number to access log files
 ! st_ctl%lpmd(5) : 
 ! st_ctl%lpmd(6) : 
 ! st_ctl%lpmd(7) : 
 ! st_ctl%lpmd(11) : 
 ! st_ctl%lpmd(12) : 
 ! st_ctl%lpmd(20) : # of threads
 ! st_ctl%lpmd(31) : MPI communicator for MPI processors in the same row when using blr structure
 ! st_ctl%lpmd(32) : # of MPI processors in the same row when using blr structure
 ! st_ctl%lpmd(33) : my rank in the same row when using blr structure
 ! st_ctl%lpmd(35) : MPI communicator for MPI processors in the same column when using blr structure
 ! st_ctl%lpmd(36) : # of MPI processors in the same column when using blr structure
 ! st_ctl%lpmd(37) : my rank in the same column when using blr structure
 
 nthr=1
!$omp parallel
  nthr = omp_get_num_threads()
  if(nthr>0) st_ctl%lpmd(20)=nthr
!$omp end parallel
 allocate(st_ctl%lod(nd),st_ctl%lthr(nthr+1),st_ctl%lnp(nrank),st_ctl%lsp(nrank),stat = ierr)
 
 call MPI_Barrier( icomm, ierr )

 if(st_ctl%param(1)>1)then
   write(logfile,'(a,i4.4,a)') 'log',irank,'.txt'
   open(st_ctl%lpmd(4),file=logfile)
 else
   st_ctl%lpmd(4)=0
 endif

 st_bemv%nd=nd
 if(st_bemv%lp61==1) then; endif
 if(ierr/=0)then; goto 9999; endif 
 
9999 continue
 HACApK_init=lrtrn
 endfunction

!***HACApK_param_set
integer function HACApK_param_set( st_ctl, param )
 implicit none
 type(st_HACApK_lcontrol) :: st_ctl
 real*8 :: param(:)
 st_ctl%param(1) = param(1) ; ! Print : 0:Only Error 1:STD 2:Dubug
 st_ctl%param(7) = param(7) ; ! 1;Make files for HACApK_view, 0; not make
 st_ctl%param(8) = param(8) ; ! 10:BLR, 20:LH, other:Leafmtx
 st_ctl%param(9) = param(9) ; ! 1:load balancer
 st_ctl%param(10)= param(10); ! 1:fulfill the matrix  0: not fulfill
 st_ctl%param(11)= param(11); ! 1:check accuracy of H-matrix 0: not check
 st_ctl%param(12)= param(12); ! 1:write down the H-matrix to file 0: not write
 st_ctl%param(21)= param(21); ! cluster : leaf size 15
 st_ctl%param(22)= param(22); ! cluster : max leaf size 1.0*nffc
 st_ctl%param(23)= param(23); ! 0:traditional clustering 1:Uniform size of block or lattice
 st_ctl%param(41)= param(41); ! p of p-by-q grid. if 0, then p=floor(sqrt(#MPI processes)). if -1 then p=#MPI processes
 st_ctl%param(42)= param(42); ! BLR:block size  LH:latice size. if 0, it is set to N/param(43)/param(41)
 st_ctl%param(43)= param(43); ! BLR, LH : Number of rows in each lattice, used only if param(42)=0 
 st_ctl%param(51)= param(51); ! H-matrix : dicision param of distance 2.0
 st_ctl%param(52)= param(52); ! H-matrix : 0:weak admissibility 1:strong
 st_ctl%param(53)= param(53); ! H-matrix : maximun depth of block cluster tree 100
 st_ctl%param(54)= param(54); ! H-matrix : 0: quad(Ver1.2), 1:BLR, 2: quad or bi
 st_ctl%param(60)= param(60); ! 1:ACA  2:ACA+
 st_ctl%param(61)= param(61); ! ACA norm 1:MREM  2:test 3:norm
 st_ctl%param(62)= param(62); ! ACA : predictive average of k
 st_ctl%param(63)= param(63); ! ACA : k-max of R_k-matrix 30
 st_ctl%param(64)= param(64); ! ACA : minimun kt
 st_ctl%param(72)= param(72); ! ACA_EPS
 st_ctl%param(83)= param(83); ! solver : maximum iterative number
 st_ctl%param(85)= param(85); ! solver : 1:BiCGSTAB 2:GCR(m)
 st_ctl%param(87)= param(87); ! solver : number of iteration for reset
 st_ctl%param(91)= param(91); ! Required accuracy iterative solver
 st_ctl%param(98)= param(98); ! Measure the time of Ax; number of trial of Ax
 st_ctl%param(99)= param(99); ! Measure the time of Ax; iterative number
 HACApK_param_set=0
end function HACApK_param_set

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

!***HACApK_generate_frame_LH
 subroutine HACApK_generate_frame_LH(st_LHp,st_bemv,st_ctl,gmid,lnmtx,nofc,nffc,ndim)
 include 'mpif.h'
 type(st_HACApK_cluster) :: st_clt
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_LHp) :: st_LHp
 type(st_HACApK_leafmtx),dimension(:), allocatable :: st_leafmtx,st_leafmtx_lcl
 real*8 :: gmid(nofc,ndim)
 integer*8 :: mem8,nlfall
 integer*4 :: lnmtx(4)
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
   print*,'HACApK_generate_frame_LH: ierr=',ierr
 endif
 do il=1,nofc
   lodfc(il)=il
 enddo 
 npgl=param(41); 
 if(npgl==0) then
   npgl = sqrt(real(nrank))
   npgt = nrank/npgl
   do while(npgl*npgt .ne. nrank)
     npgl = npgl - 1
     npgt = nrank/npgl
   enddo
   if (mpinr==0) then
     write(*,*) 'nrank=',nrank,'(',npgl,'x',npgt,')'
   endif
 elseif(npgl<0)then
   npgl=nrank
 endif
! if(param(42)==0) param(42)=sqrt(real(nd))
 if(param(42)==0) param(42)=nd/param(43)/npgl
 if(param(42)<param(21))then
  if(mpinr==0) print*, 'sub HACApK_generate_frame_LH; param(42)=',param(42),' param(21)=',param(21)
  if(mpinr==0) print*, 'param(42)(lattice size) is smaller than param(21)(leaf size)! Change to BLR!'
!  if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; param(42)(block size) must be larger than param(21)(leaf size) !!!'
! goto 9999
 endif

!!!!!!!!!!!!!!!!!! start clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nsrt=1; ndf=nofc; nclst=0; ndpth=0; ndscd=0
! call HACApK_generate_cbitree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst)
 if(st_ctl%param(23)==1)then
   call HACApK_generate_tctree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst,nblall)
 else
   call HACApK_generate_cbitree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst)
   ndpth=0; lnmtx(1:4)=0
   call HACApK_count_blrnmb(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
   nblall=lnmtx(4)
 endif
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'No. of cluster=',nclst
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'No. of cluster=',nclst

 call HACApK_bndbox(st_clt,gmid,lodfc,nofc)
 do il=1,nofc
   do ig=1,nffc
     is=ig+(il-1)*nffc
     lod(is)=lodfc(il)
   enddo
 enddo
!!!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 allocate(st_leafmtx(nblall))
 nlfalt=sqrt(real(nblall)); st_LHp%nlfalt=nlfalt
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'Number of MPI_Blocks=',nblall,'; sqrt(nblall)=',nlfalt
 ndpth=0; lnmtx(1:4)=0
 if(st_ctl%param(23)==1)then
   call HACApK_count_LH(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
   lnmtx(4)=nblall
 else
   call HACApK_count_blrleaf(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 endif

 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'No. of nsmtx',lnmtx(1:4)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'   1:Rk-matrix 2: dense-mat 3:H-matrix 4:MPI_Block'
 
 if(st_ctl%param(23)==1)then
   ndpth=0; lnmtx(1:4)=0
   call HACApK_generate_LH(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nblall,ndpth)
 else
   nblall=0; ndpth=0; lnmtx(1:4)=0
   call HACApK_generate_blrleaf(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nblall,ndpth)
 endif
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_LH; HACApK_generate_leafmtx end'
 call HACApK_sort_leafmtx(st_leafmtx,nblall)
 do ip=1,nblall
   if(st_leafmtx(ip)%ltmtx==4) call HACApK_sort_leafmtx(st_leafmtx(ip)%st_lf,st_leafmtx(ip)%nlf)
 enddo
 call HACApK_free_st_clt(st_clt)
  
!!!!!!!!!!!!!!!!!! start MPI load balance  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 npgt=nrank/npgl
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) ' npgl=',npgl,' npgt=',npgt
 if(st_ctl%param(1)>1) write(mpilog,1000) ' npgl=',npgl,' npgt=',npgt
 if(npgt>nlfalt .or. npgl>nlfalt)then
   call MPI_Barrier( icomm, ierr )
   if(mpinr==0) write(6,*) 'Error: HACApK_generate_frame_LH; Too few blocks compared with #MPI !!!'; goto 9999
 endif
 if(npgt*npgl/=nrank)then
   call MPI_Barrier( icomm, ierr )
   if(mpinr==0) write(6,*) 'Error: HACApK_generate_frame_LH; Invalid processor grid!!!'; goto 9999
 endif
 
! Split MPI communicator
  ikey=0; iclr=mpinr/npgt
  call MPI_COMM_SPLIT(icomm, iclr, ikey, icommn, ierr); st_ctl%lpmd(31)=icommn
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; MPI_COMM_SPLIT failed !!!'; goto 9999
  endif
  call MPI_Comm_size ( icommn, nrank, ierr ); st_ctl%lpmd(32)=nrank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; MPI_Comm_size failed !!!'; goto 9999
  endif
  call MPI_Comm_rank ( icommn, irank, ierr ); st_ctl%lpmd(33)=irank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; MPI_Comm_rank failed !!!'; goto 9999
  endif
  ikey=0; iclr=mod(mpinr,npgt)
  call MPI_COMM_SPLIT(icomm, iclr, ikey, icommn, ierr); st_ctl%lpmd(35)=icommn
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; MPI_COMM_SPLIT failed !!!'; goto 9999
  endif
  call MPI_Comm_size ( icommn, nrank, ierr ); st_ctl%lpmd(36)=nrank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; MPI_Comm_size failed !!!'; goto 9999
  endif
  call MPI_Comm_rank ( icommn, irank, ierr ); st_ctl%lpmd(37)=irank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_LH; MPI_Comm_rank failed !!!'; goto 9999
  endif

 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,'; irank_t=',st_ctl%lpmd(33),'; irank_l=',st_ctl%lpmd(37)
 if(st_ctl%param(1)>1) write(mpilog,1000) 'irank_t=',st_ctl%lpmd(33),'; nrank_t=',st_ctl%lpmd(32)
 if(st_ctl%param(1)>1) write(mpilog,1000) 'irank_l=',st_ctl%lpmd(37),'; nrank_l=',st_ctl%lpmd(36)

! stop
 
 nlft=nlfalt/npgt; nlfth=mod(nlfalt,npgt); mpinrth=mod(mpinr,npgt)
 if(mpinrth<nlfth) nlft=nlft+1
 nlfl=nlfalt/npgl; nlflh=mod(nlfalt,npgl); mpinrlh=mpinr/npgt
 if(mpinrlh<nlflh) nlfl=nlfl+1
 nbl=nlfl*nlft; st_LHp%nlfl=nlfl; st_LHp%nlft=nlft
 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,'; nbl=',nbl,'; nblall=',nblall
 if(st_ctl%param(1)>1) write(mpilog,1000) 'No. of blocks; nbl=',nbl,'; row=',nlfl,'; column=',nlft,'; global nbl=',nblall
 
 
   call MPI_Barrier( icomm, ierr )
! stop
 
 nrank_t=st_ctl%lpmd(32); nrank_l=st_ctl%lpmd(36)
  
 allocate(st_LHp%lbndlfs(0:nrank_l-1),st_LHp%lbndtfs(0:nrank_t-1)); st_LHp%lbndlfs=0; st_LHp%lbndtfs=0
 do il=0,nlfalt-1; is=nlfalt*il+1
   ilh=mod(il,npgl)
   st_LHp%lbndlfs(ilh)=st_LHp%lbndlfs(ilh)+st_leafmtx(is)%ndl
 enddo 
 do it=0,nlfalt-1; is=it+1
   ith=mod(it,npgt)
   st_LHp%lbndtfs(ith)=st_LHp%lbndtfs(ith)+st_leafmtx(is)%ndt
 enddo 
 if(.false.) then
! if(mpinr==0) then
   print*,'lbstrtl='
   print*,st_LHp%lbstrtl
   print*,'lbstrtt='
   print*,st_LHp%lbstrtt
   print*,'lbndlfs='
   print*,st_LHp%lbndlfs
   print*,'lbndtfs='
   print*,st_LHp%lbndtfs
 endif
 
  allocate(st_LHp%lbstrtl(nlfl),st_LHp%lbstrtt(nlft))
  allocate(st_LHp%latticel(nlfl),st_LHp%latticet(nlft))
  allocate(st_LHp%lbndcsl(nlfl+1),st_LHp%lbndcst(nlft+1))
  allocate(st_LHp%lnlfl2g(nlft,nlfl),st_LHp%st_lfp(nlfl,nlft)); 
  ip=0; nlf=0
  ndlfs=0; st_LHp%lbndcsl(1)=1
  ndtfs=0; st_LHp%lbndcst(1)=1
  do il=0,nlfalt-1; do it=0,nlfalt-1; is=it+nlfalt*il+1
    ilh=mod(il,npgl); ith=mod(it,npgt)
    ipgclr=ith+ilh*npgt
    if(ipgclr==mpinr)then
      ilf=ip/nlft+1; itf=mod(ip,nlft)+1; ip=ip+1; 
      
      
!!!      st_LHp%lnlfl2g(itf,ilf)=is
      st_LHp%lbstrtl(ilf)=st_leafmtx(is)%nstrtl
      st_LHp%lbstrtt(itf)=st_leafmtx(is)%nstrtt
      st_LHp%latticel(ilf)=il+1
      st_LHp%latticet(itf)=it+1
      st_LHp%st_lfp(ilf,itf)%ndl=st_leafmtx(is)%ndl
      st_LHp%st_lfp(ilf,itf)%ndt=st_leafmtx(is)%ndt
!!!      write(mpilog,*) 'ilf=',ilf,'; itf=',itf,'; lbndcsl_i=',st_LHp%lbndcsl(ilf), '; ndl=',st_LHp%st_lfp(ilf,1)%ndl
      st_LHp%lbndcsl(ilf+1)=st_LHp%lbndcsl(ilf)+st_LHp%st_lfp(ilf,1)%ndl
      st_LHp%lbndcst(itf+1)=st_LHp%lbndcst(itf)+st_LHp%st_lfp(1,itf)%ndt
      st_LHp%st_lfp(ilf,itf)%nlf=0
      if(st_leafmtx(is)%ltmtx==1)then
        nlf=nlf+1
        allocate(st_LHp%st_lfp(ilf,itf)%st_lf(1))
      elseif(st_leafmtx(is)%ltmtx==2)then
        nlf=nlf+1
        allocate(st_LHp%st_lfp(ilf,itf)%st_lf(1))
      else
        nlf=nlf+st_leafmtx(is)%nlf
        allocate(st_LHp%st_lfp(ilf,itf)%st_lf(st_leafmtx(is)%nlf))
      endif
    endif
  enddo; enddo
  st_LHp%nlf=nlf; st_LHp%ndlfs=st_LHp%lbndcsl(nlfl+1)-1; st_LHp%ndtfs=st_LHp%lbndcst(nlft+1)-1
  
!  print*,'mpinr=',mpinr,'; nlf=',nlf
  if(nlf<1) then
    print*, 'Error: sub HACApK_generate_frame_blrleaf; nlf<1 !!!; mpinr=',mpinr; goto 9999
  endif

 if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,'; ndlfs=',st_LHp%ndlfs,'; ndtfs=',st_LHp%ndtfs
 if(st_ctl%param(1)>1) write(mpilog,1000) 'Vector sizes; nd=',nd,'; ndlfs=',st_LHp%ndlfs,'; ndtfs=',st_LHp%ndtfs
 
 return
 9999 continue
 stop
 end subroutine

!***HACApK_generate_frame_blrleaf
 subroutine HACApK_generate_frame_blrleaf(st_leafmtxp,st_bemv,st_ctl,gmid,lnmtx,nofc,nffc,ndim)
 include 'mpif.h'
 type(st_HACApK_cluster) :: st_clt
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_leafmtx),dimension(:), allocatable :: st_leafmtx,st_leafmtx_lcl
 real*8 :: gmid(nofc,ndim)
 integer*8 :: mem8,nlfall
 integer*4 :: lnmtx(4)
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
   print*,'HACApK_generate_frame_blrleaf: ierr=',ierr
 endif
 do il=1,nofc
   lodfc(il)=il
 enddo 
 npgl=param(41); 
 if(npgl==0) then
   npgl = sqrt(real(nrank))
   npgt = nrank/npgl
   do while(npgl*npgt .ne. nrank)
     npgl = npgl - 1
     npgt = nrank/npgl
   enddo
   if (mpinr==0) then
     write(*,*) 'nrank=',nrank,'(',npgl,'x',npgt,')'
   endif
 elseif(npgl<0)then
   npgl=nrank
 endif
! if(param(42)==0) param(42)=sqrt(real(nd))
 if(param(42)==0) param(42)=nd/param(43)/npgl
 if(param(42)<param(21))then
  if(mpinr==0) print*, 'sub HACApK_generate_frame_blrleaf; param(42)=',param(42),' param(21)=',param(21)
  if(mpinr==0) print*, 'param(42)(lattice size) is smaller than param(21)(leaf size)! Change to BLR!'
!  if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; param(42)(block size) must be larger than param(21)(leaf size) !!!'
! goto 9999
 endif

!!!!!!!!!!!!!!!!!! start clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nsrt=1; ndf=nofc; nclst=0; ndpth=0; ndscd=0
! call HACApK_generate_cbitree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst)
 if(st_ctl%param(23)==1)then
   call HACApK_generate_tctree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst,nblall)
 else
   call HACApK_generate_cbitree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst)
   ndpth=0; lnmtx(1:4)=0
   call HACApK_count_blrnmb(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
   nblall=lnmtx(4)
 endif
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'No. of cluster=',nclst
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'No. of cluster=',nclst

 call HACApK_bndbox(st_clt,gmid,lodfc,nofc)
 do il=1,nofc
   do ig=1,nffc
     is=ig+(il-1)*nffc
     lod(is)=lodfc(il)
   enddo
 enddo
!!!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 allocate(st_leafmtx(nblall))
 nlfalt=sqrt(real(nblall)); st_leafmtxp%nlfalt=nlfalt
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'Number of MPI_Blocks=',nblall,'; sqrt(nblall)=',nlfalt
 ndpth=0; lnmtx(1:4)=0
 if(st_ctl%param(23)==1)then
   call HACApK_count_LH(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
   lnmtx(4)=nblall
 else
   call HACApK_count_blrleaf(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 endif

 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'No. of nsmtx',lnmtx(1:4)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'   1:Rk-matrix 2: dense-mat 3:H-matrix 4:MPI_Block'
 st_leafmtxp%nlfkt=lnmtx(1)
 nlfall=lnmtx(1)+lnmtx(2); 
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'nlf global=',nlfall
 if(nlfall<nthr)then
   print*,'Error; HACApK_generate_frame_blrleaf; # of leaves must be larger than # of threads.'
   stop
 endif
 
 if(st_ctl%param(23)==1)then
   ndpth=0; lnmtx(1:4)=0
   call HACApK_generate_LH(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nblall,ndpth)
 else
   nblall=0; ndpth=0; lnmtx(1:4)=0
   call HACApK_generate_blrleaf(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nblall,ndpth)
 endif
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_blrleaf; HACApK_generate_leafmtx end'
 call HACApK_sort_leafmtx(st_leafmtx,nblall)
 do ip=1,nblall
   if(st_leafmtx(ip)%ltmtx==4) call HACApK_sort_leafmtx(st_leafmtx(ip)%st_lf,st_leafmtx(ip)%nlf)
 enddo
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_blrleaf; HACApK_sort_leafmtx end' 
 call HACApK_free_st_clt(st_clt)
  
!!!!!!!!!!!!!!!!!! start MPI load balance  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 npgt=nrank/npgl
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) ' npgl=',npgl,' npgt=',npgt
 if(st_ctl%param(1)>1) write(mpilog,1000) ' npgl=',npgl,' npgt=',npgt
 if(npgt>nlfalt .or. npgl>nlfalt)then
   call MPI_Barrier( icomm, ierr )
   if(mpinr==0) write(6,*) 'Error: HACApK_generate_frame_blrleaf; Too few blocks compared with #MPI !!!'; goto 9999
 endif
 if(npgt*npgl/=nrank)then
   call MPI_Barrier( icomm, ierr )
   if(mpinr==0) write(6,*) 'Error: HACApK_generate_frame_blrleaf; Invalid processor grid!!!'; goto 9999
 endif
 
! Split MPI communicator
  ikey=0; iclr=mpinr/npgt
  call MPI_COMM_SPLIT(icomm, iclr, ikey, icommn, ierr); st_ctl%lpmd(31)=icommn
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; MPI_COMM_SPLIT failed !!!'; goto 9999
  endif
  call MPI_Comm_size ( icommn, nrank, ierr ); st_ctl%lpmd(32)=nrank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_size failed !!!'; goto 9999
  endif
  call MPI_Comm_rank ( icommn, irank, ierr ); st_ctl%lpmd(33)=irank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_rank failed !!!'; goto 9999
  endif
  ikey=0; iclr=mod(mpinr,npgt)
  call MPI_COMM_SPLIT(icomm, iclr, ikey, icommn, ierr); st_ctl%lpmd(35)=icommn
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; MPI_COMM_SPLIT failed !!!'; goto 9999
  endif
  call MPI_Comm_size ( icommn, nrank, ierr ); st_ctl%lpmd(36)=nrank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_size failed !!!'; goto 9999
  endif
  call MPI_Comm_rank ( icommn, irank, ierr ); st_ctl%lpmd(37)=irank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrleaf; MPI_Comm_rank failed !!!'; goto 9999
  endif

 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,'; irank_t=',st_ctl%lpmd(33),'; irank_l=',st_ctl%lpmd(37)
 if(st_ctl%param(1)>1) write(mpilog,1000) 'irank_t=',st_ctl%lpmd(33),'; nrank_t=',st_ctl%lpmd(32)
 if(st_ctl%param(1)>1) write(mpilog,1000) 'irank_l=',st_ctl%lpmd(37),'; nrank_l=',st_ctl%lpmd(36)

! stop
 
 nlft=nlfalt/npgt; nlfth=mod(nlfalt,npgt); mpinrth=mod(mpinr,npgt)
 if(mpinrth<nlfth) nlft=nlft+1
 nlfl=nlfalt/npgl; nlflh=mod(nlfalt,npgl); mpinrlh=mpinr/npgt
 if(mpinrlh<nlflh) nlfl=nlfl+1
 nbl=nlfl*nlft; st_leafmtxp%nbl=nbl; st_leafmtxp%nlfl=nlfl; st_leafmtxp%nlft=nlft
 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,'; nbl=',nbl,'; nblall=',nblall
 if(st_ctl%param(1)>1) write(mpilog,1000) 'No. of blocks; nbl=',nbl,'; row=',nlfl,'; column=',nlft,'; global nbl=',nblall
 
 
   call MPI_Barrier( icomm, ierr )
! stop
 
 nrank_t=st_ctl%lpmd(32); nrank_l=st_ctl%lpmd(36)
 
 allocate(st_leafmtxp%lbl2t(0:npgl-1)); st_leafmtxp%lbl2t(:)=0
 irank_t=st_ctl%lpmd(33); irank_l=st_ctl%lpmd(37)
 do in=0,nlfalt-1
   inml=mod(in,npgl); inmt=mod(in,npgt)
   if(inmt==irank_t) st_leafmtxp%lbl2t(inml)=1
 enddo
 if(st_ctl%param(1)>1) write(mpilog,*) 'st_leafmtxp%lbl2t'
 if(st_ctl%param(1)>1) write(mpilog,*) st_leafmtxp%lbl2t(:) 
 
 allocate(st_leafmtxp%lbstrtl(nlfalt+1),st_leafmtxp%lbstrtt(nlfalt+1))
 allocate(st_leafmtxp%lbndl(nlfalt),st_leafmtxp%lbndt(nlfalt))
 allocate(st_leafmtxp%lbndlfs(0:nrank_l-1),st_leafmtxp%lbndtfs(0:nrank_t-1)); st_leafmtxp%lbndlfs=0; st_leafmtxp%lbndtfs=0
 do il=0,nlfalt-1; is=nlfalt*il+1
   ilh=mod(il,npgl)
   st_leafmtxp%lbndlfs(ilh)=st_leafmtxp%lbndlfs(ilh)+st_leafmtx(is)%ndl
   st_leafmtxp%lbstrtl(il+1)=st_leafmtx(is)%nstrtl
   st_leafmtxp%lbndl(il+1)=st_leafmtx(is)%ndl
 enddo 
 st_leafmtxp%lbstrtl(nlfalt+1)=nd+1
 do it=0,nlfalt-1; is=it+1
   ith=mod(it,npgt)
   st_leafmtxp%lbndtfs(ith)=st_leafmtxp%lbndtfs(ith)+st_leafmtx(is)%ndt
   st_leafmtxp%lbstrtt(it+1)=st_leafmtx(is)%nstrtt
   st_leafmtxp%lbndt(it+1)=st_leafmtx(is)%ndt
 enddo 
 st_leafmtxp%lbstrtt(nlfalt+1)=nd+1
 if(.false.) then
! if(mpinr==0) then
   print*,'lbstrtl='
   print*,st_leafmtxp%lbstrtl
   print*,'lbstrtt='
   print*,st_leafmtxp%lbstrtt
   print*,'lbndlfs='
   print*,st_leafmtxp%lbndlfs
   print*,'lbndtfs='
   print*,st_leafmtxp%lbndtfs
 endif
 
  allocate(st_leafmtx_lcl(nbl),st_leafmtxp%lnlfl2g(nlft,nlfl)); 
  ip=0; nlf=0
  do il=0,nlfalt-1; do it=0,nlfalt-1; is=it+nlfalt*il+1
    ilh=mod(il,npgl); ith=mod(it,npgt)
    ipgclr=ith+ilh*npgt
    if(ipgclr==mpinr)then
      ilf=ip/nlft+1; itf=mod(ip,nlft)+1; ip=ip+1; 
      st_leafmtxp%lnlfl2g(itf,ilf)=is
      st_leafmtx_lcl(ip)=st_leafmtx(is)
      if(st_leafmtx(is)%ltmtx<=2)then
        nlf=nlf+1
        st_leafmtx(is)%lttcl=ilf; st_leafmtx(is)%lttct=itf
      else
        nlf=nlf+st_leafmtx(is)%nlf
        st_leafmtx(is)%st_lf(1:st_leafmtx(is)%nlf)%lttcl=ilf
        st_leafmtx(is)%st_lf(1:st_leafmtx(is)%nlf)%lttct=itf
      endif
    endif
  enddo; enddo
  
   ndlfs=0
   do ilf=1,nlfl
     ip=(ilf-1)*nlft+1
     ndlfs=ndlfs+st_leafmtx_lcl(ip)%ndl
   enddo
   st_leafmtxp%ndlfs=ndlfs

   ndtfs=0
   do ip=1,nlft
     ndtfs=ndtfs+st_leafmtx_lcl(ip)%ndt
   enddo
   st_leafmtxp%ndtfs=ndtfs
  
!  print*,'mpinr=',mpinr,'; nlf=',nlf
  if(nlf<1) then
    print*, 'Error: sub HACApK_generate_frame_blrleaf; nlf<1 !!!; mpinr=',mpinr; goto 9999
  endif
  st_leafmtxp%nlf=nlf
  allocate(st_leafmtxp%st_lf(nlf)); 
  ip=0; ndlfs=0; ndtfs=0
  do il=0,nlfalt-1; do it=0,nlfalt-1; is=it+nlfalt*il+1
    ilh=mod(il,npgl); ith=mod(it,npgt)
    ipgclr=ith+ilh*npgt
    if(ipgclr==mpinr)then
      if(st_leafmtx(is)%ltmtx<=2)then
        ip=ip+1; 
        st_leafmtxp%st_lf(ip)=st_leafmtx(is)
      else
        isnlf=st_leafmtx(is)%nlf
        st_leafmtxp%st_lf(ip+1:ip+isnlf)=st_leafmtx(is)%st_lf(1:isnlf)
        ip=ip+isnlf
!!!        if(st_leafmtx(is)%nstrtt==1) ndlfs=ndlfs+st_leafmtx(is)%ndl
!!!        if(st_leafmtx(is)%nstrtl==1) ndtfs=ndtfs+st_leafmtx(is)%ndt
        
      endif
    endif
  enddo; enddo

  deallocate(st_leafmtx,st_leafmtx_lcl)
  
!!! print*,'mpinr=',mpinr,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
!!! call MPI_Barrier( icomm, ierr )
! stop
  
 if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
 if(st_ctl%param(1)>1) write(mpilog,1000) 'Vector sizes; nd=',nd,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
  
 if(st_ctl%param(1)>1) then
   write(mpilog,*) 'lnlfl2g='
   do il=1,nlfl
     write(mpilog,'(10(i9))') st_leafmtxp%lnlfl2g(:,il)
   enddo
 endif

   call MPI_Barrier( icomm, ierr )
    
  lnmtx(:)=0; mem8=0; ktp=param(62)
  do ip=1,nlf
    ltmtx=st_leafmtxp%st_lf(ip)%ltmtx; ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
    if(ltmtx==1)then
      lnmtx(1)=lnmtx(1)+1; mem8=mem8+(ndt+ndl)*ktp
    else
      lnmtx(2)=lnmtx(2)+1; mem8=mem8+ns
    endif
  enddo
 if(st_ctl%param(1)>1)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:2)
  call HACApK_setcutthread(lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
 return
 9999 continue
 stop
 end subroutine HACApK_generate_frame_blrleaf

!***HACApK_generate_frame_blrmtx
 subroutine HACApK_generate_frame_blrmtx(st_leafmtxp,st_bemv,st_ctl,gmid,lnmtx,nofc,nffc,ndim)
 include 'mpif.h'
 type(st_HACApK_cluster) :: st_clt
 type(st_HACApK_lcontrol) :: st_ctl
 type(st_HACApK_calc_entry) :: st_bemv
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_leafmtx),dimension(:), allocatable :: st_leafmtx
 real*8 :: gmid(nofc,ndim)
 integer*8 :: mem8,nlfall
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
   print*,'HACApK_generate_frame_blrmtx: ierr=',ierr
 endif
 do il=1,nofc
   lodfc(il)=il
 enddo
 npgl=param(41); 
 ! if(npgl==0) npgl=sqrt(real(nrank))
 ! npgt=nrank/npgl
 if(npgl==0) then
   npgl = sqrt(real(nrank))
   npgt = nrank/npgl
   do while(npgl*npgt .ne. nrank)
     npgl = npgl - 1
     npgt = nrank/npgl
   enddo
   if (mpinr==0) then
     write(*,*) 'nrank=',nrank,'(',npgl,'x',npgt,')'
   endif
 elseif(npgl<0) then
   npgl=nrank
 endif
! if(param(42)==0) param(42)=sqrt(real(nd))
 if(param(42)==0) param(42)=nd/param(43)/npgl
 if(param(23)==1) param(53)=2
 
!!!!!!!!!!!!!!!!!! start clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nsrt=1; ndf=nofc; nclst=0; ndpth=0; ndscd=0
 if(st_ctl%param(23)==1)then
   call HACApK_generate_tctree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst,nblall)
 else
   param(21)=param(42)
   call HACApK_generate_cbitree(st_clt,gmid,param,lpmd,lodfc,ndpth,ndscd,nsrt,ndf,nofc,ndim,nclst)
 endif
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'No. of cluster=',nclst
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'No. of cluster=',nclst

 call HACApK_bndbox(st_clt,gmid,lodfc,nofc)
 do il=1,nofc
   do ig=1,nffc
     is=ig+(il-1)*nffc
     lod(is)=lodfc(il)
   enddo
 enddo
!!!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ndpth=0; lnmtx(1:3)=0
 call HACApK_count_blr(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'No. of nsmtx',lnmtx(1:3)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'   1:Rk-matrix 2: dense-mat 3:H-matrix'
 st_leafmtxp%nlfkt=lnmtx(1)
 nlfall=lnmtx(1)+lnmtx(2); nlfalt=sqrt(real(nlfall)); st_leafmtxp%nlfalt=nlfalt
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'nlf global=',nlfall,'; sqrt(nlf)=',nlfalt
 allocate(st_leafmtx(nlfall))
 if(nlfall<nthr)then
   print*,'Error; HACApK_generate_frame_leafmtx; # of leaves must be larger than # of threads.'
   stop
 endif
 
 nlf=0; ndpth=0
 call HACApK_generate_blr(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_generate_leafmtx end'
 call HACApK_sort_leafmtx(st_leafmtx,nlf)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_qsort_row_leafmtx end'
 call HACApK_free_st_clt(st_clt)
 
!!!!!!!!!!!!!!!!!! start MPI load balance  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if(.false.)then
! if(mpinr==0)then
   do ip=1,nlf
     write(*,1000) 'ip=',ip,'; nstrtl=',st_leafmtx(ip)%nstrtl,'; nstrtt=',st_leafmtx(ip)%nstrtt
   enddo
 endif

 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) ' npgl=',npgl,' npgt=',npgt
 if(st_ctl%param(1)>1) write(mpilog,1000) ' npgl=',npgl,' npgt=',npgt
 if(npgt>nlfalt .or. npgl>nlfalt)then
   call MPI_Barrier( icomm, ierr )
   if(mpinr==0) write(6,*) 'Error: HACApK_generate_frame_blrmtx; Too few blocks compared with #MPI !!!'; goto 9999
 endif
 if(npgt*npgl/=nrank)then
   call MPI_Barrier( icomm, ierr )
   if(mpinr==0) write(6,*) 'Error: HACApK_generate_frame_blrmtx; Invalid processor grid!!!'; goto 9999
 endif
 
! Split MPI communicator
  ikey=0; iclr=mpinr/npgt
  call MPI_COMM_SPLIT(icomm, iclr, ikey, icommn, ierr); st_ctl%lpmd(31)=icommn
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrmtx; MPI_COMM_SPLIT failed !!!'; goto 9999
  endif
  call MPI_Comm_size ( icommn, nrank, ierr ); st_ctl%lpmd(32)=nrank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrmtx; MPI_Comm_size failed !!!'; goto 9999
  endif
  call MPI_Comm_rank ( icommn, irank, ierr ); st_ctl%lpmd(33)=irank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrmtx; MPI_Comm_rank failed !!!'; goto 9999
  endif
  ikey=0; iclr=mod(mpinr,npgt)
  call MPI_COMM_SPLIT(icomm, iclr, ikey, icommn, ierr); st_ctl%lpmd(35)=icommn
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrmtx; MPI_COMM_SPLIT failed !!!'; goto 9999
  endif
  call MPI_Comm_size ( icommn, nrank, ierr ); st_ctl%lpmd(36)=nrank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrmtx; MPI_Comm_size failed !!!'; goto 9999
  endif
  call MPI_Comm_rank ( icommn, irank, ierr ); st_ctl%lpmd(37)=irank
  if(ierr.ne.0) then
    if(mpinr==0) print*, 'Error: sub HACApK_generate_frame_blrmtx; MPI_Comm_rank failed !!!'; goto 9999
  endif

 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,'; irank_t=',st_ctl%lpmd(33),'; irank_l=',st_ctl%lpmd(37)
 if(st_ctl%param(1)>1) write(mpilog,1000) 'irank_t=',st_ctl%lpmd(33),'; nrank_t=',st_ctl%lpmd(32)
 if(st_ctl%param(1)>1) write(mpilog,1000) 'irank_l=',st_ctl%lpmd(37),'; nrank_l=',st_ctl%lpmd(36)

! stop
 
 nlft=nlfalt/npgt; nlfth=mod(nlfalt,npgt); mpinrth=mod(mpinr,npgt)
 if(mpinrth<nlfth) nlft=nlft+1
 nlfl=nlfalt/npgl; nlflh=mod(nlfalt,npgl); mpinrlh=mpinr/npgt
 if(mpinrlh<nlflh) nlfl=nlfl+1
 nlf=nlfl*nlft; st_leafmtxp%nlf=nlf; st_leafmtxp%nlfl=nlfl; st_leafmtxp%nlft=nlft
 if(st_ctl%param(1)>1) write(*,1000) 'irank=',mpinr,'; nlf=',nlf,'; nlfall=',nlfall
 if(st_ctl%param(1)>1) write(mpilog,1000) 'No. of blocks; nlf=',nlf,'; row=',nlfl,'; column=',nlft,'; global nlf=',nlfall
 
   call MPI_Barrier( icomm, ierr )
! stop
 
 nrank_t=st_ctl%lpmd(32); nrank_l=st_ctl%lpmd(36)
 
 allocate(st_leafmtxp%lbl2t(0:npgl-1)); st_leafmtxp%lbl2t(:)=0
 irank_t=st_ctl%lpmd(33); irank_l=st_ctl%lpmd(37)
 do in=0,nlfalt-1
   inml=mod(in,npgl); inmt=mod(in,npgt)
   if(inmt==irank_t) st_leafmtxp%lbl2t(inml)=1
 enddo
 if(st_ctl%param(1)>1) write(mpilog,*) 'st_leafmtxp%lbl2t'
 if(st_ctl%param(1)>1) write(mpilog,*) st_leafmtxp%lbl2t(:) 
 
 allocate(st_leafmtxp%lbstrtl(nlfalt+1),st_leafmtxp%lbstrtt(nlfalt+1))
 allocate(st_leafmtxp%lbndl(nlfalt),st_leafmtxp%lbndt(nlfalt))
 allocate(st_leafmtxp%lbndlfs(0:nrank_l-1),st_leafmtxp%lbndtfs(0:nrank_t-1)); st_leafmtxp%lbndlfs=0; st_leafmtxp%lbndtfs=0
 do il=0,nlfalt-1; is=nlfalt*il+1
   ilh=mod(il,npgl)
   st_leafmtxp%lbndlfs(ilh)=st_leafmtxp%lbndlfs(ilh)+st_leafmtx(is)%ndl
   st_leafmtxp%lbstrtl(il+1)=st_leafmtx(is)%nstrtl
   st_leafmtxp%lbndl(il+1)=st_leafmtx(is)%ndl
 enddo 
 st_leafmtxp%lbstrtl(nlfalt+1)=nd+1
 do it=0,nlfalt-1; is=it+1
   ith=mod(it,npgt)
   st_leafmtxp%lbndtfs(ith)=st_leafmtxp%lbndtfs(ith)+st_leafmtx(is)%ndt
   st_leafmtxp%lbstrtt(it+1)=st_leafmtx(is)%nstrtt
   st_leafmtxp%lbndt(it+1)=st_leafmtx(is)%ndt
 enddo 
 st_leafmtxp%lbstrtt(nlfalt+1)=nd+1
 if(.false.) then
! if(mpinr==0) then
   print*,'lbstrtl='
   print*,st_leafmtxp%lbstrtl
   print*,'lbstrtt='
   print*,st_leafmtxp%lbstrtt
   print*,'lbndlfs='
   print*,st_leafmtxp%lbndlfs
   print*,'lbndtfs='
   print*,st_leafmtxp%lbndtfs
 endif
 
  allocate(st_leafmtxp%st_lf(nlf),st_leafmtxp%lnlfl2g(nlft,nlfl)); 
  ip=0
  do il=0,nlfalt-1; do it=0,nlfalt-1; is=it+nlfalt*il+1
    ilh=mod(il,npgl); ith=mod(it,npgt)
    ipgclr=ith+ilh*npgt
    if(ipgclr==mpinr)then
! if(mpinr==0) print*,'is=',is,'; mpinr=',mpinr
!  print*,'is=',is,'; mpinr=',mpinr
      ilf=ip/nlft+1; itf=mod(ip,nlft)+1; ip=ip+1; 
!      st_leafmtxp%lnlfl2g(ip)=is
      st_leafmtxp%lnlfl2g(itf,ilf)=is
      st_leafmtxp%st_lf(ip)=st_leafmtx(is)
    endif
  enddo; enddo
  deallocate(st_leafmtx)

   ndlfs=0
   do ilf=1,nlfl
     ip=(ilf-1)*nlft+1
     ndlfs=ndlfs+st_leafmtxp%st_lf(ip)%ndl
   enddo
   st_leafmtxp%ndlfs=ndlfs

   ndtfs=0
   do ip=1,nlft
     ndtfs=ndtfs+st_leafmtxp%st_lf(ip)%ndt
   enddo
   st_leafmtxp%ndtfs=ndtfs
  
!!! print*,'mpinr=',mpinr,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
!!! call MPI_Barrier( icomm, ierr )
! stop
  
 if(st_ctl%param(1)>1) write(*,*) 'irank=',mpinr,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
 if(st_ctl%param(1)>1) write(mpilog,1000) 'Vector sizes; nd=',nd,'; ndlfs=',st_leafmtxp%ndlfs,'; ndtfs=',st_leafmtxp%ndtfs
  
 if(st_ctl%param(1)>1) then
   write(mpilog,*) 'lnlfl2g='
   do il=1,nlfl
     write(mpilog,'(10(i9))') st_leafmtxp%lnlfl2g(:,il)
   enddo
 endif

   call MPI_Barrier( icomm, ierr )
  
  lnmtx(:)=0; mem8=0; ktp=param(62)
  do ip=1,nlf
    ltmtx=st_leafmtxp%st_lf(ip)%ltmtx; ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
    if(ltmtx==1)then
      lnmtx(1)=lnmtx(1)+1; mem8=mem8+(ndt+ndl)*ktp
    else
      lnmtx(2)=lnmtx(2)+1; mem8=mem8+ns
    endif
  enddo
 if(st_ctl%param(1)>1)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:2)
  call HACApK_setcutthread(lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
 return
 9999 continue
 stop
 end subroutine HACApK_generate_frame_blrmtx

!***HACApK_generate_frame_leafmtx_
 subroutine HACApK_generate_frame_leafmtx_(st_leafmtxp,st_bemv,st_ctl,gmid,lnmtx,nofc,nffc,ndim)
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
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'No. of cluster=',nclst

 call HACApK_bndbox(st_clt,gmid,lodfc,nofc)
!!!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ndpth=0; lnmtx(1:3)=0
 if(st_ctl%param(54)==0)then
   call HACApK_count_lntmx(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
!   call HACApK_count_olntmx(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc)
 elseif(st_ctl%param(54)==2)then
   call HACApK_count_zlntmx(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 elseif(st_ctl%param(54)==1)then
   call HACApK_count_blr(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 else
   print*,'Error; HACApK_generate_frame_leafmtx; Set param(54)=0,1 or 2.'
   stop
 endif
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'No. of nsmtx',lnmtx(1:3)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'   1:Rk-matrix 2: dense-mat 3:H-matrix'
 st_leafmtxp%nlfkt=lnmtx(1)
 nlf=lnmtx(1)+lnmtx(2)
 allocate(st_leafmtx(nlf))
 st_leafmtxp%nlf=nlf; nlf=0
 if(st_leafmtxp%nlf<nthr)then
   print*,'Error; HACApK_generate_frame_leafmtx; # of leaves must be larger than # of threads.'
   stop
 endif
 
 ndpth=0
 if(st_ctl%param(54)==0)then
   call HACApK_generate_leafmtx(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
!   call HACApK_generate_oleafmtx(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf)
 elseif(st_ctl%param(54)==2)then
   call HACApK_generate_zleafmtx(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 elseif(st_ctl%param(54)==1)then
   call HACApK_generate_blr(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 endif
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_generate_leafmtx end'
 call HACApK_sort_leafmtx(st_leafmtx,nbl)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'HACApK_generate_frame_leafmtx; HACApK_qsort_row_leafmtx end' 
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
 if(st_ctl%param(1)>1)  write(mpilog,*) 'No. of nsmtx',lnmtx(1:2)
  deallocate(st_leafmtx)
  call HACApK_setcutthread(lthr,st_leafmtxp,st_ctl,mem8,nthr,ktp)
 9999 continue
! stop
 end subroutine 
 
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
 if(st_ctl%param(1)>1)  write(mpilog,1000) 'No. of cluster=',nclst

 call HACApK_bndbox(st_clt,gmid,lodfc,nofc)
!!!!!!!!!!!!!!!!!! end clustering !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!! start construction of H-matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ndpth=0; lnmtx(1:3)=0
 call HACApK_count_olntmx(st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'No. of nsmtx',lnmtx(1:3)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'   1:Rk-matrix 2: dense-mat 3:H-matrix'
 st_leafmtxp%nlfkt=lnmtx(1)
 nlf=lnmtx(1)+lnmtx(2)
 allocate(st_leafmtx(nlf))
 st_leafmtxp%nlf=nlf; nlf=0
 
 call HACApK_generate_oleafmtx(st_leafmtx,st_clt,st_clt,param,lpmd,lnmtx,nofc,nffc,nlf)
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
! if(st_ctl%param(1)>1) print*,'HACApK_setcutthread; nlf=',nlf,' mem8=',mem8,' nthr=',nthr
 lthr(0)=1; lthr(nthr)=nlf+1
 imem=0; ith=1; kt=ktp
 do il=1,nlf
   ltmtx=st_leafmtxp%st_lf(il)%ltmtx
   ndl=st_leafmtxp%st_lf(il)%ndl; ndt=st_leafmtxp%st_lf(il)%ndt
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
! if(st_ctl%param(1)>1) print*,'HACApK_setcutthread; lthr=',lthr(0:nthr)
 end subroutine HACApK_setcutthread

!***HACApK_accuracy_leafmtx_body
 subroutine HACApK_accuracy_leafmtx_body(zhnrm,zanrm,st_leafmtxp,st_bemv,lodl,lodt,lpmd,nofc,nffc)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_calc_entry) :: st_bemv
 integer*4 :: lodl(nofc*nffc),lodt(nofc*nffc),lpmd(*)
 1000 format(5(a,i12)/)
 2000 format(5(a,1pe15.7)/)
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
! write(mpilog,1000) 'sub HACApK_accuracy_leafmtx_body; nlf=',st_leafmtxp%nlf,'; mpinr',mpinr
 do ip=1,st_leafmtxp%nlf
! write(mpilog,1000) 'sub HACApK_accuracy_leafmtx_body; ip=',ip,'; mpinr',mpinr
   ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
! write(mpilog,1000) 'sub HACApK_accuracy_leafmtx_body; ip=',ip,'; ltmtx',st_leafmtxp%st_lf(ip)%ltmtx
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
! write(mpilog,1000) 'sub HACApK_accuracy_leafmtx_body; ip=',ip,'; end'
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
 3000 format(a,i12,5(a,1pe15.8)/)
 
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 if(st_ctl%param(1)>1 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx start'
 zhnrm=0.0d0; zanrm=0.0d0; ktmax=0
 zhnrm1=0.0d0; zanrm1=0.0d0
 call HACApK_accuracy_leafmtx_body(zhnrm1,zanrm1,st_leafmtxp,st_bemv,lodl,lodt,lpmd,nofc,nffc)
! if(st_ctl%param(1)>1) write(*,3000) 'irank=',mpinr,'; zhnrm1=',zhnrm1,'; zanrm1=',zanrm1
 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( zanrm1, zanrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 call MPI_reduce( zhnrm1, zhnrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 zanrm=dsqrt(zanrm); zhnrm=dsqrt(zhnrm)
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx; zanrm=',zanrm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx; zhnrm=',zhnrm
 if(st_ctl%param(1)>0 .and. mpinr==0) print*,'sub HACApK_accuracy_leafmtx; zhnrm/zanrm=',zhnrm/zanrm
 if(st_ctl%param(1)>1)  write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zanrm1=',dsqrt(zanrm1)
 if(st_ctl%param(1)>1)  write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm1=',dsqrt(zhnrm1)
 if(st_ctl%param(1)>1)  write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm1/zanrm1=',dsqrt(zhnrm1)/dsqrt(zanrm1)
  if(mpinr==0 .and. st_ctl%param(1)>1) write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zanrm=',zanrm
  if(mpinr==0 .and. st_ctl%param(1)>1) write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm=',zhnrm
  if(mpinr==0 .and. st_ctl%param(1)>1) write(mpilog,2000) 'sub HACApK_accuracy_leafmtx; zhnrm/zanrm=',zhnrm/zanrm

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
!  print*,'k=',k,' zeps=',zeps,' ist=',ist
   pcol => zaa(:,k); prow => zab(:,k)
    call HACApK_calc_vec(zab, zaa, ndt, k-1, ist, prow, nstrtl, nstrtt,lod, st_bemv, lcol_msk,0) 
!   print*, 'prow=',prow
   call HACApK_maxabsvallocm_d(prow,row_maxval,jst,ndt,lcol_msk)
!   print*,'jst=',jst,' row_maxval=',row_maxval
! write(6,1000) 'ist=',ist,'  jst=',jst
     zdltinv=1.0d0/prow(jst); prow(:)=prow(:)*zdltinv
    call HACApK_calc_vec(zaa, zab, ndl, k-1, jst, pcol, nstrtl, nstrtt,lod, st_bemv, lrow_msk,1)
   lrow_msk(ist)=1
   call HACApK_maxabsvallocm_d(pcol,col_maxval,istn,ndl,lrow_msk)
!   print*, 'pcol=',pcol
!   print *,'zeps=',zeps
!   print*,'istn=',istn,' col_maxval=',col_maxval
   lcol_msk(jst)=1
   ist=istn; nrow_done=nrow_done+1; ncol_done=ncol_done+1
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

!***HACApK_acaplus
 integer function HACApK_acaplus(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lod,st_bemv,kmax,eps,znrmmat,pACA_EPS)
 type(st_HACApK_calc_entry) :: st_bemv
 real*8 :: param(:)
 real*8,target :: zaa(ndl,kmax),zab(ndt,kmax)
 integer*4 :: lod(:)
 integer*4,dimension(:), allocatable :: lrow_msk,lcol_msk
 real*8,dimension(:), allocatable :: pa_ref, pb_ref
 real*8,pointer :: prow(:),pcol(:)
 1000 format(5(a,i12)/)
 2000 format(5(a,1pe15.8)/)

  za_ACA_EPS=1.0e-30
! write(6,1000) 'nstrtl=',nstrtl,' nstrtt=',nstrtt,' ndl=',ndl,' ndt=',ndt
 znrm=znrmmat*sqrt(real(ndl)*real(ndt))
 if(param(61)==2 .or. param(61)==1) ACA_EPS=pACA_EPS
 if(param(61)==3) ACA_EPS=pACA_EPS*znrm
 
 HACApK_acaplus=0; ntries = max(ndl,ndt)+1; ntries_row = 6; ntries_col = 6;
 allocate(lrow_msk(ndl),lcol_msk(ndt))
 k = 0; lrow_msk(:)=0; lcol_msk(:)=0
 
  j_ref=1; allocate(pa_ref(ndl)) ! arbitrary j_ref
  call HACApK_calc_vec(zaa, zab, ndl, k, j_ref, pa_ref, nstrtl, nstrtt,lod, st_bemv, lrow_msk, 1)
!!!  print*,'pa_ref=',pa_ref
  colnorm = HACApK_unrm_d(ndl,pa_ref)

  call HACApK_minabsvalloc_d(pa_ref,rownorm,i_ref,ndl) ! determine i_ref:=argmin ||pa_ref(1:ndl)||
!!!    print*,'i_ref=',i_ref
  allocate(pb_ref(ndt))
  call HACApK_calc_vec(zab, zaa, ndt, k, i_ref, pb_ref, nstrtl, nstrtt,lod, st_bemv, lcol_msk,0)
!!!  print*,'pb_ref=',pb_ref
  rownorm=HACApK_unrm_d(ndt,pb_ref)

  apxnorm = 0.0; lstop_aca = 0;
    
  do while((k<kmax) .and. (ntries_row>0 .or. ntries_col>0) .and. (ntries>0))
    ntries=ntries-1
    pcol => zaa(:,k+1); prow => zab(:,k+1)
    col_maxval = 0.0; call HACApK_maxabsvalloc_d(pa_ref,col_maxval,i,ndl)
    row_maxval = 0.0; call HACApK_maxabsvalloc_d(pb_ref,row_maxval,j,ndt)
    
!!!    write(6,1000) 'i=',i,' i_ref=',i_ref,' j=',j,' j_ref=',j_ref
    
    if(row_maxval>col_maxval)then
      if(j/=j_ref)then; call HACApK_calc_vec(zaa, zab, ndl, k, j, pcol, nstrtl, nstrtt,lod, st_bemv, lrow_msk, 1)
      else; pcol(:)=pa_ref(:)
      endif
      call HACApK_maxabsvalloc_d(pcol,col_maxval,i,ndl)
            
      if(col_maxval < ACA_EPS .and. k>=param(64))then; lstop_aca = 1; 
!         print*,'2***************lstop_aca==1***********************2'
      else
        call HACApK_calc_vec(zab, zaa, ndt, k, i, prow, nstrtl, nstrtt,lod, st_bemv, lcol_msk,0)
        if(abs(pcol(i))>1.0e-20) then
          zinvmax=1.0/pcol(i)
        else
          k=max(k-1,0); exit
        endif
!        if(isnan(zinvmax))then
!          print*,'1.0/pcol(i)=NaN',' k=',k
!          exit
!          stop
!        endif
        pcol(:)=pcol(:)*zinvmax
      endif
    else
      if(i/=i_ref)then; call HACApK_calc_vec(zab, zaa, ndt, k, i, prow, nstrtl, nstrtt,lod, st_bemv, lcol_msk,0)
      else;  prow(:)=pb_ref(:)
      endif
      call HACApK_maxabsvalloc_d(prow,row_maxval,j,ndt)
      
      if(row_maxval < ACA_EPS .and. k>=param(64))then; lstop_aca = 1
!         print*,'3***************lstop_aca==1***********************3'
      else
        call HACApK_calc_vec(zaa, zab, ndl, k, j, pcol, nstrtl, nstrtt,lod, st_bemv, lrow_msk, 1)
        if(abs(prow(j))>1.0e-20) then
          zinvmax=1.0/prow(j)
        else
          k=max(k-1,0); exit
        endif
!        if(isnan(zinvmax))then
!          print*,'1.0/prow(j)=NaN',' k=',k
!          exit
!          stop
!        endif
        prow(:)=prow(:)*zinvmax
      endif
    endif
    lrow_msk(i) = 1; lcol_msk(j) = 1
!!!    write(6,1000) 'i=',i,' i_ref=',i_ref,' j=',j,' j_ref=',j_ref
    
    if(i/=i_ref)then
      zinvmax = -pcol(i_ref)
      pb_ref(:)=pb_ref(:)+prow(:)*zinvmax
      rownorm = HACApK_unrm_d(ndt,pb_ref)
    endif
    if(i==i_ref .or. rownorm<ACA_EPS)then
      if(i==i_ref) ntries_row=ntries_row+1
      if(ntries_row>0)then
        rownorm = 0.0; i=i_ref
!        print*,'lrow_msk',lrow_msk
        do while(i/=mod((i_ref+ndl-2),ndl)+1 .and. rownorm<za_ACA_EPS .and. ntries_row>0)
!          print*,'i=',i,' ii=',mod((i_ref+ndl-2),ndl)+1
          if(lrow_msk(i)==0)then
!            write(6,1000) 'i=',i
            call HACApK_calc_vec(zab, zaa, ndt, k+1, i, pb_ref, nstrtl, nstrtt,lod, st_bemv, lcol_msk,0)
            rownorm = HACApK_unrm_d(ndt,pb_ref)
            if(rownorm<ACA_EPS) lrow_msk(i) = 1
            ntries_row=ntries_row-1
          else
            rownorm = 0.0;
          endif
          i=mod(i,ndl)+1
        enddo
        i_ref=mod((i+ndl-2),ndl)+1
      endif
    endif
!!!    print*,'i_ref=',i_ref
    
    if(j/=j_ref)then
      zinvmax = -prow(j_ref)
      pa_ref(:)=pa_ref(:)+pcol(:)*zinvmax
      colnorm = HACApK_unrm_d(ndl,pa_ref)
    endif
    if(j==j_ref .or. colnorm<ACA_EPS)then
      if(j==j_ref) ntries_col=ntries_col+1
      if(ntries_col>0)then
        colnorm = 0.0; j=j_ref
!        print*,'lcol_msk',lcol_msk
        do while(j/=mod((j_ref+ndt-2),ndt)+1 .and. colnorm<za_ACA_EPS .and. ntries_col>0)
          if(lcol_msk(j)==0)then
            call HACApK_calc_vec(zaa, zab, ndl, k+1, j, pa_ref, nstrtl, nstrtt,lod, st_bemv, lrow_msk, 1)
            colnorm = HACApK_unrm_d(ndl,pa_ref)
            if(colnorm<ACA_EPS) lcol_msk(j)=1
            ntries_col=ntries_col-1
          else
            colnorm = 0.0
          endif
          j=mod(j,ndt)+1
        enddo
        j_ref=mod((j+ndt-2),ndt)+1
      endif
    endif

!    write(6,2000) 'colnorm=',colnorm,' rownorm=',rownorm
    if(colnorm<ACA_EPS .and. rownorm<ACA_EPS .and. k>=param(64))then
      lstop_aca=1; k=k+1;
!       print*,'1***************lstop_aca==1***********************1'
    endif

    if(lstop_aca==0)then
      blknorm = (HACApK_unrm_d(ndl,pcol)*HACApK_unrm_d(ndt,prow))
      if(k == 0)then
        if(param(61)==1)then
          apxnorm = blknorm
        elseif(param(61)==2 .or. param(61)==3)then
          apxnorm =znrm
        else
!$omp critical
          print*,'ERROR!:: invalid param(61)=',param(61)
!$omp end critical
          stop
        endif
      else
        if(    blknorm < apxnorm * eps &
         .and. rownorm < apxnorm * eps &
         .and. colnorm < apxnorm * eps &
         .and. k>=param(64)) lstop_aca = 1;  
      endif
    endif
  if(.false.)then
!$omp critical
    print*,'pcol'; print*,pcol
    print*,'prow'; print*,prow
!$omp end critical
  endif
    if(lstop_aca==1 .and. k>=param(64)) exit
    k=k+1
  enddo
!  if(k==kmax .or. ntries_row==0 .or. ntries_col==0 .or. ntries==0)then
!    k=k-1
!  endif
  
  if(k<param(64))then
!$omp critical
    print*, 'colnorm=',colnorm,' rownorm=',rownorm,'ACA_EPS=',ACA_EPS
    print*, 'col_maxval=',col_maxval,' row_maxval=',row_maxval
    print*, 'ntries_row=',ntries_row,' ntries_col=',ntries_col,' ntries=',ntries
    print*, 'k=',k
!    k=k-1; if(k<1) stop
!    stop
!$omp end critical
  endif
  deallocate(lrow_msk,lcol_msk,pa_ref,pb_ref)
  HACApK_acaplus=k
!!!  print*,'HACApK_acaplus=',HACApK_acaplus
!!!  write(6,2000) 'blknorm=',blknorm/apxnorm,' colnorm=',colnorm/apxnorm,' rownorm=',rownorm/apxnorm
!!!  if(nstrtt==         113) stop
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
 mpinr=lpmd(3); mpilog=lpmd(4); nrank=lpmd(2); icomm=lpmd(1)
 eps=param(71); ACA_EPS=param(72)*eps; kparam=param(63)
!$OMP parallel default(none) &
!$OMP          shared(st_lf,lodl,st_bemv,znrmmat,lodt,lthr,param,mpilog,mpinr) &
!$OMP          private(zab,zaa,kt,ith,ith1,nths,nthr,nthe,ltmtx,ierr,ndl,ndt,ns,nstrtl,nstrtt,ip,il,it,ill,itt,kndl,kndt) &
!$OMP          firstprivate(eps, ACA_EPS, kparam)
 ith = omp_get_thread_num()
 nthr = omp_get_num_threads()
! if(nthr == 0) write(* ,*) nthr ,ith
 !$OMP barrier
 ith1 = ith+1
 nths=lthr(ith); nthe=lthr(ith1)-1
 if(param(1)>1) then
!$omp critical
   if(mpilog>0) write(mpilog,1000) 'sub HACApK_fill_leafmtx_hyp; nths=',nths,'; nthe=',nthe
!$omp end critical
 endif
 ierr=0
 do ip=nths,nthe
   ndl   =st_lf(ip)%ndl   ; ndt   =st_lf(ip)%ndt   ; ns=ndl*ndt
   nstrtl=st_lf(ip)%nstrtl; nstrtt=st_lf(ip)%nstrtt; ltmtx=st_lf(ip)%ltmtx
! write(mpilog,1000) 'sub HACApK_fill_leafmtx_hyp; ip=',ip,'; ndl=',ndl,'; ndt=',ndt,'; mpinr',mpinr

   if(ltmtx==1)then
!     kndl=ndl/2+1; kndt=ndt/2+1
!     kparam=min(kparam,kndl,kndt)
     allocate(zab(ndt,kparam),zaa(ndl,kparam),stat=ierr)
     if(ierr.ne.0) then
!$omp critical
        write(*,*) 'sub HACApK_fill_leafmtx_hyp; zab,zaa Memory allocation failed !'
        write(*,1000) 'ip=ip',ip,' ierr=',ierr,' ndt=',ndt,' ndl=',ndl,' kparam=',kparam
!$omp end critical
        stop 10
     endif
     if(param(60)==1)then
       kt=HACApK_aca(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
     elseif(param(60)==2)then
       kt=HACApK_acaplus(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
     else
       print*,'Only ACA and ACA+ is avairable! Set param(60)=1 or 2.'
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
        write(*,*) 'sub HACApK_fill_leafmtx_hyp; a1,a2 Memory allocation failed !'
        write(*,*) 'ip=ip',ip,' ierr=',ierr,' ndt=',ndt,' ndl=',ndl,' kt=',kt
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
        write(*,*) 'sub HACApK_fill_leafmtx_hyp; a2 Memory allocation failed !'
        write(*,*) 'ip=ip',ip,' ierr=',ierr,' ndt=',ndt,' ndl=',ndl
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
     elseif(param(60)==2)then
       kt=HACApK_acaplus(zaa,zab,param,ndl,ndt,nstrtl,nstrtt,lodl,st_bemv,kparam,eps,znrmmat,ACA_EPS)
     else
       print*,'Only ACA and ACA+ is avairable! Set param(60)=1 or 2.'
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
        print*, 'sub HACApK_fill_leafmtx; Memory allocation failed !'
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
 integer :: lnmtx(:)
 integer*8 :: mem8,memh
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10))
 2000 format(5(a,f12.2))
 
 param => st_ctl%param(:)
 lpmd => st_ctl%lpmd(:)
 mpinr=lpmd(3); mpilog=lpmd(4); icomm=lpmd(1); nrank=lpmd(2)
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
 if(st_ctl%param(1)>1)  write(mpilog,1000) ' ktmax=',ktmax
 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( ktmax, ktmaxmax, 1, MPI_INTEGER, MPI_MAX,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) ' Maximun of kt=',ktmaxmax
 call MPI_reduce( ktav, ktavav, 1, MPI_INTEGER, MPI_SUM,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0 .and. st_leafmtxp%nlfkt/=0) write(*,1000) ' Average of kt=',ktavav/st_leafmtxp%nlfkt

 znndense=real(nd)*real(nd)*8.0/1024.0/1024.0
! if(st_ctl%param(1)>0)  write(mpilog,2000) 'dense matrix memory=',znndense,'(Mbyte)'
 memh=(mem8+sum(lnmtx(:)))
 znn=real(memh)*8.0/1024.0/1024.0
! write(*,*) 'mpinr=',mpinr,'H matrix memory=',znn,'(Mbyte)'
 if(st_ctl%param(1)>1)  write(mpilog,2000) '    H matrix memory=',znn,'(Mbyte)'
 
 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( znn, znnmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX,0, icomm, ierr );
 call MPI_reduce( znn, znnmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN,0, icomm, ierr );
 call MPI_reduce( znn, znnall, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Memory of the H-matrix=',znnall,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Memory compression v.s. dense matrix=',znnall/znndense*100,'(%)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Maximun memory of sub-matrices for a MPI=',znnmax,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Minimun memory of sub-matrices for a MPI=',znnmin,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Minimun memory/Maximun memory=',znnmin/znnmax
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Whole/(Maximun*Nrank)=',znnall/(znnmax*nrank)
 
 end subroutine HACApK_chk_leafmtx
 
!***HACApK_chk_blrmtx
 subroutine HACApK_chk_blrmtx(st_leafmtxp,st_ctl,lnmtx,nd,mem8)
 include 'mpif.h'
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 integer :: lnmtx(:)
 integer*8 :: mem8,memh
 real*8,pointer :: param(:)
 integer*4,pointer :: lpmd(:)
 1000 format(5(a,i10))
 2000 format(5(a,f12.2))
 
 param => st_ctl%param(:)
 lpmd => st_ctl%lpmd(:)
 mpinr=lpmd(3); mpilog=lpmd(4); icomm=lpmd(1); nrank=lpmd(2)
 kparam=param(63)
 ktofl=0 ;ktoft=0 ;ktmax=0; mem8=0; ktav=0

 nlfalt=st_leafmtxp%nlfalt; nlfalt2=nlfalt/2
 il=st_leafmtxp%lbstrtl(nlfalt2)
 call MPI_Barrier( icomm, ierr )
 do ip=1,st_leafmtxp%nlf
     ndl=st_leafmtxp%st_lf(ip)%ndl; ndt=st_leafmtxp%st_lf(ip)%ndt; ns=ndl*ndt
     nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
     if(st_leafmtxp%st_lf(ip)%ltmtx==1)then
       kt=st_leafmtxp%st_lf(ip)%kt
       if(kt>ktmax) ktmax=kt
       if((nstrtl+ndl)>nd .and. nstrtt==1) then
          ktofl=kt
          if(st_ctl%param(1)>1 ) write(*,1000) 'leftlower  value of kt=',ktofl,' ndl=',ndl,' ndt=',ndt
       endif
       if((nstrtt+ndt)>nd .and. nstrtl==1) then
          ktoft=kt
          if(st_ctl%param(1)>1 ) write(*,1000) 'rightupper value of kt=',ktoft,' ndl=',ndl,' ndt=',ndt
       endif
       mem8=mem8+(ndt+ndl)*kt
       ktav=ktav+kt
       if(il==nstrtl) then
!       write(*,*) nstrtl,nstrtt,kt,mpinr
       endif
     elseif(st_leafmtxp%st_lf(ip)%ltmtx==2)then
       if(il==nstrtl) then
!       write(*,*) nstrtl,nstrtt,0,mpinr
       endif
       mem8=mem8+ns
     endif
 enddo
 st_leafmtxp%ktmax=ktmax
 if(st_ctl%param(1)>1)  write(mpilog,1000) ' ktmax=',ktmax
 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( ktmax, ktmaxmax, 1, MPI_INTEGER, MPI_MAX,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,1000) 'Maximun of kt=',ktmaxmax
 call MPI_reduce( ktav, ktavav, 1, MPI_INTEGER, MPI_SUM,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0 .and. st_leafmtxp%nlfkt/=0) write(*,1000) 'Average of kt=',ktavav/st_leafmtxp%nlfkt

 znndense=real(nd)*real(nd)*8.0/1024.0/1024.0
! if(st_ctl%param(1)>0)  write(mpilog,2000) 'dense matrix memory=',znndense,'(Mbyte)'
 memh=(mem8+sum(lnmtx(:)))
 znn=real(memh)*8.0/1024.0/1024.0
! write(*,*) 'mpinr=',mpinr,'H matrix memory=',znn,'(Mbyte)'
 if(st_ctl%param(1)>1)  write(mpilog,2000) '    H matrix memory=',znn,'(Mbyte)'
 
 call MPI_Barrier( icomm, ierr )
 call MPI_reduce( znn, znnmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX,0, icomm, ierr );
 call MPI_reduce( znn, znnmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN,0, icomm, ierr );
 call MPI_reduce( znn, znnall, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, icomm, ierr );
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Memory of the H-matrix=',znnall,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Memory compression v.s. dense matrix=',znnall/znndense*100,'(%)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Maximun memory of sub-matrices for a MPI=',znnmax,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Minimun memory of sub-matrices for a MPI=',znnmin,'(Mbyte)'
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Minimun memory/Maximun memory=',znnmin/znnmax
 if(st_ctl%param(1)>0 .and. mpinr==0) write(*,*) 'Whole/(Maximun*Nrank)=',znnall/(znnmax*nrank)
 
 end subroutine HACApK_chk_blrmtx
 
!***HACApK_count_blrnmb
 RECURSIVE subroutine HACApK_count_blrnmb(st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
 nleaf=param(42)+1; nlmax=param(22)*nofc
 
 ndpth=ndpth+1
 mdpth=param(53)
 
 if(ndpth==mdpth .or. (ndl<nleaf .and. ndt<nleaf))then
   lnmtx(4)=lnmtx(4)+1
    return
 endif

   lnmtx(3)=lnmtx(3)+1

 if(ndl<nleaf)then
   do it=1,nnsont
     call HACApK_count_blrnmb(st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nleaf)then
   do il=1,nnsonl
     call HACApK_count_blrnmb(st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_count_blrnmb(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine HACApK_count_blrnmb

!***HACApK_count_blrleaf
 RECURSIVE subroutine HACApK_count_blrleaf(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(:)
 integer*4 :: lnmtx(:),lpmd(*),lnmtx2(3)
 real*8 :: param(*)
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
 nlttcblk=param(42)+1; nlmax=param(22)*nofc
 nleaf=param(21)+1
 
 ndpth=ndpth+1
 mdpth=param(53)
 
! print*,'ndl=',ndl,'; ndt=',ndt,'; nlttcblk=',nlttcblk,'; nleaf=',nleaf
 
 if(ndpth==mdpth .or. (ndl<nlttcblk .and. ndt<nlttcblk))then
   lnmtx(4)=lnmtx(4)+1
   ibl=lnmtx(4)

 zs=0.0d0
 do id=1,st_cltl%ndim
  if(st_cltl%bmax(id)/=st_cltt%bmax(id)) zs=1.0d-16
  if(st_cltl%bmin(id)/=st_cltt%bmin(id)) zs=1.0d-16
 enddo

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
 
 if(st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt)then
    st_leafmtx(ibl)%nlf=1
    lnmtx(1)=lnmtx(1)+1
    return
 endif
!!! if(ndl<nleaf .and. ndt<nleaf)then
 if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. (ndl<=nleaf .and. ndt<=nleaf)))then
    st_leafmtx(ibl)%nlf=1
    lnmtx(2)=lnmtx(2)+1
    return
 endif
   lnmtx(3)=lnmtx(3)+1
   iblnlf=0
   do il=1,nnsonl
     do it=1,nnsont
       lnmtx2(:)=0
       call HACApK_count_lntmx(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx2,nofc,nffc,ndpth)
       lnmtx(1)=lnmtx(1)+lnmtx2(1); lnmtx(2)=lnmtx(2)+lnmtx2(2); lnmtx(3)=lnmtx(3)+lnmtx2(3)
       iblnlf=iblnlf+lnmtx2(1)+lnmtx2(2)
       ndpth=ndpth-1
     enddo
   enddo
   if(iblnlf==0)then
     print*,'Error; HACApK_count_blrleaf; iblnlf==0'
     print*,'ndl=',ndl,'; ndt=',ndt,'; nstrtl=',nstrtl,'; nstrtt=',nstrtt
     stop
   endif
   st_leafmtx(ibl)%nlf=iblnlf
   allocate(st_leafmtx(ibl)%st_lf(iblnlf))
   return
 endif

   lnmtx(3)=lnmtx(3)+1

 if(ndl<nlttcblk)then
   do it=1,nnsont
     call HACApK_count_blrleaf(st_leafmtx,st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nlttcblk)then
   do il=1,nnsonl
     call HACApK_count_blrleaf(st_leafmtx,st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_count_blrleaf(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine HACApK_count_blrleaf

!***HACApK_generate_blrleaf
 RECURSIVE subroutine HACApK_generate_blrleaf(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt
 nlttcblk=param(42)+1; nlmax=param(22)*nofc
! print*,'nlttcblk=',nlttcblk; stop
 nleaf=param(21)+1
 
 ndpth=ndpth+1
 mdpth=param(53)
 
 if(ndpth==mdpth .or. (ndl<nlttcblk .and. ndt<nlttcblk))then
   lnmtx(4)=lnmtx(4)+1
   ibl=lnmtx(4)

   zs=0.0d0
   do id=1,st_cltl%ndim
    if(st_cltl%bmax(id)/=st_cltt%bmax(id)) zs=1.0d-16
    if(st_cltl%bmin(id)/=st_cltt%bmin(id)) zs=1.0d-16
   enddo

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
 
   nlf=nlf+1
   st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
   st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
   if(st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt)then
      st_leafmtx(nlf)%kt=0
      st_leafmtx(nlf)%ltmtx=1
!      print*,'ibl=',ibl,'; iblnlf=',1,'; ltmtx=',1
      return
!!!   elseif(ndl<nleaf .and. ndt<nleaf)then
   elseif(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. (ndl<=nleaf .and. ndt<=nleaf)))then
      st_leafmtx(nlf)%kt=0
      st_leafmtx(nlf)%ltmtx=2
!      print*,'ibl=',ibl,'; iblnlf=',1,'; ltmtx=',2
      return     
   else
!!!     allocate(st_leafmtx(ibl)%st_lf(nnsonl*nnsont))
     iblnlf=0
     do il=1,nnsonl
       do it=1,nnsont
         call HACApK_generate_leafmtx(st_leafmtx(ibl)%st_lf,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,iblnlf,ndpth)
         ndpth=ndpth-1
       enddo
     enddo
     st_leafmtx(nlf)%kt=0
     st_leafmtx(nlf)%ltmtx=4
!     print*,'ibl=',ibl,'; iblnlf=',st_leafmtx(ibl)%nlf,'; ltmtx=',4
   return
 
 endif
 endif

 if(ndl<nlttcblk)then
   do it=1,nnsont
     call HACApK_generate_blrleaf(st_leafmtx,st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nlttcblk)then
   do il=1,nnsonl
     call HACApK_generate_blrleaf(st_leafmtx,st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_generate_blrleaf(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine
 
!***HACApK_count_LH
 RECURSIVE subroutine HACApK_count_LH(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_cluster),pointer :: st_ls,st_ts
 type(st_HACApK_leafmtx) :: st_leafmtx(:)
 integer*4 :: lnmtx(:),lpmd(*),lnmtx2(3)
 real*8 :: param(*)
 
 ndpth=ndpth+1
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt
 nlttcblk=st_cltl%pc_sons(1)%nsize
 nleaf=param(21)
 mdpth=param(53)
 
 if(mdpth==1 .or. nleaf>=nofc)then
   st_leafmtx(1)%nlf=1; lnmtx(2)=1; return
 else
   lnmtx(3)=lnmtx(3)+1
   iblnlf=0
 endif
 
 do il=1,nnsonl
   st_ls=>st_cltl%pc_sons(il)
   do it=1,nnsont
     st_ts=>st_cltt%pc_sons(it)
     ibl=it+(il-1)*nnsont
     zs=0.0d0
     do id=1,st_cltl%ndim
       if(st_ls%bmax(id)/=st_ts%bmax(id)) zs=1.0d-16
       if(st_ls%bmin(id)/=st_ts%bmin(id)) zs=1.0d-16
     enddo

     do id=1,st_cltl%ndim
       if(st_ls%bmax(id)<st_ts%bmin(id))then
         zs=zs+(st_ts%bmin(id)-st_ls%bmax(id))*(st_ts%bmin(id)-st_ls%bmax(id))
       elseif(st_ts%bmax(id)<st_ls%bmin(id))then
         zs=zs+(st_ls%bmin(id)-st_ts%bmax(id))*(st_ls%bmin(id)-st_ts%bmax(id))
       else
       endif
     enddo
     zdistlt=dsqrt(zs)
     zeta=param(51)
     
     if(st_ls%zwdth<=zeta*zdistlt)then
       st_leafmtx(ibl)%nlf=1
       lnmtx(1)=lnmtx(1)+1
     elseif(mdpth==2 .or. nleaf>=nlttcblk)then
       st_leafmtx(ibl)%nlf=1
       lnmtx(2)=lnmtx(2)+1
     else
       lnmtx(3)=lnmtx(3)+1
       lnmtx2(:)=0
       call HACApK_count_lntmx(st_ls,st_ts,param,lpmd,lnmtx2,nofc,nffc,ndpth)
       lnmtx(1)=lnmtx(1)+lnmtx2(1); lnmtx(2)=lnmtx(2)+lnmtx2(2); lnmtx(3)=lnmtx(3)+lnmtx2(3)
       iblnlf=lnmtx2(1)+lnmtx2(2)
       st_leafmtx(ibl)%nlf=iblnlf
       allocate(st_leafmtx(ibl)%st_lf(iblnlf))
       ndpth=ndpth-1
     endif
   enddo
 enddo

 end subroutine

!***HACApK_generate_LH
 RECURSIVE subroutine HACApK_generate_LH(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_cluster),pointer :: st_ls,st_ts
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 
 ndpth=ndpth+1
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt
 nlttcblk=st_cltl%pc_sons(1)%nsize
 nleaf=param(21)
 mdpth=param(53)
 
 if(mdpth==1 .or. nleaf>=nofc)then
   ibl=1
   st_leafmtx(ibl)%nstrtl=nstrtl; st_leafmtx(ibl)%ndl=ndl;
   st_leafmtx(ibl)%nstrtt=nstrtt; st_leafmtx(ibl)%ndt=ndt;
   st_leafmtx(ibl)%kt=0
   st_leafmtx(ibl)%ltmtx=2
   return
 endif
 
 do il=1,nnsonl
   st_ls=>st_cltl%pc_sons(il)
   do it=1,nnsont
     st_ts=>st_cltt%pc_sons(it)
     ibl=it+(il-1)*nnsont
     st_leafmtx(ibl)%nstrtl=st_ls%nstrt; st_leafmtx(ibl)%ndl=st_ls%nsize*nffc;
     st_leafmtx(ibl)%nstrtt=st_ts%nstrt; st_leafmtx(ibl)%ndt=st_ts%nsize*nffc;
     st_leafmtx(ibl)%kt=0

     zs=0.0d0
     do id=1,st_cltl%ndim
       if(st_ls%bmax(id)/=st_ts%bmax(id)) zs=1.0d-16
       if(st_ls%bmin(id)/=st_ts%bmin(id)) zs=1.0d-16
     enddo

     do id=1,st_cltl%ndim
       if(st_ls%bmax(id)<st_ts%bmin(id))then
         zs=zs+(st_ts%bmin(id)-st_ls%bmax(id))*(st_ts%bmin(id)-st_ls%bmax(id))
       elseif(st_ts%bmax(id)<st_ls%bmin(id))then
         zs=zs+(st_ls%bmin(id)-st_ts%bmax(id))*(st_ls%bmin(id)-st_ts%bmax(id))
       else
       endif
     enddo
     zdistlt=dsqrt(zs)
     zeta=param(51)
     
     if(st_ls%zwdth<=zeta*zdistlt)then
       st_leafmtx(ibl)%ltmtx=1
!       print*,'ibl=',ibl,'; ltmtx=',1
     elseif(mdpth==2 .or. nleaf>=nlttcblk)then
       st_leafmtx(ibl)%ltmtx=2
!       print*,'ibl=',ibl,'; ltmtx=',2
     else
       iblnlf=0
       st_leafmtx(ibl)%ltmtx=3
!       print*,'ibl=',ibl,'; ltmtx=',3
       call HACApK_generate_leafmtx(st_leafmtx(ibl)%st_lf,st_ls,st_ts,param,lpmd,lnmtx,nofc,nffc,iblnlf,ndpth)
       ndpth=ndpth-1
     endif
   enddo
 enddo
 end subroutine 

!***HACApK_count_blr
 RECURSIVE subroutine HACApK_count_blr(st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
 nleaf=param(21)+1; nlmax=param(22)*nofc
 
 ndpth=ndpth+1
 mdpth=param(53)
 
 if(ndpth==mdpth .or. (ndl<nleaf .and. ndt<nleaf))then
 
 zs=0.0d0
 do id=1,st_cltl%ndim
  if(st_cltl%bmax(id)/=st_cltt%bmax(id)) zs=1.0d-16
  if(st_cltl%bmin(id)/=st_cltt%bmin(id)) zs=1.0d-16
 enddo
 
 do id=1,st_cltl%ndim
   if(st_cltl%bmax(id)<st_cltt%bmin(id))then
     zs=zs+(st_cltt%bmin(id)-st_cltl%bmax(id))*(st_cltt%bmin(id)-st_cltl%bmax(id))
   elseif(st_cltt%bmax(id)<st_cltl%bmin(id))then
     zs=zs+(st_cltl%bmin(id)-st_cltt%bmax(id))*(st_cltl%bmin(id)-st_cltt%bmax(id))
   else
   endif
 enddo
 zdistlt=dsqrt(zs)
! print*,'zdistlt=',zdistlt
 zeta=param(51)
 
 if(st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt)then
!    if((nstrtl+ndl)/=nstrtt .and. (nstrtt+ndt)/=nstrtl)then
!    if((nstrtt+ndt)/=nstrtl)then
!    if((nstrtl+ndl)/=nstrtt)then
       lnmtx(1)=lnmtx(1)+1
       return
!    endif
 endif
 lnmtx(2)=lnmtx(2)+1
 return

 endif

   lnmtx(3)=lnmtx(3)+1

 if(ndl<nleaf)then
   do it=1,nnsont
     call HACApK_count_blr(st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nleaf)then
   do il=1,nnsonl
     call HACApK_count_blr(st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_count_blr(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine HACApK_count_blr

!***HACApK_generate_blr
 RECURSIVE subroutine HACApK_generate_blr(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt
 nleaf=param(21)+1; nlmax=param(22)*nofc
! print*,'nlmax=',nlmax; stop
 
 ndpth=ndpth+1
 mdpth=param(53)
 
 if(ndpth==mdpth .or. (ndl<nleaf .and. ndt<nleaf))then

 zs=0.0d0
 do id=1,st_cltl%ndim
  if(st_cltl%bmax(id)/=st_cltt%bmax(id)) zs=1.0d-16
  if(st_cltl%bmin(id)/=st_cltt%bmin(id)) zs=1.0d-16
 enddo
 
 do id=1,st_cltl%ndim
   if(st_cltl%bmax(id)<st_cltt%bmin(id))then
     zs=zs+(st_cltt%bmin(id)-st_cltl%bmax(id))*(st_cltt%bmin(id)-st_cltl%bmax(id))
   elseif(st_cltt%bmax(id)<st_cltl%bmin(id))then
     zs=zs+(st_cltl%bmin(id)-st_cltt%bmax(id))*(st_cltl%bmin(id)-st_cltt%bmax(id))
   else
   endif
 enddo
 zdistlt=dsqrt(zs)
! print*,'zdistlt=',zdistlt
 zeta=param(51)
 
 if(st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt)then
!    if((nstrtl+ndl)/=nstrtt .and. (nstrtt+ndt)/=nstrtl)then
!    if((nstrtt+ndt)/=nstrtl)then
!    if((nstrtl+ndl)/=nstrtt)then
      nlf=nlf+1
      st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
      st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
      st_leafmtx(nlf)%kt=0
      st_leafmtx(nlf)%ltmtx=1
      return
!    endif
 endif
   nlf=nlf+1
   st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
   st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
   st_leafmtx(nlf)%ltmtx=2
!     allocate(st_leafmtx(nlf)%a1(ndt,ndl))
   return
 endif

 if(ndl<nleaf)then
   do it=1,nnsont
     call HACApK_generate_blr(st_leafmtx,st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nleaf)then
   do il=1,nnsonl
     call HACApK_generate_blr(st_leafmtx,st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_generate_blr(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine HACApK_generate_blr
 
!***HACApK_count_zlntmx
 RECURSIVE subroutine HACApK_count_zlntmx(st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
 nleaf=param(21)+1; nlmax=param(22)*nofc
 
 ndpth=ndpth+1
 mdpth=param(53)
 
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
 
 if((st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt) &
    .and.(ndl<=nlmax .and. ndt<=nlmax) &
!    .and.(ndl>=nleaf .or. ndt>=nleaf) &
  )then
    if(param(52)==0 &
      .or. param(52)==1 .and.((nstrtl+ndl)/=nstrtt .and. (nstrtt+ndt)/=nstrtl))then
       lnmtx(1)=lnmtx(1)+1
       return
    elseif(param(52)/=1)then
       print*,'Invalid admissiblity!; Set param(52)=0 or 1.'
       stop
    endif
 endif
 if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. (ndl<=nleaf .and. ndt<=nleaf)))then
   lnmtx(2)=lnmtx(2)+1
   return
 endif
 lnmtx(3)=lnmtx(3)+1
 if(ndl<nleaf)then
   do it=1,nnsont
     call HACApK_count_zlntmx(st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nleaf)then
   do il=1,nnsonl
     call HACApK_count_zlntmx(st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_count_zlntmx(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine HACApK_count_zlntmx
 
!***HACApK_generate_zleafmtx
 RECURSIVE subroutine HACApK_generate_zleafmtx(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt
 nleaf=param(21)+1; nlmax=param(22)*nofc
! print*,'nlmax=',nlmax; stop

 ndpth=ndpth+1
 mdpth=param(53)
 
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
 
 if((st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt) &
    .and.(ndl<=nlmax .and. ndt<=nlmax) &
!    .and.(ndl>=nleaf .or. ndt>=nleaf) &
  )then
    if(param(52)==0 &
      .or. param(52)==1 .and.((nstrtl+ndl)/=nstrtt .and. (nstrtt+ndt)/=nstrtl))then
      nlf=nlf+1
      st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
      st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
      st_leafmtx(nlf)%kt=0
      st_leafmtx(nlf)%ltmtx=1
      return
    endif
 endif
 if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. (ndl<=nleaf .and. ndt<=nleaf)))then
   nlf=nlf+1
   st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
   st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
   st_leafmtx(nlf)%ltmtx=2
   return
 endif
 if(ndl<nleaf)then
   do it=1,nnsont
     call HACApK_generate_zleafmtx(st_leafmtx,st_cltl,st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 elseif(ndt<nleaf)then
   do il=1,nnsonl
     call HACApK_generate_zleafmtx(st_leafmtx,st_cltl%pc_sons(il),st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 else
   do il=1,nnsonl
     do it=1,nnsont
       call HACApK_generate_zleafmtx(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
       ndpth=ndpth-1
     enddo
   enddo
 endif
 end subroutine HACApK_generate_zleafmtx

!***HACApK_count_lntmx
 RECURSIVE subroutine HACApK_count_lntmx(st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
 nleaf=param(21)+1; nlmax=param(22)*nofc
 
! print*,'ndl=',ndl,'; nleaf=',nleaf; stop
 
 ndpth=ndpth+1
 mdpth=param(53)
 
 zs=0.0d0
 do id=1,st_cltl%ndim
  if(st_cltl%bmax(id)/=st_cltt%bmax(id)) zs=1.0d-16
  if(st_cltl%bmin(id)/=st_cltt%bmin(id)) zs=1.0d-16
 enddo

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
 
 if((st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt) &
    .and.(ndl>=nleaf .and. ndt>=nleaf .and. ndl<=nlmax .and. ndt<=nlmax) &
  )then
    if(param(52)==0 &
      .or. param(52)==1 .and.((nstrtl+ndl)/=nstrtt .and. (nstrtt+ndt)/=nstrtl))then
       lnmtx(1)=lnmtx(1)+1
       return
    elseif(param(52)/=1)then
       print*,'Invalid admissiblity!; Set param(52)=0 or 1.'
       stop
    endif
 endif
! if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. (ndl<=nleaf .and. ndt<=nleaf)))then
 if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. ndl<=nleaf .or. ndt<=nleaf))then
! if((nnsonl==0 .or. nnsont==0 .or. ndl<=nleaf .or. ndt<=nleaf))then
   lnmtx(2)=lnmtx(2)+1
   return
 endif
 lnmtx(3)=lnmtx(3)+1
 do il=1,nnsonl
   do it=1,nnsont
     call HACApK_count_lntmx(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,ndpth)
     ndpth=ndpth-1
   enddo
 enddo
 end subroutine HACApK_count_lntmx
 
!***HACApK_generate_leafmtx
 RECURSIVE subroutine HACApK_generate_leafmtx(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(:),lpmd(*)
 real*8 :: param(*)
 
 ndl=st_cltl%nsize*nffc; ndt=st_cltt%nsize*nffc
 nstrtl=st_cltl%nstrt; nstrtt=st_cltt%nstrt
 nnsonl=st_cltl%nnson; nnsont=st_cltt%nnson
! print*, nnsonl,ndl,nnsont,ndt
 nleaf=param(21)+1; nlmax=param(22)*nofc
! print*,'nlmax=',nlmax; stop

 ndpth=ndpth+1
 mdpth=param(53)
 
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
 
 if((st_cltl%zwdth<=zeta*zdistlt .or. st_cltt%zwdth<=zeta*zdistlt) &
    .and.(ndl>=nleaf .and. ndt>=nleaf .and. ndl<=nlmax .and. ndt<=nlmax) &
  )then
    if(param(52)==0 &
      .or. param(52)==1 .and.((nstrtl+ndl)/=nstrtt .and. (nstrtt+ndt)/=nstrtl))then
      nlf=nlf+1
      st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
      st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
      st_leafmtx(nlf)%kt=0
      st_leafmtx(nlf)%ltmtx=1
      return
    endif
 endif
! if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. (ndl<=nleaf .and. ndt<=nleaf)))then
 if(ndpth==mdpth .or. (nnsonl==0 .or. nnsont==0 .or. ndl<=nleaf .or. ndt<=nleaf))then
! if((nnsonl==0 .or. nnsont==0 .or. ndl<=nleaf .or. ndt<=nleaf))then
   nlf=nlf+1
   st_leafmtx(nlf)%nstrtl=nstrtl; st_leafmtx(nlf)%ndl=ndl;
   st_leafmtx(nlf)%nstrtt=nstrtt; st_leafmtx(nlf)%ndt=ndt;
   st_leafmtx(nlf)%ltmtx=2
!     allocate(st_leafmtx(nlf)%a1(ndt,ndl))
   return
 endif
 do il=1,nnsonl
   do it=1,nnsont
     call HACApK_generate_leafmtx(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf,ndpth)
     ndpth=ndpth-1
   enddo
 enddo
 end subroutine HACApK_generate_leafmtx
 
!***HACApK_count_olntmx
 RECURSIVE subroutine HACApK_count_olntmx(st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 integer*4 :: lnmtx(:),lpmd(*)
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
         call HACApK_count_olntmx(st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc)
       enddo
     enddo
   endif
 endif
 end subroutine HACApK_count_olntmx
 
!***HACApK_generate_oleafmtx
 RECURSIVE subroutine HACApK_generate_oleafmtx(st_leafmtx,st_cltl,st_cltt,param,lpmd,lnmtx,nofc,nffc,nlf)
 type(st_HACApK_cluster) :: st_cltl,st_cltt
 type(st_HACApK_leafmtx) :: st_leafmtx(*)
 integer*4 :: lnmtx(:),lpmd(*)
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
         call HACApK_generate_oleafmtx(st_leafmtx,st_cltl%pc_sons(il),st_cltt%pc_sons(it),param,lpmd,lnmtx,nofc,nffc,nlf)
       enddo
     enddo
   endif
 endif
 end subroutine HACApK_generate_oleafmtx
 
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
 
!***HACApK_sort_leafmtx
 subroutine HACApK_sort_leafmtx(st_leafmtx,nlf)
 type(st_HACApK_leafmtx) :: st_leafmtx(:)
 call HACApK_qsort_row_leafmtx(st_leafmtx,1,nlf)
 ilp=1; ips=1
 do ip=1,nlf
   il=st_leafmtx(ip)%nstrtl
   if(il<ilp)then    ; print *,'Error!; HACApK_sort_leafmtx row_sort';
   elseif(il>ilp)then; 
     call HACApK_qsort_col_leafmtx(st_leafmtx,ips,ip-1)
     ilp=il; ips=ip
   endif
 enddo
 call HACApK_qsort_col_leafmtx(st_leafmtx,ips,nlf)
 
 end subroutine HACApK_sort_leafmtx
 
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

!***HACApK_generate_tctree
 RECURSIVE subroutine HACApK_generate_tctree(st_clt,zgmid,param,lpmd,lod,ndpth,ndscd,nsrt,nd,md,ndim,nclst,nblall)
   type(st_HACApK_cluster) :: st_clt,st_cltf
   real*8 :: zgmid(md,ndim)
   real*8,dimension(:),allocatable :: zlmin,zlmax
   integer*4 :: lod(md),lpmd(*)
   real*8 :: param(*),param1(100)
   
   param1(:)=param(:100); 
   param1(21)=param(42)/2+1
!   param1(21)=2
   nsrt1=1; nclst1=0; ndpth1=0; ndscd1=0
   call HACApK_generate_cbitree(st_cltf,zgmid,param1,lpmd,lod,ndpth1,ndscd1,nsrt1,nd,md,ndim,nclst1)
   nbsize=param(42)
   nson=md/param(42); if(nson==1 .or. param(53)<=1) nson=0
   
   ndpth=ndpth+1
!   print*,'nbsize=',nbsize, 'nson=',nson, ' md=',md
   st_clt=HACApK_generate_cluster(nclst,ndpth,nsrt,nd,ndim,nson)
   st_clt%ndscd=nd
   if(nson==0)then
     nblall=1; return
   else
     nblall=nson*nson
   endif
   
!   nson=2
!   print*,'nbsize=',nbsize,'; nblcoks=',nson
   ngson=0; ndpthgson=2
   ndgson=0; 
   do ic=1,nson
     ndgson=(nd+ic-1)/nson
!     print*,'nsrt=',nsrt,'; ndgson=',ndgson
!     st_clt%pc_sons(ic)=HACApK_generate_cluster(nclst,ndpthgson,nsrt,ndgson,ndim,ngson)
     call HACApK_generate_cbitree(st_clt%pc_sons(ic),zgmid,param,lpmd,lod,ndpth,ndscd,nsrt,ndgson,md,ndim,nclst)
     ndpth=ndpth-1
     nsrt=nsrt+ndgson
   enddo
 end subroutine

!***HACApK_generate_cbitree
 RECURSIVE subroutine HACApK_generate_cbitree(st_clt,zgmid,param,lpmd,lod,ndpth,ndscd,nsrt,nd,md,ndim,nclst)
   type(st_HACApK_cluster) :: st_clt
   real*8 :: zgmid(md,ndim)
   real*8,dimension(:),allocatable :: zlmin,zlmax
   integer*4 :: lod(md),lpmd(*)
   real*8 :: param(*)
   minsz=param(21)
!   minsz=param(21)/4+1
   ndpth=ndpth+1
!   ndscd=ndscd+1
!  if(i>26) stop
!   print*,''
!   print*,'    nsrt=',nsrt,' nd=',nd
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


!***HACApK_gen_mat_plot
    subroutine HACApK_gen_mat_plot(st_leafmtxp,lpmd,lthr)
    type(st_HACApK_leafmtxp) :: st_leafmtxp
    integer*4 :: lpmd(1:*),lthr(0:*),kt
    character*32 :: plot
    mpinr=lpmd(3); mpilog=lpmd(4); nthr=lpmd(20)
    nlf=st_leafmtxp%nlf

    write(plot,'(a,i4.4,a)')'plot',mpinr
    open(114,file=plot)
    
    !write .txt file for h-mat assign
    do ith=0,nthr-1
       is=lthr(ith); ie=lthr(ith+1)-1
       do ip=is,ie
          nstrtl=st_leafmtxp%st_lf(ip)%nstrtl;  nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
          ndl=st_leafmtxp%st_lf(ip)%ndl;        ndt=st_leafmtxp%st_lf(ip)%ndt
          nsub=st_leafmtxp%st_lf(ip)%ltmtx
          write(114,'(i4,a,i10,a,i10,a,i10,a,i10,a,i2)')&
               mpinr*100+ith,',',nstrtt,',',nstrtl,',',nstrtt+ndt-1,',',nstrtl+ndl-1,',',nsub
       enddo
    enddo
    close(114)
  end subroutine HACApK_gen_mat_plot

!***HACApK_write_leafmtx
 subroutine HACApK_write_leafmtx(st_leafmtxp,st_ctl)
 type(st_HACApK_leafmtxp) :: st_leafmtxp
 type(st_HACApK_lcontrol) :: st_ctl
 
 nd=st_leafmtxp%nd; nlf=st_leafmtxp%nlf
 open( 11, file="hmatrix.txt", action='write', iostat=ierr )
 write(11,*) nd

 do il=1,nd
   write(11,*) st_ctl%lod(il)
 enddo

 write(11,*) nlf 

 do ip=1,nlf
   ltmtx=st_leafmtxp%st_lf(ip)%ltmtx  ; kt    =st_leafmtxp%st_lf(ip)%kt
   ndl   =st_leafmtxp%st_lf(ip)%ndl   ; ndt   =st_leafmtxp%st_lf(ip)%ndt
   nstrtl=st_leafmtxp%st_lf(ip)%nstrtl; nstrtt=st_leafmtxp%st_lf(ip)%nstrtt
   write(11,*) ltmtx,ndl,ndt,nstrtl,nstrtt,kt

   if(ltmtx==1)then ! Low-rank matrix
     do il=1,kt
       do it=1,ndt
!         write(11,'(es,$)') st_leafmtxp%st_lf(ip)%a1(it,il)
       enddo
       write(11,*)
       do it=1,ndl
!         write(11,'(es,$)') st_leafmtxp%st_lf(ip)%a2(it,il)
       enddo
       write(11,*)
     enddo
   elseif(ltmtx==2)then ! Dense matrix
     do il=1,ndl
       do it=1,ndt
!         write(11,'(es,$)') st_leafmtxp%st_lf(ip)%a1(it,il)
       enddo
       write(11,*)
     enddo
   endif

 enddo

 close( 11 )

 end subroutine
 
!***HACApK_gen_lattice_vector
 integer function HACApK_gen_lattice_vector(st_vec,st_leafmtxp,st_ctl)
   type(st_HACApK_leafmtxp) :: st_leafmtxp
   type(st_HACApK_lcontrol) :: st_ctl
   type(st_HACApK_latticevec) :: st_vec
 1000 format(5(a,i10)/)
 2000 format(5(a,f10.4)/)

 lrtrn=0
 mpinr=st_ctl%lpmd(3); mpilog=st_ctl%lpmd(4); nrank=st_ctl%lpmd(2); icomm=st_ctl%lpmd(1); nthr=st_ctl%lpmd(20)
 if(st_ctl%param(8)/=10.and.st_ctl%param(8)/=20)then
   if(mpinr==0) then
     print*,'ERROR!; sub HACApK_gen_lattice_vector; param(8)=',st_ctl%param(8)
     print*,'ERROR!; param(8) must be 10(BLR) or 20(Lattic H)'
   endif
   stop
 endif
 if(st_ctl%lpmd(32)/=st_ctl%lpmd(36))then
   if(mpinr==0) then
     print*,'ERROR!; sub HACApK_gen_lattice_vector; Process grid must have square shape!'
   endif
   stop
 endif
!!! ndc=st_leafmtxp%ndlfs; nlfc=st_leafmtxp%nlft
 ndc=st_leafmtxp%ndtfs; nlfc=st_leafmtxp%nlft
 st_vec%ndc=ndc
 st_vec%nlfc=nlfc
 allocate(st_vec%lbstrtc(nlfc),st_vec%lbndc(nlfc),st_vec%lodc(ndc),st_vec%vs(ndc))
 st_vec%vs(:)=0.0d0
 nlfalt=st_leafmtxp%nlfalt
 ip=1
 do il=1,nlfc
   ibl=st_leafmtxp%lnlfl2g(il,1); ibl=mod(ibl-1,nlfalt)+1
   ndilb=st_leafmtxp%lbndt(ibl); lstrtlibl=st_leafmtxp%lbstrtt(ibl)
   st_vec%lbndc(il)=ndilb
   st_vec%lbstrtc(il)=lstrtlibl
   st_vec%lodc(ip:ip+ndilb-1)=st_ctl%lod(lstrtlibl:lstrtlibl+ndilb-1)
   ip=ip+ndilb
 enddo
 
  
 if(st_ctl%param(1)>1) then
   write(mpilog,1000) 'local vec size; ndc=',ndc,'; #lattice blocks=',nlfc
   write(mpilog,*) 'st_vec%lbstrtc='
   write(mpilog,'(10(i9))') st_vec%lbstrtc(:)
   write(mpilog,*) 'st_vec%lbndc='
   write(mpilog,'(10(i9))') st_vec%lbndc(:)
!   write(mpilog,*) 'st_vec%lodc='
!   write(mpilog,'(10(i9))') st_vec%lodc(:)
 endif
! stop
 HACApK_gen_lattice_vector=lrtrn
endfunction
 
endmodule m_HACApK_base


