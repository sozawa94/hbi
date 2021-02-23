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
!***dotp_d
 real*8 function HACApK_dotp_d(nd,za,zb)
 real*8 :: za(nd),zb(nd)
 HACApK_dotp_d=sum(za(:nd)*zb(:nd))
 end function
!***unrm_d
 real*8 function HACApK_unrm_d(nd,za)
 implicit real*8(a-h,o-z)
 real*8 :: za(:)
 zscl=0.0d0; zzz=1.0d0
 do il=1,nd
   zza=abs(za(il))
   if(zza<1.0d-30) cycle
   if(zscl<zza)then
     zzz=1.0d0+zzz*(zscl/zza)*(zscl/zza)
     zscl=zza
   else
     zzz=zzz+(zza/zscl)*(zza/zscl)
   endif
 enddo
 HACApK_unrm_d=zscl*dsqrt(zzz)
 end function
!*** maxabsvalloc_d
 subroutine HACApK_maxabsvalloc_d(za,zz,il,nd)
 implicit real*8(a-h,o-z)
  real*8 :: za(:)
  il = 1; zz = dabs(za(1))
  do it=2,nd
    if(dabs(za(it)) > zz)then
      il = it; zz = dabs(za(it))
    endif
  enddo
 endsubroutine
!*** maxabsvallocm_d
 subroutine HACApK_maxabsvallocm_d(za,zz,il,nd,lmask)
 implicit real*8(a-h,o-z)
 real*8 :: za(:)
 integer lmask(:)
  il = 0; zz = 0.0d0
  do it=1,nd
    if(lmask(it)==1) cycle
    if(dabs(za(it)) > zz)then
      il = it; zz = dabs(za(it))
    endif
  enddo
 endsubroutine
!*** minabsvalloc_d
 subroutine HACApK_minabsvalloc_d(za,zz,il,nd)
 implicit real*8(a-h,o-z)
  real*8 za(:)
  il = 1; zz = dabs(za(1))
  do it=2,nd
    if(dabs(za(it)) < zz)then
      il = it; zz = dabs(za(it))
    endif
  enddo
 endsubroutine

!***adot_dsm
 subroutine HACApK_adot_dsm(zau,zaa,zu,ndl,ndt,mdl)
 implicit real*8(a-h,o-z)
 real*8 :: zaa(:,:)
 real*8 :: zu(:),zau(:)

 zau=0.0d0
 do it=1,ndt
   zau(1:ndl)=zau(1:ndl)+zaa(1:ndl,it)*zu(it)
 enddo

 end subroutine HACApK_adot_dsm

!***adotsub_dsm
 subroutine HACApK_adotsub_dsm(zr,zaa,zu,ndl,ndt,mdl)
 implicit real*8(a-h,o-z)
 real*8 :: zaa(:,:)
 real*8 :: zu(:),zr(:)
 real*8,dimension(:),allocatable :: zau

 interface
    subroutine HACApK_adot_dsm(zau,zaa,zu,ndl,ndt,mdl)
      implicit real*8(a-h,o-z)
      real*8 :: zaa(:,:)
      real*8 :: zu(:),zau(:)
    end subroutine HACApK_adot_dsm
 end interface

 allocate(zau(ndl))
 call HACApK_adot_dsm(zau,zaa,zu,ndl,ndt,mdl)
 zr(1:ndl)=zr(1:ndl)-zau(1:ndl)
 deallocate(zau)

 end subroutine HACApK_adotsub_dsm

!***adotprt_dsm
 subroutine HACApK_adotprt_dsm(zau,zaa,zu,ndl,ndt,mdl,ils,ile)
 implicit real*8(a-h,o-z)
 real*8 :: zaa(mdl,*)
 real*8 :: zu(ndt),zau(ndl)
 zau(:)=0.0d0
 do it=1,ndt
   zau(ils:ile)=zau(ils:ile)+zaa(ils:ile,it)*zu(it)
 enddo
 end subroutine HACApK_adotprt_dsm
!***adotsubprt_dsm
 subroutine HACApK_adotsubprt_dsm(zr,zaa,zu,ndl,ndt,mdl,ils,ile)
 implicit real*8(a-h,o-z)
 real*8 :: zaa(mdl,*)
 real*8 :: zu(ndt),zr(ndl)
 real*8,dimension(:),allocatable :: zau
 allocate(zau(ndl))
! call adot_dsm(zau,zaa,zu,ndl,ndt,mdl,ils,ile)
 zr(ils:ile)=zr(ils:ile)-zau(ils:ile)
 deallocate(zau)
 end subroutine HACApK_adotsubprt_dsm
!***strtend_omp
 subroutine HACApK_strtend_omp(j0,j1,ith,nth,ndim)
 implicit real*8(a-h,o-z)
  inth = ith + nth;     nth2 = nth + nth
  ntdim = ndim / nth;    mth = ndim - ntdim*nth
  ntdim1 = ntdim + 1;     j0 = 0
  if ( ith .lt. mth ) then
    do i=1, ith; j0 = j0 + ntdim1; enddo
    j1 = j0 + ntdim1
  else
    do i=1, mth; j0 = j0 + ntdim1; enddo
    do i=mth+1, ith; j0 = j0 + ntdim; enddo
    j1 = j0 + ntdim
  endif
  j0=j0+1
 end subroutine HACApK_strtend_omp
!***HACApK_med3
 integer function HACApK_med3(nl,nr,nlr2)
    if(nl < nr)then
      if (nr < nlr2)then;  HACApK_med3=nr; elseif (nlr2 < nl)then; HACApK_med3=nl; else;  HACApK_med3=nlr2; endif
    else
      if (nlr2 < nr)then;  HACApK_med3=nr; elseif (nl < nlr2)then; HACApK_med3=nl; else; HACApK_med3=nlr2; endif
    endif
 endfunction
