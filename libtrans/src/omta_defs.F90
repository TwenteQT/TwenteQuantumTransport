 module omta_defs
 use sparselib
 implicit none
!!$ Global parameters:
!!$ prec      -  kind value for double-precision variables
 integer, parameter :: prec=8,plen=4

!!$ File identifiers
 integer, parameter ::  ctrlf=21, ologf=22, ifi=30
      
 type tclass
!!$ Atomic class, defined by the size of accosiated sphere and potential within.
!!$ Uniquely identified by the "label" field
   character(len=20) :: label
   real(kind=prec)  :: s, ra, rb  ! the radii of the potential sphere, mesh parameters
   ! radial mesh, potential, screening radius, sigma, screening const., t(1:8,l) ,(aD,aD^dot)
   real(kind=prec), pointer  :: r(:), v(:,:), a(:), sigma(:), alpha(:), alphadot(:), pp(:,:,:)
   real(kind=prec), pointer :: t(:,:)
   integer :: at_numb, lmx, nr, nsp    !  atomic number, max_l, number. of mesh points, number of spin
   real(kind=prec) :: so_scale !strength of the SOC on this atom type
   integer :: label_len, nm
   real(kind=prec), pointer :: a_au(:) !!$ for this, the scaled 'a'
   integer, pointer :: idxdn(:)
   real(kind=prec), pointer :: xi0(:,:,:), xi1(:,:,:), xi2(:,:,:) ! SOC parameter (the radial integrals)
   real(kind=prec), pointer :: split(:)
   real(kind=prec) :: magmom 
   !!$ nonspherical stuff
   integer :: nonsph_lmx !the l quantum number of the highest order correction the potential
   real(kind=prec), pointer :: nonsph_phiVphi(:,:,:,:) !the radial integrals
   real(kind=prec), pointer :: nonsph_V(:,:) !the nonspherical corrections to the potential, specifically the radial part of the product V(r) Y_lm
   type (zcsrmat) :: zcsr_nonsph(3)
 end type tclass

 type tsite
!!$ Site, defined by its position and occupancy (atomic class)
   real(kind=prec) ::  coord(3)
   character(len=20) :: label  ! Atomic occupancy of a site
   integer :: iclass            
   type(tclass), pointer :: ptr
   integer, pointer :: idxdn_m(:),idxsh(:),idxlead(:)
 end type tsite

 Type t_atoms_set_EMTO
    Integer :: num, nspin
    real(kind=prec) :: fermi
    real(kind=prec), Pointer :: vmtz(:)
    Type (tclass), Pointer :: at (:)
 End Type t_atoms_set_EMTO

 type t_ctrl_opt
!!$ Additional options in the CTRL.omta file
   integer :: verbos
   real(kind=prec) :: facvol, ommax1(3), ommax2(3), mix
   logical :: fixlat, sclwsr, set_sigma
 end type t_ctrl_opt
 
 Type t_omta_struct
   Integer :: ns, neighm, ldim, idim
   integer        , allocatable :: npr (:)
   real(kind=prec), allocatable :: s   (:)
   integer        , allocatable :: iax(:,:,:)
   real(kind=prec) :: base(2,2),perp_trans(3)
 End Type t_omta_struct

 type t_omta_logder
   integer :: lidim
   real(kind=prec), pointer :: pph(:,:)
   type(zcsrmat) :: logder, hsoc
   type(zcsrmat) :: logderdot
 end type t_omta_logder

 real(kind=prec), pointer :: eband(:,:,:)

!!$ This module contains some basic matrix/vector operations and/or solutions.
!!$ Do NOT change anything before you completely understand every detail !!!

 CONTAINS

   integer function ll(ilm)
!- Returns l, given lm index
   implicit none
! Passed variables:
   integer ,intent(in) :: ilm
! Local variables:
   integer, parameter :: nlmax=20
   integer lla(0:(2*nlmax-1)**2)
! Data statements:
   data lla/01*-1,01*00,03*01,05*02,07*03,09*04,                     &
                  11*05,13*06,15*07,17*08,19*09,                     &
                  21*10,23*11,25*12,27*13,29*14,                     &
                  31*15,33*16,35*17,37*18,39*19,                     &
                  41*20,43*21,45*22,47*23,49*24,                     &
                  51*25,53*26,55*27,57*28,59*29,                     &
                  61*30,63*31,65*32,67*33,69*34,                     &
                  71*35,73*36,75*37,77*38/
   ll = lla(ilm)
   end function ll
  
   subroutine cross(x1,x2,cr12)
!- Cross product cr12 = X1 cross X2
   implicit none
! Passed variables:
   real(kind=prec) ,intent( in) :: x1(3),x2(3)
   real(kind=prec) ,intent(out) :: cr12(3)

   cr12(1) = x1(2)*x2(3) - x1(3)*x2(2)
   cr12(2) = x1(3)*x2(1) - x1(1)*x2(3)
   cr12(3) = x1(1)*x2(2) - x1(2)*x2(1)
   return
   end subroutine cross

   real(kind=prec)  function ddet33(matrix)
!- calculates the determinant of a 3X3 matrix
   implicit none
! Passed variables:
   real(kind=prec)  ,intent(in) :: matrix(3,3)
! Local variables:
   real(kind=prec)  :: m1cm2(3)

   call cross(matrix(1:3,2),matrix(1:3,3),m1cm2)
   ddet33 = DOT_PRODUCT(matrix(:,1),m1cm2)
   return
   end function ddet33

   subroutine dspfa(ap,n,kpvt,info)
   implicit none
   integer          ,intent(in   ) ::  n
   integer          ,intent(  out) ::  kpvt(n)
   integer          ,intent(  out) ::  info
   real(kind=prec)  ,intent(inout) ::  ap(n*(n+1)/2)
!!$dspfa factors a real(kind=prec)  symmetric matrix stored in
!!$packed form by elimination with symmetric pivoting.
!!$
!!$to solve  a*x = b , follow dspfa by dspsl.
!!$to compute  inverse(a)*c , follow dspfa by dspsl.
!!$to compute  determinant(a) , follow dspfa by dspdi.
!!$to compute  inertia(a) , follow dspfa by dspdi.
!!$to compute  inverse(a) , follow dspfa by dspdi.
!!$
!!$on entry
!!$
!!$   ap      real(kind=prec)  (n*(n+1)/2)
!!$           the packed form of a symmetric matrix  a .  the
!!$           columns of the upper triangle are stored sequentially
!!$           in a one-dimensional array of length  n*(n+1)/2 .
!!$           see comments below for details.
!!$
!!$   n       integer
!!$           the order of the matrix  a .
!!$
!!$output
!!$
!!$   ap      a block diagonal matrix and the multipliers which
!!$           were used to obtain it stored in packed form.
!!$           the factorization can be written  a = u*d*trans(u)
!!$           where  u  is a product of permutation and unit
!!$           upper triangular matrices , trans(u) is the
!!$           transpose of  u , and  d  is block diagonal
!!$           with 1 by 1 and 2 by 2 blocks.
!!$
!!$   kpvt    integer(n)
!!$           an integer vector of pivot indices.
!!$
!!$   info    integer
!!$           = 0  normal value.
!!$           = k  if the k-th pivot block is singular. this is
!!$                not an error condition for this subroutine,
!!$                but it does indicate that dspsl or dspdi may
!!$                divide by zero if called.
!!$
!!$packed storage
!!$
!!$     the following program segment will pack the upper
!!$     triangle of a symmetric matrix.
!!$
!!$           k = 0
!!$           do 20 j = 1, n
!!$              do 10 i = 1, j
!!$                 k = k + 1
!!$                 ap(k)  = a(i,j)
!!$        10    continue
!!$        20 continue
!!$
!!$linpack. this version dated 08/14/78 .
!!$james bunch, univ. calif. san diego, argonne nat. lab.
!!$
!!$subroutines and functions
!!$
!!$blas dswap
!!$fortran ABS,MAX,SQRT
!!$
!  Local variables:
   real(kind=prec)  ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
   real(kind=prec)  absakk,alpha,colmax,rowmax
   integer ij,ijl,l,ik,ikm1,im,imax,imaxp1,imim,imj,imk
   integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep,l1(1)
   logical swap

!  initialize
!  alpha is used in choosing pivot block size.
   alpha = (1.0d0 + SQRT(17.0d0))/8.0d0
   info = 0

!  main loop on k, which goes from n to 1.
   k = n
   ik = (n*(n - 1))/2
   do
! ---- leave the loop if k=0 or k=1.
     if (k  ==  0) exit
     if (k  ==  1) then
        kpvt(1) = 1
        if (ap(1)  ==  0.0d0) info = 1
        exit
     endif
! -- this section of code determines the kind of
! -- elimination to be performed.  when it is completed,
! -- kstep will be set to the size of the pivot block, and
! -- swap will be set to .true. if an interchange is required.
     km1 = k - 1
     kk = ik + k
     absakk = ABS(ap(kk))
! -- determine the largest off-diagonal element in column k.
     l1=MAXLOC(ABS(ap(ik+1:ik+k-1))); imax = l1(1)
     
     imk = ik + imax
     colmax = ABS(ap(imk))
     if (absakk  >=  alpha*colmax) then
       kstep = 1
       swap = .false.
     else
! -- determine the largest off-diagonal element in row imax.
       rowmax = 0.0d0
       imaxp1 = imax + 1
       im = imax*(imax - 1)/2
       imj = im + 2*imax
       do j = imaxp1, k
          rowmax = MAX(rowmax,ABS(ap(imj)))
          imj = imj + j
       enddo
       if (imax  /=  1) then
          l1=MAXLOC(ABS(ap(im+1:im+imax-1))); jmax = l1(1)
          jmim = jmax + im
          rowmax = MAX(rowmax,ABS(ap(jmim)))
       endif
       imim = imax + im
       if (ABS(ap(imim))  >=  alpha*rowmax) then
          kstep = 1
          swap = .true.
       elseif (absakk  >=  alpha*colmax*(colmax/rowmax)) then
            kstep = 1
            swap = .false.
       else
          kstep = 2
          swap = imax  /=  km1
       endif
     endif
     if (MAX(absakk,colmax)  ==  0.0d0) then
!! --- column k is zero.  set info and iterate the loop.
       kpvt(k) = k
       info = k
     elseif (kstep  ==  1) then
!! --- 1 x 1 pivot block.
       if (swap) then
!! ----- perform an interchange.
         call dswap(imax,ap(im+1),1,ap(ik+1),1)
         imj = ik + imax
         do jj = imax, k
           j = k + imax - jj
           jk = ik + j
           t = ap(jk)
           ap(jk) = ap(imj)
           ap(imj) = t
           imj = imj - (j - 1)
         enddo
       endif
!! --- perform the elimination.
       ij = ik - (k - 1)
       do jj = 1, km1
         j = k - jj
         jk = ik + j
         mulk = -ap(jk)/ap(kk)
         t = mulk
         do l=1,j
           ijl = ij + l
           ap(ijl)=ap(ijl) + t * ap(ik+l) 
         enddo 
         ap(jk) = mulk
         ij = ij - (j - 1)
       enddo
!! --- set the pivot array.
       kpvt(k) = k
       if (swap) kpvt(k) = imax
     else
!! --  2 x 2 pivot block.
       km1k = ik + k - 1
       ikm1 = ik - (k - 1)
       if (swap) then
!! ------ perform an interchange.
         call dswap(imax,ap(im+1),1,ap(ikm1+1),1)
         imj = ikm1 + imax
         do jj = imax, km1
           j = km1 + imax - jj
           jkm1 = ikm1 + j
           t = ap(jkm1)
           ap(jkm1) = ap(imj)
           ap(imj) = t
           imj = imj - (j - 1)
         enddo
         t = ap(km1k)
         ap(km1k) = ap(imk)
         ap(imk) = t
       endif
!! --- perform the elimination.
       km2 = k - 2
       if (km2  /=  0) then
         ak = ap(kk)/ap(km1k)
         km1km1 = ikm1 + k - 1
         akm1 = ap(km1km1)/ap(km1k)
         denom = 1.0d0 - ak*akm1
         ij = ik - (k - 1) - (k - 2)
         do jj = 1, km2
           j = km1 - jj
           jk = ik + j
           bk = ap(jk)/ap(km1k)
           jkm1 = ikm1 + j
           bkm1 = ap(jkm1)/ap(km1k)
           mulk = (akm1*bk - bkm1)/denom
           mulkm1 = (ak*bkm1 - bk)/denom
           
           t = mulk                   
           do l=1,j
             ijl = ij + l
             ap(ijl)=ap(ijl) + t * ap(ik+l) 
           enddo

           t = mulkm1
           do l=1,j
             ijl = ij + l
             ap(ijl)=ap(ijl) + t * ap(ikm1+l) 
           enddo

           ap(jk) = mulk
           ap(jkm1) = mulkm1
           ij = ij - (j - 1)
         enddo
       endif
!! --- set the pivot array.
       kpvt(k) = 1 - k
       if (swap) kpvt(k) = -imax
       kpvt(k-1) = kpvt(k)
     endif
     ik = ik - (k - 1)
     if (kstep  ==  2) ik = ik - (k - 2)
     k = k - kstep
   end do  
   return
   end subroutine dspfa 

   subroutine zspfa(ap,n,kpvt,info)
   implicit none
   integer           ,intent(in   ) ::  n
   integer           ,intent(  out) ::  kpvt(n)
   integer           ,intent(  out) ::  info
   complex(kind=prec),intent(inout) ::  ap(n*(n+1)/2)
!!$See descriptions in dspfa
!!$Local variables:
   complex(kind=prec) ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
   real(kind=prec)  :: absakk,colmax,rowmax
   integer ij,ijl,l,ik,ikm1,im,imax,imaxp1,imim,imj,imk
   integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep,l1(1)
   logical swap
   real(kind=prec) , parameter :: zero =0.d0
   real(kind=prec) , parameter :: alpha=(1.0d0 + SQRT(17.0d0))/8.0d0

!  initialize
   info = 0

!  main loop on k, which goes from n to 1.
   k = n
   ik = (n*(n - 1))/2
   do
! --- leave the loop if k=0 or k=1.
     if (k  ==  0) exit
     if (k  ==  1) then
       kpvt(1) = 1
       if (abs(ap(1))  ==  zero) info = 1
       exit
     endif
! -- this section of code determines the kind of
! -- elimination to be performed.  when it is completed,
! -- kstep will be set to the size of the pivot block, and
! -- swap will be set to .true. if an interchange is required.
     km1 = k - 1
     kk = ik + k
     absakk = ABS(ap(kk))
! -- determine the largest off-diagonal element in column k.
     l1=MAXLOC(ABS(ap(ik+1:ik+k-1))); imax = l1(1)
     imk = ik + imax
     colmax = ABS(ap(imk))
     if (absakk  >=  alpha*colmax) then
       kstep = 1
       swap = .false.
     else
! ---- determine the largest off-diagonal element in row imax.
       rowmax = zero
       imaxp1 = imax + 1
       im = imax*(imax - 1)/2
       imj = im + 2*imax
       do j = imaxp1, k
         rowmax = MAX(rowmax,ABS(ap(imj)))
         imj = imj + j
       enddo
       if (imax  /=  1) then
         l1=MAXLOC(ABS(ap(im+1:im+imax-1))); jmax = l1(1)
         jmim = jmax + im
         rowmax = MAX(rowmax,ABS(ap(jmim)))
       endif
       imim = imax + im
       if (ABS(ap(imim))  >=  alpha*rowmax) then
         kstep = 1
         swap = .true.
       elseif (absakk  >=  alpha*colmax*(colmax/rowmax)) then
         kstep = 1
         swap = .false.
       else
         kstep = 2
         swap = imax  /=  km1
       endif
     endif
     if (MAX(absakk,colmax)  ==  zero) then
! ---- column k is zero.  set info and iterate the loop.
       kpvt(k) = k
       info = k
     elseif (kstep  ==  1) then
! ---- 1 x 1 pivot block.
       if (swap) then
! ------ perform an interchange.
         call zswap(imax,ap(im+1),1,ap(ik+1),1)
         imj = ik + imax
         do jj = imax, k
           j = k + imax - jj
           jk = ik + j
           t = ap(jk)
           ap(jk) = ap(imj)
           ap(imj) = t
           imj = imj - (j - 1)
         enddo
       endif
! ---- perform the elimination.
       ij = ik - (k - 1)
       do jj = 1, km1
         j = k - jj
         jk = ik + j
         mulk = -ap(jk)/ap(kk)
         t = mulk
         do l=1,j
           ijl = ij + l
           ap(ijl)=ap(ijl) + t * ap(ik+l) 
         enddo 
         ap(jk) = mulk
         ij = ij - (j - 1)
       enddo
! ---- set the pivot array.
       kpvt(k) = k
       if (swap) kpvt(k) = imax
     else
! ---  2 x 2 pivot block.
       km1k = ik + k - 1
       ikm1 = ik - (k - 1)
       if (swap) then
! ------ perform an interchange.
         call zswap(imax,ap(im+1),1,ap(ikm1+1),1)
         imj = ikm1 + imax
         do jj = imax, km1
           j = km1 + imax - jj
           jkm1 = ikm1 + j
           t = ap(jkm1)
           ap(jkm1) = ap(imj)
           ap(imj) = t
           imj = imj - (j - 1)
         enddo
         t = ap(km1k)
         ap(km1k) = ap(imk)
         ap(imk) = t
       endif
! ---- perform the elimination.
       km2 = k - 2
       if (km2  /=  0) then
         ak = ap(kk)/ap(km1k)
         km1km1 = ikm1 + k - 1
         akm1 = ap(km1km1)/ap(km1k)
         denom = DCMPLX(1.0d0,0.d0) - ak*akm1
         ij = ik - (k - 1) - (k - 2)
         do jj = 1, km2
           j = km1 - jj
           jk = ik + j
           bk = ap(jk)/ap(km1k)
           jkm1 = ikm1 + j
           bkm1 = ap(jkm1)/ap(km1k)
           mulk = (akm1*bk - bkm1)/denom
           mulkm1 = (ak*bkm1 - bk)/denom
           
           t = mulk                   
           do l=1,j
             ijl = ij + l
             ap(ijl)=ap(ijl) + t * ap(ik+l) 
           enddo

           t = mulkm1
           do l=1,j
             ijl = ij + l
             ap(ijl)=ap(ijl) + t * ap(ikm1+l) 
           enddo

           ap(jk) = mulk
           ap(jkm1) = mulkm1
           ij = ij - (j - 1)
         enddo
       endif
! ---- set the pivot array.
       kpvt(k) = 1 - k
       if (swap) kpvt(k) = -imax
       kpvt(k-1) = kpvt(k)
     endif
     ik = ik - (k - 1)
     if (kstep  ==  2) ik = ik - (k - 2)
     k = k - kstep
   end do  
   return
   end subroutine zspfa 

   subroutine dspsl(ap,n,kpvt,b)
   implicit none
   integer          ,intent(in   ) :: n
   integer          ,intent(in   ) :: kpvt(n)
   real(kind=prec)  ,intent(in   ) :: ap(n*(n+1)/2)
   real(kind=prec)  ,intent(inout) :: b(n)
!!$dsisl solves the real(kind=prec)  symmetric system
!!$a * x = b
!!$using the factors computed by dspfa.
!!$
!!$on entry
!!$
!!$   ap      real(kind=prec) (n*(n+1)/2)
!!$           the output from dspfa.
!!$
!!$   n       integer
!!$           the order of the matrix  a .
!!$
!!$   kpvt    integer(n)
!!$           the pivot vector from dspfa.
!!$
!!$   b       real(kind=prec) (n)
!!$           the right hand side vector.
!!$
!!$on return
!!$
!!$   b       the solution vector  x .
!!$
!!$error condition
!!$
!!$   a division by zero may occur if  dspco  has set rcond  ==  0.0
!!$   or  dspfa  has set info  /=  0  .
!!$
!!$to compute  inverse(a) * c  where  c  is a matrix
!!$with  p  columns
!!$      call dspfa(ap,n,kpvt,info)
!!$      if (info  /=  0) goto ...
!!$      do j = 1, p
!!$         call dspsl(ap,n,kpvt,c(:,j))
!!$      enddo
!!$
!!$linpack. this version dated 08/14/78 .
!!$james bunch, univ. calif. san diego, argonne nat. lab.

!!$Local variables:
   real(kind=prec)  :: ak,akm1,bk,bkm1,denom,temp,t,t1,br
   integer          :: ik,ikm1,ikp1,k,kk,km1k,km1km1,kp,l

!  loop backward applying the transformations and d inverse to b.
   k = n
   ik = (n*(n - 1))/2

   do while (k  /=  0) 
     kk = ik + k
     if (kpvt(k)  >=  0) then
! ---- 1 x 1 pivot block.
       if (k  /=  1) then
         kp = kpvt(k)
         if (kp  /=  k) then
! -------- interchange.
           temp = b(k)
           b(k) = b(kp)
           b(kp) = temp
         endif
!        apply the transformation.
         t=b(k)
         do l = 1,k-1
           b(l) = b(l) + t*ap(ik+l)
         enddo
       endif
! ---- apply d inverse.
       b(k) = b(k)/ap(kk)
       k = k - 1
       ik = ik - k
     else
! ---- 2 x 2 pivot block.
       ikm1 = ik - (k - 1)
       if (k  /=  2) then
         kp = ABS(kpvt(k))
         if (kp  /=  k-1) then
! -------- interchange.
           temp = b(k-1)
           b(k-1) = b(kp)
           b(kp)  = temp
         endif
!        apply the transformation.

         t =b(k  )
         t1=b(k-1)
         do l = 1, k-2
           b(l) = b(l) + t*ap(ik+l) + t1*ap(ikm1+l)  
         enddo
       endif
!      apply d inverse.
       km1k = ik + k - 1
       kk = ik + k
       ak = ap(kk)/ap(km1k)
       km1km1 = ikm1 + k - 1
       akm1 = ap(km1km1)/ap(km1k)
       bk = b(k)/ap(km1k)
       bkm1 = b(k-1)/ap(km1k)
       denom = ak*akm1 - 1.0d0
       b(k) = (akm1*bk - bkm1)/denom
       b(k-1) = (ak*bkm1 - bk)/denom
       k = k - 2
       ik = ik - (k + 1) - k
     endif
   enddo

!  loop forward applying the transformations.
   k = 1
   ik = 0
   
   do while (k  <=  n)
     if (kpvt(k)  >=  0) then
!      1 x 1 pivot block.
       if (k  /=  1) then
!        apply the transformation.
         br = 0.d0
         do l = 1, k-1
           br = br + ap(ik+l)*b(l)
         enddo
         b(k) = b(k) + br
         kp = kpvt(k)
         if (kp  /=  k) then  ! interchange.
           temp = b(k)
           b(k) = b(kp)
           b(kp) = temp
         endif
       endif
       ik = ik + k
       k = k + 1
     else
!      2 x 2 pivot block.
       if (k  /=  1) then
!        apply the transformation.
         br = 0.d0
         do l = 1, k-1
           br = br + ap(ik+l)*b(l)
         enddo
         b(k) = b(k) + br
 
         ikp1 = ik + k
 
         br = 0.d0
         do l = 1, k-1
           br = br + ap(ikp1+l)*b(l)
         enddo
         b(k+1) = b(k+1) + br
         kp = ABS(kpvt(k))
         if (kp  /=  k) then  ! interchange.
           temp = b(k)
           b(k) = b(kp)
           b(kp) = temp
         endif
       endif
       ik = ik + k + k + 1
       k = k + 2
     endif
   enddo
   return
   end subroutine dspsl 

   subroutine zspsl(ap,n,kpvt,b0)
   implicit none
   integer          ,intent(in   ) :: n
   integer          ,intent(in   ) :: kpvt(n)
   complex(kind=prec),intent(in   ) :: ap(n*(n+1)/2)
   real(kind=prec)       ,intent(inout) :: b0(n)
!!$dsisl solves the real(kind=prec)  symmetric system
!!$a * x = b
!!$using the factors computed by dspfa.
!!$
!!$on entry
!!$
!!$   ap      real(kind=prec) (n*(n+1)/2)
!!$           the output from dspfa.
!!$
!!$   n       integer
!!$           the order of the matrix  a .
!!$
!!$   kpvt    integer(n)
!!$           the pivot vector from dspfa.
!!$
!!$   b       real(kind=prec) (n)
!!$           the right hand side vector.
!!$
!!$on return
!!$
!!$   b       the solution vector  x .
!!$
!!$error condition
!!$
!!$   a division by zero may occur if  dspco  has set rcond  ==  0.0
!!$   or  dspfa  has set info  /=  0  .
!!$
!!$to compute  inverse(a) * c  where  c  is a matrix
!!$with  p  columns
!!$      call dspfa(ap,n,kpvt,info)
!!$      if (info  /=  0) goto ...
!!$      do j = 1, p
!!$         call dspsl(ap,n,kpvt,c(:,j))
!!$      enddo
!!$
!!$linpack. this version dated 08/14/78 .
!!$james bunch, univ. calif. san diego, argonne nat. lab.
!!$
!!$Local variables:
   complex(kind=prec)   :: ak,akm1,bk,bkm1,denom,temp,t,t1,br
   real(kind=prec) , parameter :: zero=0.d0
   complex(kind=prec),allocatable :: b(:)
   integer          :: ik,ikm1,ikp1,k,kk,km1k,km1km1,kp,l

!  loop backward applying the transformations and d inverse to b.
   allocate(b(n))
   b=b0
   k = n
   ik = (n*(n - 1))/2

   do while (k  /=  0) 
     kk = ik + k
     if (kpvt(k)  >=  0) then
! ---- 1 x 1 pivot block.
       if (k  /=  1) then
         kp = kpvt(k)
         if (kp  /=  k) then
! -------- interchange.
           temp = b(k)
           b(k) = b(kp)
           b(kp) = temp
         endif
!        apply the transformation.
         t=b(k)
         do l = 1,k-1
           b(l) = b(l) + t*ap(ik+l)
         enddo
       endif
! ---- apply d inverse.
       b(k) = b(k)/ap(kk)
       k = k - 1
       ik = ik - k
     else
! ---- 2 x 2 pivot block.
       ikm1 = ik - (k - 1)
       if (k  /=  2) then
         kp = ABS(kpvt(k))
         if (kp  /=  k-1) then
! -------- interchange.
           temp = b(k-1)
           b(k-1) = b(kp)
           b(kp)  = temp
         endif
!        apply the transformation.
         t =b(k  )
         t1=b(k-1)
         do l = 1, k-2
           b(l) = b(l) + t*ap(ik+l) + t1*ap(ikm1+l)  
         enddo
       endif
!      apply d inverse.
       km1k = ik + k - 1
       kk = ik + k
       ak = ap(kk)/ap(km1k)
       km1km1 = ikm1 + k - 1
       akm1 = ap(km1km1)/ap(km1k)
       bk = b(k)/ap(km1k)
       bkm1 = b(k-1)/ap(km1k)
       denom = ak*akm1 - DCMPLX(1.0d0,0.d0)
       b(k) = (akm1*bk - bkm1)/denom
       b(k-1) = (ak*bkm1 - bk)/denom
       k = k - 2
       ik = ik - (k + 1) - k
     endif
   enddo

!  loop forward applying the transformations.
   k = 1
   ik = 0
   
   do while (k  <=  n)
     if (kpvt(k)  >=  0) then
!      1 x 1 pivot block.
       if (k  /=  1) then
!        apply the transformation.
         br = zero
         do l = 1, k-1
           br = br + ap(ik+l)*b(l)
         enddo
         b(k) = b(k) + br
         kp = kpvt(k)
         if (kp  /=  k) then  ! interchange.
           temp = b(k)
           b(k) = b(kp)
           b(kp) = temp
         endif
       endif
       ik = ik + k
       k = k + 1
     else
!      2 x 2 pivot block.
       if (k  /=  1) then
!        apply the transformation.
         br = zero
         do l = 1, k-1
           br = br + ap(ik+l)*b(l)
         enddo
         b(k) = b(k) + br
 
         ikp1 = ik + k
         br = zero
         do l = 1, k-1
           br = br + ap(ikp1+l)*b(l)
         enddo
         b(k+1) = b(k+1) + br
         kp = ABS(kpvt(k))
         if (kp  /=  k) then  ! interchange.
           temp = b(k)
           b(k) = b(kp)
           b(kp) = temp
         endif
       endif
       ik = ik + k + k + 1
       k = k + 2
     endif
   enddo
   
   b0=dreal(b)
   deallocate(b)
   return
   end subroutine zspsl 

   subroutine dswap (n,dx,incx,dy,incy)
!-Interchanges two vectors
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :lenght of dx and dy
!io  dx    :vector
!i   incx  :incrementation for x
!io  dy    :vector
!i   incy  :incrementation for y
!o Outputs:
!io  dx    :vector
!io  dy    :vector
!r Remarks:
!r Adapted from:  jack dongarra, linpack, 3/11/78.
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer incx,incy,n
   real(kind=prec)  dx(*),dy(*)
! Local variables:
   real(kind=prec) ,allocatable:: dtemp(:)

   if(n <= 0)return
   allocate(dtemp(n))
   if(incx /= 1.or.incy /= 1) then !increments not equal to 1
     dtemp(:)=dx(1:incx*n:incx)
     dx(1:incx*n:incx)=dy(1:incy*n:incy)
     dy(1:incy*n:incy)=dtemp(:)
   else !both increments equal to 1
     dtemp(:)=dx(1:n)
     dx(1:n)=dy(1:n)
     dy(1:n)=dtemp(:)
   endif
   deallocate(dtemp)
   return
   end subroutine dswap

   subroutine  zswap (n,zx,incx,zy,incy)
!  interchanges two vectors.
!  jack dongarra, 3/11/78.
!  modified 12/3/93, array(1) declarations changed to array(*)
   complex(kind=prec) :: zx(*),zy(*),ztemp
   integer i,incx,incy,ix,iy,n

   if(n.le.0)return
   if(incx.eq.1.and.incy.eq.1)go to 20

!  code for unequal increments or equal increments not equal to 1

   ix = 1
   iy = 1
   if(incx.lt.0)ix = (-n+1)*incx + 1
   if(incy.lt.0)iy = (-n+1)*incy + 1
   do 10 i = 1,n
     ztemp = zx(ix)
     zx(ix) = zy(iy)
     zy(iy) = ztemp
     ix = ix + incx
     iy = iy + incy
10 continue
   return

!  code for both increments equal to 1
20 do 30 i = 1,n
     ztemp = zx(i)
     zx(i) = zy(i)
     zy(i) = ztemp
30 continue
   return
   end subroutine zswap

 end module omta_defs
