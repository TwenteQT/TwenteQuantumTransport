!!$ $Id: sparselib.F90 1909 2013-06-17 15:48:26Z antst $
#include "math_def.h"

Module sparselib
   Implicit None

   Type zdensemat
!!$    Auxiliary type (complx. matrix) to be used while
!!$    eg. allocating the vector of matrices of different sizes.
      Complex(Kind=DEF_DBL_PREC), Allocatable :: bl(:, :)
      Integer :: alloc = 0, nrow = 0, ncol = 0
   End Type zdensemat

   Type rdensemat
!!$    Auxiliary type (real. matrix) to be used while
!!$    eg. allocating the vector of matrices of different sizes.
      Real(Kind=DEF_DBL_PREC), Allocatable :: bl(:, :)
      Integer :: alloc = 0, nrow = 0, ncol = 0
   End Type rdensemat

   Type rvec
!!$    Auxiliary type (real vector) to be used while
!!$    eg. allocating the vector of matrices of different sizes.
      Real(Kind=DEF_DBL_PREC), Allocatable :: vec(:)
      Integer :: size
   End Type rvec

   Type matblk
      Integer :: fcol, frow
      Complex(Kind=DEF_DBL_PREC), Allocatable :: bl(:, :)
   End Type matblk

   Type rdiagmat
      Real(Kind=DEF_DBL_PREC), Allocatable :: el(:)
      Integer :: ncol, nrow, alloc = 0
   End Type rdiagmat

   Type zcoomat
      Integer :: nnz, ncol, nrow, nnzmax, alloc = 0
      Complex(Kind=DEF_DBL_PREC), Allocatable :: a(:)
      Integer, Allocatable :: icol(:), irow(:)
   End Type zcoomat

   Type rcoomat
      Integer :: nnz, ncol, nrow, nnzmax, alloc = 0
      Real(Kind=DEF_DBL_PREC), Allocatable :: a(:)
      Integer, Allocatable :: icol(:), irow(:)
   End Type rcoomat

   Type zdiagmat
      Complex(Kind=DEF_DBL_PREC), Allocatable :: el(:)
      Integer :: ncol, nrow, alloc = 0
   End Type zdiagmat

   Type zcsrmat
      Integer :: nnz, ncol, nnzmax, nrow, alloc = 0
      Complex(Kind=DEF_DBL_PREC), Allocatable :: a(:)
      Integer, Allocatable :: ir(:), jc(:)
   End Type zcsrmat

   Type rcsrmat
      Integer :: nnz, ncol, nnzmax, nrow, alloc = 0
      Real(Kind=DEF_DBL_PREC), Allocatable :: a(:)
      Integer, Allocatable :: ir(:), jc(:)
   End Type rcsrmat

   Type rmatblk
      Integer :: nrow, ncol
      Real(Kind=DEF_DBL_PREC), Allocatable :: bl(:, :)
   End Type rmatblk

   Type zmatblk
      Integer :: nrow, ncol
      Complex(Kind=DEF_DBL_PREC), Allocatable :: bl(:, :)
   End Type zmatblk

   Type rblockmat
      Integer :: ncol, nrow, nblocks, alloc = 0, maxblocks
      Type(rmatblk), Allocatable :: a(:)
      Integer, Allocatable :: startrow(:), startcol(:), endrow(:), endcol(:)
   End Type rblockmat

   Type zblockmat
      Integer :: ncol, nrow, nblocks, alloc = 0, maxblocks
      Type(rmatblk), Allocatable :: a(:)
      Integer, Allocatable :: startrow(:), startcol(:), endrow(:), endcol(:)
   End Type zblockmat

   Interface init
      Module Procedure zinitcsr, rinitcsr
   End Interface init

   Interface free
      Module Procedure rfreecsr, zfreecsr, rfreediag, zfreediag, rfreecoo, zfreecoo, zfreedense !,rfreeblockmat,zfreeblockmat
   End Interface

   Interface alloc
      Module Procedure ralloccsr, zalloccsr, rallocdiag, zallocdiag, zalloccoo, ralloccoo, zallocdense, &
     & rallocdense !,rallocblockmat,zallocblockmat
   End Interface

   Interface mrealloc
      Module Procedure rrealloccsr, zrealloccsr
   End Interface

   Interface getsubmat
      Module Procedure zgetcsrsubmat, rgetcsrsubmat
!!$            rgetdiagsumbat
   End Interface

   Interface getrowind
      Module Procedure zgetrowind, rgetrowind
   End Interface

   Interface inspmatmul
      Module Procedure mmcsrdiag_zz, mmcsrdiag_zr, mmcsrdiag_rr, mmdiagcsr_zz, mmdiagcsr_rz, mmdiagcsr_rr
!!$     $         fmmdensecsr_zz,fmmdensecsr_zr,fmmdensecsr_rz,fmmdensecsr_rr,
!!$     $     mmcsrcsr_zz,mmcsrcsr_zr,mmcsrcsr_rz,mmcsrcsr_rr
   End Interface

   Interface spmatmul
      Module Procedure fmmcsrdiag_zz, &
     & fmmcsrdiag_zr, fmmcsrdiag_rz, fmmcsrdiag_rr, fmmdiagcsr_zz, fmmdiagcsr_zr, fmmdiagcsr_rz, &
     & fmmdiagcsr_rr, fmmcsrcsr_zz, fmmcsrcsr_zr, fmmcsrcsr_rz
!!$     fmmcsrdense_zz, fmmcsrdense_zr, fmmcsrdense_rz, fmmcsrdense_rr, fmmcsrvector_zz
   End Interface

   Interface spmatmul2
      Module Procedure mmcsrdense_zz, mmcsrdense_zz_base, mmcsrdense_zz_s, mmcsrvector_zz_s
   End Interface

   Interface spmatmul3
      Module Procedure fmmcsrdiagcsr_zrz
   End Interface

   Interface spherm
      Module Procedure csrherm_z, csrherm_r
   End Interface

   Interface transp
      Module Procedure csrtransp_z
   End Interface

   Interface spmatadd
      Module Procedure fcsraddcsr_zz, fcsraddcsr_zr
   End Interface

!!$    interface spadddiag
!!$       module procedure fcsradddiag_zz,fcsradddiag_zr,fcsradddiag_rr,fcsradddiag_rz
!!$    end interface

   Interface spshiftdiag
      Module Procedure zrcsrshiftdiag, zzcsrshiftdiag
   End Interface

   Interface inspmatadd
      Module Procedure csradddiag_zz, csradddiag_zr, csradddiag_rr, denseaddcsr_zz
   End Interface

   Interface spcopy
      Module Procedure spcopy_csr_z, spcopy_csr_r, spcopy_csr_rz, spcopy_dense_z
   End Interface

   Interface dump_matrix
      Module Procedure dumpzcsr, dumprcsr, dumpzdiag, dumprdiag, dumpzdense, dumprdense, dumprvector, &
     & dumpzvector, dumpzcoo, dumpzdensemat
   End Interface

   Interface coo2csr
      Module Procedure z_coo2csr, r_coo2csr
   End Interface

   Interface conv2csr
      Module Procedure z_conv2csr, r_conv2csr
   End Interface

   Interface ordercsr
      Module Procedure z_ordercsr, r_ordercsr
   End Interface

   Interface cooaddelms
      Module Procedure z_cooaddel, r_cooaddel
   End Interface

   Interface sptofull
      Module Procedure zcsr_sptofull
   End Interface

   Interface sp_add_dupes
      Module Procedure zcsr_add_dupes
   End Interface

   Type tspblpos
      Integer, Pointer :: pos(:, :)
   End Type tspblpos

   Integer, Parameter :: dumping_unit = 12234

Contains

   Subroutine zalloccoo(matr, maxnnz, nrow, ncol_in)
      Implicit None
      Type(zcoomat) :: matr
      Integer, Intent(In) :: nrow, maxnnz
      Integer, Optional :: ncol_in
      Integer :: ncol

      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in

      If (matr%alloc /= 0) Then
         If ((matr%ncol /= ncol) .Or. (matr%nnzmax /= maxnnz)) Then
            Deallocate (matr%a)
            Deallocate (matr%irow)
            Deallocate (matr%icol)
            matr%alloc = 0
            matr%nnzmax = 0
         End If
      End If
      If (matr%alloc == 0) Then
         Allocate (matr%a(maxnnz))
         Allocate (matr%irow(maxnnz))
         Allocate (matr%icol(maxnnz))
      End If
      matr%a = 0.0d0
      matr%irow = 0
      matr%icol = 0
      matr%alloc = 1
      matr%nnz = 0
      matr%nnzmax = maxnnz
      matr%ncol = ncol
      matr%nrow = nrow
   End Subroutine zalloccoo

   Subroutine ralloccoo(matr, maxnnz, nrow, ncol_in)
      Implicit None
      Type(rcoomat) :: matr
      Integer, Intent(In) :: nrow, maxnnz
      Integer, Optional :: ncol_in
      Integer :: ncol, mnnz
      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in
      mnnz = 0
      If (matr%alloc /= 0) Then
         If ((matr%ncol /= ncol) .Or. (matr%nnzmax < maxnnz)) Then
            Deallocate (matr%a)
            Deallocate (matr%irow)
            Deallocate (matr%icol)
            matr%alloc = 0
            matr%nnzmax = 0
         Else
            mnnz = matr%nnzmax
         End If

      End If
      If (matr%alloc == 0) Then
         Allocate (matr%a(maxnnz))
         Allocate (matr%irow(maxnnz))
         Allocate (matr%icol(maxnnz))
      End If
      matr%a = 0.0d0
      matr%irow = 0
      matr%icol = 0
      matr%alloc = 1
      matr%nnz = 0
      matr%nnzmax = Max(maxnnz, mnnz)
      matr%ncol = ncol
      matr%nrow = nrow
   End Subroutine ralloccoo

   Subroutine zinitcsr(matr)
      Implicit None
      Type(zcsrmat) :: matr
      matr%alloc = 0
!!$    matr%nrow=0
!!$    matr%ncol=0
!!$    matr%nnzmax=0
!!$    matr%nnz=0
!!$    nullify(matr%a,matr%ir,matr%jc)
   End Subroutine zinitcsr

   Subroutine rinitcsr(matr)
      Implicit None
      Type(rcsrmat) :: matr
      matr%alloc = 0
!!$    matr%nrow=0
!!$    matr%ncol=0
!!$    matr%nnz=0
   End Subroutine rinitcsr

   Subroutine zalloccsr(matr, maxnnz, nrow, ncol_in)
      Implicit None
      Type(zcsrmat) :: matr
      Integer, Intent(In) :: nrow, maxnnz
      Integer, Optional :: ncol_in
      Integer :: ncol, mnnz
      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in
      mnnz = 0
      If (matr%alloc /= 0) Then
         If ((matr%nrow /= nrow) .Or. (matr%nnzmax < maxnnz)) Then
          !!$write(*,*) 'deall'
            Deallocate (matr%a)
            Deallocate (matr%ir)
            Deallocate (matr%jc)
            matr%alloc = 0
            matr%nnzmax = 0
            matr%nnz = 0
         Else
            mnnz = matr%nnzmax
         End If
      End If
      If (matr%alloc == 0) Then
         Allocate (matr%a(maxnnz))
         Allocate (matr%ir(nrow + 1))
         Allocate (matr%jc(maxnnz))
      End If
!!$    matr%a = 0.0d0
      matr%ir = 1
      matr%jc = 0
      matr%alloc = 1
      matr%nnz = 0
      matr%nnzmax = Max(maxnnz, mnnz)
      matr%ncol = ncol
      matr%nrow = nrow
   End Subroutine zalloccsr

   Subroutine ralloccsr(matr, maxnnz, nrow, ncol_in)
      Implicit None
      Type(rcsrmat) :: matr
      Integer, Intent(In) :: nrow, maxnnz
      Integer, Optional :: ncol_in
      Integer :: ncol
      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in

      If (matr%alloc /= 0) Then
         If ((matr%nrow /= nrow) .Or. (matr%nnzmax /= maxnnz)) Then
            Deallocate (matr%a)
            Deallocate (matr%ir)
            Deallocate (matr%jc)
            matr%alloc = 0
            matr%nnzmax = 0
         End If
      End If
      If (matr%alloc == 0) Then
         Allocate (matr%a(maxnnz))
         Allocate (matr%ir(nrow + 1))
         Allocate (matr%jc(maxnnz))
      End If
      matr%a = 0.0d0
      matr%ir = 1
      matr%jc = 0
      matr%alloc = 1
      matr%nnz = 0
      matr%nnzmax = maxnnz
      matr%ncol = ncol
      matr%nrow = nrow
   End Subroutine ralloccsr

   Subroutine zrealloccsr(matr, maxnnz, nrow_in, ncol_in)
      Implicit None
      Type(zcsrmat) :: matr
      Integer, Intent(In) :: maxnnz
      Integer, Optional :: nrow_in, ncol_in
      Complex(Kind=DEF_DBL_PREC), Allocatable :: b(:)
      Integer, Allocatable :: jc1(:), ir1(:)
      Integer :: mnr
      If (matr%alloc /= 0) Then
         If (maxnnz /= matr%nnzmax) Then
            matr%nnz = Min(maxnnz, matr%nnz)
            Allocate (b(maxnnz))
            Allocate (jc1(maxnnz))
            b(1:matr%nnz) = matr%a(1:matr%nnz)
            jc1(1:matr%nnz) = matr%jc(1:matr%nnz)
            Deallocate (matr%a)
            Deallocate (matr%jc)
            matr%a = b
            matr%jc = jc1
            matr%nnzmax = maxnnz
         End If
         If (present(nrow_in)) Then
            If (nrow_in /= matr%nrow) Then
               mnr = Min(nrow_in, matr%nrow) + 1
               Allocate (ir1(nrow_in + 1))
               ir1 = 1
               ir1(1:mnr) = matr%ir(1:mnr)
               Deallocate (matr%ir)
               matr%ir = ir1
               matr%nrow = nrow_in
            End If
         End If
         If (present(ncol_in)) matr%ncol = ncol_in
      Else
         If (present(nrow_in)) Then
!!$                mnr = nrow_in
!!$                If (present(ncol_in)) mnr = ncol_in
!!$                Call alloc (matr, maxnnz, nrow_in, mnr)
            Call alloc(matr, maxnnz, nrow_in, ncol_in)
         Else
            Write (*, *) 'you are in troubles....'
         End If
      End If
   End Subroutine zrealloccsr

   Subroutine rrealloccsr(matr, maxnnz, nrow_in, ncol_in)
      Implicit None
      Type(rcsrmat) :: matr
      Integer, Intent(In) :: maxnnz
      Integer, Optional :: nrow_in, ncol_in
      Real(Kind=DEF_DBL_PREC), Allocatable :: b(:)
      Integer, Allocatable :: jc1(:), ir1(:)
      Integer :: mnr
      If (matr%alloc /= 0) Then
         If (maxnnz /= matr%nnzmax) Then
            matr%nnz = Min(maxnnz, matr%nnz)
            Allocate (b(maxnnz))
            Allocate (jc1(maxnnz))
            b(1:matr%nnz) = matr%a(1:matr%nnz)
            jc1(1:matr%nnz) = matr%jc(1:matr%nnz)
            Deallocate (matr%a)
            Deallocate (matr%jc)
            matr%a = b
            matr%jc = jc1
            matr%nnzmax = maxnnz
         End If
         If (present(nrow_in)) Then
            If (nrow_in /= matr%nrow) Then
               mnr = Min(nrow_in, matr%nrow) + 1
               Allocate (ir1(nrow_in + 1))
               ir1 = 1
               ir1(1:mnr) = matr%ir(1:mnr)
               Deallocate (matr%ir)
               matr%ir = ir1
               matr%nrow = nrow_in
            End If
         End If
         If (present(ncol_in)) matr%ncol = ncol_in
      Else
         If (present(nrow_in)) Then
!!$                mnr = nrow_in
!!$                If (present(ncol_in)) mnr = ncol_in
!!$                Call alloc (matr, maxnnz, mnr)
            Call alloc(matr, maxnnz, nrow_in, ncol_in)
         Else
            Write (*, *) 'you are in troubles....'
         End If
      End If
   End Subroutine rrealloccsr

   Subroutine zfreecoo(matr)
      Implicit None
      Type(zcoomat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%a)
         Deallocate (matr%irow)
         Deallocate (matr%icol)
         matr%alloc = 0
         matr%nnzmax = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine zfreecoo

   Subroutine rfreecoo(matr)
      Implicit None
      Type(rcoomat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%a)
         Deallocate (matr%irow)
         Deallocate (matr%icol)
         matr%alloc = 0
         matr%nnzmax = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine rfreecoo

   Subroutine zfreecsr(matr)
      Implicit None
      Type(zcsrmat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%a)
         Deallocate (matr%ir)
         Deallocate (matr%jc)
         matr%alloc = 0
         matr%nnzmax = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine zfreecsr

   Subroutine rfreecsr(matr)
      Implicit None
      Type(rcsrmat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%a)
         Deallocate (matr%ir)
         Deallocate (matr%jc)
         matr%alloc = 0
         matr%nnzmax = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine rfreecsr

   Subroutine zallocdense(mat, nrow, ncol_in)
      Implicit None
      Type(zdensemat) :: mat
      Integer, Intent(In) :: nrow
      Integer, Optional :: ncol_in
      Integer :: ncol

      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in

      If ((mat%alloc /= 0) .And. ((mat%ncol /= ncol) .Or. (mat%nrow /= nrow))) Then
         Deallocate (mat%bl)
         mat%alloc = 0
      End If
      If (mat%alloc == 0) Then
         Allocate (mat%bl(nrow, ncol))
      End If
      mat%bl = dcmplx(0.0d0, 0.0d0)
      mat%ncol = ncol
      mat%nrow = nrow
      mat%alloc = 1
   End Subroutine zallocdense

   Subroutine rallocdense(mat, nrow, ncol_in)
      Implicit None
      Type(rdensemat) :: mat
      Integer, Intent(In) :: nrow
      Integer, Optional :: ncol_in
      Integer :: ncol
      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in
      If ((mat%alloc /= 0) .And. ((mat%ncol /= ncol) .Or. (mat%nrow /= nrow))) Then
         Deallocate (mat%bl)
         mat%alloc = 0
      End If
      If (mat%alloc == 0) Then
         Allocate (mat%bl(nrow, ncol))
      End If
      mat%bl = 0.0d0
      mat%ncol = ncol
      mat%nrow = nrow
      mat%alloc = 1
   End Subroutine rallocdense

   Subroutine zallocdiag(mat, nrow, ncol_in)
      Implicit None
      Type(zdiagmat) :: mat
      Integer, Intent(In) :: nrow
      Integer, Optional :: ncol_in
      Integer :: ncol
      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in
      If ((mat%alloc /= 0) .And. (mat%ncol /= ncol)) Then
         Deallocate (mat%el)
         mat%alloc = 0
      End If

      If (mat%alloc == 0) Then
         Allocate (mat%el(ncol))
      End If
      mat%el = 0.0d0
      mat%ncol = ncol
      mat%nrow = nrow
      mat%alloc = 1
   End Subroutine zallocdiag

   Subroutine rallocdiag(mat, nrow, ncol_in)
      Implicit None
      Type(rdiagmat) :: mat
      Integer, Intent(In) :: nrow
      Integer, Optional :: ncol_in
      Integer :: ncol
      ncol = nrow
      If (present(ncol_in)) ncol = ncol_in
      If ((mat%alloc /= 0) .And. (mat%ncol /= ncol)) Then
         Deallocate (mat%el)
         mat%alloc = 0
      End If

      If (mat%alloc == 0) Then
         Allocate (mat%el(ncol))
      End If
      mat%el = 0.0d0
      mat%ncol = ncol
      mat%nrow = nrow
      mat%alloc = 1
   End Subroutine rallocdiag

   Subroutine zfreedense(matr)
      Implicit None
      Type(zdensemat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%bl)
         matr%alloc = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine zfreedense

   Subroutine rfreedense(matr)
      Implicit None
      Type(rdensemat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%bl)
         matr%alloc = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine rfreedense

   Subroutine zfreediag(matr)
      Implicit None
      Type(zdiagmat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%el)
         matr%alloc = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine zfreediag

   Subroutine rfreediag(matr)
      Implicit None
      Type(rdiagmat) :: matr
      If (matr%alloc /= 0) Then
         Deallocate (matr%el)
         matr%alloc = 0
         matr%ncol = 0
         matr%nrow = 0
      End If
   End Subroutine rfreediag

!!$ convert sparce matrix to full matrix
   Function zcsr_sptofull(amat) Result(b)
      Implicit None
      Type(zcsrmat) :: amat
      Complex(Kind=DEF_DBL_PREC) :: b(amat%nrow, amat%ncol)
      Integer :: i, j

      b = DEF_cmplx_zero

      Do i = 1, amat%nrow
         Do j = amat%ir(i), amat%ir(i + 1) - 1
            b(i, amat%jc(j)) = amat%a(j)
         End Do
      End Do

   End Function zcsr_sptofull

!!$************* INPLACE (!!!) CSR matrix  SHIFT DIAGONAL matrix elements routines *************

   Subroutine zrcsrshiftdiag(amat, b)
      Implicit None
      Type(zcsrmat) :: amat
      Real(Kind=DEF_DBL_PREC) :: b
      Integer :: count, idiag(amat%nrow + 1), nen, i1
      Logical :: test
      Integer :: k1old, k2old, j, k, ii, k1, k2, kold

      Call csr_diapos(amat%nrow, amat%ir, amat%jc, idiag, count)
      nen = amat%nrow - count
      If (nen > 0) Then
         If (amat%nnzmax < (amat%nnz + nen)) Call mrealloc(amat, amat%nnz + nen)

!!$         ko = ir(amat%nrow+1)+nen
         Do ii = amat%nrow, 1, -1
!!$     go through  row ii
            k1old = amat%ir(ii)
            k2old = amat%ir(ii + 1) - 1
            k2 = k2old + nen
            amat%ir(ii + 1) = k2 + 1
            test = (idiag(ii) == 0)

            If (test) Then
               Do kold = k2old, k1old, -1
                  k = kold + nen
                  j = amat%jc(kold)
                  If (test .And. (j < ii)) Then
                     test = .False.
                     nen = nen - 1
                     amat%a(k) = 0.0d0
                     amat%jc(k) = ii
                     idiag(ii) = k
                     k = k - 1
                  End If
                  amat%a(k) = amat%a(kold)
                  amat%jc(k) = j
               End Do
!!$     diagonal element has not been added yet.
               If (test) Then
                  nen = nen - 1
                  k = k - 1
                  amat%a(k) = 0.0d0
                  amat%jc(k) = ii
                  idiag(ii) = k
               End If
            Else
               k1 = k1old + nen
               amat%ir(ii) = k1
               amat%a(k1:k2) = amat%a(k1old:k2old)
               amat%jc(k1:k2) = amat%jc(k1old:k2old)
               idiag(ii) = idiag(ii) + nen
            End If
         End Do
         Write (*, *) 'I still not sure how this "csraddiag" works ! :)'
      Else
!$omp parallel do
         Do i1 = 1, amat%nrow
            amat%a(idiag(i1)) = amat%a(idiag(i1)) + b
         End Do
!$omp end parallel do
      End If

   End Subroutine zrcsrshiftdiag

   Subroutine zzcsrshiftdiag(amat, b)
      Implicit None
      Type(zcsrmat) :: amat
      Complex(Kind=DEF_DBL_PREC) :: b
!!$          Locals
      Integer :: count, idiag(amat%nrow + 1), nen, i1
!!$          Logical :: test
      Integer :: j, k
      Type(zcsrmat) :: tmat, tmat2

      Call csr_diapos(amat%nrow, amat%ir, amat%jc, idiag, count)
      nen = amat%nrow - count
      If (nen > 0) Then
         k = Min(amat%nrow, amat%ncol)
         Call alloc(tmat, amat%nrow, amat%ncol, k)
         Do j = 1, k, 1
            tmat%jc(j) = j
            tmat%ir(j) = j
         End Do
         tmat%ir((k + 1):(tmat%nrow + 1)) = k + 1
         tmat%nnz = k

         tmat%a(:) = b
         tmat2 = spmatadd(amat, tmat)

         Call free(tmat)
         Call free(amat)
         Call spcopy(tmat2, amat)
      Else
         Do i1 = 1, amat%nrow
            amat%a(idiag(i1)) = amat%a(idiag(i1)) + b
         End Do
      End If

   End Subroutine zzcsrshiftdiag

!!$************* INPLACE (!!!) CSR matrix  add/subs DIAGONAL matrix routines *************

   Subroutine csradddiag_zz(amat, bmat, sign)
      Implicit None
      Type(zcsrmat) :: amat
      Type(zdiagmat) :: bmat
      Integer :: sign
      Integer :: count, idiag(amat%nrow + 1), nen, i1
!!$          Logical :: test
      Integer :: j, k
      Type(zcsrmat) :: tmat2, tmat

      Call csr_diapos(amat%nrow, amat%ir, amat%jc, idiag, count)
      nen = amat%nrow - count
      If (nen > 0) Then
         k = Min(amat%nrow, amat%ncol)
         Call alloc(tmat, amat%nrow, amat%ncol, k)
         Do j = 1, k, 1
            tmat%jc(j) = j
            tmat%ir(j) = j
         End Do
         tmat%ir((k + 1):(tmat%nrow + 1)) = k + 1
         tmat%nnz = k

         tmat%a(1:k) = bmat%el(1:k)
         tmat2 = spmatadd(amat, tmat)
!!$             Call free (tmat)
!!$             Call free (amat)
!!$             Call spcopy (tmat2, amat)
         Deallocate (amat%a, amat%ir, amat%jc)
         amat%a = tmat2%a
         amat%ir = tmat2%ir
         amat%jc = tmat2%jc
         amat%nnz = tmat2%nnz
         amat%nnzmax = tmat2%nnzmax
         tmat2%alloc = 0

      Else
         If (sign > 0) Then
            Do i1 = 1, amat%nrow
               amat%a(idiag(i1)) = amat%a(idiag(i1)) + bmat%el(i1)
            End Do
         Else
            Do i1 = 1, amat%nrow
               amat%a(idiag(i1)) = amat%a(idiag(i1)) - bmat%el(i1)
            End Do
         End If
      End If

   End Subroutine csradddiag_zz

   Subroutine csradddiag_rr(amat, bmat, sign)
      Implicit None
      Type(rcsrmat) :: amat
      Type(rdiagmat) :: bmat
      Integer :: sign
      Integer, Pointer, Dimension(:) :: jc, ir
      Integer :: count, idiag(amat%nrow + 1), nen, i1
      Logical :: test
      Integer :: k1old, k2old, j, k, ii, k1, k2, kold

      Call csr_diapos(amat%nrow, amat%ir, amat%jc, idiag, count)
      nen = amat%nrow - count
      If (nen > 0) Then
         If (amat%nnzmax < (amat%nnz + nen)) Call mrealloc(amat, amat%nnz + nen)

!!$         ko = ir(amat%nrow+1)+nen
         Do ii = amat%nrow, 1, -1
!!$     go through  row ii
            k1old = amat%ir(ii)
            k2old = amat%ir(ii + 1) - 1
            k2 = k2old + nen
            amat%ir(ii + 1) = k2 + 1
            test = (idiag(ii) == 0)

            If (test) Then
               Do kold = k2old, k1old, -1
                  k = kold + nen
                  j = amat%jc(kold)
                  If (test .And. (j < ii)) Then
                     test = .False.
                     nen = nen - 1
                     amat%a(k) = 0.0d0
                     amat%jc(k) = ii
                     idiag(ii) = k
                     k = k - 1
                  End If
                  amat%a(k) = amat%a(kold)
                  amat%jc(k) = j
               End Do
!!$     diagonal element has not been added yet.
               If (test) Then
                  nen = nen - 1
                  k = k - 1
                  amat%a(k) = 0.0d0
                  amat%jc(k) = ii
                  idiag(ii) = k
               End If
            Else
               k1 = k1old + nen
               amat%ir(ii) = k1
               amat%a(k1:k2) = amat%a(k1old:k2old)
               amat%jc(k1:k2) = amat%jc(k1old:k2old)
               idiag(ii) = idiag(ii) + nen
            End If
         End Do
         Write (*, *) 'I still not sure how this "csraddiag" works ! :)'
      Else
         If (sign > 0) Then
            Do i1 = 1, amat%nrow
               amat%a(idiag(i1)) = amat%a(idiag(i1)) + bmat%el(i1)
            End Do
         Else
            Do i1 = 1, amat%nrow
               amat%a(idiag(i1)) = amat%a(idiag(i1)) - bmat%el(i1)
            End Do
         End If
      End If

   End Subroutine csradddiag_rr

   Subroutine csradddiag_zr(amat, bmat, sign)
      Implicit None
      Type(zcsrmat) :: amat
      Type(rdiagmat) :: bmat
      Integer :: sign
      Integer, Pointer, Dimension(:) :: jc, ir
      Integer :: count, idiag(amat%nrow + 1), nen, i1
      Logical :: test
      Integer :: k1old, k2old, j, k, ii, k1, k2, kold

      Call csr_diapos(amat%nrow, amat%ir, amat%jc, idiag, count)
      nen = amat%nrow - count
      If (nen > 0) Then
         If (amat%nnzmax < (amat%nnz + nen)) Call mrealloc(amat, amat%nnz + nen)

!!$         ko = ir(amat%nrow+1)+nen
         Do ii = amat%nrow, 1, -1
!!$     go through  row ii
            k1old = amat%ir(ii)
            k2old = amat%ir(ii + 1) - 1
            k2 = k2old + nen
            amat%ir(ii + 1) = k2 + 1
            test = (idiag(ii) == 0)

            If (test) Then
               Do kold = k2old, k1old, -1
                  k = kold + nen
                  j = amat%jc(kold)
                  If (test .And. (j < ii)) Then
                     test = .False.
                     nen = nen - 1
                     amat%a(k) = 0.0d0
                     amat%jc(k) = ii
                     idiag(ii) = k
                     k = k - 1
                  End If
                  amat%a(k) = amat%a(kold)
                  amat%jc(k) = j
               End Do
!!$     diagonal element has not been added yet.
               If (test) Then
                  nen = nen - 1
                  k = k - 1
                  amat%a(k) = 0.0d0
                  amat%jc(k) = ii
                  idiag(ii) = k
               End If
            Else
               k1 = k1old + nen
               amat%ir(ii) = k1
               amat%a(k1:k2) = amat%a(k1old:k2old)
               amat%jc(k1:k2) = amat%jc(k1old:k2old)
               idiag(ii) = idiag(ii) + nen
            End If
         End Do
         Write (*, *) 'I still not sure how this "csraddiag" works ! :)'
      Else
         If (sign > 0) Then
            Do i1 = 1, amat%nrow
               amat%a(idiag(i1)) = amat%a(idiag(i1)) + bmat%el(i1)
            End Do
         Else
            Do i1 = 1, amat%nrow
               amat%a(idiag(i1)) = amat%a(idiag(i1)) - bmat%el(i1)
            End Do
         End If
      End If

   End Subroutine csradddiag_zr

!!$**** get position of diagonal elements in CSR matrix ******
   Subroutine csr_diapos(nrow, ir, jc, idiag, count)
      Implicit None
      Integer :: k, i, nrow, count
      Integer :: ir(nrow + 1), jc(:), idiag(nrow)

      idiag = 0
      count = 0
      Do i = 1, nrow
         Do k = ir(i), ir(i + 1) - 1
            If (jc(k) == i) Then
               idiag(i) = k
               count = count + 1
            End If

         End Do
      End Do
   End Subroutine csr_diapos

!!$************* INPLACE (!!!) CSR matrix - DIAGONAL matrix multimlication routines *************
   Subroutine mmcsrdiag_zz(amat, bmat)
      Implicit None
      Type(zcsrmat) :: amat
      Type(zdiagmat) :: bmat
      Integer, Pointer, Dimension(:) :: jc

      Integer :: i1, nnz

      nnz = amat%nnz

      Do i1 = 1, nnz
         amat%a(i1) = amat%a(i1)*bmat%el(amat%jc(i1))
      End Do

   End Subroutine mmcsrdiag_zz

   Subroutine mmcsrdiag_zr(amat, bmat)
      Implicit None
      Type(zcsrmat) :: amat
      Type(rdiagmat) :: bmat
      Integer, Pointer, Dimension(:) :: jc

      Integer :: i1, nnz

      nnz = amat%nnz

      Do i1 = 1, nnz
         amat%a(i1) = amat%a(i1)*bmat%el(amat%jc(i1))
      End Do

   End Subroutine mmcsrdiag_zr

   Subroutine mmcsrdiag_rr(amat, bmat)
      Implicit None
      Type(rcsrmat) :: amat
      Type(rdiagmat) :: bmat
      Integer, Pointer, Dimension(:) :: jc

      Integer :: i1, nnz

      nnz = amat%nnz

      Do i1 = 1, nnz
         amat%a(i1) = amat%a(i1)*bmat%el(amat%jc(i1))
      End Do

   End Subroutine mmcsrdiag_rr

!!$************* INPLACE (!!!) DIAGONAL matrix - CSR matrix  multimplication routines *************

   Subroutine mmdiagcsr_rr(amat, bmat)
      Implicit None
      Type(rcsrmat) :: bmat
      Type(rdiagmat) :: amat
      Integer, Pointer, Dimension(:) :: ir

      Integer :: i1, nr

      nr = bmat%nrow

      Do i1 = 1, nr
         bmat%a(bmat%ir(i1):(bmat%ir(i1 + 1) - 1)) = bmat%a(bmat%ir(i1):(bmat%ir(i1 + 1) - 1))*amat%el(i1)
      End Do

   End Subroutine mmdiagcsr_rr

   Subroutine mmdiagcsr_zz(amat, bmat)
      Implicit None
      Type(zcsrmat) :: bmat
      Type(zdiagmat) :: amat
      Integer, Pointer, Dimension(:) :: ir

      Integer :: i1, nr

      nr = bmat%nrow

      Do i1 = 1, nr
         bmat%a(bmat%ir(i1):(bmat%ir(i1 + 1) - 1)) = bmat%a(bmat%ir(i1):(bmat%ir(i1 + 1) - 1))*amat%el(i1)
      End Do

   End Subroutine mmdiagcsr_zz

   Subroutine mmdiagcsr_rz(amat, bmat)
      Implicit None
      Type(zcsrmat) :: bmat
      Type(rdiagmat) :: amat
      Integer, Pointer, Dimension(:) :: ir

      Integer :: i1, nr

      nr = bmat%nrow

      Do i1 = 1, nr
         bmat%a(bmat%ir(i1):(bmat%ir(i1 + 1) - 1)) = bmat%a(bmat%ir(i1):(bmat%ir(i1 + 1) - 1))*amat%el(i1)
      End Do

   End Subroutine mmdiagcsr_rz

!!$************* NON-INPLACE (!!!) DIAGONAL matrix - CSR matrix  multimplication functions, return new CSR matrix *********
   Function fmmdiagcsr_rz(amat, bmat)
      Implicit None
      Type(rdiagmat) :: amat
      Type(zcsrmat) :: bmat, fmmdiagcsr_rz
      Call alloc(fmmdiagcsr_rz, bmat%nnz, bmat%nrow, bmat%ncol)
      Call spcopy(bmat, fmmdiagcsr_rz)
      Call inspmatmul(amat, fmmdiagcsr_rz)
      Return
   End Function fmmdiagcsr_rz

   Function fmmdiagcsr_rr(amat, bmat)
      Implicit None
      Type(rdiagmat) :: amat
      Type(rcsrmat) :: bmat, fmmdiagcsr_rr
      Call alloc(fmmdiagcsr_rr, bmat%nnz, bmat%nrow, bmat%ncol)
      Call spcopy(bmat, fmmdiagcsr_rr)
      Call inspmatmul(amat, fmmdiagcsr_rr)
      Return
   End Function fmmdiagcsr_rr

   Function fmmdiagcsr_zr(amat, bmat)
      Implicit None
      Type(zdiagmat) :: amat
      Type(rcsrmat) :: bmat
      Type(zcsrmat) :: fmmdiagcsr_zr
      Call alloc(fmmdiagcsr_zr, bmat%nnz, bmat%nrow, bmat%ncol)
      Call spcopy(bmat, fmmdiagcsr_zr)
      Call inspmatmul(amat, fmmdiagcsr_zr)
      Return
   End Function fmmdiagcsr_zr

   Function fmmdiagcsr_zz(amat, bmat)
      Implicit None
      Type(zdiagmat) :: amat
      Type(zcsrmat) :: bmat, fmmdiagcsr_zz
      Call alloc(fmmdiagcsr_zz, bmat%nnz, bmat%nrow, bmat%ncol)
      Call spcopy(bmat, fmmdiagcsr_zz)
      Call inspmatmul(amat, fmmdiagcsr_zz)
      Return
   End Function fmmdiagcsr_zz

!!$************* NON-INPLACE (!!!) CSR matrix - DIAGONAL matrix  multiplication functions, return new CSR matrix *****

   Function fmmcsrdiag_zz(amat, bmat)
      Implicit None
      Type(zcsrmat) :: amat
      Type(zdiagmat) :: bmat
      Type(zcsrmat) :: fmmcsrdiag_zz

      Call alloc(fmmcsrdiag_zz, amat%nnz, amat%nrow, amat%ncol)
      Call spcopy(amat, fmmcsrdiag_zz)
      Call inspmatmul(fmmcsrdiag_zz, bmat)
      Return
   End Function fmmcsrdiag_zz

   Function fmmcsrdiag_zr(amat, bmat)
      Implicit None
      Type(zcsrmat) :: amat
      Type(rdiagmat) :: bmat
      Type(zcsrmat) :: fmmcsrdiag_zr

      Call alloc(fmmcsrdiag_zr, amat%nnz, amat%nrow, amat%ncol)
      Call spcopy(amat, fmmcsrdiag_zr)
      Call inspmatmul(fmmcsrdiag_zr, bmat)
      Return
   End Function fmmcsrdiag_zr

   Function fmmcsrdiag_rr(amat, bmat)
      Implicit None
      Type(rcsrmat) :: amat
      Type(rdiagmat) :: bmat
      Type(rcsrmat) :: fmmcsrdiag_rr

      Call alloc(fmmcsrdiag_rr, amat%nnz, amat%nrow, amat%ncol)
      Call spcopy(amat, fmmcsrdiag_rr)
      Call inspmatmul(fmmcsrdiag_rr, bmat)
      Return
   End Function fmmcsrdiag_rr

   Function fmmcsrdiag_rz(amat, bmat)
      Implicit None
      Type(rcsrmat) :: amat
      Type(zdiagmat) :: bmat
      Type(zcsrmat) :: fmmcsrdiag_rz

      Call alloc(fmmcsrdiag_rz, amat%nnz, amat%nrow, amat%ncol)
      Call spcopy(amat, fmmcsrdiag_rz)
      Call inspmatmul(fmmcsrdiag_rz, bmat)
      Return
   End Function fmmcsrdiag_rz

!!$************* CSR matrix - DENSE VECTOR multiplication functions, return DENSE VECTOR *************
   Function fmmcsrvector_zz(amat, bvec) Result(rvec)
      Implicit None
      Type(zcsrmat) :: amat
      Complex(Kind=DEF_DBL_PREC) :: bvec(:)
      Complex(Kind=DEF_DBL_PREC) :: rvec(amat%nrow)
      Integer :: i1, nr, j1
      Integer, Pointer, Dimension(:) :: jc

      rvec = 0.0d0
      nr = amat%nrow
      Do i1 = 1, nr
         Do j1 = amat%ir(i1), amat%ir(i1 + 1) - 1
            rvec(i1) = rvec(i1) + bvec(jc(j1))*amat%a(j1)
         End Do
      End Do
      Return
   End Function fmmcsrvector_zz

   Subroutine mmcsrvector_zz_s(amat, bvec, rvec)
      Implicit None
      Type(zcsrmat) :: amat
      Complex(Kind=DEF_DBL_PREC) :: bvec(:)
      Complex(Kind=DEF_DBL_PREC) :: rvec(amat%nrow)
      Integer :: i1, nr, j1
      Integer, Pointer, Dimension(:) :: jc

      rvec = 0.0d0
      nr = amat%nrow
      Do i1 = 1, nr
         Do j1 = amat%ir(i1), amat%ir(i1 + 1) - 1
            rvec(i1) = rvec(i1) + bvec(jc(j1))*amat%a(j1)
         End Do
      End Do
      Return
   End Subroutine mmcsrvector_zz_s

!!$************* CSR matrix - DENSE matrix multiplication functions, return DENSE matrix *************
   Subroutine mmcsrdense_zz(amat, bmat, cmat)
      Implicit None
      Type(zcsrmat), Intent(In) :: amat
      Type(zdensemat), Intent(In) :: bmat
      Type(zdensemat), Intent(Inout) :: cmat

      Call alloc(cmat, amat%nrow, bmat%ncol)
      Call mmcsrdense_zz_base(amat%jc, amat%ir, amat%a, bmat%bl, cmat%bl, amat%nrow, bmat%ncol)

   End Subroutine mmcsrdense_zz

   Subroutine mmcsrdense_zz_s(amat, bmat, cmat)
      Implicit None
      Type(zcsrmat), Intent(In) :: amat
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: bmat(:, :)
      Complex(Kind=DEF_DBL_PREC) :: cmat(amat%nrow, size(bmat, 2))

      cmat = dcmplx(0.0d0, 0.0d0)
      Call mmcsrdense_zz_base(amat%jc, amat%ir, amat%a, bmat, cmat, amat%nrow, size(bmat, 2))
      Return
   End Subroutine mmcsrdense_zz_s

   Subroutine mmcsrdense_zz_base(jc, ir, a, b, c, nr, nc)
      Implicit None
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: b(:, :), a(:)
      Complex(Kind=DEF_DBL_PREC), Intent(Out) :: c(:, :)
      Integer, Intent(In) :: jc(:), ir(:), nr, nc
!!$ Local vars
      Integer :: i1, i2, j1
      Do i1 = 1, nr
         Do i2 = 1, nc
            Do j1 = ir(i1), ir(i1 + 1) - 1
               c(i1, i2) = c(i1, i2) + b(jc(j1), i2)*a(j1)
            End Do
         End Do
      End Do
   End Subroutine mmcsrdense_zz_base

!!$************* CSR matrix - DENSE matrix multiplication functions, return DENSE matrix *************
   Function fmmcsrdense_zz(amat, bmat)
      Implicit None
      Type(zcsrmat), Intent(In) :: amat
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: bmat(:, :)
      Complex(Kind=DEF_DBL_PREC) :: fmmcsrdense_zz(amat%nrow, size(bmat, 2))

      fmmcsrdense_zz = 0.0d0
      Call mmcsrdense_zz_base(amat%jc, amat%ir, amat%a, bmat, fmmcsrdense_zz, amat%nrow, size(bmat, 2))
      Return
   End Function fmmcsrdense_zz

   Function fmmcsrdense_zr(amat, bmat)
      Implicit None
      Type(zcsrmat) :: amat
      Real(Kind=DEF_DBL_PREC) :: bmat(:, :)
      Complex(Kind=DEF_DBL_PREC) :: fmmcsrdense_zr(amat%nrow, size(bmat, 2))
      Integer :: i1, i2, nc, nr, j1

      fmmcsrdense_zr = 0.0d0
      nr = amat%nrow
      nc = size(bmat, 2)
      Do i1 = 1, nr
         Do i2 = 1, nc
            Do j1 = amat%ir(i1), amat%ir(i1 + 1) - 1
               fmmcsrdense_zr(i1, i2) = fmmcsrdense_zr(i1, i2) + bmat(amat%jc(j1), i2)*amat%a(j1)
            End Do
         End Do
      End Do
      Return
   End Function fmmcsrdense_zr

   Function fmmcsrdense_rz(amat, bmat)
      Implicit None
      Type(rcsrmat) :: amat
      Complex(Kind=DEF_DBL_PREC) :: bmat(:, :)
      Complex(Kind=DEF_DBL_PREC) :: fmmcsrdense_rz(amat%nrow, size(bmat, 2))
      Integer :: i1, i2, nc, nr, j1

      fmmcsrdense_rz = 0.0d0
      nr = amat%nrow
      nc = size(bmat, 2)
      Do i1 = 1, nr
         Do i2 = 1, nc
            Do j1 = amat%ir(i1), amat%ir(i1 + 1) - 1
               fmmcsrdense_rz(i1, i2) = fmmcsrdense_rz(i1, i2) + bmat(amat%jc(j1), i2)*amat%a(j1)
            End Do
         End Do
      End Do
      Return
   End Function fmmcsrdense_rz

   Function fmmcsrdense_rr(amat, bmat)
      Implicit None
      Type(rcsrmat) :: amat
      Real(Kind=DEF_DBL_PREC) :: bmat(:, :)
      Real(Kind=DEF_DBL_PREC) :: fmmcsrdense_rr(amat%nrow, size(bmat, 2))
      Integer :: i1, i2, nc, nr, j1

      fmmcsrdense_rr = 0.0d0
      nr = amat%nrow
      nc = size(bmat, 2)
      Do i1 = 1, nr
         Do i2 = 1, nc
            Do j1 = amat%ir(i1), amat%ir(i1 + 1) - 1
               fmmcsrdense_rr(i1, i2) = fmmcsrdense_rr(i1, i2) + bmat(amat%jc(j1), i2)*amat%a(j1)
            End Do
         End Do
      End Do
      Return
   End Function fmmcsrdense_rr

   Function zgetcsrsubmat(amat, i1, i2, j1, j2)
      Implicit None
      Integer i1, i2, j1, j2
      Type(zcsrmat) :: amat, zgetcsrsubmat
      Integer :: nr, nc, in1, in2, nnz, sc, cj, sr, mnnz

      nr = i2 - i1 + 1
      nc = j2 - j1 + 1
      sc = j1 - 1
      sr = i1 - 2
      mnnz = amat%ir(j2 + 1) - amat%ir(j1)
      Call alloc(zgetcsrsubmat, mnnz, nr, nc)

      zgetcsrsubmat%ir(1) = 1
      nnz = 0

      Do in1 = i1, i2
         Do in2 = amat%ir(in1), amat%ir(in1 + 1) - 1
            cj = amat%jc(in2) - sc
            If (cj > 0) Then
               If (cj > nc) Then
                  Exit
               End If
               nnz = nnz + 1
               zgetcsrsubmat%jc(nnz) = cj
               zgetcsrsubmat%a(nnz) = amat%a(in2)
            End If
         End Do
         zgetcsrsubmat%ir(in1 - sr) = nnz + 1
      End Do
      zgetcsrsubmat%nnz = nnz
      Call mrealloc(zgetcsrsubmat, nnz)
   End Function zgetcsrsubmat

   Function rgetcsrsubmat(amat, i1, i2, j1, j2)
      Implicit None
      Integer i1, i2, j1, j2
      Type(rcsrmat) :: amat, rgetcsrsubmat
      Integer :: nr, nc, in1, in2, nnz, sc, cj, sr, mnnz

      nr = i2 - i1 + 1
      nc = j2 - j1 + 1
      sc = j1 - 1
      sr = i1 - 2
      mnnz = amat%ir(j2 + 1) - amat%ir(j1)
      Call alloc(rgetcsrsubmat, mnnz, nr, nc)

      rgetcsrsubmat%ir(1) = 1
      nnz = 0

      Do in1 = i1, i2
         Do in2 = amat%ir(in1), amat%ir(in1 + 1) - 1
            cj = amat%jc(in2) - sc
            If (cj > 0) Then
               If (cj > nc) Then
                  Exit
               End If
               nnz = nnz + 1
               rgetcsrsubmat%jc(nnz) = cj
               rgetcsrsubmat%a(nnz) = amat%a(in2)
            End If
         End Do
         rgetcsrsubmat%ir(in1 - sr) = nnz + 1
      End Do
      rgetcsrsubmat%nnz = nnz
      Call mrealloc(rgetcsrsubmat, nnz)
   End Function rgetcsrsubmat

   Function zgetrowind(amat)
      Implicit None
      Type(zcsrmat) :: amat
      Integer, pointer :: zgetrowind(:)
      Integer :: i, nrow, k1, k2, k
      allocate (zgetrowind(amat%nnz))
      nrow = amat%nrow
      Do i = nrow, 1, -1
         k1 = amat%ir(i + 1) - 1
         k2 = amat%ir(i)
         Do k = k1, k2, -1
            zgetrowind(k) = i
         End Do
      End Do
      Return
   End Function zgetrowind

   Function rgetrowind(amat)
      Implicit None
      Type(rcsrmat) :: amat
      Integer, pointer :: rgetrowind(:)
      Integer :: i, nrow, k1, k2, k
      allocate (rgetrowind(amat%nnz))
      nrow = amat%nrow
      Do i = nrow, 1, -1
         k1 = amat%ir(i + 1) - 1
         k2 = amat%ir(i)
         Do k = k1, k2, -1
            rgetrowind(k) = i
         End Do
      End Do
      Return
   End Function rgetrowind

   Subroutine addblocks(ham, blocks, num)
!!$ Add blocks to sparse matrix.
!!$ I assume that in EVERY row, which will be added
!!$ There is at least one nonzero element in original matrix
      Implicit None
      Type(zcsrmat) :: ham
      Type(matblk) :: blocks(num)
      Integer :: num
!!$ Local
      Integer, Allocatable :: newjc(:), newir(:)
      Integer, Pointer :: pos(:), inds(:)
      Complex(Kind=DEF_DBL_PREC), Allocatable ::  newa(:), bptr(:, :)
      Integer :: ib, i, j, nrow, ncol, bnr, bnc, nnz, newnnz
      Integer :: frow, fcol, lrow, offs, nel, ind, ind1
!!$ Code
      Type(tspblpos) :: bnds(num)

      nrow = ham%nrow
      ncol = ham%ncol
      offs = 0
      nnz = ham%nnz
!!$ Calculate amount of space we need to hold result
      nel = 0
      Do ib = 1, num
         bnr = size(blocks(ib)%bl, 1)
         bnc = size(blocks(ib)%bl, 2)
         nel = Max(nel, bnc)
         Allocate (bnds(ib)%pos(2, bnr))
         bnds(ib)%pos = 0
         offs = offs + bnc*bnr
         ind = blocks(ib)%fcol
         ind1 = ind + bnc
         Do j = blocks(ib)%frow, blocks(ib)%frow + bnr - 1
            pos => bnds(ib)%pos(:, j - blocks(ib)%frow + 1)
            Do i = ham%ir(j), ham%ir(j + 1) - 1
               If (ham%jc(i) >= ind .And. ham%jc(i) < ind1) Then
                  If (pos(1) == 0) Then
                     pos(1) = i
                  Else
                     pos(2) = i
                  End If
                  offs = offs - 1
               End If
            End Do
            If (pos(2) == 0) pos(2) = pos(1)
         End Do
      End Do
!!$ Allocate matrix
      newnnz = ham%nnz + offs
      Allocate (newa(newnnz))
      Allocate (newjc(newnnz))
      Allocate (newir(nrow + 1))
      Allocate (inds(nel))
      Do i = 1, nel
         inds(i) = i - 1
      End Do

      lrow = 1
      newir(1) = 1
      nnz = 0
      Do ib = 1, num + 1

         If (ib > num) Then
            frow = nrow + 1
         Else
            frow = blocks(ib)%frow
         End If
!!$ Add elements before firt row of the block or before end of matrix
         offs = newir(lrow) - ham%ir(lrow)
         newir(lrow + 1:frow) = offs + ham%ir(lrow + 1:frow)
         nel = ham%ir(frow) - ham%ir(lrow)
         ind = cpdata(nel, ham%ir(lrow), ham%jc, ham%a, newir(lrow), newjc, newa)
         nnz = nnz + nel

         If (ib > num) Exit

         bptr = blocks(ib)%bl
         bnr = size(bptr, 1)
         bnc = size(bptr, 2)
         fcol = blocks(ib)%fcol
         Do j = frow, frow + bnr - 1
            pos => bnds(ib)%pos(:, j - frow + 1)
!!$ Add elements before block
            nel = pos(1) - ham%ir(j)
            ind = cpdata(nel, ham%ir(j), ham%jc, ham%a, newir(j), newjc, newa)
            nnz = nnz + nel
!!$ Add block
            ind1 = ind + bnc - 1
            newa(ind:ind1) = bptr(j - frow + 1, :)
            newjc(ind:ind1) = fcol + inds(1:bnc)
            nnz = nnz + bnc
!!$ Add old elements to block
            ind = ind - fcol
            Do i = pos(1), pos(2), 1
               newa(ind + ham%jc(i)) = newa(ind + ham%jc(i)) + ham%a(i)
            End Do
!!$ Add elements after block
            nel = ham%ir(j + 1) - pos(2) - 1
            ind1 = ind1 + 1
            ind = cpdata(nel, pos(2) + 1, ham%jc, ham%a, ind1, newjc, newa)
            nnz = nnz + nel
            newir(j + 1) = ind
         End Do
         lrow = frow + bnr
         Deallocate (bnds(ib)%pos)
      End Do
      If (newnnz /= nnz) Then
         Write (*, *) "!!!", nnz, newnnz
         Stop
      End If
      ham%nnz = nnz
      ham%nnzmax = nnz
      Deallocate (ham%a)
      Deallocate (ham%ir)
      Deallocate (ham%jc)
      Deallocate (inds)

      Allocate (ham%a(size(newa)))
      Allocate (ham%ir(size(newir)))
      Allocate (ham%jc(size(newjc)))
      ham%jc = newjc
      ham%ir = newir
      ham%a = newa
   Contains
      Function cpdata(nel, i1, jc1, a1, j1, jc2, a2) Result(j2)
         Implicit None
         Integer :: jc1(:), jc2(:)
         Complex(Kind=DEF_DBL_PREC) :: a1(:), a2(:)
         Integer :: nel, j1, j2, i1
!!$ Local
         Integer :: i2
!!$ Code
         If (nel > 0) Then
            i2 = i1 + nel - 1
            j2 = j1 + nel - 1
            jc2(j1:j2) = jc1(i1:i2)
            a2(j1:j2) = a1(i1:i2)
            j2 = j2 + 1
         Else
            j2 = j1
         End If
      End Function cpdata
   End Subroutine addblocks

   Subroutine dumpzcoo(mat, fname)
      Implicit None
      Type(zcoomat) :: mat
      Character(Len=*) :: fname
      Integer :: s1

      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, mat%nnz
         Write (dumping_unit, *) mat%irow(s1), mat%icol(s1), real(mat%a(s1)), imag(mat%a(s1))
      End Do
      Close (dumping_unit)
   End Subroutine dumpzcoo

   Subroutine dumpzcsr(mat, fname)
      Implicit None
      Type(zcsrmat) :: mat
      Character(Len=*) :: fname
      Integer ::  s1
      Integer, pointer ::  irow(:)
      irow => getrowind(mat)

      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, mat%nnz
         Write (dumping_unit, *) irow(s1), mat%jc(s1), real(mat%a(s1)), imag(mat%a(s1))
      End Do
      Write (dumping_unit, *) mat%nrow, mat%ncol, 0.0d0, 0.0d0
      Close (dumping_unit)
      deallocate (irow)
   End Subroutine dumpzcsr

   Subroutine dumprcsr(mat, fname)
      Implicit None
      Type(rcsrmat) :: mat
      Character(Len=*) :: fname
      Integer ::  s1
      Integer, pointer ::  irow(:)
      irow => getrowind(mat)
      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, mat%nnz
         Write (dumping_unit, *) irow(s1), mat%jc(s1), mat%a(s1)
      End Do
      Close (dumping_unit)
      deallocate (irow)
   End Subroutine dumprcsr

   Subroutine dumpzdiag(mat, fname)
      Implicit None
      Type(zdiagmat) :: mat
      Character(Len=*) :: fname
      Integer :: s1
      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, mat%ncol
         Write (dumping_unit, *) s1, s1, real(mat%el(s1)), imag(mat%el(s1))
      End Do
      Close (dumping_unit)
   End Subroutine dumpzdiag

   Subroutine dumprdiag(mat, fname)
      Implicit None
      Type(rdiagmat) :: mat
      Character(Len=*) :: fname
      Integer :: s1
      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, mat%ncol
         Write (dumping_unit, *) s1, s1, mat%el(s1)
      End Do
      Close (dumping_unit)
   End Subroutine dumprdiag

   Subroutine dumpzdense(mat, fname)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: mat(:, :)
      Character(Len=*) :: fname
      Integer :: s1, s2

      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s2 = 1, size(mat, 2)
         Do s1 = 1, size(mat, 1)
            Write (dumping_unit, *) s1, s2, real(mat(s1, s2)), imag(mat(s1, s2))
         End Do
      End Do
      Close (dumping_unit)
   End Subroutine dumpzdense

   Subroutine dumpzdensemat(mat, fname)
      Implicit None
      Type(zdensemat) :: mat
      Character(Len=*) :: fname
      Integer :: s1, s2

      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s2 = 1, mat%ncol
         Do s1 = 1, mat%nrow
            Write (dumping_unit, *) s1, s2, real(mat%bl(s1, s2)), imag(mat%bl(s1, s2))
         End Do
      End Do
      Close (dumping_unit)
   End Subroutine dumpzdensemat

   Subroutine dumprdense(mat, fname)
      Implicit None
      Real(Kind=DEF_DBL_PREC) :: mat(:, :)
      Character(Len=*) :: fname
      Integer :: s1, s2
      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s2 = 1, size(mat, 2)
         Do s1 = 1, size(mat, 1)
            Write (dumping_unit, *) s1, s2, mat(s1, s2)
         End Do
      End Do
      Close (dumping_unit)
   End Subroutine dumprdense

   Subroutine dumprvector(mat, fname)
      Implicit None
      Real(Kind=DEF_DBL_PREC) :: mat(:)
      Character(Len=*) :: fname
      Integer :: s1
      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, size(mat)
         Write (dumping_unit, *) s1, 1, mat(s1)
      End Do
      Close (dumping_unit)
   End Subroutine dumprvector

   Subroutine dumpzvector(mat, fname)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: mat(:)
      Character(Len=*) :: fname
      Integer :: s1
      Open (Unit=dumping_unit, File=fname, Action='write')
      Do s1 = 1, size(mat)
         Write (dumping_unit, *) s1, 1, real(mat(s1)), imag(mat(s1))
      End Do
      Close (dumping_unit)
   End Subroutine dumpzvector

   Subroutine z_coo2csr(amat, bmat)
      Implicit None
      Type(zcsrmat) :: bmat
      Type(zcoomat) :: amat

      If (bmat%nnzmax < amat%nnz) Then
         Call mrealloc(bmat, amat%nnz)
      End If

      Call zcoocsr(amat%nrow, amat%nnz, amat%a, amat%irow, amat%icol, bmat%a, bmat%jc, bmat%ir)
      bmat%nnz = amat%nnz
      bmat%nrow = amat%nrow
      bmat%ncol = amat%ncol
   End Subroutine z_coo2csr

   Subroutine r_coo2csr(amat, bmat)
      Implicit None
      Type(rcsrmat) :: bmat
      Type(rcoomat) :: amat

      If ((bmat%nnzmax < amat%nnz) .Or. bmat%nrow /= amat%nrow) Then
         Call mrealloc(bmat, amat%nnz, amat%nrow)
      End If

      Call rcoocsr(amat%nrow, amat%nnz, amat%a, amat%irow, amat%icol, bmat%a, bmat%jc, bmat%ir)
      bmat%nnz = amat%nnz
      bmat%nrow = amat%nrow
      bmat%ncol = amat%ncol
   End Subroutine r_coo2csr

   Subroutine rcoocsr(nrow, nnz, a, ir, jc, ao, jao, iao)
      Implicit None
!!$-----------------------------------------------------------------------
      Real(Kind=DEF_DBL_PREC) :: a(*), ao(*), x
      Integer :: ir(*), jc(*), jao(*), iao(*)
!!$-----------------------------------------------------------------------
!!$  Coordinate     to   Compressed Sparse Row
!!$-----------------------------------------------------------------------
!!$ converts a matrix that is stored in coordinate format
!!$  a, ir, jc into a row general sparse ao, jao, iao format.
!!$
!!$ on entry:
!!$---------
!!$ nrow        = dimension of the matrix
!!$ nnz        = number of nonzero elements in matrix
!!$ a,
!!$ ir,
!!$ jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!!$         nonzero elements of the matrix with a(k) = actual real value of
!!$           the elements, ir(k) = its row number and jc(k) = its column
!!$          number. The order of the elements is arbitrary.
!!$
!!$ on return:
!!$-----------
!!$ ir         is destroyed
!!$
!!$ ao, jao, iao = matrix in general sparse matrix format with ao
!!$         continung the real values, jao containing the column indices,
!!$        and iao being the pointer to the beginning of the row,
!!$        in arrays ao, jao.
!!$
!!$ Notes:
!!$------ This routine is NOT in place.  See coicsr
!!$
!!$------------------------------------------------------------------------
      Integer :: k, j, i, nnz, nrow, k0, iad

      Do k = 1, nrow + 1
         iao(k) = 0
      End Do

!!$ determine row-lengths.
      Do k = 1, nnz
         iao(ir(k)) = iao(ir(k)) + 1
      End Do
!!$ starting position of each row..
      k = 1
      Do j = 1, nrow + 1
         k0 = iao(j)
         iao(j) = k
         k = k + k0
      End Do

!!$ go through the structure  once more. Fill in output matrix.
      Do k = 1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) = x
         jao(iad) = j
         iao(i) = iad + 1
      End Do

!!$ shift back iao
      Do j = nrow, 1, -1
         iao(j + 1) = iao(j)
      End Do
      iao(1) = 1
      Return
   End Subroutine rcoocsr

   Subroutine zcoocsr(nrow, nnz, a, ir, jc, ao, jao, iao)
      Implicit None
!!$-----------------------------------------------------------------------
      Complex(Kind=DEF_DBL_PREC) :: a(*), ao(*), x
      Integer :: ir(*), jc(*), jao(*), iao(*)
!!$-----------------------------------------------------------------------
!!$  Coordinate     to   Compressed Sparse Row
!!$-----------------------------------------------------------------------
!!$ converts a matrix that is stored in coordinate format
!!$  a, ir, jc into a row general sparse ao, jao, iao format.
!!$
!!$ on entry:
!!$---------
!!$ nrow        = dimension of the matrix
!!$ nnz        = number of nonzero elements in matrix
!!$ a,
!!$ ir,
!!$ jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!!$         nonzero elements of the matrix with a(k) = actual real value of
!!$           the elements, ir(k) = its row number and jc(k) = its column
!!$          number. The order of the elements is arbitrary.
!!$
!!$ on return:
!!$-----------
!!$ ir         is destroyed
!!$
!!$ ao, jao, iao = matrix in general sparse matrix format with ao
!!$         continung the real values, jao containing the column indices,
!!$        and iao being the pointer to the beginning of the row,
!!$        in arrays ao, jao.
!!$
!!$ Notes:
!!$------ This routine is NOT in place.  See coicsr
!!$
!!$------------------------------------------------------------------------
      Integer :: k, j, i, nnz, nrow, k0, iad

      Do k = 1, nrow + 1
         iao(k) = 0
      End Do

!!$ determine row-lengths.
      Do k = 1, nnz
         iao(ir(k)) = iao(ir(k)) + 1
      End Do
!!$ starting position of each row..
      k = 1
      Do j = 1, nrow + 1
         k0 = iao(j)
         iao(j) = k
         k = k + k0
      End Do

!!$ go through the structure  once more. Fill in output matrix.
      Do k = 1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) = x
         jao(iad) = j
         iao(i) = iad + 1
      End Do

!!$ shift back iao
      Do j = nrow, 1, -1
         iao(j + 1) = iao(j)
      End Do
      iao(1) = 1
      Return
   End Subroutine zcoocsr

   Function r_conv2csr(amat)
      Implicit None
      Type(rcsrmat) :: r_conv2csr
      Type(rcoomat) :: amat
!!$    integer,pointer::ir1(:)
      Integer :: iwk(amat%nrow + 1)
      Call rcoicsr(amat%nrow, amat%nnz, 1, amat%a, amat%icol, amat%irow, iwk)
      r_conv2csr%alloc = 1
      r_conv2csr%nnz = amat%nnz
      r_conv2csr%nrow = amat%nrow
      r_conv2csr%ncol = amat%ncol
      r_conv2csr%nnzmax = amat%nnzmax
      r_conv2csr%a = amat%a
      r_conv2csr%jc = amat%icol
      Allocate (r_conv2csr%ir(amat%nrow + 1))
      r_conv2csr%ir(1:amat%nrow + 1) = amat%irow(1:amat%nrow + 1)
      Deallocate (amat%irow)
      Deallocate (amat%icol)
      Deallocate (amat%a)
      amat%alloc = 0
      amat%nnzmax = 0
      amat%nnz = 0
      amat%ncol = 0
      amat%nrow = 0
   End Function r_conv2csr

   Function z_conv2csr(amat)
      Implicit None
      Type(zcsrmat) :: z_conv2csr
      Type(zcoomat) :: amat
!!$    integer,pointer::ir1(:)
      Integer :: iwk(amat%nrow + 1)
      Call zcoicsr(amat%nrow, amat%nnz, 1, amat%a, amat%icol, amat%irow, iwk)
      z_conv2csr%alloc = 1
      z_conv2csr%nnz = amat%nnz
      z_conv2csr%nrow = amat%nrow
      z_conv2csr%ncol = amat%ncol
      z_conv2csr%nnzmax = amat%nnzmax
      z_conv2csr%a = amat%a
      z_conv2csr%jc = amat%icol
      Allocate (z_conv2csr%ir(amat%nrow + 1))
      z_conv2csr%ir(1:amat%nrow + 1) = amat%irow(1:amat%nrow + 1)
      Deallocate (amat%irow)
      Deallocate (amat%icol)
      Deallocate (amat%a)
      amat%alloc = 0
      amat%nnzmax = 0
      amat%nnz = 0
      amat%ncol = 0
      amat%nrow = 0
   End Function z_conv2csr

   Subroutine rcoicsr(n, nnz, job, a, ja, ia, iwk)
      Implicit None
      Integer :: nnz, n, ia(nnz), ja(nnz), iwk(n + 1), job
      Real(Kind=DEF_DBL_PREC) a(*)
!!$------------------------------------------------------------------------
!!$ IN-PLACE coo-csr conversion routine.
!!$-----------------------------------------------------------------------
!!$ this subroutine converts a matrix stored in coordinate format into
!!$ the csr format. The conversion is done in place in that the arrays
!!$ a,ja,ia of the result are overwritten onto the original arrays.
!!$------------------------------------------------------------------------
!!$ on entry:
!!$---------
!!$ n        = integer. row dimension of A.
!!$ nnz        = integer. number of nonzero elements in A.
!!$ job   = integer. Job indicator. when job=1, the real values in a are
!!$         filled. Otherwise a is not touched and the structure of the
!!$         array only (i.e. ja, ia)  is obtained.
!!$ a        = real array of size nnz (number of nonzero elements in A)
!!$         containing the nonzero elements
!!$ ja        = integer array of length nnz containing the column positions
!!$           of the corresponding elements in a.
!!$ ia        = integer array of length nnz containing the row positions
!!$           of the corresponding elements in a.
!!$ iwk        = integer work array of length n+1
!!$ on return:
!!$----------
!!$ a
!!$ ja
!!$ ia        = contains the compressed sparse row data structure for the
!!$         resulting matrix.
!!$ Note:
!!$-------
!!$         the entries of the output matrix are not sorted (the column
!!$         indices in each are not in increasing order) use coocsr
!!$         if you want them sorted.
!!$----------------------------------------------------------------------c
!!$  Coded by Y. Saad, Sep. 26 1989                                      c
!!$----------------------------------------------------------------------c
      Real(Kind=DEF_DBL_PREC) :: t, tnext
      Logical :: values
      Integer :: k, i, init, j, ipos, inext, jnext
!!$-----------------------------------------------------------------------
      values = (job .Eq. 1)
!!$ find pointer array for resulting matrix.
      iwk = 0

      Do k = 1, nnz
         i = ia(k)
         iwk(i + 1) = iwk(i + 1) + 1
      End Do
!!$------------------------------------------------------------------------
      iwk(1) = 1
      Do i = 2, n
         iwk(i) = iwk(i - 1) + iwk(i)
      End Do

!!$     loop for a cycle in chasing process.

      init = 1
      k = 0
5     If (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
!!$------------------------------------------------------------------------
6     k = k + 1
!!$     current row number is i.  determine  where to go.
      ipos = iwk(i)
!!$     save the chased element.
      If (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
!!$     then occupy its location.
      If (values) a(ipos) = t
      ja(ipos) = j
!!$     update pointer information for next element to come in row i.
      iwk(i) = ipos + 1
!!$     determine  next element to be chased,
      If (ia(ipos) .Lt. 0) Go To 65
      t = tnext
      i = inext
      j = jnext
      ia(ipos) = -1
      If (k .Lt. nnz) Go To 6
      Go To 70
65    init = init + 1
      If (init .Gt. nnz) Go To 70
      If (ia(init) .Lt. 0) Go To 65
!!$     restart chasing --
      Go To 5
70    ia(2:n + 1) = iwk(1:n)
      ia(1) = 1
      Return
   End Subroutine rcoicsr

   Subroutine zcoicsr(n, nnz, job, a, ja, ia, iwk)
      Implicit None
      Integer n, nnz, ia(nnz), ja(nnz), iwk(n + 1), job
      Complex(Kind=DEF_DBL_PREC) a(*)
!!$------------------------------------------------------------------------
!!$ IN-PLACE coo-csr conversion routine.
!!$-----------------------------------------------------------------------
!!$ this subroutine converts a matrix stored in coordinate format into
!!$ the csr format. The conversion is done in place in that the arrays
!!$ a,ja,ia of the result are overwritten onto the original arrays.
!!$------------------------------------------------------------------------
!!$ on entry:
!!$---------
!!$ n        = integer. row dimension of A.
!!$ nnz        = integer. number of nonzero elements in A.
!!$ job   = integer. Job indicator. when job=1, the real values in a are
!!$         filled. Otherwise a is not touched and the structure of the
!!$         array only (i.e. ja, ia)  is obtained.
!!$ a        = real array of size nnz (number of nonzero elements in A)
!!$         containing the nonzero elements
!!$ ja        = integer array of length nnz containing the column positions
!!$           of the corresponding elements in a.
!!$ ia        = integer array of length nnz containing the row positions
!!$           of the corresponding elements in a.
!!$ iwk        = integer work array of length n+1
!!$ on return:
!!$----------
!!$ a
!!$ ja
!!$ ia        = contains the compressed sparse row data structure for the
!!$         resulting matrix.
!!$ Note:
!!$-------
!!$         the entries of the output matrix are not sorted (the column
!!$         indices in each are not in increasing order) use coocsr
!!$         if you want them sorted.
!!$----------------------------------------------------------------------c
!!$  Coded by Y. Saad, Sep. 26 1989                                      c
!!$----------------------------------------------------------------------c
      Complex(Kind=DEF_DBL_PREC) :: t, tnext
      Logical :: values
      Integer :: k, i, init, j, ipos, inext, jnext
!!$-----------------------------------------------------------------------
      values = (job .Eq. 1)
!!$ find pointer array for resulting matrix.
      iwk = 0

      Do k = 1, nnz
         i = ia(k)
         iwk(i + 1) = iwk(i + 1) + 1
      End Do
!!$------------------------------------------------------------------------
      iwk(1) = 1
      Do i = 2, n
         iwk(i) = iwk(i - 1) + iwk(i)
      End Do

!!$     loop for a cycle in chasing process.

      init = 1
      k = 0
5     If (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
!!$------------------------------------------------------------------------
6     k = k + 1
!!$     current row number is i.  determine  where to go.
      ipos = iwk(i)
!!$     save the chased element.
      If (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
!!$     then occupy its location.
      If (values) a(ipos) = t
      ja(ipos) = j
!!$     update pointer information for next element to come in row i.
      iwk(i) = ipos + 1
!!$     determine  next element to be chased,
      If (ia(ipos) .Lt. 0) Go To 65
      t = tnext
      i = inext
      j = jnext
      ia(ipos) = -1
      If (k .Lt. nnz) Go To 6
      Go To 70
65    init = init + 1
      If (init .Gt. nnz) Go To 70
      If (ia(init) .Lt. 0) Go To 65
!!$     restart chasing --
      Go To 5
70    ia(2:n + 1) = iwk(1:n)
      ia(1) = 1
      Return
   End Subroutine zcoicsr

   Subroutine zdense2coo(a, b, th)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: a(:, :)
      Type(zcoomat) :: b
      Real(Kind=DEF_DBL_PREC) :: th
!!$ Local vars
      Integer :: nr, nc, jr, jc
      nr = size(a, dim=1)
      nc = size(a, dim=2)
      Call alloc(b, nr*nc, nr, nc)
      Do jr = 1, nr
         Do jc = 1, nc
            If (Abs(a(jr, jc)) > th) Then
               b%nnz = b%nnz + 1
               b%a(b%nnz) = a(jr, jc)
               b%icol(b%nnz) = jr
               b%irow(b%nnz) = jc
            End If
         End Do
      End Do
   End Subroutine zdense2coo

   Subroutine r_ordercsr(amat)
      Implicit None
      Type(rcsrmat) :: amat
      Integer :: i, k1, k, kn, j
      Integer, Allocatable :: iwork(:)

      kn = 0
      Do i = 1, amat%nrow
         k = amat%ir(i + 1) - amat%ir(i)
         If (k > kn) kn = k
      End Do
!$omp parallel private(iwork,k,k1,i) shared (kn,amat)
      Allocate (iwork(kn))
!$omp do
      Do i = 1, amat%nrow
         k = amat%ir(i)
         k1 = amat%ir(i + 1) - 1
         kn = k1 - k + 1
         Call iqsort(kn, amat%jc(k:k1), iwork)
         amat%a(k:k1) = amat%a(k - 1 + iwork(:kn))
         amat%jc(k:k1) = amat%jc(k - 1 + iwork(:kn))
      End Do
!$omp end do
      Deallocate (iwork)
!$omp end parallel
   End Subroutine r_ordercsr

   Subroutine z_ordercsr(amat)
      Implicit None
      Type(zcsrmat) :: amat
      Integer :: i, k, k1, kn, kn0, j
      Integer, Allocatable :: iwork(:)

      kn0 = 0
      Do i = 1, amat%nrow
         k = amat%ir(i + 1) - amat%ir(i)
         If (k > kn0) kn0 = k
      End Do

         !!!!!$omp parallel private(iwork,k,k1,i) shared(kn,amat)
!$omp parallel private(iwork,k,k1,kn,i) shared(amat)
      Allocate (iwork(kn0))
!$omp do
      Do i = 1, amat%nrow
         k = amat%ir(i)
         k1 = amat%ir(i + 1) - 1
         kn = k1 - k + 1
         Call iqsort(kn, amat%jc(k:k1), iwork)
         amat%a(k:k1) = amat%a(k - 1 + iwork(:kn))
         amat%jc(k:k1) = amat%jc(k - 1 + iwork(:kn))
      End Do
!$omp end do
      Deallocate (iwork)
!$omp end parallel
   End Subroutine z_ordercsr

   Subroutine z_cooaddel(amat, vals, ic, ir, nel)
      Implicit None
      Type(zcoomat) :: amat
      Complex(Kind=DEF_DBL_PREC) :: vals(:)
      Integer :: ic(:), ir(:), nel, nnz1, nnzn

      nnz1 = amat%nnz + 1
      nnzn = amat%nnz + nel
      amat%a(nnz1:nnzn) = vals(1:nel)
      amat%icol(nnz1:nnzn) = ic(1:nel)
      amat%irow(nnz1:nnzn) = ir(1:nel)
      amat%nnz = nnzn
   End Subroutine z_cooaddel

   Subroutine r_cooaddel(amat, vals, ic, ir, nel)
      Implicit None
      Type(rcoomat) :: amat
      Real(Kind=DEF_DBL_PREC) :: vals(:)
      Integer :: ic(:), ir(:), nel, nnz1, nnzn

      nnz1 = amat%nnz + 1
      nnzn = amat%nnz + nel
      amat%a(nnz1:nnzn) = vals(1:nel)
      amat%icol(nnz1:nnzn) = ic(1:nel)
      amat%irow(nnz1:nnzn) = ir(1:nel)
      amat%nnz = nnzn
   End Subroutine r_cooaddel

   Subroutine spmatmul_size(nrow1, ncol2, ir1, jc1, ir2, jc2, nz, nir)
      Implicit None
!!$ find the size of the new matrix which is a product of a1 and a2
      Integer, Intent(In) :: nrow1, ncol2
      Integer, Intent(In), Dimension(*) :: ir1, jc1, ir2, jc2
      Integer, Intent(Inout) :: nz, nir(:)
!!$    integer , dimension (ncol2) :: mask
      Integer, Pointer :: mask(:)
      Integer :: i, j, k, neigh, icol_add

!!$ initialise the mask array which is an array that has
!!$ value i if column index already exist in row i of the new matrix

      nz = 0
      nir = 0
!$omp parallel private(mask,j,k,icol_add,neigh) shared(nrow1) reduction(+:nz)
      Allocate (mask(ncol2))
      mask = 0
!$omp do
      Do i = 1, nrow1
         Do j = ir1(i), ir1(i + 1) - 1
            neigh = jc1(j)
            Do k = ir2(neigh), ir2(neigh + 1) - 1
               icol_add = jc2(k)
               If (mask(icol_add) /= i) Then
                  nz = nz + 1
                  nir(i + 1) = nir(i + 1) + 1
                  mask(icol_add) = i ! add mask
               End If
            End Do
         End Do
      End Do
!$omp end do
      Deallocate (mask)
!$omp end parallel
      nir(1) = 1
      Do i = 1, nrow1
         nir(i + 1) = nir(i + 1) + nir(i)
      End Do
   End Subroutine spmatmul_size

!!$************* NON-INPLACE (!!!) CSR matrix - CSR matrix  multiplication functions, return new CSR matrix *****

   Function fmmcsrcsr_zz(amat, bmat, newnnz, addnnz) Result(cmat)
      Implicit None
      Type(zcsrmat) :: amat, bmat, cmat
      Integer :: nz
      Integer, Optional :: newnnz, addnnz
      Integer :: nnzm

      Call alloc(cmat, 0, amat%nrow, bmat%ncol)

      Call spmatmul_size(amat%nrow, bmat%ncol, amat%ir, amat%jc, bmat%ir, bmat%jc, nz, cmat%ir)

      If (present(newnnz)) Then
         nnzm = Max(newnnz, nz)
      Else
         If (present(addnnz)) Then
            nnzm = nz + addnnz
         Else
            nnzm = nz
         End If
      End If
      Deallocate (cmat%jc, cmat%a)
      Allocate (cmat%jc(nnzm), cmat%a(nnzm))
      cmat%nnzmax = nnzm
!!$    cmat%a=DEF_cmplx_zero
      Call matmul_normal_zz(amat%nrow, bmat%ncol, amat%ir, amat%jc, amat%a, bmat%ir, bmat%jc, bmat%a, &
     & cmat%ir, cmat%jc, cmat%a)
      cmat%nnz = nz

      Call ordercsr(cmat)

      Return
   End Function fmmcsrcsr_zz

   Function fmmcsrcsr_zr(amat, bmat, newnnz, addnnz) Result(cmat)
      Implicit None
      Type(zcsrmat) :: amat, cmat
      Type(rcsrmat) :: bmat
      Integer :: nz
      Integer, Optional :: newnnz, addnnz
      Integer :: nnzm

      Call alloc(cmat, 0, amat%nrow, bmat%ncol)
      Call spmatmul_size(amat%nrow, bmat%ncol, amat%ir, amat%jc, bmat%ir, bmat%jc, nz, cmat%ir)

      If (present(newnnz)) Then
         nnzm = Max(newnnz, nz)
      Else
         If (present(addnnz)) Then
            nnzm = nz + addnnz
         Else
            nnzm = nz
         End If
      End If

      Deallocate (cmat%jc, cmat%a)
      Allocate (cmat%jc(nnzm), cmat%a(nnzm))
      cmat%nnzmax = nnzm
      Call matmul_normal_zr(amat%nrow, bmat%ncol, amat%ir, amat%jc, amat%a, bmat%ir, bmat%jc, bmat%a, &
     & cmat%ir, cmat%jc, cmat%a)
      cmat%nnz = nz

      Call ordercsr(cmat)

      Return
   End Function fmmcsrcsr_zr

   Function fmmcsrcsr_rz(amat, bmat, newnnz, addnnz) Result(cmat)
      Implicit None
      Type(zcsrmat) :: bmat, cmat
      Type(rcsrmat) :: amat
      Integer :: nz
      Integer, Optional :: newnnz, addnnz
      Integer :: nnzm

      Call alloc(cmat, 0, amat%nrow, bmat%ncol)
      Call spmatmul_size(amat%nrow, bmat%ncol, amat%ir, amat%jc, bmat%ir, bmat%jc, nz, cmat%ir)

      If (present(newnnz)) Then
         nnzm = Max(newnnz, nz)
      Else
         If (present(addnnz)) Then
            nnzm = nz + addnnz
         Else
            nnzm = nz
         End If
      End If

      Deallocate (cmat%jc, cmat%a)
      Allocate (cmat%jc(nnzm), cmat%a(nnzm))
      cmat%nnzmax = nnzm
      Call matmul_normal_rz(amat%nrow, bmat%ncol, amat%ir, amat%jc, amat%a, bmat%ir, bmat%jc, bmat%a, &
     & cmat%ir, cmat%jc, cmat%a)
      cmat%nnz = nz

      Call ordercsr(cmat)

      Return
   End Function fmmcsrcsr_rz

   Subroutine matmul_normal_zz(n, m, ia, ja, a, ib, jb, b, ic, jc, c)
      Implicit None
!!$ C = A*B
!!$    usage:  matmul_normal(nrow1,ncol2,ir1,jc1,val1,ir2,jc2,val2,ir3,jc3,val3)
      Integer, Intent(In) :: n, m
      Integer, Intent(In) :: ia(:), ib(:), ja(:), jb(:), ic(:)
      Integer, Intent(Out) :: jc(:)
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: a(:), b(:)
      Complex(Kind=DEF_DBL_PREC), Intent(Out) :: c(:)
      Integer :: nz, i, j, k, icol, icol_add, neigh
      Integer, Pointer :: mask(:)
      Complex(Kind=DEF_DBL_PREC) :: aij

!!$ initialise the mask array which is an array that has non-zero value if the
!!$ column index already exist, in which case
!!$ the value is the index of that column

!$omp parallel private(mask,j,k,icol_add,icol,aij,nz,neigh) shared(n)
      Allocate (mask(m))
      mask = 0
!$omp do
      Do i = 1, n
         If (ic(i + 1) /= ic(i)) Then
            nz = ic(i)
            Do j = ia(i), ia(i + 1) - 1
               aij = a(j)
               neigh = ja(j)
               Do k = ib(neigh), ib(neigh + 1) - 1
                  icol_add = jb(k)
                  icol = mask(icol_add)
                  If (icol == 0) Then
                     jc(nz) = icol_add
                     c(nz) = aij*b(k)
                     mask(icol_add) = nz ! add mask
                     nz = nz + 1
                  Else
                     c(icol) = c(icol) + aij*b(k)
                  End If
               End Do
            End Do
            mask(jc(ic(i):nz - 1)) = 0 ! done this row i, so set mask to zero again
         End If
      End Do
!$omp end do
      Deallocate (mask)
!$omp end parallel
   End Subroutine matmul_normal_zz

   Subroutine matmul_normal_zr(n, m, ia, ja, a, ib, jb, b, ic, jc, c)
      Implicit None
!!$ C = A*B
!!$    usage:  matmul_normal(nrow1,ncol2,ir1,jc1,val1,ir2,jc2,val2,ir3,jc3,val3)
      Integer, Intent(In) :: n, m
      Integer, Intent(In) :: ia(:), ib(:), ja(:), jb(:), ic(:)
      Integer, Intent(Out) :: jc(:)
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: a(:)
      Real(Kind=DEF_DBL_PREC), Intent(In) :: b(:)
      Complex(Kind=DEF_DBL_PREC), Intent(Out) :: c(:)
      Integer :: nz, i, j, k, icol, icol_add, neigh
      Integer, Pointer :: mask(:)
      Complex(Kind=DEF_DBL_PREC) :: aij

!!$ initialise the mask array which is an array that has non-zero value if the
!!$ column index already exist, in which case
!!$ the value is the index of that column

!$omp parallel private(mask,j,k,icol_add,icol,aij,nz,neigh)  shared(n)
      Allocate (mask(m))
      mask = 0
!$omp do
      Do i = 1, n
         If (ic(i + 1) /= ic(i)) Then
            nz = ic(i)
            Do j = ia(i), ia(i + 1) - 1
               aij = a(j)
               neigh = ja(j)
               Do k = ib(neigh), ib(neigh + 1) - 1
                  icol_add = jb(k)
                  icol = mask(icol_add)
                  If (icol == 0) Then
                     jc(nz) = icol_add
                     c(nz) = aij*b(k)
                     mask(icol_add) = nz ! add mask
                     nz = nz + 1
                  Else
                     c(icol) = c(icol) + aij*b(k)
                  End If
               End Do
            End Do
            mask(jc(ic(i):nz - 1)) = 0 ! done this row i, so set mask to zero again
         End If
      End Do
!$omp end do
      Deallocate (mask)
!$omp end parallel
   End Subroutine matmul_normal_zr

   Subroutine matmul_normal_rz(n, m, ia, ja, a, ib, jb, b, ic, jc, c)
      Implicit None
!!$   C = A*B
!!$   usage:  matmul_normal(nrow1,ncol2,ir1,jc1,val1,ir2,jc2,val2,ir3,jc3,val3)
      Integer, Intent(In) :: n, m
      Integer, Intent(In) :: ia(:), ib(:), ja(:), jb(:), ic(:)
      Integer, Intent(Out) :: jc(:)
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: b(:)
      Real(Kind=DEF_DBL_PREC), Intent(In) :: a(:)
      Complex(Kind=DEF_DBL_PREC), Intent(Out) :: c(:)
      Integer :: nz, i, j, k, icol, icol_add, neigh
      Integer, Pointer :: mask(:)
      Real(Kind=DEF_DBL_PREC) :: aij

!!$ initialise the mask array which is an array that has non-zero value if the
!!$ column index already exist, in which case
!!$ the value is the index of that column

!$omp parallel private(mask,j,k,icol_add,icol,aij,nz,neigh)  shared(n)
      Allocate (mask(m))
      mask = 0
!$omp do
      Do i = 1, n
         If (ic(i + 1) /= ic(i)) Then
            nz = ic(i)
            Do j = ia(i), ia(i + 1) - 1
               aij = a(j)
               neigh = ja(j)
               Do k = ib(neigh), ib(neigh + 1) - 1
                  icol_add = jb(k)
                  icol = mask(icol_add)
                  If (icol == 0) Then
                     jc(nz) = icol_add
                     c(nz) = aij*b(k)
                     mask(icol_add) = nz ! add mask
                     nz = nz + 1
                  Else
                     c(icol) = c(icol) + aij*b(k)
                  End If
               End Do
            End Do
            mask(jc(ic(i):nz - 1)) = 0 ! done this row i, so set mask to zero again
         End If
      End Do
!$omp end do
      Deallocate (mask)
!$omp end parallel
   End Subroutine matmul_normal_rz

!!$************* NON-INPLACE (!!!) 3 matrix product, CSR matrix * diag * CSR matrix  functions, return new CSR matrix *****

   Function fmmcsrdiagcsr_zrz(amat, dmat, bmat, newnnz, addnnz) Result(rmat)
      Implicit None
      Type(zcsrmat) :: amat, bmat, rmat
      Type(rdiagmat) :: dmat
      Integer :: nz
      Integer, Optional :: newnnz, addnnz
      Integer :: nnzm

      Call alloc(rmat, 0, amat%nrow, bmat%ncol)
      Call spmatmul_size(amat%nrow, bmat%ncol, amat%ir, amat%jc, bmat%ir, bmat%jc, nz, rmat%ir)

      If (present(newnnz)) Then
         nnzm = Max(newnnz, nz)
      Else
         If (present(addnnz)) Then
            nnzm = nz + addnnz
         Else
            nnzm = nz
         End If
      End If

!!$     write(*,*) 'stp-2'
      Deallocate (rmat%jc, rmat%a)
      Allocate (rmat%jc(nnzm), rmat%a(nnzm))
      rmat%nnzmax = nnzm
!!$     write(*,*) 'stp-3'
      Call matmul_withdiag_zrz(amat%nrow, bmat%ncol, amat%ir, amat%jc, amat%a, bmat%ir, bmat%jc, bmat%a, &
     & rmat%ir, rmat%jc, rmat%a, dmat%el)
      rmat%nnz = nz

      Call ordercsr(rmat)
!!$     write(*,*) 'stp-4'

      Return
   End Function fmmcsrdiagcsr_zrz

   Subroutine matmul_withdiag_zrz(n, m, ia, ja, a, ib, jb, b, ic, jc, c, diag)
      Implicit None
!!$ C = A*B
!!$    usage:  matmul_normal(nrow1,ncol2,ir1,jc1,val1,ir2,jc2,val2,ir3,jc3,val3)
      Integer, Intent(In) :: n, m
      Integer, Intent(In) :: ia(:), ib(:), ja(:), jb(:), ic(:)
      Integer, Intent(Out) :: jc(:)
      Complex(Kind=DEF_DBL_PREC), Intent(In) :: a(:), b(:)
      Complex(Kind=DEF_DBL_PREC), Intent(Out) :: c(:)
      Integer :: nz, i, j, k, icol, icol_add, neigh
      Integer, Pointer :: mask(:)
      Complex(Kind=DEF_DBL_PREC) :: aij
      Real(Kind=DEF_DBL_PREC) :: diag(:)

!!$ initialise the mask array which is an array that has non-zero value if the
!!$ column index already exist, in which case
!!$ the value is the index of that column

!$omp parallel private(mask,j,k,icol_add,icol,aij,nz,neigh)  shared(n)
      Allocate (mask(m))
      mask = 0
!$omp do
      Do i = 1, n
         If (ic(i + 1) /= ic(i)) Then
            nz = ic(i)
            Do j = ia(i), ia(i + 1) - 1
               neigh = ja(j)
               aij = a(j)*diag(neigh)
!!$            aij = a(j)
               Do k = ib(neigh), ib(neigh + 1) - 1
                  icol_add = jb(k)
                  icol = mask(icol_add)
                  If (icol == 0) Then
                     jc(nz) = icol_add
                     c(nz) = aij*b(k)
                     mask(icol_add) = nz ! add mask
                     nz = nz + 1
                  Else
                     c(icol) = c(icol) + aij*b(k)
                  End If
               End Do
            End Do
            mask(jc(ic(i):nz - 1)) = 0 ! done this row i, so set mask to zero again
         End If
      End Do
!$omp end do
      Deallocate (mask)
!$omp end parallel
   End Subroutine matmul_withdiag_zrz

!!$**************** COPY functions, create new CSR matrix A, copy of old one B
   Subroutine spcopy_dense_z(amat, bmat, sign)
      Implicit None
      Type(zdensemat) :: amat, bmat
      Integer, Optional :: sign
      Integer :: sg

      sg = 1
      If (present(sign)) sg = sign

      Call alloc(bmat, amat%nrow, amat%ncol)

      If (sg > 0) Then
         bmat%bl(:, :) = amat%bl(1:amat%nrow, 1:amat%ncol)
      Else
         bmat%bl(:, :) = -amat%bl(1:amat%nrow, 1:amat%ncol)
      End If
   End Subroutine spcopy_dense_z

   Subroutine spcopy_csr_z(amat, bmat, newnnz, addnnz, sign)
      Implicit None
      Type(zcsrmat) :: amat, bmat
      Integer, Optional :: newnnz, addnnz, sign
      Integer :: nnz, nnzm, sg

      sg = 1
      If (present(sign)) sg = sign
      If (present(newnnz)) Then
         nnzm = Max(newnnz, amat%nnz)
      Else
         If (present(addnnz)) Then
            nnzm = amat%nnz + addnnz
         Else
            nnzm = amat%nnz
         End If
      End If
!!$    if (bmat%alloc==1) call free(bmat)
      Call alloc(bmat, nnzm, amat%nrow, amat%ncol)

      nnz = amat%nnz

      If (sg > 0) Then
         bmat%a(1:nnz) = amat%a(1:nnz)
      Else
         bmat%a(1:nnz) = -amat%a(1:nnz)
      End If
      bmat%ir = amat%ir
      bmat%jc(1:nnz) = amat%jc(1:nnz)
      bmat%nnz = amat%nnz
   End Subroutine spcopy_csr_z

   Subroutine spcopy_csr_r(amat, bmat, newnnz, addnnz)
      Implicit None
      Type(rcsrmat) :: amat, bmat
      Integer, Optional :: newnnz, addnnz
      Integer :: nnz, nnzm

      If (present(newnnz)) Then
         nnzm = Max(newnnz, amat%nnz)
      Else
         If (present(addnnz)) Then
            nnzm = amat%nnz + addnnz
         Else
            nnzm = amat%nnz
         End If
      End If
!!$    if (bmat%alloc==1) call free(bmat)
      Call alloc(bmat, nnzm, amat%nrow, amat%ncol)

      nnz = amat%nnz

      bmat%a(1:nnz) = amat%a(1:nnz)
      bmat%ir = amat%ir
      bmat%jc(1:nnz) = amat%jc(1:nnz)
      bmat%nnz = amat%nnz
   End Subroutine spcopy_csr_r

   Subroutine spcopy_csr_rz(amat, bmat, newnnz, addnnz)
      Implicit None
      Type(rcsrmat) :: amat
      Type(zcsrmat) :: bmat
      Integer, Optional :: newnnz, addnnz
      Integer :: nnz, nnzm

      If (present(newnnz)) Then
         nnzm = Max(newnnz, amat%nnz)
      Else
         If (present(addnnz)) Then
            nnzm = amat%nnz + addnnz
         Else
            nnzm = amat%nnz
         End If
      End If
!!$    if (bmat%alloc==1) call free(bmat)
      Call alloc(bmat, nnzm, amat%nrow, amat%ncol)

      nnz = amat%nnz

      bmat%a(1:nnz) = amat%a(1:nnz)
      bmat%ir = amat%ir
      bmat%jc(1:nnz) = amat%jc(1:nnz)
      bmat%nnz = amat%nnz
   End Subroutine spcopy_csr_rz

   Function csrherm_z(amat) Result(bmat)
      Implicit None
      Type(zcsrmat) :: amat, bmat
!!$          Call init(bmat)
      Call alloc(bmat, amat%nnz, amat%nrow, amat%ncol)
      Call csrcsc2z(amat%nrow, amat%ncol, amat%a, amat%jc, amat%ir, bmat%a, bmat%jc, bmat%ir)
      bmat%a(1:amat%nnz) = conjg(bmat%a(1:amat%nnz))
      bmat%nnz = amat%nnz
      Call ordercsr(bmat)
      Return
   End Function csrherm_z

   Function csrtransp_z(amat) Result(bmat)
      Implicit None
      Type(zcsrmat) :: amat, bmat
      Call alloc(bmat, amat%nnz, amat%nrow, amat%ncol)
      Call csrcsc2z(amat%nrow, amat%ncol, amat%a, amat%jc, amat%ir, bmat%a, bmat%jc, bmat%ir)
      bmat%nnz = amat%nnz
      Call ordercsr(bmat)
      Return
   End Function csrtransp_z

   Function csrherm_r(amat) Result(bmat)
      Implicit None
      Type(rcsrmat) :: amat, bmat
      Call alloc(bmat, amat%nnz, amat%nrow, amat%ncol)
      Call ordercsr(amat)
      Call csrcsc2r(amat%nrow, amat%ncol, amat%a, amat%jc, amat%ir, bmat%a, bmat%jc, bmat%ir)
!!$      bmat%a=conjg(bmat%a)
      bmat%nnz = amat%nnz
      Call ordercsr(bmat)
      Return
   End Function csrherm_r

   Subroutine csrcsc2z(n, n2, a, ja, ia, ao, jao, iao)
      Implicit None
      Integer :: i, j, k, next, n, n2
      Integer :: ia(n + 1), iao(n2 + 1), ja(:), jao(:)
      Complex(Kind=DEF_DBL_PREC) :: a(:), ao(:)

      Do i = 1, n2 + 1
         iao(i) = 0
      End Do
      Do i = 1, n
         Do k = ia(i), ia(i + 1) - 1, 1
            j = ja(k) + 1
            iao(j) = iao(j) + 1
         End Do
      End Do

!!$---------- compute pointers from lengths ------------------------------
      iao(1) = 1
      Do i = 1, n2
         iao(i + 1) = iao(i) + iao(i + 1)
      End Do
!!$--------------- now do the actual copying -----------------------------

      Do i = 1, n
         Do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            next = iao(j)
            ao(next) = a(k)
            jao(next) = i
            iao(j) = next + 1
         End Do
      End Do

!!$-------------------------- reshift iao and leave ----------------------
      Do i = n2, 1, -1
         iao(i + 1) = iao(i)
      End Do
      iao(1) = 1
   End Subroutine csrcsc2z

   Subroutine csrcsc2r(n, n2, a, ja, ia, ao, jao, iao)
      Implicit None
      Integer :: i, j, k, next, n, n2
      Integer :: ia(n + 1), iao(n2 + 1), ja(:), jao(:)
      Real(Kind=DEF_DBL_PREC) :: a(:), ao(:)

      Do i = 1, n2 + 1
         iao(i) = 0
      End Do
      Do i = 1, n
         Do k = ia(i), ia(i + 1) - 1
            j = ja(k) + 1
            iao(j) = iao(j) + 1
         End Do
      End Do

!!$---------- compute pointers from lengths ------------------------------
      iao(1) = 1
      Do i = 1, n2
         iao(i + 1) = iao(i) + iao(i + 1)
      End Do
!!$--------------- now do the actual copying -----------------------------
      Do i = 1, n
         Do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            next = iao(j)
            ao(next) = a(k)
            jao(next) = i
            iao(j) = next + 1
         End Do
      End Do
!!$-------------------------- reshift iao and leave ----------------------
      Do i = n2, 1, -1
         iao(i + 1) = iao(i)
      End Do
      iao(1) = 1
   End Subroutine csrcsc2r

   Function fcsraddcsr_zz(amat, bmat, newnnz, addnnz, sign) Result(rmat)
      Implicit None
      Type(zcsrmat) :: amat, bmat, rmat
      Integer, Optional :: newnnz, addnnz, sign
!!$ Local
      Integer :: nnzm, nz, sg
!!$Code

      Call alloc(rmat, 0, amat%nrow, amat%ncol)

      Call matadd_size(amat%nrow, amat%ncol, amat%jc, amat%ir, bmat%jc, bmat%ir, nz, rmat%ir)

      If (present(sign)) Then
         sg = sign
      Else
         sg = 1
      End If
      If (present(newnnz)) Then
         nnzm = Max(newnnz, nz)
      Else
         If (present(addnnz)) Then
            nnzm = nz + addnnz
         Else
            nnzm = nz
         End If
      End If

      Deallocate (rmat%jc, rmat%a)
      Allocate (rmat%jc(nnzm), rmat%a(nnzm))
      rmat%nnzmax = nnzm

      If (sg > 0) Then
         Call aplb1_zz(amat%nrow, amat%ncol, amat%a, amat%jc, amat%ir, bmat%a, bmat%jc, bmat%ir, rmat%a, &
        & rmat%jc, rmat%ir)
      Else
         Call amib1_zz(amat%nrow, amat%ncol, amat%a, amat%jc, amat%ir, bmat%a, bmat%jc, bmat%ir, rmat%a, &
        & rmat%jc, rmat%ir)
      End If
      rmat%nnz = nz
   End Function fcsraddcsr_zz

   Function fcsraddcsr_zr(amat, bmat, newnnz, addnnz) Result(rmat)
      Implicit None
      Type(zcsrmat) :: amat, rmat
      Type(rcsrmat) :: bmat
      Integer, Optional :: newnnz, addnnz
      Integer :: nnzm, nz

!!$      nz=amat%nnz+bmat%nnz           !must be fixed here....FIXED

      Call alloc(rmat, 0, amat%nrow, amat%ncol)
      Call matadd_size(amat%nrow, amat%ncol, amat%jc, amat%ir, bmat%jc, bmat%ir, nz, rmat%ir)

      If (present(newnnz)) Then
         nnzm = Max(newnnz, nz)
      Else
         If (present(addnnz)) Then
            nnzm = nz + addnnz
         Else
            nnzm = nz
         End If
      End If

      Deallocate (rmat%jc, rmat%a)
      Allocate (rmat%jc(nnzm), rmat%a(nnzm))
      rmat%nnzmax = nnzm

      Call aplb1_zr(amat%nrow, amat%ncol, amat%a, amat%jc, amat%ir, bmat%a, bmat%jc, bmat%ir, rmat%a, &
     & rmat%jc, rmat%ir)
      rmat%nnz = nz
   End Function fcsraddcsr_zr

   Subroutine matadd_size(nrow, ncol, ja, ia, jb, ib, nz, nir)
      Implicit None
      Integer :: nrow, i, ka, kb, j1, j2, kamax, kbmax, ncol
      Integer :: ja(:), jb(:), ia(nrow + 1), ib(nrow + 1)
      Integer, Intent(Inout) :: nz, nir(:)
      Logical :: flag
      nz = 0
      nir(:) = 0
!$omp parallel private(ka,kb,kamax,kbmax,j1,j2,flag) shared(nir,ia,ib,ja,jb,nrow,ncol) reduction(+:nz)
!$omp do
      Do i = 1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i + 1) - 1
         kbmax = ib(i + 1) - 1
         flag = ((ka <= kamax) .Or. (kb <= kbmax))
!!!$             Do while ((ka <= kamax) .Or. (kb <= kbmax))
         Do while (flag)
            If (ka <= kamax) Then
               j1 = ja(ka)
            Else
               j1 = ncol + 1
            End If

            If (kb <= kbmax) Then
               j2 = jb(kb)
            Else
               j2 = ncol + 1
            End If

            If (j1 == j2) Then
               ka = ka + 1
               kb = kb + 1
               nz = nz + 1
            Else If (j1 < j2) Then
               ka = ka + 1
               nz = nz + 1
            Else If (j1 > j2) Then
               kb = kb + 1
               nz = nz + 1
            End If
            nir(i + 1) = nir(i + 1) + 1
            flag = ((ka <= kamax) .Or. (kb <= kbmax))
         End Do
      End Do
!$omp end do
!$omp end parallel

      nir(1) = 1
      Do i = 1, nrow
         nir(i + 1) = nir(i + 1) + nir(i)
      End Do
      Return

   End Subroutine matadd_size

   Subroutine aplb1_zz(nrow, ncol, a, ja, ia, b, jb, ib, c, jc, ic)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: a(:), b(:), c(:)
      Integer :: i, nrow, ncol, ka, kb, kc, j1, j2, kamax, kbmax
      Integer :: ja(:), jb(:), jc(:), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1)
      Logical :: flag

!$omp parallel private(ka,kb,kamax,kbmax,j1,j2,kc,flag)
!$omp do
      Do i = 1, nrow
         kc = ic(i)
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i + 1) - 1
         kbmax = ib(i + 1) - 1
         flag = ((ka <= kamax) .Or. (kb <= kbmax))
!!$             Do while ((ka <= kamax) .Or. (kb <= kbmax))
         Do while (flag)
            If (ka <= kamax) Then
               j1 = ja(ka)
            Else
               j1 = ncol + 1
            End If
            If (kb <= kbmax) Then
               j2 = jb(kb)
            Else
               j2 = ncol + 1
            End If

!!$     three cases

            If (j1 == j2) Then
               c(kc) = a(ka) + b(kb)
               jc(kc) = j1
               ka = ka + 1
               kb = kb + 1
            Else If (j1 < j2) Then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka + 1
            Else If (j1 > j2) Then
               jc(kc) = j2
               c(kc) = b(kb)
               kb = kb + 1
            End If
            kc = kc + 1
            flag = ((ka <= kamax) .Or. (kb <= kbmax))
         End Do
!!$            ic (i+1) = kc
      End Do
!$omp end do
!$omp end parallel

      Return
   End Subroutine aplb1_zz

   Subroutine amib1_zz(nrow, ncol, a, ja, ia, b, jb, ib, c, jc, ic)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: a(:), b(:), c(:)
      Integer :: i, nrow, ncol, ka, kb, kc, j1, j2, kamax, kbmax
      Integer :: ja(:), jb(:), jc(:), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1)
      Logical :: flag

!$omp parallel private(ka,kb,kamax,kbmax,j1,j2,kc,flag)
!$omp do
      Do i = 1, nrow
         kc = ic(i)
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i + 1) - 1
         kbmax = ib(i + 1) - 1
         flag = ((ka <= kamax) .Or. (kb <= kbmax))
!!$             Do while ((ka <= kamax) .Or. (kb <= kbmax))
         Do while (flag)
            If (ka <= kamax) Then
               j1 = ja(ka)
            Else
               j1 = ncol + 1
            End If
            If (kb <= kbmax) Then
               j2 = jb(kb)
            Else
               j2 = ncol + 1
            End If

!!$     three cases

            If (j1 == j2) Then
               c(kc) = a(ka) - b(kb)
               jc(kc) = j1
               ka = ka + 1
               kb = kb + 1
            Else If (j1 < j2) Then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka + 1
            Else If (j1 > j2) Then
               jc(kc) = j2
               c(kc) = -b(kb)
               kb = kb + 1
            End If
            kc = kc + 1
            flag = ((ka <= kamax) .Or. (kb <= kbmax))
         End Do
      End Do
!$omp end do
!$omp end parallel

      Return
   End Subroutine amib1_zz

   Subroutine amib1_zr(nrow, ncol, a, ja, ia, b, jb, ib, c, jc, ic)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: a(:), c(:)
      Real(Kind=DEF_DBL_PREC) :: b(:)
      Integer :: i, nrow, ncol, ka, kb, kc, j1, j2, kamax, kbmax
      Integer :: ja(:), jb(:), jc(:), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1)
      Logical :: flag

!$omp parallel private(ka,kb,kamax,kbmax,j1,j2,kc,flag)
!$omp do
      Do i = 1, nrow
         kc = ic(i)
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i + 1) - 1
         kbmax = ib(i + 1) - 1
         flag = ((ka <= kamax) .Or. (kb <= kbmax))
!!$             Do while ((ka <= kamax) .Or. (kb <= kbmax))
         Do while (flag)
            If (ka <= kamax) Then
               j1 = ja(ka)
            Else
               j1 = ncol + 1
            End If
            If (kb <= kbmax) Then
               j2 = jb(kb)
            Else
               j2 = ncol + 1
            End If

!!$     three cases

            If (j1 == j2) Then
               c(kc) = a(ka) - b(kb)
               jc(kc) = j1
               ka = ka + 1
               kb = kb + 1
            Else If (j1 < j2) Then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka + 1
            Else If (j1 > j2) Then
               jc(kc) = j2
               c(kc) = -b(kb)
               kb = kb + 1
            End If
            kc = kc + 1
            flag = ((ka <= kamax) .Or. (kb <= kbmax))
         End Do
      End Do
!$omp end do
!$omp end parallel

      Return
   End Subroutine amib1_zr

   Subroutine aplb1_zr(nrow, ncol, a, ja, ia, b, jb, ib, c, jc, ic)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: a(:), c(:)
      Real(Kind=DEF_DBL_PREC) :: b(:)
      Integer :: i, nrow, ncol, ka, kb, kc, j1, j2, kamax, kbmax
      Integer :: ja(:), jb(:), jc(:), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1)
      Logical :: flag

!$omp parallel private(ka,kb,kamax,kbmax,j1,j2,kc,flag)
!$omp do
      Do i = 1, nrow
         kc = ic(i)
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i + 1) - 1
         kbmax = ib(i + 1) - 1
         flag = ((ka <= kamax) .Or. (kb <= kbmax))
!!$             Do while ((ka <= kamax) .Or. (kb <= kbmax))
         Do while (flag)
            If (ka <= kamax) Then
               j1 = ja(ka)
            Else
               j1 = ncol + 1
            End If
            If (kb <= kbmax) Then
               j2 = jb(kb)
            Else
               j2 = ncol + 1
            End If

!!$     three cases

            If (j1 == j2) Then
               c(kc) = a(ka) + b(kb)
               jc(kc) = j1
               ka = ka + 1
               kb = kb + 1
            Else If (j1 < j2) Then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka + 1
            Else If (j1 > j2) Then
               jc(kc) = j2
               c(kc) = b(kb)
               kb = kb + 1
            End If
            kc = kc + 1
            flag = ((ka <= kamax) .Or. (kb <= kbmax))
         End Do
      End Do
!$omp end do
!$omp end parallel

      Return
   End Subroutine aplb1_zr

   Subroutine denseaddcsr_zz(amat, bmat)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: amat(:, :)
      Type(zcsrmat) :: bmat
      Integer :: i1, i2, j1

!!$      write(*,*) size(amat,1), bmat%nrow

      Do i1 = 1, bmat%nrow
         Do i2 = bmat%ir(i1), bmat%ir(i1 + 1) - 1
            j1 = bmat%jc(i2)
            amat(i1, j1) = amat(i1, j1) + bmat%a(i2)
         End Do
      End Do

      Return
   End Subroutine denseaddcsr_zz

   Subroutine iqsort(n, DATA, INDEX)
!!$Quick-Sort algorithm with improvements.
!!$Leaves "DATA" unchanged and returns array "INDEX" with proper ordering
      Implicit None
      Integer :: n, INDEX(n)
      Integer :: DATA(n)

      Integer :: LSTK(31), RSTK(31), ISTK
      Integer :: L, R, i, j, P, INDEXP, INDEXT
      Integer :: DATAP

      Integer, Parameter :: m = 9

      Do i = 1, n
         INDEX(i) = i
      End Do

      If (n .Le. m) Go To 910

      ISTK = 0
      L = 1
      R = n

200   Continue

      i = L
      j = R

      P = (L + R)/2
      INDEXP = INDEX(P)
      DATAP = DATA(INDEXP)

      If (DATA(INDEX(L)) .Gt. DATAP) Then
         INDEX(P) = INDEX(L)
         INDEX(L) = INDEXP
         INDEXP = INDEX(P)
         DATAP = DATA(INDEXP)
      End If

      If (DATAP .Gt. DATA(INDEX(R))) Then
         If (DATA(INDEX(L)) .Gt. DATA(INDEX(R))) Then
            INDEX(P) = INDEX(L)
            INDEX(L) = INDEX(R)
         Else
            INDEX(P) = INDEX(R)
         End If
         INDEX(R) = INDEXP
         INDEXP = INDEX(P)
         DATAP = DATA(INDEXP)
      End If

300   Continue

      i = i + 1
      If (DATA(INDEX(i)) .Lt. DATAP) Go To 300

400   Continue

      j = j - 1
      If (DATA(INDEX(j)) .Gt. DATAP) Go To 400

      If (i .Lt. j) Then
         INDEXT = INDEX(i)
         INDEX(i) = INDEX(j)
         INDEX(j) = INDEXT
         Go To 300
      Else

         If (R - j .Ge. i - L .And. i - L .Gt. m) Then
            ISTK = ISTK + 1
            LSTK(ISTK) = j + 1
            RSTK(ISTK) = R
            R = i - 1
         Else If (i - L .Gt. R - j .And. R - j .Gt. m) Then
            ISTK = ISTK + 1
            LSTK(ISTK) = L
            RSTK(ISTK) = i - 1
            L = j + 1
         Else If (R - j .Gt. m) Then
            L = j + 1
         Else If (i - L .Gt. m) Then
            R = i - 1
         Else
            If (ISTK .Lt. 1) Go To 900
            L = LSTK(ISTK)
            R = RSTK(ISTK)
            ISTK = ISTK - 1
         End If
         Go To 200
      End If

900   Continue
!!$      return
910   Continue

      Do i = 2, n
         If (DATA(INDEX(i - 1)) .Gt. DATA(INDEX(i))) Then
            INDEXP = INDEX(i)
            DATAP = DATA(INDEXP)
            P = i - 1
920         Continue
            INDEX(P + 1) = INDEX(P)
            P = P - 1
            If (P .Gt. 0) Then
               If (DATA(INDEX(P)) .Gt. DATAP) Go To 920
            End If
            INDEX(P + 1) = INDEXP
         End If
      End Do

   End Subroutine iqsort

   Subroutine dqsort(n, DATA, INDEX)
!!$Quick-Sort algorithm with improvements.
!!$Leaves "DATA" unchanged and returns array "INDEX" with proper ordering
      Implicit None
      Integer :: n, INDEX(n)
      Real(Kind=DEF_DBL_PREC) :: DATA(n)

      Integer :: LSTK(31), RSTK(31), ISTK
      Integer :: L, R, i, j, P, INDEXP, INDEXT
      Real(Kind=DEF_DBL_PREC) :: DATAP

      Integer, Parameter :: m = 9

      Do i = 1, n
         INDEX(i) = i
      End Do

      If (n .Le. m) Go To 910

      ISTK = 0
      L = 1
      R = n

200   Continue

      i = L
      j = R

      P = (L + R)/2
      INDEXP = INDEX(P)
      DATAP = DATA(INDEXP)

      If (DATA(INDEX(L)) .Gt. DATAP) Then
         INDEX(P) = INDEX(L)
         INDEX(L) = INDEXP
         INDEXP = INDEX(P)
         DATAP = DATA(INDEXP)
      End If

      If (DATAP .Gt. DATA(INDEX(R))) Then
         If (DATA(INDEX(L)) .Gt. DATA(INDEX(R))) Then
            INDEX(P) = INDEX(L)
            INDEX(L) = INDEX(R)
         Else
            INDEX(P) = INDEX(R)
         End If
         INDEX(R) = INDEXP
         INDEXP = INDEX(P)
         DATAP = DATA(INDEXP)
      End If

300   Continue

      i = i + 1
      If (DATA(INDEX(i)) .Lt. DATAP) Go To 300

400   Continue

      j = j - 1
      If (DATA(INDEX(j)) .Gt. DATAP) Go To 400

      If (i .Lt. j) Then
         INDEXT = INDEX(i)
         INDEX(i) = INDEX(j)
         INDEX(j) = INDEXT
         Go To 300
      Else

         If (R - j .Ge. i - L .And. i - L .Gt. m) Then
            ISTK = ISTK + 1
            LSTK(ISTK) = j + 1
            RSTK(ISTK) = R
            R = i - 1
         Else If (i - L .Gt. R - j .And. R - j .Gt. m) Then
            ISTK = ISTK + 1
            LSTK(ISTK) = L
            RSTK(ISTK) = i - 1
            L = j + 1
         Else If (R - j .Gt. m) Then
            L = j + 1
         Else If (i - L .Gt. m) Then
            R = i - 1
         Else
            If (ISTK .Lt. 1) Go To 900
            L = LSTK(ISTK)
            R = RSTK(ISTK)
            ISTK = ISTK - 1
         End If
         Go To 200
      End If

900   Continue
!!$      return
910   Continue

      Do i = 2, n
         If (DATA(INDEX(i - 1)) .Gt. DATA(INDEX(i))) Then
            INDEXP = INDEX(i)
            DATAP = DATA(INDEXP)
            P = i - 1
920         Continue
            INDEX(P + 1) = INDEX(P)
            P = P - 1
            If (P .Gt. 0) Then
               If (DATA(INDEX(P)) .Gt. DATAP) Go To 920
            End If
            INDEX(P + 1) = INDEXP
         End If
      End Do

   End Subroutine dqsort

!!$ This routine create new matrix "ham" where matrices "ham1" and "ham2" are diagonal blocks
   Function collectspmat(ham1, ham2) Result(ham)
      Implicit None
      Type(zcsrmat) :: ham1, ham2, ham
      Integer :: nnz1, nnz2, nnz
      nnz1 = ham1%nnz
      nnz2 = ham2%nnz
      nnz = nnz1 + nnz2
      Call alloc(ham, nnz, ham1%nrow + ham2%nrow)
      ham%a(1:nnz1) = ham1%a(1:nnz1)
      ham%jc(1:nnz1) = ham1%jc(1:nnz1)
      ham%ir(1:ham1%nrow + 1) = ham1%ir(1:ham1%nrow + 1)
      ham%a(nnz1 + 1:nnz) = ham2%a(1:nnz2)
      ham%jc(nnz1 + 1:nnz) = ham2%jc(1:nnz2) + ham1%ncol
      ham%ir(ham1%nrow + 2:ham1%nrow + ham2%nrow + 1) = ham2%ir(2:ham2%nrow + 1) + nnz1
      ham%nnz = nnz
   End Function collectspmat

!!$ This routine create new matrix "ham" where matrices "ham1" and "ham2" are OFFdiagonal blocks
   Function collectspmat_offdiag(ham1, ham2) Result(ham)
      Implicit None
      Type(zcsrmat) :: ham1, ham2, ham
      Integer :: nnz1, nnz2, nnz
      nnz1 = ham1%nnz
      nnz2 = ham2%nnz
      nnz = nnz1 + nnz2
      Call alloc(ham, nnz, ham1%nrow + ham2%nrow)
      ham%a(1:nnz1) = ham1%a(1:nnz1)
      ham%jc(1:nnz1) = ham1%jc(1:nnz1) + ham1%ncol
      ham%ir(1:ham1%nrow) = ham1%ir(1:ham1%nrow)
      ham%a(nnz1 + 1:nnz) = ham2%a(1:nnz2)
      ham%jc(nnz1 + 1:nnz) = ham2%jc(1:nnz2)
      ham%ir(ham1%nrow + 1:ham1%nrow + ham2%nrow + 1) = ham2%ir(1:ham2%nrow + 1) + nnz1
      ham%nnz = nnz
   End Function collectspmat_offdiag

   Subroutine zcsr_add_dupes(amat)
!!$ This routine adds up any duplicate entries in an ORDERED zcsr matrix
      Implicit none
      Type(zcsrmat) :: amat
       !!$ local
      Type(zcsrmat) :: bmat
      Integer :: i, j, bmat_nnz

      call alloc(bmat, amat%nnz, amat%nrow)
      bmat_nnz = 1
      do i = 1, amat%nrow
         do j = amat%ir(i), amat%ir(i + 1) - 1
            if (amat%jc(j) == amat%jc(j - 1)) then
               bmat%a(bmat_nnz - 1) = bmat%a(bmat_nnz - 1) + amat%a(j)
            else
               bmat%a(bmat_nnz) = amat%a(j)
               bmat%jc(bmat_nnz) = amat%jc(j)
               bmat%nnz = bmat_nnz
               bmat_nnz = bmat_nnz + 1
            end if
         end do
         bmat%ir(i + 1) = bmat_nnz
      end do

      deallocate (amat%a, amat%jc, amat%ir)
      allocate (amat%a(bmat%nnz), amat%jc(bmat%nnz), amat%ir(amat%nrow + 1))
      amat%a = bmat%a
      amat%jc = bmat%jc
      amat%ir = bmat%ir
      amat%nnz = bmat%nnz
      call free(bmat)

   End subroutine zcsr_add_dupes

   Subroutine blockdiag(amat, bmat, cmat)
      Implicit None
      Type(zcsrmat) :: amat, bmat
      Type(zcsrmat) :: cmat
      integer :: i, j
      !!$ Local

      call alloc(cmat, amat%nnz + bmat%nnz, amat%nrow + bmat%nrow, amat%ncol + bmat%ncol)

      cmat%nnz = amat%nnz + bmat%nnz

      cmat%a = (/amat%a, bmat%a/)
      cmat%ir = (/amat%ir, bmat%ir(2:) + amat%nnz/)
      cmat%jc = (/amat%jc, bmat%jc + amat%ncol/)

   End Subroutine blockdiag

   Function get3MMM(A, B, c) Result(D)
      Implicit None
      Type(zcsrmat) :: A, B, c, D
      Type(zcsrmat) :: t1
      t1 = spmatmul(A, B)
      D = spmatmul(t1, c)
      Call free(t1)
   End Function

      !!$ Inserts a dense matrix into a zcsr matrix on the diagonal
      !!$ at the position given by the column index of the last nonzero element
      !!$ plus one. i.e. if the last nonzero element is at (i,j), the dense mat
      !!$ will be inserted at (j+1, j+1)
   Subroutine insert_densemat_into_zcsr(amat, bmat)
      Implicit None
      Complex(Kind=DEF_DBL_PREC) :: amat(:, :)
      Type(zcsrmat) :: bmat
         !!$ Local
      Integer :: i, j, k, pos, pos0
      Integer, pointer ::  irow(:)
      irow => getrowind(bmat)

      If (bmat%alloc == 0) Then
         stop 'we require the zcsrmat to be preallocated in order to insert densemats into it'
      End if

      k = bmat%nnz + 1 !!$ index of the to-be-added nonzero element
      If (bmat%nnz == 0) Then
         pos = 1
      Else
         pos = irow(bmat%nnz) + 1
      End if
      pos0 = pos - 1
      Do i = 1, size(amat, 1)
         bmat%ir(pos+1) = bmat%ir(pos)
         Do j = 1, size(amat, 2)
            if (ABS(amat(i, j)) > 1D-20) then
               bmat%a(k) = amat(i, j)
               bmat%jc(k) = pos0 + j
               bmat%ir(pos+1) = bmat%ir(pos+1) + 1
               bmat%nnz = bmat%nnz + 1
               k = k + 1
            end if
         End do
         pos = pos + 1
      End do

   End Subroutine

End Module sparselib

