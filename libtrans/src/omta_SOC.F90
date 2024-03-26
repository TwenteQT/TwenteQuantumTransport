!!$ Module for SOC in EMTO formalism. Used in bandomtafull.x and transomtafull.x to add the SOC Hamiltonian to the complete Hamiltonian.

#include "math_def.h"

Module omta_SOC
   use omta_defs
   use omta_strrs
   use geometry_module
   use sparselib
   use logging
   implicit none

   Type t_Lop_mats
      Complex(Kind=prec), Pointer :: Lz(:, :), Ly(:, :), Lx(:, :)
      Complex(Kind=prec), Pointer :: Lp(:, :), Lm(:, :)
   End Type t_Lop_mats

   Type t_LS_mats
      Complex(Kind=prec), Pointer :: LSmat(:, :)
   End Type t_LS_mats

Contains

   Subroutine prep_so_EMTO(geom, hsoc) !!$, bsc)
  !!$ This is the subroutine which outputs the full SOC Hamiltonian matrix.
      Implicit None
      Type(t_geometry_EMTO) :: geom
      Type(zcsrmat) :: hsoc, hsoc2
      type(zcsrmat) :: temp, temp2
  !!$ Type (t_omta_struct) :: bsc(2)
  !!$ local
      Type(t_Lop_mats) :: Lop(geom%nlmax)
      Type(t_LS_mats) :: LS(2, 2, geom%nlmax)
      Type(t_LS_mats) :: Malpha(2, 2, geom%nlmax)
      Type(t_LS_mats) :: alphaMalpha(2, 2, geom%nlmax)
      Type(tclass), Pointer :: at
      Integer :: i, is, shift, j, no, il
      Integer :: ilm, jlm

      Type(zcoomat) :: sA1
      Type(zcoomat) :: aip

      j = 0
      !SMP$ DO SERIAL
      Do i = 1, geom%num
         j = j + 16 - 36*(geom%atoms(i)%ptr%lmx+1) + 20*(geom%atoms(i)%ptr%lmx+1)**2
      End Do
      j = j*4
      no = geom%norbit

      Call alloc(sA1, j, 2*no)


      Call make_Lmats(geom%nlmax, Lop)
      Call mk_LS(geom%nlmax, Lop, LS)

      shift = 0

      Do i = 1, geom%num
         at => geom%atoms(i)%ptr

         call mk_so_EMTO(aip, LS, geom%nlmax, at%xi0, at%xi1, at%xi2, i)
  !!$
         call cooaddelms(sA1, aip%A, aip%icol + shift, aip%irow + shift, aip%nnz)

         shift = shift + 2*geom%nlmax**2

      End Do

      hsoc = conv2csr(sA1)
      call ordercsr(hsoc)
      !call do_log(1,'converted sA1')

      call free(sA1)

      Call free_Lmats(Lop)
      call free_LS(LS)

      return
   End Subroutine prep_so_EMTO

   Subroutine mk_so_EMTO(aip, LS, lmax, xi0, xi1, xi2, iatom)
      Implicit None
      Type(zcoomat) :: aip
      Type(t_LS_mats), Intent(In) :: LS(2, 2, lmax)
      Real(Kind=prec), Intent(In), Pointer :: xi0(:, :, :), xi1(:, :, :), xi2(:, :, :)
      Integer, Intent(In) :: lmax, iatom
  !!$ type (t_omta_struct), Intent (In) :: bsc
  !!$ Local vars
      Integer :: no, il, is1, is2, s1, s2, ilm, jlm, il2
      Complex(Kind=prec), Pointer :: DM(:, :)
      Real(Kind=prec), pointer :: Slocal(:, :)

      call free(aip)

      Allocate (DM(2*(lmax)**2, 2*(lmax)**2))
  !!$ Allocate (Slocal(bsc%iax(7,iatom,1),bsc%iax(6,iatom,1)))

  !!$ do jlm = 1, bsc%iax(6,iatom,1)
  !!$   do ilm = 1, bsc%iax(7,iatom,1)

  !!$     Slocal(ilm,jlm) = bsc%s(bsc%iax(9,iatom,1)+(jlm-1)*bsc%iax(7,iatom,1)+ilm)

  !!$   end do
  !!$ end do

      no = (lmax)**2
      DM = 0.0d0

      do il = 2, lmax
         do is1 = 1, 2
            s1 = (is1 - 1)*no
            do is2 = 1, 2
               s2 = (is2 - 1)*no
  !!$ the xi0 terms, (onsite - onsite)
               DM(((il - 1)**2 + 1 + s1):(il**2 + s1), ((il - 1)**2 + 1 + s2):(il**2 + s2)) =  &
                       & xi0(is1, is2, il - 1)*LS(is1, is2, il)%Lsmat(1:2*il - 1, 1:2*il - 1)

  !!$ now the xi1 terms (onsite - offsite)
               !                DM (((il-1)**2+1+s1) :(il**2+s1), ((il-1)**2+1+s2) :(il**2+s2)) = DM (((il-1)**2+1+s1) :(il**2+s1), ((il-1)**2+1+s2) :(il**2+s2)) - xi1(is1, is2, il-1) * matmul ( LS(is1, is2, il)%Lsmat ,  Slocal((il-1)**2+1 : il**2, (il-1)**2+1 : il**2)) - xi1(is1, is2, il-1) * matmul ( Slocal((il-1)**2+1 : il**2, (il-1)**2+1 : il**2) , LS(is1, is2, il)%Lsmat )
  !!$ finally, the xi2 terms (offsite - offsite)
               !                do il2 = 1, lmax
               !                        DM (((il-1)**2+1+s1) :(il**2+s1), ((il-1)**2+1+s2) :(il**2+s2)) = DM (((il-1)**2+1+s1) :(il**2+s1), ((il-1)**2+1+s2) :(il**2+s2)) + xi2(is1, is2, il2-1) * matmul(matmul(transpose(Slocal((il2-1)**2+1 : il2**2, (il-1)**2+1 : il**2)), LS(is1, is2, il2)%Lsmat) ,  Slocal((il2-1)**2+1 : il2**2, (il-1)**2+1 : il**2)  )
               !                end do

            end do
         end do
      end do

      call zdense2coo(DM, aip, 1.0d-20)

      Deallocate (DM)
  !!$ Deallocate (Slocal)

      return
   End Subroutine mk_so_EMTO

   Subroutine mk_LS(nl, Lop, LS)
!!$ Here we generate the matrix elements of L.S for real spherical harmonics.
      implicit none
      Type(t_LS_mats), Intent(Out) :: LS(:, :, :)
      Type(t_Lop_mats), Intent(In) :: Lop(:)
      Integer :: nl
!!$ local vars
      Complex(Kind=prec), Parameter :: Sz(2, 2) = reshape((/0.5d0, 0.0d0, 0.0d0, -0.5d0/), (/ &
     & 2, 2/))
      Complex(Kind=prec), Parameter :: Sy(2, 2) = DEF_cmplx_Ione*reshape((/0.0d0, &
     & .5d0, -.5d0, 0.0d0/), (/2, 2/))
      Complex(Kind=prec), Parameter :: Sx(2, 2) = reshape((/0.0d0, .5d0, .5d0, 0.0d0/), (/ &
     & 2, 2/))
      Integer :: il, is1, is2
      Integer :: i1, i2

      Do il = 1, nl
         Do is1 = 1, 2
            Do is2 = 1, 2
               Allocate (LS(is1, is2, il)%Lsmat(1:2*il - 1, 1:2*il - 1))
               LS(is1, is2, il)%Lsmat(:, :) = Sx(is1, is2)*Lop(il)%Lx + Sy(is1, is2)*Lop(il)%Ly + Sz(is1, is2)*Lop(il)%Lz
            End do
         End do
      End do

! open(unit=111,file='Lsmat')
! Do il = 1, nl
!    Do is1 = 1, 2
!        Do is2 = 1, 2
!        do i1 = 1,2*il-1
!        do i2 = 1,2*il-1
!            write(111,"(I2, I2, I2, I2, I2, F8.4, F8.4)") il, is1, is2, i1, i2, real(LS(is1, is2, il)%Lsmat(i1, i2)), aimag(LS(is1, is2, il)%Lsmat(i1, i2))
!        end do
!        end do
!        End do
!    End do
! End do
! close(unit=111)

      return
   End Subroutine mk_LS

   Subroutine free_LS(LS)
      Implicit none
      Type(t_LS_mats) :: LS(:, :, :)
!!$ Local
      Integer :: is1, is2, il
      Do is1 = 1, 2
         Do is2 = 1, 2
            Do il = 1, size(LS(is1, is2, :))
               Deallocate (LS(is1, is2, il)%Lsmat)
               Nullify (LS(is1, is2, il)%Lsmat)
            End Do
         End Do
      End Do
   End Subroutine free_LS

   Subroutine make_Lmats(lmax, Lop)
      Implicit None
      Type(t_Lop_mats) :: Lop(lmax)
      Integer :: lmax
!!$ Local vars
      Complex(Kind=prec), Parameter :: ci = DEF_cmplx_Ione
      Integer :: i
      Do i = 1, lmax
         Allocate (Lop(i)%Lz(2*i - 1, 2*i - 1))
         Allocate (Lop(i)%Lm(2*i - 1, 2*i - 1))
         Allocate (Lop(i)%Lp(2*i - 1, 2*i - 1))
         Allocate (Lop(i)%Lx(2*i - 1, 2*i - 1))
         Allocate (Lop(i)%Ly(2*i - 1, 2*i - 1))
         Call CalcLInts(Lop(i)%Lp, Lop(i)%Lm, Lop(i)%Lz, i)
         Lop(i)%Lx = .5d0*(Lop(i)%Lp + Lop(i)%Lm)
         Lop(i)%Ly = -ci*.5d0*(Lop(i)%Lp - Lop(i)%Lm)
      End Do
   End Subroutine make_Lmats

   Subroutine make_bigLmat(lmax,Lop)
     Implicit None
     Type (t_Lop_mats) :: Lop(lmax)
     Integer :: lmax, i, j
     Complex (Kind=prec), Allocatable :: Lp(:,:), Lm(:,:)
     Complex(Kind=prec), Parameter :: ci = DEF_cmplx_Ione
     Integer :: ilm, lm1, lm2

     do ilm = 1,lmax

      Allocate(Lop(ilm)%Lz(2*ilm**2,2*ilm**2), Lop(ilm)%Ly(2*ilm**2,2*ilm**2), Lop(ilm)%Lx(2*ilm**2,2*ilm**2))
      Allocate(Lop(ilm)%Lp(1,1), Lop(ilm)%Lm(1,1)) !needed to free the struct without segfault
      Lop(ilm)%Lz = 0d0
      Lop(ilm)%Lx = 0d0
      Lop(ilm)%Ly = 0d0
      Allocate(Lm(ilm**2, ilm**2), Lp(ilm**2, ilm**2))
      do i=1, ilm
       j = 2*i-1
       lm1 = i**2-j+1
       lm2 = i**2

       call CalcLInts(Lp(1:j, 1:j), Lm(1:j, 1:j), Lop(ilm)%Lz(lm1:lm2, lm1:lm2), i)
       Lop(ilm)%Lx(lm1:lm2, lm1:lm2) =     .5d0*(Lp(1:j, 1:j) + Lm(1:j, 1:j))
       Lop(ilm)%Ly(lm1:lm2, lm1:lm2) = -ci*.5d0*(Lp(1:j, 1:j) - Lm(1:j, 1:j))

       Lop(ilm)%Lz(lm1+ilm**2:lm2+ilm**2, lm1+ilm**2:lm2+ilm**2) = Lop(ilm)%Lz(lm1:lm2, lm1:lm2)
       Lop(ilm)%Lx(lm1+ilm**2:lm2+ilm**2, lm1+ilm**2:lm2+ilm**2) = Lop(ilm)%Lx(lm1:lm2, lm1:lm2)
       Lop(ilm)%Ly(lm1+ilm**2:lm2+ilm**2, lm1+ilm**2:lm2+ilm**2) = Lop(ilm)%Ly(lm1:lm2, lm1:lm2)
      enddo

      Deallocate(Lm, Lp)

     enddo

   End Subroutine make_bigLmat

   Subroutine free_Lmats(Lop)
      Type(t_Lop_mats) :: Lop(:)
!!$ Local
      Integer :: i
      Do i = 1, size(Lop)
         Deallocate (Lop(i)%Lz, Lop(i)%Ly, Lop(i)%Lx, Lop(i)%Lm, Lop(i)%Lp)
         Nullify (Lop(i)%Lz, Lop(i)%Lp, Lop(i)%Lm, Lop(i)%Ly, Lop(i)%Lx)
      End Do
   End Subroutine free_Lmats

   Subroutine CalcLInts(Lp, Lm, Lz, il)
      Implicit None
      Integer :: il
      Complex(Kind=prec) :: Lp(:, :), Lm(:, :), Lz(:, :)
!!$ Local vars
      Complex(Kind=prec) :: R(-(il - 1):(il - 1), -(il - 1):(il - 1))
      Complex(Kind=prec) :: Ri(-(il - 1):(il - 1), -(il - 1):(il - 1))
      Complex(Kind=prec) :: Lb(-(il - 1):(il - 1), -(il - 1):(il - 1))
      Integer :: l, M
      l = il - 1

      Call transmat(R, l)
      Ri = transpose(conjg(R))

      Lz(:, :) = 0.0d0
      Lp(:, :) = 0.0d0
      Lm(:, :) = 0.0d0

      If (il > 1) Then

         Lb(:, :) = DEF_cmplx_zero
         Forall (M=-l:l) Lb(M, M) = cmplx(M, 0.0d0, kind=prec)
         Lz(1:(2*l + 1), 1:(2*l + 1)) = .5d0*matmul(R, matmul(Lb, Ri))
         Lb(:, :) = DEF_cmplx_zero
         Forall (M=-l:(l - 1)) Lb(M + 1, M) = Sqrt(cmplx(l*(l + 1) - M*(M + 1), 0.0d0, kind=prec))
         Lp(1:(2*l + 1), 1:(2*l + 1)) = .5d0*matmul(R, matmul(Lb, Ri))
         Lm = transpose(conjg(Lp))
      End If
   Contains

      Subroutine transmat(R, l)
!!$ Transformation matrix complex->real harmonics (without 1/sqrt(2) prefactor!!!)
         Implicit None
         Integer :: l
         Complex(Kind=prec) :: R(-l:l, -l:l)
!!$
         Integer :: M
         R(:, :) = DEF_cmplx_zero
         R(0, 0) = DEF_cmplx_one*Sqrt(2.0d0)
         Do M = 1, l, 1
            R(-M, -M) = DEF_cmplx_Ione
            R(M, M) = (-1.0d0)**M
            R(M, -M) = DEF_cmplx_one
            R(-M, M) = -DEF_cmplx_Ione*(-1.0d0)**M
         End Do
      End Subroutine transmat
   End Subroutine CalcLInts

   Subroutine spherm_coo(coo_in, coo_out)
      implicit none
      type(zcoomat), Intent(In)  :: coo_in
      type(zcoomat), Intent(Out) :: coo_out

      call alloc(coo_out, coo_in%nnz, coo_in%nrow, coo_in%ncol)
      coo_out%a = conjg(coo_in%a)
      coo_out%irow = coo_in%icol
      coo_out%icol = coo_in%irow
      coo_out%nrow = coo_in%nrow
      coo_out%ncol = coo_in%ncol
      coo_out%nnz = coo_in%nnz
   end subroutine spherm_coo

   Subroutine restructure(zcsr_in, zcsr_out, geom)
      Implicit none
      Type(zcsrmat) :: zcsr_in, zcsr_out
      Type(t_geometry_EMTO) :: geom
        !!$ local
      Integer :: i, n, ldim, atomNumberRow, atomNumberCol
      Type(zcoomat) :: zcoo_temp, zcoo_temp2

      ldim = zcsr_in%ncol/geom%num/2

      call csr2coo(zcsr_in, zcoo_temp)

      do i = 1, zcoo_temp%nnz

         atomNumberRow = (zcoo_temp%irow(i) - 1)/ldim + 1
         atomNumberCol = (zcoo_temp%icol(i) - 1)/ldim + 1

         if (atomNumberRow > geom%num) then
            atomNumberRow = atomNumberRow - geom%num
         end if
         if (atomNumberCol > geom%num) then
            atomNumberCol = atomNumberCol - geom%num
         end if

         if (zcoo_temp%icol(i) .le. zcoo_temp%ncol/2) then
            zcoo_temp%icol(i) = zcoo_temp%icol(i) + (atomNumberCol - 1)*ldim
         else
            zcoo_temp%icol(i) = zcoo_temp%icol(i) + (atomNumberCol - geom%num)*ldim
         end if

         if (zcoo_temp%irow(i) .le. zcoo_temp%nrow/2) then
            zcoo_temp%irow(i) = zcoo_temp%irow(i) + (atomNumberRow - 1)*ldim
         else
            zcoo_temp%irow(i) = zcoo_temp%irow(i) + (atomNumberRow - geom%num)*ldim
         end if
      end do

      call alloc(zcsr_out, zcoo_temp%nnz, zcoo_temp%nrow, zcoo_temp%ncol)
      call coo2csr(zcoo_temp, zcsr_out)

      call free(zcsr_in)

      return
   End Subroutine restructure

   Subroutine csr2coo(zcsr, zcoo)
      implicit none
      Type(zcsrmat) :: zcsr
      Type(zcoomat) :: zcoo
!local
      Integer :: i, jfirst, jlast

      call alloc(zcoo, zcsr%nnz, zcsr%nrow, zcsr%ncol)

      zcoo%nnz = zcsr%nnz
      zcoo%ncol = zcoo%ncol
      zcoo%nrow = zcoo%nrow

      zcoo%a = zcsr%a
      zcoo%icol = zcsr%jc
      jlast = 0
      do i = 2, size(zcsr%ir)
         jfirst = jlast + 1
         jlast = jlast + zcsr%ir(i) - zcsr%ir(i - 1)
         zcoo%irow(jfirst:jlast) = i - 1
      end do

   end subroutine csr2coo

End Module omta_SOC
