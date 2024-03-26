!!$ $Id: postprocess.F90 1664 2013-01-09 19:48:04Z antst $
#include "math_def.h"
Module postprocess
      Use sparselib
      Use geometry_module
      Use structure
      Implicit None
      Type torqsigmas
         Integer :: alloc = 0
         Integer :: nl = 0, nr = 0, nm = 0
         Real (Kind=DEF_DBL_PREC), Allocatable :: Lir (:, :), Lt (:, :), Rir (:, :), Rt (:, :)
         Real (Kind=DEF_DBL_PREC), Allocatable :: Lm (:, :), Rm (:, :)
      End Type torqsigmas
      

Contains

      Subroutine init_torque (ts)
         Implicit None
         Type (torqsigmas), Intent (Inout) :: ts
         ts%Lm (:, :) = 0.0d0
         ts%Rm (:, :) = 0.0d0
         ts%Rir (:, :) = 0.0d0
         ts%Rt (:, :) = 0.0d0
         ts%Lir (:, :) = 0.0d0
         ts%Lt (:, :) = 0.0d0
      End Subroutine init_torque

      Subroutine alloc_torque (ts, mna, lna, rna)
         Implicit None
         Integer, Intent (In) :: mna, lna, rna
         Type (torqsigmas), Intent (Inout) :: ts

         If (ts%alloc /= 0 .And. (ts%nm /= mna .Or. ts%nl /= lna .Or. ts%nr /= rna)) Then
             call free_torque(ts)
         End If

         If (ts%alloc == 0) Then
            Allocate (ts%Lm(7, mna))
            Allocate (ts%Rm(7, mna))


            Allocate (ts%Lir(7, lna))
            Allocate (ts%Lt(7, rna))

            Allocate (ts%Rir(7, rna))
            Allocate (ts%Rt(7, lna))
            ts%alloc = 1
            ts%nl = lna
            ts%nr = rna
            ts%nm = mna
         End If

         Call init_torque (ts)
      End Subroutine alloc_torque


      Subroutine free_torque (ts)
         Implicit None
         Type (torqsigmas), Intent (Inout) :: ts
             if (ts%alloc==0) return
             Deallocate (ts%Rm, ts%Lm, ts%Lir, ts%Lt, ts%Rir, ts%Rt)
             ts%nm = 0
             ts%nl = 0
             ts%nr = 0
             ts%alloc = 0
      End Subroutine free_torque

      Subroutine bzsum_torq (ts, ts1, bzw)
         Implicit None
         Type (torqsigmas), Intent (Inout) :: ts, ts1
         Real (Kind=DEF_DBL_PREC) :: bzw

         ts%Lm (:, :) = ts%Lm(:, :) + ts1%Lm(:, :) * bzw
         ts%Rm (:, :) = ts%Rm(:, :) + ts1%Rm(:, :) * bzw

         ts%Rir (:, :) = ts%Rir(:, :) + ts1%Rir(:, :) * bzw
         ts%Rt (:, :) = ts%Rt(:, :) + ts1%Rt(:, :) * bzw

         ts%Lir (:, :) = ts%Lir(:, :) + ts1%Lir(:, :) * bzw
         ts%Lt (:, :) = ts%Lt(:, :) + ts1%Lt(:, :) * bzw
      End Subroutine bzsum_torq


      Subroutine calc_torque_EMTO (rhs, trgeo, rotm, ts, nl, nr)
         Implicit None
         Type (t_geometry_EMTO), Intent (In) :: trgeo
         Type (zdensemat), target, Intent (In) :: rhs
         Real (Kind=DEF_DBL_PREC), Intent (In) :: rotm (:, :)
         Integer, Intent (In) :: nl, nr
         Type (torqsigmas), Intent (Inout) :: ts
!!$      Local vars
         Integer :: imd, inc
         Complex (Kind=DEF_DBL_PREC), Pointer :: phi1 (:)
         Integer :: ofs1o, ofs2o, ofs1a, ofs2a
!!$ **************************************************
!!$ The wave functions are already flux-normalized
!!$ but there is a 2*pi factor to be considered later
!!$ **************************************************

         Call init_torque (ts)
         inc = 1
         ofs1o = trgeo%ltno * 2 - trgeo%num_emb_l 
         ofs2o = (trgeo%tno+trgeo%ltno) * 2 - trgeo%num_emb_l - trgeo%num_emb_m
         ofs1a = trgeo%ltna
         ofs2a = trgeo%ltna + trgeo%tna
!!$          rofs = trgeo%ltna * inc + 1
         Do imd = 1, nl + nr
            phi1 => rhs%bl (:, imd)
            If (imd > nl) Then
               Call st_helper (phi1, ts%Rt(:, :), trgeo, rotm, 0, 0, trgeo%ltna)
               Call st_helper (phi1, ts%Rm(:, :), trgeo, rotm, ofs1o, ofs1a, trgeo%tna)
               Call st_helper (phi1, ts%Rir(:, :), trgeo, rotm, ofs2o, ofs2a, trgeo%rtna)
            Else
               Call st_helper (phi1, ts%Lir(:, :), trgeo, rotm, 0, 0, trgeo%ltna)
               Call st_helper (phi1, ts%Lm(:, :), trgeo, rotm, ofs1o, ofs1a, trgeo%tna)
               Call st_helper (phi1, ts%Lt(:, :), trgeo, rotm, ofs2o, ofs2a, trgeo%rtna)
            End If

         End Do

      Contains
         Subroutine st_helper (phi, ts, trgeo, rotm, ofso, ofsa, tna)
            Implicit None
            Real (Kind=DEF_DBL_PREC), Intent (In) :: rotm (:, :)
            Real (Kind=DEF_DBL_PREC), Intent (Inout) :: ts (:, :)
            Type (t_geometry_EMTO), Intent (In) :: trgeo
            Complex (Kind=DEF_DBL_PREC) :: phi (:)
            Integer, Intent (In) :: ofsa, ofso, tna
!!$ Local vars
            Real (Kind=DEF_DBL_PREC) :: s (7), xx, yy, zz, s1, s2, s3
            Integer :: ia, ofs, no, il, i1
            Type (tclass), Pointer :: at
            Integer :: numls(5) = [1,3,5,7,9]

            ofs = ofso
            Do ia = 1, tna

               xx = Sin (-rotm(1, ia+ofsa)) * Cos (-rotm(2, ia+ofsa))
               yy = Sin (-rotm(1, ia+ofsa)) * Sin (-rotm(2, ia+ofsa))
               zz = Cos (-rotm(1, ia+ofsa))
               at => trgeo%atoms(ia+ofsa)%ptr
               no = (at%lmx+1) ** 2 - dot_product(at%idxdn-1,numls(1:at%lmx+1))

               s (1:7) = 0.0d0
               Do il = 1, at%lmx+1
                  If (at%idxdn(il-1)==2) Cycle !skip this if downfolded 
                  Do i1 = 1, 2 * (il-1) + 1
                     s1 = real (Conjg(phi(ofs+i1))*phi(ofs+no+i1))
                     s2 = Imag (Conjg(phi(ofs+i1))*phi(ofs+no+i1))
                     s3 = 0.5d0 * (Abs(phi(ofs+i1))**2-Abs(phi(ofs+no+i1))**2)
                     s (7) = s (7) + 0.5d0 * (Abs(phi(ofs+i1))**2+Abs(phi(ofs+no+i1))**2)
                     s (4) = s (4) + s1
                     s (5) = s (5) + s2
                     s (6) = s (6) + s3

                     s (1) = s (1) + s1 * at%split(il-1)
                     s (2) = s (2) + s2 * at%split(il-1)
                     s (3) = s (3) + s3 * at%split(il-1)
                  End Do
                  ofs = ofs + 2 * (il-1) + 1
               End Do

               ts (4:7, ia) = ts (4:7, ia) + s (4:7)
               ts (1, ia) = ts (1, ia) + (yy*s(3)-zz*s(2))
               ts (2, ia) = ts (2, ia) + (zz*s(1)-xx*s(3))
               ts (3, ia) = ts (3, ia) + (xx*s(2)-yy*s(1))
               ofs = ofs + no
            End Do
         End Subroutine st_helper

      End Subroutine calc_torque_EMTO


      Subroutine calc_torque (rhs, trgeo, rotm, ts, nl, nr)
         Implicit None
         Type (t_geometry), Intent (In) :: trgeo
         Type (zdensemat), target, Intent (In) :: rhs
         Real (Kind=DEF_DBL_PREC), Intent (In) :: rotm (:, :)
         Integer, Intent (In) :: nl, nr
         Type (torqsigmas), Intent (Inout) :: ts
!!$      Local vars
         Integer :: imd, inc
         Complex (Kind=DEF_DBL_PREC), Pointer :: phi1 (:)
         Integer :: ofs1o, ofs2o, ofs1a, ofs2a
!!$ **************************************************
!!$ The wave functions are already flux-normalized
!!$ but there is a 2*pi factor to be considered later
!!$ **************************************************

         Call init_torque (ts)
         inc = 1
         ofs1o = trgeo%ltno * 2
         ofs2o = (trgeo%tno+trgeo%ltno) * 2
         ofs1a = trgeo%ltna
         ofs2a = trgeo%ltna + trgeo%tna
!!$          rofs = trgeo%ltna * inc + 1
         Do imd = 1, nl + nr
            phi1 => rhs%bl (:, imd)
            If (imd > nl) Then
               Call st_helper (phi1, ts%Rt(:, :), trgeo, rotm, 0, 0, trgeo%ltna)
               Call st_helper (phi1, ts%Rm(:, :), trgeo, rotm, ofs1o, ofs1a, trgeo%tna)
               Call st_helper (phi1, ts%Rir(:, :), trgeo, rotm, ofs2o, ofs2a, trgeo%rtna)
            Else
               Call st_helper (phi1, ts%Lir(:, :), trgeo, rotm, 0, 0, trgeo%ltna)
               Call st_helper (phi1, ts%Lm(:, :), trgeo, rotm, ofs1o, ofs1a, trgeo%tna)
               Call st_helper (phi1, ts%Lt(:, :), trgeo, rotm, ofs2o, ofs2a, trgeo%rtna)
            End If

         End Do

      Contains
         Subroutine st_helper (phi, ts, trgeo, rotm, ofso, ofsa, tna)
            Implicit None
            Real (Kind=DEF_DBL_PREC), Intent (In) :: rotm (:, :)
            Real (Kind=DEF_DBL_PREC), Intent (Inout) :: ts (:, :)
            Type (t_geometry), Intent (In) :: trgeo
            Complex (Kind=DEF_DBL_PREC) :: phi (:)
            Integer, Intent (In) :: ofsa, ofso, tna
!!$ Local vars
            Real (Kind=DEF_DBL_PREC) :: s (7), xx, yy, zz, s1, s2, s3
            Integer :: ia, ofs, no, il, i1
            Type (t_atom_defenition), Pointer :: at
            ofs = ofso
            Do ia = 1, tna

               xx = Sin (-rotm(1, ia+ofsa)) * Cos (-rotm(2, ia+ofsa))
               yy = Sin (-rotm(1, ia+ofsa)) * Sin (-rotm(2, ia+ofsa))
               zz = Cos (-rotm(1, ia+ofsa))
               at => trgeo%atoms(ia+ofsa)%ptr
               no = at%nl ** 2

               s (1:7) = 0.0d0
               Do il = 1, at%nl
                  Do i1 = 1, 2 * (il-1) + 1
                     s1 = real (Conjg(phi(ofs+i1))*phi(ofs+no+i1))
                     s2 = Imag (Conjg(phi(ofs+i1))*phi(ofs+no+i1))
                     s3 = 0.5d0 * (Abs(phi(ofs+i1))**2-Abs(phi(ofs+no+i1))**2)
                     s (7) = s (7) + 0.5d0 * (Abs(phi(ofs+i1))**2+Abs(phi(ofs+no+i1))**2)
                     s (4) = s (4) + s1
                     s (5) = s (5) + s2
                     s (6) = s (6) + s3

                     s (1) = s (1) + s1 * at%split(il)
                     s (2) = s (2) + s2 * at%split(il)
                     s (3) = s (3) + s3 * at%split(il)
                  End Do
                  ofs = ofs + 2 * (il-1) + 1
               End Do

               ts (4:7, ia) = ts (4:7, ia) + s (4:7)
               ts (1, ia) = ts (1, ia) + (yy*s(3)-zz*s(2))
               ts (2, ia) = ts (2, ia) + (zz*s(1)-xx*s(3))
               ts (3, ia) = ts (3, ia) + (xx*s(2)-yy*s(1))
               ofs = ofs + no
            End Do
         End Subroutine st_helper

      End Subroutine calc_torque

      Subroutine add_zcurrent (rhs, ham, trgeo, zcurr, wgt, nl, nr)
         Implicit None
         Type (t_geometry), Intent (In) :: trgeo
         Type (zdensemat), Intent (In) :: rhs
         Real (Kind=DEF_DBL_PREC), Intent (In) :: wgt
         Type (zcsrmat), Intent (In) :: ham
         Real (Kind=DEF_DBL_PREC), Intent (Inout) :: zcurr (:, :, :)
         Integer, Intent (In) :: nl, nr
!!$           Local vars
         Real (Kind=DEF_DBL_PREC), Allocatable :: zc (:, :, :)
         Allocate (zc(size(zcurr, 1), size(zcurr, 2), size(zcurr, 3)))

         Call calc_zcurrent (rhs, ham, trgeo, zc, nl, nr)

         zcurr = zcurr + zc * wgt
      End Subroutine add_zcurrent

      Subroutine calc_zcurrent (rhs, ham, trgeo, zcurr, nl, nr)
         Implicit None
         Type (t_geometry), Intent (In) :: trgeo
         Type (zdensemat), Intent (In) :: rhs
         Type (zcsrmat), Intent (In) :: ham
         Real (Kind=DEF_DBL_PREC), Intent (Inout) :: zcurr (:, :, :)
         Integer, Intent (In) :: nl, nr
!!$           Local vars
         Complex (Kind=DEF_DBL_PREC), Allocatable :: wfr (:, :), fc (:, :)
         Integer :: i1, ofs, no, imd, ld, ia
!!$          Type (t_atom_defenition), Pointer :: at
         Real (Kind=DEF_DBL_PREC) :: jzu, jzd
         Allocate (wfr(rhs%nrow, rhs%ncol), fc(rhs%nrow, rhs%ncol))

         zcurr (:, :, :) = 0.0d0
         ofs = 0
         Do i1 = 1, trgeo%num
            no = trgeo%atoms(i1)%ptr%nl ** 2
            wfr (ofs+1:ofs+2*no, :) = rhs%bl(ofs+1:ofs+2*no, :) * trgeo%atoms(i1)%coord(3)
            ofs = ofs + no * 2
         End Do
         call spmatmul2 (ham, wfr,fc)
         Deallocate (wfr)
         Do imd = 1, nl + nr
            ld = 1
            If (imd > nl) ld = 2
            ofs = trgeo%ltno * 2
            Do ia = 1, trgeo%tna
               no = trgeo%atoms(ia+trgeo%ltna)%ptr%nl ** 2
               jzu = sum (Imag(fc(ofs+1:ofs+no, imd)*Conjg(rhs%bl(ofs+1:ofs+no, imd))))
               jzd = sum (Imag(fc(ofs+no+1:ofs+2*no, imd)*Conjg(rhs%bl(ofs+no+1:ofs+2*no, imd))))
               ofs = ofs + 2 * no
               zcurr (ia, 1, ld) = zcurr (ia, 1, ld) + jzu
               zcurr (ia, 2, ld) = zcurr (ia, 2, ld) + jzd
            End Do
         End Do
         Deallocate (fc)
      End Subroutine calc_zcurrent

      Function geom_get_cm (geo) Result (coord)
         Implicit None
         Type (t_geometry), Intent (In) :: geo
         Real (Kind=DEF_DBL_PREC), Pointer :: coord (:, :)
!!$  Local vars
         Integer :: i1, j, k
         Allocate (coord(geo%norbit*2, 3))
         k = 0
         Do i1 = 1, geo%num
            Do j = 1, geo%atoms(i1)%ptr%nm * 2
               k = k + 1
               coord (k, :) = geo%atoms(i1)%coord(1:3)
            End Do
         End Do
         If (k /= geo%norbit*2) Then
            Write (*,*) 'WTF?'
            Stop
         End If
      End Function geom_get_cm


End
