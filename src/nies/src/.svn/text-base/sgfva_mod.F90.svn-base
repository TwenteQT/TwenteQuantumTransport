Module sgfva_mod
      Use definitions
Contains

!*******************
!XXX    SGFVA   ****
!*******************
      Subroutine sgfva (opt, system, nkgrid, zcn, iiter)
!***********************************************
!!$    UPDATE OF SGF FOR SEMiiNFINITE VACUUM -
!!$     -  FOR EACH K||-POINT AND ENERGY NODE
!***********************************************
         Use green_func
         Use str_const
         Use inversion ! matrix inversion from LAPACK
         Use matmul_mod ! matrix multiplication from LAPACK
	 Use green_mod
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Integer :: ie, il, ista, ifin, ibz, ii0, i, ib, iiter, nkgrid
         Complex (Kind=prec) :: ZVPF (opt%nlsq), zpsi1 (opt%nblsq, opt%nblsq), zskf1 (opt%nblsq, opt%nblsq), &
        & zsko1 (opt%nblsq, opt%nblsq), zwa (opt%nblsq, opt%nblsq), zwb (opt%nblsq, opt%nblsq), zz, zpfx, &
        & zlax, zmux, ze, zcn (opt%ne)
         Real (Kind=prec), Parameter :: rcz = 0.0d0
         Real (Kind=prec) :: ppcx, ppdx, ppqx, alfx

         zz = dcmplx (rcz, rcz)

         If (iiter == 1) zvgam = zz

!!$  !$omp parallel do default(private)&
!!$  !$omp& shared(opt,system,zcn,nkgrid,zz,zbgam,zvgam,zbski,zbskf,zbsko,zvski,zvskf,zvsko,zskf,zsko,iiter)
! ! !SMP$ DO SERIAL
         Do ie = 1, opt%ne !--------------------- ENERGY LOOP
            ze = zcn (ie)
! ! !SMP$ DO SERIAL
            Do il = 1, opt%nl
               ppcx = system%vat(1)%ppc(il, 1)
               ppdx = system%vat(1)%ppd(il, 1)
               ppqx = system%vat(1)%ppq(il, 1)
               alfx = system%qscr (il)
               ista = (il-1) ** 2 + 1
               ifin = il ** 2

               Call plmz (ppcx, ppdx, ppqx, alfx, ze, zpfx, zlax, zmux)

!SMP$ DO SERIAL
               Do i = ista, ifin
                  ZVPF (i) = zpfx ! Pot. function
               End Do
            End Do

            Do ibz = 1, nkgrid !-----------------------ibz-LOOP
               zpsi1 (:, :) = - zvski (:, :, ibz)!       -S  structure const
               Do ib = 1, system%nb
                  ii0 = (ib-1) * opt%nlsq
!SMP$ DO SERIAL
                  Do i = 1, opt%nlsq
                     zpsi1 (ii0+i, ii0+i) = zpsi1 (ii0+i, ii0+i) + ZVPF (i)! P-S
                  End Do
               End Do
               zsko1 (:, :) = zvsko (:, :, ibz)
               zskf1 (:, :) = zvskf (:, :, ibz)
               zwa (:, :) = zvgam (:, :, ibz, ie, 1)

               zwb = mmatmul_c (zwa, zsko1)
               zwa = mmatmul_c (zskf1, zwb)

               zwb = zpsi1 - zwa

               Call lcinv (zwb)
               
!--------------STORAGE OF SGF
               zvgam (:, :, ibz, ie, 1) = zwb (:, :)! vacuum Green's function
            End Do
         End Do
!!$  !$omp end parallel do

         Return
      End Subroutine sgfva

End module
