Module green_mod
      Use definitions
Contains

!*******************
!XXX    gamma    ****
!*******************
      Subroutine gamma (opt, system, nkgrid, kweight, zcn, ncase, ns, num)

!***************************************************
!!$    GAMMA FOR THE SEMIINFINITE LEFT&RIGHT BULK RANDOM ALLOY -
!!$     -  FOR EACH K||-POINT AND ENERGY NODE
!***************************************************
         Use green_func
         Use str_const
         Use inversion ! matrix inversion from LAPACK
         Use matmul_mod ! matrix multiplication from LAPACK
         Implicit None
         Real (Kind=prec), External :: dzasum ! blas function
         Type (params) :: opt
         Type (atoms_set) :: system
         Real (Kind=prec), Parameter :: rcz = 0.0d0, rc1 = 1.0d0
         Integer :: i, is, ie, its, ia, il, ib, nonc, nons, ibz, ii0, nkgrid, ncase, ns, num, stopf
         Complex (Kind=prec) :: zcn (opt%ne)
         Complex (Kind=prec) :: zz, ze, zpfx, zla, zmu, zdum
         Complex (Kind=prec), Pointer :: ZPF (:, :), ZCP (:, :, :), zom (:, :, :), zgf (:, :, :), zfi1 (:, &
        & :), ZCP1 (:, :), zom1 (:, :)
         Complex (Kind=prec), Pointer :: zw (:, :), zskf1 (:, :), zsko1 (:, :)
         Complex (Kind=prec), Pointer :: zgmb1 (:, :), zgmb (:, :, :), zgmv1 (:, :), zgmv (:, :, :)
         Complex (Kind=prec), Pointer :: zwc (:, :), zwa (:, :), zwb (:, :), ztm (:, :)
         Real (Kind=prec) :: ppc, ppd, ppq, alf, sum, sum1, kweight (nkgrid)
         Type (bulk_atom_definition), Pointer :: ab

         If (ncase == 0) write (IW6, 111)
111      Format (/ / 5 X, '**** CONVERGENCE OF LEFT BULK GAMMA:' /)
         If (ncase == 1 .And. opt%lident == 0) write (IW6, 112)
112      Format (/ / 5 X, '**** CONVERGENCE OF RIGHT BULK GAMMA:' /)
         flush (IW6)
         If (ncase /= 1 .And. ncase /= 0) Then
            Write (IW6, 116)
            Stop
         End If
116      Format (/ / 5 X, '**** ncase must be 0 or 1 :' /)

113      Format (/ / ' ****  LEFT BULK:  NON-CONVERGED ', 'AFTER ', I4, ' ITERATIONS' / '      ---  FOR IE=', &
        & I3, '    IS=', I2)
114      Format (/ / ' ****  RIGHT BULK:  NON-CONVERGED ', 'AFTER ', I4, ' ITERATIONS' / '      ---  FOR IE=',&
        &  I3, '    IS=', I2)
115      Format (/ / 5 X, '**** TRY TO INCREASE NUMBER OF ITERATIONS: opt%nmits ' /)

         zz = dcmplx (rcz, rcz)
         stopf = 0
!------------------------------- LOOP OVER SPIN AND ENERGY
!SMP$ DO SERIAL
         Do is = 1, ns
            Allocate (ZPF(opt%nlsq, system%bnum))
            Allocate (ZCP(opt%nlsq, opt%nlsq, system%nb))
            Allocate (zom(opt%nlsq, opt%nlsq, system%nb))
            Allocate (zgf(opt%nlsq, opt%nlsq, system%nb))
            Allocate (zw(opt%nlsq, opt%nlsq))
            Allocate (zgmv(opt%nblsq, opt%nblsq, nkgrid))
            Allocate (zgmb(opt%nblsq, opt%nblsq, nkgrid))
            Do ie = 1, opt%ne
               ze = zcn (ie)
               its = 0
               zgmb = zz
               zgmv = zz
210            its = its + 1 ! iteration counter
               If (its .Gt. opt%nmits) Then
                  If (ncase == 0) write (IW6, 113) opt%nmits, ie, is
                  If (ncase == 1) write (IW6, 114) opt%nmits, ie, is
                  Write (IW6, 115)
                  flush (IW6)
                  stopf = 1
               End If
!------------------------------------ CPA-PART
               nonc = system%nb
               If (its == 1) Then ! 1ST ITERATION
                  Do ia = 1, num
                     If (ncase == 0) ab => system%bat(ia)
                     If (ncase == 1) ab => system%vat(ia)
                     Do il = 1, opt%nl
                        ppc = ab%ppc (il, is)
                        ppd = ab%ppd (il, is)
                        ppq = ab%ppq (il, is)
                        alf = system%qscr (il)
                        Call plmz (ppc, ppd, ppq, alf, ze, zpfx, zla, zmu)
                        ZPF ((il-1)**2+1:il**2, ia) = zpfx
                     End Do
                  End Do
                  zom = zz
                  ZCP = zz
                  Do i = 1, opt%nlsq
                     zdum = zz
                     Do ia = 1, num
                        zdum = zdum + ab%con / ZPF (i, ia)
                     End Do
                     zdum = rc1 / zdum
                     Do ib = 1, system%nb
                        ZCP (i, i, ib) = zdum
                     End Do
                  End Do
               Else
                  Do ib = 1, system%nb
                     zfi1 => zgf (:, :, ib)
                     ZCP1 => ZCP (:, :, ib)
                     zom1 => zom (:, :, ib)

                     Call lcinv (zfi1)

                     zw = ZCP1 - zfi1 - zom1
                     sum = dzasum (opt%nlsq*opt%nlsq, zw, 1)
                     sum1 = dzasum (opt%nlsq*opt%nlsq, zom1, 1)
                     If (sum .Lt. (opt%difmc*sum1)) nonc = nonc - 1
                     zom1 = zom1 + opt%qmixc * zw

                     Call symgf (zom1, opt%bz%nsym, opt%nlsq)
                     zfi1 = zz
                     Do ia = 1, num
                        zw = - zom1
                        Do i = 1, opt%nlsq
                           zw (i, i) = zw (i, i) + ZPF (i, ia)
                        End Do
                        Call lcinv (zw)
                        zfi1 = zfi1 + ab%con * zw
                     End Do
                     Call lcinv (zfi1)
                     ZCP1 = zfi1 + zom1
                     Call symgf (ZCP1, opt%bz%nsym, opt%nlsq)
                  End Do
               End If
!---------------------------------------- SGF - PART
               zgf = zz
               nons = 0
!$omp parallel default(shared)&
!$omp& private(ibz,ztm,zwa,zwb,zskf1,zsko1,zgmb1,zgmv1,ib,i,ii0,sum,sum1) reduction(+:NONS)
               Allocate (ztm(opt%nblsq, opt%nblsq))
               Allocate (zwa(opt%nblsq, opt%nblsq), zwb(opt%nblsq, opt%nblsq))
!$omp do
               Do ibz = 1, nkgrid !-------------------------ibz-LOOP
                  If (ncase == 0) Then
                     ztm (:, :) = - zbski (:, :, ibz)!   ---  IN-LAYER QUANTITIES LEFT
                     zskf1 => zbskf (:, :, ibz)!              LEFT BULK TRANSFER MATRICES
                     zsko1 => zbsko (:, :, ibz)
                  Else
                     ztm (:, :) = - zvski (:, :, ibz)!   ---  IN-LAYER QUANTITIES RIGHT
                     zskf1 => zvskf (:, :, ibz)!              RIGHT BULK TRANSFER MATRICES
                     zsko1 => zvsko (:, :, ibz)
                  End If
                  zgmb1 => zgmb (:, :, ibz)
                  zgmv1 => zgmv (:, :, ibz)

                  Do ib = 1, system%nb
                     ii0 = (ib-1) * opt%nlsq
                     ztm (ii0+1:ii0+opt%nlsq, ii0+1:ii0+opt%nlsq) = ztm (ii0+1:ii0+opt%nlsq, &
                    & ii0+1:ii0+opt%nlsq) + ZCP (:, :, ib)
                  End Do

!!$               ITERATION OF BULK-SIDE GAMMA
                  zwa = ztm - zgmb1
                  Call lcinv (zwa)
                  zwb = mmatmul_c (zwa, zskf1)
                  zwa = mmatmul_c (zsko1, zwb)
                  zwa = zwa - zgmb1
                  sum = dzasum (opt%nblsq*opt%nblsq, zwa, 1)
                  sum1 = dzasum (opt%nblsq*opt%nblsq, zgmb1, 1)
                  If (sum .Lt. (opt%difms*sum1)) nons = nons + 1
                  zgmb1 = zgmb1 + opt%qmixs * zwa

!!$    ITERATION OF VACUUM-SIDE GAMMA
                  zwa = ztm - zgmv1
                  Call lcinv (zwa)
                  zwb = mmatmul_c (zwa, zsko1)
                  zwa = mmatmul_c (zskf1, zwb)
                  zwa = zwa - zgmv1
                  sum = dzasum (opt%nblsq*opt%nblsq, zwa, 1)
                  sum1 = dzasum (opt%nblsq*opt%nblsq, zgmv1, 1)
                  If (sum .Lt. (opt%difms*sum1)) nons = nons + 1
                  zgmv1 = zgmv1 + opt%qmixs * zwa

                  zwa = ztm - zgmv1 - zgmb1 !                                  ON-SITE GF-BLOCKS
                  Call lcinv (zwa)
!$OMP CRITICAL
                  Do ib = 1, system%nb
                     ii0 = (ib-1) * opt%nlsq
                     zgf (:, :, ib) = zgf (:, :, ib) + kweight (ibz) * zwa (ii0+1:ii0+opt%nlsq, &
                    & ii0+1:ii0+opt%nlsq)
                  End Do
!$OMP END CRITICAL
               End Do !............................... END OF ibz-LOOP
!$omp end do
               Deallocate (ztm, zwa, zwb)
!$omp end parallel
               Do ib = 1, system%nb
                  zw (:, :) = zgf (:, :, ib)
                  Call symgf (zw, opt%bz%nsym, opt%nlsq)
                  zgf (:, :, ib) = zw (:, :)
               End Do

               If ((2*nkgrid-nons+nonc) .Gt. 0) Go To 210
               Write (IW6, 150) is, ie, its
               flush (IW6)
150            Format (10 X, 'IS=', I2, 5 X, 'IE=', I3, 8 X, I4, ' ITER.')

!------- CONNECTION TO THE INTERMEDIATE REGION
!$omp parallel default(shared) private(ztm,zwb,zwc,zwa,ib,ii0,i,ibz)
               Allocate (ztm(opt%nblsq, opt%nblsq))
               Allocate (zwb(opt%nblsq, opt%nblsq))
!$omp do
               Do ibz = 1, nkgrid
                  If (ncase == 0) Then
!-----------IN-LAYER QUANTITIES FOR LEFT BULK AND RIGHT ONE IN CASE LEADS ARE IDENTICAL (lidet==0)
                     ztm (:, :) = - zbski (:, :, ibz) - zgmb (:, :, ibz)
                     zwc => zskf (:, :, system%np, ibz)
                     zwa => zbgam (:, :, ibz, ie, is)
                  Else
!-----------IN-LAYER QUANTITIES RIGHT BULK NONEQUIVALENT LEADS
                     ztm (:, :) = - zvski (:, :, ibz) - zgmv (:, :, ibz)
                     zwc => zsko (:, :, 1, ibz)
                     zwa => zvgam (:, :, ibz, ie, is)
                  End If

                  Do ib = 1, system%nb
                     ii0 = (ib-1) * opt%nlsq
                     ztm (ii0+1:ii0+opt%nlsq, ii0+1:ii0+opt%nlsq) = ztm (ii0+1:ii0+opt%nlsq, &
                    & ii0+1:ii0+opt%nlsq) + ZCP (:, :, ib)
                  End Do

                  Call lcinv (ztm)
                  Call mmatmul_z (zwb, ztm, zwc)
                  Call mmatmul_hz (zwa, zwc, zwb)
               End Do
!$omp end do
               Deallocate (ztm, zwb)
!$omp end parallel
!------- CONNECTION TO THE INTERMEDIATE REGION IN CASE LEADS ARE IDENTICAL
               If (opt%lident == 1) Then
!$omp parallel default(shared) private(ztm,zwb,zwc,zwa,ib,ii0,i,ibz)
                  Allocate (ztm(opt%nblsq, opt%nblsq))
                  Allocate (zwb(opt%nblsq, opt%nblsq))
!$omp do
                  Do ibz = 1, nkgrid
                     ztm (:, :) = - zbski (:, :, ibz) - zgmv (:, :, ibz)
                     Do ib = 1, system%nb
                        ii0 = (ib-1) * opt%nlsq
                        ztm (ii0+1:ii0+opt%nlsq, ii0+1:ii0+opt%nlsq) = ztm (ii0+1:ii0+opt%nlsq, &
                       & ii0+1:ii0+opt%nlsq) + ZCP (:, :, ib)
                     End Do
                     zwc => zsko (:, :, 1, ibz)!---------------TRANSFER MATRICES
!!$       INVERSION AND  MULTIPLICATIONS
                     Call lcinv (ztm)
                     Call mmatmul_z (zwb, ztm, zwc)
                     Call mmatmul_hz (zvgam(:, :, ibz, ie, is), zwc, zwb)
                  End Do
!$omp end do
                  Deallocate (ztm, zwb)
!$omp end parallel
               End If
            End Do !----------------------------- END OF ENERGY AND SPIN LOOP
            Deallocate (ZPF, ZCP, zom, zgf, zw, zgmb, zgmv)
         End Do
         If (stopf == 1) Stop
         If (ncase == 0) Then
            Deallocate (zbski, zbsko, zbskf)
         Else
            Deallocate (zvski, zvsko, zvskf)
         End If


         If (ncase == 0 .And. opt%lident == 1) write (IW6, 117)
117      Format (/ / '**** LEADS ARE IDENTICAL-CONVERGENCE FOR LEFT AND RIGHT BULK GAMMA IS THE SAME' /)
         flush (IW6)

         Return
      End Subroutine gamma



!*******************
!XXX    scpf    ****
!*******************
      Subroutine scpf (opt, system, zcn, iiter)
!**************************************************
!!$   SET UP THE COHERENT POTENTIAL FUNCTIONS
!**************************************************
         Use green_func
         Use str_const
         Use inversion ! matrix inversion from LAPACK
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Complex (Kind=prec) :: zp (opt%nlsq), zfi1 (opt%nlsq, opt%nlsq), ZCP1 (opt%nlsq, opt%nlsq), zom1 &
        & (opt%nlsq, opt%nlsq), zw (opt%nlsq, opt%nlsq), ze, zz, zpfx, zlax, zmux, zcn (opt%ne)
         Real (Kind=prec), Parameter :: rcz = 0.0d0, rch = 0.5d0, TOL = 1.0d-6
         Real (Kind=prec) :: ppcx, ppdx, ppqx, alfx
         Integer :: is, ie, il, ig, i, j, ista, ifin, ic, ia, iiter
         Type (atom_definition), Pointer :: aa

         zz = dcmplx (rcz, rcz)

100      Format (5 A16)
101      Format (1 X, 10 I5)
104      Format (1 X, 4 G15.7)

         If (iiter == 1 .And. .Not. (opt%rwgamma == 1 .And. opt%icont == 1)) zomg = zz !******************************** THE FIRST ITERATION

         Do is = 1, opt%ns !-------------------------------- LOOP OVER SPIN AND ENERGY
!$omp parallel do default(private) &
!$omp& shared(is,opt,system,zomg,zcn,zz,zbgam,zvgam,zbski,zbskf,zbsko,zvski,zvskf,zvsko,zskf,zsko,iiter,zcpf,IW6)
            Do ie = 1, opt%ne
               ze = zcn (ie)
               ia = 0

               Do ig = 1, system%nums !------------------------------- LOOP OVER SITES
                  zom1 (:, :) = zomg (:, :, ig, ie, is)
                  zfi1 (:, :) = zz

                  Do ic = 1, system%nc(ig)!------------------------------ LOOP OVER COMPONENTS
                     ia = ia + 1 !------------------------------ atoms counter
                     aa => system%at (ia)
                     Do il = 1, opt%nl
                        ppcx = aa%ppc (il, is)
                        ppdx = aa%ppd (il, is)
                        ppqx = aa%ppq (il, is)
                        alfx = system%qscr (il)
                        ista = (il-1) ** 2 + 1
                        ifin = il ** 2

                        Call plmz (ppcx, ppdx, ppqx, alfx, ze, zpfx, zlax, zmux)

!SMP$ DO SERIAL
                        Do i = ista, ifin
                           zp (i) = zpfx
                        End Do
                     End Do

                     zw = - zom1
!SMP$ DO SERIAL
                     Do i = 1, opt%nlsq
                        zw (i, i) = zw (i, i) + zp (i)
                     End Do

                     Call lcinv (zw)

                     Do j = 1, opt%nlsq
!SMP$ DO SERIAL
                        Do i = 1, opt%nlsq
                           zfi1 (i, j) = zfi1 (i, j) + aa%con * zw (i, j)
                        End Do
                     End Do
                  End Do !------------------------------------ END OF LOOP OVER COMPONENTS

                  Call lcinv (zfi1)

                  ZCP1 = zfi1 + zom1

                  Call symgf (ZCP1, opt%bz%nsym, opt%nlsq)

                  zcpf (:, :, ig, ie, is) = ZCP1 (:, :)! potential function

               End Do !-------------------------- end of loop over the sites
            End Do !-------------------------- over the energy
!$omp end parallel do
         End Do !-------------------------  and a spin

!!$       open(15,file='zcpf1',status='unknown')
!!$       do i=1,opt%NLSQ
!!$        do j=1,opt%NLSQ
!!$         do ig=1,system%nums
!!$          do ie=1,opt%NE
!!$           do is=1,opt%NS
!!$            write(15,*) zcpf(I,J,IG,IE,IS)
!!$           enddo
!!$          enddo
!!$         enddo
!!$        enddo
!!$       enddo
!!$       close(15)

         Return
      End Subroutine scpf

!*******************
!XXX    REPA    ****
!*******************
      Subroutine repa (opt, system, nkgrid, kweight)
!************************************************************
!!$ RECURSION PARTITIONING FOR ON-LAYER BLOCKS OF THE GF MATRIX
!************************************************************
         Use green_func
         Use str_const
         Use inversion ! matrix inversion from LAPACK
         Use matmul_mod ! matrix multiplication from LAPACK
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Integer :: isb, isv, is, ibz, ie, ip, ib, ii0, ip1, ig, nkgrid
         Complex (Kind=prec) :: zz
         Complex (Kind=prec), Pointer :: zpsi (:, :, :), zwa (:, :), zwb (:, :)
         Complex (Kind=prec), Pointer :: zgmv (:, :, :), zgmb (:, :, :), zgfg1 (:, :, :)
         Real (Kind=prec), Parameter :: rcz = 0.0d0
         Real (Kind=prec) :: kweight (nkgrid)
         zz = dcmplx (rcz, rcz)

         zgfg = zz

! ! !SMP$ DO SERIAL
         Do is = 1, opt%ns !---------------- SPIN LOOP AND ENERGY LOOP
            isb = Min (is, system%bns)
            isv = Min (is, system%vns)
!!$  !!$omp parallel do default(private)&
!!$  !!$omp& shared(kweight,is,nkgrid,opt,system,zomg,zz,ISB,ISV,IW6)&
!!$  !!$omp& shared(zbgam,zvgam,zbski,zbskf,zbsko,zvski,zvskf,zvsko,zskf,zsko,zcpf,zski,zgfg)
! ! !SMP$ DO SERIAL
            Do ie = 1, opt%ne
               zgfg1 => zgfg (:, :, :, ie, is)
!$omp parallel default(shared)&
!$omp& private(ibz,zpsi,ig,ip,ib,ii0,zgmv,ip1,zwa,zwb,zgmb)
!!$ ! !$reduction(+:zgfg1)
               Allocate (zgmv(opt%nblsq, opt%nblsq, system%np))
               Allocate (zgmb(opt%nblsq, opt%nblsq, system%np))
               Allocate (zwa(opt%nblsq, opt%nblsq))
               Allocate (zwb(opt%nblsq, opt%nblsq))
               Allocate (zpsi(opt%nblsq, opt%nblsq, system%np))
!$omp do
               Do ibz = 1, nkgrid !--------------- ibz LOOP
                  zpsi (:, :, :) = - zski (:, :, :, ibz)! -S
                  ig = 0
                  Do ip = 1, system%np
                     Do ib = 1, system%nb
                        ig = ig + 1
                        ii0 = (ib-1) * opt%nlsq
                        zpsi (ii0+1:ii0+opt%nlsq, ii0+1:ii0+opt%nlsq, ip) = zpsi (ii0+1:ii0+opt%nlsq, &
                       & ii0+1:ii0+opt%nlsq, ip) + zcpf (:, :, ig, ie, is)
                     End Do
                  End Do
!!$      RECURSION FROM VACUUM (embedding) for the first principal layer
!!$      in the intermediate region (see documentation)
                  zgmv (:, :, 1) = zvgam (:, :, ibz, ie, isv)

                  If (opt%ivac == 0) Then
                     Call mmatmul_z (zwb, zgmv(:, :, 1), zsko(:, :, 1, ibz))
                     Call mmatmul_hz (zgmv(:, :, 1), zsko(:, :, 1, ibz), zwb)
                  End If

                  If (system%np .Gt. 1) Then
                     Do ip = 2, system%np
                        ip1 = ip - 1
                        zwa (:, :) = zpsi (:, :, ip1) - zgmv (:, :, ip1)
                        Call lcinv (zwa)
                        Call mmatmul_z (zwb, zwa, zsko(:, :, ip, ibz))
                        Call mmatmul_z (zgmv(:, :, ip), zskf(:, :, ip1, ibz), zwb)
                     End Do
                  End If
!!$                              RECURSION FROM BULK (embedding)
                  zgmb (:, :, system%np) = zbgam (:, :, ibz, ie, isb)
                  If (system%np .Gt. 1) Then
                     Do ip = system%np - 1, 1, - 1
                        ip1 = ip + 1
                        zwa (:, :) = zpsi (:, :, ip1) - zgmb (:, :, ip1)
                        Call lcinv (zwa)
                        Call mmatmul_z (zwb, zwa, zskf(:, :, ip, ibz))
                        Call mmatmul_z (zgmb(:, :, ip), zsko(:, :, ip1, ibz), zwb)
                     End Do
                  End If

                  ig = 0 !                              ON-LAYER BLOCKS  AND ON-SITE BLOCKS
                  Do ip = 1, system%np
!!$ blocks of the Green's function for the scattering region
                     zwa (:, :) = zpsi (:, :, ip) - zgmv (:, :, ip) - zgmb (:, :, ip)
                     Call lcinv (zwa)
!$omp critical		
                     Do ib = 1, system%nb
                        ig = ig + 1
                        ii0 = (ib-1) * opt%nlsq
                        zgfg1 (:, :, ig) = zgfg1 (:, :, ig) + kweight (ibz) * zwa (ii0+1:ii0+opt%nlsq, &
                       & ii0+1:ii0+opt%nlsq)
                     End Do
!$omp end critical		     		
                  End Do
               End Do !------------- END OF ibz-LOOP
!$omp end do
               Deallocate (zwa, zwb, zgmv, zgmb, zpsi)
!$omp end parallel

            End Do !------------- ENERGY   LOOP
         End Do !------------- AND SPIN LOOP



!------- SYMMETRIZATION OF THE ON-SITE BLOCKS
         Do is = 1, opt%ns
!$omp parallel do default(shared) private(ig)
            Do ie = 1, opt%ne
               Do ig = 1, system%nums
                  Call symgf (zgfg(:, :, ig, ie, is), opt%bz%nsym, opt%nlsq)
               End Do
            End Do
!$omp end parallel do
         End Do
         Return
      End Subroutine repa

!*******************
!XXX    CPAIT   ****
!*******************
      Subroutine cpait (opt, system, zcn, iiter)
!****************************************************
!!$   PERFORMS ONE CPA-ITERATION, CALCULATES CONFIGURATIONALLY  AVERAGED
!!$   GREENS FUNCTIONS, AND WRITES THE RESULTING CPF FUNCTIONS
!****************************************************
         Use green_func
         Use str_const
         Use inversion ! matrix inversion from LAPACK
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Real (Kind=prec), External :: dzasum ! BLAS func. calculate norm of complex matrix used instead of CMANO (see old code)
         Integer :: nonc (system%nums, opt%ns), ia, is, ie, il, i, j, ig, ic, ista, ifin, eft, modpr, iiter
         Complex (Kind=prec) :: zp (opt%nlsq), zfi1 (opt%nlsq, opt%nlsq), ZCP1 (opt%nlsq, opt%nlsq), zom1 &
        & (opt%nlsq, opt%nlsq), zw (opt%nlsq, opt%nlsq), ze, zpfx, zlax, zmux, zcn (opt%ne)
         Real (Kind=prec), Parameter :: rc1 = 1.0d0, rch = 0.5d0
         Real (Kind=prec) :: qmixi, ppcx, ppdx, ppqx, alfx, sum, sum1
         Type (atom_definition), Pointer :: aa

         qmixi = opt%qmixc
!!$          If (iiter == 1 .and. .Not. (opt%rwgamma == 1 .And. opt%icont == 1)) qmixi = rc1
         If (iiter == 1) qmixi = rc1
         nonc = 0


! ! !SMP$ DO SERIAL
         Do is = 1, opt%ns !-------------------------------- LOOP OVER SPIN AND ENERGY
!$omp parallel do default(private)&
!$omp& shared(is,opt,system,zomg,QMIXI,nonc,IW6)&
!$omp& shared(zcn,zbgam,zvgam,zbski,zbskf,zbsko,zvski,zvskf,zvsko,zskf,zsko,zcpf,zski,zgfg,zcagf)
            Do ie = 1, opt%ne
               ze = zcn (ie)
               ia = 0
               Do ig = 1, system%nums !------------------------------ LOOP OVER SITES

                  zfi1 (:, :) = zgfg (:, :, ig, ie, is)
                  ZCP1 (:, :) = zcpf (:, :, ig, ie, is)
                  zom1 (:, :) = zomg (:, :, ig, ie, is)
                  Call lcinv (zfi1)

                  zw = ZCP1 - zfi1 - zom1

                  sum = dzasum (opt%nlsq*opt%nlsq, zw, 1)
                  sum1 = dzasum (opt%nlsq*opt%nlsq, zom1, 1)

                  If (sum .Ge. (opt%difmc*sum1)) nonc (ig, is) = nonc (ig, is) + 1

!SMP$ DO SERIAL
                  Do j = 1, opt%nlsq
!SMP$ DO SERIAL
                     Do i = 1, opt%nlsq
                        zom1 (i, j) = zom1 (i, j) + qmixi * zw (i, j)
                     End Do
                  End Do

                  Call symgf (zom1, opt%bz%nsym, opt%nlsq)

                  Do ic = 1, system%nc(ig)!--------------------------- LOOP OVER COMPONENTS
                     ia = ia + 1
                     aa => system%at (ia)
                     Do il = 1, opt%nl
                        ppcx = aa%ppc (il, is)
                        ppdx = aa%ppd (il, is)
                        ppqx = aa%ppq (il, is)
                        alfx = system%qscr (il)
                        ista = (il-1) ** 2 + 1
                        ifin = il ** 2
                        Call plmz (ppcx, ppdx, ppqx, alfx, ze, zpfx, zlax, zmux)
!SMP$ DO SERIAL
                        Do i = ista, ifin
                           zp (i) = zpfx
                        End Do
                     End Do

                     zw = - zom1
                     Do i = 1, opt%nlsq
                        zw (i, i) = zw (i, i) + zp (i)
                     End Do

                     Call lcinv (zw)

                     zcagf (:, :, ia, ie, is) = zw (:, :)
                  End Do !------------------------------------------ END OF LOOP OVER COMPONENTS
                  zomg (:, :, ig, ie, is) = zom1 (:, :)

               End Do !-------------------------- END OF LOOP OVER SITES
            End Do !-------------------------- ENERGY
!$omp end parallel do
         End Do !-------------------------- AND SPIN

         Write (IW6, 111) iiter, iiter
         flush (IW6)
111      Format (/ / 5 X, '-----------  ITERATION:  ITER=', I4, '  -----------' / / '  ****  CONVERGENCE OF C&
        &OHERENT INTERACTOR :' / '        IG,      NON-CONV. OMEGA:', 20 X, 'ITER=', I4)
         opt%lnonc = 0
         Do ig = 1, system%nums
            Write (IW6, 120) ig, (nonc(ig, is), is=1, opt%ns)
            opt%lnonc = opt%lnonc + nonc (ig, 1) + nonc (ig, opt%ns)
         End Do
         flush (IW6)
120      Format (5 X, I5, 10 X, 2 I5)

!!$                                  OUTPUT TO UNIT IW8
101      Format (1 X, 10 I5)
104      Format (1 X, 4 G15.7)

         modpr = Mod (iiter, opt%NFUPR)
         If (modpr /= 0 .And. iiter /= 1 .And. iiter /= opt%NITER) Return! concerns full output

         If (opt%ivac == 0) eft = system%bef
         If (opt%ivac == 1) eft = rch * (system%bef+system%vef)

!!$       REWIND IW8
!!$       WRITE(IW8,181)
!!$ 181   FORMAT(2X,' ----------  CPF-FILE  ----------')
!!$       WRITE(IW8,101) system%nums,opt%nl,opt%ns,opt%ne
!!$       WRITE(IW8,104) opt%diam,EFT
!!$       do is = 1, opt%ns
!!$         do ie = 1, opt%ne
!!$           do ig = 1, system%nums
!!$             write (IW8,104) ((zcpf(i,j,ig,ie,is),i=j,opt%nlsq),j=1,opt%nlsq)
!!$           enddo
!!$         enddo
!!$       enddo
!!$       WRITE(IW8,182)
!!$ 182   FORMAT(2X,' --------------------------------')
         Return
      End Subroutine cpait

!*******************
!XXX    plmz    ****
!*******************
      Subroutine plmz (C, DELTA, gamma, ALPHA, ze, zp, ZL, ZM)
!******************************************************
!!$          SCREENED POTENTIAL FUNCTION P(ZE) AND
!!$          RELATED FUNCTIONS LAMBDA(ZE), MU(ZE)
!******************************************************
         Implicit None
         Real (Kind=prec) :: gma, gamma, ALPHA, sqd, DELTA, C
         Complex (Kind=prec) :: zemc, ze, zden, zp, ZL, ZM

         gma = gamma - ALPHA
         sqd = dsqrt (DELTA)
         zemc = ze - C
         zden = DELTA + gma * zemc
         zp = zemc / zden
         ZL = gma / zden
         ZM = sqd / zden
         Return
      End Subroutine plmz


!*******************
!XXX    symgf   ****
!*******************
      Subroutine symgf (za, nsym, nlsq)
!----------------------------------------------------
!!$   SYMMETRIZATION OF LOCAL GF MATRIX ACCORDING TO
!!$   THE POINT GROUP SYMMETRY
!----------------------------------------------------
         Implicit None
         Complex (Kind=prec) :: za (nlsq, nlsq), zg (nlsq, nlsq), zp (nlsq, nlsq), zz
         Real (Kind=prec), Parameter :: rcz = 0.0d0, rch = 0.5d0
         Integer :: i, j, nsym, nlsq

!!$    nl = sqrt(nlsq)
         zz = dcmplx (rcz, rcz)
!SMP$ DO SERIAL
         Do j = 1, nlsq
!SMP$ DO SERIAL
            Do i = 1, nlsq
               zp (i, j) = rch * (za(i, j)+za(j, i))
            End Do
         End Do
!----------------------------------  GENERAL CASE
         If (nsym .Le. 0) Then
!SMP$ DO SERIAL
            Do j = 1, nlsq
!SMP$ DO SERIAL
               Do i = 1, nlsq
                  za (i, j) = zp (i, j)
               End Do
            End Do
            Return
         End If

!----------------------------------  SPECIAL CASES
         zg = zz

         If (nsym == 1) Go To 201
         If (nsym == 2) Go To 202
         If (nsym == 3) Go To 203
         If (nsym == 4) Go To 204

!!$     FCC(001) SYMMETRY - GENERATED BY FOURFOLD ROTATION
!!$     (AXIS Z) AND MIRROR REFLECTION (PLANE X-Z)

201      zg (1, 1) = zp (1, 1)
         zg (3, 1) = zp (3, 1)
         zg (7, 1) = zp (7, 1)
         zg (2, 2) = rch * (zp(2, 2)+zp(4, 4))
         zg (6, 2) = rch * (zp(6, 2)+zp(8, 4))
         zg (3, 3) = zp (3, 3)
         zg (7, 3) = zp (7, 3)
         zg (4, 4) = rch * (zp(2, 2)+zp(4, 4))
         zg (8, 4) = rch * (zp(6, 2)+zp(8, 4))
         zg (5, 5) = zp (5, 5)
         zg (6, 6) = rch * (zp(6, 6)+zp(8, 8))
         zg (7, 7) = zp (7, 7)
         zg (8, 8) = rch * (zp(6, 6)+zp(8, 8))
         zg (9, 9) = zp (9, 9)
         If (nlsq == 9) Go To 270
         zg (13, 1) = zp (13, 1)
         zg (10, 2) = rch * (zp(10, 2)-zp(16, 4))
         zg (12, 2) = rch * (zp(12, 2)+zp(14, 4))
         zg (13, 3) = zp (13, 3)
         zg (14, 4) = rch * (zp(12, 2)+zp(14, 4))
         zg (16, 4) = rch * (-zp(10, 2)+zp(16, 4))
         zg (11, 5) = zp (11, 5)
         zg (10, 6) = rch * (zp(10, 6)-zp(16, 8))
         zg (12, 6) = rch * (zp(12, 6)+zp(14, 8))
         zg (13, 7) = zp (13, 7)
         zg (14, 8) = rch * (zp(12, 6)+zp(14, 8))
         zg (16, 8) = rch * (-zp(10, 6)+zp(16, 8))
         zg (15, 9) = zp (15, 9)
         zg (10, 10) = rch * (zp(10, 10)+zp(16, 16))
         zg (12, 10) = rch * (zp(12, 10)-zp(14, 16))
         zg (11, 11) = zp (11, 11)
         zg (12, 12) = rch * (zp(12, 12)+zp(14, 14))
         zg (13, 13) = zp (13, 13)
         zg (14, 14) = rch * (zp(12, 12)+zp(14, 14))
         zg (16, 14) = rch * (-zp(10, 12)+zp(16, 14))
         zg (15, 15) = zp (15, 15)
         zg (16, 16) = rch * (zp(10, 10)+zp(16, 16))
         Go To 270

!!$     BCC(110) SYMMETRY - GENERATED BY TWO
!!$     MIRROR REFLECTIONS (PLANES X-Z AND Y-Z)

202      zg (1, 1) = zp (1, 1)
         zg (3, 1) = zp (3, 1)
         zg (7, 1) = zp (7, 1)
         zg (9, 1) = zp (9, 1)
         zg (2, 2) = zp (2, 2)
         zg (6, 2) = zp (6, 2)
         zg (3, 3) = zp (3, 3)
         zg (7, 3) = zp (7, 3)
         zg (9, 3) = zp (9, 3)
         zg (4, 4) = zp (4, 4)
         zg (8, 4) = zp (8, 4)
         zg (5, 5) = zp (5, 5)
         zg (6, 6) = zp (6, 6)
         zg (7, 7) = zp (7, 7)
         zg (9, 7) = zp (9, 7)
         zg (8, 8) = zp (8, 8)
         zg (9, 9) = zp (9, 9)
         If (nlsq == 9) Go To 270
         zg (13, 1) = zp (13, 1)
         zg (15, 1) = zp (15, 1)
         zg (10, 2) = zp (10, 2)
         zg (12, 2) = zp (12, 2)
         zg (13, 3) = zp (13, 3)
         zg (15, 3) = zp (15, 3)
         zg (14, 4) = zp (14, 4)
         zg (16, 4) = zp (16, 4)
         zg (11, 5) = zp (11, 5)
         zg (10, 6) = zp (10, 6)
         zg (12, 6) = zp (12, 6)
         zg (13, 7) = zp (13, 7)
         zg (15, 7) = zp (15, 7)
         zg (14, 8) = zp (14, 8)
         zg (16, 8) = zp (16, 8)
         zg (13, 9) = zp (13, 9)
         zg (15, 9) = zp (15, 9)
         zg (10, 10) = zp (10, 10)
         zg (12, 10) = zp (12, 10)
         zg (11, 11) = zp (11, 11)
         zg (12, 12) = zp (12, 12)
         zg (13, 13) = zp (13, 13)
         zg (15, 13) = zp (15, 13)
         zg (14, 14) = zp (14, 14)
         zg (16, 14) = zp (16, 14)
         zg (15, 15) = zp (15, 15)
         zg (16, 16) = zp (16, 16)
         Go To 270

!!$    FCC(111) SYMMETRY - GENERATED BY THREEFOLD ROTATION (AXIS Z)
!!$    AND MIRROR REFLECTION (PLANE Y-Z)

203      zg (1, 1) = zp (1, 1)
         zg (3, 1) = zp (3, 1)
         zg (7, 1) = zp (7, 1)
         zg (2, 2) = rch * (zp(2, 2)+zp(4, 4))
         zg (6, 2) = rch * (zp(6, 2)+zp(8, 4))
         zg (9, 2) = rch * (zp(9, 2)+zp(5, 4))
         zg (3, 3) = zp (3, 3)
         zg (7, 3) = zp (7, 3)
         zg (4, 4) = rch * (zp(2, 2)+zp(4, 4))
         zg (5, 4) = rch * (zp(9, 2)+zp(5, 4))
         zg (8, 4) = rch * (zp(6, 2)+zp(8, 4))
         zg (5, 5) = rch * (zp(5, 5)+zp(9, 9))
         zg (8, 5) = rch * (zp(8, 5)+zp(6, 9))
         zg (6, 6) = rch * (zp(6, 6)+zp(8, 8))
         zg (9, 6) = rch * (zp(9, 6)+zp(5, 8))
         zg (7, 7) = zp (7, 7)
         zg (8, 8) = rch * (zp(6, 6)+zp(8, 8))
         zg (9, 9) = rch * (zp(5, 5)+zp(9, 9))
         If (nlsq == 9) Go To 270
         zg (10, 1) = zp (10, 1)
         zg (13, 1) = zp (13, 1)
         zg (12, 2) = rch * (zp(12, 2)+zp(14, 4))
         zg (15, 2) = rch * (zp(15, 2)+zp(11, 4))
         zg (10, 3) = zp (10, 3)
         zg (13, 3) = zp (13, 3)
         zg (11, 4) = rch * (zp(15, 2)+zp(11, 4))
         zg (14, 4) = rch * (zp(12, 2)+zp(14, 4))
         zg (11, 5) = rch * (zp(11, 5)+zp(15, 9))
         zg (14, 5) = rch * (zp(14, 5)+zp(12, 9))
         zg (12, 6) = rch * (zp(12, 6)+zp(14, 8))
         zg (15, 6) = rch * (zp(15, 6)+zp(11, 8))
         zg (10, 7) = zp (10, 7)
         zg (13, 7) = zp (13, 7)
         zg (11, 8) = rch * (zp(15, 6)+zp(11, 8))
         zg (14, 8) = rch * (zp(12, 6)+zp(14, 8))
         zg (12, 9) = rch * (zp(14, 5)+zp(12, 9))
         zg (15, 9) = rch * (zp(11, 5)+zp(15, 9))
         zg (10, 10) = zp (10, 10)
         zg (13, 10) = zp (13, 10)
         zg (11, 11) = rch * (zp(11, 11)+zp(15, 15))
         zg (14, 11) = rch * (zp(14, 11)+zp(12, 15))
         zg (12, 12) = rch * (zp(12, 12)+zp(14, 14))
         zg (15, 12) = rch * (zp(15, 12)+zp(11, 14))
         zg (13, 13) = zp (13, 13)
         zg (14, 14) = rch * (zp(12, 12)+zp(14, 14))
         zg (15, 15) = rch * (zp(11, 11)+zp(15, 15))
         zg (16, 16) = zp (16, 16)
         Go To 270

!!$    GB (N01) SYMMETRY - GENERATED BY ONE
!!$    MIRROR REFLECTION (PLANE X-Z)

204      zg (1, 1) = zp (1, 1)
         zg (3, 1) = zp (3, 1)
         zg (4, 1) = zp (4, 1)
         zg (7, 1) = zp (7, 1)
         zg (8, 1) = zp (8, 1)
         zg (9, 1) = zp (9, 1)
         zg (2, 2) = zp (2, 2)
         zg (5, 2) = zp (5, 2)
         zg (6, 2) = zp (6, 2)
         zg (3, 3) = zp (3, 3)
         zg (4, 3) = zp (4, 3)
         zg (7, 3) = zp (7, 3)
         zg (8, 3) = zp (8, 3)
         zg (9, 3) = zp (9, 3)
         zg (4, 4) = zp (4, 4)
         zg (7, 4) = zp (7, 4)
         zg (8, 4) = zp (8, 4)
         zg (9, 4) = zp (9, 4)
         zg (5, 5) = zp (5, 5)
         zg (6, 5) = zp (6, 5)
         zg (6, 6) = zp (6, 6)
         zg (7, 7) = zp (7, 7)
         zg (8, 7) = zp (8, 7)
         zg (9, 7) = zp (9, 7)
         zg (8, 8) = zp (8, 8)
         zg (9, 8) = zp (9, 8)
         zg (9, 9) = zp (9, 9)
         If (nlsq == 9) Go To 270
         zg (13, 1) = zp (13, 1)
         zg (14, 1) = zp (14, 1)
         zg (15, 1) = zp (15, 1)
         zg (16, 1) = zp (16, 1)
         zg (10, 2) = zp (10, 2)
         zg (11, 2) = zp (11, 2)
         zg (12, 2) = zp (12, 2)
         zg (13, 3) = zp (13, 3)
         zg (14, 3) = zp (14, 3)
         zg (15, 3) = zp (15, 3)
         zg (16, 3) = zp (16, 3)
         zg (13, 4) = zp (13, 4)
         zg (14, 4) = zp (14, 4)
         zg (15, 4) = zp (15, 4)
         zg (16, 4) = zp (16, 4)
         zg (10, 5) = zp (10, 5)
         zg (11, 5) = zp (11, 5)
         zg (12, 5) = zp (12, 5)
         zg (10, 6) = zp (10, 6)
         zg (11, 6) = zp (11, 6)
         zg (12, 6) = zp (12, 6)
         zg (13, 7) = zp (13, 7)
         zg (14, 7) = zp (14, 7)
         zg (15, 7) = zp (15, 7)
         zg (16, 7) = zp (16, 7)
         zg (13, 8) = zp (13, 8)
         zg (14, 8) = zp (14, 8)
         zg (15, 8) = zp (15, 8)
         zg (16, 8) = zp (16, 8)
         zg (13, 9) = zp (13, 9)
         zg (14, 9) = zp (14, 9)
         zg (15, 9) = zp (15, 9)
         zg (16, 9) = zp (16, 9)
         zg (10, 10) = zp (10, 10)
         zg (11, 10) = zp (11, 10)
         zg (12, 10) = zp (12, 10)
         zg (11, 11) = zp (11, 11)
         zg (12, 11) = zp (12, 11)
         zg (12, 12) = zp (12, 12)
         zg (13, 13) = zp (13, 13)
         zg (14, 13) = zp (14, 13)
         zg (15, 13) = zp (15, 13)
         zg (16, 13) = zp (16, 13)
         zg (14, 14) = zp (14, 14)
         zg (15, 14) = zp (15, 14)
         zg (16, 14) = zp (16, 14)
         zg (15, 15) = zp (15, 15)
         zg (16, 15) = zp (16, 15)
         zg (16, 16) = zp (16, 16)
         Go To 270
!SMP$ DO SERIAL
270      Do i = 1, nlsq - 1
!SMP$ DO SERIAL
            Do j = i + 1, nlsq
               zg (i, j) = zg (j, i)
            End Do
         End Do

         za = zg
         Return
      End Subroutine symgf

!*******************
!XXX    CMOM    ****
!*******************
      Subroutine cmom (opt, system, zcn, zcw)
!***********************************************************
!!$ CALCULATES MOMENTS OF DOS VIA COMPLEX CONTOUR INTEGRATION
!!$ OF CONFIGURATIONALLY AVERAGED PHYSICAL GREENS FUNCTIONS
!***********************************************************
         Use green_func
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Complex (Kind=prec) :: ZGDI (opt%ne, opt%nl, opt%ns, system%num), ZGOF (opt%ne, opt%pair, opt%ns, &
        & system%num), zw0 (opt%ne), zw1 (opt%ne), zw2 (opt%ne), zw00 (opt%ne), zw01 (opt%ne), zw10 (opt%ne), &
        & zw11 (opt%ne), zw02 (opt%ne), zw20 (opt%ne), zz, ze, zpf1, zla1, zmu1, zdp1, zeps1, zeps2, zpf2, &
        & zla2, zmu2, zcn (opt%ne), zcw (opt%ne)
         Real (Kind=prec), Parameter :: rcz = 0.0d0
         Real (Kind=prec) :: ppc1, ppd1, ppq1, alf1, ppc2, ppd2, ppq2, alf2, eny1, eny2
         Integer :: i, ia, is, il, ista, ifin, ie, ipair, il1, il2, l1, l2, icen1, icen2, m, i1, I2
         Type (atom_definition), Pointer :: aa

         zz = dcmplx (rcz, rcz)
         ZGDI = zz
         ZGOF = zz
!--------------------------------- DIAGONAL ELEMENTS
!$omp parallel do default(private)&
!$omp& shared(opt,system,zcn,zcagf,zgdi,IW6)
         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%ns
               Do il = 1, opt%nl
                  ppc1 = aa%ppc (il, is)
                  ppd1 = aa%ppd (il, is)
                  ppq1 = aa%ppq (il, is)
                  alf1 = system%qscr (il)
                  ista = (il-1) ** 2 + 1
                  ifin = il ** 2
                  Do ie = 1, opt%ne
                     ze = zcn (ie)
                     Call plmz (ppc1, ppd1, ppq1, alf1, ze, zpf1, zla1, zmu1)
                     zdp1 = zmu1 ** 2
!SMP$ DO SERIAL
                     Do i = ista, ifin
                        ZGDI (ie, il, is, ia) = ZGDI (ie, il, is, ia) + zla1 + zdp1 * zcagf (i, i, ia, ie, &
                       & is)
                     End Do
                  End Do
               End Do
            End Do
         End Do
!$omp end parallel do

!----------------------------- OFF-DIAGONAL ELEMENTS
!$omp parallel do default(private)&
!$omp& shared(opt,system,zcn,zcagf,ZGOF,IW6)
         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%ns
               ipair = 0
               Do il1 = 1, opt%nl - 1
                  il2 = il1 + 1
                  l1 = il1 - 1
                  l2 = il1
                  icen1 = l1 ** 2 + l1 + 1
                  icen2 = l2 ** 2 + l2 + 1
                  ppc1 = aa%ppc (il1, is)
                  ppd1 = aa%ppd (il1, is)
                  ppq1 = aa%ppq (il1, is)
                  alf1 = system%qscr (il1)
                  ppc2 = aa%ppc (il2, is)
                  ppd2 = aa%ppd (il2, is)
                  ppq2 = aa%ppq (il2, is)
                  alf2 = system%qscr (il2)
                  Do m = - l1, l1
                     ipair = ipair + 1
                     i1 = icen1 + m
                     I2 = icen2 + m
                     Do ie = 1, opt%ne
                        Call plmz (ppc1, ppd1, ppq1, alf1, zcn(ie), zpf1, zla1, zmu1)
                        Call plmz (ppc2, ppd2, ppq2, alf2, zcn(ie), zpf2, zla2, zmu2)
                        ZGOF (ie, ipair, is, ia) = zmu1 * zcagf (i1, I2, ia, ie, is) * zmu2
                     End Do
                  End Do
               End Do
            End Do
         End Do
!$omp end parallel do

!---------------------------------- INTEGRATIONS
!$omp parallel do default(private)&
!$omp& shared(opt,system,zcn,zcw,zgdi,IW6)
         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%ns
               Do il = 1, opt%nl
                  eny1 = aa%eny (il, is)
!SMP$ DO SERIAL
                  Do ie = 1, opt%ne
                     zw0 (ie) = ZGDI (ie, il, is, ia)
                     zw1 (ie) = zw0 (ie) * (zcn(ie)-eny1)
                     zw2 (ie) = zw1 (ie) * (zcn(ie)-eny1)
                  End Do
                  aa%emdi0 (il, is) = rcint (zw0, zcw, opt%ne)
                  aa%emdi1 (il, is) = rcint (zw1, zcw, opt%ne)
                  aa%emdi2 (il, is) = rcint (zw2, zcw, opt%ne)
               End Do
            End Do
         End Do
!$omp end parallel do

!$omp parallel do default(private)&
!$omp& shared(opt,system,zcn,zcw,zgof,IW6)
         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%ns
               ipair = 0
               Do il1 = 1, opt%nl - 1
                  il2 = il1 + 1
                  l1 = il1 - 1
                  eny1 = aa%eny (il1, is)
                  eny2 = aa%eny (il2, is)
                  Do m = - l1, l1
                     ipair = ipair + 1
                     Do ie = 1, opt%ne
                        zeps1 = zcn (ie) - eny1
                        zeps2 = zcn (ie) - eny2
                        zw00 (ie) = ZGOF (ie, ipair, is, ia)
                        zw10 (ie) = zeps1 * zw00 (ie)
                        zw01 (ie) = zeps2 * zw00 (ie)
                        zw11 (ie) = zeps1 * zw01 (ie)
                        zw20 (ie) = zeps1 * zw10 (ie)
                        zw02 (ie) = zeps2 * zw01 (ie)
                     End Do
                     aa%emof00 (ipair, is) = rcint (zw00, zcw, opt%ne)
                     aa%emof10 (ipair, is) = rcint (zw10, zcw, opt%ne)
                     aa%emof01 (ipair, is) = rcint (zw01, zcw, opt%ne)
                     aa%emof11 (ipair, is) = rcint (zw11, zcw, opt%ne)
                     aa%emof20 (ipair, is) = rcint (zw20, zcw, opt%ne)
                     aa%emof02 (ipair, is) = rcint (zw02, zcw, opt%ne)
                  End Do
               End Do
            End Do
         End Do
!$omp end parallel do

         Return
      End Subroutine cmom

!*******************
!XXX    rcint   ****
!*******************
      Function rcint (zf, zcw, ne)
!********************************************************
!!$ INTEGRAL (DIVIDED BY 2*pi*i) OF A FUNCTION OVER THE
!!$  CLOSED COMPLEX CONTOUR WITH NODES AND WEIGHTS
!********************************************************
         Implicit None
         Integer :: ie, ne
         Complex (Kind=prec) :: zf (ne), zcw (ne), zd
         Real (Kind=prec) :: rcint, rcz = 0.0d0
         rcint = rcz
!SMP$ DO SERIAL
         Do ie = 1, ne
            zd = zcw (ie) * zf (ie)
            rcint = rcint + DBLE (zd)
         End Do
         Return
      End Function rcint

!*******************
!XXX   IDAMA1   ****
!*******************
      Function IDAMA1 (DX, ND, N1)
!-----------------------------------------------
!!$     FINDS THE INDEX OF MAX. ELEMENT OF DX(i)
!!$     THE MAXIMUM is SEARCHED OVER N1.LE.i.LE.ND
!-----------------------------------------------
         Implicit None
         Integer :: i, IDAMA1, N1, ND
         Real (Kind=8) :: DX (ND), wmax

         IDAMA1 = N1
         If (N1 == ND) Return
         wmax = DX (N1)
!SMP$ DO SERIAL
         Do i = N1 + 1, ND
            If (DX(i) .Le. wmax) Cycle
            IDAMA1 = i
            wmax = DX (i)
         End Do
         Return
      End Function IDAMA1

End Module green_mod

