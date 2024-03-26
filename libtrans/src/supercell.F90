#include "math_def.h"
Module supercell
      Use geometry_module
      Use omta_defs
      Use omta_SOC
      Use df
      Implicit None

      Type t_supercell
         Integer :: size = 0, no = 0, need_split = 0
         Real (Kind=DEF_DBL_PREC), Pointer :: kbase (:, :), tr_vec (:, :)
         Real (Kind=DEF_DBL_PREC) :: min_diff
         Integer, Pointer :: nreal (:, :), idx (:, :)
         Complex (Kind=DEF_DBL_PREC) :: rot (2, 2) = reshape ( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/ 2, 2 /))
      End Type t_supercell

      Type t_decomp_data
         Integer :: csz = 0, num = 0, ldwf = 0
         Integer :: alloc = 0
         Integer, Allocatable :: nin (:)
         Real (Kind=DEF_DBL_PREC), Allocatable :: vel (:), kp (:, :)
         Complex (Kind=DEF_DBL_PREC), Allocatable :: wf (:, :)

      End Type t_decomp_data
Contains


      Subroutine decomp_free (decdat)
         Implicit None
         Type (t_decomp_data) :: decdat
         If (decdat%alloc == 0) Return
         If (decdat%num > 0) Deallocate (decdat%kp, decdat%vel, decdat%wf)
         decdat%num = 0
         Deallocate (decdat%nin)
         decdat%alloc = 0
         decdat%csz = 0
         decdat%ldwf = 0
      End Subroutine decomp_free

      Subroutine decomp_realloc (csz, num, decdat, ldwf)
         Implicit None
         Type (t_decomp_data) :: decdat
         Integer :: csz, num, ldwf
         If ((decdat%alloc /= 0) .And. (decdat%num /= num .Or. decdat%csz /= csz)) Then
            Call decomp_free (decdat)
         End If

         If (decdat%alloc == 0) Then
            Allocate (decdat%nin(csz))
            decdat%csz = csz
            decdat%num = num
            decdat%ldwf = ldwf
            decdat%alloc = 1
            If (num == 0) Return
            Allocate (decdat%vel(num), decdat%kp(2, num), decdat%wf(ldwf, num))
         End If

      End Subroutine decomp_realloc

      Subroutine fill_decomp (csz, n, s_ando, kbase, kpar, decdat)
         Use ando_module
         Implicit None
         Type (t_decomp_data) :: decdat
         Integer :: csz, n
         Real (Kind=DEF_DBL_PREC) :: kbase (:, :), kpar (:)
         Type (t_ando_sollution) :: s_ando (:)
!!$ Local vars
         Integer :: i, j, num, k

         num = 0
         Do i = 1, csz
            num = num + s_ando(i)%nin
         End Do

         Call decomp_realloc (csz, num, decdat, n)
!!$          If (num == 0) Return

         k = 1
         Do i = 1, csz
            decdat%nin (i) = s_ando(i)%nin
            Do j = 1, decdat%nin(i)
               decdat%kp (1:2, k) = kbase (1:2, i) + kpar (:)
               decdat%vel (k) = s_ando(i)%Vin(j)
               decdat%wf (:, k) = s_ando(i)%Uin(1:n, j)
               k = k + 1
            End Do
         End Do
      End Subroutine fill_decomp


      Subroutine prep_dec_idx (geo, idx)
         Use geometry_module
         Implicit None
         Type (t_geometry) :: geo
         Integer, Pointer :: idx (:, :)
!!$ Local vars
         Integer :: i, no, j, ind, ofs

         Allocate (idx(2, geo%norbit))
         ofs = 1
         ind = 1
         Do i = 1, geo%num
            no = geo%atoms(i)%ptr%nm
            Do j = 1, no
               idx (1, ind) = ofs
               idx (2, ind) = ofs + no
               ind = ind + 1
               ofs = ofs + 1
            End Do
            ofs = ofs + no
         End Do
      End Subroutine prep_dec_idx

      Subroutine prep_dec_idx_EMTO (geo, idx)
         Use geometry_module
         Implicit None
         Type (t_geometry_EMTO) :: geo
         Integer, Pointer :: idx (:, :)
!!$ Local vars
         Integer :: i, no, j, ind, ofs!!$, spinoffset
         Integer :: numls(5) = [1,3,5,7,9]


         !num=size(geo%to_downfold)/geo%nspin /geo%num

         Allocate (idx(2, geo%norbit-size(geo%to_downfold)/geo%nspin))
         ofs = 1
         ind = 1
         Do i = 1, geo%num
         no = geo%atoms(i)%ptr%nm - dot_product(geo%atoms(i)%ptr%idxdn-1,numls(1:geo%atoms(i)%ptr%lmx+1)) 
            Do j = 1, no
               idx (1, ind) = ofs
               idx (2, ind) = ofs + no
               ind = ind + 1
               ofs = ofs + 1
            End Do
            ofs = ofs + no
         End Do

      End Subroutine prep_dec_idx_EMTO


      Subroutine prep_ham_idx (ngr, deg, idx)
         Implicit None
         Integer :: ngr, deg
         Integer, Pointer :: idx (:, :)
         Integer :: i, no, j, ind, ofs
         Allocate (idx(2, ngr*deg))
         ofs = 1
         ind = 1
         Do i = 1, ngr
            no = deg
            Do j = 1, no
               idx (1, ind) = ofs
               idx (2, ind) = ofs + no
               ind = ind + 1
               ofs = ofs + 1
            End Do
            ofs = ofs + no
         End Do
      End Subroutine prep_ham_idx

      Subroutine prepare_supercell (geom, scell, irel, split, rot, ngrid)
         Use geometry_module
         Implicit None
         Type (t_supercell) :: scell
         Type (t_geometry) :: geom
         Integer, Optional :: irel, split, ngrid
         Complex (Kind=DEF_DBL_PREC), Optional :: rot (2, 2)
!!$ Local
         Integer :: size, ink, n1, n2, lirel
         Real (Kind=DEF_DBL_PREC) :: cnorm, det
         Real (Kind=DEF_DBL_PREC) :: scbase (2, 2), vbg (2, 2)


         lirel = 0
         If (present(irel)) lirel = irel

         If (lirel /= 0) Then
            scell%need_split = 0
            If (present(split)) scell%need_split = split

            If ( .Not. present(rot)) Then
               scell%rot = reshape ( (/ 1.0d0, 0.0d0, 0.0d0, 1.0d0 /), (/ 2, 2 /))
            Else
               scell%rot = rot
            End If
         End If

         If ( .Not. present(ngrid)) Then
            If (scell%need_split /= 0) Call prep_dec_idx (geom, scell%idx)
         Else
            If (scell%need_split /= 0) Call prep_ham_idx (ngrid, 2, scell%idx)
         End If

         size = geom%sc_size (1) * geom%sc_size(2)
         scell%size = size
         Allocate (scell%kbase(2, size))
         Allocate (scell%tr_vec(2, size))
         Allocate (scell%nreal(size, 2))


         scbase (:, 1) = geom%base(:, 1) * geom%sc_size(1)
         scbase (:, 2) = geom%base(:, 2) * geom%sc_size(2)
         scell%no = real (size, kind=DEF_DBL_PREC) * geom%norbit * (lirel+1)
         det = scbase (1, 1) * scbase (2, 2) - scbase (1, 2) * scbase (2, 1)
         cnorm = 2.0d0 * DEF_M_PI / det
         vbg (1, 1) = cnorm * scbase (2, 2)
         vbg (2, 1) = - cnorm * scbase (1, 2)
         vbg (1, 2) = - cnorm * scbase (2, 1)
         vbg (2, 2) = cnorm * scbase (1, 1)

         ink = 1
         Do n1 = 0, geom%sc_size(1) - 1
            Do n2 = 0, geom%sc_size(2) - 1
               scell%kbase (:, ink) = n1 * vbg (:, 1) + n2 * vbg (:, 2)
               ink = ink + 1
            End Do
         End Do

         ink = 1
         Do n1 = 0, geom%sc_size(2) - 1
            Do n2 = 0, geom%sc_size(1) - 1
               scell%tr_vec (:, ink) = n2 * geom%base(:, 1) + n1 * geom%base(:, 2)
               ink = ink + 1
            End Do
         End Do

!!$ Setting the minimum distance between supercell-equivalent k|| points
         If (size > 1) Then
            scell%min_diff = Min (Sqrt(dot_product(vbg(:, 1), vbg(:, 1))), Sqrt(dot_product(vbg(:, 2), vbg(:, &
           & 2))))
         Else
            scell%min_diff = 1.0d0
         End If

      End Subroutine prepare_supercell

      Subroutine prepare_supercell_EMTO (geom, scell, split, nspin)
         Implicit None
         Type (t_supercell) :: scell
         Type (t_geometry_EMTO) :: geom
         Integer, Optional :: split
         Integer :: nspin
!!$ Local
         Integer :: sc_size, ink, n1, n2
         Real (Kind=DEF_DBL_PREC) :: cnorm, det
         Real (Kind=DEF_DBL_PREC) :: scbase (2, 2), vbg (2, 2)

         scell%need_split = split
         If (split==1) Call prep_dec_idx_EMTO (geom, scell%idx)

         sc_size = geom%sc_size (1) * geom%sc_size(2)
         scell%size = sc_size
         Allocate (scell%kbase(2, sc_size))
         Allocate (scell%tr_vec(2, sc_size))
         Allocate (scell%nreal(sc_size, 2))

         scbase (:, 1) = geom%base(:, 1) * geom%sc_size(1)
         scbase (:, 2) = geom%base(:, 2) * geom%sc_size(2)
         scell%no = real (sc_size, kind=DEF_DBL_PREC) * (geom%norbit - size(geom%to_downfold)/nspin) * (nspin)
         !scell%no = real (sc_size, kind=DEF_DBL_PREC) * geom%norbit * (nspin)

         

         det = scbase (1, 1) * scbase (2, 2) - scbase (1, 2) * scbase (2, 1)
         cnorm = 2.0d0 * DEF_M_PI / det
         vbg (1, 1) = cnorm * scbase (2, 2)
         vbg (2, 1) = - cnorm * scbase (1, 2)
         vbg (1, 2) = - cnorm * scbase (2, 1)
         vbg (2, 2) = cnorm * scbase (1, 1)

         ink = 1
         Do n1 = 0, geom%sc_size(1) - 1
            Do n2 = 0, geom%sc_size(2) - 1
               scell%kbase (:, ink) = n1 * vbg (:, 1) + n2 * vbg (:, 2)
               ink = ink + 1
            End Do
         End Do

         ink = 1
         Do n1 = 0, geom%sc_size(2) - 1
            Do n2 = 0, geom%sc_size(1) - 1
               scell%tr_vec (:, ink) = n2 * geom%base(:, 1) + n1 * geom%base(:, 2)
               ink = ink + 1
            End Do
         End Do

!!$ Setting the minimum distance between supercell-equivalent k|| points
         If (sc_size > 1) Then
            scell%min_diff = Min (Sqrt(dot_product(vbg(:, 1), vbg(:, 1))), Sqrt(dot_product(vbg(:, 2), vbg(:, &
           & 2))))
         Else
            scell%min_diff = 1.0d0
         End If


      End Subroutine prepare_supercell_EMTO

      Function get_sc_ando (sp, ando, kpar, ppar, sc, scell, po, side, actl, rotm, decd) Result (Nin_out)
!!$ Returns Ando sollution for supercell
         Use ando_module
         Use hamiltonian
         Use structure
         Implicit None
         External zgetrf, zgetfi

         Type (t_supercell) :: scell
         Type (t_ando_sollution) :: ando
         Real (Kind=DEF_DBL_PREC) :: kpar (2)
         Integer :: sp, side
         Type (t_strconst) :: sc
         Type (t_potpar_mats) :: ppar
         Type (t_pot_opts) :: po
         Type (t_ando_ctl) :: actl
         Type (t_decomp_data), Optional :: decd
         Real (Kind=DEF_DBL_PREC), Optional :: rotm (:, :)
!!$ Local
         Type (t_mathop), Target :: sys
         Type (t_ando_options) :: aopts
         Type (zcsrmat), Pointer :: hop_ptr
         Type (t_ando_sollution), Pointer :: s_ando (:)
         Type (t_ando_sollution), Pointer :: t_ando
         Complex (Kind=DEF_DBL_PREC), Pointer :: ph (:, :)
         Complex (Kind=DEF_DBL_PREC), Pointer :: ph_i (:, :)
         Complex (Kind=DEF_DBL_PREC), Pointer :: work (:)

         Real (Kind=DEF_DBL_PREC), Pointer :: kbase (:, :), tr_vec (:, :)
         Complex (Kind=DEF_DBL_PREC) :: k1 (2), ph1, phi1
         Complex (Kind=DEF_DBL_PREC), Pointer :: ons (:, :), ofs (:, :), hops (:, :)
         Integer :: csz, i, j, info, ipiv (scell%size), nin, ind1, ind2, jnd1, jnd2, sz1, l, n, nh, Nin_out
         Integer :: n1, n0, pohops, ngc (3), ngr (3)
         Complex (Kind=DEF_DBL_PREC), Pointer :: emb (:, :), t_emb (:, :)
         Real (Kind=DEF_DBL_PREC) :: scl

         Allocate (ph(scell%size, scell%size), ph_i(scell%size, scell%size), work(scell%size), &
        & s_ando(scell%size))

         pohops = 1
         If (po%kind == 3) pohops = 2
         If ((po%kind <= 3) .Or. (po%kind == 7)) Then
            n = sc%nrows / pohops
         Else If (po%kind == 5) Then
            n = scell%no / scell%size
         End If

         nh = pohops
         aopts%use_ev_solver = actl%gensol
         aopts%dir = - side
         aopts%usestates = actl%specpart
         aopts%refine = actl%refine
         If (po%need_current_pol == 1) Then
            aopts%need_pol = 1
            If (po%kind == 5) Then
               aopts%gamnas = 1
            Else
               aopts%gamnas = 0
            End If
         End If
         If (scell%need_split /= 0) Then
            aopts%need_split = 1
            aopts%split_idx => scell%idx
            aopts%rotm = scell%rot
         End If

         If ((po%kind <= 3) .Or. (po%kind == 7)) Then
            sz1 = sc%nrows
         Else If (po%kind == 5) Then
            sz1 = n
         End If

!!$ Prepare k_par and phase factors for supercell
         csz = scell%size
         kbase => scell%kbase
         tr_vec => scell%tr_vec
         Do i = 1, csz
            k1 = kbase (:, i) + kpar
            Do j = 1, csz
               ph (j, i) = Exp (DEF_cmplx_Ione*dot_product(k1, tr_vec(:, j)))
            End Do
         End Do

         ph_i = ph

!!$ Invert the ph_i a matrix
         Call zgetrf (csz, csz, ph_i, csz, ipiv, info)
         If (info == 0) Then
            Call zgetri (csz, ph_i, csz, ipiv, work, csz, info)
         End If

         If (info /= 0) Then
            Write (*,*) 'Problems with phase inversion in get_sc_ando, info =', info
            Stop
         End If

         nin = 0
!$omp parallel default(shared) private(sys,ons,ofs,hops,hop_ptr)  reduction(+:Nin)
         Call init_mathop (sys)
         !hop_ptr => sys%r
         Allocate (ons(sz1, sz1), ofs(sz1, sz1), hops(n, n*(pohops+1)))
!$omp do
         Do i = 1, csz
!!$ Prepare "hamiltonian" for small cell
            If ((po%kind <= 3) .Or. (po%kind == 7)) Then
               Call prep_system (sys, sc, kbase(:, i)+kpar, ppar, po, 1)

            Else If (po%kind == 5) Then
               If (side ==-1) Then
                  ngc (1) = po%nlx
                  ngc (2) = po%nly
                  ngc (3) = po%nlz
                  ngr (:) = ngc (:)
               Else
                  ngc (1) = po%nrx
                  ngc (2) = po%nry
                  ngc (3) = po%nrz
                  ngr (:) = ngc (:)
               End If
               Call prep_system (sys, po, ngc, ngr, kbase(:, i)+kpar, rotm)
            End If

            ons = sptofull (sys%c)

            ofs = sptofull (sys%r)
            hops (1:n, 1:n) = ons (1:n, 1:n)
            hops (1:n, n+1:n+nh*n) = ofs (n*nh-n+1:n*nh, 1:nh*n) 

!!$ this just fix for some nonstandard compilers
            s_ando(i)%alloc = 0
!!$ Solve ando problem
            Call solve_ando (ofs, hops, n, nh, s_ando(i), aopts)
            nin = nin + s_ando(i)%nin
            scell%nreal (i, sp) = s_ando(i)%nin
         End Do
!$omp end do
         Call free_mathop (sys)
         Deallocate (ons, ofs, hops)
!$omp end parallel

!!$ for chde

         If (actl%writedecomp > 0 .And. present(decd)) Then
            Call fill_decomp (csz, n, s_ando, kbase, kpar, decd)
         End If
!!$ for chde


!!$ Allocate memory for "ando" structure for supercell.
         Call alloc_ando (ando, scell%no, nin, nin, aopts)
         ando%emb = DEF_cmplx_zero
         ando%bound = DEF_cmplx_zero
         ando%Uout_i = DEF_cmplx_zero
         ando%Uin = DEF_cmplx_zero
         ando%SNin = 0
         ando%SNout = 0
         If (aopts%need_split == 1) Then
            ando%mask_in = 0
            ando%mask_out = 0
         End If
         ando%ppj = 0.d0
!!$ Construct ando for supercell
         If (nin > 0) Then

!!$ ! !******$omp parallel private(t_ando,ind1,ind2,jnd1,jnd2,ph1,phi1,j) shared(s_ando,ando,ph,ph_i,csz,sz1)
            ind1 = 0
            ind2 = 0
!!$ ! !******$omp do
            Do i = 1, csz
               If (s_ando(i)%nin == 0) Cycle
               t_ando => s_ando (i)
               ando%ppj = ando%ppj + t_ando%ppj * dble (t_ando%nin)
               If (aopts%need_split /= 0) Then
                  n1 = t_ando%SNin (1)
                  If (n1 > 0) Then
                     n0 = ando%SNin (1)
                     ando%mask_in (n0+1:n0+n1:1, 1) = ind2 + t_ando%mask_in(1:n1:1, 1)
                  End If
                  n1 = t_ando%SNin (2)
                  If (n1 > 0) Then
                     n0 = ando%SNin (2)
                     ando%mask_in (n0+1:n0+n1:1, 2) = ind2 + t_ando%mask_in(1:n1:1, 2)
                  End If
                  ando%SNin = ando%SNin + t_ando%SNin

                  n1 = t_ando%SNout (1)
                  If (n1 > 0) Then
                     n0 = ando%SNout (1)
                     ando%mask_out (n0+1:n0+n1:1, 1) = ind2 + t_ando%mask_out(1:n1:1, 1)
                  End If

                  n1 = t_ando%SNout (2)
                  If (n1 > 0) Then
                     n0 = ando%SNout (2)
                     ando%mask_out (n0+1:n0+n1:1, 2) = ind2 + t_ando%mask_out(1:n1:1, 2)
                  End If
                  ando%SNout = ando%SNout + t_ando%SNout
               End If
               ind1 = ind2 + 1
               ind2 = ind2 + t_ando%nin
               jnd1 = 0
               jnd2 = 0
               Do j = 1, csz
                  ph1 = ph (j, i)
                  phi1 = ph_i (i, j)
                  jnd1 = jnd2 + 1
                  jnd2 = jnd2 + sz1
                  ando%Uin (jnd1:jnd2, ind1:ind2) = t_ando%Uin * ph1
                  ando%bound (jnd1:jnd2, ind1:ind2) = t_ando%bound * ph1
                  ando%Uout_i (ind1:ind2, jnd1:jnd2) = t_ando%Uout_i * phi1
               End Do
               ando%Vin (ind1:ind2) = t_ando%Vin(:)
               ando%Vout (ind1:ind2) = t_ando%Vout(:)
            End Do
            ando%ppj = ando%ppj / dble(nin)
!!$ ! !******$omp end do
!!$ ! !******$omp end parallel
         End If
!$omp parallel do private(jnd1,jnd2,t_emb,emb,j,l,i) firstprivate(csz,sz1) shared(ando,ph,ph_i) default(shared)
         Do i = 1, csz
            emb => ando%emb (:, 1+sz1*(i-1) :sz1*i)
            do j = 1, csz
               t_emb => s_ando(j)%emb
               jnd1 = 0
               jnd2 = 0

               do l = 1, csz
                  jnd1 = jnd2 + 1
                  jnd2 = jnd2 + sz1
                  emb (jnd1:jnd2, :) = emb (jnd1:jnd2, :) + ph (l, j) * t_emb * ph_i (j, i)
               end do
            end do
         End do
!$omp end parallel do
         scl = sqrt (1.0d0*csz)
         i = scell%no * nin
         Call zdscal (i, 1.0d0/scl, ando%Uin, 1)
         Call zdscal (i, 1.0d0/scl, ando%bound, 1)
         Call zdscal (i, scl, ando%Uout_i, 1)
!!$ free ando structures for small cells
         Do i = 1, csz
            call free_ando (s_ando(i))
         End do
         Nin_out = nin

         Deallocate (ph, ph_i, work, s_ando)
         Nullify (aopts%split_idx)
      End Function get_sc_ando

!      Subroutine build_vshift(shift, vshift, n)
!        Implicit none
!        Integer :: n
!        Real (kind=DEF_DBL_PREC) :: shift
!        Complex (kind=DEF_DBL_PREC) :: vshift(n,n)
!
!        vshift = 0d0
!        vshift(2,2) =  -1 * shift
!        vshift(3,3) =  0 * shift
!        vshift(4,4) =  1 * shift 
!        If (n > 4) Then
!          vshift(5,5) =  -2 * shift
!          vshift(6,6) =  -1 * shift
!          vshift(7,7) =  0 * shift
!          vshift(8,8) =  1 * shift 
!          vshift(9,9) =  2 * shift 
!          If (n > 9) Then
!            vshift(10,10) = -1 * shift
!            vshift(11,11) = 0 * shift
!            vshift(12,12) =  1 * shift
!            vshift(13,13) =  2 * shift
!            vshift(14,14) =  3 * shift 
!            vshift(15,15) =  4 * shift 
!            vshift(16,16) =  5 * shift 
!          Endif
!        Endif
!      End Subroutine
!
!      Subroutine convert_vshift_rylm(vshift, n)
!       Implicit None
!       Integer :: n
!       Complex (Kind=DEF_DBL_PREC) :: vshift(n,n)
!! local
!       Complex (Kind=DEF_DBL_PREC) :: R(n,n), Ri(n,n)
!       Complex (Kind=DEF_DBL_PREC) :: new_vshift(n,n)
!       Integer :: l, nlm, lmax
!  
!       new_vshift = 0d0
!
!       lmax = sqrt(real(int(n)/2))-1
!
!       Do l=1, lmax
!        nlm = 2*(l+1)-1
!        call transmat(R(1:nlm, 1:nlm), l)
!        Ri(1:nlm, 1:nlm) = transpose(conjg(R(1:nlm, 1:nlm)))
!        new_vshift(l**2+1:(l+1)**2,l**2+1:(l+1)**2) = 0.5d0*matmul(R(1:nlm, 1:nlm), &
!& matmul(vshift(l**2+1:(l+1)**2,l**2+1:(l+1)**2), Ri(1:nlm, 1:nlm)))
!       Enddo
!
!       vshift = new_vshift
!       !finally, duplicate for downspin
!       vshift(n/2+1:n, n/2+1:n) = vshift(1:n/2, 1:n/2)
!
!       Contains
!        Subroutine transmat(R,l)
!!!$ Transformation matrix complex->real harmonics (without 1/sqrt(2) prefactor!!!)
!         Implicit None
!         Integer :: l
!         Complex(Kind=prec) :: R(-l:l, -l:l)
!!!$
!         Integer :: M
!         R(:, :) = DEF_cmplx_zero
!         R(0, 0) = DEF_cmplx_one*Sqrt(2.0d0)
!         Do M = 1, l, 1
!            R(-M, -M) = DEF_cmplx_Ione
!            R(M, M) = (-1.0d0)**M
!            R(M, -M) = DEF_cmplx_one
!            R(-M, M) = -DEF_cmplx_Ione*(-1.0d0)**M
!         End Do
!      End Subroutine transmat
!
!      End Subroutine

      Function get_sc_ando_EMTO (po, nspin, ando, kpar, ppar, geo, sc, scell, side, gens, frac) Result (Nin_out)
!!$ Returns FULL (both spins) Ando sollution for supercell
         Use omta_kink
         Use ando_module
         use structure
         use hamiltonian
         Implicit None
         External zgetf, zgetri

         Type (t_supercell) :: scell
         Type (t_ando_sollution) :: ando
         Type (t_pot_opts) :: po
         Real (Kind=DEF_DBL_PREC) :: kpar (2)
         Real (Kind=DEF_DBL_PREC) :: frac 
         Integer :: sp, side, gens
         type (t_geometry_EMTO)    :: geo
         Type (t_strconst) :: sc
         Type (t_omta_logder) :: ppar(2)
!!$ Local
         Type (t_mathop) :: sys
         Type (t_ando_options) :: aopts
         Type (t_ando_sollution), Pointer :: s_ando (:)
         Type (t_ando_sollution), Pointer :: t_ando
         Type (zcsrmat), Pointer :: hop_ptr
         Complex (Kind=DEF_DBL_PREC), Pointer :: ph (:, :)
         Complex (Kind=DEF_DBL_PREC), Pointer :: ph_i (:, :)
         Complex (Kind=DEF_DBL_PREC), Pointer :: work (:)

         Real (Kind=DEF_DBL_PREC), Pointer :: kbase (:, :), tr_vec (:, :)
         Complex (Kind=DEF_DBL_PREC) :: k1 (2), ph1, phi1
         Complex (Kind=DEF_DBL_PREC), Pointer :: ons (:, :), ofs (:, :), hops (:, :)
         Integer :: csz, i, j, info, ipiv (scell%size), Nin, ind1, ind2, jnd1, jnd2, l, n, ndim, nh, Nin_out
         Integer :: n1, n0, pohops
         Complex (Kind=DEF_DBL_PREC), Pointer :: emb (:, :), t_emb (:, :)
         Real (Kind=DEF_DBL_PREC) :: scl
         Integer :: temp1, temp2
!!$these are related to SOC
         Integer :: is, i1, i2, i3, i4, nspin, ii,jj
         Type (zcsrmat) :: hsoc, op, oph
         Type (zcsrmat) :: zcsr_c, zcsr_temp, zcsr_temp2, zcsr_temp3
         Type (zcsrmat) :: zcsr_r
         Type (zcsrmat) :: zcsr_r2, zcsr_l2
         Type (zcsrmat) :: zcsr_off, zcsr_off2, zcsr_off3, zcsr_null
         Type (zcsrmat) :: zcsr2_temp, zcsr_bigtemp, zcsr_bigtemp2, zcsr_rl2
!!$----------------------------------
         Allocate (ph(scell%size, scell%size), ph_i(scell%size, scell%size), work(scell%size), &
        & s_ando(scell%size))
         nh = 1
         pohops = 1
         if (po%df > 0) then
          n = sc%nrows !- size(geo%to_downfold)
          ndim = sc%nrows
         else
          n = sc%nrows
         endif

         aopts%use_ev_solver = gens
         aopts%dir = - side
         aopts%usestates = frac 
         aopts%refine = 1


         If (scell%need_split /= 0) Then
            aopts%need_split = 1
            aopts%split_idx => scell%idx
            aopts%rotm = scell%rot
         End if

!!$ Prepare k_par and phase factors for supercell
         csz = scell%size
         kbase => scell%kbase
         tr_vec => scell%tr_vec

         Do i = 1, csz
            k1 = kbase (:, i) + kpar
            Do j = 1, csz
               ph (j, i) = Exp (DEF_cmplx_Ione*dot_product(k1, tr_vec(:, j)))
            End Do
         End Do

         ph_i = ph

!!$ Invert the ph_i a matrix
         Call zgetrf (csz, csz, ph_i, csz, ipiv, info)
         If (info == 0) Then
            Call zgetri (csz, ph_i, csz, ipiv, work, csz, info)
         End If
         If (info /= 0) Then
            Write (*,*) 'Problems with phase inversion in get_sc_ando, info =', info
            Stop
         End If

         Nin = 0

         Call init_mathop (sys)

         Allocate (ons(n, n), ofs(n, n), hops(n, 2*n))
         Do i = 1, csz
            
            Call prep_system (sys, sc, kbase(:, i)+kpar, ppar(1), po, 1)

            ons = sptofull (sys%c)
            ofs = sptofull (sys%r)
        
            hops (1:n, 1:n) = ons(1:n, 1:n)
            hops (1:n, n+1:2*n) = ofs(1:n, 1:n)

!!$ this just fix for some nonstandard compilers
            s_ando(i)%alloc = 0

!!$ Solve ando problem
            Call solve_ando (ofs, hops, n, 1, s_ando(i), aopts)

            Nin = Nin + s_ando(i)%Nin

            do is = 1, nspin 
            scell%nreal (i, is) = s_ando(i)%Nin
            end do

         End Do

         Call free_mathop (sys)
         Deallocate (ons, ofs, hops)

!!$ Allocate memory for "ando" structure for supercell
         Call alloc_ando (ando, scell%no, Nin, Nin, aopts)
         ando%emb = DEF_cmplx_zero
         ando%bound = DEF_cmplx_zero
         ando%Uout_i = DEF_cmplx_zero
         ando%Uin = DEF_cmplx_zero
         ando%SNin = 0
         ando%SNout = 0
         If (aopts%need_split == 1) Then
            ando%mask_in = 0
            ando%mask_out = 0
         End If

!!$ Construct ando for supercell
         If (nin > 0) Then

!!$ ! !******$omp parallel private(t_ando,ind1,ind2,jnd1,jnd2,ph1,phi1,j) shared(s_ando,ando,ph,ph_i,csz,sz1)
            ind1 = 0
            ind2 = 0
!!$ ! !******$omp do
            Do i = 1, csz
               If (s_ando(i)%nin == 0) Cycle
               t_ando => s_ando (i)

               If (aopts%need_split /= 0) Then
                  n1 = t_ando%SNin (1)
             
                  If (n1 > 0) Then
                     n0 = ando%SNin (1)
                     ando%mask_in (n0+1:n0+n1:1, 1) = ind2 + t_ando%mask_in(1:n1:1, 1)
                  End If
                  n1 = t_ando%SNin (2)
              
                  If (n1 > 0) Then
                     n0 = ando%SNin (2)
                     ando%mask_in (n0+1:n0+n1:1, 2) = ind2 + t_ando%mask_in(1:n1:1, 2)
                  End If
                  ando%SNin = ando%SNin + t_ando%SNin

                  n1 = t_ando%SNout (1)
               
                  If (n1 > 0) Then
                     n0 = ando%SNout (1)
                     ando%mask_out (n0+1:n0+n1:1, 1) = ind2 + t_ando%mask_out(1:n1:1, 1)
                  End If

                  n1 = t_ando%SNout (2)
                
                  If (n1 > 0) Then
                     n0 = ando%SNout (2)
                     ando%mask_out (n0+1:n0+n1:1, 2) = ind2 + t_ando%mask_out(1:n1:1, 2)
                  End If
                  ando%SNout = ando%SNout + t_ando%SNout
               End If

               ind1 = ind2 + 1
               ind2 = ind2 + t_ando%nin
               jnd1 = 0
               jnd2 = 0
               Do j = 1, csz
                  ph1 = ph (j, i)
                  phi1 = ph_i (i, j)
                  jnd1 = jnd2 + 1
                  jnd2 = jnd2 + n
                  ando%Uin (jnd1:jnd2, ind1:ind2) = t_ando%Uin * ph1
                  ando%bound (jnd1:jnd2, ind1:ind2) = t_ando%bound * ph1
                  ando%Uout_i (ind1:ind2, jnd1:jnd2) = t_ando%Uout_i * phi1
               End Do
               ando%Vin (ind1:ind2) = t_ando%Vin(:)
               ando%Vout (ind1:ind2) = t_ando%Vout(:)
            End Do
!!$ ! !******$omp end do
!!$ ! !******$omp end parallel
         End If
!!$omp parallel do private(jnd1,jnd2,t_emb,emb,j,l,i) firstprivate(csz,n*nspin) shared(ando,ph,ph_i) default(shared)
         Do i = 1, csz
            emb => ando%emb (:, 1+n*(i-1) :n*i) 
            Do j = 1, csz
               t_emb => s_ando(j)%emb
               jnd1 = 0
               jnd2 = 0

               Do l = 1, csz
                  jnd1 = jnd2 + 1
                  jnd2 = jnd2 + n
                  emb (jnd1:jnd2, :) = emb (jnd1:jnd2, :) + ph (l, j) * t_emb * ph_i (j, i)
               End Do
            End Do
         End Do
!!$omp end parallel do
         scl = Sqrt (1.0d0*csz)
         i = scell%no * nin
         Call zdscal (i, 1.0d0/scl, ando%Uin, 1)
         Call zdscal (i, 1.0d0/scl, ando%bound, 1)
         Call zdscal (i, scl, ando%Uout_i, 1)
!!$ free ando structures for small cells
         Do i = 1, csz
            Call free_ando (s_ando(i))
         End Do
         Nin_out = nin !!$/nspin
         Deallocate (ph, ph_i, work, s_ando)
         Nullify (aopts%split_idx)
      End Function get_sc_ando_EMTO


      Subroutine write_ando (ando, kp, s, eq, c)
         Use ando_module
         Implicit None
         Type (t_ando_sollution), Intent (In) :: ando
         Integer, Intent (In) :: kp, s, eq
         Character (Len=*), Intent (In) :: c
!!$          Local vars
         Integer :: pars (14)
         Integer * 8 :: fid
         Integer * 8, External :: c_fopen
         Character (Len=50) :: fn
         Integer :: sz1, sz2, n
         Integer :: csz, rsz
         csz = 2 * DEF_DBL_PREC
         rsz = DEF_DBL_PREC
         pars (1) = ando%n
         pars (2) = ando%nin
         pars (3) = ando%nout
         pars (4) = ando%haveEmb
         pars (5) = ando%haveBound
         pars (6) = ando%haveEvanecent
         pars (7) = ando%haveF
         pars (8) = ando%have_split
         pars (9) = ando%haveAllFU
         pars (10) = ando%dir
         pars (11:12) = ando%SNin(1:2)
         pars (13:14) = ando%SNout(1:2)
         n = ando%n
         Write (fn, '("leads/k",i3.3,"-s",i1.1,"-sg",i3.3,"-",A,".dat")') kp, s, eq, c
         fid = c_fopen (trim(fn), len(trim(fn)), 'w')
         Call c_fwrite (pars, 4, 14, fid)

         If (ando%haveEvanecent == 0) Then
            sz1 = ando%nin
            sz2 = ando%nout
         Else
            sz1 = n
            sz2 = n
         End If

         Call c_fwrite (ando%lin, csz, sz1, fid)
         Call c_fwrite (ando%lout, csz, sz2, fid)
         Call c_fwrite (ando%Uout_i, csz, sz2*n, fid)
         Call c_fwrite (ando%Uin, csz, sz1*n, fid)

         If (ando%haveAllFU /= 0) Then
            Call c_fwrite (ando%Uin_i, csz, sz1*n, fid)
            Call c_fwrite (ando%Uout, csz, sz2*n, fid)
         End If

         If (ando%haveF /= 0) Then
            Call c_fwrite (ando%Fout_i, csz, (n*n), fid)
            Call c_fwrite (ando%Fin, csz, (n*n), fid)
            If (ando%haveAllFU /= 0) Then
               Call c_fwrite (ando%Fin_i, csz, (n*n), fid)
               Call c_fwrite (ando%Fout, csz, (n*n), fid)
            End If
         End If
         If (ando%haveBound /= 0) Call c_fwrite (ando%bound, csz, (n*ando%nin), fid)

         If (ando%haveEmb /= 0) Call c_fwrite (ando%emb, csz, (n*n), fid)
         Call c_fwrite (ando%Vin, rsz, (ando%nin), fid)
         Call c_fwrite (ando%Vout, rsz, (ando%nout), fid)
         If (ando%have_split /= 0) Then
            Call c_fwrite (ando%mask_in, 4, (ando%nin*2), fid)
            Call c_fwrite (ando%mask_out, 4, (ando%nout*2), fid)
         End If
         Call c_fclose (fid)

      End Subroutine write_ando

      Function read_ando (ando, kp, s, eq, c) Result (nin)
         Use ando_module
         Implicit None
         Type (t_ando_sollution), Intent (Inout) :: ando
         Integer, Intent (In) :: kp, s, eq
         Character (Len=*), Intent (In) :: c
!!$          Local vars
         Integer :: pars (14)
         Integer * 8 :: fid
         Integer * 8, External :: c_fopen
         Character (Len=50) :: fn
         Integer :: sz1, sz2, n, nin, nout
         Integer :: csz, rsz
         Type (t_ando_options) :: aopts

         csz = 2 * DEF_DBL_PREC
         rsz = DEF_DBL_PREC
         Write (fn, '("leads/k",i3.3,"-s",i1.1,"-sg",i3.3,"-",A,".dat")') kp, s, eq, c
         fid = c_fopen (trim(fn), len(trim(fn)), 'r')
         Call c_fread (pars, 4, 14, fid)

         n = pars (1)
         nin = pars (2)
         nout = pars (3)
         aopts%needEmb = pars (4)
         aopts%needBound = pars (5)
         aopts%needEvanecent = pars (6)
         aopts%needF = pars (7)
         aopts%need_split = pars (8)
         aopts%needAllFU = pars (9)
         aopts%dir = pars (10)
!!$          n = ando%n
!!$          write(*,*) nin, nout

         Call alloc_ando (ando, n, nin, nout, aopts)

         ando%SNin (1:2) = pars (11:12)
         ando%SNout (1:2) = pars (13:14)
         If (ando%haveEvanecent == 0) Then
            sz1 = ando%nin
            sz2 = ando%nout
         Else
            sz1 = n
            sz2 = n
         End If

         Call c_fread (ando%lin, csz, sz1, fid)
         Call c_fread (ando%lout, csz, sz2, fid)
         Call c_fread (ando%Uout_i, csz, sz2*n, fid)
         Call c_fread (ando%Uin, csz, sz1*n, fid)

         If (ando%haveAllFU /= 0) Then
            Call c_fread (ando%Uin_i, csz, sz1*n, fid)
            Call c_fread (ando%Uout, csz, sz2*n, fid)
         End If

         If (ando%haveF /= 0) Then
            Call c_fread (ando%Fout_i, csz, (n*n), fid)
            Call c_fread (ando%Fin, csz, (n*n), fid)
            If (ando%haveAllFU /= 0) Then
               Call c_fread (ando%Fin_i, csz, (n*n), fid)
               Call c_fread (ando%Fout, csz, (n*n), fid)
            End If
         End If
         If (ando%haveBound /= 0) Call c_fread (ando%bound, csz, (n*ando%nin), fid)

         If (ando%haveEmb /= 0) Call c_fread (ando%emb, csz, (n*n), fid)
         Call c_fread (ando%Vin, rsz, (ando%nin), fid)
         Call c_fread (ando%Vout, rsz, (ando%nout), fid)
         If (ando%have_split /= 0) Then
            Call c_fread (ando%mask_in, 4, (ando%nin*2), fid)
            Call c_fread (ando%mask_out, 4, (ando%nout*2), fid)
         End If
         Call c_fclose (fid)

      End Function read_ando
End Module supercell
