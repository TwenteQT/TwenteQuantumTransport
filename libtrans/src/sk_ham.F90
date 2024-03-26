#include "math_def.h"

Module sk_ham
      Use atoms_module
      Use structure


Contains

      Subroutine make_onsite (spin, atclass, hblock, oblock)
         Implicit None
         Type (t_atom_defenition), Intent (In) :: atclass
         Integer, Intent (In) :: spin
         Real (Kind=DEF_DBL_PREC) :: hblock (:, :), oblock (:, :)
!!$      real(kind=DEF_DBL_PREC) :: hblock(4,4), oblock(4,4)
!!$ Local
         Integer :: ioff

         hblock (1, 1) = atclass%s(spin)
         oblock (1, 1) = 1.0d0

         If (atclass%nl == 1) Return

         hblock (2, 2) = atclass%p(spin)
         hblock (3, 3) = atclass%p(spin)
         hblock (4, 4) = atclass%p(spin)

         oblock (2, 2) = 1.0d0
         oblock (3, 3) = 1.0d0
         oblock (4, 4) = 1.0d0

         If (atclass%nl == 2) Return

         hblock (5, 5) = atclass%d1(spin)
         hblock (6, 6) = atclass%d1(spin)
         hblock (8, 8) = atclass%d1(spin)

         hblock (7, 7) = atclass%d2(spin)
         hblock (9, 9) = atclass%d2(spin)

         oblock (5, 5) = 1.0d0
         oblock (6, 6) = 1.0d0
         oblock (8, 8) = 1.0d0
         oblock (7, 7) = 1.0d0
         oblock (9, 9) = 1.0d0

      End Subroutine make_onsite

      Subroutine slater_ham_ij (iclass, jclass, diff, spin, ortho, skh, hblock, sko, oblock)
!!$ Calculates blocks of Hamiltonian(H) and Overlap (O) equation, when ortho==F the oblock is
!! returned as either diagonal (for i=j) or zero (for i/=j)  matrix
         Use spherical
         Implicit None
         Type (t_atom_defenition) :: iclass, jclass
         Real (Kind=DEF_DBL_PREC) :: diff (3)
         Integer :: spin
         Type (t_sk_param) :: skh
         Type (t_sk_param) :: sko
         Real (Kind=DEF_DBL_PREC) :: hblock (:, :), oblock (:, :)
         Logical :: ortho
!!$       Local
         Integer :: ishell, ilmx, jlmx, lspin
         Real (Kind=DEF_DBL_PREC) :: dist, rot (3, 3), eul (3), dist2d
         Real (Kind=DEF_DBL_PREC), Pointer :: rotA (:, :), rotB (:, :)

         hblock = 0.0d0
         oblock = 0.0d0

         dist = Sqrt (sum(diff**2))
!!$      dist2d=sqrt(diff(1)**2+diff(3)**2)

         If ( .Not. which_shell(dist, skh, ishell)) Return

         lspin = Min (spin, skh%ns)

         If (ishell == 0) Then
            Call make_onsite (lspin, iclass, hblock, oblock)
         Else
            ilmx = iclass%nl - 1
            jlmx = jclass%nl - 1

            Allocate (rotA((ilmx+1)**2, (ilmx+1)**2))
            Allocate (rotB((jlmx+1)**2, (jlmx+1)**2))

!!$  Defining H and O for system of coordinates with z||Ri-Rj
            Call def_skblock (iclass, jclass, dist, ishell, lspin, skh, hblock, 'and')
            If ( .Not. ortho) Call def_skblock (iclass, jclass, dist, ishell, lspin, sko, oblock, 'and')

!!$ Rotation matrix to the global s.o.c
            Call rot_mat (diff, rot)
            eul = euler_ang (rot)

!!$ Wigner coefficients for A and B atoms ( can differ in rank)
            Call wigner_mat (ilmx, eul, rotA)
            Call wigner_mat (jlmx, eul, rotB)
            rotA = transpose (rotA)

!!$  Rotating Hamiltonian back to global s.o.c,  H'=DA^T*H*DB
            hblock = matmul (rotA, matmul(hblock, rotB))
            If ( .Not. ortho) oblock = matmul (rotA, matmul(oblock, rotB))

            Deallocate (rotA)
            Deallocate (rotB)

         End If
         
         hblock = hblock + oblock * (iclass%sk_shift(min(lspin,iclass%ns))+jclass%sk_shift(min(lspin,jclass%ns))) / 2

      End Subroutine slater_ham_ij

      Function scale_sk (hop, dist, rshell, l, l1, m, scaling)
         Implicit None
         Real (Kind=DEF_DBL_PREC) :: scale_sk, hop, dist, rshell
         Integer :: l, l1, m
         Character (Len=3) :: scaling

         Select Case (trim(scaling))
         Case ('and')
            scale_sk = hop * (dist/rshell) ** (-l-l1-1)
         Case ('sin')
            scale_sk = hop * (dist/rshell) ** (-3)
         Case Default
            scale_sk = hop
         End Select

      End Function scale_sk

      Subroutine def_skblock (iclass, jclass, dist, ishell, spin, sk, hblock, scaling)
         Implicit None
         Integer, Intent (In) :: spin, ishell
         Type (t_atom_defenition) :: iclass, jclass
         Real (Kind=DEF_DBL_PREC) :: dist
         Type (t_sk_param) :: sk
         Real (Kind=DEF_DBL_PREC) :: hblock (:, :)
!!$       real(kind=DEF_DBL_PREC) :: hblock(4,4)
         Character (Len=3) :: scaling
!!$ Local
!!$      character(len=1) :: mode
!!$      Local
!!$      real(kind=DEF_DBL_PREC), pointer :: hop(:,:,:)
         Real (Kind=DEF_DBL_PREC) :: rshell
         Character (Len=3) :: mode
         Integer :: ilmx, jlmx

         rshell = sk%rshell (ishell)
         ilmx = iclass%nl - 1
         jlmx = jclass%nl - 1

         hblock = 0.0d0
         If (sk%nohopping .Eqv. .True.) Return

         hblock (1, 1) = scale_sk (sk%ss_s(ishell, spin), dist, rshell, 0, 0, 0, scaling)

         If (jlmx > 0) Then
            hblock (1, 3) = scale_sk (sk%sp_s(ishell, spin), dist, rshell, 0, 1, 0, scaling)
         End If

         If (ilmx > 0) Then
            hblock (3, 1) = scale_sk (sk%ps_s(ishell, spin), dist, rshell, 1, 0, 0, scaling)
         End If

         If (jlmx > 1) Then
            hblock (1, 7) = scale_sk (sk%sd_s(ishell, spin), dist, rshell, 0, 2, 0, scaling)
         End If

         If (ilmx > 1) Then
            hblock (7, 1) = scale_sk (sk%ds_s(ishell, spin), dist, rshell, 2, 0, 0, scaling)
         End If

         If ((ilmx > 0) .And. (jlmx > 0)) Then
            hblock (2, 2) = scale_sk (sk%pp_p(ishell, spin), dist, rshell, 1, 1,-1, scaling)
            hblock (3, 3) = scale_sk (sk%pp_s(ishell, spin), dist, rshell, 1, 1, 0, scaling)
            hblock (4, 4) = scale_sk (sk%pp_p(ishell, spin), dist, rshell, 1, 1, 1, scaling)
         End If

         If ((ilmx > 0) .And. (jlmx > 1)) Then
            hblock (2, 6) = scale_sk (sk%pd_p(ishell, spin), dist, rshell, 1, 2,-1, scaling)
            hblock (3, 7) = scale_sk (sk%pd_s(ishell, spin), dist, rshell, 1, 2, 0, scaling)
            hblock (4, 8) = scale_sk (sk%pd_p(ishell, spin), dist, rshell, 1, 2, 1, scaling)
         End If

         If ((ilmx > 1) .And. (jlmx > 0)) Then
            hblock (6, 2) = scale_sk (sk%dp_p(ishell, spin), dist, rshell, 2, 1,-1, scaling)
            hblock (7, 3) = scale_sk (sk%dp_s(ishell, spin), dist, rshell, 2, 1, 0, scaling)
            hblock (8, 4) = scale_sk (sk%dp_p(ishell, spin), dist, rshell, 2, 1, 1, scaling)
         End If

         If ((ilmx > 1) .And. (jlmx > 1)) Then
            hblock (5, 5) = scale_sk (sk%dd_d(ishell, spin), dist, rshell, 2, 2,-2, scaling)
            hblock (6, 6) = scale_sk (sk%dd_p(ishell, spin), dist, rshell, 2, 2,-1, scaling)
            hblock (7, 7) = scale_sk (sk%dd_s(ishell, spin), dist, rshell, 2, 2, 0, scaling)
            hblock (8, 8) = scale_sk (sk%dd_p(ishell, spin), dist, rshell, 2, 2, 1, scaling)
            hblock (9, 9) = scale_sk (sk%dd_d(ishell, spin), dist, rshell, 2, 2, 2, scaling)
         End If

      End Subroutine def_skblock

      Function which_shell (dist, sk, ishell)
         Implicit None
         Logical :: which_shell
         Real (Kind=DEF_DBL_PREC), Intent (In) :: dist
         Type (t_sk_param), Intent (In) :: sk
         Integer, Intent (Out) :: ishell
!!$ Local
         Integer :: nshell, i
         Real (Kind=DEF_DBL_PREC) :: lim_shell

         nshell = sk%nshell

         which_shell = .False.

         If (dist <= 1.0d-3) Then
            ishell = 0
            which_shell = .True.
            Return
         End If

         Do i = 1, nshell
            lim_shell = (sk%rshell(i)+sk%rshell(i+1)) / 2
            If (dist <= lim_shell) Then
               ishell = i
               which_shell = .True.
               Return
            End If
         End Do

      End Function which_shell


      Subroutine wigner_mat (lmx, eul, mat)
         Use spherical
         Implicit None
         Integer :: lmx
         Real (Kind=DEF_DBL_PREC) :: eul (3)
         Real (Kind=DEF_DBL_PREC), Pointer :: mat (:, :)
!!$      Local
         Integer :: lm1, lm2, m1, m2

         mat = 0.0d0

         mat (1, 1) = 1

         If (lmx == 0) Return

         Do m1 = - 1, 1
            Do m2 = - 1, 1
               mat (3+m1, 3+m2) = rwigner (1, m1, m2, eul(1), eul(2), eul(3))
            End Do
         End Do

         If (lmx == 1) Return

         Do m1 = - 2, 2
            Do m2 = - 2, 2
               mat (7+m1, 7+m2) = rwigner (2, m1, m2, eul(1), eul(2), eul(3))
            End Do
         End Do

      End Subroutine wigner_mat


      Function calc_slater_ham (geom, atoms, ene, ortho, sk_2d) Result (scalpha)
         Use geometry_module
         Use structure
         Use logging
         Implicit None
         Type (t_geometry) :: geom
         Type (t_atoms_set) :: atoms
         Type (t_strconst), Pointer :: scalpha (:)
         Type (t_strconst), Pointer :: scr
!!$       Type (t_strconst) :: scalpha
         Real (Kind=DEF_DBL_PREC) :: ene
         Logical :: ortho !, sk_2d
         Logical, Optional :: sk_2d
!!$ local
         Type (t_neighbour_list), Pointer :: neibs (:)
         Type (t_cluster_list) :: clslist
         Integer :: numat, lirel, i, j, isp

         Allocate (scalpha(1:atoms%nsmax))

!!$ find neighbours for each site, the blocks of H-EO are attached
!! to the elements of the neighbour list
         Call do_log (1, '   Calculating neighbours list')

         neibs => find_slater_neibs (geom, atoms, sk_2d)
         Call rep_neighbours (neibs, geom)

         Call do_log (1, '   Preparing Hamiltonian:')
         Do isp = 1, atoms%nsmax
            scr => scalpha (isp)

            Call prep_ham_slater (geom, neibs, atoms, isp, scr, ortho, ene)

         End Do
         Call do_log (1, '   Done')
!!$ deallocate memory
         Call free_neibs (neibs)
         Deallocate (neibs)
         Call do_log (1, '')

         Nullify (scr)

      End Function calc_slater_ham


      Function find_slater_neibs (geom, atoms, sk_2d) Result (neibs)
!!$ this subroutine find neighbours on distance less that "cutrat"
         Use geometry_module
         Implicit None
         Type (t_atoms_set) :: atoms
         Type (t_geometry) :: geom
         Type (t_neighbour_list), Pointer :: neibs (:)
         Logical, Optional :: sk_2d
!!$ local
         Integer, Parameter :: initial_maxneib = 50
         Integer :: i, jx, jy, nat, k, maxneib, t, tb, te, skip
         Real (Kind=DEF_DBL_PREC) :: dist, base (2, 2), cutrat, wsr1, mdist, dist2d
         Real (Kind=DEF_DBL_PREC) :: trans (3,-1:1), scale, dir (3)
         Real (Kind=DEF_DBL_PREC), Pointer :: coord (:), coord0 (:)
         Type (t_atom_pointer), Pointer :: atpi, atpk
         Type (t_neighbour), Pointer :: t_neib
         Integer :: icls, kcls, ishell, nspin, isp, ilmx, klmx, iorb, korb, ntr1, ntr2
         Type (t_sk_param), Pointer :: skh

         ntr1 = geom%ntrpar
         ntr2 = geom%ntrpar
         If (present(sk_2d)) Then
            If (sk_2d) ntr2 = 0
         End If

         nspin = atoms%nsmax
         scale = 1.0d0

         maxneib = initial_maxneib
         nat = geom%num

!!$ allocate memory for clusters
         Allocate (neibs(nat))

!!$ put atoms itself as a first neighbours

         Allocate (coord(3), coord0(3))
         coord (:) = 0.0d0

!!$ !$omp parallel default(shared) firstprivate(scale,nat)
!$omp parallel do private(i) default(shared)
         Do i = 1, nat
            Call add_neib (neibs(i), geom%atoms(i), i, coord, 0, 0, 0)
            neibs(i)%basenl = geom%atoms(i)%ptr%nl
         End Do
!$omp end parallel do

         base = geom%base
         trans = 0.0d0
         trans (:,-1) = - geom%l_perp_trans
         trans (:, 1) = geom%r_perp_trans


!!$ loop over atoms.
         Do i = 1, nat

            atpi => geom%atoms (i)
            icls = atpi%ptr%index

            tb = 0
            te = 0

            If (i <= geom%l_transnum) Then
               tb = - 1
            End If
            If (i > (nat-geom%r_transnum)) Then
               te = 1
            End If

            Do t = tb, te
               coord0 (:) = atpi%coord + trans (:, t)
               coord (:) = coord0 (:)
!!$ loops over parallel translation indexes (jx,jy)

               Do jy = - ntr2, ntr2
                  Do jx = - ntr1, ntr1
                     skip = Abs (jx) + Abs (jy) + Abs (t)
                     coord (1:2) = coord0 (1:2) + jx * base (:, 1) + jy * base (:, 2)
!!$ loop over atoms again, if atom wchich defined by three previous loops-indexes
!!$ is in neighbourhood, then we will add it.

!$omp parallel default(shared)&
!$omp& private(atpk,t_neib,dist,dir,mdist,k)

!$omp do ordered
                     Do k = 1, nat
                        If (skip == 0 .And. k == i) Cycle
                        atpk => geom%atoms (k)
                        kcls = atpk%ptr%index

                        skh => atoms%skh (icls, kcls)

                        dist = get_dist (atpk%coord, coord, dir)
!!$                     dist2d=sqrt(dir(1)**2+dir(3)**2)

                        If (which_shell(dist, skh, ishell)) Then
                           Call add_neib (neibs(k), atpi, i, dir, t, jx, jy)
                        End If


                     End Do
!$omp end do

!$omp end parallel
                  End Do
               End Do

            End Do
         End Do
         Deallocate (coord, coord0)


!!$ convert list to array and calculate blocks of H-EO
         Do i = 1, nat
            Allocate (neibs(i)%nbr(neibs(i)%num))
            Allocate (neibs(i)%cl_dir_idx(neibs(i)%num))
            neibs(i)%nbr(1) = neibs(i)%first
            Do k = 1, neibs(i)%num
               If (k /= 1) neibs(i)%nbr(k) = neibs(i)%nbr(k-1)%next
            End Do
         End Do

         Return
      Contains
         Subroutine add_neib (neibs, atp, n, dir, t, ix, iy)
            Implicit None
            Type (t_atom_pointer) :: atp
            Type (t_neighbour_list) :: neibs
            Integer :: n, t, ix, iy
            Real (Kind=DEF_DBL_PREC) :: dir (3), scale
!!$ local
            Type (t_neighbour), Pointer :: t_neib

            Allocate (t_neib)
            t_neib%nl = atp%ptr%nl
            t_neib%atnum = n
            t_neib%clsnum = atp%ptr%index
            t_neib%nrm = 0
            t_neib%alpha => atp%ptr%alpha
            t_neib%rvec = dir
            t_neib%trpar = (/ ix, iy /)
            t_neib%trperp = t

!!$ !$omp ordered
            If (neibs%num /= 0) Then
               neibs%last%next => t_neib
            Else
               neibs%first => t_neib
            End If

            neibs%last => t_neib
            neibs%num = neibs%num + 1
!!$ !$omp end ordered
            Nullify (t_neib)
         End Subroutine add_neib

         Function get_dist (a, b, r) Result (dist)
            Implicit None
            Real (Kind=DEF_DBL_PREC) :: a (3), b (3), r (3), dist
            r = b - a
            dist = Sqrt (sum(r*r))
            Return
         End Function get_dist
      End Function find_slater_neibs



      Subroutine rep_neighbours (neibs, geom)
         Use geometry_module
         Use logging
         Implicit None
         Type (t_neighbour_list), Pointer :: neibs (:)
         Type (t_geometry) :: geom

!!$ Local
         Integer :: nat, iat, nneighb, ineighb, repf, icls
         Character (Len=1024) :: cwork, cwork1

         repf = log_stream_number

         Do iat = 1, geom%num
            Write (cwork, '("     Atom ",i5,": Cluster of ",i3," sites")') iat, neibs(iat)%num
            Call do_log (1, trim(cwork))

!!$         if(Log_Level<6) cycle

            Do ineighb = 1, neibs(iat)%num
               icls = neibs(iat)%nbr(ineighb)%atnum
               Write (cwork, '("     ",i5,". AT=",a4," DIR=",3(f7.4,1x))') ineighb, &
              & geom%atoms(icls)%ptr%label, neibs(iat)%nbr(ineighb)%rvec
               Call do_log (6, trim(cwork))
            End Do

         End Do

      End Subroutine rep_neighbours


      Subroutine scale_geom (geo, alat)
         Use geometry_module
         Implicit None
         Type (t_geometry) :: geo
         Real (Kind=DEF_DBL_PREC) :: alat
!!$ local
         Integer :: iat

         geo%base = geo%base * alat
         geo%perp_trans = geo%perp_trans * alat
         geo%l_perp_trans = geo%l_perp_trans * alat
         geo%r_perp_trans = geo%r_perp_trans * alat

         Do iat = 1, geo%num
            geo%atoms(iat)%coord = geo%atoms(iat)%coord * alat
         End Do

      End Subroutine scale_geom


      Subroutine prep_ham_slater (geom, neibs, atoms, isp, scr, ortho, eoff)
!!$ Finally construct datatype "t_strconst" which holds all information which
!!$ we need to calculate S(k) for region and hopings to outside
         Use geometry_module
         Use logging
         Implicit None

         Type (t_neighbour_list), Target :: neibs (:)
!!$         Type (t_cluster_list) :: clslist
         Type (t_strconst), Pointer :: scr
         Type (t_geometry) :: geom
         Type (t_atoms_set) :: atoms
         Integer, Intent (In) :: isp
         Logical :: ortho
         Real (Kind=DEF_DBL_PREC) :: eoff
!!$ local
         Integer :: nat, lnat, rnat, il, ir, atcls
         Integer :: i, j, norb, rowidx (atoms%nlmax**2), startrow (geom%num)
         Integer :: nlhop, nrhop, numb (3, geom%num), lnumb (3, geom%num), rnumb (3, geom%num)
         Type (t_neighbour_list), Pointer :: t_neib
         Type (t_str_site), Pointer :: t_ss
         Real (Kind=DEF_DBL_PREC) :: ldev, dev, rdev
         Character (Len=100) :: cwork

         nat = geom%num
         lnat = geom%l_transnum
         rnat = geom%r_transnum

         scr%nm = nat
         scr%base = geom%base
         Do i = 1, atoms%nlmax ** 2
            rowidx (i) = i
         End Do
         nlhop = 0
         nrhop = 0
         norb = 0
         Do i = 1, nat
            startrow (i) = norb
            t_neib => neibs (i)
            Call sort_neibs_num (t_neib)
            Call calc_numbers (t_neib, numb(:, i), lnumb(:, i), rnumb(:, i))

            If (lnumb(1, i) > 0) nlhop = nlhop + 1
            If (rnumb(1, i) > 0) nrhop = nrhop + 1
            norb = norb + t_neib%basenl ** 2

         End Do

         If (lnat < nlhop .Or. rnat < nrhop) Then
            Write (*,*) 'Definetelly something wrong with defenitions left-,right- edges'
            Write (*,*) 'Num.atoms with left-,right- hopping is ', nlhop, nrhop
            Write (*,*) 'Size of left-,right- layers is', lnat, rnat
            Stop
         End If

         If (scr%alloc /= 1) Then
            Allocate (scr%main(nat))
            Allocate (scr%lhop(nlhop))
            Allocate (scr%rhop(nrhop))
            scr%alloc = 1
         End If

         scr%nl = nlhop
         scr%nr = nrhop

         scr%nrows = norb
         scr%lnrows = norb
         scr%rnrows = norb
         scr%nnz = 0
         il = 0
         ir = 0
         ldev = 0.0d0
         rdev = 0.0d0
         dev = 0.0d0

         Do i = 1, nat
            atcls = geom%atoms(i)%ptr%index
            t_ss => scr%main (i)

            scr%nnz = scr%nnz + fillstr_slater (0, t_ss, neibs, numb(:, i), rowidx, startrow, i, dev, isp, &
           & atoms, ortho, atcls)

            If (lnumb(1, i) > 0) Then

               il = il + 1
               t_ss => scr%lhop (il)
               scr%lnnz = scr%lnnz + fillstr_slater (-1, t_ss, neibs, lnumb(:, i), rowidx, startrow, i, ldev, &
              & isp, atoms, ortho, atcls)

            End If

            If (rnumb(1, i) > 0) Then

               ir = ir + 1
               t_ss => scr%rhop (ir)
               scr%rnnz = scr%rnnz + fillstr_slater (1, t_ss, neibs, rnumb(:, i), rowidx, startrow, i, rdev, &
              & isp, atoms, ortho, atcls)

            End If

         End Do

         Do i = 1, nat
            Do j = 1, neibs(i)%num
               neibs(i)%nbr(j)%nrm = 0
               Nullify (neibs(i)%nbr(j)%sr)
            End Do
         End Do

      Contains

#define E_SAME_COORD 1.0d-8

         Function fillstr_slater (t, t_ss, neibs, numb, rowidx, startrow, sn, dev, isp, atoms, ortho, icls) &
        & Result (nnz)
            Use logging
            Implicit None
            Type (t_neighbour_list), Target :: neibs (:)
            Type (t_str_site), Target :: t_ss
            Integer :: numb (3), t, rowidx (:), startrow (:), nnz, sn
            Integer, Intent (In) :: isp, icls
            Type (t_atoms_set) :: atoms
            Logical :: ortho
!!$ Local
            Type (t_atom_defenition), Pointer :: iatom, jatom
            Type (t_neighbour_list), Pointer :: t_neib
            Integer :: i, inorb, atn, nhop, jnorb, dir, nnb, st, st1, j, jat, jdir, iat, jcls
            Type (t_neighbour), Pointer :: t_nbr, t_jnbr
            Real (Kind=DEF_DBL_PREC) :: rvec (3), dev, ef, pshift
            Real (Kind=DEF_DBL_PREC), Pointer :: hblock (:, :), oblock (:, :)
            Type (t_sk_param), Pointer :: skh, sko


            t_neib => neibs (sn)

            inorb = atoms%nlmax ** 2
            Allocate (hblock(inorb, inorb))
            Allocate (oblock(inorb, inorb))

            inorb = (t_neib%basenl) ** 2

            iatom => atoms%at (icls)
            ef = iatom%ef + eoff

            Call allocstr (t_ss, numb)

            t_ss%nrows = inorb
            atn = 0
            nhop = 0
            nnb = 0
            st1 = 1
            t_ss%srow = startrow (sn) + 1
            t_ss%ia = sn

            Do i = 1, t_neib%num
               t_nbr => t_neib%nbr (i)

               If (t_nbr%trperp == t) Then
                  nhop = nhop + 1
                  jnorb = (t_nbr%nl) ** 2 !t_cls%site(dir)%nl ** 2
                  jat = t_nbr%atnum

                  If (jat > atn) Then
                     atn = jat
                     nnb = nnb + 1
                     st = st1
                     t_ss%rowind (st:st+jnorb-1) = rowidx (1:jnorb) + startrow (atn)
                     st1 = st1 + jnorb
                  End If

                  Allocate (t_ss%hop(nhop)%bl(inorb, jnorb))

!!$ Find symmetric element
                  If (t_nbr%nrm == 0) Then

                     jcls = t_nbr%clsnum
                     jatom => atoms%at (jcls)

                     skh => atoms%skh (icls, jcls)
                     If ( .Not. ortho) sko => atoms%sko(icls, jcls)

                     rvec = t_nbr%rvec

                     hblock = 0.0d0
                     oblock = 0.0d0

!!$ Calculating H(i,j) and O(i,j) blocks
                     Call slater_ham_ij (iatom, jatom, rvec, isp, ortho, skh, hblock(1:inorb, 1:jnorb), sko, &
                    & oblock(1:inorb, 1:jnorb))
!!$ Putting together equation of motion: H-EO
                     t_ss%hop(nhop)%bl = hblock (1:inorb, 1:jnorb) - ef * oblock (1:inorb, 1:jnorb)

                     Do j = 1, neibs(jat)%num ! Symmetrization loop
                        t_jnbr => neibs(jat)%nbr(j)

                        If (sum(Abs(rvec+t_jnbr%rvec)) < E_SAME_COORD) Then

!!$                           if (t_jnbr%clsnum/=icls) then
!!$                              call do_log(1,'Problems with symmetrization in fillstr_slater')
!!$                              write(*,*) t_jnbr%clsnum,icls,t,t_jnbr%trperp
!!$                              stop
!!$                           end if

                           skh => atoms%skh (jcls, icls)
                           If ( .Not. ortho) sko => atoms%sko(jcls, icls)

                           rvec = t_jnbr%rvec
                           Call slater_ham_ij (jatom, iatom, rvec, isp, ortho, skh, hblock(1:jnorb, 1:inorb), &
                          & sko, oblock(1:jnorb, 1:inorb))

                           hblock (1:jnorb, 1:inorb) = hblock (1:jnorb, 1:inorb) - ef * oblock (1:jnorb, &
                          & 1:inorb)

                           dev = Max (dev, maxval(Abs(t_ss%hop(nhop)%bl(:, :)-transpose(hblock(1:jnorb, &
                          & 1:inorb)))))
                           t_ss%hop(nhop)%bl = 0.5d0 * (t_ss%hop(nhop)%bl+transpose(hblock(1:jnorb, &
                          & 1:inorb)))

                           If (t == 0) Then
                              neibs(jat)%nbr(j)%nrm = 1
                              neibs(jat)%nbr(j)%sr => t_ss%hop(nhop)%bl
                           End If

                           Exit
                        End If

                     End Do

                  Else
                     t_ss%hop(nhop)%bl = transpose (t_nbr%sr)
                  End If

                  t_ss%hop(nhop)%ncol = jnorb
                  t_ss%hop(nhop)%nrow = norb
                  t_ss%trpar (:, nhop) = t_nbr%trpar
                  t_ss%startn (nhop) = st
               End If
            End Do

            nnz = inorb * numb (3)

            Deallocate (hblock)
            Deallocate (oblock)

            Nullify (skh)
            Nullify (sko)

            Nullify (iatom)
            Nullify (jatom)

            Return
         End Function fillstr_slater

         Subroutine allocstr (t_ss, numb)
            Implicit None
            Type (t_str_site) :: t_ss
            Integer :: numb (3)
            Integer :: nhop, nnb, norb

            t_ss%nhop = numb (2)
            t_ss%nnb = numb (1)
            t_ss%ncols = numb (3)
            If (t_ss%alloc /= 1) Then
               t_ss%alloc = 1
               Allocate (t_ss%hop(t_ss%nhop))
               Allocate (t_ss%startn(t_ss%nhop))
               Allocate (t_ss%trpar(2, t_ss%nhop))
               Allocate (t_ss%rowind(t_ss%ncols))
            End If
         End Subroutine allocstr


         Subroutine calc_numbers (t_neib, numb, lnumb, rnumb)
            Implicit None
            Type (t_neighbour_list) :: t_neib
            Integer :: numb (3), lnumb (3), rnumb (3)
!!$ Local
            Type (t_neighbour), Pointer :: t_nbr
            Integer :: i, latn, atn, ratn, atnum, norb

            numb = 0
            lnumb = 0
            rnumb = 0
            latn = 0
            atn = 0
            ratn = 0
            Do i = 1, t_neib%num
               t_nbr => t_neib%nbr (i)
               atnum = t_nbr%atnum
               norb = t_nbr%nl ** 2
               If (t_nbr%trperp == 0) Then
                  numb (2) = numb (2) + 1
                  If (atnum > atn) Then
                     numb (1) = numb (1) + 1
                     atn = atnum
                     numb (3) = numb (3) + norb
                  End If
               Else
                  If (t_nbr%trperp < 0) Then
                     lnumb (2) = lnumb (2) + 1
                     If (atnum > latn) Then
                        lnumb (1) = lnumb (1) + 1
                        latn = atnum
                        lnumb (3) = lnumb (3) + norb
                     End If
                  End If
                  If (t_nbr%trperp > 0) Then
                     rnumb (2) = rnumb (2) + 1
                     If (atnum > ratn) Then
                        rnumb (1) = rnumb (1) + 1
                        ratn = atnum
                        rnumb (3) = rnumb (3) + norb
                     End If
                  End If
               End If
            End Do

         End Subroutine calc_numbers
      End Subroutine prep_ham_slater

      Subroutine gibz_2d (bz, vbr, bzopt)
         Use logging
         Use bzgrid
!!$****************************************************************
!!$  GENERATES NETWORK OF K||-POINTS IN IRREDUCIBLE BRILLOUIN ZONE
!!$****************************************************************
         Implicit None
         Include 'mpif.h'
         Real (Kind=DEF_DBL_PREC) :: vbr (2, 2)
         Type (bzone) :: bz
         Type (t_bz_opts), Intent (In) :: bzopt

!!$ Local
         Real (Kind=DEF_DBL_PREC) :: scale, bzdef (4), kpoint (2)
         Integer :: nk, nsym, inve, iftr

         Character (Len=30) :: cwork
         Integer :: istart !,  ext=0 nbk1 (2*(2*nk+1)**2), nbk2 (2*(2*nk+1)**2),
         Real (Kind=DEF_DBL_PREC) :: vbg (2, 2)!, k_pos(2) akbz (2, 2*(2*nk+1)**2), wkbz (2*(2*nk+1)**2),
         Real (Kind=DEF_DBL_PREC) :: romega, gomega, ajm, rnk, anu3, dky, dkx, q, dk, lf, bt, hg, wd, dx, dy
         Real (Kind=DEF_DBL_PREC) :: p, fx2, fx1, tw, cnorm, det, stx, sty, pofs (2)
         Integer :: nbz, iy, nu, ndy, ndx, nmy, nmx, i, j, i2, i1, nk1, nk2, ibz, mnbz


         Call do_log (1, ' Generating k mesh ...')

         scale = bzopt%BZpart
         nk = bzopt%nk
         nsym = bzopt%nsym
         inve = bzopt%inve
         iftr = bzopt%iftr
         bzdef = bzopt%bzone
         kpoint = bzopt%kpoint


         mnbz = (2*nk+1) ** 2

         det = vbr (1, 1) * vbr (2, 2) - vbr (1, 2) * vbr (2, 1)
         romega = Abs (det)
         gomega = (2.0d0*DEF_M_PI) ** 2 / romega
         cnorm = 2.0d0 * DEF_M_PI / det

         vbg (1, 1) = cnorm * vbr (2, 2)
         vbg (2, 1) = - cnorm * vbr (1, 2)
         vbg (1, 2) = - cnorm * vbr (2, 1)
         vbg (2, 2) = cnorm * vbr (1, 1)


         Write (cwork, '(g15.7)') romega
         Call do_log (1, '   AREA OF 2D-PRIMITIVE CELL: REAL= '//trim(cwork))
         Write (cwork, '(g15.7)') gomega
         Call do_log (1, '                        RECIPROCAL= '//trim(cwork))
         Write (cwork, '(2 g15.7)') vbg (:, 1)
         Call do_log (1, '   RECIP. BASIS: 1. VECTOR = '//trim(cwork))
         Write (cwork, '(2 g15.7)') vbg (:, 2)
         Call do_log (1, '                 2. VECTOR = '//trim(cwork))


         If (nsym .Eq.-10) Go To 1188


!!$-------------------------- NSYM=0: GENERAL CASE

!!$ SHIFT THE K|| MESH A LITTER BIT IN ORDER TO AVOID
!!$ THE KX=0 OR KY=0


         ibz = 0
         tw = 0.0d0
         nk1 = nk
         nk2 = 0

         If (iftr == 0) Then
            istart = - nk1 + 1
            nbz = 2 * nk
         Else
            istart = 1
            nbz = nk
         End If

         Allocate (bz%k(2, nbz))
         Allocate (bz%ik(2, nbz))
         Allocate (bz%kweight(nbz))

         stx = real ((2*nk1), kind=DEF_DBL_PREC)

         Do i1 = istart, nk1
            fx1 = (real(i1, kind=DEF_DBL_PREC)-0.5d0) / stx

            ibz = ibz + 1

            bz%ik (1, ibz) = i1
            bz%ik (2, ibz) = 0

            bz%kweight (ibz) = 1.0d0
            tw = tw + 1.0d0

            Do j = 1, 2
               bz%k (j, ibz) = fx1 * vbg (j, 1)
            End Do

         End Do

         Go To 290


1188     ibz = 1
         bz%k (1, 1) = kpoint (1)
         bz%k (2, 1) = kpoint (2)
         bz%ik (1, 1) = 1
         bz%ik (2, 1) = 1
         bz%kweight (1) = 1.0d0
         tw = 1.0d0
!!$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         Go To 290


290      If (ibz /= nbz) Go To 299

         Do ibz = 1, nbz
            bz%kweight (ibz) = bz%kweight(ibz) / tw
         End Do

         bz%nkgrid = nbz

         If (scale > 0.0d0 .And. scale /= 1.0d0) Then
            bz%k = bz%k * Sqrt (scale)
            bz%kweight = bz%kweight * scale
         End If


         Write (cwork, '(i8)') nbz
         Call do_log (1, '   NUMBER OF K||-POINTS:   NBZ='//trim(cwork))
         Call do_log (1, ' Done!')
         Call do_log (1, '')

         Return

299      Write (cwork, '(i6)') ibz
         Call do_log (1, '   **** ERROR IN GIBZ:  IBZ /= NBZ,  IBZ='//trim(cwork))
         Stop
      End Subroutine gibz_2d


End Module sk_ham


