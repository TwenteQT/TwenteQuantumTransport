!!$ $Id: structure.F90 2000 2013-11-20 15:11:33Z antst $
#include "math_def.h"

!!$ #define _FORCE_SC_ALPHA_HERMICITY_ 1

Module structure
      use omta_defs
      Use sparselib
      Implicit None

      Type t_str_block
         Real (Kind=DEF_DBL_PREC), Allocatable :: bl (:, :)
         Integer :: alloc = 0, ncol, nrow, nrm
      End Type t_str_block

      Type t_mathop
         Type (zcsrmat) :: c, l, r
         Integer :: alloc = 0
         Integer :: n, nl, nr, havehops = 0
      End Type t_mathop

      Type t_str_site
!!$ nnb - number of real neighbours (trpar=0)
!!$ nnsite - number of real neighbour which will be presented in sk
!!$ startn - start index of hopping to this real neighbour
!!$ nhop - number of hopping elements
!!$ trpar - translations index to hopping elements
!!$ idx - indexes of atom in matrix (indexes of diagonal element in resulting matrix)
!!$ ncols - number of columns in the rows corresponding to this atom
!!$ norb - number of orbitals in this atom
         Integer :: nnb, nrows, ncols, alloc = 0, nhop, srow, ia
         Type (t_str_block), Allocatable :: hop (:)
         Integer, Allocatable :: rowind (:), startn (:), trpar (:, :)
      End Type t_str_site

      Type t_strconst
         Integer :: nbl, ndata = 0, alloc = 0
         Integer :: nm = 0, nnz = 0, nrows = 0
         Integer :: nl = 0, lnnz = 0, lnrows = 0
         Integer :: nr = 0, rnnz = 0, rnrows = 0
         Type (t_str_site), Allocatable :: main (:), lhop (:), rhop (:)
         Real (Kind=DEF_DBL_PREC) :: base (2, 2), ltr (2), rtr (2)
      End Type t_strconst

      Type t_neighbour
!!$ datatype which keep description of neighbour
!!$ field "next" is a pointer to next neighbour
         Integer :: nl, atnum, trpar (2), trv (3), trperp, nrm = 0, clsnum
         Real (Kind=DEF_DBL_PREC) :: rvec (3), scale
         Real (Kind=DEF_DBL_PREC), Pointer :: alpha (:)
         Real (Kind=DEF_DBL_PREC), Pointer :: sr(:,:)
         Real (Kind=DEF_DBL_PREC), Pointer :: srd(:,:)
         Real (Kind=DEF_DBL_PREC), Pointer :: sdotr(:,:)
         Real (Kind=DEF_DBL_PREC), Pointer :: t1(:), t2(:), t3(:), t4(:), a(:)
         Real (Kind=DEF_DBL_PREC), Pointer :: t5(:), t6(:), t7(:), t8(:)
         Real (Kind=DEF_DBL_PREC), Pointer :: alphadot(:)
         Real (Kind=DEF_DBL_PREC), Pointer :: logder(:)
         Integer , Pointer :: idxdn(:)
         Integer :: numdf = 0
         Type (t_neighbour), Pointer :: next
      End Type t_neighbour

      Type t_neighbour_list
!!$ datatype which keep list of neighbours
!!$ initially it is organized as a list, where "first" points to first neighbour
!!$ and in each "t_neighbour" field "next" is a pointer to next "t_neighbour"
         Integer :: num = 0, cl_idx, maxnl, basenl
         Type (t_neighbour), Pointer :: nbr (:)
         Type (t_neighbour), Pointer :: first, last
         Integer, Pointer :: cl_dir_idx (:)
      End Type t_neighbour_list

      Type t_cluster_site
         Real (Kind=DEF_DBL_PREC) :: rvec (3), scale
         Real (Kind=DEF_DBL_PREC), Pointer :: alpha (:)
         Real (Kind=DEF_DBL_PREC), Pointer :: sr (:, :)
         Real (Kind=DEF_DBL_PREC), Pointer :: srd (:, :)
         Real (Kind=DEF_DBL_PREC), Pointer :: sdotr (:, :)
         Real (Kind=DEF_DBL_PREC), Pointer :: t1(:), t2(:), t3(:), t4(:) 
         Real (Kind=DEF_DBL_PREC), Pointer :: t5(:), t6(:), t7(:), t8(:) !!$ the
!!$ transformation matrix
         Real (Kind=DEF_DBL_PREC), Pointer :: a(:) !!$ as in K = a*(D-S)
         Real (Kind=DEF_DBL_PREC), Pointer :: alphadot(:)
         Real (Kind=DEF_DBL_PREC), Pointer :: logder(:)
         Integer , Pointer :: idxdn(:)
         Integer :: nl, nrm, atnum, numdf = 0
      End Type t_cluster_site

      Type t_cluster
         Integer :: nsites, nl, numorb, maxorb, numnrm
         Real (Kind=DEF_DBL_PREC), Pointer :: sr (:, :)
         Real (Kind=DEF_DBL_PREC), Pointer :: srd (:, :)
         Real (Kind=DEF_DBL_PREC), Pointer :: sdotr (:, :)
         Type (t_cluster_site), Pointer :: site (:)
         Type (t_cluster), Pointer :: next
      End Type t_cluster

      Type t_cluster_list
         Integer :: ncls = 0, maxnl = 0
         Type (t_cluster), Pointer :: first, last, cls (:)
      End Type t_cluster_list

      Type t_remember_neibs
         Integer, Pointer :: ai(:), tr(:, :)
         Integer :: n, ofs, alloc = 0
      End Type t_remember_neibs


#define E_SAME_COORD 1.0d-8

Contains

      Subroutine alloc_rneib(rneib,n)
         Implicit None
         Type(t_remember_neibs) :: rneib
         Integer :: n
         If (rneib%alloc /= 0) Then
            Deallocate(rneib%ai)
            Deallocate(rneib%tr)
         End If
         allocate(rneib%ai(n))
         allocate(rneib%tr(2,n))
         rneib%n = n
         rneib%alloc = 1
      End Subroutine alloc_rneib
      
      Subroutine make_rneibs(neibs, geom, rneibs)
         Use geometry_module
         Implicit None
         Type (t_remember_neibs), Intent(InOut) :: rneibs(:)
         Type (t_neighbour_list), Pointer, Intent(In) :: neibs (:)
         Type (t_geometry), Intent(In) :: geom
         Integer :: i, k, n, ofs, m
         
         ofs = 0
         Do i = 1, geom%num
            rneibs(i)%ofs = ofs
            n = neibs(i)%num
            Do k = 1, neibs(i)%num
               If (neibs(i)%nbr(k)%trperp /= 0) n = n - 1
            End Do
            Call alloc_rneib (rneibs(i), n)
            m = 0
            Do k = 1, neibs(i)%num
               If (neibs(i)%nbr(k)%trperp == 0) Then
                  m = m + 1
                  rneibs(i)%ai(m) = neibs(i)%nbr(k)%atnum
                  rneibs(i)%tr(:, m) = neibs(i)%nbr(k)%trpar(:)
               End If
            End Do
            ofs = ofs + 2 * geom%atoms(i)%ptr%nl**2
         End Do
      End Subroutine make_rneibs

      Subroutine make_rneibs_EMTO(neibs, geom, rneibs)
         Use geometry_module
         Implicit None
         Type (t_remember_neibs), Intent(InOut) :: rneibs(:)
         Type (t_neighbour_list), Pointer, Intent(In) :: neibs (:)
         Type (t_geometry_EMTO), Intent(In) :: geom
         Integer :: i, k, n, ofs, m
         Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

         ofs = 0
         Do i = 1, geom%num
            rneibs(i)%ofs = ofs
            n = neibs(i)%num
            Do k = 1, neibs(i)%num
               If (neibs(i)%nbr(k)%trperp /= 0) n = n - 1
            End Do
            Call alloc_rneib (rneibs(i), n)
            m = 0
            Do k = 1, neibs(i)%num
               If (neibs(i)%nbr(k)%trperp == 0) Then
                  m = m + 1
                  rneibs(i)%ai(m) = neibs(i)%nbr(k)%atnum
                  rneibs(i)%tr(:, m) = neibs(i)%nbr(k)%trpar(:)
               End If
            End Do
            ofs = ofs + 2 * (geom%atoms(i)%ptr%lmx+1)**2 - 2*dot_product(geom%atoms(i)%ptr%idxdn-1,numls(1:geom%atoms(i)%ptr%lmx+1))
         End Do
      End Subroutine make_rneibs_EMTO

      Subroutine init_mathop (m)
         Implicit None
         Type (t_mathop) :: m
         m%alloc = 0
         m%havehops = 0
         Call init (m%l)
         Call init (m%r)
         Call init (m%c)
      End Subroutine

! This function composes a new mathop struct out of four incoming mathop structs
! each matrix will be of the form
! ( 1  3 )
! ( 4  2 )
! where the number is the nth component of the array m of mathop structs
      Function compose_mathop (m, num) Result (s)
        Implicit None
        Integer :: num
        Type (t_mathop) :: m(1:num), s
        type (zcsrmat) :: temp, temp2
        
        s%alloc = 1
        s%havehops = 1  
        s%nl = m(1)%nl
        s%nr = m(1)%nr
        s%n =  m(1)%n 
        if (num == 4) then
          temp = collectspmat(m(1)%l, m(2)%l)
          temp2= collectspmat_offdiag(m(3)%l, m(4)%l)
          s%l = spmatadd(temp,temp2)

          temp = collectspmat(m(1)%r, m(2)%r)
          temp2= collectspmat_offdiag(m(3)%r, m(4)%r)
          s%r = spmatadd(temp,temp2)

          temp = collectspmat(m(1)%c, m(2)%c)
          temp2= collectspmat_offdiag(m(3)%c, m(4)%c)
          s%c = spmatadd(temp,temp2)
        else if (num == 2) then
          !s%l = collectspmat(m(1)%l, m(2)%l)
          !s%r = collectspmat(m(1)%r, m(2)%r)
          !s%c = collectspmat(m(1)%c, m(2)%c)
          call blockdiag(m(1)%l, m(2)%l, s%l)
          call blockdiag(m(1)%r, m(2)%r, s%r)
          call blockdiag(m(1)%c, m(2)%c, s%c)
          !s%nl = m(1)%nl + m(2)%nl
          !s%nr = m(1)%nr + m(2)%nr
          !s%n =  m(1)%n + m(2)%n
        endif

      End Function

      Function calc_screal (geom, irel, cr, isym_in, rneibs) Result (scalpha)
         Use geometry_module
         Use logging
         Implicit None
         Type (t_geometry) :: geom
         Type (t_strconst) :: scalpha
         Integer, Optional :: irel, isym_in
         Real (Kind=DEF_DBL_PREC), Optional :: cr
         Type (t_remember_neibs), Optional, Intent(InOut) :: rneibs(:)
!!$ local

         Type (t_neighbour_list), Pointer :: neibs (:)
         Type (t_cluster_list) :: clslist
         Integer :: numat, lirel

         If ( .Not. present(irel)) Then
            lirel = 0
         Else
            lirel = irel
         End If

         Call do_log (1, 'S(k) calculation:')

         numat = geom%num

!!$ find neighbours for each site
         neibs => find_neibs (geom, cr)
         If (present(rneibs)) Then
            Call make_rneibs(neibs, geom, rneibs)
         End If

!!$ find number of unique cluster in whole structure
         clslist = find_unique_clusters (numat, neibs)

!!$ calculate s^(\alpha)(r,r') for each cluster
         Call do_calc_strconst (clslist, isym_in)

!!$ report results
         Call do_rep_cls (clslist)

         Call do_rep_clsgem (geom, neibs, clslist)

!!$ construct data structure with s^(\alpha)(r,r') for whole structure
         If (lirel == 0) Then
            scalpha = prep_sc_real (geom, neibs, clslist, iop=1)
         Else
            scalpha = prep_sc_real_rel (geom, neibs, clslist, iop=1)
         End If

!!$    call symmetrize_sca (norbit,neibs,clslist)

!!$ deallocate memory
         Call free_neibs (neibs)
         Deallocate (neibs)
         Call free_cluster_list (clslist,iop=1)
         Call do_log (1, '')
         Return
      End Function calc_screal

      Function calc_screal_EMTO (geom, kap2, en, irel, needcur, df, ppar, cr, isym_in, rneibs, sdotca) Result (sca)
         Use geometry_module
         Use logging
         Implicit None
         Type (t_geometry_EMTO) :: geom
         Type (t_strconst) :: sca
         Real (Kind=DEF_DBL_PREC) :: kap2
         Real (Kind=DEF_DBL_PREC) :: en 
         Integer :: needcur ! 1 if calculating sdot, 0 otherwise
         Integer, Optional :: irel, isym_in, df
         Type (t_omta_logder), Optional :: ppar !contains logder and hsoc, needed for downfolding
         Real (Kind=DEF_DBL_PREC), Optional :: cr
         Type (t_remember_neibs), Optional, Intent(InOut) :: rneibs(:)
         Type (t_strconst), Optional :: sdotca ! if needcur=1, Sdot is returned in sdotca
!!$ local

         Type (t_neighbour_list), Pointer :: neibs (:)
         Type (t_cluster_list) :: clslist
         Integer :: numat, lirel
         Type (t_geometry) :: geom_LMTO
         Integer :: i
         Integer :: iop
                

         If ( .Not. present(irel)) Then
            lirel = 0
         Else
            lirel = irel
         End If

         Call do_log (1, 'S(k) calculation:')

         numat = geom%num

!!$ find neighbours for each site
         neibs => find_neibs_EMTO (geom, cr)
         If (present(rneibs)) Then
            Call make_rneibs_EMTO(neibs, geom, rneibs)
         End If

!!$ find number of unique cluster in whole structure
         clslist = find_unique_clusters (numat, neibs)
!!$ calculate s^(a)(r,r') for each cluster
         if (present(ppar)) then
           iop = 1
           if (needcur > 0) iop = 2
           Call do_calc_strconst_EMTO (clslist, geom%dawsr, kap2, isym_in, en, ppar, iop)
         else
           Call do_calc_strconst_EMTO (clslist, geom%dawsr, kap2, isym_in, en)
         endif

!!$ report results
         Call do_rep_cls (clslist)

!!$ REALLY HACKY FIX - right now, there are two geometry types, one for LMTO and
!one for EMTO. These should be merged in the future, so that we can get rid of
!all of the duplicate subroutines (there are plenty that are trivially
!duplicated, just so that we can input the _EMTO type.) Instead of duplicating
!the prep_sc_real and prep_sc_real_rel subroutines (they are way too long for
!this), I opted to do this:

         geom_LMTO%num = geom%num
         geom_LMTO%r_transnum = geom%r_transnum
         geom_LMTO%l_transnum = geom%l_transnum
         geom_LMTO%base = geom%base !old
         !geom_LMTO%base = geom%base*geom%alat
         geom_LMTO%l_perp_trans = geom%l_perp_trans 
         geom_LMTO%r_perp_trans = geom%r_perp_trans 
         allocate(geom_LMTO%atoms(size(geom%atoms)))
         do i=1,size(geom%atoms)
          geom_LMTO%atoms(i)%coord = geom%atoms(i)%coord 
          allocate(geom_LMTO%atoms(i)%ptr)
          geom_LMTO%atoms(i)%ptr%label = geom%atoms(i)%ptr%label
         enddo
!!$ A truly disgusting fix. But it should work.
         Call do_rep_clsgem (geom_LMTO, neibs, clslist)


!!$ construct data structure with s^(\alpha)(r,r') for whole structure
         If (lirel == 0) Then
            If (needcur == 0) sca = prep_sc_real (geom_LMTO, neibs, clslist, iop=1)
            If (needcur > 0) sca = prep_sc_real (geom_LMTO, neibs, clslist, iop=2, sdotcr=sdotca)
         Else
            If (needcur == 0) sca = prep_sc_real_rel (geom_LMTO, neibs, clslist, iop=1)
            If (needcur > 0) sca = prep_sc_real_rel (geom_LMTO, neibs, clslist, iop=2, sdotcr=sdotca)
         End If
!!$    call symmetrize_sca (norbit,neibs,clslist)

!!$ deallocate memory
         Call free_neibs (neibs)
         Deallocate (neibs)
         Call free_cluster_list (clslist, iop=needcur+1)
         Call do_log (1, '')
         do i=1,size(geom%atoms)
          deallocate(geom_LMTO%atoms(i)%ptr)
         enddo
         deallocate(geom_LMTO%atoms)

         Return
      End Function calc_screal_EMTO


      Function prep_sc_real (geom, neibs, clslist, iop, sdotcr) Result (scr)
!!$ Finally construct datatype "t_strconst" which holds all information which
!!$ we need to calculate S(k) for region and hopings to outside
         Use geometry_module
         Use logging
         Implicit None
         Type (t_neighbour_list), Target :: neibs (:)
         Type (t_cluster_list) :: clslist
         Type (t_strconst), Target :: scr
         Type (t_geometry) :: geom
         Integer :: iop ! 1 -> S, 2 -> Sdot
         Type (t_strconst), Target, Optional :: sdotcr !!$ returns if iop=2
!!$ local
         Integer :: nat, lnat, rnat, il, ir
         Integer :: i, norb, rowidx (clslist%maxnl**2), startrow (geom%num)
         Integer :: nlhop, nrhop, numb (3, geom%num), lnumb (3, geom%num), rnumb (3, geom%num)
         Type (t_neighbour_list), Pointer :: t_neib
         Type (t_str_site), Pointer :: t_ss
         Real (Kind=DEF_DBL_PREC) :: ldev, dev, rdev
         Character (Len=100) :: cwork

         nat = geom%num
         lnat = geom%l_transnum
         rnat = geom%r_transnum

         Allocate (scr%main(nat))
         If (iop==2) Allocate(sdotcr%main(nat))

         scr%nm = nat
         scr%base = geom%base
         scr%ltr = geom%l_perp_trans(1:2)
         scr%rtr = geom%r_perp_trans(1:2)
         If (iop==2) Then
          sdotcr%nm = nat
          sdotcr%base = geom%base
          sdotcr%ltr = geom%l_perp_trans(1:2)
          sdotcr%rtr = geom%r_perp_trans(1:2)
         Endif

         Do i = 1, clslist%maxnl ** 2
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

         Allocate (scr%lhop(nlhop))
         Allocate (scr%rhop(nrhop))
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
            t_ss => scr%main (i)
            scr%nnz = scr%nnz + fillstr (0, t_ss, neibs, clslist, numb(:, i), rowidx, startrow, i, dev, 1)
            If (lnumb(1, i) > 0) Then
               il = il + 1
               t_ss => scr%lhop (il)
               scr%lnnz = scr%lnnz + fillstr (-1, t_ss, neibs, clslist, lnumb(:, i), rowidx, startrow, i, &
              & ldev, 1)
            End If
            If (rnumb(1, i) > 0) Then
               ir = ir + 1
               t_ss => scr%rhop (ir)
               scr%rnnz = scr%rnnz + fillstr (1, t_ss, neibs, clslist, rnumb(:, i), rowidx, startrow, i, &
              & rdev, 1)
            End If
         End Do

         Call do_log (1, 'Maximal deviations in S:')
         Write (cwork, '(3(3xf8.6))') ldev, dev, rdev
         Call do_log (1, '  '//trim(cwork))

         If (iop==2) Then

          Allocate (sdotcr%lhop(nlhop))
          Allocate (sdotcr%rhop(nrhop))
          sdotcr%nl = nlhop
          sdotcr%nr = nrhop
          sdotcr%nrows = norb
          sdotcr%lnrows = norb
          sdotcr%rnrows = norb
          sdotcr%nnz = 0
          il = 0
          ir = 0
          ldev = 0.0d0
          rdev = 0.0d0
          dev = 0.0d0
          Do i = 1, nat
             t_ss => sdotcr%main (i)
             sdotcr%nnz = sdotcr%nnz + fillstr (0, t_ss, neibs, clslist, numb(:, i), rowidx, startrow, i, dev, 2)
             If (lnumb(1, i) > 0) Then
                il = il + 1
                t_ss => sdotcr%lhop (il)
                sdotcr%lnnz = sdotcr%lnnz + fillstr (-1, t_ss, neibs, clslist, lnumb(:, i), rowidx, startrow, i, &
               & ldev, 2)
             End If
             If (rnumb(1, i) > 0) Then
                ir = ir + 1
                t_ss => sdotcr%rhop (ir)
                sdotcr%rnnz = sdotcr%rnnz + fillstr (1, t_ss, neibs, clslist, rnumb(:, i), rowidx, startrow, i, &
               & rdev, 2)
             End If
          End Do

          Call do_log (1, 'Maximal deviations in Sdot:')
          Write (cwork, '(3(3xf8.6))') ldev, dev, rdev
          Call do_log (1, '  '//trim(cwork))


         Endif

      Contains

         Function fillstr (t, t_ss, neibs, clslist, numb, rowidx, startrow, sn, dev, iop) Result (nnz)
            Implicit None
            Type (t_neighbour_list), Target :: neibs (:)
            Type (t_cluster_list) :: clslist
            Type (t_str_site), Target :: t_ss
            Integer :: numb (3), t, rowidx (:), startrow (:), nnz, sn, iop
!!$ Local
            Type (t_cluster), Pointer :: t_cls, j_cls
            Type (t_neighbour_list), Pointer :: t_neib
            Integer :: i, norb, atn, nhop, jnorb, dir, nnb, st, st1, j, jat, jdir
            Type (t_neighbour), Pointer :: t_nbr
            Real (Kind=DEF_DBL_PREC) :: rvec (3), dev


            t_neib => neibs (sn)
            t_cls => clslist%cls (t_neib%cl_idx)

            Call allocstr (t_ss, numb)
            norb = t_cls%nl ** 2
            t_ss%nrows = norb
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
                  dir = t_neib%cl_dir_idx (i)
                  jnorb = t_cls%site(dir)%nl ** 2
                  jat = t_nbr%atnum
                  If (jat > atn) Then
                     atn = jat
                     nnb = nnb + 1
                     st = st1
                     t_ss%rowind (st:st+jnorb-1) = rowidx (1:jnorb) + startrow (atn)
                     st1 = st1 + jnorb
                  End If
                  Allocate (t_ss%hop(nhop)%bl(norb, jnorb))
!!$ Find symmetric element

                  If (t_nbr%nrm == 0 .and. iop == 1) Then
                     j_cls => clslist%cls (neibs(jat)%cl_idx)
                     rvec = t_cls%site(dir)%rvec
                     Do j = 1, j_cls%nsites
                        If (sum(Abs(rvec+j_cls%site(j)%rvec)) < E_SAME_COORD) Then
                           If (iop==1) then
                            dev = Max (dev, maxval(Abs(t_cls%site(dir)%sr-(transpose(j_cls%site(j)%sr)))))
                            t_ss%hop(nhop)%bl = 0.5d0 * (t_cls%site(dir)%sr+(transpose(j_cls%site(j)%sr)))
                           Else if (iop==2) then
                            dev = Max (dev, maxval(Abs(t_cls%site(dir)%sdotr-(transpose(j_cls%site(j)%sdotr)))))
                            t_ss%hop(nhop)%bl = 0.5d0 * (t_cls%site(dir)%sdotr+(transpose(j_cls%site(j)%sdotr)))
                           Else
                            write(*,*) 'iop has faulty value in prep_sc_real, abort!'
                            stop
                           Endif
                           jdir = j
                           Exit
                        End If
                     End Do
                     If (t == 0) Then
                        Do j = 1, neibs(jat)%num
                           If (neibs(jat)%cl_dir_idx(j) == jdir) Then
                              If (neibs(jat)%nbr(j)%nrm == 1 .and. iop == 1) Then
                                 Write (*,*) 'Sure that everything is OK with structure and cutrat?'
                                 Write (*,*) 'It isn''t really error, but...'
                              End If
                              neibs(jat)%nbr(j)%nrm = 1
!!$                               neibs(jat)%nbr(j)%nrm
                                  if (iop==1) neibs(jat)%nbr(j)%sr => t_ss%hop(nhop)%bl
                                  if (iop==2) neibs(jat)%nbr(j)%sdotr => t_ss%hop(nhop)%bl
!!$                               t_nbr%nrm=1
                              Exit
                           End If
                        End Do
                     End If
                  Else
                    if (iop==1) t_ss%hop(nhop)%bl = (transpose (t_nbr%sr))
                    if (iop==2) t_ss%hop(nhop)%bl = (transpose (t_nbr%sdotr))

!!$                write(*,*) 'got it!'
                  End If
                  t_ss%hop(nhop)%ncol = jnorb
                  t_ss%hop(nhop)%nrow = norb
                  t_ss%trpar (:, nhop) = t_nbr%trpar
                  t_ss%startn (nhop) = st
               End If
            End Do

            nnz = norb * numb (3)
!!$             write(*,*) norb,numb(3)
            Return
         End Function fillstr

         Subroutine allocstr (t_ss, numb)
            Implicit None
            Type (t_str_site) :: t_ss
            Integer :: numb (3)
            Integer :: nhop, nnb, norb

            t_ss%nhop = numb (2)
            t_ss%nnb = numb (1)
            t_ss%ncols = numb (3)
            t_ss%alloc = 1
            Allocate (t_ss%hop(t_ss%nhop))
            Allocate (t_ss%startn(t_ss%nhop))
            Allocate (t_ss%trpar(2, t_ss%nhop))
            Allocate (t_ss%rowind(t_ss%ncols))

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
      End Function prep_sc_real


      Function prep_sc_real_rel (geom, neibs, clslist, iop, sdotcr) Result (scr)
!!$ Finally construct datatype "t_strconst" which holds all information which
!!$ we need to calculate S(k) for region and hopings to outside
         Use geometry_module
         Use logging
         Implicit None
         Type (t_neighbour_list), Target :: neibs (:)
         Type (t_cluster_list) :: clslist
         Type (t_strconst), Target :: scr
         Type (t_geometry) :: geom
         Integer :: iop ! 1 -> S, 2 -> Sdot
         Type (t_strconst), Target, Optional :: sdotcr
!!$ local
         Integer :: nat, lnat, rnat, il, ir
         Integer :: i, norb, rowidx (clslist%maxnl**2), startrow (geom%num)
         Integer :: nlhop, nrhop, numb (3, geom%num), lnumb (3, geom%num), rnumb (3, geom%num)
         Type (t_neighbour_list), Pointer :: t_neib
         Type (t_str_site), Pointer :: t_ss (:)
         Real (Kind=DEF_DBL_PREC) :: ldev, dev, rdev
         Character (Len=100) :: cwork
         Integer :: numls(5) = [1,3,5,7,9]

         nat = geom%num
         lnat = geom%l_transnum
         rnat = geom%r_transnum

         Allocate (scr%main(nat*2))
         If (iop == 2) Allocate(sdotcr%main(nat*2))

         scr%nm = nat * 2
         scr%base = geom%base
         scr%ltr = geom%l_perp_trans(1:2)
         scr%rtr = geom%r_perp_trans(1:2)

         If (iop == 2) Then
          sdotcr%nm = nat * 2
          sdotcr%base = geom%base
          sdotcr%ltr = geom%l_perp_trans(1:2)
          sdotcr%rtr = geom%r_perp_trans(1:2)
         Endif

         Do i = 1, clslist%maxnl ** 2
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
            if (t_neib%first%numdf > 1) then
             norb = norb + t_neib%basenl ** 2 - t_neib%first%numdf 
            else
             norb = norb + t_neib%basenl ** 2
            endif
         End Do

         If (lnat < nlhop .Or. rnat < nrhop) Then
            Write (*,*) 'Definitely something wrong with defenitions left-,right- edges'
            Write (*,*) 'Num.atoms with left-,right- hopping is ', nlhop, nrhop
            Write (*,*) 'Size of left-,right- layers is', lnat, rnat
            Write (*,*) 'Check the geometry of your leads w.r.t. the scattering region!'
            Stop
         End If

         Allocate (scr%lhop(nlhop*2))
         Allocate (scr%rhop(nrhop*2))
         scr%nl = nlhop * 2
         scr%nr = nrhop * 2
         scr%nrows = norb * 2
         scr%lnrows = norb * 2
         scr%rnrows = norb * 2
         scr%nnz = 0
         il = 0
         ir = 0
         ldev = 0.0d0
         rdev = 0.0d0
         dev = 0.0d0

         Do i = 1, nat
            t_ss => scr%main ((i-1)*2+1:(i-1)*2+2)
            scr%nnz = scr%nnz + fillstr (0, t_ss, neibs, clslist, numb(:, i), rowidx, startrow, i, dev, 1)
            If (lnumb(1, i) > 0) Then
               il = il + 1
               t_ss => scr%lhop (2*(il-1)+1:2*(il-1)+2)
               scr%lnnz = scr%lnnz + fillstr (-1, t_ss, neibs, clslist, lnumb(:, i), rowidx, startrow, i, &
              & ldev, 1)
            End If
            If (rnumb(1, i) > 0) Then
               ir = ir + 1
               t_ss => scr%rhop (2*(ir-1)+1:2*(ir-1)+2)
               scr%rnnz = scr%rnnz + fillstr (1, t_ss, neibs, clslist, rnumb(:, i), rowidx, startrow, i, &
              & rdev, 1)
            End If
         End Do

         Call do_log (1, 'Maximal deviations in S:')
         Write (cwork, '(3(3xf9.6))') ldev, dev, rdev
         Call do_log (1, '  '//trim(cwork))

         If (iop == 2) Then

         Allocate (sdotcr%lhop(nlhop*2))
         Allocate (sdotcr%rhop(nrhop*2))
         sdotcr%nl = nlhop * 2
         sdotcr%nr = nrhop * 2
         sdotcr%nrows = norb * 2
         sdotcr%lnrows = norb * 2
         sdotcr%rnrows = norb * 2
         sdotcr%nnz = 0
         il = 0
         ir = 0
         ldev = 0.0d0
         rdev = 0.0d0
         dev = 0.0d0

         Do i = 1, nat
            t_ss => sdotcr%main ((i-1)*2+1:(i-1)*2+2)
            sdotcr%nnz = sdotcr%nnz + fillstr (0, t_ss, neibs, clslist, numb(:, i), rowidx, startrow, i, dev, 2)
            If (lnumb(1, i) > 0) Then
               il = il + 1
               t_ss => sdotcr%lhop (2*(il-1)+1:2*(il-1)+2)
               sdotcr%lnnz = sdotcr%lnnz + fillstr (-1, t_ss, neibs, clslist, lnumb(:, i), rowidx, startrow, i, &
              & ldev, 2)
            End If
            If (rnumb(1, i) > 0) Then
               ir = ir + 1
               t_ss => sdotcr%rhop (2*(ir-1)+1:2*(ir-1)+2)
               sdotcr%rnnz = sdotcr%rnnz + fillstr (1, t_ss, neibs, clslist, rnumb(:, i), rowidx, startrow, i, &
              & rdev, 2)
            End If
         End Do

         Call do_log (1, 'Maximal deviations in Sdot:')
         Write (cwork, '(3(3xf9.6))') ldev, dev, rdev
         Call do_log (1, '  '//trim(cwork))


         Endif


      Contains

         Function fillstr (t, t_ss2, neibs, clslist, numb, rowidx, startrow, sn, dev, iop) Result (nnz)
            Implicit None
            Type (t_neighbour_list), Target :: neibs (:)
            Type (t_cluster_list) :: clslist
            Type (t_str_site), Target :: t_ss2 (:)
            Integer :: numb (3), t, rowidx (:), startrow (:), nnz, sn, is, iop
!!$ Local
            Type (t_cluster), Pointer :: t_cls, j_cls
            Type (t_neighbour_list), Pointer :: t_neib
            Integer :: i, norb, atn, nhop, jnorb, dir, nnb, st, st1, j, jat, jdir
            Type (t_neighbour), Pointer :: t_nbr
            Real (Kind=DEF_DBL_PREC) :: rvec (3), dev
            Type (t_str_site), Pointer :: t_ss, t_ssd
            Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

            t_ss => t_ss2 (1)
            t_ssd => t_ss2 (2)
            t_neib => neibs (sn)
            t_cls => clslist%cls (t_neib%cl_idx)
            Call allocstr (t_ss, numb)
            Call allocstr (t_ssd, numb)

            if (size(t_cls%site(1)%idxdn) > 1) then
             norb = t_cls%nl ** 2 - t_cls%site(1)%numdf
            else 
             norb = t_cls%nl ** 2
            endif

            t_ss%nrows = norb
            t_ssd%nrows = norb

            atn = 0
            nhop = 0
            nnb = 0
            st1 = 1
            t_ss%srow = 2 * startrow (sn) + 1
            t_ssd%srow = 2 * startrow (sn) + 1 + norb

            t_ss%ia = sn
            t_ssd%ia = sn

            Do i = 1, t_neib%num
               t_nbr => t_neib%nbr (i)
               If (t_nbr%trperp == t) Then
                  nhop = nhop + 1
                  dir = t_neib%cl_dir_idx (i)
                  if (size(t_cls%site(dir)%idxdn) > 1) then
                  jnorb = t_cls%site(dir)%nl ** 2 - t_cls%site(dir)%numdf 
                  else
                   jnorb = t_cls%site(dir)%nl **2
                  endif
                  jat = t_nbr%atnum
                  If (jat > atn) Then
                     atn = jat
                     nnb = nnb + 1
                     st = st1
                     t_ss%rowind (st:st+jnorb-1) = rowidx (1:jnorb) + 2 * startrow (atn)
                     t_ssd%rowind (st:st+jnorb-1) = rowidx (1:jnorb) + 2 * startrow (atn) + jnorb
                     st1 = st1 + jnorb
                  End If

                  Allocate (t_ss%hop(nhop)%bl(norb, jnorb))
                  Allocate (t_ssd%hop(nhop)%bl(norb, jnorb))
!!$ Find symmetric element
                  If ((t_nbr%nrm == 0 .and. iop == 1) .or. iop == 2) Then
                     j_cls => clslist%cls (neibs(jat)%cl_idx)
                     rvec = t_cls%site(dir)%rvec
                     Do j = 1, j_cls%nsites
                        If (sum(Abs(rvec+j_cls%site(j)%rvec)) < E_SAME_COORD) Then
                          If (iop == 1) Then

                            dev = Max (dev, maxval(Abs(t_cls%site(dir)%sr-(transpose(j_cls%site(j)%sr)))))
                            t_ss%hop(nhop)%bl = 0.5d0 * (t_cls%site(dir)%sr+(transpose(j_cls%site(j)%sr)))
                            t_ssd%hop(nhop)%bl = 0.5d0 * (t_cls%site(dir)%srd+(transpose(j_cls%site(j)%srd)))
                           Elseif (iop == 2) Then
                            dev = Max (dev, maxval(Abs(t_cls%site(dir)%sdotr-(transpose(j_cls%site(j)%sdotr)))))
                            t_ss%hop(nhop)%bl = 0.5d0 * (t_cls%site(dir)%sdotr+(transpose(j_cls%site(j)%sdotr)))
                           Else
                            write(*,*) 'iop has faulty value in prep_sc_real, abort!'
                            stop
                           Endif
                           jdir = j
                           Exit
                        End If
                     End Do

                     If (t == 0) Then
                        Do j = 1, neibs(jat)%num
                           If (neibs(jat)%cl_dir_idx(j) == jdir) Then
                              If (neibs(jat)%nbr(j)%nrm == 1 .and. iop == 1) Then
                                 Write (*,*) 'Sure that everything is OK with structure and cutrat?'
                                 Write (*,*) 'It is not really error, but...'
                              End If
                              neibs(jat)%nbr(j)%nrm = 1
                              if (iop==1) then
                                neibs(jat)%nbr(j)%sr => t_ss%hop(nhop)%bl
                                neibs(jat)%nbr(j)%srd => t_ssd%hop(nhop)%bl
                              endif
                              if (iop==2) neibs(jat)%nbr(j)%sdotr => t_ss%hop(nhop)%bl
                              Exit
                           End If
                        End Do
                     End If
                  Else
                    if (iop==1) then
                      t_ss%hop(nhop)%bl = (transpose (t_nbr%sr))
                      t_ssd%hop(nhop)%bl = (transpose (t_nbr%srd))
                    endif
                    if (iop==2) t_ss%hop(nhop)%bl = (transpose (t_nbr%sdotr))
!!$               write(*,*) 'got it!'
                  End If
                  t_ss%hop(nhop)%ncol = jnorb
                  t_ss%hop(nhop)%nrow = norb
                  t_ss%trpar (:, nhop) = t_nbr%trpar
                  t_ss%startn (nhop) = st

                  if (iop==2) t_ssd%hop(nhop)%bl = t_ss%hop(nhop)%bl !old method, just copy the spin 1 to spin 2 
                  t_ssd%hop(nhop)%ncol = t_ss%hop(nhop)%ncol
                  t_ssd%hop(nhop)%nrow = t_ss%hop(nhop)%nrow
                  t_ssd%trpar (:, nhop) = t_ss%trpar(:, nhop)
                  t_ssd%startn (nhop) = t_ss%startn(nhop)

               End If
            End Do

            nnz = 2 * norb * numb (3)
            Return
         End Function fillstr

         Subroutine allocstr (t_ss, numb)
            Implicit None
            Type (t_str_site) :: t_ss
            Integer :: numb (3)
            Integer :: nhop, nnb, norb
            t_ss%nhop = numb (2)
            t_ss%nnb = numb (1)
            t_ss%ncols = numb (3)
            t_ss%alloc = 1
            Allocate (t_ss%hop(t_ss%nhop))
            Allocate (t_ss%startn(t_ss%nhop))
            Allocate (t_ss%trpar(2, t_ss%nhop))
            Allocate (t_ss%rowind(t_ss%ncols))
         End Subroutine allocstr


         Subroutine calc_numbers (t_neib, numb, lnumb, rnumb)
            Implicit None
            Type (t_neighbour_list) :: t_neib
            Integer :: numb (3), lnumb (3), rnumb (3)
!!$ Local
            Type (t_neighbour), Pointer :: t_nbr
            Integer :: i, latn, atn, ratn, atnum, norb
            Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

            numb = 0
            lnumb = 0
            rnumb = 0
            latn = 0
            atn = 0
            ratn = 0
            Do i = 1, t_neib%num
               t_nbr => t_neib%nbr (i)
               atnum = t_nbr%atnum
               !if (SIZE(t_nbr%idxdn)>1) then
               if (t_nbr%numdf>0) then
                !norb = t_nbr%nl ** 2 - dot_product(t_neib%nbr(i)%idxdn-1,numls(1:t_neib%nbr(i)%nl)) 
                norb = t_nbr%nl ** 2 - t_neib%nbr(i)%numdf 
               else
                norb = t_nbr%nl **2
               endif

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
      End Function prep_sc_real_rel


      Subroutine do_rep_cls (clslist)
         Use logging
         Implicit None
         Type (t_cluster_list) :: clslist
!!$ local
         Integer :: i, l
         !Character (Len=1024) :: cwork
         Character (Len=2048) :: cwork
         Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

         Write (cwork, '(i4)') clslist%ncls
         Call do_log (1, 'number of unique clusters is='//trim(cwork))
         Do i = 1, clslist%ncls
            Write (cwork, '("cluster: ",i3," , num.sites=",i3)') i, clslist%cls(i)%nsites
            Call do_log (2, trim(cwork))
            Call do_log (10, "Positions:")
            Do l = 1, clslist%cls(i)%nsites
               Write (cwork, '(" ",3(1x,f7.4))') clslist%cls(i)%site(l)%rvec(1:3)
               Call do_log (10, trim(cwork))
            End Do
            Write (cwork,*) 'diag. elements of s^alpha(0,0):'
            Call do_log (3, '   '//trim(cwork))
            !if (size(clslist%cls(i)%site(1)%idxdn) > 1) then
            if (clslist%cls(i)%site(1)%numdf > 0) then
             Write (cwork, '(" ",1000(1x,f9.4))') (clslist%cls(i)%site(1)%sr(l, l), l=1, &
              & clslist%cls(i)%site(1)%nl**2 - dot_product(clslist%cls(i)%site(1)%idxdn-1,numls(1:clslist%cls(i)%site(1)%nl)))
            else
             Write (cwork, '(" ",1000(1x,f9.4))') (clslist%cls(i)%site(1)%sr(l, l), l=1, &
              & clslist%cls(i)%site(1)%nl**2 )
            endif
            Call do_log (3, '  '//trim(cwork))
            !Write (cwork,*) 'diag. elements of sdot^alpha(0,0):'
            !Call do_log (3, '   '//trim(cwork))
            !Write (cwork, '(" ",100(1x,f9.4))') (clslist%cls(i)%site(1)%sdotr(l, l), l=1, &
           !& clslist%cls(i)%site(1)%nl**2)
            !Call do_log (3, '  '//trim(cwork))
         End Do
      End Subroutine do_rep_cls



      Subroutine do_rep_clsgem (geom, neibs, clslist)
         Use logging
         Use geometry_module
         Use pgsym
         Implicit None
         Type (t_cluster_list), Intent (In) :: clslist
         Type (t_geometry), Intent (In) :: geom
         Type (t_neighbour_list), Pointer, Intent (In) :: neibs (:)

!!$ local
         Integer :: s, i
         Character (Len=1024) :: cwork
         Character (Len=3) :: pgr
         Integer :: els (1024)
         Integer, Parameter :: tll = 10
         Real (Kind=DEF_DBL_PREC) :: pos (3, 1024)
         If (log_level >= tll) Then
            els (:) = 1
            els (1) = 100
            pos (:, 1) = 0.0d0
            Call do_log (tll, 'Geometry cluster report:')
            Do s = 1, geom%num
               Do i = 2, clslist%cls(neibs(s)%cl_idx)%nsites
                  pos (:, i) = clslist%cls(neibs(s)%cl_idx)%site(i)%rvec(:)
               End Do
!!$                write(*,*) s
               Call SymMol (clslist%cls(neibs(s)%cl_idx)%nsites, 1.0d-2, pos, els, pgr)
               If (1 == 0) Then
                  Write (cwork, '("cls",i3.3,".dat")') s
                  Call dumpcls (pos(1:3, 1:clslist%cls(neibs(s)%cl_idx)%nsites), trim(cwork))
               End If
               Write (cwork, '("pos: ",3(1x,g15.9)," at=",a10," CL=",i3," NNB=",i3," SYM=",a4)') &
              & geom%atoms(s)%coord, trim (geom%atoms(s)%ptr%label), neibs(s)%cl_idx, &
              & clslist%cls(neibs(s)%cl_idx)%nsites, trim (pgr)
               Call do_log (tll, trim(cwork))
            End Do
            Call do_log (tll, '-------------------------------------------------')
         End If
      Contains
         Subroutine dumpcls (mat, fname)
            Implicit None
            Real (Kind=DEF_DBL_PREC) :: mat (:, :)
            Character (Len=*) :: fname
            Character (Len=100) :: fmt
            Integer :: s1
            Open (Unit=dumping_unit, File=fname, Action='write')
            Write (Fmt, '("(",i2.2,"(1x,f13.9))")') size (mat, 2)
            Do s1 = 1, size (mat, 2)
               Write (dumping_unit, trim(Fmt)) mat (:, s1) - mat (:, 1)
            End Do
            Close (dumping_unit)
         End Subroutine dumpcls

      End Subroutine do_rep_clsgem


      Subroutine free_mathop (sk)
         Use sparselib
         Implicit None
         Type (t_mathop) :: sk
         If (sk%alloc /= 0) Then
            Call free (sk%c)
            If (sk%havehops /= 0) Then
               Call free (sk%l)
               Call free (sk%r)
               sk%havehops = 0
            End If
            sk%alloc = 0
         End If
      End Subroutine free_mathop


      Subroutine make_drho (sk, scr, kpar)
         Use sparselib
         Implicit None
         Type (t_strconst) :: scr
         Type (zcsrmat) :: sk
         Real (Kind=DEF_DBL_PREC) :: kpar (2)
!!$ local

         Call fillsk (sk, scr%nnz, scr%nrows, scr%nm, scr%main, scr%base, kpar, 1.0d0)
      Contains

         Subroutine fillsk (sk, tnnzi, nrow, na, sites, base, kpar, cfs)
!!$ function which is actually allocate and fill sk
            Implicit None
            Integer :: tnnzi, nrow, na
            Type (zcsrmat) :: sk
            Type (t_str_site), Target :: sites (:)
            Real (Kind=DEF_DBL_PREC) :: base (2, 2)
            Real (Kind=DEF_DBL_PREC) :: kpar (2), cfs
!!$ local
            Type (t_str_site), Pointer :: jsite, isite
            Integer :: ia, nc, nr, i, nnz, j, lr, tnnz
            Integer :: ii

            Call alloc (sk, tnnzi, nrow)
            tnnz = 0
            lr = 1
            ii = 0
            Do ia = 1, na

               isite => sites (ia)
               nr = isite%nrows
               nc = isite%ncols
               nnz = fillskrows (isite, kpar, sk%a(tnnz+1:), base, cfs)
               j = isite%srow
               sk%ir (lr:j) = tnnz + 1
               Do i = 1, nr
                  sk%ir (j+1) = sk%ir (j) + nc
                  j = j + 1
                  sk%jc (ii+1:ii+nc) = isite%rowind
                  ii = ii + nc
               End Do
               lr = j
               tnnz = tnnz + nnz
            End Do
            sk%nnz = tnnz
         End Subroutine fillsk

         Function fillskrows (site, k, val, base, cfs) Result (nnz)
!!$ function which fill rows of sk which correspondent to "site"
            Implicit None
            Type (t_str_site), Target :: site
            Real (Kind=DEF_DBL_PREC) :: k (2), base (2, 2), cfs
            Complex (Kind=DEF_DBL_PREC) :: val (:)
            Integer :: nnz
!!$ local
            Integer :: in, jnc, ir, jnc1
            Complex (Kind=DEF_DBL_PREC) :: els (site%nrows, site%ncols), fact
            Integer, Pointer :: tp (:, :)
            Type (t_str_block), Pointer :: hop (:)

            tp => site%trpar
            hop => site%hop

            els = DEF_cmplx_zero
            Do in = 1, site%nhop
               jnc = site%startn (in)
               jnc1 = jnc + hop(in)%ncol - 1
               fact = cfs * Exp (DEF_cmplx_Ione*Dot_product(k, dble(tp(1, in))*base(:, 1)+dble(tp(2, &
              & in))*base(:, 2)))
               If ((site%srow == site%rowind(jnc)) .And. tp(1, in) == 0 .And. tp(2, in) == 0) Then
!!$         if((site%srow==site%rowind(jnc))) then
                  fact = DEF_cmplx_zero
               End If
               els (:, jnc:jnc1) = els (:, jnc:jnc1) + fact * hop(in)%bl
            End Do

            in = site%ncols
            Do ir = 1, site%nrows
               jnc = (ir-1) * in + 1
               jnc1 = jnc + in - 1
               val (jnc:jnc1) = els (ir, :)
            End Do
            nnz = site%nrows * in
         End Function fillskrows
      End Subroutine make_drho

      Subroutine make_sk (sk, scr, kpar, needhops, cfs)
         Use sparselib
         Implicit None
         Type (t_strconst) :: scr
         Type (t_mathop) :: sk
         Real (Kind=DEF_DBL_PREC), Optional :: cfs
         Integer, Optional :: needhops
         Real (Kind=DEF_DBL_PREC) :: kpar (2)
!!$ local
         Integer :: nh

         If ( .Not. present(needhops)) Then
            nh = 0
         Else
            nh = needhops
         End If

         Call fillsk (sk%c, scr%nnz, scr%nrows, scr%nm, scr%main, scr%base, kpar, cfs)
         If (nh > 0) Then
            Call fillsk (sk%l, scr%lnnz, scr%lnrows, scr%nl, scr%lhop, scr%base, kpar, cfs)
            Call fillsk (sk%r, scr%rnnz, scr%rnrows, scr%nr, scr%rhop, scr%base, kpar, cfs)
            sk%nl = sk%l%ncol
            sk%nr = sk%r%ncol
            sk%havehops = 1
            sk%r%a = sk%r%a * Exp (DEF_cmplx_Ione*Dot_product(kpar, scr%rtr))
            sk%l%a = sk%l%a * Exp (DEF_cmplx_Ione*Dot_product(kpar,-scr%ltr))
         End If
         sk%alloc = 1
         sk%n = sk%c%ncol

!!$ Force to be hermitian
!!$    t_s=spherm(sk%c)
!!$    write(*,*) 'maxdev=',maxval(abs(sk%c%a-t_s%a))
!!$    sk%c%a=0.5d0*(sk%c%a+t_s%a)
!!$    call free(t_s)
      Contains

         Subroutine fillsk (sk, tnnzi, nrow, na, sites, base, kpar, cfs)
!!$ function which is actually allocate and fill sk
            Implicit None
            Integer :: tnnzi, nrow, na
            Type (zcsrmat), Target :: sk
            Type (t_str_site), Target :: sites (:)
            Real (Kind=DEF_DBL_PREC) :: base (2, 2)
            Real (Kind=DEF_DBL_PREC) :: kpar (2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs
!!$ local
            Type (t_str_site), Pointer :: jsite, isite
            Integer :: ia, nc, nr, i, nnz, j, lr, tnnz
            Integer, Pointer :: jc (:), ir (:)
            Complex (Kind=DEF_DBL_PREC), Pointer :: val (:)


            Call alloc (sk, tnnzi, nrow)

            jc => sk%jc (1:)
            ir => sk%ir (1:)

            val => sk%a (1:)
            tnnz = 0
            lr = 1
!!$   !$omp parallel do private(isite,nc,nr,j,) reduction(+:tnnz)
            Do ia = 1, na
               isite => sites (ia)
               nc = isite%ncols
               nr = isite%nrows
               nnz = fillskrows (isite, kpar, val, base, cfs)
               j = isite%srow
               ir (lr:j) = tnnz + 1
               Do i = 1, nr
                  ir (j+1) = ir (j) + nc
                  j = j + 1
                  jc (1:nc) = isite%rowind
                  jc => jc (nc+1:)
               End Do
               lr = j
               val => val (nnz+1:)
               tnnz = tnnz + nnz
            End Do
!!$  !$omp parallel end do
            sk%nnz = tnnz

         End Subroutine fillsk

         Function fillskrows (site, k, val, base, cfs_in) Result (nnz)
!!$ function which fill rows of sk which correspondent to "site"
            Implicit None
            Type (t_str_site), Target :: site
            Real (Kind=DEF_DBL_PREC) :: k (2), base (2, 2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs_in
            Complex (Kind=DEF_DBL_PREC) :: val (:)
            Integer :: nnz
!!$ local
            Integer :: in, jnc, ir, jnc1
            Complex (Kind=DEF_DBL_PREC) :: els (site%nrows, site%ncols), fact
            Integer, Pointer :: tp (:, :)
            Type (t_str_block), Pointer :: hop (:)
            Real (Kind=DEF_DBL_PREC) :: cfs !,b(2)

!!$             b(1)=0.0d0
!!$             b(2)=0.4082482904639
            cfs = 1.0d0
            If (present(cfs_in)) cfs = cfs_in

            tp => site%trpar
            hop => site%hop
!!$             write(*,*) '----', base(:,1)
            els = DEF_cmplx_zero
!!$ !$omp parallel do private(jnc,jnc1,fact) shared(cfs,k,base)
            Do in = 1, site%nhop
               jnc = site%startn (in)
               jnc1 = jnc + hop(in)%ncol - 1
!!$                write(*,*) k*(tp(1, in)*base(:, 1)+tp(2, in)*base(:, 2)+b)
               fact = cfs * Exp (DEF_cmplx_Ione*Dot_product(k, dble(tp(1, in))*base(:, 1)+dble(tp(2, &
              & in))*base(:, 2)))
               els (:, jnc:jnc1) = els (:, jnc:jnc1) + fact * hop(in)%bl
            End Do
!!$ !$omp end parallel do
            in = site%ncols
!!$ !$omp parallel do private(jnc,jnc1) firstprivate(in)
            Do ir = 1, site%nrows
               jnc = (ir-1) * in + 1
               jnc1 = jnc + in - 1
               val (jnc:jnc1) = els (ir, :)
            End Do
!!$ !$omp end parallel do
            nnz = site%nrows * in
         End Function fillskrows
      End Subroutine make_sk

      Subroutine make_sk_3D (sk, scr, kpar, kz, rz, needhops, cfs)
         Use sparselib
         Implicit None
         Type (t_strconst) :: scr
         Type (t_mathop) :: sk
         Type (zcsrmat) :: temp
         Real (Kind=DEF_DBL_PREC), Optional :: cfs
         Integer, Optional :: needhops
         Real (Kind=DEF_DBL_PREC) :: kpar (2), kz, rz(3)
!!$ local
         Integer :: nh

         If ( .Not. present(needhops)) Then
            nh = 0
         Else
            nh = needhops
         End If

         Call fillsk (sk%c, scr%nnz, scr%nrows, scr%nm, scr%main, scr%base, kpar, cfs)
         If (nh > 0) Then
            Call fillsk (sk%l, scr%lnnz, scr%lnrows, scr%nl, scr%lhop, scr%base, kpar, cfs)
            Call fillsk (sk%r, scr%rnnz, scr%rnrows, scr%nr, scr%rhop, scr%base, kpar, cfs)
            sk%nl = sk%l%ncol
            sk%nr = sk%r%ncol
            sk%havehops = 1
            sk%r%a = sk%r%a*Exp(DEF_cmplx_Ione*(Dot_product(kpar,scr%rtr) + kpar(1)*rz(1) + kpar(2)*rz(2) + kz*rz(3)))
            !sk%l%a = sk%l%a * Exp (-1*DEF_cmplx_Ione*(Dot_product(kpar, scr%ltr) + kz*rz))
            sk%l = spherm(sk%r)
         End If
         sk%alloc = 1
         sk%n = sk%c%ncol

         temp = spmatadd(sk%c, sk%r)
         sk%c = spmatadd(temp, sk%l)
         call free(temp)

      Contains

         Subroutine fillsk (sk, tnnzi, nrow, na, sites, base, kpar, cfs)
!!$ function which is actually allocate and fill sk
            Implicit None
            Integer :: tnnzi, nrow, na
            Type (zcsrmat), Target :: sk
            Type (t_str_site), Target :: sites (:)
            Real (Kind=DEF_DBL_PREC) :: base (2, 2)
            Real (Kind=DEF_DBL_PREC) :: kpar (2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs
!!$ local
            Type (t_str_site), Pointer :: jsite, isite
            Integer :: ia, nc, nr, i, nnz, j, lr, tnnz
            Integer, Pointer :: jc (:), ir (:)
            Complex (Kind=DEF_DBL_PREC), Pointer :: val (:)


            Call alloc (sk, tnnzi, nrow)

            jc => sk%jc (1:)
            ir => sk%ir (1:)

            val => sk%a (1:)
            tnnz = 0
            lr = 1
!!$   !$omp parallel do private(isite,nc,nr,j,) reduction(+:tnnz)
            Do ia = 1, na
               isite => sites (ia)
               nc = isite%ncols
               nr = isite%nrows
               nnz = fillskrows (isite, kpar, val, base, cfs)
               j = isite%srow
               ir (lr:j) = tnnz + 1
               Do i = 1, nr
                  ir (j+1) = ir (j) + nc
                  j = j + 1
                  jc (1:nc) = isite%rowind
                  jc => jc (nc+1:)
               End Do
               lr = j
               val => val (nnz+1:)
               tnnz = tnnz + nnz
            End Do
!!$  !$omp parallel end do
            sk%nnz = tnnz

         End Subroutine fillsk

         Function fillskrows (site, k, val, base, cfs_in) Result (nnz)
!!$ function which fill rows of sk which correspondent to "site"
            Implicit None
            Type (t_str_site), Target :: site
            Real (Kind=DEF_DBL_PREC) :: k (2), base (2, 2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs_in
            Complex (Kind=DEF_DBL_PREC) :: val (:)
            Integer :: nnz
!!$ local
            Integer :: in, jnc, ir, jnc1
            Complex (Kind=DEF_DBL_PREC) :: els (site%nrows, site%ncols), fact
            Integer, Pointer :: tp (:, :)
            Type (t_str_block), Pointer :: hop (:)
            Real (Kind=DEF_DBL_PREC) :: cfs !,b(2)

!!$             b(1)=0.0d0
!!$             b(2)=0.4082482904639
            cfs = 1.0d0
            If (present(cfs_in)) cfs = cfs_in

            tp => site%trpar
            hop => site%hop
!!$             write(*,*) '----', base(:,1)
            els = DEF_cmplx_zero
!!$ !$omp parallel do private(jnc,jnc1,fact) shared(cfs,k,base)
            Do in = 1, site%nhop
               jnc = site%startn (in)
               jnc1 = jnc + hop(in)%ncol - 1
!!$                write(*,*) k*(tp(1, in)*base(:, 1)+tp(2, in)*base(:, 2)+b)
               fact = cfs * Exp (DEF_cmplx_Ione*Dot_product(k, dble(tp(1, in))*base(:, 1)+dble(tp(2, &
              & in))*base(:, 2)))
               els (:, jnc:jnc1) = els (:, jnc:jnc1) + fact * hop(in)%bl
            End Do
!!$ !$omp end parallel do
            in = site%ncols
!!$ !$omp parallel do private(jnc,jnc1) firstprivate(in)
            Do ir = 1, site%nrows
               jnc = (ir-1) * in + 1
               jnc1 = jnc + in - 1
               val (jnc:jnc1) = els (ir, :)
            End Do
!!$ !$omp end parallel do
            nnz = site%nrows * in
         End Function fillskrows
      End Subroutine make_sk_3D

      Subroutine do_calc_strconst (clslist, isym)
         Implicit None
         Type (t_cluster_list) :: clslist
         Integer, Optional :: isym
!!$ local
         Integer :: ic
         Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :), gfrh (:, :, :)
         cnh => harm_nc (clslist%maxnl)
         gfrh => gaunty (clslist%maxnl)

!$omp parallel do
         Do ic = 1, clslist%ncls
            Call calc_sc_alpha (clslist%cls(ic), cnh, gfrh, isym)
         End Do
!$omp end parallel do

         Deallocate (cnh, gfrh)
      End Subroutine do_calc_strconst

      Subroutine do_calc_strconst_EMTO (clslist, avw, kap2, isym, en, ppar, iop)
         use omta_strrs
         Implicit None
         Type (t_cluster_list) :: clslist
         Real (Kind=DEF_DBL_PREC) :: avw
         Real (Kind=DEF_DBL_PREC) :: kap2
         Real (Kind=DEF_DBL_PREC) :: en 
         Integer, Optional :: isym, iop
         Type (t_omta_logder), Optional :: ppar !contains logder and hsoc, needed for downfolding
!!$ local
         Integer :: ic
!!$ emto stuff
         Real(kind=prec), allocatable :: cg    (:)
         Integer        , allocatable :: jcg   (:)
         Integer        , allocatable :: indxcg(:)
         Real(kind=prec), allocatable :: cy    (:)
         !Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :), gfrh (:, :, :)
         !cnh => harm_nc (clslist%maxnl)
         !gfrh => gaunty (clslist%maxnl)

         allocate(cg(clslist%maxnl**5/2+1),jcg(clslist%maxnl**5/2+1),cy(clslist%maxnl*clslist%maxnl*4),&
            indxcg(clslist%maxnl*clslist%maxnl*(clslist%maxnl*clslist%maxnl+1)/2+1))

         call sylmnc(cy,2*clslist%maxnl-1)
         call scg(clslist%maxnl-1,cg,indxcg,jcg)

   !!!$omp parallel do
         Do ic = 1, clslist%ncls
            if (present(ppar)) then
             Call calc_sc_a (clslist%cls(ic), avw, kap2, cg, indxcg, jcg, cy, isym, en, ppar, iop)
            else
             Call calc_sc_a (clslist%cls(ic), avw, kap2, cg, indxcg, jcg, cy, isym, en)
            endif
         End Do
   !!!!$omp end parallel do

         Deallocate (cg,jcg,cy,indxcg)


         !!$Deallocate (cnh, gfrh)
      End Subroutine do_calc_strconst_EMTO


      Subroutine calc_sc_alpha (cls, cnh, gfrh, isym_in)
!!$    performs s^alpha=((1-alpha*s)^(-1)-1)*alpha^-1 for small cluster
         Use sparselib
         Implicit None
!!$    arguments
         Type (t_cluster) :: cls
         Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :), gfrh (:, :, :)
         Integer, Optional :: isym_in
!!$     local variables
!!$    Real (kind=DEF_DBL_PREC)::sc_alpha(cls%numorb,cls%numorb)
!!$    Real (kind=DEF_DBL_PREC)::sc_can (cls%maxorb,cls%maxorb)
!!$    Real (kind=DEF_DBL_PREC)::sc_can1 (cls%maxorb,cls%maxorb)
         Real (Kind=DEF_DBL_PREC), Pointer :: sc_alpha (:, :)
         Real (KIND=DEF_DBL_PREC), Pointer :: sc_can (:, :), sc_can1 (:, :)
         Integer, Pointer :: symfield (:)

         Real (Kind=DEF_DBL_PREC) :: alpha
         Type (t_cluster_site), Pointer :: isite, jsite
         Integer :: icl, jcl, row, col, inl, jnl, l, inorb, jnorb
         Integer :: lup, ldown, nsym, srows, isym

         srows = cls%site(1)%nl ** 2
         nsym = srows * (cls%numorb-srows)
         isym = 0
         If (present(isym_in)) isym = isym_in

         Allocate (sc_alpha(cls%numorb, cls%numorb))
         Allocate (sc_can(cls%maxorb, cls%maxorb), sc_can1(cls%maxorb, cls%maxorb))
         Allocate (symfield(2*nsym))
!!$          Allocate (sc_alpha1(cls%numorb, cls%numorb))

         sc_alpha = 0.0d0

         row = 0
         Do icl = 1, cls%nsites
            isite => cls%site (icl)
            inl = isite%nl
            inorb = inl * inl
            col = row + inorb
            Do jcl = icl + 1, cls%nsites
               jsite => cls%site (jcl)
               jnl = jsite%nl
               jnorb = jnl * jnl
!!$ calc. s_{canonical}(i,j)
               Call cansc (isite, jsite, sc_can, cnh, gfrh)
               Call cansc (jsite, isite, sc_can1, cnh, gfrh)
               sc_can (1:inorb, 1:jnorb) = 0.5d0 * (sc_can(1:inorb, 1:jnorb)+transpose(sc_can1(1:jnorb, &
              & 1:inorb)))
               sc_alpha (row+1:row+inorb, col+1:col+jnorb) = sc_can (1:inorb, 1:jnorb)
               sc_alpha (col+1:col+jnorb, row+1:row+inorb) = transpose (sc_can(1:inorb, 1:jnorb))
               col = col + jnorb
            End Do
            row = row + inorb
         End Do
         If (isym == 1) Then
            Where (Abs(sc_alpha) < 1.0d-9)
               sc_alpha = 0.0d0
            End Where
            Call getscsym (nsym, sc_alpha(1:srows, (srows+1) :cls%numorb), symfield, 1.0d-9)
         End If
!!$         Write (*,*) 'trans-pre--'
!!$         Write (*, '(9(1xg22.14))') (sc_alpha(1:srows,icl), icl=1, cls%numorb)
!!$         Write (*,*) '-----'

         row = 0
         Do icl = 1, cls%nsites
            inl = cls%site(icl)%nl
            Do l = 1, inl
               alpha = cls%site(icl)%alpha(l)
               ldown = row + 1 + (l-1) * (l-1)
               lup = row + l * l
               sc_alpha (ldown:lup, :) = - alpha * sc_alpha (ldown:lup, :)
            End Do
            row = row + inl * inl
         End Do

         Do l = 1, cls%numorb
            sc_alpha (l, l) = 1.0d0
         End Do

!!$ (1-alpha*s_{can})^(-1)
         Call invert_rematrix (cls%numorb, sc_alpha, cls%numorb)

!!$ (1-alpha*s_{can})^(-1)-1
         Do l = 1, cls%numorb
            sc_alpha (l, l) = sc_alpha (l, l) - 1.0d0
         End Do

!!$ alpha^(-1)*((1-alpha*s_{can})^(-1)-1)
!!$ we do it only for first first rows
!!$ if you need it for whole s_{alpha}, then uncomment "cls%nsites" below
         row = 0
!!$          Do icl = 1, 1 !!$ cls%nsites
         Do icl = 1, cls%nsites
            inl = cls%site(icl)%nl
            Do l = 1, inl
               alpha = cls%site(icl)%alpha(l)
               ldown = row + 1 + (l-1) * (l-1)
               lup = row + l * l
               sc_alpha (ldown:lup, :) = sc_alpha (ldown:lup, :) / alpha
            End Do
            row = row + inl * inl
         End Do

#if defined(_FORCE_SC_ALPHA_HERMICITY_)
!!$          Where (Abs(sc_alpha) < 1.0d-2)
!!$             sc_alpha = 0.0d0
!!$          End Where
         sc_alpha (:, :) = .5d0 * (sc_alpha(:, :)+transpose(sc_alpha(:, :)))
#endif

         If (isym == 1) Then
            Call setscsym (nsym, sc_alpha(1:srows, (srows+1) :cls%numorb), symfield)
            Call symscons (srows*srows, sc_alpha(1:srows, 1:srows), 1.0d-5)

            Where (Abs(sc_alpha) < 1.0d-9)
               sc_alpha = 0.0d0
            End Where
         End If

         inorb = cls%site(1)%nl ** 2
         col = 0
         Do jcl = 1, cls%nsites
            jnorb = cls%site(jcl)%nl ** 2
            Allocate (cls%site(jcl)%sr(inorb, jnorb))
            Allocate (cls%site(jcl)%srd(inorb, jnorb))
            cls%site(jcl)%sr  = sc_alpha (1:inorb, col+1:col+jnorb)
            cls%site(jcl)%srd = sc_alpha (1:inorb, col+1:col+jnorb)
            col = col + jnorb
         End Do
         Deallocate (sc_alpha, sc_can, sc_can1, symfield)
      End Subroutine calc_sc_alpha

      Subroutine calc_sc_a (cls, avw, kap2, cg, indxcg, jcg, cy, isym_in, en, ppar, iop)
!!$  calculates s^a = t1/t3 + 1/t3*(-t4/t3-S0)^{-1}*1/t3*(t1*t4-t2*t3)
         Use sparselib
         Use omta_strrs
         Use omta_df
         Implicit None
!!$    arguments
         Type (t_cluster) :: cls
         !Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :), gfrh (:, :, :)
         Real (Kind=DEF_DBL_PREC) :: kap2
         Real (Kind=DEF_DBL_PREC) :: avw
         Integer, Optional :: isym_in
         Real(kind=prec), allocatable :: cg    (:)
         Integer        , allocatable :: jcg   (:)
         Integer        , allocatable :: indxcg(:)
         Real(kind=prec), allocatable :: cy    (:)
         Real(kind=prec), optional :: en
         Type (t_omta_logder), Optional :: ppar !contains logder and hsoc, needed for downfolding
         Integer, Optional :: iop
!!$     local variables
         Complex (Kind=DEF_DBL_PREC), Allocatable :: sc_a_c (:, :), sc_can (:, :), sc_can1 (:, :)
         !Real (Kind=DEF_DBL_PREC), Allocatable :: sc_a (:,:), sc_ad(:,:)
         Complex (Kind=DEF_DBL_PREC), Allocatable :: sc_a (:,:), sc_ad(:,:)
         Complex (Kind=DEF_DBL_PREC), Allocatable :: sdotc_a_c (:, :), sdotc_can (:, :), sdotc_can1 (:, :)
         Real (Kind=DEF_DBL_PREC), Allocatable :: sdotc_a (:,:)
         Integer, Pointer :: symfield (:)

         Type (t_cluster_site), Pointer :: isite, jsite
         Integer :: icl, jcl, row, col, inl, jnl, l, inorb, jnorb
         Integer :: nsym, srows, isym
         Integer :: i,j, idx
         Real (Kind=DEF_DBL_PREC) :: dt
         Real (Kind=DEF_DBL_PREC), allocatable :: t1(:), t2(:), t3(:), t4(:), a(:,:)
         Real (Kind=DEF_DBL_PREC), allocatable :: t5(:), t6(:), t7(:), t8(:)
         Real (Kind=DEF_DBL_PREC), allocatable :: alphadot(:), logder(:,:)
         Complex (Kind=DEF_DBL_PREC), allocatable :: hsoc(:,:) 
         Integer, Allocatable :: to_df(:)
         Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

         srows = cls%site(1)%nl ** 2
         nsym = srows * (cls%numorb-srows)
         isym = 0
         If (present(isym_in)) isym = isym_in

         ! for S
         Allocate (sc_a(cls%numorb, cls%numorb))
         Allocate (sc_ad(cls%numorb, cls%numorb)) !d for downspin
         Allocate (sc_a_c(cls%numorb, cls%numorb))
         Allocate (sc_can(cls%maxorb, cls%maxorb), sc_can1(cls%maxorb, cls%maxorb))
         ! for Sdot
         Allocate (sdotc_a(cls%numorb, cls%numorb))
         Allocate (sdotc_a_c(cls%numorb, cls%numorb))
         Allocate (sdotc_can(cls%maxorb, cls%maxorb), sdotc_can1(cls%maxorb, cls%maxorb))

         Allocate (symfield(2*nsym))

         Allocate (t1(cls%numorb), t2(cls%numorb), t3(cls%numorb), t4(cls%numorb), a(cls%numorb,cls%numorb))
         Allocate (t5(cls%numorb), t6(cls%numorb), t7(cls%numorb), t8(cls%numorb))
         Allocate (alphadot(cls%numorb))
         Allocate (logder(cls%numorb, cls%numorb), to_df(cls%numorb))
         !If (present(ppar)) Allocate(hsoc(2*cls%numorb,2*cls%numorb))

         sc_a = 0.0d0
         sc_ad = 0.0d0
         sc_a_c = 0.0d0
         sc_can = 0.0d0
         

         logder=0
         to_df=0

         row = 0
         idx = 1
         Do icl = 1, cls%nsites
            isite => cls%site (icl)
            inl = isite%nl
            inorb = inl * inl
            col = row + inorb

            do i=1,inorb 
            ! Ugly, I know... 
             t1(row+i) = isite%t1(ll(i)+1) 
             t2(row+i) = isite%t2(ll(i)+1) 
             t3(row+i) = isite%t3(ll(i)+1)
             t4(row+i) = isite%t4(ll(i)+1)
             t5(row+i) = isite%t5(ll(i)+1) 
             t6(row+i) = isite%t6(ll(i)+1) 
             t7(row+i) = isite%t7(ll(i)+1)
             t8(row+i) = isite%t8(ll(i)+1)
            !------------------------------
             a(row+i,row+i) = isite%a(ll(i)+1)
             alphadot(row+i) = isite%alphadot(ll(i)+1)
             if (isite%numdf > 0) then
              if (isite%idxdn(ll(i)+1) == 2) then
                to_df(idx) = row+i
                idx = idx + 1 
                logder(row+i,row+i) = isite%logder(ll(i)+1)
              endif
             endif
            enddo
           
!            if (kap2>0.d0) then
              ! for S
             !call KKRsc(isite, isite, avw, sc_can, cg, indxcg, jcg, cy, kap2, only_mstrx3=1, iop=1)
             !call KKRsc(isite, isite, avw, sc_can1, cg, indxcg, jcg, cy, kap2, only_mstrx3=1, iop=1)
             call KKRsc(isite, isite, avw, sc_can, cg, indxcg, jcg, cy, kap2, iop=1)
             !!call KKRsc(isite, isite, avw, sc_can1, cg, indxcg, jcg, cy, kap2, iop=1)
!             do i=1,inorb
!              sc_a_c (row+i,col-inorb+i) = sc_can(i,i)
!             enddo             
             !!sc_can (1:inorb, 1:inorb) = 0.5d0 * (sc_can(1:inorb, 1:inorb)+(transpose(sc_can1(1:inorb, &
             sc_can (1:inorb, 1:inorb) = 0.5d0 * (sc_can(1:inorb, 1:inorb)+(transpose(sc_can(1:inorb, &
             & 1:inorb))))
             !sc_a_c (row+1:row+inorb, row+1:row+inorb) = sc_can (1:inorb, 1:inorb)
             sc_a_c (row+1:row+inorb, row+1:row+inorb) = (transpose (sc_can(1:inorb, 1:inorb)))
             
              ! for Sdot
             call KKRsc(isite, isite, avw, sdotc_can, cg, indxcg, jcg, cy, kap2, iop=2)
             sdotc_can (1:inorb, 1:inorb) = 0.5d0 * (sdotc_can(1:inorb, 1:inorb)+transpose(sdotc_can(1:inorb,1:inorb)))
             sdotc_a_c (row+1:row+inorb, row+1:row+inorb) = (transpose (sc_can(1:inorb, 1:inorb)))


             !call KKRsc(isite, isite, avw, sdotc_can, cg, indxcg, jcg, cy, kap2, only_mstrx3=1, iop=2)
             !do i=1,inorb
             ! sdotc_a_c (row+i,col-inorb+i) = sdotc_can(i,i)
             !enddo             

!            endif

 
            Do jcl = icl + 1, cls%nsites
               jsite => cls%site (jcl)
               jnl = jsite%nl
               jnorb = jnl * jnl 
!!$ calc. s_{KKR}(i,j)
               Call KKRsc (isite, jsite, avw, sc_can, cg, indxcg, jcg, cy, kap2, iop=1)
               Call KKRsc (jsite, isite, avw, sc_can1, cg, indxcg, jcg, cy, kap2, iop=1)
               sc_can (1:inorb, 1:jnorb) = 0.5d0 * (sc_can(1:inorb, 1:jnorb)+(transpose(sc_can1(1:jnorb, &
               & 1:inorb))))
               sc_a_c (row+1:row+inorb, col+1:col+jnorb) = sc_can (1:inorb, 1:jnorb)
               sc_a_c (col+1:col+jnorb, row+1:row+inorb) = (transpose (sc_can(1:inorb, 1:jnorb)))
               
!!$ calc. sdot_{KKR}(i,j)
               Call KKRsc (isite, jsite, avw, sdotc_can, cg, indxcg, jcg, cy, kap2, iop=2)
               Call KKRsc (jsite, isite, avw, sdotc_can1, cg, indxcg, jcg, cy, kap2, iop=2)
               sdotc_can (1:inorb, 1:jnorb) = 0.5d0 * (sdotc_can(1:inorb, 1:jnorb)+(transpose(sdotc_can1(1:jnorb, &
               & 1:inorb))))
               sdotc_a_c (row+1:row+inorb, col+1:col+jnorb) = sdotc_can (1:inorb, 1:jnorb)
               sdotc_a_c (col+1:col+jnorb, row+1:row+inorb) = (transpose (sdotc_can(1:inorb, 1:jnorb)))

               col = col + jnorb
            End Do
            row = row + inorb
         End Do

        
         If (isym == 1) Then
            Where (Abs(sc_a) < 1.0d-9)
               sc_a = 0.0d0
            End Where
            Call getscsym (nsym, sc_a_c(1:srows, (srows+1) :cls%numorb), symfield, 1.0d-9)
            Call getscsym (nsym, sdotc_a_c(1:srows, (srows+1) :cls%numorb), symfield, 1.0d-9)
         End If
!!$         Write (*,*) 'trans-pre--'
!!$         Write (*, '(9(1xg22.14))') (sc_a(1:srows,icl), icl=1, cls%numorb)
!!$         Write (*,*) '-----'
 
         Do l = 1, cls%numorb
            sc_a_c (l, l) = sc_a_c(l, l) - t4(l)/t3(l)
            sdotc_a_c (l, l) = sdotc_a_c(l, l) - alphadot(l)/(t4(l)/t3(l))**2
         End Do
!!$ (-t4/t3-s_{can})^(-1)

         Call invert_Complexmatrix (cls%numorb, sc_a_c, cls%numorb, sc_a)

!!$ calculate (a^-1 - S^0)^-1 (-adot/a^2 - Sdot^0) (a^-1 - S^0)^-1
         sdotc_a = matmul(sc_a, matmul(sdotc_a_c, sc_a)) !maybe don't use matmul? We should use BLAS routine

!!$ scale up to correct s^a:
!!$ t1/t3 + 1/t3*(-t4/t3-s_{can})^(-1)*1/t3*[t1*t4-t2*t3]
         do i=1,cls%numorb
          do j=1,cls%numorb
           dt = t1(j)*t4(j) - t2(j)*t3(j)
           sc_a(i,j) = 1d0/t3(i)*sc_a(i,j)*1d0/t3(j)*dt
           ! now scale sdot:
            sdotc_a(i,j) = -1d0/t3(i)*sdotc_a(i,j)*1d0/t3(j)*dt &
                          -sc_a(i,j)*t7(j)/t3(j) - t7(i)/t3(i)*sc_a(i,j)
            if (i==j) sdotc_a(i,i) = sdotc_a(i,i) + (-t7(i)*t1(i)/t3(i) + t5(i))/t3(i)   
           ! done

           if (i==j) sc_a(i,i) = sc_a(i,i) + t1(i)/t3(i)
          enddo
         enddo

#if defined(_FORCE_SC_ALPHA_HERMICITY_)
!!$          Where (Abs(sc_a) < 1.0d-2)
!!$             sc_a = 0.0d0
!!$          End Where
!!$         sc_a (:, :) = .5d0 * (sc_a(:, :)+CONJG(transpose(sc_a(:, :))))
#endif

         !!$If (isym == 1) Then
         !!$   Call setscsym (nsym, sc_a(1:srows, (srows+1) :cls%numorb), symfield)
         !!$   Call symscons (srows*srows, sc_a(1:srows, 1:srows), 1.0d-5)

         !!$   Where (Abs(sc_a) < 1.0d-9)
         !!$      sc_a = 0.0d0
         !!$   End Where
         !!$End If

         if (present(ppar)) then
         ! DOWNFOLDING !
         sc_a = sc_a + logder
        ! If (present(ppar)) call pull_downfolded_states(hsoc, ppar%hsoc, cls%site(1)%idxdn, &
        !  &  cls%site(1)%atnum, cls%site(1)%nl**2)

         !call Downfold_iter(sc_a,to_df(1:idx-1),en)
         if (idx > 1) call Downfold_matrix(sc_a,to_df(1:idx-1)) 
        ! ==============================  This is only Lz*Sz =========================================
         endif
           sc_ad = sc_a

         !inorb = cls%site(1)%nl ** 2 - dot_product(cls%site(1)%idxdn-1, numls(1:cls%site(1)%nl))
         inorb = cls%site(1)%nl ** 2 - cls%site(1)%numdf 
         col = 0
         Do jcl = 1, cls%nsites
            !jnorb = cls%site(jcl)%nl ** 2 - dot_product(cls%site(jcl)%idxdn-1, numls(1:cls%site(jcl)%nl))
            jnorb = cls%site(jcl)%nl ** 2 - cls%site(jcl)%numdf 
            ! S
            Allocate (cls%site(jcl)%sr(inorb, jnorb))
            cls%site(jcl)%sr = sc_a (1:inorb, col+1:col+jnorb) 

            Allocate (cls%site(jcl)%srd(inorb, jnorb))
            cls%site(jcl)%srd = sc_ad (1:inorb, col+1:col+jnorb)  

            if (present(iop)) then
            if (iop == 2) then
            ! Sdot
            Allocate (cls%site(jcl)%sdotr(inorb, jnorb))
            cls%site(jcl)%sdotr = sdotc_a (1:inorb, col+1:col+jnorb)
            endif
            endif

            col = col + jnorb 
         End Do
         Deallocate (sc_a, sc_a_c, sc_can, sc_can1, symfield)
         Deallocate (sdotc_a, sdotc_a_c, sdotc_can, sdotc_can1)
         Deallocate (t1, t2, t3, t4)!, a)
         Deallocate (t5, t6, t7, t8)
         Deallocate (alphadot)
         Deallocate (logder, to_df)
         !If (present(ppar)) Deallocate (hsoc)

      end subroutine calc_sc_a 

      Function find_unique_clusters (nat, neibs) Result (clslist)
!!$ this routine find unique clasters in "neibs" - "array of neighbour lists"
!!$ by the way it will set to all "neibs(i)" index of correspondent cluster
         Implicit None
         Type (t_neighbour_list) :: neibs (nat)
         Type (t_cluster_list) :: clslist
         Integer :: nat
!!$ local
         Type (t_cluster), Pointer :: t_cls
         Integer :: i, j, isold

         clslist%ncls = 0
!!$omp parallel&
!!$omp& private(isold,t_cls,j)&
!!$omp& shared(clslist,neibs)
!!$omp& shared(nat)&
!!$omp do
         Do i = 1, nat
            isold = 0
            t_cls => clslist%first
            Do j = 1, clslist%ncls
!!$ check, if we have already this cluster
               If (is_same_cluster(t_cls, neibs(i), j) /= 0) Then
                  isold = 1
                  Exit
               End If
               t_cls => t_cls%next
            End Do

            If (isold == 0) Then
!!$ add new one
!!$omp critical
               Call add_cluster (clslist, neibs(i))
!!$omp end critical	         	
            End If
         End Do
!!$omp end do
!!$omp end parallel


         Allocate (clslist%cls(clslist%ncls))
         clslist%cls (1) = clslist%first
         Do i = 2, clslist%ncls
            clslist%cls (i) = clslist%cls(i-1)%next
         End Do

         Return

      Contains
         Subroutine add_cluster (clslist, neib)
            Implicit None
            Type (t_cluster_list) :: clslist
            Type (t_neighbour_list) :: neib
!!$ local
            Type (t_cluster), Pointer :: p_cluster
            Integer :: nn, i, mnl, norb, nl
            Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

!!$ If we want to have sorted cluster
!!$             Call sort_neibs_pos(neib)
            Allocate (p_cluster)
            If (clslist%ncls /= 0) Then
               clslist%last%next => p_cluster
            Else
               clslist%first => p_cluster
            End If
            clslist%last => p_cluster
            clslist%ncls = clslist%ncls + 1

            nn = neib%num
            Allocate (p_cluster%site(nn))

            p_cluster%nsites = nn
            p_cluster%nl = neib%nbr(1)%nl

            Allocate (neib%cl_dir_idx(nn))
            norb = 0
            mnl = 0
            Do i = 1, nn
               nl = neib%nbr(i)%nl
               p_cluster%site(i)%nl = nl
               norb = norb + nl * nl
               mnl = Max (nl, mnl)
               p_cluster%site(i)%alpha => neib%nbr(i)%alpha
               p_cluster%site(i)%rvec = neib%nbr(i)%rvec
               p_cluster%site(i)%scale = neib%nbr(i)%scale
               p_cluster%site(i)%nrm = 0

               ! ugly workaround
               p_cluster%site(i)%atnum = neib%nbr(i)%atnum
               p_cluster%site(i)%t1 => neib%nbr(i)%t1
               p_cluster%site(i)%t2 => neib%nbr(i)%t2
               p_cluster%site(i)%t3 => neib%nbr(i)%t3
               p_cluster%site(i)%t4 => neib%nbr(i)%t4
               p_cluster%site(i)%t5 => neib%nbr(i)%t5
               p_cluster%site(i)%t6 => neib%nbr(i)%t6
               p_cluster%site(i)%t7 => neib%nbr(i)%t7
               p_cluster%site(i)%t8 => neib%nbr(i)%t8
               ! -----
               p_cluster%site(i)%a => neib%nbr(i)%a
               p_cluster%site(i)%alphadot => neib%nbr(i)%alphadot
               !downfolding stuff
               p_cluster%site(i)%logder => neib%nbr(i)%logder
               if (neib%nbr(i)%numdf > 0) then
               p_cluster%site(i)%idxdn => neib%nbr(i)%idxdn
               p_cluster%site(i)%numdf = neib%nbr(i)%numdf
               endif

               neib%cl_dir_idx (i) = i
            End Do
            clslist%maxnl = Max (mnl, clslist%maxnl)
            p_cluster%numorb = norb
            p_cluster%numnrm = 0
            p_cluster%maxorb = mnl ** 2
            Nullify (p_cluster)
            neib%cl_idx = clslist%ncls
         End Subroutine add_cluster

         Function is_same_cluster (cluster, neib, ncl) Result (ret)
            Implicit None
            Type (t_cluster) :: cluster
            Type (t_neighbour_list) :: neib
            Integer :: ret, ncl
!!$ local
            Integer :: iseq, ic, in, nl, ndir (neib%num)
            ret = 0
            ndir = 0
            iseq = cluster%nsites
            If ((iseq == neib%num)) Then
               Do ic = 1, cluster%nsites
                  Do in = 1, neib%num
                     If (sum(Abs(cluster%site(ic)%rvec-neib%nbr(in)%rvec)) < E_SAME_COORD) Then
                        nl = cluster%site(ic)%nl - cluster%site(ic)%numdf 
                        If (nl == neib%nbr(in)%nl - neib%nbr(in)%numdf ) Then
                           If (sum(Abs(cluster%site(ic)%alpha-neib%nbr(in)%alpha)) < E_SAME_COORD) Then
                              If (cluster%site(ic)%scale == neib%nbr(in)%scale) Then
                                 iseq = iseq - 1
                                 ndir (in) = ic
                              End If
                           End If
                        End If
                     End If
                  End Do
               End Do
            End If
            If (iseq == 0) Then
               ret = 1
               Allocate (neib%cl_dir_idx(neib%num))
               neib%cl_dir_idx = ndir
               neib%cl_idx = ncl
            End If
         End Function is_same_cluster
      End Function find_unique_clusters


      Subroutine free_cluster_list (clslist,iop)
!!$ free cluster list
         Implicit None
         Type (t_cluster_list) :: clslist
         Integer :: iop
!!$ local
         Integer :: i, j
         Type (t_cluster), Pointer :: cls_ptr, cls_ptr1
         cls_ptr => clslist%first
         Do i = 1, clslist%ncls
            Do j = 1, clslist%cls(i)%nsites
               Deallocate (clslist%cls(i)%site(j)%sr)
               Deallocate (clslist%cls(i)%site(j)%srd)
               If (iop==2) Deallocate (clslist%cls(i)%site(j)%sdotr)
            End Do
            Deallocate (clslist%cls(i)%site)
            cls_ptr1 => cls_ptr%next
            Deallocate (cls_ptr)
            cls_ptr => cls_ptr1
         End Do
         Nullify (clslist%first, clslist%last)
         Deallocate (clslist%cls)
         clslist%ncls = 0
      End Subroutine free_cluster_list

!!$$ EMTO ADDITION: FIND_NEIBS_EMTO !!$$

      Function find_neibs_EMTO (geom, cr) Result (neibs)
!!$ this subroutine find neighbours on distance less that "cutrat"
         Use geometry_module
         Implicit None
         Type (t_geometry_EMTO) :: geom
         Type (t_neighbour_list), Pointer :: neibs (:)
         Real (Kind=DEF_DBL_PREC), Optional :: cr
!!$ local
         Integer, Parameter :: initial_maxneib = 50
         Integer :: i, jx, jy, nat, k, maxneib, t, tb, te, skip
         Real (Kind=DEF_DBL_PREC) :: dist, base (2, 2), cutrat, wsr1, mdist, mdst1
         Real (Kind=DEF_DBL_PREC) :: trans (3,-1:1), scale
         Real (Kind=DEF_DBL_PREC), Pointer :: coord (:), coord0 (:), dir (:)
         Type (tsite), Pointer :: atpi, atpk
!!$          Type (t_neighbour), Pointer :: t_neib
         Integer, Pointer :: slice (:)
         Integer :: slsz, ksl
         Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

         scale = geom%scale 
         maxneib = initial_maxneib
         nat = geom%num
!!$ allocate memory for clusters
         Allocate (neibs(nat))
    
!!$ put atoms itself as a first neighbours
    
         Allocate (coord(3), coord0(3), slice(nat))
         coord (:) = 0.0d0
    
!!$ !$omp parallel default(shared) firstprivate(scale,nat)
!$omp parallel do private(i) default(shared)
         Do i = 1, nat
            Call add_neib_EMTO (neibs(i), geom%atoms(i), i, coord, scale, 0, 0, 0)
            neibs(i)%basenl = geom%atoms(i)%ptr%lmx+1  
         End Do
!$omp end parallel do

         cutrat = geom%cutrat * scale 
         If (present(cr)) cutrat = cr * scale
         base = geom%base
         trans = 0.0d0
         trans (:,-1) = - geom%l_perp_trans
         trans (:, 1) = geom%r_perp_trans

!!$ loop over atoms.
         Do i = 1, nat

            atpi => geom%atoms (i)
            !wsr1 = atpi%ptr%s 
            wsr1 = geom%rawsr 
            tb = 0
            te = 0
            If (i <= geom%l_transnum) Then
               tb = - 1
            End If
            If (i > (nat-geom%r_transnum)) Then
               te = 1
            End If

            !mdist = cutrat * Sqrt (wsr1) * 2.0d0
            mdist = cutrat * Sqrt (geom%rawsr) * 2.0d0
            !mdist = cutrat * 2.0d0 * Sqrt((3/(16*DEF_M_PI))**(1.d0/3))
            slsz = 0
            Do k = 1, nat
               atpk => geom%atoms (k)
               !If (Abs(atpk%coord(3)-atpi%coord(3)) > mdist*Sqrt(atpk%ptr%s)) Cycle
               If (Abs(atpk%coord(3)-atpi%coord(3)) > mdist*Sqrt(geom%rawsr)) Cycle
               !If (Abs(atpk%coord(3)-atpi%coord(3)) > &
               !     & mdist*Sqrt((3/(16*DEF_M_PI))**(1.d0/3))) Cycle
               slsz = slsz + 1
               slice (slsz) = k
            End Do
            !mdst1 = cutrat * Sqrt (wsr1)
            mdst1 = cutrat * Sqrt (geom%rawsr)
            Do t = tb, te
               coord0 (:) = atpi%coord + trans (:, t)
               coord (:) = coord0 (:)
!!$ loops over parallel translation indexes (jx,jy)

               Do jy = - geom%ntrpar, geom%ntrpar
                  Do jx = - geom%ntrpar, geom%ntrpar
                     skip = Abs (jx) + Abs (jy) + Abs (t)
                     coord (1:2) = coord0 (1:2) + jx * base (:, 1) + jy * base (:, 2)
!!$ loop over atoms again, if atom wchich defined by three previous
!loops-indexes
!!$ is in neighbourhood, then we will add it.

!!$  !!$omp parallel default(shared)&
!!$  !!$omp& private(atpk,dist,dir,mdist,k)
                     Allocate (dir(3))
!$omp do ordered
                     Do ksl = 1, slsz
                        k = slice (ksl)
                        If (skip == 0 .And. k == i) Cycle
                        atpk => geom%atoms (k)
                        !mdist = mdst1 * Sqrt (atpk%ptr%s )
                        mdist = mdst1 * Sqrt (geom%rawsr)
                        If (dabs(atpk%coord(2)-coord(2)) > mdist .Or. dabs(atpk%coord(1)-coord(1)) > mdist &
                       & .Or. dabs(atpk%coord(3)-coord(3)) > mdist) Cycle
                        dist = get_dist (atpk%coord, coord, dir)
!$omp ordered
                        If (dist < mdist) Then
                           Call add_neib_EMTO (neibs(k), atpi, i, dir, scale, t, jx, jy)
                        End If
!$omp end ordered
                     End Do
!!$ !!$omp end do
                     Deallocate (dir)
!!$ !!$omp end parallel
                  End Do
               End Do

            End Do
         End Do

         Deallocate (coord, coord0, slice)

!!$ convert list to array
         Do i = 1, nat
            Allocate (neibs(i)%nbr(neibs(i)%num))
            neibs(i)%nbr(1) = neibs(i)%first
            Do k = 2, neibs(i)%num
               neibs(i)%nbr(k) = neibs(i)%nbr(k-1)%next
            End Do
!!$             Call sort_neibs_pos(neibs(i))
         End Do

         Return

      Contains

         Subroutine add_neib_EMTO (neibs, atp, n, dir, scale, t, ix, iy)
            Implicit None
            Type (tsite) :: atp
            Type (t_neighbour_list) :: neibs
            Integer :: n, t, ix, iy
            Real (Kind=DEF_DBL_PREC) :: dir (3), scale
!!$ local
            Type (t_neighbour), Pointer :: t_neib
            Integer :: i
            Integer :: numls(15) = [1,3,5,7,9, 11,13,15,17,19, 21,23,25,27,29]

            Allocate (t_neib)
            t_neib%nl = atp%ptr%lmx+1
            t_neib%atnum = n
            t_neib%nrm = 0
            t_neib%alpha => atp%ptr%alpha
            t_neib%rvec = dir
            t_neib%trpar = (/ ix, iy /)
            t_neib%trperp = t

            ! This is nasty code... but i can't make a rank-2 vector and do the same thing
            ! I get a weird error about 'rank remapping target must be rank 1 or simply contiguous'
            ! So i just make 8 rank-1 vectors... stupid but it works
            t_neib%t1(1:atp%ptr%lmx+1) => atp%ptr%t(1,0:atp%ptr%lmx)
            t_neib%t2(1:atp%ptr%lmx+1) => atp%ptr%t(2,0:atp%ptr%lmx)
            t_neib%t3(1:atp%ptr%lmx+1) => atp%ptr%t(3,0:atp%ptr%lmx)
            t_neib%t4(1:atp%ptr%lmx+1) => atp%ptr%t(4,0:atp%ptr%lmx)
            t_neib%t5(1:atp%ptr%lmx+1) => atp%ptr%t(5,0:atp%ptr%lmx)
            t_neib%t6(1:atp%ptr%lmx+1) => atp%ptr%t(6,0:atp%ptr%lmx)
            t_neib%t7(1:atp%ptr%lmx+1) => atp%ptr%t(7,0:atp%ptr%lmx)
            t_neib%t8(1:atp%ptr%lmx+1) => atp%ptr%t(8,0:atp%ptr%lmx)
            t_neib%a(1:atp%ptr%lmx+1) => atp%ptr%a(0:atp%ptr%lmx)
            t_neib%alphadot(1:atp%ptr%lmx+1) => atp%ptr%alphadot(0:atp%ptr%lmx)
            ! Downfolding stuff
            t_neib%logder(1:atp%ptr%lmx+1) => atp%ptr%pp(1,0:atp%ptr%lmx,1) !note that it's the first spin
            t_neib%numdf = dot_product(atp%ptr%idxdn(0:atp%ptr%lmx)-1, numls(1:atp%ptr%lmx+1))
            if (t_neib%numdf > 0) then 
             t_neib%idxdn(1:atp%ptr%lmx+1)  => atp%ptr%idxdn(0:atp%ptr%lmx) !to df or not to df
            endif 
 
            t_neib%scale = atp%ptr%s * scale
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
         End Subroutine add_neib_EMTO

         Function get_dist (a, b, r) Result (dist)
            Implicit None
            Real (Kind=DEF_DBL_PREC) :: a (3), b (3), r (3), dist
            r = b - a
            dist = Sqrt (sum(r*r))
            Return
         End Function get_dist
      End Function find_neibs_EMTO



!!$$ END OF EMTO ADDITION 'FIND_NEIBS_EMTO' !!$$
   
  Subroutine downfold_logder(geom, logder)
    ! pull the non-downfolded part out of the 'logder' matrix
     Use geometry_module 
     Type(t_geometry_EMTO) :: geom
     Type (zcsrmat)  :: logder 
     !local
     integer :: i, j, jn, df, df_c
     Type (zcsrmat)  :: new_logder 
     logical :: flag
     
     call alloc(new_logder, logder%nnz-size(geom%to_downfold), logder%ncol-size(geom%to_downfold))

     df = 0
     jn=1
     new_logder%ir(1) = 1
     new_logder%nnz = 0
     do i=1,logder%nrow
       flag=.false.
       if ( ANY( geom%to_downfold == i)) then
          flag=.true. 
        else 
          new_logder%ir(i-df+1) = new_logder%ir(i-df)
        do j=logder%ir(i),logder%ir(i+1)-1
          if ( ANY( geom%to_downfold == logder%jc(j))) then
          else
           new_logder%a(jn) = logder%a(j)
           new_logder%jc(jn) = logder%jc(j) - COUNT(geom%to_downfold < logder%jc(j)) 
           new_logder%ir(i-df+1) = new_logder%ir(i-df+1) + 1
           new_logder%nnz = new_logder%nnz+1
           jn = jn+1
          endif
         enddo
        endif
        if (flag) df = df+1
     enddo
     

     call spcopy(new_logder, logder)
     call free(new_logder)

   End Subroutine downfold_logder

!   Subroutine pull_downfolded_states(densmat, amat, idxdn, atnum, norb)
!     !pulls matrix elements from amat into densmat according to idxdn and atnum
!     !meant for spin-orbit coupling part of the hamiltonian, but could be usefule elsewhere
!     ! passed
!     Complex (Kind=DEF_DBL_PREC) :: densmat(:,:)
!     Type (zcsrmat) :: amat
!     Integer :: idxdn(4)
!     Integer :: atnum, norb
!     ! local
!     Integer :: isp,i,j,ir,jc,nz
!     Integer :: irow, icol
!
!     densmat = 0
!
!     !we have to loop over the elements of amat which correspond to atnum
!     !and put them into the densmat matrix, but only if idxdn==2
!     do isp=1,2
!      do i=1,norb
!       ir = i + (atnum-1)*norb*2 + norb*(isp-1) !shift to find the correct place in the large matrix
!       do nz=amat%ir(ir), amat%ir(ir+1)-1 !select elements of amat
!        jc = amat%jc(nz) !column index
!        j = jc - (atnum-1)*norb*2 ! shift to put it in the correct place in densmat
!
!        irow = i+(isp-1)*norb
!        icol = j
!        if ( idxdn(ll(mod(irow-1,norb)+1)+1) == 2 .or. idxdn(ll(mod(icol-1,norb)+1)+1) == 2) then
!         densmat(irow,icol) = amat%a(nz)
!        endif
!       enddo
!      enddo
!     enddo
!
!   end subroutine pull_downfolded_states

      Function find_neibs (geom, cr) Result (neibs)
!!$ this subroutine find neighbours on distance less that "cutrat"
         Use geometry_module
         Implicit None
         Type (t_geometry) :: geom
         Type (t_neighbour_list), Pointer :: neibs (:)
         Real (Kind=DEF_DBL_PREC), Optional :: cr
!!$ local
         Integer, Parameter :: initial_maxneib = 50
         Integer :: i, jx, jy, nat, k, maxneib, t, tb, te, skip
         Real (Kind=DEF_DBL_PREC) :: dist, base (2, 2), cutrat, wsr1, mdist, mdst1
         Real (Kind=DEF_DBL_PREC) :: trans (3,-1:1), scale
         Real (Kind=DEF_DBL_PREC), Pointer :: coord (:), coord0 (:), dir (:)
         Type (t_atom_pointer), Pointer :: atpi, atpk
!!$          Type (t_neighbour), Pointer :: t_neib
         Integer, Pointer :: slice (:)
         Integer :: slsz, ksl
         scale = geom%scale

         maxneib = initial_maxneib
         nat = geom%num
!!$ allocate memory for clusters
         Allocate (neibs(nat))

!!$ put atoms itself as a first neighbours

         Allocate (coord(3), coord0(3), slice(nat))
         coord (:) = 0.0d0

!!$ !$omp parallel default(shared) firstprivate(scale,nat)
!$omp parallel do private(i) default(shared)
         Do i = 1, nat
            Call add_neib (neibs(i), geom%atoms(i), i, coord, scale, 0, 0, 0)
            neibs(i)%basenl = geom%atoms(i)%ptr%nl
         End Do
!$omp end parallel do

         cutrat = geom%cutrat * scale
         If (present(cr)) cutrat = cr * scale
         base = geom%base
         trans = 0.0d0
         trans (:,-1) = - geom%l_perp_trans
         trans (:, 1) = geom%r_perp_trans

!!$ loop over atoms.
         Do i = 1, nat

            atpi => geom%atoms (i)
            wsr1 = atpi%ptr%wsr
            tb = 0
            te = 0
            If (i <= geom%l_transnum) Then
               tb = - 1
            End If
            If (i > (nat-geom%r_transnum)) Then
               te = 1
            End If

            mdist = cutrat * Sqrt (wsr1) * 2.0d0
            slsz = 0
            Do k = 1, nat
               atpk => geom%atoms (k)
               If (Abs(atpk%coord(3)-atpi%coord(3)) > mdist*Sqrt(atpk%ptr%wsr)) Cycle
               slsz = slsz + 1
               slice (slsz) = k
            End Do
            mdst1 = cutrat * Sqrt (wsr1)
            Do t = tb, te
               coord0 (:) = atpi%coord + trans (:, t)
               coord (:) = coord0 (:)
!!$ loops over parallel translation indexes (jx,jy)

               Do jy = - geom%ntrpar, geom%ntrpar
                  Do jx = - geom%ntrpar, geom%ntrpar
                     skip = Abs (jx) + Abs (jy) + Abs (t)
                     coord (1:2) = coord0 (1:2) + jx * base (:, 1) + jy * base (:, 2)
!!$ loop over atoms again, if atom wchich defined by three previous loops-indexes
!!$ is in neighbourhood, then we will add it.

!!$  !!$omp parallel default(shared)&
!!$  !!$omp& private(atpk,dist,dir,mdist,k)
                     Allocate (dir(3))
!$omp do ordered
                     Do ksl = 1, slsz
                        k = slice (ksl)
                        If (skip == 0 .And. k == i) Cycle
                        atpk => geom%atoms (k)
                        mdist = mdst1 * Sqrt (atpk%ptr%wsr)
                        If (dabs(atpk%coord(2)-coord(2)) > mdist .Or. dabs(atpk%coord(1)-coord(1)) > mdist &
                       & .Or. dabs(atpk%coord(3)-coord(3)) > mdist) Cycle
                        dist = get_dist (atpk%coord, coord, dir)
!$omp ordered
                        If (dist < mdist) Then
                           Call add_neib (neibs(k), atpi, i, dir, scale, t, jx, jy)
                        End If
!$omp end ordered
                     End Do
!!$ !!$omp end do
                     Deallocate (dir)
!!$ !!$omp end parallel
                  End Do
               End Do

            End Do
         End Do

         Deallocate (coord, coord0, slice)

!!$ convert list to array
         Do i = 1, nat
            Allocate (neibs(i)%nbr(neibs(i)%num))
            neibs(i)%nbr(1) = neibs(i)%first
            Do k = 2, neibs(i)%num
               neibs(i)%nbr(k) = neibs(i)%nbr(k-1)%next
            End Do
!!$             Call sort_neibs_pos(neibs(i))
         End Do

         Return

      Contains

         Subroutine add_neib (neibs, atp, n, dir, scale, t, ix, iy)
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
            t_neib%nrm = 0
            t_neib%alpha => atp%ptr%alpha
            t_neib%rvec = dir
            t_neib%trpar = (/ ix, iy /)
            t_neib%trperp = t
            t_neib%numdf = 0

            t_neib%scale = scale * atp%ptr%wsr
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
      End Function find_neibs


      Subroutine free_neibs (neibs)
         Implicit None
         Type (t_neighbour_list), Pointer :: neibs (:)
!!$ local
         Integer :: i, j
         Type (t_neighbour), Pointer :: nbr_ptr, nbr_ptr1


         Do i = 1, size (neibs)
            nbr_ptr => neibs(i)%first
            Do j = 1, neibs(i)%num
               nbr_ptr1 => nbr_ptr%next
               Deallocate (nbr_ptr)
               nbr_ptr => nbr_ptr1
            End Do
            Deallocate (neibs(i)%nbr)
            Deallocate (neibs(i)%cl_dir_idx)
!!$            Deallocate (neibs(i)%cl_dir_idx)
            neibs(i)%num = 0
         End Do

      End Subroutine free_neibs


      Subroutine invert_rematrix (n, a, lda)
!!$    a wrapper for the lapack routines for inverting the real
!!$    square matrix
         Implicit None
         External dgetrf, dgetri
         
         Integer, Intent (In) :: n, lda
         Real (Kind=DEF_DBL_PREC), Intent (Inout) :: a (:, :)

         Real (Kind=DEF_DBL_PREC), Pointer :: work (:)
         Integer :: info
         Integer, Pointer :: ipiv (:)

         Allocate (ipiv(n), work(n))

!!$ DAMMMIT! Valgrind complains here!
         Call dgetrf (n, n, a, lda, ipiv, info)! computes lu factorization (lapack
!!$


         If (info /= 0) Then
            Write (*,*) 'problems encountered while in dgetrf routine, info =', info
            Stop
         End If
         Call dgetri (n, a, lda, ipiv, work, n, info)! inverses the a matrix
         If (info /= 0) Then
            Write (*,*) 'problems encountered while in dgetrfi routine, info =', info
            Stop
         End If
         Deallocate (ipiv, work)
      End Subroutine invert_rematrix

      Subroutine invert_Complexmatrix (n, a, lda, aout)
!!$    a wrapper for the lapack routines for inverting a complex
!!$    square matrix
         Implicit None
         External zgetrf, zgetri

         Integer, Intent(In) :: n, lda
         Complex (Kind=DEF_DBL_PREC), Intent (In) :: a (:, :)
         Complex (Kind=DEF_DBL_PREC), Intent (Out) :: aout (:, :)

         Complex (Kind=DEF_DBL_PREC), Allocatable :: work (:)
         Integer :: info
         Integer, Allocatable :: ipiv (:)

         Allocate (ipiv(n), work(n))

!!$ DAMMMIT! Valgrind complains here!
         Call zgetrf (n, n, a, lda, ipiv, info)! computes lu factorization (lapack

         If (info /= 0) Then
            Write (*,*) 'problems encountered while in zgetrf routine, info =', info
            Stop
         End If
         Call zgetri (n, a, lda, ipiv, work, n, info)! inverses the a matrix
         If (info /= 0) Then
            Write (*,*) 'problems encountered while in zgetri routine, info =', info
            Stop
         End If
         Deallocate (ipiv, work)
        
         aout = REAL(a)

      End Subroutine invert_Complexmatrix

!!$ Bellow two modified QSort routines
!!$ Can somebody tell me finally when qsort will be standart in fortran and
!!$ i will not need to include tens of qsort functions for different datatypes
!!$ into the code?????

      Recursive Subroutine sort_neibs_pos (neibs)
!!$ sort sites in t_neighbour_list accroding to their coordinates
!!$ priority of coordinates is (z,y,x)
!!$ i beleive that this is increase accuracy
!!$ it uses  modified quicksort
         Implicit None
         Type (t_neighbour_list) :: neibs
!!$ local
         Type (t_neighbour), Pointer :: neib_p, dneib (:)
         Integer :: n
         Integer :: Index (neibs%num-1)
         Integer :: lstk (31), rstk (31), istk
         Integer :: l, r, i, j, p, indexp, indext, fl1
         Integer, Parameter :: m = 9

         n = neibs%num - 1
         dneib => neibs%nbr (2:neibs%num)
         Do i = 1, n
            Index (i) = i
         End Do
         If (n .Le. m) Go To 900

!!$ qsort itself
         istk = 0
         l = 1
         r = n
200      Continue
         i = l
         j = r
         p = (l+r) / 2
         indexp = Index (p)
         neib_p => dneib (indexp)
         If (cmpr_fun(dneib(Index(l)), neib_p)) Then
            Index (p) = Index (l)
            Index (l) = indexp
            indexp = Index (p)
            neib_p => dneib (indexp)
         End If
         If (cmpr_fun(neib_p, dneib(Index(r)))) Then
            If (cmpr_fun(dneib(Index(l)), dneib(Index(r)))) Then
               Index (p) = Index (l)
               Index (l) = Index (r)
            Else
               Index (p) = Index (r)
            End If
            Index (r) = indexp
            indexp = Index (p)
            neib_p => dneib (indexp)
         End If
300      Continue
         i = i + 1
         If (cmpr_fun(neib_p, dneib(Index(i)))) Go To 300
400      Continue
         j = j - 1
         If (cmpr_fun(dneib(Index(j)), neib_p)) Go To 400
         If (i .Lt. j) Then
            indext = Index (i)
            Index (i) = Index (j)
            Index (j) = indext
            Go To 300
         Else
            If (r-j .Ge. i-l .And. i-l .Gt. m) Then
               istk = istk + 1
               lstk (istk) = j + 1
               rstk (istk) = r
               r = i - 1
            Else If (i-l .Gt. r-j .And. r-j .Gt. m) Then
               istk = istk + 1
               lstk (istk) = l
               rstk (istk) = i - 1
               l = j + 1
            Else If (r-j .Gt. m) Then
               l = j + 1
            Else If (i-l .Gt. m) Then
               r = i - 1
            Else
               If (istk .Lt. 1) Go To 900
               l = lstk (istk)
               r = rstk (istk)
               istk = istk - 1
            End If
            Go To 200
         End If

900      Continue
!!$ "buble sort". used if we have datasize<9 and for last iteration after "qsort"
         Do i = 2, n
            If (cmpr_fun(dneib(Index(i-1)), dneib(Index(i)))) Then
               indexp = Index (i)
               neib_p => dneib (indexp)
               p = i - 1
               fl1 = 0
               Do While (fl1 ==  0)
                  Index (p+1) = Index (p)
                  p = p - 1
                  If (p .Le. 0) Then
                     fl1 = 1
                  Else
                     If (cmpr_fun(neib_p, dneib(Index(p)))) fl1 = 1
                  End If
               End Do
               Index (p+1) = indexp
            End If
         End Do
         dneib = dneib (Index)

         If ((-dneib(1)%rvec(3)) > (dneib(n)%rvec(3))) Then
            Write (*,*) 'inversion!'
            Do i = 1, n
               Index (i) = n - i + 1
               dneib(i)%rvec(:) = - dneib(i)%rvec(:)
            End Do
            dneib = dneib (Index)
            Call sort_neibs_pos (neibs)
         End If

         Nullify (dneib, neib_p)
      Contains
         Function cmpr_fun (nb1, nb2) Result (res)
!!$ helper function. compare two datatypes (positions of sites here)
            Implicit None
            Type (t_neighbour) :: nb1, nb2
            Logical :: res
!!$ local
            Real (Kind=DEF_DBL_PREC) :: r (3)
            res = .False.
            r = nb1%rvec - nb2%rvec

            If (r(3) == 0.0d0) Then
               If (r(2) == 0.0d0) Then
                  If (r(1) > 0.0d0) Then
                     res = .True.
                     Return
                  End If
               Else
                  If (r(2) > 0.0d0) Then
                     res = .True.
                     Return
                  End If
               End If
            Else
               If (r(3) > 0.0d0) Then
                  res = .True.
                  Return
               End If
            End If
         End Function cmpr_fun
      End Subroutine sort_neibs_pos

      Subroutine sort_neibs_num (neibs)
!!$ sort sites in t_neighbour_list accroding to their number
!!$ it uses  modified quicksort
         Implicit None
         Type (t_neighbour_list) :: neibs
!!$ local
         Type (t_neighbour), Pointer :: neib_p, dneib (:)
         Integer :: n
         Integer :: Index (neibs%num)
         Integer :: lstk (31), rstk (31), istk
         Integer :: l, r, i, j, p, indexp, indext, fl1
         Integer, Parameter :: m = 9

         n = neibs%num
         dneib => neibs%nbr
         Do i = 1, n
            Index (i) = i
         End Do
         If (n .Le. m) Go To 900

!!$ qsort itself
         istk = 0
         l = 1
         r = n
200      Continue
         i = l
         j = r
         p = (l+r) / 2
         indexp = Index (p)
         neib_p => dneib (indexp)
         If (cmpr_fun(dneib(Index(l)), neib_p)) Then
            Index (p) = Index (l)
            Index (l) = indexp
            indexp = Index (p)
            neib_p => dneib (indexp)
         End If
         If (cmpr_fun(neib_p, dneib(Index(r)))) Then
            If (cmpr_fun(dneib(Index(l)), dneib(Index(r)))) Then
               Index (p) = Index (l)
               Index (l) = Index (r)
            Else
               Index (p) = Index (r)
            End If
            Index (r) = indexp
            indexp = Index (p)
            neib_p => dneib (indexp)
         End If
300      Continue
         i = i + 1
         If (cmpr_fun(neib_p, dneib(Index(i)))) Go To 300
400      Continue
         j = j - 1
         If (cmpr_fun(dneib(Index(j)), neib_p)) Go To 400
         If (i .Lt. j) Then
            indext = Index (i)
            Index (i) = Index (j)
            Index (j) = indext
            Go To 300
         Else
            If (r-j .Ge. i-l .And. i-l .Gt. m) Then
               istk = istk + 1
               lstk (istk) = j + 1
               rstk (istk) = r
               r = i - 1
            Else If (i-l .Gt. r-j .And. r-j .Gt. m) Then
               istk = istk + 1
               lstk (istk) = l
               rstk (istk) = i - 1
               l = j + 1
            Else If (r-j .Gt. m) Then
               l = j + 1
            Else If (i-l .Gt. m) Then
               r = i - 1
            Else
               If (istk .Lt. 1) Go To 900
               l = lstk (istk)
               r = rstk (istk)
               istk = istk - 1
            End If
            Go To 200
         End If

900      Continue
!!$ "buble sort". used if we have datasize<9 and for last iteration after "qsort"
         Do i = 2, n
            If (cmpr_fun(dneib(Index(i-1)), dneib(Index(i)))) Then
               indexp = Index (i)
               neib_p => dneib (indexp)
               p = i - 1
               fl1 = 0
               Do While (fl1 ==  0)
                  Index (p+1) = Index (p)
                  p = p - 1
                  If (p .Le. 0) Then
                     fl1 = 1
                  Else
                     If (cmpr_fun(neib_p, dneib(Index(p)))) fl1 = 1
                  End If
               End Do
               Index (p+1) = indexp
            End If
         End Do
         dneib = dneib (Index)
         neibs%cl_dir_idx = neibs%cl_dir_idx(Index)
         Nullify (dneib, neib_p)
      Contains
         Function cmpr_fun (nb1, nb2) Result (res)
!!$ helper function. compare two datatypes (positions of sites here)
            Implicit None
            Type (t_neighbour) :: nb1, nb2
            Logical :: res
!!$ local
            If (nb1%atnum > nb2%atnum) Then
               res = .True.
            Else
               res = .False.
            End If
         End Function cmpr_fun
      End Subroutine sort_neibs_num


!!$ heritage from Ilja bellow, fixed slightly

      Subroutine cansc (isite, jsite, sc, cnh, gfrh)
!!$---------------------------------------------------------------
!!$   canonical structure constant connecting 2 sites "isite" and "jsite"
!!$       sc(i1,i2) - the structure constant matrix
!!$---------------------------------------------------------------
!!$   remark: "cnh" should be calculated before by "harm_nc" function
!!$           "gfrh" should be calculated before by "gaunty" function
!!$---------------------------------------------------------------
         Implicit None
         Type (t_cluster_site) :: isite, jsite
         Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :), gfrh (:, :, :)
         Real (Kind=DEF_DBL_PREC) :: sc (:, :)
!!$ local
!!$         Real (Kind=DEF_DBL_PREC) :: yps ((isite%nl+jsite%nl-1)**2)
         Real (Kind=DEF_DBL_PREC), Pointer :: yps (:)
         Real (Kind=DEF_DBL_PREC) :: tx, ty, tz, r, rmin, w1, w2
         Integer :: nl1, nl2, l, i2, l2, i1, l1, ista, ifin, m1, i, m2
         Integer :: lmaxh, lmax, lmaxh1, lmax2, lmax1, lmaxh2
         Real (Kind=DEF_DBL_PREC) :: dfac (0:(isite%nl+jsite%nl-1)), am1lp1 (0:(isite%nl+jsite%nl-1))
         Real (Kind=DEF_DBL_PREC) :: pow1 (0:isite%nl), pow2 (0:jsite%nl)
         Real (Kind=DEF_DBL_PREC) :: pref (0:isite%nl, 0:jsite%nl)
         Real (Kind=DEF_DBL_PREC) :: dy, dx, dz, sum, rat2, rat1, cit, twolm1, ajm

!!$    mnl = nlmax, mnlsq = mnl ** 2, mlmax = 2 * (nlmax-1)

         w1 = isite%scale
         w2 = jsite%scale
         sc = 0.0d0
         nl1 = isite%nl
         nl2 = jsite%nl

         tx = jsite%rvec (1) - isite%rvec(1)
         ty = jsite%rvec (2) - isite%rvec(2)
         tz = jsite%rvec (3) - isite%rvec(3)

         r = Sqrt (tx**2+ty**2+tz**2)
         rmin = 0.01 * (w1+w2)
         If (r .Lt. rmin) Return
         lmax1 = nl1 - 1
         lmaxh1 = 2 * lmax1

         lmax2 = nl2 - 1
         lmaxh2 = 2 * lmax2

         lmax = Max (lmax1, lmax2)
         lmaxh = Max (lmaxh1, lmaxh2)


         dfac (0) = 1.0d0
         Do l = 1, lmaxh
            twolm1 = Real (2*l-1, kind=DEF_DBL_PREC)
            dfac (l) = twolm1 * dfac (l-1)
         End Do

         am1lp1 (0) = - 1.0d0
         Do l = 1, lmax
            am1lp1 (l) = - am1lp1 (l-1)
         End Do

         rat1 = w1 / r
         pow1 (0) = Sqrt (rat1)
         Do l = 1, lmax1
            pow1 (l) = rat1 * pow1 (l-1)
         End Do

         rat2 = w2 / r
         pow2 (0) = Sqrt (rat2)
         Do l = 1, lmax2
            pow2 (l) = rat2 * pow2 (l-1)
         End Do

         Do l2 = 0, lmax2
            Do l1 = 0, lmax1
               l = l1 + l2
               cit = 8.0d0 * DEF_M_PI * am1lp1 (l2) * dfac (l) * pow1 (l1) * pow2 (l2)
               ajm = dfac (l1) * dfac (l2)
               pref (l1, l2) = cit / ajm
            End Do
         End Do

         dx = tx / r
         dy = ty / r
         dz = tz / r

         yps => harm (lmaxh, dx, dy, dz, cnh)

         i2 = 0
         Do l2 = 0, lmax2
            Do m2 = - l2, l2
               i2 = i2 + 1

               i1 = 0
               Do l1 = 0, lmax1
                  Do m1 = - l1, l1
                     i1 = i1 + 1

                     l = l1 + l2
                     ista = l ** 2 + 1
                     ifin = (l+1) ** 2
                     sum = 0.0d0
                     Do i = ista, ifin
                        sum = sum + gfrh (i, i1, i2) * yps (i)
                        sc (i1, i2) = pref (l1, l2) * sum

                     End Do
                  End Do
               End Do
            End Do
         End Do
         Deallocate (yps)
      End Subroutine cansc

      Subroutine KKRsc (isite, jsite, avw, sc, cg, indxcg, jcg, cy, kap2, only_mstrx3, iop)
!!$---------------------------------------------------------------
!!$   KKR structure constant connecting 2 sites "isite" and "jsite"
!!$       sc(i1,i2) - the KKR structure constant matrix
!!$   OR sc(i1,i2) is the ENERGY DERIVATIVE, if iop=2
!!$---------------------------------------------------------------
!!$   remark: cg, indxcg, jcg and cy are pre-calculated with the  
!!$           symlnc() and scg() subroutines in omta_strrs.o
!!$---------------------------------------------------------------
         Use omta_strrs
         Implicit None
         Type (t_cluster_site) :: isite, jsite
         Real (Kind=DEF_DBL_PREC) :: avw
         Integer :: indxcg(*), jcg(*)
         Real (Kind=DEF_DBL_PREC) :: cy(*), cg(*)
         Complex (Kind=DEF_DBL_PREC) :: sc (:, :)
         Real (Kind=DEF_DBL_PREC) :: kap2
         Integer, Optional :: only_mstrx3
         Integer :: iop
!!$ local
         Real (Kind=DEF_DBL_PREC) :: tx, ty, tz, r
         Integer :: nl1, nl2, nlboth
         Real (Kind=DEF_DBL_PREC) :: dr(3)
         Real (Kind=DEF_DBL_PREC), dimension(isite%nl**2,jsite%nl**2) :: s0a, s0j
         Real (Kind=DEF_DBL_PREC) :: bohr=0.529177210903d0
!!$    mnl = nlmax, mnlsq = mnl ** 2, mlmax = 2 * (nlmax-1)
         

         sc = 0.0d0
         nl1 = isite%nl**2
         nl2 = jsite%nl**2
         nlboth = isite%nl-1 + jsite%nl-1

         tx = (jsite%rvec (1) - isite%rvec(1)) / avw 
         ty = (jsite%rvec (2) - isite%rvec(2)) / avw 
         tz = (jsite%rvec (3) - isite%rvec(3)) / avw 

         r = tx**2+ty**2+tz**2
         
         dr = (/ tx, ty, tz /)

         s0a=0.d0

         if (.not. present(only_mstrx3)) then
         call mstrx2(kap2, dr, nl1, nl2, nlboth, cg, indxcg, jcg, cy, s0a, iop)
 
          if (kap2 <= 0.d0) then
           sc(1:size(s0a,1), 1:size(s0a,2)) = DCMPLX(s0a,0d0)
           return
          else 
           s0j=0.d0
           call mstrx3(kap2, dr, nl1, nl2, nlboth, cg, indxcg, jcg, cy, s0j, iop)
           call renorm_sJ(nl1,nl2,s0j,kap2)
           sc(1:size(s0a,1), 1:size(s0a,2)) = DCMPLX(s0a, s0j)
           return
          endif
         else
          if (kap2 <= 0.d0) then
           sc = 0.d0 
           return
          else
           s0j=0.d0
           call mstrx3(kap2, dr, nl1, nl2, nlboth, cg, indxcg, jcg, cy, s0j, iop)
           call renorm_sJ(nl1,nl2,s0j,kap2)
           sc(1:size(s0a,1), 1:size(s0a,2)) = DCMPLX(s0a, s0j)
           return
          endif
         endif 
      End Subroutine KKRsc


      Function harm_nc (lmax) Result (cnh)
!!$                   calculation of normal. constants for spherical harmonics
         Implicit None
         Integer, Intent (In) :: lmax
         Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :)
!!$ local
         Integer :: m, l, i, mlmax
         Real (Kind=DEF_DBL_PREC) :: twolp1, prod

         If (lmax .Lt. 2) Stop ' *** harm: lmax.lt.2 !'

         mlmax = 2 * (lmax-1)
         Allocate (cnh(0:mlmax, 0:mlmax))
!!$    write(*,*) 'size=', size(cnh,dim=1),size(cnh,dim=2)
         cnh (0, 0) = 0.25d0 * DEF_M_1_PI
         Do l = 1, mlmax
            twolp1 = Real (2*l+1, kind=DEF_DBL_PREC)
            cnh (l, 0) = twolp1 * 0.25d0 * DEF_M_1_PI
            Do m = 1, l
               prod = 1.0d0
               Do i = l - m + 1, l + m
                  prod = prod * Real (i, kind=DEF_DBL_PREC)
               End Do
               cnh (l, m) = twolp1 / (2.0d0*DEF_M_PI*prod)
            End Do
         End Do

         Do l = 0, mlmax
            Do m = 0, l
               cnh (l, m) = Sqrt (cnh(l, m))
            End Do
         End Do
      End Function harm_nc



!!$*******************
!!$xxx    harm    ****
!!$*******************
      Function harm (lmax, d1, d2, d3, cnh) Result (yps)
!!$---------------------------------------------------------
!!$   real(kind=prec) spherical harmonics normalized to unity
!!$     y(l,m;d1,d2,d3) = c(l,|m|)*p(l,|m|;theta)*t(m;phi)
!!$         where: cos(theta)=d3,
!!$                sin(theta)*exp(i*phi)=d1+i*d2,
!!$                c(l,|m|) is positive normalizing factor,
!!$                p(l,|m|;theta) is legendre polynomial,
!!$                t(m;phi) = cos(m*phi)   for m.ge.0,
!!$                t(m;phi) = sin(|m|*phi) for m.lt.0.
!!$---------------------------------------------------------
!!$          input:
!!$    lmax - maximal angular number (2.le.lmax.le.mlmax)
!!$    d1,d2,d3 - direction cosines
!!$    chn - constants, should be precalculated by harm_nc
!!$         output:
!!$    yps(.) - real(kind=prec) spherical harmonics
!!$             the harmonics y(l,m) are stored in array yps(i)
!!$             in the sequence ((y(l,m),m=-l,l),l=0,lmax)
!!$---------------------------------------------------------
!!$  remark:  if (d1,d2,d3) is not a unit vector,
!!$           then y(l,m) are equal to the corresponding
!!$           homogeneous harmonic polynomials:
!!$      y(l,m;d1,d2,d3) = r**l * y(l,m;e1,e2,e3) ,
!!$            where (e1,e2,e3) is a unit vector and
!!$            (d1,d2,d3) = r * (e1,e2,e3) .
!!$----------------------------------------------------------
         Implicit None
!!$ arguments
         Integer, Intent (In) :: lmax
         Real (Kind=DEF_DBL_PREC), Intent (In) :: d1, d2, d3
         Real (Kind=DEF_DBL_PREC), Pointer :: cnh (:, :)
         Real (Kind=DEF_DBL_PREC), Pointer :: yps (:)
!!$ local
         Integer :: m, l, ma, i
         Real (Kind=DEF_DBL_PREC) :: u (0:2*(lmax-1),-2*(lmax-1) :2*(lmax-1))
         Real (Kind=DEF_DBL_PREC) :: dsq, flp1, flm1, twolm1, twolp1
         Complex (Kind=DEF_DBL_PREC) :: zp, zq


         If (lmax .Lt. 2) Stop ' *** harm: lmax.lt.2 !'
!!$      mlmax= 2 * (lmax-1)
!!$      if (lmax .gt. mlmax) stop ' *** harm: lmax.gt.mlmax !'

!!$                     unnormalized harmonics

         dsq = d1 ** 2 + d2 ** 2 + d3 ** 2
         zp = Cmplx (d1, d2, kind=DEF_DBL_PREC)
         zq = DEF_cmplx_one
         u (0, 0) = 1.0d0

         Do l = 1, lmax
            twolm1 = Real (2*l-1, kind=DEF_DBL_PREC)
            zq = twolm1 * zp * zq
            u (l,-l) = imag (zq)
            u (l, l) = Real (zq, kind=DEF_DBL_PREC)
         End Do

         Do m = - lmax + 1, lmax - 1
            l = Abs (m)
            twolp1 = Real (2*l+1, kind=DEF_DBL_PREC)
            u (l+1, m) = twolp1 * d3 * u (l, m)
         End Do

         Do m = - lmax + 2, lmax - 2
            ma = Abs (m)
            Do l = ma + 1, lmax - 1
               twolp1 = Real (2*l+1, kind=DEF_DBL_PREC)
               flm1 = Real (l+ma, kind=DEF_DBL_PREC)
               flp1 = Real (l-ma+1, kind=DEF_DBL_PREC)
               u (l+1, m) = (twolp1*d3*u(l, m)-flm1*dsq*u(l-1, m)) / flp1
            End Do
         End Do

!!$                            normalized harmonics
         Allocate (yps((2*(lmax-1)+1)**2))
         i = 0
         Do l = 0, lmax
            Do m = - l, l
               i = i + 1
               ma = Abs (m)
               yps (i) = cnh (l, ma) * u (l, m)
            End Do
         End Do

      End Function harm


!!$*******************
!!$xxx   gaunty   ****
!!$*******************
      Function gaunty (nl) Result (gfrh)
         Implicit None
         Integer :: nl
         Real (Kind=DEF_DBL_PREC), Pointer :: gfrh (:, :, :)
!!$ local
         Integer :: lmax, lmaxh, i2, l2, i1, l1, i, l, m, m1, m2

!!$------------------------------------------------------
!!$   gaunt's coefficients for real(kind=prec) spherical harmonics:
!!$      input - nl,  output - common/ghar/, where
!!$       gfrh(i,i1,i2) = integral over the unit sphere
!!$                    of  y(l,m) * y(l1,m1) * y(l2,m2).
!!$   here the simple indices i,i1,i2 correspond to the
!!$   composed (l,m), (l1,m1), (l2,m2), respectively.
!!$   the maximal value of l1 and l2 is (nl-1),
!!$   the maximal value of l is two times greater.
!!$------------------------------------------------------

         lmax = nl - 1
         lmaxh = 2 * lmax


         Allocate (gfrh((lmaxh+1)**2, nl**2, nl**2))

         i2 = 0
         Do l2 = 0, lmax
            Do m2 = - l2, l2
               i2 = i2 + 1

               i1 = 0
               Do l1 = 0, lmax
                  Do m1 = - l1, l1
                     i1 = i1 + 1

                     i = 0
                     Do l = 0, lmaxh
                        Do m = - l, l
                           i = i + 1

                           gfrh (i, i1, i2) = gacor (l, m, l1, m1, l2, m2)

                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Function gaunty

!!$*******************
!!$xxx   gacor    ****
!!$*******************
      Function gacor (kl1, km1, kl2, km2, kl3, km3)
!!$----------------------------------------------------
!!$     gaunt coefficient for real(kind=prec) spherical harmonics
!!$----------------------------------------------------
!!$     arguments must satisfy natural inequalities
!!$   kl1.ge.abs(km1), kl2.ge.abs(km2), kl3.ge.abs(km3)
!!$----------------------------------------------------
         Implicit None
         Real (Kind=DEF_DBL_PREC) :: gacor
         Integer :: kl1, km1, kl2, km2, kl3, km3
!!$ local
         Integer :: l, m, l1, m1, is, ir, it, iy, ix, ip, itmax, itmin, m2, l2
         Integer :: ig2, ig, mtest, kma1, kma2, kma3, mamax, mamin
         Real (Kind=DEF_DBL_PREC) :: delta, pref2, a, b, pref1, gamma, sum


         gacor = 0.0d0

         ig2 = kl1 + kl2 + kl3
         If (Mod(ig2, 2) .Ne. 0) Return
         ig = ig2 / 2
         If (Max(kl1, kl2, kl3) .Gt. ig) Return

         mtest = (2*km1+1) * (2*km2+1) * (2*km3+1)
         If (mtest .Lt. 0) Return
         kma1 = Abs (km1)
         kma2 = Abs (km2)
         kma3 = Abs (km3)
         mamax = Max (kma1, kma2, kma3)
         If (2*mamax .Ne. kma1+kma2+kma3) Return
!!$                                        phi-integral
         mamin = Min (kma1, kma2, kma3)
         If (mamin .Eq. 0) Then
            delta = 1.0d0 / Sqrt (2.0d0*DEF_M_PI)
         Else
            delta = 0.5d0 / Sqrt (DEF_M_PI)
            If (km1+km2+km3 .Eq. 0) delta = - delta
         End If
!!$                                      theta-integral
!!$   see e.u. condon, g.h. shortley:
!!$      the theory of atomic spectra (cambridge 1957)
!!$      page 176, eq.(11)

         If (kma1 .Eq. mamax) Then
            l = kl2
            m = kma2
            l1 = kl3
            m1 = kma3
            l2 = kl1
            m2 = kma1
         Else
            If (kma2 .Eq. mamax) Then
               l = kl1
               m = kma1
               l1 = kl3
               m1 = kma3
               l2 = kl2
               m2 = kma2
            Else
               If (kma3 .Eq. mamax) Then
                  l = kl1
                  m = kma1
                  l1 = kl2
                  m1 = kma2
                  l2 = kl3
                  m2 = kma3
               End If
            End If
         End If

         ix = l2 + m2
         iy = l + l1 - m2
         ip = l2 - m2
         ir = l - l1 + m2
         is = l1 - m1
         itmin = Max (0,-ir)
         itmax = Min (iy, ip, is)
         sum = 0.0d0
         Do it = itmin, itmax
            a = (-1.0d0) ** it * tuff (ix+it) * tuff (iy-it)
            b = tuff (ip-it) * tuff (ir+it) * tuff (is-it) * tuff (it)
            sum = sum + a / b
         End Do

         a = Real ((2*l+1)*(2*l1+1)*(2*l2+1)) * tuff (ip) * tuff (l+m) * tuff (l1+m1) * tuff (is)
         b = 2.0d0 * tuff (ix) * tuff (l-m)
         pref2 = Sqrt (a/b)

         a = (-1.0d0) ** (ig-l-m1) * tuff (ig2-2*l1) * tuff (ig)
         b = tuff (ig-l) * tuff (ig-l1) * tuff (ig-l2) * tuff (ig2+1)
         pref1 = a / b

         gamma = pref1 * pref2 * sum

!!$                                  gaunt coefficient
         gacor = gamma * delta
         Return
      End Function gacor

!!$*******************
!!$xxx    tuff    ****
!!$*******************
      Function tuff (n)
!!$----------------------------------------------
!!$        factorial of a non-negative integer
!!$----------------------------------------------
         Implicit None
         Integer :: n
         Real (Kind=DEF_DBL_PREC) :: tuff
!!$ local
         Integer :: i

         tuff = 1.0
         If (n .Le. 1) Return
         Do i = 2, n
            tuff = tuff * Real (i, kind=DEF_DBL_PREC)
         End Do
      End Function tuff

End Module structure

Subroutine getscsym (n, sc, f, crit)
      Implicit None
      Integer, Intent (In) :: n
      Real (Kind=DEF_DBL_PREC), Intent (In) :: crit
      Real (Kind=DEF_DBL_PREC), Intent (Inout) :: sc (n)
      Integer, Intent (Out) :: f (2*n)
!!$  Local
      Integer, Pointer :: mask (:)
      Integer :: i, j, mptr, num
      Real (Kind=DEF_DBL_PREC) :: base
      Allocate (mask(n))
      mask (1:n) = 0
      f (1) = 0
      mptr = 1
      Do i = 1, n, 1
         base = sc (i)
!!$          If (Abs(base) > crit .And. mask(i) == 0) Then
         If (mask(i) == 0) Then
            num = 1
            Do j = i + 1, n, 1
!!$                If ((mask(j) == 0) .And. (Abs(sc(j)) > crit)) Then
               If ((mask(j) == 0)) Then
                  If (Abs(base-sc(j)) < crit) Then
                     mask (j) = 1
                     num = num + 1
                     f (mptr+num) = j
                  Else
                     If (Abs(base+sc(j)) < crit) Then
                        mask (j) = 1
                        num = num + 1
                        f (mptr+num) = - j
                     End If
                  End If
               End If
            End Do
            If (num > 1) Then
               f (mptr) = num
               f (mptr+1) = i
               mptr = mptr + num + 1
               f (mptr) = 0
            End If
         End If
      End Do
      Deallocate (mask)
End Subroutine getscsym

Subroutine setscsym (n, sc, f)
      Implicit None
      Integer, Intent (In) :: n
      Real (Kind=DEF_DBL_PREC), Intent (Inout) :: sc (n)
      Integer, Intent (In) :: f (2*n)
!!$  Local
      Integer :: i, mptr, cntr
      Real (Kind=DEF_DBL_PREC) :: base
      mptr = 1
      Do While (f(mptr) /=  0)
         base = 0.0d0
         cntr = f (mptr)
!!$           write(*,*)
         Do i = mptr + 1, mptr + cntr, 1
            base = base + Real (sign(1, f(i)), kind=DEF_DBL_PREC) * sc (Abs(f(i)))
!!$               write(*,*) sc(abs(f(i)))
         End Do
         base = base / Real (cntr, kind=DEF_DBL_PREC)
         Do i = mptr + 1, mptr + cntr, 1
            sc (Abs(f(i))) = Real (sign(1, f(i)), kind=DEF_DBL_PREC) * base
         End Do
         mptr = mptr + cntr + 1
      End Do
End Subroutine setscsym

Subroutine symscons (n, sc, crit)
      Implicit None
      Integer, Intent (In) :: n
      Real (Kind=DEF_DBL_PREC), Intent (Inout) :: sc (n)
      Real (Kind=DEF_DBL_PREC) :: crit
!!$ Locals
      Integer :: i, j, num
      Integer, Pointer :: mask (:), refs (:)
      Real (Kind=DEF_DBL_PREC) :: base, dev, abase, aref, acc, ref
      Allocate (mask(n))
      Allocate (refs(n))
      mask (1:n) = 0
!!$       write(*,*)crit
      Do i = 1, n, 1
         If (mask(i) == 0) Then
            base = sc (i)
            abase = Abs (base)
            acc = abase
            num = 1
            refs (1) = i * sign (1.0d0, base)
            num = 1
            Do j = i + 1, n, 1
               If ((mask(j) == 0)) Then
                  ref = sc (j)
                  aref = Abs (ref)
                  dev = Abs (2.0d0*(abase-aref)/(abase+aref))
!!$                   write(*,*)'dev=',dev
                  If (dev <= crit) Then
                     num = num + 1
                     refs (num) = j * sign (1.0d0, ref)
!!$                      write(*,*) base, ref
                     acc = acc + aref
                     mask (j) = 1
                  End If
               End If
            End Do

            If (num > 1) Then
               acc = acc / Real (num, kind=DEF_DBL_PREC)
               Do j = 1, num, 1
                  sc (Abs(refs(j))) = Real (sign(1, refs(j)), kind=DEF_DBL_PREC) * acc
               End Do
            End If
         End If
      End Do
      Deallocate (mask, refs)
End Subroutine symscons



