!!$   $Id:trans.F90 425 2007-05-23 18:18:58Z antst $
!!$    This program calculates the condactance of
!!$    "ideal lead|scattering region|ideal lead" system using Landauer
!!$    formula. The transmission coefficients are calculated according to
!!$    the method described in T.Ando PRB 44, 8017 (1991) modified for
!!$    KKR (P-S) equation in screened representation (K.Xia et al, unpublished).
!!$    The current incarnation is strongly influenced by the previous
!!$    implementation by K.Xia and "enriched" by the outright theft of the
!!$    parts of the self-consistent LMTO code by I.Turek.
#include "math_def.h"

Program htrans_eigvals
   Use transport
   Use ando_module
   Use atoms_module
   Use geometry_module
   Use structure
   Use supercell
   Use hamiltonian
   Use logging
   Use sparse_solvers
   Use bzgrid
   Use rotations
   Use readcfg
   Use postprocess
   Use hdf5io
   Use compat
   Use helpers
   Use interatcur
   Use timer
!!$ EMTO modules
   Use omta_kink
   Use omta_defs
   Use omta_strrs
   Use omta_pots
   Use omta_SOC
   Use omta_nonsph
   Use df
   Use slepc_solvers
!!$ PETSC
   !Use petsc_solvers
   Use mpi
   Implicit None
!   Include 'mpif.h'
!!$ Local
   Integer, Parameter :: clen = 200
   Type(t_options) :: opt
   Type(t_atoms_set) :: atoms
   Type(t_geometry) :: mgeo, lgeo, rgeo, trgeo, llgeo, rlgeo
   Type(t_strconst) :: msc, lsc, rsc
   Type(t_strconst) :: msdotc
   Type(t_ando_sollution) :: l_ando, r_ando
   Type(t_potpar_mats) :: mpar, lpar, rpar
   Type(t_mathop) :: msys, sk
   Type(t_supercell) :: lscell, rscell
   Type(bzone) :: bz
   Type(t_tranrefl) :: tran(2)
   Type(zdensemat) :: rhs, wavefn_phase, temp_sol
   Real(Kind=DEF_DBL_PREC) :: iitime, ttime !, ptime
   Integer :: ik, havelr, haverl, needsolve, nkloops, loop, hasdata = 0
   Real(Kind=DEF_DBL_PREC) :: k(2), cons, consr(2) = 0.0d0
   Complex(Kind=DEF_DBL_PREC), ALlocatable :: eiv(:)
!   Type(zcsrmat) :: testing_matrix
   Real(Kind=DEF_DBL_PREC) :: t_lr(2, 2) = 0.0d0, r_lr(2, 2) = 0.0d0, t_lr_spec(2, 2) = 0.0d0, &
  & r_lr_spec(2, 2) = 0.0d0
   Real(Kind=DEF_DBL_PREC) :: t_rl(2, 2) = 0.0d0, r_rl(2, 2) = 0.0d0, t_rl_spec(2, 2) = 0.0d0, &
  & r_rl_spec(2, 2) = 0.0d0
   Real(Kind=DEF_DBL_PREC) :: r_shc(2) = 0.0d0, l_shc(2) = 0.0d0, r_polj = 0.d0, l_polj = 0.d0
   Complex(Kind=DEF_DBL_PREC) :: mix(2) = dcmplx(0.0d0, 0.0d0)
   Real(Kind=DEF_DBL_PREC), Pointer :: rotm(:, :)
   Real(Kind=DEF_DBL_PREC), Allocatable :: zcurr(:, :, :), zcurr1(:, :, :)
   Integer :: rnm = 0, lnm = 0, hops, ngc(3), ngr(3), ntr, nmodl(4), nmodr(4), is1, is2
   Character(Len=clen) :: cwork
   Real(Kind=DEF_DBL_PREC), External :: dznrm2
   Type(t_mask) :: mask
   Type(zdiagmat) :: tempmat1
!!$       Type (zcsrmat) :: tempmat2
   Integer(HID_T) :: h5out
   Type(h5_iterout) :: h5tro
   Character(Len=2) :: slbl = 'ud'
   Type(torqsigmas) :: tors, tors1
   Type(t_neighbour_list), Pointer :: neibs(:)
   Type(t_remember_neibs), Allocatable :: rneibs(:)
   Type(t_interat_currents) :: iatcu, iatcu1
   Integer :: itr, itr2
!!$ MPI vars
   Integer, Allocatable :: jobs(:), nsols(:)
   Integer :: ierr, my_mpi_id, local_comm, solve_comm, mpi_sz, boss, havejob, root = 0
   Integer :: mpi_loc_id, error
   Integer :: solsize, ngrps, grpid, nth
   Type(t_decomp_data) :: ldecd, rdecd
   Real(Kind=DEF_DBL_PREC) :: magvals(3), magmom(3)
   Integer, External :: OMP_GET_MAX_THREADS
!!$ EMTO vars
   Type(t_atoms_set_EMTO) :: atoms_EMTO
   Type(t_geometry_EMTO) :: mgeo_EMTO, lgeo_EMTO, rgeo_EMTO, trgeo_EMTO, llgeo_EMTO, rlgeo_EMTO
   Type(t_geometry_EMTO) :: trgeo_EMTO_nodf
   Type(t_omta_struct) :: msc_EMTO(2), lsc_EMTO(2), rsc_EMTO(2)
   Type(t_omta_logder) :: mpar_EMTO(2), lpar_EMTO(2), rpar_EMTO(2)
   Type(t_mathop_EMTO) :: tempsys(2)
   Real(kind=DEF_DBL_PREC) :: kappa2
   Integer :: ldimidim(2), rdimidim(2), mdimidim(2), is
   Type(zcsrmat) :: zcsr_c, zcsr_temp, hsoc, op, oph, zcsr_df, zcsr_df2
   Integer :: i, l_ldim, m_ldim, r_ldim, nrow_i, iii
!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------

   Call mpi_init(ierr)
   Call init_random_seed()
   solve_comm = mpi_comm_world
   Call mpi_comm_rank(solve_comm, my_mpi_id, ierr)
   Call mpi_comm_size(solve_comm, mpi_sz, ierr)
   nth = -1
#if defined(_OPENMP) || defined(HASOMPF)
   nth = OMP_GET_MAX_THREADS()
#endif

   iitime = MPI_Wtime()
   ttime = iitime

   Log_Level = -1
   If (my_mpi_id == root) Then
      Call unlink('andout')
      Log_Level = 1
   End If
#ifdef _VERS_
   Call do_log(1, 'Transport code v'//_VERS_//' by Antst (HDF)')
#else
   Call do_log(1, 'Transport code vXXXX by Antst (HDF)')
#endif

   open (file='eigvals_withleads', unit=818, status='replace')
   close (818)
   open (file='eigvals_noleads', unit=818, status='replace')
   close (818)

!!$ Read  file with parameters
   Write (cwork, '("TASKS=",i4,",  NTHREADS/TASK=",i3)') mpi_sz, nth
   Call do_log(1, trim(cwork))

   Call trans_config_init(opt, ispar_in=(mpi_sz > 1))

   opt%po%full = 1

   Call redistcomm(opt%par_nprock, solve_comm, local_comm, ngrps, grpid)
   Call mpi_comm_rank(local_comm, mpi_loc_id, ierr)
   Call mpi_comm_size(local_comm, solsize, ierr)
   If (mpi_loc_id == 0) boss = 1
   Allocate (nsols(solsize))

!!$ Set LogLevel and open output files
   Allocate (jobs(mpi_sz))
   If (my_mpi_id == root) Then
      Log_Level = opt%loglvl
      Call hdf5lib_init()
      Call h5fcreate_f('transout.h5', H5F_ACC_TRUNC_F, h5out, error)
      Call hdfsetversion(h5out)

      compat_flags(cpf_base) = 1
      compat_flags(cpf_no_spin_separation) = 1
      compat_flags(cpf_reliabledecompose) = 1
      compat_flags(cpf_haslr) = opt%needlr
      compat_flags(cpf_hasrl) = opt%needrl

      !Open (Unit=timef, File='times', Action='write')
   Else
!!$ uncomment line bellow if you want output from rest of CPUs
!!$     Log_Level=opt%loglvl
      Write (cwork, '(i2.2)') my_mpi_id
      Write (logfilename, '("log.",A)') trim(cwork)
   End If
   timef = open_timef(my_mpi_id, opt%alltimes)

   Call init_mathop(msys)
   Call init_mathop(sk)
   hops = 1
   If (opt%po%kind == 3) hops = 2

   If (mpi_loc_id == 0) Then
      hasdata = 1
      If (opt%po%kind <= 3) Then
!!$ Read atomic potentials
         atoms = read_atoms('atomlist')
!!$ Read geometries
         lgeo = read_geom('geom_l', atoms)
         rgeo = read_geom('geom_r', atoms)
         mgeo = read_geom('geom_m', atoms, lgeo%scale)
         If (my_mpi_id == root) Then
            Call hdfwrite(h5out, 'geom/mgeo', mgeo)
            Call hdfwrite(h5out, 'geom/rgeo', rgeo)
            Call hdfwrite(h5out, 'geom/lgeo', lgeo)
         End If

!!$ Prepare geometry for leads
         llgeo = make_leadgeom(lgeo, hops, -1)
         rlgeo = make_leadgeom(rgeo, hops, 1)
!!$ Free original geometry for leads, we don't need it anymore
         Call free_geom(lgeo)
         Call free_geom(rgeo)
!!$ Prepare geometry for transport region
         trgeo = make_transportgeom(llgeo, rlgeo, mgeo)
         If (my_mpi_id == root) Then
            Call hdfwrite(h5out, 'geom/trgeo', trgeo)
            Call hdfwrite(h5out, 'geom/llgeo', llgeo)
            Call hdfwrite(h5out, 'geom/rlgeo', rlgeo)
         End If

!!$ Free geometry for 'm' region, we don't need it anymore
         Call free_geom(mgeo)
!!$ Prepare parameters for supercell calculation
         Call prepare_supercell(llgeo, lscell, 1, split=1)
         Call prepare_supercell(rlgeo, rscell, 1, split=1)

         rotm => prep_rot_mask(trgeo, opt%po)
         If (calc_magmoms(trgeo, rotm, magvals, magmom)) compat_flags(cpf_hasmagmom) = 1
!!$ Calculate S(R)
         lsc = calc_screal(llgeo, 1)
         rsc = calc_screal(rlgeo, 1)

         if (opt%do_iatcu == 0 .and. opt%do_iatLcu == 0) msc = calc_screal(trgeo, 1)

         If (opt%do_iatcu /= 0 .or. opt%do_iatLcu /= 0) Then
            allocate (rneibs(trgeo%num))
            msc = calc_screal(trgeo, 1, rneibs=rneibs)
            Call init_iatcu(rneibs, iatcu, trgeo, opt%po%kind, 4*opt%do_iatcu + 3*opt%do_iatLcu)
            iatcu1 = copy_iatcu(iatcu)
            itr = iatcu%ntr
            Call MPI_Bcast(itr, 1, mpi_integer, 0, local_comm, ierr)
         End If

!!$ If you want to have loop over Energy  do it here
!!$ Just change opt%Eoffset

!!$ Calculate potential parameters
!         opt%atopt%ei = opt%po%ei ! the following routine doesn't have access to the opt%ei field,
         ! so we manually add it
         Call calc_potpar(atoms, opt%atopt, opt%Eoffset)

!!$ Prepare potential parameters matrices
!!$ for leads and transport region

         lpar = make_ppar_full(llgeo, opt%po, rotm(:, 1:1))
         rpar = make_ppar_full(rlgeo, opt%po, rotm(:, size(rotm, 2):size(rotm, 2)))
         mpar = make_ppar_full(trgeo, opt%po, rotm)
         If (opt%sdcdir /= 0) Then
            lscell%rot = spinrm(rotm(:, 1:1))
            rscell%rot = spinrm(rotm(:, size(rotm, 2):size(rotm, 2)))
         End If
      Else If (opt%po%kind == 5) Then
!!$ Do some preparing parts before constructing Hamiltonian
!!$ Brillouin zone information to be fixed
         If (Mod(opt%po%ngridx, opt%po%nlx) /= 0 .Or. Mod(opt%po%ngridy, opt%po%nly) /= 0) Then
            Call do_log(1, 'Model Hamiltonian setting; supercell error')
            Stop
         Else
            llgeo%sc_size(1) = opt%po%ngridx/opt%po%nlx
            llgeo%sc_size(2) = opt%po%ngridy/opt%po%nly
         End If
         llgeo%base(1, 1) = dble(opt%po%nlx)*opt%po%dx
         llgeo%base(1, 2) = 0.d0
         llgeo%base(2, 2) = dble(opt%po%nly)*opt%po%dy
         llgeo%base(2, 1) = 0.d0
         llgeo%norbit = 2*opt%po%nlx*opt%po%nly*opt%po%nlz
         If (Mod(opt%po%ngridx, opt%po%nrx) /= 0 .Or. Mod(opt%po%ngridy, opt%po%nry) /= 0) Then
            Call do_log(1, 'Model Hamiltonian setting; supercell error')
            Stop
         Else
            rlgeo%sc_size(1) = opt%po%ngridx/opt%po%nrx
            rlgeo%sc_size(2) = opt%po%ngridy/opt%po%nry
         End If
         trgeo%base(1, 1) = dble(llgeo%sc_size(1))*llgeo%base(1, 1)
         trgeo%base(1, 2) = 0.d0
         trgeo%base(2, 1) = 0.d0
         trgeo%base(2, 2) = dble(llgeo%sc_size(2))*llgeo%base(2, 2)

         rlgeo%base(1, 1) = dble(opt%po%nrx)*opt%po%dx
         rlgeo%base(1, 2) = 0.d0
         rlgeo%base(2, 2) = dble(opt%po%nry)*opt%po%dy
         rlgeo%base(2, 1) = 0.d0
         rlgeo%norbit = 2*opt%po%nrx*opt%po%nry*opt%po%nrz

         Call prepare_supercell(llgeo, lscell, 1, split=1, ngrid=opt%po%nlx*opt%po%nly*opt%po%nlz)
         Call prepare_supercell(rlgeo, rscell, 1, split=1, ngrid=opt%po%nrx*opt%po%nry*opt%po%nrz)

        !!$ ================ Added option: EMTO ==================

      Else if (opt%po%kind == 6) then

         call read_omta_pot('atomlist_EMTO', atoms_EMTO)

         mgeo_EMTO = read_geom_EMTO('geom_m', atoms_EMTO, .false.)
         lgeo_EMTO = read_geom_EMTO('geom_l', atoms_EMTO, .true.)
         rgeo_EMTO = read_geom_EMTO('geom_r', atoms_EMTO, .true.)

         If (my_mpi_id == root) Then
            Call hdfwrite(h5out, 'geom/mgeo', mgeo_EMTO)
            Call hdfwrite(h5out, 'geom/rgeo', rgeo_EMTO)
            Call hdfwrite(h5out, 'geom/lgeo', lgeo_EMTO)
         End If

         hops = 1

         llgeo_EMTO = make_leadgeom_EMTO(lgeo_EMTO, hops, -1)
         rlgeo_EMTO = make_leadgeom_EMTO(rgeo_EMTO, hops, 1)
         call free_geom_EMTO(lgeo_EMTO)
         call free_geom_EMTO(rgeo_EMTO)

         !if (opt%po%df > 0) then
         !   trgeo_EMTO = make_transportgeom_EMTO_df(llgeo_EMTO, rlgeo_EMTO, mgeo_EMTO, atoms_EMTO)
         !   trgeo_EMTO_nodf = make_transportgeom_EMTO(llgeo_EMTO, rlgeo_EMTO, mgeo_EMTO)
         !else
         trgeo_EMTO = make_transportgeom_EMTO(llgeo_EMTO, rlgeo_EMTO, mgeo_EMTO)
         !end if

         If (my_mpi_id == root) Then
            Call hdfwrite(h5out, 'geom/trgeo', trgeo_EMTO)
            Call hdfwrite(h5out, 'geom/llgeo', llgeo_EMTO)
            Call hdfwrite(h5out, 'geom/rlgeo', rlgeo_EMTO)
         End If

         call makidx_lead(llgeo_EMTO, l_ldim)
         call makidx_lead(rlgeo_EMTO, r_ldim)
         call makidx_m(trgeo_EMTO, m_ldim)

         call free_geom_EMTO(mgeo_EMTO)

         call prepare_supercell_EMTO(llgeo_EMTO, lscell, 1, atoms_EMTO%nspin)
         call prepare_supercell_EMTO(rlgeo_EMTO, rscell, 1, atoms_EMTO%nspin)

         rotm => prep_rot_mask_EMTO(trgeo_EMTO, opt%po)
         If (calc_magmoms_EMTO(trgeo_EMTO, rotm, magvals, magmom)) compat_flags(cpf_hasmagmom) = 1

         kappa2 = (atoms_EMTO%fermi - atoms_EMTO%vmtz(1))*trgeo_EMTO%dawsr*trgeo_EMTO%dawsr*trgeo_EMTO%alat*trgeo_EMTO%alat

         !call mkalph(kappa2, trgeo_EMTO%nlmax, trgeo_EMTO%dawsr, atoms_EMTO)
         !call mktral(kappa2, trgeo_EMTO%nlmax, trgeo_EMTO%dawsr, atoms_EMTO)
         call mkalph(kappa2, trgeo_EMTO%nlmax, trgeo_EMTO%dawsr*trgeo_EMTO%alat, atoms_EMTO)
         call mktral(kappa2, trgeo_EMTO%nlmax, trgeo_EMTO%dawsr*trgeo_EMTO%alat, atoms_EMTO)

         call atompp_nmto(atoms_EMTO)
         if (opt%po%nonsph > 0) call calc_omta_nonsph(atoms_EMTO)

         do is = 1, atoms_EMTO%nspin
            mdimidim(is) = m_ldim
            ldimidim(is) = l_ldim
            rdimidim(is) = r_ldim
         end do

         call makpph_nmto(atoms_EMTO%nspin, llgeo_EMTO, ldimidim, lpar_EMTO, .true., opt%po%so, opt%po%nonsph)
         call makpph_nmto(atoms_EMTO%nspin, rlgeo_EMTO, rdimidim, rpar_EMTO, .true., opt%po%so, opt%po%nonsph)
         call makpph_nmto(atoms_EMTO%nspin, trgeo_EMTO, mdimidim, mpar_EMTO, .false., opt%po%so, opt%po%nonsph)

         lsc = calc_screal_EMTO(llgeo_EMTO, kappa2, 0d0, 1, needcur=0, df=opt%po%df, ppar=lpar_EMTO(1))
         rsc = calc_screal_EMTO(rlgeo_EMTO, kappa2, 0d0, 1, needcur=0, df=opt%po%df, ppar=rpar_EMTO(1))

         if (opt%do_iatcu == 0 .and. opt%do_iatLcu == 0) msc = &
           & calc_screal_EMTO(trgeo_EMTO, kappa2, 0d0, 1, needcur=0, df=opt%po%df, ppar=mpar_EMTO(1))

         If (opt%do_iatcu /= 0 .or. opt%do_iatLcu /= 0) Then
            allocate (rneibs(trgeo_EMTO%num))
            msc = calc_screal_EMTO(trgeo_EMTO, kappa2, 0d0, 1, needcur=0, rneibs=rneibs, df=opt%po%df, ppar=mpar_EMTO(1)) !, sdotca=msdotc)
            call init_iatcu_EMTO(rneibs(:), iatcu, trgeo_EMTO, opt%po%kind, 4*opt%do_iatcu + 3*opt%do_iatLcu)
            iatcu1 = copy_iatcu(iatcu)
            itr = iatcu%ntr
            Call MPI_Bcast(itr, 1, mpi_integer, 0, local_comm, ierr)
         End If

         if (opt%po%df > 0) then
            call downfold_logder(lgeo_EMTO, lpar_EMTO(1)%logder)
            call downfold_logder(lgeo_EMTO, lpar_EMTO(1)%hsoc)

            call downfold_logder(rgeo_EMTO, rpar_EMTO(1)%logder)
            call downfold_logder(rgeo_EMTO, rpar_EMTO(1)%hsoc)

            call downfold_logder(trgeo_EMTO, mpar_EMTO(1)%logder)
            call downfold_logder(trgeo_EMTO, mpar_EMTO(1)%hsoc)
         end if

         if (opt%po%nonsph > 0) then
            call mk_nonsphblock(llgeo_EMTO, lpar_EMTO(1)%logder, lsc, k)
            call mk_nonsphblock(rlgeo_EMTO, rpar_EMTO(1)%logder, rsc, k)
            call mk_nonsphblock(trgeo_EMTO, mpar_EMTO(1)%logder, msc, k)
         end if
        !!$ ======================================================

      End If !opt%po%kind
      if (opt%po%kind == 6) then
         call alloc_torque(tors, trgeo_EMTO%tna, trgeo_EMTO%ltna, trgeo_EMTO%rtna)
         call alloc_torque(tors1, trgeo_EMTO%tna, trgeo_EMTO%ltna, trgeo_EMTO%rtna)
      else
         Call alloc_torque(tors, trgeo%tna, trgeo%ltna, trgeo%rtna)
         Call alloc_torque(tors1, trgeo%tna, trgeo%ltna, trgeo%rtna)
      end if

      if (opt%po%kind == 6) then
         Allocate (zcurr(trgeo_EMTO%tna, 2, 2))
         Allocate (zcurr1(trgeo_EMTO%tna, 2, 2))
         zcurr(:, :, :) = 0.0d0
      else

         Allocate (zcurr(trgeo%tna, 2, 2))
         Allocate (zcurr1(trgeo%tna, 2, 2))
         zcurr(:, :, :) = 0.0d0

      end if
!!$       Complex E fix
      If ((opt%po%kind > 1) .And. (opt%po%ei /= 0.0d0)) Then
         Call load_mask(mask, trgeo, lscell%size*llgeo%num, rscell%size*rlgeo%num, 1)
         Call alloc(tempmat1, trgeo%norbit*2)
         tempmat1%el(:) = Cmplx(0.0d0, opt%po%ei, kind=DEF_DBL_PREC)*mask%o(:)
         Call inspmatadd(mpar%c0, tempmat1, 1)
         Call free(tempmat1)
      End If
         !!$ Prepare grid for BZ integration
      if (opt%po%kind == 6) then
         Call gibz(bz, trgeo_EMTO%base, opt%bzo)
      else
         Call gibz(bz, trgeo%base, opt%bzo)
      end if
      !!$write (*, *) 'gibz done! ', mpi_loc_id
   Else !mpi_loc_id .ne. 0
      Call alloc_torque(tors, 1, 1, 1)
      Call alloc_torque(tors1, 1, 1, 1)
      Allocate (zcurr(1, 2, 2))
      Allocate (zcurr1(1, 2, 2))
      If (opt%do_iatcu /= 0 .or. opt%do_iatLcu /= 0) Then
         Call MPI_Bcast(itr, 1, mpi_integer, 0, local_comm, ierr)
         Call make_empty_iatcu(iatcu, itr)
         Call make_empty_iatcu(iatcu1, itr)
      End If
   End If

   If (my_mpi_id == 0) Then
      Call hdfwrite(h5out, '/opts', opt)
      Call hdfwrite(h5out, '/geom/trgeo/rotm', rotm)
      If (opt%po%kind == 6) Then
         Call hdfwrite(h5out, '/geom/mgeo/rotm', rotm(:, (trgeo_EMTO%ltna + 1):(trgeo_EMTO%ltna + trgeo_EMTO%tna)))
      Else
         Call hdfwrite(h5out, '/geom/mgeo/rotm', rotm(:, (trgeo%ltna + 1):(trgeo%ltna + trgeo%tna)))
      End If
      Call hdfwrite(h5out, '/damp/magvals', magvals)
      Call hdfwrite(h5out, '/damp/magmom', magmom)

      Call hdfwrite(h5out, 'bz', bz)
      Call init_iterout(h5out, h5tro, bz%nkgrid)
      nkloops = ceiling(bz%nkgrid/real(ngrps, kind=DEF_DBL_PREC))
!! TIMES CHANGE HERE
      Write (cwork, '("   Init Time = ", f12.4," seconds")') MPI_Wtime() - iitime
      Call do_log(1, trim(cwork))
      Call do_log(1, 'Starting transport calculation...')
   End If
   consr(1) = 0.0d0

!!$       Call do_log (1, '     IBZ, spin,  kx,      ky,   cons,nr_lm,nr_rm,ikx,iky,   weight')
   Call mpi_bcast(nkloops, 1, mpi_integer, 0, solve_comm, ierr)

   If (mpi_loc_id == 0 .And. opt%po%kind == 5) Then
!!$ rotate magnetization for eq_type = 5
      ntr = dble(opt%po%ngridx*opt%po%ngridy*(opt%po%nlz + opt%po%nrz + opt%po%ngridz))
!!$         rotm => prep_rot_mask (ntr, llgeo%norbit/2, rlgeo%norbit/2, opt%po)

      If (opt%sdcdir /= 0) Then
         lscell%rot = spinrm(rotm(:, 1:1))
         rscell%rot = spinrm(rotm(:, size(rotm, 2):size(rotm, 2)))
      End If
   End If

!!$ Workaround of bug in Fortran MPI interface
   Call mrealloc(msys%c, 1, 1)

!SMP$ DO SERIAL
   Do loop = 0, nkloops - 1 !!$ loop over k-points
      Call restart_times(times, MPI_Wtime())

      ik = 1 + loop*(ngrps) + grpid

      havejob = 0
      If (ik <= bz%nkgrid .and. mpi_loc_id == 0) Then
         k = bz%k(:, ik)
         havejob = 1
      End If
      If ((mpi_loc_id /= 0) .And. (solsize > 1)) havejob = 0
      If (havejob /= 0) Then
         If (opt%po%kind <= 3) Call make_sk(sk, msc, k, hops - 1)
         If (opt%po%kind == 6) Call make_sk(sk, msc, k, hops - 1)
      End If
      Call mpi_gather(havejob, 1, mpi_integer, jobs, 1, mpi_integer, root, solve_comm, ierr)
      Call sub_time(times, MPI_Wtime(), ti_mpi)

      cons = 0.0d0
      havelr = 0
      haverl = 0
      needsolve = 0
      If (havejob /= 0) Then
         If (opt%po%kind <= 3) Then
!!$ Prepare Ando for left lead
            lnm = get_sc_ando(1, l_ando, k, lpar, lsc, lscell, opt%po, -1, opt%actl, decd=ldecd)
!!$ Prepare Ando for Right lead
            rnm = get_sc_ando(1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl, decd=rdecd)
         Else If (opt%po%kind == 5) Then
!!$ Prepare Ando for left lead
            lnm = get_sc_ando(1, l_ando, k, lpar, lsc, lscell, opt%po, -1, opt%actl, rotm(:, 1:1), &
           & decd=ldecd)
!!$ Prepare Ando for Right lead
            rnm = get_sc_ando(1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl, rotm(:, size(rotm, 2) &
           & :size(rotm, 2)), decd=rdecd)

!!$ ================ Added option: EMTO ==================

         Else if (opt%po%kind == 6) then

!                lnm = get_sc_ando_EMTO (opt%po, atoms_EMTO%nspin, l_ando, k, lpar_EMTO, llgeo_EMTO, lsc_EMTO, lscell, -1, rotm(:, 1:1), opt%actl%gensol, opt%actl%specpart)
!                rnm = get_sc_ando_EMTO (opt%po, atoms_EMTO%nspin, r_ando, k, rpar_EMTO, rlgeo_EMTO, rsc_EMTO, rscell,  1, rotm(:, size(rotm, 2) : size(rotm, 2)), opt%actl%gensol, opt%actl%specpart)

            lnm = get_sc_ando_EMTO(opt%po, atoms_EMTO%nspin, l_ando, k, lpar_EMTO, llgeo_EMTO,&
                   & lsc, lscell, -1, gens=opt%actl%gensol, frac=opt%actl%specpart)

            rnm = get_sc_ando_EMTO(opt%po, atoms_EMTO%nspin, r_ando, k, rpar_EMTO, rlgeo_EMTO,&
                   & rsc, rscell, 1, gens=opt%actl%gensol, frac=opt%actl%specpart)

!!$ ======================================================

         End If

         nmodl(1:2) = l_ando%SNin
         nmodr(1:2) = r_ando%SNin
         nmodl(3:4) = l_ando%SNout
         nmodr(3:4) = r_ando%SNout

         Call sub_time(times, MPI_Wtime(), ti_ando)

         If (l_ando%Nin*opt%needlr > 0) havelr = 1
         If (r_ando%Nin*opt%needrl > 0) haverl = 1
         needsolve = havelr + haverl
         If (needsolve /= 0) Then
!!$ Prepare P-S or H-HOH system for whole device
            If (opt%po%kind <= 3) Then
               Call prep_system(msys, sk, mpar, opt%po, 0)
            Else If (opt%po%kind == 5) Then
!!$             calculated the dimension of the Hamiltonian
               ngc(1) = opt%po%ngridx
               ngc(2) = opt%po%ngridy
               ngc(3) = opt%po%ngridz + opt%po%nlz + opt%po%nrz
               ngr(1) = ngc(1)
               ngr(2) = ngc(2)
               ngr(3) = opt%po%ngridz
               Call prep_system(msys, opt%po, ngc, ngr, k, rotm)

!!$ ================ Added option: EMTO ==================

            Else if (opt%po%kind == 6) then

               Call prep_system(msys, sk, mpar_EMTO(1), opt%po, 0)

               !!$ Call GetRotOp_EMTO (op, trgeo_EMTO, rotm)
               !!$ oph = spherm (op)
               !!$ Call RotateZCSR (oph, msys%c, op, opt%po%rnegl)
               !!$ Call free (op)
               !!$ Call free (oph)

               !if (opt%po%df > 0) then
               !   call downfold_kink_scattering(msys, trgeo_EMTO)
               !end if

!!$ ======================================================
            End If
            Call sub_time(times, MPI_Wtime(), ti_prep)

            Allocate (eiv(msys%c%nrow))
            Call slepc_solve(msys%c, root, solve_comm, msys%c%nrow, eiv)
            open (file='eigvals_noleads', unit=818, position='append')
            do iii = 1, size(eiv)
               write (818, *) k, real(eiv(iii)), aimag(eiv(iii))
            end do
            close (unit=818)
            Deallocate (eiv)
!!$ Embed boundary conditions and prepare RHS
            If (opt%embed_only_left == 0) then
              Call embed_boundary(msys%c, rhs, l_ando, r_ando, opt%needlr, opt%needrl)
            Else
              Call embed_boundary_only_left(msys%c, rhs, l_ando, opt%needlr)
            Endif

            Call sub_time(times, MPI_Wtime(), ti_embed)
         End If !needsolve
      End If !havejob

      If (opt%writeveloc == 1) Then
         Call hdfwrtset_par(h5out, '/kset', ik, 'leadvel/l_in', solve_comm, root, l_ando%Nin, l_ando%vin)
         Call hdfwrtset_par(h5out, '/kset', ik, 'leadvel/l_out', solve_comm, root, l_ando%Nout, &
        & l_ando%vout)
         Call hdfwrtset_par(h5out, '/kset', ik, 'leadvel/r_in', solve_comm, root, r_ando%Nin, r_ando%vin)
         Call hdfwrtset_par(h5out, '/kset', ik, 'leadvel/r_out', solve_comm, root, r_ando%Nout, &
        & r_ando%vout)
      End If

      If (opt%lowmem /= 0) Then
         Call free_ando_emb(l_ando)
         Call free_ando_emb(r_ando)
      End If

      If (opt%actl%writedecomp == 1) Then
         Call hdfwrtset_par(h5out, '/kset', ik, 'ldecomp', solve_comm, root, needsolve, ldecd)
         Call hdfwrtset_par(h5out, '/kset', ik, 'rdecomp', solve_comm, root, needsolve, rdecd)
      End If

      If (opt%writeham == 1) Then
         Call hdfwrtset_par(h5out, '/kset', ik, 'ham', solve_comm, root, needsolve, msys%c)
      End If

      If (opt%writerhs == 1) Then
         Call hdfwrtset_par(h5out, '/kset', ik, 'rhs', solve_comm, root, needsolve, rhs%bl)
      End If
      Call sub_time(times, MPI_Wtime(), ti_hdf)

      !!$Call MPI_Barrier (local_comm, ierr)
      Call mpi_allgather(needsolve, 1, mpi_integer, nsols, 1, mpi_integer, local_comm, ierr)
      Call sub_time(times, MPI_Wtime(), ti_mpi)

!!$ Solve system
!SMP$ DO SERIAL

      Allocate (eiv(msys%c%nrow))
      Call slepc_solve(msys%c, root, solve_comm, msys%c%nrow, eiv)
      open (file='eigvals_withleads', unit=818, position='append')
      do iii = 1, size(eiv)
         write (818, *) k, real(eiv(iii)), aimag(eiv(iii))
      end do
      close (unit=818)
      Deallocate (eiv)

      !If (opt%leq%solver<2) Then
      !  Do boss = 0, solsize - 1
      !     If (nsols(boss + 1) /= 0) Then
      !        Call sp_solve(msys%c, rhs, boss, local_comm, opt%leq)
      !     End If
      !  End Do

      !Else ! i.e. opt%leq%solver .ge. 2, this is PETSc
      !  !Do boss = 0, solsize - 1
      !  !   If (nsols(boss + 1) /= 0) Then
      !  !    Call petsc_solve(msys%c, rhs, boss, local_comm, opt%petsc_opt)
      !  !   Endif
      !  !Enddo
      !Endif

      Call sub_time(times, MPI_Wtime(), ti_solve)

   End Do !!$ end of loop over k points

!      call alloc(testing_matrix, 3, 3)
!      testing_matrix%a = (/ DEF_cmplx_Ione, -1*DEF_cmplx_Ione, 2*DEF_cmplx_Ione /)
!      testing_matrix%ir = (/ 1, 2, 3, 4 /)
!      testing_matrix%jc = (/ 1, 2, 3 /)
!      testing_matrix%nnz = 3

!      Allocate(eiv(3))
!      Call slepc_solve(testing_matrix, root, solve_comm, 3, eiv)
!      do iii=1,size(eiv)
!      write(*,*) iii, real(eiv(iii)), aimag(eiv(iii))
!      enddo
!      Deallocate(eiv)

   ttime = MPI_Wtime() - ttime
   If (timef /= -1) Then
      Write (timef, '("   Total Time = ", f12.4," seconds")') ttime
      Close (timef)
   End If

   If (my_mpi_id == 0) Then
      Call close_iterout(h5tro)

      Call hdfwrite(h5out, '/compat', compat_flags)
      Call hdfwrite(h5out, '/time', ttime)
      Call h5fclose_f(h5out, error)
      Call h5close_f(error)
      Write (cwork, '("  Total Time = ", f12.4," seconds")') ttime
      !Write (timef, '("   Total Time = ", f12.4," seconds")') ttime
      !Close (timef)

      Call do_log(1, trim(cwork))
   End If

   Call mpi_finalize(ierr)

!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
   Stop
End Program htrans_eigvals
