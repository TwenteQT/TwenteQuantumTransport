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

Program hdamp
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
      Use mpi
      Implicit None
      !Include 'mpif.h'

!!$ Local
      Integer, Parameter :: clen = 200, timef = 201
      Type (t_options) :: opt
      Type (t_atoms_set) :: atoms
      Type (t_geometry) :: mgeo, lgeo, rgeo, trgeo, llgeo, rlgeo
      Type (t_strconst) :: msc, lsc, rsc
      Type (t_ando_sollution) :: l_ando, r_ando
      Type (t_potpar_mats) :: mpar, lpar, rpar
      Type (t_mathop) :: msys, sk
      Type (t_supercell) :: lscell, rscell
      Type (bzone) :: bz
      Type (t_tranrefl) :: tran (2)
      Type (zdensemat) :: rhs, xdiffs (2)
      Real (Kind=DEF_DBL_PREC) :: iitime, ttime, ptime
      Integer :: ik, havelr, haverl, needsolve, nkloops, loop, hasdata = 0
      Real (Kind=DEF_DBL_PREC) :: k (2), cons, consr (2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: t_lr (2, 2) = 0.0d0, r_lr (2, 2) = 0.0d0, t_lr_spec (2, 2) = 0.0d0, &
     & r_lr_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: t_rl (2, 2) = 0.0d0, r_rl (2, 2) = 0.0d0, t_rl_spec (2, 2) = 0.0d0, &
     & r_rl_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: r_shc (2) = 0.0d0, l_shc (2) = 0.0d0, r_polj = 0.d0, l_polj = 0.d0
      Complex (Kind=DEF_DBL_PREC) :: mix (2) = dcmplx (0.0d0, 0.0d0)
      Complex (Kind=DEF_DBL_PREC) :: damp (2, 2) = dcmplx (0.0d0, 0.0d0), damp1 (2, 2) = dcmplx (0.0d0, &
     & 0.0d0)
      Real (Kind=DEF_DBL_PREC), Pointer :: rotm (:, :)
      Real (Kind=DEF_DBL_PREC), Allocatable :: zcurr (:, :, :), zcurr1 (:, :, :)
      Integer :: rnm = 0, lnm = 0, hops, ngc (3), ngr (3), ntr, nmodl (4), nmodr (4), is1, is2
      Character (Len=clen) :: cwork
      Real (Kind=DEF_DBL_PREC), External :: dznrm2
      Type (t_mask) :: mask
      Type (zcsrmat) :: hdiffs (2)
      Integer (HID_T) :: h5out
      Type (h5_iterout) :: h5tro
      Type (t_h5_openarray) :: h5damp
      Character (Len=2) :: slbl = 'ud'
      Type (torqsigmas) :: tors, tors1
!!$ MPI vars
      Integer, Allocatable :: jobs (:), nsols (:)
      Integer :: ierr, my_mpi_id, local_comm, solve_comm, mpi_sz, boss, havejob, root = 0
      Integer :: mpi_loc_id, error
      Integer :: solsize, ngrps, grpid, i, sl1 (2), sl2 (2), nth
      Type (t_decomp_data) :: ldecd, rdecd
      Real (Kind=DEF_DBL_PREC) :: magvals (3), magmom (3), mdir (2), zr2 (2) = 0.0d0
      Integer, External :: OMP_GET_MAX_THREADS


!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------

      Call init_random_seed ()
      Call mpi_init (ierr)
      solve_comm = mpi_comm_world
      Call mpi_comm_rank (solve_comm, my_mpi_id, ierr)
      Call mpi_comm_size (solve_comm, mpi_sz, ierr)
      nth = - 1
#if defined(_OPENMP) || defined(HASOMPF)
      nth = OMP_GET_MAX_THREADS ()
#endif

      Log_Level = - 1
      If (my_mpi_id == root) Then
         Call unlink ('andout')
         Log_Level = 1
      End If
#ifdef _VERS_
      Call do_log (1, 'Transport code v'//_VERS_//' by Antst (HDF)')
#else
      Call do_log (1, 'Transport code vXXXX by Antst (HDF)')
#endif
!!$ Read  file with parameters

      Write (cwork, '("TASKS=",i4,",  NTHREADS/TASK=",i3)') mpi_sz, nth
      Call do_log (1, trim(cwork))
      Call trans_config_init (opt, ispar_in=(mpi_sz > 1))

      opt%po%full = 1
      opt%needlr = 1
      opt%needrl = 1

      Call redistcomm (opt%par_nprock, solve_comm, local_comm, ngrps, grpid)
      Call mpi_comm_rank (local_comm, mpi_loc_id, ierr)
      Call mpi_comm_size (local_comm, solsize, ierr)
      If (mpi_loc_id == 0) boss = 1
      Allocate (nsols(solsize))

!!$ Set LogLevel and open output files
      Allocate (jobs(mpi_sz))
      If (my_mpi_id == root) Then
         Log_Level = opt%loglvl
         Call hdf5lib_init ()
         Call h5fcreate_f ('transout.h5', H5F_ACC_TRUNC_F, h5out, error)
         Call hdfsetversion (h5out)

         compat_flags (cpf_base) = 1
         compat_flags (cpf_no_spin_separation) = 1
         compat_flags (cpf_reliabledecompose) = 1
         compat_flags (cpf_hasdamping) = 1
         compat_flags (cpf_haslr) = opt%needlr
         compat_flags (cpf_hasrl) = opt%needrl

         Open (Unit=timef, File='times', Action='write')
      Else
!!$ uncomment line bellow if you want output from rest of CPUs
!!$     Log_Level=opt%loglvl
         Write (cwork, '(i2.2)') my_mpi_id
         Write (logfilename, '("log.",A)') trim (cwork)
      End If


      ttime = MPI_Wtime ()
      iitime = MPI_Wtime ()
      ptime = 0
      Call init_mathop (msys)
      Call init_mathop (sk)
      hops = 1
      If (opt%po%kind == 3) hops = 2

      If (mpi_loc_id == 0) Then
         hasdata = 1
         If (opt%po%kind <= 3) Then
!!$ Read atomic potentials
            atoms = read_atoms ('atomlist')
!!$ Read geometries
            lgeo = read_geom ('geom_l', atoms)
            rgeo = read_geom ('geom_r', atoms)
            mgeo = read_geom ('geom_m', atoms, lgeo%scale)
            If (my_mpi_id == root) Then
               Call hdfwrite (h5out, 'geom/mgeo', mgeo)
               Call hdfwrite (h5out, 'geom/rgeo', rgeo)
               Call hdfwrite (h5out, 'geom/lgeo', lgeo)
            End If

!!$ Prepare geometry for leads
            llgeo = make_leadgeom (lgeo, hops,-1)
            rlgeo = make_leadgeom (rgeo, hops, 1)
!!$ Free original geometry for leads, we don't need it anymore
            Call free_geom (lgeo)
            Call free_geom (rgeo)
!!$ Prepare geometry for transport region
            trgeo = make_transportgeom (llgeo, rlgeo, mgeo)
            If (my_mpi_id == root) Then
               Call hdfwrite (h5out, 'geom/trgeo', trgeo)
               Call hdfwrite (h5out, 'geom/llgeo', llgeo)
               Call hdfwrite (h5out, 'geom/rlgeo', rlgeo)
            End If

!!$ Free geometry for 'm' region, we don't need it anymore
            Call free_geom (mgeo)
!!$ Prepare parameters for supercell calculation
            Call prepare_supercell (llgeo, lscell, 1, split=1)
            Call prepare_supercell (rlgeo, rscell, 1, split=1)

            rotm => prep_rot_mask (trgeo, opt%po)
            Call load_mask (mask, trgeo, lscell%size*llgeo%num, rscell%size*rlgeo%num, 1)
            mdir = rotspheric (zr2, opt%po%rot_angle)
            If (calc_magmoms(trgeo, rotm, magvals, magmom)) Then
               compat_flags (cpf_hasmagmom) = 1
               If (opt%damp%axesdir /= 0) Then
                  mdir (1) = Acos (magmom(3)/magvals(2))
                  mdir (2) = Atan2 (magmom(2)/magvals(2), magmom(1)/magvals(2))
               End If
            End If

!!$ Calculate S(R)
            msc = calc_screal (trgeo, 1)
            lsc = calc_screal (llgeo, 1)
            rsc = calc_screal (rlgeo, 1)
!!$ If you want to have loop over Energy  do it here
!!$ Just change opt%Eoffset

!!$ Calculate potential parameters
            Call calc_potpar (atoms, opt%atopt, opt%Eoffset)
!!$             write(*,*) (atoms%at(ierr)%split(:),ierr=1,atoms%num)
!!$ Prepare potential parameters matrices
!!$ for leads and transport region
            lpar = make_ppar_full (llgeo, opt%po, rotm(:, 1:1))
            rpar = make_ppar_full (rlgeo, opt%po, rotm(:, size(rotm, 2) :size(rotm, 2)))
            mpar = make_ppar_full (trgeo, opt%po, rotm)
            If (opt%sdcdir /= 0) Then
               lscell%rot = spinrm (rotm(:, 1:1))
               rscell%rot = spinrm (rotm(:, size(rotm, 2) :size(rotm, 2)))
            End If
         Else If (opt%po%kind == 5) Then
!!$ Do some preparing parts before constructing Hamiltonian
!!$ Brillouin zone information to be fixed
            If (Mod(opt%po%ngridx, opt%po%nlx) /= 0 .Or. Mod(opt%po%ngridy, opt%po%nly) /= 0) Then
               Call do_log (1, 'Model Hamiltonian setting; supercell error')
               Stop
            Else
               llgeo%sc_size (1) = opt%po%ngridx / opt%po%nlx
               llgeo%sc_size (2) = opt%po%ngridy / opt%po%nly
            End If
            llgeo%base (1, 1) = dble (opt%po%nlx) * opt%po%dx
            llgeo%base (1, 2) = 0.d0
            llgeo%base (2, 2) = dble (opt%po%nly) * opt%po%dy
            llgeo%base (2, 1) = 0.d0
            llgeo%norbit = 2 * opt%po%nlx * opt%po%nly * opt%po%nlz
            If (Mod(opt%po%ngridx, opt%po%nrx) /= 0 .Or. Mod(opt%po%ngridy, opt%po%nry) /= 0) Then
               Call do_log (1, 'Model Hamiltonian setting; supercell error')
               Stop
            Else
               rlgeo%sc_size (1) = opt%po%ngridx / opt%po%nrx
               rlgeo%sc_size (2) = opt%po%ngridy / opt%po%nry
            End If
            trgeo%base (1, 1) = dble (llgeo%sc_size(1)) * llgeo%base(1, 1)
            trgeo%base (1, 2) = 0.d0
            trgeo%base (2, 1) = 0.d0
            trgeo%base (2, 2) = dble (llgeo%sc_size(2)) * llgeo%base(2, 2)

            rlgeo%base (1, 1) = dble (opt%po%nrx) * opt%po%dx
            rlgeo%base (1, 2) = 0.d0
            rlgeo%base (2, 2) = dble (opt%po%nry) * opt%po%dy
            rlgeo%base (2, 1) = 0.d0
            rlgeo%norbit = 2 * opt%po%nrx * opt%po%nry * opt%po%nrz

            Call prepare_supercell (llgeo, lscell, 1, split=1, ngrid=opt%po%nlx*opt%po%nly*opt%po%nlz)
            Call prepare_supercell (rlgeo, rscell, 1, split=1, ngrid=opt%po%nrx*opt%po%nry*opt%po%nrz)

         End If !opt%po%kind

         Call alloc_torque (tors, trgeo%tna, trgeo%ltna, trgeo%rtna)
         Call alloc_torque (tors1, trgeo%tna, trgeo%ltna, trgeo%rtna)

         Allocate (zcurr(trgeo%tna, 2, 2))
         Allocate (zcurr1(trgeo%tna, 2, 2))
         zcurr (:, :, :) = 0.0d0

!!$ Prepare grid for BZ integration

         Call gibz (bz, trgeo%base, opt%bzo)
      Else
         Call alloc_torque (tors, 1, 1, 1)
         Call alloc_torque (tors1, 1, 1, 1)
         Allocate (zcurr(1, 2, 2))
         Allocate (zcurr1(1, 2, 2))
      End If


      If (my_mpi_id == 0) Then
         Call hdfwrite (h5out, '/damp/mdir', mdir)
         Call hdfwrite (h5out, '/opts', opt)
         Call hdfwrite (h5out, '/geom/trgeo/rotm', rotm)
         Call hdfwrite (h5out, '/geom/mgeo/rotm', rotm(:, (trgeo%ltna+1) :(trgeo%ltna+trgeo%tna)))

         Call hdfwrite (h5out, '/damp/magvals', magvals)
         Call hdfwrite (h5out, '/damp/magmom', magmom)
         Call hdfwrite (h5out, 'bz', bz)
         Call init_iterout (h5out, h5tro, bz%nkgrid)
         Call hdf_mkopen2 (h5out, '/damp/kdep', h5damp, 4, bz%nkgrid, hdf_complex_t(H5T_NATIVE_DOUBLE))
         nkloops = ceiling (bz%nkgrid/real(ngrps, kind=DEF_DBL_PREC))
!! TIMES CHANGE HERE
         Write (cwork, '("   Init Time = ", f12.4," seconds")') MPI_Wtime () - iitime
         Call do_log (1, trim(cwork))
         Call do_log (1, 'Starting transport calculation...')
      End If
      consr (1) = 0.0d0

!!$       Call do_log (1, '     IBZ, spin,  kx,      ky,   cons,nr_lm,nr_rm,ikx,iky,   weight')


      Call mpi_bcast (nkloops, 1, mpi_integer, 0, solve_comm, ierr)

      If (mpi_loc_id == 0 .And. opt%po%kind == 5) Then
!!$ rotate magnetization for eq_type = 5
         ntr = dble (opt%po%ngridx*opt%po%ngridy*(opt%po%nlz+opt%po%nrz+opt%po%ngridz))
!!$         rotm => prep_rot_mask (ntr, llgeo%norbit/2, rlgeo%norbit/2, opt%po)

         If (opt%sdcdir /= 0) Then
            lscell%rot = spinrm (rotm(:, 1:1))
            rscell%rot = spinrm (rotm(:, size(rotm, 2) :size(rotm, 2)))
         End If
      End If

!!$ Workaround of bug in Fortran MPI interface
      Call mrealloc (msys%c, 1, 1)
      sl2 = (/ 4, 1 /)
      sl1 = (/ 1, 1 /)

!SMP$ DO SERIAL
      Do loop = 0, nkloops - 1
         iitime = MPI_Wtime ()
         ik = 1 + loop * (ngrps) + grpid

         havejob = 0
         If (ik <= bz%nkgrid) Then
            k = bz%k (:, ik)
            havejob = 1
         End If
         If ((mpi_loc_id /= 0) .And. (solsize > 1)) havejob = 0
         If (havejob /= 0) Then
            If (opt%po%kind <= 3) Call make_sk (sk, msc, k, hops-1)
         End If
         Call mpi_gather (havejob, 1, mpi_integer, jobs, 1, mpi_integer, root, solve_comm, ierr)

         cons = 0.0d0
         havelr = 0
         haverl = 0
         needsolve = 0
         If (havejob /= 0) Then

            If (opt%po%kind <= 3) Then
!!$ Prepare Ando for left lead
               lnm = get_sc_ando (1, l_ando, k, lpar, lsc, lscell, opt%po,-1, opt%actl, decd=ldecd)
!!$ Prepare Ando for Right lead
               rnm = get_sc_ando (1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl, decd=rdecd)
!!$                Call hdfwrite (h5out, '/l_ando', l_ando)
!!$                Call hdfwrite (h5out, '/r_ando', r_ando)

            Else If (opt%po%kind == 5) Then
!!$ Prepare Ando for left lead
               lnm = get_sc_ando (1, l_ando, k, lpar, lsc, lscell, opt%po,-1, opt%actl, rotm(:, 1:1), &
              & decd=ldecd)
!!$ Prepare Ando for Right lead
               rnm = get_sc_ando (1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl, rotm(:, size(rotm, 2) &
              & :size(rotm, 2)), decd=rdecd)
            End If

            nmodl (1:2) = l_ando%SNin
            nmodr (1:2) = r_ando%SNin
            nmodl (3:4) = l_ando%SNout
            nmodr (3:4) = r_ando%SNout

            If (my_mpi_id == 0) Then
               Rewind (timef)
               Write (timef, '("   Prev Time = ", f12.4," seconds")') ptime
               Write (timef, '("   Ando Time = ", f12.4," seconds, kpt=",i3)') MPI_Wtime () - iitime, ik
               flush (timef)
            End If
            If (l_ando%Nin*opt%needlr > 0) havelr = 1
            If (r_ando%Nin*opt%needrl > 0) haverl = 1
            needsolve = havelr + haverl
            If (needsolve /= 0) Then
!!$ Prepare P-S or H-HOH system for whole device
               If (opt%po%kind <= 3) Then
                  Call prep_system (msys, sk, mpar, opt%po, 0)

                  Call get_ham_diffs (trgeo, sk, rotm, mask, hdiffs, opt%po, opt%damp)
!!$                   Call hdfwrite (h5out, '/hd/1', hdiffs(1))
!!$                   Call hdfwrite (h5out, '/hd/2', hdiffs(2))
               Else If (opt%po%kind == 5) Then
!!$             calculated the dimension of the Hamiltonian
                  ngc (1) = opt%po%ngridx
                  ngc (2) = opt%po%ngridy
                  ngc (3) = opt%po%ngridz + opt%po%nlz + opt%po%nrz
                  ngr (1) = ngc (1)
                  ngr (2) = ngc (2)
                  ngr (3) = opt%po%ngridz
                  Call prep_system (msys, opt%po, ngc, ngr, k, rotm)
               End If


               If (my_mpi_id == 0) Then
                  Write (timef, '("   Prep Time = ", f12.4," seconds")') MPI_Wtime () - iitime
                  flush (timef)
               End If
!!$  call dump_matrix(msys%c,'h0')
!!$ Embed boundary conditions and prepare RHS
               Call embed_boundary (msys%c, rhs, l_ando, r_ando, opt%needlr, opt%needrl)
               If (my_mpi_id == 0) Then
                  Write (timef, '("  Embed Time = ", f12.4," seconds")') MPI_Wtime () - iitime
                  flush (timef)
               End If
            End If
         End If

         If (opt%writeveloc == 1) Then
            Call hdfwrtset_par (h5out, '/kset', ik, 'leadvel/l_in', solve_comm, root, l_ando%Nin, l_ando%vin)
            Call hdfwrtset_par (h5out, '/kset', ik, 'leadvel/l_out', solve_comm, root, l_ando%Nout, &
           & l_ando%vout)
            Call hdfwrtset_par (h5out, '/kset', ik, 'leadvel/r_in', solve_comm, root, r_ando%Nin, r_ando%vin)
            Call hdfwrtset_par (h5out, '/kset', ik, 'leadvel/r_out', solve_comm, root, r_ando%Nout, &
           & r_ando%vout)
         End If

         If (opt%lowmem /= 0) Then
            Call free_ando_emb (l_ando)
            Call free_ando_emb (r_ando)
         End If

         If (opt%actl%writedecomp == 1) Then
            Call hdfwrtset_par (h5out, '/kset', ik, 'ldecomp', solve_comm, root, needsolve, ldecd)
            Call hdfwrtset_par (h5out, '/kset', ik, 'rdecomp', solve_comm, root, needsolve, rdecd)
         End If

         If (opt%writeham == 1) Then
            Call hdfwrtset_par (h5out, '/kset', ik, 'ham', solve_comm, root, needsolve, msys%c)
         End If

         If (opt%writerhs == 1) Then
            Call hdfwrtset_par (h5out, '/kset', ik, 'rhs', solve_comm, root, needsolve, rhs%bl)
         End If
!!$          Call MPI_Barrier (local_comm, ierr)
         Call mpi_allgather (needsolve, 1, mpi_integer, nsols, 1, mpi_integer, local_comm, ierr)

!!$ Solve system
!SMP$ DO SERIAL
         Do boss = 0, solsize - 1
            If (nsols(boss+1) /= 0) Then
               Call sp_solve_diffs (msys%c, hdiffs, rhs, xdiffs, boss, local_comm, opt%leq)
            End If
         End Do

         If (opt%writewf == 1) Then
            Call hdfwrtset_par (h5out, '/kset', ik, 'wavefunc', solve_comm, root, needsolve, rhs%bl)
            Call hdfwrtset_par (h5out, '/kset', ik, 'wfdiff1', solve_comm, root, needsolve, xdiffs(1)%bl)
            Call hdfwrtset_par (h5out, '/kset', ik, 'wfdiff2', solve_comm, root, needsolve, xdiffs(2)%bl)
         End If

         damp1 = dcmplx (0.0d0, 0.0d0)
!!$ All postprocessing/output, which depends on solution for local K-poin comes here
         If (needsolve > 0) Then
            Call calc_torque (rhs, trgeo, rotm, tors1, l_ando%Nin*havelr, r_ando%Nin*haverl)

            Call bzsum_torq (tors, tors1, bz%kweight(ik))

            Call calc_zcurrent (rhs, msys%c, trgeo, zcurr1, l_ando%Nin*havelr, r_ando%Nin*haverl)
            zcurr (:, :, :) = zcurr (:, :, :) + zcurr1 (:, :, :) * bz%kweight(ik)

            damp1 = calc_damp (xdiffs, l_ando, r_ando, opt%damp%dstep)
            damp = damp + bz%kweight(ik) * damp1
         End If

         sl1 (2) = ik
!!$ Write k-resolved damping
         Call hdf_wrtslab2z_par (solve_comm, root, havejob, h5damp, damp1, sl1, sl2)

         If (opt%writetorque == 1) Then
!!$ WARNING! Memory layout must be changed to avoid creation of temporary arrays here
            Call hdfwrtset_par (h5out, '/kset', ik, 'torsig', solve_comm, root, needsolve, tors1)
         End If

         If (opt%do_oldcurr > 1) Then
            Call hdfwrtset_par (h5out, '/kset', ik, '/current/Lz', solve_comm, root, needsolve, zcurr1(:, &
           & 1:2, 1))
            Call hdfwrtset_par (h5out, '/kset', ik, '/current/Rz', solve_comm, root, needsolve, zcurr1(:, &
           & 1:2, 2))
         End If

!!$ Calculate conductance
         If (havejob /= 0) Then
            cons = cond_full (tran, rhs, l_ando, r_ando, havelr, haverl)
            consr (1) = consr (1) + 0.5d0 * cons * bz%kweight(ik)

            Write (cwork, '(i7,1x,i4,2(1x,f8.4),5(1x,i4),3x,f8.4,3x,g17.8)') ik, 0, k (1), k (2), floor &
           & (Log10(cons)), lnm, rnm, bz%ik(1, ik), bz%ik(2, ik), bz%kweight(ik), cons

            If (consr(2) < cons) consr (2) = cons
!!$             Call cond_spec (is1, is2, tran, lscell, rscell)

            t_lr = t_lr + bz%kweight(ik) * tran(1)%T%val
            r_lr = r_lr + bz%kweight(ik) * tran(1)%R%val
            t_rl = t_rl + bz%kweight(ik) * tran(2)%T%val
            r_rl = r_rl + bz%kweight(ik) * tran(2)%R%val

            l_shc = l_shc + bz%kweight(ik) * l_ando%SNin
            r_shc = r_shc + bz%kweight(ik) * r_ando%SNin

            l_polj = l_polj + bz%kweight(ik) * l_ando%ppj * (l_ando%SNin(1)+l_ando%SNin(2))
            r_polj = r_polj + bz%kweight(ik) * r_ando%ppj * (r_ando%SNin(1)+r_ando%SNin(2))


            t_lr_spec = t_lr_spec + bz%kweight(ik) * tran(1)%T%spec
            r_lr_spec = r_lr_spec + bz%kweight(ik) * tran(1)%R%spec
            t_rl_spec = t_rl_spec + bz%kweight(ik) * tran(2)%T%spec
            r_rl_spec = r_rl_spec + bz%kweight(ik) * tran(2)%R%spec

            Call free (rhs)
         End If
         If (opt%writetm /= 0) Then
            Do is1 = 1, 2
               Do is2 = 1, 2
                  If (opt%needlr > 0) Then
                     Call hdfwrtset_par (h5out, '/kset', ik, 'tmx/lr/'//slbl(is1:is1)//slbl(is2:is2), &
                    & solve_comm, root, tran(1)%T%exists(is1, is2), tran(1)%T%mat(is1, is2)%bl)
                     Call hdfwrtset_par (h5out, '/kset', ik, 'tmx/ll/'//slbl(is1:is1)//slbl(is2:is2), &
                    & solve_comm, root, tran(1)%R%exists(is1, is2), tran(1)%R%mat(is1, is2)%bl)
                  End If
                  If (opt%needrl > 0) Then
                     Call hdfwrtset_par (h5out, '/kset', ik, 'tmx/rl/'//slbl(is1:is1)//slbl(is2:is2), &
                    & solve_comm, root, tran(2)%T%exists(is1, is2), tran(2)%T%mat(is1, is2)%bl)
                     Call hdfwrtset_par (h5out, '/kset', ik, 'tmx/rr/'//slbl(is1:is1)//slbl(is2:is2), &
                    & solve_comm, root, tran(2)%R%exists(is1, is2), tran(2)%R%mat(is1, is2)%bl)
                  End If
               End Do
            End Do

         End If

         Call MPI_Barrier (solve_comm, ierr)

         Call write_trans_parallel (h5tro, opt, solve_comm, root, jobs, bz, ik, tran, cons, nmodl, nmodr)
         If (opt%loglvl > 1) Call log_parallel (solve_comm, root, havejob, jobs, cwork)


         If (my_mpi_id == 0) Then
            Write (timef, '(i3.3," x KP Time = ", f12.4," seconds")') ngrps, MPI_Wtime () - iitime
            ptime = MPI_Wtime () - iitime
            flush (timef)
         End If
      End Do

      Call MPI_Sum_conditional_d (tors%lm, tors1%lm, 7*trgeo%tna, hasdata, 0, solve_comm)
      Call MPI_Sum_conditional_d (tors%rm, tors1%rm, 7*trgeo%tna, hasdata, 0, solve_comm)

      Call MPI_Sum_conditional_d (tors%Lt, tors1%Lt, 7*trgeo%rtna, hasdata, 0, solve_comm)
      Call MPI_Sum_conditional_d (tors%Lir, tors1%Lir, 7*trgeo%ltna, hasdata, 0, solve_comm)

      Call MPI_Sum_conditional_d (tors%Rt, tors1%Rt, 7*trgeo%ltna, hasdata, 0, solve_comm)
      Call MPI_Sum_conditional_d (tors%Rir, tors1%Rir, 7*trgeo%rtna, hasdata, 0, solve_comm)

      Call MPI_Sum_conditional_d (zcurr, zcurr1, 2*trgeo%tna*2, hasdata, 0, solve_comm)

      If (my_mpi_id == 0) Then
         Call hdfwrite (h5out, 'torsig', tors1)
         Call hdfwrite (h5out, 'torsig/L-R', tors1%lm(:, :)-tors1%rm(:, :))
         Call hdfwrite (h5out, 'current/Lz', transpose(zcurr1(:, 1:2, 1)))
         Call hdfwrite (h5out, 'current/Rz', transpose(zcurr1(:, 1:2, 2)))
         Call hdfwrite (h5out, 'current/consz', zcurr1(:, 1, 1)+zcurr1(:, 1, 2)+zcurr1(:, 2, 1)+zcurr1(:, 2, &
        & 2))
      End If

      Call do_log (1, 'Done!')


      Call MPI_summ_res_d (solve_comm, 0, t_lr, 4, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, t_rl, 4, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, r_lr, 4, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, r_rl, 4, mpi_sum)

      Call MPI_summ_res_d (solve_comm, 0, t_lr_spec, 4, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, t_rl_spec, 4, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, r_lr_spec, 4, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, r_rl_spec, 4, mpi_sum)

      Call MPI_summ_res_d (solve_comm, 0, l_shc, 2, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, r_shc, 2, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, consr(1), 1, mpi_sum)
      Call MPI_summ_res_d (solve_comm, 0, consr(2), 1, mpi_max)

      Call MPI_summ_res_z (solve_comm, 0, damp, 4, mpi_sum)


      If (my_mpi_id == 0) Then
         Call write_legacy_ando (opt, t_lr, t_rl, t_lr_spec, t_rl_spec, l_shc, r_shc, consr, damp)
         Call hdf_write_conds (h5out, t_lr, t_rl, r_lr, r_rl, t_lr_spec, t_rl_spec, r_lr_spec, r_rl_spec, &
        & mix, mix, l_shc, r_shc, consr)

         Call hdfwrite (h5out, '/damp/total', damp)
      End If

      Call MPI_Barrier (solve_comm, ierr)
!!$ End here you loop over E

      If (my_mpi_id == 0) Then
         Call close_iterout (h5tro)
         Call h5dclose_f (h5damp%dset, error)

         Call hdfwrite (h5out, '/compat', compat_flags)
         Call hdfwrite (h5out, '/time', MPI_Wtime()-ttime)
         Call h5fclose_f (h5out, error)
         Call h5close_f (error)

         Write (cwork, '("  Total Time = ", f12.4," seconds")') MPI_Wtime () - ttime
         Write (timef, '("  Total Time = ", f12.4," seconds")') MPI_Wtime () - ttime
         Close (timef)

         Call do_log (1, trim(cwork))
      End If

      Call mpi_finalize (ierr)


!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop

Contains
      Subroutine write_legacy_ando (opt, t_lr, t_rl, t_lr_spec, t_rl_spec, l_shc, r_shc, consr, damp)
!!$ Write results
         Implicit None
         Type (t_options) :: opt
         Real (Kind=DEF_DBL_PREC) :: t_lr (2, 2), t_lr_spec (2, 2)
         Real (Kind=DEF_DBL_PREC) :: t_rl (2, 2), t_rl_spec (2, 2)
         Real (Kind=DEF_DBL_PREC) :: r_shc (2), l_shc (2)
         Real (Kind=DEF_DBL_PREC) :: consr (2)
         Complex (Kind=DEF_DBL_PREC) :: damp (2, 2)

!!$Local
         Integer, Parameter :: fl = 200

!!$ Output results into "andout"
         Open (Unit=fl, File='andout', Action='write')

         Write (fl,*) ' Conservation: integrated,      worst'
         Write (fl, '("         ",2(5x,g13.6))') consr

         Write (fl,*) ' Left Sharvin conductance :'
         Write (fl, '("     by spin :",2(3x,f15.10))') l_shc
         Write (fl, '("     total   :",3x,f15.10)') sum (l_shc)
!!$         If (abs(l_polj1)>1.d-6) &
!!$         Write (fl, '(" Current pol.:",3x,f15.10)') l_polj1/sum(l_shc1)
         Write (fl,*) ' Right Sharvin conductance :'
         Write (fl, '("     by spin :",2(3x,f15.10))') r_shc
         Write (fl, '("     total   :",3x,f15.10)') sum (r_shc)
!!$         If (abs(r_polj1)>1.d-6) &
!!$         Write (fl, '(" Current pol.:",3x,f15.10)') r_polj1/sum(r_shc1)
         Write (fl,*) '---------------------------------------------'
         Write (fl,*) ''
         Write (fl,*) ' Damping constants, Re part:'
         Write (fl, '("   X    :",2(2x,g17.10))') real (damp(1, :))
         Write (fl, '("   Y    :",2(2x,g17.10))') real (damp(2, :))
         Write (fl,*) ''
         Write (fl,*) ' Damping constants, Im part:'
         Write (fl, '("   X    :",2(2x,g17.10))') imag (damp(1, :))
         Write (fl, '("   Y    :",2(2x,g17.10))') imag (damp(2, :))
         Write (fl,*) '---------------------------------------------'
         Write (fl,*) ''

         If (opt%needlr /= 0) Then
            Write (fl,*) ' L->R Conductance(carefull with decomposition!):'
            Write (fl,*) '              Up:                Down:'
            Write (fl, '("   Up   :",2(3x,g17.10))') t_lr (1, :)
            Write (fl, '("   Down :",2(3x,g17.10))') t_lr (2, :)
            Write (fl,*) ''
            Write (fl,*) ' Conductance total:'
            Write (fl, '("         ",2(3x,g17.10))') sum (t_lr(1, :)), sum (t_lr(2, :))
            Write (fl, '("  Grand total:                  ",g17.10)') sum (t_lr(:, :))
            Write (fl,*) ''
            Write (fl,*) ' Specular:'
            Write (fl,*) '              Up:                Down:'
            Write (fl, '("   Up   :",2(3x,g17.10))') t_lr_spec (1, :)
            Write (fl, '("   Down :",2(3x,g17.10))') t_lr_spec (2, :)
            Write (fl,*) ''
         End If

         If (opt%needrl /= 0) Then
            Write (fl,*) ' R->L Conductance(carefull with decomposition!):'
            Write (fl,*) '              Up                Down'
            Write (fl, '("   Up   :",2(3x,g17.10))') t_rl (1, :)
            Write (fl, '("   Down :",2(3x,g17.10))') t_rl (2, :)
            Write (fl,*) ''
            Write (fl,*) ' Conductance total:'
            Write (fl, '("         ",2(3x,g17.10))') sum (t_rl(1, :)), sum (t_rl(2, :))
            Write (fl, '("  Grand total:                  ",g17.10)') sum (t_rl(:, :))
            Write (fl,*) ''
            Write (fl,*) ' Specular:'
            Write (fl,*) '              Up:                Down:'
            Write (fl, '("   Up   :",2(3x,g17.10))') t_rl_spec (1, :)
            Write (fl, '("   Down :",2(3x,g17.10))') t_rl_spec (2, :)
            Write (fl,*) ''
         End If
         Close (fl)
      End Subroutine write_legacy_ando

      Function calc_damp (xdiffs, l_ando, r_ando, dstep) Result (damp)
         Use transport
         Use ando_module
         Implicit None
         Type (zdensemat) :: xdiffs (2)
         Type (t_ando_sollution) :: l_ando, r_ando
         Complex (Kind=DEF_DBL_PREC) :: damp (2, 2)
         Real (Kind=DEF_DBL_PREC) :: dstep
!!$ Local
         Type (zdensemat) :: sder (2)
         Integer :: i1, i2

         Do i1 = 1, 2
            Call cond_Sdiff (sder(i1), xdiffs(i1), l_ando, r_ando)
         End Do


         Do i1 = 1, 2
            Do i2 = 1, 2
               damp (i1, i2) = sum (sder(i1)%bl*conjg(sder(i2)%bl))
            End Do
         End Do
         damp = damp / (2*Sin(dstep/2.d0)) ** 2

         Do i1 = 1, 2
            Call free (sder(i1))
         End Do
      End Function calc_damp


      Subroutine cond_Sdiff (sd, psid, l_ando, r_ando)
!!$    Calculates derivative of scattering matrix
!!$    psid - derivatives of wave-functions
!!$    sd - matrix with the derivative

         Use ando_module
         Implicit None
!!$    Arguments
         Type (t_ando_sollution) :: l_ando, r_ando
         Type (zdensemat) :: psid
         Type (zdensemat) :: sd
!!$    Local variables
         Complex (Kind=DEF_DBL_PREC), Parameter :: z1 = DEF_cmplx_one, z0 = DEF_cmplx_zero
         Integer :: n, l, m, ldr, n1, i1

         ldr = psid%nrow
         l = psid%ncol
         n1 = l_ando%Nout + r_ando%Nout
         Call alloc (sd, n1, l)

         i1 = 0
         If (l_ando%Nout > 0) Then
            n = l_ando%Nout
            m = l_ando%n
            Call zgemm ('N', 'N', n, l, m, z1, l_ando%Uout_i, n, psid%bl(1, 1), ldr, z0, sd%bl(1, 1), n1)
            i1 = i1 + l_ando%Nout
         End If

         If (r_ando%Nout > 0) Then
            n = r_ando%Nout
            m = r_ando%n
            Call zgemm ('N', 'N', n, l, m, z1, r_ando%Uout_i, n, psid%bl(ldr-m+1, 1), ldr, z0, sd%bl(i1+1, &
           & 1), n1)
         End If
      End Subroutine cond_Sdiff

      Subroutine get_ham_diffs (trgeo, sk, rotm0, mask, hdiffs, po, dopt)
         Use sparselib
         Use hamiltonian
         Use structure
         Implicit None
         Type (zcsrmat), Intent (Inout) :: hdiffs (:)
         Real (Kind=DEF_DBL_PREC), Intent (In) :: rotm0 (:, :)
         Type (t_pot_opts), Intent (In) :: po
         Type (t_damp_opts), Intent (In) :: dopt
         Type (t_mathop), Intent (In) :: sk
         Type (t_geometry), Intent (In) :: trgeo
         Type (t_mask), Intent (In) :: mask
!!$ Local vars
         Real (Kind=DEF_DBL_PREC), Allocatable :: rotm (:, :)
         Integer :: i, ddir, dind
         Real (Kind=DEF_DBL_PREC) :: aa (2, 1)
         Type (t_potpar_mats) :: mpar
         Type (t_mathop) :: ms (2)

         Allocate (rotm(2, trgeo%num))

         Do ddir = 1, 2
            Do dind = 1, 2
               Do i = 1, trgeo%num
                  aa (:, 1) = getderangle (rotm0(:, i), mdir, dopt%dstep, ddir, dind)
                  rotm (:, i) = mask%a(i) * aa (:, 1) + (1.0d0-mask%a(i)) * rotm0 (:, i)
               End Do
               mpar = make_ppar_full (trgeo, po, rotm)
               Call prep_system (ms(dind), sk, mpar, po, 0)
               Call free_ppar_f (mpar)
            End Do
            Call free (hdiffs(ddir))
            hdiffs (ddir) = spmatadd (ms(2)%c, ms(1)%c, sign=-1)
            Call free_mathop (ms(1))
            Call free_mathop (ms(2))
         End Do
         Deallocate (rotm)
      End Subroutine get_ham_diffs

End Program hdamp
