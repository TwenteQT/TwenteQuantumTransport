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

Program makemat
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
      Use helpers
      Implicit None
      Include 'mpif.h'

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
      Type (zcsrmat) :: rhs
      Type (zdensemat) :: wf
      Real (Kind=DEF_DBL_PREC) :: iitime, ttime, ptime
      Integer :: ik, havelr, haverl, needsolve
      Real (Kind=DEF_DBL_PREC) :: k (2), cons, consr (2) = 0.0d0
      Real (Kind=DEF_DBL_PREC), Pointer :: rotm (:, :)
      Integer :: rnm = 0, lnm = 0, hops
      Character (Len=clen) :: cwork
      Real (Kind=DEF_DBL_PREC), External :: dznrm2
      Type (t_mask) :: mask
      Type (zdiagmat) :: tempmat1
!!$       Type (zcsrmat) :: tempmat2
      Integer (HID_T) :: h5out

!!$ MPI vars
      Integer, Allocatable :: jobs (:), nsols (:)
      Integer :: ierr, my_mpi_id, local_comm, solve_comm, mpi_sz, havejob, root = 0, solsize
      Integer :: mpi_loc_id, error
!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------

      Call init_random_seed ()
      Call mpi_init (ierr)
      solve_comm = mpi_comm_world
      Call mpi_comm_rank (solve_comm, my_mpi_id, ierr)
      Call mpi_comm_size (solve_comm, mpi_sz, ierr)

      Log_Level = - 1
      If (my_mpi_id == root) Then
         Call unlink ('andout')
         Log_Level = 1
      End If
#ifdef _VERS_
      Call do_log (1, 'Transport code v'//_VERS_//' by Antst')
#else
      Call do_log (1, 'Transport code vXXXX by Antst')
#endif
!!$ Read  file with parameters

      Call trans_config_init (opt)
      opt%po%full = 1

      local_comm = solve_comm
      solsize = mpi_sz
      Allocate (nsols(mpi_sz))

      Call mpi_comm_rank (local_comm, mpi_loc_id, ierr)

!!$ Set LogLevel and open output files
      Allocate (jobs(mpi_sz))
      If (my_mpi_id == root) Then
         Log_Level = opt%loglvl
         Call hdf5lib_init ()
         Call h5fcreate_f ('transout.h5', H5F_ACC_TRUNC_F, h5out, error)
         Call hdfsetversion (h5out)
         Call hdfwrite (h5out, '/opts', opt)
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
!!$ Read atomic potentials
         atoms = read_atoms ('atomlist')
!!$ Read geometries
         lgeo = read_geom ('geom_l', atoms)
         rgeo = read_geom ('geom_r', atoms)
         mgeo = read_geom ('geom_m', atoms,lgeo%scale)
         
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

!!$ Calculate S(R)
         msc = calc_screal (trgeo, 1)
         lsc = calc_screal (llgeo, 1)
         rsc = calc_screal (rlgeo, 1)
!!$ If you want to have loop over Energy  do it here
!!$ Just change opt%Eoffset

!!$ Calculate potential parameters
         Call calc_potpar (atoms, opt%atopt, opt%Eoffset)
!!$ Prepare potential parameters matrices
!!$ for leads and transport region
         lpar = make_ppar_full (llgeo, opt%po, rotm(:, 1:1))
         rpar = make_ppar_full (rlgeo, opt%po, rotm(:, size(rotm, 2) :size(rotm, 2)))
         mpar = make_ppar_full (trgeo, opt%po, rotm)
         If (opt%sdcdir /= 0) Then
            lscell%rot = spinrm (rotm(:, 1:1))
            rscell%rot = spinrm (rotm(:, size(rotm, 2) :size(rotm, 2)))
         End If
      End If

!!$       Complex E fix
      If ((opt%po%kind > 1) .And. (opt%po%ei /= 0.0d0)) Then
         Call load_mask (mask, trgeo, lscell%size*llgeo%num, rscell%size*rlgeo%num, 1)
         Call alloc (tempmat1, trgeo%norbit*2)
         tempmat1%el (:) = Cmplx (0.0d0, opt%po%ei, kind=DEF_DBL_PREC) * mask%o(:)
         Call inspmatadd (mpar%c0, tempmat1, 1)
         Call free (tempmat1)
      End If


!!$ Prepare grid for BZ integration

      Call gibz (bz, trgeo%base, opt%bzo)


      consr (1) = 0.0d0

      If (my_mpi_id == root) Then
         Call hdfwrite (h5out, 'bz', bz)
         Write (timef, '("   Init Time = ", f12.4," seconds")') MPI_Wtime () - iitime
         flush (timef)
      End If
      Call MPI_Barrier (solve_comm, ierr)

!!$ Workaround of bug in Fortran MPI interface
      Call mrealloc (msys%c, 1, 1)



      iitime = MPI_Wtime ()
      ik = 1 + my_mpi_id
      havejob = 0
      If (ik <= bz%nkgrid) Then
         k = bz%k (:, ik)
         havejob = 1
      End If
      If ((mpi_loc_id /= 0) .And. (solsize > 1)) havejob = 0
      If (havejob /= 0) Then
         Call make_sk (sk, msc, k, hops-1)
      End If
      Call mpi_gather (havejob, 1, mpi_integer, jobs, 1, mpi_integer, root, solve_comm, ierr)

      cons = 0.0d0
      havelr = 0
      haverl = 0
      needsolve = 0
      If (havejob /= 0) Then
         If (opt%po%kind <= 3) Then
!!$ Prepare Ando for left lead
            lnm = get_sc_ando (1, l_ando, k, lpar, lsc, lscell, opt%po,-1, opt%actl)
!!$ Prepare Ando for Right lead
            rnm = get_sc_ando (1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl)
         Else If (opt%po%kind == 5) Then
!!$ Prepare Ando for left lead
            lnm = get_sc_ando (1, l_ando, k, lpar, lsc, lscell, opt%po,-1, opt%actl, rotm(:, 1:1))
!!$ Prepare Ando for Right lead
            rnm = get_sc_ando (1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl, rotm(:, size(rotm, 2) &
           & :size(rotm, 2)))
         End If

         If (my_mpi_id == 0) Then
            Write (timef, '("   Ando Time = ", f12.4," seconds, kpt=",i)') MPI_Wtime () - iitime, ik
            flush (timef)
         End If
         If (l_ando%Nin*opt%needlr > 0) havelr = 1
         If (r_ando%Nin*opt%needrl > 0) haverl = 1
         needsolve = havelr + haverl
         If (needsolve /= 0) Then
!!$ Prepare P-S or H-HOH system for whole device
            Call prep_system (msys, sk, mpar, opt%po, 0)
            Call hdfwrtset_par (h5out, '/kset', ik, 'ham0', solve_comm, root, needsolve, msys%c)
            Call hdfwrtset_par (h5out, '/kset', ik, 'lemb', solve_comm, root, needsolve, l_ando%emb)
            Call hdfwrtset_par (h5out, '/kset', ik, 'remb', solve_comm, root, needsolve, r_ando%emb)

            If (my_mpi_id == 0) Then
               Write (timef, '("   Prep Time = ", f12.4," seconds")') MPI_Wtime () - iitime
               flush (timef)
            End If
!!$  call dump_matrix(msys%c,'h0')
!!$ Embed boundary conditions and prepare RHS
            Call embed_boundary_srhs (msys%c, rhs, l_ando, r_ando, opt%needlr, opt%needrl)
            If (my_mpi_id == 0) Then
               Write (timef, '("  Embed Time = ", f12.4," seconds")') MPI_Wtime () - iitime
               flush (timef)
            End If
         End If
      End If

      iitime = MPI_Wtime ()
      Call hdfwrtset_par (h5out, '/kset', ik, 'ham', solve_comm, root, needsolve, msys%c)
      Call hdfwrtset_par (h5out, '/kset', ik, 'rhs', solve_comm, root, needsolve, rhs)

      If (my_mpi_id == 0) Then
         Write (timef, '("  Write Time = ", f12.4," seconds")') MPI_Wtime () - iitime
         flush (timef)
         Call h5fclose_f (h5out, error)
         Call h5close_f (error)
      End If
      Call MPI_Barrier (local_comm, ierr)
      Call mpi_finalize (ierr)

!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop

End Program makemat
