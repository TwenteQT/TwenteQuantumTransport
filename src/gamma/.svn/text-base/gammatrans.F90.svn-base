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

Program gammatrans
      Use transport
      Use ando_module
      Use atoms_module
      Use geometry_module
      Use structure
      Use supercell
      Use hamiltonian
      Use logging
      Use sparse_solvers
      Use rotations
      Use readcfg
      Use postprocess
      Use hdf5io

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
      Type (t_tranrefl) :: tran (2)
      Type (zdensemat) :: rhs
      Real (Kind=DEF_DBL_PREC) :: iitime, ttime, ptime
      Integer :: havelr, haverl
      Real (Kind=DEF_DBL_PREC) :: k (2), cons
      Real (Kind=DEF_DBL_PREC) :: t_lr (2, 2) = 0.0d0, r_lr (2, 2) = 0.0d0, t_lr_spec (2, 2) = 0.0d0, &
     & r_lr_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: t_rl (2, 2) = 0.0d0, r_rl (2, 2) = 0.0d0, t_rl_spec (2, 2) = 0.0d0, &
     & r_rl_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: r_shc (2) = 0.0d0, l_shc (2) = 0.0d0, r_polj = 0.d0, l_polj = 0.d0
      Real (Kind=DEF_DBL_PREC), Allocatable :: torqsigma (:, :, :)
      Real (Kind=DEF_DBL_PREC), Pointer :: rotm (:, :)
      Real (Kind=DEF_DBL_PREC), Allocatable :: zcurr (:, :, :)
      Integer :: rnm = 0, lnm = 0, hops, nmod (4), is1, is2
      Character (Len=clen) :: cwork
      Real (Kind=DEF_DBL_PREC), External :: dznrm2
      Integer (HID_T) :: h5out
      Character (Len=2) :: slbl = 'ud'

!!$ MPI vars
      Integer :: ierr, my_mpi_id, solve_comm, mpi_sz, havejob, root = 0
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

      Call trans_config_init (opt, ispar_in=(mpi_sz > 1))
      opt%po%full = 1
      opt%par_nprock = 0

!!$ Set LogLevel and open output files
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

      If (my_mpi_id == 0) Then
         k = opt%bzo%kpoint
         Call hdfwrite (h5out, 'kpoint', k)
!! TIMES CHANGE HERE
         Write (cwork, '("   Init Time = ", f12.4," seconds")') MPI_Wtime () - iitime
         Call do_log (1, trim(cwork))
         Call do_log (1, 'Starting transport calculation...')
      End If

      Call MPI_Barrier (solve_comm, ierr)

!!$ Workaround of bug in Fortran MPI interface
      Call mrealloc (msys%c, 1, 1)


!SMP$ DO SERIAL
      iitime = MPI_Wtime ()

      havejob = 0
      If (my_mpi_id == 0) havejob = 1

      If (havejob /= 0) Then
         Call make_sk (sk, msc, k, hops-1)
      End If

      cons = 0.0d0
      If (my_mpi_id == root) Then
         havelr = 0
         haverl = 0
!!$ Prepare Ando for left lead
         lnm = get_sc_ando (1, l_ando, k, lpar, lsc, lscell, opt%po,-1, opt%actl)
!!$ Prepare Ando for Right lead
         rnm = get_sc_ando (1, r_ando, k, rpar, rsc, rscell, opt%po, 1, opt%actl)

         nmod (1:2) = l_ando%SNin
         nmod (3:4) = r_ando%SNin

         Rewind (timef)
         Write (timef, '("   Prev Time = ", f12.4," seconds")') ptime
         Write (timef, '("   Ando Time = ", f12.4," seconds, kpt=",i)') MPI_Wtime () - iitime, 1
         flush (timef)
         If (l_ando%Nin*opt%needlr > 0) havelr = 1
         If (r_ando%Nin*opt%needrl > 0) haverl = 1

         If ((havelr+haverl) < 1) Then
            Write (*,*) 'nothing to do!'
            Stop
         End If

!!$ Prepare P-S or H-HOH system for whole device
         Call prep_system (msys, sk, mpar, opt%po, 0)
         Write (timef, '("   Prep Time = ", f12.4," seconds")') MPI_Wtime () - iitime
         flush (timef)
!!$  call dump_matrix(msys%c,'h0')
!!$ Embed boundary conditions and prepare RHS
         Call embed_boundary (msys%c, rhs, l_ando, r_ando, opt%needlr, opt%needrl)
         Write (timef, '("  Embed Time = ", f12.4," seconds")') MPI_Wtime () - iitime
         flush (timef)
         If (opt%writeham == 1) Then
            Call hdfwrite (h5out, '/ham', msys%c)
         End If

         If (opt%writerhs == 1) Then
            Call hdfwrite (h5out, '/rhs', rhs%bl)
         End If
      End If
      If (opt%lowmem /= 0) Then
         Call free_ando_emb (l_ando)
         Call free_ando_emb (r_ando)
      End If
      opt%leq%keepmtx = 1

!!$ Solve system
      Call sp_solve_drhs (msys%c, rhs, root, solve_comm, opt%leq)

!!$ All postprocessing/output, which depends on solution for local K-poin comes here
      If (my_mpi_id == root) Then
         If (opt%writewf == 1) Then
            Call hdfwrite (h5out, '/wavefunc', rhs%bl)
         End If

         If (msys%c%alloc /= 1) Call prep_system (msys, sk, mpar, opt%po, 0)

         Allocate (torqsigma(7, trgeo%tna, 2))
         Allocate (zcurr(trgeo%tna, 2, 2))
         Call calc_torque (rhs, trgeo, rotm, torqsigma, l_ando%Nin*havelr, r_ando%Nin*haverl)
         Call calc_zcurrent (rhs, msys%c, trgeo, zcurr, l_ando%Nin*havelr, r_ando%Nin*haverl)

!!$ Calculate conductance
         cons = cond_full (tran, rhs, l_ando, r_ando, havelr, haverl)

         Write (cwork, '(i7,1x,i4,2(1x,f8.4),5(1x,i4),3x,f8.4,3x,g17.8)') 1, 0, k (1), k (2), floor &
        & (Log10(cons)), lnm, rnm, 1, 1, 1.0d0, cons
         Call do_log (1, trim(cwork))


         t_lr = tran(1)%T%val
         r_lr = tran(1)%R%val
         t_rl = tran(2)%T%val
         r_rl = tran(2)%R%val
         l_shc = l_ando%SNin
         r_shc = r_ando%SNin
         l_polj = l_ando%ppj * (l_ando%SNin(1)+l_ando%SNin(2))
         r_polj = r_ando%ppj * (r_ando%SNin(1)+r_ando%SNin(2))


         t_lr_spec = tran(1)%T%spec
         r_lr_spec = tran(1)%R%spec
         t_rl_spec = tran(2)%T%spec
         r_rl_spec = tran(2)%R%spec

         Write (timef, '(i3.3," x KP Time = ", f12.4," seconds")') 1, MPI_Wtime () - iitime
         ptime = MPI_Wtime () - iitime
         flush (timef)

         Call free (rhs)
         If (opt%writetm /= 0) Then
            Do is1 = 1, 2
               Do is2 = 1, 2
                  If (opt%needlr > 0) Then
                     If (tran(1)%T%exists(is1, is2) /= 0) Call hdfwrite (h5out, '/tmx/lr/'//slbl(is1:&
                    & is1)//slbl(is2:is2), tran(1)%T%mat(is1, is2)%bl)

                     If (tran(1)%R%exists(is1, is2) /= 0) Call hdfwrite (h5out, '/tmx/ll/'//slbl(is1:&
                    & is1)//slbl(is2:is2), tran(1)%R%mat(is1, is2)%bl)
                  End If
                  If (opt%needrl > 0) Then
                     If (tran(2)%T%exists(is1, is2) /= 0) Call hdfwrite (h5out, '/tmx/rl/'//slbl(is1:&
                    & is1)//slbl(is2:is2), tran(2)%T%mat(is1, is2)%bl)
                     If (tran(2)%R%exists(is1, is2) /= 0) Call hdfwrite (h5out, '/tmx/rr/'//slbl(is1:&
                    & is1)//slbl(is2:is2), tran(2)%R%mat(is1, is2)%bl)
                  End If
               End Do
            End Do

         End If
!!$          Call write_trans_parallel (h5tro, opt, solve_comm, root, jobs, bz, ik, tran, cons, nmod)
         Call hdfwrite (h5out, 'sd/Ltorque', torqsigma(1:3, :, 1))
         Call hdfwrite (h5out, 'sd/Rtorque', torqsigma(1:3, :, 2))
         Call hdfwrite (h5out, 'sd/Lsigma', torqsigma(4:7, :, 1))
         Call hdfwrite (h5out, 'sd/Rsigma', torqsigma(4:7, :, 2))
         Call hdfwrite (h5out, 'sd/torque', torqsigma(1:3, :, 1)-torqsigma(1:3, :, 2))
         Call hdfwrite (h5out, 'sd/sigma', torqsigma(4:7, :, 1)-torqsigma(4:7, :, 2))
         Call hdfwrite (h5out, '/current/Lz', transpose(zcurr(:, 1:2, 1)))
         Call hdfwrite (h5out, '/current/Rz', transpose(zcurr(:, 1:2, 2)))
         Call hdfwrite (h5out, '/current/consz', zcurr(:, 1, 1)+zcurr(:, 1, 2)+zcurr(:, 2, 1)+zcurr(:, 2, 2))
         Call hdfwrite (h5out, '/cond/lr', t_lr)
         Call hdfwrite (h5out, '/cond/rl', t_rl)
         Call hdfwrite (h5out, '/cond/ll', r_lr)
         Call hdfwrite (h5out, '/cond/rr', r_rl)
         Call hdfwrite (h5out, '/cond/lr_spec', t_lr_spec)
         Call hdfwrite (h5out, '/cond/rl_spec', t_rl_spec)
         Call hdfwrite (h5out, '/cond/consint', cons)
         Call hdfwrite (h5out, '/cond/conspeak', cons)
         Call hdfwrite (h5out, '/cond/l_shc', l_shc)
         Call hdfwrite (h5out, '/cond/r_shc', r_shc)

         Call write_legacy_ando (opt, t_lr, t_rl, t_lr_spec, t_rl_spec, l_shc, r_shc, cons)
         Call h5fclose_f (h5out, error)
         Call h5close_f (error)
         Write (cwork, '("  Total Time = ", f12.4," seconds")') MPI_Wtime () - ttime
         Call do_log (1, trim(cwork))
      End If

      Call do_log (1, 'Done!')


      Call MPI_Barrier (solve_comm, ierr)

      Call mpi_finalize (ierr)


!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop

Contains


      Subroutine write_legacy_ando (opt, t_lr, t_rl, t_lr_spec, t_rl_spec, l_shc, r_shc, cons)
!!$ Write results
         Implicit None
         Type (t_options) :: opt
         Real (Kind=DEF_DBL_PREC) :: t_lr (2, 2), t_lr_spec (2, 2)
         Real (Kind=DEF_DBL_PREC) :: t_rl (2, 2), t_rl_spec (2, 2)
         Real (Kind=DEF_DBL_PREC) :: r_shc (2), l_shc (2)
         Real (Kind=DEF_DBL_PREC) :: cons
!!$Local
         Integer, Parameter :: fl = 200

!!$ Output results into "andout"
         Open (Unit=fl, File='andout', Action='write')

         Write (fl,*) ' Conservation: integrated,      worst'
         Write (fl, '("         ",2(5x,g13.6))') cons, cons

         Write (fl,*) ' Left Sharvin conductance :'
         Write (fl, '("     by spin :",2(3x,f15.10))') l_shc
         Write (fl, '("     total   :",3x,f15.10)') sum (l_shc)
         Write (fl,*) ' Right Sharvin conductance :'
         Write (fl, '("     by spin :",2(3x,f15.10))') r_shc
         Write (fl, '("     total   :",3x,f15.10)') sum (r_shc)
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
         End If
         Close (fl)
      End Subroutine write_legacy_ando

End Program gammatrans
