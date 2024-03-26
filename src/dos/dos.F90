!!#   $Id:trans.F90 425 2007-05-23 18:18:58Z antst $
!!$    This program calculates the condactance of
!!$    "ideal lead|scattering region|ideal lead" system using Landauer
!!$    formula. The transmission coefficients are calculated according to
!!$    the method described in T.Ando PRB 44, 8017 (1991) modified for
!!$    KKR (P-S) equation in screened representation (K.Xia et al, unpublished).
!!$    The current incarnation is strongly influenced by the previous
!!$    implementation by K.Xia and "enriched" by the outright theft of the
!!$    parts of the self-consistent LMTO code by I.Turek.
#include "math_def.h"

Program trans
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
      Use readcfg
      Implicit None
      Include 'mpif.h'

!!$ Local
      Integer, Parameter :: clen = 200, timef = 201
      Integer, Parameter :: trfiles (2) = (/ 202, 203 /)
      Type (t_options) :: opt
      Type (t_atoms_set) :: atoms
      Type (t_potparp_options) :: optspp
      Type (t_geometry) :: mgeo, lgeo, rgeo, trgeo, llgeo, rlgeo
      Type (t_strconst) :: msc, lsc, rsc
      Type (t_ando_sollution) :: l_ando, r_ando
      Type (t_potpar_mats) :: mpar (2), lpar (2), rpar (2)
      Type (t_mathop) :: msys, sk
      Type (t_supercell) :: lscell, rscell
      Type (bzone) :: bz
      Type (t_tranrefl) :: tran (2)
      Type (zdensemat) :: rhs
      Real (Kind=DEF_DBL_PREC) :: iitime, ttime
      Integer :: is, ik, havelr, haverl, needsolve, nkloops, loop, is1 = 1, is2 = 2
      Real (Kind=DEF_DBL_PREC) :: k (2), consv, cons1, cons, consw = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: t_lr (2, 2) = 0.0d0, r_lr (2, 2) = 0.0d0, t_lr_spec (2, 2) = 0.0d0, &
     & r_lr_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: t_rl (2, 2) = 0.0d0, r_rl (2, 2) = 0.0d0, t_rl_spec (2, 2) = 0.0d0, &
     & r_rl_spec (2, 2) = 0.0d0
      Real (Kind=DEF_DBL_PREC) :: r_shc (2) = 0.0d0, l_shc (2) = 0.0d0
      Complex (Kind=DEF_DBL_PREC) :: tmix (2) = 0.0d0, gmix (2) = 0.0d0
      Integer :: rnm (2) = 0, lnm (2) = 0, hops
      Character (Len=clen) :: cwork
!!$ MPI vars
      Integer, Pointer :: jobs (:), nsols (:)
      Integer :: ierr, my_mpi_id, local_comm, solve_comm, mpi_sz, boss, havejob, root = 0, solsize
      Integer :: mpi_loc_id
!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------

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

      If (opt%par_kind == 1 .And. mpi_sz > 1) Then
         Allocate (nsols(mpi_sz))
         local_comm = solve_comm
         solsize = mpi_sz
      Else
         Call MPI_COMM_SPLIT (solve_comm, my_mpi_id, 0, local_comm, ierr)
         solsize = 1
         boss = 1
         Allocate (nsols(1))
      End If

      Call mpi_comm_rank (local_comm, mpi_loc_id, ierr)

!!$ Set LogLevel and open output files
      Allocate (jobs(mpi_sz))
      If (my_mpi_id == root) Then
         Log_Level = opt%loglvl
         Call open_trans_fles (trfiles, opt)
         Open (Unit=timef, File='times', Action='write')
!$$ Create directories if need
         If (opt%writewf /= 0) Call c_mkdir ('wf', 2)
         If (opt%writeham /= 0) Call c_mkdir ('ham', 3)
         If (opt%writerhs /= 0) Call c_mkdir ('rhs', 3)
         If (opt%writetm /= 0) Call c_mkdir ('tmx', 3)
         If (opt%extlead == 1) Call c_mkdir ('leads', 5)
      Else
!!$ uncomment line bellow if you want output from rest of CPUs
!!$     Log_Level=opt%loglvl
         Write (cwork, '(i2.2)') my_mpi_id
         Write (logfilename, '("log.",A)') trim (cwork)
      End If


      ttime = MPI_Wtime ()
      iitime = MPI_Wtime ()

      Call init_mathop (msys)
      Call init_mathop (sk)
      hops = 1
      If (opt%po%kind > 2) hops = 2


      If (mpi_loc_id == 0) Then
!!$ Read atomic potentials
         atoms = read_atoms ('atomlist')
!!$ Read geometries
         mgeo = read_geom ('geom_m', atoms)
         lgeo = read_geom ('geom_l', atoms)
         rgeo = read_geom ('geom_r', atoms)
!!$ Prepare geometry for leads
         llgeo = make_leadgeom (lgeo, hops,-1)
         rlgeo = make_leadgeom (rgeo, hops, 1)
!!$ Free original geometry for leads, we don't need it anymore
         Call free_geom (lgeo)
         Call free_geom (rgeo)
!!$ Prepare geometry for transport region
         trgeo = make_transportgeom (llgeo, rlgeo, mgeo)
!!$ Free geometry for 'm' region, we don't need it anymore
         Call free_geom (mgeo)
!!$ Prepare parameters for supercell calculation
         Call prepare_supercell (llgeo, lscell)
         Call prepare_supercell (rlgeo, rscell)
!!$ Calculate S(R)
         msc = calc_screal (trgeo)
         lsc = calc_screal (llgeo)
         rsc = calc_screal (rlgeo)
!!$ If you want to have loop over Energy  do it here
!!$ Just change opt%Eoffset

!!$ Calculate potential parameters
         optspp%irel = opt%irel
         Call calc_potpar (atoms, optspp, opt%Eoffset)
!!$ Prepare potential parameters matrices
!!$ for leads and transport region
         mpar = make_ppar (trgeo, opt%po)
         lpar = make_ppar (llgeo, opt%po)
         rpar = make_ppar (rlgeo, opt%po)
      End If
!!$ Prepare grid for BZ integration

      Call gibz (bz, opt%BZpart, trgeo%base, opt%nk, opt%nsym, opt%inve, opt%iftr, opt%bzone, opt%kpoint)

      If (my_mpi_id == 0) Then
         Write (timef, '("   Init Time = ", f12.4," seconds"/)') MPI_Wtime () - iitime
         flush (timef)
      End If
      consv = 0.0d0

      Call do_log (1, 'Starting transport calculation...')
      Call do_log (1, '     IBZ, spin,  kx,      ky,   cons,nr_lm,nr_rm,ikx,iky,   weight')


!!$  if ((mpi_sz>1) .and. solsize>1 .and. opt%lowmem/=0) then
!!$     nkloops=bz%nkgrid
!!$  else
      nkloops = ceiling (bz%nkgrid/real(mpi_sz/solsize, kind=DEF_DBL_PREC))
!!$  end if


      Call MPI_Barrier (solve_comm, ierr)
      If (opt%onlyspin == 0) Then
         is1 = 1
         is2 = 2
      Else
         is1 = opt%onlyspin
         is2 = opt%onlyspin
      End If

!SMP$ DO SERIAL
      Do loop = 0, nkloops - 1
         iitime = MPI_Wtime ()
         ik = 1 + loop * (mpi_sz/solsize) + my_mpi_id
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
!SMP$ DO SERIAL
         Do is = is1, is2
            havelr = 0
            haverl = 0
            needsolve = 0
            If (havejob /= 0) Then
!!$ Prepare Ando for left lead
               If ((is-opt%lpol) < 2 .Or. opt%lowmem /= 0) Then
                  If (opt%extlead == 2) Then
                     lnm (is) = read_ando (l_ando, ik, is, opt%signat, 'l')
                  Else
                     lnm (is) = get_sc_ando (is, l_ando, k, lpar(is), lsc, lscell, opt%po,-1, opt%gensol, &
                    & opt%specpart)
                     If (opt%extlead == 1) Call write_ando (l_ando, ik, is, opt%signat, 'l')
                  End If
               Else
                  lnm (is) = lnm (is-1)
                  lscell%nreal (:, is) = lscell%nreal(:, is-1)
               End If
!!$ Prepare Ando for Right lead
               If ((is-opt%rpol) < 2 .Or. opt%lowmem /= 0) Then
                  If (opt%extlead == 2) Then
                     rnm (is) = read_ando (r_ando, ik, is, opt%signat, 'r')
                  Else
                     rnm (is) = get_sc_ando (is, r_ando, k, rpar(is), rsc, rscell, opt%po, 1, opt%gensol, &
                    & opt%specpart)
                     If (opt%extlead == 1) Call write_ando (r_ando, ik, is, opt%signat, 'r')
                  End If
               Else
                  rnm (is) = rnm (is-1)
                  rscell%nreal (:, is) = rscell%nreal(:, is-1)
               End If
               If (my_mpi_id == 0) Then
                  Write (timef, '("   Ando Time = ", f12.4," seconds")') MPI_Wtime () - iitime
                  flush (timef)
               End If
               If (l_ando%Nin*opt%needlr > 0) havelr = 1
               If (r_ando%Nin*opt%needrl > 0) haverl = 1
               needsolve = havelr + haverl
               If (needsolve /= 0) Then
!!$ Prepare P-S or H-HOH system for whole device
                  Call prep_system (msys, sk, mpar(is), opt%po, 0)
                  If (my_mpi_id == 0) Then
                     Write (timef, '("   Prep Time = ", f12.4," seconds")') MPI_Wtime () - iitime
                     flush (timef)
                  End If
!!$ Embed boundary conditions and prepare RHS
                  Call embed_boundary (msys%c, rhs, l_ando, r_ando, opt%needlr, opt%needrl)
                  If (my_mpi_id == 0) Then
                     Write (timef, '("  Embed Time = ", f12.4," seconds")') MPI_Wtime () - iitime
                     flush (timef)
                  End If
                  If (opt%lowmem /= 0) Then
                     Call free_ando_emb (l_ando)
                     Call free_ando_emb (r_ando)
                  End If
                  If (opt%writeham == 1) Then
                     Write (cwork, '("ham/",i3.3,"_",i1,".dat")') ik, is
                     Call dump_matrix (msys%c, trim(cwork))
                  End If
                  If (opt%writerhs == 1) Then
                     Write (cwork, '("rhs/",i3.3,"_",i1,".dat")') ik, is
                     Call dump_matrix (rhs%bl, trim(cwork))
                  End If
               End If
            End If
            Call MPI_Barrier (local_comm, ierr)
            Call mpi_allgather (needsolve, 1, mpi_integer, nsols, 1, mpi_integer, local_comm, ierr)
!!$ Solve system
!SMP$ DO SERIAL
            Do boss = 0, solsize - 1
               If (nsols(boss+1) /= 0) Then
                  Call sp_solve (msys%c, rhs, boss, local_comm, opt%leq)
               End If
            End Do
            If (opt%writewf == 1 .And. needsolve > 0) Then
               Write (cwork, '("wf/",i3.3,"_",i1,".dat")') ik, is
               Call dump_matrix (rhs%bl, trim(cwork))
            End If

!!$ Calculate conductance
            If (havejob /= 0) Then
               cons1 = cond_singlespin (is, tran, rhs, l_ando, r_ando, havelr, haverl)
               Write (cwork, '(1x,i6,1x,i4,2(1x,f8.4),5(1x,i4),3x,f8.4,3x,g17.8)') ik, is, k (1), k (2), &
              & floor (Log10(cons1)), lnm (is), rnm (is), bz%ik(1, ik), bz%ik(2, ik), bz%kweight(ik), cons1

               cons = cons + cons1
               If (consw < cons1) consw = cons1
               Call free (rhs)
            End If
            Call MPI_Barrier (solve_comm, ierr)
            Call log_parallel (solve_comm, root, havejob, jobs, cwork)
         End Do

         If (havejob /= 0) Then
!!$ Calculate specular components
            If (opt%do_spec /= 0) Call cond_spec (is1, is2, tran, lscell, rscell, opt%needrl, .True.)

!!$ Write results etc
            If (opt%needlr /= 0) Call write_trans_parallel (trfiles(1), solve_comm, root, havejob, jobs, bz, &
           & ik, tran(1))
            If (opt%needrl /= 0) Call write_trans_parallel (trfiles(2), solve_comm, root, havejob, jobs, bz, &
           & ik, tran(2))

            If (opt%writetm /= 0) Then
               If (opt%needlr /= 0) Call write_trans_matrices ('l', tran(1), ik)
               If (opt%needrl /= 0) Call write_trans_matrices ('r', tran(2), ik)
            End If

            If (opt%lpol == 0 .And. lnm(1) /= 0) Then
               gmix (1) = gmix (1) + bz%kweight(ik) * (lnm(1)-sum(tran(1)%R%mat(1, &
              & 1)%bl*conjg(tran(1)%R%mat(2, 2)%bl)))
!!$                gmix (1) = gmix (1) + bz%kweight(ik) * (lnm(1)-sum(tran(1)%R%mat(1, &
!!$               & 1)%bl*conjg(tran(1)%R%mat(1, 1)%bl)))
!!$                gmix (2) = gmix (1) + bz%kweight(ik) * (lnm(1)-sum(tran(1)%R%mat(2, &
!!$               & 2)%bl*conjg(tran(1)%R%mat(2, 2)%bl)))

               If (opt%rpol == 0 .And. rnm(1) /= 0) Then
                  tmix (1) = tmix (1) + bz%kweight(ik) * (sum(tran(1)%T%mat(1, 1)%bl*conjg(tran(1)%T%mat(2, &
                 & 2)%bl)))
               End If
            End If


            If (opt%rpol == 0 .And. rnm(2) /= 0 .and. haverl/=0) Then
               gmix (2) = gmix (2) + bz%kweight(ik) * (rnm(1)-sum(tran(2)%R%mat(1, &
              & 1)%bl*conjg(tran(2)%R%mat(2, 2)%bl)))
               If (opt%lpol == 0 .And. lnm(2) /= 0) Then
                  tmix (2) = tmix (2) + bz%kweight(ik) * (sum(tran(2)%T%mat(1, 1)%bl*conjg(tran(2)%T%mat(2, &
                 & 2)%bl)))
               End If
            End If

            consv = consv + 0.5d0 * cons * bz%kweight(ik)
            t_lr = t_lr + bz%kweight(ik) * tran(1)%T%val
            r_lr = r_lr + bz%kweight(ik) * tran(1)%R%val
            t_rl = t_rl + bz%kweight(ik) * tran(2)%T%val
            r_rl = r_rl + bz%kweight(ik) * tran(2)%R%val
            l_shc = l_shc + bz%kweight(ik) * lnm
            r_shc = r_shc + bz%kweight(ik) * rnm

            t_lr_spec = t_lr_spec + bz%kweight(ik) * tran(1)%T%spec
            r_lr_spec = r_lr_spec + bz%kweight(ik) * tran(1)%R%spec
            t_rl_spec = t_rl_spec + bz%kweight(ik) * tran(2)%T%spec
            r_rl_spec = r_rl_spec + bz%kweight(ik) * tran(2)%R%spec
         End If

         If (my_mpi_id == 0) Then
            Write (timef, '(i3.3," x KP Time = ", f12.4," seconds")') mpi_sz / solsize, MPI_Wtime () - iitime
            flush (timef)
         End If
      End Do
      Call do_log (1, 'Done!')
      l_shc (2) = l_shc (1+opt%lpol)
      r_shc (2) = r_shc (1+opt%rpol)

      Call write_results (solve_comm, opt, t_lr, r_lr, t_rl, r_rl, t_lr_spec, t_rl_spec, l_shc, r_shc, consv, &
     & consw, tmix, gmix)

      Call MPI_Barrier (solve_comm, ierr)
!!$ End here you loop over E

      If (my_mpi_id == 0) write (timef, '(/"  Total Time = ", f12.4," seconds")') MPI_Wtime () - ttime

      Call mpi_finalize (ierr)


!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop

Contains


      Subroutine write_trans_parallel (trfile, comm, root, havejob, jobs, bz, ki, tm)
         Implicit None
         Type (t_tranrefl) :: tm, tm1
         Integer :: comm, root, trfile
         Integer :: havejob, jobs (:), ki
         Type (bzone) :: bz
!!$ Local
         Integer :: ierr, i, mpi_sz, my_id, ki1
         Integer :: status (MPI_STATUS_SIZE)

         Call mpi_comm_rank (comm, my_id, ierr)
         Call mpi_comm_size (comm, mpi_sz, ierr)
         If (my_id /= root) Then
            If (havejob /= 0) Then
               Call MPI_Send (tm%T%exists, 4, mpi_integer, root, my_id+2001, comm, ierr)
               Call MPI_Send (tm%T%val, 4, mpi_double_precision, root, my_id+3001, comm, ierr)
               Call MPI_Send (tm%R%exists, 4, mpi_integer, root, my_id+4001, comm, ierr)
               Call MPI_Send (tm%R%val, 4, mpi_double_precision, root, my_id+5001, comm, ierr)
               Call MPI_Send (ki, 1, mpi_integer, root, my_id+6001, comm, ierr)
            End If
         Else
            Call write_trans (trfile, bz%k(:, ki), bz%ik(:, ki), tm)
            Do i = 2, mpi_sz
               If (jobs(i) /= 0) Then

                  Call MPI_Recv (tm1%T%exists, 4, mpi_integer, i-1, i+2000, comm, status, ierr)
                  Call MPI_Recv (tm1%T%val, 4, mpi_double_precision, i-1, i+3000, comm, status, ierr)
                  Call MPI_Recv (tm1%R%exists, 4, mpi_integer, i-1, i+4000, comm, status, ierr)
                  Call MPI_Recv (tm1%R%val, 4, mpi_double_precision, i-1, i+5000, comm, status, ierr)
                  Call MPI_Recv (ki1, 4, mpi_integer, i-1, i+6000, comm, status, ierr)

                  Call write_trans (trfile, bz%k(:, ki1), bz%ik(:, ki1), tm1)
               End If
            End Do
            flush (trfile)
         End If
      End Subroutine write_trans_parallel


      Subroutine write_trans (trfile, k, ki, tm)
         Implicit None
         Type (t_tranrefl) :: tm
         Integer :: i, j, ki (:), trfile
         Real (Kind=DEF_DBL_PREC) :: tval (2, 2), rval (2, 2), k (:)

         Do i = 1, 2
            Do j = 1, 2
               If (tm%T%exists(i, j) /= 0) Then
                  tval (i, j) = tm%T%val(i, j)
               Else
#if defined(HAVE_NAN)
                  tval (i, j) = 0.0d0 / 0.0d0
#else
                  tval (i, j) = - 100
#endif
               End If

               If (tm%R%exists(i, j) /= 0) Then
                  rval (i, j) = tm%R%val(i, j)
               Else
#if defined(HAVE_NAN)
                  rval (i, j) = 0.0d0 / 0.0d0
#else
                  rval (i, j) = - 100
#endif
               End If
            End Do
         End Do
         Write (trfile, '(2(1x,i4),2(1x,f13.8),8(1x,g17.8))') ki, k, tval (:, :), rval (:, :)
      End Subroutine write_trans

      Subroutine open_trans_fles (trfiles, opt)
         Implicit None
         Type (t_options) :: opt
         Integer :: trfiles (2)
1111     Format (2 x, "ikx", 2 x, "iky", 8 x, "kx", 12 x, "ky", 11 x, "T_up_up", 9 x, "T_up_down", 9 x, "T_do&
        &wn_up", 8 x, "T_down_down", 10 x, "R_up_up", 9 x, "R_up_down", 9 x, "R_down_up", 8 x, "R_down_down")
         Open (Unit=timef, File='times', Action='write')
         If (opt%needlr /= 0) Then
            Open (Unit=trfiles(1), File='tran_lr', Action='write')
            Write (trfiles(1), 1111)
         End If
         If (opt%needrl /= 0) Then
            Open (Unit=trfiles(2), File='tran_rl', Action='write')
            Write (trfiles(2), 1111)
         End If
      End Subroutine open_trans_fles

      Subroutine write_results (comm, opt, t_lr, r_lr, t_rl, r_rl, t_lr_spec, t_rl_spec, l_shc, r_shc, consv, &
     & consw, tmix, gmix)
!!$ Write results
         Implicit None
         Integer :: ierr, comm, my_id
         Type (t_options) :: opt
         Real (Kind=DEF_DBL_PREC) :: t_lr (2, 2), r_lr (2, 2), t_lr_spec (2, 2)
         Real (Kind=DEF_DBL_PREC) :: t_rl (2, 2), r_rl (2, 2), t_rl_spec (2, 2)
         Real (Kind=DEF_DBL_PREC) :: r_shc (2), l_shc (2)
         Complex (Kind=DEF_DBL_PREC) :: tmix (2), gmix (2)
         Real (Kind=DEF_DBL_PREC) :: consv, consw
!!$Local
         Real (Kind=DEF_DBL_PREC) :: t_lr1 (2, 2) = 0.0d0, r_lr1 (2, 2) = 0.0d0, t_lr1_spec (2, 2) = 0.0d0
         Real (Kind=DEF_DBL_PREC) :: t_rl1 (2, 2) = 0.0d0, r_rl1 (2, 2) = 0.0d0, t_rl1_spec (2, 2) = 0.0d0
         Real (Kind=DEF_DBL_PREC) :: r_shc1 (2) = 0.0d0, l_shc1 (2) = 0.0d0
         Complex (Kind=DEF_DBL_PREC) :: tmix1 (2) = 0.0d0, gmix1 (2) = 0.0d0
         Real (Kind=DEF_DBL_PREC) :: consv1, consw1
         Integer, Parameter :: fl = 200

         Call mpi_comm_rank (comm, my_id, ierr)

         Call mpi_reduce (t_lr, t_lr1, 4, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (r_lr, r_lr1, 4, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (r_rl, r_rl1, 4, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (t_rl, t_rl1, 4, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (l_shc, l_shc1, 2, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (r_shc, r_shc1, 2, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (consv, consv1, 1, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (consw, consw1, 1, mpi_double_precision, mpi_max, 0, comm, ierr)

         Call mpi_reduce (t_lr_spec, t_lr1_spec, 4, mpi_double_precision, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (t_rl_spec, t_rl1_spec, 4, mpi_double_precision, mpi_sum, 0, comm, ierr)

         Call mpi_reduce (tmix, tmix1, 2, MPI_DOUBLE_COMPLEX, mpi_sum, 0, comm, ierr)
         Call mpi_reduce (gmix, gmix1, 2, MPI_DOUBLE_COMPLEX, mpi_sum, 0, comm, ierr)

!!$ Output results into "andout"
         If (my_id == 0) Then
            Open (Unit=fl, File='andout', Action='write')

            Write (fl,*) ' Conservation: integrated,      worst'
            Write (fl, '("         ",2(5x,g13.6))') consv1, consw1

            Write (fl,*) ' Left Sharvin conductance :'
            Write (fl, '("     by spin :",2(3x,f15.10))') l_shc1
            Write (fl, '("     total   :",3x,f15.10)') sum (l_shc1)
            Write (fl,*) ' Right Sharvin conductance :'
            Write (fl, '("     by spin :",2(3x,f15.10))') r_shc1
            Write (fl, '("     total   :",3x,f15.10)') sum (r_shc1)
            Write (fl,*) '---------------------------------------------'
            Write (fl,*) ''

            If (opt%needlr /= 0) Then
               Write (fl,*) ' L->R Conductance:'
               Write (fl,*) '              Up:                Down:'
               Write (fl, '("   Up   :",2(3x,g17.10))') t_lr1 (1, :)
               Write (fl, '("   Down :",2(3x,g17.10))') t_lr1 (2, :)
               Write (fl,*) ''
               Write (fl,*) ' Conductance total:'
               Write (fl, '("         ",2(3x,g17.10))') sum (t_lr1(1, :)), sum (t_lr1(2, :))
               Write (fl, '("  Grand total:                  ",g17.10)') sum (t_lr1(:, :))
               If (opt%do_spec /= 0) Then
                  Write (fl,*) ''
                  Write (fl,*) ' Specular:'
                  Write (fl,*) '              Up:                Down:'
                  Write (fl, '("   Up   :",2(3x,g17.10))') t_lr1_spec (1, :)
                  Write (fl, '("   Down :",2(3x,g17.10))') t_lr1_spec (2, :)
               End If
               Write (fl,*) ''
               Write (fl,*) ' Mixing conductance (Re, Im):'
               Write (fl, '("         ",2(3x,g17.10))') real (gmix1(1)), imag (gmix1(1))
               Write (fl,*) ' Mixing transmission (Re, Im):'
               Write (fl, '("         ",2(3x,g17.10))') real (tmix1(1)), imag (tmix1(1))
               Write (fl,*) '---------------------------------------------'
               Write (fl,*) ''
            End If

            If (opt%needrl /= 0) Then
               Write (fl,*) ' R->L Conductance:'
               Write (fl,*) '              Up                Down'
               Write (fl, '("   Up   :",2(3x,g17.10))') t_rl1 (1, :)
               Write (fl, '("   Down :",2(3x,g17.10))') t_rl1 (2, :)
               Write (fl,*) ''
               Write (fl,*) ' Conductance total:'
               Write (fl, '("         ",2(3x,g17.10))') sum (t_rl1(1, :)), sum (t_rl1(2, :))
               Write (fl, '("  Grand total:                  ",g17.10)') sum (t_rl1(:, :))
               If (opt%do_spec /= 0) Then
                  Write (fl,*) ''
                  Write (fl,*) ' Specular:'
                  Write (fl,*) '              Up:                Down:'
                  Write (fl, '("   Up   :",2(3x,g17.10))') t_rl1_spec (1, :)
                  Write (fl, '("   Down :",2(3x,g17.10))') t_rl1_spec (2, :)
               End If
               Write (fl,*) ''
               Write (fl,*) ' Mixing conductance (Re, Im):'
               Write (fl, '("         ",2(3x,g17.10))') real (gmix1(2)), imag (gmix1(2))
               Write (fl,*) ' Mixing transmission (Re, Im):'
               Write (fl, '("         ",2(3x,g17.10))') real (tmix1(2)), imag (tmix1(2))
               Write (fl,*) '---------------------------------------------'
               Write (fl,*) ''
            End If

            Close (fl)
         End If
      End Subroutine write_results

      Subroutine log_parallel (comm, root, havejob, jobs, cwork)
!!$ Log in parallel something
         Implicit None
         Integer :: root, havejob, jobs (:)
         Character (Len=clen) :: cwork
!!$Local
         Integer :: ierr, procnum, comm, mpi_sz, my_mpi_id
         Integer :: status (MPI_STATUS_SIZE)
         Call mpi_comm_rank (comm, my_mpi_id, ierr)
         Call mpi_comm_size (comm, mpi_sz, ierr)
         If (my_mpi_id /= root) Then
            If (havejob /= 0) Then
               Call MPI_Send (cwork, clen, MPI_CHARACTER, root, my_mpi_id+mpi_sz+1000, comm, ierr)
            End If
         Else
            Call do_log (1, trim(cwork))
            Do procnum = 2, mpi_sz
               If (jobs(procnum) > 0) Then
                  Call MPI_Recv (cwork, clen, MPI_CHARACTER, procnum-1, procnum+mpi_sz+999, comm, status, &
                 & ierr)
                  Call do_log (1, trim(cwork))
               End If
            End Do
         End If
      End Subroutine log_parallel
End Program trans
!!$  call dump_matrix(sk%m,'sk.dat')
