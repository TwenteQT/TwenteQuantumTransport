#include "math_def.h"

Program fermproj
      Use ando_module
      Use atoms_module
      Use geometry_module
      Use structure
      Use supercell
      Use hamiltonian
      Use logging
      Use bzgrid
      Use rotations
      Use readcfg
      Use mpi
      Implicit None
!      Include 'mpif.h'
!!$ Local
      Integer, Parameter :: clen = 200, timef = 201, fermf = 284645
      Integer, Parameter :: trfiles (2) = (/ 202, 203 /)
      Type (t_options) :: opt
      Type (t_atoms_set) :: atoms
      Type (t_geometry) :: bgeo, fbgeo
      Type (t_strconst) :: fsc
      Type (t_ando_sollution) :: ando
      Type (t_potpar_mats) :: fpar (2)
      Type (t_supercell) :: scell
      Type (bzone) :: bz
      Real (Kind=DEF_DBL_PREC) :: iitime, ttime
      Integer :: ik, nkloops, loop
      Real (Kind=DEF_DBL_PREC) :: k (2)
      Real (Kind=DEF_DBL_PREC) :: shc (2) = 0.0d0
      Real (Kind=DEF_DBL_PREC), Pointer :: rotm (:, :)
      Integer :: fnm (2) = 0, hops !, i, j, s
      Character (Len=clen) :: cwork, cwork1
      Real (Kind=DEF_DBL_PREC), External :: dznrm2
      Integer, Pointer :: jobs (:)
!!$ MPI vars
      Integer :: ierr, my_mpi_id, local_comm, solve_comm, mpi_sz, havejob, root = 0, solsize
      Integer :: mpi_loc_id
      Integer :: ifull = 0
!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------


!!$      Call check_li_size()

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

      Call MPI_COMM_SPLIT (solve_comm, my_mpi_id, 0, local_comm, ierr)
      solsize = 1

      Call mpi_comm_rank (local_comm, mpi_loc_id, ierr)

!!$ Set LogLevel and open output files
      Allocate (jobs(mpi_sz))
      If (my_mpi_id == root) Then
         Log_Level = opt%loglvl
         Open (Unit=timef, File='times', Action='write')
         Open (Unit=fermf, File='fermproj', Action='write')
      Else
!!$ uncomment line bellow if you want output from rest of CPUs
!!$     Log_Level=opt%loglvl
         Write (cwork, '(i2.2)') my_mpi_id
         Write (logfilename, '("log.",A)') trim (cwork)
      End If


      ttime = MPI_Wtime ()
      iitime = MPI_Wtime ()

      hops = 1
      If (opt%po%kind > 2) hops = 2
      If ((opt%po%so /= 0) .Or. (Abs(opt%po%rot_angle(1)) > 1.0d-8) .Or. (Abs(opt%po%rot_angle(2)) > 1.0d-8)) &
     & Then
         ifull = 1
         opt%po%full = 1
      End If

      If (mpi_loc_id == 0) Then
!!$ Read atomic potentials
         atoms = read_atoms ('atomlist')
!!$ Read geometries
         bgeo = read_geom ('geom_f', atoms)
!!$ Prepare geometry for leads
         fbgeo = make_leadgeom (bgeo, hops,-1)
!!$ Free original geometry for leads, we don't need it anymore
         Call free_geom (bgeo)
!!$ Prepare parameters for supercell calculation
         Call prepare_supercell (fbgeo, scell, 1, split=ifull)

         rotm => prep_rot_mask (fbgeo, opt%po)

!!$ Calculate S(R)
         fsc = calc_screal (fbgeo, ifull)
!!$ If you want to have loop over Energy  do it here
!!$ Just change opt%Eoffset

!!$ Calculate potential parameters
         Call calc_potpar (atoms, opt%atopt, opt%Eoffset)
!!$ Prepare potential parameters matrices
!!$ for leads and transport region
         If (ifull /= 0) Then
            fpar (1) = make_ppar_full (fbgeo, opt%po)
         Else
            fpar = make_ppar (fbgeo, opt%po)
         End If
      End If
!!$ Prepare grid for BZ integration

      Call gibz (bz, fbgeo%base, opt%bzo)

      If (my_mpi_id == 0) Then
         Write (timef, '("   Init Time = ", f12.4," seconds"/)') MPI_Wtime () - iitime
         flush (timef)
      End If
      Write (cwork, '("Starting fermi proj. calculation  for FS=",i1)') ifull
      Call do_log (1, trim(cwork))
      Call do_log (1, '      IBZ ,    kx,     ky,    nm, ikx , iky ,   weight')

      nkloops = ceiling (bz%nkgrid/real(mpi_sz/solsize, kind=DEF_DBL_PREC))

      Call MPI_Barrier (solve_comm, ierr)
!SMP$ DO SERIAL
      Do loop = 0, nkloops - 1
         iitime = MPI_Wtime ()
         ik = 1 + loop * (mpi_sz/solsize) + my_mpi_id
         havejob = 0
         If (ik <= bz%nkgrid) Then
            k = bz%k (:, ik)
            havejob = 1
            If (ifull == 1) Then
               fnm = get_sc_ando (1, ando, k, fpar(1), fsc, scell, opt%po,-1, opt%actl)
               fnm = ando%SNin
            Else
               fnm (1) = get_sc_ando (1, ando, k, fpar(1), fsc, scell, opt%po,-1, opt%actl)
               fnm (2) = get_sc_ando (2, ando, k, fpar(2), fsc, scell, opt%po,-1, opt%actl)
            End If

            shc = shc + bz%kweight(ik) * fnm
            Write (cwork, '(1x,i8,2(1x,f8.4),1x,i4,2(1x,i5),3x,f8.4)') ik, k (1), k (2), fnm (1) + fnm (2), &
           & bz%ik(1, ik), bz%ik(2, ik), bz%kweight(ik)
            Write (cwork1, '(2(1x,f12.8),2(1x,i4))') k (1), k (2), fnm (1), fnm (2)
         End If

         Call MPI_Barrier (solve_comm, ierr)
         Call mpi_gather (havejob, 1, mpi_integer, jobs, 1, mpi_integer, root, solve_comm, ierr)
         Call log_bz (solve_comm, root, havejob, jobs, cwork, cwork1, fermf)

         If (my_mpi_id == 0) Then
            Write (timef, '(i3.3," x KP Time = ", f12.4," seconds")') mpi_sz / solsize, MPI_Wtime () - iitime
            flush (timef)
         End If
      End Do
      Call do_log (1, 'Done!')
!!$       dens(:)=sqrt(dens(:))
      Call write_results (solve_comm, shc)

      Call MPI_Barrier (solve_comm, ierr)
!!$ End here you loop over E

      If (my_mpi_id == 0) Then
         Write (timef, '(/"  Total Time = ", f12.4," seconds")') MPI_Wtime () - ttime
         Close (Unit=timef)
         Close (Unit=fermf)
      End If

      Call mpi_finalize (ierr)


!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop
Contains
      Subroutine write_results (comm, shc)
!!$ Write results
         Implicit None
         Integer :: ierr, comm, my_id
         Real (Kind=DEF_DBL_PREC) :: shc (2)
!!$Local
         Real (Kind=DEF_DBL_PREC) :: shc1 (2) = 0.0d0
         Integer, Parameter :: fl = 200

         Call mpi_comm_rank (comm, my_id, ierr)

         Call mpi_reduce (shc, shc1, 2, mpi_double_precision, mpi_sum, 0, comm, ierr)

!!$ Output results into "andout"
         If (my_id == 0) Then
            Open (Unit=fl, File='andout', Action='write')

            Write (fl,*) ' Conservation: integrated,      worst'
            Write (fl, '("         ",2(5x,g13.6))') 0.0d0, 0.0d0

            Write (fl,*) ' Bulk Sharvin conductance :'
            Write (fl, '("         ",2(3x,f15.10))') shc1

            Close (fl)
         End If
      End Subroutine write_results

      Subroutine log_bz (comm, root, havejob, jobs, cwork, cwork1, fermf)
!!$ Log in parallel something
         Implicit None
         Integer :: root, havejob, jobs (:), fermf
         Character (Len=clen) :: cwork, cwork1
!!$Local
         Integer :: ierr, procnum, comm, mpi_sz, my_mpi_id
         Integer :: status (MPI_STATUS_SIZE)
         Call mpi_comm_rank (comm, my_mpi_id, ierr)
         Call mpi_comm_size (comm, mpi_sz, ierr)
         If (my_mpi_id /= root) Then
            If (havejob /= 0) Then
               Call MPI_Send (cwork, clen, MPI_CHARACTER, root, my_mpi_id+mpi_sz+1000, comm, ierr)
               Call MPI_Send (cwork1, clen, MPI_CHARACTER, root, my_mpi_id+mpi_sz+2000, comm, ierr)
            End If
         Else
            Call do_log (1, trim(cwork))
            Write (fermf, '(A)') trim (cwork1)
            Do procnum = 2, mpi_sz
               If (jobs(procnum) > 0) Then
                  Call MPI_Recv (cwork, clen, MPI_CHARACTER, procnum-1, procnum+mpi_sz+999, comm, status, &
                 & ierr)
                  Call MPI_Recv (cwork1, clen, MPI_CHARACTER, procnum-1, procnum+mpi_sz+1999, comm, status, &
                 & ierr)
                  Call do_log (1, trim(cwork))
                  Write (fermf, '(A)') trim (cwork1)
               End If
            End Do
         End If
      End Subroutine log_bz

End Program fermproj
