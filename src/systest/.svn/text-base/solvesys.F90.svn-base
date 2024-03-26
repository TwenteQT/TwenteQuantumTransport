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

Program testh
      Use sparse_solvers
      Use HDF5io
      Use sparselib
      Use sparse_solvers
      Use readcfg
      Implicit None
      Include 'mpif.h'
!!$
      Type (zcsrmat) :: ham
      Type (zcsrmat) :: rhs
      Type (zdensemat) :: wf
      Integer (HID_T) :: h5in, h5out
      Integer :: error
      Type (t_leq_param) :: opt
      Integer :: ierr, my_mpi_id, solve_comm, mpi_sz, root = 0
      Real (Kind=DEF_DBL_PREC) :: tm1, tm2
      Character (Len=1024) :: cwork


      Call h5open_f (error)
      Call h5fopen_f ('transout.h5', H5F_ACC_RDONLY_F, h5in, error)

      Call mpi_init (ierr)
      solve_comm = mpi_comm_world
      Call mpi_comm_rank (solve_comm, my_mpi_id, ierr)
      Call mpi_comm_size (solve_comm, mpi_sz, ierr)

      If (my_mpi_id == root) Then
         Log_Level = 1
      End If
      Call readconf (opt)
      If (my_mpi_id == root) Then
         Log_Level = 1
         Call h5open_f (error)
         Call h5fopen_f ('transout.h5', H5F_ACC_RDONLY_F, h5in, error)
         Call hdfread (h5in, '/kset/1/ham', ham)
         Call hdfread (h5in, '/kset/1/rhs', rhs)
         write(*,*) rhs%ncol,rhs%nrow
         Call do_log (1, 'Data is loaded!', .True.)
      End If

      opt%solver = 0
      opt%debug = 2
      Call MPI_Barrier (solve_comm, ierr)
      tm1 = MPI_Wtime ()
      Call sp_solve_srhs (ham, rhs, wf, root, solve_comm, opt)
      Write (cwork, '("Solve time=",g,"sec"))') MPI_Wtime () - tm1
      If (my_mpi_id == root) Then
         Call h5fcreate_f ('matout.h5', H5F_ACC_TRUNC_F, h5out, error)
         Call hdfsetversion (h5out)
         Call hdfwrite (h5out, 'wavefunc', wf%bl)
         Call h5fclose_f (h5out, error)
         Call h5close_f (error)
      End If
      Call do_log (1, trim(cwork), .True.)
      Call MPI_Barrier (solve_comm, ierr)
      Call mpi_finalize (ierr)
Contains
      Subroutine readconf (opt)
         Type (t_leq_param), Intent (Inout) :: opt
!!$
         Type (t_config_option), Pointer :: optlist
         Call read_config_file ('test.conf', optlist)

         opt%memover = getopt_i (optlist, 'mumps_memory_overhead',-1)
         opt%partitioner = getopt_i (optlist, 'mumps_partitioner', 4)
         opt%redistmtx = getopt_i (optlist, 'mumps_redist_mtx', 0)
         opt%permscale = getopt_i (optlist, 'mumps_permscale', 7)
         opt%workhost = getopt_i (optlist, 'mumps_workhost', 1)
         opt%debug = getopt_i (optlist, 'mumps_debug', 1)
         Call do_log (1, '# For OOC solver set variable MUMPS_TMP_DIR. Then OOC solver will be selected')
         opt%tmpdir = getopt_s (optlist, 'MUMPS_TMP_DIR', '')
      End Subroutine readconf

End Program testh
