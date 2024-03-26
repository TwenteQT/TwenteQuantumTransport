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

Program hdfaddmoms
      Use transport
      Use atoms_module
      Use geometry_module
      Use supercell
      Use structure
      Use logging
      Use rotations
      Use readcfg
      Use hdf5io
      Use compat
      Use helpers

      Implicit None
      Include 'mpif.h'

!!$ Local
      Integer, Parameter :: clen = 200
      Type (t_options) :: opt
      Type (t_atoms_set) :: atoms
      Type (t_geometry) :: mgeo, lgeo, rgeo, trgeo, llgeo, rlgeo
      Type (t_supercell) :: lscell, rscell
      Real (Kind=DEF_DBL_PREC), Pointer :: rotm (:, :)
      Integer :: hops
      Real (Kind=DEF_DBL_PREC), External :: dznrm2
!!$       Type (t_mask) :: mask
!!$       Type (zcsrmat) :: tempmat2
      Integer (HID_T) :: h5out
      Integer, Parameter :: ifile = 3944, ofile = 34345, timef = 339
!!$ MPI vars
      Integer :: ierr, my_mpi_id, solve_comm, mpi_sz, root = 0
      Integer :: error
      Real (Kind=DEF_DBL_PREC) :: magvals (3), magmom (3),mdir(2),zr2(2)=0.0d0
      
!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------

      Call init_random_seed ()
      Call mpi_init (ierr)
      solve_comm = mpi_comm_world
      Call mpi_comm_rank (solve_comm, my_mpi_id, ierr)
      Call mpi_comm_size (solve_comm, mpi_sz, ierr)

      Log_Level = - 1
!!$ Read  file with parameters


!!$ Set LogLevel and open output files
      If (my_mpi_id == root) Then
          Call trans_config_init (opt)
          opt%po%full = 1

         Call hdf5lib_init ()
         Call h5fopen_f ('transout.h5', H5F_ACC_RDWR_F, h5out, error)
         if (error/=0) then
             write(*,*) 'transout.h5 does not exist!'
             Call h5close_f (error)
             stop
         end if
         call hdfread(h5out, '/compat', compat_flags,100)
         if (compat_flags (cpf_hasmagmom) == 1) then
             write(*,*) 'file already has moments!'
             Call h5fclose_f (h5out, error)
             Call h5close_f (error)
             stop
         end if

         hops = 1
         If (opt%po%kind == 3) hops = 2

!!$ Read atomic potentials
         atoms = read_atoms ('atomlist')
!!$ Read geometries
         lgeo = read_geom ('geom_l', atoms)
         rgeo = read_geom ('geom_r', atoms)
         mgeo = read_geom ('geom_m', atoms,lgeo%scale)
         Call hdfwrite (h5out, 'geom/mgeo', mgeo)
         Call hdfwrite (h5out, 'geom/rgeo', rgeo)
         Call hdfwrite (h5out, 'geom/lgeo', lgeo)

!!$ Prepare geometry for leads
         llgeo = make_leadgeom (lgeo, hops,-1)
         rlgeo = make_leadgeom (rgeo, hops, 1)
!!$ Free original geometry for leads, we don't need it anymore
         Call free_geom (lgeo)
         Call free_geom (rgeo)
!!$ Prepare geometry for transport region
         trgeo = make_transportgeom (llgeo, rlgeo, mgeo)
         Call hdfwrite (h5out, 'geom/trgeo', trgeo)
         Call hdfwrite (h5out, 'geom/llgeo', llgeo)
         Call hdfwrite (h5out, 'geom/rlgeo', rlgeo)

!!$ Free geometry for 'm' region, we don't need it anymore
         Call free_geom (mgeo)
!!$ Prepare parameters for supercell calculation
         Call prepare_supercell (llgeo, lscell, 1, split=1)
         Call prepare_supercell (rlgeo, rscell, 1, split=1)

         rotm => prep_rot_mask (trgeo, opt%po)
         
         mdir = rotspheric(zr2,opt%po%rot_angle)
         
         If (calc_magmoms(trgeo, rotm, magvals, magmom)) Then
            compat_flags (cpf_hasmagmom) = 1
            If (opt%damp%axesdir /= 0) Then
               mdir (1) = acos(magmom(3)/magvals(2))
               mdir (2) = atan2(magmom(2)/magvals(2),magmom(1)/magvals(2))
            End If
            call hdfwrite(h5out,'/damp/mdir',mdir)
            Call hdfwrite (h5out, '/geom/trgeo/rotm', rotm)
            Call hdfwrite (h5out, '/geom/mgeo/rotm', rotm(:,(trgeo%ltna+1):(trgeo%ltna+trgeo%tna)))
            Call hdfwrite (h5out, '/damp/magvals', magvals)
            Call hdfwrite (h5out, '/damp/magmom', magmom)
            Call hdfwrite (h5out, '/compat', compat_flags)
         else
             write(*,*) 'did not find moments in atomic files!'
         End If
         

         Call h5fclose_f (h5out, error)
         Call h5close_f (error)

      End If
      Call MPI_Barrier (solve_comm, ierr)
      Call mpi_finalize (ierr)

!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop
End Program hdfaddmoms
