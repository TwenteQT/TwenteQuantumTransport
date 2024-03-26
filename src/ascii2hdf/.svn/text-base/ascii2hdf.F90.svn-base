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

Program ascii2hdf
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

      Implicit None
      Include 'mpif.h'

!!$ Local
      Integer, Parameter :: clen = 200
      Type (t_options) :: opt
      Type (t_atoms_set) :: atoms
      Type (t_geometry) :: mgeo, lgeo, rgeo, trgeo, llgeo, rlgeo
      Type (t_mathop) :: msys, sk, msys1
      Type (t_supercell) :: lscell, rscell
      Type (bzone) :: bz
      Real (Kind=DEF_DBL_PREC), Pointer :: rotm (:, :)
      Integer :: hops
      Character (Len=clen) :: cwork
      Real (Kind=DEF_DBL_PREC), External :: dznrm2
!!$       Type (t_mask) :: mask
!!$       Type (zcsrmat) :: tempmat2
      Integer (HID_T) :: h5out
      Integer, Parameter :: ifile = 3944, ofile = 34345, timef = 339
      Real (Kind=DEF_DBL_PREC), Allocatable :: tt (:, :), rr (:, :), cons (:), tt1 (:, :), rr1 (:, :), cons1 &
     & (:)
      Complex (Kind=DEF_DBL_PREC) :: dc (2, 2)
      Integer, Allocatable :: nmodes (:, :)
      Integer :: nk, i, st, nk1
      Real (Kind=DEF_DBL_PREC) :: cond (2, 2), tm (4), r1, r2, dweight
      Logical :: flag, lnew,exist
      Integer :: islr, isrl, isdamp, ddir, dd1, dd2, dd0, dind
!!$ MPI vars
      Integer :: ierr, my_mpi_id, solve_comm, mpi_sz, root = 0
      Integer :: error
      Real (Kind=DEF_DBL_PREC) :: magvals (3), magmom (3),mdir(2),zr2(2)=0.0d0
      
!!$ Porgram starts here
!!$ -----------------------------------------------------------------------------

      inquire(file='andout', exist=exist) 
      if (.not. exist) then
          write(*,*) 'No andout file, nothing to do'
          stop
      end if
    
      Call init_random_seed ()
      Call mpi_init (ierr)
      solve_comm = mpi_comm_world
      Call mpi_comm_rank (solve_comm, my_mpi_id, ierr)
      Call mpi_comm_size (solve_comm, mpi_sz, ierr)

      Log_Level = - 1
!!$ Read  file with parameters

      Call trans_config_init (opt)
      opt%po%full = 1
      
!!$ Set LogLevel and open output files
      If (my_mpi_id == root) Then

         Open (Unit=ifile, File='output', Action='read')
         Read (ifile, '(A)') cwork
         If (index(cwork, 'HDF') > 0) Then
            Write (*,*) 'Results produced by HDF version of code!'
            Call mpi_finalize (ierr)
            Stop
         End If

         Open (Unit=ofile, File='output.hdf', Action='write')
         cwork = trim (cwork) // ' (HDF)'
         Write (ofile, '(A)') trim (cwork)

         isrl = 0
         islr = 0
         ddir = 0
         flag = .False.
         Do while ( .Not. flag)
            Read (ifile, '(A)', IoStat=st) cwork
            Write (ofile, '(A)') trim (cwork)
            If (get_i(cwork, 'DO_LEFT_TO_RIGHT', i)) islr = i
            If (get_i(cwork, 'DO_RIGHT_TO_LEFT', i)) isrl = i
            If (get_i(cwork, 'DAMP_DIR', i)) ddir = i
            If (st /= 0) Then
               Write (*,*) 'End of output file reached, wrong input!'
               Stop
            End If
            flag = (index(cwork, 'Starting transport calculation') > 0)
         End Do
         Close (ifile)

         isdamp = 0
         Open (Unit=ifile, File='andout', Action='read')
         flag = .True.
         Do while (flag)
            Read (ifile, '(A)', IoStat=st) cwork
            isdamp = index (cwork, 'Damping')
            flag = ((st == 0) .And. (isdamp < 1))
         End Do
         Close (ifile)
         If (isdamp > 0) isdamp = 1

         compat_flags (cpf_base) = 1
         compat_flags (cpf_no_spin_separation) = 1
         compat_flags (cpf_hasdamping) = isdamp
         compat_flags (cpf_converted) = 1

!!$          Write (*,*) islr, isrl
         Call hdf5lib_init ()
         Call h5fcreate_f ('transout.h5', H5F_ACC_TRUNC_F, h5out, error)
         Call hdfsetversion (h5out)
         Call hdfwrite (h5out, '/opts', opt)

         Call init_mathop (msys)
         Call init_mathop (msys1)
         Call init_mathop (sk)
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
         End If
         call hdfwrite(h5out,'/damp/mdir',mdir)
         Call hdfwrite (h5out, '/geom/trgeo/rotm', rotm)
         Call hdfwrite (h5out, '/geom/mgeo/rotm', rotm(:,(trgeo%ltna+1):(trgeo%ltna+trgeo%tna)))
         Call hdfwrite (h5out, '/damp/magvals', magvals)
         Call hdfwrite (h5out, '/damp/magmom', magmom)

         Call gibz (bz, trgeo%base, opt%bzo)

         If (ddir == 0) Then
            dd1 = 1
            dd2 = 2
            dweight = 0.25d0
         Else
            dd1 = ddir
            dd2 = ddir
            dweight = 0.5d0
         End If


         Call hdfwrite (h5out, 'bz', bz)

         nk = bz%nkgrid
         Write (*, '("starting conversion for nk=",i)') nk
         nk1 = bz%nkgrid * (1+isdamp*(1+2*(dd2-dd1)))
         Allocate (rr(4, nk1), tt(4, nk1), cons(nk1), nmodes(4, nk1))
         If (isdamp /= 0) Allocate (rr1(4, nk), tt1(4, nk), cons1(nk))


         nmodes (:, :) = 0

         If (islr > 0) Then
            If (rtr('tran_lr', nk1, tt, rr, nmodes(1:2, :))) Then
               If (isdamp == 0) Then
                  Call hdfwrite (h5out, '/trans/lr', tt)
                  Call hdfwrite (h5out, '/trans/ll', rr)
               Else
                  i = 0
                  tt1 (:, :) = 0.0d0
                  rr1 (:, :) = 0.0d0
                  Do dd0 = dd1, dd2
                     Do dind = 1, 2
                        Write (cwork(1:1), '(i1)') dd0
                        Write (cwork(2:2), '(i1)') dind
                        Call hdfwrite (h5out, '/damp/'//cwork(1:1)//'/'//cwork(2:2)//'/trans/lr', tt(:, &
                       & i*nk+1:i*nk+nk))
                        Call hdfwrite (h5out, '/damp/'//cwork(1:1)//'/'//cwork(2:2)//'/trans/ll', rr(:, &
                       & i*nk+1:i*nk+nk))
                        tt1 (:, :) = tt1 (:, :) + tt (:, i*nk+1:i*nk+nk)
                        rr1 (:, :) = rr1 (:, :) + rr (:, i*nk+1:i*nk+nk)
                        i = i + 1
                     End Do
                  End Do
                  tt1 = tt1 * dweight
                  rr1 = rr1 * dweight
                  Call hdfwrite (h5out, '/trans/lr', tt1)
                  Call hdfwrite (h5out, '/trans/ll', rr1)

               End If
               Write (*,*) 'tran_lr'
            Else
               Write (*,*) 'error reading tran_lr!'
               Stop
            End If
         End If

         If (isrl > 0) Then
            If (rtr('tran_rl', nk1, tt, rr, nmodes(3:4, :))) Then
               If (isdamp == 0) Then
                  Call hdfwrite (h5out, '/trans/rl', tt)
                  Call hdfwrite (h5out, '/trans/rr', rr)
               Else
                  i = 0
                  tt1 (:, :) = 0.0d0
                  rr1 (:, :) = 0.0d0
                  Do dd0 = dd1, dd2
                     Do dind = 1, 2
                        Write (cwork(1:1), '(i1)') dd0
                        Write (cwork(2:2), '(i1)') dind
                        Call hdfwrite (h5out, '/damp/'//cwork(1:1)//'/'//cwork(2:2)//'/trans/rl', tt(:, &
                       & i*nk+1:i*nk+nk))
                        Call hdfwrite (h5out, '/damp/'//cwork(1:1)//'/'//cwork(2:2)//'/trans/rr', rr(:, &
                       & i*nk+1:i*nk+nk))
                        tt1 (:, :) = tt1 (:, :) + tt (:, i*nk+1:i*nk+nk)
                        rr1 (:, :) = rr1 (:, :) + rr (:, i*nk+1:i*nk+nk)
                        i = i + 1
                     End Do
                  End Do
                  tt1 = tt1 * dweight
                  rr1 = rr1 * dweight
                  Call hdfwrite (h5out, '/trans/rl', tt1)
                  Call hdfwrite (h5out, '/trans/rr', rr1)
               End If
               Write (*,*) 'tran_rl'
            Else
               Write (*,*) 'error reading tran_rl!'
               Stop
            End If
         End If

         If (rcons(cons, nk1, nmodes)) Then
            If (isdamp == 0) Then
               Call hdfwrite (h5out, '/trans/cons', cons)
            Else
               i = 0
               cons1 (:) = 0.0d0
               Do dd0 = dd1, dd2
                  Do dind = 1, 2
                     Write (cwork(1:1), '(i1)') dd0
                     Write (cwork(2:2), '(i1)') dind
                     Call hdfwrite (h5out, '/damp/'//cwork(1:1)//'/'//cwork(2:2)//'/trans/cons', &
                    & cons(i*nk+1:i*nk+nk))
                     cons1 (:) = cons1 (:) + cons (i*nk+1:i*nk+nk)
                     i = i + 1
                  End Do
               End Do
               cons1 = cons1 * dweight
               Call hdfwrite (h5out, '/trans/cons', cons1)
            End If
            Write (*,*) 'conservation'
         Else
            Write (*,*) 'problem with conservation data!'
            Stop
         End If


         Call hdfwrite (h5out, '/trans/nmodes', nmodes(:, 1:nk))

         If (isdamp /= 0) Then
            Do dd0 = dd1, dd2
               Do dind = 1, 2
                  Write (cwork(1:1), '(i1)') dd0
                  Write (cwork(2:2), '(i1)') dind
                  Call hdfwrite (h5out, '/damp/'//cwork(1:1)//'/'//cwork(2:2)//'/trans/nmodes', nmodes(:, &
                 & 1:nk))
               End Do
            End Do
         End If

         Open (Unit=ifile, File='andout', Action='read')
         Read (ifile, '(A)') cwork
         Read (ifile, '("         ",2(5x,g13.6))') r1, r2
         Call hdfwrite (h5out, '/cond/consint', r1)
         Call hdfwrite (h5out, '/cond/conspeak', r2)

         If (rsharv(ifile, 'Left Sharvin', r1, r2)) Call hdfwrite (h5out, '/cond/l_shc', (/ r1, r2 /))
         If (rsharv(ifile, 'Right Sharvin', r1, r2)) Call hdfwrite (h5out, '/cond/r_shc', (/ r1, r2 /))

         If (rcond(ifile, 'L->R', cond, lnew)) Then
            compat_flags (cpf_haslr) = 1
         End If

         Call hdfwrite (h5out, '/cond/lr', cond)
         Call hdfwrite (h5out, '/cond/totlr', sum(cond))

         If (rcond(ifile, 'R->L', cond, lnew)) Then
            compat_flags (cpf_hasrl) = 1
         End If

         Call hdfwrite (h5out, '/cond/rl', cond)
         Call hdfwrite (h5out, '/cond/totrl', sum(cond))

         If (rcond(ifile, 'L->L', cond, lnew)) Then
         End If

         Call hdfwrite (h5out, '/cond/ll', cond)
         Call hdfwrite (h5out, '/cond/totll', sum(cond))

         If (rcond(ifile, 'R->R', cond, lnew)) Then
         End If

         Call hdfwrite (h5out, '/cond/rr', cond)
         Call hdfwrite (h5out, '/cond/totrr', sum(cond))

         If (lnew) Then
            compat_flags (cpf_reliabledecompose) = 1
         End If
         cond (:, :) = 0.0d0

         Call hdfwrite (h5out, '/cond/lr_spec', cond)
         Call hdfwrite (h5out, '/cond/rl_spec', cond)
         Call hdfwrite (h5out, '/cond/ll_spec', cond)
         Call hdfwrite (h5out, '/cond/rr_spec', cond)

         dc (:, :) = dcmplx (0.0d0, 0.0d0)

         Call hdfwrite (h5out, '/cond/gmix', dc(1:2, 1))
         Call hdfwrite (h5out, '/cond/tmix', dc(1:2, 1))


         If (isdamp /= 0) Then
            If (rdamp(ifile, dc)) Then
               Call hdfwrite (h5out, '/damp/total', dc)
            Else
               Write (*,*) 'Can not read damping from andout!'
               Stop
            End If
         End If

         Close (ifile)




         If (rtms(tm, r1, r2, i)) Then
            Write (ofile, '("   Init Time = ", f12.4," seconds")') r1
            Write (ofile, '("  Total Time = ", f12.4," seconds")') r2

            Open (Unit=timef, File='times.hdf', Action='write')
            Write (timef, '("   Ando Time = ", f12.4," seconds, kpt=",i)') tm (1), 0
            Write (timef, '("   Prep Time = ", f12.4," seconds")') tm (2)
            Write (timef, '("  Embed Time = ", f12.4," seconds")') tm (3)
            Write (timef, '(i3.3," x KP Time = ", f12.4," seconds")') i, tm (4)
            Close (timef)
         Else
            Write (*,*) 'error reading times!'
            Stop
         End If


         Write (ofile, '(A)') 'Done!'
         Close (ofile)

         Call hdfwrite (h5out, '/compat', compat_flags)

         Call h5fclose_f (h5out, error)
         Call h5close_f (error)

         If (1 == 1) Then
            Call unlink ('output')
            Call unlink ('times')
            Call c_rename ('output.hdf', 10, 'output', 6)
            Call c_rename ('times.hdf', 9, 'times', 5)

            Call unlink ('tran_lr')
            Call unlink ('tran_rl')
         End If
      End If
      Call MPI_Barrier (solve_comm, ierr)
      Call mpi_finalize (ierr)

!!$ End of Program!
!!$ ---------------------------------------------------------------------------------------
      Stop
Contains
      Function rsharv (ifile, lbl, r1, r2) Result (ret)
         Implicit None
         Integer :: ifile
         Character (Len=*) :: lbl
         Real (Kind=DEF_DBL_PREC) :: r1, r2
         Logical :: ret
!!$           local
         Integer :: st
         Logical :: flag
         ret = .False.
         Rewind (ifile)
         flag = .False.
         Do while ( .Not. flag)
            Read (ifile, '(A)', IoStat=st) cwork
            If (st /= 0) Return
            flag = (index(cwork, lbl) > 0)
         End Do
         Write (*,*) trim (cwork)
         Read (ifile, '(A)') cwork
         Read (cwork(index(cwork, ':')+1:len(cwork)), '(2(3x,f15.10))') r1, r2
         ret = .True.
         Return
      End Function rsharv

      Function rtr (fname, nk, tt, rr, nmodes) Result (ret)
         Implicit None
         Character (Len=*) :: fname
         Real (Kind=DEF_DBL_PREC) :: tt (:, :), rr (:, :)
         Integer :: nmodes (:, :)
         Logical :: ret
         Integer :: nk
!!$ local
         Integer, Parameter :: ifile = 8474
         Integer :: i, ik (2)
         Real (Kind=DEF_DBL_PREC) :: k (2)
         Integer :: st
         Real (Kind=DEF_DBL_PREC) :: r1, r2
         ret = .False.
         Open (Unit=ifile, File=fname, Action='read', IoStat=st)
         If (st /= 0) Return
         Read (ifile, '(A)', IoStat=st) cwork
         If (st /= 0) Return
         Do i = 1, nk
            Read (ifile, '(2(1x,i4),2(1x,f13.8),8(1x,g17.8))', IoStat=st) ik, k, tt (1:4, i), rr (1:4, i)
            If (st /= 0) Return
         End Do

         Close (ifile)
         Do i = 1, nk
            r1 = 0.0d0
            r2 = 0.0d0
            If (tt(1, i) > 0) r1 = r1 + tt (1, i)
            If (tt(3, i) > 0) r1 = r1 + tt (3, i)
            If (rr(1, i) > 0) r1 = r1 + rr (1, i)
            If (rr(3, i) > 0) r1 = r1 + rr (3, i)

            If (tt(2, i) > 0) r2 = r2 + tt (2, i)
            If (tt(4, i) > 0) r2 = r2 + tt (4, i)
            If (rr(2, i) > 0) r2 = r2 + rr (2, i)
            If (rr(4, i) > 0) r2 = r2 + rr (4, i)
            nmodes (1, i) = Nint (r1)
            nmodes (2, i) = Nint (r2)
         End Do
         ret = .True.
      End Function rtr

      Function rcons (cons, nk, nmod) Result (ret)
         Implicit None
         Real (Kind=DEF_DBL_PREC) :: cons (:)
         Logical :: ret
         Integer :: nk, nmod (:, :)
!!$           local
         Integer :: st
         Integer, Parameter :: ifile = 47464
         Integer :: n1, n2, n3, n4, n5, ik (2)
         Real (Kind=DEF_DBL_PREC) :: k (2), r1, cn1, cn2
         Logical :: flag
         ret = .False.

         Open (Unit=ifile, File='output', Action='read', IoStat=st)
         If (st /= 0) Return

         flag = .False.
         Do while ( .Not. flag)
            Read (ifile, '(A)', IoStat=st) cwork
            If (st /= 0) Then
               Write (*,*) 'End of output file reached, wrong input!'
               Return
            End If
            flag = (index(cwork, 'IBZ, spin,  kx,      ky,   cons,nr_lm,nr_rm,ikx,iky,   weight') > 0)
         End Do

         Do i = 1, nk
            Read (ifile, '(i7,1x,i4,2(1x,f8.4),5(1x,i4),3x,f8.4,3x,g17.8)') n3, n4, k (1:2), n5, n1, n2, ik &
           & (1:2), r1, cn1
            If (n4 == 0) Then
               cons (i) = cn1
            Else
               compat_flags (cpf_no_spin_separation) = 0
               nmod (1, i) = n1
               nmod (3, i) = n2
               Read (ifile, '(i7,1x,i4,2(1x,f8.4),5(1x,i4),3x,f8.4,3x,g17.8)') n3, n4, k (1:2), n5, n1, n2, &
              & ik (1:2), r1, cn2
               nmod (2, i) = n1
               nmod (4, i) = n2
               cons (i) = .5d0 * (cn1+cn2)
            End If
         End Do
         Close (ifile)

         ret = .True.
      End Function rcons
      Function rdamp (ifile, dc) Result (ret)
         Implicit None
         Integer :: ifile
         Complex (Kind=DEF_DBL_PREC) :: dc (2, 2)
         Logical :: ret
!!$           local
         Integer :: st, i, j, k
         Logical :: flag (2)
         Real (Kind=DEF_DBL_PREC) :: dct (2, 2, 2)
         ret = .False.
         Rewind (ifile)
         dc = dcmplx (0.0d0, 0.0d0)
         flag (:) = .True.
         Do while (flag( 1) .Or. flag(2))
            Read (ifile, '(A)', IoStat=st) cwork
            If (st /= 0) Return
            If (index(cwork, 'Damping') > 0) Then
               Write (*,*) trim (cwork)
               If (index(cwork, 'Re') > 0) Then
                  i = 1
               Else
                  If (index(cwork, 'Im') > 0) Then
                     i = 2
                  Else
                     Return
                  End If
               End If
               Do j = 1, 2
                  Read (ifile, '(A)', IoStat=st) cwork
                  If (st /= 0) Return
                  k = index (cwork, ':')
                  If (k < 1) Return
                  Read (cwork(k+1:len(cwork)), '(2(2x,g17.10))', IoStat=st) dct (j, :, i)
                  If (st /= 0) Return
               End Do
               flag (i) = .False.
            End If
         End Do

         dc (:, :) = dcmplx (dct(:, :, 1), dct(:, :, 2))
         ret = .True.
      End Function rdamp

      Function rcond (ifile, lbl, cond, lnew) Result (ret)
         Implicit None
         Integer :: ifile
         Character (Len=*) :: lbl
         Real (Kind=DEF_DBL_PREC) :: cond (2, 2)
         Logical :: ret, lnew
!!$           local
         Integer :: st
         Logical :: flag
         ret = .False.
         cond (:, :) = 0.0d0
         Rewind (ifile)
         flag = .False.
         Do while ( .Not. flag)
            Read (ifile, '(A)', IoStat=st) cwork
            If (st /= 0) Return
            flag = (index(cwork, lbl) > 0)
         End Do
         Write (*,*) trim (cwork)
         lnew = .False.
         If (index(cwork, 'decomp') > 0) lnew = .True.

         Read (ifile, '(A)') cwork

         Read (ifile, '(A)') cwork
         Read (cwork(index(cwork, ':')+1:len(cwork)), '(2(3x,g17.10))') cond (1, :)
         Read (ifile, '(A)') cwork
         Read (cwork(index(cwork, ':')+1:len(cwork)), '(2(3x,g17.10))') cond (2, :)
         If (cond(1, 2) < 1e-15 .And. cond(2, 1) < 1e-15) lnew = .True.

         If ( .Not. lnew) Then
            cond (1, 1) = cond (1, 1) + cond (1, 2)
            cond (1, 2) = 0.0d0
            cond (2, 2) = cond (2, 2) + cond (2, 1)
            cond (2, 1) = 0.0d0
         End If
         ret = .True.
         Return
      End Function rcond

      Function rtms (tm, t1, t2, np) Result (ret)
         Implicit None
         Real (Kind=DEF_DBL_PREC) :: tm (:), t1, t2
         Logical :: ret
         Integer :: np
!!$           local
         Integer :: st
         Integer, Parameter :: ifile = 47464
         Integer :: i, nt
         Logical :: flag
         Real (Kind=DEF_DBL_PREC) :: tm1 (4)
         ret = .False.

         flag = .False.


         Open (Unit=ifile, File='times', Action='read', IoStat=st)
         If (st /= 0) Return
         Read (ifile, '(A)', IoStat=st) cwork
         If (st /= 0) Return
         Read (cwork(index(cwork, '=')+2:len(cwork)), '(f12.4)') t1

         flag = .True.
         tm1 (:) = 0.0d0
         tm (:) = 0.0d0
         nt = 0
         np = 0
         Do while (flag)
            i = 0
            Read (ifile, '(A)', IoStat=st) cwork
            If (st /= 0) Return
            If (index(cwork, 'Ando') > 0) i = 1
            If (index(cwork, 'Prep') > 0) i = 2
            If (index(cwork, 'Embed') > 0) i = 3
            If (index(cwork, 'KP') > 0) i = 4

            flag = .Not. (index(cwork, 'Total') > 0)
            If (i > 0 .And. flag .And. st == 0) Then
               Read (cwork(index(cwork, '=')+2:len(cwork)), '(f12.4)', IoStat=st) tm1 (i)
               If (st /= 0) Return
            End If

            If (i == 4) Then
               nt = nt + 1
               tm (:) = tm (:) + tm1 (:)
               tm1 (:) = 0.0d0
               If (np == 0) Then
                  If (index(cwork(1:3), '*') > 0) Then
                     np = 999
                  Else
                     Read (cwork, '(i3.3)') np
                  End If
               End If
            End If
         End Do
         tm = tm / dble (nt)

         t2 = 0.0d0
         If ( .Not. flag) Read (cwork(index(cwork, '=')+2:len(cwork)), '(f12.4)') t2

         ret = .True.
         Return
      End Function rtms

      Function get_i (txt, optname, val) Result (res)
         Implicit None
         Character (Len=*) :: txt
         Character (Len=*) :: optname
         Integer :: val
         Logical :: res
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Integer :: p1, l2, p2
         Character (Len=256) :: cw
         res = .False.
         topt = optname
         Call s_cap (topt)
         p1 = index (txt, topt)
         If (p1 < 1) Return
         l2 = len_trim (txt)
         p1 = p1 + len_trim (optname)
         p2 = index (txt(p1:l2), '=')
         If (p1 < 1) Return
         cw = trim (adjustl(txt(p1+p2:l2)))
         Read (cw(1:l2),*) val
         res = .True.
      End Function get_i

      Function get_d (txt, optname, val) Result (res)
         Implicit None
         Character (Len=*) :: txt
         Character (Len=*) :: optname
         Real (Kind=DEF_DBL_PREC) :: val
         Logical :: res
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Integer :: p1, l2, p2
         Character (Len=256) :: cw
         res = .False.
         topt = optname
         Call s_cap (topt)
         p1 = index (txt, topt)
         If (p1 < 1) Return
         l2 = len_trim (txt)
         p1 = p1 + len_trim (optname)
         p2 = index (txt(p1:l2), '=')
         If (p1 < 1) Return
         cw = trim (adjustl(txt(p1+p2:l2)))
         Read (cw(1:l2),*) val
         res = .True.
      End Function get_d

End Program ascii2hdf
