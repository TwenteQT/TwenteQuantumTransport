!!$ $Id: readcfg.F90 2008 2014-04-01 14:51:28Z rien $
#include "math_def.h"
Module readcfg
      Use logging
      Use sparse_solvers
      !Use petsc_solvers
      Use hamiltonian
      Use potpars
      Use iter_solver
      Use atoms_module
      Use ando_module
      Use bzgrid
      Use adaptive_bzgrid
      Implicit None


      Type t_config_option
         Character (Len=100) :: name
         Character (Len=200) :: value
         Integer :: vl, nl
         Type (t_config_option), Pointer :: next
      End Type t_config_option

      Type t_damp_opts
         Real (Kind=DEF_DBL_PREC) :: dstep
         Integer :: axesdir
      End Type t_damp_opts

      Type t_ddw_opts
         Real (Kind=DEF_DBL_PREC) :: dr,dphi
         Integer :: nmodes
      End Type t_ddw_opts

      Type t_options
         Real (Kind=DEF_DBL_PREC) :: Eoffset, drshift, spop_angle
         Integer :: par_nprock
         Integer :: onlyspin
         Integer :: loglvl, lowmem
         Integer :: sdcdir, signat
         Integer :: leads_ident, write_tr, icont
         Integer :: needlr, needrl, writewf, writeham, writetm, writerhs
         Integer :: writetorque, writedens, writeveloc
         Integer :: do_iatcu, do_oldcurr
         Integer :: do_OAMdens !orbital angular momentum density
         Integer :: do_iatLcu
         Integer :: alltimes
         Integer :: embed_only_left
         Type (t_leq_param) :: leq
         !Type (t_petsc_opt) :: petsc_opt
         Type (t_pot_opts) :: po
         Type (t_ando_ctl) :: actl
         Type (t_potparp_options) :: atopt
         Type (t_bz_opts) :: bzo
         Type (t_damp_opts) :: damp
         Type (t_ddw_opts) :: ddw
      End Type t_options

      Type t_bsfp_options
         Real (Kind=DEF_DBL_PREC) :: Eoffset
         Integer :: nkz, isdeconv, symsc, isav,isdrop
         Integer :: loglvl
         Integer :: signat
         Type (t_pot_opts) :: po
         Type (t_potparp_options) :: atopt
         Type (t_bz_opts) :: bzo
      End Type t_bsfp_options

      Type t_bandopts
         Integer :: loglvl, nE
         Real (Kind=DEF_DBL_PREC) :: bE, eE, fermiE, sk_alat = 0d0 !sk_alat added for sk
         Real (Kind=DEF_DBL_PREC) :: kpar (2), ucrit
         Integer :: sdcdir
         Integer :: par_nprock
         Type (t_pot_opts) :: po
         Type (t_leq_param) :: leq
         Type (t_ando_ctl) :: actl
         Type (t_potparp_options) :: atopt
         logical :: sk_ortho = .TRUE. !added for sk
         Integer :: gensol       !added for sk
      End Type t_bandopts

      Type t_fsopts !options for Fermi Surface program
         Integer :: loglvl
         Real (Kind=DEF_DBL_PREC) :: en, ucrit
         Type (t_adapt_bz_opts) :: bzo
         Integer :: sdcdir
         Integer :: par_nprock
         Type (t_pot_opts) :: po
         Type (t_leq_param) :: leq
         Type (t_ando_ctl) :: actl
         Type (t_potparp_options) :: atopt
      End Type t_fsopts

Contains

      Subroutine free_option_list (optlist)
         Implicit None
         Type (t_config_option), Pointer :: optlist
!!$ Local vars
         Type (t_config_option), Pointer :: ptr, ptr1
         ptr => optlist
         Do while (associated(ptr))
            ptr1 => ptr
            ptr => ptr%next
            Deallocate (ptr1)
            Nullify (ptr1)
         End Do
         Nullify (optlist, ptr)
      End Subroutine free_option_list


      Subroutine solver_config_init (optlist, opt, ispar, np)
         Implicit None
         Type (t_leq_param), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist
         Logical, Optional :: ispar
         Integer :: np
!!$ Local vars
         Integer :: mpar (4)
         mpar = (/ - 1, 2, 0, 7 /)
         If ((np /= 1) .And. ispar .Eqv. .True.) mpar = (/ 100, 4, 0, 6 /)

         Call do_log (1, '### Options to control solvers behaviour:')
         Call do_log (1, '### ')

         Call do_log (1, '### Direct linear solver options:')
         opt%solver = getopt_i (optlist, 'lineq_solver', 0)

         opt%memover = getopt_i (optlist, 'mumps_memory_overhead', mpar(1))
         opt%partitioner = getopt_i (optlist, 'mumps_partitioner', mpar(2))
         opt%redistmtx = getopt_i (optlist, 'mumps_redist_mtx', mpar(3))
         opt%permscale = getopt_i (optlist, 'mumps_permscale', mpar(4))

         Call do_log (1, '### if par_nprock/=1, then workhost controls occupancy of root-process in solver co&
        &mm:')
         opt%workhost = getopt_i (optlist, 'mumps_workhost', 1)
         opt%debug = getopt_i (optlist, 'mumps_debug',-1)
         Call do_log (1, '# For OOC solver set variable MUMPS_TMP_DIR. Then OOC solver will be selected')
         opt%tmpdir = getopt_s (optlist, 'MUMPS_TMP_DIR', '')
         Call do_log (1, '# For half-iterative solver (LINEQ_SOLVER==1)')
         If (opt%solver == 1) Then
            opt%maxit = getopt_i (optlist, 'lineq_maxiter', 4)
            opt%tol = getopt_d (optlist, 'lineq_accuracy', 1.0d-12)
         Else
            Call log_opt_i ('lineq_maxiter', 4)
            Call log_opt_d ('lineq_accuracy', 1.0d-12)
         End If


#if defined(USE_SCSL) || defined(USE_WSMP) || defined(USE_PARDISO)
         opt%par_nprock = 1
#endif
      End Subroutine solver_config_init


      Subroutine atopt_config_init (optlist, opt)
         Implicit None
         Type (t_potparp_options), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist

         Call do_log (1, '### Options to control calculation of atomic parameters:')
         Call do_log (1, '### ')
         Call do_log (1, '### Scalar relativistic calculations')
         opt%irel = getopt_i (optlist, 'atopt_irel', 1)
         opt%nsirk = getopt_i (optlist, 'atopt_nsirk', 2)
         opt%dE = getopt_d (optlist, 'atopt_dE_step', 0.0001d0)
!!$  call do_log(1,'######################################')
         Call do_log (1, '')

      End Subroutine atopt_config_init

      Subroutine bz_config_init (optlist, opt)
         Implicit None
         Type (t_bz_opts), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist
!!$ Local
         Real (Kind=8) :: t1 (4) = 0.0d0
         Integer :: needsym, needar
         Call do_log (1, '### Parameters of BZ sampling:')
         Call do_log (1, '#')

         needsym = 1
         needar = 1

         Call do_log (1, '# Perform calculation only for one K-point')
         Call do_log (1, '# Ignore rest of BZ options')
         If (is_opt_exists(optlist, 'bz_kpoint') == 1) Then
            opt%kpoint = getopt_da (optlist, 'bz_kpoint', 2, t1)
            opt%nsym = - 10
            opt%nk = 1
            needsym = 0
            needar = 0
         Else
            Call do_log (1, '# BZ_KPOINT = [0.0,0.0]')
         End If
         Call do_log (1, '#')

         Call do_log (1, '### Perform calculation for user specified area of BZ')
         Call do_log (1, '### Ignore symmetry related BZ options')
         Call do_log (1, '### Performs sampling only in selected BZ_AREA instead of full B.Zone')
         If (needar == 1 .And. is_opt_exists(optlist, 'bz_area') == 1) Then
            opt%bzone = getopt_da (optlist, 'bz_area', 4, t1)
            needsym = 0
            If (is_opt_exists(optlist, 'bz_area') == 1) Then
               needsym = getopt_i (optlist, 'bz_area_sym', 0)
            End If
            opt%nsym = - 1 - needsym
            needsym = 0

         Else
            Call do_log (1, '# BZ_AREA = [0.0,1.0,0.0,1.0] # [x0,x1,y0,y1]')
         End If
         Call do_log (1, '#')
         Call do_log (1, '### Defines number of K-points ib B.Zone for sampling')
         Call do_log (1, '### NK=(2*bz_num_kpoints)^2')
         If (needar == 1) Then
            opt%nk = getopt_i (optlist, 'bz_num_kpoints', 30)
         Else
            Call log_opt_i ('bz_num_kpoints', 1)
         End If
         Call do_log (1, '#')

         Call do_log (1, '### bz_sym and bz_inve defines symmetry in B.Zone')
         If (needsym == 1) Then
            opt%nsym = getopt_i (optlist, 'bz_sym', 0)
            opt%inve = getopt_i (optlist, 'bz_inve', 0)
            Call do_log (1, '### Set if deal with system with time-reversal symmetry')
            opt%iftr = getopt_i (optlist, 'bz_timereversal', 0)
         Else
            Call log_opt_i ('bz_sym', 0)
            Call log_opt_i ('bz_inve', 0)
            Call do_log (1, '### Set if deal with system with time-reversal symmetry')
            Call log_opt_i ('bz_timereversal', 0)
         End If

         Call do_log (1, '### If "bz_part" not equal to 1.0d0, then only bz_part*100%')
         Call do_log (1, '### in center of BZ will be used for sampling (bz_part>0)')

         If (needsym == 1) Then
            opt%BZpart = getopt_d (optlist, 'bz_part', 1.0d0)
         Else
            Call log_opt_d ('bz_part', 1.0d0)
         End If
         Call do_log (1, '###')
!!$  call do_log(1,'######################################')
         Call do_log (1, '')

      End Subroutine bz_config_init

      Subroutine adapt_bz_config_init(optlist, opt)
        Implicit None
         Type (t_adapt_bz_opts), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist

         Call do_log (1, '### Parameters of BZ sampling:')
         Call do_log (1, '#')

         Call do_log (1, '### Defines number of K-points ib B.Zone for sampling')
         Call do_log (1, '### NK=(2*bz_num_kpoints)^2')
         opt%nk = getopt_i (optlist, 'bz_num_kpoints', 30)

         Call do_log (1, '### Defines the tolerance on the refinement of the B.Zone grid')
         opt%tol = getopt_d (optlist, 'bz_tol', 1d-3)

         !! More options based on symmetry etc.
          
      End Subroutine adapt_bz_config_init 

      Subroutine ando_config_init (optlist, opt, so, nh, ap_in)
         Implicit None
         Type (t_ando_ctl), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist
         Real (Kind=DEF_DBL_PREC) :: ap_in
         Integer :: so, nh
!!$ Local vars
         Real (Kind=DEF_DBL_PREC) :: ap
         Call do_log (1, '### Options to control Ando solver:')
         Call do_log (1, '### ')

         Call do_log (1, '### Chooses eigenvalue solver used by Ando solver.')
         Call do_log (1, '### 0 - simple(fastest), 10 - general(optimal choice), 17 - quadro(slow, most stabl&
        &e)')
         Call do_log (1, '### in some cases when accuracy not enough, code will automaticaly fail-back to 17')
         opt%gensol = getopt_i (optlist, 'eigval_solver', 10)
         If (so /= 0 .And. opt%gensol < 15) Then
            opt%gensol = 15
            Call do_log (1, 'EIGVAL_SOLVER = 15 # required by spin-orbit')
         End If

         opt%refine = getopt_i (optlist, 'eigval_refine', 1)

         ap = 1.0d0 / nh
         If (ap_in >= 0.0d0) ap = ap_in
         Call do_log (1, '### in some cases')
         Call do_log (1, '### use ando_part_of_spectrum*100% of eigenvalues for Ando')

         opt%specpart = getopt_d (optlist, 'ando_part_of_spectrum', ap)
         Call do_log (1, '### Leads parameters:')
         Call do_log (1, '### ')
         Call do_log (1, '### Polarization in left/right leads is now automatically analyzed ')
         Call do_log (1, '### left_pol and right_pol are not used anymore')
         Call do_log (1, '### ')
         Call do_log (1, '### Identity of leads is to be automatically decided')
         Call do_log (1, '### leads_ident is not used anymore')
         opt%extlead = getopt_i (optlist, 'ext_leads', 0)
!!$  call do_log(1,'######################################')
         Call do_log (1, '')

         Call do_log (1, '### Write data for channel decomposition:')
         opt%writedecomp = getopt_i (optlist, 'ando_decomp_data', 0)
!!$ Call do_log(1,'######################################')
         Call do_log (1, '')
      End Subroutine ando_config_init

      Subroutine pot_config_init (optlist, opt)
         Implicit None
         Type (t_pot_opts), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist
!!$ Local vars
         Real (Kind=8) :: t1 (4) = 0.0d0
         Logical :: flag
         Character (Len=300) :: cwork

         Call do_log (1, '### Kind of system to solve: 1=P-S, 2=H, 3=H-HOH, 5=Model Hamiltonian, 6=EMTO')
         opt%kind = getopt_i (optlist, 'equation_type', 1)
         Call do_log (1, '### ')
         Call do_log (1, '### Spin-orbit interaction parameters:')
         Call do_log (1, '### ')
         Call do_log (1, '### Switch on spin-orbit. Valid only if "trans_full" binary used')
         If (opt%kind == 3) opt%nhops = 2

         opt%so = getopt_i (optlist, 'SO_ON', 0)
         Call do_log (1, '### Scaling of spin-orbit')
         opt%so_scale = getopt_d (optlist, 'SO_scale', 1.0d0)
         Call do_log (1, '### Complex part of energy')
         opt%ei = getopt_d (optlist, 'EI', 0.0d0)

         Call do_log (1, '######################################')

         Call do_log (1, '### Rotation of magnetization parameters:')
         Call do_log (1, '### ')
         Call do_log (1, '### Valid only if "trans_full" binary used')
         opt%rot_mask = getopt_i (optlist, 'rot_mask', 0)
         If (opt%rot_mask == 1) Then
            Call do_log (1, '### Using spin-rotation mask from external file "rot_mask"!')
            Call do_log (1, '### For left/right lead rotations will be')
            Call do_log (1, '### used same values as for first/last atoms!')
         End If
         Call do_log (1, '### Angle for stiff rotation of magnetization (after rot_mask is applied)')
         opt%rot_angle = getopt_da (optlist, 'rot_angle', 3, t1, flag)
         If ( .Not. flag) Then
            Call do_log (1, '### Could not read rot_angle, will try to read old-style 2-element vector')
            opt%rot_angle = getopt_da (optlist, 'rot_angle', 2, t1)
            opt%rot_angle (3) = 0.0d0
            Write (cwork, '(A,A,2(G14.6,","),G14.6,A1)') 'ROT_ANGLE', ' = [', opt%rot_angle(1:3), ']'
            Call do_log (1, cwork)
         End If
         opt%rot_angle = opt%rot_angle * DEF_M_PI
         opt%rotcf = getopt_i (optlist, 'rotcf', 0)
         Call do_log (1, '### Values in rotated potential parameters')
         Call do_log (1, '###  smaller than "rnegl" will be neglected.')
         opt%rnegl = getopt_d (optlist, 'rnegl', 1.0d-20)
         Call do_log (1, '')

         If (opt%kind == 5) Then
            Call do_log (1, '')
            Call do_log (1, '### Parameters for model Hamiltonian')
            opt%ngridx = getopt_i (optlist, 'ngridx', 0)
            Call do_log (1, '### Number of x-grid in scattering region:')
            opt%ngridy = getopt_i (optlist, 'ngridy', 0)
            Call do_log (1, '### Number of y-grid in scattering region:')
            opt%ngridz = getopt_i (optlist, 'ngridz', 0)
            Call do_log (1, '### Number of z-grid in scattering region:')
            opt%ngrid = opt%ngridx * opt%ngridy * opt%ngridz

            opt%nlx = getopt_i (optlist, 'nlx', 0)
            Call do_log (1, '### Number of x-grid in the left lead:')
            opt%nly = getopt_i (optlist, 'nly', 0)
            Call do_log (1, '### Number of y-grid in the left lead:')
            opt%nlz = getopt_i (optlist, 'nlz', 0)
            Call do_log (1, '### Number of z-grid in the left lead:')
            opt%nrx = getopt_i (optlist, 'nrx', 0)
            Call do_log (1, '### Number of x-grid in the right lead:')
            opt%nry = getopt_i (optlist, 'nry', 0)
            Call do_log (1, '### Number of y-grid in the right lead:')
            opt%nrz = getopt_i (optlist, 'nrz', 0)
            Call do_log (1, '### Number of z-grid in the right lead:')

            opt%periodic = getopt_i (optlist, 'periodic', 0)
            Call do_log (1, '### Periodic boundary condition in model Hamiltonian system')

            opt%dx = getopt_d (optlist, 'lengthx', 0.d0)
            Call do_log (1, '### Grid length in x direction')
            opt%dy = getopt_d (optlist, 'lengthy', 0.d0)
            Call do_log (1, '### Grid length in y direction')
            opt%dz = getopt_d (optlist, 'lengthz', 0.d0)
            Call do_log (1, '### Grid length in z direction')

            opt%energy = getopt_d (optlist, 'energy', 0.d0)
            Call do_log (1, '### Fermi energy in the model Hamiltonian system')
            opt%gamma1 = getopt_d (optlist, 'gamma1', 1.d0)
            Call do_log (1, '### Parameter in the model Hamiltonian system')
            opt%gamma2 = getopt_d (optlist, 'gamma2', 0.d0)
            Call do_log (1, '### Parameter in the model Hamiltonian system')
            opt%h0 = getopt_d (optlist, 'h0', 0.d0)
            Call do_log (1, '### Local exchange strength in the model Hamiltonian system')
            opt%v0 = getopt_d (optlist, 'v0', 0.d0)
            Call do_log (1, '### Impurity scattering strength in the model Hamiltonian system')
         End If

      End Subroutine pot_config_init

      Subroutine damp_config_init (optlist, opt)
         Type (t_damp_opts), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist
         Call do_log (1, '### Damping options:')
         opt%dstep = getopt_d (optlist, 'damping_diff_step', 0.0001d0)
         opt%dstep = opt%dstep * DEF_M_PI
         opt%axesdir = getopt_i (optlist, 'damping_axis_orientation', 0)
         Call do_log (1, '### If damping_axis_orientation==0 then differentiation about axis defined by rot_a&
        &ngle')
         Call do_log (1, '### If damping_axis_orientation==1 then differentiation about axis defined by macro&
        &-spin orientation')
         Call do_log (1, '')
      End Subroutine damp_config_init

      Subroutine ddw_config_init (optlist, opt)
         Type (t_ddw_opts), Intent (Inout) :: opt
         Type (t_config_option), Pointer, Intent (In) :: optlist
         Call do_log (1, '### Domain-wall options:')
         opt%dr = getopt_d (optlist, 'dw_center_step', 0.01d0)
         opt%dphi = getopt_d (optlist, 'dw_flapping_step', 0.0001d0)
         opt%dphi = opt%dphi * DEF_M_PI
         opt%nmodes = getopt_i (optlist, 'dw_nmodes', 0)
         Call do_log (1, '### If dw==0, calculate both alpha_theta (only if dw==1) &
        & and alpha_phi (only if dw==2)')
         Call do_log (1, '')
      End Subroutine ddw_config_init

      Subroutine trans_config_init (opt, ap_in, ispar_in)
         Implicit None
         Type (t_options) :: opt
         Real (Kind=8), Optional :: ap_in
         Logical, Optional :: ispar_in
!!$ Local vars
         Type (t_config_option), Pointer :: optlist
         Real (Kind=8) :: ap
         Logical :: ispar
         Integer :: dummy
         Logical :: def
         
         Call read_config_file ('trans.conf', optlist)

         Call do_log (1, '######################################')
         Call do_log (1, '#     Equivalent of given config     #')
         Call do_log (1, '######################################')
         Call do_log (1, '')

         Call do_log (1, '### Output verbosity level:')
         opt%loglvl = getopt_i (optlist, 'LogLevel', 1)
!!$  call do_log(1,'######################################')
         opt%alltimes = getopt_i (optlist, 'write_all_times', 0)
         Call do_log (1, '')

         Call pot_config_init (optlist, opt%po)
         Call atopt_config_init (optlist, opt%atopt)
         ap = - 1.0d0
         If (present(ap_in)) ap = ap_in
         Call ando_config_init (optlist, opt%actl, opt%po%so, opt%po%nhops, ap)
         Call bz_config_init (optlist, opt%bzo)

         Call do_log (1, '### when performs spin decomposition, if')
         Call do_log (1, '### "spin_dec_dir"= 0 , then Sz=Lz, else spin states')
         Call do_log (1, '### correspond to Sz direction in lead')
         opt%sdcdir = getopt_i (optlist, 'spin_dec_dir', 0)

         Call do_log (1, '### for rotated magnetic lead "spin_dec_dir"= 0 will lie.')
         Call do_log (1, '### for rotated non=magnetic lead "spin_dec_dir"= 1 will lie.')

         Call do_log (1, '######################################')
         Call do_log (1, '')

         Call do_log (1, '### General options:')
         Call do_log (1, '### ')
         Call do_log (1, '### Shift from Fermy level')
         opt%Eoffset = getopt_d (optlist, 'energy_offset', 0.0d0)

         Call do_log (1, '### Minimal number of MPI processes per k-point:')
         Call do_log (1, '### if par_nprock == 0, then effectively par_nprock == mpi_comm_size (mpi_comm_worl&
        &d)')
         opt%par_nprock = getopt_i (optlist, 'par_nprock', 1)

         Call do_log (1, '### Option to only embed the left lead, not the right')
         Call do_log (1, '### (only relevant to the htrans_eigvals program)')
         opt%embed_only_left = getopt_i (optlist, 'embed_only_left', 0)

         ispar = .True.
         If (present(ispar_in)) ispar = ispar_in
         Call solver_config_init (optlist, opt%leq, ispar, opt%par_nprock)


!!$  call do_log(1,'######################################')
         Call do_log (1, '')

         Call do_log (1, '### Try to reduce memory usage on a price of performance')
         opt%lowmem = getopt_i (optlist, 'conservative_memory_usage', 0)

!!$  call do_log(1,'######################################')
         Call do_log (1, '')


         Call do_log (1, '### Work to be performed:')
         Call do_log (1, '### Preform calculations only for selected spin')
         Call do_log (1, '### 0 - both, 1- Up, 2-Down')
         Call do_log (1, '### valid only for "trans_ps" binary')
         opt%onlyspin = getopt_i (optlist, 'do_only_spin', 0)
         If (opt%onlyspin > 2) opt%onlyspin = 2
         If (opt%onlyspin < 0) opt%onlyspin = 0

         Call do_log(1, '### Perform downfolded calculations')
         opt%po%df = getopt_i (optlist, 'downfold', 0)

         Call do_log(1, '### Use non-spherical corrections in the Hamiltonian')
         opt%po%nonsph = getopt_i (optlist, 'nonsph', 0) 

         Call do_log (1, '### Calculations for left-to-right going electrons:')
         opt%needlr = getopt_i (optlist, 'do_Left_to_Right', 1)
         Call do_log (1, '### Calculations for right-to-left going electrons:')
         opt%needrl = getopt_i (optlist, 'do_Right_to_Left', 1)
         If (opt%needlr == 0 .And. opt%needrl == 0) Then
            Call do_log (1, '# Give me a job! :)')
            Stop
         End If
         Call do_log (1, '### Specular transmission is calculated by default now')
         Call do_log (1, '### do_spec is not used anymore')
!!$  call do_log(1,'######################################')
         Call do_log (1, '')

         Call damp_config_init (optlist, opt%damp)
         Call ddw_config_init (optlist, opt%ddw)

         Call do_log (1, '### Write extra data:')
         Call do_log (1, '### Write wave functions:')
         opt%writewf = getopt_i (optlist, 'write_wavefuncs', 0)
         Call do_log (1, '### Write hamiltionians (or P-S) used for calculations:')
         opt%writeham = getopt_i (optlist, 'write_ham', 0)
         Call do_log (1, '### Write right-hand-side of linear equations')
         opt%writerhs = getopt_i (optlist, 'write_rhs', 0)
         Call do_log (1, '### Write transmission matrices')
         opt%writetm = getopt_i (optlist, 'write_trans_matrix', 0)
         Call do_log (1, '### Write torques')
         opt%writetorque = getopt_i (optlist, 'write_torque', 0)
         Call do_log (1, '### Write lead velocities')
         opt%writeveloc = getopt_i (optlist, 'write_veloc', 0)
         Call do_log (1, '### Do old currents (commutator [H, z])')
         opt%do_oldcurr = getopt_i (optlist, 'do_oldcurr', 0, def)
         If (def) Then ! this block can be deleted if compatibility with write_curr is not needed anymore
            Call do_log (1, '### Write currents (old option kept for compatibility if do_oldcurr undefined)')
            dummy = getopt_i (optlist, 'write_current', 0, def)
            If (def .eqv. .FALSE.) opt%do_oldcurr = dummy + 1
         End If
         Call do_log (1, '### Calculate interatomic currents (only for eq type 2, 3 or 6)')
         opt%do_iatcu = getopt_i (optlist, 'do_iatcu', 0)
         If (opt%do_iatcu /= 0 .and. all(opt%po%kind .ne. (/2, 3, 6/))) Then
            Call do_log (1, '# Option do_iatcu incompatible with equation_type')
            stop
         End If

         Call do_log (1, '### Calculate interatomic Orbital currents (only for eq type 2, 3 or 6)')
         opt%do_iatLcu = getopt_i (optlist, 'do_iatlcu', 0)
         If (opt%do_iatcu /= 0 .and. all(opt%po%kind .ne. (/2, 3, 6/))) Then
            Call do_log (1, '# Option do_iatLcu incompatible with equation_type')
            stop
         End If

         Call do_log (1, '### Calculate interatomic Orbital density (only for eq type 2, 3 or 6)')
         opt%do_OAMdens = getopt_i (optlist, 'do_OAMdens', 0)
         If (opt%do_OAMdens /= 0 .and. all(opt%po%kind .ne. (/2, 3, 6/))) Then
            Call do_log (1, '# Option do_OAMdens incompatible with equation_type')
            stop
         End If


         opt%writedens = getopt_i (optlist, 'write_dens', 0)

         Call do_log (1, '### Continue interrupted calculation')
         opt%icont = getopt_i (optlist, 'icont', 0)

         !Call do_log (1, '### PETSc linear solver options:')
         !Call do_log (1, '### only used if LINEQ_SOLVER >= 2')
         !opt%petsc_opt%solver = getopt_i (optlist, 'ksp_solver', 1)
         !opt%petsc_opt%preco  = getopt_i (optlist, 'ksp_preco', 1)
         !opt%petsc_opt%tol  = getopt_d (optlist, 'ksp_tol', 1d-12)
      


         Call do_log (1, '### Movement of domain wall center')
         opt%drshift = getopt_d (optlist, 'drshift', 0.d0)

         Call do_log (1, '')
         Call do_log (1, '######################################')
         Call do_log (1, '#     End of config file is here!    #')
         Call do_log (1, '######################################')
         Call free_option_list (optlist)
         opt%signat = opt%po%kind * 100 + opt%po%so

      End Subroutine trans_config_init

      Subroutine bsfp_config_init (opt, ap_in, ispar_in)
         Implicit None
         Type (t_bsfp_options) :: opt
         Real (Kind=8), Optional :: ap_in
         Logical, Optional :: ispar_in
!!$ Local vars
         Type (t_config_option), Pointer :: optlist
         Real (Kind=8) :: ap
         Logical :: ispar
         Call read_config_file ('trans.conf', optlist)

         Call do_log (1, '######################################')
         Call do_log (1, '#     Equivalent of given config     #')
         Call do_log (1, '######################################')
         Call do_log (1, '')

         Call do_log (1, '### Output verbosity level:')
         opt%loglvl = getopt_i (optlist, 'LogLevel', 1)
!!$  call do_log(1,'######################################')
         Call do_log (1, '')

         Call pot_config_init (optlist, opt%po)
         Call atopt_config_init (optlist, opt%atopt)
         ap = - 1.0d0
         If (present(ap_in)) ap = ap_in
         Call bz_config_init (optlist, opt%bzo)


         opt%nkz = getopt_i (optlist, 'NKZ', 100)
         opt%isdrop = getopt_i (optlist, 'write_zdep', 0)
         opt%isdeconv = getopt_i (optlist, 'DO_DECONV', 1)
         opt%symsc = getopt_i (optlist, 'DO_SC_SYM', 0)
         opt%isav = getopt_i (optlist, 'AVERAGE_ZCPF', 0)

         Call do_log (1, '')
         Call do_log (1, '######################################')
         Call do_log (1, '#     End of config file is here!    #')
         Call do_log (1, '######################################')
         Call free_option_list (optlist)
         opt%signat = opt%po%kind * 100 + opt%po%so

      End Subroutine bsfp_config_init

      Subroutine band_config_init (opt, ap_in, ispar_in)
         Implicit None
         Type (t_bandopts) :: opt
         Logical, Optional :: ispar_in
         Real (Kind=DEF_DBL_PREC), Optional :: ap_in
!!$ Local vars
         Type (t_config_option), Pointer :: optlist
         Real (Kind=8) :: t1 (2) = (/ 0.0d0, 0.0d0 /), ap
         Logical :: ispar

         Call read_config_file ('band.conf', optlist)
         Call do_log (1, '######################################')
         Call do_log (1, '#     Equivalent of given config     #')
         Call do_log (1, '######################################')
         Call do_log (1, '')

         Call do_log (1, '### Output verbosity level:')
         opt%loglvl = getopt_i (optlist, 'LogLevel', 1)
         Call do_log (1, '')
         opt%fermiE = getopt_d (optlist, 'Fermi_energy',0.0d0)
         opt%bE = getopt_d (optlist, 'energy_from',-0.2d0)
         opt%eE = getopt_d (optlist, 'energy_to', 0.2d0)

         opt%ucrit = getopt_d (optlist, 'ucrit', 1.0d-6)

         opt%nE = getopt_i (optlist, 'num_points', 100)


         opt%sk_alat = getopt_d(optlist, 'sk_alat', 1.0d0)
         opt%sk_ortho = (getopt_i(optlist, 'sk_ortho', 1) == 1)

         Call pot_config_init (optlist, opt%po)
         Call do_log(1, '### Use downfolding (see atomlist_EMTO file')
         opt%po%df = getopt_i (optlist, 'downfold', 0)

         Call do_log(1, '### Use non-spherical corrections in the Hamiltonian')
         opt%po%nonsph = getopt_i (optlist, 'nonsph', 0) 

         Call atopt_config_init (optlist, opt%atopt)
         ap = - 1.0d0
         If (present(ap_in)) ap = ap_in
         Call ando_config_init (optlist, opt%actl, opt%po%so, opt%po%nhops, ap)

         Call do_log (1, '### when performs spin decomposition, if')
         Call do_log (1, '### "spin_dec_dir"= 0 , then Sz=Lz, else spin states')
         Call do_log (1, '### correspond to Sz direction in lead')
         opt%sdcdir = getopt_i (optlist, 'spin_dec_dir', 0)

         Call do_log (1, '### for rotated magnetic lead "spin_dec_dir"= 0 will lie.')
         Call do_log (1, '### for rotated non=magnetic lead "spin_dec_dir"= 1 will lie.')

         Call do_log (1, '######################################')
         Call do_log (1, '')
         opt%kpar = getopt_da (optlist, 'kpar', 2, t1)

         Call do_log (1, '### Minimal number of MPI processes per k-point:')
         Call do_log (1, '### if par_nprock == 0, then effectively par_nprock == mpi_comm_size (mpi_comm_worl&
        &d)')
         opt%par_nprock = getopt_i (optlist, 'par_nprock', 1)

         ispar = .True.
         If (present(ispar_in)) ispar = ispar_in
         Call solver_config_init (optlist, opt%leq, ispar, opt%par_nprock)


         Call do_log (1, '######################################')
         Call free_option_list (optlist)
      End Subroutine band_config_init

      Subroutine fs_config_init (opt, ap_in, ispar_in)
         Implicit None
         Type (t_fsopts) :: opt
         Logical, Optional :: ispar_in
         Real (Kind=DEF_DBL_PREC), Optional :: ap_in
!!$ Local vars
         Type (t_config_option), Pointer :: optlist
         Real (Kind=8) :: t1 (2) = (/ 0.0d0, 0.0d0 /), ap
         Logical :: ispar

         Call read_config_file ('fermi.conf', optlist)
         Call do_log (1, '######################################')
         Call do_log (1, '#     Equivalent of given config     #')
         Call do_log (1, '######################################')
         Call do_log (1, '')

         Call do_log (1, '### Output verbosity level:')
         opt%loglvl = getopt_i (optlist, 'LogLevel', 1)
         Call do_log (1, '')

         opt%en = getopt_d (optlist, 'energy',-0.2d0)
         opt%ucrit = getopt_d (optlist, 'ucrit', 1.0d-6)

         Call adapt_bz_config_init (optlist, opt%bzo)

         Call pot_config_init (optlist, opt%po)
         Call do_log(1, '### Use downfolding (see atomlist_EMTO file')
         opt%po%df = getopt_i (optlist, 'downfold', 0)

         Call do_log(1, '### Use non-spherical corrections in the Hamiltonian')
         opt%po%nonsph = getopt_i (optlist, 'nonsph', 0) 

         Call atopt_config_init (optlist, opt%atopt)
         ap = - 1.0d0
         If (present(ap_in)) ap = ap_in
         Call ando_config_init (optlist, opt%actl, opt%po%so, opt%po%nhops, ap)

         Call do_log (1, '### when performs spin decomposition, if')
         Call do_log (1, '### "spin_dec_dir"= 0 , then Sz=Lz, else spin states')
         Call do_log (1, '### correspond to Sz direction in lead')
         opt%sdcdir = getopt_i (optlist, 'spin_dec_dir', 0)

         Call do_log (1, '### for rotated magnetic lead "spin_dec_dir"= 0 will lie.')
         Call do_log (1, '### for rotated non=magnetic lead "spin_dec_dir"= 1 will lie.')

         Call do_log (1, '######################################')

         Call do_log (1, '### Minimal number of MPI processes per k-point:')
         Call do_log (1, '### if par_nprock == 0, then effectively par_nprock == mpi_comm_size (mpi_comm_worl&
        &d)')
         opt%par_nprock = getopt_i (optlist, 'par_nprock', 1)

         ispar = .True.
         If (present(ispar_in)) ispar = ispar_in
         Call solver_config_init (optlist, opt%leq, ispar, opt%par_nprock)


         Call do_log (1, '######################################')
         Call free_option_list (optlist)
      End Subroutine fs_config_init

!!$  subroutine do_log(wantlevel,text)
!!$    implicit none
!!$    integer::wantlevel
!!$    character(len=*)::text
!!$    write(*,'(A)') text
!!$  end subroutine do_log


      Function is_opt_exists (optlist, optname) Result (res)
         Implicit None
         Type (t_config_option), Pointer :: optlist
         Character (Len=*) :: optname
         Integer :: res
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Type (t_config_option), Pointer :: ptr
         topt = optname
         ptr => optlist
         Call s_cap (topt)
         res = 0
         Do while (associated(ptr))
            If (ptr%name(1:ptr%nl) == topt) Then
               res = 1
               Exit
            End If
            ptr => ptr%next
         End Do
      End Function is_opt_exists

      Subroutine log_opt_i (optname, val)
         Implicit None
         Character (Len=*) :: optname
         Integer :: val
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Character (Len=300) :: slog
         topt = optname
         Call s_cap (topt)

         Write (slog, '(A,A,A,I6)') "# ", topt, ' = ', val
         Call do_log (1, trim(slog))

      End Subroutine log_opt_i

      Subroutine log_opt_s (optname, val)
         Implicit None
         Character (Len=*) :: optname
         Character (Len=*) :: val
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Character (Len=300) :: slog
         topt = optname
         Call s_cap (topt)

         Write (slog, '(A,A,A,I6)') "# ", topt, ' = ', trim (val)
         Call do_log (1, trim(slog))

      End Subroutine log_opt_s

      Subroutine log_opt_d (optname, val)
         Implicit None
         Character (Len=*) :: optname
         Real (Kind=8) :: val
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Character (Len=300) :: slog
         topt = optname
         Call s_cap (topt)

         Write (slog, '(A,A,A,G15.7)') "# ", topt, ' = ', val
         Call do_log (1, trim(slog))

      End Subroutine log_opt_d

      Function getopt_s (optlist, optname, defval, defout) Result (res)
         Implicit None
         Type (t_config_option), Pointer :: optlist
         Character (Len=*) :: optname
         Character (Len=*) :: defval
         Character (Len=256) :: res
         Logical, Optional, Intent(out) :: defout
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Character (Len=300) :: slog, cmnt
         Type (t_config_option), Pointer :: ptr
         Integer :: def
         def = 1
         topt = optname
         ptr => optlist
         Call s_cap (topt)
         res = defval
         Do while (associated(ptr))
            If (ptr%name(1:ptr%nl) == topt) Then
               res = trim (adjustl(ptr%value(1:ptr%vl)))
               def = 0
               Exit
            End If
            ptr => ptr%next
         End Do
         cmnt = ''
         If (def == 1) cmnt = '   # default value'
         Write (slog, '(A,A,A,A)') topt, ' = ', trim (res), trim (cmnt)
         Call do_log (1, trim(slog))
         If (present(defout)) defout = (def .eq. 1)
      End Function getopt_s


      Function getopt_i (optlist, optname, defval, defout) Result (res)
         Implicit None
         Type (t_config_option), Pointer :: optlist
         Character (Len=*) :: optname
         Integer :: defval, res
         Logical, Optional, Intent(out) :: defout
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Character (Len=300) :: slog, cmnt
         Type (t_config_option), Pointer :: ptr
         Integer :: def
         def = 1
         topt = optname
         ptr => optlist
         Call s_cap (topt)
         res = defval
         Do while (associated(ptr))
            If (ptr%name(1:ptr%nl) == topt) Then
               Read (ptr%value(1:ptr%vl),*) res
               def = 0
               Exit
            End If
            ptr => ptr%next
         End Do
         cmnt = ''
         If (def == 1) cmnt = '   # default value'
         Write (slog, '(A,A,I6,A)') topt, ' = ', res, trim (cmnt)
         Call do_log (1, trim(slog))
         If (present(defout)) defout = (def .eq. 1)
      End Function getopt_i

      Function getopt_d (optlist, optname, defval, defout) Result (res)
         Implicit None
         Type (t_config_option), Pointer :: optlist
         Character (Len=*) :: optname
         Real (Kind=8) :: defval, res
         Logical, Optional, Intent(out) :: defout
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Type (t_config_option), Pointer :: ptr
         Integer :: def
         Character (Len=300) :: slog, cmnt
         topt = optname
         ptr => optlist
         Call s_cap (topt)
         res = defval
         def = 1
         Do while (associated(ptr))
            If (ptr%name(1:ptr%nl) == topt) Then
               If (index(ptr%value(1:ptr%vl), '.') == 0) Then
                  ptr%value (ptr%vl+1:ptr%vl+2) = '.0'
                  ptr%vl = ptr%vl + 2
               End If
               Read (ptr%value(1:ptr%vl),*) res
               def = 0
               Exit
            End If
            ptr => ptr%next
         End Do
         cmnt = ''
         If (def == 1) cmnt = ' # default value'
         Write (slog, '(A,A,G15.7,A)') topt, ' = ', res, trim (cmnt)
         Call do_log (1, trim(slog))
         If (present(defout)) defout = (def .eq. 1)
      End Function getopt_d

      Function getopt_da (optlist, optname, sz, defval, flag) Result (res)
         Implicit None
         Type (t_config_option), Pointer :: optlist
         Character (Len=*) :: optname
         Integer :: sz
         Real (Kind=8) :: defval (sz), res (sz)
         Logical, Optional :: flag
!!$ Local vars
         Character (Len=len_trim(optname)) :: topt
         Type (t_config_option), Pointer :: ptr
         Character (Len=30) :: strs (sz)
         Character (Len=2) :: fmt1
         Character (Len=30) :: fmt2
         Character (Len=300) :: slog
         Integer :: i, j
         topt = optname
         ptr => optlist
         Call s_cap (topt)
         res = defval
         Do while (associated(ptr))
            If (ptr%name(1:ptr%nl) == topt) Then
               If (split_array_opt(ptr%value(1:ptr%vl), sz, strs) == sz) Then
                  Do i = 1, sz
                     j = len_trim (strs(i))
                     If (index(strs(i) (1:j), '.') == 0) Then
                        strs (i) (j+1:j+2) = '.0'
                        j = j + 2
                     End If
                     Read (strs(i) (1:j),*) res (i)
                  End Do
               Else
                  Write (*,*) 'Wrong value of ', topt, ' option!'
                  If (present(flag)) Then
                     flag = .False.
                     Return
                  Else
                     Stop
                  End If
               End If
               Exit
            End If
            ptr => ptr%next
         End Do
         Write (fmt1, '(I2.2)') sz - 1
         fmt2 = '(A,A,' // fmt1 // '(G14.6,","),G14.6,A1)'
         Write (slog, fmt2) topt, ' = [', (res(i), i=1, sz), ']'

!!$    strange behavior with fmt1
!!$    write(slog,'(A,A,'//fmt1//'(G14.6","),G14.6,A1)') topt,' = [',(res(i),i=1,sz),']'

         Call do_log (1, trim(slog))
         If (present(flag)) Then
            flag = .True.
         End If
      End Function getopt_da


      Function split_array_opt (str, n, ostr) Result (res)
         Implicit None
         Character (Len=*) :: str, ostr (n)
         Integer :: n, res
!!$ Local vars
         Integer :: ll, le, pb, pe, i, m
         ll = len_trim (str)

         res = 0
         If (str(1:1) /= '[' .Or. str(ll:ll) /= ']') Return
         If (process_str(str(2:ll-1), pb, pe) == 0) Return
         pb = pb + 1
         pe = pe + 1
         i = 0
         m = index (str(pb:pe), ',') + pb - 2
         If (m == 0) Return
         Do while (pb <= pe)
            If (process_str(str(pb:m), ll, le) == 0) Exit
            i = i + 1
            If (i > n) Return
            ostr (i) = str (pb+ll-1:pb+le-1)
            pb = m + 2
            m = index (str(pb:pe), ',')
            If (m == 0) Then
               m = pe
            Else
               m = m + pb - 2
            End If
         End Do
         If (i == n) res = n
      End Function split_array_opt

      Subroutine dump_config (optlist)
         Implicit None
         Type (t_config_option), Pointer :: optlist
!!$ Local vars
         Type (t_config_option), Pointer :: ptr
         ptr => optlist

         Do while (associated(ptr))
            Write (*,*) ptr%name(1:ptr%nl), ' = ', ptr%value(1:ptr%vl)
            ptr => ptr%next
         End Do
      End Subroutine dump_config

      Subroutine read_config_file (filename, optlist)
         Implicit None
         Type (t_config_option), Pointer :: optlist
         Character (*) :: filename
!!$ Local vars
         Type (t_config_option), Pointer :: ptr
         Integer :: fp = 1234
         Character (Len=255) :: cwork
         Integer :: pb, pe, del, pb1, pe1, pb2, pe2, ef

         Allocate (ptr)
         optlist => ptr
         Open (Unit=fp, File=filename, Action='read')
         Do while (.True.)
            Read (fp, Fmt='(A)', IoStat=ef) cwork
            If (ef /= 0) Exit
            If (process_str(cwork, pb, pe) > 0) Then
               del = index (cwork(pb:pe), '=')
               If (del > 0) Then
                  If (process_str(cwork(pb:pb+del-2), pb1, pe1) /= 0 .And. process_str(cwork(pb+del:pe), pb2, &
                 & pe2) > 0) Then
                     Allocate (ptr%next)
                     ptr => ptr%next
                     Nullify (ptr%next)
                     ptr%name = cwork (pb+pb1-1:pb+pe1-1)
                     Call s_cap (ptr%name)
                     ptr%value = cwork (pb+del+pb2-1:pb+del+pe2-1)
                     ptr%nl = pe1 - pb1 + 1
                     ptr%vl = pe2 - pb2 + 1
                  End If
               End If
            End If
         End Do
         Close (Unit=fp)
         ptr => optlist
         optlist => ptr%next
         Deallocate (ptr)
         Return

      End Subroutine read_config_file

      Function process_str (str, pb, pe) Result (isc)
         Implicit None
         Integer :: pe, pb, isc
         Character (*) :: str
!!$Local vars
         Integer :: i

         isc = 0
         pb = 1
         pe = len_trim (str)
         Do i = 1, pe
            Select Case (str(i:i))
            Case (' ')
               pb = pb + 1
            Case (char(9))
               pb = pb + 1
            Case ('#')
               isc = 0
               Return
            Case Default
               Exit
            End Select
         End Do

         If (pb <= pe) Then
            i = index (str(pb:pe), '#')
            If (i /= 0) pe = pb + i - 2
            Do i = pe, pb, - 1
               Select Case (str(i:i))
               Case (' ')
                  pe = pe - 1
               Case (char(9))
                  pe = pe - 1
               Case Default
                  Exit
               End Select
            End Do
            If (pb <= pe) isc = 1
         End If
      End Function process_str

      Subroutine s_cap (s)
!!$ S_CAP replaces any lowercase letters by uppercase ones in a string.
!!$  Parameters:
!!$    Input/output, character ( len = * ) S, the string to be transformed.
         Implicit None
         Character (Len=*) :: s
!!$Local vars
         Integer i, itemp
         Integer nchar
         nchar = len_trim (s)
         Do i = 1, nchar
            itemp = ichar (s(i:i))
            If (97 <= itemp .And. itemp <= 122) Then
               s (i:i) = char (itemp-32)
            End If
         End Do
         Return
      End Subroutine s_cap

End Module readcfg
