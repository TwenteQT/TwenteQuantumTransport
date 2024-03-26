!!#   $Id: band.F90 149 2006-10-26 16:29:08Z antst $

#include "math_def.h"

Program fermisurf
   Use transport
   Use ando_module
   Use atoms_module
   Use geometry_module
   Use structure
   Use supercell
   Use hamiltonian
   Use logging
   Use sparse_solvers
   Use adaptive_bzgrid
   Use readcfg
   Use rotations
!!$ EMTO mods
   Use omta_kink
   Use omta_defs
   Use omta_strrs
   Use omta_pots
   Use omta_SOC
   Use df

   Use mpi
   Implicit None
   !Include 'mpif.h'

   Integer, Parameter :: clen = 200

   Type(t_atoms_set) :: atoms
   Type(t_fsopts) :: opt
   Type(t_geometry) :: bgeo, tbgeo
   Type(t_strconst) :: bsc
   Type(t_mathop) :: msys, bsk
   Type(t_potpar_mats) :: bpar
   Type(t_ando_options) :: aopts
   Type(t_ando_sollution) :: ando
   Type(adapt_bzone) :: bz
   Real(Kind=DEF_DBL_PREC), Pointer :: k(:, :), ek(:, :)
   Real(Kind=DEF_DBL_PREC) :: kpar(2)
   Real(Kind=DEF_DBL_PREC) :: en, pi, rydberg
   Integer :: is, n, nh, j, nmod, i
   Integer :: bandf(2) = (/33, 34/)
   Complex(Kind=DEF_DBL_PREC), Pointer :: ons(:, :)
   Complex(Kind=DEF_DBL_PREC), Pointer :: ofs(:, :)
   Complex(Kind=DEF_DBL_PREC), Pointer :: hops(:, :)
   Complex(Kind=DEF_DBL_PREC), Pointer :: ham(:, :)
   Complex(Kind=DEF_DBL_PREC) :: ph
   Real(Kind=DEF_DBL_PREC), Pointer :: rotm(:, :)
   Integer :: npars, ik, igp, iloop
   Real(Kind=DEF_DBL_PREC) :: kk, k1, dffact
   Integer :: nb, nc, ngc(3), ngr(3)
   Complex(Kind=DEF_DBL_PREC) :: fact
   Integer :: lwork, info, mdf, nd, pohops
   Complex(Kind=DEF_DBL_PREC), Pointer :: work(:)
   Real(Kind=DEF_DBL_PREC), Pointer :: rwork(:), ev(:), kdf(:)
   Integer, Pointer :: iwork(:)
   Character(Len=clen) :: cwork
!!$ MPI vars
   Integer :: ierr, my_mpi_id, local_comm, solve_comm, mpi_sz, root = 0, solsize
   Integer :: mpi_loc_id
!!$ EMTO vars
   real(kind=DEF_DBL_PREC) :: bzk(2)
   Type(t_geometry_EMTO) :: bgeo_EMTO, tbgeo_EMTO
   Type(t_omta_struct) :: bsc_EMTO(2)
   Type(t_omta_logder) :: bpar_EMTO(2)
   Complex(kind=DEF_DBL_PREC), Pointer :: ofs_sc(:, :)
   Complex(kind=DEF_DBL_PREC), Pointer :: ons_soc(:, :)
   Type(zcsrmat) :: hsoc, htemp, htemp2
   Type(zcsrmat) :: zcsr_temp, zcsr_temp2, zcsr_temp3
   Type(zcsrmat) :: zcsr2_temp, zcsr2_temp2
   Type(zcsrmat) :: zcsr_off, zcsr_off2, zcsr_off3, zcsr_null
   Type(t_mathop_EMTO) :: tempsys(2)
   Type(t_atoms_set_EMTO) :: atoms_EMTO
   Integer :: ldimidim(2)
   Real(kind=DEF_DBL_PREC) :: kappa2
   Type(zcsrmat) :: temp
!!$ janky array: indices of states that are to be downfolded
   integer :: ndim, lidim

   pi = 4.d0*datan(1.d0)
   rydberg = 0.5d0*27.211396d0

   Call mpi_init(ierr)
   solve_comm = mpi_comm_world
   Call mpi_comm_rank(solve_comm, my_mpi_id, ierr)
   Call mpi_comm_size(solve_comm, mpi_sz, ierr)

   Log_Level = -1
   If (my_mpi_id == root) Then
!!$          Call unlink ('outp')
      Log_Level = 1
   End If

#ifdef _VERS_
   Call do_log(1, 'Fermi surface code v'//_VERS_//' by Msrang')
#else
   Call do_log(1, 'Fermi surface code vXXXX by Msrang')
#endif

   Call fs_config_init(opt)
!!$       Log_Level = opt%loglvl
   Call MPI_COMM_SPLIT(solve_comm, my_mpi_id, 0, local_comm, ierr)
   solsize = 1

   Call mpi_comm_rank(local_comm, mpi_loc_id, ierr)

!!$       Allocate (jobs(mpi_sz))
   If (my_mpi_id == root) Then
      Log_Level = opt%loglvl
      Open (Unit=bandf(1), File='plotFermi_up.b', Action='write')
      Open (Unit=bandf(2), File='plotFermi_down.b', Action='write')
   Else
!!$ uncomment line bellow if you want output from rest of CPUs
!!$     Log_Level=opt%loglvl
      Write (cwork, '(i2.2)') my_mpi_id
      Write (logfilename, '("log.",A)') trim(cwork)
   End If

   If (mpi_loc_id == 0) Then

      If (opt%po%kind <= 3) Then
         atoms = read_atoms('atomlist')
         tbgeo = read_geom('geom_f', atoms)


         pohops = 1
         If (opt%po%kind == 3) pohops = 2

         bgeo = make_leadgeom(tbgeo, pohops, 1)

         Call free_geom(tbgeo)

         rotm => prep_rot_mask(bgeo, opt%po)

         bsc = calc_screal(bgeo, 1)

         !!$ Prepare grid for BZ integration
         !Call gibz(bz, bgeo%base, opt%bzo)
         Call init_bz(bz, bgeo%base, opt%bzo)
      End If

!!$ ================ Added option: EMTO ==================

      If (opt%po%kind == 6) then
         pohops = 1
         call read_omta_pot('atomlist_EMTO', atoms_EMTO)
         tbgeo_EMTO = read_geom_EMTO('geom_f', atoms_EMTO, .true.)
         bgeo_EMTO = make_leadgeom_EMTO(tbgeo_EMTO, 1, 1)
         call free_geom_EMTO(tbgeo_EMTO)

         call makidx_lead(bgeo_EMTO, lidim)

         !!$ Prepare grid for BZ integration
         !Call gibz(bz, bgeo_EMTO%base, opt%bzo)
         Call init_bz(bz, bgeo_EMTO%base, opt%bzo)

      End if

!!$ ======================================================

      en = opt%en

      aopts%use_ev_solver = opt%actl%gensol
      aopts%dir = 1
      aopts%usestates = 1
      aopts%needEmb = 0
      aopts%needBound = 0
      aopts%refine = opt%actl%refine
      aopts%usestates = opt%actl%specpart
      aopts%need_split = 1

      If (opt%po%kind .Le. 3) Then

         If (opt%sdcdir /= 0) Then
            aopts%rotm = spinrm(rotm(:, 1:1))
         End If

         Call prep_dec_idx(bgeo, aopts%split_idx)

         n = bsc%nrows/pohops
         nh = pohops

         Allocate (ons(bsc%nrows, bsc%nrows))
         Allocate (ofs(bsc%nrows, bsc%nrows))
         Allocate (hops(n, n*(nh + 1)))
         Allocate (ham(n, n))
         lwork = 4*n
!!$        lwork = 2 * n ** 2 + 5 * n + 1
         Allocate (work(lwork))
         Allocate (rwork(lwork))
         Allocate (iwork(lwork))
         Allocate (ev(n))

         Allocate (k(bsc%nrows + 2, 2))
         Allocate (ek(bsc%nrows + 1, 2))
         ph = Exp(DEF_cmplx_Ione*DEF_M_PI)
         npars = 0

         Call calc_potpar(atoms, opt%atopt, en)
!!$ Prepare potential parameters matrices for leads
         bpar = make_ppar_full(bgeo, opt%po, rotm)

      Else If (opt%po%kind .Eq. 5) Then

         pohops = 1
         aopts%usestates = 1.0d0
         nb = 2
         nc = 2*nb
         n = nc*opt%po%ngrid
         nh = pohops

         If (aopts%need_split == 1) Call prep_ham_idx(opt%po%ngrid, nb, aopts%split_idx)

         If (opt%po%rot_mask /= 0) Then
            rotm => read_rot_mask(n)
         Else
            Allocate (rotm(2, n))
            rotm(:, :) = 0.d0
         End If

         Allocate (ons(n, n))
         Allocate (ofs(n, n))
         Allocate (hops(n, n*(nh + 1)))

         lwork = 4*n
         Allocate (work(lwork))
         Allocate (rwork(lwork))
         Allocate (iwork(lwork))
         Allocate (ev(n))
         Allocate (k(n + 2, 2))

!!$ ================ Added option: EMTO ==================
      Else if (opt%po%kind .Eq. 6) Then

         Call prep_dec_idx_EMTO(bgeo_EMTO, aopts%split_idx)
         n = bgeo_EMTO%nspin*bgeo_EMTO%norbit - size(bgeo_EMTO%to_downfold)
         nh = 1

         ndim = bgeo_EMTO%nspin*bgeo_EMTO%norbit

         Allocate (ons(n, n))
         Allocate (ofs(n, n))
         Allocate (hops(n, 2*n))
         !Allocate (k(n/bgeo_EMTO%nspin + 2 - size(bgeo_EMTO%to_downfold)/2, 2))
         Allocate (k(n + 2, 2))
         ph = Exp(DEF_cmplx_Ione*DEF_M_PI)

!!$ ======================================================

      End If

      mdf = pohops
      dffact = 1.0d0/real(mdf, kind=8)
!!$      Allocate (mine(mdf))
      Allocate (kdf(mdf))

   End If

   !nkloops = ceiling(bz%nkgrid/real(ngrps, kind=DEF_DBL_PREC)) for MPI parallization over k-points, later (maybe)

   en = opt%en

   iloop = 0
   Do While (.not. bz%conv)
    iloop = iloop + 1
   ! if (iloop == 3) stop
    Do igp = 1, size(bz%nc_gp) 
       ik = bz%nc_gp(igp)

       kpar = bz%bz(ik)%kpar
       write(cwork,'("Calculating igp, ik=",i7,1x,i7,1x,2f8.4)') igp,ik,kpar
       call do_log(2, trim(cwork))

       If (opt%po%kind .Le. 3) Then

          Call make_sk(bsk, bsc, kpar, 1)
          Call prep_system(msys, bsk, bpar, opt%po, 1)

       Else If (opt%po%kind == 5) Then

          opt%po%energy = en
          ngc(1) = opt%po%ngridx
          ngc(2) = opt%po%ngridy
          ngc(3) = opt%po%ngridz
          ngr(:) = ngc(:)

          Call make_sk(bsk, bsc, kpar, 1)
          Call prep_system(msys, opt%po, ngc, ngr, kpar, rotm)

!!$  ================ Added option: EMTO ==================

       Else if (opt%po%kind == 6) then

          atoms_EMTO%fermi = en
          ldimidim(1) = lidim
          ldimidim(2) = lidim

          kappa2 = (en - atoms_EMTO%vmtz(1))*bgeo_EMTO%dawsr*bgeo_EMTO%dawsr*bgeo_EMTO%alat*bgeo_EMTO%alat
          call mkalph(kappa2, bgeo_EMTO%nlmax, bgeo_EMTO%dawsr*bgeo_EMTO%alat, atoms_EMTO)
          call mktral(kappa2, bgeo_EMTO%nlmax, bgeo_EMTO%dawsr*bgeo_EMTO%alat, atoms_EMTO)
          call atompp_nmto(atoms_EMTO)
          call makpph_nmto(atoms_EMTO%nspin, bgeo_EMTO, ldimidim, bpar_EMTO, .true., opt%po%so, opt%po%nonsph)

          bsc = calc_screal_EMTO(bgeo_EMTO, kappa2, en, irel=1, needcur=0, df=opt%po%df, ppar=bpar_EMTO(1))

          call make_sk(bsk, bsc, kpar, 1)

          if (opt%po%nonsph > 0) call calc_omta_nonsph(atoms_EMTO)
          ! remove downfolded states from bpar(1)%logder
          if (opt%po%df > 0) call downfold_logder(bgeo_EMTO, bpar_EMTO(1)%logder)
          if (opt%po%df > 0) call downfold_logder(bgeo_EMTO, bpar_EMTO(1)%hsoc)

          call prep_system(msys, bsk, bpar_EMTO(1), opt%po, 1)
!!$  ======================================================

       End If

       ons = sptofull(msys%c)
       ofs = sptofull(msys%r)

       Call free_mathop(msys)

       hops(1:n, 1:n) = ons(1:n, 1:n)
       hops(1:n, n + 1:n + nh*n) = ofs(n*nh - n + 1:n*nh, 1:nh*n)

!!$  Solve ando problem
       Call solve_ando(ofs, hops, n, nh, ando, aopts)

       Do is = 1, 2

          nmod = ando%Snin(is) + ando%Snout(is)
          !if (nmod > 1) STOP 'This code does not support multiple modes per k-point yet'
          if (nmod==0) Then 
            bz%bz(ik)%kz = 0
          else
            bz%bz(ik)%kz =(imag(Log(ando%lin(ando%mask_in(1, is))))) 
          Endif

          k(1, is) = en
          k(2, is) = real(nmod, kind=DEF_DBL_PREC)
          Do i = 1, nmod, 1
             if (i > ando%Snin(is)) then
                k(i + 2, is) = (imag(Log(ando%lout(ando%mask_out(i-ando%Snin(is), is)))))
             else 
                k(i + 2, is) = (imag(Log(ando%lin(ando%mask_in(i, is)))))
             endif

             kk = k(i + 2, is)*dffact
          End Do

          If (opt%po%kind .Le. 3) Then
             k(nmod + 3:bsc%nrows + 2, is) = 0.0 !!$/ 0.0
          Else If (opt%po%kind > 4) Then
             k(nmod + 3:bsc_EMTO(1)%ldim + 2 - size(bgeo_EMTO%to_downfold)/2, is) = 0.0 !!$/ 0.0
          End If

       End Do

       If (opt%po%kind .Le. 3) Then
          Call log_bands(solve_comm, root, kpar, k, bsc%nrows + 2, bandf)
       Else If (opt%po%kind > 4) Then
          Call log_bands(solve_comm, root, kpar, k, n + 2 , bandf)
       End If

    End Do ! k-grid loop

    Call refine_bz(bz, opt%bzo%tol)

   End Do ! refinement loop

   If (opt%po%kind .Le. 3) Call free_ppar_f(bpar)

   If (my_mpi_id == root) Then
      Close (Unit=bandf(1))
      Close (Unit=bandf(2))
   End If

   Call MPI_Barrier(solve_comm, ierr)

   Call mpi_finalize(ierr)

   Call do_log(1, 'Done!')

Contains
   Subroutine log_bands(comm, root, kpar, k, n, bandf)
      Implicit None
      Integer :: root, bandf(2), n, comm
      Real(Kind=DEF_DBL_PREC) :: kpar(2), k(n, 2)
!!$ Locals
      Integer :: ierr, procnum, mpi_sz, my_mpi_id
      Integer :: status(MPI_STATUS_SIZE)
      Integer :: is, j, nmod
      Call mpi_comm_rank(comm, my_mpi_id, ierr)
      Call mpi_comm_size(comm, mpi_sz, ierr)
      If (my_mpi_id /= root) Then
         Call MPI_Send(k, n*2, MPI_DOUBLE_PRECISION, root, my_mpi_id + mpi_sz + 1000, comm, ierr)
      Else
         Do is = 1, 2
            nmod = Int(k(2, is))
            Do j = 1, nmod
                Write (bandf(is), '(3(1x,f12.8))') kpar(1), kpar(2),  k(j + 2, is)
            End Do
         End Do
         Do procnum = 2, mpi_sz
            Call MPI_Recv(k, n*2, MPI_DOUBLE_PRECISION, procnum - 1, procnum + mpi_sz + 999, comm, status, &
           & ierr)
            Do is = 1, 2
               nmod = Int(k(2, is))
               Do j = 1, nmod
                 Write (bandf(is), '(3(1x,f12.8))') kpar(1), kpar(2),  k(j + 2, is)
               End Do
            End Do
         End Do
         flush (bandf(1))
         flush (bandf(2))
      End If

   End Subroutine log_bands

End Program fermisurf
