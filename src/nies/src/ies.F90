Program ies

      Use definitions ! types and parameters
!!$    use str_const   ! matrixes for str. constants
      Use green_func ! matrixes for Green's functions
      Use read_write ! reads input and write output files
      Use bzgrid ! 2D BZ
      Use ies_structure ! structure constants
      Use maco_mod ! calculate Madelung constants
      Use pot_mod ! potentials and potential parameters stuff (mixing, transformation et cet)
      Use green_mod ! stuff for all the Green's functions + CPA
      Use dens_mod ! calculate densities of core and valense states
      Use sgfva_mod
      Use logging
!!$   use quad3_mod   ! interpolation and other stuff used for solving Schroedinger and Dirac equations and so on

      Implicit None
!!$ global variables
      Type (params) :: opt
      Type (atoms_set) :: system
      Type (complex_nodes) :: cnodes
      Type (bzone) :: bz
      Type (sharm) :: harmon
      Type (madelung) :: madcon
      Type (andmix) :: mixing
      Integer :: iiter, ncase
      Integer, Parameter :: gflst = 84643
      Integer :: grecl, precl

!!$ PREPERATION STAGE
      Call rall (opt, system, cnodes)!   reading of all data
      flush (IW6)
      Call begi (opt, system, harmon, cnodes)!   start calculations
      flush (IW6)
!!$      Call gibz (bz, system%vbr, system%vbg, opt%nk, opt%nsym, opt%inve)!   2DBZ stuff

      call gibz (bz, system%vbr, opt%bz)
      flush (IW6)
      Call tbk (opt, system, bz, harmon)!   calculation of structure constants for the system
      flush (IW6)

      Call maco (opt, system, madcon)!   Madelung constants
      flush (IW6)
      Call slda (opt, system, mixing)!   start of LDA calculations
      flush (IW6)
!!$  allocation of GREEN'S FUNCTIONS and their parts for bulk(gabu),vacuum(gava) and scat.region
!!$  major memory bottle neck. Kept in memory as a multidimesional arrays while only 2D blocks
!!$  defined by the first 2 dimensions are used for all matrix operations: inversion, matrix
!!$  multiplications, matrix symmetrization (see symgf) et cetera.
      Allocate (zbgam(opt%nblsq, opt%nblsq, bz%nkgrid, opt%ne, system%bns))
      Allocate (zvgam(opt%nblsq, opt%nblsq, bz%nkgrid, opt%ne, system%vns))
      Allocate (zgfg(opt%nlsq, opt%nlsq, system%nums, opt%ne, opt%ns))
      Allocate (zcpf(opt%nlsq, opt%nlsq, system%nums, opt%ne, opt%ns))
      Allocate (zomg(opt%nlsq, opt%nlsq, system%nums, opt%ne, opt%ns))
      Allocate (zcagf(opt%nlsq, opt%nlsq, system%num, opt%ne, opt%ns))
!!$  end of allocation
      ncase = 0

      grecl = size (zbgam) * prec *2

      If (opt%icont == 1 .And. opt%rwgamma == 1) Then
         Open (Unit=gflst, File='gamma.dat', Action='read', Form='UNFORMATTED', Access='direct', Recl=grecl)
         Read (gflst, Rec=1) zvgam
         Read (gflst, Rec=2) zbgam
         Close (Unit=gflst)
         precl = size (zgfg) * prec * 2
         Open (Unit=gflst, File='rest1.dat', Action='read', Form='UNFORMATTED', Access='direct', Recl=precl)
         Read (gflst, Rec=1) zgfg
         Read (gflst, Rec=2) zcpf
         Read (gflst, Rec=3) zomg
         Close (Unit=gflst)

         precl = size (zcagf) * prec * 2
         Open (Unit=gflst, File='rest2.dat', Action='read', Form='UNFORMATTED', Access='direct', Recl=precl)
         Read (gflst, Rec=1) zcagf
         Close (Unit=gflst)

      Else
         Call gamma (opt, system, bz%nkgrid, bz%kweight, cnodes%zcn, ncase, system%bns, system%bnum)!----- GAMMA OF LEFT BULK
         flush (IW6)
!!$      stop
         If (opt%ivac == 1 .And. opt%lident == 0) Then
            ncase = 1
            Call gamma (opt, system, bz%nkgrid, bz%kweight, cnodes%zcn, ncase, system%vns, system%vnum)!---- GAMMA OF RIGHT BULK
            flush (IW6)
         End If
         If (opt%rwgamma == 1) Then
            Open (Unit=gflst, File='gamma.dat', Action='write', Form='UNFORMATTED', Access='direct', &
           & Recl=grecl)
            Write (gflst, Rec=1) zvgam
            Write (gflst, Rec=2) zbgam
         End If
      End If


!!$ END OF PREPARATION STAGE
!!$ INTERATIVE STAGE
!SMP$ DO SERIAL 
      Do iiter = 1, opt%niter
         If (opt%ivac == 0) Call sgfva (opt, system, bz%nkgrid, cnodes%zcn, iiter)!----- UPDATE OF VACUUM SGF
         Call scpf (opt, system, cnodes%zcn, iiter)!----- COHERENT POTENTIAL FUNCTION
         flush (IW6)
         Call repa (opt, system, bz%nkgrid, bz%kweight)!----- RECURSION PARTITIONING-MOST EXPENSIVE PART OF CALC.
         flush (IW6)
         Call cpait (opt, system, cnodes%zcn, iiter)!----- CPA STEP
         flush (IW6)
         Call cmom (opt, system, cnodes%zcn, cnodes%zcw)!----- ENERGY MOMENTS
         flush (IW6)
         Call code (opt%irel, opt%nsirk, opt%thresh, opt%ns, system)!----- CORE DENSITIES
         flush (IW6)
         Call vade (opt%ns, opt%nl, system)!----- VALENCE DENSITIES
         flush (IW6)
         Call dimo (opt, system, harmon)!----- DIPOLE MOMENTS
         flush (IW6)
         Call nepo (opt, system, madcon, mixing, opt%nam, iiter)!----- NEW POTENTIALS
         flush (IW6)
         Call tisk (opt, system, iiter)!----- PRINT OF RESULTS
         flush (IW6)
      End Do
!!$ END OF INTERATIVE STAGE
      Call tisk (opt, system, iiter)!----- PRINT OF RESULTS

      Write (IW6, 191)
      flush (IW6)
191   Format (/ '        * * *   END OF IES   * * * ' /)

      Stop
End Program ies
