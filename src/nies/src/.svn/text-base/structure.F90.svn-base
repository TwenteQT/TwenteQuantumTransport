Module ies_structure
      Use definitions
      Use structure
Contains


      Subroutine begi (opt, system, harmon, cnodes)
!************************************
!!$   TESTS OF INPUT DATA AND
!!$   STARTING AUXILIARY CALCULATIONS
!************************************
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Type (sharm) :: harmon
         Type (complex_nodes) :: cnodes
         Type (atom_definition), Pointer :: aa
         Type (bulk_atom_definition), Pointer :: ab
         Real (Kind=prec), Parameter :: TOL = 1.0d-6, EPSCOR = 1.0D-9
         Real (Kind=prec), Parameter :: QSC3 (3) = (/ 0.3484850d0, 0.05303d0, 0.010714d0 /)! screened constants for spd basis set
         Real (Kind=prec), Parameter :: QSC4 (4) = (/ 0.3850570d0, 0.073209D0, 0.022481D0, 0.006069D0 /)! screened constants for spdf basis set
         Integer :: i, j, ia, is, ip, ib, ic, ig, il, ie, l
         Real (Kind=prec) :: sum, pi, det, prod, dif, fact, hdba, s, s2, s4, sqs, sqs3, tl1, tl3, tl5, tl7, &
        & sqtl3, cnorm

         pi = 4.0d0 * datan (1.0d0)

!----------------------------------- SCREENING CONSTANTS
         Allocate (system%qscr(opt%nl))
         If (opt%nl == 3) Then
            system%qscr = QSC3
         Else
            system%qscr = QSC4
         End If

         Write (IW6, 114) (system%qscr(il), il=1, opt%nl)
114      Format (/ 3 X, 'SCREENING CONSTANTS: ' / 5 X, 4 G15.7)

!---------------------------------- PRECISION OF CORES
         opt%THRESH = EPSCOR
         Write (IW6, 115) opt%THRESH
115      Format (/ 9 X, ' PRECISION OF CORE ENERGIES: ', G15.7)

!------------------------- ACCURACY OF RUNGE-KUTTA METHOD (1.GT.NSIRK.LT.5)
         opt%NSIRK = 2
         Write (IW6, 116) opt%NSIRK
116      Format (/ 9 X, ' ACCURACY OF RUNGE-KUTTA METHOD: ', '  NSIRK=', I2)

!---------------------------------- MAX. NO. OF SGF ITERATIONS
         opt%nmits = 500
         Write (IW6, 117) opt%nmits
117      Format (/ 3 X, ' MAX. NUMBER OF SGF ITERATIONS: ', I5)

!---------------------------------- SHifT OF ENY'S
         Do ia = 1, system%num
            aa => system%at (ia)
            Allocate (aa%isheny(opt%nl, opt%ns))
            Do is = 1, opt%ns
               Do il = 1, opt%nl
                  aa%isheny (il, is) = 1
               End Do
            End Do
         End Do

         Write (IW6, 118)
118      Format (/ 2 X, '  LABEL,             IA,  IS,', '    ISHENY - FOR SHIFT OF ENY:')

         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%ns
               Write (IW6, 119) aa%otxta, ia, is, (aa%isheny(il, is), il=1, opt%nl)
            End Do
         End Do
119      Format (2 X, A16, 2 X, 2 I5, 8 X, 4 I2)


!---------------------------------- GAUNT FACTORS and normal constants
         harmon%cnh => harm_nc (opt%nl)
         harmon%gfrh => gaunty (opt%nl)

!------------------ SYSTEM FERMI ENERGY AND SCALING OF COMPLEX CONTOUR  NODES AND WEIGHTS

         If (opt%ivac == 0) Then
            system%ef = system%bef
         Else
            system%ef = (system%bef+system%vef) / 2.0d0
         End If

         Do ie = 1, opt%ne
            cnodes%zcn (ie) = opt%diam * cnodes%zcn(ie) + system%ef
            cnodes%zcw (ie) = opt%diam * cnodes%zcw(ie)
         End Do
         Write (IW6, 120) system%ef
120      Format (/ 5 X, '****  SYSTEM FERMI ENERGY= ', G15.7)
         Write (IW6, 121) cnodes%zcn(1), cnodes%zcn(opt%ne)
121      Format (/ 5 X, 'ENERGY CONTOUR: FIRST NODE= ', 2 G15.7 / 5 X, '                 LAST NODE= ', 2 &
        & G15.7)


!---------------------------- CHECK OF CONSISTENCY OF CONCENTRATIONS
         ia = 0
180      Format (/ ' **** INPUT ERROR: ', ' INCONSISTENCY IN CONCENTRATIONS ')
         Do ig = 1, system%nums
            sum = 0.0d0
            Do ic = 1, system%nc(ig)
               ia = ia + 1
               aa => system%at (ia)
               If (aa%con < 0.0d0 .Or. aa%con > 1.0d0) Then
                  Write (IW6, 180)
                  Stop
               End If
               sum = sum + aa%con
            End Do
            If (dabs(sum-1.0d0) > TOL) Then
               Write (IW6, 180)
               Stop
            End If
         End Do

         sum = 0.0d0
181      Format (/ ' **** INPUT ERROR: ' / 3 X, 'INCONSISTENCY IN CONCENTRATIONS OF BULK SUBSTRATE')
         Do ia = 1, system%bnum
            ab => system%bat (ia)
            If (ab%con < 0.0d0 .Or. ab%con > 1.0d0) Then
               Write (IW6, 181)
               Stop
            End If
            sum = sum + ab%con
         End Do
         If (dabs(sum-1.0d0) > TOL) Then
            Write (IW6, 181)
            Stop
         End If


182      Format (/ ' **** INPUT ERROR: ' / 3 X, 'INCONSISTENCY IN CONCENTRATIONS OF VAC. SUBSTRATE')
         If (opt%ivac == 1) Then
            sum = 0.0d0
            Do ia = 1, system%vnum
               ab => system%vat (ia)
               If (ab%con < 0.0d0 .Or. ab%con > 1.0d0) Then
                  Write (IW6, 182)
                  Stop
               End If
               sum = sum + ab%con
            End Do
            If (dabs(sum-1.0d0) > TOL) Then
               Write (IW6, 182)
               Stop
            End If
         End If


!----------------------------------- BASIS IN RECIP. SPACE

         det = system%vbr (1, 1) * system%vbr(2, 2) - system%vbr(1, 2) * system%vbr(2, 1)
         system%romega = dabs (det)
         system%gomega = (2.0d0*pi) ** 2 / system%romega
         cnorm = 2.0d0 * pi / det

         system%VBG (1, 1) = cnorm * system%vbr(2, 2)
         system%VBG (2, 1) = - cnorm * system%vbr(1, 2)
         system%VBG (1, 2) = - cnorm * system%vbr(2, 1)
         system%VBG (2, 2) = cnorm * system%vbr(1, 1)

         Write (IW6, 124) system%romega, system%gomega
124      Format (/ ' AREA OF 2D-PRIMITIVE CELL: REAL= ', G15.7 / '                      RECIPROCAL= ', G15.7)
         Write (IW6, 125) ((system%VBG(i, j), i=1, 2), j=1, 2)
125      Format (/ ' RECIP. BASIS: 1. VECTOR = ', 2 G15.7 / '               2. VECTOR = ', 2 G15.7)

!------------- CHECK OF CONSISTENCY OF DIRECTIONS OF BULK AND VACUUM TRANSL. VECTORS
191      Format (/ ' **** INPUT ERROR: ' / 3 X, 'VACUUM AND BULK TRANSLATION VECTORS ARE INCONSISTENT')
         prod = system%BTRV (3) * system%VTRV(3)
         If (prod >= 0.0d0) Then
            Write (IW6, 191)
            Stop
         End If

!------------------------ DIMENSIONLESS WS-RADII OF BULK AND VACUUM ----------------
         prod = system%romega * dabs (system%BTRV(3)) / dble (system%NB)
         system%abws = ((3.0d0*prod)/(4.0d0*pi)) ** (1.0d0/3.0d0)
         prod = system%romega * dabs (system%VTRV(3)) / dble (system%NB)
         system%avws = ((3.0d0*prod)/(4.0d0*pi)) ** (1.0d0/3.0d0)

         Write (IW6, 126) system%abws, system%avws
126      Format (/ ' DIMENSIONLESS WS-RADII:  BULK= ', G15.7 / '                        VACUUM= ', G15.7)

!------------------------ TRUE AVERAGE WS-RADII OF BULK AND VACUUM  ----------------
         sum = 0.0d0
         Do ia = 1, system%bnum
            ab => system%bat (ia)
            sum = sum + ab%con * ab%wsav ** 3
         End Do
         system%bwst = sum ** (1.0d0/3.0d0)

         If (opt%ivac == 0) Then
            system%vwst = system%bwst * system%avws / system%abws
         Else
            sum = 0.0d0
            Do ia = 1, system%vnum
               ab => system%vat (ia)
               sum = sum + ab%con * ab%wsav ** 3
            End Do
            system%vwst = sum ** (1.0d0/3.0d0)
         End If


         Write (IW6, 127) system%bwst, system%vwst
127      Format (/ ' TRUE AVERAGE WS-RADII:  BULK= ', G15.7 / '                       VACUUM= ', G15.7)

         dif = (system%bwst*system%avws) / (system%vwst*system%abws) - 1.0d0

         If (dabs(dif) > TOL*10.0d0) Then
            Write (IW6, 192)
            Stop
         End If
192      Format (/ ' **** INPUT ERROR: ' / 'VACUUM AND BULK AVERAGE DIMENSIONLESS AND TRUE ' / '   WS-RADII A&
        &RE INCONSISTENT')

!-----------IMENSIONLESS AVERAGE WS-RADII IN INTERMEDIATE REGION

         fact = system%abws / system%bwst
         ia = 0
         ig = 0
         Allocate (system%aws(system%NB, system%np))
         Do ip = 1, system%np
            Do ib = 1, system%NB
               ig = ig + 1
               sum = 0.0d0
               Do ic = 1, system%nc(ig)
                  ia = ia + 1
                  aa => system%at (ia)
                  sum = sum + aa%con * aa%wsav ** 3
               End Do
               system%aws (ib, ip) = fact * sum ** (1.0d0/3.0d0)
            End Do
         End Do

         Write (IW6, 128)
128      Format (/ 5 X, 'DIMENSIONLESS AVERAGE WS-RADII:' /)
         Do ip = 1, system%np
            Do ib = 1, system%NB
               Write (IW6, 129) ip, ib, system%aws(ib, ip)
            End Do
         End Do
129      Format (8 X, 'IP=', I4, 5 X, 'IB=', I2, 5 X, 'AWS=', G15.7)

!-------------- POTENTIAL PARAMETERS OF THE SUBSTRATES CASE OF SURFACE

         If (opt%ivac == 0) Then
            system%vnum = 1
            system%vns = 1
            Allocate (system%vat(1))
            ab => system%vat (1)
            Allocate (ab%eny(opt%nl, opt%ns))
            Allocate (ab%ppc(opt%nl, opt%ns))
            Allocate (ab%ppd(opt%nl, opt%ns))
            Allocate (ab%ppq(opt%nl, opt%ns))
            Allocate (ab%ppp(opt%nl, opt%ns))
            Allocate (ab%dny(opt%nl, opt%ns))
            Allocate (ab%finy(opt%nl, opt%ns))
            Allocate (ab%finyd(opt%nl, opt%ns))

            ab%con = 1.0d0
            ab%az = 0.0d0
            ab%ws = system%vwst
            ab%wsav = system%vwst
            s = system%vwst
            s2 = s ** 2
            s4 = s2 ** 2
            sqs = dsqrt (s)
            sqs3 = s * sqs

            Do il = 1, opt%nl
               l = il - 1
               tl1 = real (2*l+1, kind=prec)
               tl3 = real (2*l+3, kind=prec)
               tl5 = real (2*l+5, kind=prec)
               tl7 = real (2*l+7, kind=prec)
               sqtl3 = dsqrt (tl3)
               ab%eny (il, 1) = system%DBA
               ab%ppc (il, 1) = system%DBA + tl1 * tl5 / (2.0d0*s2)
               ab%ppd (il, 1) = tl5 ** 2 / (2.0d0*4.0d0*tl3*s2)
               ab%ppq (il, 1) = tl5 / (4.0d0*tl1*tl3)
               ab%ppp (il, 1) = s4 / (tl3*tl5**2*tl7)
               ab%dny (il, 1) = real (l, kind=prec)
               ab%finy (il, 1) = sqtl3 / sqs3
               ab%finyd (il, 1) = - sqs / (tl5*sqtl3)
            End Do
         End If

!!$                                     CASE OF INTERFACE
         If (opt%ivac == 1) Then
            system%DBA = system%bef - system%vef
            hdba = 0.5d0 * system%DBA
            Do ia = 1, system%vnum
               ab => system%vat (ia)
               Do is = 1, system%vns
                  Do il = 1, opt%nl
                     ab%eny (il, is) = ab%eny(il, is) + hdba
                     ab%ppc (il, is) = ab%ppc(il, is) + hdba
                  End Do
               End Do
            End Do

            Do ia = 1, system%bnum
               ab => system%bat (ia)
               Do is = 1, system%bns
                  Do il = 1, opt%nl
                     ab%eny (il, is) = ab%eny(il, is) - hdba
                     ab%ppc (il, is) = ab%ppc(il, is) - hdba
                  End Do
               End Do
            End Do
         End If


!---------------- CHECK OF CONSISTENCY OF NSYM, INVE AND 2D-TRANSLATION VECTORS

         If (opt%bz%NSYM < 0 .Or. opt%bz%NSYM > 4) Then
            Write (IW6, 195)
            Stop
         End If
195      Format (/ ' **** INPUT ERROR: NSYM MUST BE 0,1,2,3,4')


         If (opt%bz%NSYM == 0) Go To 222
         If (opt%bz%NSYM == 1) Go To 201
         If (opt%bz%NSYM == 2) Go To 202
         If (opt%bz%NSYM == 3) Go To 203
         If (opt%bz%NSYM == 4) Go To 204

!-------------  NSYM = 1 - SYMMETRY OF CUBIC (001) SURFACES:
!!$   SYMMETRY GROUP GIVEN BY 4-FOLD ROTATION AXIS Z AND BY MIRROR PLANE X-Z

201      If (opt%bz%INVE < 1 .Or. opt%bz%INVE > 2) Go To 296
         If (opt%bz%INVE == 1) Go To 2011
         If (opt%bz%INVE == 2) Go To 2012

!!$         INVE = 1: THE 2D-TRANSL. VECTORS MUST HAVE THE FORM:
!!$                                      A1 = (P,0)
!!$                                      A2 = (0,P)
!!$                               WHERE P is POSITIVE.

2011     If (dabs(system%vbr(2, 1)) > 0.0d0) Go To 297
         If (dabs(system%vbr(1, 2)) > 0.0d0) Go To 297
         If (system%vbr(1, 1) <= 0.0d0) Go To 297
         If (system%vbr(2, 2) <= 0.0d0) Go To 297
         dif = system%vbr (1, 1) / system%vbr(2, 2) - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         Go To 222

!!$                   INVE = 2: THE 2D-TRANSL. VECTORS  MUST HAVE THE FORM:
!!$                                      A1 = (P,P)
!!$                                      A2 = (-P,P)*U
!!$                   WHERE P is POSITIVE AND U=1 OR U=-1.

2012     If (system%vbr(1, 1) <= 0.0d0) Go To 297
         dif = system%vbr (2, 1) / system%vbr(1, 1) - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         If (dabs(system%vbr(1, 2)) <= 0.0d0) Go To 297
         dif = system%vbr (2, 2) / system%vbr(1, 2) + 1.0d0
         If (dabs(dif) > TOL) Go To 297
         dif = dabs (system%vbr(2, 2)) / system%vbr(1, 1) - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         Go To 222

!---  NSYM = 2 - SYMMETRY OF CUBIC (110) SURFACES: SYMMETRY GROUP GIVEN BY TWO MIRROR PLANES X-Z AND Y-Z

!---  NSYM = 4 - SYMMETRY OF CUBIC (N10) SURFACES AND TILT GRAIN BOUNDARIES: SYMMETRY GROUP GIVEN BY MIRROR PLANE X-Z

204      Continue
202      If (opt%bz%INVE < 1 .Or. opt%bz%INVE > 2) Go To 296
         If (opt%bz%INVE == 1) Go To 2021
         If (opt%bz%INVE == 2) Go To 2022

!!$               INVE = 1: THE 2D-TRANSL. VECTORS MUST HAVE THE FORM:
!!$                                      A1 = (P,0)
!!$                                      A2 = (0,Q)
!!$                         WHERE P AND Q ARE POSITIVE.

2021     If (dabs(system%vbr(2, 1)) > 0.0d0) Go To 297
         If (dabs(system%vbr(1, 2)) > 0.0d0) Go To 297
         If (system%vbr(1, 1) <= 0.0d0) Go To 297
         If (system%vbr(2, 2) <= 0.0d0) Go To 297
         Go To 222


!!$           INVE = 2: THE 2D-TRANSL. VECTORS MUST HAVE THE FORM:
!!$                                      A1 = (P,Q)
!!$                                      A2 = (-P,Q)*U
!!$           WHERE P AND Q ARE POSITIVE AND U=1 OR U=-1.

2022     If (system%vbr(1, 1) <= 0.0d0) Go To 297
         If (system%vbr(2, 1) <= 0.0d0) Go To 297
         dif = dabs (system%vbr(1, 2)) / system%vbr(1, 1) - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         dif = dabs (system%vbr(2, 2)) / system%vbr(2, 1) - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         dif = (system%vbr(1, 2)*system%vbr(2, 2)) / (system%vbr(1, 1)*system%vbr(2, 1)) + 1.0d0
         If (dabs(dif) > TOL) Go To 297
         Go To 222

!--------  NSYM = 3 - SYMMETRY OF CUBIC (111)  AND  HEXAGONAL (0001) SURFACES:
!!$                     SYMMETRY GROUP GIVEN BY 3-FOLD ROTATION AXIS Z AND BY MIRROR PLANE Y-Z
!!$          INVE - NOT ACTIVE. THE 2D-TRANSL. VECTORS MUST HAVE THE FORM:
!!$                                   A1 = (P,0)
!!$                                   A2 = (P/2)*(1,S)*U
!!$          WHERE P is POSITIVE,S=SQRT(3) OR S=-SQRT(3),AND U=1 OR U=-1.

203      If (dabs(system%vbr(2, 1)) > 0.0d0) Go To 297
         If (system%vbr(1, 1) <= 0.0d0) Go To 297
         dif = 2.0d0 * dabs (system%vbr(1, 2)) / system%vbr(1, 1) - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         dif = (system%vbr(1, 2)**2+system%vbr(2, 2)**2) / system%vbr(1, 1) ** 2 - 1.0d0
         If (dabs(dif) > TOL) Go To 297
         Go To 222

222      Continue

!!$     DEFINING GLOBAL PARAMETERS AND DIMESIONS FOR DYNAMICALLY ALLOCATED ARRAYS
         opt%NPP1 = system%np + 1
         opt%NPP2 = system%np + 2
         opt%nlsq = opt%nl ** 2
         opt%lmax = 2 * (opt%nl-1)! MAX. ANGULAR QUANTUM NUMBER OF SPHERICAL HARMONICS IN CANONICAL STRUCTURE CONSTANTS
         opt%nharm = (opt%lmax+1) ** 2 ! MAX. NUMBER OF SPHERICAL HARMONICS IN CANONICAL STRUCTURE CONSTANTS
         opt%nblsq = system%NB * opt%nlsq ! DIMENSION OF THE ON-LAYER BLOCKS OF THE K||-DEPENDENT GREEN'S FUNCTION
         opt%mtbcl = MNCL * opt%nlsq
         opt%pair = (opt%nl-1) ** 2 !PARAMETER: MAX. NUMBER OF GF-ELEMENTS FOR THE DIPOLE MOMENTS
         opt%nam = system%num * opt%ns * system%nrmax !PARAMETER: MAX. DIMENSION OF VECTORS FOR  ANDERSON MIXING

         Return

!--------------------------------- ERROR MESSAGES
296      Write (IW6, 196)
196      Format (/ ' **** INPUT ERROR:  NSYM AND INVE ', 'ARE INCONSISTENT')
         Stop
297      Write (IW6, 197)
197      Format (/ ' **** INPUT ERROR:  NSYM, INVE ' / 3 X, '  AND 2D-TRANSLATION VECTORS ARE INCONSISTENT')
         Stop

      End Subroutine begi

!*******************
!XXX    TBK     ****
!*******************

      Subroutine tbk (opt, system, bz, harmon)
!***********************************************
!!$   BLOCH TRANSFORM OF TB-STRUCTURE CONSTANTS
!***********************************************
         Use str_const
         Use bzgrid
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Type (bzone) :: bz
         Type (sharm) :: harmon
         Real (Kind=prec) :: TRAN (2), arg, xefa, yefa, w1x, w2x, w1y, w2y, DEV, DEV1
         Real (Kind=prec), Parameter :: rcz = 0.0d0
         Integer :: i, ibz, ip, ib, ii, jj, ii0, jj0, jcl, jdir, jb, jtra1, jtra2, iq, jq, ipp1
         Complex (Kind=prec) :: zz, zw1, zw2, zefa
         Real (Kind=prec), Pointer :: str (:, :, :, :, :), bstr (:, :, :, :), vstr (:, :, :, :)
         Integer, Pointer :: jdpp (:, :, :), jbva (:, :, :), jtra (:, :, :, :), JDRPB (:, :), JBVAB (:, :), &
        & JTRAB (:, :, :), JDRPV (:, :), JBVAV (:, :), JTRAV (:, :, :), NSCL (:, :), NSCLB (:), NSCLV (:)
         zz = DCMPLX (rcz, rcz)

!!$     first we calculate structure constant for bulk, vaccum, and scattering region
!!$  scattering region
         Allocate (str(opt%nlsq, opt%nlsq, MNCL, system%NB, 0:opt%NPP1))
         str = 0.0d0
         Allocate (jdpp(MNCL, system%NB, 0:opt%NPP1))
         jdpp = 0
         Allocate (jbva(MNCL, system%NB, 0:opt%NPP1))
         jbva = 0
         Allocate (jtra(2, MNCL, system%NB, 0:opt%NPP1))
         jtra = 0
         Allocate (NSCL(system%NB, 0:opt%NPP1))
         NSCL = 0

         Call tbri (opt, system, harmon, str, jdpp, jbva, jtra, NSCL)!                      STRUCTURE CONSTANTS FOR SCATTERING REGION


!!$  bulk
         Allocate (bstr(opt%nlsq, opt%nlsq, MNCL, system%NB))
         bstr = 0.0d0
         Allocate (JBVAB(MNCL, system%NB))
         JBVAB = 0
         Allocate (JDRPB(MNCL, system%NB))
         JDRPB = 0
         Allocate (JTRAB(2, MNCL, system%NB))
         JTRAB = 0
         Allocate (NSCLB(system%NB))
         NSCLB = 0

         Call tbrb (opt, system, harmon, bstr, JBVAB, JDRPB, JTRAB, NSCLB)!                      STRUCTURE CONSTANTS FOR BULK


!!$  vacuum
         Allocate (vstr(opt%nlsq, opt%nlsq, MNCL, system%NB))
         vstr = 0.0d0
         Allocate (JBVAV(MNCL, system%NB))
         JBVAV = 0
         Allocate (JDRPV(MNCL, system%NB))
         JDRPV = 0
         Allocate (JTRAV(2, MNCL, system%NB))
         JTRAV = 0
         Allocate (NSCLV(system%NB))
         NSCLV = 0

         Call tbrv (opt, system, harmon, vstr, JBVAV, JDRPV, JTRAV, NSCLV)!                      STRUCTURE CONSTANTS FOR VACUUM

!!$                   BLOCH TRANSFORMATION OF STRUCTURE CONSTANTS
!!$  scattering region
         Write (IW6, 111)
111      Format (/ / 3 X, ' **********  BLOCH TRANSFORM ', 'OF TB-CONSTANTS  ********** ' / / 10 X, '    CHEC&
        &K OF HERMITICITY : ')

         Allocate (zski(opt%nblsq, opt%nblsq, system%np, bz%nkgrid))
         Allocate (zsko(opt%nblsq, opt%nblsq, system%np, bz%nkgrid))
         Allocate (zskf(opt%nblsq, opt%nblsq, system%np, bz%nkgrid))
         zski = zz
         zsko = zz
         zskf = zz

         Do ibz = 1, bz%nkgrid
            Do ip = 1, system%np
               Do ib = 1, system%NB
                  ii0 = (ib-1) * opt%nlsq
                  Do jcl = 1, NSCL (ib, ip)
                     jdir = jdpp (jcl, ib, ip)
                     jb = jbva (jcl, ib, ip)
                     jtra1 = jtra (1, jcl, ib, ip)
                     jtra2 = jtra (2, jcl, ib, ip)
                     jj0 = (jb-1) * opt%nlsq
                     Do i = 1, 2
                        TRAN (i) = jtra1 * system%vbr(i, 1) + jtra2 * system%vbr(i, 2)
                     End Do
                     arg = rcz
                     Do i = 1, 2
                        arg = arg + bz%k(i, ibz) * TRAN (i)
                     End Do
                     xefa = dCOS (arg)
                     yefa = dSIN (arg)
                     zefa = DCMPLX (xefa, yefa)
                     If (jdir == 0) Then
                        Do iq = 1, opt%nlsq
                           Do jq = 1, opt%nlsq
                              zski (jj0+jq, ii0+iq, ip, ibz) = zski (jj0+jq, ii0+iq, ip, ibz) + zefa * str &
                             & (jq, iq, jcl, ib, ip)
                           End Do
                        End Do
                     End If
                     If (jdir == 1) Then
                        Do iq = 1, opt%nlsq
                           Do jq = 1, opt%nlsq
                              zskf (jj0+jq, ii0+iq, ip, ibz) = zskf (jj0+jq, ii0+iq, ip, ibz) + zefa * str &
                             & (jq, iq, jcl, ib, ip)
                           End Do
                        End Do
                     End If
                     If (jdir ==-1) Then
                        Do iq = 1, opt%nlsq
                           Do jq = 1, opt%nlsq
                              zsko (jj0+jq, ii0+iq, ip, ibz) = zsko (jj0+jq, ii0+iq, ip, ibz) + zefa * str &
                             & (jq, iq, jcl, ib, ip)
                           End Do
                        End Do
                     End If
                  End Do
               End Do
            End Do
         End Do

!!$                       CHECK OF HERMITICITY
         Do ip = 1, system%np
            DEV = rcz
            Do ibz = 1, bz%nkgrid
               Do ii = 1, opt%nblsq
                  Do jj = ii, opt%nblsq
                     zw1 = zski (ii, jj, ip, ibz)
                     w1x = dble (zw1)
                     w1y = DIMAG (zw1)
                     zw2 = zski (jj, ii, ip, ibz)
                     w2x = dble (zw2)
                     w2y = DIMAG (zw2)
                     DEV1 = Abs (w1x-w2x) + Abs (w1y+w2y)
                     DEV = Max (DEV, DEV1)
                  End Do
               End Do
            End Do
            Write (IW6, 120) ip, DEV
120         Format (/ 4 X, '  LAYER:  IP=', I4, 6 X, 'DEVIATION=', G12.4)
         End Do

         If (system%np .Gt. 1) Then
            Do ip = 1, system%np - 1
               ipp1 = ip + 1
               DEV = rcz
               Do ibz = 1, bz%nkgrid
                  Do ii = 1, opt%nblsq
                     Do jj = 1, opt%nblsq
                        zw1 = zskf (ii, jj, ip, ibz)
                        w1x = dble (zw1)
                        w1y = DIMAG (zw1)
                        zw2 = zsko (jj, ii, ipp1, ibz)
                        w2x = dble (zw2)
                        w2y = DIMAG (zw2)
                        DEV1 = dabs (w1x-w2x) + dabs (w1y+w2y)
                        DEV = Max (DEV, DEV1)
                     End Do
                  End Do
               End Do
               Write (IW6, 125) ip, ipp1, DEV
125            Format (/ 4 X, 'LAYERS:  IP=', I4, '  JP=', I4, 6 X, 'DEVIATION=', G12.4)
            End Do
         End If

!----------------------------------- BULK REGION
         Allocate (zbski(opt%nblsq, opt%nblsq, bz%nkgrid))
         Allocate (zbsko(opt%nblsq, opt%nblsq, bz%nkgrid))
         Allocate (zbskf(opt%nblsq, opt%nblsq, bz%nkgrid))
         zbski = zz
         zbsko = zz
         zbskf = zz

         Do ibz = 1, bz%nkgrid
            Do ib = 1, system%NB
               ii0 = (ib-1) * opt%nlsq
               Do jcl = 1, NSCLB (ib)
                  jdir = JDRPB (jcl, ib)
                  jb = JBVAB (jcl, ib)
                  jtra1 = JTRAB (1, jcl, ib)
                  jtra2 = JTRAB (2, jcl, ib)
                  jj0 = (jb-1) * opt%nlsq
                  Do i = 1, 2
                     TRAN (i) = jtra1 * system%vbr(i, 1) + jtra2 * system%vbr(i, 2)
                  End Do
                  arg = rcz
                  Do i = 1, 2
                     arg = arg + bz%k(i, ibz) * TRAN (i)
                  End Do
                  xefa = Cos (arg)
                  yefa = Sin (arg)
                  zefa = DCMPLX (xefa, yefa)
                  If (jdir == 0) Then
                     Do iq = 1, opt%nlsq
                        Do jq = 1, opt%nlsq
                           zbski (jj0+jq, ii0+iq, ibz) = zbski (jj0+jq, ii0+iq, ibz) + zefa * bstr (jq, iq, &
                          & jcl, ib)
                        End Do
                     End Do
                  End If
                  If (jdir == 1) Then
                     Do iq = 1, opt%nlsq
                        Do jq = 1, opt%nlsq
                           zbskf (jj0+jq, ii0+iq, ibz) = zbskf (jj0+jq, ii0+iq, ibz) + zefa * bstr (jq, iq, &
                          & jcl, ib)
                        End Do
                     End Do
                  End If
                  If (jdir .Eq.-1) Then
                     Do iq = 1, opt%nlsq
                        Do jq = 1, opt%nlsq
                           zbsko (jj0+jq, ii0+iq, ibz) = zbsko (jj0+jq, ii0+iq, ibz) + zefa * bstr (jq, iq, &
                          & jcl, ib)
                        End Do
                     End Do
                  End If
               End Do
            End Do
         End Do

!!$                       CHECK OF HERMITICITY
         DEV = rcz
         Do ibz = 1, bz%nkgrid
            Do ii = 1, opt%nblsq
               Do jj = ii, opt%nblsq
                  zw1 = zbski (ii, jj, ibz)
                  w1x = dble (zw1)
                  w1y = DIMAG (zw1)
                  zw2 = zbski (jj, ii, ibz)
                  w2x = dble (zw2)
                  w2y = DIMAG (zw2)
                  DEV1 = Abs (w1x-w2x) + Abs (w1y+w2y)
                  DEV = Max (DEV, DEV1)
               End Do
            End Do
         End Do

         Write (IW6, 140) DEV
140      Format (/ 4 X, ' LEFT BULK - IN-LAYER', 6 X, 'DEVIATION=', G12.4)

         DEV = rcz
         Do ibz = 1, bz%nkgrid
            Do ii = 1, opt%nblsq
               Do jj = 1, opt%nblsq
                  zw1 = zbskf (ii, jj, ibz)
                  w1x = dble (zw1)
                  w1y = DIMAG (zw1)
                  zw2 = zbsko (jj, ii, ibz)
                  w2x = dble (zw2)
                  w2y = DIMAG (zw2)
                  DEV1 = Abs (w1x-w2x) + Abs (w1y+w2y)
                  DEV = Max (DEV, DEV1)
               End Do
            End Do
         End Do
         Write (IW6, 145) DEV
145      Format (/ 4 X, ' LEFT  BULK - OFF-LAYER', 5 X, 'DEVIATION=', G12.4)

!----------------------------------- VACUUM REGION
         If (opt%lident == 0) Then ! if leads are identical than we use the same str.const for bulk and vacuum
            Allocate (zvski(opt%nblsq, opt%nblsq, bz%nkgrid))
            Allocate (zvsko(opt%nblsq, opt%nblsq, bz%nkgrid))
            Allocate (zvskf(opt%nblsq, opt%nblsq, bz%nkgrid))
            zvski = zz
            zvsko = zz
            zvskf = zz

            Do ibz = 1, bz%nkgrid
               Do ib = 1, system%NB
                  ii0 = (ib-1) * opt%nlsq
                  Do jcl = 1, NSCLV (ib)
                     jdir = JDRPV (jcl, ib)
                     jb = JBVAV (jcl, ib)
                     jtra1 = JTRAV (1, jcl, ib)
                     jtra2 = JTRAV (2, jcl, ib)
                     jj0 = (jb-1) * opt%nlsq
                     Do i = 1, 2
                        TRAN (i) = jtra1 * system%vbr(i, 1) + jtra2 * system%vbr(i, 2)
                     End Do
                     arg = rcz
                     Do i = 1, 2
                        arg = arg + bz%k(i, ibz) * TRAN (i)
                     End Do
                     xefa = dCOS (arg)
                     yefa = dSIN (arg)
                     zefa = DCMPLX (xefa, yefa)
                     If (jdir == 0) Then
                        Do iq = 1, opt%nlsq
                           Do jq = 1, opt%nlsq
                              zvski (jj0+jq, ii0+iq, ibz) = zvski (jj0+jq, ii0+iq, ibz) + zefa * vstr (jq, &
                             & iq, jcl, ib)
                           End Do
                        End Do
                     End If
                     If (jdir == 1) Then
                        Do iq = 1, opt%nlsq
                           Do jq = 1, opt%nlsq
                              zvskf (jj0+jq, ii0+iq, ibz) = zvskf (jj0+jq, ii0+iq, ibz) + zefa * vstr (jq, &
                             & iq, jcl, ib)
                           End Do
                        End Do
                     End If
                     If (jdir .Eq.-1) Then
                        Do iq = 1, opt%nlsq
                           Do jq = 1, opt%nlsq
                              zvsko (jj0+jq, ii0+iq, ibz) = zvsko (jj0+jq, ii0+iq, ibz) + zefa * vstr (jq, &
                             & iq, jcl, ib)
                           End Do
                        End Do
                     End If
                  End Do
               End Do
            End Do

!!$                       CHECK OF HERMITICITY
            DEV = rcz
            Do ibz = 1, bz%nkgrid
               Do ii = 1, opt%nblsq
                  Do jj = ii, opt%nblsq
                     zw1 = zvski (ii, jj, ibz)
                     w1x = dble (zw1)
                     w1y = DIMAG (zw1)
                     zw2 = zvski (jj, ii, ibz)
                     w2x = dble (zw2)
                     w2y = DIMAG (zw2)
                     DEV1 = Abs (w1x-w2x) + Abs (w1y+w2y)
                     DEV = Max (DEV, DEV1)
                  End Do
               End Do
            End Do

            Write (IW6, 160) DEV
160         Format (/ 4 X, 'RIGHT BULK(VACUUM) - IN-LAYER', 6 X, 'DEVIATION=', G12.4)

            DEV = rcz
            Do ibz = 1, bz%nkgrid
               Do ii = 1, opt%nblsq
                  Do jj = 1, opt%nblsq
                     zw1 = zvskf (ii, jj, ibz)
                     w1x = dble (zw1)
                     w1y = DIMAG (zw1)
                     zw2 = zvsko (jj, ii, ibz)
                     w2x = dble (zw2)
                     w2y = DIMAG (zw2)
                     DEV1 = Abs (w1x-w2x) + Abs (w1y+w2y)
                     DEV = Max (DEV, DEV1)
                  End Do
               End Do
            End Do

            Write (IW6, 165) DEV
165         Format (/ 4 X, 'RIGHT BULK(VACUUM) - OFF-LAYER', 5 X, 'DEVIATION=', G12.4)
         Else
            Write (IW6, 146)
146         Format (/ 4 X, 'LEADS ARE IDENTICAL- THUS DEVIATIONS FOR RIGHT AND LEFT BULK ARE THE SAME ')
         End If

!!$     DEALLOCATION BLOCK FOR ALL THE GARBAGE WE WILL NOT NEED
         Deallocate (str)
         Deallocate (jdpp)
         Deallocate (jbva)
         Deallocate (jtra)
         Deallocate (NSCL)
         Deallocate (bstr)
         Deallocate (JDRPB)
         Deallocate (JBVAB)
         Deallocate (JTRAB)
         Deallocate (NSCLB)
         Deallocate (vstr)
         Deallocate (JDRPV)
         Deallocate (JBVAV)
         Deallocate (JTRAV)
         Deallocate (NSCLV)
!!$     END OF DEALLOCATION BLOCK

!!$       open(15,file='zski',status='old')
!!$       do ii = 1, opt%nblsq
!!$        do jj = 1, opt%nblsq
!!$         do ip = 1, system%np
!!$          do ibz = 1, bz%nkgrid
!!$           read(15,*) zski(ii,jj,ip,ibz)
!!$          enddo
!!$         enddo
!!$        enddo
!!$       enddo
!!$       close(15)

!!$       open(15,file='zskf',status='old')
!!$       do ii = 1, opt%nblsq
!!$        do jj = 1, opt%nblsq
!!$         do ip = 1, system%np
!!$          do ibz = 1, bz%nkgrid
!!$           read(15,*) zskf(ii,jj,ip,ibz)
!!$          enddo
!!$         enddo
!!$        enddo
!!$       enddo
!!$       close(15)

!!$       open(15,file='zsko',status='old')
!!$       do ii = 1, opt%nblsq
!!$        do jj = 1, opt%nblsq
!!$         do ip = 1, system%np
!!$          do ibz = 1, bz%nkgrid
!!$           read(15,*) zsko(ii,jj,ip,ibz)
!!$          enddo
!!$         enddo
!!$        enddo
!!$       enddo
!!$       close(15)

         Return
      End Subroutine tbk


!*******************
!XXX    TBRI    ****
!*******************
      Subroutine tbri (opt, system, harmon, str, jdpp, jbva, jtra, NSCL)
!----------------------------------------------------
!!$   TB-STRUCTURE CONSTANTS FOR INTERMEDIATE REGION
!----------------------------------------------------
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Type (sharm) :: harmon
         Real (Kind=prec) :: POW (3, system%NB,-1:opt%NPP2), POL (3, MNCL), WSW (system%NB,-1:opt%NPP2), WSL &
        & (MNCL)
         Real (Kind=prec) :: POLT (3), DPOL (3), WWMA (opt%nlsq, opt%nlsq), WWMB (opt%nlsq, opt%nlsq)
         Real (Kind=prec), Parameter :: rcz = 0.0d0, SMALL = 0.01d0
         Integer :: ib, i, ip, jdir, jtra1, jtra2, jcl, iq, jq, jp, icl, isum, jb, ncl
         Real (Kind=prec) :: DMINSQ, CUTRSQ, DMAXSQ, sum, DEV, DEV1, ssc (opt%nlsq, opt%nlsq, MNCL)
         Real (Kind=prec), Pointer :: str (:, :, :, :, :)
         Integer, Pointer :: jdpp (:, :, :), jbva (:, :, :), jtra (:, :, :, :), NSCL (:, :)

         Write (IW6, 111)
111      Format (/ / ' *********  TB-CONSTANTS FOR', ' THE INTERMEDIATE REGION  *********')

!-------------- POSITIONS AND LOCAL AVERAGE WS-RADII INCLUDING TWO  NEIGHBORING LAYERS ON BOTH
!!$                                SIDES OF INERMEDIATE REGION

         Do ib = 1, system%NB
            Do i = 1, 3
               POW (i, ib,-1) = system%vpos(i, ib) + system%VTRV(i)
            End Do
            WSW (ib,-1) = system%avws
         End Do

         Do ib = 1, system%NB
            Do i = 1, 3
               POW (i, ib, 0) = system%vpos(i, ib)
            End Do
            WSW (ib, 0) = system%avws
         End Do

         Do ip = 1, system%np
            Do ib = 1, system%NB
               Do i = 1, 3
                  POW (i, ib, ip) = system%POS(i, ib, ip)
               End Do
               WSW (ib, ip) = system%aws(ib, ip)
            End Do
         End Do

         Do ib = 1, system%NB
            Do i = 1, 3
               POW (i, ib, opt%NPP1) = system%BPOS(i, ib)
            End Do
            WSW (ib, opt%NPP1) = system%abws
         End Do

         Do ib = 1, system%NB
            Do i = 1, 3
               POW (i, ib, opt%NPP2) = system%BPOS(i, ib) + system%BTRV(i)
            End Do
            WSW (ib, opt%NPP2) = system%abws
         End Do
         DMINSQ = SMALL ** 2 * system%abws * system%avws
         CUTRSQ = system%cutrat ** 2

!--------------------- LOOP OVER CENTRAL ATOMS FOR SMALL CLUSTERS
         Do ip = 0, opt%NPP1
            Do ib = 1, system%NB
               jcl = 1 !------------------------- SELECTION OF THE SMALL CLUSTER
               Do i = 1, 3
                  POL (i, 1) = POW (i, ib, ip)
               End Do
               jdpp (1, ib, ip) = 0
               jbva (1, ib, ip) = ib
               jtra (1, 1, ib, ip) = 0
               jtra (2, 1, ib, ip) = 0
               WSL (1) = WSW (ib, ip)

               Do jdir = - 1, 1
                  jp = ip + jdir
                  Do jb = 1, system%NB
                     Do jtra1 = - system%nmtr, system%nmtr
                        Do jtra2 = - system%nmtr, system%nmtr
                           Do i = 1, 2
                              POLT (i) = POW (i, jb, jp) + jtra1 * system%vbr(i, 1) + jtra2 * system%vbr(i, &
                             & 2)
                           End Do
                           POLT (3) = POW (3, jb, jp)
                           Do i = 1, 3
                              DPOL (i) = POLT (i) - POL (i, 1)
                           End Do
                           DMAXSQ = CUTRSQ * WSW (ib, ip) * WSW (jb, jp)
!c$$$             DMAXSQ=CUTRSQ*AVWS*ABWS
                           sum = DPOL (1) ** 2 + DPOL (2) ** 2 + DPOL (3) ** 2
                           If (sum .Lt. DMINSQ .Or. sum .Gt. DMAXSQ) Cycle
                           jcl = jcl + 1
                           If (jcl .Gt. MNCL) Go To 291
                           Do i = 1, 3
                              POL (i, jcl) = POLT (i)
                           End Do
                           jdpp (jcl, ib, ip) = jdir
                           jbva (jcl, ib, ip) = jb
                           jtra (1, jcl, ib, ip) = jtra1
                           jtra (2, jcl, ib, ip) = jtra2
                           WSL (jcl) = WSW (jb, jp)
                        End Do
                     End Do
                  End Do
               End Do
               ncl = jcl
               NSCL (ib, ip) = ncl
               Write (IW6, 115) ip, ib, ncl
115            Format (/ '    SITE:  IP=', I4, '  IB=', I3, '      CLUSTER SIZE= ', I3)

!---------- TB-CONSTANTS FOR THE SMALL CLUSTER AND THEIR STORAGE
               Call TBCL (ncl, POL, WSL, ssc, opt%nl, opt%nlsq, opt%mtbcl, system%qscr, harmon)

               Do jcl = 1, ncl
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        str (iq, jq, jcl, ib, ip) = ssc (iq, jq, jcl)
                     End Do
                  End Do
               End Do

            End Do
         End Do


!---------------------------- SYMMETRIZATION
!!$                -------- IN-LAYER CONSTANTS
         Do ip = 1, system%np
            jp = ip
            DEV = rcz
            Do ib = 1, system%NB
               Do jcl = 1, NSCL (ib, ip)
                  If (jdpp(jcl, ib, ip) .Ne. 0) Cycle
                  jb = jbva (jcl, ib, ip)
                  Do icl = 1, NSCL (jb, jp)
                     If (jdpp(icl, jb, jp) .Ne. 0) Cycle
                     If (jbva(icl, jb, jp) .Ne. ib) Cycle
                     isum = jtra (1, jcl, ib, ip) + jtra (1, icl, jb, jp)
                     If (isum .Ne. 0) Cycle
                     isum = jtra (2, jcl, ib, ip) + jtra (2, icl, jb, jp)
                     If (isum .Ne. 0) Cycle
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           WWMA (iq, jq) = str (iq, jq, jcl, ib, ip)
                        End Do
                     End Do
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           WWMB (iq, jq) = str (iq, jq, icl, jb, jp)
                        End Do
                     End Do
                     Call MASY (WWMA, WWMB, DEV1, opt%nlsq)
                     DEV = Max (DEV, DEV1)

                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           str (iq, jq, jcl, ib, ip) = WWMA (iq, jq)
                        End Do
                     End Do
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           str (iq, jq, icl, jb, jp) = WWMB (iq, jq)
                        End Do
                     End Do
                  End Do
               End Do
            End Do

            Write (IW6, 125) ip, DEV
125         Format (/ 4 X, '  LAYER:  IP=', I4, 6 X, 'DEVIATION=', G12.4)

         End Do
!!$                       -------- OFF-LAYER CONSTANTS
         Do ip = 0, system%np
            jp = ip + 1
            DEV = rcz
            Do ib = 1, system%NB
               Do jcl = 1, NSCL (ib, ip)
                  If (jdpp(jcl, ib, ip) .Ne. 1) Cycle
                  jb = jbva (jcl, ib, ip)
                  Do icl = 1, NSCL (jb, jp)
                     If (jdpp(icl, jb, jp) .Ne.-1) Cycle
                     If (jbva(icl, jb, jp) .Ne. ib) Cycle
                     isum = jtra (1, jcl, ib, ip) + jtra (1, icl, jb, jp)
                     If (isum .Ne. 0) Cycle
                     isum = jtra (2, jcl, ib, ip) + jtra (2, icl, jb, jp)
                     If (isum .Ne. 0) Cycle
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           WWMA (iq, jq) = str (iq, jq, jcl, ib, ip)
                        End Do
                     End Do
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           WWMB (iq, jq) = str (iq, jq, icl, jb, jp)
                        End Do
                     End Do
                     Call MASY (WWMA, WWMB, DEV1, opt%nlsq)
                     DEV = Max (DEV, DEV1)
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           str (iq, jq, jcl, ib, ip) = WWMA (iq, jq)
                        End Do
                     End Do
                     Do jq = 1, opt%nlsq
                        Do iq = 1, opt%nlsq
                           str (iq, jq, icl, jb, jp) = WWMB (iq, jq)
                        End Do
                     End Do
                  End Do
               End Do
            End Do

            Write (IW6, 135) ip, jp, DEV
135         Format (/ 4 X, 'LAYERS:  IP=', I4, '   JP=', I4, 4 X, 'DEVIATION=', G12.4)

         End Do

!------------ PRINT OF DIAGONAL ELEMENTS
         Write (IW6, 150)
150      Format (/ 5 X, '***  DIAGONAL ELEMENTS OF ', 'TB-STRUCTURE CONSTANTS :')
151      Format (/ '      SITE:  IP=', I4, '  IB=', I3)
104      Format (1 X, 4 G15.7)
         Do ip = 0, opt%NPP1
            Do ib = 1, system%NB
               Write (IW6, 151) ip, ib
               Write (IW6, 104) (str(iq, iq, 1, ib, ip), iq=1, opt%nlsq)
            End Do
         End Do

         Return

291      Write (IW6, 191) jcl
191      Format (/ ' **** ERROR IN TBRI :' / 10 X, ' JCL GREATER THAN mncl,  JCL=', I3)

         Stop
      End Subroutine tbri

!*******************
!XXX    MASY    ****
!*******************
      Subroutine MASY (A, B, DEV, nlsq)
         Implicit None
         Integer :: i, j, nlsq
         Real (Kind=prec) :: DEV, eps, val, A (nlsq, nlsq), B (nlsq, nlsq)
         Real (Kind=prec), Parameter :: rcz = 0.0d0, rch = 0.5d0
!------------------------------------------------
!!$   SPECIAL SYMMETRIZATION OF MATRICES A AND B :
!!$  - CALCULATES DEVIATION OF A(TRANSPOSE) AND B
!!$  - CORRECTS A AND B TO SATISFY A(TRANSPOSE) = B
!------------------------------------------------
         DEV = rcz

         Do i = 1, nlsq
            Do j = 1, nlsq
               eps = dabs (A(i, j)-B(j, i))
               DEV = Max (DEV, eps)
               val = rch * (A(i, j)+B(j, i))
               A (i, j) = val
               B (j, i) = val
            End Do
         End Do
         Return
      End Subroutine MASY

!*******************
!XXX    TBCL    ****
!*******************
      Subroutine TBCL (ncl, POL, WSA, ssc, nl, nlsq, mtbcl, qscr, harmon)
!--------------------------------------------------------
!!$     TB-STRUCTURE CONSTANTS FOR A SMALL CLUSTER
!!$     CALCULATED BY MATRIX INVERSION IN real SPACE
!--------------------------------------------------------
!!$   INPUT:  NL - NUMBER OF PARTIAL WAVES,
!!$           ncl - SIZE OF THE SMALL CLUSTER,
!!$           POL(K,ICL) (K=1,2,3, ICL=1,2,...,ncl) -
!!$            - COORDINATES OF THE ATOMS OF THE SMALL CLUSTER,
!!$           ICL=1 CORRESPONDS TO THE CENTRAL ATOM OF THE CLUSTER,
!!$           WSA(ICL) (ICL=1,2,...,ncl) - LOCAL AVERAGE WS-RADII.
!!$  OUTPUT:  ssc(i,J,ICL) (i,J - ORBITAL INDICES,
!!$                        ICL=1,2,...,ncl) -
!!$           - TB-STRUCTURE CONSTANTS S(i,ICL,J,1) CONNECTING
!!$           THE ICL-TH ATOM WITH THE FIRST (CENTRAL) ATOM
!--------------------------------------------------------
!!$   REMARK: THE GAUNT FACTORS MUST BE CALCULATED BEFORE !
!!$           NCALLH MUST BE SET ZERO BEFORE !
!---------------------------------------------------------

         Implicit None
         Integer :: nl, nlsq, mtbcl
         Real (Kind=prec) :: POL (3, ncl), WSA (ncl), RI (3), RJ (3), QR (nlsq), aa (mtbcl, mtbcl), BB &
        & (mtbcl, nlsq), WW (mtbcl), SC (nlsq, nlsq), SC1 (nlsq, nlsq), ssc (nlsq, nlsq, MNCL), qscr (nl)
         Real (Kind=prec), Parameter :: rc1 = 1.0d0
         Integer :: INDW (nlsq, ncl), il, ista, ifin, ii, jj, icl, i, j, jcl, k, nn, ncl, lmax
         Type (sharm) :: harmon
#if defined(DO_SC_SYM_TUREK)
         Integer, Pointer :: symfield (:)
         Integer :: NSYMs, srows
#endif

         lmax = 2 * (nl-1)
         aa = 0.0d0
         BB = 0.0d0
         QR = 0.0d0
         Do il = 1, nl
            ista = (il-1) ** 2 + 1
            ifin = il ** 2
            Do i = ista, ifin
               QR (i) = rc1 / qscr (il)
            End Do
         End Do

         ii = 0
         Do icl = 1, ncl
            Do i = 1, nlsq
               ii = ii + 1
               INDW (i, icl) = ii
            End Do
         End Do
         nn = ii

#if defined(DO_SC_SYM_TUREK)
         srows = nlsq
         NSYMs = srows * (nn-srows)
         Allocate (symfield(2*NSYMs))
#endif

!!$    FORMATION OF THE MATRIX (QR - S) (=AA)  AND OF THE RIGHT-HAND SIDES (=BB)

         Do jcl = 1, ncl
            Do k = 1, 3
               RJ (k) = POL (k, jcl)
            End Do
            Do icl = jcl, ncl
               Do k = 1, 3
                  RI (k) = POL (k, icl)
               End Do
               Call CANSC_L (WSA(icl), WSA(jcl), RI, RJ, SC, nl, lmax, nlsq, harmon)
               Call CANSC_L (WSA(jcl), WSA(icl), RJ, RI, SC1, nl, lmax, nlsq, harmon)
           
               SC = 0.5d0 * (SC+transpose(SC1))
#if defined(DO_SC_SYM_TUREK)
               Where (Abs(SC) < 1.0d-10)
                  SC = 0.0d0
               End Where
#endif               
               Do j = 1, nlsq
                  jj = INDW (j, jcl)
                  Do i = 1, nlsq
                     ii = INDW (i, icl)
                     aa (ii, jj) = - SC (i, j)
                     If (jcl == 1) BB (ii, j) = SC (i, j)
                  End Do
               End Do
            End Do
         End Do


#if defined(DO_SC_SYM_TUREK)
         Call getscsym (NSYMs, aa((srows+1) :nn, 1:srows), symfield, 1.0d-10)
#endif

!!$          Write (*,*) 'iface-pre--'
!!$          Write (*, '(9(1xg22.14))') (aa(i, 1:srows), i=1, nn)
!!$          Write (*,*) '-----'


         Do icl = 1, ncl
            Do i = 1, nlsq
               ii = INDW (i, icl)
               aa (ii, ii) = QR (i)
            End Do
         End Do

!!$                      SOLUTION OF THE LINEAR EQUATIONS
         Call SOPO (aa, mtbcl, nn, BB, mtbcl, nlsq, WW)
!!$                      THE TB-STRUCTURE CONSTANTS


         Do j = 1, nlsq
            Do icl = 1, ncl
               Do i = 1, nlsq
                  ii = INDW (i, icl)
                  BB (ii, j) = QR (i) * BB (ii, j)
               End Do
            End Do
         End Do

#if defined(DO_SC_SYM_TUREK)
         Call setscsym (NSYMs, BB((srows+1) :nn, 1:srows), symfield)
         Call symscons (srows*srows, BB(1:srows, 1:srows), 1.0d-5)
         Where (Abs(BB) < 1.0d-8)
            BB = 0.0d0
         End Where
         Deallocate (symfield)
#endif
!!$          Write (*,*) 'iface--'
!!$          Write (*, '(9(1xg22.14))') (BB(i, 1:srows), i=1, nn)
!!$          Write (*,*) '-----'

         Do j = 1, nlsq
            Do icl = 1, ncl
               Do i = 1, nlsq
                  ii = INDW (i, icl)
                  ssc (i, j, icl) = BB (ii, j)
               End Do
            End Do
         End Do


         Return
      End Subroutine TBCL

!*******************
!XXX   CANSC    ****
!*******************
      Subroutine CANSC_L (W1, W2, R1, R2, SC, nl, lmax, nlsq, harmon)
         Implicit None

!---------------------------------------------------------------
!!$   CANONICAL STRUCTURE CONSTANT CONNECTING 2 POINTS R1 AND R2
!!$       NL - NUMBER OF PARTIAL WAVES
!!$       W1,W2 -  LOCAL AVERAGE W.S.-RADII
!!$       R1(i),R2(i) - THE COORDINATES OF THE 2 POINTS
!!$       SC(I1,I2) - THE STRUCTURE CONSTANT MATRIX
!---------------------------------------------------------------
!!$   REMARK: THE GAUNT FACTORS MUST BE CALCULATED BEFORE !
!!$           NCALLH MUST BE SET ZERO BEFORE !
!---------------------------------------------------------------

         Type (sharm) :: harmon
         Real (Kind=prec) :: DFAC (0:lmax), AM1LP1 (0:nl), POW1 (0:nl), POW2 (0:nl), PREF (0:nl, 0:nl)
         Real (Kind=prec), Pointer :: yps (:)
         Real (Kind=prec) :: SC (nlsq, nlsq), dx, dy, dz, sum
         Real (Kind=prec), Parameter :: rcz = 0.0d0, rc1 = 1.0d0, rc2 = 2.0d0, rc4 = 4.0d0, RMIN = 0.01d0
         Integer :: i, i1, I2, lmax1, lmaxh, l, l2, l1, m2, m1, ista, ifin, nl, lmax, nlsq
         Real (Kind=prec) :: tx, ty, tz, r, rbw, asig, pi, pi8, rat1, W1, W2, R1 (3), R2 (3), RAT2, TWOLM1, &
        & CIT, AJM

         SC (:, :) = rcz

         tx = R2 (1) - R1 (1)
         ty = R2 (2) - R1 (2)
         tz = R2 (3) - R1 (3)
         r = dsqrt (tx**2+ty**2+tz**2)
         rbw = r / (W1+W2)
         If (rbw .Lt. RMIN) Return



         lmax1 = nl - 1
         lmaxh = 2 * lmax1

         DFAC (0) = 1.0d0
         Do l = 1, lmaxh
            TWOLM1 = dble (2*l-1)
            DFAC (l) = TWOLM1 * DFAC (l-1)
         End Do

         asig = rc1
         Do l = 0, lmax1
            asig = - asig
            AM1LP1 (l) = asig
         End Do

         pi = rc4 * datan (rc1)
         pi8 = rc2 * rc4 * pi

         rat1 = W1 / r
         POW1 (0) = dsqrt (rat1)
         Do l = 1, lmax1
            POW1 (l) = rat1 * POW1 (l-1)
         End Do

         RAT2 = W2 / r
         POW2 (0) = dsqrt (RAT2)
         Do l = 1, lmax1
            POW2 (l) = RAT2 * POW2 (l-1)
         End Do

         Do l2 = 0, lmax1
            Do l1 = 0, lmax1
               l = l1 + l2
               CIT = AM1LP1 (l2) * pi8 * DFAC (l) * POW1 (l1) * POW2 (l2)
               AJM = DFAC (l1) * DFAC (l2)
               PREF (l1, l2) = CIT / AJM
            End Do
         End Do

         dx = tx / r
         dy = ty / r
         dz = tz / r

         yps => harm (lmaxh, dx, dy, dz, harmon%cnh)

         I2 = 0
         Do l2 = 0, lmax1
            Do m2 = - l2, l2
               I2 = I2 + 1
               i1 = 0
               Do l1 = 0, lmax1
                  Do m1 = - l1, l1
                     i1 = i1 + 1
                     l = l1 + l2
                     ista = l ** 2 + 1
                     ifin = (l+1) ** 2
                     sum = rcz
                     Do i = ista, ifin
                        sum = sum + harmon%gfrh(i, i1, I2) * yps (i)
                     End Do
                     SC (i1, I2) = PREF (l1, l2) * sum
                  End Do
               End Do
            End Do
         End Do
         Deallocate (yps)
         Return
      End Subroutine CANSC_L

!*******************
!XXX    SOPO    ****
!*******************
      Subroutine SOPO (A, LDA, N, B, LDB, M, W)
         Implicit None

!---------------------------------------------------------------
!!$      SOLUTION OF A SET OF LINEAR EQUATIONS WITH A POSITIVE
!!$      DEFINITE real SYMMETRIC MATRIX AND WITH MULTIPLE R.H.S.
!---------------------------------------------------------------
!!$     NA VSTUPU :
!!$        N ... DIMENZE MATICE
!!$        M ... POCET PRAVYCH STRAN
!!$        A(i,J) ... MATICE SOUSTAVY - JEN DOLNI TROJUHELNIK JE
!!$                   TREBA   (TJ.  i = 1, ... N,   J = 1, ... i)
!!$        B(i,K) ... PRAVE STRANY (K = 1, ... M,    i = 1, ... N)
!!$        W(i)  ...  PRACOVNI POLE DELKY N
!!$     NA VYSTUPU : A(i,J) JE PREPSANA SVYM CHOLESKYHO FAKTOREM
!!$                  B(i,K) .... PRISLUSNA RESENI

         Integer :: i, j, l, M, N, LDA, LDB, mq
         Real (Kind=prec), Parameter :: rc1 = 1.0d0
         Real (Kind=prec) :: A (LDA, N), B (LDB, M), W (N), DUM
!!$                                   ROZKLAD
         Do l = 1, N - 1
            A (l, l) = dsqrt (A(l, l))
            DUM = rc1 / A (l, l)
            Do i = l + 1, N
               W (i) = A (i, l) * DUM
               A (i, l) = W (i)
            End Do
            Do j = l + 1, N
               Do i = j, N
                  A (i, j) = A (i, j) - W (i) * W (j)
               End Do
            End Do
         End Do
         A (N, N) = dsqrt (A(N, N))

!!$                                 CYKLUS PRES PRAVE STRANY
         Do mq = 1, M
            Do l = 1, N - 1 ! INVERZE DOLNI TROJUH. MATICE
               W (l) = B (l, mq) / A (l, l)
               Do i = l + 1, N
                  B (i, mq) = B (i, mq) - A (i, l) * W (l)
               End Do
            End Do
            W (N) = B (N, mq) / A (N, N)
!!$                                 INVERZE HORNI TROJUH. MATICE
            B (N, mq) = W (N) / A (N, N)
            Do l = N - 1, 1, - 1
               Do i = N, l + 1, - 1
                  W (l) = W (l) - B (i, mq) * A (i, l)
               End Do
               B (l, mq) = W (l) / A (l, l)
            End Do
         End Do
         Return
      End Subroutine SOPO

!*******************
!XXX    TBRB    ****
!*******************
      Subroutine tbrb (opt, system, harmon, bstr, JBVAB, JDRPB, JTRAB, NSCLB)
         Implicit None
!----------------------------------------------------
!!$   TB-STRUCTURE CONSTANTS FOR SEMiiNFINITE BULK
!----------------------------------------------------
         Type (params) :: opt
         Type (atoms_set) :: system
         Type (sharm) :: harmon
         Real (Kind=prec), Parameter :: rcz = 0.0d0, SMALL = 0.01d0
         Real (Kind=prec) :: POL (3, MNCL), WSL (MNCL), POLT (3), DPOL (3), WWMA (opt%nlsq, opt%nlsq), WWMB &
        & (opt%nlsq, opt%nlsq)
         Real (Kind=prec) :: DMAXSQ, DMINSQ, sum, DEV, DEV1, ssc (opt%nlsq, opt%nlsq, MNCL)
         Integer :: ib, jcl, jdir, jb, jtra1, jtra2, i, jq, iq, ncl, icl, isum
         Real (Kind=prec), Pointer :: bstr (:, :, :, :)
         Integer, Pointer :: JDRPB (:, :), JTRAB (:, :, :), JBVAB (:, :), NSCLB (:)


         ssc = 0.0d0

         Write (IW6, 111)
111      Format (/ / ' ************  TB-CONSTANTS FOR', ' THE BULK REGION  ************')

!------------------------------- CUT-OFF DISTANCE
         DMAXSQ = (system%cutrat*system%abws) ** 2
         DMINSQ = (SMALL*system%abws) ** 2


         Do ib = 1, system%NB !------------------ LOOP OVER CENTRAL ATOMS FOR SMALL CLUSTERS
            jcl = 1 !-------------------- SELECTION OF THE SMALL CLUSTER
            Do i = 1, 3
               POL (i, 1) = system%BPOS(i, ib)
            End Do
            JDRPB (1, ib) = 0
            JBVAB (1, ib) = ib
            JTRAB (1, 1, ib) = 0
            JTRAB (2, 1, ib) = 0
            WSL (1) = system%abws
            Do jdir = - 1, 1
               Do jb = 1, system%NB
                  Do jtra1 = - system%nmtr, system%nmtr
                     Do jtra2 = - system%nmtr, system%nmtr
                        Do i = 1, 2
                           POLT (i) = system%BPOS(i, jb) + jdir * system%BTRV(i) + jtra1 * system%vbr(i, 1) + &
                          & jtra2 * system%vbr(i, 2)
                        End Do
                        POLT (3) = system%BPOS(3, jb) + jdir * system%BTRV(3)
                        Do i = 1, 3
                           DPOL (i) = POLT (i) - POL (i, 1)
                        End Do
                        sum = DPOL (1) ** 2 + DPOL (2) ** 2 + DPOL (3) ** 2
                        If (sum .Lt. DMINSQ .Or. sum .Gt. DMAXSQ) Cycle
                        jcl = jcl + 1
                        If (jcl .Gt. MNCL) Go To 291
                        Do i = 1, 3
                           POL (i, jcl) = POLT (i)
                        End Do
                        JDRPB (jcl, ib) = jdir
                        JBVAB (jcl, ib) = jb
                        JTRAB (1, jcl, ib) = jtra1
                        JTRAB (2, jcl, ib) = jtra2
                        WSL (jcl) = system%abws
                     End Do
                  End Do
               End Do
            End Do

            ncl = jcl
            NSCLB (ib) = ncl
            Write (IW6, 115) ib, ncl
115         Format (/ '    SITE:    IB=', I3, '        CLUSTER SIZE= ', I3)

!------------------------ TB-CONSTANTS FOR THE SMALL CLUSTER AND  THEIR STORAGE
            Call TBCL (ncl, POL, WSL, ssc, opt%nl, opt%nlsq, opt%mtbcl, system%qscr, harmon)

            Do jcl = 1, ncl
               Do jq = 1, opt%nlsq
                  Do iq = 1, opt%nlsq
                     bstr (iq, jq, jcl, ib) = ssc (iq, jq, jcl)
                  End Do
               End Do
            End Do

         End Do

!----------------------------------------- SYMMETRIZATION
!!$                       -------- IN-LAYER CONSTANTS
         DEV = rcz
         Do ib = 1, system%NB
            Do jcl = 1, NSCLB (ib)
               If (JDRPB(jcl, ib) .Ne. 0) Cycle
               jb = JBVAB (jcl, ib)
               Do icl = 1, NSCLB (jb)
                  If (JDRPB(icl, jb) .Ne. 0) Cycle
                  If (JBVAB(icl, jb) .Ne. ib) Cycle
                  isum = JTRAB (1, jcl, ib) + JTRAB (1, icl, jb)
                  If (isum .Ne. 0) Cycle
                  isum = JTRAB (2, jcl, ib) + JTRAB (2, icl, jb)
                  If (isum .Ne. 0) Cycle

                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMA (iq, jq) = bstr (iq, jq, jcl, ib)
                     End Do
                  End Do

                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMB (iq, jq) = bstr (iq, jq, icl, jb)
                     End Do
                  End Do

                  Call MASY (WWMA, WWMB, DEV1, opt%nlsq)

                  DEV = Max (DEV, DEV1)
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        bstr (iq, jq, jcl, ib) = WWMA (iq, jq)
                     End Do
                  End Do

                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        bstr (iq, jq, icl, jb) = WWMB (iq, jq)
                     End Do
                  End Do
               End Do
            End Do
         End Do

         Write (IW6, 125) DEV
125      Format (/ 4 X, ' IN-LAYER CONSTANTS : DEVIATION=', G12.4)

!!$                       -------- OFF-LAYER CONSTANTS
         DEV = rcz

         Do ib = 1, system%NB
            Do jcl = 1, NSCLB (ib)
               If (JDRPB(jcl, ib) .Ne. 1) Cycle
               jb = JBVAB (jcl, ib)
               Do icl = 1, NSCLB (jb)
                  If (JDRPB(icl, jb) .Ne.-1) Cycle
                  If (JBVAB(icl, jb) .Ne. ib) Cycle
                  isum = JTRAB (1, jcl, ib) + JTRAB (1, icl, jb)
                  If (isum .Ne. 0) Cycle
                  isum = JTRAB (2, jcl, ib) + JTRAB (2, icl, jb)
                  If (isum .Ne. 0) Cycle
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMA (iq, jq) = bstr (iq, jq, jcl, ib)
                     End Do
                  End Do
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMB (iq, jq) = bstr (iq, jq, icl, jb)
                     End Do
                  End Do

                  Call MASY (WWMA, WWMB, DEV1, opt%nlsq)

                  DEV = Max (DEV, DEV1)
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        bstr (iq, jq, jcl, ib) = WWMA (iq, jq)
                     End Do
                  End Do

                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        bstr (iq, jq, icl, jb) = WWMB (iq, jq)
                     End Do
                  End Do
               End Do
            End Do
         End Do


135      Format (/ 4 X, 'OFF-LAYER CONSTANTS : DEVIATION=', G12.4)
150      Format (/ 5 X, '***  DIAGONAL ELEMENTS OF ', 'TB-STRUCTURE CONSTANTS :')
151      Format (/ '        SITE:       IB=', I3)
104      Format (1 X, 4 G15.7)
191      Format (/ ' **** ERROR IN TBRB :' / 10 X, ' JCL GREATER THAN mncl,  JCL=', I3)

         Write (IW6, 135) DEV
!---------------------- PRINT OF DIAGONAL  ELEMENTS
         Write (IW6, 150)
         Do ib = 1, system%NB
            Write (IW6, 151) ib
            Write (IW6, 104) (bstr(iq, iq, 1, ib), iq=1, opt%nlsq)
         End Do

         Return

291      Write (IW6, 191) jcl
         Stop
      End Subroutine tbrb

!*******************
!XXX    TBRV    ****
!*******************
      Subroutine tbrv (opt, system, harmon, vstr, JBVAV, JDRPV, JTRAV, NSCLV)
!----------------------------------------------------
!!$   TB-STRUCTURE CONSTANTS FOR SEMiiNFINITE VACUUM
!----------------------------------------------------
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Type (sharm) :: harmon
         Real (Kind=prec) :: POL (3, MNCL), WSL (MNCL), POLT (3), DPOL (3), WWMA (opt%nlsq, opt%nlsq), WWMB &
        & (opt%nlsq, opt%nlsq)
         Real (Kind=prec) :: DMAXSQ, DMINSQ, sum, DEV, DEV1, ssc (opt%nlsq, opt%nlsq, MNCL)
         Real (Kind=prec), Parameter :: rcz = 0.0d0, SMALL = 0.01d0
         Integer :: ib, jcl, icl, jdir, jb, jtra1, jtra2, i, iq, jq, ncl, isum
         Real (Kind=prec), Pointer :: vstr (:, :, :, :)
         Integer, Pointer :: JDRPV (:, :), JTRAV (:, :, :), JBVAV (:, :), NSCLV (:)

         ssc = 0.0d0

         Write (IW6, 111)
111      Format (/ / ' ***********  TB-CONSTANTS FOR', ' THE VACUUM REGION  ***********')

!------------------------------- CUT-OFF DISTANCE
         DMAXSQ = (system%cutrat*system%avws) ** 2
         DMINSQ = (SMALL*system%avws) ** 2

         Do ib = 1, system%NB !--------------------- LOOP OVER CENTRAL ATOMS FOR SMALL CLUSTERS
            jcl = 1 !--------------------- SELECTION OF THE SMALL CLUSTER
            Do i = 1, 3
               POL (i, 1) = system%vpos(i, ib)
            End Do
            JDRPV (1, ib) = 0
            JBVAV (1, ib) = ib
            JTRAV (1, 1, ib) = 0
            JTRAV (2, 1, ib) = 0
            WSL (1) = system%avws
            Do jdir = - 1, 1
               Do jb = 1, system%NB
                  Do jtra1 = - system%nmtr, system%nmtr
                     Do jtra2 = - system%nmtr, system%nmtr
                        Do i = 1, 2
                           POLT (i) = system%vpos(i, jb) - jdir * system%VTRV(i) + jtra1 * system%vbr(i, 1) + &
                          & jtra2 * system%vbr(i, 2)
                        End Do
                        POLT (3) = system%vpos(3, jb) - jdir * system%VTRV(3)
                        Do i = 1, 3
                           DPOL (i) = POLT (i) - POL (i, 1)
                        End Do
                        sum = DPOL (1) ** 2 + DPOL (2) ** 2 + DPOL (3) ** 2
                        If (sum .Lt. DMINSQ .Or. sum .Gt. DMAXSQ) Cycle
                        jcl = jcl + 1
                        If (jcl .Gt. MNCL) Go To 291
                        Do i = 1, 3
                           POL (i, jcl) = POLT (i)
                        End Do
                        JDRPV (jcl, ib) = jdir
                        JBVAV (jcl, ib) = jb
                        JTRAV (1, jcl, ib) = jtra1
                        JTRAV (2, jcl, ib) = jtra2
                        WSL (jcl) = system%avws
                     End Do
                  End Do
               End Do
            End Do

            ncl = jcl
            NSCLV (ib) = ncl
            Write (IW6, 115) ib, ncl
115         Format (/ '    SITE:    IB=', I3, '        CLUSTER SIZE= ', I3)

!------ TB-CONSTANTS FOR THE SMALL CLUSTER AND THEIR STORAGE

            Call TBCL (ncl, POL, WSL, ssc, opt%nl, opt%nlsq, opt%mtbcl, system%qscr, harmon)
!SMP$ DO SERIAL
            Do jcl = 1, ncl
               Do jq = 1, opt%nlsq
                  Do iq = 1, opt%nlsq
                     vstr (iq, jq, jcl, ib) = ssc (iq, jq, jcl)
                  End Do
               End Do
            End Do
         End Do

!----------------------------------------- SYMMETRIZATION
!!$                       -------- IN-LAYER CONSTANTS
         DEV = rcz
         Do ib = 1, system%NB
            Do jcl = 1, NSCLV (ib)
               If (JDRPV(jcl, ib) .Ne. 0) Cycle
               jb = JBVAV (jcl, ib)
               Do icl = 1, NSCLV (jb)
                  If (JDRPV(icl, jb) .Ne. 0) Cycle
                  If (JBVAV(icl, jb) .Ne. ib) Cycle
                  isum = JTRAV (1, jcl, ib) + JTRAV (1, icl, jb)
                  If (isum .Ne. 0) Cycle
                  isum = JTRAV (2, jcl, ib) + JTRAV (2, icl, jb)
                  If (isum .Ne. 0) Cycle
!SMP$ DO SERIAL
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMA (iq, jq) = vstr (iq, jq, jcl, ib)
                     End Do
                  End Do
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMB (iq, jq) = vstr (iq, jq, icl, jb)
                     End Do
                  End Do
                  Call MASY (WWMA, WWMB, DEV1, opt%nlsq)
                  DEV = Max (DEV, DEV1)
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        vstr (iq, jq, jcl, ib) = WWMA (iq, jq)
                     End Do
                  End Do
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        vstr (iq, jq, icl, jb) = WWMB (iq, jq)
                     End Do
                  End Do
               End Do
            End Do
         End Do

         Write (IW6, 125) DEV
125      Format (/ 4 X, ' IN-LAYER CONSTANTS : DEVIATION=', G12.4)
!!$                       -------- OFF-LAYER CONSTANTS
         DEV = rcz

         Do ib = 1, system%NB
            Do jcl = 1, NSCLV (ib)
               If (JDRPV(jcl, ib) .Ne. 1) Cycle
               jb = JBVAV (jcl, ib)
               Do icl = 1, NSCLV (jb)
                  If (JDRPV(icl, jb) .Ne.-1) Cycle
                  If (JBVAV(icl, jb) .Ne. ib) Cycle
                  isum = JTRAV (1, jcl, ib) + JTRAV (1, icl, jb)
                  If (isum .Ne. 0) Cycle
                  isum = JTRAV (2, jcl, ib) + JTRAV (2, icl, jb)
                  If (isum .Ne. 0) Cycle
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMA (iq, jq) = vstr (iq, jq, jcl, ib)
                     End Do
                  End Do
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        WWMB (iq, jq) = vstr (iq, jq, icl, jb)
                     End Do
                  End Do

                  Call MASY (WWMA, WWMB, DEV1, opt%nlsq)

                  DEV = Max (DEV, DEV1)
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        vstr (iq, jq, jcl, ib) = WWMA (iq, jq)
                     End Do
                  End Do
                  Do jq = 1, opt%nlsq
                     Do iq = 1, opt%nlsq
                        vstr (iq, jq, icl, jb) = WWMB (iq, jq)
                     End Do
                  End Do
               End Do
            End Do
         End Do


135      Format (/ 4 X, 'OFF-LAYER CONSTANTS : DEVIATION=', G12.4)
150      Format (/ 5 X, '***  DIAGONAL ELEMENTS OF ', 'TB-STRUCTURE CONSTANTS :')
151      Format (/ '        SITE:       IB=', I3)
104      Format (1 X, 4 G15.7)

         Write (IW6, 135) DEV
!------------------------ PRINT OF DIAGONAL ELEMENTS
         Write (IW6, 150)
         Do ib = 1, system%NB
            Write (IW6, 151) ib
            Write (IW6, 104) (vstr(iq, iq, 1, ib), iq=1, opt%nlsq)
         End Do

         Return

291      Write (IW6, 191) jcl
191      Format (/ ' **** ERROR IN TBRV :' / 10 X, ' JCL GREATER THAN mncl,  JCL=', I3)
         Stop
      End Subroutine tbrv
#if defined(DO_SC_SYM_TUREK)
      Subroutine getscsym (N, SC, f, crit)
         Implicit None
         Integer, Intent (In) :: N
         Real (Kind=prec), Intent (In) :: crit
         Real (Kind=prec), Intent (Inout) :: SC (N)
         Integer, Intent (Out) :: f (2*N)
!!$  Local
         Integer, Pointer :: mask (:)
         Integer :: i, j, mptr, num
         Real (Kind=prec) :: base
         Allocate (mask(N))
         mask (1:N) = 0
         f (1) = 0
         mptr = 1
         Do i = 1, N, 1
            base = SC (i)
!!$          If (Abs(base) > crit .And. mask(i) == 0) Then
            If (mask(i) == 0) Then
               num = 1
               Do j = i + 1, N, 1
!!$                If ((mask(j) == 0) .And. (Abs(sc(j)) > crit)) Then
                  If ((mask(j) == 0)) Then
                     If (Abs(base-SC(j)) < crit) Then
                        mask (j) = 1
                        num = num + 1
                        f (mptr+num) = j
                     Else
                        If (Abs(base+SC(j)) < crit) Then
                           mask (j) = 1
                           num = num + 1
                           f (mptr+num) = - j
                        End If
                     End If
                  End If
               End Do
               If (num > 1) Then
                  f (mptr) = num
                  f (mptr+1) = i
                  mptr = mptr + num + 1
                  f (mptr) = 0
               End If
            End If
         End Do
         Deallocate (mask)
      End Subroutine getscsym

      Subroutine setscsym (N, SC, f)
         Implicit None
         Integer, Intent (In) :: N
         Real (Kind=prec), Intent (Inout) :: SC (N)
         Integer, Intent (In) :: f (2*N)
!!$  Local
         Integer :: i, mptr, cntr
         Real (Kind=prec) :: base
         mptr = 1
         Do While (f(mptr) /=  0)
            base = 0.0d0
            cntr = f (mptr)
!!$           write(*,*)
            Do i = mptr + 1, mptr + cntr, 1
               base = base + real (sign(1, f(i)), kind=prec) * SC (Abs(f(i)))
!!$               write(*,*) sc(abs(f(i)))
            End Do
            base = base / real (cntr, kind=prec)
            Do i = mptr + 1, mptr + cntr, 1
               SC (Abs(f(i))) = real (sign(1, f(i)), kind=prec) * base
            End Do
            mptr = mptr + cntr + 1
         End Do
      End Subroutine setscsym

      Subroutine symscons (N, SC, crit)
         Implicit None
         Integer, Intent (In) :: N
         Real (Kind=prec), Intent (Inout) :: SC (N)
         Real (Kind=prec) :: crit
!!$ Locals
         Integer :: i, j, num
         Integer, Pointer :: mask (:), refs (:)
         Real (Kind=prec) :: base, DEV, abase, aref, acc, ref
         Allocate (mask(N))
         Allocate (refs(N))
         mask (1:N) = 0
!!$       write(*,*)crit
         Do i = 1, N, 1
            If (mask(i) == 0) Then
               base = SC (i)
               abase = Abs (base)
               acc = abase
               num = 1
               refs (1) = i * sign (1.0d0, base)
               num = 1
               Do j = i + 1, N, 1
                  If ((mask(j) == 0)) Then
                     ref = SC (j)
                     aref = Abs (ref)
                     DEV = Abs (2.0d0*(abase-aref)/(abase+aref))
!!$                   write(*,*)'dev=',dev
                     If (DEV <= crit) Then
                        num = num + 1
                        refs (num) = j * sign (1.0d0, ref)
!!$                      write(*,*) base, ref
                        acc = acc + aref
                        mask (j) = 1
                     End If
                  End If
               End Do

               If (num > 1) Then
                  acc = acc / real (num, kind=prec)
                  Do j = 1, num, 1
                     SC (Abs(refs(j))) = real (sign(1, refs(j)), kind=prec) * acc
                  End Do
               End If
            End If
         End Do
         Deallocate (mask, refs)
      End Subroutine symscons
#endif

End Module ies_structure
