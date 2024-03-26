!!$   $Id: definitions.f,v 1.39 2003/11/03 19:04:43 maciej Exp $
Module read_write
!!$    This module contains the definitions of data structures used elsewhere
!!$    in the program and declarations od some constants.
!!$    prec   - kind value for double precision
!!$    lprec  - kind value for double double (extended double) precision
!!$    clmax  - max. size of the cluster for S^alpha calculations.
      Use definitions
      Use logging
      Implicit None

Contains

!*******************
!XXX    RALL    ****
!*******************
      Subroutine rall (opt, system, cnodes)
!***************************************
!!$   INPUT OF ALL DATA
!***************************************
         Implicit None
         Character (Len=16) :: OWORK (5), OCASE (0:1), ORELA (2, 0:2), OXCP (2, 0:2)
         Integer :: i, j, ip, ib, ia, ic, il, is, ie, ne1, ix, nc, nsv, nrad1
         Real (Kind=prec) :: dummy1, dummy2, skipPOT
         Character (Len=5) :: owork1
         Type (atom_definition), Pointer :: aa
         Type (bulk_atom_definition), Pointer :: ab
         Type (params) :: opt
         Type (atoms_set) :: system
         Type (complex_nodes) :: cnodes
         Character (Len=16) :: cwork, cwork1
         Data OCASE / '    SURFACE     ', '   INTERFACE    ' /
         Data ORELA / '      NONE      ', '                ', '      SCALAR    ', '                ', '      &
        &SCALAR    ', '( CORES - FULL )' /
         Data OXCP / '    VON BARTH - ', 'HEDIN           ', '  CEPERLEY-ALDER', ' (PERDEW-ZUNGER)', '    VOS&
        &KO - WILK', ' - NUSAIR       ' /

         system%nrmax = 0
100      Format (5 A16)
101      Format (1 X, 10 I5)
104      Format (1 X, 4 G15.9)

         Open (Unit=IR1, File='inpge', Action='read')
         Open (Unit=IR2, File='inpch', Action='read')
!!$      open(unit=IR4,file='inpcp',action='read')
         Open (Unit=IR5, File='input', Action='read')
         Open (Unit=IR11, File='inpbu', Action='read')
         Open (Unit=IR13, File='cci.nw', Action='read')
         IW6 = log_stream_number
         Open (Unit=IW6, File='outit', Action='write')
         have_log_open = 1
!!$      open(unit=IW8,file='outcp',action='write')

         Write (IW6, 111)
111      Format (/ '       ****  OUTPUT OF IES  **** ' // '     ** READING OF ALL DATA:')

!-------------------------  CONTROL DATA - UNIT IR5

1001     Format (/ / 5 A16)
102      Format (1 X, F10.5, I5)
107      Format (1 X, 2 F10.5)

         Read (IR5, 100) OWORK
         Write (IW6, 1001) OWORK
         Read (IR5, 100) OWORK
         Write (IW6, 1001) OWORK

         Read (IR5, 101) opt%ivac, opt%nl, opt%NS
         If (opt%ivac < 0 .Or. opt%ivac > 1) Then
            Write (IW6, '(" **** INPUT ERROR:  ivac MUST BE 0 OR 1 ")')
            Stop
         End If
         If (opt%nl < 3 .Or. opt%nl > 4) Then
            Write (IW6, '(" **** INPUT ERROR:  nl MUST BE 3 OR 4 ")')
            Stop
         End If

         Read (IR5, 101) opt%IREL, opt%IVXC

         If (opt%IREL < 0 .Or. opt%IREL > 2) Then
            Write (IW6, '(" **** INPUT ERROR:  IREL MUST BE 0, 1 OR 2 ")')
            Stop
         End If
         If (opt%NS == 2 .And. opt%IREL == 2) Then
            Write (IW6, '(" **** INPUT ERROR:  NS=2 AND IREL=2 ")')
            Stop
         End If
         If (opt%IVXC < 0 .Or. opt%IVXC > 2) Then
            Write (IW6, '(" **** INPUT ERROR:  IVXC MUST BE 0, 1 OR 2 ")')
            Stop
         End If

         Read (IR5, 101) opt%bz%NSYM, opt%bz%NK, opt%bz%INVE
         Read (IR5, 101) opt%NE, MNCL
         If (opt%NE /= 7 .And. opt%NE /= 10 .And. opt%NE /= 14 .And. opt%NE /= 20 .And. opt%NE /= 28 .And. &
        & opt%NE /= 40 .And. opt%NE /= 56 .And. opt%NE /= 80) Then
            Write (IW6, '(" *** INPUT ERROR: NE MUST BE 7,10,14,20,28,40,56, OR 80 ")')
            Stop
         End If

         Read (IR5, 107) opt%DIAM

         Read (IR5, 101) opt%icont, opt%rwgamma
         If (opt%icont == 0) Then
            Write (IW6, 205) opt%icont
205         Format (5 X, 'icont=', I5, ' new calculations ')
         Else If (opt%icont == 1) Then
            Write (IW6, 206) opt%icont
206         Format (5 X, 'icont=', I5, ' continue with calculations ')
         Else
            Write (IW6, 207) opt%icont
207         Format (5 X, 'icont=', I5, ' but must be 0 or 1')
            Stop
         End If

         Read (IR5, 101) opt%NFUPR
         Read (IR5, 107) opt%DifMC, opt%QMIXC
         Read (IR5, 107) opt%DifMS, opt%QMIXS
         Read (IR5, 102) opt%ALFA, opt%NITER
         Read (IR5, 102) opt%BETA, opt%NITERA
         opt%NITERA = Max (opt%NITERA, 2)

         Read (IR5, 102) opt%W0AM, opt%nuh

         opt%nuh = Max (opt%nuh, 2)

         Read (IR5, 101) opt%lident
         Read (IR5,*) opt%tol

         If (opt%lident < 0 .Or. opt%lident > 1) Then
            Write (IW6, '(" **** INPUT ERROR:  LIDENT MUST BE 0 OR 1 ")')
            Stop
         End If

         Write (IW6, 115) OCASE (opt%ivac), opt%nl, opt%NS
115      Format (/ 4 X, 'CASE:', A16, 7 X, 'NL=', I2, 7 X, 'NS=', I2)
         Write (IW6, 1161) (ORELA(i, opt%IREL), i=1, 2)
1161     Format (/ 8 X, 'RELATIVITY:   ', 2 A16)
         Write (IW6, 1162) (OXCP(i, opt%IVXC), i=1, 2)
1162     Format (/ 8 X, 'XC-POTENTIAL:   ', 2 A16)
         Write (IW6, 117) opt%bz%NSYM, opt%bz%NK, opt%bz%INVE
117      Format (/ '       LOCAL SYMMETRY:    NSYM=', I2 // '   K||-MESH:   NK=', I3, 10 X, 'INVE=', I2)
         Write (IW6, 118) opt%NE, opt%DIAM
118      Format (/ '  ENERGY CONTOUR:    NE=', I3, 12 X, 'DIAM=', F12.5)
         Write (IW6, 119) opt%lident
119      Format (/ '  LEADS IDENTICAL:  LIDENT=', I1)
         Write (IW6, 120) opt%NFUPR
120      Format (/ '          FULL PRINT AFTER', I4, ' ITERATIONS')
         Write (IW6, 121) opt%DifMC, opt%QMIXC
121      Format (/ '  CPA ITERATIONS:  DIFMAX=', G12.4, 8 X, 'MIXING=', F10.5)
         Write (IW6, 122) opt%DifMS, opt%QMIXS
122      Format (/ '  SGF ITERATIONS:  DIFMAX=', G12.4, 8 X, 'MIXING=', F10.5)
         Write (IW6, 125) opt%ALFA, opt%NITER
125      Format (/ '  LDA ITERATIONS: ' / 5 X, 'STARTING MIXING=', F10.5, '    TOTAL:', I5, ' ITERATIONS')
         Write (IW6, 126) opt%BETA, opt%NITERA
126      Format (5 X, 'ANDERSON MIXING=', F10.5, '     AFTER', I5, ' ITERATIONS')
         Write (IW6, 127) opt%W0AM, opt%nuh
127      Format (5 X, 'QUANTITY W0=', F10.5, 4 X, '  PREVIOUS', I5, ' ITERATIONS')
         Write (IW6, 1281) opt%tol
1281     Format (/ '  STOP CRITERIA:  TOL=', G12.4)
         Write (IW6, 1282) MNCL
1282     Format (/ '  MAXIMAL CLUSTER SIZE :  MNCL=', I3)




!-------------------------  GEOMETRY - UNIT IR1

         Read (IR1, 100) OWORK
         Write (IW6, 1001) OWORK

         Read (IR1, 101) system%NP, system%NB, system%NMTR
         system%nums = system%NP * system%NB
         Allocate (system%nc(system%nums))
         Allocate (system%pos(3, system%NB, system%NP))
         Allocate (system%vpos(3, system%NB))
         Allocate (system%bpos(3, system%NB))
!!$ 10414   format(1X,4G15.9)
         Read (IR1, 104) system%CUTRAT
         Read (IR1, 104) (system%SCX(i), i=1, 3)
         Read (IR1, 104) system%VBR(1, 1), system%VBR(2, 1)
         Read (IR1, 104) system%VBR(1, 2), system%VBR(2, 2)
         Read (IR1, 104) (system%VTRV(i), i=1, 3)
         Read (IR1, 104) (system%BTRV(i), i=1, 3)
         Do ib = 1, system%NB
            Read (IR1, 104) (system%vpos(i, ib), i=1, 3)
         End Do

         Do ip = 1, system%NP
            Do ib = 1, system%NB
               Read (IR1, 104) (system%pos(i, ib, ip), i=1, 3)
            End Do
         End Do

         Do ib = 1, system%NB
            Read (IR1, 104) (system%bpos(i, ib), i=1, 3)
         End Do
!!$                          SCALING OF COORDINATES
         Do j = 1, 2
            Do i = 1, 2
               system%VBR (i, j) = system%SCX(i) * system%VBR(i, j)
            End Do
         End Do

         Do i = 1, 3
            system%VTRV (i) = system%SCX(i) * system%VTRV(i)
            system%BTRV (i) = system%SCX(i) * system%BTRV(i)
         End Do

         Do ib = 1, system%NB
            Do i = 1, 3
               system%vpos (i, ib) = system%SCX(i) * system%vpos(i, ib)
               system%bpos (i, ib) = system%SCX(i) * system%bpos(i, ib)
            End Do
         End Do

         Do ip = 1, system%NP
            Do ib = 1, system%NB
               Do i = 1, 3
                  system%pos (i, ib, ip) = system%SCX(i) * system%pos(i, ib, ip)
               End Do
            End Do
         End Do

         Write (IW6, 130) system%NP, system%NB, system%NMTR, system%CUTRAT
130      Format (/ 8 X, 'NP=', I4, 8 X, 'NB=', I2 // '   MAX. ', 'COEFFICIENT OF 2D-TRANSL. VECTORS:    NMTR=&
        &', I2 / / '   CUT-OFF DISTANCE / WS-RADIUS:  CUTRAT=', F18.12)
         Write (IW6, 131) (system%SCX(i), i=1, 3)
131      Format (/ '  SCALING FACTORS: ', 3 G22.14)
         Write (IW6, 132) ((system%VBR(i, j), i=1, 2), j=1, 2)
132      Format (/ '  TRANSL. VECTORS: 1. VECTOR = ', 2 G22.14 / '                   2. VECTOR = ', 2 G22.14)
         Write (IW6, 133) (system%VTRV(i), i=1, 3)
133      Format (/ '    TRANSL. VECTOR OF SEMIINFINITE VACUUM:' / 4 X, 3 G22.14)
         Write (IW6, 134) (system%BTRV(i), i=1, 3)
134      Format (/ '    TRANSL. VECTOR OF SEMIINFINITE BULK:' / 4 X, 3 G22.14)
1041     Format (3 X, 'IB=', I2, 5 X, 3 G22.14)
         Write (IW6, 135)
135      Format (/ '       SITES OF VACUUM PRINCIPAL LAYER:')
         Do ib = 1, system%NB
            Write (IW6, 1041) ib, (system%vpos(i, ib), i=1, 3)
         End Do

         Do ip = 1, system%NP
            Write (IW6, 136) ip
            Do ib = 1, system%NB
               Write (IW6, 1041) ib, (system%pos(i, ib, ip), i=1, 3)
            End Do
         End Do
136      Format (/ '       SITES OF PRINCIPAL LAYER NO.: ', I4)
         Write (IW6, 137)
137      Format (/ '       SITES OF BULK PRINCIPAL LAYER: ')

         Do ib = 1, system%NB
            Write (IW6, 1041) ib, (system%bpos(i, ib), i=1, 3)
         End Do

!------------------------- CHEMICAL OCCUPATION - UNIT IR2
108      Format (A16)
1081     Format (A16, A16)
         Read (IR2, 100) OWORK
         is = 0 ! sites
         ia = 0 ! atoms
         Do ip = 1, system%NP
            Do ib = 1, system%NB
               is = is + 1 ! site counter
               Read (IR2, 101) nc
               Do ic = 1, nc
                  ia = ia + 1 ! atom counter
                  Read (IR2, 108) owork1
                  Read (IR2, 104) dummy1, dummy2
               End Do
            End Do
         End Do
         system%num = ia ! total number of atoms

         Allocate (system%at(system%num))

         Rewind IR2
         Read (IR2, 100) OWORK
         Write (IW6, 1001) OWORK
         is = 0 ! sites
         ia = 0 ! atoms
         Do ip = 1, system%NP
            Do ib = 1, system%NB
               is = is + 1 ! site counter
               Read (IR2, 101) system%nc(is)
               Do ic = 1, system%nc(is)
                  ia = ia + 1 ! atom counter
                  aa => system%at (ia)
                  Read (IR2, 1081) aa%otxta, aa%otxta_in
                  Read (IR2, 104) aa%con, aa%valz
               End Do
            End Do
         End Do

141      Format (/ '  SITE:  IP=', I4, '   IB=', I2, '   IG=', I4, '   OCCUPIED BY ', I2, ' ATOMIC TYPE(S): ')
1421     Format (' IA=', I4, ' - LABEL,LABEL_IN,CONC.,VALENCY: ', A16, A16, F10.5, F10.3)
143      Format (/ 4 X, ' TOTAL NUMBER OF DIFFERENT ', 'ATOMIC TYPES:  NA= ', I4)

         is = 0
         ia = 0
         Do ip = 1, system%NP
            Do ib = 1, system%NB
               is = is + 1
               Write (IW6, 141) ip, ib, is, system%nc(is)
               Do ic = 1, system%nc(is)
                  ia = ia + 1
                  aa => system%at (ia)
                  Write (IW6, 1421) ia, aa%otxta, aa%otxta_in, aa%con, aa%valz
               End Do
            End Do
         End Do
         Write (IW6, 143) system%num

!!$ ---------------------------- LSDA FILE - UNIT IR3

         Do ia = 1, system%num
            aa => system%at (ia)
            Allocate (aa%eny(opt%nl, opt%NS))
            Allocate (aa%ppc(opt%nl, opt%NS))
            Allocate (aa%ppd(opt%nl, opt%NS))
            Allocate (aa%ppq(opt%nl, opt%NS))
            Allocate (aa%ppp(opt%nl, opt%NS))
            Allocate (aa%dny(opt%nl, opt%NS))
            Allocate (aa%finy(opt%nl, opt%NS))
            Allocate (aa%finyd(opt%nl, opt%NS))
            Allocate (aa%ynam(opt%nl, opt%NS))
         End Do

         opt%pair = (opt%nl-1) ** 2

         Do ia = 1, system%num
            aa => system%at (ia)
            If (opt%icont == 0) Then
               Read (aa%otxta_in,*) cwork
               Read (aa%otxta,*) cwork1
               Write (IW6,*) ' Read defaults for:', trim (cwork1), ' from:', trim (cwork)
            Else
               Read (aa%otxta,*) cwork
               Write (IW6,*) ' Read previous resuts for:', trim (cwork)
            End If
            Open (Unit=atomf, File='atoms/'//trim(cwork), Action='read')
            Do is = 1, 3
               Read (atomf,*) cwork
            End Do

            Read (atomf, 104) aa%az, aa%ws, aa%wsav
            Do is = 1, 4
               Read (atomf,*) cwork
            End Do
            Do is = 1, opt%NS
               If (is == 1) Then
                  Read (atomf, 101) aa%nszrad
                  system%nrmax = Max (system%nrmax, aa%nszrad)
                  Allocate (aa%pot(aa%nszrad, opt%NS))!  potential
                  Allocate (aa%fpot(aa%nszrad, opt%NS))!  potential
                  Allocate (aa%phi(aa%nszrad, opt%nl, opt%NS))!  wave function
                  Allocate (aa%phid(aa%nszrad, opt%nl, opt%NS))!  first derivative of wave function
                  Allocate (aa%phidd(aa%nszrad, opt%nl, opt%NS))!  second derivative
                  Allocate (aa%rhocor(aa%nszrad, opt%NS))!  core density
                  Allocate (aa%rhoval(aa%nszrad, opt%NS))!  valence density
                  Allocate (aa%emdi0(opt%nl, opt%NS))!  zero diagonal energy moments
                  Allocate (aa%emdi1(opt%nl, opt%NS))!  first diagonal energy moments
                  Allocate (aa%emdi2(opt%nl, opt%NS))!  second diagonal energy moments
                  Allocate (aa%emof00(opt%pair, opt%NS))!  off-diagonal energy moments
                  Allocate (aa%emof10(opt%pair, opt%NS))
                  Allocate (aa%emof01(opt%pair, opt%NS))
                  Allocate (aa%emof11(opt%pair, opt%NS))
                  Allocate (aa%emof20(opt%pair, opt%NS))
                  Allocate (aa%emof02(opt%pair, opt%NS))
               Else
                  Read (atomf,*) cwork
               End If
               Read (atomf, 104) (aa%pot(i, is), i=1, aa%nszrad)
            End Do
            Read (atomf,*) cwork
            Do is = 1, opt%NS
               Read (atomf, 104) cwork
               Read (atomf, 104) (aa%eny(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%ppc(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%ppd(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%ppq(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%ppp(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%dny(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%finy(il, is), il=1, opt%nl)
               Read (atomf, 104) (aa%finyd(il, is), il=1, opt%nl)
            End Do
            Read (atomf, 101) aa%numcor
            If (aa%numcor /= 0) Then
               Allocate (aa%ncor(aa%numcor))
               Allocate (aa%lcor(aa%numcor))
               Allocate (aa%nobc(aa%numcor))
               Allocate (aa%ecor(aa%numcor, 2))
               Do j = 1, aa%numcor
                  Read (atomf, 101) aa%ncor(j), aa%lcor(j), aa%nobc(j)
                  Read (atomf, 104) (aa%ecor(j, ix), ix=1, 2)
               End Do
            End If
            If (opt%ivac == 0) read (atomf, 104) system%dba
            Close (atomf)
         End Do

         If (opt%ivac == 0) write (IW6, 149) system%dba
149      Format (/ 6 X, ' DIPOLE BARRIER=', F12.5)
         Write (IW6, 150)
150      Format (/ 6 X, '  LABEL,         IA,', '      Z,         WS,         WSAV: ')

         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 151) aa%otxta, ia, aa%az, aa%ws, aa%wsav
         End Do

151      Format (4 X, A16, I5, 3 F12.5)
         Write (IW6, 153)
153      Format (/ 6 X, '  LABEL,         IA,', '   NSZRAD, NUMCOR:')

         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 154) aa%otxta, ia, aa%nszrad, aa%numcor
         End Do

154      Format (4 X, A16, I5, 2 I8)

!!$ ---------------  DATA FOR 2ND (VAC.) SUBSTRATE - UNIT IR12
         If (opt%ivac == 1) Then
            Read (IR11, 100) OWORK
            Write (IW6, 1001) OWORK
            Read (IR11, 101) system%vnum, system%vns
            Allocate (system%vat(system%vnum))
            Do ia = 1, system%vnum
               ab => system%vat (ia)
               Allocate (ab%eny(opt%nl, system%vns))
               Allocate (ab%ppc(opt%nl, system%vns))
               Allocate (ab%ppd(opt%nl, system%vns))
               Allocate (ab%ppq(opt%nl, system%vns))
               Allocate (ab%ppp(opt%nl, system%vns))
               Allocate (ab%dny(opt%nl, system%vns))
               Allocate (ab%finy(opt%nl, system%vns))
               Allocate (ab%finyd(opt%nl, system%vns))
               Read (IR11, 104) ab%con
               Read (IR11, 108) ab%otxta
               Read (ab%otxta,*) cwork
               Write (IW6,*) ' Read left lead(vac.) from:', cwork
               Open (Unit=atomf, File='atoms/'//trim(cwork), Action='read')
               Read (atomf,*) cwork
               Read (atomf,*) cwork
               Read (atomf,*) cwork
               Read (atomf, 104) ab%az, ab%ws, ab%wsav
               Read (atomf, 104) system%vef
               Read (atomf,*) cwork
               Read (atomf,*) cwork
               Read (atomf,*) nsv
               If (nsv < system%vns) Then
                  Write (*,*) 'Error vac01'
                  Stop
               End If
               Do is = 1, nsv
                  Read (atomf, 101) nrad1
                  Read (atomf, 104) (skipPOT, i=1, nrad1)
               End Do
               Read (atomf,*) cwork
               Do is = 1, system%vns
                  Read (atomf,*) cwork
                  Read (atomf, 104) (ab%eny(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%ppc(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%ppd(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%ppq(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%ppp(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%dny(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%finy(il, is), il=1, opt%nl)
                  Read (atomf, 104) (ab%finyd(il, is), il=1, opt%nl)
               End Do
               Close (atomf)
            End Do
            Write (IW6, 170) system%vef, system%vnum, system%vns
170         Format (/ 4 X, 'VAC. SUBSTRATE:  EF=', G15.7, 8 X, 'NAV=', I2, 8 X, 'NSV=', I2)
            Write (IW6, 161)
            Do ia = 1, system%vnum
               ab => system%vat (ia)
               Write (IW6, 162) ia, ab%con, ab%az, ab%ws, ab%wsav
            End Do
         End If
!!$ ---------------  DATA FOR 1ST (BULK) SUBSTRATE - UNIT IR11

         Read (IR11, 100) OWORK
         Write (IW6, 1001) OWORK
         Read (IR11, 101) system%bnum, system%bns
         Allocate (system%bat(system%bnum))
         Do ia = 1, system%bnum
            ab => system%bat (ia)
            Allocate (ab%eny(opt%nl, system%bns))
            Allocate (ab%ppc(opt%nl, system%bns))
            Allocate (ab%ppd(opt%nl, system%bns))
            Allocate (ab%ppq(opt%nl, system%bns))
            Allocate (ab%ppp(opt%nl, system%bns))
            Allocate (ab%dny(opt%nl, system%bns))
            Allocate (ab%finy(opt%nl, system%bns))
            Allocate (ab%finyd(opt%nl, system%bns))
            Read (IR11, 104) ab%con
            Read (IR11, 108) ab%otxta
            Read (ab%otxta,*) cwork
            Write (IW6,*) ' Read right lead(bulk.) from:', cwork
            Open (Unit=atomf, File='atoms/'//trim(cwork), Action='read')
            Read (atomf,*) cwork
            Read (atomf,*) cwork
            Read (atomf,*) cwork
            Read (atomf, 104) ab%az, ab%ws, ab%wsav
            Read (atomf, 104) system%bef
            Read (atomf,*) cwork
            Read (atomf,*) cwork
            Read (atomf,*) nsv
            If (nsv < system%bns) Then
               Write (*,*) 'Error vac01'
               Stop
            End If
            Do is = 1, nsv
               Read (atomf, 101) nrad1
               Read (atomf, 104) (skipPOT, i=1, nrad1)
            End Do
            Read (atomf,*) cwork
            Do is = 1, system%bns
               Read (atomf,*) cwork
               Read (atomf, 104) (ab%eny(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%ppc(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%ppd(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%ppq(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%ppp(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%dny(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%finy(il, is), il=1, opt%nl)
               Read (atomf, 104) (ab%finyd(il, is), il=1, opt%nl)
            End Do
            Close (atomf)
         End Do

         Write (IW6, 160) system%bef, system%bnum, system%bns
160      Format (/ 4 X, 'BULK SUBSTRATE:  EF=', G15.7, 8 X, 'NAB=', I2, 8 X, 'NSB=', I2)
         Write (IW6, 161)
161      Format (/ 4 X, 'COMPONENT,   CONC.,      Z,       WS,        WSAV: ')
         Do ia = 1, system%bnum
            ab => system%bat (ia)
            Write (IW6, 162) ia, ab%con, ab%az, ab%ws, ab%wsav
         End Do
162      Format (10 X, I2, F12.6, F9.3, 2 F12.6)

         If (system%vef /= system%bef .And. opt%ivac == 1) write (IW6, '(A50)') 'WARNING LEFT AND RIGHT FERMI&
        & LEVELS ARE DIFFERENT'


!!$ ------------------------ COMPLEX NODES AND WEIGHTS - UNIT IR13
         Read (IR13, 100) OWORK
         Write (IW6, 1001) OWORK

         Allocate (cnodes%zcn(opt%NE))
         Allocate (cnodes%zcw(opt%NE))
         Do j = 1, 8
            Read (IR13, 101) ne1
            Do ie = 1, ne1
               Read (IR13, 104) cnodes%zcn(ie), cnodes%zcw(ie)
            End Do
            If (ne1 == opt%NE) Exit
         End Do

         Write (IW6, 195)
195      Format (/ '     **  END OF READING ')

      End Subroutine rall

!*******************
!XXX    TISK    ****
!*******************
      Subroutine tisk (opt, system, iiter)
!****************************************
!!$  SHORT AND FULL PRINT OF THE RESULTS
!****************************************
         Use green_func
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Real (Kind=prec) :: QT (system%num), QP (opt%nl, system%num), AMP (opt%nl, system%num), BHFC &
        & (system%num), BHFV (system%num), BHFT (system%num), DNC (system%num), DNV (system%num), DNT &
        & (system%num)
         Real (Kind=prec), Parameter :: rcz = 0.0d0, rch = 0.5d0, DLIM = 0.02d0, BCONST = 52.4308d0
         Real (Kind=prec) :: qtot, amtot, etot, wrkf, pole, modpr
         Integer :: ia, is, il, j, ix, iiter
         Logical :: do_end
         Type (atom_definition), Pointer :: aa
         Integer :: precl
         Integer, Parameter :: gflst = 84643

!------------------------------ SHORT PRINT IN EACH ITERATION
!!$                               ATOMIC CHARGES AND MOMENTS
         Write (IW6, 110) iiter
110      Format (/ 1 X, '----- ATOM, CHARGE TRANSF.,', ' MAGN. MOMENT, DIPOLE MOMENT:        ITER=', I4)
         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 111) ia, aa%chatra, aa%amgmom, aa%dipmom
         End Do
111      Format (5 X, I5, 3 F15.5)

!!$     TOTAL CHARGE, MAGNETIC MOMENT, DIPOLE BARRIER, WORK FUNCTION, AND TOTAL ENERGY
         qtot = rcz
         amtot = rcz
         etot = rcz
         Do ia = 1, system%num
            aa => system%at (ia)
            qtot = qtot + aa%con * aa%chatra
            amtot = amtot + aa%con * aa%amgmom
            etot = etot + aa%con * aa%ACTE
         End Do

         If (opt%ivac == 0) wrkf = system%dba - system%bef
         If (opt%ivac == 1) wrkf = rcz

         Write (IW6, 116) qtot, amtot, iiter
116      Format (/ 2 X, 'TOTAL CHARGE=', F12.6, 5 X, 'TOTAL MAGMOM=', F11.5, 5 X, 'ITER=', I4)
         Write (IW6, 117) system%dba, wrkf, iiter
117      Format (/ 2 X, 'DIP. BARRIER=', F11.5, 5 X, 'WORK FUNCTION=', F11.5, 5 X, 'ITER=', I4)
         Write (IW6, 118) etot, iiter
118      Format (/ 2 X, 'TOTAL ENERGY=', G20.12, 20 X, 'ITER=', I4)

!------------------ POLES OF THE SCREENED POTENTIAL FUNCTIONS

127      Format (1 X, '  * POLE IN E= ', G15.7, '    FOR', '  IA=', I4, '  IS=', I1, '  IL=', I1)

         If (opt%ivac == 0) system%ef = system%bef
         If (opt%ivac == 1) system%ef = rch * (system%bef+system%vef)
         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%NS
               Do il = 1, opt%nl
                  pole = aa%ppc (il, is) - aa%ppd(il, is) / (aa%ppq(il, is)-system%qscr(il))
                  If (dabs(pole-system%ef) .Gt. DLIM) Cycle
                  Write (IW6, 127) pole, ia, is, il
               End Do
            End Do
         End Do

!---------------- FULL PRINT IN SELECTED ITERATIONS
         do_end = (opt%lnonc == 0 .And. opt%wtol <= opt%tol)

         modpr = Mod (iiter, opt%NFUPR)
         If (modpr .Ne. 0 .And. iiter .Ne. 1 .And. iiter .Ne. opt%NITER .And. ( .Not. do_end)) Return

         Write (IW6, 129) iiter
129      Format (/ 3 X, '**********  FULL PRINT IN ITERATION:', ' ITER=', I4, '  **********')

!!$       ATOMIC MADELUNG TERMS AND conTRIBUTIONS TO TOTAL ENERGY
         Write (IW6, 112) iiter
112      Format (/ 1 X, '----- ATOM, MADELUNG TERM, DIP. MAD. TERM,', '  ENERGY TERM:      ITER=', I4)
         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 113) ia, aa%VMAD, aa%DMAD, aa%ACTE
         End Do
113      Format (5 X, I5, 2 X, 2 G15.7, 2 X, G20.12)
!!$                                        POT. PARAMETERS C
         Write (IW6, 114) iiter
114      Format (/ 1 X, '----- ATOM, SPIN, POTENTIAL PARAMETERS C:', 15 X, 'ITER=', I4)
         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%NS
               Write (IW6, 115) ia, is, (aa%ppc(il, is), il=1, opt%nl)
            End Do
         End Do
115      Format (5 X, 2 I5, 4 F15.5)
!!$                                       CORE ENERGIES
         Write (IW6, 130) iiter
130      Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 13 X, '- CORES:   N,   L,  NOBC,', 8 X, 'EIGEN&
        &VALUES:' /)
131      Format (4 X, A16)
132      Format (20 X, 3 I5, 3 X, 2 G15.7)

         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 131) aa%otxta
            If (aa%numcor == 0) Cycle
            Do j = 1, aa%numcor
               Write (IW6, 132) aa%ncor(j), aa%lcor(j), aa%nobc(j), (aa%ecor(j, ix), ix=1, 2)
            End Do
         End Do
!!$              ENERGY MOMENTS - L-DIAGONAL

         Write (IW6, 135) iiter
135      Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 6 X, '- IS,  IL, ISHENY,   ENERGY MOMENTS (0,1&
        &,2):' /)
136      Format (4 X, A16)
137      Format (5 X, 3 I5, 3 X, 3 G15.7)
         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 136) aa%otxta
            Do is = 1, opt%NS
               Do il = 1, opt%nl
                  Write (IW6, 137) is, il, aa%ISHENY(il, is), aa%emdi0(il, is), aa%emdi1(il, is), &
                 & aa%emdi2(il, is)
               End Do
            End Do
         End Do
!!$                                       VALENCE CHARGES
         Do ia = 1, system%num
            aa => system%at (ia)
            QT (ia) = rcz
            Do il = 1, opt%nl
               QP (il, ia) = aa%emdi0(il, 1) + aa%emdi0(il, opt%NS)
               QT (ia) = QT (ia) + QP (il, ia)
            End Do
         End Do

         Write (IW6, 140) iiter
140      Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 5 X, '  - TOTAL  AND  PARTIAL VALENCE CHARGES:&
        &' /)
142      Format (4 X, A16)
143      Format (1 X, G15.7, 3 X, 4 G15.7)
         Do ia = 1, system%num
            aa => system%at (ia)
            Write (IW6, 142) aa%otxta
            Write (IW6, 143) QT (ia), (QP(il, ia), il=1, opt%nl)
         End Do

!!$                                       MAGNETIC MOMENTS
         If (opt%NS == 2) Then
            Do ia = 1, system%num
               aa => system%at (ia)
               Do il = 1, opt%nl
                  AMP (il, ia) = aa%emdi0(il, 1) - aa%emdi0(il, opt%NS)
               End Do
            End Do

            Write (IW6, 145) iiter
145         Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 5 X, '  - TOTAL  AND  PARTIAL MAGNETIC MOME&
           &NTS:' /)
147         Format (4 X, A16)
148         Format (1 X, G15.7, 3 X, 4 G15.7)
            Do ia = 1, system%num
               aa => system%at (ia)
               Write (IW6, 147) aa%otxta
               Write (IW6, 148) aa%amgmom, (AMP(il, ia), il=1, opt%nl)
            End Do
         End If

!!$                                       POTENTIAL PARAMETERS
         Write (IW6, 150) iiter
150      Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 15 X, '  ---- POTENTIAL PARAMETERS :')
1501     Format (/ 4 X, A16, 15 X, 'SPIN=', I2, 15 X, 'ITER=', I4 / 5 X, 'AZ=', F7.3, 8 X, 'WS=', F10.6, 8 X, &
        & 'WSAV=', F10.6)
151      Format (1 X, ' ENY : ', 4 G15.7)
152      Format (1 X, '  C  : ', 4 G15.7)
153      Format (1 X, 'DELTA: ', 4 G15.7)
154      Format (1 X, '  Q  : ', 4 G15.7)
155      Format (1 X, '  P  : ', 4 G15.7)
156      Format (1 X, ' DNY : ', 4 G15.7)
157      Format (1 X, 'FINY : ', 4 G15.7)
158      Format (1 X, 'FINYD: ', 4 G15.7)

         Do ia = 1, system%num
            aa => system%at (ia)
            Do is = 1, opt%NS
               Write (IW6, 1501) aa%otxta, is, iiter, aa%az, aa%ws, aa%wsav
               Write (IW6, 151) (aa%eny(il, is), il=1, opt%nl)
               Write (IW6, 152) (aa%ppc(il, is), il=1, opt%nl)
               Write (IW6, 153) (aa%ppd(il, is), il=1, opt%nl)
               Write (IW6, 154) (aa%ppq(il, is), il=1, opt%nl)
               Write (IW6, 155) (aa%ppp(il, is), il=1, opt%nl)
               Write (IW6, 156) (aa%dny(il, is), il=1, opt%nl)
               Write (IW6, 157) (aa%finy(il, is), il=1, opt%nl)
               Write (IW6, 158) (aa%finyd(il, is), il=1, opt%nl)
            End Do
         End Do

!!$                                       NUCLEAR QUANTITIES
         If (opt%IREL == 0) Then
!!$                                       HYPERFINE FIELDS
            If (opt%NS == 2) Then
               Do ia = 1, system%num
                  aa => system%at (ia)
                  BHFC (ia) = (aa%rhocor(1, 1)-aa%rhocor(1, opt%NS)) * BCONST
                  BHFV (ia) = (aa%rhoval(1, 1)-aa%rhoval(1, opt%NS)) * BCONST
                  BHFT (ia) = BHFC (ia) + BHFV (ia)
               End Do

170            Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 20 X, '- TOTAL, CORE, AND VALENCE HYPERF&
              &INE FIELDS (T):' /)
172            Format (1 X, A16, 4 X, G15.7, 3 X, G15.7, G15.7)

               Write (IW6, 170) iiter
               Do ia = 1, system%num
                  aa => system%at (ia)
                  Write (IW6, 172) aa%otxta, BHFT (ia), BHFC (ia), BHFV (ia)
               End Do
            End If
!!$                                       DENSITIES AT NUCLEI
            Do ia = 1, system%num
               aa => system%at (ia)
               DNC (ia) = aa%rhocor(1, 1) + aa%rhocor(1, opt%NS)
               DNV (ia) = aa%rhoval(1, 1) + aa%rhoval(1, opt%NS)
               DNT (ia) = DNC (ia) + DNV (ia)
            End Do

            Write (IW6, 175) iiter
175         Format (/ 1 X, ' ****  ATOMS, ', 45 X, ' ITER=', I4 / 20 X, '- TOTAL, CORE, AND VALENCE DENSITY A&
           &T NUCLEUS:' /)
177         Format (1 X, A16, 1 X, 2 G18.10, G15.7)

            Do ia = 1, system%num
               aa => system%at (ia)
               Write (IW6, 177) aa%otxta, DNT (ia), DNC (ia), DNV (ia)
            End Do
         End If

         Call write_atoms (opt, system)
         If (opt%rwgamma == 1) Then
            precl = size (zgfg) * prec * 2
            Open (Unit=gflst, File='rest1.dat', Action='write', Form='UNFORMATTED', Access='direct', &
           & Recl=precl)
            Write (gflst, Rec=1) zgfg
            Write (gflst, Rec=2) zcpf
            Write (gflst, Rec=3) zomg
            Close (Unit=gflst)
            precl = size (zcagf) * prec * 2
            Open (Unit=gflst, File='rest2.dat', Action='write', Form='UNFORMATTED', Access='direct', &
           & Recl=precl)
            Write (gflst, Rec=1) zcagf
            Close (Unit=gflst)
         End If

         If (do_end) Then
            Write (*,*) ' '
            Write (IW6,*) '  ***  Required accuracy reached at ITER=', iiter
            Write (IW6, 116) qtot, amtot, iiter
            Write (IW6, 117) system%dba, wrkf, iiter
            Write (IW6, 118) etot, iiter
            Write (*,*) ' '
            Write (IW6,*) '  *** Terminating XIES ***'
            Stop
         End If
         Return
      End Subroutine tisk

      Subroutine write_atoms (opt, system)
!!$  writes atoms in a new format
!!$  taken from Anton
         Use definitions
         Implicit None
         Type (params) :: opt
         Type (atoms_set) :: system
         Character (Len=16) :: cwork
         Integer :: ia, il, is, i, ix, atomf1
         Type (atom_definition), Pointer :: aa

992      Format (1 X, I5, 54 X, '# ', A10, ' IS=', I1)
993      Format (1 X, I5, 42 X, A16)
994      Format (1 X, 4 G15.7)
995      Format (1 X, G15.7, 42 X, A16)
!!$ 996      Format (1 X, E10.1, 44 X, A16)
99511 FORMAT(1X,G15.7,G15.7,29X,A20)
996   FORMAT(1X,E15.7,   39X,A16)
997      Format (1 X, 3 G15.7, 12 X, A22)
1101     Format (1 X, 10 I5)
1104     Format (1 X, 4 G15.7)
1190     Format (1 X, '  --------  LSDA-FILE  -------- ')
1191     Format (1 X, 49 X, A16, ' IS=', I1)
1192     Format (1 X, I5, 44 X, A16, ' IS=', I1)
1193     Format (1 X, I5, 44 X, A16)
1195     Format (1 X, '  ----------------------------- ')
!$omp parallel do default(private) shared(opt,system)
         Do ia = 1, system%num
            atomf1 = atomf + ia
            aa => system%at (ia)
            Read (aa%otxta,*) cwork
            Open (Unit=atomf1, File='atoms/'//trim(cwork), Action='write')
            Write (atomf1, '(A,A,A,I1,A,I5)') '#  Interface results for: ', trim (cwork), ', IVXC=', &
           & opt%IVXC, ', NK=', opt%bz%NK
            Write (atomf1, 993) opt%nl, '# NL'
            Write (atomf1, 994) (system%qscr(il), il=1, opt%nl)
            Write (atomf1, 997) aa%az, aa%ws, aa%wsav, '# AZ, WSR, WSR(avrg)'
            Write (atomf1, 99511) system%ef, aa%amgmom, '# Fermi.En, Mag.Mom'
            Write (atomf1, 996) 0.0d0, '# Pot.Shift'
            Write (atomf1, 993) 0, '# SW'
            Write (atomf1, 993) opt%NS, '# NS'
            Do is = 1, opt%NS
               Write (atomf1, 992) aa%nszrad, aa%otxta, is
               Write (atomf1, 994) (aa%pot(i, is), i=1, aa%nszrad)
            End Do
            Write (atomf1,*) 'XXXXXXXXXXXX CORE STATES XXXXXXXXXXXXXXXXXXXXXX'
            Do is = 1, opt%NS
               Write (atomf1, 1191) aa%otxta, is
               Write (atomf1, 1104) (aa%eny(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%ppc(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%ppd(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%ppq(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%ppp(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%dny(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%finy(il, is), il=1, opt%nl)
               Write (atomf1, 1104) (aa%finyd(il, is), il=1, opt%nl)
            End Do
            Write (atomf1, 1193) aa%numcor, aa%otxta
            If (aa%numcor /= 0) Then
               Do i = 1, aa%numcor
                  Write (atomf1, 1101) aa%ncor(i), aa%lcor(i), aa%nobc(i)
                  Write (atomf1, 1104) (aa%ecor(i, ix), ix=1, 2)
               End Do
            End If
            If (opt%ivac == 0) write (atomf1, 1104) system%dba
            Close (atomf1)
         End Do
!$omp end parallel do
         Return
      End Subroutine write_atoms
End Module read_write
