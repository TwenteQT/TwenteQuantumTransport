!!$ $Id: atoms.F90 1964 2013-06-28 09:01:11Z yuanz $
#include "math_def.h"

!!$
Module atoms_module
!!$    Procedures for reading atoms defenitions
   Use sparselib
   Use heritage
   Implicit None

   Type t_sk_param
      Integer :: nshell, ns
      Logical :: nohopping, defined
!!$ *_s, *_p, *_d - sigma, pi and delta S-K parameters for Hamiltonian
      Real(Kind=DEF_DBL_PREC), Pointer :: rshell(:), ss_s(:, :), sp_s(:, :), pp_s(:, :), pp_p(:, :), &
     & sd_s(:, :), pd_s(:, :), pd_p(:, :), dd_s(:, :), dd_p(:, :), dd_d(:, :), ps_s(:, :), ds_s(:, &
     & :), dp_p(:, :), dp_s(:, :)
   End Type t_sk_param

   Type t_potpar
      Real(Kind=DEF_DBL_PREC), Pointer :: c(:), d(:), q(:), p(:), dny(:), ovl(:), finy(:), finyd &
     & (:), pf(:), o(:), sqrtd(:), eny(:)
      Complex (Kind=DEF_DBL_PREC), Pointer :: bdotl(:, :) !this is an effective B-field on the atomic basis
      !functions. It needs to be a matrix because it is
      !not diagonal in the complex spherical harmonics
      Real(Kind=DEF_DBL_PREC), Pointer :: PRP(:)
   End Type t_potpar

   Type t_potpar_el
      Type(t_potpar), Pointer :: el
   End Type t_potpar_el

   Type t_dbl_ptr
      Real(Kind=DEF_DBL_PREC), Pointer :: val(:)
   End Type t_dbl_ptr

   Type t_atom_defenition
      Character(Len=20) :: label
      Integer :: nl, ns, label_len, nopol
      Real(Kind=DEF_DBL_PREC) :: az, pot_shift, Ef, magmom
      Real(Kind=DEF_DBL_PREC) :: bz_shift ! effective B-field used to orbitally polarize in the z-direction
      Integer :: nm, haveSO = 0
      Real(Kind=DEF_DBL_PREC) :: wsr, eny
      Complex(Kind=DEF_DBL_PREC) :: zny
      Type(rvec) :: pot(2)
      Real(Kind=DEF_DBL_PREC), Pointer :: alpha(:)
      Type(t_potpar_el) :: pp(2)
      Real(Kind=DEF_DBL_PREC), Pointer :: PvP(:, :), PvPD(:, :), PDvPD(:, :)
      Real(Kind=DEF_DBL_PREC), Pointer :: split(:)
      Real(Kind=DEF_DBL_PREC) :: socscaling
      Integer, Pointer :: idxdn(:)
!!$ Slater-Koster  data
      Integer :: index
      Real(Kind=DEF_DBL_PREC), Pointer :: s(:), p(:), d1(:), d2(:), d0(:), pd(:), opd(:)
      Real(Kind=DEF_DBL_PREC), Pointer :: sk_shift(:)
   End Type t_atom_defenition

   Type t_atoms_set
      Integer :: num, nlmax, nrmax, nsmax
      Type(t_atom_defenition), Pointer :: at(:)
!!$ Slater-Koster  data
      Type(t_sk_param), Pointer :: skh(:, :)
      Type(t_sk_param), Pointer :: sko(:, :)
   End Type t_atoms_set

   Type t_potparp_options
      Integer :: irel = 1, nsirk = 2
      Real(Kind=DEF_DBL_PREC) :: dE = 0.0001
      Real(Kind=DEF_DBL_PREC) :: ei
   End Type t_potparp_options

Contains
!>Read atom definitions from set of files.
!!Atomic files must be in "atoms" directory
!!@param filename name of file with atomlist
   Function read_atoms(filename) Result(atoms)
      Use logging
      Implicit None
      Type(t_atoms_set) :: atoms !< structure with atomic parameters returned as result
      Character(Len=*) :: filename
!!$ Local
      Integer :: numb
      Character(Len=100) :: cwork !< temporary string variable
      Character(Len=20) :: cwork1 !< temporary string variable
      Integer :: i, ins, nr, nat, nlmax = 0, nrmax = 0, nsmax = 0, swapspins, nl, ios
      Integer, Parameter :: listf = 100, atomf = 101
      Integer, Parameter :: socf = 102
      Logical            :: socsca
      Character(Len=100) :: cworksoc
      Type(t_atom_defenition), Pointer :: at
      Type(t_potpar), Pointer :: pp
      Type(rvec) :: rv

      Call do_log(1, 'Reading atom defenitions')
      Open (Unit=listf, File=filename, Action='read')
      Inquire (File='soc_mask', exist=socsca)
      If (socsca) Then
         Call do_log(1, 'SOC strength for each potential is scaled from soc_mask')
         Open (Unit=socf, File='soc_mask', Action='read')
      End If
100   Format(5 a16)
101   Format(1 x, 10 i5)
104   Format(1 x, 4 g15.7)
10414 Format(1 x, 10 g15.7)

!!$ Skip first two lines left for comments
      Read (listf, 100) cwork
      Read (listf, 100) cwork
!!$    read (listf, '(i)') numb
      Read (listf, *) numb

      atoms%num = numb
      Write (cwork, *) numb
      Call do_log(2, ' Number of atoms = '//trim(cwork))

      Allocate (atoms%at(numb))
      Do nat = 1, numb
         at => atoms%at(nat)
         Read (listf, *) cwork1, cwork
!!$       allocate(at%label(len_trim(cwork1)))
         at%label = trim(cwork1)
         at%label_len = len_trim(cwork1)
         Call do_log(2, ' Reading atom: '//at%label(1:10)//' from file: '//trim(cwork))
         If (socsca) Then
            Read (socf, *) at%socscaling
            write (cworksoc, '(a,f8.3)') 'SOC for '//trim(at%label)//' is scaled to be:', at%socscaling
            Call do_log(1, trim(cworksoc))
         Else
            at%socscaling = 1.d0
         End If

         Open (Unit=atomf, File='atoms/'//trim(cwork), Action='read')
         Read (atomf, 100) cwork
         Read (atomf, *) nl ! No. of partial waves
         at%nl = nl
         nlmax = Max(nlmax, nl)
         Allocate (at%alpha(nl))
         Read (atomf, 10414) at%alpha ! screening
         Read (atomf, 104) at%az, at%wsr ! At. number, WS radii
         Read (atomf, *, IoStat=ios) at%Ef, at%magmom ! Fermi Energy
         If (ios > 0) at%magmom = 0.0d0
!!$             If (ios < 0)
!!$       pshift=(ef-at%Ef)
         Read (atomf, *, IoStat=ios) at%pot_shift, at%bz_shift ! global potential shift, B.Lz potential shift
         If (ios > 0) at%bz_shift = 0d0
         swapspins = 0
         Read (atomf, *) swapspins ! Do we need swap spins?
         Read (atomf, *) at%ns ! No. of partial waves
         if (at%ns == 1) at%magmom = 0.0d0
         nsmax = Max(nsmax, at%ns)
!!$ Read potential
         Do ins = 1, at%ns
            Read (atomf, *) nr ! No. of radial points
            nrmax = Max(nrmax, nr)
            at%nm = nl**2
            at%pot(ins)%size = nr
            Allocate (at%pot(ins)%vec(nr))
            Read (atomf, 104) (at%pot(ins)%vec(i), i=1, nr)! potenti
            at%pot(ins)%vec = at%pot(ins)%vec + at%pot_shift
            Allocate (at%pp(ins)%el)
            pp => at%pp(ins)%el
            Allocate (pp%c(nl))
            Allocate (pp%d(nl))
            Allocate (pp%q(nl))
            Allocate (pp%p(nl))
            Allocate (pp%dny(nl))
            Allocate (pp%ovl(nl))
            Allocate (pp%bdotl(at%nm, at%nm))
            Allocate (pp%o(nl))
            Allocate (pp%sqrtd(nl))
            Allocate (pp%finy(nl))
            Allocate (pp%finyd(nl))
            Allocate (pp%pf(nl))
            Allocate (pp%PRP(nl))
         End Do
         Allocate (at%PvP(3, nl))
         Allocate (at%PvPD(3, nl))
         Allocate (at%PDvPD(3, nl))
         Allocate (at%split(nl))
!!$ Set potential for spin=2 to be equal potential for spin=1, if we have spinless input file
         If (at%ns < 2) Then
            at%pot(2) = at%pot(1)
            at%pp(2)%el => at%pp(1)%el
         End If
         If (swapspins /= 0) Then
            at%magmom = -at%magmom
            pp => at%pp(2)%el
            at%pp(2)%el => at%pp(1)%el
            at%pp(1)%el => pp
            rv = at%pot(1)
            at%pot(1) = at%pot(2)
            at%pot(2) = rv
         End If
         at%nopol = 0
         Close (atomf)
      End Do
      atoms%nrmax = nrmax
      atoms%nlmax = nlmax
      atoms%nsmax = nsmax
      Close (listf)
      Call do_log(1, 'Reading atom defenitions ....Done!')
      Call do_log(1, '')
      Return
   End Function read_atoms

   Subroutine rep_pfn(atoms, repf)
!!$    Writes the calculated P values to the "andout" file
      Use logging
      Implicit None
      Type(t_atoms_set) :: atoms
      Integer :: repf
!!$    Arguments
      Type(t_atom_defenition), Pointer :: at
!!$    Local variables
      Integer :: inp, s, l

      If (Log_Level < 3) Return
!!$ 38506    Format (12 x, 10(a1, G10.4, 2 x, G10.4, 2 x, G10.4, a1, 3 x))
38506 Format(12 x, 10(a1, G10.4, 2 x, G10.4, a1, 3 x))
      Write (repf, *) '---------------------------------------------------------------'
      Write (repf, *) 'Potential parameters report:'
      Do inp = 1, atoms%num
         at => atoms%at(inp)
         Do s = 1, at%ns
            Write (repf, '(a12,": sp=",i1,"  PFN=",10(1x,f12.6))') trim(at%label), s, &
           & at%pp(s)%el%pf(1:at%nl)
         End Do
         Write (repf, '("            alpha-A:")')
         Write (repf, 38506) ("|", at%PvP(1:2, l), "|", l=2, at%nl)
         Write (repf, '("            alpha-B:")')
         Write (repf, 38506) ("|", at%PvPD(1:2, l), "|", l=2, at%nl)
         Write (repf, '("            alpha-C:")')
         Write (repf, 38506) ("|", at%PDvPD(1:2, l), "|", l=2, at%nl)
         Write (repf, *)
      End Do
      Write (repf, *) '---------------------------------------------------------------'
   End Subroutine rep_pfn

   Subroutine calc_potpar(atoms, opts, eoffs)
!!$    Calculates potential function for the part of the system.
!!$    P is determined exactely (no linearization) by calling
!!$    an atomic solver and determining psi.
      Use logging
      Use sparselib
      Implicit None
!!$    Arguments
      Type(t_atoms_set) :: atoms
      Type(t_potparp_options) :: opts
      Real(Kind=DEF_DBL_PREC) :: eoffs
!!$    Local variables
!!$ "2/c^2"
      Real(Kind=DEF_DBL_PREC), Parameter :: scf0 = 2.663967179924343e-05
      Integer :: inp, nr, irel1, il, l, is, nsirk
      Real(Kind=DEF_DBL_PREC), Target :: fid(atoms%nrmax, atoms%nlmax, 2), fi(atoms%nrmax, atoms%nlmax, &
     & 2)
      Real(Kind=DEF_DBL_PREC), Pointer :: fip(:, :), fidp(:, :)
      Real(Kind=DEF_DBL_PREC) :: vi(opts%nsirk*2 - 1, atoms%nrmax), r(atoms%nrmax), deltaV(atoms%nrmax)
      Real(Kind=DEF_DBL_PREC), Target :: dV(atoms%nrmax, 2), dVm(atoms%nrmax)
      Real(Kind=DEF_DBL_PREC), Pointer :: v(:)
      Type(t_potpar), Pointer :: pp
      Real(Kind=DEF_DBL_PREC) :: amg, cfag
      Real(Kind=DEF_DBL_PREC) :: az1, ws1, wsav1
      Real(Kind=DEF_DBL_PREC), Parameter :: azl = 2.5d0
      Type(t_atom_defenition), Pointer :: at
      Type(t_dbl_ptr) :: dv_ptr(2)
      Call do_log(1, 'Calculating potential parameters for atom set...')

      nsirk = opts%nsirk

      dv_ptr(1)%val => dV(:, 1)

      Do inp = 1, atoms%num
         at => atoms%at(inp)
         az1 = at%az
         ws1 = at%wsr
         wsav1 = at%wsr

         Do is = 1, at%ns

            irel1 = Min(opts%irel, 1)
            If (az1 .Lt. azl) irel1 = 0
            nr = at%pot(is)%size
            v => at%pot(is)%vec
            pp => at%pp(is)%el
            Call rapo(nr, ws1, r)! Preparing the radial mesh
            Call lipo(az1, nsirk, nr, r, v, vi)! Lagrange interpolation for poten
!!$ We will need fi and fidot later, so ,we set pointers to array for appropriate spin

            fip => fi(:, :, is)
            fidp => fid(:, :, is)
            dv_ptr(2)%val => dV(:, is)
            Do il = 1, at%nl
               l = il - 1
               at%eny = eoffs + at%Ef

               Call rsel(az1, at%eny, opts%dE, wsav1, il, irel1, nsirk, nr, r, v, vi, fip, fidp, pp)! Atomic
!!$ Screen pot.fun. and calculate sqrt(\delta) and overlap
               pp%pf(il) = pp%pf(il)/(1.0d0 - at%alpha(il)*pp%pf(il))
               amg = (at%alpha(il) - pp%q(il))
               cfag = 1.0d0 + amg/pp%d(il)*(pp%c(il) - at%eny)
               pp%sqrtd(il) = cfag*Sqrt(pp%d(il))
               pp%ovl(il) = -amg/(cfag*pp%d(il))

               pp%PRP(il) = quad3(nr, r, fi(:, il, is)**2*r(:)**3)
            End Do

!!$ Generate B.L potential shift
            Call build_vshift(at%bz_shift, pp%bdotl, at%nm)
            Call convert_vshift_rylm(pp%bdotl, at%nm)

!!$ Calculate radial derivative of potential
            dv_ptr(is)%val(1) = 0
            Call dlangr(5, nr - 1, r(2:nr), v(2:nr), dv_ptr(is)%val(2:nr))
!// dV -> dV/dR*R
            dv_ptr(is)%val(1:nr) = dv_ptr(is)%val(1:nr)*r(1:nr)

         End Do

!!$             write(*,*) inp
!!$             write(*,*) pp%d(:)
!!$             write(*,*) pp%sqrtd (:)
!!$              write(*,*) pp%ovl (:)
!!$             write(*,*) pp%pf (:)

         If (at%ns > 1) Then
            dVm(1:nr) = 0.5d0*(dv_ptr(1)%val(1:nr) + dv_ptr(2)%val(1:nr))
            deltaV(1:nr) = at%pot(2)%vec(1:nr) - at%pot(1)%vec(1:nr)
         End If
!!$    //  Calculate expectation values for 1/R*dV/dR
         Do il = 1, at%nl
!!$    // Calculate diagonal(in spin-space) elements
            Do is = 1, at%ns
!// Int(phi^2*dV/dR*R)
               at%PvP(is, il) = scf0*quad3(nr, r, dv_ptr(is)%val(1:nr)*fi(1:nr, il, is)*fi(1:nr, il, &
              & is))
!// Int(phi_dot*dV/dR*R*phi)
               at%PvPD(is, il) = scf0*quad3(nr, r, dv_ptr(is)%val(1:nr)*fi(1:nr, il, is)*fid(1:nr, il, &
              & is))
!// Int(phi_dot^2*dV/dR*R)
               at%PDvPD(is, il) = scf0*quad3(nr, r, dv_ptr(is)%val(1:nr)*fid(1:nr, il, is)*fid(1:nr, &
              & il, is))
            End Do

            If (at%ns < 2) Then
               at%PvP(2, il) = at%PvP(1, il)
               at%PvPD(2, il) = at%PvPD(1, il)
               at%PDvPD(2, il) = at%PDvPD(1, il)
               at%split(il) = 0.0d0
            Else
               at%split(il) = 0.5d0*quad3(nr, r, ((fi(1:nr, il, 1)**2) + (fi(1:nr, il, &
              & 2)**2))*(r(1:nr)**2)*deltaV(1:nr))
            End If

            at%PvP(3, il) = .5d0*(at%PvP(2, il) + at%PvP(1, il))
            at%PvPD(3, il) = .5d0*(at%PvPD(2, il) + at%PvPD(1, il))
            at%PDvPD(3, il) = .5d0*(at%PDvPD(2, il) + at%PDvPD(1, il))

!!$                If (at%ns > 1) Then
!!$ !!$            // Calculate off-diagonal(in spin-space) elements
!!$ !// Int(phi*dV/dR*R*phi_prime)
!!$                   at%PvP (3, il) = scf0 * quad3 (nr, r, dVm(1:nr)*fi(1:nr, il, 1)*fi(1:nr, il, 2))
!!$ !// Int(phi_dot*dV/dR*R*phi_dot_prime)
!!$                   at%PDvPD (3, il) = scf0 * quad3 (nr, r, dVm(1:nr)*fid(1:nr, il, 1)*fid(1:nr, il, 2))
!!$ !// Int(phi*dV/dR*R*phi_dot_prime)
!!$                   at%PvPD (3, il) = scf0 * quad3 (nr, r, .5d0*dVm(1:nr)*(fi(1:nr, il, 1)*fid(1:nr, il, &
!!$                  & 2)+fi(1:nr, il, 2)*fid(1:nr, il, 1)))
!!$                End If

         End Do
!!$         SOC scaling for each atom potential according to soc_mask
         at%PvP = at%PvP*at%socscaling
         at%PvPD = at%PvPD*at%socscaling
         at%PDvPD = at%PDvPD*at%socscaling

!!$             If (at%ns < 2) Then
!!$                at%PvP (2, :) = at%PvP(1,:)
!!$                at%PvP (3, :) = at%PvP(1,:)
!!$                at%PvPD (2, :) = at%PvPD(1, :)
!!$                at%PvPD (3, :) = at%PvPD(1, :)
!!$                at%PDvPD (2, :) = at%PDvPD(1, :)
!!$                at%PDvPD (3, :) = at%PDvPD(1, :)
!!$             End If
         If (maxval(dabs(at%pp(1)%el%pf(1:at%nl) - at%pp(2)%el%pf(1:at%nl))) < 1e-15) at%nopol = 1
      End Do
      Call rep_pfn(atoms, log_stream_number)
      Call do_log(1, 'Calculating potential parameters for atom set...Done!')
      Call do_log(1, '')
   End Subroutine calc_potpar

   Subroutine build_vshift(shift, vshift, n)
      Implicit none
      Integer :: n
      Real(kind=DEF_DBL_PREC) :: shift
      Complex(kind=DEF_DBL_PREC) :: vshift(n, n)
!!$ local
      Integer :: il, lmx, offset, M

      lmx = sqrt(real(int(n))) - 1
      vshift = 0d0
      Do il = 1,lmx
        offset = il**2 + il + 1
        ! To build Lz
        Forall (M=-il:il) vshift(offset+M, offset+M) = shift*cmplx(M, 0.0d0, kind=DEF_DBL_PREC)
        ! To build Lx
        !Forall (M=-il  :il-1) vshift(offset+M, offset+M+1) = sqrt(1.*il*(il+1)-(M+1)*(M))*shift*cmplx(0.5d0, 0.0d0, kind=DEF_DBL_PREC)
        !Forall (M=-il+1:il  ) vshift(offset+M, offset+M-1) = sqrt(1.*il*(il+1)-(M-1)*(M))*shift*cmplx(0.5d0, 0.0d0, kind=DEF_DBL_PREC)
        ! To build Ly
        !Forall (M=-il  :il-1) vshift(offset+M, offset+M+1) = sqrt(1.*il*(il+1)-(M+1)*(M))*shift*cmplx(0.0d0, 0.5d0, kind=DEF_DBL_PREC)
        !Forall (M=-il+1:il  ) vshift(offset+M, offset+M-1) = sqrt(1.*il*(il+1)-(M-1)*(M))*shift*cmplx(0.0d0,-0.5d0, kind=DEF_DBL_PREC)
      Enddo
   End Subroutine build_vshift

   Subroutine convert_vshift_rylm(vshift, n)
      Implicit None
      Integer :: n
      Complex(Kind=DEF_DBL_PREC) :: vshift(n, n)
! local
      Complex(Kind=DEF_DBL_PREC) :: R(n, n), Ri(n, n)
      Complex(Kind=DEF_DBL_PREC) :: new_vshift(n, n)
      Integer :: l, nlm, lmax

      new_vshift = 0d0

      lmax = sqrt(real(int(n))) - 1

      Do l = 1, lmax
         nlm = 2*(l + 1) - 1
         call transmat(R(1:nlm, 1:nlm), l)
         Ri(1:nlm, 1:nlm) = transpose(conjg(R(1:nlm, 1:nlm)))
         new_vshift(l**2 + 1:(l + 1)**2, l**2 + 1:(l + 1)**2) = 0.5d0*matmul(R(1:nlm, 1:nlm), &
 & matmul(vshift(l**2 + 1:(l + 1)**2, l**2 + 1:(l + 1)**2), Ri(1:nlm, 1:nlm)))
      End do

      vshift = new_vshift

   Contains
      Subroutine transmat(R, l)
!!$ Transformation matrix complex->real harmonics (without 1/sqrt(2) prefactor!!!)
         Implicit None
         Integer :: l
         Complex(Kind=DEF_DBL_PREC) :: R(-l:l, -l:l)
!!$
         Integer :: M
         R(:, :) = DEF_cmplx_zero
         R(0, 0) = DEF_cmplx_one*Sqrt(2.0d0)
         Do M = 1, l, 1
            R(-M, -M) = DEF_cmplx_Ione
            R(M, M) = (-1.0d0)**M
            R(M, -M) = DEF_cmplx_one
            R(-M, M) = -DEF_cmplx_Ione*(-1.0d0)**M
         End Do
      End Subroutine transmat

   End Subroutine convert_vshift_rylm

   Subroutine dlangr(ord, n, x, y, yn)
      Implicit None
      Integer :: ord, n
      Real(Kind=DEF_DBL_PREC) :: x(:), y(:), yn(:)
!!$   Local vars
      Integer :: n1, n2, i
      n1 = ord/2
      n2 = n - (ord - n1)
      Do i = 1, n1
         yn(i) = auxdlangr(ord, 0, x, y, x(i))
      End Do
      Do i = n1 + 1, n2
         yn(i) = auxdlangr(ord, i - n1 - 1, x, y, x(i))
      End Do
      Do i = n2 + 1, n
         yn(i) = auxdlangr(ord, n2 - n1 - 1, x, y, x(i))
      End Do

      Return
   Contains

      Function auxdlangr(order, offs, x, y, xn) Result(yd)
         Implicit None
         Integer :: order, offs
         Real(Kind=DEF_DBL_PREC) :: x(:), y(:), xn, yd
!!$   Local vars
         Real(Kind=DEF_DBL_PREC) :: t1, t2, t3
         Integer :: i, j, k
         yd = 0.0d0
         Do i = 1 + offs, order + offs
            t1 = 1
            Do j = 1 + offs, order + offs
               If (j /= i) t1 = t1*(x(i) - x(j))
            End Do
            t1 = y(i)/t1
            t2 = 0
            Do k = 1 + offs, order + offs
               If (k /= i) Then
                  t3 = 1
                  Do j = 1 + offs, order + offs
                     If ((j /= k) .And. (j /= i)) t3 = t3*(xn - x(j))
                  End Do
                  t2 = t2 + t3
               End If
            End Do
            yd = yd + t1*t2
         End Do
      End Function auxdlangr
   End Subroutine dlangr

!!$*******************
!!$XXX    RSEL    ****
!!$*******************
   Subroutine rsel(az, eny, eh, wsa, il, irel, nsirk, nr, r, v, vi, fi, fid, pp)
!!$---------------------------------------------------------
!!$     LINEARIZATION OF ENERGY DEPENDENCE
!!$     OF SOLUTION OF RADIAL SCHROEDINGER EQUATION
!!$     FOR VALENCE ELECTRONS:
!!$      BOTH  NON-RELATIVISTIC (IREL=0) AND
!!$      SCALAR-RELATIVISTIC (IREL=1) VERSION
!!$---------------------------------------------------------
!!$  INPUT:
!!$     AZ - ATOMIC NUMBER
!!$     EH - ENERGY STEP FOR NUMERICAL DERIVATIVE
!!$     WSA - AVERAGE WIGNER-SEITZ RADIUS
!!$           (TO SCALE THE POTENTIAL PARAMETERS)
!!$     L - ORBITAL QUANTUM NUMBER
!!$     IREL - RELATIVITY
!!$     NR - SIZE OF RADIAL MESH
!!$     R(.) - RADIAL MESH
!!$     V(.) - POTENTIAL
!!$     VI(.,.) - INTERPOLATED POTENTIAL
!!$     ENY - ENERGY VALUE
!!$  OUTPUT:
!!$     PPC, ..., FINYD - POTENTIAL PARAMETERS
!!$     FI(.),FID(.),FIDD(.) - WAVEFUNCTION (NORMALIZED TO
!!$                  UNITY) AND ITS TWO ENERGY DERIVATIVES
!!$---------------------------------------------------------
      Implicit None
      Type(t_potpar), Pointer :: pp
      Real(Kind=DEF_DBL_PREC) :: eny, r(:), v(:), vi(:, :), fi(:, :), fid(:, :), wsa, eh, az
      Integer :: nsirk, nr, il, irel

!!$Local
      Real(Kind=DEF_DBL_PREC) :: ws, al, alp1, fak, dny, dnyd, wf(nr)
      Integer :: l, i, ie
      Real(Kind=DEF_DBL_PREC) :: wb(nr, -2:2), finy, finyd
      Real(Kind=DEF_DBL_PREC) :: h1, r1, r2, ajm
!!$ variables bellow used  for "fidd" calculation
!!$      real (kind=DEF_DBL_PREC) :: fidd(nr),a0,s1,s2,h2,rc16=16.0D0,rc30-30.0d0

      Real(Kind=DEF_DBL_PREC) :: rc1, rc2, rc8, rc12, e

      Data rc1/1.0D0/, rc2/2.0D0/, rc8/8.0D0/, rc12/12.0D0/

      l = il - 1
      al = real(l, kind=DEF_DBL_PREC)
      alp1 = al + rc1
      ws = r(nr)
      fak = (ws/wsa)**(2*l + 1)
      Do ie = -2, 2
         e = eny + real(ie)*eh
         Call rsev(az, e, eny, l, irel, nsirk, nr, r, v, vi, wf, wb(:, ie))
         If (ie .Eq. 0) dny = ws*wf(nr)/wb(nr, 0)
      End Do

      pp%pf(il) = rc2*real(2*l + 1, kind=DEF_DBL_PREC)*fak*(dny + al + rc1)/(dny - al)
      pp%dny(il) = dny

      h1 = rc1/(rc12*eh)
!!$ This is commented because we dont need fidd
!!$    h2 = h1 / eh

      fi(1:nr, il) = wb(1:nr, 0)

      Do i = 1, nr
         r1 = wb(i, 1) - wb(i, -1)
         r2 = wb(i, 2) - wb(i, -2)
         fid(i, il) = (rc8*r1 - r2)*h1
!!$ This is commented because we dont need fidd
!!$       a0 = wb (i, 0)
!!$       s1 = wb (i,-1) + wb (i, 1)
!!$       s2 = wb (i,-2) + wb (i, 2)
!!$         fidd (i) = (rc16*s1-s2-rc30*a0) * h2
      End Do

      finy = fi(nr, il)
      finyd = fid(nr, il)
      pp%finy(il) = finy
      pp%finyd(il) = finyd
      pp%o(il) = quad3(nr, r, fi(:, il)*fid(:, il))

      dnyd = dny - rc1/(ws*finy*finyd)
      ajm = dnyd + alp1

      pp%c(il) = eny - finy*(dny + alp1)/(finyd*ajm)
      pp%d(il) = fak/(rc2*ws*(finyd*ajm)**2)
      pp%q(il) = fak*(dnyd - al)/(rc2*(al + alp1)*ajm)
      pp%p(il) = quad3(nr, r, (r(1:nr)*fid(1:nr, il))**2)
!!$  write(*,*) pp%c (il),pp%d (il),pp%q (il),pp%p (il),pp%o (il),finy,finyd
      Return
   End Subroutine rsel

End Module atoms_module
