!   $Id: definitions.f,v 1.39 2003/11/03 19:04:43 maciej Exp $
module pot_mod
  use definitions   
  use heritage
  use ies_structure
contains

!*******************
!XXX    SLDA    ****
!*******************
  subroutine slda(opt,system,mixing)
!**********************************************
!   PREPARES THE LDA-ITERATIONS:
!   - DEFINES THE WEIGHTS FOR SCALAR PRODUCT
!   - CALCULATES THE VALENCE WAVE FUNCTIONS
!   - RECALCULATES THE POTENTIAL PARAMETERS
!**********************************************
   implicit none
   type (params) :: opt 
   type (atoms_set) :: system
   type (andmix):: mixing
   real (kind=prec),target :: fi(system%nrmax,opt%nl,2),fid(system%nrmax,opt%nl,2),fidd(system%nrmax,opt%nl,2)
   real (kind=prec), pointer :: fip(:,:),fidp(:,:),fiddp(:,:),v(:)
   real(kind=prec),parameter :: rcz=0.0D0,rc1=1.0D0,rc3=3.0D0,CLIM=0.01D0,PLIM=0.01D0,AZL=2.5D0,DMUL=0.03D0
   real(kind=prec) :: pref,vi(9,system%nrmax),r(system%nrmax)
   real (kind=prec):: ws1,wsav1,eny1,arg,az1,de
   integer :: i,j,is,nr1,irel1,ia,il,nam
   type(atom_definition), pointer :: aa

!     allocation of arrays needed for mixings
      allocate(mixing%SP(opt%NAM))
      allocate(mixing%DXP(opt%NAM,opt%NUH))
      allocate(mixing%DFP(opt%NAM,opt%NUH))
      allocate(mixing%XL(opt%NAM))
      allocate(mixing%FL(opt%NAM))
      allocate(mixing%XN(opt%NAM))
      allocate(mixing%VOMA(opt%NUH,opt%NUH))
!     end of allocation

   mixing%dxp  = rcz
   mixing%dfp  = rcz
   mixing%voma = rcz
   vi   = rcz
   
   j = 0 !                            WEIGHTS FOR SCALAR PRODUCT we need sp for mixing 
   do ia = 1, system%num
     aa=>system%at(ia)
     nr1 = aa%nszrad
     ws1 = aa%ws
     call rapo(nr1,ws1,r)
     pref = max (clim,aa%con)     !CC      PREF=rc1
     do is = 1, opt%ns
       do i = 2, nr1
         j = j + 1
         mixing%sp(j) = pref * r(i)**2 * (r(i) - r(i-1))
       enddo
     enddo
   enddo

   if (opt%ivac==0) then
     j = j+1
     mixing%sp(j) = real(opt%ns,kind=prec) * system%vwst**3 / rc3
   endif
   nam=j

   write(IW6,150) nam
150 format(//5X,'*** LENGTH OF VECTORS FOR ANDERSON', ' MIXING= ',I8)

   do ia = 1, system%num           !     LOOP OVER ATOMS
     aa=>system%at(ia)
     az1 = aa%az
     nr1 = aa%nszrad
     irel1 = min(opt%irel,1)
     if (az1.LT.azl) irel1 = 0
     ws1 = aa%ws
     wsav1 = aa%wsav
     do is = 1, opt%ns             !     LOOP OVER SPINS
      v=>aa%pot(:,is) 
      call rapo (nr1,ws1,r)        !  Preparing the radial mesh
      call lipo (az1,opt%nsirk,nr1,r,v,vi)   !  Lagrange interpolation for potentials
      fip=>fi(:,:,is)
      fidp=>fid(:,:,is)
      fiddp=>fidd(:,:,is)
      do il = 1, opt%nl            !     LOOP OVER L
        eny1 = aa%eny(il,is)
        arg  = MAX(aa%ppp(il,is),plim)
        de   = dmul/dSQRT(ARG)
        call rsel_t (az1,eny1,de,wsav1,il,is,irel1,nr1,r,v,vi,fip,fidp,fiddp,aa,opt%nsirk)
      enddo
     enddo
   enddo 

!---------------------------------------    PRINT
      write(IW6,160)
160   format(//5X,' ******  INPUT POTENTIAL PARAMETERS  : ')
161   format(/4X,'IA=',I4,8X,A16,8X,'IS=',I2/ 5X,'AZ=',F7.3,8X,'WS=',F10.6,8X,'WSAV=',F10.6)
162   format(1X,' ENY : ',4G15.7)
163   format(1X,'  C  : ',4G15.7)
164   format(1X,'DELTA: ',4G15.7)
165   format(1X,'  Q  : ',4G15.7)
166   format(1X,'  P  : ',4G15.7)
167   format(1X,' DNY : ',4G15.7)
168   format(1X,'FINY : ',4G15.7)
169   format(1X,'FINYD: ',4G15.7)
      do ia = 1, system%num
        aa=>system%at(ia)
        do is = 1, opt%ns
          write(iw6,161)ia,aa%otxta,is,aa%az,aa%ws,aa%wsav
          write(iw6,162)(aa%eny(il,is),il = 1,opt%nl)
          write(iw6,163)(aa%ppc(il,is),il = 1,opt%nl)
          write(iw6,164)(aa%ppd(il,is),il = 1,opt%nl)
          write(iw6,165)(aa%ppq(il,is),il = 1,opt%nl)
          write(iw6,166)(aa%ppp(il,is),il = 1,opt%nl)
          write(iw6,167)(aa%dny(il,is),il = 1,opt%nl)
          write(iw6,168)(aa%finy(il,is),il = 1,opt%nl)
          write(iw6,169)(aa%finyd(il,is),il = 1,opt%nl)
        enddo
      enddo
   return
 end subroutine slda

!!$*******************
!!$XXX    RSEL    ****
!!$*******************
 subroutine rsel_t (az,eny,eh,wsa,il,is,irel,nr,r,v,vi,fi,fid,fidd,aa,nsirk)
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
    implicit none
    type (atom_definition),pointer:: aa
    real (kind=prec) :: eny,r(:),v(:),vi(:,:),fi(:,:),fid(:,:),fidd(:,:),wsa, eh, az
    integer::nr,il,is,irel

!!$Local
    real (kind=prec) :: ws, al, alp1, fak,  dny, dnyd,wf(nr)
    integer :: l, i, ie, nsirk
    real (kind=prec) :: wb (nr,-2:2), finy, finyd
    real (kind=prec) :: h1, h2, r1, r2, ajm
!!$ variables bellow used  for "fidd" calculation
    real (kind=prec) :: a0,s1,s2,rc16=16.0D0,rc30=30.0d0
    real (kind=prec) :: rc1, rc2, rc8, rc12, e

    data rc1 / 1.0D0 /, rc2 / 2.0D0 /, rc8 / 8.0D0 /, rc12 / 12.0D0 /

    l = il - 1
    al = real(l,kind=prec)
    alp1 = al+rc1
    ws = r(nr)
    fak = (ws/wsa) ** (2*l+1)
    do ie = - 2, 2
     e = eny + real (ie, kind=prec) * eh
     call rsev (az, e, eny, l, irel,nsirk, nr,r,v,vi,wf,wb(:,ie))
     if (ie .eq. 0) dny = ws * wf (nr) / wb (nr,0)
    end do
    
    aa%dny (il,is) = dny

    h1 = rc12 * eh  !                                     phi, phi-DOT, ...
    h2 = h1 * eh

    fi (1:nr,il)=wb(1:nr,0)
    do i = 1, nr
       r1 = wb (i, 1) - wb (i,-1)
       r2 = wb (i, 2) - wb (i,-2)
       fid (i,il) = (rc8*r1-r2) / h1
       a0 = wb (i, 0)       
       s1 = wb (i,-1) + wb (i, 1)
       s2 = wb (i,-2) + wb (i, 2)       
       fidd (i,il) = (rc16*s1-s2-rc30*a0) / h2
    end do


    finy = fi (nr,il)
    finyd = fid (nr,il)
    dnyd = dny-rc1 / (ws*finy*finyd)

    aa%finy(il,is) = finy
    aa%finyd(il,is) = finyd

    ajm = dnyd + alp1

    aa%ppc (il,is) = eny - finy * (dny+alp1) / (finyd*ajm)
    aa%ppd (il,is) = fak / (rc2*ws*(finyd*ajm)**2)
    aa%ppq (il,is) = fak * (dnyd-al) / (rc2*(al+alp1)*ajm)   
    aa%ppp (il,is) = quad3 (nr, r, (r(1:nr)*fid(1:nr, il)) ** 2)

    aa%phi (1:nr,il,is)   = fi (1:nr,il)
    aa%phid (1:nr,il,is)  = fid (1:nr,il)
    aa%phidd (1:nr,il,is) = fidd (1:nr,il)

   return
 end subroutine rsel_t


!*******************
!XXX    DIMO    ****
!*******************
 subroutine dimo(opt,system,harmon)
!***********************************************************************
!   CALCULATES THE DIPOLE MOMENTS FROM THE L-OFF-DIAGONAL MOMENTS OF DOS
!***********************************************************************
  implicit none
  type(params):: opt
  type(atoms_set)::system
  type(sharm)::harmon
  real(kind=prec):: UINT(opt%pair),DIPMS(opt%NS),RCUB(system%nrmax),QQ(system%nrmax)
  real(kind=prec),parameter :: rcz=0.0d0, rc1=1.0d0, rc2=2.0d0, rc4=4.0d0, rch=0.5d0, rc3=3.0d0
  real(kind=prec):: pi,pref,arg,wm00,wm10,wm01,wm11,hm20,hm02,sum,rint,ws1,r(system%nrmax)
  integer:: ia,il,is,il1,il2,ipair,l1,l2,i1,i2,m,icen1,icen2,nr1,i
  type(atom_definition),pointer :: aa 

!$omp parallel default(private)&
!$omp& shared(system,opt,harmon,IW6)
   pi=rc4*dATAN(rc1)
!   ANGULAR INTEGRALS INCLUDING FACTOR of 2
   PREF = rc2 * dSQRT( rc4 * pi / rc3 )
   ipair=0
   do IL1=1,opt%nl-1
     IL2=IL1+1
     L1=IL1-1
     L2=IL1
     ICEN1=L1**2+L1+1
     ICEN2=L2**2+L2+1
     do M=-L1,L1
      ipair=ipair+1
      I1=ICEN1+M
      I2=ICEN2+M
      UINT(ipair)=PREF*harmon%gfrh(I2,I1,3)
     enddo
    enddo

!                           SIGNS FOR SQRT(PDOT) FUNCTIONS
!$omp do
   do ia = 1, system%num
     aa => system%at(ia)
     do is = 1, opt%ns
       do il = 1, opt%nl
         ARG = rc1/aa%FINY(il,is) - aa%ws * aa%FINYD(il,is)*(aa%DNY(il,is)+real(il,kind=prec))
         aa%ynam(il,is) = SIGN(rc1,ARG)
       enddo
     enddo
   enddo
!$omp end do

!$omp do
   do ia = 1, system%num!                        LOOP OVER ATOMS
     aa => system%at(ia)
     NR1 = aa%nszrad    !                        RADIAL MESH
     WS1 = aa%WS
     call RAPO(NR1,WS1,R)
     do i = 1, nr1
       RCUB(i) = R(i)**3
     enddo
     do is = 1, opt%ns !                         LOOP OVER SPINS
       DIPMS(is) = rcz
       ipair = 0
       do IL1=1,opt%nl-1 !                       LOOP OVER (L1,L2)-opt%pairS
         IL2=IL1+1
         L1=IL1-1
         do M=-L1,L1
          ipair = ipair + 1

          WM00 = aa%emof00(ipair,is)
          WM10 = aa%emof10(ipair,is)
          WM01 = aa%emof01(ipair,is)
          WM11 = aa%emof11(ipair,is)
          HM20 = rch * aa%emof20(ipair,is)
          HM02 = rch * aa%emof02(ipair,is)

          do i = 1, NR1
           sum = aa%phi(i,IL1,is) * aa%phi(i,IL2,is) * wm00 + aa%phid(i,IL1,is)  * aa%phid(i,IL2,is) * wm11 &
           &   + aa%phi(i,IL1,is) * aa%phid(i,IL2,is) * wm01 + aa%phid(i,IL1,is)  * aa%phi(i,IL2,is)  * wm10 &
           &   + aa%phi(i,IL1,is) * aa%phidd(i,IL2,is) * HM02 + aa%phidd(i,IL1,is) * aa%phi(i,IL2,is)  * HM20
           QQ(i) = RCUB(i) * sum
          enddo

          RINT = QUAD3(NR1,R,QQ)
          DIPMS(is) = DIPMS(is) + RINT * UINT(ipair) * aa%YNAM(IL1,is) * aa%YNAM(IL2,is)
         enddo
        enddo
       enddo
       aa%dipmom =  DIPMS(1) + DIPMS(opt%ns)
      enddo
!$omp end do
!$omp end parallel
   return
 end subroutine dimo


!*******************
!XXX    NEPO    ****
!*******************
 subroutine nepo(opt,system,madcon,mixing,nam,iiter)
!**********************************************************************
! CALCULATES NEW ONE-ELECTRON POTENTIALS AND TOTAL ENERGY CONTRIBUTIONS
!**********************************************************************
  implicit none
  type (params) :: opt 
  type (atoms_set) :: system
  type (madelung) :: madcon
  type (andmix):: mixing
  real(kind=prec):: YMAX(2,system%num),CHTRG(system%nums),DIPMG(system%nums),VMADG(system%nums), &
     &  DMADG(system%nums),WDUM(system%nrmax),RHOT(system%nrmax),RHOS(system%nrmax,opt%ns),VXC(2),VHART(system%nrmax)
  real(kind=prec),parameter :: rcz=0.0d0, rc1=1.0d0, rc2=2.0d0, rc4=4.0d0, rch=0.5d0, AZL=2.5d0, DMUL=0.03d0
  real (kind=prec),target :: fi(system%nrmax,opt%nl,2),fid(system%nrmax,opt%nl,2),fidd(system%nrmax,opt%nl,2)
  real (kind=prec),pointer :: fip(:,:),fidp(:,:),fiddp(:,:),v(:)
  real(kind=prec):: pi,pi4,az1,ws1,RHOUP,RHODO,vnuc,vcoul,dum,scor,amul1,amul2,sval,cote,elst,fdba,exc,rdvr,s2, &
     &  tl1,tl5,wsav1,r(system%nrmax),vi(9,system%nrmax),de,eny1
  integer:: ia,il,ic,ig,jg,is,i,j,nr1,iadd,l2p1,nuh1,ih,jh,l,irel1,ish1,IMAX(2,system%num),nam,iiter
  type(atom_definition),pointer ::aa

!$omp parallel default(private)&
!$omp& shared(opt,system,madcon,mixing,nam,iiter,CHTRG,DIPMG,VMADG,DMADG,FDBA,YMAX,IMAX,IW6)

!---------------------------------------------------------------------------------
!    WITH ANDERSON MIXING FOR POTENTIALS (V.EYERT, J.COMPUT.PHYS. 124 (1996) 271)
!---------------------------------------------------------------------------------
   pi = rc4*datan(rc1)
   pi4 = rc4*pi
   

!---------------------------- INTERATOMIC MADELUNG TERMS
!            ATOMIC CHARGE TRANSFERS AND MAGNETIC MOMENTS
!$omp do
    do ia = 1, system%num
     aa => system%at(ia) 
     aa%chatra = -aa%valz
     aa%amgmom = rcz
     do il = 1, opt%nl
      aa%chatra = aa%chatra + aa%EMDI0(il,1) + aa%EMDI0(il,opt%ns)
      aa%amgmom = aa%amgmom + aa%EMDI0(il,1) - aa%EMDI0(il,opt%ns)
     enddo
    enddo
!$omp end do
!        SITE-RESOLVED CHARGE TRANSFERS AND DIPOLE MOMENTS
    ia=0
!   !$omp do shared(CHTRG,DIPMG)

!$omp sections
    do ig = 1, system%nums
      CHTRG(ig)=rcz
      DIPMG(ig)=rcz
      do IC=1,system%nc(ig)
        ia = ia + 1
        aa => system%at(ia)
        CHTRG(ig)=CHTRG(ig) + aa%CON * aa%chatra
        DIPMG(ig) = DIPMG(ig) + aa%CON * aa%DIPMOM
      enddo
    enddo
!   !$omp end do

!                           SITE-RESOLVED MADELUNG TERMS

      VMADG(1:system%nums)=rcz
      DMADG(1:system%nums)=rcz
!$omp end sections

!$omp do reduction(+:VMADG,DMADG)
    do jg = 1, system%nums
      do ig = 1, system%nums
        VMADG(ig) = VMADG(ig) + madcon%mmm(ig,jg)*CHTRG(jg) + madcon%mmd(ig,jg)*DIPMG(jg)
        DMADG(ig) = DMADG(ig) + madcon%mdm(ig,jg)*CHTRG(jg) + madcon%mdd(ig,jg)*DIPMG(jg)
      enddo
    enddo
!$omp end do
!                                   ATOMIC MADELUNG TERMS
!$omp sections
    ia = 0
    do ig = 1, system%nums
      do IC = 1, system%nc(ig)
       ia = ia + 1
       aa=> system%at(ia)
       aa%VMAD = VMADG(ig)
       aa%DMAD = DMADG(ig)
      enddo
     enddo
!$omp end sections
!-------------  OUT-POTENTIALS AND ENERGY CONTRIBUTIONS

!$omp do
   do ia = 1, system%num !                    LOOP OVER ATOMS
      aa => system%at(ia)
      NR1 = aa%nszrad     !                     RADIAL MESH
      AZ1 = aa%az
      WS1 = aa%ws
      call RAPO(NR1,WS1,R)
!                                   SPIN DEPENDENT DENSITY
      do is = 1, opt%ns
        do i = 1, nr1
          RHOS(i,is)=aa%rhocor(i,is)+aa%rhoval(i,is)   ! full DOS = core + valence densities of states
        enddo
      enddo
      do i = 1, nr1
        RHOT(i) = RHOS(i,1) + RHOS(i,opt%ns)!    TOTAL CHARGE DENSITY
      enddo

      call HAPO(NR1,R,RHOT,VHART,system%nrmax)!   HARTREE POTENTIAL

      WDUM(1)=rcz
      do i = 2, nr1 !             LOOP OVER RADIAL POINTS
       RHOUP=RHOS(i,1)
       RHODO=RHOS(i,opt%ns)
       if(opt%IVXC==0) CALL XCVBH(RHOUP,RHODO,VXC(1),VXC(2),EXC)  ! Different exch-cor.pots: von Bart - Hadin
       if(opt%IVXC==1) CALL XCCAPZ(RHOUP,RHODO,VXC(1),VXC(2),EXC) !                            Ceperley-Alder
       if(opt%IVXC==2) CALL XCVWN(RHOUP,RHODO,VXC(1),VXC(2),EXC)  !                        Vosko-Wilk-Nusair 

       VNUC = -rc2 * AZ1 / R(i)
       VCOUL = VNUC + VHART(i) + aa%VMAD  ! nuclear pot. + Hartree pot. + Maduleng term = Coulomb interactions
       do is = 1, opt%ns
         aa%fpot(i,is) = VCOUL + VXC(is)  ! Coulomg part + exchange-correlation = part of a new potential
       enddo

       DUM = - RHOUP * aa%pot(i,1) - RHODO*aa%POT(i,opt%ns) + RHOT(i) * (VNUC + rch * VHART(i) + EXC)
       WDUM(i) = DUM * R(i)**2
     enddo

!           Calculate different contribution into an atomic energy
!                 CORE TERM: 
     SCOR = rcz
     if (aa%numcor/=0) then
      if (opt%irel==2) then 
       IADD = 1
      else
       IADD = 0 
      endif
      do J=1, aa%numcor
        L2P1 = 2 * aa%lcor(J) + 1
        AMUL1 = real(L2P1+IADD,kind=prec)
        AMUL2 = real(L2P1-IADD,kind=prec)
        SCOR = SCOR+AMUL1*aa%ECOR(J,1)+AMUL2*aa%ECOR(J,2)
      enddo
     endif

!              VALENCE TERM
     SVAL = rcz
     do il = 1, opt%nl
      SVAL = SVAL + aa%EMDI0(il,1) * aa%ENY(il,1) + aa%EMDI1(il,1) + aa%EMDI0(il,opt%ns) * aa%ENY(il,opt%ns) + aa%EMDI1(il,opt%ns)
     enddo
     COTE = pi4 * QUAD3(NR1,R,WDUM)!       CORRECTION TERM
     ELST = rch * aa%chatra * aa%VMAD + rch * aa%DIPMOM * aa%DMAD!           ELECTROSTATIC TERM

     aa%acte = SCOR + SVAL + COTE + ELST   ! total energy per atom
   enddo
!$omp end do

!                                    DIPOLE BARRIER
   FDBA = rcz
!$omp do reduction(+:FDBA) 
   do ig = 1, system%nums
     FDBA = FDBA + madcon%mdb(ig) * CHTRG(ig) + madcon%ddb(ig) * DIPMG(ig)
   enddo 
!$omp end do
!---------- DEVIATION OF THE POTENTIALS AND THE DIPOLE BARRIER
!$omp sections
   FDBA = FDBA - system%DBA
   
   WRITE(IW6,125) FDBA, iiter
125   FORMAT(/' ***  BARRIER DEVIATION=',G13.3,10X,'ITER=',I4)

!$omp end sections

!$omp do
   do ia = 1, system%num
     aa => system%at(ia)
     NR1 = aa%nszrad
     WS1 = aa%ws
     call RAPO(NR1,WS1,R)
     do is = 1, opt%ns
      aa%fpot(1,is) = rcz
      do i = 2, nr1
        aa%FPOT(i,is) = aa%FPOT(i,is) - aa%POT(i,is)
      enddo
     enddo

     do is = 1, opt%ns
       YMAX(is,ia) = rcz
       IMAX(is,ia) = 0
       do i = 1, nr1
         RDVR = R(i) * dABS(aa%FPOT(i,is))
         if(RDVR.LT.YMAX(is,ia)) cycle
         YMAX(is,ia) = RDVR
         IMAX(is,ia) = i
       enddo
     enddo
   enddo
!$omp end do
!                                           PRINT
!$omp sections
      WRITE(IW6,130) iiTER
130   FORMAT(/' **** ATOM, ',' MAX. DEVIATION(S), MESH POINT(S):',10X,'ITER=',I4)
   opt%wtol=0.0d0
   do ia = 1, system%num
      WRITE(IW6,131) ia,(YMAX(is,ia),IMAX(is,ia),is=1,opt%ns)
      opt%wtol=max(opt%wtol,max(YMAX(1,ia),YMAX(opt%ns,ia)))
   enddo
131   FORMAT(5X,I5,G17.3,I6,G17.3,I6)
!$omp end sections

!---------------   MIXING OF POTENTIALS AND DIPOLE BARRIER
!$omp sections
   J=0
   do ia = 1, system%num
     aa=> system%at(ia) 
     NR1 = aa%nszrad
     do is = 1, opt%ns
      mixing%XL(J+1:J+NR1-1) = aa%pot(2:NR1,is)        
      mixing%FL(J+1:J+NR1-1) = aa%fpot(2:NR1,is)      
  	  J=J+NR1-1
     enddo     
   enddo

   if(opt%ivac==0) then 
     mixing%XL(NAM) = system%DBA
     mixing%FL(NAM) = FDBA
   endif 

   if(iiter.GT.1)  then
      do J = 1, NAM
        mixing%DXP(J,1) = mixing%DXP(J,1) - mixing%XL(J)
        mixing%DFP(J,1) = mixing%DFP(J,1) - mixing%FL(J)
      enddo
      CALL OM1C(NAM,opt%NUH,mixing)

   endif

   if (iiTER.LE.opt%NITERA) then
     do J = 1, NAM
       mixing%XN(J) = mixing%XL(J) + opt%alfa * mixing%FL(J)
     enddo 
   else  
!CCC      NUH1=MIN(NUH,iiTER-NITERA+1)
     NUH1=MIN(opt%NUH,iiTER-1)
     CALL AMST(NAM,NUH1, opt%w0am, opt%beta, mixing)
   endif

   do JH = opt%NUH, 2 , -1
     do IH = JH , opt%NUH
       mixing%VOMA(IH,JH) = mixing%VOMA(IH-1,JH-1)
     enddo
     do J = 1, NAM
       mixing%DXP(J,JH) = mixing%DXP(J,JH-1)
       mixing%DFP(J,JH) = mixing%DFP(J,JH-1)
     enddo
   enddo

   mixing%DXP(:,1) = mixing%XL(:)
   mixing%DFP(:,1) = mixing%FL(:)

   J=0
   do ia = 1, system%num
     aa=>system%at(ia)
     NR1 = aa%nszrad
     do is = 1, opt%ns
       do i = 2, nr1
         J = J + 1
         aa%pot(i,is) = mixing%xn(j)
       enddo
     enddo
   enddo

   if(opt%ivac==0) system%DBA = mixing%xn(NAM)

!----------------- POTENTIAL PARAMETERS OF THE NEW POTENTIALS
!                  VACUUM REGION
   if(opt%ivac==0) then
     S2 = system%vwst**2
     do il = 1, opt%NL
      L   = il-1
      TL1 = real(2*L+1,kind=prec)
      TL5 = real(2*L+5,kind=prec)
      system%vat(1)%eny(il,1) = system%DBA
      system%vat(1)%ppc(il,1) = system%DBA + TL1 * TL5 / ( rc2 * S2 )
     enddo
   endif
!$omp end sections

!$omp do
   do ia = 1, system%num           !     LOOP OVER ATOMS
     aa=>system%at(ia)
     az1 = aa%az
     nr1 = aa%nszrad
     irel1 = min(opt%irel,1)
     if (az1.LT.azl) irel1 = 0
     ws1 = aa%ws
     wsav1 = aa%wsav
     do is = 1, opt%ns             !     LOOP OVER SPINS
      v=>aa%pot(:,is) 
      call rapo (nr1,ws1,r)        !  Preparing the radial mesh
      call lipo (az1,opt%nsirk,nr1,r,v,vi)   !  Lagrange interpolation for potentials
      fip=>fi(:,:,is)
      fidp=>fid(:,:,is)
      fiddp=>fidd(:,:,is)
      do il = 1, opt%nl            !     LOOP OVER L
        ISH1 = aa%ISHENY(il,is)    !                NEW ENY
        if(ISH1==0) ENY1 = aa%ENY(il,is)
        if(ISH1==1.AND.aa%EMDI0(il,is).GT.rcz)  ENY1 = aa%ENY(il,is) + aa%EMDI1(il,is)/aa%EMDI0(il,is)
        if(ISH1==2) ENY1 = aa%PPC(il,is)
        de = dmul / dsqrt(aa%ppp(il,is))
!        write(*,*) az1,eny1,de,wsav1,il,is,irel1,nr1,r,v,vi,fip,fidp,fiddp,opt%nsirk
        call rsel_t (az1,eny1,de,wsav1,il,is,irel1,nr1,r,v,vi,fip,fidp,fiddp,aa,opt%nsirk)
        aa%ENY(il,is)   = eny1
      enddo
     enddo
   enddo 
!$omp end do
!-------------------------------     OUTPUT - UNIT IW7
!$omp end parallel
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)

190   FORMAT(1X,'  --------  LSDA-FILE  -------- ')
191   FORMAT(1X,3G15.7,4X,A16,' IS=',I1)
192   FORMAT(1X,I5,   44X,A16,' IS=',I1)
193   FORMAT(1X,I5,   44X,A16)
195   FORMAT(1X,'  ----------------------------- ')

   return
 end subroutine nepo

!*******************
!XXX    OM1C    ****
!*******************
  subroutine om1c(NAM,NUH,mixing)
!***********************************************************************
! THE FIRST COLUMN OF OVERLAP MATRIX FOR  THE ANDERSON MIXING PROCEDURE
!***********************************************************************
!     NAM - NUMBER OF VARIABLES
!     NUH - NUMBER OF PREVIOUS VECTORS
!-------------------------------------------------
!   ON INPUT:
!     SP(i), i=1,...,NAM,  -  WEIGHTS FOR THE SCALAR PRODUCT
!     DFP(i,IH), i=1,...,NAM, IH=1,...,NUH,  -
!              - DIFFERENCES OF PREVIOUS VECTORS F
!   ON OUTPUT:
!     VOMA(IH,1),  IH=1,...,NUH,  -1ST COLUMN OF OVERLAP MATRIX
!**********************************************************************
    implicit none
    integer:: i,ih,nam,nuh
    real(kind=prec):: VECT(NAM)
    type(andmix):: mixing

    do IH=1,nuh
      mixing%VOMA(IH,1)= 0.0d0
    enddo

    do i=1,NAM
      VECT(i) = mixing%SP(i) * mixing%DFP(i,1)
    enddo
    do IH=1,nuh
      do i=1,NAM
        mixing%VOMA(IH,1) = mixing%VOMA(IH,1) + mixing%DFP(i,IH) * VECT(i)
      enddo
    enddo

   return
  end subroutine om1c

!*******************
!XXX    AMST    ****
!*******************
 subroutine AMST (NAM, NUH, w0am, beta, mixing)
!*********************************************************************************
!  ONE STEP OF THE ANDERSON MIXING PROCEDURE TO SOLVE NON-LINEAR EQUATIONS  F(X)=0
!*********************************************************************************
!    opt%beta - MIXING PARAMETER
!    opt%w0am - A SMALL QUANTITY TO IMPROVE STABILITY
!     NAM - NUMBER OF VARIABLES
!     NUH - NUMBER OF PREVIOUS VECTORS (NUH.GE.2)
!----------------------------------------------------------------------------------
! ON INPUT:
!     DXP(i,IH), i=1,...,NAM, IH=1,...,NUH, - DIFFERENCES OF PREVIOUS VECTORS X
!     DFP(i,IH), i=1,...,NAM, IH=1,...,NUH, - DIFFERENCES OF PREVIOUS VECTORS F
!     SP(i), i=1,...,NAM,  -  WEIGHTS FOR THE SCALAR PRODUCT
!     XL(i), i=1,...,NAM,  -  THE LAST VECTOR X
!     FL(i), i=1,...,NAM,  -  THE LAST VECTOR F
!     VOMA(IH,JH),  IH=1,...,NUH, JH=1,...,IH, - LOWER TRIANGLE OF OVERLAP MATRIX
! ON OUTPUT:
!     XN(i), i=1,...,NAM,   -  THE NEW VECTOR X
!**********************************************************************************
     implicit none
     integer:: ih,jh,i,nuh,nam
     real(kind=prec),parameter:: rcz=0.0d0, rc1=1.0d0
     real(kind=prec):: VECT(NAM),A(NUH,NUH),T(NUH,1),WORK(NUH),dum,bt,w0am,beta
     type(andmix):: mixing

      do IH=1, NUH
       T(IH,1)=rcz
      enddo 
      do i=1, NAM
       mixing%XN(i)=rcz
      enddo
!                                    OVERLAP MATRIX
      do JH = 1, NUH
        do IH = JH, NUH
          A(IH,JH) = mixing%VOMA(IH,JH)
        enddo
      enddo
      DUM=rc1+w0am**2
      do IH = 1, NUH
        A(IH,IH) = DUM * A(IH,IH)
      enddo
!                                     R.H.S.
      do i = 1, NAM
        VECT(i) = mixing%SP(i) * mixing%FL(i)
      enddo

      do IH = 1, NUH
        do i = 1, NAM
         T(IH,1) = T(IH,1) + mixing%DFP(i,IH) * VECT(i)
        enddo
      enddo

!                            SOLUTION OF LINEAR SYSTEM
      CALL SOPO(A,NUH,NUH,T,NUH,1,WORK)
!                                    NEW VECTOR X
      do IH = 1, NUH
        BT = beta * T(IH,1)
        do i = 1, NAM
          mixing%XN(i) = mixing%XN(i) + mixing%DFP(i,IH)*BT
        enddo
      enddo

      do IH = 1, NUH
         do i = 1, NAM
           mixing%XN(i) = mixing%XN(i) + mixing%DXP(i,IH) * T(IH,1)
         enddo
      enddo

      do i = 1, NAM
        mixing%XN(i) = - mixing%XN(i) + beta * mixing%FL(i)
        mixing%XN(i) = mixing%XN(i) + mixing%XL(i)
      enddo

   return
 end subroutine amst


!*******************
!XXX    HAPO    ****
!*******************
 subroutine hapo (n, r, rho, v, nrmax)
!-----------------------------------------------
!   HARTREE POTENTIAL INSIDE ONE ATOMIC SPHERE
!------------------------------------------------
   implicit none
   integer:: i,n,nrmax
   real(kind=prec),parameter:: rc1=1.0d0,rc4=4.0d0,RC8=8.0d0
   real(kind=prec):: R(N),RHO(N),V(N),W1(nrmax),W2(nrmax),pi,pi8,d

   pi=rc4*ATAN(rc1)
   PI8=RC8*pi

   do i = 1, N
     W1(i)=R(i)*RHO(i)*PI8
   enddo
   call PRIM3(N,R,W1,W2)

   d = W2(N)
   do i = 1, N
     V(i)=D-W2(i)
   enddo

   do i = 1, N
     W1(i)=W1(i)*R(i)
   enddo
   call PRIM3(N,R,W1,W2)

   do i = 2, N
     V(i)=V(i)+W2(i)/R(i)
   enddo
  return
 end subroutine hapo


!*******************
!XXX    XCVBH   ****
!*******************
 subroutine XCVBH(RHO1,RHO2,VXC1,VXC2,EXC)
!---------------------------------------------------------------
! XC-POTENTIAL AND ENERGY ACCORDING TO VON BARTH AND HEDIN (J. PHYS. C 5 (1972) 1629)
! WITH PARAMETERS ACCORDING TO HEDIN AND LUNDQVIST (J. PHYS. C 4 (1971) 2064)
! AND JANAK (SOLID STATE COMMUN. 25 (1978) 53)
!---------------------------------------------------------------
  implicit none
  real(kind=prec),parameter:: rcz=0.0d0,rc1=1.0d0,rc2=2.0d0,rc3=3.0d0,rc4=4.0d0, &
  & C238=0.2387324146d0, C916=0.9163305866d0, CP=0.045d0, RP=21.0d0, CF=0.0225d0, RF=52.916684d0
  real(kind=prec)::rho,rho1,rho2,vxc1,vxc2,exc,r,s,a1ps,a1ms,dum,fs,dsfs,epx,efx,ex,rdrex,dsex, & 
  & x,g,epc,rdrepc,xdxg,ec,rdrec,dsec,rdrefc,T13,T43,T213,efc,rdrexc,dsexc

  if(RHO1.LE.rcz) STOP ' ***  RHO1.LE.0  *** '
  if(RHO2.LE.rcz) STOP ' ***  RHO2.LE.0  *** '

!                 CONSTANTS
  T13 = rc1 / rc3
  T43 = rc4 / rc3
  T213 = rc2**T13

!                 VARIABLE RS
   RHO = RHO1 + RHO2
   R = ( C238 / RHO )**T13

!                 VARIABLE DZETA
   S = ( RHO1 - RHO2 ) / RHO
   A1PS = rc1 + S
   A1MS = rc1 - S

!                 FUNCTION F(DZETA)
   DUM  = rc2 * ( T213 - rc1 )
   FS   = ( A1PS**T43 + A1MS**T43 - rc2 ) / DUM
   DSFS = T43 * (A1PS**T13 - A1MS**T13 ) / DUM

!                 EXCHANGE PART
   EPX = - C916 / R
   EFX = T213 * EPX
   EX  = EPX + ( EFX - EPX ) * FS
   RDREX = - EX
   DSEX  = ( EFX - EPX ) * DSFS

!                 CORRELATION PART
   X = R / RP
   CALL AUXVBH(X,G,XDXG)
   EPC    = - CP * G
   RDREPC = - CP * XDXG

   X = R / RF
   CALL AUXVBH(X,G,XDXG)
   EFC    = - CF * G
   RDREFC = - CF * XDXG
   EC     = EPC + ( EFC - EPC ) * FS
   RDREC  = RDREPC + ( RDREFC - RDREPC ) * FS
   DSEC   = ( EFC - EPC ) * DSFS

!                 XC QUANTITIES
   EXC    = EX + EC
   RDREXC = RDREX + RDREC
   DSEXC  = DSEX + DSEC
   DUM    = EXC - T13 * RDREXC
   VXC1   = DUM + A1MS * DSEXC
   VXC2   = DUM - A1PS * DSEXC

  return
 end subroutine XCVBH

!*******************
!XXX   AUXVBH   ****
!*******************
 subroutine AUXVBH(X,G,XDXG)
!-------------------------------------------------------------------------
! AUXILIARY FUNCTION FOR VON BARTH - HEDIN XC. FUNCTION G(X) is DEFINED BY 
! EQ. (5.11) OF THEIR ARTICLE. XDXG MEANS X TIMES X-DERIVATIVE OF G(X).
!-------------------------------------------------------------------------
 implicit none
 real(kind=prec),parameter:: rcz=0.0d0, rc1=1.0d0, rch=0.5d0, rc3=3.0D0, XMAX=20.0d0, RELERR=1.0d-7
 real(kind=prec):: x,g,xdxg,gol,bra,t,tk,suma,sumb,clena,clenb
 integer::k

   if (X.LE.XMAX) then
    GOL = dLOG( rc1 + rc1 / X )
    BRA = X * ( X * ( X * GOL - rc1 ) + rch ) - rc1 / rc3
    G=GOL+BRA
    XDXG=rc3*BRA
   else
    T  = - rc1 / X
    TK = rc1
    K  = 0
    SUMA = rcz
    SUMB = rcz
1   K    = K + 1
    TK   = TK * T
    CLENB = TK / real(K+3)
    CLENA = CLENB / real(K)
    SUMA  = SUMA + CLENA
    SUMB  = SUMB + CLENB
    if(dABS(CLENB).GT.RELERR*dABS(SUMB)) GO TO 1
    G    = -rc3 * SUMA
    XDXG = rc3 * SUMB
   endif
  return
 end subroutine auxvbh

!*******************
!XXX   XCCAPZ   ****
!*******************
  subroutine XCCAPZ(RHO1,RHO2,VXC1,VXC2,EXC)
!-----------------------------------------------------------------------------------
! XC-POTENTIAL AND ENERGY BASED ON WORK OF CEPERLEY AND ALDER (PHYS. REV. LETT. 45 
! (1980) 566)  AS PARAMETRIZED BY PERDEW AND ZUNGER (PHYS. REV. B 23 (1981) 5048)
!-----------------------------------------------------------------------------------
   implicit none
   real(kind=prec),parameter:: rcz=0.0d0, rc1=1.0d0, rc2=2.0d0, rc3=3.0d0, rc4=4.0d0, &
    & C238=0.2387324146d0, C916=0.9163305866d0, GAMP=-0.2846d0, BETP1=1.0529d0, BETP2=0.3334d0, &
    & AP=0.0622d0, BP=-0.096d0, CP=0.0040d0, DP=-0.0232d0, GAMF=-0.1686d0, BETF1=1.3981d0, &
    & BETF2=0.2611d0, AF=0.0311D0, BF=-0.0538d0, CF=0.0014d0, DF=-0.0096d0
   real(kind=prec):: rho,rho1,rho2,r,s,a1ps,a1ms,dum,fs,dsfs,epx,efx,ex,rdrex,dsex,ec,rdrec,dsec,&
    &  exc,rdrexc,dsexc,vxc1,vxc2,t13,t43,t213,epc,efc,rdrepc,rdrefc

   if(RHO1.LE.rcz) STOP ' ***  RHO1.LE.0  *** '
   if(RHO2.LE.rcz) STOP ' ***  RHO2.LE.0  *** '

!                             CONSTANTS
   T13  = rc1 / rc3
   T43  = rc4 / rc3
   T213 = rc2**T13
!                             VARIABLE RS
   RHO = RHO1 + RHO2
   R = ( C238 / RHO )**T13
!                             VARIABLE DZETA
   S = ( RHO1 - RHO2 ) / RHO
   A1PS = rc1 + S
   A1MS = rc1 - S
!                             FUNCTION F(DZETA)
   DUM = rc2 * ( T213 - rc1 )
   FS  = ( A1PS**T43 + A1MS**T43 - rc2 ) / DUM
   DSFS = T43 * ( A1PS**T13 - A1MS**T13 ) / DUM
!                              EXCHANGE PART
   EPX = -C916 / R
   EFX = T213 * EPX
   EX  = EPX + ( EFX - EPX ) * FS
   RDREX = -EX
   DSEX  = ( EFX - EPX ) * DSFS

!                                  CORRELATION PART
   CALL AUXPZ(R,GAMP,BETP1,BETP2,AP,BP,CP,DP,EPC,RDREPC)
   CALL AUXPZ(R,GAMF,BETF1,BETF2,AF,BF,CF,DF,EFC,RDREFC)

   EC=EPC+(EFC-EPC)*FS
   RDREC=RDREPC+(RDREFC-RDREPC)*FS
   DSEC=(EFC-EPC)*DSFS

!                                  XC QUANTITIES
   EXC=EX+EC
   RDREXC=RDREX+RDREC
   DSEXC=DSEX+DSEC
   DUM=EXC-T13*RDREXC
   VXC1=DUM+A1MS*DSEXC
   VXC2=DUM-A1PS*DSEXC

  return 
 end subroutine xccapz

!*******************
!XXX   AUXPZ    ****
!*******************
 subroutine auxpz (RS,GAM,BET1,BET2,A,B,C,D,EC,RDREC)
!-----------------------------------------------------------------------
!  AUXILIARY FUNCTION FOR PERDEW-ZUNGER CORRELATION ENERGY EC DEPENDENT 
!  ON DENSITY VARIABLE RS. RDREC MEANS RS TIMES RS-DERIVATIVE OF EC.
!-----------------------------------------------------------------------
  implicit none
  real(kind=prec),parameter:: rc1=1.0d0,rch=0.5d0
  real(kind=prec):: clen1,clen2,suma,ec,rdrec,gol,a,b,c,d,bet1,bet2,rs,gam

  if(RS.GE.rc1) then
   CLEN1 = BET1 * dSQRT(RS)
   CLEN2 = BET2 * RS
   SUMA  = rc1 + CLEN1 + CLEN2
   EC    = GAM / SUMA
   RDREC = -EC * ( rch * CLEN1 + CLEN2 ) / SUMA
  else
   GOL = LOG(RS)
   EC    = ( A + C * RS ) * GOL + B + D * RS
   RDREC = A + C * RS * GOL + ( C + D ) * RS
  endif
  return
 end subroutine auxpz 

!*******************
!XXX   XCVWN    ****
!*******************
 subroutine xcvwn (RHO1,RHO2,VXC1,VXC2,EXC)
!----------------------------------------------------------------
!      XC-POTENTIAL AND ENERGY ACCORDING TO VOSKO, WILK, NUSAIR
!         (CAN. J. PHYS. 58 (1980) 1200)
!----------------------------------------------------------------
  implicit none
  real(kind=prec),parameter::rcz=0.0d0, rc1=1.0d0, rc2=2.0d0, rc3=3.0d0, rc4=4.0d0, C238=0.2387324146d0, &
  &   C916=0.9163305866d0, AP=0.0621814d0, XP0=-0.10498d0, BP=3.72744d0, CP=12.9352d0, AF=0.0310907d0,   & 
  &   XF0=-0.32500D0, BF=7.06042D0, CF=18.0578D0, AA=-0.0197516d0, XA0=-0.0047584d0, BA=1.13107d0, CA=13.0045D0
  real(kind=prec):: s,rho1,rho2,vxc1,vxc2,exc,rho,r,a1ps,a1ms,s3,s4,dum,fs,dsfs,epx,efx,ex,rdrex,dsex,bra, &
              &  rdrbra,dsbra,ec,rdrec,dsec,rdrexc,dsexc,t13,t43,t213,alt,rdralt,rdrepc,rdrefc,efc,epc

  if(RHO1.LE.rcz) STOP ' ***  RHO1.LE.0  *** '
  if(RHO2.LE.rcz) STOP ' ***  RHO2.LE.0  *** '

!                                   CONSTANTS
  T13  = rc1 / rc3
  T43  = rc4 / rc3
  T213 = rc2**T13
!                                  VARIABLE RS
  RHO = RHO1 + RHO2
  R   = ( C238 / RHO )**T13
!                                  VARIABLE DZETA
  S = ( RHO1 - RHO2 ) / RHO
  A1PS = rc1 + S
  A1MS = rc1 - S
  S3 = S**3
  S4 = S * S3
!                                  FUNCTION F(DZETA)
  DUM = rc2 * ( T213 - rc1 )
  FS  = ( A1PS**T43 + A1MS**T43 - rc2 ) / DUM
  DSFS= T43 * ( A1PS**T13 - A1MS**T13 ) / DUM
!                                  EXCHANGE PART
  EPX = -C916 / R
  EFX = T213 * EPX
!
  EX = EPX + ( EFX - EPX ) * FS
  RDREX = - EX
  DSEX  = ( EFX - EPX ) * DSFS
!                                  CORRELATION PART
  CALL AUXVWN(R,AP,XP0,BP,CP,EPC,RDREPC)
  CALL AUXVWN(R,AF,XF0,BF,CF,EFC,RDREFC)
  CALL AUXVWN(R,AA,XA0,BA,CA,ALT,RDRALT)
!
   BRA = ALT * ( rc1 - S4 ) + ( EFC - EPC ) * S4
   RDRBRA = RDRALT * ( rc1 - S4 ) + ( RDREFC - RDREPC ) * S4
   DSBRA  = rc4 * S3 * ( -ALT + EFC - EPC )
!
   EC = EPC + BRA * FS
   RDREC = RDREPC + RDRBRA * FS
   DSEC  = DSBRA * FS + BRA * DSFS
!                                  XC QUANTITIES
   EXC = EX + EC
   RDREXC = RDREX + RDREC
   DSEXC  = DSEX + DSEC
   DUM  = EXC - T13 * RDREXC
   VXC1 = DUM + A1MS * DSEXC
   VXC2 = DUM - A1PS * DSEXC
!
   return
 end subroutine xcvwn

!*******************
!XXX   AUXVWN   ****
!*******************
 subroutine AUXVWN(RS,A,X0,B,C,Y,RDRY)
!-------------------------------------------------
!  AUXILIARY FUNCTION FOR VOSKO-WILK-NUSAIR CORRELATION ENERGIES
!  (EQS. (4.3) AND (4.4) OF THEIR ARTICLE)
!    RS  ............... DENSITY VARIABLE
!    A, X0, B, C  ...... CONSTANTS
!    Y  ................ DEFINED BY EQ. (4.4)
!    RDRY  ............. DEFINED BY EQ. (4.3)
!-------------------------------------------------
  implicit none
  real(kind=prec),parameter:: rc2=2.0d0, rc4=4.0d0
  real(kind=prec):: a,b,c,rs,x0,y,rdry,x,wx,wx0,xmx0,q,arg1,coef2,arg2,coef3,arg3

  X    = dSQRT(RS)
  WX   = X * ( X + B ) + C
  WX0  = X0 * ( X0 + B ) + C
  XMX0 = X - X0
  Q    = dSQRT(rc4*C-B**2)

  ARG1  = RS / WX
  COEF2 = - B * X0 / WX0
  ARG2  = XMX0**2 / WX
  COEF3 = rc2 * B * ( C - X0**2 ) / ( Q * WX0 )
  ARG3  = Q / ( rc2 * X + B )

  Y = A * (dLOG(ARG1) + COEF2 * dLOG(ARG2) + COEF3 * dATAN(ARG3))

  RDRY = A * ( C - B * X0 * X / XMX0 ) / WX

  return
 end subroutine auxvwn

 
end module pot_mod

