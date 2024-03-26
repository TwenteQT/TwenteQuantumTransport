module dens_mod
 use definitions
 use heritage

contains

!*******************
!XXX    CODE    ****
!*******************
  subroutine code(irel,nsirk,thresh,ns,system)
!********************************************************
! CALCULATES THE CORE DENSITIES AND THE CORE EIGENVALUES
!********************************************************
   implicit none
   type (atoms_set):: system
   real(kind=prec),parameter :: C274=274.074d0,AZL=2.5d0,rcz=0.0d0,rc1=1.0d0,rc4=4.0d0
   real(kind=prec)::pi,pi4,az1,ws1,r(system%nrmax),uc,ucsq,vi(9,system%nrmax),e,cmul,blam,ame,frac,dens, &
   &             wf(system%nrmax),wg(system%nrmax),thresh
   real (kind=prec), pointer ::v(:)
   integer:: is,ia,nr1,irel1,numc1,n,l,nobc1,ik,kappa,i,j,irel,ns,nsirk
   type(atom_definition),pointer::aa

   pi = rc4 * dATAN(rc1)
   pi4 = rc4 * pi

!$omp parallel do default(private) firstprivate(pi,pi4)&
!$omp& shared(system,irel,nsirk,thresh,ns,IW6)
   do ia = 1, system%num !                LOOP OVER ATOMS
    aa=>system%at(ia)
    aa%rhocor=rcz 
    numc1 = aa%numcor
    if(numc1==0) cycle
    nr1 = aa%nszrad
    az1 = aa%az
    ws1 = aa%ws
    call RAPO(nr1, ws1, r)!               RADIAL MESH
    if (az1.LT.AZL) then 
      irel1 = 0
     else
      irel1 = irel
    endif 

    if (IREL1/=2) then!----- NON-RELATIVISTIC AND SCALAR-RELATIVISTIC CASE (INCLUDING SPIN POLARIZATION)
     UC = real(IREL1) / C274
     UCSQ = UC**2
     do is = 1, ns  !                 LOOP OVER SPIN
      v => aa%pot(:,is)

      call lipo(az1,nsirk,nr1,r,v,vi)

      do j = 1, NUMC1 !               LOOP OVER CORE ORBITALS
       E = aa%ECOR(j,is)
       N = aa%NCOR(j)
       L = aa%LCOR(j)
       NOBC1 = aa%NOBC(j)

       call RSEC(AZ1,E,THRESH,N,L,IREL1,NR1,NOBC1,r,v,vi,wf,wg,nsirk)!   SCHROEDINGER EQUATION
       CMUL = real ( 2 * L + 1 ) / pi4                         !   ADDITION TO CORE DENSITY
       BLAM = real ( L * ( L + 1 ) )

       aa%rhocor(1,is) =  aa%rhocor(1,is) + CMUL * WG(1)**2
       do i = 2, NR1
         AME = rc1 + UCSQ * ( E - V(i) )
         FRAC = WG(i) / ( R(i) * AME )
         DENS = WG(i)**2 + UCSQ * (WF(i)**2 + BLAM*FRAC**2)
         aa%rhocor(i,is) = aa%rhocor(i,is) + CMUL*DENS
       enddo
       aa%ecor(j,is) = e !                                   EIGENVALUE
      enddo
     enddo

      if(ns == 1) then  ! even in case of 1 spin we have ECOR for both spins but in this case they are identical
       do j = 1, NUMC1
         aa%ecor(j,2) = aa%ecor(j,1)
       enddo
      endif
    endif

    if(IREL1==2) then!---------------FULLY RELATIVISTIC CASE (WITHOUT SPIN POLARIZATION)
      UC = rc1/C274
      UCSQ = UC**2
      v = aa%pot(:,is)

      call lipo(az1,nsirk,nr1,r,v,vi)

      do j=1,NUMC1 !         LOOP OVER CORE ORBITALS
       N=aa%NCOR(j)
       L=aa%LCOR(j)
       NOBC1=aa%NOBC(j)
       do IK = 1, 2 !               LOOP OVER KAPPA
        if(IK==1) KAPPA=-L-1
        if(IK==2) KAPPA=L
        if(KAPPA.NE.0) then
         E=aa%ecor(j,IK)

         call RDEC(AZ1,E,THRESH,N,KAPPA,NR1,NOBC1,r,v,vi,wf,wg,nsirk) !              DIRAC EQUATION

         CMUL=real(abs(KAPPA))/pi4 !  ADDITION TO CORE DENSITY
         do i = 1, nr1
          DENS = WG(i)**2 + UCSQ * WF(i)**2
          aa%rhocor(i,1) = aa%rhocor(i,1) + CMUL * DENS
         enddo
        endif 
        aa%ECOR(j,IK)=E  !             EIGENVALUE
       enddo
      enddo
    endif
   enddo
!$omp end parallel do

  return
 end subroutine code

!*******************
!XXX    RSEC    ****
!*******************
  subroutine rsec (AZ,E,THRESH,N,L,IREL,NR,NOBC,r,v,vi,wf,wg,nsirk)
!-------------------------------------------------------------------
! SOLUTION OF RADIAL SCHROEDINGER EQUATION FOR CORE ELECTRONS: BOTH
! NON-RELATIVISTIC (IREL=0) AND SCALAR-RELATIVISTIC (IREL=1) VERSION
!-------------------------------------------------------------------
!  INPUT:
!     AZ - ATOMIC NUMBER
!     E - ENERGY (ESTIMATE)
!     THRESH - ABSOLUTE ACCURACY OF THE EIGENVALUE
!     N - PRINCIPAL QUANTUM NUMBER
!     L - ORBITAL QUANTUM NUMBER
!     IREL - RELATIVITY
!     NR - SIZE OF RADIAL MESH
!     NOBC - OUTER BOUNDARY CONDITION
!           (NOBC=0 - DEEP LEVEL, NOBC=1 - SHALLOW LEVEL)
!     R(.) - RADIAL MESH
!     V(.) - POTENTIAL
!     VI(.,.) - INTERPOLATED POTENTIAL
!     NSIRK - NO. OF SUBINTERVALS FOR RUNGE-KUTTA METHOD
!  OUTPUT:
!     E - ENERGY
!     WG(.) - WAVE FUNCTION (LARGE COMPONENT)
!     WF(.) - FUNCTION RELATED TO RADIAL DERIVATIVE OF WG
!             (SMALL COMPONENT)
!---------------------------------------------------------
      implicit none
      real(kind=prec):: AME(NR),U(NR),P(NR),Q(NR),WRK(NR),TAPP(0:10),TAPQ(0:10),TAQP(0:10),TAQQ(0:10),TR(0:10),TV(0:10)
      real(kind=prec),parameter :: ARGMAX=120.0d0,DSG=1.5d0,C274=274.074d0,rcz=0.0d0,rc1=1.0d0,rc2=2.0d0,rc3=3.0d0,rch=0.5d0
      integer,parameter:: NRST=10,MAXIT=200
      real(kind=prec)::al,alp1,blam,twoz,uc,ucsq,un2,u0,a1,a2,b1,b2,beta,a1pb,a1mb,th,tp,tq,rtb,asig,bsig,sum,cnorm,emin,&
                  & emax,sigde,frac,anil,xld1,dld1,az,e,thresh,v(:),vi(:,:),r(:),abse,amgde,pmal,xldl,vme,pmar,xldr,  &
                  & anir,xsi,exfa,de,dum,wf(:),wg(:)
      integer :: i,j,numnod,nsi2,it,imin,imax,nrma,nuno,n,l,irel,nr,nobc,nrii,nsirk

!                             SET UP CONSTANTS
      AL   = real(L)
      ALP1 = real(L + 1)
      BLAM = AL * ALP1
      TWOZ = rc2 * AZ
      NUMNOD = N - L - 1
      UC = real(IREL) / C274
      UCSQ = UC**2
      NSI2 = 2*NSIRK
      UN2 = rc1 / real(NSI2)
      IT = 0
      IMIN = 0
      IMAX = 0

200   IT = IT + 1
      if(IT.GT.MAXIT) STOP ' *** ig>MAXIT, RSEC - NO CONVERGENCY '
      P = rcz
      Q = rcz
      do i=2, NR
        AME(i) = rc1 + UCSQ * ( E - V(i) )
        U(i)   = BLAM / ( R(i)**2 * AME(i) ) + V(i) - E
      enddo

!                            CLASSICAL TURNING POINT
      i = NR + 1
205   i = i - 1
      if (U(i).GT.rcz.AND.i.GT.NRST+10) GO TO 205
      NRMA = MIN(i+1,NR)

!---------------   OUTWARD INTEGRATION

      if (IREL==0) then !--------------START -  NON-RELATIVISTIC
       U0 = V(2) + TWOZ/R(2) - E
       A1 = -AZ / ALP1
       A2 = ( U0 - TWOZ * A1 ) / ( rc2 * ( rc2 * AL + rc3 ) )
       B1 = -AZ
       B2 = A2 * ( AL + rc2 )
       P(2) = rc1 + R(2) * ( A1 + R(2) * A2 )
       Q(2) = B1 + R(2) * B2
       BETA = rc1
       if(L.GT.0) then
         P(2)=R(2)*P(2)
         Q(2)=AL+R(2)*Q(2)
         BETA=AL
       endif
      endif

      if(IREL==1) then !----------------START -  SCALAR-RELATIVISTIC
       BETA = dSQRT( BLAM + rc1 - ( TWOZ * UC )**2 )
       P(2) = rc1
       Q(2) = ( BETA - rc1 ) / ( TWOZ * UCSQ )
      endif

!                       STARTING RUNGE-KUTTA INTEGRATION
      if(NRST==2) GO TO 208
       A1PB = rc1 + BETA
       A1MB = rc1 - BETA
       do i = 3, NRST
        do j = 0, NSI2
          TR(j) = ( R(i-1) * real(NSI2 - j) + R(i) * real(j) ) * UN2
        enddo
        TV(0) = V(i-1)
        do j = 1, NSI2 - 1
         TV(j) = VI(j,i)
        enddo
        TV(NSI2) = V(i)

        do j = 0, NSI2
         TAPP(j) = A1MB / TR(j)
         TAPQ(j) = rc1 + UCSQ * (E - TV(j) )
         TAQP(j) = BLAM / (TAPQ(j) * TR(j)**2 ) + TV(j) - E
         TAQQ(j) = -A1PB / TR(j)
        enddo
        TH = R(i)-R(i-1)
        TP = P(i-1)
        TQ = Q(i-1)

        call rukust(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)

        P(i) = TP
        Q(i) = TQ
       enddo

208    do i = 2, NRST
        RTB = R(i)**BETA
        P(i) = RTB * P(i)
        Q(i) = RTB * Q(i)
       enddo

!                             RUNGE-KUTTA INTEGRATION
      do i = NRST + 1, NRMA
       do j=0,NSI2
        TR(j) = (R(i-1) * real( NSI2 - j ) + R(i) * real( j ) ) * UN2
       enddo
       TV(0) = V(i-1)
       do j = 1,NSI2 - 1
        TV(j) = VI(j,i)
       enddo
       TV(NSI2) = V(i)

       do j = 0, NSI2
        TAPP(j) = rc1 / TR(j)
        TAPQ(j) = rc1 + UCSQ * ( E - TV( j ) )
        TAQP(j) = BLAM / ( TAPQ( j ) * TR( j )**2 ) + TV( j ) -E
        TAQQ(j) = -rc1 / TR( j )
       enddo

       TH = R( i ) - R( i - 1 )
       TP = P( i - 1 )
       TQ = Q( i - 1 )
       call RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
       P( i ) = TP
       Q( i ) = TQ
      enddo
!                                    NUMBER OF NODES
      NUNO=0
      ASIG=SIGN(rc1,P(2))
      do i=3,NRMA
       BSIG=SIGN(rc1,P(i))
       if(abs(ASIG-BSIG).LT.DSG) cycle
       NUNO=NUNO+1
       ASIG=BSIG
      enddo

      if(NUNO.LT.NUMNOD) then
       if(IMIN==0) EMIN=E !                                   ENERGY TOO LOW
       if(IMIN==1) EMIN=MAX(EMIN,E)
       IMIN=1
       SIGDE=rc1
      endif

      if(NUNO.GT.NUMNOD) then
       if(IMAX==0) EMAX=E !                                   ENERGY TOO HIGH
       if(IMAX==1) EMAX=MIN(EMAX,E)
       IMAX=1
       SIGDE=-rc1
      endif

      if(NUNO==NUMNOD) GO TO 250

      if(IMIN*IMAX==0) then
       ABSE=abs(E)
       AMGDE=MAX(ABSE,rc1)/real(10)
       E=E+SIGDE*AMGDE
      else
       E=rch*(EMAX+EMIN)
      endif
      GO TO 200

!                            ENERGY IN THE CORRECT WINDOW
250   PMAL = P( NRMA )
      XLDL = Q( NRMA ) / PMAL
      WRK(1) = rcz
      do i = 2, NRMA
       FRAC = P(i) / ( R( i ) * AME( i ) )
       WRK(i) = P(i)**2 + UCSQ * (Q(i)**2 + BLAM * FRAC**2)
      enddo
      ANIL=rcz
      do i = 2, NRMA
       ANIL = ANIL + ( WRK( i - 1 ) + WRK( i ) ) * (R( i ) - R( i - 1 ) )
      enddo
      ANIL = rch * ANIL

!                    STARTING POINT FOR INWARD INTEGRATION
      if(NOBC==0) then
       XLD1 = -dSQRT( U( NR ) / AME( NR ) )
       DLD1 = -rch / XLD1
      endif
      if(NOBC==1) then
       VME = V( NR ) - E
       CALL COBC(XLD1,DLD1,VME,R(NR),L,-1)
      END if

      if(NRMA.LT.NR) GO TO 253
      PMAR=rc1
      XLDR=XLD1
      ANIR=rcz
      GO TO 260

253   do i=NRMA,NR
        WRK(i)=dSQRT(U(i)*AME(i))
      enddo

      sum=rcz
      i=NRMA
255   i=i+1
      sum=sum+(WRK(i-1)+WRK(i))*(R(i)-R(i-1))
      if(sum.LT.ARGMAX.AND.i.LT.NR) GO TO 255
      NRii=i

!------------------------ INWARD INTEGRATION
!        START
      if(NRii.LT.NR) THEN
       XLD1=-SQRT(U(NRii)/AME(NRii))
       DLD1=-rch/XLD1
      END if
      P(NRii)=rc1
      Q(NRii)=XLD1
!                                RUNGE-KUTTA INTEGRATION
      do i=NRii-1,NRMA,-1
       do j=0,NSI2
         TR(j) = (R(i+1)*real(NSI2-j)+R(i)*real(j))*UN2
       enddo
       TV(0)=V(i+1)
       do j=1,NSI2-1
         TV(j) = VI(NSI2-j,i+1)
       enddo
       TV(NSI2)=V(i)

       XSI=WRK(i)
       do j=0,NSI2
        TAPP(j) = XSI + rc1 / TR( j )
        TAPQ(j) = rc1 + UCSQ * ( E - TV( j ) )
        TAQP(j) = BLAM / (TAPQ( j ) * TR( j )**2 ) + TV( j ) - E
        TAQQ(j) = XSI - rc1 / TR( j )
       enddo
       TH=R(i)-R(i+1)
       TP=P(i+1)
       TQ=Q(i+1)

       call RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)

       EXFA = dEXP( -XSI * TH )
       P(i) = EXFA * TP
       Q(i) = EXFA * TQ
      enddo

      PMAR = P(NRMA)
      XLDR = Q(NRMA) / PMAR
      do i = NRMA, NRii
       FRAC = P( i ) / ( R( i ) * AME( i ) )
       WRK( i ) = P(i)**2 + UCSQ *( Q( i )**2 + BLAM * FRAC**2 )
      enddo
      ANIR = rcz
      do i = NRMA+1, NRii
       ANIR = ANIR + ( WRK( i - 1 ) + WRK( i ) ) * ( R( i ) - R( i - 1 ) )
      enddo
      ANIR = rch * ANIR

!                   ENERGY CHANGE FROM PERTURBATION THEORY
260   DE = ( XLDL - XLDR ) / ( ANIL / PMAL**2 + ( DLD1 + ANIR ) / PMAR**2 )
      if(abs(DE).LT.THRESH) GO TO 280

      if(DE.GT.rcz) THEN
       if(IMIN==0) EMIN=E
       if(IMIN==1) EMIN=MAX(EMIN,E)
       IMIN=1
      ELSE
       if(IMAX==0) EMAX=E
       if(IMAX==1) EMAX=MIN(EMAX,E)
       IMAX=1
      END if
      E=E+DE
      GO TO 200
!                               MATCHING OF BOTH PARTS
280   if(NRMA==NR) GO TO 290
      DUM=PMAL/PMAR
      do i=NRMA,NRii
       P(i)=P(i)*DUM
       Q(i)=Q(i)*DUM
      enddo
!                                  NORMALIZATION
290   WRK(1)=rcz
      do i=2,NR
       FRAC=P(i)/(R(i)*AME(i))
       WRK(i)=P(i)**2 + UCSQ*(Q(i)**2+BLAM*FRAC**2)
      enddo

      sum=QUAD3(NR,R,WRK)
      CNORM=SQRT(sum)
!
      WG(1)=rcz
      WF(1)=rcz
      if(IREL==0) THEN
       if(L==0) THEN
         WG(1)=rc1/CNORM
         WF(1)=-AZ/CNORM
       END if
       if(L==1) WF(1)=rc1/CNORM
      END if
!
      do i=2,NR
       DUM=CNORM*R(i)
       WG(i)=P(i)/DUM
       WF(i)=Q(i)/DUM
      enddo
!
     return
   end subroutine rsec

!*******************
!XXX    RDEC    ****
!*******************
  subroutine rdec (AZ,E,THRESH,N,KAPPA,NR,NOBC,r,v,vi,wf,wg,nsirk)
!---------------------------------------------------------
! SOLUTION OF RADIAL DIRAC EQUATION FOR CORE ELECTRONS (NON-MAGNETIC CASE)
!---------------------------------------------------------
!  INPUT:
!     AZ - ATOMIC NUMBER
!     E - ENERGY (ESTIMATE)
!     THRESH - ABSOLUTE ACCURACY OF THE EIGENVALUE
!     N - PRINCIPAL QUANTUM NUMBER
!     KAPPA - RELATIVISTIC QUANTUM NUMBER
!     NR - SIZE OF RADIAL MESH
!     NOBC - OUTER BOUNDARY CONDITION
!           (NOBC=0 - DEEP LEVEL, NOBC=1 - SHALLOW LEVEL)
!     R(.) - RADIAL MESH
!     V(.) - POTENTIAL
!     VI(.,.) - INTERPOLATED POTENTIAL
!     NSIRK - NO. OF SUBINTERVALS FOR RUNGE-KUTTA METHOD
!  OUTPUT:
!     E - ENERGY
!     WG(.) - LARGE COMPONENT (G)
!     WF(.) - SMALL COMPONENT (C*F) OF WAVE FUNCTION NORMALIZED TO UNITY
!---------------------------------------------------------
    implicit none
    real(kind=prec)::AME(NR),U(NR),P(NR),Q(NR),WRK(NR),TAPP(0:10),TAPQ(0:10),TAQP(0:10),TAQQ(0:10),TR(0:10),TV(0:10), &
    &     ak,twoz,uc,ucsq,un2,beta,akpb,akmb,th,tp,tq,dum,cnorm,sum,asig,bsig,sigde,abse,amgde,e,az,thresh,        &
    &     r(:),v(:),vi(:,:),rtb,emin,emax,anil,xld1,dld1,de,pmal,xldl,vme,pmar,xldr,anir,xsi,exfa,wf(:),wg(:)
    real(kind=prec),parameter::ARGMAX=120.0d0,DSG=1.5D0,C274=274.074D0,rcz=0.0D0,rc1=1.0D0,rc2=2.0D0,rch=0.5D0
    integer,parameter:: NRST=10, MAXIT=200
    integer :: kappa,l,numnod,nsi2,it,imin,imax,i,j,nrma,nobc,nr,n,nuno,nrii,nsirk

!                             SET UP CONSTANTS
      AK = real(KAPPA)
      if(KAPPA.GT.0) L = KAPPA
      if(KAPPA.LT.0) L = -KAPPA - 1
      NUMNOD=N-L-1
      TWOZ=rc2*AZ
      UC=rc1/C274
      UCSQ=UC**2
      NSI2=2*NSIRK
      UN2=rc1/real(NSI2)
      IT=0
      IMIN=0
      IMAX=0

200   IT=IT+1
      if(IT.GT.MAXIT) STOP 'IT<MAXIT *** RDEC - NO CONVERGENCY '

      P = rcz
      Q = rcz

      do i = 2, NR
       AME(i) = rc1 + UCSQ * ( E - V( i ) )
       U(i)   = V( i ) - E
      enddo

!                                   MATCHING RADIUS
      i = NR + 1
205   i = i - 1
      if(U(i).GT.rcz.AND.i.GT.NRST+10) GO TO 205
      NRMA = MIN(i+1,NR)

!----------------- OUTWARD INTEGRATION INITIAL CONDITION

      BETA = dSQRT( AK**2 - ( TWOZ * UC )**2 )
      P(2) = rc1
      Q(2) = ( BETA + AK ) / ( TWOZ * UCSQ )
!                       STARTING RUNGE-KUTTA INTEGRATION
      if(NRST==2) GO TO 208
      AKPB = AK + BETA
      AKMB = AK - BETA

      do i = 3, NRST
       do j = 0, NSI2
         TR(j) = (R(i-1)*real(NSI2-j)+R(i)*real(j))*UN2
       enddo
       TV(0) = V(i-1)
       do j = 1,NSI2-1
         TV(j) = VI(j,i)
       enddo 
       TV(NSI2) = V(i)

       do j = 0, NSI2
         TAPP(j) = -AKPB / TR( j )
         TAPQ(j) = rc1 + UCSQ * ( E - TV( j ) )
         TAQP(j) = TV( j ) - E
         TAQQ(j) = AKMB / TR(j)
       enddo

       TH=R(i)-R(i-1)
       TP=P(i-1)
       TQ=Q(i-1)
       CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
       P(i)=TP
       Q(i)=TQ
      enddo


208   do i = 2, NRST
       RTB    = R( i )**BETA
       P( i ) = RTB * P( i )
       Q( i ) = RTB * Q( i )
      enddo

!                             RUNGE-KUTTA INTEGRATION
      do i = NRST + 1, NRMA
        do j = 0, NSI2
          TR( j ) = ( R( i-1 ) * real( NSI2 - j ) + R( i ) * real( j ) ) * UN2
        enddo
        TV( 0 ) = V( i - 1 )
        do j=1,NSI2-1
          TV( j ) = VI( j , i )
        enddo
        TV( NSI2 ) = V( i )

        do j = 0, NSI2
         TAPP( j ) = -AK / TR( j )
         TAPQ( j ) = rc1 + UCSQ * ( E - TV( j ) )
         TAQP( j ) = TV( j ) - E
         TAQQ( j ) = AK / TR( j )
        enddo

        TH=R(i)-R(i-1)
        TP=P(i-1)
        TQ=Q(i-1)

        CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)

        P(i)=TP
        Q(i)=TQ
      enddo

!                                    NUMBER OF NODES
      NUNO = 0
      ASIG = SIGN(rc1,P(2))
      do i=3,NRMA
        BSIG=SIGN(rc1,P(i))
        if(abs(ASIG-BSIG).LT.DSG) cycle
        NUNO=NUNO+1
        ASIG=BSIG
      enddo

      if(NUNO.NE.NUMNOD) then  ! energy is not in the correct window

       if (NUNO.LT.NUMNOD) then  !               ENERGY TOO LOW
        if(IMIN==0) EMIN=E
        if(IMIN==1) EMIN=MAX(EMIN,E)
        IMIN=1
        SIGDE=rc1
       elseif (NUNO.GT.NUMNOD) then   !          ENERGY TOO HIGH
        if(IMAX==0) EMAX=E
        if(IMAX==1) EMAX=MIN(EMAX,E)
        IMAX=1
        SIGDE=-rc1
       endif

       if(IMIN*IMAX==0) THEN
        ABSE=abs(E)
        AMGDE=MAX(ABSE,rc1)/10.0d0
        E=E+SIGDE*AMGDE
       ELSE
        E=rch*(EMAX+EMIN)
       END if
       goto 200

      endif


      PMAL = P(NRMA)              !         ENERGY IN THE CORRECT WINDOW
      XLDL = Q(NRMA) / PMAL
      do i = 1, NRMA
        WRK( i ) = P( i )**2 + UCSQ * Q( i )**2
      enddo
      ANIL=rcz
      do i = 2, NRMA
        ANIL = ANIL + ( WRK( i - 1 ) + WRK( i ) ) * ( R( i ) - R( i - 1 ) )
      enddo
      ANIL = rch * ANIL

!                    STARTING POINT FOR INWARD INTEGRATION
      if(NOBC==0) THEN
       XLD1=-SQRT(U(NR)/AME(NR))
       DLD1=-rch/XLD1
      END if
      if(NOBC==1) THEN
       VME=V(NR)-E
       CALL COBC(XLD1,DLD1,VME,R(NR),L,KAPPA)
      END if

      if(NRMA.LT.NR) GO TO 253
      PMAR=rc1
      XLDR=XLD1
      ANIR=rcz
      GO TO 260

253   do i=NRMA,NR
        WRK(i)=SQRT(U(i)*AME(i))
      enddo
      sum = rcz
      i = NRMA
255   i = i + 1
      sum = sum + (WRK ( i - 1 ) + WRK( i ) ) * ( R( i ) - R( i - 1 ) )
      if(sum.LT.ARGMAX.AND.i.LT.NR) GO TO 255
      NRii=i

! ------------------------------------ INWARD INTEGRATION
!                                   START
     if(NRii.LT.NR) THEN
       XLD1 = -dSQRT( U( NRii ) / AME( NRii ) )
       DLD1 = -rch / XLD1
     END if
     P(NRii) = rc1
     Q(NRii) = XLD1
!                                RUNGE-KUTTA INTEGRATION
     do i = NRii - 1, NRMA , -1

       do j = 0, NSI2
        TR( j ) = ( R( i + 1 ) * real( NSI2 - j ) + R( i ) * real( j ) ) * UN2
       enddo
       TV( 0 ) = V( i + 1 )
       do j = 1, NSI2 - 1
        TV( j ) = VI( NSI2 - j , i + 1 )
       enddo
       TV(NSI2)=V(i)

       XSI=WRK(i)
       do j=0,NSI2
         TAPP(j) = XSI - AK / TR( j )
         TAPQ(j) = rc1 + UCSQ * ( E - TV( j ) )
         TAQP(j) = TV(j) - E
         TAQQ(j) = XSI + AK / TR(j)
       enddo

       TH=R(i)-R(i+1)
       TP=P(i+1)
       TQ=Q(i+1)

       CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)

       EXFA=EXP(-XSI*TH)
       P(i)=EXFA*TP
       Q(i)=EXFA*TQ
     enddo

     PMAR = P( NRMA )
     XLDR = Q( NRMA ) / PMAR
     do i = NRMA, NRii
      WRK(i)=P(i)**2+UCSQ*Q(i)**2
     enddo
     ANIR=rcz
     do i=NRMA+1,NRii
      ANIR=ANIR+(WRK(i-1)+WRK(i))*(R(i)-R(i-1))
     enddo
     ANIR=rch*ANIR

!                   ENERGY CHANGE FROM PERTURBATION THEORY
260   DE=(XLDL-XLDR)/(ANIL/PMAL**2 + (DLD1+ANIR)/PMAR**2)
      if(abs(DE).LT.THRESH) GO TO 280

      if(DE.GT.rcz) THEN
       if(IMIN==0) EMIN=E
       if(IMIN==1) EMIN=MAX(EMIN,E)
       IMIN=1
      ELSE
       if(IMAX==0) EMAX=E
       if(IMAX==1) EMAX=MIN(EMAX,E)
       IMAX=1
      END if
      E=E+DE
      GO TO 200
!                               MATCHING OF BOTH PARTS
280   if(NRMA==NR) GO TO 290
      DUM=PMAL/PMAR
      do i=NRMA,NRii
        P(i)=P(i)*DUM
        Q(i)=Q(i)*DUM
      enddo
!                                  NORMALIZATION
290   do i=1,NR
        WRK(i)=P(i)**2+UCSQ*Q(i)**2
      enddo

      sum=QUAD3(NR,R,WRK)
      CNORM=SQRT(sum)

      WG(1)=rcz
      WF(1)=rcz
      do i=2,NR
       DUM=CNORM*R(i)
       WG(i)=P(i)/DUM
       WF(i)=Q(i)/DUM
      enddo
   return
  end subroutine rdec

!*******************
!XXX    COBC    ****
!*******************
 subroutine cobc(XLD,DLD,VME,R,L,KAPPA)
!-------------------------------------------------------------------------------
! OUTER BOUNDARY CONDITION FOR CORE ORBITALS BASED ON SPHERICAL HANKEL FUNCTIONS
!-------------------------------------------------------------------------------
  implicit none
  real(kind=prec),parameter :: rcz=0.0d0,rc1=1.0d0,rc2=2.0d0,rc3=3.0d0,RC6=6.0d0,RC15=15.0d0
  real(kind=prec):: vme,u,xsi,r,t,y,yp,w,ak,bl,xld,h,dld
  integer:: l,kappa 

   if(VME.LE.rcz) STOP ' ***  COBC: VME.LE.0  *** '
   if(L.LT.0.OR.L.GT.3) STOP ' ***  COBC: WRONG L  *** '
   XSI = dSQRT(VME)
   T   = XSI * R
   U   = rc1 / T

   if(L==0) then
     Y  = rc1
     YP = rcz
   elseif(L==1) THEN
     Y  = rc1 + U
     YP = -U**2
   elseif(L==2) THEN
     Y  = rc1 + rc3 * U * ( rc1 + U )
     YP = -U**2 * ( rc3 + RC6 * U )
   elseif(L==3) THEN
     Y  = rc1 + U * ( RC6 + RC15 * U * ( rc1 + U ) )
     YP = -U**2 * ( RC6 + RC15 * U * ( rc2 + rc3 * U ) )
   endif
     W  = YP / Y
     AK = real(KAPPA)
     BL = real(L) * real( L + 1 )
     XLD = XSI * ( -rc1 + AK * U + W )
     H = -rc1 + W * ( rc1 + T * ( rc2 - W ) ) + BL * U
     DLD = -H / ( rc2 * XSI )
   return
 end subroutine cobc

!*******************
!XXX    VADE    ****
!*******************
 subroutine vade(ns,nl,system)
!********************************************************************
! CALCULATES THE VALENCE DENSITIES FROM THE L-DIAGONAL MOMENTS OF DOS
!********************************************************************
  implicit none
  type (atoms_set):: system
  real (kind=prec),parameter:: rcz=0.0d0,rc1=1.0d0,rc2=2.0d0,rc4=4.0d0
  real (kind=prec):: pi,pi4,wm0,tm1,wm2,f,fd,fdd
  integer:: ia,is,il,i,nl,ns 
  type(atom_definition),pointer :: aa

   pi = rc4*dATAN(rc1)
   pi4 = rc4*pi

!$omp parallel do default(private) firstprivate(pi,pi4)&
!$omp& shared(system,ns,nl,IW6)
   do ia = 1, system%num
     aa => system%at(ia)
     aa%rhoval = rcz 
     do is = 1, ns
       do il = 1, nl
        WM0 = aa%emdi0(il,is) / pi4
        TM1 = rc2 * aa%emdi1(il,is) / pi4
        WM2 = aa%emdi2(il,is) / pi4

        do i = 1, aa%nszrad
         f   = aa%phi(i,il,is)
         fd  = aa%phid(i,il,is)
         fdd = aa%phidd(i,il,is)
         aa%rhoval(i,is) = aa%rhoval(i,is) + WM0*F**2 + TM1 * F * FD + WM2 * ( FD**2 + F * FDD )
        enddo
       enddo
      enddo
     enddo
!$omp end parallel do
    return
 end subroutine vade

end module dens_mod
