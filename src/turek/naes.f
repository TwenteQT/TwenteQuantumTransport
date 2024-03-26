      PROGRAM AES
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C************************************************************
C
C      CALCULATES POTENTIAL PARAMETERS, POTENTIALS,
C      AND CORE ENERGIES FOR SELF-CONSISTENT CALCULATION
C      FOR A SINGLE ATOMIC SPHERE
C
C       NON-RELATIVISTIC (IREL=0) AND
C       SCALAR-RELATIVISTIC (IREL=1) VERSION
C
C************************************************************
C        VERSION  APRIL 1998
C************************************************************
C
C  INPUT: IR5 (=5) - GENERAL DATA
C OUTPUT: IW6 (=6) - CURRENT OUTPUT
C         IW8 (=8) - RADIAL DEPENDENCE OF CHARGE DENSITY
C
C************************************************************
C
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
      COMMON/IDAT/ AZ,WS,VOCC(4,2),NL,NS,IREL,NR,NVAL(4)
      COMMON/PP/ ENY(4,2),PPC(4,2),PPD(4,2),PPQ(4,2),
     &           PPP(4,2),DNY(4,2),FINY(4,2),FINYD(4,2)
      COMMON/RW/ IR5,IW6,IW8
      COMMON/COR/ ECOR(20,2),THRESH,NCOR(20),LCOR(20),MCOR
      COMMON/PTO/ POT(MNR,2)
      COMMON/ITER/ DEV(2),DEVM,ALFA,IITER,NITER
      COMMON/ROH/ RHOCOR(MNR,2),RHOVAL(MNR,2)
      COMMON/LIN/ ENY1,PPC1,PPD1,PPQ1,PPP1,DNY1,FINY1,FINYD1,
     &            FI(MNR),FID(MNR),FIDD(MNR)
C
C
C                                       READING AND START
      CALL ATOMA
C                                      LOOP
      IITER=0
C
222   IITER=IITER+1
C                                       CORE DENSITIES
      CALL ATOMC
C                                       VALENCE DENSITIES
      CALL ATOMD
C                                NEW POTENTIALS & OUTPUT
      CALL ATOMF
C
      DEVY= MAX(DEV(1),DEV(NS))
      IF(IITER.LT.NITER.AND.DEVY.GT.DEVM) GO TO 222
C
      WRITE(IW6,150)
150   FORMAT(//1X,'   ***   END OF RUN   *** ')
C
      STOP
      END

      subroutine write_atoms(NL,NS,NSZRAD,AZ,WS,POT)
      implicit none           
      integer,PARAMETER::MNR=512

      COMMON/PP/ ENY(4,2),PPC(4,2),PPD(4,2),PPQ(4,2),
     &           PPP(4,2),DNY(4,2),FINY(4,2),FINYD(4,2)
      COMMON/COR/ ECOR(20,2),THRESH,NCOR(20),LCOR(20),MCOR

      integer:: LCOR,NCOR,MCOR,NOBC
      real(kind=8)::ECOR,ENY,PPC,PPD,PPP,FINY,DNY,PPQ,FINYD
      real(kind=8)::THRESH
      
      integer::NL,NS
      integer::NSZRAD
      REAL(kind=8)::AZ,WS
      REAL(kind=8)::POT(MNR,2)
      CHARACTER*16::OTXTA='AES'
      
      integer::IS,IL,I,J
      integer,parameter::atomf=999
      character(len=16):: cwork

992   FORMAT(1X,I5,   54X,'# ',A10,' IS=',I1)
993   FORMAT(1X,I5,   42X,A16)      
994   FORMAT(1X,4G15.7)      
995   FORMAT(1X,G15.7,   42X,A16)
996   FORMAT(1X,E10.1,   44X,A16)
997   FORMAT(1X,3G15.7,7X,A16)

1101   FORMAT(1X,10I5)
1104   FORMAT(1X,4G15.7)
C
1190   FORMAT(1X,'  --------  LSDA-FILE  -------- ')
1191   FORMAT(1X,49X,A16,' IS=',I1)
1192   FORMAT(1X,I5,   44X,A16,' IS=',I1)
1193   FORMAT(1X,I5,   44X,A16)
1195   FORMAT(1X,'  ----------------------------- ')
      
       read(OTXTA,*) cwork
       open (unit=atomf, file='atoms/'//trim(cwork),
     &      action='write')
       write(atomf,*) '#  AES results:'
       write(atomf,993) NL, '# NL'
       write(atomf,*) '#  No QSCR at this stage'
       WRITE(atomf,997) AZ,WS,WS,'# AZ, WSR'
       write(atomf,995) 0.0d0, '# Fermi Energy'         
       write(atomf,996) 0.0d0,'# Pot.Shift'
       write(atomf,993) 0,'# SW'
       write(atomf,993) NS,'# NS'
       DO IS=1,NS
          WRITE(atomf,992) NSZRAD,OTXTA,IS
          WRITE(atomf,994) (POT(I,IS),I=1,NSZRAD)
       end do
       write(atomf,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
     $      
       DO IS=1,NS
          WRITE(atomf,1191) OTXTA,IS
          WRITE(atomf,1104) (ENY(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (PPC(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (PPD(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (PPQ(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (PPP(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (DNY(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (FINY(IL,IS),IL=1,NL)
          WRITE(atomf,1104) (FINYD(IL,IS),IL=1,NL)
       end do
       NOBC=0
       WRITE(atomf,1193) MCOR,OTXTA
       IF(MCOR.NE.0) then
          DO J=1,MCOR
             WRITE(atomf,1101) NCOR(J),LCOR(J),NOBC
             WRITE(atomf,1104) ECOR(J,1),ECOR(J,NS)
          end do
       end if
       close (unit=atomf)       
       return
       end
      
      
C********************
CXXX    ATOMA    ****
C********************
      SUBROUTINE ATOMA
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT CHARACTER*4 (O)
C
      PARAMETER(MNR=512)
C
C--------------------------------------------------
C        READING OF INPUT DATA AND INITIALIZATION
C--------------------------------------------------
C
      DIMENSION OTXT(20)
C
      COMMON/IDAT/ AZ,WS,VOCC(4,2),NL,NS,IREL,NR,NVAL(4)
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
      COMMON/COR/ ECOR(20,2),THRESH,NCOR(20),LCOR(20),MCOR
      COMMON/PP/ ENY(4,2),PPC(4,2),PPD(4,2),PPQ(4,2),
     &           PPP(4,2),DNY(4,2),FINY(4,2),FINYD(4,2)
      COMMON/PTO/ POT(MNR,2)
      COMMON/RW/ IR5,IW6,IW8
      COMMON/ITER/ DEV(2),DEVM,ALFA,IITER,NITER
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RCH/0.5D0/
C
      DATA EPSSCF/1.0D-7/,EPSCOR/1.0D-8/
C
100   FORMAT(20A4)
101   FORMAT(1X,10I5)
105   FORMAT(1X,5F10.5)
C                                 I/O UNITS
      IR5=5
      IW6=6
      IW8=8
C
      OPEN(UNIT=IR5,FILE='aes.inp',FORM='FORMATTED')
      OPEN(UNIT=IW6,FILE='aes.ot6',FORM='FORMATTED')
      OPEN(UNIT=IW8,FILE='aes.ot8',FORM='FORMATTED')
C
C                                 PRECISION
      DEVM=EPSSCF
      THRESH=EPSCOR
C                                        READING
      READ(IR5,100) OTXT
      READ(IR5,105) AZ,WS
      READ(IR5,101) MCOR
      IF(MCOR.GT.20) STOP ' *** MCOR.GT.20 *** '
      IF(MCOR.EQ.0) GO TO 299
      DO 301 J=1,MCOR
      READ(IR5,101) NCOR(J),LCOR(J)
301   CONTINUE
299   CONTINUE
      READ(IR5,101) NL,NS,IREL,NR
      IF(NL.GT.4) STOP ' *** NL.GT.4 *** '
      IF(NR.GT.MNR) STOP ' *** NR.GT.MNR *** '
      DO 302 IL=1,NL
      READ(IR5,105) (VOCC(IL,IS),IS=1,NS)
302   CONTINUE
      READ(IR5,101) NITER
      READ(IR5,105) ALFA
C                                         PRINT
      WRITE(IW6,110)
110   FORMAT(//5X,'   AES   -  INPUT DATA : '//)
      WRITE(IW6,100) OTXT
      WRITE(IW6,111) AZ,WS
111   FORMAT(/9X,' AT. NUMBER=',F6.2,'     W.S. RADIUS=',F8.5)
      WRITE(IW6,112) MCOR
112   FORMAT(/1X,' --------- CORE ORBITALS ----------'/
     &       /9X,' NUMBER OF CORE ORBITALS=',I5/)
      IF(MCOR.EQ.0) GO TO 298
      DO 310 J=1,MCOR
      WRITE(IW6,113) J,NCOR(J),LCOR(J)
310   CONTINUE
298   CONTINUE
113   FORMAT(9X,'  SHELL NO.',I2,'   N=',I2,'   L=',I2)
      WRITE(IW6,114) NL,NS,IREL,NR
114   FORMAT(/1X,' --------- VALENCE ORBITALS ---------'/
     & /9X,'NL(=LMAX+1)=',I2,8X,'NO. OF SPIN DIRECTIONS=',I2/
     & /9X,'IREL=',I2,6X,'(0: NON-RELAT., 1: SCALAR-RELAT.)'/
     & /9X,' SIZE OF RADIAL MESH UP TO W.S. RADIUS=',I5/)
      DO 311 IL=1,NL
      WRITE(IW6,115) IL,(VOCC(IL,IS),IS=1,NS)
311   CONTINUE
115   FORMAT(1X,' IL(=L+1)=',I2,
     &          '   VALENCE OCCUPATION PER SPIN=',2F10.5)
      WRITE(IW6,119) NITER
119   FORMAT(/1X,' -------- SELFCONSISTENCY ----------'/
     &  /9X,'  MAX. NUMBER OF ITER.=',I5)
      WRITE(IW6,116) ALFA
116   FORMAT(/9X,'  MIXING PARAMETER ALFA=',F10.5)
      WRITE(IW6,117) DEVM
117   FORMAT(/9X,'  MAX. DEVIATION OF R*V(R)=',E10.3)
      WRITE(IW6,118) THRESH
118   FORMAT(/9X,' ABSOLUTE PRECISION OF ENERGIES=',E10.3)
C
C                                     TEST OF NEUTRALITY
      AION=AZ
      IF(MCOR.EQ.0) GO TO 297
      DO 317 J=1,MCOR
      AL=REAL(LCOR(J))
      AION=AION-RC2*(RC2*AL+RC1)
317   CONTINUE
297   CONTINUE
      DO 318 IL=1,NL
      AION=AION-VOCC(IL,1)-VOCC(IL,NS)
318   CONTINUE
C
      WRITE(IW6,129) AION
129   FORMAT(//1X,'  **** IONICITY=',F12.5)
C
C      ----------------------------     STARTING VALUES
C
C                         PRINCIPAL QUANTUM
C                         NUMBERS OF VALENCE ORBITALS
      DO 350 IL=1,NL
      L=IL-1
      NVAL(IL)=IL
      IF(MCOR.EQ.0) GO TO 350
      DO 351 J=1,MCOR
      IF(LCOR(J).NE.L) GO TO 351
      IF(NCOR(J).GE.NVAL(IL)) NVAL(IL)=NCOR(J)+1
351   CONTINUE
350   CONTINUE
C                                PRINT
      WRITE(IW6,160)
160   FORMAT(/1X,' *** VALENCE ORBITALS -',
     &           ' PRINCIPAL QUANTUM NUMBERS:')
      DO 355 IL=1,NL
      WRITE(IW6,161) IL,NVAL(IL)
355   CONTINUE
161   FORMAT(9X,'  IL(=L+1)=',I2,'     NVAL(IL)=',I2)
C
C                            VALENCE ENERGIES
      DO 340 IS=1,NS
      DO 341 IL=1,NL
      ENY(IL,IS)=RCZ
341   CONTINUE
340   CONTINUE
C                                  RADIAL MESH
      CALL RAPO(NR,WS,R)
C                                   POTENTIALS
      ANV=REAL(NVAL(1))
      AHU=RCH*ANV**2/AZ
      BHU=RC2*AZ/AHU
C
      POT(1,1)=RCZ
      POT(1,NS)=RCZ
      DO 320 I=2,NR
      DUM=-BHU/(EXP(R(I)/AHU)-RC1)
      DO 321 IS=1,NS
      POT(I,IS)=DUM
321   CONTINUE
320   CONTINUE
C                                   CORE ENERGIES
      IF(MCOR.EQ.0) GO TO 295
      DO 330 J=1,MCOR
      AN=REAL(NCOR(J))
      DUM=AZ/AN-RCH*AN/AHU
      DO 331 IS=1,NS
      ECOR(J,IS)=-DUM**2
331   CONTINUE
330   CONTINUE
C                                   PRINT
      WRITE(IW6,150)
150   FORMAT(/1X,' ***  CORE ORBITALS - ESTIMATED EIGENVALUES:')
      DO 338 J=1,MCOR
      WRITE(IW6,151) NCOR(J),LCOR(J),(ECOR(J,IS),IS=1,NS)
338   CONTINUE
151   FORMAT(5X,' N=',I2,'    L=',I2,'     E= ',2G15.7)
295   CONTINUE
C
      RETURN
      END
C********************
CXXX    ATOMC    ****
C********************
      SUBROUTINE ATOMC
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C---------------------------------------
C        CORE DENSITIES
C---------------------------------------
C
      COMMON/IDAT/ AZ,WS,VOCC(4,2),NL,NS,IREL,NR,NVAL(4)
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
      COMMON/COR/ ECOR(20,2),THRESH,NCOR(20),LCOR(20),MCOR
      COMMON/ROH/ RHOCOR(MNR,2),RHOVAL(MNR,2)
      COMMON/PTO/ POT(MNR,2)
      COMMON/ITER/ DEV(2),DEVM,ALFA,IITER,NITER
      COMMON/RW/ IR5,IW6,IW8
C
      DATA C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      DO 300 IS=1,NS
      DO 301 I=1,NR
      RHOCOR(I,IS)=RCZ
301   CONTINUE
300   CONTINUE
C
      IF(MCOR.EQ.0) RETURN
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
      UC= REAL(IREL)/C274
      UCSQ=UC**2
C
C                                          LOOP OVER SPIN
      DO 302 IS=1,NS
C                                     POTENTIAL
      DO 303 I=2,NR
      V(I)=POT(I,IS)
303   CONTINUE
      CALL LIPO(AZ,NR)
C                                  LOOP OVER CORE ORBITALS
      DO 305 J=1,MCOR
C
      E=ECOR(J,IS)
      L=LCOR(J)
      N=NCOR(J)
      CONO= REAL(2*L+1)/PI4
C
      CALL RSEC(AZ,E,THRESH,N,L,IREL,NR)
C
      ECOR(J,IS)=E
      BLAM=REAL(L*(L+1))
      RHOCOR(1,IS)=RHOCOR(1,IS) + CONO*WG(1)**2
      DO 310 I=2,NR
      AME=RC1 + UCSQ*(E-V(I))
      FRAC=WG(I)/(R(I)*AME)
      DENS=WG(I)**2 + UCSQ * (WF(I)**2 + BLAM*FRAC**2)
      RHOCOR(I,IS)=RHOCOR(I,IS) + CONO*DENS
310   CONTINUE
C
305   CONTINUE
302   CONTINUE
C
C                                              PRINT
      WRITE(IW6,110) IITER
110   FORMAT(/1X,'ITER=',I3,' - CORE ENERGIES :')
      DO 350 J=1,MCOR
      WRITE(IW6,112) NCOR(J),LCOR(J),(ECOR(J,IS),IS=1,NS)
350   CONTINUE
112   FORMAT(1X,'  N=',I2,'    L=',I2,'     E= ',2G15.7)
C
      RETURN
      END
C********************
CXXX    ATOMD    ****
C********************
      SUBROUTINE ATOMD
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C----------------------------------
C         VALENCE DENSITIES
C----------------------------------
C
      COMMON/IDAT/ AZ,WS,VOCC(4,2),NL,NS,IREL,NR,NVAL(4)
      COMMON/ROH/ RHOCOR(MNR,2),RHOVAL(MNR,2)
      COMMON/COR/ ECOR(20,2),THRESH,NCOR(20),LCOR(20),MCOR
      COMMON/PP/ ENY(4,2),PPC(4,2),PPD(4,2),PPQ(4,2),
     &           PPP(4,2),DNY(4,2),FINY(4,2),FINYD(4,2)
      COMMON/PTO/ POT(MNR,2)
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
      COMMON/RW/ IR5,IW6,IW8
      COMMON/ITER/ DEV(2),DEVM,ALFA,IITER,NITER
C
      DATA C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
      UC= REAL(IREL)/C274
      UCSQ=UC**2
C
C                                     LOOP OVER SPINS
      DO 300 IS=1,NS
C
      DO 301 I=1,NR
      RHOVAL(I,IS)=RCZ
301   CONTINUE
C                                      POTENTIAL
      DO 302 I=2,NR
      V(I)=POT(I,IS)
302   CONTINUE
      CALL LIPO(AZ,NR)
C                                     LOOP OVER L
      DO 304 IL=1,NL
C
      E=ENY(IL,IS)
      L=IL-1
      N=NVAL(IL)
      CONO=VOCC(IL,IS)/PI4
C
      CALL RSEC(AZ,E,THRESH,N,L,IREL,NR)
C
      ENY(IL,IS)=E
      BLAM=REAL(L*(L+1))
      RHOVAL(1,IS)=RHOVAL(1,IS) + CONO*WG(1)**2
      DO 308 I=2,NR
      AME=RC1 + UCSQ*(E-V(I))
      FRAC=WG(I)/(R(I)*AME)
      DENS=WG(I)**2 + UCSQ * (WF(I)**2 + BLAM*FRAC**2)
      RHOVAL(I,IS)=RHOVAL(I,IS) + CONO*DENS
308   CONTINUE
C
304   CONTINUE
300   CONTINUE
C
C                                           PRINT
      WRITE(IW6,110) IITER
110   FORMAT(/1X,'ITER=',I3,' - VALENCE ENERGIES :')
      DO 350 IL=1,NL
      WRITE(IW6,112) IL,(ENY(IL,IS),IS=1,NS)
350   CONTINUE
112   FORMAT(5X,' IL=',I2,'     E= ',2G15.7)
C
      RETURN
      END
C********************
CXXX    ATOMF    ****
C********************
      SUBROUTINE ATOMF
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C---------------------------------------------------
C    NEW POTENTIALS FOR SELF-CONSISTENT CALCULATION
C           AND OUTPUT
C---------------------------------------------------
C
      DIMENSION VCOUL(MNR),VTOT(MNR,2),RHO(MNR,2),Y(MNR),
     &          RHOT(MNR),VXC(2),IRMAX(2)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
      COMMON/IDAT/ AZ,WS,VOCC(4,2),NL,NS,IREL,NR,NVAL(4)
      COMMON/PP/ ENY(4,2),PPC(4,2),PPD(4,2),PPQ(4,2),
     &           PPP(4,2),DNY(4,2),FINY(4,2),FINYD(4,2)
      COMMON/RW/ IR5,IW6,IW8
      COMMON/COR/ ECOR(20,2),THRESH,NCOR(20),LCOR(20),MCOR
      COMMON/PTO/ POT(MNR,2)
      COMMON/ITER/ DEV(2),DEVM,ALFA,IITER,NITER
      COMMON/ROH/ RHOCOR(MNR,2),RHOVAL(MNR,2)
      COMMON/LIN/ ENY1,PPC1,PPD1,PPQ1,PPP1,DNY1,FINY1,FINYD1,
     &            FI(MNR),FID(MNR),FIDD(MNR)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
C                                 SPIN DEPENDENT DENSITY
      DO 310 IS=1,NS
      DO 311 I=1,NR
      RHO(I,IS)=RHOCOR(I,IS)+RHOVAL(I,IS)
311   CONTINUE
310   CONTINUE
C                                  TOTAL CHARGE DENSITY
!SMP$ DO SERIAL
      DO 312 I=1,NR
      RHOT(I)=RHO(I,1)+RHO(I,NS)
312   CONTINUE
C                                      POISSON EQUATION
      CALL POIS1(AZ,NR,R,RHOT,VCOUL)
C                                       TOTAL POTENTIAL
      DO 315 I=2,NR
      RHOUP=RHO(I,1)
      RHODO=RHO(I,NS)
      CALL XCVBH(RHOUP,RHODO,VXC(1),VXC(2),EXC)
      DO 316 IS=1,NS
      VTOT(I,IS)=VCOUL(I)+VXC(IS)
316   CONTINUE
315   CONTINUE
C
C                                           DEVIATION
      DO 325 IS=1,NS
      Y(1)=RCZ
!SMP$ DO SERIAL      
      DO 326 I=2,NR
      Y(I)=R(I)*ABS(VTOT(I,IS)-POT(I,IS))
326   CONTINUE
      YMX=RCZ
      IMX=0
      DO 327 I=1,NR
      IF(Y(I).LT.YMX) GO TO 327
      YMX=Y(I)
      IMX=I
327   CONTINUE
      DEV(IS)=YMX
      IRMAX(IS)=IMX
325   CONTINUE
C                                          PRINT
      WRITE(IW6,130)
130   FORMAT(1X)
      DO 331 IS=1,NS
      WRITE(IW6,131) IITER,IS,DEV(IS),IRMAX(IS)
331   CONTINUE
131   FORMAT(1X,'-----ITER=',I3,'   IS=',I2,
     &      '  MAX. DEV.=',G11.3,'   RAD. POINT=',I4)
C
      DEVY= MAX(DEV(1),DEV(NS))
      IF(IITER.EQ.NITER.OR.DEVY.LE.DEVM) GO TO 244
C
C                                    MIXING AND RETURN
      BETA=RC1-ALFA
      DO 335 IS=1,NS
      POT(1,IS)=RCZ
      DO 336 I=2,NR
      POT(I,IS) = ALFA*VTOT(I,IS) + BETA*POT(I,IS)
336   CONTINUE
335   CONTINUE
C
      RETURN
C
C********************************* LINEARIZATION
C
244   CONTINUE
C                                  LOOP OVER SPIN
      DO 340 IS=1,NS
C                                      POTENTIAL
      DO 341 I=2,NR
      V(I)=POT(I,IS)
341   CONTINUE
      CALL LIPO(AZ,NR)
C                                  LOOP OVER L
      DO 345 IL=1,NL
C
      L=IL-1
      ENY1=ENY(IL,IS)
      DE=RC1/(REAL(10)*REAL(L+2))
C
      CALL RSEL(AZ,DE,WS,L,IREL,NR)
C
      PPC(IL,IS)=PPC1
      PPD(IL,IS)=PPD1
      PPQ(IL,IS)=PPQ1
      PPP(IL,IS)=PPP1
      DNY(IL,IS)=DNY1
      FINY(IL,IS)=FINY1
      FINYD(IL,IS)=FINYD1
C
345   CONTINUE
340   CONTINUE
C
C ********************************   PRINT OF RESULTS
C
      WRITE(IW6,160)
160   FORMAT(//15X,' ****** POTENTIAL PARAMETERS  : ')
      DO 355 IS=1,NS
      WRITE(IW6,161) IS
161   FORMAT(/1X,'  SPIN : ',I3/)
      WRITE(IW6,162)(ENY(IL,IS),IL=1,NL)
162   FORMAT(1X,' ENY : ',4G15.7)
      WRITE(IW6,163)(PPC(IL,IS),IL=1,NL)
163   FORMAT(1X,'  C  : ',4G15.7)
      WRITE(IW6,164)(PPD(IL,IS),IL=1,NL)
164   FORMAT(1X,'DELTA: ',4G15.7)
      WRITE(IW6,165)(PPQ(IL,IS),IL=1,NL)
165   FORMAT(1X,'  Q  : ',4G15.7)
      WRITE(IW6,166)(PPP(IL,IS),IL=1,NL)
166   FORMAT(1X,'  P  : ',4G15.7)
      WRITE(IW6,167)(DNY(IL,IS),IL=1,NL)
167   FORMAT(1X,' DNY : ',4G15.7)
      WRITE(IW6,168)(FINY(IL,IS),IL=1,NL)
168   FORMAT(1X,'FINY : ',4G15.7)
      WRITE(IW6,169)(FINYD(IL,IS),IL=1,NL)
169   FORMAT(1X,'FINYD: ',4G15.7)
355   CONTINUE
C
      call write_atoms(NL,NS,NR,AZ,WS,POT)
C
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
175   FORMAT(1X,'---------  LSDA-FILE  ----------')
177   FORMAT(1X,'--------------------------------')
C
C -----------------------------    OUTPUT -  UNIT IW8
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
!SMP$ DO SERIAL      
      DO 370 I=1,NR
      Y(I)=PI4*R(I)**2*RHOT(I)
370   CONTINUE
C
      WRITE(IW8,180) AZ,IREL
180   FORMAT(1X,' DENSITY*(4*PI*R**2) FOR  Z= ',F6.2,
     &          '   IREL= ',I1)
      DO 372 I=1,NR
      WRITE(IW8,104) R(I),Y(I)
372   CONTINUE
      WRITE(IW8,177)
C
      RETURN
      END
C*******************
CXXX    PRIM2   ****
C*******************
      SUBROUTINE PRIM2(N,X,Y,F)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
C---------------------------------------------------
C     CALCULATES INDEFINITE INTEGRAL F(I)
C     FOR INTEGRAND Y(I) OF VARIABLE X(I)
C     (WHERE I=1,2,... N, AND N.GE.4)
C     USING CUBIC LAGRANGE INTERPOLATION POLYNOMIAL
C---------------------------------------------------
      DIMENSION X(N),Y(N),F(N)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/,RC4/4.0D0/,RC6/6.0D0/
C
      AUXF(T0,T1,T2,T3,TT)=
     &    (TT-T1)*(TT-T2)*(TT-T3)/((T0-T1)*(T0-T2)*(T0-T3))
C
      F(1)=RCZ
      DO 301 I=2,N
      I0=I-3
      IF(I.EQ.2) I0=0
      IF(I.EQ.N) I0=N-4
      X1=X(I0+1)
      X2=X(I0+2)
      X3=X(I0+3)
      X4=X(I0+4)
C                             LAGRANGE INTERPOLATION
      XT=RCH*(X(I-1)+X(I))
      YT=Y(I0+1)*AUXF(X1,X2,X3,X4,XT)
     &  +Y(I0+2)*AUXF(X2,X3,X4,X1,XT)
     &  +Y(I0+3)*AUXF(X3,X4,X1,X2,XT)
     &  +Y(I0+4)*AUXF(X4,X1,X2,X3,XT)
C                                       SIMPSON RULE
      F(I)=F(I-1)
     &    +(Y(I-1)+RC4*YT+Y(I))*(X(I)-X(I-1))/RC6
301   CONTINUE
      RETURN
      END
C*******************
CXXX    QUAD2   ****
C*******************
      REAL*8 FUNCTION QUAD2(K,N,X,Y)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
C---------------------------------------------------
C     QUADRATURE OF INTEGRAND Y(I) OF VARIABLE X(I)
C       (WHERE K.LE.I.LE.N  AND  N.GE.K+3)
C     USING CUBIC LAGRANGE INTERPOLATION POLYNOMIAL
C---------------------------------------------------
      DIMENSION X(N),Y(N)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/,RC4/4.0D0/,RC6/6.0D0/
C
      AUXF(T0,T1,T2,T3,TT)=
     &    (TT-T1)*(TT-T2)*(TT-T3)/((T0-T1)*(T0-T2)*(T0-T3))
C
      QUAD2=RCZ
      DO 301 I=K+1,N
      I0=I-3
      IF(I.EQ.K+1) I0=K-1
      IF(I.EQ.N) I0=N-4
      X1=X(I0+1)
      X2=X(I0+2)
      X3=X(I0+3)
      X4=X(I0+4)
C                             LAGRANGE INTERPOLATION
      XT=RCH*(X(I-1)+X(I))
      YT=Y(I0+1)*AUXF(X1,X2,X3,X4,XT)
     &  +Y(I0+2)*AUXF(X2,X3,X4,X1,XT)
     &  +Y(I0+3)*AUXF(X3,X4,X1,X2,XT)
     &  +Y(I0+4)*AUXF(X4,X1,X2,X3,XT)
C                                       SIMPSON RULE
      QUAD2=QUAD2
     &    +(Y(I-1)+RC4*YT+Y(I))*(X(I)-X(I-1))/RC6
301   CONTINUE
      RETURN
      END
C*******************
CXXX    RAPO    ****
C*******************
      SUBROUTINE RAPO(NR,WS,R)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
C-------------------------------------
C   GENERATES RADIAL POINTS
C-------------------------------------
C     R(1)=0.0 , R(2)=H , R(NR)=WS ,
C     R(I) = H * (I-1) * Q**(I-2)
C     H - A SMALL QUANTITY, H=2.0E-5
C     NR - NUMBER OF POINTS
C     WS - WIGNER-SEITZ RADIUS
C-------------------------------------
      DIMENSION R(NR)
C
      DATA RCZ/0.0D0/,H/2.0D-5/
C
      ARG=WS/(H*REAL(NR-1))
      DUM=LOG(ARG)/REAL(NR-2)
      Q=EXP(DUM)
C
      R(1)=RCZ
      R(2)=H
      HQ=H
      DO 301 I=3,NR-1
      HQ=HQ*Q
      R(I)=REAL(I-1)*HQ
301   CONTINUE
      R(NR)=WS
      RETURN
      END
C*******************
CXXX    LIPO    ****
C*******************
      SUBROUTINE LIPO(AZ,NR)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C-------------------------------------------------------------
C   LAGRANGE INTERPOLATION OF POTENTIAL FOR THE RUNGE-KUTTA
C   INTEGRATION OF THE RADIAL SCHROEDINGER EQUATION
C-------------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     NR - SIZE OF RADIAL MESH
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C  OUTPUT:
C     VI(.) - INTERPOLATED POTENTIAL: THE VALUE VI(I) REFERS
C     TO THE RADIUS   RM = 0.5*(R(I-1)+R(I)), I = 2,... NR.
C-------------------------------------------------------------
C
      DIMENSION X(4),Y(4)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
C
      DATA RCZ/0.0D0/,RC2/2.0D0/,RCH/0.5D0/
C
      Y(1)=-RC2*AZ
      VI(2)=(Y(1)+V(2)*R(2))/R(2)
      DO 310 I=3,NR
      XM=RCH*(R(I)+R(I-1))
      I0=I-3
      IF(I.EQ.NR) I0=I-4
      X(1)=R(I0+1)
      X(2)=R(I0+2)
      X(3)=R(I0+3)
      X(4)=R(I0+4)
      IF(I.GT.3) Y(1)=X(1)*V(I0+1)
      Y(2)=X(2)*V(I0+2)
      Y(3)=X(3)*V(I0+3)
      Y(4)=X(4)*V(I0+4)
C                             LAGRANGE INTERPOLATION
      FM=RCZ
      DO 321 J=1,4
      P=Y(J)
      DO 322 K=1,4
      IF(K.EQ.J) GO TO 322
      P=P*(XM-X(K))/(X(J)-X(K))
322   CONTINUE
      FM=FM+P
321   CONTINUE
C
      VI(I)=FM/XM
310   CONTINUE
      RETURN
      END
C*******************
CXXX    RSEV    ****
C*******************
      SUBROUTINE RSEV(AZ,E,E0,L,IREL,NR)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C-----------------------------------------------------------------
C   SOLUTION OF RADIAL SCHROEDINGER EQUATION FOR VALENCE ELECTRONS
C        BOTH  NON-RELATIVISTIC (IREL=0) AND
C        SCALAR-RELATIVISTIC (IREL=1) VERSION
C-----------------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     E - ENERGY
C     E0 - ENERGY FOR DOWNFOLDING OF SMALL COMPONENTS
C     L - ORBITAL QUANTUM NUMBER
C     IREL - RELATIVITY
C     NR - SIZE OF RADIAL MESH
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     VI(.) - INTERPOLATED POTENTIAL
C  OUTPUT:
C     WG(.) - WAVE FUNCTION NORMALIZED TO UNITY
C             (GREAT COMPONENT)
C     WF(.) - FUNCTION RELATED TO RADIAL DERIVATIVE OF WG
C             (SMALL COMPONENT)
C-----------------------------------------------------------------
C
      DIMENSION AME(MNR),U(MNR),P(MNR),Q(MNR)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
C
      DATA NRST/10/,C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/,
     &     RCH/0.5D0/,RC6/6.0D0/
C
C                             SET UP CONSTANTS
      LP1=L+1
      AL= REAL(L)
      ALP1= REAL(LP1)
      BLAM=AL*ALP1
      TWOZ=RC2*AZ
C
      UC= REAL(IREL)/C274
      UCSQ=UC**2
      DO 303 I=2,NR
      AME(I)=RC1 + UCSQ*(E0-V(I))
303   CONTINUE
      DO 306 I=2,NR
      U(I)=BLAM/(R(I)**2*AME(I)) + V(I) - E
306   CONTINUE
C
C                            START OF INTEGRATION
      P(1)=RCZ
      Q(1)=RCZ
C                              NON-RELATIVISTIC
      IF(IREL.EQ.0) THEN
      U0=V(2)+TWOZ/R(2)-E
      A1=-AZ/ALP1
      A2=(TWOZ*AZ/ALP1+U0)/(RC4*AL+RC6)
      B1=-AZ
      B2=A2*(AL+RC2)
      P(2)=RC1+R(2)*(A1+R(2)*A2)
      Q(2)=B1+R(2)*B2
      BETA=RC1
       IF(L.GT.0) THEN
       P(2)=R(2)*P(2)
       Q(2)=AL+R(2)*Q(2)
       BETA=AL
       END IF
      END IF
C                           SCALAR-RELATIVISTIC
      IF(IREL.EQ.1) THEN
      BETA=SQRT(BLAM+RC1-(TWOZ*UC)**2)
      P(2)=RC1
      Q(2)=(BETA-RC1)/(TWOZ*UCSQ)
      END IF
C
C                       STARTING RUNGE-KUTTA INTEGRATION
      IF(NRST.EQ.2) GO TO 222
      A1PB=RC1+BETA
      A1MB=RC1-BETA
      DO 301 I=3,NRST
      H=R(I)-R(I-1)
      H2=H/RC2
      H6=H/RC6
C
      RM=RCH*(R(I-1)+R(I))
      VM=VI(I)
      AMEM=RC1 + UCSQ*(E0-VM)
      UM=BLAM/(RM**2*AMEM) + VM - E
C
      PK1=A1MB*P(I-1)/R(I-1) + AME(I-1)*Q(I-1)
      QK1=U(I-1)*P(I-1) - A1PB*Q(I-1)/R(I-1)
      PDUM=P(I-1)+H2*PK1
      QDUM=Q(I-1)+H2*QK1
      PK2=A1MB*PDUM/RM  + AMEM *QDUM
      QK2=UM *PDUM - A1PB*QDUM/RM
      PDUM=P(I-1)+H2*PK2
      QDUM=Q(I-1)+H2*QK2
      PK3=A1MB*PDUM/RM  + AMEM *QDUM
      QK3=UM *PDUM - A1PB*QDUM/RM
      PDUM=P(I-1)+H*PK3
      QDUM=Q(I-1)+H*QK3
      PK4=A1MB*PDUM/R(I) + AME(I)*QDUM
      QK4=U(I)*PDUM - A1PB*QDUM/R(I)
      P(I)=P(I-1)+H6*(PK1+RC2*(PK2+PK3)+PK4)
      Q(I)=Q(I-1)+H6*(QK1+RC2*(QK2+QK3)+QK4)
301   CONTINUE
222   DO 310 I=2,NRST
      RTB=R(I)**BETA
      P(I)=RTB*P(I)
      Q(I)=RTB*Q(I)
310   CONTINUE
C
C                             RUNGE-KUTTA INTEGRATION
      DO 320 I=NRST+1,NR
      H=R(I)-R(I-1)
      H2=H/RC2
      H6=H/RC6
C
      RM=RCH*(R(I-1)+R(I))
      VM=VI(I)
      AMEM=RC1 + UCSQ*(E0-VM)
      UM=BLAM/(RM**2*AMEM) + VM - E
C
      PK1=P(I-1)/R(I-1) + AME(I-1)*Q(I-1)
      QK1=U(I-1)*P(I-1) - Q(I-1)/R(I-1)
      PDUM=P(I-1)+H2*PK1
      QDUM=Q(I-1)+H2*QK1
      PK2=PDUM/RM  + AMEM *QDUM
      QK2=UM *PDUM - QDUM/RM
      PDUM=P(I-1)+H2*PK2
      QDUM=Q(I-1)+H2*QK2
      PK3=PDUM/RM  + AMEM *QDUM
      QK3=UM *PDUM - QDUM/RM
      PDUM=P(I-1)+H*PK3
      QDUM=Q(I-1)+H*QK3
      PK4=PDUM/R(I) + AME(I)*QDUM
      QK4=U(I)*PDUM - QDUM/R(I)
      P(I)=P(I-1)+H6*(PK1+RC2*(PK2+PK3)+PK4)
      Q(I)=Q(I-1)+H6*(QK1+RC2*(QK2+QK3)+QK4)
320   CONTINUE
C
C                                  NORMALIZATION
      DO 330 I=1,NR
      U(I)=P(I)**2
330   CONTINUE
      SUM=QUAD2(1,NR,R,U)
      CNORM=SQRT(SUM)
C
      WG(1)=RCZ
      WF(1)=RCZ
      IF(IREL.EQ.0) THEN
       IF(L.EQ.0) THEN
         WG(1)=RC1/CNORM
         WF(1)=-AZ/CNORM
       END IF
       IF(L.EQ.1) WF(1)=RC1/CNORM
      END IF
C
      DO 333 I=2,NR
      DUM=CNORM*R(I)
      WG(I)=P(I)/DUM
      WF(I)=Q(I)/DUM
333   CONTINUE
C
      RETURN
      END
C*******************
CXXX    RSEL    ****
C*******************
      SUBROUTINE RSEL(AZ,EH,WSA,L,IREL,NR)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C-----------------------------------------------------------------
C     LINEARIZATION OF ENERGY DEPENDENCE OF SOLUTION
C     OF RADIAL SCHROEDINGER EQUATION FOR VALENCE ELECTRONS
C        BOTH  NON-RELATIVISTIC (IREL=0) AND
C        SCALAR-RELATIVISTIC (IREL=1) VERSION
C-----------------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     EH - ENERGY STEP FOR NUMERICAL DERIVATIVE
C     WSA - AVERAGE WIGNER-SEITZ RADIUS
C           (TO SCALE THE POTENTIAL PARAMETERS)
C     L - ORBITAL QUANTUM NUMBER
C     IREL - RELATIVITY
C     NR - SIZE OF RADIAL MESH
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     VI(.) - INTERPOLATED POTENTIAL
C     ENY - ENERGY VALUE
C  OUTPUT:
C     PPC, ..., FINYD - POTENTIAL PARAMETERS
C     FI(.),FID(.),FIDD(.) - WAVEFUNCTION (NORMALIZED TO UNITY)
C                            AND ITS TWO ENERGY DERIVATIVES
C-----------------------------------------------------------------
C
      DIMENSION SN(MNR,5),B(5),DUM(MNR)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
      COMMON/LIN/ ENY,PPC,PPD,PPQ,PPP,DNY,FINY,FINYD,
     &            FI(MNR),FID(MNR),FIDD(MNR)

     
C
      DATA RC1/1.0D0/,RC2/2.0D0/,RC8/8.0D0/,
     &     RC12/12.0D0/,RC16/16.0D0/,RC30/30.0D0/
C
      AL= REAL(L)
      ALP1=AL+RC1
      WS=R(NR)
      FAK=(WS/WSA)**(2*L+1)
C
      E0=ENY
C                                  CYKLUS PRES 5 ENERGII
      DO 305 IE=1,5
      E=ENY+REAL(IE-3)*EH
      CALL RSEV(AZ,E,E0,L,IREL,NR)
      IF(IE.EQ.3) DNY=WS*WF(NR)/WG(NR)
      DO 306 I=1,NR
      SN(I,IE)=WG(I)
306   CONTINUE
305   CONTINUE

C                                   VYPOCET FI,FID,FIDD
      H1=RC12*EH
      H2=H1*EH
C
      DO 313 I=1,NR
      DO 317 IE=1,5
      B(IE)=SN(I,IE)
317   CONTINUE
C
      A0=B(3)
      S1=B(2)+B(4)
      S2=B(1)+B(5)
      R1=B(4)-B(2)
      R2=B(5)-B(1)
C
      FI(I)=A0
      FID(I)=(RC8*R1-R2)/H1
      FIDD(I)=(RC16*S1-S2-RC30*A0)/H2
313   CONTINUE
C                            VYPOCET POTENCIAL. PARAMETRU
C                                  FINY,FINYD,DNY,DNYD
      FINY=FI(NR)
      FINYD=FID(NR)
      DNYD=DNY-RC1/(WS*FINY*FINYD)
C                                     PPC,PPD,PPQ
      AJM=DNYD+ALP1
C
      PPC=ENY - FINY*(DNY+ALP1)/(FINYD*AJM)
      PPD=FAK/(RC2*WS*(FINYD*AJM)**2)
      PPQ=FAK*(DNYD-AL)/(RC2*(AL+ALP1)*AJM)
C            
C                                PPP

!SMP$ DO SERIAL
       DO I=1,NR
       DUM(I)=(R(I)*FID(I))**2
       END DO
      PPP=QUAD2(1,NR,R,DUM)
C
      RETURN
      END
C*******************
CXXX    RSEC    ****
C*******************
      SUBROUTINE RSEC(AZ,E,THRESH,N,L,IREL,NR)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C-----------------------------------------------------------------
C   SOLUTION OF RADIAL SCHROEDINGER EQUATION FOR CORE ELECTRONS
C        BOTH  NON-RELATIVISTIC (IREL=0) AND
C        SCALAR-RELATIVISTIC (IREL=1) VERSION
C-----------------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     E - ENERGY (ESTIMATE)
C     THRESH - ABSOLUTE ACCURACY OF THE EIGENVALUE
C     N - PRINCIPAL QUANTUM NUMBER
C     L - ORBITAL QUANTUM NUMBER
C     IREL - RELATIVITY
C     NR - SIZE OF RADIAL MESH
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     VI(.) - INTERPOLATED POTENTIAL
C  OUTPUT:
C     E - ENERGY
C     WG(.) - WAVE FUNCTION (GREAT COMPONENT)
C     WF(.) - FUNCTION RELATED TO RADIAL DERIVATIVE OF WG
C             (SMALL COMPONENT)
C-----------------------------------------------------------------
C
      DIMENSION AME(MNR),U(MNR),P(MNR),Q(MNR),WRK(MNR)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(MNR),WG(MNR),WF(MNR)
C
      DATA NRST/10/,MAXIT/200/,ARGMAX/120.0D0/,
     &     DSG/1.5D0/,C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/,
     &     RCH/0.5D0/,RC6/6.0D0/
C
C                             SET UP CONSTANTS
      LP1=L+1
      AL= REAL(L)
      ALP1= REAL(LP1)
      BLAM=AL*ALP1
      TWOZ=RC2*AZ
      NUMNOD=N-LP1
      AKLSQ=(ALP1/R(NR))**2
      UC= REAL(IREL)/C274
      UCSQ=UC**2
      IT=0
      IMIN=0
      IMAX=0
C
200   IT=IT+1
C>>>>
      IF(2*IT.EQ.MAXIT) WRITE(6,100) N,L
100   FORMAT(/3X,'-----RSEC: PROBLEMS FOR  N=',I2,'    L=',I2)
      IF(2*IT.GT.MAXIT) WRITE(6,101) IT,E,IMIN,IMAX
101   FORMAT(3X,'---RSEC:  IT=',I3,'     E=',G20.10/
     &        10X,'     IMIN=',I2,'       IMAX=',I2)
C>>>>
      IF(IT.GT.MAXIT) STOP ' *** RSEC - NO CONVERGENCY '
C
      DO 302 I=1,NR
      P(I)=RCZ
      Q(I)=RCZ
302   CONTINUE
      DO 303 I=2,NR
      AME(I)=RC1 + UCSQ*(E-V(I))
      U(I)=BLAM/(R(I)**2*AME(I)) + V(I) - E
303   CONTINUE
C
C                      SEARCH FOR CLASSICAL TURNING POINT
      I=NR+1
205   I=I-1
      IF(U(I).GT.RCZ.AND.I.GT.NRST+10) GO TO 205
      NRMA=I+1
      IF(NRMA+3.GT.NR) NRMA=NR
C>>>>
      IF(2*IT.GT.MAXIT) WRITE(6,105) NRMA
105   FORMAT(10X,'   NRMA=',I3)
C>>>>
C
C -------------------------   OUTWARD INTEGRATION
C                          START -  NON-RELATIVISTIC
      IF(IREL.EQ.0) THEN
      U0=V(2)+TWOZ/R(2)-E
      A1=-AZ/ALP1
      A2=(TWOZ*AZ/ALP1+U0)/(RC4*AL+RC6)
      B1=-AZ
      B2=A2*(AL+RC2)
      P(2)=RC1+R(2)*(A1+R(2)*A2)
      Q(2)=B1+R(2)*B2
      BETA=RC1
       IF(L.GT.0) THEN
       P(2)=R(2)*P(2)
       Q(2)=AL+R(2)*Q(2)
       BETA=AL
       END IF
      END IF
C                         START -  SCALAR-RELATIVISTIC
      IF(IREL.EQ.1) THEN
      BETA=SQRT(BLAM+RC1-(TWOZ*UC)**2)
      P(2)=RC1
      Q(2)=(BETA-RC1)/(TWOZ*UCSQ)
      END IF
C                       STARTING RUNGE-KUTTA INTEGRATION
      IF(NRST.EQ.2) GO TO 208
      A1PB=RC1+BETA
      A1MB=RC1-BETA
      DO 311 I=3,NRST
      H=R(I)-R(I-1)
      H2=H/RC2
      H6=H/RC6
C
      RM=RCH*(R(I-1)+R(I))
      VM=VI(I)
      AMEM=RC1 + UCSQ*(E-VM)
      UM=BLAM/(RM**2*AMEM) + VM - E
C
      PK1=A1MB*P(I-1)/R(I-1) + AME(I-1)*Q(I-1)
      QK1=U(I-1)*P(I-1) - A1PB*Q(I-1)/R(I-1)
      PDUM=P(I-1)+H2*PK1
      QDUM=Q(I-1)+H2*QK1
      PK2=A1MB*PDUM/RM  + AMEM *QDUM
      QK2=UM *PDUM - A1PB*QDUM/RM
      PDUM=P(I-1)+H2*PK2
      QDUM=Q(I-1)+H2*QK2
      PK3=A1MB*PDUM/RM  + AMEM *QDUM
      QK3=UM *PDUM - A1PB*QDUM/RM
      PDUM=P(I-1)+H*PK3
      QDUM=Q(I-1)+H*QK3
      PK4=A1MB*PDUM/R(I) + AME(I)*QDUM
      QK4=U(I)*PDUM - A1PB*QDUM/R(I)
      P(I)=P(I-1)+H6*(PK1+RC2*(PK2+PK3)+PK4)
      Q(I)=Q(I-1)+H6*(QK1+RC2*(QK2+QK3)+QK4)
311   CONTINUE
208   DO 312 I=2,NRST
      RTB=R(I)**BETA
      P(I)=RTB*P(I)
      Q(I)=RTB*Q(I)
312   CONTINUE
C                             RUNGE-KUTTA INTEGRATION
      DO 314 I=NRST+1,NRMA
      H=R(I)-R(I-1)
      H2=H/RC2
      H6=H/RC6
C
      RM=RCH*(R(I-1)+R(I))
      VM=VI(I)
      AMEM=RC1 + UCSQ*(E-VM)
      UM=BLAM/(RM**2*AMEM) + VM - E
C
      PK1=P(I-1)/R(I-1) + AME(I-1)*Q(I-1)
      QK1=U(I-1)*P(I-1) - Q(I-1)/R(I-1)
      PDUM=P(I-1)+H2*PK1
      QDUM=Q(I-1)+H2*QK1
      PK2=PDUM/RM  + AMEM *QDUM
      QK2=UM *PDUM - QDUM/RM
      PDUM=P(I-1)+H2*PK2
      QDUM=Q(I-1)+H2*QK2
      PK3=PDUM/RM  + AMEM *QDUM
      QK3=UM *PDUM - QDUM/RM
      PDUM=P(I-1)+H*PK3
      QDUM=Q(I-1)+H*QK3
      PK4=PDUM/R(I) + AME(I)*QDUM
      QK4=U(I)*PDUM - QDUM/R(I)
      P(I)=P(I-1)+H6*(PK1+RC2*(PK2+PK3)+PK4)
      Q(I)=Q(I-1)+H6*(QK1+RC2*(QK2+QK3)+QK4)
314   CONTINUE
C                                    NUMBER OF NODES
      NUNO=0
      ASIG=SIGN(RC1,P(2))
      DO 316 I=3,NRMA
      BSIG=SIGN(RC1,P(I))
      IF(ABS(ASIG-BSIG).LT.DSG) GO TO 316
      NUNO=NUNO+1
      ASIG=BSIG
316   CONTINUE
      PRODPQ=P(NRMA)*Q(NRMA)
      IF(NUNO.EQ.NUMNOD.AND.PRODPQ.LT.RCZ) GO TO 250
C>>>>
      IF(2*IT.GT.MAXIT) WRITE(6,112) NUNO,PRODPQ
112   FORMAT(10X,'   NUNO=',I2,'     PRODPQ=',G15.7)
C>>>>
C
C                             ENERGY OUT OF CORRECT WINDOW
      IF(NUNO.GT.NUMNOD) GO TO 230
C                                   ENERGY TOO LOW
      IF(IMIN.EQ.0) EMIN=E
      IF(IMIN.EQ.1) EMIN=MAX(EMIN,E)
      IMIN=1
      SIGDE=RC1
      GO TO 240
C                                   ENERGY TOO HIGH
230   IF(IMAX.EQ.0) EMAX=E
      IF(IMAX.EQ.1) EMAX=MIN(EMAX,E)
      IMAX=1
      SIGDE=-RC1
C
240   IPROD=IMIN*IMAX
      IF(IPROD.EQ.0) THEN
      ABSE=ABS(E)
      AMGDE= MAX(ABSE,RC1)/ REAL(10)
      E=E+SIGDE*AMGDE
      END IF
      IF(IPROD.EQ.1) E=RCH*(EMAX+EMIN)
      GO TO 200
C                            ENERGY IN THE CORRECT WINDOW
250   PMAL=P(NRMA)
      XLDL=Q(NRMA)/PMAL
      WRK(1)=RCZ
      DO 322 I=2,NRMA
      FRAC=P(I)/(R(I)*AME(I))
      WRK(I)=P(I)**2 + UCSQ*(Q(I)**2+BLAM*FRAC**2)
322   CONTINUE
      ANIL=QUAD2(1,NRMA,R,WRK)
C
C                      STARTING POINT FOR INWARD INTEGRATION
      XLD1SQ= MAX( U(NR)/AME(NR) , AKLSQ )
      XLD1=-SQRT(XLD1SQ)
      IF(NRMA.LT.NR) GO TO 253
      PMAR=RC1
      DLD1=-RCH/XLD1
      XLDR=XLD1
      ANIR=RCZ
      GO TO 260
253   DO 325 I=NRMA,NR
      WRK(I)=SQRT(U(I)*AME(I))
325   CONTINUE
      SUM=RCZ
      I=NRMA
255   I=I+1
      SUM=SUM+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
      IF(SUM.LT.ARGMAX.AND.I.LT.NR) GO TO 255
      NRII=MAX(I,NRMA+3)
C>>>>
      IF(2*IT.GT.MAXIT) WRITE(6,117) NRII
117   FORMAT(10X,'     NRII=',I3)
C>>>>
C
C ------------------------------------ INWARD INTEGRATION
C                                   START
      P(NRII)=RC1
      IF(NRII.EQ.NR) THEN
       Q(NRII)=XLD1
       DLD1=-RCH/XLD1
      END IF
      IF(NRII.LT.NR) THEN
       Q(NRII)=-WRK(NRII)/AME(NRII)
       DLD1=RCZ
      END IF
C                                    RUNGE-KUTTA INTEGRATION
      DO 330 I=NRII-1,NRMA,-1
      H=R(I)-R(I+1)
      H2=H/RC2
      H6=H/RC6
C
      RM=RCH*(R(I)+R(I+1))
      VM=VI(I+1)
      AMEM=RC1 + UCSQ*(E-VM)
      UM=BLAM/(RM**2*AMEM) + VM - E
C
      XSI=WRK(I)
      EXPO=EXP(-XSI*H)
C
      PK1=(XSI+RC1/R(I+1))*P(I+1) + AME(I+1)*Q(I+1)
      QK1=U(I+1)*P(I+1) + (XSI-RC1/R(I+1))*Q(I+1)
      PDUM=P(I+1)+H2*PK1
      QDUM=Q(I+1)+H2*QK1
      PK2=(XSI+RC1/RM)*PDUM + AMEM *QDUM
      QK2=UM *PDUM + (XSI-RC1/RM)* QDUM
      PDUM=P(I+1)+H2*PK2
      QDUM=Q(I+1)+H2*QK2
      PK3=(XSI+RC1/RM)*PDUM  + AMEM *QDUM
      QK3=UM *PDUM + (XSI-RC1/RM)* QDUM
      PDUM=P(I+1)+H*PK3
      QDUM=Q(I+1)+H*QK3
      PK4=(XSI+RC1/R(I))*PDUM + AME(I)*QDUM
      QK4=U(I)*PDUM + (XSI-RC1/R(I))* QDUM
      P(I)=EXPO*( P(I+1) + H6*(PK1+RC2*(PK2+PK3)+PK4) )
      Q(I)=EXPO*( Q(I+1) + H6*(QK1+RC2*(QK2+QK3)+QK4) )
330   CONTINUE
C
      PMAR=P(NRMA)
      XLDR=Q(NRMA)/PMAR
      DO 332 I=NRMA,NRII
      FRAC=P(I)/(R(I)*AME(I))
      WRK(I)=P(I)**2 + UCSQ*(Q(I)**2+BLAM*FRAC**2)
332   CONTINUE
      ANIR=QUAD2(NRMA,NRII,R,WRK)
C
C                         ENERGY CHANGE FROM PERTURBATION THEORY
260   DE = (XLDL-XLDR) / (ANIL/PMAL**2 + (DLD1+ANIR)/PMAR**2)
C>>>>
      IF(2*IT.GT.MAXIT) WRITE(6,130) XLDL,XLDR,ANIL,
     &                               ANIR,PMAL,PMAR,DE
130   FORMAT(10X,' XLDL=',G20.10,'  XLDR=',G20.10/
     &       10X,' ANIL=',G20.10,'  ANIR=',G20.10/
     &       10X,' PMAL=',G20.10,'  PMAR=',G20.10/
     &       10X,'        DE=',G15.7)
C>>>>
      IF(ABS(DE).LT.THRESH) GO TO 280
C
      IF(DE.GT.RCZ) THEN
       IF(IMIN.EQ.0) EMIN=E
       IF(IMIN.EQ.1) EMIN=MAX(EMIN,E)
       IMIN=1
      END IF
      IF(DE.LT.RCZ) THEN
       IF(IMAX.EQ.0) EMAX=E
       IF(IMAX.EQ.1) EMAX=MIN(EMAX,E)
       IMAX=1
      END IF
      E=E+DE
      GO TO 200
C                               MATCHING OF BOTH PARTS
280   IF(NRMA.EQ.NR) GO TO 290
      DUM=PMAL/PMAR
      DO 370 I=NRMA,NRII
      P(I)=P(I)*DUM
      Q(I)=Q(I)*DUM
370   CONTINUE
C                                  NORMALIZATION
290   WRK(1)=RCZ
      DO 380 I=2,NR
      FRAC=P(I)/(R(I)*AME(I))
      WRK(I)=P(I)**2 + UCSQ*(Q(I)**2+BLAM*FRAC**2)
380   CONTINUE
      SUM=QUAD2(1,NR,R,WRK)
      CNORM=SQRT(SUM)
C
      WG(1)=RCZ
      WF(1)=RCZ
      IF(IREL.EQ.0) THEN
       IF(L.EQ.0) THEN
         WG(1)=RC1/CNORM
         WF(1)=-AZ/CNORM
       END IF
       IF(L.EQ.1) WF(1)=RC1/CNORM
      END IF
C
      DO 385 I=2,NR
      DUM=CNORM*R(I)
      WG(I)=P(I)/DUM
      WF(I)=Q(I)/DUM
385   CONTINUE
C
      RETURN
      END
C*******************
CXXX    POIS1   ****
C*******************
      SUBROUTINE POIS1(AZ,N,R,RHO,V)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
      PARAMETER(MNR=512)
C
C------------------------------------------------------------
C  INTEGRATION OF POISSON'S EQUATION INSIDE 1 ATOMIC SPHERE
C------------------------------------------------------------
C
      DIMENSION R(N),RHO(N),V(N),W1(MNR),W2(MNR)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
C
      AZ2=RC2*AZ
      PI=RC4*ATAN(RC1)
      PI8=RC2*RC4*PI
C
!SMP$ DO SERIAL
      DO 300 I=1,N
      W1(I)=R(I)*RHO(I)*PI8
300   CONTINUE
      CALL PRIM2(N,R,W1,W2)
C
      D=W2(N)
      DO 301 I=1,N
      V(I)=D-W2(I)
301   CONTINUE
C
      DO 302 I=1,N
      W1(I)=W1(I)*R(I)
302   CONTINUE
      CALL PRIM2(N,R,W1,W2)
C
      V(1)=RCZ
      DO 303 I=2,N
      V(I)=V(I)+(-AZ2+W2(I))/R(I)
303   CONTINUE
C
      RETURN
      END
C*******************
CXXX    XCVBH   ****
C*******************
      SUBROUTINE XCVBH(RHO1,RHO2,VXC1,VXC2,EXC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
C---------------------------------------------------------------
C      XC-POTENTIAL AND ENERGY ACCORDING TO VON BARTH AND HEDIN
C         (J. PHYS. C 5 (1972) 1629)
C      WITH PARAMETERS ACCORDING TO HEDIN AND LUNDQVIST
C         (J. PHYS. C 4 (1971) 2064)
C      AND JANAK
C         (SOLID STATE COMMUN. 25 (1978) 53)
C---------------------------------------------------------------
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/,
     &     RC4/4.0D0/
C
      DATA C238/0.2387324146D0/, C916/0.9163305866D0/
C
      DATA CP/0.045D0/,RP/21.0D0/,CF/0.0225D0/,RF/52.916684D0/
C
      IF(RHO1.LE.RCZ) STOP ' ***  RHO1.LE.0  *** '
      IF(RHO2.LE.RCZ) STOP ' ***  RHO2.LE.0  *** '
C
C                                   CONSTANTS
      T13=RC1/RC3
      T43=RC4/RC3
      T213=RC2**T13
C
C                                  VARIABLE RS
      RHO=RHO1+RHO2
      R=(C238/RHO)**T13
C
C                                  VARIABLE DZETA
      S=(RHO1-RHO2)/RHO
      A1PS=RC1+S
      A1MS=RC1-S
C
C                                  FUNCTION F(DZETA)
      DUM=RC2*(T213-RC1)
      FS=(A1PS**T43+A1MS**T43-RC2)/DUM
      DSFS=T43*(A1PS**T13-A1MS**T13)/DUM
C
C                                  EXCHANGE PART
      EPX=-C916/R
      EFX=T213*EPX
C
      EX=EPX+(EFX-EPX)*FS
      RDREX=-EX
      DSEX=(EFX-EPX)*DSFS
C
C                                  CORRELATION PART
      X=R/RP
      CALL AUXVBH(X,G,XDXG)
      EPC=-CP*G
      RDREPC=-CP*XDXG
C
      X=R/RF
      CALL AUXVBH(X,G,XDXG)
      EFC=-CF*G
      RDREFC=-CF*XDXG
C
      EC=EPC+(EFC-EPC)*FS
      RDREC=RDREPC+(RDREFC-RDREPC)*FS
      DSEC=(EFC-EPC)*DSFS
C
C                                  XC QUANTITIES
      EXC=EX+EC
      RDREXC=RDREX+RDREC
      DSEXC=DSEX+DSEC
      DUM=EXC-T13*RDREXC
      VXC1=DUM+A1MS*DSEXC
      VXC2=DUM-A1PS*DSEXC
C
      RETURN
      END
C*******************
CXXX   AUXVBH   ****
C*******************
      SUBROUTINE AUXVBH(X,G,XDXG)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
C
C--------------------------------------------------
C     AUXILIARY FUNCTION FOR VON BARTH - HEDIN XC.
C     FUNCTION G(X) IS DEFINED BY EQ. (5.11) OF
C     THEIR ARTICLE.
C     XDXG MEANS X TIMES X-DERIVATIVE OF G(X).
C--------------------------------------------------
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RCH/0.5D0/,RC3/3.0D0/
C
      DATA XMAX/20.0D0/, RELERR/1.0D-7/
C
       IF(X.LE.XMAX) THEN
      GOL=LOG(RC1+RC1/X)
      BRA=X*(X*(X*GOL-RC1)+RCH)-RC1/RC3
      G=GOL+BRA
      XDXG=RC3*BRA
       ELSE
      T=-RC1/X
      TK=RC1
      K=0
      SUMA=RCZ
      SUMB=RCZ
222   K=K+1
      TK=TK*T
      CLENB=TK/REAL(K+3)
      CLENA=CLENB/REAL(K)
      SUMA=SUMA+CLENA
      SUMB=SUMB+CLENB
      IF(ABS(CLENB).GT.RELERR*ABS(SUMB)) GO TO 222
      G=-RC3*SUMA
      XDXG=RC3*SUMB
       END IF
C
      RETURN
      END
