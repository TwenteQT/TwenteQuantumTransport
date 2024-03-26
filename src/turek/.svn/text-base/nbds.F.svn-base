C*******************
CXXX    BDS     ****
C*******************
      PROGRAM BDS
C
C*******************************************************
C      PROGRAM FOR CALCULATION OF
C       BULK DENSITIES OF STATES AND
C       THE BLOCH SPECTRAL FUNCTIONS
C       BY MEANS OF THE TB-LMTO-CPA METHOD
C*******************************************************
C      VERSION:    JANUARY 1999
C*******************************************************
C   BASIC FEATURES OF THE PROGRAM ARE:
C     - ELECTRONIC STRUCTURE OF ORDERED AND DISORDERED
C       BULK SYSTEMS
C     - CASES: CUBIC, HEXAGONAL, ...
C     - GEOMETRY: QUITE GENERAL
C       (WITH 3D-TRANSLATIONAL SYMMETRY)
C     - NON-MAGNETIC AND SPIN-POLARIZED
C       (WITH COLLINEAR SPIN STRUCTURES ONLY)
C     - SUBSTITUTIONAL CHEMICAL DISORDER
C     - ATOMIC SPHERE APPROXIMATION
C     - SPD AND SPDF CASES
C*******************************************************
C   INPUT:   IR1  (= 1)  -  GEOMETRY OF THE PROBLEM
C            IR2  (= 2)  -  CHEMICAL OCCUPATIONS
C                           OF THE SITES
C            IR3  (= 3)  -  LSDA FILE (POT. PARAMETERS)
C            IR4  (= 4)  -      NOT USED
C            IR5  (= 5)  -  CONTROL PARAMETERS
C   OUTPUT:  IW6  (= 6)  -  CURRENT OUTPUT
C            IW7  (= 7)  -  ORBITAL-RESOLVED DOS
C            IW8  (= 8)  -  ATOM-RESOLVED DOS (SPIN UP)
C            IW9  (= 9)  -  ATOM-RESOLVED DOS (SPIN DOWN)
C            IW17 (=17)  -  ORBITAL-RESOLVED BSF
C            IW18 (=18)  -  SITE-RESOLVED BSF (SPIN UP)
C            IW19 (=19)  -  SITE-RESOLVED BSF (SPIN DOWN)
C*******************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C*******************************************************
C         IMPORTANT DIMENSIONS AND PARAMETERS:
C
C    NB ...... SIZE OF THE BASIS
C              (NUMBER OF ATOMIC SITES IN THE BASIS)
C   MNB ...... PARAMETER: MAX. NB
C    NA ...... NUMBER OF CHEMICALLY DIFFERENT
C              ATOMS FOR THE LSDA PART
C   MNA ...... PARAMETER: MAX. NA
C    NL ...... NUMBER OF VALENCE RADIAL WAVE FUNCTIONS
C              (NL=3 FOR SPD-CASE, NL=4 FOR SPDF-CASE)
C   MNL ...... PARAMETER: MAX. NL
C    NS ...... NUMBER OF SPIN ORIENTATIONS
C              (NS=1 FOR NON-MAGNETIC CASE,
C               NS=2 FOR SPIN-POLARIZED CASE)
C   MNS ...... PARAMETER: MAX. NS
C    NE ...... NUMBER OF ENERGY POINTS
C   MNE ...... PARAMETER: MAX. NE
C    NBZ ..... NUMBER OF K-POINTS IN THE IBZ
C   MNBZ ..... PARAMETER: MAX. NBZ
C    NBZA .... NUMBER OF K-POINTS FOR THE BSF
C   MNBZA .... PARAMETER: MAX. NBZA
C    NLSQ .... DIMENSION OF THE ON-SITE BLOCKS
C              OF GREEN'S FUNCTION (NLSQ=NL**2)
C   MNLSQ .... PARAMETER: MAX. NLSQ
C   MLMAX .... PARAMETER: MAX. ANGULAR QUANTUM NUMBER
C              OF SPHERICAL HARMONICS IN CANONICAL
C              STRUCTURE CONSTANTS ( MLMAX=2*(MNL-1) )
C   MHARM .... PARAMETER: MAX. NUMBER OF SPHERICAL
C              HARMONICS IN CANONICAL STRUCTURE
C              CONSTANTS ( MHARM=(MLMAX+1)**2 )
C    NBLSQ ... DIMENSION OF THE K-DEPENDENT GREEN'S
C              FUNCTION (NBLSQ=NB*NLSQ)
C   MNBLSQ ... PARAMETER: MAX. NBLSQ
C   MNCL ..... PARAMETER: MAX. SIZE OF SMALL CLUSTERS
C              FOR CALCULATION OF TB-STRUCTURE CONSTANTS
C   MTBCL .... PARAMETER: MAX. DIMENSION OF MATRICES
C              FOR CALCULATION OF TB-STRUCTURE CONSTANTS
C              (MTBCL=MNCL*MNLSQ)
C******************************************************
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNBZA=300)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=32)
C
C*******************************************************
C             ALL COMMON BLOCKS  :
C
C-----------------------------------------  GEOMETRY
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/KMESHA/ AKBZA(3,MNBZA),NBZA
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB),JBVA(MNCL,MNB),
     &            JTRA(3,MNCL,MNB),NSCL(MNB)
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/SCKA/ ZSKA(MNBLSQ,MNBLSQ,MNBZA)
C
C-------------------------------------  GREEN'S FUNCTIONS
C
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/GDOS/ ZGDOS(MNL,MNA,MNE,MNS)
      COMMON/GBLO/ ZGBLO(MNL,MNB,MNBZA,MNE,MNS)
C
C-------------------------------------  LSDA
C
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
C
C-------------------------------------  GENERAL
C
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
C
C-------------------------------------  SPECIAL
C
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
C
C******************************************************
C
C                             INITIALIZATION
      CALL PINI
C                             READING OF ALL DATA
      CALL RALL
      flush(IW6)      
C                             STARTING CALCULATIONS
      CALL BEGI
      flush(IW6)      
C                   REAL SPACE TB-STRUCTURE CONSTANTS
      CALL TBR
      flush(IW6)      
C                             GENERATING K-MESH
      CALL GIBZ
      flush(IW6)      
C                             LFT OF TB-CONSTANTS
      CALL TBK
      flush(IW6)      
C                             LFT OF TB-CONSTANTS FOR BSF
      IF (NBZA.GT.0) CALL ATBK
      flush(IW6)      
C
C----------------------------------  ENERGY LOOP
      DO 310 IE=1,NE
C                               POTENTIAL FUNCTIONS
      CALL SETPF(IE)
C                               CPA SELFCONSISTENCY
C                               AND GF FOR DOS
      DO 320 IS=1,NS
      CALL CPAIT(IE,IS)
      
C       write(*,*) ZCPF(3,3,1,IS)
320   CONTINUE
C                               GF FOR BSF
      IF (NBZA.GT.0) THEN
      DO 330 IS=1,NS
      CALL GBSF(IE,IS)
330   CONTINUE
      flush(IW6)      
      END IF
C------------------------------  END OF ENERGY LOOP
310   CONTINUE
C                                DOS EVALUATION
      CALL XDOS
C                                BSF EVALUATION
      IF (NBZA.GT.0) CALL XBSF
C
290   WRITE(IW6,191)
191   FORMAT(/8X,'* * *   END OF BDS   * * * '/)
C
      STOP
      END
C*******************
CXXX    PINI    ****
C*******************
      SUBROUTINE PINI
C
C************************************
C   INITIALIZATION OF THE PROGRAM
C************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
C
C                     INPUT UNITS
      IR1=1
      IR2=2      
      IR4=4
      IR5=5
C                    OUTPUT UNITS
      IW6=6
      IW7=7
      IW8=8
      IW9=9
      IW17=17
      IW18=18
      IW19=19
C
      OPEN(UNIT=IR1,FILE='inpge',FORM='FORMATTED')
      OPEN(UNIT=IR2,FILE='inpch',FORM='FORMATTED')
C      OPEN(UNIT=IR3,FILE='inpld',FORM='FORMATTED')
      OPEN(UNIT=IR5,FILE='inpds',FORM='FORMATTED')
      OPEN(UNIT=IW6,FILE='otdit',FORM='FORMATTED')
      OPEN(UNIT=IW7,FILE='otdor',FORM='FORMATTED')
      OPEN(UNIT=IW8,FILE='otda1',FORM='FORMATTED')
      OPEN(UNIT=IW9,FILE='otda2',FORM='FORMATTED')
      OPEN(UNIT=IW17,FILE='otbor',FORM='FORMATTED')
      OPEN(UNIT=IW18,FILE='otba1',FORM='FORMATTED')
      OPEN(UNIT=IW19,FILE='otba2',FORM='FORMATTED')
      OPEN(UNIT=355,FILE='abk1',FORM='FORMATTED')
      OPEN(UNIT=356,FILE='abk2',FORM='FORMATTED')
C
      RETURN
      END
C*******************
CXXX    RALL    ****
C*******************
      SUBROUTINE RALL
C
C***************************************
C   INPUT OF ALL DATA
C***************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNBZA=300)
      integer,parameter::atomf=999            
C
      DIMENSION OWORK(5),SCX(3),SCKV(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/KMESHA/ AKBZA(3,MNBZA),NBZA
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
      character(len=16)::cwork       
C
100   FORMAT(5A16)
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
C
      WRITE(IW6,111)
111   FORMAT(/'       ****  OUTPUT OF BDS  **** '/
     &     /'     ** READING OF ALL DATA:')
C
C-------------------------  CONTROL DATA - UNIT IR5
C
1001  FORMAT(//5A16)
107   FORMAT(1X,7F10.5)
C
      READ (IR5,100) OWORK
      WRITE(IW6,1001) OWORK
      READ (IR5,100) OWORK
      WRITE(IW6,1001) OWORK
C
      READ (IR5,101) NL,NS
      IF(NL.LT.3.OR.NL.GT.4) GO TO 210
      IF(NL.GT.MNL) GO TO 211
      IF(NS.GT.MNS) GO TO 212
      READ (IR5,101) NSYM,NK,INVE
      READ (IR5,101) NE,NDC
      IF(NE.GT.MNE) GO TO 217
      READ (IR5,104) EBOT,ETOP
      READ (IR5,107) DIFMC,QMIXC
      READ (IR5,101) NBZA
      IF(NBZA.GT.MNBZA) GO TO 218
      IF(NBZA.GT.0) THEN
      READ (IR5,104) (SCKV(I),I=1,3)
      DO 307 IBZA=1,NBZA
      READ (IR5,1041) (AKBZA(I,IBZA),I=1,3),OTXTK(IBZA)
307   CONTINUE
1041  FORMAT(1X,3G15.7,A16)
      END IF
C                                   SCALING
      IF(NBZA.GT.0) THEN
      DO 3081 IBZA=1,NBZA
      DO 3082 I=1,3
      AKBZA(I,IBZA)=SCKV(I)*AKBZA(I,IBZA)
3082  CONTINUE
3081  CONTINUE
      END IF
C
      WRITE(IW6,115) NL,NS
115   FORMAT(/17X,'NL=',I2,17X,'NS=',I2)
      WRITE(IW6,117) NSYM,NK,INVE
117   FORMAT(/'       LOCAL SYMMETRY:    NSYM=',I2/
     &      /'     K-MESH:   NK=',I3,10X,'INVE=',I2)
      WRITE(IW6,118) NE,NDC,EBOT,ETOP
118   FORMAT(/'    ENERGY MESH:    NE=',I4,'       NDC=',I2/
     &    /13X,' EBOT=',F12.5,'          ETOP=',F12.5)
      WRITE(IW6,121) DIFMC,QMIXC
121   FORMAT(/'  CPA ITERATIONS:  DIFMAX=',G12.4,
     &        8X,'MIXING=',F10.5)
      WRITE(IW6,124) NBZA
124   FORMAT(/'     BLOCH SPECTRAL FUNCTIONS:    NBZA=',I5)
      IF(NBZA.GT.0) THEN
      WRITE(IW6,126) (SCKV(I),I=1,3)
126   FORMAT(/'  SCALING FACTORS: ',3G15.7)
      WRITE(IW6,128)
128   FORMAT(/'          SELECTED K-VECTORS:')
      DO 309 IBZA=1,NBZA
      WRITE(IW6,1042) IBZA,OTXTK(IBZA),(AKBZA(I,IBZA),I=1,3)
309   CONTINUE
1042  FORMAT(2X,I2,2X,A16,2X,3G15.7)
      END IF
C
C-------------------------  GEOMETRY - UNIT IR1
C
      READ (IR1,100) OWORK
      WRITE(IW6,1001) OWORK
C
      READ (IR1,101) NB,NMTR
      IF(NB.GT.MNB) GO TO 222
      READ (IR1,104) CUTRAT
      READ (IR1,104) (SCX(I),I=1,3)
      DO 310 J=1,3
      READ (IR1,104) (VBR(I,J),I=1,3)
310   CONTINUE
      DO 311 IB=1,NB
      READ (IR1,104) (POS(I,IB),I=1,3)
311   CONTINUE
C
C                          SCALING OF COORDINATES
      DO 3151 J=1,3
      DO 3152 I=1,3
      VBR(I,J)=SCX(I)*VBR(I,J)
3152  CONTINUE
3151  CONTINUE
      DO 3191 IB=1,NB
      DO 3192 I=1,3
      POS(I,IB)=SCX(I)*POS(I,IB)
3192  CONTINUE
3191  CONTINUE
C
      WRITE(IW6,130) NB,NMTR,CUTRAT
130   FORMAT(/18X,'NB=',I2//'   MAX. ',
     & 'COEFFICIENT OF TRANSL. VECTORS:    NMTR=',I2/
     & /'   CUT-OFF DISTANCE / WS-RADIUS:  CUTRAT=',F10.5)
      WRITE(IW6,131) (SCX(I),I=1,3)
131   FORMAT(/'  SCALING FACTORS: ',3G15.7)
      WRITE(IW6,132)((VBR(I,J),I=1,3),J=1,3)
132   FORMAT(/8X,'TRANSLATION VECTORS IN REAL SPACE:'/
     &        '    1. VECTOR = ',3G15.7/
     &        '    2. VECTOR = ',3G15.7/
     &        '    3. VECTOR = ',3G15.7)
1044  FORMAT(3X,'IB=',I2,5X,3G15.7)
      WRITE(IW6,135)
135   FORMAT(/'       SITES OF THE ATOMIC BASIS: ')
      DO 320 IB=1,NB
      WRITE(IW6,1044) IB,(POS(I,IB),I=1,3)
320   CONTINUE
C
C------------------------- CHEMICAL OCCUPATION - UNIT IR2
C
      READ (IR2,100) OWORK
      WRITE(IW6,1001) OWORK
C
108   FORMAT(A16)
      IA=0
      DO 340 IB=1,NB
      READ (IR2,101) NC(IB)
      DO 341 IC=1,NC(IB)
      IA=IA+1
      IF(IA.GT.MNA) GO TO 241
      READ (IR2,108) OTXTA(IA)
      READ (IR2,104) CON(IA),VALZ(IA)
341   CONTINUE
340   CONTINUE
      NA=IA
C
      IA=0
      DO 344 IB=1,NB
      WRITE(IW6,141) IB,NC(IB)
141   FORMAT(/'    SITE:  IB=',I2,
     &     '      OCCUPIED BY ',I2,' ATOMIC TYPE(S): ')
      DO 345 IC=1,NC(IB)
      IA=IA+1
      WRITE(IW6,142) IA,OTXTA(IA),CON(IA),VALZ(IA)
142   FORMAT(' IA=',I2,' - LABEL,CONC.,VALENCY: ',
     &                       A16,F10.5,F10.3)
345   CONTINUE
344   CONTINUE
      WRITE(IW6,143) NA
143   FORMAT(/4X,' TOTAL NUMBER OF DIFFERENT ',
     &    'ATOMIC TYPES:  NA=',I3)
C
C ---------------------------- LSDA FILE - UNIT IR3
C

      do IA=1,NA
         read(OTXTA(IA),*) cwork
         WRITE(IW6,*)' Read defs for:', CWORK
         open (unit=atomf, file='atoms/'//trim(cwork)//'_B',
     &        action='read')                   
         
         do is=1,3
            read(atomf,*) cwork
         end do
         READ (atomf,104) AZ(IA),WS(IA),WSAV(IA)
      
         do is=1,4
            read(atomf,*) cwork
         end do         

         do IS=1,NS
            READ (atomf,101) NRAD1
            READ (atomf,104) (skipPOT,I=1,NRAD1)
         end do         

         read(atomf,*) cwork
         do IS=1,NS
            READ (atomf,104) cwork
            READ (atomf,104) (ENY(IL,IS,IA),IL=1,NL)
C             write (*,104) (ENY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPC(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPD(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPQ(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPP(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (DNY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (FINY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (FINYD(IL,IS,IA),IL=1,NL)
         end do
         close(atomf)
      end do
    
      WRITE(IW6,150)
150   FORMAT(/6X,'  LABEL,         IA,',
     &        '      Z,         WS,         WSAV: ')
      DO 355 IA=1,NA
      WRITE(IW6,151) OTXTA(IA),IA,AZ(IA),WS(IA),WSAV(IA)
355   CONTINUE
151   FORMAT(4X,A16,I5,3F12.5)
C
      WRITE(IW6,195)
195   FORMAT(/'     **  END OF READING ')
      RETURN
C
C-------------------------------- ERROR MESSAGES
C
210   WRITE(IW6,910)
910   FORMAT(/' **** INPUT ERROR:  NL MUST BE 3 OR 4 ')
      STOP
211   WRITE(IW6,911)
911   FORMAT(/' **** INPUT ERROR:  NL GREATER THAN MNL ')
      STOP
212   WRITE(IW6,912)
912   FORMAT(/' **** INPUT ERROR:  NS GREATER THAN MNS ')
      STOP
217   WRITE(IW6,917)
917   FORMAT(/' **** INPUT ERROR:  NE GREATER THAN MNE ')
      STOP
218   WRITE(IW6,918)
918   FORMAT(/' **** INPUT ERROR: NBZA GREATER THAN MNBZA ')
      STOP
222   WRITE(IW6,922)
922   FORMAT(/' **** INPUT ERROR: NB GREATER THAN MNB  ')
      STOP
241   WRITE(IW6,941)
941   FORMAT(/' **** INPUT ERROR: IA EXCEEDED MNA  ')
      STOP
      END
C*******************
CXXX    BEGI    ****
C*******************
      SUBROUTINE BEGI
C
C************************************
C   TESTS OF INPUT DATA AND
C   STARTING AUXILIARY CALCULATIONS
C************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
C
      DIMENSION CR23(3),CR31(3),CR12(3),QSC3(3),QSC4(4)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/
      DATA RC4/4.0D0/,TOL/1.0D-6/
C
      DATA QSC3/0.348485D0,0.053030D0,0.010714D0/
      DATA QSC4/0.385057D0,0.073209D0,0.022481D0,0.006069D0/
C
      PI=RC4*ATAN(RC1)
C
C----------------------------------- SCREENING CONSTANTS
      IF(NL.EQ.3) THEN
      DO 303 IL=1,3
      QSCR(IL)=QSC3(IL)
303   CONTINUE
      END IF
      IF(NL.EQ.4) THEN
      DO 304 IL=1,4
      QSCR(IL)=QSC4(IL)
304   CONTINUE
      END IF
      WRITE(IW6,115) (QSCR(IL),IL=1,NL)
115   FORMAT(/3X,'SCREENING CONSTANTS: '/
     &     5X,4G15.7)
C
C------------------------------- MAX. NUMBER OF CPA ITERATIONS
      NMITC= 150
      WRITE(IW6,118) NMITC
118   FORMAT(/3X,' MAX. NUMBER OF CPA ITERATIONS: ',I5)
C
C---------------------------------- NCALLH FOR SPH. HARMONICS
      NCALLH=0
C
C---------------------------------- GAUNT FACTORS
      CALL GAUNTY(NL)
C
C---------------------------------- COMPLEX ENERGIES
C
      EH=(ETOP-EBOT)/DBLE(NE-1)   
C       EH=1.0d-7   
C       write(*,*) ETOP-EBOT,EH
      EI= DBLE(NDC)*EH
      DO 310 IE=1,NE
      ER=EBOT+DBLE(IE-1)*EH
      ZCN(IE)=DCMPLX(ER,EI)
310   CONTINUE
      WRITE(IW6,121) EH,EI
121   FORMAT(/9X,'ENERGY MESH:   STEP= ',G15.7/
     &        9X,'     IMAGINARY PART= ',G15.7)
C
C----------------------------------- BASIS IN RECIP. SPACE
C
      CR23(1)=VBR(2,2)*VBR(3,3)-VBR(3,2)*VBR(2,3)
      CR23(2)=VBR(3,2)*VBR(1,3)-VBR(1,2)*VBR(3,3)
      CR23(3)=VBR(1,2)*VBR(2,3)-VBR(2,2)*VBR(1,3)
C
      CR31(1)=VBR(2,3)*VBR(3,1)-VBR(3,3)*VBR(2,1)
      CR31(2)=VBR(3,3)*VBR(1,1)-VBR(1,3)*VBR(3,1)
      CR31(3)=VBR(1,3)*VBR(2,1)-VBR(2,3)*VBR(1,1)
C
      CR12(1)=VBR(2,1)*VBR(3,2)-VBR(3,1)*VBR(2,2)
      CR12(2)=VBR(3,1)*VBR(1,2)-VBR(1,1)*VBR(3,2)
      CR12(3)=VBR(1,1)*VBR(2,2)-VBR(2,1)*VBR(1,2)
C
      DET=RCZ
      DO 340 I=1,3
      DET=DET+VBR(I,1)*CR23(I)
340   CONTINUE
      ROMEGA=ABS(DET)
      GOMEGA=(RC2*PI)**3/ROMEGA
      CNORM=RC2*PI/DET
C
      DO 341 I=1,3
      VBG(I,1)=CNORM*CR23(I)
341   CONTINUE
      DO 342 I=1,3
      VBG(I,2)=CNORM*CR31(I)
342   CONTINUE
      DO 343 I=1,3
      VBG(I,3)=CNORM*CR12(I)
343   CONTINUE
C
      WRITE(IW6,124) ROMEGA,GOMEGA
124   FORMAT(/' VOLUME OF PRIMITIVE CELL: REAL= ',G15.7/
     &        '                     RECIPROCAL= ',G15.7)
      WRITE(IW6,125)((VBG(I,J),I=1,3),J=1,3)
125   FORMAT(/'              RECIPROCAL BASIS: '/
     &      '    1. VECTOR = ',3G15.7/
     &      '    2. VECTOR = ',3G15.7/
     &      '    3. VECTOR = ',3G15.7)
C
C---------------------------- CHECK OF CONSISTENCY
C                             OF CONCENTRATIONS
      IA=0
      DO 325 IB=1,NB
      SUM=RCZ
      DO 326 IC=1,NC(IB)
      IA=IA+1
      IF(CON(IA).LT.RCZ.OR.CON(IA).GT.RC1) GO TO 280
      SUM=SUM+CON(IA)
326   CONTINUE
      IF(ABS(SUM-RC1).GT.TOL) GO TO 280
325   CONTINUE
C
C--------------------------- RANDOM/NON-RANDOM CASE
      IRAND=0
      DO 328 IA=1,NA
      PROD=CON(IA)*(RC1-CON(IA))
      IF(PROD.GT.RCZ) IRAND=1
328   CONTINUE
      IF(IRAND.EQ.0) WRITE(IW6,126)
126   FORMAT(/5X,'  --- SYSTEM WITHOUT DISORDER --- ')
C
C----------------------------------- AVERAGE WS-RADII:
C
C                                 DIMENSIONLESS
      AVOL=ROMEGA/dble(NB)
      AWS=((RC3*AVOL)/(RC4*PI))**(RC1/RC3)
      WRITE(IW6,127) AWS
127   FORMAT(/'    DIMENSIONLESS AVER. WS-RADIUS= ',G15.7)
C
C                                IN ATOMIC UNITS
      SUM=RCZ
      DO 330 IA=1,NA
      SUM=SUM+CON(IA)*WS(IA)**3
330   CONTINUE
      TWS=(SUM/dble(NB))**(RC1/RC3)
      WRITE(IW6,128) TWS
128   FORMAT(/'  AVER. WS-RADIUS IN ATOMIC UNITS= ',G15.7)
C
C---------------------------------- CHECK OF CONSISTENCY
C                                    OF NSYM, INVE AND
C                                   TRANSLATION VECTORS
C
      IF(NSYM.LT.0.OR.NSYM.GT.5) GO TO 295
      IF(NSYM.EQ.0) GO TO 222
      IF(NSYM.EQ.1) GO TO 201
      IF(NSYM.EQ.2) GO TO 202
      IF(NSYM.EQ.3) GO TO 203
      IF(NSYM.EQ.4) GO TO 204
      IF(NSYM.EQ.5) GO TO 205
C
C--------  NSYM = 1 - FULL CUBIC SYMMETRY
C
201   IF(INVE.LT.1.OR.INVE.GT.3) GO TO 296
C
C                INVE = 1: SIMPLE CUBIC LATTICE
C                INVE = 2: FACE-CENTRED CUBIC LATTICE
C                INVE = 3: BODY-CENTRED CUBIC LATTICE
C
C--------  NSYM = 3  -  TETRAHEDRAL SYMMETRY OF ORDERED
C                       DO3 (AND SIMILAR) STRUCTURES
203   CONTINUE
C
C                   INVE - NOT ACTIVE
C
C     ALL 3 VECTORS MUST HAVE THE SAME MAGNITUDE
C
      SUM1=VBR(1,1)**2+VBR(2,1)**2+VBR(3,1)**2
      SUM2=VBR(1,2)**2+VBR(2,2)**2+VBR(3,2)**2
      SUM3=VBR(1,3)**2+VBR(2,3)**2+VBR(3,3)**2
      DIF=SUM2/SUM1-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=SUM3/SUM1-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      GO TO 222
C
C--------  NSYM = 2 - HEXAGONAL SYMMETRY
C                SYMMETRY GROUP GIVEN BY 3-FOLD ROTATION
C                AXIS Z AND 2 MIRROR PLANES X-Y AND Y-Z
C
C          INVE - NOT ACTIVE.
C     THE TRANSLATION VECTORS MUST HAVE THE FORM:
C                          A1 = (A,0,0)
C                          A2 = (A/2)*(1,S,0)*U
C                          A3 = (0,0,C)
C       WHERE A AND C ARE POSITIVE,
C       S=SQRT(3) OR S=-SQRT(3), AND U=1 OR U=-1.
C
202   IF (ABS(VBR(2,1)).GT.RCZ) GO TO 297
      IF (ABS(VBR(3,1)).GT.RCZ) GO TO 297
      IF (VBR(1,1).LE.RCZ) GO TO 297
      IF (ABS(VBR(3,2)).GT.RCZ) GO TO 297
      DIF=RC2*ABS(VBR(1,2))/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=(VBR(1,2)**2+VBR(2,2)**2)/VBR(1,1)**2-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      IF (ABS(VBR(1,3)).GT.RCZ) GO TO 297
      IF (ABS(VBR(2,3)).GT.RCZ) GO TO 297
      IF (VBR(3,3).LE.RCZ) GO TO 297
      GO TO 222
C
C--------  NSYM = 4 - TETRAGONAL SYMMETRY
C
204   IF(INVE.LT.1.OR.INVE.GT.2) GO TO 296
      IF(INVE.EQ.1) GO TO 2041
      IF(INVE.EQ.2) GO TO 2042
C
C            INVE = 1 :  THE TRANSLATION VECTORS
C                        MUST HAVE THE FORM:
C                          A1 = (A,0,0)
C                          A2 = (0,A,0)
C                          A3 = (0,0,C)
C                WHERE A AND C ARE POSITIVE.
C
2041  IF (ABS(VBR(2,1)).GT.RCZ) GO TO 297
      IF (ABS(VBR(3,1)).GT.RCZ) GO TO 297
      IF (VBR(1,1).LE.RCZ) GO TO 297
      IF (ABS(VBR(1,2)).GT.RCZ) GO TO 297
      IF (ABS(VBR(3,2)).GT.RCZ) GO TO 297
      IF (VBR(2,2).LE.RCZ) GO TO 297
      DIF=VBR(1,1)/VBR(2,2)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      IF (ABS(VBR(1,3)).GT.RCZ) GO TO 297
      IF (ABS(VBR(2,3)).GT.RCZ) GO TO 297
      IF (VBR(3,3).LE.RCZ) GO TO 297
      GO TO 222
C
C            INVE = 2 :  THE TRANSLATION VECTORS
C                        MUST HAVE THE FORM:
C                          A1 =  (P,P,0)
C                          A2 = (-P,P,0)*U
C                          A3 =  (0,0,C)
C                  WHERE P AND C ARE POSITIVE
C                           AND U=1 OR U=-1.
C
2042  IF (VBR(1,1).LE.RCZ) GO TO 297
      DIF=VBR(2,1)/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      IF (ABS(VBR(2,2)).LE.RCZ) GO TO 297
      DIF=VBR(1,2)/VBR(2,2)+RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=ABS(VBR(2,2))/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      IF (ABS(VBR(1,3)).GT.RCZ) GO TO 297
      IF (ABS(VBR(2,3)).GT.RCZ) GO TO 297
      IF (VBR(3,3).LE.RCZ) GO TO 297
      GO TO 222
C
C--------  NSYM = 5 - ORTHORHOMBIC SYMMETRY
C
C            INVE - NOT ACTIVE.
C     THE TRANSLATION VECTORS MUST HAVE THE FORM:
C                          A1 = (A,0,0)
C                          A2 = (0,B,0)
C                          A3 = (0,0,C)
C                WHERE A, B, C ARE POSITIVE.
C
205   IF (ABS(VBR(2,1)).GT.RCZ) GO TO 297
      IF (ABS(VBR(3,1)).GT.RCZ) GO TO 297
      IF (VBR(1,1).LE.RCZ) GO TO 297
      IF (ABS(VBR(1,2)).GT.RCZ) GO TO 297
      IF (ABS(VBR(3,2)).GT.RCZ) GO TO 297
      IF (VBR(2,2).LE.RCZ) GO TO 297
      IF (ABS(VBR(1,3)).GT.RCZ) GO TO 297
      IF (ABS(VBR(2,3)).GT.RCZ) GO TO 297
      IF (VBR(3,3).LE.RCZ) GO TO 297
      GO TO 222
C
222   CONTINUE
C
      RETURN
C
C--------------------------------- ERROR MESSAGES
C
280   WRITE(IW6,180)
180   FORMAT(/' **** INPUT ERROR: ',
     &    ' INCONSISTENCY IN CONCENTRATIONS ')
      STOP
295   WRITE(IW6,195)
195   FORMAT(/' **** INPUT ERROR: NSYM MUST BE ',
     &                            '0,1,2,3,4,5')
      STOP
296   WRITE(IW6,196)
196   FORMAT(/' **** INPUT ERROR:  NSYM AND INVE ',
     &  'ARE INCONSISTENT')
      STOP
297   WRITE(IW6,197)
197   FORMAT(/' **** INPUT ERROR:  NSYM, INVE '/3X,
     &  '  AND TRANSLATION VECTORS ARE INCONSISTENT')
      STOP
C
      END
C*******************
CXXX    TBCL    ****
C*******************
      SUBROUTINE TBCL(NL,NCL,POL,WSA)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=32)
      PARAMETER(MTBCL=MNCL*MNLSQ)
C
C--------------------------------------------------------
C     TB-STRUCTURE CONSTANTS FOR A SMALL CLUSTER
C     CALCULATED BY MATRIX INVERSION IN REAL SPACE
C--------------------------------------------------------
C   INPUT:  NL - NUMBER OF PARTIAL WAVES,
C           NCL - SIZE OF THE SMALL CLUSTER,
C           POL(K,ICL) (K=1,2,3, ICL=1,2,...,NCL) -
C            - COORDINATES OF THE ATOMS OF THE SMALL CLUSTER,
C           ICL=1 CORRESPONDS TO THE CENTRAL ATOM OF THE CLUSTER,
C           WSA - AVERAGE W.S. RADIUS.
C  OUTPUT:  SSC(I,J,ICL) (I,J - ORBITAL INDICES,
C                        ICL=1,2,...,NCL) -
C           - TB-STRUCTURE CONSTANTS S(I,ICL,J,1) CONNECTING
C           THE ICL-TH ATOM WITH THE FIRST (CENTRAL) ATOM
C--------------------------------------------------------
C   REMARK: THE GAUNT FACTORS MUST BE CALCULATED BEFORE !
C           NCALLH MUST BE SET ZERO BEFORE !
C---------------------------------------------------------
C
      DIMENSION POL(3,MNCL),RI(3),RJ(3),QR(MNLSQ),
     &          AA(MTBCL,MTBCL),BB(MTBCL,MNLSQ),WW(MTBCL),
     &          INDW(MNLSQ,MNCL),SC(MNLSQ,MNLSQ)
C
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
#if defined(DO_SC_SYM_TUREK)
       Integer, Pointer :: symfield (:)
       Integer :: NSYMs, srows
#endif
      
C       do i=1,NCL
C           write(*,'(4g11.4)') pol(:,i)
C       end do
C       write(*,*) '-----cls----', NCL
C
      DATA RC1/1.0D0/
C
      IF(NCL.GT.MNCL) STOP ' TBCL-ERROR:  NCL.GT.MNCL ! '
C
      NLSQ=NL**2
C
      DO 300 IL=1,NL
      ISTA=(IL-1)**2+1
      IFIN=IL**2
!SMP$ DO SERIAL      
      DO 301 I=ISTA,IFIN
      QR(I)=RC1/QSCR(IL)
301   CONTINUE
300   CONTINUE
C
      II=0
      DO 304 ICL=1,NCL
      DO 305 I=1,NLSQ
      II=II+1
      INDW(I,ICL)=II
305   CONTINUE
304   CONTINUE
      NN=II
#if defined(DO_SC_SYM_TUREK)
       srows = nlsq
       NSYMs = srows * (nn-srows)
       Allocate (symfield(2*NSYMs))
#endif

C
C                   FORMATION OF THE MATRIX (QR - S) (=AA)
C                     AND OF THE RIGHT-HAND SIDES (=BB)
C
      DO 310 JCL=1,NCL
      DO 3101 K=1,3
      RJ(K)=POL(K,JCL)
3101  CONTINUE
      DO 312 ICL=JCL,NCL
      DO 3121 K=1,3
      RI(K)=POL(K,ICL)
3121  CONTINUE
      CALL CANSC(NL,WSA,RI,RJ,SC)
      DO 314 J=1,NLSQ
      JJ=INDW(J,JCL)
      DO 315 I=1,NLSQ
      II=INDW(I,ICL)
      AA(II,JJ)=-SC(I,J)
      IF(JCL.EQ.1) BB(II,J)=SC(I,J)
315   CONTINUE
314   CONTINUE
312   CONTINUE
310   CONTINUE
#if defined(DO_SC_SYM_TUREK)
      Call getscsym (NSYMs, aa((srows+1) :nn, 1:srows), symfield, 
     &1.0d-12)
#endif

      DO 318 ICL=1,NCL
      DO 319 I=1,NLSQ
      II=INDW(I,ICL)
      AA(II,II)=QR(I)
319   CONTINUE
318   CONTINUE
C
C                      SOLUTION OF THE LINEAR EQUATIONS
C
      CALL SOPO(AA,MTBCL,NN,BB,MTBCL,NLSQ,WW)
C
C                      THE TB-STRUCTURE CONSTANTS
C


      DO J=1,NLSQ
      DO ICL=1,NCL
      DO I=1,NLSQ
      II=INDW(I,ICL)
      BB(II,J)=QR(I)*BB(II,J)
      end do
      end do
      end do

#if defined(DO_SC_SYM_TUREK)
      Call setscsym (NSYMs, BB((srows+1) :nn, 1:srows), symfield)
      Call symscons (srows*srows, BB(1:srows, 1:srows), 1.0d-5)
      Where (Abs(BB) < 1.0d-14)
         BB = 0.0d0
      End Where
      Deallocate (symfield)
#endif
      DO J=1,NLSQ
      DO ICL=1,NCL
      DO I=1,NLSQ
      II=INDW(I,ICL)
      SSC(I,J,ICL)=BB(II,J)
      end do
      end do
      end do



C
      RETURN
#if defined(DO_SC_SYM_TUREK)

      contains
      Subroutine getscsym (N, SC, f, crit)
      Implicit None
      Integer, Intent (In) :: N
      Real (kind=8), Intent (In) :: crit
      Real (kind=8), Intent (Inout) :: SC (N)
      Integer, Intent (Out) :: f (2*N)
C  Local
      Integer, Pointer :: mask (:)
      Integer :: i, j, mptr, num
      Real (kind=8) :: base
      Allocate (mask(N))
      mask (1:N) = 0
      f (1) = 0
      mptr = 1
      Do i = 1, N, 1
        base = SC (i)
C         If (Abs(base) > crit .And. mask(i) == 0) Then
        If (mask(i) == 0) Then
           num = 1
           Do j = i + 1, N, 1
C              If ((mask(j) == 0) .And. (Abs(sc(j)) > crit)) Then
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
      Real (kind=8), Intent (Inout) :: SC (N)
      Integer, Intent (In) :: f (2*N)
C      Local
      Integer :: i, mptr, cntr
      Real (kind=8) :: base
      mptr = 1
      Do While (f(mptr) /=  0)
        base = 0.0d0
        cntr = f (mptr)
C           write(*,*)
        Do i = mptr + 1, mptr + cntr, 1
           base = base + dble(sign(1, f(i)), kind=8) * SC (Abs(f(i)))
C               write(*,*) sc(abs(f(i)))
        End Do
        base = base / dble(cntr, kind=8)
        Do i = mptr + 1, mptr + cntr, 1
           SC (Abs(f(i))) = dble(sign(1, f(i)), kind=8) * base
        End Do
        mptr = mptr + cntr + 1
      End Do
      End Subroutine setscsym

      Subroutine symscons (N, SC, crit)
      Implicit None
      Integer, Intent (In) :: N
      Real (kind=8), Intent (Inout) :: SC (N)
      Real (kind=8) :: crit
C Locals
      Integer :: i, j, num
      Integer, Pointer :: mask (:), refs (:)
      Real (kind=8) :: base, DEV, abase, aref, acc, ref
      Allocate (mask(N))
      Allocate (refs(N))
      mask (1:N) = 0
C       write(*,*)crit
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
C                   write(*,*)'dev=',dev
                 If (DEV <= crit) Then
                    num = num + 1
                    refs (num) = j * sign (1.0d0, ref)
C                      write(*,*) base, ref
                    acc = acc + aref
                    mask (j) = 1
                 End If
              End If
           End Do

           If (num > 1) Then
              acc = acc / dble(num, kind=8)
              Do j = 1, num, 1
              SC(Abs(refs(j))) = dble(sign(1,refs(j)), kind=8)*acc
              End Do
           End If
        End If
      End Do
      Deallocate (mask, refs)
      End Subroutine symscons
#endif
      
      END
C*******************
CXXX   CANSC    ****
C*******************
      SUBROUTINE CANSC(NL,WSA,R1,R2,SC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
C
C---------------------------------------------------------------
C   CANONICAL STRUCTURE CONSTANT CONNECTING 2 POINTS R1 AND R2
C       NL - NUMBER OF PARTIAL WAVES
C       WSA - AVERAGE W.S. RADIUS
C       R1(I),R2(I) - THE COORDINATES OF THE 2 POINTS
C       SC(I1,I2) - THE STRUCTURE CONSTANT MATRIX
C---------------------------------------------------------------
C   REMARK: THE GAUNT FACTORS MUST BE CALCULATED BEFORE !
C           NCALLH MUST BE SET ZERO BEFORE !
C---------------------------------------------------------------
C
      DIMENSION R1(3),R2(3),SC(MNLSQ,MNLSQ)
      DIMENSION DFAC(0:MLMAX),AM1LP1(0:MNL)
      DIMENSION RBWLP1(0:MLMAX),PREF(0:MNL,0:MNL)
C
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
      DATA RMIN/0.01D0/
C
      DO 302 I2=1,MNLSQ
      DO 301 I1=1,MNLSQ
      SC(I1,I2)=RCZ
301   CONTINUE
302   CONTINUE
C
      TX=R2(1)-R1(1)
      TY=R2(2)-R1(2)
      TZ=R2(3)-R1(3)
      R= SQRT(TX**2+TY**2+TZ**2)
      RBW=R/WSA
      IF(RBW.LT.RMIN) RETURN
C
      LMAX=NL-1
      LMAXH=2*LMAX
C
      DFAC(0)= dble(1)
      DO 310 L=1,LMAXH
      TWOLM1= dble(2*L-1)
      DFAC(L)=TWOLM1*DFAC(L-1)
310   CONTINUE
C
      ASIG=RC1
      DO 311 L=0,LMAX
      ASIG=-ASIG
      AM1LP1(L)=ASIG
311   CONTINUE
C
      PI=RC4*ATAN(RC1)
      PI8=RC2*RC4*PI
C
      RBWLP1(0)=RBW
      DO 313 L=1,LMAXH
      RBWLP1(L)=RBWLP1(L-1)*RBW
313   CONTINUE
C
      DO 315 L2=0,LMAX
      DO 316 L1=0,LMAX
      L=L1+L2
      CIT=AM1LP1(L2)*PI8*DFAC(L)
      AJM=DFAC(L1)*DFAC(L2)*RBWLP1(L)
      PREF(L1,L2)=CIT/AJM
316   CONTINUE
315   CONTINUE
C
      DX=TX/R
      DY=TY/R
      DZ=TZ/R
C
      CALL HARM(LMAXH,DX,DY,DZ)
C
      I2=0
      DO 3228 L2=0,LMAX
      DO 3229 M2=-L2,L2
      I2=I2+1
C
      I1=0
      DO 3218 L1=0,LMAX
      DO 3219 M1=-L1,L1
      I1=I1+1
C
      L=L1+L2
      ISTA=L**2+1
      IFIN=(L+1)**2
      SUM=RCZ
      DO 320 I=ISTA,IFIN
      SUM=SUM+GFRH(I,I1,I2)*YPS(I)
320   CONTINUE
      SC(I1,I2)=PREF(L1,L2)*SUM
C
3219  CONTINUE
3218  CONTINUE
3229  CONTINUE
3228  CONTINUE
C
      RETURN
      END
C*******************
CXXX    HARM    ****
C*******************
      SUBROUTINE HARM(LMAX,D1,D2,D3)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNL=3)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
C
C---------------------------------------------------------
C   REAL SPHERICAL HARMONICS NORMALIZED TO UNITY
C     Y(L,M;D1,D2,D3) = C(L,|M|)*P(L,|M|;THETA)*T(M;PHI)
C         WHERE: COS(THETA)=D3,
C                SIN(THETA)*EXP(I*PHI)=D1+I*D2,
C                C(L,|M|) IS POSITIVE NORMALIZING FACTOR,
C                P(L,|M|;THETA) IS LEGENDRE POLYNOMIAL,
C                T(M;PHI) = COS(M*PHI)   FOR M.GE.0,
C                T(M;PHI) = SIN(|M|*PHI) FOR M.LT.0.
C---------------------------------------------------------
C          INPUT:
C    NCALLH - MUST BE ZERO BEFORE THE FIRST CALL !!
C    LMAX - MAXIMAL ANGULAR NUMBER (2.LE.LMAX.LE.MLMAX)
C    D1,D2,D3 - DIRECTION COSINES
C         OUTPUT:
C    YPS(.) - REAL SPHERICAL HARMONICS
C             THE HARMONICS Y(L,M) ARE STORED IN ARRAY YPS(I)
C             IN THE SEQUENCE ((Y(L,M),M=-L,L),L=0,LMAX)
C---------------------------------------------------------
C  REMARK:  IF (D1,D2,D3) IS NOT A UNIT VECTOR,
C           THEN Y(L,M) ARE EQUAL TO THE CORRESPONDING
C           HOMOGENEOUS HARMONIC POLYNOMIALS:
C      Y(L,M;D1,D2,D3) = R**L * Y(L,M;E1,E2,E3) ,
C            WHERE (E1,E2,E3) IS A UNIT VECTOR AND
C            (D1,D2,D3) = R * (E1,E2,E3) .
C----------------------------------------------------------
C
      DIMENSION U(0:MLMAX,-MLMAX:MLMAX)
C
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
C
      IF(LMAX.LT.2) STOP ' ***  HARM: LMAX.LT.2 !'
      IF(LMAX.GT.MLMAX) STOP ' ***  HARM: LMAX.GT.MLMAX !'
C
      IF(NCALLH.GT.0) GO TO 222
C
C                   CALCULATION OF NORMAL. CONSTANTS
C
      NCALLH=1
      PI=RC4*ATAN(RC1)
      CNH(0,0)=RC1/(RC4*PI)
      DO 301 L=1,MLMAX
      TWOLP1=dble(2*L+1)
      CNH(L,0)=TWOLP1/(RC4*PI)
      DO 302 M=1,L
      PROD=RC1
      DO 303 I=L-M+1,L+M
      PROD=PROD*dble(I)
303   CONTINUE
      CNH(L,M)=TWOLP1/(RC2*PI*PROD)
302   CONTINUE
301   CONTINUE
      DO 304 L=0,MLMAX
      DO 305 M=0,L
      CNH(L,M)= SQRT(CNH(L,M))
305   CONTINUE
304   CONTINUE
C
C                       UNNORMALIZED HARMONICS
C
222   CONTINUE
      DSQ=D1**2+D2**2+D3**2
      ZP=DCMPLX(D1,D2)
      ZQ=DCMPLX(RC1,RCZ)
      U(0,0)=RC1
      DO 310 L=1,LMAX
      TWOLM1=dble(2*L-1)
      ZQ=TWOLM1*ZP*ZQ
      U(L,-L)=DIMAG(ZQ)
      U(L,L) =DREAL(ZQ)
310   CONTINUE
      DO 312 M=-LMAX+1,LMAX-1
      L= ABS(M)
      TWOLP1=dble(2*L+1)
      U(L+1,M)=TWOLP1*D3*U(L,M)
312   CONTINUE
      DO 314 M=-LMAX+2,LMAX-2
      MA= ABS(M)
      DO 315 L=MA+1,LMAX-1
      TWOLP1=dble(2*L+1)
      FLM1=dble(L+MA)
      FLP1=dble(L-MA+1)
      U(L+1,M)=(TWOLP1*D3*U(L,M)-FLM1*DSQ*U(L-1,M))/FLP1
315   CONTINUE
314   CONTINUE
C
C                              NORMALIZED HARMONICS
C
      I=0
      DO 320 L=0,LMAX
      DO 321 M=-L,L
      I=I+1
      MA= ABS(M)
      YPS(I)=CNH(L,MA)*U(L,M)
321   CONTINUE
320   CONTINUE
C
      RETURN
      END
C*******************
CXXX   GACOR    ****
C*******************
      REAL*8 FUNCTION GACOR(KL1,KM1,KL2,KM2,KL3,KM3)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------------
C     GAUNT COEFFICIENT FOR REAL SPHERICAL HARMONICS
C----------------------------------------------------
C     ARGUMENTS MUST SATISFY NATURAL INEQUALITIES
C   KL1.GE.ABS(KM1), KL2.GE.ABS(KM2), KL3.GE.ABS(KM3)
C----------------------------------------------------
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
C
      GACOR=RCZ
C
      IG2=KL1+KL2+KL3
      IF(MOD(IG2,2).NE.0) RETURN
      IG=IG2/2
      IF(MAX(KL1,KL2,KL3).GT.IG) RETURN
C
      MTEST=(2*KM1+1)*(2*KM2+1)*(2*KM3+1)
      IF(MTEST.LT.0) RETURN
      KMA1=ABS(KM1)
      KMA2=ABS(KM2)
      KMA3=ABS(KM3)
      MAMAX=MAX(KMA1,KMA2,KMA3)
      IF(2*MAMAX.NE.KMA1+KMA2+KMA3) RETURN
C
C                                        PHI-INTEGRAL
      PI=RC4*ATAN(RC1)
      MAMIN=MIN(KMA1,KMA2,KMA3)
      IF(MAMIN.EQ.0) THEN
       DELTA=RC1/SQRT(RC2*PI)
      ELSE
       DELTA=RC1/(RC2*SQRT(PI))
       IF(KM1+KM2+KM3.EQ.0) DELTA=-DELTA
      END IF
C
C                                      THETA-INTEGRAL
C   SEE E.U. CONDON, G.H. SHORTLEY:
C      THE THEORY OF ATOMIC SPECTRA (CAMBRIDGE 1957)
C      PAGE 176, EQ.(11)
C
      IF(KMA1.EQ.MAMAX) GO TO 221
      IF(KMA2.EQ.MAMAX) GO TO 222
      IF(KMA3.EQ.MAMAX) GO TO 223
221   L =KL2
      M =KMA2
      L1=KL3
      M1=KMA3
      L2=KL1
      M2=KMA1
       GO TO 230
222   L =KL1
      M =KMA1
      L1=KL3
      M1=KMA3
      L2=KL2
      M2=KMA2
       GO TO 230
223   L =KL1
      M =KMA1
      L1=KL2
      M1=KMA2
      L2=KL3
      M2=KMA3
       GO TO 230
C
230   IX=L2+M2
      IY=L+L1-M2
      IP=L2-M2
      IR=L-L1+M2
      IS=L1-M1
      ITMIN=MAX(0,-IR)
      ITMAX=MIN(IY,IP,IS)
      SUM=RCZ
      DO 301 IT=ITMIN,ITMAX
      A=(-RC1)**IT * TUFF(IX+IT)*TUFF(IY-IT)
      B=TUFF(IP-IT)*TUFF(IR+IT)*TUFF(IS-IT)*TUFF(IT)
      SUM=SUM+A/B
301   CONTINUE
C
      A=dble((2*L+1)*(2*L1+1)*(2*L2+1))
     &   *TUFF(IP)*TUFF(L+M)*TUFF(L1+M1)*TUFF(IS)
      B=RC2*TUFF(IX)*TUFF(L-M)
      PREF2=SQRT(A/B)
C
      A=(-RC1)**(IG-L-M1) * TUFF(IG2-2*L1)*TUFF(IG)
      B=TUFF(IG-L)*TUFF(IG-L1)*TUFF(IG-L2)*TUFF(IG2+1)
      PREF1=A/B
C
      GAMMA=PREF1*PREF2*SUM
C
C                                  GAUNT COEFFICIENT
      GACOR=GAMMA*DELTA
      RETURN
      END
C*******************
CXXX   GAUNTY   ****
C*******************
      SUBROUTINE GAUNTY(NL)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
C
C------------------------------------------------------
C   GAUNT'S COEFFICIENTS FOR REAL SPHERICAL HARMONICS:
C      INPUT - NL,  OUTPUT - COMMON/GHAR/, WHERE
C       GFRH(I,I1,I2) = INTEGRAL OVER THE UNIT SPHERE
C                    OF  Y(L,M) * Y(L1,M1) * Y(L2,M2).
C   HERE THE SIMPLE INDICES I,I1,I2 CORRESPOND TO THE
C   COMPOSED (L,M), (L1,M1), (L2,M2), RESPECTIVELY.
C   THE MAXIMAL VALUE OF L1 AND L2 IS (NL-1),
C   THE MAXIMAL VALUE OF L IS TWO TIMES GREATER.
C------------------------------------------------------
C
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
C
      LMAX=NL-1
      LMAXH=2*LMAX
C
      I2=0
      DO 328 L2=0,LMAX
      DO 329 M2=-L2,L2
      I2=I2+1
C
      I1=0
      DO 318 L1=0,LMAX
      DO 319 M1=-L1,L1
      I1=I1+1
C
      I=0
      DO 308 L=0,LMAXH
      DO 309 M=-L,L
      I=I+1
C
      GFRH(I,I1,I2)=GACOR(L,M,L1,M1,L2,M2)
C
309   CONTINUE
308   CONTINUE
319   CONTINUE
318   CONTINUE
329   CONTINUE
328   CONTINUE
C
      RETURN
      END
C*******************
CXXX    TUFF    ****
C*******************
      REAL*8 FUNCTION TUFF(N)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------
C        FACTORIAL OF A NON-NEGATIVE INTEGER
C----------------------------------------------
C
      TUFF= dble(1)
      IF(N.LE.1) RETURN
      DO 301 I=2,N
      TUFF=TUFF* dble(I)
301   CONTINUE
      RETURN
      END
C*******************
CXXX    SOPO    ****
C*******************
      SUBROUTINE SOPO(A,LDA,N,B,LDB,M,W)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C---------------------------------------------------------------
C      SOLUTION OF A SET OF LINEAR EQUATIONS WITH A POSITIVE
C      DEFINITE REAL SYMMETRIC MATRIX AND WITH MULTIPLE R.H.S.
C---------------------------------------------------------------
C
C     NA VSTUPU :
C        N ... DIMENZE MATICE
C        M ... POCET PRAVYCH STRAN
C        A(I,J) ... MATICE SOUSTAVY - JEN DOLNI TROJUHELNIK JE
C                   TREBA   (TJ.  I = 1, ... N,   J = 1, ... I)
C        B(I,K) ... PRAVE STRANY (K = 1, ... M,    I = 1, ... N)
C        W(I)  ...  PRACOVNI POLE DELKY N
C
C     NA VYSTUPU : A(I,J) JE PREPSANA SVYM CHOLESKYHO FAKTOREM
C                  B(I,K) .... PRISLUSNA RESENI
C
C
      DIMENSION A(LDA,N),B(LDB,M),W(N)
C
      DATA RC1/1.0D0/
C
C                                   ROZKLAD
      DO 310 L=1,N-1
      A(L,L)=SQRT(A(L,L))
      DUM=RC1/A(L,L)
      DO 311 I=L+1,N
      W(I)=A(I,L)*DUM
      A(I,L)=W(I)
311   CONTINUE
      DO 312 J=L+1,N
      DO 313 I=J,N
      A(I,J)=A(I,J)-W(I)*W(J)
313   CONTINUE
312   CONTINUE
310   CONTINUE
      A(N,N)=SQRT(A(N,N))
C
C                                 CYKLUS PRES PRAVE STRANY
      DO 320 MQ=1,M
C                                 INVERZE DOLNI TROJUH. MATICE
      DO 321 L=1,N-1
      W(L)=B(L,MQ)/A(L,L)
      DO 322 I=L+1,N
      B(I,MQ)=B(I,MQ)-A(I,L)*W(L)
322   CONTINUE
321   CONTINUE
      W(N)=B(N,MQ)/A(N,N)
C                                 INVERZE HORNI TROJUH. MATICE
      B(N,MQ)=W(N)/A(N,N)
      DO 331 L=N-1,1,-1
      DO 332 I=N,L+1,-1
      W(L)=W(L)-B(I,MQ)*A(I,L)
332   CONTINUE
      B(L,MQ)=W(L)/A(L,L)
331   CONTINUE
C
320   CONTINUE
C
      RETURN
      END
C*******************
CXXX    MASY    ****
C*******************
      SUBROUTINE MASY(N,ND,A,B,DEV)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C------------------------------------------------
C   SPECIAL SYMMETRIZATION OF MATRICES A AND B :
C  - CALCULATES DEVIATION OF A(TRANSPOSE) AND B
C  - CORRECTS A AND B TO SATISFY A(TRANSPOSE) = B
C------------------------------------------------
C
      DIMENSION A(ND,ND),B(ND,ND)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/
C
      DEV=RCZ
      DO 301 I=1,N
      DO 302 J=1,N
      EPS=ABS(A(I,J)-B(J,I))
      DEV=MAX(DEV,EPS)
      VAL=RCH*(A(I,J)+B(J,I))
      A(I,J)=VAL
      B(J,I)=VAL
302   CONTINUE
301   CONTINUE
      RETURN
      END
C*******************
CXXX    TBR     ****
C*******************
      SUBROUTINE TBR
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------
C   TB-STRUCTURE CONSTANTS IN THE REAL SPACE
C-----------------------------------------------
C
      PARAMETER(MNB=2)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=32)
C
      DIMENSION POL(3,MNCL),POLT(3),DPOL(3),
     &          WWMA(MNLSQ,MNLSQ),WWMB(MNLSQ,MNLSQ)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB),JBVA(MNCL,MNB),
     &            JTRA(3,MNCL,MNB),NSCL(MNB)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
C
      DATA RCZ/0.0D0/,SMALL/0.01D0/
C
      WRITE(IW6,111)
111   FORMAT(//' ************  TB-CONSTANTS IN',
     &    ' THE REAL SPACE  *************')
C
      NLSQ=NL**2
C
C------------------------------- CUT-OFF DISTANCE
      WMAX=CUTRAT*AWS
      DMAXSQ=WMAX**2
      WMIN=SMALL*AWS
      DMINSQ=WMIN**2
C
C------------------------------- LOOP OVER CENTRAL ATOMS
C                                FOR SMALL CLUSTERS
      DO 310 IB=1,NB
C
C------------------------- SELECTION OF THE SMALL CLUSTER
      JCL=1
      DO 312 I=1,3
      POL(I,1)=POS(I,IB)
312   CONTINUE
      JBVA(1,IB)=IB
      JTRA(1,1,IB)=0
      JTRA(2,1,IB)=0
      JTRA(3,1,IB)=0
C
      DO 313 JB=1,NB
      DO 3131 JTRA1=-NMTR,NMTR
      DO 3132 JTRA2=-NMTR,NMTR
      DO 3133 JTRA3=-NMTR,NMTR
C
      DO 314 I=1,3
      POLT(I)=POS(I,JB)
     &       +JTRA1*VBR(I,1)+JTRA2*VBR(I,2)+JTRA3*VBR(I,3)
314   CONTINUE
      DO 315 I=1,3
      DPOL(I)=POLT(I)-POL(I,1)
      IF (ABS(DPOL(I)).GT.WMAX) GO TO 3133
315   CONTINUE
      SUM=DPOL(1)**2+DPOL(2)**2+DPOL(3)**2
      IF (SUM.LT.DMINSQ.OR.SUM.GT.DMAXSQ) GO TO 3133
C
      JCL=JCL+1
      IF (JCL.GT.MNCL) GO TO 291
      DO 316 I=1,3
      POL(I,JCL)=POLT(I)
316   CONTINUE
      JBVA(JCL,IB)=JB
      JTRA(1,JCL,IB)=JTRA1
      JTRA(2,JCL,IB)=JTRA2
      JTRA(3,JCL,IB)=JTRA3
C
3133  CONTINUE
3132  CONTINUE
3131  CONTINUE
313   CONTINUE
C
      NCL=JCL
      NSCL(IB)=NCL
      WRITE(IW6,115) IB,NCL
115   FORMAT(/'    SITE:    IB=',I3,
     &       '        CLUSTER SIZE= ',I3)
C
C--------------------------------- TB-CONSTANTS FOR THE
C                                  SMALL CLUSTER AND
C                                  THEIR STORAGE
      CALL TBCL(NL,NCL,POL,AWS)
C
      DO 320 JCL=1,NCL
      DO 321 JQ=1,NLSQ
      DO 322 IQ=1,NLSQ
      STR(IQ,JQ,JCL,IB)=SSC(IQ,JQ,JCL)
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
310   CONTINUE
C
C----------------------------------------- SYMMETRIZATION
C
      DEV=RCZ
C
      DO 330 IB=1,NB
      DO 331 JCL=1,NSCL(IB)
      JB=JBVA(JCL,IB)
      DO 332 ICL=1,NSCL(JB)
      IF (JBVA(ICL,JB).NE.IB) GO TO 332
      ISUM=JTRA(1,JCL,IB)+JTRA(1,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
      ISUM=JTRA(2,JCL,IB)+JTRA(2,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
      ISUM=JTRA(3,JCL,IB)+JTRA(3,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
C
      DO 3358 JQ=1,NLSQ
      DO 3359 IQ=1,NLSQ
      WWMA(IQ,JQ)=STR(IQ,JQ,JCL,IB)
3359  CONTINUE
3358  CONTINUE
      DO 3368 JQ=1,NLSQ
      DO 3369 IQ=1,NLSQ
      WWMB(IQ,JQ)=STR(IQ,JQ,ICL,JB)
3369  CONTINUE
3368  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3378 JQ=1,NLSQ
      DO 3379 IQ=1,NLSQ
      STR(IQ,JQ,JCL,IB)=WWMA(IQ,JQ)
3379  CONTINUE
3378  CONTINUE
      DO 3388 JQ=1,NLSQ
      DO 3389 IQ=1,NLSQ
      STR(IQ,JQ,ICL,JB)=WWMB(IQ,JQ)
3389  CONTINUE
3388  CONTINUE
C
332   CONTINUE
331   CONTINUE
330   CONTINUE
C
      WRITE(IW6,125) DEV
125   FORMAT(/15X,'    ------ DEVIATION=',G12.4)
C
C--------------------------------- PRINT OF DIAGONAL
C                                  ELEMENTS
      WRITE(IW6,150)
150   FORMAT(/5X,'***  DIAGONAL ELEMENTS OF ',
     &    'TB-STRUCTURE CONSTANTS :')
      DO 350 IB=1,NB
      WRITE(IW6,151) IB
151   FORMAT(/'        SITE:       IB=',I3)
      WRITE(IW6,104) (STR(IQ,IQ,1,IB),IQ=1,NLSQ)
104   FORMAT(1X,4G15.7)
350   CONTINUE
C
      RETURN
C
291   WRITE(IW6,191) JCL
191   FORMAT(/' **** ERROR IN TBR :'/10X,
     &   ' JCL GREATER THAN MNCL,  JCL=',I3)
      STOP
      END
C*******************
CXXX    GIBZ    ****
C*******************
      SUBROUTINE GIBZ
C
C**************************************************************
C  GENERATES NETWORK OF K-POINTS IN IRREDUCIBLE BRILLOUIN ZONE
C**************************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNBZ=160002)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
      DATA RCH/0.5D0/,RC3/3.0D0/,RC6/6.0D0/
      DATA DUM1/0.8D0/,DUM2/1.2D0/
C
      PI=RC4*ATAN(RC1)
C
      IF(NSYM.EQ.0) GO TO 200
      IF(NSYM.EQ.1) GO TO 201
      IF(NSYM.EQ.2) GO TO 202
      IF(NSYM.EQ.3) GO TO 2012
      IF(NSYM.EQ.4) GO TO 204
      IF(NSYM.EQ.5) GO TO 205
C
C-------------------------- NSYM=0: GENERAL CASE
C
200   IBZ=0
      TW=RCZ
      NK1=NK
      NK2=NK
      NK3=NK
      DO 301 I1=-NK1+1,NK1
      FX1=(dble(I1)-RCH)/dble(2*NK1)
      DO 302 I2=-NK2+1,NK2
      FX2=(dble(I2)-RCH)/dble(2*NK2)
      DO 303 I3=1,NK3
      FX3=(dble(I3)-RCH)/dble(2*NK3)
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      DO 305 J=1,3
      AKBZ(J,IBZ)=FX1*VBG(J,1)+FX2*VBG(J,2)+FX3*VBG(J,3)
305   CONTINUE
303   CONTINUE
302   CONTINUE
301   CONTINUE
      GO TO 290
C
201   CONTINUE
      IF(INVE.EQ.1) GO TO 2011
      IF(INVE.EQ.2) GO TO 2012
      IF(INVE.EQ.3) GO TO 2013
C
C-------------------------- NSYM=1, INVE=1: SC LATTICE
2011  IBZ=0
      TW=RCZ
      A=ABS(VBR(1,1))+ABS(VBR(2,1))+ABS(VBR(3,1))
      DK=PI/(A*dble(NK))
      DO 3101 I1=1,NK
      DO 3102 I2=1,I1
      DO 3103 I3=1,I2
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC6
      IF(I1.EQ.I2.OR.I2.EQ.I3) WKBZ(IBZ)=RC3
      IF(I1.EQ.I3) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(dble(I1)-RCH)*DK
      AKBZ(2,IBZ)=(dble(I2)-RCH)*DK
      AKBZ(3,IBZ)=(dble(I3)-RCH)*DK
3103  CONTINUE
3102  CONTINUE
3101  CONTINUE
      GO TO 290
C
C--------------------------- NSYM=1, INVE=2: FCC LATTICE
C                                  (ALSO FOR NSYM=3)
2012  IF(MOD(NK,2).NE.0) GO TO 295
      IBZ=0
      TW=RCZ
      A=ABS(VBR(1,1))+ABS(VBR(2,1))+ABS(VBR(3,1))
      DK=RC2*PI/(A*dble(NK))
      DO 3141 I1=1,NK
      DO 3142 I2=1,I1
      DO 3143 I3=1,I2
      IF(2*(I1+I2+I3).GT.3*(NK+1)) GO TO 3143
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC6
      IF(I1.EQ.I2.OR.I2.EQ.I3) WKBZ(IBZ)=RC3
      IF(I1.EQ.I3) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(dble(I1)-RCH)*DK
      AKBZ(2,IBZ)=(dble(I2)-RCH)*DK
      AKBZ(3,IBZ)=(dble(I3)-RCH)*DK
3143  CONTINUE
3142  CONTINUE
3141  CONTINUE
      GO TO 290
C
C--------------------------- NSYM=1, INVE=3: BCC LATTICE
2013  IF(MOD(NK,2).NE.0) GO TO 295
      IBZ=0
      TW=RCZ
      A=ABS(VBR(1,1))+ABS(VBR(2,1))
      DK=RC2*PI/(A*dble(NK))
      DO 3181 I1=1,NK
      DO 3182 I2=1,I1
      IF(I1+I2.GT.NK+1) GO TO 3182
      DO 3183 I3=1,I2
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC6
      IF(I1.EQ.I2.OR.I2.EQ.I3.OR.I1+I2.EQ.NK+1) WKBZ(IBZ)=RC3
      IF(I2.EQ.I3.AND.I1+I2.EQ.NK+1) WKBZ(IBZ)=RC1
      IF(I1.EQ.I3) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(dble(I1)-RCH)*DK
      AKBZ(2,IBZ)=(dble(I2)-RCH)*DK
      AKBZ(3,IBZ)=(dble(I3)-RCH)*DK
3183  CONTINUE
3182  CONTINUE
3181  CONTINUE
      GO TO 290
C
C-------------------------- NSYM=2 :  HEXAGONAL LATTICE
202   IBZ=0
      TW=RCZ
      A=VBR(1,1)
      C=VBR(3,3)
      NKX=2*NK
      NKZ=NK
      DKX=RC2*PI/(A*dble(NKX))
      DKY=DKX/SQRT(RC3)
      DKZ=PI/(C*dble(NKZ))
      DO 321 IZ=1,NKZ
      DO 322 NU=1,2
      ANU3=dble(NU)/RC3
      DO 323 I=1,NKX
      DO 324 J=I,NKX
      IY=I+2*J+NU-3
      IF(IY.GT.NKX) GO TO 324
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC6
      IF(I.EQ.J.OR.IY.EQ.NKX) WKBZ(IBZ)=RC3
      IF(I.EQ.J.AND.IY.EQ.NKX) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(dble(I-1)+ANU3)*DKX
      AKBZ(2,IBZ)=dble(IY)*DKY
      AKBZ(3,IBZ)=(dble(IZ)-RCH)*DKZ
324   CONTINUE
323   CONTINUE
322   CONTINUE
321   CONTINUE
      GO TO 290
C
204   CONTINUE
      IF(INVE.EQ.1) GO TO 2041
      IF(INVE.EQ.2) GO TO 2042
C
C-------------------------- NSYM=4, INVE=1 :
C                                   SIMPLE TETRAGONAL (I)
2041  IBZ=0
      TW=RCZ
      A=VBR(1,1)
      C=VBR(3,3)
      COA=C/A
      NKZ=NK
      IF(COA.LT.DUM1.OR.COA.GT.DUM2) NKZ=INT(dble(NK)/COA)
      DK=PI/(A*dble(NK))
      DKZ=PI/(C*dble(NKZ))
      DO 341 I1=1,NK
      DO 342 I2=1,I1
      DO 343 I3=1,NKZ
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC2
      IF(I1.EQ.I2) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(dble(I1)-RCH)*DK
      AKBZ(2,IBZ)=(dble(I2)-RCH)*DK
      AKBZ(3,IBZ)=(dble(I3)-RCH)*DKZ
343   CONTINUE
342   CONTINUE
341   CONTINUE
      GO TO 290
C
C-------------------------- NSYM=4, INVE=2 :
C                                   SIMPLE TETRAGONAL (II)
2042  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      C=VBR(3,3)
      SR2=SQRT(RC2)
      NKZ=INT(dble(NK)*P*SR2/C)
      DK=PI/(RC2*P*dble(NK))
      DKZ=PI/(C*dble(NKZ))
      DO 345 I=1,NK
      DO 346 J=1,I
      DO 347 I3=1,NKZ
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC2
      IF(I.EQ.J) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=dble(I-J)*DK
      AKBZ(2,IBZ)=dble(I+J-1)*DK
      AKBZ(3,IBZ)=(dble(I3)-RCH)*DKZ
347   CONTINUE
346   CONTINUE
345   CONTINUE
      GO TO 290
C
C-------------------------- NSYM=5 : SIMPLE ORTHORHOMBIC
205   IBZ=0
      TW=RCZ
      A=VBR(1,1)
      B=VBR(2,2)
      C=VBR(3,3)
      BOA=B/A
      COA=C/A
      NKX=NK
      NKY=NK
      IF(BOA.LT.DUM1.OR.BOA.GT.DUM2) NKY=INT(dble(NK)/BOA)
      NKZ=NK
      IF(COA.LT.DUM1.OR.COA.GT.DUM2) NKZ=INT(dble(NK)/COA)
      DKX=PI/(A*dble(NKX))
      DKY=PI/(B*dble(NKY))
      DKZ=PI/(C*dble(NKZ))
      DO 351 I1=1,NKX
      DO 352 I2=1,NKY
      DO 353 I3=1,NKZ
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(dble(I1)-RCH)*DKX
      AKBZ(2,IBZ)=(dble(I2)-RCH)*DKY
      AKBZ(3,IBZ)=(dble(I3)-RCH)*DKZ
353   CONTINUE
352   CONTINUE
351   CONTINUE
      GO TO 290
C
290   NBZ=IBZ
      DO 390 IBZ=1,NBZ
      WKBZ(IBZ)=WKBZ(IBZ)/TW
390   CONTINUE
      WRITE(IW6,190) NBZ
190   FORMAT(//3X,'*****  NUMBER OF K-POINTS:   NBZ=',I6)
      RETURN
C
295   WRITE(IW6,195) NK
195   FORMAT(/3X,' **** ERROR IN GIBZ: ',
     &    ' NK MUST BE EVEN,    NK=',I4)
      STOP
C
299   WRITE(IW6,199) IBZ
199   FORMAT(/3X,' **** ERROR IN GIBZ: ',
     &    ' IBZ EXCEEDED MNBZ,  IBZ=',I6)
      STOP
      END
C*******************
CXXX    TBK     ****
C*******************
      SUBROUTINE TBK
C
C***********************************************
C   BLOCH TRANSFORM OF TB-STRUCTURE CONSTANTS
C***********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNL=3)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=32)
C
      DIMENSION TRAN(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB),JBVA(MNCL,MNB),
     &            JTRA(3,MNCL,MNB),NSCL(MNB)
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      WRITE(IW6,111)
111   FORMAT(//3X,' **********  BLOCH TRANSFORM ',
     &   'OF TB-CONSTANTS  ********** '/
     &   /10X,'    CHECK OF HERMITICITY : ')
C
      DO 301 IBZ=1,MNBZ
      DO 302 II=1,MNBLSQ
      DO 303 JJ=1,MNBLSQ
      ZSK(JJ,II,IBZ)=ZZ
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
      DO 331 IBZ=1,NBZ
      DO 332 IB=1,NB
      II0=(IB-1)*NLSQ
C
      DO 333 JCL=1,NSCL(IB)
      JB=JBVA(JCL,IB)
      JTRA1=JTRA(1,JCL,IB)
      JTRA2=JTRA(2,JCL,IB)
      JTRA3=JTRA(3,JCL,IB)
C
      JJ0=(JB-1)*NLSQ
      DO 334 I=1,3
      TRAN(I)=JTRA1*VBR(I,1)+JTRA2*VBR(I,2)+JTRA3*VBR(I,3)
334   CONTINUE
      ARG=RCZ
      DO 335 I=1,3
      ARG=ARG+AKBZ(I,IBZ)*TRAN(I)
335   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      DO 336 IQ=1,NLSQ
      DO 337 JQ=1,NLSQ
      ZSK(JJ0+JQ,II0+IQ,IBZ) = ZSK(JJ0+JQ,II0+IQ,IBZ)
     &                       + ZEFA * STR(JQ,IQ,JCL,IB)
337   CONTINUE
336   CONTINUE
C
333   CONTINUE
332   CONTINUE
331   CONTINUE
C
C                                   CHECK OF HERMITICITY
      DEV=RCZ
      DO 341 IBZ=1,NBZ
      DO 342 II=1,NBLSQ
      DO 343 JJ=II,NBLSQ
      ZW1=ZSK(II,JJ,IBZ)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSK(JJ,II,IBZ)
      W2X=DREAL(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
343   CONTINUE
342   CONTINUE
341   CONTINUE
      WRITE(IW6,140) DEV
140   FORMAT(/7X,'            ------ DEVIATION=',G12.4)
C
      RETURN
      END
C*******************
CXXX    SETPF   ****
C*******************
      SUBROUTINE SETPF(IE)
C
C***************************************************
C    CALCULATION OF POTENTIAL FUNCTIONS AND RELATED
C    QUANTITIES FOR ONE ENERGY
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNLSQ=MNL**2)
C
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/SCREEN/ QSCR(MNL)
C
      DATA RC1/1.0D0/
C
      ZE1=ZCN(IE)
C       write(*,*) ZE1
C
      DO 311 IS=1,NS
      DO 312 IA=1,NA
      DO 313 IL=1,NL
      PPC1=PPC(IL,IS,IA)
      PPD1=PPD(IL,IS,IA)
      PPQ1=PPQ(IL,IS,IA)
      ALF1=QSCR(IL)
      CALL PLMZ(PPC1,PPD1,PPQ1,ALF1,ZE1,ZPF1,ZLA1,ZMU1)
      ARG = RC1/FINY(IL,IS,IA)
     &    - WS(IA)*FINYD(IL,IS,IA)*(DNY(IL,IS,IA)+dble(IL))
      YNAM1=SIGN(RC1,ARG)
      ZMU1=YNAM1*ZMU1
C
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 316 I=ISTA,IFIN
      ZPF(I,IA,IS)=ZPF1
      ZMU(I,IA,IS)=ZMU1
      ZLA(I,IA,IS)=ZLA1
316   CONTINUE
313   CONTINUE
312   CONTINUE
311   CONTINUE
C
      RETURN
      END
C*******************
CXXX    ONGF    ****
C*******************
      SUBROUTINE ONGF(IS)
C
C***************************************************
C    CALCULATION OF ON-SITE BLOCKS OF THE GF MATRIX
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZWA(MNBLSQ,MNBLSQ),ZWG(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 IB=1,MNB
      DO 302 J=1,MNLSQ
      DO 303 I=1,MNLSQ
      ZGFG(I,J,IB,IS)=ZZ
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
C-------------------------------- IBZ LOOP
      DO 312 IBZ=1,NBZ
C                               MATRIX P(Z)-S(K)
      DO 314 JJ=1,NBLSQ
      DO 315 II=1,NBLSQ
      ZWA(II,JJ)=-ZSK(II,JJ,IBZ)
315   CONTINUE
314   CONTINUE
      DO 316 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 317 J=1,NLSQ
      DO 318 I=1,NLSQ
      ZWA(II0+I,II0+J)=ZWA(II0+I,II0+J)+ZCPF(I,J,IB,IS)
318   CONTINUE
317   CONTINUE
316   CONTINUE
C                                      INVERSION
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
C                                    ON-SITE BLOCKS
      DO 324 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 326 J=1,NLSQ
      DO 327 I=1,NLSQ
      ZGFG(I,J,IB,IS) = ZGFG(I,J,IB,IS)
     &            + WKBZ(IBZ) * ZWA(II0+I,II0+J)
327   CONTINUE
326   CONTINUE
324   CONTINUE
C---------------------------------- END OF IBZ-LOOP
312   CONTINUE
C
C---------------------------------- SYMMETRIZATION OF THE
C                                   ON-SITE BLOCKS
      DO 350 IB=1,NB
      DO 352 J=1,NLSQ
      DO 353 I=1,NLSQ
      ZWG(I,J)=ZGFG(I,J,IB,IS)
353   CONTINUE
352   CONTINUE
C
      CALL SYMGF(ZWG,NSYM,NL)
C
      DO 356 J=1,NLSQ
      DO 357 I=1,NLSQ
      ZGFG(I,J,IB,IS)=ZWG(I,J)
357   CONTINUE
356   CONTINUE
350   CONTINUE
C
      RETURN
      END
C*******************
CXXX    CPAIT   ****
C*******************
      SUBROUTINE CPAIT(IE,IS)
C
C******************************************
C   PERFORMS CPA-SELFCONSISTENCY LOOP AND
C    CALCULATES GREEN FUNCTIONS FOR DOS
C******************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/GDOS/ ZGDOS(MNL,MNA,MNE,MNS)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2
      ZZ=DCMPLX(RCZ,RCZ)
C
      QMIXI=QMIXC
      IF(IRAND.EQ.0) QMIXI=RC1
C
      IITER=0
C                             STARTING OMEGA
      IF(IE.EQ.1) THEN
      DO 3001 IB=1,MNB
      DO 3002 J=1,MNLSQ
      DO 3003 I=1,MNLSQ
      ZOMG(I,J,IB,IS)=ZZ
3003  CONTINUE
3002  CONTINUE
3001  CONTINUE
      END IF
C
210   IITER=IITER+1
      IF(IITER.GT.NMITC) GO TO 290
C
      NONC=0
C
C::::::::::::::::::::::::::::::::::::   STEP 1
C
      IA=0
      DO 3100 IB=1,NB
C
      DO 3116 J=1,NLSQ
      DO 3117 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IB,IS)
3117  CONTINUE
3116  CONTINUE
      DO 3124 J=1,NLSQ
      DO 3125 I=1,NLSQ
      ZFI1(I,J)=ZZ
3125  CONTINUE
3124  CONTINUE
C
      DO 3130 IC=1,NC(IB)
      IA=IA+1
C
      DO 3134 J=1,NLSQ
      DO 3135 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3135  CONTINUE
3134  CONTINUE
      DO 3136 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZPF(I,IA,IS)
3136  CONTINUE
C
      CALL CINVC(ZW,NLSQ,MNLSQ)
C
      DO 3141 J=1,NLSQ
      DO 3142 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+CON(IA)*ZW(I,J)
3142  CONTINUE
3141  CONTINUE
C
3130  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
C
      DO 3150 J=1,NLSQ
      DO 3151 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3151  CONTINUE
3150  CONTINUE
C
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3154 J=1,NLSQ
      DO 3155 I=1,NLSQ
      ZCPF(I,J,IB,IS)=ZCP1(I,J)
3155  CONTINUE
3154  CONTINUE
C
3100  CONTINUE
C
C::::::::::::::::::::::::::::::::::::   STEP 2
C
      CALL ONGF(IS)
C
C::::::::::::::::::::::::::::::::::::   STEP 3
C
      DO 3300 IB=1,NB
C
      DO 3304 J=1,NLSQ
      DO 3305 I=1,NLSQ
      ZFI1(I,J)=ZGFG(I,J,IB,IS)
3305  CONTINUE
3304  CONTINUE
      DO 3308 J=1,NLSQ
      DO 3309 I=1,NLSQ
      ZCP1(I,J)=ZCPF(I,J,IB,IS)
3309  CONTINUE
3308  CONTINUE
      DO 3312 J=1,NLSQ
      DO 3313 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IB,IS)
3313  CONTINUE
3312  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
C
      DO 3316 J=1,NLSQ
      DO 3317 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
3317  CONTINUE
3316  CONTINUE
C
      SUM=RCZ
      SUM1=RCZ
      DO 3320 J=1,NLSQ
      DO 3321 I=1,NLSQ
      XW=DREAL(ZW(I,J))
      YW=DIMAG(ZW(I,J))
      SUM=SUM+ ABS(XW)+ ABS(YW)
      XW1=DREAL(ZOM1(I,J))
      YW1=DIMAG(ZOM1(I,J))
      SUM1=SUM1+ ABS(XW1)+ ABS(YW1)
3321  CONTINUE
3320  CONTINUE
      IF(SUM.GT.(DIFMC*SUM1)) NONC=NONC+1
C
      DO 3324 J=1,NLSQ
      DO 3325 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXI*ZW(I,J)
3325  CONTINUE
3324  CONTINUE
C
      CALL SYMGF(ZOM1,NSYM,NL)
C
      DO 3328 J=1,NLSQ
      DO 3329 I=1,NLSQ
      ZOMG(I,J,IB,IS)=ZOM1(I,J)
3329  CONTINUE
3328  CONTINUE
C
3300  CONTINUE
C                             TEST OF CONVERGENCE
C
      IF(IRAND.EQ.1.AND.NONC.GT.0) GO TO 210
      WRITE(IW6,111) IE,IS,IITER
111   FORMAT(9X,'IE=',I4,'    IS=',I1,
     &      '         CPA : ',I4,' ITER.')
C
C------------------------------------ GF FOR DOS
C
      IA=0
      DO 350 IB=1,NB
C
      DO 3516 J=1,NLSQ
      DO 3517 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IB,IS)
3517  CONTINUE
3516  CONTINUE
C
      DO 353 IC=1,NC(IB)
      IA=IA+1
C
      DO 3534 J=1,NLSQ
      DO 3535 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3535  CONTINUE
3534  CONTINUE
      DO 3536 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZPF(I,IA,IS)
3536  CONTINUE
C
      CALL CINVC(ZW,NLSQ,MNLSQ)
C
      DO 355 IL=1,NL
      ZSUM=ZZ
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 3552 I=ISTA,IFIN
      ZSUM = ZSUM + ZLA(I,IA,IS)
     &            + ZMU(I,IA,IS)**2 * ZW(I,I)
3552  CONTINUE
      ZGDOS(IL,IA,IE,IS)=ZSUM
355   CONTINUE
C
353   CONTINUE
350   CONTINUE
C
      RETURN
C
C------------------------------ NO CONVERGENCY
C
290   WRITE(IW6,190) NMITC,IE,IS
190   FORMAT(//' **** CPAIT:  NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   '      ---  FOR    IE=',I4,'    IS=',I2)
      STOP
      END
C*******************
CXXX    PLMZ    ****
C*******************
      SUBROUTINE PLMZ(C,DELTA,GAMMA,ALPHA,ZE,ZP,ZL,ZM)
C
C******************************************************
C          SCREENED POTENTIAL FUNCTION P(ZE) AND
C          RELATED FUNCTIONS LAMBDA(ZE), MU(ZE)
C******************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      GMA=GAMMA-ALPHA
      SQD=SQRT(DELTA)
      ZEMC=ZE-C
      ZDEN=DELTA+GMA*ZEMC
      ZP=ZEMC/ZDEN
      ZL=GMA/ZDEN
      ZM=SQD/ZDEN
      RETURN
      END
C*******************
CXXX    CINVC   ****
C*******************
      SUBROUTINE CINVC(ZA,NDIM,NTOT)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C  **************************************************************
C  * INVERT COMPLEX MATRIX ZA BY GAUSS ELIMINATION WITH PARTIAL *
C  *       PIVOTING   -   USING COMPLEX ARITHMETIC              *
C  **************************************************************
C  *    IMPORTANT: NDIM <= NTOT, NDIM <= NMAX, WHERE NMAX IS    *
C  *    DIMENSION OF ARRAYS USED FOR EVALUATION OF INV(ZA).     *
C  **************************************************************
C  *        P. MARKSTEINER, UNIVERSITY OF VIENNA, 1990          *
C  **************************************************************
C---------------------------------
      PARAMETER(MNB=2)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C---------------------------------
      PARAMETER(NMAX=MNBLSQ)
C
      DIMENSION ZA(NTOT,NDIM)
      DIMENSION AUX(NMAX),INDA(NMAX)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      ZZ=DCMPLX(RCZ,RCZ)
      Z1=DCMPLX(RC1,RCZ)
C
      N=NDIM
      IF(N.GT.NMAX) STOP '*** CINVC: NMAX TOO SMALL'
C
      DO 31 I=1,N
!SMP$ DO SERIAL      
      DO 32 J=I,N
      AUX(J)=ZA(I,J)*DCONJG(ZA(I,J))
32    CONTINUE
      K=IDAMA1(AUX,N,I)
      INDA(I)=K
      IF(K.NE.I) THEN
      DO 33 J=1,N
      ZQ=ZA(J,I)
      ZA(J,I)=ZA(J,K)
      ZA(J,K)=ZQ
33    CONTINUE
      END IF
      ZQ=Z1/ZA(I,I)
      ZA(I,I)=Z1
      DO 35 J=1,N
      ZA(J,I)=ZQ*ZA(J,I)
35    CONTINUE
      DO 36 J=1,N
      IF(J.NE.I) THEN
      ZQ=-ZA(I,J)
      ZA(I,J)=ZZ
      DO 37 K=1,N
      ZA(K,J)=ZA(K,J)+ZQ*ZA(K,I)
37    CONTINUE
      END IF
36    CONTINUE
31    CONTINUE
      DO 38 I=N,1,-1
      J=INDA(I)
      IF(I.NE.J) THEN
      DO 39 K=1,N
      ZQ=ZA(I,K)
      ZA(I,K)=ZA(J,K)
      ZA(J,K)=ZQ
39    CONTINUE
      END IF
38    CONTINUE
      RETURN
      END
C*******************
CXXX   IDAMA1   ****
C*******************
      INTEGER FUNCTION IDAMA1(DX,ND,N1)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------
C     FINDS THE INDEX OF MAX. ELEMENT OF DX(I)
C     THE MAXIMUM IS SEARCHED OVER N1.LE.I.LE.ND
C-----------------------------------------------
C
      DIMENSION DX(ND)
C
      IDAMA1=N1
      IF(N1.EQ.ND) RETURN
      WMAX=DX(N1)
      DO 30 I=N1+1,ND
      IF(DX(I).LE.WMAX) GO TO 30
      IDAMA1=I
      WMAX=DX(I)
30    CONTINUE
      RETURN
      END
C*******************
CXXX    SYMGF   ****
C*******************
      SUBROUTINE SYMGF(ZA,NSYM,NL)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
C
C----------------------------------------------------
C   SYMMETRIZATION OF LOCAL GF MATRIX ACCORDING TO
C   THE POINT GROUP SYMMETRY
C----------------------------------------------------
C
      DIMENSION ZA(MNLSQ,MNLSQ),ZG(16,16),ZP(16,16)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/
C
      DATA        C333/0.333333333333333D0/,
     &            C208/0.208333333333333D0/,
     &            C161/0.161374306091976D0/,
     &            C263/0.263523138347365D0/,
     &            C125/0.125000000000000D0/,
     &            C204/0.204124145231932D0/,
     &            C177/0.177083333333333D0/,
     &            C040/0.040343576522994D0/,
     &            C156/0.156250000000000D0/
C
      NLSQ=NL**2
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 J=1,NLSQ
      DO 302 I=1,NLSQ
      ZP(I,J)=RCH*(ZA(I,J)+ZA(J,I))
302   CONTINUE
301   CONTINUE
C
C----------------------------------  GENERAL CASE
C
      IF(NSYM.GT.0) GO TO 2000
C
      DO 304 J=1,NLSQ
      DO 305 I=1,NLSQ
      ZA(I,J)=ZP(I,J)
305   CONTINUE
304   CONTINUE
      RETURN
C
C----------------------------------  SPECIAL CASES
C
2000  CONTINUE
      DO 307 J=1,NLSQ
      DO 308 I=1,NLSQ
      ZG(I,J)=ZZ
308   CONTINUE
307   CONTINUE
C
      IF(NSYM.EQ.1) GO TO 201
      IF(NSYM.EQ.2) GO TO 202
      IF(NSYM.EQ.3) GO TO 203
      IF(NSYM.EQ.4) GO TO 204
      IF(NSYM.EQ.5) GO TO 205
C
C                              CUBIC SYMMETRY
C
201   ZG(1,1)=ZP(1,1)
      ZG(2,2)=C333*(ZP(2,2)+ZP(3,3)+ZP(4,4))
      ZG(3,3)=ZG(2,2)
      ZG(4,4)=ZG(2,2)
      ZG(5,5)=C333*(ZP(5,5)+ZP(6,6)+ZP(8,8))
      ZG(6,6)=ZG(5,5)
      ZG(7,7)=RCH*(ZP(7,7)+ZP(9,9))
      ZG(8,8)=ZG(5,5)
      ZG(9,9)=ZG(7,7)
       IF(NL.EQ.3) GO TO 270
      ZG(10,2)=   C208*ZP(10,2) + C161*ZP(12,2)
     &          - C263*ZP(13,3) + C161*ZP(14,4)
     &          - C208*ZP(16,4)
      ZG(12,2)=   C161*ZP(10,2) + C125*ZP(12,2)
     &          - C204*ZP(13,3) + C125*ZP(14,4)
     &          - C161*ZP(16,4)
      ZG(13,3)= - C263*ZP(10,2) - C204*ZP(12,2)
     &          + C333*ZP(13,3) - C204*ZP(14,4)
     &          + C263*ZP(16,4)
      ZG(14,4)=ZG(12,2)
      ZG(16,4)=-ZG(10,2)
      ZG(10,10)=   C177*ZP(10,10) + C040*ZP(12,10)
     &           + C040*ZP(10,12) + C156*ZP(12,12)
     &           + C208*ZP(13,13) + C156*ZP(14,14)
     &           - C040*ZP(16,14) + C125*ZP(15,15)
     &           - C040*ZP(14,16) + C177*ZP(16,16)
      ZG(12,10)=   C040*ZP(10,10) + C156*ZP(12,10)
     &           + C156*ZP(10,12) - C040*ZP(12,12)
     &           + C161*ZP(13,13) - C040*ZP(14,14)
     &           - C156*ZP(16,14) - C161*ZP(15,15)
     &           - C156*ZP(14,16) + C040*ZP(16,16)
      ZG(11,11)=ZP(11,11)
      ZG(12,12)=   C156*ZP(10,10) - C040*ZP(12,10)
     &           - C040*ZP(10,12) + C177*ZP(12,12)
     &           + C125*ZP(13,13) + C177*ZP(14,14)
     &           + C040*ZP(16,14) + C208*ZP(15,15)
     &           + C040*ZP(14,16) + C156*ZP(16,16)
      ZG(13,13)=   C208*ZP(10,10) + C161*ZP(12,10)
     &           + C161*ZP(10,12) + C125*ZP(12,12)
     &           + C333*ZP(13,13) + C125*ZP(14,14)
     &           - C161*ZP(16,14) - C161*ZP(14,16)
     &           + C208*ZP(16,16)
      ZG(14,14)=ZG(12,12)
      ZG(16,14)=-ZG(12,10)
      ZG(15,15)=   C125*ZP(10,10) - C161*ZP(12,10)
     &           - C161*ZP(10,12) + C208*ZP(12,12)
     &           + C208*ZP(14,14) + C161*ZP(16,14)
     &           + C333*ZP(15,15) + C161*ZP(14,16)
     &           + C125*ZP(16,16)
      ZG(16,16)=ZG(10,10)
       GO TO 270
C
C   HEXAGONAL SYMMETRY - GENERATED BY 3-FOLD ROTATION (AXIS Z)
C    AND TWO MIRROR REFLECTIONS (PLANES X-Y AND Y-Z)
C
202   ZG(1,1)=ZP(1,1)
      ZG(7,1)=ZP(7,1)
      ZG(2,2)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(9,2)=RCH*(ZP(9,2)+ZP(5,4))
      ZG(3,3)=ZP(3,3)
      ZG(4,4)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(5,4)=RCH*(ZP(9,2)+ZP(5,4))
      ZG(5,5)=RCH*(ZP(5,5)+ZP(9,9))
      ZG(6,6)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(7,7)=ZP(7,7)
      ZG(8,8)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(9,9)=RCH*(ZP(5,5)+ZP(9,9))
       IF(NL.EQ.3) GO TO 270
      ZG(10,1)=ZP(10,1)
      ZG(12,2)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(13,3)=ZP(13,3)
      ZG(14,4)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(14,5)=RCH*(ZP(14,5)+ZP(12,9))
      ZG(15,6)=RCH*(ZP(15,6)+ZP(11,8))
      ZG(10,7)=ZP(10,7)
      ZG(11,8)=RCH*(ZP(15,6)+ZP(11,8))
      ZG(12,9)=RCH*(ZP(14,5)+ZP(12,9))
      ZG(10,10)=ZP(10,10)
      ZG(11,11)=RCH*(ZP(11,11)+ZP(15,15))
      ZG(12,12)=RCH*(ZP(12,12)+ZP(14,14))
      ZG(13,13)=ZP(13,13)
      ZG(14,14)=RCH*(ZP(12,12)+ZP(14,14))
      ZG(15,15)=RCH*(ZP(11,11)+ZP(15,15))
      ZG(16,16)=ZP(16,16)
       GO TO 270
C
C                          TETRAHEDRAL SYMMETRY
C
203   ZG(1,1)=ZP(1,1)
      ZG(2,2)=C333*(ZP(2,2)+ZP(3,3)+ZP(4,4))
      ZG(8,2)=C333*(ZP(8,2)+ZP(5,3)+ZP(6,4))
      ZG(3,3)=ZG(2,2)
      ZG(5,3)=ZG(8,2)
      ZG(4,4)=ZG(2,2)
      ZG(6,4)=ZG(8,2)
      ZG(5,5)=C333*(ZP(5,5)+ZP(6,6)+ZP(8,8))
      ZG(6,6)=ZG(5,5)
      ZG(7,7)=RCH*(ZP(7,7)+ZP(9,9))
      ZG(8,8)=ZG(5,5)
      ZG(9,9)=ZG(7,7)
       IF(NL.EQ.3) GO TO 270
      ZG(11,1)=ZP(11,1)
      ZG(10,2)=   C208*ZP(10,2) + C161*ZP(12,2)
     &          - C263*ZP(13,3) + C161*ZP(14,4)
     &          - C208*ZP(16,4)
      ZG(12,2)=   C161*ZP(10,2) + C125*ZP(12,2)
     &          - C204*ZP(13,3) + C125*ZP(14,4)
     &          - C161*ZP(16,4)
      ZG(13,3)= - C263*ZP(10,2) - C204*ZP(12,2)
     &          + C333*ZP(13,3) - C204*ZP(14,4)
     &          + C263*ZP(16,4)
      ZG(14,4)=ZG(12,2)
      ZG(16,4)=-ZG(10,2)
      ZG(13,5)=   C333*ZP(13,5) - C204*ZP(14,6)
     &          + C263*ZP(16,6) - C263*ZP(10,8)
     &          - C204*ZP(12,8)
      ZG(14,6)= - C204*ZP(13,5) + C125*ZP(14,6)
     &          - C161*ZP(16,6) + C161*ZP(10,8)
     &          + C125*ZP(12,8)
      ZG(16,6)=   C263*ZP(13,5) - C161*ZP(14,6)
     &          + C208*ZP(16,6) - C208*ZP(10,8)
     &          - C161*ZP(12,8)
      ZG(10,8)=-ZG(16,6)
      ZG(12,8)=ZG(14,6)
      ZG(10,10)=   C177*ZP(10,10) + C040*ZP(12,10)
     &           + C040*ZP(10,12) + C156*ZP(12,12)
     &           + C208*ZP(13,13) + C156*ZP(14,14)
     &           - C040*ZP(16,14) + C125*ZP(15,15)
     &           - C040*ZP(14,16) + C177*ZP(16,16)
      ZG(12,10)=   C040*ZP(10,10) + C156*ZP(12,10)
     &           + C156*ZP(10,12) - C040*ZP(12,12)
     &           + C161*ZP(13,13) - C040*ZP(14,14)
     &           - C156*ZP(16,14) - C161*ZP(15,15)
     &           - C156*ZP(14,16) + C040*ZP(16,16)
      ZG(11,11)=ZP(11,11)
      ZG(12,12)=   C156*ZP(10,10) - C040*ZP(12,10)
     &           - C040*ZP(10,12) + C177*ZP(12,12)
     &           + C125*ZP(13,13) + C177*ZP(14,14)
     &           + C040*ZP(16,14) + C208*ZP(15,15)
     &           + C040*ZP(14,16) + C156*ZP(16,16)
      ZG(13,13)=   C208*ZP(10,10) + C161*ZP(12,10)
     &           + C161*ZP(10,12) + C125*ZP(12,12)
     &           + C333*ZP(13,13) + C125*ZP(14,14)
     &           - C161*ZP(16,14) - C161*ZP(14,16)
     &           + C208*ZP(16,16)
      ZG(14,14)=ZG(12,12)
      ZG(16,14)=-ZG(12,10)
      ZG(15,15)=   C125*ZP(10,10) - C161*ZP(12,10)
     &           - C161*ZP(10,12) + C208*ZP(12,12)
     &           + C208*ZP(14,14) + C161*ZP(16,14)
     &           + C333*ZP(15,15) + C161*ZP(14,16)
     &           + C125*ZP(16,16)
      ZG(16,16)=ZG(10,10)
       GO TO 270
C
C                          TETRAGONAL SYMMETRY
C
204   ZG(1,1)=ZP(1,1)
      ZG(7,1)=ZP(7,1)
      ZG(2,2)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(3,3)=ZP(3,3)
      ZG(4,4)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(5,5)=ZP(5,5)
      ZG(6,6)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(7,7)=ZP(7,7)
      ZG(8,8)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(9,9)=ZP(9,9)
       IF(NL.EQ.3) GO TO 270
      ZG(10,2)=RCH*(ZP(10,2)-ZP(16,4))
      ZG(12,2)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(13,3)=ZP(13,3)
      ZG(14,4)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(16,4)=RCH*(-ZP(10,2)+ZP(16,4))
      ZG(10,10)=RCH*(ZP(10,10)+ZP(16,16))
      ZG(12,10)=RCH*(ZP(12,10)-ZP(14,16))
      ZG(11,11)=ZP(11,11)
      ZG(12,12)=RCH*(ZP(12,12)+ZP(14,14))
      ZG(13,13)=ZP(13,13)
      ZG(14,14)=RCH*(ZP(12,12)+ZP(14,14))
      ZG(16,14)=RCH*(-ZP(10,12)+ZP(16,14))
      ZG(15,15)=ZP(15,15)
      ZG(16,16)=RCH*(ZP(10,10)+ZP(16,16))
       GO TO 270
C
C                          ORTHORHOMBIC SYMMETRY
C
205   ZG(1,1)=ZP(1,1)
      ZG(7,1)=ZP(7,1)
      ZG(9,1)=ZP(9,1)
      ZG(2,2)=ZP(2,2)
      ZG(3,3)=ZP(3,3)
      ZG(4,4)=ZP(4,4)
      ZG(5,5)=ZP(5,5)
      ZG(6,6)=ZP(6,6)
      ZG(7,7)=ZP(7,7)
      ZG(9,7)=ZP(9,7)
      ZG(8,8)=ZP(8,8)
      ZG(9,9)=ZP(9,9)
       IF(NL.EQ.3) GO TO 270
      ZG(10,2)=ZP(10,2)
      ZG(12,2)=ZP(12,2)
      ZG(13,3)=ZP(13,3)
      ZG(15,3)=ZP(15,3)
      ZG(14,4)=ZP(14,4)
      ZG(16,4)=ZP(16,4)
      ZG(10,10)=ZP(10,10)
      ZG(12,10)=ZP(12,10)
      ZG(11,11)=ZP(11,11)
      ZG(12,12)=ZP(12,12)
      ZG(13,13)=ZP(13,13)
      ZG(15,13)=ZP(15,13)
      ZG(14,14)=ZP(14,14)
      ZG(16,14)=ZP(16,14)
      ZG(15,15)=ZP(15,15)
      ZG(16,16)=ZP(16,16)
       GO TO 270
C
270   DO 370 I=1,NLSQ-1
      DO 371 J=I+1,NLSQ
      ZG(I,J)=ZG(J,I)
371   CONTINUE
370   CONTINUE
      DO 374 J=1,NLSQ
      DO 375 I=1,NLSQ
      ZA(I,J)=ZG(I,J)
375   CONTINUE
374   CONTINUE
C
      RETURN
      END
C*******************
CXXX    XDOS    ****
C*******************
      SUBROUTINE XDOS
C
C************************************************
C   CALCULATES LOCAL DOS BY MEANS OF
C   ANALYTIC CONTINUATION OF CONDITIONALLY
C   AVERAGED PHYSICAL GREENS FUNCTIONS
C************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZA=300)
C
      DIMENSION ZWG(MNE),EW(MNE),
     &          DOSP(MNE,MNL,MNS,MNA),
     &          DOST(MNE,MNS,MNA)
C
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/GDOS/ ZGDOS(MNL,MNA,MNE,MNS)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      PI=RC4*ATAN(RC1)
C
      DO 3031 IA=1,MNA
      DO 3032 IS=1,MNS
      DO 3033 IE=1,MNE
      DOST(IE,IS,IA)=RCZ
3033  CONTINUE
3032  CONTINUE
3031  CONTINUE
C
!SMP$ DO SERIAL
      DO 304 IE=1,NE
      EW(IE)=DREAL(ZCN(IE))
304   CONTINUE
C
C------------------------------- PARTIAL DOS BY
C                             ANALYTIC CONTINUATION
      DO 320 IA=1,NA
      DO 321 IS=1,NS
      DO 322 IL=1,NL
C
!SMP$ DO SERIAL
      DO 324 IE=1,NE
      ZWG(IE)=ZGDOS(IL,IA,IE,IS)
324   CONTINUE
C
      CALL DECO(ZWG,NE,NDC)
C
      DO 325 IE=1,NE
      DOSP(IE,IL,IS,IA) = - DIMAG(ZWG(IE))/PI
325   CONTINUE
C
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
C-------------------------------  TOTAL DOS
C                            (ATOM- AND SPIN-RESOLVED)
      DO 326 IA=1,NA
      DO 327 IS=1,NS
      DO 328 IL=1,NL
      DO 329 IE=1,NE
      DOST(IE,IS,IA)=DOST(IE,IS,IA)+DOSP(IE,IL,IS,IA)
329   CONTINUE
328   CONTINUE
327   CONTINUE
326   CONTINUE
C
C-------------------------------- PRINT OF DOS
101   FORMAT(1X,10I5)
106   FORMAT(6F12.5)
C.......................... ORBITAL-RESOLVED DOS
C                                   - UNIT IW7
      WRITE(IW7,130)
130   FORMAT(2X,'------ ORBITAL-RESOLVED DOS ------')
      WRITE(IW7,101) NE,NL,NS,NA
      DO 330 IA=1,NA
      DO 331 IS=1,NS
      WRITE(IW7,116) OTXTA(IA),IS
116   FORMAT(4X,A16,4X,'IS=',I1)
      DO 332 IE=1,NE
      WRITE(IW7,106) EW(IE),DOST(IE,IS,IA),
     &               (DOSP(IE,IL,IS,IA),IL=1,NL)
332   CONTINUE
331   CONTINUE
330   CONTINUE
      WRITE(IW7,135)
135   FORMAT(2X,'---------------------------------')
C
C........................ ATOM AND SPIN RESOLVED DOS
C                          - UNITS IW8 AND IW9
      DO 341 IS=1,NS
      IF(IS.EQ.1) IWW=IW8
      IF(IS.EQ.2) IWW=IW9
      WRITE(IWW,140) IS
140   FORMAT(2X,'------  ATOM-RESOLVED DOS FOR SPIN IS=',
     &     I1,'  ------')
      WRITE(IWW,101) NE,NA
      DO 342 IE=1,NE
      WRITE(IWW,106) EW(IE),(DOST(IE,IS,IA),IA=1,NA)
342   CONTINUE
      WRITE(IWW,145)
145   FORMAT(2X,'---------------------------------')
341   CONTINUE
C
      RETURN
      END
C*******************
CXXX    DECO    ****
C*******************
      SUBROUTINE DECO(ZG,NE,NDC)
C
C****************************************************
C      DECONVOLUTION (ANALYTIC CONTINUATION)
C      OF A GREEN'S FUNCTION FROM THE UPPER
C      COMPLEX HALF-PLANE TO THE REAL ENERGY AXIS
C****************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNE=1000)
C
      DIMENSION ZG(MNE),ZW(MNE)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RCH/0.5D0/
C
      ZZ=DCMPLX(RCZ,RCZ)
C
C                                QUADRATIC CONTINUATION
C      ZM=DCMPLX(-RCH,RCH)
C      AA=RC2
C      ZP=DCMPLX(-RCH,-RCH)
C                                LINEAR CONTINUATION
      ZM=DCMPLX(RCZ,RCH)
      AA=RC1
      ZP=DCMPLX(RCZ,-RCH)
C
      DO 310 IDC=1,NDC
!SMP$ DO SERIAL      
      DO 312 IE=1,NE
      ZW(IE)=ZG(IE)      
312   CONTINUE
      ZG(IDC)=ZZ
      DO 316 IE=IDC+1,NE-IDC
      ZG(IE) = ZM*ZW(IE-1) + AA*ZW(IE) + ZP*ZW(IE+1)
316   CONTINUE
      ZG(NE-IDC+1)=ZZ
310   CONTINUE
C
      RETURN
      END
C*******************
CXXX    ATBK    ****
C*******************
      SUBROUTINE ATBK
C
C***********************************************
C   BLOCH TRANSFORM OF TB-STRUCTURE CONSTANTS
C    ( FOR THE BLOCH SPECTRAL FUNCTIONS )
C***********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNL=3)
      PARAMETER(MNBZA=300)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=32)
C
      DIMENSION TRAN(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/KMESHA/ AKBZA(3,MNBZA),NBZA
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB),JBVA(MNCL,MNB),
     &            JTRA(3,MNCL,MNB),NSCL(MNB)
      COMMON/SCKA/ ZSKA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
C       WRITE(*,*) AKBZA(1:3,1)
C       write(*,*) VBR(1:3,1)
C       write(*,*) VBR(1:3,2)
C       write(*,*) VBR(1:3,3)
C       write(*,*)'------'
      
      WRITE(IW6,111)
111   FORMAT(//3X,' **********  BLOCH TRANSFORM ',
     &   'OF TB-CONSTANTS  ********** '/
     &    14X,'( FOR THE BLOCH SPECTRAL FUNCTIONS )'/
     &   /10X,'    CHECK OF HERMITICITY : ')
C
      DO 301 IBZA=1,MNBZA
      DO 302 II=1,MNBLSQ
      DO 303 JJ=1,MNBLSQ
      ZSKA(JJ,II,IBZA)=ZZ
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
      DO 331 IBZA=1,NBZA
      DO 332 IB=1,NB
      II0=(IB-1)*NLSQ
C
      DO 333 JCL=1,NSCL(IB)
      JB=JBVA(JCL,IB)
      JTRA1=JTRA(1,JCL,IB)
      JTRA2=JTRA(2,JCL,IB)
      JTRA3=JTRA(3,JCL,IB)
C
      JJ0=(JB-1)*NLSQ
      
      DO 334 I=1,3
      TRAN(I)=JTRA1*VBR(I,1)+JTRA2*VBR(I,2)+JTRA3*VBR(I,3)
334   CONTINUE
C       write(*,'(6g11.4)') tran(1:3)
      ARG=RCZ
C       write(*,'(6g11.4)') AKBZA(1:3,IBZA)*TRAN(1:3)
C        write(*,'(4g12.5)') TRAN(1:3)*AKBZA(1:3,IBZA),sum(TRAN(1:3)
C      * *AKBZA(1:3,IBZA))
C       write(*,*)'------'
      DO 335 I=1,3
      ARG=ARG+AKBZA(I,IBZA)*TRAN(I)
335   CONTINUE
C       write(*,*) 'F:',ARG
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      DO 336 IQ=1,NLSQ
      DO 337 JQ=1,NLSQ
      ZSKA(JJ0+JQ,II0+IQ,IBZA) = ZSKA(JJ0+JQ,II0+IQ,IBZA)
     &                         + ZEFA * STR(JQ,IQ,JCL,IB)
337   CONTINUE
336   CONTINUE
C
333   CONTINUE
332   CONTINUE
331   CONTINUE
C
C                                   CHECK OF HERMITICITY
      DEV=RCZ
      DO 341 IBZA=1,NBZA
      DO 342 II=1,NBLSQ
      DO 343 JJ=II,NBLSQ
      ZW1=ZSKA(II,JJ,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKA(JJ,II,IBZA)
      W2X=DREAL(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
343   CONTINUE
342   CONTINUE
341   CONTINUE
      WRITE(IW6,140) DEV
140   FORMAT(/7X,'            ------ DEVIATION=',G12.4)
C
      RETURN
      END
C*******************
CXXX    GBSF    ****
C*******************
      SUBROUTINE GBSF(IE,IS)
C
C********************************************
C   CALCULATES THE GREEN FUNCTION QUANTITIES
C   FOR THE BLOCH SPECTRAL FUNCTIONS
C********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZA=300)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZWD(MNLSQ,MNLSQ),ZWF(MNLSQ,MNLSQ),
     &          ZWG(MNLSQ,MNLSQ),ZWA(MNBLSQ,MNBLSQ),
     &          ZCMU(MNLSQ,MNLSQ,MNB),
     &          ZCLA(MNLSQ,MNLSQ,MNB)
C
      COMMON/KMESHA/ AKBZA(3,MNBZA),NBZA
      COMMON/SCKA/ ZSKA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNS)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/GBLO/ ZGBLO(MNL,MNB,MNBZA,MNE,MNS)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
      Z1=DCMPLX(RC1,RCZ)
C
      DO 3003 IB=1,MNB
      DO 3004 J=1,MNLSQ
      DO 3005 I=1,MNLSQ
      ZCMU(I,J,IB)=ZZ
3005  CONTINUE
3004  CONTINUE
3003  CONTINUE
      DO 3013 IB=1,MNB
      DO 3014 J=1,MNLSQ
      DO 3015 I=1,MNLSQ
      ZCLA(I,J,IB)=ZZ
3015  CONTINUE
3014  CONTINUE
3013  CONTINUE
C
C---------------------- CALIGRAPHIC QUANTITIES MU AND LAMBDA
C
C.................................. SITE-LOOP
      IA=0
      DO 3106 IB=1,NB
C
      DO 3120 J=1,NLSQ
      DO 3121 I=1,NLSQ
      ZWG(I,J)=ZGFG(I,J,IB,IS)
3121  CONTINUE
3120  CONTINUE
C.................................. SPECIES-LOOP
      DO 3124 IC=1,NC(IB)
      IA=IA+1
C                                    MATRIX F
      DO 3127 J=1,NLSQ
      DO 3128 I=1,NLSQ
      ZWD(I,J)=-ZCPF(I,J,IB,IS)
3128  CONTINUE
3127  CONTINUE
      DO 3129 I=1,NLSQ
      ZWD(I,I)=ZWD(I,I)+ZPF(I,IA,IS)
3129  CONTINUE
C
      CALL MAMU(ZWD,ZWG,ZWF,NLSQ,MNLSQ)
C
      DO 3132 I=1,NLSQ
      ZWF(I,I)=ZWF(I,I)+Z1
3132  CONTINUE
C
      CALL CINVC(ZWF,NLSQ,MNLSQ)
C                                  CONTRIBUTION TO
C                                  CALIGRAPHIC MU
      DO 3136 J=1,NLSQ
      DO 3137 I=1,NLSQ
      ZCMU(I,J,IB) = ZCMU(I,J,IB)
     &        +  CON(IA) * ZWF(I,J) * ZMU(J,IA,IS)
3137  CONTINUE
3136  CONTINUE
C                                  CONTRIBUTION TO
C                                  CALIGRAPHIC LAMBDA
      CALL MAMU(ZWG,ZWF,ZWD,NLSQ,MNLSQ)
C
      DO 3140 I=1,NLSQ
      DO 3141 J=1,NLSQ
      ZWD(I,J)=ZMU(I,IA,IS)*ZWD(I,J)
3141  CONTINUE
3140  CONTINUE
      DO 3144 J=1,NLSQ
      DO 3145 I=1,NLSQ
      ZWD(I,J)=ZWD(I,J)*ZMU(J,IA,IS)
3145  CONTINUE
3144  CONTINUE
      DO 3148 I=1,NLSQ
      ZWD(I,I)=ZWD(I,I)+ZLA(I,IA,IS)
3148  CONTINUE
C
      DO 3152 J=1,NLSQ
      DO 3153 I=1,NLSQ
      ZCLA(I,J,IB) = ZCLA(I,J,IB)
     &             + CON(IA) * ZWD(I,J)
3153  CONTINUE
3152  CONTINUE
C.............................. END OF SPECIES-LOOP
3124  CONTINUE
C                                  CALIGRAPHIC LAMBDA
      DO 3160 J=1,NLSQ
      DO 3161 I=1,NLSQ
      ZWF(I,J)=ZCMU(I,J,IB)
3161  CONTINUE
3160  CONTINUE
C
      CALL MAMU(ZWG,ZWF,ZWD,NLSQ,MNLSQ)
C
      DO 3164 J=1,NLSQ
      DO 3165 I=1,NLSQ
      ZWF(J,I)=ZCMU(I,J,IB)
3165  CONTINUE
3164  CONTINUE
C
      CALL MAMU(ZWF,ZWD,ZWG,NLSQ,MNLSQ)
C
      DO 3168 J=1,NLSQ
      DO 3169 I=1,NLSQ
      ZCLA(I,J,IB)=ZCLA(I,J,IB)-ZWG(I,J)
3169  CONTINUE
3168  CONTINUE
C.............................. END OF SITE-LOOP
3106  CONTINUE
C
C------------------------------ GREEN FUNCTION QUANTITIES
C
C................................. IBZA-LOOP
C       do i=1,9
C       write(*,*) dreal(ZSKA(i,i,1)),dimag(ZSKA(i,i,1))
C       end do

      DO 3202 IBZA=1,NBZA
C                               MATRIX P(Z)-S(K)
      DO 3214 JJ=1,NBLSQ
      DO 3215 II=1,NBLSQ
      ZWA(II,JJ)=-ZSKA(II,JJ,IBZA)
3215  CONTINUE
3214  CONTINUE
      DO 3216 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3217 J=1,NLSQ
      DO 3218 I=1,NLSQ
      ZWA(II0+I,II0+J)=ZWA(II0+I,II0+J)+ZCPF(I,J,IB,IS)
3218  CONTINUE
3217  CONTINUE
3216  CONTINUE

C       write(*,*) dreal(ZCPF(1:9,1:9,1,IS))
C       write(*,*) dimag(ZCPF(1:9,1:9,1,IS))

C                                      INVERSION
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
C       write(*,*) 'h:',-dimag(ZWA(3,3))
C
C.................................... SITE-LOOP
      DO 3221 IB=1,NB
C
      II0=(IB-1)*NLSQ
      DO 3224 J=1,NLSQ
      DO 3225 I=1,NLSQ
      ZWG(I,J)=ZWA(II0+I,II0+J)
3225  CONTINUE
3224  CONTINUE
C
      DO 3228 J=1,NLSQ
      DO 3229 I=1,NLSQ
      ZWF(J,I)=ZCMU(I,J,IB)
3229  CONTINUE
3228  CONTINUE
C
      CALL MAMU(ZWF,ZWG,ZWD,NLSQ,MNLSQ)
C
      DO 3242 J=1,NLSQ
      DO 3243 I=1,NLSQ
      ZWF(I,J)=ZCMU(I,J,IB)
3243  CONTINUE
3242  CONTINUE
C
      CALL MAMU(ZWD,ZWF,ZWG,NLSQ,MNLSQ)
C
      DO 325 IL=1,NL
      ZSUM=ZZ
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 3257 I=ISTA,IFIN
      ZSUM = ZSUM + ZWG(I,I) + ZCLA(I,I,IB)
3257  CONTINUE
      ZGBLO(IL,IB,IBZA,IE,IS)=ZSUM
325   CONTINUE
C::::::::::::::::::::::::::::::::::: END OF THE LOOPS
3221  CONTINUE
3202  CONTINUE
C
      RETURN
      END
C*******************
CXXX    XBSF    ****
C*******************
      SUBROUTINE XBSF
C
C********************************************
C   CALCULATES THE BLOCH SPECTRAL FUNCTIONS
C   BY ANALYTIC CONTINUATION
C********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=2)
      PARAMETER(MNA=4)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=1000)
      PARAMETER(MNBZA=300)
C
      DIMENSION ZWX(MNE),EW(MNE),
     &          ABLOP(MNE,MNL,MNS,MNB,MNBZA),
     &          ABLOT(MNE,MNS,MNB,MNBZA)
C
      COMMON/KMESHA/ AKBZA(3,MNBZA),NBZA
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/GBLO/ ZGBLO(MNL,MNB,MNBZA,MNE,MNS)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      PI=RC4*ATAN(RC1)
C
      DO 3041 IBZA=1,MNBZA
      DO 3042 IB=1,MNB
      DO 3043 IS=1,MNS
      DO 3044 IE=1,MNE
      ABLOT(IE,IS,IB,IBZA)=RCZ
3044  CONTINUE
3043  CONTINUE
3042  CONTINUE
3041  CONTINUE
C
!SMP$ DO SERIAL
      DO 308 IE=1,NE
      EW(IE)=DREAL(ZCN(IE))
308   CONTINUE
      
C       DO  IS=1,NS
C       DO  IE=1,NE
C        write(*,*) -dimag(sum(ZGBLO(1:NL,1,1,IE,IS)))/PI
C       end do
C       end do
C
C------------------------------ BLOCH SPECTRAL FUNCTIONS
C
C                                      PARTIAL
      DO 3302 IBZA=1,NBZA
      DO 3303 IB=1,NB
      DO 3304 IS=1,NS
      DO 3305 IL=1,NL
C
!SMP$ DO SERIAL
      DO 3314 IE=1,NE
      ZWX(IE)=ZGBLO(IL,IB,IBZA,IE,IS)
3314  CONTINUE
C
      CALL DECO(ZWX,NE,NDC)
C
      DO 3316 IE=1,NE
      ABLOP(IE,IL,IS,IB,IBZA) = - DIMAG(ZWX(IE))/PI
3316  CONTINUE
C
3305  CONTINUE
3304  CONTINUE
3303  CONTINUE
3302  CONTINUE
C
C                                       TOTAL
      DO 3321 IBZA=1,NBZA
      DO 3322 IB=1,NB
      DO 3323 IS=1,NS
      DO 3324 IL=1,NL
      DO 3325 IE=1,NE
      ABLOT(IE,IS,IB,IBZA) = ABLOT(IE,IS,IB,IBZA)
     &                     + ABLOP(IE,IL,IS,IB,IBZA)
3325  CONTINUE
3324  CONTINUE
3323  CONTINUE
3322  CONTINUE
3321  CONTINUE
C
C-------------------------------- PRINT OF BSF
101   FORMAT(1X,10I5)
106   FORMAT(6F12.5)
C.......................... ORBITAL-RESOLVED BSF
C                                   - UNIT IW17
      WRITE(IW17,130)
130   FORMAT(2X,'------ ORBITAL-RESOLVED BSF ------')
      WRITE(IW17,101) NE,NL,NS,NB,NBZA
      DO 340 IBZA=1,NBZA
      DO 341 IB=1,NB
      DO 342 IS=1,NS
      WRITE(IW17,132) OTXTK(IBZA),IBZA,IB,IS
132   FORMAT(4X,A16,4X,'(IBZA=',I2,')',
     &              4X,'IB=',I2,4X,'IS=',I1)
      DO 343 IE=1,NE
      WRITE(IW17,106) EW(IE),ABLOT(IE,IS,IB,IBZA),
     &               (ABLOP(IE,IL,IS,IB,IBZA),IL=1,NL)
343   CONTINUE
342   CONTINUE
341   CONTINUE
340   CONTINUE
      WRITE(IW17,135)
135   FORMAT(2X,'---------------------------------')
C
C........................ SITE AND SPIN RESOLVED BSF
C                          - UNITS IW18 AND IW19
      DO 351 IS=1,NS
      IF(IS.EQ.1) IWW=IW18
      IF(IS.EQ.2) IWW=IW19
      WRITE(IWW,140) IS
140   FORMAT(2X,'------  SITE-RESOLVED BSF FOR SPIN IS=',
     &     I1,'  ------')
      WRITE(IWW,101) NE,NB,NBZA
      DO 352 IBZA=1,NBZA
      WRITE(IWW,142) OTXTK(IBZA),IBZA,IS
142   FORMAT(4X,A16,4X,'(IBZA=',I2,')',4X,'IS=',I1)
      DO 353 IE=1,NE
      WRITE(IWW,106) EW(IE),(ABLOT(IE,IS,IB,IBZA),IB=1,NB)
353   CONTINUE
      IE=NE/2+1
      WRITE(354+IS,106) EW(IE),(ABLOT(IE,IS,IB,IBZA),IB=1,NB)
352   CONTINUE
      WRITE(IWW,145)
145   FORMAT(2X,'---------------------------------')
351   CONTINUE
C
      RETURN
      END
C*******************
CXXX    MAMU    ****
C*******************
      SUBROUTINE MAMU(ZA,ZB,ZC,N,ND)
C
C----------------------------------------------------
C   MULTIPLICATION OF COMPLEX MATRICES: ZC = ZA * ZB
C----------------------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      DIMENSION ZA(ND,N),ZB(ND,N),ZC(ND,N)
C
      DATA RCZ/0.0D0/
C
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 J=1,N
      DO 302 I=1,N
      ZC(I,J)=ZZ
302   CONTINUE
      DO 303 K=1,N
      DO 304 I=1,N
      ZC(I,J)=ZC(I,J)+ZA(I,K)*ZB(K,J)
304   CONTINUE
303   CONTINUE
301   CONTINUE
C
      RETURN
      END

