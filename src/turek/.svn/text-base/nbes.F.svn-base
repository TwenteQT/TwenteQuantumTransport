C*******************
CXXX    BES     ****
C*******************
      PROGRAM BES
C
C*******************************************************
C      PROGRAM FOR SELFCONSISTENT ELECTRONIC
C        STRUCTURE CALCULATION OF
C        BULK ORDERED AND DISORDERED ALLOYS
C       BY MEANS OF THE TB-LMTO-CPA METHOD
C*******************************************************
C      VERSION  OCTOBER 1999
C*******************************************************
C   BASIC FEATURES OF THE PROGRAM ARE:
C     - ELECTRONIC STRUCTURE OF ORDERED AND DISORDERED
C       BULK SYSTEMS
C     - CASES: CUBIC, HEXAGONAL, ...
C     - GEOMETRY: QUITE GENERAL
C       (WITH 3D-TRANSLATIONAL SYMMETRY)
C     - NON-RELATIVISTIC AND SCALAR-RELATIVISTIC VERSION
C       (OPTIONALLY WITH FULLY RELATIVISTIC CORES)
C     - NON-MAGNETIC AND SPIN-POLARIZED
C       (WITH COLLINEAR SPIN STRUCTURES ONLY)
C     - SUBSTITUTIONAL CHEMICAL DISORDER
C     - ATOMIC SPHERE APPROXIMATION
C     - SPD AND SPDF CASES
C     - EVALUATION OF TOTAL ENERGIES
C*******************************************************
C   INPUT:   IR1  (= 1)  -  GEOMETRY OF THE PROBLEM
C            IR2  (= 2)  -  CHEMICAL OCCUPATIONS
C                           OF THE SITES
C            IR4  (= 4)  -  COHERENT POTENTIAL FUNCTIONS
C            IR5  (= 5)  -  CONTROL PARAMETERS
C            IR12 (=12)  -  NORMALIZED COMPLEX NODES AND
C                           WEIGHTS FOR CONTOUR INTEGRATION
C   OUTPUT:  IW6  (= 6)  -  CURRENT OUTPUT
C            IW7  (= 7)  -  LSDA FILE - SELFCONSISTENT
C                           POT. PARAMETERS AND POTENTIALS
C            IW8  (= 8)  -  COHERENT POTENTIAL FUNCTIONS
C                           (SELFCONSISTENT)
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
C    NE ...... NUMBER OF NODES ON THE COMPLEX
C              ENERGY CONTOUR
C   MNE ...... PARAMETER: MAX. NE
C    NBZ ..... NUMBER OF K-POINTS IN THE IBZ
C   MNBZ ..... PARAMETER: MAX. NBZ
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
C     MNR .... PARAMETER: MAX. SIZE OF RADIAL MESH
C   MNAM ..... PARAMETER: MAX. DIMENSION OF VECTORS FOR
C              ANDERSON MIXING (MNAM=MNA*MNS*MNR)
C   MNUH ..... PARAMETER: MAX. NUMBER OF PREVIOUS
C              VECTORS FOR ANDERSON MIXING
C******************************************************
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=30)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
C*******************************************************
C             ALL COMMON BLOCKS  :
C
C-----------------------------------------  GEOMETRY
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB),JBVA(MNCL,MNB),
     &            JTRA(3,MNCL,MNB),NSCL(MNB)
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/AMAD/ AMCMM(MNB,MNB)
C
C-------------------------------------  GREEN'S FUNCTIONS
C
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CAGF/ ZCAGF(MNLSQ,MNLSQ,MNA,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
C
C-------------------------------------  LSDA
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/WAF/ PHI(MNR,MNL,MNS,MNA),
     &            PHID(MNR,MNL,MNS,MNA),
     &            PHIDD(MNR,MNL,MNS,MNA)
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/LIN/ ENY1,PPC1,PPD1,PPQ1,PPP1,DNY1,FINY1,FINYD1,
     &            FI(MNR),FID(MNR),FIDD(MNR)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA)
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),
     &            VMAD(MNA),ACTE(MNA)
C
C-------------------------------------  GENERAL
C
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
C
C-------------------------------------  SPECIAL
C
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/ITEF/ DIFME,QMIXE
      COMMON/AMXP/ DXP(MNAM,MNUH)
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
      COMMON/TEXT/ OTXTA(MNA)
      common/bshift/ ishift,iashift,ilshift,eshift,
     &               ivshift,iavshift,vshift
      integer::iconv
      integer,external::cpait
C
C********************************************************
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
C                       BLOCH TRANSFORM OF TB-CONSTANTS
      CALL TBK
      flush(IW6)      
C                              MADELUNG CONSTANTS
      CALL MACO
      flush(IW6)      
C                          START OF LDA ITERATIONS
      CALL SLDA
      flush(IW6)

C
C----------------------------------  ITERATION LOOP
      DO 310 IITER=1,NITER
C                            COHERENT POTENTIAL FUNCTION
      CALL SCPF
      flush(IW6)      
C                            ON-SITE GREEN FUNCTIONS
      CALL ONGF
      flush(IW6)      
C                               CPA STEP
      ICONV=CPAIT(1)
C      CALL CPAIT
C                             ENERGY MOMENTS
      CALL CMOM
      flush(IW6)      
C                              CORE DENSITIES
      CALL CODE
      flush(IW6)      
C                              VALENCE DENSITIES
      CALL VADE
      flush(IW6)      
C                              NEW POTENTIALS
      CALL NEPO
      flush(IW6)      
C                              PRINT OF RESULTS
      CALL TISK
      flush(IW6)      
C------------------------------  END OF ITERATION LOOP
310   CONTINUE

      if (ICONV/=0) then
         write(IW6,*) ' No convergence!!! '
         flush(IW6)         
C         stop
      end if      
      call write_atoms(NA,NL,NS,NSZRAD,EF,AZ,WS,WSAV,
     & QSCR,POT,OTXTA,AMGMOM)
C
      WRITE(IW6,191)
      flush(IW6)
191   FORMAT(/'        * * *   END OF BES   * * * '/)
C
      STOP
      END
C  
      subroutine write_atoms(NA,NL,NS,NSZRAD,EF,AZ,WS,WSAV,
     & QSCR,POT,OTXTA,AMGMOM)
      implicit none
      integer,PARAMETER::MNA=4
      integer,PARAMETER::MNL=4
      integer,PARAMETER::MNS=2
      integer,PARAMETER::MNR=400
      integer,PARAMETER::MNBZ=160002
      
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      
      integer:: LCOR,NCOR,NUMCOR,NOBC,IREL,IVXC,ISHENY
      integer:: NSYM,NK,INVE,NBZ
      real(kind=8)::ECOR,ENY,PPC,PPD,PPP,FINY,DNY,PPQ,FINYD
      real(kind=8)::THRESH,AKBZ,WKBZ
      Real (Kind=8) :: AMGMOM(NA)
      
      integer::NA,NL,NS
      integer::NSZRAD(NA)
      REAL(kind=8)::AZ(NA),WS(NA),QSCR(NL),EF,WSAV(NA)
      REAL(kind=8)::POT(MNR,MNS,MNA)
      CHARACTER*16::OTXTA(NA)
      
      integer::IA,IS,IL,I,J,IX
      integer,parameter::atomf=999
      character(len=16):: cwork

992   FORMAT(1X,I5,   54X,'# ',A10,' IS=',I1)
993   FORMAT(1X,I5,   42X,A16)      
994   FORMAT(1X,4G15.7)      
995   FORMAT(1X,G15.7,   42X,A16)
99511 FORMAT(1X,G15.7,G15.7,29X,A20)
996   FORMAT(1X,E15.7,   39X,A16)
997   FORMAT(1X,3G15.7,7X,A16)

1101   FORMAT(1X,10I5)
1104   FORMAT(1X,4G15.7)
C
1190   FORMAT(1X,'  --------  LSDA-FILE  -------- ')
1191   FORMAT(1X,49X,A16,' IS=',I1)
1192   FORMAT(1X,I5,   44X,A16,' IS=',I1)
1193   FORMAT(1X,I5,   44X,A16)
1195   FORMAT(1X,'  ----------------------------- ')
      
      do IA=1,NA
         read(OTXTA(IA),*) cwork
         open (unit=atomf, file='atoms/'//trim(cwork)//'_B',
     &    action='write')
         write(atomf,'(A,A,A,I1,A,I5)') '#  Bulk results for: ',
     &   trim(cwork),', IVXC=',IVXC,', NK=',NK
         write(atomf,993) NL, '# NL'
         write(atomf,994) (QSCR(IL),IL=1,NL)
         WRITE(atomf,997) AZ(IA),WS(IA),WSAV(IA),'# AZ, WSR'
         write(atomf,99511) EF, AMGMOM(IA), '# Fermi.En, Mag.Mom'         
         write(atomf,996) 0.0d0,'# Pot.Shift'
         write(atomf,993) 0,'# SW'
         write(atomf,993) NS,'# NS'
         DO IS=1,NS
            WRITE(atomf,992) NSZRAD(IA),OTXTA(IA),IS
            WRITE(atomf,994) (POT(I,IS,IA),I=1,NSZRAD(IA))
         end do
         write(atomf,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'         
         DO IS=1,NS
            WRITE(atomf,1191) OTXTA(IA),IS
            WRITE(atomf,1104) (ENY(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (PPC(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (PPD(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (PPQ(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (PPP(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (DNY(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (FINY(IL,IS,IA),IL=1,NL)
            WRITE(atomf,1104) (FINYD(IL,IS,IA),IL=1,NL)
         end do
         WRITE(atomf,1193) NUMCOR(IA),OTXTA(IA)
         IF(NUMCOR(IA).NE.0) then
            DO J=1,NUMCOR(IA)
               WRITE(atomf,1101) NCOR(J,IA),LCOR(J,IA),NOBC(J,IA)
               WRITE(atomf,1104) (ECOR(J,IX,IA),IX=1,2)
            end do
         end if
         close (unit=atomf)
      end do      
      return
      end
      
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
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
C
C                     INPUT UNITS
      IR1=1
      IR2=2
      IR4=4
      IR5=5
      IR11=11
      IR12=12
C                    OUTPUT UNITS
      IW6=6
      IW7=7
      IW8=8
      IW9=9
C
      OPEN(UNIT=IR1,FILE='inpge',FORM='FORMATTED')
      OPEN(UNIT=IR2,FILE='inpch',FORM='FORMATTED')
!      OPEN(UNIT=IR4,FILE='inpcp',FORM='FORMATTED')
      OPEN(UNIT=IR5,FILE='input',FORM='FORMATTED')
      OPEN(UNIT=IR12,FILE='cci.nw',FORM='FORMATTED')
      OPEN(UNIT=IW6,FILE='outit',FORM='FORMATTED')
C      OPEN(UNIT=IW7,FILE='outld',FORM='FORMATTED')
      OPEN(UNIT=IW8,FILE='outcp',FORM='FORMATTED')
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
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNR=400)
      PARAMETER(MNUH=48)
      integer,parameter::atomf=999      
C
      DIMENSION OWORK(5),ORELA(2,0:2),OXCP(2,0:2),SCX(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/ITEF/ DIFME,QMIXE
      COMMON/TEXT/ OTXTA(MNA)
      DIMENSION OTXTA1(MNA)
      common/bshift/ ishift,iashift,ilshift,eshift,
     &               ivshift,iavshift,vshift
   
      character(len=16)::cwork
C
      DATA ORELA/'      NONE      ','                ',
     &           '      SCALAR    ','                ',
     &           '      SCALAR    ','( CORES - FULL )'/
      DATA OXCP/'    VON BARTH - ','HEDIN           ',
     &          '  CEPERLEY-ALDER',' (PERDEW-ZUNGER)',
     &          '    VOSKO - WILK',' - NUSAIR       '/
C
100   FORMAT(5A16)
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
C
      WRITE(IW6,111)
      flush(IW6)
111   FORMAT(/'       ****  OUTPUT OF BES  **** '/
     &     /'     ** READING OF ALL DATA:')
C
C-------------------------  CONTROL DATA - UNIT IR5
C
1001  FORMAT(//5A16)
102   FORMAT(1X,F10.5,I5)
107   FORMAT(1X,7F10.5)
11011 FORMAT(1X,3I5,G15.7)
11012 FORMAT(1X,2I5,G15.7)
C
      READ (IR5,100) OWORK
      WRITE(IW6,1001) OWORK
      READ (IR5,100) OWORK
      WRITE(IW6,1001) OWORK
      flush(IW6)      
C
      READ (IR5,101) NL,NS
      IF(NL.LT.3.OR.NL.GT.4) GO TO 210
      IF(NL.GT.MNL) GO TO 211
      IF(NS.GT.MNS) GO TO 212
      READ (IR5,101) IREL,IVXC
      IF(IREL.LT.0.OR.IREL.GT.2) GO TO 214
      IF(NS.EQ.2.AND.IREL.EQ.2) GO TO 215
      IF(IVXC.LT.0.OR.IVXC.GT.2) GO TO 216
      READ (IR5,101) NSYM,NK,INVE
      READ (IR5,101) NE
      IF(NE.GT.MNE) GO TO 217
      IF(NE.NE.7.AND.NE.NE.10.AND.NE.NE.14.AND.NE.NE.20.AND.
     &  NE.NE.28.AND.NE.NE.40.AND.NE.NE.56.AND.NE.NE.80)
     &  GO TO 218
      READ (IR5,107) DIAM,EF
      READ (IR5,101) ICONT
      IF(ICONT.LT.0.OR.ICONT.GT.1) GO TO 219
      READ (IR5,101) NFUPR
      READ (IR5,107) DIFMC,QMIXC
      READ (IR5,107) DIFME,QMIXE
      READ (IR5,102) ALFA,NITER
      READ (IR5,102) BETA,NITERA
      NITERA=MAX(NITERA,2)
      READ (IR5,102) W0AM,NUH
      NUH=MAX(NUH,2)
      IF(NUH.GT.MNUH) GO TO 220
      read (ir5,11011) ishift,iashift,ilshift,eshift
      read (ir5,11012) ivshift,iavshift,vshift
      if (ishift.ne.0) ishift=1
      if (ivshift.ne.0) ivshift=1
C
      WRITE(IW6,115) NL,NS
115   FORMAT(/17X,'NL=',I2,17X,'NS=',I2)
      WRITE(IW6,1161) (ORELA(I,IREL),I=1,2)
1161  FORMAT(/8X,'RELATIVITY:   ',2A16)
      WRITE(IW6,1162) (OXCP(I,IVXC),I=1,2)
1162  FORMAT(/8X,'XC-POTENTIAL:   ',2A16)
      WRITE(IW6,117) NSYM,NK,INVE
117   FORMAT(/'       LOCAL SYMMETRY:    NSYM=',I2/
     &      /'     K-MESH:   NK=',I3,10X,'INVE=',I2)
      WRITE(IW6,118) NE,DIAM,EF
118   FORMAT(/'  ENERGY CONTOUR:    NE=',I3/
     &    /10X,' DIAM=',F12.5,'        STARTING EF=',F12.5)
      WRITE(IW6,119) ICONT
119   FORMAT(/'    CONTINUATION:  ICONT=',I1)
      WRITE(IW6,120) NFUPR
120   FORMAT(/'          FULL PRINT AFTER',I4,' ITERATIONS')
      WRITE(IW6,121) DIFMC,QMIXC
121   FORMAT(/'  CPA ITERATIONS:  DIFMAX=',G12.4,
     &        8X,'MIXING=',F10.5)
      WRITE(IW6,122) DIFME,QMIXE
122   FORMAT(/'   EF ITERATIONS:  DIFMAX=',G12.4,
     &        8X,'MIXING=',F10.5)
      WRITE(IW6,125) ALFA,NITER
125   FORMAT(/'  LDA ITERATIONS: '/5X,
     & 'STARTING MIXING=',F10.5,'    TOTAL:',I5,' ITERATIONS')
      WRITE(IW6,126) BETA,NITERA
126   FORMAT(5X,'ANDERSON MIXING=',F10.5,
     &           '     AFTER',I5,' ITERATIONS')
      WRITE(IW6,127) W0AM,NUH
127   FORMAT(5X,'QUANTITY W0=',F10.5,4X,
     &           '  PREVIOUS',I5,' ITERATIONS')
      flush(IW6)
      WRITE(IW6,11271) ishift,iashift,ilshift,eshift
11271 FORMAT(5X,'ISHIFT=',I2,4X,
     &  'IASHIFT=',I2,4x,'ILSHIFT=',I2,4x,'ESHIFT=',G15.7)
      flush(IW6)
      WRITE(IW6,11272) ivshift,iavshift,vshift
11272 FORMAT(5X,'IVSHIFT=',I2,4X,
     &  'IAVSHIFT=',I2,4x,'VSHIFT=',G15.7)
      flush(IW6)
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
      DO 315 J=1,3
      DO 316 I=1,3
      VBR(I,J)=SCX(I)*VBR(I,J)
316   CONTINUE
315   CONTINUE
      DO 318 IB=1,NB
      DO 319 I=1,3
      POS(I,IB)=SCX(I)*POS(I,IB)
319   CONTINUE
318   CONTINUE
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
1041  FORMAT(3X,'IB=',I2,5X,3G15.7)
      WRITE(IW6,135)
135   FORMAT(/'       SITES OF THE ATOMIC BASIS: ')
      DO 320 IB=1,NB
      WRITE(IW6,1041) IB,(POS(I,IB),I=1,3)
320   CONTINUE
      flush(IW6)
C
C------------------------- CHEMICAL OCCUPATION - UNIT IR2
C
      READ (IR2,100) OWORK
      WRITE(IW6,1001) OWORK
C
108   FORMAT(A16)
1108   FORMAT(A16,A16)
      IA=0
      DO 340 IB=1,NB
      READ (IR2,101) NC(IB)
      DO 341 IC=1,NC(IB)
      IA=IA+1
      IF(IA.GT.MNA) GO TO 241
      READ (IR2,1108) OTXTA(IA),OTXTA1(IA)
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
      WRITE(IW6,142) IA,OTXTA(IA),OTXTA1(IA),CON(IA),VALZ(IA)
142   FORMAT(' IA=',I2,' - LABEL,INIT,CONC.,VALENCY: ',
     &                       A16,A16,F10.5,F10.3)
345   CONTINUE
344   CONTINUE
      WRITE(IW6,143) NA
      flush(IW6)      
143   FORMAT(/4X,' TOTAL NUMBER OF DIFFERENT ',
     &    'ATOMIC TYPES:  NA=',I3)
C
C ---------------------------- LSDA FILE - UNIT IR3
C
C ----------------------------------------------------
      
      do IA=1,NA
         if (ICONT.EQ.0) then
            read(OTXTA1(IA),*) cwork
            WRITE(IW6,*)' Read defs for:', CWORK
            open (unit=atomf, file='atoms/'//trim(cwork),
     &           action='read')                   
         else
            read(OTXTA(IA),*) cwork
            WRITE(IW6,*)' Read old defs for:', CWORK
            open (unit=atomf, file='atoms/'//trim(cwork)//'_B',
     &           action='read')       
            
         end if            

         do is=1,3
            read(atomf,*) cwork
         end do
         READ (atomf,104) AZ(IA),WS(IA),WSAV(IA)
         if (ICONT.EQ.1) then
            READ (atomf,104) EF
         else
            read(atomf,*) cwork
         end if
      
         do is=1,3
            read(atomf,*) cwork
         end do         
         
         do IS=1,NS
            READ (atomf,101) NSZRAD(IA)
            IF(NSZRAD(IA).GT.MNR) GO TO 245
            READ (atomf,104) (POT(I,IS,IA),I=1,NSZRAD(IA))
         end do         
         read(atomf,*) cwork
         do IS=1,NS
            READ (atomf,104) cwork
            READ (atomf,104) (ENY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPC(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPD(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPQ(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (PPP(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (DNY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (FINY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (FINYD(IL,IS,IA),IL=1,NL)
         end do
         READ (atomf,101) NUMCOR(IA)
         IF(NUMCOR(IA).NE.0) then
            DO J=1,NUMCOR(IA)
               READ (atomf,101) NCOR(J,IA),LCOR(J,IA),NOBC(J,IA)
               READ (atomf,104) (ECOR(J,IX,IA),IX=1,2)
            end do
         end if
         close(atomf)
      end do
      
C ----------------------------------------------------      
C
      WRITE(IW6,150)
150   FORMAT(/6X,'  LABEL,         IA,',
     &        '      Z,         WS,         WSAV: ')
      DO 355 IA=1,NA
      WRITE(IW6,151) OTXTA(IA),IA,AZ(IA),WS(IA),WSAV(IA)
355   CONTINUE
151   FORMAT(4X,A16,I5,3F12.5)
      WRITE(IW6,153)
153   FORMAT(/6X,'  LABEL,         IA,',
     &        '   NSZRAD, NUMCOR:')
      DO 356 IA=1,NA
      WRITE(IW6,154) OTXTA(IA),IA,
     &         NSZRAD(IA),NUMCOR(IA)
356   CONTINUE
154   FORMAT(4X,A16,I5,2I8)
C
C --------------------- NORMALIZED COMPLEX NODES AND WEIGHTS
C                       FOR CONTOUR INTEGRATION - UNIT IR12
C
      READ (IR12,100) OWORK
      WRITE(IW6,1001) OWORK
C
      DO 370 J=1,8
      READ (IR12,101) NE1
      DO 371 IE=1,NE1
      READ (IR12,104) ZNCN(IE),ZNCW(IE)
371   CONTINUE
      IF(NE1.EQ.NE) GO TO 270
370   CONTINUE
270   CONTINUE
C
      WRITE(IW6,195)
195   FORMAT(/'     **  END OF READING ')
      flush(IW6)
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
214   WRITE(IW6,914)
914   FORMAT(/' **** INPUT ERROR:  IREL MUST BE 0, 1 OR 2 ')
      STOP
215   WRITE(IW6,915)
915   FORMAT(/' **** INPUT ERROR:  NS=2 AND IREL=2 ')
      STOP
216   WRITE(IW6,916)
916   FORMAT(/' **** INPUT ERROR:  IVXC MUST BE 0, 1 OR 2 ')
      STOP
217   WRITE(IW6,917)
917   FORMAT(/' **** INPUT ERROR:  NE GREATER THAN MNE ')
      STOP
218   WRITE(IW6,918)
918   FORMAT(/' *** INPUT ERROR: NE MUST BE 7,10,14,20,',
     &                   '28,40,56, OR 80')
      STOP
219   WRITE(IW6,919)
919   FORMAT(/' **** INPUT ERROR:  ICONT MUST BE 0 OR 1 ')
      STOP
220   WRITE(IW6,920)
920   FORMAT(/' **** INPUT ERROR:  NUH GREATER THAN MNUH ')
      STOP
222   WRITE(IW6,922)
922   FORMAT(/' **** INPUT ERROR: NB GREATER THAN MNB  ')
      STOP
241   WRITE(IW6,941)
941   FORMAT(/' **** INPUT ERROR: IA EXCEEDED MNA  ')
      STOP
245   WRITE(IW6,945)
945   FORMAT(/' **** INPUT ERROR: NSZRAD(IA) GREATER THAN MNR ')
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
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNR=400)
C
      DIMENSION CR23(3),CR31(3),CR12(3),QSC3(3),QSC4(4)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/TEXT/ OTXTA(MNA)
      common/bshift/ ishift,iashift,ilshift,eshift,
     &               ivshift,iavshift,vshift
         
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/
      DATA RC4/4.0D0/,TOL/1.0D-6/,EPSCOR/1.0D-9/
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
C---------------------------------- PRECISION OF CORES
      THRESH=EPSCOR
      WRITE(IW6,116) THRESH
116   FORMAT(/9X,' PRECISION OF CORE ENERGIES: ',G15.7)
C
C------------------------- ACCURACY OF RUNGE-KUTTA METHOD
C                            (1.LE.NSIRK.LE.5)
      NSIRK=2
      WRITE(IW6,117) NSIRK
117   FORMAT(/9X,' ACCURACY OF RUNGE-KUTTA METHOD: ',
     &      '  NSIRK=',I2)
C
C---------------------------------- SHIFT OF ENY'S
      DO 305 IA=1,NA
      DO 306 IS=1,NS
      DO 307 IL=1,NL
      ISHENY(IL,IS,IA)=1

      IF((ishift.eq.1).and.(IA.EQ.iashift).and.(IL.EQ.ilshift) ) 
     & ISHENY(IL,IS,IA)=3
C      IF((ishift.eq.1)) 
C     & ISHENY(IL,IS,IA)=3
          
307   CONTINUE
306   CONTINUE
305   CONTINUE
      WRITE(IW6,118)
118   FORMAT(/2X,'  LABEL,             IA,  IS,',
     &       '    ISHENY - FOR SHIFT OF ENY:')
      DO 308 IA=1,NA
      DO 309 IS=1,NS
      WRITE(IW6,119) OTXTA(IA),IA,IS,
     &            (ISHENY(IL,IS,IA),IL=1,NL)
309   CONTINUE
308   CONTINUE
119   FORMAT(2X,A16,2X,2I5,8X,4I2)
C
C---------------------------------- NCALLH FOR SPH. HARMONICS
      NCALLH=0
C
C---------------------------------- GAUNT FACTORS
      CALL GAUNTY(NL)
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
C----------------------------------- AVERAGE WS-RADII:
C
C                                 DIMENSIONLESS
      AVOL=ROMEGA/REAL(NB)
      AWS=((RC3*AVOL)/(RC4*PI))**(RC1/RC3)
      WRITE(IW6,127) AWS
127   FORMAT(/'    DIMENSIONLESS AVER. WS-RADIUS= ',G15.7)
C
C                                IN ATOMIC UNITS
      SUM=RCZ
      DO 330 IA=1,NA
      SUM=SUM+CON(IA)*WS(IA)**3
330   CONTINUE
      TWS=(SUM/REAL(NB))**(RC1/RC3)
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
      PARAMETER(MNL=4)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=30)
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
     &          INDW(MNLSQ,MNCL),SC(MNLSQ,MNLSQ),SC1(MNLSQ,MNLSQ)
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
      CALL CANSC(NL,WSA,RJ,RI,SC1)
      SC(1:nlsq,1:nlsq)=0.5d0*(SC(1:nlsq,1:nlsq)+
     &transpose(SC1(1:nlsq,1:nlsq)))
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
C       Write (*,*) 'bulk-pre--'
C       Write (*, '(9(1xg22.14))') (aa(i, 1:srows), i=1, nn)
C       Write (*,*) '-----'

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
C       write(*,*)'bulk--'
C       write(*, '(9(1xg22.14))') (BB(i, 1:srows),i=1,nn)      
C       write(*,*)'-----'      
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
           base = base + real (sign(1, f(i)), kind=8) * SC (Abs(f(i)))
C               write(*,*) sc(abs(f(i)))
        End Do
        base = base / real (cntr, kind=8)
        Do i = mptr + 1, mptr + cntr, 1
           SC (Abs(f(i))) = real (sign(1, f(i)), kind=8) * base
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
              acc = acc / real (num, kind=8)
              Do j = 1, num, 1
              SC(Abs(refs(j))) = real(sign(1,refs(j)), kind=8)*acc
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
      PARAMETER(MNL=4)
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
      DFAC(0)= REAL(1)
      DO 310 L=1,LMAXH
      TWOLM1= REAL(2*L-1)
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
      PARAMETER(MNL=4)
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
      TWOLP1=REAL(2*L+1)
      CNH(L,0)=TWOLP1/(RC4*PI)
      DO 302 M=1,L
      PROD=RC1
      DO 303 I=L-M+1,L+M
      PROD=PROD*REAL(I)
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
      TWOLM1=REAL(2*L-1)
      ZQ=TWOLM1*ZP*ZQ
      U(L,-L)=DIMAG(ZQ)
      U(L,L) =DREAL(ZQ)
310   CONTINUE
      DO 312 M=-LMAX+1,LMAX-1
      L= ABS(M)
      TWOLP1=REAL(2*L+1)
      U(L+1,M)=TWOLP1*D3*U(L,M)
312   CONTINUE
      DO 314 M=-LMAX+2,LMAX-2
      MA= ABS(M)
      DO 315 L=MA+1,LMAX-1
      TWOLP1=REAL(2*L+1)
      FLM1=REAL(L+MA)
      FLP1=REAL(L-MA+1)
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
      A=REAL((2*L+1)*(2*L1+1)*(2*L2+1))
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
      PARAMETER(MNL=4)
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
C----------------------------------------------
C        FACTORIAL OF A NON-NEGATIVE INTEGER
C----------------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      TUFF= REAL(1)
      IF(N.LE.1) RETURN
      DO 301 I=2,N
      TUFF=TUFF* REAL(I)
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
C-----------------------------------------------
C   TB-STRUCTURE CONSTANTS IN THE REAL SPACE
C-----------------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNL=4)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=30)
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
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
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
!SMP$ DO SERIAL
      DO 310 IB=1,NB
C
C------------------------- SELECTION OF THE SMALL CLUSTER
      JCL=1
!SMP$ DO SERIAL      
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
!SMP$ DO SERIAL
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
!SMP$ DO SERIAL      
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
c      do jcl=2,NCL
c         write(IW6,*) 'JCL=',jcl
c         do jq=1,nlsq
c            WRITE(IW6,'(9(G10.4,1x))') (STR(jQ,IQ,jcl,IB),IQ=1,NLSQ)
c         end do
c      end do
c
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
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
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
      FX1=(REAL(I1)-RCH)/REAL(2*NK1)
      DO 302 I2=-NK2+1,NK2
      FX2=(REAL(I2)-RCH)/REAL(2*NK2)
      DO 303 I3=1,NK3
      FX3=(REAL(I3)-RCH)/REAL(2*NK3)
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
      DK=PI/(A*REAL(NK))
      DO 3101 I1=1,NK
      DO 3102 I2=1,I1
      DO 3103 I3=1,I2
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC6
      IF(I1.EQ.I2.OR.I2.EQ.I3) WKBZ(IBZ)=RC3
      IF(I1.EQ.I3) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I1)-RCH)*DK
      AKBZ(2,IBZ)=(REAL(I2)-RCH)*DK
      AKBZ(3,IBZ)=(REAL(I3)-RCH)*DK
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
      DK=RC2*PI/(A*REAL(NK))
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
      AKBZ(1,IBZ)=(REAL(I1)-RCH)*DK
      AKBZ(2,IBZ)=(REAL(I2)-RCH)*DK
      AKBZ(3,IBZ)=(REAL(I3)-RCH)*DK
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
      DK=RC2*PI/(A*REAL(NK))
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
      AKBZ(1,IBZ)=(REAL(I1)-RCH)*DK
      AKBZ(2,IBZ)=(REAL(I2)-RCH)*DK
      AKBZ(3,IBZ)=(REAL(I3)-RCH)*DK
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
      DKX=RC2*PI/(A*REAL(NKX))
      DKY=DKX/SQRT(RC3)
      DKZ=PI/(C*REAL(NKZ))
      DO 321 IZ=1,NKZ
      DO 322 NU=1,2
      ANU3=REAL(NU)/RC3
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
      AKBZ(1,IBZ)=(REAL(I-1)+ANU3)*DKX
      AKBZ(2,IBZ)=REAL(IY)*DKY
      AKBZ(3,IBZ)=(REAL(IZ)-RCH)*DKZ
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
      IF(COA.LT.DUM1.OR.COA.GT.DUM2) NKZ=INT(REAL(NK)/COA)
      DK=PI/(A*REAL(NK))
      DKZ=PI/(C*REAL(NKZ))
      DO 341 I1=1,NK
      DO 342 I2=1,I1
      DO 343 I3=1,NKZ
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC2
      IF(I1.EQ.I2) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I1)-RCH)*DK
      AKBZ(2,IBZ)=(REAL(I2)-RCH)*DK
      AKBZ(3,IBZ)=(REAL(I3)-RCH)*DKZ
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
      NKZ=INT(REAL(NK)*P*SR2/C)
      DK=PI/(RC2*P*REAL(NK))
      DKZ=PI/(C*REAL(NKZ))
      DO 345 I=1,NK
      DO 346 J=1,I
      DO 347 I3=1,NKZ
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC2
      IF(I.EQ.J) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=REAL(I-J)*DK
      AKBZ(2,IBZ)=REAL(I+J-1)*DK
      AKBZ(3,IBZ)=(REAL(I3)-RCH)*DKZ
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
      IF(BOA.LT.DUM1.OR.BOA.GT.DUM2) NKY=INT(REAL(NK)/BOA)
      NKZ=NK
      IF(COA.LT.DUM1.OR.COA.GT.DUM2) NKZ=INT(REAL(NK)/COA)
      DKX=PI/(A*REAL(NKX))
      DKY=PI/(B*REAL(NKY))
      DKZ=PI/(C*REAL(NKZ))
      DO 351 I1=1,NKX
      DO 352 I2=1,NKY
      DO 353 I3=1,NKZ
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I1)-RCH)*DKX
      AKBZ(2,IBZ)=(REAL(I2)-RCH)*DKY
      AKBZ(3,IBZ)=(REAL(I3)-RCH)*DKZ
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
CHARACTER
      PARAMETER(MNB=1)
      PARAMETER(MNL=4)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=30)
C
      DIMENSION TRAN(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB),JBVA(MNCL,MNB),
     &            JTRA(3,MNCL,MNB),NSCL(MNB)
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
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
CXXX    MACO    ****
C*******************
      SUBROUTINE MACO
C
C***************************************************
C           MADELUNG CONSTANTS
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
C
      DIMENSION X(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB),AWS,TWS
      COMMON/AMAD/ AMCMM(MNB,MNB)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
C
C--------------------------- EWALD'S TECHNIQUE
C
      DO 320 JB=1,NB
      DO 321 IB=1,NB
      DO 324 I=1,3
      X(I)=POS(I,IB)-POS(I,JB)
324   CONTINUE
      CALL BEWA(X,POT)
      AMCMM(IB,JB)=POT
321   CONTINUE
320   CONTINUE
C
C---------------------------- SCALING OF CONSTANTS BY THE
C                           DIMENSIONLESS AVERAGE WS-RADIUS
      DO 340 JB=1,NB
      DO 341 IB=1,NB
      AMCMM(IB,JB)=AWS*AMCMM(IB,JB)
341   CONTINUE
340   CONTINUE
C
      WRITE(IW6,121)
121   FORMAT(//2X,' **********  DIMENSIONLESS MADELUNG ',
     &    'CONSTANTS  ********** '/
     &   /'        IB,  JB,       CONSTANT: ')
      DO 345 IB=1,NB
      DO 346 JB=1,NB
      WRITE(IW6,123) IB,JB,AMCMM(IB,JB)
346   CONTINUE
345   CONTINUE
123   FORMAT(5X,2I5,5X,G15.7)
C
C---------------------------- SCALING OF CONSTANTS BY THE
C                                TRUE AVERAGE WS-RADIUS
      DO 350 JB=1,NB
      DO 351 IB=1,NB
      AMCMM(IB,JB)=AMCMM(IB,JB)/TWS
351   CONTINUE
350   CONTINUE
C
      RETURN
      END
C*******************
CXXX    BEWA    ****
C*******************
      SUBROUTINE BEWA(XI,POT)
C
C----------------------------------------------
C   BULK EWALD'S TECHNIQUE
C----------------------------------------------
C    INPUT: XI(.) - COORDINATES OF THE POINT
C           + COMMON/BUV/
C   OUTPUT: POT - VALUE OF ELSTAT. POTENTIAL
C----------------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      DIMENSION X(3),XI(3),CX(3),GV(3),RV(3)
C
      COMMON/BUV/ VBR(3,3),VBG(3,3),ROMEGA,GOMEGA
C
      DATA ARGMAX/7.0D0/, CMUL/0.282D0/, TINY/1.0D-4/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/,
     &     RC3/3.0D0/
C
C                                     CONSTANTS
      PI=RC4*ATAN(RC1)
      TWOPI=RC2*PI
      SIG=CMUL*ROMEGA**(RC1/RC3)
      TWOSIG=RC2*SIG
      SIGSQ=SIG**2
C
      CONSTG=RC4*TWOPI/ROMEGA
      CONSTR=RC2/(SIG*SQRT(PI))
C
      GMAX=ARGMAX/SIG
      RMAX=ARGMAX*TWOSIG
C
C                                REDUCED COORDINATES
      DO 302 J=1,3
      COOR=XI(1)*VBG(1,J)+XI(2)*VBG(2,J)+XI(3)*VBG(3,J)
      COOR=COOR/TWOPI
      CX(J)= MOD(COOR,RC1)
302   CONTINUE
      DO 305 I=1,3
      X(I)=CX(1)*VBR(I,1)+CX(2)*VBR(I,2)+CX(3)*VBR(I,3)
305   CONTINUE
C
C                                 G-SUM
C
      DUM=SQRT(VBR(1,1)**2+VBR(2,1)**2+VBR(3,1)**2)/TWOPI
      LG1=1+INT(DUM*GMAX)
      DUM=SQRT(VBR(1,2)**2+VBR(2,2)**2+VBR(3,2)**2)/TWOPI
      LG2=1+INT(DUM*GMAX)
      DUM=SQRT(VBR(1,3)**2+VBR(2,3)**2+VBR(3,3)**2)/TWOPI
      LG3=1+INT(DUM*GMAX)
C
      SUMG=-SIGSQ
C
      DO 311 I1=-LG1,LG1
      DO 312 I2=-LG2,LG2
      DO 313 I3=0,LG3
      IF(I3.EQ.0.AND.I2.LT.0) GO TO 313
      IF(I3.EQ.0.AND.I2.EQ.0.AND.I1.LE.0) GO TO 313
C
      DO 315 J=1,3
      GV(J)=I1*VBG(J,1)+I2*VBG(J,2)+I3*VBG(J,3)
315   CONTINUE
      GSQ=GV(1)**2+GV(2)**2+GV(3)**2
      G=SQRT(GSQ)
      IF(G.GT.GMAX) GO TO 313
C
      ALFA=X(1)*GV(1)+X(2)*GV(2)+X(3)*GV(3)
      ARGEXP=SIGSQ*GSQ
      SUMG = SUMG + RC2*EXP(-ARGEXP)*COS(ALFA)/GSQ
C
313   CONTINUE
312   CONTINUE
311   CONTINUE
C
      SUMG = CONSTG*SUMG
C
C                                 R-SUM
C
      DUM=SQRT(VBG(1,1)**2+VBG(2,1)**2+VBG(3,1)**2)/TWOPI
      LR1=2+INT(DUM*RMAX)
      DUM=SQRT(VBG(1,2)**2+VBG(2,2)**2+VBG(3,2)**2)/TWOPI
      LR2=2+INT(DUM*RMAX)
      DUM=SQRT(VBG(1,3)**2+VBG(2,3)**2+VBG(3,3)**2)/TWOPI
      LR3=2+INT(DUM*RMAX)
C
      SUMR=RCZ
C
      DO 321 I1=-LR1,LR1
      DO 322 I2=-LR2,LR2
      DO 323 I3=-LR3,LR3
C
      DO 325 J=1,3
      RV(J)=I1*VBR(J,1)+I2*VBR(J,2)+I3*VBR(J,3)-X(J)
325   CONTINUE
      R=SQRT(RV(1)**2+RV(2)**2+RV(3)**2)
      IF(R.GT.RMAX) GO TO 323
C
      ARGER=R/TWOSIG
      IF(ARGER.LE.TINY) CLEN=-CONSTR
      IF(ARGER.GT.TINY) CLEN=RC2*ERRFC(ARGER)/R
      SUMR = SUMR + CLEN
C
323   CONTINUE
322   CONTINUE
321   CONTINUE
C
      POT=SUMG+SUMR
C
      RETURN
      END
C*******************
CXXX   ERRFC    ****
C*******************
      REAL*8 FUNCTION ERRFC(X)
C
C-------------------------------------------------------
C   COMPLEMENTARY ERROR FUNCTION
C-------------------------------------------------------
C    ERRFC = ( 2 / SQRT(PI) ) *
C    * INTEGRAL OD X DO NEKONECNA Z FUNKCE EXP(-T**2) DT
C-------------------------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      DATA SQRTPI/1.772453850905516D0/, EPS/1.0D-14/,
     &     ALIM/13.0D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RCH/0.5D0/
C
      ABSX=ABS(X)
      IF(ABSX.LT.ALIM)GO TO 201
C
C                                      EXTREMNI ABSX
      IF(X.GT.RCZ)ERRFC=RCZ
      IF(X.LT.RCZ)ERRFC=RC2
      RETURN
C
201   X2=X**2
      EMX2=EXP(-X2)
C
      IF(ABSX.LT.RC2)GO TO 280
C
C                                         VELKE ABSX
C
      N=0
      PN=RC1
      PNP1=ABSX
      QN=RCZ
      QNP1=RC1
      YLOM=RCZ
      CIT=RCZ
C
222   N=N+1
      YLOMW=YLOM
C
      SUM=ABS(PNP1)+ABS(QNP1)
C
      PNM1=PN/SUM
      PN=PNP1/SUM
      QNM1=QN/SUM
      QN=QNP1/SUM
C
      CIT=CIT+RCH
C
      PNP1=ABSX*PN+CIT*PNM1
      QNP1=ABSX*QN+CIT*QNM1
C
      YLOM=QNP1/PNP1
C
      IF(ABS(YLOM-YLOMW).GT.EPS*ABS(YLOM))GO TO 222
C
      ERRFC=EMX2*YLOM/SQRTPI
C
      IF(X.LT.RCZ)ERRFC=RC2-ERRFC
C
      RETURN
C
C                                          MALE ABSX
C
280   N=0
      PN=RC1
      PNP1=RC1
      QN=RCZ
      QNP1=RC1
      YLOM=RCZ
      ALFA=RC1
      ABSBET=RCZ
      SIGBET=RC1
      TWOX2=RC2*X2
C
233   N=N+1
      YLOMW=YLOM
C
      SUM=ABS(PNP1)+ABS(QNP1)
C
      PNM1=PN/SUM
      PN=PNP1/SUM
      QNM1=QN/SUM
      QN=QNP1/SUM
C
      ALFA=ALFA+RC2
      ABSBET=ABSBET+TWOX2
      SIGBET=-SIGBET
      BETA=SIGBET*ABSBET
C
      PNP1=ALFA*PN+BETA*PNM1
      QNP1=ALFA*QN+BETA*QNM1
C
      YLOM=QNP1/PNP1
C
      IF(ABS(YLOM-YLOMW).GT.EPS*ABS(YLOM))GO TO 233
C
      ERRFC=RC1-RC2*X*EMX2*YLOM/SQRTPI
C
      RETURN
      END
C*******************
CXXX    SLDA    ****
C*******************
      SUBROUTINE SLDA
C
C**********************************************
C   PREPARES THE LDA-ITERATIONS:
C   - DEFINES THE WEIGHTS FOR SCALAR PRODUCT
C   - CALCULATES THE VALENCE WAVE FUNCTIONS
C   - RECALCULATES THE POTENTIAL PARAMETERS
C**********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/WAF/ PHI(MNR,MNL,MNS,MNA),
     &            PHID(MNR,MNL,MNS,MNA),
     &            PHIDD(MNR,MNL,MNS,MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/LIN/ ENY1,PPC1,PPD1,PPQ1,PPP1,DNY1,FINY1,FINYD1,
     &            FI(MNR),FID(MNR),FIDD(MNR)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/AMXP/ DXP(MNAM,MNUH)
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
      COMMON/TEXT/ OTXTA(MNA)
C
      DATA RCZ/0.0D0/
C
      DATA CLIM/0.01D0/,PLIM/0.01D0/,AZL/2.5D0/,DMUL/0.03D0/
C
      DO 301 IH=1,MNUH
      DO 302 J=1,MNAM
      DXP(J,IH)=RCZ
302   CONTINUE
301   CONTINUE
      DO 303 IH=1,MNUH
      DO 304 J=1,MNAM
      DFP(J,IH)=RCZ
304   CONTINUE
303   CONTINUE
      DO 305 JH=1,MNUH
      DO 306 IH=1,MNUH
      VOMA(IH,JH)=RCZ
306   CONTINUE
305   CONTINUE
C
C                            WEIGHTS FOR SCALAR PRODUCT
      J=0
      DO 310 IA=1,NA
      NR1=NSZRAD(IA)
      WS1=WS(IA)
      CALL RAPO(NR1,WS1,R)
CCC      PREF=RC1
      PREF=MAX(CLIM,CON(IA))
C
      DO 311 IS=1,NS
      DO 312 I=2,NR1
      J=J+1
      SP(J)=PREF*R(I)**2*(R(I)-R(I-1))
312   CONTINUE
311   CONTINUE
310   CONTINUE
      NAM=J
C
      WRITE(IW6,105) NAM
105   FORMAT(//5X,'*** LENGTH OF VECTORS FOR ANDERSON',
     &            ' MIXING= ',I8)
C
C                                       ENERGY CONTOUR
C
      DO 315 IE=1,NE
      ZCN(IE)=DIAM*ZNCN(IE)+EF
      ZCW(IE)=DIAM*ZNCW(IE)
315   CONTINUE
      WRITE(IW6,115) ZCN(1),ZCN(NE)
115   FORMAT(//3X,'1ST ENERGY CONTOUR: 1ST NODE= ',2G15.7/
     &        3X,'                    LAST NODE= ',2G15.7)
C
C
C                                         LOOP OVER ATOMS
      DO 320 IA=1,NA
C                                         RADIAL MESH
      AZ1=AZ(IA)
      NR1=NSZRAD(IA)
      IREL1=MIN(IREL,1)
      IF(AZ1.LT.AZL) IREL1=0
      WS1=WS(IA)
      WSAV1=WSAV(IA)
      CALL RAPO(NR1,WS1,R)
C                                         LOOP OVER SPINS
      DO 321 IS=1,NS
C                                          POTENTIAL
      DO 322 I=1,NR1
      V(I)=POT(I,IS,IA)
322   CONTINUE
      CALL LIPO(AZ1,NR1)
C                                           LOOP OVER L
      DO 323 IL=1,NL
C
      L=IL-1
      ENY1=ENY(IL,IS,IA)
      ARG=MAX(PPP(IL,IS,IA),PLIM)
      DE=DMUL/SQRT(ARG)
C
      CALL RSEL(AZ1,DE,WSAV1,L,IREL1,NR1)
C
      PPC(IL,IS,IA)=PPC1
      PPD(IL,IS,IA)=PPD1
      PPQ(IL,IS,IA)=PPQ1
      PPP(IL,IS,IA)=PPP1
      DNY(IL,IS,IA)=DNY1
      FINY(IL,IS,IA)=FINY1
      FINYD(IL,IS,IA)=FINYD1
C                              STORAGE OF PHI,PHID,PHIDD
      DO 325 I=1,NR1
      PHI(I,IL,IS,IA)=FI(I)
      PHID(I,IL,IS,IA)=FID(I)
      PHIDD(I,IL,IS,IA)=FIDD(I)
325   CONTINUE
C
323   CONTINUE
321   CONTINUE
320   CONTINUE
C
C---------------------------------------    PRINT
      WRITE(IW6,160)
160   FORMAT(//5X,' ******  INPUT POTENTIAL PARAMETERS  : ')
      DO 330 IA=1,NA
      DO 331 IS=1,NS
      WRITE(IW6,161)IA,OTXTA(IA),IS,AZ(IA),WS(IA),WSAV(IA)
161   FORMAT(/4X,'IA=',I3,8X,A16,8X,'IS=',I2/
     &  5X,'AZ=',F7.3,8X,'WS=',F10.6,8X,'WSAV=',F10.6)
      WRITE(IW6,162)(ENY(IL,IS,IA),IL=1,NL)
162   FORMAT(1X,' ENY : ',4G15.7)
      WRITE(IW6,163)(PPC(IL,IS,IA),IL=1,NL)
163   FORMAT(1X,'  C  : ',4G15.7)
      WRITE(IW6,164)(PPD(IL,IS,IA),IL=1,NL)
164   FORMAT(1X,'DELTA: ',4G15.7)
      WRITE(IW6,165)(PPQ(IL,IS,IA),IL=1,NL)
165   FORMAT(1X,'  Q  : ',4G15.7)
      WRITE(IW6,166)(PPP(IL,IS,IA),IL=1,NL)
166   FORMAT(1X,'  P  : ',4G15.7)
      WRITE(IW6,167)(DNY(IL,IS,IA),IL=1,NL)
167   FORMAT(1X,' DNY : ',4G15.7)
      WRITE(IW6,168)(FINY(IL,IS,IA),IL=1,NL)
168   FORMAT(1X,'FINY : ',4G15.7)
      WRITE(IW6,169)(FINYD(IL,IS,IA),IL=1,NL)
169   FORMAT(1X,'FINYD: ',4G15.7)
331   CONTINUE
330   CONTINUE
C
      RETURN
      END
C*******************
CXXX    SCPF    ****
C*******************
      SUBROUTINE SCPF
C
C**************************************************
C        SETS THE COHERENT POTENTIAL FUNCTIONS
C**************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
C
      DIMENSION OWORK(5),ZP(MNLSQ),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
C
      DATA RCZ/0.0D0/,TOL/1.0D-5/
C
      NLSQ=NL**2
      ZZ=DCMPLX(RCZ,RCZ)
C
      WRITE(IW6,111) IITER
111   FORMAT(//5X,'-----------  ITERATION:  ITER=',
     &       I4,'  -----------')
C
      IF(IITER.GT.1) GO TO 250
C
C******************************** FOR THE FIRST ITERATION
C
      DO 310 IS=1,MNS
      DO 311 IE=1,MNE
      DO 312 IB=1,MNB
      DO 313 J=1,MNLSQ
      DO 314 I=1,MNLSQ
      ZOMG(I,J,IB,IE,IS)=ZZ
314   CONTINUE
313   CONTINUE
312   CONTINUE
311   CONTINUE
310   CONTINUE
C
C      IF(ICONT.EQ.0) GO TO 250
      GO TO 250
C
C-------------------------------  CONTINUATION OF
C                                 A PREVIOUS RUN (ICONT=1)
100   FORMAT(5A16)
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
C
!      READ (IR4,100) OWORK
!      READ (IR4,101) NBX,NLX,NSX,NEX
!      READ (IR4,104) DIAMX,EFX
!      IDIF= ABS(NBX-NB)+ ABS(NLX-NL)
!     &    + ABS(NSX-NS)+ ABS(NEX-NE)
!      DIF = ABS(DIAMX-DIAM)+ ABS(EFX-EF)
!      IF(IDIF.GT.0.OR.DIF.GT.TOL) GO TO 295
!      DO 340 IS=1,NS
!      DO 341 IE=1,NE
!      DO 342 IB=1,NB
!      READ (IR4,104) ((ZCP1(I,J),I=J,NLSQ),J=1,NLSQ)
!      DO 343 J=1,NLSQ
!      DO 344 I=J,NLSQ
!      ZCP1(J,I)=ZCP1(I,J)
!344   CONTINUE
!343   CONTINUE
!      DO 345 J=1,NLSQ
!      DO 346 I=1,NLSQ
!      ZCPF(I,J,IB,IE,IS)=ZCP1(I,J)
!346   CONTINUE
!345   CONTINUE
!342   CONTINUE
!341   CONTINUE
!340   CONTINUE
!      WRITE(IW6,140)
140   FORMAT(/5X,' *****  CPA STARTED FROM',
     &       ' A PREVIOUS RUN  *****'/)
!      WRITE(IW6,100) OWORK
!C
!      RETURN
!C
C********************************* GENERAL PART
C
250   CONTINUE
C
C-------------------------------- LOOP OVER SPIN AND ENERGY
      DO 360 IS=1,NS
      DO 361 IE=1,NE
      ZE=ZCN(IE)
C
C------------------------------- LOOP OVER SITES
      IA=0
      DO 362 IB=1,NB
C
      DO 364 J=1,NLSQ
      DO 365 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IB,IE,IS)
365   CONTINUE
364   CONTINUE
C
      DO 367 J=1,NLSQ
      DO 368 I=1,NLSQ
      ZFI1(I,J)=ZZ
368   CONTINUE
367   CONTINUE
C
C..............................  LOOP OVER COMPONENTS
      DO 370 IC=1,NC(IB)
      IA=IA+1
C
      DO 372 IL=1,NL
      PPCX=PPC(IL,IS,IA)
      PPDX=PPD(IL,IS,IA)
      PPQX=PPQ(IL,IS,IA)
      ALFX=QSCR(IL)
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      CALL PLMZ(PPCX,PPDX,PPQX,ALFX,ZE,ZPFX,ZLAX,ZMUX)
!SMP$ DO SERIAL
      DO 373 I=ISTA,IFIN
      ZP(I)=ZPFX
373   CONTINUE
372   CONTINUE
C
      DO 375 J=1,NLSQ
      DO 376 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
376   CONTINUE
375   CONTINUE
      DO 377 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZP(I)
377   CONTINUE
C
      CALL CINVC(ZW,NLSQ,MNLSQ)
C
      DO 378 J=1,NLSQ
      DO 379 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+CON(IA)*ZW(I,J)
379   CONTINUE
378   CONTINUE
C
C.......................... END OF LOOP OVER COMPONENTS
370   CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
C
      DO 381 J=1,NLSQ
      DO 382 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
382   CONTINUE
381   CONTINUE
C
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 384 J=1,NLSQ
      DO 385 I=1,NLSQ
      ZCPF(I,J,IB,IE,IS)=ZCP1(I,J)
385   CONTINUE
384   CONTINUE
C
C-------------------------- END OF LOOP OVER SITES
362   CONTINUE
C
C-------------------------- END OF LOOP OVER
C                                   ENERGY AND SPIN
361   CONTINUE
360   CONTINUE
C
      RETURN
C
C------------------------------ ERROR MESSAGE
C
295   WRITE(IW6,195)
195   FORMAT(/3X,' **** INCONSISTENCY IN CPF-DATA ')
      STOP
      END
C*******************
CXXX    ONGF    ****
C*******************
      SUBROUTINE ONGF
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
      PARAMETER(MNB=1)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZWA(MNBLSQ,MNBLSQ),ZWG(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSK(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 IS=1,MNS
      DO 302 IE=1,MNE
      DO 303 IB=1,MNB
      DO 304 J=1,MNLSQ
      DO 305 I=1,MNLSQ
      ZGFG(I,J,IB,IE,IS)=ZZ
305   CONTINUE
304   CONTINUE
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
C--------------------------------- SPIN LOOP AND
C                                  ENERGY LOOP
      DO 310 IS=1,NS
      DO 311 IE=1,NE
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
      ZWA(II0+I,II0+J)=ZWA(II0+I,II0+J)+ZCPF(I,J,IB,IE,IS)
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
      ZGFG(I,J,IB,IE,IS) = ZGFG(I,J,IB,IE,IS)
     &            + WKBZ(IBZ) * ZWA(II0+I,II0+J)
327   CONTINUE
326   CONTINUE
324   CONTINUE
C---------------------------------- END OF IBZ-LOOP, ENERGY
C                                   LOOP, AND SPIN LOOP
312   CONTINUE
311   CONTINUE
310   CONTINUE
C
C---------------------------------- SYMMETRIZATION OF THE
C                                   ON-SITE BLOCKS
      DO 350 IS=1,NS
      DO 351 IE=1,NE
      DO 352 IB=1,NB
C
      DO 354 J=1,NLSQ
      DO 355 I=1,NLSQ
      ZWG(I,J)=ZGFG(I,J,IB,IE,IS)
355   CONTINUE
354   CONTINUE
C
      CALL SYMGF(ZWG,NSYM,NL)
C
      DO 357 J=1,NLSQ
      DO 358 I=1,NLSQ
      ZGFG(I,J,IB,IE,IS)=ZWG(I,J)
358   CONTINUE
357   CONTINUE
C
352   CONTINUE
351   CONTINUE
350   CONTINUE
C
      RETURN
      END
C*******************
CXXX    CPAIT   ****
C*******************
      FUNCTION CPAIT(inarg) result(ICONV)
C      subroutine CPAIT      
C
C****************************************************
C   PERFORMS ONE CPA-ITERATION, CALCULATES
C   CONFIGURATIONALLY AVERAGED GREENS FUNCTIONS,
C   AND WRITES THE RESULTING CPF FUNCTIONS
C****************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=160002)
      PARAMETER(MNLSQ=MNL**2)
C
      DIMENSION NONC(MNB,MNS),ZP(MNLSQ),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(3,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNB,MNE,MNS)
      COMMON/CAGF/ ZCAGF(MNLSQ,MNLSQ,MNA,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR      
C     
      integer::ICONV
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2      
C
      QMIXI=QMIXC
      IF(IITER.EQ.1) QMIXI=RC1
C
      DO 301 IS=1,MNS
      DO 302 IB=1,MNB
      NONC(IB,IS)=0

302   CONTINUE
301   CONTINUE
C
C-------------------------------- LOOP OVER SPIN AND ENERGY
C
      DO 310 IS=1,NS
      DO 311 IE=1,NE
      ZE=ZCN(IE)
C
C------------------------------- LOOP OVER SITES
C
      IA=0
      DO 312 IB=1,NB
C
      DO 314 J=1,NLSQ
      DO 315 I=1,NLSQ
      ZFI1(I,J)=ZGFG(I,J,IB,IE,IS)
315   CONTINUE
314   CONTINUE
      DO 316 J=1,NLSQ
      DO 317 I=1,NLSQ
      ZCP1(I,J)=ZCPF(I,J,IB,IE,IS)
317   CONTINUE
316   CONTINUE
      DO 318 J=1,NLSQ
      DO 319 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IB,IE,IS)
319   CONTINUE
318   CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
C
      DO 320 J=1,NLSQ
      DO 321 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
321   CONTINUE
320   CONTINUE
C
      SUM=RCZ
      SUM1=RCZ
      DO 324 J=1,NLSQ
      DO 325 I=1,NLSQ
      XW=DREAL(ZW(I,J))
      YW=DIMAG(ZW(I,J))
      SUM=SUM+ ABS(XW)+ ABS(YW)
      XW1=DREAL(ZOM1(I,J))
      YW1=DIMAG(ZOM1(I,J))
      SUM1=SUM1+ ABS(XW1)+ ABS(YW1)
325   CONTINUE
324   CONTINUE
      IF(SUM.GT.(DIFMC*SUM1)) NONC(IB,IS)=NONC(IB,IS)+1
C
      DO 327 J=1,NLSQ
      DO 328 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXI*ZW(I,J)
328   CONTINUE
327   CONTINUE
C
      CALL SYMGF(ZOM1,NSYM,NL)
C
C..............................  LOOP OVER COMPONENTS
C
      DO 330 IC=1,NC(IB)
      IA=IA+1
C
      DO 332 IL=1,NL
      PPCX=PPC(IL,IS,IA)
      PPDX=PPD(IL,IS,IA)
      PPQX=PPQ(IL,IS,IA)
      ALFX=QSCR(IL)
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      CALL PLMZ(PPCX,PPDX,PPQX,ALFX,ZE,ZPFX,ZLAX,ZMUX)
!SMP$ DO SERIAL
      DO 333 I=ISTA,IFIN
      ZP(I)=ZPFX
333   CONTINUE
332   CONTINUE
C
      DO 335 J=1,NLSQ
      DO 336 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
336   CONTINUE
335   CONTINUE
      DO 337 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZP(I)
337   CONTINUE
C
      CALL CINVC(ZW,NLSQ,MNLSQ)
C
      DO 338 J=1,NLSQ
      DO 339 I=1,NLSQ
      ZCAGF(I,J,IA,IE,IS)=ZW(I,J)
339   CONTINUE
338   CONTINUE
C
C.......................... END OF LOOP OVER COMPONENTS
330   CONTINUE
C
      DO 342 J=1,NLSQ
      DO 343 I=1,NLSQ
      ZOMG(I,J,IB,IE,IS)=ZOM1(I,J)
343   CONTINUE
342   CONTINUE
C
C-------------------------- END OF LOOP OVER SITES
312   CONTINUE
C
C-------------------------- END OF LOOP OVER
C                                   ENERGY AND SPIN
311   CONTINUE
310   CONTINUE
C
      WRITE(IW6,111) IITER
111   FORMAT(/'  ****  CONVERGENCE OF COHERENT INTERACTOR :'/
     &  '        IB,      NON-CONV. OMEGA:',
     &   20X,'ITER=',I4)
      ICONV=0
      DO 360 IB=1,NB
      WRITE(IW6,120) IB,(NONC(IB,IS),IS=1,NS)
!SMP$ DO SERIAL      
      do IS=1,NS
         ICONV=ICONV+NONC(IB,IS)
      end do
360   CONTINUE

120   FORMAT(5X,I5,10X,2I5)
C
C                                  OUTPUT TO UNIT IW8
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
C
      MODPR= MOD(IITER,NFUPR)
      IF(MODPR.NE.0.AND.IITER.NE.1
     &             .AND.IITER.NE.NITER) RETURN
C
      
      REWIND IW8
      WRITE(IW8,181)
181   FORMAT(2X,' ----------  CPF-FILE  ----------')
      WRITE(IW8,101) NB,NL,NS,NE
      WRITE(IW8,104) DIAM,EF
      DO 370 IS=1,NS
      DO 371 IE=1,NE
      DO 372 IB=1,NB
      WRITE(IW8,104) ((ZCPF(I,J,IB,IE,IS),I=J,NLSQ),J=1,NLSQ)
372   CONTINUE
371   CONTINUE
370   CONTINUE
      WRITE(IW8,182)
182   FORMAT(2X,' --------------------------------')
C
      RETURN
      END FUNCTION CPAIT
C*******************
CXXX    CMOM    ****
C*******************
      SUBROUTINE CMOM
C
C************************************************
C   CALCULATES MOMENTS OF DOS VIA COMPLEX
C   CONTOUR INTEGRATION OF CONFIGURATIONALLY
C   AVERAGED PHYSICAL GREENS FUNCTIONS
C   AND ITERATES THE FERMI ENERGY
C************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNLSQ=MNL**2)
C
      DIMENSION ZGDI(MNE,MNL,MNS,MNA),DOSP(MNL,MNS,MNA),
     &          ZW0(MNE),ZW1(MNE),ZW2(MNE)
C
      COMMON/CAGF/ ZCAGF(MNLSQ,MNLSQ,MNA,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITEF/ DIFME,QMIXE
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      PI=RC4*ATAN(RC1)
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 IA=1,MNA
      DO 302 IS=1,MNS
      DO 303 IL=1,MNL
      DO 304 IE=1,MNE
      ZGDI(IE,IL,IS,IA)=ZZ
304   CONTINUE
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
C--------------------------------- DIAGONAL ELEMENTS
C
      DO 310 IA=1,NA
      DO 311 IS=1,NS
      DO 312 IL=1,NL
      PPC1=PPC(IL,IS,IA)
      PPD1=PPD(IL,IS,IA)
      PPQ1=PPQ(IL,IS,IA)
      ALF1=QSCR(IL)
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 315 IE=1,NE
      ZE=ZCN(IE)
      CALL PLMZ(PPC1,PPD1,PPQ1,ALF1,ZE,ZPF1,ZLA1,ZMU1)

      ZDP1=ZMU1**2
      DO 317 I=ISTA,IFIN
      ZGDI(IE,IL,IS,IA)=ZGDI(IE,IL,IS,IA)
     &        + ZLA1 + ZDP1 * ZCAGF(I,I,IA,IE,IS)
317   CONTINUE
315   CONTINUE
312   CONTINUE
311   CONTINUE
310   CONTINUE
C
C---------------------------------- INTEGRATIONS
C
      DO 330 IA=1,NA
      DO 331 IS=1,NS
      DO 332 IL=1,NL
      ENY1=ENY(IL,IS,IA)
!SMP$ DO SERIAL      
      DO 334 IE=1,NE
      ZEPS=ZCN(IE)-ENY1
      ZW0(IE)=ZGDI(IE,IL,IS,IA)
      ZW1(IE)=ZW0(IE)*ZEPS
      ZW2(IE)=ZW1(IE)*ZEPS
334   CONTINUE
      EMDI0(IL,IS,IA)=RCINT(ZW0)
      EMDI1(IL,IS,IA)=RCINT(ZW1)
      EMDI2(IL,IS,IA)=RCINT(ZW2)
332   CONTINUE
331   CONTINUE
330   CONTINUE
C
C-------------------------------- CHARGE DEVIATION
      CHDEV=RCZ
      DO 350 IA=1,NA
      CHTR=-VALZ(IA)
      DO 351 IL=1,NL
      CHTR=CHTR+EMDI0(IL,1,IA)+EMDI0(IL,NS,IA)
351   CONTINUE
      CHDEV=CHDEV+CON(IA)*CHTR
350   CONTINUE
C
C-------------------------------- ESTIMATION OF DOS(EF)
      DO 353 IA=1,NA
      DO 354 IS=1,NS
      DO 355 IL=1,NL
      DOSP(IL,IS,IA)=DIMAG(ZGDI(NE,IL,IS,IA))/PI
355   CONTINUE
354   CONTINUE
353   CONTINUE
      DOST=RCZ
      DO 357 IA=1,NA
      DO 358 IL=1,NL
      DOST = DOST
     &     + CON(IA)*(DOSP(IL,1,IA)+DOSP(IL,NS,IA))
358   CONTINUE
357   CONTINUE
C
C-------------------------------- CORRECTION TO EF AND
C                                 TO THE ENERGY MOMENTS
      EFDEV=-CHDEV/DOST
      DO 360 IA=1,NA
      DO 361 IS=1,NS
      DO 362 IL=1,NL
      ADD=EFDEV*DOSP(IL,IS,IA)
      FEPS=EF-ENY(IL,IS,IA)
      EMDI0(IL,IS,IA)=EMDI0(IL,IS,IA)+ADD
      EMDI1(IL,IS,IA)=EMDI1(IL,IS,IA)+ADD*FEPS
      EMDI2(IL,IS,IA)=EMDI2(IL,IS,IA)+ADD*FEPS**2
362   CONTINUE
361   CONTINUE
360   CONTINUE
C
C-------------------------------- ESTIMATION OF NEW EF
C                               AND NEW ENERGY CONTOUR
      DEF=EFDEV*QMIXE
      DEF= MIN(DEF,DIFME)
      DEF= MAX(DEF,-DIFME)
      EF=EF+DEF
C
      DO 365 IE=1,NE
      ZCN(IE)=DIAM*ZNCN(IE)+EF
365   CONTINUE
C
      WRITE(IW6,161) CHDEV,EFDEV,IITER
161   FORMAT(/2X,'CHARGE DEV.=',G12.3,6X,'EF DEV.=',G12.3,
     &   6X,'ITER=',I4)
      WRITE(IW6,162) EF,IITER
162   FORMAT(/2X,'FERMI LEVEL=',F13.7,26X,'ITER=',I4)
C
      RETURN
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
      PARAMETER(MNB=1)
      PARAMETER(MNL=4)
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
CXXX    RCINT   ****
C*******************
      REAL*8 FUNCTION RCINT(ZF)
C
C********************************************************
C    INTEGRAL (DIVIDED BY 2*PI*I) OF A FUNCTION OVER
C    THE CLOSED COMPLEX CONTOUR WITH
C    NODES AND WEIGHTS GIVEN IN COMMON/CNW/
C    (WHOLE-CIRCLE CONTOUR)
C********************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNE=14)
C
      DIMENSION ZF(MNE)
C
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
C
      DATA RCZ/0.0D0/
C
      RCINT=RCZ
      DO 310 IE=1,NE
      ZD=ZCW(IE)*ZF(IE)
      RCINT=RCINT+DREAL(ZD)
310   CONTINUE
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
      PARAMETER(MNL=4)
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
CXXX    CODE    ****
C*******************
      SUBROUTINE CODE
C
C************************************
C   CALCULATES THE CORE DENSITIES AND
C   THE CORE EIGENVALUES
C************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
C
      DATA C274/274.074D0/, AZL/2.5D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
C
      DO 300 IA=1,MNA
      DO 301 IS=1,MNS
      DO 302 I=1,MNR
      RHOCOR(I,IS,IA)=RCZ
302   CONTINUE
301   CONTINUE
300   CONTINUE
C
C                                        LOOP OVER ATOMS
      DO 310 IA=1,NA
C
      NUMC1=NUMCOR(IA)
      IF(NUMC1.EQ.0) GO TO 310
C
      NR1=NSZRAD(IA)
      AZ1=AZ(IA)
      WS1=WS(IA)
C                                          RADIAL MESH
      CALL RAPO(NR1,WS1,R)
C
      IREL1=IREL
      IF(AZ1.LT.AZL) IREL1=0
      IF(IREL1.EQ.2) GO TO 220
C.................................... NON-RELATIVISTIC AND
C                                  SCALAR-RELATIVISTIC CASE
C                              (INCLUDING SPIN POLARIZATION)
      UC=REAL(IREL1)/C274
      UCSQ=UC**2
C                                      LOOP OVER SPIN
      DO 312 IS=1,NS
C                                          POTENTIAL
      DO 313 I=1,NR1
      V(I)=POT(I,IS,IA)
313   CONTINUE
      CALL LIPO(AZ1,NR1)
C                                 LOOP OVER CORE ORBITALS
      DO 315 J=1,NUMC1
C
      E=ECOR(J,IS,IA)
      N=NCOR(J,IA)
      L=LCOR(J,IA)
      NOBC1=NOBC(J,IA)
C                                   SCHROEDINGER EQUATION
      CALL RSEC(AZ1,E,THRESH,N,L,IREL1,NR1,NOBC1)
C
C                                ADDITION TO CORE DENSITY
      CMUL=REAL(2*L+1)/PI4
      BLAM=REAL(L*(L+1))
      RHOCOR(1,IS,IA)=RHOCOR(1,IS,IA) + CMUL*WG(1)**2
      DO 317 I=2,NR1
      AME=RC1+UCSQ*(E-V(I))
      FRAC=WG(I)/(R(I)*AME)
      DENS = WG(I)**2 + UCSQ * (WF(I)**2 + BLAM*FRAC**2)
      RHOCOR(I,IS,IA)=RHOCOR(I,IS,IA) + CMUL*DENS
317   CONTINUE
C                                   EIGENVALUE
      ECOR(J,IS,IA)=E
C
315   CONTINUE
312   CONTINUE
C
      IF(NS.EQ.1) THEN
      DO 319 J=1,NUMC1
      ECOR(J,2,IA)=ECOR(J,1,IA)
319   CONTINUE
      END IF
      GO TO 310
C.................................. FULLY RELATIVISTIC CASE
C                                (WITHOUT SPIN POLARIZATION)
220   UC=RC1/C274
      UCSQ=UC**2
C                                          POTENTIAL
      DO 323 I=1,NR1
      V(I)=POT(I,1,IA)
323   CONTINUE
      CALL LIPO(AZ1,NR1)
C                                 LOOP OVER CORE ORBITALS
      DO 325 J=1,NUMC1
C
      N=NCOR(J,IA)
      L=LCOR(J,IA)
      NOBC1=NOBC(J,IA)
C                                      LOOP OVER KAPPA
      DO 326 IK=1,2
C
      IF(IK.EQ.1) KAPPA=-L-1
      IF(IK.EQ.2) KAPPA=L
      IF(KAPPA.EQ.0) GO TO 222
C
      E=ECOR(J,IK,IA)
C                                        DIRAC EQUATION
      CALL RDEC(AZ1,E,THRESH,N,KAPPA,NR1,NOBC1)
C
C                                ADDITION TO CORE DENSITY
      CMUL=REAL(ABS(KAPPA))/PI4
      DO 327 I=1,NR1
      DENS = WG(I)**2 + UCSQ * WF(I)**2
      RHOCOR(I,1,IA)=RHOCOR(I,1,IA) + CMUL*DENS
327   CONTINUE
C                                   EIGENVALUE
222   ECOR(J,IK,IA)=E
C
326   CONTINUE
325   CONTINUE
C
310   CONTINUE
C
      RETURN
      END
C*******************
CXXX    VADE    ****
C*******************
      SUBROUTINE VADE
C
C*****************************************
C   CALCULATES THE VALENCE DENSITIES
C   FROM THE L-DIAGONAL MOMENTS OF DOS
C*****************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
C
      COMMON/WAF/ PHI(MNR,MNL,MNS,MNA),
     &            PHID(MNR,MNL,MNS,MNA),
     &            PHIDD(MNR,MNL,MNS,MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
C
      DO 300 IA=1,MNA
      DO 301 IS=1,MNS
      DO 302 I=1,MNR
      RHOVAL(I,IS,IA)=RCZ
302   CONTINUE
301   CONTINUE
300   CONTINUE
C
      DO 330 IA=1,NA
      DO 331 IS=1,NS
      DO 332 IL=1,NL
C
      WM0=EMDI0(IL,IS,IA)/PI4
       TM1=RC2*EMDI1(IL,IS,IA)/PI4
      WM2=EMDI2(IL,IS,IA)/PI4
C
      DO 335 I=1,NSZRAD(IA)
C
      F=PHI(I,IL,IS,IA)
      FD=PHID(I,IL,IS,IA)
      FDD=PHIDD(I,IL,IS,IA)
C
      RHOVAL(I,IS,IA) =  RHOVAL(I,IS,IA)
     &     +  WM0*F**2 + TM1*F*FD + WM2*(FD**2+F*FDD)
C
335   CONTINUE
332   CONTINUE
331   CONTINUE
330   CONTINUE
C
      RETURN
      END
C*******************
CXXX    NEPO    ****
C*******************
      SUBROUTINE NEPO
C
C**********************************************
C   CALCULATES NEW ONE-ELECTRON POTENTIALS
C   AND TOTAL ENERGY CONTRIBUTIONS
C**********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
      PARAMETER(MNE=14)

C     
      DIMENSION YMAX(2,MNA),IMAX(2,MNA),
     &          CHTRB(MNB),VMADB(MNB),WDUM(MNR),
     &          RHOT(MNR),RHOS(MNR,MNS),VXC(2),
     &          VHART(MNR),FPOT(MNR,MNS,MNA)
C
      COMMON/AMAD/ AMCMM(MNB,MNB)
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/WAF/ PHI(MNR,MNL,MNS,MNA),
     &            PHID(MNR,MNL,MNS,MNA),
     &            PHIDD(MNR,MNL,MNS,MNA)
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/LIN/ ENY1,PPC1,PPD1,PPQ1,PPP1,DNY1,FINY1,FINYD1,
     &            FI(MNR),FID(MNR),FIDD(MNR)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA)
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),
     &            VMAD(MNA),ACTE(MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/AMXP/ DXP(MNAM,MNUH)
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
      COMMON/TEXT/ OTXTA(MNA)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/SCREEN/ QSCR(MNL)
      common/bshift/ ishift,iashift,ilshift,eshift,
     &               ivshift,iavshift,vshift

C     
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/,
     &     RCH/0.5D0/
C
      DATA AZL/2.5D0/, DMUL/0.03D0/
C
C------------------------------------------------------
C    WITH ANDERSON MIXING FOR THE POTENTIALS
C    (V. EYERT, J. COMPUT. PHYS. 124 (1996) 271)
C------------------------------------------------------
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
C
C---------------------------- INTERATOMIC MADELUNG TERMS
C
C                            ATOMIC CHARGE TRANSFERS
C                            AND MAGNETIC MOMENTS
      DO 301 IA=1,NA
      CHATRA(IA)=-VALZ(IA)
      AMGMOM(IA)=RCZ
      DO 302 IL=1,NL
      CHATRA(IA)=CHATRA(IA)+EMDI0(IL,1,IA)+EMDI0(IL,NS,IA)
      AMGMOM(IA)=AMGMOM(IA)+EMDI0(IL,1,IA)-EMDI0(IL,NS,IA)
302   CONTINUE
301   CONTINUE
C                           SITE-RESOLVED CHARGE TRANSFERS
      IA=0
      DO 303 IB=1,NB
      CHTRB(IB)=RCZ
      DO 304 IC=1,NC(IB)
      IA=IA+1
      CHTRB(IB)=CHTRB(IB)+CON(IA)*CHATRA(IA)
304   CONTINUE
303   CONTINUE
C                           SITE-RESOLVED MADELUNG TERMS
      DO 305 IB=1,NB
      VMADB(IB)=RCZ
305   CONTINUE
      DO 306 JB=1,NB
      DO 307 IB=1,NB
      VMADB(IB)=VMADB(IB) + AMCMM(IB,JB)*CHTRB(JB)
307   CONTINUE
306   CONTINUE
C
C                                 ATOMIC MADELUNG TERMS
      IA=0
      DO 308 IB=1,NB
      DO 309 IC=1,NC(IB)
      IA=IA+1
      VMAD(IA)=VMADB(IB)
309   CONTINUE
308   CONTINUE
C
C-----------------------------------  OUT-POTENTIALS AND
C                                   ENERGY CONTRIBUTIONS
C
C                                    LOOP OVER ATOMS
      DO 310 IA=1,NA
C                                        RADIAL MESH
      NR1=NSZRAD(IA)
      AZ1=AZ(IA)
      WS1=WS(IA)
      CALL RAPO(NR1,WS1,R)
C                                   SPIN DEPENDENT DENSITY
      DO 311 IS=1,NS
      DO 312 I=1,NR1
      RHOS(I,IS)=RHOCOR(I,IS,IA)+RHOVAL(I,IS,IA)
312   CONTINUE
311   CONTINUE
C                                     TOTAL CHARGE DENSITY
      DO 313 I=1,NR1
      RHOT(I)=RHOS(I,1)+RHOS(I,NS)
313   CONTINUE
C                                      HARTREE POTENTIAL
      CALL HAPO(NR1,R,RHOT,VHART)
C                                    LOOP OVER RADIAL POINTS
      WDUM(1)=RCZ
      DO 314 I=2,NR1
      RHOUP=RHOS(I,1)
      RHODO=RHOS(I,NS)
      IF(IVXC.EQ.0) CALL XCVBH(RHOUP,RHODO,VXC(1),VXC(2),EXC)
      IF(IVXC.EQ.1) CALL XCCAPZ(RHOUP,RHODO,VXC(1),VXC(2),EXC)
      IF(IVXC.EQ.2) CALL XCVWN(RHOUP,RHODO,VXC(1),VXC(2),EXC)
C
      VNUC=-RC2*AZ1/R(I)
      VCOUL=VNUC+VHART(I)+VMAD(IA)
      DO 315 IS=1,NS
      FPOT(I,IS,IA)=VCOUL+VXC(IS)   
C       if ((ivshift.eq.1).and.(iavshift.eq.IA)) then
C        FPOT(I,IS,IA)=FPOT(I,IS,IA)+vshift
C       end if
         
315   CONTINUE

C
      DUM = - RHOUP*POT(I,1,IA) - RHODO*POT(I,NS,IA)
     &      + RHOT(I)*(VNUC+RCH*VHART(I)+EXC)
      WDUM(I)=DUM*R(I)**2
314   CONTINUE
      if ((ivshift.eq.1).and.(iavshift.eq.IA)) then
       FPOT(:,:,IA)=FPOT(:,:,IA)+vshift
      end if

C                                        CORE TERM
      SCOR=RCZ
      IF(NUMCOR(IA).EQ.0) GO TO 201
      IADD=0
      IF(IREL.EQ.2) IADD=1
      DO 316 J=1,NUMCOR(IA)
      L2P1=2*LCOR(J,IA)+1
      AMUL1=REAL(L2P1+IADD)
      AMUL2=REAL(L2P1-IADD)
      SCOR=SCOR+AMUL1*ECOR(J,1,IA)+AMUL2*ECOR(J,2,IA)
316   CONTINUE
C
C                                        VALENCE TERM
201   SVAL=RCZ
      DO 317 IL=1,NL
      SVAL=SVAL + EMDI0(IL,1,IA)*ENY(IL,1,IA)+EMDI1(IL,1,IA)
     &       + EMDI0(IL,NS,IA)*ENY(IL,NS,IA)+EMDI1(IL,NS,IA)
317   CONTINUE
C
C                                      CORRECTION TERM
      COTE=PI4*QUAD3(NR1,R,WDUM)
C                                      ELECTROSTATIC TERM
      ELST=RCH*CHATRA(IA)*VMAD(IA)
C
      ACTE(IA)=SCOR+SVAL+COTE+ELST


C
310   CONTINUE



C
C------------------------------ DEVIATION OF THE POTENTIALS
C
      DO 320 IA=1,NA
C
      NR1=NSZRAD(IA)
      WS1=WS(IA)
      CALL RAPO(NR1,WS1,R)
C
      DO 321 IS=1,NS
      FPOT(1,IS,IA)=RCZ
      DO 322 I=2,NR1
      FPOT(I,IS,IA)=FPOT(I,IS,IA)-POT(I,IS,IA)
322   CONTINUE
321   CONTINUE
C
      DO 323 IS=1,NS
      YMAX(IS,IA)=RCZ
      IMAX(IS,IA)=0
      DO 324 I=1,NR1
      RDVR=R(I)*ABS(FPOT(I,IS,IA))
      IF(RDVR.LT.YMAX(IS,IA)) GO TO 324
      YMAX(IS,IA)=RDVR
      IMAX(IS,IA)=I
324   CONTINUE
323   CONTINUE
C
320   CONTINUE
C                                           PRINT
      WRITE(IW6,130) IITER
130   FORMAT(/' **** ATOM, ',
     &' MAX. DEVIATION(S), MESH POINT(S):',10X,'ITER=',I4)
      DO 325 IA=1,NA
      WRITE(IW6,131) IA,(YMAX(IS,IA),IMAX(IS,IA),IS=1,NS)
325   CONTINUE
131   FORMAT(5X,I5,G17.3,I6,G17.3,I6)
C
C---------------------------- MIXING OF POTENTIALS
C
      J=0
      DO 3301 IA=1,NA
      NR1=NSZRAD(IA)
      DO 3302 IS=1,NS
      DO 3303 I=2,NR1
      J=J+1
      XL(J)=POT(I,IS,IA)
3303  CONTINUE
3302  CONTINUE
3301  CONTINUE
C
      J=0
      DO 3305 IA=1,NA
      NR1=NSZRAD(IA)
      DO 3306 IS=1,NS
      DO 3307 I=2,NR1
      J=J+1
      FL(J)=FPOT(I,IS,IA)
3307  CONTINUE
3306  CONTINUE
3305  CONTINUE
C
      IF(IITER.GT.1) THEN
      DO 3311 J=1,NAM
      DXP(J,1)=DXP(J,1)-XL(J)
3311  CONTINUE
      DO 3312 J=1,NAM
      DFP(J,1)=DFP(J,1)-FL(J)
3312  CONTINUE
      CALL OM1C(NAM,NUH)
      END IF
C
      IF(IITER.LE.NITERA) THEN
      DO 332 J=1,NAM
      XN(J)=XL(J)+ALFA*FL(J)
332   CONTINUE
      ELSE
CCC      NUH1=MIN(NUH,IITER-NITERA+1)
      NUH1=MIN(NUH,IITER-1)
      CALL AMST(BETA,W0AM,NAM,NUH1)
      END IF
C
      DO 333 JH=NUH,2,-1
      DO 334 IH=JH,NUH
      VOMA(IH,JH)=VOMA(IH-1,JH-1)
334   CONTINUE
333   CONTINUE
C
      DO 3351 IH=NUH,2,-1
      DO 3361 J=1,NAM
      DXP(J,IH)=DXP(J,IH-1)
3361  CONTINUE
3351  CONTINUE
      DO 3371 J=1,NAM
      DXP(J,1)=XL(J)
3371  CONTINUE
C
      DO 3352 IH=NUH,2,-1
      DO 3362 J=1,NAM
      DFP(J,IH)=DFP(J,IH-1)
3362  CONTINUE
3352  CONTINUE
      DO 3372 J=1,NAM
      DFP(J,1)=FL(J)
3372  CONTINUE
C
      J=0
      DO 3381 IA=1,NA
      NR1=NSZRAD(IA)
      DO 3382 IS=1,NS
      DO 3383 I=2,NR1
      J=J+1
      POT(I,IS,IA)=XN(J)
3383  CONTINUE
3382  CONTINUE
3381  CONTINUE
C
C--------------------------------- POTENTIAL PARAMETERE
C                                  OF THE NEW POTENTIALS
C
C                                       LOOP OVER ATOMS
      DO 350 IA=1,NA
C                                          RADIAL MESH
      AZ1=AZ(IA)
      IREL1=MIN(IREL,1)
      IF(AZ1.LT.AZL) IREL1=0
      NR1=NSZRAD(IA)
      WS1=WS(IA)
      WSAV1=WSAV(IA)
      CALL RAPO(NR1,WS1,R)
C                                       LOOP OVER SPINS
      DO 351 IS=1,NS
C                                           POTENTIAL
      DO 352 I=1,NR1
      V(I)=POT(I,IS,IA)
352   CONTINUE
      CALL LIPO(AZ1,NR1)
C                                          LOOP OVER L
      DO 353 IL=1,NL
C                                           NEW ENY
      ISH1=ISHENY(IL,IS,IA)
      IF(ISH1.EQ.0) ENY1=ENY(IL,IS,IA)
      IF(ISH1.EQ.1.AND.EMDI0(IL,IS,IA).GT.RCZ)
     &  ENY1 = ENY(IL,IS,IA) + EMDI1(IL,IS,IA)/EMDI0(IL,IS,IA)
      IF(ISH1.EQ.2) ENY1=PPC(IL,IS,IA)
      IF(ISH1.EQ.3) ENY1=eshift
C
      L=IL-1
      DE=DMUL/SQRT(PPP(IL,IS,IA))
C                                         LINEARIZATION
      CALL RSEL(AZ1,DE,WSAV1,L,IREL1,NR1)
C
      ENY(IL,IS,IA)=ENY1
      PPC(IL,IS,IA)=PPC1
      PPD(IL,IS,IA)=PPD1
      PPQ(IL,IS,IA)=PPQ1
      PPP(IL,IS,IA)=PPP1
      DNY(IL,IS,IA)=DNY1
      FINY(IL,IS,IA)=FINY1
      FINYD(IL,IS,IA)=FINYD1
C                              STORAGE OF PHI,PHID,PHIDD
      DO 355 I=1,NR1
      PHI(I,IL,IS,IA)=FI(I)
      PHID(I,IL,IS,IA)=FID(I)
      PHIDD(I,IL,IS,IA)=FIDD(I)
355   CONTINUE
C
353   CONTINUE
351   CONTINUE
350   CONTINUE
C
C-------------------------------     OUTPUT - UNIT IW7
C
      MODPR= MOD(IITER,NFUPR)
      IF(MODPR.NE.0.AND.IITER.NE.1
     &             .AND.IITER.NE.NITER) RETURN
C
      
      call write_atoms(NA,NL,NS,NSZRAD,EF,AZ,WS,WSAV,
     & QSCR,POT,OTXTA,amgmom)

C
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
C
190   FORMAT(1X,'  --------  LSDA-FILE  -------- ')
191   FORMAT(1X,3G15.7,4X,A16,' IS=',I1)
192   FORMAT(1X,I5,   44X,A16,' IS=',I1)
193   FORMAT(1X,I5,   44X,A16)
195   FORMAT(1X,'  ----------------------------- ')
C
      RETURN
      END
C*******************
CXXX    TISK    ****
C*******************
      SUBROUTINE TISK
C
C****************************************
C  SHORT AND FULL PRINT OF THE RESULTS
C****************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNB=1)
      PARAMETER(MNA=4)
      PARAMETER(MNL=4)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNR=400)
C
      DIMENSION QT(MNA),QP(MNL,MNA),AMP(MNL,MNA),
     &          BHFC(MNA),BHFV(MNA),BHFT(MNA),
     &          DNC(MNA),DNV(MNA),DNT(MNA)
C
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),ZNCN(MNE),ZNCW(MNE),
     &            DIAM,EF,NE
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA)
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),
     &            VMAD(MNA),ACTE(MNA)
      COMMON/DATDIM/ NB,NC(MNB),NA,NL,NS
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR12,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/TEXT/ OTXTA(MNA)
C
      DATA RCZ/0.0D0/, DLIM/0.02D0/, BCONST/52.4308D0/
C
C------------------------------ SHORT PRINT IN EACH ITERATION
C
C                                ATOMIC CHARGES AND MOMENTS
      WRITE(IW6,110) IITER
110   FORMAT(/1X,'----- ATOM, CHARGE TRANSF.,',
     &  ' MAGN. MOMENT:',15X,' ITER=',I4)
      DO 311 IA=1,NA
      WRITE(IW6,111) IA,CHATRA(IA),AMGMOM(IA)
311   CONTINUE
111   FORMAT(5X,I5,2F15.5)
C                                ATOMIC MADELUNG TERMS AND
C                             CONTRIBUTIONS TO TOTAL ENERGY
      WRITE(IW6,112) IITER
112   FORMAT(/1X,'----- ATOM, MADELUNG TERM,',
     &  '  ENERGY CONTRIBUTION:',8X,' ITER=',I4)
      DO 313 IA=1,NA
      WRITE(IW6,113) IA,VMAD(IA),ACTE(IA)
313   CONTINUE
113   FORMAT(5X,I5,2X,G15.7,2X,G20.12)
C                                        POT. PARAMETERS C
      WRITE(IW6,114) IITER
114   FORMAT(/1X,'----- ATOM, SPIN, POTENTIAL PARAMETERS C:',
     &                    15X,'ITER=',I4)
      DO 314 IA=1,NA
      DO 315 IS=1,NS
      WRITE(IW6,115) IA,IS,(PPC(IL,IS,IA),IL=1,NL)
315   CONTINUE
314   CONTINUE
115   FORMAT(5X,2I5,4F15.5)
C                                 TOTAL MAGMOM AND ENERGY
      AMTOT=RCZ
      ETOT=RCZ
      DO 316 IA=1,NA
      AMTOT=AMTOT+CON(IA)*AMGMOM(IA)
      ETOT=ETOT+CON(IA)*ACTE(IA)
316   CONTINUE
      IF(NS.EQ.2) WRITE(IW6,116) AMTOT,IITER
116   FORMAT(/2X,'TOTAL MAGMOM=',F12.6,28X,'ITER=',I4)
      WRITE(IW6,117) ETOT,IITER
117   FORMAT(/2X,'TOTAL ENERGY=',G20.12,20X,'ITER=',I4)
      ETOTPA=ETOT/REAL(NB)
      WRITE(IW6,118) ETOTPA,IITER
118   FORMAT(2X,'ETOT PER ATOM=',G20.12,20X,'ITER=',I4)
C
C------------------------------- POLES OF THE SCREENED
C                                 POTENTIAL FUNCTIONS
      DO 320 IA=1,NA
      DO 321 IS=1,NS
      DO 322 IL=1,NL
      POLE=PPC(IL,IS,IA)-PPD(IL,IS,IA)/(PPQ(IL,IS,IA)-QSCR(IL))
      IF(ABS(POLE-EF).GT.DLIM) GO TO 322
      WRITE(IW6,120) POLE,IA,IS,IL
120   FORMAT(1X,'  * POLE IN E= ',G15.7,'    FOR',
     &          '  IA=',I2,'  IS=',I1,'  IL=',I1)
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
C--------------------------------- FULL PRINT IN SELECTED
C                                    ITERATIONS
      MODPR= MOD(IITER,NFUPR)
      IF(MODPR.NE.0.AND.IITER.NE.1
     &             .AND.IITER.NE.NITER) RETURN
C
      WRITE(IW6,129) IITER
129   FORMAT(/3X,'**********  FULL PRINT IN ITERATION:',
     &   ' ITER=',I4,'  **********')
C                                       CORE ENERGIES
      WRITE(IW6,130) IITER
130   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     &   13X,'- CORES:   N,   L,  NOBC,',8X,'EIGENVALUES:'/)
      DO 331 IA=1,NA
      WRITE(IW6,131) OTXTA(IA)
131   FORMAT(4X,A16)
      IF(NUMCOR(IA).EQ.0) GO TO 331
      DO 332 J=1,NUMCOR(IA)
      WRITE(IW6,132) NCOR(J,IA),LCOR(J,IA),NOBC(J,IA),
     &                     (ECOR(J,IX,IA),IX=1,2)
332   CONTINUE
132   FORMAT(20X,3I5,3X,2G15.7)
331   CONTINUE
C                                       ENERGY MOMENTS -
C                                       - L-DIAGONAL
      WRITE(IW6,135) IITER
135   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     & 6X,'- IS,  IL, ISHENY,   ENERGY MOMENTS (0,1,2):'/)
      DO 336 IA=1,NA
      WRITE(IW6,136) OTXTA(IA)
136   FORMAT(4X,A16)
      DO 337 IS=1,NS
      DO 338 IL=1,NL
      WRITE(IW6,137) IS,IL,ISHENY(IL,IS,IA),
     &     EMDI0(IL,IS,IA),EMDI1(IL,IS,IA),EMDI2(IL,IS,IA)
338   CONTINUE
337   CONTINUE
137   FORMAT(5X,3I5,3X,3G15.7)
336   CONTINUE
C                                       VALENCE CHARGES
      DO 340 IA=1,NA
      QT(IA)=RCZ
      DO 341 IL=1,NL
      QP(IL,IA)=EMDI0(IL,1,IA)+EMDI0(IL,NS,IA)
      QT(IA)=QT(IA)+QP(IL,IA)
341   CONTINUE
340   CONTINUE
C
      WRITE(IW6,140) IITER
140   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     & 5X,'  - TOTAL  AND  PARTIAL VALENCE CHARGES:'/)
      DO 342 IA=1,NA
      WRITE(IW6,142) OTXTA(IA)
142   FORMAT(4X,A16)
      WRITE(IW6,143) QT(IA),(QP(IL,IA),IL=1,NL)
143   FORMAT(1X,G15.7,3X,4G15.7)
342   CONTINUE
C                                       MAGNETIC MOMENTS
       IF(NS.EQ.2) THEN
      DO 345 IA=1,NA
      DO 346 IL=1,NL
      AMP(IL,IA)=EMDI0(IL,1,IA)-EMDI0(IL,NS,IA)
346   CONTINUE
345   CONTINUE
C
      WRITE(IW6,145) IITER
145   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     & 5X,'  - TOTAL  AND  PARTIAL MAGNETIC MOMENTS:'/)
      DO 347 IA=1,NA
      WRITE(IW6,147) OTXTA(IA)
147   FORMAT(4X,A16)
      WRITE(IW6,148) AMGMOM(IA),(AMP(IL,IA),IL=1,NL)
148   FORMAT(1X,G15.7,3X,4G15.7)
347   CONTINUE
       END IF
C                                       POTENTIAL PARAMETERS
      WRITE(IW6,150) IITER
150   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     & 15X,'  ---- POTENTIAL PARAMETERS :')
      DO 351 IA=1,NA
      DO 352 IS=1,NS
      WRITE(IW6,1501) OTXTA(IA),IS,IITER,AZ(IA),WS(IA),WSAV(IA)
1501  FORMAT(/4X,A16,15X,'SPIN=',I2,15X,'ITER=',I4/
     &  5X,'AZ=',F7.3,8X,'WS=',F10.6,8X,'WSAV=',F10.6)
      WRITE(IW6,151)(ENY(IL,IS,IA),IL=1,NL)
151   FORMAT(1X,' ENY : ',4G15.7)
      WRITE(IW6,152)(PPC(IL,IS,IA),IL=1,NL)
152   FORMAT(1X,'  C  : ',4G15.7)
      WRITE(IW6,153)(PPD(IL,IS,IA),IL=1,NL)
153   FORMAT(1X,'DELTA: ',4G15.7)
      WRITE(IW6,154)(PPQ(IL,IS,IA),IL=1,NL)
154   FORMAT(1X,'  Q  : ',4G15.7)
      WRITE(IW6,155)(PPP(IL,IS,IA),IL=1,NL)
155   FORMAT(1X,'  P  : ',4G15.7)
      WRITE(IW6,156)(DNY(IL,IS,IA),IL=1,NL)
156   FORMAT(1X,' DNY : ',4G15.7)
      WRITE(IW6,157)(FINY(IL,IS,IA),IL=1,NL)
157   FORMAT(1X,'FINY : ',4G15.7)
      WRITE(IW6,158)(FINYD(IL,IS,IA),IL=1,NL)
158   FORMAT(1X,'FINYD: ',4G15.7)
352   CONTINUE
351   CONTINUE
C                                       NUCLEAR QUANTITIES
      IF(IREL.EQ.0) THEN
C                                       HYPERFINE FIELDS
       IF(NS.EQ.2) THEN
!SMP$ DO SERIAL       
      DO 370 IA=1,NA
      BHFC(IA)=(RHOCOR(1,1,IA)-RHOCOR(1,NS,IA))*BCONST
      BHFV(IA)=(RHOVAL(1,1,IA)-RHOVAL(1,NS,IA))*BCONST
      BHFT(IA)=BHFC(IA)+BHFV(IA)
370   CONTINUE
C
      WRITE(IW6,170) IITER
170   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     & 20X,'- TOTAL, CORE, AND VALENCE HYPERFINE FIELDS (T):'/)
      DO 372 IA=1,NA
      WRITE(IW6,172) OTXTA(IA),BHFT(IA),BHFC(IA),BHFV(IA)
172   FORMAT(1X,A16,4X,G15.7,3X,G15.7,G15.7)
372   CONTINUE
       END IF
C                                       DENSITIES AT NUCLEI
!SMP$ DO SERIAL
      DO 375 IA=1,NA
      DNC(IA)=RHOCOR(1,1,IA)+RHOCOR(1,NS,IA)
      DNV(IA)=RHOVAL(1,1,IA)+RHOVAL(1,NS,IA)
      DNT(IA)=DNC(IA)+DNV(IA)
375   CONTINUE
C
      WRITE(IW6,175) IITER
175   FORMAT(/1X,' ****  ATOMS, ',45X,' ITER=',I4/
     & 20X,'- TOTAL, CORE, AND VALENCE DENSITY AT NUCLEUS:'/)
      DO 377 IA=1,NA
      WRITE(IW6,177) OTXTA(IA),DNT(IA),DNC(IA),DNV(IA)
177   FORMAT(1X,A16,1X,2G18.10,G15.7)
377   CONTINUE
C
      END IF
C
      RETURN
      END
C*******************
CXXX    PRIM3   ****
C*******************
      SUBROUTINE PRIM3(N,X,Y,F)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C---------------------------------------------------
C     CALCULATES INDEFINITE INTEGRAL F(I)
C     FOR INTEGRAND Y(I) OF VARIABLE X(I)
C     (WHERE I=1,2,... N, AND N.GE.6)
C     USING 5TH DEGREE INTERPOLATION POLYNOMIAL
C---------------------------------------------------
      DIMENSION X(N),Y(N),F(N)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/,RC1/1.0D0/,RC5/5.0D0/,
     &     RC12/12.0D0/
C
      AUXF(T0,T1,T2,T3,T4,T5,TT)=
     &    (TT-T1)*(TT-T2)*(TT-T3)*(TT-T4)*(TT-T5)/
     &   ((T0-T1)*(T0-T2)*(T0-T3)*(T0-T4)*(T0-T5))
C
      SQ15=SQRT(RC1/RC5)
      CP=RCH*(RC1+SQ15)
      CM=RCH*(RC1-SQ15)
C
      F(1)=RCZ
      DO 301 I=2,N
      I0=MAX(I-4,0)
      I0=MIN(I0,N-6)
      X1=X(I0+1)
      X2=X(I0+2)
      X3=X(I0+3)
      X4=X(I0+4)
      X5=X(I0+5)
      X6=X(I0+6)
C                             LAGRANGE INTERPOLATION
      XT1=CP*X(I-1)+CM*X(I)
      YT1=Y(I0+1)*AUXF(X1,X2,X3,X4,X5,X6,XT1)
     &   +Y(I0+2)*AUXF(X2,X3,X4,X5,X6,X1,XT1)
     &   +Y(I0+3)*AUXF(X3,X4,X5,X6,X1,X2,XT1)
     &   +Y(I0+4)*AUXF(X4,X5,X6,X1,X2,X3,XT1)
     &   +Y(I0+5)*AUXF(X5,X6,X1,X2,X3,X4,XT1)
     &   +Y(I0+6)*AUXF(X6,X1,X2,X3,X4,X5,XT1)
C
      XT2=CM*X(I-1)+CP*X(I)
      YT2=Y(I0+1)*AUXF(X1,X2,X3,X4,X5,X6,XT2)
     &   +Y(I0+2)*AUXF(X2,X3,X4,X5,X6,X1,XT2)
     &   +Y(I0+3)*AUXF(X3,X4,X5,X6,X1,X2,XT2)
     &   +Y(I0+4)*AUXF(X4,X5,X6,X1,X2,X3,XT2)
     &   +Y(I0+5)*AUXF(X5,X6,X1,X2,X3,X4,XT2)
     &   +Y(I0+6)*AUXF(X6,X1,X2,X3,X4,X5,XT2)
C                                       LOBATTO RULE
      F(I)=F(I-1)+
     &  (Y(I-1)+RC5*(YT1+YT2)+Y(I))*(X(I)-X(I-1))/RC12
301   CONTINUE
      RETURN
      END
C*******************
CXXX    QUAD3   ****
C*******************
      REAL*8 FUNCTION QUAD3(N,X,Y)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C---------------------------------------------------
C     QUADRATURE OF INTEGRAND Y(I) OF VARIABLE X(I)
C     (WHERE I=1,2,... N, AND N.GE.6)
C     USING 5TH DEGREE INTERPOLATION POLYNOMIAL
C---------------------------------------------------
      DIMENSION X(N),Y(N)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/,RC1/1.0D0/,RC5/5.0D0/,
     &     RC12/12.0D0/
C
      AUXF(T0,T1,T2,T3,T4,T5,TT)=
     &    (TT-T1)*(TT-T2)*(TT-T3)*(TT-T4)*(TT-T5)/
     &   ((T0-T1)*(T0-T2)*(T0-T3)*(T0-T4)*(T0-T5))
C
      SQ15=SQRT(RC1/RC5)
      CP=RCH*(RC1+SQ15)
      CM=RCH*(RC1-SQ15)
C
      QUAD3=RCZ
      DO 301 I=2,N
      I0=MAX(I-4,0)
      I0=MIN(I0,N-6)
      X1=X(I0+1)
      X2=X(I0+2)
      X3=X(I0+3)
      X4=X(I0+4)
      X5=X(I0+5)
      X6=X(I0+6)
C                             LAGRANGE INTERPOLATION
      XT1=CP*X(I-1)+CM*X(I)
      YT1=Y(I0+1)*AUXF(X1,X2,X3,X4,X5,X6,XT1)
     &   +Y(I0+2)*AUXF(X2,X3,X4,X5,X6,X1,XT1)
     &   +Y(I0+3)*AUXF(X3,X4,X5,X6,X1,X2,XT1)
     &   +Y(I0+4)*AUXF(X4,X5,X6,X1,X2,X3,XT1)
     &   +Y(I0+5)*AUXF(X5,X6,X1,X2,X3,X4,XT1)
     &   +Y(I0+6)*AUXF(X6,X1,X2,X3,X4,X5,XT1)
C
      XT2=CM*X(I-1)+CP*X(I)
      YT2=Y(I0+1)*AUXF(X1,X2,X3,X4,X5,X6,XT2)
     &   +Y(I0+2)*AUXF(X2,X3,X4,X5,X6,X1,XT2)
     &   +Y(I0+3)*AUXF(X3,X4,X5,X6,X1,X2,XT2)
     &   +Y(I0+4)*AUXF(X4,X5,X6,X1,X2,X3,XT2)
     &   +Y(I0+5)*AUXF(X5,X6,X1,X2,X3,X4,XT2)
     &   +Y(I0+6)*AUXF(X6,X1,X2,X3,X4,X5,XT2)
C                                       LOBATTO RULE
      QUAD3=QUAD3+
     &  (Y(I-1)+RC5*(YT1+YT2)+Y(I))*(X(I)-X(I-1))/RC12
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
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
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
      DATA RCZ/0.0D0/, H/2.0D-5/
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
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNR=400)
C
C-----------------------------------------------------------
C   LAGRANGE INTERPOLATION OF POTENTIAL FOR THE RUNGE-KUTTA
C   INTEGRATION OF THE RADIAL SCHROEDINGER EQUATION
C-----------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     NR - SIZE OF RADIAL MESH
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     NSIRK - NO. OF SUBINTERVALS (1.LE.NSIRK.LE.5)
C  OUTPUT:
C     VI(.,.) - INTERPOLATED POTENTIAL:  THE VALUES VI(J,I),
C               WHERE J=1,2,... 2*NSIRK-1, AND I=3,4,... NR,
C               REFER TO RADIAL POINTS SAMPLING UNIFORMLY
C               THE INTERVAL (R(I-1),R(I)).
C-----------------------------------------------------------
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
C
      DATA RC1/1.0D0/,RC2/2.0D0/
C
      AUXF(T0,T1,T2,T3,TT)=
     &    (TT-T1)*(TT-T2)*(TT-T3)/((T0-T1)*(T0-T2)*(T0-T3))
C
      N2=2*NSIRK
      UN2=RC1/REAL(N2)
C
      DO 310 I=3,NR
      I0=I-3
      IF(I.EQ.NR) I0=I-4
      X1=R(I0+1)
      X2=R(I0+2)
      X3=R(I0+3)
      X4=R(I0+4)
      IF(I.EQ.3) THEN
       Y1=-RC2*AZ
      ELSE
       Y1=X1*V(I0+1)
      END IF
      Y2=X2*V(I0+2)
      Y3=X3*V(I0+3)
      Y4=X4*V(I0+4)
C                             LAGRANGE INTERPOLATION
      DO 320 J=1,N2-1
      XT=(R(I-1)*REAL(N2-J)+R(I)*REAL(J))*UN2
      YT=Y1*AUXF(X1,X2,X3,X4,XT)
     &  +Y2*AUXF(X2,X3,X4,X1,XT)
     &  +Y3*AUXF(X3,X4,X1,X2,XT)
     &  +Y4*AUXF(X4,X1,X2,X3,XT)
      VI(J,I)=YT/XT
320   CONTINUE
310   CONTINUE
      RETURN
      END

C*******************
CXXX   RUKUST   ****
C*******************
      SUBROUTINE RUKUST(P,Q,H,APP,APQ,AQP,AQQ,NSI)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------------
C    ONE STEP OF RUNGE-KUTTA INTEGRATION OF A SYSTEM
C    OF TWO LINEAR DIFFERENTIAL EQUATIONS:
C        P'(X) = APP(X) * P(X) + APQ(X) * Q(X)
C        Q'(X) = AQP(X) * P(X) + AQQ(X) * Q(X)
C-----------------------------------------------------
C  ON INPUT: P,Q - INITIAL VALUES
C            H  -  INCREMENT OF X
C            NSI - NO. OF SUBINTERVALS (1.LE.NSI.LE.5)
C            APP(J), APQ(J), AQP(J), AQQ(J)
C              (J = 0, ... , 2*NSI) - VALUES OF THE
C              COEFFICIENTS FOR ARGUMENTS X SAMPLING
C              UNIFORMLY THE WHOLE STEP OF LENGTH H
C-----------------------------------------------------
C  ON OUTPUT: P,Q - FINAL VALUES
C-----------------------------------------------------
C
      DIMENSION APP(0:10),APQ(0:10),AQP(0:10),AQQ(0:10)
C
      DATA RC2/2.0D0/,RC6/6.0D0/
C
      H1=H/REAL(NSI)
      H2=H1/RC2
      H6=H1/RC6
C
      J=0
      DO 301 ISI=1,NSI
      PK1=APP(J)*P+APQ(J)*Q
      QK1=AQP(J)*P+AQQ(J)*Q
      J=J+1
      PDUM=P+H2*PK1
      QDUM=Q+H2*QK1
      PK2=APP(J)*PDUM+APQ(J)*QDUM
      QK2=AQP(J)*PDUM+AQQ(J)*QDUM
      PDUM=P+H2*PK2
      QDUM=Q+H2*QK2
      PK3=APP(J)*PDUM+APQ(J)*QDUM
      QK3=AQP(J)*PDUM+AQQ(J)*QDUM
      J=J+1
      PDUM=P+H1*PK3
      QDUM=Q+H1*QK3
      PK4=APP(J)*PDUM+APQ(J)*QDUM
      QK4=AQP(J)*PDUM+AQQ(J)*QDUM
      P=P+H6*(PK1+RC2*(PK2+PK3)+PK4)
      Q=Q+H6*(QK1+RC2*(QK2+QK3)+QK4)
301   CONTINUE
C
      RETURN
      END
C*******************
CXXX    RSEV    ****
C*******************
      SUBROUTINE RSEV(AZ,E,E0,L,IREL,NR)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNR=400)
C
C---------------------------------------------------------
C     SOLUTION OF RADIAL SCHROEDINGER EQUATION
C     FOR VALENCE ELECTRONS:
C     BOTH NON-RELATIVISTIC (IREL=0) AND
C     SCALAR-RELATIVISTIC (IREL=1) VERSION
C---------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     E - ENERGY
C     E0 - ENERGY FOR DOWNFOLDING OF SMALL COMPONENTS
C     L - ORBITAL QUANTUM NUMBER
C     IREL - RELATIVITY
C     NR - SIZE OF RADIAL MESH
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     VI(.,.) - INTERPOLATED POTENTIAL
C     NSIRK - NO. OF SUBINTERVALS FOR RUNGE-KUTTA METHOD
C  OUTPUT:
C     WG(.) - WAVE FUNCTION NORMALIZED TO UNITY
C             (LARGE COMPONENT)
C     WF(.) - FUNCTION RELATED TO RADIAL DERIVATIVE OF WG
C             (SMALL COMPONENT)
C---------------------------------------------------------
C
      DIMENSION P(MNR),Q(MNR),WRK(MNR)
      DIMENSION TAPP(0:10),TAPQ(0:10),TAQP(0:10),TAQQ(0:10),
     &          TR(0:10),TV(0:10)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
C
      DATA NRST/10/, C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/
C
C                             SET UP CONSTANTS
      AL=REAL(L)
      ALP1=REAL(L+1)
      BLAM=AL*ALP1
      TWOZ=RC2*AZ
      UC=REAL(IREL)/C274
      UCSQ=UC**2
      NSI2=2*NSIRK
      UN2=RC1/REAL(NSI2)
C                            START OF INTEGRATION
      P(1)=RCZ
      Q(1)=RCZ
C                              NON-RELATIVISTIC
      IF(IREL.EQ.0) THEN
      U0=V(2)+TWOZ/R(2)-E
      A1=-AZ/ALP1
      A2=(U0-TWOZ*A1)/(RC2*(RC2*AL+RC3))
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
C
      DO 301 I=3,NRST
C
      DO 3011 J=0,NSI2
      TR(J)=(R(I-1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3011  CONTINUE
      TV(0)=V(I-1)
      DO 3012 J=1,NSI2-1
      TV(J)=VI(J,I)
3012  CONTINUE
      TV(NSI2)=V(I)
C
      DO 3014 J=0,NSI2
      TAPP(J)=A1MB/TR(J)
3014  CONTINUE
      DO 3015 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E0-TV(J))
3015  CONTINUE
      DO 3016 J=0,NSI2
      TAQP(J)=BLAM/(TAPQ(J)*TR(J)**2)+TV(J)-E
3016  CONTINUE
      DO 3017 J=0,NSI2
      TAQQ(J)=-A1PB/TR(J)
3017  CONTINUE
C
      TH=R(I)-R(I-1)
      TP=P(I-1)
      TQ=Q(I-1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      P(I)=TP
      Q(I)=TQ
301   CONTINUE
C
222   DO 310 I=2,NRST
      RTB=R(I)**BETA
      P(I)=RTB*P(I)
      Q(I)=RTB*Q(I)
310   CONTINUE
C
C                             RUNGE-KUTTA INTEGRATION
      DO 320 I=NRST+1,NR
C
      DO 3201 J=0,NSI2
      TR(J)=(R(I-1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3201  CONTINUE
      TV(0)=V(I-1)
      DO 3202 J=1,NSI2-1
      TV(J)=VI(J,I)
3202  CONTINUE
      TV(NSI2)=V(I)
C
      DO 3204 J=0,NSI2
      TAPP(J)=RC1/TR(J)
3204  CONTINUE
      DO 3205 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E0-TV(J))
3205  CONTINUE
      DO 3206 J=0,NSI2
      TAQP(J)=BLAM/(TAPQ(J)*TR(J)**2)+TV(J)-E
3206  CONTINUE
      DO 3207 J=0,NSI2
      TAQQ(J)=-RC1/TR(J)
3207  CONTINUE
C
      TH=R(I)-R(I-1)
      TP=P(I-1)
      TQ=Q(I-1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      P(I)=TP
      Q(I)=TQ
320   CONTINUE
C
C                                  NORMALIZATION
      DO 330 I=1,NR
      WRK(I)=P(I)**2
330   CONTINUE
      SUM=QUAD3(NR,R,WRK)
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
CXXX    RSEC    ****
C*******************
      SUBROUTINE RSEC(AZ,E,THRESH,N,L,IREL,NR,NOBC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNR=400)
C
C---------------------------------------------------------
C     SOLUTION OF RADIAL SCHROEDINGER EQUATION
C     FOR CORE ELECTRONS:
C     BOTH  NON-RELATIVISTIC (IREL=0) AND
C     SCALAR-RELATIVISTIC (IREL=1) VERSION
C---------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     E - ENERGY (ESTIMATE)
C     THRESH - ABSOLUTE ACCURACY OF THE EIGENVALUE
C     N - PRINCIPAL QUANTUM NUMBER
C     L - ORBITAL QUANTUM NUMBER
C     IREL - RELATIVITY
C     NR - SIZE OF RADIAL MESH
C     NOBC - OUTER BOUNDARY CONDITION
C           (NOBC=0 - DEEP LEVEL, NOBC=1 - SHALLOW LEVEL)
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     VI(.,.) - INTERPOLATED POTENTIAL
C     NSIRK - NO. OF SUBINTERVALS FOR RUNGE-KUTTA METHOD
C  OUTPUT:
C     E - ENERGY
C     WG(.) - WAVE FUNCTION (LARGE COMPONENT)
C     WF(.) - FUNCTION RELATED TO RADIAL DERIVATIVE OF WG
C             (SMALL COMPONENT)
C---------------------------------------------------------
C
      DIMENSION AME(MNR),U(MNR),P(MNR),Q(MNR),WRK(MNR)
      DIMENSION TAPP(0:10),TAPQ(0:10),TAQP(0:10),TAQQ(0:10),
     &          TR(0:10),TV(0:10)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
C
      DATA NRST/10/,MAXIT/200/, ARGMAX/120.0D0/,
     &     DSG/1.5D0/, C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/,
     &     RCH/0.5D0/
C
C                             SET UP CONSTANTS
      AL=REAL(L)
      ALP1=REAL(L+1)
      BLAM=AL*ALP1
      TWOZ=RC2*AZ
      NUMNOD=N-L-1
      UC=REAL(IREL)/C274
      UCSQ=UC**2
      NSI2=2*NSIRK
      UN2=RC1/REAL(NSI2)
      IT=0
      IMIN=0
      IMAX=0
C
200   IT=IT+1
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
C                            CLASSICAL TURNING POINT
      I=NR+1
205   I=I-1
      IF(U(I).GT.RCZ.AND.I.GT.NRST+10) GO TO 205
      NRMA=MIN(I+1,NR)
C
C -------------------------   OUTWARD INTEGRATION
C                          START -  NON-RELATIVISTIC
      IF(IREL.EQ.0) THEN
      U0=V(2)+TWOZ/R(2)-E
      A1=-AZ/ALP1
      A2=(U0-TWOZ*A1)/(RC2*(RC2*AL+RC3))
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
C
      DO 311 I=3,NRST
C
      DO 3111 J=0,NSI2
      TR(J)=(R(I-1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3111  CONTINUE
      TV(0)=V(I-1)
      DO 3112 J=1,NSI2-1
      TV(J)=VI(J,I)
3112  CONTINUE
      TV(NSI2)=V(I)
C
      DO 3114 J=0,NSI2
      TAPP(J)=A1MB/TR(J)
3114  CONTINUE
      DO 3115 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E-TV(J))
3115  CONTINUE
      DO 3116 J=0,NSI2
      TAQP(J)=BLAM/(TAPQ(J)*TR(J)**2)+TV(J)-E
3116  CONTINUE
      DO 3117 J=0,NSI2
      TAQQ(J)=-A1PB/TR(J)
3117  CONTINUE
C
      TH=R(I)-R(I-1)
      TP=P(I-1)
      TQ=Q(I-1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      P(I)=TP
      Q(I)=TQ
311   CONTINUE
C
208   DO 312 I=2,NRST
      RTB=R(I)**BETA
      P(I)=RTB*P(I)
      Q(I)=RTB*Q(I)
312   CONTINUE
C                             RUNGE-KUTTA INTEGRATION
      DO 314 I=NRST+1,NRMA
C
      DO 3141 J=0,NSI2
      TR(J)=(R(I-1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3141  CONTINUE
      TV(0)=V(I-1)
      DO 3142 J=1,NSI2-1
      TV(J)=VI(J,I)
3142  CONTINUE
      TV(NSI2)=V(I)
C
      DO 3144 J=0,NSI2
      TAPP(J)=RC1/TR(J)
3144  CONTINUE
      DO 3145 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E-TV(J))
3145  CONTINUE
      DO 3146 J=0,NSI2
      TAQP(J)=BLAM/(TAPQ(J)*TR(J)**2)+TV(J)-E
3146  CONTINUE
      DO 3147 J=0,NSI2
      TAQQ(J)=-RC1/TR(J)
3147  CONTINUE
C
      TH=R(I)-R(I-1)
      TP=P(I-1)
      TQ=Q(I-1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      P(I)=TP
      Q(I)=TQ
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
      IF(NUNO.LT.NUMNOD) GO TO 220
      IF(NUNO.GT.NUMNOD) GO TO 230
      IF(NUNO.EQ.NUMNOD) GO TO 250
C                                   ENERGY TOO LOW
220   IF(IMIN.EQ.0) EMIN=E
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
240   IF(IMIN*IMAX.EQ.0) THEN
       ABSE=ABS(E)
       AMGDE=MAX(ABSE,RC1)/REAL(10)
       E=E+SIGDE*AMGDE
      ELSE
       E=RCH*(EMAX+EMIN)
      END IF
      GO TO 200
C                            ENERGY IN THE CORRECT WINDOW
250   PMAL=P(NRMA)
      XLDL=Q(NRMA)/PMAL
      WRK(1)=RCZ
      DO 322 I=2,NRMA
      FRAC=P(I)/(R(I)*AME(I))
      WRK(I)=P(I)**2 + UCSQ*(Q(I)**2+BLAM*FRAC**2)
322   CONTINUE
      ANIL=RCZ
      DO 323 I=2,NRMA
      ANIL=ANIL+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
323   CONTINUE
      ANIL=RCH*ANIL
C
C                    STARTING POINT FOR INWARD INTEGRATION
      IF(NOBC.EQ.0) THEN
       XLD1=-SQRT(U(NR)/AME(NR))
       DLD1=-RCH/XLD1
      END IF
      IF(NOBC.EQ.1) THEN
       VME=V(NR)-E
       CALL COBC(XLD1,DLD1,VME,R(NR),L,-1)
      END IF
C
      IF(NRMA.LT.NR) GO TO 253
      PMAR=RC1
      XLDR=XLD1
      ANIR=RCZ
      GO TO 260
C
253   DO 325 I=NRMA,NR
      WRK(I)=SQRT(U(I)*AME(I))
325   CONTINUE
      SUM=RCZ
      I=NRMA
255   I=I+1
      SUM=SUM+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
      IF(SUM.LT.ARGMAX.AND.I.LT.NR) GO TO 255
      NRII=I
C
C ------------------------------------ INWARD INTEGRATION
C                                   START
      IF(NRII.LT.NR) THEN
       XLD1=-SQRT(U(NRII)/AME(NRII))
       DLD1=-RCH/XLD1
      END IF
      P(NRII)=RC1
      Q(NRII)=XLD1
C                                RUNGE-KUTTA INTEGRATION
      DO 330 I=NRII-1,NRMA,-1
C
      DO 3301 J=0,NSI2
      TR(J)=(R(I+1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3301  CONTINUE
      TV(0)=V(I+1)
      DO 3302 J=1,NSI2-1
      TV(J)=VI(NSI2-J,I+1)
3302  CONTINUE
      TV(NSI2)=V(I)
C
      XSI=WRK(I)
      DO 3304 J=0,NSI2
      TAPP(J)=XSI+RC1/TR(J)
3304  CONTINUE
      DO 3305 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E-TV(J))
3305  CONTINUE
      DO 3306 J=0,NSI2
      TAQP(J)=BLAM/(TAPQ(J)*TR(J)**2)+TV(J)-E
3306  CONTINUE
      DO 3307 J=0,NSI2
      TAQQ(J)=XSI-RC1/TR(J)
3307  CONTINUE
C
      TH=R(I)-R(I+1)
      TP=P(I+1)
      TQ=Q(I+1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      EXFA=EXP(-XSI*TH)
      P(I)=EXFA*TP
      Q(I)=EXFA*TQ
330   CONTINUE
C
      PMAR=P(NRMA)
      XLDR=Q(NRMA)/PMAR
      DO 332 I=NRMA,NRII
      FRAC=P(I)/(R(I)*AME(I))
      WRK(I)=P(I)**2 + UCSQ*(Q(I)**2+BLAM*FRAC**2)
332   CONTINUE
      ANIR=RCZ
      DO 333 I=NRMA+1,NRII
      ANIR=ANIR+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
333   CONTINUE
      ANIR=RCH*ANIR
C
C                   ENERGY CHANGE FROM PERTURBATION THEORY
260   DE=(XLDL-XLDR)/
     &   (ANIL/PMAL**2 + (DLD1+ANIR)/PMAR**2)
      IF(ABS(DE).LT.THRESH) GO TO 280
C
      IF(DE.GT.RCZ) THEN
       IF(IMIN.EQ.0) EMIN=E
       IF(IMIN.EQ.1) EMIN=MAX(EMIN,E)
       IMIN=1
      ELSE
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
      SUM=QUAD3(NR,R,WRK)
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
CXXX    RDEC    ****
C*******************
      SUBROUTINE RDEC(AZ,E,THRESH,N,KAPPA,NR,NOBC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNR=400)
C
C---------------------------------------------------------
C     SOLUTION OF RADIAL DIRAC EQUATION
C     FOR CORE ELECTRONS (NON-MAGNETIC CASE)
C---------------------------------------------------------
C  INPUT:
C     AZ - ATOMIC NUMBER
C     E - ENERGY (ESTIMATE)
C     THRESH - ABSOLUTE ACCURACY OF THE EIGENVALUE
C     N - PRINCIPAL QUANTUM NUMBER
C     KAPPA - RELATIVISTIC QUANTUM NUMBER
C     NR - SIZE OF RADIAL MESH
C     NOBC - OUTER BOUNDARY CONDITION
C           (NOBC=0 - DEEP LEVEL, NOBC=1 - SHALLOW LEVEL)
C     R(.) - RADIAL MESH
C     V(.) - POTENTIAL
C     VI(.,.) - INTERPOLATED POTENTIAL
C     NSIRK - NO. OF SUBINTERVALS FOR RUNGE-KUTTA METHOD
C  OUTPUT:
C     E - ENERGY
C     WG(.) - LARGE COMPONENT (G)
C     WF(.) - SMALL COMPONENT (C*F)
C             OF WAVE FUNCTION NORMALIZED TO UNITY
C---------------------------------------------------------
C
      DIMENSION AME(MNR),U(MNR),P(MNR),Q(MNR),WRK(MNR)
      DIMENSION TAPP(0:10),TAPQ(0:10),TAQP(0:10),TAQQ(0:10),
     &          TR(0:10),TV(0:10)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
C
      DATA NRST/10/,MAXIT/200/, ARGMAX/120.0D0/,
     &     DSG/1.5D0/, C274/274.074D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RCH/0.5D0/
C
C                             SET UP CONSTANTS
      AK=REAL(KAPPA)
      IF(KAPPA.GT.0) L=KAPPA
      IF(KAPPA.LT.0) L=-KAPPA-1
      NUMNOD=N-L-1
      TWOZ=RC2*AZ
      UC=RC1/C274
      UCSQ=UC**2
      NSI2=2*NSIRK
      UN2=RC1/REAL(NSI2)
      IT=0
      IMIN=0
      IMAX=0
C
200   IT=IT+1
      IF(IT.GT.MAXIT) STOP ' *** RDEC - NO CONVERGENCY '
C
      DO 302 I=1,NR
      P(I)=RCZ
      Q(I)=RCZ
302   CONTINUE
      DO 303 I=2,NR
      AME(I)=RC1+UCSQ*(E-V(I))
      U(I)=V(I)-E
303   CONTINUE
C
C                                   MATCHING RADIUS
      I=NR+1
205   I=I-1
      IF(U(I).GT.RCZ.AND.I.GT.NRST+10) GO TO 205
      NRMA=MIN(I+1,NR)
C
C -------------------------   OUTWARD INTEGRATION
C                                     INITIAL CONDITION
      BETA=SQRT(AK**2-(TWOZ*UC)**2)
      P(2)=RC1
      Q(2)=(BETA+AK)/(TWOZ*UCSQ)
C                       STARTING RUNGE-KUTTA INTEGRATION
      IF(NRST.EQ.2) GO TO 208
      AKPB=AK+BETA
      AKMB=AK-BETA
C
      DO 311 I=3,NRST
C
      DO 3111 J=0,NSI2
      TR(J)=(R(I-1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3111  CONTINUE
      TV(0)=V(I-1)
      DO 3112 J=1,NSI2-1
      TV(J)=VI(J,I)
3112  CONTINUE
      TV(NSI2)=V(I)
C
      DO 3114 J=0,NSI2
      TAPP(J)=-AKPB/TR(J)
3114  CONTINUE
      DO 3115 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E-TV(J))
3115  CONTINUE
      DO 3116 J=0,NSI2
      TAQP(J)=TV(J)-E
3116  CONTINUE
      DO 3117 J=0,NSI2
      TAQQ(J)=AKMB/TR(J)
3117  CONTINUE
C
      TH=R(I)-R(I-1)
      TP=P(I-1)
      TQ=Q(I-1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      P(I)=TP
      Q(I)=TQ
311   CONTINUE
C
208   DO 312 I=2,NRST
      RTB=R(I)**BETA
      P(I)=RTB*P(I)
      Q(I)=RTB*Q(I)
312   CONTINUE
C                             RUNGE-KUTTA INTEGRATION
      DO 314 I=NRST+1,NRMA
C
      DO 3141 J=0,NSI2
      TR(J)=(R(I-1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3141  CONTINUE
      TV(0)=V(I-1)
      DO 3142 J=1,NSI2-1
      TV(J)=VI(J,I)
3142  CONTINUE
      TV(NSI2)=V(I)
C
      DO 3144 J=0,NSI2
      TAPP(J)=-AK/TR(J)
3144  CONTINUE
      DO 3145 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E-TV(J))
3145  CONTINUE
      DO 3146 J=0,NSI2
      TAQP(J)=TV(J)-E
3146  CONTINUE
      DO 3147 J=0,NSI2
      TAQQ(J)=AK/TR(J)
3147  CONTINUE
C
      TH=R(I)-R(I-1)
      TP=P(I-1)
      TQ=Q(I-1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      P(I)=TP
      Q(I)=TQ
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
      IF(NUNO.LT.NUMNOD) GO TO 220
      IF(NUNO.GT.NUMNOD) GO TO 230
      IF(NUNO.EQ.NUMNOD) GO TO 250
C                                   ENERGY TOO LOW
220   IF(IMIN.EQ.0) EMIN=E
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
240   IF(IMIN*IMAX.EQ.0) THEN
       ABSE=ABS(E)
       AMGDE=MAX(ABSE,RC1)/REAL(10)
       E=E+SIGDE*AMGDE
      ELSE
       E=RCH*(EMAX+EMIN)
      END IF
      GO TO 200
C                            ENERGY IN THE CORRECT WINDOW
250   PMAL=P(NRMA)
      XLDL=Q(NRMA)/PMAL
      DO 322 I=1,NRMA
      WRK(I)=P(I)**2+UCSQ*Q(I)**2
322   CONTINUE
      ANIL=RCZ
      DO 323 I=2,NRMA
      ANIL=ANIL+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
323   CONTINUE
      ANIL=RCH*ANIL
C
C                    STARTING POINT FOR INWARD INTEGRATION
      IF(NOBC.EQ.0) THEN
       XLD1=-SQRT(U(NR)/AME(NR))
       DLD1=-RCH/XLD1
      END IF
      IF(NOBC.EQ.1) THEN
       VME=V(NR)-E
       CALL COBC(XLD1,DLD1,VME,R(NR),L,KAPPA)
      END IF
C
      IF(NRMA.LT.NR) GO TO 253
      PMAR=RC1
      XLDR=XLD1
      ANIR=RCZ
      GO TO 260
C
253   DO 325 I=NRMA,NR
      WRK(I)=SQRT(U(I)*AME(I))
325   CONTINUE
      SUM=RCZ
      I=NRMA
255   I=I+1
      SUM=SUM+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
      IF(SUM.LT.ARGMAX.AND.I.LT.NR) GO TO 255
      NRII=I
C
C ------------------------------------ INWARD INTEGRATION
C                                   START
      IF(NRII.LT.NR) THEN
       XLD1=-SQRT(U(NRII)/AME(NRII))
       DLD1=-RCH/XLD1
      END IF
      P(NRII)=RC1
      Q(NRII)=XLD1
C                                RUNGE-KUTTA INTEGRATION
      DO 330 I=NRII-1,NRMA,-1
C
      DO 3301 J=0,NSI2
      TR(J)=(R(I+1)*REAL(NSI2-J)+R(I)*REAL(J))*UN2
3301  CONTINUE
      TV(0)=V(I+1)
      DO 3302 J=1,NSI2-1
      TV(J)=VI(NSI2-J,I+1)
3302  CONTINUE
      TV(NSI2)=V(I)
C
      XSI=WRK(I)
      DO 3304 J=0,NSI2
      TAPP(J)=XSI-AK/TR(J)
3304  CONTINUE
      DO 3305 J=0,NSI2
      TAPQ(J)=RC1+UCSQ*(E-TV(J))
3305  CONTINUE
      DO 3306 J=0,NSI2
      TAQP(J)=TV(J)-E
3306  CONTINUE
      DO 3307 J=0,NSI2
      TAQQ(J)=XSI+AK/TR(J)
3307  CONTINUE
C
      TH=R(I)-R(I+1)
      TP=P(I+1)
      TQ=Q(I+1)
      CALL RUKUST(TP,TQ,TH,TAPP,TAPQ,TAQP,TAQQ,NSIRK)
      EXFA=EXP(-XSI*TH)
      P(I)=EXFA*TP
      Q(I)=EXFA*TQ
330   CONTINUE
C
      PMAR=P(NRMA)
      XLDR=Q(NRMA)/PMAR
      DO 332 I=NRMA,NRII
      WRK(I)=P(I)**2+UCSQ*Q(I)**2
332   CONTINUE
      ANIR=RCZ
      DO 333 I=NRMA+1,NRII
      ANIR=ANIR+(WRK(I-1)+WRK(I))*(R(I)-R(I-1))
333   CONTINUE
      ANIR=RCH*ANIR
C
C                   ENERGY CHANGE FROM PERTURBATION THEORY
260   DE=(XLDL-XLDR)/
     &   (ANIL/PMAL**2 + (DLD1+ANIR)/PMAR**2)
      IF(ABS(DE).LT.THRESH) GO TO 280
C
      IF(DE.GT.RCZ) THEN
       IF(IMIN.EQ.0) EMIN=E
       IF(IMIN.EQ.1) EMIN=MAX(EMIN,E)
       IMIN=1
      ELSE
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
290   DO 380 I=1,NR
      WRK(I)=P(I)**2+UCSQ*Q(I)**2
380   CONTINUE
      SUM=QUAD3(NR,R,WRK)
      CNORM=SQRT(SUM)
C
      WG(1)=RCZ
      WF(1)=RCZ
      DO 385 I=2,NR
      DUM=CNORM*R(I)
      WG(I)=P(I)/DUM
      WF(I)=Q(I)/DUM
385   CONTINUE
C
      RETURN
      END
C*******************
CXXX    COBC    ****
C*******************
      SUBROUTINE COBC(XLD,DLD,VME,R,L,KAPPA)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------------
C       OUTER BOUNDARY CONDITION FOR CORE ORBITALS
C       BASED ON SPHERICAL HANKEL FUNCTIONS
C----------------------------------------------------
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/,
     &     RC6/6.0D0/,RC15/15.0D0/
C
      IF(VME.LE.RCZ) STOP ' ***  COBC: VME.LE.0  *** '
      IF(L.LT.0.OR.L.GT.3) STOP ' ***  COBC: WRONG L  *** '
C
      XSI=SQRT(VME)
      T=XSI*R
      U=RC1/T
C
      IF(L.EQ.0) THEN
       Y=RC1
       YP=RCZ
      ELSE IF(L.EQ.1) THEN
       Y=RC1+U
       YP=-U**2
      ELSE IF(L.EQ.2) THEN
       Y=RC1+RC3*U*(RC1+U)
       YP=-U**2*(RC3+RC6*U)
      ELSE IF(L.EQ.3) THEN
       Y=RC1+U*(RC6+RC15*U*(RC1+U))
       YP=-U**2*(RC6+RC15*U*(RC2+RC3*U))
      END IF
C
      W=YP/Y
      AK=REAL(KAPPA)
      BL=REAL(L)*REAL(L+1)
      XLD=XSI*(-RC1+AK*U+W)
      H=-RC1+W*(RC1+T*(RC2-W))+BL*U
      DLD=-H/(RC2*XSI)
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
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNR=400)
C
C---------------------------------------------------------
C     LINEARIZATION OF ENERGY DEPENDENCE
C     OF SOLUTION OF RADIAL SCHROEDINGER EQUATION
C     FOR VALENCE ELECTRONS:
C      BOTH  NON-RELATIVISTIC (IREL=0) AND
C      SCALAR-RELATIVISTIC (IREL=1) VERSION
C---------------------------------------------------------
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
C     VI(.,.) - INTERPOLATED POTENTIAL
C     NSIRK - NO. OF SUBINTERVALS FOR RUNGE-KUTTA METHOD
C     ENY - ENERGY VALUE
C  OUTPUT:
C     PPC, ..., FINYD - POTENTIAL PARAMETERS
C     FI(.),FID(.),FIDD(.) - WAVEFUNCTION (NORMALIZED TO
C                  UNITY) AND ITS TWO ENERGY DERIVATIVES
C---------------------------------------------------------
C
      DIMENSION WB(MNR,-2:2),WDU(MNR)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/LIN/ ENY,PPC,PPD,PPQ,PPP,DNY,FINY,FINYD,
     &            FI(MNR),FID(MNR),FIDD(MNR)
C
      DATA RC1/1.0D0/,RC2/2.0D0/,RC8/8.0D0/,
     &     RC12/12.0D0/,RC16/16.0D0/,RC30/30.0D0/
C
      AL=REAL(L)
      ALP1=AL+RC1
      WS=R(NR)
      FAK=(WS/WSA)**(2*L+1)
      E0=ENY
C                                   LOOP OVER 5 ENERGIES
      DO 305 IE=-2,2
      E=ENY+REAL(IE)*EH
      CALL RSEV(AZ,E,E0,L,IREL,NR)
      IF(IE.EQ.0) DNY=WS*WF(NR)/WG(NR)
      DO 306 I=1,NR
      WB(I,IE)=WG(I)
306   CONTINUE
305   CONTINUE
C                                     PHI, PHI-DOT, ...
      H1=RC12*EH
      H2=H1*EH
C
      DO 310 I=1,NR
      A0=WB(I,0)
      S1=WB(I,-1)+WB(I,1)
      S2=WB(I,-2)+WB(I,2)
      R1=WB(I,1)-WB(I,-1)
      R2=WB(I,2)-WB(I,-2)
C
      FI(I)=A0
      FID(I)=(RC8*R1-R2)/H1
      FIDD(I)=(RC16*S1-S2-RC30*A0)/H2
310   CONTINUE
C                                  POTENTIAL PARAMETERS
C                                    FINY,FINYD,DNYD
      FINY=FI(NR)
      FINYD=FID(NR)
      DNYD=DNY-RC1/(WS*FINY*FINYD)
C                                     PPC,PPD,PPQ
      AJM=DNYD+ALP1
C
      PPC=ENY - FINY*(DNY+ALP1)/(FINYD*AJM)
      PPD=FAK/(RC2*WS*(FINYD*AJM)**2)
      PPQ=FAK*(DNYD-AL)/(RC2*(AL+ALP1)*AJM)
C                                            PPP
!SMP$ DO SERIAL
      DO 320 I=1,NR
      WDU(I)=(R(I)*FID(I))**2
320   CONTINUE
      PPP=QUAD3(NR,R,WDU)
C
      RETURN
      END
C*******************
CXXX    HAPO    ****
C*******************
      SUBROUTINE HAPO(N,R,RHO,V)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNR=400)
C
C------------------------------------------------
C     HARTREE POTENTIAL INSIDE ONE ATOMIC SPHERE
C------------------------------------------------
C
      DIMENSION R(N),RHO(N),V(N),W1(MNR),W2(MNR)
C
      DATA RC1/1.0D0/,RC4/4.0D0/,RC8/8.0D0/
C
      PI=RC4*ATAN(RC1)
      PI8=RC8*PI
C
!SMP$ DO SERIAL
      DO 300 I=1,N
      W1(I)=R(I)*RHO(I)*PI8
300   CONTINUE
      CALL PRIM3(N,R,W1,W2)
C
      D=W2(N)
      DO 301 I=1,N
      V(I)=D-W2(I)
301   CONTINUE
C
      DO 302 I=1,N
      W1(I)=W1(I)*R(I)
302   CONTINUE
      CALL PRIM3(N,R,W1,W2)
C
      DO 303 I=2,N
      V(I)=V(I)+W2(I)/R(I)
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
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
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
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
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
C*******************
CXXX   XCCAPZ   ****
C*******************
      SUBROUTINE XCCAPZ(RHO1,RHO2,VXC1,VXC2,EXC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------------------------
C    XC-POTENTIAL AND ENERGY BASED ON WORK OF CEPERLEY AND ALDER
C         (PHYS. REV. LETT. 45 (1980) 566)
C      AS PARAMETRIZED BY PERDEW AND ZUNGER
C         (PHYS. REV. B 23 (1981) 5048)
C----------------------------------------------------------------
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/,
     &     RC4/4.0D0/
C
      DATA C238/0.2387324146D0/, C916/0.9163305866D0/
C
      DATA GAMP/-0.2846D0/, BETP1/1.0529D0/, BETP2/0.3334D0/,
     &     AP/0.0622D0/,BP/-0.096D0/,CP/0.0040D0/,DP/-0.0232D0/
C
      DATA GAMF/-0.1686D0/, BETF1/1.3981D0/, BETF2/0.2611D0/,
     &     AF/0.0311D0/,BF/-0.0538D0/,CF/0.0014D0/,DF/-0.0096D0/
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
C
      CALL AUXPZ(R,GAMP,BETP1,BETP2,AP,BP,CP,DP,EPC,RDREPC)
      CALL AUXPZ(R,GAMF,BETF1,BETF2,AF,BF,CF,DF,EFC,RDREFC)
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
CXXX   AUXPZ    ****
C*******************
      SUBROUTINE AUXPZ(RS,GAM,BET1,BET2,A,B,C,D,EC,RDREC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------------
C    AUXILIARY FUNCTION FOR PERDEW-ZUNGER CORRELATION
C    ENERGY EC DEPENDENT ON DENSITY VARIABLE RS.
C    RDREC MEANS RS TIMES RS-DERIVATIVE OF EC.
C-----------------------------------------------------
C
      DATA RC1/1.0D0/,RCH/0.5D0/
C
       IF(RS.GE.RC1) THEN
      CLEN1=BET1*SQRT(RS)
      CLEN2=BET2*RS
      SUMA=RC1+CLEN1+CLEN2
      EC=GAM/SUMA
      RDREC=-EC*(RCH*CLEN1+CLEN2)/SUMA
       ELSE
      GOL=LOG(RS)
      EC=(A+C*RS)*GOL+B+D*RS
      RDREC=A+C*RS*GOL+(C+D)*RS
       END IF
C
      RETURN
      END
C*******************
CXXX   XCVWN    ****
C*******************
      SUBROUTINE XCVWN(RHO1,RHO2,VXC1,VXC2,EXC)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------------------------
C      XC-POTENTIAL AND ENERGY ACCORDING TO VOSKO, WILK, NUSAIR
C         (CAN. J. PHYS. 58 (1980) 1200)
C----------------------------------------------------------------
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/,
     &     RC4/4.0D0/
C
      DATA C238/0.2387324146D0/, C916/0.9163305866D0/
C
      DATA AP/0.0621814D0/, XP0/-0.10498D0/, BP/3.72744D0/,
     &     CP/12.9352D0/
C
      DATA AF/0.0310907D0/, XF0/-0.32500D0/, BF/7.06042D0/,
     &     CF/18.0578D0/
C
      DATA AA/-0.0197516D0/, XA0/-0.0047584D0/, BA/1.13107D0/,
     &     CA/13.0045D0/
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
      S3=S**3
      S4=S*S3
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
C
      CALL AUXVWN(R,AP,XP0,BP,CP,EPC,RDREPC)
      CALL AUXVWN(R,AF,XF0,BF,CF,EFC,RDREFC)
      CALL AUXVWN(R,AA,XA0,BA,CA,ALT,RDRALT)
C
      BRA=ALT*(RC1-S4)+(EFC-EPC)*S4
      RDRBRA=RDRALT*(RC1-S4)+(RDREFC-RDREPC)*S4
      DSBRA=RC4*S3*(-ALT+EFC-EPC)
C
      EC=EPC+BRA*FS
      RDREC=RDREPC+RDRBRA*FS
      DSEC=DSBRA*FS+BRA*DSFS
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
CXXX   AUXVWN   ****
C*******************
      SUBROUTINE AUXVWN(RS,A,X0,B,C,Y,RDRY)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-------------------------------------------------
C    AUXILIARY FUNCTION FOR VOSKO-WILK-NUSAIR
C    CORRELATION ENERGIES
C         (EQS. (4.3) AND (4.4) OF THEIR ARTICLE)
C    RS  ............... DENSITY VARIABLE
C    A, X0, B, C  ...... CONSTANTS
C    Y  ................ DEFINED BY EQ. (4.4)
C    RDRY  ............. DEFINED BY EQ. (4.3)
C-------------------------------------------------
C
      DATA RC2/2.0D0/,RC4/4.0D0/
C
      X=SQRT(RS)
      WX=X*(X+B)+C
      WX0=X0*(X0+B)+C
      XMX0=X-X0
      Q=SQRT(RC4*C-B**2)
C
      ARG1=RS/WX
      COEF2=-B*X0/WX0
      ARG2=XMX0**2/WX
      COEF3=RC2*B*(C-X0**2)/(Q*WX0)
      ARG3=Q/(RC2*X+B)
C
      Y=A*(LOG(ARG1)+COEF2*LOG(ARG2)+COEF3*ATAN(ARG3))
C
      RDRY=A*(C-B*X0*X/XMX0)/WX
C
      RETURN
      END
C*******************
CXXX    OM1C    ****
C*******************
      SUBROUTINE OM1C(NAM,NUH)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNA=4)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
C*************************************************
C     THE FIRST COLUMN OF OVERLAP MATRIX FOR
C     THE ANDERSON MIXING PROCEDURE
C*************************************************
C     NAM - NUMBER OF VARIABLES
C     NUH - NUMBER OF PREVIOUS VECTORS
C-------------------------------------------------
C   ON INPUT:
C     SP(I), I=1,...,NAM,  -  WEIGHTS FOR THE
C                             SCALAR PRODUCT
C     DFP(I,IH), I=1,...,NAM, IH=1,...,NUH,  -
C              - DIFFERENCES OF PREVIOUS VECTORS F
C
C   ON OUTPUT:
C     VOMA(IH,1),  IH=1,...,NUH,  -
C                 - 1ST COLUMN OF OVERLAP MATRIX
C*************************************************
C
      DIMENSION VECT(MNAM)
C
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
C
      DATA RCZ/0.0D0/
C
!SMP$ DO SERIAL
      DO 310 IH=1,MNUH
      VOMA(IH,1)=RCZ
310   CONTINUE
!SMP$ DO SERIAL
      DO 312 I=1,NAM
      VECT(I) = SP(I) * DFP(I,1)
312   CONTINUE
      DO 320 IH=1,NUH
      DO 325 I=1,NAM
      VOMA(IH,1) = VOMA(IH,1) + DFP(I,IH) * VECT(I)
325   CONTINUE
320   CONTINUE
C
      RETURN
      END
C*******************
CXXX    AMST    ****
C*******************
      SUBROUTINE AMST(BETA,W0,NAM,NUH)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNA=4)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
C*************************************************
C     ONE STEP OF THE ANDERSON MIXING PROCEDURE
C     TO SOLVE NON-LINEAR EQUATIONS  F(X)=0
C*************************************************
C    BETA - MIXING PARAMETER
C      W0 - A SMALL QUANTITY TO IMPROVE STABILITY
C     NAM - NUMBER OF VARIABLES
C     NUH - NUMBER OF PREVIOUS VECTORS (NUH.GE.2)
C-------------------------------------------------
C   ON INPUT:
C     DXP(I,IH), I=1,...,NAM, IH=1,...,NUH,  -
C              - DIFFERENCES OF PREVIOUS VECTORS X
C     DFP(I,IH), I=1,...,NAM, IH=1,...,NUH,  -
C              - DIFFERENCES OF PREVIOUS VECTORS F
C     SP(I), I=1,...,NAM,  -  WEIGHTS FOR THE
C                             SCALAR PRODUCT
C     XL(I), I=1,...,NAM,  -  THE LAST VECTOR X
C     FL(I), I=1,...,NAM,  -  THE LAST VECTOR F
C     VOMA(IH,JH),  IH=1,...,NUH, JH=1,...,IH, -
C              - LOWER TRIANGLE OF OVERLAP MATRIX
C
C   ON OUTPUT:
C     XN(I), I=1,...,NAM,   -  THE NEW VECTOR X
C*************************************************
C
      DIMENSION VECT(MNAM),A(MNUH,MNUH),T(MNUH,1),
     &          WORK(MNUH)
C
      COMMON/AMXP/ DXP(MNAM,MNUH)
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      DO 301 IH=1,MNUH
      T(IH,1)=RCZ
301   CONTINUE
      DO 302 I=1,MNAM
      XN(I)=RCZ
302   CONTINUE
C                                    OVERLAP MATRIX
      DO 310 JH=1,NUH
      DO 311 IH=JH,NUH
      A(IH,JH)=VOMA(IH,JH)
311   CONTINUE
310   CONTINUE
      DUM=RC1+W0**2
      DO 314 IH=1,NUH
      A(IH,IH)=DUM*A(IH,IH)
314   CONTINUE
C                                        R.H.S.
!SMP$ DO SERIAL
      DO 320 I=1,NAM
      VECT(I) = SP(I) * FL(I)
320   CONTINUE
      DO 322 IH=1,NUH
      DO 324 I=1,NAM
      T(IH,1) = T(IH,1) + DFP(I,IH) * VECT(I)
324   CONTINUE
322   CONTINUE
C                            SOLUTION OF LINEAR SYSTEM
      CALL SOPO(A,MNUH,NUH,T,MNUH,1,WORK)
C                                    NEW VECTOR X
      DO 340 IH=1,NUH
      BT=BETA*T(IH,1)
      DO 341 I=1,NAM
      XN(I) = XN(I) + DFP(I,IH)*BT
341   CONTINUE
340   CONTINUE
      DO 342 IH=1,NUH
      DO 343 I=1,NAM
      XN(I) = XN(I) + DXP(I,IH)*T(IH,1)
343   CONTINUE
342   CONTINUE
      DO 345 I=1,NAM
      XN(I) = - XN(I) + BETA*FL(I)
345   CONTINUE
      DO 346 I=1,NAM
      XN(I) = XN(I) + XL(I)
346   CONTINUE
C
      RETURN
      END
