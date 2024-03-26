C*******************
CXXX    IES     ****
C*******************
      PROGRAM IES
C
C*******************************************************
C      PROGRAM FOR SELFCONSISTENT CALCULATION OF
C        INTERFACE ELECTRONIC STRUCTURE
C     BY MEANS OF THE TB-LMTO-SGF-CPA METHOD
C*******************************************************
C      VERSION  DECEMBER 1999
C*******************************************************
C   BASIC FEATURES OF THE PROGRAM ARE:
C     - ELECTRONIC STRUCTURE OF ORDERED AND DISORDERED
C       SURFACES AND INTERFACES OF PURE METALS AND
C       DISORDERED ALLOYS
C     - CASES: FCC(001), FCC(111), FCC(110), HCP(0001),
C              BCC(110), BCC(001) AND SIMILAR
C     - GEOMETRY: INTERFACE NORMAL LIES IN Z AXIS,
C       FROZEN RELAXATIONS ARE ALLOWED
C     - NON-RELATIVISTIC AND SCALAR-RELATIVISTIC VERSION
C       (OPTIONALLY WITH FULLY RELATIVISTIC CORES)
C     - NON-MAGNETIC AND SPIN-POLARIZED
C       (WITH COLLINEAR SPIN STRUCTURES ONLY)
C     - SUBSTITUTIONAL CHEMICAL DISORDER
C     - ATOMIC SPHERE APPROXIMATION BUT WITH DIPOLE
C       MOMENTS PERPENDICULAR TO THE INTERFACE INCLUDED
C       INTO INTERATOMIC MADELUNG TERMS
C     - SPD AND SPDF CASES
C     - EVALUATION OF TOTAL ENERGY
C*******************************************************
C   INPUT:   IR1  (= 1)  -  GEOMETRY OF THE INTERFACE
C            IR2  (= 2)  -  CHEMICAL OCCUPATIONS
C                           OF THE SITES
C            IR3  (= 3)  -  LSDA FILE (POT. PARAMETERS
C                           AND 1-EL. POTENTIALS)
C            IR4  (= 4)  -  COHERENT POTENTIAL FUNCTIONS
C            IR5  (= 5)  -  CONTROL PARAMETERS
C            IR11 (=11)  -  1ST SUBSTRATE POT. PARAMETERS,
C                           COMPOSITION, AND FERMI ENERGY
C            IR13 (=13)  -  NORMALIZED COMPLEX NODES AND
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
C    NP ...... NUMBER OF PRINCIPAL LAYERS
C              IN THE INTERMEDIATE REGION
C   MNP ...... PARAMETER: MAX. NP
C   MNPP1..... PARAMETER: MAX. NP + 1
C    NB ...... SIZE OF BASIS IN 1 PRINCIPAL LAYER
C              (NUMBER OF ATOMIC SITES IN THE BASIS)
C   MNB ...... PARAMETER: MAX. NB
C    NG ...... TOTAL NUMBER OF SITES IN THE
C              INTERMEDIATE REGION (NG=NP*NB)
C   MNG ...... PARAMETER: MAX. NG
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
C    NBZ ..... NUMBER OF K||-POINTS IN THE
C              IRREDUCIBLE PART OF 2D-BZ
C   MNBZ ..... PARAMETER: MAX. NBZ
C    NAB ..... NUMBER OF CHEMICALLY DIFFERENT
C              ATOMS IN THE 1ST (BULK) SUBSTRATE
C   MNAB ..... PARAMETER: MAX. NAB
C    NAV ..... NUMBER OF CHEMICALLY DIFFERENT
C              ATOMS IN THE 2ND (VACUUM) SUBSTRATE
C   MNAV ..... PARAMETER: MAX. NAV
C    NSB ..... NUMBER OF SPIN ORIENTATIONS
C              OF THE 1ST (BULK) SUBSTRATE
C   MNSB ..... PARAMETER: MAX. NSB
C    NSV ..... NUMBER OF SPIN ORIENTATIONS
C              OF THE 2ND (VACUUM) SUBSTRATE
C   MNSV ..... PARAMETER: MAX. NSV
C    NLSQ .... DIMENSION OF THE ON-SITE BLOCKS
C              OF GREEN'S FUNCTION (NLSQ=NL**2)
C   MNLSQ .... PARAMETER: MAX. NLSQ
C   MLMAX .... PARAMETER: MAX. ANGULAR QUANTUM NUMBER
C              OF SPHERICAL HARMONICS IN CANONICAL
C              STRUCTURE CONSTANTS ( MLMAX=2*(MNL-1) )
C   MHARM .... PARAMETER: MAX. NUMBER OF SPHERICAL
C              HARMONICS IN CANONICAL STRUCTURE
C              CONSTANTS ( MHARM=(MLMAX+1)**2 )
C    NBLSQ ... DIMENSION OF THE ON-LAYER BLOCKS OF
C              THE K||-DEPENDENT GREEN'S FUNCTION
C              (NBLSQ=NB*NLSQ)
C   MNBLSQ ... PARAMETER: MAX. NBLSQ
C   MNCL ..... PARAMETER: MAX. SIZE OF SMALL CLUSTERS
C              FOR CALCULATION OF TB-STRUCTURE CONSTANTS
C   MTBCL .... PARAMETER: MAX. DIMENSION OF MATRICES
C              FOR CALCULATION OF TB-STRUCTURE CONSTANTS
C              (MTBCL=MNCL*MNLSQ)
C   MPAIR .... PARAMETER: MAX. NUMBER OF GF-ELEMENTS
C              FOR THE DIPOLE MOMENTS (MPAIR=(MNL-1)**2)
C    MNR ..... PARAMETER: MAX. SIZE OF RADIAL MESH
C   MNAM ..... PARAMETER: MAX. DIMENSION OF VECTORS FOR
C              ANDERSON MIXING (MNAM=MNA*MNS*MNR)
C   MNUH ..... PARAMETER: MAX. NUMBER OF PREVIOUS
C              VECTORS FOR ANDERSON MIXING
C******************************************************
C
      PARAMETER(MNP=30) 
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNAB=20)
      PARAMETER(MNAV=1)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=45)
      PARAMETER(MPAIR=(MNL-1)**2)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
C*******************************************************
C             ALL COMMON BLOCKS  :
C
C-----------------------------------------  GEOMETRY
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB,0:MNPP1),
     &            JDRP(MNCL,MNB,0:MNPP1),JBVA(MNCL,MNB,0:MNPP1),
     &            JTRA(2,MNCL,MNB,0:MNPP1),NSCL(MNB,0:MNPP1)
      COMMON/BSCR/ BSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPB(MNCL,MNB),JBVAB(MNCL,MNB),
     &             JTRAB(2,MNCL,MNB),NSCLB(MNB)
      COMMON/VSCR/ VSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPV(MNCL,MNB),JBVAV(MNCL,MNB),
     &             JTRAV(2,MNCL,MNB),NSCLV(MNB)
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/BSCK/ ZBSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/VSCK/ ZVSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/AMAD/ AMCMM(MNG,MNG),AMCMD(MNG,MNG),
     &             AMCDM(MNG,MNG),AMCDD(MNG,MNG),
     &             BARCM(MNG),BARCD(MNG)
C
C-------------------------------------  GREEN'S FUNCTIONS
C
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSB)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSV)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CAGF/ ZCAGF(MNLSQ,MNLSQ,MNA,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
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
      COMMON/BPOPA/ BENY(MNL,MNSB,MNAB),BPPC(MNL,MNSB,MNAB),
     &              BPPD(MNL,MNSB,MNAB),BPPQ(MNL,MNSB,MNAB),
     &              BPPP(MNL,MNSB,MNAB),BDNY(MNL,MNSB,MNAB),
     &              BFINY(MNL,MNSB,MNAB),BFINYD(MNL,MNSB,MNAB)
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA),
     &             EMOF00(MPAIR,MNS,MNA),EMOF10(MPAIR,MNS,MNA),
     &             EMOF01(MPAIR,MNS,MNA),EMOF11(MPAIR,MNS,MNA),
     &             EMOF20(MPAIR,MNS,MNA),EMOF02(MPAIR,MNS,MNA)
      COMMON/DIBA/ DBA
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),DIPMOM(MNA),
     &            VMAD(MNA),DMAD(MNA),ACTE(MNA)
C
C-------------------------------------  GENERAL
C
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
C
C-------------------------------------  SPECIAL
C
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/AMXP/ DXP(MNAM,MNUH)
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
      COMMON/TEXT/ OTXTA(MNA)     
C
C******************************************************
C
C                              INITIALIZATION
      CALL PINI
C                             READING OF ALL DATA
      CALL RALL
C                             STARTING CALCULATIONS
      CALL BEGI
C      CALL BEGI(1)      
C                   REAL SPACE TB-STRUCTURE CONSTANTS
      CALL TBRI
      CALL TBRB
      CALL TBRV
C                             GENERATING K||-MESH
      CALL GIBZ
C                       BLOCH TRANSFORM OF TB-CONSTANTS
      CALL TBK
C                           MADELUNG CONSTANTS
      CALL MACO
C                          START OF LDA ITERATIONS
      CALL SLDA
C                           GAMMA OF BULK SUBSTRATE
      CALL GABU
C                           GAMMA OF VACUUM SUBSTRATE
      IF(IVAC.EQ.1) CALL GAVA
C
C----------------------------------  ITERATION LOOP
      DO 310 IITER=1,NITER
C                               UPDATE OF VACUUM SGF
      IF(IVAC.EQ.0) CALL SGFVA
C                            COHERENT POTENTIAL FUNCTION
      CALL SCPF
C                          RECURSION PARTITIONING
      CALL REPA
C                               CPA STEP
      CALL CPAIT
C                             ENERGY MOMENTS
      CALL CMOM
C                              CORE DENSITIES
      CALL CODE
C                              VALENCE DENSITIES
      CALL VADE
C                              DIPOLE MOMENTS
      CALL DIMO
C                              NEW POTENTIALS
      CALL NEPO
C                              PRINT OF RESULTS
      CALL TISK
C------------------------------  END OF ITERATION LOOP
310   CONTINUE
C

      
      WRITE(IW6,191)
191   FORMAT(/'        * * *   END OF IES   * * * '/)
C
      STOP
      END

      subroutine write_atoms(NA,NL,NS,NSZRAD,EF,AZ,WS,WSAV,
     & QSCR,POT,OTXTA,DBA,IVAC)
      implicit none
      integer,PARAMETER::MNA=100
      integer,PARAMETER::MNL=3
      integer,PARAMETER::MNS=2
      integer,PARAMETER::MNR=400
      integer,PARAMETER::MNBZ=700
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      
      integer:: LCOR,NCOR,NUMCOR,NOBC,IVAC
      integer:: NSYM,NK,INVE,NBZ,IREL,IVXC,ISHENY
      real(kind=8)::ECOR,ENY,PPC,PPD,PPP,FINY,DNY,PPQ,FINYD
      real(kind=8)::THRESH,DBA,AKBZ,WKBZ
      
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
      
      do IA=1,NA
         read(OTXTA(IA),*) cwork
         open (unit=atomf, file='atoms/'//trim(cwork),
     &    action='write')
         write(atomf,'(A,A,A,I1,A,I5)') '#  Interface results for: ',
     &        trim(cwork),', IVXC=',IVXC,', NK=',NK         
         write(atomf,993) NL, '# NL'
         write(atomf,994) (QSCR(IL),IL=1,NL)
         WRITE(atomf,997) AZ(IA),WS(IA),WSAV(IA),'# AZ, WSR'         
         write(atomf,995) EF, '# Fermi Energy'         
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
      IF(IVAC.EQ.0) WRITE(atomf,1104) DBA         
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
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
C
C                     INPUT UNITS
      IR1=1
      IR2=2      
      IR4=4
      IR5=5
      IR11=11
      IR13=13
C                    OUTPUT UNITS
      IW6=6
      IW7=7
      IW8=8
      IW9=9
C
      OPEN(UNIT=IR1,FILE='inpge',FORM='FORMATTED')
      OPEN(UNIT=IR2,FILE='inpch',FORM='FORMATTED')
      OPEN(UNIT=IR4,FILE='inpcp',FORM='FORMATTED')
      OPEN(UNIT=IR5,FILE='input',FORM='FORMATTED')
      OPEN(UNIT=IR11,FILE='inpbu',FORM='FORMATTED')
      OPEN(UNIT=IR13,FILE='cci.nw',FORM='FORMATTED')
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNAB=20)
      PARAMETER(MNAV=1)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNR=400)
      PARAMETER(MNUH=48)
      integer,parameter::atomf=999      
C
      DIMENSION OWORK(5),OCASE(0:1),ORELA(2,0:2),OXCP(2,0:2),
     &          SCX(3)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/BPOPA/ BENY(MNL,MNSB,MNAB),BPPC(MNL,MNSB,MNAB),
     &              BPPD(MNL,MNSB,MNAB),BPPQ(MNL,MNSB,MNAB),
     &              BPPP(MNL,MNSB,MNAB),BDNY(MNL,MNSB,MNAB),
     &              BFINY(MNL,MNSB,MNAB),BFINYD(MNL,MNSB,MNAB)
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/DIBA/ DBA
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/TEXT/ OTXTA(MNA)
      DIMENSION OTXTA1(MNA)
      character(len=16)::cwork
      integer:: nrad1
      real(kind=8)::skippot
C
      DATA OCASE/'    SURFACE     ','   INTERFACE    '/
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
111   FORMAT(/'       ****  OUTPUT OF IES  **** '/
     &     /'     ** READING OF ALL DATA:')
C
C-------------------------  CONTROL DATA - UNIT IR5
C
1001  FORMAT(//5A16)
102   FORMAT(1X,F10.5,I5)
107   FORMAT(1X,2F10.5)
C
      READ (IR5,100) OWORK
      WRITE(IW6,1001) OWORK
      READ (IR5,100) OWORK
      WRITE(IW6,1001) OWORK
C
      READ (IR5,101) IVAC,NL,NS
      IF(IVAC.LT.0.OR.IVAC.GT.1) GO TO 209
      IF(NL.LT.3.OR.NL.GT.4) GO TO 210
      IF(NL.GT.MNL) GO TO 211
      IF(NS.GT.MNS) GO TO 212
      READ (IR5,101) IREL,IVXC
      IF(IREL.LT.0.OR.IREL.GT.2) GO TO 215
      IF(NS.EQ.2.AND.IREL.EQ.2) GO TO 2152
      IF(IVXC.LT.0.OR.IVXC.GT.2) GO TO 216
      READ (IR5,101) NSYM,NK,INVE
      READ (IR5,101) NE
      IF(NE.GT.MNE) GO TO 217
      IF(NE.NE.7.AND.NE.NE.10.AND.NE.NE.14.AND.NE.NE.20.
     &  AND.NE.NE.28.AND.NE.NE.40.AND.NE.NE.56.AND.NE.NE.80)
     &   GO TO 218
      READ (IR5,107) DIAM
      READ (IR5,101) ICONT
      IF(ICONT.LT.0.OR.ICONT.GT.1) GO TO 219
      READ (IR5,101) NFUPR
      READ (IR5,107) DIFMC,QMIXC
      READ (IR5,107) DIFMS,QMIXS
      READ (IR5,102) ALFA,NITER
      READ (IR5,102) BETA,NITERA
      NITERA=MAX(NITERA,2)
      READ (IR5,102) W0AM,NUH
      NUH=MAX(NUH,2)
      IF(NUH.GT.MNUH) GO TO 220
C
      WRITE(IW6,115) OCASE(IVAC),NL,NS
115   FORMAT(/4X,'CASE:',A16,7X,'NL=',I2,7X,'NS=',I2)
      WRITE(IW6,1161) (ORELA(I,IREL),I=1,2)
1161  FORMAT(/8X,'RELATIVITY:   ',2A16)
      WRITE(IW6,1162) (OXCP(I,IVXC),I=1,2)
1162  FORMAT(/8X,'XC-POTENTIAL:   ',2A16)
      WRITE(IW6,117) NSYM,NK,INVE
117   FORMAT(/'       LOCAL SYMMETRY:    NSYM=',I2/
     &      /'   K||-MESH:   NK=',I3,10X,'INVE=',I2)
      WRITE(IW6,118) NE,DIAM
118   FORMAT(/'  ENERGY CONTOUR:    NE=',I3,12X,'DIAM=',F12.5)
      WRITE(IW6,119) ICONT
119   FORMAT(/'    CONTINUATION:  ICONT=',I1)
      WRITE(IW6,120) NFUPR
120   FORMAT(/'          FULL PRINT AFTER',I4,' ITERATIONS')
      WRITE(IW6,121) DIFMC,QMIXC
121   FORMAT(/'  CPA ITERATIONS:  DIFMAX=',G12.4,
     &        8X,'MIXING=',F10.5)
      WRITE(IW6,122) DIFMS,QMIXS
122   FORMAT(/'  SGF ITERATIONS:  DIFMAX=',G12.4,
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
C
C-------------------------  GEOMETRY - UNIT IR1
C
      READ (IR1,100) OWORK
      WRITE(IW6,1001) OWORK
C
      READ (IR1,101) NP,NB,NMTR
      IF(NP.GT.MNP) GO TO 221
      IF(NB.GT.MNB) GO TO 222
      READ (IR1,104) CUTRAT
      READ (IR1,104) (SCX(I),I=1,3)
      READ (IR1,104) VBR(1,1),VBR(2,1)
      READ (IR1,104) VBR(1,2),VBR(2,2)
      READ (IR1,104) (VTRV(I),I=1,3)
      READ (IR1,104) (BTRV(I),I=1,3)
      DO 310 IB=1,NB
      READ (IR1,104) (VPOS(I,IB),I=1,3)
310   CONTINUE
      DO 311 IP=1,NP
      DO 312 IB=1,NB
      READ (IR1,104) (POS(I,IB,IP),I=1,3)
312   CONTINUE
311   CONTINUE
      DO 313 IB=1,NB
      READ (IR1,104) (BPOS(I,IB),I=1,3)
313   CONTINUE
C
C                          SCALING OF COORDINATES
      DO 314 J=1,2
      DO 315 I=1,2
      VBR(I,J)=SCX(I)*VBR(I,J)
315   CONTINUE
314   CONTINUE
      DO 316 I=1,3
      VTRV(I)=SCX(I)*VTRV(I)
316   CONTINUE
      DO 317 I=1,3
      BTRV(I)=SCX(I)*BTRV(I)
317   CONTINUE
      DO 318 IB=1,NB
      DO 319 I=1,3
      VPOS(I,IB)=SCX(I)*VPOS(I,IB)
319   CONTINUE
318   CONTINUE
      DO 320 IP=1,NP
      DO 321 IB=1,NB
      DO 322 I=1,3
      POS(I,IB,IP)=SCX(I)*POS(I,IB,IP)
322   CONTINUE
321   CONTINUE
320   CONTINUE
      DO 323 IB=1,NB
      DO 324 I=1,3
      BPOS(I,IB)=SCX(I)*BPOS(I,IB)
324   CONTINUE
323   CONTINUE
C
      WRITE(IW6,130) NP,NB,NMTR,CUTRAT
130   FORMAT(/8X,'NP=',I4,8X,'NB=',I2//'   MAX. ',
     & 'COEFFICIENT OF 2D-TRANSL. VECTORS:    NMTR=',I2/
     & /'   CUT-OFF DISTANCE / WS-RADIUS:  CUTRAT=',F10.5)
      WRITE(IW6,131) (SCX(I),I=1,3)
131   FORMAT(/'  SCALING FACTORS: ',3G15.7)
      WRITE(IW6,132)((VBR(I,J),I=1,2),J=1,2)
132   FORMAT(/'  TRANSL. VECTORS: 1. VECTOR = ',2G15.7/
     &        '                   2. VECTOR = ',2G15.7)
      WRITE(IW6,133) (VTRV(I),I=1,3)
133   FORMAT(/'    TRANSL. VECTOR OF SEMIINFINITE VACUUM:'/
     &      4X,3G15.7)
      WRITE(IW6,134) (BTRV(I),I=1,3)
134   FORMAT(/'    TRANSL. VECTOR OF SEMIINFINITE BULK:'/
     &      4X,3G15.7)
1041  FORMAT(3X,'IB=',I2,5X,3G15.7)
      WRITE(IW6,135)
135   FORMAT(/'       SITES OF VACUUM PRINCIPAL LAYER:')
      DO 3101 IB=1,NB
      WRITE(IW6,1041) IB,(VPOS(I,IB),I=1,3)
3101  CONTINUE
      DO 3111 IP=1,NP
      WRITE(IW6,136) IP
136   FORMAT(/'       SITES OF PRINCIPAL LAYER NO.: ',I4)
      DO 3121 IB=1,NB
      WRITE(IW6,1041) IB,(POS(I,IB,IP),I=1,3)
3121  CONTINUE
3111  CONTINUE
      WRITE(IW6,137)
137   FORMAT(/'       SITES OF BULK PRINCIPAL LAYER: ')
      DO 3131 IB=1,NB
      WRITE(IW6,1041) IB,(BPOS(I,IB),I=1,3)
3131  CONTINUE
C
C------------------------- CHEMICAL OCCUPATION - UNIT IR2
C
      READ (IR2,100) OWORK
      WRITE(IW6,1001) OWORK
C
108   FORMAT(A16)
1108   FORMAT(A16,A16)      
      IG=0
      IA=0
      DO 340 IP=1,NP
      DO 341 IB=1,NB
      IG=IG+1
      READ (IR2,101) NC(IG)
      DO 342 IC=1,NC(IG)
      IA=IA+1
      IF(IA.GT.MNA) GO TO 241
      READ (IR2,1108) OTXTA(IA),OTXTA1(IA)
      READ (IR2,104) CON(IA),VALZ(IA)
342   CONTINUE
341   CONTINUE
340   CONTINUE
      NA=IA
C
      IG=0
      IA=0
      DO 343 IP=1,NP
      DO 344 IB=1,NB
      IG=IG+1
      WRITE(IW6,141) IP,IB,IG,NC(IG)
141   FORMAT(/'  SITE:  IP=',I4,'   IB=',I2,'   IG=',I4,
     &     '   OCCUPIED BY ',I2,' ATOMIC TYPE(S): ')
      DO 345 IC=1,NC(IG)
      IA=IA+1
      WRITE(IW6,142) IA,OTXTA(IA),OTXTA1(IA),CON(IA),VALZ(IA)
142   FORMAT(' IA=',I4,' - LABEL,INIT,CONC.,VALENCY: ',
     &                       A16,A16,F10.5,F10.3)
345   CONTINUE
344   CONTINUE
343   CONTINUE
      WRITE(IW6,143) NA
143   FORMAT(/4X,' TOTAL NUMBER OF DIFFERENT ',
     &    'ATOMIC TYPES:  NA= ',I4)
C
C ---------------------------- LSDA FILE - UNIT IR3
C

C---------------------------------------------------      

      do IA=1,NA
         if (ICONT.EQ.0) then
            read(OTXTA1(IA),*) cwork
            WRITE(IW6,*)' Read defs for:', CWORK
         else
            read(OTXTA(IA),*) cwork
            WRITE(IW6,*)' Read old defs for:', CWORK
         end if            
         open (unit=atomf, file='atoms/'//trim(cwork),
     &        action='read')       

         do is=1,3
            read(atomf,*) cwork
         end do
         READ (atomf,104) AZ(IA),WS(IA),WSAV(IA)

         do is=1,4
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
         IF(IVAC.EQ.0) READ (atomf,104) DBA                  
         close(atomf)
      end do

         
C---------------------------------------------------      
C
      IF(IVAC.EQ.0) WRITE(IW6,149) DBA
149   FORMAT(/6X,' DIPOLE BARRIER=',F12.5)
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
     &          NSZRAD(IA),NUMCOR(IA)
356   CONTINUE
154   FORMAT(4X,A16,I5,2I8)


C ---------------  DATA FOR LEFT (VAC.) SUBSTRATE - UNIT IR11
C
      IF(IVAC.EQ.1) THEN
      READ (IR11,100) OWORK
      WRITE(IW6,1001) OWORK
C
      READ (IR11,101) NAV,NSV
      IF(NAV.GT.MNAV) GO TO 243
      IF(NSV.GT.MNSV) GO TO 2122
      IF(NS.LT.NSV) GO TO 214
      READ (IR11,104) (VCON(IA),IA=1,NAV)
      do IA=1,NAV
         READ (IR11,108) OTXTA1(IA)
         read(OTXTA1(IA),*) cwork
         WRITE(IW6,*)' Read vac. defs for:', CWORK
         open (unit=atomf, file='atoms/'//trim(cwork),
     &        action='read')
         read(atomf,*) cwork
         read(atomf,*) cwork
         read(atomf,*) cwork            
         READ (atomf,104) VAZ(IA),VWS(IA),VWSAV(IA)
         READ (atomf,104) EFV
         read(atomf,*) cwork
         read(atomf,*) cwork
         read(atomf,*) cwork          
         do IS=1,NSV
            READ (atomf,101) NRAD1
            READ (atomf,104) (skipPOT,I=1,NRAD1)
         end do         
         read(atomf,*) cwork
         do IS=1,NSV
            READ (atomf,104) cwork
            READ (atomf,104) (VENY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VPPC(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VPPD(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VPPQ(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VPPP(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VDNY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VFINY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (VFINYD(IL,IS,IA),IL=1,NL)
         end do
         close(atomf)
      end do
C
      WRITE(IW6,170) EFV,NAV,NSV
170   FORMAT(/4X,'VAC. SUBSTRATE:  EF=',G15.7,
     &        8X,'NAV=',I2,8X,'NSV=',I2)
      WRITE(IW6,161)
      DO 372 IA=1,NAV
      WRITE(IW6,162) IA,VCON(IA),VAZ(IA),VWS(IA),VWSAV(IA)
372   CONTINUE
      END IF
C
      
      
C ---------------  DATA FOR 1ST (BULK) SUBSTRATE - UNIT IR11
C
      READ (IR11,100) OWORK
      WRITE(IW6,1001) OWORK
C
      READ (IR11,101) NAB,NSB
      IF(NAB.GT.MNAB) GO TO 242
      IF(NSB.GT.MNSB) GO TO 2121
      IF(NS.LT.NSB) GO TO 213
      READ (IR11,104) (BCON(IA),IA=1,NAB)
      do IA=1,NAB
         READ (IR11,108) OTXTA1(IA)
         read(OTXTA1(IA),*) cwork
         WRITE(IW6,*)' Read bulk defs for:', CWORK
         open (unit=atomf, file='atoms/'//trim(cwork),
     &        action='read')
         read(atomf,*) cwork
         read(atomf,*) cwork
         read(atomf,*) cwork         
         READ (atomf,104) BAZ(IA),BWS(IA),BWSAV(IA)
         READ (atomf,104) EFB
         read(atomf,*) cwork
         read(atomf,*) cwork
         read(atomf,*) cwork         
         do IS=1,NSB
            READ (atomf,101) NRAD1
            READ (atomf,104) (skipPOT,I=1,NRAD1)
         end do         
         read(atomf,*) cwork
         do IS=1,NSB
            READ (atomf,104) cwork
            READ (atomf,104) (BENY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BPPC(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BPPD(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BPPQ(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BPPP(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BDNY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BFINY(IL,IS,IA),IL=1,NL)
            READ (atomf,104) (BFINYD(IL,IS,IA),IL=1,NL)
         end do
         close(atomf)
      end do
C
      WRITE(IW6,160) EFB,NAB,NSB
160   FORMAT(/4X,'BULK SUBSTRATE:  EF=',G15.7,
     &        8X,'NAB=',I2,8X,'NSB=',I2)
      WRITE(IW6,161)
161   FORMAT(/4X,'COMPONENT,   CONC.,      Z,',
     &           '       WS,        WSAV: ')
      DO 362 IA=1,NAB
      WRITE(IW6,162) IA,BCON(IA),BAZ(IA),BWS(IA),BWSAV(IA)
362   CONTINUE
162   FORMAT(10X,I2,F12.6,F9.3,2F12.6)
C
C ------------------------ COMPLEX NODES AND WEIGHTS - UNIT IR13
C
      READ (IR13,100) OWORK
      WRITE(IW6,1001) OWORK
C
      DO 380 J=1,8
      READ (IR13,101) NE1
      DO 381 IE=1,NE1
      READ (IR13,104) ZCN(IE),ZCW(IE)
381   CONTINUE
      IF(NE1.EQ.NE) GO TO 280
380   CONTINUE
280   CONTINUE
C
      WRITE(IW6,195)
195   FORMAT(/'     **  END OF READING ')
      RETURN
C
C-------------------------------- ERROR MESSAGES
C
209   WRITE(IW6,909)
909   FORMAT(/' **** INPUT ERROR:  IVAC MUST BE 0 OR 1 ')
      STOP
210   WRITE(IW6,910)
910   FORMAT(/' **** INPUT ERROR:  NL MUST BE 3 OR 4 ')
      STOP
211   WRITE(IW6,911)
911   FORMAT(/' **** INPUT ERROR:  NL GREATER THAN MNL ')
      STOP
212   WRITE(IW6,912)
912   FORMAT(/' **** INPUT ERROR:  NS GREATER THAN MNS ')
      STOP
2121  WRITE(IW6,9121)
9121  FORMAT(/' **** INPUT ERROR:  NSB GREATER THAN MNSB ')
      STOP
2122  WRITE(IW6,9122)
9122  FORMAT(/' **** INPUT ERROR:  NSV GREATER THAN MNSV ')
      STOP
213   WRITE(IW6,913)
913   FORMAT(/' **** INPUT ERROR:  NSB GREATER THAN NS ')
      STOP
214   WRITE(IW6,914)
914   FORMAT(/' **** INPUT ERROR:  NSV GREATER THAN NS ')
      STOP
215   WRITE(IW6,915)
915   FORMAT(/' **** INPUT ERROR:  IREL MUST BE 0, 1 OR 2 ')
      STOP
2152  WRITE(IW6,9152)
9152  FORMAT(/' **** INPUT ERROR:  NS=2 AND IREL=2 ')
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
221   WRITE(IW6,921)
921   FORMAT(/' **** INPUT ERROR:  NP GREATER THAN MNP  ')
      STOP
222   WRITE(IW6,922)
922   FORMAT(/' **** INPUT ERROR:  NB GREATER THAN MNB  ')
      STOP
241   WRITE(IW6,941)
941   FORMAT(/' **** INPUT ERROR:  IA EXCEEDED MNA  ')
      STOP
242   WRITE(IW6,942)
942   FORMAT(/' **** INPUT ERROR:  NAB GREATER THAN MNAB ')
      STOP
243   WRITE(IW6,943)
943   FORMAT(/' **** INPUT ERROR:  NAV GREATER THAN MNAV ')
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNAB=20)
      PARAMETER(MNSB=2)
      PARAMETER(MNAV=1)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNR=400)
C
      DIMENSION QSC3(3),QSC4(4)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/BPOPA/ BENY(MNL,MNSB,MNAB),BPPC(MNL,MNSB,MNAB),
     &              BPPD(MNL,MNSB,MNAB),BPPQ(MNL,MNSB,MNAB),
     &              BPPP(MNL,MNSB,MNAB),BDNY(MNL,MNSB,MNAB),
     &              BFINY(MNL,MNSB,MNAB),BFINYD(MNL,MNSB,MNAB)
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/DIBA/ DBA
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/TEXT/ OTXTA(MNA)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/
      DATA RCH/0.5D0/,RC4/4.0D0/,TOL/1.0D-6/,EPSCOR/1.0D-9/
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
      WRITE(IW6,114) (QSCR(IL),IL=1,NL)
114   FORMAT(/3X,'SCREENING CONSTANTS: '/
     &     5X,4G15.7)
C
C---------------------------------- PRECISION OF CORES
      THRESH=EPSCOR
      WRITE(IW6,115) THRESH
115   FORMAT(/9X,' PRECISION OF CORE ENERGIES: ',G15.7)
C
C------------------------- ACCURACY OF RUNGE-KUTTA METHOD
C                            (1.LE.NSIRK.LE.5)
      NSIRK=2
      WRITE(IW6,116) NSIRK
116   FORMAT(/9X,' ACCURACY OF RUNGE-KUTTA METHOD: ',
     &      '  NSIRK=',I2)
C
C---------------------------------- MAX. NO. OF SGF ITERATIONS
      NMITS=500
      WRITE(IW6,117) NMITS
117   FORMAT(/3X,' MAX. NUMBER OF SGF ITERATIONS: ',I5)
C
C---------------------------------- SHIFT OF ENY'S
      DO 305 IA=1,NA
      DO 306 IS=1,NS
      DO 307 IL=1,NL
      ISHENY(IL,IS,IA)=1
307   CONTINUE
306   CONTINUE
305   CONTINUE
      WRITE(IW6,118)
118   FORMAT(/2X,'  LABEL,             IA,  IS,',
     &       '    ISHENY - FOR SHIFT OF ENY:')
      DO 308 IA=1,NA
      DO 309 IS=1,NS
      WRITE(IW6,119) OTXTA(IA),IA,IS,
     &             (ISHENY(IL,IS,IA),IL=1,NL)
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
C---------------------------------- SYSTEM FERMI ENERGY AND
C                                   SCALING OF COMPLEX CONTOUR
C                                   NODES AND WEIGHTS
      IF(IVAC.EQ.0) EFT=EFB
      IF(IVAC.EQ.1) EFT=RCH*(EFB+EFV)
      DO 310 IE=1,NE
      ZCN(IE)=DIAM*ZCN(IE)+EFT
      ZCW(IE)=DIAM*ZCW(IE)
310   CONTINUE
      EFW=EFT
      WRITE(IW6,120) EFT
120   FORMAT(/5X,'****  SYSTEM FERMI ENERGY= ',G15.7)
      WRITE(IW6,121) ZCN(1),ZCN(NE)
121   FORMAT(/5X,'ENERGY CONTOUR: FIRST NODE= ',2G15.7/
     &        5X,'                 LAST NODE= ',2G15.7)
C
C---------------------------- CHECK OF CONSISTENCY
C                             OF CONCENTRATIONS
      NG=NB*NP
      IA=0
      DO 312 IG=1,NG
      SUM=RCZ
      DO 313 IC=1,NC(IG)
      IA=IA+1
      IF(CON(IA).LT.RCZ.OR.CON(IA).GT.RC1) GO TO 280
      SUM=SUM+CON(IA)
313   CONTINUE
      IF(ABS(SUM-RC1).GT.TOL) GO TO 280
312   CONTINUE
C
      SUM=RCZ
      DO 3141 IA=1,NAB
      IF(BCON(IA).LT.RCZ.OR.BCON(IA).GT.RC1) GO TO 281
      SUM=SUM+BCON(IA)
3141  CONTINUE
      IF(ABS(SUM-RC1).GT.TOL) GO TO 281
C
      IF(IVAC.EQ.1) THEN
      SUM=RCZ
      DO 3142 IA=1,NAV
      IF(VCON(IA).LT.RCZ.OR.VCON(IA).GT.RC1) GO TO 282
      SUM=SUM+VCON(IA)
3142  CONTINUE
      IF(ABS(SUM-RC1).GT.TOL) GO TO 282
      END IF
C
C----------------------------------- BASIS IN RECIP. SPACE
C
      DET=VBR(1,1)*VBR(2,2)-VBR(1,2)*VBR(2,1)
      ROMEGA=ABS(DET)
      GOMEGA=(RC2*PI)**2/ROMEGA
      CNORM=RC2*PI/DET
C
      VBG(1,1)= CNORM*VBR(2,2)
      VBG(2,1)=-CNORM*VBR(1,2)
      VBG(1,2)=-CNORM*VBR(2,1)
      VBG(2,2)= CNORM*VBR(1,1)
C
      WRITE(IW6,124) ROMEGA,GOMEGA
124   FORMAT(/' AREA OF 2D-PRIMITIVE CELL: REAL= ',G15.7/
     &        '                      RECIPROCAL= ',G15.7)
      WRITE(IW6,125)((VBG(I,J),I=1,2),J=1,2)
125   FORMAT(/' RECIP. BASIS: 1. VECTOR = ',2G15.7/
     &        '               2. VECTOR = ',2G15.7)
C
C--------------------------- CHECK OF CONSISTENCY OF DIRECTIONS
C                            OF BULK AND VACUUM TRANSL. VECTORS
      PROD=BTRV(3)*VTRV(3)
      IF (PROD.GE.RCZ) GO TO 291
C
C----------------------------------- DIMENSIONLESS WS-RADII
C                                        OF BULK AND VACUUM
      AVOLB=ROMEGA*ABS(BTRV(3))/REAL(NB)
      ABWS=((RC3*AVOLB)/(RC4*PI))**(RC1/RC3)
      AVOLV=ROMEGA*ABS(VTRV(3))/REAL(NB)
      AVWS=((RC3*AVOLV)/(RC4*PI))**(RC1/RC3)
C
      WRITE(IW6,126) ABWS,AVWS
126   FORMAT(/' DIMENSIONLESS WS-RADII:  BULK= ',G15.7/
     &        '                        VACUUM= ',G15.7)
C
C--------------------------------- TRUE AVERAGE WS-RADII
C                                     OF BULK AND VACUUM
      SUM=RCZ
      DO 321 IA=1,NAB
      SUM=SUM+BCON(IA)*BWSAV(IA)**3
321   CONTINUE
      BWST=SUM**(RC1/RC3)
C
      IF(IVAC.EQ.0) THEN
      VWST=BWST*AVWS/ABWS
      END IF
C
      IF(IVAC.EQ.1) THEN
      SUM=RCZ
      DO 322 IA=1,NAV
      SUM=SUM+VCON(IA)*VWSAV(IA)**3
322   CONTINUE
      VWST=SUM**(RC1/RC3)
      END IF
C
      WRITE(IW6,127) BWST,VWST
127   FORMAT(/' TRUE AVERAGE WS-RADII:  BULK= ',G15.7/
     &        '                       VACUUM= ',G15.7)
C
      DIF=(BWST*AVWS)/(VWST*ABWS)-RC1

      IF(ABS(DIF).GT.TOL*10.0d0) GO TO 292
c$$$      IF(ABS(DIF).GT.TOL) GO TO 292
C
C------------------------- DIMENSIONLESS AVERAGE WS-RADII
C                          IN INTERMEDIATE REGION
      FACT=ABWS/BWST
      IA=0
      IG=0
      DO 324 IP=1,NP
      DO 325 IB=1,NB
      IG=IG+1
      SUM=RCZ
      DO 326 IC=1,NC(IG)
      IA=IA+1
      SUM=SUM+CON(IA)*WSAV(IA)**3
326   CONTINUE
      AWS(IB,IP)=FACT*SUM**(RC1/RC3)
325   CONTINUE
324   CONTINUE
C
      WRITE(IW6,128)
128   FORMAT(/5X,'DIMENSIONLESS AVERAGE WS-RADII:'/)
      DO 328 IP=1,NP
      DO 329 IB=1,NB
      WRITE(IW6,129) IP,IB,AWS(IB,IP)
329   CONTINUE
328   CONTINUE
129   FORMAT(8X,'IP=',I4,5X,'IB=',I2,5X,'AWS=',G15.7)
C
C------------------------ POTENTIAL PARAMETERS
C                         OF THE SUBSTRATES
C                                    CASE OF SURFACE
      IF(IVAC.EQ.0) THEN
      NAV=1
      NSV=1
      VCON(1)=RC1
      VAZ(1)=RCZ
      VWS(1)=VWST
      VWSAV(1)=VWST
C
      S=VWST
      S2=S**2
      S4=S2**2
      SQS=SQRT(S)
      SQS3=S*SQS
C
      DO 331 IL=1,NL
      L=IL-1
      TL1=REAL(2*L+1)
      TL3=REAL(2*L+3)
      TL5=REAL(2*L+5)
      TL7=REAL(2*L+7)
      SQTL3=SQRT(TL3)
C
      VENY(IL,1,1)=DBA
      VPPC(IL,1,1)=DBA+TL1*TL5/(RC2*S2)
      VPPD(IL,1,1)=TL5**2/(RC2*RC4*TL3*S2)
      VPPQ(IL,1,1)=TL5/(RC4*TL1*TL3)
      VPPP(IL,1,1)=S4/(TL3*TL5**2*TL7)
      VDNY(IL,1,1)=REAL(L)
      VFINY(IL,1,1)=SQTL3/SQS3
      VFINYD(IL,1,1)=-SQS/(TL5*SQTL3)
331   CONTINUE
      END IF
C                                     CASE OF INTERFACE
      IF(IVAC.EQ.1) THEN
      DBA=EFB-EFV
      HDBA=RCH*DBA
      DO 333 IA=1,NAV
      DO 334 IS=1,NSV
      DO 335 IL=1,NL
      VENY(IL,IS,IA)=VENY(IL,IS,IA)+HDBA
      VPPC(IL,IS,IA)=VPPC(IL,IS,IA)+HDBA
335   CONTINUE
334   CONTINUE
333   CONTINUE
      DO 337 IA=1,NAB
      DO 338 IS=1,NSB
      DO 339 IL=1,NL
      BENY(IL,IS,IA)=BENY(IL,IS,IA)-HDBA
      BPPC(IL,IS,IA)=BPPC(IL,IS,IA)-HDBA
339   CONTINUE
338   CONTINUE
337   CONTINUE
      END IF
C
C---------------------------------- CHECK OF CONSISTENCY
C                                    OF NSYM, INVE AND
C                                 2D-TRANSLATION VECTORS
C
      IF(NSYM.LT.0.OR.NSYM.GT.4) GO TO 295
      IF(NSYM.EQ.0) GO TO 222
      IF(NSYM.EQ.1) GO TO 201
      IF(NSYM.EQ.2) GO TO 202
      IF(NSYM.EQ.3) GO TO 203
      IF(NSYM.EQ.4) GO TO 204
C
C--------  NSYM = 1 - SYMMETRY OF CUBIC (001) SURFACES:
C                     SYMMETRY GROUP GIVEN BY 4-FOLD ROTATION
C                     AXIS Z AND BY MIRROR PLANE X-Z
C
201   IF(INVE.LT.1.OR.INVE.GT.2) GO TO 296
      IF(INVE.EQ.1) GO TO 2011
      IF(INVE.EQ.2) GO TO 2012
C
C                     INVE = 1: THE 2D-TRANSL. VECTORS
C                               MUST HAVE THE FORM:
C                                      A1 = (P,0)
C                                      A2 = (0,P)
C                               WHERE P IS POSITIVE.
C
2011  IF (ABS(VBR(2,1)).GT.RCZ) GO TO 297
      IF (ABS(VBR(1,2)).GT.RCZ) GO TO 297
      IF (VBR(1,1).LE.RCZ) GO TO 297
      IF (VBR(2,2).LE.RCZ) GO TO 297
      DIF=VBR(1,1)/VBR(2,2)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      GO TO 222
C
C                     INVE = 2: THE 2D-TRANSL. VECTORS
C                               MUST HAVE THE FORM:
C                                      A1 = (P,P)
C                                      A2 = (-P,P)*U
C                               WHERE P IS POSITIVE
C                               AND U=1 OR U=-1.
C
2012  IF (VBR(1,1).LE.RCZ) GO TO 297
      DIF=VBR(2,1)/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      IF (ABS(VBR(1,2)).LE.RCZ) GO TO 297
      DIF=VBR(2,2)/VBR(1,2)+RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=ABS(VBR(2,2))/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      GO TO 222
C
C--------  NSYM = 2 - SYMMETRY OF CUBIC (110) SURFACES:
C                     SYMMETRY GROUP GIVEN BY TWO
C                     MIRROR PLANES X-Z AND Y-Z
C
C--------  NSYM = 4 - SYMMETRY OF CUBIC (N10) SURFACES
C                     AND TILT GRAIN BOUNDARIES:
C                     SYMMETRY GROUP GIVEN BY MIRROR PLANE X-Z
C
204   CONTINUE
202   IF(INVE.LT.1.OR.INVE.GT.2) GO TO 296
      IF(INVE.EQ.1) GO TO 2021
      IF(INVE.EQ.2) GO TO 2022
C
C                     INVE = 1: THE 2D-TRANSL. VECTORS
C                               MUST HAVE THE FORM:
C                                      A1 = (P,0)
C                                      A2 = (0,Q)
C                               WHERE P AND Q ARE POSITIVE.
C
2021  IF (ABS(VBR(2,1)).GT.RCZ) GO TO 297
      IF (ABS(VBR(1,2)).GT.RCZ) GO TO 297
      IF (VBR(1,1).LE.RCZ) GO TO 297
      IF (VBR(2,2).LE.RCZ) GO TO 297
      GO TO 222
C
C                     INVE = 2: THE 2D-TRANSL. VECTORS
C                               MUST HAVE THE FORM:
C                                      A1 = (P,Q)
C                                      A2 = (-P,Q)*U
C                               WHERE P AND Q ARE POSITIVE
C                               AND U=1 OR U=-1.
C
2022  IF (VBR(1,1).LE.RCZ) GO TO 297
      IF (VBR(2,1).LE.RCZ) GO TO 297
      DIF=ABS(VBR(1,2))/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=ABS(VBR(2,2))/VBR(2,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=(VBR(1,2)*VBR(2,2))/(VBR(1,1)*VBR(2,1))+RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      GO TO 222
C
C--------  NSYM = 3 - SYMMETRY OF CUBIC (111)  AND
C                     HEXAGONAL (0001) SURFACES:
C                     SYMMETRY GROUP GIVEN BY 3-FOLD ROTATION
C                     AXIS Z AND BY MIRROR PLANE Y-Z
C
C          INVE - NOT ACTIVE.   THE 2D-TRANSL. VECTORS
C                               MUST HAVE THE FORM:
C                                      A1 = (P,0)
C                                   A2 = (P/2)*(1,S)*U
C                               WHERE P IS POSITIVE,
C                               S=SQRT(3) OR S=-SQRT(3),
C                               AND U=1 OR U=-1.
C
203   IF (ABS(VBR(2,1)).GT.RCZ) GO TO 297
      IF (VBR(1,1).LE.RCZ) GO TO 297
      DIF=RC2*ABS(VBR(1,2))/VBR(1,1)-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
      DIF=(VBR(1,2)**2+VBR(2,2)**2)/VBR(1,1)**2-RC1
      IF (ABS(DIF).GT.TOL) GO TO 297
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
281   WRITE(IW6,181)
181   FORMAT(/' **** INPUT ERROR: '/3X,
     &'INCONSISTENCY IN CONCENTRATIONS OF BULK SUBSTRATE')
      STOP
282   WRITE(IW6,182)
182   FORMAT(/' **** INPUT ERROR: '/3X,
     &'INCONSISTENCY IN CONCENTRATIONS OF VAC. SUBSTRATE')
      STOP
291   WRITE(IW6,191)
191   FORMAT(/' **** INPUT ERROR: '/3X,
     &'VACUUM AND BULK TRANSLATION VECTORS ARE INCONSISTENT')
      STOP
292   WRITE(IW6,192)
192   FORMAT(/' **** INPUT ERROR: '/
     & '   VACUUM AND BULK AVERAGE DIMENSIONLESS AND TRUE '/
     & '   WS-RADII ARE INCONSISTENT')
      STOP
295   WRITE(IW6,195)
195   FORMAT(/' **** INPUT ERROR: NSYM MUST BE 0,1,2,3,4')
      STOP
296   WRITE(IW6,196)
196   FORMAT(/' **** INPUT ERROR:  NSYM AND INVE ',
     &  'ARE INCONSISTENT')
      STOP
297   WRITE(IW6,197)
197   FORMAT(/' **** INPUT ERROR:  NSYM, INVE '/3X,
     &  '  AND 2D-TRANSLATION VECTORS ARE INCONSISTENT')
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
      PARAMETER(MNCL=45)
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
C           WSA(ICL) (ICL=1,2,...,NCL) - LOCAL AVERAGE WS-RADII.
C  OUTPUT:  SSC(I,J,ICL) (I,J - ORBITAL INDICES,
C                        ICL=1,2,...,NCL) -
C           - TB-STRUCTURE CONSTANTS S(I,ICL,J,1) CONNECTING
C           THE ICL-TH ATOM WITH THE FIRST (CENTRAL) ATOM
C--------------------------------------------------------
C   REMARK: THE GAUNT FACTORS MUST BE CALCULATED BEFORE !
C           NCALLH MUST BE SET ZERO BEFORE !
C---------------------------------------------------------
C
      DIMENSION POL(3,MNCL),WSA(MNCL),RI(3),RJ(3),QR(MNLSQ),
     &          AA(MTBCL,MTBCL),BB(MTBCL,MNLSQ),WW(MTBCL),
     &          INDW(MNLSQ,MNCL),SC(MNLSQ,MNLSQ)
C
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
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
      CALL CANSC(NL,WSA(ICL),WSA(JCL),RI,RJ,SC)
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
C       Write (*,*) 'oies-pre--'
C       Write (*, '(9(1xg22.14))') (aa(i, 1:nlsq), i=1, nn)
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
      DO 320 J=1,NLSQ
      DO 321 ICL=1,NCL
      DO 322 I=1,NLSQ
      II=INDW(I,ICL)
      SSC(I,J,ICL)=QR(I)*BB(II,J)
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
      RETURN
      END
C*******************
CXXX   CANSC    ****
C*******************
      SUBROUTINE CANSC(NL,W1,W2,R1,R2,SC)
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
C       W1,W2 -  LOCAL AVERAGE W.S.-RADII
C       R1(I),R2(I) - THE COORDINATES OF THE 2 POINTS
C       SC(I1,I2) - THE STRUCTURE CONSTANT MATRIX
C---------------------------------------------------------------
C   REMARK: THE GAUNT FACTORS MUST BE CALCULATED BEFORE !
C           NCALLH MUST BE SET ZERO BEFORE !
C---------------------------------------------------------------
C
      DIMENSION R1(3),R2(3),SC(MNLSQ,MNLSQ)
      DIMENSION DFAC(0:MLMAX),AM1LP1(0:MNL)
      DIMENSION POW1(0:MNL),POW2(0:MNL),PREF(0:MNL,0:MNL)
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
      RBW=R/(W1+W2)
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
      RAT1=W1/R
      POW1(0)= SQRT(RAT1)
      DO 3131 L=1,LMAX
      POW1(L)=RAT1*POW1(L-1)
3131  CONTINUE
C
      RAT2=W2/R
      POW2(0)= SQRT(RAT2)
      DO 3132 L=1,LMAX
      POW2(L)=RAT2*POW2(L-1)
3132  CONTINUE
C
      DO 315 L2=0,LMAX
      DO 316 L1=0,LMAX
      L=L1+L2
      CIT=AM1LP1(L2)*PI8*DFAC(L)*POW1(L1)*POW2(L2)
      AJM=DFAC(L1)*DFAC(L2)
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
      U(L,L) =DBLE(ZQ)
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
CXXX    TBRI    ****
C*******************
      SUBROUTINE TBRI
C
C----------------------------------------------------
C   TB-STRUCTURE CONSTANTS FOR INTERMEDIATE REGION
C----------------------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNPP2=MNP+2)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=45)
C
      DIMENSION POW(3,MNB,-1:MNPP2),POL(3,MNCL),
     &          WSW(MNB,-1:MNPP2),WSL(MNCL),
     &          POLT(3),DPOL(3),
     &          WWMA(MNLSQ,MNLSQ),WWMB(MNLSQ,MNLSQ)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB,0:MNPP1),
     &            JDRP(MNCL,MNB,0:MNPP1),JBVA(MNCL,MNB,0:MNPP1),
     &            JTRA(2,MNCL,MNB,0:MNPP1),NSCL(MNB,0:MNPP1)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
111   FORMAT(//' *********  TB-CONSTANTS FOR',
     &    ' THE INTERMEDIATE REGION  *********')
C
      NLSQ=NL**2
C
C------------------------------- POSITIONS AND LOCAL AVERAGE
C                                WS-RADII INCLUDING TWO
C                                NEIGHBORING LAYERS ON BOTH
C                                SIDES OF INERMEDIATE REGION
      NPP1=NP+1
      NPP2=NP+2
      DO 301 IB=1,NB
      DO 3011 I=1,3
      POW(I,IB,-1)=VPOS(I,IB)+VTRV(I)
3011  CONTINUE
      WSW(IB,-1)=AVWS
301   CONTINUE
      DO 302 IB=1,NB
      DO 3021 I=1,3
      POW(I,IB,0)=VPOS(I,IB)
3021  CONTINUE
      WSW(IB,0)=AVWS
302   CONTINUE
      DO 3039 IP=1,NP
      DO 303 IB=1,NB
      DO 3031 I=1,3
      POW(I,IB,IP)=POS(I,IB,IP)
3031  CONTINUE
      WSW(IB,IP)=AWS(IB,IP)
303   CONTINUE
3039  CONTINUE
      DO 304 IB=1,NB
      DO 3041 I=1,3
      POW(I,IB,NPP1)=BPOS(I,IB)
3041  CONTINUE
      WSW(IB,NPP1)=ABWS
304   CONTINUE
      DO 305 IB=1,NB
      DO 3051 I=1,3
      POW(I,IB,NPP2)=BPOS(I,IB)+BTRV(I)
3051  CONTINUE
      WSW(IB,NPP2)=ABWS
305   CONTINUE
C
      DMINSQ=SMALL**2*ABWS*AVWS
      CUTRSQ=CUTRAT**2
C
C------------------------------- LOOP OVER CENTRAL ATOMS
C                                FOR SMALL CLUSTERS
      DO 310 IP=0,NPP1
      DO 311 IB=1,NB
C
C------------------------- SELECTION OF THE SMALL CLUSTER
      JCL=1
      DO 312 I=1,3
      POL(I,1)=POW(I,IB,IP)
312   CONTINUE
      JDRP(1,IB,IP)=0
      JBVA(1,IB,IP)=IB
      JTRA(1,1,IB,IP)=0
      JTRA(2,1,IB,IP)=0
      WSL(1)=WSW(IB,IP)
C
      DO 3139 JDIR=-1,1
      JP=IP+JDIR
      DO 313 JB=1,NB
      DO 3131 JTRA1=-NMTR,NMTR
      DO 3132 JTRA2=-NMTR,NMTR
C
      DO 314 I=1,2
      POLT(I)=POW(I,JB,JP)+JTRA1*VBR(I,1)+JTRA2*VBR(I,2)
314   CONTINUE
      POLT(3)=POW(3,JB,JP)
      DO 315 I=1,3
      DPOL(I)=POLT(I)-POL(I,1)
315   CONTINUE
C
      DMAXSQ=CUTRSQ*WSW(IB,IP)*WSW(JB,JP)
c$$$      DMAXSQ=CUTRSQ*AVWS*ABWS
      SUM=DPOL(1)**2+DPOL(2)**2+DPOL(3)**2
      IF (SUM.LT.DMINSQ.OR.SUM.GT.DMAXSQ) GO TO 3132
C
      JCL=JCL+1
      IF (JCL.GT.MNCL) GO TO 291
      DO 316 I=1,3
      POL(I,JCL)=POLT(I)
316   CONTINUE
      JDRP(JCL,IB,IP)=JDIR
      JBVA(JCL,IB,IP)=JB
      JTRA(1,JCL,IB,IP)=JTRA1
      JTRA(2,JCL,IB,IP)=JTRA2
      WSL(JCL)=WSW(JB,JP)
C
3132  CONTINUE
3131  CONTINUE
313   CONTINUE
3139  CONTINUE
C
      NCL=JCL
      NSCL(IB,IP)=NCL
      WRITE(IW6,115) IP,IB,NCL
115   FORMAT(/'    SITE:  IP=',I4,'  IB=',I3,
     &       '      CLUSTER SIZE= ',I3)
C
C--------------------------------- TB-CONSTANTS FOR THE
C                                  SMALL CLUSTER AND
C                                  THEIR STORAGE
      CALL TBCL(NL,NCL,POL,WSL)
C
      DO 320 JCL=1,NCL
      DO 321 JQ=1,NLSQ
      DO 322 IQ=1,NLSQ
      STR(IQ,JQ,JCL,IB,IP)=SSC(IQ,JQ,JCL)
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
311   CONTINUE
310   CONTINUE
C
C----------------------------------------- SYMMETRIZATION
C
C                       -------- IN-LAYER CONSTANTS
      DO 330 IP=1,NP
      JP=IP
      DEV=RCZ
C
      DO 331 IB=1,NB
      DO 332 JCL=1,NSCL(IB,IP)
      IF (JDRP(JCL,IB,IP).NE.0) GO TO 332
      JB=JBVA(JCL,IB,IP)
      DO 333 ICL=1,NSCL(JB,JP)
      IF (JDRP(ICL,JB,JP).NE.0) GO TO 333
      IF (JBVA(ICL,JB,JP).NE.IB) GO TO 333
      ISUM=JTRA(1,JCL,IB,IP)+JTRA(1,ICL,JB,JP)
      IF (ISUM.NE.0) GO TO 333
      ISUM=JTRA(2,JCL,IB,IP)+JTRA(2,ICL,JB,JP)
      IF (ISUM.NE.0) GO TO 333
C
      DO 3358 JQ=1,NLSQ
      DO 3359 IQ=1,NLSQ
      WWMA(IQ,JQ)=STR(IQ,JQ,JCL,IB,IP)
3359  CONTINUE
3358  CONTINUE
      DO 3368 JQ=1,NLSQ
      DO 3369 IQ=1,NLSQ
      WWMB(IQ,JQ)=STR(IQ,JQ,ICL,JB,JP)
3369  CONTINUE
3368  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3378 JQ=1,NLSQ
      DO 3379 IQ=1,NLSQ
      STR(IQ,JQ,JCL,IB,IP)=WWMA(IQ,JQ)
3379  CONTINUE
3378  CONTINUE
      DO 3388 JQ=1,NLSQ
      DO 3389 IQ=1,NLSQ
      STR(IQ,JQ,ICL,JB,JP)=WWMB(IQ,JQ)
3389  CONTINUE
3388  CONTINUE
C
333   CONTINUE
332   CONTINUE
331   CONTINUE
C
      WRITE(IW6,125) IP,DEV
125   FORMAT(/4X,'  LAYER:  IP=',I4,6X,'DEVIATION=',G12.4)
C
330   CONTINUE
C
C                       -------- OFF-LAYER CONSTANTS
      DO 340 IP=0,NP
      JP=IP+1
      DEV=RCZ
C
      DO 341 IB=1,NB
      DO 342 JCL=1,NSCL(IB,IP)
      IF (JDRP(JCL,IB,IP).NE.1) GO TO 342
      JB=JBVA(JCL,IB,IP)
      DO 343 ICL=1,NSCL(JB,JP)
      IF (JDRP(ICL,JB,JP).NE.-1) GO TO 343
      IF (JBVA(ICL,JB,JP).NE.IB) GO TO 343
      ISUM=JTRA(1,JCL,IB,IP)+JTRA(1,ICL,JB,JP)
      IF (ISUM.NE.0) GO TO 343
      ISUM=JTRA(2,JCL,IB,IP)+JTRA(2,ICL,JB,JP)
      IF (ISUM.NE.0) GO TO 343
C
      DO 3458 JQ=1,NLSQ
      DO 3459 IQ=1,NLSQ
      WWMA(IQ,JQ)=STR(IQ,JQ,JCL,IB,IP)
3459  CONTINUE
3458  CONTINUE
      DO 3468 JQ=1,NLSQ
      DO 3469 IQ=1,NLSQ
      WWMB(IQ,JQ)=STR(IQ,JQ,ICL,JB,JP)
3469  CONTINUE
3468  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3478 JQ=1,NLSQ
      DO 3479 IQ=1,NLSQ
      STR(IQ,JQ,JCL,IB,IP)=WWMA(IQ,JQ)
3479  CONTINUE
3478  CONTINUE
      DO 3488 JQ=1,NLSQ
      DO 3489 IQ=1,NLSQ
      STR(IQ,JQ,ICL,JB,JP)=WWMB(IQ,JQ)
3489  CONTINUE
3488  CONTINUE
C
343   CONTINUE
342   CONTINUE
341   CONTINUE
C
      WRITE(IW6,135) IP,JP,DEV
135   FORMAT(/4X,'LAYERS:  IP=',I4,'   JP=',I4,
     &        4X,'DEVIATION=',G12.4)
C
340   CONTINUE
C
C--------------------------------- PRINT OF DIAGONAL
C                                  ELEMENTS
      WRITE(IW6,150)
150   FORMAT(/5X,'***  DIAGONAL ELEMENTS OF ',
     &    'TB-STRUCTURE CONSTANTS :')
      DO 350 IP=0,NPP1
      DO 351 IB=1,NB
      WRITE(IW6,151) IP,IB
151   FORMAT(/'      SITE:  IP=',I4,'  IB=',I3)
      WRITE(IW6,104) (STR(IQ,IQ,1,IB,IP),IQ=1,NLSQ)
104   FORMAT(1X,4G15.7)
351   CONTINUE
350   CONTINUE
C
      RETURN
C
291   WRITE(IW6,191) JCL
191   FORMAT(/' **** ERROR IN TBRI :'/10X,
     &   ' JCL GREATER THAN MNCL,  JCL=',I3)
      STOP
      END
C*******************
CXXX    TBRB    ****
C*******************
      SUBROUTINE TBRB
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------------
C   TB-STRUCTURE CONSTANTS FOR SEMIINFINITE BULK
C----------------------------------------------------
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=45)
C
      DIMENSION POL(3,MNCL),WSL(MNCL),POLT(3),DPOL(3),
     &          WWMA(MNLSQ,MNLSQ),WWMB(MNLSQ,MNLSQ)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/BSCR/ BSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPB(MNCL,MNB),JBVAB(MNCL,MNB),
     &             JTRAB(2,MNCL,MNB),NSCLB(MNB)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
111   FORMAT(//' ************  TB-CONSTANTS FOR',
     &    ' THE BULK REGION  ************')
C
      NLSQ=NL**2
C
C------------------------------- CUT-OFF DISTANCE
      DMAXSQ=(CUTRAT*ABWS)**2
      DMINSQ=(SMALL*ABWS)**2
C
C------------------------------- LOOP OVER CENTRAL ATOMS
C                                FOR SMALL CLUSTERS
      DO 310 IB=1,NB
C
C------------------------- SELECTION OF THE SMALL CLUSTER
      JCL=1
      DO 312 I=1,3
      POL(I,1)=BPOS(I,IB)
312   CONTINUE
      JDRPB(1,IB)=0
      JBVAB(1,IB)=IB
      JTRAB(1,1,IB)=0
      JTRAB(2,1,IB)=0
      WSL(1)=ABWS
C
      DO 3139 JDIR=-1,1
      DO 313 JB=1,NB
      DO 3131 JTRA1=-NMTR,NMTR
      DO 3132 JTRA2=-NMTR,NMTR
C
      DO 314 I=1,2
      POLT(I)=BPOS(I,JB)+JDIR*BTRV(I)
     &       +JTRA1*VBR(I,1)+JTRA2*VBR(I,2)
314   CONTINUE
      POLT(3)=BPOS(3,JB)+JDIR*BTRV(3)
      DO 315 I=1,3
      DPOL(I)=POLT(I)-POL(I,1)
315   CONTINUE
C
      SUM=DPOL(1)**2+DPOL(2)**2+DPOL(3)**2
      IF (SUM.LT.DMINSQ.OR.SUM.GT.DMAXSQ) GO TO 3132
C
      JCL=JCL+1
      IF (JCL.GT.MNCL) GO TO 291
      DO 316 I=1,3
      POL(I,JCL)=POLT(I)
316   CONTINUE
      JDRPB(JCL,IB)=JDIR
      JBVAB(JCL,IB)=JB
      JTRAB(1,JCL,IB)=JTRA1
      JTRAB(2,JCL,IB)=JTRA2
      WSL(JCL)=ABWS
C
3132  CONTINUE
3131  CONTINUE
313   CONTINUE
3139  CONTINUE
C
      NCL=JCL
      NSCLB(IB)=NCL
      WRITE(IW6,115) IB,NCL
115   FORMAT(/'    SITE:    IB=',I3,
     &       '        CLUSTER SIZE= ',I3)
C
C--------------------------------- TB-CONSTANTS FOR THE
C                                  SMALL CLUSTER AND
C                                  THEIR STORAGE
      CALL TBCL(NL,NCL,POL,WSL)
C
      DO 320 JCL=1,NCL
      DO 321 JQ=1,NLSQ
      DO 322 IQ=1,NLSQ
      BSTR(IQ,JQ,JCL,IB)=SSC(IQ,JQ,JCL)
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
310   CONTINUE
C
C----------------------------------------- SYMMETRIZATION
C
C                       -------- IN-LAYER CONSTANTS
      DEV=RCZ
C
      DO 330 IB=1,NB
      DO 331 JCL=1,NSCLB(IB)
      IF (JDRPB(JCL,IB).NE.0) GO TO 331
      JB=JBVAB(JCL,IB)
      DO 332 ICL=1,NSCLB(JB)
      IF (JDRPB(ICL,JB).NE.0) GO TO 332
      IF (JBVAB(ICL,JB).NE.IB) GO TO 332
      ISUM=JTRAB(1,JCL,IB)+JTRAB(1,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
      ISUM=JTRAB(2,JCL,IB)+JTRAB(2,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
C
      DO 3358 JQ=1,NLSQ
      DO 3359 IQ=1,NLSQ
      WWMA(IQ,JQ)=BSTR(IQ,JQ,JCL,IB)
3359  CONTINUE
3358  CONTINUE
      DO 3368 JQ=1,NLSQ
      DO 3369 IQ=1,NLSQ
      WWMB(IQ,JQ)=BSTR(IQ,JQ,ICL,JB)
3369  CONTINUE
3368  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3378 JQ=1,NLSQ
      DO 3379 IQ=1,NLSQ
      BSTR(IQ,JQ,JCL,IB)=WWMA(IQ,JQ)
3379  CONTINUE
3378  CONTINUE
      DO 3388 JQ=1,NLSQ
      DO 3389 IQ=1,NLSQ
      BSTR(IQ,JQ,ICL,JB)=WWMB(IQ,JQ)
3389  CONTINUE
3388  CONTINUE
C
332   CONTINUE
331   CONTINUE
330   CONTINUE
C
      WRITE(IW6,125) DEV
125   FORMAT(/4X,' IN-LAYER CONSTANTS : DEVIATION=',G12.4)
C
C                       -------- OFF-LAYER CONSTANTS
      DEV=RCZ
C
      DO 340 IB=1,NB
      DO 341 JCL=1,NSCLB(IB)
      IF (JDRPB(JCL,IB).NE.1) GO TO 341
      JB=JBVAB(JCL,IB)
      DO 342 ICL=1,NSCLB(JB)
      IF (JDRPB(ICL,JB).NE.-1) GO TO 342
      IF (JBVAB(ICL,JB).NE.IB) GO TO 342
      ISUM=JTRAB(1,JCL,IB)+JTRAB(1,ICL,JB)
      IF (ISUM.NE.0) GO TO 342
      ISUM=JTRAB(2,JCL,IB)+JTRAB(2,ICL,JB)
      IF (ISUM.NE.0) GO TO 342
C
      DO 3458 JQ=1,NLSQ
      DO 3459 IQ=1,NLSQ
      WWMA(IQ,JQ)=BSTR(IQ,JQ,JCL,IB)
3459  CONTINUE
3458  CONTINUE
      DO 3468 JQ=1,NLSQ
      DO 3469 IQ=1,NLSQ
      WWMB(IQ,JQ)=BSTR(IQ,JQ,ICL,JB)
3469  CONTINUE
3468  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3478 JQ=1,NLSQ
      DO 3479 IQ=1,NLSQ
      BSTR(IQ,JQ,JCL,IB)=WWMA(IQ,JQ)
3479  CONTINUE
3478  CONTINUE
      DO 3488 JQ=1,NLSQ
      DO 3489 IQ=1,NLSQ
      BSTR(IQ,JQ,ICL,JB)=WWMB(IQ,JQ)
3489  CONTINUE
3488  CONTINUE
C
342   CONTINUE
341   CONTINUE
340   CONTINUE
C
      WRITE(IW6,135) DEV
135   FORMAT(/4X,'OFF-LAYER CONSTANTS : DEVIATION=',G12.4)
C
C--------------------------------- PRINT OF DIAGONAL
C                                  ELEMENTS
      WRITE(IW6,150)
150   FORMAT(/5X,'***  DIAGONAL ELEMENTS OF ',
     &    'TB-STRUCTURE CONSTANTS :')
      DO 350 IB=1,NB
      WRITE(IW6,151) IB
151   FORMAT(/'        SITE:       IB=',I3)
      WRITE(IW6,104) (BSTR(IQ,IQ,1,IB),IQ=1,NLSQ)
104   FORMAT(1X,4G15.7)
350   CONTINUE
C
      RETURN
C
291   WRITE(IW6,191) JCL
191   FORMAT(/' **** ERROR IN TBRB :'/10X,
     &   ' JCL GREATER THAN MNCL,  JCL=',I3)
      STOP
      END
C*******************
CXXX    TBRV    ****
C*******************
      SUBROUTINE TBRV
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C----------------------------------------------------
C   TB-STRUCTURE CONSTANTS FOR SEMIINFINITE VACUUM
C----------------------------------------------------
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=45)
C
      DIMENSION POL(3,MNCL),WSL(MNCL),POLT(3),DPOL(3),
     &          WWMA(MNLSQ,MNLSQ),WWMB(MNLSQ,MNLSQ)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/TBSC/ SSC(MNLSQ,MNLSQ,MNCL)
      COMMON/VSCR/ VSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPV(MNCL,MNB),JBVAV(MNCL,MNB),
     &             JTRAV(2,MNCL,MNB),NSCLV(MNB)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
111   FORMAT(//' ***********  TB-CONSTANTS FOR',
     &    ' THE VACUUM REGION  ***********')
C
      NLSQ=NL**2
C
C------------------------------- CUT-OFF DISTANCE
      DMAXSQ=(CUTRAT*AVWS)**2
      DMINSQ=(SMALL*AVWS)**2
C
C------------------------------- LOOP OVER CENTRAL ATOMS
C                                FOR SMALL CLUSTERS
      DO 310 IB=1,NB
C
C------------------------- SELECTION OF THE SMALL CLUSTER
      JCL=1
      DO 312 I=1,3
      POL(I,1)=VPOS(I,IB)
312   CONTINUE
      JDRPV(1,IB)=0
      JBVAV(1,IB)=IB
      JTRAV(1,1,IB)=0
      JTRAV(2,1,IB)=0
      WSL(1)=AVWS
C
      DO 3139 JDIR=-1,1
      DO 313 JB=1,NB
      DO 3131 JTRA1=-NMTR,NMTR
      DO 3132 JTRA2=-NMTR,NMTR
C
      DO 314 I=1,2
      POLT(I)=VPOS(I,JB)-JDIR*VTRV(I)
     &       +JTRA1*VBR(I,1)+JTRA2*VBR(I,2)
314   CONTINUE
      POLT(3)=VPOS(3,JB)-JDIR*VTRV(3)
      DO 315 I=1,3
      DPOL(I)=POLT(I)-POL(I,1)
315   CONTINUE
C
      SUM=DPOL(1)**2+DPOL(2)**2+DPOL(3)**2
      IF (SUM.LT.DMINSQ.OR.SUM.GT.DMAXSQ) GO TO 3132
C
      JCL=JCL+1
      IF (JCL.GT.MNCL) GO TO 291
      DO 316 I=1,3
      POL(I,JCL)=POLT(I)
316   CONTINUE
      JDRPV(JCL,IB)=JDIR
      JBVAV(JCL,IB)=JB
      JTRAV(1,JCL,IB)=JTRA1
      JTRAV(2,JCL,IB)=JTRA2
      WSL(JCL)=AVWS
C
3132  CONTINUE
3131  CONTINUE
313   CONTINUE
3139  CONTINUE
C
      NCL=JCL
      NSCLV(IB)=NCL
      WRITE(IW6,115) IB,NCL
115   FORMAT(/'    SITE:    IB=',I3,
     &       '        CLUSTER SIZE= ',I3)
C
C--------------------------------- TB-CONSTANTS FOR THE
C                                  SMALL CLUSTER AND
C                                  THEIR STORAGE
      CALL TBCL(NL,NCL,POL,WSL)
C
      DO 320 JCL=1,NCL
      DO 321 JQ=1,NLSQ
      DO 322 IQ=1,NLSQ
      VSTR(IQ,JQ,JCL,IB)=SSC(IQ,JQ,JCL)
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
310   CONTINUE
C
C----------------------------------------- SYMMETRIZATION
C
C                       -------- IN-LAYER CONSTANTS
      DEV=RCZ
C
      DO 330 IB=1,NB
      DO 331 JCL=1,NSCLV(IB)
      IF (JDRPV(JCL,IB).NE.0) GO TO 331
      JB=JBVAV(JCL,IB)
      DO 332 ICL=1,NSCLV(JB)
      IF (JDRPV(ICL,JB).NE.0) GO TO 332
      IF (JBVAV(ICL,JB).NE.IB) GO TO 332
      ISUM=JTRAV(1,JCL,IB)+JTRAV(1,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
      ISUM=JTRAV(2,JCL,IB)+JTRAV(2,ICL,JB)
      IF (ISUM.NE.0) GO TO 332
C
      DO 3358 JQ=1,NLSQ
      DO 3359 IQ=1,NLSQ
      WWMA(IQ,JQ)=VSTR(IQ,JQ,JCL,IB)
3359  CONTINUE
3358  CONTINUE
      DO 3368 JQ=1,NLSQ
      DO 3369 IQ=1,NLSQ
      WWMB(IQ,JQ)=VSTR(IQ,JQ,ICL,JB)
3369  CONTINUE
3368  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3378 JQ=1,NLSQ
      DO 3379 IQ=1,NLSQ
      VSTR(IQ,JQ,JCL,IB)=WWMA(IQ,JQ)
3379  CONTINUE
3378  CONTINUE
      DO 3388 JQ=1,NLSQ
      DO 3389 IQ=1,NLSQ
      VSTR(IQ,JQ,ICL,JB)=WWMB(IQ,JQ)
3389  CONTINUE
3388  CONTINUE
C
332   CONTINUE
331   CONTINUE
330   CONTINUE
C
      WRITE(IW6,125) DEV
125   FORMAT(/4X,' IN-LAYER CONSTANTS : DEVIATION=',G12.4)
C
C                       -------- OFF-LAYER CONSTANTS
      DEV=RCZ
C
      DO 340 IB=1,NB
      DO 341 JCL=1,NSCLV(IB)
      IF (JDRPV(JCL,IB).NE.1) GO TO 341
      JB=JBVAV(JCL,IB)
      DO 342 ICL=1,NSCLV(JB)
      IF (JDRPV(ICL,JB).NE.-1) GO TO 342
      IF (JBVAV(ICL,JB).NE.IB) GO TO 342
      ISUM=JTRAV(1,JCL,IB)+JTRAV(1,ICL,JB)
      IF (ISUM.NE.0) GO TO 342
      ISUM=JTRAV(2,JCL,IB)+JTRAV(2,ICL,JB)
      IF (ISUM.NE.0) GO TO 342
C
      DO 3458 JQ=1,NLSQ
      DO 3459 IQ=1,NLSQ
      WWMA(IQ,JQ)=VSTR(IQ,JQ,JCL,IB)
3459  CONTINUE
3458  CONTINUE
      DO 3468 JQ=1,NLSQ
      DO 3469 IQ=1,NLSQ
      WWMB(IQ,JQ)=VSTR(IQ,JQ,ICL,JB)
3469  CONTINUE
3468  CONTINUE
C
      CALL MASY(NLSQ,MNLSQ,WWMA,WWMB,DEV1)
      DEV= MAX(DEV,DEV1)
C
      DO 3478 JQ=1,NLSQ
      DO 3479 IQ=1,NLSQ
      VSTR(IQ,JQ,JCL,IB)=WWMA(IQ,JQ)
3479  CONTINUE
3478  CONTINUE
      DO 3488 JQ=1,NLSQ
      DO 3489 IQ=1,NLSQ
      VSTR(IQ,JQ,ICL,JB)=WWMB(IQ,JQ)
3489  CONTINUE
3488  CONTINUE
C
342   CONTINUE
341   CONTINUE
340   CONTINUE
C
      WRITE(IW6,135) DEV
135   FORMAT(/4X,'OFF-LAYER CONSTANTS : DEVIATION=',G12.4)
C
C--------------------------------- PRINT OF DIAGONAL
C                                  ELEMENTS
      WRITE(IW6,150)
150   FORMAT(/5X,'***  DIAGONAL ELEMENTS OF ',
     &    'TB-STRUCTURE CONSTANTS :')
      DO 350 IB=1,NB
      WRITE(IW6,151) IB
151   FORMAT(/'        SITE:       IB=',I3)
      WRITE(IW6,104) (VSTR(IQ,IQ,1,IB),IQ=1,NLSQ)
104   FORMAT(1X,4G15.7)
350   CONTINUE
C
      RETURN
C
291   WRITE(IW6,191) JCL
191   FORMAT(/' **** ERROR IN TBRV :'/10X,
     &   ' JCL GREATER THAN MNCL,  JCL=',I3)
      STOP
      END
C*******************
CXXX    GIBZ    ****
C*******************
      SUBROUTINE GIBZ
C
C****************************************************************
C  GENERATES NETWORK OF K||-POINTS IN IRREDUCIBLE BRILLOUIN ZONE
C****************************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNBZ=700)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
      DATA RCH/0.5D0/,RC3/3.0D0/,RC6/6.0D0/
      DATA DUM1/0.75D0/,DUM2/0.4D0/,DUM3/1.2D0/
C
      PI=RC4*ATAN(RC1)
C
      IF(NSYM.EQ.0) GO TO 200
      IF(NSYM.EQ.1) GO TO 201
      IF(NSYM.EQ.2) GO TO 202
      IF(NSYM.EQ.3) GO TO 203
      IF(NSYM.EQ.4) GO TO 204
C
C------------------------------- NSYM=0: GENERAL CASE
200   IBZ=0
      TW=RCZ
      NK1=NK
      NK2=NK
      DO 301 I1=-NK1+1,NK1
      FX1=(REAL(I1)-RCH)/REAL(2*NK1)
      DO 302 I2=1,NK2
      FX2=(REAL(I2)-RCH)/REAL(2*NK2)
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      DO 305 J=1,2
      AKBZ(J,IBZ)=FX1*VBG(J,1)+FX2*VBG(J,2)
305   CONTINUE
302   CONTINUE
301   CONTINUE
      GO TO 290
C
201   CONTINUE
      IF(INVE.EQ.1) GO TO 2011
      IF(INVE.EQ.2) GO TO 2012
C----------------------------------- NSYM=1, INVE=1
2011  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      DK=PI/(P*REAL(NK))
      DO 311 I=1,NK
      DO 312 J=1,I
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      IF(I.EQ.J) WKBZ(IBZ)=RCH
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I)-RCH)*DK
      AKBZ(2,IBZ)=(REAL(J)-RCH)*DK
312   CONTINUE
311   CONTINUE
      GO TO 290
C----------------------------------- NSYM=1, INVE=2
2012  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      DK=PI/(RC2*P*REAL(NK))
      DO 314 I=1,NK
      DO 315 J=1,I
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      IF(I.EQ.J) WKBZ(IBZ)=RCH
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=REAL(I-J)*DK
      AKBZ(2,IBZ)=REAL(I+J-1)*DK
315   CONTINUE
314   CONTINUE
      GO TO 290
C
202   CONTINUE
      IF(INVE.EQ.1) GO TO 2021
      IF(INVE.EQ.2) GO TO 2022
C----------------------------------- NSYM=2, INVE=1
2021  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      Q=VBR(2,2)
      NMX=1
      NMY=1
      IF(P.LT.(DUM1*Q)) NMX=2
      IF(Q.LT.(DUM1*P)) NMY=2
      IF(P.LT.(DUM2*Q)) WRITE(IW6,115)
      IF(Q.LT.(DUM2*P)) WRITE(IW6,115)
115   FORMAT(/' **** GIBZ - WARNING: ',
     &    '  P DEVIATES TOO MUCH FROM Q ')
      NDX=NK*NMX
      NDY=NK*NMY
      DKX=PI/(P*REAL(NDX))
      DKY=PI/(Q*REAL(NDY))
      DO 321 I=1,NDX
      DO 322 J=1,NDY
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I)-RCH)*DKX
      AKBZ(2,IBZ)=(REAL(J)-RCH)*DKY
322   CONTINUE
321   CONTINUE
      GO TO 290
C----------------------------------- NSYM=2, INVE=2
2022  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      Q=VBR(2,1)
      NMX=1
      NMY=1
      IF(P.LT.(DUM1*Q)) NMX=2
      IF(Q.LT.(DUM1*P)) NMY=2
      IF(P.LT.(DUM2*Q)) WRITE(IW6,125)
      IF(Q.LT.(DUM2*P)) WRITE(IW6,125)
125   FORMAT(/' **** GIBZ - WARNING: ',
     &    '  P DEVIATES TOO MUCH FROM Q ')
      NDX=NK*NMX
      NDY=NK*NMY
      DKX=PI/(P*REAL(NDX))
      DKY=PI/(Q*REAL(NDY))
      DO 325 I=1,NDX
      DO 326 J=1,NDY
      IF((2*J-1)*NMX.GT.(2*I-1)*NMY) GO TO 326
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      IF(I.EQ.J.AND.NMX.EQ.NMY) WKBZ(IBZ)=RCH
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I)-RCH)*DKX
      AKBZ(2,IBZ)=(REAL(J)-RCH)*DKY
326   CONTINUE
325   CONTINUE
      GO TO 290
C----------------------------------- NSYM=3
203   IBZ=0
      TW=RCZ
      P=VBR(1,1)
      DKX=RC2*PI/(P*REAL(NK))
      DKY=DKX/SQRT(RC3)
      DO 330 NU=1,2
      ANU3=REAL(NU)/RC3
      DO 331 I=1,NK
      DO 332 J=I,NK
      IY=I+2*J+NU-3
      IF(IY.GT.NK) GO TO 332
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC6
      IF(I.EQ.J.OR.IY.EQ.NK) WKBZ(IBZ)=RC3
      IF(I.EQ.J.AND.IY.EQ.NK) WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I-1)+ANU3)*DKX
      AKBZ(2,IBZ)=REAL(IY)*DKY
332   CONTINUE
331   CONTINUE
330   CONTINUE
      GO TO 290
C
204   CONTINUE
      IF(INVE.EQ.1) GO TO 2041
      IF(INVE.EQ.2) GO TO 2042
C----------------------------------- NSYM=4, INVE=1
2041  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      Q=VBR(2,2)
      RNK=REAL(NK)
      NDX=NK
      NDY=NK
      IF(P.GT.(DUM3*Q)) NDY=INT(P*RNK/Q)
      IF(Q.GT.(DUM3*P)) NDX=INT(Q*RNK/P)
      DKX=PI/(P*REAL(NDX))
      DKY=PI/(Q*REAL(NDY))
      DO 341 I=1,NDX
      DO 342 J=1,NDY
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I)-RCH)*DKX
      AKBZ(2,IBZ)=(REAL(J)-RCH)*DKY
342   CONTINUE
341   CONTINUE
      GO TO 290
C----------------------------------- NSYM=4, INVE=2
2042  IBZ=0
      TW=RCZ
      P=VBR(1,1)
      Q=VBR(2,1)
      RNK=REAL(NK)
      NDX=NK
      NDY=NK
      IF(P.GT.(DUM3*Q)) NDY=INT(P*RNK/Q)
      IF(Q.GT.(DUM3*P)) NDX=INT(Q*RNK/P)
      DKX=PI/(P*REAL(NDX))
      DKY=PI/(Q*REAL(NDY))
      DO 345 I=1,NDX
      DO 346 J=1,NDY
      IF((2*J-1)*NDX.GT.(2*I-1)*NDY) GO TO 346
      IBZ=IBZ+1
      IF(IBZ.GT.MNBZ) GO TO 299
      WKBZ(IBZ)=RC1
      IF((2*J-1)*NDX.EQ.(2*I-1)*NDY) WKBZ(IBZ)=RCH
      TW=TW+WKBZ(IBZ)
      AKBZ(1,IBZ)=(REAL(I)-RCH)*DKX
      AKBZ(2,IBZ)=(REAL(J)-RCH)*DKY
346   CONTINUE
345   CONTINUE
      GO TO 290
C
290   NBZ=IBZ
      DO 395 IBZ=1,NBZ
      WKBZ(IBZ)=WKBZ(IBZ)/TW
395   CONTINUE
      WRITE(IW6,190) NBZ
190   FORMAT(//3X,'*****  NUMBER OF K||-POINTS:   NBZ=',I6)
      RETURN
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
      PARAMETER(MNP=30)
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNBZ=700)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=45)
C
      DIMENSION TRAN(2)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB,0:MNPP1),
     &            JDRP(MNCL,MNB,0:MNPP1),JBVA(MNCL,MNB,0:MNPP1),
     &            JTRA(2,MNCL,MNB,0:MNPP1),NSCL(MNB,0:MNPP1)
      COMMON/BSCR/ BSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPB(MNCL,MNB),JBVAB(MNCL,MNB),
     &             JTRAB(2,MNCL,MNB),NSCLB(MNB)
      COMMON/VSCR/ VSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPV(MNCL,MNB),JBVAV(MNCL,MNB),
     &             JTRAV(2,MNCL,MNB),NSCLV(MNB)
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/BSCK/ ZBSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/VSCK/ ZVSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
C----------------------------------- INTERMEDIATE REGION
C
      DO 3101 IBZ=1,MNBZ
      DO 3102 IP=1,MNP
      DO 3103 II=1,MNBLSQ
      DO 3104 JJ=1,MNBLSQ
      ZSKI(JJ,II,IP,IBZ)=ZZ
      ZSKO(JJ,II,IP,IBZ)=ZZ
      ZSKF(JJ,II,IP,IBZ)=ZZ
3104  CONTINUE
3103  CONTINUE
3102  CONTINUE
3101  CONTINUE
C
      DO 311 IBZ=1,NBZ
      DO 312 IP=1,NP
      DO 3125 IB=1,NB
      II0=(IB-1)*NLSQ
C
      DO 313 JCL=1,NSCL(IB,IP)
      JDIR=JDRP(JCL,IB,IP)
      JB=JBVA(JCL,IB,IP)
      JTRA1=JTRA(1,JCL,IB,IP)
      JTRA2=JTRA(2,JCL,IB,IP)
C
      JJ0=(JB-1)*NLSQ
      DO 314 I=1,2
      TRAN(I)=JTRA1*VBR(I,1)+JTRA2*VBR(I,2)
314   CONTINUE
      ARG=RCZ
      DO 315 I=1,2
      ARG=ARG+AKBZ(I,IBZ)*TRAN(I)
315   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      IF(JDIR.EQ.0) THEN
      DO 3161 IQ=1,NLSQ
      DO 3162 JQ=1,NLSQ
      ZSKI(JJ0+JQ,II0+IQ,IP,IBZ) = ZSKI(JJ0+JQ,II0+IQ,IP,IBZ)
     &                           + ZEFA * STR(JQ,IQ,JCL,IB,IP)
3162  CONTINUE
3161  CONTINUE
      END IF
      IF(JDIR.EQ.1) THEN
      DO 3171 IQ=1,NLSQ
      DO 3172 JQ=1,NLSQ
      ZSKF(JJ0+JQ,II0+IQ,IP,IBZ) = ZSKF(JJ0+JQ,II0+IQ,IP,IBZ)
     &                           + ZEFA * STR(JQ,IQ,JCL,IB,IP)
3172  CONTINUE
3171  CONTINUE
      END IF
      IF(JDIR.EQ.-1) THEN
      DO 3181 IQ=1,NLSQ
      DO 3182 JQ=1,NLSQ
      ZSKO(JJ0+JQ,II0+IQ,IP,IBZ) = ZSKO(JJ0+JQ,II0+IQ,IP,IBZ)
     &                           + ZEFA * STR(JQ,IQ,JCL,IB,IP)
3182  CONTINUE
3181  CONTINUE
      END IF
C
313   CONTINUE
3125  CONTINUE
312   CONTINUE
311   CONTINUE
C                       CHECK OF HERMITICITY
      DO 320 IP=1,NP
      DEV=RCZ
      DO 321 IBZ=1,NBZ
      DO 322 II=1,NBLSQ
      DO 323 JJ=II,NBLSQ
      ZW1=ZSKI(II,JJ,IP,IBZ)
      W1X=DBLE(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKI(JJ,II,IP,IBZ)
      W2X=DBLE(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
323   CONTINUE
322   CONTINUE
321   CONTINUE
      WRITE(IW6,120) IP,DEV
120   FORMAT(/4X,'  LAYER:  IP=',I4,6X,'DEVIATION=',G12.4)
320   CONTINUE
C
      IF(NP.GT.1) THEN
      DO 325 IP=1,NP-1
      IPP1=IP+1
      DEV=RCZ
      DO 326 IBZ=1,NBZ
      DO 327 II=1,NBLSQ
      DO 328 JJ=1,NBLSQ
      ZW1=ZSKF(II,JJ,IP,IBZ)
      W1X=DBLE(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKO(JJ,II,IPP1,IBZ)
      W2X=DBLE(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
328   CONTINUE
327   CONTINUE
326   CONTINUE
      WRITE(IW6,125) IP,IPP1,DEV
125   FORMAT(/4X,'LAYERS:  IP=',I4,'  JP=',I4,
     &                       6X,'DEVIATION=',G12.4)
325   CONTINUE
      END IF
C
C----------------------------------- BULK REGION
C
      DO 3301 IBZ=1,MNBZ
      DO 3302 II=1,MNBLSQ
      DO 3303 JJ=1,MNBLSQ
      ZBSKI(JJ,II,IBZ)=ZZ
      ZBSKO(JJ,II,IBZ)=ZZ
      ZBSKF(JJ,II,IBZ)=ZZ
3303  CONTINUE
3302  CONTINUE
3301  CONTINUE
C
      DO 331 IBZ=1,NBZ
      DO 332 IB=1,NB
      II0=(IB-1)*NLSQ
C
      DO 333 JCL=1,NSCLB(IB)
      JDIR=JDRPB(JCL,IB)
      JB=JBVAB(JCL,IB)
      JTRA1=JTRAB(1,JCL,IB)
      JTRA2=JTRAB(2,JCL,IB)
C
      JJ0=(JB-1)*NLSQ
      DO 334 I=1,2
      TRAN(I)=JTRA1*VBR(I,1)+JTRA2*VBR(I,2)
334   CONTINUE
      ARG=RCZ
      DO 335 I=1,2
      ARG=ARG+AKBZ(I,IBZ)*TRAN(I)
335   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      IF(JDIR.EQ.0) THEN
      DO 3361 IQ=1,NLSQ
      DO 3362 JQ=1,NLSQ
      ZBSKI(JJ0+JQ,II0+IQ,IBZ) = ZBSKI(JJ0+JQ,II0+IQ,IBZ)
     &                         + ZEFA * BSTR(JQ,IQ,JCL,IB)
3362  CONTINUE
3361  CONTINUE
      END IF
      IF(JDIR.EQ.1) THEN
      DO 3371 IQ=1,NLSQ
      DO 3372 JQ=1,NLSQ
      ZBSKF(JJ0+JQ,II0+IQ,IBZ) = ZBSKF(JJ0+JQ,II0+IQ,IBZ)
     &                         + ZEFA * BSTR(JQ,IQ,JCL,IB)
3372  CONTINUE
3371  CONTINUE
      END IF
      IF(JDIR.EQ.-1) THEN
      DO 3381 IQ=1,NLSQ
      DO 3382 JQ=1,NLSQ
      ZBSKO(JJ0+JQ,II0+IQ,IBZ) = ZBSKO(JJ0+JQ,II0+IQ,IBZ)
     &                         + ZEFA * BSTR(JQ,IQ,JCL,IB)
3382  CONTINUE
3381  CONTINUE
      END IF
C
333   CONTINUE
332   CONTINUE
331   CONTINUE
C                       CHECK OF HERMITICITY
      DEV=RCZ
      DO 341 IBZ=1,NBZ
      DO 342 II=1,NBLSQ
      DO 343 JJ=II,NBLSQ
      ZW1=ZBSKI(II,JJ,IBZ)
      W1X=DBLE(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZBSKI(JJ,II,IBZ)
      W2X=DBLE(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
343   CONTINUE
342   CONTINUE
341   CONTINUE
      WRITE(IW6,140) DEV
140   FORMAT(/4X,'  BULK - IN-LAYER',6X,'DEVIATION=',G12.4)
C
      DEV=RCZ
      DO 346 IBZ=1,NBZ
      DO 347 II=1,NBLSQ
      DO 348 JJ=1,NBLSQ
      ZW1=ZBSKF(II,JJ,IBZ)
      W1X=DBLE(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZBSKO(JJ,II,IBZ)
      W2X=DBLE(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
348   CONTINUE
347   CONTINUE
346   CONTINUE
      WRITE(IW6,145) DEV
145   FORMAT(/4X,'  BULK - OFF-LAYER',5X,'DEVIATION=',G12.4)
C
C----------------------------------- VACUUM REGION
C
      DO 3501 IBZ=1,MNBZ
      DO 3502 II=1,MNBLSQ
      DO 3503 JJ=1,MNBLSQ
      ZVSKI(JJ,II,IBZ)=ZZ
      ZVSKO(JJ,II,IBZ)=ZZ
      ZVSKF(JJ,II,IBZ)=ZZ
3503  CONTINUE
3502  CONTINUE
3501  CONTINUE
C
      DO 351 IBZ=1,NBZ
      DO 352 IB=1,NB
      II0=(IB-1)*NLSQ
C
      DO 353 JCL=1,NSCLV(IB)
      JDIR=JDRPV(JCL,IB)
      JB=JBVAV(JCL,IB)
      JTRA1=JTRAV(1,JCL,IB)
      JTRA2=JTRAV(2,JCL,IB)
C
      JJ0=(JB-1)*NLSQ
      DO 354 I=1,2
      TRAN(I)=JTRA1*VBR(I,1)+JTRA2*VBR(I,2)
354   CONTINUE
      ARG=RCZ
      DO 355 I=1,2
      ARG=ARG+AKBZ(I,IBZ)*TRAN(I)
355   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      IF(JDIR.EQ.0) THEN
      DO 3561 IQ=1,NLSQ
      DO 3562 JQ=1,NLSQ
      ZVSKI(JJ0+JQ,II0+IQ,IBZ) = ZVSKI(JJ0+JQ,II0+IQ,IBZ)
     &                         + ZEFA * VSTR(JQ,IQ,JCL,IB)
3562  CONTINUE
3561  CONTINUE
      END IF
      IF(JDIR.EQ.1) THEN
      DO 3571 IQ=1,NLSQ
      DO 3572 JQ=1,NLSQ
      ZVSKF(JJ0+JQ,II0+IQ,IBZ) = ZVSKF(JJ0+JQ,II0+IQ,IBZ)
     &                         + ZEFA * VSTR(JQ,IQ,JCL,IB)
3572  CONTINUE
3571  CONTINUE
      END IF
      IF(JDIR.EQ.-1) THEN
      DO 3581 IQ=1,NLSQ
      DO 3582 JQ=1,NLSQ
      ZVSKO(JJ0+JQ,II0+IQ,IBZ) = ZVSKO(JJ0+JQ,II0+IQ,IBZ)
     &                         + ZEFA * VSTR(JQ,IQ,JCL,IB)
3582  CONTINUE
3581  CONTINUE
      END IF
C
353   CONTINUE
352   CONTINUE
351   CONTINUE
C                       CHECK OF HERMITICITY
      DEV=RCZ
      DO 361 IBZ=1,NBZ
      DO 362 II=1,NBLSQ
      DO 363 JJ=II,NBLSQ
      ZW1=ZVSKI(II,JJ,IBZ)
      W1X=DBLE(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZVSKI(JJ,II,IBZ)
      W2X=DBLE(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
363   CONTINUE
362   CONTINUE
361   CONTINUE
      WRITE(IW6,160) DEV
160   FORMAT(/4X,'VACUUM - IN-LAYER',6X,'DEVIATION=',G12.4)
C
      DEV=RCZ
      DO 366 IBZ=1,NBZ
      DO 367 II=1,NBLSQ
      DO 368 JJ=1,NBLSQ
      ZW1=ZVSKF(II,JJ,IBZ)
      W1X=DBLE(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZVSKO(JJ,II,IBZ)
      W2X=DBLE(ZW2)
      W2Y=DIMAG(ZW2)
      DEV1= ABS(W1X-W2X) +  ABS(W1Y+W2Y)
      DEV= MAX(DEV,DEV1)
368   CONTINUE
367   CONTINUE
366   CONTINUE
      WRITE(IW6,165) DEV
165   FORMAT(/4X,'VACUUM - OFF-LAYER',5X,'DEVIATION=',G12.4)
C
      RETURN
      END
C*******************
CXXX    MACO    ****
C*******************
      SUBROUTINE MACO
C
C***************************************************
C   SURFACE MADELUNG CONSTANTS
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
C
      DIMENSION POW(3,MNG),X(2),
     &          BREF(MNG),BREFD(MNG),
     &          VREF(MNG),VREFD(MNG)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/AMAD/ AMCMM(MNG,MNG),AMCMD(MNG,MNG),
     &             AMCDM(MNG,MNG),AMCDD(MNG,MNG),
     &             BARCM(MNG),BARCD(MNG)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
C
      DATA RCZ/0.0D0/,RCH/0.5D0/
C
      RNB=REAL(NB)
      NG=NB*NP
C
      IG=0
      DO 300 IP=1,NP
      DO 301 IB=1,NB
      IG=IG+1
      DO 302 I=1,3
      POW(I,IG)=POS(I,IB,IP)
302   CONTINUE
301   CONTINUE
300   CONTINUE
C
C------------------------------ AUXILIARY QUANTITIES
      DO 310 IG=1,NG
      BREF(IG)=RCZ
      BREFD(IG)=RCZ
      VREF(IG)=RCZ
      VREFD(IG)=RCZ
310   CONTINUE
C
      DO 311 IG=1,NG
      DO 312 IB=1,NB
C
      DO 314 I=1,2
      X(I)=BPOS(I,IB)-POW(I,IG)
314   CONTINUE
      XZ  =BPOS(3,IB)-POW(3,IG)
      CALL SEW(X,XZ,P,PD,PDD)
      BREF(IG)=BREF(IG)+P
      BREFD(IG)=BREFD(IG)+PD
C
      DO 315 I=1,2
      X(I)=VPOS(I,IB)-POW(I,IG)
315   CONTINUE
      XZ  =VPOS(3,IB)-POW(3,IG)
      CALL SEW(X,XZ,P,PD,PDD)
      VREF(IG)=VREF(IG)+P
      VREFD(IG)=VREFD(IG)+PD
C
312   CONTINUE
311   CONTINUE
C
      DO 318 IG=1,NG
      BREF(IG)=BREF(IG)/RNB
      BREFD(IG)=BREFD(IG)/RNB
      VREF(IG)=VREF(IG)/RNB
      VREFD(IG)=VREFD(IG)/RNB
318   CONTINUE
C
      FV=RCH*REAL(IVAC)
      FB=RCH*REAL(2-IVAC)
C
C------------------------------ MADELUNG CONSTANTS
      DO 320 JG=1,NG
      DO 321 IG=1,NG
C
      DO 324 I=1,2
      X(I)=POW(I,IG)-POW(I,JG)
324   CONTINUE
      XZ  =POW(3,IG)-POW(3,JG)
      CALL SEW(X,XZ,P,PD,PDD)
C
      AMCMM(IG,JG)= P - FB*BREF(JG) - FV*VREF(JG)
      AMCMD(IG,JG)= FB*BREFD(JG) + FV*VREFD(JG) - PD
      AMCDM(IG,JG)= PD
      AMCDD(IG,JG)= - PDD
321   CONTINUE
320   CONTINUE
C
C---------------------------- CONSTANTS FOR DIPOLE BARRIER
      DO 330 IG=1,NG
      BARCM(IG) = VREF(IG) - BREF(IG)
      BARCD(IG) = BREFD(IG) - VREFD(IG)
330   CONTINUE
C
C---------------------------- SCALING OF CONSTANTS BY THE
C                            DIMENSIONLESS BULK WS-RADIUS
      ABWS2=ABWS**2
      ABWS3=ABWS*ABWS2
      DO 340 JG=1,NG
      DO 341 IG=1,NG
      AMCMM(IG,JG)=ABWS*AMCMM(IG,JG)
      AMCMD(IG,JG)=ABWS2*AMCMD(IG,JG)
      AMCDM(IG,JG)=ABWS2*AMCDM(IG,JG)
      AMCDD(IG,JG)=ABWS3*AMCDD(IG,JG)
341   CONTINUE
340   CONTINUE
      DO 342 IG=1,NG
      BARCM(IG)=ABWS*BARCM(IG)
      BARCD(IG)=ABWS2*BARCD(IG)
342   CONTINUE
C
      WRITE(IW6,121)
121   FORMAT(//2X,' **********  DIMENSIONLESS MADELUNG ',
     &    'CONSTANTS  ********** '//'   IG,  JG,',
     & '     M-M,           M-D,',
     & '          D-M,           D-D:')
      DO 345 IG=1,NG
         DO 346 JG=1,NG
      WRITE(IW6,123) IG,JG,AMCMM(IG,JG),AMCMD(IG,JG),
     &                     AMCDM(IG,JG),AMCDD(IG,JG)
346   CONTINUE
345   CONTINUE
123   FORMAT(1X,I4,1X,I4,1X,4G15.7)
C
      WRITE(IW6,125)
125   FORMAT(/1X,' ******  DIMENSIONLESS CONSTANTS',
     &    ' FOR DIPOLE BARRIER  ******'/
     & /'       IG,     MONOPOLE,      DIPOLE: ')
      DO 348 IG=1,NG
      WRITE(IW6,127) IG,BARCM(IG),BARCD(IG)
348   CONTINUE
127   FORMAT(4X,I5,2X,2G15.7)
C
C---------------------------- SCALING OF CONSTANTS BY THE
C                             ACTUAL BULK WS-RADIUS
      BWST2=BWST**2
      BWST3=BWST*BWST2
      DO 350 JG=1,NG
      DO 351 IG=1,NG
      AMCMM(IG,JG)=AMCMM(IG,JG)/BWST
      AMCMD(IG,JG)=AMCMD(IG,JG)/BWST2
      AMCDM(IG,JG)=AMCDM(IG,JG)/BWST2
      AMCDD(IG,JG)=AMCDD(IG,JG)/BWST3
351   CONTINUE
350   CONTINUE
      DO 352 IG=1,NG
      BARCM(IG)=BARCM(IG)/BWST
      BARCD(IG)=BARCD(IG)/BWST2
352   CONTINUE
C
      RETURN
      END
C*******************
CXXX    SEW     ****
C*******************
      SUBROUTINE SEW(X,XZ,POT,POTD,POTDD)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------------------
C   ELECTROSTATIC POTENTIAL GENERATED BY A 2D-LATTICE
C   OF UNIT MONOPOLES. THE 2D-LATTICE LIES IN THE XY-PLANE.
C   THE POTENTIAL IS EVALUATED:
C   AT A POINT WELL OUTSIDE THE XY-PLANE USING SUBR. SEWOF,
C   AT A POINT NEAR THE XY-PLANE USING SUBR. SEWON.
C-----------------------------------------------------------
C   INPUT:
C     X(1),X(2),XZ - COORDINATES OF THE POINT:
C        X(1),X(2) - THE PROJECTION ONTO THE XY-PLANE
C               XZ - COORDINATE ALONG Z-AXIS
C   OUTPUT:
C      POT - THE ELECTROSTATIC POTENTIAL
C      POTD - 1ST DERIVATIVE OF POT ACCORDING TO XZ
C      POTDD - 2ND DERIVATIVE OF POT ACCORDING TO XZ
C----------------------------------------------------------
C  REMARK: COMMON/SUV/ MUST BE SET BEFORE !!!
C----------------------------------------------------------
C
      DIMENSION X(2),XR(2),C(2)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
C
      DATA RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
C
C                                      REDUCED X||
      PI=RC4*ATAN(RC1)
      PI2=RC2*PI
      DO 312 J=1,2
      COOR=(X(1)*VBG(1,J)+X(2)*VBG(2,J))/PI2
      C(J)= MOD(COOR,RC1)
312   CONTINUE
      DO 315 I=1,2
      XR(I)=C(1)*VBR(I,1)+C(2)*VBR(I,2)
315   CONTINUE
C
      IF(XZ**2.GE.ROMEGA) THEN
      CALL SEWOF(XR,XZ,POT,POTD,POTDD)
      ELSE
      CALL SEWON(XR,XZ,POT,POTD,POTDD)
      END IF
C
      RETURN
      END
C*******************
CXXX   SEWOF    ****
C*******************
      SUBROUTINE SEWOF(X,XZ,POT,POTD,POTDD)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------------------
C   ELECTROSTATIC POTENTIAL GENERATED BY A 2D-LATTICE
C   OF UNIT MONOPOLES:
C   FOR A POINT WELL OUTSIDE THE PLANE OF THE 2D-LATTICE
C-----------------------------------------------------------
C
      DIMENSION X(2),G(2)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
C
      DATA ARGMAX/50.0D0/
C
      ABSZ=ABS(XZ)
      GMAX=ARGMAX/ABSZ
C
      DUM=GMAX*SQRT(VBG(1,2)**2+VBG(2,2)**2)/GOMEGA
      MAXG1=1+INT(DUM)
      DUM=GMAX*SQRT(VBG(1,1)**2+VBG(2,1)**2)/GOMEGA
      MAXG2=1+INT(DUM)
C
      SUM=RCZ
      SUMD=RCZ
      SUMDD=RCZ
C
      DO 311 IG1=-MAXG1,MAXG1
      DO 312 IG2=0,MAXG2
      IF (IG1.LE.0.AND.IG2.EQ.0) GO TO 312
C
      DO 315 J=1,2
      G(J) = IG1*VBG(J,1) + IG2*VBG(J,2)
315   CONTINUE
C
      GG=SQRT(G(1)**2+G(2)**2)
      IF (GG.GT.GMAX) GO TO 312
C
      ALFA=G(1)*X(1)+G(2)*X(2)
      ARG=GG*ABSZ
      PROD=EXP(-ARG)*COS(ALFA)
C
      SUM   = SUM    + PROD/GG
      SUMD  = SUMD   + PROD
      SUMDD = SUMDD  + PROD*GG
C
312   CONTINUE
311   CONTINUE
C
      SUM   = RC2*SUM
      SUMD  = RC2*SUMD
      SUMDD = RC2*SUMDD
C
C
      PI=RC4*ATAN(RC1)
      FACTOR=RC4*PI/ROMEGA
C
      POT   = FACTOR * (-ABSZ + SUM)
      POTD  = FACTOR * (RC1 + SUMD)
      IF (XZ.GT.RCZ) POTD=-POTD
      POTDD = FACTOR * SUMDD
C
      RETURN
      END
C*******************
CXXX   SEWON    ****
C*******************
      SUBROUTINE SEWON(X,XZ,POT,POTD,POTDD)
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
C-----------------------------------------------------------
C   ELECTROSTATIC POTENTIAL GENERATED BY A 2D-LATTICE
C   OF UNIT MONOPOLES:
C   FOR A POINT NEAR OR IN THE PLANE OF THE 2D-LATTICE
C-----------------------------------------------------------
C
      DIMENSION X(2),R(2),G(2)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
C
      DATA ARGMAX/7.0D0/, CMUL/0.282D0/, TINY/1.0D-4/,
     &     BIG/150.0D0/
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/,
     &     RC3/3.0D0/,RC8/8.0D0/
C
C                                         CONSTANTS
      PI=RC4*ATAN(RC1)
      SPI=SQRT(PI)
      SIGMA=CMUL*SQRT(ROMEGA)
      TWOSIG=RC2*SIGMA
      SIGSPI=SIGMA*SPI
      TAU=XZ/TWOSIG
      TAU2=TAU**2
C
C-------------------------------- LONG-RANGE POTENTIAL
C
      QERF=RC1-ERRFC(TAU)
      QEXP=RCZ
      IF(TAU2.LT.BIG) QEXP=EXP(-TAU2)
C
      PLR = -(RC4*PI*XZ*QERF + RC8*SIGSPI*QEXP)/ROMEGA
      PLRD = -RC4*PI*QERF/ROMEGA
      PLRDD = -RC4*SPI*QEXP/(SIGMA*ROMEGA)
C
C--------------------------------- R-SUMMATION
C
      RMAX=ARGMAX*TWOSIG
      RMIN=TINY*SIGMA
      XZ2=XZ**2
      USIG2=RC1/SIGMA**2
C
      DUM=RMAX*SQRT(VBR(1,2)**2+VBR(2,2)**2)/ROMEGA
      MAXR1=2+INT(DUM)
      DUM=RMAX*SQRT(VBR(1,1)**2+VBR(2,1)**2)/ROMEGA
      MAXR2=2+INT(DUM)
C
      SUMR=RCZ
      SUMRD=RCZ
      SUMRDD=RCZ
C
      DO 301 IR1=-MAXR1,MAXR1
      DO 302 IR2=-MAXR2,MAXR2
C
      DO 305 J=1,2
      R(J) = IR1*VBR(J,1) + IR2*VBR(J,2) - X(J)
305   CONTINUE
C
      RR2=R(1)**2+R(2)**2+XZ2
      RR=SQRT(RR2)
      IF (RR.GT.RMAX) GO TO 302
      IF (RR.GT.RMIN) GO TO 220
C
      SUMR = SUMR - RC2/SIGSPI
      SUMRDD = SUMRDD + USIG2/(RC3*SIGSPI)
      GO TO 302
C
220   ARG=RR/TWOSIG
      ARG2=ARG**2
      QERFC=ERRFC(ARG)
      QEXP=EXP(-ARG2)
C
      DRDX=XZ/RR
      DRDX2=(RC1-XZ2/RR2)/RR
      DFDR=-RC2*(QERFC/RR2 + QEXP/(RR*SIGSPI))
      DFDR2=RC4*QERFC/(RR*RR2) + (RC4/RR2+USIG2)*QEXP/SIGSPI
      F=RC2*QERFC/RR
      DFDX=DFDR*DRDX
      DFDX2=DFDR2*DRDX**2 + DFDR*DRDX2
C
      SUMR = SUMR + F
      SUMRD = SUMRD + DFDX
      SUMRDD = SUMRDD + DFDX2
C
302   CONTINUE
301   CONTINUE
C
C--------------------------------- G-SUMMATION
C
      GMAX=(ARGMAX+ABS(TAU))/SIGMA
C
      DUM=GMAX*SQRT(VBG(1,2)**2+VBG(2,2)**2)/GOMEGA
      MAXG1=1+INT(DUM)
      DUM=GMAX*SQRT(VBG(1,1)**2+VBG(2,1)**2)/GOMEGA
      MAXG2=1+INT(DUM)
C
      SUMG=RCZ
      SUMGD=RCZ
      SUMGDD=RCZ
C
      DO 311 IG1=-MAXG1,MAXG1
      DO 312 IG2=0,MAXG2
      IF (IG1.LE.0.AND.IG2.EQ.0) GO TO 312
C
      DO 315 J=1,2
      G(J) = IG1*VBG(J,1) + IG2*VBG(J,2)
315   CONTINUE
C
      GG=SQRT(G(1)**2+G(2)**2)
      IF (GG.GT.GMAX) GO TO 312
C
      ALFA=G(1)*X(1)+G(2)*X(2)
      ARGEXP=GG*XZ
      SIGGG=SIGMA*GG
      ARGERP=SIGGG+TAU
      ARGERM=SIGGG-TAU
      ARGEXQ=SIGGG**2+TAU2
C
      COSALF=COS(ALFA)
      QEXPP=EXP(ARGEXP)
      QEXPM=EXP(-ARGEXP)
      QERFCP=ERRFC(ARGERP)
      QERFCM=ERRFC(ARGERM)
      QEXPQ=EXP(-ARGEXQ)
C
      PRODUP=QEXPP*QERFCP
      PRODUM=QEXPM*QERFCM
      SOUCET=PRODUP+PRODUM
      ROZDIL=PRODUP-PRODUM
C
      SUMG = SUMG + COSALF*SOUCET/GG
      SUMGD = SUMGD + COSALF*ROZDIL
      SUMGDD = SUMGDD + COSALF*(GG*SOUCET-RC2*QEXPQ/SIGSPI)
C
312   CONTINUE
311   CONTINUE
C
      SUMG = RC2*SUMG
      SUMGD = RC2*SUMGD
      SUMGDD = RC2*SUMGDD
C
      FACT = RC2*PI/ROMEGA
C
      SUMG = FACT*SUMG
      SUMGD = FACT*SUMGD
      SUMGDD = FACT*SUMGDD
C
C
      POT = PLR + SUMR + SUMG
      POTD = PLRD + SUMRD + SUMGD
      POTDD = PLRDD + SUMRDD + SUMGDD
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
CXXX    SCPF    ****
C*******************
      SUBROUTINE SCPF
C
C**************************************************
C   SET UP THE COHERENT POTENTIAL FUNCTIONS
C**************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNLSQ=MNL**2)
C
      DIMENSION OWORK(5),ZP(MNLSQ),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
C
      DATA RCZ/0.0D0/,RCH/0.5D0/,TOL/1.0D-6/
C
      NG=NP*NB
      NLSQ=NL**2
      ZZ=DCMPLX(RCZ,RCZ)
C
      IF(IITER.GT.1) GO TO 250
C
C******************************** FOR THE FIRST ITERATION
C
      DO 310 IS=1,MNS
      DO 311 IE=1,MNE
      DO 312 IG=1,MNG
      DO 313 J=1,MNLSQ
      DO 314 I=1,MNLSQ
      ZOMG(I,J,IG,IE,IS)=ZZ
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
      IF(IVAC.EQ.0) EFT=EFB
      IF(IVAC.EQ.1) EFT=RCH*(EFB+EFV)
C
      READ (IR4,100) OWORK
      READ (IR4,101) NGX,NLX,NSX,NEX
      READ (IR4,104) DIAMX,EFTX
      IDIF= ABS(NGX-NG)+ ABS(NLX-NL)
     &    + ABS(NSX-NS)+ ABS(NEX-NE)
      DIF = ABS(DIAMX-DIAM)+ ABS(EFTX-EFT)
      IF(IDIF.GT.0.OR.DIF.GT.TOL) GO TO 295
      DO 340 IS=1,NS
      DO 341 IE=1,NE
      DO 342 IG=1,NG
      READ (IR4,104) ((ZCP1(I,J),I=J,NLSQ),J=1,NLSQ)
      DO 343 J=1,NLSQ
      DO 344 I=J,NLSQ
      ZCP1(J,I)=ZCP1(I,J)
344   CONTINUE
343   CONTINUE
      DO 345 J=1,NLSQ
      DO 346 I=1,NLSQ
      ZCPF(I,J,IG,IE,IS)=ZCP1(I,J)
346   CONTINUE
345   CONTINUE
342   CONTINUE
341   CONTINUE
340   CONTINUE
      WRITE(IW6,140)
140   FORMAT(//5X,' *****  CPA STARTED FROM',
     &       ' A PREVIOUS RUN  *****'/)
      WRITE(IW6,100) OWORK
C
      RETURN
C
C********************************* FOR THE NEXT ITERATIONS
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
      DO 362 IG=1,NG
C
      DO 364 J=1,NLSQ
      DO 365 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IG,IE,IS)
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
      DO 370 IC=1,NC(IG)
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
      ZCPF(I,J,IG,IE,IS)=ZCP1(I,J)
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
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
C
      DATA RCZ/0.0D0/
C
      RCINT=RCZ
      DO 310 IE=1,NE
      ZD=ZCW(IE)*ZF(IE)
      RCINT=RCINT+DBLE(ZD)
310   CONTINUE
      RETURN
      END
C*******************
CXXX    SGFVA   ****
C*******************
      SUBROUTINE SGFVA
C
C***********************************************
C    UPDATE OF SGF FOR SEMIINFINITE VACUUM -
C     -  FOR EACH K||-POINT AND ENERGY NODE
C***********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNAV=1)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZVPF(MNLSQ),ZPSI1(MNBLSQ,MNBLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/VSCK/ ZVSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSV)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      IF (IITER.EQ.1) THEN
      DO 301 ISV=1,MNSV
      DO 302 IE=1,MNE
      DO 303 IBZ=1,MNBZ
      DO 304 JJ=1,MNBLSQ
      DO 305 II=1,MNBLSQ
      ZVGAM(II,JJ,IBZ,IE,ISV)=ZZ
305   CONTINUE
304   CONTINUE
303   CONTINUE
302   CONTINUE
301   CONTINUE
      END IF
C
C-------------------------------- ENERGY LOOP
      DO 320 IE=1,NE
      ZE=ZCN(IE)
C
      DO 310 IL=1,NL
      PPCX=VPPC(IL,1,1)
      PPDX=VPPD(IL,1,1)
      PPQX=VPPQ(IL,1,1)
      ALFX=QSCR(IL)
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      CALL PLMZ(PPCX,PPDX,PPQX,ALFX,ZE,ZPFX,ZLAX,ZMUX)
!SMP$ DO SERIAL
      DO 313 I=ISTA,IFIN
      ZVPF(I)=ZPFX
313   CONTINUE
310   CONTINUE
C                                   IBZ-LOOP
      DO 321 IBZ=1,NBZ
C
      DO 322 JJ=1,NBLSQ
      DO 323 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZVSKI(II,JJ,IBZ)
323   CONTINUE
322   CONTINUE
      DO 324 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 325 I=1,NLSQ
      ZPSI1(II0+I,II0+I)=ZPSI1(II0+I,II0+I) + ZVPF(I)
325   CONTINUE
324   CONTINUE
      DO 331 JJ=1,NBLSQ
      DO 332 II=1,NBLSQ
      ZSKO1(II,JJ)=ZVSKO(II,JJ,IBZ)
332   CONTINUE
331   CONTINUE
      DO 335 JJ=1,NBLSQ
      DO 336 II=1,NBLSQ
      ZSKF1(II,JJ)=ZVSKF(II,JJ,IBZ)
336   CONTINUE
335   CONTINUE
C
      DO 341 JJ=1,NBLSQ
      DO 342 II=1,NBLSQ
      ZWA(II,JJ)=ZVGAM(II,JJ,IBZ,IE,1)
342   CONTINUE
341   CONTINUE
C
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 345 JJ=1,NBLSQ
      DO 346 II=1,NBLSQ
      ZWB(II,JJ)=ZPSI1(II,JJ)-ZWA(II,JJ)
346   CONTINUE
345   CONTINUE
C
      CALL CINVC(ZWB,NBLSQ,MNBLSQ)
C
C                                  STORAGE OF SGF
      DO 348 JJ=1,NBLSQ
      DO 349 II=1,NBLSQ
      ZVGAM(II,JJ,IBZ,IE,1)=ZWB(II,JJ)
349   CONTINUE
348   CONTINUE
C
321   CONTINUE
320   CONTINUE
C
      RETURN
      END
C*******************
CXXX    REPA    ****
C*******************
      SUBROUTINE REPA
C
C***************************************************
C    RECURSION PARTITIONING FOR ON-LAYER BLOCKS
C                 OF THE GF MATRIX
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZPSI(MNBLSQ,MNBLSQ,MNP),ZWG(MNLSQ,MNLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZGMV(MNBLSQ,MNBLSQ,MNP),ZGMB(MNBLSQ,MNBLSQ,MNP),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSB)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSV)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      NG=NB*NP
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 IS=1,MNS
      DO 302 IE=1,MNE
      DO 303 IG=1,MNG
      DO 304 J=1,MNLSQ
      DO 305 I=1,MNLSQ
      ZGFG(I,J,IG,IE,IS)=ZZ
305   CONTINUE
304   CONTINUE
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
C--------------------------------- SPIN LOOP AND
C                                  ENERGY LOOP
      DO 310 IS=1,NS
      ISB= MIN(IS,NSB)
      ISV= MIN(IS,NSV)
      DO 311 IE=1,NE
C
C-------------------------------- IBZ LOOP
C
      DO 312 IBZ=1,NBZ
C                               IN-LAYER QUANTITIES
      DO 314 IP=1,NP
      DO 3158 JJ=1,NBLSQ
      DO 3159 II=1,NBLSQ
      ZPSI(II,JJ,IP)=-ZSKI(II,JJ,IP,IBZ)
3159  CONTINUE
3158  CONTINUE
314   CONTINUE
      IG=0
      DO 316 IP=1,NP
      DO 317 IB=1,NB
      IG=IG+1
      II0=(IB-1)*NLSQ
      DO 3181 J=1,NLSQ
      DO 3182 I=1,NLSQ
      ZPSI(II0+I,II0+J,IP) = ZPSI(II0+I,II0+J,IP)
     &                     + ZCPF(I,J,IG,IE,IS)
3182  CONTINUE
3181  CONTINUE
317   CONTINUE
316   CONTINUE
C
C                              RECURSION FROM VACUUM
      DO 3208 JJ=1,NBLSQ
      DO 3209 II=1,NBLSQ
      ZGMV(II,JJ,1)=ZVGAM(II,JJ,IBZ,IE,ISV)
3209  CONTINUE
3208  CONTINUE
C
      IF (IVAC.EQ.0) THEN
      DO 3801 JJ=1,NBLSQ
      DO 3802 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKO(II,JJ,1,IBZ)
3802  CONTINUE
3801  CONTINUE
      DO 3813 JJ=1,NBLSQ
      DO 3814 II=1,NBLSQ
      ZSKF1(II,JJ)=DCONJG(ZSKO1(JJ,II))
3814  CONTINUE
3813  CONTINUE
      DO 3825 JJ=1,NBLSQ
      DO 3826 II=1,NBLSQ
      ZWA(II,JJ)=ZGMV(II,JJ,1)
3826  CONTINUE
3825  CONTINUE
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3849 JJ=1,NBLSQ
      DO 3840 II=1,NBLSQ
      ZGMV(II,JJ,1)=ZWA(II,JJ)
3840  CONTINUE
3849  CONTINUE
      END IF
C
      IF(NP.GT.1) THEN
      DO 321 IP=2,NP
      IP1=IP-1
      DO 3228 JJ=1,NBLSQ
      DO 3229 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI(II,JJ,IP1)-ZGMV(II,JJ,IP1)
3229  CONTINUE
3228  CONTINUE
      DO 3248 JJ=1,NBLSQ
      DO 3249 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKF(II,JJ,IP1,IBZ)
3249  CONTINUE
3248  CONTINUE
      DO 3268 JJ=1,NBLSQ
      DO 3269 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKO(II,JJ,IP,IBZ)
3269  CONTINUE
3268  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3281 JJ=1,NBLSQ
      DO 3282 II=1,NBLSQ
      ZGMV(II,JJ,IP)=ZWA(II,JJ)
3282  CONTINUE
3281  CONTINUE
321   CONTINUE
      END IF
C
C                              RECURSION FROM BULK
      DO 3308 JJ=1,NBLSQ
      DO 3309 II=1,NBLSQ
      ZGMB(II,JJ,NP)=ZBGAM(II,JJ,IBZ,IE,ISB)
3309  CONTINUE
3308  CONTINUE
C
      IF(NP.GT.1) THEN
      DO 331 IP=NP-1,1,-1
      IP1=IP+1
      DO 3328 JJ=1,NBLSQ
      DO 3329 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI(II,JJ,IP1)-ZGMB(II,JJ,IP1)
3329  CONTINUE
3328  CONTINUE
      DO 3348 JJ=1,NBLSQ
      DO 3349 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKF(II,JJ,IP,IBZ)
3349  CONTINUE
3348  CONTINUE
      DO 3368 JJ=1,NBLSQ
      DO 3369 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKO(II,JJ,IP1,IBZ)
3369  CONTINUE
3368  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3381 JJ=1,NBLSQ
      DO 3382 II=1,NBLSQ
      ZGMB(II,JJ,IP)=ZWA(II,JJ)
3382  CONTINUE
3381  CONTINUE
331   CONTINUE
      END IF
C                                  ON-LAYER BLOCKS
C                               AND ON-SITE BLOCKS
      IG=0
      DO 340 IP=1,NP
      DO 3428 JJ=1,NBLSQ
      DO 3429 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI(II,JJ,IP)-ZGMV(II,JJ,IP)-ZGMB(II,JJ,IP)
3429  CONTINUE
3428  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 344 IB=1,NB
      IG=IG+1
      II0=(IB-1)*NLSQ
      DO 3461 J=1,NLSQ
      DO 3462 I=1,NLSQ
      ZGFG(I,J,IG,IE,IS) = ZGFG(I,J,IG,IE,IS)
     &            + WKBZ(IBZ) * ZWA(II0+I,II0+J)
3462  CONTINUE
3461  CONTINUE
344   CONTINUE
340   CONTINUE
C
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
      DO 352 IG=1,NG
C
      DO 3548 J=1,NLSQ
      DO 3549 I=1,NLSQ
      ZWG(I,J)=ZGFG(I,J,IG,IE,IS)
3549  CONTINUE
3548  CONTINUE
C
      CALL SYMGF(ZWG,NSYM,NL)
C
      DO 3568 J=1,NLSQ
      DO 3569 I=1,NLSQ
      ZGFG(I,J,IG,IE,IS)=ZWG(I,J)
3569  CONTINUE
3568  CONTINUE
C
352   CONTINUE
351   CONTINUE
350   CONTINUE
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
C*******************
CXXX    CMANO   ****
C*******************
      REAL*8 FUNCTION CMANO(ZA,N,ND)
C
C----------------------------------------
C       NORM OF A COMPLEX MATRIX
C----------------------------------------
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      DIMENSION ZA(ND,N)
C
      DATA RCZ/0.0D0/
C
      CMANO=RCZ
      DO 301 J=1,N
      DO 302 I=1,N
      XW=DBLE(ZA(I,J))
      YW=DIMAG(ZA(I,J))
      CMANO=CMANO+ ABS(XW)+ ABS(YW)
302   CONTINUE
301   CONTINUE
C
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
C
C                 FCC(001) SYMMETRY - GENERATED BY
C                 FOURFOLD ROTATION (AXIS Z) AND
C                 MIRROR REFLECTION (PLANE X-Z)
C
201   ZG(1,1)=ZP(1,1)
      ZG(3,1)=ZP(3,1)
      ZG(7,1)=ZP(7,1)
      ZG(2,2)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(6,2)=RCH*(ZP(6,2)+ZP(8,4))
      ZG(3,3)=ZP(3,3)
      ZG(7,3)=ZP(7,3)
      ZG(4,4)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(8,4)=RCH*(ZP(6,2)+ZP(8,4))
      ZG(5,5)=ZP(5,5)
      ZG(6,6)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(7,7)=ZP(7,7)
      ZG(8,8)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(9,9)=ZP(9,9)
       IF(NL.EQ.3) GO TO 270
      ZG(13,1)=ZP(13,1)
      ZG(10,2)=RCH*(ZP(10,2)-ZP(16,4))
      ZG(12,2)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(13,3)=ZP(13,3)
      ZG(14,4)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(16,4)=RCH*(-ZP(10,2)+ZP(16,4))
      ZG(11,5)=ZP(11,5)
      ZG(10,6)=RCH*(ZP(10,6)-ZP(16,8))
      ZG(12,6)=RCH*(ZP(12,6)+ZP(14,8))
      ZG(13,7)=ZP(13,7)
      ZG(14,8)=RCH*(ZP(12,6)+ZP(14,8))
      ZG(16,8)=RCH*(-ZP(10,6)+ZP(16,8))
      ZG(15,9)=ZP(15,9)
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
C                 BCC(110) SYMMETRY - GENERATED BY TWO
C                 MIRROR REFLECTIONS (PLANES X-Z AND Y-Z)
C
202   ZG(1,1)=ZP(1,1)
      ZG(3,1)=ZP(3,1)
      ZG(7,1)=ZP(7,1)
      ZG(9,1)=ZP(9,1)
      ZG(2,2)=ZP(2,2)
      ZG(6,2)=ZP(6,2)
      ZG(3,3)=ZP(3,3)
      ZG(7,3)=ZP(7,3)
      ZG(9,3)=ZP(9,3)
      ZG(4,4)=ZP(4,4)
      ZG(8,4)=ZP(8,4)
      ZG(5,5)=ZP(5,5)
      ZG(6,6)=ZP(6,6)
      ZG(7,7)=ZP(7,7)
      ZG(9,7)=ZP(9,7)
      ZG(8,8)=ZP(8,8)
      ZG(9,9)=ZP(9,9)
       IF(NL.EQ.3) GO TO 270
      ZG(13,1)=ZP(13,1)
      ZG(15,1)=ZP(15,1)
      ZG(10,2)=ZP(10,2)
      ZG(12,2)=ZP(12,2)
      ZG(13,3)=ZP(13,3)
      ZG(15,3)=ZP(15,3)
      ZG(14,4)=ZP(14,4)
      ZG(16,4)=ZP(16,4)
      ZG(11,5)=ZP(11,5)
      ZG(10,6)=ZP(10,6)
      ZG(12,6)=ZP(12,6)
      ZG(13,7)=ZP(13,7)
      ZG(15,7)=ZP(15,7)
      ZG(14,8)=ZP(14,8)
      ZG(16,8)=ZP(16,8)
      ZG(13,9)=ZP(13,9)
      ZG(15,9)=ZP(15,9)
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
C                 FCC(111) SYMMETRY - GENERATED BY
C                 THREEFOLD ROTATION (AXIS Z) AND
C                 MIRROR REFLECTION (PLANE Y-Z)
C
203   ZG(1,1)=ZP(1,1)
      ZG(3,1)=ZP(3,1)
      ZG(7,1)=ZP(7,1)
      ZG(2,2)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(6,2)=RCH*(ZP(6,2)+ZP(8,4))
      ZG(9,2)=RCH*(ZP(9,2)+ZP(5,4))
      ZG(3,3)=ZP(3,3)
      ZG(7,3)=ZP(7,3)
      ZG(4,4)=RCH*(ZP(2,2)+ZP(4,4))
      ZG(5,4)=RCH*(ZP(9,2)+ZP(5,4))
      ZG(8,4)=RCH*(ZP(6,2)+ZP(8,4))
      ZG(5,5)=RCH*(ZP(5,5)+ZP(9,9))
      ZG(8,5)=RCH*(ZP(8,5)+ZP(6,9))
      ZG(6,6)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(9,6)=RCH*(ZP(9,6)+ZP(5,8))
      ZG(7,7)=ZP(7,7)
      ZG(8,8)=RCH*(ZP(6,6)+ZP(8,8))
      ZG(9,9)=RCH*(ZP(5,5)+ZP(9,9))
       IF(NL.EQ.3) GO TO 270
      ZG(10,1)=ZP(10,1)
      ZG(13,1)=ZP(13,1)
      ZG(12,2)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(15,2)=RCH*(ZP(15,2)+ZP(11,4))
      ZG(10,3)=ZP(10,3)
      ZG(13,3)=ZP(13,3)
      ZG(11,4)=RCH*(ZP(15,2)+ZP(11,4))
      ZG(14,4)=RCH*(ZP(12,2)+ZP(14,4))
      ZG(11,5)=RCH*(ZP(11,5)+ZP(15,9))
      ZG(14,5)=RCH*(ZP(14,5)+ZP(12,9))
      ZG(12,6)=RCH*(ZP(12,6)+ZP(14,8))
      ZG(15,6)=RCH*(ZP(15,6)+ZP(11,8))
      ZG(10,7)=ZP(10,7)
      ZG(13,7)=ZP(13,7)
      ZG(11,8)=RCH*(ZP(15,6)+ZP(11,8))
      ZG(14,8)=RCH*(ZP(12,6)+ZP(14,8))
      ZG(12,9)=RCH*(ZP(14,5)+ZP(12,9))
      ZG(15,9)=RCH*(ZP(11,5)+ZP(15,9))
      ZG(10,10)=ZP(10,10)
      ZG(13,10)=ZP(13,10)
      ZG(11,11)=RCH*(ZP(11,11)+ZP(15,15))
      ZG(14,11)=RCH*(ZP(14,11)+ZP(12,15))
      ZG(12,12)=RCH*(ZP(12,12)+ZP(14,14))
      ZG(15,12)=RCH*(ZP(15,12)+ZP(11,14))
      ZG(13,13)=ZP(13,13)
      ZG(14,14)=RCH*(ZP(12,12)+ZP(14,14))
      ZG(15,15)=RCH*(ZP(11,11)+ZP(15,15))
      ZG(16,16)=ZP(16,16)
       GO TO 270
C
C                 GB (N01) SYMMETRY - GENERATED BY ONE
C                 MIRROR REFLECTION (PLANE X-Z)
C
204   ZG(1,1)=ZP(1,1)
      ZG(3,1)=ZP(3,1)
      ZG(4,1)=ZP(4,1)
      ZG(7,1)=ZP(7,1)
      ZG(8,1)=ZP(8,1)
      ZG(9,1)=ZP(9,1)
      ZG(2,2)=ZP(2,2)
      ZG(5,2)=ZP(5,2)
      ZG(6,2)=ZP(6,2)
      ZG(3,3)=ZP(3,3)
      ZG(4,3)=ZP(4,3)
      ZG(7,3)=ZP(7,3)
      ZG(8,3)=ZP(8,3)
      ZG(9,3)=ZP(9,3)
      ZG(4,4)=ZP(4,4)
      ZG(7,4)=ZP(7,4)
      ZG(8,4)=ZP(8,4)
      ZG(9,4)=ZP(9,4)
      ZG(5,5)=ZP(5,5)
      ZG(6,5)=ZP(6,5)
      ZG(6,6)=ZP(6,6)
      ZG(7,7)=ZP(7,7)
      ZG(8,7)=ZP(8,7)
      ZG(9,7)=ZP(9,7)
      ZG(8,8)=ZP(8,8)
      ZG(9,8)=ZP(9,8)
      ZG(9,9)=ZP(9,9)
       IF(NL.EQ.3) GO TO 270
      ZG(13,1)=ZP(13,1)
      ZG(14,1)=ZP(14,1)
      ZG(15,1)=ZP(15,1)
      ZG(16,1)=ZP(16,1)
      ZG(10,2)=ZP(10,2)
      ZG(11,2)=ZP(11,2)
      ZG(12,2)=ZP(12,2)
      ZG(13,3)=ZP(13,3)
      ZG(14,3)=ZP(14,3)
      ZG(15,3)=ZP(15,3)
      ZG(16,3)=ZP(16,3)
      ZG(13,4)=ZP(13,4)
      ZG(14,4)=ZP(14,4)
      ZG(15,4)=ZP(15,4)
      ZG(16,4)=ZP(16,4)
      ZG(10,5)=ZP(10,5)
      ZG(11,5)=ZP(11,5)
      ZG(12,5)=ZP(12,5)
      ZG(10,6)=ZP(10,6)
      ZG(11,6)=ZP(11,6)
      ZG(12,6)=ZP(12,6)
      ZG(13,7)=ZP(13,7)
      ZG(14,7)=ZP(14,7)
      ZG(15,7)=ZP(15,7)
      ZG(16,7)=ZP(16,7)
      ZG(13,8)=ZP(13,8)
      ZG(14,8)=ZP(14,8)
      ZG(15,8)=ZP(15,8)
      ZG(16,8)=ZP(16,8)
      ZG(13,9)=ZP(13,9)
      ZG(14,9)=ZP(14,9)
      ZG(15,9)=ZP(15,9)
      ZG(16,9)=ZP(16,9)
      ZG(10,10)=ZP(10,10)
      ZG(11,10)=ZP(11,10)
      ZG(12,10)=ZP(12,10)
      ZG(11,11)=ZP(11,11)
      ZG(12,11)=ZP(12,11)
      ZG(12,12)=ZP(12,12)
      ZG(13,13)=ZP(13,13)
      ZG(14,13)=ZP(14,13)
      ZG(15,13)=ZP(15,13)
      ZG(16,13)=ZP(16,13)
      ZG(14,14)=ZP(14,14)
      ZG(15,14)=ZP(15,14)
      ZG(16,14)=ZP(16,14)
      ZG(15,15)=ZP(15,15)
      ZG(16,15)=ZP(16,15)
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
CXXX    CPAIT   ****
C*******************
      SUBROUTINE CPAIT
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNLSQ=MNL**2)
C
      DIMENSION NONC(MNG,MNS),ZP(MNLSQ),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNG,MNE,MNS)
      COMMON/CAGF/ ZCAGF(MNLSQ,MNLSQ,MNA,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
C
      DATA RC1/1.0D0/,RCH/0.5D0/
C
      NLSQ=NL**2
      NG=NP*NB
C
      QMIXI=QMIXC
      IF(IITER.EQ.1) QMIXI=RC1
C
      DO 301 IS=1,MNS
      DO 302 IG=1,MNG
      NONC(IG,IS)=0
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
      DO 312 IG=1,NG
C
      DO 314 J=1,NLSQ
      DO 315 I=1,NLSQ
      ZFI1(I,J)=ZGFG(I,J,IG,IE,IS)
315   CONTINUE
314   CONTINUE
      DO 317 J=1,NLSQ
      DO 318 I=1,NLSQ
      ZCP1(I,J)=ZCPF(I,J,IG,IE,IS)
318   CONTINUE
317   CONTINUE
      DO 320 J=1,NLSQ
      DO 321 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IG,IE,IS)
321   CONTINUE
320   CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
C
      DO 323 J=1,NLSQ
      DO 324 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
324   CONTINUE
323   CONTINUE
C
      SUM=CMANO(ZW,NLSQ,MNLSQ)
      SUM1=CMANO(ZOM1,NLSQ,MNLSQ)
      IF(SUM.GE.(DIFMC*SUM1)) NONC(IG,IS)=NONC(IG,IS)+1
C
      DO 326 J=1,NLSQ
      DO 327 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXI*ZW(I,J)
327   CONTINUE
326   CONTINUE
C
      CALL SYMGF(ZOM1,NSYM,NL)
C
C..............................  LOOP OVER COMPONENTS
C
      DO 330 IC=1,NC(IG)
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
      ZOMG(I,J,IG,IE,IS)=ZOM1(I,J)
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
      WRITE(IW6,111) IITER,IITER
111   FORMAT(//5X,'-----------  ITERATION:  ITER=',
     &       I4,'  -----------'/
     &  /'  ****  CONVERGENCE OF COHERENT INTERACTOR :'/
     &  '        IG,      NON-CONV. OMEGA:',
     &   20X,'ITER=',I4)
      DO 360 IG=1,NG
      WRITE(IW6,120) IG,(NONC(IG,IS),IS=1,NS)
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
      IF(IVAC.EQ.0) EFT=EFB
      IF(IVAC.EQ.1) EFT=RCH*(EFB+EFV)
C
      REWIND IW8
      WRITE(IW8,181)
181   FORMAT(2X,' ----------  CPF-FILE  ----------')
      WRITE(IW8,101) NG,NL,NS,NE
      WRITE(IW8,104) DIAM,EFT
      DO 380 IS=1,NS
      DO 381 IE=1,NE
      DO 382 IG=1,NG
      WRITE(IW8,104) ((ZCPF(I,J,IG,IE,IS),I=J,NLSQ),J=1,NLSQ)
382   CONTINUE
381   CONTINUE
380   CONTINUE
      WRITE(IW8,182)
182   FORMAT(2X,' --------------------------------')
C
      RETURN
      END
C*******************
CXXX    CMOM    ****
C*******************
      SUBROUTINE CMOM
C
C************************************************
C   CALCULATES MOMENTS OF DOS VIA COMPLEX
C   CONTOUR INTEGRATION OF CONFIGURATIONALLY
C   AVERAGED PHYSICAL GREENS FUNCTIONS
C************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE=14)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MPAIR=(MNL-1)**2)
C
      DIMENSION ZGDI(MNE,MNL,MNS,MNA),
     &          ZGOF(MNE,MPAIR,MNS,MNA),
     &          ZW0(MNE),ZW1(MNE),ZW2(MNE),
     &          ZW00(MNE),ZW01(MNE),ZW10(MNE),
     &          ZW11(MNE),ZW02(MNE),ZW20(MNE)
C
      COMMON/CAGF/ ZCAGF(MNLSQ,MNLSQ,MNA,MNE,MNS)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA),
     &             EMOF00(MPAIR,MNS,MNA),EMOF10(MPAIR,MNS,MNA),
     &             EMOF01(MPAIR,MNS,MNA),EMOF11(MPAIR,MNS,MNA),
     &             EMOF20(MPAIR,MNS,MNA),EMOF02(MPAIR,MNS,MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/SCREEN/ QSCR(MNL)
C
      DATA RCZ/0.0D0/
C
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
C----------------------------- OFF-DIAGONAL ELEMENTS
C
      DO 320 IA=1,NA
      DO 321 IS=1,NS
      IPAIR=0
      DO 322 IL1=1,NL-1
      IL2=IL1+1
      L1=IL1-1
      L2=IL1
      ICEN1=L1**2+L1+1
      ICEN2=L2**2+L2+1
      PPC1=PPC(IL1,IS,IA)
      PPD1=PPD(IL1,IS,IA)
      PPQ1=PPQ(IL1,IS,IA)
      ALF1=QSCR(IL1)
      PPC2=PPC(IL2,IS,IA)
      PPD2=PPD(IL2,IS,IA)
      PPQ2=PPQ(IL2,IS,IA)
      ALF2=QSCR(IL2)
      DO 323 M=-L1,L1
      IPAIR=IPAIR+1
      I1=ICEN1+M
      I2=ICEN2+M
      DO 325 IE=1,NE
      ZE=ZCN(IE)
      CALL PLMZ(PPC1,PPD1,PPQ1,ALF1,ZE,ZPF1,ZLA1,ZMU1)
      CALL PLMZ(PPC2,PPD2,PPQ2,ALF2,ZE,ZPF2,ZLA2,ZMU2)
      ZGOF(IE,IPAIR,IS,IA) =
     &         ZMU1 * ZCAGF(I1,I2,IA,IE,IS) * ZMU2
325   CONTINUE
323   CONTINUE
322   CONTINUE
321   CONTINUE
320   CONTINUE
C
C---------------------------------- INTEGRATIONS
C
      DO 350 IA=1,NA
      DO 351 IS=1,NS
      DO 352 IL=1,NL
      ENY1=ENY(IL,IS,IA)
      DO 355 IE=1,NE
      ZEPS=ZCN(IE)-ENY1
      ZW0(IE)=ZGDI(IE,IL,IS,IA)
      ZW1(IE)=ZW0(IE)*ZEPS
      ZW2(IE)=ZW1(IE)*ZEPS
355   CONTINUE
      EMDI0(IL,IS,IA)=RCINT(ZW0)
      EMDI1(IL,IS,IA)=RCINT(ZW1)
      EMDI2(IL,IS,IA)=RCINT(ZW2)
352   CONTINUE
351   CONTINUE
350   CONTINUE
C
      DO 360 IA=1,NA
      DO 361 IS=1,NS
      IPAIR=0
      DO 362 IL1=1,NL-1
      IL2=IL1+1
      L1=IL1-1
      ENY1=ENY(IL1,IS,IA)
      ENY2=ENY(IL2,IS,IA)
      DO 363 M=-L1,L1
      IPAIR=IPAIR+1
      DO 365 IE=1,NE
      ZEPS1=ZCN(IE)-ENY1
      ZEPS2=ZCN(IE)-ENY2
      ZW00(IE)=ZGOF(IE,IPAIR,IS,IA)
      ZW10(IE)=ZEPS1*ZW00(IE)
      ZW01(IE)=ZEPS2*ZW00(IE)
      ZW11(IE)=ZEPS1*ZW01(IE)
      ZW20(IE)=ZEPS1*ZW10(IE)
      ZW02(IE)=ZEPS2*ZW01(IE)
365   CONTINUE
      EMOF00(IPAIR,IS,IA)=RCINT(ZW00)
      EMOF10(IPAIR,IS,IA)=RCINT(ZW10)
      EMOF01(IPAIR,IS,IA)=RCINT(ZW01)
      EMOF11(IPAIR,IS,IA)=RCINT(ZW11)
      EMOF20(IPAIR,IS,IA)=RCINT(ZW20)
      EMOF02(IPAIR,IS,IA)=RCINT(ZW02)
363   CONTINUE
362   CONTINUE
361   CONTINUE
360   CONTINUE
C
      RETURN
      END
C*******************
CXXX    GABU    ****
C*******************
      SUBROUTINE GABU
C
C***************************************************
C    GAMMA FOR THE SEMIINFINITE BULK RANDOM ALLOY -
C     -  FOR EACH K||-POINT AND ENERGY NODE
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNAB=20)
      PARAMETER(MNSB=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZPF(MNLSQ,MNAB),ZCP(MNLSQ,MNLSQ,MNB),
     &          ZOM(MNLSQ,MNLSQ,MNB),ZGF(MNLSQ,MNLSQ,MNB),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      DIMENSION ZPSI1(MNBLSQ,MNBLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZGMB1(MNBLSQ,MNBLSQ),ZGMB(MNBLSQ,MNBLSQ,MNBZ),
     &          ZGMV1(MNBLSQ,MNBLSQ),ZGMV(MNBLSQ,MNBLSQ,MNBZ),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/BSCK/ ZBSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSB)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/BPOPA/ BENY(MNL,MNSB,MNAB),BPPC(MNL,MNSB,MNAB),
     &              BPPD(MNL,MNSB,MNAB),BPPQ(MNL,MNSB,MNAB),
     &              BPPP(MNL,MNSB,MNAB),BDNY(MNL,MNSB,MNAB),
     &              BFINY(MNL,MNSB,MNAB),BFINYD(MNL,MNSB,MNAB)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      WRITE(IW6,111)
111   FORMAT(//5X,'**** CONVERGENCE OF BULK GAMMA:'/)
C
C------------------------------- LOOP OVER SPIN AND ENERGY
C
      DO 310 IS=1,NSB
      DO 311 IE=1,NE
C
      ZE=ZCN(IE)
C
      ITS=0
C
      DO 301 IBZ=1,MNBZ
      DO 3018 JJ=1,MNBLSQ
      DO 3019 II=1,MNBLSQ
      ZGMB(II,JJ,IBZ)=ZZ
3019  CONTINUE
3018  CONTINUE
301   CONTINUE
      DO 302 IBZ=1,MNBZ
      DO 3028 JJ=1,MNBLSQ
      DO 3029 II=1,MNBLSQ
      ZGMV(II,JJ,IBZ)=ZZ
3029  CONTINUE
3028  CONTINUE
302   CONTINUE
C
210   ITS=ITS+1
      IF(ITS.GT.NMITS) GO TO 291
C
C------------------------------------ CPA-PART
C
      NONC=NB
      IF(ITS.GT.1) GO TO 212
C..................................  1ST ITERATION
C
      DO 321 IA=1,NAB
      DO 322 IL=1,NL
      PPCX=BPPC(IL,IS,IA)
      PPDX=BPPD(IL,IS,IA)
      PPQX=BPPQ(IL,IS,IA)
      ALFX=QSCR(IL)
      CALL PLMZ(PPCX,PPDX,PPQX,ALFX,ZE,ZPFX,ZLAX,ZMUX)
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 323 I=ISTA,IFIN
      ZPF(I,IA)=ZPFX
323   CONTINUE
322   CONTINUE
321   CONTINUE
C
      DO 305 IB=1,MNB
      DO 3058 J=1,MNLSQ
      DO 3059 I=1,MNLSQ
      ZOM(I,J,IB)=ZZ
3059  CONTINUE
3058  CONTINUE
305   CONTINUE
      DO 306 IB=1,MNB
      DO 3068 J=1,MNLSQ
      DO 3069 I=1,MNLSQ
      ZCP(I,J,IB)=ZZ
3069  CONTINUE
3068  CONTINUE
306   CONTINUE
C
      DO 325 I=1,NLSQ
      ZDUM=ZZ
      DO 326 IA=1,NAB
      ZDUM=ZDUM+BCON(IA)/ZPF(I,IA)
326   CONTINUE
      ZDUM=RC1/ZDUM
      DO 327 IB=1,NB
      ZCP(I,I,IB)=ZDUM
327   CONTINUE
325   CONTINUE
      GO TO 220
C................................... NEXT ITERATIONS
C
212   DO 330 IB=1,NB
C
      DO 3318 J=1,NLSQ
      DO 3319 I=1,NLSQ
      ZFI1(I,J)=ZGF(I,J,IB)
3319  CONTINUE
3318  CONTINUE
      DO 3328 J=1,NLSQ
      DO 3329 I=1,NLSQ
      ZCP1(I,J)=ZCP(I,J,IB)
3329  CONTINUE
3328  CONTINUE
      DO 3338 J=1,NLSQ
      DO 3339 I=1,NLSQ
      ZOM1(I,J)=ZOM(I,J,IB)
3339  CONTINUE
3338  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3358 J=1,NLSQ
      DO 3359 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
3359  CONTINUE
3358  CONTINUE
C
      SUM=CMANO(ZW,NLSQ,MNLSQ)
      SUM1=CMANO(ZOM1,NLSQ,MNLSQ)
      IF(SUM.LT.(DIFMC*SUM1)) NONC=NONC-1
      DO 3371 J=1,NLSQ
      DO 3372 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXC*ZW(I,J)
3372  CONTINUE
3371  CONTINUE
      CALL SYMGF(ZOM1,NSYM,NL)
C
      DO 3381 J=1,NLSQ
      DO 3382 I=1,NLSQ
      ZFI1(I,J)=ZZ
3382  CONTINUE
3381  CONTINUE
C
      DO 340 IA=1,NAB
      DO 3418 J=1,NLSQ
      DO 3419 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3419  CONTINUE
3418  CONTINUE
      DO 3425 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZPF(I,IA)
3425  CONTINUE
      CALL CINVC(ZW,NLSQ,MNLSQ)
      DO 3438 J=1,NLSQ
      DO 3439 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+BCON(IA)*ZW(I,J)
3439  CONTINUE
3438  CONTINUE
340   CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3458 J=1,NLSQ
      DO 3459 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3459  CONTINUE
3458  CONTINUE
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3471 J=1,NLSQ
      DO 3472 I=1,NLSQ
      ZCP(I,J,IB)=ZCP1(I,J)
3472  CONTINUE
3471  CONTINUE
      DO 3481 J=1,NLSQ
      DO 3482 I=1,NLSQ
      ZOM(I,J,IB)=ZOM1(I,J)
3482  CONTINUE
3481  CONTINUE
C
330   CONTINUE
C
C---------------------------------------- SGF-PART
C
220   NONS=2*NBZ
C
      DO 308 IB=1,MNB
      DO 3081 J=1,MNLSQ
      DO 3082 I=1,MNLSQ
      ZGF(I,J,IB)=ZZ
3082  CONTINUE
3081  CONTINUE
308   CONTINUE
C...................................... IBZ-LOOP
      DO 350 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3528 JJ=1,NBLSQ
      DO 3529 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZBSKI(II,JJ,IBZ)
3529  CONTINUE
3528  CONTINUE
      DO 354 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3558 J=1,NLSQ
      DO 3559 I=1,NLSQ
      ZPSI1(II0+I,II0+J)=ZPSI1(II0+I,II0+J) + ZCP(I,J,IB)
3559  CONTINUE
3558  CONTINUE
354   CONTINUE
C                              TRANSFER MATRICES
      DO 3571 JJ=1,NBLSQ
      DO 3572 II=1,NBLSQ
      ZSKF1(II,JJ)=ZBSKF(II,JJ,IBZ)
3572  CONTINUE
3571  CONTINUE
      DO 3581 JJ=1,NBLSQ
      DO 3582 II=1,NBLSQ
      ZSKO1(II,JJ)=ZBSKO(II,JJ,IBZ)
3582  CONTINUE
3581  CONTINUE
C
C                              ITERATION OF BULK-SIDE GAMMA
      DO 3618 JJ=1,NBLSQ
      DO 3619 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB(II,JJ,IBZ)
3619  CONTINUE
3618  CONTINUE
C
      DO 3628 JJ=1,NBLSQ
      DO 3629 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMB1(II,JJ)
3629  CONTINUE
3628  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3661 JJ=1,NBLSQ
      DO 3662 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMB1(II,JJ)
3662  CONTINUE
3661  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMB1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3681 JJ=1,NBLSQ
      DO 3682 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB1(II,JJ)+QMIXS*ZWA(II,JJ)
3682  CONTINUE
3681  CONTINUE
      DO 3691 JJ=1,NBLSQ
      DO 3692 II=1,NBLSQ
      ZGMB(II,JJ,IBZ)=ZGMB1(II,JJ)
3692  CONTINUE
3691  CONTINUE
C
C                           ITERATION OF VACUUM-SIDE GAMMA
      DO 3718 JJ=1,NBLSQ
      DO 3719 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV(II,JJ,IBZ)
3719  CONTINUE
3718  CONTINUE
C
      DO 3728 JJ=1,NBLSQ
      DO 3729 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)
3729  CONTINUE
3728  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3761 JJ=1,NBLSQ
      DO 3762 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMV1(II,JJ)
3762  CONTINUE
3761  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMV1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3781 JJ=1,NBLSQ
      DO 3782 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV1(II,JJ)+QMIXS*ZWA(II,JJ)
3782  CONTINUE
3781  CONTINUE
      DO 3791 JJ=1,NBLSQ
      DO 3792 II=1,NBLSQ
      ZGMV(II,JJ,IBZ)=ZGMV1(II,JJ)
3792  CONTINUE
3791  CONTINUE
C
C                                  ON-SITE GF-BLOCKS
      DO 3828 JJ=1,NBLSQ
      DO 3829 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)-ZGMB1(II,JJ)
3829  CONTINUE
3828  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 384 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3861 J=1,NLSQ
      DO 3862 I=1,NLSQ
      ZGF(I,J,IB) = ZGF(I,J,IB) + WKBZ(IBZ) * ZWA(II0+I,II0+J)
3862  CONTINUE
3861  CONTINUE
384   CONTINUE
C............................... END OF IBZ-LOOP
350   CONTINUE
C
      DO 387 IB=1,NB
      DO 3881 J=1,NLSQ
      DO 3882 I=1,NLSQ
      ZW(I,J)=ZGF(I,J,IB)
3882  CONTINUE
3881  CONTINUE
      CALL SYMGF(ZW,NSYM,NL)
      DO 3891 J=1,NLSQ
      DO 3892 I=1,NLSQ
      ZGF(I,J,IB)=ZW(I,J)
3892  CONTINUE
3891  CONTINUE
387   CONTINUE
C
      IF((NONS+NONC).GT.0) GO TO 210
C
      WRITE(IW6,150) IS,IE,ITS
150   FORMAT(10X,'IS=',I2,5X,'IE=',I3,8X,I4,' ITER.')
C
C----------------------------- CONNECTION TO THE
C                              INTERMEDIATE REGION
      DO 390 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3918 JJ=1,NBLSQ
      DO 3919 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB(II,JJ,IBZ)
3919  CONTINUE
3918  CONTINUE
      DO 3928 JJ=1,NBLSQ
      DO 3929 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZBSKI(II,JJ,IBZ)
3929  CONTINUE
3928  CONTINUE
      DO 394 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3958 J=1,NLSQ
      DO 3959 I=1,NLSQ
      ZPSI1(II0+I,II0+J)=ZPSI1(II0+I,II0+J) + ZCP(I,J,IB)
3959  CONTINUE
3958  CONTINUE
394   CONTINUE
C                              TRANSFER MATRICES
      DO 3961 JJ=1,NBLSQ
      DO 3962 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKF(II,JJ,NP,IBZ)
3962  CONTINUE
3961  CONTINUE
      DO 3971 JJ=1,NBLSQ
      DO 3972 II=1,NBLSQ
      ZSKO1(II,JJ)=DCONJG(ZSKF1(JJ,II))
3972  CONTINUE
3971  CONTINUE
C                                     INVERSION AND
C                                    MULTIPLICATIONS
      DO 3981 JJ=1,NBLSQ
      DO 3982 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMB1(II,JJ)
3982  CONTINUE
3981  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3991 JJ=1,NBLSQ
      DO 3992 II=1,NBLSQ
      ZBGAM(II,JJ,IBZ,IE,IS)=ZWA(II,JJ)
3992  CONTINUE
3991  CONTINUE
C
390   CONTINUE
C
C----------------------------- END OF ENERGY AND SPIN LOOP
311   CONTINUE
310   CONTINUE
C
      RETURN
C
C------------------------------ NO CONVERGENCY
C
291   WRITE(IW6,191) NMITS,IE,IS
191   FORMAT(//' ****  GABU:  NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   '      ---  FOR IE=',I3,'    IS=',I2)
      STOP
      END
C*******************
CXXX    GAVA    ****
C*******************
      SUBROUTINE GAVA
C
C***************************************************
C    GAMMA FOR THE SEMIINFINITE BULK RANDOM ALLOY -
C     -  FOR EACH K||-POINT AND ENERGY NODE
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNE=14)
      PARAMETER(MNBZ=700)
      PARAMETER(MNAV=1)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZPF(MNLSQ,MNAV),ZCP(MNLSQ,MNLSQ,MNB),
     &          ZOM(MNLSQ,MNLSQ,MNB),ZGF(MNLSQ,MNLSQ,MNB),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      DIMENSION ZPSI1(MNBLSQ,MNBLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZGMB1(MNBLSQ,MNBLSQ),ZGMB(MNBLSQ,MNBLSQ,MNBZ),
     &          ZGMV1(MNBLSQ,MNBLSQ),ZGMV(MNBLSQ,MNBLSQ,MNBZ),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/VSCK/ ZVSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNE,MNSV)
      COMMON/CNW/ ZCN(MNE),ZCW(MNE),NE
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,ICONT
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      WRITE(IW6,111)
111   FORMAT(//5X,'**** CONVERGENCE OF VACUUM GAMMA:'/)
C
C------------------------------- LOOP OVER SPIN AND ENERGY
C
      DO 310 IS=1,NSV
      DO 311 IE=1,NE
C
      ZE=ZCN(IE)
C
      ITS=0
C
      DO 301 IBZ=1,MNBZ
      DO 3018 JJ=1,MNBLSQ
      DO 3019 II=1,MNBLSQ
      ZGMB(II,JJ,IBZ)=ZZ
3019  CONTINUE
3018  CONTINUE
301   CONTINUE
      DO 302 IBZ=1,MNBZ
      DO 3028 JJ=1,MNBLSQ
      DO 3029 II=1,MNBLSQ
      ZGMV(II,JJ,IBZ)=ZZ
3029  CONTINUE
3028  CONTINUE
302   CONTINUE
C
210   ITS=ITS+1
      IF(ITS.GT.NMITS) GO TO 291
C
C------------------------------------ CPA-PART
C
      NONC=NB
      IF(ITS.GT.1) GO TO 212
C..................................  1ST ITERATION
C
      DO 321 IA=1,NAV
      DO 322 IL=1,NL
      PPCX=VPPC(IL,IS,IA)
      PPDX=VPPD(IL,IS,IA)
      PPQX=VPPQ(IL,IS,IA)
      ALFX=QSCR(IL)
      CALL PLMZ(PPCX,PPDX,PPQX,ALFX,ZE,ZPFX,ZLAX,ZMUX)
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 323 I=ISTA,IFIN
      ZPF(I,IA)=ZPFX
323   CONTINUE
322   CONTINUE
321   CONTINUE
C
      DO 305 IB=1,MNB
      DO 3058 J=1,MNLSQ
      DO 3059 I=1,MNLSQ
      ZOM(I,J,IB)=ZZ
3059  CONTINUE
3058  CONTINUE
305   CONTINUE
      DO 306 IB=1,MNB
      DO 3068 J=1,MNLSQ
      DO 3069 I=1,MNLSQ
      ZCP(I,J,IB)=ZZ
3069  CONTINUE
3068  CONTINUE
306   CONTINUE
C
      DO 325 I=1,NLSQ
      ZDUM=ZZ
      DO 326 IA=1,NAV
      ZDUM=ZDUM+VCON(IA)/ZPF(I,IA)
326   CONTINUE
      ZDUM=RC1/ZDUM
      DO 327 IB=1,NB
      ZCP(I,I,IB)=ZDUM
327   CONTINUE
325   CONTINUE
      GO TO 220
C................................... NEXT ITERATIONS
C
212   DO 330 IB=1,NB
C
      DO 3318 J=1,NLSQ
      DO 3319 I=1,NLSQ
      ZFI1(I,J)=ZGF(I,J,IB)
3319  CONTINUE
3318  CONTINUE
      DO 3328 J=1,NLSQ
      DO 3329 I=1,NLSQ
      ZCP1(I,J)=ZCP(I,J,IB)
3329  CONTINUE
3328  CONTINUE
      DO 3338 J=1,NLSQ
      DO 3339 I=1,NLSQ
      ZOM1(I,J)=ZOM(I,J,IB)
3339  CONTINUE
3338  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3358 J=1,NLSQ
      DO 3359 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
3359  CONTINUE
3358  CONTINUE
C
      SUM=CMANO(ZW,NLSQ,MNLSQ)
      SUM1=CMANO(ZOM1,NLSQ,MNLSQ)
      IF(SUM.LT.(DIFMC*SUM1)) NONC=NONC-1
      DO 3371 J=1,NLSQ
      DO 3372 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXC*ZW(I,J)
3372  CONTINUE
3371  CONTINUE
      CALL SYMGF(ZOM1,NSYM,NL)
C
      DO 3381 J=1,NLSQ
      DO 3382 I=1,NLSQ
      ZFI1(I,J)=ZZ
3382  CONTINUE
3381  CONTINUE
C
      DO 340 IA=1,NAV
      DO 3418 J=1,NLSQ
      DO 3419 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3419  CONTINUE
3418  CONTINUE
      DO 3425 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZPF(I,IA)
3425  CONTINUE
      CALL CINVC(ZW,NLSQ,MNLSQ)
      DO 3438 J=1,NLSQ
      DO 3439 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+VCON(IA)*ZW(I,J)
3439  CONTINUE
3438  CONTINUE
340   CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3458 J=1,NLSQ
      DO 3459 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3459  CONTINUE
3458  CONTINUE
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3471 J=1,NLSQ
      DO 3472 I=1,NLSQ
      ZCP(I,J,IB)=ZCP1(I,J)
3472  CONTINUE
3471  CONTINUE
      DO 3481 J=1,NLSQ
      DO 3482 I=1,NLSQ
      ZOM(I,J,IB)=ZOM1(I,J)
3482  CONTINUE
3481  CONTINUE
C
330   CONTINUE
C
C---------------------------------------- SGF-PART
C
220   NONS=2*NBZ
C
      DO 308 IB=1,MNB
      DO 3081 J=1,MNLSQ
      DO 3082 I=1,MNLSQ
      ZGF(I,J,IB)=ZZ
3082  CONTINUE
3081  CONTINUE
308   CONTINUE
C...................................... IBZ-LOOP
      DO 350 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3528 JJ=1,NBLSQ
      DO 3529 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZVSKI(II,JJ,IBZ)
3529  CONTINUE
3528  CONTINUE
      DO 354 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3558 J=1,NLSQ
      DO 3559 I=1,NLSQ
      ZPSI1(II0+I,II0+J)=ZPSI1(II0+I,II0+J) + ZCP(I,J,IB)
3559  CONTINUE
3558  CONTINUE
354   CONTINUE
C                              TRANSFER MATRICES
      DO 3571 JJ=1,NBLSQ
      DO 3572 II=1,NBLSQ
      ZSKF1(II,JJ)=ZVSKF(II,JJ,IBZ)
3572  CONTINUE
3571  CONTINUE
      DO 3581 JJ=1,NBLSQ
      DO 3582 II=1,NBLSQ
      ZSKO1(II,JJ)=ZVSKO(II,JJ,IBZ)
3582  CONTINUE
3581  CONTINUE
C
C                              ITERATION OF BULK-SIDE GAMMA
      DO 3618 JJ=1,NBLSQ
      DO 3619 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB(II,JJ,IBZ)
3619  CONTINUE
3618  CONTINUE
C
      DO 3628 JJ=1,NBLSQ
      DO 3629 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMB1(II,JJ)
3629  CONTINUE
3628  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3661 JJ=1,NBLSQ
      DO 3662 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMB1(II,JJ)
3662  CONTINUE
3661  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMB1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3681 JJ=1,NBLSQ
      DO 3682 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB1(II,JJ)+QMIXS*ZWA(II,JJ)
3682  CONTINUE
3681  CONTINUE
      DO 3691 JJ=1,NBLSQ
      DO 3692 II=1,NBLSQ
      ZGMB(II,JJ,IBZ)=ZGMB1(II,JJ)
3692  CONTINUE
3691  CONTINUE
C
C                           ITERATION OF VACUUM-SIDE GAMMA
      DO 3718 JJ=1,NBLSQ
      DO 3719 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV(II,JJ,IBZ)
3719  CONTINUE
3718  CONTINUE
C
      DO 3728 JJ=1,NBLSQ
      DO 3729 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)
3729  CONTINUE
3728  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3761 JJ=1,NBLSQ
      DO 3762 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMV1(II,JJ)
3762  CONTINUE
3761  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMV1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3781 JJ=1,NBLSQ
      DO 3782 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV1(II,JJ)+QMIXS*ZWA(II,JJ)
3782  CONTINUE
3781  CONTINUE
      DO 3791 JJ=1,NBLSQ
      DO 3792 II=1,NBLSQ
      ZGMV(II,JJ,IBZ)=ZGMV1(II,JJ)
3792  CONTINUE
3791  CONTINUE
C
C                                  ON-SITE GF-BLOCKS
      DO 3828 JJ=1,NBLSQ
      DO 3829 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)-ZGMB1(II,JJ)
3829  CONTINUE
3828  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 384 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3861 J=1,NLSQ
      DO 3862 I=1,NLSQ
      ZGF(I,J,IB) = ZGF(I,J,IB) + WKBZ(IBZ) * ZWA(II0+I,II0+J)
3862  CONTINUE
3861  CONTINUE
384   CONTINUE
C............................... END OF IBZ-LOOP
350   CONTINUE
C
      DO 387 IB=1,NB
      DO 3881 J=1,NLSQ
      DO 3882 I=1,NLSQ
      ZW(I,J)=ZGF(I,J,IB)
3882  CONTINUE
3881  CONTINUE
      CALL SYMGF(ZW,NSYM,NL)
      DO 3891 J=1,NLSQ
      DO 3892 I=1,NLSQ
      ZGF(I,J,IB)=ZW(I,J)
3892  CONTINUE
3891  CONTINUE
387   CONTINUE
C
      IF((NONS+NONC).GT.0) GO TO 210
C
      WRITE(IW6,150) IS,IE,ITS
150   FORMAT(10X,'IS=',I2,5X,'IE=',I3,8X,I4,' ITER.')
C
C----------------------------- CONNECTION TO THE
C                              INTERMEDIATE REGION
      DO 390 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3918 JJ=1,NBLSQ
      DO 3919 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV(II,JJ,IBZ)
3919  CONTINUE
3918  CONTINUE
      DO 3928 JJ=1,NBLSQ
      DO 3929 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZVSKI(II,JJ,IBZ)
3929  CONTINUE
3928  CONTINUE
      DO 394 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3958 J=1,NLSQ
      DO 3959 I=1,NLSQ
      ZPSI1(II0+I,II0+J)=ZPSI1(II0+I,II0+J) + ZCP(I,J,IB)
3959  CONTINUE
3958  CONTINUE
394   CONTINUE
C                              TRANSFER MATRICES
      DO 3961 JJ=1,NBLSQ
      DO 3962 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKO(II,JJ,1,IBZ)
3962  CONTINUE
3961  CONTINUE
      DO 3971 JJ=1,NBLSQ
      DO 3972 II=1,NBLSQ
      ZSKF1(II,JJ)=DCONJG(ZSKO1(JJ,II))
3972  CONTINUE
3971  CONTINUE
C                                     INVERSION AND
C                                    MULTIPLICATIONS
      DO 3981 JJ=1,NBLSQ
      DO 3982 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)
3982  CONTINUE
3981  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3991 JJ=1,NBLSQ
      DO 3992 II=1,NBLSQ
      ZVGAM(II,JJ,IBZ,IE,IS)=ZWA(II,JJ)
3992  CONTINUE
3991  CONTINUE
C
390   CONTINUE
C
C----------------------------- END OF ENERGY AND SPIN LOOP
311   CONTINUE
310   CONTINUE
C
      RETURN
C
C------------------------------ NO CONVERGENCY
C
291   WRITE(IW6,191) NMITS,IE,IS
191   FORMAT(//' ****  GAVA:  NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   '      ---  FOR IE=',I3,'    IS=',I2)
      STOP
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
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
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
      DATA RCZ/0.0D0/,RC3/3.0D0/
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
C
      IF(IVAC.EQ.0) THEN
      J=J+1
      SP(J)=REAL(NS)*VWST**3/RC3
      END IF
      NAM=J
C
      WRITE(IW6,150) NAM
150   FORMAT(//5X,'*** LENGTH OF VECTORS FOR ANDERSON',
     &            ' MIXING= ',I8)
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
161   FORMAT(/4X,'IA=',I4,8X,A16,8X,'IS=',I2/
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNA=100)
      PARAMETER(MNS=2)
      PARAMETER(MNR=400)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/OPT/ POT(MNR,MNS,MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MPAIR=(MNL-1)**2)
      PARAMETER(MNR=400)
C
      COMMON/WAF/ PHI(MNR,MNL,MNS,MNA),
     &            PHID(MNR,MNL,MNS,MNA),
     &            PHIDD(MNR,MNL,MNS,MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA),
     &             EMOF00(MPAIR,MNS,MNA),EMOF10(MPAIR,MNS,MNA),
     &             EMOF01(MPAIR,MNS,MNA),EMOF11(MPAIR,MNS,MNA),
     &             EMOF20(MPAIR,MNS,MNA),EMOF02(MPAIR,MNS,MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
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
CXXX    DIMO    ****
C*******************
      SUBROUTINE DIMO
C
C*****************************************
C   CALCULATES THE DIPOLE MOMENTS
C   FROM THE L-OFF-DIAGONAL MOMENTS OF DOS
C*****************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MPAIR=(MNL-1)**2)
      PARAMETER(MNR=400)
C
      DIMENSION UINT(MPAIR),DIPMS(MNS),YNAM(MNL,MNS,MNA),
     &          RCUB(MNR),QQ(MNR)
C
      COMMON/RVW/ R(MNR),V(MNR),VI(9,MNR),WG(MNR),WF(MNR),
     &            NSIRK
      COMMON/WAF/ PHI(MNR,MNL,MNS,MNA),
     &            PHID(MNR,MNL,MNS,MNA),
     &            PHIDD(MNR,MNL,MNS,MNA)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA),
     &             EMOF00(MPAIR,MNS,MNA),EMOF10(MPAIR,MNS,MNA),
     &             EMOF01(MPAIR,MNS,MNA),EMOF11(MPAIR,MNS,MNA),
     &             EMOF20(MPAIR,MNS,MNA),EMOF02(MPAIR,MNS,MNA)
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),DIPMOM(MNA),
     &            VMAD(MNA),DMAD(MNA),ACTE(MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/
      DATA RCH/0.5D0/,RC3/3.0D0/
C
      PI=RC4*ATAN(RC1)
C                                  ANGULAR INTEGRALS
C                                  INCLUDING FACTOR 2
      PREF=RC2*SQRT(RC4*PI/RC3)
      IPAIR=0
      DO 302 IL1=1,NL-1
      IL2=IL1+1
      L1=IL1-1
      L2=IL1
      ICEN1=L1**2+L1+1
      ICEN2=L2**2+L2+1
!SMP$ DO SERIAL
      DO 303 M=-L1,L1
      IPAIR=IPAIR+1
      I1=ICEN1+M
      I2=ICEN2+M
      UINT(IPAIR)=PREF*GFRH(I2,I1,3)
303   CONTINUE
302   CONTINUE
C                           SIGNS FOR SQRT(PDOT) FUNCTIONS
      DO 306 IA=1,NA
      DO 307 IS=1,NS
      DO 308 IL=1,NL
      ARG = RC1/FINY(IL,IS,IA)
     &    - WS(IA)*FINYD(IL,IS,IA)*(DNY(IL,IS,IA)+REAL(IL))
      YNAM(IL,IS,IA)=SIGN(RC1,ARG)
308   CONTINUE
307   CONTINUE
306   CONTINUE
C
C                                    LOOP OVER ATOMS
      DO 310 IA=1,NA
C                                     RADIAL MESH
      NR1=NSZRAD(IA)
      WS1=WS(IA)
      CALL RAPO(NR1,WS1,R)
C
!SMP$ DO SERIAL
      DO 312 I=1,NR1
      RCUB(I)=R(I)**3
312   CONTINUE
C                                   LOOP OVER SPINS
      DO 315 IS=1,NS
C
      DIPMS(IS)=RCZ
      IPAIR=0
C                              LOOP OVER (L1,L2)-PAIRS
      DO 320 IL1=1,NL-1
      IL2=IL1+1
      L1=IL1-1
      DO 321 M=-L1,L1
      IPAIR=IPAIR+1
C
      WM00=EMOF00(IPAIR,IS,IA)
      WM10=EMOF10(IPAIR,IS,IA)
      WM01=EMOF01(IPAIR,IS,IA)
      WM11=EMOF11(IPAIR,IS,IA)
       HM20=RCH*EMOF20(IPAIR,IS,IA)
       HM02=RCH*EMOF02(IPAIR,IS,IA)
C
!SMP$ DO SERIAL
      DO 325 I=1,NR1
C
      SUM = PHI(I,IL1,IS,IA)   * PHI(I,IL2,IS,IA)   * WM00
     &    + PHID(I,IL1,IS,IA)  * PHID(I,IL2,IS,IA)  * WM11
     &    + PHI(I,IL1,IS,IA)   * PHID(I,IL2,IS,IA)  * WM01
     &    + PHID(I,IL1,IS,IA)  * PHI(I,IL2,IS,IA)   * WM10
     &    + PHI(I,IL1,IS,IA)   * PHIDD(I,IL2,IS,IA) *  HM02
     &    + PHIDD(I,IL1,IS,IA) * PHI(I,IL2,IS,IA)   *  HM20
C
      QQ(I)=RCUB(I)*SUM
325   CONTINUE
C
      RINT=QUAD3(NR1,R,QQ)
C
      DIPMS(IS) = DIPMS(IS)
     &   +  RINT*UINT(IPAIR)*YNAM(IL1,IS,IA)*YNAM(IL2,IS,IA)
C
321   CONTINUE
320   CONTINUE
C
315   CONTINUE
      DIPMOM(IA) = DIPMS(1) + DIPMS(NS)
310   CONTINUE
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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNAV=1)
      PARAMETER(MNSV=2)
      PARAMETER(MPAIR=(MNL-1)**2)
      PARAMETER(MNR=400)
      PARAMETER(MNAM=MNA*MNS*MNR)
      PARAMETER(MNUH=48)
C
      DIMENSION YMAX(2,MNA),IMAX(2,MNA),
     &          CHTRG(MNG),DIPMG(MNG),VMADG(MNG),DMADG(MNG),
     &          WDUM(MNR),RHOT(MNR),RHOS(MNR,MNS),VXC(2),
     &          VHART(MNR),FPOT(MNR,MNS,MNA)
C
      COMMON/AMAD/ AMCMM(MNG,MNG),AMCMD(MNG,MNG),
     &             AMCDM(MNG,MNG),AMCDD(MNG,MNG),
     &             BARCM(MNG),BARCD(MNG)
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
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA),
     &             EMOF00(MPAIR,MNS,MNA),EMOF10(MPAIR,MNS,MNA),
     &             EMOF01(MPAIR,MNS,MNA),EMOF11(MPAIR,MNS,MNA),
     &             EMOF20(MPAIR,MNS,MNA),EMOF02(MPAIR,MNS,MNA)
      COMMON/DIBA/ DBA
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),DIPMOM(MNA),
     &            VMAD(MNA),DMAD(MNA),ACTE(MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/AMXP/ DXP(MNAM,MNUH)
      COMMON/AMFP/ DFP(MNAM,MNUH)
      COMMON/AMSP/ SP(MNAM),XL(MNAM),FL(MNAM),XN(MNAM),
     &             VOMA(MNUH,MNUH)
      COMMON/TEXT/ OTXTA(MNA)
      COMMON/SCREEN/ QSCR(MNL)      
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC4/4.0D0/,
     &     RCH/0.5D0/
C
      DATA AZL/2.5D0/, DMUL/0.03D0/
C
C------------------------------------------------------------
C    WITH ANDERSON MIXING FOR THE POTENTIALS
C    (V. EYERT, J. COMPUT. PHYS. 124 (1996) 271)
C------------------------------------------------------------
C
      PI=RC4*ATAN(RC1)
      PI4=RC4*PI
      NG=NB*NP
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
C                        SITE-RESOLVED CHARGE TRANSFERS
C                              AND DIPOLE MOMENTS
      IA=0
      DO 303 IG=1,NG
      CHTRG(IG)=RCZ
      DIPMG(IG)=RCZ
      DO 304 IC=1,NC(IG)
      IA=IA+1
      CHTRG(IG)=CHTRG(IG)+CON(IA)*CHATRA(IA)
      DIPMG(IG)=DIPMG(IG)+CON(IA)*DIPMOM(IA)
304   CONTINUE
303   CONTINUE
C                           SITE-RESOLVED MADELUNG TERMS
      DO 305 IG=1,NG
      VMADG(IG)=RCZ
      DMADG(IG)=RCZ
305   CONTINUE
      DO 306 JG=1,NG
      DO 307 IG=1,NG
      VMADG(IG)=VMADG(IG) + AMCMM(IG,JG)*CHTRG(JG)
     &                    + AMCMD(IG,JG)*DIPMG(JG)
      DMADG(IG)=DMADG(IG) + AMCDM(IG,JG)*CHTRG(JG)
     &                    + AMCDD(IG,JG)*DIPMG(JG)
307   CONTINUE
306   CONTINUE
C                                   ATOMIC MADELUNG TERMS
      IA=0
      DO 308 IG=1,NG
      DO 309 IC=1,NC(IG)
      IA=IA+1
      VMAD(IA)=VMADG(IG)
      DMAD(IA)=DMADG(IG)
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
315   CONTINUE
C
      DUM = - RHOUP*POT(I,1,IA) - RHODO*POT(I,NS,IA)
     &      + RHOT(I)*(VNUC+RCH*VHART(I)+EXC)
      WDUM(I)=DUM*R(I)**2
314   CONTINUE
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
     &    +RCH*DIPMOM(IA)*DMAD(IA)
C
      ACTE(IA)=SCOR+SVAL+COTE+ELST
C
310   CONTINUE
C
C                                    DIPOLE BARRIER
      FDBA=RCZ
      DO 319 IG=1,NG
      FDBA=FDBA+BARCM(IG)*CHTRG(IG)+BARCD(IG)*DIPMG(IG)
319   CONTINUE
C
C------------------------------ DEVIATION OF THE POTENTIALS
C                                    AND THE DIPOLE BARRIER
      FDBA=FDBA-DBA
      WRITE(IW6,125) FDBA,IITER
125   FORMAT(/' ***  BARRIER DEVIATION=',G13.3,10X,'ITER=',I4)
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
C----------------------------   MIXING OF POTENTIALS
C                                AND DIPOLE BARRIER
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
      IF(IVAC.EQ.0) XL(NAM)=DBA
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
      IF(IVAC.EQ.0) FL(NAM)=FDBA
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
      IF(IVAC.EQ.0) DBA=XN(NAM)
C
C--------------------------------- POTENTIAL PARAMETERS
C                                 OF THE NEW POTENTIALS
C
C                                 VACUUM REGION
      IF(IVAC.EQ.0) THEN
      S2=VWST**2
      DO 340 IL=1,NL
      L=IL-1
      TL1=REAL(2*L+1)
      TL5=REAL(2*L+5)
      VENY(IL,1,1)=DBA
      VPPC(IL,1,1)=DBA+TL1*TL5/(RC2*S2)
340   CONTINUE
      END IF
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
C  call write here!!!
      call write_atoms(NA,NL,NS,NSZRAD,EFW,AZ,WS,WSAV,
     & QSCR,POT,OTXTA,DBA,IVAC)

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
      PARAMETER(MNP=30)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=100)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MPAIR=(MNL-1)**2)
      PARAMETER(MNR=400)
C
      DIMENSION QT(MNA),QP(MNL,MNA),AMP(MNL,MNA),
     &          BHFC(MNA),BHFV(MNA),BHFT(MNA),
     &          DNC(MNA),DNV(MNA),DNT(MNA)
C
      COMMON/COR/ ECOR(20,2,MNA),THRESH,NCOR(20,MNA),
     &            LCOR(20,MNA),NOBC(20,MNA),NUMCOR(MNA)
      COMMON/RHO/ RHOCOR(MNR,MNS,MNA),RHOVAL(MNR,MNS,MNA)
      COMMON/POPA/ ENY(MNL,MNS,MNA),PPC(MNL,MNS,MNA),
     &             PPD(MNL,MNS,MNA),PPQ(MNL,MNS,MNA),
     &             PPP(MNL,MNS,MNA),DNY(MNL,MNS,MNA),
     &             FINY(MNL,MNS,MNA),FINYD(MNL,MNS,MNA)
      COMMON/EMOM/ EMDI0(MNL,MNS,MNA),EMDI1(MNL,MNS,MNA),
     &             EMDI2(MNL,MNS,MNA),
     &             EMOF00(MPAIR,MNS,MNA),EMOF10(MPAIR,MNS,MNA),
     &             EMOF01(MPAIR,MNS,MNA),EMOF11(MPAIR,MNS,MNA),
     &             EMOF20(MPAIR,MNS,MNA),EMOF02(MPAIR,MNS,MNA)
      COMMON/DIBA/ DBA
      COMMON/CMD/ CHATRA(MNA),AMGMOM(MNA),DIPMOM(MNA),
     &            VMAD(MNA),DMAD(MNA),ACTE(MNA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA),NSZRAD(MNA)
      COMMON/DATSUB/ DIAM,EFB,EFV,BWST,VWST,EFW
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/ITLDA/ ALFA,BETA,W0AM,NITER,NITERA,NUH,NAM,
     &              IITER,NFUPR
      COMMON/ITLDI/ ISHENY(MNL,MNS,MNA),IREL,IVXC
      COMMON/TEXT/ OTXTA(MNA)
C
      DATA RCZ/0.0D0/,RCH/0.5D0/,DLIM/0.02D0/,BCONST/52.4308D0/
C
C------------------------------ SHORT PRINT IN EACH ITERATION
C
C                                ATOMIC CHARGES AND MOMENTS
      WRITE(IW6,110) IITER
110   FORMAT(/1X,'----- ATOM, CHARGE TRANSF.,',
     &  ' MAGN. MOMENT, DIPOLE MOMENT:        ITER=',I4)
      DO 311 IA=1,NA
      WRITE(IW6,111) IA,CHATRA(IA),AMGMOM(IA),DIPMOM(IA)
311   CONTINUE
111   FORMAT(5X,I5,3F15.5)
C
C                              TOTAL CHARGE, MAGNETIC MOMENT,
C                             DIPOLE BARRIER, WORK FUNCTION,
C                               AND TOTAL ENERGY
      QTOT=RCZ
      AMTOT=RCZ
      ETOT=RCZ
      DO 316 IA=1,NA
      QTOT=QTOT+CON(IA)*CHATRA(IA)
      AMTOT=AMTOT+CON(IA)*AMGMOM(IA)
      ETOT=ETOT+CON(IA)*ACTE(IA)
316   CONTINUE
      IF(IVAC.EQ.0) WRKF=DBA-EFB
      IF(IVAC.EQ.1) WRKF=RCZ
C
      WRITE(IW6,116) QTOT,AMTOT,IITER
116   FORMAT(/2X,'TOTAL CHARGE=',F12.6,
     &        5X,'TOTAL MAGMOM=',F11.5,5X,'ITER=',I4)
      WRITE(IW6,117) DBA,WRKF,IITER
117   FORMAT(/2X,'DIP. BARRIER=',F11.5,
     &        5X,'WORK FUNCTION=',F11.5,5X,'ITER=',I4)
      WRITE(IW6,118) ETOT,IITER
118   FORMAT(/2X,'TOTAL ENERGY=',G20.12,20X,'ITER=',I4)
C
C------------------------------- POLES OF THE SCREENED
C                                 POTENTIAL FUNCTIONS
      IF(IVAC.EQ.0) EFT=EFB
      IF(IVAC.EQ.1) EFT=RCH*(EFB+EFV)
      DO 320 IA=1,NA
      DO 321 IS=1,NS
      DO 322 IL=1,NL
      POLE=PPC(IL,IS,IA)-PPD(IL,IS,IA)/(PPQ(IL,IS,IA)-QSCR(IL))
      IF(ABS(POLE-EFT).GT.DLIM) GO TO 322
      WRITE(IW6,127) POLE,IA,IS,IL
127   FORMAT(1X,'  * POLE IN E= ',G15.7,'    FOR',
     &          '  IA=',I4,'  IS=',I1,'  IL=',I1)
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
C
C                                ATOMIC MADELUNG TERMS AND
C                              CONTRIBUTIONS TO TOTAL ENERGY
      WRITE(IW6,112) IITER
112   FORMAT(/1X,'----- ATOM, MADELUNG TERM, DIP. MAD. TERM,',
     &  '  ENERGY TERM:      ITER=',I4)
      DO 313 IA=1,NA
      WRITE(IW6,113) IA,VMAD(IA),DMAD(IA),ACTE(IA)
313   CONTINUE
113   FORMAT(5X,I5,2X,2G15.7,2X,G20.12)
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
      PARAMETER(MNA=100)
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
      PARAMETER(MNA=100)
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
