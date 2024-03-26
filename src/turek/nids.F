C*******************
CXXX    IDS     ****
C*******************
      PROGRAM IDS
C
C*******************************************************
C      PROGRAM FOR CALCULATION OF
C       INTERFACE DENSITY OF STATES AND
C       THE BLOCH SPECTRAL FUNCTIONS
C       BY MEANS OF THE TB-LMTO-SGF-CPA METHOD
C*******************************************************
C      VERSION:  DECEMBER 1999
C*******************************************************
C   BASIC FEATURES OF THE PROGRAM ARE:
C     - ELECTRONIC STRUCTURE OF ORDERED AND DISORDERED
C       SURFACES AND INTERFACES OF PURE METALS AND
C       DISORDERED ALLOYS
C     - CASES: FCC(001), FCC(111), FCC(110), HCP(0001),
C              BCC(110), BCC(001) AND SIMILAR
C     - GEOMETRY: INTERFACE NORMAL LIES IN Z AXIS,
C       FROZEN RELAXATIONS ARE ALLOWED
C     - NON-MAGNETIC AND SPIN-POLARIZED
C       (WITH COLLINEAR SPIN STRUCTURES ONLY)
C     - SUBSTITUTIONAL CHEMICAL DISORDER
C     - ATOMIC SPHERE APPROXIMATION
C     - SPD AND SPDF CASES
C*******************************************************
C   INPUT:   IR1  (= 1)  -  GEOMETRY OF THE INTERFACE
C            IR2  (= 2)  -  CHEMICAL OCCUPATIONS
C                           OF THE SITES
C            IR5  (= 5)  -  CONTROL PARAMETERS
C            IR11 (=11)  -  1ST SUBSTRATE POT. PARAMETERS,
C                           COMPOSITION, AND FERMI ENERGY
C   OUTPUT:  IW6  (= 6)  -  CURRENT OUTPUT
C            IW7  (= 7)  -  ORBITAL-RESOLVED DOS
C            IW8  (= 8)  -  ATOM-RESOLVED DOS FOR SPIN UP
C            IW9  (= 9)  -  ATOM-RESOLVED DOS FOR SPIN DOWN
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
C    NE ...... NUMBER OF ENERGY POINTS
C   MNE ...... PARAMETER: MAX. NE
C    NBZ ..... NUMBER OF K||-POINTS IN THE
C              IRREDUCIBLE PART OF 2D-BZ
C   MNBZ ..... PARAMETER: MAX. NBZ
C    NBZA .... NUMBER OF K||-POINTS FOR THE BSF
C   MNBZA .... PARAMETER: MAX. NBZA
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
C******************************************************
C
      PARAMETER(MNP=12)
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZ=128)
      PARAMETER(MNBZA=4)
      PARAMETER(MNAB=2)
      PARAMETER(MNAV=2)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
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
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/KMESHA/ AKBZA(2,MNBZA),NBZA
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
      COMMON/SCKA/ ZSKIA(MNBLSQ,MNBLSQ,MNP,MNBZA),
     &             ZSKOA(MNBLSQ,MNBLSQ,MNP,MNBZA),
     &             ZSKFA(MNBLSQ,MNBLSQ,MNP,MNBZA)
      COMMON/BSCKA/ ZBSKIA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZBSKOA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZBSKFA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/VSCKA/ ZVSKIA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZVSKOA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZVSKFA(MNBLSQ,MNBLSQ,MNBZA)
C
C-------------------------------------  GREEN'S FUNCTIONS
C
      COMMON/BCPF/ ZBCPF(MNLSQ,MNLSQ,MNB,MNSB)
      COMMON/VCPF/ ZVCPF(MNLSQ,MNLSQ,MNB,MNSV)
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMB(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMV(MNBLSQ,MNBLSQ,MNBZ,MNSB)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMB(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMV(MNBLSQ,MNBLSQ,MNBZ,MNSV)
      COMMON/BOMG/ ZBOMG(MNLSQ,MNLSQ,MNB,MNSB)
      COMMON/VOMG/ ZVOMG(MNLSQ,MNLSQ,MNB,MNSV)
      COMMON/BPF/ ZBPF(MNLSQ,MNAB,MNSB)
      COMMON/VPF/ ZVPF(MNLSQ,MNAV,MNSV)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/GDOS/ ZGDOS(MNL,MNA,MNE,MNS)
      COMMON/GBLO/ ZGBLO(MNL,MNG,MNBZA,MNE,MNS)
      COMMON/BSURA/ ZBSURA(MNBLSQ,MNBLSQ,MNBZA,MNSB)
      COMMON/VSURA/ ZVSURA(MNBLSQ,MNBLSQ,MNBZA,MNSV)
C
C-------------------------------------  LSDA
C
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
C
C-------------------------------------  GENERAL
C
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/DATSUB/ DBA,EFB,EFV,BWST,VWST
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
C
C-------------------------------------  SPECIAL
C
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/YHAR/ YPS(MHARM)
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND,IRANDB,IRANDV
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
C
C******************************************************
C
C                             INITIALIZATION
      CALL PINI
C                             READING OF ALL DATA
      CALL RALL
C                             STARTING CALCULATIONS
      CALL BEGI
C                   REAL SPACE TB-STRUCTURE CONSTANTS
      CALL TBRI
      CALL TBRB
      CALL TBRV
C                             GENERATING K||-MESH
      CALL GIBZ
C                             LFT OF TB-CONSTANTS
      CALL TBK
C                             LFT OF TB-CONSTANTS FOR BSF
      IF (NBZA.GT.0) CALL ATBK
C
C----------------------------------  ENERGY LOOP
      DO 310 IE=1,NE
C                               POTENTIAL FUNCTIONS
      CALL SETPF(IE)
C                            GAMMA OF THE BULK PART
      DO 312 ISB=1,NSB
      CALL GABU(IE,ISB)
312   CONTINUE
C                            GAMMA OF THE VACUUM PART
      DO 314 ISV=1,NSV
      CALL GAVA(IE,ISV)
314   CONTINUE
C                               CPA SELFCONSISTENCY
C                               AND GF FOR DOS
      DO 320 IS=1,NS
      CALL CPAIT(IE,IS)
320   CONTINUE
C                                  GF FOR BSF
      IF (NBZA.GT.0) CALL GBSF(IE)
C--------------------------------  END OF ENERGY LOOP
310   CONTINUE
C                                DOS EVALUATION
      CALL XDOS
C                                BSF EVALUATION
      IF (NBZA.GT.0) CALL XBSF
C
      WRITE(IW6,191)
191   FORMAT(/8X,'* * *   END OF IDS   * * * '/)
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
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
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
      IW17=17
      IW18=18
      IW19=19
C
      OPEN(UNIT=IR1,FILE='inpge',FORM='FORMATTED')
      OPEN(UNIT=IR2,FILE='inpch',FORM='FORMATTED')
      OPEN(UNIT=IR5,FILE='inpds',FORM='FORMATTED')
      OPEN(UNIT=IR11,FILE='inpbu',FORM='FORMATTED')
      OPEN(UNIT=IW6,FILE='otdit',FORM='FORMATTED')
      OPEN(UNIT=IW7,FILE='otdor',FORM='FORMATTED')
      OPEN(UNIT=IW8,FILE='otda1',FORM='FORMATTED')
      OPEN(UNIT=IW9,FILE='otda2',FORM='FORMATTED')
      OPEN(UNIT=IW17,FILE='otbor',FORM='FORMATTED')
      OPEN(UNIT=IW18,FILE='otba1',FORM='FORMATTED')
      OPEN(UNIT=IW19,FILE='otba2',FORM='FORMATTED')
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZ=128)
      PARAMETER(MNBZA=4)
      PARAMETER(MNAB=2)
      PARAMETER(MNAV=2)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      integer,parameter::atomf=999            
C
      DIMENSION OCASE(0:1),OWORK(5),SCX(3),SCKV(2)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/KMESHA/ AKBZA(2,MNBZA),NBZA
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
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
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/DATSUB/ DBA,EFB,EFV,BWST,VWST
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/CUT/ CUTRAT,NMTR
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND,IRANDB,IRANDV
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
      character(len=16)::cwork
      DIMENSION OTXTA1(MNA)      
C
      DATA OCASE/'    SURFACE     ','   INTERFACE    '/
C
100   FORMAT(5A16)
101   FORMAT(1X,10I5)
104   FORMAT(1X,4G15.7)
C
      WRITE(IW6,111)
111   FORMAT(/'       ****  OUTPUT OF IDS  **** '/
     &     /'     ** READING OF ALL DATA:')
C
C-------------------------  CONTROL DATA - UNIT IR5
C
1001  FORMAT(//5A16)
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
      READ (IR5,101) NSYM,NK,INVE
      READ (IR5,101) NE,NDC
      IF(NE.GT.MNE) GO TO 217
      READ (IR5,107) EBOT,ETOP
      READ (IR5,107) DIFMC,QMIXC
      READ (IR5,107) DIFMS,QMIXS
      READ (IR5,101) NBZA
      IF(NBZA.GT.MNBZA) GO TO 218
      IF(NBZA.GT.0) THEN
      READ (IR5,104) (SCKV(I),I=1,2)
      DO 307 IBZA=1,NBZA
      READ (IR5,1041) (AKBZA(I,IBZA),I=1,2),OTXTK(IBZA)
307   CONTINUE
1041  FORMAT(1X,2G15.7,A16)
      END IF
C                                   SCALING
      IF(NBZA.GT.0) THEN
      DO 3081 IBZA=1,NBZA
      DO 3082 I=1,2
      AKBZA(I,IBZA)=SCKV(I)*AKBZA(I,IBZA)
3082  CONTINUE
3081  CONTINUE
      END IF
C
      WRITE(IW6,115) OCASE(IVAC),NL,NS
115   FORMAT(/4X,'CASE:',A16,7X,'NL=',I2,7X,'NS=',I2)
      WRITE(IW6,117) NSYM,NK,INVE
117   FORMAT(/'       LOCAL SYMMETRY:    NSYM=',I2/
     &      /'   K||-MESH:   NK=',I3,10X,'INVE=',I2)
      WRITE(IW6,118) NE,NDC,EBOT,ETOP
118   FORMAT(/'  ENERGY MESH:     NE=',I4,'       NDC=',I2/
     &       /15X,'EBOT=',F12.5,'       ETOP=',F12.5)
      WRITE(IW6,121) DIFMC,QMIXC
121   FORMAT(/'  CPA ITERATIONS:  DIFMAX=',G12.4,
     &        8X,'MIXING=',F10.5)
      WRITE(IW6,122) DIFMS,QMIXS
122   FORMAT(/'  SGF ITERATIONS:  DIFMAX=',G12.4,
     &        8X,'MIXING=',F10.5)
      WRITE(IW6,124) NBZA
124   FORMAT(/'     BLOCH SPECTRAL FUNCTIONS:    NBZA=',I3)
      IF(NBZA.GT.0) THEN
      WRITE(IW6,126) (SCKV(I),I=1,2)
126   FORMAT(/'  SCALING FACTORS: ',2G15.7)
      WRITE(IW6,128)
128   FORMAT(/'        SELECTED K||-VECTORS:')
      DO 309 IBZA=1,NBZA
      WRITE(IW6,1042) IBZA,OTXTK(IBZA),(AKBZA(I,IBZA),I=1,2)
309   CONTINUE
1042  FORMAT(2X,I2,2X,A16,2X,2G15.7)
      END IF
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
1044  FORMAT(3X,'IB=',I2,5X,3G15.7)
      WRITE(IW6,135)
135   FORMAT(/'       SITES OF VACUUM PRINCIPAL LAYER:')
      DO 3101 IB=1,NB
      WRITE(IW6,1044) IB,(VPOS(I,IB),I=1,3)
3101  CONTINUE
      DO 3111 IP=1,NP
      WRITE(IW6,136) IP
136   FORMAT(/'       SITES OF PRINCIPAL LAYER NO.: ',I4)
      DO 3121 IB=1,NB
      WRITE(IW6,1044) IB,(POS(I,IB,IP),I=1,3)
3121  CONTINUE
3111  CONTINUE
      WRITE(IW6,137)
137   FORMAT(/'       SITES OF BULK PRINCIPAL LAYER: ')
      DO 3131 IB=1,NB
      WRITE(IW6,1044) IB,(BPOS(I,IB),I=1,3)
3131  CONTINUE
C
C------------------------- CHEMICAL OCCUPATION - UNIT IR2
C
      READ (IR2,100) OWORK
      WRITE(IW6,1001) OWORK
C
108   FORMAT(A16)
      IG=0
      IA=0
      DO 340 IP=1,NP
      DO 341 IB=1,NB
      IG=IG+1
      READ (IR2,101) NC(IG)
      DO 342 IC=1,NC(IG)
      IA=IA+1
      IF(IA.GT.MNA) GO TO 241
      READ (IR2,108) OTXTA(IA)
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
      WRITE(IW6,142) IA,OTXTA(IA),CON(IA),VALZ(IA)
142   FORMAT(' IA=',I4,' - LABEL,CONC.,VALENCY: ',
     &                       A16,F10.5,F10.3)
345   CONTINUE
344   CONTINUE
343   CONTINUE
      WRITE(IW6,143) NA
143   FORMAT(/4X,' TOTAL NUMBER OF DIFFERENT ',
     &    'ATOMIC TYPES:  NA=',I4)
C
C ---------------------------- LSDA FILE - UNIT IR3
C


      do IA=1,NA
         read(OTXTA(IA),*) cwork
         WRITE(IW6,*)' Read defs for:', CWORK
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
            READ (atomf,101) NRAD1
            READ (atomf,104) (skipPOT,I=1,NRAD1)
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
         READ (atomf,101) NUMCOR
         IF(NUMCOR.NE.0) then
            DO J=1,NUMCOR
               READ (atomf,*) cwork
               READ (atomf,*) cwork
            end do
         end if
         IF(IVAC.EQ.0) READ (atomf,104) DBA                  
         close(atomf)
      end do      
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
C------------------------------------------------------
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
217   WRITE(IW6,917)
917   FORMAT(/' **** INPUT ERROR:  NE GREATER THAN MNE ')
      STOP
218   WRITE(IW6,918)
918   FORMAT(/' **** INPUT ERROR: NBZA GREATER THAN MNBZA ')
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZ=128)
      PARAMETER(MNAB=2)
      PARAMETER(MNAV=2)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
C
      DIMENSION QSC3(3),QSC4(4)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/GEPO/ POS(3,MNB,MNP),AWS(MNB,MNP)
      COMMON/BGEPO/ BPOS(3,MNB),BTRV(3),ABWS
      COMMON/VGEPO/ VPOS(3,MNB),VTRV(3),AVWS
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/BPOPA/ BENY(MNL,MNSB,MNAB),BPPC(MNL,MNSB,MNAB),
     &              BPPD(MNL,MNSB,MNAB),BPPQ(MNL,MNSB,MNAB),
     &              BPPP(MNL,MNSB,MNAB),BDNY(MNL,MNSB,MNAB),
     &              BFINY(MNL,MNSB,MNAB),BFINYD(MNL,MNSB,MNAB)
      COMMON/VPOPA/ VENY(MNL,MNSV,MNAV),VPPC(MNL,MNSV,MNAV),
     &              VPPD(MNL,MNSV,MNAV),VPPQ(MNL,MNSV,MNAV),
     &              VPPP(MNL,MNSV,MNAV),VDNY(MNL,MNSV,MNAV),
     &              VFINY(MNL,MNSV,MNAV),VFINYD(MNL,MNSV,MNAV)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/DATSUB/ DBA,EFB,EFV,BWST,VWST
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/SCREEN/ QSCR(MNL)
      COMMON/NHAR/ CNH(0:MLMAX,0:MLMAX),NCALLH
      COMMON/GHAR/ GFRH(MHARM,MNLSQ,MNLSQ)
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND,IRANDB,IRANDV
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC2/2.0D0/,RC3/3.0D0/
      DATA RCH/0.5D0/,RC4/4.0D0/,TOL/1.0D-6/
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
C---------------------------------- MAX. NUMBER OF SGF
C                                   AND CPA ITERATIONS
      NMITS=500
      WRITE(IW6,117) NMITS
117   FORMAT(/3X,' MAX. NUMBER OF SGF ITERATIONS: ',I5)
      NMITC=50
      WRITE(IW6,118) NMITC
118   FORMAT(/3X,' MAX. NUMBER OF CPA ITERATIONS: ',I5)
C
C---------------------------------- NCALLH FOR SPH. HARMONICS
      NCALLH=0
C
C---------------------------------- GAUNT FACTORS
      CALL GAUNTY(NL)
C
C---------------------------------- SYSTEM FERMI ENERGY
      IF(IVAC.EQ.0) EFT=EFB
      IF(IVAC.EQ.1) EFT=RCH*(EFB+EFV)
      WRITE(IW6,120) EFT
120   FORMAT(/5X,'****  SYSTEM FERMI ENERGY= ',G15.7)
C
C---------------------------------- COMPLEX ENERGIES
C
      EH=(ETOP-EBOT)/REAL(NE-1)
      EI= REAL(NDC)*EH
      DO 310 IE=1,NE
      ER=EBOT+REAL(IE-1)*EH
      ZCN(IE)=DCMPLX(ER,EI)
310   CONTINUE
      WRITE(IW6,121) EH,EI
121   FORMAT(/9X,'ENERGY MESH:   STEP= ',G15.7/
     &        9X,'     IMAGINARY PART= ',G15.7)
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
      VBG(1,1)=CNORM*VBR(2,2)
      VBG(2,1)=-CNORM*VBR(1,2)
      VBG(1,2)=-CNORM*VBR(2,1)
      VBG(2,2)=CNORM*VBR(1,1)
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
      IF(ABS(DIF).GT.TOL) GO TO 292
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
C---------------------------   RANDOM/NON-RANDOM CASE
      IRAND=0
      DO 370 IA=1,NA
      PROD=CON(IA)*(RC1-CON(IA))
      IF(PROD.GT.RCZ) IRAND=1
370   CONTINUE
      IF(IRAND.EQ.0) WRITE(IW6,170)
170   FORMAT(/5X,'--- INTERMEDIATE REGION ',
     &           'WITHOUT DISORDER ---')
C
      IRANDB=0
      DO 371 IA=1,NAB
      PROD=BCON(IA)*(RC1-BCON(IA))
      IF(PROD.GT.RCZ) IRANDB=1
371   CONTINUE
      IF(IRANDB.EQ.0) WRITE(IW6,171)
171   FORMAT(/5X,'--- BULK SUBSTRATE ',
     &           'WITHOUT DISORDER ---')
C
      IRANDV=0
      DO 372 IA=1,NAV
      PROD=VCON(IA)*(RC1-VCON(IA))
      IF(PROD.GT.RCZ) IRANDV=1
372   CONTINUE
      IF(IRANDV.EQ.0) WRITE(IW6,172)
172   FORMAT(/5X,'--- VAC. SUBSTRATE ',
     &           'WITHOUT DISORDER ---')
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
      PARAMETER(MNP=24)
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNPP2=MNP+2)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=32)
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=32)
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MLMAX=2*(MNL-1))
      PARAMETER(MHARM=(MLMAX+1)**2)
      PARAMETER(MNCL=32)
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
      PARAMETER(MNBZ=128)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
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
      PARAMETER(MNP=24)
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNBZ=128)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=32)
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
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKI(JJ,II,IP,IBZ)
      W2X=DREAL(ZW2)
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
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKO(JJ,II,IPP1,IBZ)
      W2X=DREAL(ZW2)
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
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZBSKI(JJ,II,IBZ)
      W2X=DREAL(ZW2)
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
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZBSKO(JJ,II,IBZ)
      W2X=DREAL(ZW2)
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
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZVSKI(JJ,II,IBZ)
      W2X=DREAL(ZW2)
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
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZVSKO(JJ,II,IBZ)
      W2X=DREAL(ZW2)
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
      XW=DREAL(ZA(I,J))
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNAB=2)
      PARAMETER(MNAV=2)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
C
      COMMON/BPF/ ZBPF(MNLSQ,MNAB,MNSB)
      COMMON/VPF/ ZVPF(MNLSQ,MNAV,MNSV)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
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
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/SCREEN/ QSCR(MNL)
C
      DATA RC1/1.0D0/
C
      ZE1=ZCN(IE)
C
      DO 341 IS=1,NS
      DO 342 IA=1,NA
      DO 343 IL=1,NL
      PPC1=PPC(IL,IS,IA)
      PPD1=PPD(IL,IS,IA)
      PPQ1=PPQ(IL,IS,IA)
      ALF1=QSCR(IL)
      CALL PLMZ(PPC1,PPD1,PPQ1,ALF1,ZE1,ZPF1,ZLA1,ZMU1)
      ARG = RC1/FINY(IL,IS,IA)
     &    - WS(IA)*FINYD(IL,IS,IA)*(DNY(IL,IS,IA)+REAL(IL))
      YNAM1=SIGN(RC1,ARG)
      ZMU1=YNAM1*ZMU1
C
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 346 I=ISTA,IFIN
      ZPF(I,IA,IS)=ZPF1
      ZMU(I,IA,IS)=ZMU1
      ZLA(I,IA,IS)=ZLA1
346   CONTINUE
343   CONTINUE
342   CONTINUE
341   CONTINUE
C
      DO 351 IS=1,NSB
      DO 352 IA=1,NAB
      DO 353 IL=1,NL
      PPC1=BPPC(IL,IS,IA)
      PPD1=BPPD(IL,IS,IA)
      PPQ1=BPPQ(IL,IS,IA)
      ALF1=QSCR(IL)
      CALL PLMZ(PPC1,PPD1,PPQ1,ALF1,ZE1,ZPF1,ZLA1,ZMU1)
C
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 356 I=ISTA,IFIN
      ZBPF(I,IA,IS)=ZPF1
356   CONTINUE
353   CONTINUE
352   CONTINUE
351   CONTINUE
C
      DO 361 IS=1,NSV
      DO 362 IA=1,NAV
      DO 363 IL=1,NL
      PPC1=VPPC(IL,IS,IA)
      PPD1=VPPD(IL,IS,IA)
      PPQ1=VPPQ(IL,IS,IA)
      ALF1=QSCR(IL)
      CALL PLMZ(PPC1,PPD1,PPQ1,ALF1,ZE1,ZPF1,ZLA1,ZMU1)
C
      ISTA=(IL-1)**2+1
      IFIN=IL**2
      DO 366 I=ISTA,IFIN
      ZVPF(I,IA,IS)=ZPF1
366   CONTINUE
363   CONTINUE
362   CONTINUE
361   CONTINUE
C
      RETURN
      END
C*******************
CXXX    GABU    ****
C*******************
      SUBROUTINE GABU(IE,IS)
C
C***************************************************
C    GAMMA FOR THE SEMIINFINITE BULK RANDOM ALLOY -
C     -  FOR EACH K||-POINT
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNBZ=128)
      PARAMETER(MNAB=2)
      PARAMETER(MNSB=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZGF(MNLSQ,MNLSQ,MNB),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      DIMENSION ZPSI1(MNBLSQ,MNBLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZGMB1(MNBLSQ,MNBLSQ),ZGMV1(MNBLSQ,MNBLSQ),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/BSCK/ ZBSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZBSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/BCPF/ ZBCPF(MNLSQ,MNLSQ,MNB,MNSB)
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMB(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMV(MNBLSQ,MNBLSQ,MNBZ,MNSB)
      COMMON/BOMG/ ZBOMG(MNLSQ,MNLSQ,MNB,MNSB)
      COMMON/BPF/ ZBPF(MNLSQ,MNAB,MNSB)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHB/ BCON(MNAB),BAZ(MNAB),BWS(MNAB),BWSAV(MNAB)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND,IRANDB,IRANDV
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      ITS=0
C
      IF(IE.EQ.1) THEN
      DO 3001 IBZ=1,MNBZ
      DO 3002 JJ=1,MNBLSQ
      DO 3003 II=1,MNBLSQ
      ZBGMB(II,JJ,IBZ,IS)=ZZ
3003  CONTINUE
3002  CONTINUE
3001  CONTINUE
      DO 3004 IBZ=1,MNBZ
      DO 3005 JJ=1,MNBLSQ
      DO 3006 II=1,MNBLSQ
      ZBGMV(II,JJ,IBZ,IS)=ZZ
3006  CONTINUE
3005  CONTINUE
3004  CONTINUE
      DO 3007 IB=1,MNB
      DO 3008 J=1,MNLSQ
      DO 3009 I=1,MNLSQ
      ZBOMG(I,J,IB,IS)=ZZ
3009  CONTINUE
3008  CONTINUE
3007  CONTINUE
      END IF
C
210   ITS=ITS+1
      IF(ITS.GT.NMITS) GO TO 291
C
C------------------------------------ CPA-PART
C
      NONC=IRANDB*NB
      IF(ITS.GT.1) GO TO 212
C..................................  1ST ITERATION
C
      DO 3020 IB=1,NB
C
      DO 3021 J=1,NLSQ
      DO 3022 I=1,NLSQ
      ZOM1(I,J)=ZBOMG(I,J,IB,IS)
3022  CONTINUE
3021  CONTINUE
      DO 3024 J=1,NLSQ
      DO 3025 I=1,NLSQ
      ZFI1(I,J)=ZZ
3025  CONTINUE
3024  CONTINUE
C
      DO 3027 IA=1,NAB
      DO 3030 J=1,NLSQ
      DO 3031 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3031  CONTINUE
3030  CONTINUE
      DO 3032 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZBPF(I,IA,IS)
3032  CONTINUE
      CALL CINVC(ZW,NLSQ,MNLSQ)
      DO 3034 J=1,NLSQ
      DO 3035 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+BCON(IA)*ZW(I,J)
3035  CONTINUE
3034  CONTINUE
3027  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3037 J=1,NLSQ
      DO 3038 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3038  CONTINUE
3037  CONTINUE
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3041 J=1,NLSQ
      DO 3042 I=1,NLSQ
      ZBCPF(I,J,IB,IS)=ZCP1(I,J)
3042  CONTINUE
3041  CONTINUE
C
3020  CONTINUE
C
      GO TO 220
C................................... NEXT ITERATIONS
C
212   IF (IRANDB.EQ.1) THEN
      DO 3100 IB=1,NB
C
      DO 3110 J=1,NLSQ
      DO 3111 I=1,NLSQ
      ZFI1(I,J)=ZGF(I,J,IB)
3111  CONTINUE
3110  CONTINUE
      DO 3113 J=1,NLSQ
      DO 3114 I=1,NLSQ
      ZCP1(I,J)=ZBCPF(I,J,IB,IS)
3114  CONTINUE
3113  CONTINUE
      DO 3116 J=1,NLSQ
      DO 3117 I=1,NLSQ
      ZOM1(I,J)=ZBOMG(I,J,IB,IS)
3117  CONTINUE
3116  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3120 J=1,NLSQ
      DO 3121 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
3121  CONTINUE
3120  CONTINUE
C
      SUM=CMANO(ZW,NLSQ,MNLSQ)
      SUM1=CMANO(ZOM1,NLSQ,MNLSQ)
      IF(SUM.LT.(DIFMC*SUM1)) NONC=NONC-1
      DO 3123 J=1,NLSQ
      DO 3124 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXC*ZW(I,J)
3124  CONTINUE
3123  CONTINUE
      CALL SYMGF(ZOM1,NSYM,NL)
C
      DO 3126 J=1,NLSQ
      DO 3127 I=1,NLSQ
      ZFI1(I,J)=ZZ
3127  CONTINUE
3126  CONTINUE
C
      DO 3130 IA=1,NAB
      DO 3132 J=1,NLSQ
      DO 3133 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3133  CONTINUE
3132  CONTINUE
      DO 3134 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZBPF(I,IA,IS)
3134  CONTINUE
      CALL CINVC(ZW,NLSQ,MNLSQ)
      DO 3136 J=1,NLSQ
      DO 3137 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+BCON(IA)*ZW(I,J)
3137  CONTINUE
3136  CONTINUE
3130  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3141 J=1,NLSQ
      DO 3142 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3142  CONTINUE
3141  CONTINUE
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3144 J=1,NLSQ
      DO 3145 I=1,NLSQ
      ZBCPF(I,J,IB,IS)=ZCP1(I,J)
3145  CONTINUE
3144  CONTINUE
      DO 3147 J=1,NLSQ
      DO 3148 I=1,NLSQ
      ZBOMG(I,J,IB,IS)=ZOM1(I,J)
3148  CONTINUE
3147  CONTINUE
C
3100  CONTINUE
      END IF
C
C---------------------------------------- SGF-PART
C
220   NONS=(1+IRANDB)*NBZ
C
      DO 3901 IB=1,MNB
      DO 3902 J=1,MNLSQ
      DO 3903 I=1,MNLSQ
      ZGF(I,J,IB)=ZZ
3903  CONTINUE
3902  CONTINUE
3901  CONTINUE
C...................................... IBZ-LOOP
      DO 320 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3228 JJ=1,NBLSQ
      DO 3229 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZBSKI(II,JJ,IBZ)
3229  CONTINUE
3228  CONTINUE
      DO 324 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3258 J=1,NLSQ
      DO 3259 I=1,NLSQ
      ZPSI1(II0+I,II0+J) = ZPSI1(II0+I,II0+J)
     &                   + ZBCPF(I,J,IB,IS)
3259  CONTINUE
3258  CONTINUE
324   CONTINUE
C                              TRANSFER MATRICES
      DO 3271 JJ=1,NBLSQ
      DO 3272 II=1,NBLSQ
      ZSKF1(II,JJ)=ZBSKF(II,JJ,IBZ)
3272  CONTINUE
3271  CONTINUE
      DO 3281 JJ=1,NBLSQ
      DO 3282 II=1,NBLSQ
      ZSKO1(II,JJ)=ZBSKO(II,JJ,IBZ)
3282  CONTINUE
3281  CONTINUE
C
C                              ITERATION OF BULK-SIDE GAMMA
      DO 3318 JJ=1,NBLSQ
      DO 3319 II=1,NBLSQ
      ZGMB1(II,JJ)=ZBGMB(II,JJ,IBZ,IS)
3319  CONTINUE
3318  CONTINUE
C
      DO 3328 JJ=1,NBLSQ
      DO 3329 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMB1(II,JJ)
3329  CONTINUE
3328  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3361 JJ=1,NBLSQ
      DO 3362 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMB1(II,JJ)
3362  CONTINUE
3361  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMB1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3381 JJ=1,NBLSQ
      DO 3382 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB1(II,JJ)+QMIXS*ZWA(II,JJ)
3382  CONTINUE
3381  CONTINUE
      DO 3391 JJ=1,NBLSQ
      DO 3392 II=1,NBLSQ
      ZBGMB(II,JJ,IBZ,IS)=ZGMB1(II,JJ)
3392  CONTINUE
3391  CONTINUE
C
C                           ITERATION OF VACUUM-SIDE GAMMA
      IF (IRANDB.EQ.1) THEN
      DO 3418 JJ=1,NBLSQ
      DO 3419 II=1,NBLSQ
      ZGMV1(II,JJ)=ZBGMV(II,JJ,IBZ,IS)
3419  CONTINUE
3418  CONTINUE
C
      DO 3428 JJ=1,NBLSQ
      DO 3429 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)
3429  CONTINUE
3428  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3461 JJ=1,NBLSQ
      DO 3462 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMV1(II,JJ)
3462  CONTINUE
3461  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMV1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3481 JJ=1,NBLSQ
      DO 3482 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV1(II,JJ)+QMIXS*ZWA(II,JJ)
3482  CONTINUE
3481  CONTINUE
      DO 3491 JJ=1,NBLSQ
      DO 3492 II=1,NBLSQ
      ZBGMV(II,JJ,IBZ,IS)=ZGMV1(II,JJ)
3492  CONTINUE
3491  CONTINUE
      END IF
C
C                                  ON-SITE GF-BLOCKS
      IF (IRANDB.EQ.1) THEN
      DO 3518 JJ=1,NBLSQ
      DO 3519 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)-ZGMB1(II,JJ)
3519  CONTINUE
3518  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 353 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3551 J=1,NLSQ
      DO 3552 I=1,NLSQ
      ZGF(I,J,IB) = ZGF(I,J,IB) + WKBZ(IBZ) * ZWA(II0+I,II0+J)
3552  CONTINUE
3551  CONTINUE
353   CONTINUE
      END IF
C............................... END OF IBZ-LOOP
320   CONTINUE
C
      DO 357 IB=1,NB
      DO 3581 J=1,NLSQ
      DO 3582 I=1,NLSQ
      ZW(I,J)=ZGF(I,J,IB)
3582  CONTINUE
3581  CONTINUE
      CALL SYMGF(ZW,NSYM,NL)
      DO 3591 J=1,NLSQ
      DO 3592 I=1,NLSQ
      ZGF(I,J,IB)=ZW(I,J)
3592  CONTINUE
3591  CONTINUE
357   CONTINUE
C
      IF((NONS+NONC).GT.0) GO TO 210
C
      WRITE(IW6,150) IE,IS,ITS
150   FORMAT(9X,'IE=',I4,'    IS=',I1,
     &   '    GAMMA-BULK : ',I4,' ITER.')
C
C----------------------------- CONNECTION TO THE
C                              INTERMEDIATE REGION
      DO 360 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3618 JJ=1,NBLSQ
      DO 3619 II=1,NBLSQ
      ZGMB1(II,JJ)=ZBGMB(II,JJ,IBZ,IS)
3619  CONTINUE
3618  CONTINUE
      DO 3628 JJ=1,NBLSQ
      DO 3629 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZBSKI(II,JJ,IBZ)
3629  CONTINUE
3628  CONTINUE
      DO 364 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3651 J=1,NLSQ
      DO 3652 I=1,NLSQ
      ZPSI1(II0+I,II0+J) = ZPSI1(II0+I,II0+J)
     &                   + ZBCPF(I,J,IB,IS)
3652  CONTINUE
3651  CONTINUE
364   CONTINUE
C                              TRANSFER MATRICES
      DO 3671 JJ=1,NBLSQ
      DO 3672 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKF(II,JJ,NP,IBZ)
3672  CONTINUE
3671  CONTINUE
      DO 3681 JJ=1,NBLSQ
      DO 3682 II=1,NBLSQ
      ZSKO1(II,JJ)=DCONJG(ZSKF1(JJ,II))
3682  CONTINUE
3681  CONTINUE
C                                     INVERSION AND
C                                    MULTIPLICATIONS
      DO 3718 JJ=1,NBLSQ
      DO 3719 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMB1(II,JJ)
3719  CONTINUE
3718  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3738 JJ=1,NBLSQ
      DO 3739 II=1,NBLSQ
      ZBGAM(II,JJ,IBZ,IS)=ZWA(II,JJ)
3739  CONTINUE
3738  CONTINUE
C
360   CONTINUE
C
      RETURN
C
C------------------------------ NO CONVERGENCY
C
291   WRITE(IW6,191) NMITS,IE,IS
191   FORMAT(//' ****  GABU:  NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   '      ---  FOR IE=',I4,'    IS=',I2)
      STOP
      END
C*******************
CXXX    GAVA    ****
C*******************
      SUBROUTINE GAVA(IE,IS)
C
C***************************************************
C    GAMMA FOR THE SEMIINFINITE BULK RANDOM ALLOY -
C     -  FOR EACH K||-POINT
C***************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNBZ=128)
      PARAMETER(MNAV=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZGF(MNLSQ,MNLSQ,MNB),
     &          ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      DIMENSION ZPSI1(MNBLSQ,MNBLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZGMB1(MNBLSQ,MNBLSQ),ZGMV1(MNBLSQ,MNBLSQ),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/VSCK/ ZVSKI(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKO(MNBLSQ,MNBLSQ,MNBZ),
     &             ZVSKF(MNBLSQ,MNBLSQ,MNBZ)
      COMMON/VCPF/ ZVCPF(MNLSQ,MNLSQ,MNB,MNSV)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMB(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMV(MNBLSQ,MNBLSQ,MNBZ,MNSV)
      COMMON/VOMG/ ZVOMG(MNLSQ,MNLSQ,MNB,MNSV)
      COMMON/VPF/ ZVPF(MNLSQ,MNAV,MNSV)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHV/ VCON(MNAV),VAZ(MNAV),VWS(MNAV),VWSAV(MNAV)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND,IRANDB,IRANDV
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      ZZ=DCMPLX(RCZ,RCZ)
C
      ITS=0
C
      IF(IE.EQ.1) THEN
      DO 3001 IBZ=1,MNBZ
      DO 3002 JJ=1,MNBLSQ
      DO 3003 II=1,MNBLSQ
      ZVGMB(II,JJ,IBZ,IS)=ZZ
3003  CONTINUE
3002  CONTINUE
3001  CONTINUE
      DO 3004 IBZ=1,MNBZ
      DO 3005 JJ=1,MNBLSQ
      DO 3006 II=1,MNBLSQ
      ZVGMV(II,JJ,IBZ,IS)=ZZ
3006  CONTINUE
3005  CONTINUE
3004  CONTINUE
      DO 3007 IB=1,MNB
      DO 3008 J=1,MNLSQ
      DO 3009 I=1,MNLSQ
      ZVOMG(I,J,IB,IS)=ZZ
3009  CONTINUE
3008  CONTINUE
3007  CONTINUE
      END IF
C
210   ITS=ITS+1
      IF(ITS.GT.NMITS) GO TO 291
C
C------------------------------------ CPA-PART
C
      NONC=IRANDV*NB
      IF(ITS.GT.1) GO TO 212
C..................................  1ST ITERATION
C
      DO 3020 IB=1,NB
C
      DO 3021 J=1,NLSQ
      DO 3022 I=1,NLSQ
      ZOM1(I,J)=ZVOMG(I,J,IB,IS)
3022  CONTINUE
3021  CONTINUE
      DO 3024 J=1,NLSQ
      DO 3025 I=1,NLSQ
      ZFI1(I,J)=ZZ
3025  CONTINUE
3024  CONTINUE
C
      DO 3027 IA=1,NAV
      DO 3030 J=1,NLSQ
      DO 3031 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3031  CONTINUE
3030  CONTINUE
      DO 3032 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZVPF(I,IA,IS)
3032  CONTINUE
      CALL CINVC(ZW,NLSQ,MNLSQ)
      DO 3034 J=1,NLSQ
      DO 3035 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+VCON(IA)*ZW(I,J)
3035  CONTINUE
3034  CONTINUE
3027  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3037 J=1,NLSQ
      DO 3038 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3038  CONTINUE
3037  CONTINUE
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3041 J=1,NLSQ
      DO 3042 I=1,NLSQ
      ZVCPF(I,J,IB,IS)=ZCP1(I,J)
3042  CONTINUE
3041  CONTINUE
C
3020  CONTINUE
C
      GO TO 220
C................................... NEXT ITERATIONS
C
212   IF (IRANDV.EQ.1) THEN
      DO 3100 IB=1,NB
C
      DO 3110 J=1,NLSQ
      DO 3111 I=1,NLSQ
      ZFI1(I,J)=ZGF(I,J,IB)
3111  CONTINUE
3110  CONTINUE
      DO 3113 J=1,NLSQ
      DO 3114 I=1,NLSQ
      ZCP1(I,J)=ZVCPF(I,J,IB,IS)
3114  CONTINUE
3113  CONTINUE
      DO 3116 J=1,NLSQ
      DO 3117 I=1,NLSQ
      ZOM1(I,J)=ZVOMG(I,J,IB,IS)
3117  CONTINUE
3116  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3120 J=1,NLSQ
      DO 3121 I=1,NLSQ
      ZW(I,J)=ZCP1(I,J)-ZFI1(I,J)-ZOM1(I,J)
3121  CONTINUE
3120  CONTINUE
C
      SUM=CMANO(ZW,NLSQ,MNLSQ)
      SUM1=CMANO(ZOM1,NLSQ,MNLSQ)
      IF(SUM.LT.(DIFMC*SUM1)) NONC=NONC-1
      DO 3123 J=1,NLSQ
      DO 3124 I=1,NLSQ
      ZOM1(I,J)=ZOM1(I,J)+QMIXC*ZW(I,J)
3124  CONTINUE
3123  CONTINUE
      CALL SYMGF(ZOM1,NSYM,NL)
C
      DO 3126 J=1,NLSQ
      DO 3127 I=1,NLSQ
      ZFI1(I,J)=ZZ
3127  CONTINUE
3126  CONTINUE
C
      DO 3130 IA=1,NAV
      DO 3132 J=1,NLSQ
      DO 3133 I=1,NLSQ
      ZW(I,J)=-ZOM1(I,J)
3133  CONTINUE
3132  CONTINUE
      DO 3134 I=1,NLSQ
      ZW(I,I)=ZW(I,I)+ZVPF(I,IA,IS)
3134  CONTINUE
      CALL CINVC(ZW,NLSQ,MNLSQ)
      DO 3136 J=1,NLSQ
      DO 3137 I=1,NLSQ
      ZFI1(I,J)=ZFI1(I,J)+VCON(IA)*ZW(I,J)
3137  CONTINUE
3136  CONTINUE
3130  CONTINUE
C
      CALL CINVC(ZFI1,NLSQ,MNLSQ)
      DO 3141 J=1,NLSQ
      DO 3142 I=1,NLSQ
      ZCP1(I,J)=ZFI1(I,J)+ZOM1(I,J)
3142  CONTINUE
3141  CONTINUE
      CALL SYMGF(ZCP1,NSYM,NL)
C
      DO 3144 J=1,NLSQ
      DO 3145 I=1,NLSQ
      ZVCPF(I,J,IB,IS)=ZCP1(I,J)
3145  CONTINUE
3144  CONTINUE
      DO 3147 J=1,NLSQ
      DO 3148 I=1,NLSQ
      ZVOMG(I,J,IB,IS)=ZOM1(I,J)
3148  CONTINUE
3147  CONTINUE
C
3100  CONTINUE
      END IF
C
C---------------------------------------- SGF-PART
C
220   NONS=(1+IRANDV)*NBZ
C
      DO 3901 IB=1,MNB
      DO 3902 J=1,MNLSQ
      DO 3903 I=1,MNLSQ
      ZGF(I,J,IB)=ZZ
3903  CONTINUE
3902  CONTINUE
3901  CONTINUE
C...................................... IBZ-LOOP
      DO 320 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3228 JJ=1,NBLSQ
      DO 3229 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZVSKI(II,JJ,IBZ)
3229  CONTINUE
3228  CONTINUE
      DO 324 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3258 J=1,NLSQ
      DO 3259 I=1,NLSQ
      ZPSI1(II0+I,II0+J) = ZPSI1(II0+I,II0+J)
     &                   + ZVCPF(I,J,IB,IS)
3259  CONTINUE
3258  CONTINUE
324   CONTINUE
C                              TRANSFER MATRICES
      DO 3271 JJ=1,NBLSQ
      DO 3272 II=1,NBLSQ
      ZSKF1(II,JJ)=ZVSKF(II,JJ,IBZ)
3272  CONTINUE
3271  CONTINUE
      DO 3281 JJ=1,NBLSQ
      DO 3282 II=1,NBLSQ
      ZSKO1(II,JJ)=ZVSKO(II,JJ,IBZ)
3282  CONTINUE
3281  CONTINUE
C
C                              ITERATION OF BULK-SIDE GAMMA
      IF (IRANDV.EQ.1) THEN
      DO 3318 JJ=1,NBLSQ
      DO 3319 II=1,NBLSQ
      ZGMB1(II,JJ)=ZVGMB(II,JJ,IBZ,IS)
3319  CONTINUE
3318  CONTINUE
C
      DO 3328 JJ=1,NBLSQ
      DO 3329 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMB1(II,JJ)
3329  CONTINUE
3328  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3361 JJ=1,NBLSQ
      DO 3362 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMB1(II,JJ)
3362  CONTINUE
3361  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMB1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3381 JJ=1,NBLSQ
      DO 3382 II=1,NBLSQ
      ZGMB1(II,JJ)=ZGMB1(II,JJ)+QMIXS*ZWA(II,JJ)
3382  CONTINUE
3381  CONTINUE
      DO 3391 JJ=1,NBLSQ
      DO 3392 II=1,NBLSQ
      ZVGMB(II,JJ,IBZ,IS)=ZGMB1(II,JJ)
3392  CONTINUE
3391  CONTINUE
      END IF
C
C                           ITERATION OF VACUUM-SIDE GAMMA
      DO 3418 JJ=1,NBLSQ
      DO 3419 II=1,NBLSQ
      ZGMV1(II,JJ)=ZVGMV(II,JJ,IBZ,IS)
3419  CONTINUE
3418  CONTINUE
C
      DO 3428 JJ=1,NBLSQ
      DO 3429 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)
3429  CONTINUE
3428  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3461 JJ=1,NBLSQ
      DO 3462 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZGMV1(II,JJ)
3462  CONTINUE
3461  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZGMV1,NBLSQ,MNBLSQ)
      IF(SUM.LT.(DIFMS*SUM1)) NONS=NONS-1
      DO 3481 JJ=1,NBLSQ
      DO 3482 II=1,NBLSQ
      ZGMV1(II,JJ)=ZGMV1(II,JJ)+QMIXS*ZWA(II,JJ)
3482  CONTINUE
3481  CONTINUE
      DO 3491 JJ=1,NBLSQ
      DO 3492 II=1,NBLSQ
      ZVGMV(II,JJ,IBZ,IS)=ZGMV1(II,JJ)
3492  CONTINUE
3491  CONTINUE
C
C                                  ON-SITE GF-BLOCKS
      IF (IRANDV.EQ.1) THEN
      DO 3518 JJ=1,NBLSQ
      DO 3519 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)-ZGMB1(II,JJ)
3519  CONTINUE
3518  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 353 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3551 J=1,NLSQ
      DO 3552 I=1,NLSQ
      ZGF(I,J,IB) = ZGF(I,J,IB) + WKBZ(IBZ) * ZWA(II0+I,II0+J)
3552  CONTINUE
3551  CONTINUE
353   CONTINUE
      END IF
C............................... END OF IBZ-LOOP
320   CONTINUE
C
      DO 357 IB=1,NB
      DO 3581 J=1,NLSQ
      DO 3582 I=1,NLSQ
      ZW(I,J)=ZGF(I,J,IB)
3582  CONTINUE
3581  CONTINUE
      CALL SYMGF(ZW,NSYM,NL)
      DO 3591 J=1,NLSQ
      DO 3592 I=1,NLSQ
      ZGF(I,J,IB)=ZW(I,J)
3592  CONTINUE
3591  CONTINUE
357   CONTINUE
C
      IF((NONS+NONC).GT.0) GO TO 210
C
      WRITE(IW6,150) IE,IS,ITS
150   FORMAT(9X,'IE=',I4,'    IS=',I1,
     &   '    GAMMA-VAC. : ',I4,' ITER.')
C
C----------------------------- CONNECTION TO THE
C                              INTERMEDIATE REGION
      DO 360 IBZ=1,NBZ
C                             IN-LAYER QUANTITIES
      DO 3618 JJ=1,NBLSQ
      DO 3619 II=1,NBLSQ
      ZGMV1(II,JJ)=ZVGMV(II,JJ,IBZ,IS)
3619  CONTINUE
3618  CONTINUE
      DO 3628 JJ=1,NBLSQ
      DO 3629 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZVSKI(II,JJ,IBZ)
3629  CONTINUE
3628  CONTINUE
      DO 364 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3651 J=1,NLSQ
      DO 3652 I=1,NLSQ
      ZPSI1(II0+I,II0+J) = ZPSI1(II0+I,II0+J)
     &                   + ZVCPF(I,J,IB,IS)
3652  CONTINUE
3651  CONTINUE
364   CONTINUE
C                              TRANSFER MATRICES
      DO 3671 JJ=1,NBLSQ
      DO 3672 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKO(II,JJ,1,IBZ)
3672  CONTINUE
3671  CONTINUE
      DO 3681 JJ=1,NBLSQ
      DO 3682 II=1,NBLSQ
      ZSKF1(II,JJ)=DCONJG(ZSKO1(JJ,II))
3682  CONTINUE
3681  CONTINUE
C                                     INVERSION AND
C                                    MULTIPLICATIONS
      DO 3718 JJ=1,NBLSQ
      DO 3719 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZGMV1(II,JJ)
3719  CONTINUE
3718  CONTINUE
C
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
C
      DO 3738 JJ=1,NBLSQ
      DO 3739 II=1,NBLSQ
      ZVGAM(II,JJ,IBZ,IS)=ZWA(II,JJ)
3739  CONTINUE
3738  CONTINUE
C
360   CONTINUE
C
      RETURN
C
C------------------------------ NO CONVERGENCY
C
291   WRITE(IW6,191) NMITS,IE,IS
191   FORMAT(//' ****  GAVA:  NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   '      ---  FOR IE=',I4,'    IS=',I2)
      STOP
      END
C*******************
CXXX    REPA    ****
C*******************
      SUBROUTINE REPA(IS)
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNBZ=128)
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
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMB(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMV(MNBLSQ,MNBLSQ,MNBZ,MNSB)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMB(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMV(MNBLSQ,MNBLSQ,MNBZ,MNSV)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
C
      DATA RCZ/0.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      NG=NB*NP
      ZZ=DCMPLX(RCZ,RCZ)
C
      DO 301 IG=1,MNG
      DO 302 J=1,MNLSQ
      DO 303 I=1,MNLSQ
      ZGFG(I,J,IG,IS)=ZZ
303   CONTINUE
302   CONTINUE
301   CONTINUE
C
      ISB= MIN(IS,NSB)
      ISV= MIN(IS,NSV)
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
     &                     + ZCPF(I,J,IG,IS)
3182  CONTINUE
3181  CONTINUE
317   CONTINUE
316   CONTINUE
C
C                              RECURSION FROM VACUUM
      DO 3208 JJ=1,NBLSQ
      DO 3209 II=1,NBLSQ
      ZGMV(II,JJ,1)=ZVGAM(II,JJ,IBZ,ISV)
3209  CONTINUE
3208  CONTINUE
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
      ZGMB(II,JJ,NP)=ZBGAM(II,JJ,IBZ,ISB)
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
      ZGFG(I,J,IG,IS) = ZGFG(I,J,IG,IS)
     &            + WKBZ(IBZ) * ZWA(II0+I,II0+J)
3462  CONTINUE
3461  CONTINUE
344   CONTINUE
340   CONTINUE
C
C---------------------------------- END OF IBZ-LOOP
312   CONTINUE
C
C---------------------------------- SYMMETRIZATION OF THE
C                                   ON-SITE BLOCKS
      DO 350 IG=1,NG
      DO 3528 J=1,NLSQ
      DO 3529 I=1,NLSQ
      ZWG(I,J)=ZGFG(I,J,IG,IS)
3529  CONTINUE
3528  CONTINUE
      CALL SYMGF(ZWG,NSYM,NL)
      DO 3548 J=1,NLSQ
      DO 3549 I=1,NLSQ
      ZGFG(I,J,IG,IS)=ZWG(I,J)
3549  CONTINUE
3548  CONTINUE
350   CONTINUE
C
      RETURN
      END
C*******************
CXXX    CPAIT   ****
C*******************
      SUBROUTINE CPAIT(IE,IS)
C
C*****************************************
C   PERFORMS CPA-SELFCONSISTENCY LOOP AND
C   CALCULATES GREEN FUNCTIONS FOR DOS
C*****************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZ=128)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZFI1(MNLSQ,MNLSQ),ZCP1(MNLSQ,MNLSQ),
     &          ZOM1(MNLSQ,MNLSQ),ZW(MNLSQ,MNLSQ)
C
      COMMON/KMESH/ AKBZ(2,MNBZ),WKBZ(MNBZ),NSYM,NK,INVE,NBZ
      COMMON/SCK/ ZSKI(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKO(MNBLSQ,MNBLSQ,MNP,MNBZ),
     &            ZSKF(MNBLSQ,MNBLSQ,MNP,MNBZ)
      COMMON/BGAM/ ZBGAM(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMB(MNBLSQ,MNBLSQ,MNBZ,MNSB),
     &             ZBGMV(MNBLSQ,MNBLSQ,MNBZ,MNSB)
      COMMON/VGAM/ ZVGAM(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMB(MNBLSQ,MNBLSQ,MNBZ,MNSV),
     &             ZVGMV(MNBLSQ,MNBLSQ,MNBZ,MNSV)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/OMG/ ZOMG(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/GDOS/ ZGDOS(MNL,MNA,MNE,MNS)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/ITCPA/ DIFMC,QMIXC,NMITC,IRAND,IRANDB,IRANDV
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2
      NG=NP*NB
      ZZ=DCMPLX(RCZ,RCZ)
C
      QMIXI=QMIXC
      IF(IRAND.EQ.0) QMIXI=RC1
C
      IITER=0
C                             STARTING OMEGA
      IF(IE.EQ.1) THEN
      DO 3001 IG=1,MNG
      DO 3002 J=1,MNLSQ
      DO 3003 I=1,MNLSQ
      ZOMG(I,J,IG,IS)=ZZ
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
      DO 3100 IG=1,NG
C
      DO 3116 J=1,NLSQ
      DO 3117 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IG,IS)
3117  CONTINUE
3116  CONTINUE
      DO 3124 J=1,NLSQ
      DO 3125 I=1,NLSQ
      ZFI1(I,J)=ZZ
3125  CONTINUE
3124  CONTINUE
C
      DO 3130 IC=1,NC(IG)
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
      ZCPF(I,J,IG,IS)=ZCP1(I,J)
3155  CONTINUE
3154  CONTINUE
C
3100  CONTINUE
C
C::::::::::::::::::::::::::::::::::::   STEP 2
C
      CALL REPA(IS)
C
C::::::::::::::::::::::::::::::::::::   STEP 3
C
      DO 3300 IG=1,NG
C
      DO 3304 J=1,NLSQ
      DO 3305 I=1,NLSQ
      ZFI1(I,J)=ZGFG(I,J,IG,IS)
3305  CONTINUE
3304  CONTINUE
      DO 3308 J=1,NLSQ
      DO 3309 I=1,NLSQ
      ZCP1(I,J)=ZCPF(I,J,IG,IS)
3309  CONTINUE
3308  CONTINUE
      DO 3312 J=1,NLSQ
      DO 3313 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IG,IS)
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
      SUM=CMANO(ZW,NLSQ,MNLSQ)
      SUM1=CMANO(ZOM1,NLSQ,MNLSQ)
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
      ZOMG(I,J,IG,IS)=ZOM1(I,J)
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
      DO 350 IG=1,NG
C
      DO 3516 J=1,NLSQ
      DO 3517 I=1,NLSQ
      ZOM1(I,J)=ZOMG(I,J,IG,IS)
3517  CONTINUE
3516  CONTINUE
C
      DO 353 IC=1,NC(IG)
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
CXXX    XDOS    ****
C*******************
      SUBROUTINE XDOS
C
C************************************************
C      CALCULATES LOCAL DOS BY MEANS OF
C      ANALYTIC CONTINUATION OF CONDITIONALLY
C      AVERAGED PHYSICAL GREEN'S FUNCTIONS
C************************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZA=4)
C
      DIMENSION ZWG(MNE),EW(MNE),
     &          DOSP(MNE,MNL,MNS,MNA),
     &          DOST(MNE,MNS,MNA)
C
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/GDOS/ ZGDOS(MNL,MNA,MNE,MNS)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
      PARAMETER(MNE= 512)
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
C     ( FOR THE BLOCH SPECTRAL FUNCTIONS )
C***********************************************
C
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,P-Y)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT CHARACTER*16 (O)
C
      PARAMETER(MNP=24)
      PARAMETER(MNPP1=MNP+1)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNL=3)
      PARAMETER(MNBZA=4)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
      PARAMETER(MNCL=32)
C
      DIMENSION TRAN(2)
C
      COMMON/SUV/ VBR(2,2),VBG(2,2),ROMEGA,GOMEGA
      COMMON/KMESHA/ AKBZA(2,MNBZA),NBZA
      COMMON/SCR/ STR(MNLSQ,MNLSQ,MNCL,MNB,0:MNPP1),
     &            JDRP(MNCL,MNB,0:MNPP1),JBVA(MNCL,MNB,0:MNPP1),
     &            JTRA(2,MNCL,MNB,0:MNPP1),NSCL(MNB,0:MNPP1)
      COMMON/BSCR/ BSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPB(MNCL,MNB),JBVAB(MNCL,MNB),
     &             JTRAB(2,MNCL,MNB),NSCLB(MNB)
      COMMON/VSCR/ VSTR(MNLSQ,MNLSQ,MNCL,MNB),
     &             JDRPV(MNCL,MNB),JBVAV(MNCL,MNB),
     &             JTRAV(2,MNCL,MNB),NSCLV(MNB)
      COMMON/SCKA/ ZSKIA(MNBLSQ,MNBLSQ,MNP,MNBZA),
     &             ZSKOA(MNBLSQ,MNBLSQ,MNP,MNBZA),
     &             ZSKFA(MNBLSQ,MNBLSQ,MNP,MNBZA)
      COMMON/BSCKA/ ZBSKIA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZBSKOA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZBSKFA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/VSCKA/ ZVSKIA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZVSKOA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZVSKFA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
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
     &    14X,'( FOR THE BLOCH SPECTRAL FUNCTIONS )'/
     &   /10X,'    CHECK OF HERMITICITY : ')
C
C----------------------------------- INTERMEDIATE REGION
C
      DO 3101 IBZA=1,MNBZA
      DO 3102 IP=1,MNP
      DO 3103 II=1,MNBLSQ
      DO 3104 JJ=1,MNBLSQ
      ZSKIA(JJ,II,IP,IBZA)=ZZ
      ZSKOA(JJ,II,IP,IBZA)=ZZ
      ZSKFA(JJ,II,IP,IBZA)=ZZ
3104  CONTINUE
3103  CONTINUE
3102  CONTINUE
3101  CONTINUE
C
      DO 311 IBZA=1,NBZA
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
      ARG=ARG+AKBZA(I,IBZA)*TRAN(I)
315   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      IF(JDIR.EQ.0) THEN
      DO 3161 IQ=1,NLSQ
      DO 3162 JQ=1,NLSQ
      ZSKIA(JJ0+JQ,II0+IQ,IP,IBZA)=ZSKIA(JJ0+JQ,II0+IQ,IP,IBZA)
     &                            + ZEFA * STR(JQ,IQ,JCL,IB,IP)
3162  CONTINUE
3161  CONTINUE
      END IF
      IF(JDIR.EQ.1) THEN
      DO 3171 IQ=1,NLSQ
      DO 3172 JQ=1,NLSQ
      ZSKFA(JJ0+JQ,II0+IQ,IP,IBZA)=ZSKFA(JJ0+JQ,II0+IQ,IP,IBZA)
     &                            + ZEFA * STR(JQ,IQ,JCL,IB,IP)
3172  CONTINUE
3171  CONTINUE
      END IF
      IF(JDIR.EQ.-1) THEN
      DO 3181 IQ=1,NLSQ
      DO 3182 JQ=1,NLSQ
      ZSKOA(JJ0+JQ,II0+IQ,IP,IBZA)=ZSKOA(JJ0+JQ,II0+IQ,IP,IBZA)
     &                            + ZEFA * STR(JQ,IQ,JCL,IB,IP)
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
      DO 321 IBZA=1,NBZA
      DO 322 II=1,NBLSQ
      DO 323 JJ=II,NBLSQ
      ZW1=ZSKIA(II,JJ,IP,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKIA(JJ,II,IP,IBZA)
      W2X=DREAL(ZW2)
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
      DO 326 IBZA=1,NBZA
      DO 327 II=1,NBLSQ
      DO 328 JJ=1,NBLSQ
      ZW1=ZSKFA(II,JJ,IP,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZSKOA(JJ,II,IPP1,IBZA)
      W2X=DREAL(ZW2)
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
      DO 3301 IBZA=1,MNBZA
      DO 3302 II=1,MNBLSQ
      DO 3303 JJ=1,MNBLSQ
      ZBSKIA(JJ,II,IBZA)=ZZ
      ZBSKOA(JJ,II,IBZA)=ZZ
      ZBSKFA(JJ,II,IBZA)=ZZ
3303  CONTINUE
3302  CONTINUE
3301  CONTINUE
C
      DO 331 IBZA=1,NBZA
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
      ARG=ARG+AKBZA(I,IBZA)*TRAN(I)
335   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      IF(JDIR.EQ.0) THEN
      DO 3361 IQ=1,NLSQ
      DO 3362 JQ=1,NLSQ
      ZBSKIA(JJ0+JQ,II0+IQ,IBZA) = ZBSKIA(JJ0+JQ,II0+IQ,IBZA)
     &                           + ZEFA * BSTR(JQ,IQ,JCL,IB)
3362  CONTINUE
3361  CONTINUE
      END IF
      IF(JDIR.EQ.1) THEN
      DO 3371 IQ=1,NLSQ
      DO 3372 JQ=1,NLSQ
      ZBSKFA(JJ0+JQ,II0+IQ,IBZA) = ZBSKFA(JJ0+JQ,II0+IQ,IBZA)
     &                           + ZEFA * BSTR(JQ,IQ,JCL,IB)
3372  CONTINUE
3371  CONTINUE
      END IF
      IF(JDIR.EQ.-1) THEN
      DO 3381 IQ=1,NLSQ
      DO 3382 JQ=1,NLSQ
      ZBSKOA(JJ0+JQ,II0+IQ,IBZA) = ZBSKOA(JJ0+JQ,II0+IQ,IBZA)
     &                           + ZEFA * BSTR(JQ,IQ,JCL,IB)
3382  CONTINUE
3381  CONTINUE
      END IF
C
333   CONTINUE
332   CONTINUE
331   CONTINUE
C                       CHECK OF HERMITICITY
      DEV=RCZ
      DO 341 IBZA=1,NBZA
      DO 342 II=1,NBLSQ
      DO 343 JJ=II,NBLSQ
      ZW1=ZBSKIA(II,JJ,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZBSKIA(JJ,II,IBZA)
      W2X=DREAL(ZW2)
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
      DO 346 IBZA=1,NBZA
      DO 347 II=1,NBLSQ
      DO 348 JJ=1,NBLSQ
      ZW1=ZBSKFA(II,JJ,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZBSKOA(JJ,II,IBZA)
      W2X=DREAL(ZW2)
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
      DO 3501 IBZA=1,MNBZA
      DO 3502 II=1,MNBLSQ
      DO 3503 JJ=1,MNBLSQ
      ZVSKIA(JJ,II,IBZA)=ZZ
      ZVSKOA(JJ,II,IBZA)=ZZ
      ZVSKFA(JJ,II,IBZA)=ZZ
3503  CONTINUE
3502  CONTINUE
3501  CONTINUE
C
      DO 351 IBZA=1,NBZA
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
      ARG=ARG+AKBZA(I,IBZA)*TRAN(I)
355   CONTINUE
      XEFA=COS(ARG)
      YEFA=SIN(ARG)
      ZEFA=DCMPLX(XEFA,YEFA)
C
      IF(JDIR.EQ.0) THEN
      DO 3561 IQ=1,NLSQ
      DO 3562 JQ=1,NLSQ
      ZVSKIA(JJ0+JQ,II0+IQ,IBZA) = ZVSKIA(JJ0+JQ,II0+IQ,IBZA)
     &                           + ZEFA * VSTR(JQ,IQ,JCL,IB)
3562  CONTINUE
3561  CONTINUE
      END IF
      IF(JDIR.EQ.1) THEN
      DO 3571 IQ=1,NLSQ
      DO 3572 JQ=1,NLSQ
      ZVSKFA(JJ0+JQ,II0+IQ,IBZA) = ZVSKFA(JJ0+JQ,II0+IQ,IBZA)
     &                           + ZEFA * VSTR(JQ,IQ,JCL,IB)
3572  CONTINUE
3571  CONTINUE
      END IF
      IF(JDIR.EQ.-1) THEN
      DO 3581 IQ=1,NLSQ
      DO 3582 JQ=1,NLSQ
      ZVSKOA(JJ0+JQ,II0+IQ,IBZA) = ZVSKOA(JJ0+JQ,II0+IQ,IBZA)
     &                           + ZEFA * VSTR(JQ,IQ,JCL,IB)
3582  CONTINUE
3581  CONTINUE
      END IF
C
353   CONTINUE
352   CONTINUE
351   CONTINUE
C                       CHECK OF HERMITICITY
      DEV=RCZ
      DO 361 IBZA=1,NBZA
      DO 362 II=1,NBLSQ
      DO 363 JJ=II,NBLSQ
      ZW1=ZVSKIA(II,JJ,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZVSKIA(JJ,II,IBZA)
      W2X=DREAL(ZW2)
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
      DO 366 IBZA=1,NBZA
      DO 367 II=1,NBLSQ
      DO 368 JJ=1,NBLSQ
      ZW1=ZVSKFA(II,JJ,IBZA)
      W1X=DREAL(ZW1)
      W1Y=DIMAG(ZW1)
      ZW2=ZVSKOA(JJ,II,IBZA)
      W2X=DREAL(ZW2)
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
CXXX    GBSF    ****
C*******************
      SUBROUTINE GBSF(IE)
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZA=4)
      PARAMETER(MNSB=2)
      PARAMETER(MNSV=2)
      PARAMETER(MNLSQ=MNL**2)
      PARAMETER(MNBLSQ=MNB*MNLSQ)
C
      DIMENSION ZWD(MNLSQ,MNLSQ),ZWF(MNLSQ,MNLSQ),
     &          ZWG(MNLSQ,MNLSQ)
C
      DIMENSION ZPSI1(MNBLSQ,MNBLSQ),ZWS(MNBLSQ,MNBLSQ),
     &          ZSKF1(MNBLSQ,MNBLSQ),ZSKO1(MNBLSQ,MNBLSQ),
     &          ZWA(MNBLSQ,MNBLSQ),ZWB(MNBLSQ,MNBLSQ)
C
      DIMENSION ZCMU(MNLSQ,MNLSQ,MNG,MNS),
     &          ZCLA(MNLSQ,MNLSQ,MNG,MNS)
C
      DIMENSION ZPSI(MNBLSQ,MNBLSQ,MNP),
     &          ZGMB(MNBLSQ,MNBLSQ,MNP),ZGMV(MNBLSQ,MNBLSQ,MNP)
C
      COMMON/KMESHA/ AKBZA(2,MNBZA),NBZA
      COMMON/SCKA/ ZSKIA(MNBLSQ,MNBLSQ,MNP,MNBZA),
     &             ZSKOA(MNBLSQ,MNBLSQ,MNP,MNBZA),
     &             ZSKFA(MNBLSQ,MNBLSQ,MNP,MNBZA)
      COMMON/BSCKA/ ZBSKIA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZBSKOA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZBSKFA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/VSCKA/ ZVSKIA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZVSKOA(MNBLSQ,MNBLSQ,MNBZA),
     &              ZVSKFA(MNBLSQ,MNBLSQ,MNBZA)
      COMMON/BCPF/ ZBCPF(MNLSQ,MNLSQ,MNB,MNSB)
      COMMON/VCPF/ ZVCPF(MNLSQ,MNLSQ,MNB,MNSV)
      COMMON/GFG/ ZGFG(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/CPF/ ZCPF(MNLSQ,MNLSQ,MNG,MNS)
      COMMON/PFALL/ ZPF(MNLSQ,MNA,MNS),
     &              ZMU(MNLSQ,MNA,MNS),
     &              ZLA(MNLSQ,MNA,MNS)
      COMMON/GBLO/ ZGBLO(MNL,MNG,MNBZA,MNE,MNS)
      COMMON/BSURA/ ZBSURA(MNBLSQ,MNBLSQ,MNBZA,MNSB)
      COMMON/VSURA/ ZVSURA(MNBLSQ,MNBLSQ,MNBZA,MNSV)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/DATCHE/ CON(MNA),AZ(MNA),WS(MNA),WSAV(MNA),
     &               VALZ(MNA)
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13,
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/ITSGF/ DIFMS,QMIXS,NMITS
C
      DATA RCZ/0.0D0/,RC1/1.0D0/
C
      NLSQ=NL**2
      NBLSQ=NB*NLSQ
      NG=NB*NP
      ZZ=DCMPLX(RCZ,RCZ)
      Z1=DCMPLX(RC1,RCZ)
C
      DO 3001 IS=1,MNS
      DO 3003 IG=1,MNG
      DO 3004 J=1,MNLSQ
      DO 3005 I=1,MNLSQ
      ZCMU(I,J,IG,IS)=ZZ
3005  CONTINUE
3004  CONTINUE
3003  CONTINUE
3001  CONTINUE
      DO 3011 IS=1,MNS
      DO 3013 IG=1,MNG
      DO 3014 J=1,MNLSQ
      DO 3015 I=1,MNLSQ
      ZCLA(I,J,IG,IS)=ZZ
3015  CONTINUE
3014  CONTINUE
3013  CONTINUE
3011  CONTINUE
C
C---------------------- CALIGRAPHIC QUANTITIES MU AND LAMBDA
C
C.................................. SPIN-LOOP
      DO 3102 IS=1,NS
C.................................. SITE-LOOP
      IA=0
      DO 3106 IG=1,NG
C
      DO 3120 J=1,NLSQ
      DO 3121 I=1,NLSQ
      ZWG(I,J)=ZGFG(I,J,IG,IS)
3121  CONTINUE
3120  CONTINUE
C.................................. SPECIES-LOOP
      DO 3124 IC=1,NC(IG)
      IA=IA+1
C                                    MATRIX F
      DO 3127 J=1,NLSQ
      DO 3128 I=1,NLSQ
      ZWD(I,J)=-ZCPF(I,J,IG,IS)
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
      ZCMU(I,J,IG,IS) = ZCMU(I,J,IG,IS)
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
      ZCLA(I,J,IG,IS) = ZCLA(I,J,IG,IS)
     &                + CON(IA) * ZWD(I,J)
3153  CONTINUE
3152  CONTINUE
C.............................. END OF SPECIES-LOOP
3124  CONTINUE
C                                  CALIGRAPHIC LAMBDA
      DO 3160 J=1,NLSQ
      DO 3161 I=1,NLSQ
      ZWF(I,J)=ZCMU(I,J,IG,IS)
3161  CONTINUE
3160  CONTINUE
C
      CALL MAMU(ZWG,ZWF,ZWD,NLSQ,MNLSQ)
C
      DO 3164 J=1,NLSQ
      DO 3165 I=1,NLSQ
      ZWF(J,I)=ZCMU(I,J,IG,IS)
3165  CONTINUE
3164  CONTINUE
C
      CALL MAMU(ZWF,ZWD,ZWG,NLSQ,MNLSQ)
C
      DO 3168 J=1,NLSQ
      DO 3169 I=1,NLSQ
      ZCLA(I,J,IG,IS)=ZCLA(I,J,IG,IS)-ZWG(I,J)
3169  CONTINUE
3168  CONTINUE
C.............................. END OF SITE-LOOP
3106  CONTINUE
C.............................. END OF SPIN-LOOP
3102  CONTINUE
C
C------------------------------ GREEN FUNCTION QUANTITIES
C
C::::::::::::::::::::::::::::::::::: BULK-SIDE SGF
      DO 3691 ISB=1,NSB
      DO 3692 IBZA=1,NBZA
C
      ITS=0
      IF(IE.EQ.1) THEN
      DO 3681 JJ=1,MNBLSQ
      DO 3682 II=1,MNBLSQ
      ZBSURA(II,JJ,IBZA,ISB)=ZZ
3682  CONTINUE
3681  CONTINUE
      END IF
      DO 3601 JJ=1,NBLSQ
      DO 3602 II=1,NBLSQ
      ZWS(II,JJ)=ZBSURA(II,JJ,IBZA,ISB)
3602  CONTINUE
3601  CONTINUE
      DO 3604 JJ=1,NBLSQ
      DO 3605 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZBSKIA(II,JJ,IBZA)
3605  CONTINUE
3604  CONTINUE
      DO 3607 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3608 J=1,NLSQ
      DO 3609 I=1,NLSQ
      ZPSI1(II0+I,II0+J) = ZPSI1(II0+I,II0+J)
     &                   + ZBCPF(I,J,IB,ISB)
3609  CONTINUE
3608  CONTINUE
3607  CONTINUE
      DO 3611 JJ=1,NBLSQ
      DO 3612 II=1,NBLSQ
      ZSKF1(II,JJ)=ZBSKFA(II,JJ,IBZA)
3612  CONTINUE
3611  CONTINUE
      DO 3614 JJ=1,NBLSQ
      DO 3615 II=1,NBLSQ
      ZSKO1(II,JJ)=ZBSKOA(II,JJ,IBZA)
3615  CONTINUE
3614  CONTINUE
C
201   ITS=ITS+1
      IF(ITS.GT.NMITS) GO TO 291
      CALL MAMU(ZWS,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3621 JJ=1,NBLSQ
      DO 3622 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZWA(II,JJ)
3622  CONTINUE
3621  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 3624 JJ=1,NBLSQ
      DO 3625 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZWS(II,JJ)
3625  CONTINUE
3624  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZWS,NBLSQ,MNBLSQ)
      DO 3628 JJ=1,NBLSQ
      DO 3629 II=1,NBLSQ
      ZWS(II,JJ)=ZWS(II,JJ)+QMIXS*ZWA(II,JJ)
3629  CONTINUE
3628  CONTINUE
      IF(SUM.GT.(DIFMS*SUM1)) GO TO 201
C
      DO 3631 JJ=1,NBLSQ
      DO 3632 II=1,NBLSQ
      ZBSURA(II,JJ,IBZA,ISB)=ZWS(II,JJ)
3632  CONTINUE
3631  CONTINUE
C
3692  CONTINUE
3691  CONTINUE
C
C::::::::::::::::::::::::::::::::::: VACUUM-SIDE SGF
      DO 3791 ISV=1,NSV
      DO 3792 IBZA=1,NBZA
C
      ITS=0
      IF(IE.EQ.1) THEN
      DO 3781 JJ=1,MNBLSQ
      DO 3782 II=1,MNBLSQ
      ZVSURA(II,JJ,IBZA,ISV)=ZZ
3782  CONTINUE
3781  CONTINUE
      END IF
      DO 3701 JJ=1,NBLSQ
      DO 3702 II=1,NBLSQ
      ZWS(II,JJ)=ZVSURA(II,JJ,IBZA,ISV)
3702  CONTINUE
3701  CONTINUE
      DO 3704 JJ=1,NBLSQ
      DO 3705 II=1,NBLSQ
      ZPSI1(II,JJ)=-ZVSKIA(II,JJ,IBZA)
3705  CONTINUE
3704  CONTINUE
      DO 3707 IB=1,NB
      II0=(IB-1)*NLSQ
      DO 3708 J=1,NLSQ
      DO 3709 I=1,NLSQ
      ZPSI1(II0+I,II0+J) = ZPSI1(II0+I,II0+J)
     &                   + ZVCPF(I,J,IB,ISV)
3709  CONTINUE
3708  CONTINUE
3707  CONTINUE
      DO 3711 JJ=1,NBLSQ
      DO 3712 II=1,NBLSQ
      ZSKF1(II,JJ)=ZVSKFA(II,JJ,IBZA)
3712  CONTINUE
3711  CONTINUE
      DO 3714 JJ=1,NBLSQ
      DO 3715 II=1,NBLSQ
      ZSKO1(II,JJ)=ZVSKOA(II,JJ,IBZA)
3715  CONTINUE
3714  CONTINUE
C
202   ITS=ITS+1
      IF(ITS.GT.NMITS) GO TO 292
      CALL MAMU(ZWS,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3721 JJ=1,NBLSQ
      DO 3722 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI1(II,JJ)-ZWA(II,JJ)
3722  CONTINUE
3721  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      DO 3724 JJ=1,NBLSQ
      DO 3725 II=1,NBLSQ
      ZWA(II,JJ)=ZWA(II,JJ)-ZWS(II,JJ)
3725  CONTINUE
3724  CONTINUE
      SUM=CMANO(ZWA,NBLSQ,MNBLSQ)
      SUM1=CMANO(ZWS,NBLSQ,MNBLSQ)
      DO 3728 JJ=1,NBLSQ
      DO 3729 II=1,NBLSQ
      ZWS(II,JJ)=ZWS(II,JJ)+QMIXS*ZWA(II,JJ)
3729  CONTINUE
3728  CONTINUE
      IF(SUM.GT.(DIFMS*SUM1)) GO TO 202
C
      DO 3731 JJ=1,NBLSQ
      DO 3732 II=1,NBLSQ
      ZVSURA(II,JJ,IBZA,ISV)=ZWS(II,JJ)
3732  CONTINUE
3731  CONTINUE
C
3792  CONTINUE
3791  CONTINUE
C
C:::::::::::::::::::::::::::::::: RECURSION PARTITIONING
      DO 3891 IS=1,NS
      ISB=MIN(IS,NSB)
      ISV=MIN(IS,NSV)
      DO 3892 IBZA=1,NBZA
C
      DO 3801 IP=1,NP
      DO 3802 JJ=1,NBLSQ
      DO 3803 II=1,NBLSQ
      ZPSI(II,JJ,IP)=-ZSKIA(II,JJ,IP,IBZA)
3803  CONTINUE
3802  CONTINUE
3801  CONTINUE
      IG=0
      DO 3805 IP=1,NP
      DO 3806 IB=1,NB
      IG=IG+1
      II0=(IB-1)*NLSQ
      DO 3808 J=1,NLSQ
      DO 3809 I=1,NLSQ
      ZPSI(II0+I,II0+J,IP) = ZPSI(II0+I,II0+J,IP)
     &                     + ZCPF(I,J,IG,IS)
3809  CONTINUE
3808  CONTINUE
3806  CONTINUE
3805  CONTINUE
C                             RECURSION FROM VACUUM
      DO 3811 JJ=1,NBLSQ
      DO 3812 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKOA(II,JJ,1,IBZA)
3812  CONTINUE
3811  CONTINUE
      DO 3813 JJ=1,NBLSQ
      DO 3814 II=1,NBLSQ
      ZSKF1(II,JJ)=DCONJG(ZSKO1(JJ,II))
3814  CONTINUE
3813  CONTINUE
      DO 3815 JJ=1,NBLSQ
      DO 3816 II=1,NBLSQ
      ZWA(II,JJ)=ZVSURA(II,JJ,IBZA,ISV)
3816  CONTINUE
3815  CONTINUE
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3817 JJ=1,NBLSQ
      DO 3818 II=1,NBLSQ
      ZGMV(II,JJ,1)=ZWA(II,JJ)
3818  CONTINUE
3817  CONTINUE
C
      IF(NP.GT.1) THEN
      DO 382 IP=2,NP
      IP1=IP-1
      DO 3821 JJ=1,NBLSQ
      DO 3822 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI(II,JJ,IP1)-ZGMV(II,JJ,IP1)
3822  CONTINUE
3821  CONTINUE
      DO 3823 JJ=1,NBLSQ
      DO 3824 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKFA(II,JJ,IP1,IBZA)
3824  CONTINUE
3823  CONTINUE
      DO 3825 JJ=1,NBLSQ
      DO 3826 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKOA(II,JJ,IP,IBZA)
3826  CONTINUE
3825  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKO1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKF1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3827 JJ=1,NBLSQ
      DO 3828 II=1,NBLSQ
      ZGMV(II,JJ,IP)=ZWA(II,JJ)
3828  CONTINUE
3827  CONTINUE
382   CONTINUE
      END IF
C                             RECURSION FROM BULK
      DO 3831 JJ=1,NBLSQ
      DO 3832 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKFA(II,JJ,NP,IBZA)
3832  CONTINUE
3831  CONTINUE
      DO 3833 JJ=1,NBLSQ
      DO 3834 II=1,NBLSQ
      ZSKO1(II,JJ)=DCONJG(ZSKF1(JJ,II))
3834  CONTINUE
3833  CONTINUE
      DO 3835 JJ=1,NBLSQ
      DO 3836 II=1,NBLSQ
      ZWA(II,JJ)=ZBSURA(II,JJ,IBZA,ISB)
3836  CONTINUE
3835  CONTINUE
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3837 JJ=1,NBLSQ
      DO 3838 II=1,NBLSQ
      ZGMB(II,JJ,NP)=ZWA(II,JJ)
3838  CONTINUE
3837  CONTINUE
C
      IF(NP.GT.1) THEN
      DO 384 IP=NP-1,1,-1
      IP1=IP+1
      DO 3841 JJ=1,NBLSQ
      DO 3842 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI(II,JJ,IP1)-ZGMB(II,JJ,IP1)
3842  CONTINUE
3841  CONTINUE
      DO 3843 JJ=1,NBLSQ
      DO 3844 II=1,NBLSQ
      ZSKF1(II,JJ)=ZSKFA(II,JJ,IP,IBZA)
3844  CONTINUE
3843  CONTINUE
      DO 3845 JJ=1,NBLSQ
      DO 3846 II=1,NBLSQ
      ZSKO1(II,JJ)=ZSKOA(II,JJ,IP1,IBZA)
3846  CONTINUE
3845  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
      CALL MAMU(ZWA,ZSKF1,ZWB,NBLSQ,MNBLSQ)
      CALL MAMU(ZSKO1,ZWB,ZWA,NBLSQ,MNBLSQ)
      DO 3847 JJ=1,NBLSQ
      DO 3848 II=1,NBLSQ
      ZGMB(II,JJ,IP)=ZWA(II,JJ)
3848  CONTINUE
3847  CONTINUE
384   CONTINUE
      END IF
C
      IG=0
      DO 321 IP=1,NP
      DO 3214 JJ=1,NBLSQ
      DO 3215 II=1,NBLSQ
      ZWA(II,JJ)=ZPSI(II,JJ,IP)-ZGMV(II,JJ,IP)-ZGMB(II,JJ,IP)
3215  CONTINUE
3214  CONTINUE
      CALL CINVC(ZWA,NBLSQ,MNBLSQ)
C.................................... SITE-LOOP
      DO 322 IB=1,NB
      IG=IG+1
      II0=(IB-1)*NLSQ
      DO 3224 J=1,NLSQ
      DO 3225 I=1,NLSQ
      ZWG(I,J)=ZWA(II0+I,II0+J)
3225  CONTINUE
3224  CONTINUE
C
      DO 3228 J=1,NLSQ
      DO 3229 I=1,NLSQ
      ZWF(J,I)=ZCMU(I,J,IG,IS)
3229  CONTINUE
3228  CONTINUE
C
      CALL MAMU(ZWF,ZWG,ZWD,NLSQ,MNLSQ)
C
      DO 3242 J=1,NLSQ
      DO 3243 I=1,NLSQ
      ZWF(I,J)=ZCMU(I,J,IG,IS)
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
      ZSUM = ZSUM + ZWG(I,I) + ZCLA(I,I,IG,IS)
3257  CONTINUE
      ZGBLO(IL,IG,IBZA,IE,IS)=ZSUM
325   CONTINUE
C.................................. END OF SITE LOOP
322   CONTINUE
321   CONTINUE
C
3892  CONTINUE
3891  CONTINUE
C
      RETURN
C
C------------------------------ NO CONVERGENCY
C
291   WRITE(IW6,191) NMITS,IE,ISB,IBZA
191   FORMAT(//' **** GBSF: BULK SGF NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   5X,'---  FOR IE=',I4,5X,'ISB=',I1,5X,'IBZA=',I2)
      STOP
292   WRITE(IW6,192) NMITS,IE,ISV,IBZA
192   FORMAT(//' **** GBSF: VAC. SGF NON-CONVERGED ',
     &  'AFTER ',I4,' ITERATIONS'/
     &   5X,'---  FOR IE=',I4,5X,'ISV=',I1,5X,'IBZA=',I2)
      STOP
C
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
      PARAMETER(MNP=24)
      PARAMETER(MNB=2)
      PARAMETER(MNG=MNP*MNB)
      PARAMETER(MNA=40)
      PARAMETER(MNL=3)
      PARAMETER(MNS=2)
      PARAMETER(MNE= 512)
      PARAMETER(MNBZA=4)
C
      DIMENSION ZWX(MNE),EW(MNE),
     &          ABLOP(MNE,MNL,MNS,MNG,MNBZA),
     &          ABLOT(MNE,MNS,MNG,MNBZA)
C
      COMMON/KMESHA/ AKBZA(2,MNBZA),NBZA
      COMMON/EMESH/ ZCN(MNE),EBOT,ETOP,EH,NE,NDC
      COMMON/GBLO/ ZGBLO(MNL,MNG,MNBZA,MNE,MNS)
      COMMON/DATDIM/ IVAC,NP,NB,NC(MNG),NA,NAB,NAV,NL,NS,NSB,NSV
      COMMON/REWR/ IR1,IR2,IR4,IR5,IR11,IR13, 
     &             IW6,IW7,IW8,IW9,IW17,IW18,IW19
      COMMON/TEXT/ OTXTA(MNA),OTXTK(MNBZA)
C
      DATA RCZ/0.0D0/,RC1/1.0D0/,RC4/4.0D0/
C
      NG=NB*NP
      PI=RC4*ATAN(RC1)
C
      DO 3041 IBZA=1,MNBZA
      DO 3042 IG=1,MNG
      DO 3043 IS=1,MNS
      DO 3044 IE=1,MNE
      ABLOT(IE,IS,IG,IBZA)=RCZ
3044  CONTINUE
3043  CONTINUE
3042  CONTINUE
3041  CONTINUE
C
!SMP$ DO SERIAL
      DO 308 IE=1,NE
      EW(IE)=DREAL(ZCN(IE))
308   CONTINUE
C
C------------------------------ BLOCH SPECTRAL FUNCTIONS
C
C                                      PARTIAL
      DO 3302 IBZA=1,NBZA
      DO 3303 IG=1,NG
      DO 3304 IS=1,NS
      DO 3305 IL=1,NL
C
!SMP$ DO SERIAL
      DO 3314 IE=1,NE
      ZWX(IE)=ZGBLO(IL,IG,IBZA,IE,IS)
3314  CONTINUE
C
      CALL DECO(ZWX,NE,NDC)
C
      DO 3316 IE=1,NE
      ABLOP(IE,IL,IS,IG,IBZA) = - DIMAG(ZWX(IE))/PI
3316  CONTINUE
C
3305  CONTINUE
3304  CONTINUE
3303  CONTINUE
3302  CONTINUE
C
C                                       TOTAL
      DO 3321 IBZA=1,NBZA
      DO 3322 IG=1,NG
      DO 3323 IS=1,NS
      DO 3324 IL=1,NL
      DO 3325 IE=1,NE
      ABLOT(IE,IS,IG,IBZA) = ABLOT(IE,IS,IG,IBZA)
     &                     + ABLOP(IE,IL,IS,IG,IBZA)
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
      WRITE(IW17,101) NE,NL,NS,NG,NBZA
      DO 340 IBZA=1,NBZA
      DO 341 IG=1,NG
      DO 342 IS=1,NS
      WRITE(IW17,132) OTXTK(IBZA),IBZA,IG,IS
132   FORMAT(4X,A16,4X,'(IBZA=',I2,')',
     &              4X,'IG=',I4,4X,'IS=',I1)
      DO 343 IE=1,NE
      WRITE(IW17,106) EW(IE),ABLOT(IE,IS,IG,IBZA),
     &               (ABLOP(IE,IL,IS,IG,IBZA),IL=1,NL)
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
      WRITE(IWW,101) NE,NG,NBZA
      DO 352 IBZA=1,NBZA
      WRITE(IWW,142) OTXTK(IBZA),IBZA,IS
142   FORMAT(4X,A16,4X,'(IBZA=',I2,')',4X,'IS=',I1)
      DO 353 IE=1,NE
      WRITE(IWW,106) EW(IE),(ABLOT(IE,IS,IG,IBZA),IG=1,NG)
353   CONTINUE
352   CONTINUE
      WRITE(IWW,145)
145   FORMAT(2X,'---------------------------------')
351   CONTINUE
C
      RETURN
      END

