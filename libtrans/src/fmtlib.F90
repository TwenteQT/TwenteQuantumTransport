!+******************************************************************************
! FMTLIB.F90
!
! 16-MAR-2000 GK created
!*******************************************************************************

 MODULE FMTLIB

 IMPLICIT NONE

 PRIVATE
 
 CHARACTER*5 ,PARAMETER,PUBLIC :: FMTLIB_REV  = '01.03'       !   
 CHARACTER*11,PARAMETER,PUBLIC :: FMTLIB_DATE = '16-MAR-2000' ! 

 INTEGER,PARAMETER    ,PRIVATE :: LF = 8  ! Kind param. for 8-byte real numbers 


 ! IF DEBUG_MESSAGES_FMTLIB id true some debug information is displayed (on standard output)

 LOGICAL      ,PUBLIC :: DEBUG_MESSAGES_FMTLIB = .TRUE. 
 INTEGER      ,PUBLIC :: LU_FMTLIB = 6


 INTEGER      ,PUBLIC :: DIGITS_FMTLIB = 8        ! default maximum number of digits     
 REAL(KIND=LF),PUBLIC :: MAXF_FMTLIB = 1.0E12_LF  ! largest real with F edit descriptor
                                                  ! instead with an E edit descriptor
 REAL(KIND=LF),PUBLIC :: MINF_FMTLIB = 1.0E-5_LF  ! smallest real with F edit descriptor
                                                  ! instead with an E edit descriptor
 PUBLIC GETFMT     ! 

! The following INTERFACE command is used to be able to use a common procedure
! name, CALL GETFMT(..) for all possible data types: 


 INTERFACE GETFMT

      MODULE PROCEDURE GETFMT_REAL_SCALAR, GETFMT_REAL_VECTOR, GETFMT_REAL_MATRIX, &
                       GETFMT_INT__SCALAR, GETFMT_INT__VECTOR,                     &
                       GETFMT_CHAR_SCALAR, GETFMT_CHAR_VECTOR,                     &
                       GETFMT_LOG__SCALAR, GETFMT_LOG__VECTOR

 END INTERFACE


 CONTAINS

!+******************************************************************************
!
! DECIMAL_DIGITS 
!
! Determines the number of decimal digits necessary to represent a real number
! to a precision EPS
!
! Example DECIMAL_DIGITS(0.123400009,1.0E-6) returns 4 
!
! 10-MAR-1999 GK created
!*******************************************************************************

 INTEGER FUNCTION  DECIMAL_DIGITS(X,EPS) 

! ------------- ARGUMENTS-------------------------------------------------------

   REAL(KIND=LF),INTENT(IN) ::  X   ! real number
   REAL(KIND=LF),INTENT(IN) ::  EPS ! smallest number which should not written as 0

!*------------- LOCAL DATA -----------------------------------------------------


   REAL(KIND=LF) :: EPSC,XC !
   INTEGER       :: D       !

     XC    = X
     EPSC  = EPS
     
     D=0
     DO
        IF (ABS(XC-DNINT(XC)) < EPSC) EXIT
        XC   = XC    * 10.0_LF
        EPSC = EPSC  * 10.0_LF
        D    = D + 1
     ENDDO

     DECIMAL_DIGITS = D
  
  END FUNCTION DECIMAL_DIGITS


!+******************************************************************************
!
! GETFMT_REAL_VECTOR, GETFMT_INT__VECTOR, GETFMT_CHAR_VECTOR, GETFMT_LOG__VECTOR
! GETFMT_REAL_SCALAR,GETFMT_INT__SCALAR,GETFMT_CHAR_SCALAR,GETFMT_LOG__SCALAR
!
! grouped together into generic procedure GETFMT
!
!+******************************************************************************
!
! GETFMT_REAL_SCALAR
!
! Determines the format specifier FMT for a real number X
!
! The F edit descriptor (simple decimal floating-point representation) is used 
! if the absolute value of X lies between MINF_FMTLIB and MAXF_FMTLIB. 
! If ABS(X) is smaller than MINF_FMTLIB or larger than MAXF_FMTLIB, the E edit
! descriptor (exponential notation) is used. 
! 
! For the F edit descriptor, depending on the numerical value of X,  the width of 
! the integer and decimal part is determined: there will be no trailing zeros 
! nor leading blanks.
!
! With the optional argument DIGITS, it is possible to specify the maximum number 
! of decimal digits to be displayed in F or E edit description. 
! For LF=4 (REAL*4) DIGITS should lie beteen 0 and 7 and for LF=8 DIGITS should 
! lie beteen 0 and 16.    
!                                                         
! The length of FMT is 8 for the F edit descriptor and 11 for the E edit descriptor

!
! Example#1: for the real 5/4, the returned format is '(F04.02)' 
!            and the real written with this format will give '1.25' 
!           
! Example#2: for the real X=9.8765 and DIGITS=2, the returned format is '(F04.02)' 
!            and the real written with this format will give '9.88' 
!
! 06-MAR-1999 GK created
!*******************************************************************************

  
 SUBROUTINE GETFMT_REAL_SCALAR(X,FMT,DIGITS)

! ------------- ARGUMENTS-------------------------------------------------------

   REAL(KIND=LF),INTENT(IN )          ::  X      ! format specifier of this number
                                                 ! will be determined
   INTEGER      ,INTENT(IN ),OPTIONAL ::  DIGITS ! maximum number of decimal digits
                                                 ! to be displayed
   CHARACTER*(*),INTENT(OUT)          ::  FMT    ! format specifier

!*------------- LOCAL DATA -----------------------------------------------------

   INTEGER :: D          ! width of decimal part
   INTEGER :: I          ! width of integer part (including sign and decimal point)
   INTEGER :: W          !  
   INTEGER :: LOC_DIGITS ! maximum number of digits  
   REAL(KIND=LF) :: XA   ! absolute value of X    
   REAL(KIND=LF) :: XN   ! "normalized X" 1<=XM<10
   REAL(KIND=LF) :: EPS  ! smallest number which will not written as 0

	   IF (PRESENT(DIGITS)) THEN
	     LOC_DIGITS=DIGITS
	   ELSE
	     LOC_DIGITS=DIGITS_FMTLIB
	   END IF
	   
       EPS = 10.0_LF**(-LOC_DIGITS)

       XA = ABS(X)

       IF ((XA==0.0d0).OR.(MINF_FMTLIB < XA .AND. XA < MAXF_FMTLIB)) THEN

         ! F edit descriptor

         D = DECIMAL_DIGITS(X,EPS)
         I = 2 + LOG10(MAX(XA+EPS,1.0_LF)) ! including decimal point
         IF (X < -EPS) I = I + 1           ! for minus sign
         WRITE(FMT,301) I+D,D
       
       ELSE 
       
         ! E edit descriptor

         XN = XA*10.0_LF**(-INT(LOG10(XA))) 
         IF  (XA<1) XN = XN*10.0_LF            ! 1<= XN < 10
         D = DECIMAL_DIGITS(XN,EPS) 
         W = D + 6
         IF (X < 0) W = W + 1                  ! for minus sign
         WRITE(FMT,302) W,D
       
       ENDIF

 301   FORMAT('(F',i2.2,'.',i2.2,')')
 302   FORMAT('(1P,E',i2.2,'.',i2.2,')')
       
  END SUBROUTINE  GETFMT_REAL_SCALAR


!+******************************************************************************
!
! GETFMT_REAL_VECTOR 
!
! Determines the format specifier FMT for an array X of real numbers. 
! X must have less than 100 elements.
!
! The F edit descriptor (simple decimal floating-point representation) is used 
! if the absolute value of X lies between MINF_FMTLIB and MAXF_FMTLIB. 
! If ABS(X) is smaller than MINF_FMTLIB or larger than MAXF_FMTLIB, the E edit
! descriptor (exponential notation) is used. 
! For MODE=COMMON or MAX, the F edit descriptor is only used, if all elements
! can be written in that format.
!
! With the optional argument DIGITS, it is possible to specify the maximum number 
! of decimal digits to be displayed in F or E edit description. 
! For LF=4 (REAL*4) DIGITS should lie beteen 0 and 7 and for LF=8 DIGITS should 
! lie beteen 0 and 16.  
!
! argument MODE can take one of the following values:
!
! 'INDIV'  each element of X has it individual edit descriptor:
!          Example: (F04.01,X,1P,E06.00,0P,X,F07.02) for X=15.5,1.E+12,-344.66
!
! 'COMMON' all elements of X have the same edit descriptor (default)
!          Example: (03(F08.03,X))                  for X=15.5,2.777,-344.66
!                   
! 'MAX'    same edit descriptor as for 'COMMON', but without repeat specification
!          Example: (F08.03)                    for X=15.5,2.777,-344.66
!                
!
! The length LEN of FMT depends on MODE: 
!
!  MODE=MAX    -> LEN=8 
!  MODE=COMMON -> LEN=14 for F edit descriptor and LEN=17 for E edit descriptor
!  MODE=INDIV  -> LEN depends on the values of X but is smaller than 15*N-1
!
! 06-MAR-1999 GK created
!*******************************************************************************
 
 SUBROUTINE GETFMT_REAL_VECTOR(X,FMT,DIGITS,MODE)

! ------------- ARGUMENTS-------------------------------------------------------

   REAL(KIND=LF),INTENT(IN )          ::  X(:)   ! array 
   CHARACTER*(*),INTENT(OUT)          ::  FMT    ! format specifier to write X
   INTEGER      ,INTENT(IN ),OPTIONAL ::  DIGITS ! maximum number of decimal digits
                                                 ! to be displayed
   CHARACTER*(*),INTENT(IN ),OPTIONAL ::  MODE   ! can be either 'INDIV' 'COMMON','MAX'

!*------------- LOCAL DATA -----------------------------------------------------

   INTEGER       :: IR         ! dummy
   INTEGER       :: D          ! width of decimal part
   INTEGER       :: DMAX       ! maximum D
   INTEGER       :: I          ! width of integer part (including sign and decimal point)
   INTEGER       :: IMAX       ! maximum I
   INTEGER       :: W          !
   INTEGER       :: WMAX       ! maximum W
   INTEGER       :: N          ! size of X
   INTEGER       :: LOC_DIGITS !
   CHARACTER*6   :: LOC_MODE   !
   CHARACTER*12  :: FMT0       ! elementary F edit descriptor: Ex.: '1P,E08.02'
   LOGICAL       :: ALL_F      ! all elements can be written with the F edit descriptor 
   LOGICAL       :: THIS_F     ! actual elements can be written with the F edit descriptor 
   REAL(KIND=LF) :: XA         ! absolute value of X    
   REAL(KIND=LF) :: XN         ! "normalized X" 1<=XM<10
   REAL(KIND=LF) :: EPS        ! smallest number which will not written as 0


	   IF (PRESENT(MODE)) THEN
	     LOC_MODE=MODE
	   ELSE
	     LOC_MODE='COMMON'
	   END IF

	   IF (PRESENT(DIGITS)) THEN
	     LOC_DIGITS=DIGITS
	   ELSE
	     LOC_DIGITS=DIGITS_FMTLIB  ! default value
	   END IF  

       N    = SIZE(X)
       FMT  = ''
       IMAX = 0
       DMAX = 0
       WMAX = 0
       EPS = 10.0_LF**(-LOC_DIGITS)

       ALL_F = ALL( (X==0.0d0) .OR. (MINF_FMTLIB < ABS(X) .AND. ABS(X) < MAXF_FMTLIB) ) 

       DO IR=1,N

         XA = ABS(X(IR))

         THIS_F = (XA==0.0d0) .OR. (MINF_FMTLIB < XA .AND. XA < MAXF_FMTLIB)

         IF ( (LOC_MODE=='INDIV' .AND. THIS_F)  .OR. ALL_F) THEN

           ! F edit descriptor

           D = DECIMAL_DIGITS(X(IR),EPS)
           I = 2 + LOG10(MAX(XA+EPS,1.0_LF)) ! including decimal point
           IF (X(IR) < -EPS) I = I + 1           ! for minus sign
           
           IF (LOC_MODE=='INDIV')  THEN
             WRITE(FMT0,301) I+D,D
             FMT=TRIM(FMT)//FMT0
             IF (IR/=N)  FMT=TRIM(FMT)//',X,'
           ELSE
             DMAX = MAX(DMAX,D)
             IMAX = MAX(IMAX,I)
           ENDIF
 
         ELSE 
       
           ! E edit descriptor
           IF (XA==0.0_LF) THEN
             D = 0 
           ELSE
             XN = XA*10.0_LF**(-INT(LOG10(XA))) 
             IF (XA<1) XN = XN*10.0_LF            ! 1<= XN < 10
             D = DECIMAL_DIGITS(XN,EPS) 
           ENDIF
           W = D + 6
                      
           IF (LOC_MODE=='INDIV')  THEN
             IF (X(IR) < 0) W = W + 1              ! for minus sign
             WRITE(FMT0,302) W,D
             FMT=TRIM(FMT)//FMT0
             IF (IR/=N)  FMT=TRIM(FMT)//',X,'
           ELSE
             DMAX = MAX(DMAX,D)
             IF (ANY(X < 0.0_LF)) THEN
               WMAX = MAX(WMAX,W+1)
             ELSE
               WMAX = MAX(WMAX,W)  
             END IF            
           ENDIF      
         ENDIF

       ENDDO

       IF (LOC_MODE=='INDIV') THEN  
         FMT='('//TRIM(FMT)//')'
       ELSEIF(LOC_MODE=='COMMON') THEN
         IF (ALL_F) THEN 
           WRITE(FMT,303) N,IMAX+DMAX,DMAX
         ELSE
           WRITE(FMT,304) N,     WMAX,DMAX
         ENDIF
       ELSEIF (LOC_MODE=='MAX') THEN  
         IF (ALL_F) THEN 
           WRITE(FMT,305)   IMAX+DMAX,DMAX
         ELSE
           WRITE(FMT,306)        WMAX,DMAX
         ENDIF
       ELSE
         WRITE(LU_FMTLIB,*)'GETFMT_REAL_VECTOR: UNKNOWN VALUE FOR ARGUMENT MODE: ',LOC_MODE
         STOP
       ENDIF     

 301   FORMAT('F',i2.2,'.',i2.2)
 302   FORMAT('1P,E',i2.2,'.',i2.2,',0P')
       
 303   FORMAT('(',i2.2,'(F'   ,i2.2,'.',i2.2,',X))')
 304   FORMAT('(',i2.2,'(1P,E',i2.2,'.',i2.2,',X))')

 305   FORMAT('(F'   ,i2.2,'.',i2.2,')')
 306   FORMAT('(1P,E',i2.2,'.',i2.2,')')
       		
END SUBROUTINE GETFMT_REAL_VECTOR

!+******************************************************************************
!
! GETFMT_REAL_MATRIX 
!
! see comments of GETFMT_REAL_VECTOR 
!-------------------------------------------------------------------------------

SUBROUTINE GETFMT_REAL_MATRIX(X,FMT,DIGITS,MODE)

! ------------- ARGUMENTS-------------------------------------------------------

   REAL(KIND=LF),INTENT(IN )          ::  X(:,:) ! array 
   CHARACTER*(*),INTENT(OUT)          ::  FMT    ! format specifier to write X
   INTEGER      ,INTENT(IN ),OPTIONAL ::  DIGITS ! maximum number of decimal digits
                                                 ! to be displayed
   CHARACTER*(*),INTENT(IN ),OPTIONAL ::  MODE   ! can be either 'INDIV' 'COMMON','MAX'

!-------------------------------------------------------

    CALL GETFMT_REAL_VECTOR(RESHAPE(X,(/SIZE(X)/)),FMT,DIGITS,MODE)

END SUBROUTINE GETFMT_REAL_MATRIX

!+******************************************************************************
! GETFMT_LOG__SCALAR 
!
! Returns the format specifier for a logical: '(L1)'
!
! 06-MAR-1999 GK created
!*******************************************************************************
  
 SUBROUTINE GETFMT_LOG__SCALAR(X,FMT)

! ------------- ARGUMENTS-------------------------------------------------------
   
   LOGICAL      ,INTENT(IN )          ::  X      !
   CHARACTER*(*),INTENT(OUT)          ::  FMT    ! format specifier to write X

!*------------- LOCAL DATA -----------------------------------------------------
   
   LOGICAL :: LDUM  ! dummy
  
       LDUM = X     ! to prevent compiler warnings
       FMT ='(L1)'
 
 END SUBROUTINE  GETFMT_LOG__SCALAR


!+******************************************************************************
! GETFMT_LOG__VECTOR 
!
! Determines the format specifier FMT for an array X of logiccals. 
! X must have less than 100 elements.
! The argument MODE can take one of the following values:
!
! MODE='COMMON' all elements of X have the same edit descriptor (default)
!               Ex: (03(L1,X))
!
! MODE='MAX'    same edit descriptor as for MODE='COMMON', but without repeat 
!               specification: (L1)
!
! 06-MAR-1999 GK created
!*******************************************************************************
  
 SUBROUTINE GETFMT_LOG__VECTOR(X,FMT,MODE)

! ------------- ARGUMENTS-------------------------------------------------------
   
   LOGICAL      ,INTENT(IN )          ::  X(:)   ! 
   CHARACTER*(*),INTENT(OUT)          ::  FMT    ! format specifier to write X
   CHARACTER*(*),INTENT(IN ),OPTIONAL ::  MODE   ! can be either 'COMMON','MAX'

!*------------- LOCAL DATA -----------------------------------------------------

   INTEGER       :: N         !
   CHARACTER*6   :: LOC_MODE  !
 

	   IF (PRESENT(MODE)) THEN
	     LOC_MODE=MODE
	   ELSE
	     LOC_MODE='COMMON'
	   END IF

	   N=SIZE(X)

       IF (LOC_MODE=='COMMON') THEN
         WRITE(FMT,300) N
       ELSEIF (LOC_MODE=='MAX') THEN  
         FMT ='(L1)'
       ELSE
         WRITE(LU_FMTLIB,*)'GETFMT_LOG__VECTOR: UNKNOWN VALUE FOR ARGUMENT MODE: ',LOC_MODE
         STOP
       ENDIF 

 300   FORMAT('(',i2.2,'L1)')
        
 END SUBROUTINE  GETFMT_LOG__VECTOR

!+******************************************************************************
! GETFMT_INT__SCALAR 
!
! Determines the format specifier for an integer number.
!
! The character width depends on the numerical value of X; there will be no 
! leading blanks.
!                                                         
! Example: for -565  the returned format is '(I04)' 
!          and the number written with this format gives '-565' 
!
! The length of FMT is 5 

! 06-MAR-1999 GK created
!*******************************************************************************
  
 SUBROUTINE GETFMT_INT__SCALAR(X,FMT)

! ------------- ARGUMENTS-------------------------------------------------------

   INTEGER   	,INTENT(IN )          ::  X   !
   CHARACTER*(*),INTENT(OUT)          ::  FMT ! format specifier to write X

!*------------- LOCAL DATA -----------------------------------------------------
   
   INTEGER   :: W   ! character width
 
       IF (X == 0)  THEN
           W = 1
	   ELSE
		   W = 1 + LOG10(0.5+ABS(X))
           IF (X < 0)  W = W + 1        ! for the minus sign
	   ENDIF       

       WRITE(FMT,300) W
       
 300   FORMAT('(I',i2.2,')')
      
 END SUBROUTINE  GETFMT_INT__SCALAR

!+******************************************************************************
! GETFMT_INT__VECTOR 
!
! Determines the format specifier FMT for an array X of integer numbers.. 
! X must have less than 100 elements.
!
! argument MODE can take one of the following values:
!
! 'INDIV'  each element of X has it individual edit descriptor:
!          Example: (I04,X,I05,X,I07)  for X=4444, 55555, -123456
!                   
! 'COMMON' all elements of X have the same edit descriptor (default):
!          Example: (03(I07,X))        for X=4444, 55555, -123456
!                   
! 'MAX'    same edit descriptor as for 'COMMON', but without repeat specification
!          Example: (I07)              for X=4444, 55555, -123456
!
! The length LEN of FMT depends on MODE:  'MAX'    -> LEN=5 
!                                         'COMMON' -> LEN=11
!                                         'INDIV'  -> LEN=6*N-1
! 06-MAR-1999 GK created
!*******************************************************************************
  
 SUBROUTINE GETFMT_INT__VECTOR(X,FMT,MODE)

! ------------- ARGUMENTS-------------------------------------------------------

   INTEGER   	,INTENT(IN )          :: X(:) !
   CHARACTER*(*),INTENT(OUT)          :: FMT  ! format specifier to write X
   CHARACTER*(*),INTENT(IN ),OPTIONAL :: MODE ! can be either 'INDIV' 'COMMON','MAX'

!*------------- LOCAL DATA -----------------------------------------------------
  
   INTEGER     :: W        ! character width 
   INTEGER     :: WMAX     ! maximum of W
   INTEGER     :: IR       ! dummy
   INTEGER     :: N        ! size of X
   CHARACTER*3 :: FMT0     ! elementary I edit descriptor. Ex.: 'I06'
   CHARACTER*6 :: LOC_MODE !
 

	   IF (PRESENT(MODE)) THEN
	     LOC_MODE=MODE
	   ELSE
	     LOC_MODE='COMMON'
	   END IF

       N    = SIZE(X)
       FMT  = ''
       WMAX = 0

       DO IR=1,N
         IF (X(IR) == 0)  THEN
           W = 1
	     ELSE
		   W = 1 + LOG10(0.5+ABS(X(IR)))
           IF (X(IR) < 0)  W = W + 1    ! for the minus sign
	     ENDIF         
         WMAX = MAX (WMAX,W)
         
         IF (LOC_MODE=='INDIV')  THEN
           WRITE(FMT0,300)W
           FMT=TRIM(FMT)//FMT0
           IF (IR/=N)  FMT=TRIM(FMT)//',X,'
         ENDIF 

       ENDDO

       IF (LOC_MODE=='INDIV') THEN  
         FMT='('//TRIM(FMT)//')'
       ELSEIF(LOC_MODE=='COMMON') THEN
         WRITE(FMT,301) N,WMAX
       ELSEIF (LOC_MODE=='MAX') THEN  
         WRITE(FMT,302)   WMAX
       ELSE
         WRITE(LU_FMTLIB,*)'GETFMT_INT__VECTOR: UNKNOWN VALUE FOR ARGUMENT MODE: ',LOC_MODE
         STOP
       ENDIF     

 300   FORMAT(          'I',i2.2)
 301   FORMAT('(',i2.2,'(I',i2.2,',X))')
 302   FORMAT('(I'         ,i2.2,')')
    
       		     
 END SUBROUTINE  GETFMT_INT__VECTOR

!+******************************************************************************
! GETFMT_CHAR_SCALAR 
!
! Determines the format specifier FMT for a character string X
!
! The character width is the length of X, without counting trailing blanks.  
! If the optional argument APOSTROPHES is true and if blanks are included in the 
! string, additional apostrophes will be printed if X is written with format FMT.
!
! If the optional argument APOSTROPHES is true and if X contains blanks the entire
! character string written with format FMT will be enclosed by apostrophes.
!
! The length LEN of FMT is either 6 (no apostrophes) or 16 (with apostrophes)
!                                                        
! Example: for X='Hello world ' and APOSTROPHES=.TRUE.,
!          the returned format is:  '('''',A011,'''')' and the string written with 
!          this format gives 'Hello world' (the apostrophes are really printed.)
!
! 
! 06-MAR-1999 GK created
!*******************************************************************************
  
 SUBROUTINE GETFMT_CHAR_SCALAR(X,FMT,APOSTROPHES)

! ------------- ARGUMENTS-------------------------------------------------------

   CHARACTER*(*),INTENT(IN )          :: X           ! 
   CHARACTER*(*),INTENT(OUT)          :: FMT         ! format specifier to write X
   LOGICAL      ,INTENT(IN ),OPTIONAL :: APOSTROPHES !

!*------------- LOCAL DATA -----------------------------------------------------

   LOGICAL     :: BLANKS_INCLUDED,LOC_APOSTROPHES
 

	   IF (PRESENT(APOSTROPHES)) THEN
	     LOC_APOSTROPHES=APOSTROPHES
	   ELSE
	     LOC_APOSTROPHES=.FALSE.
	   END IF

 
       BLANKS_INCLUDED= INDEX(TRIM(X),' ')/=0
       
       IF (LOC_APOSTROPHES .AND. BLANKS_INCLUDED) THEN  
         WRITE(FMT,301)LEN_TRIM(X) 
       ELSE
         WRITE(FMT,300)LEN_TRIM(X) 
       ENDIF  
      
 300   FORMAT('(A',I3.3,')')
 301   FORMAT('('''''''',A',I3.3,','''''''')')

 END SUBROUTINE  GETFMT_CHAR_SCALAR


!+******************************************************************************
! GETFMT_CHAR_VECTOR 
!
! argument MODE can take one of the following values:
!
! 'INDIV'  each element of X has it individual edit descriptor 
! 'COMMON' all elements of X have the same edit descriptor (default)
! 'MAX'    same edit descriptor as for 'COMMON', but without repeat specification
!
! The character width is the length of X, without counting trailing blanks.  
!
! The optional argument APOSTROPHES has an effect only when MODE equals 'INDIV'.
! In this case, if the i'th element of X contains blanks, it will will be enclosed 
! by apostrophes when written with format FMT.
!
! The length LEN of FMT is either 6 (APOSTROPHES=false) or 16 (APOSTROPHES=true)
! The length LEN of FMT depends on MODE:  'MAX'    -> LEN=6 
!                                         'COMMON' -> LEN=12
!                                         'INDIV'  -> LEN at most 17*N-1
! 06-MAR-1999 GK created
!*******************************************************************************
  
 SUBROUTINE GETFMT_CHAR_VECTOR(X,FMT,APOSTROPHES,MODE)

! ------------- ARGUMENTS-------------------------------------------------------

   CHARACTER*(*),INTENT(IN )          :: X(:) !
   CHARACTER*(*),INTENT(OUT)          :: FMT  ! format specifier to write X
   LOGICAL      ,INTENT(IN ),OPTIONAL :: APOSTROPHES !
   CHARACTER*(*),INTENT(IN ),OPTIONAL :: MODE ! can be either 'INDIV' 'COMMON','MAX'

!*------------- LOCAL DATA -----------------------------------------------------

   INTEGER     :: WMAX,W,IR!
   INTEGER     :: N        ! size of X
   CHARACTER*14:: FMT0     ! elementary I edit descriptor. Ex.: 'A125'
   CHARACTER*6 :: LOC_MODE !
   LOGICAL     :: BLANKS_INCLUDED,LOC_APOSTROPHES
 

	   IF (PRESENT(APOSTROPHES)) THEN
	     LOC_APOSTROPHES=APOSTROPHES
	   ELSE
	     LOC_APOSTROPHES=.FALSE.
	   END IF 

	   IF (PRESENT(MODE)) THEN
	     LOC_MODE=MODE
	   ELSE
	     LOC_MODE='COMMON'
	   END IF
    
	   N = SIZE(X)
       WMAX=0
       FMT=''

       DO IR=1,N
	     W = LEN_TRIM(X(IR))  
         WMAX = MAX (WMAX,W)
         IF (LOC_MODE=='INDIV')  THEN
           BLANKS_INCLUDED= INDEX(TRIM(X(IR)),' ')/=0
           IF (LOC_APOSTROPHES .AND. BLANKS_INCLUDED) THEN  
             WRITE(FMT0,301)W 
           ELSE
             WRITE(FMT0,300)W 
           ENDIF
           FMT=TRIM(FMT)//FMT0
           IF (IR/=N)  FMT=TRIM(FMT)//',X,'
         ENDIF
       ENDDO
       
       IF (LOC_MODE=='INDIV') THEN  
         FMT='('//TRIM(FMT)//')'
       ELSEIF(LOC_MODE=='COMMON') THEN
         WRITE(FMT,302) N,WMAX
       ELSEIF (LOC_MODE=='MAX') THEN  
         WRITE(FMT,303)   WMAX
       ELSE
         WRITE(LU_FMTLIB,*)'GETFMT_INT__VECTOR: UNKNOWN VALUE FOR ARGUMENT MODE: ',LOC_MODE
         STOP
       ENDIF 

 300   FORMAT('A',I3.3)
 301   FORMAT(''''''''',A',I3.3,',''''''''')
 302   FORMAT('(',i2.2,'(A',i3.3,',X))')
 303   FORMAT('(A',i3.3,')')
      
 END SUBROUTINE  GETFMT_CHAR_VECTOR


 END MODULE FMTLIB
