!+******************************************************************************
! SETUPLIB.F90
!
! Purpose: module for enhanced free format reading from a setup file.
!
! Idea: Instead of assuming a fixed sequence of variables to be read from a setup 
!       file, the user can arrange the variables in an order and a format he likes.
!
!
! The setup file is a free format ASCII-file with the following features:
!
! * Each element is identified by a keyword (=token) followed by one or more input 
!   data of either type real, integer, character or logical. These data are read 
!   by list-directed input (FORMAT=*)
!
! * Any number of spaces, comments and blank lines are allowed to structure the 
!   setup file. Comments are introduced by the exclamation sign '!' 
!   (= variable COMMENT_SETUPLIB) and continue to the end of the line. 
!   Lines longer than 200 (=constant LINE_LENGTH) characters will be truncated.
!
! * The setup file is subdivided into different sections (categories,subcategories). 
!   A nonblank character in the first column defines the beginning of a new category.  
!   The end of a category is defined by the first occurrence of a nonblank character 
!   in column one (=start of next category).     
!   Caution: As a consequence it follows that the first column of the setup file is 
!   for category keywords only, a token erroneously positioned in the first column is 
!   interpreted as new category.  
!  
! * Free order of the categories within the setup file. If multiple categories with 
!   the same name exist, the category that occurs first is considered.
!
! * Free order of the tokens within each category. If multiple tokens with the same 
!   name exist, the token that occurs first is considered.
!
! * Categories and tokens may be either mandatory or optional. If a mandatory 
!   category or a mandatory token is missing, the program aborts with an error message.
!
!     
!
! * On output the variable may have either:
!
!   a) the value read from the setup file (if the token is found and has a value mapped)
!
!   If the token is optional and not found 
!      or a null is mapped to the variable (see list-directed input) then:
!   b) a default value              (if such a default value is given)
!   c) the value as it had on input (if no default value was given)
!  
!
!   Notice that in the following examples all tokens end with a "=", 
!   but this is not mandatory.
!
!
! Example #1: read the following Kepler elements:
! ----------------------------------------------   
! the next 3 lines represent a possible setup-file:
!
!   A= 42166770    ECC=0.000340014  INC=0.062872   ! enter any comment
!   LON=9.3434                                     ! not needed here
!   M= 249.897611  AOP=36.413842    RAN=6.044634  
!
! you can:
!
! * freely choose the order of the elements in the setup file (automatic)
!
! * allow any number of blanks and linefeeds (automatic)
!
! * enter comments in setup file (comment character is '!' by default) 
!
! * use simple syntax in the program: 
!   -> CALL READ_SETUP('A=',sma)
!   the "=" character behind each name of a variable is only a suggestion, 
!   not a must. The name of the FORTRAN variable (sma) and the token in the
!   setup-file ('A=') need not to be identical.
!
! * avoid mixing up different variables ( an additional blank must precede 
!   each token (X=7Y=8 will raise an error) and each token should end 
!   with a "=".)
!
! * ignore other items in the setup-file.  (here longitude LON) 
! 
! * expect meaningful error messages 
!   Example : "category KEPLER, token AOP is missing"
!             "category MANEUVER, token AOP not all 3 elements correcly filled"
!
!
! Example #2: read the following propagator model information:
! -------------------------------------------------------------
!
!  GM_SUN=1.327D11          ! not present: GM_MOON 
!  POTENTIAL_ORDER= 7
!
! you can:
!
! * allow variables to be optional (optional keyword: MANDATORY=.FALSE) and 
!   assign default values for that case (optional keyword: DEFAULT=):
!   -> CALL READ_SETUP('GM_MOON=',gm_moon,MANDATORY=.FALSE.,   &
!                      DEFAULT=4902.799d0,FOUND=exist)
!   an optional output variable (keyword FOUND=) tells you IF the token was
!   found. 
!   IF the token (here 'GM_MOON=') was found, a value must be given, or the
!   program will abort. (GM_MOON=?   -> abort)
!
!
!
! Example #3: read maneuver information:
! ----------------------------------------------
!
!MANEUVER
!   MAN_ID=EW112 BCT=21-NOV-1997:05:15:00  DV=0.01389  0.10354 0.00000    
!   MAN_ID=WE112 BCT=21-NOV-1997:16:15:00  DV=0.00888 -0.07837 0.00000
!
! you can:
!                                        
! * read REAL,INTEGER,CHARACTER*(*) and LOGICAL typed variables:
!   the common subroutine name is READ (it is a generic procedure, the compiler
!    decided which specific subroutine will be called).
!
! * read arrays 
!   -> CALL READ_SETUP('DV=',DELTAV(1:3))
!
! * group total setup file into smaller parts (called here a CATEGORY) 
!   a new category is assumed IF the first character of a line is nonblank:
!   -> CALL CATEGORY('MANEUVER',MANDATORY=.FALSE.,FOUND=CONTINUE)   
!
! * divide a category into sub-categories (and even lower with keyword LEVEL)
!   -> CALL SUB_CATEGORY('MAN_ID',FOUND=FOUND)
!   
! * either start the search from the beginning  of the category (optional 
!   keyword: REWIND=.TRUE.) or continue after the last token read.
!   -> CALL READ_SETUP('MAN_ID=',man_id(N_MAN),MANDATORY=FALSE,.REWIND=.FALSE)
!
! the following default values for the optional arguments are set:
!
!
! ARGUMENT  I/O    DEFAULT VALUE
! ------------------------------------------------------------------------------
! REWIND    IN     .TRUE.
! DEFAULT   IN     none
! MANDATORY IN     .FALSE. if argument "DEFAULT" or argument "FOUND" present, 
!                  DEFAULT_MANDATORY_SETUPLIB else
! FOUND     OUT    none
!
!
! SUBROUTINES used: NONE
!
! Modifications:
! 25-NOV-1997 GK created
! 06-APR-1998 GK default for mandatory = false if keyword "FOUND" present
! 09-APR-1998 GK SUB_CATEGORY created
! 17-APR-1998 GK ISIZE_1 and ISIZE_2 initialized to ISIZE_0 in OPEN_SETUP
! 23-APR-1998 GK ISIZE_0,ISIZE_1, ISIZE_2 -> ISIZE(0:MXLEV) ;ISIZE -> ISIZE(LEV)
! 10-JUN-1998 GK migration to PC 
! 29-JUN-1998 GK additional argument MANDATORY= for subroutine SUB_CATEGORY 
! 04-AUG-1998 GK Kind paramater LF
! 04-SEP-1998 GK category keyword now included in category
! 10-SEP-1998 GK Tabs are replaced by single spaces.
! 05-OCT-1998 GK new subsubcategory (level=3)
! 07-OCT-1998 GK V01.14 new optional argument FILE to OPEN_SETUP
! 25-FEB-1999 GK V02.01 Major revision: added WRITE_SETUP,PARTOK,GETFMT
! 08-MAR-1999 GK V02.02 added real matrix part
! 05-JAN-2000 GK V02.03 introduced output unit LU_DEBUG_SETUPLIB
! 19-MAR-2000 GK V02.04 strong privacy,renamed PARTOK to IO_SETUP, 
!                       new global variable DEFAULT_MANDATORY_SETUPLIB
!                       moved getfmt section into dedicated module FMTLIB
! 14-APR-2000 GK V02.05 new: LIST_DIRECTED_IO_SETUPLIB, N_READ, ERRMSG
!*******************************************************************************

 MODULE SETUPLIB 

 USE FMTLIB

 IMPLICIT NONE 

 PRIVATE ! strong privacy
 
 CHARACTER*5 ,PARAMETER,PUBLIC :: SETUPLIB_REV  = '02.05'       !   
 CHARACTER*11,PARAMETER,PUBLIC :: SETUPLIB_DATE = '14-APR-2000' ! 

 INTEGER,PARAMETER,PRIVATE :: LF = 8  ! Kind param. for 8-byte real numbers 
 
!*******************************************************************************
! Block#1
! Reading from setup file
!*******************************************************************************

 PUBLIC OPEN_SETUP     ! reads SETUP-file into RAM
 PUBLIC CATEGORY       ! restricts search to category
 PUBLIC SUB_CATEGORY   ! restricts search to sub_category
 PUBLIC READ_SETUP     ! Searchs token in [sub_]category and reads data behind
 PUBLIC SHOW_CONTENT   ! not documented, do only use in exceptional cases
 PUBLIC FOUND_TOKEN    ! returns true if the last token was found
 PUBLIC ELEMENTS_READ  ! returns the number of elements read
 PUBLIC FOUND_CATEGORY ! returns true if the last category was found

! public variables:

! if DEBUG_MESSAGES_SETUPLIB is true, debug informations 
! are written to unit LU_DEBUG_SETUPLIB

 LOGICAL,PUBLIC::DEBUG_MESSAGES_SETUPLIB = .FALSE. 

 INTEGER,PUBLIC::LU_DEBUG_SETUPLIB = 6

 ! The following logical specifies whether a token is required when the 
 ! optional argument MANDATORY is missing in the call to READ_SETUP:

 LOGICAL,PUBLIC::DEFAULT_MANDATORY_SETUPLIB=.TRUE.
 
 ! Severity level when optional token found, but not filled correctly 
 ! can be either [M]ute, [W]arning, [E]rror or [F]atal

 CHARACTER,PUBLIC:: SEV_OPT_IO_ERROR_SETUPLIB = 'W'
 
 ! if the next logical is true, pure list-directed IO is implemented 

 LOGICAL,PUBLIC::LIST_DIRECTED_IO_SETUPLIB =.FALSE.

 
 ! Comment character; all characters behind this character are ignored:

 CHARACTER,PUBLIC::COMMENT_SETUPLIB='!' 
   
 
! private constants for input SETUP file
! ------------------------------------------------------------------------------
 INTEGER,PARAMETER::MAX_CHARS=299999 ! maximum number of character to be read
                                    ! into memory (trailing blanks
                                    ! in each line are not count)
 INTEGER,PARAMETER::LINE_LENGTH=200 ! maximum SETUP-file line size
 INTEGER,PARAMETER::MXLEV   = 3     ! max. number of levels (file,category,...)

! private variables for input SETUP file
 
 CHARACTER(LEN=MAX_CHARS) :: FSTR ! at begin, entire input SETUP-file
                                  ! is copied into RAM (string FSTR)  

 LOGICAL :: TOKEN_FOUND   =.FALSE.! set by call to READ_SETUP or SUBCATEGORY
                                  ! true if the last token was found
 LOGICAL :: CATEGORY_FOUND=.FALSE.! true if the last category was found     

 INTEGER :: LEV=0  ! level of actual section: 0 for entire setupfile
                   ! 1 for category ,2 for subcategory, 3 for sub-subcategory
 INTEGER :: OFFSET ! offset from ISTART (characters)
 INTEGER :: N_READ ! number of fields read 
  INTEGER :: IOS    ! last IOSTAT 

 CHARACTER*100::LABEL (0:MXLEV) ! label of different levels
 INTEGER      ::ISTART(0:MXLEV) ! start of section 
 INTEGER      ::ISIZE (0:MXLEV) ! size  of section (in characters)  

! special characters

 CHARACTER,PARAMETER :: LINE_FEED = CHAR(13) ! marks end_of_line  
 CHARACTER,PARAMETER :: TAB       = CHAR(09) ! will be replaced by blanks   
 CHARACTER,PARAMETER :: BLANK     = ' '      !
   

! The following INTERFACE command is used to be able to use a common procedure
! name, CALL READ_SETUP(..) for all possible data types: 

 INTERFACE READ_SETUP

      MODULE PROCEDURE READ_REAL_SCALAR, READ_REAL_VECTOR, READ_REAL_MATRIX,   &
                       READ_INT__SCALAR, READ_INT__VECTOR,                     &
                       READ_CHAR_SCALAR, READ_CHAR_VECTOR,                     &
                       READ_LOG__SCALAR, READ_LOG__VECTOR
 END INTERFACE   




!*******************************************************************************
! Block#2
! Writing to setup file
!*******************************************************************************

 PUBLIC OPEN_SETUP_OUT  ! assigns a logical unit and opens the new SETUP-file
 PUBLIC WRITE_SETUP     ! Writes token and its value into the new SETUP-file
 PUBLIC WRITE_CATEGORY  ! Writes category into the new SETUP-file

! public variables for SETUP_OUT file:
 
 INTEGER,PUBLIC :: COLCAT_OUT_SETUPLIB = 10 ! first .. columns reserved for 
                                            ! categories, must be at least 1
 INTEGER,PUBLIC :: MAXCOL_OUT_SETUPLIB = 80 ! maximum SETUP_OUT line size
 
 INTEGER,PUBLIC :: VERBOS_SETUPLIB     = 30 !

! private variables for SETUP_OUT file:
! ------------------------------------------------------------------------------
    
 INTEGER :: UNIT_OUT  ! file unit of SETUP_OUT
 INTEGER :: ICOL = 0  ! actual position of cursor in SETUP_OUT

! The following INTERFACE command is used to be able to use a common procedure
! name, WRITE_SETUP(..) for all possible data types: 


 INTERFACE WRITE_SETUP
      
      MODULE PROCEDURE WRITE_REAL_SCALAR, WRITE_REAL_VECTOR, WRITE_REAL_MATRIX,   &
                       WRITE_INT__SCALAR, WRITE_INT__VECTOR,    &
                       WRITE_CHAR_SCALAR, WRITE_CHAR_VECTOR,    &
                       WRITE_LOG__SCALAR, WRITE_LOG__VECTOR

 END INTERFACE


!*******************************************************************************
! Block#3
! combined Reading/Writing to setup file
!*******************************************************************************

 INTEGER,PUBLIC :: OPTIO_SETUPLIB = 1  !

 PUBLIC IO_CATEGORY !
 PUBLIC IO_SETUP    !

! The following INTERFACE command is used to be able to use a common procedure
! name, IO_SETUP(..) for all possible data types: 

INTERFACE IO_SETUP

      MODULE PROCEDURE IO_REAL_SCALAR, IO_REAL_VECTOR, IO_REAL_MATRIX,         &
                       IO_INT__SCALAR, IO_INT__VECTOR,                         &
                       IO_CHAR_SCALAR, IO_CHAR_VECTOR,                         &
                       IO_LOG__SCALAR, IO_LOG__VECTOR

END INTERFACE








 CONTAINS

!+******************************************************************************
! IEND
! 15-NOV-1998 GK created
!*******************************************************************************
  INTEGER FUNCTION IEND (L)     
  ! ------------- ARGUMENTS-------------------------------------------------------

  INTEGER ,INTENT( IN) :: L  ! either 0 (SETUPFILE),1(CATEGORY),2(SUBCATEGORY),...
    IEND = ISTART(L)+ISIZE(L)-1
  END FUNCTION IEND

!+******************************************************************************
! TO_CHAR
!
! Purpose:   transforms an integer into a character string
!*******************************************************************************
  CHARACTER*13 FUNCTION TO_CHAR(I)
  
   INTEGER,INTENT(IN) :: I
  
     WRITE(TO_CHAR,*)I
     TO_CHAR=ADJUSTL(TO_CHAR)

  END FUNCTION TO_CHAR

!+******************************************************************************
! ERRMSG
!
! Purpose:   Displays error message and stops
!*******************************************************************************
  
  SUBROUTINE ERRMSG (SEV,MESSAGE)

!--------------------------- ARGUMENTS -------------------------------------

  CHARACTER*1,   INTENT(IN) :: SEV       ! 'M','W','E','F'
  CHARACTER*(*), INTENT(IN) :: MESSAGE   ! message
  
     
     IF      ( SEV  ==  'M') THEN
        CONTINUE 
     ELSE IF ( SEV  ==  'W') THEN
        WRITE (LU_DEBUG_SETUPLIB,'(A)')' Warning: '//TRIM(MESSAGE)
     ELSE IF ( SEV  ==  'E') THEN
        WRITE (LU_DEBUG_SETUPLIB,'(A)')' Error  : '//TRIM(MESSAGE)
     ELSE IF ( SEV  ==  'F') THEN
        WRITE (LU_DEBUG_SETUPLIB,'(A)')' Fatal  : '//TRIM(MESSAGE)
     END IF

     IF (SEV=='E' .OR. SEV=='F') STOP ' *** Program aborted ***'

  END SUBROUTINE ERRMSG

!+******************************************************************************
! SEVERITY
!
! Purpose:   returns severity for the case where a token is not correctly filled
!*******************************************************************************

   CHARACTER FUNCTION SEVERITY(MANDATORY)
     
    LOGICAL :: MANDATORY
    
     IF (MANDATORY) THEN
       SEVERITY='E'
     ELSE
       SEVERITY=SEV_OPT_IO_ERROR_SETUPLIB
     END IF
  
   END FUNCTION SEVERITY

!+******************************************************************************
! SPLIT
!
! Splits a string into substrings (separator = blank) with a special handling 
! of apostrophes and quotes.
!
! If a substring contains blanks it must be embedded into apostrophes or 
! quotation marks.  
!*******************************************************************************
 SUBROUTINE SPLIT(S_IN,M,S_OUT,N)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*),INTENT(IN ) :: S_IN
  INTEGER      ,INTENT(IN ) :: M
  INTEGER      ,INTENT(OUT) :: N
  CHARACTER*(*),INTENT(OUT) :: S_OUT(M)

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER                    :: I,I0,L            !
  CHARACTER,PARAMETER        :: APOSTROPHE=''''   !
  CHARACTER,PARAMETER        :: QUOTES    ='"'    !
  CHARACTER(LEN=LEN(S_IN)+1) :: S                 !

   N=0
   S=S_IN
   S_OUT(1:M)=''
   L = LEN_TRIM(S_IN)
   I=1

   DO WHILE (I<=L .AND. N<M) 
     DO WHILE (S(I:I)==BLANK .AND. I<=L)
       I=I+1
     END DO
     N=N+1
     I0=I
     IF (S(I:I)==APOSTROPHE) THEN
       I=I+1
       DO WHILE ((S(I:I)/=APOSTROPHE .OR. S(I+1:I+1)==APOSTROPHE) .AND. I<=L)
         IF (S(I:I)==APOSTROPHE .AND. S(I+1:I+1)==APOSTROPHE) S = S(:I)//S(I+2:)
         I=I+1
       END DO
       S_OUT(N)=S(I0+1:I-1)
       IF (I<=L) I=I+1 
     ELSEIF (S(I:I)==QUOTES) THEN
       I=I+1
       DO WHILE ((S(I:I)/=QUOTES .OR. S(I+1:I+1)==QUOTES) .AND. I<=L)
         IF (S(I:I)==QUOTES .AND. S(I+1:I+1)==QUOTES) S = S(:I)//S(I+2:)
         I=I+1
       END DO
       S_OUT(N)=S(I0+1:I-1)
       IF (I<=L) I=I+1
     ELSE      
       DO WHILE (S(I:I)/=BLANK .AND. I<=L)
         I=I+1
       END DO
       S_OUT(N)=S(I0:I-1)
     END IF
   END DO  
  END SUBROUTINE SPLIT

!+******************************************************************************
! SHOW_CONTENT
!
! 19-MAR-2000 GK created
!*******************************************************************************
  SUBROUTINE SHOW_CONTENT(LEV,STRING)

! ------------- ARGUMENTS-------------------------------------------------------
  
  INTEGER      ,INTENT(IN)  :: LEV
  CHARACTER*(*),INTENT(OUT) :: STRING

  STRING=FSTR(ISTART(LEV):IEND(LEV))
  END SUBROUTINE SHOW_CONTENT

!+******************************************************************************
! ELEMENTS_READ
!
! Purpose:   Returns nuber of elements found
!*******************************************************************************
  INTEGER FUNCTION ELEMENTS_READ()
    IF (LIST_DIRECTED_IO_SETUPLIB) THEN
      CALL ERRMSG('F','ELEMENTS_READ cannot be called for LIST_DIRECTED_IO')
    ELSE
      ELEMENTS_READ=N_READ
    ENDIF
  END FUNCTION ELEMENTS_READ

    
!+******************************************************************************
! FOUND_TOKEN
!
! Purpose:   Returns true if the the last token was found
!*******************************************************************************
  LOGICAL FUNCTION FOUND_TOKEN()
    FOUND_TOKEN=TOKEN_FOUND
  END FUNCTION FOUND_TOKEN

!+******************************************************************************
! FOUND_CATEGORY
!
! Purpose:  Returns true if the the last category was found
!*******************************************************************************
  LOGICAL FUNCTION  FOUND_CATEGORY() 
  FOUND_CATEGORY = CATEGORY_FOUND
  END FUNCTION  FOUND_CATEGORY

!+******************************************************************************
! OPEN_SETUP
!
!   reads SETUP-file into RAM.
!
! Arguments:
!
!   UNIT: logical unit    (OPTIONAL)
!         if not specified the next free unit starting from 10 is chosen.
!
!   FILE: SETUP file name (OPTIONAL)
!         if FILE is specified, an OPEN is performed.
!
!
! At least one of the two arguments UNIT and FILE should be specified. There are
! typically two different situation:
!
! * If the setup-file is not yet opened,  OPEN_SETUP should be called with
!   argument FILE only (recommended), the unit will be determined automatically,
!   the file will be opened, read and closed before leaving OPEN_SETUP.
!
! * The setup-file is already opened, and has unit LU assigned to it. In this
!   case OPEN_SETUP must be called with argument UNIT only. Be aware that no
!   rewind will be performed before reading the file. In addition the file
!   will not be closed before leaving OPEN_SETUP. (the file was opened outside
!   OPEN_SETUP and must be closed outside OPEN_SETUP)
!
!
! The following global variables are set:
!
!   *   ISTART(0..M)   := 1
!   *   ISIZE (0..M)   := size of entire setup-file
!   *   LEV            := 0 <=> whole setup
!   *   OFFSET         := 0
!   *   LABEL (0..M)   := ''
!   *   FSTR           := content of entire setup-file
!   *   OPTIO_SETUPLIB := 1
!
!
! The following parameters are used:
!
!   * maximum length of setup-file: MAX_CHARS
!   * maximum line length         : LINE_LENGTH
!
! External calls:
!
! 15-DEC-1997 GK created
! 10-SEP-1998 GK Tabs are replaced by single spaces.
!                use of intrinsic functions LEN_TRIM,TRIM
! 07-OCT-1998 GK new optional argument:FILE=, more exceptions 											
!*******************************************************************************

  SUBROUTINE OPEN_SETUP(UNIT,FILE)

! ------------- ARGUMENTS-------------------------------------------------------

  INTEGER       ,INTENT( IN),OPTIONAL :: UNIT   ! file handle (5 = standard input)
  CHARACTER*(*) ,INTENT( IN),OPTIONAL :: FILE   ! file name

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER :: LENGTH   !
  INTEGER :: SIZE     !
  INTEGER :: I        !
  INTEGER :: LOC_UNIT !
  LOGICAL :: IS_OPEN  ! true if unit already opened

  CHARACTER(LEN=LINE_LENGTH):: LINE ! content of one line of setup file


   IF (PRESENT(UNIT))   THEN
     LOC_UNIT = UNIT
   ELSE
     DO LOC_UNIT= 10,100
       INQUIRE (UNIT = LOC_UNIT, OPENED = IS_OPEN,ERR=900)
	   IF (.NOT. IS_OPEN) EXIT
	   IF (LOC_UNIT == 100) GOTO 910
     END DO
   END IF

   IF (PRESENT(FILE))   THEN
	 OPEN(LOC_UNIT,FILE=FILE,ACTION='READ',STATUS='OLD',ERR=920)
   END IF

   SIZE=1

   DO
      IF  (SIZE > MAX_CHARS - LINE_LENGTH -4)  GOTO 930

! --  read single line
      READ(LOC_UNIT,'(A)',END=20,ERR=940) LINE

! --  Ignore all behind  comment CHARACTER
      I = INDEX(LINE,COMMENT_SETUPLIB)
      IF (I /= 0) LINE = LINE(1:I-1)

! --  Ignore all trailing blanks
      LENGTH = LEN_TRIM(LINE)

! --  Replace tabs by blanks
      DO I=1,LENGTH
        IF (LINE(I:I)==TAB)  LINE(I:I) = BLANK
      END DO

! --  identify new category
      IF (LINE(1:1) /= BLANK) THEN
        FSTR(SIZE:SIZE+1) = LINE_FEED//BLANK
        SIZE=SIZE+2
      END IF

! --  append content of line
      FSTR(SIZE:SIZE+LENGTH)=TRIM(LINE)//BLANK
      SIZE=SIZE+LENGTH+1

   END DO

20 FSTR(SIZE:SIZE) = LINE_FEED        ! LAST CHARACTER

   LEV            = 0
   ISTART(0:MXLEV)= 1
   ISIZE (0:MXLEV)= SIZE
   LABEL (0:MXLEV)= ''
   OFFSET         = 0

! read keyword DEBUG_SETUPLIB=

   CALL READ_SETUP('DEBUG_SETUPLIB=',DEBUG_MESSAGES_SETUPLIB,                  &
                    MANDATORY=.FALSE.,DEFAULT=.FALSE.)

   OFFSET = 0

   IF (PRESENT(FILE))   THEN
  	  CLOSE(LOC_UNIT)
   END IF 	
	
   RETURN

! . EXCEPTIONS

900 CALL ERRMSG('E','OPEN_SETUP: Error during INQUIRE statement')
910 CALL ERRMSG('E','OPEN_SETUP: No free unit found')
920 CALL ERRMSG('E','OPEN_SETUP: Error opening file')
930 CALL ERRMSG('E','OPEN_SETUP: SETUP-file too large, increase MAX_CHARS')
940 CALL ERRMSG('E','OPEN_SETUP: Error reading SETUP file ')

  END SUBROUTINE OPEN_SETUP


!+******************************************************************************
! CATEGORY
!
! Finds a begin and size of a category
!
! following global variables are set:
!
!   ISTART(0)   not modified
!   ISTART(1)   := first character after CATEG (IEND(0) IF category not found)
!   ISTART(2:M) := ISTART(1)
!
!   ISIZE(0)    not modified
!   ISIZE(1)    := size of category (0 IF category not found)
!   ISIZE(2:M)  := 0
!
!   LEV         := 1 <=> category
!   OFFSET      := 0
!   LABEL(1)    := CATEG
!
! External calls: NONE
!
! 15-DEC-1997 GK created
! 06-APR-1998 GK IF keyword "FOUND" present => default for mandatory is .FALSE.
! 08-APR-1998 GK ISTART now corRESponds to first character after CATEG.
! 04-SEP-1998 GK category keyword itseft now included in category
!             -> this corrects the bug that there were needed at least 2 spaced
!                between category keyword and the first subcategory keyword.!
! 05-OCT-1998 GK Undo the following modification from 04-SEP-1998:
!                to be sure that the category itself cannot be recognized as token
!                the two characters preceeding the category LINE_FEED,BLANK have been
!                replaced by a LINE_FEED. As a blank is always expected in front of a
!                token and as now a LINE_FEED preceeds the category, there is no
!                risk of mixing up the category and the token.
!*******************************************************************************

  SUBROUTINE CATEGORY(CATEG,MANDATORY,FOUND)

! --------------------- ARGUMENTS----------------------------------------------

  CHARACTER*(*),   INTENT( IN) :: CATEG      ! name of category
  LOGICAL,OPTIONAL,INTENT( IN) :: MANDATORY  ! aborts if unmatched category
  LOGICAL,OPTIONAL,INTENT(OUT) :: FOUND      ! true if category found

!*--------------------- LOCAL DATA --------------------------------------------

  LOGICAL :: LOC_MANDATORY        !
  INTEGER :: CTMP                 !



    LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND. .NOT. PRESENT(FOUND)
    IF ( PRESENT(MANDATORY) )  LOC_MANDATORY = MANDATORY

! --find the category
    CTMP = INDEX(FSTR(ISTART(0):IEND(0)),LINE_FEED//BLANK//CATEG)
    CATEGORY_FOUND = (CTMP /= 0)


    IF (CATEGORY_FOUND) THEN
! ---  Found a category
       ISTART(1) = 2 + CTMP ! + LEN(CATEG)

! ---  determine the size -> begin of next category = end of actual category
       ISIZE (1) = INDEX( FSTR(ISTART(1):IEND(0)),LINE_FEED ) -1
    ELSE
       ISTART(1) = IEND(0)
       ISIZE (1) = 0
    END IF

    IF (LOC_MANDATORY .AND. .NOT.CATEGORY_FOUND) THEN
      CALL ERRMSG('E','SETUPLIB: unmatched category: '//TRIM(CATEG))
    END IF

    IF (PRESENT(FOUND)) FOUND=CATEGORY_FOUND

    OFFSET          = 0
    LEV             = 1
    LABEL (1)       = CATEG
    ISIZE (2:MXLEV) = 0
    ISTART(2:MXLEV) = ISTART(1)
    
    IF (DEBUG_MESSAGES_SETUPLIB) WRITE(LU_DEBUG_SETUPLIB,*)CATEG

   END SUBROUTINE CATEGORY

!+******************************************************************************
! SUB_CATEGORY
!
! Find starting point and size of a sub_category
!
! The following global variables are set:
!
!      ISTART(0)      not modified
!      ISTART(1..L-1) not modified
!      ISTART(L)      := first character after SUBCATEG (IEND(1) IF not found)
!      ISTART(L+1:M)  := ISTART(L)
!
!      ISIZE(0)       not modified
!      ISIZE(1:L-1)   not modified
!      ISIZE(L)       := size of subcategory (0 IF not found)
!      ISIZE(L+1:M)   := 0
!
!      OFFSET         := 0
!      LEV            := 2
!      LABEL(L)       := SUBCATEG
!
! External calls: NONE
!
! 08-APR-1998 GK created
!*******************************************************************************

  SUBROUTINE SUB_CATEGORY(SUBCATEG,MANDATORY,LEVEL,FOUND)


! --------------------- ARGUMENTS----------------------------------------------

  CHARACTER*(*),   INTENT( IN) :: SUBCATEG   ! name of subcategory
  LOGICAL,OPTIONAL,INTENT( IN) :: MANDATORY  ! aborts if unmatched subcategory
  INTEGER,OPTIONAL,INTENT( IN) :: LEVEL      ! level: 2 for subcategory    (default)
!        3 for subsubcategory
  LOGICAL,OPTIONAL,INTENT(OUT) :: FOUND      ! true if subcategory found

!*--------------------- LOCAL DATA --------------------------------------------

  LOGICAL :: LOC_MANDATORY        !
  INTEGER :: CTMP                 !

  IF ( PRESENT(LEVEL) )  THEN
    IF (LEVEL > MXLEV .OR. LEVEL < 2) THEN
      CALL ERRMSG('F','SETUPLIB: wrong LEVEL, must lie between 2 and MXLEV')
    END IF
    LEV = LEVEL
  ELSE
    LEV = 2
  END IF


  LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND. .NOT. PRESENT(FOUND)
  IF ( PRESENT(MANDATORY) )  LOC_MANDATORY = MANDATORY

! --- find a subcategory

   CTMP = INDEX(FSTR(ISTART(LEV)+1:IEND(LEV-1)),BLANK//SUBCATEG)
   TOKEN_FOUND = (CTMP /= 0)

   IF (TOKEN_FOUND) THEN
! ----Found a subcategory , determine the size
      ISTART   (LEV) = ISTART(LEV) + CTMP

! subcategory ends at next appearence of the same toke
      ISIZE    (LEV) = INDEX(FSTR(ISTART(LEV)+1:IEND(LEV-1)),BLANK//SUBCATEG) + 1

! IF same token was not found  subcategory ends at category end:
      IF (ISIZE(LEV) == 1) ISIZE(LEV) = ISTART(LEV-1)+ISIZE(LEV-1)-ISTART(LEV)
   ELSE
      ISTART   (LEV) = IEND(LEV-1)
      ISIZE    (LEV) = 0
   END IF


   IF (LOC_MANDATORY .AND. .NOT.TOKEN_FOUND) THEN
     CALL ERRMSG('E','SETUPLIB: unmatched sub_category: '//TRIM(SUBCATEG))
   END IF

   IF (PRESENT(FOUND)) FOUND=TOKEN_FOUND

   IF (LEV<MXLEV) ISIZE (1+LEV:MXLEV) = 0
   IF (LEV<MXLEV) ISTART(1+LEV:MXLEV) = ISTART(LEV)
   OFFSET    = 0
   LABEL(LEV)= SUBCATEG

  END SUBROUTINE SUB_CATEGORY


!------------------------------------------------------------------------------
!
! SEARCH_TOKEN
!
! Purpose:
!
!   Returns true if the the last category was found
!
!---------------------------------------------------------------------------

SUBROUTINE SEARCH_TOKEN(TOKEN,MANDATORY,REWIND,FOUND) 

  CHARACTER*(*)   ,INTENT(IN) ::TOKEN     ! string, identifies item in SETUP
  LOGICAL         ,INTENT(IN) ::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL         ,INTENT(IN) ::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL,OPTIONAL,INTENT(OUT)::FOUND     ! true if token found

  INTEGER :: I,I0 
  
! --rewind to category begin
   IF (REWIND) OFFSET = 0
   I0 = ISTART(LEV)+OFFSET

! --Seek TOKEN within a category (each token must be preceded by a blank)
   I =INDEX(FSTR(I0:IEND(LEV)),BLANK//TOKEN)
   TOKEN_FOUND = ( I/=0 )
   IF (PRESENT(FOUND)) FOUND = TOKEN_FOUND

   IF (.NOT.TOKEN_FOUND) THEN
     IF(MANDATORY) CALL ERRMSG('E','READ_SETUP: category '//TRIM(LABEL(1))//   &
                               ', token '//TRIM(TOKEN)//' is missing')
   ELSE
! -- TOKEN has been found, try reading next element
     OFFSET = OFFSET + I + LEN(TOKEN)
   ENDIF

END SUBROUTINE SEARCH_TOKEN



!+******************************************************************************
! READ_REAL_VECTOR, READ_INT__VECTOR, READ_CHAR_VECTOR, READ_LOG__VECTOR
! READ_REAL_SCALAR, READ_INT__SCALAR, READ_CHAR_SCALAR, READ_LOG__SCALAR
!
! grouped together into generic procedure READ_SETUP
!
!+******************************************************************************
!
! Searchs token TOKEN and reads an array or scalar of either real, integer,logical
! or character RES located behind TOKEN.
!
! The read of the data RES is performed as list-directed input:
! For arrays, each element can be either separated by a comma ',' or a blank.
! each elements  be either a value or a null. (Ex. N=1,,3  the second element
! will not be changed).  A slash '/' terminates the input stream.
!
! If the optional argument MANDATORY is true, the token must be present, if it is
! false and the token is not found,  RES will be returned with the default
! value DEFAULT. The optional argument FOUND contains the information if TOKEN has
! been found.
!
! The search of TOKEN starts either from the beginning of the section, or continues
! after the last item read if the optional argument REWIND is present and false.
!
!
! The following global variables are used:
!
!  LEV                      level of section (entire file, category, subcateg, ...)
!  ISTART(LEV)              starting position of section
!  ISIZE(LEV)               size of section
!  LABEL(1)                 category name
!  DEBUG_MESSAGES_SETUPLIB  if true some debug messages are displayed
!
! The following global variables are modified
!
!  OFFSET   at end of READ_SETUP will be the position behind last item read
!
! External calls: NONE
!
! 15-DEC-1997 GK created
! 06-APR-1998 GK default for mandatory = false if keyword "FOUND" present
!*******************************************************************************

! ------------------------------------------------------------------------------

  SUBROUTINE READ_REAL_SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  REAL(KIND=LF)         ,INTENT(INOUT)::RES       ! scalar into which data is read
  REAL(KIND=LF),OPTIONAL,INTENT(IN   )::DEFAULT   ! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*--------- LOCAL DATA ---------------------------------------------------------

  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(1),ERRTOK            !
  CHARACTER     :: SEV                       ! 

   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY) 
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       SRES(1)=FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(SRES(1),*,IOSTAT=IOS)RES
      IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//'not correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),1,SRES,N_READ)
       IF (N_READ==0) CALL ERRMSG(SEV,TRIM(ERRTOK)//'empty')
       READ(SRES(1),*,IOSTAT=IOS) RES
       IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//TRIM(SRES(1))//' cannot be read')
     END IF
   END IF
  
   IF (DEBUG_MESSAGES_SETUPLIB) WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES

  END SUBROUTINE READ_REAL_SCALAR

! -----------------------------------------------------------------------------

  SUBROUTINE READ_INT__SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  INTEGER               ,INTENT(INOUT)::RES       ! scalar into which data is read
  INTEGER      ,OPTIONAL,INTENT(IN   )::DEFAULT   ! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*--------- LOCAL DATA ---------------------------------------------------------

  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(1),ERRTOK            !
  CHARACTER     :: SEV                       ! 

   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY) 
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       SRES(1)=FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(SRES(1),*,IOSTAT=IOS)RES
      IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//'not correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),1,SRES,N_READ)
       IF (N_READ==0) CALL ERRMSG(SEV,TRIM(ERRTOK)//'empty')
       READ(SRES(1),*,IOSTAT=IOS) RES
       IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//TRIM(SRES(1))//' cannot be read')
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES

  END SUBROUTINE READ_INT__SCALAR

! ------------------------------------------------------------------------------

  SUBROUTINE READ_LOG__SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  LOGICAL               ,INTENT(INOUT)::RES       ! scalar into which data is read
  LOGICAL      ,OPTIONAL,INTENT(IN   )::DEFAULT   ! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*--------- LOCAL DATA ---------------------------------------------------------

  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(1),ERRTOK            !
  CHARACTER     :: SEV                       ! 

   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY) 
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       SRES(1)=FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(SRES(1),*,IOSTAT=IOS)RES
      IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//'not correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),1,SRES,N_READ)
       IF (N_READ==0) CALL ERRMSG(SEV,TRIM(ERRTOK)//'empty')
       READ(SRES(1),*,IOSTAT=IOS) RES
       IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//TRIM(SRES(1))//' cannot be read')
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES

  END SUBROUTINE READ_LOG__SCALAR

! ------------------------------------------------------------------------------

  SUBROUTINE READ_CHAR_SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  CHARACTER*(*)         ,INTENT(INOUT)::RES       ! scalar into which data is read
  CHARACTER*(*),OPTIONAL,INTENT(IN   )::DEFAULT   ! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*--------- LOCAL DATA ---------------------------------------------------------

  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(1),ERRTOK            !
  CHARACTER     :: SEV                       ! 

   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY) 
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       SRES(1)=FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(SRES(1),*,IOSTAT=IOS)RES
      IF(IOS/=0)CALL ERRMSG(SEV,TRIM(ERRTOK)//'not correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),1,SRES,N_READ)
       IF (N_READ==0) CALL ERRMSG(SEV,TRIM(ERRTOK)//'empty')
       RES=SRES(1)
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES

  END SUBROUTINE READ_CHAR_SCALAR

! ------------------------------------------------------------------------------

  SUBROUTINE READ_REAL_VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  REAL(KIND=LF)         ,INTENT(INOUT)::RES    (:)! array into which data is read
  REAL(KIND=LF),OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER       :: I,N                       !
  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(SIZE(RES)),ERRTOK    !
  CHARACTER*1024:: READ_STRING               !
  CHARACTER     :: SEV                       ! 

   N             = SIZE(RES)
   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY)
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       READ_STRING = FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(READ_STRING,*,IOSTAT=IOS) RES
       IF (IOS/=0) CALL ERRMSG(SEV,TRIM(ERRTOK)//' not all '//     &
                               TRIM(TO_CHAR(N))//' elements correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),N,SRES,N_READ)
       IF (N_READ<N) CALL ERRMSG(SEV,TRIM(ERRTOK)//', only '//     &
                                 TRIM(TO_CHAR(N_READ))//' from '//           &
                                 TRIM(TO_CHAR(N     ))//' elements found')
       DO I=1,N_READ
         READ(SRES(I),*,IOSTAT=IOS)RES(I)
         IF (IOS/=0) THEN
           CALL ERRMSG(SEV,TRIM(ERRTOK)//', cannot read element#'//&
                       TRIM(TO_CHAR(I))//': '//TRIM(SRES(I)))
           N_READ=I-1
           EXIT
         END IF
       END DO
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) THEN
     WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES(1)
     DO I=2,N
       WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV+LEN_TRIM(TOKEN)),RES(I)
     ENDDO     
   END IF

  END SUBROUTINE READ_REAL_VECTOR


! ------------------------------------------------------------------------------

  SUBROUTINE READ_INT__VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  INTEGER               ,INTENT(INOUT)::RES    (:)! array into which data is read
  INTEGER      ,OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER       :: I,N                       !
  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(SIZE(RES)),ERRTOK    !
  CHARACTER*1024:: READ_STRING               !
  CHARACTER     :: SEV                       ! 

   N             = SIZE(RES)
   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY)
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       READ_STRING = FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(READ_STRING,*,IOSTAT=IOS) RES
       IF (IOS/=0) CALL ERRMSG(SEV,TRIM(ERRTOK)//' not all '//     &
                               TRIM(TO_CHAR(N))//' elements correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),N,SRES,N_READ)
       IF (N_READ<N) CALL ERRMSG(SEV,TRIM(ERRTOK)//', only '//     &
                                 TRIM(TO_CHAR(N_READ))//' from '//           &
                                 TRIM(TO_CHAR(N     ))//' elements found')
       DO I=1,N_READ
         READ(SRES(I),*,IOSTAT=IOS)RES(I)
         IF (IOS/=0) THEN
           CALL ERRMSG(SEV,TRIM(ERRTOK)//', cannot read element#'//&
                       TRIM(TO_CHAR(I))//': '//TRIM(SRES(I)))
           N_READ=I-1
           EXIT
         END IF
       END DO
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) THEN
     WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES(1)
     DO I=2,N
       WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV+LEN_TRIM(TOKEN)),RES(I)
     ENDDO     
   END IF

  END SUBROUTINE READ_INT__VECTOR

! ------------------------------------------------------------------------------

  SUBROUTINE READ_LOG__VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  LOGICAL               ,INTENT(INOUT)::RES    (:)! array into which data is read
  LOGICAL      ,OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER       :: I,N                       !
  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: SRES(SIZE(RES)),ERRTOK    !
  CHARACTER*1024:: READ_STRING               !
  CHARACTER     :: SEV                       ! 

   N             = SIZE(RES)
   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY)
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       READ_STRING = FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(READ_STRING,*,IOSTAT=IOS) RES
       IF (IOS/=0) CALL ERRMSG(SEV,TRIM(ERRTOK)//' not all '//     &
                               TRIM(TO_CHAR(N))//' elements correcty filled')
     ELSE
       CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),N,SRES,N_READ)
       IF (N_READ<N) CALL ERRMSG(SEV,TRIM(ERRTOK)//', only '//     &
                                 TRIM(TO_CHAR(N_READ))//' from '//           &
                                 TRIM(TO_CHAR(N     ))//' elements found')
       DO I=1,N_READ
         READ(SRES(I),*,IOSTAT=IOS)RES(I)
         IF (IOS/=0) THEN
           CALL ERRMSG(SEV,TRIM(ERRTOK)//', cannot read element#'//&
                       TRIM(TO_CHAR(I))//': '//TRIM(SRES(I)))
           N_READ=I-1
           EXIT
         END IF
       END DO
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) THEN
     WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES(1)
     DO I=2,N
       WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV+LEN_TRIM(TOKEN)),RES(I)
     ENDDO     
   END IF

  END SUBROUTINE READ_LOG__VECTOR

! -----------------------------------------------------------------------------

  SUBROUTINE READ_CHAR_VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  CHARACTER*(*)         ,INTENT(INOUT)::RES    (:)! array into which data is read
  CHARACTER*(*),OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER       :: I,N                       !
  LOGICAL       :: LOC_MANDATORY,LOC_REWIND  !
  CHARACTER*256 :: ERRTOK                    !
  CHARACTER*1024:: READ_STRING               !
  CHARACTER     :: SEV                       ! 

   N             = SIZE(RES)
   LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND.                            &
                  .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
   IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY
   LOC_REWIND    =.TRUE.
   IF (PRESENT(REWIND   )) LOC_REWIND    = REWIND
   
   IF (PRESENT(DEFAULT  )) RES = DEFAULT ! filling with default values

   CALL SEARCH_TOKEN(TOKEN,LOC_MANDATORY,LOC_REWIND,FOUND)

   IF (TOKEN_FOUND)THEN ! token has been found, try reading next element
     SEV=SEVERITY(LOC_MANDATORY)
     ERRTOK='READ_SETUP: category '//TRIM(LABEL(1))//', token '//TRIM(TOKEN)//' '
     IF (LIST_DIRECTED_IO_SETUPLIB)THEN
       READ_STRING = FSTR(ISTART(LEV)+OFFSET:IEND(LEV))
       READ(READ_STRING,*,IOSTAT=IOS) RES
       IF (IOS/=0) CALL ERRMSG(SEV,TRIM(ERRTOK)//' not all '//     &
                               TRIM(TO_CHAR(N))//' elements correcty filled')
     ELSE
         CALL SPLIT(FSTR(ISTART(LEV)+OFFSET:IEND(LEV)),N,RES,N_READ)
         IF (N_READ<N) CALL ERRMSG(SEV,TRIM(ERRTOK)//', only '//     &
                                   TRIM(TO_CHAR(N_READ))//' from '//           &
                                   TRIM(TO_CHAR(N     ))//' elements found')
     END IF
   END IF
   
   IF (DEBUG_MESSAGES_SETUPLIB) THEN
     WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV)//TRIM(TOKEN),RES(1)
     DO I=2,N
       WRITE(LU_DEBUG_SETUPLIB,*)REPEAT(' ',9*LEV+LEN_TRIM(TOKEN)),RES(I)
     ENDDO     
   END IF

  END SUBROUTINE READ_CHAR_VECTOR

! ------------------------------------------------------------------------------

  SUBROUTINE READ_REAL_MATRIX(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN       ! string, identifies item in SETUP
  REAL(KIND=LF)         ,INTENT(INOUT)::RES    (:,:)! array into which data is read
  REAL(KIND=LF),OPTIONAL,INTENT(IN   )::DEFAULT(:,:)! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY   ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND      ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND       ! true if TOKEN found

!*------------- LOCAL DATA -----------------------------------------------------

   REAL(KIND=LF) :: RESV(SIZE(RES))

      RESV=RESHAPE(RES,(/SIZE(RES)/))

      IF (PRESENT(DEFAULT)) THEN
        CALL READ_REAL_VECTOR(TOKEN,RESV,MANDATORY,                            &
                              RESHAPE(DEFAULT,(/SIZE(DEFAULT)/)),              &
                              REWIND=REWIND,FOUND=FOUND)
      ELSE
        CALL READ_REAL_VECTOR(TOKEN,RESV,MANDATORY,REWIND=REWIND,FOUND=FOUND)
      END IF

      RES=RESHAPE(RESV,(/SIZE(RES,1),SIZE(RES,2)/))

  END SUBROUTINE READ_REAL_MATRIX






!*******************************************************************************
! output to new setup file section
!*******************************************************************************


!+******************************************************************************
! OPEN_SETUP_OUT
!
! assigns a logical unit and opens the new SETUP-file 
!
! UNIT: logical unit    (OPTIONAL)    
!       if not specified the next free unit starting from 10 is chosen.
!
! FILE: SETUP file name (OPTIONAL)
!       only if FILE is specified, an OPEN is performed.
!
! The following global variables are set: 
!
! *  UNIT_OUT       := logical unit of new SETUP-file 
! *  OPTIO_SETUPLIB := 2
!
! External calls: NONE
!
! 27-FEB-1999 GK created
!*******************************************************************************
    
  SUBROUTINE OPEN_SETUP_OUT(UNIT,FILE)

! ------------- ARGUMENTS-------------------------------------------------------

  INTEGER       ,INTENT( IN),OPTIONAL :: UNIT   ! file handle (6 = standard output) 
  CHARACTER*(*) ,INTENT( IN),OPTIONAL :: FILE   ! file name

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL              :: IS_OPEN

  IF (PRESENT(UNIT))   THEN
    UNIT_OUT = UNIT 
  ELSE
    DO UNIT_OUT= 10,100
       INQUIRE (UNIT = UNIT_OUT, OPENED = IS_OPEN,ERR=900)
	   IF (.NOT. IS_OPEN) EXIT
	   IF (UNIT_OUT == 100) GOTO 910
    END DO
  END IF
  
  IF (PRESENT(FILE))   THEN
	OPEN(UNIT_OUT,FILE=FILE,ACTION='WRITE',ERR=920)
  END IF
  
  OPTIO_SETUPLIB = 2
   
  RETURN

! EXCEPTIONS
900 CALL ERRMSG('E','OPEN_SETUP_OUT: Error during INQUIRE statement')
910 CALL ERRMSG('E','OPEN_SETUP_OUT: No free unit found')
920 CALL ERRMSG('E','OPEN_SETUP_OUT: Error opening SETUP_OUT file ')

  END SUBROUTINE OPEN_SETUP_OUT


!+******************************************************************************
! WRITE_TOKEN 
!
! writes TOKEN into the new SETUP-file
!
! if optional argument COLUMN is present, the output should be written at position 
! COLUMN (counted from MAXCOL_OUT_SETUPLIB) else behind last output
!
! 26-FEB-1998 GK created
!*******************************************************************************

 SUBROUTINE WRITE_TOKEN(TOKEN,COLUMN)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN  ! string, identifies item in SETUP
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN ! output should be written at this position 

!*------------- LOCAL DATA -----------------------------------------------------

   INTEGER             :: ICOL_NEW   ! column at which output will be written
   INTEGER             :: LENTOK     ! length of token
       
     IF (PRESENT(COLUMN)) THEN  ! ----- Put cursor at position COLUMN
       ICOL_NEW = MAX(COLCAT_OUT_SETUPLIB + COLUMN - 1, 0)
       IF (ICOL_NEW > ICOL) THEN
         WRITE(UNIT_OUT,'(a,$)')REPEAT(BLANK,ICOL_NEW-ICOL)   
       ELSEIF (ICOL_NEW < ICOL) THEN 
         WRITE(UNIT_OUT,'(/,a,$)')REPEAT(BLANK,ICOL_NEW)
       END IF
       ICOL=ICOL_NEW
     END IF
       
! - write TOKEN at position ICOL or into next line if not enough space
     LENTOK=LEN_TRIM(TOKEN)                  

     IF (ICOL + LENTOK + 1 <= MAXCOL_OUT_SETUPLIB) THEN       
       WRITE(UNIT_OUT,'(a,$)') TRIM(TOKEN)  
       ICOL = ICOL  + LENTOK
     ELSE
       WRITE(UNIT_OUT,'(/,2a,$)')REPEAT(BLANK,COLCAT_OUT_SETUPLIB),TRIM(TOKEN)
       ICOL = COLCAT_OUT_SETUPLIB + LENTOK
     END IF

  END SUBROUTINE WRITE_TOKEN 
      
!+******************************************************************************
! WRITE_CATEGORY
!
! writes CATEG at beginning of a new line
!
! global variables set: ICOL := column of last character written
!
! 27-FEB-1998 GK created
!*******************************************************************************

 SUBROUTINE WRITE_CATEGORY(CATEG)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*),   INTENT( IN) :: CATEG      ! name of category 

!*------------- LOCAL DATA -----------------------------------------------------

  INTEGER             :: LENCAT  ! length of string CATEG
  INTEGER             :: NBLANKS ! number of blanks behind category
       
    LENCAT = LEN_TRIM(CATEG)
    NBLANKS=COLCAT_OUT_SETUPLIB - LENCAT

    IF (ICOL/=0)  WRITE(UNIT_OUT,'(1x)')  ! -> newline 
                                            
    IF (NBLANKS<1) THEN
      WRITE(UNIT_OUT,'(a,x,$)')          TRIM(CATEG)
      ICOL = LENCAT+1
    ELSE
      WRITE(UNIT_OUT,'(a,a,$)') TRIM(CATEG),REPEAT(BLANK,NBLANKS)
      ICOL = COLCAT_OUT_SETUPLIB
    END IF
    
END SUBROUTINE WRITE_CATEGORY   

!+******************************************************************************
! WRITE_REAL_VECTOR, WRITE_INT__VECTOR, WRITE_CHAR_VECTOR, WRITE_LOG__VECTOR
! WRITE_REAL_SCALAR,WRITE_INT__SCALAR,WRITE_CHAR_SCALAR,WRITE_LOG__SCALAR
!
! grouped together into generic procedure WRITE_SETUP
!
! Writes TOKEN and its value RES into the new SETUP-file. 
!
! If the optional argument COLUMN is present, then the output is written at 
! column COLUMN (counted from MAXCOL_OUT_SETUPLIB), if it is not present, then the output is 
! written behind then last output. If there is not enough space the output is 
! written into the next line- 
!
! The output is written with the format descriptor FMT, if it is not present, 
! then it is automatically determined. For arrays the same descriptor is used for
! all elements.
!
! Examples of FMT: '(A12)', '(''X'',I2)', '(f9.1)', '(l1)'
!
! For real numbers it is possible to specify with the optional argument DIGITS
! the maximum number of decimal digits to be displayed in F edit description.
!
! In most cases it should be not nessecary to specify any optional keyword.
!
! global variables set: ICOL := column behind last character written
!
! parameters: max. output line size                      : MAXCOL_OUT_SETUPLIB
!             number of columns reserved for CATEGORIES  : COLCAT_OUT_SETUPLIB
!
! 25-FEB-1998 GK created
!*******************************************************************************


SUBROUTINE WRITE_REAL_SCALAR(TOKEN,RES,COLUMN,DIGITS,FMT)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN  ! string, identifies item in SETUP
  REAL(KIND=LF)          ,INTENT(IN) :: RES    ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN ! output should start at this column
  INTEGER      ,OPTIONAL ,INTENT(IN) :: DIGITS ! maximum number of decimal digits
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT    ! format descriptor for RES

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*20        :: LOC_FMT ! format descriptor 
   CHARACTER*99        :: CHR_RES ! RES written with format FMT
   INTEGER             :: LENRES  ! length of string CHR_RES
   INTEGER             :: I       ! dummy
       
!  write TOKEN at position COLUMN

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
     ELSE
         CALL GETFMT(RES,LOC_FMT,DIGITS=DIGITS)
     END IF

!  write result into CHR_RES

     WRITE(CHR_RES,LOC_FMT) RES     
     
     I=INDEX(CHR_RES,' .')  ; IF (I/=0) CHR_RES(I:I+1)='0.'
     I=INDEX(CHR_RES,' -.') ; IF (I/=0) CHR_RES(I:I+2)='-0.'
     I=LEN_TRIM(CHR_RES)    ; IF (CHR_RES(I:I)=='.') CHR_RES(I:I)=' ' 
 
!  write CHR_RES at column ICOL or into next line if not enough space

     LENRES=LEN_TRIM(CHR_RES)
     IF (ICOL + LENRES + 1 <= MAXCOL_OUT_SETUPLIB) THEN
       WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
     ELSE
       ICOL = COLCAT_OUT_SETUPLIB
       WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
     END IF
       
     ICOL = ICOL + LENRES + 1

 END SUBROUTINE WRITE_REAL_SCALAR

! ------------------------------------------------------------------------------

 SUBROUTINE WRITE_REAL_VECTOR(TOKEN,RES,COLUMN,DIGITS,FMT,LINEFEED)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN    ! string, identifies item in SETUP
  REAL(KIND=LF)          ,INTENT(IN) :: RES(:)   ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN   ! output should start at this column
  INTEGER      ,OPTIONAL ,INTENT(IN) :: DIGITS   ! maximum number of decimal digits
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT      ! format descriptor for RES
  INTEGER      ,OPTIONAL ,INTENT(IN) :: LINEFEED ! new line at least each LINEFEED elements of RES 

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*200       :: LOC_FMT    ! format descriptor 
   CHARACTER*99        :: CHR_RES    ! RES written with format FMT
   INTEGER             :: LENRES     ! length of string CHR_RES
   INTEGER             :: N          ! size of RES
   INTEGER             :: I          ! dummy
   INTEGER             :: IR         ! runs over range of RES
   INTEGER             :: ICOL_START ! column at which first value of RES is written 
   INTEGER             :: LOC_LF     ! new line at least each LINEFEED elements of RES     

     N = SIZE(RES)

!  write TOKEN at position COLUMN (column =COLCAT_OUT_SETUPLIB+COLUMN)

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

     IF (PRESENT(LINEFEED)) THEN
         LOC_LF = LINEFEED 
     ELSE
         LOC_LF = N
     END IF

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT  
     ELSE
         CALL GETFMT(RES,LOC_FMT,DIGITS=DIGITS,MODE='MAX') 
     END IF
     
     ICOL_START=ICOL
     
     DO IR=1,N

!  write result into CHR_RES

       WRITE(CHR_RES,LOC_FMT) RES(IR) ! same format for all elements
           
       I=INDEX(CHR_RES,' .')  ; IF (I/=0) CHR_RES(I:I+1)='0.'
       I=INDEX(CHR_RES,' -.') ; IF (I/=0) CHR_RES(I:I+2)='-0.'
       I=LEN_TRIM(CHR_RES)    ; IF (CHR_RES(I:I)=='.') CHR_RES(I:I)=' '
 
!  write CHR_RES at column ICOL or into next line if not enough space

       LENRES=LEN_TRIM(CHR_RES)

       IF ((ICOL + LENRES + 1 > MAXCOL_OUT_SETUPLIB) .OR. (MOD(IR,LOC_LF)==1 .AND. IR/=1)) THEN
         ICOL = ICOL_START
         WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
       ELSE
         WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
       END IF

       ICOL = ICOL + LENRES + 1
     
     END DO 

 END SUBROUTINE WRITE_REAL_VECTOR

! ------------------------------------------------------------------------------

 SUBROUTINE WRITE_REAL_MATRIX(TOKEN,RES,COLUMN,DIGITS,FMT,LINEFEED)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN    ! string, identifies item in SETUP
  REAL(KIND=LF)          ,INTENT(IN) :: RES(:,:) ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN   ! output should start at this column
  INTEGER      ,OPTIONAL ,INTENT(IN) :: DIGITS   ! maximum number of decimal digits
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT      ! format descriptor for RES
  INTEGER      ,OPTIONAL ,INTENT(IN) :: LINEFEED ! new line at least each LINEFEED elements of RES 

!*------------- LOCAL DATA -----------------------------------------------------
  INTEGER LOC_LF 
    
    IF (PRESENT(LINEFEED)) THEN
      LOC_LF = LINEFEED 
    ELSE
      LOC_LF = SIZE(RES,1)
    END IF

    CALL WRITE_REAL_VECTOR(TOKEN,RESHAPE(RES,(/SIZE(RES)/)),COLUMN,DIGITS,FMT,LOC_LF)

  END SUBROUTINE WRITE_REAL_MATRIX

! ------------------------------------------------------------------------------
 
 SUBROUTINE WRITE_INT__SCALAR(TOKEN,RES,COLUMN,FMT)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN  ! string, identifies item in SETUP
  INTEGER                ,INTENT(IN) :: RES    ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN ! output should start at this column
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT    ! format descriptor for RES

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*20        :: LOC_FMT ! format descriptor 
   CHARACTER*50        :: CHR_RES ! RES written with format FMT
   INTEGER             :: LENRES  ! length of string CHR_RES
       
!  write TOKEN at position COLUMN

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
     ELSE
         CALL GETFMT(RES,LOC_FMT)
     END IF

!  write result into CHR_RES

     WRITE(CHR_RES,LOC_FMT) RES     
 
!  write CHR_RES at column ICOL or into next line if not enough space

     LENRES=LEN_TRIM(CHR_RES)
     IF (ICOL + LENRES + 1 <= MAXCOL_OUT_SETUPLIB) THEN
       WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
     ELSE
       ICOL = COLCAT_OUT_SETUPLIB
       WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
     END IF
       
     ICOL = ICOL + LENRES + 1

  END SUBROUTINE WRITE_INT__SCALAR

! ------------------------------------------------------------------------------

 SUBROUTINE WRITE_INT__VECTOR(TOKEN,RES,COLUMN,FMT,LINEFEED)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN    ! string, identifies item in SETUP
  INTEGER                ,INTENT(IN) :: RES(:)   ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN   ! output should start at this column
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT      ! format descriptor for RES
  INTEGER      ,OPTIONAL ,INTENT(IN) :: LINEFEED ! new line at least each LINEFEED elements of RES 

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*200       :: LOC_FMT    ! format descriptor 
   CHARACTER*50        :: CHR_RES    ! RES written with format FMT
   INTEGER             :: LENRES     ! length of string CHR_RES
   INTEGER             :: N          ! size of RES
   INTEGER             :: IR         ! runs over range of RES
   INTEGER             :: ICOL_START ! column at which first value of RES is written 
   INTEGER             :: LOC_LF     ! new line at least each LINEFEED elements of RES     

     N = SIZE(RES)

!  write TOKEN at position COLUMN (column =COLCAT_OUT_SETUPLIB+COLUMN)

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

     IF (PRESENT(LINEFEED)) THEN
         LOC_LF = LINEFEED 
     ELSE
         LOC_LF = N
     END IF

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
     ELSE
         CALL GETFMT(RES,LOC_FMT,MODE='MAX')
     END IF
     
     ICOL_START=ICOL
     
     DO IR=1,N

!  write result into CHR_RES

       WRITE(CHR_RES,LOC_FMT) RES(IR) ! same format for all elements
 
!  write CHR_RES at column ICOL or into next line if not enough space

       LENRES=LEN_TRIM(CHR_RES)

       IF ((ICOL + LENRES + 1 > MAXCOL_OUT_SETUPLIB) .OR. (MOD(IR,LOC_LF)==1 .AND. IR/=1)) THEN
         ICOL = ICOL_START
         WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
       ELSE
         WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
       END IF

       ICOL = ICOL + LENRES + 1
     
     END DO 


 END SUBROUTINE WRITE_INT__VECTOR

! ------------------------------------------------------------------------------

 SUBROUTINE WRITE_LOG__SCALAR(TOKEN,RES,COLUMN,FMT)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN  ! string, identifies item in SETUP
  LOGICAL                ,INTENT(IN) :: RES    ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN ! output should start at this column
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT    ! format descriptor for RES

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*20        :: LOC_FMT ! format descriptor 
   CHARACTER*50        :: CHR_RES ! RES written with format FMT
   INTEGER             :: LENRES  ! length of string CHR_RES
       
!  write TOKEN at position COLUMN

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
     ELSE
         CALL GETFMT(RES,LOC_FMT)
     END IF

!  write result into CHR_RES

     WRITE(CHR_RES,LOC_FMT) RES     
 
!  write CHR_RES at column ICOL or into next line if not enough space

     LENRES=LEN_TRIM(CHR_RES)
     IF (ICOL + LENRES + 1 <= MAXCOL_OUT_SETUPLIB) THEN
       WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
     ELSE
       ICOL = COLCAT_OUT_SETUPLIB
       WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
     END IF
       
     ICOL = ICOL + LENRES + 1

 END SUBROUTINE WRITE_LOG__SCALAR

! ------------------------------------------------------------------------------

 SUBROUTINE WRITE_LOG__VECTOR(TOKEN,RES,COLUMN,FMT,LINEFEED)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN    ! string, identifies item in SETUP
  LOGICAL                ,INTENT(IN) :: RES(:)   ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN   ! output should start at this column
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT      ! format descriptor for RES
  INTEGER      ,OPTIONAL ,INTENT(IN) :: LINEFEED ! new line at least each LINEFEED elements of RES 

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*20        :: LOC_FMT    ! format descriptor 
   CHARACTER*50        :: CHR_RES    ! RES written with format FMT
   INTEGER             :: LENRES     ! length of string CHR_RES
   INTEGER             :: N          ! size of RES
   INTEGER             :: IR         ! runs over range of RES
   INTEGER             :: ICOL_START ! column at which first value of RES is written 
   INTEGER             :: LOC_LF     ! new line at least each LINEFEED elements of RES     

     N = SIZE(RES)

!  write TOKEN at position COLUMN (column =COLCAT_OUT_SETUPLIB+COLUMN)

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

     IF (PRESENT(LINEFEED)) THEN
         LOC_LF = LINEFEED 
     ELSE
         LOC_LF = N
     END IF

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
     ELSE
         CALL GETFMT(RES,LOC_FMT,MODE='MAX') 
     END IF
     
     ICOL_START=ICOL
     
     DO IR=1,N

!  write result into CHR_RES

       WRITE(CHR_RES,LOC_FMT) RES(IR) ! same format for all elements
 
!  write CHR_RES at column ICOL or into next line if not enough space

       LENRES=LEN_TRIM(CHR_RES)

       IF ((ICOL + LENRES + 1 > MAXCOL_OUT_SETUPLIB) .OR. (MOD(IR,LOC_LF)==1 .AND. IR/=1)) THEN
         ICOL = ICOL_START
         WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
       ELSE
         WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
       END IF

       ICOL = ICOL + LENRES + 1
     
     END DO 

END SUBROUTINE WRITE_LOG__VECTOR

! ------------------------------------------------------------------------------

  SUBROUTINE WRITE_CHAR_SCALAR(TOKEN,RES,COLUMN,FMT)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN  ! string, identifies item in SETUP
  CHARACTER*(*)          ,INTENT(IN) :: RES    ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN ! output should start at this column
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT    ! format descriptor for RES

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*100       :: LOC_FMT ! format descriptor 
   CHARACTER*50        :: CHR_RES ! RES written with format FMT
   INTEGER             :: LENRES  ! length of string CHR_RES
       
!  write TOKEN at position COLUMN

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

!  eventually determines the format of RES

     IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
     ELSE
         CALL GETFMT(RES,LOC_FMT,APOSTROPHES=.TRUE.)
     END IF

!  write result into CHR_RES and determine length
     
     WRITE(CHR_RES,LOC_FMT) RES     
     LENRES=LEN_TRIM(CHR_RES)
      
!  write RES at column ICOL or into next line if not enough space


     IF (ICOL + LENRES + 1 <= MAXCOL_OUT_SETUPLIB) THEN
       WRITE(UNIT_OUT,'('//TRIM(LOC_FMT)//',a,$)')RES,BLANK
     ELSE
       ICOL = COLCAT_OUT_SETUPLIB
       WRITE(UNIT_OUT,'(/,a,'//LOC_FMT//',a,$)')REPEAT(BLANK,ICOL),RES,BLANK
     END IF
       
     ICOL = ICOL + LENRES + 1

 END SUBROUTINE WRITE_CHAR_SCALAR
 
! ------------------------------------------------------------------------------

 SUBROUTINE WRITE_CHAR_VECTOR(TOKEN,RES,COLUMN,FMT,LINEFEED)

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)          ,INTENT(IN) :: TOKEN  ! string, identifies item in SETUP
  CHARACTER*(*)          ,INTENT(IN) :: RES(:) ! value of TOKEN
  INTEGER      ,OPTIONAL ,INTENT(IN) :: COLUMN ! output should start at this column
  CHARACTER*(*),OPTIONAL ,INTENT(IN) :: FMT    ! format descriptor for RES
  INTEGER      ,OPTIONAL ,INTENT(IN) :: LINEFEED ! new line at least each LINEFEED elements of RES 

!*------------- LOCAL DATA -----------------------------------------------------

   CHARACTER*1000      :: LOC_FMT    ! format descriptor 
   CHARACTER*1000      :: CHR_RES    ! RES written with format FMT
   INTEGER             :: LENRES     ! length of string CHR_RES
   INTEGER             :: N          ! size of RES
   INTEGER             :: IR         ! runs over range of RES
   INTEGER             :: ICOL_START ! column at which first value of RES is written 
   INTEGER             :: LOC_LF     ! new line at least each LINEFEED elements of RES     

     N = SIZE(RES)

!  write TOKEN at position COLUMN (column =COLCAT_OUT_SETUPLIB+COLUMN)

     CALL WRITE_TOKEN(TOKEN,COLUMN)       

     IF (PRESENT(LINEFEED)) THEN
         LOC_LF = LINEFEED 
     ELSE
         LOC_LF = N
     END IF
     
     ICOL_START=ICOL
     
     DO IR=1,N

!  eventually determines the format of RES

       IF (PRESENT(FMT)) THEN
         LOC_FMT = FMT 
       ELSE
         CALL GETFMT(RES(IR),LOC_FMT,APOSTROPHES=.TRUE.)
       END IF

!  write result into CHR_RES

       WRITE(CHR_RES,LOC_FMT) RES(IR) ! maybe different format for different elements    
 
!  write CHR_RES at column ICOL or into next line if not enough space

       LENRES=LEN_TRIM(CHR_RES)

       IF ((ICOL + LENRES + 1 > MAXCOL_OUT_SETUPLIB) .OR. (MOD(IR,LOC_LF)==1 .AND. IR/=1)) THEN
         ICOL = ICOL_START
         WRITE(UNIT_OUT,'(/,3a,$)') REPEAT(BLANK,ICOL),TRIM(CHR_RES),BLANK
       ELSE
         WRITE(UNIT_OUT,'(  2a,$)')                    TRIM(CHR_RES),BLANK
       END IF

       ICOL = ICOL + LENRES + 1
     
     END DO 

  END SUBROUTINE WRITE_CHAR_VECTOR

 


!*******************************************************************************
! 
!*******************************************************************************






 SUBROUTINE IO_REAL_SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                           COLUMN,DIGITS,FMT,LPRMOD,LPR)
  
! ------------- ARGUMENTS--------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  REAL(KIND=LF)         ,INTENT(INOUT)::RES       ! array from/into which data is read
  REAL(KIND=LF),OPTIONAL,INTENT(IN   )::DEFAULT   ! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  INTEGER      ,OPTIONAL,INTENT(IN)   ::DIGITS    ! maximum number of decimal digits
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL       :: LOC_LPR       !
  CHARACTER*3   :: LOC_LPRMOD    !
  REAL(KIND=LF) :: RES0          !


  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN

    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT
    
    CALL READ_REAL_SCALAR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF             
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. RES/=RES0 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. RES/=DEFAULT 
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123    

    IF (LOC_LPR) CALL WRITE_REAL_SCALAR(TOKEN,RES,COLUMN,DIGITS,FMT)

  END IF

 END SUBROUTINE IO_REAL_SCALAR

! ------------------------------------------------------------------------------

 SUBROUTINE IO_REAL_VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                           COLUMN,DIGITS,FMT,LINEFEED,LPRMOD,LPR)  

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  REAL(KIND=LF)         ,INTENT(INOUT)::RES(:)    ! array from/into which data is read
  REAL(KIND=LF),OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  INTEGER      ,OPTIONAL,INTENT(IN)   ::DIGITS    ! maximum number of decimal digits
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  INTEGER      ,OPTIONAL,INTENT(IN)   ::LINEFEED  ! new line at least each LINEFEED elements of RES 
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

  LOGICAL       :: LOC_LPR        !
  CHARACTER*3   :: LOC_LPRMOD     !
  REAL(KIND=LF) :: RES0(SIZE(RES))!

  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN
    
    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT

    CALL READ_REAL_VECTOR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF              
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. ANY(RES/=RES0) 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. ANY(RES/=DEFAULT)
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_REAL_VECTOR(TOKEN,RES,COLUMN,DIGITS,FMT,LINEFEED)

  END IF

 END SUBROUTINE IO_REAL_VECTOR

! ------------------------------------------------------------------------------

 SUBROUTINE IO_REAL_MATRIX(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                               COLUMN,DIGITS,FMT,LINEFEED,LPRMOD,LPR)
  
! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  REAL(KIND=LF)         ,INTENT(INOUT)::RES(:,:)  ! array from/into which data is read
  REAL(KIND=LF),OPTIONAL,INTENT(IN   )::DEFAULT(:,:)! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  INTEGER      ,OPTIONAL,INTENT(IN)   ::DIGITS    ! maximum number of decimal digits
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  INTEGER      ,OPTIONAL,INTENT(IN)   ::LINEFEED  ! new line at least each LINEFEED elements of RES 
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL       :: LOC_LPR        !
  CHARACTER*3   :: LOC_LPRMOD     !
  REAL(KIND=LF) :: RES0(SIZE(RES,1),SIZE(RES,2))!

  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN
    
    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT

    CALL READ_REAL_MATRIX (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF                          
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. ANY(RES/=RES0) 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. ANY(RES/=DEFAULT)
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_REAL_MATRIX(TOKEN,RES,COLUMN,DIGITS,FMT,LINEFEED)

  END IF

 END SUBROUTINE IO_REAL_MATRIX

! ------------------------------------------------------------------------------

 SUBROUTINE IO_INT__SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                           COLUMN,FMT,LPRMOD,LPR)
  
! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  INTEGER               ,INTENT(INOUT)::RES       ! array from/into which data is read
  INTEGER      ,OPTIONAL,INTENT(IN   )::DEFAULT   ! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL       :: LOC_LPR       !
  CHARACTER*3   :: LOC_LPRMOD    !
  INTEGER       :: RES0          !


  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN

    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT
    
    CALL READ_INT__SCALAR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF             
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. RES/=RES0 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. RES/=DEFAULT 
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_INT__SCALAR(TOKEN,RES,COLUMN,FMT)

  END IF

 END SUBROUTINE IO_INT__SCALAR
  
! ------------------------------------------------------------------------------

SUBROUTINE IO_INT__VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                              COLUMN,FMT,LINEFEED,LPRMOD,LPR)
  
! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  INTEGER               ,INTENT(INOUT)::RES(:)    ! array from/into which data is read
  INTEGER      ,OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  INTEGER      ,OPTIONAL,INTENT(IN)   ::LINEFEED  ! new line at least each LINEFEED elements of RES 
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL     :: LOC_LPR        !
  CHARACTER*3 :: LOC_LPRMOD     !
  INTEGER     :: RES0(SIZE(RES))!

  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN
    
    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT
    
    CALL READ_INT__VECTOR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF                         
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. ANY(RES/=RES0) 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. ANY(RES/=DEFAULT)
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_INT__VECTOR(TOKEN,RES,COLUMN,FMT,LINEFEED)

  END IF

 END SUBROUTINE IO_INT__VECTOR 
 
! ------------------------------------------------------------------------------

 
 SUBROUTINE IO_CHAR_SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                               COLUMN,FMT,LPRMOD,LPR)

   
  

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  CHARACTER*(*)         ,INTENT(INOUT)::RES       ! array into which data is read
  CHARACTER*(*),OPTIONAL,INTENT(IN   )::DEFAULT   ! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------



  LOGICAL              :: LOC_LPR        !
  CHARACTER*3          :: LOC_LPRMOD     !
  CHARACTER*(LEN(RES)) :: RES0           !


  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN
    
    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT
    
    CALL READ_CHAR_SCALAR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF                           
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. RES/=RES0 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. RES/=DEFAULT 
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_CHAR_SCALAR(TOKEN,RES,COLUMN,FMT)

  END IF

 END SUBROUTINE IO_CHAR_SCALAR

! ------------------------------------------------------------------------------

 SUBROUTINE IO_CHAR_VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,           &
                               COLUMN,FMT,LINEFEED,LPRMOD,LPR)
 
  

! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  CHARACTER*(*)         ,INTENT(INOUT)::RES    (:)! array into which data is read
  CHARACTER*(*),OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  INTEGER      ,OPTIONAL,INTENT(IN)   ::LINEFEED  ! new line at least each LINEFEED elements of RES 
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL     :: LOC_MANDATORY ! needed due to internal compiler error cannot pass DEFAULT to READ_CHAR_VECTOR
  LOGICAL     :: LOC_LPR       !
  CHARACTER*3 :: LOC_LPRMOD    !

  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN


    IF (PRESENT(DEFAULT)) RES=DEFAULT

    LOC_MANDATORY = DEFAULT_MANDATORY_SETUPLIB .AND. .NOT.(PRESENT(DEFAULT).OR.PRESENT(FOUND))
    IF (PRESENT(MANDATORY)) LOC_MANDATORY = MANDATORY    
    IF (PRESENT(DEFAULT)) RES=DEFAULT
    
    CALL READ_CHAR_VECTOR (TOKEN,RES,MANDATORY=LOC_MANDATORY,REWIND=REWIND,FOUND=FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF                          
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y' .AND. PRESENT(DEFAULT)) LPR = LPR .OR. ANY(RES/=DEFAULT) 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. ANY(RES/=DEFAULT)
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_CHAR_VECTOR(TOKEN,RES,COLUMN,FMT,LINEFEED)

  END IF

 END SUBROUTINE IO_CHAR_VECTOR

! ------------------------------------------------------------------------------

 SUBROUTINE IO_LOG__SCALAR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,  &
                              COLUMN,FMT,LPRMOD,LPR)
  
! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  LOGICAL               ,INTENT(INOUT)::RES       ! array from/into which data is read
  LOGICAL      ,OPTIONAL,INTENT(IN   )::DEFAULT   ! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL       :: LOC_LPR       !
  CHARACTER*3   :: LOC_LPRMOD    !
  LOGICAL       :: RES0          !


  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN

    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT
    
    CALL READ_LOG__SCALAR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF            
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. (RES.NEQV.RES0) 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. (RES.NEQV.DEFAULT)  
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_LOG__SCALAR(TOKEN,RES,COLUMN,FMT)

  END IF

 END SUBROUTINE IO_LOG__SCALAR

! ------------------------------------------------------------------------------

SUBROUTINE IO_LOG__VECTOR(TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND,            &
                              COLUMN,FMT,LINEFEED,LPRMOD,LPR)
  
! ------------- ARGUMENTS-------------------------------------------------------

  CHARACTER*(*)         ,INTENT(IN)   ::TOKEN     ! string, identifies item in SETUP
  LOGICAL               ,INTENT(INOUT)::RES(:)    ! array from/into which data is read
  LOGICAL      ,OPTIONAL,INTENT(IN   )::DEFAULT(:)! default value for RES
  LOGICAL      ,OPTIONAL,INTENT(IN   )::MANDATORY ! TRUE -> aborts if unmatched TOKEN
  LOGICAL      ,OPTIONAL,INTENT(IN   )::REWIND    ! TRUE -> reads from start of cat.
  LOGICAL      ,OPTIONAL,INTENT(  OUT)::FOUND     ! true if TOKEN found
  INTEGER      ,OPTIONAL,INTENT(IN)   ::COLUMN    ! output should start at this column
  CHARACTER*(*),OPTIONAL,INTENT(IN)   ::FMT       ! format descriptor for RES
  INTEGER      ,OPTIONAL,INTENT(IN)   ::LINEFEED  ! new line at least each LINEFEED elements of RES 
  CHARACTER*3  ,OPTIONAL,INTENT(IN)   ::LPRMOD    ! controls whether token is written
  LOGICAL      ,OPTIONAL,INTENT(INOUT)::LPR       ! TRUE -> token will be be written

!*------------- LOCAL DATA -----------------------------------------------------

  LOGICAL     :: LOC_LPR       !
  CHARACTER*3 :: LOC_LPRMOD    !
  LOGICAL     :: RES0(SIZE(RES))!

  IF     (OPTIO_SETUPLIB==0) THEN
    
    CONTINUE

  ELSEIF (OPTIO_SETUPLIB==1) THEN
    
    RES0=RES
    IF (PRESENT(DEFAULT)) RES0=DEFAULT

    CALL READ_LOG__VECTOR (TOKEN,RES,MANDATORY,DEFAULT,REWIND,FOUND)
    
    IF (PRESENT(LPR)) THEN    
      IF (PRESENT(LPRMOD)) THEN
         LOC_LPRMOD = LPRMOD 
      ELSE
         LOC_LPRMOD = 'YYY'
      END IF                       
      IF (LOC_LPRMOD(1:1)=='Y') LPR = LPR .OR. (TOKEN_FOUND .AND. VERBOS_SETUPLIB >= 30) 
      IF (LOC_LPRMOD(2:2)=='Y') LPR = LPR .OR. ANY(RES.NEQV.RES0) 
      IF (LOC_LPRMOD(3:3)=='Y') LPR = LPR .OR. VERBOS_SETUPLIB >= 50
    END IF

  ELSEIF (OPTIO_SETUPLIB==2) THEN

    IF (PRESENT(LPR)) THEN
      LOC_LPR  = LPR
    ELSEIF (PRESENT(DEFAULT)) THEN
      LOC_LPR  = VERBOS_SETUPLIB >= 50 .OR. ANY(RES.NEQV.RES0) 
    ELSE
      LOC_LPR = .TRUE.       
    END IF
    LOC_LPR  = LOC_LPR .OR. VERBOS_SETUPLIB==123

    IF (LOC_LPR) CALL WRITE_LOG__VECTOR(TOKEN,RES,COLUMN,FMT,LINEFEED)

  END IF

 END SUBROUTINE IO_LOG__VECTOR

! ------------------------------------------------------------------------------

 SUBROUTINE IO_CATEGORY(CATEG,MANDATORY,FOUND,LPR)

! --------------------- ARGUMENTS-----------------------------------------------

  CHARACTER*(*),   INTENT( IN) :: CATEG      ! name of category
  LOGICAL,OPTIONAL,INTENT( IN) :: MANDATORY  ! aborts if unmatched category
  LOGICAL,OPTIONAL,INTENT(OUT) :: FOUND      ! true if category found
  LOGICAL,OPTIONAL,INTENT( IN) :: LPR        ! TRUE -> category will be be written
    
  LOGICAL     :: LOC_LPR       !
  

  IF     (OPTIO_SETUPLIB==0) THEN
    CONTINUE
  ELSEIF (OPTIO_SETUPLIB==1) THEN
    CALL CATEGORY (CATEG,MANDATORY,FOUND)
  ELSEIF (OPTIO_SETUPLIB==2) THEN
    IF (PRESENT(LPR)) THEN
      LOC_LPR =  LPR 
    ELSE
      LOC_LPR = .TRUE.
    END IF
    IF (LOC_LPR) CALL WRITE_CATEGORY(CATEG)
  END IF

 END SUBROUTINE IO_CATEGORY

 END MODULE SETUPLIB
