!!$ $Id: math_def.h 101 2006-10-19 10:31:50Z antst $
!!$-*-F90-*- !This line forces emacs to handle it as fortran source file
!!$ This header contains the definitions of data structures used elsewhere
!!$ in the program and declarations od some constants.


#if !defined(_antst_math_defs_)

!!$ Size of default float datatypes
#define DEF_SING_PREC 4
#define DEF_DBL_PREC 8


!!$ "0" in default precision complex variable. (Im tired to define it everywhere)
#define DEF_cmplx_zero cmplx(0.0d0,0.0d0,kind=DEF_DBL_PREC)

!!$ "1" in default precision complex variable. (Im tired to define it everywhere)
#define DEF_cmplx_one cmplx(1.0d0,0.0d0,kind=DEF_DBL_PREC)

!!$ "i" in default precision complex variable. (Im tired to define it everywhere)
#define DEF_cmplx_Ione cmplx(0.0d0,1.0d0,kind=DEF_DBL_PREC)


#if defined(HAVE_QUADPREC)
!!$ Size of extended float datatype
#define DEF_QUAD_PREC 16

!!$ "0" in long precision real variable. (Im tired to define it everywhere)
#define DEF_long_zero 0.0q0

!!$ "1" in long precision real variable. (Im tired to define it everywhere)
#define DEF_long_one 1.0q0

!!$ "2" in long precision real variable. (Im tired to define it everywhere)
#define DEF_long_two 2.0q0

!!$ "0" in default precision complex variable. (Im tired to define it everywhere)
#define DEF_cmplx_zeroL cmplx(0.0q0,0.0q0,kind=DEF_QUAD_PREC)

!!$ "1" in default precision complex variable. (Im tired to define it everywhere)
#define DEF_cmplx_oneL cmplx(1.0q0,0.0q0,kind=DEF_QUAD_PREC)

!!$ "i" in default precision complex variable. (Im tired to define it everywhere)
#define DEF_cmplx_IoneL cmplx(0.0q0,1.0q0,kind=DEF_QUAD_PREC)

#endif

!!$Defention of PI and some derevated constants in default precision floats variables
!!$  pi strange...fortran doesnt like normal PI
#define DEF_M_PI real(3.1415926535897932384626433832795029d0,kind=DEF_DBL_PREC)
!!$ pi/2 
#define DEF_M_PI_2 real(1.5707963267948966192313216916397514d0,kind=DEF_DBL_PREC)
!!$ pi/4 
#define DEF_M_PI_4 real(0.7853981633974483096156608458198757d0,kind=DEF_DBL_PREC)
!!$ 1/pi
#define DEF_M_1_PI real(0.3183098861837906715377675267450287d0,kind=DEF_DBL_PREC)
!!$ 2/pi
#define DEF_M_2_PI real(0.6366197723675813430755350534900574d0,kind=DEF_DBL_PREC) 
#define _antst_math_defs_ 1
#endif
