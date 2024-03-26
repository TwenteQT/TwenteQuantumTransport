#include "math_def.h"

Program eig_test
      Use sparselib
      Use potpars
      Implicit None
      Integer, Parameter :: n = 144
      Complex (Kind=DEF_DBL_PREC) :: Lmat (n, n), Rmat (n, n), Lmat0 (n, n), Rmat0 (n, n)
      Complex (Kind=DEF_DBL_PREC) :: alpha (n), beta (n), evec (n, n)
      Integer :: info
      Integer :: i
      Real (Kind=DEF_DBL_PREC) :: res (n)

      Call readzdensemat (Lmat, 'lmatr.dat')
      Call readzdensemat (Rmat, 'rmatr.dat')
      Lmat0 = Lmat
      Rmat0 = Rmat
      Call newzggevx_driver ('B', n, Lmat, Rmat, evec, alpha, beta, info)
!!$       Call r16tomsqz_driver ( n, Lmat, Rmat, evec, alpha, beta, info)

      Write (*,*) info

      Write (*,*) alpha (1) / beta (1)
      Write (*,*) alpha (n) / beta (n)

!!$       Write (*,*) alpha (8) / beta (8)
      Do i = 1, n
             write(*,*) beta(i)
         res (i) = maxval (Abs(matmul(Lmat0(:, :), evec(1:n, i))-(matmul(Rmat0(:, :), evec(1:n, &
        & i))*(alpha(i)/beta(i)))))
      End Do
      write(*,*)maxval(res)
      Write (*,*) 'done!'
Contains
         Subroutine r16tomsqz_driver (n, Lm, Rm, wr, alpha, beta, info)
!!$ driver for TOMSQZ routine in quadro.prec.
            Implicit None
            Integer :: n, info
            Complex (Kind=DEF_DBL_PREC) :: wr (n, n), alpha (n), beta (n)
            Complex (Kind=DEF_DBL_PREC) :: Lm (n, n), Rm (n, n)
!!$ Local
            Real (Kind=DEF_QUAD_PREC), Pointer :: LmR (:, :), LmI (:, :), RmR (:, :), RmI (:, :)
            Real (Kind=DEF_QUAD_PREC), Pointer :: alphaR (:), alphaI (:), betaQ (:), evR (:, :), evI (:, :)
write(*,*) 'quad'
            Allocate (LmR(n, n), LmI(n, n), RmR(n, n), RmI(n, n), alphaR(n), alphaI(n), betaQ(n), evR(n, n), &
           & evI(n, n))
            LmR = DEF_long_zero
            LmI = DEF_long_zero
            RmR = DEF_long_zero
            RmI = DEF_long_zero
            LmR = real (Lm)
            LmI = aimag (Lm)
            RmR = real (Rm)
            RmI = aimag (Rm)

            Call r16tomsqz (n, n, LmR, LmI, RmR, RmI, alphaR, alphaI, betaQ, evR, evI, info)

            wr = cmplx (evR, evI, kind=DEF_DBL_PREC)
            alpha = cmplx (alphaR, alphaI, kind=DEF_DBL_PREC)
            beta = real (betaQ, kind=DEF_DBL_PREC)
            Deallocate (LmR, LmI, RmR, RmI, alphaR, alphaI, betaQ, evR, evI)
         End Subroutine r16tomsqz_driver
      Subroutine newzggevx_driver (bal, n, Lm, Rm, wr, alpha, beta, info)
!!$ ZGGEVX driver
         Implicit None
         Integer :: n, info
         Character :: bal
         Complex (Kind=DEF_DBL_PREC) :: beta (n), wr (n, n), alpha (n)
         Complex (Kind=DEF_DBL_PREC) :: Lm (n, n), Rm (n, n)
!!$ Local
         Integer :: lwork, ilo, ihi
         Complex (Kind=DEF_DBL_PREC) :: zerop
         Real (Kind=DEF_DBL_PREC) :: abnrm, bbnrm
         Complex (Kind=DEF_DBL_PREC), Pointer :: work (:)
         Real (Kind=DEF_DBL_PREC), Pointer :: rwork (:)
         Logical, Pointer :: bwork (:)
         Integer, Pointer :: lscale (:), rscale (:), iwork (:)
         lwork = 5 * (2*n*n+2*n)

         Allocate (work(lwork), lscale(n), rscale(n), iwork(n+2), rwork(6*n), bwork(n))
#if defined(NTS_LPK)	
!$OMP critical
#endif
         Call ZGGEVX (bal, 'N', 'V', 'N', n, Lm, n, Rm, n, alpha, beta, zerop, n, wr, n, ilo, ihi, lscale, &
        & rscale, abnrm, bbnrm, zerop, zerop, work, lwork, rwork, iwork, bwork, info)

!!$
!!$            SUBROUTINE ZGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA,
!!$                                    B, LDB,  ALPHA,  BETA,  VL,  LDVL,  VR,
!!$                                    LDVR,  ILO, IHI, LSCALE, RSCALE, ABNRM,
!!$                                    BBNRM,  RCONDE,  RCONDV,  WORK,  LWORK,
!!$                                    RWORK, IWORK, BWORK, INFO )


#if defined(NTS_LPK)	
!$OMP end critical
#endif
         Deallocate (work, lscale, rscale, iwork, rwork, bwork)
      End Subroutine newzggevx_driver

      Subroutine readzdensemat (mat, fname)
         Implicit None
         Complex (Kind=DEF_DBL_PREC) :: mat (:, :)
         Character (Len=*) :: fname
         Integer :: s1, s2, nn, i
         Real (Kind=DEF_DBL_PREC) :: ai, ar
         Integer :: dumping_unit = 22823

         Open (Unit=dumping_unit, File=fname, Action='read')

         nn = size (mat, dim=1) * size (mat, dim=2)

         Do i = 1, nn
            Read (dumping_unit,*) s1, s2, ar, ai
            mat (s1, s2) = dcmplx (ar, ai)
         End Do
         Close (dumping_unit)
      End Subroutine readzdensemat
End Program eig_test
