!!$   $Id: matmul_blas.f,v 1.4 2003/11/03 19:04:52 maciej Exp $
Module matmul_mod
!!$    This module defines the wrappers for matrix multiplications
!!$    calling BLAS *gemm functions.
!!$    Uncomment the entries for *_re and *_ce functions
!!$    if prec.NE.lprec in definitions.f
!!$    use parameters

      Interface mmatmul
         Module Procedure mmatmul_r, mmatmul_c
      End Interface


Contains


      Function mmatmul_r (a, b)
         Use definitions

         Implicit None
         Real (Kind=prec), Intent (In) :: a (:, :), b (:, :)
         Real (Kind=prec) :: mmatmul_r (size(a, 1), size(b, 2))

         Real (Kind=prec) :: alpha, beta

         alpha = 1.0d0
         beta = 0.0d0

         Call dgemm ('N', 'N', size(a, 1), size(b, 2), size(a, 2), alpha, a, size(a, 1), b, size(b, 1), beta, &
        & mmatmul_r, size(a, 1))


      End Function mmatmul_r




      Function mmatmul_c (a, b)
         Use definitions

         Implicit None
         Complex (Kind=prec), Intent (In) :: a (:, :), b (:, :)
         Complex (Kind=prec) :: mmatmul_c (size(a, 1), size(b, 2))

         Complex (Kind=prec) :: alpha, beta

         alpha = 1.0d0
         beta = 0.0d0

         Call zgemm ('N', 'N', size(a, 1), size(b, 2), size(a, 2), alpha, a, size(a, 1), b, size(b, 1), beta, &
        & mmatmul_c, size(a, 1))

      End Function mmatmul_c

      Subroutine mmatmul_hz (c, a, b)
         Use definitions
         Implicit None
         Complex (Kind=prec) :: a (:, :), b (:, :)
         Complex (Kind=prec) :: c (size(a, 2), size(b, 2))
         Complex (Kind=prec) :: alpha, beta
         alpha = 1.0d0
         beta = 0.0d0
         Call zgemm ('C', 'N', size(a, 2), size(b, 2), size(a, 1), alpha, a, size(a, 2), b, size(b, 1), beta, &
        & c, size(a, 2))
      End Subroutine mmatmul_hz

      Subroutine mmatmul_z (c, a, b)
         Use definitions
         Implicit None
         Complex (Kind=prec) :: a (:, :), b (:, :)
         Complex (Kind=prec) :: c (size(a, 1), size(b, 2))
         Complex (Kind=prec) :: alpha, beta
         alpha = 1.0d0
         beta = 0.0d0
         Call zgemm ('N', 'N', size(a, 1), size(b, 2), size(a, 2), alpha, a, size(a, 1), b, size(b, 1), beta, &
        & c, size(a, 1))
      End Subroutine mmatmul_z




End Module matmul_mod
