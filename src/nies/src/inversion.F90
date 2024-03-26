!!$   $Id: inversion.f,v 1.2 2003/06/17 16:01:37 maciej Exp $
Module inversion
!!$    This module defines the wrapper for complex*16 matrix inversion

Contains
      Subroutine lcinv (A)
         Use definitions
!!$     Inversion of complex matrix (wrapper of Lapack routines)
         Implicit None
         Complex (Kind=prec) :: a (:, :)
         Complex (Kind=prec) :: work (size(a, 1))
         Integer :: ipiv (size(a, 1)), info, n

         n = size (a, 1)
         Call zgetrf (n, n, a, n, ipiv, info)! computes LU factorization
         If (info /= 0) Then
            Write (*,*) 'Problems encountered while in zgetrf routine, info =', info
            Stop
         End If
         Call zgetri (n, a, n, ipiv, work, n, info)! inverses the a matrix
         
         If (info /= 0) Then
            Write (*,*) 'Problems encountered while in zgetrfi routine, info =', info
            Stop
         End If
      End Subroutine lcinv
End Module inversion
