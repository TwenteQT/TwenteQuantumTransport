!!$ Module which contains the subroutines related to the non-spherical corrections.
#include "math_def.h"

Module omta_df
 Use sparselib
 Implicit None

 Contains

   Subroutine Downfold_matrix(A,to_df)
   !!$ Downfold the states given in to_df in the matrix A in-place
   !!$ Note that A still has the larger dimensions...
   Implicit None
   Complex (Kind=DEF_DBL_PREC) :: A(:,:) !!$ to-downfold matrix
   Integer :: to_df(:)        !!$ the array of to-downfold entries of vector x
   !!$ local
   Complex (Kind=DEF_DBL_PREC) :: B(size(A,1),size(A,2)) !!$ intermediate matrix
   Integer :: to_df2(size(to_df))  !!$ don't mess up the original array

   to_df2 = Reorder(A,to_df)
   B = 0d0
   B(1:size(A,1)-size(to_df), 1:size(A,2)-size(to_df)) = ActuallyDownfold(A,to_df2)
   A = B

   End Subroutine

   Function Reorder(A,to_df) result (to_df2)
     Implicit None 
     Complex (Kind=DEF_DBL_PREC) :: A(:,:) !!$ to-downfold matrix
     Integer :: to_df(:)
     !!$ local
     Real (Kind=DEF_DBL_PREC) :: B(size(A,1),size(A,2))
     Integer :: to_df2(size(to_df))
     Integer :: i,j
     Integer :: n,ndf,idf,jdf
     Integer :: i_n, j_n
     n = size(A,1)
     ndf = size(to_df)
     idf = 0
     i_n = 0
 
     !reorder the rows
     do i=1,n
      if (ANY(to_df==i)) then
        idf = idf+1
        B(n-ndf+idf,:) = A(i,:)
      else
        i_n = i_n+1
        B(i_n,:) = A(i,:)
      endif
     enddo
 
     jdf = 0
     j_n = 0
     !reorder the columns
     do j=1,n
      if (ANY(to_df==j)) then
        jdf = jdf+1
        A(:,n-ndf+jdf) = B(:,j)
      else
        j_n = j_n + 1
        A(:,j_n) = B(:,j)
      endif
     enddo
 
     to_df2 = [(n-ndf+j,j=1,ndf)]
 
   End function

   Function ActuallyDownfold(A,to_df2) result (B)
    Implicit None
    Complex (Kind=DEF_DBL_PREC) :: A(:,:) !!$ to-downfold matrix
    Integer :: to_df2(:)
    Complex (Kind=DEF_DBL_PREC) :: B(size(A,1)-size(to_df2),size(A,2)-size(to_df2))
    !!$local
    Integer :: n, ndf
    n = size(A,1)
    ndf = size(to_df2)

    B = A(1:n-ndf, 1:n-ndf) - MATMUL(A(1:n-ndf, n-ndf+1:n),MATMUL( &
      & invert_cmatrix(ndf,A(n-ndf+1:n, n-ndf+1:n),ndf),A(n-ndf+1:n,1:n-ndf)))

   End Function

      Function invert_cmatrix (n, a, lda) result(aout)
!!$    a wrapper for the lapack routines for inverting a complex
!!$    square matrix
         Implicit None
         !External zgetrf, zgetri

         Integer, Intent (In) :: n, lda
         Complex (Kind=DEF_DBL_PREC), Intent (In) :: a (:, :)
         Complex (Kind=DEF_DBL_PREC) :: aout (size(a,1),size(a,2))

         Complex (Kind=DEF_DBL_PREC), Pointer :: work (:)
         Integer :: info
         Integer, Pointer :: ipiv (:)

         Allocate (ipiv(n), work(n))

!!$ DAMMMIT! Valgrind complains here!
         Call zgetrf (n, n, a, lda, ipiv, info)! computes lu factorization (lapack
!!$


         If (info /= 0) Then
            Write (*,*) 'problems encountered while in zgetrf routine, info =', info
            Write (*,*) 'in omta_DF module' 
            Stop
         End If
         Call zgetri (n, a, lda, ipiv, work, n, info)! inverses the a matrix
         If (info /= 0) Then
            Write (*,*) 'problems encountered while in zgetrfi routine, info =', info
            Stop
         End If
         Deallocate (ipiv, work)
        
         aout = REAL(a)

      End Function invert_cmatrix

   Subroutine Downfold_iter(A,to_df,en)
   !!$ iteratively downfold the states given in to_df in the matrix A in-place
   !!$ Note that A still has the larger dimensions...
   Implicit None
   Complex (Kind=DEF_DBL_PREC) :: A(:,:) !!$ to-downfold matrix
   Real (Kind=DEF_DBL_PREC) :: en
   Integer :: to_df(:)        !!$ the array of to-downfold entries of vector x
   !!$ local
   Integer :: n, i, k, j
   Complex (Kind=DEF_DBL_PREC) :: C(size(A,1),size(A,2)) !!$ intermediate matrix
   Integer :: df(size(to_df))  !!$ don't mess up the original array

   n = size(A,1)

   df = to_df

   do i=1,size(df)
     k = n-i
     call df_once(k+1,A(1:k+1, 1:k+1),C(1:k, 1:k),df(i),en)
     df = df - 1 !!$ reduce indices by one
     A = 0d0
     A(1:k, 1:k) = C(1:k, 1:k)
   enddo

   End Subroutine

   Subroutine df_once(n,A,B,df_i, en)
   !!$ downfolds one variable of the output vector into the matrix.
   !!$ 'u' means upper, 'l' means lower, which is the to-downfold variable
   !!$ that means 'll' is the diagonal value of the variable, which is inverted
   !as 1 / ll

   Implicit None
   Complex (Kind=DEF_DBL_PREC) :: A(n,n) !!$ to-downfold matrix
   Complex (Kind=DEF_DBL_PREC) :: B(n-1,n-1) !!$ downfolded matrix
   Real (Kind=DEF_DBL_PREC) :: en
   Integer :: n, df_i !!$ index of to-downfold entry
   !!$ local
   Real (Kind=DEF_DBL_PREC) :: uu(n-1,n-1)
   Real (Kind=DEF_DBL_PREC) :: ul(n-1,1), lu(1,n-1)
   Real (Kind=DEF_DBL_PREC) :: ll
   Integer :: lo, up, m, i, j

   lo = df_i-1
   up = df_i+1
   m  = n-1

   if (df_i == 1) then
     ll = A(df_i,df_i)
     ul(:,1) = A(df_i,up:n)
     lu(1,:) = A(up:n,df_i)
     uu = A(df_i:n,df_i:n)
   else if (df_i == n) then
     ll = A(df_i,df_i)
     ul(:,1) = A(df_i,1:lo)
     lu(1,:) = A(1:lo, df_i)
     uu = A(1:lo, 1:lo)
   else
    ll = A(df_i,df_i)
    ul(:,1) = [A(df_i, 1:lo), A(df_i, up:n)]
    lu(1,:) = [A(1:lo, df_i), A(up:n, df_i)]

    uu(1:lo  , 1:lo) = A(1:lo, 1:lo)
    uu(lo+1:m, 1:lo) = A(up:n, 1:lo)
    uu(1:lo  , lo+1:m) = A(1:lo, up:n)
    uu(lo+1:m, lo+1:m) = A(up:n, up:n)
   endif

   B = uu - MATMUL(ul,lu)/(ll-en)

!   do i=1,m
!    do j=1,m
!     B(i,j) = uu(i,j) - ul(j) * 1d0/(ll-en) * lu(i)
!    enddo
!   enddo
   End Subroutine



End

