#include "math_def.h"
 Module omta_nonsph
 Use omta_defs
 Use sparselib
 Use geometry_module
 Use structure
 Implicit None
 
 Contains

   Subroutine mk_nonsphblock(geo,amat,sc,kpar)
    !This subroutine generates a giant matrix with the dimensions of the scattering problem from subblocks
    !given by a matrix for each atom and adds it to the input matrix amat
    Implicit None
    Type (t_geometry_EMTO) :: geo
    Type (zcsrmat) :: amat
    Type (t_strconst) :: sc
    Real (Kind=DEF_DBL_PREC) :: kpar(2)
    !local
    Type (zcsrmat) :: V(3), zcsr_tmp
    Type (t_mathop) :: sk
    Integer :: i

    do i=1,3
     call mk_V(geo, i, V(i)) 
    enddo

!    call make_sk_diag(sk, sc, kpar, 0, 1d0)

    zcsr_tmp = spmatadd(amat, V(1)) 
    call spcopy(zcsr_tmp, amat)
!    amat = spmatadd(zcsr_tmp, spmatmul(V(2),sk%c), sign=-1)
!    zcsr_tmp = spmatadd(amat, spmatmul(sk%c,V(2)), sign=-1)
!    amat = spmatadd(zcsr_tmp, spmatmul(sk%c,spmatmul(V(3),sk%c)))
!
!    call free(zcsr_tmp)
!    call free(sk%c)
!    do i=1,3
!     call free(V(i)) 
!    enddo

   End Subroutine mk_nonsphblock

   Subroutine mk_V(geo, ind, zcsr_out)
    ! generates the atom-diagonal part of the nonspherical corrections
    Implicit None
    Type (t_geometry_EMTO) :: geo
    Integer :: ind
    Type (zcsrmat) :: zcsr_out
    !local
    Integer :: ibas
    Type (zcsrmat) :: zcsr_2sp, zcsr_sum, zcsr_tmp

    call blockdiag(geo%atoms(1)%ptr%zcsr_nonsph(ind), geo%atoms(1)%ptr%zcsr_nonsph(ind), zcsr_sum)
    Do ibas=2,geo%num
     call blockdiag(geo%atoms(ibas)%ptr%zcsr_nonsph(ind), geo%atoms(ibas)%ptr%zcsr_nonsph(ind), zcsr_2sp)
     call blockdiag(zcsr_sum, zcsr_2sp, zcsr_tmp)
     call spcopy(zcsr_tmp, zcsr_sum)
    Enddo 

    call spcopy(zcsr_sum, zcsr_out)
    
    call free(zcsr_2sp)
    call free(zcsr_sum)
    call free(zcsr_tmp)

   End Subroutine mk_V

      Subroutine make_sk_diag (sk, scr, kpar, needhops, cfs)
         Use sparselib
         Implicit None
         Type (t_strconst) :: scr
         Type (t_mathop) :: sk
         Real (Kind=DEF_DBL_PREC), Optional :: cfs
         Integer, Optional :: needhops
         Real (Kind=DEF_DBL_PREC) :: kpar (2)
!!$ local
         Integer :: nh

         If ( .Not. present(needhops)) Then
            nh = 0
         Else
            nh = needhops
         End If

         Call fillsk (sk%c, scr%nnz, scr%nrows, scr%nm, scr%main, scr%base, kpar, cfs)
         If (nh > 0) Then
            Call fillsk (sk%l, scr%lnnz, scr%lnrows, scr%nl, scr%lhop, scr%base, kpar, cfs)
            Call fillsk (sk%r, scr%rnnz, scr%rnrows, scr%nr, scr%rhop, scr%base, kpar, cfs)
            sk%nl = sk%l%ncol
            sk%nr = sk%r%ncol
            sk%havehops = 1
            sk%r%a = sk%r%a * Exp (DEF_cmplx_Ione*Dot_product(kpar, scr%rtr))
            sk%l%a = sk%l%a * Exp (DEF_cmplx_Ione*Dot_product(kpar,-scr%ltr))
         End If
         sk%alloc = 1
         sk%n = sk%c%ncol

!!$ Force to be hermitian
!!$    t_s=spherm(sk%c)
!!$    write(*,*) 'maxdev=',maxval(abs(sk%c%a-t_s%a))
!!$    sk%c%a=0.5d0*(sk%c%a+t_s%a)
!!$    call free(t_s)
      Contains

         Subroutine fillsk (sk, tnnzi, nrow, na, sites, base, kpar, cfs)
!!$ function which is actually allocate and fill sk
            Implicit None
            Integer :: tnnzi, nrow, na
            Type (zcsrmat), Target :: sk
            Type (t_str_site), Target :: sites (:)
            Real (Kind=DEF_DBL_PREC) :: base (2, 2)
            Real (Kind=DEF_DBL_PREC) :: kpar (2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs
!!$ local
            Type (t_str_site), Pointer :: jsite, isite
            Integer :: ia, nc, nr, i, nnz, j, lr, tnnz
            Integer, Pointer :: jc (:), ir (:)
            Complex (Kind=DEF_DBL_PREC), Pointer :: val (:)


            Call alloc (sk, tnnzi, nrow)

            jc => sk%jc (1:)
            ir => sk%ir (1:)

            val => sk%a (1:)
            tnnz = 0
            lr = 1
!!$   !$omp parallel do private(isite,nc,nr,j,) reduction(+:tnnz)
            Do ia = 1, na
               isite => sites (ia)
               nc = isite%nrows
               nr = isite%nrows
               j = isite%srow
               nnz = fillskrows (isite, kpar, val, base, cfs, j)
               ir (lr:j) = tnnz + 1
               Do i = 1, nr
                  ir (j+1) = ir (j) + nc
                  j = j + 1
                  jc (1:nc) = [( i,i=isite%srow,isite%srow+nc )] 
                  jc => jc (nc+1:)
               End Do
               lr = j
               val => val (nnz+1:)
               tnnz = tnnz + nnz
            End Do
!!$  !$omp parallel end do
            sk%nnz = tnnz

         End Subroutine fillsk

         Function fillskrows (site, k, val, base, cfs_in, j) Result (nnz)
!!$ function which fill rows of sk which correspondent to "site"
            Implicit None
            Type (t_str_site), Target :: site
            Real (Kind=DEF_DBL_PREC) :: k (2), base (2, 2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs_in
            Complex (Kind=DEF_DBL_PREC) :: val (:)
            Integer :: nnz, j
!!$ local
            Integer :: in, jnc, ir, jnc1
            Complex (Kind=DEF_DBL_PREC) :: els (site%nrows, site%ncols), fact
            Integer, Pointer :: tp (:, :)
            Type (t_str_block), Pointer :: hop (:)
            Real (Kind=DEF_DBL_PREC) :: cfs !,b(2)

!!$             b(1)=0.0d0
!!$             b(2)=0.4082482904639
            cfs = 1.0d0
            If (present(cfs_in)) cfs = cfs_in

            tp => site%trpar
            hop => site%hop
!!$             write(*,*) '----', base(:,1)
            els = DEF_cmplx_zero
!!$ !$omp parallel do private(jnc,jnc1,fact) shared(cfs,k,base)
            Do in = 1, site%nhop
               jnc = site%startn (in)
               if (jnc == j) then
                jnc1 = jnc + hop(in)%ncol - 1
!!$                 write(*,*) k*(tp(1, in)*base(:, 1)+tp(2, in)*base(:, 2)+b)
                fact = cfs * Exp (DEF_cmplx_Ione*Dot_product(k, dble(tp(1, in))*base(:, 1)+dble(tp(2, &
              &  in))*base(:, 2)))
                els (:, jnc:jnc1) = els (:, jnc:jnc1) + fact * hop(in)%bl
               endif
            End Do
!!$ !$omp end parallel do
            in = site%nrows
!!$ !$omp parallel do private(jnc,jnc1) firstprivate(in)
            Do ir = 1, site%nrows
               jnc = (ir-1) * in + 1
               jnc1 = jnc + in - 1
               val (jnc:jnc1) = els (ir, :)
            End Do
!!$ !$omp end parallel do
            nnz = site%nrows * in
         End Function fillskrows
      End Subroutine make_sk_diag

 End Module
