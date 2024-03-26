!!$ Module which contains the subroutines related to the non-spherical corrections.
#include "math_def.h"

Module df
 Use sparselib
 Implicit None

 Contains

 Subroutine delete_nearzeros(zcsr_in,tol)
 Implicit none
 type (zcsrmat) :: zcsr_in
 real (kind=DEF_DBL_PREC) :: tol
 !!$ local
 integer :: i,j, temp_nnz
 type (zcsrmat) :: zcsr_temp

 call alloc(zcsr_temp,zcsr_in%nnz,zcsr_in%ncol)

 temp_nnz = 1

 do j=1,zcsr_in%ncol
  do i = zcsr_in%ir(j), zcsr_in%ir(j+1)-1
   if (abs(zcsr_in%a(i)) > tol) then
    zcsr_temp%a(temp_nnz) = zcsr_in%a(i)
    zcsr_temp%jc(temp_nnz) = zcsr_in%jc(i)
    zcsr_temp%nnz = temp_nnz
    temp_nnz = temp_nnz + 1
   endif        
  enddo
  zcsr_temp%ir(j+1) = temp_nnz
 enddo

 call free(zcsr_in)

 call alloc(zcsr_in,temp_nnz-1,zcsr_temp%ncol)
 zcsr_in%a = zcsr_temp%a
 zcsr_in%jc = zcsr_temp%jc
 zcsr_in%ir = zcsr_temp%ir
 zcsr_in%nnz = zcsr_temp%nnz

 call free(zcsr_temp)

 End subroutine

 Subroutine extract_upper_lower(n,m,ia,ja,a,ul_ir,ul_jc,ul_a,ul_nnz,lu_ir,lu_jc,lu_a,lu_nnz,uu_ir,uu_jc,uu_a,uu_nnz,to_downfold,en)
 Implicit None
 integer, intent(in) :: n,m
 integer, intent(in) :: ia(:), ja(:)
 complex (kind=DEF_DBL_PREC), intent(in) :: a(:)
 integer :: to_downfold !!$ This is the index of the to-downfold state 
 real (kind=DEF_DBL_PREC) :: en
 integer, intent(out) :: uu_jc(:), ul_jc(:), lu_jc(:)
 integer, intent(out) :: uu_ir(:), ul_ir(:), lu_ir(:)
 complex (kind=DEF_DBL_PREC), intent(out) :: uu_a(:), ul_a(:), lu_a(:)
!!$ local
 integer :: i,j, uu_nnz, ul_nnz, lu_nnz
 complex (kind=DEF_DBL_PREC) :: ll

 uu_nnz=1
 ul_nnz=1
 lu_nnz=1

 do j=1,n
  if (j==to_downfold) then
   do i=ia(j),ia(j+1)-1

    if (ja(i) == to_downfold) then
     ll = a(i)
    else
     if (abs(a(i)) > 1E-8) then
      if (ja(i) < to_downfold) then
       lu_a(lu_nnz) = a(i)
       lu_jc(lu_nnz) = ja(i)
       lu_nnz = lu_nnz+1
      else
       lu_a(lu_nnz) = a(i)
       lu_jc(lu_nnz) = ja(i) - 1
       lu_nnz = lu_nnz+1
      endif
     endif
    endif

   enddo
   lu_ir(2) = lu_nnz 
  else
   do i=ia(j),ia(j+1)-1

    if (ja(i) == to_downfold) then
     if (abs(a(i)) > 1E-8) then
      ul_a(ul_nnz) = a(i)
      ul_jc(ul_nnz) = 1 !!$ja(i)
      ul_nnz = ul_nnz+1
     endif
    else
     if (abs(a(i)) > 1E-8) then
      if (ja(i) < to_downfold) then
       uu_a(uu_nnz) = a(i) 
       uu_jc(uu_nnz) = ja(i)
      else
       uu_a(uu_nnz) = a(i)
       uu_jc(uu_nnz) = ja(i)-1
       uu_nnz = uu_nnz+1
      end if
     end if
    endif
   enddo
   if (j < to_downfold) then
    uu_ir(j+1) = uu_nnz
    !!$lu_ir(j+1) = lu_nnz 
    ul_ir(j+1) = ul_nnz
   else
    uu_ir(j)   = uu_nnz 
    !!$lu_ir(j+1) = lu_nnz
    ul_ir(j)   = ul_nnz
   endif
  endif
 enddo

 ul_a = ul_a/(en-ll)

 End Subroutine 
 
 Subroutine matmul_plus (n, m, ia, ja, a, ib, jb, b, ic, jc, c, id, jd, d)
 Implicit none
 integer, intent(in) :: n, m
 integer, intent(in) :: ia(:), ib(:), ic(:) 
 integer, intent(in) :: ja(:), jb(:), jc(:) 
 Complex (Kind=DEF_DBL_PREC), Intent(In) :: a(:), b(:), c(:)
 Integer, intent(out) :: id(:), jd(:)
 Complex (Kind=DEF_DBL_PREC), Intent(Out) :: d(:)
 !! Local
 integer :: ie(n+1), je(size(a)*size(b))
 Complex (Kind=DEF_DBL_PREC) :: e(size(a)*size(b))


 call matmul_normal_zz (n, m, ia, ja, a, ib, jb, b, ie, je, e)
 
 call aplb1_zz(n, m, c, jc, ic, e, je, ie, d, jd, id)

 End Subroutine matmul_plus 



 Subroutine downfolding(zcsr_in,zcsr_out,to_downfold,en)
 Implicit none
 !!$ passed
 type (zcsrmat) :: zcsr_in
 type (zcsrmat) ::  zcsr_out !!$ matrix comes in, matrix goes out. Can't explain that.
 integer, dimension(:) :: to_downfold(:) !!$ This is a vector with the indices of the states that are downfolded, 
 real (kind=DEF_DBL_PREC) :: en
 !!$ local
 integer :: i
 character(len=15) :: str
 integer :: to_d(size(to_downfold)) !!$ copy of to_downfold
 integer, allocatable :: AB_ir(:), BA_ir(:), AA_ir(:)
 integer, allocatable :: AB_jc(:), BA_jc(:), AA_jc(:)
 complex (kind=DEF_DBL_PREC), allocatable :: AB_a(:), BA_a(:), AA_a(:)
 integer :: AB_nnz, BA_nnz, AA_nnz

 to_d = to_downfold

 allocate(AB_ir(zcsr_in%nrow),BA_ir(2),AA_ir(zcsr_in%nrow))
 allocate(AB_jc(zcsr_in%ncol-1),BA_jc(zcsr_in%ncol-1),AA_jc(zcsr_in%nnz))
 allocate(AB_a(zcsr_in%nrow-1),BA_a(zcsr_in%ncol-1),AA_a(zcsr_in%nnz))

 call alloc(zcsr_out, zcsr_in%nnz, zcsr_in%ncol)

 do i=1,size(to_d)
 
  call extract_upper_lower(zcsr_in%nrow, zcsr_in%ncol, zcsr_in%ir, zcsr_in%jc, zcsr_in%a, &
 & AB_ir, AB_jc, AB_a, AB_nnz, BA_ir, BA_jc, BA_a, BA_nnz, AA_ir, AA_jc, AA_a, AA_nnz, to_d(i), en)
  call matmul_plus(zcsr_in%nrow-1, zcsr_in%ncol-1, AB_ir(1:zcsr_in%nrow), AB_jc(1:AB_nnz), AB_a(1:AB_nnz), BA_ir(1:2),& 
 & BA_jc(1:zcsr_in%nrow-1), BA_a(1:BA_nnz), AA_ir(1:zcsr_in%nrow), AA_jc(1:AA_nnz), AA_a(1:AA_nnz), zcsr_out%ir, zcsr_out%jc, zcsr_out%a, zcsr_out%nnz)

  zcsr_in%nrow = zcsr_in%nrow-1 
  zcsr_in%ncol = zcsr_in%ncol-1 
  zcsr_in%a(1:zcsr_out%nnz) = zcsr_out%a(1:zcsr_out%nnz)
  zcsr_in%jc(1:zcsr_out%nnz) = zcsr_out%jc(1:zcsr_out%nnz)
  zcr_in%ir(1:zcsr_in%nrow) = zcsr_out%ir(1:zcsr_in%nrow)

  to_d = to_d - 1
 enddo

 deallocate(AB_ir, BA_ir, AA_ir)
 deallocate(AB_jc, BA_jc, AA_jc)
 deallocate(AB_a, BA_a, AA_a)

 call free(zcsr_in)

 End Subroutine

 Subroutine offdiag(amat, bmat, cmat)
  Implicit None
  Type (zcsrmat) :: amat, bmat
  Type (zcsrmat) :: cmat
  !!$ Local

  if (cmat%alloc==1) call free(cmat)
  call alloc (cmat, amat%nnz + bmat%nnz, amat%nrow + bmat%nrow, amat%ncol + bmat%ncol)

  cmat%nnz = amat%nnz+bmat%nnz
  cmat%ncol = amat%ncol+bmat%ncol
  cmat%nrow = amat%nrow+bmat%nrow

  cmat%a = (/ amat%a, bmat%a /)
  cmat%ir = (/ amat%ir, bmat%ir(2:size(bmat%ir))+amat%nnz /)
  cmat%jc = (/ amat%jc+amat%ncol, bmat%jc /)

 End Subroutine offdiag

 Subroutine extract(amat,bmat,rows,cols,gap_in_col,clear)
  Implicit None
  !!$ Passed
  Type (zcsrmat), intent(in) :: amat
  Type (zcsrmat) :: bmat
  Integer, dimension(:) :: rows(:)
  Integer, dimension(:) :: cols(:)
  Integer, Optional :: gap_in_col
  Integer, Optional :: clear
  !!$ Local
  Integer :: i, j, bmat_nnz, ip


  if (clear>0) then
   if (bmat%alloc==1) call free(bmat)
   call alloc(bmat,amat%nnz,size(rows),size(cols))
  else
   if (bmat%alloc==0) call alloc(bmat,amat%nnz,size(rows),size(cols))
  endif

  bmat_nnz=1

  if (present(gap_in_col)) then

   do ip=1,size(rows)
    i=rows(ip)
    do j=amat%ir(i),amat%ir(i+1)-1
     if (any(cols==amat%jc(j))) then
      bmat%a(bmat_nnz) = amat%a(j)
      if (amat%jc(j) < gap_in_col) bmat%jc(bmat_nnz) = amat%jc(j)
      if (amat%jc(j) > gap_in_col) bmat%jc(bmat_nnz) = amat%jc(j) - 1 
      bmat%nnz = bmat_nnz
      bmat_nnz = bmat_nnz+1
     endif
    enddo
    bmat%ir(ip+1) = bmat_nnz
   enddo

  else

   do ip=1,size(rows)
    i=rows(ip) 
    do j=amat%ir(i),amat%ir(i+1)-1
     if (any(cols==amat%jc(j))) then
      bmat%a(bmat_nnz) = amat%a(j)
      bmat%jc(bmat_nnz) = amat%jc(j) - cols(1) + 1
      bmat%nnz = bmat_nnz
      bmat_nnz = bmat_nnz+1
     endif
    enddo
    bmat%ir(ip+1) = bmat_nnz
   enddo       

  endif

!!$ call ordercsr(bmat)

 End Subroutine extract



end Module df 
