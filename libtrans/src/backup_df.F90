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

 do j=1,zcsr_in%nrow
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

 Function downfold(zcsr_in,to_downfold,en) result (zcsr_out)
 Implicit None
 type (zcsrmat) :: zcsr_in
 type (zcsrmat) ::  zcsr_out !!$ matrix comes in, matrix goes out. Can't explain that.
 integer :: to_downfold !!$ This is the index of the to-downfold state 
 real (kind=DEF_DBL_PREC) :: en
 !!$ Local
 type (zcsrmat) :: AA, AB, BA
 complex (kind=DEF_DBL_PREC) :: BB
 integer :: i,j, aa_nnz, ab_nnz, ba_nnz

 call alloc(AA,zcsr_in%nnz,zcsr_in%nrow-1)
 call alloc(AB,zcsr_in%ncol,zcsr_in%nrow-1)
 call alloc(BA,zcsr_in%ncol,zcsr_in%nrow-1)

 aa_nnz=1
 ab_nnz=1
 ba_nnz=1

 do j=1,zcsr_in%nrow
  if (j==to_downfold) then
   do i=zcsr_in%ir(j),zcsr_in%ir(j+1)-1

    if (zcsr_in%jc(i) == to_downfold) then
     BB = zcsr_in%a(i)
    else
     if (abs(zcsr_in%a(i)) > 1E-8) then
      if (zcsr_in%jc(i) < to_downfold) then
       BA%a(ba_nnz) = zcsr_in%a(i)
       BA%jc(ba_nnz) = zcsr_in%jc(i)
       BA%nnz = ba_nnz
       ba_nnz = ba_nnz+1
      else
       BA%a(ba_nnz) = zcsr_in%a(i)
       BA%jc(ba_nnz) = zcsr_in%jc(i) - 1
       BA%nnz = ba_nnz
       ba_nnz = ba_nnz+1
      endif
     endif
    endif

   enddo
   BA%ir(2) = ba_nnz 
  else
   do i=zcsr_in%ir(j),zcsr_in%ir(j+1)-1

    if (zcsr_in%jc(i) == to_downfold) then
     if (abs(zcsr_in%a(i)) > 1E-8) then
      AB%a(ab_nnz) = zcsr_in%a(i)
      AB%jc(ab_nnz) = 1 !!$zcsr_in%jc(i)
      AB%nnz = ab_nnz
      ab_nnz = ab_nnz+1
     endif
    else
     if (abs(zcsr_in%a(i)) > 1E-8) then
      if (zcsr_in%jc(i) < to_downfold) then
       AA%a(aa_nnz) = zcsr_in%a(i) 
       AA%jc(aa_nnz) = zcsr_in%jc(i)
       AA%nnz = aa_nnz
       aa_nnz = aa_nnz+1
      else
       AA%a(aa_nnz) = zcsr_in%a(i)
       AA%jc(aa_nnz) = zcsr_in%jc(i)-1
       AA%nnz = aa_nnz
       aa_nnz = aa_nnz+1
      end if
     end if
    endif
   enddo
   if (j < to_downfold) then
    AA%ir(j+1) = aa_nnz
    !!$BA%ir(j+1) = ba_nnz 
    AB%ir(j+1) = ab_nnz
   else
    AA%ir(j)   = aa_nnz 
    !!$BA%ir(j+1) = ba_nnz
    AB%ir(j)   = ab_nnz
   endif
  endif
 enddo

 AB%a = AB%a/(en-BB)


 zcsr_out = spmatadd(AA,spmatmul(AB,BA)) 
 
 call free(AA)
 call free(AB)
 call free(BA)

 End Function

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

 to_d = to_downfold

 do i=1,size(to_d)

  zcsr_out = downfold(zcsr_in,to_d(i),en)
  call delete_nearzeros(zcsr_out,1D-10)
  call spcopy(zcsr_out,zcsr_in) 
 
  to_d = to_d - 1
 enddo

 !!$call alloc(zcsr_out,zcsr_in%nnz,zcsr_in%ncol-size(to_d))
 !!$zcsr_out%a(:) = zcsr_in%a(:)
 !!$zcsr_out%jc(:) = zcsr_in%jc(:)
 !!$zcsr_out%ir(:) = zcsr_in%ir(:)

 !!$call dump_matrix(zcsr_out,'zcsr_finalout')

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
