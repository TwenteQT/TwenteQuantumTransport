!!$ Module which contains the subroutines related to the non-spherical corrections.
#include "math_def.h"

Module df
 Use sparselib
 Use potpars
 Use structure
 Use geometry_module
 Implicit None

 Contains

 Subroutine downfold_kink_scattering(sys, geo, en)
 Implicit None
 Type (t_mathop)        :: sys
 Type (t_geometry_EMTO) :: geo
 Real (Kind=DEF_DBL_PREC), Optional :: en
 !!$ local
 Type (zcsrmat)         :: zcsr_temp, zcsr_df2
 Integer                :: i

 if (present(en)) then
  call downfolding_alt(sys%c, zcsr_temp, geo%to_downfold,en)
 else
  call downfolding_alt(sys%c, zcsr_temp, geo%to_downfold,0d0)
 endif

 call extract(zcsr_temp, zcsr_df2,[(i,i=geo%num_emb_l+1,zcsr_temp%nrow-geo%num_emb_r )],&
  &[(i,i=geo%num_emb_l+1,zcsr_temp%ncol-geo%num_emb_r)],clear=1)

 call spcopy(zcsr_df2,sys%c)

 !call spcopy(zcsr_temp, sys%c)
 call free(zcsr_temp)
 call free(zcsr_df2)

 End Subroutine

 Subroutine downfold_kink_lead(sys, geo, n, ndim, en)
 Implicit None
 Type (t_mathop)        :: sys
 Type (t_geometry_EMTO) :: geo
 Integer                :: n, ndim
 Real (Kind=DEF_DBL_PREC), Optional :: en
 !!$ local
 Type (zcsrmat)         :: zcsr_temp, zcsr2_temp, zcsr_temp2, zcsr2_temp2, zcsr_temp3
 Type (zcsrmat)         :: zcsr_null, zcsr_off, zcsr_off2, zcsr_off3
 Integer                :: i
 Integer                :: to_df_cr(4*size(geo%to_downfold))


  !!$ Construct {{sys%c, sys%r},{sys%r^*, sys%c}} matrix
  !!$ for downfolding the interactions properly
  call blockdiag(sys%c, sys%c, zcsr_temp)
  call blockdiag(zcsr_temp, zcsr_temp, zcsr2_temp) !!!!$

  !!$ offdiagonal component, sys%r and its hermitian conjugate
  sys%l = spherm(sys%r)
  call offdiag(sys%r,sys%l,zcsr_temp2)

  call blockdiag(zcsr_temp2, zcsr_temp2, zcsr2_temp2) !!!!$

  call spcopy(sys%r, zcsr_null)
  zcsr_null%a = 0d0

  call offdiag(zcsr_null,sys%r,zcsr_off)
  call offdiag(sys%l,zcsr_null,zcsr_off2)
  call offdiag(zcsr_off,zcsr_off2,zcsr_off3)

  zcsr_temp2 = spmatadd(zcsr2_temp2,zcsr_off3)

  !!!$zcsr_temp3 = spmatadd(zcsr_temp, zcsr_temp2)
  zcsr_temp3 = spmatadd(zcsr2_temp, zcsr_temp2) !!!!$

  call free(zcsr_off)
  call free(zcsr_off2)
  call free(zcsr_off3)


  to_df_cr = (/ geo%to_downfold, geo%to_downfold+ndim, &
             & geo%to_downfold+2*ndim, geo%to_downfold+3*ndim /)

  call free(zcsr_temp)
  if (present(en)) then
  call downfolding_alt(zcsr_temp3,zcsr_temp,to_df_cr,en)
  else
  call downfolding_alt(zcsr_temp3,zcsr_temp,to_df_cr,0d0)
  endif

  call free_mathop(sys)

  !!$ extract new sys%c and sys%r matrices from downfolded matrix
  call extract(zcsr_temp, sys%c,[(i,i=n+1,2*n )],[(i,i=n+1  ,2*n)] ,clear=1)
  call extract(zcsr_temp, sys%r,[(i,i=n+1,2*n )],[(i,i=2*n+1,3*n)] ,clear=1)

  call free(zcsr_temp)
  call free(zcsr2_temp)
  call free(zcsr_temp2)
  call free(zcsr2_temp2)
  call free(zcsr_temp3)
 End Subroutine

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

 call spcopy(zcsr_temp,zcsr_in)

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
 
 ul_ir(1) = 1
 lu_ir(1) = 1
 uu_ir(1) = 1 

 do j=1,n
  if (j==to_downfold) then
   do i=ia(j),ia(j+1)-1

    if (ja(i) == to_downfold) then
     ll = a(i)
    else
!!$     if (abs(a(i)) > 1D-10) then
      if (ja(i) < to_downfold) then
       lu_a(lu_nnz) = a(i)
       lu_jc(lu_nnz) = ja(i)
       lu_nnz = lu_nnz+1
      else
       lu_a(lu_nnz) = a(i)
       lu_jc(lu_nnz) = ja(i) - 1
       lu_nnz = lu_nnz+1
      endif
!!$     endif
    endif

   enddo
   lu_ir(2) = lu_nnz 
  else
   do i=ia(j),ia(j+1)-1

    if (ja(i) == to_downfold) then
!!$     if (abs(a(i)) > 1D-10) then
      ul_a(ul_nnz) = a(i)
      ul_jc(ul_nnz) = 1 !!$ja(i)
      ul_nnz = ul_nnz+1
!!$     endif
    else
!!$     if (abs(a(i)) > 1D-10) then
      if (ja(i) < to_downfold) then
       uu_a(uu_nnz) = a(i) 
       uu_jc(uu_nnz) = ja(i)
       uu_nnz = uu_nnz+1
      else
       uu_a(uu_nnz) = a(i)
       uu_jc(uu_nnz) = ja(i)-1
       uu_nnz = uu_nnz+1
      end if
!!$     end if
    endif
   enddo
   if (j < to_downfold) then
    uu_ir(j+1) = uu_nnz
    ul_ir(j+1) = ul_nnz
   else
    uu_ir(j)   = uu_nnz 
    ul_ir(j)   = ul_nnz
   endif
  endif
 enddo

 ul_a = ul_a/(en-ll)

 End Subroutine 
 
 Subroutine matmul_plus (n, m, ia, ja, a, ib, jb, b, ic, jc, c, id, jd, d, d_nnz)
 Implicit none
 integer, intent(in) :: n, m
 integer, intent(in) :: ia(:), ib(:), ic(:) 
 integer, intent(in) :: ja(:), jb(:), jc(:) 
 Complex (Kind=DEF_DBL_PREC), Intent(In) :: a(:), b(:), c(:)
 Integer, intent(out) :: id(n+1), jd(:), d_nnz
 Complex (Kind=DEF_DBL_PREC), Intent(Out) :: d(:)
 !! Local
 integer :: ie(n+1)
 !integer, allocatable :: je(:)
 !Complex (Kind=DEF_DBL_PREC), allocatable :: e(:)
 integer :: je(n*n)
 Complex (Kind=DEF_DBL_PREC) :: e(n*n)
 integer :: nz
 integer :: i,j

 ie = 0
 je = 0
 e = 0d0
 id = 0
 jd = 0
 d = 0d0

 call spmatmul_size(n, m, ia, ja, ib, jb, nz, ie)
 !allocate(je(nz),e(nz))

 call matmul_normal_zz (n, m, ia, ja, a, ib, jb, b, ie, je(1:nz), e(1:nz))

 call z_ordercsr_nonptr (n, ie, je, e)

 call matadd_size(n, m, jc, ic, je, ie, nz, id) 
 
 call aplb1_zz(n, m, c, jc, ic, e, je, ie, d(1:nz), jd(1:nz), id)

 !deallocate(je,e)

 d_nnz = nz

 End Subroutine matmul_plus 

      Subroutine z_ordercsr_nonptr (n, ia, ja, a)
         Implicit None
         Integer, Intent(In) :: n
         Integer, Intent(in), target :: ia(:), ja(:)
         Complex (kind=DEF_DBL_PREC), intent(in), target :: a(:)
         Integer :: i, k, k1, kn, kn0
         Integer, Pointer :: iwork (:)
         Complex (Kind=DEF_DBL_PREC), Pointer :: aa (:)
         Integer, Pointer :: jc (:)

         kn0 = 0
         Do i = 1, n
            k = ia (i+1) - ia(i)
            If (k > kn0) kn0 = k
         End Do

!$omp parallel private(iwork,k,k1,aa,jc,kn) shared(kn0)
         Allocate (iwork(kn0))
!$omp do
         Do i = 1, n
            k = ia (i)
            k1 = ia (i+1) - 1
            kn = k1 - k + 1
            aa => a (k:k1)
            jc => ja (k:k1)
            Call iqsort (kn, jc, iwork)
            aa = aa (iwork(:kn))
            jc = jc (iwork(:kn))
         End Do
!$omp end do
         Deallocate (iwork)
!$omp end parallel

      End Subroutine z_ordercsr_nonptr

 Subroutine downfolding_alt(zcsr_in,zcsr_out,to_downfold,en)
 Implicit None
 Type (zcsrmat) :: zcsr_in !!$ this will be K, the Kink matrix
 Type (zcsrmat), Intent(Out) :: zcsr_out
 Integer, Dimension(:) :: to_downfold(:) !!$ This is a vector with the indices of the states that are downfolded, 
 Real (Kind=DEF_DBL_PREC) :: en
 !!$ local 
 Integer :: to_d(size(to_downfold)) !!$ copy of to_downfold
 Integer :: i, j, u_nz, ndf, idx
 Type (zcsrmat) :: I_u, I_l !!$ partial identity matrices 
 Type (zcsrmat) :: K_uu, K_ul, K_lu !!$upperupper,  upperlower and lowerupper parts of K 
 Type (zcsrmat) :: K_temp, zcsr_temp 
 Real (Kind=DEF_DBL_PREC) :: K_ll
 Integer :: non_df(zcsr_in%nrow)

 call delete_nearzeros(zcsr_in,1D-15)

 to_d = to_downfold
 ndf = size(to_downfold)

 Allocate(I_u%a(zcsr_in%nrow), I_u%ir(zcsr_in%nrow+1), I_u%jc(zcsr_in%nrow))
 I_u%alloc = 1
 I_u%nrow = zcsr_in%nrow
 I_u%ncol = zcsr_in%ncol
 I_u%ir(1) = 1

 Allocate(I_l%a(zcsr_in%nrow), I_l%ir(zcsr_in%nrow+1), I_l%jc(zcsr_in%nrow))
 I_l%alloc = 1
 I_l%nrow = zcsr_in%nrow
 I_l%ncol = zcsr_in%ncol
 I_l%ir(1) = 1


 do i=1,ndf

  u_nz = 1

  I_u%a = 0
  I_u%ir(1) = 1

  I_l%a = 0
  I_l%ir(1) = 1
  !!$ Generate Id matrices I_u and I_l of dimension zcsr_in%nrow
  do j=1,zcsr_in%nrow
   if (j .ne. to_d(i)) then
    I_u%a(u_nz) = 1
    I_u%jc(u_nz) = j
    u_nz = u_nz + 1

    I_u%ir(j+1) = I_u%ir(j) + 1
    I_l%ir(j+1) = I_l%ir(j)
   else
    I_l%a(1) = 1
    I_l%jc(1) = j

    I_u%ir(j+1) = I_u%ir(j)
    I_l%ir(j+1) = I_l%ir(j) + 1
   endif
  enddo

  I_u%nnz = u_nz - 1 
  I_l%nnz = 1

  !!$ Get K_ul = I_u * K * I_l
  K_ul = get3MMM(I_u, zcsr_in, I_l)
  !!$ and K_lu = I_l * K * I_u
  K_lu = get3MMM(I_l, zcsr_in, I_u)

  !!$ and K_uu
  K_uu = get3MMM(I_u, zcsr_in, I_u)

  !!$ extract K_ll from K
  do j=zcsr_in%ir(to_d(i)), zcsr_in%ir(to_d(i)+1)-1
   if (zcsr_in%jc(j) == to_d(i)) K_ll = zcsr_in%a(j)
  enddo

  !!$ evaluate K_df = K - K_ul * 1/K_ll * K_lu
  K_temp = spmatmul(K_ul,K_lu)
  K_temp%a = K_temp%a/K_ll

  zcsr_temp = spmatadd(K_uu,K_temp,sign=-1)
  call spcopy(zcsr_temp,zcsr_in)

 enddo

 Deallocate(I_u%a, I_u%ir, I_u%jc)
 Deallocate(I_l%a, I_l%ir, I_l%jc)
 I_u%alloc = 0
 I_l%alloc = 0

 idx = 1
 do i=1,zcsr_in%nrow
  if (any(to_d==i)) cycle
  non_df(idx) = i
  idx = idx + 1 
 enddo

 call delete_nearzeros(zcsr_in,1D-15)
 call extract(zcsr_in,zcsr_out,non_df(1:idx-1),non_df(1:idx-1),clear=1)

 call free(zcsr_temp)
 call free(K_uu)
 call free(K_ul)
 call free(K_lu)
 call free(K_temp)
 End subroutine


 Subroutine downfolding(zcsr_in,zcsr_out,to_downfold,en)
 Implicit none
 !!$ passed
 type (zcsrmat), Intent(In) :: zcsr_in
 type (zcsrmat), Intent(Out) ::  zcsr_out !!$ matrix comes in, matrix goes out. Can't explain that.
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
 type (zcsrmat) :: zcsr_work

 to_d = to_downfold

 allocate(AB_ir(zcsr_in%nrow),AB_jc(zcsr_in%ncol-1),AB_a(zcsr_in%nrow-1))
 allocate(BA_ir(2),BA_jc(zcsr_in%ncol-1),BA_a(zcsr_in%ncol-1))
 allocate(AA_ir(zcsr_in%nrow),AA_jc(30*zcsr_in%nnz),AA_a(30*zcsr_in%nnz))


 call alloc(zcsr_work,30*zcsr_in%nnz,zcsr_in%nrow) !!$ fixme this is terrible, a
!!$ magic number to ensure enough memory allocation....

 zcsr_work%a(1:zcsr_in%nnz) = zcsr_in%a(1:zcsr_in%nnz)
 zcsr_work%jc(1:zcsr_in%nnz) = zcsr_in%jc(1:zcsr_in%nnz)
 zcsr_work%ir(1:zcsr_in%nrow+1) = zcsr_in%ir(1:zcsr_in%nrow+1)
 zcsr_work%nnz = zcsr_in%nnz

 do i=1, size(to_d)
 
  call extract_upper_lower(zcsr_work%nrow, zcsr_work%ncol, &
&zcsr_work%ir(1:zcsr_work%nrow+1), zcsr_work%jc(1:zcsr_work%nnz), zcsr_work%a(1:zcsr_work%nnz), &
 & AB_ir, AB_jc, AB_a, AB_nnz, BA_ir, BA_jc, BA_a, BA_nnz, AA_ir, AA_jc, AA_a, AA_nnz, to_d(i), en)
 
  call matmul_plus(zcsr_work%nrow-1, zcsr_work%ncol-1, AB_ir(1:zcsr_work%nrow), AB_jc(1:AB_nnz), &
&AB_a(1:AB_nnz), BA_ir(1:2), BA_jc(1:BA_nnz), BA_a(1:BA_nnz), AA_ir(1:zcsr_work%nrow), &
&AA_jc(1:AA_nnz), AA_a(1:AA_nnz), zcsr_work%ir(1:zcsr_work%nrow), zcsr_work%jc, zcsr_work%a, zcsr_work%nnz)


  zcsr_work%nrow = zcsr_work%nrow-1 
  zcsr_work%ncol = zcsr_work%ncol-1 
  call delete_nearzeros(zcsr_work,1D-10)

  to_d = to_d - 1
 enddo

 
 !!$call delete_nearzeros(zcsr_work,1D-10)
 call spcopy(zcsr_work,zcsr_out)
 call zcsrcompress(zcsr_out,1D-10)

 deallocate(AB_ir, AB_jc, AB_a)
 deallocate(BA_ir, BA_jc, BA_a)
 deallocate(AA_ir, AA_jc, AA_a)

!!$ write(*,*) 'some success! DF almost worked'
 call free(zcsr_work)
!!$ call free(zcsr_in)

 !!$write(*,*) 'great success! DF worked'

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
  Type (zcsrmat) :: amat
  Type (zcsrmat) :: bmat
  Integer :: rows(:)
  Integer :: cols(:)
  Integer, Optional :: gap_in_col
  Integer, Optional :: clear
  !!$ Local
  Integer :: i, j, bmat_nnz, ip, new_col


  if (clear>0) then
   if (bmat%alloc==1) call free(bmat)
   call alloc(bmat,amat%nnz,size(rows),size(cols))
  else
   if (bmat%alloc==0) call alloc(bmat,amat%nnz,size(rows),size(cols))
  endif

  bmat_nnz=1

  do ip=1,size(rows)
   i=rows(ip) 
   do j=amat%ir(i),amat%ir(i+1)-1
    new_col = Findloc(cols,amat%jc(j), 1)
    if (new_col>0) then
     bmat%a(bmat_nnz) = amat%a(j)
     bmat%jc(bmat_nnz) = new_col 
     bmat%nnz = bmat_nnz
     bmat_nnz = bmat_nnz+1
    endif
   enddo
   bmat%ir(ip+1) = bmat_nnz
  enddo       

 call ordercsr(bmat)

 End Subroutine extract



end Module df 
