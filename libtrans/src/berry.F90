
#include "math_def.h"
Module berry
 Implicit None

 Enum, bind(c)
  Enumerator :: ANOMALOUS=1, SPIN=2, ORBITAL=3
 Endenum

Contains

Function compute_HC(bz, cond) Result (HC)
 use adaptive_bzgrid
 Implicit None
 Type (adapt_bzone) :: bz
 Integer :: cond
 Real (Kind=DEF_DBL_PREC) :: HC(3)
 !local
 Integer :: ik

 HC = 0d0

 Open(file='fermiS_up.stl',unit=919)
 Open(file='fermiS_down.stl',unit=920)
 Write(919,*) 'solid Fermi surface'
 Write(920,*) 'solid Fermi surface'

 Do ik = 1,bz%nk
  If (bz%bz(ik)%bp > 0) Cycle
  HC = HC + calc_HC_bzp(bz, ik, cond)
 Enddo

 Write(919,*) 'endsolid Fermi surface'
 Write(920,*) 'endsolid Fermi surface'
 Close(unit=919)
 Close(unit=920)
 
End Function

Function calc_HC_bzp(bz, ik, cond) Result (HC)
 use adaptive_bzgrid
 Implicit None
 Type (adapt_bzone) :: bz
 Integer :: ik
 Integer :: cond
 Real (Kind=DEF_DBL_PREC) :: HC(3)
 !local
 Integer :: ib

 HC = 0d0
 
 Do ib = 1,bz%bz(ik)%nmod
  HC = HC + calc_HC_bzp_band(bz, ik, cond, ib) 
 Enddo

End Function

Function calc_HC_bzp_band(bz, ik, cond, ib) Result (HC)
 use adaptive_bzgrid
 Implicit None
 Type (adapt_bzone) :: bz
 Integer :: ik, cond, ib
 Real (Kind=DEF_DBL_PREC) :: HC(3)
 !local
 Integer :: inb, inb2
 Logical :: connected

 HC = 0d0

 !print*,'num',8-count(bz%bz(ik)%connect(2:9,ib)>0)

 connected=.true.
 inb = 2
 Do while (connected .and. inb < 10)
  inb2 = inb+1
  If (inb == 9) inb2 = 2
  If (inb .ne. 5 .and. inb .ne. 9) connected = connect_triangle(bz, ik, ib, 1, inb, inb2)
  inb = inb + 1
 Enddo

 
 If (connected) HC = calc_berry(bz, ik, ib, cond, bz%bz(ik)%nmod)
End Function

Function connect_triangle(bz, ik, ib, nb1, nb2, nb3) Result (conn)
 use adaptive_bzgrid
 Implicit None
 Type (adapt_bzone), Target :: bz
 Integer :: ik, ib, nb1, nb2, nb3
 Logical :: conn
 !local
 Type (bz_point), Pointer :: bzp
 Integer :: i12, i13, i23, nb
 Real (kind=DEF_DBL_PREC) :: n(3), r1(3), r2(3), r3(3)

 conn = .false. !by default, not connected
 bzp => bz%bz(ik)

 i12 = bzp%connect(nb2,ib) ! check 1 -> 2
 i13 = bzp%connect(nb3,ib) ! check 1 -> 3
 If (i12 == 0 .or. i13 == 0) Then
  Return
 Endif

 If (nb2 == 2 .or. nb2 == 9) nb = 4
 If (nb2 == 3 .or. nb2 == 4) nb = 6
 If (nb2 == 5 .or. nb2 == 6) nb = 8
 If (nb2 == 7 .or. nb2 == 8) nb = 2

 i23 = bz%bz(bzp%nbs(nb2))%connect(nb,i12) !check 1 -> 2 -> 3

 If (i13 == i23) Then
  conn = .true. !if 1 -> 2 -> 3 is the same as 1 -> 3, we have a triangle
 Else
 Endif

 If (conn) Then
  call find_triangle(bzp, bz%bz(bzp%nbs(nb2)), bz%bz(bzp%nbs(nb3)), ib, i12, i13, r1, r2, r3, n)

  If (bzp%s_z(ib) == 1) Then
   If (bzp%kz(ib) < 0) call write_face(919, n, r1, r2, r3)
   If (bzp%kz(ib) > 0) call write_face(919, n, r1, r3, r2)
  Endif
  If (bzp%s_z(ib) == -1) Then
   If (bzp%kz(ib) < 0) call write_face(920, n, r1, r2, r3)
   If (bzp%kz(ib) > 0) call write_face(920, n, r1, r3, r2)
  Endif
 Endif

End Function

Function calc_berry(bz, ik, ib, cond, nmod) Result (berry_c)
 use adaptive_bzgrid
 Implicit None
 Type (adapt_bzone), Target :: bz
 Integer :: ik, ib, cond, nmod !bzp id, nmod id, cond type (i.e. anom, spin, orb), number of modes at ik
 Real (Kind=DEF_DBL_PREC) :: berry_c(3)
 !local
 Integer :: inb, inb2, ib1, ib2
 Type (bz_point), Pointer :: bzp, bzp1, bzp2
 Complex (Kind=DEF_DBL_PREC) :: prod, Aexp
 Complex (Kind=DEF_DBL_PREC) :: A(nmod,nmod)
 Real (kind=DEF_DBL_PREC) :: n(3), r1(3), r2(3), r3(3)

 bzp => bz%bz(ik)

 berry_c = 0d0
 !if (bzp%kz(ib) > 0) return

 prod = 1d0
 Do inb = 2,9
  inb2 = inb + 1
  if (inb == 9 ) inb2 = 2

  bzp1 => bz%bz(bzp%nbs(inb))
  bzp2 => bz%bz(bzp%nbs(inb2))
  ib1 = bzp%connect(inb, ib)
  ib2 = bzp%connect(inb2, ib)

! draw faces here, such that only faces which are part of Wilson loops are drawn
!!  call find_triangle(bzp, bzp1, bzp2, ib, ib1, ib2, r1, r2, r3, n)
!!  If (bzp%s_z(ib) == 1) Then
!!   If (bzp%kz(ib) < 0) call write_face(919, n, r1, r2, r3)
!!   If (bzp%kz(ib) > 0) call write_face(919, n, r1, r3, r2)
!!  Endif
!!  If (bzp%s_z(ib) == -1) Then
!!   If (bzp%kz(ib) < 0) call write_face(920, n, r1, r2, r3)
!!   If (bzp%kz(ib) > 0) call write_face(920, n, r1, r3, r2)
!!  Endif
!!
! end draw

  prod = prod * dot_product(bzp1%vec(:,ib1),bzp2%vec(:,ib2))
 Enddo

 A = create_operator(cond, nmod)

 Aexp = dot_product(bzp%vec(1:nmod,ib),matmul(A(1:nmod, 1:nmod),bzp%vec(1:nmod,ib)))
 !print*,'Aimag(log(prod))',Aimag(log(prod))
 !print*,'Aexp',Aexp

 berry_c(:) = -0.33333333d0 * (/ bz%bz(ik)%kpar(1), bz%bz(ik)%kpar(2), bz%bz(ik)%kz(ib) /) * Aimag(Log(prod)) * Aexp

End Function

Function create_operator(cond, nmod) Result (A)
 Implicit None
 Integer :: cond, nmod
 Complex (Kind=DEF_DBL_PREC) :: A(nmod, nmod)
 !local
 Integer :: i

 A = 0d0

 If (cond == ANOMALOUS) Then
  Do i = 1,nmod
   A(i,i) = 1d0
  Enddo 
 Elseif (cond == SPIN) Then
  !s_z
  !up
  Do i = 1,nmod/2
   A(i,i) = 1d0
  Enddo 
  !and down
  Do i = nmod/2+1,nmod
   A(i,i) = -1d0
  Enddo 

  !s_x
  !Do i = 1,nmod/2
  ! A(i,i+nmod/2) = 1d0
  !Enddo
  !Do i = nmod/2+1,nmod
  ! A(i,i-nmod/2) = 1d0
  !Enddo
 Elseif (cond == ORBITAL) Then
  !todo
 Endif


End Function

End Module
