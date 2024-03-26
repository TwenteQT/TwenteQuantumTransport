
! This module contains an alternative BZ grid sampling, specifically written for the Fermi surface
! program. It is more involved than the 'bzgrid' module, allowing for adaptive refinement at slight
! memory cost.

#include "math_def.h"
Module adaptive_bzgrid
  Implicit None

  Type t_adapt_bz_opts
    Integer :: nk=10
    Real (Kind=DEF_DBL_PREC) :: tol=1d-3 !refinement tolerance
  End Type t_adapt_bz_opts

  Type bz_point
    Integer :: id ! index of the point on the grid
    Integer :: nbs(9) ! index of its neighbours
    Real (Kind=DEF_DBL_PREC) :: kpar(2) ! kx,ky values of this gridpoint
    Real (Kind=DEF_DBL_PREC), Allocatable :: kz(:) !kz values associated with this kx,ky pair (might be multiple...)
    Complex (Kind=DEF_DBL_PREC), Allocatable :: vec(:,:) !eigenvectors associated with this kx,ky pair (might be multiple
! first index goes over basis functions, second index over nmod
    Integer, Allocatable :: s_z(:) !s_z character of each eigvec
    Integer, Allocatable :: connect(:,:) !connectivity array: first index goes over neighbours, second over nmod
! if the value at (i,j) is k, that means that bz%bz(id)%kz(j) is connected to bz%bz(bz%bz(id)%nbs(i))%kz(k)
! importantly, the first row contains the internal connectivity at this kx,ky pair. A negative number here means
! that the points are degenerate
! Note that this strategy implies that every point at this kx,ky is connect to at most one point at kx2,ky2
    Integer :: nmod ! no. of kz values associated with this kx,ky pair
    Integer :: bp=0 ! 0-interior, 1-bottom, 2-top, 3-left, 4-right, 5-botleft, 6-botright, 7-topleft, 8-topright
  End Type  bz_point

  Type adapt_bzone
    Integer :: nk=0
    Type (bz_point) :: bz(100000)
    Integer, Allocatable :: nc_gp(:) ! list of indices of non-converged gridpoints
    Logical :: conv=.FALSE. ! flag for convergence
    Real (Kind=DEF_DBL_PREC) :: b1(2), b2(2) !translational vectors
  End Type adapt_bzone

  Public :: init_bz, refine_bz, compute_bz_connect, find_triangle, write_face, find_bz_degenerate

  Private :: setup_bz, generate_grid

  Contains

  Subroutine init_bz(bz, vbr, bzopt)
    Use logging
    Implicit None
    Type (t_adapt_bz_opts), Intent (In) :: bzopt
    Type (adapt_bzone) :: bz
    Real (Kind=DEF_DBL_PREC), Intent (In) :: vbr(2, 2)
    ! local
    Real (Kind=DEF_DBL_PREC) :: vbg(2,2)

    Call setup_bz(vbr,vbg)
    Call generate_grid(bz, bzopt%nk, vbg)
    bz%b1(:) = vbg(:,1)
    bz%b2(:) = vbg(:,2)

  End Subroutine init_bz
  
  Subroutine setup_bz(vbr,vbg) 
    Use logging
    Implicit None
    Real (Kind=DEF_DBL_PREC), Intent (In) :: vbr(2, 2)
    ! local
    Real (Kind=DEF_DBL_PREC) :: det, cnorm, vbg(2,2)
    Real (Kind=DEF_DBL_PREC) :: romega, gomega
    Character (Len=30) :: cwork

    ! set some parameters
    det = vbr(1,1)*vbr(2,2) - vbr(1,2)*vbr(2,1)
    cnorm = 2.0d0 * DEF_M_PI / det
    romega = Abs(det)
    gomega = (2.0d0*DEF_M_PI)**2 / romega

    ! compute reciprocial vectors
    vbg(1,1) = cnorm * vbr(2,2)
    vbg(2,1) = - cnorm * vbr(1,2)
    vbg(1,2) = - cnorm * vbr(2,1)
    vbg(2,2) = cnorm * vbr(1,1)

    ! some logging
      Write (cwork, '(g15.7)') romega
      Call do_log (1, '   AREA OF 2D-PRIMITIVE CELL: REAL= '//trim(cwork))
      Write (cwork, '(g15.7)') gomega
      Call do_log (1, '                        RECIPROCAL= '//trim(cwork))
      Write (cwork, '(2 g15.7)') vbg (:, 1)
      Call do_log (1, '   RECIP. BASIS: 1. VECTOR = '//trim(cwork))
      Write (cwork, '(2 g15.7)') vbg (:, 2)
      Call do_log (1, '                 2. VECTOR = '//trim(cwork))
  End Subroutine setup_bz

  Subroutine generate_grid(bz, nk, vbg)
    Use logging
    Implicit None
    Type (adapt_bzone), Target :: bz
    Integer :: nk
    Real (Kind=DEF_DBL_PREC), Intent (In) :: vbg(2, 2)
    ! local
    Integer :: ibz, i1, i2, nk_pr
    Real (Kind=DEF_DBL_PREC) :: part
    Real (Kind=DEF_DBL_PREC) :: fx1, fx2
    Character (Len=30) :: cwork

    part = 2d0

    nk_pr = 2*nk-2 !number of kpoints per row (or column)
    bz%nk = nk_pr**2 !set number of kpoints in grid derived type
    Allocate (bz%nc_gp(bz%nk)) ! allocate array of non-converged gridpoints (nc_gp)
    ibz = 0

    Do i1 = -nk+1, nk-2 !the nk-1'th element would describe the exact same point as -nk+1 due to period bc
      fx1 = 1.0d0* i1 / (part*(nk-1))
      Do i2 = -nk+1, nk-2
        fx2 = 1.0d0* i2 / (part*(nk-1))
        ibz = ibz + 1
        bz%nc_gp(ibz) = ibz ! at this point, every gridpoint is nonconverged
        bz%bz(ibz)%id = ibz ! set id values

        bz%bz(ibz)%kpar(:) = fx1*vbg(:,1) + fx2*vbg(:,2)
        
        bz%bz(ibz)%nbs(1) = ibz ! first element points at itself
        ! set neighbours
        bz%bz(ibz)%nbs(2) = ibz+1
        bz%bz(ibz)%nbs(3) = ibz+nk_pr+1
        bz%bz(ibz)%nbs(4) = ibz+nk_pr
        !If (bz%bz(ibz)%bp == 0 ) Then
         bz%bz(ibz)%nbs(5) = ibz+nk_pr-1
         bz%bz(ibz)%nbs(6) = ibz-1
         bz%bz(ibz)%nbs(7) = ibz-nk_pr-1
         bz%bz(ibz)%nbs(8) = ibz-nk_pr
         bz%bz(ibz)%nbs(9) = ibz-nk_pr+1
        !Endif
        If (i2 == -nk+1) Then !bottom edge
         bz%bz(ibz)%bp = 1
         bz%bz(ibz)%nbs(5) = ibz + nk_pr - 1 + nk_pr
         bz%bz(ibz)%nbs(6) = ibz + nk_pr - 1
         bz%bz(ibz)%nbs(7) = ibz + nk_pr - 1 - nk_pr
        Endif
        If (i2 == nk-2) Then !top edge
         bz%bz(ibz)%bp = 2
         bz%bz(ibz)%nbs(2) = ibz - nk_pr + 1
         bz%bz(ibz)%nbs(3) = ibz - nk_pr + 1 + nk_pr
         bz%bz(ibz)%nbs(9) = ibz - nk_pr + 1 - nk_pr
        Endif
        If (i1 == -nk+1) Then !left edge
         bz%bz(ibz)%bp = 3
         bz%bz(ibz)%nbs(7) = ibz + nk_pr*(nk_pr-1) - 1
         bz%bz(ibz)%nbs(8) = ibz + nk_pr*(nk_pr-1) 
         bz%bz(ibz)%nbs(9) = ibz + nk_pr*(nk_pr-1) + 1
        Endif
        If (i1 == nk-2) Then !right edge
         bz%bz(ibz)%bp = 4
         bz%bz(ibz)%nbs(3) = ibz - nk_pr*(nk_pr-1) + 1
         bz%bz(ibz)%nbs(4) = ibz - nk_pr*(nk_pr-1)
         bz%bz(ibz)%nbs(5) = ibz - nk_pr*(nk_pr-1) - 1
        Endif
        If (i2 == -nk+1 .and. i1 == -nk+1) Then !bottom left
         bz%bz(ibz)%bp = 5
         bz%bz(ibz)%nbs(7) = bz%nk
        Endif
        If (i2 == -nk+1 .and. i1 ==  nk-2) Then !bottom right
         bz%bz(ibz)%bp = 6
         bz%bz(ibz)%nbs(5) = ibz - (nk_pr-1)*(nk_pr-1)
        Endif
        If (i2 ==  nk-2 .and. i1 == -nk+1) Then !top left
         bz%bz(ibz)%bp = 7
         bz%bz(ibz)%nbs(9) = ibz + (nk_pr-1)*(nk_pr-1)
        Endif
        If (i2 ==  nk-2 .and. i1 ==  nk-2) Then !top right
         bz%bz(ibz)%bp = 8
         bz%bz(ibz)%nbs(3) = 1
        Endif

      Enddo
    Enddo

    Write (cwork, '(i8)') bz%nk 
    Call do_log (1, '   NUMBER OF K||-POINTS:   NBZ='//trim(cwork))
    Call do_log (1, '')

  End Subroutine
!  3 6 9
!  2 5 8   ids
!  1 4 7
    ! #############
    ! #           #
    ! #  9  2  3  #
    ! #  8  1  4  # nbs
    ! #  7  6  5  #
    ! #           #
    ! #############

  Subroutine refine_bz(bz, tol)
    Use logging
    Implicit None
    Type (adapt_bzone), Target :: bz
    Real (Kind=DEF_DBL_PREC) :: tol
    ! local
    Integer :: igp, ik, inb, nb
    Real (Kind=DEF_DBL_PREC) :: avk
    Logical :: conv
    Integer :: nc(size(bz%nc_gp)*5), num_nc 
    Integer :: nc_tmp(size(bz%nc_gp)*5), num_nc_tmp
    Real (Kind=DEF_DBL_PREC) :: kpar(4,2)
    Character (Len=30) :: cwork

    num_nc = 0 ! we are optimistic people: assume that all gridpoints will be converged :)

    ! loop over non-converged gridpoints
    ! 'non-converged gridpoint' really means that the parallelogram of which ik is the bottom-left corner
    ! is not converged.. really every gridpoint has one unique parallelogram associated with it. This happens
    ! because of the periodic boundary conditions (very nice)
    Do igp = 1, size(bz%nc_gp)
      ik = bz%nc_gp(igp) 

      if (bz%bz(ik)%bp > 0) cycle

      conv = .true.
      Do inb = 2,4
        nb = bz%bz(ik)%nbs(inb)
        if (Abs(bz%bz(nb)%kpar(1) - bz%bz(ik)%kpar(1)) < 2*DEF_M_PI/100) cycle
        If (Abs(bz%bz(nb)%kz(1) - bz%bz(ik)%kz(1)) > tol) then
         conv = .FALSE.
        endif
      Enddo
      
      If (Abs(Abs(bz%bz(ik)%kz(1) - bz%bz(ik)%kz(2)) - 2*DEF_M_PI) < tol) conv = .TRUE.
      Do inb = 2,4
        nb = bz%bz(ik)%nbs(inb)
        If (Abs(bz%bz(nb)%kz(1) - bz%bz(ik)%kz(1)) < 0.01*tol) conv = .TRUE.
      Enddo


      If (.not. conv) Then
        num_nc = num_nc + 1
        nc(num_nc) = ik !save the index of the non-conv grid point in this array
      Endif
    Enddo

    If (num_nc == 0) Then
      Call do_log(1, 'Converged k-grid!')      
      Deallocate(bz%nc_gp)
      bz%conv = .TRUE.
      return
    Endif
  
    ! sketch of the gridpoints
    !   X are the inb=1,4 of the old grid
    !   O are the to-be-added gridpoints
    ! ######################
    ! #                    #
    ! #  X (2) 0     X (3) #
    ! #                    #
    ! #  0     0     0     #
    ! #                    #
    ! #  X (1) 0     X (4) #
    ! #                    #
    ! ######################
    ! the new neighbours of (1)
    ! ######################
    ! #                    #
    ! #  X     0     X     #
    ! #                    #
    ! #  0 (2) 0 (3) 0     #
    ! #                    #
    ! #  X (1) 0 (4) X     #
    ! #                    #
    ! ######################
    ! The new gridpoints are given indices
    ! ############################
    ! #                          #
    ! #  X       0(nk+4) X       #
    ! #                          #
    ! #  0(nk+1) 0(nk+3) 0(nk+5) #
    ! #                          #
    ! #  X       0(nk+2) X       #
    ! #                          #
    ! ############################
    ! Only nk+1, nk+3 and nk+2 have an associated parallelogram, nk+4 and nk+5 don't.

    num_nc_tmp = 0 !to keep track of the number of gridpoints we have generated
    ! now loop over non-converged gridpoints and generate new gridpoints
    Do igp = 1, num_nc
      ik = nc(igp)
      ! set kpar values
      ! define some shorthands
      kpar(1,:) = bz%bz(ik)%kpar(:)
      kpar(2,:) = bz%bz(bz%bz(ik)%nbs(2))%kpar(:)
      kpar(3,:) = bz%bz(bz%bz(ik)%nbs(3))%kpar(:)
      kpar(4,:) = bz%bz(bz%bz(ik)%nbs(4))%kpar(:)

      bz%bz(bz%nk+1)%id = bz%nk+1
      bz%bz(bz%nk+2)%id = bz%nk+2
      bz%bz(bz%nk+3)%id = bz%nk+3
      bz%bz(bz%nk+4)%id = bz%nk+4
      bz%bz(bz%nk+5)%id = bz%nk+5

      ! 1,2,3 are halfway between the current gridpoint and its 4 neighbours
      bz%bz(bz%nk+1)%kpar(:) = kpar(1,:) + 0.5d0 * (kpar(2,:) - kpar(1,:))
      bz%bz(bz%nk+3)%kpar(:) = kpar(1,:) + 0.5d0 * (kpar(3,:) - kpar(1,:))
      bz%bz(bz%nk+2)%kpar(:) = kpar(1,:) + 0.5d0 * (kpar(4,:) - kpar(1,:))
      ! 4,5 are in between the neighbours
      bz%bz(bz%nk+4)%kpar(:) = kpar(2,:) + 0.5d0 * (kpar(3,:) - kpar(2,:))
      bz%bz(bz%nk+5)%kpar(:) = kpar(4,:) + 0.5d0 * (kpar(3,:) - kpar(4,:))

      ! set neighbours, only necessary for 1,2,3
      bz%bz(bz%nk+1)%nbs(1) = bz%nk+1
      bz%bz(bz%nk+1)%nbs(2) = bz%bz(ik)%nbs(2)
      bz%bz(bz%nk+1)%nbs(3) = bz%nk+4
      bz%bz(bz%nk+1)%nbs(4) = bz%nk+3

      bz%bz(bz%nk+3)%nbs(1) = bz%nk+3
      bz%bz(bz%nk+3)%nbs(2) = bz%nk+4
      bz%bz(bz%nk+3)%nbs(3) = bz%bz(ik)%nbs(3)
      bz%bz(bz%nk+3)%nbs(4) = bz%nk+5

      bz%bz(bz%nk+2)%nbs(1) = bz%nk+2
      bz%bz(bz%nk+2)%nbs(2) = bz%nk+3
      bz%bz(bz%nk+2)%nbs(3) = bz%nk+5
      bz%bz(bz%nk+2)%nbs(4) = bz%bz(ik)%nbs(4)

      ! finally update the neighbours of the current gridpoint
      bz%bz(ik)%nbs(2) = bz%nk+1
      bz%bz(ik)%nbs(3) = bz%nk+3
      bz%bz(ik)%nbs(4) = bz%nk+2

      ! Finally, add gridpoints to list of non-converged gridpoints
      ! We can't just add them to 'nc' as we are currently looping over those...
      nc_tmp(num_nc_tmp+1) = bz%nk+1
      nc_tmp(num_nc_tmp+2) = bz%nk+2
      nc_tmp(num_nc_tmp+3) = bz%nk+3
      nc_tmp(num_nc_tmp+4) = bz%nk+4
      nc_tmp(num_nc_tmp+5) = bz%nk+5
      num_nc_tmp = num_nc_tmp+5
      ! We don't want to test gridpoints 4 and 5 for convergence.. so flag them as if they were boundary points
      bz%bz(bz%nk+4)%bp = 99
      bz%bz(bz%nk+5)%bp = 99

      bz%nk = bz%nk + 5 ! we have added 5 gridpoints
    Enddo

    !Finally add the newly generated gridpoints to the list of non-converged gridpoints
    Do ik = 1, num_nc_tmp 
      num_nc = num_nc + 1
      nc(num_nc) = nc_tmp(ik)
    Enddo

    Deallocate(bz%nc_gp)

    Allocate(bz%nc_gp(num_nc))

    Do ik = 1, num_nc
      bz%nc_gp(ik) = nc(ik)
    Enddo

    Write (cwork, '(i8)') num_nc
    Call do_log (1, '   NUMBER OF K||-POINTS:   NBZ='//trim(cwork))
    Call do_log (1, '')

    End Subroutine refine_bz
 
    Subroutine diagonalize_degenerate_subspace(bz)
     Implicit None
     Type (adapt_bzone), Target :: bz
     !local
     Type (bz_point), Pointer :: bzp
     Integer :: ik, ib, ib2, flags(100), n !assuming nmod <= 100
     Complex (Kind=DEF_DBL_PREC) :: work_vec(100,2) !assuming 2*nl <= 100

     Do ik = 1,bz%nk
      bzp => bz%bz(ik)
      flags = 0
      Do ib = 1,bzp%nmod
       !check for degenerate states, it'll be ib2 if it exists
       If (flags(ib) > 0) cycle

       If (bzp%connect(1,ib) >= 0) Then
        cycle
       Else
        ib2 = -bzp%connect(1,ib)
        !eliminate ib2 from the loop
        flags(ib2) = 1
       Endif 

       work_vec(:,:) = 0d0
       n = size(bzp%vec,1)
       !diagonalize wrt s_z
       work_vec(1:n,1) = bzp%vec(1:n,ib)
       work_vec(1:n,2) = bzp%vec(1:n,ib2)
       bzp%vec(1:n,1:2) = diag_sz(work_vec(1:n,1:2), n)

       !write to bzp%vec(1:n,ib) and bzp%vec(1:n,ib2)

       bzp%s_z(ib) = 1
       bzp%s_z(ib2) = -1

      Enddo !ib

     Enddo !ik

    Contains 

     Function diag_sz(vec_i, n) Result (vec_o)
      Implicit None
      Integer :: N
      Complex (Kind=DEF_DBL_PREC) :: vec_i(n,2), vec_o(n,2)
      !local
      Real (Kind=DEF_DBL_PREC) :: s_z(n,n)
      Integer :: i
      Complex (Kind=DEF_DBL_PREC) :: u(n), v(n), B(2,2), alpha(2), beta(2)
      Real (Kind=DEF_DBL_PREC) :: det

      u(1:n) = vec_i(1:n,1)
      v(1:n) = vec_i(1:n,2)

      s_z = 0d0

      Do i = 1,n/2
       s_z(i,i) = 1d0
      Enddo
      Do i = n/2+1,n
       s_z(i,i) = -1d0
      Enddo

      B(1,1) = Dot_product(u(1:n),Matmul(s_z(1:n, 1:n),u(1:n)))
      B(1,2) = Dot_product(u(1:n),Matmul(s_z(1:n, 1:n),v(1:n)))
      B(2,1) = Dot_product(v(1:n),Matmul(s_z(1:n, 1:n),u(1:n)))
      B(2,2) = Dot_product(v(1:n),Matmul(s_z(1:n, 1:n),v(1:n)))

      det = sqrt(4*Abs(B(2,1))**2 + (B(1,1)-B(2,2))**2)
      !print*,'======================================================='
      !print*,'B',B
   
      alpha(1) = (B(1,1)-B(2,2)+det)/(2*B(2,1))
      alpha(2) = 1d0
      beta(1) = (B(1,1)-B(2,2)-det)/(2*B(2,1))
      beta(2) = 1d0

      alpha = alpha / sqrt(alpha(1)**2+alpha(2)**2)
      beta = beta / sqrt(beta(1)**2+beta(2)**2)

      !Print*,'eigval1',dot_product(alpha(1:2), matmul(B, alpha(1:2)))
      !Print*,'eigval2',dot_product(beta(1:2), matmul(B, beta(1:2)))


      vec_o(1:n,1) = alpha(1)*u(1:n) + alpha(2)*v(1:n)
      vec_o(1:n,2) = beta(1)*u(1:n) + beta(2)*v(1:n)
 
      !Print*,'expval1',dot_product(vec_o(1:n,1), matmul(s_z, vec_o(1:n,1)))
      !Print*,'expval2',dot_product(vec_o(1:n,2), matmul(s_z, vec_o(1:n,2)))
      
     End Function

    End Subroutine diagonalize_degenerate_subspace

    Subroutine find_bz_degenerate(bz)
     Implicit None
     Type (adapt_bzone), Target :: bz
     !local
     Integer :: ik, ib1, ib2
     Type (bz_point), Pointer :: bzp

     Do ik = 1,bz%nk ! for every kpoint
      bzp => bz%bz(ik)

      !check internal connectivity
      Do ib1 = 1,bzp%nmod
       Do ib2 = 1,bzp%nmod
        If (Abs(bzp%kz(ib1) - bzp%kz(ib2)) < 1d-10) Then !degenerate points
         bzp%connect(1,ib2) = -ib1
        Endif
       Enddo
      Enddo
     Enddo

    End Subroutine find_bz_degenerate

    Subroutine compute_bz_connect(bz, tol)
     Implicit None
     Type (adapt_bzone), Target :: bz
     Real (Kind=DEF_DBL_PREC) :: tol
     !local
     Integer :: ik, inb, ib1, ib2, closest_ib2
     Type (bz_point), Pointer :: bzp, bzp_nb
     Real (Kind=DEF_DBL_PREC) :: dist

     Do ik = 1,bz%nk ! for every kpoint
      bzp => bz%bz(ik)

      !check internal connectivity
      Do ib1 = 1,bzp%nmod
       Do ib2 = 1,bzp%nmod
        If (Abs(bzp%kz(ib1) - bzp%kz(ib2)) < tol) Then !connected points
         bzp%connect(1,ib2) = ib1
        Endif
       Enddo
      Enddo

      Do inb = 6,8,1 ! check all neighbours to the bottom, bottom left and left
       bzp_nb => bz%bz(bzp%nbs(inb))

       Do ib1 = 1,bzp_nb%nmod ! loop over meshpoints to connect to
        dist = 100d0
        closest_ib2 = 0
        Do ib2 = 1,bzp%nmod   !                      to connect
         ! find the ib2 closest to ib1
         ! and only allow same-spin bands to connect
         If (Abs(bzp_nb%kz(ib1)-bzp%kz(ib2)) < dist .and. bzp_nb%s_z(ib1) == bzp%s_z(ib2)) Then
          dist = Abs(bzp_nb%kz(ib1)-bzp%kz(ib2))
          closest_ib2 = ib2
         Endif
        Enddo
        If (dist < tol) Then
         bzp%connect(inb,closest_ib2) = ib1
         bzp_nb%connect(inb-4,ib1) = closest_ib2
        Endif
       Enddo

      Enddo
     Enddo
 
    End Subroutine compute_bz_connect

    Subroutine print_bz_as_mesh(bz)
     Implicit None
     Type (adapt_bzone), Target :: bz
     !local
     Type (bz_point), Pointer :: bzp, bzp_top, bzp_topright, bzp_right
     Integer :: ik, ib, i_t, i_tr, i_r, i_t_r, i_r_t, ib_n, i_t_n
     Real (Kind=DEF_DBL_PREC) :: n(3), r1(3), r2(3), r3(3), a1(3), a2(3)
     Integer :: f1, f2, f3

     Open(file='fermiS.stl',unit=919)

     Write(919,*) 'solid Fermi surface'

     Do ik = 1,bz%nk
      bzp => bz%bz(ik)
      if (bzp%bp > 0) cycle

      bzp_top => bz%bz(bzp%nbs(2))
      bzp_topright => bz%bz(bzp%nbs(3))
      bzp_right => bz%bz(bzp%nbs(4))

      Do ib = 1,bzp%nmod
       i_t  = bzp%connect(2,ib)
       i_tr = bzp%connect(3,ib)
       i_r  = bzp%connect(4,ib)

       f1 = 0
       f2 = 0
       f3 = 0

       If (i_t /= 0) i_t_r = bzp_top%connect(4,i_t)
       If (i_r /= 0) i_r_t = bzp_right%connect(2,i_r)

       !triangle within first part of quad defined by kx,ky and its 2,3,4, nbs
       If (i_tr /= 0 .and. i_t /= 0 .and. i_t_r /= 0) Then
       If (i_tr == i_t_r .or. i_tr == -bzp_topright%connect(1,i_t_r) .or. i_t_r == -bzp_topright%connect(1,i_tr)) Then
        call find_triangle(bzp,bzp_top,bzp_topright, ib,i_t,i_tr, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
        f1 = 1
       Endif
       Endif

       !triangle within second part of quad defined by kx,ky and its 2,3,4, nbs
       If (i_tr /= 0 .and. i_r /= 0 .and. i_r_t /= 0) Then
       If (i_tr == i_r_t .or. i_tr == -bzp_topright%connect(1,i_r_t) .or. i_r_t == -bzp_topright%connect(1,i_tr)) Then
        call find_triangle(bzp,bzp_topright,bzp_right, ib,i_tr,i_r, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
        f2 = 1
       Endif
       Endif

       If (i_t /= 0) Then
       If (f1+f2==2 .and. i_r == bzp_top%connect(5,i_t) .and. i_r /= 0) Then
        call find_triangle(bzp,bzp_top,bzp_right, ib,i_t,i_r, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
        f3 = 1
       Endif
       Endif

       If (f1+f2+f3 > 1) cycle

       !connectivity between the kz(:) in connect(:) is defined such that every kz(:) is connected to at most
       !one point at kx+1,ky+1, so there can exist no triangles parralel to kz. quads can exist, though
       !quads between id and nbs(2)
       If (i_t /= 0) i_t_n = bzp_top%connect(1,i_t)
       ib_n = bzp%connect(1,ib)
       If (ib_n == 0) cycle

       If (i_t_n > 0) Then
       If (ib_n == i_t_n .and. i_t /= 0) Then
        call find_triangle(bzp,bzp_top,bzp_top, ib,i_t_n,i_t, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
        call find_triangle(bzp,bzp,bzp_top, ib,ib_n,i_t_n, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
       Endif
       Endif
       !quads between id and nbs(4)
       If (i_r /= 0) i_t_n = bzp_right%connect(1,i_r)
       If (i_t_n > 0) then
       If (ib_n == i_t_n .and. i_r /= 0) Then
        call find_triangle(bzp,bzp_right,bzp_right, ib,i_t_n,i_r, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
        call find_triangle(bzp,bzp,bzp_right, ib,ib_n,i_t_n, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
       Endif
       Endif
       !quads between id and nbs(3)

       If (i_tr /= 0) i_t_n = bzp_topright%connect(1,i_tr)
       If (i_t_n > 0) Then
       If (ib_n == i_t_n .and. i_tr /= 0) Then
        call find_triangle(bzp,bzp_topright,bzp_topright, ib, i_t_n, i_tr, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
        call find_triangle(bzp,bzp,bzp_topright, ib, ib_n, i_t_n, r1,r2,r3, n)
        call write_face(919, n, r1, r2, r3)
       Endif
       Endif



      Enddo
     Enddo 
     Write(919,*) 'endsolid Fermi surface'

     Close(unit=919)
  End subroutine

     Subroutine find_triangle(bzp1, bzp2, bzp3, ib1, ib2, ib3, r1, r2, r3, n)
      Implicit None
      Type (bz_point) :: bzp1, bzp2, bzp3
      Integer :: ib1, ib2, ib3
      Real (Kind=DEF_DBL_PREC) :: n(3), r1(3), r2(3), r3(3)
      !local
      Real (Kind=DEF_DBL_PREC) :: a1(3), a2(3)

      r1 = (/ bzp1%kpar(1), bzp1%kpar(2), bzp1%kz(ib1) /)
      r2 = (/ bzp2%kpar(1), bzp2%kpar(2), bzp2%kz(ib2) /)
      r3 = (/ bzp3%kpar(1), bzp3%kpar(2), bzp3%kz(ib3) /)

      a1 = r2 - r1
      a2 = r3 - r1

      n = sign(1d0,bzp1%kz(ib1))*(/ a1(2)*a2(3) - a1(3)*a1(2), a1(3)*a2(1) - a1(1)*a2(3), a1(1)*a2(2) - a1(2)*a2(1) /)
      
     End Subroutine

     Subroutine write_face(fileunit, n, r1, r2, r3)
      Implicit None
      Integer :: fileunit
      Real (Kind=DEF_DBL_PREC) :: n(3), r1(3), r2(3), r3(3)

        Write(fileunit,"(A12, 1X, ES14.7, 1X, ES14.7, 1X, ES14.7)") 'facet normal', n(1), n(2), n(3)
        Write(fileunit,*) 'outer loop'
        Write(fileunit,"(A6, 1X, ES14.7, 1X, ES14.7, 1X, ES14.7)") 'vertex', r1(1), r1(2), r1(3)
        Write(fileunit,"(A6, 1X, ES14.7, 1X, ES14.7, 1X, ES14.7)") 'vertex', r2(1), r2(2), r2(3)
        Write(fileunit,"(A6, 1X, ES14.7, 1X, ES14.7, 1X, ES14.7)") 'vertex', r3(1), r3(2), r3(3)
        Write(fileunit,*) 'endloop'
        Write(fileunit,*) 'endfacet'

     End Subroutine


! 9 2 3
! 8 1 4
! 7 6 5


End
