
#include "math_def.h"
Module OAMdensity
  Implicit None

  Type OAM_density
    Integer :: alloc = 0
    Integer :: nl = 0, nr = 0, nm = 0
    Real (Kind=DEF_DBL_PREC), Allocatable :: Lir (:, :), Lt (:, :), Rir (:, :), Rt (:, :)
    Real (Kind=DEF_DBL_PREC), Allocatable :: Lm (:, :), Rm (:, :)
  End Type OAM_density

Contains
  
  Subroutine alloc_oam_density (oam_dens, mna, lna, rna)
    Implicit None
    Type (OAM_density), Intent (Inout) :: oam_dens
    Integer, Intent(In) :: mna, lna, rna

    If (oam_dens%alloc /= 0 .And. (oam_dens%nm /= mna .Or. oam_dens%nl /= lna .Or. oam_dens%nr /= rna)) Then
        call free_oam_density(oam_dens)
    End If

    If (oam_dens%alloc == 0) Then
       Allocate (oam_dens%Lm(3, mna))
       Allocate (oam_dens%Rm(3, mna))

       Allocate (oam_dens%Lir(3, lna))
       Allocate (oam_dens%Lt(3, rna))

       Allocate (oam_dens%Rir(3, rna))
       Allocate (oam_dens%Rt(3, lna))
       oam_dens%alloc = 1
       oam_dens%nl = lna
       oam_dens%nr = rna
       oam_dens%nm = mna
    End If

    Call init_oam_density (oam_dens)

  End Subroutine alloc_oam_density

  Subroutine init_oam_density (oam_dens)
    Implicit None
    Type (OAM_density), Intent (Inout) :: oam_dens
    oam_dens%Lm (:, :) = 0.0d0
    oam_dens%Rm (:, :) = 0.0d0
    oam_dens%Rir (:, :) = 0.0d0
    oam_dens%Rt (:, :) = 0.0d0
    oam_dens%Lir (:, :) = 0.0d0
    oam_dens%Lt (:, :) = 0.0d0
  End Subroutine init_oam_density

  Subroutine free_oam_density (oam_dens)
    Implicit None
    Type (OAM_density), Intent (Inout) :: oam_dens
    if (oam_dens%alloc==0) return
    Deallocate (oam_dens%Rm, oam_dens%Lm, oam_dens%Lir, oam_dens%Lt, oam_dens%Rir, oam_dens%Rt)
    oam_dens%nm = 0
    oam_dens%nl = 0
    oam_dens%nr = 0
    oam_dens%alloc = 0
  End Subroutine free_oam_density

  Subroutine bzsum_oam_density (oam_dens, oam_dens1, bzw)
    Implicit None
    Type (OAM_density), Intent (Inout) :: oam_dens, oam_dens1
    Real (Kind=DEF_DBL_PREC), Intent(In) :: bzw

    oam_dens%Lm (:, :) = oam_dens%Lm(:, :) + oam_dens1%Lm(:, :) * bzw
    oam_dens%Rm (:, :) = oam_dens%Rm(:, :) + oam_dens1%Rm(:, :) * bzw
 
    oam_dens%Rir (:, :) = oam_dens%Rir(:, :) + oam_dens1%Rir(:, :) * bzw
    oam_dens%Rt (:, :) = oam_dens%Rt(:, :) + oam_dens1%Rt(:, :) * bzw
 
    oam_dens%Lir (:, :) = oam_dens%Lir(:, :) + oam_dens1%Lir(:, :) * bzw
    oam_dens%Lt (:, :) = oam_dens%Lt(:, :) + oam_dens1%Lt(:, :) * bzw
  End Subroutine bzsum_oam_density

  Subroutine calc_oam_density_EMTO (rhs, trgeo, oam_dens, nl, nr)
    use omta_SOC !needed for make_bigLmat and the t_Lop_mats type
    Implicit None
    Type (t_geometry_EMTO), Intent (In) :: trgeo
    Type (zdensemat), target, Intent (In) :: rhs
    Integer, Intent (In) :: nl, nr
    Type (OAM_density), Intent (Inout) :: oam_dens
!!$ Local vars
    Integer :: imd, inc
    Complex (Kind=DEF_DBL_PREC), Pointer :: phi1 (:)
    Integer :: ofs1o, ofs2o, ofs1a, ofs2a
    Type(t_Lop_mats) :: Lop(trgeo%nlmax)
!!$ **************************************************
!!$ The wave functions are already flux-normalized
!!$ but there is a 2*pi factor to be considered later
!!$ **************************************************
    call make_bigLmat(trgeo%nlmax, Lop)

    Call init_oam_density (oam_dens)
    inc = 1
    ofs1o = trgeo%ltno * 2 - trgeo%num_emb_l
    ofs2o = (trgeo%tno+trgeo%ltno) * 2 - trgeo%num_emb_l - trgeo%num_emb_m
    ofs1a = trgeo%ltna
    ofs2a = trgeo%ltna + trgeo%tna
    Do imd = 1, nl + nr
       phi1 => rhs%bl (:, imd)
       If (imd > nl) Then
          Call oam_helper (phi1, oam_dens%Rt(:, :), trgeo, 0, 0, trgeo%ltna, Lop)
          Call oam_helper (phi1, oam_dens%Rm(:, :), trgeo, ofs1o, ofs1a, trgeo%tna, Lop)
          Call oam_helper (phi1, oam_dens%Rir(:, :), trgeo, ofs2o, ofs2a, trgeo%rtna, Lop)
       Else
          Call oam_helper (phi1, oam_dens%Lir(:, :), trgeo, 0, 0, trgeo%ltna, Lop)
          Call oam_helper (phi1, oam_dens%Lm(:, :), trgeo, ofs1o, ofs1a, trgeo%tna, Lop)
          Call oam_helper (phi1, oam_dens%Lt(:, :), trgeo, ofs2o, ofs2a, trgeo%rtna, Lop)
       End If

    End Do

  Contains
    Subroutine oam_helper (phi, oam_dens, trgeo, ofso, ofsa, tna, Lop)
      Implicit none
      Real (Kind=DEF_DBL_PREC), Intent (Inout) :: oam_dens (:, :)
      Type (t_geometry_EMTO), Intent (In) :: trgeo
      Complex (Kind=DEF_DBL_PREC) :: phi (:)
      Integer, Intent (In) :: ofsa, ofso, tna
      Type(t_Lop_mats) :: Lop(:)
    !!$ Local vars
      Integer :: ia, ofs, no
      Real (Kind=DEF_DBL_PREC) :: ldens(3)
      Type (tclass), Pointer :: at
      ! this routine computes <psi|L_a|psi> for each atom
      ! where a = {x,y,z}, and adds the values found for up and down spins
      ofs = ofso + 1
      Do ia = 1, tna
        at => trgeo%atoms(ia+ofsa)%ptr
        no = (at%lmx+1) **2
        ldens = 0d0
        ! first spin up
        ldens(1) = dot_product(phi(ofs:ofs+no), matmul(Lop(at%lmx+1)%Lx(1:no, 1:no), phi(ofs:ofs+no))) 
        ldens(2) = dot_product(phi(ofs:ofs+no), matmul(Lop(at%lmx+1)%Ly(1:no, 1:no), phi(ofs:ofs+no)))  
        ldens(3) = dot_product(phi(ofs:ofs+no), matmul(Lop(at%lmx+1)%Lz(1:no, 1:no), phi(ofs:ofs+no)))  
        ! then spin down
        ldens(1) = ldens(1) + dot_product(phi(ofs+no+1:ofs+2*no), matmul(Lop(at%lmx+1)%Lx(1:no, 1:no), phi(ofs+no+1:ofs+2*no))) 
        ldens(2) = ldens(2) + dot_product(phi(ofs+no+1:ofs+2*no), matmul(Lop(at%lmx+1)%Ly(1:no, 1:no), phi(ofs+no+1:ofs+2*no))) 
        ldens(3) = ldens(3) + dot_product(phi(ofs+no+1:ofs+2*no), matmul(Lop(at%lmx+1)%Lz(1:no, 1:no), phi(ofs+no+1:ofs+2*no)))  

        oam_dens(:, ia) = oam_dens(:, ia) + ldens(:)

        ofs = ofs + 2*no
      End Do

    End Subroutine oam_helper

  End Subroutine calc_oam_density_EMTO

  Subroutine calc_oam_density (rhs, trgeo, oam_dens, nl, nr)
    use omta_SOC !needed for make_bigLmat and the t_Lop_mats type
    Implicit None
    Type (t_geometry), Intent (In) :: trgeo
    Type (zdensemat), target, Intent (In) :: rhs
    Integer, Intent (In) :: nl, nr
    Type (OAM_density), Intent (Inout) :: oam_dens
!!$ Local vars
    Integer :: imd, inc
    Complex (Kind=DEF_DBL_PREC), Pointer :: phi1 (:)
    Integer :: ofs1o, ofs2o, ofs1a, ofs2a
    Type(t_Lop_mats) :: Lop(trgeo%nlmax)
!!$ **************************************************
!!$ The wave functions are already flux-normalized
!!$ but there is a 2*pi factor to be considered later
!!$ **************************************************
    call make_bigLmat(trgeo%nlmax, Lop)

    Call init_oam_density (oam_dens)
    inc = 1
    ofs1o = trgeo%ltno * 2
    ofs2o = (trgeo%tno+trgeo%ltno) * 2
    ofs1a = trgeo%ltna
    ofs2a = trgeo%ltna + trgeo%tna
    Do imd = 1, nl + nr
       phi1 => rhs%bl (:, imd)
       If (imd > nl) Then
          Call oam_helper (phi1, oam_dens%Rt(:, :), trgeo, 0, 0, trgeo%ltna, Lop)
          Call oam_helper (phi1, oam_dens%Rm(:, :), trgeo, ofs1o, ofs1a, trgeo%tna, Lop)
          Call oam_helper (phi1, oam_dens%Rir(:, :), trgeo, ofs2o, ofs2a, trgeo%rtna, Lop)
       Else
          Call oam_helper (phi1, oam_dens%Lir(:, :), trgeo, 0, 0, trgeo%ltna, Lop)
          Call oam_helper (phi1, oam_dens%Lm(:, :), trgeo, ofs1o, ofs1a, trgeo%tna, Lop)
          Call oam_helper (phi1, oam_dens%Lt(:, :), trgeo, ofs2o, ofs2a, trgeo%rtna, Lop)
       End If

    End Do

  Contains
    Subroutine oam_helper (phi, oam_dens, trgeo, ofso, ofsa, tna, Lop)
      Implicit none
      Real (Kind=DEF_DBL_PREC), Intent (Inout) :: oam_dens (:, :)
      Type (t_geometry), Intent (In) :: trgeo
      Complex (Kind=DEF_DBL_PREC) :: phi (:)
      Integer, Intent (In) :: ofsa, ofso, tna
      Type(t_Lop_mats) :: Lop(:)
    !!$ Local vars
      Integer :: ia, ofs, no
      Real (Kind=DEF_DBL_PREC) :: ldens(3)
      Type (t_atom_defenition), Pointer :: at
      ! this routine computes <psi|L_a|psi> for each atom
      ! where a = {x,y,z}, and adds the values found for up and down spins
      ofs = ofso + 1
      Do ia = 1, tna
        at => trgeo%atoms(ia+ofsa)%ptr
        no = at%nl **2
        ldens = 0d0
        ! first spin up
        ldens(1) = dot_product(phi(ofs:ofs+no-1), matmul(Lop(at%nl)%Lx(1:no, 1:no), phi(ofs:ofs+no-1))) 
        ldens(2) = dot_product(phi(ofs:ofs+no-1), matmul(Lop(at%nl)%Ly(1:no, 1:no), phi(ofs:ofs+no-1)))  
        ldens(3) = dot_product(phi(ofs:ofs+no-1), matmul(Lop(at%nl)%Lz(1:no, 1:no), phi(ofs:ofs+no-1)))  
        ! then spin down
        ldens(1) = ldens(1) + dot_product(phi(ofs+no:ofs+2*no-1), matmul(Lop(at%nl)%Lx(1:no, 1:no), phi(ofs+no:ofs+2*no-1))) 
        ldens(2) = ldens(2) + dot_product(phi(ofs+no:ofs+2*no-1), matmul(Lop(at%nl)%Ly(1:no, 1:no), phi(ofs+no:ofs+2*no-1))) 
        ldens(3) = ldens(3) + dot_product(phi(ofs+no:ofs+2*no-1), matmul(Lop(at%nl)%Lz(1:no, 1:no), phi(ofs+no:ofs+2*no-1)))  

        oam_dens(:, ia) = oam_dens(:, ia) + ldens(:)

        ofs = ofs + 2*no
      End Do

    End Subroutine oam_helper

  End Subroutine calc_oam_density

End Module
