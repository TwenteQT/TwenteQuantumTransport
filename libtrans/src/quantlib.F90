#include "math_def.h"
Module quantlib
    Implicit None
   !-----------discursh types-----------!
    Type planes_region
        integer, allocatable :: ind(:) ! plane index of each atom
        integer, allocatable :: napp(:) ! number of atoms per plane
        integer :: np, nat ! number of planes, atoms
        integer, allocatable :: order(:)
    end type planes_region

    Type planes_system
        type(planes_region) :: m !midle
        type(planes_region) :: l !left
        type(planes_region) :: r !right
        Double precision, allocatable :: th(:)
    end type planes_system

    Type spatialquantdis
        Real (Kind=DEF_DBL_PREC), Allocatable :: Lm (:,:), Rm (:,:)
        Integer :: alloc = 0
    End Type spatialquantdis

    Type confdisdata
        Type(spatialquantdis) :: xcur,ycur,zcur,ts,ts_av,spinmux,spinmuy,spinmuz,OAMdens,OAMdens_av
        Integer :: nat, ldim
        Integer :: alloc = 0
    End Type confdisdata

    Type confdisav
        Type(spatialquantdis) :: xcur,ycur,zcur
    End Type confdisav

    Type grid_info
        Integer,allocatable :: wt(:,:,:),ia(:,:,:,:),count(:,:,:)
        Real (Kind=DEF_DBL_PREC),Allocatable :: ctr(:,:,:,:)
    End Type grid_info

    Type gridquant
        Real (Kind=DEF_DBL_PREC), Allocatable :: Lm (:,:,:), Rm (:,:,:)
    End Type gridquant

    Type confgrid
        Type(gridquant) :: xcur,ycur,zcur,ts,ts_av,spinmux,spinmuy,spinmuz,OAMdens,OAMdens_av
        Real (Kind=DEF_DBL_PREC), Allocatable :: coord(:,:,:),wt(:,:)
        Integer :: alloc = 0
    End Type confgrid

    Type cell_system
        Real (Kind=DEF_DBL_PREC), Allocatable :: cell(:,:),pts(:,:),tr(:,:)
        Real (Kind=DEF_DBL_PREC), Allocatable :: zlr(:,:),wt(:)
        Integer :: info(3),sc(2),nl=0,nm=0,nr=0,nch=0
        Integer, Allocatable :: ind(:),pl(:)
        Real (Kind=DEF_DBL_PREC) :: base(2,2)
        Real (Kind=DEF_DBL_PREC), Allocatable :: dx(:),dy(:),dz(:)
    End Type cell_system

!-------------qofzlib types-----------------
    Type spatialquant ! general local quantity
        Integer :: alloc = 0
        Integer :: nl = 0, nr = 0, nm = 0 ! number of left/right/middle entries
        Integer :: ldim ! leading dimension
        Real (Kind=DEF_DBL_PREC), Allocatable :: m (:, :), l(:, :), r(:, :)
    End Type spatialquant
    
    Type spadirquant ! general local quantity that depends on direction of transport
        Integer :: alloc = 0
        Integer :: nl = 0, nr = 0, nm = 0 ! number of left/right/middle entries
        Integer :: ldim ! leading dimension
        Real (Kind=DEF_DBL_PREC), Allocatable :: Lir (:, :), Lt (:, :), Rir (:, :), Rt (:, :)
        Real (Kind=DEF_DBL_PREC), Allocatable :: Lm (:, :), Rm (:, :)
    End Type spadirquant
    
    Type confdata
      Integer :: writetorque, writeveloc, do_iatcu, do_iatlcu=0, do_oldcurr, do_OAMdens
        double precision :: cond, area, area_nm, base(2,2)
        Type (spatialquant) :: z, z_si, spindos, mdir
        Type (spadirquant) :: spinmu, ts, ts_pscr, xcur, ycur, zcur, cons, OAMdens, OAMdens_pscr
        Type (spadirquant) :: xcur_ps, ycur_ps, zcur_ps
        Type (spatialquant) :: SHang
    End Type confdata
    
    Type t_dirlist
        Type(t_dirlist), Pointer :: next
        Character(len=:), allocatable :: path
        Integer :: newgeom
    End Type
end Module quantlib
