!   $Id: definitions.f,v 1.39 2003/11/03 19:04:43 maciej Exp $
module definitions
!    This module contains the definitions of data structures used elsewhere
!    in the program and declarations od some constants.
!    prec   - kind value for double precision
!    lprec  - kind value for double double (extended double) precision
!    clmax  - max. size of the cluster for S^alpha calculations.
    Use bzgrid
    implicit none
!
  integer, parameter :: IR1=1,IR2=2,IR4=4,IR5=5,IR11=11,IR13=13,IW8=8,atomf=10
  integer, parameter :: prec = 8, lprec = 16
  integer:: IW6,MNCL=0 ! MNCL ..... PARAMETER: MAX. SIZE OF SMALL CLUSTERS FOR CALCULATION OF TB-STRUCTURE CONSTANTS
                                ! the only parameter that needs to be allocated dynamically but after rewriting structre const stuff
!
!   real (kind=prec), parameter :: img_max = 1.0d-4, min_dist = 0.01d0
!
 type params
     real (kind=prec)::difmc,qmixc,difms,qmixs,diam,alfa,beta,w0am,thresh
     integer::ivac,nl,ns,irel,ivxc,lident,nfupr,lmax,nharm, rwgamma,ne
     integer::niter,nitera,nuh,nsirk,nmits,npp1,npp2,nlsq,nblsq,mtbcl,pair,nam,icont
     integer::lnonc
     real (kind=prec)::tol,wtol 
     Type (t_bz_opts) :: bz
  end type params

  type atom_definition ! scattering region atom
     character(len=16):: otxta,otxta_in
     integer :: nm,nszrad,numcor
     integer, pointer :: ncor(:),lcor(:),nobc(:),isheny(:,:)
     real (kind=prec) :: con,valz,az,ws,wsav,dipmom,amgmom,chatra,vmad,dmad,acte
     real (kind=prec),pointer:: pot(:,:),ecor(:,:),phi(:,:,:),phid(:,:,:),phidd(:,:,:),rhocor(:,:),rhoval(:,:)
     real (kind=prec),pointer:: eny(:,:),ppc(:,:), ppd(:,:), ppq(:,:),ppp(:,:),dny(:,:),finy(:,:),finyd(:,:)
     real (kind=prec),pointer:: ynam(:,:),fpot(:,:),emdi0(:,:),emdi1(:,:),emdi2(:,:)
     real (kind=prec),pointer:: emof00(:,:),emof10(:,:),emof01(:,:),emof11(:,:),emof20(:,:),emof02(:,:)
  end type atom_definition

  type bulk_atom_definition ! bulk atom
     character(len=16):: otxta
     real (kind=prec) :: ws,wsav,az,ef,con
     real (kind=prec),pointer:: eny(:,:), ppc(:,:), ppd(:,:), ppq(:,:), ppp(:,:), dny(:,:),finy(:,:), finyd(:,:)
  end type bulk_atom_definition

  type atoms_set  !   our system
     integer::np,nb,nums,num,bnum,bns,vnum,vns,nmtr,nrmax              !num,nlmax,nrmax,nsmax
     real (kind=prec):: cutrat,scx(3),vbr(2,2),vbg(2,2),vtrv(3),btrv(3),bef,vef,ef,dba,bwst,vwst,abws,avws,romega,gomega
     real (kind=prec), pointer :: pos(:,:,:),vpos(:,:),bpos(:,:),qscr(:),aws(:,:)
     integer,pointer:: nc(:)
     type (atom_definition),pointer:: at(:) 
     type (bulk_atom_definition),pointer:: bat(:),vat(:) 
  end type atoms_set

  type complex_nodes  ! complex nodes & weights
     complex(kind=prec),pointer :: zcn(:), zcw(:) 
  end type complex_nodes

  type madelung
     real(kind=prec),pointer:: mmm(:,:),mmd(:,:),mdm(:,:),mdd(:,:),mdb(:),ddb(:)
  end type madelung

    type sharm  ! spherical harminics stuff
     real(kind=prec), pointer:: cnh(:,:),gfrh(:,:,:)
  end type sharm

  type andmix
     real(kind=prec), pointer:: dxp(:,:),dfp(:,:),sp(:),xl(:),fl(:),xn(:),voma(:,:) 
  end type andmix
end module definitions

module str_const
  use definitions
  complex(kind=prec),pointer::  zski(:,:,:,:),zsko(:,:,:,:),zskf(:,:,:,:),zbski(:,:,:), &
                              & zbsko(:,:,:),zbskf(:,:,:),zvski(:,:,:),zvsko(:,:,:),zvskf(:,:,:)
end module str_const

module green_func
  use definitions
  complex(kind=prec),pointer::zbgam(:,:,:,:,:),zvgam(:,:,:,:,:),zgfg(:,:,:,:,:),zcpf(:,:,:,:,:), &
                             & zomg(:,:,:,:,:),zcagf(:,:,:,:,:)
end module green_func

