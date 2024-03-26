!!$ Module for calculation localized structure matrix in real space
!!$ code extracted from Stuttgart-NMTO v4.7, modified by
!!$ Zhe Yuan@UTwente
module omta_strrs
  use omta_defs
  use geometry_module
  use logging
  implicit none

contains

  subroutine calcsc(atoms,geo,kappa2,scalpha)
  implicit none
  type(t_atoms_set_EMTO)   :: atoms
  type(t_geometry_EMTO)    :: geo
  type(t_omta_struct) :: scalpha
  real(kind=prec)     :: kappa2

  scalpha%base=geo%base
  scalpha%perp_trans=geo%perp_trans
  call gtdims(geo, scalpha%neighm, scalpha%ns)
  allocate(scalpha%s(scalpha%ns),scalpha%iax(9,geo%num,scalpha%neighm),scalpha%npr(geo%num))
  call strrs(atoms,geo,scalpha,kappa2)

  end subroutine calcsc 

  subroutine strrs(atoms,geo,scalpha,ki)
!- Calculates localized structure constants for all spheres
!-----------------------------------------------------------------------
!i Inputs:
!i   alat  :length scale
!i   avw   :average Wigner-Seitz sphere radius
!i   nbas  :number of atoms in the basis
!i   nclass:number of different classes
!i   nl    :number of l's
!i   plat  :primitive lattice vectors (scaled by alat)
!i   rmaxs :cluster radius (scaled by avw)
!o Outputs:
!o   s     :screened structure constant matrix
!r Remarks:
!r These subroutines calculate tight-binding structure
!r constants and sdot for all inequivalent atoms in a
!r lattice.  Structure constants S for the various inequivalent atoms
!r are strung together in one long vector s of length npair, the
!r total number of inequivalent pairs of atoms with nonvanishing
!r structure constants.  Finally S is stored in a onedimensional array.
!r
!r iax contains information for the mapping of S, as stored in the
!r machine to the structure constant array S(R'L',RL).
!r iax is an integer array of the form iax(9,neighm,nbas).
!r R and R' are separated by some primitive translation vectors,
!r plus some difference of basis vectors within the unit cell.
!r iax(1) is an index to which of the atoms in bas R corresponds to;
!r iax(2) is an index to which of the atoms in bas R' corresponds to;
!r iax(3..5) are the number of primitive lattice translation vectors
!r           separating the unit cells corresponding to R and R'.
!r iax(6)  is nlsqr(R) ; nlsqr=(lmx+1)**2
!r iax(7)  is nlsqr(R')
!r iax(8)  is an offset to R'-block in s0a
!r iax(9)  is an offset to RR'-block in s
!r
!r STRRS returns with npr(ibas) which contains the number of atoms
!r within the cluster surrounding the ith atom.
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   type(t_atoms_set_EMTO)   :: atoms
   type(t_geometry_EMTO)    :: geo
   type(t_omta_struct) :: scalpha
   real(kind=prec)     :: ki
! Local variables:
   logical          :: lsymst=.true.
   integer          :: nlmax
   real(kind=prec), allocatable :: cg    (:)
   integer        , allocatable :: jcg   (:)
   integer        , allocatable :: indxcg(:)
   real(kind=prec), allocatable :: cy    (:)

!!$prepare Clebsch-Gordon and spherical harmonics to be used 
!!$in structure constant, the maximum angular-momentum is needed
   nlmax=geo%nlmax
   allocate(cg(nlmax**5/2+1),jcg(nlmax**5/2+1),cy(nlmax*nlmax*4),&
            indxcg(nlmax*nlmax*(nlmax*nlmax+1)/2+1))
   call sylmnc(cy,2*nlmax-1)
   call scg(nlmax-1,cg,indxcg,jcg)
   call nghbr0(geo,scalpha)
   call salph0(atoms,geo,scalpha,ki,cg,cy,indxcg,jcg,lsymst)
   deallocate(cg,jcg,cy,indxcg)
   return
   end subroutine strrs

   subroutine salph0(atoms,geo,scalpha,ki,cg,cy,indxcg,jcg,lsymst)
!- Calculates localized structure constants for all spheres
   implicit none
!  Passed variables:
   type(t_atoms_set_EMTO)   :: atoms
   type(t_geometry_EMTO)    :: geo
   type(t_omta_struct) :: scalpha
   real(kind=prec)     :: ki
   real(kind=prec)     :: cg(*)
   real(kind=prec)     :: cy(*)
   integer             :: indxcg(*)
   integer             :: jcg(*)
   logical             :: lsymst
!  integer          ,intent(in ) :: nbas
!  integer          ,intent(in ) :: nclass
!  integer          ,intent(in ) :: nlmax
!  integer         , intent(in ) :: iax(9,nbas,*)
!  integer          ,intent(in ) :: npr(nbas)
!  integer          ,intent(in ) :: ns
!  real(kind=prec)  ,intent(in ) :: avw
!  real(kind=prec)  ,intent(in ) :: plat(3,3)
!  real(kind=prec)  ,intent(out) :: s(ns)
!  type(tclass)                  :: aclass(nclass)
!  type(tsite)                   :: site(nbas)
!  Local variables:
   integer          :: nbas,jbas,ndiml,nlsqr,nlmax,jc,npr
   real(kind=prec)  :: kap2i,avw
   real(kind=prec) ,allocatable :: sc(:,:)

   nbas=geo%num
   nlmax=geo%nlmax
   avw  =geo%dawsr
   nlsqr=geo%nlmax**2
   kap2i=ki
   ! ----- calculate screening constants and transformation matrix
   call mkalph(kap2i,nlmax,avw,atoms)
   call mktral(kap2i,nlmax,avw,atoms)
   ! ------- Make the structure constants for all sites
   do jbas = 1, nbas
      npr=scalpha%npr(jbas)
      ndiml=0
      do jc=1,npr
        ndiml=ndiml+scalpha%iax(7,jbas,jc)
      enddo
      ! --------- Work array for S0_RR',salpha
      allocate(sc(ndiml,nlmax*nlmax))
      ! --------- Calculate salpha 
      call salph1(avw,ndiml,nlmax,nlsqr,jbas,cg,cy,indxcg,jcg,kap2i,geo,scalpha,atoms,sc)
      ! --------- Append salpha and/or salpha-dot to structure constant s
      call addtos(jbas,ndiml,sc,scalpha)
      ! --------- convergence check for LMTO-ASA structure constants
!     call ccSrs(iax,jbas,nbas,site,nclass,aclass,npr,s)
      deallocate(sc)
   enddo
   ! ------- Symmetrize structure constants
   if (lsymst) call symstr(nbas,scalpha)

303 format(' SALPH0: cluster',i3,', atom ',a4,i3,' sites; lmx=',i1,'; ndiml=',i4)
   return
   end subroutine salph0

!  subroutine salph1(avw,cg,cy,iax,indxcg,jbas,jcg,kap2,nclass,nbas,aclass,site,&
!                    ndiml,nlmax,nlsqr,npr,plat,sc)
   subroutine salph1(avw,ndiml,nlmax,nlsqr,jbas,cg,cy,indxcg,jcg,kap2,geo,scalpha,atoms,sc)
!- Calculate single-kappa s for one atom
   implicit none
   ! Passed variables:
   integer             :: jbas
   integer             :: ndiml
   integer             :: nlmax
   integer             :: nlsqr
   integer             :: indxcg(*)
   integer             :: jcg(*)
   real(kind=prec)     :: cg(*)
   real(kind=prec)     :: cy(*)
   real(kind=prec)     :: avw
   real(kind=prec)     :: kap2
   real(kind=prec)     :: sc(ndiml,nlmax*nlmax) 
   type(t_geometry_EMTO)    :: geo
   type(t_omta_struct) :: scalpha
   type(t_atoms_set_EMTO)   :: atoms
   ! Local variables:
   integer                        :: ierr,jlm,nlsqrj,npr
   integer           ,allocatable :: kpvt(:)
   real(kind=prec)   ,allocatable :: diag(:,:)
   real(kind=prec)   ,allocatable :: s0a(:), s0j(:)
   complex(kind=prec)   ,allocatable :: s0c(:)

   nlsqrj=scalpha%iax(6,jbas,1)
   npr=scalpha%npr(jbas)
   allocate(kpvt(ndiml),diag(nlsqr,npr),s0a(ndiml*(ndiml+1)/2))
   if(kap2>0.d0) then
     allocate (s0J(ndiml*(ndiml+1)/2),s0c(ndiml*(ndiml+1)/2)) 
   else
     allocate (s0J(1))
   end if ! 
   ! --- s0a -> alpha^-1-S^0
   diag = 0.d0
   call mkdiag(diag,jbas,nlsqr,npr,geo,scalpha)
   s0a = 0.0d0 ; s0J = 0.d0
   call maksp0(avw,cg,cy,diag,indxcg,jbas,jcg,kap2,nlsqr,npr,geo,scalpha,s0a,s0j,1)
   ! --- Add term delta_RL-,R'L' to rhs
   sc(:,:) = 0.0d0
   do jlm=1,nlsqrj
      sc(jlm,jlm)=1.0d0
   enddo
   ! --- Beta for inner block from back substitution
   ! --- sc -> (alpha^-1-S^0)^-1

   if(kap2>0.d0) then
     s0c = DCMPLX(s0a,s0J)
!!$     call invert_complex_packed(ndiml,nlmax*nlmax,s0c,sc)
     call zspfa(s0c,ndiml,kpvt,ierr)
   else
!!$     call invert_real_packed(ndiml,nlmax*nlmax,s0a,sc)
     call dspfa(s0a,ndiml,kpvt,ierr)
   end if ! 

   if (ierr /= 0)then
      write(*,*) ' SALPH1: matrix singular'
      stop
   endif

   sc(:,:) = 0.0d0
   do jlm=1,nlsqrj
      sc(jlm,jlm)=1.0d0
   enddo
   ! --- Beta for inner block from back substitution
   do jlm = 1, nlsqrj
     if(kap2>0.d0) then
       call zspsl(s0c,ndiml,kpvt,sc(1,jlm))
     else
       call dspsl(s0a,ndiml,kpvt,sc(1,jlm))
     end if
   end do


   deallocate(s0J)
   if (kap2>0.d0) deallocate(s0c)

   call scals(jbas,ndiml,nlmax,npr,geo,scalpha,atoms,sc)

   deallocate(kpvt,diag,s0a)
   return
   end subroutine salph1

!  subroutine mkdiag(diag,iax,jbas,nbas,site,nclass,aclass,nl,nlsqr,npr)
   subroutine mkdiag(diag,jbas,nlsqr,npr,geo,scalpha)
!-- Makes diagonal part added to S^0 or Sdot^0
   implicit none
!  Passed variables:
   integer             :: npr
   integer             :: jbas
   integer             :: nlsqr
   real(kind=prec)     :: diag(nlsqr,npr)
   type(t_geometry_EMTO)    :: geo
   type(t_omta_struct) :: scalpha
!  Local variables:
   integer :: ic,ilm,ipr,l,nlsqri

   do ipr=1,npr
      nlsqri=scalpha%iax(7,jbas,ipr)
      if(ipr /= npr) ic=scalpha%iax(2,jbas,ipr)
      do ilm=1,nlsqri
         if (ipr == npr) stop 'mkdiag error!'
         l=ll(ilm)
         diag(ilm,ipr)=1.d0/geo%atoms(ic)%ptr%alpha(l)
      enddo
   enddo
   return
   end subroutine mkdiag

!  subroutine maksp0(avw,nbas,site,cg,cy,diag,iax,indxcg,jbas,jcg,kap2, &
!                    nlsqr,npr,plat,sp0,s0J)
   subroutine maksp0(avw,cg,cy,diag,indxcg,jbas,jcg,kap2,nlsqr,npr,geo,scalpha,sp0,s0J,iop)
!- Make real space unscreened structure constants OR its energy derivative
! ----------------------------------------------------------------------
!o Outputs:
!o   sp0    :alpha^-1 - S0, stored in packed form,
!r Remarks:
!r   Constructs structure matrix for neighbor table iax at energy kap2.
!r   A diagonal term diag is added.
!r   sp0 is assumed initialized to zero.
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer             :: indxcg(*)
   integer             :: jcg(*)
   integer             :: jbas
   integer             :: npr
   integer             :: nlsqr
   real(kind=prec)     :: cg(*)
   real(kind=prec)     :: cy(*)
   real(kind=prec)     :: avw
   real(kind=prec)     :: kap2
   real(kind=prec)     :: diag(nlsqr,npr)
   real(kind=prec)     :: sp0(*),s0J(*)
   type(t_geometry_EMTO)    :: geo
   type(t_omta_struct) :: scalpha
   integer             :: iop
! Local variables:
   logical         :: islead
   integer         :: ipr,jpr,ij,ilm,jlm,idxd,lmax,offs,nlsqri,nlsqrj
   integer         :: i1,i2,j1,j2,k1,k2
   real(kind=prec) :: dr(3),dsqr,plat(3,3),coord1(3),coord2(3)
   real(kind=prec),allocatable :: s0(:,:),sJ(:,:)

   if(geo%num==geo%l_transnum)then
     islead=.true.
   else
     islead=.false.
   endif

   plat=0.d0
   plat(1:2,1)=geo%base(1:2,1)/avw
   plat(1:2,2)=geo%base(1:2,2)/avw
   plat(1:3,3)=geo%perp_trans(1:3)/avw

   idxd = 0
   do ipr = 1, npr
     nlsqri = scalpha%iax(7,jbas,ipr)
     coord1 = geo%atoms(scalpha%iax(2,jbas,ipr))%coord/avw
     i1=scalpha%iax(3,jbas,ipr)
     j1=scalpha%iax(4,jbas,ipr)
     k1=scalpha%iax(5,jbas,ipr)
! ----- Add diag to sp0
! ----- idxd(idxd+3)/2 = idxd(idxd+1)/2+idxd is offset to
! ----- diagonal matrix element s0(idxd,idxd):
     do ilm = 1, nlsqri
       sp0(1+idxd*(idxd+3)/2) = diag(ilm,ipr)
       idxd = idxd+1
     end do ! ilm
     do jpr = ipr+1, npr
       nlsqrj = scalpha%iax(7,jbas,jpr)
       coord2 = geo%atoms(scalpha%iax(2,jbas,jpr))%coord/avw
       i2=scalpha%iax(3,jbas,jpr)
       j2=scalpha%iax(4,jbas,jpr)
       k2=scalpha%iax(5,jbas,jpr)

       if(.not.islead) then
         if(k1*k2.lt.0) then
           write(*,*) 'hopping error in maksp0!'
           stop
         endif
         if(k1+k2.lt.0)then
           plat(1:3,3)=geo%l_perp_trans(1:3)/avw
         elseif(k1+k2.gt.0)then
           plat(1:3,3)=geo%r_perp_trans(1:3)/avw
         else
           plat(1:3,3)=0.d0
         endif
       endif

       dsqr=drr2(plat,coord1,coord2,i2-i1,j2-j1,k2-k1,dr)

       allocate (s0(nlsqri,nlsqrj))
       lmax = ll(nlsqri)+ll(nlsqrj)
       if (jpr /= npr) then
         call mstrx2(kap2,dr,nlsqri,nlsqrj,lmax,cg,indxcg,jcg,cy,s0,iop)
       else if (nlsqrj /= 0) then
         write(*,*) 'maksp0 strange'
         stop
       end if
       do jlm = 1, nlsqrj
! --------- Calculate packed equivalent of address to:
! ---------       s0_{R'L',RL}
! ---------  <==> s0(R_ipr,ilm;R_jpr,jlm)
! ---------  <==> sp0(iax(8,jbas,ipr)+ilm,iax(8,jbas,jpr)+jlm)
         offs = scalpha%iax(8,jbas,jpr)+jlm-1
         offs = (offs*(offs+1))/2 + scalpha%iax(8,jbas,ipr)
         do ilm = 1, nlsqri
           ij = offs+ilm
           sp0(ij) = s0(ilm,jlm)
         end do ! ilm
       end do ! jlm

       deallocate (s0)
     end do ! jpr
   end do ! ipr

   if(kap2<=0.d0) return

! =============== case kap2 > 0 ===========
   do ipr = 1, npr-1
     nlsqri = scalpha%iax(7,jbas,ipr)
     coord1 = geo%atoms(scalpha%iax(2,jbas,ipr))%coord/avw
     i1=scalpha%iax(3,jbas,ipr)
     j1=scalpha%iax(4,jbas,ipr)
     k1=scalpha%iax(5,jbas,ipr)
     do jpr = ipr, npr-1
       nlsqrj = scalpha%iax(7,jbas,jpr)
       coord2 = geo%atoms(scalpha%iax(2,jbas,jpr))%coord/avw
       i2=scalpha%iax(3,jbas,jpr)
       j2=scalpha%iax(4,jbas,jpr)
       k2=scalpha%iax(5,jbas,jpr)

       if(.not.islead) then
         if(k1*k2.lt.0) then
           write(*,*) 'hopping error in maksp0!'
           stop
         endif
         if(k1+k2.lt.0)then
           plat(1:3,3)=geo%l_perp_trans(1:3)/avw
         elseif(k1+k2.gt.0)then
           plat(1:3,3)=geo%r_perp_trans(1:3)/avw
         else
           plat(1:3,3)=0.d0
         endif
       endif

       dsqr=drr2(plat,coord1,coord2,i2-i1,j2-j1,k2-k1,dr)
!      dsqr = drr2(plat,site(iax(2,jbas,ipr))%coord/avw,site(iax(2,jbas,jpr))%coord/avw, &
!                  iax(3,jbas,jpr)-iax(3,jbas,ipr),                    &
!                  iax(4,jbas,jpr)-iax(4,jbas,ipr),                    &
!                  iax(5,jbas,jpr)-iax(5,jbas,ipr),dr)

       allocate (sJ(nlsqri,nlsqrj))
       sJ = 0.d0
       lmax = ll(nlsqri) + ll(nlsqrj)
       call mstrx3(kap2,dr,nlsqri,nlsqrj,lmax,cg,indxcg,jcg,cy,sJ,iop)
       call renorm_sJ(nlsqri,nlsqrj,sJ,kap2)
       if(ipr /= jpr) then
         do jlm = 1, nlsqrj
           offs = scalpha%iax(8,jbas,jpr)+jlm-1
           offs = (offs*(offs+1))/2 + scalpha%iax(8,jbas,ipr)
           do ilm = 1, nlsqri
             ij = offs + ilm
             s0J(ij) = sJ(ilm,jlm)
           end do ! ilm
         end do ! jlm
       else
         do jlm = 1, nlsqrj
           offs = scalpha%iax(8,jbas,jpr)+jlm-1
           offs = (offs*(offs+1))/2 + scalpha%iax(8,jbas,ipr)
           ij = offs + jlm
           s0J(ij) = sJ(jlm,jlm)
         end do ! jlm
       end if ! ipr /= jpr
       deallocate (sJ)
     end do ! jpr
   end do ! ipr

   return
   end subroutine maksp0

   real(kind=prec) function drr2(plat,bas1,bas2,i,j,k,dr)
!- Calculates the vector connecting two sites in a solid
! ----------------------------------------------------------------------
!i Inputs:
!i   plat  :primitive lattice vectors(scaled with avw)
!i   bas1  :basis vector of first site
!i   bas2  :basis vector of second site
!i   i,j,k :the number of primitive lattice vectors separating sites
!o Outputs:
!o   dr    :connecting vector bas2 - bas1
!o   drr2  :square of the length of this vector
!r Remarks:
!r   Using the TB package and a table of indices iax, the connecting
!r   vector and the square of the distance is obtained by:
!r   rsqr=drr2(plat,bas(1,iax(1)),bas(1,iax(2)),iax(3),iax(4),iax(5),dr)
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer         , intent(in ) :: i
   integer         , intent(in ) :: j
   integer         , intent(in ) :: k
   real(kind=prec) , intent(in ) :: plat(3,3)
   real(kind=prec) , intent(in ) :: bas1(3)
   real(kind=prec) , intent(in ) :: bas2(3)
   real(kind=prec) , intent(out) :: dr(3)
! Local variables:
   integer ix

   drr2 = 0.d0
   do ix = 1, 3
     dr(ix)=bas2(ix)-bas1(ix)+plat(ix,1)*i+plat(ix,2)*j+plat(ix,3)*k
     drr2 = drr2 + dr(ix)*dr(ix)
   enddo

   return
   end function drr2

   subroutine scals(jbas,ndiml,nl,npr,geo,scalpha,atoms,sc)
!- Scales S
   implicit none
! Passed variables:
   integer             :: jbas
   integer             :: ndiml
   integer             :: nl
   integer             :: npr
   real(kind=prec)     :: sc(ndiml,nl*nl)
   type(t_geometry_EMTO)    :: geo 
   type(t_omta_struct) :: scalpha
   type(t_atoms_set_EMTO)   :: atoms
! Local variables:
   integer :: ibas,ic,ii,ilm,ipr,kbas,kc,kk,klm,kpr,l,li,lk,nlsqri,nlsqrk,nclass
   real(kind=prec)  :: dt
   real(kind=prec),allocatable :: wkscl(:,:,:)

   nclass=atoms%num
   allocate(wkscl(2,0:nl-1,0:nclass))
   wkscl(:,:,:)=0.d0
   do ic=1,nclass
      do l=0,atoms%at(ic)%lmx
         wkscl(1,l,ic)=atoms%at(ic)%t(1,l)/atoms%at(ic)%t(3,l)
         wkscl(2,l,ic)=1.d0/atoms%at(ic)%t(3,l)
      enddo
   enddo
   ii = 1
   do ipr = 1, npr
      ibas=scalpha%iax(2,jbas,ipr)
      ic=geo%atoms(ibas)%iclass
      if (ipr == npr) ic=0
      nlsqri=scalpha%iax(7,jbas,ipr)
      do ilm = 1, nlsqri
         li=ll(ilm)
         kk = 1
         do kpr = 1, 1
            kbas=scalpha%iax(2,jbas,kpr)
            kc=geo%atoms(kbas)%iclass
            if (kpr == npr) kc=0
            nlsqrk=scalpha%iax(7,jbas,kpr)
            do klm = 1, nlsqrk
               lk=ll(klm)
               dt=atoms%at(kc)%t(1,lk)*atoms%at(kc)%t(4,lk)-atoms%at(kc)%t(2,lk)*atoms%at(kc)%t(3,lk)
               sc(ii,kk)= wkscl(2,li,ic)*sc(ii,kk)*wkscl(2,lk,kc)*dt
               if (kk == ii) sc(kk,kk)=sc(kk,kk)+wkscl(1,lk,kc)
               kk = kk+1
            enddo
         enddo
         ii = ii+1
      enddo
   enddo
   deallocate(wkscl)

   return
   end subroutine scals

!  subroutine addtos(iax,jbas,nbas,ndiml,npr,s,sc)
   subroutine addtos(jbas,ndiml,sc,scalpha)
!- Appends structure constant vector for single sphere to s
! ----------------------------------------------------------------------
!i Inputs:
!i   jbas  :center of current cluster
!i   sc    :structure constants in a two-dimensional array ndiml by nl*2
!o Outputs:
!o   s     :sc is appended
!r Remarks:
!r   Converts s to format s(*), the indices corresponding
!r   to S_ R'L', RL.  Indices RL correspond to sphere and L channel
!r   of the head of the screened orbital and R'L' to the sphere and
!r   L channel of which S is the one-center expansion for the screened
!r   function.
!r   The ordering of the indices is L',L,R',R i.e., L' is the fastest
!r   changing, and R the slowest changing variable.
!r   But, attention S cannot be writen as an
!r   array S(nl*nl,nl*nl,neighm,neighm). Instead:
!r        S(i,j) = S(offs)
!r        with: offs=iax(9,jbas,ipr)+(jlm-1)*nlsqri+ilm
!r        and  j: head index      i: tail index
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer             :: jbas
   integer             :: ndiml
   real(kind=prec)     :: sc(ndiml,*)
   type(t_omta_struct) :: scalpha
! Local variables:
   integer :: npr,i,ij,ipr,ilm,jlm,nlsqri,nlsqrj,offs,offsc

   npr=scalpha%npr(jbas)
   do ipr = 1, npr
      nlsqrj=scalpha%iax(6,jbas,ipr)
      nlsqri=scalpha%iax(7,jbas,ipr)
      offsc =scalpha%iax(8,jbas,ipr)
      do jlm = 1, nlsqrj
         offs  =scalpha%iax(9,jbas,ipr)+(jlm-1)*nlsqri
         do ilm = 1, nlsqri
            ij = offs  + ilm
            i  = offsc + ilm
            scalpha%s(ij)=sc(i,jlm)
         enddo
      enddo
   enddo

   return
   end subroutine addtos

!  subroutine symstr(iax,nbas,npr,s)
   subroutine symstr(nbas,scalpha)
!- Symmetrize structure constants
! ----------------------------------------------------------------------
! Remarks:
!  Force s(r1,l1,r2,l2) = s(r2,l2,r1,l1) by averaging the two.
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer             :: nbas
   type(t_omta_struct) :: scalpha
! Local variables:
   integer :: i1,i2,i3,ibas,ij,ilm,ipr,j1,j2,j3,jbas,  &
              ji,kbas,jlm,jpr,nlsqri,nlsqrj,offsi,offsj
   real(kind=prec) :: temp,tol
   logical :: lfirst
   character*9,parameter :: sordot='STR      '

   lfirst=.true.
   tol=0.10d0

   do jbas=1,nbas
     do jpr=1,scalpha%npr(jbas)-1
       ibas  =scalpha%iax(2,jbas,jpr)
       j1    =scalpha%iax(3,jbas,jpr)
       j2    =scalpha%iax(4,jbas,jpr)
       j3    =scalpha%iax(5,jbas,jpr)
       nlsqrj=scalpha%iax(6,jbas,jpr)
       offsj =scalpha%iax(9,jbas,jpr)
       do ipr=1,scalpha%npr(ibas)-1
         kbas  =scalpha%iax(2,ibas,ipr)
         i1    =scalpha%iax(3,ibas,ipr)
         i2    =scalpha%iax(4,ibas,ipr)
         i3    =scalpha%iax(5,ibas,ipr)
         nlsqri=scalpha%iax(6,ibas,ipr)
         offsi =scalpha%iax(9,ibas,ipr)
         if(kbas == jbas.and.i1 == -j1.and.i2 == -j2.and.i3 == -j3) then
           do jlm = 1, nlsqrj
             do ilm = 1, nlsqri
               ji=offsj+(jlm-1)*nlsqri+ilm
               ij=offsi+(ilm-1)*nlsqrj+jlm
               temp = (scalpha%s(ji)+ scalpha%s(ij))*0.5d0
!              if(ABS(scalpha%s(ij)-temp) > tol)then
!                 if (lfirst) write(*,400)
!                 lfirst=.false.
!                 write(*,401)sordot,jbas,jpr,jlm,ibas,ipr,&
!                                    ilm,scalpha%s(ij),scalpha%s(ji)
!                 if(ABS(scalpha%s(ij)-temp) > 2.d0*tol) then
!                   write(*,402)4.d0*tol
!                 endif ! ABS(s(ij)-temp) > 2.d0*tol
!              endif ! ABS(s(ij)-temp) > tol
               scalpha%s(ij) = temp
               scalpha%s(ji) = temp
             enddo ! ilm
           enddo ! jlm
           goto 100
         endif ! kbas == jbas.and.i1 == -j1.and.i2 == -j2.and.i3 == -j3
       enddo ! ipr
       scalpha%s(scalpha%iax(9,jbas,jpr)+1: &
      & scalpha%iax(9,jbas,jpr)+scalpha%iax(6,jbas,jpr)*scalpha%iax(7,jbas,jpr))=0.d0
   100  continue
     enddo ! jpr
   enddo ! jbas

300 format(/' SYMSTR: Symmetrize ',a)
400 format(' SYMSTR: Large difference between S_ij and S_ji:',        &
          /9x,'TYPE jbas jpr jlm   ibas ipr ilm     s_ij      s_ji',  &
          /9x,53('-'))
401 format(9x,a4,3(1x,i3),3x,3(1x,i3),2f10.5)
402 format(' SYMSTR: difference larger than',f6.3,                    &
           '|        increase RMAXS or delete small empty spheres.$')

   return 
   end subroutine symstr

   subroutine ccSrs(iax,jbas,nbas,site,nclass,aclass,npr,s)
!----------------------------------------------------------!
!    convergence check for LMTO-ASA structure constants    !
!    date: 03.07.02 author: mkorotin@optics.imp.uran.ru    !
!----------------------------------------------------------!
! The programm checks the simple charge sum-rule involving !
! the ss-elements [Varenna Eq.(123) and Kanpur Eq.(3.34)]: !
!             R_max            alpha                       !
!           1+ Sum  aplpha  * S       <= tolerance         !
!              R=R'       R0   R0,R'0                      !
!----------------------------------------------------------!
   implicit none
   integer          ,intent(in) :: jbas,nbas,nclass,npr(*),iax(9,nbas,*)
   real(kind=prec)  ,intent(in) :: s(*)
   type(tclass)                 :: aclass(nclass)
   type(tsite)                  :: site(nbas)
   integer          :: jpr
   real(kind=prec)  :: cc, tolerance=1.0d-3

   if(jbas == 1) then
     write(*,'(/52h Convergence check for LMTO-ASA structure constants:)')
   else
     if(site(jbas)%iclass==site(jbas-1)%iclass) return
   end if

   cc = 0.d0
   do jpr = 1, npr(jbas)-1
     cc = cc + aclass(site(iax(2,jbas,jpr))%iclass)%alpha(0)*s(iax(9,jbas,jpr)+1)
   end do

   cc=ABS(1.d0+cc)
   write(*,'(2x,7h class ,a4,26h convergency of s-channel=,ES12.3)') &
         aclass(site(jbas)%iclass)%label, cc

   if(cc > tolerance) then
     write(*,'(/a/a)')' Convergency of real space structure constants &
                    &is bad',' Recommendation: increase NDIMIN by 100 &
                    &in CTRL file.'
     stop
   endif

   end subroutine ccSrs

   subroutine nghbr0(geo,scalpha)
!- Create a table of all neighbors within spec'd range of a spec'd atom
! ----------------------------------------------------------------------
!o Outputs:
!o   iax   :structure-constant information relating to s (see remarks)
!o   npr   :returns number of neighbors within rmaxs
!r Remarks:
!r   Creates a neighour list for a specified atom, generating iax
!r   Table is sorted according to increasing length from central atom.
!o   for specified atom; it contains:
!r    iax(1): jbas = center of cluster
!r    iax(2): ibas = atom in cluster
!r    iax(3): i
!r    iax(4): j
!r    iax(5): k
!r    iax(6): nlsqr(jbas)
!r    iax(7): nlsqr(ibas)
!r    iax(8): offset to sc or s0
!r    iax(9): offset to s
!r   Sort of iax is according to following priority:
!r   (1)distance, (2)ibas, (3)i, (4)j, (5)k ...
!r   iax(8)  and iax(9) are used to have a faster access to elements
!r   of s, sc, and S0.
!r   if sc is declared sc(ndiml,nlsqr) and with ipr=R' the
!r   L'L-block of S_{R'0} starts at sc(iax(8,jbas,ipr)+1,1)
!r
!r   if s is declared s(*) and with ipr=(R',R) the L'L-block of
!r   S_{R'R} starts at s(iax(9,jbas,ipr)+1)
!r   and ends at s(iax(9,jbas,ipr+1))
!r 
!r  Attention: NEIGHM0 must perform exactly the same calculation 
!r  as subroutime GTDIMS does. GTDIMS is first called to determine 
!r  the size of IAX and S. In NEIGHM0 iax is the populated with 
!r  values, which themselves determine the size of S.   
! ----------------------------------------------------------------
   implicit none
! Passed variables:
   type(t_omta_struct) :: scalpha
   type(t_geometry_EMTO)    :: geo
! Local variables:
   logical :: islead
   integer :: i,i1,i2,i3_l,i3_r,ibas,im1,ipr,   &
              j,jbas,jc,k,nlsqrj,inx,nbas
   real(kind=prec) :: dsqr,rmx,rmxsqr,avw
   real(kind=prec) :: plat(3,3),coord1(3),coord2(3),dr(3)
   integer,allocatable :: iwk(:,:)

   rmx=geo%cutrat
   rmxsqr = rmx*rmx

   nbas=geo%num
   avw=geo%dawsr

   if(geo%num==geo%l_transnum)then
     islead=.true.
   else
     islead=.false.
   endif

   i1=geo%ntrpar
   i2=geo%ntrpar
   i3_l=geo%l_ntrperp
   i3_r=geo%r_ntrperp

   allocate (iwk(5,(2*i1+1)*(2*i2+1)*(i3_l+i3_r+1)*nbas))

   plat=0.d0
   plat(1:2,1)=geo%base(1:2,1)/avw
   plat(1:2,2)=geo%base(1:2,2)/avw
   plat(1:3,3)=geo%perp_trans(1:3)/avw

   jc=0
   do jbas = 1,nbas
     nlsqrj=(geo%atoms(jbas)%ptr%lmx+1)**2
     ipr=0
     coord1(1:3)=geo%atoms(jbas)%coord(1:3)/avw
     do  ibas = 1, nbas
       coord2(1:3)=geo%atoms(ibas)%coord(1:3)/avw
!------- Sweep lattice translations to find all neighbors within rmaxs.
       do i = -i1, i1
         do j = -i2, i2
           do k = -i3_l, i3_r

             if((.not.islead).and.(k.ne.0)) then
               if(k.lt.0)then
                 if(ibas.gt.geo%l_transnum) cycle
                 plat(1:3,3)=geo%l_perp_trans(1:3)/avw
               elseif(k.gt.0)then
                 if(ibas.le.geo%num-geo%r_transnum) cycle
                 plat(1:3,3)=geo%r_perp_trans(1:3)/avw
               else
                 plat(1:3,3)=0.d0
               endif
             endif

             dsqr = drr2(plat,coord1,coord2,i,j,k,dr)
             if (dsqr  <=  rmxsqr) then
               ipr=ipr+1
               iwk(1:5,ipr) = (/INT(100000*dsqr),ibas,i,j,k/)
             endif
           enddo
         enddo
       enddo
     enddo

     call ishell(5,ipr,iwk)

     scalpha%iax(1:5,jbas,1:ipr) = iwk(1:5,1:ipr)
     scalpha%iax(8,jbas,1)=0
     scalpha%iax(9,1   ,1)=0
     if (jbas /= 1)                                                   &
       scalpha%iax(9,jbas,1)=scalpha%iax(9,jbas-1,scalpha%npr(jbas-1))&
                    +scalpha%iax(6,jbas-1,scalpha%npr(jbas-1))        &
                    *scalpha%iax(7,jbas-1,scalpha%npr(jbas-1))
     do i = 1,ipr
       scalpha%iax(1,jbas,i)=jbas
       scalpha%iax(6,jbas,i)=nlsqrj
       inx=scalpha%iax(2,jbas,i)
       scalpha%iax(7,jbas,i)=(geo%atoms(inx)%ptr%lmx+1)**2
       if (i /= 1) then
         im1=i-1
         scalpha%iax(8,jbas,i)=scalpha%iax(8,jbas,im1)+scalpha%iax(7,jbas,im1)
         scalpha%iax(9,jbas,i)=scalpha%iax(9,jbas,im1) &
                       +scalpha%iax(6,jbas,im1)*scalpha%iax(7,jbas,im1)
       endif
     enddo

!---- Add information for Watson orbitals
     ipr=ipr+1
     im1=ipr-1
     scalpha%iax(1,jbas,ipr)=jbas
     scalpha%iax(2,jbas,ipr)=jbas
     scalpha%iax(3,jbas,ipr)=0
     scalpha%iax(4,jbas,ipr)=0
     scalpha%iax(5,jbas,ipr)=0
     scalpha%iax(6,jbas,ipr)=nlsqrj
     scalpha%iax(7,jbas,ipr)=0
     scalpha%iax(8,jbas,ipr)=scalpha%iax(8,jbas,im1)+scalpha%iax(7,jbas,im1)
     scalpha%iax(9,jbas,ipr)=scalpha%iax(9,jbas,im1)+scalpha%iax(6,jbas,im1)*scalpha%iax(7,jbas,im1)
     scalpha%npr(jbas) = ipr

   enddo
   deallocate (iwk)

   return
   end subroutine nghbr0

   subroutine gtdims(geo, neighm, ns)
!- Determines maximum number of neighbors and dimension of s
! ----------------------------------------------------------------------
!o  neighm:maximum number of neighbors around each atom
!o  ns    :dimension of s (considering atoms up to jbas)
!r Remarks:
!r  Attention: gtdims must perform exactly the same calculation 
!r  as subroutime neighm0 does. 
! ----------------------------------------------------------------
   implicit none
   type(t_geometry_EMTO) :: geo

   integer          :: nbas
   real(kind=prec)  :: plat(3,3)
   real(kind=prec)  :: avw 
   integer, intent(out) :: neighm
   integer, intent(out) :: ns
! Local variables:
   logical :: islead
   integer i,i1,i2,i3_l,i3_r,ibas,ipr,j,jbas,k,nlsqri,nlsqrj
   real(kind=prec) :: dsqr,dr(3),rmx,rmxsqr,coord1(3),coord2(3)
   integer,allocatable :: iwk(:,:)

   ns     = 0
   neighm = 0
   rmx=geo%cutrat
   rmxsqr = rmx*rmx

   nbas=geo%num
   avw=geo%dawsr

   if(geo%num==geo%l_transnum)then
     islead=.true.
   else
     islead=.false.
   endif

   i1=geo%ntrpar
   i2=geo%ntrpar
   i3_l=geo%l_ntrperp
   i3_r=geo%r_ntrperp

   allocate (iwk(5,(2*i1+1)*(2*i2+1)*(i3_l+i3_r+1)*nbas))

   plat=0.d0
   plat(1:2,1)=geo%base(1:2,1)/avw
   plat(1:2,2)=geo%base(1:2,2)/avw
   plat(1:3,3)=geo%perp_trans(1:3)/avw

   do jbas = 1,nbas
     nlsqrj=(geo%atoms(jbas)%ptr%lmx+1)**2
     ipr=0
     coord1(1:3)=geo%atoms(jbas)%coord(1:3)/avw
     do  ibas = 1, nbas
       coord2(1:3)=geo%atoms(ibas)%coord(1:3)/avw
!------ Sweep lattice translations to find all neighbors within rmx.
       do i = -i1, i1
         do j = -i2, i2
           do k = -i3_l, i3_r
             if((.not.islead).and.(k.ne.0)) then
               if(k.lt.0)then
                 if(ibas.gt.geo%l_transnum) cycle
                 plat(1:3,3)=geo%l_perp_trans(1:3)/avw
               elseif(k.gt.0)then
                 if(ibas.le.geo%num-geo%r_transnum) cycle
                 plat(1:3,3)=geo%r_perp_trans(1:3)/avw
               else
                 plat(1:3,3)=0.d0
               endif
             endif
             dsqr = drr2(plat,coord1,coord2,i,j,k,dr)
             if (dsqr  <=  rmxsqr) then
               ipr=ipr+1
               iwk(1:5,ipr) = (/INT(100000*dsqr),ibas,i,j,k/)
             endif
           enddo
         enddo
       enddo
     enddo

     call ishell(5,ipr,iwk)

     do i=1, ipr
       ibas=iwk(2,i)
       nlsqri=(geo%atoms(ibas)%ptr%lmx+1)**2
       ns = ns + nlsqrj * nlsqri
     enddo
! ----- Add information for Watson orbitals
! ----- although we do not use Watson,
! ----- the dimension of matrix is still kept
! ----- as in STUTTGART code and neighm=neighm+1
     ipr=ipr+1
     neighm = MAX(neighm,ipr)
   enddo

   deallocate (iwk)
   return
   end subroutine gtdims

   subroutine dinv33(matrix,iopt,invers,det)
!- Inverts 3X3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   matrix:input matrix
!i   iopt  :if 0, usual inverse
!i             1, transpose of inverse
!i             2, 2*pi*inverse
!i             3, 2*pi*transpose of inverse
!o Outputs:
!o   invers:as modified according to iopt
!o   det   :determinant, (or det/2*pi if iopt=2,3)
!r Remarks:
!r  To generate reciprocal lattice vectors, call dinv33(plat,3,plat)
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer         ,intent(in)  :: iopt
   real(kind=prec) ,intent(in)  :: matrix(3,3)
   real(kind=prec) ,intent(out) :: invers(3,3)
   real(kind=prec) ,intent(out) :: det
! Local variables:
   real(kind=prec) ,parameter :: twopi=6.28318530717958648d0

   call cross(matrix(1:3,2),matrix(1:3,3),invers(1:3,1))
   call cross(matrix(1:3,3),matrix(1:3,1),invers(1:3,2))
   call cross(matrix(1:3,1),matrix(1:3,2),invers(1:3,3))
   
   det = DOT_PRODUCT(matrix(1:3,1),invers(1:3,1))
   
   if (det  ==  0.d0) then
     write(*,400)matrix
     stop 'dinv33 error'
   endif
    
   if (iopt >= 2) det = det/twopi
   if (MOD(iopt,2)  ==  0) then
     invers=TRANSPOSE(invers)
   endif
   
   invers=invers/det

400 format(' DINV33: vanishing determinant of matrix:',&
           '|',3x,3f12.5,'|',3x,3f12.5,'|',3x,3f12.5,'$')
   return
   end subroutine dinv33

   real(kind=prec)  function d3nrm2(r)
!- 'norm'-square
! ----------------------------------------------------------------------
!i Inputs:
!i   r     :vector
!o Outputs:
!o   d3nrm2:norm-square
!r Remarks:
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   real(kind=prec) ,intent(in) :: r(3)

   d3nrm2=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)

   return
   end function d3nrm2

   subroutine ishell(m,n,iarray)
!- shell sort of a array of integer vectors
! ----------------------------------------------------------------------
!i Inputs:
!i   m     :number of components in iarray
!i   n     :number of elements in iarray
!i   iarray:array to be sorted
!o Outputs:
!o   iarray:array to be sorted
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer,intent(in   ) :: m
   integer,intent(in   ) :: n
   integer,intent(inout) :: iarray(m,0:n-1)
! Local variables:
   integer lognb2,i,j,k,l,n2,nn,it,mm,mmm

   lognb2 = int(log(float(n+1))*1.4426950)
   n2 = n
   do nn = 1, lognb2
     n2 = n2/2
     k = n - n2
     do  11  j = 1, k
       i = j - 1
 3     continue
       l = i + n2
       do  15  mm = 1, m
         if (iarray(mm,l) - iarray(mm,i)) 16,15,11
16       continue
         do mmm = 1, m
           it = iarray(mmm,i)
           iarray(mmm,i) = iarray(mmm,l)
           iarray(mmm,l) = it
         enddo
         i = i - n2
         if (i  >=  0) goto 3
         goto 11
15     continue
11   continue
   enddo

   return
   end subroutine ishell

   subroutine bessl(y,lmin,lmax,fi,gi)
!- Radial part of Bessel functions:
! ----------------------------------------------------------------------
!i Inputs:
!i   y     :y = e*r**2 = (kappa*r)**2 = x**2
!i   lmin  :minimum l
!i   lmax  :maximum l
!o Outputs:
!o   fi    :fi(l)=  j_l(x)/x^l    *(2l-1)!!/2; l=0,1,...,lmax,
!o   gi    :gi(l)= -n_l(x)*x^(l+1)/(2l-1)!!  ; l=0,1,...,lmax, y>0
!o                =i*h_l(x)*x^(l+1)/(2l-1)!!  ; l=0,1,...,lmax, y<0
!r Remarks:
!r   j,n,h=j+i*n are spherical Bessel, Neumann and Hankel functions.
!r   Bessel functions are calculated for lmax and lmax-1 by a expansion
!r   in powers of x^2=y:
!r                (-x^2/2)^k
!r  fi =  Sum_k  --------------  = dum(lmx+1-l)
!r               k! (2l+1+2k)!!
!r
!r   The remaining functions are obtained by the recursion formula:
!r
!r      j_{l-2}(x)=(2l-1)j_{l-1}(x)/x-j{l}(x)    l=lmx-2,lmx-3,..,-lmx-1
!r
!r  <==> dum(k) = (2*lmx+5-2*k)*dum(k-1)-y*dum(k-2)   k=3,4,...,2*lmx+2
!r
!r   and the Neumann function are given through:
!r
!r         n_l(x)=j_{-l-1}*(-1)^{l+1}
!r
!r   Andersen's definition of Bessel and Neumann functions is:
!r   J(kappa^2,r) = fi(y=kappa^2*r^2) * (r/w)^l
!r   K(kappa^2,r) = gi(y=kappa^2*r^2) / (r/w)^(l+1)
!r   fi(y) -> 1/(4l+2)  for y->0
!r   gi(y) -> 1         for y->0
!r   where w is the average Wigner Seitz radius
!r
!r   Warning:  For x > 40 this algorithm is numerically unstable !!!
!r   Warning:  for lmax < -1 program must be checked
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer          ,intent(in ) :: lmin
   integer          ,intent(in ) :: lmax
   real(kind=prec)  ,intent(in ) :: y
   real(kind=prec)  ,intent(out) :: fi(lmin:lmax+1)
   real(kind=prec)  ,intent(out) :: gi(lmin:lmax+1)
! Local variables:
   integer i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,nf,tlp1,ll1,ll2
   real(kind=prec)  :: dt,dt2,exppr,my,srmy,t
   real(kind=prec)  :: dum     (2*MAX(lmax,2)+2)
   real(kind=prec)  :: fac2l(lmin:MAX(lmax,2)+1)
   logical         , parameter :: lhank=.true.

   lmx = MAX(lmax,2)

   ! --- A table of fac2l(l)=(2l-1)!!
   !     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
   fac2l(0) = 1.d0
   do l = 1, lmx+1
      fac2l(l) = fac2l(l-1) * (l+l-1)
   enddo
   do l = -1,lmin,-1
      fac2l(l) = fac2l(l+1) / (l+l+1)
   enddo

   ! --- Case kappa=0
   if (y == 0.d0)then
      do l = lmin, lmx
         fi(l) = 1.d0/(4*l+2)
         gi(l) = 1.d0
      enddo
      return
   endif

   my = -y

   ! --- get dum(1)=j_{lmx}(x)/x^{lmx}=fi(lmx)
   tlp1 = lmx+lmx+1
   dt = 1.d0
   t = 1.d0
   i = 0
   do while ( ABS(dt) >= EPSILON(t) )
      i   = i + 2
      dt2 = i + tlp1
      dt  = dt *my/(i*dt2)
      t   = t + dt
   enddo
   dum(1) = t/fac2l(lmx+1)

   ! --- get dum(2)=j_{lmx-1}(x)/x^{lmx-1}=fi(lmx-1)
   tlp1 =  tlp1-2
   dt = 1.d0
   t = 1.d0
   i = 0
   do while ( ABS(dt) >= EPSILON(t) )
      i   = i + 2
      dt2 = i + tlp1
      dt  = dt*my/(i*dt2)
      t   = t + dt
   enddo
   dum(2) = t/fac2l(lmx)

   ! --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
   ll1 = lmx + lmx + 1
   ll2 = ll1 + 1
   nf = ll1
   do k = 3, ll2
      nf = nf-2
      dum(k) = nf*dum(k-1) - y*dum(k-2)
   enddo

   ! --- Get fi and gi from dum
   lmxp1 = lmx+1
   lmxp2 = lmx+2
   isn = (-1)**lmin
   do k = lmin, lmx
      j1 = lmxp1-k
      j2 = lmxp2+k
      fi(k) = dum(j1)
      ! ----- n_l(x)=j_{-l-1}*(-1)^{l+1}
      gi(k) = dum(j2)*isn
      isn = -isn
   enddo

   ! --- below, a part is provided which defines spherical Hankel functions
   ! --- as basis functions for EKAP < 0.
   if (lhank .and. y < 0d0) then
      srmy = SQRT(-y)
      gi(0) = 1.d0
      gi(1) = 1.d0+srmy
      if (lmx >= 2) then
         tlp1 = 1
         do l = 2, lmx
            tlp1 = tlp1+2
            gi(l) = tlp1*gi(l-1) - y*gi(l-2)
         enddo
      endif
      tlp1 = 3
      if (lmin <= -1) then
         do l = -1,lmin,-1
            tlp1 = tlp1-2
            gi(l)=((l+l+3)*gi(l+1)-gi(l+2))/y
         enddo
      endif
      exppr = 1.d0/EXP(srmy)
      do l = lmin, lmx
         gi(l) = gi(l)*exppr
      enddo
   endif

   do l = lmin, lmx
      fi(l) = fi(l)*fac2l(l)*0.5d0
      gi(l) = gi(l)/fac2l(l)
   enddo

   return
   end subroutine bessl

   real(kind=prec)  function f100(fac,j1,j2,j3,m1,m2,m3)
!- calculating factors for Clebsch-Gordon
   implicit none
! Passed variables:
   integer          ,intent(in) :: j1,j2,j3,m1,m2,m3
   real(kind=prec)  ,intent(in) :: fac(*)
! Local variables:
   integer n,n1,n2,np1,nm1
   real(kind=prec)  t,t1

   f100 = 0.d0
   if (m3 /= m1+m2) return
   n1 = MAX(j2-j3-m1,j1-j3+m2,0) + 1
   n2 = MIN(j1+j2-j3,j1-m1,j2+m2) + 1
   if (n1 > n2) return
   t = (2*j3+1)*fac(j1+j2-j3+1)*fac(j3+j1-j2+1)*fac(j3+j2-j1+1)/fac(j1+j2+j3+2)
   t = SQRT(t*fac(j1+m1+1)*fac(j1-m1+1)*fac(j2+m2+1)*fac(j2-m2+1)*  &
        fac(j3+m3+1)*fac(j3-m3+1))
   t1 = 0.d0
   do np1 = n1, n2
      n = np1 - 1
      nm1 = n - 1
      t1 = t1 + dble(1+4*(n/2)-2*n)/                                  &
           (fac(np1)*fac(j1+j2-j3-nm1)*fac(j1-m1-nm1)*                    &
           fac(j2+m2-nm1)*fac(j3-j2+m1+np1)*fac(j3-j1-m2+np1))
   enddo
   f100 = t*t1

   return
   end function f100

   real(kind=prec)  function f102(fac,l1,l2,l3)
!- calculating factors for Clebsch-Gordon
   implicit none
! Passed variables:
   integer          ,intent(in) :: l1,l2,l3
   real(kind=prec)  ,intent(in) :: fac(*)
! Local variables:
   integer :: lt,p,x

   lt = l1 + l2 + l3
   p = lt/2
   if (p+p  /=  lt) then
      f102 = 0.d0
   else
      f102 = SQRT(DBLE(l3+l3+1)/DBLE(lt+1))
      f102 = f102*fac(p+1)/SQRT(fac(p+p+1))
      x = p-l1
      f102 = f102*SQRT(fac(x+x+1))/fac(x+1)
      x = p-l2
      f102 = f102*SQRT(fac(x+x+1))/fac(x+1)
      x = p-l3
      f102 = f102*SQRT(fac(x+x+1))/fac(x+1)
      if (x  >  2*(x/2)) f102 = -f102
   endif

   return
   end function f102

   subroutine sylmnc(cy,lmax)
!- Normalization constants for the spherical harmonics
! ----------------------------------------------------------------------
!i Inputs:
!i   lmax  :maximum l
!o Outputs:
!o   cy    :normalization constants
!o   for given (l,m!=0) sqrt((2l+1)/2pi)*sqrt((l-|m|)!/(l+|m|)!) 
!o   for given (l,0   ) sqrt((2l+1)/4pi)
!r Remarks:
!r   use together with SYLM
!r   The true spherical harmonic is:  Y_L = yl(lm,r) / r^l
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer         ,intent(in ) :: lmax
   real(kind=prec) ,intent(out) :: cy((lmax+1)**2)
! Local variables:
   integer :: i,i1,i2,l,lav,lp1,m,n1,n2,n3
   real(kind=prec)  :: fn2,tlp1,y0
   real(kind=prec) , parameter :: sqrtpi = 1.7724538509055160d0
   real(kind=prec) , parameter :: twopi  = 6.28318530717958648d0
   real(kind=prec) , parameter :: fpi    =12.5663706143591725d0

   y0 = 0.5d0/sqrtpi
   cy(1) = y0
   do l = 1, lmax
      lp1 = l+1
      tlp1 = l+lp1
      lav = l*lp1 + 1
      cy(lav) = SQRT(tlp1/fpi)
      do m = 1, l
         n2 = lp1-m
         n1 = n2+1
         n3 = l+m
         fn2 = n2
         do i = n1, n3
            fn2 = fn2*i
         enddo
         i1 = lav+m
         i2 = lav-m
         cy(i1) = SQRT(tlp1/(fn2*twopi))
         cy(i2) = cy(i1)
      enddo
   enddo

   return
   end subroutine sylmnc

   subroutine sylm(cy,dr,ylm,lmax,r)
!- Generate normalized spherical harmonic polynomials
! ----------------------------------------------------------------------
!i Inputs:
!i   cy    :normalization constants
!i   dr    :vector with 3 components
!i   lmax  :maximum l for which ylm is calculated
!o Outputs:
!o   ylm   :normalized spherical harmonic polynomials
!o   r     :length of dr
!r Remarks:
!r   polar axis along 001 axis. (adapted from asw programs)
!r   use together with sylmnc:
!r   The true spherical harmonic is:  Y_L = ylm(lm,r) / r^l
!r   The first 9 values are for ylm/cy:
!r   l=0:            1
!r   l=1:       y    z    x
!r   l=2:  xy  yz  3zz-rr  xz  xx-yy
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer         , intent(in ) :: lmax
   real(kind=prec) , intent(in ) :: cy( (lmax+1)**2)
   real(kind=prec) , intent(in ) :: dr(3)
   real(kind=prec) , intent(out) :: ylm((lmax+1)**2)
   real(kind=prec) , intent(out) :: r
! Local variables:
   real(kind=prec) , parameter :: r2min=1.d-16
   integer :: i,l,lav,lavml,lavmm,lavpl,lavpm,lm1,lmm,lp1,m,mp1,nlsqr,nt
   real(kind=prec)  :: x,y,z,r2,zz,xxyy
   real(kind=prec)  :: c(lmax+1),s(lmax+1),p(lmax+1,lmax+1)

   nlsqr = (lmax+1)*(lmax+1)
   ylm(1) = cy(1)
   x = dr(1)
   y = dr(2)
   z = dr(3)
   zz=z*z
   xxyy=x*x+y*y
   r2 = xxyy+zz
   r = SQRT(r2)

   if (lmax  >  0) then
      if (r2  <=  r2min) then
         do i = 2, nlsqr
            ylm(i) = 0.d0
         enddo
         return
      endif
      ylm(2) = y*cy(2)
      ylm(3) = z*cy(3)
      ylm(4) = x*cy(4)
   endif

   if (lmax  >  1) then
      c(1)=1.d0
      c(2)=x
      s(1)=0.d0
      s(2)=y
      p(1,1)=1.d0
      p(2,1)=z
      p(2,2)=1.d0
      nt = 1
      do l = 2, lmax
         lp1 = l+1
         lm1 = l-1
         lav = l*lp1 + 1
         p(lp1,1) = ((l+lm1)*z*p(l,1) - lm1*r2*p(lm1,1)) / l
         ylm(lav) = p(lp1,1)*cy(lav)
         nt = nt+2
         p(lp1,lp1) = p(l,l)*nt
         c(lp1) = x*c(l) - y*s(l)
         s(lp1) = x*s(l) + y*c(l)
         lavpl = lav+l
         ylm(lavpl) = p(lp1,lp1)*c(lp1)*cy(lavpl)
         lavml = lav-l
         ylm(lavml) = p(lp1,lp1)*s(lp1)*cy(lavml)
         ! ----- In order to avoid numerical instabilities distinguish
         if (xxyy  >  zz) then
            do m = 1, lm1
               mp1 = m+1
               lavpm = lav+m
               lavmm = lav-m
               p(lp1,mp1) = ((lm1+m)*r2*p(l,m)-(lp1-m)*z*p(lp1,m))/(xxyy)
               ylm(lavpm) = p(lp1,mp1)*c(mp1)*cy(lavpm)
               ylm(lavmm) = p(lp1,mp1)*s(mp1)*cy(lavmm)
            enddo
         else
            do lmm = 1, lm1
               m = l-lmm
               lavpm = lav+m
               lavmm = lav-m
               mp1 = m+1
               p(lp1,mp1) = (r2*(l+m)*p(l,mp1)-xxyy*p(lp1,mp1+1))/(z*lmm)
               ylm(lavpm) = p(lp1,mp1)*c(mp1)*cy(lavpm)
               ylm(lavmm) = p(lp1,mp1)*s(mp1)*cy(lavmm)
            enddo
         endif
      enddo
   endif

   return
   end subroutine sylm

   subroutine scg(lmax,cg,indxcg,jcg)
!- Computes Clebsch-Gordan coefficients
! ----------------------------------------------------------------------
!i Inputs:
!i   lmax  :maximum l
!o Outputs:
!o   cg    :Clebsch Gordon coefficients
!o   indxcg:gives index for Clebsch Gordon coefficients
!o   jcg   :lm" see remarks
!r Remarks:
!r   (FORMERLY S104 IN ASW)
!r   Notation: L=(l,m),  L'=(l',m'),  L"=(l",m")
!r             lm:=1+m+l(1+l); lm':=1+m'+l'(1+l') ;lm":=1+m"+l"(1+l")
!r   for fixed  lm and lm' with lm'<=lm
!r   the first nonvanishing C_{LL'L"} is : cg(indxcg(lm*(lm-1)/2+lm'))
!r   and the corresponding lm" is jcg(indxcg(lm*(lm-1)/2+lm'))
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer         , intent(in ) :: lmax
   integer         , intent(out) :: indxcg(*)
   integer         , intent(out) :: jcg(*)
   real(kind=prec) , intent(out) :: cg (*)
! Local variables:
   integer :: i,i1,i2,i3,i31,i32,ic,j1,j1s,j2,j2s,k2,l1,l2,l3,lmindx,   &
              m1,m2,m3,mb,n1,n2,n3,nl,nm3,s1,s2,s3,t1,t2,t3
   real(kind=prec)  :: q1,t,fs
   real(kind=prec)  :: fac(4*lmax+2)
   real(kind=prec) ,parameter :: sqrtpi=1.7724538509055160d0
   real(kind=prec) ,parameter :: sr2   =1.41421356237309515d0

   fs(i) = 1 + 4*(i/2) - i - i
   nl = lmax+1

   fac(1) = 1.d0
   do i = 1, 4*nl-3
      fac(i+1) = i*fac(i)
   enddo

   ic = 0
   lmindx = 0
   do i1 = 1, nl
      l1 = i1-1
      j1s = l1+l1+1
      do j1 = 1, j1s
         m1 = j1-i1
         n1 = ABS(m1)
         s1 = 0
         if (m1  <  0) s1 = 1
         t1 = 0
         if (m1  ==  0) t1 = 1
         do i2 = 1, i1
            l2 = i2-1
            i31 = l1 - l2 + 1
            i32 = l1 + l2 + 1
            j2s = l2+l2 + 1
            k2 = j1s*j2s
            if (i2  ==  i1) j2s = j1
            do j2 = 1, j2s
               lmindx = lmindx + 1
               indxcg(lmindx) = ic + 1
               m2 = j2-i2
               n2 = ABS(m2)
               s2 = 0
               if (m2  <  0) s2 = 1
               t2 = 0
               if (m2  ==  0) t2 = 1
               if (m1*m2 < 0) then
                  m3 = -n1 - n2
                  mb = -ABS(n1-n2)
                  if (mb  ==  0) then
                     nm3 = 1
                  else
                     nm3 = 2
                  endif
               elseif (m1*m2 == 0) then
                  m3 = m1+m2
                  nm3 = 1
               else
                  m3 = n1+n2
                  mb = ABS(n1-n2)
                  nm3 = 2
               endif
               do while (nm3 > 0)
                  n3 = ABS(m3)
                  s3 = 0
                  if (m3  <  0) s3 = 1
                  t3 = 0
                  if (m3  ==  0) t3 = 1
                  q1=SQRT(dble(k2))*fs(n3+(s1+s2+s3)/2)/ (2*sr2**(1+t1+t2+t3))
                  do i3 = i31, i32,2
                     l3 = i3-1
                     if (n3 <= l3) then
                        t = 0.d0
                        if (n1+n2  == -n3) t = t + f102(fac,l1,l2,l3)
                        if (n1+n2  ==  n3) t = t + f100(fac,l1,l2,l3, n1,n2,n3)*fs(n3+s3)
                        if (n1-n2  == -n3) t = t + f100(fac,l1,l2,l3,n1,-n2,-n3)*fs(n2+s2)
                        if (n1-n2  ==  n3) t = t + f100(fac,l1,l2,l3,-n1,n2,-n3)*fs(n1+s1)
                        ic = ic + 1
                        cg(ic) = q1*t*f102(fac,l1,l2,l3)/(sqrtpi*SQRT(dble(l3+l3+1)))
                        jcg(ic) = l3*(l3+1) + m3 + 1
                     endif
                  enddo
                  nm3 = nm3-1
                  m3 = mb
               enddo
            enddo
         enddo
      enddo
   enddo
   indxcg(lmindx+1) = ic + 1

   return
   end subroutine scg

   subroutine mkalph(kap2,nlmax,avw,atoms)
!!$--- Generates screening parameters alpha_{Rl}
   implicit none
! Passed variables:
   integer           :: nlmax
   real(kind=prec)   :: kap2,avw
   type(t_atoms_set_EMTO) :: atoms 
! Local variables
   integer         :: l,ic,nclass,lmx
   !real(kind=prec) :: fi(-1:nlmax),gi(-1:nlmax),hcr
   real(kind=prec), Allocatable :: fi(:),gi(:)
   real(kind=prec) :: hcr

   Allocate(fi(-1:max(nlmax,2)), gi(-1:max(nlmax,2)))

   nclass=atoms%num
   do ic = 1, nclass
      lmx=atoms%at(ic)%lmx
      allocate(atoms%at(ic)%alpha(0:lmx))
      allocate(atoms%at(ic)%alphadot(0:lmx))
      do  l = 0, lmx
         hcr=atoms%at(ic)%a(l) /avw
         call bessl(kap2*hcr*hcr,-1,max(l+1,2),fi(-1),gi(-1))
         atoms%at(ic)%alpha(l)=fi(l)/gi(l)*hcr**(l+l+1)
         atoms%at(ic)%alphadot(l)= -0.5d0*hcr**(l+l+3)*                 &
                 (fi(l+1)*gi(l)/(l+l+1)+fi(l)*gi(l-1)/(l+l-1))/(gi(l)*gi(l))
      enddo
   enddo

   Deallocate(fi, gi)

   return
   end subroutine mkalph

   subroutine mktral(kap2,nlmax,avw,atoms)
!- Makes the transformation matrix
! ----------------------------------------------------------------------
!r Remarks: The transformation matrix is:
!r          |N^a>=t(1)*|N>+t(2)*|J>
!r          |J^a>=t(3)*|N>+t(4)*|J>
!r      and ||N^a>=|N^a>-|J^a>S^a
!r      |N^a(kappa)> has value 1 and derivative 0
!r      |J^a(kappa)> has value 0 and derivative 1/a*a
!r                              at hard core radius a=sigma*wsr
!r      alpha = - t(3) / t(4)
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   type(t_atoms_set_EMTO) :: atoms 
   real(kind=prec)   :: kap2
   integer           :: nlmax
   real(kind=prec)   :: avw
! Local variables:
   integer          :: ic,l,nclass,lmx
   real(kind=prec)  :: al,ad,dj,djdot,dn,dndot,j,jdot,  &
                       n,ndot,r2,rfac,hcr,w1,w2,w3,w4,w5,t(8),fac
   real(kind=prec), allocatable  :: fi(:), gi(:)

   Allocate(fi(-1:max(nlmax,2)), gi(-1:max(nlmax,2)))


   nclass=atoms%num
   do ic=1,nclass
      lmx=atoms%at(ic)%lmx
      allocate(atoms%at(ic)%t(8,0:lmx))
      do l=0,lmx
         hcr=atoms%at(ic)%a(l)  /avw
         r2=hcr*hcr
         !call bessl(kap2*r2,-1,max(l+1,2),fi(-1),gi(-1))
         call bessl(kap2*r2,-1,max(l+1,1),fi(-1),gi(-1))
         rfac=hcr**(-l-1)

         !!$ Bessel functions
         n    = rfac*gi(l)
         ndot = rfac*gi(l-1)/(4*l-2)*r2
         dn   = l-(l+l+1)*gi(l+1)/gi(l)
         dndot= l-(l+l-1)*gi(l  )/gi(l-1) !!$ energy derivative
         j    = fi(l)/(rfac*hcr)
         jdot = -fi(l+1)/(4*l+2)/(rfac/hcr) !!$ ed
         dj   = -(l+1)+(l+l-1)*fi(l-1)/fi(l)
         djdot= -(l+1)+(l+l+1)*fi(l  )/fi(l+1) !!$ ed

         !!$ shorthands
         al=atoms%at(ic)%alpha(l)
         ad=atoms%at(ic)%alphadot(l)

         w1=n*j*hcr*(dj-dn)   *avw
         w2=  j*hcr*(dj   )  *avw
         w3=n  *hcr*(  -dn)  *avw
         w4=jdot*hcr*avw*( djdot) !!$ needed for ed
         w5=ndot*hcr*avw*(-dndot) !!$ idem

         !!$ actually calculate t parameters
         t(1)=w2/w1
         t(2)=w3/w1
         t(4)=1.d0/(t(1)+t(2)*al)
         t(3)=-al*t(4)

         !!$ and those for the energy derivative
         t(5)=w4/w1
         t(6)=w5/w1
         t(7)=-(t(1)*ad/al/al-t(5)/al-t(6))*t(3)*t(3)
         t(8)=-(t(5)+al*t(6)+ad*t(2))*t(4)*t(4)

         !fac=1         ! leads to itrans==3 in NMTO, i.e. |J^a(kappa)> has derivative w/(2*a*a) 
         fac=2.d0/avw ! leads to itrans==4 in NMTO, i.e. |J^a(kappa)> has derivative 1/a*a
         !fac=-2.d0*hcr ! leads to itrans==5 in NMTO, i.e. |J^a(kappa)> has derivative -1/a

         t(3)=t(3)*fac
         t(4)=t(4)*fac
         t(7)=t(7)*fac
         t(8)=t(8)*fac

         !!$ save the parameters
         atoms%at(ic)%t(:,l)=t
      end do
   end do

   deallocate(fi,gi)

   return
   end subroutine mktral

   subroutine mstrx2(kap2,dr,nlsqri,nlsqrj,lmax,cg,indxcg,jcg,cy,s,iop)
!- Calculate one-center expansion coefficents to j of h (struc. consts)
!  and their energy derivative
! ----------------------------------------------------------------------
!i Inputs:
!i   kap2  :(kappa*average WS-radius)**2
!i   dr    :vector connecting sphere at R' to that at R (units of avw)
!i   nlsqri:Generate coefficients S_R'L',RL for L' < nlsqri
!i   nlsqrj:Generate coefficients S_R'L',RL for L  < nlsqrj
!i   nlsqr :row dimension of structure constant matrix s
!i   cg    :Clebsch Gordon coefficients
!i   indxcg:index for cg
!i   jcg   :jcg(icg)=lm"
!i   cy    :normalization constants for spherical harmonics (see SYLMNC)
!o Outputs: (depending on iop)
!o   s     :structure constant matrix
!o          multiplied with a factor -1
!r Remarks:
!r  Conventions for structure constants: (in the following equations
!r  the lengths and energies are given in atomic units, whereas
!r  in the program r is scaled by w (=avw) and kap2 means {kw}^2.
!r  The solutions of the Helmholtz equation are defined as (OKA)
!r  (k=kappa):
!r      H_L(k^2,r)  = - (kw)^(l+1)/(2l-1)!!   n_l(kr) Y_L(r)
!r      J_L(k^2,r)  =   (kw)^(-l) *(2l-1)!!/2 j_l(kr) Y_L(r)
!r  they depend on kappa only through kappa^2 and apply both for
!r  positive and negative kappa^2.
!r  n_0(kr) = - cos(kr)/kr, j_0(kr) = sin(kr)/kr, ... are the
!r  usual Bessel and Neumann Functions of fractional order.
!r
!r  These satisfy:
!r    (1)  cases kap2  /=  0 and kap2  ==  0 same equations
!r    (2)  functions are real
!r    (3)  S_{R'L',RL}  =  S^*_{RL,R'L'}
!r
!r  Expansion Theorem: H_RL(r) = H_L(r-R)
!r
!r  (1)    H_{RL}(k^2,r) = - J_{R'L'}(k^2,r) * S_{R'L',RL}
!r  where
!r  (2)  S_R'L',RL  = - Sum_l" 4 pi C_{LL'L"} i^{-l+l'-l"}  *
!r
!r                           l+l'-l"
!r                     2 (kw)       (2l"-1)!!
!r                  * -----------------------  H_L"(k^2,R-R')
!r                      (2l'-1)!!  (2l-1)!!
!r
!r  Indices for L: middle is m=0; (z axis for m=0)
!r          after is real part of spherical harmonics: m=1,2, ...
!r          before is imaginary part m=1,2, ... (* +/- sqrt(2))
!r  For s,p,d the indices 1-9 correspond to:
!r         1     2     3     4     5     6       7     8     9
!r         s     y     z     x     xy    yz  3z^2-r^2  xz  x^2-y^2
!r
!r   the structure constants differ by the ratio
!r   S_RL,R'L'(MSM) / S_RL,R'L'(OKA) = -(2l-1)!!(2l'-1)!!/2
!r
!r  Note: the kap2-> zero limit of the ss(sigma) sdot diverges as
!r  1/sqrt(kap2).  This singularity is removed for case kap2=0, making
!r  a discontinuity in Sdot at kap2=0.
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer          ,intent(in ) :: lmax
   integer          ,intent(in ) :: nlsqri
   integer          ,intent(in ) :: nlsqrj
   integer          ,intent(in ) :: indxcg(*)
   integer          ,intent(in ) :: jcg(*)
   real(kind=prec)  ,intent(in ) :: kap2
   real(kind=prec)  ,intent(out) :: s   (nlsqri,nlsqrj)
   real(kind=prec)  ,intent(in ) :: cy(*)
   real(kind=prec)  ,intent(in ) :: cg(*)
   real(kind=prec)  ,intent(in ) :: dr(3)
   integer          ,intent(in ) :: iop
! Local variables:
   integer icg,icg1,icg2,ii,ilm,indx,ipow,jlm,jlm0,klm,l,li,lj,lk,mk
   real(kind=prec)  :: edot (0:lmax)
   real(kind=prec)  :: efac (0:lmax)
   real(kind=prec)  :: fac2l(0:lmax)
   real(kind=prec)  :: sig  (0:lmax)
   real(kind=prec), allocatable  :: fi(:), gi(:)
   real(kind=prec)  :: hl  ((lmax+1)**2)
   real(kind=prec)  :: hdot((lmax+1)**2)
   real(kind=prec)  :: ylm ((lmax+1)**2) 
   real(kind=prec)  :: gdot,r,rfac,sum,sumd
   real(kind=prec) , parameter :: fpi   =12.5663706143591725d0 !!$ MAGIC NUMBER, BOO! it's 4*pi 
   real(kind=prec) , parameter :: rmin  =1.d-10

   allocate(fi(-1:max(lmax,2)), gi(-1:max(lmax,2)))

   ! --- fac2l(l) = (2l-1)!!
   ! --- efac(l)=(-kap2)^l
   ! --- sig(l)=(-)^l
   efac(0) = 1.d0
   edot(0) = 0.d0
   sig (0) = 1.d0
   fac2l(0) = 1.d0
   do l = 1, lmax
      efac (l) = -kap2*efac(l-1)
      edot (l) = -l   *efac(l-1)
      sig  (l) = -sig(l-1)
      fac2l(l) = fac2l(l-1) * (2*l-1)
   enddo

   call sylm(cy,dr,ylm,lmax,r)
   if (r  <  rmin) return
   call bessl(kap2*r*r,-1,max(lmax,2),fi(-1),gi(-1))

   ! --- Generate H_L(DR) and HDOT
   rfac = r
   klm  = 0
   do lk = 0, lmax
      gdot= gi(lk-1)*r*r/(4*lk-2)
      ! ----- rfac = r^-{1+2*lk}
      rfac = rfac/r/r
      do mk = -lk, lk
         klm = klm+1
         hl  (klm) = gi(lk)*ylm(klm)*rfac*fac2l(lk)
         hdot(klm) = gdot  *ylm(klm)*rfac*fac2l(lk)
      enddo
   enddo

   if (iop == 1) then
    ! ----- Generate S_MK
    do ilm = 1, nlsqri
       li = ll(ilm)
       jlm0 = ilm
       ! ------- Handles subblock nlsqrj < ilm <= nlsqri, if there is one
       if (ilm  >  nlsqrj) jlm0 = 1
       do jlm = jlm0, nlsqrj
          lj = ll(jlm)
          ii = MAX(ilm,jlm)
          indx = (ii*(ii-1))/2 + MIN(ilm,jlm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          sum = 0.d0
          do icg = icg1, icg2
             klm  = jcg(icg)
             lk   = ll(klm)
             ipow = (li+lj-lk)/2
             sum  = sum + cg(icg)*efac(ipow)*hl(klm)
          enddo
          s(ilm,jlm) = fpi*sum*2.d0/fac2l(lj)/fac2l(li)*sig(lj)
          if (jlm <= nlsqri .and. ilm <= nlsqrj)                      &
               s(jlm,ilm) = s(ilm,jlm)*sig(lj)*sig(li)
       enddo
    enddo
   elseif (iop == 2) then
     ! ----- Generate Sdot_MK
         do ilm = 1, nlsqri
           li = ll(ilm)
           jlm0 = ilm
     ! ------- Handles subblock nlsqrj < ilm <= nlsqri, if there is one
           if (ilm  >  nlsqrj) jlm0 = 1
           do jlm = jlm0, nlsqrj
             lj = ll(jlm)
             ii = MAX(ilm,jlm)
             indx = (ii*(ii-1))/2 + MIN(ilm,jlm)
             icg1 = indxcg(indx)
             icg2 = indxcg(indx+1) - 1
             sumd = 0d0
             do icg = icg1, icg2
               klm  = jcg(icg)
               lk   = ll(klm)
               ipow = (li+lj-lk)/2
               sumd = sumd + cg(icg)*(efac(ipow)*hdot(klm)+edot(ipow)*hl(klm))
             enddo
             s(ilm,jlm) = fpi*sumd*2.d0/fac2l(lj)/fac2l(li)*sig(lj)
             if (jlm  <=  nlsqri .and. ilm  <=  nlsqrj)                  &
               s(jlm,ilm) = s(ilm,jlm)*sig(lj)*sig(li)
           enddo
         enddo
   else
    write(*,*) 'iop set to ', iop, ' which is not allowed. For reference, iop=1 is S and iop=2 is Sdot'
    stop
   endif

   return
   end subroutine mstrx2

   subroutine mstrx3(kap2,dr,nlsqri,nlsqrj,lmax,cg,indxcg,jcg,cy,s,iop)
!- Calculate one-center expansion coefficents to J of J (struc. consts)
! ----------------------------------------------------------------------
!i Inputs:
!i   kap2  :(kappa*average WS-radius)**2
!i   dr    :vector connecting sphere at R' to that at R (units of avw)
!i   nlsqri:Generate coefficients S_R'L',RL for L' < nlsqri
!i   nlsqrj:Generate coefficients S_R'L',RL for L  < nlsqrj
!i   lmax  :
!i   cg    :Clebsch Gordon coefficients
!i   indxcg:index for cg
!i   jcg   :jcg(icg)=lm"
!i   cy    :normalization constants for spherical harmonics (see SYLMNC)
!i   iop   : determines if S or Sdot is calculated. (S -> iop=1, Sdot -> iop=2)
!o Outputs: (depending on iop)
!o   s     :structure constant matrix
!o          multiplied with a factor -1
!o         OR its energy derivative, depending on iop
!r Remarks:
!r  Conventions for structure constants: (in the following equations
!r  the lengths and energies are given in atomic units, whereas
!r  in the program r is scaled by w (=avw) and kap2 means {kw}^2.
!r  The solutions of the Helmholtz equation are defined as (OKA)
!r  (k=kappa):
!r
!r  H_L(k^2,r)  = - (kw)^(l+1)/(2l-1)!!   n_l(kr) Y_L(r)
!r  J_L(k^2,r)  =   (kw)^(-l) *(2l-1)!!/2 j_l(kr) Y_L(r)
!r
!r  they depend on kappa only through kappa^2 and apply both for
!r  positive and negative kappa^2.
!r  n_0(kr) = - cos(kr)/kr, j_0(kr) = sin(kr)/kr, ... are the
!r  usual Bessel and Neumann Functions of fractional order.
!r
!r  These satisfy:
!r    (1)  cases kap2  /=  0 and kap2  ==  0 same equations
!r    (2)  functions are real
!r    (3)  S_{R'L',RL}  =  S^*_{RL,R'L'}
!r
!r  Expansion Theorem: J_RL(r) = J_L(r-R)
!r
!r  (1)    J_{RL}(k^2,r) = - J_{R'L'}(k^2,r) * S_{R'L',RL}
!r  where
!r  (2)  S_R'L',RL  = - Sum_l" 4 pi C_{LL'L"} i^{-l+l'-l"}  *
!r
!r                          -l+l'+l"
!r                     2 (kw)       (2l-1)!!
!r                  * -----------------------  J_L"(k^2,R-R')
!r                      (2l'-1)!!  (2l"-1)!!
!r
!r  Indices for L: middle is m=0; (z axis for m=0)
!r          after is real part of spherical harmonics: m=1,2, ...
!r          before is imaginary part m=1,2, ... (* +/- sqrt(2))
!r  For s,p,d the indices 1-9 correspond to:
!r         1     2     3     4     5     6       7     8     9
!r         s     y     z     x     xy    yz  3z^2-r^2  xz  x^2-y^2
!r   the structure constants differ by the ratio
!r   S_RL,R'L'(MSM) / S_RL,R'L'(OKA) = -(2l'-1)!!/(2l-1)!!
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer          ,intent(in ) :: lmax
   integer          ,intent(in ) :: nlsqri
   integer          ,intent(in ) :: nlsqrj
   integer          ,intent(in ) :: indxcg(*)
   integer          ,intent(in ) :: jcg(*)
   real(kind=prec)  ,intent(in ) :: kap2
   real(kind=prec)  ,intent(in ) :: cy(*)
   real(kind=prec)  ,intent(in ) :: cg(*)
   real(kind=prec)  ,intent(in ) :: dr(3)
   real(kind=prec)  ,intent(out) :: s   (nlsqri,nlsqrj)
   integer          ,intent(in ) :: iop
! Local variables:
   integer          :: icg,icg1,icg2,ii,ilm,indx,ipow,jlm,klm,l,li,lj,lk,mk
   real(kind=prec)  :: fdot,r,sum,sumd
   real(kind=prec)  :: edot (0:lmax)
   real(kind=prec)  :: efac (0:lmax)
   real(kind=prec)  :: fac2l(0:lmax)
   real(kind=prec)  :: sig  (0:lmax)
   real(kind=prec)  :: fi   (0:lmax+1)
   real(kind=prec)  :: gi   (0:lmax+1)
   real(kind=prec)  :: bl  ((lmax+1)**2)
   real(kind=prec)  :: bdot((lmax+1)**2)
   real(kind=prec)  :: ylm ((lmax+1)**2) 
   real(kind=prec) , parameter :: fpi   =12.5663706143591725d0

   ! --- fac2l(l) = (2l-1)!!
   ! --- efac(l)=(-kap2)^l
   ! --- sig(l)=(-)^l
   efac (0) = 1.d0
   edot (0) = 0.d0
   sig  (0) = 1.d0
   fac2l(0) = 1.d0
   do l = 1, lmax
      efac (l) = -kap2*efac(l-1)
      edot (l) = -l   *efac(l-1)
      sig  (l) = -sig(l-1)
      fac2l(l) = fac2l(l-1) * (2*l-1)
   enddo

   call sylm(cy,dr,ylm,lmax,r)
   call bessl(kap2*r*r,0,lmax+1,fi(0),gi(0))
   klm = 0
   do lk = 0, lmax
      fdot  = -fi(lk+1)*r*r/(4*lk+2)
      do mk = -lk, lk
         klm = klm+1
         bl  (klm) = fi(lk)*ylm(klm)/fac2l(lk)
         bdot(klm) = fdot  *ylm(klm)/fac2l(lk)
      enddo
   enddo

   if (iop == 1) then
   ! ----- Generate S_IJ
   do ilm = 1, nlsqri
      li = ll(ilm)
      do jlm = 1, nlsqrj
         lj = ll(jlm)
         ii = MAX(ilm,jlm)
         indx = (ii*(ii-1))/2 + MIN(ilm,jlm)
         icg1 = indxcg(indx)
         icg2 = indxcg(indx+1) - 1
         sum = 0.d0
         do icg = icg1, icg2
            klm  = jcg(icg)
            lk   = ll(klm)
            ipow = (li-lj+lk)/2
            sum  = sum + cg(icg)*efac(ipow)*bl(klm)*sig(lk)
         enddo
         s(ilm,jlm) = fpi*sum*2.d0*fac2l(lj)/fac2l(li)
      enddo
   enddo
   elseif (iop == 2) then
 ! ----- Generate Sdot_IJ
         do ilm = 1, nlsqri
           li = ll(ilm)
           do jlm = 1, nlsqrj
             lj = ll(jlm)
             ii = MAX(ilm,jlm)
             indx = (ii*(ii-1))/2 + MIN(ilm,jlm)
             icg1 = indxcg(indx)
             icg2 = indxcg(indx+1) - 1
             sumd = 0d0
             do icg = icg1, icg2
               klm  = jcg(icg)
               lk   = ll(klm)
               ipow = (li-lj+lk)/2
               sumd = sumd + cg(icg)*sig(lk)*                            &
                             (efac(ipow)*bdot(klm)+edot(ipow)*bl(klm))
             enddo
             s(ilm,jlm) = fpi*sumd*2.d0*fac2l(lj)/fac2l(li)
           enddo
         enddo
   else
     write(*,*) 'iop set to ', iop,' in mstrx3 which is an illegal value. iop=1 is S and iop=2 is Sdot'
     stop
   endif
   return
   end subroutine mstrx3

   subroutine renorm_sJ(nlsqri,nlsqrj,sJ,kap2)
   implicit none
! Passed variables
   integer,          intent(in   ) :: nlsqri,nlsqrj
   real(kind=prec) , intent(in   ) :: kap2
   real(kind=prec) , intent(inout) :: sJ(nlsqri,nlsqrj)
! Local variables
   integer          :: l,m,lm
   real(kind=prec)  :: efac,fac2l 

   lm    = 0
   efac  = sqrt(kap2)
   fac2l = 1.d0
   do l = 0, ll(nlsqrj)
      do m = -l, l
         lm = lm+1
         sJ(:,lm) = 2.d0*efac/fac2l/fac2l*sJ(:,lm)
      end do ! m
      efac  = efac*kap2           ! kap2**(l+0.5)
      fac2l = fac2l*dble(2*l+1)   ! (2l-1)!!
   end do ! l

   return
   end subroutine renorm_sJ

!  real(kind=prec)  function avwsr(plat,vol,nbas)
!  implicit none
!  integer          ,intent(in ) :: nbas
!  real(kind=prec)  ,intent(in ) :: plat(3,3)
!  real(kind=prec)  ,intent(out) :: vol
!  real(kind=prec) ,parameter :: fpi3=4.18879020478639053d0
!  vol = ABS(ddet33(plat))
!  avwsr = (vol/fpi3/nbas)**(1.d0/3.d0)
!  return
!  end function avwsr

      Subroutine invert_real_packed (n,m,ap,c)
!!$    a wrapper for the lapack routines for inverting a real matrix stored in
!packed format.
         Implicit None
         External dgetrf, dgetri
         
         integer :: n,m
         Real (Kind=prec) :: ap (:)
         Real (Kind=prec) :: c (n,m)
!!$ local
         Real (Kind=prec), Allocatable :: a (:,:)
         Real (Kind=prec), Pointer :: work (:)
         Integer :: info
         Integer, Pointer :: ipiv (:)
         Integer :: i,j,k
!!$ unpack
         Allocate(a(n,n))
         do i=1,n
          k=i*(i-1)/2 
          a(1:i,i) = ap(1+k:i+k)
          a(i,1:i-1) = ap(1+k:i-1+k)
         enddo
!!$

         Allocate (ipiv(n), work(n))

         Call dgetrf (n, n, a, n, ipiv, info)! computes lu factorization (lapack
 
         If (info /= 0) Then
            Write (*,*) 'problems encountered while in dgetrf routine, info =', info
            Stop
         End If
         Call dgetri (n, a, n, ipiv, work, n, info)! inverses the a matrix
         If (info /= 0) Then
            Write (*,*) 'problems encountered while in dgetrfi routine, info =', info
            Stop
         End If
         Deallocate (ipiv, work)

      c = a(1:n,1:m)

      End Subroutine invert_real_packed

      Subroutine invert_complex_packed (n,m,ap,c)
!!$    a wrapper for the lapack routines for inverting a real matrix stored in
!packed format.
         Implicit None
         External zgetrf, zgetri

         integer :: n,m
         Complex (Kind=prec) :: ap (:)
         Real  (Kind=prec) :: c (n,m)
!!$ local
         Complex (Kind=prec), Allocatable :: a (:,:)
         Complex (Kind=prec), Pointer :: work (:)
         Integer :: info
         Integer, Pointer :: ipiv (:)
         Integer :: i,j,k
!!$ unpack
         Allocate(a(n,n))
         do i=1,n
          k=i*(i-1)/2
          a(1:i,i) = ap(1+k:i+k)
          a(i,1:i-1) = ap(1+k:i-1+k)
         enddo
!!$

         Allocate (ipiv(n), work(n))

         Call zgetrf (n, n, a, n, ipiv, info)! computes lu factorization (lapack

         If (info /= 0) Then
            Write (*,*) 'problems encountered while in zgetrf routine, info =', info
            Write (*,*) ' in omta_strr.F90'
            Stop
         End If
         Call zgetri (n, a, n, ipiv, work, n, info)! inverses the a matrix
         If (info /= 0) Then
            Write (*,*) 'problems encountered while in zgetrfi routine, info =', info
            Stop
         End If
         Deallocate (ipiv, work)

         c = a(1:n,1:m)
 

      End Subroutine invert_complex_packed



end module omta_strrs

