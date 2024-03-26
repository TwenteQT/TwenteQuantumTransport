#include "math_def.h"

Module omta_kink
 use omta_defs
 use omta_SOC 
 use geometry_module
 use sparselib
 use logging
 implicit none
 
 Type t_mathop_EMTO
    Type (zcsrmat) :: c, l, r
    Integer :: alloc = 0
    Integer :: n, nl, nr, havehops = 0
 End Type t_mathop_EMTO

Contains

   Subroutine init_mathop_EMTO (m)
      Implicit None
      Type (t_mathop_EMTO) :: m
      m%alloc = 0
      m%havehops = 0
      Call init (m%l)
      Call init (m%r)
      Call init (m%c)
   End Subroutine

   Subroutine free_mathop_EMTO (sk)
      Implicit None
      Type (t_mathop_EMTO) :: sk
      If (sk%alloc /= 0) Then
         Call free (sk%c)
         If (sk%havehops /= 0) Then
            Call free (sk%l)
            Call free (sk%r)
            sk%havehops = 0
         End If
         sk%alloc = 0
      End If
   End Subroutine free_mathop_EMTO

   subroutine kinkmat(sys,geo,scalpha,bki,par,need_hop,tr)
!- Set up Hamiltonian and Overlap and diagonalize secular matrix.
   implicit none
! Passed variables:
   integer        ,    intent(in ) :: need_hop
   real(kind=DEF_DBL_PREC),    intent(in ) :: bki(2)
   type(t_mathop_EMTO)                  :: sys
   type(t_geometry_EMTO)                :: geo
   type(t_omta_struct)             :: scalpha
   type(t_omta_logder)             :: par
   integer, optional :: tr(2) !!$ This is for the currents scheme.
! Local variables:
   type(zcsrmat) :: struct,hopping, logder
   integer :: i


   call alloc (logder,scalpha%ldim,scalpha%ldim)

   do i=1,scalpha%ldim
     logder%A(i)  = cmplx(par%pph(1,i),0.d0,Kind=DEF_DBL_PREC)
        
     logder%jc(i) = i
     logder%ir(i) = i
   enddo
   logder%nnz=scalpha%ldim
   logder%ir(scalpha%ldim+1)=scalpha%ldim+1

   If (present(tr)) Then
   call bloch(bki,geo,scalpha,struct,hopping,need_hop,par,tr)
   Else
   call bloch(bki,geo,scalpha,struct,hopping,need_hop,par)
   End if

   call alloc(sys%c,logder%nnz+struct%nnz,scalpha%ldim)
   if (present(tr)) then
       if (tr(1)==0 .and. tr(2)==0) then
                sys%c=spmatadd(struct,logder)
        else
                call spcopy(struct, sys%c)
                
        end if
   else
        sys%c = spmatadd(struct,logder)
   end if
   sys%n=sys%c%ncol
   call free(logder)
   call free(struct)

   if(need_hop==1)then
!    sys%r=spherm(hopping)
     call  spcopy(hopping,sys%r)
     sys%nr=sys%r%ncol 
     sys%havehops=1
     call free(hopping)
   endif
   sys%alloc=1

   end subroutine kinkmat


   subroutine bloch(bk,geo,scalpha,s0,j0,need_hop,par,tr)
!- Bloch transform of real-space matrix
! ----------------------------------------------------------------------
!r   s(k;r1,l1,r2,l2) = sum_T s(r1,l1,T+r2,l2) * exp(-i k . T)
!r   where r1 and r2 are basis vectors and T = t2-t1 is the difference
!r   in primitive lattice translation vectors.
!r   If not folding down, sil and sii may point to dummy addresses.
! ----------------------------------------------------------------------
   implicit none
!! Passed variables:
   real(kind=DEF_DBL_PREC)     :: bk(2)
   type(t_omta_struct) :: scalpha
   type(t_geometry_EMTO)    :: geo
   type(zcsrmat)       :: s0, j0
   integer             :: need_hop
   type(t_omta_logder)   :: par
   integer, optional :: tr(2)
!! Local variables:
   logical :: hop=.false.
   integer :: i,ibas,ilm,offs,ipr,j,jbas,jlm,k,nlsqri,nlsqrj,nbas,nr
   integer :: ns,snnz,snz,nh,hnnz,hnz
   real(kind=DEF_DBL_PREC)  :: tdotk,cost,sint
   real(kind=DEF_DBL_PREC), parameter :: nzero=1.d-20
   complex(kind=DEF_DBL_PREC), parameter :: czero=(0.d0,0.d0)
   type(zcoomat)                  :: sf, hf
   integer :: dum_i, hf_skip, sf_skip, hf_nnz, sf_nnz
   type(zcsrmat) :: sf_csr, logder

   hop=.false.
   nr=geo%nlmax*geo%nlmax
   nbas=geo%num

   call alloc(sf,scalpha%ns,scalpha%ldim)

   if(need_hop==1)then
     call alloc(hf,scalpha%ns,scalpha%ldim)
     hop=.true.
   endif
   sf_nnz = 1
   hf_nnz = 1 

   do jbas = 1,nbas  ! loop over atoms
    do ipr=1,scalpha%npr(jbas)-1  ! loop over neighbours within cut off
      if (present(tr)) then
       if (tr(1) .ne. scalpha%iax(3,jbas,ipr) &
       .or. tr(2) .ne. scalpha%iax(4,jbas,ipr)) then
         !!$write(*,*) 'cycling'
         cycle
       end if
      end if
        tdotk = 0.d0
        do j = 1,2
          do k = 1,2
              tdotk = tdotk - bk(j)*scalpha%base(j,k)*scalpha%iax(2+k,jbas,ipr)
          enddo
        enddo

       ibas = scalpha%iax(2,jbas,ipr)
       nlsqrj=scalpha%iax(6,jbas,ipr)
       nlsqri=scalpha%iax(7,jbas,ipr)
       offs = scalpha%iax(9,jbas,ipr)
       cost = COS(tdotk)
       sint = SIN(tdotk)

       do jlm = 1, nlsqrj
         if(hop)then
           j = geo%atoms(jbas)%idxlead(jlm) 
         else
           j = geo%atoms(jbas)%idxsh(jlm)
         endif
         do ilm = 1,nlsqri
            offs=offs+1

            if (abs(scalpha%s(offs)) > nzero) then

             if(hop)then
              i = geo%atoms(ibas)%idxlead(ilm)
             else
              i = geo%atoms(ibas)%idxsh(ilm)
             endif
             if(scalpha%iax(5,jbas,ipr)==0)then
                sf_skip =0
                if (sf_nnz>1) then
                 do dum_i = 1, sf_nnz-1 
                        if (sf%icol(dum_i) == i .and. sf%irow(dum_i) == j) then
                                sf%a(dum_i) = sf%a(dum_i) + cmplx(scalpha%s(offs)*cost,-scalpha%s(offs)*sint,Kind=DEF_DBL_PREC)
                                sf_skip =1
                        end if
                 end do
                end if
                if (sf_skip==0) then
                 sf%a(sf_nnz) = cmplx(scalpha%s(offs)*cost,-scalpha%s(offs)*sint,Kind=DEF_DBL_PREC)
                 sf%icol(sf_nnz) = i 
                 sf%irow(sf_nnz) = j
                 sf%nnz = sf_nnz 
                 sf_nnz = sf_nnz + 1
                end if

             elseif(scalpha%iax(5,jbas,ipr)==1.and.hop)then
                hf_skip =0
                if (hf_nnz>1) then
                 do dum_i = 1, hf_nnz-1
                        if (hf%icol(dum_i) == i .and. hf%irow(dum_i) == j) then
                                hf%a(dum_i) = hf%a(dum_i) + cmplx(scalpha%s(offs)*cost,-scalpha%s(offs)*sint,Kind=DEF_DBL_PREC)
                                hf_skip =1    
                        end if
                 end do
                end if
                if (hf_skip==0) then
                 hf%a(hf_nnz) = cmplx(scalpha%s(offs)*cost,-scalpha%s(offs)*sint,Kind=DEF_DBL_PREC)
                 hf%icol(hf_nnz) = i 
                 hf%irow(hf_nnz) = j
                 hf%nnz = hf_nnz 
                 hf_nnz = hf_nnz + 1
                end if
             endif
            endif
         enddo
       enddo
     enddo  
   enddo

   s0 = conv2csr(sf)
   call ordercsr(s0)
   !!$call sp_add_dupes(s0) 
 
   if (hop) j0 = conv2csr(hf)
   if (hop) call ordercsr(j0)
   !!$if (hop) call sp_add_dupes(j0)
  
   return    
   end subroutine bloch

   subroutine makidx_lead(geo,out_ldim)
!- Makes row and column indices for struc. const. and ham. matrices
   implicit none
!  Passed variables:
   type(t_geometry_EMTO)    :: geo
!  Local variables:
   integer :: nl2,ldim,out_ldim,idim
   integer :: ibas,idxnew,idxold,inxsh(3),iprior,ip,ilm

   idxnew = 0
   do iprior = 1,3
      idxold = 0
      do ibas = 1, geo%num
         nl2=size(geo%atoms(ibas)%idxdn_m)
         do ilm = 1, nl2
            idxold = idxold+1
            ip=geo%atoms(ibas)%idxdn_m(ilm) !idx_lm,R for m-dependent downfoding
            if (ip <= 0.or.ip >= 4) then
               write(*,403)ibas,ilm,ip,3
               geo%atoms(ibas)%idxdn_m(ilm)=3
               ip=3
            endif
            if(ip == iprior) then
               idxnew = idxnew+1
               geo%atoms(ibas)%idxlead(ilm) = idxnew
            endif
         enddo
      enddo
      inxsh(iprior) = idxnew
   enddo
   if (idxnew  /=  idxold) then
      write(*,401)idxnew,idxold
      stop
   endif
   ldim=inxsh(1)
   idim=inxsh(2)-inxsh(1)
  
   if(ldim == 0) then
      write(*,402)
      stop
   endif

   out_ldim=ldim
  !!$ scalpha%idim=idim THIS IS OLD

401 format( ' MAKIDX: bad input, idxnew=',i3,'idxold=',i3)
402 format( ' MAKIDX: ldim = 0 does not work.')
403 format( ' WARNING in MAKIDX: site',i2,', L=',i1,', bad IDXDN=',i2,', replaced by ',i1)
   return
   end subroutine makidx_lead

   subroutine makidx_m(geo,out_ldim)
!- Makes row and column indices for struc. const. and ham. matrices
   implicit none
!  Passed variables:
   type(t_geometry_EMTO)    :: geo
!  Local variables:
   integer :: nl2,ldim,idim,out_ldim
   integer :: ibas,idxnew,idxold,inxsh(3),iprior,ip,ilm

   idxnew = 0
   do iprior = 1,3
      idxold = 0
      do ibas = 1, geo%num
         nl2=size(geo%atoms(ibas)%idxdn_m)
         do ilm = 1, nl2
            idxold = idxold+1
            ip=geo%atoms(ibas)%idxdn_m(ilm) !idx_lm,R for m-dependent downfoding
            if (ip <= 0.or.ip >= 4) then
               write(*,403)ibas,ilm,ip,3
               geo%atoms(ibas)%idxdn_m(ilm)=3
               ip=3
            endif
            if(ip == iprior) then
               idxnew = idxnew+1
               geo%atoms(ibas)%idxsh(ilm) = idxnew
            endif
         enddo
      enddo
      inxsh(iprior) = idxnew
   enddo
   if (idxnew  /=  idxold) then
      write(*,401)idxnew,idxold
      stop
   endif
   ldim=inxsh(1)
   idim=inxsh(2)-inxsh(1)
  
   if(ldim == 0) then
      write(*,402)
      stop
   endif

   out_ldim=ldim
 !!$  scalpha%idim=idim

401 format( ' MAKIDX: bad input, idxnew=',i3,'idxold=',i3)
402 format( ' MAKIDX: ldim = 0 does not work.')
403 format( ' WARNING in MAKIDX: site',i2,', L=',i1,', bad IDXDN=',i2,&
         ', replaced by ',i1)
   return
   end subroutine makidx_m

   subroutine makpph_nmto(nspin,geo,lidim, par,lead,so,nonsph)
   use sparselib
!- Create a vector of site-dependent potential parameters:
! ----------------------------------------------------------------------
!o   pph   :potential parameters in downfolding order
!r Remarks:
!r   Mapping of composite RL index to this vector is same as that of
!r   eigenvectors permutated according to downfolding rules in idxsh
!r   pph(1) : aD
!r   pph(2) : aDdot
! TODO: clean up this subroutine, it's truly terrible right now
! ----------------------------------------------------------------------
! Passed variables:
   integer             :: nspin,lidim(2)
   type(t_geometry_EMTO)    :: geo
   type(t_omta_logder) :: par(2)
   logical             :: lead
   integer, optional   :: so, nonsph
! Local variables:
   integer         :: nbas,isp,ibas,idx,ilm,l,m
   real(kind=DEF_DBL_PREC) :: aD1,aD2
   real(kind=DEF_DBL_PREC) :: aDdot1,aDdot2
   integer :: i, j, ilm1, ilm2
   type (zcsrmat) :: zcsr_temp, zcsr_temp2, zcsr_nonsph

   nbas=geo%num
   !do isp=1,nspin
   do isp=1,1
      par(isp)%lidim=lidim(isp)*2
      allocate(par(isp)%pph(2,2*lidim(isp)))
      do ibas = 1, nbas
         ilm=1
         do l = 0,geo%atoms(ibas)%ptr%lmx
            aD1   =geo%atoms(ibas)%ptr%pp(1,l,1)
            aD2   =geo%atoms(ibas)%ptr%pp(1,l,2)
            aDdot1=geo%atoms(ibas)%ptr%pp(2,l,1)
            aDdot2=geo%atoms(ibas)%ptr%pp(2,l,2)
            do m = -l, l
               if(lead)then
                 idx=geo%atoms(ibas)%idxlead(ilm) 
               else
                 idx=geo%atoms(ibas)%idxsh(ilm)
               endif
               if (idx <= lidim(isp)) then
                  par(isp)%pph(1,idx * 2 - ilm)=aD1
                  par(isp)%pph(1,idx * 2 - ilm+(geo%atoms(ibas)%ptr%lmx+1)**2)=aD2 !the other spin
                  par(isp)%pph(2,idx * 2 - ilm)=aDdot1
                  par(isp)%pph(2,idx * 2 - ilm+(geo%atoms(ibas)%ptr%lmx+1)**2)=aDdot2 !the other spin
               endif
               ilm=ilm+1
            enddo
         enddo
      enddo
   enddo

!!$ kind of hacky... we put both up and down logder into the first spin. We need
!this later

      if (nonsph == 0) call alloc (par(1)%logder,2*lidim(1),2*lidim(1))
      if (nonsph > 0 ) call alloc (par(1)%logder,10*lidim(1),2*lidim(1))
       
      i=1
      ! do j=1,2
        do idx=1,2*lidim(1)

         par(1)%logder%a(i) = cmplx(par(1)%pph(1,idx),0d0,kind=DEF_DBL_PREC)
         par(1)%logder%jc(i) = i
         par(1)%logder%ir(i) = i
         i = i+1
          
        enddo
      ! enddo
        
      par(1)%logder%ir(i) = i
      par(1)%logder%nnz = i-1

      if (present(so) .and. so>0) then
         call prep_so_EMTO (geo, par(1)%hsoc)
      else
         call alloc(par(1)%hsoc,0,par(1)%logder%nrow) 
      endif      



      !!$ do the same hacky thing for the energy derivative of logder
      call alloc (par(1)%logderdot,2*lidim(1),2*lidim(1))
       
      i=1
      !do ibas=1,nbas
       !do j=1,2
        do idx=1,2*lidim(1)

        par(1)%logderdot%a(i) = cmplx(par(1)%pph(2,idx),0d0) 
        par(1)%logderdot%jc(i) = i
        par(1)%logderdot%ir(i) = i
        i = i+1
         
       enddo
     ! enddo
     !enddo
      par(1)%logderdot%ir(i) = i
      par(1)%logderdot%nnz = i-1


      !do isp=1,nspin
      do isp=1,1
       deallocate(par(isp)%pph) !don't need this anymore
      enddo
   return
   end subroutine makpph_nmto

 End Module omta_kink
