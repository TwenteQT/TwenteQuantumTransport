#include "math_def.h"

 module omta_pots
 use omta_defs
 use omta_strrs
 use logging
 use sparselib
 implicit none

 Contains

   subroutine read_omta_pot(filename,atoms)
   implicit none
   integer :: nclass
   character(len=*) :: filename
   type(t_atoms_set_EMTO) :: atoms
   !!$ local
   integer :: ic, nl, listf=73001
   character(len=100) :: cwork
   Call do_log (1, 'Reading omta atom definitions')
   Open (Unit=listf, File=filename, Action='read')

!!$ Skip first two lines left for comments
   Read (listf,*)
   Read (listf,*)
   Read (listf,*) nclass

   Write (cwork,*) nclass
   Call do_log (1, ' Number of potentials = '//trim(cwork))

   atoms%num=nclass
   Allocate (atoms%at(nclass))

   do ic=1, nclass

     read (listf,*) atoms%at(ic)%label,atoms%at(ic)%at_numb,atoms%at(ic)%lmx, &
        &atoms%at(ic)%magmom, atoms%at(ic)%nonsph_lmx, atoms%at(ic)%so_scale
     nl=atoms%at(ic)%lmx+1
     atoms%at(ic)%nm=nl*nl
     allocate(atoms%at(ic)%idxdn(0:nl-1))
     read (listf,*) atoms%at(ic)%idxdn(0:nl-1)
     allocate(atoms%at(ic)%a(0:nl-1))

     !allocate(atoms%at(ic)%a_au(0:nl-1))

     allocate(atoms%at(ic)%sigma(0:nl-1))
     read (listf,*) atoms%at(ic)%sigma(0:nl-1)
     atoms%at(ic)%label_len=len_trim(atoms%at(ic)%label)
     call aiopot(atoms%at(ic))
   enddo
   read(listf,*) atoms%nspin
   allocate(atoms%vmtz(atoms%nspin))
   read(listf,*) atoms%fermi
   read(listf,*) atoms%vmtz(1:atoms%nspin)
   close(listf)

   write(cwork,'(i2)') atoms%nspin
   call do_log (1, ' Number of spins = '//trim(cwork))
   write(cwork,'(f18.10)') atoms%fermi
   call do_log (1, ' Fermi level at = '//trim(cwork))
   write(cwork,'(2(1x,f18.10))') atoms%vmtz(:)
   call do_log (1, ' Muffin-tin zero at = '//trim(cwork))
   call do_log (1, 'Reading potentials done!')
   call do_log (1, '')

   return
   end subroutine read_omta_pot

   subroutine mkrofi(a,b,nr,rofi)
!- Makes radial mesh rofi(i) = b [e^(a(i-1)) -1]
   implicit none
! Passed variables:
   integer        , intent(in ) :: nr
   real(kind=prec), intent(in ) :: a
   real(kind=prec), intent(in ) :: b
   real(kind=prec), intent(out) :: rofi(nr)
   ! Local variables:
   integer          :: ir
   real(kind=prec)  :: ea,rpb

   ea = EXP(a)
   rpb = b
   do ir = 1, nr
      rofi(ir) = rpb - b
      rpb = rpb*ea
   enddo

   end subroutine mkrofi

   logical function lscat(ifi,categ,lrewnd)   
!- Scan file for a category
   implicit none
! Passed variables:    
   integer       ,intent(in) :: ifi
   logical       ,intent(in) :: lrewnd
   character*(*) ,intent(in) :: categ
! Local variables:
   integer,parameter :: recl0=80
   integer           :: len_cat
   character*(recl0) :: line

   if (lrewnd) rewind ifi

   lscat = .true.
   len_cat=len(categ)

   do 
      read(ifi,'(a)',end=20) line
      if (line(1:len_cat)==categ) return
   enddo

20 lscat = .false.

   return
   end function lscat

   subroutine aiopot(aclass)
!- File I/O for cell potential.
   implicit none
   integer, parameter :: atf=45
! Passed variables:
   type(tclass) :: aclass
! Local variables:
   integer nr,nsp
   real(kind=prec) :: a,wsr
   integer i,iost,isp
   real(kind=prec) b
   real(kind=prec), parameter :: tiny=1.d-5
   character(len=1024) :: filename

   open(unit=atf,file='atoms_EMTO/'//trim(aclass%label),action='read')

   if (.not. lscat(atf,'POT:',.true.)) stop 'lscat error' 
   read(atf,301,err=15,end=15,iostat=iost) nr,nsp,a,wsr
   
   aclass%nsp=nsp   
   aclass%s  =wsr   
   do i=0,aclass%lmx
     aclass%a(i)=aclass%sigma(i)*wsr
  enddo

   allocate(aclass%r(1:nr))
   allocate(aclass%v(1:nr,nsp))

   b = wsr/(EXP(a*(nr - 1)) - 1.d0)
   aclass%ra=a
   aclass%rb=b
   aclass%nr=nr
   call mkrofi(a,b,nr,aclass%r(1:nr))
   do isp = 1, nsp
      !read(atf,300,err=15,end=15,iostat=iost)(aclass%v(i,isp),i = 1,nr)
      read(atf,*,err=15,end=15,iostat=iost)(aclass%v(i,isp),i = 1,nr)
   enddo
   close(atf)
  
!!$ open(191,file='l1m0.dat')
!!$  do i=1,nr
!!$        write(191,*) aclass%r(i)*2*0.3402/sqrt(3/(4*3.1415))/(2*wsr)
!!$  end do
!!$ close(191)
   return

! --- handle read exception ---
15 continue
   write(*,401) iost
   stop

300 format(1p,5d17.9)
301 format(2i5,2f12.5)
401 format(' AIOPOT: reading error, iostat=',i4,'$')

   return
   end subroutine aiopot

   subroutine aiopot_nonsph(aclass)
!- File I/O for cell potential.
   implicit none
   integer, parameter :: atf=45
! Passed variables:
   type(tclass) :: aclass
! Local variables:
   integer nr,nsp
   real(kind=prec) :: a,wsr
   integer i, iost, l2, m2, lm2
   real(kind=prec) b
   real(kind=prec), parameter :: tiny=1.d-5
   character(len=1024) :: filename
   character (len=8) :: fmt2, fmt3 !format descriptor
   character (len=3) :: x2,x3
   fmt2='(I2)'
   fmt3='(I3)'

   nr = size(aclass%r)
 
   allocate(aclass%nonsph_v(1:nr,(aclass%nonsph_lmx+1)**2))
  
   do l2 = 1,aclass%nonsph_lmx
    do m2 = -l2,l2

    lm2 = (l2+1)**2 - l2 + m2

    if (abs(l2)<10) then
     write(x2,fmt2) l2
    else
     write(x2,fmt3) l2
    endif

    if (abs(m2)<10) then
     write(x3,fmt2) m2
    else
     write(x3,fmt3) m2  
    endif  
       
    open(unit=atf+lm2,file='atoms_EMTO/'//trim(aclass%label)//'_'//&
&trim(adjustl(trim(x2)))//'_'//trim(adjustl(trim(x3))),action='read')
    if (.not. lscat(atf+lm2,'POT:',.true.)) stop 'lscat error'
     !read(atf+lm2,301,err=15,end=15,iostat=iost) nr,nsp,a,wsr
     read(atf+lm2,*,err=15,end=15,iostat=iost) nr,nsp,a,wsr

     !read(atf+lm2,300,err=15,end=15,iostat=iost)(aclass%nonsph_v(i,lm2),i = 1,nr)
     read(atf+lm2,*,err=15,end=15,iostat=iost)(aclass%nonsph_v(i,lm2),i = 1,nr)
     close(atf+lm2)

    end do
   end do

   return

! --- handle read exception ---
15 continue
   write(*,401) iost
   stop

300 format(1p,5d16.9)
301 format(2i5,2f12.5)
401 format(' AIOPOT: reading error, iostat=',i4,'$')
   end subroutine aiopot_nonsph

   subroutine atompp_nmto(atoms)
!- Makes potential parameters for all atoms
   implicit none
! Passed variables:
   type(t_atoms_set_EMTO) :: atoms
! Local variables:
   integer           :: ic

   do ic = 1,atoms%num
     call potpar50(atoms%nspin,atoms%vmtz,atoms%at(ic),atoms%fermi)
   enddo

   return
   end subroutine atompp_nmto

   subroutine potpar50(nsp,vmtz,aclass,globalEnu)
!- Generates potential parameters: for given potential and energies
   implicit none
! Passed variables:
   type(tclass) :: aclass
   integer, intent(in) :: nsp
   real(kind=prec) ,intent(in) :: vmtz(nsp), globalEnu
! Local variables:
   real(kind=prec)  ::a,b
   integer          :: isp,l,nod,nr,lmx,ir, imin
   real(kind=prec)  :: e,slo,val,wsr,hcr,s,s0,dla,ea,rpb,z
   real(kind=prec) , pointer :: g(:,:,:),gfac(:,:),ph0(:,:)

   Real (Kind=prec), Pointer :: f1(:,:), f2(:,:), dV (:,:), &
           & xi0(:,:), xi1(:,:), xi2(:,:), nl(:,:,:), jl(:,:,:), Bess_nl(:,:), Bess_jl(:,:)
   Real (Kind=prec), Pointer :: deltaV(:)
       Character (Len=200) :: cwork
   real (Kind=prec), Pointer :: radialphi(:,:,:)
   integer :: l1, m1, lm1, l2, m2, lm2, l3, m3, lm3, i

   nr=aclass%nr
   wsr=aclass%s
   a=aclass%ra
   b=aclass%rb
   lmx=aclass%lmx
   z=1.0d0*aclass%at_numb
  
 
   allocate(g(nr,2,nsp),gfac(nr,nsp),ph0(nr,nsp),dV(nr,3),f1(nr,nsp),f2(nr,nsp))
   allocate(nl(0:lmx+1,nr,nsp), jl(0:lmx+1,nr,nsp))
   !allocate(xi0(2,2), xi1(2,2), xi2(2,2), Bess_nl(nr,2), Bess_jl(nr,2))
   allocate(Bess_nl(nr,2), Bess_jl(nr,2))
   allocate(aclass%pp(2,0:lmx,nsp),aclass%xi0(nsp,nsp,0:lmx))!,aclass%xi1(nsp,nsp,0:lmx),aclass%xi2(nsp,nsp,0:lmx))
   allocate(aclass%split(0:lmx),deltaV(nr))
   allocate(aclass%nonsph_phiVphi((lmx+1)**2,(aclass%nonsph_lmx+1)**2,(lmx+1)**2,3))
   allocate(radialphi(nr,0:lmx,2))

   !xi0 = 0d0
   !xi1 = 0d0
   !xi2 = 0d0

   ea  = exp(a)
   rpb = b*ea
   e   = globalEnu

   do isp = 1, nsp
   !!$
!!$ radial derivative of the potential
       call dlangr (5, nr-1, aclass%r(2:nr),-(z+z)/aclass%r(2:nr) + aclass%v(2:nr,isp), dV(2:nr, isp) )
       dV(1,isp) = 0
       dV(1:nr, isp) = dV(1:nr, isp) * aclass%r(1:nr)
!!$ ----------------------------------
!call do_log(1,'dlangr successful')

   end do
   
    deltaV(1:nr) = aclass%v(1:nr,2) - aclass%v(1:nr,1)

!!$ For the spin-offdiagonal components, we need the average of dV(:,1) and dV(:,2)
     dV(1:nr, 3) = 0.5d0 * (dV(1:nr, 1) + dV(1:nr, 2))
!!$--------------------------------------------------------------------------------

!call do_log(1,'attempting to generate the bessel functions')
   do isp = 1, nsp
     do ir = 1, nr
       call bessl((e-vmtz(isp))*aclass%r(ir)**2, 0, max(lmx,2), jl(:,ir,isp), nl(:,ir,isp))
     end do
   end do
!call do_log(1,'successful')

   do l = 0, aclass%lmx
     do isp = 1, nsp
       hcr=wsr*aclass%sigma(l)
       !--- find phi by integrating SE outwards to S sphere
       call rsq1(a,b,e,g(:,:,isp),l,nr,nod,nr,aclass%r,slo,aclass%v(:,isp),val,z)
       call mkgfac(e,gfac(:,isp),l,nr,aclass%r,aclass%v(:,isp),z)
       !--- normalize phi to phi0(a)=1
       call norm2g(g(:,:,isp),e-vmtz(isp),l,lmx,nr,aclass%sigma,slo,val,wsr)
       !--- make D (log.deriv. of phi0 at hard sphere)
       call makdla(dla,e-vmtz(isp),l,lmx,aclass%sigma,slo,val,wsr)
       !--- calculate phi0 on a mesh
       call mkphi0(a,b,e-vmtz(isp),l,lmx,nr,ph0(:,isp),aclass%sigma,slo,val,wsr,aclass%r)
       !--- integrate outwards from origin to s-phere (<phi|phi>)
       call gintsr(a,b,g(:,:,isp),g(:,:,isp),gfac(:,isp),nr,aclass%r,s)
       !--- integrate backwards in a-s region (<phi0|phi0>)
       call gints2(a,b,ph0(:,isp),ph0(:,isp),hcr,nr,aclass%r,s0) 
       aclass%pp(1:2,l,isp)=(/hcr*dla,s0-s/)   !aD,a*Ddot

       imin=int(1.d0+dlog(hcr/b+1.d0)/a)-1

       Bess_nl(1:imin-1,isp) = 0
       Bess_nl(imin:nr,isp) = nl(l,imin:nr,isp) !*(aclass%r(imin:nr)/wsr)**l
       Bess_jl(1:imin-1,isp) = 0
       Bess_jl(imin:nr,isp) = jl(l,imin:nr,isp) !*(aclass%r(imin:nr)/wsr)**(-l-1)

     end do

!!$ Calculating nonspherical components
     if (aclass%nonsph_lmx > 0) then
         radialphi(1:nr,l,:) = 0.0d0
         radialphi(2:nr,l,1) = g(2:nr,1,1) /aclass%r(2:nr) !- ph0(2:nr,1) + &
                               ! &Bess_nl(2:nr,1)*aclass%t(1,l) + Bess_jl(2:nr,1)*aclass%t(2,l) 
         radialphi(2:nr,l,2) =   Bess_nl(2:nr,1)*aclass%t(3,l) + Bess_jl(2:nr,1)*aclass%t(4,l)
     end if

!call do_log(1,'attempting to set full radial functions')
     !--- radial integrals for SOC ----
     if (nsp==2) then
!!$ There are two radial functions to define

       f1(1:nr,1) = 0
       f1(1:nr,2) = 0

       f1(2:nr,1) = g(2:nr,1,1) /aclass%r(2:nr) !- ph0(2:nr,1) !+ Bess_nl(2:nr,1)*aclass%t(1,l) + Bess_jl(2:nr,1)*aclass%t(2,l)
       f1(2:nr,2) = g(2:nr,1,2) /aclass%r(2:nr) !- ph0(2:nr,2) !+ Bess_nl(2:nr,2)*aclass%t(1,l) + Bess_jl(2:nr,2)*aclass%t(2,l)

       !f2(1:nr,1) = 0
       !f2(1:nr,2) = 0

       !f2(2:nr,1) = Bess_nl(2:nr,1)*aclass%t(3,l) + Bess_jl(2:nr,1)*aclass%t(4,l)
       !f2(2:nr,2) = Bess_nl(2:nr,2)*aclass%t(3,l) + Bess_jl(2:nr,2)*aclass%t(4,l)


       call calc_radialints(f1, f1, dV(1:nr,1:3), aclass%r(1:nr), nr, xi0)
       !call calc_radialints(f1, f2, dV(1:nr,1:3), aclass%r(1:nr), nr, xi1)
       !call calc_radialints(f2, f2, dV(1:nr,1:3), aclass%r(1:nr), nr, xi2)

       aclass%split (l) = 0.5d0 * quad3 (nr, aclass%r, ((f1(1:nr, 1)**2)+(f1(1:nr, 2)**2))*(aclass%r(1:nr)**2)*deltaV(1:nr))                  
       aclass%xi0(:,:,l) = aclass%so_scale * xi0
       !aclass%xi1(:,:,l) = aclass%so_scale * xi1
       !aclass%xi2(:,:,l) = aclass%so_scale * xi2
     else
     aclass%split(l) = 0d0

     end if
   enddo

      call aiopot_nonsph(aclass) !!$ read in the nonspherical potentials

     aclass%nonsph_phiVphi = 0.0d0
     if (aclass%nonsph_lmx > 0) then
       do l1 = 0, aclass%lmx
        do m1 = -l1,l1

         lm1 = (l1+1)**2-l1+m1

         do l2 = 1, aclass%nonsph_lmx
          do m2 = -l2, l2

           lm2 = (l2+1)**2-l2+m2

           do l3 = 0, aclass%lmx
            do m3 = -l3, l3

             lm3 = (l3+1)**2-l3+m3

             aclass%nonsph_phiVphi(lm1,lm2,lm3,1) = quad3(nr, aclass%r, radialphi(1:nr,l1,1) &
              &* 4*aclass%nonsph_V(1:nr,lm2) * radialphi(1:nr,l3,1))

             aclass%nonsph_phiVphi(lm1,lm2,lm3,2) = quad3(nr, aclass%r, radialphi(1:nr,l1,1) &
              &* 4*aclass%nonsph_V(1:nr,lm2) * radialphi(1:nr,l3,2))

             aclass%nonsph_phiVphi(lm1,lm2,lm3,3) = quad3(nr, aclass%r, radialphi(1:nr,l1,2) &
              &* 4*aclass%nonsph_V(1:nr,lm2) * radialphi(1:nr,l3,2))

            end do
           end do

          end do
         end do

        end do
       end do

      end if
   
   deallocate(aclass%nonsph_V) !!$ and deallocate the nonspherical potentials.
   deallocate(radialphi)


!call do_log(1,'attempting to deallocate')
   deallocate(g,gfac,ph0,dV,f1,f2,nl,jl,xi0)
   Deallocate(deltaV)
   deallocate(Bess_nl, Bess_jl)
   return
   end subroutine potpar50

   Subroutine calc_omta_nonsph(atoms)
   use structure
   Implicit None
   type(t_atoms_set_EMTO) :: atoms
   Integer :: lmx
   !!$ local
   Integer :: maxnum, ic, i, j
   Integer :: lm1, lm2, lm3
   Real (Kind=DEF_DBL_PREC), Pointer :: gfrh(:,:,:)
   Complex (Kind=DEF_DBL_PREC), Allocatable :: nonsph_mat(:,:)
   Type (zcoomat) :: zcoo_tmp
   Integer :: ind

   lmx = 0
   do ic=1,atoms%num
    lmx = max(lmx, atoms%at(ic)%lmx)
    lmx = max(lmx, ceiling(atoms%at(ic)%nonsph_lmx/2d0))
   enddo

   gfrh => gaunty(lmx+1)
   maxnum = size(gfrh)

   !!$ generate the matrix of nonsph matrix elements for each atom type
  allocate(nonsph_mat(100,100)) 

  do ic = 1,atoms%num
   do ind = 1,3 !fi1.fi1, fi1.fi2, fi2.fi2 ( see Rang MSc thesis (5.1.27) )
    nonsph_mat = DEF_cmplx_zero 

    do lm1=1,(atoms%at(ic)%lmx+1)**2
     do lm2=1,(atoms%at(ic)%nonsph_lmx+1)**2
      do lm3=1,(atoms%at(ic)%lmx+1)**2
       nonsph_mat(lm1,lm3) = nonsph_mat(lm1,lm3) - &
       & gfrh(lm1,lm2,lm3)*atoms%at(ic)%nonsph_phiVphi(lm1,lm2,lm3,ind) 
      enddo
     enddo
    enddo

    if (atoms%at(ic)%nonsph_lmx > 0) then
     call zdense2coo(nonsph_mat(1:(atoms%at(ic)%lmx+1)**2,1:(atoms%at(ic)%lmx+1)**2), zcoo_tmp, 1d-15) 
     atoms%at(ic)%zcsr_nonsph(ind) = conv2csr(zcoo_tmp)
    else
     call alloc(atoms%at(ic)%zcsr_nonsph(ind), 1, (atoms%at(ic)%lmx+1)**2) ! alloc empty zcsr matrix
    endif

   enddo
  enddo

  nullify(gfrh)

  deallocate(nonsph_mat)
   

   return
   End Subroutine calc_omta_nonsph



Subroutine calc_radialints(f, g, dV, r, nr, h)
 Implicit none
 Real (Kind=prec), Intent(in) :: r(:), dV(:,:)
 Integer, Intent(in) :: nr
 Real (Kind=prec), Intent(in), Pointer :: f(:,:), g(:,:)
 Real (Kind=prec), Intent(Out), Pointer :: h(:,:)
 !!$local vars
 Integer :: is1, is2
 Real (Kind=prec), Parameter :: scf0 = 2.663967179924343d-5 !!$ this is 2/c^2, c=speed of light = 274

 allocate(h(2,2))

 do is1 = 1,2
   do is2 = 1,2
     if (is1==1 .and. is2==1) then
       h(is1, is2) = scf0 * quad3(nr, r, dV(1:nr, 1) * f(1:nr, is1) * g(1:nr, is2))
     else if (is1==2 .and. is2==2) then
       h(is1, is2) = scf0 * quad3(nr, r, dV(1:nr, 2) * f(1:nr, is1) * g(1:nr, is2))
     else
       h(is1, is2) = scf0 * quad3(nr, r, dV(1:nr, 3) * f(1:nr, is1) * g(1:nr, is2))
     end if
   end do
 end do

 return
 End Subroutine calc_radialints



   subroutine rsq1(a,b,e,g,l,nc,nod,nr,rofi,slonc,v,valnc,z)
!- integrate the scalar relativistic eqn outward from 0 to rofi(nc)
! ----------------------------------------------------------------------
!i Inputs:
!i   nc    :integration from 0.0 to rofi(nc)
!o Outputs:
!o   g     :wave function times r
!o   nod   :number of nodes
!o   slonc :slope of g at rofi(nc)
!o   valnc :value of g at rofi(nc)
!r Remarks:
!r   Boundary condition does not fix value of wave function near the
!r   origin; integration can be scaled by an arbritrary factor.
!r   Scalar relativistic version
!r
!r   g'_l =  2M  * f_l + 1/r * g_l
!r   f'_l = -1/r * f_l +  y  * g_l   where y=l(l+1)/(2Mr^2)+V-Eq

!r   see Koelling and Harmon J.Phys.C 10,3107(1977) for details 
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer          ,intent(in ) :: l
   integer          ,intent(in ) :: nc
   integer          ,intent(in ) :: nr
   real(kind=prec)  ,intent(in ) :: a
   real(kind=prec)  ,intent(in ) :: b
   real(kind=prec)  ,intent(in ) :: e
   real(kind=prec)  ,intent(in ) :: rofi(nr)
   real(kind=prec)  ,intent(in ) :: v(nr)
   real(kind=prec)  ,intent(in ) :: z
   integer          ,intent(out) :: nod
   real(kind=prec)  ,intent(out) :: g(nr,2)
   real(kind=prec)  ,intent(out) :: slonc
   real(kind=prec)  ,intent(out) :: valnc
! Local variables:
   integer :: ir,irm1
   real(kind=prec)  :: aa,b1,b2,c,d(2,3),det,df1,df2,df3,dg1,dg2,dg3,   &
        drdi,f0,fllp1,g0,h83,phi,r,r1,r2,r3,r83sq,s,sf,  &
        u,x,x2,y,zz,z1

   c=274.071979d0

   z1=z
   nod   = 0
   zz    = z1 + z1
   fllp1 = l*(l + 1)
   r83sq = 64.d0/9.d0
   r1    =  1.d0/9.d0
   r2    = -5.d0/9.d0
   r3    = 19.d0/9.d0
   h83   =  8.d0/3.d0
   ! --- Approximate g,f by leading term near zero
   g0   = 1.d0
   if (z1 < 0.9d0) then
      s  = l+1
      sf = l
      f0 = l/c
   else
      aa = zz/c
      s  = SQRT(fllp1 + 1.d0 - aa*aa)
      sf = s
      f0 = g0*(s - 1.d0)/aa
   endif
   g(1,1:2)= 0.d0
   do ir = 2, 4
      r     = rofi(ir)
      drdi  = a*(r + b)
      g(ir,1) = g0*r**s
      g(ir,2) = f0*r**sf
      d(1,ir-1) = drdi*g(ir,1)*s/r
      d(2,ir-1) = drdi*g(ir,2)*sf/r
   enddo
   ! --- integrate over rest of points ------
   dg1 = d(1,1)
   dg2 = d(1,2)
   dg3 = d(1,3)
   df1 = d(2,1)
   df2 = d(2,2)
   df3 = d(2,3)
   do ir = 5, nc
      irm1  = ir - 1
      r     = rofi(ir)
      drdi  = a*(r + b)
      phi   = (e + zz/r - v(ir))*drdi/c
      u     = drdi*c + phi
      x     = -drdi/r
      x2    = x*x
      y     = -fllp1*x2/u + phi
      det   = r83sq-x2+u*y
      b1    = g(irm1,1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
      b2    = g(irm1,2)*h83 + r1*df1 + r2*df2 + r3*df3
      g(ir,1) = (b1*(h83-x) + b2*u)/det
      g(ir,2) = (b2*(h83+x) - b1*y)/det
      if    (g(ir,1)*g(irm1,1) < 0.d0) nod=nod+1
      dg1   = dg2
      dg2   = dg3
      dg3   = u*g(ir,2) - x*g(ir,1)
      df1   = df2
      df2   = df3
      df3   = x*g(ir,2) - y*g(ir,1)
   enddo
   valnc = g(nc,1)
   drdi  = a*(rofi(nc) + b)
   slonc = dg3/drdi

   return
   end subroutine rsq1

   subroutine mkgfac(e,gfac,l,nr,rofi,v,z)
!-Makes normalization factor for first component of g
   implicit none
! Passed variables:
   integer          ,intent(in ) :: nr
   integer          ,intent(in ) :: l
   real(kind=prec)  ,intent(in ) :: e
   real(kind=prec)  ,intent(in ) :: v(nr)
   real(kind=prec)  ,intent(in ) :: rofi(nr)
   real(kind=prec)  ,intent(in ) :: z
   real(kind=prec)  ,intent(out) :: gfac(nr)
! Local variables:
   integer          :: ir
   real(kind=prec)  :: c,r,tmc,tmcr,zz,fllp1

   c=274.071979d0

   fllp1 = l*(l+1)
   zz=z+z

   do  ir = 2, nr
      r = rofi(ir)
      tmc = c + (e + zz/r - v(ir))/c
      tmcr= tmc*r
      gfac(ir)= 1.d0+fllp1/tmcr/tmcr
   enddo

! extrapolate to r=0
   gfac(1) = (gfac(2)*rofi(3) - gfac(3)*rofi(2))/(rofi(3) - rofi(2))

   return
   end subroutine mkgfac

   subroutine norm2g(g,kappa2,l,lmx,nr,sigma,slo,val,wsr)
!- Normalizes wave function phi so that phi^0(a)=1
!- The normalization constant is determined by matching
!- partial wave and back-extrapolated partial wave at r=s
!- with both values and slopes, and the back-extrapolated
!- partial wave equals to 1 at r=a (hard sphere) so that
!- it match the screened spherical wave automatically at
!- the hard sphere.
   implicit none
! Passed variables:
   integer        , intent(in)   :: l,lmx,nr
   real(kind=prec), intent(in)   :: sigma(0:lmx),kappa2,wsr
   real(kind=prec), intent(inout):: g(2*nr),slo,val
! Local variables:
   real(kind=prec) :: er2,fi(0:lmx+2),gi(0:lmx+2),wn,wj,dphi,dn,dj,hcr,&
                      phi,phio

   er2 = kappa2*wsr*wsr
   call bessl(er2,0,l+1,fi(0),gi(0))
   phi   = val/wsr
   dphi  = wsr*slo/val-1.d0
   dn    = l-gi(l+1)/gi(l)*(l+l+1)
   dj    = l-fi(l+1)/fi(l)/(l+l+1)*er2
   wj    = (dphi-dj)*fi(l)*phi
   wn    = (dphi-dn)*gi(l)*phi
   hcr   = wsr*sigma(l)
   er2   = kappa2*hcr*hcr
   call bessl(er2,0,l+1,fi(0),gi(0))
   phio  =2.d0*(wn*fi(l)*sigma(l)**l-wj*gi(l)*sigma(l)**(-l-1))
!---normilize phi
   g(:)=g(:)/phio
   val=val/phio
   slo=slo/phio

   return
   end subroutine norm2g

   subroutine makdla(dla,kappa2,l,lmx,sigma,slo,val,wsr)
!- Gets logarithmic derivative of wave function at screening radius r=s
! ----------------------------------------------------------------------
!    dla   :a * logarithmic derivative at a
!           a = hcr = hard-core radius
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer        , intent(in ) :: l,lmx
   real(kind=prec), intent(in ) :: sigma(0:lmx),slo,val,kappa2,wsr
   real(kind=prec), intent(out) :: dla
! Local variables:
   real(kind=prec) :: er2,fi(0:lmx+2),gi(0:lmx+2),wn,wj,dphi,dj,dn,hcr,&
                      phi,rdfi,rdgi

   er2 = kappa2*wsr*wsr
   call bessl(er2,0,l+1,fi(0),gi(0))
   phi   =  val/wsr
   dphi  =  wsr*slo / val - 1.d0
   dn    = l-gi(l+1)/gi(l)*(l+l+1)
   dj    = l-fi(l+1)/fi(l)/(l+l+1)*er2

   !free dphi=dj
   wj    = (dphi-dj)*fi(l)*phi
   wn    = (dphi-dn)*gi(l)*phi
   hcr   = wsr*sigma(l)
   er2   = kappa2*hcr*hcr
   call bessl(er2,0,l+1,fi(0),gi(0))
   rdgi  = l*gi(l)-gi(l+1)*(l+l+1)
   rdfi  = l*fi(l)-fi(l+1)/(l+l+1)*er2

   dla=2.d0*(wn* rdfi*sigma(l)**l-wj* rdgi*sigma(l)**(-l-1))

   return
   end subroutine makdla

   subroutine mkphi0(a,b,kappa2,l,lmx,nr,phi0,sigma,slo,val,wsr,rofi)
!- calculates phi^0 on a radial mesh 
!-----------------------------------------------------------------------
!r- phi is solution of Schroedinger-equation for energy E_nu and 
!r- potential v(r). 
!r- phi^0 is then integrated back for flat potential vmtz and energy
!r- ebak=E_nu.
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer,intent(in):: l,lmx,nr
   real(kind=prec) ,intent(in )::a,b,kappa2,sigma(0:lmx),slo,val,wsr,rofi(nr)
   real(kind=prec) ,intent(out)::phi0(nr)
! Local variables:
   integer i,imin
   real(kind=prec)  er2,fi(0:lmx+2),gi(0:lmx+2),wn,wj,dphi,dn,dj,hcr,phi,ri

   er2 = kappa2*wsr*wsr
   call bessl(er2,0,l+1,fi(0),gi(0))
   phi   = val/wsr
   dphi  = wsr*slo/val-1.d0
   dn    = l-gi(l+1)/gi(l)*(l+l+1)
   dj    = l-fi(l+1)/fi(l)/(l+l+1)*er2
   wj    = (dphi-dj)*fi(l)*phi
   wn    = (dphi-dn)*gi(l)*phi
   hcr   = wsr*sigma(l)
   imin=int(1.d0+dlog(hcr/b+1.d0)/a)-1
   phi0(1:imin-1)=0.d0
   do i=imin,nr
      ri=rofi(i)
      er2=kappa2*ri*ri
      call bessl(er2,0,l+1,fi(0),gi(0))
      phi0(i)=(2.d0*(wn*fi(l)*(ri/wsr)**l-wj*gi(l)*(ri/wsr)**(-l-1)))
   enddo
   return
   end subroutine mkphi0

   subroutine gintsr(a,b,ga,gb,gfac,nr,rofi,sum)
!- Integrate inner product of two wave equations
! ----------------------------------------------------------------------
!r   Uses Simpson's rule
!r   Attention: an integration from 0 to nre cannot be done by putting
!r   nr=nre
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer          ,intent(in ) :: nr
   real(kind=prec)  ,intent(in ) :: a
   real(kind=prec)  ,intent(in ) :: b
   real(kind=prec)  ,intent(in ) :: ga(nr,2)
   real(kind=prec)  ,intent(in ) :: gb(nr,2)
   real(kind=prec)  ,intent(in ) :: gfac(nr)
   real(kind=prec)  ,intent(in ) :: rofi(nr)
   real(kind=prec)  ,intent(out) :: sum
! Local variables:
   real(kind=prec)  :: fi
   integer          :: ir
   real(kind=prec)  :: wgt

   sum = 0.d0
   do  ir = 2, nr-1, 2
      fi  = gfac(ir)*ga(ir,1)*gb(ir,1) + ga(ir,2)*gb(ir,2)
      wgt = rofi(ir) + b
      sum = sum + wgt * fi
   enddo
   sum = sum + sum

   do  ir = 3, nr-2, 2
      fi  = gfac(ir)*ga(ir,1)*gb(ir,1) + ga(ir,2)*gb(ir,2)
      wgt = rofi(ir) + b
      sum = sum + wgt * fi
   enddo

   fi  = gfac(1 )*ga(1 ,1)*gb(1 ,1) + ga(1 ,2)*gb(1 ,2)
   wgt = (rofi(1 ) + b)*0.5d0
   sum = sum + wgt * fi

   fi  = gfac(nr)*ga(nr,1)*gb(nr,1) + ga(nr,2)*gb(nr,2)
   wgt = (rofi(nr) + b)*0.5d0
   sum = sum + wgt * fi

   sum = sum * a * 2.d0/3.d0

   end subroutine gintsr

   real(kind=prec)  function di3int(ix,ya,x)
!- Interpolates y = f(x) for given xa and ya=f(xa)
! ----------------------------------------------------------------------
!i Inputs:
!i   ix    :xa(1) is integer and equals ix
!i          and xa(1),xa(2) and xa(3) differ exactly by 1.d0
!i   ya    :value of f at xa
!i   x     :x-value at which f is interpolated
!o Outputs:
!o   di3int:interpolated value
!r Remarks:
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer ix
   real(kind=prec)  x,ya(3)
! Local variables:
   real(kind=prec)  xa(3)

   xa(1)=dble(ix)
   xa(2)=dble(ix+1)
   xa(3)=dble(ix+2)

   di3int=0.5d0*(x-xa(2))*((x-xa(3))*ya(1)+(x-xa(1))*ya(3))          &
        -(x-xa(1))* (x-xa(3))*ya(2)

   end function di3int

   subroutine gints2(a,b,ga,gb,hcr,nr,rofi,sum0)
!- Integrate inner product of two back extrapolated wave equations
! ----------------------------------------------------------------------
!i Inputs:
!i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
!i   b     :                 -//-
!i   ga    :First wave function (times r)
!i   gb    :Second wave function (times r)
!i   nr    :number of mesh points
!i   rofi  :radial mesh points
!o Outputs:
!o   sum0  :inner product
!r Remarks:
!r   Uses Simpson's rule
!r   Attention: an integration from 0 to nre cannot be done by putting
!r   nr=nre
! ----------------------------------------------------------------------
   implicit none
! Passed variables:
   integer,intent(in)::nr
   real(kind=prec) ,intent(in )::a,b,ga(nr),gb(nr),hcr,rofi(nr)
   real(kind=prec) ,intent(out)::sum0
! Local variables:
   real(kind=prec)  xx,gaa(2),gba(2),sum1,sum2
   integer na

   xx=1.d0+dlog(hcr/b+1.d0)/a
   na=int(xx)
   ! -- must have an odd number of points for integration
   if(mod((nr-na),2)/=0) na=na+1
   na=min0(na,nr-2)
   sum1=SUM(ga(1+na:nr-1:2)*gb(1+na:nr-1:2)*&
        (rofi(1+na:nr-1:2)+b)*rofi(1+na:nr-1:2)*rofi(1+na:nr-1:2))
   sum2=SUM(ga(2+na:nr-2:2)*gb(2+na:nr-2:2)*&
        (rofi(2+na:nr-2:2)+b)*rofi(2+na:nr-2:2)*rofi(2+na:nr-2:2))
   sum0=sum1+sum1+sum2+&
        ga(na)*gb(na)*(rofi(na)+b)*0.5d0*rofi(na)*rofi(na)+&
        ga(nr)*gb(nr)*(rofi(nr)+b)*0.5d0*rofi(nr)*rofi(nr)
   sum0 = sum0*a*2.d0/3.d0
   gaa(1)=di3int(na,ga(na),xx)
   gba(1)=di3int(na,gb(na),xx)
   sum0=sum0+&
        ((gaa(1)*gba(1))*hcr*hcr+(ga(na)*gb(na))*rofi(na)*rofi(na))*&
        (rofi(na)-hcr)/2.d0

   end subroutine gints2

!      Subroutine dlangr (ord, n, x, y, yn)
!         Implicit None
!         Integer :: ord, n
!         Real (Kind=prec) :: x (:), y (:), yn (:)
!!$   Local vars
!         Integer :: n1, n2, i
!         n1 = ord / 2
!         n2 = n - (ord-n1)
!         Do i = 1, n1
!            yn (i) = auxdlangr (ord, 0, x, y, x(i))
!         End Do
!         Do i = n1 + 1, n2
!            yn (i) = auxdlangr (ord, i-n1-1, x, y, x(i))
!         End Do
!         Do i = n2 + 1, n
!            yn (i) = auxdlangr (ord, n2-n1-1, x, y, x(i))
!         End Do

!         Return
!      Contains

!         Function auxdlangr (order, offs, x, y, xn) Result (yd)
!            Implicit None
!            Integer :: order, offs
!            Real (Kind=prec) :: x (:), y (:), xn, yd
!!$   Local vars        
!            Real (Kind=prec) :: t1, t2, t3
!            Integer :: i, j, k
!            yd = 0.0d0
!            Do i = 1 + offs, order + offs
!               t1 = 1
!               Do j = 1 + offs, order + offs
!                  If (j /= i) t1 = t1 * (x(i)-x(j))
!               End Do
!               t1 = y (i) / t1
!               t2 = 0
!               Do k = 1 + offs, order + offs
!                  If (k /= i) Then
!                     t3 = 1
!                     Do j = 1 + offs, order + offs
!                        If ((j /= k) .And. (j /= i)) t3 = t3 * (xn-x(j))
!                     End Do
!                     t2 = t2 + t3
!                  End If
!               End Do
!               yd = yd + t1 * t2
!            End Do
!         End Function auxdlangr
!      End Subroutine dlangr

!!$*******************
!!$XXX    QUAD3   ****
!!$*******************
!      Function quad3 (n, x, y)
!         Implicit None
!         Integer :: n
!         Real (Kind=prec), intent(in) :: x (:), y (:)
!         Real (Kind=prec) :: quad3
!!!$ Local
!        Integer :: i, i0
!         Real (Kind=prec) :: rc1, rcz, rch, rc5, rc12
!         Real (Kind=prec) :: cm, cp, sq15
!         Real (Kind=prec) :: x1, x2, x3, x4, x5, x6, xt1, yt1, xt2, yt2
!!!$---------------------------------------------------
!!!$     QUADRATURE OF INTEGRAND Y(I) OF VARIABLE X(I)
!!!$     (WHERE I=1,2,... N, AND N.GE.6)
!!!$     USING 5TH DEGREE INTERPOLATION POLYNOMIAL
!!!$---------------------------------------------------
!
!         Data rcz / 0.0D0 /, rch / 0.5D0 /, rc1 / 1.0D0 /, rc5 / 5.0D0 /, rc12 / 12.0D0 /
!
!         sq15 = Sqrt (rc1/rc5)
!         cp = rch * (rc1+sq15)
!         cm = rch * (rc1-sq15)
!!
!         quad3 = rcz
!         Do i = 2, n
!            i0 = Max (i-4, 0)
!            i0 = Min (i0, n-6)
!            x1 = x (i0+1)
!            x2 = x (i0+2)
!            x3 = x (i0+3)
!            x4 = x (i0+4)
!            x5 = x (i0+5)
!            x6 = x (i0+6)
!!$                             LAGRANGE INTERPOLATION
!            xt1 = cp * x (i-1) + cm * x (i)
!            yt1 = y (i0+1) * auxf (x1, x2, x3, x4, x5, x6, xt1) + y (i0+2) * auxf (x2, x3, x4, x5, x6, x1, &
!           & xt1) + y (i0+3) * auxf (x3, x4, x5, x6, x1, x2, xt1) + y (i0+4) * auxf (x4, x5, x6, x1, x2, x3, &
!           & xt1) + y (i0+5) * auxf (x5, x6, x1, x2, x3, x4, xt1) + y (i0+6) * auxf (x6, x1, x2, x3, x4, x5, &
!           & xt1)
!!$
!            xt2 = cm * x (i-1) + cp * x (i)
!            yt2 = y (i0+1) * auxf (x1, x2, x3, x4, x5, x6, xt2) + y (i0+2) * auxf (x2, x3, x4, x5, x6, x1, &
!           & xt2) + y (i0+3) * auxf (x3, x4, x5, x6, x1, x2, xt2) + y (i0+4) * auxf (x4, x5, x6, x1, x2, x3, &
!           & xt2) + y (i0+5) * auxf (x5, x6, x1, x2, x3, x4, xt2) + y (i0+6) * auxf (x6, x1, x2, x3, x4, x5, &
!           & xt2)
!!$                                       LOBATTO RULE
!            quad3 = quad3 + (y(i-1)+rc5*(yt1+yt2)+y(i)) * (x(i)-x(i-1)) / rc12
!         End Do

!!$        RETURN
!      Contains
!         Function auxf (t0, t1, t2, t3, t4, t5, tt)
!            Implicit None
!            Real (Kind=prec) :: auxf, tt, t1, t2, t3, t0, t4, t5
!            auxf = (tt-t1) * (tt-t2) * (tt-t3) * (tt-t4) * (tt-t5) / &
!           & ((t0-t1)*(t0-t2)*(t0-t3)*(t0-t4)*(t0-t5))
!         End Function auxf
     ! End Function quad3


 end module omta_pots
