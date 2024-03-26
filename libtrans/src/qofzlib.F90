#include "math_def.h"
Module qofzlib
Use quantlib
   Implicit None

   Interface paverage
      Module procedure paverage_mat, paverage_spaqu, paverage_spadirqu, paverage_rdensemat
   End Interface
   
   Interface getweights
      Module procedure getweights_mat, getweights_spaqu, getweights_spadirqu
   End Interface
   
   Interface putin1mat
      Module procedure putin1mat_spaqu, putin1mat_spadirqu
   End Interface
   
   Interface writeplot
      Module procedure writeplot_mat, writeplot_rdense, writeplot_spadirqu, writeplot_spaqu
   End Interface
    
contains
   
   Subroutine print_help_qofz()
      Implicit None
      write(*,'(x,A)') 'DESCRIPTION:'
      write(*,'(xxx,A)') ' Firstly, this program evaluates the following quantities using transout.h5 :'
      write(*,'(xxx,A)') '  - spin dependent density of states'
      write(*,'(xxx,A)') '  - spin dependent chemical potential'
      write(*,'(xxx,A)') '  - screened spin density'
      write(*,'(xxx,A)') '  - screened spin torque density'
      write(*,'(xxx,A)') ' The last 3 require the use of a screening method.'
      write(*,'(xxx,A)') ' By default this is perfect screening.'
      write(*,'(xxx,A)') ' No other method has been implemented.'
      write(*,*)
      write(*,'(xxx,A)') ' Secondly, all above data is averaged over xy-planes.'
      write(*,'(xxx,A)') ' If the required data is present, interatomic currents are converted to continuum currents'
      write(*,'(xxx,A)') ' and the currents calculated by the ''old'' method are plane-averaged and scaled properly.'
      write(*,*)
      write(*,'(xxx,A)') ' Output can be found in the group "/qofz" inside transout.h5 .'
      write(*,'(xxx,A)') ' In the folder "plots" ascii files are written to be used with grace.'
      write(*,'(xxx,A)') ' The first column in every file is the z coordinate (col 1).'
      write(*,'(xxx,A)') ' Some abbreviations used in the file names of these plots:'
      write(*,'(xxx,A)') ' nopav: no planar average'
      write(*,'(xxx,A)') ' L:     transport from left to right'
      write(*,'(xxx,A)') ' R:     transport from right to left'
      write(*,'(xxx,A)') ' pscr:  perfectly screened'
      write(*,'(xxx,A)') ' torq:  xyz torquances (col 234), units: dimensionless (Hartree/Hartree)'
      write(*,'(xxx,A)') ' dens:  xyz spin and total densities (col 2345), units: number/Hartree'
      write(*,'(xxx,A)') ' ldos:  local (sphere integrated) spin up and down density of states, units: number/Hartree'
      write(*,'(xxx,A)') ' mu:    spin up and down normalized chemical potentials, dimensionless'
      write(*,'(xxx,A)') ' icur:  current in i-direction evaluated from interatomic currents'
      write(*,'(xxx,A)') '        xyz spin and total currents (col 2345)'
      write(*,'(xxx,A)') '        units: dimensionless (fraction of the total conductance of the system)'
      write(*,'(xxx,A)') ' zold:  current in z-direction evaluated from ''old'' method'
      write(*,'(xxx,A)') '        xyz spin and total currents (col 2345)'
      write(*,'(xxx,A)') '        units: dimensionless (fraction of the total conductance of the system)'
      write(*,*)
      write(*,'(x,A)') 'OPTIONS:'
      write(*,'(xxx,A)') '-help'
      write(*,'(xxx,A)') '    Display help info.'
      write(*,'(xxx,A)') '-log'
      write(*,'(xxx,A)') '    Write progress to the terminal.'
      write(*,'(xxx,A)') '-pav [value]'
      write(*,'(xxx,A)') '    Specify an other-than-default maximum z-distance between atoms within an ''xy-plane''.'
      write(*,'(xxx,A)') '    The default value is half of the average wigner-seitz radius.'
      write(*,'(xxx,A)') '    If [value] is 0 , no planar averaging will be carried out.'
      write(*,'(xxx,A)') '    [value] should be in bohr radii.'
      write(*,'(xxx,A)') '-nocfav'
      write(*,'(xxx,A)') '    Do not do averaging over configurations.'
      write(*,'(xxx,A)') '-doleads'
      write(*,'(xxx,A)') '    Calculate spin-dos in the leads.'
      write(*,'(xxx,A)') '-plotall'
      write(*,'(xxx,A)') '    Plot all data of every configuration.'
      write(*,'(xxx,A)') '-oldunit'
      write(*,'(xxx,A)') '    Use the old way to plot spin current.'
      write(*,'(xxx,A)') '-fix'
      write(*,'(xxx,A)') '    Fix the position of planes, in this situation the value of Pav should be Integer '
      write(*,'(xxx,A)') '    and means the number of planes now. ( in test by Lei).'
      write(*,'(xxx,A)') '-dis'
      write(*,'(xxx,A)') '    invoke discrete current scheme'
      write(*,'(xxx,A)') '-dpi'
      write(*,'(xxx,A)') '    integer tuning factor to increase the grid density for discrete current scheme'
      write(*,'(xxx,A)') '-noscr'
      write(*,'(xxx,A)') '    switch off L->R and R->L smoothing'
      write(*,'(xxx,A)') '-cscr'
      write(*,'(xxx,A)') '    print particle accumulations n^{L,R} instead of chemical potentials'
      write(*,'(xxx,A)') '-ofs'
      write(*,'(xxx,A)') '    offset x and y coordinates in the discrete current scheme to plot for actual width of& 
 &the SC, useful when perp tr vector has a component in x or y'

   End Subroutine
   
   Subroutine logscreen(log, adv, msg)
      Implicit None
      Logical :: log
      Integer :: adv
      Character(*) :: msg
      If (log) Then
         If (adv==0) Then
            write(*,'(x,A)', advance='no') msg
         Else 
            write(*,'(x,A)') msg
         End If
      End If
   End Subroutine
   
   Recursive Subroutine get_dirlist(dirlist, path, f, count)
      Implicit None
      Type(t_dirlist), Pointer, Intent(InOut) :: dirlist
      Character(*), Intent(In) :: path
      Integer, Intent(In) :: f
      Integer, Intent(InOut) :: count
      
      Logical :: flag
      Character(40) :: str
      Integer :: i, n, level
      
      !!$ gfortran and ifort treat directories differently. 
      Inquire(DIRECTORY=path, exist=flag) !!$ ifort
      !!$Inquire(file=path//"/.", exist=flag) !!$ gfortran
      If (flag) Then
         Inquire(FILE=path//"/dirlist", exist=flag)
         If (flag) Then
            open(unit=f, file=path//"/dirlist")
               read(f,*) level
               read(f,*) n
               If (level == 0) Then 
                  dirlist%newgeom = n
                  count = count + 1
               End If
               Do i = 1, n
                  read(f,*) str
                  Call get_dirlist(dirlist, path//"/"//trim(str), f+1, count)
               End Do
            close(f)
         Else
            allocate(character(len=len(path)) :: dirlist%path)
            dirlist%path = path
            allocate(dirlist%next)
            dirlist => dirlist%next
            dirlist%newgeom = 0
         End If
      Else
         Write(*,'("Cannot find ",A)') path
         Write(*,'("Skipping,A")') path
      End If
   End Subroutine
   
   Subroutine make_cfdir_name(cfav,icf,cfdir)
      Implicit None
      Character (Len=8), Intent(Out) :: cfdir
      Integer, Intent(In) :: icf
      Logical, Intent(In) :: cfav
      
      Character (Len=5) :: cwork
      
      If (cfav) then
         Write (cwork, '(i5)') icf
         cfdir = 'cf-'//trim(adjustl(cwork))//'/'
      Else
         cfdir = ' '
      End If
   End Subroutine
   
   Subroutine add_confdata(cf1, cf2, factor)
      Implicit None
      Type (confdata), Intent (Inout) :: cf2
      Type (confdata), Intent (Inout) :: cf1
      Double precision, Optional, Intent (In) :: factor
      
      Double precision :: f
      f = 1.0d0
      If (present(factor)) f = factor
      cf1%cond = cf1%cond + cf2%cond*f
      cf1%area = cf2%area
      cf1%area_nm = cf2%area_nm
      cf1%base = cf2%base
      Call add_spaqu(cf1%z, cf2%z, factor)
      Call add_spaqu(cf1%z_si, cf2%z_si, factor)
      Call add_spaqu(cf1%spindos, cf2%spindos, factor)
      Call add_spaqu(cf1%mdir, cf2%mdir, factor)
      Call add_spadirqu(cf1%spinmu, cf2%spinmu, factor)
      Call add_spadirqu(cf1%ts, cf2%ts, factor)
      Call add_spadirqu(cf1%ts_pscr, cf2%ts_pscr, factor)
      If (cf2%do_OAMdens /= 0) Then
        Call add_spadirqu(cf1%OAMdens, cf2%OAMdens, factor)
        Call add_spadirqu(cf1%OAMdens_pscr, cf2%OAMdens_pscr, factor)
      Endif
      If (cf2%do_iatcu /= 0) Then
         Call add_spadirqu(cf1%xcur, cf2%xcur, factor)
         Call add_spadirqu(cf1%ycur, cf2%ycur, factor)
         Call add_spadirqu(cf1%zcur, cf2%zcur, factor)
         Call add_spadirqu(cf1%cons, cf2%cons, factor)
         Call add_spadirqu(cf1%xcur_ps, cf2%xcur_ps, factor)
         Call add_spadirqu(cf1%ycur_ps, cf2%ycur_ps, factor)
         Call add_spadirqu(cf1%zcur_ps, cf2%zcur_ps, factor)
         
         Call add_spaqu(cf1%SHang, cf2%SHang, factor)
      End If
   End Subroutine
   
   Subroutine reset_conf(cf)
      Implicit None
      Type(confdata), Intent(Inout) :: cf
      
      cf%cond = 0.0d0
      Call free_spaqu(cf%z)
      Call free_spaqu(cf%z_si)
      Call free_spaqu(cf%spindos)
      Call free_spaqu(cf%mdir)
      
      Call free_spadirqu(cf%spinmu)
      Call free_spadirqu(cf%ts)
      Call free_spadirqu(cf%ts_pscr)
      Call free_spadirqu(cf%OAMdens)
      Call free_spadirqu(cf%OAMdens_pscr)
      Call free_spadirqu(cf%xcur)
      Call free_spadirqu(cf%ycur)
      Call free_spadirqu(cf%zcur)
      Call free_spadirqu(cf%cons)
      Call free_spadirqu(cf%xcur_ps)
      Call free_spadirqu(cf%ycur_ps)
      Call free_spadirqu(cf%zcur_ps)
      
      Call free_spaqu(cf%SHang)
      cf%cond = 0.0d0
   End Subroutine
   
   Subroutine alloc_spadirqu (squ, mna, lna, rna, ld)
      Implicit None
      Integer, Intent (In) :: mna, lna, rna, ld
      Type (spadirquant), Intent (Inout) :: squ

      If (squ%alloc /= 0 .And. (squ%nm /= mna .Or. squ%nl /= lna .Or. & 
                           &  squ%nr /= rna .Or. squ%ldim /= ld )) Then
          call free_spadirqu(squ)
      End If

      If (squ%alloc == 0) Then
         Allocate (squ%Lm(ld, mna))
         Allocate (squ%Rm(ld, mna))
         Allocate (squ%Lir(ld, lna))
         Allocate (squ%Rt(ld, lna))
         Allocate (squ%Rir(ld, rna))
         Allocate (squ%Lt(ld, rna))
         squ%alloc = 1
         squ%nl = lna
         squ%nr = rna
         squ%nm = mna
         squ%ldim = ld
      End If
   End Subroutine alloc_spadirqu
   
   Subroutine init_spadirqu (squ)
      Implicit None
      Type (spadirquant), Intent (Inout) :: squ
      
      If (squ%alloc /= 0) Then
         squ%Lm = 0.0d0
         squ%Rm = 0.0d0
         squ%Lir = 0.0d0
         squ%Rt = 0.0d0
         squ%Rir = 0.0d0
         squ%Lt = 0.0d0
      Else
         write(*,*) 'Cannot call init_spadirqu(squ) when squ is not allocated.'
         stop
      End If
   End Subroutine
   
   Subroutine add_spadirqu (q1, q2, factor)
      Implicit None
      Type (spadirquant), Intent (Inout) :: q1
      Type (spadirquant), Intent (In) :: q2
      Double precision, Optional, Intent (In) :: factor
      Double precision :: f
      
      If (q1%alloc == 0) Then
         Call alloc_spadirqu(q1, q2%nm, q2%nl, q2%nr, q2%ldim)
         Call init_spadirqu(q1)
      End If
      If (present(factor)) Then
         f = factor
      Else
         f = 1.0d0
      End If
      q1%Lm = q1%Lm + q2%Lm*f
      q1%Rm = q1%Rm + q2%Rm*f
      q1%Lir = q1%Lir + q2%Lir*f
      q1%Rir = q1%Rir + q2%Rir*f
      q1%Lt = q1%Lt + q2%Lt*f
      q1%Rt = q1%Rt + q2%Rt*f
   End Subroutine

   Subroutine free_spadirqu (squ)
      Implicit None
      Type (spadirquant), Intent (Inout) :: squ
          if (squ%alloc==0) return
          Deallocate (squ%Lm, squ%Lir, squ%Lt, squ%Rm, squ%Rir, squ%Rt)
          squ%nm = 0
          squ%nl = 0
          squ%nr = 0
          squ%ldim = 0
          squ%alloc = 0
   End Subroutine free_spadirqu

   Subroutine alloc_spaqu (squ, mna, lna, rna, ld)
      Implicit None
      Integer, Intent (In) :: mna, lna, rna, ld
      Type (spatialquant), Intent (Inout) :: squ

      If (squ%alloc /= 0 .And. (squ%nm /= mna .Or. squ%nl /= lna .Or. & 
                           &  squ%nr /= rna .Or. squ%ldim /= ld )) Then
          call free_spaqu(squ)
      End If

      If (squ%alloc == 0) Then
         Allocate (squ%m(ld, mna))
         Allocate (squ%l(ld, lna))
         Allocate (squ%r(ld, rna))
         squ%alloc = 1
         squ%nl = lna
         squ%nr = rna
         squ%nm = mna
         squ%ldim = ld
      End If
   End Subroutine alloc_spaqu
   
   Subroutine init_spaqu (squ)
      Implicit None
      Type (spatialquant), Intent (Inout) :: squ
      
      If (squ%alloc /= 0) Then
         squ%m = 0.0d0
         squ%l = 0.0d0
         squ%r = 0.0d0
      Else
         write(*,*) 'Cannot call init_spaqu(squ) when squ is not allocated.'
         stop
      End If
   End Subroutine
   
   Subroutine add_spaqu (q1, q2, factor)
      Implicit None
      Type (spatialquant), Intent (Inout) :: q1
      Type (spatialquant), Intent (In) :: q2
      Double precision, Optional, Intent (In) :: factor
      Double precision :: f
      
      If (q1%alloc == 0) Then
         Call alloc_spaqu(q1, q2%nm, q2%nl, q2%nr, q2%ldim)
         Call init_spaqu(q1)
      End If
      If (present(factor)) Then
         f = factor
      Else
         f = 1.0d0
      End If
      q1%m = q1%m + q2%m*f
      q1%l = q1%l + q2%l*f
      q1%r = q1%r + q2%r*f
   End Subroutine

   Subroutine free_spaqu (squ)
      Implicit None
      Type (spatialquant), Intent (Inout) :: squ
          if (squ%alloc==0) return
          Deallocate (squ%m, squ%l, squ%r)
          squ%nm = 0
          squ%nl = 0
          squ%nr = 0
          squ%ldim = 0
          squ%alloc = 0
   End Subroutine free_spaqu
   
   Subroutine remake_oldc(oldc, pl, scaling, cond)
      Use sparselib
      Implicit None
      type(rdensemat) :: oldc
      type(planes_system) :: pl
      Double Precision :: scaling, cond
      
      Integer :: i
      Double precision :: new(2, 1)
      
      Do i = 1, pl%m%np
         new(1, 1) = oldc%bl(1, i) - oldc%bl(2, i)
         new(2, 1) = oldc%bl(1, i) + oldc%bl(2, i)
         oldc%bl(:, i) = new(:, 1)/scaling/pl%th(i)*dble(pl%m%napp(i))/cond
         ! In principle this scaling is not correct. The layer thickness should be determind such that
         ! the volume of every layer corresponds to the sum of volumes of the asa-spheres. Still even
         ! then the old currents would violate current conservation if there are hopping currents that
         ! pass a layer while both atoms are not in that layer.
      End Do
   End Subroutine remake_oldc
   
   subroutine process_iatcu(iatcu, coords, pl, cf, ifoldunit,dis,cl)
      Use interatcur
      Implicit None
      Type (confdata), Intent(Inout), Target :: cf
      type(t_interat_currents), Intent(In) :: iatcu
      type(planes_system), Intent(In) :: pl
      type(spatialquant), Intent(In) :: coords
      logical, intent(in) :: ifoldunit , dis
      Type(cell_system), Intent(InOut) :: cl
      type(spadirquant), Pointer :: xcur, ycur, zcur
      double precision :: base(2, 2), cond
      Integer :: itr, i1, j, ip, ia1, ia2, ip1, ip2, ia
      Type(t_curtensor), Pointer :: ctp
      Double precision :: xyof(2), factor, z1, z2, beta1, beta2, area, vol, r1(3), r2(3)
      Double precision, allocatable :: cmat(:, :)
      Double precision, allocatable :: zl(:), zr(:)
      type(spadirquant) :: cons
      Call putin1mat_spaqu(coords, cmat)
      Call get_zlr(pl, cf%z, zl, zr)
!-------------setting up grid essentials for discrete current scheme-----------------------
      if(dis) then
        if (allocated(cl%pts) .or. allocated(cl%zlr) ) then
        Deallocate(cl%pts,cl%zlr)
        end if
         allocate(cl%pts(coords%ldim,coords%nl+coords%nr+coords%nm))
         allocate(cl%zlr(2,size(zl)))
         cl%pts =cmat
        cl%zlr(1,:)=zl
        cl%zlr(2,:)=zr
        cl%info(3)=pl%m%np
        cl%nl=coords%nl
        cl%nm=coords%nm
        cl%nr=coords%nr
        cl%nch=iatcu%nch
      end if
!-----------------------------------------
      base = cf%base
      cond = cf%cond
      xcur => cf%xcur
      ycur => cf%ycur
      zcur => cf%zcur
      
      Call alloc_spadirqu(zcur, pl%m%np, pl%l%np, pl%r%np, iatcu%nch)
      zcur%Lm = 0.0d0
      zcur%Rm = 0.0d0
      Call alloc_spadirqu(ycur, pl%m%np, pl%l%np, pl%r%np, iatcu%nch)
      ycur%Lm = 0.0d0
      ycur%Rm = 0.0d0
      Call alloc_spadirqu(xcur, pl%m%np, pl%l%np, pl%r%np, iatcu%nch)
      xcur%Lm = 0.0d0
      xcur%Rm = 0.0d0
      Call alloc_spadirqu(cons, pl%m%nat, pl%l%nat, pl%r%nat, 1)
      Call alloc_spadirqu(cf%cons, pl%m%np, pl%l%np, pl%r%np, 1)
      Do ia = 1, pl%l%nat
         cons%Lir(1, ia) = iatcu%cons(1, ia)
         cons%Rt(1, ia) = iatcu%cons(2, ia)
      End Do
      Do ia = 1, pl%m%nat
         cons%Lm(1, ia) = iatcu%cons(1, ia + pl%l%nat)
         cons%Rm(1, ia) = iatcu%cons(2, ia + pl%l%nat)
      End Do
      Do ia = 1, pl%r%nat
         cons%Lt(1, ia) = iatcu%cons(1, ia + pl%l%nat + pl%m%nat)
         cons%Rir(1, ia) = iatcu%cons(2, ia + pl%l%nat + pl%m%nat)
      End Do
      Call paverage (cons, cf%cons, pl)
      area = abs(base(1,1)*base(2,2) - base(1,2)*base(2,1))
      Do itr = 1, iatcu%ntr
         ctp => iatcu%curten(itr)
         xyof = base(:,1)*iatcu%trlist(1, itr) + base(:,2)*iatcu%trlist(2, itr)
         Do i1 = 1, iatcu%nat
            ia1 = i1 - pl%l%nat
            z1 = cmat(3, i1)
            r1 = cmat(:, i1)
            If ( ia1 <= 0 ) Then
               ip1 = 0
            Elseif ( ia1 > pl%m%nat) Then
               ip1 = pl%m%nat + 1
            Else
               ip1 = pl%m%ind(ia1)
            End If
            Do j = ctp%ir(i1), ctp%ir(i1+1)-1
               ia2 = ctp%jc(j) - pl%l%nat
               z2 = cmat(3, ctp%jc(j))
               r2 = cmat(:, ctp%jc(j))
               r2(1:2) = r2(1:2) + xyof
               If ( ia2 <= 0 ) Then
                  ip2 = 0
               Elseif ( ia2 > pl%m%nat) Then
                  ip2 = pl%m%nat + 1
               Else
                  ip2 = pl%m%ind(ia2)
               End If
!               if (z1==z2) ip2=ip1
               Do ip = 1, pl%m%np
                  If ( ip >= min(ip1, ip2) .and. ip <= max(ip1, ip2)) Then
                     If ( ip == ip1) Then
                        beta1 = 0.0d0
                     Else
                        beta1 = min(abs(zl(ip)-z1), abs(zr(ip)-z1)) / abs(z1-z2)
                     End If
                     If ( ip == ip2) Then
                        beta2 = 0.0d0
                     Else
                        beta2 = min(abs(zl(ip)-z2), abs(zr(ip)-z2)) / abs(z1-z2)
                     End If
                     vol = (zr(ip)-zl(ip))
                     if (ifoldunit) then
                       factor = 0.5d0*( (1.0d0-beta1)**2 - beta2**2 ) / vol 
                     else
                       factor = 0.5d0*( (1.0d0-beta1)**2 - beta2**2 ) / vol / cond
                     end if

                     zcur%Lm(:,ip) = zcur%Lm(:,ip) + factor*(r1(3)-r2(3)) *ctp%a(j, :, 1)
                     zcur%Rm(:,ip) = zcur%Rm(:,ip) + factor*(r1(3)-r2(3)) *ctp%a(j, :, 2)
                     ycur%Lm(:,ip) = ycur%Lm(:,ip) + factor*(r1(2)-r2(2)) *ctp%a(j, :, 1)
                     ycur%Rm(:,ip) = ycur%Rm(:,ip) + factor*(r1(2)-r2(2)) *ctp%a(j, :, 2)
                     xcur%Lm(:,ip) = xcur%Lm(:,ip) + factor*(r1(1)-r2(1)) *ctp%a(j, :, 1)
                     xcur%Rm(:,ip) = xcur%Rm(:,ip) + factor*(r1(1)-r2(1)) *ctp%a(j, :, 2)
                        !------------------check NAN-----------------
                        if (isnan(1/vol)) stop '"vol" is a NaN'
                        !--------------------------------------------
                  End If
               End Do
            End Do
         End Do
      End Do
   End subroutine process_iatcu

   Subroutine pscreen_cur(cf)
      Implicit None
      Type (confdata), Intent (Inout) :: cf
      Integer :: ip
      
      call alloc_spadirqu(cf%xcur_ps, cf%xcur%nm, cf%xcur%nl, cf%xcur%nr, cf%xcur%ldim)
      call alloc_spadirqu(cf%ycur_ps, cf%ycur%nm, cf%ycur%nl, cf%ycur%nr, cf%ycur%ldim)
      call alloc_spadirqu(cf%zcur_ps, cf%zcur%nm, cf%zcur%nl, cf%zcur%nr, cf%zcur%ldim)
      do ip = 1, cf%ts%nm
         call pscr( cf%xcur%Lm(:,ip), cf%xcur%Rm(:,ip), cf%ts%Lm(7,ip), cf%ts%Rm(7,ip), cf%xcur_ps%Lm(:,ip), cf%xcur_ps%Rm(:,ip) )
         call pscr( cf%ycur%Lm(:,ip), cf%ycur%Rm(:,ip), cf%ts%Lm(7,ip), cf%ts%Rm(7,ip), cf%ycur_ps%Lm(:,ip), cf%ycur_ps%Rm(:,ip) )
         call pscr( cf%zcur%Lm(:,ip), cf%zcur%Rm(:,ip), cf%ts%Lm(7,ip), cf%ts%Rm(7,ip), cf%zcur_ps%Lm(:,ip), cf%zcur_ps%Rm(:,ip) )
      end do
      call alloc_spaqu(cf%SHang, cf%xcur%nm, cf%xcur%nl, cf%xcur%nr, 1)
      do ip = 1, cf%SHang%nm
         cf%SHang%m(1,ip) = ( cf%xcur_ps%Lm(2,ip) - cf%xcur_ps%Rm(2,ip) - cf%ycur_ps%Lm(1,ip) + cf%ycur_ps%Rm(1,ip) )/8.0d0
      end do
      
   End Subroutine pscreen_cur

   Subroutine cscreen_cur(cf)
      Implicit None
      Type (confdata), Intent (Inout) :: cf
      Integer :: ip
      
      call alloc_spadirqu(cf%xcur_ps, cf%xcur%nm, cf%xcur%nl, cf%xcur%nr, cf%xcur%ldim)
      call alloc_spadirqu(cf%ycur_ps, cf%ycur%nm, cf%ycur%nl, cf%ycur%nr, cf%ycur%ldim)
      call alloc_spadirqu(cf%zcur_ps, cf%zcur%nm, cf%zcur%nl, cf%zcur%nr, cf%zcur%ldim)
      do ip = 1, cf%ts%nm
         call pscr( cf%xcur%Lm(:,ip), cf%xcur%Rm(:,ip), cf%ts%Lm(7,ip), cf%ts%Rm(7,ip), cf%xcur_ps%Lm(:,ip), cf%xcur_ps%Rm(:,ip) )
         call pscr( cf%ycur%Lm(:,ip), cf%ycur%Rm(:,ip), cf%ts%Lm(7,ip), cf%ts%Rm(7,ip), cf%ycur_ps%Lm(:,ip), cf%ycur_ps%Rm(:,ip) )
         call pscr( cf%zcur%Lm(:,ip), cf%zcur%Rm(:,ip), cf%ts%Lm(7,ip), cf%ts%Rm(7,ip), cf%zcur_ps%Lm(:,ip), cf%zcur_ps%Rm(:,ip) )
      end do
      call alloc_spaqu(cf%SHang, cf%xcur%nm, cf%xcur%nl, cf%xcur%nr, 1)
      do ip = 1, cf%SHang%nm
         cf%SHang%m(1,ip) = ( cf%xcur_ps%Lm(2,ip) - cf%xcur_ps%Rm(2,ip) - cf%ycur_ps%Lm(1,ip) + cf%ycur_ps%Rm(1,ip) )/8.0d0
      end do
      
   End Subroutine cscreen_cur
   
   Subroutine get_zlr(pl, zav, zl, zr)
      Implicit None
      type(planes_system) :: pl
      type(spatialquant) :: zav
      double precision, allocatable :: zl(:), zr(:)
      
      Integer :: ip
      
      Allocate (zl(pl%m%np), zr(pl%m%np))
      zl(1) = 0.5d0*(zav%m(1, 1) + maxval(zav%l))
      Do ip = 1, pl%m%np-1
         zr(ip) = 0.5d0*(zav%m(1, ip) + zav%m(1, ip+1))
         zl(ip+1) = zr(ip)
      End Do
      zr(pl%m%np) = 0.5d0*(zav%m(1, pl%m%np) + minval(zav%r))
   End Subroutine get_zlr
   
   Subroutine make_all_plots(cf, cfdir)
      Implicit None
      type(confdata), Intent(In) :: cf
      Character(*), Intent(In) :: cfdir
      integer :: i 
      call writeplot(cf%z_si, cf%mdir, 'magdir', cfdir)
      call writeplot(cf%z_si, cf%spindos, 'ldos', cfdir)
      call writeplot(cf%z_si, cf%spinmu, 'mu', cfdir)
      call writeplot(cf%z_si, cf%ts, 'torq', cfdir, [(i,i=1,3)])
      call writeplot(cf%z_si, cf%ts, 'dens', cfdir, [(i,i=4,7)])
      call writeplot(cf%z_si, cf%ts_pscr, 'torq_pscr', cfdir, [(i,i=1,3)])
      call writeplot(cf%z_si, cf%ts_pscr, 'dens_pscr', cfdir, [(i,i=4,7)])
      if (cf%do_OAMdens /= 0) then
        call writeplot(cf%z_si, cf%OAMdens, 'OAM_dens', cfdir, [(i,i=1,3)])
        call writeplot(cf%z_si, cf%OAMdens_pscr, 'OAM_dens_pscr', cfdir, [(i,i=1,3)])
      endif
      if (cf%do_iatcu /= 0) then
         call writeplot(cf%z_si, cf%xcur%Lm, 'Lxcur', cfdir)
         call writeplot(cf%z_si, cf%xcur%Rm, 'Rxcur', cfdir)
         call writeplot(cf%z_si, cf%ycur%Lm, 'Lycur', cfdir)
         call writeplot(cf%z_si, cf%ycur%Rm, 'Rycur', cfdir)
         call writeplot(cf%z_si, cf%zcur%Lm, 'Lzcur', cfdir)
         call writeplot(cf%z_si, cf%zcur%Rm, 'Rzcur', cfdir)
         call writeplot(cf%z_si, cf%cons%Lm, 'Lcons', cfdir)
         call writeplot(cf%z_si, cf%cons%Rm, 'Rcons', cfdir)
         
         call writeplot(cf%z_si, cf%xcur_ps%Lm, 'Lxcur_ps', cfdir)
         call writeplot(cf%z_si, cf%xcur_ps%Rm, 'Rxcur_ps', cfdir)
         call writeplot(cf%z_si, cf%ycur_ps%Lm, 'Lycur_ps', cfdir)
         call writeplot(cf%z_si, cf%ycur_ps%Rm, 'Rycur_ps', cfdir)
         call writeplot(cf%z_si, cf%zcur_ps%Lm, 'Lzcur_ps', cfdir)
         call writeplot(cf%z_si, cf%zcur_ps%Rm, 'Rzcur_ps', cfdir)
         
         call writeplot(cf%z_si, cf%SHang%m, 'SHang', cfdir)
      end if
      open(100, file=trim(cfdir)//'/qfin')
      write(100,'(3I10)') cf%z_si%nl, cf%z_si%nm, cf%z_si%nr
      write(100,'(g17.8e3)') cf%cond
      write(100,'(g17.8e3)') cf%area_nm
      close(100)
   End Subroutine
   
   
   Subroutine writeplot_spaqu(z, q, name, cfdir, select_in)
      Implicit None
      type(spatialquant) :: z
      type(spatialquant) :: q
      Character(*) :: name, cfdir
      Integer, Optional :: select_in(:)
      
      integer :: i, nq
      Integer, allocatable :: select(:)
      Character(len=20) :: formatstring
      
      If (present(select_in)) Then
         nq = size(select_in)
         allocate(select(nq))
         select(:) = select_in(:)
      Else
         nq = q%ldim
         allocate(select(nq))
         select(:) = (/(I, I=1, nq)/)
      End If
      write(formatstring, '("(", I4, "g17.8e3)")') nq + 1
      open(100, file=trim(cfdir)//'/plots/'//name//'.dat')
      do i=1,z%nl
         write(100, formatstring) z%l(1, i), q%l(select, i)
      end do
      do i=1, z%nm
         write(100, formatstring) z%m(1, i), q%m(select, i)
      end do
      do i=1,z%nr
         write(100, formatstring) z%r(1, i), q%r(select, i)
      end do
      close(100)
   End Subroutine writeplot_spaqu
   
   Subroutine writeplot_spadirqu(z, q, name, cfdir, select_in)
      Implicit None
      type(spatialquant) :: z
      type(spadirquant) :: q
      Character(*) :: name, cfdir
      Integer, Optional :: select_in(:)
      
      integer :: i, nq
      Integer, allocatable :: select(:)
      Character(len=20) :: formatstring
      
      If (present(select_in)) Then
         nq = size(select_in)
         allocate(select(nq))
         select(:) = select_in(:)
      Else
         nq = q%ldim
         allocate(select(nq))
         select(:) = (/(I, I=1, nq)/)
      End If
      write(formatstring, '("(", I4, "g17.8e3)")') nq + 1
      open(100, file=trim(cfdir)//'/plots/L'//name//'.dat')
      open(200, file=trim(cfdir)//'/plots/R'//name//'.dat')
      do i=1,z%nl
         write(100, formatstring) z%l(1, i), q%Lir(select, i)
         write(200, formatstring) z%l(1, i), q%Rt(select, i)
      end do
      do i=1, z%nm
         write(100, formatstring) z%m(1, i), q%Lm(select, i)
         write(200, formatstring) z%m(1, i), q%Rm(select, i)
      end do
      do i=1,z%nr
         write(100, formatstring) z%r(1, i), q%Lt(select, i)
         write(200, formatstring) z%r(1, i), q%Rir(select, i)
      end do
      close(100)
      close(200)
   End Subroutine writeplot_spadirqu
   
   Subroutine writeplot_rdense(z, rdmat, name, cfdir)
      Use sparselib
      Implicit None
      type(spatialquant) :: z
      type(rdensemat) :: rdmat
      Character(*) :: name, cfdir
      
      call writeplot_mat(z, rdmat%bl, name, cfdir)
   End subroutine writeplot_rdense
   
   subroutine writeplot_mat(z, mat, name, cfdir)
      Implicit None
      type(spatialquant) :: z
      Double Precision :: mat(:, :)
      Character(*) :: name, cfdir
      
      integer :: i
      Character(len=20) :: formatstring
      
      write(formatstring, '("(", I4, "g17.8e3)")') size(mat, 1) + 1
      open(100, file=trim(cfdir)//'/plots/'//name//'.dat')
      do i=1, z%nm
         write(100, formatstring) z%m(1, i), mat(:, i)
      end do
      close(100)  
   end subroutine writeplot_mat
   
   subroutine putin1mat_spaqu(q,mat)
      Implicit None
      double precision, allocatable :: mat(:,:)
      type(spatialquant) :: q
      
      allocate(mat(q%ldim,q%nl+q%nr+q%nm))
      mat(:,1:q%nl) = q%l
      mat(:,q%nl+1:q%nl+q%nm) = q%m
      mat(:,q%nl+q%nm+1:q%nl+q%nr+q%nm) = q%r
   end subroutine putin1mat_spaqu


   subroutine putin1mat_spadirqu(q, Lmat, Rmat)
      Implicit None
      double precision, allocatable :: Lmat(:, :)
      double precision, allocatable, optional :: Rmat(:, :)
      type(spadirquant) :: q
      
      allocate(Lmat(q%ldim, q%nl+q%nr+q%nm))
      Lmat(:, 1:q%nl) = q%Lir
      Lmat(:, q%nl+1:q%nl+q%nm) = q%Lm
      Lmat(:, q%nl+q%nm+1:q%nl+q%nr+q%nm) = q%Lt
      If (present(Rmat)) Then
         allocate(Rmat(q%ldim, q%nl+q%nr+q%nm))
         Rmat(:, 1:q%nl) = q%Rt
         Rmat(:, q%nl+1:q%nl+q%nm) = q%Rm
         Rmat(:, q%nl+q%nm+1:q%nl+q%nr+q%nm) = q%Rir
      End If
   end subroutine putin1mat_spadirqu
   
   
   function var(array,mean)
      implicit none
      double precision :: var
      double precision :: array(:,:), mean
      integer :: i1, i2
      var=0.0d0
      do i1=1,size(array,1)
         do i2=1,size(array,2)
            var = var+(array(i1,i2)-mean)**2
         end do
      end do
   end function var
   
   
   subroutine getweights_mat(qin,qav,p,we)
      Implicit none
      Double Precision,intent(in) :: qin(:,:)
      Double Precision, intent(in) :: qav(:,:)
      Double Precision, intent(out) :: we(:,:)
      type(planes_region), intent(in) :: p
      
      Integer :: ip, ia
      Double Precision :: vari(size(we,1),size(we,2))
   
      vari = 0.0d0
      do ia = 1, p%nat
         ip = p%ind(ia)
         vari(:,ip) = vari(:,ip) + (qin(:,ia)-qav(:,ip))**2
      end do
      do ip = 1, p%np
         if (p%napp(ip) > 1) then
            we(:,ip) = dble(p%napp(ip)-1)/vari(:,ip)
         else
            we(:,ip) = 1.0d0
         end if
      end do
   end subroutine getweights_mat
   
   
   subroutine getweights_spaqu(q,qav,pl,we)
      Implicit None
      type(spatialquant), intent(in) :: q
      type(spatialquant), intent(in) :: qav
      type(spatialquant), intent(out) :: we
      type(planes_system), intent(in) :: pl
         
      call alloc_spaqu(we,pl%m%np,pl%l%np,pl%r%np,q%ldim)
      call getweights(q%m,qav%m,pl%m,we%m)
      call getweights(q%l,qav%l,pl%l,we%l)
      call getweights(q%r,qav%r,pl%r,we%r)
   end subroutine getweights_spaqu
   
   
   subroutine getweights_spadirqu(q,qav,pl,we)
      Implicit None
      type(spadirquant), intent(in) :: q
      type(spadirquant), intent(in) :: qav
      type(spadirquant), intent(out) :: we
      type(planes_system), intent(in) :: pl
         
      call alloc_spadirqu(we,pl%m%np,pl%l%np,pl%r%np,q%ldim)
      call getweights(q%Lm,qav%Lm,pl%m,we%Lm)
      call getweights(q%Rm,qav%Rm,pl%m,we%Rm)
      call getweights(q%Lir,qav%Lir,pl%l,we%Lir)
      call getweights(q%Lt,qav%Lt,pl%r,we%Lt)
      call getweights(q%Rir,qav%Rir,pl%r,we%Rir)
      call getweights(q%Rt,qav%Rt,pl%l,we%Rt)
   end subroutine getweights_spadirqu
   
   
    subroutine paverage_mat(qin,qout,p,we)
      Implicit none
      Double Precision, intent(in) :: qin(:,:)
      Double Precision, intent(out) :: qout(:,:)
      Double Precision, allocatable, optional, intent(out) :: we(:,:)
      type(planes_region), intent(in) :: p
      
      Integer :: ia, ip
      
      qout = 0.0d0
      do ia = 1, p%nat
         ip = p%ind(ia)
         qout(:,ip) = qout(:,ip) + qin(:,ia)
      end do
      do ip = 1, p%np
         qout(:,ip) = qout(:,ip)/dble(p%napp(ip))
      end do
      
      if (present(we)) then
         call getweights(qin, qout, p, we)
      end if
    end subroutine paverage_mat

   subroutine paverage_rdensemat(q,qav,pl,we)
      Use sparselib
      Implicit None
      type(rdensemat), intent(in) :: q
      type(rdensemat), intent(out) :: qav
      type(planes_system), intent(in) :: pl
      type(rdensemat), optional, intent(out) :: we
      
      Call alloc (qav, q%nrow, pl%m%np)
      Call paverage(q%bl, qav%bl, pl%m)
      If (present(we)) Then
         Call alloc (we, q%nrow, pl%m%np)
         Call getweights(q%bl, qav%bl, pl%m, we%bl)
      End If
   end subroutine paverage_rdensemat

   subroutine paverage_spaqu(q,qav,pl,we)
      Implicit None
      type(spatialquant), intent(in) :: q
      type(spatialquant), intent(out) :: qav
      type(planes_system), intent(in) :: pl
      type(spatialquant), optional, intent(out) :: we
      
      call alloc_spaqu(qav,pl%m%np,pl%l%np,pl%r%np,q%ldim)
      call paverage(q%m,qav%m,pl%m)
      call paverage(q%l,qav%l,pl%l)
      call paverage(q%r,qav%r,pl%r)
      if (present(we)) call getweights(q,qav,pl,we)
   end subroutine paverage_spaqu
   
   
   subroutine paverage_spadirqu(q,qav,pl,we)
      Implicit None
      type(spadirquant), intent(in) :: q
      type(spadirquant), intent(out) :: qav
      type(planes_system), intent(in) :: pl
      type(spadirquant), optional, intent(out) :: we
         
      call alloc_spadirqu(qav,pl%m%np,pl%l%np,pl%r%np,q%ldim)
      call paverage(q%Lm,qav%Lm,pl%m)
      call paverage(q%Rm,qav%Rm,pl%m)
      call paverage(q%Lir,qav%Lir,pl%l)
      call paverage(q%Lt,qav%Lt,pl%r)
      call paverage(q%Rir,qav%Rir,pl%r)
      call paverage(q%Rt,qav%Rt,pl%l)
      if (present(we)) call getweights(q,qav,pl,we)
   end subroutine paverage_spadirqu


   Subroutine getplanes2(z,maxvar,pl)
      Implicit None
      double precision, intent(in) :: z(:, :), maxvar
      type(planes_region), intent(out) :: pl
      
      integer :: ip, ia, np
      double precision :: dum(2,size(z,2)), d
      
      pl%nat = size(z,2)
      allocate(pl%ind(pl%nat))
      np = 1
      pl%ind(1) = 1
      dum(1,1) = z(1,1)
      dum(2,1) = 0.0d0
      do ia=2, pl%nat
         ip = minloc( abs( dum(1,1:np)-z(1,ia) ), 1)
         d = z(1,ia)-dum(1,ip)
         if (abs(d) > maxvar-dum(2,ip)) then
            np = np + 1
            pl%ind(ia) = np
            dum(1,np) = z(1,ia)
            dum(2,np) = 0.0d0
         else
            pl%ind(ia) = ip
            if (abs(d) > dum(2,ip)) then
               dum(1,ip) = dum(1,ip) + 0.5*d - 0.5*sign(dum(2,ip),d)
               dum(2,ip) = 0.5*dum(2,ip) + 0.5*abs(d)
            end if
         end if
      end do
      pl%np = np
      allocate(pl%napp(np))
      pl%napp = 0
      do ia = 1, pl%nat
         pl%napp(pl%ind(ia)) = pl%napp(pl%ind(ia)) + 1
      end do
      allocate(pl%order(np))
   end subroutine getplanes2
   
   Subroutine getplanes3(z,maxvar,pl)
! modified by Icer
! date : 2014/09/28
! fix plane
      Implicit None
      double precision, intent(in) :: z(:, :), maxvar
      type(planes_region), intent(out) :: pl
      
      integer ::  ia, np, ip, newmaxvar
      double precision :: dum1, d, dum2, dplane
      
      if (maxvar==0) then
        write(*,'(xxx,A)')'Worning! the value of Pav can not be Zero in Fix plane situation!!!'
        stop
      end if

      if (abs(maxvar-nint(maxvar))>0.000001) then
        write(*,'(xxx,A)')'Input value of Pav is not integer!'
      end if

      newmaxvar=nint(maxvar)+1

      write(*,*)'The value of Pav will be treat as Pav= ',newmaxvar

      if (newmaxvar<10) then
        write(*,'(xxx,A)')'Worning! The number of planes is LESS than 10, please careful with this.'
      end if

      pl%nat = size(z,2)
      allocate(pl%ind(pl%nat))
      dum1=minval(z(1,:))-0.001d0
      dum2=maxval(z(1,:))+0.001d0
      dplane=(dum2-dum1)/real(newmaxvar)
      do ia=1, pl%nat
         d = z(1,ia)-dum1
         ip=floor(d/dplane)
         pl%ind(ia) = ip+1
      end do
      pl%np = maxval(pl%ind(:))
      np = pl%np
      allocate(pl%napp(np))
      pl%napp = 0
      do ia = 1, pl%nat
         pl%napp(pl%ind(ia)) = pl%napp(pl%ind(ia)) + 1
      end do
      allocate(pl%order(np))
   end subroutine getplanes3
   
   Subroutine getplanes(z, zav, maxvar, pl)
      Implicit None
      type(spatialquant), intent(in) :: z
      type(spatialquant), intent(out) :: zav
      double precision, intent(in) :: maxvar
      type(planes_system), intent(out) :: pl
      
      call getplanes2(z%m,maxvar,pl%m)
      call getplanes2(z%l,maxvar,pl%l)
      call getplanes2(z%r,maxvar,pl%r)
      call paverage(z,zav,pl)
      call orderplanes(pl%m, zav%m(1, :))
      call orderplanes(pl%l, zav%l(1, :))
      call orderplanes(pl%r, zav%r(1, :))
      call getthicknesses(pl, zav)
   end subroutine getplanes
   
   Subroutine fixplanes(z, zav, maxvar, pl)
! modified by Icer
! date : 2014/09/28
! fix plane
      Implicit None
      type(spatialquant), intent(in) :: z
      type(spatialquant), intent(out) :: zav
      double precision, intent(in) :: maxvar
      type(planes_system), intent(out) :: pl
      
      call getplanes3(z%m,maxvar,pl%m)
      call getplanes2(z%l,maxvar,pl%l)
      call getplanes2(z%r,maxvar,pl%r)
      call paverage(z,zav,pl)
      call orderplanes(pl%m, zav%m(1, :))
      call orderplanes(pl%l, zav%l(1, :))
      call orderplanes(pl%r, zav%r(1, :))
      call getthicknesses(pl, zav)
   end subroutine fixplanes
   
   Subroutine getthicknesses(p, zav)
      Implicit None
      type(spatialquant), intent(in) :: zav
      type(planes_system), intent(inout) :: p
         
      Integer :: i, np
      np = p%m%np
      allocate(p%th(np))
      p%th(:) = 0.0d0
      p%th(1) = zav%m(1, 1) - zav%l(1, zav%nl)
      Do i = 1, np-1
         p%th(i:i+1) = p%th(i:i+1) + (zav%m(1, i+1) - zav%m(1, i))
      End Do
      p%th(np) = p%th(np) + zav%r(1, 1) - zav%m(1, np)
      p%th(:) = 0.5d0*p%th(:)
   End Subroutine
   
   Subroutine orderplanes(p, zav)
      Use sparselib
      Implicit None
      double precision, intent(inout) :: zav(:)
      type(planes_region), intent(inout) :: p
      
      Integer, allocatable :: order(:), newind(:), newnapp(:)
      Integer :: ia, ip
      double precision, allocatable :: newzav(:)
      
      allocate(order(p%np))
      allocate(newind(p%nat))
      !allocate(newzav(p%nat))
      allocate(newzav(p%np))
      allocate(newnapp(p%nat))
      call dqsort(p%np, zav, order)
      Do ip = 1, p%np
         newnapp(order(ip)) = p%napp(ip)
         newzav(order(ip)) = zav(ip)
         Do ia = 1, p%nat
            If (order(ip) == p%ind(ia)) Then
               newind(ia) = ip
            End If
         End Do
      End Do
      p%ind = newind
      p%napp = newnapp
      zav = newzav
   End Subroutine

   Function suminv(array)
      Implicit None
      double precision :: suminv
      double precision :: array(:)
      Integer ::  sz, i

      sz = size(array)
      suminv = 0.0
      do i=1,sz
         suminv = suminv + abs(1.0/array(i))
      end do
   End Function suminv
   
   Subroutine get_magnetization(magnet, rotm, magmom, q)
      Use sparselib
      Implicit None
      Type(spatialquant), Intent(inout) :: magnet
      Type(spatialquant), Intent(in) :: q
      Type(rdensemat), Intent(in) :: rotm
      type(rdensemat), Intent(in) :: magmom
         
      Double precision :: m(3)
      Integer :: ia
     
      Call alloc_spaqu(magnet, q%nm, q%nl, q%nr, 3)
      magnet%l(3, :) = magmom%bl(1, 1:q%nl)
      magnet%m(3, :) = magmom%bl(1, q%nl+1:q%nl+q%nm)
      magnet%r(3, :) = magmom%bl(1, q%nl+q%nm+1:q%nl+q%nm+q%nr)
      
      call getm(rotm, 1, m)
      do ia = 1, magnet%nl
         magnet%l(:, ia) = m*magnet%l(3, ia)
      end do
      do ia = 1, magnet%nm
         call getm(rotm,ia,m)
         magnet%m(:, ia) = m*magnet%m(3, ia)
      end do
      call getm(rotm, magnet%nm, m)
      do ia = 1, magnet%nr
         magnet%r(:, ia) = m*magnet%r(3, ia)
      end do
   End Subroutine

   Subroutine mudos(ts, rotm, spinmu, spindos)
      Use sparselib
      Implicit None
      Type(spadirquant), Intent (in) :: ts
      Type(rdensemat), Intent (in) :: rotm
      Type(spadirquant), Intent (out) :: spinmu
      Type(spatialquant), Intent (out) :: spindos

      Double precision :: m(3)
      Integer :: ia
     
      Call alloc_spaqu(spindos, ts%nm, ts%nl, ts%nr, 2)
      Call alloc_spadirqu(spinmu, ts%nm, ts%nl, ts%nr, 2)
      
      call getm(rotm,1,m)
      do ia = 1,ts%nl
         call mudos1(ts%Lir(:,ia),ts%Rt(:,ia),m,spinmu%Lir(:,ia),spinmu%Rt(:,ia),spindos%l(:,ia))
      end do
      do ia = 1,ts%nm
         call getm(rotm,ia,m)
         call mudos1(ts%Lm(:,ia),ts%Rm(:,ia),m,spinmu%Lm(:,ia),spinmu%Rm(:,ia),spindos%m(:,ia))
      end do
      call getm(rotm,ts%nm,m)
      do ia = 1,ts%nr
         call mudos1(ts%Lt(:,ia),ts%Rir(:,ia),m,spinmu%Lt(:,ia),spinmu%Rir(:,ia),spindos%r(:,ia))
      end do
   End Subroutine mudos

   Subroutine mudos1(Lts,Rts,m,Lspinmu,Rspinmu,spindos)
      Implicit None
      Double precision, Intent (in) :: Lts(7), Rts(7), m(3)
      Double precision, Intent (out) :: Lspinmu(2), Rspinmu(2), spindos(2)

      Double precision :: Lsdotm, Rsdotm, sdotm, dens

      Lsdotm = Lts(4)*m(1)+Lts(5)*m(2)+Lts(6)*m(3)
      Rsdotm = Rts(4)*m(1)+Rts(5)*m(2)+Rts(6)*m(3)
      sdotm = Lsdotm + Rsdotm
      dens = Lts(7)+Rts(7)
      spindos(1) = 0.5*(dens + sdotm)
      spindos(2) = 0.5*(dens - sdotm)
      Lspinmu(1) = (Lts(7) + Lsdotm)/(dens+sdotm)
      Lspinmu(2) = (Lts(7) - Lsdotm)/(dens-sdotm)
      Rspinmu(1) = (Rts(7) + Rsdotm)/(dens+sdotm)
      Rspinmu(2) = (Rts(7) - Rsdotm)/(dens-sdotm)
   End subroutine mudos1

   Subroutine mudos_cscr(ts, rotm, spinmu, spindos)
      Use sparselib
      Implicit None
      Type(spadirquant), Intent (in) :: ts
      Type(rdensemat), Intent (in) :: rotm
      Type(spadirquant), Intent (out) :: spinmu
      Type(spatialquant), Intent (out) :: spindos

      Double precision :: m(3)
      Integer :: ia
     
      Call alloc_spaqu(spindos, ts%nm, ts%nl, ts%nr, 2)
      Call alloc_spadirqu(spinmu, ts%nm, ts%nl, ts%nr, 2)
      
      call getm(rotm,1,m)
      do ia = 1,ts%nl
         call mudos1_cscr(ts%Lir(:,ia),ts%Rt(:,ia),m,spinmu%Lir(:,ia),spinmu%Rt(:,ia),spindos%l(:,ia))
      end do
      do ia = 1,ts%nm
         call getm(rotm,ia,m)
         call mudos1_cscr(ts%Lm(:,ia),ts%Rm(:,ia),m,spinmu%Lm(:,ia),spinmu%Rm(:,ia),spindos%m(:,ia))
      end do
      call getm(rotm,ts%nm,m)
      do ia = 1,ts%nr
         call mudos1_cscr(ts%Lt(:,ia),ts%Rir(:,ia),m,spinmu%Lt(:,ia),spinmu%Rir(:,ia),spindos%r(:,ia))
      end do
   End Subroutine mudos_cscr

   Subroutine mudos1_cscr(Lts,Rts,m,Lspinmu,Rspinmu,spindos)
      Implicit None
      Double precision, Intent (in) :: Lts(7), Rts(7), m(3)
      Double precision, Intent (out) :: Lspinmu(2), Rspinmu(2), spindos(2)

      Double precision :: Lsdotm, Rsdotm, sdotm, dens

      Lsdotm = Lts(4)*m(1)+Lts(5)*m(2)+Lts(6)*m(3)
      Rsdotm = Rts(4)*m(1)+Rts(5)*m(2)+Rts(6)*m(3)
      sdotm = Lsdotm + Rsdotm
      dens = Lts(7)+Rts(7)
      spindos(1) = 0.5*(dens + sdotm)
      spindos(2) = 0.5*(dens - sdotm)
      Lspinmu(1) = (Lts(7) + Lsdotm)
      Lspinmu(2) = (Lts(7) - Lsdotm)
      Rspinmu(1) = (Rts(7) + Rsdotm)
      Rspinmu(2) = (Rts(7) - Rsdotm)
   End subroutine mudos1_cscr

   Subroutine mudos_noscr(ts, rotm, spinmu, spindos)
      Use sparselib
      Implicit None
      Type(spadirquant), Intent (in) :: ts
      Type(rdensemat), Intent (in) :: rotm
      Type(spadirquant), Intent (out) :: spinmu
      Type(spatialquant), Intent (out) :: spindos

      Double precision :: m(3)
      Integer :: ia
     
      Call alloc_spaqu(spindos, ts%nm, ts%nl, ts%nr, 2)
      Call alloc_spadirqu(spinmu, ts%nm, ts%nl, ts%nr, 2)
      
      call getm(rotm,1,m)
      do ia = 1,ts%nl
         call mudos1_noscr(ts%Lir(:,ia),ts%Rt(:,ia),m,spinmu%Lir(:,ia),spinmu%Rt(:,ia),spindos%l(:,ia))
      end do
      do ia = 1,ts%nm
         call getm(rotm,ia,m)
         call mudos1_noscr(ts%Lm(:,ia),ts%Rm(:,ia),m,spinmu%Lm(:,ia),spinmu%Rm(:,ia),spindos%m(:,ia))
      end do
      call getm(rotm,ts%nm,m)
      do ia = 1,ts%nr
         call mudos1_noscr(ts%Lt(:,ia),ts%Rir(:,ia),m,spinmu%Lt(:,ia),spinmu%Rir(:,ia),spindos%r(:,ia))
      end do
   End Subroutine mudos_noscr

   Subroutine mudos1_noscr(Lts,Rts,m,Lspinmu,Rspinmu,spindos)
      Implicit None
      Double precision, Intent (in) :: Lts(7), Rts(7), m(3)
      Double precision, Intent (out) :: Lspinmu(2), Rspinmu(2), spindos(2)

      Double precision :: Lsdotm, Rsdotm, sdotm, dens

      Lsdotm = Lts(4)*m(1)+Lts(5)*m(2)+Lts(6)*m(3)
      Rsdotm = Rts(4)*m(1)+Rts(5)*m(2)+Rts(6)*m(3)
      sdotm = Lsdotm + Rsdotm
      dens = Lts(7)+Rts(7)
      spindos(1) = 0.5*(dens + sdotm)
      spindos(2) = 0.5*(dens - sdotm)
      Lspinmu(1) = (Lts(7) + Lsdotm)
      Lspinmu(2) = (Lts(7) - Lsdotm)
      Rspinmu(1) = (Rts(7) + Rsdotm)
      Rspinmu(2) = (Rts(7) - Rsdotm)
   End subroutine mudos1_noscr
   
   Subroutine scale_ts(q)
      Implicit None
      Type(spadirquant), Intent(inout) :: q
      ! this makes all quantities have proper atomic units
      
      ! The output of the code in case of spin densities is 1/2 sigma.
      ! Since in the code h = 1 instead of hbar = 1, there is a factor of 2pi missing in the velocities.
      ! This becomes 1/2pi in the flux normalization.
      ! Combining this, we should scale the densities by 1/pi to get the number of spins and electrons per Hartree.
      q%Lir(4:7, :) = q%Lir(4:7, :)*DEF_M_1_PI
      q%Rt(4:7, :)  = q%Rt(4:7, :) *DEF_M_1_PI
      q%Lm(4:7, :)  = q%Lm(4:7, :) *DEF_M_1_PI
      q%Rm(4:7, :)  = q%Rm(4:7, :) *DEF_M_1_PI
      q%Rir(4:7, :) = q%Rir(4:7, :)*DEF_M_1_PI
      q%Lt(4:7, :)  = q%Lt(4:7, :) *DEF_M_1_PI
      
      ! In case of the torques, the factor 1/2 in sigma is cancelled by its omission in 'split'.
      ! The factor 1/2pi should still be put in to get the torque per unit energy (torquance) in atomic units.
      q%Lir(1:3, :) = 0.5d0*q%Lir(1:3, :)*DEF_M_1_PI
      q%Rt(1:3, :)  = 0.5d0*q%Rt(1:3, :) *DEF_M_1_PI
      q%Lm(1:3, :)  = 0.5d0*q%Lm(1:3, :) *DEF_M_1_PI
      q%Rm(1:3, :)  = 0.5d0*q%Rm(1:3, :) *DEF_M_1_PI
      q%Rir(1:3, :) = 0.5d0*q%Rir(1:3, :)*DEF_M_1_PI
      q%Lt(1:3, :)  = 0.5d0*q%Lt(1:3, :) *DEF_M_1_PI
   End Subroutine

   Subroutine getm(rotm,ia,m)
      Use sparselib
      Implicit None
      Integer, Intent (in) :: ia
      Type(rdensemat), Intent(in) :: rotm
      Double precision, Intent (out) :: m(3)
      
      Integer :: i

      m(1)=Sin(-rotm%bl(1,ia))*Cos(-rotm%bl(2,ia))
      
      m(2)=Sin(-rotm%bl(1,ia))*Sin(-rotm%bl(2,ia))
      m(3)=Cos(-rotm%bl(1,ia))
      Do i = 1, 3
         if(abs(m(i)) < 1.0e-14) m(i) = 0.0d0
      End Do
   End subroutine getm


   Subroutine pscreen(ts,ts_pscr)
      Implicit None
      type(spadirquant), intent(in) :: ts
      type(spadirquant), intent(out) :: ts_pscr
      integer :: ia

      Call alloc_spadirqu(ts_pscr,ts%nm,ts%nl,ts%nr,7)
      do ia = 1,ts%nl
         call pscr(ts%Lir(:,ia),ts%Rt(:,ia),ts%Lir(7,ia),ts%Rt(7,ia),ts_pscr%Lir(:,ia),ts_pscr%Rt(:,ia))
         !call pscr1(ts%Lir(:,ia),ts%Rt(:,ia),ts_pscr%Lir(:,ia),ts_pscr%Rt(:,ia))
      end do
      do ia = 1,ts%nm
         call pscr(ts%Lm(:,ia),ts%Rm(:,ia),ts%Lm(7,ia),ts%Rm(7,ia),ts_pscr%Lm(:,ia),ts_pscr%Rm(:,ia))
         !call pscr1(ts%Lm(:,ia),ts%Rm(:,ia),ts_pscr%Lm(:,ia),ts_pscr%Rm(:,ia))
      end do
      do ia = 1,ts%nr
         call pscr(ts%Lt(:,ia),ts%Rir(:,ia),ts%Lt(7,ia),ts%Rir(7,ia),ts_pscr%Lt(:,ia),ts_pscr%Rir(:,ia))
         !call pscr1(ts%Lt(:,ia),ts%Rir(:,ia),ts_pscr%Lt(:,ia),ts_pscr%Rir(:,ia))
      end do
   end Subroutine pscreen

   
   Subroutine pscr(L, R, nL, nR, Lscr, Rscr)
      Implicit None
      double precision, intent (in) :: L(:), R(:)
      double precision, intent (out) :: Lscr(:), Rscr(:)
      double precision, intent (in) :: nL, nR
      Integer :: dimL, i
      Double precision :: muL, muR
      
      muL = nL/(nL+nR)
      muR = nR/(nL+nR)

      dimL = size(L)
      Do i = 1, dimL
         Lscr(i) = L(i) - muL*(L(i)+R(i))
         Rscr(i) = R(i) - muR*(L(i)+R(i))
      End Do
   End Subroutine pscr

   Subroutine pscreen_oam_density(od,od_pscr,ts)
      Implicit None
      type(spadirquant), intent(in) :: od
      type(spadirquant), intent(out) :: od_pscr
      type(spadirquant), intent(in) :: ts
      integer :: ia
! od is shorthand for Orbital angular momentum Density
! ts are the torques and densities; the index 7 is the charge density
! (why it is 'ts' I don't know)

      Call alloc_spadirqu(od_pscr,od%nm,od%nl,od%nr,3)
      do ia = 1,od%nl
         call pscr(od%Lir(:,ia),od%Rt(:,ia),ts%Lir(7,ia),ts%Rt(7,ia),od_pscr%Lir(:,ia),od_pscr%Rt(:,ia))
         !call pscr1(od%Lir(:,ia),od%Rt(:,ia),od_pscr%Lir(:,ia),od_pscr%Rt(:,ia))
      end do
      do ia = 1,od%nm
         call pscr(od%Lm(:,ia),od%Rm(:,ia),ts%Lm(7,ia),ts%Rm(7,ia),od_pscr%Lm(:,ia),od_pscr%Rm(:,ia))
         !call pscr1(od%Lm(:,ia),od%Rm(:,ia),od_pscr%Lm(:,ia),od_pscr%Rm(:,ia))
      end do
      do ia = 1,od%nr
         call pscr(od%Lt(:,ia),od%Rir(:,ia),ts%Lt(7,ia),ts%Rir(7,ia),od_pscr%Lt(:,ia),od_pscr%Rir(:,ia))
         !call pscr1(od%Lt(:,ia),od%Rir(:,ia),od_pscr%Lt(:,ia),od_pscr%Rir(:,ia))
      end do
   end Subroutine pscreen_oam_density

   Subroutine cscreen(ts,ts_pscr)
      Implicit None
      type(spadirquant), intent(in) :: ts
      type(spadirquant), intent(out) :: ts_pscr
      integer :: ia

      Call alloc_spadirqu(ts_pscr,ts%nm,ts%nl,ts%nr,7)
      do ia = 1,ts%nl
         call cscr(ts%Lir(:,ia),ts%Rt(:,ia),ts%Lir(7,ia),ts%Rt(7,ia),ts_pscr%Lir(:,ia),ts_pscr%Rt(:,ia))
      end do
      do ia = 1,ts%nm
         call cscr(ts%Lm(:,ia),ts%Rm(:,ia),ts%Lm(7,ia),ts%Rm(7,ia),ts_pscr%Lm(:,ia),ts_pscr%Rm(:,ia))
      end do
      do ia = 1,ts%nr
         call cscr(ts%Lt(:,ia),ts%Rir(:,ia),ts%Lt(7,ia),ts%Rir(7,ia),ts_pscr%Lt(:,ia),ts_pscr%Rir(:,ia))
      end do
   end Subroutine cscreen

   
   Subroutine cscr(L, R, nL, nR, Lscr, Rscr)
      Implicit None
      double precision, intent (in) :: L(:), R(:)
      double precision, intent (out) :: Lscr(:), Rscr(:)
      double precision, intent (in) :: nL, nR
      Integer :: dimL, i
      Double precision :: muL, muR
      
      muL = nL
      muR = nR

      dimL = size(L)
      Do i = 1, dimL
         Lscr(i) = L(i) - muL*(L(i)+R(i))
         Rscr(i) = R(i) - muR*(L(i)+R(i))
      End Do
   End Subroutine cscr

!   Subroutine pscr1(Lts,Rts,Lts_pscr,Rts_pscr)
!      Implicit None
!      double precision, intent (in) :: Lts(7), Rts(7)
!      double precision, intent (out) :: Lts_pscr(7), Rts_pscr(7)
!      double precision :: ts0(7), mu
!
!      ts0 = Lts+Rts
!      mu = Lts(7)/ts0(7)
!      Lts_pscr(:) = Lts(:) - mu*ts0(:)
!      mu = Rts(7)/ts0(7)
!      Rts_pscr = Rts - mu*ts0
!   End Subroutine pscr1
   
end Module qofzlib
