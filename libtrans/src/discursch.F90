#include "math_def.h"
Module discursch
   Use quantlib
   Implicit None
contains
!------------------------------
   subroutine transform(cl, base, ptr)
      Implicit none
      Type(cell_system), Intent(InOut) :: cl
      Real(Kind=DEF_DBL_PREC), Intent(In)  :: base(:, :), ptr(:)
      Real(Kind=DEF_DBL_PREC)  :: sxy, sxz, syx, syz, x1, y1, scale
      Real(Kind=DEF_DBL_PREC), Allocatable  :: coord(:, :)
      Integer :: i
      Allocate (coord(2, cl%nl + cl%nm + cl%nr))
      if (allocated(cl%tr)) then
         deallocate (cl%tr)
      end if
      Allocate (cl%tr(3, cl%nl + cl%nm + cl%nr))
      sxy = base(1, 2)/base(2, 2)
      sxz = ptr(1)/ptr(3)
      syx = base(2, 1)/base(1, 1)
      syz = ptr(2)/ptr(3)
      scale = base(1, 1)/(base(1, 1) - sxy*base(2, 1))
      Do i = 1, cl%nl + cl%nm + cl%nr
         x1 = cl%pts(1, i) - sxz*cl%pts(3, i)
         y1 = cl%pts(2, i) - syz*cl%pts(3, i)
         coord(1, i) = (x1 - sxy*y1)*scale
         coord(2, i) = y1 - syx*coord(1, i)
      End Do
      cl%tr(1:2, :) = coord(1:2, :)
      cl%tr(3, :) = cl%pts(3, :)
      cl%base(1, 1) = base(1, 1)
      cl%base(2, 1) = 0.0
      cl%base(1, 2) = 0.0
      cl%base(2, 2) = base(2, 2)
      Deallocate (coord)
   end subroutine transform
!------------------------------

   subroutine process_gridcur(iatcu, ifoldunit, base, cond, cl, pl, cf_d)
      Use interatcur
      Implicit None
      logical, intent(in) :: ifoldunit
      type(t_interat_currents), Intent(In) :: iatcu
      Type(confdisdata), Intent(Inout) :: cf_d
      type(cell_system), Intent(InOut):: cl
      type(planes_system), Intent(In) :: pl
      Real(Kind=DEF_DBL_PREC), Intent(In) :: base(:, :), cond
      Type(t_curtensor), Pointer :: ctp
      Double precision :: xyof(2), xyof2(2), factor, vol, r1(3), r2(3), r11(3), r22(3)
      Double precision, allocatable :: xl(:), xr(:), yl(:), yr(:), zl(:), zr(:)
      Double precision :: beta1, beta2, t1, t2, area
      Integer :: itr, i1, j, ic, ia1, ia2, ic1, ic2, ip1, ip2, ip
      Call get_cell(cl, pl, xl, xr, yl, yr, zl, zr)
      call alloc_confdisdata(cf_d, size(zl), iatcu%nch)

      area = abs(base(1, 1)*base(2, 2) - base(1, 2)*base(2, 1))
      Do itr = 1, iatcu%ntr
         ctp => iatcu%curten(itr)
         xyof = base(:, 1)*iatcu%trlist(1, itr) + base(:, 2)*iatcu%trlist(2, itr)
         xyof2 = cl%base(:, 1)*iatcu%trlist(1, itr) + cl%base(:, 2)*iatcu%trlist(2, itr)
         Do i1 = 1, iatcu%nat
            ia1 = i1 - pl%l%nat

            If (ia1 <= 0) Then
               ip1 = 0
            Elseif (ia1 > pl%m%nat) Then
               ip1 = pl%m%nat + 1
            Else
               ip1 = pl%m%ind(ia1)
            End If

            r1 = cl%pts(:, i1)
            r11 = cl%tr(:, i1)
            ic1 = cl%ind(i1)
            Do j = ctp%ir(i1), ctp%ir(i1 + 1) - 1
               ia2 = ctp%jc(j) - pl%l%nat

               If (ia2 <= 0) Then
                  ip2 = 0
               Elseif (ia2 > pl%m%nat) Then
                  ip2 = pl%m%nat + 1
               Else
                  ip2 = pl%m%ind(ia2)
               End If

               r2 = cl%pts(:, ctp%jc(j))
               r22 = cl%tr(:, ctp%jc(j))
               r2(1:2) = r2(1:2) + xyof
               r22(1:2) = r22(1:2) + xyof2
               ic2 = cl%ind(ctp%jc(j))
               Do ic = 1, product(cl%info)
                  ip = cl%pl(ic)
                  If (ip >= min(ip1, ip2) .and. ip <= max(ip1, ip2)) Then
                     call inbox(xl(ic), xr(ic), yl(ic), yr(ic), zl(ic), zr(ic), r11, r22, t1, t2)
                     If (t1 >= 0.0d0 .and. t2 <= 1.0d0) Then
                        beta1 = t1
                        beta2 = 1.0d0 - t2
                        vol = abs(zr(ic) - zl(ic))
                        if (ifoldunit) then
                           factor = 0.5d0*((1.0d0 - beta1)**2 - beta2**2)/vol
                        else
                           factor = 0.5d0*((1.0d0 - beta1)**2 - beta2**2)/vol/cond
                        end if
                        cf_d%zcur%Lm(:, ic) = cf_d%zcur%Lm(:, ic) + factor*(r1(3) - r2(3))*ctp%a(j, :, 1)
                        cf_d%zcur%Rm(:, ic) = cf_d%zcur%Rm(:, ic) + factor*(r1(3) - r2(3))*ctp%a(j, :, 2)
                        cf_d%ycur%Lm(:, ic) = cf_d%ycur%Lm(:, ic) + factor*(r1(2) - r2(2))*ctp%a(j, :, 1)
                        cf_d%ycur%Rm(:, ic) = cf_d%ycur%Rm(:, ic) + factor*(r1(2) - r2(2))*ctp%a(j, :, 2)
                        cf_d%xcur%Lm(:, ic) = cf_d%xcur%Lm(:, ic) + factor*(r1(1) - r2(1))*ctp%a(j, :, 1)
                        cf_d%xcur%Rm(:, ic) = cf_d%xcur%Rm(:, ic) + factor*(r1(1) - r2(1))*ctp%a(j, :, 2)
                     End If
                  End If
!            if (isnan(factor)) then
!                print*,r1,r2,ip1,ip2,ip
!                stop '"factor" is a NaN'
!            end if
               End Do
            End Do
         End Do
      End Do
   end subroutine process_gridcur
!--------------
   subroutine alloc_confdisdata(qin, tan, ld)
      Implicit None
      Integer, Intent(In) :: tan, ld
      Type(confdisdata), Intent(Inout) :: qin
      if (qin%alloc /= 0) then
         call free_confdisdata(qin)
      end if
      if (qin%alloc == 0) then
         Allocate (qin%xcur%Lm(ld, tan))
         Allocate (qin%ycur%Lm(ld, tan))
         Allocate (qin%zcur%Lm(ld, tan))
         Allocate (qin%xcur%Rm(ld, tan))
         Allocate (qin%ycur%Rm(ld, tan))
         Allocate (qin%zcur%Rm(ld, tan))
         Allocate (qin%ts_av%Lm(1, tan))
         Allocate (qin%ts_av%Rm(1, tan))
         Allocate (qin%ts%Lm(ld, tan))
         Allocate (qin%ts%Rm(ld, tan))
         Allocate (qin%OAMdens%Lm(3, tan))
         Allocate (qin%OAMdens%Rm(3, tan))
         qin%xcur%Lm = 0.0d0
         qin%ycur%Lm = 0.0d0
         qin%zcur%Lm = 0.0d0
         qin%xcur%Rm = 0.0d0
         qin%ycur%Rm = 0.0d0
         qin%zcur%Rm = 0.0d0
         qin%ts%Lm = 0.0d0
         qin%ts%Rm = 0.0d0
         qin%OAMdens%Lm = 0.0d0
         qin%OAMdens%Rm = 0.0d0
         qin%ts_av%Lm = 0.0d0
         qin%ts_av%Rm = 0.0d0

         qin%ts%alloc = 1
         qin%OAMdens%alloc = 1
         qin%ts_av%alloc = 1
         qin%alloc = 1

         qin%nat = tan
         qin%ldim = ld
      end if
   end subroutine alloc_confdisdata
!--------------
   subroutine free_confdisdata(qin)
      Implicit None
      Type(confdisdata), Intent(Inout) :: qin
      if (qin%alloc == 0) return
      Deallocate (qin%xcur%Lm, qin%ycur%Lm, qin%zcur%Lm, qin%ts%Lm, qin%ts_av%Lm, qin%OAMdens%Lm)
      Deallocate (qin%spinmuz%Lm, qin%spinmux%Lm, qin%spinmuy%Lm)
      Deallocate (qin%xcur%Rm, qin%ycur%Rm, qin%zcur%Rm, qin%ts%Rm, qin%ts_av%Rm, qin%OAMdens%Rm)
      Deallocate (qin%spinmuz%Rm, qin%spinmux%Rm, qin%spinmuy%Rm)
      qin%nat = 0
      qin%ldim = 0
      qin%alloc = 0
      qin%ts%alloc = 0
      qin%ts_av%alloc = 0
      qin%OAMdens%alloc = 0
      qin%spinmux%alloc = 0
      qin%spinmuy%alloc = 0
      qin%spinmuz%alloc = 0
   end subroutine free_confdisdata
!-----------------
   subroutine free_confgrid(qin)
      Implicit None
      Type(confgrid), Intent(Inout) :: qin
      if (qin%alloc == 0) return
      Deallocate (qin%xcur%Lm, qin%ycur%Lm, qin%zcur%Lm)
      Deallocate (qin%xcur%Rm, qin%ycur%Rm, qin%zcur%Rm)
      Deallocate (qin%ts%Lm, qin%ts%Rm)
      Deallocate (qin%ts_av%Lm, qin%ts_av%Rm)
      Deallocate (qin%OAMdens%Lm, qin%OAMdens%Rm)
      Deallocate (qin%spinmuz%Lm, qin%spinmuz%Rm)
      Deallocate (qin%spinmux%Lm, qin%spinmux%Rm)
      Deallocate (qin%spinmuy%Lm, qin%spinmuy%Rm)
      Deallocate (qin%coord)
      Deallocate (qin%wt)
      qin%alloc = 0
   end subroutine free_confgrid
!-------------

   subroutine alloc_gridata(cf_dx, cf_dy, cf_dxdy, cl)
      Implicit None
      Type(cell_system), Intent(in) :: cl
      Type(confgrid), Intent(InOut):: cf_dx, cf_dy, cf_dxdy
      Integer :: n, m, l
      if (cf_dx%alloc == 0 .and. cf_dy%alloc == 0) then

         n = cl%info(1)
         m = cl%info(2)
         l = cl%info(3)

         call alloc_cf(cf_dx, n, l, cl%nch)
         call set_cf(cf_dx)

         call alloc_cf(cf_dy, m, l, cl%nch)
         call set_cf(cf_dy)

         call alloc_cf(cf_dxdy, n, m, cl%nch)
         call set_cf(cf_dxdy)

      end if

   Contains
    Subroutine alloc_cf(cf_d, i, j, nch)
      Implicit None
      Type(confgrid), Intent(InOut):: cf_d
      Integer :: i,j, nch
      Allocate (cf_d%xcur%Lm   (nch, i, j))
      Allocate (cf_d%xcur%Rm   (nch, i, j))
      Allocate (cf_d%ycur%Lm   (nch, i, j))
      Allocate (cf_d%ycur%Rm   (nch, i, j))
      Allocate (cf_d%zcur%Lm   (nch, i, j))
      Allocate (cf_d%zcur%Rm   (nch, i, j))
      Allocate (cf_d%ts%Lm     (7, i, j))
      Allocate (cf_d%ts%Rm     (7, i, j))
      Allocate (cf_d%ts_av%Lm  (1, i, j))
      Allocate (cf_d%ts_av%Rm  (1, i, j))
      Allocate (cf_d%OAMdens%Lm(3, i, j))
      Allocate (cf_d%OAMdens%Rm(3, i, j))
      Allocate (cf_d%spinmuz%Lm(2, i, j))
      Allocate (cf_d%spinmuz%Rm(2, i, j))
      Allocate (cf_d%spinmux%Lm(2, i, j))
      Allocate (cf_d%spinmux%Rm(2, i, j))
      Allocate (cf_d%spinmuy%Lm(2, i, j))
      Allocate (cf_d%spinmuy%Rm(2, i, j))
      Allocate (cf_d%coord     (2, i, j))
      Allocate (cf_d%wt           (i, j))
      cf_d%alloc = 1
    End Subroutine alloc_cf

    Subroutine set_cf(cf_d)
      Implicit None
      Type(confgrid), Intent(InOut):: cf_d
      cf_d%xcur%Lm = 0.0d0
      cf_d%ycur%Lm = 0.0d0
      cf_d%zcur%Lm = 0.0d0
      cf_d%xcur%Rm = 0.0d0
      cf_d%ycur%Rm = 0.0d0
      cf_d%zcur%Rm = 0.0d0
      cf_d%ts%Lm = 0.0d0
      cf_d%ts%Rm = 0.0d0
      cf_d%ts_av%Lm = 0.0d0
      cf_d%ts_av%Rm = 0.0d0
      cf_d%OAMdens%Lm = 0.0d0
      cf_d%OAMdens%Rm = 0.0d0
      cf_d%spinmuz%Lm = 0.0d0
      cf_d%spinmuz%Rm = 0.0d0
      cf_d%spinmux%Lm = 0.0d0
      cf_d%spinmux%Rm = 0.0d0
      cf_d%spinmuy%Lm = 0.0d0
      cf_d%spinmuy%Rm = 0.0d0
      cf_d%coord = 0.0d0
      cf_d%wt = 0
    End Subroutine set_cf

   end subroutine alloc_gridata

   subroutine free_confdisav(qin)
      Implicit None
      Type(confdisav), Intent(Inout) :: qin
      if (allocated(qin%xcur%Lm)) then
         Deallocate (qin%xcur%Lm, qin%ycur%Lm, qin%zcur%Lm)
         Deallocate (qin%xcur%Rm, qin%ycur%Rm, qin%zcur%Rm)
         qin%xcur%alloc = 0
         qin%ycur%alloc = 0
         qin%zcur%alloc = 0
      end if
   end subroutine free_confdisav
   subroutine add_gridata(cf1, cf2, factor)
      Implicit None
      Type(confgrid), Intent(InOut) :: cf1
      Type(confgrid), Intent(In) :: cf2
      Double precision, Optional, Intent(In) :: factor
      Double precision :: f
      f = 1.0d0
      If (present(factor)) f = factor
      cf1%xcur%Lm = cf1%xcur%Lm + cf2%xcur%Lm*f
      cf1%xcur%Rm = cf1%xcur%Rm + cf2%xcur%Rm*f
      cf1%ycur%Lm = cf1%ycur%Lm + cf2%ycur%Lm*f
      cf1%ycur%Rm = cf1%ycur%Rm + cf2%ycur%Rm*f
      cf1%zcur%Lm = cf1%zcur%Lm + cf2%zcur%Lm*f
      cf1%zcur%Rm = cf1%zcur%Rm + cf2%zcur%Rm*f
      cf1%ts%Lm = cf1%ts%Lm + cf2%ts%Lm*f
      cf1%ts%Rm = cf1%ts%Rm + cf2%ts%Rm*f
      cf1%ts_av%Lm = cf1%ts_av%Lm + cf2%ts_av%Lm*f
      cf1%ts_av%Rm = cf1%ts_av%Rm + cf2%ts_av%Rm*f
      cf1%OAMdens%Lm = cf1%OAMdens%Lm + cf2%OAMdens%Lm*f
      cf1%OAMdens%Rm = cf1%OAMdens%Rm + cf2%OAMdens%Rm*f
      cf1%spinmuz%Lm = cf1%spinmuz%Lm + cf2%spinmuz%Lm*f
      cf1%spinmuz%Rm = cf1%spinmuz%Rm + cf2%spinmuz%Rm*f
      cf1%spinmux%Lm = cf1%spinmux%Lm + cf2%spinmux%Lm*f
      cf1%spinmux%Rm = cf1%spinmux%Rm + cf2%spinmux%Rm*f
      cf1%spinmuy%Lm = cf1%spinmuy%Lm + cf2%spinmuy%Lm*f
      cf1%spinmuy%Rm = cf1%spinmuy%Rm + cf2%spinmuy%Rm*f
      cf1%coord = cf1%coord + cf2%coord*f
      cf1%wt = cf1%wt + cf2%wt*f
   end subroutine add_gridata
!------------------------------
   subroutine reset_confdis(cf)
      Implicit None
      Type(confdisdata), Intent(InOut) :: cf
      call free_confdisdata(cf)
   end subroutine reset_confdis
   subroutine reset_confdisav(cf)
      Implicit None
      Type(confdisav), Intent(InOut) :: cf
      call free_confdisav(cf)
   end subroutine reset_confdisav
!-------------------------------
   subroutine reset_confgrid(cf)
      Implicit None
      Type(confgrid), Intent(InOut) :: cf
      call free_confgrid(cf)
   end subroutine reset_confgrid
!------------------------------------
   subroutine write_dis(dis, cl, scl, path)
      Implicit None
      Type(confdisdata), Intent(In) :: dis
      Type(cell_system), Intent(In) :: cl
      Real(Kind=DEF_DBL_PREC), Intent(in) ::scl
      Character(*), Intent(In) :: path
      Integer :: i
      open (unit=370, file=trim(path))
      Do i = 1, product(cl%info)
         write (370, 1001) cl%cell(:, i)*scl, dis%zcur%Lm(:, i)
      End Do
      close (370)
1001  format(7g17.8e3)
   end subroutine write_dis

   subroutine writegridcur(cf, path)
      Implicit None
      Type(confgrid), Intent(In) :: cf
      Character(*), Intent(In) :: path
      Integer :: i, j, k
      ! pretty disgusting routine, should be cleaned up... 
      open (unit=371, file=trim(path)//'/Lxcur.dat')
      open (unit=372, file=trim(path)//'/Lycur.dat')
      open (unit=373, file=trim(path)//'/Lzcur.dat')
      open (unit=374, file=trim(path)//'/Rxcur.dat')
      open (unit=375, file=trim(path)//'/Rycur.dat')
      open (unit=376, file=trim(path)//'/Rzcur.dat')
      open (unit=377, file=trim(path)//'/Lmuz.dat')
      open (unit=378, file=trim(path)//'/Rmuz.dat')
      open (unit=379, file=trim(path)//'/Lmux.dat')
      open (unit=380, file=trim(path)//'/Rmux.dat')
      open (unit=381, file=trim(path)//'/Lmuy.dat')
      open (unit=382, file=trim(path)//'/Rmuy.dat')
      open (unit=383, file=trim(path)//'/wts.dat')
      open (unit=384, file=trim(path)//'/Rdens.dat')
      open (unit=385, file=trim(path)//'/Ldens.dat')
      open (unit=386, file=trim(path)//'/ROAMdens.dat')
      open (unit=387, file=trim(path)//'/LOAMdens.dat')
      Do i = 1, size(cf%xcur%Lm, DIM=3)
      Do j = 1, size(cf%xcur%Lm, DIM=2)
      k = size(cf%xcur%Lm, DIM=1)
      If (cf%wt(j, i) /= 0) then
         if (k == 4) then
         write (371, 1002) cf%coord(:, j, i), cf%xcur%Lm(:, j, i)
         write (372, 1002) cf%coord(:, j, i), cf%ycur%Lm(:, j, i)
         write (373, 1002) cf%coord(:, j, i), cf%zcur%Lm(:, j, i)
         write (374, 1002) cf%coord(:, j, i), cf%xcur%Rm(:, j, i)
         write (375, 1002) cf%coord(:, j, i), cf%ycur%Rm(:, j, i)
         write (376, 1002) cf%coord(:, j, i), cf%zcur%Rm(:, j, i)
         else if (k == 7) then
         write (371, 1001) cf%coord(:, j, i), cf%xcur%Lm(:, j, i)
         write (372, 1001) cf%coord(:, j, i), cf%ycur%Lm(:, j, i)
         write (373, 1001) cf%coord(:, j, i), cf%zcur%Lm(:, j, i)
         write (374, 1001) cf%coord(:, j, i), cf%xcur%Rm(:, j, i)
         write (375, 1001) cf%coord(:, j, i), cf%ycur%Rm(:, j, i)
         write (376, 1001) cf%coord(:, j, i), cf%zcur%Rm(:, j, i)
         else 
          STOP 'number of currents is not equal to 4 (c + js) or 7 (c + js + ls)'
         endif
         write (377, 1003) cf%coord(:, j, i), cf%spinmuz%Lm(:, j, i)
         write (378, 1003) cf%coord(:, j, i), cf%spinmuz%Rm(:, j, i)
         write (379, 1003) cf%coord(:, j, i), cf%spinmux%Lm(:, j, i)
         write (380, 1003) cf%coord(:, j, i), cf%spinmux%Rm(:, j, i)
         write (381, 1003) cf%coord(:, j, i), cf%spinmuy%Lm(:, j, i)
         write (382, 1003) cf%coord(:, j, i), cf%spinmuy%Rm(:, j, i)
         write (383, 1004) cf%coord(:, j, i), cf%wt(j, i)
         write (384, 1002) cf%coord(:, j, i), cf%ts%Lm(:, j, i)
         write (385, 1002) cf%coord(:, j, i), cf%ts%Rm(:, j, i)
         write (386, 1005) cf%coord(:, j, i), cf%OAMdens%Lm(:, j, i)
         write (387, 1005) cf%coord(:, j, i), cf%OAMdens%Rm(:, j, i)
      End if
      End Do
      End Do
      close (371)
      close (372)
      close (373)
      close (374)
      close (375)
      close (376)
      close (377)
      close (378)
      close (379)
      close (380)
      close (381)
      close (382)
      close (383)
      close (384)
      close (385)
      close (386)
      close (387)
1001  format(9g17.8e3)
1002  format(9g17.8e3)
1003  format(4g17.8e3)
1004  format(3g17.8e3)
1005  format(5g17.8e3)
   end subroutine writegridcur
!--------------------perfect screening--------------------
   subroutine pscr_grid(cf)
      Implicit None
      Type(confgrid), Intent(InOut) :: cf
      Integer :: i, j, n, m
      n = size(cf%ts_av%Lm, DIM=2)
      m = size(cf%ts_av%Lm, DIM=3)
      Do i = 1, n
      Do j = 1, m
      if (all(cf%zcur%Lm(:, i, j) /= 0.0d0)) then
         call pscr_dis(cf%xcur%Lm(:, i, j), cf%xcur%Rm(:, i, j), cf%ts_av%Lm(1, i, j), cf%ts_av%Rm(1, i, j))
         call pscr_dis(cf%ycur%Lm(:, i, j), cf%ycur%Rm(:, i, j), cf%ts_av%Lm(1, i, j), cf%ts_av%Rm(1, i, j))
         call pscr_dis(cf%zcur%Lm(:, i, j), cf%zcur%Rm(:, i, j), cf%ts_av%Lm(1, i, j), cf%ts_av%Rm(1, i, j))
         call pscr_dis(cf%ts%Lm(:, i, j), cf%ts%Rm(:, i, j), cf%ts_av%Lm(1, i, j), cf%ts_av%Rm(1, i, j))
         call pscr_dis(cf%OAMdens%Lm(:, i, j), cf%OAMdens%Rm(:, i, j), cf%ts_av%Lm(1, i, j), cf%ts_av%Rm(1, i, j))
      end if
      End Do
      End Do
   end subroutine pscr_grid
   subroutine pscr_dis(L, R, nL, nR)
      Implicit None
      Real(Kind=DEF_DBL_PREC), Intent(in) :: nL, nR
      Real(Kind=DEF_DBL_PREC), Intent(inout) :: L(:), R(:)
      Real(Kind=DEF_DBL_PREC), allocatable :: Lscr(:), Rscr(:)
      Integer :: i, dim
      double precision :: muL, muR
      dim = size(L)
      Allocate (Lscr(dim), Rscr(dim))
      Lscr = 0.0
      Rscr = 0.0
      muL = nL/(nL + nR)
      muR = nR/(nL + nR)
      do i = 1, dim
         Lscr(i) = L(i) - muL*(L(i) + R(i))
         Rscr(i) = R(i) - muR*(L(i) + R(i))
      end do
      L = Lscr
      R = Rscr
      Deallocate (Lscr, Rscr)
   end subroutine pscr_dis
!----------------------------------------
!--------------------Coulomb screening--------------------
   subroutine cscr_grid(cf)
      Implicit None
      Type(confgrid), Intent(InOut) :: cf
      Integer :: i, j, n, m
      n = size(cf%ts%Lm, DIM=2)
      m = size(cf%ts%Lm, DIM=3)
      Do i = 1, n
      Do j = 1, m
      if (all(cf%zcur%Lm(:, i, j) /= 0.0d0)) then
         call pscr_dis(cf%xcur%Lm(:, i, j), cf%xcur%Rm(:, i, j), cf%ts%Lm(1, i, j), cf%ts%Rm(1, i, j))
         call pscr_dis(cf%ycur%Lm(:, i, j), cf%ycur%Rm(:, i, j), cf%ts%Lm(1, i, j), cf%ts%Rm(1, i, j))
         call pscr_dis(cf%zcur%Lm(:, i, j), cf%zcur%Rm(:, i, j), cf%ts%Lm(1, i, j), cf%ts%Rm(1, i, j))
      end if
      End Do
      End Do
   end subroutine cscr_grid
   subroutine cscr_dis(L, R, nL, nR)
      Implicit None
      Real(Kind=DEF_DBL_PREC), Intent(in) :: nL, nR
      Real(Kind=DEF_DBL_PREC), Intent(inout) :: L(:), R(:)
      Real(Kind=DEF_DBL_PREC), allocatable :: Lscr(:), Rscr(:)
      Integer :: i, dim
      double precision :: muL, muR
      dim = size(L)
      Allocate (Lscr(dim), Rscr(dim))
      Lscr = 0.0
      Rscr = 0.0
      muL = nL+nR
      muR = nR+nL
      do i = 1, dim
         Lscr(i) = L(i) - muL*(L(i) + R(i))
         Rscr(i) = R(i) - muR*(L(i) + R(i))
      end do
      L = Lscr
      R = Rscr
      Deallocate (Lscr, Rscr)
   end subroutine cscr_dis
!----------------------------------------

!-------------generating grid-------------------
   subroutine grid_dis(cl, dpi)
      Implicit None
      Integer, Intent(In) :: dpi !increase this factor to make mesh denser
!        Real (Kind=DEF_DBL_PREC),Intent(In)  :: base(:,:)
      Real(Kind=DEF_DBL_PREC), allocatable:: zl(:), zr(:), divx(:), divy(:), divz(:)
      Type(cell_system), Intent(InOut) :: cl
      Real(Kind=DEF_DBL_PREC), Allocatable :: xyc(:, :), xf(:), xb(:), yu(:), yd(:)
      Integer :: i, j, l, n, m, k, h
      Allocate (zl(size(cl%zlr(1, :))), zr(size(cl%zlr(2, :))))
!        coord=cl%pts
      zl = cl%zlr(1, :)
      zr = cl%zlr(2, :)
      n = cl%sc(1)*dpi
      m = cl%sc(2)*dpi
      
      l = cl%info(3)
      Allocate (xyc(2, cl%sc(1)*cl%sc(2)))
      Allocate (xf(l), xb(l), yu(l), yd(l))
      Allocate (divx(l), divy(l), divz(l))
      !----experimental----
      Do i = 1, l
         k = 0
         xyc = 0
         Do j = 1, size(cl%tr(1, :))
            if ((cl%tr(3, j) > zl(i)) .and. (cl%tr(3, j) < zr(i))) then
               k = k + 1
               if (k > n*m)  then
                  write(*,*) 'number of atoms in a layer: ',k,' exceeds SC: ', n*m
                  !stop
               endif
               xyc(1:2, k) = cl%tr(1:2, j)
            end if
         End Do
         if (k < n*m) then
            write(*,*) 'number of atoms in a layer: ',k,' is less than SC: ', n*m
            !stop
         end if
         xf(i) = MINVAL(xyc(1, :))
         xb(i) = MAXVAL(xyc(1, :))
         yu(i) = MINVAL(xyc(2, :))
         yd(i) = MAXVAL(xyc(2, :))
         divx(i) = abs(xb(i) - xf(i))/dble(n - 1)
         divy(i) = abs(yd(i) - yu(i))/dble(m - 1)
         divz(i) = abs(zr(i) - zl(i))
         !------adding extra cells on the edges to prevent current leakage
         xf(i) = xf(i) - divx(i)
         xb(i) = xb(i) + divx(i)
         yu(i) = yu(i) - divy(i)
         yd(i) = yd(i) + divy(i)
         !----------------------------------------------------------------
      End Do
      !---increasing n,m to account for the extra cells
      n = n + 2
      m = m + 2
      !------------------
      h = 0
      if (allocated(cl%cell)) Deallocate (cl%cell, cl%dx, cl%dy, cl%dz)
      Allocate (cl%cell(3, n*m*l))
      Allocate (cl%dx(n*m*l), cl%dy(n*m*l), cl%dz(n*m*l))
      Do i = 1, l
         Do j = 1, n
            Do k = 1, m
               h = h + 1
               cl%cell(1, h) = xf(i) + (j - 1)*divx(i)
               cl%cell(2, h) = yu(i) + (k - 1)*divy(i)
               cl%cell(3, h) = 0.5d0*(zl(i) + zr(i))
               cl%dx(h) = divx(i)
               cl%dy(h) = divy(i)
               cl%dz(h) = divz(i)
            End Do
         End Do
      End Do
      cl%info(1) = n
      cl%info(2) = m
      Deallocate (xf, xb, yu, yd, zl, zr)
      Deallocate (divx, divy, divz)
      Deallocate (xyc)
!        cl%dx=divx
!        cl%dy=divy
!        cl%dz=divz
!--------------------------
!        xmin=MINVAL(coord(1,:))
!        xmax=MAXVAL(coord(1,:))
!        ymin=MINVAL(coord(2,:))
!        ymax=MAXVAL(coord(2,:))
!        zmin=MINVAL(zl)
!        zmax=MAXVAL(zr)
!        divx=(xmax-xmin)/dble(n-1)
!        divy=(ymax-ymin)/dble(m-1)
!!        divx=abs(base(1,1))/dble(n-1)
!!        divy=abs(base(1,2))/dble(m-1)
!        divz=(zmax-zmin)/dble(l)
      !adjusting to the global bounding box for non-orthogonal systems
! !       n=NINT((xmax-xmin)/divx)+1
! !      m=NINT((ymax-ymin)/divy)+1
!        Do i=1,l-1
!        zr(i)=zl(i)+ divz
!        zl(i+1)=zr(i)
!        End Do
!        zr(l)=zl(l)+divz
!        Allocate(mesh(n+1,m+1))
!        xmin=xmin-0.50d0*divx
!        ymin=ymin-0.50d0*divy
!        DO i=1,n+1
!        mesh(i,1)=xmin+dble(i-1)*divx !i+2 , to have positive grid indexing
!        END DO
!        DO j=1,m+1
!        mesh(j,2)=ymin+dble(j-1)*divy !j+2 , to have positive grid indexing
!        END DO
!        h=0
!        if (allocated(cl%cell)) then
!        Deallocate(cl%cell)
!        end if
!        Allocate(cl%cell(3,n*m*l))
!        DO i=1,n
!        DO j=1,m
!        DO k=1,l
!        h=h+1
!        cl%cell(1,h)=0.5d0*(mesh(i+1,1)+mesh(i,1))
!        cl%cell(2,h)=0.5d0*(mesh(j+1,2)+mesh(j,2))
!        cl%cell(3,h)=0.5d0*(zl(k)+zr(k))
!        END DO
!        END DO
!        END DO
!        cl%info(1)=n
!        cl%info(2)=m
!        cl%dx=divx
!        cl%dy=divy
!        cl%dz=divz
!        Deallocate(mesh)
   end subroutine grid_dis
!------------------------------------

!----------averaging currents in the grid----------
   subroutine makegridcur(cf_d, cl, cf_dx, cf_dy, cf_dxdy, scl)
      Implicit None
      Type(confdisdata), Intent(In) :: cf_d
      Type(cell_system), Intent(In) :: cl
      Real(Kind=DEF_DBL_PREC), Intent(In) ::scl
      Type(confgrid), Intent(Out) :: cf_dx, cf_dy, cf_dxdy
      Real(Kind=DEF_DBL_PREC), allocatable :: x(:, :), y(:, :), z(:), xtmp(:), ytmp(:)
      Real(Kind=DEF_DBL_PREC), allocatable ::xav(:), yav(:)
      Integer i, j, k, h, n, m, l
      if (cf_dx%alloc /= 0 .and. cf_dy%alloc /= 0) then
         call free_confgrid(cf_dx)
         call free_confgrid(cf_dy)
      end if
      call alloc_gridata(cf_dx, cf_dy, cf_dxdy, cl)
      n = cl%info(1)
      m = cl%info(2)
      l = cl%info(3)
      allocate (x(n, l), y(m, l), z(l), xtmp(n*m), ytmp(n*m))
!----------------------------------------
      Do k = 1, l
         j = 0
         Do h = 1, l*m*n
            if (cl%pl(h) == k) then
               j = j + 1
               xtmp(j) = cl%cell(1, h)
               ytmp(j) = cl%cell(2, h)
            end if
         End Do
         call unique(xtmp, x(:, k), n)
         call unique(ytmp, y(:, k), m)
      End Do
      call unique(cl%cell(3, :), z, l)
!-----------Averaging coordinates ----------------
      Allocate (xav(n), yav(m))
      Do i = 1, n
         xav(i) = sum(x(i, :))/l
      End Do
      Do j = 1, m
         yav(j) = sum(y(j, :))/l
      End Do
! Note from MSR: these are very naive implementations.. there is a huge redundancy in
! these double loops. The 3D grid is generated through code and should be consist of parallelipeds, hence
! it should be possible to find the correct cells for the xz, yz and xy averaging through mathematical 
! operations (i.e. mod and integer division). TODO refactor this
!--------------------------yaverage---------------
      Do k = 1, l
      Do i = 1, n
         cf_dx%coord(2, i, k) = z(k)
         cf_dx%coord(1, i, k) = xav(i)
         Do h = 1, n*m*l
            if (cl%cell(1, h) == x(i, k) .and. cl%cell(3, h) == z(k)) then
               cf_dx%xcur%Lm(:, i, k) = cf_dx%xcur%Lm(:, i, k) + cf_d%xcur%Lm(:, h)
               cf_dx%ycur%Lm(:, i, k) = cf_dx%ycur%Lm(:, i, k) + cf_d%ycur%Lm(:, h)
               cf_dx%zcur%Lm(:, i, k) = cf_dx%zcur%Lm(:, i, k) + cf_d%zcur%Lm(:, h)
               cf_dx%xcur%Rm(:, i, k) = cf_dx%xcur%Rm(:, i, k) + cf_d%xcur%Rm(:, h)
               cf_dx%ycur%Rm(:, i, k) = cf_dx%ycur%Rm(:, i, k) + cf_d%ycur%Rm(:, h)
               cf_dx%zcur%Rm(:, i, k) = cf_dx%zcur%Rm(:, i, k) + cf_d%zcur%Rm(:, h)
               cf_dx%ts%Lm(:, i, k) = cf_dx%ts%Lm(:, i, k) + cf_d%ts%Lm(:, h)
               cf_dx%ts%Rm(:, i, k) = cf_dx%ts%Rm(:, i, k) + cf_d%ts%Rm(:, h)
               cf_dx%ts_av%Lm(:, i, k) = cf_dx%ts_av%Lm(:, i, k) + cf_d%ts_av%Lm(:, h)
               cf_dx%ts_av%Rm(:, i, k) = cf_dx%ts_av%Rm(:, i, k) + cf_d%ts_av%Rm(:, h)
               cf_dx%OAMdens%Lm(:, i, k) = cf_dx%OAMdens%Lm(:, i, k) + cf_d%OAMdens%Lm(:, h)
               cf_dx%OAMdens%Rm(:, i, k) = cf_dx%OAMdens%Rm(:, i, k) + cf_d%OAMdens%Rm(:, h)
               cf_dx%spinmuz%Lm(:, i, k) = cf_dx%spinmuz%Lm(:, i, k) + cf_d%spinmuz%Lm(:, h)
               cf_dx%spinmuz%Rm(:, i, k) = cf_dx%spinmuz%Rm(:, i, k) + cf_d%spinmuz%Rm(:, h)
               cf_dx%spinmux%Lm(:, i, k) = cf_dx%spinmux%Lm(:, i, k) + cf_d%spinmux%Lm(:, h)
               cf_dx%spinmux%Rm(:, i, k) = cf_dx%spinmux%Rm(:, i, k) + cf_d%spinmux%Rm(:, h)
               cf_dx%spinmuy%Lm(:, i, k) = cf_dx%spinmuy%Lm(:, i, k) + cf_d%spinmuy%Lm(:, h)
               cf_dx%spinmuy%Rm(:, i, k) = cf_dx%spinmuy%Rm(:, i, k) + cf_d%spinmuy%Rm(:, h)
               cf_dx%wt(i, k) = cf_dx%wt(i, k) + cl%wt(h)
            End if
         End Do
      End Do
      End Do
!--------------------------------xaverage
      Do k = 1, l
      Do j = 1, m
         cf_dy%coord(2, j, k) = z(k)
         cf_dy%coord(1, j, k) = yav(j)
         Do h = 1, n*m*l
            if (cl%cell(2, h) == y(j, k) .and. cl%cell(3, h) == z(k)) then
               cf_dy%xcur%Lm(:, j, k) = cf_dy%xcur%Lm(:, j, k) + cf_d%xcur%Lm(:, h)
               cf_dy%ycur%Lm(:, j, k) = cf_dy%ycur%Lm(:, j, k) + cf_d%ycur%Lm(:, h)
               cf_dy%zcur%Lm(:, j, k) = cf_dy%zcur%Lm(:, j, k) + cf_d%zcur%Lm(:, h)
               cf_dy%xcur%Rm(:, j, k) = cf_dy%xcur%Rm(:, j, k) + cf_d%xcur%Rm(:, h)
               cf_dy%ycur%Rm(:, j, k) = cf_dy%ycur%Rm(:, j, k) + cf_d%ycur%Rm(:, h)
               cf_dy%zcur%Rm(:, j, k) = cf_dy%zcur%Rm(:, j, k) + cf_d%zcur%Rm(:, h)
               cf_dy%ts%Lm(:, j, k) = cf_dy%ts%Lm(:, j, k) + cf_d%ts%Lm(:, h)
               cf_dy%ts%Rm(:, j, k) = cf_dy%ts%Rm(:, j, k) + cf_d%ts%Rm(:, h)
               cf_dy%ts_av%Lm(:, j, k) = cf_dy%ts_av%Lm(:, j, k) + cf_d%ts_av%Lm(:, h)
               cf_dy%ts_av%Rm(:, j, k) = cf_dy%ts_av%Rm(:, j, k) + cf_d%ts_av%Rm(:, h)
               cf_dy%OAMdens%Lm(:, j, k) = cf_dy%OAMdens%Lm(:, j, k) + cf_d%OAMdens%Lm(:, h)
               cf_dy%OAMdens%Rm(:, j, k) = cf_dy%OAMdens%Rm(:, j, k) + cf_d%OAMdens%Rm(:, h)
               cf_dy%spinmuz%Lm(:, j, k) = cf_dy%spinmuz%Lm(:, j, k) + cf_d%spinmuz%Lm(:, h)
               cf_dy%spinmuz%Rm(:, j, k) = cf_dy%spinmuz%Rm(:, j, k) + cf_d%spinmuz%Rm(:, h)
               cf_dy%spinmux%Lm(:, j, k) = cf_dy%spinmux%Lm(:, j, k) + cf_d%spinmux%Lm(:, h)
               cf_dy%spinmux%Rm(:, j, k) = cf_dy%spinmux%Rm(:, j, k) + cf_d%spinmux%Rm(:, h)
               cf_dy%spinmuy%Lm(:, j, k) = cf_dy%spinmuy%Lm(:, j, k) + cf_d%spinmuy%Lm(:, h)
               cf_dy%spinmuy%Rm(:, j, k) = cf_dy%spinmuy%Rm(:, j, k) + cf_d%spinmuy%Rm(:, h)
               cf_dy%wt(j, k) = cf_dy%wt(j, k) + cl%wt(h)
            end if
         End Do
      End Do
      End Do
!--------------------------------xy plane, zaverage
! MSR: this routine is especially bad, see note above
      Do i = 1, n
      Do j = 1, m
        cf_dxdy%coord(1, i, j) = xav(i)
        cf_dxdy%coord(2, i, j) = yav(j)
      Do k = 1, l
        Do h = 1, n*m*l
          If (cl%cell(1, h) == x(i, k) .and. cl%cell(2, h) == y(j, k)) then
            cf_dxdy%xcur%Lm   (:, i, j) = cf_dxdy%xcur%Lm   (:, i, j) + cf_d%xcur%Lm   (:, h)
            cf_dxdy%ycur%Lm   (:, i, j) = cf_dxdy%ycur%Lm   (:, i, j) + cf_d%ycur%Lm   (:, h)
            cf_dxdy%zcur%Lm   (:, i, j) = cf_dxdy%zcur%Lm   (:, i, j) + cf_d%zcur%Lm   (:, h)
            cf_dxdy%xcur%Rm   (:, i, j) = cf_dxdy%xcur%Rm   (:, i, j) + cf_d%xcur%Rm   (:, h)
            cf_dxdy%ycur%Rm   (:, i, j) = cf_dxdy%ycur%Rm   (:, i, j) + cf_d%ycur%Rm   (:, h)
            cf_dxdy%zcur%Rm   (:, i, j) = cf_dxdy%zcur%Rm   (:, i, j) + cf_d%zcur%Rm   (:, h)
            cf_dxdy%ts%Lm     (:, i, j) = cf_dxdy%ts%Lm     (:, i, j) + cf_d%ts%Lm     (:, h)
            cf_dxdy%ts%Rm     (:, i, j) = cf_dxdy%ts%Rm     (:, i, j) + cf_d%ts%Rm     (:, h)
            cf_dxdy%ts_av%Lm     (:, i, j) = cf_dxdy%ts_av%Lm     (:, i, j) + cf_d%ts_av%Lm     (:, h)
            cf_dxdy%ts_av%Rm     (:, i, j) = cf_dxdy%ts_av%Rm     (:, i, j) + cf_d%ts_av%Rm     (:, h)
            cf_dxdy%OAMdens%Lm     (:, i, j) = cf_dxdy%OAMdens%Lm     (:, i, j) + cf_d%OAMdens%Lm     (:, h)
            cf_dxdy%OAMdens%Rm     (:, i, j) = cf_dxdy%OAMdens%Rm     (:, i, j) + cf_d%OAMdens%Rm     (:, h)
            cf_dxdy%spinmuz%Lm(:, i, j) = cf_dxdy%spinmuz%Lm(:, i, j) + cf_d%spinmuz%Lm(:, h)
            cf_dxdy%spinmuz%Rm(:, i, j) = cf_dxdy%spinmuz%Rm(:, i, j) + cf_d%spinmuz%Rm(:, h)
            cf_dxdy%spinmux%Lm(:, i, j) = cf_dxdy%spinmux%Lm(:, i, j) + cf_d%spinmux%Lm(:, h)
            cf_dxdy%spinmux%Rm(:, i, j) = cf_dxdy%spinmux%Rm(:, i, j) + cf_d%spinmux%Rm(:, h)
            cf_dxdy%spinmuy%Lm(:, i, j) = cf_dxdy%spinmuy%Lm(:, i, j) + cf_d%spinmuy%Lm(:, h)
            cf_dxdy%spinmuy%Rm(:, i, j) = cf_dxdy%spinmuy%Rm(:, i, j) + cf_d%spinmuy%Rm(:, h)
            cf_dxdy%wt           (i, j) = cf_dxdy%wt           (i, j) + cl%wt             (h)
          Endif
        Enddo
      Enddo
      Enddo
      Enddo

!    !--shiftig currents from empty cells
      call fix_cur(cf_dx, cl%info(1), cl%info(3))
      call fix_cur(cf_dy, cl%info(2), cl%info(3))
      call fix_cur(cf_dxdy, cl%info(1), cl%info(2))
      !scaling to nm
      cf_dx%coord = cf_dx%coord*scl
      cf_dy%coord = cf_dy%coord*scl
      cf_dxdy%coord = cf_dxdy%coord*scl
      deallocate (x, y, z, xtmp, ytmp, xav, yav)
   end subroutine makegridcur

   subroutine get_cell(cl, pl, xl, xr, yl, yr, zl, zr)
      Implicit None
      type(planes_system), Intent(In) :: pl
      type(cell_system), Intent(InOut) :: cl
      double precision, allocatable, Intent(Out) :: xl(:), xr(:), yl(:), yr(:), zl(:), zr(:)
!    double precision, allocatable :: zpl(:),zpr(:)
      Integer :: i, ic, n, ip
      n = product(cl%info)
      Allocate (xl(n), xr(n), yl(n), yr(n), zl(n), zr(n))
      Do ic = 1, n
         xr(ic) = cl%cell(1, ic) + 0.5d0*cl%dx(ic)
         xl(ic) = cl%cell(1, ic) - 0.5d0*cl%dx(ic)
         yr(ic) = cl%cell(2, ic) + 0.5d0*cl%dy(ic)
         yl(ic) = cl%cell(2, ic) - 0.5d0*cl%dy(ic)
         zr(ic) = cl%cell(3, ic) + 0.5d0*cl%dz(ic)
         zl(ic) = cl%cell(3, ic) - 0.5d0*cl%dz(ic)
      End Do
!    Allocate(zpl(cl%info(3)),zpr(cl%info(3)))
!    zpl(1)=zl(1)
!    zpr(1)=zr(1)
!    Do i=2,cl%info(3)
!    zpl(i)=zpr(i-1)
!    zpr(i)=zpl(i)+cl%dz(i-1)
!    End Do
      if (allocated(cl%ind) .and. allocated(cl%pl)) then
         Deallocate (cl%ind, cl%pl)
      end if
      Allocate (cl%ind(size(cl%tr(1, :))))
      Allocate (cl%pl(n))
      Do i = 1, size(cl%ind)
         ip = i - pl%l%nat
         if (ip <= 0) then
            cl%ind(i) = 0
         else if (ip > pl%m%nat) then
            cl%ind(i) = n + 100
         else
            Do ic = 1, n
               if ((cl%tr(1, i) >= xl(ic) .and. cl%tr(1, i) <= xr(ic)) &
                   .and. (cl%tr(2, i) >= yl(ic) .and. cl%tr(2, i) <= yr(ic)) &
                   .and. (cl%tr(3, i) >= zl(ic) .and. cl%tr(3, i) <= zr(ic))) then
                  cl%ind(i) = ic
               end if
            End Do
         end if
      End Do
      Do i = 1, cl%info(3)
      Do ic = 1, n
         if (cl%cell(3, ic) >= cl%zlr(1, i) .and. cl%cell(3, ic) <= cl%zlr(2, i)) then
            cl%pl(ic) = i
         end if
      End Do
      End Do
   end subroutine get_cell
!    subroutine inbox(xl,xr,yl,yr,zl,zr,r1,r2,t1,t2)
!        !note! this works only for axis aligned cells
!        Implicit none
!        double precision, Intent(In) :: xl,xr,yl,yr,zl,zr,r1(:),r2(:)
!        double precision, Intent(Out) :: t1,t2
!        double precision ::d(3),t(2),n(3,6),P(6,3,4),nd,np
!        double precision ::ctr(3),r3(3),hit
!        !logical :: inbox
!        Integer :: i,j
!        t1=100
!        t2=100
!        t=100
!        !inbox =.false.
!        ctr=(/ xl+xr, yl+yr, zl+zr/)*0.5d0
!        d=r2-r1
!        If(xl<=min(r1(1),r2(1)).and.xr>=max(r1(1),r2(1)) .and. &
!           yl<=min(r1(2),r2(2)).and.yr>=max(r1(2),r2(2)) .and. &
!           zl<=min(r1(3),r2(3)).and.zr>=max(r1(3),r2(3))) Then
!            t1=0.0d0
!            t2=0.0d0
!                  !inbox=.true.
!        Else If(All(d/=0.0d0)) Then
!            P(1,:,:)=RESHAPE ((/ xl,yl,zl,xl,yl,zr,xl,yr,zl,xl,yr,zr/), (/3, 4/))
!            P(2,:,:)=RESHAPE ((/ xr,yl,zl,xr,yl,zr,xr,yr,zl,xr,yr,zr/), (/3, 4/))
!            P(3,:,:)=RESHAPE ((/ xl,yl,zl,xl,yl,zr,xr,yl,zl,xr,yl,zr/), (/3, 4/))
!            P(4,:,:)=RESHAPE ((/ xl,yr,zl,xl,yr,zr,xr,yr,zl,xr,yr,zr/), (/3, 4/))
!            P(5,:,:)=RESHAPE ((/ xl,yl,zl,xr,yl,zl,xl,yr,zl,xr,yr,zl/), (/3, 4/))
!            P(6,:,:)=RESHAPE ((/ xl,yl,zr,xr,yl,zr,xl,yr,zr,xr,yr,zr/), (/3, 4/))
!            n(:,1)= (/ 1,0,0/)
!            n(:,2)= (/ 1,0,0/)
!            n(:,3)= (/ 0,1,0/)
!            n(:,4)= (/ 0,1,0/)
!            n(:,5)= (/ 0,0,1/)
!            n(:,6)= (/ 0,0,1/)
!            j=0
!            Do i=1,6
!            nd =DOT_PRODUCT(n(:,i),d)
!            if (nd /=0.0d0) then
!            np=DOT_PRODUCT(n(:,i),P(i,:,1)-r1)
!            hit=np/nd
!            if (hit >=0.0d0 .and. hit<=1.0d0 ) then
!            r3 = r1+hit*d
!            if (r3(1)>=MINVAL(P(i,1,:))  .and. r3(1)<= MAXVAL(P(i,1,:)) &
!            .and. r3(2)>=MINVAL(P(i,2,:)) .and. r3(2)<= MAXVAL(P(i,2,:)) &
!            .and. r3(3)>=MINVAL(P(i,3,:)) .and. r3(3)<= MAXVAL(P(i,3,:))) then
!            j=j+1
!            t(j)=hit
!!inbox = .true.
!            end if
!            end if
!            end if
!            End Do
!        End If
!        if(j==1)then
!            if (norm2(ctr-r1) > norm2(ctr-r2))then
!            t1=0.0d0
!            t2=t(j)
!            else
!            t1=t(j)
!            t2=0.0d0
!            end if
!        else if(j==2)then
!            t1=min(t(1),t(2))
!            t2=max(t(1),t(2))
!        end if
!!        End If
!    End subroutine inbox
   subroutine makescr(cf_d, nL, nR, cl)
      Implicit None
      Type(cell_system), Intent(In) :: cl
      Type(confdisdata), Intent(InOut) :: cf_d
      Double precision, Intent(In) :: nL(:), nR(:)
      Integer :: i, j
      Do i = 1, product(cl%info)
         j = cl%pl(i)
         !cf_d%ts_av%Lm(:, i) = nL(j)
         !cf_d%ts_av%Rm(:, i) = nR(j)
         cf_d%ts_av%Lm(1, i) = nL(j)
         cf_d%ts_av%Rm(1, i) = nR(j)
      End Do
!     n=product(cl%info)
!     cf_d%ts%LM=0.0d0
!     cf_d%ts%RM=0.0d0
!!     Do i=1,size(cl%pts(1,:))
!!         j=cl%ind(i)
!!        if( j>0 .and. j<n+100 ) then
!!            cf_d%ts%Lm(:,j)=cf_d%ts%Lm(:,j)+nL(i)
!!            cf_d%ts%Rm(:,j)=cf_d%ts%Rm(:,j)+nR(i)
!!        end if
!!     End Do
   end subroutine makescr
   subroutine inbox(xl, xr, yl, yr, zl, zr, r1, r2, t1, t2)
      Implicit None
      double precision, Intent(In) :: xl, xr, yl, yr, zl, zr, r1(:), r2(:)
      double precision, Intent(Out) :: t1, t2
      double precision ::dx, dy, dz, p, q, r, inf
      integer :: i
      logical :: hit
      hit = .true.
      inf = huge(inf)
      t1 = 0.0d0
      t2 = 1.0d0
      dx = r2(1) - r1(1)
      dy = r2(2) - r1(2)
      dz = r2(3) - r1(3)
      do i = 1, 6
         select case (i)
         case (1)
            p = -dx
            q = -(xl - r1(1))
         case (2)
            p = dx
            q = xr - r1(1)
         case (3)
            p = -dy
            q = -(yl - r1(2))
         case (4)
            p = dy
            q = yr - r1(2)
         case (5)
            p = -dz
            q = -(zl - r1(3))
         case (6)
            p = dz
            q = zr - r1(3)
         end select
         r = q/p
         if (p == 0.0d0 .and. q < 0.0d0) then!no crossing
            hit = .false.
            exit
         else if (p < 0.0d0) then
            if (r > t2) then
               hit = .false.
               exit
            elseif (r >= t1) then
               t1 = r
            end if
         else if (p > 0.0d0) then
            if (r < t1) then
               hit = .false.
               exit
            elseif (r <= t2) then
               t2 = r
            end if
         end if
      end do
      if (.not. (hit)) then
         t1 = -inf
         t2 = inf
      end if
   end subroutine inbox
   subroutine offset_cell(cf_d)
      Implicit None
      Type(confgrid) :: cf_d
      Double Precision :: ofst
      Double Precision, allocatable::ref(:)
      Integer :: i, j, count
      allocate (ref(size(cf_d%coord, DIM=2)))
      Do i = 1, size(cf_d%coord, DIM=3)
         count = 0
         ref = 0.0d0
         Do j = 1, size(cf_d%coord, DIM=2)
            if (All(cf_d%zcur%Lm(:, j, i) /= 0.0d0)) then
               count = count + 1
               ref(count) = cf_d%coord(1, j, i)
            end if
         End Do
         ofst = MINVAL(ref(1:count))
         cf_d%coord(1, :, i) = cf_d%coord(1, :, i) - ofst
      End Do
   end subroutine offset_cell
!!!!!!!!!!! interpolating chemical potential (and other spatial quantities in the future)!!!!!!!!!!!!
   subroutine dis_quant(qout, qin, cl, p)
      Implicit None
      Type(spatialquantdis), Intent(InOut) :: qout
      Type(cell_system), Intent(In) :: cl
      Type(planes_system), Intent(In) :: p
      Type(spadirquant), Intent(In)::qin
!Real (Kind=DEF_DBL_PREC), Intent(In) ::qinl(:,:),qinr(:,:)
      Real(Kind=DEF_DBL_PREC) :: dx, dy, dz
      Integer i, j, k, ip, ndim

      ndim = SIZE(qin%Lir(:,:),1)

      !If (qout%alloc == 1) Then
      !  Deallocate (qout%Lm)
      !  Deallocate (qout%Rm)
      !  qout%alloc = 0
      !Endif
      If (qout%alloc == 0) Then
        Allocate (qout%Lm(ndim, product(cl%info)))
        Allocate (qout%Rm(ndim, product(cl%info)))
        qout%alloc=1
      Endif

      qout%Lm = 0.0d0 !initializing
      qout%Rm = 0.0d0

      Do i = 1, product(cl%info)
      Do j = 1, cl%nm + cl%nl + cl%nr
         dx = abs(cl%cell(1, i) - cl%tr(1, j))
         dy = abs(cl%cell(2, i) - cl%tr(2, j))
         dz = abs(cl%cell(3, i) - cl%tr(3, j))

         if (dx < 0.5*cl%dx(i) .and. dy < 0.5*cl%dy(i) .and. dz < 0.5*cl%dz(i)) then
            if (j < (cl%nl + 1)) then
               ip = p%l%ind(j)
               qout%Lm(:, i) = qout%Lm(:, i) + qin%Lir(:, j)/dble(p%l%napp(ip))
               qout%Rm(:, i) = qout%Rm(:, i) + qin%Rt(:, j)/dble(p%l%napp(ip))
            else if (j > cl%nl .and. j < (cl%nl + cl%nm + 1)) then
               k = j - cl%nl
               ip = p%m%ind(k)
               qout%Lm(:, i) = qout%Lm(:, i) + qin%Lm(:, k)/dble(p%m%napp(ip))
               qout%Rm(:, i) = qout%Rm(:, i) + qin%Rm(:, k)/dble(p%m%napp(ip))
            else if (j > (cl%nl + cl%nm)) then
               k = j - (cl%nl + cl%nm)
               ip = p%r%ind(k)
               qout%Lm(:, i) = qout%Lm(:, i) + qin%Lt(:, k)/dble(p%r%napp(ip))
               qout%Rm(:, i) = qout%Rm(:, i) + qin%Rir(:, k)/dble(p%r%napp(ip))
            end if
         end if
      End Do
      End Do
   end subroutine dis_quant
!!!!!!!!!!!
   subroutine makerotxy(rotx, roty, pl)
      Use sparselib
      Implicit None
      type(planes_system), Intent(In) :: pl
      type(rdensemat), Intent(out) ::rotx, roty
      REAL, PARAMETER :: Pi = 3.1415927
      Allocate (rotx%bl(2, pl%l%nat + pl%m%nat + pl%r%nat), roty%bl(2, pl%l%nat + pl%m%nat + pl%r%nat))
      rotx%bl(1, :) = -Pi/2
      rotx%bl(2, :) = 0.0d0
      roty%bl(1, :) = -Pi/2
      roty%bl(2, :) = -Pi/2
   end subroutine makerotxy
!fixing any atoms on the edge of a cell.
!subroutine fixoverlap(cl)
!Implicit none
!Type(cell_system), Intent(InOut) :: cl
!Real (Kind=DEF_DBL_PREC)::eps,dx,dy,dz
!Integer :: i,j
!eps=0.00000001
!Do i=1,product(cl%info)
!Do j=1,cl%nl+cl%nm+cl%nr
!dx=cl%cell(1,i)-cl%pts(1,j)
!dy=cl%cell(2,i)-cl%pts(2,j)
!dz=cl%cell(3,i)-cl%pts(3,j)
!if(abs(dx)==0.5*cl%dx)then
!cl%pts(1,j)=cl%pts(1,j)+eps*dx
!else if(abs(dy)==0.5*cl%dy)then
!cl%pts(2,j)=cl%pts(2,j)+eps*dy
!else if(abs(dz)==0.5*cl%dz)then
!cl%pts(3,j)=cl%pts(3,j)+eps*dz
!end if
!End Do
!end Do
!end subroutine fixoverlap
!---------------------------------------
   subroutine getclwts(cl)
      Type(cell_system), Intent(InOut) :: cl
      Real(Kind=DEF_DBL_PREC) :: dx, dy, dz
      Integer i, j
      If (allocated(cl%wt)) then
         Deallocate (cl%wt)
      End if
      Allocate (cl%wt(product(cl%info)))
      cl%wt = 0
!call fixoverlap(cl)
      Do i = 1, product(cl%info)
      Do j = 1, cl%nl + cl%nm + cl%nr
         dx = abs(cl%cell(1, i) - cl%tr(1, j))
         dy = abs(cl%cell(2, i) - cl%tr(2, j))
         dz = abs(cl%cell(3, i) - cl%tr(3, j))
         if ((dx < 0.5d0*cl%dx(i)) .and. (dy < 0.5d0*cl%dy(i)) .and. (dz < 0.5d0*cl%dz(i))) then
            cl%wt(i) = cl%wt(i) + 1.0d0
         end if
      End Do
      End Do
   end subroutine getclwts
!-----------------------------
   subroutine unique(vec, vec_unique, n)
! Return only the unique values from vec.
      implicit none
      Real(Kind=DEF_DBL_PREC), intent(in) :: vec(:)
      Real(Kind=DEF_DBL_PREC), intent(out) :: vec_unique(:)
      integer, intent(in):: n
      integer :: i, num, index(n)
      logical, dimension(size(vec)) :: mask

      mask = .false.

      do i = 1, size(vec)

!count the number of occurrences of this element:
         num = count(vec(i) == vec)

         if (num == 1) then
!there is only one, flag it:
            mask(i) = .true.
         else
!flag this value only if it hasn't already been flagged:
            if (.not. any(vec(i) == vec .and. mask)) mask(i) = .true.
         end if

      end do

!return only flagged elements:
!allocate( vec_unique(count(mask)) )
      vec_unique = pack(vec, mask)
      index = 0
      call hpsort_eps_epw(n, vec_unique, index, 0.0d0)
   end subroutine unique
!-----------
!---------------------------------------------------------------------
   subroutine hpsort_eps_epw(n, ra, ind, eps)
!---------------------------------------------------------------------
! sort an array ra(1:n) into ascending order using heapsort algorithm,
! and considering two elements being equal if their values differ
! for less than "eps".
! n is input, ra is replaced on output by its sorted rearrangement.
! create an index table (ind) by making an exchange in the index array
! whenever an exchange is made on the sorted data array (ra).
! in case of equal values in the data array (ra) the values in the
! index array (ind) are used to order the entries.
! if on input ind(1)  = 0 then indices are initialized in the routine,
! if on input ind(1) != 0 then indices are assumed to have been
!                initialized before entering the routine and these
!                indices are carried around during the sorting process
!
! adapted from Numerical Recipes pg. 329 (new edition)
!
      implicit none
!-input/output variables
      integer, intent(in)   :: n
      Real(Kind=DEF_DBL_PREC), intent(in)  :: eps
      integer :: ind(n)
      Real(Kind=DEF_DBL_PREC) :: ra(n)
!-local variables
      integer :: i, ir, j, l, iind
      Real(Kind=DEF_DBL_PREC) :: rra
!
! initialize index array
      IF (ind(1) .eq. 0) then
      DO i = 1, n
         ind(i) = i
      END DO
      END IF
! nothing to order
      IF (n .lt. 2) return
! initialize indices for hiring and retirement-promotion phase
      l = n/2 + 1

      ir = n

      sorting: do

! still in hiring phase
         IF (l .gt. 1) then
            l = l - 1
            rra = ra(l)
            iind = ind(l)
! in retirement-promotion phase.
         ELSE
! clear a space at the end of the array
            rra = ra(ir)
!
            iind = ind(ir)
! retire the top of the heap into it
            ra(ir) = ra(1)
!
            ind(ir) = ind(1)
! decrease the size of the corporation
            ir = ir - 1
! done with the last promotion
            IF (ir .eq. 1) then
! the least competent worker at all !
               ra(1) = rra
!
               ind(1) = iind
               exit sorting
            END IF
         END IF
! whether in hiring or promotion phase, we
         i = l
! set up to place rra in its proper level
         j = l + l
!
         DO while (j .le. ir)
         IF (j .lt. ir) then
! compare to better underling
            IF (hslt(ra(j), ra(j + 1))) then
               j = j + 1
!else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
! this means ra(j) == ra(j+1) within tolerance
!  if (ind (j) .lt.ind (j + 1) ) j = j + 1
            END IF
         END IF
! demote rra
         IF (hslt(rra, ra(j))) then
            ra(i) = ra(j)
            ind(i) = ind(j)
            i = j
            j = j + j
!else if ( .not. hslt ( ra(j) , rra ) ) then
!this means rra == ra(j) within tolerance
! demote rra
! if (iind.lt.ind (j) ) then
!    ra (i) = ra (j)
!    ind (i) = ind (j)
!    i = j
!    j = j + j
! else
! set j to terminate do-while loop
!    j = ir + 1
! endif
! this is the right place for rra
         ELSE
! set j to terminate do-while loop
            j = ir + 1
         END IF
         END DO
         ra(i) = rra
         ind(i) = iind

      END DO sorting
   contains

!  internal function
!  compare two real number and return the result

      logical function hslt(a, b)
         Real(Kind=DEF_DBL_PREC) :: a, b
         IF (abs(a - b) < eps) then
            hslt = .false.
         ELSE
            hslt = (a < b)
         end if
      end function hslt

!
   end subroutine hpsort_eps_epw
!!-----------------
   subroutine fix_cur(cf_d, n1, n2)
      Integer, Intent(In) :: n1, n2
      Type(confgrid), Intent(InOut) :: cf_d
      Integer :: i, k
      Do k = 1, n2
      Do i = 1, n1
      if (cf_d%wt(i, k) == 0) then
         if (i == 1) then
            !moving currents to periodic image cell
            cf_d%xcur%Lm(:, n1 - 1, k) = cf_d%xcur%Lm(:, n1 - 1, k) + cf_d%xcur%Lm(:, i, k)
            cf_d%ycur%Lm(:, n1 - 1, k) = cf_d%ycur%Lm(:, n1 - 1, k) + cf_d%ycur%Lm(:, i, k)
            cf_d%zcur%Lm(:, n1 - 1, k) = cf_d%zcur%Lm(:, n1 - 1, k) + cf_d%zcur%Lm(:, i, k)
            cf_d%xcur%Rm(:, n1 - 1, k) = cf_d%xcur%Rm(:, n1 - 1, k) + cf_d%xcur%Rm(:, i, k)
            cf_d%ycur%Rm(:, n1 - 1, k) = cf_d%ycur%Rm(:, n1 - 1, k) + cf_d%ycur%Rm(:, i, k)
            cf_d%zcur%Rm(:, n1 - 1, k) = cf_d%zcur%Rm(:, n1 - 1, k) + cf_d%zcur%Rm(:, i, k)
            !----------------------setting currents to zero in extra cells
            cf_d%xcur%Lm(:, i, k) = 0.0d0
            cf_d%ycur%Lm(:, i, k) = 0.0d0
            cf_d%zcur%Lm(:, i, k) = 0.0d0
            cf_d%xcur%Rm(:, i, k) = 0.0d0
            cf_d%ycur%Rm(:, i, k) = 0.0d0
            cf_d%zcur%Rm(:, i, k) = 0.0d0
         else if (i == n1) then
            !moving currents to periodic image cell
            cf_d%xcur%Lm(:, 2, k) = cf_d%xcur%Lm(:, 2, k) + cf_d%xcur%Lm(:, i, k)
            cf_d%ycur%Lm(:, 2, k) = cf_d%ycur%Lm(:, 2, k) + cf_d%ycur%Lm(:, i, k)
            cf_d%zcur%Lm(:, 2, k) = cf_d%zcur%Lm(:, 2, k) + cf_d%zcur%Lm(:, i, k)
            cf_d%xcur%Rm(:, 2, k) = cf_d%xcur%Rm(:, 2, k) + cf_d%xcur%Rm(:, i, k)
            cf_d%ycur%Rm(:, 2, k) = cf_d%ycur%Rm(:, 2, k) + cf_d%ycur%Rm(:, i, k)
            cf_d%zcur%Rm(:, 2, k) = cf_d%zcur%Rm(:, 2, k) + cf_d%zcur%Rm(:, i, k)
            !----------------------setting currents to zero in extra cells
            cf_d%xcur%Lm(:, i, k) = 0.0d0
            cf_d%ycur%Lm(:, i, k) = 0.0d0
            cf_d%zcur%Lm(:, i, k) = 0.0d0
            cf_d%xcur%Rm(:, i, k) = 0.0d0
            cf_d%ycur%Rm(:, i, k) = 0.0d0
            cf_d%zcur%Rm(:, i, k) = 0.0d0
         end if
      end if
      End Do
      End Do
   end subroutine fix_cur
end Module discursch
