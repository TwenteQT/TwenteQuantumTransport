#include "math_def.h"
program qfit
Use HDF5io
Use qfitlib

implicit none

Character (Len=100) :: cwork
logical :: log, ignore_l, ignore_r
logical :: do_sh
integer :: iar, i, j
!integer :: statvals(13)
integer :: error
integer :: lvl

integer :: nint
double precision, allocatable :: zint(:)

integer :: n, nl, nr
Double precision :: cond, area, cpera
double precision, allocatable :: curpol(:,:), z(:), magmom(:,:)
double precision, allocatable :: mus(:,:), coldat(:,:)
double precision, allocatable :: spindens(:,:), shang(:,:)

Type(t_parameters), allocatable :: pars(:)

integer :: ns
Integer, allocatable :: nlay(:), layind(:)

Type(t_slabdata), allocatable, target :: slabdata(:)
Type(t_slabdata), Pointer :: sd
Type(t_fitdata), allocatable, target :: colfit(:), shfit(:)
Type(t_dependmat) :: coldep, shdep
   
Integer :: ns_in, sym
Character, allocatable :: materials(:)
Double precision :: qa(3)
Logical :: magflag
Logical :: pu_fm(5,4), pu_nm(5,4), pu_0(5,4), p1_specific(4)

Type(t_fit_options) :: fopt

Integer :: ny

Double precision, parameter :: condquant = 7.7480917346d-5

! All data is fitted to the following form within a layer:
! f = p(1) + p(2)*exp(-(z-zl)/p(4)) + p(3)*exp((z-zr)/p(4)) + p(5)*z
! where zl and zr are the z coordinates of the left and right end of a layer
! Some of the terms above are not neccesary in certain cases.
! This is determined by the arrays with flags below:

ny = 4 ! number of quantities that need to be fitted

! For magnetic layers:
pu_fm(:,1) = (/ .true.,.true.,.true.,.true.,.false. /) ! spin current
pu_fm(:,2) = (/ .true.,.true.,.true.,.true.,.true. /) ! mu up
pu_fm(:,3) = (/ .true.,.true.,.true.,.true.,.true. /) ! mu down
pu_fm(:,4) = (/ .false.,.true.,.true.,.true.,.false. /) ! spin density

! For non-magnetic layers in a structure with a magnetic layer somewhere:
pu_nm(:,1) = (/ .false.,.true.,.true.,.true.,.false. /)
pu_nm(:,2) = (/ .true.,.true.,.true.,.true.,.true. /)
pu_nm(:,3) = (/ .true.,.true.,.true.,.true.,.true. /)
pu_nm(:,4) = (/ .false.,.true.,.true.,.true.,.false. /)

! For non-magnetic layers in a structure without any magnetic layers:
pu_0(:,1) = (/ .false.,.false.,.false.,.false.,.false. /)
pu_0(:,2) = (/ .true.,.false.,.false.,.false.,.true. /)
pu_0(:,3) = (/ .true.,.false.,.false.,.false.,.true. /)
pu_0(:,4) = (/ .false.,.false.,.false.,.false.,.false. /)

! determines wether the constant (p1) is material specific:
! (for spin current yes, chemical potentials no, density not relevant)
p1_specific = (/ .true.,.false.,.false.,.false./)

log = .false.
ignore_l = .false.
ignore_r = .false.
ns_in = -1
sym = 0
do_sh = .false.

iar = iargc()
i = 1
if (iar >= 1) then
   do
      call getarg(i,cwork)
      select case (trim(cwork))
      
         case ('-log')
         log = .true.
         
         case ('-l')
         ignore_l = .true.
         
         case ('-r')
         ignore_r = .true.
         
         case('-m')
         i = i + 1
         call getarg(i,cwork)
         select case (cwork(1:1))
         
            case ('-')
            write(*,'(x,A)') 'No input after -m'
            stop
               
            case default
            ns_in = len_trim(cwork)
            allocate(materials(ns_in))
            Do j = 1, ns_in
               materials(j) = cwork(j:j)
            End Do
    
         end select
         
         case('-e')
         sym = 1
         
         case('-o')
         sym = 2
         
         case('-s')
         i = i + 1
         call getarg(i,cwork)
         select case (trim(cwork))
            case ('independent')
            fopt%lsf_opt = 1
            case ('coupled')
            fopt%lsf_opt = 2
            case ('mu')
            fopt%lsf_opt = 3
            case ('current')
            fopt%lsf_opt = 4
            case ('dens')
            fopt%lsf_opt = 5
            case default
            write(*,'(x,A)') 'Invalid input after option -s'
         end select
         
         case('-n')
         i = i + 1
         call getarg(i,cwork)
         select case (cwork(1:1))
            case ('-')
            write(*,'(x,A)') 'No input after -n'
            stop
               
            case default
            read(cwork,*) fopt%niter
         end select
         
         case('-t')
         i = i + 1
         call getarg(i,cwork)
         select case (cwork(1:1))
            case ('-')
            write(*,'(x,A)') 'No input after -t'
            stop
               
            case default
            read(cwork,*) fopt%tol
         end select
         
         case('-d')
         i = i + 1
         call getarg(i,cwork)
         select case (cwork(1:1))
            case ('-')
            write(*,'(x,A)') 'No input after -d'
            stop
               
            case default
            read(cwork,*) fopt%damp
         end select
         
         case('-v')
         i = i + 1
         call getarg(i,cwork)
         select case (cwork(1:1))
            case ('-')
            write(*,'(x,A)') 'No input after -v'
            stop
               
            case default
            read(cwork,*) fopt%v
         end select
         
         case ('-w')
         fopt%weigh = 1
         i = i + 1
         call getarg(i,cwork)
         select case (cwork(1:1))
            case ('-')
            i = i - 1
               
            case default
            read(cwork,*) fopt%tune
         end select
         
         case ('-sh')
         !secret option
         do_sh = .true.
         
         case ('-help')
         Call print_help_qfit()
         stop
            
         case default
         write(*,'(x,3A)') 'Invalid option "', trim(cwork), '"'
         stop
            
      end select
      if (i == iar) exit
      i = i + 1
   end do
end if

open(200, file="geom_i", action="read", IOSTAT=error)
If (error /= 0) Then
   open(unit=201, file="dirlist", action="read", IOSTAT=error)
   If (error /=0) Then
      write(*,*) "Cannot find geom_i; no dirlist file."
      stop
   End If
   read(201,*) lvl
   If (lvl /= 0) Then
      write(*,*) "Level in dirlist is not 0: working directory does not contain configurations of the same geometry."
      stop
   End If
   read(201,*)
   read(201,*) cwork
   close(201)
   open(unit=200, file=trim(adjustl(cwork))//"/geom_i", action="read", IOSTAT=error)
   If (error /=0) Then
      write(*,*) "Cannot find geom_i file."
      stop
   End If
End If

Read(200,*)
Read(200,*) nint
If (ignore_l) Then
   Read(200,*)
   nint = nint - 1
End If
If (ignore_r) nint = nint - 1
If (nint <=1) Then
   write(*,*) "There must be at least 2 interfaces."
   stop
End If
Allocate(zint(nint))
Do i = 1, nint
   Read(200,*) zint(i)
   zint(i) = zint(i)*0.52917721092e-1
End Do
Close(200)

write(*,*) "Interfaces at the following z-positions (nm):"
write(*,*) zint

ns = nint - 1
If (ns_in /= -1) Then
   If (ns_in /= ns) Then
      write(*,'( "number of slabs = n_int-1 = ", I3)') ns
      write(*,'( "not consistent with user input: ", I3)') ns_in
      stop
   End If
Else
   allocate(materials(ns))
   Do i = 1, ns
      materials(i) = char(i+32)
   End Do
End If
allocate(nlay(ns), slabdata(ns), colfit(ns), pars(ns))

!n = getnumlines("plots/Lzcur.dat")
open(100, file='qfin')
read(100, *) nl, n, nr
read(100, *) cond
read(100, *) area ! in nm
close(100)
cpera = cond*condquant/2.0d0/area ! the conductance quantum is 2*e^2/h, we need e^2/h

Allocate(curpol(3, n), z(n), mus(2, n), spindens(4, n))
Call read_plot("plots/Lzcur.dat", curpol, n, 0, z)
Call read_plot("plots/Lmu.dat", mus, n, nl)
Call read_plot("plots/Ldens_pscr.dat", spindens, n, nl)


Allocate(magmom(3, n))
Call read_plot("plots/magdir.dat", magmom, n, nl)


Allocate(layind(n))
layind(:) = -1
nlay(:) = 0
Do i = 1, n
   Do j = 1, ns
      If (z(i) > zint(j) .and. z(i) < zint(j+1)) Then
         layind(i) = j
         nlay(j) = nlay(j) + 1
         Exit
      End if
   End Do
End Do

magflag = .false.
Do i = 1, ns
   sd => slabdata(i)
   Call alloc_slabdata(sd, nlay(i))
   sd%zl = zint(i)
   sd%zr = zint(i+1)
   sd%zm = (sd%zl+sd%zr)/2.0d0
   Call transfer_data(1, z, sd%z, layind, i)
   Call transfer_data(3, transpose(magmom), sd%magmom, layind, i)
   If (ismagnetic(sd)) magflag = .true.
End Do

Call get_quant_axis(slabdata, qa)
Allocate(coldat(n, 4))
Do i = 1, n
   coldat(i, 1) = dot_product(curpol(:,i), qa)
   coldat(i, 2:3) = mus(1:2, i)
   coldat(i, 4) = dot_product(spindens(1:3,i), qa)
End Do
!coldat(:,2:3) = transpose(mus)

Do i = 1, ns
   sd => slabdata(i)
   If (magflag) Then
      If (sd%magnetic) Then
         Call alloc_fitdata(colfit(i), sd, pu_fm, ny)
         Call Init_pars(pars(i), fopt%lsf_opt, do_sh, pu_fm, ny)
      Else
         Call alloc_fitdata(colfit(i), sd, pu_nm, ny)
         Call Init_pars(pars(i), fopt%lsf_opt, do_sh, pu_nm, ny)
      End If
   Else
      Call alloc_fitdata(colfit(i), sd, pu_0, ny)
      Call Init_pars(pars(i), fopt%lsf_opt, do_sh, pu_0, ny)
   End If
   Call transfer_data(ny, coldat, colfit(i)%y, layind, i)
End Do


Call make_dep(colfit, coldep, materials, p1_specific)
select case (sym)
   case(1) ! delta mu is odd, spin current is even
   Call force_sym_col(colfit, coldep, 1)
   case(2) ! delta mu is even, spin current is odd
   Call force_sym_col(colfit, coldep, -1)
   case default
end select

Do i = 1, ns
   If (fopt%lsf_opt == 2) Then
      Call Couple_pars(coldep, i, 4, 1, i, 4, 2, 1, 1)
      Call Couple_pars(coldep, i, 4, 1, i, 4, 4, 1, 1)
   Else If (fopt%lsf_opt == 3) Then
      Call Couple_pars(coldep, i, 4, 2, i, 4, 1, 1, 0)
      Call Couple_pars(coldep, i, 4, 2, i, 4, 4, 1, 0)
   Else If (fopt%lsf_opt == 4) Then
      Call Couple_pars(coldep, i, 4, 1, i, 4, 2, 1, 0)
      Call Couple_pars(coldep, i, 4, 1, i, 4, 4, 1, 0)
   Else If (fopt%lsf_opt == 5) Then
      Call Couple_pars(coldep, i, 4, 4, i, 4, 1, 1, 0)
      Call Couple_pars(coldep, i, 4, 4, i, 4, 2, 1, 0)
   End If
   
   ! Make sure that the mu_u and mu_d (row 2 and 3 of the data)
   ! have the same offset, l_sf and slope (parameter 1, 4, 5) :
   Call Couple_pars(coldep, i, 1, 2, i, 1, 3, 1, 1)
   Call Couple_pars(coldep, i, 4, 2, i, 4, 3, 1, 1)
   Call Couple_pars(coldep, i, 5, 2, i, 5, 3, 1, 1)
End Do

!write(*,'(45I5)') coldep%irow
!write(*,'(45I5)') coldep%ip
!write(*,'(45I5)') coldep%iy
!write(*,'(45F5.1)') coldep%val

Call direct_search(colfit, coldep, 50)
Call fit_linexp(colfit, coldep, fopt)
Call write_fitdata(colfit, "plots/Ljs_fit.dat", (/1/))
Call write_fitdata(colfit, "plots/Lmu_fit.dat", (/2, 3/))
Call write_fitdata(colfit, "plots/Ldens_fit.dat", (/4/))
Call convert_pars(colfit, coldep)
Call get_vfpars(colfit, coldep, pars, cpera)

If (do_sh) Then
   write(*,*) "Doing SH angle"

   allocate(shang(1,n), shfit(ns))
   Call read_plot("plots/SHang.dat", shang, n, 0)
   Do i = 1, ns
      sd => slabdata(i)
      Call alloc_fitdata(shfit(i), sd, (/ .true.,.true.,.true.,.true.,.false. /), 1)
      Call transfer_data(1, shang, shfit(i)%y, layind, i)
   End Do
   Call make_dep(shfit, shdep, materials, (/.true./))
   Call force_sym_sh(shfit, shdep)
   Call direct_search(shfit, shdep, 50)
   Call fit_linexp(shfit, shdep, fopt)
   Call write_fitdata(shfit, "plots/SHang_fit.dat", (/1/))
   Call convert_pars(shfit, shdep)
   Call get_shpars(shfit, shdep, pars)
End If

Call write_pars(pars, fopt, materials, zint, qa)
write(*,*) 'Finished'

End program qfit