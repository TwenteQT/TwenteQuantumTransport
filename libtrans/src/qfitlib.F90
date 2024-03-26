#include "math_def.h"
Module qfitlib
   Implicit None
   
   Type t_fitdata
      Integer :: ip_ofs, np, n, ny
      Integer, allocatable :: ip(:), iy(:)
      Type(t_slabdata), Pointer :: sd
      Double precision, allocatable :: var(:), mad(:), var_w(:), yw(:)
      Double precision, allocatable :: p(:,:), y(:,:), f(:,:), w(:,:), jac(:,:)
      Double precision :: tune, tol, LMfactor, LMdamp
      Double precision, allocatable :: A(:), B(:), C(:), D(:), l(:)
      Double precision, allocatable :: dA(:), dB(:), dC(:), dD(:), dl(:)
   End Type
   
   Type t_slabdata
      Integer :: n
      Logical :: magnetic
      Double precision :: zl, zr, zm, mdir(3)
      Double Precision, allocatable :: z(:), magmom(:, :)
      Type(t_fitdata) :: pol
   End Type
   
   Type t_dependmat
      Integer :: nr, nc, nz, ny
      Integer, allocatable :: irow(:), is(:), ip(:), iy(:)
      Double Precision, allocatable :: val(:), fac(:)
   End Type
   
   Type t_one_parameter
      Double precision :: val, d
      Logical :: need = .true.
      Character (Len=10) :: units, name
   End Type

   Type t_parameters
      Type(t_one_parameter) :: lsf_j, rho, lsf_mu, lsf_n, lsf_sh, shang, beta
   End Type
   
   Type t_fit_options
      Integer :: weigh=0, niter=1000, lsf_opt=4
      Double precision :: damp = 0.0d0, v = 2.0d0 ! Levenberg-Marquardt parameters
      Double precision :: tune = 4.685, tol = 1.0d-8
   End Type
      
   
#define MAXNPARS 5
    
contains
   
   Subroutine print_help_qfit()
      Implicit None
      write(*,'(x,A)') 'DESCRIPTION:'
      write(*,'(xxx,A)') 'This program fits these quantities:'
      write(*,'(xxx,A)') ' - Spin current'
      write(*,'(xxx,A)') ' - Spin dependent chemical potential'
      write(*,'(xxx,A)') ' - Spin density'
      write(*,'(xxx,A)') 'to the form A + B*exp(-(z-zm)/l) + C*exp((z-zm)/l) + D*z'
      write(*,'(xxx,A)') 'and thereby determines the parameters A, B, C, D and l,'
      write(*,'(xxx,A)') 'where zm is the middle of a layer.'
      write(*,'(xxx,A)') 'Then some useful material parameters can be calculated.'
      write(*,'(xxx,A)') 'The following are implemented:'
      write(*,'(xxx,A)') ' - Beta (the bulk spin polarization)'
      write(*,'(xxx,A)') ' - Resistivity'
      write(*,'(xxx,A)') ' - Spin diffusion length l_sf.'
      write(*,'(xxx,A)') 'The parameters are of course determined for each layer separately,'
      write(*,'(xxx,A)') 'but they will be coupled to each other or put to 0 depending on'
      write(*,'(xxx,A)') 'which materials and symmetries are in the geometry. These things'
      write(*,'(xxx,A)') 'are specified using the options -m, -e, -o, -s'
      write(*,'(xxx,A)') ''
      write(*,'(x,A)') 'OPTIONS:'
      write(*,'(xxx,A)') ' -l        Ignore one layer on the left side'
      write(*,'(xxx,A)') ' -r        Ignore one layer on the right side'
      write(*,'(xxx,A)') ' -m <str>  Specify which layers are of the same material: <str> is'
      write(*,'(xxx,A)') '           a sequence of characters, for example ABA, where each'
      write(*,'(xxx,A)') '           character represents a layer. If 2 layers have the same'
      write(*,'(xxx,A)') '           character, it means they are made from the same'
      write(*,'(xxx,A)') '           material and parameters such as resistivity and l_sf'
      write(*,'(xxx,A)') '           are coupled during the fitting.'
      write(*,'(xxx,A)') ' -e        The magnetization has even symmetry around the middle.'
      write(*,'(xxx,A)') ' -o        The magnetization has odd symmetry around the middle.'
      write(*,'(xxx,A)') ' -s <str>  The spin diffusion length can be fitted to each'
      write(*,'(xxx,A)') '           of the quantities. The possibilities are:'
      write(*,'(xxx,A)') '           <str> = independent'
      write(*,'(xxx,A)') '           <str> = coupled'
      write(*,'(xxx,A)') '           <str> = current'
      write(*,'(xxx,A)') '           <str> = mu'
      write(*,'(xxx,A)') '           <str> = dens'
      write(*,'(xxx,A)') ' -n <num>  Change the maximum number of iterations in the fitting'
      write(*,'(xxx,A)') '           procedure. Should be unnecessary: the default is 1000.'
      write(*,'(xxx,A)') ' -t <num>  Change the tolerance in the fitting. Default is 1e-8.'
      write(*,'(xxx,A)') '           This determines when the fitting is converged.'
      write(*,'(xxx,A)') ' -d <num>  Damp the updating of the parameters in the fitting'
      write(*,'(xxx,A)') '           procedure using the Levenberg-Marquardt (LM) algorithm. The'
      write(*,'(xxx,A)') '           initial damping is <num>. Default is 0.'
      write(*,'(xxx,A)') ' -v <num>  Change the Marquardt factor ''v'' that is used for changing the'
      write(*,'(xxx,A)') '           LM damping parameter ''d'' during the fitting procedure in the'
      write(*,'(xxx,A)') '           following way:'
      write(*,'(xxx,A)') '           If the variance increases, d_new = d*v,'
      write(*,'(xxx,A)') '           If the variance decreases, d_new = d/v.'
      write(*,'(xxx,A)') '           This is necessary to get good convergence.'
      write(*,'(xxx,A)') '           The default is v = 2.0 .'
      write(*,'(xxx,A)') ' -w <num>  Use robust bisquare weights to exclude outliers from'
      write(*,'(xxx,A)') '           the data. The argument is optional and determines'
      write(*,'(xxx,A)') '           the tuning parameter. The default is 4.685. If it is'
      write(*,'(xxx,A)') '           increased, less points are considered to be outliers.'
      
   End Subroutine
   
   Function getnumlines(name) Result(n)
      Implicit None
      Integer :: n
      Character(*) :: name
      Integer :: err
      
      n = 0
      open(200, file=name)
      Do
         read(200, *, IOSTAT=err)
         if (err/=0) exit
         n = n + 1
      End Do
      close(200)
   End Function
   
   Subroutine read_plot(name, arr, n, nskip, z)
      Implicit None
      Character(*) :: name
      Double precision :: arr(:,:)
      Integer :: n, nskip
      Double precision, optional :: z(:)
      
      Integer :: i
      Double precision :: dum
      
      open(200, file=name)
      Do i = 1, nskip
         read(200, *)
      End Do
      Do i = 1, n
         read(200, *) dum, arr(:, i)
         If (present(z)) z(i) = dum
      end Do
   End Subroutine
   
   Subroutine write_fitdata(fd, name, iy)
      Implicit None
      Type(t_fitdata), Target :: fd(:)
      Character(*) :: name
      Integer :: iy(:)
      
      Integer :: i, ns, j, k, ny
      Type(t_fitdata), Pointer :: fdp
      
      ny = size(iy)
      
      Open(200, file=name, action='write')
      ns = size(fd)
      Do i = 1, ns
         fdp => fd(i)
         write(200, '(g17.8)', advance='no') fdp%sd%zl
         Do j = 1, ny
            write(200, '(g17.8)', advance='no') f_linexp_2(fdp%p(:,iy(j)), fdp%sd%zl, fdp%sd%zr, fdp%sd%zl)
         End Do
         write(200,*)
         
         Do j = 1, fdp%n
            write(200, '(g17.8)', advance='no') fdp%sd%z(j)
            Do k = 1, ny
               write(200, '(g17.8)', advance='no') fdp%f(j, iy(k))
            End Do
            write(200,*)
         End Do
         
         write(200, '(g17.8)', advance='no') fdp%sd%zr
         Do j = 1, ny
            write(200, '(g17.8)', advance='no') f_linexp_2(fdp%p(:,iy(j)), fdp%sd%zl, fdp%sd%zr, fdp%sd%zr)
         End Do
         write(200,*)
      End Do
      Close(200)
   End Subroutine
   
   
   
   Subroutine transfer_data(ld, ar_in, ar_out, layind, slab)
      Implicit None
      double precision :: ar_in(*), ar_out(*)
      Integer :: layind(:), slab, ld
      
      Integer i, j, l, n, ofs
      
      n = size(layind)
      j=1
      ofs = 0
      Do l = 1, ld
         Do i = 1, n
            If (layind(i) == slab) Then
               ar_out(j) = ar_in(ofs + i)
               j = j + 1
            End If
         End Do
         ofs = ofs + n
      End Do
   End Subroutine
   
   
   Subroutine Init_pars(pars, lsf_opt, do_sh, par_use, ny)
      Implicit None
      Type(t_parameters) :: pars
      Integer :: lsf_opt, ny
      Logical :: par_use(MAXNPARS, ny)
      Logical :: do_sh
      
      pars%rho%units = "uOhm.cm"
      pars%lsf_j%units = "nm"
      pars%lsf_mu%units = "nm"
      pars%lsf_n%units = "nm"
      pars%lsf_sh%units = "nm"
      pars%beta%units = ""
      pars%SHang%units = ""
      
      pars%rho%name = "rho"
      pars%lsf_j%name = "lsf_j"
      pars%lsf_mu%name = "lsf_mu"
      pars%lsf_n%name = "lsf_n"
      pars%lsf_sh%name = "l_SH"
      pars%beta%name = "beta"
      pars%SHang%name = "SH-angle"
      
      If (lsf_opt == 3) Then
         pars%lsf_j%need = .false.
         pars%lsf_n%need = .false.
      Else If (lsf_opt == 4) Then
         pars%lsf_mu%need = .false.
         pars%lsf_n%need = .false.
      Else If (lsf_opt == 5) Then
         pars%lsf_j%need = .false.
         pars%lsf_mu%need = .false.
      End If
      If (.not. par_use(4,1)) pars%lsf_j%need = .false.
      If (.not. par_use(4,2)) pars%lsf_mu%need = .false.
      If (.not. par_use(4,4)) pars%lsf_n%need = .false.
      pars%beta%need = par_use(1,1)
      pars%shang%need = do_sh
      pars%lsf_sh%need = do_sh
   End Subroutine Init_pars
   
   
   Subroutine alloc_fitdata(fd, sd, par_use, ny)
      Implicit None
      Type(t_fitdata) :: fd
      Type(t_slabdata), Target :: sd
      Logical :: par_use(MAXNPARS, ny)
      Integer :: ny
      
      Integer :: i, np, iy
      
      allocate(fd%ip(MAXNPARS*ny), fd%p(MAXNPARS, ny), fd%iy(MAXNPARS*ny))
      allocate(fd%mad(ny), fd%var(ny), fd%var_w(ny), fd%yw(ny))
      allocate(fd%A(ny), fd%B(ny), fd%C(ny), fd%D(ny), fd%l(ny))
      allocate(fd%dA(ny), fd%dB(ny), fd%dC(ny), fd%dD(ny), fd%dl(ny))
      fd%ny = ny
      np = 0
      Do iy = 1, ny
         Do i = 1, MAXNPARS
            If (par_use(i, iy)) Then
               np = np + 1
               fd%ip(np) = i
               fd%iy(np) = iy
            End If
         End Do
      End Do
      fd%np = np
      fd%sd => sd
      allocate(fd%y(sd%n, ny), fd%f(sd%n, ny), fd%w(sd%n, ny))
      allocate(fd%jac(np, sd%n))
      fd%n = sd%n
      fd%p = 0.0d0
      fd%w = 1.0d0
      fd%yw = 1.0d0
      fd%A = 0.0d0
      fd%B = 0.0d0
      fd%C = 0.0d0
      fd%D = 0.0d0
      fd%l = 0.0d0
      fd%dA = 0.0d0
      fd%dB = 0.0d0
      fd%dC = 0.0d0
      fd%dD = 0.0d0
      fd%dl = 0.0d0
   End Subroutine
            
   Subroutine alloc_slabdata(sd, n)
      Implicit None
      Integer :: n
      Type(t_slabdata) :: sd
      
      allocate(sd%z(n), sd%magmom(n, 3))
      sd%n = n
   End Subroutine
   
   
   Function ismagnetic(sd) Result(flag)
      Implicit None
      Type(t_slabdata), Intent(inout) :: sd
      Logical :: flag
      
      Integer :: n, j
      Double precision :: m(3), a
      
      n = sd%n
      m = 0.0d0
      Do j = 1, n
         m = m + sd%magmom(j, :)
      End Do
      m = m/dble(n)
      a = sqrt(sum(m**2))
      flag = (a > 1.0e-12)
      sd%magnetic = flag
      If (flag) Then
         sd%mdir(:) = m(:)/a
      Else
         sd%mdir = 0.0d0
      End If
   End Function
   
   
   Subroutine get_quant_axis(sd, qa)
      Implicit None
      Type(t_slabdata) :: sd(:)
      Double precision :: qa(3)
      
      Integer :: is, ns, nmag, idx(1)
      Double precision :: s
      
      ns = size(sd)
      nmag = 0
      Do is = 1, ns
         If (sd(is)%magnetic) Then
            nmag = nmag + 1
            idx = maxloc(abs(sd(is)%mdir))
         End If
      End Do
      qa = 0.0d0
      If (nmag > 0) Then
         Do is = 1, ns
            If (sd(is)%magnetic) Then
               s = sign(dble(sd(is)%n), sd(is)%mdir(idx(1)))
               qa = qa + sd(is)%mdir*s
            End If
         End Do
         qa = qa/sqrt(sum(qa**2))
         Do is = 1, ns
            If (sd(is)%magnetic) Then
               sd(is)%mdir = sign(qa, sd(is)%mdir(idx(1)))
            End If
         End Do
      Else
         qa(3) = 1.0d0
      End if
      write(*,'("Quant. axis: ", 3F8.3)') qa
   End Subroutine
      
   
   
   Subroutine make_dep(fd, dep, ar, p1_spec)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Character, Intent(in) :: ar(:)
      Logical :: p1_spec(:)
      
      Integer :: ns, np, i, j, k, ofs, nz
      
      ns = size(fd)
      np = 0
      nz = 0
      Do i = 1, ns
         np=np+fd(i)%np
         nz = nz + fd(i)%n
      End Do
      dep%nz = nz
      allocate(dep%irow(np), dep%val(np), dep%is(np), dep%ip(np), dep%iy(np), dep%fac(np))
      dep%val = 1.0d0
      dep%fac = 1.0d0
      dep%nc = np
      dep%nr = np
      dep%ny = fd(1)%ny
      ofs=0
      Do i = 1, ns
         Do j = 1, fd(i)%np
            k = j + ofs
            dep%is(k) = i
            dep%ip(k) = fd(i)%ip(j)
            dep%iy(k) = fd(i)%iy(j)
            dep%irow(k) = k
         End Do
         ofs = ofs + fd(i)%np
      End Do
      
      Do i = 1, ns
         Do j = 1, i-1
            If (ar(j) == ar(i)) Then
               Do k = 1, fd(i)%ny
                  Call Couple_pars(dep, j, 4, k, i, 4, k, 1, 1)
                  Call Couple_pars(dep, j, 5, k, i, 5, k, 1, 1)
                  If (p1_spec(k)) Call Couple_pars(dep, j, 1, k, i, 1, k, 1, 1)
               End Do
               exit
            End If
         End Do
      End Do
   End Subroutine
   
   Subroutine force_sym_col(fd, dep, fact)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Integer :: fact
      
      integer :: ns, is, is2
      
      ns = size(fd)
      If (mod(ns, 2) /= 1) Then
         write(*,*) 'Cannot enforce symmetry. Number of slabs must be odd.'
         return
      End If
      
      Do is = 1, ns
         is2 = ns+1-is
         Call Couple_pars(dep, is, 2, 1, is2, 3, 1, fact, 1)
         Call Couple_pars(dep, is, 2, 2, is2, 3, 2, -fact, 1)
         Call Couple_pars(dep, is, 2, 3, is2, 3, 3, -fact, 1)
         Call Couple_pars(dep, is, 2, 4, is2, 3, 4, -fact, 1)
      End Do
   End Subroutine
   
   Subroutine force_sym_sh(fd, dep)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      
      integer :: ns, is, is2
      
      ns = size(fd)
      If (mod(ns, 2) /= 1) Then
         write(*,*) 'Cannot enforce symmetry. Number of slabs must be odd.'
         return
      End If
      
      Do is = 1, ns
         is2 = ns+1-is
         Call Couple_pars(dep, is, 2, 1, is2, 3, 1, 1, 1)
      End Do
   End Subroutine
   
   Subroutine couple_pars(dep, is1, ip1, iy1, is2, ip2, iy2, factor, w)
      ! Parameter 2 is coupled to parameter 1
      ! This means its value is <factor> times that of parameter 1
      ! Only put values of -1, 0, 1 for <factor>, other values will cause problems.
      ! A parameter is specified by its slab index is, parameter index ip and data column iy to which it is relevant
      ! <w> determines the weight of parameter 2 in the fitting:
      ! if w=0, the value of parameter 2 is determined by fitting parameter 1 to the data in column iy1
      ! if w=1, the value is determined by fitting to all columns of data that are relevant
      Implicit None
      Type(t_dependmat) :: dep
      Integer :: factor, w
      Integer :: is1, is2, ip1, ip2, iy1, iy2
      
      Integer :: i, j, k, l, ik, il, m
      Integer :: imin, jmin
      Double precision :: f1, f2, v1, v2
      
      If (abs(factor) > 1) Then
         write(*,*) 'Invalid call of couple_pars, factor should be one of -1, 0, 1'
         stop
      End If
      
      Do i = 1, dep%nc
         If (dep%ip(i) == ip1 .and. dep%is(i) == is1 .and. dep%iy(i) == iy1) Then
            imin = i
            f1 = dep%val(i)
            Do m = 1, dep%nc
               If(dep%irow(m) == dep%irow(i)) imin = min(imin, m)
            End Do
            f1 = f1*dep%val(imin)
            Do j = 1, dep%nc
               If (dep%ip(j) == ip2 .and. dep%is(j) == is2 .and. dep%iy(j) == iy2) Then
                  jmin = j
                  f2 = dep%val(j)
                  Do m = 1, dep%nc
                     If(dep%irow(m) == dep%irow(j)) jmin = min(jmin, m)
                  End Do
                  f2 = f2*dep%val(jmin)
                  If (jmin == imin) return
                  k = min(imin, jmin)
                  l = max(imin, jmin)
                  ik = dep%irow(k)
                  il = dep%irow(l)
                  
                  ik = dep%irow(i)
                  il = dep%irow(j)
                  v1 = dep%val(i)
                  v2 = dep%val(j)
                  f1 = dep%fac(i)
                  f2 = dep%fac(j)
                  
                  Do m = 1, dep%nc
                     If (dep%irow(m) == il) Then
                        dep%irow(m) = ik
                        dep%val(m) = dble(factor*w)*dep%val(m)*v1*v2
                        dep%fac(m) = dble(factor)*dep%fac(m)*f1*f2
                     End If
                  End Do
                  Do m = 1, dep%nc
                     If (dep%irow(m) > il) Then
                        dep%irow(m) = dep%irow(m) - 1
                     End If
                  End Do
                  dep%nr = dep%nr - 1
                  return ! DONE HERE
               End If
            End Do
         End If
      End Do
   End Subroutine
   
   
   
   Subroutine solve_Ab(A,b,x,inv)
      ! For symmetric, positive definite matrices
      Implicit None
      Double precision, intent(in) :: A(:,:), b(:)
      Double precision, intent(out) :: x(:)
      Double precision, optional, intent(out) :: inv(:,:)
      
      Double precision, allocatable :: work(:), adum(:,:)
      Double precision :: work1(1)
      Integer :: n, info, lwork, j, k
      Integer, allocatable :: ipiv(:)
      
      n = size(b)
      If (n /= 0) Then
         If (present(inv)) Then
            inv(:,:) = A(:,:)
            Call getinv_A(inv)
            x = 0.0d0
            Do j = 1, n
               Do k = 1, n
                  x(j) = x(j) + inv(j,k)*b(k)
               End Do
            End Do
         Else
            allocate(ipiv(n), adum(n,n))
            call dsysv('L', n, 1, adum, n, ipiv, x, n, work1, -1, info)
            lwork = work1(1)
            adum(:,:) = A(:,:)
            x(:) = b(:)
            allocate(work(lwork))
            call dsysv('L', n, 1, adum, n, ipiv, x, n, work, lwork, info)
            deallocate(ipiv, adum, work)
            If(info/=0) Then
               write(*,*) 'Error in solver dsysv.'
               stop
            End If
         End If
      End If
   End Subroutine
   
   Subroutine getinv_A(A)
      Implicit None
      External dgetrf, dgetri

      Double precision, intent(inout) :: A(:,:)
      
      Integer :: n, info, lwork
      Double precision, allocatable :: ipiv(:), work(:)
      Double precision :: work1(1)
      
      n = size(A,1)
      allocate(ipiv(n))
      
      Call dgetrf(n, n, A, n, ipiv, info)
      If(info/=0) Then
         write(*,*) 'Error in solver detrf.'
         stop
      End If
      
      Call dgetri(n, A, n, ipiv, work1, -1, info)
      lwork = work1(1)
      allocate(work(lwork))
      
      Call dgetri(n, A, n, ipiv, work, lwork, info)
      If(info/=0) Then
         write(*,*) 'Error in solver dgetri.'
         stop
      End If
      
      deallocate(ipiv, work)
   End Subroutine
   
   
   
   
   Subroutine Ab_system(fd, dep, w_flag, A, b, damp)
      Implicit None
      Type(t_fitdata), Target :: fd(:)
      Type(t_dependmat) :: dep
      Double precision :: A(dep%nr, dep%nr)
      Double precision, Optional :: b(dep%nr)
      Double precision, Optional :: damp
      Integer :: w_flag
      
      Integer :: is, ofs, i, j, iz, iy
      Double precision :: jac(dep%nr, dep%nz, dep%ny)
      Double precision :: w(dep%nz, dep%ny), d(dep%nz, dep%ny)
      
      Call jac_system(fd, dep, jac)
      If (w_flag == 1) Then
         ofs = 0
         Do is = 1, size(fd)
            Do iz = 1, fd(is)%n
               Do iy = 1, fd(is)%ny
                  w(ofs+iz, iy) = fd(is)%w(iz, iy)*fd(is)%yw(iy)
                  d(ofs+iz, iy) = fd(is)%y(iz, iy) - fd(is)%f(iz, iy)
               End Do
            End Do
            ofs = ofs + fd(is)%n
         End Do
      Else If (w_flag == 0) Then
         ofs = 0
         Do is = 1, size(fd)
            Do iz = 1, fd(is)%n
               Do iy = 1, fd(is)%ny
                  w(ofs+iz, iy) = fd(is)%yw(iy)
                  d(ofs+iz, iy) = fd(is)%y(iz, iy) - fd(is)%f(iz, iy)
               End Do
            End Do
            ofs = ofs + fd(is)%n
         End Do
      Else
         ofs = 0
         Do is = 1, size(fd)
            Do iz = 1, fd(is)%n
               Do iy = 1, fd(is)%ny
                  w(ofs+iz, iy) = 1.0d0
                  d(ofs+iz, iy) = fd(is)%y(iz, iy) - fd(is)%f(iz, iy)
               End Do
            End Do
            ofs = ofs + fd(is)%n
         End Do
      End If
      A = 0.0d0
      Do iz = 1, dep%nz
         Do iy = 1, dep%ny
            Do i = 1, dep%nr
               Do j = 1, dep%nr
                  A(i, j) = A(i, j) + jac(i, iz, iy)*jac(j, iz, iy)*w(iz, iy)
               End Do
            End Do
         End Do
      End Do
      If (present(b)) Then
         b = 0.0d0
         Do iz = 1, dep%nz
            Do iy = 1, dep%ny
               b(:) = b(:) + jac(:, iz, iy)*w(iz, iy)*d(iz, iy)
            End Do
         End Do
      End If
      If (present(damp)) Then
         Do i = 1, dep%nr
            A(i,i) = A(i,i)*(1.0d0+damp)
         End Do
      End If
      
   End Subroutine
   
   
   
   Subroutine jac_system(fd, dep, jac)
      Implicit None
      Type(t_fitdata), Target :: fd(:)
      Type(t_dependmat) :: dep
      Double precision :: jac(dep%nr, dep%nz, dep%ny)
      
      Integer :: ns, is, ofs, i, iz, ofsz, ip, iy
      Type(t_fitdata), Pointer :: fdp
      Double Precision :: fac, A(dep%nr, dep%nr)
      
      ns = size(fd)
      jac = 0.0d0
      ofs = 0
      ofsz = 0
      Do is = 1, ns
         fdp => fd(is)
         Call f_linexp_sl(fdp)
         Call jac_linexp_sl(fdp)
         Do iz = 1, fdp%n
            Do ip = 1, fdp%np
               i = dep%irow(ip + ofs)
               iy = dep%iy(ip + ofs)
               fac = dep%val(ip + ofs)
               jac(i, iz+ofsz, iy) = jac(i, iz+ofsz, iy) + fdp%jac(ip, iz) *fac
            End Do
         End Do
         ofsz = ofsz + fdp%n
         ofs = ofs + fdp%np
      End Do
      A = 0.0d0
      Do iy = 1, dep%ny
         Do iz = 1, dep%nz
            Do ip = 1, dep%nr
               Do i = 1, dep%nr
                  A(ip, i) = A(ip, i) + jac(ip, iz, iy)*jac(i, iz, iy)
               End Do
            End Do
         End Do
      End Do
   End Subroutine
   
   
   
   
   Subroutine f_linexp_sl(fd)
      Implicit None
      Type(t_fitdata) :: fd
      
      Integer :: i, iy
      
      Do iy = 1, fd%ny
         Do i = 1, fd%n
            fd%f(i, iy) = f_linexp_2(fd%p(:, iy), fd%sd%zl, fd%sd%zr, fd%sd%z(i))
         End Do
      End Do
   End Subroutine
   
   Subroutine jac_linexp_sl(fd)
      Implicit None
      Type(t_fitdata) :: fd
      
      Integer :: ip, i, iz, iy
      Double Precision :: jac(MAXNPARS, fd%ny)
      
      Do iz = 1, fd%n
         Do iy = 1, fd%ny
            jac(:, iy) = jac_linexp_2(fd%p(:, iy), fd%sd%zl, fd%sd%zr, fd%sd%z(iz))
         End Do
         Do i = 1, fd%np
            ip = fd%ip(i)
            iy = fd%iy(i)
            fd%jac(i, iz) = jac(ip, iy)
         End Do
      End Do
   End Subroutine
   
   
   
   Function f_linexp_1(p, zl, zr, z) Result(f)
      Implicit None
      double precision :: p(MAXNPARS), zl, zr, z, f, f1, f2
      
      If (z == zl) Then
         f1 = 1.0d0
         f2 = exp((z-zr)/p(4))
      Else if (z == zr) Then
         f1 = exp(-(z-zl)/p(4))
         f2 = 1.0d0
      Else
         f1 = exp(-(z-zl)/p(4))
         f2 = exp((z-zr)/p(4))
      End If
      f = p(1) + p(2)*f1 + p(3)*f2 + p(5)*z
   End Function
   
   Function jac_linexp_1(p, zl, zr, z) Result(jac)
      Implicit None
      double precision :: p(MAXNPARS), zl, zr, z, jac(MAXNPARS)
      
      jac(1) = 1.0d0
      jac(2) = exp(-(z-zl)/p(4))
      jac(3) = exp((z-zr)/p(4))
      jac(4) = p(2)/p(4)**2*(z-zl)*jac(2) - p(3)/p(4)**2*(z-zr)*jac(3)
      jac(5) = z
   End Function
   
   
   
   
   Function f_linexp_2(p, zl, zr, z) Result(f)
      Implicit None
      double precision :: p(MAXNPARS), zl, zr, z, f, f1, f2, zm
      
      
      If (p(4) == 0.0d0) Then
         f1 = 0.0d0
         f2 = 0.0d0
      Else
         zm = (zr+zl)/2.0d0
         f1 = exp(-(z-zm)/p(4))
         f2 = exp((z-zm)/p(4))
         If (p(1) /= 0.0d0) Then
            f1 = f1 - 1.0d0
            f2 = f2 - 1.0d0
         End If
      End If
      f = p(1) + p(2)*f1 + p(3)*f2 + p(5)*z
   End Function
   
   Function jac_linexp_2(p, zl, zr, z) Result(jac)
      Implicit None
      double precision :: p(MAXNPARS), zl, zr, z, jac(MAXNPARS)
      
      Double precision :: zm, f1, f2
      
      zm = (zr+zl)/2.0d0
      f1 = exp(-(z-zm)/p(4))
      f2 = exp((z-zm)/p(4))
      
      jac(1) = 1.0d0
      If (p(1) == 0.0d0) Then
         jac(2) = f1
         jac(3) = f2
      Else
         jac(2) = f1-1.0d0
         jac(3) = f2-1.0d0
      End If
      jac(4) = (p(2)*f1*(z-zm) - p(3)*f2*(z-zm)) / p(4)**2
      jac(5) = z
   End Function
   
   
   
   Function f_linexp_3(p, zl, zr, z) Result(f)
      Implicit None
      double precision :: p(MAXNPARS), zl, zr, z, f, f1, f2, f3
      
      
      If (p(4) == 0.0d0) Then
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
      Else
         f3 = exp((zl-zr)/p(4))
         f1 = exp(-(z-zl)/p(4))
         f2 = exp((z-zr)/p(4))
      End If
      f = p(1) + p(2)*(f1-f3*f2) + p(3)*(f2-f3*f1) + p(5)*z
   End Function
   
   Function jac_linexp_3(p, zl, zr, z) Result(jac)
      Implicit None
      double precision :: p(MAXNPARS), zl, zr, z, jac(MAXNPARS)
      
      Double precision :: f1, f2, f3
      
      f3 = exp((zl-zr)/p(4))
      f1 = exp(-(z-zl)/p(4))
      f2 = exp((z-zr)/p(4))
      
      jac(1) = 1.0d0
      jac(2) = f1-f3*f2
      jac(3) = f2-f3*f1
      jac(4) = ( p(2)*(f1*(z-zl)+f3*f2*(z+zl-zr-zr) ) - p(3)*(f2*(z-zr)+f3*f1*(z+zr-zl-zl) ) ) / p(4)**2
      jac(5) = z
   End Function
   

   
   Subroutine var_mad(fd)
      Implicit None
      Type(t_fitdata), Target :: fd
      
      Integer :: i, iy
      Double precision :: res, mad, var, var_w, w_tot
      
      Do iy = 1, fd%ny
         var = 0.0d0
         mad = 0.0d0 ! median absolute deviation
         var_w = 0.0d0 ! weighted variance
         w_tot = 0.0d0 ! total of weights, used to normalize weighted variance
         Do i = 1, fd%n
            res = fd%y(i, iy)-fd%f(i, iy)
            mad = mad + abs(res)
            var = var + res**2
            var_w = var_w + res**2*fd%w(i, iy)
            w_tot = w_tot + fd%w(i, iy)
         End Do
         fd%mad(iy) = mad/dble(fd%n)
         fd%var(iy) = var/dble(fd%n)
         fd%var_w(iy) = var_w/w_tot
      End Do
   End Subroutine
   
   Subroutine fit_linexp(fd, dep, fopt)
      Implicit None
      Type(t_fitdata), Target, intent(inout) :: fd(:)
      Type(t_dependmat), intent(in) :: dep
      Type(t_fit_options), intent(in) :: fopt
      !Double precision, Optional :: damp
      
      Integer :: i, j
      Double precision :: A(dep%nr, dep%nr), b(dep%nr), dp(dep%nr)
      !Double precision, parameter :: tol = 1.0d-9, v = 1.1
      Double precision :: var, pvar, d
      Logical :: begindamp
      
      d = fopt%damp
      begindamp = .true.
      Call w_ydependent(fd)
      pvar = var_wav(fd)
      j = 0
      Do i = 1, fopt%niter
         j = j + 1
         Call Ab_system(fd, dep, fopt%weigh, A, b, d)
         Call solve_Ab(A, b, dp)
         Call update(fd, dep, dp)
         Call w_ydependent(fd)
         If (fopt%weigh /= 0) Call w_bisquare(fd, dep, fopt)
         var = var_wav(fd)
         If (var > pvar .and. begindamp) Then ! make sure the sum of squares decreases
            If (begindamp) Call update(fd, dep, -dp)
            d = d*fopt%v
         Else ! decrease the damping factor to have better convergence
            begindamp = .false.
            d = d/fopt%v
         End If
         !write(*,*) i, d, var
         If( abs(pvar-var)/pvar < fopt%tol ) exit
         pvar = var
      End Do
      
      If (j == fopt%niter) Then
         write(*,'(A)') 'Max number of iterations hit. Check output and possibly repeat with option -n <max. num. it.>'
      Else
         write(*,'(" Fitting converged in ", I6, " steps.")') j
      End If
   End Subroutine
   
   Function var_wav(fd)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Double precision :: var_wav
         
      Integer :: ns, is, ntot, ny, iy
      
      ns = size(fd)
      ny = fd(1)%ny
      ntot = 0
      var_wav = 0.0d0
      Do is = 1, ns
         Do iy = 1, ny
            ntot = ntot + fd(is)%n
            var_wav = var_wav + fd(is)%var_w(iy)*dble(fd(is)%n)*fd(is)%yw(iy)
         End Do
      End Do
      var_wav = var_wav/dble(ntot)
   End Function
   
   Function var_av(fd)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Double precision :: var_av
         
      Integer :: ns, is, ntot, iy
      
      ns = size(fd)
      ntot = 0
      var_av = 0.0d0
      Do is = 1, ns
         ntot = ntot + fd(is)%n
         Do iy = 1, fd(is)%ny
            var_av = var_av + fd(is)%var(iy)*dble(fd(is)%n)*fd(is)%yw(iy)
         End Do
      End Do
      var_av = var_av/dble(ntot)
   End Function
   
   Function mad_av(fd)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Double precision :: mad_av
         
      Integer :: ns, is, ntot, iy
      
      ns = size(fd)
      ntot = 0
      mad_av = 0.0d0
      Do is = 1, ns
         ntot = ntot + fd(is)%n
         Do iy = 1, fd(is)%ny
            mad_av = mad_av + fd(is)%mad(iy)*dble(fd(is)%n)*fd(is)%yw(iy)
         End Do
      End Do
      mad_av = mad_av/dble(ntot)
   End Function
   
   Subroutine update(fd, dep, dp)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Double precision :: dp(dep%nr)
      
      Integer :: is, ns, i, ir, ofs, ip, iy
      Double precision :: fr
      
      ns = size(fd)
      ofs = 0
      Do is = 1, ns
         Do i = 1, fd(is)%np
            ir = dep%irow(ofs + i)
            fr = dep%fac(ofs + i)
            ip = fd(is)%ip(i)
            iy = fd(is)%iy(i)
            fd(is)%p(ip, iy) = fd(is)%p(ip, iy) + dp(ir)*fr
         End Do
         Call f_linexp_sl(fd(is))
         Call var_mad(fd(is))
         Call jac_linexp_sl(fd(is))
         ofs = ofs + fd(is)%np
      End Do
   End Subroutine
   
   
   Subroutine w_ydependent(fd)
      Implicit None
      Type(t_fitdata) :: fd(:)
      
      Integer :: is, ns, iy
      
      ns = size(fd)
      Do is = 1, ns
         Do iy = 1, fd(is)%ny
            fd(is)%yw = sum(fd(is)%mad)/fd(is)%mad(iy)/dble(fd(is)%ny)
         End Do
      End Do
   End Subroutine
   
   
   Subroutine w_bisquare(fd, dep, fopt)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Type(t_fit_options) :: fopt
      Double precision :: A(dep%nr, dep%nr)
      
      Integer :: is, ns, i, j, ir, ofsi, ic, iz, iy
      Double precision, parameter :: mad_normal=0.6745
      Double precision :: scale, h(fd(1)%ny), htot(fd(1)%ny), r, fr, fc
      
      
      ns = size(fd)
      !mad = mad_av(fd)
      !scale = mad/mad_normal
      Call Ab_system(fd, dep, 0, A)
      Call getinv_A(A)
      htot = 0.0d0
      ofsi = 0
      Do is = 1, ns
         Do iz = 1, fd(is)%n
            h = 0.0d0
            Do i = 1, fd(is)%np
               iy = fd(is)%iy(i)
               ir = dep%irow(ofsi+i)
               fr = dep%val(ofsi+i)
               Do j = 1, fd(is)%np
                  If (fd(is)%iy(j) == iy) Then
                     ic = dep%irow(ofsi+j)
                     fc = dep%val(ofsi+j)
                     h(iy) = h(iy) + fd(is)%jac(i, iz) * fr * A(ir, ic) * fc * fd(is)%jac(j, iz) * fd(is)%yw(iy)
                     ! h are the leverages that are needed to calculate robust weights (diagonal values of the Hat-matrix).
                     ! these values are calculated without the full weighing, only the y-dependent weigth is considered.
                  End If
               End Do
            End Do
            htot = htot + h
            Do iy = 1, fd(is)%ny
               scale = fd(is)%mad(iy)/mad_normal
               r = (fd(is)%y(iz, iy)-fd(is)%f(iz, iy))/(fopt%tune*scale*sqrt(1.0d0-h(iy))) ! adjusted residuals
               If ( abs(r) < 1.0d0 ) Then
                  fd(is)%w(iz, iy) = (1.0d0-r**2)**2 ! robust bisquare weights
               Else
                  fd(is)%w(iz, iy) = 0.0d0
               End If
            End Do
         End Do
         ofsi = ofsi + fd(is)%np
      End Do
      !write(*,*) sum(htot)
   End Subroutine
   
   Subroutine covar_pars(fd, dep, covar)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Double precision :: covar(dep%nr, dep%nr)
      
      Double precision :: JxAinv(dep%nz, dep%ny, dep%nr)
      Double precision :: A(dep%nr, dep%nr)
      Double precision :: jac(dep%nr, dep%nz, dep%ny)
      Double precision :: dum(dep%nz, dep%ny, dep%nr), w(dep%nz, dep%ny)
      Double precision :: d
      Integer :: iz, iy, i, j, is, ofs
      
      
      Call jac_system(fd, dep, jac)
      ofs = 0
      Do is = 1, size(fd)
         Do iz = 1, fd(is)%n
            Do iy = 1, fd(is)%ny
               w(ofs+iz, iy) = fd(is)%w(iz, iy)*fd(is)%yw(iy)
            End Do
         End Do
         ofs = ofs + fd(is)%n
      End Do
      
      Call Ab_system(fd, dep, 1, A)
      Call getinv_A(A)
      
      JxAinv = 0.0d0
      Do iy = 1, dep%ny
         Do iz = 1, dep%nz
            Do i = 1, dep%nr
               Do j = 1, dep%nr
                  JxAinv(iz, iy, i) = JxAinv(iz, iy, i) + jac(j, iz, iy)*A(j, i)*w(iz, iy)
               End Do
            End Do
         End Do
      End Do
      dum = 0.0d0
      ofs = 0
      Do is = 1, size(fd)
         Do iz = 1, fd(is)%n
            Do iy = 1, fd(is)%ny
               d = fd(is)%y(iz, iy) - fd(is)%f(iz, iy)
               Do i = 1, dep%nr
                  dum(iz+ofs, iy, i) = dum(iz+ofs, iy, i) + d*d*JxAinv(iz+ofs, iy, i)
               End Do
            End Do
         End Do
         ofs = ofs + fd(is)%n
      End Do
      covar = 0.0d0
      Do iy = 1, dep%ny
         Do iz = 1, dep%nz
            Do i = 1, dep%nr
               Do j = 1, dep%nr
                  covar(i, j) = covar(i, j) + JxAinv(iz, iy, i)*dum(iz, iy, j)
               End Do
            End Do
         End Do
      End Do
      !Do i = 1, dep%nr
      !   write(*,'(20g12.4)') covar(i,:)
      !End Do
      
   End Subroutine
   
   
   Subroutine convert_pars(fd, dep)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
         
      Double precision :: covar(dep%nr, dep%nr)
      Integer :: i, is, iy, ip, ir, iy2
      
      Call covar_pars(fd, dep, covar)
      Do i = 1, dep%nc
         ir = dep%irow(i)
         is = dep%is(i)
         iy = dep%iy(i)
         ip = dep%ip(i)
         Do iy2 = 1, fd(is)%ny
            If (iy == iy2) Then
               select case (ip)
                  case(1)
                  fd(is)%A(iy) = fd(is)%A(iy) + fd(is)%p(ip, iy)
                  fd(is)%dA(iy) = fd(is)%A(iy) + sqrt(covar(ir, ir))
                  case(2)
                  fd(is)%A(iy) = fd(is)%A(iy) - fd(is)%p(ip, iy)
                  fd(is)%dA(iy) = fd(is)%A(iy) + sqrt(covar(ir, ir))
                  fd(is)%B(iy) = fd(is)%p(ip, iy)
                  fd(is)%dB(iy) = sqrt(covar(ir, ir))
                  case(3)
                  fd(is)%A(iy) = fd(is)%A(iy) - fd(is)%p(ip, iy)
                  fd(is)%dA(iy) = fd(is)%A(iy) + sqrt(covar(ir, ir))
                  fd(is)%C(iy) = fd(is)%p(ip, iy)
                  fd(is)%dC(iy) = sqrt(covar(ir, ir))
                  case(5)
                  fd(is)%D(iy) = fd(is)%p(ip, iy)
                  fd(is)%dD(iy) = sqrt(covar(ir, ir))
                  case(4)
                  fd(is)%l(iy) = fd(is)%p(ip, iy)
                  fd(is)%dl(iy) = sqrt(covar(ir, ir))
               end select
            End If
         End Do
      End Do
   End Subroutine
   
   
   Subroutine get_vfpars(fd, dep, pars, cpera)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Type(t_parameters), Target :: pars(:)
      Double Precision :: cpera
         
      Double precision :: covar(dep%nr, dep%nr)
      Integer :: i, is, iy, ip, ir
      Type(t_parameters), Pointer :: p
      
      Call covar_pars(fd, dep, covar)
      Do i = 1, dep%nc
         ir = dep%irow(i)
         is = dep%is(i)
         iy = dep%iy(i)
         ip = dep%ip(i)
         p => pars(is)
         If (ip == 1 .and. iy == 1 .and. p%lsf_j%need) Then
            p%beta%val = fd(is)%p(ip, iy)
            p%beta%d = sqrt(covar(ir, ir))
            !write(*,'(2I3, "   beta =", g15.5, "std", g15.5)') is, iy, p%beta%val, p%beta%d
         End If
         If (ip == 4 .and. iy == 1 .and. p%lsf_j%need) Then
            p%lsf_j%val = fd(is)%p(ip, iy)
            p%lsf_j%d = sqrt(covar(ir, ir))
            !write(*,'(2I3, "    lsf =", g15.5, "std", g15.5)') is, iy, p%lsf_j%val, p%lsf_j%d
         End If
         If (ip == 4 .and. iy == 2 .and. p%lsf_mu%need) Then
            p%lsf_mu%val = fd(is)%p(ip, iy)
            p%lsf_mu%d = sqrt(covar(ir, ir))
            !write(*,'(2I3, "    lsf =", g15.5, "std", g15.5)') is, iy, p%lsf_mu%val, p%lsf_mu%d
         End If
         If (ip == 4 .and. iy == 4 .and. p%lsf_n%need) Then
            p%lsf_n%val = fd(is)%p(ip, iy)
            p%lsf_n%d = sqrt(covar(ir, ir))
            !write(*,'(2I3, "    lsf =", g15.5, "std", g15.5)') is, iy, p%lsf_n%val, p%lsf_n%d
         End If
         If (ip == 5 .and. iy == 2 .and. p%rho%need) Then
            p%rho%val = abs(fd(is)%p(ip, iy))/cpera*0.1d0
            p%rho%d = sqrt(covar(ir, ir))/cpera*0.1d0
            !write(*,'(2I3, "    rho =", g15.5, "std", g15.5)') is, iy, p%rho%val, p%rho%d
         End If
         !write(*,'(3I5, 2g15.4)') is, iy, ip, fd(is)%p(ip, iy), sqrt(covar(ir, ir))
      End Do
   End Subroutine
   
   Subroutine get_shpars(fd, dep, pars)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Type(t_parameters), Target :: pars(:)
         
      Double precision :: covar(dep%nr, dep%nr)
      Integer :: i, is, iy, ip, ir
      Type(t_parameters), Pointer :: p
      
      !Call covar_pars(fd, dep, covar)
      Do i = 1, dep%nc
         ir = dep%irow(i)
         is = dep%is(i)
         iy = dep%iy(i)
         ip = dep%ip(i)
         p => pars(is)
         If (ip == 1) Then
            p%shang%val = fd(is)%A(iy)
            p%shang%d = fd(is)%dA(iy)
            !write(*,'(2I3, "    SH_ang =", g15.5, "std", g15.5)') is, iy, p%shang%val, p%shang%d
         End If
         If (ip == 4) Then
            p%lsf_sh%val = fd(is)%l(iy)
            p%lsf_sh%d = fd(is)%dl(iy)
            !write(*,'(2I3, "      l_SH =", g15.5, "std", g15.5)') is, iy, p%lsf_sh%val, p%lsf_sh%d
         End If
      End Do
   End Subroutine
   
   Subroutine write_pars(pars, fopt, materials, zint, qa)
      Implicit None
      Type(t_parameters), Target :: pars(:)
      Character :: materials(:)
      Type(t_fit_options) :: fopt
      Double Precision :: zint(:), qa(3)
      
      Integer :: i, ns
      character (Len=100) :: work
      character (Len=60) :: line
      Type(t_parameters), Pointer :: p
      
      line = '("--------------------------------------------------------")'
      
      open(205, file="qfout", action="write")
      write(205,*) "Fitting options:"
      write(205,'(" Tolerance:  ", g15.5)') fopt%tol
      If (fopt%damp == 0.0d0) Then
         write(205,*) "No Levenberg-Marquardt damping."
      Else
         write(205,*) "Levenberg-Marquardt damping with:"
         write(205,'("    Damping parameter: ", f9.5)') fopt%damp
         write(205,'("    Marquardt factor:  ", f9.5)') fopt%v
      End If
      
      If (fopt%weigh == 0) Then
         write(205,*) "No robust bisquare weighing."
      Else
         write(205,'(" Robust bisquare weighing with tuning parameter ", f7.4)') fopt%tune
      End If
      
      If (fopt%lsf_opt == 1) work = "independently."
      If (fopt%lsf_opt == 2) work = "coupled."
      If (fopt%lsf_opt == 3) work = "to the chemical potentials."
      If (fopt%lsf_opt == 4) work = "to the spin current."
      If (fopt%lsf_opt == 3) work = "to the spin density."
      write(205,'(x,A,A)') "The spin diffusion length is fitted ", trim(work)
      write(205,line)
      write(205,'(" Quantization axis: ", 3f9.5)') qa
      write(205,line)
      
      ns = size(pars)
      Do i = 1, ns
         write(205,'("Slab ", I4, ",    Material ", A1)') i, materials(i)
         write(205,'(A10," = ", f13.5, 20x, A10)') 'zl', zint(i), 'nm'
         write(205,'(A10," = ", f13.5, 20x, A10)') 'zr', zint(i+1), 'nm'
         p => pars(i)
         Call wp_helper(205, p%rho)
         Call wp_helper(205, p%beta)
         Call wp_helper(205, p%lsf_j)
         Call wp_helper(205, p%lsf_mu)
         Call wp_helper(205, p%lsf_n)
         
         Call wp_helper(205, p%SHang)
         Call wp_helper(205, p%lsf_sh)
         
         write(205,line)
      End Do
      close(205)
      
   End Subroutine
   
   Subroutine wp_helper(file, par)
      Implicit None
      Type(t_one_parameter) :: par
      Integer :: file
      
      If (par%need) write(file, '(A10," = ", f13.5, "  std", g15.5, A10)') trim(par%name), par%val, par%d, trim(par%units)
   End Subroutine
   
   Subroutine direct_search(fd, dep, nsearch)
      Implicit None
      Type(t_fitdata) :: fd(:)
      Type(t_dependmat) :: dep
      Integer :: nsearch
      
      Integer :: ns, ir, ic, i, jr, besti
      Double precision :: A(dep%nr, dep%nr), b(dep%nr)
      Double precision :: dl, dp(dep%nr)
      Double precision :: mad, bestmad
      Logical :: dosearch, lsf_flag(dep%nr)
      
      ns = size(fd)
      dl = 0.0d0
      Do i = 1, ns
         dl = max(dl, (fd(i)%sd%zr - fd(i)%sd%zl) )
      End Do
      dp = 0.0d0
      lsf_flag = .false.
      Do ic = 1, dep%nc
         If (dep%ip(ic) == 4) Then
            ir = dep%irow(ic)
            dp(ir) = dl
            lsf_flag(ir) = .true.
         End If
      End Do
      Call update(fd, dep, dp)
      bestmad = mad_av(fd)
      dl = dl/dble(nsearch)
      
      Do ir = 1, dep%nr
         dp = 0.0d0
         dosearch = .false.
         Do ic = 1, dep%nc
            If (dep%irow(ic) == ir .and. dep%ip(ic) == 4) Then
               dosearch = .true.
               exit
            End If
         End Do
         If (dosearch) Then
            Do i = 1, nsearch
               Call Ab_system(fd, dep, 0, A, b)
               Do jr = 1, dep%nr
                  If (lsf_flag(jr)) Then
                     A(:,jr) = 0.0d0
                     A(jr,:) = 0.0d0
                     A(jr,jr) = 1.0d0
                  End If
               End Do
               
               Call solve_Ab(A, b, dp)
               Do jr = 1, dep%nr
                  If (lsf_flag(jr)) Then
                     dp(jr) = 0.0d0
                  End If
               End Do
               
               Call update(fd, dep, dp)
               
               mad = mad_av(fd)
               If (mad < bestmad) Then
                  bestmad = mad
                  besti = i
               End If
               
               dp = 0.0d0
               dp(ir) = -dl
               Call update(fd, dep, dp)
            End Do
            
            dp = 0.0d0
            dp(ir) = dble(nsearch-besti+1)*dl
            Call update(fd, dep, dp)
            
            Call Ab_system(fd, dep, 0, A, b)
            Do jr = 1, dep%nr
               If (lsf_flag(jr)) Then
                  A(:,jr) = 0.0d0
                  A(jr,:) = 0.0d0
                  A(jr,jr) = 1.0d0
               End If
            End Do
            
            Call solve_Ab(A, b, dp)
            Do jr = 1, dep%nr
               If (lsf_flag(jr)) Then
                  dp(jr) = 0.0d0
               End If
            End Do
            
            Call update(fd, dep, dp)
         End If
      End Do
   End Subroutine
   
end Module qfitlib
