#include "math_def.h"
Module iter_solver
      Implicit None
      Include 'zmumps_struc.h'

      Type t_iter_opts
         Integer :: debug = 1
         Integer :: restart = 100, maxiter = 10000, opt1 = 0, opt2 = 0
         Real (Kind=DEF_DBL_PREC) :: tol = 1.0d-8
      End Type t_iter_opts

      Type t_gmres_data
         Integer :: icntl (8), info (3), irc (5)
         Real (Kind=DEF_DBL_PREC) :: cntl (5), rinfo (2)
         Integer :: restart, lwork, n, maxit, nloc
         Complex (Kind=DEF_DBL_PREC), Allocatable :: work (:)
      End Type t_gmres_data

      Type t_iter_local_data
         Integer :: alloc
         Complex (Kind=DEF_DBL_PREC), Pointer :: a (:)
         Integer, Pointer :: ir (:), jc (:)
         Integer :: locnrow, locncol, locnnz
         Integer :: startcol, startrow, gn
         Integer, Allocatable :: rowoffs (:), nrows (:), coloffs (:), ncols (:)
         Complex (Kind=DEF_DBL_PREC), Pointer :: mmbuf1 (:, :), mmbuf2 (:, :)
         Complex (Kind=DEF_DBL_PREC), Pointer :: mmsin (:, :), mmsout (:, :)
         Complex (Kind=DEF_DBL_PREC), Pointer :: mmdin (:, :), mmdout (:, :)

         Complex (Kind=DEF_DBL_PREC), Allocatable :: lemb (:, :), remb (:, :), rhs (:, :)
         Integer, Allocatable :: lembofs (:), lembgn (:), rembofs (:), rembgn (:)
         Integer, Allocatable :: rhsmap (:)
         Logical, Allocatable :: rhsflag (:)
         Integer :: nrhs, globrhs, lrhs
         Integer :: les, Len, res, ren, ldmi, ldmo
         Integer :: gsz, comm, myid, myidf
         Logical :: boss
         Complex (Kind=DEF_DBL_PREC), Allocatable :: p1 (:), p2 (:)

         Integer :: rsz, lsz, dbsz
         Integer :: dbg
         Type (t_gmres_data), Pointer :: gmres
      End Type t_iter_local_data

      Interface mdtype
         Module Procedure mdtype_int, mdtype_int1, mdtype_dbl, mdtype_dbl1
      End Interface


Contains

      Subroutine sp_solve_iter (ham, lando, rando, task, wf, root, solvecomm, opt)
         Use sparselib
         Use ando_module
         use mpi
         Implicit None
         !Include 'mpif.h'
         Type (t_iter_opts), Intent (In) :: opt
         Type (zcsrmat), Intent (Inout), Target :: ham
         Type (zdensemat), Intent (Inout) :: wf
         Type (t_ando_sollution), Intent (Inout) :: lando, rando
         Integer, Intent (In) :: root, solvecomm, task
!!$ Local vars
         Type (t_iter_local_data) :: dat
         Integer :: ierr, comm, gsz, key, myid, i1
         Integer :: i,j
!!$          Complex (Kind=DEF_DBL_PREC), Pointer :: v1 (:), v2 (:), v3 (:)


         Call mpi_comm_rank (solvecomm, i1, ierr)
         If (i1 == root) key = 0
         key = 1

         Call MPI_COMM_SPLIT (solvecomm, 1, key, comm, ierr)
         Call mpi_comm_size (comm, gsz, ierr)
         Call mpi_comm_rank (comm, myid, ierr)

         If ((i1 == root) .And. (myid /= 0)) Then
            Write (*,*) 'New comm rank assign error1!'
            Write (*,*) 'myid=', myid
            Write (*,*) 'orig_id=', i1
            Write (*,*) 'root=', root
         End If

         dat%comm = comm
         dat%gsz = gsz
         dat%myid = myid
         dat%myidf = myid + 1
         dat%dbg = opt%debug
         dat%boss = (myid == 0)

!!$          Call testmul_init (dat, ham, lando%emb, rando%emb, v1, v2, v3, 1)

!!$           set preconditioning matrix
call dump_matrix(ham,'ham')
stop
         Allocate(dat%p1(ham%nrow))
         do i=1,ham%nrow
           dat%p1(i) = 1d0
           do j=ham%ir(i), ham%ir(i+1)-1
             if (ham%jc(j) == i) dat%p1(i) = ham%a(j)
           enddo
         enddo 
         
         Call distribute_mtxA (ham, dat)
         Call distribute_emb (lando%emb, rando%emb, dat)
         Call distribute_rhs (lando%bound, rando%bound, dat, task)
         Call multiply_alloc (dat)
         lando%haveEmb = 0
         lando%haveBound = 0
         rando%haveEmb = 0
         rando%haveBound = 0

!!$          Call multiply (dat, v1, v3, .True.)
!!$          Call mdtype (dat%comm, maxval(Abs(v2(:)-v3(:))), 'res')
         Call gmres_solver_init (dat, opt)
         Call alloc (wf, dat%gmres%n, dat%nrhs)
         Call gmres_solver_exec (dat,wf%bl)
         Call MPI_Barrier (dat%comm, ierr)
         !Call mpi_finalize (ierr)

         Deallocate(dat%gmres%work)
         Deallocate(dat%gmres)
      End Subroutine sp_solve_iter


      Subroutine gmres_solver_init (dat, opt)
         Use sparselib
         use mpi
         Implicit None
         !Include 'mpif.h'
         Type (t_iter_local_data) :: dat
         Type (t_iter_opts), Intent (In) :: opt
!!$           Local
         Integer :: m, n
!!$          Real (Kind=DEF_DBL_PREC) :: scl
         Integer :: scl, ierr
         Integer, Pointer :: icntl (:)
         Real (Kind=DEF_DBL_PREC), Pointer :: cntl (:)
         Character (Len=256) :: cwork
         Integer, Parameter :: debf = 17
         Allocate (dat%gmres)

         icntl => dat%gmres%icntl
         cntl => dat%gmres%cntl

         Call init_zgmres (icntl, cntl)
!!$      Tune some parameters
         icntl (2) = 0
         If (opt%debug > 0) icntl (2) = 6

         icntl (3) = 0

         If (opt%debug > 1) Then
            Write (cwork, '(i3.3)') dat%myid
            Open (Unit=debf, File='isolv_log_'//trim(cwork), Action='write')
            icntl (3) = debf
            icntl (2) = debf
            icntl (1) = debf
         End If

!!$      Save the convergence history on standard output
!!$          icntl (3) = 40

!!$           Maximum number of iterations
         icntl (7) = opt%maxiter
!!$           preconditioner location
         icntl (4) = opt%opt1

!!$           orthogonalization scheme
         icntl (5) = opt%opt2
!!$          opt%opt2
!!$           initial guess
         icntl (6) = 0
!!$           residual calculation strategy at restart
         icntl (8) = 1


         cntl (1) = opt%tol

         m = opt%restart
         n = dat%gn

         dat%gmres%n = n
         dat%gmres%nloc = n
         dat%gmres%restart = m
         Call MPI_Bcast (dat%gmres%restart, 1, MPI_INTEGER, 0, dat%comm, ierr)
         Allocate(dat%gmres%work(n))

!!$          scl = opt%dopt1
         scl = 1
!!$ dat%ist1 == lwork
         If (icntl(5) == 0 .Or. icntl(5) == 1) Then
            If (icntl(8) == 1) dat%gmres%lwork = scl * (m*m+m*(n+5)+5*n+2)
            If (icntl(8) == 0) dat%gmres%lwork = scl * (m*m+m*(n+5)+6*n+2)
         Else
            If (icntl(8) == 1) dat%gmres%lwork = scl * (m*m+m*(n+5)+5*n+m+1)
            If (icntl(8) == 0) dat%gmres%lwork = scl * (m*m+m*(n+5)+6*n+m+1)
         End If

         Call MPI_Bcast (icntl, 8, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Bcast (cntl, 5, MPI_DOUBLE_PRECISION, 0, dat%comm, ierr)
         dat%gmres%maxit = icntl (7)
      End Subroutine gmres_solver_init

      Subroutine gmres_solver_exec (dat, sols)
         use mpi
         Implicit None
         !Include 'mpif.h'
         Type (t_iter_local_data) :: dat
         Complex (Kind=DEF_DBL_PREC) :: sols(:,:)
!!$           Local
         Integer :: n, lwork, m, info (3), irc (5)

         Complex (Kind=DEF_DBL_PREC), Allocatable, Target :: work (:)
         Real (Kind=DEF_DBL_PREC) :: rinfo (2), tm1
         Complex (Kind=DEF_DBL_PREC), Parameter :: rone = dcmplx (1.0d0, 0.0d0), rzero = dcmplx (0.0d0, &
        & 0.0d0)
         Logical :: hw
         Integer :: revcom, colx, coly, colz, nbscal, ierr, nit, num
         Integer, Parameter :: matvec = 1, precondLeft = 2, precondRight = 3, dotProd = 4
         Logical :: flag
         Integer :: i

         lwork = dat%gmres%lwork
         Allocate (work(lwork))

         n = dat%gn

         m = dat%gmres%restart
         colx = 1
         coly = 1
         colz = 1

        do num=1,dat%nrhs
          hw = (isolv_set_rhs(dat, num, work(n+1:2*n)) > 0)
          !work(1:n) = .5d0
          work(1:n) = 0d0

          tm1 = MPI_Wtime ()
!!$           write(*,*) hw
          nit = 0
          info(1)=1
          revcom = 99
          Do While (revcom .ne. 0 .and. info(1) .ne. 0)
          !print*,'num',num
             nit = nit + 1
             If (hw) Then
                Call drive_zgmres (n, n, m, lwork, work, irc, dat%gmres%icntl, dat%gmres%cntl, info, rinfo)
                revcom = irc (1)
                colx = irc (2)
                coly = irc (3)
                colz = irc (4)
                nbscal = irc (5)
             End If

!!$              write(*,*) 'hhh',revcom,maxval(Abs(work(colx:colx+n-1))),maxval(Abs(work(colz:colz+n-1)))
             Call multiply (dat, work(colx:colx+n-1), work(colz:colz+n-1), hw .And. (revcom .Eq. matvec))
!!$              write(*,*) 'ppp',maxval(Abs(work(colz:colz+n-1))),colz,colx

             If (revcom .Eq. precondLeft) Then
!!$                   calc non-preconditioned solution
!!$                   work(colz) <-- M^{-1} * work(colx) 
!!$                   For now the preconditioner will just be the Jacobi one, i.e. diagonal
                do i=1,n
                   work(colz+i-1) = work(colx+i-1)/dat%p1(i)
                enddo
             End If
             If (revcom .Eq. precondRight) Then
                Write (*,*) 'Right preconditioning is not implemented yet'
                Stop
             End If

             If (revcom .Eq. dotProd .And. hw) Then
!!$                   perform the scalar product
!!$                   work(colz) <-- work(colx) work(coly)
                Call zgemv ('C', n, nbscal, rone, work(colx:lwork), n, work(coly:lwork), 1, rzero, &
               & work(colz:lwork), 1)
             End If
!!$              Call MPI_Barrier (dat%comm, i)

!!$              exit

             if (revcom .eq. 0 .and. info(1) .eq. 0)  then
               sols(:,num) = work(1:n)
             endif

             If ( info(1) .ne. 0 .and. revcom < 1 .or. revcom > 4) hw = .False.
             Call MPI_AllReduce ( .Not. hw, flag, 1, mpi_logical, mpi_land, dat%comm, ierr)

             If (flag) Then
                Write (*,*) 'I need to handle this!', n
                Write (*,*) hw, dat%myid, revcom, info
                Write (*,*) rinfo, MPI_Wtime () - tm1
                Stop 'DID NOT CONVERGE ITERATIVE SOLUTION OF AX=B'
                hw = .False.
             End If
!!$              if (mod(nit,100)==0) then
!!$                  Write (*,*) dat%myid,maxval(abs(x-b))
!!$              endif
          End Do
!!$           write(*,*)     MPI_Wtime ()-tm1

        End Do
        Deallocate(work)

      End Subroutine gmres_solver_exec


      Subroutine testmul_init (dat, ham, lemb, remb, v1, v2, v3, flag)
         Use sparselib
         !use mpi
         Implicit None
         Include 'mpif.h'
         Type (zcsrmat), Intent (In), Target :: ham
         Complex (Kind=DEF_DBL_PREC), Pointer, Intent (Inout) :: v1 (:), v2 (:), v3 (:)
         Complex (Kind=DEF_DBL_PREC), Intent (In) :: lemb (:, :), remb (:, :)
         Type (t_iter_local_data) :: dat
         Integer :: flag
!!$           Local
         Integer :: key, i1, ierr
         Integer :: status (MPI_STATUS_SIZE)

         key = ham%nrow
         Call MPI_Bcast (key, 1, MPI_INTEGER, 0, dat%comm, ierr)

         Allocate (v1(key), v2(key), v3(key))
         Do i1 = 1, key
            v1 (i1) = dcmplx (i1, dat%myid)
         End Do

         If (dat%boss) Then
            Do i1 = 1, dat%gsz - 1
               Call MPI_recv (v3, key, mpi_double_complex, i1, i1, dat%comm, status, ierr)
               call spmatmul2 (ham, v3,v2)
               If (flag /= 0) Then
                  v2 (1:size(lemb, 1)) = v2 (1:size(lemb, 1)) + matmul (lemb, v3(1:size(lemb, 1)))
                  v2 (key-size(remb, 1)+1:key) = v2 (key-size(remb, 1)+1:key) + matmul (remb, &
                 & v3(key-size(remb, 1)+1:key))
               End If
               Call MPI_send (v2, key, mpi_double_complex, i1, i1+1024, dat%comm, status, ierr)
            End Do
            call spmatmul2 (ham, v1,v2)
            If (flag /= 0) Then
               v2 (1:size(lemb, 1)) = v2 (1:size(lemb, 1)) + matmul (lemb, v1(1:size(lemb, 1)))
               v2 (key-size(remb, 1)+1:key) = v2 (key-size(remb, 1)+1:key) + matmul (remb, v1(key-size(remb, &
              & 1)+1:key))
            End If
         Else
            Call MPI_send (v1, key, mpi_double_complex, 0, dat%myid, dat%comm, status, ierr)
            Call MPI_recv (v2, key, mpi_double_complex, 0, dat%myid+1024, dat%comm, status, ierr)
         End If

      End Subroutine testmul_init

      Subroutine multiply (dat, v, v1, hj)
         !use mpi
         Implicit None
         Include 'mpif.h'
         Type (t_iter_local_data) :: dat
         Complex (Kind=DEF_DBL_PREC), Intent (In), Target :: v (:)
         Complex (Kind=DEF_DBL_PREC), Intent (Out) :: v1 (:)
         Logical :: hj
!!$
         Logical :: jobs (dat%gsz)
         Integer :: i, nr, i1, j1, n1, n2, ierr, i2
         Complex (Kind=DEF_DBL_PREC), Allocatable :: tv2 (:)
         Complex (Kind=DEF_DBL_PREC), Pointer :: a (:), mat (:, :), mat2 (:, :)
         Integer, Pointer :: jc (:)
         Complex (Kind=DEF_DBL_PREC), Parameter :: cone = dcmplx (1.0d0, 0.0d0), czero = dcmplx (0.0d0, &
        & 0.0d0)
         Integer :: lsz, rsz, gn, nc, myc
         Call mpi_allgather (hj, 1, mpi_logical, jobs, 1, mpi_logical, 0, dat%comm, ierr)

!!$ Distr data for sparse-dense op
         nc = 0
!SMP$ DO SERIAL
         Do i = 1, dat%gsz
            If (jobs(i)) Then
               nc = nc + 1
               Call MPI_Scatterv (v, dat%ncols, dat%coloffs, mpi_double_complex, dat%mmsin(:, nc), &
              & dat%locncol, mpi_double_complex, i-1, dat%comm, ierr)
            End If
         End Do
         If (nc == 0) Return

         rsz = dat%rsz
         lsz = dat%lsz
         gn = dat%gn
         a => dat%a
         jc => dat%jc

!!$ SparseMM
         mat => dat%mmsin
         mat2 => dat%mmsout
         nr = dat%locnrow
!$omp parallel do
         Do i2 = 1, nc
            Do i1 = 1, nr
               mat2 (i1, i2) = czero
            End Do
         End Do
!$omp end parallel do
!!$ !$omp parallel do
         Do i2 = 1, nc
            Do i1 = 1, nr
               n1 = dat%ir (i1)
               n2 = dat%ir (i1+1) - 1
               Do j1 = n1, n2
                  mat2 (i1, i2) = mat2 (i1, i2) + mat (jc(j1), i2) * a (j1)
               End Do
            End Do
         End Do
!!$ !$omp end parallel do

         nc = 0
         myc = 1
         Do i = 1, dat%gsz
            If (jobs(i)) Then
               nc = nc + 1
               Call mpi_gatherv (mat2(:, nc), dat%locnrow, mpi_double_complex, v1, dat%nrows, dat%rowoffs, &
              & mpi_double_complex, i-1, dat%comm, ierr)
               If (i == dat%myidf) myc = nc
            End If
         End Do


!!$ Dense blocks MV
         Allocate (tv2(dat%dbsz))
!!$ Left
         mat => dat%mmdin
         mat2 => dat%mmdout
         mat (1:lsz, myc) = v (1:lsz)
         nc = 0
!SMP$ DO SERIAL
         Do i = 1, dat%gsz
            If (jobs(i)) Then
               nc = nc + 1
               Call MPI_Bcast (mat(:, nc), lsz, mpi_double_complex, i-1, dat%comm, ierr)
            End If
         End Do

         Call zgemm ('N', 'N', dat%len, nc, lsz, cone, dat%lemb, dat%len, mat, dat%ldmi, czero, mat2, &
        & dat%ldmo)

         nc = 0
!SMP$ DO SERIAL
         Do i = 1, dat%gsz
            If (jobs(i)) Then
               nc = nc + 1
               Call mpi_gatherv (mat2(:, nc), dat%len, mpi_double_complex, tv2, dat%lembgn, dat%lembofs, &
              & mpi_double_complex, i-1, dat%comm, ierr)
            End If
         End Do
         If (hj) v1 (1:lsz) = v1 (1:lsz) + tv2 (1:lsz)


!!$ Right
         n1 = gn - rsz + 1
         mat (1:rsz, myc) = v (n1:gn)
         nc = 0
!SMP$ DO SERIAL
         Do i = 1, dat%gsz
            If (jobs(i)) Then
               nc = nc + 1
               Call MPI_Bcast (mat(:, nc), rsz, mpi_double_complex, i-1, dat%comm, ierr)
            End If
         End Do

         Call zgemm ('N', 'N', dat%ren, nc, rsz, cone, dat%remb, dat%ren, mat, dat%ldmi, czero, mat2, &
        & dat%ldmo)

         nc = 0
!SMP$ DO SERIAL
         Do i = 1, dat%gsz
            If (jobs(i)) Then
               nc = nc + 1
               Call mpi_gatherv (mat2(:, nc), dat%ren, mpi_double_complex, tv2, dat%rembgn, dat%rembofs, &
              & mpi_double_complex, i-1, dat%comm, ierr)
            End If
         End Do
         If (hj) v1 (n1:gn) = v1 (n1:gn) + tv2 (1:rsz)

         Deallocate (tv2)
         Nullify (a, jc, mat, mat2)
      End Subroutine multiply


      Subroutine multiply_alloc (dat)
         use mpi
         Implicit None
         !Include 'mpif.h'
         Type (t_iter_local_data) :: dat
!!$ Local
         Logical :: set
         Integer :: in1, in2, on1, on2, n, m

         in1 = Max (dat%lsz, dat%rsz)
         on1 = dat%len
         in2 = dat%locncol
         on2 = dat%locnrow
         set = ((Max(in1, in2)+Max(on1, on2)) < (Max(in1, on2)+Max(on1, in2)))
         If (set) Then
            n = Max (in1, in2)
            m = Max (on1, on2)
         Else
            n = Max (in1, on2)
            m = Max (on1, in2)
         End If

         Allocate (dat%mmbuf1(n, dat%gsz), dat%mmbuf2(m, dat%gsz))
         dat%mmdin => dat%mmbuf1
         dat%mmdout => dat%mmbuf2
         dat%ldmi = n
         dat%ldmo = m

         If (set) Then
            dat%mmsin => dat%mmbuf1
            dat%mmsout => dat%mmbuf2
         Else
            dat%mmsin => dat%mmbuf2
            dat%mmsout => dat%mmbuf1
         End If
      End Subroutine multiply_alloc

      Subroutine distribute_rhs (lbnd, rbnd, dat, task)
         Use helpers
         use mpi
         Implicit None
         !Include 'mpif.h'
         Complex (Kind=DEF_DBL_PREC), Pointer, Intent (Inout) :: lbnd (:, :), rbnd (:, :)
         Type (t_iter_local_data) :: dat
         Integer :: task
!!$          Local
         Integer :: needl, needr, i, nrhs, tag, j, nt, nr, nl, ierr, n1
         Integer, Allocatable :: tmpint1 (:), tmpint2 (:), tmpint3 (:), tmpint4 (:)
         Logical, Allocatable :: tmplog (:)
         Complex (Kind=DEF_DBL_PREC), Pointer :: lb (:, :), rb (:, :)
         Complex (Kind=DEF_DBL_PREC), Allocatable :: zbuf (:, :)
         Allocate (tmpint1(dat%gsz), tmpint2(dat%gsz))

         If (dat%boss) Then
            needl = Mod (task, 2)
            needr = Mod (task, 4) / 2
!!$             Write (*,*) task, needr, needl
            nr = size (rbnd, 2) * needr
            nl = size (lbnd, 2) * needl
            nrhs = nl + nr
            Call calcredist (dat%gsz, nrhs, tmpint1, tmpint2)
            dat%globrhs = nrhs
            lb => lbnd
            rb => rbnd
         Else
            Allocate (lb(1, 1), rb(1, 1))
         End If

         Call MPI_Bcast (dat%globrhs, 1, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Scatter (tmpint2, 1, MPI_INTEGER, dat%nrhs, 1, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Scatter (tmpint1, 1, MPI_INTEGER, tag, 1, MPI_INTEGER, 0, dat%comm, ierr)

         dat%lrhs = Max (dat%lsz, dat%rsz)
         Allocate (dat%rhsmap(dat%nrhs))
         Allocate (dat%rhsflag(dat%nrhs))
         Allocate (dat%rhs(dat%lrhs, dat%nrhs))

         Do i = 1, dat%nrhs
            dat%rhsmap (i) = tag + i
         End Do

         If (dat%boss) Then
            Allocate (tmplog(nrhs))
            j = nl + 1
            Do i = 1, nrhs
               tmplog (i) = i < j
            End Do
         Else
            Allocate (tmplog(1))
         End If

         Call MPI_Scatterv (tmplog, tmpint2, tmpint1, mpi_logical, dat%rhsflag, dat%nrhs, mpi_logical, 0, &
        & dat%comm, ierr)
         Deallocate (tmplog)

         Allocate (tmpint3(dat%gsz), tmpint4(dat%gsz))
         If (dat%boss) Then
            tmpint3 (:) = 0
            j = nl
            n1 = 0
            Do i = 1, dat%gsz
               tmpint3 (i) = Min (tmpint2(i), j)
               j = j - tmpint3 (i)
               n1 = n1 + tmpint3 (i)
            End Do
!!$             Write (*,*) tmpint2 (:)
!!$             Write (*,*) tmpint3 (:)
!!$             Write (*,*) tmpint1 (:)
            tmpint1 (:) = tmpint1 (:) * dat%lsz
            tmpint4 (:) = tmpint3 * dat%lsz
!!$             Write (*,*) tmpint1 (:)
!!$             Write (*,*) tmpint4 (:)
!!$             Write (*,*) '---------'
         End If
         Call MPI_Scatter (tmpint3, 1, MPI_INTEGER, nt, 1, MPI_INTEGER, 0, dat%comm, ierr)


         Allocate (zbuf(dat%lsz, nt))
         Call MPI_Scatterv (lb(:, :), tmpint4, tmpint1, mpi_double_complex, zbuf(:, :), nt*dat%lsz, &
        & mpi_double_complex, 0, dat%comm, ierr)
         If (nt > 0) dat%rhs(1:dat%lsz, 1:nt) = zbuf (1:dat%lsz, 1:nt)
         Deallocate (zbuf, lb)

         If (dat%boss) Then
            Do i = 1, dat%gsz
               tmpint1 (i) = tmpint1 (i) - nl * dat%lsz
               If (tmpint1(i) < 0) tmpint1 (i) = 0
            End Do

            tmpint3 = tmpint2 (:) - tmpint3 (:)
            tmpint4 (:) = tmpint3 * dat%rsz
!!$             Write (*,*) tmpint3 (:)
!!$             Write (*,*) tmpint1 (:) / dat%rsz
!!$             Write (*,*) tmpint4 (:)
         End If

         tag = nt
         Call MPI_Scatter (tmpint3, 1, MPI_INTEGER, nt, 1, MPI_INTEGER, 0, dat%comm, ierr)

         Allocate (zbuf(dat%rsz, nt))
         Call MPI_Scatterv (rb, tmpint4, tmpint1, mpi_double_complex, zbuf, nt*dat%rsz, mpi_double_complex, &
        & 0, dat%comm, ierr)
!!$          Call mdtype (dat%comm, tag, 'tag')
!!$          Call mdtype (dat%comm, nt, 'nt')
         If (nt > 0) dat%rhs(1:dat%rsz, tag+1:tag+nt) = zbuf (1:dat%rsz, 1:nt)
         Deallocate (zbuf)

         Deallocate (tmpint1, tmpint2, tmpint3, tmpint4, rb)
!!$          Stop
!!$    If (dat%nrhs > 0) allocate (dat%vrhs(dat%gn))
!!$          Do i = 1, dat%nrhs
!!$             Write (*,*) 'SGN:', dat%myid, i, maxval (Abs(dat%rhs(:, i)))
!!$          End Do
      End Subroutine distribute_rhs

      Function isolv_set_rhs (dat, num, v) Result (resok)
         Implicit None
         Type (t_iter_local_data) :: dat
         Integer :: num, resok, n1, n2
         Complex (Kind=DEF_DBL_PREC) :: v (1:dat%gn)

         resok = 0
         If (num > dat%nrhs) Return
!!$          call mdtype(dat%comm,dat%nrhs,'NRHS')
!!$          stop
         If (dat%rhsflag(num)) Then
            n1 = 1
            n2 = dat%lsz
         Else
            n1 = dat%gn - dat%rsz + 1
            n2 = dat%rsz
         End If
!!$          Call mdtype (dat%comm, n1, 'n1')
!!$          Call mdtype (dat%comm, n2, 'n2')
         v (:) = dcmplx (0.0d0, 0.0d0)
         v (n1:n1+n2-1) = dat%rhs(1:n2, num)
!!$          Call mdtype (dat%comm, maxval(Abs(v)), 'mB:')
         resok = 1
      End Function isolv_set_rhs


      Subroutine distribute_emb (lemb, remb, dat)
         Use helpers
         use mpi
         Implicit None
         !Include 'mpif.h'
         Complex (Kind=DEF_DBL_PREC), Pointer, Intent (Inout) :: lemb (:, :), remb (:, :)
         Type (t_iter_local_data) :: dat
!!$ local
         Integer :: ierr
         Integer :: i
         Integer, Allocatable :: tmpint1 (:), tmpint2 (:)
         Allocate (tmpint1(dat%gsz+1), tmpint2(dat%gsz+1))

         If (dat%boss) Then
            tmpint1 (1) = size (lemb, 1)
            tmpint1 (2) = size (remb, 1)
         End If

         Call MPI_Bcast (tmpint1, 2, MPI_INTEGER, 0, dat%comm, ierr)
         dat%lsz = tmpint1 (1)
         dat%rsz = tmpint1 (2)
         Allocate (dat%lembofs(dat%gsz), dat%lembgn(dat%gsz))
         Allocate (dat%rembofs(dat%gsz), dat%rembgn(dat%gsz))

         If (dat%boss) Call calcredist (dat%gsz, dat%lsz, dat%lembofs, dat%lembgn)
         If (dat%boss) Call calcredist (dat%gsz, dat%rsz, dat%rembofs, dat%rembgn)

         Call MPI_Bcast (dat%lembofs, dat%gsz, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Bcast (dat%lembgn, dat%gsz, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Bcast (dat%rembofs, dat%gsz, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Bcast (dat%rembgn, dat%gsz, MPI_INTEGER, 0, dat%comm, ierr)
         dat%les = dat%lembofs(dat%myidf) + 1
         dat%len = dat%lembgn(dat%myidf)
         dat%res = dat%rembofs(dat%myidf) + 1
         dat%ren = dat%rembgn(dat%myidf)
         Allocate (dat%lemb(dat%len, dat%lsz))
         Allocate (dat%remb(dat%ren, dat%rsz))
         If ( .Not. dat%boss) allocate (lemb(1, dat%lsz), remb(1, dat%rsz))

         Do i = 1, dat%lsz
            Call MPI_Scatterv (lemb(1, i), dat%lembgn, dat%lembofs, mpi_double_complex, dat%lemb(1, i), &
           & dat%len, mpi_double_complex, 0, dat%comm, ierr)
         End Do
         Deallocate (lemb)
         Do i = 1, dat%rsz
            Call MPI_Scatterv (remb(1, i), dat%rembgn, dat%rembofs, mpi_double_complex, dat%remb(1, i), &
           & dat%ren, mpi_double_complex, 0, dat%comm, ierr)
         End Do
         Deallocate (remb)

         dat%dbsz = Max (dat%lsz, dat%rsz)

      End Subroutine distribute_emb

      Subroutine distribute_mtxA (a, dat)
         Use helpers
         Use sparselib
         !use mpi
         Implicit None
         Include 'mpif.h'
         Type (zcsrmat), Intent (Inout) :: a
         Type (t_iter_local_data) :: dat
!!$ local
         Integer :: ierr
         Integer :: i, i1, i2
         Integer, Allocatable :: tmpint1 (:), tmpint2 (:)

         Allocate (dat%nrows(dat%gsz), dat%rowoffs(dat%gsz))
         Allocate (dat%ncols(dat%gsz), dat%coloffs(dat%gsz))

         Allocate (tmpint1(dat%gsz), tmpint2(dat%gsz))
         If (dat%boss) Call calcredist (dat%gsz, a%nrow, dat%rowoffs, dat%nrows)
         dat%gn = a%nrow
         Call MPI_Bcast (dat%rowoffs, dat%gsz, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Bcast (dat%nrows, dat%gsz, MPI_INTEGER, 0, dat%comm, ierr)
         Call MPI_Bcast (dat%gn, 1, MPI_INTEGER, 0, dat%comm, ierr)

         dat%locnrow = dat%nrows(dat%myidf)
         dat%startrow = dat%rowoffs(dat%myidf) + 1

!!$          Call mdtype (dat%comm, dat%locnrow, 'nrow')
!!$          Call mdtype (dat%comm, dat%startrow, 'srow')
         Allocate (dat%ir(dat%locnrow+1))
         tmpint1 (1:dat%gsz) = dat%nrows(1:dat%gsz) + 1

         Call MPI_Scatterv (a%ir, tmpint1, dat%rowoffs, MPI_INTEGER, dat%ir, dat%locnrow+1, MPI_INTEGER, 0, &
        & dat%comm, ierr)

         dat%locnnz = dat%ir(dat%locnrow+1) - dat%ir(1)

         Allocate (dat%jc(dat%locnnz))
         Allocate (dat%a(dat%locnnz))


         If (dat%boss) Then
            Do i = 1, dat%gsz
               tmpint1 (i) = a%ir(dat%rowoffs(i)+1) - 1
               tmpint2 (i) = a%ir(dat%rowoffs(i)+dat%nrows(i)+1) - a%ir(dat%rowoffs(i)+1)
            End Do
         End If
         If (dat%boss) deallocate (a%ir)

         Call MPI_Scatterv (a%jc, tmpint2, tmpint1, MPI_INTEGER, dat%jc, dat%locnnz, MPI_INTEGER, 0, &
        & dat%comm, ierr)
         If (dat%boss) deallocate (a%jc)

         Call MPI_Scatterv (a%a, tmpint2, tmpint1, mpi_double_complex, dat%a, dat%locnnz, mpi_double_complex, &
        & 0, dat%comm, ierr)
         If (dat%boss) deallocate (a%a)
         a%alloc = 0

         dat%startcol = minval (dat%jc)
         dat%locncol = maxval (dat%jc) - dat%startcol + 1
         Call MPI_Allgather (dat%startcol-1, 1, MPI_INTEGER, dat%coloffs, 1, MPI_INTEGER, 0, dat%comm, ierr)
         Call mpi_allgather (dat%locncol, 1, MPI_INTEGER, dat%ncols, 1, MPI_INTEGER, 0, dat%comm, ierr)

         i1 = dat%ir (1) - 1
         i2 = dat%locnrow + 1
!$omp parallel do
         Do i = 1, i2
            dat%ir (i) = dat%ir(i) - i1
         End Do
!$omp end parallel do
         i1 = dat%startcol - 1
         i2 = dat%locnnz
!$omp parallel do
         Do i = 1, i2
            dat%jc (i) = dat%jc(i) - i1
         End Do
!$omp end parallel do
         Deallocate (tmpint1, tmpint2)
!!$          Call mdtype (comm, dat%jc, 'locjc')



!!$          Call mdtype (comm, a%nnz, 'nnz0')
!!$          Call mdtype (comm, dat%locnnz, 'nnz')
!!$          Call mdtype (comm, dat%locncol, 'locncols')
!!$          Call mdtype (comm, dat%startcol, 'startcol')
!!$          Call mdtype (comm, dat%ncols, 'ncols')

      End Subroutine distribute_mtxA

      Subroutine mdtype_int1 (comm, val, msg)
         use mpi
         Implicit None
         !Include 'mpif.h'
         Integer :: myid, i, comm, ierr, gsz
         Integer :: val (:)
         Character (Len=*) :: msg
         Call mpi_comm_rank (comm, myid, ierr)
         Call mpi_comm_size (comm, gsz, ierr)
         Do i = 1, gsz
            If (myid == (i-1)) Then
               Write (*, '(A," data from id=",i3)') msg, i - 1
               Write (*,*) val (:)
!!$                flush (6)
            End If
            Call MPI_Barrier (comm, ierr)
         End Do

      End Subroutine mdtype_int1

      Subroutine mdtype_dbl1 (comm, val, msg)
         use mpi
         Implicit None
         !Include 'mpif.h'
         Integer :: myid, i, comm, ierr, gsz
         Real (Kind=DEF_DBL_PREC) :: val (:)
         Character (Len=*) :: msg
         Call mpi_comm_rank (comm, myid, ierr)
         Call mpi_comm_size (comm, gsz, ierr)
         Do i = 1, gsz
            If (myid == (i-1)) Then
               Write (*, '(A," data from id=",i3)') msg, i - 1
               Write (*,*) val (:)
!!$                flush (6)
            End If
            Call MPI_Barrier (comm, ierr)
         End Do

      End Subroutine mdtype_dbl1

      Subroutine mdtype_int (comm, val, msg)
         use mpi
         Implicit None
         !Include 'mpif.h'
         Integer :: myid, i, comm, ierr, gsz
         Integer :: val
         Character (Len=*) :: msg
         Call mpi_comm_rank (comm, myid, ierr)
         Call mpi_comm_size (comm, gsz, ierr)
         Do i = 1, gsz
            If (myid == (i-1)) Then
               Write (*, '(A," data from id=",i3," : ",i6)') msg, i - 1, val
!!$                flush (6)
            End If
            Call MPI_Barrier (comm, ierr)
         End Do
      End Subroutine mdtype_int

      Subroutine mdtype_dbl (comm, val, msg)
         use mpi
         Implicit None
         !Include 'mpif.h'
         Integer :: myid, i, comm, ierr, gsz
         Real (Kind=DEF_DBL_PREC) :: val
         Character (Len=*) :: msg
         Call mpi_comm_rank (comm, myid, ierr)
         Call mpi_comm_size (comm, gsz, ierr)
         Do i = 1, gsz
            If (myid == (i-1)) Then
               Write (*, '(A," data from id=",i3," : ",g15.7)') msg, i - 1, val
!!$                flush (6)
            End If
            Call MPI_Barrier (comm, ierr)
         End Do
      End Subroutine mdtype_dbl


End Module iter_solver
