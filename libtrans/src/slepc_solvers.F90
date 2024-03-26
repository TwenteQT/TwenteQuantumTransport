#include "math_def.h"
#include <slepc/finclude/slepceps.h>

  Module slepc_solvers
  Use mpi
  Use sparselib
  Use slepceps

  Implicit None

  Type t_petsc_opt
    Integer :: solver = 1 ! 1) GMRES
    Integer :: preco  = 1 ! 1) ILU 
    Real (Kind=DEF_DBL_PREC) :: tol=1d-8
  End Type t_petsc_opt

  Contains

    Subroutine slepc_solve(ham, root, comm, k, ev)
      Implicit None
      Type (zcsrmat) :: ham
      Integer :: root, comm
      Integer :: k !the number of eigenvalues
      Complex (Kind=DEF_DBL_PREC) :: ev(k)
      !!$ local
      Integer :: commsize, err 
      Integer :: sols 
      Integer :: ei
      Complex (Kind=DEF_DBL_PREC), Pointer :: temp(:)
      Mat               mat
      EPS               eps
      !PetscScalar       eigr, eigi
      PetscComplex      eigr, eigi
      Vec               vr, vi
      PetscErrorCode    ierr
      PetscInt          m, n, i(ham%nrow+1), j(ham%nnz)
      PetscInt          nev, iev
      !PetscScalar       a(ham%nnz)
      PetscComplex      a(ham%nnz)


      ! Set PETSC_COMM_WORLD to subset of MPI_COMM_WORLD
      PETSC_COMM_WORLD = comm

      ! Check comm size: at this time we don't support parallel matrices
      call MPI_COMM_SIZE(comm, commsize, err)
      if (commsize > 1) stop 'We currently do not support parallel matrices'

      ! Initialize PETSc
      Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      Call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      If (ierr .ne. 0) then
        stop 'Unable to initialize PETSc'
      endif

      ! create PETSc matrix
      i(:) = ham%ir(:)-1 ! 0-based index for C
      j(:) = ham%jc(:)-1 ! idem
      a(:) = ham%a(:)
      m = ham%nrow
      n = ham%ncol

      !Call MatCreateSeqAIJWithArrays(comm, m, n, i(:), j(:), ham%a(:), mat, ierr)
      Call MatCreateSeqAIJWithArrays(comm, m, n, i(:), j(:), a(:), mat, ierr)
      !Call MatView(mat, PETSC_VIEWER_STDOUT_WORLD, ierr)

      ! Create eigenvectors
      Call VecCreate(PETSC_COMM_SELF, vr, ierr)
      Call VecSetSizes(vr, PETSC_DECIDE, ham%nrow, ierr)
      Call VecSetFromOptions(vr, ierr)
      Call VecDuplicate(vr, vi, ierr)

      ! Create eigensolver context
      Call EPSCreate(PETSC_COMM_WORLD, eps, ierr)

      ! Set operator
      Call EPSSetOperators(eps, mat, PETSC_NULL_MAT, ierr)
      Call EPSSetProblemType(eps, EPS_NHEP, ierr)

      ! Set eigenvalue selection method
      Call EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL, ierr)

      ! Set number of eigenvalues to search for
      nev = k
      Call EPSSetDimensions(eps, nev, PETSC_DEFAULT_INTEGER, PETSC_DEFAULT_INTEGER, ierr)
      
      ! Set (additional) parameters at runtime (i.e. from the command line)
      Call EPSSetFromOptions(eps, ierr)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !           Solve the eigensystem
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Call EPSSolve(eps, ierr)

      ! Extract eigenpairs
      !do ei = 1,k-1
      do ei = 0,k-1
        iev = ei
        Call EPSGetEigenpair(eps, iev, eigr, eigi, vr, vi, ierr)
        !ev(ei+1) = -1.0*cmplx(REAL(eigr),REAL(eigi))
        ev(ei+1) = -1.0*cmplx(REAL(eigr),AIMAG(eigr))
      enddo

      ! cleanup
      Call MatDestroy(mat, ierr)
      Call VecDestroy(vr, ierr)
      Call VecDestroy(vi, ierr)
      Call EPSDestroy(eps, ierr)
      call PetscFinalize(ierr)

      Return
    End Subroutine slepc_solve

  End Module
