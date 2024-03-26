#include "math_def.h"
#include <petsc/finclude/petscksp.h>

  Module petsc_solvers
  Use mpi
  Use petscksp
  Use sparselib
  Implicit None

  Type t_petsc_opt
    Integer :: solver = 1 ! 1) GMRES
    Integer :: preco  = 1 ! 1) ILU 
    Real (Kind=DEF_DBL_PREC) :: tol=1d-8
  End Type t_petsc_opt

  Contains

    Subroutine petsc_solve(ham, rhs, root, comm, opt)
      Implicit None
      Type (t_petsc_opt) :: opt
      Type (zcsrmat)     :: ham
      Type (zdensemat)   :: rhs
      Integer :: root, comm
      !!$ local
      Integer :: commsize, err 
      Integer :: sols 
      Complex (Kind=DEF_DBL_PREC), Pointer :: temp(:)
      Mat   mat
      Vec   x,b
      KSP   ksp
      PC    pc
      PetscErrorCode ierr
      PetscInt          m, n, i(ham%nrow+1), j(ham%nnz)
      PetscScalar       a(ham%nnz)


      ! Set PETSC_COMM_WORLD to subset of MPI_COMM_WORLD
      PETSC_COMM_WORLD = comm

      ! Check comm size: at this time we don't support parallel matrices
      call MPI_COMM_SIZE(comm, commsize, err)
      if (commsize > 1) stop 'We currently do not support parallel matrices'

      ! Initialize PETSc
      Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      If (ierr .ne. 0) then
        stop 'Unable to initialize PETSc'
      endif

      ! create PETSc matrix
      i(:) = ham%ir(:)-1 ! 0-based index for C
      j(:) = ham%jc(:)-1 ! idem
      a(:) = ham%a(:)
      m = ham%nrow
      n = ham%ncol

      Call MatCreateSeqAIJWithArrays(comm, m, n, i(:), j(:), ham%a(:), mat, ierr)
      !Call MatView(mat, PETSC_VIEWER_STDOUT_WORLD, ierr)

      ! create rhs and solution vectors
      Call VecCreate(PETSC_COMM_SELF, x, ierr)
      Call VecSetSizes(x,PETSC_DECIDE, rhs%nrow, ierr)
      Call VecSetFromOptions(x, ierr)
      Call VecDuplicate(x, b, ierr)

      ! create linear solver context
      Call KSPCreate(PETSC_COMM_SELF, ksp, ierr)

      ! set operators. For now, the preconditioner will be based on mat
      Call KSPSetOperators(ksp, mat, mat, ierr)

      ! set preconditioner
      Call KSPGetPC(ksp,pc,ierr)
      If (opt%preco == 1) Then !!$ ILU
        Call PCSetType(pc, PCILU, ierr)
      ElseIf (opt%preco == 2) Then !!$ ASM 
        Call PCSetType(pc, PCASM, ierr)
      Else
        stop 'I do not understand what preconditioner you want to use'
      Endif

      Call KSPSetTolerances(ksp, opt%tol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL,&
                              & PETSC_DEFAULT_INTEGER, ierr)

      ! set solving method
      If (opt%solver == 1) Then !!$ GMRES
        Call KSPSetType(ksp, KSPGMRES, ierr)
      ElseIf (opt%solver == 2) Then !!$ TFQMR
        Call KSPSetType(ksp, KSPTFQMR, ierr)
      ElseIf (opt%solver == 3) Then !!$ BiCGSTAB
        Call KSPSetType(ksp, KSPBICG, ierr)
      Else
        stop 'I do not understand what solver you want to use'
      Endif

      Call KSPSetFromOptions(ksp, ierr)

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !             Solve the linear system of equations
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      Do sols=1,rhs%ncol !solve different RHS
        ! set rhs
        Call VecGetArrayF90(b, temp, ierr)
        temp(:) = rhs%bl(:,sols)
        Call VecRestoreArrayF90(b, temp, ierr)

        ! solve lineq
        Call KSPSolve(ksp, b, x, ierr)

        ! save solution in rhs, this is what the code expects
        Call VecGetArrayReadF90(x, temp, ierr)
        rhs%bl(:,sols) = temp(:)
        Call VecRestoreArrayF90(x, temp, ierr)
      Enddo

      ! cleanup
      Call VecDestroy(x, ierr)
      Call VecDestroy(b, ierr)
      Call MatDestroy(mat, ierr)
      Call KSPDestroy(ksp, ierr)
      call PetscFinalize(ierr)

      Return
    End Subroutine petsc_solve

  End Module
