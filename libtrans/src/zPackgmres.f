        subroutine drive_zgmres(n,nloc,m,lwork,work,
     &                         irc,icntl,cntl,info,rinfo)
*
*  Purpose
*  =======
*    drive_zgmres is the driver routine for solving the linear system 
*  Ax = b using the *  Generalized Minimal Residual iterative method 
*  with preconditioning.
*  This solver is implemented with a reverse communication scheme: control
*  is returned to the user for computing the (preconditioned) 
*  matrix-vector product.
*  See the User's Guide for an example of use.
*
*
* Written : June 1996
* Authors : Luc Giraud, Serge Gratton, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : June 1998
* Authors : Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
* Purpose : Make clear that the warning and error messages come from the
*           zgmres modules.
*
* Updated : December 2002 - L. Giraud, J.Langou
* Purpose : Add the capability to avoid explicit residual calculation at restart
*
* Updated : March 2005 - L. Giraud
* Purpose : small bug in the computation of the restart if too large vs the
*           size of the workspace
*
*  Arguments
*  =========
*
*  n      (input) INTEGER.
*          On entry, the dimension of the problem.
*          Unchanged on exit.
*
*  nloc   (input) INTEGER.
*          On entry, the dimension of the local problem.
*          In a parallel distributed envirionment, this corresponds
*          to the size of the subset of entries of the right hand side
*          and solution allocated to the calling process.
*          Unchanged on exit.
*
*
*  m      (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix H (see WORK and H).
*          Unchanged on exit.
*
*  lwork  (input) INTEGER
*          size of the workspace
*          if (icntl(5) = 0 or 1 )
*            lwork >= m*m + m*(n+5) + 5*n+2, if icntl(8) = 1
*            lwork >= m*m + m*(n+5) + 6*n+2, if icntl(8) = 0
*          if (icntl(5) = 2 or 3 )
*            lwork >= m*m + m*(n+5) + 5*n+m+1, if icntl(8) = 1
*            lwork >= m*m + m*(n+5) + 6*n+m+1, if icntl(8) = 0
*
*  work   (workspace) double precision/double complex array, length lwork
*          work contains the required vector and matrices stored in the 
*          following order :
*            x  (n,1)       : computed solution.
*            b  (n,1)       : right hand side.
*            r0 (n,1)       : vector workspace.
*            w  (n,1)       : vector workspace.
*            V  (n,m)       : Krylov basis.
*            H  (m+1,m+1)   : Hessenberg matrix (full storage).
*            yCurrent (m,1) : solution of the current LS
*            xCurrent (n,1) : current iterate
*            rotSin (m,1)   : Sine of the Givens rotation
*            rotCos (m,1)   : Cosine of the Givens rotation
*
*  irc     (input/output) INTEGER array. length 5
*            irc(1) : REVCOM   used for reverse communication
*                             (type of external operation)
*            irc(2) : COLX     used for reverse communication
*            irc(3) : COLY     used for reverse communication
*            irc(4) : COLZ     used for reverse communication
*            irc(5) : NBSCAL   used for reverse communication
*
*  icntl   (input) INTEGER array. length 7
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       2 - right preconditioning
*                       3 - double side preconditioning
*                       4 - error, default set in Init
*            icntl(5) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(7) : maximum number of iterations
*            icntl(8) : 0 - use recurence formula at restart
*                       1 - default compute the true residual at each restart
*
*  cntl    (input) double precision array, length 5
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*            cntl(4) : scaling factor for normwise perturbation on the
*                      preconditioned matrix
*            cntl(5) : scaling factor for normwise perturbation on 
*                      preconditioned right hand side
*
*  info    (output) INTEGER array, length 3
*            info(1) :  0 - normal exit
*                      -1 - n < 1
*                      -2 - m < 1
*                      -3 - lwork too small
*                      -4 - convergence not achieved after icntl(7) iterations
*                      -5 - precondition type not set by user
*            info(2) : if info(1)=0 - number of iteration to converge
*                      if info(1)=-3 - minimum workspace size necessary
*            info(3) : optimal size for the workspace
*
*  rinfo   (output) double precision array, length 2
*            if info(1)=0 
*              rinfo(1) : backward error for the preconditioned system
*              rinfo(2) : backward error for the unpreconditioned system
*
* Input variables
* ---------------
       integer n, nloc, lwork, icntl(*)
       double precision   cntl(*)
       double precision   sA, sb, sPA, sPb
* Output variables
* ----------------
       integer  info(*)
       double precision    rinfo(*)
* Input/Output variables
* ----------------------
       integer  m, irc(*)
       double complex work(*)
* Local variables
* ---------------
       integer xptr, bptr, wptr, r0ptr, Vptr, Hptr, dotptr
       integer yCurrent,rotSin, rotCos, xCurrent
       integer sizeWrk, newRestart
       integer iwarn, ierr, ihist, compRsd
       real    rn,rx,rc
       double precision DZRO
       parameter (DZRO = 0.0d0)
*
       integer icheck
       save icheck
       DATA icheck /0/
*
       intrinsic ifix, float
*
*       Executable statements :
*
       ierr    = icntl(1)
       iwarn   = icntl(2)
       ihist   = icntl(3)
       compRsd = icntl(8)
*
       if (ierr.lt.0) ierr = 6
*
       if (compRsd.eq.1) then
          sizeWrk  = m*m + m*(nloc+5) + 5*nloc+ 1
       else
          sizeWrk  = m*m + m*(nloc+5) + 6*nloc+ 1
       endif
*
       if (icheck.eq.0) then
* Check the value of the arguments
         if ((n.lt.1).or.(nloc.lt.1)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     N < 1 '
            write(ierr,*)
            info(1) = -1
            irc(1)  = 0
            return
         endif
         if (m.lt.1) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES :'
            write(ierr,*)'     M < 1 '
            write(ierr,*)
            info(1) = -2
            irc(1)  = 0
            return
         endif
         if ((icntl(4).ne.0).and.(icntl(4).ne.1).and.
     &     (icntl(4).ne.2).and.(icntl(4).ne.3)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     Undefined preconditioner '
            write(ierr,*)
            info(1) = -5
            irc(1)  = 0
            return
         endif
*
         if ((icntl(5).lt.0).or.(icntl(5).gt.3)) then
           icntl(5) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING  GMRES : '
             write(iwarn,*) '       Undefined orthogonalisation '  
             write(iwarn,*) '       Default MGS '
             write(iwarn,*)
           endif
         endif
*
         if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
* the workspace should be large enough to store the m dot-products
            sizeWrk  = sizeWrk  + m
         else
            sizeWrk  = sizeWrk  + 1
         endif
*
C          if (iwarn.ne.0) then
C            write(iwarn,*)
C            write(iwarn,*) ' WARNING GMRES : '
C            write(iwarn,*) '       For M = ',m,' optimal value '     
C            write(iwarn,*) '       for LWORK =  ', sizeWrk
C            write(iwarn,*)
C          endif
*
         if ((icntl(6).ne.0).and.(icntl(6).ne.1)) then
           icntl(6) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING GMRES : '
             write(iwarn,*) '       Undefined intial guess '
             write(iwarn,*) '       Default x0 = 0 '
             write(iwarn,*)
           endif
         endif
         if (icntl(7).le.0) then
           icntl(7) = n
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING GMRES :'
             write(iwarn,*) '       Negative max number of iterations'  
             write(iwarn,*) '       Default N '
             write(iwarn,*)
           endif
         endif
         if ((icntl(8).ne.0).and.(icntl(8).ne.1)) then
           icntl(8) = 1
           write(iwarn,*)
           write(iwarn,*) ' WARNING GMRES :'
           write(iwarn,*) '       Undefined strategy for the residual'
           write(iwarn,*) '       at restart'
           write(iwarn,*) '       Default 1 '
           write(iwarn,*)
         endif
* Check if the restart parameter is correct and if the size of the
*  workspace is big enough for the restart.
* If not try to fix correctly the parameters
*
         if ((m .gt. n).or.(lwork.lt.sizeWrk)) then
           if (m .gt. n) then
             m = n
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING GMRES : '
               write(iwarn,*) '       Parameter M bigger than N'  
               write(iwarn,*) '       New value for M ',m
               write(iwarn,*)
             endif
             if (compRsd.eq.1) then
               sizeWrk = m*m + m*(nloc+5) + 5*nloc+1
             else
               sizeWrk = m*m + m*(nloc+5) + 6*nloc+1
             endif
             if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
* the workspace should be large enough to store the m dot-products
                sizeWrk  = sizeWrk  + m
             else
                sizeWrk  = sizeWrk  + 1
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.eq.nloc)) then
* Compute the maximum size of the m according to the memory space
             rn         = float(n)
             rx         = rn + 5.0
             rc         = 5.0*rn + 1 - float(lwork)
*
* Update the linear part of the second order equation to be solved
             if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
               rx = rx + 1
             endif
* Update the constant part of the second order equation to be solved
*             
             if (icntl(8).eq.0) then
               rc = rc + rn
             endif
* Solution of the the second order equation
             newRestart = ifix((-rx+sqrt(rx**2-4.0*rc))/2.0)
             if (newRestart.gt.0) then
               m = newRestart
               if (iwarn.ne.0) then
                 write(iwarn,*)
                 write(iwarn,*)' WARNING GMRES : '
                 write(iwarn,*)'       Workspace too small for M'  
                 write(iwarn,*)'       New value for M ',m
                 write(iwarn,*)
               endif
             else
               write(ierr,*)
               write(ierr,*)' ERROR GMRES : '
               write(ierr,*)'     Not enough space for the problem'
               write(ierr,*)'     the space does not permit any m'
               write(ierr,*)
               info(1) = -3
               irc(1)  = 0
               return
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.ne.nloc)) then
              write(ierr,*)
              write(ierr,*)' ERROR GMRES : '
              write(ierr,*)'     Not enough space for the problem'
              write(ierr,*)
              info(1) = -3
              irc(1)  = 0
              return
           endif
         endif
*
         info(3) = sizeWrk
         icheck = 1
*
* save the parameters the the history file
*
         if (ihist.ne.0) then
           write(ihist,'(10x,A39)') 'CONVERGENCE HISTORY FOR GMRES'
           write(ihist,*)
           write(ihist,'(A30,I2)') 'Errors are displayed in unit: ',ierr
           if (iwarn.eq.0) then
             write(ihist,'(A27)') 'Warnings are not displayed:'
           else
             write(ihist,'(A32,I2)') 'Warnings are displayed in unit: ',
     &                               iwarn
           endif 
           write(ihist,'(A13,I7)') 'Matrix size: ',n
           write(ihist,'(A19,I7)') 'Local matrix size: ',nloc
           write(ihist,'(A9,I7)') 'Restart: ',m
           if (icntl(4).eq.0) then
             write(ihist,'(A18)') 'No preconditioning'
           elseif (icntl(4).eq.1) then
             write(ihist,'(A20)') 'Left preconditioning'
           elseif (icntl(4).eq.2) then
             write(ihist,'(A21)') 'Right preconditioning'
           elseif (icntl(4).eq.3) then
             write(ihist,'(A30)') 'Left and right preconditioning'
           endif
           if (icntl(5).eq.0) then
             write(ihist,'(A21)') 'Modified Gram-Schmidt'
           elseif (icntl(5).eq.1) then
             write(ihist,'(A31)') 'Iterative modified Gram-Schmidt'
           elseif (icntl(5).eq.2) then
             write(ihist,'(A22)') 'Classical Gram-Schmidt'
           else
             write(ihist,'(A32)') 'Iterative classical Gram-Schmidt'
           endif
           if (icntl(6).eq.0) then
             write(ihist,'(A29)') 'Default initial guess x_0 = 0'
           else
             write(ihist,'(A27)') 'User supplied initial guess'
           endif
           if (icntl(8).eq.1) then
             write(ihist,'(A33)') 'True residual computed at restart'
           else
             write(ihist,'(A30)') 'Recurrence residual at restart'
           endif
           write(ihist,'(A30,I5)') 'Maximum number of iterations: ',
     &                              icntl(7)
           write(ihist,'(A27,E8.2)') 'Tolerance for convergence: ', 
     &                                cntl(1) 
* 
           write(ihist,'(A53)') 
     &       'Backward error on the unpreconditioned system Ax = b:'
           sA       = cntl(2)
           sb       = cntl(3)
           if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
             write(ihist,'(A39)') 
     &       '    the residual is normalised by ||b||'
           else
             write(ihist,1) sA,sb
1     format('    the residual is normalised by         ',E8.2,
     &       ' * ||x|| + ',E8.2)
           endif
           sPA      = cntl(4)
           sPb      = cntl(5)
             write(ihist,2)
2     format('Backward error on the preconditioned system',
     &        ' (P1)A(P2)y = (P1)b:')
           if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
             write(ihist,3)
3     format('    the preconditioned residual is normalised ',
     &       'by ||(P1)b||')
           else
             write(ihist,4) sPA, sPb
4     format('    the preconditioned residual is normalised by ', E8.2,
     &        ' * ||(P2)y|| + ',E8.2)
           endif
*
           write(ihist,5) info(3)
5     format('Optimal size for the local workspace:',I7)
           write(ihist,*) 
           write(ihist,6)
6     format('Convergence history: b.e. on the preconditioned system')
           write(ihist,7)
7     format(' Iteration   Arnoldi b.e.    True b.e.')
         endif
*
       endif
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + nloc
       r0ptr    = bptr + nloc
       wptr     = r0ptr + nloc
       Vptr     = wptr + nloc
       if (compRsd.eq.1) then
          Hptr     = Vptr + m*nloc
       else
          Hptr     = Vptr + (m+1)*nloc
       endif
       dotptr   = Hptr + (m+1)*(m+1)
       if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
          yCurrent = dotptr + m
       else
         yCurrent = dotptr + 1
       endif
       xCurrent = yCurrent + m
       rotSin   = xCurrent + nloc
       rotCos   = rotSin + m
*
       call zgmres(nloc,m,work(bptr),work(xptr),
     &            work(Hptr),work(wptr),work(r0ptr),work(Vptr),
     &            work(dotptr),work(yCurrent),work(xCurrent),
     &            work(rotSin),work(rotCos),irc,icntl,cntl,info,rinfo)
*
       if (irc(1).eq.0) then
         icheck = 0
       endif
*
       return
       end
*
        subroutine zgmres(n,m,b,x,H,w,r0,V,dot,yCurrent,xCurrent,rotSin,
     &                   rotCos,irc,icntl,cntl,info,rinfo)
*
*
*  Purpose
*  =======
*  zgmres solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method
*
* When preconditioning is used we solve :
*     M_1^{-1} A M_2^{-1} y = M_1^{-1} b
*     x = M_2^{-1} y
*
*   Convergence test based on the normwise backward error for
*  the preconditioned system
*
* Written : June 1996
* Authors : Luc Giraud, Serge Gratton, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : March 1998
* Purpose : Pb with F90 on DEC ws
*           cure : remove "ZDSCAL" when used to initialize vectors to zero
*
* Updated : May 1998
* Purpose : r0(1) <-- r0'r0 : pb when used with DGEMV for the dot product
*           cure : w(1) <--  r0'r0
*
* Updated : June 1998
* Purpose : Make clear that the warning and error messages come from the
*           zgmres modules.
*
* Updated : February 2001 - L. Giraud
* Purpose : In complex version, initializations to zero performed  in complex
*           arithmetic to avoid implicit conversion by the compiler.
*
* Updated : July 2001 - L. Giraud, J. Langou
* Purpose : Avoid to compute the approximate solution at each step of
*           the Krylov space construction when spA is zero.
*
* Updated : November 2002 - S. Gratton
* Purpose : Use Givens rotations conform to the classical definition.
*           No impact one the convergence history.
*
* Updated : November 2002 - L. Giraud
* Purpose : Properly handle the situation when the convergence is obtained
*           exactly at the "IterMax" iteration
*
* Updated : December 2002 - L. Giraud, J.Langou
* Purpose : Add the capability to avoid explicit residual calculation at restart
*
* Updated : January  2003 - L. Giraud, S. Gratton
* Purpose : Use Givens rotations from BLAS.
*
* Updated : March    2003 - L. Giraud
* Purpose : Set back retlbl to zero, if initial guess is solution
*           or right-hand side is zero
*
* Updated : September 2003 - L. Giraud
* Purpose : Include room in the workspace to store the results of the dot products
*           Fix the bugs that appeared when M > Nloc 
*
*  Arguments
*  =========
*
*  n       (input) INTEGER.
*           On entry, the dimension of the problem.
*           Unchanged on exit.
*
*  m        (input) INTEGER
*           Restart parameter, <= N. This parameter controls the amount
*           of memory required for matrix H (see WORK and H).
*           Unchanged on exit.
*
*  b        (input) double precision/double complex
*           Right hand side of the linear system.
*
*  x        (output) double precision/double complex
*           Computed solution of the linear system.
*
*  H        (workspace)  double precision/double complex
*           Hessenberg matrix built within dgmres
*
*  w        (workspace)  double precision/double complex
*           Vector used as temporary storage
*
*  r0       (workspace)  double precision/double complex
*           Vector used as temporary storage
*
*  V        (workspace)  double precision/double complex
*           Basis computed by the Arnoldi's procedure.
*  
*  dot      (workspace) double precision/double complex
*           Store the results of the dot product calculation
*
*  yCurrent (workspace) double precision/double complex
*           solution of the current LS
*
*  xCurrent (workspace) double precision/double complex
*           current iterate
*
*  rotSin   (workspace) double precision/double complex
*           Sine of the Givens rotation
*
*  rotCos   (workspace) double precision
*           Cosine of the Givens rotation
*
*  irc      (input/output) INTEGER array. length 3
*             irc(1) : REVCOM   used for reverse communication
*                              (type of external operation)
*             irc(2) : COLX     used for reverse communication
*             irc(3) : COLY     used for reverse communication
*             irc(4) : COLZ     used for reverse communication
*             irc(5) : NBSCAL   used for reverse communication
*
*  icntl    (input) INTEGER array. length 7
*             icntl(1) : stdout for error messages
*             icntl(2) : stdout for warnings
*             icntl(3) : stdout for convergence history
*             icntl(4) : 0 - no preconditioning
*                        1 - left preconditioning
*                        2 - right preconditioning
*                        3 - double side preconditioning
*                        4 - error, default set in Init
*             icntl(5) : 0 - modified Gram-Schmidt
*                        1 - iterative modified Gram-Schmidt
*                        2 - classical Gram-Schmidt
*                        3 - iterative classical Gram-Schmidt
*             icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
*                        1 - user supplied initial guess
*             icntl(7) : maximum number of iterations
*             icntl(8) : 1 - default compute the true residual at each restart
*                        0 - use recurence formula at restart
*
*  cntl     (input) double precision array, length 5
*             cntl(1) : tolerance for convergence
*             cntl(2) : scaling factor for normwise perturbation on A
*             cntl(3) : scaling factor for normwise perturbation on b
*             cntl(4) : scaling factor for normwise perturbation on the
*                       preconditioned matrix
*             cntl(5) : scaling factor for normwise perturbation on
*                       preconditioned right hand side
*
*  info     (output) INTEGER array, length 2
*             info(1) :  0 - normal exit
*                       -1 - n < 1
*                       -2 - m < 1
*                       -3 - lwork too small
*                       -4 - convergence not achieved after icntl(7) iterations
*                       -5 - precondition type not set by user
*             info(2) : if info(1)=0 - number of iteration to converge
*                       if info(1)=-3 - minimum workspace size necessary
*             info(3) : optimal size for the workspace
*
* rinfo     (output) double precision array, length 2
*             if info(1)=0 
*               rinfo(1) : backward error for the preconditioned system
*               rinfo(2) : backward error for the unpreconditioned system
*
* Input variables
* ---------------
        integer  n, m, icntl(*)
        double complex b(*)
        double precision    cntl(*)
*
* Output variables
* ----------------
       integer  info(*)
       double precision    rinfo(*)
*
* Input/Output variables
* ----------------------
       integer  irc(*)
       double complex x(*), H(m+1,*), w(*), r0(*) 
       double complex V(n,*),dot(*),yCurrent(*)  
       double complex xCurrent(*), rotSin(*)
       double complex rotCos(*)
*
* Local variables
* ---------------
       integer  j, jH, iterOut, nOrtho, iterMax, initGuess, iOrthog
       integer  xptr, bptr, wptr, r0ptr, Vptr, Hptr, yptr, xcuptr
       integer  dotptr
       integer  typePrec, leftPrec, rightPrec, dblePrec, noPrec
       integer  iwarn, ihist
       integer  compRsd
       double precision    beta, bn, sA, sb, sPA, sPb, bea, be, temp
       double precision    dloo, dnormw, dnormx, dnormres, trueNormRes
       double complex dVi, aux
       double complex auxHjj, auxHjp1j
*
       parameter (noPrec = 0, leftPrec = 1)
       parameter (rightPrec = 2, dblePrec = 3)  
*
       double complex ZERO, ONE
       parameter (ZERO = (0.0d0, 0.0d0), ONE = (1.0d0, 0.0d0))
       double precision DZRO,DONE
       parameter (DZRO = 0.0d0, DONE = 1.0d0)
*
*
* External functions
* ------------------
       double precision    dznrm2
       external dznrm2
*
* Saved variables
* ---------------
       save iterOut, jH, beta, bn, dnormres, retlbl, j
       save sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be
       save dloo, nOrtho, compRsd
*
* Intrinsic function
* ------------------
       intrinsic abs, dsqrt, dreal, dconjg
*
*
* Reverse communication variables
* -------------------------------
       integer matvec, precondLeft, precondRight, prosca
       parameter(matvec=1, precondLeft=2, precondRight=3, prosca=4)
       integer retlbl
       DATA retlbl /0/
*
*       Executable statements
*
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + n
       r0ptr    = bptr + n
       wptr     = r0ptr + n
       Vptr     = wptr + n
       if (icntl(8).eq.1) then
         Hptr     = Vptr + m*n
       else
         Hptr     = Vptr + (m+1)*n
       endif
       dotptr   = Hptr + (m+1)*(m+1)
       if ((icntl(5).eq.2).or.(icntl(5).eq.3)) then
         yptr = dotptr + m
       else
         yptr = dotptr + 1
       endif
       xcuptr   = yptr + m
*
       iwarn      = icntl(2)
       ihist      = icntl(3)
       typePrec   = icntl(4)
       iOrthog    = icntl(5)
       initGuess  = icntl(6)
       iterMax    = icntl(7)
*
       if (retlbl.eq.0) then
         compRsd    = icntl(8)
       endif
*
      if (retlbl.ne.0) then
        if (retlbl.eq.5) then
          goto 5
        else if (retlbl.eq.6) then
          goto 6
        else if (retlbl.eq.8) then
          goto 8
        else if (retlbl.eq.11) then
          goto 11
        else if (retlbl.eq.16) then
          goto 16
        else if (retlbl.eq.18) then
          goto 18
        else if (retlbl.eq.21) then
          goto 21
        else if (retlbl.eq.26) then
          goto 26
        else if (retlbl.eq.31) then
          goto 31
        else if (retlbl.eq.32) then
          goto 32
        else if (retlbl.eq.33) then
          goto 33
        else if (retlbl.eq.34) then
          goto 34 
        else if (retlbl.eq.36) then
          goto 36
        else if (retlbl.eq.37) then
          goto 37
        else if (retlbl.eq.38) then
          goto 38
        else if (retlbl.eq.41) then
          goto 41
        else if (retlbl.eq.43) then
          goto 43
        else if (retlbl.eq.46) then
          goto 46
        else if (retlbl.eq.48) then
          goto 48
        else if (retlbl.eq.51) then
          goto 51
        else if (retlbl.eq.52) then
          goto 52
        else if (retlbl.eq.61) then
          goto 61
        else if (retlbl.eq.66) then
          goto 66
        else if (retlbl.eq.68) then
          goto 68
        endif
      endif
*
*
* intialization of various variables
*
      iterOut  = 0
      beta     = DZRO
*
      if (initGuess.eq.0) then
        do j=1,n
          x(j) = ZERO
        enddo
      endif
*
*        bn = dznrm2(n,b,1)
*
      irc(1) = prosca
      irc(2) = bptr
      irc(3) = bptr
      irc(4) = dotptr
      irc(5) = 1
      retlbl = 5
      return
 5    continue
      bn = dsqrt(dreal(dot(1)))
*
      if (bn.eq.DZRO) then
        do j=1,n
          x(j) = ZERO
        enddo  
        if (iwarn.ne.0) then
          write(iwarn,*)
          write(iwarn,*) ' WARNING GMRES : '
          write(iwarn,*) '       Null right hand side'
          write(iwarn,*) '       solution set to zero'
          write(iwarn,*)
        endif
        jH = 0
        bea = DZRO
        be  = DZRO
        write(ihist,'(I5,11x,E8.2,$)') jH,bea
        write(ihist,'(7x,E8.2)') be
        info(1)  = 0
        info(2)  = 0
        rinfo(1) = DZRO
        rinfo(2) = DZRO
        irc(1)   = 0
        retlbl = 0
        return
      endif
*
* Compute the scaling factor for the backward error on the 
*  unpreconditioned sytem
*
      sA       = cntl(2)
      sb       = cntl(3)
      if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
        sb = bn
      endif
* Compute the scaling factor for the backward error on the
*  preconditioned sytem
*
       sPA      = cntl(4)
       sPb      = cntl(5)
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
           sPb = bn
         else
           irc(1) = precondLeft
           irc(2) = bptr
           irc(4) = r0ptr
           retlbl = 6
           return
         endif
       endif
 6     continue
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.dblePrec).or.(typePrec.eq.leftPrec)) then
*
*           sPb = dznrm2(n,r0,1)
*
           irc(1) = prosca
           irc(2) = r0ptr
           irc(3) = r0ptr
           irc(4) = dotptr
           irc(5) = 1
           retlbl = 8
           return
         endif
       endif
 8     continue
       if ((sPA.eq.DZRO).and.(sPb.eq.DZRO)) then
         if ((typePrec.eq.dblePrec).or.(typePrec.eq.leftPrec)) then
           sPb = dsqrt(dreal(dot(1)))
* 
         endif
       endif
*
*
* Compute the first residual
*           Y = AX : r0 <-- A x
*
* The residual is computed only if the initial guess is not zero
*
       if (initGuess.ne.0) then
         irc(1) = matvec
         irc(2) = xptr
         irc(4) = r0ptr
         retlbl = 11
         return
       endif
 11    continue
       if (initGuess.ne.0) then
         do j=1,n
           r0(j) = b(j)-r0(j)
         enddo
       else
         call zcopy(n,b,1,r0,1)
       endif 
*
* Compute the preconditioned residual if necessary
*      M_1Y = X : w <-- M_1^{-1} r0
*
       if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
         call zcopy(n,r0,1,w,1)
       else
         irc(1) = precondLeft
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 16
         return
       endif
 16    continue
*
*
*       beta = dznrm2(n,w,1)
*
*
       irc(1) = prosca
       irc(2) = wptr
       irc(3) = wptr
       irc(4) = dotptr
       irc(5) = 1
       retlbl = 18
       return
 18    continue
       beta = dsqrt(dreal(dot(1)))
*
       if (beta .eq. DZRO) then
*  The residual is exactly zero : x is the exact solution
         info(1) = 0
         info(2) = 0
         rinfo(1) = DZRO
         rinfo(2) = DZRO
         irc(1)   = 0
         retlbl = 0
         jH = 0
         bea = DZRO
         be  = DZRO
         write(ihist,'(I5,11x,E8.2,$)') jH,bea
         write(ihist,'(7x,E8.2)') be
         if (iwarn.ne.0) then
          write(iwarn,*)
          write(iwarn,*) ' WARNING GMRES : '
          write(iwarn,*) '       Intial residual is zero'
          write(iwarn,*) '       initial guess is solution'
          write(iwarn,*)
         endif
         return
       endif
*
       aux = ONE/beta
       do j=1,n
         V(j,1) = ZERO
       enddo
       call zaxpy(n,aux,w,1,V(1,1),1)
*
*       Most outer loop : zgmres iteration
*
*       REPEAT
 7     continue
*
*
       H(1,m+1)=beta
       do j=1,m
         H(j+1,m+1) = ZERO
       enddo
*
*        Construction of the hessenberg matrix WORK and of the orthogonal
*        basis V such that AV=VH 
*
       jH = 1
 10    continue
* Remark : this  do loop has been written with a while do
*          because the
*               " do jH=1,restart "
*         fails with the reverse communication.
*      do  jH=1,restart
*
*
* Compute the preconditioned residual if necessary
*
       if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then  
*
*           Y = M_2^{-1}X : w <-- M_2^{-1} V(1,jH)
*
         irc(1) = precondRight
         irc(2) = vptr + (jH-1)*n
         irc(4) = wptr
         retlbl = 21
         return
       else
         call zcopy(n,V(1,jH),1,w,1)
       endif
 21    continue
*
*           Y = AX : r0 <-- A w
*
       irc(1) = matvec
       irc(2) = wptr
       irc(4) = r0ptr
       retlbl = 26
       return
 26    continue
*
*      MY = X : w <-- M_1^{-1} r0
*
       if ((typePrec.eq.noPrec).or.(typePrec.eq.rightPrec)) then
         call zcopy(n,r0,1,w,1)
       else
         irc(1) = precondLeft
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 31
         return
       endif
 31    continue
*
* Orthogonalization using either MGS or IMGS
*  
* initialize the Hessenberg matrix to zero in order to be able to use
*     IMGS as orthogonalization procedure.
       do j=1,jH
         H(j,jH) = ZERO
       enddo
       nOrtho = 0
 19    continue
       nOrtho = nOrtho +1
       dloo   = DZRO
*
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
* MGS
*
*           do j=1,jH
*
         j = 1
*           REPEAT
       endif
 23    continue
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
*
*             dVi     = zdotc(n,V(1,j),1,w,1)
*
         irc(1) = prosca
         irc(2) = vptr + (j-1)*n
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 32
         return
       endif
 32    continue
       if ((iOrthog.eq.0).or.(iOrthog.eq.1)) then
         dVi     = dot(1)
         H(j,jH) = H(j,jH) + dVi
         dloo    = dloo + abs(dVi)**2
         aux = -ONE*dVi
         call zaxpy(n,aux,V(1,j),1,w,1)
         j = j + 1
         if (j.le.jH) goto 23
*          enddo_j
       else
* CGS
* Gathered dot product calculation
*
*           call zgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1)
*
         irc(1) = prosca
         irc(2) = vptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = jH
         retlbl = 34
         return
       endif
 34    continue
       if ((iOrthog.eq.2).or.(iOrthog.eq.3)) then
*
         call zaxpy(jH,ONE,dot,1,H(1,jH),1)
         call zgemv('N',n,jH,-ONE,V(1,1),n,dot,1,ONE,w,1)
         dloo = dznrm2(jH,dot,1)**2
       endif
*
*         dnormw = dznrm2(n,w,1)
*
       irc(1) = prosca
       irc(2) = wptr
       irc(3) = wptr
       irc(4) = dotptr
       irc(5) = 1
       retlbl = 33
       return
 33    continue
       dnormw = dsqrt(dreal(dot(1)))
*
       if ((iOrthog.eq.1).or.(iOrthog.eq.3)) then
* IMGS / CGS orthogonalisation
         dloo = dsqrt(dloo)
* check the orthogonalization quality
         if (((2.0*dnormw).le.dloo).and.(nOrtho.lt.3)) then
           goto 19
         endif
       endif
*
       H(jH+1,jH) = dnormw
       if ((jH.lt.m).or.(icntl(8).eq.0)) then
         aux = ONE/dnormw
         do j=1,n
           V(j,jH+1) = ZERO
         enddo
         call zaxpy(n,aux,w,1,V(1,jH+1),1)
       endif
* Apply previous Givens rotations to the new column of H
       do j=1,jH-1
         call zrot(1, H(j,jH), 1, H(j+1,jH), 1,dreal(rotCos(j)), 
     &             rotSin(j))
       enddo
       auxHjj = H(jH,jH)
       auxHjp1j= H(jH+1,jH)
       call zrotg(auxHjj, auxHjp1j,temp,rotSin(jH))
       rotCos(jH)= dcmplx(temp, DZRO)
* Apply current rotation to the rhs of the least squares problem
       call zrot(1, H(jH,m+1), 1, H(jH+1,m+1), 1, dreal(rotCos(jH)),
     &            rotSin(jH))
*
* zabs(H(jH+1,m+1)) is the residual computed using the least squares
*          solver
* Complete the QR factorisation of the Hessenberg matrix by apply the current
* rotation to the last entry of the collumn
       call zrot(1, H(jH,jH), 1, H(jH+1,jH), 1, dreal(rotCos(jH)),
     &           rotSin(jH))
       H(jH+1,jH) = ZERO
*
* Get the Least square residual
*
       dnormres = abs(H(jH+1,m+1))
       if (sPa.ne.DZRO) then
*
* Compute the solution of the current linear least squares problem
*
         call zcopy(jH,H(1,m+1),1,yCurrent,1)
         call ztrsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Compute the value of the new iterate 
*
         call zgemv('N',n,jH,ONE,v,n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
         if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then  
*
*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
*
           irc(1) = precondRight
           irc(2) = xcuptr
           irc(4) = r0ptr
           retlbl = 36
           return
         else
           call zcopy(n,xCurrent,1,r0,1)
         endif
       endif
 36    continue
*
*
       if (sPa.ne.DZRO) then
* Update the current solution
         call zcopy(n,x,1,xCurrent,1)
         call zaxpy(n,ONE,r0,1,xCurrent,1)
*
*         dnormx = dznrm2(n,xCurrent,1)
*
         irc(1) = prosca
         irc(2) = xcuptr
         irc(3) = xcuptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 38
         return
       else
         dnormx    = DONE
       endif
 38    continue
       if (sPa.ne.DZRO) then
         dnormx = dsqrt(dreal(dot(1)))
       endif
*
       bea = dnormres/(sPA*dnormx+sPb)
*
* Check the convergence based on the Arnoldi Backward error for the
* preconditioned system
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then  
* 
* The Arnoldi Backward error indicates that zgmres might have converge
* enforce the calculation of the true residual at next restart
         compRsd = 1
*
*  If the update of X has not yet been performed
         if (sPA.eq.DZRO) then
*
* Compute the solution of the current linear least squares problem
*
           call zcopy(jH,H(1,m+1),1,yCurrent,1)
           call ztrsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Compute the value of the new iterate 
*
           call zgemv('N',n,jH,ONE,v,n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
           if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then
* 
*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
*
             irc(1) = precondRight
             irc(2) = xcuptr
             irc(4) = r0ptr 
             retlbl = 37
             return
           else
             call zcopy(n,xCurrent,1,r0,1)
           endif
         endif
       endif
 37    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         if (sPA.eq.DZRO) then
* Update the current solution
            call zcopy(n,x,1,xCurrent,1)
            call zaxpy(n,ONE,r0,1,xCurrent,1)
         endif
*
         call zcopy(n,xCurrent,1,r0,1)
* Compute the true residual, the Arnoldi one may be unaccurate
*
*           Y = AX : w  <-- A r0
*
         irc(1) = matvec
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 41
         return
       endif
 41    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
*
         do j=1,n
           w(j) = b(j) - w(j)
         enddo
* Compute the norm of the unpreconditioned residual
*
*        trueNormRes = dznrm2(n,w,1)
*
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 43
         return
       endif
 43    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         trueNormRes = dsqrt(dreal(dot(1)))
*
         if ((typePrec.eq.leftPrec).or.(typePrec.eq.dblePrec)) then  
*
*      MY = X : r0 <-- M_1^{-1} w 
*
           irc(1) = precondLeft
           irc(2) = wptr
           irc(4) = r0ptr
           retlbl = 46
           return
         else 
           call zcopy(n,w,1,r0,1)
         endif
       endif
 46    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
*
*        dnormres = dznrm2(n,r0,1) 
*
         irc(1) = prosca
         irc(2) = r0ptr
         irc(3) = r0ptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 48
         return
       endif
 48    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         dnormres = dsqrt(dreal(dot(1)))
*
         be = dnormres/(sPA*dnormx+sPb)
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,1000)iterOut*m+jH,bea,be
1000  format(I5,11x,E8.2,7x,E8.2)
         endif
*
       endif
*
*
* Check again the convergence
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
         if ((be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
* The convergence has been achieved, we restore the solution in x
* and compute the two backward errors.
           call zcopy(n,xCurrent,1,x,1)
*
           if (sA.ne.DZRO) then
*
*            dnormx = dznrm2(n,x,1)
*
             irc(1) = prosca
             irc(2) = xptr
             irc(3) = xptr
             irc(4) = dotptr
             irc(5) = 1
             retlbl = 51
             return
           endif
         endif
       endif
 51    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         if ((be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
           if (sA.ne.DZRO) then
             dnormx = dsqrt(dreal(dot(1)))
*
           else
             dnormx = DONE
           endif
* Return the backward errors
           rinfo(1) = be
           rinfo(2) = trueNormRes/(sA*dnormx+sb)
           if (be.le.cntl(1)) then
             info(1) = 0
             if (ihist.ne.0) then
               write(ihist,*)
               write(ihist,'(A20)') 'Convergence achieved'
             endif
           else if (be.gt.cntl(1)) then
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING GMRES : '
               write(iwarn,*) '       No convergence after '
               write(iwarn,*) iterOut*m+jH,' iterations '
               write(iwarn,*)
             endif
             if (ihist.ne.0) then
               write(ihist,*)
               write(ihist,*) ' WARNING GMRES :'
               write(ihist,*) '       No convergence after '
               write(ihist,*) iterOut*m+jH,' iterations '
               write(ihist,*)
             endif
             info(1) = -4
           endif
           if (ihist.ne.0) then
             write(ihist,1010) rinfo(1)
             write(ihist,1011) rinfo(2)
1010  format('B.E. on the preconditioned system:   ',E8.2)
1011  format('B.E. on the unpreconditioned system: ',E8.2)
           endif
           info(2) = iterOut*m+jH
           if (ihist.ne.0) then
             write(ihist,'(A10,I2)') 'info(1) = ',info(1)
             write(ihist,'(A32,I5)') 
     &                'Number of iterations (info(2)): ',info(2)  
           endif
           irc(1)  = 0
           retlbl  = 0
           return
         endif
       else
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,1001)iterOut*m+jH,bea
1001  format(I5,11x,E8.2,7x,'--')
         endif
*
       endif  
*
       jH = jH + 1
       if (jH.le.m) then
         goto 10
       endif
*
       iterOut = iterOut + 1
*
* we have completed the Krylov space construction, we restart if
* we have not yet exceeded the maximum number of iterations allowed.
*
       if ((sPa.eq.DZRO).and.(bea.gt.cntl(1))) then
*
* Compute the solution of the current linear least squares problem
*
         jH = jH - 1
         call zcopy(jH,H(1,m+1),1,yCurrent,1)
         call ztrsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Compute the value of the new iterate
*
         call zgemv('N',n,jH,ONE,v,n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
         if ((typePrec.eq.rightPrec).or.(typePrec.eq.dblePrec)) then
*
*         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent
*
           irc(1) = precondRight
           irc(2) = xcuptr
           irc(4) = r0ptr
           retlbl = 52
           return
         else
           call zcopy(n,xCurrent,1,r0,1)
         endif
       endif
 52    continue
       if ((sPa.eq.DZRO).and.(bea.gt.cntl(1))) then
* Update the current solution
         call zcopy(n,x,1,xCurrent,1)
         call zaxpy(n,ONE,r0,1,xCurrent,1)
       endif
*
       call zcopy(n,xCurrent,1,x,1)
*
       if (compRsd.eq.1) then
*
* Compute the true residual
*
         call zcopy(n,x,1,w,1)
         irc(1) = matvec
         irc(2) = wptr
         irc(4) = r0ptr
         retlbl = 61
         return
       endif
 61    continue
       if (compRsd.eq.1) then
         do j=1,n
           r0(j) = b(j) - r0(j)
         enddo
*
* Precondition the new residual if necessary
*
         if ((typePrec.eq.leftPrec).or.(typePrec.eq.dblePrec)) then
*
*      MY = X : w <-- M_1^{-1} r0
*
           irc(1) = precondLeft
           irc(2) = r0ptr
           irc(4) = wptr
           retlbl = 66
           return
         else
           call zcopy(n,r0,1,w,1)
         endif
       endif
 66    continue
*
*           beta = dznrm2(n,w,1)
*
       if (compRsd.eq.1) then
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 68
         return
       endif
 68    continue
       if (compRsd.eq.1) then
         beta = dsqrt(dreal(dot(1)))
*
       else
* Use recurrence to approximate the residual at restart
         beta = abs(H(m+1,m+1))
* Apply the Givens rotation is the reverse order
         do j=m,1,-1
           H(j,m+1)   = ZERO
           call zrot(1, H(j,m+1), 1, H(j+1,m+1), 1,
     &               dreal(rotCos(j)), -rotSin(j))
         enddo
*
* On applique les vecteurs V
*
         call zgemv('N',n,m+1,ONE,v,n,H(1,m+1),1,ZERO,w,1)
*
       endif
       do j=1,n
         V(j,1) = ZERO
       enddo
       aux = ONE/beta
       call zaxpy(n,aux,w,1,V(1,1),1)
*
       goto 7
*
       end
*
*
        subroutine init_zgmres(icntl,cntl)
*
*  Purpose
*  =======
*    Set default values for the parameters defining the characteristics
* of the Gmres algorithm.
*  See the User's Guide for an example of use.
*
*
* Written : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
*
*  Arguments
*  =========
*
* icntl    (input) INTEGER array. length 6
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - no preconditioning
*                       1 - left preconditioning
*                       2 - right preconditioning
*                       3 - double side preconditioning
*                       4 - error, default set in Init
*            icntl(5) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(7) : maximum number of iterations
*            icntl(8) : 1 - default compute the true residual at each restart
*                       0 - use recurence formaula at restart
*
* cntl     (input) double precision array, length 5
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*            cntl(4) : scaling factor for normwise perturbation on the
*                      preconditioned matrix
*            cntl(5) : scaling factor for normwise perturbation on 
*                      preconditioned right hand side
*
* Output variables
* ----------------
       integer icntl(*)
       double precision   cntl(*)
*
       icntl(1) = 6
       icntl(2) = 6
       icntl(3) = 0
       icntl(4) = 4
       icntl(5) = 0
       icntl(6) = 0
       icntl(7) = -1
       icntl(8) = 1
* 
       cntl(1) = 1.0 d -5
       cntl(2) = 0.0 d 0
       cntl(3) = 0.0 d 0
       cntl(4) = 0.0 d 0
       cntl(5) = 0.0 d 0
*
       return
       end
