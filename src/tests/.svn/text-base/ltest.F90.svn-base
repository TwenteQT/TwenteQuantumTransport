#include "math_def.h"

Program ltest
!!$       Use sparselib
!!$       Use potpars
!!$       Use rotations
!!$       Implicit None
!!$       Integer, Parameter :: l = 2
!!$ !!$       Complex (Kind=DEF_DBL_PREC) :: Lp (2*l-1, 2*l-1), Lm (2*l-1, 2*l-1), Lz (2*l-1, 2*l-1)
!!$ !!$       Real (Kind=DEF_DBL_PREC) :: Rm(3,3)
!!$       Complex (Kind=8) :: qq (2, 2)
!!$       Complex (Kind=DEF_DBL_PREC), Parameter :: sz (2, 2) = reshape ( (/ 0.5d0, 0.0d0, 0.0d0,-0.5d0 /), (/ 2, &
!!$      & 2 /))
!!$       Complex (Kind=DEF_DBL_PREC), Parameter :: Sm (2, 2) = reshape ( (/ 0.0d0, 1.0d0, 0.0d0, 0.0d0 /), (/ 2, &
!!$      & 2 /))
!!$       Complex (Kind=DEF_DBL_PREC), Parameter :: Sp (2, 2) = reshape ( (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0 /), (/ 2, &
!!$      & 2 /))
!!$       Complex (Kind=DEF_DBL_PREC), Parameter :: ci = DEF_cmplx_Ione
!!$ 
!!$       Real (Kind=8) :: theta (2) = (/ 0.75d0, .6d0 /)
!!$       Complex (Kind=DEF_DBL_PREC) :: rm (2, 2), sx (2, 2), Sy (2, 2), sn (2, 2), sx1 (2, 2), sy1 (2, 2)
!!$       Complex (Kind=DEF_DBL_PREC) :: r3 (3, 3), Lp (5, 5), Lm (5, 5), Lz (5, 5), Lx (5, 5), Ly (5, 5), Lx1 &
!!$      & (5, 5)
!!$       Complex (Kind=DEF_DBL_PREC) :: LR (5, 5), LRi (5, 5), Lx2 (5, 5), Ly1 (5, 5), Ly2 (5, 5)
!!$       Type (t_Lop_mats) :: Lop (3)
!!$       integer :: m
!!$ 
!!$ 101   Format (4(f10.5))
!!$ 
!!$       Call make_Lmats (3, Lop)
!!$ 
!!$       sx = .5d0 * (Sp+Sm)
!!$       Sy = - ci * .5d0 * (Sp-Sm)
!!$ 
!!$       rm = spinrm (theta)
!!$       r3 = rm3d (theta)
!!$ 
!!$       sn = rotspin (Sp, rm)
!!$ 
!!$       Write (*, 101) sn
!!$       Write (*,*)
!!$       sx1 = r3 (1, 1) * sx + r3 (1, 2) * Sy + r3 (1, 3) * sz
!!$       sy1 = r3 (2, 1) * sx + r3 (2, 2) * Sy + r3 (2, 3) * sz
!!$ 
!!$       sn = sx1 + ci * sy1
!!$       Write (*, 101) sn
!!$       Write (*,*)
!!$ 
!!$ !!$ 
!!$ !!$       Lp = Lop(3)%Lp
!!$ !!$       Lm = Lop(3)%Lm
!!$ !!$       Lz = Lop(3)%Lz
!!$ !!$ 
!!$ !!$ 
!!$ !!$       Lx = .5d0 * (Lp+Lm)
!!$ !!$       Ly = - ci * .5d0 * (Lp-Lm)
!!$ !!$ !!$       LR=matmul(expm(ci*Ly,theta(1)),expm(ci*Lz,theta(2)))
!!$ !!$       LR = rmLop (Ly, Lz, theta)
!!$ !!$       LRi = transpose (conjg(LR))
!!$ !!$ 
!!$ !!$ !!$       Lx1=r3(1,1)*Lx+r3(1,2)*Ly+r3(1,3)*Lz
!!$ !!$       Lx1 = rot3D (Lx, Ly, Lz, r3, 1)
!!$ !!$       Ly1 = r3 (2, 1) * Lx + r3 (2, 2) * Ly + r3 (2, 3) * Lz
!!$ !!$ 
!!$ !!$ !!$       Lx2=matmul(LRi,matmul(Lx,LR))
!!$ !!$       Lx2 = RotL (Lx, LR)
!!$ !!$       Ly2 = matmul (LRi, matmul(Ly, LR))
!!$ !!$ 
!!$ !!$       Write (*,*) maxval (Abs(Lx1-Lx2))
!!$ !!$       Write (*,*) maxval (Abs(Ly1-Ly2))
!!$ !!$ 102   Format (5(f10.5))
!!$ !!$ 
!!$ !!$       Write (*,102) Abs (matmul(Lm, Lp))
!!$ !!$       write(*,*)
!!$ !!$       Write (*,102) Abs (matmul(Lp, Lm))
!!$ 103   Format (5("("f10.5" "f10.5") "))
!!$       Write (*, 103) (Lp(m,:),m=1,5)
!!$       Write (*,*)
!!$       Write (*, 103) (Lz(m,:),m=1,5)
!!$ !!$       Write (*,*)
!!$ !!$       Write (*, 103) Lm
!!$ !!$       Write (*,*)
!!$ !!$       Write (*, 103) Lm-transpose(conjg(Lp))
!!$ 
!!$ !!$
!!$ !!$       rm = Srot3D (theta1)
!!$ !!$       Write (*, 101) rm
!!$ !!$       Write (*,*)
!!$ !!$
!!$ !!$       rm = Srot3D (theta3)
!!$ !!$       Write (*, 101) rm
!!$ !!$
!!$ !!$       Stop
!!$ !!$       rm = Srot3D (theta1)
!!$ !!$       Write (*, 101) rm
!!$ !!$       Write (*, 101) rotSpin (qq, rm)
!!$ !!$       Write (*,*)
!!$ !!$       rm = Srot3D (theta2)
!!$ !!$       Write (*, 101) rm
!!$ !!$       Write (*, 101) rotSpin (qq, rm)
!!$ !!$       Write (*,*)
!!$ !!$      rm=Srot3D (theta3)
!!$ !!$      write(*,101) rotSpin(qq,rm)
!!$ !!$      write(*,*)
!!$ !!$      rm=Srot3D (theta4)
!!$ !!$      write(*,101) rotSpin(qq,rm)
!!$ !!$      write(*,*)
!!$ !!$      rm=Srot3D (theta5)
!!$ !!$      write(*,101) rotSpin(qq,rm)
!!$ 
!!$ 
!!$ !!$       Call CalcLInts (Lp, Lm, Lz, l)
!!$ !!$
!!$ !!$       Call dump_matrix (Lz, 'Lz')
!!$ !!$       Call dump_matrix (Lp, 'Lp')
!!$ !!$       Call dump_matrix (Lm, 'Lm')
!!$ !!$       Call rotmat3d (0.0d0,0.0d0,DEF_M_PI/6.0d0,Rmz)
!!$ !!$       Call rotmat3d (0.0d0,DEF_M_PI/13.0d0,0.0d0,Rmy)
!!$ !!$       Call rotmat3d (DEF_M_PI/5.0d0,0.0d0,0.0d0,Rmx)
!!$ !!$       Call rotmat3d (0.0d0,DEF_M_PI/2.0d0,0.0d0,Rm)
!!$ !!$       Rm=matmul(Rmz,matmul(Rmy,Rmx))
!!$ !!$      Call dump_matrix (Rm, 'rm')        
!!$ !!$      call make_angints (ai, 5, theta )
!!$ 
!!$ !     call rotateLints(Lp,Lm,Lz,Rm,l)
!!$ !!$      Call dump_matrix (ai%ud(l), 'ud')
!!$ !!$           Call dump_matrix (ai%du(l), 'du')
!!$ !!$      Call dump_matrix (Lp, 'nqqq1')
!!$ !!$      Call dump_matrix (Lm, 'nqqq2')
!!$ !!$      Call dump_matrix (q, 'lpe')
!!$ !!$      call prep_so (sopar, geom, soopt)
!!$ 
!!$       Write (*,*) 'done!'
End Program ltest
