!!$ $Id: $
#include "math_def.h"
!#define DEF_DBL_PREC 8
Module rotations
      Implicit None
Contains
      Function read_rot_mask (num) Result (rm)
         Use logging
         Implicit None
         Integer :: num
         Real (Kind=DEF_DBL_PREC), Pointer :: rm (:, :)
!!$ Local vars
         Integer :: fl = 32898, ios
         Integer :: i

!!$      Call do_log (1, 'Reading rotation mask...')
         Allocate (rm(2, num))
         rm = 0.0d0
         Open (Unit=fl, File='rot_mask', Action='read', IoStat=ios)
         If (ios /= 0) Then
            Write (*,*) 'Error opening file "rot_mask".'
            Stop
         End If
         Read (fl,*)
         Read (fl, '(2g17.10)') (rm(1:2, i), i=1, num)
         Close (Unit=fl)
         rm (:, :) = rm (:, :) * DEF_M_PI
!!$      Call do_log (1, 'Done!')
      End Function read_rot_mask

      Function spinrm (thetai) Result (rm)
!!$ Generate rotation matrix in spin space
         Implicit None
         Real (Kind=DEF_DBL_PREC), Intent (In) :: thetai (2)
         Complex (Kind=DEF_DBL_PREC) :: rm (2, 2)
!!$ Local vars
         Complex (Kind=DEF_QUAD_PREC), Parameter :: ci = DEF_cmplx_IoneL
         Complex (Kind=DEF_QUAD_PREC), Parameter :: sz (2, 2) = reshape ( (/ 0.5q0, 0.0q0, 0.0q0,-0.5q0 &
        & /), (/ 2, 2 /))
         Complex (Kind=DEF_QUAD_PREC), Parameter :: Sy (2, 2) = ci * reshape ( (/ 0.0q0, .5q0,-.5q0, 0.0q0 &
        & /), (/ 2, 2 /))
         Real (Kind=DEF_QUAD_PREC) :: theta (2)
         Complex (Kind=DEF_QUAD_PREC) :: rmq (2, 2)
!!$          Integer :: i, j


         theta = thetai
         rmq (1, 1) = Exp (-ci*theta(2)/2.d0) * Cos (theta(1)/2.d0)
         rmq (1, 2) = - Exp (ci*theta(2)/2.d0) * Sin (theta(1)/2.d0)
         rmq (2, 1) = - conjg (rmq(1, 2))
         rmq (2, 2) = conjg (rmq(1, 1))

         rm = rmq
         rm = matmul (expm(dcmplx(-ci*Sy), thetai(1)), expm(dcmplx(-ci*sz), thetai(2)))
      End Function spinrm

!!$       Function rmLop (Ly, Lz, theta) Result (rm)
!!$          Implicit None
!!$          Real (Kind=DEF_DBL_PREC), Intent (In) :: theta (2)
!!$          Complex (Kind=DEF_DBL_PREC), Intent (In) :: Ly (:, :), Lz (:, :)
!!$          Complex (Kind=DEF_DBL_PREC) :: rm (size(Ly, 1), size(Ly, 2))
!!$          Complex (Kind=DEF_DBL_PREC), Parameter :: ci = DEF_cmplx_IoneL
!!$          Integer :: i, j
!!$
!!$         rm = matmul (expm(-ci*Ly, theta(1)),expm(-ci*Lz, theta(2)))
!!$       End Function rmLop
!!$
!!$
!!$       Function rm3D (thetai) Result (rm)
!!$          Implicit None
!!$          Real (Kind=DEF_DBL_PREC), Intent (In) :: thetai (2)
!!$          Complex (Kind=DEF_DBL_PREC) :: rm (3, 3)
!!$          Real (Kind=DEF_QUAD_PREC) :: theta (2)
!!$          Complex (Kind=DEF_QUAD_PREC) :: rmq (3, 3)
!!$          theta = thetai
!!$
!!$          rmq (1, :) = (/Cos(theta(1))*Cos(theta(2)),-(Cos(theta(1))*Sin(theta(2))),Sin(theta(1))/)
!!$          rmq (2, :) = (/Sin(theta(2)),Cos(theta(2)),0/)
!!$          rmq (3, :) = (/ -(Cos(theta(2))*Sin(theta(1))),Sin(theta(1))*Sin(theta(2)),Cos(theta(1))/)
!!$          rm = rmq
!!$       End Function rm3D

!!$       Function rotSpin (s, R) Result (SR)
!!$          Implicit None
!!$          Complex (Kind=DEF_DBL_PREC), Intent (In) :: s (2, 2), R (2, 2)
!!$          Complex (Kind=DEF_DBL_PREC) :: SR (2, 2)
!!$ !!$ Locals
!!$          Complex (Kind=DEF_QUAD_PREC) :: t1, t3, t5, t9, t15, t17, t19, t23
!!$          t1 = conjg (R(1, 1))
!!$          t3 = conjg (R(2, 1))
!!$          t5 = t1 * s (1, 1) + t3 * s (2, 1)
!!$          t9 = t1 * s (1, 2) + t3 * s (2, 2)
!!$          t15 = conjg (R(1, 2))
!!$          t17 = conjg (R(2, 2))
!!$          t19 = t15 * s (1, 1) + t17 * s (2, 1)
!!$          t23 = t15 * s (1, 2) + t17 * s (2, 2)
!!$          SR (1, 1) = t5 * R (1, 1) + t9 * R (2, 1)
!!$          SR (1, 2) = t5 * R (1, 2) + t9 * R (2, 2)
!!$          SR (2, 1) = t19 * R (1, 1) + t23 * R (2, 1)
!!$          SR (2, 2) = t19 * R (1, 2) + t23 * R (2, 2)
!!$       End Function rotSpin

!!$       Function rotL (L, R) Result (LR)
!!$          Implicit None
!!$          Complex (Kind=DEF_DBL_PREC), Intent (In) :: L (:, :), R (:, :)
!!$          Complex (Kind=DEF_DBL_PREC) :: LR (size(L, 1), size(L, 2))
!!$
!!$          LR = matmul (transpose(conjg(R)), matmul(L, R))
!!$       End Function rotL
!!$
!!$       Function Rot3D (mx, my, mz, rm, ind) Result (nm)
!!$          Implicit None
!!$          Complex (Kind=DEF_DBL_PREC), Intent (In) :: mx (:, :), my (:, :), mz (:, :), rm (3, 3)
!!$          Complex (Kind=DEF_DBL_PREC) :: nm (size(mx, 1), size(mx, 2))
!!$          Integer, Intent (In) :: ind
!!$          nm = rm (ind, 1) * mx + rm (ind, 2) * my + rm (ind, 3) * mz
!!$       End Function Rot3D

      Function getderangle (base, ax, step, dir, ind, f1) Result (ang)
         Implicit None
         Real (Kind=DEF_DBL_PREC) :: base (2), ang (2), step, ax (2)
         Real (Kind=DEF_DBL_PREC), Optional :: f1
         Integer :: dir, ind
         Real (Kind=DEF_DBL_PREC) :: ing (3)


         ing (1) = step / 2.0d0
         ing (2) = DEF_M_PI * (0.5d0*real(dir-1, kind=8)-real(ind-1, kind=8))
         If (present(f1)) ing (2) = ing (2) + f1
         ing (3) = 0

         ang(1:2) = rotspheric (base, ing, ax)

      End Function



      Recursive Function rotspheric (veci, rotang, ax) Result (v)
         Real (Kind=DEF_DBL_PREC), Intent (In) :: veci (2), rotang (3)
         Real (Kind=DEF_DBL_PREC), Optional, Intent (In) :: ax (2)
         Real (Kind=DEF_DBL_PREC) :: v (2)
!!$ local
         Real (Kind=DEF_DBL_PREC) :: vec (2), tv, rf (3), v1 (2)
         Logical :: drot

         drot = .False.
         vec(:) = veci(:)

         If (present(ax)) Then
            If (maxval(Abs(ax)) > 1.0d-12) Then
               drot = .True.
               rf = (/ ax (1), ax (2) + DEF_M_PI, 0.0d0 /)
               vec(1:2) = rotspheric (veci, rf)
            End If
         End If
         
         v1 (1) = Acos &
        & (Cos(rotang(1))*Cos(vec(1))-Cos(rotang(2)-rotang(3)-vec(2))*Sin(rotang(1))*Sin(vec(1)))
         tv = rotang (2) + Atan2 (-(Sin(vec(1))*Sin(rotang(2)-rotang(3)-vec(2))), &
        & Cos(vec(1))*Sin(rotang(1))+Cos(rotang(1))*Cos(rotang(2)-rotang(3)-vec(2))*Sin(vec(1)))
         v1 (2) = Atan2 (Sin(tv), Cos(tv))

         If (drot) Then
            rf (2) = ax (2)
            v(1:2) = rotspheric (v1, rf)
         Else
            v(1:2) = v1(1:2)
         End If

      End Function rotspheric

!!$       Function rotspheric (basei, anglei) Result (nang)
!!$          Implicit None
!!$          Real (Kind=DEF_DBL_PREC) :: basei (2), anglei (2), nango (2)
!!$          Real (Kind=DEF_QUAD_PREC), Parameter :: pi = DEF_M_PI
!!$          Real (Kind=DEF_QUAD_PREC) :: base (2), angle (2), nang (2)
!!$          base = - basei
!!$          angle = anglei
!!$
!!$          nang (1) = Acos (Cos(angle(1))*Cos(base(1))+Cos(angle(2))*Sin(angle(1))*Sin(base(1)))
!!$
!!$          nang (2) = Atan2 (Cos(base(2))*Sin(angle(1))*Sin(angle(2))+(-&
!!$         & (Cos(angle(2))*Cos(base(1))*Sin(angle(1)))+Cos(angle(1))*Sin(base(1)))*Sin(base(2)), &
!!$         & Cos(angle(2))*Cos(base(1))*Cos(base(2))*Sin(angle(1))-Cos(angle(1))*Cos(base(2))*Sin(base(1))+&
!!$         & Sin(angle(1))*Sin(angle(2))*Sin(base(2)))
!!$          nango = nang
!!$       End Function

      Function zmatmul3 (R1, R2, R3) Result (rm)
         Implicit None
         Complex (Kind=DEF_DBL_PREC) :: R1 (:, :), R2 (:, :), R3 (:, :)
         Complex (Kind=DEF_DBL_PREC) :: rm (size(R1, 1), size(R3, 2))
         rm = matmul (R1, matmul(R2, R3))
      End Function zmatmul3

!!$----------------------------------------------------------------------|
      Function expm (H, t_in, iflag_out) Result (rexp)
         Implicit None
         Real (Kind=DEF_DBL_PREC), Intent (In), Optional :: t_in
         Integer, Intent (Out), Optional :: iflag_out
         Complex (Kind=DEF_DBL_PREC) :: H (:, :), rexp (size(H, 1), size(H, 2))


         Integer, Parameter :: ideg = 6
         Integer :: lwsp, ipiv (size(H, 2)), ns
         Complex (Kind=DEF_DBL_PREC) :: wsp (4*size(H, 2)**2+ideg+1)
         Integer :: m, iflag
         Real (Kind=DEF_DBL_PREC) :: T

!!$-----Purpose----------------------------------------------------------|
!!$
!!$     Computes exp(t*H), the matrix exponential of a general complex
!!$     matrix in full, using the irreducible rational Pade approximation
!!$     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
!!$     combined with scaling-and-squaring.
!!$
!!$-----Arguments--------------------------------------------------------|
!!$
!!$     ideg      : (input) the degre of the diagonal Pade to be used.
!!$                 a value of 6 is generally satisfactory.
!!$
!!$     H(ldh,m)  : (input) argument matrix.
!!$
!!$     t         : (input) time-scale (can be < 0).
!!$
!!$     ns        : (output) number of scaling-squaring used.
!!$
!!$     iflag     : (output) exit flag.
!!$                       0 - no problem
!!$                      <0 - problem
!!$----------------------------------------------------------------------|
         Integer :: i, j, k, icoef, mm, ih2, iodd, iused, ifree, iq, ip, iput, iget
         Real (Kind=DEF_DBL_PREC) :: hnorm
         Complex (Kind=DEF_DBL_PREC) :: cp, cq, scale, scale2, ZERO, ONE

         Parameter (ZERO=(0.0d0, 0.0d0), ONE=(1.0d0, 0.0d0))

         m = size (H, 2)
         lwsp = 4 * m * m + ideg + 1
         T = 1.0d0
         If (present(t_in)) T = t_in
!!$---  check restrictions on input parameters ...
         mm = m * m
         iflag = 0
         If (size(H, 1) /= m) Stop 'matrix should be square'
         If (lwsp .Lt. 4*mm+ideg+1) iflag = - 2
         If (iflag .Ne. 0) Stop 'bad sizes (in input of ZGPADM)'

!!$---  initialise pointers ...
         icoef = 1
         ih2 = icoef + (ideg+1)
         ip = ih2 + mm
         iq = ip + mm
         ifree = iq + mm

!!$---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
!!$     and set scale = t/2^ns ...

         wsp (:) = ZERO
         Do j = 1, m
            Do i = 1, m
               wsp (i) = wsp (i) + Abs (H(i, j))
            End Do
         End Do
         hnorm = 0.0d0
         Do i = 1, m
            hnorm = Max (hnorm, dble(wsp(i)))
         End Do
         hnorm = dAbs (T*hnorm)
         If ((dAbs(hnorm) .Le. 1.0d-8)) Then
            iput = 1
            wsp (:) = ZERO
            Do k = 1, m
               wsp ((k-1)*m+k) = 1.0d0
            End Do
            Go To 200
         End If

         ns = Max (0, Int(dLog(hnorm)/dLog(2.0d0))+2)
         scale = dcmplx (T/dble(2**ns), 0.0d0)
         scale2 = scale * scale

!!$---  compute Pade coefficients ...

         i = ideg + 1
         j = 2 * ideg + 1
         wsp (icoef) = ONE
         Do k = 1, ideg
            wsp (icoef+k) = (wsp(icoef+k-1)*dble(i-k)) / dble (k*(j-k))
         End Do

!!$---  H2 = scale2*H*H ...

         Call ZGEMM ('n', 'n', m, m, m, scale2, H, m, H, m, ZERO, wsp(ih2), m)

!!$---  initialise p (numerator) and q (denominator) ...

         cp = wsp (icoef+ideg-1)
         cq = wsp (icoef+ideg)
         Do j = 1, m
            Do i = 1, m
               wsp (ip+(j-1)*m+i-1) = ZERO
               wsp (iq+(j-1)*m+i-1) = ZERO
            End Do
            wsp (ip+(j-1)*(m+1)) = cp
            wsp (iq+(j-1)*(m+1)) = cq
         End Do

!!$---  Apply Horner rule ...

         iodd = 1
         k = ideg - 1
100      Continue
         iused = iodd * iq + (1-iodd) * ip
         Call ZGEMM ('n', 'n', m, m, m, ONE, wsp(iused), m, wsp(ih2), m, ZERO, wsp(ifree), m)
         Do j = 1, m
            wsp (ifree+(j-1)*(m+1)) = wsp (ifree+(j-1)*(m+1)) + wsp (icoef+k-1)
         End Do
         ip = (1-iodd) * ifree + iodd * ip
         iq = iodd * ifree + (1-iodd) * iq
         ifree = iused
         iodd = 1 - iodd
         k = k - 1
         If (k .Gt. 0) Go To 100

!!$---  Obtain (+/-)(I + 2*(p\q)) ...

         If (iodd .Ne. 0) Then
            Call ZGEMM ('n', 'n', m, m, m, scale, wsp(iq), m, H, m, ZERO, wsp(ifree), m)
            iq = ifree
         Else
            Call ZGEMM ('n', 'n', m, m, m, scale, wsp(ip), m, H, m, ZERO, wsp(ifree), m)
            ip = ifree
         End If
         Call ZAXPY (mm,-ONE, wsp(ip), 1, wsp(iq), 1)
         Call ZGESV (m, m, wsp(iq), m, ipiv, wsp(ip), m, iflag)
         If (iflag .Ne. 0) Stop 'Problem in ZGESV (within ZGPADM)'
         Call ZDSCAL (mm, 2.0d0, wsp(ip), 1)
         Do j = 1, m
            wsp (ip+(j-1)*(m+1)) = wsp (ip+(j-1)*(m+1)) + ONE
         End Do
         iput = ip
         If (ns .Eq. 0 .And. iodd .Ne. 0) Then
            Call ZDSCAL (mm,-1.0d0, wsp(ip), 1)
            Go To 200
         End If

!!$--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...

         iodd = 1
         Do k = 1, ns
            iget = iodd * ip + (1-iodd) * iq
            iput = (1-iodd) * ip + iodd * iq
            Call ZGEMM ('n', 'n', m, m, m, ONE, wsp(iget), m, wsp(iget), m, ZERO, wsp(iput), m)
            iodd = 1 - iodd
         End Do
200      Continue

         Do k = 1, m
            rexp (1:m, k) = wsp (iput+(k-1)*m:iput+k*m-1)
         End Do

         If (present(iflag_out)) iflag_out = iflag
      End Function expm

!!$       Function expmq (H, t_in, ideg, iflag_out) Result (rexp)
!!$          Implicit None
!!$          Real (Kind=DEF_QUAD_PREC), Intent (In), Optional :: t_in
!!$          Integer, Intent (Out), Optional :: iflag_out!,ideg_in
!!$          Complex (Kind=DEF_QUAD_PREC) :: H (:, :), rexp (size(H, 1), size(H, 2))
!!$
!!$
!!$          Integer :: ideg
!!$          Integer :: lwsp, ipiv (size(H, 2)), ns
!!$          Complex (Kind=DEF_QUAD_PREC) :: wsp (4*size(H, 2)**2+ideg+1)
!!$          Integer :: m, iflag
!!$          Real (Kind=DEF_QUAD_PREC) :: T
!!$
!!$ !!$-----Purpose----------------------------------------------------------|
!!$ !!$
!!$ !!$     Computes exp(t*H), the matrix exponential of a general complex
!!$ !!$     matrix in full, using the irreducible rational Pade approximation
!!$ !!$     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
!!$ !!$     combined with scaling-and-squaring.
!!$ !!$
!!$ !!$-----Arguments--------------------------------------------------------|
!!$ !!$
!!$ !!$     ideg      : (input) the degre of the diagonal Pade to be used.
!!$ !!$                 a value of 6 is generally satisfactory.
!!$ !!$
!!$ !!$     H(ldh,m)  : (input) argument matrix.
!!$ !!$
!!$ !!$     t         : (input) time-scale (can be < 0).
!!$ !!$
!!$ !!$     ns        : (output) number of scaling-squaring used.
!!$ !!$
!!$ !!$     iflag     : (output) exit flag.
!!$ !!$                       0 - no problem
!!$ !!$                      <0 - problem
!!$ !!$----------------------------------------------------------------------|
!!$          Integer :: i, j, k, icoef, mm, ih2, iodd, iused, ifree, iq, ip, iput, iget
!!$          Real (Kind=DEF_QUAD_PREC) :: hnorm
!!$          Complex (Kind=DEF_QUAD_PREC) :: cp, cq, scale, scale2, ZERO, ONE
!!$
!!$          Parameter (ZERO=(0.0q0, 0.0q0), ONE=(1.0q0, 0.0q0))
!!$
!!$          m = size (H, 2)
!!$          lwsp = 4 * m * m + ideg + 1
!!$          T = 1.0q0
!!$          If (present(t_in)) T = t_in
!!$ !!$          ideg=6
!!$ !!$          If (present(ideg_in)) ideg = ideg_in
!!$ !!$---  check restrictions on input parameters ...
!!$          mm = m * m
!!$          iflag = 0
!!$          If (size(H, 1) /= m) Stop 'matrix should be square'
!!$          If (lwsp .Lt. 4*mm+ideg+1) iflag = - 2
!!$          If (iflag .Ne. 0) Stop 'bad sizes (in input of ZGPADM)'
!!$
!!$ !!$---  initialise pointers ...
!!$          icoef = 1
!!$          ih2 = icoef + (ideg+1)
!!$          ip = ih2 + mm
!!$          iq = ip + mm
!!$          ifree = iq + mm
!!$
!!$ !!$---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
!!$ !!$     and set scale = t/2^ns ...
!!$
!!$          wsp (:) = ZERO
!!$          Do j = 1, m
!!$             Do i = 1, m
!!$                wsp (i) = wsp (i) + Abs (H(i, j))
!!$             End Do
!!$          End Do
!!$          hnorm = 0.0q0
!!$          Do i = 1, m
!!$             hnorm = Max (hnorm, real(wsp(i), kind=16))
!!$          End Do
!!$          hnorm = Abs (T*hnorm)
!!$          If ((Abs(hnorm) .Le. 1.0q-16)) Then
!!$             iput = 1
!!$             wsp (:) = ZERO
!!$             Do k = 1, m
!!$                wsp ((k-1)*m+k) = 1.0q0
!!$             End Do
!!$             Go To 200
!!$          End If
!!$          ns = Max (0, Int(Log(hnorm)/Log(2.0q0))+2)
!!$          scale = CMPLX (T/real(2**ns, kind=16), 0.0q0)
!!$          scale2 = scale * scale
!!$
!!$ !!$---  compute Pade coefficients ...
!!$
!!$          i = ideg + 1
!!$          j = 2 * ideg + 1
!!$          wsp (icoef) = ONE
!!$          Do k = 1, ideg
!!$             wsp (icoef+k) = (wsp(icoef+k-1)*real(i-k, kind=16)) / real (k*(j-k), kind=16)
!!$          End Do
!!$
!!$
!!$ !!$---  H2 = scale2*H*H ...
!!$
!!$          Call ZGEMMQ ('n', 'n', m, m, m, scale2, H, m, H, m, ZERO, wsp(ih2), m)
!!$
!!$ !!$---  initialise p (numerator) and q (denominator) ...
!!$
!!$          cp = wsp (icoef+ideg-1)
!!$          cq = wsp (icoef+ideg)
!!$          Do j = 1, m
!!$             Do i = 1, m
!!$                wsp (ip+(j-1)*m+i-1) = ZERO
!!$                wsp (iq+(j-1)*m+i-1) = ZERO
!!$             End Do
!!$             wsp (ip+(j-1)*(m+1)) = cp
!!$             wsp (iq+(j-1)*(m+1)) = cq
!!$          End Do
!!$
!!$
!!$ !!$---  Apply Horner rule ...
!!$
!!$          iodd = 1
!!$          k = ideg - 1
!!$ 100      Continue
!!$          iused = iodd * iq + (1-iodd) * ip
!!$          Call ZGEMMQ ('n', 'n', m, m, m, ONE, wsp(iused), m, wsp(ih2), m, ZERO, wsp(ifree), m)
!!$          Do j = 1, m
!!$             wsp (ifree+(j-1)*(m+1)) = wsp (ifree+(j-1)*(m+1)) + wsp (icoef+k-1)
!!$          End Do
!!$          ip = (1-iodd) * ifree + iodd * ip
!!$          iq = iodd * ifree + (1-iodd) * iq
!!$          ifree = iused
!!$          iodd = 1 - iodd
!!$          k = k - 1
!!$          If (k .Gt. 0) Go To 100
!!$
!!$ !!$---  Obtain (+/-)(I + 2*(p\q)) ...
!!$
!!$          If (iodd .Ne. 0) Then
!!$             Call ZGEMMQ ('n', 'n', m, m, m, scale, wsp(iq), m, H, m, ZERO, wsp(ifree), m)
!!$             iq = ifree
!!$          Else
!!$             Call ZGEMMQ ('n', 'n', m, m, m, scale, wsp(ip), m, H, m, ZERO, wsp(ifree), m)
!!$             ip = ifree
!!$          End If
!!$          Call ZAXPYQ (mm,-ONE, wsp(ip), 1, wsp(iq), 1)
!!$
!!$          Call ZGESVQ (m, m, wsp(iq), m, ipiv, wsp(ip), m, iflag)
!!$ !!$          Call QGESV (m, m, wsp(iq), m, ipiv, wsp(ip), m, iflag)
!!$
!!$          If (iflag .Ne. 0) Stop 'Problem in ZGESV (within ZGPADM)'
!!$          Call ZDSCALQ (mm, 2.0q0, wsp(ip), 1)
!!$
!!$          Do j = 1, m
!!$             wsp (ip+(j-1)*(m+1)) = wsp (ip+(j-1)*(m+1)) + ONE
!!$          End Do
!!$          iput = ip
!!$
!!$          If (ns .Eq. 0 .And. iodd .Ne. 0) Then
!!$             Call ZDSCALQ (mm,-1.0q0, wsp(ip), 1)
!!$             Go To 200
!!$          End If
!!$
!!$ !!$--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
!!$
!!$          iodd = 1
!!$          Do k = 1, ns
!!$             iget = iodd * ip + (1-iodd) * iq
!!$             iput = (1-iodd) * ip + iodd * iq
!!$             Call ZGEMMQ ('n', 'n', m, m, m, ONE, wsp(iget), m, wsp(iget), m, ZERO, wsp(iput), m)
!!$             iodd = 1 - iodd
!!$          End Do
!!$ 200      Continue
!!$
!!$          Do k = 1, m
!!$             rexp (1:m, k) = wsp (iput+(k-1)*m:iput+k*m-1)
!!$          End Do
!!$
!!$          If (present(iflag_out)) iflag_out = iflag
!!$       End Function expmq


End Module rotations
!!$
!!$       Subroutine ZGEMMQ (TRANSA, TRANSB, m, N, k, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!!$          Complex (Kind=DEF_QUAD_PREC) :: ALPHA, BETA
!!$          Integer k, LDA, LDB, LDC, m, N
!!$          Character TRANSA, TRANSB
!!$          Complex (Kind=DEF_QUAD_PREC) :: A (LDA,*), B (LDB,*), C (LDC,*)
!!$          Logical LSAME
!!$          External LSAME
!!$          External XERBLA
!!$ !!$          Intrinsic DCONJG, Max
!!$          Complex (Kind=DEF_QUAD_PREC) :: TEMP
!!$          Integer i, INFO, j, L, NCOLA, NROWA, NROWB
!!$          Logical CONJA, CONJB, NOTA, NOTB
!!$          Complex (Kind=DEF_QUAD_PREC) :: ONE
!!$          Parameter (ONE=(1.0q+0, 0.0q+0))
!!$          Complex (Kind=DEF_QUAD_PREC) :: ZERO
!!$          Parameter (ZERO=(0.0q+0, 0.0q+0))
!!$          NOTA = LSAME (TRANSA, 'N')
!!$          NOTB = LSAME (TRANSB, 'N')
!!$          CONJA = LSAME (TRANSA, 'C')
!!$          CONJB = LSAME (TRANSB, 'C')
!!$          If (NOTA) Then
!!$             NROWA = m
!!$             NCOLA = k
!!$          Else
!!$             NROWA = k
!!$             NCOLA = m
!!$          End If
!!$          If (NOTB) Then
!!$             NROWB = k
!!$          Else
!!$             NROWB = N
!!$          End If
!!$          INFO = 0
!!$          If (( .Not. NOTA) .And. ( .Not. CONJA) .And. (.NOT. LSAME(TRANSA, 'T'))) THEN
!!$          INFO = 1
!!$       Else If (( .Not. NOTB) .And. ( .Not. CONJB) .And.(.NOT. LSAME(TRANSB, 'T'))) THEN
!!$          INFO = 2
!!$       Else If (m .Lt. 0) Then
!!$          INFO = 3
!!$       Else If (N .Lt. 0) Then
!!$          INFO = 4
!!$       Else If (k .Lt. 0) Then
!!$          INFO = 5
!!$       Else If (LDA .Lt. Max(1, NROWA)) Then
!!$          INFO = 8
!!$       Else If (LDB .Lt. Max(1, NROWB)) Then
!!$          INFO = 10
!!$       Else If (LDC .Lt. Max(1, m)) Then
!!$          INFO = 13
!!$       End If
!!$       If (INFO .Ne. 0) Then
!!$          Call XERBLA ('ZGEMM ', INFO)
!!$          Return
!!$       End If
!!$       If ((m .Eq. 0) .Or. (N .Eq. 0) .Or. (((ALPHA .Eq. ZERO) .Or. (k .Eq. 0)) .And. (BETA .Eq. ONE))) RETURN
!!$       If (ALPHA .Eq. ZERO) Then
!!$          If (BETA .Eq. ZERO) Then
!!$             Do 20 j = 1, N
!!$                Do 10 i = 1, m
!!$                   C (i, j) = ZERO
!!$ 10             Continue
!!$ 20          Continue
!!$          Else
!!$             Do 40 j = 1, N
!!$                Do 30 i = 1, m
!!$                   C (i, j) = BETA * C (i, j)
!!$ 30             Continue
!!$ 40          Continue
!!$          End If
!!$          Return
!!$       End If
!!$       If (NOTB) Then
!!$          If (NOTA) Then
!!$             Do 90 j = 1, N
!!$                If (BETA .Eq. ZERO) Then
!!$                   Do 50 i = 1, m
!!$                      C (i, j) = ZERO
!!$ 50                Continue
!!$                Else If (BETA .Ne. ONE) Then
!!$                   Do 60 i = 1, m
!!$                      C (i, j) = BETA * C (i, j)
!!$ 60                Continue
!!$                End If
!!$                Do 80 L = 1, k
!!$                   If (B(L, j) .Ne. ZERO) Then
!!$                      TEMP = ALPHA * B (L, j)
!!$                      Do 70 i = 1, m
!!$                         C (i, j) = C (i, j) + TEMP * A (i, L)
!!$ 70                   Continue
!!$                   End If
!!$ 80             Continue
!!$ 90          Continue
!!$          Else If (CONJA) Then
!!$             Do 120 j = 1, N
!!$                Do 110 i = 1, m
!!$                   TEMP = ZERO
!!$                   Do 100 L = 1, k
!!$                      TEMP = TEMP + QCONJG (A(L, i)) * B (L, j)
!!$ 100               Continue
!!$                   If (BETA .Eq. ZERO) Then
!!$                      C (i, j) = ALPHA * TEMP
!!$                   Else
!!$                      C (i, j) = ALPHA * TEMP + BETA * C (i, j)
!!$                   End If
!!$ 110            Continue
!!$ 120         Continue
!!$          Else
!!$             Do 150 j = 1, N
!!$                Do 140 i = 1, m
!!$                   TEMP = ZERO
!!$                   Do 130 L = 1, k
!!$                      TEMP = TEMP + A (L, i) * B (L, j)
!!$ 130               Continue
!!$                   If (BETA .Eq. ZERO) Then
!!$                      C (i, j) = ALPHA * TEMP
!!$                   Else
!!$                      C (i, j) = ALPHA * TEMP + BETA * C (i, j)
!!$                   End If
!!$ 140            Continue
!!$ 150         Continue
!!$          End If
!!$       Else If (NOTA) Then
!!$          If (CONJB) Then
!!$             Do 200 j = 1, N
!!$                If (BETA .Eq. ZERO) Then
!!$                   Do 160 i = 1, m
!!$                      C (i, j) = ZERO
!!$ 160               Continue
!!$                Else If (BETA .Ne. ONE) Then
!!$                   Do 170 i = 1, m
!!$                      C (i, j) = BETA * C (i, j)
!!$ 170               Continue
!!$                End If
!!$                Do 190 L = 1, k
!!$                   If (B(j, L) .Ne. ZERO) Then
!!$                      TEMP = ALPHA * QCONJG (B(j, L))
!!$                      Do 180 i = 1, m
!!$                         C (i, j) = C (i, j) + TEMP * A (i, L)
!!$ 180                  Continue
!!$                   End If
!!$ 190            Continue
!!$ 200         Continue
!!$          Else
!!$             Do 250 j = 1, N
!!$                If (BETA .Eq. ZERO) Then
!!$                   Do 210 i = 1, m
!!$                      C (i, j) = ZERO
!!$ 210               Continue
!!$                Else If (BETA .Ne. ONE) Then
!!$                   Do 220 i = 1, m
!!$                      C (i, j) = BETA * C (i, j)
!!$ 220               Continue
!!$                End If
!!$                Do 240 L = 1, k
!!$                   If (B(j, L) .Ne. ZERO) Then
!!$                      TEMP = ALPHA * B (j, L)
!!$                      Do 230 i = 1, m
!!$                         C (i, j) = C (i, j) + TEMP * A (i, L)
!!$ 230                  Continue
!!$                   End If
!!$ 240            Continue
!!$ 250         Continue
!!$          End If
!!$       Else If (CONJA) Then
!!$          If (CONJB) Then
!!$             Do 280 j = 1, N
!!$                Do 270 i = 1, m
!!$                   TEMP = ZERO
!!$                   Do 260 L = 1, k
!!$                      TEMP = TEMP + CONJG (A(L, i)) * CONJG (B(j, L))
!!$ 260               Continue
!!$                   If (BETA .Eq. ZERO) Then
!!$                      C (i, j) = ALPHA * TEMP
!!$                   Else
!!$                      C (i, j) = ALPHA * TEMP + BETA * C (i, j)
!!$                   End If
!!$ 270            Continue
!!$ 280         Continue
!!$          Else
!!$             Do 310 j = 1, N
!!$                Do 300 i = 1, m
!!$                   TEMP = ZERO
!!$                   Do 290 L = 1, k
!!$                      TEMP = TEMP + CONJG (A(L, i)) * B (j, L)
!!$ 290               Continue
!!$                   If (BETA .Eq. ZERO) Then
!!$                      C (i, j) = ALPHA * TEMP
!!$                   Else
!!$                      C (i, j) = ALPHA * TEMP + BETA * C (i, j)
!!$                   End If
!!$ 300            Continue
!!$ 310         Continue
!!$          End If
!!$       Else
!!$          If (CONJB) Then
!!$             Do 340 j = 1, N
!!$                Do 330 i = 1, m
!!$                   TEMP = ZERO
!!$                   Do 320 L = 1, k
!!$                      TEMP = TEMP + A (L, i) * CONJG (B(j, L))
!!$ 320               Continue
!!$                   If (BETA .Eq. ZERO) Then
!!$                      C (i, j) = ALPHA * TEMP
!!$                   Else
!!$                      C (i, j) = ALPHA * TEMP + BETA * C (i, j)
!!$                   End If
!!$ 330            Continue
!!$ 340         Continue
!!$          Else
!!$             Do 370 j = 1, N
!!$                Do 360 i = 1, m
!!$                   TEMP = ZERO
!!$                   Do 350 L = 1, k
!!$                      TEMP = TEMP + A (L, i) * B (j, L)
!!$ 350               Continue
!!$                   If (BETA .Eq. ZERO) Then
!!$                      C (i, j) = ALPHA * TEMP
!!$                   Else
!!$                      C (i, j) = ALPHA * TEMP + BETA * C (i, j)
!!$                   End If
!!$ 360            Continue
!!$ 370         Continue
!!$          End If
!!$       End If
!!$       Return
!!$    End Subroutine ZGEMMQ
!!$
!!$
!!$       SUBROUTINE ZDSCALQ(N,DA,ZX,INCX)
!!$       real (Kind=DEF_QUAD_PREC) :: DA
!!$       INTEGER INCX,N
!!$       Complex (Kind=DEF_QUAD_PREC) ::ZX(*)
!!$       INTEGER I,IX
!!$       IF (N.LE.0 .OR. INCX.LE.0) RETURN
!!$       IF (INCX.EQ.1) GO TO 20
!!$       IX = 1
!!$       DO 10 I = 1,N
!!$           ZX(IX) = CMPLX(DA,0.0q0)*ZX(IX)
!!$           IX = IX + INCX
!!$    10 CONTINUE
!!$       RETURN
!!$    20 DO 30 I = 1,N
!!$           ZX(I) = CMPLX(DA,0.0q0)*ZX(I)
!!$    30 CONTINUE
!!$       RETURN
!!$       END SUBROUTINE ZDSCALQ
!!$
!!$
!!$          SUBROUTINE ZAXPYQ(N,ZA,ZX,INCX,ZY,INCY)
!!$          Complex (Kind=DEF_QUAD_PREC) ::ZA
!!$          INTEGER INCX,INCY,N
!!$          Complex (Kind=DEF_QUAD_PREC) :: ZX(*),ZY(*)
!!$          INTEGER I,IX,IY
!!$          real(Kind=DEF_QUAD_PREC),external::DCABS1Q
!!$
!!$          IF (N.LE.0) RETURN
!!$          IF (DCABS1Q(ZA) .EQ. 0.0q0) RETURN
!!$          IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!!$          IX = 1
!!$          IY = 1
!!$          IF (INCX.LT.0) IX = (-N+1)*INCX + 1
!!$          IF (INCY.LT.0) IY = (-N+1)*INCY + 1
!!$          DO 10 I = 1,N
!!$              ZY(IY) = ZY(IY) + ZA*ZX(IX)
!!$              IX = IX + INCX
!!$              IY = IY + INCY
!!$       10 CONTINUE
!!$          RETURN
!!$       20 DO 30 I = 1,N
!!$              ZY(I) = ZY(I) + ZA*ZX(I)
!!$       30 CONTINUE
!!$          RETURN
!!$          END SUBROUTINE ZAXPYQ
!!$
!!$           FUNCTION DCABS1Q(Z)
!!$               real (Kind=DEF_QUAD_PREC)::DCABS1Q
!!$               Complex (Kind=DEF_QUAD_PREC):: Z
!!$               DCABS1Q = ABS(real(Z,kind=16)) + ABS(IMAG(Z))
!!$               RETURN
!!$         END FUNCTION DCABS1Q
!!$
!!$         Subroutine ZGESVQ ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
!!$         INTEGER            INFO, LDA, LDB, N, NRHS
!!$          INTEGER            IPIV( * )
!!$          Complex (Kind=DEF_QUAD_PREC) ::  A( LDA, * ), B( LDB, * )
!!$
!!$          Complex (Kind=DEF_DBL_PREC) :: AR(N,N),BR(N,NRHS)
!!$
!!$          AR=A(1:N,1:N)
!!$          BR=B(1:N,1:NRHS)
!!$         call ZGESV ( N, NRHS, AR, LDA, IPIV, BR, LDB, INFO)
!!$         A(1:N,1:N)=AR
!!$         B(1:N,1:NRHS)=BR
!!$         End subroutine
