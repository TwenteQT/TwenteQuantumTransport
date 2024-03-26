!!$ $Id:$
#include "math_def.h"

!!$
Module heritage
      Implicit None
Contains
!!$ Herritage from Ilja bellow. Removed global variables, common blocks and static memory allocations

!!$*******************
!!$XXX    RAPO    ****
!!$*******************
      Subroutine rapo (nr, ws, r)
!!$
!!$-------------------------------------
!!$   GENERATES RADIAL POINTS
!!$-------------------------------------
!!$     R(1)=0.0 , R(2)=H , R(NR)=WS ,
!!$     R(I) = H * (I-1) * Q**(I-2)
!!$     H - A SMALL QUANTITY, H=2.0E-5
!!$     NR - NUMBER OF POINTS
!!$     WS - WIGNER-SEITZ RADIUS
!!$-------------------------------------
         Implicit None
         Real (Kind=DEF_DBL_PREC) :: r (:), ws
         Integer :: nr
!!$ Local
         Integer :: i
         Real (Kind=DEF_DBL_PREC) :: rcz, h, hq, alpha

         Data rcz / 0.0D0 /, h / 2.0D-5 /

         alpha = Exp (Log(ws/(h*real(nr-1, kind=DEF_DBL_PREC)))/real(nr-2, kind=DEF_DBL_PREC))

         r (1) = rcz
         r (2) = h
         hq = h
         Do i = 3, nr - 1
            hq = hq * alpha
            r (i) = real (i-1, kind=DEF_DBL_PREC) * hq
         End Do
         r (nr) = ws
         Return
      End Subroutine rapo

!!$*******************
!!$XXX    LIPO    ****
!!$*******************
      Subroutine lipo (az, nsirk, nr, r, v, vi)
         Implicit None
         Integer :: nr, nsirk
         Real (Kind=DEF_DBL_PREC) :: r (:), v (:), vi (:, :), az
!!$ Local
         Integer :: i, j, n2, i0
         Real (Kind=DEF_DBL_PREC) :: rc1, un2, rc2
         Real (Kind=DEF_DBL_PREC) :: x1, x2, x3, x4, y1, y2, y3, y4, xt, yt

!!$-----------------------------------------------------------
!!$   LAGRANGE INTERPOLATION OF POTENTIAL FOR THE RUNGE-KUTTA
!!$   INTEGRATION OF THE RADIAL SCHROEDINGER EQUATION
!!$-----------------------------------------------------------
!!$  INPUT:
!!$     AZ - ATOMIC NUMBER
!!$     NR - SIZE OF RADIAL MESH
!!$     R(.) - RADIAL MESH
!!$     V(.) - POTENTIAL
!!$  OUTPUT:
!!$     VI(.,.) - INTERPOLATED POTENTIAL:  THE VALUES VI(J,I),
!!$               WHERE J=1,2,... 2*NSIRK-1, AND I=3,4,... NR,
!!$               REFER TO RADIAL POINTS SAMPLING UNIFORMLY
!!$               THE INTERVAL (R(I-1),R(I)).
!!$-----------------------------------------------------------
         Data rc1 / 1.0D0 /, rc2 / 2.0D0 /

         n2 = 2 * nsirk
         un2 = rc1 / real (n2, kind=DEF_DBL_PREC)

         Do i = 3, nr
            i0 = i - 3
            If (i .Eq. nr) i0 = i - 4
            x1 = r (i0+1)
            x2 = r (i0+2)
            x3 = r (i0+3)
            x4 = r (i0+4)
            If (i .Eq. 3) Then
               y1 = - rc2 * az
            Else
               y1 = x1 * v (i0+1)
            End If
            y2 = x2 * v (i0+2)
            y3 = x3 * v (i0+3)
            y4 = x4 * v (i0+4)
!!$                             LAGRANGE INTERPOLATION
            Do j = 1, n2 - 1
               xt = (r(i-1)*real(n2-j, kind=DEF_DBL_PREC)+r(i)*real(j, kind=DEF_DBL_PREC)) * un2
               yt = y1 * auxf (x1, x2, x3, x4, xt) + y2 * auxf (x2, x3, x4, x1, xt) + y3 * auxf (x3, x4, x1, &
              & x2, xt) + y4 * auxf (x4, x1, x2, x3, xt)
               vi (j, i) = yt / xt
            End Do
         End Do
         Return
      Contains
         Function auxf (t0, t1, t2, t3, tt)
            Implicit None
            Real (Kind=DEF_DBL_PREC) :: auxf, tt, t1, t2, t3, t0
            auxf = (tt-t1) * (tt-t2) * (tt-t3) / ((t0-t1)*(t0-t2)*(t0-t3))
         End Function auxf
      End Subroutine lipo

!!$*******************
!!$XXX    QUAD3   ****
!!$*******************
      Function quad3 (n, x, y)
         Implicit None
         Integer :: n
         Real (Kind=DEF_DBL_PREC), intent(in) :: x (:), y (:)
         Real (Kind=DEF_DBL_PREC) :: quad3
!!$ Local
         Integer :: i, i0
         Real (Kind=DEF_DBL_PREC) :: rc1, rcz, rch, rc5, rc12
         Real (Kind=DEF_DBL_PREC) :: cm, cp, sq15
         Real (Kind=DEF_DBL_PREC) :: x1, x2, x3, x4, x5, x6, xt1, yt1, xt2, yt2
!!$---------------------------------------------------
!!$     QUADRATURE OF INTEGRAND Y(I) OF VARIABLE X(I)
!!$     (WHERE I=1,2,... N, AND N.GE.6)
!!$     USING 5TH DEGREE INTERPOLATION POLYNOMIAL
!!$---------------------------------------------------

         Data rcz / 0.0D0 /, rch / 0.5D0 /, rc1 / 1.0D0 /, rc5 / 5.0D0 /, rc12 / 12.0D0 /

         sq15 = Sqrt (rc1/rc5)
         cp = rch * (rc1+sq15)
         cm = rch * (rc1-sq15)

         quad3 = rcz
         Do i = 2, n
            i0 = Max (i-4, 0)
            i0 = Min (i0, n-6)
            x1 = x (i0+1)
            x2 = x (i0+2)
            x3 = x (i0+3)
            x4 = x (i0+4)
            x5 = x (i0+5)
            x6 = x (i0+6)
!!$                             LAGRANGE INTERPOLATION
            xt1 = cp * x (i-1) + cm * x (i)
            yt1 = y (i0+1) * auxf (x1, x2, x3, x4, x5, x6, xt1) + y (i0+2) * auxf (x2, x3, x4, x5, x6, x1, &
           & xt1) + y (i0+3) * auxf (x3, x4, x5, x6, x1, x2, xt1) + y (i0+4) * auxf (x4, x5, x6, x1, x2, x3, &
           & xt1) + y (i0+5) * auxf (x5, x6, x1, x2, x3, x4, xt1) + y (i0+6) * auxf (x6, x1, x2, x3, x4, x5, &
           & xt1)
!!$
            xt2 = cm * x (i-1) + cp * x (i)
            yt2 = y (i0+1) * auxf (x1, x2, x3, x4, x5, x6, xt2) + y (i0+2) * auxf (x2, x3, x4, x5, x6, x1, &
           & xt2) + y (i0+3) * auxf (x3, x4, x5, x6, x1, x2, xt2) + y (i0+4) * auxf (x4, x5, x6, x1, x2, x3, &
           & xt2) + y (i0+5) * auxf (x5, x6, x1, x2, x3, x4, xt2) + y (i0+6) * auxf (x6, x1, x2, x3, x4, x5, &
           & xt2)
!!$                                       LOBATTO RULE
            quad3 = quad3 + (y(i-1)+rc5*(yt1+yt2)+y(i)) * (x(i)-x(i-1)) / rc12
         End Do

!!$        RETURN
      Contains
         Function auxf (t0, t1, t2, t3, t4, t5, tt)
            Implicit None
            Real (Kind=DEF_DBL_PREC) :: auxf, tt, t1, t2, t3, t0, t4, t5
            auxf = (tt-t1) * (tt-t2) * (tt-t3) * (tt-t4) * (tt-t5) / &
           & ((t0-t1)*(t0-t2)*(t0-t3)*(t0-t4)*(t0-t5))
         End Function auxf
      End Function quad3

!*******************
!XXX    PRIM3   ****
!*******************
      Subroutine PRIM3 (n, x, y, F)
!---------------------------------------------------------------------------
!!$ CALCULATES INDEFINITE INTEGRAL F(i) FOR INTEGRAND Y(i) OF VARIABLE X(i)
!!$ (WHERE i=1,2,... N, AND N.GE.6) USING 5TH DEGREE INTERPOLATION POLYNOMIAL
!---------------------------------------------------------------------------
         Implicit None
         Real (Kind=DEF_DBL_PREC), Parameter :: rcz = 0.0d0, rch = 0.5d0, rc1 = 1.0d0, rc5 = 5.0d0, rc12 = &
        & 12.0d0
         Real (Kind=DEF_DBL_PREC) :: x (n), y (n), F (n), sq15, cp, cm, xt1, yt1, xt2, yt2, x1, x2, x3, x4, &
        & x5, x6
         Integer :: i, i0, n

         sq15 = dSQRT (rc1/rc5)
         cp = rch * (rc1+sq15)
         cm = rch * (rc1-sq15)

         F (1) = rcz
         Do i = 2, n
            i0 = Max (i-4, 0)
            i0 = Min (i0, n-6)
            x1 = x (i0+1)
            x2 = x (i0+2)
            x3 = x (i0+3)
            x4 = x (i0+4)
            x5 = x (i0+5)
            x6 = x (i0+6)

!!$                  LAGRANGE INTERPOLATION
            xt1 = cp * x (i-1) + cm * x (i)
            yt1 = y (i0+1) * AUXF2 (x1, x2, x3, x4, x5, x6, xt1) + y (i0+2) * AUXF2 (x2, x3, x4, x5, x6, x1, &
           & xt1) + y (i0+3) * AUXF2 (x3, x4, x5, x6, x1, x2, xt1) + y (i0+4) * AUXF2 (x4, x5, x6, x1, x2, &
           & x3, xt1) + y (i0+5) * AUXF2 (x5, x6, x1, x2, x3, x4, xt1) + y (i0+6) * AUXF2 (x6, x1, x2, x3, &
           & x4, x5, xt1)

            xt2 = cm * x (i-1) + cp * x (i)
            yt2 = y (i0+1) * AUXF2 (x1, x2, x3, x4, x5, x6, xt2) + y (i0+2) * AUXF2 (x2, x3, x4, x5, x6, x1, &
           & xt2) + y (i0+3) * AUXF2 (x3, x4, x5, x6, x1, x2, xt2) + y (i0+4) * AUXF2 (x4, x5, x6, x1, x2, &
           & x3, xt2) + y (i0+5) * AUXF2 (x5, x6, x1, x2, x3, x4, xt2) + y (i0+6) * AUXF2 (x6, x1, x2, x3, &
           & x4, x5, xt2)
!!$                 LOBATTO RULE
            F (i) = F (i-1) + (y(i-1)+rc5*(yt1+yt2)+y(i)) * (x(i)-x(i-1)) / rc12
         End Do
         Return
      Contains
         Function AUXF2 (t0, t1, t2, t3, t4, t5, tt)
            Implicit None
            Real (Kind=DEF_DBL_PREC) :: AUXF2, tt, t1, t2, t3, t0, t4, t5
            AUXF2 = (tt-t1) * (tt-t2) * (tt-t3) * (tt-t4) * (tt-t5) / &
           & ((t0-t1)*(t0-t2)*(t0-t3)*(t0-t4)*(t0-t5))
         End Function AUXF2
      End Subroutine PRIM3


!!$*******************
!!$XXX   RUKUST   ****
!!$*******************
      Subroutine rukust (p, q, h, app, apq, aqp, aqq, nsi)
         Implicit None
         Integer :: nsi
         Real (Kind=DEF_DBL_PREC) :: p, q, h, app (0:10), apq (0:10), aqp (0:10), aqq (0:10)
!!$ Local
         Integer :: isi, j
         Real (Kind=DEF_DBL_PREC) :: h1, h2, h6, pk1, pk2, pk3, pk4, pdum, qdum, qk1, qk2, qk3, qk4
         Real (Kind=DEF_DBL_PREC) :: rc2, rc6
!!$-----------------------------------------------------
!!$    ONE STEP OF RUNGE-KUTTA INTEGRATION OF A SYSTEM
!!$    OF TWO LINEAR DIFFERENTIAL EQUATIONS:
!!$        P^prime(X) = APP(X) * P(X) + APQ(X) * Q(X)
!!$        Q^prime(X) = AQP(X) * P(X) + AQQ(X) * Q(X)
!!$-----------------------------------------------------
!!$  ON INPUT: P,Q - INITIAL VALUES
!!$            H  -  INCREMENT OF X
!!$            NSI - NO. OF SUBINTERVALS (1.LE.NSI.LE.5)
!!$            APP(J), APQ(J), AQP(J), AQQ(J)
!!$              (J = 0, ... , 2*NSI) - VALUES OF THE
!!$              COEFFICIENTS FOR ARGUMENTS X SAMPLING
!!$              UNIFORMLY THE WHOLE STEP OF LENGTH H
!!$-----------------------------------------------------
!!$  ON OUTPUT: P,Q - FINAL VALUES
!!$-----------------------------------------------------

         Data rc2 / 2.0D0 /, rc6 / 6.0D0 /

         h1 = h / real (nsi, kind=DEF_DBL_PREC)
         h2 = h1 / rc2
         h6 = h1 / rc6

         j = 0
         Do isi = 1, nsi
            pk1 = app (j) * p + apq (j) * q
            qk1 = aqp (j) * p + aqq (j) * q
            j = j + 1
            pdum = p + h2 * pk1
            qdum = q + h2 * qk1
            pk2 = app (j) * pdum + apq (j) * qdum
            qk2 = aqp (j) * pdum + aqq (j) * qdum
            pdum = p + h2 * pk2
            qdum = q + h2 * qk2
            pk3 = app (j) * pdum + apq (j) * qdum
            qk3 = aqp (j) * pdum + aqq (j) * qdum
            j = j + 1
            pdum = p + h1 * pk3
            qdum = q + h1 * qk3
            pk4 = app (j) * pdum + apq (j) * qdum
            qk4 = aqp (j) * pdum + aqq (j) * qdum
            p = p + h6 * (pk1+rc2*(pk2+pk3)+pk4)
            q = q + h6 * (qk1+rc2*(qk2+qk3)+qk4)
         End Do
         Return
      End Subroutine rukust

!!$*******************
!!$XXX    RSEV    ****
!!$*******************
      Subroutine rsev (az, e, e0, l, irel, nsirk, nr, r, v, vi, wf, wg)
         Implicit None
         Integer :: nr, l, irel, nsirk
         Real (Kind=DEF_DBL_PREC) :: az, e, e0
         Real (Kind=DEF_DBL_PREC) :: r (:), v (:), vi (:, :), wf (:), wg (:)
!!$ Local
         Integer :: i, j, nrst, nsi2
         Real (Kind=DEF_DBL_PREC) :: rcz, rc1, rc2, rc3
         Real (Kind=DEF_DBL_PREC) :: p (nr), q (nr), wrk (nr)
         Real (Kind=DEF_DBL_PREC) :: tapp (0:10), tapq (0:10), taqp (0:10), taqq (0:10), tr (0:10), tv (0:10)
         Real (Kind=DEF_DBL_PREC) :: al, alp1, blam, a1, a2, sum, cnorm, dum, c274, rtb, tq, tp, th, a1mb, &
        & a1pb, beta
         Real (Kind=DEF_DBL_PREC) :: b2, b1, u0, un2, ucsq, uc, twoz

!!$---------------------------------------------------------
!!$     SOLUTION OF RADIAL SCHROEDINGER EQUATION
!!$     FOR VALENCE ELECTRONS:
!!$     BOTH NON-RELATIVISTIC (IREL=0) AND
!!$     SCALAR-RELATIVISTIC (IREL=1) VERSION
!!$---------------------------------------------------------
!!$  INPUT:
!!$     AZ - ATOMIC NUMBER
!!$     E - ENERGY
!!$     E0 - ENERGY FOR DOWNFOLDING OF SMALL COMPONENTS
!!$     L - ORBITAL QUANTUM NUMBER
!!$     IREL - RELATIVITY
!!$     NR - SIZE OF RADIAL MESH
!!$     R(.) - RADIAL MESH
!!$     V(.) - POTENTIAL
!!$     VI(.,.) - INTERPOLATED POTENTIAL
!!$  OUTPUT:
!!$     WG(.) - WAVE FUNCTION NORMALIZED TO UNITY
!!$             (LARGE COMPONENT)
!!$     WF(.) - FUNCTION RELATED TO RADIAL DERIVATIVE OF WG
!!$             (SMALL COMPONENT)
!!$---------------------------------------------------------

         Data nrst / 10 /, c274 / 274.074D0 /

         Data rcz / 0.0D0 /, rc1 / 1.0D0 /, rc2 / 2.0D0 /, rc3 / 3.0D0 /

!!$                             SET UP CONSTANTS
         al = real (l, kind=DEF_DBL_PREC)
         alp1 = real (l+1, kind=DEF_DBL_PREC)
         blam = al * alp1
         twoz = rc2 * az
         uc = real (irel, kind=DEF_DBL_PREC) / c274
         ucsq = uc ** 2
!!$          write(*,*) uc, ucsq
         nsi2 = 2 * nsirk
         un2 = rc1 / real (nsi2, kind=DEF_DBL_PREC)
!!$                            START OF INTEGRATION
         p (1) = rcz
         q (1) = rcz
!!$                              NON-RELATIVISTIC
         If (irel .Eq. 0) Then
            u0 = v (2) + twoz / r (2) - e
            a1 = - az / alp1
            a2 = (u0-twoz*a1) / (rc2*(rc2*al+rc3))
            b1 = - az
            b2 = a2 * (al+rc2)
            p (2) = rc1 + r (2) * (a1+r(2)*a2)
            q (2) = b1 + r (2) * b2
            beta = rc1
            If (l .Gt. 0) Then
               p (2) = r (2) * p (2)
               q (2) = al + r (2) * q (2)
               beta = al
            End If
         End If
!!$                           SCALAR-RELATIVISTIC
         If (irel .Eq. 1) Then
            beta = Sqrt (blam+rc1-(twoz*uc)**2)
            p (2) = rc1
            q (2) = (beta-rc1) / (twoz*ucsq)
         End If

!!$                       STARTING RUNGE-KUTTA INTEGRATION
         If (nrst /= 2) Then
            a1pb = rc1 + beta
            a1mb = rc1 - beta

            Do i = 3, nrst

               Do j = 0, nsi2
                  tr (j) = (r(i-1)*real(nsi2-j, kind=DEF_DBL_PREC)+r(i)*real(j, kind=DEF_DBL_PREC)) * un2
               End Do
               tv (0) = v (i-1)
               Do j = 1, nsi2 - 1
                  tv (j) = vi (j, i)
               End Do
               tv (nsi2) = v (i)

               Do j = 0, nsi2
                  tapp (j) = a1mb / tr (j)
               End Do
               Do j = 0, nsi2
                  tapq (j) = rc1 + ucsq * (e0-tv(j))
               End Do
               Do j = 0, nsi2
                  taqp (j) = blam / (tapq(j)*tr(j)**2) + tv (j) - e
               End Do
               Do j = 0, nsi2
                  taqq (j) = - a1pb / tr (j)
               End Do

               th = r (i) - r (i-1)
               tp = p (i-1)
               tq = q (i-1)
               Call rukust (tp, tq, th, tapp, tapq, taqp, taqq, nsirk)
               p (i) = tp
               q (i) = tq
            End Do
         End If

         Do i = 2, nrst
            rtb = r (i) ** beta
            p (i) = rtb * p (i)
            q (i) = rtb * q (i)
         End Do

!!$                             RUNGE-KUTTA INTEGRATION
         Do i = nrst + 1, nr

            Do j = 0, nsi2
               tr (j) = (r(i-1)*real(nsi2-j, kind=DEF_DBL_PREC)+r(i)*real(j, kind=DEF_DBL_PREC)) * un2
            End Do
            tv (0) = v (i-1)
            Do j = 1, nsi2 - 1
               tv (j) = vi (j, i)
            End Do

            tv (nsi2) = v (i)

            Do j = 0, nsi2
               tapp (j) = rc1 / tr (j)
            End Do

            Do j = 0, nsi2
               tapq (j) = rc1 + ucsq * (e0-tv(j))
            End Do

            Do j = 0, nsi2
               taqp (j) = blam / (tapq(j)*tr(j)**2) + tv (j) - e
            End Do

            Do j = 0, nsi2
               taqq (j) = - rc1 / tr (j)
            End Do

            th = r (i) - r (i-1)
            tp = p (i-1)
            tq = q (i-1)
            Call rukust (tp, tq, th, tapp, tapq, taqp, taqq, nsirk)
            p (i) = tp
            q (i) = tq
         End Do

!!$                                  NORMALIZATION
         Do i = 1, nr
            wrk (i) = p (i) ** 2
         End Do

         sum = quad3 (nr, r, wrk)
         cnorm = Sqrt (sum)

         wg (1) = rcz
         wf (1) = rcz
         If (irel .Eq. 0) Then
            If (l .Eq. 0) Then
               wg (1) = rc1 / cnorm
               wf (1) = - az / cnorm
            End If
            If (l .Eq. 1) wf (1) = rc1 / cnorm
         End If
         Do i = 2, nr
            dum = cnorm * r (i)
            wg (i) = p (i) / dum
            wf (i) = q (i) / dum
         End Do
         Return
      End Subroutine rsev
End Module heritage
