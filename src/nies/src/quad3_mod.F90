module quad3_mod
 use definitions
 contains

!!$*******************
!!$XXX    QUAD3   ****
!!$*******************
  function quad3 (n, x, y) 
    implicit none
    integer::n
    real (kind=prec)::x(:),y(:)
    real (kind=prec)::quad3    
!!$ Local
    integer::i,i0
    real (kind=prec)::rc1,rcz,rch,rc5,rc12
    real (kind=prec)::cm,cp,sq15
    real (kind=prec)::x1,x2,x3,x4,x5,x6,xt1,yt1,xt2,yt2
!!$---------------------------------------------------
!!$     QUADRATURE OF INTEGRAND Y(i) OF VARIABLE X(i)
!!$     (WHERE i=1,2,... N, AND N.GE.6)
!!$     USING 5TH DEGREE INTERPOLATION POLYNOMIAL
!!$---------------------------------------------------

    data rcz / 0.0D0 /, rch / 0.5D0 /, rc1 / 1.0D0 /, rc5 / 5.0D0 /, rc12 / 12.0D0 /

    sq15 = sqrt (rc1/rc5)
    cp = rch * (rc1+sq15)
    cm = rch * (rc1-sq15)

    quad3 = rcz
    do i = 2, n
       i0 = max (i-4, 0)
       i0 = min (i0, n-6)
       x1 = x (i0+1)
       x2 = x (i0+2)
       x3 = x (i0+3)
       x4 = x (i0+4)
       x5 = x (i0+5)
       x6 = x (i0+6)
!!$                             LAGRANGE INTERPOLATION
       xt1 = cp * x (i-1) + cm * x (i)
       yt1 = y (i0+1) * auxf2 (x1, x2, x3, x4, x5, x6, xt1) &
  &+ y (i0+2) * auxf2 (x2, x3, x4, x5, x6, x1, xt1) + y (i0+3) &
  & * auxf2 (x3,x4, x5, x6, x1, x2, xt1) + y (i0+4) * &
  & auxf2 (x4, x5, x6, x1, x2, x3, xt1) + y (i0+5) * &
  & auxf2 (x5, x6, x1, x2, x3, x4, xt1) + y (i0+6) * &
  & auxf2 (x6, x1, x2, x3, x4, x5, xt1)
!!$
       xt2 = cm * x (i-1) + cp * x (i)
       yt2 = y (i0+1) * auxf2 (x1, x2, x3, x4, x5, x6, xt2) + &
   & y (i0+2) * auxf2 (x2, x3, x4, x5, x6, x1, xt2) + y (i0+3) * &
   & auxf2 (x3, x4, x5, x6, x1, x2, xt2) + y (i0+4) * &
   & auxf2 (x4, x5, x6, x1, x2, x3, xt2) + y (i0+5) * & 
   & auxf2 (x5, x6, x1, x2, x3, x4, xt2) + y (i0+6) * &
   & auxf2 (x6, x1, x2, x3, x4, x5, xt2)
!!$                                       LOBATTO RULE
       quad3 = quad3 + (y(i-1)+rc5*(yt1+yt2)+y(i)) * (x(i)-x(i-1)) / rc12
    end do

!!$        RETURN
  end function quad3

   function auxf2 (t0, t1, t2, t3, t4, t5, tt)
      implicit none
      real (kind=prec)::auxf2,tt,t1,t2,t3,t0,t4,t5
      auxf2 = (tt-t1) * (tt-t2) * (tt-t3) * (tt-t4) * (tt-t5) / ((t0-t1)*(t0-t2)*(t0-t3)*(t0-t4)*(t0-t5))
   end function auxf2
 

!!$*******************
!!$XXX   RUKUST   ****
!!$*******************
 subroutine rukust (p, q, h, app, apq, aqp, aqq, nsi)
    implicit none
    integer::nsi
    real (kind=prec)::p,q,h,app (0:10), apq (0:10), aqp (0:10), aqq (0:10)
!!$ Local
    integer::isi,j
    real (kind=prec)::h1,h2,h6,pk1,pk2,pk3,pk4,pdum,qdum,qk1,qk2,qk3,qk4
    real (kind=prec)::rc2,rc6
!!$-----------------------------------------------------
!!$    ONE STEP OF RUNGE-KUTTA INTEGRATION OF A SYSTEM
!!$    OF TWO LINEAR DIFFERENTIAL EQUATIONS:
!!$        P'(X) = APP(X) * P(X) + APQ(X) * Q(X)
!!$        Q'(X) = AQP(X) * P(X) + AQQ(X) * Q(X)
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

    data rc2 / 2.0D0 /, rc6 / 6.0D0 /

    h1 = h / real (nsi,kind=prec)
    h2 = h1 / rc2
    h6 = h1 / rc6

    j = 0
    do isi = 1, nsi
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
    end do
    return
  end subroutine rukust

!!$*******************
!!$XXX    rapo    ****
!!$*******************
 subroutine rapo (nr, ws, r)
!!$
!!$-------------------------------------
!!$   GENERATES RADIAL POINTS
!!$-------------------------------------
!!$     R(1)=0.0 , R(2)=H , R(NR)=WS ,
!!$     R(i) = H * (i-1) * Q**(i-2)
!!$     H - A SMALL QUANTITY, H=2.0E-5
!!$     NR - NUMBER OF POINTS
!!$     WS - WIGNER-SEITZ RADIUS
!!$------------------------------------- 
    implicit none
    real (kind=prec)::r(:),ws
    integer::nr
!!$ Local    
    integer::i
    real (kind=prec)::rcz,h,hq,alpha

    data rcz / 0.0D0 /, h / 2.0D-5 /

    alpha= exp(log (ws / (h*real(nr-1,kind=prec))) / real (nr-2,kind=prec))

    r (1) = rcz
    r (2) = h
    hq = h
    do i = 3, nr - 1
       hq = hq * alpha
       r (i) = real (i-1,kind=prec) * hq
    end do
    r (nr) = ws
    return
  end subroutine rapo

!!$*******************
!!$XXX    LIPO    ****
!!$*******************
  subroutine lipo (az,nr,r,v,vi,nsirk)
    implicit none
    integer:: nr,nsirk
    real (kind=prec)::r(:),v(:),vi(:,:),az
!!$ Local
    integer::i,j,n2,i0
    real (kind=prec)::rc1,un2,rc2
    real (kind=prec)::x1,x2,x3,x4,y1,y2,y3,y4,xt,yt

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
!!$     VI(.,.) - INTERPOLATED POTENTIAL:  THE VALUES VI(J,i),
!!$               WHERE J=1,2,... 2*NSIRK-1, AND i=3,4,... NR,
!!$               REFER TO RADIAL POINTS SAMPLING UNIFORMLY
!!$               THE INTERVAL (R(i-1),R(i)).
!!$-----------------------------------------------------------
    data rc1 / 1.0D0 /, rc2 / 2.0D0 /

    n2 = 2 * nsirk
    un2 = rc1 / real (n2,kind=prec)
    
    do i = 3, nr
       i0 = i - 3
       if (i .eq. nr) i0 = i - 4
       x1 = r (i0+1)
       x2 = r (i0+2)
       x3 = r (i0+3)
       x4 = r (i0+4)
       if (i .eq. 3) then
          y1 = - rc2 * az
       else
          y1 = x1 * v (i0+1)
       end if
       y2 = x2 * v (i0+2)
       y3 = x3 * v (i0+3)
       y4 = x4 * v (i0+4)
!!$                             LAGRANGE INTERPOLATION
       do  j = 1, n2 - 1
          xt = (r(i-1)*real(n2-j,kind=prec)+r(i)*real(j,kind=prec)) * un2
          yt = y1 * auxf1 (x1, x2, x3, x4, xt) + y2 * auxf1 (x2, x3, x4, x1, xt) + y3 * auxf1 (x3, x4, x1, x2, xt) + & 
      &  y4 * auxf1 (x4,  x1, x2, x3, xt)
          vi (j, i) = yt / xt
       end do
    end do
    return
  end subroutine lipo

   function auxf1 (t0, t1, t2, t3, tt)
      implicit none
      real (kind=prec)::auxf1,tt,t1,t2,t3,t0    
      auxf1 = (tt-t1) * (tt-t2) * (tt-t3) / ((t0-t1)*(t0-t2)*(t0-t3))
    end function auxf1


end module quad3_mod