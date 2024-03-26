Module bzgrid
!!$    This module contains procedures for generating k-mesh
!!$    for BZ integration copied from I.Turek LMTO code.
      Implicit None

Contains

!!$*******************
!!$XXX    GIBZ    ****
!!$*******************
      Subroutine gibz (bz, vbr, vbg, nk, nsym, inve)
         Use definitions
!!$****************************************************************
!!$  GENERATES NETWORK OF K||-POINTS IN IRREDUCIBLE BRILLOUIN ZONE
!!$****************************************************************
         Implicit None

         Type (bzone) :: bz
         Integer :: nk, nsym, inve, iftr
         Integer :: nbk1 ((2*nk+1)**2), nbk2 ((2*nk+1)**2), istart !,  ext=0
         Real (Kind=8) :: akbz (2, (2*nk+1)**2), wkbz ((2*nk+1)**2), vbr (2, 2), vbg (2, 2)!, k_pos(2)
         Real (Kind=8) :: ajm, rnk, anu3, dky, dkx, q, dk, lf, bt, hg, wd, dx, dy
         Real (Kind=8) :: p, fx2, fx1, tw, pi
         Integer :: nbz, iy, nu, ndy, ndx, nmy, nmx, i, j, i2, i1, nk1, nk2, ibz, mnbz

         pi = 4.0d0 * datan (1.0d0)
         mnbz = (2*nk+1) ** 2

         If (nsym .Eq. 10) Go To 888
         If (nsym .Eq. 11) Go To 988
         If (nsym .Eq.-1) Go To 1088
         If (nsym .Gt. 0) Go To 2000

!!$-------------------------- NSYM=0: GENERAL CASE

!!$ SHifT THE K|| MESH A LITTER BIT IN ORDER TO AVOID
!!$ THE KX=0 OR KY=0
         ibz = 0
         tw = 0.0d0
         nk1 = nk
         nk2 = nk

         If (iftr == 0) Then
            istart = - nk1 + 1
         Else
            istart = 1
         End If

         Do i1 = istart, nk1
!!$      DO 300 I1=1,NK1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - nk2 + 1, nk2
!!$      DO 300 I2=1,NK2
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Go To 290

2000     Go To (201, 202, 203, 204), nsym

201      Go To (2011, 2012), inve

!!$--------  NSYM = 1 - SYMMETRY OF CUBIC (001) SURFACES:
!!$                     SYMMETRY GROUP GIVEN BY 4-FOLD ROTATION
!!$                     AXIS Z AND BY MIRROR PLANE X-Z
!!$
!!$------------------------------ NSYM=1, INVE=1
!!$                     INVE = 1: THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,0)
!!$                                      A2 = (0,P)
!!$                               WHERE P is POSITIVE.
2011     ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         dk = pi / (p*real(nk))
         Do i = 1, nk
            Do j = 1, i
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               wkbz (ibz) = 1.0d0
               If (i .Eq. j) wkbz (ibz) = 0.5d0
               tw = tw + wkbz (ibz)
               nbk1 (ibz) = i
               nbk2 (ibz) = j
               akbz (1, ibz) = (real(i)-0.5d0) * dk
               akbz (2, ibz) = (real(j)-0.5d0) * dk
            End Do
         End Do

         Go To 290
!!$------------------------------ NSYM=1, INVE=2
!!$
!!$                     INVE = 2: THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,P)
!!$                                      A2 = (-P,P)*U
!!$                               WHERE P is POSITIVE
!!$                               AND U=1 OR U=-1.
2012     ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         dk = pi / (2.0d0*p*real(nk))
         Do i = 1, nk
            Do j = 1, i
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               wkbz (ibz) = 1.0d0
               If (i .Eq. j) wkbz (ibz) = 0.5d0
               tw = tw + wkbz (ibz)
               nbk1 (ibz) = i
               nbk2 (ibz) = j
               akbz (1, ibz) = real (i-j) * dk !+1.0D-4
               akbz (2, ibz) = real (i+j-1) * dk
            End Do
         End Do

         Go To 290

202      Go To (2021, 2022), inve

!!$--------  NSYM = 2 - SYMMETRY OF CUBIC (110) SURFACES:
!!$                     SYMMETRY GROUP GIVEN BY TWO
!!$                     MIRROR PLANES X-Z AND Y-Z
!!$
!!$------------------------------ NSYM=2, INVE=1
!!$
!!$                    INVE = 1: THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,0)
!!$                                      A2 = (0,Q)
!!$                               WHERE P AND Q ARE POSITIVE.

2021     ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         q = vbr (2, 2)
         nmx = 1
         nmy = 1
         If (p .Lt. (0.75d0*q)) nmx = 2
         If (q .Lt. (0.75d0*p)) nmy = 2
         If (p .Lt. (0.4d0*q)) write (IW6, '(" **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ")')
         If (q .Lt. (0.4d0*p)) write (IW6, '(" **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ")')
         flush(IW6)
         ndx = nk * nmx
         ndy = nk * nmy
         dkx = pi / (p*real(ndx))
         dky = pi / (q*real(ndy))
         Do i = 1, ndx
            Do j = 1, ndy
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               nbk1 (ibz) = i
               nbk2 (ibz) = j
               akbz (1, ibz) = (real(i)-0.5d0) * dkx
               akbz (2, ibz) = (real(j)-0.5d0) * dky
            End Do
         End Do

         Go To 290
!!$------------------------------ NSYM=2, INVE=2
!!$                     INVE = 2: THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,Q)
!!$                                      A2 = (-P,Q)*U
!!$                               WHERE P AND Q ARE POSITIVE
!!$                               AND U=1 OR U=-1.

2022     ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         q = vbr (2, 1)
         nmx = 1
         nmy = 1
         If (p .Lt. (0.75d0*q)) nmx = 2
         If (q .Lt. (0.75d0*p)) nmy = 2
         If (p .Lt. (0.4d0*q)) write (IW6, '(" **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ")')
         If (q .Lt. (0.4d0*p)) write (IW6, '(" **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ")')
         flush(IW6)
         ndx = nk * nmx
         ndy = nk * nmy
         dkx = pi / (p*real(ndx))
         dky = pi / (q*real(ndy))
         Do i = 1, ndx
            Do j = 1, ndy
               If ((2*j-1)*nmx .Gt. (2*i-1)*nmy) Cycle
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               wkbz (ibz) = 1.0d0
               If (i .Eq. j .And. nmx .Eq. nmy) wkbz (ibz) = 0.5d0
               tw = tw + wkbz (ibz)
               nbk1 (ibz) = i
               nbk2 (ibz) = j
               akbz (1, ibz) = (real(i)-0.56d0) * dkx
               akbz (2, ibz) = (real(j)-0.47d0) * dky
               If (i .Eq. j) Then
                  akbz (1, ibz) = (real(i)-0.42d0) * dkx
                  akbz (2, ibz) = (real(j)-0.58d0) * dky
               End If
               If (i .Eq. ndx) Then
                  akbz (1, ibz) = (real(i)-0.52d0) * dkx
               End If
               If (j .Eq. ndy) Then
                  akbz (2, ibz) = (real(j)-0.43d0) * dkx
               End If
            End Do
         End Do

         Go To 290
!!$------------------------------ NSYM=3
!!$--------  NSYM = 3 - SYMMETRY OF CUBIC (111)  AND
!!$                     HEXAGONAL (0001) SURFACES:
!!$                     SYMMETRY GROUP GIVEN BY 3-FOLD ROTATION
!!$                     AXIS Z AND BY MIRROR PLANE Y-Z
!!$
!!$          INVE - NOT ACTIVE.   THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,0)
!!$                                   A2 = (P/2)*(1,S)*U
!!$                               WHERE P is POSITIVE,
!!$                               S=SQRT(3) OR S=-SQRT(3),
!!$                               AND U=1 OR U=-1.
!!$
203      ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         dkx = 2.0d0 * pi / (p*real(nk))
         dky = dkx / Sqrt (3.0d0)
         Do nu = 1, 2
            anu3 = real (nu) / 3.0d0
            Do i = 1, nk
               Do j = i, nk
                  iy = i + 2 * j + nu - 3
                  If (iy .Gt. nk) Cycle
                  ibz = ibz + 1
!!$              if (ibz .gt. mnbz) goto 299
                  ajm = 1.0d0
                  If (i .Eq. j .Or. iy .Eq. nk) ajm = 2.0d0
                  If (i .Eq. j .And. iy .Eq. nk) ajm = 6.0d0
                  wkbz (ibz) = 1.0d0 / ajm
                  tw = tw + wkbz (ibz)
                  nbk1 (ibz) = i
                  nbk2 (ibz) = j
                  akbz (1, ibz) = (real(i-1)+anu3) * dkx
                  akbz (2, ibz) = (real(iy)) * dky
               End Do
            End Do
         End Do

         Go To 290

204      Go To (2041, 2042), inve

!!$--------  NSYM = 4 - SYMMETRY OF CUBIC (N10) SURFACES
!!$                     AND TILT GRAIN BOUNDARIES:
!!$                     SYMMETRY GROUP GIVEN BY MIRROR PLANE X-Z
!!$
!!$------------------------------ NSYM=4, INVE=1
!!$                     INVE = 1: THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,0)
!!$                                      A2 = (0,Q)
!!$                               WHERE P AND Q ARE POSITIVE.

2041     ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         q = vbr (2, 2)
         rnk = real (nk)
         ndx = nk
         ndy = nk
         If (p .Gt. (1.2d0*q)) ndy = Int (p*rnk/q)
         If (q .Gt. (1.2d0*p)) ndx = Int (q*rnk/p)
         dkx = pi / (p*real(ndx))
         dky = pi / (q*real(ndy))
         Do i = 1, ndx
            Do j = 1, ndy
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               nbk1 (ibz) = i
               nbk2 (ibz) = j
               akbz (1, ibz) = (real(i)-0.5d0) * dkx
               akbz (2, ibz) = (real(j)-0.5d0) * dky
            End Do
         End Do


         Go To 290
!!$------------------------------ NSYM=4, INVE=2
!!$                     INVE = 2: THE 2D-TRANSL. VECTORS
!!$                               MUST HAVE THE FORM:
!!$                                      A1 = (P,Q)
!!$                                      A2 = (-P,Q)*U
!!$                               WHERE P AND Q ARE POSITIVE
!!$                               AND U=1 OR U=-1.

2042     ibz = 0
         tw = 0.0d0
         p = vbr (1, 1)
         q = vbr (2, 1)
         rnk = real (nk)
         ndx = nk
         ndy = nk
         If (p .Gt. (1.2d0*q)) ndy = Int (p*rnk/q)
         If (q .Gt. (1.2d0*p)) ndx = Int (q*rnk/p)
         dkx = pi / (p*real(ndx))
         dky = pi / (q*real(ndy))
         Do i = 1, ndx
            Do j = 1, ndy
               If ((2*j-1)*ndx .Gt. (2*i-1)*ndy) Cycle
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               wkbz (ibz) = 1.0d0
               If ((2*j-1)*ndx .Eq. (2*i-1)*ndy) wkbz (ibz) = 0.5d0
               tw = tw + wkbz (ibz)
               nbk1 (ibz) = i
               nbk2 (ibz) = j
               akbz (1, ibz) = (real(i)-0.5d0) * dkx
               akbz (2, ibz) = (real(j)-0.5d0) * dky
            End Do
         End Do

         Go To 290

!!$-------------------------- NSYM=10: GENERAL CASE PARTINIONED
!!$                                    FOR 8-fold parallelization
!!$ SHifT THE K|| MESH A LITTER BIT IN ORDER TO AVOID
!!$ THE KX=0 OR KY=0
888      ibz = 0
         tw = 0.0d0
         nk1 = nk
         nk2 = nk

         Do i1 = - nk1 + 1, 0
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - nk2 + 1, i1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = - nk1 + 1, - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = i1 + 1, 0
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1 - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - nk2 + 1, - i1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - i1 + 1, 0
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = - nk1 + 1, 0
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = 1, - i1 + 1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = - nk1 + 2, 0
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - i1 + 2, nk2
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = 1, i1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1 - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = i1 + 1, nk2
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Go To 290


!!$-------------------------- NSYM=11: GENERAL CASE PARTINIONED
!!$                                    FOR 8-fold parallelization
!!$ SHifT THE K|| MESH A LITTER BIT IN ORDER TO AVOID
!!$ THE KX=0 OR KY=0
988      ibz = 0
         tw = 0.0d0
         nk1 = nk
         nk2 = nk

         Do i1 = - nk1 + 1, 0
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - nk2 + 1, i1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = - nk1 + 1, - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = i1 + 1, 0
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1 - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - nk2 + 1, - i1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do


         Do i1 = nk1, 1, - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - i1 + 1, 0
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 0, - nk1 + 1, - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = 1, - i1 + 1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = - nk1 + 2, 0
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = - i1 + 2, nk2
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = 1, i1
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Do i1 = 1, nk1 - 1
            fx1 = (real(i1)-0.5d0) / real (2*nk1)
            Do i2 = i1 + 1, nk2
               fx2 = (real(i2)-0.5d0) / real (2*nk2)
               ibz = ibz + 1
!!$           if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 1.0d0
               tw = tw + wkbz (ibz)
               Do j = 1, 2
                  akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
               End Do
            End Do
         End Do

         Go To 290

!!$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$   Adaptive mesh, sort of. Only for 1x1 supercell - clean calculations

1088     Open (15, File='bzgrid', Status='old')
         Read (15,*) lf, bt, hg, wd ! read left,bottom,hight,width (like in Matlab) of the BZ part
         Close (15)
         ibz = 0
         tw = 1.0d0
         dx = (wd-lf) / nk
         dy = (hg-bt) / nk
         Do i1 = 1, nk
            fx1 = lf + (i1-0.5d0) * dx
            Do i2 = 1, nk
               fx2 = bt + (i2-0.5d0) * dy
               ibz = ibz + 1
!!$            if (ibz .gt. mnbz) goto 299
               nbk1 (ibz) = i1
               nbk2 (ibz) = i2
               wkbz (ibz) = 4 * dx * dy / ((2.0d0*pi)**2)
               akbz (1, ibz) = fx1
               akbz (2, ibz) = fx2
            End Do
         End Do

!!$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         Go To 290

290      nbz = ibz

         Do ibz = 1, nbz
            wkbz (ibz) = wkbz (ibz) / tw
         End Do


         Allocate (bz%k(2, nbz))
         Allocate (bz%ik(2, nbz))
         Allocate (bz%kweight(nbz))


         bz%nkgrid = nbz
         bz%k = akbz (:, 1:nbz)
         bz%kweight = wkbz (1:nbz)
         bz%ik (1, :) = nbk1 (1:nbz)
         bz%ik (2, :) = nbk2 (1:nbz)

         Write (IW6, 190) bz%nkgrid
         flush (IW6)
190      Format (/ / 3 X, '*****  NUMBER OF K||-POINTS:   NBZ=', I6)

         Return
      End Subroutine gibz
End Module bzgrid

