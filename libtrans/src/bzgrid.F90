!!$ $Id: bzgrid.F90 1974 2013-09-09 18:19:13Z antst $
#include "math_def.h"
Module bzgrid
!!$    This module contains procedures for generating k-mesh
!!$    for BZ integration copied from I.Turek LMTO code.
      Implicit None

      Type t_bz_opts
         Integer :: nsym=0, nk=10, inve=0, iftr=0
         Real (Kind=DEF_DBL_PREC) :: bzone (4), kpoint (2) = (/ 0.0d0, 0.0d0 /), BZpart=1.0d0
      End Type t_bz_opts

      Type bzone
         Integer :: nkgrid
         Integer, Pointer :: ik (:, :)
         Real (Kind=DEF_DBL_PREC) :: vbg(2,2),gomega
         Real (Kind=DEF_DBL_PREC), Pointer :: kweight (:), k (:, :)
      End Type bzone
Contains

  subroutine exceed_mnbz(ibz)
    Use logging
    integer :: ibz
    Character (Len=30) :: cwork
    Write (cwork, '(i6)') ibz
    Call do_log (1, '   **** ERROR IN GIBZ:  IBZ EXCEEDED MNBZ,  IBZ='//trim(cwork))
    Stop
  end subroutine exceed_mnbz

  
!!$*******************
!!$XXX    GIBZ    ****
!!$*******************
      Subroutine gibz (bz, vbr, bzopt)
        Use logging 
!!$****************************************************************
!!$  GENERATES NETWORK OF K||-POINTS IN IRREDUCIBLE BRILLOUIN ZONE
!!$****************************************************************
         Implicit None
         Type (t_bz_opts), Intent (In) :: bzopt
         Type (bzone) :: bz
         Real (Kind=DEF_DBL_PREC), Intent (In) :: vbr (2, 2)
!!$         Local vars
         Real (Kind=DEF_DBL_PREC) :: scale, bzdef (4), kpoint (2)
         Integer :: nk, nsym, inve, iftr
         Character (Len=30) :: cwork
         Integer :: nbk1 (2*(2*bzopt%nk+1)**2), nbk2 (2*(2*bzopt%nk+1)**2), istart !,  ext=0
         Real (Kind=DEF_DBL_PREC) :: akbz (2, 2*(2*bzopt%nk+1)**2), wkbz (2*(2*bzopt%nk+1)**2), vbg (2, 2)!, k_pos(2)
         Real (Kind=DEF_DBL_PREC) :: romega, gomega, ajm, rnk, anu3, dky, dkx, q, dk, lf, bt, hg, wd, dx, dy
         Real (Kind=DEF_DBL_PREC) :: p, fx2, fx1, tw, cnorm, det, stx, sty
         Integer :: nbz, iy, nu, ndy, ndx, nmy, nmx, i, j, i2, i1, nk1, nk2, ibz, mnbz


         Call do_log (1, ' Generating k mesh ...')


         scale = bzopt%BZpart
         nk = bzopt%nk
         nsym = bzopt%nsym
         inve = bzopt%inve
         iftr = bzopt%iftr
         bzdef = bzopt%bzone
         kpoint = bzopt%kpoint

         mnbz = (2*nk+1) ** 2

         det = vbr (1, 1) * vbr (2, 2) - vbr (1, 2) * vbr (2, 1)
         romega = Abs (det)
         gomega = (2.0d0*DEF_M_PI) ** 2 / romega
         cnorm = 2.0d0 * DEF_M_PI / det

         vbg (1, 1) = cnorm * vbr (2, 2)
         vbg (2, 1) = - cnorm * vbr (1, 2)
         vbg (1, 2) = - cnorm * vbr (2, 1)
         vbg (2, 2) = cnorm * vbr (1, 1)


         Write (cwork, '(g15.7)') romega
         Call do_log (1, '   AREA OF 2D-PRIMITIVE CELL: REAL= '//trim(cwork))
         Write (cwork, '(g15.7)') gomega
         Call do_log (1, '                        RECIPROCAL= '//trim(cwork))
         Write (cwork, '(2 g15.7)') vbg (:, 1)
         Call do_log (1, '   RECIP. BASIS: 1. VECTOR = '//trim(cwork))
         Write (cwork, '(2 g15.7)') vbg (:, 2)
         Call do_log (1, '                 2. VECTOR = '//trim(cwork))
         
         bz%vbg(:,:)=vbg
         bz%gomega=gomega        

!!$         Now we select from the specific in-plane symmetries of the system.
!!$         The cases are defined by the combination of BZ_SYM, BZ_INVE, and BZ_TIMEREVERSAL variables set in trans.conf.
!!$         Below, they are NSYM, INVE and IFTR respectively
         
         NSYM_SELECT: select case(nsym)
         case(0)           !!$-------------------------- NSYM=0: GENERAL CASE           
            ibz = 0         !!$ SHIFT THE K|| MESH A LITTER BIT IN ORDER TO AVOID        
            tw = 0.0d0      !!$ THE KX=0 OR KY=0                                         
            nk1 = nk
            nk2 = nk
            
            If (iftr == 0) Then
               istart = - nk1 + 1
            Else
               istart = 1
            End If
            stx = real ((2*nk1-1), kind=DEF_DBL_PREC)
            sty = real ((2*nk2-1), kind=DEF_DBL_PREC)
            Do i1 = istart, nk1-1
               fx1 = (real(i1, kind=DEF_DBL_PREC)) / stx
               Do i2 = - nk2 + 1, nk2-1
                  fx2 = (real(i2, kind=DEF_DBL_PREC)) / sty
                  ibz = ibz + 1
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                  nbk1 (ibz) = i1
                  nbk2 (ibz) = i2
                  wkbz (ibz) = 1.0d0
                  tw = tw + wkbz (ibz)
                  Do j = 1, 2
                     akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)!+pofs(j)
                  End Do
                  !                end if
               End Do
            End Do
            
         case(1)          !!$  NSYM = 1 - SYMMETRY OF CUBIC (001) SURFACES: 4-FOLD ROTATION, AXIS Z AND BY MIRROR PLANE X-Z            

            INVE_SELECT: select case(inve)
            case(1)                                                                               
               ibz = 0                       !!$---------- INVE=1                                                   
               tw = 0.0d0                    !!$ INVE = 1: THE 2D-TRANSL. VECTORS                                      
               p = vbr (1, 1)                !!$           MUST HAVE THE FORM:                                         
               dk = DEF_M_PI / (p*real(nk))  !!$                  A1 = (P,0)                                           
               Do i = 1, nk                  !!$                  A2 = (0,P)                                           
                  Do j = 1, i                !!$           WHERE P IS POSITIVE.                   
                     ibz = ibz + 1                                         
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                     wkbz (ibz) = 1.0d0
                     If (i .Eq. j) wkbz (ibz) = 0.5d0
                     tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = (real(i)-0.5d0) * dk
                     akbz (2, ibz) = (real(j)-0.5d0) * dk
                  End Do
               End Do
            case(2)                                    !!$ INVE = 2: THE 2D-TRANSL. VECTORS    
               ibz = 0                                 !!$           MUST HAVE THE FORM:      
               tw = 0.0d0                              !!$                  A1 = (P,P)        
               p = vbr (1, 1)                          !!$                  A2 = (-P,P)*U     
               dk = DEF_M_PI / (2.0d0*p*real(nk))      !!$           WHERE P IS POSITIVE      
               Do i = 1, nk                            !!$           AND U=1 OR U=-1.          
                  Do j = 1, i
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                     wkbz (ibz) = 1.0d0
                     If (i .Eq. j) wkbz (ibz) = 0.5d0
                     tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = real (i-j) * dk !+1.0D-4
                     akbz (2, ibz) = real (i+j-1) * dk
                  End Do
               End Do
            end select INVE_SELECT

         case(2)           !!$--  NSYM = 2 - SYMMETRY OF CUBIC (110) SURFACES: MIRROR PLANES X-Z AND Y-Z

            INVE_SELECT2: select case(inve) 
               
            case(1)               !!$  INVE = 1: THE 2D-TRANSL. VECTORS                     
               ibz = 0            !!$             MUST HAVE THE FORM:                       
               tw = 0.0d0         !!$                    A1 = (P,0)                         
               p = vbr (1, 1)     !!$                    A2 = (0,Q)                         
               q = vbr (2, 2)     !!$             WHERE P AND Q ARE POSITIVE.  
               nmx = 1            
               nmy = 1                           
               If (p .Lt. (0.75d0*q)) nmx = 2    
               If (q .Lt. (0.75d0*p)) nmy = 2
               If (p .Lt. (0.4d0*q)) Call do_log (0, ' **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ')
               If (q .Lt. (0.4d0*p)) Call do_log (0, ' **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ')
               ndx = nk * nmx
               ndy = nk * nmy
               dkx = DEF_M_PI / (p*real(ndx))
               dky = DEF_M_PI / (q*real(ndy))
               Do i = 1, ndx
                  Do j = 1, ndy
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                     wkbz (ibz) = 1.0d0
                      tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = (real(i)-0.5d0) * dkx
                     akbz (2, ibz) = (real(j)-0.5d0) * dky
                  End Do
               End Do

            case(2)
               ibz = 0                  !!$  INVE = 2: THE 2D-TRANSL. VECTORS        
               tw = 0.0d0               !!$            MUST HAVE THE FORM:           
               p = vbr (1, 1)           !!$                   A1 = (P,Q)             
               q = vbr (2, 1)           !!$                   A2 = (-P,Q)*U          
               nmx = 1                  !!$            WHERE P AND Q ARE POSITIVE    
               nmy = 1                  !!$            AND U=1 OR U=-1.              
               If (p .Lt. (0.75d0*q)) nmx = 2
               If (q .Lt. (0.75d0*p)) nmy = 2
               If (p .Lt. (0.4d0*q)) Call do_log (0, ' **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ')
               If (q .Lt. (0.4d0*p)) Call do_log (0, ' **** GIBZ - WARNING:   P DEVIATES TOO MUCH FROM Q ')
               ndx = nk * nmx
               ndy = nk * nmy
               dkx = DEF_M_PI / (p*real(ndx))
               dky = DEF_M_PI / (q*real(ndy))
               Do i = 1, ndx
                  Do j = 1, ndy
                     If ((2*j-1)*nmx .Gt. (2*i-1)*nmy) Cycle
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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

            end select INVE_SELECT2
            
         case(3)
            ibz = 0                                    !!$--------  NSYM = 3 - SYMMETRY OF CUBIC (111)  AND HEXAGONAL (0001) SURFACES:       
            tw = 0.0d0                                 !!$                     3-FOLD ROTATION, AXIS Z AND BY MIRROR PLANE Y-Z               
            p = vbr (1, 1)                             !!$          INVE - NOT ACTIVE.   THE 2D-TRANSL. VECTORS                              
            dkx = 2.0d0 * DEF_M_PI / (p*real(nk))      !!$                               MUST HAVE THE FORM:                                 
            dky = dkx / Sqrt (3.0d0)                   !!$                                      A1 = (P,0)                                   
            Do nu = 1, 2                               !!$                                   A2 = (P/2)*(1,S)*U                              
               anu3 = real (nu) / 3.0d0                !!$                               WHERE P IS POSITIVE,                                
               Do i = 1, nk                            !!$                               S=SQRT(3) OR S=-SQRT(3),                            
                  Do j = i, nk                         !!$                               AND U=1 OR U=-1.                                    
                     iy = i + 2 * j + nu - 3
                     If (iy .Gt. nk) Cycle
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
         case(4)                                          !!$--------  NSYM = 4 - SYMMETRY OF CUBIC (N10) SURFACES AND TILT GRAIN BOUNDARIES:  
                                                          !!$                     SYMMETRY GROUP GIVEN BY MIRROR PLANE X-Z
            INVE_SELECT3: select case(inve)
            case(1)
               ibz = 0                 !!$  INVE = 1: THE 2D-TRANSL. VECTORS     
               tw = 0.0d0              !!$            MUST HAVE THE FORM:        
               p = vbr (1, 1)          !!$                   A1 = (P,0)          
               q = vbr (2, 2)          !!$                   A2 = (0,Q)          
               rnk = real (nk)         !!$            WHERE P AND Q ARE POSITIVE.
               ndx = nk
               ndy = nk
               If (p .Gt. (1.2d0*q)) ndy = Int (p*rnk/q)
               If (q .Gt. (1.2d0*p)) ndx = Int (q*rnk/p)
               dkx = DEF_M_PI / (p*real(ndx))
               dky = DEF_M_PI / (q*real(ndy))
               Do i = 1, ndx
                  Do j = 1, ndy
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                     wkbz (ibz) = 1.0d0
                     tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = (real(i)-0.5d0) * dkx
                     akbz (2, ibz) = (real(j)-0.5d0) * dky
                  End Do
               End Do
               
            case(2)
               
               ibz = 0                 !!$   INVE = 2: THE 2D-TRANSL. VECTORS      
               tw = 0.0d0              !!$             MUST HAVE THE FORM:         
               p = vbr (1, 1)          !!$                    A1 = (P,Q)           
               q = vbr (2, 1)          !!$                    A2 = (-P,Q)*U        
               rnk = real (nk)         !!$             WHERE P AND Q ARE POSITIVE  
               ndx = nk                !!$             AND U=1 OR U=-1.            
               ndy = nk
               If (p .Gt. (1.2d0*q)) ndy = Int (p*rnk/q)
               If (q .Gt. (1.2d0*p)) ndx = Int (q*rnk/p)
               dkx = DEF_M_PI / (p*real(ndx))
               dky = DEF_M_PI / (q*real(ndy))
               Do i = 1, ndx
                  Do j = 1, ndy
                     If ((2*j-1)*ndx .Gt. (2*i-1)*ndy) Cycle
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                     wkbz (ibz) = 1.0d0
                     If ((2*j-1)*ndx .Eq. (2*i-1)*ndy) wkbz (ibz) = 0.5d0
                     tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = (real(i)-0.5d0) * dkx
                     akbz (2, ibz) = (real(j)-0.5d0) * dky
                  End Do
               End Do
            end select INVE_SELECT3

         case(10)         !!$------- NSYM=10: GENERAL CASE PARTINIONED FOR 8-fold parallelization
            ibz = 0
            tw = 0.0d0
            nk1 = nk
            nk2 = nk

            Do i1 = - nk1 + 1, 0
               fx1 = (real(i1)-0.5d0) / real (2*nk1)
               Do i2 = - nk2 + 1, i1
                  fx2 = (real(i2)-0.5d0) / real (2*nk2)
                  ibz = ibz + 1
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                  nbk1 (ibz) = i1
                  nbk2 (ibz) = i2
                  wkbz (ibz) = 1.0d0
                  tw = tw + wkbz (ibz)
                  Do j = 1, 2
                     akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
                  End Do
               End Do
            End Do
         case(11)                     !!$- NSYM=11: GENERAL CASE PARTINIONED  
                                      !!$           FOR 8-fold parallelization
            ibz = 0
            tw = 0.0d0
            nk1 = nk
            nk2 = nk

            Do i1 = - nk1 + 1, 0
               fx1 = (real(i1)-0.5d0) / real (2*nk1)
               Do i2 = - nk2 + 1, i1
                  fx2 = (real(i2)-0.5d0) / real (2*nk2)
                  ibz = ibz + 1
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
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
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                  nbk1 (ibz) = i1
                  nbk2 (ibz) = i2
                  wkbz (ibz) = 1.0d0
                  tw = tw + wkbz (ibz)
                  Do j = 1, 2
                     akbz (j, ibz) = fx1 * vbg (j, 1) + fx2 * vbg (j, 2)
                  End Do
               End Do
            End Do

         case(-1)                      !!$ NSYM=-1
                                       !!$    Sampling over a subsection of the BZ (1x1 only)
            lf = bzdef (1)
            wd = bzdef (2)
            bt = bzdef (3)
            hg = bzdef (4)
            
            ibz = 0
            tw = 1.0d0
            dx = (wd-lf) / nk
            dy = (hg-bt) / nk
            Do i1 = 1, nk
               fx1 = lf + (i1-0.5d0) * dx
               Do i2 = 1, nk
                  fx2 = bt + (i2-0.5d0) * dy
                  ibz = ibz + 1
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                  nbk1 (ibz) = i1
                  nbk2 (ibz) = i2
                  wkbz (ibz) = dx * dy / gomega
                  akbz (1, ibz) = fx1
                  akbz (2, ibz) = fx2
               End Do
            End Do
            
         case(-2)
            lf = bzdef (1)
            wd = bzdef (2)
            bt = bzdef (3)
            hg = bzdef (4)
            
            ibz = 0
            tw = 1.0d0
            dx = (wd-lf) / nk
            dy = (hg-bt) / nk
            Do i1 = 1, nk
               fx1 = lf + (i1-0.5d0) * dx
               Do i2 = 1, nk
                  fx2 = bt + (i2-0.5d0) * dy
                  ibz = ibz + 1
                  If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                  nbk1 (ibz) = i1
                  nbk2 (ibz) = i2
                  
                  wkbz (ibz) = dx * dy / gomega
                  akbz (1, ibz) = fx1
                  akbz (2, ibz) = fx2
               End Do
            End Do
            
            nbk1 ((ibz+1) :2*ibz) = - nbk1 (1:ibz)
            nbk2 ((ibz+1) :2*ibz) = - nbk2 (1:ibz)
            wkbz ((ibz+1) :2*ibz) = wkbz (1:ibz)
            akbz (1:2, (ibz+1) :2*ibz) = - akbz (1:2, 1:ibz)
            ibz = ibz * 2
            
         case(-3)           !!$--------  NSYM = 3 - SYMMETRY OF CUBIC (111)  AND HEXAGONAL (0001) SURFACES:  
                            !!$  This is a variant of NSYM=3 with -k|| points added.                         
            ibz = 0         !!$  Use when calculating T_mix (transmission mixing conductance).                                        
            tw = 0.0d0      !!$  Not necessary for regular G_mix (reflection mixing conductance).                                     
            p = vbr (1, 1)                           
            dkx = 2.0d0 * DEF_M_PI / (p*real(nk))    
            dky = dkx / Sqrt (3.0d0)
            Do nu = 1, 2
               anu3 = real (nu) / 3.0d0
               Do i = 1, nk
                  Do j = i, nk
                     iy = i + 2 * j + nu - 3
                     If (iy .Gt. nk) Cycle
                     ibz = ibz + 1
                     If (ibz .Gt. mnbz) call exceed_mnbz(ibz)
                     ajm = 1.0d0
                     If (i .Eq. j .Or. iy .Eq. nk) ajm = 2.0d0
                     If (i .Eq. j .And. iy .Eq. nk) ajm = 6.0d0
                     wkbz (ibz) = 1.0d0 / ajm
                     tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = (real(i-1)+anu3) * dkx
                     akbz (2, ibz) = (real(iy)) * dky
                     ibz = ibz + 1

                     wkbz (ibz) = wkbz (ibz-1)
                     tw = tw + wkbz (ibz)
                     nbk1 (ibz) = i
                     nbk2 (ibz) = j
                     akbz (1, ibz) = - akbz (1, ibz-1)
                     akbz (2, ibz) = - akbz (2, ibz-1)
                  End Do
               End Do
            End Do
            
         case(-10)    !!  Single k-point calculations

            ibz = 1
            akbz (1, 1) = kpoint (1)
            akbz (2, 1) = kpoint (2)
            nbk1 (1) = 1
            nbk2 (1) = 1
            wkbz (1) = 1.0d0
            tw = 1.0d0

         end select NSYM_SELECT

         nbz = ibz

         Do ibz = 1, nbz
            wkbz (ibz) = wkbz (ibz) / tw
         End Do


         Allocate (bz%k(2, nbz))
         Allocate (bz%ik(2, nbz))
         Allocate (bz%kweight(nbz))


         bz%nkgrid = nbz

         If (scale > 0.0d0 .And. scale /= 1.0d0) Then
            bz%k = akbz (:, 1:nbz) * Sqrt (scale)
            bz%kweight = wkbz (1:nbz) * scale
         Else
            bz%k = akbz (:, 1:nbz)
            bz%kweight = wkbz (1:nbz)
         End If

         bz%ik (1, :) = nbk1 (1:nbz)
         bz%ik (2, :) = nbk2 (1:nbz)

         Write (cwork, '(i8)') nbz
         
         Call do_log (1, '   NUMBER OF K||-POINTS:   NBZ='//trim(cwork))
         Call do_log (1, ' Done!')
         Call do_log (1, '')

      End Subroutine gibz
End Module bzgrid
