!!$ $Id: postprocess.F90 1664 2013-09-10 15:53:00Z rien $
#include "math_def.h"
Module interatcur
      Use sparselib
      Use geometry_module
      Use structure
      Use omta_kink
      Use omta_SOC
      Implicit None

      Type t_curtensor
         Integer :: nnz = 0
         Real (Kind=DEF_DBL_PREC), Pointer :: a (:, :, :) ! nnz x 4 x 2 array
         Integer, Pointer :: ir (:), jc (:)
      End Type
      
      Type t_tr_components
         Integer :: n
         !Integer, Pointer :: trvia(:, :)
         Integer, Allocatable :: trvia(:, :)
      End Type
      
      Type t_interat_currents
         Integer :: ntr = 0, nat = 0, nch = 0, alloc = 0
         Integer, Pointer :: trlist(:,:) ! 2 x ntr matrix containing all relevant translations of the uc
         Integer, Pointer :: direct(:)
         Type (t_curtensor), Pointer :: curten(:)
         Type (t_tr_components), Pointer :: via(:)
         Real (Kind=DEF_DBL_PREC), Pointer :: cons(:, :)
      End Type
      
Contains
     
        Subroutine dealloc_iatcu(iatcu)
           Implicit None
           Type(t_interat_currents) :: iatcu
           Integer i
           
           If (iatcu%alloc /= 0) Then
              Deallocate(iatcu%trlist)
              Deallocate(iatcu%direct)
              Do i = 1, iatcu%ntr
                 Call dealloc_curten(iatcu%curten(i))
              End Do
              Deallocate(iatcu%curten)
              Do i = 1, iatcu%ntr
                 If (Allocated(iatcu%via(i)%trvia)) Call dealloc_via(iatcu%via(i))
              End Do
              Deallocate(iatcu%via)
              Deallocate(iatcu%cons)
              iatcu%ntr = 0
              iatcu%nat = 0
              iatcu%alloc = 0
           End If
        End Subroutine dealloc_iatcu
     
        Subroutine alloc_iatcu(iatcu, ntr, nat, nch)
           Implicit None
           Type(t_interat_currents) :: iatcu
           Integer :: ntr, nat, nch ! nch is the number of channels, i.e. 
                                    !                       2 for spin (up dn)
                                    !                       6 for spin+orbital (up dn)*(lx ly lz) 
           
           Call Dealloc_iatcu(iatcu)
           Allocate(iatcu%trlist(2, ntr))
           Allocate(iatcu%direct(ntr))
           Allocate(iatcu%curten(ntr))
           Allocate(iatcu%via(ntr))
           Allocate(iatcu%cons(2, nat))
           iatcu%cons = 0.0d0
           iatcu%ntr = ntr
           iatcu%nat = nat
           iatcu%nch = nch
           iatcu%alloc = 1
        End Subroutine alloc_iatcu
        
        Subroutine make_empty_iatcu(iatcu, ntr)
           Implicit None
           Type(t_interat_currents) :: iatcu
           Integer :: ntr
           
           Integer :: itr
           Call alloc_iatcu(iatcu, ntr, 1, 1)
           Do itr = 1, ntr
              Call alloc_curten(iatcu%curten(itr), 1, 1, 1)
           End Do
        End Subroutine
     
        Function copy_iatcu(iatcu) Result(copy)
           Implicit None
           Type(t_interat_currents) :: iatcu, copy
           Integer :: itr, ntr, nnz
           ntr = iatcu%ntr
           Call alloc_iatcu(copy, ntr, iatcu%nat, iatcu%nch)
           copy%trlist = iatcu%trlist
           copy%direct = iatcu%direct
           copy%cons = iatcu%cons
           Do itr = 1, ntr
              Call alloc_via(copy%via(itr), iatcu%via(itr)%n)
              copy%via(itr)%trvia = iatcu%via(itr)%trvia
              nnz = iatcu%curten(itr)%nnz
              Call alloc_curten(copy%curten(itr), nnz, copy%nat, copy%nch)
              copy%curten(itr)%ir = iatcu%curten(itr)%ir
              copy%curten(itr)%jc = iatcu%curten(itr)%jc
           End Do
        End Function copy_iatcu
        
        Subroutine zero_iatcu(iatcu)
           Implicit None
           Type (t_interat_currents), Intent (Inout) :: iatcu
              
           Integer :: itr
           Do itr = 1, iatcu%ntr
              iatcu%curten(itr)%a = 0.0d0
           End Do
           iatcu%cons(:, :) = 0.0d0
        End Subroutine
    
  !!$ EMTO: calc_iatcu_EMTO

  Subroutine calc_iatcu_EMTO(rhs, msc, msdotc, energy, par_EMTO, trgeo, rneibs, opt, iatcu, k, nl, nr, do_iatcu, do_iatLcu, comm)
           Use hamiltonian
           Use readcfg
           Implicit None
           Type (zdensemat), Intent (In) :: rhs
           Type (t_strconst), Intent(In) :: msc
           Type (t_strconst), Intent(In) :: msdotc
           Real (Kind=DEF_DBL_PREC) :: energy
           Type (t_omta_logder), Intent (In) :: par_EMTO(2)
           Type (t_geometry_EMTO), Intent (In) :: trgeo
           Type (t_remember_neibs), Intent(In) :: rneibs(:)
           Type (t_options) :: opt
           Type (t_interat_currents), Intent (Inout) :: iatcu ! interatomic currents
           Real (Kind=DEF_DBL_PREC), Intent(In) :: k (2)
           Integer, Intent (In) :: nl, nr, do_iatcu, do_iatLcu
           Integer :: comm

           Integer :: itr, tr(2), ivia, tr1(2), tr2(2), flag, process_rank, ierror
           Type (zcsrmat) :: ham, ham_dot !, H1, dum, H2
           Type (t_mathop) :: sk
           Type (t_mathop) :: sdotk
           !mpi Integer :: alloc

           !mpi Call MPI_Comm_rank(MPI_COMM_WORLD, process_rank, ierror) 
           ! the iatcu struct only exists on the root mpi rank
           ! so we need to communicate often between ranks

           !mpi If (process_rank == 0) then
           Call zero_iatcu(iatcu)
           !mpi alloc = iatcu%alloc
           !mpi Endif
           !mpi call MPI_Bcast(alloc, 1, MPI_INT, 0, comm, ierror) !root rank broadcasts alloc status to all ranks
           !If (alloc == 0) Return

           Do itr = 1, iatcu%ntr

            !mpi If (process_rank == 0) tr = iatcu%trlist(:,itr)
            tr = iatcu%trlist(:,itr)
            !mpi Call MPI_Bcast(tr, 1, MPI_INT, 0, comm, ierror)

            !mpi If (process_rank == 0) then
             If (iatcu%direct(itr) /= 0) Then
              flag = setflag(tr)
              Call make_sk_decomp (sk, msc, k, tr, 0)
              ham = sk%c
              Call free_mathop(sk)
              If ( flag /= 0 ) Call addtoM (ham, par_EMTO(1)%logder)
              If ( flag /= 0 ) Call addtoM (ham, par_EMTO(1)%hsoc)
              ham%a = -1*ham%a !attempt
             Else
              Call alloc (ham, 0, par_EMTO(1)%logder%nrow)
             End If
            !mpi endif

            !!TODO WIP parallelization of the currents calculations
            ! we probably require the parallelization of both the construction of the Hamiltonian as well
            ! as the actual computations of the currents... It should be done on a 'row-atom' basis, as the
            ! Hamiltonian is constructed that way AND the currents are calculated that way.
            Call calc_curr_EMTO (rhs, ham, trgeo, rneibs, itr, iatcu, nl, nr, do_iatcu, do_iatLcu, comm)

            !mpi If (process_rank == 0) Call free (ham)
            Call free (ham)
           End Do
           !mpi If (process_rank == 0) Call free_mathop(sk)
           Call free_mathop(sk)
        
      Contains

         Function setflag(tr) Result(flag)
            Implicit None
            Integer :: tr(2), flag
            If ( tr(1) == 0 .and. tr(2) == 0 ) Then
               flag = 1
            Else
               flag = 0
            End If
         End Function
      End Subroutine calc_iatcu_EMTO
!!$ end EMTO
 
Subroutine calc_iatcu(rhs, msc, mpar, trgeo, rneibs, opt, iatcu, k, nl, nr, do_iatcu, do_iatLcu, comm)
         Use hamiltonian
         Use readcfg
         Implicit None
         Type (zdensemat), Intent (In) :: rhs
         Type (t_strconst), Intent(In) :: msc
         Type (t_potpar_mats), Intent(In) :: mpar
         Type (t_geometry), Intent (In) :: trgeo
         Type (t_remember_neibs), Intent(In) :: rneibs(:)
         Type (t_options) :: opt
         Type (t_interat_currents), Intent (Inout) :: iatcu ! interatomic currents
         Real (Kind=DEF_DBL_PREC), Intent(In) :: k (2)
         Integer, Intent (In) :: nl, nr, do_iatcu, do_iatLcu
         Integer :: comm
         
         Type (zcsrmat) :: ham, H1, dum, H2
         Type (t_mathop) :: sk
         Integer :: itr, tr(2), ivia, tr1(2), tr2(2), flag
         
         Call zero_iatcu(iatcu)

         If (opt%po%kind > 3) Return
         If (iatcu%alloc == 0) Return

         Call init_mathop(sk)
         Do itr = 1, iatcu%ntr
            tr = iatcu%trlist(:,itr)
            If (iatcu%direct(itr) /= 0) Then
               flag = setflag(tr)
               Call make_sk_decomp (sk, msc, k, tr, 0)
               H1 = getH(sk%c, mpar, flag)
               Call free_mathop(sk)
               ham = make1cent (H1, mpar%c1, mpar%c1conj)
               Call free (H1)
               If ( flag /= 0 ) Call addtoM (ham, mpar%c0)
            Else
               Call alloc (ham, 0, mpar%c0%nrow)
            End If
            
            If ( iatcu%via(itr)%n > 0 ) Then
               Do ivia = 1, iatcu%via(itr)%n
                  tr1 = iatcu%via(itr)%trvia(:, ivia)
                  tr2 = tr - tr1
                  Call make_sk_decomp (sk, msc, k, tr1, 0)
                  flag = setflag(tr1)
                  H1 = getH (sk%c, mpar, flag)
                  Call free_mathop(sk)
                  Call make_sk_decomp (sk, msc, k, tr2, 0)
                  flag = setflag(tr2)
                  H2 = getH (sk%c, mpar, flag)
                  Call free_mathop(sk)
                  dum = get3MMM(H2, mpar%c2, H1)
                  Call addtoM(ham, dum)
                  Call free(dum)
                  Call free(H1)
                  Call free(H2)
               End Do
            End If


            Call calc_curr (rhs, ham, trgeo, rneibs, itr, iatcu, nl, nr, do_iatcu, do_iatLcu, comm)
            call free (ham)
         End Do
         Call free_mathop(sk)
         
      Contains
         
         Function setflag(tr) Result(flag)
            Implicit None
            Integer :: tr(2), flag
            If ( tr(1) == 0 .and. tr(2) == 0 ) Then
               flag = 1
            Else
               flag = 0
            End If
         End Function
      End Subroutine calc_iatcu
      
      Subroutine addtoM(mtot, madd)
         Implicit None
         Type (zcsrmat), Intent(Inout) :: mtot
         Type (zcsrmat), Intent(In) :: madd
         Type (zcsrmat) :: dum
         dum = spmatadd (mtot, madd)
         Call free(mtot)
         Call spcopy(dum, mtot)
         Call free(dum)
      End Subroutine addtoM
      
      Subroutine bzsum_iatcu(iatcu, iatcu1, kw)
         Type (t_interat_currents), Intent (Inout) :: iatcu
         Type (t_interat_currents), Intent (In) :: iatcu1
         Real (Kind=DEF_DBL_PREC), Intent(In) :: kw
         
         Integer :: itr
         If (iatcu%alloc == 0) Return
         Do itr = 1, iatcu%ntr
            iatcu%curten(itr)%a(:, :, :) = iatcu%curten(itr)%a(:, :, :) + iatcu1%curten(itr)%a(:, :, :)*kw
         End Do
         iatcu%cons(:, :) = iatcu%cons(:, :) + iatcu1%cons(:, :)*kw
      End Subroutine bzsum_iatcu
      
      Subroutine alloc_via(via, n)
         Implicit None
         Type (t_tr_components) :: via
         Integer :: n
         Allocate(via%trvia(2, n))
         via%n = n
      End Subroutine alloc_via
      
      Subroutine dealloc_via(via)
         Implicit None
         Type (t_tr_components) :: via
         Deallocate(via%trvia)
         via%n = 0
      End Subroutine
      
      Subroutine alloc_curten(curten, nnz, nat, nch)
         Implicit None
         Type (t_curtensor) :: curten
         Integer :: nnz, nat, nch !nch is the number of current channels
                                  ! sx,sy,sz,c makes 4 (interatomic (spin) currents)
                                  ! (up dn)*(lx ly lz) makes 6 (interatomic orbital currents split into up and dn)
         
         If (nnz == 0) Then
            Allocate(curten%ir(0))
         Else
            Allocate(curten%ir(nat+1))
         End If
         Allocate(curten%jc(nnz))
         Allocate(curten%a(nnz, nch, 2))
         curten%nnz = nnz
         curten%a = 0.0d0
         curten%ir = 1
         curten%jc = 0
      End Subroutine alloc_curten
      
      Subroutine dealloc_curten(curten)
         Implicit None
         Type (t_curtensor) :: curten
         
         Deallocate(curten%ir)
         Deallocate(curten%jc)
         Deallocate(curten%a)
         curten%nnz = 0
      End Subroutine dealloc_curten
      
      Subroutine init_iatcu (rneibs, iatcu, trgeo, kind, nch)
         Implicit None
         Type (t_interat_currents), Intent (Inout), target :: iatcu
         Type (t_remember_neibs), Intent(inout) :: rneibs(:)
         Integer, Intent(in) :: kind, nch
         Type (t_geometry), Intent (In) :: trgeo
         
         Integer :: ntrmax
         Integer :: ia, ia2, inb, inb2, nat, ntr, n, i, itr
         Integer, allocatable :: list(:, :), ind(:)
         Integer :: tr1(2), tr2(2), ione(1), vec(3)
         Integer :: maxnn, nnz
         Type (t_curtensor), Pointer :: ctp
         Type (t_remember_neibs), Allocatable :: rneibs2(:)
         Integer, Allocatable :: dumvia(:,:,:)
         Integer, Allocatable :: nvia(:), dumdirect(:)
         
         ntrmax = (trgeo%ntrpar*4+1)**2 ! (the *4 is to be safe when kind == 3)
         nat = size(rneibs)
         iatcu%nat = nat
         maxnn = 0
         Do ia = 1, nat
            If (rneibs(ia)%n > maxnn) maxnn = rneibs(ia)%n
         End Do
         Allocate(list(2, ntrmax))
         Allocate(dumvia(2,maxnn,ntrmax))
         Allocate(nvia(ntrmax), dumdirect(ntrmax))
         
         ! Compose a list of all relevant translations
         ntr = 0
         nvia = 0
         dumvia = -1
         dumdirect = 0
         Do ia = 1, nat
            Do inb = 1, rneibs(ia)%n
               tr1 = rneibs(ia)%tr(:,inb)
               Call check_list(tr1, list, ntr, itr)
               dumdirect(itr) = 1
               If (kind == 3) Then
                  ia2 = rneibs(ia)%ai(inb)
                  Do inb2 = 1, rneibs(ia2)%n
                     tr2 = tr1 + rneibs(ia2)%tr(:,inb2)
                     Call check_list(tr2, list, ntr, itr)
                     Call check_list(tr1, dumvia(1:2,:,itr), nvia(itr)) ! needed for next thing
                  End Do
               End If
            End Do
         End Do
         Call alloc_iatcu(iatcu, ntr, nat, nch)
         iatcu%trlist(:,1:ntr) = list(:,1:ntr)
         
         ! Per translation, compose a list of possible combinations
         ! of translations that lead to that translation
         Do itr = 1, ntr
            Call alloc_via(iatcu%via(itr), nvia(itr))
            iatcu%via(itr)%trvia = dumvia(:, 1:nvia(itr), itr)
            iatcu%direct(itr) = dumdirect(itr)
         End Do
         Deallocate(list, dumvia, dumdirect)
         
         
         ! If kind = 3 there are more hoppings than stored in rneibs
         ! So get a new rneibs that includes 3 center terms
         If (kind == 3) Then
            Allocate(list(3, maxnn**2))
            Allocate(rneibs2(nat))
            Do ia = 1, nat
               n = rneibs(ia)%n
               list(1:2,1:n) = rneibs(ia)%tr
               list(3,1:n) = rneibs(ia)%ai
               Do inb = 1, rneibs(ia)%n
                  tr1 = rneibs(ia)%tr(:,inb)
                  ia2 = rneibs(ia)%ai(inb)
                  Do inb2 = 1, rneibs(ia2)%n
                     vec(1:2) = tr1 + rneibs(ia2)%tr(:,inb2)
                     vec(3) = rneibs(ia2)%ai(inb2)
                     Call check_list(vec, list, n)
                  End Do
               End Do
               Call alloc_rneib(rneibs2(ia), n)
               rneibs2(ia)%tr = list(1:2,1:n)
               rneibs2(ia)%ai = list(3,1:n)
               rneibs2(ia)%ofs = rneibs(ia)%ofs
            End Do
            Call copy_rnb(rneibs2, rneibs)
            Deallocate(rneibs2)
            Deallocate(list)
            Do ia = 1, nat
               If (rneibs(ia)%n > maxnn) maxnn = rneibs(ia)%n
            End Do
         End If
         
         ! We don't want to keep currents between the
         ! lead atoms included in the scattering region
         Allocate(list(3, maxnn))
         Allocate(rneibs2(nat))
         Do ia = 1, nat
            If (ia > trgeo%ltna .and. ia <= nat - trgeo%rtna) Then
               Call copy_rnb(rneibs(ia:ia), rneibs2(ia:ia))
            Else
               n = 0
               Do inb = 1, rneibs(ia)%n
                  ia2 = rneibs(ia)%ai(inb)
                  If (ia2 > trgeo%ltna .and. ia2 <= nat - trgeo%rtna) Then
                     n = n + 1
                     list(3,n) = ia2
                     list(1:2,n) = rneibs(ia)%tr(:, inb)
                  End If
               End Do
               Call alloc_rneib(rneibs2(ia), n)
               rneibs2(ia)%ai = list(3,1:n)
               rneibs2(ia)%tr = list(1:2,1:n)
               rneibs2(ia)%ofs = rneibs(ia)%ofs
            End If
         End Do
         Call copy_rnb(rneibs2, rneibs)
         Deallocate(list)
         Deallocate(rneibs2)
         
         ! construct the sparse tensors
         Allocate(list(1, maxnn))
         Allocate(ind(maxnn))
         Do itr = 1, ntr
            ctp => iatcu%curten(itr)
            nnz = 0
            Do ia = 1, nat ! first loop is to determine sizes of curten
               n = 0
               Do inb = 1, rneibs(ia)%n
                  tr1 = rneibs(ia)%tr(:,inb)
                  If ( all(tr1 .eq. iatcu%trlist(:,itr)) ) Then
                     ione(1) = rneibs(ia)%ai(inb)
                     Call check_list(ione, list, n)
                  End If
               End Do
               nnz = nnz + n
            End Do
            Call alloc_curten(ctp, nnz, nat, nch) ! allocate curten
            If (nnz > 0) Then
               Do ia = 1, nat ! second loop is to fill row and col indices ( ir and jc)
                  n = 0
                  Do inb = 1, rneibs(ia)%n
                     tr1 = rneibs(ia)%tr(:,inb)
                     If ( all(tr1 .eq. iatcu%trlist(:,itr)) ) Then
                        ione(1) = rneibs(ia)%ai(inb)
                        Call check_list(ione, list, n)
                     End If
                  End Do
                  ctp%ir(ia+1) = ctp%ir(ia) + n
                  Call iqsort(n, list(1,1:n), ind(1:n)) ! sorting routine from sparselib
                  Do i = 1, n
                     ctp%jc(ctp%ir(ia)+i-1) = list(1,ind(i))
                  End Do
               End Do
            End If
         End Do
         Deallocate(list, ind)
      End Subroutine init_iatcu
     
      Subroutine init_iatcu_EMTO (rneibs, iatcu, trgeo, kind, nch)
         Implicit None
         Type (t_interat_currents), Intent (Inout), target :: iatcu
         Type (t_remember_neibs), Intent(inout) :: rneibs(:)
         Integer, Intent(in) :: kind, nch
         Type (t_geometry_EMTO), Intent (In) :: trgeo
         
         Integer :: ntrmax
         Integer :: ia, ia2, inb, inb2, nat, ntr, n, i, itr
         Integer, allocatable :: list(:, :), ind(:)
         Integer :: tr1(2), tr2(2), ione(1), vec(3)
         Integer :: maxnn, nnz
         Type (t_curtensor), Pointer :: ctp
         Type (t_remember_neibs), Allocatable :: rneibs2(:)
         Integer, Allocatable :: dumvia(:,:,:)
         Integer, Allocatable :: nvia(:), dumdirect(:)
         
         ntrmax = (trgeo%ntrpar*4+1)**2 ! (the *4 is to be safe when kind == 3)
         nat = size(rneibs)
         iatcu%nat = nat
         maxnn = 0
         Do ia = 1, nat
            If (rneibs(ia)%n > maxnn) maxnn = rneibs(ia)%n
         End Do
         Allocate(list(2, ntrmax))
         Allocate(dumvia(2,maxnn,ntrmax))
         Allocate(nvia(ntrmax), dumdirect(ntrmax))
         
         ! Compose a list of all relevant translations
         ntr = 0
         nvia = 0
         dumvia = -1
         dumdirect = 0
         Do ia = 1, nat
            Do inb = 1, rneibs(ia)%n
               tr1 = rneibs(ia)%tr(:,inb)
               Call check_list(tr1, list, ntr, itr)
               dumdirect(itr) = 1
               If (kind == 3) Then
                  ia2 = rneibs(ia)%ai(inb)
                  Do inb2 = 1, rneibs(ia2)%n
                     tr2 = tr1 + rneibs(ia2)%tr(:,inb2)
                     Call check_list(tr2, list, ntr, itr)
                     Call check_list(tr1, dumvia(1:2,:,itr), nvia(itr)) ! needed for next thing
                  End Do
               End If
            End Do
         End Do
         Call alloc_iatcu(iatcu, ntr, nat, nch)
         iatcu%trlist(:,1:ntr) = list(:,1:ntr)
         
         ! Per translation, compose a list of possible combinations
         ! of translations that lead to that translation
         Do itr = 1, ntr
            Call alloc_via(iatcu%via(itr), nvia(itr))
            iatcu%via(itr)%trvia = dumvia(:, 1:nvia(itr), itr)
            iatcu%direct(itr) = dumdirect(itr)
         End Do
         Deallocate(list, dumvia, dumdirect)
         
         
         ! If kind = 3 there are more hoppings than stored in rneibs
         ! So get a new rneibs that includes 3 center terms
!         If (kind == 3) Then
!            Allocate(list(3, maxnn**2))
!            Allocate(rneibs2(nat))
!            Do ia = 1, nat
!               n = rneibs(ia)%n
!               list(1:2,1:n) = rneibs(ia)%tr
!               list(3,1:n) = rneibs(ia)%ai
!               Do inb = 1, rneibs(ia)%n
!                  tr1 = rneibs(ia)%tr(:,inb)
!                  ia2 = rneibs(ia)%ai(inb)
!                  Do inb2 = 1, rneibs(ia2)%n
!                     vec(1:2) = tr1 + rneibs(ia2)%tr(:,inb2)
!                     vec(3) = rneibs(ia2)%ai(inb2)
!                     Call check_list(vec, list, n)
!                  End Do
!               End Do
!               Call alloc_rneib(rneibs2(ia), n)
!               rneibs2(ia)%tr = list(1:2,1:n)
!               rneibs2(ia)%ai = list(3,1:n)
!               rneibs2(ia)%ofs = rneibs(ia)%ofs
!            End Do
!            Call copy_rnb(rneibs2, rneibs)
!            Deallocate(rneibs2)
!            Deallocate(list)
!            Do ia = 1, nat
!               If (rneibs(ia)%n > maxnn) maxnn = rneibs(ia)%n
!            End Do
!         End If
         
         ! We don't want to keep currents between the
         ! lead atoms included in the scattering region
         Allocate(list(3, maxnn))
         Allocate(rneibs2(nat))
         Do ia = 1, nat
            If (ia > trgeo%ltna .and. ia <= nat - trgeo%rtna) Then
               Call copy_rnb(rneibs(ia:ia), rneibs2(ia:ia))
            Else
               n = 0
               Do inb = 1, rneibs(ia)%n
                  ia2 = rneibs(ia)%ai(inb)
                  If (ia2 > trgeo%ltna .and. ia2 <= nat - trgeo%rtna) Then
                     n = n + 1
                     list(3,n) = ia2
                     list(1:2,n) = rneibs(ia)%tr(:, inb)
                  End If
               End Do
               Call alloc_rneib(rneibs2(ia), n)
               rneibs2(ia)%ai = list(3,1:n)
               rneibs2(ia)%tr = list(1:2,1:n)
               rneibs2(ia)%ofs = rneibs(ia)%ofs
            End If
         End Do
         Call copy_rnb(rneibs2, rneibs)
         Deallocate(list)
         Deallocate(rneibs2)
         
         ! construct the sparse tensors
         Allocate(list(1, maxnn))
         Allocate(ind(maxnn))
         Do itr = 1, ntr
            ctp => iatcu%curten(itr)
            nnz = 0
            Do ia = 1, nat ! first loop is to determine sizes of curten
               n = 0
               Do inb = 1, rneibs(ia)%n
                  tr1 = rneibs(ia)%tr(:,inb)
                  If ( all(tr1 .eq. iatcu%trlist(:,itr)) ) Then
                     ione(1) = rneibs(ia)%ai(inb)
                     Call check_list(ione, list, n)
                  End If
               End Do
               nnz = nnz + n
            End Do
            Call alloc_curten(ctp, nnz, nat, nch) ! allocate curten
            If (nnz > 0) Then
               Do ia = 1, nat ! second loop is to fill row and col indices ( ir and jc)
                  n = 0
                  Do inb = 1, rneibs(ia)%n
                     tr1 = rneibs(ia)%tr(:,inb)
                     If ( all(tr1 .eq. iatcu%trlist(:,itr)) ) Then
                        ione(1) = rneibs(ia)%ai(inb)
                        Call check_list(ione, list, n)
                     End If
                  End Do
                  ctp%ir(ia+1) = ctp%ir(ia) + n
                  Call iqsort(n, list(1,1:n), ind(1:n)) ! sorting routine from sparselib
                  Do i = 1, n
                     ctp%jc(ctp%ir(ia)+i-1) = list(1,ind(i))
                  End Do
               End Do
            End If
         End Do
         Deallocate(list, ind)
      End Subroutine init_iatcu_EMTO

 
      Subroutine copy_rnb(rnb, copy)
         Implicit None
         Type (t_remember_neibs) :: rnb(:)
         Type (t_remember_neibs) :: copy(:)
         
         Integer :: nat, ia, inb
         nat = size(rnb)
         Do ia = 1, nat
            Call alloc_rneib(copy(ia), rnb(ia)%n)
            Do inb = 1, rnb(ia)%n
               copy(ia)%ai(inb) = rnb(ia)%ai(inb)
               copy(ia)%tr(:,inb) = rnb(ia)%tr(:,inb)
            End Do
            copy(ia)%ofs = rnb(ia)%ofs
         End Do
      End Subroutine copy_rnb
      
      Subroutine check_list(vec_in, list, n, i_out)
         Implicit None
         Integer :: vec_in(:)
         Integer :: list(:,:), n
         Integer, optional :: i_out
         
         Integer :: i
         Logical :: exists
         
         exists = .false.
         Do i = 1, n
            If ( all( vec_in .eq. list(:,i) ) ) Then
               exists = .true.
               If (present(i_out)) i_out = i
               Exit
            End If
         End Do
         
         If (exists .neqv. .true.) Then
            n = n + 1
            list(:, n) = vec_in
            If (present(i_out)) i_out = n
         End If
      End Subroutine check_list

!!$ EMTO: calc_curr_EMTO

Subroutine calc_curr_EMTO (rhs, ham, trgeo, rneibs, itr, iatcu, nl, nr, do_iatcu, do_iatLcu, comm)
         Implicit None
         Type (t_geometry_EMTO), Intent (In), Target :: trgeo
         !Type (zdensemat), Intent (In), Target :: rhs
         Type (zdensemat), Intent (In) :: rhs
         Type (zcsrmat), Intent (In) :: ham
         Type (t_interat_currents), Intent (Inout) :: iatcu ! interatomic currents
         Integer, Intent (In) :: nl, nr, itr, do_iatcu, do_iatLcu
         Integer :: comm
         Type (t_remember_neibs) :: rneibs(:)
!!$           Local vars         
         Complex (Kind=DEF_DBL_PREC) :: hamp (2*trgeo%nlmax**2, 2*trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: psi_r(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: psi_r_swap(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: psi_c(2*trgeo%nlmax**2)
         !Complex (Kind=DEF_DBL_PREC), Pointer :: psi_r(:)
         !Complex (Kind=DEF_DBL_PREC), Pointer :: psi_r_swap(:)
         !Complex (Kind=DEF_DBL_PREC), Pointer :: psi_c(:)
         Complex (Kind=DEF_DBL_PREC) :: prod(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod_lx(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod_ly(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod_lz(trgeo%nlmax**2)
         Type (tclass), Pointer :: pr, pc
         Real (Kind=DEF_DBL_PREC) :: rr(3)
         Integer :: ofsr, ofsr1, ofsc, nor, noc, imd, ld, iar, iac, isr, i1, i2, j, ind
         Integer :: spinsign, swap, tr(2)
         Real (Kind=DEF_DBL_PREC) :: jval(2)
         Real (Kind=DEF_DBL_PREC) :: lxval(2), lyval(2), lzval(2)
         Complex (Kind=DEF_DBL_PREC) :: jxy(2)
         Type(t_curtensor), Pointer :: ctp
         Integer :: numls(5) = [1,3,5,7,9]
         ! for Lcurs
         Type(t_Lop_mats) :: Lop(trgeo%nlmax)

         if (do_iatLcu /= 0) call make_bigLmat(trgeo%nlmax, Lop)

         tr = iatcu%trlist(:, itr)
         ctp => iatcu%curten(itr)
         ofsr1 = 0
        
         Do iar = 1, trgeo%num ! loop over 'row atoms' of hamiltonian
            pr => trgeo%atoms(iar)%ptr
            rr = trgeo%atoms(iar)%coord(:)
            nor = (pr%lmx+1)**2 - dot_product(pr%idxdn-1,numls(1:pr%lmx+1))
       
            Do ind = ctp%ir(iar), ctp%ir(iar+1)-1
               iac = ctp%jc(ind)
               pc => trgeo%atoms(iac)%ptr
               ofsc = rneibs(iac)%ofs
               noc = 2 * ((pc%lmx+1)**2 - dot_product(pc%idxdn-1,numls(1:pr%lmx+1)))
               ! get the relevant part of the hamiltonian
               hamp = 0.0d0

               Do i1 = 1, 2*nor
                  Do j = ham%ir(ofsr1+i1), ham%ir(ofsr1+i1+1)-1
                       If (ham%jc(j) > ofsc ) Then
                        i2 = ham%jc(j)-ofsc
                        If ( i2 > noc ) Exit
                        !hamp(i2,i1) = ham%a(j) ! (transposed)
                        hamp(i1,i2) = ham%a(j) 
                       End If
                  End Do
               End Do

               Do isr = 0, 1 ! row spin (0 = up, 1 = down)
                  ofsr = ofsr1 + isr*nor
                  spinsign = sign(1,-isr)
                  swap = spinsign*nor
                  jval = 0.0d0
                  jxy = DEF_cmplx_zero
                  lxval = 0.0d0
                  lyval = 0.0d0
                  lzval = 0.0d0
                  ! Do the rest
                  Do imd = 1, nl + nr
                     psi_r(1:nor) = rhs%bl(ofsr+1:ofsr+nor, imd)
                     psi_r_swap(1:nor) = rhs%bl(ofsr+1+swap:ofsr+nor+swap, imd)
                     psi_c(1:noc) = rhs%bl(ofsc+1:ofsc+noc, imd)

                     !Calculate H|iac>
                     if (do_iatcu /= 0)  then
                      !prod(1:nor) = matmul ( hamp(1+isr*nor:nor+isr*nor, 1:noc), psi_c(1:noc) )
                      prod(1:nor) = matmul ( hamp(1:noc,1+isr*nor:nor+isr*nor), psi_c )
                     endif
                     if (do_iatLcu /= 0) then
                      prod_lx(1:nor) = matmul( hamp(1:noc, 1+isr*nor:nor+isr*nor), matmul(Lop(pc%lmx+1)%Lx, psi_c ) )
                      prod_ly(1:nor) = matmul( hamp(1:noc, 1+isr*nor:nor+isr*nor), matmul(Lop(pc%lmx+1)%Ly, psi_c ) )
                      prod_lz(1:nor) = matmul( hamp(1:noc, 1+isr*nor:nor+isr*nor), matmul(Lop(pc%lmx+1)%Lz, psi_c ) )
                     endif

                     If (imd > nl) Then
                        ld = 2
                     Else
                        ld = 1
                     End If

                     if (do_iatcu /= 0) then
                      jval(ld) = jval(ld) + 2.0d0*(sum( Dimag( prod(1:nor)*Conjg(psi_r))))
                      jxy(ld)  = jxy(ld)  + 2.0d0* sum(        prod(1:nor)*Conjg(psi_r_swap))
                      !jval(ld) = jval(ld) + 2.0d0*Dimag( dot_product( psi_r, prod(1:nor)))
                      !jxy(ld)  = jxy(ld)  + 2.0d0* dot_product( psi_r_swap, prod(1:nor))
                     endif
                     !Orbital currents here, i.e. <psi|H L|psi>
                     if (do_iatLcu /= 0) then

                      lxval(ld)=lxval(ld)+2.0d0*Dimag(dot_product(psi_r, prod_lx(1:nor)))
                      lyval(ld)=lyval(ld)+2.0d0*Dimag(dot_product(psi_r, prod_ly(1:nor)))
                      lzval(ld)=lzval(ld)+2.0d0*Dimag(dot_product(psi_r, prod_lz(1:nor)))

                     endif
                  End Do
                  if (do_iatcu /= 0) then
                   ctp%a(ind, 1, :) = ctp%a(ind, 1, :) - Dimag(jxy(:)) ! x-spin current
                   ctp%a(ind, 2, :) = ctp%a(ind, 2, :) - dble(jxy(:))*dble(spinsign) ! y-spin current
                   ctp%a(ind, 3, :) = ctp%a(ind, 3, :) - jval(:)*dble(spinsign) ! z-spin current
                   ctp%a(ind, 4, :) = ctp%a(ind, 4, :) - jval(:) ! electron current
                   iatcu%cons(:, iar) = iatcu%cons(:, iar) + jval(:)
                   if (do_iatLcu /= 0) then
                    ctp%a(ind, 5, :) = ctp%a(ind, 5, :) + lxval(:) ! these signs need to be double
                    ctp%a(ind, 6, :) = ctp%a(ind, 6, :) - lyval(:) ! checked...
                    ctp%a(ind, 7, :) = ctp%a(ind, 7, :) + lzval(:) !
                   endif
                  endif
               End Do
            End Do
            ofsr1 = ofsr1 + 2*nor
         End Do
    
         if (do_iatLcu /= 0) call free_Lmats(Lop)   


      End Subroutine calc_curr_EMTO
!!$ end of EMTO: calc_curr_EMTO     

 
Subroutine calc_curr (rhs, ham, trgeo, rneibs, itr, iatcu, nl, nr, do_iatcu, do_iatLcu, comm)
         Implicit None
         Type (t_geometry), Intent (In), Target :: trgeo
         Type (zdensemat), Intent (In) :: rhs
         Type (zcsrmat), Intent (In) :: ham
         Type (t_interat_currents), Intent (Inout) :: iatcu ! interatomic currents
         Integer, Intent (In) :: nl, nr, itr, do_iatcu, do_iatLcu
         Integer :: comm
         Type (t_remember_neibs) :: rneibs(:)
!!$           Local vars         
         Complex (Kind=DEF_DBL_PREC) :: hamp (2*trgeo%nlmax**2, 2*trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: psi_r(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: psi_r_swap(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: psi_c(2*trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod(trgeo%nlmax**2) !, prod2(2*trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod_lx(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod_ly(trgeo%nlmax**2)
         Complex (Kind=DEF_DBL_PREC) :: prod_lz(trgeo%nlmax**2)
         Type (t_atom_defenition), Pointer :: pr, pc
         Real (Kind=DEF_DBL_PREC) :: rr(3)
         Integer :: ofsr, ofsr1, ofsc, nor, noc, imd, ld, iar, iac, isr, i1, i2, j, ind
         Integer :: spinsign, swap, tr(2)
         Real (Kind=DEF_DBL_PREC) :: jval(2)
         Real (Kind=DEF_DBL_PREC) :: cons(2)
         Real (Kind=DEF_DBL_PREC) :: lxval(2), lyval(2), lzval(2)
         Complex (Kind=DEF_DBL_PREC) :: jxy(2)
         Type(t_curtensor), Pointer :: ctp
         ! for Lcurs
         Type(t_Lop_mats) :: Lop(trgeo%nlmax)
         Integer :: pc_nl

         if (do_iatLcu /= 0) call make_bigLmat(trgeo%nlmax, Lop)
         
         tr = iatcu%trlist(:, itr)
         ctp => iatcu%curten(itr)
         ofsr1 = 0

         Do iar = 1, trgeo%num ! loop over 'row atoms' of hamiltonian
            pr => trgeo%atoms(iar)%ptr
            rr = trgeo%atoms(iar)%coord(:)
            nor = pr%nl**2
            cons=0d0
!$OMP PARALLEL DO PRIVATE(iac, ofsc, pc_nl, noc, i1, j, i2, hamp, isr, ofsr, spinsign, swap,&
!$OMP& jval, jxy, lxval, lyval, lzval, imd, psi_r, psi_r_swap, psi_c, prod, prod_lx, prod_ly, prod_lz, ld)&
!$OMP& REDUCTION(+:cons)
            Do ind = ctp%ir(iar), ctp%ir(iar+1)-1
               iac = ctp%jc(ind)
               ofsc = rneibs(iac)%ofs 
               !temp pc => trgeo%atoms(iac)%ptr
               pc_nl = trgeo%atoms(iac)%ptr%nl
               !temp noc = 2*pc%nl*pc%nl
               noc = 2*pc_nl*pc_nl
               ! get the relevant part of the hamiltonian
               hamp = 0.0d0

               Do i1 = 1, 2*nor
                  Do j = ham%ir(ofsr1+i1), ham%ir(ofsr1+i1+1)-1
                     If (ham%jc(j) > ofsc ) Then
                        i2 = ham%jc(j)-ofsc
                        If ( i2 > noc ) Exit
                        !hamp(i2,i1) = ham%a(j) ! (transposed)
                        hamp(i1,i2) = ham%a(j) 
                     End If
                  End Do
               End Do

               Do isr = 0, 1 ! row spin (0 = up, 1 = down)
                ofsr = ofsr1 + isr*nor
                spinsign = sign(1,-isr)
                swap = spinsign*nor
                jval = 0.0d0
                jxy = DEF_cmplx_zero
                lxval = 0.0d0
                lyval = 0.0d0
                lzval = 0.0d0
                ! Do the rest
                Do imd = 1, nl + nr
                 !temp psi_r(1:nor)      = rhs%bl(ofsr+1     :ofsr+nor     , imd)
                 !temp psi_r_swap(1:nor) = rhs%bl(ofsr+1+swap:ofsr+nor+swap, imd)
                 !temp psi_c(1:noc)      = rhs%bl(ofsc+1     :ofsc+noc     , imd)
                 psi_r(1:nor)      = rhs%bl(ofsr+1     :ofsr+nor     , imd)
                 psi_r_swap(1:nor) = rhs%bl(ofsr+1+swap:ofsr+nor+swap, imd)
                 psi_c(1:noc)      = rhs%bl(ofsc+1     :ofsc+noc     , imd)

                 !Calculate H|iac>
                 if (do_iatcu /= 0) prod(1:nor) = matmul ( hamp(1+isr*nor:nor+isr*nor, 1:noc), psi_c(1:noc) )
                 !if (do_iatcu /= 0) then
                 !        !temp prod(1:nor) = matmul ( hamp(1+isr*nor:nor+isr*nor, 1:noc), psi_c(1:noc) )
                 !        prod(1:nor) = matmul ( hamp(1+isr*nor:nor+isr*nor, 1:noc), rhs%bl(ofsc+1:ofsc+noc     , imd))
                 !endif
                 if (do_iatLcu /= 0) then
                  !temp prod2(1:noc) = matmul ( conjg(transpose(hamp(1+isr*nor:nor+isr*nor, 1:noc))), psi_r(1:nor) )
                  !temp prod_lx(1:noc) = matmul(Lop(pc%nl)%Lx(1:noc,1:noc), prod2(1:noc) ) 
                  !temp prod_ly(1:noc) = matmul(Lop(pc%nl)%Ly(1:noc,1:noc), prod2(1:noc) ) 
                  !temp prod_lz(1:noc) = matmul(Lop(pc%nl)%Lz(1:noc,1:noc), prod2(1:noc) ) 
                  !temp prod2(1:noc) = matmul ( conjg(transpose(hamp(1+isr*nor:nor+isr*nor,1:noc))), rhs%bl(ofsr+1:ofsr+nor, imd) )
                  !temp prod_lx(1:noc) = matmul(Lop(pc%nl)%Lx(1:noc,1:noc), prod2(1:noc) ) 
                  !temp prod_ly(1:noc) = matmul(Lop(pc%nl)%Ly(1:noc,1:noc), prod2(1:noc) ) 
                  !temp prod_lz(1:noc) = matmul(Lop(pc%nl)%Lz(1:noc,1:noc), prod2(1:noc) ) 
                  !prod_lx(1:noc) = matmul(Lop(pc_nl)%Lx(1:noc,1:noc), prod2(1:noc) ) 
                  !prod_ly(1:noc) = matmul(Lop(pc_nl)%Ly(1:noc,1:noc), prod2(1:noc) ) 
                  !prod_lz(1:noc) = matmul(Lop(pc_nl)%Lz(1:noc,1:noc), prod2(1:noc) ) 
                  prod_lx(1:nor) = matmul(hamp(1+isr*nor:nor+isr*nor, 1:noc), matmul(Lop(pc_nl)%Lx(1:noc,1:noc), psi_c(1:noc) ) )
                  prod_ly(1:nor) = matmul(hamp(1+isr*nor:nor+isr*nor, 1:noc), matmul(Lop(pc_nl)%Ly(1:noc,1:noc), psi_c(1:noc) ) )
                  prod_lz(1:nor) = matmul(hamp(1+isr*nor:nor+isr*nor, 1:noc), matmul(Lop(pc_nl)%Lz(1:noc,1:noc), psi_c(1:noc) ) )
                 endif

                 If (imd > nl) Then
                    ld = 2
                 Else
                    ld = 1
                 End If

                 if (do_iatcu /= 0) then
                  jval(ld) = jval(ld) + 2.0d0* Dimag(dot_product(psi_r(1:nor), prod(1:nor)))
                  jxy(ld)  = jxy(ld)  + 2.0d0*  dot_product(psi_r_swap(1:nor), prod(1:nor))
                  !jval(ld) = jval(ld) + 2.0d0* Dimag(dot_product(rhs%bl(ofsr+1:ofsr+nor, imd), prod(1:nor)))
                  !jxy(ld)  = jxy(ld)  + 2.0d0*  dot_product(rhs%bl(ofsr+1+swap:ofsr+nor+swap, imd), prod(1:nor))
                 endif
                 !Orbital currents here, i.e. <psi|H L|psi>
                 if (do_iatLcu /= 0) then

                  lxval(ld)=lxval(ld)+2.0d0*Dimag(dot_product(psi_r(1:nor), prod_lx(1:nor)))
                  lyval(ld)=lyval(ld)+2.0d0*Dimag(dot_product(psi_r(1:nor), prod_ly(1:nor)))
                  lzval(ld)=lzval(ld)+2.0d0*Dimag(dot_product(psi_r(1:nor), prod_lz(1:nor)))
                  !lxval(ld)=lxval(ld)-2.0d0*Dimag(dot_product(rhs%bl(ofsc+1:ofsc+noc, imd), prod_lx(1:noc)))
                  !lyval(ld)=lyval(ld)-2.0d0*Dimag(dot_product(rhs%bl(ofsc+1:ofsc+noc, imd), prod_ly(1:noc)))
                  !lzval(ld)=lzval(ld)-2.0d0*Dimag(dot_product(rhs%bl(ofsc+1:ofsc+noc, imd), prod_lz(1:noc)))

                 endif
                End Do
                if (do_iatcu /= 0) then
                 ctp%a(ind, 1, :) = ctp%a(ind, 1, :) - Dimag(jxy(:)) ! x-spin current
                 ctp%a(ind, 2, :) = ctp%a(ind, 2, :) - dble(jxy(:))*dble(spinsign) ! y-spin current
                 ctp%a(ind, 3, :) = ctp%a(ind, 3, :) - jval(:)*dble(spinsign) ! z-spin current
                 ctp%a(ind, 4, :) = ctp%a(ind, 4, :) - jval(:) ! electron current
!                 iatcu%cons(:, iar) = iatcu%cons(:, iar) + jval(:)
                 cons = cons + jval(:)
                 if (do_iatLcu /= 0) then
                  ctp%a(ind, 5, :) = ctp%a(ind, 5, :) + lxval(:) ! somewhere, something is wrong with the 
                  ctp%a(ind, 6, :) = ctp%a(ind, 6, :) - lyval(:) ! signs...
                  ctp%a(ind, 7, :) = ctp%a(ind, 7, :) + lzval(:) !
                 endif
                endif
               End Do
            End Do
!$OMP END PARALLEL DO
            iatcu%cons(:, iar) = iatcu%cons(:, iar) + cons

            ofsr1 = ofsr1 + 2*nor
         End Do

         if (do_iatLcu /= 0) call free_Lmats(Lop)   
      End Subroutine calc_curr
     

 
      Subroutine make_sk_decomp (sk, scr, kpar, tr, needhops, cfs)
         Use sparselib
         Implicit None
         Type (t_strconst) :: scr
         Type (t_mathop) :: sk
         Integer :: tr(2)
         Real (Kind=DEF_DBL_PREC), Optional :: cfs
         Integer, Optional :: needhops
         Real (Kind=DEF_DBL_PREC) :: kpar (2)
!!$ local
         Integer :: nh

         If ( .Not. present(needhops)) Then
            nh = 0
         Else
            nh = needhops
         End If

         Call fillsk (sk%c, scr%nnz, scr%nrows, scr%nm, scr%main, scr%base, kpar, tr, cfs)
         If (nh > 0) Then
            Call fillsk (sk%l, scr%lnnz, scr%lnrows, scr%nl, scr%lhop, scr%base, kpar, tr, cfs)
            Call fillsk (sk%r, scr%rnnz, scr%rnrows, scr%nr, scr%rhop, scr%base, kpar, tr, cfs)
            sk%nl = sk%l%ncol
            sk%nr = sk%r%ncol
            sk%havehops = 1
         End If
         sk%alloc = 1
         sk%n = sk%c%ncol
         
      Contains
         
         Subroutine fillsk (sk, tnnzi, nrow, na, sites, base, kpar, tr, cfs)
!!$ function that allocates and fills sk
!!$ BUT uses only the hoppings corresponding to translation tr
            Implicit None
            Integer :: tnnzi, nrow, na
            Type (zcsrmat) :: sk
            Type (t_str_site), Target :: sites (:)
            Real (Kind=DEF_DBL_PREC) :: base (2, 2)
            Real (Kind=DEF_DBL_PREC) :: kpar (2)
            Integer :: tr(2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs
!!$ local
            Type (t_str_site), Pointer :: jsite, isite
            Integer :: ia, i, nnz, j, lr, tnnz


            Call alloc (sk, tnnzi, nrow)

            tnnz = 0
            lr = 1
            Do ia = 1, na
               isite => sites (ia)
               j = isite%srow
               sk%ir (lr:j) = tnnz + 1
               lr = j + isite%nrows
               nnz = fillskrows (isite, kpar, sk%a(tnnz+1:), sk%ir(j:j+isite%nrows), sk%jc(tnnz+1:), base, tr, cfs)
               tnnz = tnnz + nnz
            End Do
            sk%nnz = tnnz
            
         End Subroutine fillsk
         
         Function fillskrows (site, k, val, ir, jc, base, tr, cfs_in) Result (nnz)
!!$ function that fills the rows of sk that correspondent to "site"
!!$ BUT uses only the hoppings corresponding to translation tr
            Implicit None
            Type (t_str_site), Target :: site
            Real (Kind=DEF_DBL_PREC) :: k (2), base (2, 2)
            Real (Kind=DEF_DBL_PREC), Optional :: cfs_in
            Complex (Kind=DEF_DBL_PREC) :: val (:)
            Integer :: tr(2)
            Integer :: nnz
            Integer :: jc (:), ir (:)
!!$ local
            Integer :: i, fii, lai, ih
            Integer :: nh, nc, nr, ind(site%nhop), startn(site%nhop+1)
            Integer, allocatable :: colind(:)
            Complex (Kind=DEF_DBL_PREC) :: fact
            Complex (Kind=DEF_DBL_PREC), allocatable :: els (:, :)
            Integer, Pointer :: tp (:, :)
            Type (t_str_block), Pointer :: hop (:)
            Real (Kind=DEF_DBL_PREC) :: cfs

            cfs = 1.0d0
            If (present(cfs_in)) cfs = cfs_in

            tp => site%trpar
            hop => site%hop
            nh = 0
            startn(1) = 1
            Do i = 1, site%nhop
               If (tp(1,i) == tr(1) .and. tp(2,i) == tr(2)) Then
                  nh = nh + 1
                  ind(nh) = i
                  startn(nh+1) = startn(nh) + hop(i)%ncol
               End If
            End Do
            nc = startn(nh+1) - 1
            nr = site%nrows
            If (nc > 0) Then
               allocate(els(nr,nc))
               allocate(colind(nc))
               els = DEF_cmplx_zero
               fii = 1
               lai = 0
               Do i = 1, nh
                  ih = ind(i)
                  lai = lai + hop(ih)%ncol
                  fact = cfs * Exp (DEF_cmplx_Ione*Dot_product(k, tp(1, ih)*base(:, 1)+tp(2, ih)*base(:, 2)))
                  els (:, fii:lai) = els (:, fii:lai) + fact * hop(ind(i))%bl
                  colind(fii:lai) = site%rowind( site%startn(ih) : site%startn(ih) + hop(ih)%ncol - 1 )
                  fii = fii + hop(ih)%ncol
               End Do
               fii = 1
               lai = nc
               Do i = 1, nr
                  val (fii:lai) = els (i, :)
                  jc (fii:lai) = colind(:)
                  fii = fii + nc
                  lai = lai + nc
               End Do
               deallocate(els, colind)
            End If
            Do i = 1, nr
               ir(i+1) = ir(i) + nc
            End Do
            nnz = nr * nc
         End Function fillskrows
      End Subroutine make_sk_decomp
      
End
