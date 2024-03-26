module sk_io
#include "math_def.h"

use atoms_module

type t_char_container
   Character (Len=20) :: str
end type t_char_container

CONTAINS

  subroutine read_slater_atom(atclass,label,ef)    
!!$ Reads S-K on-site parameters and other (lmax, polarization etc.) data
    use fmtlib
    use setuplib
    use logging
    implicit none
    type(t_atom_defenition) :: atclass
    character(len=20) :: label
    real(kind=DEF_DBL_PREC) :: ef
!!$ Local 
    integer :: nspin,lmx
    integer :: is,ishell, nparm
    character(len=10) :: tag
    character(len=2) :: eunit

    call open_setup(file='atoms.sk/'//trim(label))

    call category('OPTIONS ',mandatory=.true.)
    call read_setup('LMX=',lmx,mandatory=.true.) ! Max l value
    call read_setup('NSPIN=',atclass%ns,mandatory=.false., default=1) ! number of spin directions to be read

    nspin=atclass%ns

    allocate(atclass%sk_shift(nspin))
    atclass%sk_shift=0.0d0
!!$  Energy unit
    call read_setup('EUNIT=',eunit,mandatory=.false.,default='Ry')    
!!$ Fermi level for which the parameters are defined
    call read_setup('EF=',atclass%ef,mandatory=.false., default=ef)   
!!$ Additional potential shift
    call read_setup('SK_POT_SHIFT=',atclass%sk_shift(1:nspin),mandatory=.false.)

!!$ Anton's code uses numebr of partial waves (nl=1,2,3...) instead of maximum l (lmx=0,1,2..) convention
    atclass%nl=lmx+1

    allocate(atclass%s(nspin))
    allocate(atclass%p(nspin))
    allocate(atclass%d1(nspin))
    allocate(atclass%d2(nspin))
    
    do is=1,nspin
       tag='ONSITE '
       if(is==2) tag='ONSITE.DOWN '
       call category(trim(tag),mandatory=.true.)
       call read_setup('S=',atclass%s(is),mandatory=.false.,default=0.0d0)
       call read_setup('P=',atclass%p(is),mandatory=.false.,default=0.0d0)
       call read_setup('D1=',atclass%d1(is),mandatory=.false.,default=0.0d0) ! t2g (xy,yz,xz)
       call read_setup('D2=',atclass%d2(is),mandatory=.false.,default=0.0d0) ! eg (3z2,x2-y2)
    end do

    call scale_sk_atom(atclass,eunit)

  end subroutine read_slater_atom


  subroutine scale_sk_atom(atclass,eunit)    
    type(t_atom_defenition) :: atclass
    character(len=2) :: eunit
!$    
    real(kind=DEF_DBL_PREC) :: emult
    integer :: is

    select case(eunit) 
    case('eV','EV','Ev','ev')
       emult=1.0d0/13.6056923d0
    case('Ry','RY','ry','rY')
       emult=1.0d0
    case default 
       write(*,*) 'I do not understand energy unit'
    end select

    atclass%ef=atclass%ef*emult
    atclass%sk_shift=atclass%sk_shift*emult
    do is=1,atclass%ns
       atclass%s(is) =atclass%s(is)*emult
       atclass%p(is) =atclass%p(is)*emult
       atclass%d1(is)=atclass%d1(is)*emult
       atclass%d2(is)=atclass%d2(is)*emult
    end do

  end subroutine scale_sk_atom


  subroutine rep_onsite(iat,atclass,ef,label,repf)
!!$ Report on the values of on-site parameters
    use logging
    implicit none
    integer :: iat, repf
    type(t_atom_defenition) :: atclass
    real(kind=DEF_DBL_PREC) :: ef
    type (t_char_container), pointer :: label(:)
!!$ Local
    integer :: is,i
    real(kind=DEF_DBL_PREC) :: loc_pshift
    character(len=6) :: cwork
    character(len=20) :: cwork1

    If (Log_Level < 6) Return

    do i=1,iat-1
       if(label(i)%str==label(iat)%str) then
          return
       end if
    end do

    Write (repf,*) '------------------------------------------------------'
    Write (repf,'("  OPTIONS: LMX=",i2," NSPIN=",i2," EF=",f6.4" SK_POT_SHIFT=",f6.4,1x,f6.4)') &
         &atclass%nl-1,atclass%ns,atclass%ef,atclass%sk_shift
    
    do is=1,atclass%ns
       cwork='     :'
       if (is==2) cwork='.DOWN:'
       write(repf,'("  ONSITE",a6," S=",f8.5," P=",f8.5," D1=",f8.5," D2&
            &=",f8.5)') cwork, atclass%s(is),atclass%p(is),atclass%d1(is)&
            &,atclass%d2(is)
    end do

    Write (repf,*) '--------------------------------------------------&
         &----'
    call do_log(6,'   ')


  end subroutine rep_onsite

  subroutine allocate_sk(sk,nshell,nspin)
    implicit none
    integer, intent(in) :: nshell, nspin
    type(t_sk_param) :: sk
!!$ Local

    allocate(sk%rshell(nshell+1))

    allocate(sk%ss_s(nshell,nspin))

    allocate(sk%sp_s(nshell,nspin))
    allocate(sk%pp_s(nshell,nspin))
    allocate(sk%pp_p(nshell,nspin))

    allocate(sk%ps_s(nshell,nspin))

    allocate(sk%sd_s(nshell,nspin))
    allocate(sk%pd_s(nshell,nspin))
    allocate(sk%pd_p(nshell,nspin))

    allocate(sk%dd_s(nshell,nspin))
    allocate(sk%dd_p(nshell,nspin))
    allocate(sk%dd_d(nshell,nspin))

    allocate(sk%ds_s(nshell,nspin))
    allocate(sk%dp_s(nshell,nspin))
    allocate(sk%dp_p(nshell,nspin))

  end subroutine allocate_sk

  subroutine read_slater_hopping(icls,jcls,atoms,ortho,ilabel,jlabel)    
    use fmtlib
    use setuplib
    use logging
    implicit none
    integer, intent(in) :: icls,jcls
    Type(t_atoms_set) :: atoms     
    logical :: ortho
    Character (Len=20) :: ilabel,jlabel
!!$ Local 
    type(t_atom_defenition), pointer :: iclass,jclass
    type(t_sk_param), pointer :: skh, sko
    integer :: is,ishell, nparm,nshell,nspin
    character(len=15) :: fileAB,fileBA,filename,sec_name
    character(len=1) :: cshell
    character(len=2) :: runit,eunit
    logical :: existsAB, existsBA, switch
    real(kind=DEF_DBL_PREC) :: rmult,scale(2),emult
    
    iclass=>atoms%at(icls)
    jclass=>atoms%at(jcls)

    nspin=max(iclass%ns,jclass%ns)
    skh=>atoms%skh(icls,jcls)
    skh%ns=nspin
    if(.not.ortho) then
       sko=>atoms%sko(icls,jcls)
       sko%ns=nspin
    end if

!!$ Checking which files exist, atoms/AB, atoms/BA ot both   
    fileAB='atoms.sk/'//trim(ilabel)//trim(jlabel)
    inquire(file=fileAB,EXIST=existsAB)

    fileBA='atoms.sk/'//trim(jlabel)//trim(ilabel)
    inquire(file=fileBA,EXIST=existsBA)


    if(existsAB) then
       filename=fileAB
       switch=.false.
    elseif((existsBA).AND.(icls/=jcls)) then
       filename=fileBA
       switch=.true.       
    else
       call do_log(2,"        Can't find hopping file for "//trim(iclass%label)//' and '//trim(jclass%label))
       call do_log(2,'       Assuming no hopping')
       atoms%skh(icls,jcls)%nohopping=.true.
       atoms%skh(jcls,icls)%nohopping=.true.
       if(.not.ortho) then
          atoms%sko(icls,jcls)%nohopping=.true.
          atoms%sko(jcls,icls)%nohopping=.true.
       end if
       return
    end if
    
    call open_setup(file=trim(filename))

    call category('OPTIONS ',mandatory=.true.)    
    call read_setup('NSHELL=',nshell,mandatory=.false., default=1) ! Number of shells
    call read_setup('RUNIT=',runit,mandatory=.false., default='AU')
    call read_setup('EUNIT=',eunit,mandatory=.false., default='Ry')
    scale=1.0d0
    call read_setup('SCALE=',scale(1:nspin),mandatory=.false.)

    select case(runit)
    case('AU','au','Au','aU') 
       rmult=1.0d0
    case('AN','an','aN','An') 
       rmult=1.0d0/0.52917720859d0
    case default
       write(*,*) 'I do not understand units for shell radiuses'
    end select

    skh%nshell=nshell
    call allocate_sk(skh,nshell,nspin)

    if(.not.ortho) then
       sko%nshell=nshell
       call allocate_sk(sko,nshell,nspin)
    end if

    do ishell=1,nshell+1
       write(cshell,'(i1)') ishell
       call category('SHELL'//cshell//' ',mandatory=.true.)
       call read_setup('RSHELL=',skh%rshell(ishell),mandatory=.true.) ! The radius of i-th shell
       skh%rshell(ishell)=skh%rshell(ishell)*rmult
       if(ishell==nshell+1) exit

       if(.not.ortho)  sko%rshell(ishell)=skh%rshell(ishell)

       do is=1,nspin
          sec_name='HAM:'
          if(is==2) sec_name='HAM.DOWN:'
          call sub_category(trim(sec_name),mandatory=.true.)
          call read_setup('SS_S=',skh%ss_s(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('SP_S=',skh%sp_s(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('PP_S=',skh%pp_s(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('PP_P=',skh%pp_p(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('SD_S=',skh%sd_s(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('PD_S=',skh%pd_s(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('PD_P=',skh%pd_p(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('DD_S=',skh%dd_s(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('DD_P=',skh%dd_p(ishell,is),mandatory=.false.,default=0.0d0)
          call read_setup('DD_D=',skh%dd_d(ishell,is),mandatory=.false.,default=0.0d0)

          call read_setup('PS_S=',skh%ps_s(ishell,is),mandatory=.false.,default=-skh%sp_s(ishell,is))
          call read_setup('DS_S=',skh%ds_s(ishell,is),mandatory=.false.,default= skh%sd_s(ishell,is))
          call read_setup('DP_S=',skh%dp_s(ishell,is),mandatory=.false.,default=-skh%pd_s(ishell,is))
          call read_setup('DP_P=',skh%dp_p(ishell,is),mandatory=.false.,default=-skh%pd_p(ishell,is))

          if(.not.ortho) then
             sec_name='OVP:'
             if(is==2) sec_name='OVP.DOWN:'
             call sub_category(trim(sec_name),mandatory=.true.)
             call read_setup('SS_S=',sko%ss_s(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('SP_S=',sko%sp_s(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('PP_S=',sko%pp_s(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('PP_P=',sko%pp_p(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('SD_S=',sko%sd_s(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('PD_S=',sko%pd_s(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('PD_P=',sko%pd_p(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('DD_S=',sko%dd_s(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('DD_P=',sko%dd_p(ishell,is),mandatory=.false.,default=0.0d0)
             call read_setup('DD_D=',sko%dd_d(ishell,is),mandatory=.false.,default=0.0d0)

             call read_setup('PS_S=',sko%ps_s(ishell,is),mandatory=.false.,default=-sko%sp_s(ishell,is))
             call read_setup('DS_S=',sko%ds_s(ishell,is),mandatory=.false.,default= sko%sd_s(ishell,is))
             call read_setup('DP_S=',sko%dp_s(ishell,is),mandatory=.false.,default=-sko%pd_s(ishell,is))
             call read_setup('DP_P=',sko%dp_p(ishell,is),mandatory=.false.,default=-sko%pd_p(ishell,is))
          end if
       end do
    end do

    call scale_sk_hopping(nspin,skh,eunit,scale(1:nspin))
    if(.not.ortho)     call scale_sk_hopping(nspin,sko,'Ry',scale(1:nspin))

!!$ If the parameters were read from BA file we must switch some
    !! parameters (eg sp and ps) and adjust the signs accordingly
    if(switch) then
       call hswitch(nspin,skh)
          if(.not.ortho) then
             call hswitch(nspin,sko)
          end if
    end if

    skh%nohopping=.false. ! Hoppings are defined
    if(.not.ortho) then
       sko%nohopping=.false.
    end if

  end subroutine read_slater_hopping


  subroutine scale_sk_hopping(nspin,sk,eunit,scale)
    implicit none
    integer :: nspin
    type(t_sk_param) :: sk
    character(len=2) :: eunit
    real(kind=DEF_DBL_PREC) :: scale(1:nspin)
!!$ Local
    real(kind=DEF_DBL_PREC) :: emult(1:nspin)
    integer :: is,ish

    select case(eunit) 
    case('eV','EV','Ev','ev')
       emult=1.0d0/ 13.605698d0*scale
    case('Ry','RY','ry','rY')
       emult=1.0d0*scale
    case default 
       write(*,*) 'I do not understand energy unit'
    end select

    do is=1,nspin
       do ish=1,sk%nshell          
          sk%ss_s(ish,is)=sk%ss_s(ish,is)*emult(is)
          sk%sp_s(ish,is)=sk%sp_s(ish,is)*emult(is)
          sk%pp_s(ish,is)=sk%pp_s(ish,is)*emult(is)
          sk%pp_p(ish,is)=sk%pp_p(ish,is)*emult(is)
          sk%sd_s(ish,is)=sk%sd_s(ish,is)*emult(is)
          sk%pd_s(ish,is)=sk%pd_s(ish,is)*emult(is)
          sk%pd_p(ish,is)=sk%pd_p(ish,is)*emult(is)
          sk%dd_s(ish,is)=sk%dd_s(ish,is)*emult(is)
          sk%dd_p(ish,is)=sk%dd_p(ish,is)*emult(is)
          sk%dd_d(ish,is)=sk%dd_d(ish,is)*emult(is)
          sk%ps_s(ish,is)=sk%ps_s(ish,is)*emult(is)
          sk%ds_s(ish,is)=sk%ds_s(ish,is)*emult(is)
          sk%dp_s(ish,is)=sk%dp_s(ish,is)*emult(is)
          sk%dp_p(ish,is)=sk%dp_p(ish,is)*emult(is)
       end do
    end do

  end subroutine scale_sk_hopping


  subroutine rep_hopping(icls,jcls,atoms,ortho,label,repf)    
    use logging
    implicit none
    integer, intent(in) :: icls,jcls,repf
    Type(t_atoms_set) :: atoms     
    logical :: ortho
    type (t_char_container), pointer :: label(:)
!!$ Local
    type(t_sk_param), pointer :: skh, sko
    integer :: i,j,jend,ns,is,nlmax,ishell,iorth,north
    character(len=20) :: cwork,cwork1

    cwork=trim(trim(label(icls)%str)//trim(label(jcls)%str))

    If (Log_Level < 6) Return

    jend=atoms%num
    do i=1,icls
       if(i==icls) jend= jcls-1
       do j=1,jend
          cwork1=trim(trim(label(i)%str)//trim(label(j)%str))
          if(cwork==cwork1) then
             return
          end if
       end do
    end do

    ns=max(atoms%at(icls)%ns,atoms%at(jcls)%ns)
    nlmax=max(atoms%at(icls)%nl,atoms%at(jcls)%nl)

    north=1
    if(.not.ortho) then
       north=2
    end if
    
    Write (repf,*) '------------------------------------------------------'
    Write (repf,*) ' NSHELL=',atoms%skh(icls,jcls)%nshell
    
    do ishell=1,atoms%skh(icls,jcls)%nshell
       Write (repf,'("  SHELL",i1,"  RSHELL=",f8.5)') ishell,atoms%skh(icls,jcls)%rshell(ishell)
       do is=1,ns
          do iorth=1,north    
             skh=>atoms%skh(icls,jcls)                 
             cwork='HAM'
             if (iorth==2) then
                skh=>atoms%sko(icls,jcls)                 
                cwork='OVP'
             end if
             cwork1='     :'             
             if (is==2) cwork1='.DOWN:'

             Write (repf,'("   ",a9," SS_S=",f8.5)')cwork(1:3)//cwork1(1:6),skh&
                  &%ss_s(ishell,is)
             if(nlmax>1) then
                Write (repf,'("             SP_S=",f8.5," PS_S=",f8.5," PP&
                     &_S",f8.5," PP_P=",f8.5)') skh%sp_s(ishell,is),&
                     & skh%ps_s(ishell,is),skh%pp_s(ishell,is),skh%pp_p(ishell,is)
             end if
             if(nlmax>2) then
                Write (repf,'("             SD_S=",f8.5," DS_S=",f8.5," PD&
                     &_S",f8.5," DP_S",f8.5," PD_P",f8.5," DP_P",f8.5)') skh%sd_s(ishell,is),skh%ds_s(ishell&
                     &,is),skh%pd_s(ishell,is),skh%dp_s(ishell,is),skh%pd_p(ishell,is),skh%dp_p(ishell,is)
                Write (repf,'("             DD_S=",f8.5," DD_P=",f8.5," DD&
                     &_D",f8.5)') skh%dd_s(ishell,is),skh%dd_p(ishell&
                     &,is),skh%dd_d(ishell,is)
             end if
          end do
       end do
    end do
    Write (repf,*) '------------------------------------------------------'
    Write (repf,*) ''
    
  end subroutine rep_hopping

  subroutine hswitch(nspin,sk)
    implicit none
    integer :: nspin
    type(t_sk_param) :: sk
!!$    type(ctrl_opt) :: opt
!!$ Local
    integer :: is,ish

    do is=1,nspin
       do ish=1,sk%nshell          
          call rswitch(sk%sp_s(ish,is),sk%ps_s(ish,is),sign=-1)
          call rswitch(sk%sd_s(ish,is),sk%ds_s(ish,is),sign=+1)
          call rswitch(sk%pd_s(ish,is),sk%dp_s(ish,is),sign=-1)
          call rswitch(sk%pd_p(ish,is),sk%dp_p(ish,is),sign=-1)
       end do
    end do
  end subroutine hswitch

  subroutine rswitch(a,b,sign)
    implicit none
    real(kind=DEF_DBL_PREC) :: a,b,tmp
    integer, optional :: sign
    integer :: lsign

    if(present(sign)) then
       lsign=sign
    else
       lsign=1
    end if

    tmp=a
    a=b*sign
    b=tmp*sign

  end subroutine rswitch


  Function read_slater_atoms (filename,ortho,ef) Result (atoms)
    Use logging
    Implicit None
    Type (t_atoms_set) :: atoms !< structure with atomic parameters returned as result
    Character (Len=*) :: filename
    logical :: ortho
    real(kind=DEF_DBL_PREC) :: ef
!!$ Local
    Integer :: numb
    Character (Len=20) :: cwork !< temporary string variable
    Character (Len=20) :: cwork1 !< temporary string variable
    Integer :: i, ins, nr, iat, jat, nlmax = 0, nrmax = 0, nsmax = 0, swapspins, nl
    Integer, Parameter :: listf = 100, atomf = 101
    Type (t_atom_defenition), Pointer :: at,at1
    type (t_char_container), pointer :: label(:)

    Call do_log (1, 'Reading atom defenitions')
    Open (Unit=listf, File=filename, Action='read')
100 Format (5 a16)
101 Format (1 x, 10 i5)
104 Format (1 x, 4 g15.7)
10414 Format (1 x, 10 g15.7)

!!$ Skip first two lines left for comments
    Read (listf, 100) cwork
    Read (listf, 100) cwork
    Read (listf,*) numb

    atoms%num = numb
    Write (cwork,*) numb
    Call do_log (2, ' Number of atoms = '//trim(cwork))

    Allocate (atoms%at(numb))
    Allocate (atoms%skh(numb,numb))
    if(.not.ortho) Allocate (atoms%sko(numb,numb))

    Allocate (label(numb))

    atoms%nsmax=1
    atoms%nlmax=1

    Do iat = 1, numb

       at => atoms%at (iat)
       Read (listf,*) cwork1, cwork

       label(iat)%str=trim(cwork)
       at%label = trim (cwork1)
       at%label_len = len_trim (cwork1)
       at%index=iat
       at%wsr=1.0d0

       Call do_log (2, ' Reading on site parameters for: '//trim(at%label)//' from file: '//trim(cwork))
       call read_slater_atom(at,cwork,ef)
       call rep_onsite(iat,at,ef,label,log_stream_number) ! Report on the on-site parameters
!!$    Shift on-site parameters by EF-AT%EF+AT%POT_SHIFT
       at%sk_shift=ef-at%ef+at%sk_shift
       at%ef=ef

       atoms%nlmax=max(atoms%nlmax,at%nl)
       atoms%nsmax=max(atoms%nsmax,at%ns)

    end Do

    Call do_log (2, ' Reading hopping parameters:')

    do iat=1,atoms%num
       do jat=1,atoms%num
          at => atoms%at (iat)
          at1 => atoms%at (jat)
          Call do_log (2, '      between '//trim(at%label)//' and '&
               &//trim(at1%label)//' ('//trim(label(iat)%str)&
               &//trim(label(jat)%str)//')') 
          call read_slater_hopping(iat,jat,atoms,ortho,label(iat)%str&
               &,label(jat)%str)
          call rep_hopping(iat,jat,atoms,ortho,label,log_stream_number) ! Report on the hopping S-K parameters
       end do
    end do

    deallocate(label)

    Close (listf)
    Call do_log (1, 'Reading atom defenitions ....Done!')
    Call do_log (1, '')
    Return
  End Function read_slater_atoms


  subroutine sk_shift_pot(at,ef)
!!$ This is wrong for non-orthogonal basis!!!!
    use logging
    implicit none
    Type (t_atom_defenition), Pointer :: at
    real(kind=DEF_DBL_PREC) :: ef
!!$     Local
    integer :: iat,is
    real(kind=DEF_DBL_PREC) :: loc_pshift

    do is=1,at%ns
       loc_pshift=ef-at%ef+at%pot_shift
       at%s(is)=at%s(is)+loc_pshift
       at%p(is)=at%p(is)+loc_pshift
       at%d1(is)=at%d1(is)+loc_pshift
       at%d2(is)=at%d2(is)+loc_pshift       
    end do
    at%ef=ef    
    
  end subroutine sk_shift_pot

end module sk_io
