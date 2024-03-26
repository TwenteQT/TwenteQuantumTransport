#include "math_def.h"

program anmfix

    Use hdf5io
    implicit none
    Integer (HID_T) :: h5out
    Integer :: error
    
    Call hdf5lib_init ()
    Call h5fopen_f ('transout.h5', H5F_ACC_RDWR_F, h5out, error)
    if (error/=0) then
        write(*,*) 'transout.h5 does not exist!'
        Call h5close_f (error)
        stop
    end if
    
    call fixit(h5out,'lgeo')
    call fixit(h5out,'llgeo')
    call fixit(h5out,'mgeo')
    call fixit(h5out,'rgeo') 
    call fixit(h5out,'rlgeo')  
    call fixit(h5out,'trgeo')
    Call h5fclose_f (h5out, error)
    Call h5close_f (error)
    write(*,*) 'fixed!'
    
   stop
contains
    subroutine fixit(h5out,geonm)
        implicit none
        Integer (HID_T) :: h5out
        Real (Kind=DEF_DBL_PREC) :: anm
        character(len=*) :: geonm
        call hdfread(h5out, '/geom/'//geonm//'/anm', anm)
        if (anm<0.2d0) then
            anm=1/(anm/5.2917720859D-2)*5.2917720859D-2
        end if
        call hdfwrite(h5out, '/geom/'//geonm//'/anm', anm)
        
    end subroutine fixit
end program anmfix
