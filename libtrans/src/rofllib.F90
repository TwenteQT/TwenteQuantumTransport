! Module to compute resistance (R) as a function of length (L)
! used for calculating the resistivity rho
Module rofllib
  Use qofzlib
  Use HDF5io
Contains
  Subroutine area_from_hdf5(area, filename)
    Implicit None
    !passed
    Double Precision :: area
    Character(len=*) :: filename
    !local
    Integer (HID_T) :: hid
    Integer :: error
    Type (t_dirlist), Target :: dirlist
    Type (t_dirlist), Pointer :: dp
    Integer :: igeo, ngeo
    Logical :: flag
    Double Precision :: base(2,2), sc_size(2), rscale

    ! following code is directly lifted from qofz.F90
    ! should be in qofzlib.F90 as a function in the future
    ! TODO when adding tests to qofzlib.F90
    !----

    dp => dirlist
    ngeo = 0
    Call get_dirlist(dp, ".", 200, ngeo)
    deallocate (dp)
    If (ngeo == 0) Then
       !TODO change this to proper logging 
       !Write (*, '(x,"No ''dirlist'' file found in working directory.")')
       !Write (*, '(x,"Not doing multiple configurations.")')
       Inquire (FILE="./"//filename//".h5", exist=flag)
       If (.not. flag) then
          Write (*, '(x,"Also no '//filename//'.h5 found.")')
          Write (*, '(x,"Stopping.")')
          Stop
       End If
       dirlist%newgeom = 1
       ngeo = 1
    End if
    !---

    ! try to open the first HDF5 file and read the geometry to get the surface
    ! if it fails, it fails. One should change dirlist to get rid of the faulty geometry
    ! Again, these lines of code are more or less identical to the beginning of qofz.F90.
    ! This could (should) be refactored.
    dp => dirlist
    do igeo = 1, ngeo
      call h5open_f(error)
      call h5fopen_f(dp%path//'/'//filename//'.h5', H5F_ACC_RDWR_F, hid, error)
      If (error /= 0) then
        Write (*, '(x,"Error opening ",A)') dp%path//'/'//filename//'.h5'
        Write (*, '(x,"Stopping.")')
        Stop
      Endif

      Call hdfread(hid, '/geom/llgeo/base', base, 2, 2)
      Call hdfread(hid, '/geom/llgeo/sc_size', sc_size, 2)
      Call hdfread(hid, '/geom/llgeo/scale', rscale)

      area = sc_size(1)*sc_size(2)*abs(base(1,1)*base(2,2) - base(1,2)*base(2,1)) / rscale**2

    enddo

  End Subroutine

End Module rofllib
