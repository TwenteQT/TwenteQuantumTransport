
Program rofl
  Implicit None
  Double Precision :: area

  call area_from_hdf5(area) !pulls the area of the lead geometry from the first available HDF5 file 
  
  !open L*/cf-*/transout.h5 file to find area
  !loop over lengths
  !!loop over confgs
  !!!read andout to find resistance
  !!!read conductance from hdf5 file?
  !!endloop
  !compute average resistance for this length
  !endloop
  !
  !compute rho via R(L) = R_contact + rho*L
  !write rho to hdf5
  !write rho to terminal (opt)
  !write rho to file (opt)


End Program
