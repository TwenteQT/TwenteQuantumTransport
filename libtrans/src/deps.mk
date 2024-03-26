omta_defs.o: sparselib.o
ando.o: sparselib.o
atoms.o: heritage.o sparselib.o logging.o
bzgrid.o: logging.o
adaptive_bzgrid.o: logging.o
geometry.o: omta_defs.o atoms.o logging.o
hamiltonian.o: omta_SOC.o geometry.o structure.o potpars.o
hdf5io.o:  sparselib.o geometry.o structure.o bzgrid.o spsolvers.o readcfg.o potpars.o ando.o transport.o supercell.o postprocess.o interatcur.o qofzlib.o OAMdensity.o
postprocess.o : sparselib.o geometry.o structure.o
potpars.o: sparselib.o structure.o rotations.o logging.o geometry.o
readcfg.o: logging.o spsolvers.o hamiltonian.o potpars.o itersolv.o ando.o atoms.o bzgrid.o #petsc_solvers.o
rotations.o: logging.o
setuplib.o: fmtlib.o 
sk_ham.o: bzgrid.o atoms.o geometry.o logging.o spherical.o structure.o
sk_io.o: atoms.o fmtlib.o logging.o setuplib.o
spsolvers.o: logging.o sparselib.o
structure.o: omta_defs.o omta_strrs.o omta_df.o logging.o geometry.o pgsym.o sparselib.o
supercell.o: geometry.o ando.o hamiltonian.o structure.o
transport.o: ando.o supercell.o sparselib.o
itersolv.o: ando.o sparselib.o helpers.o
compat.o:
helpers.o: logging.o
interatcur.o: sparselib.o geometry.o structure.o hamiltonian.o readcfg.o
discursch.o:quantlib.o interatcur.o
qofzlib.o: quantlib.o discursch.o sparselib.o interatcur.o postprocess.o
df.o: sparselib.o potpars.o structure.o geometry.o
omta_kink.o: omta_SOC.o omta_defs.o geometry.o sparselib.o logging.o
omta_strrs.o: omta_defs.o geometry.o logging.o
omta_pots.o: structure.o omta_defs.o omta_strrs.o logging.o
omta_SOC.o: omta_defs.o omta_strrs.o geometry.o sparselib.o logging.o
omta_df.o: sparselib.o potpars.o structure.o geometry.o
omta_nonsph.o: structure.o
rofllib.o: qofzlib.o hdf5io.o
OAMdensity.o: omta_SOC.o
#slepc_solvers.o: sparselib.o
#berry.o: adaptive_bzgrid.o
#petsc_solvers.o: sparselib.o




