include ../../make.sys

OBJ= omta_defs.o  ando.o adaptive_bzgrid.o  bzgrid.o    hamiltonian.o  sparselib.o df.o \
     atoms.o  geometry.o omta_kink.o omta_strrs.o logging.o omta_SOC.o supercell.o \
     omta_pots.o omta_nonsph.o omta_df.o structure.o  transport.o r16tomsqz.o \
     r8tomsqz.o spsolvers.o readcfg.o potpars.o rotations.o heritage.o praxis.o \
     pgsym.o postprocess.o hdf5io.o spherical.o sk_io.o sk_ham.o fmtlib.o setuplib.o \
     itersolv.o zPackgmres.o compat.o helpers.o interatcur.o quantlib.o discursch.o qofzlib.o qfitlib.o timer.o \
     OAMdensity.o
     #rofllib.o slepc_solvers.o berry.o OAMdensity.o #petsc_solvers.o 
    
 

LIB=../../lib/libtrans.a

INC1= -I. -I../src

all: $(LIB)

$(OBJ): ../src/math_def.h

clean:
	@$(RM) *.o *.mod *.MOD *.lst  $(LIB) *~

$(LIB): $(OBJ)
	$(AR) $@ $(OBJ)
	@$(RANLIB) $@  

tags:
	@$(TAGS)  *.F90 

#.SUFFIXES :
.SUFFIXES : .o .F .F90 .c

%.o: ../src/%.F90
	$(F90) $(F90FLAGS) $(INC1) -c ../src/$*.F90
%.o: ../src/%.F
	$(F77) $(F90FLAGS) $(INC1) -c ../src/$*.F
%.o: ../src/%.f
	$(F77) $(F77FLAGS) $(INC1) -c ../src/$*.f
%.o: ../src/%.c
	$(CC) $(CFLAGS)  $(INC1) $(CDEFS)  -c ../src/$*.c

include ../src/deps.mk
