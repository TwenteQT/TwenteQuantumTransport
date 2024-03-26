include ../../../make.sys

BINDIR_=../../../bin/
PROGS=nies.x 
QBINDIR=../../../bin/
INC1= -I../../../libtrans/obj -I../../../libtrans/src
LFLAGS1=  -L../../../lib -ltrans -lcutils $(LFLAGS)

BINPROGS=$(patsubst %,$(BINDIR_)%,$(PROGS))

OBJ=       definitions.o structure.o inversion.o matmul_blas.o \
           green_mod.o read_write.o dens_mod.o maco_mod.o pot_mod.o ies.o \
	   sgfva_mod.o


all:$(BINPROGS)

clean: 
	@rm -f *.o *.mod *.MOD  *.lst $(BINPROGS)


.SUFFIXES : .o .F .F90 .c .x

%.o: ../src/%.F90
	$(F90) $(F90FLAGS)  $(INC1) -c ../src/$*.F90
%.o: ../src/%.F
	$(F77) $(F77FLAGS)  $(INC1) -c ../src/$*.F
%.o: ../src/%.c
	$(CC) $(CFLAGS)  $(INC1)  -c ../src/$*.c

$(BINDIR_)%.x:$(OBJ)
	$(F90)  $(OBJ) $(LFLAGS1)  -o $@ 

include ../src/deps.mk
	