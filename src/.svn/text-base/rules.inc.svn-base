include ../../VERSION

BINDIR=../../bin/
BINPROGS=$(patsubst %,$(BINDIR)%.x,$(PROGS)) 

all: $(OBJ) $(BINPROGS)

.SUFFIXES :
.SUFFIXES : .o .F .F90 .f .x .cxx .c

INC1=-I../../libtrans/obj -I../../libtrans/src -I../../libcutils/src -I./

DEPS1=../../lib/*.a

ifdef C_CODE
    LINKER=$(CL)
else
    LINKER=$(FL)
endif

LFLAGS1= -L../../lib -ltrans -lcutils $(LFLAGS)
F90FLAGS1=$(F90FLAGS) $(INC1) $(FDEFS) $(FDP)_VERS_=\'$(VERSION)\' 
F77FLAGS1=$(F77FLAGS) $(INC1) $(FDEFS) $(FDP)_VERS_=\'$(VERSION)\' 
CFLAGS1=$(CFLAGS)  $(INC1) $(CDEFS) $(CDP)_VERS_=\'$(VERSION)\' 

$(BINDIR)%.x:%.F90 %.o  $(DEPS1)
	$(FL) $< $(OBJ)  $(F90FLAGS1) $(LFLAGS1) -o $@

$(BINDIR)%.x:%.F  %.o $(DEPS1)
	$(FL) $< $(OBJ) $(F90FLAGS1) $(LFLAGS1) -o $@

$(BINDIR)%.x:%.f  %.o $(DEPS1) $(DEPS)
	$(FL) $< $(OBJ) $(F77FLAGS1) $(LFLAGS1) -o $@

$(BINDIR)%.x:%.c  %.o $(DEPS1)
	$(CL) $< $(OBJ) $(CFLAGS1) $(LFLAGS1) -o $@

.F90.o: $(DEPS1) $(DEPS)
	$(F90) $(F90FLAGS1) -c $<

.F.o: $(DEPS1) $(DEPS)
	$(F90) $(F90FLAGS1) -c $<

.f.o: $(DEPS1) $(DEPS)
	$(F77) $(F77FLAGS1) -c $<

.c.o:   $(DEPS1) $(DEPS)
	$(CC) $(CFLAGS1) -c $<

.cxx.o:   $(DEPS1) $(DEPS)
	$(CC) $(CFLAGS1) -c $<

clean:
	@$(RM) *.o *.mod *.MOD   *.lst *~  $(BINPROGS)

include deps.mk

