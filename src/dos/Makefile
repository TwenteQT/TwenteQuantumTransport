include ../../make.sys
include LASTS

OBJ= dos.o 
PROGS = dos

BINPROGS=$(patsubst %,$(BINDIR)%$(SUF).x,$(PROGS))

DEPLIB=../../lib/*$(SUF).a
do_clean=$(subst $(subst $(LSUF),clean_o,$(SUF)),,clean_o)

all: $(do_clean) $(BINPROGS)
	@echo LSUF=$(SUF) > LASTS

#all:
#	echo $(BINPROGS)

clean_o:
	@$(RM) *.o *.mod *.MOD  *~  *.lst

clean:
	@$(RM) *.o *.mod *.MOD  *~ *.lst  $(patsubst %,$(BINDIR)%*,$(PROGS))

.SUFFIXES :
.SUFFIXES : .o .F .F90 .x

%.x: $(OBJ) $(DEPLIB)
	$(FL) $(OBJ) $(F90OPTS)  $(LOPTS)  -o $@


.F90.o: $(DEPLIB)
	$(F90) $(F90OPTS) -c $<

    
include deps.mk
