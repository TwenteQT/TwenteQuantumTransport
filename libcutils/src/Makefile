include ../../make.sys

OBJ= c_mkdir.o memstat.o c_fileio.o c_conv_ptr.o
LIB= ../../lib/libcutils.a

INCS= -I. -I../src 


all: $(LIB)

clean:
	@$(RM) *.o *.mod *.MOD *.lst $(LIB) *~

$(OBJ):../src/*.h

$(LIB): $(OBJ)
	@$(AR) $@ $?
	@$(RANLIB) $@ 

tags:
	@$(TAGS)  *.h *.c

.SUFFIXES :
.SUFFIXES : .o .c


%.o: ../src/%.c
	$(CC) $(CFLAGS)  $(INC) $(CDEFS) $(INCS)  -c ../src/$*.c

include ../src/deps.mk
