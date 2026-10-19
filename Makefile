.PHONY: clean libcutils libtrans doc bin

all: bin

bin: libcutils libtrans
	+@(cd src && make)

tags: 
#	+@(cd src && make tags)#
#	+@(cd libtrans && make tags)
	find . -name *.F90 -o -name *.[ch] -o -name *.[fF] | /usr/bin/xargs etags

libcutils:
	+(cd libcutils && make)

libtrans: 
	+@(cd libtrans && make)

tests:
	+@(cd libtrans && make tests)

doc:
	doxygen 

clean: 
	+@(cd libtrans && make clean)
	+@(cd libcutils && make clean)	
	+@(cd src && make clean)
