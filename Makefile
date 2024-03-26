.PHONY: clean libcutils libtrans doc bin

all: bin

bin: libcutils libtrans
	+@(cd src && make)

tags: 
	+@(cd src && make tags)
	+@(cd libtrans && make tags)

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
