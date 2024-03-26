include ../make.sys

DIRS= obj

.PHONY: pord par seq clean all

par: all
seq: all

all:
	+@(cd obj && make -f ../src/Makefile )

clean:
	+@for dir in $(DIRS); do \
	(cd $$dir && make -f ../src/Makefile clean;) \
	done \

tags:
	@(cd src && etags *.F90 *.F *.h;)
