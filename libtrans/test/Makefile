USEMPI=YES
include $(PFUNIT_DIR)/PFUNIT.mk
include ../../make.sys

all: tests

tests: test_rofllib test_ando test_omta_df run

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../obj

LFLAGS1= -L../../lib -ltrans -lcutils $(LFLAGS)

test_rofllib_TESTS := test_rofllib.pf
test_rofllib_OTHER_LIBRARIES := -L../../lib -ltrans $(F90FLAGS) $(LFLAGS1)
$(eval $(call make_pfunit_test,test_rofllib))

test_ando_TESTS := test_ando.pf
test_ando_OTHER_LIBRARIES := -L../../lib -ltrans $(F90FLAGS) $(LFLAGS1)
$(eval $(call make_pfunit_test,test_ando))

test_omta_df_TESTS := test_omta_df.pf
test_omta_df_OTHER_LIBRARIES := -L../../lib -ltrans $(F90FLAGS) $(LFLAGS1)
$(eval $(call make_pfunit_test,test_omta_df))

run:
	./test_rofllib
	./test_ando
	./test_omta_df

clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) *.F90
