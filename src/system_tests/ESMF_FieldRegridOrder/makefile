# $Id$ 

ALL: tree_build_system_tests

run: tree_run_system_tests

LOCDIR	  = src/system_tests/ESMF_FieldRegridOrder


SYSTEM_TESTS_BUILD    = $(ESMF_TESTDIR)/ESMF_FieldRegridOrderSTest

# Object files other than SysTest%.o that the
# system tests executable will depend on.
# List objects files in the order that they
# are to be compiled/created.
SYSTEM_TESTS_OBJ      = user_model1.o user_model2.o user_coupler.o

SYSTEM_TESTS_RUN      = RUN_FieldRegridOrder

SYSTEM_TESTS_RUN_UNI  = RUN_FieldRegridOrderUNI


ESMF_FieldRegridOrderSTest.o : $(SYSTEM_TESTS_OBJ)

include $(ESMF_DIR)/makefile

DIRS = 

CLEANDIRS   =
CLEANFILES  = $(SYSTEM_TESTS_BUILD)
CLOBBERDIRS =


RUN_FieldRegridOrder:
	$(MAKE) TNAME=FieldRegridOrder NP=6 stest

RUN_FieldRegridOrderUNI:
	$(MAKE) TNAME=FieldRegridOrder NP=1 stest



