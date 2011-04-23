#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006-2010 Lehigh University. All Rights Reserved.           #
#                                                                           #
# This software is licensed under the Common Public License. Please see     #
# accompanying file for terms.                                              #
#                                                                           #
#===========================================================================#

## $Id: Makefile.am 372 2006-06-15 14:18:02Z andreasw $

AUTOMAKE_OPTIONS = foreign

########################################################################
#                         unitTest for Symphony                        #
########################################################################

AM_CPPFLAGS = $(OSITESTS_CFLAGS) $(SYMPHONYLIB_CFLAGS)
SYMCPPFLAGS =
SYMPCPPFLAGS =
APPLCPPFLAGS =
SYMPALLLDADD =
LDADD = $(OSITESTS_LIBS) $(SYMPHONYLIB_LIBS)
DEPENDENCIES = $(OSITESTS_DEPENDENCIES) $(SYMPHONYLIB_DEPENDENCIES)
ALLLDADD =

LIBS =
SYMLIB = 
SYMPLIB = 

# Add include files
AM_CPPFLAGS += -DCOIN_HAS_SYMPHONY

# List all additionally required solver and Osi libraries

if CLP_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_CLP__
endif

if CPLEX_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_CPLEX__
endif

if GLPK_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_GLPK__
endif 

if SOPLEX_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_SOPLEX__
endif

if XPRESS_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_XPRESS__
endif

########################################################################
#                            Additional Definitions                    #
########################################################################

ALLLDADD += $(LDADD)
SYMCPPFLAGS += $(AM_CPPFLAGS) $(SYMDEFS)
SYMLIB += ../src/libSym.la ../src/OsiSym/libOsiSym.la

SYMPCPPFLAGS += $(AM_CPPFLAGS) $(SYMPDEFS)
SYMPALLLDADD += $(ALLLDADD)

if SYM_PARALLEL
SYMPCPPFLAGS += -I`$(CYGPATH_W) $(PVM_ROOT)/include`
SYMPALLLDADD += -L`$(CYGPATH_W) $(PVM_ROOT)/lib/$(ARCH)` \
	-lgpvm3 -lpvm3
SYMPLIB += ../src/libSym$(MASTEREXT).la ../src/OsiSym/libOsiSym.la
else
SYMPLIB += $(SYMLIB)
endif

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(top_srcdir)/include` \
		   -I`$(CYGPATH_W) $(top_builddir)/include` \
		   -I`$(CYGPATH_W) $(top_srcdir)/src/OsiSym`

TestLDADD = $(LDADD)
TestLINK = 
if USE_CL
TestLINK += cp -f OsiColCutTest.$(OBJEXT) unitTest-OsiColCutTest.$(OBJEXT); \
	cp -f OsiCutsTest.$(OBJEXT)  unitTest-OsiCutsTest.$(OBJEXT); \
	cp -f OsiRowCutDebuggerTest.$(OBJEXT) unitTest-OsiRowCutDebuggerTest.$(OBJEXT); \
	cp -f OsiRowCutTest.$(OBJEXT) unitTest-OsiRowCutTest.$(OBJEXT); \
	cp -f OsiSolverInterfaceTest.$(OBJEXT) unitTest-OsiSolverInterfaceTest.$(OBJEXT); \
	cp -f OsiSymSolverInterfaceTest.$(OBJEXT) unitTest-OsiSymSolverInterfaceTest.$(OBJEXT); \
	cp -f unitTestSym.$(OBJEXT) unitTest-unitTestSym.$(OBJEXT); \
	$(LIBTOOL) \
	--tag=CXX --mode=link $(CXX) $(AM_CFLAGS) $(CCLCFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
TestLDADD += OsiColCutTest.$(OBJEXT) OsiCutsTest.$(OBJEXT) \
	OsiRowCutDebuggerTest.$(OBJEXT) OsiRowCutTest.$(OBJEXT) \
	OsiSolverInterfaceTest.$(OBJEXT) OsiSymSolverInterfaceTest.$(OBJEXT) \
	 unitTestSym.$(OBJEXT)
else 
TestLINK += $(CXXLINK)
endif

noinst_PROGRAMS = unitTest punitTest

unitTest_SOURCES = \
	OsiSolverInterfaceTest.cpp \
	OsiSymSolverInterfaceTest.cpp \
	unitTestSym.cpp

punitTest_SOURCES = \
	OsiSolverInterfaceTest.cpp \
	OsiSymSolverInterfaceTest.cpp \
	unitTestSym.cpp

unitTest_CPPFLAGS = $(SYMCPPFLAGS)
unitTest_LDADD = $(SYMLIB) $(ALLLDADD)
unitTest_DEPENDENCIES = $(SYMLIB) $(DEPENDENCIES)
unitTest_LDFLAGS = -rpath $(libdir)
unitTest_LINK = $(TestLINK)

punitTest_CPPFLAGS = $(SYMPCPPFLAGS)
punitTest_LDADD = $(SYMPLIB) $(SYMPALLLDADD)
punitTest_DEPENDENCIES = $(SYMPLIB) $(DEPENDENCIES)
punitTest_LDFLAGS = -rpath $(libdir)


# We are using the Cbc solver executable to do the unit test

if COIN_HAS_MIPLIB3
MIPLIBDIR=$(MIPLIB3_DATA)
else
MIPLIBDIR = $(top_srcdir)/miplib3
endif

testFILES = air03 dcmulti egout flugpl khb05250 l152lav \
	lseu mod010 p0033 p0201 stein27 vpm1

test: unitTest$(EXEEXT)
	if test x"$(NETLIB_DATA)" != x; then \
		NETLIBARG="-netlibDir $(NETLIB_DATA) -testOsiSolverInterface";\
	fi;
	./unitTest$(EXEEXT) $$NETLIBARG -mpsDir=`$(CYGPATH_W) $(SAMPLE_DATA) | sed -e 's/\\\\/\\\\\\\\/g'`

fulltest: check_files unitTest$(EXEEXT)
	if test x"$(NETLIB_DATA)" != x; then \
		NETLIBARG="-netlibDir $(NETLIB_DATA) -testOsiSolverInterface";\
	fi;
	./unitTest$(EXEEXT) $$NETLIBARG -T $(MIPLIBDIR)

ptest: punitTest$(EXEEXT)
	if test x"$(NETLIB_DATA)" != x; then \
		NETLIBARG="-netlibDir $(NETLIB_DATA) -testOsiSolverInterface";\
	fi;
	./punitTest$(EXEEXT) $$NETLIBARG -p 2

pfulltest: check_files punitTest$(EXEEXT)
	if test x"$(NETLIB_DATA)" != x; then \
		NETLIBARG="-netlibDir $(NETLIB_DATA) -testOsiSolverInterface";\
	fi;
	./punitTest$(EXEEXT) $$NETLIBARG -T $(MIPLIBDIR) -p 2

check_files:
	  if test -r $(MIPLIBDIR); then \
	   echo "MIBLIB3 directory found..."; \
	  else \
	   echo "Downloading required test files to $(MIPIBDIR)"; \
	   mkdir $(MIPLIBDIR); \
	   for file in $(testFILES); do \
	       (echo $$file; svn cat https://projects.coin-or.org/svn/Data/stable/1.0/miplib3/$$file.gz >> $(MIPLIBDIR)/$$file.gz) \
	   done; \
	  fi; \
	  echo "Uncompressing files"; \
	  cd $(MIPLIBDIR); \
	  for file in $(testFILES); do \
	   if test -r $$file.gz; then \
		(gzip -d $$file.gz) \
	   fi \
	  done

.PHONY: test ptest fulltest pfulltest check_files

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES = *.mps *.mps.gz *.lp test2out *.out.gz *.out
