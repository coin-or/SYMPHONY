#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY Branch, Cut, and Price Library.         #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006 Lehigh University. All Rights Reserved.                #
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

test_srcdir = $(OSISRCDIR)/test

AM_CPPFLAGS =
SYMCPPFLAGS =
APPLCPPFLAGS =
LDADD =
ALLLDADD =

LIBS =

ADDLIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`

if USE_CGL_CUTS
LDADD += $(CGLOBJDIR)/src/libCgl.la
ADDLIBS += \
	`cat $(CGLOBJDIR)/cgl_addlibs.txt`
endif

# Add include files
AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`

AM_CPPFLAGS += -DCOIN_HAS_SYM

# List all additionally required solver and Osi libraries

if CLP_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CLPOBJDIR)/src/libClp.la
  ADDLIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
  AM_CPPFLAGS += -D__OSI_CLP__
endif

if CPLEX_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
  AM_CPPFLAGS += -D__OSI_CPLEX__
endif

if GLPK_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la 
  AM_CPPFLAGS += -D__OSI_GLPK__
else 
if USE_GMPL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)`
endif
endif

if OSL_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
  AM_CPPFLAGS += -D__OSI_OSL__
endif

if SOPLEX_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
AM_CPPFLAGS += -D__OSI_SPX__
endif

if XPRESS_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
  AM_CPPFLAGS += -D__OSI_XPRESS__
endif

LDADD += $(OSIOBJDIR)/src/libOsi.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS += -I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglDuplicateRow` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglLiftAndProject` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglRedSplit` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglSimpleRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CGLOBJDIR)/inc`

########################################################################
#                            Additional Definitions                    #
########################################################################

ALLLDADD += $(LDADD) $(ADDLIBS)
SYMCPPFLAGS += $(AM_CPPFLAGS) $(SYMDEFS)

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

noinst_PROGRAMS = unitTest

unitTest_SOURCES = \
	$(test_srcdir)/OsiColCutTest.cpp \
	$(test_srcdir)/OsiCutsTest.cpp \
	$(test_srcdir)/OsiRowCutDebuggerTest.cpp \
	$(test_srcdir)/OsiRowCutTest.cpp \
	OsiSolverInterfaceTest.cpp \
	OsiSymSolverInterfaceTest.cpp \
	unitTestSym.cpp

unitTest_DEPENDENCIES = $(LDADD)
unitTest_CPPFLAGS = $(SYMCPPFLAGS)
unitTest_LDADD = $(SYMPHONYOBJDIR)/src/libSym.la $(ALLLDADD)
unitTest_LDFLAGS = -rpath $(libdir)


# We are using the Cbc solver executable to do the unit test

test: unitTest$(EXEEXT)
	cd ../../Data/miplib3; make uncompress
	cd ../../Data/Netlib; make uncompress
	./unitTest$(EXEEXT) -T $(DATAOBJDIR)/miplib3 \
	-mpsDir=$(DATAOBJDIR)/Sample \
	-netlibDir=$(DATAOBJDIR)/Netlib
.PHONY: test

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES = *.mps *.mps.gz *.lp test2out *.out.gz
