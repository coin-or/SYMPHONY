# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 387 2006-06-22 23:15:50Z menal $

# Authors: Ted Ralphs, Menal Guzelsoy           

AUTOMAKE_OPTIONS = foreign

SUBDIRS = GMPL

########################################################################
#                               libSym                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libSym.la

# List all source files for this library, including headers
libSym_la_SOURCES = \
	../include/SymConfig.h \
	../include/BB_constants.h \
	../include/BB_macros.h \
	../include/BB_types.h \
	../include/cg.h \
	../include/cg_params.h \
	../include/cg_u.h \
	../include/cp.h \
	../include/cp_params.h \
	../include/cp_u.h \
	../include/dg.h \
	../include/dg_params.h \
	../include/dg_u.h \
	../include/lp.h \
	../include/lp_params.h \
	../include/lp_solver.h \
	../include/lp_u.h \
	../include/master.h \
	../include/master_params.h \
	../include/master_u.h \
	../include/messages.h \
	../include/pack_array.h \
	../include/pack_cut.h \
	../include/proccomm.h \
	../include/proto.h \
	../include/qsortucb.h \
	../include/return_values.h \
	../include/symphony_api.h \
	../include/timemeas.h \
	../include/tm.h \
	../include/tm_params.h \
	Master/master.c \
	Master/master_wrapper.c \
	Master/master_io.c \
	Master/master_func.c \
	TreeManager/tm_func.c \
	TreeManager/tm_proccomm.c \
	LP/lp_solver.c \
	LP/lp_varfunc.c \
	LP/lp_rowfunc.c \
	LP/lp_genfunc.c \
	LP/lp_proccomm.c \
	LP/lp_wrapper.c \
	LP/lp_free.c \
	LP/lp_branch.c \
	CutGen/cg_func.c \
	CutGen/cg_wrapper.c \
	CutPool/cp_proccomm.c \
	CutPool/cp_func.c \
	CutPool/cp_wrapper.c \
	Common/timemeas.c \
	Common/qsortucb.c \
	Common/qsortucb_i.c \
	Common/qsortucb_ii.c \
	Common/qsortucb_id.c \
	Common/qsortucb_di.c \
	Common/qsortucb_ic.c \
	Common/proccomm.c \
	Common/pack_cut.c \
	Common/pack_array.c

# Add Gmpl library

libSym_la_LIBADD = GMPL/libGmpl.la

# This is for libtool (on Windows)
libSym_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                        symphony program                              #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = symphony

# List all source files for this executable, including headers
symphony_SOURCES = \
	Master/main.c ../include/symphony_api.h

# List all required libraries -other than the LP solver's- needed to create
# the executable

symphony_LDADD = libSym.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la \
	$(OSIOBJDIR)/src/libOsi.la \
	$(CGLOBJDIR)/src/libCgl.la

LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt` \
	`cat $(CGLOBJDIR)/cgl_addlibs.txt`

# Add include files
AM_CPPFLAGS = -I`$(CYGPATH_W) ../include`

# List all additionally required solver and Osi libraries
if COIN_HAS_CLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  symphony_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
  AM_CPPFLAGS += -D__OSI_CLP__
endif

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  symphony_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif

# if COIN_HAS_DYLP
#   AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
#         -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
#         -I`$(CYGPATH_W) $(DYLPOBJDIR)/inc` \
# 	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
#   symphony_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la \
# 	$(DYLPOBJDIR)/src/Dylp/libDylp.la
#   LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
# endif

if COIN_HAS_GLPK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  symphony_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif

if COIN_HAS_OSL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  symphony_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif

if COIN_HAS_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  symphony_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif

#if COIN_HAS_SYM
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMSRCDIR)/src` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
#  symphony_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la \
#	$(SYMOBJDIR)/src/libSym.la
#  LIBS += `cat $(SYMOBJDIR)/sym_addlibs.txt`
#endif

if COIN_HAS_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  symphony_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif

LIBS += $(ADDLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
symphony_LDFLAGS = -rpath $(libdir)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
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
	-I`$(CYGPATH_W) $(srcdir)/GMPL`

# List additional defines
AM_CPPFLAGS += -DHAS_RANDOM -DHAS_SRANDOM -DUSE_CGL_CUTS -D__NONE__ -D__LINUX \
	-DUSE_GLPMPL -DCOMPILE_IN_CG -DCOMPILE_IN_CP -DCOMPILE_IN_LP \
	-DCOMPILE_IN_TM

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc


# Compile everything with g++
CC = $(CXX)

# Get rid of C flags
CFLAGS =

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'

# They are already in include dir. 

# include_HEADERS = \
# 	../include/config_sym.h \
# 	../include/symphony_api.h

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = SymConfig.h

# It is already defined in include dir. 

# install-exec-local:
# 	echo "#ifndef HAVE_CONFIG_H" >bla
# 	echo "#define HAVE_CONFIG_H" >>bla
# 	echo "#endif" >> bla
# 	cat $(includedir)/$(ConfigHeader) >> bla
# 	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
# 	rm -f bla

# uninstall-local:
# 	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)
