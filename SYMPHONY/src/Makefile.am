# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 387 2006-06-22 23:15:50Z menal $

# Authors: Ted Ralphs, Menal Guzelsoy           

AUTOMAKE_OPTIONS = foreign

SUBDIRS = 

SYMSOURCES  = \
	$(SYMSRCDIR)/Common      :\
	$(SYMSRCDIR)/LP          :\
	$(SYMSRCDIR)/CutGen      :\
	$(SYMSRCDIR)/CutPool     :\
	$(SYMSRCDIR)/DrawGraph   :\
	$(SYMSRCDIR)/Master      :\
	$(SYMSRCDIR)/../include     :\
	$(SYMSRCDIR)/ :\
	$(SYMSRCDIR)/TreeManager

VPATH  = $(SYMSOURCES)

##############################################################################
##############################################################################
# Global source files
##############################################################################
##############################################################################

SYMINCFILES = \
	SymConfig.h \
	BB_constants.h \
	BB_macros.h \
	BB_types.h \
	cg.h \
	cg_params.h \
	cg_u.h \
	cp.h \
	cp_params.h \
	cp_u.h \
	dg.h \
	dg_params.h \
	dg_u.h \
	lp.h \
	lp_params.h \
	lp_solver.h \
	lp_u.h \
	master.h \
	master_params.h \
	master_u.h \
	messages.h \
	pack_array.h \
	pack_cut.h \
	proccomm.h \
	proto.h \
	qsortucb.h \
	return_values.h \
	symphony_api.h \
	timemeas.h \
	tm.h \
	tm_params.h

MASTER_SRC	= master.c master_wrapper.c master_io.c master_func.c

MASTER_MAIN_SRC     = main.c

DG_SRC		= draw_graph.c

if SYM_COMPILE_IN_TM
TM_SRC		= tm_func.c tm_proccomm.c
else
TM_SRC          = treemanager.c tm_func.c tm_proccomm.c
endif
if SYM_COMPILE_IN_LP
TM_SRC         += lp_solver.c lp_varfunc.c lp_rowfunc.c lp_genfunc.c
TM_SRC         += lp_proccomm.c lp_wrapper.c lp_free.c
if PSEUDO_COSTS
TM_SRC         += lp_pseudo_branch.c
else
TM_SRC         += lp_branch.c
endif
if SYM_COMPILE_IN_CG
TM_SRC         += cg_func.c cg_wrapper.c
endif
else
MASTER_SRC += lp_solver.c
endif

if SYM_COMPILE_IN_CP
TM_SRC	       += cp_proccomm.c cp_func.c cp_wrapper.c
endif
if SYM_COMPILE_IN_TM
MASTER_SRC     += $(TM_SRC)
endif

LP_SRC		= lp_solver.c lp_varfunc.c lp_rowfunc.c lp_genfunc.c \
		  lp_proccomm.c lp_wrapper.c lp.c lp_free.c
if PSEUDO_COSTS
LP_SRC         += lp_pseudo_branch.c
else
LP_SRC         += lp_branch.c
endif

if SYM_COMPILE_IN_CG
LP_SRC         += cg_func.c cg_wrapper.c
endif

CP_SRC		= cut_pool.c cp_proccomm.c cp_func.c cp_wrapper.c

CG_SRC		= cut_gen.c cg_proccomm.c cg_func.c cg_wrapper.c

QSORT_SRC	= qsortucb.c qsortucb_i.c qsortucb_ii.c qsortucb_id.c \
		  qsortucb_di.c qsortucb_ic.c
TIME_SRC	= timemeas.c
PROCCOMM_SRC	= proccomm.c
PACKCUT_SRC	= pack_cut.c
PACKARRAY_SRC	= pack_array.c

BB_SRC = $(MASTER_SRC) $(DG_SRC) $(TM_SRC) $(LP_SRC) $(CP_SRC) $(CG_SRC) \
	 $(QSORT_SRC) $(TIME_SRC) $(PROCCOMM_SRC) $(PACKCUT_SRC) \
	 $(PACKARRAY_SRC)

ALL_SRC = $(BB_SRC)

ALL_MASTER	 = $(MASTER_SRC)
ALL_MASTER 	+= $(TIME_SRC)
ALL_MASTER 	+= $(QSORT_SRC)
ALL_MASTER 	+= $(PROCCOMM_SRC)
ALL_MASTER 	+= $(PACKCUT_SRC)
ALL_MASTER 	+= $(PACKARRAY_SRC)

ALL_LP	 = $(LP_SRC)
ALL_LP 	+= $(TIME_SRC)
ALL_LP 	+= $(QSORT_SRC)
ALL_LP 	+= $(PROCCOMM_SRC)
ALL_LP 	+= $(PACKCUT_SRC)
ALL_LP 	+= $(PACKARRAY_SRC)

ALL_TM	 = $(TM_SRC)
ALL_TM 	+= $(TIME_SRC)
ALL_TM 	+= $(PROCCOMM_SRC)
ALL_TM 	+= $(PACKCUT_SRC)
ALL_TM 	+= $(PACKARRAY_SRC)
if SYM_COMPILE_IN_LP
ALL_TM  += $(QSORT_SRC)
endif

ALL_CG	 = $(CG_SRC)
ALL_CG 	+= $(TIME_SRC)
ALL_CG 	+= $(QSORT_SRC)
ALL_CG 	+= $(PROCCOMM_SRC)
ALL_CG 	+= $(PACKCUT_SRC)

ALL_CP	 = $(CP_SRC)
ALL_CP 	+= $(TIME_SRC)
ALL_CP 	+= $(QSORT_SRC)
ALL_CP 	+= $(PROCCOMM_SRC)
ALL_CP 	+= $(PACKCUT_SRC)

MAIN_SRC = $(MASTER_MAIN_SRC)
##############################################################################
#  Extensions and filenames for various configurations
##############################################################################

# if SYM_COMPILE_IN_CG
# LPEXT = _cg
# endif
# if SYM_COMPILE_IN_CP
# CPEXT = _cp
# endif
# if SYM_COMPILE_IN_LP
# TMEXT = _lp$(LPEXT)$(CPEXT)
# else
# TMEXT = $(CPEXT)
# endif
# if SYM_COMPILE_IN_TM
# MASTEREXT = _m_tm$(TMEXT)
# endif

# LIBNAME = Sym
# MASTERNAME = symphony

# if SYM_PARALLEL
# MASTERLIBNAME = $(LIBNAME)$(MASTEREXT)
# MASTERBIN = $(MASTERNAME)$(MASTEREXT)
# else
# MASTERLIBNAME = $(LIBNAME)
# MASTERBIN = $(MASTERNAME)
# endif

########################################################################
#                           Common Flags                               #
########################################################################

AM_CPPFLAGS = $(SYMDEFS)
LDADD =
LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`

if USE_CGL_CUTS
LDADD += $(CGLOBJDIR)/src/libCgl.la
LIBS += \
	`cat $(CGLOBJDIR)/cgl_addlibs.txt`
endif

# Add include files
AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMINCDIR)`

# List all additionally required solver and Osi libraries

if CLP_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
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

if SYM_PARALLEL
AM_CPPFLAGS += -I`$(CYGPATH_W) $(PVM_ROOT)/include`
LDADD +=  -L$(PVM_ROOT)/lib/$(ARCH) -lgpvm3 -lpvm3
endif

#LIBS += $(ADDLIBS)
LDADD += $(ADDLIBS)

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
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CGLOBJDIR)/inc`
#	-I`$(CYGPATH_W) $(srcdir)/GMPL`

########################################################################
#                            Additional Definitions                    #
########################################################################

if COIN_HAS_READLINE
AM_CPPFLAGS += -DHAS_READLINE
endif

if SENSITIVITY_ANALYSIS
AM_CPPFLAGS += -DSENSITIVITY_ANALYSIS
endif

if USE_CGL_CUTS
AM_CPPFLAGS += -DUSE_CGL_CUTS
endif

if USE_OSI_INTERFACE
AM_CPPFLAGS += 
endif

if ROOT_NODE_ONLY
AM_CPPFLAGS += -DROOT_NODE_ONLY
endif

#if CCMALLOC
#AM_CPPFLAGS += 
#endif

if COMPILE_FRAC_BRANCHING
AM_CPPFLAGS += -DCOMPILE_FRAC_BRANCHING
endif

if DO_TESTS
AM_CPPFLAGS += -DDO_TESTS
endif

if TM_BASIS_TESTS
AM_CPPFLAGS += -DTM_BASIS_TESTS
endif

if TRACE_PATH
AM_CPPFLAGS += -DTRACE_PATH
endif

if CHECK_CUT_VALIDITY
AM_CPPFLAGS += -DCHECK_CUT_VALIDITY
endif

if STATISTICS
AM_CPPFLAGS += -DSTATISTICS
endif

if PSEUDO_COSTS
AM_CPPFLAGS += -DPSEUDO_COSTS
endif

if DRAWGRAPH
AM_CPPFLAGS += 
endif

#FIXME For now, assume that we have gmpl files, later check for GLPK package
if USE_GMPL
AM_CPPFLAGS += -DUSE_GLPMPL

endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

LT_LDFLAGS += -rpath $(libdir)

##############################################################################
##############################################################################
# What to make ? This has to go here in case the user has any targets.
##############################################################################
##############################################################################

# lib_LTLIBRARIES =

#WHATTOMAKELIBS =
#WHATTOMAKEBINS = 

# if SYM_COMPILE_IN_LP
# lib_LTLIBRARIES += lib$(LIBNAME)_lp$(LPEXT).la
# #WHATTOMAKEBINS += $(MASTERNAME)_lp$(LPEXT)
# endif

# if SYM_COMPILE_IN_CP
# lib_LTLIBRARIES += lib$(LIBNAME)_cp.la
# #WHATTOMAKEBINS += $(MASTERNAME)_cp
# endif

# if SYM_COMPILE_IN_CG
# lib_LTLIBRARIES += lib$(LIBNAME)_cg.la
# #WHATTOMAKEBINS += $(MASTERNAME)_cg
# endif

# if SYM_COMPILE_IN_TM
# lib_LTLIBRARIES += lib$(LIBNAME)_tm$(TMEXT).la
# #WHATTOMAKEBINS += $(MASTERNAME)_tm$(TMEXT)
# endif

# lib_LTLIBRARIES += lib$(MASTERLIBNAME).la
# #WHATTOMAKEBINS += $(MASTERNAME)


# if !SYM_COMPILE_IN_LP
# WHATTOMAKELIBS += lib$(LIBNAME)_lp$(LPEXT).la
# WHATTOMAKEBINS += $(MASTERNAME)_lp$(LPEXT)
# endif

# if !SYM_COMPILE_IN_CP
# WHATTOMAKELIBS += lib$(LIBNAME)_cp.la
# WHATTOMAKEBINS += $(MASTERNAME)_cp
# endif

# if !SYM_COMPILE_IN_CG
# WHATTOMAKELIBS += lib$(LIBNAME)_cg.la
# WHATTOMAKEBINS += $(MASTERNAME)_cg
# endif

# if !SYM_COMPILE_IN_TM
# WHATTOMAKELIBS += lib$(LIBNAME)_tm$(TMEXT).la
# WHATTOMAKEBINS += $(MASTERNAME)_tm$(TMEXT)
# endif

# WHATTOMAKELIBS += lib$(MASTERLIBNAME).la
# WHATTOMAKEBINS += $(MASTERNAME)

########################################################################
#                               libSym                                 #
########################################################################

EXTRA_LTLIBRARIES = libSym.la libSym_cg.la libSym_lp.la libSym_lp_cg.la \
	libSym_tm.la libSym_tm_lp.la libSym_tm_cp.la libSym_tm_lp_cg.la \
	libSym_tm_lp_cg_cp.la libSym_m_tm.la libSym_m_tm_lp.la \
	libSym_m_tm_cp.la libSym_m_tm_lp_cg.la libSym_m_tm_lp_cg_cp.la \
	libSym_cp.la 

#$(WHATTOMAKELIBS) $(WHATTOMAKEBINS)
#if !SYM_PARALLEL

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory

lib_LTLIBRARIES = $(WHATTOMAKELIBS)


STATIC_LINK =  $(LIBTOOL) --tag=CXX --mode=link $(CCLD) -static\
	$(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

# List all source files for this library, including headers
libSym_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)

libSym_m_tm_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSym_m_tm_lp_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSym_m_tm_cp_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSym_m_tm_lp_cg_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSym_m_tm_lp_cg_cp_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSym_lp_la_SOURCES = $(SYMINCFILES) $(ALL_LP)
libSym_lp_cg_la_SOURCES = $(SYMINCFILES) $(ALL_LP)
libSym_cg_la_SOURCES = $(SYMINCFILES) $(ALL_CG)
libSym_cp_la_SOURCES = $(SYMINCFILES) $(ALL_CP)
libSym_tm_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSym_tm_lp_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSym_tm_lp_cg_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSym_tm_lp_cg_cp_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSym_tm_cp_la_SOURCES = $(SYMINCFILES) $(ALL_TM)

libSym_lp_la_LINK = $(STATIC_LINK)
libSym_lp_cg_la_LINK = $(STATIC_LINK)
libSym_cg_la_LINK = $(STATIC_LINK)
libSym_cp_la_LINK = $(STATIC_LINK)
libSym_tm_la_LINK = $(STATIC_LINK)
libSym_tm_lp_la_LINK = $(STATIC_LINK)
libSym_tm_lp_cg_la_LINK = $(STATIC_LINK)
libSym_tm_lp_cg_cp_la_LINK = $(STATIC_LINK)
libSym_tm_cp_la_LINK = $(STATIC_LINK)

#nodist_libSym_la_SOURCES = dummy.cpp

libSym_la_LDFLAGS = $(LT_LDFLAGS) 
libSym_m_tm_la_LDFLAGS = $(LT_LDFLAGS)
libSym_m_tm_lp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_m_tm_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_m_tm_lp_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSym_m_tm_lp_cg_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_lp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_lp_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSym_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSym_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_tm_la_LDFLAGS = $(LT_LDFLAGS)
libSym_tm_lp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_tm_lp_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSym_tm_lp_cg_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSym_tm_cp_la_LDFLAGS = $(LT_LDFLAGS)

#else 
#lib_LTLIBRARIES = libSymP.la 

# # if SYM_COMPILE_IN_LP
# # lib_LTLIBRARIES  += lplib lp
# # endif

# # if SYM_COMPILE_IN_CP
# # lib_LTLIBRARIES  += cplib cp
# # WHATTOMAKE += cplib cp
# # endif

# # ifeq ($(SYM_COMPILE_IN_CG),FALSE)
# # lib_LTLIBRARIES  += lplib lp
# # WHATTOMAKE += cglib cg
# # PWHATTOMAKE += pcg
# # QWHATTOMAKE += qcg
# # endif

# # ifeq ($(SYM_COMPILE_IN_TM),FALSE)
# # lib_LTLIBRARIES  += lplib lp
# # WHATTOMAKE += tmlib tm
# # endif


# # List all source files for this library, including headers
# libSymP_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)

# nodist_libSymP_la_SOURCES = dummy.cpp

# # This is for libtool (on Windows)
# libSymP_la_LDFLAGS = $(LT_LDFLAGS)
# endif

#dummy.cpp:
#	echo 'static void dummy();' > $@

########################################################################
#                        symphony program                              #
########################################################################

 EXTRA_PROGRAMS = symphony symphony_m_tm symphony_m_tm_lp symphony_m_tm_cp \
	symphony_m_tm_lp_cg symphony_m_tm_lp_cg_cp symphony_cp symphony_lp \
	symphony_lp_cg symphony_lp_cg_cp symphony_tm symphony_tm_cp \
	symphony_tm_lp 	symphony_tm_lp_cg symphony_tm_lp_cg_cp 

# # # Name of the executable compiled in this directory.  We want it to be
# # # installed in the 'bin' directory
  bin_PROGRAMS = $(WHATTOMAKEBINS)

# LDADD += libSym$(MASTEREXT).la
# LDFLAGS += -rpath $(libdir)

# # List all source files for this executable, including headers
symphony_SOURCES = $(MAIN_SRC) symphony_api.h
symphony_LDADD = libSym$(MASTEREXT).la $(LDADD)
symphony_LDFLAGS = -rpath $(libdir)

symphony_m_tm_SOURCES = $(MAIN_SRC) symphony_api.h
symphony_m_tm_LDADD = libSym$(MASTEREXT).la $(LDADD)
symphony_m_tm_LDFLAGS = -rpath $(libdir)

symphony_m_tm_lp_SOURCES = $(MAIN_SRC) symphony_api.h
symphony_m_tm_lp_LDADD = libSym$(MASTEREXT).la $(LDADD)
symphony_m_tm_lp_LDFLAGS = -rpath $(libdir)

symphony_m_tm_cp_SOURCES = $(MAIN_SRC) symphony_api.h
symphony_m_tm_cp_LDADD = libSym$(MASTEREXT).la $(LDADD)
symphony_m_tm_cp_LDFLAGS = -rpath $(libdir)

symphony_m_tm_lp_cg_SOURCES = $(MAIN_SRC) symphony_api.h
symphony_m_tm_lp_cg_LDADD = libSym$(MASTEREXT).la $(LDADD)
symphony_m_tm_lp_cg_LDFLAGS = -rpath $(libdir)

symphony_m_tm_lp_cg_cp_SOURCES = $(MAIN_SRC) symphony_api.h
symphony_m_tm_lp_cg_cp_LDADD = libSym$(MASTEREXT).la $(LDADD)
symphony_m_tm_lp_cg_cp_LDFLAGS = -rpath $(libdir)

symphony_tm_SOURCES =
symphony_tm_LDADD = libSym_tm$(TMEXT).la $(LDADD)
symphony_tm_LDFLAGS = -rpath $(libdir)

symphony_tm_lp_SOURCES =
symphony_tm_lp_LDADD = libSym_tm$(TMEXT).la $(LDADD)
symphony_tm_lp_LDFLAGS = -rpath $(libdir)

symphony_tm_cp_SOURCES =
symphony_tm_cp_LDADD = libSym_tm$(TMEXT).la $(LDADD)
symphony_tm_cp_LDFLAGS = -rpath $(libdir)

symphony_tm_lp_cg_SOURCES =
symphony_tm_lp_cg_LDADD = libSym_tm$(TMEXT).la $(LDADD)
symphony_tm_lp_cg_LDFLAGS = -rpath $(libdir)

symphony_tm_lp_cg_cp_SOURCES =
symphony_tm_lp_cg_cp_LDADD = libSym_tm$(TMEXT).la $(LDADD)
symphony_tm_lp_cg_cp_LDFLAGS = -rpath $(libdir)

symphony_cp_SOURCES =
symphony_cp_LDADD = libSym_cp.la $(LDADD)
symphony_cp_LDFLAGS = -rpath $(libdir)

symphony_lp_SOURCES =
symphony_lp_LDADD = libSym_lp$(LPEXT).la $(LDADD)
symphony_lp_LDFLAGS = -rpath $(libdir)

symphony_lp_cg_SOURCES =
symphony_lp_cg_LDADD = libSym_lp$(LPEXT).la $(LDADD)
symphony_lp_cg_LDFLAGS = -rpath $(libdir)

symphony_lp_cg_cp_SOURCES =
symphony_lp_cg_cp_LDADD = libSym_lp$(LPEXT).la $(LDADD)
symphony_lp_cg_cp_LDFLAGS = -rpath $(libdir)


########################################################################
#                Application Library                                   #    
########################################################################

if SYM_WITH_APPL

lib_LTLIBRARIES += libSymAppl.la

libSymAppl_la_SOURCES = $(libSym_la_SOURCES)

libSymAppl_la_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_SYM_APPLICATION

# This is for libtool (on Windows)
libSymAppl_la_LDFLAGS = $(LT_LDFLAGS)

endif 

COIN_WARNINGS = -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion
# Compile everything with g++
CC = $(CXX)
CCLD = $(CXX)
# Get rid of warning flags
#FIXME-Try to have the following later
#CFLAGS := $(subst $(COIN_WARNINGS), ,$(CFLAGS))
#CXXFLAGS := $(subst $(COIN_WARNINGS), ,$(CXXFLAGS))

CFLAGS := 

if SYM_DEBUG
CFLAGS += -g
endif

CFLAGS += -pipe
CXXFLAGS := $(CFLAGS)
# Link with g++
LINK = $(LIBTOOL) --tag=CXX --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@

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

CLEANFILES = dummy.cpp
