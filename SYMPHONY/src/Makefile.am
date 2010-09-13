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

## $Id: Makefile.am 387 2006-06-22 23:15:50Z menal $

AUTOMAKE_OPTIONS = foreign

SUBDIRS = OsiSym

##############################################################################
##############################################################################
# Global source files
##############################################################################
##############################################################################

SYMHEADERS = \
	../include/sym_constants.h \
	../include/sym_macros.h \
	../include/sym_types.h \
	../include/sym_cg.h \
	../include/sym_cg_params.h \
	../include/sym_cg_u.h \
	../include/sym_cp.h \
	../include/sym_cp_params.h \
	../include/sym_cp_u.h \
	../include/sym_dg.h \
	../include/sym_dg_params.h \
	../include/sym_dg_u.h \
	../include/sym_lp.h \
	../include/sym_lp_params.h \
	../include/sym_lp_solver.h \
	../include/sym_lp_u.h \
	../include/sym_master.h \
	../include/sym_master_params.h \
	../include/sym_master_u.h \
	../include/sym_messages.h \
	../include/sym_pack_array.h \
	../include/sym_pack_cut.h \
	../include/sym_proccomm.h \
	../include/sym_proto.h \
	../include/sym_qsort.h \
	../include/sym_return_values.h \
	../include/symphony.h \
	../include/symphony_api.h \
	../include/sym_timemeas.h \
	../include/sym_tm.h \
	../include/sym_tm_params.h \
	../include/sym_primal_heuristics.h \
	../include/sym_prep.h  \
	../include/sym_prep_params.h

if USE_CL
 SYMHEADERS += ../include/sym_win32_time.h
endif


#Add SymConfig.h
SYMINCFILES = ../include/SymConfig.h $(SYMHEADERS)

SYM_MASTER_SRC	= Master/master.c Master/master_wrapper.c Master/master_io.c \
	Master/master_func.c Master/master_prep.c Master/master_prep_base.c Master/master_prep_sr.c 

MASTER_SRC =$(SYM_MASTER_SRC)

MASTER_MAIN_SRC     = Master/main.c

DG_SRC		= DrawGraph/draw_graph.c

SYM_TM_SRC = TreeManager/tm_func.c TreeManager/tm_proccomm.c

if SYM_COMPILE_IN_TM
TM_SRC		= TreeManager/tm_func.c TreeManager/tm_proccomm.c
else
TM_SRC          = TreeManager/treemanager.c TreeManager/tm_func.c \
	TreeManager/tm_proccomm.c
endif

SYM_LP_SRC = LP/lp_solver.c LP/lp_varfunc.c LP/lp_rowfunc.c \
	LP/lp_genfunc.c LP/lp_proccomm.c LP/lp_wrapper.c LP/lp_free.c LP/lp_heuristics.c LP/lp_sp.c
if PSEUDO_COSTS
SYM_LP_SRC += LP/lp_pseudo_branch.c
else
SYM_LP_SRC += LP/lp_branch.c
endif

SYM_CG_SRC = CutGen/cg_func.c CutGen/cg_wrapper.c
SYM_CP_SRC = CutPool/cp_proccomm.c CutPool/cp_func.c CutPool/cp_wrapper.c

if SYM_COMPILE_IN_LP
TM_SRC         += LP/lp_solver.c LP/lp_varfunc.c LP/lp_rowfunc.c \
	LP/lp_genfunc.c
TM_SRC         += LP/lp_proccomm.c LP/lp_wrapper.c LP/lp_free.c
if PSEUDO_COSTS
TM_SRC         += LP/lp_pseudo_branch.c
else
TM_SRC         += LP/lp_branch.c
endif
if SYM_COMPILE_IN_CG
TM_SRC         += CutGen/cg_func.c CutGen/cg_wrapper.c
endif
else
MASTER_SRC += LP/lp_solver.c
endif

if SYM_COMPILE_IN_CP
TM_SRC	       += CutPool/cp_proccomm.c CutPool/cp_func.c CutPool/cp_wrapper.c
endif
if SYM_COMPILE_IN_TM
MASTER_SRC     += $(TM_SRC)
endif

LP_SRC		= LP/lp_solver.c LP/lp_varfunc.c LP/lp_rowfunc.c \
	LP/lp_genfunc.c LP/lp_proccomm.c LP/lp_wrapper.c LP/lp.c LP/lp_free.c LP/lp_heuristics.c LP/lp_sp.c
if PSEUDO_COSTS
LP_SRC         += LP/lp_pseudo_branch.c
else
LP_SRC         += LP/lp_branch.c
endif

if SYM_COMPILE_IN_CG
LP_SRC         += CutGen/cg_func.c CutGen/cg_wrapper.c
endif

CP_SRC		= CutPool/cut_pool.c CutPool/cp_proccomm.c CutPool/cp_func.c \
	CutPool/cp_wrapper.c

CG_SRC		= CutGen/cut_gen.c CutGen/cg_proccomm.c CutGen/cg_func.c \
	CutGen/cg_wrapper.c

QSORT_SRC	= Common/sym_qsort.c 
TIME_SRC	= Common/timemeas.c
PROCCOMM_SRC	= Common/proccomm.c
PACKCUT_SRC	= Common/pack_cut.c
PACKARRAY_SRC	= Common/pack_array.c

BB_SRC = $(MASTER_SRC) $(DG_SRC) $(TM_SRC) $(LP_SRC) $(CP_SRC) $(CG_SRC) \
	 $(QSORT_SRC) $(TIME_SRC) $(PROCCOMM_SRC) $(PACKCUT_SRC) \
	 $(PACKARRAY_SRC) 

ALL_SRC = $(BB_SRC)

ALL_SYM_MASTER = $(SYM_MASTER_SRC) $(SYM_TM_SRC) $(SYM_LP_SRC) $(SYM_CG_SRC) \
	$(SYM_CP_SRC) $(TIME_SRC) $(QSORT_SRC) $(PROCCOMM_SRC) $(PACKCUT_SRC) \
	$(PACKARRAY_SRC) 

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

SYMCPPFLAGS =
PCPPFLAGS =
APPLCPPFLAGS =
APPLPCPPFLAGS =
ALLLDADD =
PLDADD =
#ADDLIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt` \
#	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`
LDADD = $(SYMPHONYLIB_LIBS)
AM_CPPFLAGS = $(COINDEPEND_CFLAGS) -I`$(CYGPATH_W) $(SYMINCDIR)` \
              -I`$(CYGPATH_W) $(ABSBUILDDIR)/include`

if USE_CGL_CUTS
AM_CPPFLAGS += $(CGL_CFLAGS)
#ADDLIBS += \
#	`cat $(CGLOBJDIR)/cgl_addlibs.txt`
endif

# Add include files
#AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMINCDIR)`

# List all additionally required solver and Osi libraries

#if CLP_LP_SOLVER
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
#	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
#  LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
#	$(CLPOBJDIR)/src/libClp.la
#  ADDLIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
if CLP_LP_SOLVER
  AM_CPPFLAGS += $(CLP_CFLAGS) -D__OSI_CLP__
endif

if CPLEX_LP_SOLVER
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
#  LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
  AM_CPPFLAGS += $(CPX_CFLAGS) -D__OSI_CPLEX__
endif

if GLPK_LP_SOLVER
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
#  LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la 
  AM_CPPFLAGS += $(GLPK_CFLAGS) -D__OSI_GLPK__
else 
if USE_GMPL
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)`
  AM_CPPFLAGS += $(GLPK_CFLAGS) -DUSE_GLPMPL
endif
endif

if SOPLEX_LP_SOLVER
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
#  LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
  AM_CPPFLAGS += $(SPX_CFLAGS) -D__OSI_SOPLEX__
endif

if XPRESS_LP_SOLVER
#  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
#  LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
  AM_CPPFLAGS += $(XPR_CFLAGS) -D__OSI_XPRESS__
endif

#LDADD += $(OSIOBJDIR)/src/libOsi.la \
#	$(COINUTILSOBJDIR)/src/libCoinUtils.la

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
#AM_CPPFLAGS += -I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglDuplicateRow` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglLiftAndProject` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglLandP` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglRedSplit` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglSimpleRounding` \
#	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
#	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
#	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
#	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
#	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
#	-I`$(CYGPATH_W) $(CGLOBJDIR)/inc`
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

if SIGHANDLER
AM_CPPFLAGS += -DSIGHANDLER
endif

SYMCPPFLAGS += $(AM_CPPFLAGS) $(SYMDEFS)

if SYM_PARALLEL
PCPPFLAGS += $(AM_CPPFLAGS) $(SYMPDEFS) -I`$(CYGPATH_W) $(PVM_ROOT)/include`
PLDADD += $(LDADD) -L`$(CYGPATH_W) $(PVM_ROOT)/lib/$(ARCH)` \
	-lgpvm3 -lpvm3
APPLPCPPFLAGS += $(PCPPFLAGS) -DUSE_SYM_APPLICATION
endif

APPLCPPFLAGS += $(SYMCPPFLAGS) -DUSE_SYM_APPLICATION

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

LT_LDFLAGS += -rpath $(libdir)

########################################################################
########################################################################

# Compile everything with g++
CC = $(CXX)
CCLD = $(CXX)

#Static link
STATIC_LINK =  $(LIBTOOL) --tag=CXX --mode=link $(CCLD) -all-static\
	$(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

# Link with g++
LINK =
if STATICEXE
LINK +=  $(STATIC_LINK) 
else
LINK += $(LIBTOOL) --tag=CXX --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
endif

SYMLINK =
if USE_CL
SYMLINK += cp -f main.$(OBJEXT) symphony-main.$(OBJEXT); $(LIBTOOL) \
	--tag=CXX --mode=link $(CCLD) $(AM_CFLAGS) $(CCLCFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
else
SYMLINK += $(LINK)
endif

########################################################################
#                               libSym                                 #
########################################################################

EXTRA_LTLIBRARIES = libSym.la libSym_cg.la libSym_lp.la libSym_lp_cg.la \
	libSym_tm.la libSym_tm_lp.la libSym_tm_cp.la libSym_tm_lp_cg.la \
	libSym_tm_lp_cg_cp.la libSym_m_tm.la libSym_m_tm_lp.la \
	libSym_m_tm_cp.la libSym_m_tm_lp_cg.la libSym_m_tm_lp_cg_cp.la \
	libSym_cp.la libSymAppl.la libSymAppl_cg.la libSymAppl_lp.la \
	libSymAppl_lp_cg.la libSymAppl_tm.la libSymAppl_tm_lp.la \
	libSymAppl_tm_cp.la libSymAppl_tm_lp_cg.la libSymAppl_tm_lp_cg_cp.la \
	libSymAppl_m_tm.la libSymAppl_m_tm_lp.la libSymAppl_m_tm_cp.la \
	libSymAppl_m_tm_lp_cg.la libSymAppl_m_tm_lp_cg_cp.la \
	libSymAppl_cp.la libSym_m.la libSymAppl_m.la
#$(WHATTOMAKELIBS) $(WHATTOMAKEBINS)
#if !SYM_PARALLEL

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory

coinlibdir = ${libdir}/coin
coinlib_LTLIBRARIES = $(WHATTOMAKELIBS)

# List all source files for this library, including headers
libSym_la_SOURCES = $(SYMINCFILES) $(ALL_SYM_MASTER)
libSym_m_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
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

libSym_la_CPPFLAGS = $(SYMCPPFLAGS)
libSym_m_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_m_tm_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_m_tm_lp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_m_tm_cp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_m_tm_lp_cg_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_m_tm_lp_cg_cp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_lp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_lp_cg_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_cg_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_cp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_tm_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_tm_lp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_tm_lp_cg_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_tm_lp_cg_cp_la_CPPFLAGS = $(PCPPFLAGS) 
libSym_tm_cp_la_CPPFLAGS = $(PCPPFLAGS) 

libSym_lp_la_LINK = $(STATIC_LINK)
libSym_lp_cg_la_LINK = $(STATIC_LINK)
libSym_cg_la_LINK = $(STATIC_LINK)
libSym_cp_la_LINK = $(STATIC_LINK)
libSym_tm_la_LINK = $(STATIC_LINK)
libSym_tm_lp_la_LINK = $(STATIC_LINK)
libSym_tm_lp_cg_la_LINK = $(STATIC_LINK)
libSym_tm_lp_cg_cp_la_LINK = $(STATIC_LINK)
libSym_tm_cp_la_LINK = $(STATIC_LINK)

libSym_la_LDFLAGS = $(LT_LDFLAGS)
libSym_m_la_LDFLAGS = $(LT_LDFLAGS)
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

########################################################################
#                        symphony program                              #
########################################################################

EXTRA_PROGRAMS = symphony symphony_m symphony_m_tm symphony_m_tm_lp \
	symphony_m_tm_cp \
	symphony_m_tm_lp_cg symphony_m_tm_lp_cg_cp symphony_cp symphony_lp \
	symphony_lp_cg symphony_lp_cg_cp symphony_tm symphony_tm_cp \
	symphony_tm_lp 	symphony_tm_lp_cg symphony_tm_lp_cg_cp 

# # # Name of the executable compiled in this directory.  We want it to be
# # # installed in the 'bin' directory

bin_PROGRAMS = $(WHATTOMAKEBINS)

symphony_SOURCES = $(MAIN_SRC) symphony.h
#symphony_DEPENDENCIES = $(SYMLDADD) libSym.la
symphony_CPPFLAGS = $(SYMCPPFLAGS)
symphony_LDADD = libSym.la $(LDADD) $(ADDLIBS)
symphony_LDFLAGS = -rpath $(libdir)
symphony_LINK = $(SYMLINK)

symphony_m_SOURCES = $(MAIN_SRC) symphony.h
#symphony_m_DEPENDENCIES = $(LDADD) libSym$(MASTEREXT).la
symphony_m_CPPFLAGS = $(PCPPFLAGS)
symphony_m_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_LDFLAGS = -rpath $(libdir)

symphony_m_tm_SOURCES = $(MAIN_SRC) symphony.h
#symphony_m_tm_DEPENDENCIES = $(LDADD) libSym$(MASTEREXT).la
symphony_m_tm_CPPFLAGS = $(PCPPFLAGS)
symphony_m_tm_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_LDFLAGS = -rpath $(libdir)

symphony_m_tm_lp_SOURCES = $(MAIN_SRC) symphony.h
#symphony_m_tm_lp_DEPENDENCIES = $(LDADD) libSym$(MASTEREXT).la
symphony_m_tm_lp_CPPFLAGS = $(PCPPFLAGS)
symphony_m_tm_lp_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_lp_LDFLAGS = -rpath $(libdir)

symphony_m_tm_cp_SOURCES = $(MAIN_SRC) symphony.h
#symphony_m_tm_cp_DEPENDENCIES = $(LDADD) libSym$(MASTEREXT).la
symphony_m_tm_cp_CPPFLAGS = $(PCPPFLAGS)
symphony_m_tm_cp_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_cp_LDFLAGS = -rpath $(libdir)

symphony_m_tm_lp_cg_SOURCES = $(MAIN_SRC) symphony.h
#symphony_m_tm_lp_cg_DEPENDENCIES = $(LDADD) libSym$(MASTEREXT).la
symphony_m_tm_lp_cg_CPPFLAGS = $(PCPPFLAGS)
symphony_m_tm_lp_cg_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_lp_cg_LDFLAGS = -rpath $(libdir)

symphony_m_tm_lp_cg_cp_SOURCES = $(MAIN_SRC) symphony.h
#symphony_m_tm_lp_cg_cp_DEPENDENCIES = $(LDADD) libSym$(MASTEREXT).la
symphony_m_tm_lp_cg_cp_CPPFLAGS = $(PCPPFLAGS)
symphony_m_tm_lp_cg_cp_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_lp_cg_cp_LDFLAGS = -rpath $(libdir)

symphony_tm_SOURCES =
#symphony_tm_DEPENDENCIES = $(LDADD) libSym_tm$(TMEXT).la
symphony_tm_CPPFLAGS = $(PCPPFLAGS)
symphony_tm_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_LDFLAGS = -rpath $(libdir)

symphony_tm_lp_SOURCES =
#symphony_tm_lp_DEPENDENCIES = $(LDADD) libSym_tm$(TMEXT).la
symphony_tm_lp_CPPFLAGS = $(PCPPFLAGS)
symphony_tm_lp_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_lp_LDFLAGS = -rpath $(libdir)

symphony_tm_cp_SOURCES =
#symphony_tm_cp_DEPENDENCIES = $(LDADD) libSym_tm$(TMEXT).la
symphony_tm_cp_CPPFLAGS = $(PCPPFLAGS)
symphony_tm_cp_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_cp_LDFLAGS = -rpath $(libdir)

symphony_tm_lp_cg_SOURCES =
#symphony_tm_lp_cg_DEPENDENCIES = $(LDADD) libSym_tm$(TMEXT).la
symphony_tm_lp_cg_CPPFLAGS = $(PCPPFLAGS)
symphony_tm_lp_cg_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_lp_cg_LDFLAGS = -rpath $(libdir)

symphony_tm_lp_cg_cp_SOURCES =
#symphony_tm_lp_cg_cp_DEPENDENCIES = $(LDADD) libSym_tm$(TMEXT).la
symphony_tm_lp_cg_cp_CPPFLAGS = $(PCPPFLAGS)
symphony_tm_lp_cg_cp_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_lp_cg_cp_LDFLAGS = -rpath $(libdir)

symphony_cp_SOURCES =
#symphony_cp_DEPENDENCIES = $(LDADD) libSym_cp.la
symphony_cp_CPPFLAGS = $(PCPPFLAGS)
symphony_cp_LDADD = libSym_cp.la $(PLDADD)
symphony_cp_LDFLAGS = -rpath $(libdir)

symphony_lp_SOURCES =
#symphony_lp_DEPENDENCIES = $(LDADD) 
symphony_lp_CPPFLAGS = $(PCPPFLAGS)
symphony_lp_LDADD = libSym_lp$(LPEXT).la $(PLDADD)
symphony_lp_LDFLAGS = -rpath $(libdir)

symphony_lp_cg_SOURCES =
#symphony_lp_cg_DEPENDENCIES = $(LDADD) libSym_lp$(LPEXT).la 
symphony_lp_cg__CPPFLAGS = $(PCPPFLAGS)
symphony_lp_cg_LDADD = libSym_lp$(LPEXT).la $(PLDADD)
symphony_lp_cg_LDFLAGS = -rpath $(libdir)

symphony_lp_cg_cp_SOURCES =
#symphony_lp_cg_cp_DEPENDENCIES = $(LDADD) libSym_lp$(LPEXT).la
symphony_lp_cg_cp__CPPFLAGS = $(PCPPFLAGS)
symphony_lp_cg_cp_LDADD = libSym_lp$(LPEXT).la $(PLDADD)
symphony_lp_cg_cp_LDFLAGS = -rpath $(libdir)


########################################################################
#                Application Library                                   #    
########################################################################

libSymAppl_la_SOURCES = $(SYMINCFILES) $(ALL_SYM_MASTER)
libSymAppl_m_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSymAppl_m_tm_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSymAppl_m_tm_lp_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSymAppl_m_tm_cp_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSymAppl_m_tm_lp_cg_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSymAppl_m_tm_lp_cg_cp_la_SOURCES = $(SYMINCFILES) $(ALL_MASTER)
libSymAppl_lp_la_SOURCES = $(SYMINCFILES) $(ALL_LP)
libSymAppl_lp_cg_la_SOURCES = $(SYMINCFILES) $(ALL_LP)
libSymAppl_cg_la_SOURCES = $(SYMINCFILES) $(ALL_CG)
libSymAppl_cp_la_SOURCES = $(SYMINCFILES) $(ALL_CP)
libSymAppl_tm_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSymAppl_tm_lp_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSymAppl_tm_lp_cg_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSymAppl_tm_lp_cg_cp_la_SOURCES = $(SYMINCFILES) $(ALL_TM)
libSymAppl_tm_cp_la_SOURCES = $(SYMINCFILES) $(ALL_TM)

libSymAppl_lp_la_LINK = $(STATIC_LINK)
libSymAppl_lp_cg_la_LINK = $(STATIC_LINK)
libSymAppl_cg_la_LINK = $(STATIC_LINK)
libSymAppl_cp_la_LINK = $(STATIC_LINK)
libSymAppl_tm_la_LINK = $(STATIC_LINK)
libSymAppl_tm_lp_la_LINK = $(STATIC_LINK)
libSymAppl_tm_lp_cg_la_LINK = $(STATIC_LINK)
libSymAppl_tm_lp_cg_cp_la_LINK = $(STATIC_LINK)
libSymAppl_tm_cp_la_LINK = $(STATIC_LINK)

libSymAppl_la_LDFLAGS = $(LT_LDFLAGS) 
libSymAppl_m_la_LDFLAGS = $(LT_LDFLAGS) 
libSymAppl_m_tm_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_m_tm_lp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_m_tm_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_m_tm_lp_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_m_tm_lp_cg_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_lp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_lp_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_tm_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_tm_lp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_tm_lp_cg_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_tm_lp_cg_cp_la_LDFLAGS = $(LT_LDFLAGS)
libSymAppl_tm_cp_la_LDFLAGS = $(LT_LDFLAGS)

libSymAppl_la_CPPFLAGS = $(APPLCPPFLAGS) 
libSymAppl_m_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymAppl_m_tm_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_m_tm_lp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_m_tm_cp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_m_tm_lp_cg_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_m_tm_lp_cg_cp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_lp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_lp_cg_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_cg_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_cp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_tm_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_tm_lp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_tm_lp_cg_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_tm_lp_cg_cp_la_CPPFLAGS = $(APPLPCPPFLAGS)
libSymAppl_tm_cp_la_CPPFLAGS = $(APPLPCPPFLAGS)

########################################################################
#                Other stuff                                           #    
########################################################################

#COIN_WARNINGS = -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion
# Get rid of warning flags
#FIXME-Try to have the following later
#CFLAGS := $(subst $(COIN_WARNINGS), ,$(CFLAGS))
#CXXFLAGS := $(subst $(COIN_WARNINGS), ,$(CXXFLAGS))

rename_obj:
	mv main.obj symphony-main.obj

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../include/symphony.h
#	$(SYMHEADERS)

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = SymConfig.h

install-exec-local:
	cat $(top_srcdir)/include/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)


