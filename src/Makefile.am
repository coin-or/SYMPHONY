#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006-2022 Lehigh University. All Rights Reserved.           #
#                                                                           #
# This software is licensed under the Eclipse Public License. Please see    #
# accompanying file for terms.                                              #
#                                                                           #
#===========================================================================#

## $Id: Makefile.am 387 2006-06-22 23:15:50Z menal $

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

if COIN_CXX_IS_CL
 SYMHEADERS += ../include/sym_win32_time.h
endif

#Add SymConfig.h
SYMINCFILES = ../include/SymConfig.h $(SYMHEADERS)

SYM_MASTER_SRC	= Master/master.cpp Master/master_wrapper.c Master/master_io.c \
		Master/master_func.c Master/master_prep.c \
		Master/master_prep_base.c Master/master_prep_sr.c 

MASTER_SRC =$(SYM_MASTER_SRC)

MASTER_MAIN_SRC     = Master/main.cpp

DG_SRC		= DrawGraph/draw_graph.c

SYM_TM_SRC = TreeManager/tm_func.c TreeManager/tm_proccomm.c

if SYM_COMPILE_IN_TM
TM_SRC		= TreeManager/tm_func.c TreeManager/tm_proccomm.c
else
TM_SRC          = TreeManager/treemanager.c TreeManager/tm_func.c \
	TreeManager/tm_proccomm.c
endif

SYM_LP_SRC = LP/lp_solver.c LP/lp_varfunc.c LP/lp_rowfunc.c \
	LP/lp_genfunc.c LP/lp_proccomm.c LP/lp_wrapper.c LP/lp_free.c LP/lp_sp.c \
	LP/lp_heuristics.c
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

LP_SRC		= LP/lp.cpp LP/lp_solver.c LP/lp_varfunc.c LP/lp_rowfunc.c \
	LP/lp_genfunc.c LP/lp_proccomm.c LP/lp_wrapper.c LP/lp_free.c LP/lp_sp.c \
	LP/lp_heuristics.c
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

USER_APPL_DIR = ../Applications/USER/src/

USER_MASTER_SRC = $(USER_APPL_DIR)/Master/user_master.c 
USER_LP_SRC = $(USER_APPL_DIR)/LP/user_lp.c $(USER_APPL_DIR)/LP/user_lp_branch.c
USER_CG_SRC = $(USER_APPL_DIR)/CutGen/user_cg.c 
USER_CP_SRC = $(USER_APPL_DIR)/CutPool/user_cp.c 
USER_DG_SRC = $(USER_APPL_DIR)/DrawGraph/user_dg.c

USER_SRC = $(USER_APPL_DIR)/Master/user_master.c \
	$(USER_APPL_DIR)/LP/user_lp.c \
	$(USER_APPL_DIR)/LP/user_lp_branch.c \
	$(USER_APPL_DIR)/CutGen/user_cg.c \
	$(USER_APPL_DIR)/CutPool/user_cp.c \
	$(USER_APPL_DIR)/DrawGraph/user_dg.c

########################################################################
#                           Common Flags                               #
########################################################################

COMMON_FLAGS = -I$(top_srcdir)/include \
              -I$(top_builddir)/include \
              -I$(top_srcdir)/Applications/USER/include

COMMON_FLAGS += $(SYM_OPENMP_FLAGS) $(SYMPHONYLIB_CFLAGS)
if COIN_CXX_IS_CL
COMMON_FLAGS += -TP
endif

if CLP_LP_SOLVER
  COMMON_FLAGS += -D__OSI_CLP__
endif

if CPLEX_LP_SOLVER
  COMMON_FLAGS += -D__OSI_CPLEX__
endif

if GLPK_LP_SOLVER
  COMMON_FLAGS += -D__OSI_GLPK__
else 
if USE_GMPL
  COMMON_FLAGS += -DUSE_GLPMPL
endif
endif

if SOPLEX_LP_SOLVER
  COMMON_FLAGS += -D__OSI_SOPLEX__
endif

if XPRESS_LP_SOLVER
  COMMON_FLAGS += -D__OSI_XPRESS__
endif

########################################################################
#                            Additional Definitions                    #
########################################################################

if SENSITIVITY_ANALYSIS
COMMON_FLAGS += -DSENSITIVITY_ANALYSIS
endif

if USE_CGL_CUTS
COMMON_FLAGS += -DUSE_CGL_CUTS
endif

if ROOT_NODE_ONLY
COMMON_FLAGS += -DROOT_NODE_ONLY
endif

if COMPILE_FRAC_BRANCHING
COMMON_FLAGS += -DCOMPILE_FRAC_BRANCHING
endif

if DO_TESTS
COMMON_FLAGS += -DDO_TESTS
endif

if TM_BASIS_TESTS
COMMON_FLAGS += -DTM_BASIS_TESTS
endif

if TRACE_PATH
COMMON_FLAGS += -DTRACE_PATH
endif

if CHECK_CUT_VALIDITY
COMMON_FLAGS += -DCHECK_CUT_VALIDITY
endif

if STATISTICS
COMMON_FLAGS += -DSTATISTICS
endif

if PSEUDO_COSTS
COMMON_FLAGS += -DPSEUDO_COSTS
endif

if DRAWGRAPH
COMMON_FLAGS += 
endif

if SIGHANDLER
COMMON_FLAGS += -DSIGHANDLER
endif

AM_CPPFLAGS = $(COMMON_FLAGS) $(SYMDEFS)
APPLCPPFLAGS = $(AM_CPPFLAGS) -DUSE_SYM_APPLICATION

if SYM_PARALLEL
PCPPFLAGS = $(COMMON_FLAGS) $(SYMPDEFS) -I$(PVM_ROOT)/include
PLDADD = -lgpvm3 -lpvm3
APPLPCPPFLAGS = $(PCPPFLAGS) -DUSE_SYM_APPLICATION
endif

# Compile everything with g++
CC = $(CXX)
CCLD = $(CXX)
CXXLD = $(CXX)

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
	libSymAppl_cp.la libSym_m.la libSymApplStubs.la libSymAppl_m.la \
	libSymApplStubs_cg.la libSymApplStubs_lp.la libSymApplStubs_m.la \
	libSymApplStubs_lp_cg.la libSymApplStubs_tm.la libSymApplStubs_tm_lp.la \
	libSymApplStubs_tm_cp.la libSymApplStubs_tm_lp_cg.la \
	libSymApplStubs_tm_lp_cg_cp.la \
	libSymApplStubs_m_tm.la libSymApplStubs_m_tm_lp.la \
	libSymApplStubs_m_tm_cp.la \
	libSymApplStubs_m_tm_lp_cg.la libSymApplStubs_m_tm_lp_cg_cp.la \
	libSymApplStubs_cp.la
# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = $(WHATTOINSTALLLIBS)

AM_LDFLAGS = -rpath $(libdir) $(LT_LDFLAGS)

libSym_la_SOURCES = $(SYMINCFILES) $(ALL_SYM_MASTER)
libSym_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
#libSym_la_CPPFLAGS = $(AM_CPPFLAGS)

# List all source files for this library, including headers
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

libSym_m_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_m_tm_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_m_tm_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_m_tm_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_m_tm_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_m_tm_lp_cg_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_tm_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_tm_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_tm_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_tm_lp_cg_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSym_tm_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)

########################################################################
#                        symphony program                              #
########################################################################

bin_PROGRAMS = $(WHATTOMAKEBINS)

EXTRA_PROGRAMS = symphony symphony_m symphony_m_tm symphony_m_tm_lp \
	symphony_m_tm_cp symphony_m_tm_lp_cg symphony_m_tm_lp_cg_cp \
	symphony_tm symphony_tm_lp symphony_tm_cp symphony_tm_lp_cg \
	symphony_tm_lp_cg_cp symphony_cp symphony_lp symphony_lp_cg \
	symphony_lp_cg_cp

symphony_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_LDADD = libSym.la 

symphony_m_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_m_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_CPPFLAGS = $(PCPPFLAGS)

symphony_m_tm_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_m_tm_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_CPPFLAGS = $(PCPPFLAGS)

symphony_m_tm_lp_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_m_tm_lp_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_lp_CPPFLAGS = $(PCPPFLAGS)

symphony_m_tm_cp_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_m_tm_cp_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_cp_CPPFLAGS = $(PCPPFLAGS)

symphony_m_tm_lp_cg_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_m_tm_lp_cg_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_lp_cg_CPPFLAGS = $(PCPPFLAGS)

symphony_m_tm_lp_cg_cp_SOURCES = $(MASTER_MAIN_SRC) symphony.h
symphony_m_tm_lp_cg_cp_LDADD = libSym$(MASTEREXT).la $(PLDADD)
symphony_m_tm_lp_cg_cp_CPPFLAGS = $(PCPPFLAGS)

symphony_tm_SOURCES = 
symphony_tm_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_CPPFLAGS = $(PCPPFLAGS)

symphony_tm_lp_SOURCES = 
symphony_tm_lp_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_lp_CPPFLAGS = $(PCPPFLAGS)

symphony_tm_cp_SOURCES = 
symphony_tm_cp_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_cp_CPPFLAGS = $(PCPPFLAGS)

symphony_tm_lp_cg_SOURCES = 
symphony_tm_lp_cg_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_lp_cg_CPPFLAGS = $(PCPPFLAGS)

symphony_tm_lp_cg_cp_SOURCES = 
symphony_tm_lp_cg_cp_LDADD = libSym_tm$(TMEXT).la $(PLDADD)
symphony_tm_lp_cg_cp_CPPFLAGS = $(PCPPFLAGS)

symphony_cp_SOURCES = 
symphony_cp_LDADD = libSym_cp.la $(PLDADD)
symphony_cp_CPPFLAGS = $(PCPPFLAGS)

symphony_lp_SOURCES = 
symphony_lp_LDADD = libSym_lp.la $(PLDADD)
symphony_lp_CPPFLAGS = $(PCPPFLAGS)

symphony_lp_cg_SOURCES = 
symphony_lp_cg_LDADD = libSym_lp$(LPEXT).la libSym$(MASTEREXT).la $(PLDADD)
symphony_lp_cg_CPPFLAGS = $(PCPPFLAGS)

symphony_lp_cg_cp_SOURCES = 
symphony_lp_cg_cp_LDADD = libSym_lp$(LPEXT).la $(PLDADD)
symphony_lp_cg_cp_CPPFLAGS = $(PCPPFLAGS)

########################################################################
#                Application Library                                   #    
########################################################################

libSymAppl_la_SOURCES = $(SYMINCFILES) $(ALL_SYM_MASTER)
libSymAppl_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_la_CPPFLAGS = $(APPLCPPFLAGS)

# List all source files for this library, including headers
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

libSymAppl_m_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_m_tm_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_m_tm_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_m_tm_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_m_tm_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_m_tm_lp_cg_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_tm_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_tm_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_tm_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_tm_lp_cg_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymAppl_tm_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)

########################################################################
#     A Dummy Library with Stubs for User Callbacks                    #
########################################################################

libSymApplStubs_la_SOURCES  = $(SYMINCFILES) $(USER_APPL_DIR)/include/user.h
libSymApplStubs_la_SOURCES += $(USER_SRC)
libSymApplStubs_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_la_CPPFLAGS = $(APPLCPPFLAGS)

libSymApplStubs_m_la_SOURCES = $(SYMINCFILES) $(USER_MASTER_SRC)
libSymApplStubs_m_tm_la_SOURCES = $(SYMINCFILES) $(USER_MASTER_SRC)
libSymApplStubs_m_tm_lp_la_SOURCES = $(SYMINCFILES) $(USER_MASTER_SRC) $(USER_LP_SRC)
libSymApplStubs_m_tm_cp_la_SOURCES = $(SYMINCFILES) $(USER_MASTER_SRC) $(USER_CP_SRC)
libSymApplStubs_m_tm_lp_cg_la_SOURCES = $(SYMINCFILES) $(USER_MASTER_SRC) $(USER_LP_SRC) \
	$(USER_CG_SRC)
libSymApplStubs_m_tm_lp_cg_cp_la_SOURCES = $(SYMINCFILES) $(USER_MASTER_SRC) $(USER_LP_SRC) \
	$(USER_CG_SRC) $(USER_CP_SRC)
libSymApplStubs_lp_la_SOURCES = $(SYMINCFILES) $(USER_LP_SRC) 
libSymApplStubs_lp_cg_la_SOURCES = $(SYMINCFILES) $(USER_LP_SRC) $(USER_CG_SRC)
libSymApplStubs_cg_la_SOURCES = $(SYMINCFILES) $(USER_CG_SRC)
libSymApplStubs_cp_la_SOURCES = $(SYMINCFILES) $(USER_CP_SRC)
libSymApplStubs_tm_la_SOURCES = $(SYMINCFILES) 
libSymApplStubs_tm_lp_la_SOURCES = $(SYMINCFILES) $(USER_LP_SRC)
libSymApplStubs_tm_lp_cg_la_SOURCES = $(SYMINCFILES) $(USER_LP_SRC) $(USER_CG_SRC)
libSymApplStubs_tm_lp_cg_cp_la_SOURCES = $(SYMINCFILES) $(USER_LP_SRC) $(USER_CG_SRC) \
	$(USER_CP_SRC)
libSymApplStubs_tm_cp_la_SOURCES = $(SYMINCFILES) $(USER_CP_SRC)

libSymApplStubs_m_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_m_tm_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_m_tm_lp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_m_tm_cp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_m_tm_lp_cg_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_m_tm_lp_cg_cp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_lp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_lp_cg_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_cg_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_cp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_tm_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_tm_lp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_tm_lp_cg_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_tm_lp_cg_cp_la_CPPFLAGS = $(APPLPCPPFLAGS) 
libSymApplStubs_tm_cp_la_CPPFLAGS = $(APPLPCPPFLAGS) 

libSymApplStubs_m_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_m_tm_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_m_tm_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_m_tm_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_m_tm_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_m_tm_lp_cg_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_tm_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_tm_lp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_tm_lp_cg_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_tm_lp_cg_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)
libSymApplStubs_tm_cp_la_LIBADD = $(SYMPHONYLIB_LFLAGS) $(SYM_OPENMP_LIBS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = ../include/symphony.h
#	$(SYMHEADERS)

#######################################################################
# Create the Config.h file that has all public defines and install it #
#######################################################################

install-exec-local:
	$(install_sh_DATA) ../include/config_sym.h $(DESTDIR)$(includecoindir)/SymConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/SymConfig.h
