#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
#                                                                           #
# SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         #
# Laci Ladanyi (ladanyi@us.ibm.com).                                        #
#                                                                           #
# The author of this file is Menal Guzelsoy                                 #
#                                                                           #
# (c) Copyright 2006-2019 Lehigh University. All Rights Reserved.           #
#                                                                           #
# This software is licensed under the Eclipse Public License. Please see    #
# accompanying file for terms.                                              #
#                                                                           #
#===========================================================================#

## $Id: Makefile.am 387 2006-06-22 23:15:50Z menal $

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

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
SYMINCFILES = ../include/SymConfig.h ../include/OsiSymConfig.h $(SYMHEADERS)

MAIN_SRC = Master/main.c

SRC  = Master/master.c Master/master_wrapper.c Master/master_io.c Master/master_func.c 
SRC += Master/master_prep.c Master/master_prep_base.c Master/master_prep_sr.c 

#SRC += DrawGraph/draw_graph.c

SRC += TreeManager/tm_func.c TreeManager/tm_proccomm.c

SRC += LP/lp_solver.c LP/lp_varfunc.c LP/lp_rowfunc.c LP/lp_genfunc.c LP/lp_branch.c
SRC += LP/lp_proccomm.c LP/lp_wrapper.c LP/lp.c LP/lp_free.c LP/lp_sp.c LP/lp_heuristics.c

SRC += CutGen/cg_func.c CutGen/cg_wrapper.c

SRC += CutPool/cp_proccomm.c CutPool/cp_func.c CutPool/cp_wrapper.c

SRC += Common/sym_qsort.c 
SRC += Common/timemeas.c
SRC += Common/proccomm.c
SRC += Common/pack_cut.c
SRC += Common/pack_array.c

USER_APPL_DIR = ../Applications/USER/src/

USER_APPL_STUBS_SRC = $(USER_APPL_DIR)/Master/user_master.c \
	$(USER_APPL_DIR)/LP/user_lp.c \
	$(USER_APPL_DIR)/LP/user_lp_branch.c \
	$(USER_APPL_DIR)/CutGen/user_cg.c \
	$(USER_APPL_DIR)/CutPool/user_cp.c \
	$(USER_APPL_DIR)/DrawGraph/user_dg.c

########################################################################
#                           Common Flags                               #
########################################################################

AM_CPPFLAGS = -DSYMPHONY_BUILD $(SYMPHONYLIB_CFLAGS)

if CLP_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_CLP__
endif

if CPLEX_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_CPLEX__
endif

if GLPK_LP_SOLVER
  AM_CPPFLAGS += -D__OSI_GLPK__
else 
if USE_GMPL
  AM_CPPFLAGS += -DUSE_GLPMPL
endif
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

AM_CPPFLAGS += $(SYMDEFS)
AM_CPPFLAGS += -I`$(CYGPATH_W) $(top_srcdir)/include` \
               -I`$(CYGPATH_W) $(top_builddir)/include` \
               -I`$(CYGPATH_W) $(top_srcdir)/Applications/USER/include`

# Compile everything with g++
CC = $(CXX)
CCLD = $(CXX)

########################################################################
#                               libSym                                 #
########################################################################

lib_LTLIBRARIES = libSym.la libSymApplStubs.la
libSym_la_SOURCES = $(SYMINCFILES) $(SRC)
libSym_la_CPPFLAGS =  -DSYMPHONY_BUILD $(SYMPHONYLIB_CFLAGS) $(AM_CPPFLAGS)
libSym_la_LIBADD = $(SYMPHONYLIB_LFLAGS)

########################################################################
#                        symphony program                              #
########################################################################

bin_PROGRAMS = symphony
symphony_SOURCES = $(MAIN_SRC) symphony.h
symphony_CPPFLAGS = $(SYMPHONYLIB_CFLAGS) $(AM_CPPFLAGS)
symphony_LDADD = libSym.la 

########################################################################
#     A Dummy Library with Stubs for User Callbacks                    #
########################################################################

libSymApplStubs_la_SOURCES  = $(SYMINCFILES) $(USER_APPL_DIR)/include/user.h
libSymApplStubs_la_SOURCES += $(USER_APPL_STUBS_SRC)
libSymApplStubs_la_CPPFLAGS = $(SYMPHONYLIB_CFLAGS) $(AM_CPPFLAGS)

rename_obj:
	mv main.obj symphony-main.obj

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = ../include/symphony.h
#	$(SYMHEADERS)

#######################################################################
# Create the Config.h file that has all public defines and install it #
#######################################################################

install-exec-local:
	$(install_sh_DATA) ../include/config_sym.h $(DESTDIR)$(includecoindir)/SymConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/SymConfig.h
