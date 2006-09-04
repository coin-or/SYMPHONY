# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# $Id$

MASTERNAME = match

EXE = $(MASTERNAME)

# CHANGEME: Additional sources
ADDSOURCES =

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = 


# Directory with COIN header files
COININCDIR = /home/menal/Sym-Dev/include

# Directory with COIN libraries
COINLIBDIR = /home/menal/Sym-Dev/lib

# Directoy with SYMPHONY header files 
SYMLIBDIR = /home/menal/Sym-Dev/SYMPHONY/src

# Directoy with SYMPHONY library
SYMINCDIR = /home/menal/Sym-Dev/SYMPHONY/include

#Directory for LP solver
GLPKINCDIR = /home/menal/GLPK/glpk-4.11/include
CPXINCDIR = 
OSLINCDIR = 
SPXINCDIR = 
XPRINCDIR = 

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile the      #
#  COIN package.                                                         #
##########################################################################

SYM_COMPILE_IN_TM = TRUE
SYM_COMPILE_IN_LP = TRUE
SYM_COMPILE_IN_CP = TRUE
SYM_COMPILE_IN_CG = TRUE
USE_CGL_CUTS = TRUE
CLP_LP_SOLVER = #TRUE
GLPK_LP_SOLVER = TRUE
CPLEX_LP_SOLVER = #TRUE
OSL_LP_SOLVER = #TRUE
SOPLEX_LP_SOLVER = #TRUE
XPRESS_LP_SOLVER = #TRUE
SYM_PARALLEL = #TRUE
USE_GMPL = TRUE
MASTEREXT = _m
LPEXT = _cg
TMEXT = _lp_cg_cp

USER_LP_OBJS       = src/LP/user_lp_branch.o src/LP/user_lp.o 
ifeq ($(SYM_COMPILE_IN_CG),TRUE)
USER_LP_OBJS      += src/CutGen/user_cg.o
endif

USER_TM_OBJS       =
ifeq ($(SYM_COMPILE_IN_LP),TRUE)
USER_TM_OBJS      += $(USER_LP_OBJS)
endif
ifeq ($(SYM_COMPILE_IN_CP),TRUE)
USER_TM_OBJS      += src/CutPool/user_cp.o 
endif

USER_MASTER_OBJS   = src/Master/user_main.o \
	src/Master/user_master.o $(ADDSURCES)

ifeq ($(SYM_COMPILE_IN_TM),TRUE)
USER_MASTER_OBJS  += $(USER_TM_OBJS)
endif

USER_CG_OBJS       = src/CutGen/user_cg.o

USER_CP_OBJS       = src/CutPool/user_cp.o

USER_DG_OBJS       = src/DrawGraph/user_dg.o

ALL_OBJS =   src/Master/user_main.o src/Master/user_master.o\
	src/LP/user_lp_branch.o src/LP/user_lp.o \
	$(USER_CG_OBJS) $(USER_CP_OBJS) $(USER_DG_OBJS)

PUSER_MASTER_OBJS = $(addsuffix appl, $(USER_MASTER_OBJS))
PUSER_LP_OBJS = $(addsuffix appl, $(USER_LP_OBJS))
PUSER_CG_OBJS = $(addsuffix appl, $(USER_CG_OBJS))
PUSER_CP_OBJS = $(addsuffix appl, $(USER_CP_OBJS))
PUSER_TM_OBJS = $(addsuffix appl, $(USER_TM_OBJS))

PALL_OBJS = $(PUSER_MASTER_OBJS) $(PUSER_LP_OBJS) $(PUSER_CP_OBJS) \
	$(PUSER_TM_OBJS)

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -g -pipe

# additional C++ Compiler options for linking
CXXLINKFLAGS = 

#Include files
INCL =  -Iinclude -I`$(CYGPATH_W) $(COININCDIR)` \
	-I`$(CYGPATH_W) $(SYMINCDIR)`

# Libraries necessary to link
LIBS = -L$(COINLIBDIR) -L$(SYMLIBDIR)

ifeq ($(USE_CGL_CUTS),TRUE)
LIBS += -lCgl `cat /home/menal/Sym-Dev/SYMPHONY/../Cgl/cgl_addlibs.txt`
ADDINCFLAGS += -DUSE_CGL_CUTS
endif

ifeq ($(CLP_LP_SOLVER),TRUE)
LIBS += -lOsiClp -lClp `cat /clp_addlibs.txt`
ADDINCFLAGS += -D__OSI_CLP__
endif

ifeq ($(GLPK_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(GLPKINCDIR)`
LIBS += -lOsiGlpk #-lglpk
ADDINCFLAGS += -D__OSI_GLPK__
else
ifeq ($(USE_GMPL), TRUE)
INCL += -I`$(CYGPATH_W) $(GLPKINCDIR)`
ADDINCFLAGS += -DUSE_GLPMPL
endif
endif

ifeq ($(CPLEX_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(CPXINCDIR)`
LIBS += -l OsiCpx -lcplex
ADDINCFLAGS += -D__OSI_CPLEX__
endif

ifeq ($(OSL_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(OSLINCDIR)`
LIBS += -l OsiOsl -losl
ADDINCFLAGS += -D__OSI_OSL__
endif

ifeq ($(SOPLEX_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(CPXINCDIR)`
LIBS += -l
ADDINCFLAGS += -D__OSI_SPX__
endif

ifeq ($(XPRESS_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(CPXINCDIR)`
LIBS += -l
ADDINCFLAGS += -D__OSI_SPX__
endif

ADDLIBS += /home/menal/GLPK/glpk-4.11/src/libglpk.a -lreadline -lncurses  -lm

LIBS += -lOsi -lCoinUtils `cat /home/menal/Sym-Dev/SYMPHONY/../Osi/osi_addlibs.txt` \
	`cat /home/menal/Sym-Dev/SYMPHONY/../CoinUtils/coinutils_addlibs.txt` $(ADDLIBS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

WHATTOMAKE = $(EXE)

ifeq ($(SYM_PARALLEL),TRUE)
WHATTOMAKE += $(MASTERNAME)$(MASTEREXT)
endif

ifneq ($(SYM_COMPILE_IN_LP),TRUE)
WHATTOMAKE  += $(MASTERNAME)_lp$(LPEXT)
endif

ifneq ($(SYM_COMPILE_IN_CP),TRUE)
WHATTOMAKE += $(MASTERNAME)_cp
endif

ifneq ($(SYM_COMPILE_IN_CG),TRUE)
WHATTOMAKE += $(MASTERNAME)_cg
endif

ifneq ($(SYM_COMPILE_IN_TM),TRUE)
WHATTOMAKE += $(MASTERNAME)_tm$(TMEXT)
endif

ADDFLAGS += $(ADDINCFLAGS) -DHAS_RANDOM -DHAS_SRANDOM -D__NONE__ -D__LINUX -DCOMPILE_IN_CG -DCOMPILE_IN_CP -DCOMPILE_IN_LP -DCOMPILE_IN_TM
ADDPFLAGS += $(ADDINCFLAGS) -DHAS_RANDOM -DHAS_SRANDOM -D__NONE__ -D__LINUX -DCOMPILE_IN_CG -DCOMPILE_IN_CP -DCOMPILE_IN_LP -DCOMPILE_IN_TM

INCL += -I`$(CYGPATH_W) $(PVM_ROOT)/include`
SYMLIBS = -lSymAppl $(LIBS) 
PSYMLIBS = -L$(PVM_ROOT)/lib/LINUX -lgpvm3 -lpvm3 $(LIBS) 

all: $(APPL_PARALLEL) $(WHATTOMAKE)

.c.o:
	$(CXX) $(CXXFLAGS) $(INCL) $(ADDFLAGS) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.c.obj:
	$(CXX) $(CXXFLAGS) $(INCL) $(ADDFLAGS) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`

.c.oappl:
	$(CXX) $(CXXFLAGS) $(INCL) $(ADDPFLAGS) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.c.objappl:
	$(CXX) $(CXXFLAGS) $(INCL) $(ADDPFLAGS) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`

$(EXE): $(ALL_OBJS)
	bla=;\
	for file in $(ALL_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  $(SYMLIBS)

$(MASTERNAME)$(MASTEREXT) : $(PUSER_MASTER_OBJS)
	bla=;\
	for file in $(PUSER_MASTER_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  \
	-lSymAppl$(MASTEREXT) $(PSYMLIBS)

$(MASTERNAME)_lp$(LPEXT) : $(PUSER_LP_OBJS)
	bla=;\
	for file in $(PUSER_LP_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  \
	-lSymAppl_lp$(LPEXT) $(PSYMLIBS)

$(MASTERNAME)_cp : $(PUSER_CP_OBJS)
	bla=;\
	for file in $(PUSER_CP_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  \
	-lSymAppl_cp $(PSYMLIBS)

$(MASTERNAME)_cg : $(PUSER_CG_OBJS)
	bla=;\
	for file in $(PUSER_CG_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  $(PSYMLIBS) \
	-lSymAppl_cg

$(MASTERNAME)_tm$(TMEXT) : $(PUSER_TM_OBJS)
	bla=;\
	for file in $(PUSER_TM_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla  \
	-lSymAppl_tm$(TMEXT) $(PSYMLIBS)

.SUFFIXES: .cpp .c .o .obj .oappl .objappl

clean:
	rm -rf $(WHATTOMAKE) $(ALL_OBJS) $(PALL_OBJS)
