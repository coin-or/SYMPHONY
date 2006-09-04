# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

# $Id$

MASTERNAME = mckp

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

USE_CGL_CUTS = TRUE
CLP_LP_SOLVER = #TRUE
GLPK_LP_SOLVER = TRUE
CPLEX_LP_SOLVER = #TRUE
OSL_LP_SOLVER = #TRUE
SOPLEX_LP_SOLVER = #TRUE
XPRESS_LP_SOLVER = #TRUE
USE_GMPL = TRUE

USER_MASTER_OBJS   = mckp_main.o $(ADDSOURCES)

ALL_OBJS =   $(USER_MASTER_OBJS)

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
	`cat /home/menal/Sym-Dev/SYMPHONY/../CoinUtils/coinutils_addlibs.txt` -lSym \
	$(ADDLIBS)

ADDFLAGS += $(ADDINCFLAGS) -DHAS_RANDOM -DHAS_SRANDOM -D__NONE__ -D__LINUX -DCOMPILE_IN_CG -DCOMPILE_IN_CP -DCOMPILE_IN_LP -DCOMPILE_IN_TM

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

WHATTOMAKE = $(EXE)



all: $(WHATTOMAKE)

.c.o:
	$(CXX) $(CXXFLAGS) $(INCL) $(ADDFLAGS) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.c.obj:
	$(CXX) $(CXXFLAGS) $(INCL) $(ADDFLAGS) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`


$(EXE): $(ALL_OBJS)
	bla=;\
	for file in $(ALL_OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; \
	done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS)

.SUFFIXES: .cpp .c .o .obj

clean:
	rm -rf $(WHATTOMAKE) $(USER_MASTER_OBJS)
