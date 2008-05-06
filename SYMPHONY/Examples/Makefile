#===========================================================================#
#                                                                           #
# This file is part of the SYMPHONY MILP Solver Framework.                  #
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

# $Id$

EXAMPLES = milp milp2 bicriteria sensitivity \
	warm_start1 warm_start2 warm_start3

##############################################################################
# If you wish to use SYMPHONY through the SYMPHONY OSI interface, 
# set USE_OSI_INTERFACE to TRUE below. This will have the examples 
# to call the OSI SYMPHONY interface routines. 
##############################################################################

USE_OSI_INTERFACE = FALSE

ifeq ($(USE_OSI_INTERFACE),TRUE)
EXAMPLES += milpOsi milpOsi2
endif

# Directory with COIN header files
COININCDIR = /home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/include/coin

# Directory with COIN libraries
COINLIBDIR = /home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/lib

# Directoy with SYMPHONY header files 
#SYMINCDIR = /home/ted/SYMPHONY/stable/5.1/tmp_checkout/SYMPHONY/include

# Directoy with SYMPHONY library
#SYMLIBDIR = /home/ted/SYMPHONY/stable/5.1/tmp_checkout/SYMPHONY/src

# Directory with Application files
SYMEXDIR = /home/ted/SYMPHONY/stable/5.1/tmp_checkout/SYMPHONY/src/../Examples

#Directory for LP solver
GLPKINCDIR = /home/ted/SYMPHONY/stable/5.1/tmp_checkout/SYMPHONY/../ThirdParty/Glpk/glpk/include
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
CLP_LP_SOLVER = TRUE
GLPK_LP_SOLVER = #TRUE
CPLEX_LP_SOLVER = #TRUE
OSL_LP_SOLVER = #TRUE
SOPLEX_LP_SOLVER = #TRUE
XPRESS_LP_SOLVER = #TRUE

VPATH = \
	$(SYMEXDIR)

WHATTOMAKE = $(addsuffix , $(EXAMPLES))

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -fomit-frame-pointer -pipe -DNDEBUG -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion   -Wno-unknown-pragmas -Wno-write-strings -Wno-sign-compare -Wno-cast-qual

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/lib

#Include files
INCL =  -I`$(CYGPATH_W) $(COININCDIR)`

# Libraries necessary to link
LIBS = 

# Additional flags
ADDINCFLAGS =

ifeq ($(USE_CGL_CUTS),TRUE)
LIBS += -lCgl `cat /home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/lib/../share/doc/coin/Cgl/cgl_addlibs.txt`
#ADDINCFLAGS += -DUSE_CGL_CUTS
endif

ifeq ($(CLP_LP_SOLVER),TRUE)
LIBS += -lOsiClp -lClp `cat /home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/lib/../share/doc/coin/Clp/clp_addlibs.txt`
#ADDINCFLAGS += -D__OSI_CLP__
endif

ifeq ($(GLPK_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(GLPKINCDIR)`
LIBS += -lOsiGlpk #-lglpk
#ADDINCFLAGS += -D__OSI_GLPK__
else
ifeq ($(USE_GMPL), TRUE)
INCL += -I`$(CYGPATH_W) $(GLPKINCDIR)`
#ADDINCFLAGS += -DUSE_GLPMPL
endif
endif

ifeq ($(CPLEX_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(CPXINCDIR)`
LIBS += -lOsiCpx #-lcplex
#ADDINCFLAGS += -D__OSI_CPLEX__
endif

ifeq ($(OSL_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(OSLINCDIR)`
LIBS += -lOsiOsl #-losl
#ADDINCFLAGS += -D__OSI_OSL__
endif

ifeq ($(SOPLEX_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(SPXINCDIR)`
LIBS += -lOsiSpx #-l
#ADDINCFLAGS += -D__OSI_SOPLEX__
endif

ifeq ($(XPRESS_LP_SOLVER),TRUE)
INCL += -I`$(CYGPATH_W) $(XPRINCDIR)`
LIBS += -lOsiXpr #-l
#ADDINCFLAGS += -D__OSI_XPRESS__
endif

ADDLIBS += -lm 

LIBS += -lOsi -lCoinUtils `cat /home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/lib/../share/doc/coin/Osi/osi_addlibs.txt` \
	`cat /home/ted/SYMPHONY/stable/5.1/tmp_checkout/build/lib/../share/doc/coin/CoinUtils/coinutils_addlibs.txt` $(ADDLIBS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

SYMLIBS = -L$(COINLIBDIR) #-L$(SYMLIBDIR)

ifeq ($(USE_OSI_INTERFACE),TRUE)
SYMLIBS += -lOsiSym
ADDINCFLAGS += -DUSE_OSI_INTERFACE
endif

SYMLIBS += -lSym $(LIBS) 

all: $(WHATTOMAKE)

%: %.c
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) $(INCL) $(ADDINCFLAGS) $< -o $@ $(SYMLIBS)

clean:
	rm -rf $(WHATTOMAKE)
