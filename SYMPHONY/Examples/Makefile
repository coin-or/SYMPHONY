##############################################################################
##############################################################################
#                                                                            #
# This file is part of the SYMPHONY MILP Solver Framework.                   #
#                                                                            #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and     #
# Laci Ladanyi (ladanyi@us.ibm.com).                                         #
#                                                                            #
# (c) Copyright 2005-2006 Lehigh University. All Rights Reserved.            #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
# The author of this file is Menal Guzelsoy (megb@lehigh.edu)                #
#                                                                            #
##############################################################################
##############################################################################

CC = g++
OPT = -O #-fPIC
SHLINKPREFIX = -Wl,-rpath,

SYMPHONYROOT = ..
COINROOT = ../..


LP_SOLVER = OSI
OSI_INTERFACE = CLP

##############################################################################
# If you wish to use SYMPHONY through the SYMPHONY OSI interface, 
# set USE_OSI_INTERFACE to TRUE below. This will have the examples 
# to call the OSI SYMPHONY interface routines. 
##############################################################################
USE_OSI_INTERFACE = FALSE

LPINCDIR = $(COINROOT)/include
LPLIBPATHS = $(COINROOT)/lib
LPLIB = -lSym
SYMDEFS =

ifeq ($(LP_SOLVER),CPLEX)
       LPINCDIR += /usr/local/include/
       LPLIBPATHS += /usr/local/lib/ 
       LPLIB += -lcplex -lpthread
endif

ifeq ($(LP_SOLVER),OSL)
       LPINCDIR +=  $(HOME)/OSL/osllib/include/
       LPLIBPATHS += $(HOME)/OSL/osllib/lib
       LPLIB += -losl
endif

ifeq ($(LP_SOLVER),OSI)
       LPINCDIR += $(COINROOT)/include
       LPLIBPATHS += $(COINROOT)/lib
       LPLIB += -lCgl
ifeq ($(USE_OSI_INTERFACE),TRUE)
       LPLIB += -lOsiSym
       SYMDEFS += -DUSE_OSI_INTERFACE    
endif
ifeq ($(OSI_INTERFACE),CPLEX)
       LPINCDIR += /usr/local/include/
       LPLIBPATHS += /usr/local/lib/
       LPLIB += -lOsiCpx -lcplex -lpthread
endif
ifeq ($(OSI_INTERFACE),OSL)
       LPINCDIR += $(HOME)/OSL/osllib/include
       LPLIBPATHS += $(HOME)/OSL/osllib/lib
       LPLIB += -lOsiOsl -losl
endif
ifeq ($(OSI_INTERFACE),CLP)
       LPINCDIR += 
       LPLIBPATHS +=
       LPLIB += -lOsiClp -lClp
endif
ifeq ($(OSI_INTERFACE),GLPK)
       LPINCDIR += ${HOME}/GLPK/glpk-4.3/include
       LPLIBPATHS += ${HOME}/GLPK/glpk-4.3/src
       LPLIB += -lOsiGlpk -lglpk
endif
       LPLIB += -lOsi -lCoinUtils
endif

INCDIR        = $(addprefix -I,${LPINCDIR})
LDFLAGS       = $(addprefix -L,${LPLIBPATHS})
LDFLAGS      += $(SYMDEFS)

LINKFLAGS += $(addprefix ${SHLINKPREFIX},${LPLIBPATHS})

EXAMPLES = milp milp2 bicriteria sensitivity warm_start1 warm_start2 warm_start3

examples: $(EXAMPLES)
	true

milp : milp.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $(LINKFLAGS) $< -o $@ $(LPLIB)

milp2 : milp2.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $(LINKFLAGS) $< -o $@ $(LPLIB)

bicriteria : bicriteria.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR)  $(LINKFLAGS) $< -o $@ $(LPLIB)

sensitivity : sensitivity.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $(LINKFLAGS) $< -o $@ $(LPLIB)

warm_start1 : warm_start1.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR)  $(LINKFLAGS) $< -o $@ $(LPLIB)

warm_start2 : warm_start2.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR)  $(LINKFLAGS) $< -o $@ $(LPLIB)

warm_start3 : warm_start3.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR)  $(LINKFLAGS) $< -o $@ $(LPLIB)


clean : 

	rm -rf $(EXAMPLES)

.SILENT:
