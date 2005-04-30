##############################################################################
##############################################################################
#                                                                            #
# This file is part of the SYMPHONY Branch, Cut, and Price Library.          #
#                                                                            #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and     #
# Laci Ladanyi (ladanyi@us.ibm.com).                                         #
#                                                                            #
# (c) Copyright 2000-2004  Ted Ralphs. All Rights Reserved.                  #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
##############################################################################
##############################################################################

CC = /usr/local/bin/g++
OPT = -g

SYMPHONYROOT = ..
COINROOT = ../..

LP_SOLVER = OSI
OSI_INTERFACE = CLP

LPINCDIR = $(SYMPHONYROOT)/include
LPLIBPATHS = $(SYMPHONYROOT)/lib
LPLIB = -lsym

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
       LPLIB += -lCoin -lOsi -lCgl
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
endif

INCDIR        = $(addprefix -I,${LPINCDIR})
LDFLAGS       = $(addprefix -L,${LPLIBPATHS})

EXAMPLES = milp bicriteria sensitivity warm_start1 warm_start2 warm_start3

examples: $(EXAMPLES)
	true



milp : milp.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB)

bicriteria : bicriteria.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB)

sensitivity : sensitivity.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB)

warm_start1 : warm_start1.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB)

warm_start2 : warm_start2.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB)

warm_start3 : warm_start3.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB)


clean : 

	rm -rf $(EXAMPLES)

.SILENT:
