##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the          #
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for  #
# multi-criteria knapsack.                                                   #
#                                                                            #
# (c) Copyright 2005 Lehigh University. All Rights Reserved.                 #
#                                                                            #
# This application was written by Ted Ralphs (tkralphs@lehigh.edu).          #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
##############################################################################
##############################################################################

CC = g++
OPT = -O
ARCH = LINUX

SYMPHONYROOT = ../..
COINROOT = ../../..
MCKPROOT = .

LP_SOLVER = OSI
OSI_INTERFACE = CLP

LPINCDIR = $(SYMPHONYROOT)/include
LPLIB = -lsym

ifeq ($(LP_SOLVER),CPLEX)
       LPINCDIR += /usr/local/include/
       LPLIBPATHS = /usr/local/lib/ 
       LPLIBPATHS += $(SYMPHONYROOT)/lib/${ARCH}/$(LP_SOLVER)
       LPLIB += -lcplex -lpthread
endif

ifeq ($(LP_SOLVER),OSL)
       LPINCDIR +=  $(HOME)/OSL/osllib/include/
       LPLIBPATHS = $(HOME)/OSL/osllib/lib
       LPLIBPATHS += $(SYMPHONYROOT)/lib/${ARCH}/$(LP_SOLVER)
       LPLIB += -losl
endif

ifeq ($(LP_SOLVER),OSI)
       LPINCDIR += $(COINROOT)/include
       LPLIBPATHS += $(COINROOT)/lib
       LPLIBPATHS += $(SYMPHONYROOT)/lib/${ARCH}/$(LP_SOLVER)_$(OSI_INTERFACE)
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
LDFLAGS      += -DHAS_RANDOM -DHAS_SRANDOM

mckp : mckp_main.c

	$(CC) $(OPT) $(LDFLAGS) $(INCDIR) $< -o $@ $(LPLIB) $(LPLIB)

clean : 

	rm mckp
