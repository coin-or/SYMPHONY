##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the          #
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for  #
# the Vehicle Routing Problem and the Traveling Salesman Problem.            #
#                                                                            #
# (c) Copyright 2000, 2001, 2002  Ted Ralphs. All Rights Reserved.           #
#                                                                            #
# This application was developed by Ted Ralphs (tkralphs@lehigh.edu)         #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
##############################################################################
##############################################################################

##############################################################################
# First, we have to define some variables that will be used
##############################################################################

##############################################################################
# SYMPHONY's home directory
##############################################################################

SYMPHONYROOT = ${HOME}/SYMPHONY

##############################################################################
# The home directory of the user files (must be an absolute path)
##############################################################################

USERROOT = ${SYMPHONYROOT}/Template

##############################################################################
# The customized name of the master program (if desired)
##############################################################################

MASTERNAME = mip

##############################################################################
# These are only needed if parallel computation is used. Please see the
# explanation of these variables in SYMPHONY's Makefile.
##############################################################################

COMPILE_IN_TM = FALSE
COMPILE_IN_LP = FALSE
COMPILE_IN_CG = TRUE
COMPILE_IN_CP = FALSE

##############################################################################
# This section is for listing source file names
##############################################################################

USER_LP_SRC       = user_lp_branch.c user_lp.c 
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += user_cg.c
endif

USER_TM_SRC       =
ifeq ($(COMPILE_IN_LP),TRUE)
USER_TM_SRC      += $(USER_LP_SRC)
endif
ifeq ($(COMPILE_IN_CP),TRUE)
USER_TM_SRC      += user_cp.c 
endif

USER_MASTER_SRC   = user_master.c

ifeq ($(COMPILE_IN_TM),TRUE)
USER_MASTER_SRC  += $(USER_TM_SRC)
endif

USER_CG_SRC       = user_cg.c

USER_CP_SRC       = user_cp.c

USER_DG_SRC       = user_dg.c

ALL_USER_SRC      = user_lp_branch.c user_lp.c
ALL_USER_SRC     += user_master.c user_cg.c user_cp.c user_dg.c

USER_SRCDIR      += 

##############################################################################
# Include the SYMPHONY makefile after all the source file names (and before
# any used-defined targets!)
##############################################################################

include $(SYMPHONYROOT)/Makefile
