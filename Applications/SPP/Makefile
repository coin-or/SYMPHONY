##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the         */
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
# Set Partitioning Problems.                                                */
#                                                                           */
# (c) Copyright 2004 Marta Eso and Ted Ralphs. All Rights Reserved.         */
#                                                                           */
# This application was originally developed by Marta Eso and was modified   */
# Ted Ralphs (tkralphs@lehigh.edu)                                          */
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

SYMPHONYROOT = ..

##############################################################################
# The home directory of the user files (must be an absolute path)
##############################################################################

USERROOT = .

##############################################################################
# The customized name of the master program (if desired)
##############################################################################

MASTERNAME = spp

##############################################################################
# These are only needed if parallel computation is used. Please see the
# explanation of these variables in SYMPHONY's Makefile.
##############################################################################

COMPILE_IN_TM = TRUE
COMPILE_IN_LP = TRUE
COMPILE_IN_CG = TRUE
COMPILE_IN_CP = TRUE

##############################################################################
# This section is for listing SPP source file names
##############################################################################

USER_COMMON_SRC   = spp_common.c

USER_CG_SRC       = spp_cg.c 

USER_CP_SRC       = spp_cp.c 

USER_DG_SRC       = spp_dg.c 

USER_LP_SRC       = spp_lp_branch.c spp_lp.c
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += $(USER_CG_SRC)
endif

USER_TM_SRC       =
ifeq ($(COMPILE_IN_LP),TRUE)
USER_TM_SRC      += $(USER_LP_SRC)
endif
ifeq ($(COMPILE_IN_CP),TRUE)
USER_TM_SRC      += $(USER_CP_SRC)
endif

USER_MASTER_SRC   = spp_main.c spp_master.c spp_master_functions.c \
$(USER_COMMON_SRC)

ifeq ($(COMPILE_IN_TM),TRUE)
USER_MASTER_SRC  += $(USER_TM_SRC)
endif

ALL_USER_SRC      = spp_common.c spp_lp_branch.c spp_lp.c
ALL_USER_SRC     += spp_master.c spp_master_functions.c 
ALL_USER_SRC     += spp_cg.c spp_cp.c spp_dg.c

USER_SRCDIR      += 

##############################################################################
# This is the flag for the SYMPHONY Makefile to determine whether
# a SYMPHONY application is to be used. You shouldn't need to modify this
# variable. 
##############################################################################

USE_SYM_APPL = TRUE

##############################################################################
# Include the SYMPHONY makefile after all the source file names (and before
# any used-defined targets!)
##############################################################################

include $(SYMPHONYROOT)/Makefile
