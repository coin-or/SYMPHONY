##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the          #
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for  #
# the Vehicle Routing Problem and the Traveling Salesman Problem.            #
#                                                                            #
# (c) Copyright 2000-2003  Ted Ralphs. All Rights Reserved.                  #
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

USERROOT = ${SYMPHONYROOT}/CNRP

##############################################################################
# The customized name of the master program (if desired)
##############################################################################

MASTERNAME = cnrp

##############################################################################
# These are only needed if parallel computation is used. Please see the
# explanation of these variables in SYMPHONY's Makefile.
##############################################################################

COMPILE_IN_TM = TRUE
COMPILE_IN_LP = TRUE
COMPILE_IN_CG = TRUE
COMPILE_IN_CP = TRUE

##############################################################################
# Other variables to be set by the user
##############################################################################

USER_MAIN = TRUE
DO_CONCORDE_CUTS = FALSE
DIRECTED_X_VARS = TRUE
ADD_FLOW_VARS = TRUE
ADD_CAP_CUTS = TRUE
ADD_X_CUTS = TRUE

ifeq ($(DO_CONCORDE_CUTS),TRUE)
USER_BB_DEFINES += -DDO_TSP_CUTS
CONCORDE_ROOT = $(USERROOT)/Tsp
USERLIBS = -lconcorde
USERLDFLAGS = -L${HOME}/lib
endif

ifeq ($(ADD_FLOW_VARS),TRUE)
USER_BB_DEFINES += -DADD_FLOW_VARS
ifeq ($(ADD_CAP_CUTS),TRUE)
USER_BB_DEFINES += -DADD_CAP_CUTS
endif
endif

ifeq ($(DIRECTED_X_VARS),TRUE)
USER_BB_DEFINES += -DDIRECTED_X_VARS
ifeq ($(ADD_X_CUTS),TRUE)
USER_BB_DEFINES += -DADD_X_CUTS
endif
endif

##############################################################################
# This section is for listing CNRP source file names
##############################################################################

MACROS_SRC        = cnrp_macros.c
COST_SRC          = compute_cost.c
NET_SRC           = network.c
DG_FUNC_SRC       = cnrp_dg_functions.c

USER_MASTER_SRC   = cnrp_master.c cnrp_io.c cnrp_master_functions.c 
USER_MASTER_SRC  += small_graph.c
ifeq ($(USER_MAIN),TRUE)
USER_MASTER_SRC  += cnrp_main.c
endif
USER_MASTER_SRC  += $(COST_SRC)
ifeq ($(COMPILE_IN_TM),TRUE)
USER_MASTER_SRC  += $(USER_TM_SRC)
ifeq ($(COMPILE_IN_LP),FALSE)
USER_MASTER_SRC  += $(MACROS_SRC) $(DG_FUNC_SRC)
endif
else
USER_MASTER_SRC  += $(MACROS_SRC) $(DG_FUNC_SRC)
endif

ifeq ($(COMPILE_IN_LP),TRUE)
USER_TM_SRC       = $(USER_LP_SRC) 
endif
ifeq ($(COMPILE_IN_CP),TRUE)
USER_TM_SRC      += cnrp_cp.c 
endif

USER_LP_SRC       = cnrp_lp_branch.c cnrp_lp.c 
USER_LP_SRC      += $(MACROS_SRC) $(NET_SRC) $(DG_FUNC_SRC)
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += cnrp_cg.c biconnected.c shrink.c 
endif

USER_CG_SRC       = cnrp_cg.c biconnected.c 
USER_CG_SRC      += compute_cost.c shrink.c $(DG_FUNC_SRC)
USER_CG_SRC      += $(MACROS_SRC) $(NET_SRC)

USER_CP_SRC       = cnrp_cp.c $(MACROS_SRC)

USER_DG_SRC       = cnrp_dg.c cnrp_dg_network.c $(MACROS_SRC)

ALL_USER_SRC      = $(MACROS_SRC) $(COST_SRC) $(NET_SRC) $(DG_FUNC_SRC)
ALL_USER_SRC     += vrp_lp_branch.c vrp_lp.c
ALL_USER_SRC     += vrp_master.c vrp_io.c vrp_master_functions.c
ALL_USER_SRC     += start_heurs.c lower_bound.c route_heur.c exchange_heur.c
ALL_USER_SRC     += cluster_heur.c receive_rout.c small_graph.c
ALL_USER_SRC     += vrp_cg.c biconnected.c
ALL_USER_SRC     += compute_cost.c decomp_user.c vrp_dg.c heur_routines.c
ALL_USER_SRC     += vrp_cp.c vrp_dg.c vrp_dg_network.c 

USER_INCDIR       =

##############################################################################
# This section is for CONCORDE source file names
##############################################################################

##############################################################################
# This solver can now solve TSPs using separation routines from CONCORDE, the
# TSP solver of Applegate, Bixby, Chvatal, and Cook. To enable this option,
# set the variables DO_CONCORDE_CUTS to TRUE. Download the source code from
# http://www.math.princeton.edu/tsp/concorde.html and compile it as a library.
# Then make sure to add the path to the library file libconcorde.a to the
# variable USERLIBPATHS, as below (currently, it is set to ${HOME}/lib).
# Alternatively, this makefile can also be used to make the library. Just set
# the variable CONCORDE_ROOT below to the CONCORDE source directory and type
# "make concorde.a" in the SYMPHONY root directory (not in this directory).
##############################################################################

DO_CONCORDE_CUTS = FALSE

ifeq ($(DO_CONCORDE_CUTS),TRUE)

USER_BB_DEFINES += -DDO_TSP_CUTS
USER_BB_DEFINES += -DCC_ONE_ENV
CONCORDE_ROOT = $(HOME)/Concorde/
USERLIBS = -lconcorde -lcplex
USERLIBPATHS  = ${HOME}/lib
USER_INCDIR += -I$(CONCORDE_ROOT)/INCLUDE

ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += tsp.c # $(TSP_SRC)
endif

USER_CG_SRC      += tsp.c # $(TSP_SRC)

TSP_TSP_SRC       = cliqwork.c combcliq.c skeleton.c generate.c ex_price.c 
TSP_TSP_SRC      += growcomb.c teething.c bcontrol.c
TSP_TSP_SRC      += blkcomb.c blossom.c tighten.c control.c 
TSP_TSP_SRC      += prclique.c cutpool.c cutcall.c tsp_call.c 
TSP_TSP_SRC      += tsp_lp.c prob_io.c branch.c 
TSP_TSP_SRC      += cliqhash.c ddecker.c qsparse.c xtour.c
TSP_UTIL_SRC      = edgelen.c genhash.c edgemap.c getdata.c sortrus.c
TSP_UTIL_SRC     += allocrus.c edgeutil.c priority.c urandom.c bgetopt.c 
TSP_UTIL_SRC     += eunion.c profrus.c util.c dheaps_i.c fastread.c safe_io.c
TSP_UTIL_SRC     += zeit.c 
TSP_CUT_SRC       = cut_st.c shrink.c gomoryhu.c mincut.c connect.c
TSP_CUT_SRC      += segments.c
TSP_BIGGUY_SRC    = bigguy.c
TSP_COMB_SRC     += block.c dngreedy.c
TSP_LP_SRC       += lpcplex5.c
TSP_PQ_SRC       += consec1.c necklace.c cuttree.c pq.c
TSP_KD_SRC        = kdnear.c kdtwoopt.c kdbuild.c kdspan.c
TSP_EDGEGEN_SRC   = delaunay.c xnear.c edgegen.c mlinkern.c
TSP_FMATCH_SRC    = fmatch.c
TSP_LOCALCUT_SRC  = lift.c separate.c first.c tsporacl.c localcut.c
TSP_LOCALCUT_SRC += intmat.c chunks.c
TSP_VERIFY_SRC    = verify.c
TSP_HELDKARP_SRC  = heldkarp.c
TSP_LINKERN_SRC   = linkern.c flip_two.c
TSP_TINY_SRC      = benttsp.c bnbtsp.c tinytsp.c bnbmsp.c

TSP_SRC           = $(TSP_TSP_SRC) $(TSP_UTIL_SRC) $(TSP_CUT_SRC) 
TSP_SRC          += $(TSP_BIGGUY_SRC) $(TSP_COMB_SRC) $(TSP_LP_SRC)
TSP_SRC          += $(TSP_PQ_SRC) $(TSP_KD_SRC) $(TSP_EDGEGEN_SRC)
TSP_SRC          += $(TSP_FMATCH_SRC) $(TSP_LOCALCUT_SRC)
TSP_SRC          += $(TSP_VERIFY_SRC) $(TSP_HELDKARP_SRC)
TSP_SRC          += $(TSP_LINKERN_SRC) $(TSP_TINY_SRC)

ALL_USER_SRC     += tsp.c # $(TSP_SRC)

USER_SRCDIR    += $(CONCORDE_ROOT)/TSP $(CONCORDE_ROOT)/UTIL 
USER_SRCDIR    += $(CONCORDE_ROOT)/CUT $(CONCORDE_ROOT)/BIGGUY
USER_SRCDIR    += $(CONCORDE_ROOT)/COMBS $(CONCORDE_ROOT)/LP
USER_SRCDIR    += $(CONCORDE_ROOT)/PQ $(CONCORDE_ROOT)/KDTREE
USER_SRCDIR    += $(CONCORDE_ROOT)/EDGEGEN $(CONCORDE_ROOT)/FMATCH
USER_SRCDIR    += $(CONCORDE_ROOT)/LOCALCUT $(CONCORDE_ROOT)/VERIFY
USER_SRCDIR    += $(CONCORDE_ROOT)/HELDKARP $(CONCORDE_ROOT)/LINKERN
USER_SRCDIR    += $(CONCORDE_ROOT)/TINY

endif

TSP_OBJS = $(addprefix $(USER_OBJDIR)/,$(notdir $(TSP_SRC:.c=.o)))

##############################################################################
# Include the SYMPHONY makefile after all the source file names (and before
# any used-defined targets!)
##############################################################################

include $(SYMPHONYROOT)/Makefile

##############################################################################
# Targets for CONCORDE
##############################################################################

concorde:
	mv ${USERROOT}/CutGen/shrink.c ${USERROOT}/CutGen/temp.c
	rm $(USER_OBJDIR)/shrink.o

concorde.a : concorde $(TSP_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	$(AR) $(HOME)/lib/libconcorde.a $(TSP_OBJS)
	$(RANLIB) $(HOME)/lib/libconcorde.a
	mv ${USERROOT}/CutGen/temp.c ${USERROOT}/CutGen/shrink.c
	rm $(USER_OBJDIR)/shrink.o
	@echo ""

