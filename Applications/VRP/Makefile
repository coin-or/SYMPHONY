##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the          #
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for  #
# the Vehicle Routing Problem and the Traveling Salesman Problem.            #
#                                                                            #
# (c) Copyright 2000-2005  Ted Ralphs. All Rights Reserved.                  #
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

SYMPHONYROOT = ..

##############################################################################
# The home directory of the user files (must be an absolute path)
##############################################################################

USERROOT = .

##############################################################################
# The customized name of the master program (if desired)
##############################################################################

MASTERNAME = vrp

##############################################################################
# These are only needed if parallel computation is used. Please see the
# explanation of these variables in SYMPHONY's Makefile.
##############################################################################

COMPILE_IN_TM = TRUE
COMPILE_IN_LP = TRUE
COMPILE_IN_CG = TRUE
COMPILE_IN_CP = TRUE

##############################################################################
# This section is for listing VRP source file names
##############################################################################

MACROS_SRC        = vrp_macros.c
COST_SRC          = compute_cost.c
NET_SRC           = network.c
DG_FUNC_SRC       = vrp_dg_functions.c 

USER_MASTER_SRC   = vrp_main.c vrp_master.c vrp_io.c vrp_master_functions.c 
USER_MASTER_SRC  += small_graph.c
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
USER_TM_SRC      += vrp_cp.c 
endif

USER_LP_SRC       = vrp_lp_branch.c vrp_lp.c 
USER_LP_SRC      += $(MACROS_SRC) $(NET_SRC) $(DG_FUNC_SRC)
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += vrp_cg.c biconnected.c shrink.c 
endif

USER_CG_SRC       = vrp_cg.c biconnected.c
USER_CG_SRC      += compute_cost.c shrink.c $(DG_FUNC_SRC)
USER_CG_SRC      += $(MACROS_SRC) $(NET_SRC)

USER_CP_SRC       = vrp_cp.c $(MACROS_SRC)

USER_DG_SRC       = vrp_dg.c vrp_dg_network.c $(MACROS_SRC)

ALL_USER_SRC      = $(MACROS_SRC) $(COST_SRC) $(NET_SRC) $(DG_FUNC_SRC)
ALL_USER_SRC     += vrp_main.c vrp_lp_branch.c vrp_lp.c
ALL_USER_SRC     += vrp_master.c vrp_io.c vrp_master_functions.c
ALL_USER_SRC     += vrp_cg.c biconnected.c
ALL_USER_SRC     += compute_cost.c decomp_user.c vrp_dg.c heur_routines.c
ALL_USER_SRC     += vrp_cp.c vrp_dg.c vrp_dg_network.c 
#__BEGIN_EXPERIMENTAL_SECTION__#
ALL_USER_SRC     += start_heurs.c lower_bound.c route_heur.c exchange_heur.c
ALL_USER_SRC     += cluster_heur.c receive_rout.c small_graph.c
#__END_EXPERIMENTAL_SECTION__#

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
# "make concord.a" in the SYMPHONY root directory (not in this directory).
##############################################################################

DO_CONCORDE_CUTS = FALSE

ifeq ($(DO_CONCORDE_CUTS),TRUE)

USER_BB_DEFINES += -DDO_TSP_CUTS
USER_BB_DEFINES += -DCC_ONE_ENV
CONCORDE_ROOT = $(HOME)/Concorde/
USERLIBS = -lconcorde -lcplex
USERLIBPATHS  = ${HOME}/lib
USER_INCDIR += -I$(CONCORDE_ROOT)/INCLUDE

ifeq ($(APPL_COMPILE_IN_CG),TRUE)
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

#__BEGIN_EXPERIMENTAL_SECTION__#
##############################################################################
# This section is for the VRP heuristic source file names
##############################################################################

USER_MASTER_SRC  += start_heurs.c lower_bound.c route_heur.c
USER_MASTER_SRC  += exchange_heur.c cluster_heur.c
USER_MASTER_SRC  += receive_rout.c collect_solutions.c heur_routines.c
USER_MASTER_SRC  += proccomm.c timemeas.c

HEURS_SRC     = $(MACROS_SRC) sweep.c savings.c savings_rout.c savings2.c 
HEURS_SRC    += savings_rout2.c binomial.c heur_routines.c nearest_ins.c 
HEURS_SRC    += ins_routines.c farthest_ins.c farnear_ins.c compute_cost.c
HEURS_SRC    += near_cluster.c ins_routines2.c tsp_ni.c tsp_ins_rout.c 
HEURS_SRC    += s_path.c tsp_fi.c tsp_fini.c exchange.c exchange2.c mst.c 
HEURS_SRC    += mst_ins_rout.c savings3.c savings_rout3.c qsort.c
HEURS_SRC    += $(TIME_SRC) $(PVMERR_SRC) heur_main.c

HEURS_OBJS    = $(addprefix $(USER_OBJDIR)/,$(notdir $(HEURS_SRC:.c=.o)))
HEURS_DEP     = $(addprefix $(USER_DEPDIR)/,$(notdir $(HEURS_SRC:.c=.d)))

USER_SRCDIR   = $(USERROOT)/Heuristics
USER_SRCDIR  += $(USERROOT)/Master/Heuristics

USER_INCDIR  += -I$(USERROOT)/include/heurs

ALL_USER_SRC += $(ALL_HEURS_SRC)

##############################################################################
# This section is for DECOMP source file names
##############################################################################

COMPILE_OUR_DECOMP = FALSE
FARKAS_LIFTING_BOUND = TSP
FARKAS_LIFTING_METHOD = NONE
ifeq ($(COMPILE_OUR_DECOMP),TRUE)
USER_BB_DEFINES += -DCOMPILE_OUR_DECOMP
endif
ifeq ($(FARKAS_LIFTING_BOUND),TSP)
USER_BB_DEFINES += -DEXACT_LIFTING
endif
ifeq ($(FARKAS_LIFTING_METHOD),SEQUENTIAL)
USER_BB_DEFINES += -DFARKAS_LIFT_SEQUENTIAL
endif
ifeq ($(FARKAS_LIFTING_METHOD),NAIVE)
USER_BB_DEFINES += -DFARKAS_LIFT_NAIVE
endif
ifeq ($(FARKAS_LIFTING_METHOD),NONE)
USER_BB_DEFINES += -DNO_LIFTING
endif

ifeq ($(APPL_COMPILE_IN_CG),TRUE)
ifeq ($(DECOMP),TRUE)
ifeq ($(COMPILE_OUR_DECOMP),TRUE)
USER_LP_SRC      += my_decomp_user.c decomp_lower_bound.c
else
USER_LP_SRC      += decomp_user.c ind_sort.c
endif
endif
endif

ifeq ($(DECOMP),TRUE)
ifeq ($(COMPILE_OUR_DECOMP),TRUE)
USER_CG_SRC      += my_decomp_user.c decomp_lower_bound.c
else
USER_CG_SRC      += decomp_user.c ind_sort.c
endif
endif

ifeq ($(COMPILE_OUR_DECOP),TRUE)
USER_INCDIR      += -I$(USER_ROOT)/include/decomp
endif
#__END_EXPERIMENTAL_SECTION__#

##############################################################################
# This is the flag for the SYMPHONY Makefile to determine whether
# an application is to be used. You shouldn't need to modify this
# variable. 
##############################################################################

USE_SYM_APPL = TRUE

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

#__BEGIN_EXPERIMENTAL_SECTION__#
##############################################################################
# Targets for the heuristics
##############################################################################

clean_heurs:
	rm -f $(addprefix $(USER_OBJDIR)/, $(notdir $(HEURS_OBJS)))

vrp_heuristics: $(BINDIR)/vrp_heuristics
	true

$(BINDIR)/vrp_heuristics: $(HEURS_DEP) $(HEURS_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(HEURS_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

##############################################################################
# Targets for decomp
##############################################################################

clean_decomp:
	rm -f $(USER_OBJDIR)/my_decomp_user.o 
	$(USER_OBJDIR)/decomp_lower_bound.o

#__END_EXPERIMENTAL_SECTION__#


