##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the          #
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for  #
# the Vehicle Routing Problem and the Traveling Salesman Problem.            #
#                                                                            #
# (c) Copyright 2000-2004  Ted Ralphs. All Rights Reserved.                  #
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

USERROOT = ${SYMPHONYROOT}/Vrp

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
USER_MASTER_SRC  += receive_rout.c

SWEEP_SRC         = $(MACROS_SRC) sweep.c qsort.c compute_cost.c 
SWEEP_SRC        += heur_routines.c $(TIME_SRC) $(PVMERR_SRC)

SAVINGS_SRC       = $(MACROS_SRC) savings.c compute_cost.c savings_rout.c 
SAVINGS_SRC      += heur_routines.c $(TIME_SRC) $(PVMERR_SRC)

SAVINGS2_SRC      = $(MACROS_SRC) savings2.c compute_cost.c savings_rout2.c 
SAVINGS2_SRC     += heur_routines.c binomial.c $(TIME_SRC) $(PVMERR_SRC)

SAVINGS3_SRC      = $(MACROS_SRC) savings3.c compute_cost.c savings_rout3.c 
SAVINGS3_SRC     += heur_routines.c binomial.c ins_routines2.c qsort.c 
SAVINGS3_SRC     += $(TIME_SRC) $(PVMERR_SRC)

ROUTE_OPT_SRC     = $(MACROS_SRC) route_opt.c compute_cost.c heur_routines.c 
ROUTE_OPT_SRC    += vrp_routines.c receive_rout.c $(TIME_SRC) $(PVMERR_SRC)

NEAREST_INS_SRC   = $(MACROS_SRC) nearest_ins.c ins_routines.c compute_cost.c 
NEAREST_INS_SRC  += heur_routines.c $(TIME_SRC) $(PVMERR_SRC)

FARTHEST_INS_SRC  = $(MACROS_SRC) farthest_ins.c ins_routines.c compute_cost.c 
FARTHEST_INS_SRC += heur_routines.c $(TIME_SRC) $(PVMERR_SRC)

FARNEAR_INS_SRC   = $(MACROS_SRC) farnear_ins.c ins_routines.c compute_cost.c 
FARNEAR_INS_SRC  += heur_routines.c $(TIME_SRC) $(PVMERR_SRC)

NEAR_CLUSTER_SRC  = $(MACROS_SRC) near_cluster.c ins_routines2.c 
NEAR_CLUSTER_SRC += compute_cost.c heur_routines.c qsort.c $(TIME_SRC) 
NEAR_CLUSTER_SRC += $(PVMERR_SRC)

TSP_NI_SRC        = $(MACROS_SRC) tsp_ni.c tsp_ins_rout.c heur_routines.c 
TSP_NI_SRC       += compute_cost.c s_path.c $(TIME_SRC) $(PVMERR_SRC)

TSP_FI_SRC        = $(MACROS_SRC) tsp_fi.c tsp_ins_rout.c heur_routines.c 
TSP_FI_SRC       += compute_cost.c s_path.c $(TIME_SRC) $(PVMERR_SRC)

TSP_FINI_SRC      = $(MACROS_SRC) tsp_fini.c tsp_ins_rout.c heur_routines.c 
TSP_FINI_SRC     += compute_cost.c s_path.c $(TIME_SRC) $(PVMERR_SRC)

EXCHANGE_SRC      = $(MACROS_SRC) exchange.c heur_routines.c compute_cost.c 
EXCHANGE_SRC     += $(TIME_SRC) $(PVMERR_SRC)

EXCHANGE2_SRC     = $(MACROS_SRC) exchange2.c heur_routines.c compute_cost.c 
EXCHANGE2_SRC    += $(TIME_SRC) $(PVMERR_SRC)

MST_SRC           = $(MACROS_SRC) mst.c compute_cost.c mst_ins_rout.c 
MST_SRC          += heur_routines.c $(TIME_SRC) $(PVMERR_SRC)

ALL_HEURS_SRC     = $(MACROS_SRC) sweep.c savings.c savings_rout.c savings2.c 
ALL_HEURS_SRC    += savings_rout2.c binomial.c route_opt.c nearest_ins.c 
ALL_HEURS_SRC    += ins_routines.c farthest_ins.c farnear_ins.c 
ALL_HEURS_SRC    += near_cluster.c ins_routines2.c tsp_ni.c tsp_ins_rout.c 
ALL_HEURS_SRC    += s_path.c tsp_fi.c tsp_fini.c exchange.c exchange2.c mst.c 
ALL_HEURS_SRC    += mst_ins_rout.c savings3.c savings_rout3.c qsort.c
ALL_HEURS_SRC    += $(TIME_SRC) $(PVMERR_SRC)

HEUR_BINS         = sweep savings route_opt nearest_ins farthest_ins 
HEUR_BINS        += farnear_ins savings2 near_cluster tsp_ni tsp_fi
HEUR_BINS        += tsp_fini exchange exchange2 mst savings3

SWEEP_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(SWEEP_SRC:.c=.o)))
SWEEP_DEP 	  = $(addprefix $(DEPDIR)/,$(SWEEP_SRC:.c=.d))

SAVINGS_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(SAVINGS_SRC:.c=.o)))
SAVINGS_DEP 	  = $(addprefix $(DEPDIR)/,$(SAVINGS_SRC:.c=.d))

SAVINGS2_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(SAVINGS2_SRC:.c=.o)))
SAVINGS2_DEP 	  = $(addprefix $(DEPDIR)/,$(SAVINGS2_SRC:.c=.d))

SAVINGS3_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(SAVINGS3_SRC:.c=.o)))
SAVINGS3_DEP 	  = $(addprefix $(DEPDIR)/,$(SAVINGS3_SRC:.c=.d))

ROUTE_OPT_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(ROUTE_OPT_SRC:.c=.o)))
ROUTE_OPT_DEP 	  = $(addprefix $(DEPDIR)/,$(ROUTE_OPT_SRC:.c=.d))

NEAREST_INS_OBJS  = $(addprefix $(USER_OBJDIR)/,$(notdir $(NEAREST_INS_SRC:.c=.o)))
NEAREST_INS_DEP   = $(addprefix $(DEPDIR)/,$(NEAREST_INS_SRC:.c=.d))

FARTHEST_INS_OBJS = $(addprefix $(USER_OBJDIR)/,$(notdir $(FARTHEST_INS_SRC:.c=.o)))
FARTHEST_INS_DEP  = $(addprefix $(DEPDIR)/,$(FARTHEST_INS_SRC:.c=.d))

FARNEAR_INS_OBJS  = $(addprefix $(USER_OBJDIR)/,$(notdir $(FARNEAR_INS_SRC:.c=.o)))
FARNEAR_INS_DEP   = $(addprefix $(DEPDIR)/,$(FARNEAR_INS_SRC:.c=.d))

NEAR_CLUSTER_OBJS = $(addprefix $(USER_OBJDIR)/,$(notdir $(NEAR_CLUSTER_SRC:.c=.o)))
NEAR_CLUSTER_DEP  = $(addprefix $(DEPDIR)/,$(NEAR_CLUSTER_SRC:.c=.d))

TSP_NI_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(TSP_NI_SRC:.c=.o)))
TSP_NI_DEP 	  = $(addprefix $(DEPDIR)/,$(TSP_NI_SRC:.c=.d))

TSP_FI_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(TSP_FI_SRC:.c=.o)))
TSP_FI_DEP 	  = $(addprefix $(DEPDIR)/,$(TSP_FI_SRC:.c=.d))

TSP_FINI_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(TSP_FINI_SRC:.c=.o)))
TSP_FINI_DEP 	  = $(addprefix $(DEPDIR)/,$(TSP_FINI_SRC:.c=.d))

EXCHANGE_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(EXCHANGE_SRC:.c=.o)))
EXCHANGE_DEP 	  = $(addprefix $(DEPDIR)/,$(EXCHANGE_SRC:.c=.d))

EXCHANGE2_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(EXCHANGE2_SRC:.c=.o)))
EXCHANGE2_DEP 	  = $(addprefix $(DEPDIR)/,$(EXCHANGE2_SRC:.c=.d))

MST_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(MST_SRC:.c=.o)))
MST_DEP 	  = $(addprefix $(DEPDIR)/,$(MST_SRC:.c=.d))

USER_SRCDIR 	= $(USERROOT)/Heuristics
USER_SRCDIR    += $(USERROOT)/Master/Heuristics

USER_INCDIR    += -I$(USERROOT)/include/heurs

ALL_USER_SRC     += $(ALL_HEURS_SRC)

HEUR_OBJS 	  = $(addprefix $(USER_OBJDIR)/,$(notdir $(ALL_HEURS_SRC:.c=.o)))

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

heuristics : 
	$(MAKE) $(HEUR_BINS)

clean_heurs:
	rm -f $(addprefix $(USER_OBJDIR)/, $(notdir $(HEUR_OBJS)))


sweep : $(BINDIR)/sweep
	true

$(BINDIR)/sweep : $(SWEEP_DEP) $(SWEEP_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SWEEP_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

savings : $(BINDIR)/savings
	true

$(BINDIR)/savings : $(SAVINGS_DEP) $(SAVINGS_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SAVINGS_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

savings2 : $(BINDIR)/savings2
	true

$(BINDIR)/savings2 : $(SAVINGS2_DEP) $(SAVINGS2_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SAVINGS2_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

savings3 : $(BINDIR)/savings3
	true

$(BINDIR)/savings3 : $(SAVINGS3_DEP) $(SAVINGS3_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SAVINGS3_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

route_opt : $(BINDIR)/route_opt
	true

$(BINDIR)/route_opt : $(ROUTE_OPT_DEP) $(ROUTE_OPT_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(ROUTE_OPT_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

nearest_ins : $(BINDIR)/nearest_ins
	true

$(BINDIR)/nearest_ins : $(NEAREST_INS_DEP) $(NEAREST_INS_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(NEAREST_INS_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

farthest_ins : $(BINDIR)/farthest_ins
	true

$(BINDIR)/farthest_ins : $(FARTHEST_INS_DEP) $(FARTHEST_INS_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(FARTHEST_INS_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

farnear_ins : $(BINDIR)/farnear_ins
	true

$(BINDIR)/farnear_ins : $(FARNEAR_INS_DEP) $(FARNEAR_INS_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(FARNEAR_INS_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

near_cluster : $(BINDIR)/near_cluster
	true

$(BINDIR)/near_cluster : $(NEAR_CLUSTER_DEP) $(NEAR_CLUSTER_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(NEAR_CLUSTER_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

tsp_ni : $(BINDIR)/tsp_ni
	true

$(BINDIR)/tsp_ni : $(TSP_NI_DEP) $(TSP_NI_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(TSP_NI_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

tsp_fi : $(BINDIR)/tsp_fi
	true

$(BINDIR)/tsp_fi : $(TSP_FI_DEP) $(TSP_FI_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(TSP_FI_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

tsp_fini : $(BINDIR)/tsp_fini
	true

$(BINDIR)/tsp_fini : $(TSP_FINI_DEP) $(TSP_FINI_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(TSP_FINI_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

exchange : $(BINDIR)/exchange
	true

$(BINDIR)/exchange : $(EXCHANGE_DEP) $(EXCHANGE_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXCHANGE_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

exchange2 : $(BINDIR)/exchange2
	true

$(BINDIR)/exchange2 : $(EXCHANGE2_DEP) $(EXCHANGE2_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(EXCHANGE2_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

mst : $(BINDIR)/mst
	true

$(BINDIR)/mst : $(MST_DEP) $(MST_OBJS)
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MST_OBJS) $(LIBS) \
	-lmaster$(MASTEREXT)
	@echo ""

##############################################################################
# Targets for decomp
##############################################################################

clean_decomp:
	rm -f $(USER_OBJDIR)/my_decomp_user.o 
	$(USER_OBJDIR)/decomp_lower_bound.o

#__END_EXPERIMENTAL_SECTION__#


