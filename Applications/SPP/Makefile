USER_COMMON_SRC   = spp_common.c

USER_LP_SRC       = spp_lp_branch.c spp_lp.c
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += spp_cg.c
endif

USER_TM_SRC       =
ifeq ($(COMPILE_IN_LP),TRUE)
USER_TM_SRC      += $(USER_LP_SRC)
endif
ifeq ($(COMPILE_IN_CP),TRUE)
USER_TM_SRC      += spp_cp.c
endif

USER_MASTER_SRC   = spp_master.c spp_master_functions.c $(USER_COMMON_SRC)

ifeq ($(COMPILE_IN_TM),TRUE)
USER_MASTER_SRC  += $(USER_TM_SRC)
endif

USER_CG_SRC       = spp_cg.c 

USER_CP_SRC       = spp_cp.c 

USER_DG_SRC       = spp_dg.c 

ALL_USER_SRC      = spp_common.c spp_lp_branch.c spp_lp.c
ALL_USER_SRC     += spp_master.c app_master_functions.c 
ALL_USER_SRC     += spp_cg.c spp_cp.c spp_dg.c

USER_SRCDIR      += 

