USER_LP_SRC       = mpp_lp_branch.c mpp_lp.c 
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += mpp_cg.c
endif

USER_TM_SRC       =
ifeq ($(COMPILE_IN_LP),TRUE)
USER_TM_SRC      += $(USER_LP_SRC)
endif
ifeq ($(COMPILE_IN_CP),TRUE)
USER_TM_SRC      += mpp_cp.c 
endif

USER_MASTER_SRC   = mpp_master.c

ifeq ($(COMPILE_IN_TM),TRUE)
USER_MASTER_SRC  += $(USER_TM_SRC)
endif

USER_CG_SRC       = mpp_cg.c

USER_CP_SRC       = mpp_cp.c

USER_DG_SRC       = mpp_dg.c

ALL_USER_SRC      = mpp_lp_branch.c mpp_lp.c
ALL_USER_SRC     += mpp_master.c mpp_cg.c mpp_cp.c mpp_dg.c

USER_SRCDIR      += 

