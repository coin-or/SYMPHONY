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

