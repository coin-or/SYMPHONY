USER_LP_SRC       = lp_user_branch.c lp_user.c 
ifeq ($(COMPILE_IN_CG),TRUE)
USER_LP_SRC      += cg_user.c
endif

USER_TM_SRC       =
ifeq ($(COMPILE_IN_LP),TRUE)
USER_TM_SRC      += $(USER_LP_SRC)
endif
ifeq ($(COMPILE_IN_CP),TRUE)
USER_TM_SRC      += cp_user.c 
endif

USER_MASTER_SRC   = master_user.c

ifeq ($(COMPILE_IN_TM),TRUE)
USER_MASTER_SRC  += $(USER_TM_SRC)
endif

USER_CG_SRC       = cg_user.c

USER_CP_SRC       = cp_user.c

USER_DG_SRC       = dg_user.c

ALL_USER_SRC      = lp_user_branch.c lp_user.c
ALL_USER_SRC     += master_user.c cg_user.c cp_user.c dg_user.c

USER_SRCDIR      += 

