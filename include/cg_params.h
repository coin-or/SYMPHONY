/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2002 Ted Ralphs. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _CUT_GEN_PARAMS_H
#define _CUT_GEN_PARAMS_H

/*stores the parameters needed by the cut generator*/
typedef struct CG_PARAMS{
   int     verbosity;
   int     do_findcuts;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int     do_decomp;
   int     decomp_sol_pool_check_freq;
   int     decomp_wait_for_cols;
   int     decomp_max_col_num_per_iter;
   int     decomp_col_block_size;
   int     decomp_mat_block_size;
   double  decomp_initial_timeout;
   double  decomp_dynamic_timeout;
   int     decomp_complete_enum;
   /*___END_EXPERIMENTAL_SECTION___*/
}cg_params;

#endif
