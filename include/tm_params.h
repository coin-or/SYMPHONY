/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _TM_PARAMS_H
#define _TM_PARAMS_H

#include "proto.h"

typedef struct TM_STAT{
   double      root_lb;
   int         cuts_in_pool;
   int         max_depth;          /* keeps track of the deepest level reached
				      in the tree so far */
   int         chains;             /* the number of diving chains */
   int         diving_halts;       /* how many times was an already started
				      dive stopped */
   int         tree_size;          /* number of search tree nodes */
   int         created;            /* the number of created nodes (not
				      necessarily the same as tree_size
				      (trimming...) */
   int         analyzed;           /* the number of analyzed (i.e., CG-LP
				      iteration) nodes (not necessarily same
				      as created, leaves can be cut off
				      without analyzing; trimming) */
   int         leaves_before_trimming;
   int         leaves_after_trimming;
   int         vars_not_priced;    /* How many variables did not price out
				      after the first phase */
   char        nf_status;          /* nf_status of the root node after
				      repricing */
}tm_stat;

/*===========================================================================*\
 * The params structure contains all of the user-specified parameters   
 * to be read in from the parameter file.
\*===========================================================================*/

typedef struct TM_PARAMS{
   int         verbosity;
   double      granularity;
   char        lp_exe[MAX_FILE_NAME_LENGTH +1];
   char        cg_exe[MAX_FILE_NAME_LENGTH +1];
   char        cp_exe[MAX_FILE_NAME_LENGTH +1];
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   char        sp_exe[MAX_FILE_NAME_LENGTH +1];
   /*___END_EXPERIMENTAL_SECTION___*/
   int         lp_debug;
   int         cg_debug;
   int         cp_debug;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int         sp_debug;
   /*___END_EXPERIMENTAL_SECTION___*/
   int         max_active_nodes;
   int         max_cp_num;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int         max_sp_num;
   /*___END_EXPERIMENTAL_SECTION___*/

   /* if a ..._machine_num is not 0 and there MUST be that many machine
      names listed in ..._machines (one name can be listed more than once) */
   int         lp_mach_num;
   char      **lp_machs;
   int         cg_mach_num;
   char      **cg_machs;
   int         cp_mach_num;
   char      **cp_machs;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int         sp_mach_num;
   char      **sp_machs;
   /*___END_EXPERIMENTAL_SECTION___*/

   int         use_cg;

   int         random_seed;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int         do_decomp;
   /*___END_EXPERIMENTAL_SECTION___*/
   double      unconditional_dive_frac;
   int         diving_strategy;
   int         diving_k;
   double      diving_threshold;
   int         node_selection_rule;

   int         keep_description_of_pruned;
   int         vbc_emulation;
   char        vbc_emulation_file_name[MAX_FILE_NAME_LENGTH +1];
   int         warm_start;
   int         logging;
   int         logging_interval;
   int         cp_logging;
   char        pruned_node_file_name[MAX_FILE_NAME_LENGTH +1];
   char        warm_start_tree_file_name[MAX_FILE_NAME_LENGTH +1];
   char        warm_start_cut_file_name[MAX_FILE_NAME_LENGTH +1];
   char        tree_log_file_name[MAX_FILE_NAME_LENGTH +1];
   char        cut_log_file_name[MAX_FILE_NAME_LENGTH +1];
   int         price_in_root;
   int         trim_search_tree;

   int         colgen_strat[2]; /* the column generattion strategy for the LP
				   in the two phases */
   int         not_fixed_storage_size;
   double      time_limit;
}tm_params;

#endif
