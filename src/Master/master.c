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

#define COMPILING_FOR_MASTER

#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __PVM__
#include <pvmtev.h>
#endif

#include "proccomm.h"
#include "timemeas.h"
#include "messages.h"
#include "BB_macros.h"
#include "pack_cut.h"
#include "pack_array.h"
#include "master.h"
#include "master_u.h"
#include "lp_solver.h"
#ifdef COMPILE_IN_TM
#include "tm.h"
#ifdef COMPILE_IN_LP
#include "lp.h"
#endif
#endif

#ifndef TEV_INIT_MASK
/* We must have pvm3.4 where it is called TEV_MASK_INIT */
#  define TEV_INIT_MASK(m)  TEV_MASK_INIT(m)
#  define TEV_SET_MASK(m,k)  TEV_MASK_SET(m,k)
#  define TEV_MCAST0  TEV_MCAST
#  define TEV_RECV0   TEV_RECV
#  define TEV_SEND0   TEV_SEND
#  define TEV_NRECV0  TEV_NRECV
#endif

/*===========================================================================*/
/*===========================================================================*/

problem *sym_open_environment()
{
   problem *p;
#if (!defined(COMPILE_IN_LP) || !defined(COMPILE_IN_CG) || \
   !defined(COMPILE_IN_CP)) && defined(__PVM__)
   int xpvm_tid;
   Pvmtmask trace_mask;
#endif

   setvbuf(stdout, (char *)NULL, _IOLBF, 0);
   
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)
       
   register_process();   /* Enroll this process */

#ifdef __PVM__
   pvm_catchout(stdout); /* Tells PVM to treat all output from the children of
			    this process as output from this process itself*/
#endif
#endif
   
   /* This next set of commands has to be executed if we want to create a PVM
      trace file for viewing in xpvm (this is a very slow process) */

#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
   if (p->par.pvm_trace){
      if ((xpvm_tid = pvm_gettid((char *)"xpvm", 0)) > 0){
	 pvm_setopt(PvmSelfTraceTid, xpvm_tid);
	 pvm_setopt(PvmSelfTraceCode, 666);
	 pvm_setopt(PvmSelfOutputTid, xpvm_tid);
	 pvm_setopt(PvmSelfOutputCode, 667);
	 pvm_setopt(PvmTraceTid, xpvm_tid);
	 pvm_setopt(PvmTraceCode, 666);
	 pvm_setopt(PvmOutputTid, xpvm_tid);
	 pvm_setopt(PvmOutputCode, 667);
	 TEV_INIT_MASK(trace_mask);
	 TEV_SET_MASK(trace_mask, TEV_MCAST0);
	 TEV_SET_MASK(trace_mask, TEV_RECV0);
	 TEV_SET_MASK(trace_mask, TEV_SEND0);
	 TEV_SET_MASK(trace_mask, TEV_NRECV0);
	 pvm_settmask(PvmTaskSelf, trace_mask);
	 pvm_settmask(PvmTaskChild, trace_mask);
      }else{
	 PVM_ERROR(xpvm_tid);
      }
   }
#endif
  
   printf("\n");
   printf("*******************************************************\n");
   printf("*   This is SYMPHONY Version 4.0                      *\n");
   printf("*   Copyright 2000-2003 Ted Ralphs                    *\n");
   printf("*   All Rights Reserved.                              *\n");
   printf("*   Distributed under the Common Public License 1.0   *\n");
   printf("*******************************************************\n");
   printf("\n");
   
   p = (problem *) calloc(1, sizeof(problem));

   if (initialize_u(p) == FUNCTION_TERMINATED_NORMALLY){
      return(p);
   }else{
      FREE(p);
      return(NULL);
   }
}   

/*===========================================================================*/
/*===========================================================================*/

int sym_set_defaults(problem *p)
{
   int termcode = 0;
   
   tm_params *tm_par = &p->par.tm_par;
   lp_params *lp_par = &p->par.lp_par;
   cg_params *cg_par = &p->par.cg_par;
   cp_params *cp_par = &p->par.cp_par;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   sp_params *sp_par = &p->par.sp_par;
#endif
   /*___END_EXPERIMENTAL_SECTION___*/
   dg_params *dg_par = &p->par.dg_par;

   /************************* Global defaults ********************************/
   p->ub = 0;
   p->has_ub = FALSE;
   p->lb = 0;
   p->termcode = TM_NO_PROBLEM;
   p->par.verbosity = 0;
   p->par.random_seed = 17;
   p->par.tm_machine_set = FALSE;
   p->par.dg_machine_set = FALSE;
   strcpy(p->par.tm_exe, "tm");
#ifdef COMPILE_IN_LP
   strcat(p->par.tm_exe, "_lp");
#ifdef COMPILE_IN_CG
   strcat(p->par.tm_exe, "_cg");
#endif
#endif
#ifdef COMPILE_IN_CP
   strcat(p->par.tm_exe, "_cp");
#endif   
   strcpy(p->par.dg_exe, "dg");
   p->par.tm_debug = 0;
   p->par.dg_debug = 0;
   p->par.pvm_trace = 0;
   p->par.do_branch_and_cut = 1;
   p->par.do_draw_graph = FALSE;
   p->par.use_permanent_cut_pools = FALSE;

   /************************** treemanager defaults **************************/
   tm_par->verbosity = 0;
   tm_par->granularity = 0.000001;
   strcpy(tm_par->lp_exe, "lp");
#ifdef COMPILE_IN_CG
   strcat(tm_par->lp_exe, "_cg");
#endif
   strcpy(tm_par->cg_exe, "cg");
   strcpy(tm_par->cp_exe, "cp");
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   strcpy(tm_par->sp_exe, "sp");
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->lp_debug = 0;
   tm_par->cg_debug = 0;
   tm_par->cp_debug = 0;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   tm_par->sp_debug = 0;
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->max_active_nodes = 1;
   tm_par->max_cp_num = 1;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   tm_par->max_sp_num = 0;
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->lp_mach_num = 0;
   tm_par->lp_machs = NULL;
   tm_par->cg_mach_num = 0;
   tm_par->cg_machs = NULL;
   tm_par->cp_mach_num = 0;
   tm_par->cp_machs = NULL;

   tm_par->use_cg = FALSE;
   tm_par->random_seed = 17;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   tm_par->do_decomp = FALSE;
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->unconditional_dive_frac = .1;
   tm_par->diving_strategy = BEST_ESTIMATE;
   tm_par->diving_k = 1;
   tm_par->diving_threshold = 0;
   tm_par->node_selection_rule = LOWEST_LP_FIRST;

   tm_par->keep_description_of_pruned = DISCARD;
   tm_par->warm_start = FALSE;
   tm_par->logging = NO_LOGGING;
   tm_par->logging_interval = 1800;
   tm_par->vbc_emulation = NO_VBC_EMULATION;
   tm_par->price_in_root = FALSE;
   tm_par->trim_search_tree = FALSE;
   tm_par->colgen_strat[0] = (FATHOM__DO_NOT_GENERATE_COLS__DISCARD  |
			      BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   tm_par->colgen_strat[1] = (FATHOM__DO_NOT_GENERATE_COLS__DISCARD  |
			      BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   tm_par->not_fixed_storage_size = 2048;
   tm_par->time_limit = 0;
   tm_par->node_limit = 0;
   tm_par->gap_limit = 0;
   tm_par->find_first_feasible = FALSE;
   
   /************************** lp defaults ***********************************/
   lp_par->verbosity = 0;
   lp_par->granularity = tm_par->granularity;
   lp_par->use_cg = tm_par->use_cg;
   lp_par->set_obj_upper_lim = FALSE;
   lp_par->scaling = -1; /* CPLEX'ism ... don't scale */
   lp_par->fastmip = 1; /* CPLEX'ism ... set it to 1 */
   lp_par->try_to_recover_from_error = TRUE;
   lp_par->problem_type = ZERO_ONE_PROBLEM;
   lp_par->keep_description_of_pruned = tm_par->keep_description_of_pruned;
   lp_par->not_fixed_storage_size = tm_par->not_fixed_storage_size;
   lp_par->cut_pool_check_freq = 10;
   lp_par->load_balance_level = -1;
   lp_par->load_balance_iterations = -1;
   lp_par->load_balance_compare_candidates = HIGHEST_LOW_OBJ;
   lp_par->fractional_diving_ratio = 0.02;
   lp_par->fractional_diving_num = 0;
   lp_par->max_non_dual_feas_to_add_frac = 0.05;
   lp_par->max_non_dual_feas_to_add_min = 20;
   lp_par->max_non_dual_feas_to_add_max = 200;
   lp_par->max_not_fixable_to_add_frac = 0.1;
   lp_par->max_not_fixable_to_add_min = 100;
   lp_par->max_not_fixable_to_add_max = 500;
   lp_par->mat_col_compress_num = 50;
   lp_par->mat_col_compress_ratio = .05;
   lp_par->mat_row_compress_num = 20;
   lp_par->mat_row_compress_ratio = .05;
   lp_par->tailoff_gap_backsteps = 2;
   lp_par->tailoff_gap_frac = .99;
   lp_par->tailoff_obj_backsteps = 3;
   lp_par->tailoff_obj_frac = .75;
   lp_par->tailoff_absolute = 0.0001;
   lp_par->ineff_cnt_to_delete = 0;
   lp_par->eff_cnt_before_cutpool = 3;
   lp_par->ineffective_constraints = BASIC_SLACKS_ARE_INEFFECTIVE;
   lp_par->base_constraints_always_effective = TRUE;
   lp_par->branch_on_cuts = FALSE;
   lp_par->discard_slack_cuts = DISCARD_SLACKS_BEFORE_NEW_ITERATION;
   lp_par->first_lp.first_cut_time_out = 0;
   lp_par->first_lp.all_cuts_time_out = 0;
   lp_par->later_lp.first_cut_time_out = 5;
   lp_par->later_lp.first_cut_time_out = 0;
   lp_par->later_lp.all_cuts_time_out = 1;
   lp_par->later_lp.all_cuts_time_out = 0;
   lp_par->max_cut_num_per_iter = 20;
   lp_par->do_reduced_cost_fixing = TRUE;
   lp_par->gap_as_ub_frac = .1;
   lp_par->gap_as_last_gap_frac = .7;
   lp_par->do_logical_fixing = 1;
   lp_par->fixed_to_ub_before_logical_fixing = 1;
   lp_par->fixed_to_ub_frac_before_logical_fixing = .01;

   lp_par->generate_cgl_cuts = TRUE;

#ifdef __OSI_GLPK__
   lp_par->max_presolve_iter = -1;
#else
   lp_par->max_presolve_iter = 50;
#endif
   
   lp_par->is_feasible_default = TEST_INTEGRALITY;
   lp_par->send_feasible_solution_default = SEND_NONZEROS;
   lp_par->display_solution_default = DISP_NOTHING;
   lp_par->shall_we_branch_default = USER__BRANCH_IF_TAILOFF;
   lp_par->select_candidates_default = USER__CLOSE_TO_HALF;
   lp_par->strong_branching_cand_num_max = 25;
   lp_par->strong_branching_cand_num_min = 5;
   lp_par->strong_branching_red_ratio = 1;
   lp_par->compare_candidates_default = HIGHEST_LOW_OBJ;
   lp_par->select_child_default = PREFER_LOWER_OBJ_VALUE;
   lp_par->pack_lp_solution_default = SEND_NONZEROS;

   /************************** cut_gen defaults *****************************/
   cg_par->verbosity = 0;
   cg_par->do_findcuts = TRUE;

   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   cg_par->do_decomp = FALSE;
   cg_par->decomp_sol_pool_check_freq = 10;
   cg_par->decomp_wait_for_cols = TRUE;
   cg_par->decomp_max_col_num_per_iter = 1000;
   cg_par->decomp_col_block_size = 1000;
   cg_par->decomp_mat_block_size = 100000;
   cg_par->decomp_initial_timeout = 5;
   cg_par->decomp_dynamic_timeout = 5;
   cg_par->decomp_complete_enum = TRUE;

   /*___END_EXPERIMENTAL_SECTION___*/
   /************************** cutpool defaults ******************************/
   cp_par->verbosity = 0;
   cp_par->warm_start = FALSE;
   cp_par->logging = FALSE;
   cp_par->block_size = 5000;
   cp_par->max_size = 2000000;
   cp_par->max_number_of_cuts = 10000;
   cp_par->cuts_to_check = 1000;
   cp_par->delete_which = DELETE_BY_QUALITY;
   cp_par->touches_until_deletion = 10;
   cp_par->min_to_delete = 1000;
   cp_par->check_which = CHECK_ALL_CUTS;

   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   /************************** solpool defaults ******************************/
#ifdef COMPILE_DECOMP
   sp_par->verbosity = 0;
   sp_par->etol = 0.000001;
   sp_par->block_size = 1000;
   sp_par->max_size = 1000000;
   sp_par->max_number_of_sols = 10000;
   sp_par->delete_which = DELETE_DUPLICATE_COLS;
   sp_par->touches_until_deletion = 10;
   sp_par->min_to_delete = 100;
   sp_par->compress_num = 10;
   sp_par->compress_ratio = .01;
   sp_par->check_which = CHECK_COL_LEVEL_AND_TOUCHES;
#endif
   /*___END_EXPERIMENTAL_SECTION___*/
   /********************** draw_graph defaults  ******************************/
   strcpy(dg_par->source_path, ".");
   dg_par->echo_commands = FALSE;
   dg_par->canvas_width = 1000;
   dg_par->canvas_height = 700;
   dg_par->viewable_width = 600;
   dg_par->viewable_height = 400;
   dg_par->disp_nodelabels = 1;
   dg_par->disp_nodeweights = 1;
   dg_par->disp_edgeweights = 1;
   dg_par->node_dash[0] = 0;
   dg_par->edge_dash[0] = 0;
   dg_par->node_radius = 8;
   dg_par->interactive_mode = 1;
   dg_par->mouse_tracking = 1;
   dg_par->scale_factor = 1;
   strcpy(dg_par->nodelabel_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");
   strcpy(dg_par->nodeweight_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");
   strcpy(dg_par->edgeweight_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_parse_command_line(problem *p, int argc, char **argv)
{
   int termcode = 0;

   CALL_WRAPPER_FUNCTION( readparams_u(p, argc, argv) );

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_load_problem(problem *p)
{
   double t = 0;
   int termcode = 0;
 
   /*------------------------------------------------------------------------*\
    *                         start reading in problem                        
   \*------------------------------------------------------------------------*/

   (void) used_time(&t);

   /* Get the problem data */
   CALL_WRAPPER_FUNCTION( io_u(p) );

   /* Start up the graphics window*/
#ifndef WIN32
   CALL_WRAPPER_FUNCTION( init_draw_graph_u(p) );
#endif

   /*------------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*------------------------------------------------------------------------*/

   CALL_WRAPPER_FUNCTION( initialize_root_node_u(p) );

   p->comp_times.readtime = used_time(&t);

   p->termcode = TM_NO_SOLUTION;

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_find_initial_bounds(problem *p)
{
   double total_time = 0;
   int termcode = 0;
   
   /* Finds the upper and lower bounds for the problem */
   CALL_WRAPPER_FUNCTION( start_heurs_u(p) );

   if (!p->par.do_branch_and_cut){
      printf("\n****************************************************\n");
      printf(  "* Heuristics Finished!!!!!!!                       *\n");
      printf(  "* Now displaying stats and best solution....       *\n");
      printf(  "****************************************************\n\n");
      total_time += p->comp_times.ub_overhead + p->comp_times.ub_heurtime;
      total_time += p->comp_times.lb_overhead + p->comp_times.lb_heurtime;
#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
      printf( "  Problem IO     %.3f\n", p->comp_times.readtime);
      printf( "  Overhead: UB   %.3f\n", p->comp_times.ub_overhead);
      printf( "            LB   %.3f\n", p->comp_times.lb_overhead);
      printf( "  Runtime:  UB   %.3f\n", p->comp_times.ub_heurtime);
      printf( "            LB   %.3f\n", p->comp_times.lb_heurtime);
      printf( "  Total User Time    %.3f\n", total_time);
#endif
      if (p->mip->obj_sense == MAXIMIZE){
	 printf( "Upper Bound: %.3f\n", -p->ub + p->mip->obj_offset);
      }else{
 	 printf( "Upper Bound: %.3f\n", p->ub + p->mip->obj_offset);
      } 
      CALL_WRAPPER_FUNCTION( display_solution_u(p, 0) );
      if (p->par.tm_par.lp_machs)
	 FREE(p->par.tm_par.lp_machs[0]);
      FREE(p->par.tm_par.lp_machs);
   }

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_solve(problem *p)
{
   int s_bufid, r_bufid, bytes, msgtag = 0, sender, termcode = 0, temp, i;
   char lp_data_sent = FALSE, cg_data_sent = FALSE, cp_data_sent = FALSE;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   char sp_data_sent = TRUE; /*for now, we are not using this one*/
   /*___END_EXPERIMENTAL_SECTION___*/
#ifndef COMPILE_IN_TM
   char repricing, node_type;
#else
   tm_prob *tm;
#endif
   double start_time, lb;
   struct timeval timeout = {10, 0};
   double t = 0, total_time = 0;

   node_desc *rootdesc = p->rootdesc;
   base_desc *base = p->base;

   start_time = wall_clock(NULL);

#ifndef COMPILE_IN_TM
   /*------------------------------------------------------------------------*\
    * Start the tree manager and send the parameters
   \*------------------------------------------------------------------------*/

   if (p->par.tm_machine_set){
      spawn(p->par.tm_exe, (char **)NULL, p->par.tm_debug | TaskHost,
	    p->par.tm_machine, 1, &p->tm_tid);
   }else{
      spawn(p->par.tm_exe, (char **)NULL, p->par.tm_debug, (char *)NULL, 1,
	    &p->tm_tid);
   }
   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&p->par.tm_par), sizeof(tm_params));
   send_char_array(&p->has_ub, 1);
   if (p->has_ub)
      send_dbl_array(&p->ub, 1);
   send_char_array(&p->has_ub_estimate, 1);
   if (p->has_ub_estimate)
      send_dbl_array(&p->ub_estimate, 1);
   if (p->par.tm_par.lp_mach_num)
      send_char_array(p->par.tm_par.lp_machs[0],
		      p->par.tm_par.lp_mach_num*MACH_NAME_LENGTH);
   if (p->par.tm_par.cg_mach_num)
      send_char_array(p->par.tm_par.cg_machs[0],
		      p->par.tm_par.cg_mach_num*MACH_NAME_LENGTH);
   if (p->par.tm_par.cp_mach_num)
      send_char_array(p->par.tm_par.cp_machs[0],
		      p->par.tm_par.cp_mach_num*MACH_NAME_LENGTH);
   send_int_array(&base->varnum, 1);
   send_int_array(&base->cutnum, 1);
#ifdef TRACE_PATH
   {
      int feas_sol;
      int *feas_sol_size;

      if (user_send_feas_sol(p->user, &feas_sol_size, &feas_sol)==USER_NO_PP){
	 send_int_array(&feas_sol_size, 1);
	 if (feas_sol_size){
	    send_int_array(feas_sol, feas_sol_size);
	 }
      }
   }
#endif   
   send_msg(p->tm_tid, TM_DATA);
      
   /*------------------------------------------------------------------------*\
    * Send out the root node
   \*------------------------------------------------------------------------*/

   if (!p->par.warm_start){
      repricing = FALSE;
      node_type = ROOT_NODE;
      
      s_bufid = init_send(DataInPlace);
      send_char_array(&repricing, 1);
      send_char_array(&node_type, 1);
      send_dbl_array(&p->lb, 1);
      send_int_array(&rootdesc->nf_status, 1);
      pack_array_desc(&rootdesc->uind);
      if (rootdesc->nf_status == NF_CHECK_AFTER_LAST ||
	  rootdesc->nf_status == NF_CHECK_UNTIL_LAST)
	 pack_array_desc(&rootdesc->not_fixed);
      pack_array_desc(&rootdesc->cutind);
      pack_basis(&rootdesc->basis, TRUE);
      send_int_array(&rootdesc->desc_size, 1);
      if (rootdesc->desc_size)
	 send_char_array(rootdesc->desc, rootdesc->desc_size);
      if (rootdesc->cutind.size > 0){ /* Hey, we have cuts! Pack them, too. */
	 /* Pack their number again, so we can call unpack_cut_set in TM */
	 int i;
	 send_int_array(&rootdesc->cutind.size, 1);
	 for (i = 0; i < rootdesc->cutind.size; i++)
	    pack_cut(rootdesc->cuts[i]);
      }
      send_msg(p->tm_tid, TM_ROOT_DESCRIPTION);
      freebuf(s_bufid);
   }
#else
   
   /*------------------------------------------------------------------------*\
    * Create the treemanager and copy the problem data
   \*------------------------------------------------------------------------*/

   p->tm = tm = (tm_prob *) calloc(1, sizeof(tm_prob));

   tm->par = p->par.tm_par;

   if ((tm->has_ub = p->has_ub))
      tm->ub = p->ub;
   if ((tm->has_ub_estimate = p->has_ub_estimate))
      tm->ub_estimate = p->ub_estimate;

#ifdef COMPILE_IN_LP
   CALL_WRAPPER_FUNCTION( send_lp_data_u(p, 0) );
   lp_data_sent = TRUE;
#ifdef COMPILE_IN_CG
   CALL_WRAPPER_FUNCTION( send_cg_data_u(p, 0) );
   cg_data_sent = TRUE;
#endif
#endif
#ifdef COMPILE_IN_CP
   if (p->cp && p->par.use_permanent_cut_pools){
      tm->cpp = p->cp;
   }else{
      CALL_WRAPPER_FUNCTION( send_cp_data_u(p, 0) );
   }
   cp_data_sent = TRUE;
#endif

   if (p->warm_start && p->par.tm_par.warm_start){
      /* Load warm start info */
      tm->rootnode = p->warm_start->rootnode;
      tm->cuts = p->warm_start->cuts;
      tm->cut_num = p->warm_start->cut_num;
      tm->allocated_cut_num = p->warm_start->allocated_cut_num;
      tm->stat = p->warm_start->stat;
      tm->comp_times = p->warm_start->comp_times;
      tm->lb = p->warm_start->lb;
      if (p->warm_start->has_ub){
	 if (p->warm_start->ub < tm->ub || !tm->has_ub){
	    tm->ub = p->warm_start->ub;
	 }
	 tm->has_ub = TRUE;
      }
      tm->phase = p->warm_start->phase;
   }else if (p->warm_start){
      /* Otherwise, free what was saved */
      free_subtree(p->warm_start->rootnode);
      if (p->warm_start->cuts){
	 for (i = p->warm_start->cut_num - 1; i >= 0; i--)
	    if (p->warm_start->cuts[i]){
	       FREE(p->warm_start->cuts[i]->coef);
	       FREE(p->warm_start->cuts[i]);
	    }
	 FREE(p->warm_start->cuts);
      }
   }
   /* Now the tree manager owns everything */
   FREE(p->warm_start);
   
   if ((termcode = tm_initialize(tm , base, rootdesc)) < 0){
      tm_close(tm, termcode);

      if (p->par.do_draw_graph){
	 s_bufid = init_send(DataInPlace);
	 send_msg(p->dg_tid, CTOI_YOU_CAN_DIE);
	 freebuf(s_bufid);
      }
      
      if (p->par.tm_par.lp_machs)
	 FREE(p->par.tm_par.lp_machs[0]);
      FREE(p->par.tm_par.lp_machs);
      if (p->par.tm_par.cg_machs)
	 FREE(p->par.tm_par.cg_machs[0]);
      FREE(p->par.tm_par.cg_machs);
      if (p->par.tm_par.cp_machs)
	 FREE(p->par.tm_par.cp_machs[0]);
      FREE(p->par.tm_par.cp_machs);
      
      free_tm(tm);

      p->termcode = termcode;
      
      return(termcode);
   }
   

#ifdef TRACE_PATH
   {
      int feas_sol_size;
      int *feas_sol;
      
      if (user_send_feas_sol(p->user, &feas_sol_size, &feas_sol)==USER_NO_PP){
	 tm->feas_sol_size = feas_sol_size;
	 tm->feas_sol = (int *) calloc (tm->feas_sol_size, sizeof(int));
	 memcpy((char *)tm->feas_sol, (char *)feas_sol, feas_sol_size * ISIZE);
      }
   }
#endif
#endif   
   
   /*------------------------------------------------------------------------*\
    * Wait for messages
   \*------------------------------------------------------------------------*/
   
#ifdef COMPILE_IN_TM
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   while (!lp_data_sent || !cg_data_sent || !cp_data_sent || !sp_data_sent){
   /*___END_EXPERIMENTAL_SECTION___*/
   /*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
   while (!lp_data_sent || !cg_data_sent || !cp_data_sent){
#endif
#else
   do{
#endif
      r_bufid = treceive_msg(ANYONE, ANYTHING, &timeout);
      if (r_bufid == 0){
#ifndef COMPILE_IN_TM
	 if (pstat(p->tm_tid) != PROCESS_OK){
	    printf("\nThe treemanager has died :-(\n\n");
#else
	 if (!processes_alive(p->tm)){
#endif
	    termcode = msgtag = SOMETHING_DIED;
	    break;
	 }else{
	    continue;
	 }
      }
      bufinfo(r_bufid, &bytes, &msgtag, &sender);

      switch (msgtag){
       case FEASIBLE_SOLUTION_NONZEROS:
       case FEASIBLE_SOLUTION_USER:
	 CALL_WRAPPER_FUNCTION( receive_feasible_solution_u(p, msgtag) );
	 if (p->par.verbosity > 0){
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
	    CALL_WRAPPER_FUNCTION( display_solution_u(p, p->tm->opt_thread_num) );
#else
	    CALL_WRAPPER_FUNCTION( display_solution_u(p, 0) );
#endif
	 }
	 break;

       case REQUEST_FOR_LP_DATA:
	 /* An LP process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_lp_data_u(p, sender) );
	 lp_data_sent = TRUE;
	 break;

       case REQUEST_FOR_CG_DATA:
	 /* A CG process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_cg_data_u(p, sender) );
	 cg_data_sent = TRUE;
	 break;

       case REQUEST_FOR_CP_DATA:
	 /* A CP process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_cp_data_u(p, sender) );
	 cp_data_sent = TRUE;
	 break;

       /*__BEGIN_EXPERIMENTAL_SECTION__*/
       case REQUEST_FOR_SP_DATA:
	 /* An SP process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_sp_data_u(p, sender) );
	 sp_data_sent = TRUE;
	 break;

       /*___END_EXPERIMENTAL_SECTION___*/
       case TM_FIRST_PHASE_FINISHED:
	 receive_char_array((char *)(&p->comp_times.bc_time),
			     sizeof(node_times));
	 receive_dbl_array(&lb, 1);
	 if (lb > p->lb) p->lb = lb;
	 receive_char_array((char *)&p->warm_start->stat, sizeof(problem_stat));
	 printf( "\n");
	 printf( "****************************************************\n");
	 printf( "* Branch and Cut First Phase Finished!!!!          *\n");
	 printf( "* Now displaying stats and best solution...        *\n");
	 printf( "****************************************************\n\n");

	 print_statistics(&(p->comp_times.bc_time), &(p->warm_start->stat),
			  p->ub, p->lb, 0, start_time);
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
	 CALL_WRAPPER_FUNCTION( display_solution_u(p, p->tm->opt_thread_num) );
#else
	 CALL_WRAPPER_FUNCTION( display_solution_u(p, 0) );
#endif
	 break;

       case SOMETHING_DIED:
       case TM_TIME_LIMIT_EXCEEDED:
       case TM_NODE_LIMIT_EXCEEDED:
       case TM_TARGET_GAP_ACHIEVED:
       case TM_FOUND_FIRST_FEASIBLE:
       case TM_OPTIMAL_SOLUTION_FOUND:
       case TM_ERROR__NO_BRANCHING_CANDIDATE:
       case TM_ERROR__ILLEGAL_RETURN_CODE:
       case TM_ERROR__NUMERICAL_INSTABILITY:
       case TM_ERROR__COMM_ERROR:
       case TM_ERROR__USER:
	 receive_char_array((char *)(&p->comp_times.bc_time),
			    sizeof(node_times));
	 receive_dbl_array(&lb, 1);
	 if (lb > p->lb) p->lb = lb;
	 receive_char_array((char *)&p->warm_start->stat, sizeof(problem_stat));
	 break;

       default:
	 CALL_WRAPPER_FUNCTION( process_own_messages_u(p, msgtag) );
	 break;
      }
      freebuf(r_bufid);

#ifndef COMPILE_IN_TM
   }while (msgtag != TM_OPTIMAL_SOLUTION_FOUND && msgtag != SOMETHING_DIED &&
	   msgtag != TM_TIME_LIMIT_EXCEEDED &&
	   msgtag != TM_NODE_LIMIT_EXCEEDED &&
	   msgtag != TM_TARGET_GAP_ACHIEVED &&
	   msgtag != TM_FOUND_FIRST_FEASIBLE &&
	   msgtag != TM_ERROR__NO_BRANCHING_CANDIDATE &&
	   msgatg != TM_ERROR__ILLEGAL_RETURN_CODE &&
	   msgtag != TM_ERROR__NUMERICAL_INSTABLITY &&
	   msgtag != TM_ERROR__COMM_ERROR &&
	   msgtag != TM_ERROR__USER);

   termcode = msgtag;
#else
   }
   
   /*------------------------------------------------------------------------*\
    * Solve the problem and receive solutions                         
   \*------------------------------------------------------------------------*/

   tm->start_time += start_time;

   termcode = solve(tm);

   /* Save the warm start info */
   p->warm_start = (warm_start_desc *) calloc (1, sizeof(warm_start_desc));
   p->warm_start->rootnode = tm->rootnode;
   p->warm_start->cuts = p->tm->cuts;
   p->warm_start->cut_num = p->tm->cut_num;
   p->warm_start->allocated_cut_num = p->tm->allocated_cut_num;
   p->warm_start->stat = tm->stat;
   p->warm_start->phase = tm->phase;
   p->warm_start->lb = tm->lb;
   if (p->warm_start->has_ub = tm->has_ub){
      p->warm_start->ub = tm->ub;
   }
   tm->rootnode = NULL;
   tm->cuts = NULL;
   tm->cut_num = tm->allocated_cut_num = 0;
   if (p->cp && p->par.use_permanent_cut_pools){
      tm->cpp = NULL;
   }
   tm_close(tm, termcode);

#ifndef COMPILE_IN_LP
   if (termcode != SOMETHING_DIED){
      do{
	 r_bufid = receive_msg(ANYONE, ANYTHING);
	 if (r_bufid == 0){
	    printf("\nError receiving solution ...\n");
	    break;
	 }
	 bufinfo(r_bufid, &bytes, &msgtag, &sender);
	 if (msgtag == FEASIBLE_SOLUTION_NONZEROS ||
	     msgtag == FEASIBLE_SOLUTION_USER){
	    CALL_WRAPPER_FUNCTION( receive_feasible_solution_u(p, msgtag) );
	 }
      }while (msgtag != FEASIBLE_SOLUTION_NONZEROS &&
	      msgtag != FEASIBLE_SOLUTION_USER);
   }
#endif
#endif

   /*------------------------------------------------------------------------*\
    * Display the the results and solution data                               
   \*------------------------------------------------------------------------*/

   printf("\n****************************************************\n");
   if (termcode == TM_OPTIMAL_SOLUTION_FOUND){
      printf(  "* Branch and Cut Finished                          *\n");
   }else if (termcode == TM_TIME_LIMIT_EXCEEDED){
      printf(  "* Time Limit Reached                               *\n");
   }else if (termcode == TM_NODE_LIMIT_EXCEEDED){
      printf(  "* Node Limit Reached                               *\n");
   }else if (termcode == TM_TARGET_GAP_ACHIEVED){
      printf(  "* Target Gap Achieved                              *\n");
   }else if (termcode == TM_FOUND_FIRST_FEASIBLE){
      printf(  "* Stopping After Finding First Feasible Solution   *\n");
   }else if (termcode == TM_ERROR__NO_BRANCHING_CANDIDATE ||
	     termcode == TM_ERROR__ILLEGAL_RETURN_CODE ||
	     termcode == TM_ERROR__NUMERICAL_INSTABILITY ||
	     termcode == TM_ERROR__COMM_ERROR ||
	     termcode == TM_ERROR__USER){
      printf(  "* Terminated abnormally with error message %i      *\n",termcode);
   }else{
      printf("* A process has died abnormally -- halting \n\n");
   }
      printf(  "* Now displaying stats and best solution found...  *\n");
      printf(  "****************************************************\n\n");

   total_time  = p->comp_times.readtime;
   total_time += p->comp_times.ub_overhead + p->comp_times.ub_heurtime;
   total_time += p->comp_times.lb_overhead + p->comp_times.lb_heurtime;

#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
   printf( "====================== Misc Timing =========================\n");
   printf( "  Problem IO        %.3f\n", p->comp_times.readtime);
   printf( "  UB overhead:      %.3f\n", p->comp_times.ub_overhead);
   printf( "  UB runtime:       %.3f\n", p->comp_times.ub_heurtime);
   printf( "  LB overhead:      %.3f\n", p->comp_times.lb_overhead);
   printf( "  LB runtime:       %.3f\n", p->comp_times.lb_heurtime);
#endif
   
#ifdef COMPILE_IN_TM
   if (tm->lb > p->lb) p->lb = tm->lb;
   print_statistics(&(tm->comp_times), &(tm->stat), tm->ub, p->lb, total_time,
		    start_time);

   temp = termcode;
#ifdef COMPILE_IN_LP
   CALL_WRAPPER_FUNCTION( display_solution_u(p, p->tm->opt_thread_num) );
#else
   CALL_WRAPPER_FUNCTION( display_solution_u(p, 0) );
#endif
#else
   print_statistics(&(p->comp_times.bc_time), &(p->warm_start->stat), p->ub,
		    p->lb, 0, start_time);
   CALL_WRAPPER_FUNCTION( display_solution_u(p, 0) );
#endif
   termcode = temp;
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   if (p->tm && p->tm->lpp[p->tm->opt_thread_num]){
      p->best_sol = p->tm->lpp[p->tm->opt_thread_num]->best_sol;
      p->tm->lpp[p->tm->opt_thread_num]->best_sol.xlength = 0;
      p->tm->lpp[p->tm->opt_thread_num]->best_sol.xind = NULL;
      p->tm->lpp[p->tm->opt_thread_num]->best_sol.xval = NULL;
   }
#endif

   if (p->par.do_draw_graph){
      s_bufid = init_send(DataInPlace);
      send_msg(p->dg_tid, CTOI_YOU_CAN_DIE);
      freebuf(s_bufid);
   }

   if (p->par.tm_par.lp_machs)
      FREE(p->par.tm_par.lp_machs[0]);
   FREE(p->par.tm_par.lp_machs);
   if (p->par.tm_par.cg_machs)
      FREE(p->par.tm_par.cg_machs[0]);
   FREE(p->par.tm_par.cg_machs);
   if (p->par.tm_par.cp_machs)
      FREE(p->par.tm_par.cp_machs[0]);
   FREE(p->par.tm_par.cp_machs);

   free_tm(tm);

   p->termcode = termcode;
   
   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_create_permanent_cut_pools(problem *p)
{
#if !(defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP))
   return(0);
#else
   int i;
   
   if (p->par.tm_par.max_cp_num){
      p->cp =
	 (cut_pool **) malloc(p->par.tm_par.max_cp_num*sizeof(cut_pool *));
      for (i = 0; i < p->par.tm_par.max_cp_num; i++){
	 p->cp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
	 p->cp[i]->par = p->par.cp_par;
	 CALL_USER_FUNCTION( user_send_cp_data(p->user, &p->cp[i]->user) );
      }
      get_cp_ptr(p->cp, 0);
      return(p->par.tm_par.max_cp_num);
   }else{
      return(0);
   }
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int sym_close_environment(problem *p)
{
   int termcode = 0;
   
   CALL_WRAPPER_FUNCTION( free_master_u(p) );

   FREE(p);

#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
   pvm_catchout(0);
   comm_exit();
#endif

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

