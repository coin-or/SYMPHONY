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

#include "symphony_api.h"
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
   p->par.use_permanent_cut_pools = TRUE;
   p->par.multi_criteria = FALSE;
   p->par.mc_binary_search_tolerance = 0; 
   p->par.mc_compare_solution_tolerance = .001;
   p->par.mc_search_order = MC_FIFO;
   p->par.mc_warm_start = TRUE;

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

   lp_par->multi_criteria = FALSE;
   lp_par->mc_find_nondominated_solutions = TRUE;
   lp_par->mc_gamma = 1;       /* Determines the weight on objective 1 */
   lp_par->mc_tau   = 0;       /* Determines the weight on objective 2 */
   lp_par->mc_rho   = 0.00001; /* For augmented Chebyshev norm */
   
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

   if (p->par.tm_par.node_selection_rule == BEST_FIRST_SEARCH){
      switch (p->mip->obj_sense){
       case MAXIMIZE:
	  p->par.tm_par.node_selection_rule = HIGHEST_LP_FIRST;
       case MINIMIZE:
	  p->par.tm_par.node_selection_rule = LOWEST_LP_FIRST;
      }
   }
   
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
      if (p->has_ub){
	 if (p->mip->obj_sense == MAXIMIZE){
	    printf( "Lower Bound: %.3f\n", -p->ub + p->mip->obj_offset);
	 }else{
	    printf( "Upper Bound: %.3f\n", p->ub + p->mip->obj_offset);
	 } 
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
      p->best_sol = p->warm_start->best_sol;
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
      FREE(p->best_sol.xind);
      FREE(p->best_sol.xval);
      memset(&(p->best_sol), 0, sizeof(lp_sol));
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
   while (!lp_data_sent || !cg_data_sent || !cp_data_sent)
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
	 receive_char_array((char *)&p->warm_start->stat,sizeof(problem_stat));
	 printf( "\n");
	 printf( "****************************************************\n");
	 printf( "* Branch and Cut First Phase Finished!!!!          *\n");
	 printf( "* Now displaying stats and best solution...        *\n");
	 printf( "****************************************************\n\n");

	 print_statistics(&(p->comp_times.bc_time), &(p->warm_start->stat),
			  p->ub, p->lb, 0, start_time, p->mip->obj_offset,
			  p->mip->obj_sense, p->has_ub);
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
	 receive_char_array((char *)&p->warm_start->stat,sizeof(problem_stat));
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
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   if (p->tm){
      int thread_num = p->tm->opt_thread_num;
      if (p->tm->lpp[thread_num]){
	 if (p->tm->lpp[thread_num]->best_sol.xlength){
	    p->warm_start->best_sol = p->tm->lpp[thread_num]->best_sol;
	 }else if (!p->par.multi_criteria){
	    p->tm->lpp[thread_num]->best_sol = p->best_sol;
	 }
      }
   }
#else
   p->warm_start->best_sol = p->best_sol;
#endif
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
      printf(  "* Terminated abnormally with error message %i      *\n",
	       termcode);
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
		    start_time, p->mip->obj_offset, p->mip->obj_sense,
		    p->has_ub);

   temp = termcode;
#ifdef COMPILE_IN_LP
   CALL_WRAPPER_FUNCTION( display_solution_u(p, p->tm->opt_thread_num) );
#else
   CALL_WRAPPER_FUNCTION( display_solution_u(p, 0) );
#endif
#else
   print_statistics(&(p->comp_times.bc_time), &(p->warm_start->stat), p->ub,
		    p->lb, 0, start_time, p->mip->obj_offset,
		    p->mip->obj_sense, p->has_ub);
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

int sym_resolve(problem *p)
{
   int i, change_type;

   /* first check for the updates! */

   if (p->warm_start && p->par.tm_par.warm_start){
      if (!p->mip->change_num){
	 return sym_solve(p);
      }
      else{
	 if (!p->warm_start){
	    printf("sym_solve():");
	    printf("Unable to process an empty warm start description!\n");
	    return(TM_NO_SOLUTION);
	 }
	 
	 p->warm_start->has_ub = FALSE;
	 p->warm_start->ub = 0.0;

	 for(i = 0; i < p->mip->change_num; i++){
	    change_type = p->mip->change_type[i];
	    if (change_type == RHS_CHANGED){

#ifdef USE_CGL_CUTS
	       printf("sym_resolve(): SYMPHONY can not resolve for the\n");
	       printf("rhs change when cuts exist, for now!\n"); 
	       return(TM_NO_SOLUTION);
	       
#else

	       p->mip->change_num = 0;
	       update_tree_bound(p, p->warm_start->rootnode, RHS_CHANGED);
	       return sym_solve(p);
#endif
	    }
	 }

	 for(i = 0; i < p->mip->change_num; i++){
	    change_type = p->mip->change_type[i];
	    switch(change_type){
	     case OBJ_COEFF_CHANGED:
		update_tree_bound(p, p->warm_start->rootnode, 
				      OBJ_COEFF_CHANGED);
		break;
	     default:
		printf("sum_resolve():");
		printf("Unable to re-solve this change,for now!\n");
		return(TM_NO_SOLUTION); 
	    }
	 }      
      }   
   }
   
   p->mip->change_num = 0;
   return sym_solve(p);
   
}

/*===========================================================================*/
/* These data types are for multi-criteria problems and are only used here   */
/*===========================================================================*/
 
typedef struct SOLUTION_DATA{
   double  obj[2];
   double  gamma;
   double  tau;
   int     length;
   int    *indices;
   double *values;
}solution_data;

/*===========================================================================*/

typedef struct SOLUTION_PAIRS{
   int solution1;
   int solution2;
   double gamma1;
   double gamma2;
}solution_pairs;

/*===========================================================================*/

#define MAX_NUM_PAIRS 10000
#define MAX_NUM_SOLUTIONS 10000
#define MAX_NUM_INFEASIBLE 10000

/*===========================================================================*/

int sym_mc_solve(problem *p)
{
   int i;
   double gamma, gamma0, gamma1, tau, slope;
   double start_time;
   warm_start_desc *ws;

   solution_data utopia1;
   solution_data utopia2;
   solution_data solutions[MAX_NUM_PAIRS];
   int numsolutions = 0, numprobs = 0, numinfeasible = 0;
   solution_pairs pairs[MAX_NUM_PAIRS];
   int numpairs = 0, cur_position = 0, first = 0, last = 0, previous = 0;
   int *indices;
   double *values;
   int length, termcode;
   int solution1, solution2;
   double utopia[2];
   node_desc *root= NULL;
   base_desc *base = NULL;
   double compare_sol_tol, ub = 0.0;
   int binary_search = FALSE;
   
   for (i = 0; i < p->mip->n; i++){
      if (p->mip->obj2[i] != 0){
	 break;
      }
   }
   if (i == p->mip->n){
      printf("Second objective function is identically zero.\n");
      printf("Switching to standard branch and bound.\n\n");
      return(sym_solve(p));
   }

   sym_set_int_param(p, "multi_criteria", TRUE);
   memcpy((char *)p->mip->obj1, (char *)p->mip->obj, DSIZE*p->mip->n);
   if (p->par.lp_par.mc_find_nondominated_solutions){
      p->base->cutnum += 2;
      p->rootdesc->uind.size++;
      p->rootdesc->uind.list = (int *) realloc(p->rootdesc->uind.list,
					       p->rootdesc->uind.size*ISIZE);
      p->rootdesc->uind.list[p->rootdesc->uind.size-1] = p->mip->n;
   }else{
      sym_set_int_param(p, "keep_description_of_pruned", KEEP_IN_MEMORY);
   }
   
   start_time = wall_clock(NULL);

   /* Set some parameters */
   compare_sol_tol = p->par.mc_compare_solution_tolerance;
   p->par.tm_par.granularity = p->par.lp_par.granularity =
      -MAX(p->par.lp_par.mc_rho, compare_sol_tol);

   if (p->par.verbosity >= 0){
      if (p->par.mc_binary_search_tolerance > 0){
	 binary_search = TRUE;
	 printf("Using binary search with tolerance = %f...\n",
		p->par.mc_binary_search_tolerance);
      }
      if (p->par.mc_search_order == MC_LIFO){
	 printf("Using LIFO search order...\n");
      }else{
	 printf("Using FIFO search order...\n");
      }
      if (p->par.lp_par.mc_rho > 0){
	 printf("Using augmented Chebyshev weight %.8f\n",
		p->par.lp_par.mc_rho);
      }
      if (p->par.use_permanent_cut_pools){
	 printf("Saving the global cut pool between iterations...\n");
	 sym_create_permanent_cut_pools(p);
      }
      printf("\n");
   }

   /* First, calculate the utopia point */
   p->par.lp_par.mc_gamma = 1.0;
   p->par.lp_par.mc_tau = 0.0;
      
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = 1.0 tau = 0.0 \n");  
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Solve */
   if (termcode = sym_solve(p) < 0){
      p->base->cutnum -=2;
      p->rootdesc->uind.size--;
      return(termcode);
   }
   numprobs++;
   
   if (!p->par.lp_par.mc_find_nondominated_solutions){
      ws = sym_get_warm_start(p, TRUE);
   }
   
   /* Store the solution */
   length = solutions[numsolutions].length = p->best_sol.xlength;
   indices = solutions[numsolutions].indices = (int *) calloc(length, ISIZE);
   values = solutions[numsolutions].values = (double *) calloc(length, DSIZE);
   memcpy((char *) indices, p->best_sol.xind, length * ISIZE);
   memcpy((char *) values, p->best_sol.xval, length * DSIZE);
   solutions[numsolutions].gamma = 1.0;
   solutions[numsolutions].tau = 0.0;
   solutions[numsolutions].obj[0] = p->obj[0];
   solutions[numsolutions++].obj[1] = p->obj[1];
   utopia[0] = p->obj[0];
      
   p->par.lp_par.mc_gamma = 0.0;
   p->par.lp_par.mc_tau = 1.0;
      
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = 0.0 tau = 1.0 \n");  
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Resolve */
   if (!p->par.lp_par.mc_find_nondominated_solutions){
      sym_set_warm_start(p, ws);
      for (i = 0; i < p->mip->n; i++){
	 sym_set_obj_coeff(p, i, p->mip->obj2[i] +
			   p->par.lp_par.mc_rho*(p->mip->obj1[i] +
					      p->mip->obj2[i]));
      }
      if (termcode = sym_resolve(p) < 0){
	 sym_delete_warm_start(ws);
	 p->base->cutnum -=2;
	 p->rootdesc->uind.size--;
	 return(termcode);
      }
   }else{
      if (termcode = sym_solve(p) < 0){
	 p->base->cutnum -=2;
	 p->rootdesc->uind.size--;
	 return(termcode);
      }
   }      
   numprobs++;
   
   /* Store the solution */
   length = solutions[numsolutions].length = p->best_sol.xlength;
   indices = solutions[numsolutions].indices = (int *) calloc(length, ISIZE);
   values = solutions[numsolutions].values = (double *) calloc(length, DSIZE);
   memcpy((char *) indices, p->best_sol.xind, length * ISIZE);
   memcpy((char *) values, p->best_sol.xval, length * DSIZE);
   solutions[numsolutions].gamma = 0.0;
   solutions[numsolutions].tau = 1.0;
   solutions[numsolutions].obj[0] = p->obj[0];
   solutions[numsolutions++].obj[1] = p->obj[1];
   utopia[1] = p->obj[1];
   
   p->utopia[1] = utopia[1];
   p->utopia[0] = utopia[0];
   
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Utopia point has first  objective value %.3f\n", utopia[0]);
   printf("                 second objective value %.3f\n", utopia[1]);
   printf("***************************************************\n");
   printf("***************************************************\n\n");
   
   /* Add the first pair to the list */
   if (binary_search){
      pairs[first].gamma1 = 1.0;
      pairs[first].gamma2 = 0.0;
   }
   pairs[first].solution1 = 0;
   pairs[first].solution2 = 1;

   first = last = 0;
   numpairs = 1;

   /* Keep taking pairs off the list and processing them until there are none
      left */
   while (numpairs > 0 && numpairs < MAX_NUM_PAIRS &&
	  numsolutions < MAX_NUM_SOLUTIONS &&
	  numinfeasible < MAX_NUM_INFEASIBLE){

      if (p->par.mc_search_order == MC_LIFO){
	 solution1 = pairs[last].solution1;
	 solution2 = pairs[last].solution2;
	 cur_position = last;
	 if (--last < 0){
	    last = MAX_NUM_PAIRS - 1;
	 }
	 numpairs--;
      }else{
	 solution1 = pairs[first].solution1;
	 solution2 = pairs[first].solution2;
	 cur_position = first;
	 if (++first > MAX_NUM_PAIRS-1)
	    first = 0;
	 numpairs--;
      }

      if (binary_search){
	 gamma = (pairs[cur_position].gamma1 + pairs[cur_position].gamma2)/2;
      }else if (p->par.lp_par.mc_find_nondominated_solutions){
	 gamma = (utopia[1] - solutions[solution1].obj[1])/
	    (utopia[0] - solutions[solution2].obj[0] +
	     utopia[1] - solutions[solution1].obj[1]);
      }else{
	 slope = (solutions[solution1].obj[1] -
		  solutions[solution2].obj[1])/
	    (solutions[solution2].obj[0] -
	     solutions[solution1].obj[0]);
	 gamma = slope/(1+slope);
      }
      tau = 1 - gamma;
      
      p->par.lp_par.mc_gamma = gamma;
      p->par.lp_par.mc_tau = tau;

      /* Find upper bound */

      p->has_mc_ub = p->has_ub = FALSE;
      p->mc_ub = p->ub = MAXDOUBLE;
      if (!binary_search){
	 for (i = 0; i < numsolutions; i++){
	    if (p->par.lp_par.mc_find_nondominated_solutions){
	       ub = MAX(gamma*(solutions[i].obj[0] - utopia[0]),
			tau*(solutions[i].obj[1] - utopia[1]));
	    }else{
	       ub = gamma*solutions[i].obj[0] + tau*solutions[i].obj[1] +
		  p->par.lp_par.mc_rho * (solutions[i].obj[0] +
					  solutions[i].obj[1]);
	    }
	    if (ub + p->par.lp_par.mc_rho * (solutions[i].obj[0] +
					     solutions[i].obj[1]) < p->ub){
	       p->has_mc_ub = p->has_ub = TRUE;
	       p->ub = ub + p->par.lp_par.mc_rho *
		  (solutions[i].obj[0] + solutions[i].obj[1]) - compare_sol_tol;
	       p->obj[0] = solutions[i].obj[0];
	       p->obj[1] = solutions[i].obj[1];
	       p->mc_ub = ub;
	    }
	 }
      }
      
      printf("***************************************************\n");
      printf("***************************************************\n");
      printf("Now solving with gamma = %.6f tau = %.6f \n", gamma, tau);  
      printf("***************************************************\n");
      printf("***************************************************\n\n");
      
      p->obj[0] = p->obj[1] = 0.0;
      
      if (!p->par.lp_par.mc_find_nondominated_solutions){
	 sym_set_warm_start(p, ws);
	 for (i = 0; i < p->mip->n; i++){
	    sym_set_obj_coeff(p, i, gamma*p->mip->obj1[i]+tau*p->mip->obj2[i]
			      + p->par.lp_par.mc_rho*(p->mip->obj1[i] +
						      p->mip->obj2[i]));
	 }
	 if (termcode = sym_resolve(p) < 0){
	    sym_delete_warm_start(ws);
	    p->base->cutnum -=2;
	    p->rootdesc->uind.size--;
	    return(termcode);
	 }
      }else{
	 if (termcode = sym_solve(p) < 0){
	    p->base->cutnum -=2;
	    p->rootdesc->uind.size--;
	    return(termcode);
	 }
      }
      numprobs++;
      
      if (binary_search){
	 if (p->obj[0] - solutions[solution1].obj[0] <
	     compare_sol_tol &&
	     solutions[solution1].obj[1] - p->obj[1] <
	     compare_sol_tol){
	    if (pairs[cur_position].gamma1 - gamma >
		p->par.mc_binary_search_tolerance){
	       if (++last > MAX_NUM_PAIRS - 1)
		  last = 0;
	       pairs[last].solution1 = solution1;
	       pairs[last].solution2 = solution2;
	       pairs[last].gamma1 = gamma;
	       pairs[last].gamma2 = pairs[cur_position].gamma2;
	       numpairs++;
	    }
	    continue;
	 }
	 if (solutions[solution2].obj[0] - p->obj[0] < compare_sol_tol
	     && p->obj[1] - solutions[solution2].obj[1] <
	     compare_sol_tol){
	    if (gamma - pairs[cur_position].gamma2 >
		p->par.mc_binary_search_tolerance){
	       if (++last > MAX_NUM_PAIRS - 1)
		  last = 0;
	       pairs[last].solution1 = solution1;
	       pairs[last].solution2 = solution2;
	       pairs[last].gamma1 = pairs[cur_position].gamma1;
	       pairs[last].gamma2 = gamma;
	       numpairs++;
	    }
	    continue;
	 }
      }else{
	 if (p->obj[0] == 0.0 && p->obj[1] == 0.0){
	    numinfeasible++;
	    continue;
	 }else if (p->obj[0] - solutions[solution1].obj[0] <
		   compare_sol_tol &&
		   solutions[solution1].obj[1] - p->obj[1] <
		   compare_sol_tol){
	    numinfeasible++;
	    continue;
	 }else if (solutions[solution2].obj[0] - p->obj[0] <
		   compare_sol_tol &&
		   p->obj[1] - solutions[solution2].obj[1] <
		   compare_sol_tol){
	    numinfeasible++;
	    continue;
	 }
      }
      
      /* Insert new solution */
      numinfeasible = 0;
      if (last + 2 == MAX_NUM_PAIRS){
	 last = 0;
	 previous = MAX_NUM_PAIRS - 1;
      }else if (last + 2 == MAX_NUM_PAIRS + 1){
	 last = 1;
	 previous = 0;
      }else{
	 last += 2;
	 previous = last - 1;
      }
      if (binary_search){
	 pairs[previous].gamma1 = pairs[cur_position].gamma1;
	 pairs[previous].gamma2 = gamma;
	 pairs[last].gamma1 = gamma;
	 pairs[last].gamma2 = pairs[cur_position].gamma2;
      }
      pairs[previous].solution1 = solution1;
      pairs[previous].solution2 = solution2;
      pairs[last].solution1 = solution2;
      pairs[last].solution2 = solution2+1;
      numpairs += 2;
      for (i = numsolutions; i > solution2; i--){
	 solutions[i] = solutions[i-1];
      }
      numsolutions++;
      if (p->par.mc_search_order == MC_FIFO){
	 if (first < last){
	    for (i = first; i < last - 1; i++){
	       if (pairs[i].solution1 >= solution2){
		  pairs[i].solution1++;
	       }
	       if (pairs[i].solution2 >= solution2){
		  pairs[i].solution2++;
	       }
	    }
	 }else{
	    for (i = first; i < MAX_NUM_PAIRS - (last == 0 ? 1 : 0); i++){
	       if (pairs[i].solution1 >= solution2){
		  pairs[i].solution1++;
	       }
	       if (pairs[i].solution2 >= solution2){
		  pairs[i].solution2++;
	       }
	    }
	    for (i = 0; i < last - 1; i++){
	       if (pairs[i].solution1 >= solution2){
		  pairs[i].solution1++;
	       }
	       if (pairs[i].solution2 >= solution2){
		  pairs[i].solution2++;
	       }
	    }
	 }
      }

      length = solutions[solution2].length = p->best_sol.xlength;
      indices = solutions[solution2].indices = (int *) calloc(length, ISIZE);
      values = solutions[solution2].values = (double *) calloc(length, DSIZE);
      memcpy((char *) indices, p->best_sol.xind, length * ISIZE);
      memcpy((char *) values, p->best_sol.xval, length * DSIZE);
      solutions[solution2].gamma = gamma;
      solutions[solution2].tau = tau;
      solutions[solution2].obj[0] = p->obj[0];
      solutions[solution2].obj[1] = p->obj[1];
   }

   printf("\n********************************************************\n");

   if (numsolutions >= MAX_NUM_SOLUTIONS){
      printf("Maximum number of solutions (%i) reached\n\n",
	     MAX_NUM_SOLUTIONS);
   }

   if (numinfeasible >= MAX_NUM_INFEASIBLE){
      printf("Maximum number of infeasible subproblems (%i) reached\n\n",
	     MAX_NUM_INFEASIBLE);
   }
   
   if (numpairs >= MAX_NUM_PAIRS){
      printf("Maximum number of solution pairs (%i) reached\n\n",
	     MAX_NUM_PAIRS);
      printf("\n********************************************************\n");
      if (p->par.lp_par.mc_find_nondominated_solutions){
	 printf(  "* Found set of non-dominated solutions!!!!!!! *\n");
      }else{
	 printf(  "* Found set of supported solutions!!!!!!!     *\n");
      }
   }else{
      printf("\n********************************************************\n");
      if (p->par.lp_par.mc_find_nondominated_solutions){
	 printf(  "* Found complete set of non-dominated solutions!!!!!!! *\n");
      }else{
	 printf(  "* Found complete set of supported solutions!!!!!!!     *\n");
      }
   }
   printf(  "* Now displaying stats...                              *\n");
   printf(  "********************************************************\n\n");

   if (p->par.use_permanent_cut_pools){
      for (i = 0; i < p->par.tm_par.max_cp_num; i++){
	 p->comp_times.bc_time.cut_pool += p->cp[i]->cut_pool_time;
	 p->warm_start->stat.cuts_in_pool += p->cp[i]->cut_num;
      }
   }
   
   print_statistics(&(p->comp_times.bc_time), &(p->warm_start->stat), 0.0,
		    0.0, 0, start_time, p->mip->obj_offset,
		    p->mip->obj_sense, p->has_ub);

   printf("\nNumber of subproblems solved: %i\n", numprobs);
   printf("Number of solutions found: %i\n\n", numsolutions);
   
   printf("***************************************************\n");
   printf("***************************************************\n");
   if (p->par.lp_par.mc_find_nondominated_solutions){
      printf("Displaying non-dominated solution values and breakpoints\n");  
   }else{
      printf("Displaying supported solution values and breakpoints\n");  
   }
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   gamma0 = 1.0;
   for (i = 0; i < numsolutions - 1; i++){
      if (p->par.lp_par.mc_find_nondominated_solutions){
	 gamma1 = (utopia[1] - solutions[i].obj[1])/
	    (utopia[0] - solutions[i+1].obj[0] +
	     utopia[1] - solutions[i].obj[1]);
      }else{
	 slope = (solutions[i].obj[1] -
		  solutions[i+1].obj[1])/
	    (solutions[i+1].obj[0] -
	     solutions[i].obj[0]);
	 gamma1 = slope/(1+slope);
      }
      printf("First Objective: %.3f Second Objective: %.3f ",
	     solutions[i].obj[0], solutions[i].obj[1]);
      printf("Range: %.6f - %.6f\n", gamma1, gamma0);
      gamma0 = gamma1;
   }
   printf("First Objective: %.3f Second Objective: %.3f ",
	  solutions[i].obj[0], solutions[i].obj[1]);
   printf("Range: %.6f - %.6f\n", 0.0, gamma0);
   
   for (i = 0 ; i < numsolutions; i++){
      FREE(solutions[i].values);
      FREE(solutions[i].indices);
   }
   if (!p->par.lp_par.mc_find_nondominated_solutions){
      sym_delete_warm_start(ws);
   }
   p->base->cutnum -=2;
   p->rootdesc->uind.size--;

   return(TM_OPTIMAL_SOLUTION_FOUND);
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

int sym_load_problem_user(problem * p, int numcols, int numrows, int *start, 
			  int *index, double *value, double *collb,
			  double *colub, double *obj, char *rowsen, 
			  double *rowrhs, double *rowrng)
{
   int termcode = 0;   

   if (numcols == 0){
      printf("sym_load_problem_user():The given problem is empty!\n");
      return (0);
   }

   //Assuming all the pointers are always given NOT null, except rowrng!
   char free_range = FALSE;
   double t =0;
   int j;
   (void)used_time(&t);
   
   if (!rowrng){
      rowrng = (double*)calloc(numrows, DSIZE);
      free_range = TRUE;
   }
   p->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
 
   p->mip->m  = numrows;
   p->mip->n  = numcols;
   p->mip->nz = start[numcols];
   
   p->mip->obj    = (double *) malloc(DSIZE * numcols);
   p->mip->rhs    = (double *) malloc(DSIZE * numrows);
   p->mip->sense  = (char *)   malloc(CSIZE * numrows);
   p->mip->rngval = (double *) malloc(DSIZE * numrows);
   p->mip->ub     = (double *) malloc(DSIZE * numcols);
   p->mip->lb     = (double *) malloc(DSIZE * numcols);
   p->mip->is_int = (char *)   calloc(CSIZE, numcols);

   memcpy(p->mip->obj, obj, DSIZE * numcols); 
   memcpy(p->mip->rhs, rowrhs, DSIZE * numrows); 
   memcpy(p->mip->sense, rowsen, CSIZE * numrows); 
   memcpy(p->mip->rngval, rowrng, DSIZE * numrows); 	  
   memcpy(p->mip->ub, colub, DSIZE * numcols); 
   memcpy(p->mip->lb, collb, DSIZE * numcols); 

   //user defined matind, matval, matbeg--fill as column ordered
   
   p->mip->matbeg = (int *) malloc(ISIZE * (numcols + 1));
   p->mip->matval = (double *) malloc(DSIZE*start[numcols]);
   p->mip->matind = (int *)    malloc(ISIZE*start[numcols]);
   
   memcpy(p->mip->matbeg, start, ISIZE * (numcols + 1));
   memcpy(p->mip->matval, value, DSIZE * 
	  start[numcols]);  
   memcpy(p->mip->matind, index, ISIZE * 
	  start[numcols]);  

   /* Start up the graphics window*/
#ifndef WIN32
   CALL_WRAPPER_FUNCTION( init_draw_graph_u(p) );   
#endif

   /*------------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*------------------------------------------------------------------------*/

   CALL_WRAPPER_FUNCTION(initialize_root_node_u(p) ); 
   
   p->comp_times.readtime = used_time(&t);
 
   p->termcode = TM_NO_SOLUTION;

   if (free_range){
      FREE(rowrng);
   }
   
   return termcode;

}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_abandoned(problem *p)
{

   if (p){
      switch(p->termcode){
	 
       case SOMETHING_DIED:
       case TM_ERROR__NUMERICAL_INSTABILITY:    
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_abandoned():The problem is empty!\n");
      return FALSE;
   }
 
   return FALSE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_proven_optimal(problem *p)
{

   if (p){
      switch(p->termcode){
	 
       case TM_OPTIMAL_SOLUTION_FOUND:
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_proven_optimal():The problem is empty!\n");
      return FALSE;
   }
   
   return FALSE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_proven_primal_infeasible(problem *p)
{

   if (p){
      switch(p->termcode){
	 
       case TM_NO_SOLUTION:
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_proven_primal_infeasible():The problem is empty!\n");
      return FALSE;
   }
   
   return FALSE;
}

/*===========================================================================*/
/*===========================================================================*/
int sym_is_primal_objective_limit_reached(problem *p)
{
   if (p){
      switch(p->termcode){
	 
       case TM_TARGET_GAP_ACHIEVED:
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_primal_objective_limit_reached():");
      printf("The problem is empty!\n");
      return FALSE;
   }
   
   return FALSE;
}
/*===========================================================================*/
/*===========================================================================*/

int sym_is_iteration_limit_reached(problem *p)
{
      if (p){
      switch(p->termcode){
	 
       case TM_NODE_LIMIT_EXCEEDED:
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_iteration_limit_reached():The problem is empty!\n");
      return FALSE;
   }
   
   return FALSE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_time_limit_reached(problem *p)
{
      if (p){
      switch(p->termcode){
	 
       case TM_TIME_LIMIT_EXCEEDED:
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_time_limit_reached():The problem is empty!\n");
      return FALSE;
   }
   
   return FALSE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_target_gap_achieved(problem *p)
{
   if (p){
      switch(p->termcode){
	 
       case TM_TARGET_GAP_ACHIEVED:
	  return TRUE;
       default:
	  break;
      }
   }
   else{
      printf("sym_is_time_limit_reached():The problem is empty!\n");
      return FALSE;
   }
   
   return FALSE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_num_cols(problem *p)
{
   if (!p->mip){
      printf("sym_get_num_cols():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->n;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_num_rows(problem *p)
{
   if (!p->mip){
      printf("sym_get_num_rows():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->m;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_num_elements(problem *p)
{
   if (!p->mip){
      printf("sym_get_num_elements():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->nz;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_col_lower(problem *p)
{
   if (!p->mip){
      printf("sym_get_col_lower():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->lb;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_col_upper(problem *p)
{
   if (!p->mip){
      printf("sym_get_col_upper():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->ub;
}

/*===========================================================================*/
/*===========================================================================*/

char *sym_get_row_sense(problem *p)
{
   if (!p->mip){
      printf("sym_get_row_sense():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->sense;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_rhs(problem *p)
{
   if (!p->mip){
      printf("sym_get_rhs():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->rhs;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_row_range(problem *p)
{
   if (!p->mip){
      printf("sym_get_row_range():The problem description is empty!\n");
      return (0);
   }
   
   return p->mip->rngval;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_row_lower(problem *p)
{
   if (!p->mip){
      printf("sym_get_row_lower():The problem description is empty!\n");
      return (0);
   }
   
   double * lower = (double *)malloc(DSIZE*p->mip->m);
   double upper;
   int i;
   double rhs, range, inf = INFINITY;
   char sense;

   for ( i = p->mip->m - 1; i >= 0; --i )
      {

	 rhs   = p->mip->rhs[i];
	 range = p->mip->rngval[i];
	 sense = p->mip->sense[i];
	 
	 switch (sense) {
	  case 'E':
	     lower[i] = upper = rhs;
	     break;
	  case 'L':
	     lower[i] = -inf;
	     upper = rhs;
	     break;
	  case 'G':
	     lower[i] = rhs;
	     upper = inf;
	     break;
	  case 'R':
	     lower[i] = rhs - range;
	     upper = rhs;
	     break;
	  case 'N':
	     lower[i] = -inf;
	     upper = inf;
	     break;
	 }
      }

   return lower;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_row_upper(problem *p)
{
   if (!p->mip){
      printf("sym_get_row_upper():The problem description is empty!\n");
      return (0);
   }

   double * upper = (double *)malloc(DSIZE*p->mip->m);
   double lower;
   int i;
   double rhs, range, inf = INFINITY;
   char sense;
   
   for ( i = p->mip->m - 1; i >= 0; --i )
      {

	 rhs   = p->mip->rhs[i];
	 range = p->mip->rngval[i];
	 sense = p->mip->sense[i];
	 
	 switch (sense) {
	  case 'E':
	     lower = upper[i] = rhs;
	     break;
	  case 'L':
	     lower = -inf;
	     upper[i] = rhs;
	     break;
	  case 'G':
	     lower = rhs;
	     upper[i] = inf;
	     break;
	  case 'R':
	     lower = rhs - range;
	     upper[i] = rhs;
	     break;
	  case 'N':
	     lower = -inf;
	     upper[i] = inf;
	     break;
	 }
      }

   return upper;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_obj_coeff(problem *p)
{
   if (!p->mip){
      printf("sym_get_obj_coeff():The problem description is empty!\n");
      return (0);
   }

   return p->mip->obj;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_obj_sense(problem *p)
{
   if (!p->mip){
      printf("sym_get_obj_sense():The problem description is empty!\n");
      return (0);
   }

   if (p->mip->obj_sense == MINIMIZE)
      return 1;
   else if (p->mip->obj_sense == MAXIMIZE)
      return -1;
   else
      return 1; 
}   

/*===========================================================================*/
/*===========================================================================*/

int sym_is_continuous(problem *p, int index)
{
   if (!p->mip){
      printf("sym_is_continuous():The problem description is empty!\n");
      return (-1);
   }

   return(p->mip->is_int[index]);
}
/*===========================================================================*/
/*===========================================================================*/

int sym_is_binary(problem *p, int index)
{
   if (!p->mip){
      printf("sym_is_binary():The problem description is empty!\n");
      return (-1);
   }

   if (p->mip->is_int[index] && p->mip->lb[index] == 0.0 &&
      p->mip->ub[index] == 1.0) {
      return TRUE;
   }
   else{
      return FALSE;
   }
}
/*===========================================================================*/
/*===========================================================================*/

int sym_is_integer(problem *p, int index)
{
   if (!p->mip){
      printf("sym_is_integer():The problem description is empty!\n");
      return (-1);
   }

   return(p->mip->is_int[index]);
}

/*===========================================================================*/
/*===========================================================================*/

double sym_get_infinity()
{
   return INFINITY;
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_col_solution(problem *p)
{
   int i;
   double * colSol;
   lp_sol sol;

   sol = p->best_sol;

   if (!sol.xlength){
      printf("\nNo Solution Found\n\n");
      return 0;
   }
   else{
      colSol = (double*)calloc(DSIZE, p->mip->n);
      for( i = 0; i<sol.xlength; i++){
	 colSol[sol.xind[i]] = sol.xval[i];
      }
      return (colSol);
   }
}

/*===========================================================================*/
/*===========================================================================*/

double *sym_get_row_activity(problem *p)
{
   double * rowAct;
   double * colSol;  
   int i, j;

   const int * matBeg;
   const double * matVal;
   const int * matInd;

   colSol = sym_get_col_solution(p);

   if (!p->mip){
      printf("sym_get_row_activity():The problem description is empty!\n");
      return (0);
   }

   rowAct = (double*)calloc(DSIZE, p->mip->m);

   if (colSol){

      matBeg = p->mip->matbeg;
      matVal = p->mip->matval;
      matInd = p->mip->matind;

      for(i = 0; i<p->mip->n; i++){
	 for(j = matBeg[i]; j<matBeg[i+1]; j++){
	    rowAct[matInd[j]] += matVal[j] * colSol[i];
	 }
      }
      return rowAct;
   }
   else
      return 0;         
}

/*===========================================================================*/
/*===========================================================================*/

double sym_get_obj_val(problem *p)
{
   if (p->best_sol.xlength)
      return (p->best_sol.objval);
   else{
      printf("\nNo Solution Found\n\n");
      return INFINITY;
   }
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_iteration_count(problem *p)
{
   if (!p->warm_start){
      printf("sym_get_iteration_count():");
      printf("The problem warm start description is empty!\n");
      return (0);
   }

   return p->warm_start->stat.analyzed;  
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_obj_coeff(problem *p, int index, double value)
{

   int i;

   if (!p->mip){
      printf("sym_set_obj_coeff():The problem description is empty!\n");
      return FALSE;
   }
   
   p->mip->obj[index] = value;

   if (p->mip->change_num){
      for(i = p->mip->change_num - 1 ; i >=0 ; i--){
	 if (p->mip->change_type[i] == OBJ_COEFF_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 p->mip->change_type[p->mip->change_num++] = OBJ_COEFF_CHANGED;
      }
   }
   else{
      p->mip->change_type[p->mip->change_num++] = OBJ_COEFF_CHANGED;
   }
   
   
   
   return TRUE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_obj2_coeff(problem *p, int index, double value)
{

   int i;

   if (!p->mip){
      printf("sym_set_obj_coeff():The problem description is empty!\n");
      return FALSE;
   }
   
   p->mip->obj2[index] = value;

   return TRUE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_lower(problem *p, int index, double value)
{
   if (!p->mip){
      printf("sym_set_col_lower():The problem description is empty!\n");
      return FALSE;
   }

   p->mip->lb[index] = value;
   return TRUE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_upper(problem *p, int index, double value)
{
   if (!p->mip){
      printf("sym_set_col_upper():The problem description is empty!\n");
      return FALSE;
   }

   p->mip->ub[index] = value;
   return TRUE;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_row_lower(problem *p, int index, double value)
{
   double rhs, range, lower, upper, inf = INFINITY;
   char   sense;
   int i;

   if (!p->mip){
      printf("sym_set_row_lower():The problem description is empty!\n");
      return FALSE;
   }

   rhs   = p->mip->rhs[index];
   range = p->mip->rngval[index];
   sense = p->mip->sense[index];
   
   switch (sense) {
    case 'E':
       lower = upper = rhs;
       break;
    case 'L':
       lower = -inf;
       upper = rhs;
       break;
    case 'G':
       lower = rhs;
       upper = inf;
       break;
    case 'R':
       lower = rhs - range;
       upper = rhs;
       break;
    case 'N':
       lower = -inf;
       upper = inf;
       break;
   }

   if ( lower != value ) {
      lower = value;
      range = 0.0;
      if (lower > -inf) {
	 if (upper < inf) {
	    rhs = upper;
	    if (upper==lower) {
	       sense = 'E';
	    } else {
	       sense = 'R';
	       range = upper - lower;
	    }
	 } else {
	    sense = 'G';
	    rhs = lower;
	 }
      } else {
	 if (upper < inf) {
	    sense = 'L';
	    rhs = upper;
	 } else {
	    sense = 'N';
	    rhs = 0.0;
	 }
      }    

      p->mip->sense[index] = sense;   
      p->mip->rhs[index] = rhs;
      p->mip->rngval[index] = range;
   }

   if (p->mip->change_num){
      for(i = p->mip->change_num - 1 ; i >=0 ; i--){
	 if (p->mip->change_type[i] == RHS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 p->mip->change_type[p->mip->change_num] = RHS_CHANGED;
	 p->mip->change_num++;
      }
   }
   else{
      p->mip->change_type[p->mip->change_num] = RHS_CHANGED;
      p->mip->change_num++;
   }
   
   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_row_upper(problem *p, int index, double value)
{
   double rhs, range, lower, upper, inf = INFINITY;
   char   sense;
   int i;

   if (!p->mip){
      printf("sym_set_row_upper():The problem description is empty!\n");
      return FALSE;
   }

   rhs   = p->mip->rhs[index];
   range = p->mip->rngval[index];
   sense = p->mip->sense[index];

   switch (sense) {
    case 'E':
       lower = upper = rhs;
       break;
    case 'L':
       lower = -inf;
       upper = rhs;
       break;
    case 'G':
       lower = rhs;
       upper = inf;
       break;
    case 'R':
       lower = rhs - range;
       upper = rhs;
       break;
    case 'N':
       lower = -inf;
       upper = inf;
       break;
   }

   /*   convertSenseToBound( sense, rhs, range, lower, upper );*/
   
   if ( upper != value ) {
      upper = value;
      /* convertBountToSense */
      range = 0.0;
      if (lower > -inf) {
	 if (upper < inf) {
	    rhs = upper;
	    if (upper==lower) {
	       sense = 'E';
	    } else {
	       sense = 'R';
	       range = upper - lower;
	    }
	 } else {
	    sense = 'G';
	    rhs = lower;
	 }
      } else {
	 if (upper < inf) {
	    sense = 'L';
	    rhs = upper;
	 } else {
	    sense = 'N';
	    rhs = 0.0;
	 }
      } 
      
      p->mip->sense[index] = sense;   
      p->mip->rhs[index] = rhs;
      p->mip->rngval[index] = range;
   }

   if (p->mip->change_num){
      for(i = p->mip->change_num - 1 ; i >=0 ; i--){
	 if (p->mip->change_type[i] == RHS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 p->mip->change_type[p->mip->change_num] = RHS_CHANGED;
	 p->mip->change_num++;

      }
   }
   else{
      p->mip->change_type[p->mip->change_num] = RHS_CHANGED;
      p->mip->change_num++;

   }

   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_row_type(problem *p, int index, char rowsense, double rowrhs, 
		      double rowrng)
{

   int i;

   if (!p->mip){
      printf("sym_set_row_type():The problem description is empty!\n");
      return FALSE;
   }

   //FIX_ME! what if sense is changed? 

   p->mip->sense[index] = rowsense;   
   p->mip->rhs[index] = rowrhs;
   p->mip->rngval[index] = rowrng;


   if (p->mip->change_num){
      for(i = p->mip->change_num - 1 ; i >=0 ; i--){
	 if (p->mip->change_type[i] == RHS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 p->mip->change_type[p->mip->change_num] = RHS_CHANGED;
	 p->mip->change_num++;

      }
   }
   else{
      p->mip->change_type[p->mip->change_num] = RHS_CHANGED;
      p->mip->change_num++;

   }
   
   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_obj_sense(problem *p, int sense)
{
   if (!p->mip){
      printf("sym_set_obj_type():The problem description is empty!\n");
      return FALSE;
   }

   if (sense == 1){
      p->mip->obj_sense = MINIMIZE;
   }
   else if (sense==-1){
      p->mip->obj_sense = MAXIMIZE;
   }
   else{
      /* assume it to be min problem */
      p->mip->obj_sense = MINIMIZE;
   }
   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_solution(problem *p, double * colsol)
{
   int i, j, nz =0, *matBeg, *matInd;
   double value, *rowAct, *matVal; 
   char feasible;
   lp_sol * sol;

   if (!p->mip){
      printf("sym_set_col_solution():The problem description is empty!\n");
      return FALSE;
   }

   /* Feasibility Check*/

   /* step 1. check for bounds and integrality */   
   for (i = p->mip->n - 1; i >= 0; i--){
      if (colsol[i] < p->mip->lb[i] || colsol[i] > p->mip->ub[i])
	 break;
      if (colsol[i] !=0.0)
	 nz++;
      if (!p->mip->is_int[i])
	 continue; /* Not an integer variable */
      value = colsol[i];
      if (colsol[i] > p->mip->lb[i] && colsol[i] < p->mip->ub[i]
	  && colsol[i]-floor(colsol[i]) > p->par.lp_par.granularity &&
	  ceil(colsol[i])-colsol[i] > p->par.lp_par.granularity){
	 break;   //FIX_ME, can we use granularity here?
      }
   }

   feasible = i < 0 ? true : false;
   
   /* step 2. check for the constraint matrix */
   
   if (feasible){      
      rowAct = (double*) calloc(p->mip->m, DSIZE);
      matBeg = p->mip->matbeg;
      matVal = p->mip->matval;
      matInd = p->mip->matind;

      for(i = 0; i<p->mip->n; i++){
	 for(j = matBeg[i]; j<matBeg[i+1]; j++){
	    rowAct[matInd[j]] += matVal[j] * colsol[i];
	 }
      }	 
 
      for(i = 0; i<p->mip->m; i++){
	 switch(p->mip->sense[i]){
	  case 'L': 
	     if (rowAct[i] > p->mip->rhs[i])
		feasible = FALSE;
	     break;
	  case 'G':
	     if (rowAct[i] < p->mip->rhs[i])
		feasible = FALSE;
	     break;
	  case 'E':
	     if (rowAct[i] != p->mip->rhs[i])
		feasible = FALSE;
	     break;
	  case 'R':
	     if (rowAct[i] > p->mip->rhs[i] || 
		rowAct[i] < p->mip->rhs[i] - p->mip->rngval[i])
		feasible = FALSE;
	     break;
	  case 'N':
	  default:
	     break;
	 }
	 
	 if (!feasible) 
	    break;
      }
   }

   if (feasible){
      /* now, it is feasible, set the best_sol to colsol */
      //FIX_ME
      
      sol = (lp_sol*) calloc(1,sizeof(lp_sol));      
      sol->xlength = nz;
      sol->xind = (int*)malloc(ISIZE*nz);
      sol->xval = (double*)calloc(nz,DSIZE);
      
      for(i = 0, j =0 ; i<p->mip->n; i++){
	 if (colsol[i] != 0.0){
	    sol->xind[j] = i;
	    sol->xval[j] = colsol[i];
	    sol->objval += colsol[i] * p->mip->obj[i];	   
	    j++;
	 }      
      }

      if (p->has_ub_estimate){
	 if (p->ub_estimate > sol->objval)
	    p->ub_estimate = sol->objval; //no need for this, I guess.
      }
      else{
	 p->has_ub_estimate = TRUE;
	 p->ub_estimate = sol->objval; //no need for this, I guess.
      }

      if (p->has_ub){
	 if (p->ub > sol->objval)
	    p->ub = sol->objval;
      }
      else{
	 p->has_ub = TRUE;
	 p->ub = sol->objval;
      }

      if (p->best_sol.xlength){
	 if (p->best_sol.objval > sol->objval){
	    FREE(p->best_sol.xind);
	    FREE(p->best_sol.xval);
	    p->best_sol = *sol;
	 }
	 else{
	    printf("sym_set_col_solution(): The col.solution was not set:\n");
	    printf("the current solution is better!\n");
	 }
      }
      else{
	 p->best_sol = *sol;           
      }      
   }
   else{
      printf("sym_set_col_solution(): The col.solution was not set:\n");
      printf("the given solution is not feasible!\n");
   }  

   if (rowAct)
      FREE(rowAct);
   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_continuous(problem *p, int index)
{
   if (!p->mip){
      printf("sym_set_continuous():The problem description is empty!\n");
      return FALSE;
   }

   p->mip->is_int[index] = FALSE;
   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_integer(problem *p, int index)
{
   if (!p->mip){
      printf("sym_set_integer():The problem description is empty!\n");
      return FALSE;
   }

   p->mip->is_int[index] = TRUE;
   return TRUE;      
}
/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_names(problem * p, char **colname)
{
   int j;

   if (!p->mip){
      printf("sym_set_col_names():The problem description is empty!\n");
      return FALSE;
   }

   p->mip->colname = (char **)  malloc(sizeof(char *) * p->mip->n);
   
   for (j = 0; j < p->mip->n; j++){
      /* FIX_ME! Restricting col names to 20 characters! */
      p->mip->colname[j] = (char *) malloc(CSIZE * 20); 
      strncpy(p->mip->colname[j], colname[j], 20);
      p->mip->colname[j][19] = 0;  /* ??? */
   }      
   return TRUE;      
}
/*===========================================================================*/
/*===========================================================================*/

int sym_add_col(problem *p, int num_elements, int *indices, 
			double *elements, double collb, double colub,
			double obj, char *name)
{

   int i, n, nz, *matBeg, *matInd;
   double *matVal, *colLb, *colUb, *objN;
   char *isInt, **colName;

   if (!p->mip){
      printf("sym_add_col():The problem description is empty!\n");
      return FALSE;
   }

   n = p->mip->n;
   nz = p->mip->nz;

   matBeg = (int*) malloc(ISIZE*(n+2));
   matInd = (int*) malloc(ISIZE*(nz+num_elements));
   matVal = (double*) malloc(DSIZE*(nz+num_elements));
   colLb = (double*) malloc(DSIZE*(n+1));
   colUb = (double*) malloc(DSIZE*(n+1));
   objN = (double*) malloc(DSIZE*(n+1));
   isInt = (char*) calloc(CSIZE, (n+1));
   colName = (char**) malloc(sizeof(char*)*(n+1));

   memcpy(matBeg, p->mip->matbeg, ISIZE*(n+1));
   memcpy(matInd, p->mip->matind, ISIZE*nz);
   memcpy(matVal, p->mip->matval, DSIZE*nz);
   memcpy(colLb, p->mip->lb, DSIZE*n);
   memcpy(colUb, p->mip->ub, DSIZE*n);
   memcpy(objN, p->mip->obj, DSIZE*n);
   memcpy(isInt, p->mip->is_int, CSIZE*n);

   for (i = 0; i < n; i++){
      colName[i] = (char *) malloc(CSIZE * 20); 
      strncpy(colName[i], p->mip->colname[i],20);
      colName[i][19] = 0;
   }

   matBeg[n+1] = matBeg[n] + num_elements;
   memcpy(matInd + nz, indices, ISIZE*num_elements);
   memcpy(matVal + nz, elements, DSIZE*num_elements); 
   colLb[n] = collb;
   colUb[n] = colub;
   objN[n] = obj;
   
   colName[n] = (char *) malloc(CSIZE * 20); 
   strncpy(colName[n], name,20);
   colName[n][19] = 0;

   FREE(p->mip->matbeg);
   FREE(p->mip->matind);
   FREE(p->mip->matval);
   FREE(p->mip->lb);
   FREE(p->mip->ub);
   FREE(p->mip->obj);
   FREE(p->mip->is_int);
   for (i = 0; i < n; i++){
      FREE(p->mip->colname[i]);
   }
   FREE(p->mip->colname);

   p->mip->n = n+1;
   p->mip->nz = nz + num_elements;
   p->mip->matbeg = matBeg;
   p->mip->matind = matInd;
   p->mip->matval = matVal;
   p->mip->lb =  colLb;
   p->mip->ub = colUb;
   p->mip->obj = objN;
   p->mip->is_int = isInt;   
   p->mip->colname = colName;

   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/
int sym_add_row(problem *p, int num_elements, int *indices, 
			double *elements, char rowsen, double rowrhs,
			double rowrng)
{
   int i, j, m, n, nz, *matBeg, *matInd, *lengths;
   double *matVal, *rhs, *range;
   char *sense;
   
   if (!p->mip){
      printf("sym_add_row():The problem description is empty!\n");
      return FALSE;
   }

   /*FIX_ME! Put sym_add_row(numelem, indices, elements, rowsen, 
     rowrhs, rowrng) here! */

   m = p->mip->m;
   n = p->mip->n;
   nz = p->mip->nz;

   matBeg = (int*) calloc (n, ISIZE);
   matInd = (int*) malloc(ISIZE*(nz+num_elements));
   matVal = (double*) malloc(DSIZE*(nz+num_elements));  
   lengths = (int*) malloc (ISIZE*n);
   sense = (char*) malloc(CSIZE*(m+1));
   rhs = (double*) malloc(DSIZE*(m+1));
   range = (double*) malloc(DSIZE*(m+1));
   
   memcpy(sense, p->mip->sense, CSIZE*m);
   memcpy(range, p->mip->rngval, DSIZE*m);
   memcpy(rhs, p->mip->rhs, DSIZE*m);

   for(i = 0; i<n; i++){     
      lengths[i] = p->mip->matbeg[i+1] - p->mip->matbeg[i];
   }

   for(i = 0; i<num_elements; i++){
      lengths[indices[i]]++;
   }

   for(i = 0, j = 0; i<n; i++){
      matBeg[i+1] = matBeg[i] + lengths[i];
      memcpy(matInd + matBeg[i], p->mip->matind + p->mip->matbeg[i], 
	     ISIZE * (p->mip->matbeg[i+1]-p->mip->matbeg[i])); 
      memcpy(matVal + matBeg[i], p->mip->matval + p->mip->matbeg[i], 
	     DSIZE * (p->mip->matbeg[i+1]-p->mip->matbeg[i])); 
      if (indices[j] == i){
	 matInd[matBeg[i+1]-1] = m;
	 matVal[matBeg[i+1]-1] = elements[j];
	 j++;
      }
   }

   sense[m] = rowsen;
   rhs[m] = rowrhs;
   range[m] = rowrng;
   
   /*can use FREE_mip_desc???*/
   FREE(p->mip->matbeg);
   FREE(p->mip->matind);
   FREE(p->mip->matval);
   FREE(p->mip->sense);
   FREE(p->mip->rhs);
   FREE(p->mip->rngval);
   FREE(lengths);

   p->mip->m = m+1;
   p->mip->nz = nz + num_elements;
   p->mip->matbeg = matBeg;
   p->mip->matind = matInd;
   p->mip->matval = matVal;
   p->mip->sense =  sense;
   p->mip->rhs = rhs;
   p->mip->rngval = range;

   return TRUE;      
}
/*===========================================================================*/
/*===========================================================================*/

int sym_delete_cols(problem *p, int num, int * indices)
{

   int i, j, k, n, nz, temp = 0, numElements = 0, *matBeg, *matInd, *lengths;
   double *matVal, *colLb, *colUb, *objN;
   char *isInt;

   if (!p->mip){
      printf("sym_delete_cols():The problem description is empty!\n");
      return FALSE;
   }


   /* sort the indices in case they are not given sorted! */

   for(i = 0; i < num; i++){
     temp = indices[i];
     for(j = i-1; j >= 0 && temp < indices[j]; j--){
	indices[j+1] = indices[j];
     }
     indices[j+1] = temp;
  }

   n = p->mip->n;
   nz = p->mip->nz;

   if (num>n){
      printf("sym_delete_cols():The number of cols to be deleted exceeds the");
      printf("real size!\n");
      return FALSE;
   }

   lengths = (int*) malloc (ISIZE*n);

   for(i = 0; i<n; i++){     
      lengths[i] = p->mip->matbeg[i+1] - p->mip->matbeg[i];
   }

   for( i = 0; i<num; i++){
      if (indices[i]<n){
	 numElements += lengths[indices[i]];
      }
      else{
	 /*FIX_ME*/
	 printf("sym_delete_cols(): Column index is out of range!\n");
	 return FALSE;
      }
   }

   matBeg = (int*) malloc(ISIZE*(n-num));
   matInd = (int*) malloc(ISIZE*(nz-numElements));
   matVal = (double*) malloc(DSIZE*(nz-numElements));
   colLb = (double*) malloc(DSIZE*(n-num));
   colUb = (double*) malloc(DSIZE*(n-num));
   objN = (double*) malloc(DSIZE*(n-num));
   isInt = (char*) calloc(CSIZE, (n-num));

   matBeg[0] = 0;
   for(i = 0, j = 0, k = 0; i < n; i++){
      if (indices[j] == i){
	 j++;
      }
      else{
	 matBeg[k+1] = matBeg[k] + lengths[i];
	 memcpy(matInd + matBeg[k], p->mip->matind + p->mip->matbeg[i], 
		ISIZE * lengths[i]); 
	 memcpy(matVal + matBeg[k], p->mip->matval + p->mip->matbeg[i], 
		DSIZE * lengths[i]); 
	 colLb[k] = p->mip->lb[i];
	 colUb[k] = p->mip->ub[i];
	 objN[k] = p->mip->obj[i];
	 isInt[k] = p->mip->is_int[i];
	 k++;
      }
   }

   FREE(p->mip->matbeg);
   FREE(p->mip->matind);
   FREE(p->mip->matval);
   FREE(p->mip->lb);
   FREE(p->mip->ub);
   FREE(p->mip->obj);
   FREE(p->mip->is_int);
   FREE(lengths);

   p->mip->n = n-num;
   p->mip->nz = nz - numElements;
   p->mip->matbeg = matBeg;
   p->mip->matind = matInd;
   p->mip->matval = matVal;
   p->mip->lb =  colLb;
   p->mip->ub = colUb;
   p->mip->obj = objN;
   p->mip->is_int = isInt;   

   return TRUE;      

}

/*===========================================================================*/
/*===========================================================================*/
int sym_delete_rows(problem *p, int num, int * indices)
{

   int i, j = 0, k, n, m, nz, numElements = 0, numRows = 0, *matBeg, *matInd; 
   int deletedRows, deleted;
   double *matVal, *rhs, *range;
   char *sense;

   if (!p->mip){
      printf("sym_delete_rows():The problem description is empty!\n");
      return FALSE;
   }

   n = p->mip->n;
   m = p->mip->m;
   nz = p->mip->nz;

   if (num>m){
      printf("sym_delete_rows():The number of rows to be deleted exceeds the");
      printf("real row number!\n");
      return FALSE;
   }
   
   matBeg = p->mip->matbeg;
   matInd = p->mip->matind;
   matVal = p->mip->matval;
   sense = p->mip->sense;
   rhs = p->mip->rhs;
   range = p->mip->rngval;

   /* assuming that the rowIndices may not be given in order! */
   /* FIX_ME, ask to Prof. Ralphs. */
   for(i = 0; i<n; i++){
      for(; j<matBeg[i+1]; j++){
	 for( k = 0, deleted = 0, deletedRows = 0; k<num; k++){
	    if (matInd[j] == indices[k]){	    
	       deleted = 1;
	       break;
	    }
	    if (matInd[j] > indices[k]){
	       deletedRows++;
	    }
	 }
	 if (!deleted){
	    matInd[numElements] = matInd[j] - deletedRows;
	    matVal[numElements] = matVal[j];
	    numElements++;
	 }
      }
      j = matBeg[i+1];
      matBeg[i+1] = numElements;      
   }

   for(i = 0; i<m; i++){
      for( k = 0, deleted = 0; k<num; k++){
	 if (i == k){	    
	    deleted = 1;
	    break;
	 }
      }
      if (!deleted){
	 sense[numRows] = sense[i];
	 rhs[numRows] = rhs[i];
	 range[numRows] = range[i];
	 numRows++;
      }
   }

   if (numRows != m - num){
      printf("sym_delete_rows(): Unknown error!\n");
      return FALSE;
   }


   p->mip->m  = numRows;
   p->mip->nz = numElements;
   
   p->mip->rhs    = (double *) malloc(DSIZE * numRows);
   p->mip->sense  = (char *)   malloc(CSIZE * numRows);
   p->mip->rngval = (double *) malloc(DSIZE * numRows);
   
   p->mip->matval = (double *) malloc(DSIZE*matBeg[n]);
   p->mip->matind = (int *)    malloc(ISIZE*matBeg[n]);


   memcpy(p->mip->rhs, rhs, DSIZE*numRows);
   memcpy(p->mip->rngval, range, DSIZE*numRows);
   memcpy(p->mip->sense, sense, CSIZE*numRows);
   
   memcpy(p->mip->matval, matVal, DSIZE * matBeg[n]);  
   memcpy(p->mip->matind, matInd, ISIZE * matBeg[n]);


   FREE(matVal);
   FREE(matInd);
   FREE(sense);
   FREE(rhs);
   FREE(range);

   return TRUE;      
}

/*===========================================================================*/
/*===========================================================================*/

/*===========================================================================*/
/*===========================================================================*/

int sym_write_warm_start_desc(warm_start_desc *ws, char *file)
{
 
   FILE * f = NULL;
   int i, j, temp;
   cut_data ** cuts;
   problem_stat stat;
   node_times compT;

   f = fopen(file, "w");

   if (!ws){
      printf("There is no loaded warmStart to write!\n");
      fclose(f);
      return FALSE;
   }
   else{
      fprintf(f, "########################################################\n");
      fprintf(f, " BOUND INFO \n");
      fprintf(f, "########################################################\n");
      fprintf(f, " PHASE      : %i\n", ws->phase);
      fprintf(f, " LB         : %.4f\n", ws->lb);
      fprintf(f, " HAS_UB     : %i\n", (int)ws->has_ub);
      fprintf(f, " UB         : %.4f\n\n", ws->ub);

      fprintf(f, "########################################################\n");
      fprintf(f, " CUT INFO \n");
      fprintf(f, "########################################################\n");
      fprintf(f, " CUT_NUM             : %i\n", ws->cut_num);
      fprintf(f, " ALLOCATED_CUT_NUM   : %i\n\n", 
	      ws->allocated_cut_num);

      //FIX_ME! WHAT TYPE A CUT CAN BE OTHER THAN EXPLICIT_ROW

      cuts = ws->cuts;
      
      for(i=0; i<ws->cut_num; i++){
	 fprintf(f, " CUT %i : \n",i);
	 fprintf(f, " SIZE        : %i \n",cuts[i]->size);
	 fprintf(f, " ELEMENTS    : ");
	 for(j=0; j<cuts[i]->size; j++){
	    fprintf(f," %i",(int)cuts[i]->coef[j]);
	 }
	 fprintf(f, "\n");
	 fprintf(f, " RHS         : %.4f \n",cuts[i]->rhs);
	 fprintf(f, " RANGE       : %.4f \n",cuts[i]->range);
	 fprintf(f, " TYPE        : %i \n",(int)cuts[i]->type);
	 fprintf(f, " SENSE       : %c \n",cuts[i]->sense);
	 fprintf(f, " DELETABLE   : %i \n",(int)cuts[i]->deletable);
	 fprintf(f, " BRANCH      : %i \n",(int)cuts[i]->branch);
	 fprintf(f, " NAME        : %i \n\n",cuts[i]->name);
      }

      fprintf(f, "########################################################\n");
      fprintf(f, " PROBLEM STATISTICS \n");
      fprintf(f, "########################################################\n");

      stat= ws->stat;


      fprintf(f," ROOT_LB                : %.4f\n", stat.root_lb);
      fprintf(f," CUTS_IN_POOL           : %i\n", stat.cuts_in_pool);
      fprintf(f," MAXIMIM_DEPTH          : %i\n", stat.max_depth);
      fprintf(f," DIVING_CHAINS          : %i\n", stat.chains);
      fprintf(f," DIVING_STOPS           : %i\n", stat.diving_halts);
      fprintf(f," TREE_SIZE              : %i\n", stat.tree_size);
      fprintf(f," CREATED_NODES          : %i\n", stat.created);
      fprintf(f," ANALYZED_NODES         : %i\n", stat.analyzed);
      fprintf(f," LEAVES_BEFORE_TRIMMING : %i\n", stat.leaves_before_trimming);
      fprintf(f," LEAVES_BEFORE_TRIMMING : %i\n", stat.leaves_after_trimming);
      fprintf(f," NOT_FIXED_VARIABLE_NUM : %i\n", stat.vars_not_priced);
      fprintf(f," NF_STATUS_OF_ROOT      : %i\n\n", (int)stat.nf_status);
     
      fprintf(f, "########################################################\n");
      fprintf(f, " COMPUTATION TIMES \n");
      fprintf(f, "########################################################\n");

      compT = ws->comp_times;

      fprintf(f," COMMUNICATION       : %.4f\n",compT.communication);
      fprintf(f," LP                  : %.4f\n",compT.lp);
      fprintf(f," SEPARATION          : %.4f\n",compT.separation);
      fprintf(f," FIXING              : %.4f\n",compT.fixing);
      fprintf(f," PRICING             : %.4f\n",compT.pricing);
      fprintf(f," STRONG_BRANCHING    : %.4f\n",compT.strong_branching);
      fprintf(f," WALL_CLOCK_LP       : %.4f\n",compT.wall_clock_lp);
      fprintf(f," RAMP_UP_TM          : %.4f\n",compT.ramp_up_tm);
      fprintf(f," RAMP_UP_LP          : %.4f\n",compT.ramp_up_lp);
      fprintf(f," RAMP_DOWN_TIME      : %.4f\n",compT.ramp_down_time);
      fprintf(f," IDLE_DIVING         : %.4f\n",compT.idle_diving);
      fprintf(f," IDLE_NODE           : %.4f\n",compT.idle_node);
      fprintf(f," IDLE_NAMES          : %.4f\n",compT.idle_names);
      fprintf(f," IDLE_CUTS           : %.4f\n",compT.idle_cuts);
      fprintf(f," START_NODE          : %.4f\n",compT.start_node);
      fprintf(f," CUT_POOL            : %.4f\n\n",compT.cut_pool);

      fprintf(f, "########################################################\n");
      fprintf(f, " TREE DESCRIPTION \n");
      fprintf(f, "########################################################\n");

      write_tree(ws->rootnode, f);
      fclose(f);
      return TRUE;
   }   
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc *sym_read_warm_start(char *file)
{   
   FILE * f;
   char str[80], str2[80], str3[80], str4[80];
   int i=0, j=0, num=0, ch=0;
   int temp =0;
   cut_data *cut;
   problem_stat stat;
   node_times compT;
   warm_start_desc * ws;   
  
   if (!(f = fopen(file, "r"))){
      printf("sym_read_warm_start():");
      printf("Can not open the warm start file to read!\n");
      return (0);
   }
   else{
      ws = (warm_start_desc*)calloc(1,sizeof(warm_start_desc));     
      
      /* bound info */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %i", str, str, &ws->phase);
      fscanf(f,"%s %s %lf", str, str, &ws->lb);
      fscanf(f,"%s %s %i", str, str, &ch);
      ws->has_ub = (char)ch;
      fscanf(f,"%s %s %lf", str, str, &ws->ub);
      
      /* cut info */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %i", str, str, &ws->cut_num);
      fscanf(f,"%s %s %i", str, str, &temp);
      ws->allocated_cut_num = temp;
      
      if (temp){
	 ws->cuts = (cut_data **)malloc(temp *sizeof(cut_data *));
	 for(i = 0; i < ws->cut_num; i++){
	    cut = (cut_data*)malloc(sizeof(cut_data));
	    fscanf(f,"%s %i %s", str, &num, str);
	    fscanf(f,"%s %s %i", str, str, &cut->size);
	    cut->coef = (char*)malloc(CSIZE*cut->size);
	    fscanf(f,"%s %s", str, str);
	    
	    for(j=0; j<cut->size; j++){
	       fscanf(f,"%i", &ch);
	       cut->coef[j] = (char)ch;
	    } 
	    fscanf(f,"%s %s %lf", str, str, &cut->rhs);
	    fscanf(f,"%s %s %lf", str, str, &cut->range);
	    fscanf(f,"%s %s %i", str, str, &ch);
	    cut->type = (char)ch;
	    fscanf(f,"%s %s %c", str, str, &cut->sense);
	    fscanf(f,"%s %s %i", str, str, &ch);
	    cut->deletable=(char)ch;
	    fscanf(f,"%s %s %i", str, str, &ch);
	    cut->branch = (char)ch;
	    fscanf(f,"%s %s %i", str, str, &cut->name);
	    
	    ws->cuts[i] = cut;
	 }
      }
      
      /* problem stats */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %lf", str, str, &stat.root_lb);
      fscanf(f,"%s %s %i", str, str, &stat.cuts_in_pool);
      fscanf(f,"%s %s %i", str, str, &stat.max_depth);
      fscanf(f,"%s %s %i", str, str, &stat.chains);
      fscanf(f,"%s %s %i", str, str, &stat.diving_halts);
      fscanf(f,"%s %s %i", str, str, &stat.tree_size);
      fscanf(f,"%s %s %i", str, str, &stat.created);
      fscanf(f,"%s %s %i", str, str, &stat.analyzed);
      fscanf(f,"%s %s %i", str, str, &stat.leaves_before_trimming);
      fscanf(f,"%s %s %i", str, str, &stat.leaves_after_trimming);
      fscanf(f,"%s %s %i", str, str, &stat.vars_not_priced);
      fscanf(f,"%s %s %i", str, str, &ch);
      stat.nf_status = (char)ch; 
      
      ws->stat = stat;
      
      /* computation times */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %lf", str, str, &compT.communication);
      fscanf(f,"%s %s %lf", str, str, &compT.lp);
      fscanf(f,"%s %s %lf", str, str, &compT.separation);
      fscanf(f,"%s %s %lf", str, str, &compT.fixing);
      fscanf(f,"%s %s %lf", str, str, &compT.pricing);
      fscanf(f,"%s %s %lf", str, str, &compT.strong_branching);
      fscanf(f,"%s %s %lf", str, str, &compT.wall_clock_lp);
      fscanf(f,"%s %s %lf", str, str, &compT.ramp_up_tm);
      fscanf(f,"%s %s %lf", str, str, &compT.ramp_up_lp);
      fscanf(f,"%s %s %lf", str, str, &compT.ramp_down_time);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_diving);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_node);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_names);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_cuts);
      fscanf(f,"%s %s %lf", str, str, &compT.start_node);
      fscanf(f,"%s %s %lf", str, str, &compT.cut_pool);
      
      ws->comp_times = compT;
      
      /* tree description */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      ws->rootnode = (bc_node*)calloc(1,sizeof(bc_node));	 
      read_tree(ws->rootnode, f);
   }

   fclose(f);
   return ws;   
}

/*===========================================================================*/
/*===========================================================================*/

void sym_delete_warm_start(warm_start_desc *ws)
{
   int i, temp;
   if (ws) {
      if (ws->rootnode) {
	 free_subtree(ws->rootnode);
      }
      if (ws->cuts){
	 temp = ws->cut_num;
	 for(i = 0; i < temp; i++){
	    if (ws->cuts[i]){
	       if (ws->cuts[i]->coef){
		  FREE(ws->cuts[i]->coef);
	       }
	       FREE(ws->cuts[i]);
	    }
	 }
	 FREE(ws->cuts);
      }
      FREE(ws);
   }
   
   ws = 0;
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc *sym_get_warm_start(problem *p, int copy_warm_start)
{

   int i, num=0, allocated_cut_num = 0;
   warm_start_desc * ws;
   
   if (p){
      if (!p->warm_start){
	 printf("sym_get_warm_start_desc():");
	 printf("The problem warm start description is empty!\n");
	 return (0);
      }
   }
   else{
      	 printf("sym_get_warm_start_desc():");
	 printf("The problem description is empty!\n");
	 return (0);
   }

   if (copy_warm_start){
      ws = create_copy_warm_start(p->warm_start);
   }
   else{
      ws = p->warm_start;
      p->warm_start = 0;
   }  

   return ws;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_warm_start (problem *p, warm_start_desc *ws)
{

   if (!ws){
      printf("sym_set_warm_start():The warm_start desc. is empty!\n");
      return FALSE;
   }
   
   warm_start_desc * ws_copy = create_copy_warm_start(ws);
   sym_delete_warm_start(p->warm_start);
   p->warm_start = ws_copy;
   p->par.tm_par.warm_start = TRUE;
   
   return TRUE;
}

/*===========================================================================*/
/*===========================================================================*/

void sym_trim_tree(bc_node *node)
{
   free_subtree(node);
}

/*===========================================================================*/
/*===========================================================================*/

void sym_set_int_param(problem *p, char *key, int value)
{
   char *line = (char*)malloc(CSIZE*(MAX_LINE_LENGTH+1));
   sprintf(line, "%s %d", key, value);  
   set_param(p, line);
   FREE(line);
}

/*===========================================================================*/
/*===========================================================================*/

void sym_set_dbl_param(problem *p, char *key, double value)
{
   char *line = (char*)malloc(CSIZE*(MAX_LINE_LENGTH+1));
   sprintf(line, "%s %f", key, value);  
   set_param(p, line);
   FREE(line);
}

/*===========================================================================*/
/*===========================================================================*/

void sym_set_str_param(problem *p, char *key, char *value)
{
   char *line = (char*)malloc(CSIZE*(MAX_LINE_LENGTH+1));
   sprintf(line, "%s %s", key, value);  
   set_param(p, line);
   FREE(line);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_int_param(problem *p,  char *key)
{

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
   
   if (strcmp(key, "verbosity") == 0){
      return (p->par.verbosity);
   }
   else if (strcmp(key, "random_seed") == 0){
      return(p->par.random_seed);
   }
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   else if (strcmp(key, "do_decomp") == 0 ||
	    strcmp(key, "CG_do_decomp") == 0 ||
	    strcmp(key, "TM_do_decomp") == 0){
      return(tm_par->do_decomp);
   }
   /*___END_EXPERIMENTAL_SECTION___*/
   
   /***********************************************************************
    ***                    Master params                            ***
    ***********************************************************************/
   else if (strcmp(key, "M_verbosity") == 0){
      return(p->par.verbosity);
   }
   else if (strcmp(key, "M_random_seed") == 0){
      return(p->par.random_seed);
   }
   else if (strcmp(key, "tm_debug") == 0 ||
	    strcmp(key, "M_tm_debug") == 0){
      return(p->par.tm_debug);
      if (p->par.tm_debug) p->par.tm_debug = 4;
   }
   else if (strcmp(key, "dg_debug") == 0 ||
	    strcmp(key, "M_dg_debug") == 0){
      return(p->par.dg_debug);
      if (p->par.dg_debug) p->par.dg_debug = 4;
   }
   else if (strcmp(key, "pvm_trace") == 0 ||
	    strcmp(key, "M_pvm_trace") == 0){
      return(p->par.pvm_trace);
   }
   else if (strcmp(key, "do_branch_and_cut") == 0 ||
	    strcmp(key, "M_do_branch_and_cut") == 0){
      return(p->par.do_branch_and_cut);
   }
   else if (strcmp(key, "do_draw_graph") == 0 ||
	    strcmp(key, "M_do_draw_graph") == 0){
      return(p->par.do_draw_graph);
   }
   else if (strcmp(key, "use_permanent_cut_pools") == 0 ||
	    strcmp(key, "M_use_permanent_cut_pools") == 0){
      return(p->par.use_permanent_cut_pools);
   }
   
   /***********************************************************************
    ***                 DrawGraph params                            ***
    ***********************************************************************/
   else if (strcmp(key, "echo_commands") == 0 ||
	    strcmp(key, "DG_echo_commands") == 0){
      return(dg_par->echo_commands);
   }
   else if (strcmp(key, "canvas_width") == 0 ||
	    strcmp(key, "DG_canvas_width") == 0){
      return(dg_par->canvas_width);
   }
   else if (strcmp(key, "canvas_height") == 0 ||
	    strcmp(key, "DG_canvas_height") == 0){
      return(dg_par->canvas_height);
   }
   else if (strcmp(key, "viewable_width") == 0 ||
	    strcmp(key, "DG_viewable_width") == 0){
      return(dg_par->viewable_width);
   }
   else if (strcmp(key, "viewable_height") == 0 ||
	    strcmp(key, "DG_viewable_height") == 0){
      return(dg_par->viewable_width);
   }
   else if (strcmp(key, "disp_nodelabels") == 0 ||
	    strcmp(key, "DG_disp_nodelabels") == 0){
      return(dg_par->disp_nodelabels);
   }
   else if (strcmp(key, "disp_nodeweights") == 0 ||
	    strcmp(key, "DG_disp_nodeweights") == 0){
      return(dg_par->disp_nodeweights);
   }
   else if (strcmp(key, "disp_edgeweights") == 0 ||
	    strcmp(key, "DG_disp_edgeweights") == 0){
      return(dg_par->disp_edgeweights);
   }
   else if (strcmp(key, "node_radius") == 0 ||
	    strcmp(key, "DG_node_radius") == 0){
      return(dg_par->node_radius);
   }
   else if (strcmp(key, "interactive_mode") == 0 ||
	    strcmp(key, "DG_interactive_mode") == 0){
      return(dg_par->interactive_mode);
   }
   else if (strcmp(key, "mouse_tracking") == 0 ||
	    strcmp(key, "DG_mouse_tracking") == 0){
      return(dg_par->mouse_tracking);
   }

   /***********************************************************************
    ***                  Treemanager params                         ***
    ***********************************************************************/

   else if (strcmp(key, "TM_verbosity") == 0){
      return(tm_par->verbosity);
   }
   else if (strcmp(key, "lp_debug") == 0 ||
	    strcmp(key, "TM_lp_debug") == 0){
      return(tm_par->lp_debug);
   }
   else if (strcmp(key, "cg_debug") == 0 ||
	    strcmp(key, "TM_cg_debug") == 0){
      return(tm_par->cg_debug);
   }
   else if (strcmp(key, "cp_debug") == 0 ||
	    strcmp(key, "TM_cp_debug") == 0){
      return(tm_par->cp_debug);
   }
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   else if (strcmp(key, "sp_debug") == 0 ||
	    strcmp(key, "TM_sp_debug") == 0){
      return(tm_par->sp_debug);
   }
   /*___END_EXPERIMENTAL_SECTION___*/
   else if (strcmp(key, "max_active_nodes") == 0 ||
	    strcmp(key, "TM_max_active_nodes") == 0){
      return(tm_par->max_active_nodes);
   }
   else if (strcmp(key, "max_cp_num") == 0 ||
	    strcmp(key, "TM_max_cp_num") == 0){
      return(tm_par->max_cp_num);
   }
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   else if (strcmp(key, "max_sp_num") == 0 ||
	    strcmp(key, "TM_max_sp_num") == 0){
      return(tm_par->max_sp_num);
   }
   /*___END_EXPERIMENTAL_SECTION___*/
   else if (strcmp(key, "lp_mach_num") == 0 ||
	    strcmp(key, "TM_lp_mach_num") == 0){
      return(tm_par->lp_mach_num);
   }
   else if (strcmp(key, "cg_mach_num") == 0 ||
	    strcmp(key, "TM_cg_mach_num") == 0){
      return(tm_par->cg_mach_num);
   }
   else if (strcmp(key, "cp_mach_num") == 0 ||
	    strcmp(key, "TM_cp_mach_num") == 0){
      return(tm_par->cp_mach_num);
   }
#ifndef COMPILE_IN_CG
   else if (strcmp(key, "use_cg") == 0 ||
	    strcmp(key, "TM_use_cg") == 0 ||
	    strcmp(key, "LP_use_cg") == 0){
      return(tm_par->use_cg);
   }
#endif
   else if (strcmp(key, "TM_random_seed") == 0){
      return(tm_par->random_seed);
   }
   else if (strcmp(key, "diving_strategy") == 0 ||
	    strcmp(key, "TM_diving_strategy") == 0){
      return(tm_par->diving_strategy);
   }
   else if (strcmp(key, "diving_k") == 0 ||
	    strcmp(key, "TM_diving_k") == 0){
      return(tm_par->diving_k);
   }
   else if (strcmp(key, "node_selection_rule") == 0 ||
	    strcmp(key, "TM_node_selection_rule") == 0){
      return(tm_par->node_selection_rule);
   }
   else if (strcmp(key, "keep_description_of_pruned") == 0 ||
	    strcmp(key, "TM_keep_description_of_pruned") == 0){
      return(tm_par->keep_description_of_pruned);
   }
   else if (strcmp(key, "warm_start") == 0 ||
	    strcmp(key, "TM_warm_start") == 0){
      return(tm_par->warm_start);
   }
   else if (strcmp(key, "vbc_emulation") == 0 ||
	    strcmp(key, "TM_vbc_emulation") == 0){
      return(tm_par->vbc_emulation);
   }
   else if (strcmp(key, "logging_interval") == 0 ||
	    strcmp(key, "TM_logging_interval") == 0){
      return(tm_par->logging_interval);
   }
   else if (strcmp(key, "logging") == 0 ||
	    strcmp(key, "TM_logging") == 0){
      return(tm_par->logging);
   }
   else if (strcmp(key, "price_in_root") == 0 ||
	    strcmp(key, "TM_price_in_root") == 0){
      return(tm_par->price_in_root);
   }
   else if (strcmp(key, "trim_search_tree") == 0 ||
	    strcmp(key, "TM_trim_search_tree") == 0){
      return(tm_par->trim_search_tree);
   }
   else if (strcmp(key, "colgen_in_first_phase") == 0 ||
	    strcmp(key, "TM_colgen_in_first_phase") == 0){
      return(tm_par->colgen_strat[0]);
   }
   
   else if (strcmp(key, "colgen_in_second_phase") == 0 ||
	    strcmp(key, "TM_colgen_in_second_phase") == 0){
      return(tm_par->colgen_strat[1]);
   }
   else if (strcmp(key, "node_limit") == 0 ||
	    strcmp(key, "TM_node_limit") == 0){
      return(tm_par->node_limit);
   }
   else if (strcmp(key, "find_first_feasible") == 0 ||
	    strcmp(key, "TM_find_first_feasible") == 0){
      return(tm_par->find_first_feasible);
   }
   
   /***********************************************************************
    ***                      LP params                              ***
    ***********************************************************************/
   if (strcmp(key, "LP_verbosity") == 0){
      return(lp_par->verbosity);
   }
   else if (strcmp(key, "set_obj_upper_lim") == 0 ||
	    strcmp(key, "LP_set_obj_upper_lim") == 0){
      return(lp_par->set_obj_upper_lim);
   }
   
   else if (strcmp(key, "scaling") == 0 ||
	    strcmp(key, "LP_scaling") == 0){
      return(lp_par->scaling);
   }
   else if (strcmp(key, "fastmip") == 0 ||
	    strcmp(key, "LP_fastmip") == 0){
      return(lp_par->fastmip);
   }
   else if (strcmp(key, "try_to_recover_from_error") == 0 ||
	    strcmp(key, "LP_try_to_recover_from_error") == 0){
      return(lp_par->try_to_recover_from_error);
   }
   else if (strcmp(key, "problem_type") == 0 ||
	    strcmp(key, "LP_problem_type") == 0){
      return(lp_par->problem_type);
   }
   else if (strcmp(key, "not_fixed_storage_size") == 0 ||
	    strcmp(key, "LP_not_fixed_storage_size") == 0 ||
	    strcmp(key, "TM_not_fixed_storage_size") == 0 ){
      return(lp_par->not_fixed_storage_size);
   }
   else if (strcmp(key, "cut_pool_check_frequency") == 0 ||
	    strcmp(key, "LP_cut_pool_check_frequency") == 0){
      return(lp_par->cut_pool_check_freq);
   }
   else if (strcmp(key, "load_balance_level") == 0 ||
	    strcmp(key, "LP_load_balance_level") == 0){
      return(lp_par->load_balance_level);
   }
   else if (strcmp(key, "load_balance_iterations") == 0 ||
	    strcmp(key, "LP_load_balance_iterations") == 0){
      return(lp_par->load_balance_iterations);
   }
   else if (strcmp(key, "load_balance_compare_candidates") == 0 ||
	    strcmp(key, "LP_load_balance_compare_candidates") == 0){
      return(lp_par->load_balance_compare_candidates);
   }
   else if (strcmp(key, "fractional_diving_num") == 0 ||
	    strcmp(key, "LP_fractional_diving_num") == 0){
      return(lp_par->fractional_diving_num);
   }
   else if (strcmp(key, "max_cols_to_add_min") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_min") == 0){
      return(lp_par->max_non_dual_feas_to_add_min);
   }
   else if (strcmp(key, "max_non_dual_feas_to_add_max") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_max") == 0){
      return(lp_par->max_non_dual_feas_to_add_max);
   }
   else if (strcmp(key, "max_not_fixable_to_add_min") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_min") == 0){
      return(lp_par->max_not_fixable_to_add_min);
   }
   else if (strcmp(key, "max_not_fixable_to_add_max") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_max") == 0){
      return(lp_par->max_not_fixable_to_add_max);
   }
   
   else if (strcmp(key, "mat_col_compress_num") == 0 ||
	    strcmp(key, "LP_mat_col_compress_num") == 0){
      return(lp_par->mat_col_compress_num);
   }
   else if (strcmp(key, "mat_row_compress_num") == 0 ||
	    strcmp(key, "LP_mat_row_compress_num") == 0){
      return(lp_par->mat_row_compress_num);
   }
   else if (strcmp(key, "tailoff_gap_backsteps") == 0 ||
	    strcmp(key, "LP_tailoff_gap_backsteps") == 0){
      return(lp_par->tailoff_gap_backsteps);
   }
   else if (strcmp(key, "tailoff_obj_backsteps") == 0 ||
	    strcmp(key, "LP_tailoff_obj_backsteps") == 0){
      return(lp_par->tailoff_obj_backsteps);
   }
   else if (strcmp(key, "ineff_cnt_to_delete") == 0 ||
	    strcmp(key, "LP_ineff_cnt_to_delete") == 0){
      return(lp_par->ineff_cnt_to_delete);
   }
   else if (strcmp(key, "eff_cnt_before_cutpool") == 0 ||
	    strcmp(key, "LP_eff_cnt_before_cutpool") == 0){
      return(lp_par->eff_cnt_before_cutpool);
   }
   else if (strcmp(key, "ineffective_constraints") == 0 ||
	    strcmp(key, "LP_ineffective_constraints") == 0){
      return(lp_par->ineffective_constraints);
   }
   else if (strcmp(key, "base_constraints_always_effective") == 0 ||
	    strcmp(key, "LP_base_constraints_always_effective") == 0){
      return(lp_par->base_constraints_always_effective);
   }
   else if (strcmp(key, "branch_on_cuts") == 0 ||
	    strcmp(key, "LP_branch_on_cuts") == 0){
      return(lp_par->branch_on_cuts);
	     }
   else if (strcmp(key, "discard_slack_cuts") == 0 ||
	    strcmp(key, "LP_discard_slack_cuts") == 0){
      return(lp_par->discard_slack_cuts);
   }
   else if (strcmp(key, "max_cut_num_per_iter") == 0 ||
	    strcmp(key, "LP_max_cut_num_per_iter") == 0){
      return(lp_par->max_cut_num_per_iter);
   }
   
   /* variable fixing params */
   else if (strcmp(key, "do_reduced_cost_fixing") == 0 ||
	    strcmp(key, "LP_do_reduced_cost_fixing") == 0){
      return(lp_par->do_reduced_cost_fixing);
   }
   else if (strcmp(key, "do_logical_fixing") == 0 ||
	    strcmp(key, "LP_do_logical_fixing") == 0){
      return(lp_par->do_logical_fixing);
   }
   else if (strcmp(key, "fixed_to_ub_before_logical_fixing") == 0 ||
	    strcmp(key, "LP_fixed_to_ub_before_logical_fixing") == 0){
      return(lp_par->fixed_to_ub_before_logical_fixing);
   }
   else if (strcmp(key, "generate_cgl_cuts") == 0 ||
	    strcmp(key, "generate_cgl_cuts") == 0){
      return(cg_par->do_findcuts);
   }
   
   else if (strcmp(key, "max_presolve_iter") == 0 ||
	    strcmp(key, "LP_max_presolve_iter") == 0){
      return(lp_par->max_presolve_iter);
   }
   
   /* user-defined function defaults */
   else if (strcmp(key, "is_feasible_default") == 0 ||
	    strcmp(key, "LP_is_feasible_default") == 0){
      return(lp_par->is_feasible_default);
   }
   else if (strcmp(key, "send_feasible_solution_default") == 0 ||
	    strcmp(key, "LP_send_feasible_solution_default") == 0){
      return(lp_par->send_feasible_solution_default);
   }
   else if (strcmp(key, "display_solution_default") == 0 ||
	    strcmp(key, "LP_display_solution_default") == 0){
      return(lp_par->display_solution_default);
   }
   else if (strcmp(key, "shall_we_branch_default") == 0 ||
	    strcmp(key, "LP_shall_we_branch_default") == 0){
      return(lp_par->shall_we_branch_default);
   }
   else if (strcmp(key, "select_candidates_default") == 0 ||
	    strcmp(key, "LP_select_candidates_default") == 0){
      return(lp_par->select_candidates_default);
   }
   else if (strcmp(key, "strong_branching_cand_num") == 0){
      return(lp_par->strong_branching_cand_num_max);
   }
   else if (strcmp(key, "strong_branching_cand_num_max") == 0 ||
	    strcmp(key, "LP_strong_branching_cand_num_max") == 0){
      return(lp_par->strong_branching_cand_num_max);
   }
   else if (strcmp(key, "strong_branching_cand_num_min") == 0 ||
	    strcmp(key, "LP_strong_branching_cand_num_min") == 0){
      return(lp_par->strong_branching_cand_num_min);
   }
   else if (strcmp(key, "compare_candidates_default") == 0 ||
	    strcmp(key, "LP_compare_candidates_default") == 0){
      return(lp_par->compare_candidates_default);
   }

   else if (strcmp(key, "select_child_default") == 0 ||
	    strcmp(key, "LP_select_child_default") == 0){
      return(lp_par->select_child_default);
   }
   else if (strcmp(key, "pack_lp_solution_default") == 0 ||
	       strcmp(key, "LP_pack_lp_solution_default") == 0){
      return(lp_par->pack_lp_solution_default);
   }
   
   /***********************************************************************
    ***                     cut_gen params                          ***
    ***********************************************************************/
   else if (strcmp(key, "CG_verbosity") == 0){
      return(cg_par->verbosity);
   }
   else if (strcmp(key, "do_findcuts") == 0 ||
	    strcmp(key, "CG_do_findcuts") == 0){
      return(cg_par->do_findcuts);
   }
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   else if (strcmp(key, "decomp_sol_pool_check_freq") == 0 ||
	    strcmp(key, "CG_decomp_sol_pool_check_freq") == 0){
      return(cg_par->decomp_sol_pool_check_freq);
   }
   else if (strcmp(key, "decomp_wait_for_cols") == 0 ||
	    strcmp(key, "CG_decomp_wait_for_cols") == 0){
      return(cg_par->decomp_wait_for_cols);
   }
   else if (strcmp(key, "decomp_max_col_num_per_iter") == 0 ||
	    strcmp(key, "CG_decomp_max_col_num_per_iter") == 0){
     return( cg_par->decomp_max_col_num_per_iter);
   }
   else if (strcmp(key, "decomp_col_block_size") == 0 ||
	    strcmp(key, "CG_decomp_col_block_size") == 0){
      return(cg_par->decomp_col_block_size);
   }
   else if (strcmp(key, "decomp_mat_block_size") == 0 ||
	    strcmp(key, "CG_decomp_mat_block_size") == 0){
      return(cg_par->decomp_mat_block_size);
   }
   else if (strcmp(key, "decomp_complete_enum") == 0 ||
	    strcmp(key, "CG_decomp_complete_enum") == 0){
	 return(cg_par->decomp_complete_enum);
   }
   /*___END_EXPERIMENTAL_SECTION___*/
   
   /***********************************************************************
    ***                      cutpool params                         ***
    ***********************************************************************/
   else if (strcmp(key, "CP_verbosity") == 0){
      return(cp_par->verbosity);
   }
   else if (strcmp(key, "cp_warm_start") == 0 ||
	    strcmp(key, "CP_warm_start") == 0){
      return(cp_par->warm_start);
   }
   else if (strcmp(key, "cp_logging") == 0 ||
	    strcmp(key, "CP_logging") == 0){
      return(cp_par->logging);
   }
   else if (strcmp(key, "block_size") == 0 ||
	    strcmp(key, "CP_block_size") == 0){
      return(cp_par->block_size);
   }
   else if (strcmp(key, "max_size") == 0 ||
	    strcmp(key, "CP_max_size") == 0){
      return(cp_par->max_size);
   }
   else if (strcmp(key, "max_number_of_cuts") == 0 ||
	    strcmp(key, "CP_max_number_of_cuts") == 0){
      return(cp_par->max_number_of_cuts);
   }
   else if (strcmp(key, "cuts_to_check") == 0 ||
	    strcmp(key, "cuts_to_check") == 0){
      return(cp_par->cuts_to_check);
   }
   else if (strcmp(key, "delete_which") == 0 ||
	    strcmp(key, "CP_delete_which") == 0){
      return(cp_par->delete_which);
   }
   else if (strcmp(key, "touches_until_deletion") == 0 ||
	    strcmp(key, "CP_touches_until_deletion") == 0){
      return(cp_par->touches_until_deletion);
   }
   else if (strcmp(key, "min_to_delete") == 0 ||
	    strcmp(key, "CP_min_to_delete") == 0){
      return(cp_par->min_to_delete);
   }
      else if (strcmp(key, "check_which") == 0 ||
	    strcmp(key, "CP_check_which") == 0){
	 return(cp_par->check_which);
		}
/*__BEGIN_EXPERIMENTAL_SECTION__*/

   /***********************************************************************
    ***                     solpool params                          ***
    ***********************************************************************/
#ifdef COMPILE_DECOMP
   else if (strcmp(key, "SP_verbosity") == 0){
	 return(sp_par->verbosity);
   }
   else if (strcmp(key, "SP_block_size") == 0){
      return(sp_par->block_size);
   }
   else if (strcmp(key, "SP_max_size") == 0){
      return(sp_par->max_size);
   }
   else if (strcmp(key, "max_number_of_sols") == 0 ||
	    strcmp(key, "SP_max_number_of_sols") == 0){
      return(sp_par->max_number_of_sols);
   }
   else if (strcmp(key, "SP_delete_which") == 0){
      return(sp_par->delete_which);
   }
   else if (strcmp(key, "SP_touches_until_deletion") == 0){
      return(sp_par->touches_until_deletion);
   }
   else if (strcmp(key, "SP_min_to_delete") == 0){
      return(sp_par->min_to_delete);
   }
   else if (strcmp(key, "SP_compress_num") == 0){
      return(sp_par->compress_num);
   }
   else if (strcmp(key, "SP_check_which") == 0){
      return(sp_par->check_which);
   }
#endif

   /*___END_EXPERIMENTAL_SECTION___*/

   return (0);
}

/*===========================================================================*/
/*===========================================================================*/

double sym_get_dbl_param(problem *p, char *key)
{

   double timeout;

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
   
   if (strcmp(key, "granularity") == 0){
      return(tm_par->granularity);
   }
   else if (strcmp(key, "upper_bound") == 0 ||
	    strcmp(key, "M_upper_bound") == 0){
      return(p->ub);
   }
   else if (strcmp(key, "upper_bound_estimate") == 0 ||
	    strcmp(key, "M_upper_bound_estimate") == 0){
      return(p->ub_estimate);
   }
   else if (strcmp(key, "lower_bound") == 0 ||
	    strcmp(key, "M_lower_bound") == 0){
      return(p->lb);
   }
   else if (strcmp(key, "scale_factor") == 0 ||
	    strcmp(key, "DG_scale_factor") == 0){
      return(dg_par->scale_factor);
   }
   
   /***********************************************************************
    ***                  Treemanager params                         ***
    ***********************************************************************/
   else if (strcmp(key, "TM_granularity") == 0){
      return(tm_par->granularity);
   }
   else if (strcmp(key, "unconditional_dive_frac") == 0 ||
	    strcmp(key, "TM_unconditional_dive_frac") == 0){
      return(tm_par->unconditional_dive_frac);
   }
   else if (strcmp(key, "diving_threshold") == 0 ||
	    strcmp(key, "TM_diving_threshold") == 0){
     return( tm_par->diving_threshold);
   }
   else if (strcmp(key, "time_limit") == 0 ||
	    strcmp(key, "TM_time_limit") == 0){
     return( tm_par->time_limit);
   }
   else if (strcmp(key, "gap_limit") == 0 ||
	    strcmp(key, "TM_gap_limit") == 0){
      return(tm_par->gap_limit);
   }
   
   /***********************************************************************
    ***                      LP params                              ***
    ***********************************************************************/
   else if (strcmp(key, "LP_granularity") == 0){
      return(lp_par->granularity);
   }
   else if (strcmp(key, "fractional_diving_ratio") == 0 ||
	    strcmp(key, "LP_fractional_diving_ratio") == 0){
      return(lp_par->fractional_diving_ratio);
   }
   else if (strcmp(key, "max_non_dual_feas_to_add_frac") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_frac") == 0){
      return(lp_par->max_non_dual_feas_to_add_frac);
   }
   else if (strcmp(key, "max_not_fixable_to_add_frac") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_frac") == 0){
      return(lp_par->max_not_fixable_to_add_frac);
   }
   else if (strcmp(key, "mat_col_compress_ratio") == 0 ||
	    strcmp(key, "LP_mat_col_compress_ratio") == 0){
      return(lp_par->mat_col_compress_ratio);
   }
   else if (strcmp(key, "mat_row_compress_ratio") == 0 ||
	    strcmp(key, "LP_mat_row_compress_ratio") == 0){
      return(lp_par->mat_row_compress_ratio);
   }
   else if (strcmp(key, "tailoff_gap_frac") == 0 ||
	    strcmp(key, "LP_tailoff_gap_frac") == 0){
      return(lp_par->tailoff_gap_frac);
   }
   else if (strcmp(key, "tailoff_obj_frac") == 0 ||
	    strcmp(key, "LP_tailoff_obj_frac") == 0){
      return(lp_par->tailoff_obj_frac);
   }
   else if (strcmp(key, "tailoff_absolute") == 0 ||
	    strcmp(key, "LP_tailoff_absolute") == 0){
      return(lp_par->tailoff_absolute);
   }

   /* timeouts on receiving cuts */
   else if (strcmp(key, "first_lp_first_cut_time_out") == 0 ||
	    strcmp(key, "LP_first_lp_first_cut_time_out") == 0){
      return(lp_par->first_lp.first_cut_time_out);
   }
   else if (strcmp(key, "first_lp_all_cuts_time_out") == 0 ||
	    strcmp(key, "LP_first_lp_all_cuts_time_out") == 0){
      return(lp_par->first_lp.all_cuts_time_out);
   }
   else if (strcmp(key, "later_lp_first_cut_time_out") == 0 ||
	    strcmp(key, "LP_later_lp_first_cut_time_out") == 0){
      return(lp_par->later_lp.first_cut_time_out);
   }
   else if (strcmp(key, "later_lp_all_cuts_time_out") == 0 ||
	    strcmp(key, "LP_later_lp_all_cuts_time_out") == 0){
      return(lp_par->later_lp.all_cuts_time_out);
   }

   else if (strcmp(key, "gap_as_ub_frac") == 0 ||
	    strcmp(key, "LP_gap_as_ub_frac") == 0){
      return(lp_par->gap_as_ub_frac);
   }
   else if (strcmp(key, "gap_as_last_gap_frac") == 0 ||
	    strcmp(key, "LP_gap_as_last_gap_frac") == 0){
      return(lp_par->gap_as_last_gap_frac);
   }
   else if (strcmp(key, "fixed_to_ub_frac_before_logical_fixing")==0 ||
	    strcmp(key, "LP_fixed_to_ub_frac_before_logical_fixing")==0){
      return(lp_par->fixed_to_ub_frac_before_logical_fixing);
   }
   else if (strcmp(key,"strong_branching_red_ratio") == 0 ||
	    strcmp(key,"LP_strong_branching_red_ratio") == 0){
      return(lp_par->strong_branching_red_ratio);
   }
   
   /***********************************************************************
    ***                     cut_gen params                          ***
    ***********************************************************************/
   else if (strcmp(key, "decomp_initial_timeout") == 0 ||
	    strcmp(key, "CG_decomp_initial_timeout") == 0){
      return(cg_par->decomp_initial_timeout);
   }
   else if (strcmp(key, "decomp_dynamic_timeout") == 0 ||
	    strcmp(key, "CG_decomp_dynamic_timeout") == 0){
      return(cg_par->decomp_dynamic_timeout);
   }

   /***********************************************************************
    ***                     solpool params                          ***
    ***********************************************************************/
#ifdef COMPILE_DECOMP
   else if (strcmp(key, "SP_etol") == 0){
      return(sp_par->etol);
   }
   else if (strcmp(key, "SP_compress_ratio") == 0){
      return(sp_par->compress_ratio);
   }
#endif

   return (0);
}

/*===========================================================================*/
/*===========================================================================*/

char *sym_get_str_param(problem *p, char *key)
{

   int len, i;
   
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
   
   if (strcmp(key, "problem_name") == 0){      
      return (p->probname);
   }  
   else if (strcmp(key, "tm_executable_name") == 0 ||
	    strcmp(key, "tm_exe") == 0 ||
	    strcmp(key, "M_tm_exe") == 0 ||
	    strcmp(key, "M_tm_executable_name") == 0){
      return(p->par.tm_exe);
   }
   else if (strcmp(key, "dg_executable_name") == 0 ||
	    strcmp(key, "dg_exe") == 0 ||
	    strcmp(key, "M_dg_exe") == 0 ||
	    strcmp(key, "M_dg_executable_name") == 0){
      return(p->par.dg_exe); 
   }
   else if (strcmp(key, "tm_machine") == 0 ||
	    strcmp(key, "M_tm_machine") == 0){
      return(p->par.tm_machine);
   }
   else if (strcmp(key, "dg_machine") == 0 ||
	    strcmp(key, "M_dg_machine") == 0){
      return(p->par.dg_machine);
   }

   /***********************************************************************
    ***                 DrawGraph params                            ***
    ***********************************************************************/
   
   else if (strcmp(key, "source_path") == 0 ||
	    strcmp(key, "DG_source_path") == 0){
      return(dg_par->source_path);
   }
   else if (strcmp(key, "node_dash") == 0 ||
	    strcmp(key, "DG_node_dash") == 0){
      return(dg_par->node_dash);
   }
   else if (strcmp(key, "edge_dash") == 0 ||
	    strcmp(key, "DG_edge_dash") == 0){
      return(dg_par->edge_dash);
   }
   else if (strcmp(key, "nodelabel_font") == 0 ||
	    strcmp(key, "DG_nodelabel_font") == 0){
      return(dg_par->nodelabel_font);
   }
   else if (strcmp(key, "nodeweight_font") == 0 ||
	    strcmp(key, "DG_nodeweight_font") == 0){
      return(dg_par->nodeweight_font);
   }
   else if (strcmp(key, "edgeweight_font") == 0 ||
	    strcmp(key, "DG_edgeweight_font") == 0){
      return(dg_par->edgeweight_font);
   }
   
   /***********************************************************************
    ***                  Treemanager params                         ***
    ***********************************************************************/
   else if (strcmp(key, "lp_executable_name") == 0 ||
	    strcmp(key, "lp_exe") == 0 ||
	    strcmp(key, "TM_lp_exe") == 0 ||
	    strcmp(key, "TM_lp_executable_name") == 0){
      return(tm_par->lp_exe);
   }
   else if (strcmp(key, "cg_executable_name") == 0 ||
	    strcmp(key, "cg_exe") == 0 ||
	    strcmp(key, "TM_cg_exe") == 0 ||
	    strcmp(key, "TM_cg_executable_name") == 0){
      return(tm_par->cg_exe);
   }
   else if (strcmp(key, "cp_executable_name") == 0 ||
	    strcmp(key, "cp_exe") == 0 ||
	    strcmp(key, "TM_cp_exe") == 0 ||
	    strcmp(key, "TM_cp_executable_name") == 0){
      return(tm_par->cp_exe);
   }
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   else if (strcmp(key, "sp_executable_name") == 0 ||
	    strcmp(key, "sp_exe") == 0 ||
	    strcmp(key, "TM_sp_exe") == 0 ||
	    strcmp(key, "TM_sp_executable_name") == 0){
      return(tm_par->sp_exe);
   }

   return (0);
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc *sym_create_copy_warm_start(warm_start_desc *ws)
{
   return create_copy_warm_start(ws);
}

/*===========================================================================*/
/*===========================================================================*/

MIPdesc *sym_create_copy_mip_desc(problem * p)
{
   if (p){
      return create_copy_mip_desc(p->mip);
   }
   else{
      printf("sym_create_copy_mip_desc():");
      printf("An empty problem is given!\n");
      return 0;
   }
}

/*===========================================================================*/
/*===========================================================================*/

problem * sym_create_copy_problem (problem *p)
{
   return create_copy_problem(p);
}

/*===========================================================================*/
/*===========================================================================*/
double sym_get_lb_for_new_rhs(problem *p, int cnt, int *new_rhs_ind, 
			      double *new_rhs_val)
{
#ifdef USE_CGL_CUTS
   printf("sym_get_lb_for_new_rhs():\n");
   printf("SYMPHONY can not analyse the warm start - rhs change when cuts exist, for now!\n"); 
   return(0);
#else
   if (!p || !p->mip || 
      p->par.tm_par.keep_description_of_pruned != KEEP_IN_MEMORY){ 
      printf("sym_get_lb_for_new_rhs():\n");
      printf("Trying to read an empty problem, an empty problem description"); 
      printf(" or tree nodes were not kept in memory!\n");
      return 0;
   }
   else{
      if (!p->warm_start){
	 printf("sym_get_lb_for_new_rhs():\n");
	 printf("No available tree to incur sens. analysis on. \n");
	 return 0;
      }
      else{
	 return get_lb_for_new_rhs(p->warm_start->rootnode, p->mip, cnt, 
				   new_rhs_ind, new_rhs_val);
      }
   }
#endif
}
/*===========================================================================*/
/*===========================================================================*/
