/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2008 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILE_FOR_LP

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sym_proccomm.h"
#include "sym_qsort.h"
#include "sym_lp.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_types.h"
#include "sym_pack_cut.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains general LP functions.
\*===========================================================================*/

/*===========================================================================*\
 * This function receives the problem data (if we are running in parallel)   
 * and intitializes the data structures.                                     
\*===========================================================================*/

int lp_initialize(lp_prob *p, int master_tid)
{
#ifndef COMPILE_IN_LP
   int msgtag, bytes, r_bufid;
#endif
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP)
   int s_bufid;
#endif
   int i, j;
   row_data *rows;
   var_desc **vars;

#ifdef COMPILE_IN_LP

   p->master = master_tid;

#else
   
   /* set stdout to be line buffered */
   setvbuf(stdout, (char *)NULL, _IOLBF, 0);

   register_process();

   /*------------------------------------------------------------------------*\
    * Receive tid info; request and receive problem specific data
   \*------------------------------------------------------------------------*/
   r_bufid = receive_msg(ANYONE, MASTER_TID_INFO);
   bufinfo(r_bufid, &bytes, &msgtag, &p->tree_manager);
   receive_int_array(&p->master, 1);
   receive_int_array(&p->proc_index, 1);
   freebuf(r_bufid);

#endif

   p->lp_data = (LPdata *) calloc(1, sizeof(LPdata));
   p->lp_data->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
   
#pragma omp critical (lp_solver)
   open_lp_solver(p->lp_data);

   (void) used_time(&p->tt);

#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP)
   s_bufid = init_send(DataInPlace);
   send_msg(p->master, REQUEST_FOR_LP_DATA);
   freebuf(s_bufid);
   int termcode; 
   CALL_WRAPPER_FUNCTION( receive_lp_data_u(p) );
#endif
   
   if (p->par.tailoff_gap_backsteps > 0 ||
       p->par.tailoff_obj_backsteps > 1){
      i = MAX(p->par.tailoff_gap_backsteps, p->par.tailoff_obj_backsteps);
      p->obj_history = (double *) malloc((i + 1) * DSIZE);
      for (j = 0; j <= i; j++){
	 p->obj_history[j] = -DBL_MAX;
      }
   }
#ifndef COMPILE_IN_LP
   if (p->par.use_cg){
      r_bufid = receive_msg(p->tree_manager, LP__CG_TID_INFO);
      receive_int_array(&p->cut_gen, 1);
      freebuf(r_bufid);
   }
#endif
   p->lp_data->rows =
      (row_data *) malloc((p->base.cutnum + BB_BUNCH) * sizeof(row_data));
   rows = p->lp_data->rows;
   for (i = p->base.cutnum - 1; i >= 0; i--){
      ( rows[i].cut = (cut_data *) malloc(sizeof(cut_data)) )->coef = NULL;
   }

   if (p->base.varnum > 0){
      vars = p->lp_data->vars = (var_desc **)
	 malloc(p->base.varnum * sizeof(var_desc *));
      for (i = p->base.varnum - 1; i >= 0; i--){
	 vars[i] = (var_desc *) malloc( sizeof(var_desc) );
	 vars[i]->userind = p->base.userind[i];
	 vars[i]->colind = i;
      }
   }

   /* Just to make sure this array is sufficently big */
   p->lp_data->not_fixed = (int *) malloc(p->par.not_fixed_storage_size*ISIZE);
   p->lp_data->tmp.iv = (int *) malloc(p->par.not_fixed_storage_size* 2*ISIZE);
   p->lp_data->tmp.iv_size = 2*p->par.not_fixed_storage_size;
   p->lp_data->cgl = p->par.cgl;

#ifdef COMPILE_IN_CG
   if (!p->cgp){
      p->cgp = (cg_prob *) calloc(1, sizeof(cg_prob));
   }
   
   cg_initialize(p->cgp, p->master);
#endif

   return(FUNCTION_TERMINATED_NORMALLY);
}   

/*===========================================================================*/

/*===========================================================================*\
 * This function continues to dive down the current chain until told to stop
 * by the tree manager. 
\*===========================================================================*/

int process_chain(lp_prob *p)
{
   int termcode;
   
   /* Create the LP */
   if ((termcode = create_subproblem_u(p)) < 0){
      /* User had problems creating initial LP. Abandon node. */
      return(termcode);
   }

   p->last_gap = 0.0;
   p->dive = CHECK_BEFORE_DIVE;
   if (p->has_ub && p->par.set_obj_upper_lim)
      set_obj_upper_lim(p->lp_data, p->ub - p->par.granularity);
   
   if (p->colgen_strategy & COLGEN_REPRICING){
      if (p->par.verbosity > 0){
	 printf("****************************************************\n");
	 printf("* Now repricing NODE %i LEVEL %i\n",
		p->bc_index, p->bc_level);
	 printf("****************************************************\n\n");
      }
      termcode = repricing(p);
      free_node_dependent(p);
   }else{
      if (p->par.verbosity > 0){
	 printf("****************************************************\n");
	 printf("* Now processing NODE %i LEVEL %i (from TM)\n",
		p->bc_index, p->bc_level);
	 printf("****************************************************\n\n");
	 PRINT(p->par.verbosity, 4, ("Diving set to %i\n\n", p->dive));
      }
      termcode = fathom_branch(p);

#ifdef COMPILE_IN_LP
      p->tm->stat.chains++;
      p->tm->active_node_num--;
      free_node_dependent(p);
#else
      /* send_lp_is_free()  calls  free_node_dependent() */
      send_lp_is_free(p);
#endif
   }
   p->lp_data->col_set_changed = TRUE;

   p->comp_times.lp += used_time(&p->tt);

   return(termcode);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function receives information for an active node, processes that     *
 * node and then decides which one of the children of that node should be    *
 * processed next. It then recursively processes the child until the branch  *
 * is pruned at point                                                        *
\*===========================================================================*/

int fathom_branch(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   node_times *comp_times = &p->comp_times;
   char first_in_loop = TRUE;
   int iterd, termcode;
   int cuts = 0, no_more_cuts_count;
   int num_errors = 0;
   int cut_term = 0;
   double obj_before_cuts = 0;
#ifdef DO_TESTS
   double oldobjval = lp_data->objval;
#endif

   check_ub(p);
   p->iter_num = p->node_iter_num = 0;

   // TODO: replace check_bounds with a better preprocessor
   termcode = LP_OPTIMAL; // just to initialize
   check_bounds(p, &termcode);
   if (termcode == LP_D_UNBOUNDED) {
      if (fathom(p, FALSE)) {
         comp_times->communication += used_time(&p->tt);
         return(FUNCTION_TERMINATED_NORMALLY);
      }
   }

  
   /*------------------------------------------------------------------------*\
    * The main loop -- continue solving relaxations until no new cuts
    * are found
   \*------------------------------------------------------------------------*/

   while (TRUE){
      if (p->par.branch_on_cuts && p->slack_cut_num > 0){
	 switch (p->par.discard_slack_cuts){
	  case DISCARD_SLACKS_WHEN_STARTING_NEW_NODE:
	    if (p->iter_num != 0)
	       break;
	  case DISCARD_SLACKS_BEFORE_NEW_ITERATION:
	    free_cuts(p->slack_cuts, p->slack_cut_num);
	    p->slack_cut_num = 0;
	    break;
	 }
      }

      p->iter_num++;
      p->node_iter_num++;
      lp_data->lp_count++;

      PRINT(p->par.verbosity, 2,
	    ("\n\n**** Starting iteration %i ****\n\n", p->iter_num));

      p->bound_changes_in_iter = 0;
      if (p->iter_num < 2 && (p->par.should_warmstart_chain == FALSE || 
               p->bc_level < 1)) {
         if (p->bc_level < 1) {
            PRINT(p->par.verbosity, -1, ("solving root lp relaxation\n"));
         }
         termcode = initial_lp_solve(lp_data, &iterd);
         if (p->bc_level < 1) {
            PRINT(p->par.verbosity, -1, ("finished solving lp relaxation.\n"));
         }
      } else {
         termcode = dual_simplex(lp_data, &iterd);
      }
      p->lp_stat.lp_calls++;

#ifdef DO_TESTS
      if (lp_data->objval < oldobjval - .01){
	 printf ("#####Error: LP objective value decrease from %.3f to %.3f\n",
		 oldobjval, lp_data->objval);
      }
#endif
      
      /* Get relevant data */
      get_dj_pi(lp_data);
      get_slacks(lp_data);
      get_x(lp_data);
      
      /* display the current solution */
      if (p->mip->obj_sense == SYM_MAXIMIZE){
	 PRINT(p->par.verbosity, 2, ("The LP value is: %.3f [%i,%i]\n\n",
				     -lp_data->objval + p->mip->obj_offset,
				     termcode, iterd));

      }else{
	 PRINT(p->par.verbosity, 2, ("The LP value is: %.3f [%i,%i]\n\n",
				     lp_data->objval+ p->mip->obj_offset,
				     termcode, iterd));
      }
      switch (termcode){
       case LP_D_ITLIM:      /* impossible, since itlim is set to infinity */
       case LP_D_INFEASIBLE: /* this is impossible (?) as of now */
       case LP_ABANDONED:
	 printf("####### Unexpected termcode: %i \n", termcode);
	 if (p->par.try_to_recover_from_error && (++num_errors == 1)){
	    /* Try to resolve it from scratch */
	    printf("####### Trying to recover by resolving from scratch...\n");
	    
	    continue;
	 }else{
	    char name[50] = "";
	    printf("####### Recovery failed. %s%s",
		   "LP solver is having numerical difficulties :(.\n",
		   "####### Dumping current LP to MPS file and exiting.\n\n");
	    sprintf(name, "matrix.%i.%i", p->bc_index, p->iter_num);
	    write_mps(lp_data, name);
	    return(ERROR__NUMERICAL_INSTABILITY);
	 }

       case LP_D_UNBOUNDED: /* the primal problem is infeasible */
       case LP_D_OBJLIM:
       case LP_OPTIMAL:
	 if (num_errors == 1){
	    printf("####### Recovery succeeded! Continuing with node...\n\n");
	    num_errors = 0;
	 }
	 if (termcode == LP_D_UNBOUNDED){
	    PRINT(p->par.verbosity, 1, ("Feasibility lost -- "));
#if 0
	    char name[50] = "";
	    sprintf(name, "matrix.%i.%i", p->bc_index, p->iter_num);
	    write_mps(lp_data, name);
#endif
	 }else if ((p->has_ub && lp_data->objval > p->ub - p->par.granularity)
		   || termcode == LP_D_OBJLIM){
	    PRINT(p->par.verbosity, 1, ("Terminating due to high cost -- "));
	 }else{ /* optimal and not too high cost */
#ifdef COMPILE_IN_LP
            if (p->node_iter_num < 2 && p->bc_index > 0 && 
                  p->par.should_use_rel_br) {
               update_pcost(p);
            }
            if (cuts > 0) {
               p->lp_stat.cuts_added_to_lps += cuts;
            }
            if (p->node_iter_num > 0 && p->bc_level > 0) {
               if (cuts > 0) {
                  p->lp_stat.num_cuts_added_in_path += cuts;
               }
               if (p->lp_stat.avg_cuts_obj_impr_in_path > 0) {
                  p->lp_stat.avg_cuts_obj_impr_in_path = 
                     (p->lp_stat.avg_cuts_obj_impr_in_path *
                      (p->lp_stat.num_cut_iters_in_path-1) + p->lp_data->objval - 
                      obj_before_cuts)/p->lp_stat.num_cut_iters_in_path;
               }
            } 
            obj_before_cuts = lp_data->objval;
            comp_times->lp += used_time(&p->tt);
#endif
            break;
	 }
	 comp_times->lp += used_time(&p->tt);
	 if (fathom(p, (termcode != LP_D_UNBOUNDED))){
	    comp_times->communication += used_time(&p->tt);
	    return(FUNCTION_TERMINATED_NORMALLY);
	 }else{
	    first_in_loop = FALSE;
	    comp_times->communication += used_time(&p->tt);
	    continue;
	 }
      }

      /* If come to here, the termcode must have been OPTIMAL and the
       * cost cannot be too high. */
      /* is_feasible_u() fills up lp_data->x, too!! */
      if (is_feasible_u(p, FALSE, FALSE) == IP_FEASIBLE){
	 cuts = -1;
      }else{
	 /*------------------------------------------------------------------*\
	  * send the current solution to the cut generator, and also to the
	  * cut pool if we are either
	  *  - at the beginning of a chain (but not in the root in the
	  *         first phase)
	  *  - or this is the cut_pool_check_freq-th iteration.
	 \*------------------------------------------------------------------*/
	 cuts = 0;
	 no_more_cuts_count = 0;
	 if (p->cut_pool &&
	     ((first_in_loop && (p->bc_level>0 || p->phase==1)) ||
	      (p->iter_num % p->par.cut_pool_check_freq == 0)) ){
	    no_more_cuts_count += send_lp_solution_u(p, p->cut_pool);
	 }
	 if (p->cut_gen){
	    no_more_cuts_count += send_lp_solution_u(p, p->cut_gen);
	 }

	 if (p->par.verbosity > 4){
	    printf ("Now displaying the relaxed solution ...\n");
	    display_lp_solution_u(p, DISP_RELAXED_SOLUTION);
	 }

	 comp_times->lp += used_time(&p->tt);

	 tighten_bounds(p);

	 comp_times->fixing += used_time(&p->tt);

	 if (!first_in_loop){
	    cuts = check_row_effectiveness(p);
	 }

	 /*------------------------------------------------------------------*\
	  * receive the cuts from the cut generator and the cut pool
	 \*------------------------------------------------------------------*/

#ifdef USE_SYM_APPLICATION
            if ((cut_term = receive_cuts(p, first_in_loop,
                        no_more_cuts_count)) >=0 ){
               cuts += cut_term;
            }else{
               return(ERROR__USER);
            }
#else
         if (!check_tailoff(p)) {
            if ((cut_term = receive_cuts(p, first_in_loop,
                        no_more_cuts_count)) >=0 ){
               cuts += cut_term;
            }else{
               return(ERROR__USER);
            }
         }
#endif
      }

      comp_times->lp += used_time(&p->tt);
      if (cuts < 0){ /* i.e. feasible solution is found */
	 if (fathom(p, TRUE)){
	    return(FUNCTION_TERMINATED_NORMALLY);
	 }else{
	    first_in_loop = FALSE;
	    check_ub(p);
	    continue;
	 }
      }

      PRINT(p->par.verbosity, 2,
	    ("\nIn iteration %i, before calling branch()\n", p->iter_num));
      if (cuts == 0){
	 PRINT(p->par.verbosity, 2, ("... no cuts were added.\n"));
	 if (p->par.verbosity > 4){
	    printf("Now displaying final relaxed solution...\n\n");
	    display_lp_solution_u(p, DISP_FINAL_RELAXED_SOLUTION);
	 }
      }else{
	 PRINT(p->par.verbosity, 2,
	       ("... %i violated cuts were added\n", cuts));
      }
      
      comp_times->lp += used_time(&p->tt);

      switch (cuts = branch(p, cuts)){

       case NEW_NODE:
#ifndef ROOT_NODE_ONLY
	 if (p->par.verbosity > 0){
	    printf("*************************************************\n");
	    printf("* Now processing NODE %i LEVEL %i\n",
		   p->bc_index, p->bc_level);
	    printf("*************************************************\n\n");
	 }
	 p->node_iter_num = 0;
         /*
         printf("node = %d\n", p->bc_index);
         printf("cut iters = %d\n", p->lp_stat.num_cut_iters_in_path);
         printf("cuts added = %d\n", p->lp_stat.num_cuts_added_in_path);
         printf("cut removed = %d\n", p->lp_stat.num_cuts_slacked_out_in_path);
         printf("cut obj impr = %f\n", p->lp_stat.avg_cuts_obj_impr_in_path);

         printf("strong br cands = %d\n", p->lp_stat.num_str_br_cands_in_path);
         printf("str br impr = %f\n", p->lp_stat.avg_br_obj_impr_in_path);

         printf("fp calls = %d\n", p->lp_stat.num_fp_calls_in_path);
         */
	 break;
#endif
       case FATHOMED_NODE:
	 comp_times->strong_branching += used_time(&p->tt);
	 return(FUNCTION_TERMINATED_NORMALLY);

       case BRANCHING_INF_NODE:
	 comp_times->strong_branching += used_time(&p->tt);
	 if (fathom(p, FALSE)){
	    return(FUNCTION_TERMINATED_NORMALLY);
	 }else{
	    return(FUNCTION_TERMINATED_ABNORMALLY);
	 }

       case ERROR__NO_BRANCHING_CANDIDATE: /* Something went wrong */
	 return(ERROR__NO_BRANCHING_CANDIDATE);

       case FEAS_SOL_FOUND:
         PRINT(p->par.verbosity,2,("solution found before branching\n"));
       default: /* the return value is the number of cuts added */
	 if (p->par.verbosity > 2){
	    printf("Continue with this node.");
	    if (cuts > 0)
	       printf(" %i cuts added altogether in iteration %i\n",
		      cuts, p->iter_num);
            if (p->bound_changes_in_iter > 0) {
               printf(" %i bounds added altogether in iteration %i\n",
                     p->bound_changes_in_iter, p->iter_num);
            }
	    printf("\n\n");
	 }
	 break;
      }
      comp_times->strong_branching += used_time(&p->tt);

      check_ub(p);
      first_in_loop = FALSE;

#ifdef COMPILE_IN_LP
      if (p->tm->par.time_limit >= 0.0 &&
	  wall_clock(NULL) - p->tm->start_time >= p->tm->par.time_limit){
#else
      if (p->par.time_limit >= 0.0 &&
	  wall_clock(NULL) - p->start_time >= p->par.time_limit){
#if 0
      to unconfuse vi
      }
#endif
#endif
	 if (fathom(p, TRUE)){
	    return(FUNCTION_TERMINATED_NORMALLY);
	 }else{
	    return(FUNCTION_TERMINATED_ABNORMALLY);
	 }
      }
   }

   comp_times->lp += used_time(&p->tt);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

/* fathom() returns true if it has really fathomed the node, false otherwise
   (i.e., if it had added few variables) */

int fathom(lp_prob *p, int primal_feasible)
{
   LPdata *lp_data = p->lp_data;
   our_col_set *new_cols = NULL;
   int new_vars;
   int colgen = p->colgen_strategy & COLGEN__FATHOM;
   int termcode = p->lp_data->termcode;
   
   if (p->lp_data->nf_status == NF_CHECK_NOTHING){
      PRINT(p->par.verbosity, 1,
	    ("fathoming node (no more cols to check)\n\n"));
      if (primal_feasible){
	 switch (termcode){
	  case LP_OPT_FEASIBLE:
	    send_node_desc(p, FEASIBLE_PRUNED);
	    break;
	  case LP_OPTIMAL:
	    send_node_desc(p, INTERRUPTED_NODE);
	    break;
	  default:
	    send_node_desc(p, OVER_UB_PRUNED);
	    break;
	 }
      }else{
	 send_node_desc(p, INFEASIBLE_PRUNED);
      }
      return(TRUE);
   }
	 
   if (p->colgen_strategy & COLGEN_REPRICING)
      colgen = FATHOM__GENERATE_COLS__RESOLVE;

   switch (colgen){
    case FATHOM__DO_NOT_GENERATE_COLS__DISCARD:
      PRINT(p->par.verbosity, 1, ("Pruning node\n\n"));
      send_node_desc(p, termcode == LP_OPT_FEASIBLE ? FEASIBLE_PRUNED :
		     DISCARDED_NODE);
      return(TRUE);

    case FATHOM__DO_NOT_GENERATE_COLS__SEND:
      PRINT(p->par.verbosity, 1, ("Sending node for pricing\n\n"));
      send_node_desc(p, primal_feasible ? OVER_UB_HOLD_FOR_NEXT_PHASE :
		     INFEASIBLE_HOLD_FOR_NEXT_PHASE);
      return(TRUE);

    case FATHOM__GENERATE_COLS__RESOLVE:
      check_ub(p);
      /* Note that in case of COLGEN_REPRICING we must have UB. */
      if (! p->has_ub){
	 PRINT(p->par.verbosity, 1,
	       ("\nCan't generate cols before sending (no UB)\n"));
	 send_node_desc(p, primal_feasible ? OVER_UB_HOLD_FOR_NEXT_PHASE :
			INFEASIBLE_HOLD_FOR_NEXT_PHASE);
	 return(TRUE);
      }
      PRINT(p->par.verbosity, 1,
	    ("\nGenerating columns before fathoming/resolving\n"));
      new_cols = price_all_vars(p);
      p->comp_times.pricing += used_time(&p->tt);
      new_vars = new_cols->num_vars + new_cols->rel_ub + new_cols->rel_lb;
      if (new_cols->dual_feas == NOT_TDF){
	 /* Don't have total dual feasibility. The non-dual-feasible vars
	  * have already been added. Go back and resolve. */
	 PRINT(p->par.verbosity, 2,
	       ("%i variables added in price-out.\n", new_vars));
	 free_col_set(&new_cols);
	 return(FALSE);
      }
      /* Now we know that we have total dual feasibility */
      if ((p->has_ub && lp_data->objval > p->ub - p->par.granularity) ||
	  termcode == LP_D_OBJLIM || termcode == LP_OPT_FEASIBLE){
	 /* fathomable */
	 if (termcode == LP_D_OBJLIM ||
	     (p->has_ub && lp_data->objval > p->ub - p->par.granularity)){
	    PRINT(p->par.verbosity, 1,
		  ("Fathoming node (discovered tdf & high cost)\n\n"));
	 }else{
	    PRINT(p->par.verbosity, 1,
		  ("Fathoming node (discovered tdf & feasible)\n\n"));
	 }
	 send_node_desc(p, termcode == LP_OPT_FEASIBLE ? FEASIBLE_PRUNED :
			OVER_UB_PRUNED);
	 free_col_set(&new_cols);
	 return(TRUE);
      }
      /* If we ever arrive here then we must have tdf and the function
       * was called with a primal infeasible LP.
       *
       * Again, note that in case of COLGEN_REPRICING, since we do that
       * only in the root node, the lp relaxation MUST be primal feasible,
       *
       * If TDF_HAS_ALL, then whatever can be used to restore
       * primal feasibility is already in the matrix so don't bother
       * to figure out restorability, just return and resolve the problem
       * (if new_vars == 0 then even returning is unnecessary, the node
       * can be fathomed, nothing can restore feasibility).
       */
      if (new_cols->dual_feas == TDF_HAS_ALL){
	 if (new_vars == 0){
	    PRINT(p->par.verbosity, 1,
		  ("fathoming node (no more cols to check)\n\n"));
	    send_node_desc(p, INFEASIBLE_PRUNED);
	    free_col_set(&new_cols);
	    return(TRUE);
	 }else{
	    free_col_set(&new_cols);
	    return(FALSE);
	 }
      }
      /* Sigh. There were too many variables not fixable even though we have
       * proved tdf. new_cols contains a good many of the non-fixables, use
       * new_cols to start with in restore_lp_feasibility(). */
      if (! restore_lp_feasibility(p, new_cols)){
	 PRINT(p->par.verbosity, 1,
	       ("Fathoming node (discovered tdf & not restorable inf.)\n\n"));
	 send_node_desc(p, INFEASIBLE_PRUNED);
	 free_col_set(&new_cols);
	 return(TRUE);
      }
      /* So primal feasibility is restorable. Exactly one column has been
       * added (released or a new variable) to destroy the proof of
       * infeasibility */
      free_col_set(&new_cols);
      p->comp_times.pricing += used_time(&p->tt);
      return(FALSE);
   }

   return(TRUE); /* fake return */
}

/*****************************************************************************/
/*****************************************************************************/
/* NOTE: this version of repricing works ONLY for repricing in the root node */
/*****************************************************************************/
/*****************************************************************************/

int repricing(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   node_times *comp_times = &p->comp_times;
   int iterd, termcode;
   int num_errors = 0;
   our_col_set *new_cols = NULL;
   int dual_feas, new_vars, cuts, no_more_cuts_count;
   int cut_term = 0;
   
   check_ub(p);
   p->iter_num = 0;
  
   /*------------------------------------------------------------------------*\
    * The main loop -- continue solving relaxations until TDF
   \*------------------------------------------------------------------------*/

   while (TRUE){
      p->iter_num++;

      PRINT(p->par.verbosity, 2,
	    ("\n\n**** Starting iteration %i ****\n\n", p->iter_num));

      termcode = dual_simplex(lp_data, &iterd);
      p->lp_stat.lp_calls++;

      /* Get relevant data */
      get_dj_pi(lp_data);
      get_slacks(lp_data);

      /* display the current solution */
      if (p->mip->obj_sense == SYM_MAXIMIZE){
	 PRINT(p->par.verbosity, 2, ("The LP value is: %.3f [%i,%i]\n\n",
				     -lp_data->objval + p->mip->obj_offset,
				     termcode, iterd));

      }else{
	 PRINT(p->par.verbosity, 2, ("The LP value is: %.3f [%i,%i]\n\n",
				     lp_data->objval+ p->mip->obj_offset,
				     termcode, iterd));
      }
      comp_times->lp += used_time(&p->tt);

      switch (termcode){
       case LP_D_ITLIM:      /* impossible, since itlim is set to infinity */
       case LP_D_INFEASIBLE: /* this is impossible (?) as of now */
       case LP_ABANDONED:
	 printf("######## Unexpected termcode: %i \n", termcode);
	 if (p->par.try_to_recover_from_error && (++num_errors == 1)){
	    /* Try to resolve it from scratch */
	    printf("######## Trying to recover by resolving from scratch...\n");
	    
	    continue;
	 }else{
	    char name[50] = "";
	    printf("######## Recovery failed. %s%s",
		   "LP solver is having numerical difficulties :(.\n",
		   "######## Dumping current LP to MPS file and exiting.\n\n");
	    sprintf(name, "matrix.%i.%i", p->bc_index, p->iter_num);
	    write_mps(lp_data, name);
	    return(ERROR__NUMERICAL_INSTABILITY);
	 }

       case LP_D_UNBOUNDED: /* the primal problem is infeasible */
       case LP_D_OBJLIM:
       case LP_OPTIMAL:
	 if (termcode == LP_D_UNBOUNDED){
	    PRINT(p->par.verbosity, 1, ("Feasibility lost -- "));
	 }else if ((p->has_ub && lp_data->objval > p->ub - p->par.granularity)
		   || termcode == LP_D_OBJLIM){
	    PRINT(p->par.verbosity, 1, ("Terminating due to high cost -- "));
	 }else{ /* optimal and not too high cost */
	    break;
	 }
	 comp_times->lp += used_time(&p->tt);
	 if (fathom(p, (termcode != LP_D_UNBOUNDED))){
	    comp_times->communication += used_time(&p->tt);
	    return(FUNCTION_TERMINATED_NORMALLY);
	 }else{
	    comp_times->communication += used_time(&p->tt);
	    continue;
	 }
      }

      /* If come to here, the termcode must have been OPTIMAL and the
       * cost cannot be too high. */
      /* is_feasible_u() fills up lp_data->x, too!! */
      if (is_feasible_u(p, FALSE, FALSE) == IP_FEASIBLE){
	 if (p->par.verbosity > 2){
	    printf ("Now displaying the feasible solution ...\n");
	    display_lp_solution_u(p, DISP_FEAS_SOLUTION);
	 }
	 cuts = -1;
      }else{

	 /*------------------------------------------------------------------*\
	  * send the current solution to the cut generator, and also to the
	  * cut pool if this is the 1st or cut_pool_check_freq-th iteration.
	 \*------------------------------------------------------------------*/

	 no_more_cuts_count = 0;
	 if (p->cut_pool &&
	     ((p->iter_num-1) % p->par.cut_pool_check_freq == 0) ){
	    no_more_cuts_count += send_lp_solution_u(p, p->cut_pool);
	 }
	 if (p->cut_gen){
	    no_more_cuts_count += send_lp_solution_u(p, p->cut_gen);
	 }

	 if (p->par.verbosity > 4){
	    printf ("Now displaying the relaxed solution ...\n");
	    display_lp_solution_u(p, DISP_RELAXED_SOLUTION);
	 }

	 comp_times->lp += used_time(&p->tt);

	 tighten_bounds(p);

	 comp_times->fixing += used_time(&p->tt);

	 cuts = 0;
	 if (p->cut_gen || p->cut_pool){
	    cuts = check_row_effectiveness(p);
	 }

	 /*------------------------------------------------------------------*\
	  * receive the cuts from the cut generator and the cut pool
	 \*------------------------------------------------------------------*/
         if ((cut_term = receive_cuts(p, TRUE, no_more_cuts_count)) >= 0){
            cuts += cut_term;
         }else{
            return(ERROR__USER);
         }
      }

      comp_times->lp += used_time(&p->tt);
      if (cuts < 0){ /* i.e. feasible solution is found */
	 if (fathom(p, TRUE)){
	    comp_times->communication += used_time(&p->tt);
	    return(FUNCTION_TERMINATED_NORMALLY);
	 }else{
	    comp_times->communication += used_time(&p->tt);
	    check_ub(p);
	    continue;
	 }
      }

      if (cuts == 0){
	 PRINT(p->par.verbosity, 2,
	       ("\nIn iteration %i ... no cuts were added.\n", p->iter_num));
      }else{
	 /* Go back to top */
	 PRINT(p->par.verbosity, 2,
	       ("\nIn iteration %i ... %i violated cuts were added.\n",
		p->iter_num, cuts));
	 continue;
      }

      comp_times->lp += used_time(&p->tt);

      /* So no cuts were found. Price out everything */
      new_cols = price_all_vars(p);
      new_vars = new_cols->num_vars + new_cols->rel_ub + new_cols->rel_lb;
      dual_feas = new_cols->dual_feas;
      free_col_set(&new_cols);
      comp_times->pricing += used_time(&p->tt);
      if (dual_feas != NOT_TDF)
	 break;

      /* Don't have total dual feasibility. The non-dual-feasible vars
       * have already been added. Go back and resolve. */
      PRINT(p->par.verbosity, 2,
	    ("%i variables added in price-out.\n", new_vars));
   }

   /* Now we know that we have TDF, just send back the node */
   comp_times->lp += used_time(&p->tt);
   send_node_desc(p, REPRICED_NODE);
   comp_times->communication += used_time(&p->tt);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int bfind(int key, int *table, int size)
{
   int i = 0, k = size;
   int j = size >> 1;   /* the element to be probed */
   while ( i < k ){
      if (table[j] == key){
	 return(j);
      }else if (table[j] < key){
	 i = j + 1;
      }else{
	 k = j;
      }
      j = (i + k) >> 1;
   }
   return(j-1); /* key is not found and it is between the (j-1)st and j-th */
}

/*===========================================================================*/

int collect_nonzeros(lp_prob *p, double *x, int *tind, double *tx)
{
   var_desc **vars = p->lp_data->vars;
   int n = p->lp_data->n;
   int i, cnt = 0;
   double lpetol = p->lp_data->lpetol;

   if (p->par.is_userind_in_order != TRUE) {
      colind_sort_extra(p);
      for (i = 0; i < n; i++){
         if (x[i] > lpetol || x[i] < -lpetol){
            tind[cnt] = vars[i]->userind;
            tx[cnt++] = x[i];
         }
      }
      /* order indices and values according to indices */
      qsort_id(tind, tx, cnt);
   } else {
      for (i = 0; i < n; i++){
         if (x[i] > lpetol || x[i] < -lpetol){
            tind[cnt] = i;
            tx[cnt++] = x[i];
         }
      }
   }
   return(cnt);
}

/*===========================================================================*/

int collect_fractions(lp_prob *p, double *x, int *tind, double *tx)
{
   var_desc **vars = p->lp_data->vars;
   int n = p->lp_data->n;
   int i, cnt = 0;
   double lpetol = p->lp_data->lpetol, xi;

   colind_sort_extra(p);
   for (i = 0; i < n; i++){
      xi = x[i];
      if (xi - floor(xi) > lpetol && ceil(xi) - xi > lpetol){
	 tind[cnt] = vars[i]->userind;
	 tx[cnt++] = x[i];
      }
   }
   /* order indices and values according to indices */
   qsort_id(tind, tx, cnt);
   return(cnt);
}

/*===========================================================================*/

node_desc *create_explicit_node_desc(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   int m = lp_data->m, n = lp_data->n;

   int bvarnum = p->base.varnum;
   var_desc **extravars = lp_data->vars + bvarnum;
   int extravarnum = n - bvarnum;

   int bcutnum = p->base.cutnum;
   row_data *rows = lp_data->rows;
   int extrarownum = m - bcutnum;
   int cutindsize;

   node_desc *desc = (node_desc *) calloc(1, sizeof(node_desc));

   /* Will need these anyway for basis */
   int *rstat = (int *) malloc(m * ISIZE);
   int *cstat = (int *) malloc(n * ISIZE);
   int *erstat = (extrarownum == 0) ? NULL : (int *) malloc(extrarownum*ISIZE);
   int *ecstat = (extravarnum == 0) ? NULL : (int *) malloc(extravarnum*ISIZE);

   int *ulist, *clist; /* this later uses tmp.i1 */
   int cutcnt, i, j;
#ifndef COMPILE_IN_LP
   int s_bufid, r_bufid;
#endif

   get_basis(lp_data, cstat, rstat);
   if (extrarownum > 0)
      memcpy(erstat, rstat + bcutnum, extrarownum * ISIZE);
   if (extravarnum > 0)
      memcpy(ecstat, cstat + bvarnum, extravarnum * ISIZE);

   /* To start with, send the non-indexed cuts (only those which will be
      saved) to the treemanager and ask for names */
   for (cutcnt = cutindsize = 0, i = bcutnum; i < m; i++){
      if ((rows[i].cut->branch & CUT_BRANCHED_ON) ||
	  !rows[i].free || (rows[i].free && rstat[i] != SLACK_BASIC)){
	 cutindsize++;
	 if (rows[i].cut->name < 0)
	    cutcnt++;
      }
   }
   if (cutcnt > 0){
#ifdef COMPILE_IN_LP
      row_data *tmp_rows = (row_data *) malloc(cutcnt*sizeof(row_data));
      
      for (j = 0, i = bcutnum; j < cutcnt; i++){
	 if (rows[i].cut->name < 0 &&
	     (!rows[i].free || (rows[i].free && rstat[i] != SLACK_BASIC)))
	    tmp_rows[j++] = rows[i];
      }
      unpack_cut_set(p->tm, 0, cutcnt, tmp_rows);
      FREE(tmp_rows);
#else
      s_bufid = init_send(DataInPlace);
      send_int_array(&cutcnt, 1);
      for (i = bcutnum; i < m; i++){
	 if (rows[i].cut->name < 0 &&
	     (!rows[i].free || (rows[i].free && rstat[i] != SLACK_BASIC)))
	    pack_cut(rows[i].cut);
      }
      send_msg(p->tree_manager, LP__CUT_NAMES_REQUESTED);
      freebuf(s_bufid);
#endif
   }

   /* create the uind list and the extravars basis description */
   desc->uind.type = EXPLICIT_LIST;
   desc->uind.added = 0;
   desc->uind.size = extravarnum;
   desc->basis.extravars.type = EXPLICIT_LIST;
   desc->basis.extravars.size = extravarnum;
   desc->basis.extravars.list = NULL;
   if (extravarnum > 0){
      desc->uind.list = ulist = (int *) malloc(extravarnum * ISIZE);
      desc->basis.extravars.stat = ecstat;
      for (i = extravarnum - 1; i >= 0; i--)
	 ulist[i] = extravars[i]->userind;
      if (lp_data->ordering == COLIND_ORDERED)
	 qsort_ii(ulist, ecstat, extravarnum);
   }else{
      desc->uind.list = NULL;
      desc->basis.extravars.stat = NULL;
   }
   /* create the basevars basis description */
   desc->basis.basevars.type = EXPLICIT_LIST;
   desc->basis.basevars.size = bvarnum;
   desc->basis.basevars.list = NULL;
   if (bvarnum)
      desc->basis.basevars.stat = cstat;
   else
      FREE(cstat);

   /* create the not_fixed list */
   desc->nf_status = lp_data->nf_status;
   if (desc->nf_status == NF_CHECK_AFTER_LAST ||
       desc->nf_status == NF_CHECK_UNTIL_LAST){
      desc->not_fixed.type = EXPLICIT_LIST;
      desc->not_fixed.added = 0;
      if ((desc->not_fixed.size = lp_data->not_fixed_num) > 0){
	 desc->not_fixed.list = (int *) malloc(desc->not_fixed.size * ISIZE);
	 memcpy(desc->not_fixed.list, lp_data->not_fixed,
		lp_data->not_fixed_num * ISIZE);
      }else{
	 desc->not_fixed.list = NULL;
      }
   }

#ifndef COMPILE_IN_LP
   /* At this point we will need the missing names */
   if (cutcnt > 0){
      static struct timeval tout = {15, 0};
      int *names = lp_data->tmp.i1; /* m */
      double start = wall_clock(NULL);
      do{
	 r_bufid = treceive_msg(p->tree_manager, LP__CUT_NAMES_SERVED, &tout);
	 if (! r_bufid){
	    if (pstat(p->tree_manager) != PROCESS_OK){
	       printf("TM has died -- LP exiting\n\n");
	       exit(-301);
	    }
	 }
      }while (! r_bufid);
      p->comp_times.idle_names += wall_clock(NULL) - start;
      receive_int_array(names, cutcnt);
      for (j = 0, i = bcutnum; j < cutcnt; i++){
	 if (rows[i].cut->name < 0 &&
	     (!rows[i].free || (rows[i].free && rstat[i] != SLACK_BASIC)))
	    rows[i].cut->name = names[j++];
      }
   }
#endif

   /* create the cutind list and the extrarows basis description */
   desc->cutind.type = EXPLICIT_LIST;
   desc->cutind.added = 0;
   desc->cutind.size = cutindsize;
   desc->basis.extrarows.type = EXPLICIT_LIST;
   desc->basis.extrarows.list = NULL;
   desc->basis.extrarows.size = cutindsize;
   if (cutindsize > 0){
      desc->cutind.list = clist = (int *) malloc(cutindsize * ISIZE);
      desc->basis.extrarows.stat = erstat;
      for (cutindsize = 0, i = bcutnum; i < m; i++){
	 if ((rows[i].cut->branch & CUT_BRANCHED_ON) ||
	     !rows[i].free || (rows[i].free && rstat[i] != SLACK_BASIC)){
	    clist[cutindsize] = rows[i].cut->name;
	    erstat[cutindsize++] = rstat[i];
	 }
      }
      qsort_ii(clist, erstat, cutindsize);
   }else{
      desc->cutind.list = NULL;
      desc->basis.extrarows.stat = NULL;
   }
   /* create the baserows basis description */
   desc->basis.baserows.type = EXPLICIT_LIST;
   desc->basis.baserows.size = bcutnum;
   desc->basis.baserows.list = NULL;
   if (bcutnum)
      desc->basis.baserows.stat = rstat;
   else
      FREE(rstat);

   /* Mark that there is a basis */
   desc->basis.basis_exists = TRUE;

   /* Add user description */
   add_to_desc_u(p, desc);

   return(desc);
}

/*===========================================================================*/

int check_tailoff(lp_prob *p)
{
   int gap_backsteps = p->par.tailoff_gap_backsteps;
   int obj_backsteps = p->par.tailoff_obj_backsteps;
   double *obj_hist = p->obj_history;

   int i;
   double sum, ub;
   int maxsteps = MAX(gap_backsteps, obj_backsteps);

   /*
   p->has_tailoff = FALSE;
   return (FALSE);
   */

   p->has_tailoff = TRUE;
   if (gap_backsteps >= 1 || obj_backsteps >= 2) {

      /* shift the data in obj_hist by one to the right and insert the
	 most recent objval to be the 0th */
      for (i = MIN(p->node_iter_num-1, maxsteps) - 1; i >= 0; i--) {
	 obj_hist[i+1] = obj_hist[i];
      }
      obj_hist[0] = p->lp_data->objval;

      if (p->bc_index == 0) {
         /*
          * root policy: generate cuts for min_root_cut_rounds and then stop.
          * if obj value doesnt improve in last
          * tailoff_max_no_impr_iters_root, then stop.
          */
         if (obj_hist[0] <= obj_hist[1] + p->lp_data->lpetol) {
            p->obj_no_impr_iters++;
         } else {
            p->obj_no_impr_iters = 0;
         }
         if (p->obj_no_impr_iters >=
               p->par.tailoff_max_no_impr_iters_root) {
            p->has_tailoff = TRUE;
            return (TRUE);
         }
         if (p->node_iter_num < p->par.min_root_cut_rounds) {
            p->has_tailoff = FALSE;
            return (FALSE);
         } else {
            p->has_tailoff = TRUE;
            return (TRUE);
         }
      }

      /* if there is an upper bound and we want gap based tailoff:
	 tailoff_gap is false if the average of the consecutive gap ratios is
	 less than gap_frac */
      if (p->node_iter_num>gap_backsteps && p->has_ub && gap_backsteps > 0) {
	 ub = p->ub;
	 for (i = 1, sum = 0; i <= gap_backsteps; i++) {
	    sum += (ub - obj_hist[i-1]) / (ub - obj_hist[i]);
	 }
	 if (sum / gap_backsteps > p->par.tailoff_gap_frac) {
	    PRINT(p->par.verbosity, 3, ("Branching because of tailoff in gap!\n"));
	    return(TRUE); /* there is tailoff */
	 }
      }

      /* if we want objective value based tailoff:
	 tailoff_obj is true if the average of the objective difference
	 ratios is smaller than par.tailoff_obj_frac */
      if (p->node_iter_num>obj_backsteps){
	 for (i = 2, sum = 0; i <= obj_backsteps; i++){
	    if (obj_hist[i-1] - obj_hist[i] > p->lp_data->lpetol){
	       sum += (obj_hist[i-2]-obj_hist[i-1]) / (obj_hist[i-1]-obj_hist[i]);
	    }else if (obj_hist[i-2] - obj_hist[i-1] > p->lp_data->lpetol){
	       sum += obj_backsteps;
	    }
	 }
	 if (sum / (obj_backsteps - 1) < p->par.tailoff_obj_frac){
	    PRINT(p->par.verbosity, 3, ("Branching because of tailoff in objective function!\n"));
	    PRINT(p->par.verbosity, 3, ("sum/n = %f, tailoff_obj_frac = %f\n",sum / (obj_backsteps - 1) , p->par.tailoff_obj_frac));
	    return(TRUE); /* there is tailoff */
	 }
      }

      /* Another check. All other checks seem to show that there is no
       * tailoff yet. 
       */
      if (p->node_iter_num>1 && 
	    obj_hist[0] - obj_hist[1] < p->par.tailoff_absolute){
	 PRINT(p->par.verbosity, 3, ("Branching because of tailoff in value of objective function!\n"));
	 return(TRUE);
      }

   } else {
      /* Both gap_backsteps and obj_backsteps are too small to procede with
         check_tailoff. The user asks for tailoff (since we came to this
	 function) yet doesn't want to check any kind of tailoff (since this
	 condition is true). Report no tailoff. */
      p->has_tailoff=FALSE;
      return(FALSE); /* no tailoff */
   }

   p->has_tailoff=FALSE;
   return(FALSE); /* gone thru everything ==> no tailoff */
}

/*===========================================================================*/

// Adapted from COIN's BRANCH AND CUT (CBC) solver! 

// See if rounding will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers 
// Fix values if asked for
// Returns 1 if solution, 0 if not

int round_solution(lp_prob *p, double *solutionValue, double *betterSolution)
{

  LPdata *lp_data = p->lp_data;
  int numberColumns = lp_data->n;
  //  int numberRows = lp_data->m; 
  int numberRows = p->base.cutnum + p->desc->cutind.size, nz = lp_data->nz;
  int returnCode = 0, numberIntegers = 0;
  double primalTolerance = lp_data->lpetol, integerTolerance = primalTolerance;
  double *lower, *upper, *rowLower, *rowUpper, *solution, *objective;
  double direction = p->mip->obj_sense == SYM_MINIMIZE ? 1: -1 ;
  double newSolutionValue = direction*lp_data->objval;
  double *element, *elementByRow;
  int * integerVariable, *isInteger;
  int *row, *column, *columnStart, *rowStart, *columnLength, *rowLength;
  int i, j;

  get_bounds(lp_data);
  get_x(lp_data);

  lower = lp_data->lb;
  upper = lp_data->ub;
  solution = lp_data->x;

  element = new double[nz];
  row = new int[nz];
  columnStart = new int[numberColumns+1];
  columnLength = new int[numberColumns];
  objective = new double[numberColumns];     

  elementByRow = new double[nz];
  column = new int[nz];
  rowStart = new int[numberRows+1];
  rowLength = new int[numberRows];
  rowUpper = new double[numberRows];
  rowLower = new double[numberRows];

  columnStart[0] = 0;
  rowStart[0] = 0;

  for (i = 0; i < numberColumns; i++){
     get_column(lp_data, i, &element[columnStart[i]], &row[columnStart[i]], 
		&columnLength[i], &objective[i]);     
     columnStart[i+1] = columnStart[i] + columnLength[i];

     for(j = 0; j < columnLength[i]; j++){
	if(row[columnStart[i] + j] >= numberRows){
	   columnLength[i] = j;
	   break;
	}
     }     
  }

  for (i = 0; i < numberRows; i++){
     get_row(lp_data, i, &elementByRow[rowStart[i]], &column[rowStart[i]],
	     &rowLength[i], &rowUpper[i], &rowLower[i]);
     rowStart[i+1] = rowStart[i] + rowLength[i];
  }	     

  isInteger = new int[numberColumns];
  integerVariable = new int[numberColumns];

  for (i = 0; i<numberColumns; i++){
     isInteger[i] = 0;
     if (lp_data->vars[i]->is_int){
	isInteger[i] = 1;
	integerVariable[numberIntegers++] = i;
     }
  }
 
  // Get solution array for heuristic solution

  double * newSolution = new double [numberColumns];
  memcpy(newSolution,solution,numberColumns*sizeof(double));

  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));
  for (i=0;i<numberColumns;i++) {
    int j;
    double value = newSolution[i];
    if (value) {
      for (j=columnStart[i];
	   j<columnStart[i]+columnLength[i];j++) {
	int iRow=row[j];
	//	printf("rowind %i: %i \n", j, iRow);
	//	if(j < 5){
	//	printf("element %i: %f \n", j, element[j]);
	//	}
	rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  for (i=0;i<numberRows;i++) {
    if(rowActivity[i]<rowLower[i]) {
      //assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
      rowActivity[i]=rowLower[i];
    } else if(rowActivity[i]>rowUpper[i]) {
      //assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
      rowActivity[i]=rowUpper[i];
    }
  }
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    double value=newSolution[iColumn];
    if (fabs(floor(value+0.5)-value)>integerTolerance) {
      double below = floor(value);
      double newValue=newSolution[iColumn];
      double cost = direction * objective[iColumn];
      double move;
      if (cost>0.0) {
	// try up
	move = 1.0 -(value-below);
      } else if (cost<0.0) {
	// try down
	move = below-value;
      } else {
	// won't be able to move unless we can grab another variable
	// just for now go down
	move = below-value;
      }
      newValue += move;
      newSolution[iColumn] = newValue;
      newSolutionValue += move*cost;
      int j;
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	int iRow = row[j];
	rowActivity[iRow] += move*element[j];
      }
    }
  }

  double penalty=0.0;
  
  // see if feasible
  for (i=0;i<numberRows;i++) {
    double value = rowActivity[i];
    double thisInfeasibility=0.0;
    if (value<rowLower[i]-primalTolerance)
      thisInfeasibility = value-rowLower[i];
    else if (value>rowUpper[i]+primalTolerance)
      thisInfeasibility = value-rowUpper[i];
    if (thisInfeasibility) {
      // See if there are any slacks I can use to fix up
      // maybe put in coding for multiple slacks?
      double bestCost = 1.0e50;
      int k;
      int iBest=-1;
      double addCost=0.0;
      double newValue=0.0;
      double changeRowActivity=0.0;
      double absInfeasibility = fabs(thisInfeasibility);
      for (k=rowStart[i];k<rowStart[i]+rowLength[i];k++) {
	int iColumn = column[k];
	if (columnLength[iColumn]==1) {
	  double currentValue = newSolution[iColumn];
	  double elementValue = elementByRow[k];
	  double lowerValue = lower[iColumn];
	  double upperValue = upper[iColumn];
	  double gap = rowUpper[i]-rowLower[i];
	  double absElement=fabs(elementValue);
	  if (thisInfeasibility*elementValue>0.0) {
	    // we want to reduce
	    if ((currentValue-lowerValue)*absElement>=absInfeasibility) {
	      // possible - check if integer
	      double distance = absInfeasibility/absElement;
	      double thisCost = -direction*objective[iColumn]*distance;
	      if (isInteger[iColumn]) {
		distance = ceil(distance-primalTolerance);
		if (currentValue-distance>=lowerValue-primalTolerance) {
		  if (absInfeasibility-distance*absElement< -gap-primalTolerance)
		    thisCost=1.0e100; // no good
		  else
		    thisCost = -direction*objective[iColumn]*distance;
		} else {
		  thisCost=1.0e100; // no good
		}
	      }
	      if (thisCost<bestCost) {
		bestCost=thisCost;
		iBest=iColumn;
		addCost = thisCost;
		newValue = currentValue-distance;
		changeRowActivity = -distance*elementValue;
	      }
	    }
	  } else {
	    // we want to increase
	    if ((upperValue-currentValue)*absElement>=absInfeasibility) {
	      // possible - check if integer
	      double distance = absInfeasibility/absElement;
	      double thisCost = direction*objective[iColumn]*distance;
	      if (isInteger[iColumn]) {
		distance = ceil(distance-primalTolerance);
		//assert (currentValue-distance<=upperValue+primalTolerance);
		if (absInfeasibility-distance*absElement< -gap-primalTolerance)
		  thisCost=1.0e100; // no good
		else
		  thisCost = direction*objective[iColumn]*distance;
	      }
	      if (thisCost<bestCost) {
		bestCost=thisCost;
		iBest=iColumn;
		addCost = thisCost;
		newValue = currentValue+distance;
		changeRowActivity = distance*elementValue;
	      }
	    }
	  }
	}
      }
      if (iBest>=0) {
	/*printf("Infeasibility of %g on row %d cost %g\n",
	  thisInfeasibility,i,addCost);*/
	newSolution[iBest]=newValue;
	thisInfeasibility=0.0;
	newSolutionValue += addCost;
	rowActivity[i] += changeRowActivity;
      }
      penalty += fabs(thisInfeasibility);
    }
  }

  // Could also set SOS (using random) and repeat
  if (!penalty) {
    // See if we can do better
    //seed_++;
    //CoinSeedRandom(seed_);
    // Random number between 0 and 1.
    double randomNumber = CoinDrand48();
    int iPass;
    int start[2];
    int end[2];
    int iRandom = (int) (randomNumber*((double) numberIntegers));
    start[0]=iRandom;
    end[0]=numberIntegers;
    start[1]=0;
    end[1]=iRandom;
    for (iPass=0;iPass<2;iPass++) {
      int i;
      for (i=start[iPass];i<end[iPass];i++) {
	int iColumn = integerVariable[i];
	//double value=newSolution[iColumn];
	//assert (fabs(floor(value+0.5)-value)<integerTolerance);
	double cost = direction * objective[iColumn];
	double move=0.0;
	if (cost>0.0)
	  move = -1.0;
	else if (cost<0.0)
	  move=1.0;
	while (move) {
	  bool good=true;
	  double newValue=newSolution[iColumn]+move;
	  if (newValue<lower[iColumn]-primalTolerance||
	      newValue>upper[iColumn]+primalTolerance) {
	    move=0.0;
	  } else {
	    // see if we can move
	    int j;
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      double newActivity = rowActivity[iRow] + move*element[j];
	      if (newActivity<rowLower[iRow]-primalTolerance||
		  newActivity>rowUpper[iRow]+primalTolerance) {
		good=false;
		break;
	      }
	    }
	    if (good) {
	      newSolution[iColumn] = newValue;
	      newSolutionValue += move*cost;
	      int j;
	      for (j=columnStart[iColumn];
		   j<columnStart[iColumn]+columnLength[iColumn];j++) {
		int iRow = row[j];
		rowActivity[iRow] += move*element[j];
	      }
	    } else {
	      move=0.0;
	    }
	  }
	}
      }
    }
    if (newSolutionValue < *solutionValue) {
      // paranoid check
      memset(rowActivity,0,numberRows*sizeof(double));
      for (i=0;i<numberColumns;i++) {
	int j;
	double value = newSolution[i];
	if (value) {
	  for (j=columnStart[i];
	       j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    rowActivity[iRow] += value*element[j];
	  }
	}
      }
      // check was approximately feasible
      bool feasible=true;
      for (i=0;i<numberRows;i++) {
	if(rowActivity[i]<rowLower[i]) {
	  if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	    feasible = false;
	} else if(rowActivity[i]>rowUpper[i]) {
	  if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
	    feasible = false;
	}
      }
      if (feasible) {
	// new solution
	memcpy(betterSolution, newSolution, numberColumns*DSIZE);
	*solutionValue = newSolutionValue;
	//printf("** Solution of %g found by rounding\n",newSolutionValue);
	returnCode=1;
      } else {
	// Can easily happen
	//printf("Debug CbcRounding giving bad solution\n");
      }
    }
  }
  delete [] integerVariable;
  delete [] isInteger;

  delete [] element;
  delete [] row;
  delete [] columnStart;
  delete [] columnLength;
  delete [] objective;

  delete [] elementByRow;
  delete [] column;
  delete [] rowStart;
  delete [] rowLength;
  delete [] rowUpper;
  delete [] rowLower;

  delete [] newSolution;
  delete [] rowActivity;
  return returnCode;
}

/*===========================================================================*/

// Adapted from COIN's BRANCH AND CUT (CBC) solver! 

/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */

int local_search(lp_prob *p, double *solutionValue, double *colSolution,
		 double *betterSolution)
{
 
  LPdata *lp_data = p->lp_data;
  int numberColumns = lp_data->n;
  int numberRows = p->base.cutnum + p->desc->cutind.size, nz = lp_data->nz;
  int returnCode = 0, numberIntegers = 0;
  double primalTolerance = lp_data->lpetol;
  double *rowLower, *rowUpper, *solution = colSolution, *objective;
  double direction = p->mip->obj_sense == SYM_MINIMIZE ? 1: -1 ;
  double newSolutionValue = direction*p->ub;
  double *element, *elementByRow;
  int * integerVariable, *isInteger;
  int *row, *columnStart, *columnLength, *column, rowLength;
  int i, j;
  
  element = new double[nz];
  row = new int[nz];
  columnStart = new int[numberColumns+1];
  columnLength = new int[numberColumns];
  objective = new double[numberColumns];     

  rowUpper = new double[numberRows];
  rowLower = new double[numberRows];


  elementByRow = new double[numberColumns];
  column = new int[numberColumns];
  
  columnStart[0] = 0;

  for (i = 0; i < numberColumns; i++){
     get_column(lp_data, i, &element[columnStart[i]], &row[columnStart[i]], 
		&columnLength[i], &objective[i]);     
     columnStart[i+1] = columnStart[i] + columnLength[i];

     for(j = 0; j < columnLength[i]; j++){
	if(row[columnStart[i] + j] >= numberRows){
	   columnLength[i] = j;
	   break;
	}
     }     
  }

  for (i = 0; i < numberRows; i++){
     get_row(lp_data, i, elementByRow, column, &rowLength, &rowUpper[i], 
	     &rowLower[i]);
  }

  isInteger = new int[numberColumns];
  integerVariable = new int[numberColumns];

  for (i = 0; i<numberColumns; i++){
     isInteger[i] = 0;
     if (lp_data->vars[i]->is_int){
	isInteger[i] = 1;
	integerVariable[numberIntegers++] = i;
     }
  }

  // Column copy
  /* 
  const double * element = matrix.getElements();
  const int * row = matrix.getIndices();
  const CoinBigIndex * columnStart = matrix.getVectorStarts();
  const int * columnLength = matrix.getVectorLengths();
  */

  // Get solution array for heuristic solution
  double * newSolution = new double [numberColumns];
  memcpy(newSolution,solution,numberColumns*sizeof(double));

  // way is 1 if down possible, 2 if up possible, 3 if both possible
  int * way = new int[numberIntegers];
  // corrected costs
  double * cost = new double[numberIntegers];
  // for array to mark infeasible rows after iColumn branch
  char * mark = new char[numberRows];
  memset(mark,0,numberRows);
  // space to save values so we don't introduce rounding errors
  double * save = new double[numberRows];

  // clean solution
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    
    // get original bounds
    double originalLower = lp_data->vars[iColumn]->lb; //p->mip->lb[iColumn];
    double originalUpper = lp_data->vars[iColumn]->ub; //p->mip->ub[iColumn];

    //  double originalLower = lp_data->lb[iColumn];
    //double originalUpper = lp_data->ub[iColumn];

    //   double originalLower = p->mip->lb[iColumn];
    //  double originalUpper = p->mip->ub[iColumn];

    double value=newSolution[iColumn];
    double nearest=floor(value+0.5);
    //assert(fabs(value-nearest)<10.0*primalTolerance);
    value=nearest;
    newSolution[iColumn]=nearest;
    // if away from lower bound mark that fact
    if (nearest>originalLower) {
      //      used_[iColumn]=1;
    }
    cost[i] = direction*objective[iColumn];
    int iway=0;
    
    if (value>originalLower+0.5) 
      iway = 1;
    if (value<originalUpper-0.5) 
      iway |= 2;
    way[i]=iway;
  }
  // get row activities
  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));

  for (i=0;i<numberColumns;i++) {
    int j;
    double value = newSolution[i];
    if (value) {
      for (j=columnStart[i];
	   j<columnStart[i]+columnLength[i];j++) {
	int iRow=row[j];
	rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  // if very infeasible then give up
  bool tryHeuristic=true;
  for (i=0;i<numberRows;i++) {
    if(rowActivity[i]<rowLower[i]) {
      if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowLower[i];
    } else if(rowActivity[i]>rowUpper[i]) {
      if (rowActivity[i]<rowUpper[i]+10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowUpper[i];
    }
  }
  if (tryHeuristic) {
    
    // best change in objective
    double bestChange=0.0;
    
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      
      double objectiveCoefficient = cost[i];
      int k;
      int j;
      int goodK=-1;
      int wayK=-1,wayI=-1;
      if ((way[i]&1)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] -= element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try down
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (-objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=-1;
		bestChange = -objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (-objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=-1;
		bestChange = -objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if ((way[i]&2)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] += element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try up
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=1;
		bestChange = objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=1;
		bestChange = objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if (goodK>=0) {
	// we found something - update solution
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayI * element[j];
	}
	newSolution[iColumn] += wayI;
	int kColumn = integerVariable[goodK];
	for (j=columnStart[kColumn];
	     j<columnStart[kColumn]+columnLength[kColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayK * element[j];
	}
	newSolution[kColumn] += wayK;
	// See if k can go further ?
	// get original bounds
	double originalLower = p->mip->lb[kColumn];
	double originalUpper = p->mip->ub[kColumn];
	
	double value=newSolution[kColumn];
	int iway=0;
	if (value>originalLower+0.5) 
	  iway = 1;
	if (value<originalUpper-0.5) 
	  iway |= 2;
	way[goodK]=iway;
      }
    }
    if (bestChange+newSolutionValue<*solutionValue) {
      // new solution
      memcpy(betterSolution, newSolution, numberColumns*DSIZE);
      returnCode=1;
      *solutionValue = newSolutionValue + bestChange;
      if (bestChange>1.0e-12)
	printf("Local search heuristic improved solution by %g\n",
	     -bestChange);
      // paranoid check
      memset(rowActivity,0,numberRows*sizeof(double));
      
      for (i=0;i<numberColumns;i++) {
	int j;
	double value = newSolution[i];
	if (value) {
	  for (j=columnStart[i];
	       j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    rowActivity[iRow] += value*element[j];
	  }
	}
      }
      // check was approximately feasible
      for (i=0;i<numberRows;i++) {
	if(rowActivity[i]<rowLower[i]) {
	   //assert (rowActivity[i]>rowLower[i]-10.0*primalTolerance);
	} else if(rowActivity[i]>rowUpper[i]) {
	   //assert (rowActivity[i]<rowUpper[i]+10.0*primalTolerance);
	}
      }
      for (i=0;i<numberIntegers;i++) {
	int iColumn = integerVariable[i];
	double originalLower = p->mip->lb[iColumn];
	//double originalUpper = integerObject->originalUpperBound();

	double value=newSolution[iColumn];
	// if away from lower bound mark that fact
	if (value>originalLower) {
	  //	  used_[iColumn]=1;
	}
      }
    }
  }


  delete [] integerVariable;
  delete [] isInteger;

  delete [] element;
  delete [] row;
  delete [] columnStart;
  delete [] columnLength;
  delete [] objective;

  delete [] elementByRow;
  delete [] column;
  delete [] rowUpper;
  delete [] rowLower;

  delete [] newSolution;
  delete [] rowActivity;
  delete [] way;
  delete [] cost;
  delete [] save;
  delete [] mark;

  return returnCode;
}

/*===========================================================================*/

void lp_exit(lp_prob *p)
{
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_msg(p->tree_manager, SOMETHING_DIED);
   freebuf(s_bufid);
   comm_exit();
   exit(-1);
}

/*===========================================================================*/

void lp_close(lp_prob *p)
{
#ifndef COMPILE_IN_LP
   int s_bufid;
   
   /* Send back the timing data for the whole algorithm */
   s_bufid = init_send(DataInPlace);
   send_char_array((char *)&p->comp_times, sizeof(node_times));
   send_msg(p->tree_manager, LP__TIMING);
   freebuf(s_bufid);
#else
#pragma omp critical (timing_update)
{
   p->tm->comp_times.communication    += p->comp_times.communication;
   p->tm->comp_times.lp               += p->comp_times.lp;
   p->tm->comp_times.separation       += p->comp_times.separation;
   p->tm->comp_times.fixing           += p->comp_times.fixing;
   p->tm->comp_times.pricing          += p->comp_times.pricing;
   p->tm->comp_times.strong_branching += p->comp_times.strong_branching;
   p->tm->comp_times.fp               += p->comp_times.fp;
   p->tm->comp_times.primal_heur      += p->comp_times.primal_heur;

   p->tm->comp_times.cuts             += p->comp_times.cuts;
   p->tm->comp_times.gomory_cuts      += p->comp_times.gomory_cuts;
   p->tm->comp_times.knapsack_cuts    += p->comp_times.knapsack_cuts;
   p->tm->comp_times.oddhole_cuts     += p->comp_times.oddhole_cuts;
   p->tm->comp_times.clique_cuts      += p->comp_times.clique_cuts;
   p->tm->comp_times.probing_cuts     += p->comp_times.probing_cuts;
   p->tm->comp_times.mir_cuts         += p->comp_times.mir_cuts;
   p->tm->comp_times.twomir_cuts      += p->comp_times.twomir_cuts;
   p->tm->comp_times.rounding_cuts    += p->comp_times.rounding_cuts;
   p->tm->comp_times.landp_cuts       += p->comp_times.landp_cuts;
   p->tm->comp_times.flowcover_cuts   += p->comp_times.flowcover_cuts;
   p->tm->comp_times.lift_and_project_cuts += 
      p->comp_times.lift_and_project_cuts;
   p->tm->comp_times.redsplit_cuts += p->comp_times.redsplit_cuts;
   p->tm->comp_times.dupes_and_bad_coeffs_in_cuts += 
      p->comp_times.dupes_and_bad_coeffs_in_cuts;

   p->tm->lp_stat.lp_calls                += p->lp_stat.lp_calls;
   p->tm->lp_stat.str_br_lp_calls         += p->lp_stat.str_br_lp_calls;
   p->tm->lp_stat.lp_sols                 += p->lp_stat.lp_sols;
   p->tm->lp_stat.str_br_bnd_changes      += p->lp_stat.str_br_bnd_changes;
   p->tm->lp_stat.str_br_nodes_fathomed   += p->lp_stat.str_br_nodes_fathomed;

   p->tm->lp_stat.cuts_generated          += p->lp_stat.cuts_generated;
   p->tm->lp_stat.gomory_cuts             += p->lp_stat.gomory_cuts;
   p->tm->lp_stat.knapsack_cuts           += p->lp_stat.knapsack_cuts;
   p->tm->lp_stat.oddhole_cuts            += p->lp_stat.oddhole_cuts;
   p->tm->lp_stat.clique_cuts             += p->lp_stat.clique_cuts;
   p->tm->lp_stat.probing_cuts            += p->lp_stat.probing_cuts;
   p->tm->lp_stat.mir_cuts                += p->lp_stat.mir_cuts;
   p->tm->lp_stat.twomir_cuts             += p->lp_stat.twomir_cuts;
   p->tm->lp_stat.rounding_cuts           += p->lp_stat.rounding_cuts;
   p->tm->lp_stat.landp_cuts              += p->lp_stat.landp_cuts;
   p->tm->lp_stat.flowcover_cuts          += p->lp_stat.flowcover_cuts;
   p->tm->lp_stat.lift_and_project_cuts   += p->lp_stat.lift_and_project_cuts;
   p->tm->lp_stat.redsplit_cuts           += p->lp_stat.redsplit_cuts;

   p->tm->lp_stat.cuts_root               += p->lp_stat.cuts_root;
   p->tm->lp_stat.gomory_cuts_root        += p->lp_stat.gomory_cuts_root;
   p->tm->lp_stat.knapsack_cuts_root      += p->lp_stat.knapsack_cuts_root;
   p->tm->lp_stat.oddhole_cuts_root       += p->lp_stat.oddhole_cuts_root;
   p->tm->lp_stat.clique_cuts_root        += p->lp_stat.clique_cuts_root;
   p->tm->lp_stat.probing_cuts_root       += p->lp_stat.probing_cuts_root;
   p->tm->lp_stat.mir_cuts_root           += p->lp_stat.mir_cuts_root;
   p->tm->lp_stat.twomir_cuts_root        += p->lp_stat.twomir_cuts_root;
   p->tm->lp_stat.rounding_cuts_root      += p->lp_stat.rounding_cuts_root;
   p->tm->lp_stat.landp_cuts_root         += p->lp_stat.landp_cuts_root;
   p->tm->lp_stat.flowcover_cuts_root     += p->lp_stat.flowcover_cuts_root;
   p->tm->lp_stat.lift_and_project_cuts_root += 
      p->lp_stat.lift_and_project_cuts_root;
   p->tm->lp_stat.redsplit_cuts_root += 
      p->lp_stat.redsplit_cuts_root;

   p->tm->lp_stat.num_poor_cuts           += p->lp_stat.num_poor_cuts;
   p->tm->lp_stat.num_duplicate_cuts      += p->lp_stat.num_duplicate_cuts;
   p->tm->lp_stat.num_unviolated_cuts     += p->lp_stat.num_unviolated_cuts;
   p->tm->lp_stat.cuts_deleted_from_lps   += p->lp_stat.cuts_deleted_from_lps;
   p->tm->lp_stat.cuts_added_to_lps       += p->lp_stat.cuts_added_to_lps;

   p->tm->lp_stat.gomory_calls            += p->lp_stat.gomory_calls;
   p->tm->lp_stat.knapsack_calls          += p->lp_stat.knapsack_calls;
   p->tm->lp_stat.oddhole_calls           += p->lp_stat.oddhole_calls;
   p->tm->lp_stat.clique_calls            += p->lp_stat.clique_calls;
   p->tm->lp_stat.probing_calls           += p->lp_stat.probing_calls;
   p->tm->lp_stat.mir_calls               += p->lp_stat.mir_calls;
   p->tm->lp_stat.twomir_calls            += p->lp_stat.twomir_calls;
   p->tm->lp_stat.rounding_calls          += p->lp_stat.rounding_calls;
   p->tm->lp_stat.landp_calls             += p->lp_stat.landp_calls;
   p->tm->lp_stat.flowcover_calls         += p->lp_stat.flowcover_calls;
   p->tm->lp_stat.lift_and_project_calls  += p->lp_stat.lift_and_project_calls;
   p->tm->lp_stat.redsplit_calls          += p->lp_stat.redsplit_calls;

   p->tm->lp_stat.fp_calls                += p->lp_stat.fp_calls;
   p->tm->lp_stat.fp_lp_calls             += p->lp_stat.fp_lp_calls;
   p->tm->lp_stat.fp_num_sols             += p->lp_stat.fp_num_sols;
}
#endif
#ifdef COMPILE_IN_CG
   cg_close(p->cgp);
#endif
#ifndef COMPILE_IN_TM
   free_lp(p);
#endif
}

/*===========================================================================*/
/*
 * save the changes in bounds that occurred while processing the current node
 * into current-node's node_desc. These changes are available by comparing
 * vars[i]->lb and vars[i]->new_lb etc. After saving the changes, vars[i]->lb,
 * vars[i]->ub are changed to new_lb and new_ub so that the same changes are
 * not saved in the child-node's desc.
 */
int add_bound_changes_to_desc(node_desc *desc, lp_prob *p)
{
#ifdef COMPILE_IN_LP
   LPdata                *lp_data = p->lp_data; 
   var_desc             **vars = lp_data->vars;
   int                    i, num_bnd_changes, cnt;
   bounds_change_desc    *bnd_change;
   int                   *index;
   char                  *lbub;
   double                *value;

   num_bnd_changes = 0;
   for (i=0;i<lp_data->n;i++) {
      if (vars[i]->new_lb>vars[i]->lb) {
         num_bnd_changes++;
      }
      if (vars[i]->new_ub<vars[i]->ub) {
         num_bnd_changes++;
      }
   }
   if (num_bnd_changes>0) {
      bnd_change = desc->bnd_change = (bounds_change_desc *) 
         calloc (1, sizeof(bounds_change_desc));
      bnd_change->num_changes = num_bnd_changes;
      index = bnd_change->index = (int *)malloc(num_bnd_changes*ISIZE);
      lbub  = bnd_change->lbub = (char *)malloc(num_bnd_changes*CSIZE);
      value = bnd_change->value = (double *)malloc(num_bnd_changes*DSIZE);
      cnt = 0;
      for (i=0;i<lp_data->n;i++) {
         if (vars[i]->new_lb>vars[i]->lb) {
            index[cnt] = vars[i]->userind;
            lbub[cnt] = 'L';
            value[cnt] = vars[i]->new_lb;
            cnt++;
            vars[i]->lb = vars[i]->new_lb;
         }
         if (vars[i]->new_ub<vars[i]->ub) {
            index[cnt] = vars[i]->userind;
            lbub[cnt] = 'U';
            value[cnt] = vars[i]->new_ub;
            cnt++;
            vars[i]->ub = vars[i]->new_ub;
         }
      }
   } else {
      desc->bnd_change = NULL;
   }
#endif

   return 0;
}

/*===========================================================================*/
int str_br_bound_changes(lp_prob *p, int num_bnd_changes, double *bnd_val, 
      int *bnd_ind, char *bnd_sense)
{
#ifdef COMPILE_IN_LP
   bounds_change_desc    *bnd_change;
   int                   i, j;
   var_desc              **vars = p->lp_data->vars;
   int                   *index;
   double                *value;
   char                  *lbub;

   if (num_bnd_changes<1) {
      return 0;
   }
   if (p->tm->active_nodes[p->proc_index]->desc.bnd_change == NULL) {
      bnd_change = (bounds_change_desc *)calloc(1, sizeof(bounds_change_desc));
      index = bnd_change->index = (int *)malloc(num_bnd_changes*ISIZE);
      lbub = bnd_change->lbub = (char *)malloc(num_bnd_changes*CSIZE);
      value = bnd_change->value = (double *)malloc(num_bnd_changes*DSIZE);
      bnd_change->num_changes = num_bnd_changes;
      j = 0;
   } else {
      bnd_change = p->tm->active_nodes[p->proc_index]->desc.bnd_change;
      j = bnd_change->num_changes;
      bnd_change->num_changes += num_bnd_changes;
      index = bnd_change->index = (int *)realloc(bnd_change->index, 
            bnd_change->num_changes*ISIZE);
      lbub = bnd_change->lbub = (char *)realloc(bnd_change->lbub,
            bnd_change->num_changes*CSIZE);
      value = bnd_change->value = (double *)realloc(bnd_change->value,
            bnd_change->num_changes*DSIZE);
   }
   for (i = 0; i<num_bnd_changes; i++) {
      index[i+j] = vars[bnd_ind[i]]->userind;
      lbub[i+j] = (bnd_sense[i] == 'L') ? 'U' : 'L';
      value[i+j] = bnd_val[i];
   }
   p->tm->active_nodes[p->proc_index]->desc.bnd_change = bnd_change;

#endif
   return 0;
}

/*===========================================================================*/
/* this function is called after root node has been processed. we update
 * frequency of cut generation for different cuts depending upon how many cuts
 * were generated and how much time was used
 */
int update_cut_parameters(lp_prob *p)
{
#ifdef USE_CGL_CUTS
   /* TODO: check (a) time (b) if any cuts are actually in the LP */
   lp_stat_desc  lp_stat  = p->lp_stat;
   cgl_params   *par      = &(p->par.cgl);
   cgl_params   *data_par = &(p->lp_data->cgl);
   /* probing cuts */
   if (par->generate_cgl_probing_cuts == GENERATE_IF_IN_ROOT && 
       lp_stat.probing_cuts_root<1) {
      par->generate_cgl_probing_cuts_freq = -1;
   }
   if (par->generate_cgl_probing_cuts == GENERATE_DEFAULT) {
      if (lp_stat.probing_cuts_root<1) {
         data_par->generate_cgl_probing_cuts_freq = 
              par->generate_cgl_probing_cuts_freq = -1;
      } else {
         data_par->generate_cgl_probing_cuts_freq = 
              par->generate_cgl_probing_cuts_freq = 100;
      }
   }

   /* cliques cuts */
   if (par->generate_cgl_clique_cuts == GENERATE_IF_IN_ROOT && 
       lp_stat.clique_cuts_root<1) {
      par->generate_cgl_clique_cuts_freq = -1;
   }
   if (par->generate_cgl_probing_cuts == GENERATE_DEFAULT) {
      if (lp_stat.probing_cuts_root<1) {
         data_par->generate_cgl_clique_cuts_freq = 
              par->generate_cgl_clique_cuts_freq = 200;
      } else {
         data_par->generate_cgl_clique_cuts_freq = 
              par->generate_cgl_clique_cuts_freq = 10;
      }
   }

   /* flow and cover cuts */
   if (par->generate_cgl_flowcover_cuts == GENERATE_IF_IN_ROOT && 
       lp_stat.flowcover_cuts_root<1) {
      par->generate_cgl_flowcover_cuts_freq = -1;
   }
   if (par->generate_cgl_flowcover_cuts == GENERATE_DEFAULT) {
      if (lp_stat.flowcover_cuts_root<1) {
         data_par->generate_cgl_flowcover_cuts_freq = 
              par->generate_cgl_flowcover_cuts_freq = -1;
      } else {
         data_par->generate_cgl_flowcover_cuts_freq = 
              par->generate_cgl_flowcover_cuts_freq = 100;
      }
   }

   /* knapsack cuts */
   if (par->generate_cgl_knapsack_cuts == GENERATE_IF_IN_ROOT && 
       lp_stat.knapsack_cuts_root<1) {
      par->generate_cgl_knapsack_cuts_freq = -1;
   }
   if (par->generate_cgl_knapsack_cuts == GENERATE_DEFAULT) {
      if (lp_stat.knapsack_cuts_root<1) {
         data_par->generate_cgl_knapsack_cuts_freq = 
              par->generate_cgl_knapsack_cuts_freq = 200;
      } else {
         data_par->generate_cgl_knapsack_cuts_freq = 
              par->generate_cgl_knapsack_cuts_freq = 20;
      }
   }

   /* gomory cuts */
   if (par->generate_cgl_gomory_cuts == GENERATE_IF_IN_ROOT && 
       lp_stat.gomory_cuts_root<1) {
      par->generate_cgl_gomory_cuts_freq = -1;
   }
   if (par->generate_cgl_gomory_cuts == GENERATE_DEFAULT) {
      if (lp_stat.gomory_cuts_root<1) {
         data_par->generate_cgl_gomory_cuts_freq = 
              par->generate_cgl_gomory_cuts_freq = 100;
      } else {
         data_par->generate_cgl_gomory_cuts_freq = 
              par->generate_cgl_gomory_cuts_freq = 10;
      }
   }
#endif
   return 0;
}

/*===========================================================================*/
int generate_cgl_cuts_new(lp_prob *p, int *num_cuts, cut_data ***cuts, 
      int send_to_pool, int *bound_changes)
{
#ifdef USE_CGL_CUTS
   int *should_generate = (int *) malloc(CGL_NUM_GENERATORS*ISIZE);
   int i, should_stop = FALSE, repeat_with_long = TRUE, max_cut_length;
   OsiCuts cutlist;
   const int n                 = p->lp_data->n;
   OsiXSolverInterface  *si    = p->lp_data->si;
   var_desc             **vars = p->lp_data->vars;
   int                  was_tried = FALSE;
   
   if (p->iter_num < 2) {
      for (i = 0; i < n; i++) {
         if (vars[i]->is_int) { // integer or binary
            si->setInteger(i);
         }
      }  
   }

   max_cut_length = p->par.max_cut_length;
   if (p->par.tried_long_cuts == TRUE) {
      repeat_with_long = FALSE;
   }
   for (i=0; i<CGL_NUM_GENERATORS; i++) {
      generate_cgl_cut_of_type(p, i, &cutlist, &was_tried);
      check_and_add_cgl_cuts(p, i, cuts, num_cuts, bound_changes, &cutlist, 
            send_to_pool);
      should_stop_adding_cgl_cuts(p, i, &should_stop);
      if (should_stop == TRUE) {
         break;
      }
      if (i==CGL_NUM_GENERATORS-1 && p->bc_index < 1 && *num_cuts < 1 && 
            repeat_with_long == TRUE) {
         p->par.max_cut_length = 1000;
         i = 0;
         repeat_with_long = FALSE;
         p->par.tried_long_cuts = TRUE;
      }
   }
   p->par.max_cut_length = max_cut_length;

   if (was_tried == TRUE && p->bc_index > 0) {
      p->lp_stat.num_cut_iters_in_path++;
   }

   FREE(should_generate);
#endif
   return 0;
}

/*===========================================================================*/
int should_use_cgl_generator(lp_prob *p, int *should_generate, 
      int which_generator, void *generator)
{

   int bc_index = p->bc_index;
   int is_root_node = (bc_index < 1) ? TRUE : FALSE;
#ifdef USE_CGL_CUTS
   *should_generate = FALSE;
   switch (which_generator) {
    case CGL_PROBING_GENERATOR:
      {
         CglProbing *probing = (CglProbing *)generator;
         int param = p->par.cgl.generate_cgl_probing_cuts;
         int freq  = p->par.cgl.generate_cgl_probing_cuts_freq;
         if (param < 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT && (freq < 1 || 
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_ONLY_IN_ROOT && bc_index > 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_IF_IN_ROOT && (freq < 1 || 
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_PERIODICALLY && (freq < 1 || 
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT) {
            if (bc_index > 0) {
               if (p->comp_times.probing_cuts > p->comp_times.lp/5) {
                  *should_generate = FALSE;
                  break;
               }
            } else {
               if (p->lp_stat.probing_cuts > p->lp_stat.cuts_generated/2
                   && p->comp_times.probing_cuts > 3*p->comp_times.lp) {
                  p->par.cgl.probing_is_expensive = TRUE;
                  *should_generate = FALSE;
                  break;
               } else if (p->lp_stat.probing_cuts <= 
                     p->lp_stat.cuts_generated/2 && 
                     p->comp_times.probing_cuts > 3*p->comp_times.lp) {
                  p->par.cgl.probing_is_expensive = TRUE;
                  *should_generate = FALSE;
                  break;
               }
            }
         }

         probing->setRowCuts(3); 
         probing->setMode(2);
         if (p->has_ub) {
            probing->setUsingObjective(1);
         }
         if (is_root_node == TRUE && !p->par.cgl.probing_is_expensive) {
            probing->setMaxPass(10); /* default is 3 */
            probing->setMaxPassRoot(10); /* default is 3 */
            probing->setMaxElements(10000);  /* default is 1000 */
            probing->setMaxElementsRoot(10000); /* default is 10000 */
            probing->setMaxLook(500);    /* default is 50 */
            probing->setMaxLookRoot(500);    /* default is 50 */
            probing->setMaxProbe(200);   /* default is 100 */
            probing->setMaxProbeRoot(200);   /* default is 100 */
         }
         *should_generate = TRUE;
         p->lp_stat.probing_calls++;
         break;
      }
    case CGL_CLIQUE_GENERATOR:
      {
         CglClique *clique = (CglClique *)generator;
         int param = p->par.cgl.generate_cgl_clique_cuts;
         int freq  = p->par.cgl.generate_cgl_clique_cuts_freq;
         if (param < 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT && (freq < 0 ||
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_ONLY_IN_ROOT && bc_index > 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_IF_IN_ROOT && (freq < 0 ||
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_PERIODICALLY && (freq < 0 ||
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } 
         *should_generate = TRUE;
         clique->setStarCliqueReport(FALSE);
         clique->setRowCliqueReport(FALSE);
         p->lp_stat.clique_calls++;
         break;
      }
    case CGL_KNAPSACK_GENERATOR:
      {
         CglKnapsackCover *knapsack = (CglKnapsackCover *)generator;
         int param = p->par.cgl.generate_cgl_knapsack_cuts;
         int freq  = p->par.cgl.generate_cgl_knapsack_cuts_freq;
         if (param < 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT && (freq < 1 ||
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_ONLY_IN_ROOT && bc_index > 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_IF_IN_ROOT && (freq < 1 ||
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_PERIODICALLY && (freq < 1 ||
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } 
         *should_generate = TRUE;
         knapsack->setMaxInKnapsack(p->par.max_cut_length); // default is 50
         knapsack->switchOffExpensive(); // gets into infinite loop if on 
         p->lp_stat.knapsack_calls++;
         break;
      }
    case CGL_GOMORY_GENERATOR:
      {
         CglGomory *gomory = (CglGomory *)generator;
         int param = p->par.cgl.generate_cgl_gomory_cuts;
         int freq  = p->par.cgl.generate_cgl_gomory_cuts_freq;
         if (param < 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT && (freq < 1 || 
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_ONLY_IN_ROOT && bc_index > 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_IF_IN_ROOT && (freq < 1 || 
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_PERIODICALLY && (freq < 1 || 
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT) {
            if (bc_index > 0) {
               if (p->comp_times.gomory_cuts > p->comp_times.lp/5) {
                  *should_generate = FALSE;
                  break;
               }
            } else {
               if (p->lp_stat.gomory_cuts > p->lp_stat.cuts_generated/2
                   && p->comp_times.gomory_cuts > 3*p->comp_times.lp) {
                  *should_generate = FALSE;
                  break;
               } else if (p->lp_stat.gomory_cuts <= 
                     p->lp_stat.cuts_generated/2 && 
                     p->comp_times.gomory_cuts > 3*p->comp_times.lp) {
                  *should_generate = FALSE;
                  break;
               }
            }
         }
         gomory->setLimit(p->par.max_cut_length);
         *should_generate = TRUE;
         p->lp_stat.gomory_calls++;
         break;
      }
    case CGL_TWOMIR_GENERATOR:
      {
         CglTwomir *twomir = (CglTwomir *)generator;
         int param = p->par.cgl.generate_cgl_twomir_cuts;
         int freq  = p->par.cgl.generate_cgl_twomir_cuts_freq;
         if (param < 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT && (freq < 1 ||
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_ONLY_IN_ROOT && bc_index > 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_IF_IN_ROOT && (freq < 1 ||
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_PERIODICALLY && (freq < 1 ||
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } 
         *should_generate = TRUE;
         twomir->setMaxElements(p->par.max_cut_length);
         twomir->setCutTypes (TRUE, TRUE, TRUE, TRUE);
         p->lp_stat.twomir_calls++;
         break;
      }
    case CGL_FLOWCOVER_GENERATOR:
      {
         CglFlowCover *flowcover = (CglFlowCover *)generator;
         int param = p->par.cgl.generate_cgl_flowcover_cuts;
         int freq  = p->par.cgl.generate_cgl_flowcover_cuts_freq;
         if (param < 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_DEFAULT && (freq < 1 ||
                  bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_ONLY_IN_ROOT && bc_index > 0) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_IF_IN_ROOT && (freq < 1 ||
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } else if (param == GENERATE_PERIODICALLY && (freq < 1 ||
               bc_index % freq != 0)) {
            *should_generate = FALSE;
            break;
         } 
         *should_generate = TRUE;
         flowcover->setNumFlowCuts(0); //needs to be called because static
         p->lp_stat.flowcover_calls++;
         break;
      }
   }
#endif
   return 0;
}

/*===========================================================================*/
#ifdef USE_CGL_CUTS
int generate_cgl_cut_of_type(lp_prob *p, int i, OsiCuts *cutlist_p, 
      int *was_tried)
{
   OsiCuts cutlist = *cutlist_p;
   int should_generate = FALSE;
   double total_time, cut_time;

   /* two times is necessary */
   cut_time     = used_time(&total_time);
   cut_time     = used_time(&total_time);

   switch (i) {
    case CGL_PROBING_GENERATOR:
      {
         CglProbing *probing = new CglProbing;
         should_use_cgl_generator(p, &should_generate, i, (void *)probing);
         if (should_generate == TRUE) {
            probing->generateCuts(*(p->lp_data->si), cutlist);
            *was_tried = TRUE;
         }
         delete probing;
         cut_time     = used_time(&total_time);
         p->comp_times.probing_cuts += cut_time;
         break;
      }
    case CGL_CLIQUE_GENERATOR:
      {
         CglClique *clique = new CglClique;
         should_use_cgl_generator(p, &should_generate, i, (void *)clique);
         if (should_generate == TRUE) {
            clique->generateCuts(*(p->lp_data->si), cutlist);
            *was_tried = TRUE;
         }
         delete clique;
         cut_time     = used_time(&total_time);
         p->comp_times.clique_cuts += cut_time;
         break;
      }
    case CGL_KNAPSACK_GENERATOR:
      {
         CglKnapsackCover *knapsack = new CglKnapsackCover;
         should_use_cgl_generator(p, &should_generate, i, (void *)knapsack);
         if (should_generate == TRUE) {
            knapsack->generateCuts(*(p->lp_data->si), cutlist);
            *was_tried = TRUE;
         }
         delete knapsack;
         cut_time     = used_time(&total_time);
         p->comp_times.knapsack_cuts += cut_time;
         break;
      }
    case CGL_GOMORY_GENERATOR:
      {
         CglGomory *gomory = new CglGomory;
         should_use_cgl_generator(p, &should_generate, i, (void *)gomory);
         if (should_generate == TRUE) {
            gomory->generateCuts(*(p->lp_data->si), cutlist);
            *was_tried = TRUE;
         }
         delete gomory;
         cut_time     = used_time(&total_time);
         p->comp_times.gomory_cuts += cut_time;
         break;
      }
    case CGL_TWOMIR_GENERATOR:
      {
         CglTwomir *twomir = new CglTwomir;
         should_use_cgl_generator(p, &should_generate, i, (void *)twomir);
         if (should_generate == TRUE) {
            twomir->generateCuts(*(p->lp_data->si), cutlist);
            *was_tried = TRUE;
         }
         delete twomir;
         cut_time     = used_time(&total_time);
         p->comp_times.twomir_cuts += cut_time;
         break;
      }
    case CGL_FLOWCOVER_GENERATOR:
      {
         CglFlowCover *flowcover = new CglFlowCover;
         should_use_cgl_generator(p, &should_generate, i, (void *)flowcover);
         if (should_generate == TRUE) {
            flowcover->generateCuts(*(p->lp_data->si), cutlist);
            *was_tried = TRUE;
         }
         delete flowcover;
         cut_time     = used_time(&total_time);
         p->comp_times.flowcover_cuts += cut_time;
         break;
      }
   }
   *cutlist_p = cutlist;
   p->comp_times.cuts += cut_time;
   return 0;
}
#endif

/*===========================================================================*/
int check_and_add_cgl_cuts(lp_prob *p, int generator, cut_data ***cuts, 
      int *num_cuts, int *bound_changes, OsiCuts *cutlist, int send_to_pool) 
{
   int          i, j, k, num_row_cuts, num_col_cuts, *is_deleted, num_elements,
                *indices, discard_cut, num_poor_quality = 0, num_unviolated = 0,
                num_duplicate = 0, *cut_size, *matind; 
   const int    max_elements = p->par.max_cut_length, 
                verbosity = p->par.verbosity;
   LPdata       *lp_data = p->lp_data;
   int          *tmp_matind = lp_data->tmp.i1;
   double       *hashes, *elements, rhs, max_coeff, min_coeff, hash_value, 
                violation, *matval, total_time, cut_time;
   double       *random_hash = lp_data->random_hash;
   const double lpetol = lp_data->lpetol;
   const double etol1000 = lpetol * 1000;
   const double *x     = lp_data->x;
   OsiRowCut    row_cut;
   OsiColCut    col_cut;
   var_desc     **vars = lp_data->vars;
   const int    is_userind_in_order = p->par.is_userind_in_order;
   cut_data     *sym_cut;

   /* two times is necessary */
   cut_time     = used_time(&total_time);
   cut_time     = used_time(&total_time);

   num_row_cuts = cutlist->sizeRowCuts();
   num_col_cuts = cutlist->sizeColCuts();
   hashes       = (double *) malloc(num_row_cuts*DSIZE);
   is_deleted   = (int *) calloc(num_row_cuts, ISIZE);
   cut_size     = (int *) calloc(num_row_cuts, ISIZE);

   j = 0;
   for (i=0; i<num_row_cuts; i++) {
      /* check for violation, duplicacy, quality of coefficients, length */
      row_cut = cutlist->rowCut(i);
      num_elements = row_cut.row().getNumElements();
      cut_size[i] = num_elements;
      indices = const_cast<int *> (row_cut.row().getIndices());
      elements = const_cast<double *> (row_cut.row().getElements());
      rhs = row_cut.rhs();
      discard_cut = FALSE;
      max_coeff = 0;
      min_coeff = DBL_MAX;

      /* length */
      if (num_elements > max_elements) {
         PRINT(verbosity,5,("Threw out cut because its length %d is too "
                  "high.\n\n\n", num_elements));
         num_poor_quality++;
         is_deleted[i] = TRUE;
         continue;
      }

      /* hash value, min, max, violation */
      hash_value = 0;
      violation = 0;
      for (int el_num=0; el_num<num_elements; el_num++) {
         if (fabs(elements[el_num])>max_coeff) {
            max_coeff = fabs(elements[el_num]);
         }
         if (fabs(elements[el_num]) < min_coeff) {
            min_coeff = fabs(elements[el_num]);
         }
         tmp_matind[el_num] = vars[indices[el_num]]->userind;
         hash_value += elements[el_num]*random_hash[tmp_matind[el_num]];
         violation += elements[el_num]*x[tmp_matind[el_num]];
      }
      hashes[i] = hash_value;
      /* see rhs as well */
      if (fabs(rhs) > lpetol) {
         if (fabs(rhs) < min_coeff) { 
            min_coeff = fabs(rhs);
         }
         if (fabs(rhs) > max_coeff) {
            max_coeff = fabs(rhs);
         }
      }
      switch (row_cut.sense()) {
       case 'L':
         violation -= rhs;
         break;
       case 'G':
         violation = rhs - violation;
         break;
       case 'E':
         violation = fabs(rhs - violation);
         break;
      }

      /* check violation */
      if (violation < lpetol) {
         PRINT(verbosity,5,("violation = %f. Threw out cut.\n", 
                  violation));
         num_unviolated++;
         is_deleted[i] = TRUE;
         continue;
      }

      /* check quality */
      if (num_elements>0) {
         if ( (max_coeff > 0 && min_coeff/max_coeff < etol1000)||
               (min_coeff > 0 && min_coeff < etol1000) ) {
            PRINT(verbosity,5,("Threw out cut because of bad coeffs.\n"));
            num_poor_quality++;
            is_deleted[i] = TRUE;
            continue;
         }
      }

      /* check for duplicates */
      if (num_elements>0) {
         for (k=i-1; k>-1; k--) {
            if (is_deleted[k] == TRUE) {
               continue;
            }
            if (cut_size[k] != num_elements || 
                  fabs(hashes[k]-hash_value) > lpetol) {
               continue;
            } else {
               break;
            }
         }
         if (k>-1) {
            PRINT(verbosity,5,("cut #%d is same as cut #%d\n", i, k));
            num_duplicate++;
            is_deleted[i] = TRUE;
            continue;
         }
      }

      /* check if sense is 'R' */
      if (row_cut.sense()=='R') {
         PRINT(verbosity,5,("cut #%d has a range. thrown out.\n", i));
         is_deleted[i] = TRUE;
         continue;
      }

      /* cut is accepted. congratulations. */
      j++;
   }

   /* copy the accepted cuts */
   if (*cuts){
      *cuts = (cut_data **)realloc(*cuts, (*num_cuts+j)*sizeof(cut_data *));
   }else{
      *cuts = (cut_data **)malloc(j*sizeof(cut_data *));
   }
   k = *num_cuts;
   for (i=0; i<num_row_cuts; i++) {
      if (is_deleted[i] == TRUE) {
         continue;
      }
      row_cut = cutlist->rowCut(i);
      rhs = row_cut.rhs();
      num_elements = row_cut.row().getNumElements();
      //PRINT(verbosity, -1,("length = %d \n", num_elements));
      indices = const_cast<int *> (row_cut.row().getIndices());
      elements = const_cast<double *> (row_cut.row().getElements());
      (*cuts)[k] =  (cut_data *) calloc(1, sizeof(cut_data));
      sym_cut    = (*cuts)[k];
      sym_cut->type = EXPLICIT_ROW;
      sym_cut->rhs = rhs;
      sym_cut->range = row_cut.range();
      sym_cut->size = (num_elements * (ISIZE + DSIZE) + DSIZE);
      sym_cut->coef = (char *) malloc (sym_cut->size);
      sym_cut->sense = row_cut.sense();
      ((double *) (sym_cut->coef))[0] = 0; // otherwise valgrind complains.
      ((int *) (sym_cut->coef))[0] = num_elements;

      //Here, we have to pad the initial int to avoid misalignment, so we
      //add DSIZE bytes to get to a double boundary
      matval = (double *) (sym_cut->coef + DSIZE);
      matind = (int *) (sym_cut->coef + (num_elements + 1)*DSIZE);
      memcpy((char *)matval, (char *)elements, num_elements * DSIZE);
      if (is_userind_in_order == TRUE) {
         memcpy((char*)matind, (char *)indices, num_elements * ISIZE);
      } else {
         for (int i2=0; i2<num_elements; i2++) {
            tmp_matind[i2] = vars[indices[i2]]->userind;
         }
         memcpy((char*)matind, (char *)tmp_matind, num_elements * ISIZE);
      }

      qsort_id(matind, matval, num_elements);

      sym_cut->branch = DO_NOT_BRANCH_ON_THIS_ROW;

      sym_cut->deletable = TRUE;
      if (send_to_pool){
         sym_cut->name = CUT__SEND_TO_CP;
      }else{
         sym_cut->name = CUT__DO_NOT_SEND_TO_CP;
      }	    
      k++;
   }
   *num_cuts = k;
   // TODO: short circuit the copying to row data and si */
   for (i=0; i<num_row_cuts; i++) {
      cutlist->eraseRowCut(0);
   }

   for (i=0; i<num_col_cuts; i++) {
      col_cut = cutlist->colCut(i);
      if (verbosity>10) {
         col_cut.print();
      }
      indices  = const_cast<int *>(col_cut.lbs().getIndices());
      elements = const_cast<double *>(col_cut.lbs().getElements());
      for (j=0;j<col_cut.lbs().getNumElements();j++) {
         if (vars[indices[j]]->new_lb < elements[j]) {
            vars[indices[j]]->new_lb = elements[j];
            change_lbub(lp_data, indices[j], elements[j], 
                  vars[indices[j]]->new_ub);
            (*bound_changes)++;
         }
      }
      indices  = const_cast<int *>(col_cut.ubs().getIndices());
      elements = const_cast<double *>(col_cut.ubs().getElements());
      for (j=0;j<col_cut.ubs().getNumElements();j++) {
         if (vars[indices[j]]->new_ub > elements[j]) {
            vars[indices[j]]->new_ub = elements[j];
            change_lbub(lp_data, indices[j], vars[indices[j]]->new_lb,
                  elements[j]);
            (*bound_changes)++;
         }
      }
   }


   
   for (i=0; i<num_col_cuts; i++) {
      cutlist->eraseColCut(0);
   }

   FREE(hashes);
   FREE(is_deleted);
   FREE(cut_size);

   /* update statistics */
   p->lp_stat.num_duplicate_cuts += num_duplicate;
   p->lp_stat.num_poor_cuts += num_poor_quality;
   p->lp_stat.num_unviolated_cuts += num_unviolated;
   p->lp_stat.cuts_generated += num_row_cuts;
   if (p->bc_level<1) {
      p->lp_stat.cuts_root   += num_row_cuts;
   }

   switch (generator) {
    case (CGL_PROBING_GENERATOR):
      p->lp_stat.probing_cuts += num_row_cuts;
      if (p->bc_level<1) {
         p->lp_stat.probing_cuts_root += num_row_cuts;
      }
      break;
    case (CGL_CLIQUE_GENERATOR):
      p->lp_stat.clique_cuts += num_row_cuts;
      if (p->bc_level<1) {
         p->lp_stat.clique_cuts_root += num_row_cuts;
      }
      break;
    case (CGL_KNAPSACK_GENERATOR):
      p->lp_stat.knapsack_cuts += num_row_cuts;
      if (p->bc_level<1) {
         p->lp_stat.knapsack_cuts_root += num_row_cuts;
      }
      break;
    case (CGL_GOMORY_GENERATOR):
      p->lp_stat.gomory_cuts += num_row_cuts;
      if (p->bc_level<1) {
         p->lp_stat.gomory_cuts_root += num_row_cuts;
      }
      break;
    case (CGL_TWOMIR_GENERATOR):
      p->lp_stat.twomir_cuts += num_row_cuts;
      if (p->bc_level<1) {
         p->lp_stat.twomir_cuts_root += num_row_cuts;
      }
      break;
    case (CGL_FLOWCOVER_GENERATOR):
      p->lp_stat.flowcover_cuts += num_row_cuts;
      if (p->bc_level<1) {
         p->lp_stat.flowcover_cuts_root += num_row_cuts;
      }
      break;
   }

   cut_time = used_time(&total_time);
   p->comp_times.dupes_and_bad_coeffs_in_cuts += cut_time;
        
   return 0;
}

/*===========================================================================*/
int should_stop_adding_cgl_cuts(lp_prob *p, int i, int *should_stop)
{
   *should_stop = 0;
   return 0;
}

/*===========================================================================*/
int update_pcost(lp_prob *p)
{
#ifdef COMPILE_IN_LP
   bc_node *parent = p->tm->active_nodes[p->proc_index]->parent;
   char sense = parent->bobj.sense[0];
   int branch_var = parent->bobj.position;
   double *pcost_down = p->pcost_down;
   double *pcost_up = p->pcost_up;
   int *br_rel_down = p->br_rel_down;
   int *br_rel_up = p->br_rel_up;
   double objval = p->lp_data->objval;
   double oldobjval = p->tm->active_nodes[p->proc_index]->lower_bound;
   double oldx =  parent->bobj.value;
   double *x;
   get_x(p->lp_data);
   x = p->lp_data->x;
   if (parent->children[0]->bc_index != p->bc_index) {
      sense = (sense == 'L') ? 'G' : 'L';
   }
   if (sense == 'L') {
      if (oldx - x[branch_var] > 1e-5) {
         pcost_down[branch_var] = (pcost_down[branch_var]*
               br_rel_down[branch_var] + (objval - oldobjval)/
               (oldx-x[branch_var]))/(br_rel_down[branch_var] + 1);
         //printf("new pcost_down[%d] = %f\n", branch_var, pcost_down[branch_var]);
         br_rel_down[branch_var]++;
      } else {
         PRINT(p->par.verbosity, 0, ("warning: poor lpetol used while branching\n"));
      }
   } else {
      if (x[branch_var] - oldx > 1e-5) {
         pcost_up[branch_var] = (pcost_up[branch_var]*
               br_rel_up[branch_var] + (objval - oldobjval)/
               (x[branch_var]-oldx))/(br_rel_up[branch_var] + 1);
         //printf("new pcost_up[%d] = %f\n", branch_var, pcost_up[branch_var]);
         br_rel_up[branch_var]++;
      } else {
         PRINT(p->par.verbosity, 0, ("warning: poor lpetol used while branching\n"));
      }
   }

   p->lp_stat.avg_br_obj_impr_in_path = ((p->bc_level-1)*
         p->lp_stat.avg_br_obj_impr_in_path + objval - oldobjval)/p->bc_level;
#endif
   return 0;
}
/*===========================================================================*/

/* check if lb <= ub for each variable. otherwise fathom this branch. */
int check_bounds(lp_prob *p, int *termcode)
{
   int i;
   double *lb, *ub;
   const double lpetol = p->lp_data->lpetol;
   const int n = p->lp_data->n;
   LPdata *lp_data = p->lp_data;
   
   get_bounds(lp_data);
   lb = lp_data->lb;
   ub = lp_data->ub;

   for (i=0; i<n; i++) {
      if (lb[i] > ub[i]+lpetol) {
         break;
      }
   }
   if (i<n) {
      *termcode = LP_D_UNBOUNDED;
   }
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
