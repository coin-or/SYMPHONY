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

#define COMPILE_FOR_LP

#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "proccomm.h"
#include "qsortucb.h"
#include "lp.h"
#include "messages.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "BB_types.h"
#include "pack_cut.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains general LP functions.
\*===========================================================================*/

/*===========================================================================*\
 * This small function simply returns a pointer to the current LP 
 * problem structure so that it can be accessed from anywhere within the LP  
 * without explicitly passing a pointer. This is essentially only uselful    
 * to the more advanced user who wishes to have access to the internal data  
 * structures                                                                
\*===========================================================================*/

lp_prob *get_lp_ptr(lp_prob **lp_list)
{
#ifdef _OPENMP
   int thread_num = omp_get_thread_num();
#else
   int thread_num = 0;
#endif
   static lp_prob **lp;

   if (lp_list){
      lp = lp_list;
   }

   return(lp[thread_num]);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function receives the problem data (if we are running in parallel)   
 * and intitializes the data structures.                                     
\*===========================================================================*/

void lp_initialize(lp_prob *p, int master_tid)
{
#ifndef COMPILE_IN_LP
   int msgtag, bytes, r_bufid;
#endif
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP)
   int s_bufid;
#endif
   int i;
   constraint *rows;
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
#pragma omp critical (lp_solver)
/*__BEGIN_EXPERIMENTAL_SECTION__*/
   {
#if 0
      system("/bin/rm -f /usr/lib/.cpxlicense/.tko*");
#endif
      open_lp_solver(p->lp_data);
   }
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
   open_lp_solver(p->lp_data);
#endif

   (void) used_time(&p->tt);

#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP)
   s_bufid = init_send(DataInPlace);
   send_msg(p->master, REQUEST_FOR_LP_DATA);
   freebuf(s_bufid);

   receive_lp_data_u(p);
#endif
   
   if (p->par.tailoff_gap_backsteps > 0 ||
       p->par.tailoff_obj_backsteps > 1){
      i = MAX(p->par.tailoff_gap_backsteps, p->par.tailoff_obj_backsteps);
      p->obj_history = (double *) malloc((i + 1) * DSIZE);
   }
#ifndef COMPILE_IN_LP
   if (p->par.use_cg){
      r_bufid = receive_msg(p->tree_manager, LP__CG_TID_INFO);
      receive_int_array(&p->cut_gen, 1);
      freebuf(r_bufid);
   }
#endif
   p->lp_data->rows =
      (constraint *) malloc((p->base.cutnum + BB_BUNCH) * sizeof(constraint));
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

#ifdef COMPILE_IN_CG
   if (!p->cgp){
      p->cgp = (cg_prob *) calloc(1, sizeof(cg_prob));
      get_cg_ptr(&(p->cgp));
   }
   
   cg_initialize(p->cgp, p->master);
   
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   if (p->cgp->dcmp_data.lp_data)
      p->cgp->dcmp_data.lp_data->cpxenv = p->lp_data->cpxenv;
#endif
/*___END_EXPERIMENTAL_SECTION___*/
#endif
}   

/*===========================================================================*/

/*===========================================================================*\
 * This function continues to dive down the current chain until told to stop
 * by the tree manager. 
\*===========================================================================*/

int process_chain(lp_prob *p)
{
   /* Create the LP */
   if (! create_lp_u(p)){
      /* User had problems creating initial LP. Abandon node. */
      return(FALSE);
   }

   p->last_gap = 0.0;
   if (p->has_ub && p->par.set_obj_upper_lim)
      set_obj_upper_lim(p->lp_data, p->ub - p->par.granularity);
   
   if (p->colgen_strategy & COLGEN_REPRICING){
      if (p->par.verbosity > 0){
	 printf("****************************************************\n");
	 printf("* Now repricing NODE %i LEVEL %i\n",
		p->bc_index, p->bc_level);
	 printf("****************************************************\n\n");
      }
      repricing(p);
      free_node_dependent(p);
   }else{
      if (p->par.verbosity > 0){
	 printf("****************************************************\n");
	 printf("* Now processing NODE %i LEVEL %i (from TM)\n",
		p->bc_index, p->bc_level);
	 printf("****************************************************\n\n");
	 PRINT(p->par.verbosity, 4, ("Diving set to %i\n\n", p->dive));
      }
      fathom_branch(p);

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

   return(TRUE);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function receives information for an active node, processes that     *
 * node and then decides which one of the children of that node should be    *
 * processed next. It then recursively processes the child until the branch  *
 * is pruned at point                                                        *
\*===========================================================================*/

void fathom_branch(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   node_times *comp_times = &p->comp_times;
   char first_in_loop = TRUE;
   int iterd, termcode, i;
   int cuts, no_more_cuts_count;
   int num_errors = 0;

   check_ub(p);
   p->iter_num = p->node_iter_num = 0;
  
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

      PRINT(p->par.verbosity, 2,
	    ("\n\n**** Starting iteration %i ****\n\n", p->iter_num));

      /* the following call fills up:
       *    lp_data->basis, lp_data->xbzero, lp_data->objval
       */
      termcode = dual_simplex(lp_data, &iterd);
      /* Get all the data we will need */
      get_dj_pi(lp_data);
      get_slacks(lp_data);

      /* display the current solution */
      PRINT(p->par.verbosity, 2, ("The LP value is: %.3f [%i,%i]\n\n",
				  lp_data->objval, termcode, iterd));

      switch (termcode){
       case D_ITLIM:      /* impossible, since itlim is set to infinity */
       case D_INFEASIBLE: /* this is impossible (?) as of now */
       case ABANDONED:
	 printf("######## Unexpected termcode: %i ########\n\n", termcode);
	 if (p->par.try_to_recover_from_error && (++num_errors == 1)){
	    /* Try to resolve it from scratch */
	    continue;
	 }
	 MakeMPS(lp_data, p->bc_index, p->iter_num);
	 exit(-3);

       case D_UNBOUNDED: /* the primal problem is infeasible */
       case D_OBJLIM:
       case OPTIMAL:
	 if (termcode == D_UNBOUNDED){
	    PRINT(p->par.verbosity, 1, ("Feasibility lost -- "));
	 }else if ((p->has_ub && lp_data->objval > p->ub - p->par.granularity)
		   || termcode == D_OBJLIM){
	    PRINT(p->par.verbosity, 1, ("Terminating due to high cost -- "));
	 }else{ /* optimal and not too high cost */
	    break;
	 }
	 comp_times->lp += used_time(&p->tt);
	 if (fathom(p, (termcode != D_UNBOUNDED))){
	    comp_times->communication += used_time(&p->tt);
	    return;
	 }else{
	    first_in_loop = FALSE;
	    comp_times->communication += used_time(&p->tt);
	    continue;
	 }
      }

      /* If come to here, the termcode must have been OPTIMAL and the
       * cost cannot be too high. */
      /* is_feasible_u() fills up lp_data->x, too!! */
      if (is_feasible_u(p) == FEASIBLE){
	 if (p->par.verbosity > 2){
	    printf ("Now displaying the feasible solution ...\n");
	    display_lp_solution_u(p, DISP_FEAS_SOLUTION);
	 }
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

	 fix_variables(p);

	 comp_times->fixing += used_time(&p->tt);

	 if (!first_in_loop){
	    cuts = check_row_effectiveness(p);
	 }

	 /*------------------------------------------------------------------*\
	  * receive the cuts from the cut generator and the cut pool
	 \*------------------------------------------------------------------*/

	 cuts += receive_cuts(p, first_in_loop, no_more_cuts_count);
      }

      comp_times->lp += used_time(&p->tt);
      if (cuts < 0){ /* i.e. feasible solution is found */
	 if (fathom(p, TRUE)){
	    return;
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
	    p->node_iter_num = 0;
	 }
	 break;
#endif
       case FATHOMED_NODE:
	 comp_times->strong_branching += used_time(&p->tt);
	 return;

       default: /* the return value is the number of cuts added */
	 if (p->par.verbosity > 2){
	    printf("Continue with this node.");
	    if (cuts > 0)
	       printf(" %i cuts added alltogether in iteration %i",
		      cuts, p->iter_num);
	    printf("\n\n");
	 }
	 break;
      }
      comp_times->strong_branching += used_time(&p->tt);

      check_ub(p);
      first_in_loop = FALSE;
   }
   comp_times->lp += used_time(&p->tt);
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
      send_node_desc(p, primal_feasible ? (termcode == OPT_FEASIBLE ?
					   FEASIBLE_PRUNED: OVER_UB_PRUNED) :
		     INFEASIBLE_PRUNED);
      return(TRUE);
   }

   if (p->colgen_strategy & COLGEN_REPRICING)
      colgen = FATHOM__GENERATE_COLS__RESOLVE;

   switch (colgen){
    case FATHOM__DO_NOT_GENERATE_COLS__DISCARD:
      PRINT(p->par.verbosity, 1, ("Pruning node\n\n"));
      send_node_desc(p, termcode == OPT_FEASIBLE ? FEASIBLE_PRUNED :
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
	  termcode == D_OBJLIM || termcode == OPT_FEASIBLE){
	 /* fathomable */
	 if (termcode == D_OBJLIM ||
	     (p->has_ub && lp_data->objval > p->ub - p->par.granularity)){
	    PRINT(p->par.verbosity, 1,
		  ("Fathoming node (discovered tdf & high cost)\n\n"));
	 }else{
	    PRINT(p->par.verbosity, 1,
		  ("Fathoming node (discovered tdf & feasible)\n\n"));
	 }
	 send_node_desc(p, termcode == OPT_FEASIBLE ? FEASIBLE_PRUNED :
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
	 free_col_set(&new_cols); return(TRUE);
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

void repricing(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   node_times *comp_times = &p->comp_times;
   int iterd, termcode;
   int num_errors = 0;
   our_col_set *new_cols = NULL;
   int dual_feas, new_vars, cuts, no_more_cuts_count;

   check_ub(p);
   p->iter_num = 0;
  
   /*------------------------------------------------------------------------*\
    * The main loop -- continue solving relaxations until TDF
   \*------------------------------------------------------------------------*/

   while (TRUE){
      p->iter_num++;

      PRINT(p->par.verbosity, 2,
	    ("\n\n**** Starting iteration %i ****\n\n", p->iter_num));

      /* the following call fills up:
       *    lp_data->bhead, lp_data->xbzero, lp_data->objval
       */
      termcode = dual_simplex(lp_data, &iterd);
      /* Get all the data we will need */
      get_dj_pi(lp_data);
      get_slacks(lp_data);

      /* display the current solution */
      PRINT(p->par.verbosity, 2, ("The LP value is: %.3f [%i,%i]\n\n",
				  lp_data->objval, termcode, iterd));

      comp_times->lp += used_time(&p->tt);

      switch (termcode){
       case D_ITLIM:      /* impossible, since itlim is set to infinity */
       case D_INFEASIBLE: /* this is impossible (?) as of now */
       case ABANDONED:
	 printf("######## Unexpected termcode: %i ########\n\n", termcode);
	 if (p->par.try_to_recover_from_error && (++num_errors == 1)){
	    /* Try to resolve it from scratch */
	    continue;
	 }
	 exit(-2);

       case D_UNBOUNDED: /* the primal problem is infeasible */
       case D_OBJLIM:
       case OPTIMAL:
	 if (termcode == D_UNBOUNDED){
	    PRINT(p->par.verbosity, 1, ("Feasibility lost -- "));
	 }else if ((p->has_ub && lp_data->objval > p->ub - p->par.granularity)
		   || termcode == D_OBJLIM){
	    PRINT(p->par.verbosity, 1, ("Terminating due to high cost -- "));
	 }else{ /* optimal and not too high cost */
	    break;
	 }
	 comp_times->lp += used_time(&p->tt);
	 if (fathom(p, (termcode != D_UNBOUNDED))){
	    comp_times->communication += used_time(&p->tt);
	    return;
	 }else{
	    comp_times->communication += used_time(&p->tt);
	    continue;
	 }
      }

      /* If come to here, the termcode must have been OPTIMAL and the
       * cost cannot be too high. */
      /* is_feasible_u() fills up lp_data->x, too!! */
      if (is_feasible_u(p) == FEASIBLE){
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

	 fix_variables(p);

	 comp_times->fixing += used_time(&p->tt);

	 cuts = 0;
	 if (p->cut_gen || p->cut_pool){
	    cuts = check_row_effectiveness(p);
	 }

	 /*------------------------------------------------------------------*\
	  * receive the cuts from the cut generator and the cut pool
	 \*------------------------------------------------------------------*/

	 cuts += receive_cuts(p, TRUE, no_more_cuts_count);
      }

      comp_times->lp += used_time(&p->tt);
      if (cuts < 0){ /* i.e. feasible solution is found */
	 if (fathom(p, TRUE)){
	    comp_times->communication += used_time(&p->tt);
	    return;
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
}

/*===========================================================================*/

void resize_lp_arrays(LPdata *lp_data, char do_realloc, char set_max,
		      int row_num, int col_num, int nzcnt)
{
   char resize_m = FALSE;
   char resize_n = FALSE;
   int maxm, maxn, maxnz, maxmax;

   if (set_max){
      maxm = row_num;
      maxn = col_num;
      maxnz = nzcnt;
   }else{
      maxm = lp_data->m + row_num;
      maxn = lp_data->n + col_num;
      maxnz = lp_data->nz + nzcnt;
   }

   if (maxm > lp_data->maxm){
      resize_m = TRUE;
      lp_data->maxm = maxm + (set_max ? 0 : BB_BUNCH);
      if (! do_realloc){
         FREE(lp_data->dualsol);
         FREE(lp_data->bhead);
         FREE(lp_data->xbzero);
         lp_data->dualsol = (double *) malloc(lp_data->maxm * DSIZE);
         lp_data->xbzero = (double *) malloc(lp_data->maxm * DSIZE);
         lp_data->bhead = (int *) malloc(lp_data->maxm * DSIZE);
      }else{
         lp_data->dualsol = (double *) realloc((char *)lp_data->dualsol,
                                               lp_data->maxm * DSIZE);
         lp_data->xbzero = (double *) realloc((char *)lp_data->xbzero,
                                              lp_data->maxm * DSIZE);
         lp_data->bhead = (int *) realloc((char *)lp_data->bhead,
                                          lp_data->maxm * DSIZE);
      }
      /* rows is realloc'd in either case just to keep the base constr */
      lp_data->rows = (constraint *) realloc((char *)lp_data->rows,
                                             lp_data->maxm*sizeof(constraint));
   }
   if (maxn > lp_data->maxn){
      int i, oldmaxn = MAX(lp_data->maxn, lp_data->n);
      var_desc **vars;
      resize_n = TRUE;
      lp_data->maxn = maxn + (set_max ? 0 : 5 * BB_BUNCH);
      if (! do_realloc){
         FREE(lp_data->x);
         FREE(lp_data->dj);
         FREE(lp_data->status);
         lp_data->x = (double *) malloc(lp_data->maxn * DSIZE);
         lp_data->dj = (double *) malloc(lp_data->maxn * DSIZE);
         lp_data->status = (char *) malloc(lp_data->maxn * CSIZE);
      }else{
         lp_data->x = (double *) realloc((char *)lp_data->x,
                                         lp_data->maxn * DSIZE);
         lp_data->dj = (double *) realloc((char *)lp_data->dj,
                                          lp_data->maxn * DSIZE);
         lp_data->status = (char *) realloc((char *)lp_data->status,
                                            lp_data->maxn * CSIZE);
      }
      /* vers is realloc'd in either case just to keep the base vars */
      lp_data->vars = (var_desc **) realloc((char *)lp_data->vars,
                                            lp_data->maxn*sizeof(var_desc *));
      vars = lp_data->vars + oldmaxn;
      for (i = lp_data->maxn - oldmaxn - 1; i >= 0; i--)
	 vars[i] = (var_desc *) malloc( sizeof(var_desc) );
   }
   if (maxnz > lp_data->maxnz){
      lp_data->maxnz = maxnz + (set_max ? 0 : 20 * BB_BUNCH);
   }

   /* re(m)alloc the tmp arrays */
   if (resize_m || resize_n){
      temporary *tmp = &lp_data->tmp;
      maxm = lp_data->maxm;
      maxn = lp_data->maxn;
      maxmax = MAX(maxm, maxn);
      /* anything with maxm and maxn in it has to be resized */
      FREE(tmp->c);
      FREE(tmp->i1);
      FREE(tmp->d);
      tmp->c = (char *) malloc(CSIZE * maxmax);
      tmp->i1 = (int *) malloc(ISIZE * MAX(3*maxm, 2*maxn + 1));
      tmp->d = (double *) malloc(DSIZE * 2 * maxmax);
      /* These have to be resized only if maxm changes */
      if (resize_m){
	 FREE(tmp->i2);
	 FREE(tmp->p1);
	 FREE(tmp->p2);
	 tmp->i2 = (int *) malloc(maxm * ISIZE);
	 tmp->p1 = (void *) malloc(maxm * sizeof(void *));
	 tmp->p2 = (void *) malloc(maxm * sizeof(void *));
      }
   }

   resize_lp_solver_arrays(lp_data);
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

/*__BEGIN_EXPERIMENTAL_SECTION__*/
int collect_nonzeros(lp_prob *p, double *x, int *tind, double *tx,
		     char *status)
{
   var_desc **vars = p->lp_data->vars;
   int n = p->lp_data->n;
   int i, cnt = 0;
   double lpetol = p->lp_data->lpetol;

   colind_sort_extra(p);
   if (status){
      int *temp_tind;
      
      for (i=0; i<n; i++){
	 if (x[i] > lpetol){
	    tind[cnt] = vars[i]->userind;
	    tx[cnt] = x[i];
	    status[cnt++] = p->lp_data->status[i];
	 }
      }
      temp_tind = (int *) malloc(cnt * sizeof(int));
      memcpy(temp_tind, tind, cnt * sizeof(int));
      qsortucb_ic(temp_tind, status, cnt);
      FREE(temp_tind);
   }else{
      for (i=0; i<n; i++){
	 if (x[i] > lpetol){
	    tind[cnt] = vars[i]->userind;
	    tx[cnt++] = x[i];
	 }
      }
   }
   /* order indices and values according to indices */
   qsortucb_id(tind, tx, cnt);
   return(cnt);
}
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
int collect_nonzeros(lp_prob *p, double *x, int *tind, double *tx)
{
   var_desc **vars = p->lp_data->vars;
   int n = p->lp_data->n;
   int i, cnt = 0;
   double lpetol = p->lp_data->lpetol;

   colind_sort_extra(p);
   for (i=0; i<n; i++){
      if (x[i] > lpetol){
	 tind[cnt] = vars[i]->userind;
	 tx[cnt++] = x[i];
      }
   }
   /* order indices and values according to indices */
   qsortucb_id(tind, tx, cnt);
   return(cnt);
}
#endif

/*===========================================================================*/

int collect_fractions(lp_prob *p, double *x, int *tind, double *tx)
{
   var_desc **vars = p->lp_data->vars;
   int n = p->lp_data->n;
   int i, cnt = 0;
   double lpetol = p->lp_data->lpetol, xi;

   colind_sort_extra(p);
   for (i=0; i<n; i++){
      xi = x[i];
      if (xi - floor(xi) > lpetol && ceil(xi) - xi > lpetol){
	 tind[cnt] = vars[i]->userind;
	 tx[cnt++] = x[i];
      }
   }
   /* order indices and values according to indices */
   qsortucb_id(tind, tx, cnt);
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
   constraint *rows = lp_data->rows;
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
      constraint *tmp_rows=(constraint *) malloc(cutcnt*sizeof(constraint));
      
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
	 qsortucb_ii(ulist, ecstat, extravarnum);
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
	 desc->not_fixed.list = malloc(desc->not_fixed.size * ISIZE);
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
	    if (pstat(p->tree_manager) != OK){
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
      qsortucb_ii(clist, erstat, cutindsize);
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

   if (gap_backsteps < 1 && obj_backsteps < 2)
      /* The user is stupid... asks for tailoff (since we came to this
	 function) yet doesn't want to check any kind of tailoff (since this
	 condition is true). Report no tailoff. */
      return(FALSE);

   /* shift the data in obj_hist by one to the right and insert the
      most recent objval to be the 0th */
   for (i = MIN(p->iter_num-1, maxsteps) - 1; i >= 0; i--)
      obj_hist[i+1] = obj_hist[i];
   obj_hist[0] = p->lp_data->objval;

   /* If the history is not long enough then just return */
   if (p->iter_num <= maxsteps)
      return(FALSE);

   /* if there is an upper bound and we want gap based tailoff:
      tailoff_gap is false if the average of the consecutive gap ratios is
      less than gap_frac */
   if (p->has_ub && gap_backsteps > 0) {
      ub = p->ub;
      for (i = 1, sum = 0; i <= gap_backsteps; i++)
	 sum += (ub - obj_hist[i-1]) / (ub - obj_hist[i]);
      if (sum / gap_backsteps < p->par.tailoff_gap_frac)
	 return(FALSE); /* no tailoff */
   }

   /* if we want objective value based tailoff:
      tailoff_obj is false if the average of the objective difference ratios
      is greater than obj_frac */
   if (obj_backsteps > 1){
      for (i = 2, sum = 0; i <= obj_backsteps; i++)
	 if (obj_hist[i-1] - obj_hist[i] > 0)
	    sum += (obj_hist[i-2]-obj_hist[i-1]) / (obj_hist[i-1]-obj_hist[i]);
	 else
	    sum += 1;
      if (sum / (obj_backsteps - 1) > p->par.tailoff_obj_frac)
	 return(FALSE); /* no tailoff */
   }

   PRINT(p->par.verbosity, 3, ("Branching because of tailoff!\n"));

   return(TRUE); /* gone thru everything ==> tailoff */
}

/*===========================================================================*/

void MakeMPS(LPdata *lp_data, int bc_index, int iter_num)
{
   FILE *f;
   int i, j;
   char name1[50] = "";

   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   sprintf(name1,"/home/tkr/tmp/matrices/test.%i.%i.mps", bc_index, iter_num);
   /*___END_EXPERIMENTAL_SECTION___*/
   /*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
   sprintf(name1,"matrix.%i.%i.mps", bc_index, iter_num);
#endif
   
   f = fopen(name1, "w");

   fprintf(f, "NAME    %s\nROWS\n N      obj\n", name1);

   for (i = 0; i < lp_data->m; i++)
      fprintf(f, " %c  %6ir\n", lp_data->sense[i], i);

   fprintf(f, "COLUMNS\n");

   for (j = 0; j < lp_data->n; j++){
      fprintf(f, " %6ix      obj  %f\n", j, lp_data->obj[j]);
      for (i = 0; i < lp_data->matcnt[j]; i++){
	 fprintf(f, " %6ix  %6ir  %f\n",
		 j, lp_data->matind[lp_data->matbeg[j]+i],
		 lp_data->matval[lp_data->matbeg[j]+i]);
      }
   }

   fprintf(f, "RHS\n");

   for (i = 0; i < lp_data->m; i++)
      fprintf(f, "     rhs  %6ir  %f\n", i, lp_data->rhs[i]);
   
   fprintf(f, "BOUNDS\n");

   for (j = 0; j < lp_data->n; j++){
      fprintf(f, " LO  BOUND  %6ix      %f\n", j, lp_data->lb[j]);
      fprintf(f, " UP  BOUND  %6ix      %f\n", j, lp_data->ub[j]);
   }

   fprintf(f, "ENDATA\n");

   fclose(f);
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
}
#endif
#ifdef COMPILE_IN_CG
   cg_close(p->cgp);
#endif
#ifndef COMPILE_IN_TM
   free_lp(p);
#endif
}

