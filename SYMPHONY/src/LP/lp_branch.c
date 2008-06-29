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

#include <math.h>
#include <memory.h>
#include <stdlib.h>

#include "sym_lp.h"
#include "sym_qsort.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_types.h"
#include "sym_lp_solver.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP functions related to branching.
\*===========================================================================*/

void add_slacks_to_matrix(lp_prob *p, int cand_num, branch_obj **candidates)
{
   LPdata *lp_data = p->lp_data;
   int *index;
   int m = p->lp_data->m;
   int j, k;
   branch_obj *can;
   row_data *newrows;
   waiting_row **wrows;

   for (j=cand_num-1; j >= 0; j--)
      if (candidates[j]->type == CANDIDATE_CUT_NOT_IN_MATRIX)
	 break;

   if (j < 0) /* there is nothing to add */
      return;

   /* We'll create a waiting row for each cut, add them to the matrix,
      then set their status to be free */
   wrows = (waiting_row **) malloc(cand_num * sizeof(waiting_row *));
   /* can't use tmp.p, because that might get resized in add_row_set */
   for (k=0; j >= 0; j--){
      can = candidates[j];
      if (can->type == CANDIDATE_CUT_NOT_IN_MATRIX){
	 wrows[k] = can->row;
	 can->row = NULL;
	 can->position = m + k;
	 can->type = CANDIDATE_CUT_IN_MATRIX;
	 k++;
      }
   }
   add_row_set(p, wrows, k);
   /* To satisfy the size requirements in free_row_set, the following sizes
    * are needed: tmp.c:2*m   tmp.i1:3*m   tmp.d:m */
   FREE(wrows);
   index = lp_data->tmp.i1;
   for (j = 0; j < k; j++)
      index[j] = m + j;
   free_row_set(lp_data, k, index);
   newrows = lp_data->rows + m; /* m is still the old one! */
   for (j = 0; j < k; j++){
      newrows[j].ineff_cnt = (MAXINT) >> 1; /* it is slack... */
      newrows[j].free = TRUE;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Ok. So there were violated cuts, either in waiting_rows or among the
 * slacks (or both). We just have to add those among the slacks to the
 * waiting_rows; add the appropriate number of cuts to the problem and
 * continue the node (this second part is just like the end of receive_cuts().
\*===========================================================================*/
   
int add_violated_slacks(lp_prob *p, int cand_num, branch_obj **candidates)
{
   LPdata *lp_data = p->lp_data;
   waiting_row **new_rows;
   int i, new_row_num = 0;

   /* If there are any violated (former) slack, unpack them and add them
    * to the set of waiting rows. */
   if (cand_num > 0){
      new_rows = (waiting_row **) lp_data->tmp.p1; /* m (actually, candnum<m */
      for (i=0; i<cand_num; i++){
	 if (candidates[i]->type == VIOLATED_SLACK){
	    new_rows[new_row_num++] = candidates[i]->row;
	    candidates[i]->row = NULL;
	 }
      }
      if (new_row_num > 0)
	 add_new_rows_to_waiting_rows(p, new_rows, new_row_num);
   }

   return( p->waiting_row_num == 0 ? 0 : add_best_waiting_rows(p) );
}

/*===========================================================================*/

int select_branching_object(lp_prob *p, int *cuts, branch_obj **candidate)
{

   LPdata *lp_data = p->lp_data;
   var_desc **vars;
   row_data *rows;
   int m;
#ifndef MAX_CHILDREN_NUM
   int maxnum;
   double *objval, *pobj;
   int *termcode, *pterm, *feasible, *pfeas, *iterd, *piter;
#ifdef COMPILE_FRAC_BRANCHING
   int *frnum, *pfrnum, **frind, **pfrind;
   double **frval, **pfrval;
#endif
#endif
   int i, j, k, l, branch_var, branch_row, min_ind;
   double lb, ub, oldobjval, min_obj;
   cut_data *cut;
   branch_obj *can, *best_can = NULL;
#ifdef COMPILE_FRAC_BRANCHING
   int *xind;
   double *xval;
#endif

   /* These are the return values from select_candidates_u() */
   int cand_num = 0, new_vars = 0, *indices;
   double *values;
   branch_obj **candidates = NULL;
#ifdef STATISTICS
   int itlim = 0, cnum = 0;
#endif
   int should_use_hot_starts;
   double total_time = 0;
   double st_time = 0;
   int total_iters, should_continue;

   /*------------------------------------------------------------------------*\
    * First we call select_candidates_u() to select candidates. It can
    * -- return with DO_BRANCH and a bunch of candidates, or
    * -- return with DO_NOT_BRANCH along with a bunch of violated cuts
    *    in the matrix and/or among the slack_cuts, or
    * -- return with DO_NOT_BRANCH__FATHOMED, i.e., the node can be fathomed.
   \*------------------------------------------------------------------------*/

   j = select_candidates_u(p, cuts, &new_vars, &cand_num, &candidates);
   switch (j){
    case DO_NOT_BRANCH__FATHOMED:
      *candidate = NULL;
      return(DO_NOT_BRANCH__FATHOMED);
    case DO_NOT_BRANCH__FEAS_SOL:
      *candidate = NULL;
      return(DO_NOT_BRANCH__FEAS_SOL);
    case DO_NOT_BRANCH:
      if (cand_num)
	 *cuts += add_violated_slacks(p, cand_num, candidates);
#ifdef DO_TESTS
      if (*cuts == 0 && new_vars == 0){
	 printf("Error! Told not to branch, but there are no new cuts or ");
	 printf("variables!\n");
	 *candidate = NULL;
	 return(ERROR__NO_BRANCHING_CANDIDATE);
      }
#endif
      /* Free the candidates */
      if (candidates){
	 for (i=0; i<cand_num; i++){
	    free_candidate(candidates + i);
	 }
	 FREE(candidates);
      }
      *candidate = NULL;
      return(DO_NOT_BRANCH);

    case DO_BRANCH:
      break;

   case ERROR__NO_BRANCHING_CANDIDATE:
      *candidate = NULL;
      return(ERROR__NO_BRANCHING_CANDIDATE);
   }

   /* OK, now we have to branch. */

   /* First of all, send everything to the cutpool that hasn't been sent
      before and send the current node description to the TM. */
   p->comp_times.strong_branching += used_time(&p->tt);
#pragma omp critical(cut_pool)
   send_cuts_to_pool(p, -1);
   send_node_desc(p, NODE_BRANCHED_ON);
   p->comp_times.communication += used_time(&p->tt);

   /* Add all the branching cuts */
   if (p->par.branch_on_cuts)
      add_slacks_to_matrix(p, cand_num, candidates);
   m = lp_data->m;
   rows = lp_data->rows;

#ifndef MAX_CHILDREN_NUM
   /* The part below is not needed when we have MAX_CHILDREN_NUM specified */
   /* Count how many objval/termcode/feasible entry we might need
      and allocate space for it */
   for (maxnum = candidates[0]->child_num, j=0, i=1; i<cand_num; i++){
      if (maxnum < candidates[i]->child_num)
	 maxnum = candidates[i]->child_num;
   }

   objval   = (double *) malloc(maxnum * DSIZE);
   termcode = (int *) malloc(maxnum * ISIZE);
   feasible = (int *) malloc(maxnum * ISIZE);
   iterd    = (int *) malloc(maxnum * ISIZE);
#ifdef COMPILE_FRAC_BRANCHING
   frval = (double **) malloc(maxnum * sizeof(double *));
   pfrval = (double **) malloc(maxnum * sizeof(double *));
   frind = (int **) malloc(maxnum * sizeof(int *));
   pfrind = (int **) malloc(maxnum * sizeof(int *));
   frnum = (int *) malloc(maxnum * ISIZE);
   pfrnum = (int *) malloc(maxnum * ISIZE);
#endif
   pobj  = (double *) malloc(maxnum * DSIZE);
   pterm = (int *) malloc(maxnum * ISIZE);
   pfeas = (int *) malloc(maxnum * ISIZE);
   piter = (int *) malloc(maxnum * ISIZE);
#endif

   /* 
    * see if hot-starts should be used. in theory if strong branching is used
    * and only variable bounds are changed then, hot-starts should be faster
    */
   if (p->par.use_hot_starts && !p->par.branch_on_cuts) {
      should_use_hot_starts = TRUE;
   } else {
      should_use_hot_starts = FALSE;
   }

   /* Set the iteration limit */
   if (p->par.max_presolve_iter > 0)
      set_itlim(lp_data, p->par.max_presolve_iter);

   if (should_use_hot_starts) {
      mark_hotstart(lp_data);
      set_itlim_hotstart(lp_data, p->par.max_presolve_iter);
   }
   vars = lp_data->vars;

   /* Look at the candidates one-by-one and presolve them. */
   oldobjval   = lp_data->objval;

   st_time     = used_time(&total_time);
   st_time     = used_time(&total_time);
   total_iters = 0;
   //printf("number of candidates = %d\n",cand_num);
   //printf("lp time = %f\n",p->comp_times.lp);
   //printf ("used time = %f\n",p->tt);
   //printf ("str time = %f\n",p->comp_times.strong_branching);
   for (i=0; i<cand_num; i++){
      can = candidates[i];

#ifndef MAX_CHILDREN_NUM
      can->objval = pobj;
      can->termcode = pterm;
      can->feasible = pfeas;
      can->iterd = piter;
      if (p->tm->par.keep_description_of_pruned == KEEP_IN_MEMORY){
	 can->solutions = (double **) calloc(maxnum, sizeof(double *));	
	 can->sol_inds = (int **) calloc(maxnum, sizeof(int *));	
	 can->sol_size = (int *) calloc(maxnum, ISIZE);	
      }

#ifdef SENSITIVITY_ANALYSIS
      if (p->tm->par.sensitivity_analysis){      
	 can->duals = (double **) calloc(maxnum, sizeof(double *));
      }else{
	 can->duals = NULL;	 
      }
#endif
#ifdef COMPILE_FRAC_BRANCHING
      can->frac_num = pfrnum;
      can->frac_ind = pfrind;
      can->frac_val = pfrval;
#endif

#else
      if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){	 	 
	 can->solutions = (double **) calloc (MAX_CHILDREN_NUM, 
					      sizeof(double *));
	 can->sol_inds = (int **) calloc(MAX_CHILDREN_NUM, 
					 sizeof(int *));	
	 can->sol_sizes = (int *) calloc(MAX_CHILDREN_NUM, ISIZE);	
      }
#ifdef SENSITIVITY_ANALYSIS
      if (p->par.sensitivity_analysis){      
	 can->duals = (double **) calloc (MAX_CHILDREN_NUM, sizeof(double *));
      }else{
	 can->duals = NULL;	 
      }
#endif
#endif

#ifdef STATISTICS
      cnum += can->child_num;
#endif

      /* Now depending on the type, adjust ub/lb or rhs/range/sense */
      switch (can->type){
       case CANDIDATE_VARIABLE:
	 branch_var = can->position;
#if 0
	 if (lp_data->status[branch_var] & PERM_FIXED_TO_LB ||
	     lp_data->status[branch_var] & PERM_FIXED_TO_UB){
	 if (vars[branch_var]->lb == vars[branch_var]->ub){
	    printf("Warning -- branching candidate is already fixed. \n");
	    printf("SYMPHONY has encountered numerical difficulties \n");
	    printf("With the LP solver. Exiting...\n\n");
	 }
         /* } to unconfuse vi*/
#endif
	 lb = vars[branch_var]->lb;
	 ub = vars[branch_var]->ub;
	 for (j = 0; j < can->child_num; j++){
	    switch (can->sense[j]){
	     case 'E':
	       change_lbub(lp_data, branch_var, can->rhs[j], can->rhs[j]);
	       break;
	     case 'R':
	       change_lbub(lp_data, branch_var, can->rhs[j],
			   can->rhs[j] + can->range[j]);
	       break;
	     case 'L':
	       change_lbub(lp_data, branch_var, lb, can->rhs[j]);
	       break;
	     case 'G':
	       change_lbub(lp_data, branch_var, can->rhs[j], ub);
	       break;
	    }
	    check_ub(p);
	    /* The original basis is in lp_data->lpbas */
            if (should_use_hot_starts) {
               can->termcode[j] = solve_hotstart(lp_data, can->iterd+j);
               total_iters+=*(can->iterd+j);
            } else {
               can->termcode[j] = dual_simplex(lp_data, can->iterd+j);
               total_iters+=*(can->iterd+j);
            }
            p->lp_stat.lp_calls++;
	    can->objval[j] = lp_data->objval;
	    get_x(lp_data);

#ifdef SENSITIVITY_ANALYSIS
	    if (p->par.sensitivity_analysis){      
	       get_dj_pi(lp_data);
	       can->duals[j] = (double *) malloc (DSIZE*p->base.cutnum);
	       memcpy(can->duals[j], lp_data->dualsol, DSIZE*p->base.cutnum);
	    }
#endif

	    if (can->termcode[j] == LP_OPTIMAL){
	       /* is_feasible_u() fills up lp_data->x, too!! */
	       switch (is_feasible_u(p, TRUE, FALSE)){

		  /*NOTE: This is confusing but not all that citical...*/
		  /*The "feasible" field is only filled out for the
		    purposes of display (in vbctool) to keep track of
		    where in the tree the feasible solutions were
		    found. Since this may not be the actual candidate
		    branched on, we need to pass this info on to whatever
		    candidate does get branched on so the that the fact that
		    a feasible solution was found in presolve can be recorded*/

		case IP_FEASIBLE:
		  can->termcode[j] = LP_OPT_FEASIBLE;
		  can->feasible[j] = TRUE;
		  if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
		     can->solutions[j] = (double *) malloc (DSIZE*lp_data->n);
		     memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
		  }
		  break;
		  
		case IP_FEASIBLE_BUT_CONTINUE:
		  can->termcode[j] = LP_OPT_FEASIBLE_BUT_CONTINUE;
		  can->feasible[j] = TRUE;
		  if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
		     can->solutions[j] = (double *) malloc (DSIZE*lp_data->n);
		     memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
		  }
		  break;

		default:
		  break;
	       }
	    }
#ifdef COMPILE_FRAC_BRANCHING
	    else{
	       if (can->termcode[j] != LP_ABANDONED){
		  get_x(lp_data);
	       }
	    }
	    if (can->termcode[j] != LP_ABANDONED){
	       xind = lp_data->tmp.i1; /* n */
	       xval = lp_data->tmp.d; /* n */
	       can->frac_num[j] = collect_fractions(p, lp_data->x, xind, xval);
	       if (can->frac_num[j] > 0){
		  can->frac_ind[j] = (int *) malloc(can->frac_num[j] * ISIZE);
		  can->frac_val[j] = (double *) malloc(can->frac_num[j]*DSIZE);
		  memcpy(can->frac_ind[j], xind, can->frac_num[j] * ISIZE);
		  memcpy(can->frac_val[j], xval, can->frac_num[j] * DSIZE);
	       }
	    }else{
	       can->frac_num[j] = 0;
	    }
#endif
#ifdef STATISTICS
	    if (can->termcode[j] == LP_D_ITLIM)
	       itlim++;
#endif
	 }
	 change_lbub(lp_data, branch_var, lb, ub);
	 break;

       case CANDIDATE_CUT_IN_MATRIX:
	 branch_row = can->position;
	 for (j = 0; j < can->child_num; j++){
	    change_row(lp_data, branch_row,
		       can->sense[j], can->rhs[j], can->range[j]);
	    check_ub(p);
	    /* The original basis is in lp_data->lpbas */
	    can->termcode[j] = dual_simplex(lp_data, can->iterd+j);
            p->lp_stat.lp_calls++;
	    can->objval[j] = lp_data->objval;


	    get_x(lp_data);

#ifdef SENSITIVITY_ANALYSIS
	    if (p->par.sensitivity_analysis){      
	       get_dj_pi(lp_data);
	       can->duals[j] = (double *) malloc (DSIZE*p->base.cutnum);
	       memcpy(can->duals[j], lp_data->dualsol, DSIZE*p->base.cutnum);
	    }
#endif

	    if (can->termcode[j] == LP_OPTIMAL){
	       /* is_feasible_u() fills up lp_data->x, too!! */
	       switch (is_feasible_u(p, TRUE, FALSE)){

		  /*NOTE: This is confusing but not all that citical...*/
		  /*The "feasible" field is only filled out for the
		    purposes of display (in vbctool) to keep track of
		    where in the tree the feasible solutions were
		    found. Since this may not be the actual candidate
		    branched on, we need to pass this info on to whatever
		    candidate does get branched on so the that the fact that
		    a feasible solution was found in presolve can be recorded*/

		case IP_FEASIBLE:
		  can->termcode[j] = LP_OPT_FEASIBLE;
		  can->feasible[j] = TRUE;
		  if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
		     can->solutions[j] = (double *) malloc (DSIZE*
							    lp_data->n);
		     memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
		  }
		  break;
		  
		case IP_FEASIBLE_BUT_CONTINUE:
		  can->termcode[j] = LP_OPT_FEASIBLE_BUT_CONTINUE;
		  can->feasible[j] = TRUE;
		  if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
		     can->solutions[j] = (double *) malloc (DSIZE*
							    lp_data->n);
		     memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
		  }
		  break;
		  
		default:
		  break;
	       }
	    }
#ifdef COMPILE_FRAC_BRANCHING
	    else{
	       if (can->termcode[j] != LP_ABANDONED)
		  get_x(lp_data);
	    }
	    if (can->termcode[j] != LP_ABANDONED){
	       xind = lp_data->tmp.i1; /* n */
	       xval = lp_data->tmp.d; /* n */
	       can->frac_num[j] = collect_fractions(p, lp_data->x, xind, xval);
	       if (can->frac_num[j] > 0){
		  can->frac_ind[j] = (int *) malloc(can->frac_num[j] * ISIZE);
		  can->frac_val[j] = (double *) malloc(can->frac_num[j]*DSIZE);
		  memcpy(can->frac_ind[j], xind, can->frac_num[j] * ISIZE);
		  memcpy(can->frac_val[j], xval, can->frac_num[j] * DSIZE);
	       }
	    }else{
	       can->frac_num[j] = 0;
	    }
#endif
#ifdef STATISTICS
	    if (can->termcode[j] == LP_D_ITLIM)
	       itlim++;
#endif
	 }
	 cut = rows[branch_row].cut;
	 change_row(lp_data, branch_row, cut->sense, cut->rhs, cut->range);
	 free_row_set(lp_data, 1, &branch_row);
	 break;
      }

      switch ((j = compare_candidates_u(p, oldobjval, best_can, can))){
       case FIRST_CANDIDATE_BETTER:
       case FIRST_CANDIDATE_BETTER_AND_BRANCH_ON_IT:
	  if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
	     min_ind = -1;
	     for (k = can->child_num - 1; k >= 0; k--){
		if (can->feasible[k]){
		   if (min_ind < 0){
		      min_obj = SYM_INFINITY;
		      for (l = best_can->child_num - 1; l >= 0; l--){
			 if (best_can->feasible[l] && best_can->objval[k] < 
			     min_obj){
			    min_obj = best_can->objval[l]; 
			    min_ind = l;		      
			 }
		      }
		   }		   
		   if (min_ind > -1){
		      if(can->objval[k] > best_can->objval[min_ind]){
			 best_can->feasible[k] = TRUE;
			 best_can->solutions[k] = can->solutions[k];
			 can->solutions[k] = 0;
			 min_ind = -1;
		      }
		   }
		}
	     }
	  } else{
	     for (k = best_can->child_num - 1; k >= 0; k--){
		/* Again, this is only for tracking that there was a feasible
		   solution discovered in presolve for display purposes */
		if (can->feasible[k]){
		   best_can->feasible[k] = TRUE;
		}
	     }
	  }
	  free_candidate(candidates + i);
	  break;
       case SECOND_CANDIDATE_BETTER:
       case SECOND_CANDIDATE_BETTER_AND_BRANCH_ON_IT:
#ifndef MAX_CHILDREN_NUM
	 if (best_can == NULL){
	    pobj  = objval;
	    pterm = termcode;
	    pfeas = feasible;
	    piter = iterd;
#ifdef COMPILE_FRAC_BRANCHING
	    pfrnum = frnum;
	    pfrind = frind;
	    pfrval = frval;
#endif
	 }else{
	    pobj  = best_can->objval;
	    pterm = best_can->termcode;
	    pfeas = best_can->feasible;
	    piter = best_can->iterd;
#ifdef COMPILE_FRAC_BRANCHING
	    pfrnum = best_can->frac_num;
	    pfrind = best_can->frac_ind;
	    pfrval = best_can->frac_val;
#endif
	 }
#endif
	 if (best_can){
	    if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
	       min_ind = -1;
	       for (k = best_can->child_num - 1; k >= 0; k--){
		  if (best_can->feasible[k]){
		     if (min_ind < 0){
			min_obj = SYM_INFINITY;
			for (l = can->child_num - 1; l >= 0; l--){
			   if (can->feasible[l] && can->objval[k] < 
			       min_obj){
			      min_obj = can->objval[l]; 
			      min_ind = l;		      
			   }
			}
		     }		   
		     if (min_ind > -1){
			if(best_can->objval[k] > can->objval[min_ind]){
			   can->feasible[k] = TRUE;
			   can->solutions[k] = best_can->solutions[k];
			   best_can->solutions[k] = 0;
			   min_ind = -1;
			}
		     }
		  }
	       }
	    }	
	    else{
	       for (k = can->child_num - 1; k >= 0; k--){
		  /* Again, this is only for tracking that there was a feasible
		     solution discovered in presolve for display purposes */
		  if (best_can->feasible[k]){
		     can->feasible[k] = TRUE;
		  }
	       }
	    }
	    free_candidate(&best_can);
	 }
	 best_can = can;
	 candidates[i] = NULL;
	 break;
      }
      if ((j & BRANCH_ON_IT)){
	 break;
      }
      st_time += used_time(&total_time);
      should_continue_strong_branching(p,i,cand_num,st_time,total_iters,
            &should_continue);
      if (should_continue==FALSE) {
         PRINT(p->par.verbosity, 0, 
               ("too much time in strong branching, breaking\n"));
         break;
      }
   }
   //printf ("total_iters = %d \n",total_iters);
   //printf ("candidates evaluated = %d \n",i);

   if (should_use_hot_starts) {
      unmark_hotstart(lp_data);
   }

#if 0
   if (best_can->type == CANDIDATE_VARIABLE &&
       vars[best_can->position]->lb == vars[best_can->position]->ub){
      printf("Error -- branching variable is already fixed. \n");
      printf("SYMPHONY has encountered numerical difficulties \n");
      printf("with the LP solver. Exiting...\n\n");
   }
   
   if (best_can->type == CANDIDATE_VARIABLE &&
       best_can->rhs[0] < vars[best_can->position]->lb ||
       best_can->rhs[1] > vars[best_can->position]->ub){
      printf("Warning -- possible illegal branching. \n");
      printf("SYMPHONY has encountered possible numerical difficulties \n");
      printf("with the LP solver. Exiting...\n\n");
   }
#endif
   
#ifndef MAX_CHILDREN_NUM
   FREE(pobj); FREE(pterm); FREE(pfeas); FREE(piter);
#  ifdef COMPILE_FRAC_BRANCHING
   FREE(pfrnum); FREE(pfrind); FREE(pfrval);
#  endif
#endif

   if (p->par.max_presolve_iter > 0)
      set_itlim(lp_data, -1);

#ifdef STATISTICS
   PRINT(p->par.verbosity, 5,
	 ("Itlim reached %i times out of %i .\n\n", itlim, cnum));
#endif

   for (i++; i<cand_num; i++){
      /* Free the remaining candidates */
      free_candidate(candidates + i);
   }
   FREE(candidates);
   if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
      indices = lp_data->tmp.i1;
      values = lp_data->tmp.d;
      for (k = best_can->child_num - 1; k >= 0; k--){
	 if (best_can->feasible[k]){
	    best_can->sol_sizes[k] = collect_nonzeros(p, 
						      best_can->solutions[k],
						     indices, values);
	    FREE(best_can->solutions[k]);
	    best_can->sol_inds[k] = (int *) malloc(best_can->sol_sizes[k]* 
						   ISIZE);
	    best_can->solutions[k] = (double *) malloc(best_can->sol_sizes[k]* 
						       DSIZE);
	    memcpy(best_can->sol_inds[k], indices, best_can->sol_sizes[k] * 
		   ISIZE);
	    memcpy(best_can->solutions[k], values, best_can->sol_sizes[k]* 
		   DSIZE);
	    break;
	 }
      }
   }
   
   *candidate = best_can;
   
   return(DO_BRANCH);
}

/*===========================================================================*/

int branch(lp_prob *p, int cuts)
{
   LPdata *lp_data = p->lp_data;
   branch_obj *can;
   char *action;
   int branch_var, branch_row, keep;
   var_desc *var;
   cut_data *cut;
   node_desc *desc;
   int termcode;

   termcode = select_branching_object(p, &cuts, &can);

   if (termcode == ERROR__NO_BRANCHING_CANDIDATE){
      return(termcode);
   }
   
   if (can == NULL){
      if (termcode == DO_NOT_BRANCH__FEAS_SOL) {
         /* a better feasible solution was found. return to do reduced cost
          * fixing etc.
          */
         return(FEAS_SOL_FOUND);
      }
      /* We were either able to fathom the node or found violated cuts
       * In any case, send the qualifying cuts to the cutpool */
      p->comp_times.strong_branching += used_time(&p->tt);
#pragma omp critical(cut_pool)
      send_cuts_to_pool(p, p->par.eff_cnt_before_cutpool);
      p->comp_times.communication += used_time(&p->tt);
      return( termcode == DO_NOT_BRANCH__FATHOMED ? FATHOMED_NODE : cuts );
   }

   /*------------------------------------------------------------------------*\
    * Now we evaluate can, the best of the candidates.
   \*------------------------------------------------------------------------*/
   action = lp_data->tmp.c; /* n (good estimate... can->child_num) */
   if ((termcode = select_child_u(p, can, action)) < 0)
      return(termcode);
   if (p->par.verbosity > 4)
      print_branch_stat_u(p, can, action);

   for (keep = can->child_num-1; keep >= 0; keep--)
      if (action[keep] == KEEP_THIS_CHILD) break;

   /* Send the branching information to the TM and inquire whether we
      should dive */
   p->comp_times.strong_branching += used_time(&p->tt);
   send_branching_info(p, can, action, &keep);
   p->comp_times.communication += used_time(&p->tt);

   /* If we don't dive then return quickly */
   if (keep < 0 || p->dive == DO_NOT_DIVE){
      free_candidate_completely(&can);
      return(FATHOMED_NODE);
   }

   desc = p->desc;
   switch (can->type){
    case CANDIDATE_VARIABLE:
      var = lp_data->vars[branch_var = can->position];
      switch (can->sense[keep]){
       case 'E':
	 var->new_lb = var->new_ub = can->rhs[keep];
	 var->lb = var->ub = can->rhs[keep];                             break;
       case 'R':
	 var->new_lb = can->rhs[keep]; 
         var->new_ub = var->lb + can->range[keep];
	 var->lb = can->rhs[keep]; var->ub = var->lb + can->range[keep]; break;
       case 'L':
	 var->new_ub = can->rhs[keep];
	 var->ub = can->rhs[keep];                                       break;
       case 'G':
	 var->new_lb = can->rhs[keep];
	 var->lb = can->rhs[keep];                                       break;
      }
      change_col(lp_data, branch_var, can->sense[keep], var->lb, var->ub);
      lp_data->status[branch_var] |= VARIABLE_BRANCHED_ON;
      break;
    case CANDIDATE_CUT_IN_MATRIX:
      branch_row = can->position;
      cut = lp_data->rows[branch_row].cut;
      /* To maintain consistency with TM we have to fix a few more things if
	 we had a non-base, new branching cut */
      if (branch_row >= p->base.cutnum && !(cut->branch & CUT_BRANCHED_ON)){
	 /* insert cut->name into p->desc.cutind.list, and insert SLACK_BASIC
	    to the same position in p->desc.basis.extrarows.stat */
#ifdef DO_TESTS
	 if (desc->cutind.size != desc->basis.extrarows.size){
	    printf("Oops! desc.cutind.size != desc.basis.extrarows.size! \n");
	    exit(-123);
	 }
#endif
#ifdef COMPILE_IN_LP
	 /* Because these cuts are shared with the treemanager, we have to
	    make a copy before changing them if the LP is compiled in */
	 cut = (cut_data *) malloc(sizeof(cut_data));
	 memcpy((char *)cut, (char *)lp_data->rows[branch_row].cut,
		sizeof(cut_data));
	 if (cut->size){
	    cut->coef = (char *) malloc(cut->size);
	    memcpy((char *)cut->coef,
		   (char *)lp_data->rows[branch_row].cut->coef, cut->size);
	 }
	 lp_data->rows[branch_row].cut = cut;
#endif
	 if (desc->cutind.size == 0){
	    desc->cutind.size = 1;
	    desc->cutind.list = (int *) malloc(ISIZE);
	    desc->cutind.list[0] = cut->name;
	    desc->basis.extrarows.size = 1; /* this must have been 0, too */
	    desc->basis.extrarows.stat = (int *) malloc(ISIZE);
	    desc->basis.extrarows.stat[0] = SLACK_BASIC;
	 }else{
	    int i, name = cut->name;
	    int *list;
	    int *stat;
	    /* most of the time the one extra element will fit into the
	       already allocated memory chunk, so it's worth to realloc */
	    desc->cutind.size++;
	    list = desc->cutind.list =
	       (int *) realloc(desc->cutind.list, desc->cutind.size * ISIZE);
	    desc->basis.extrarows.size++;
	    stat = desc->basis.extrarows.stat =
	       (int *) realloc(desc->basis.extrarows.stat,
			       desc->cutind.size * ISIZE);
	    for (i = desc->cutind.size - 1; i > 0; i--){
#ifdef DO_TESTS
	       if (name == list[i-1]){
		  printf("Oops! name == desc.cutind.list[i] !\n");
		  exit(-124);
	       }
#endif
	       if (name < list[i-1]){
		  list[i] = list[i-1];
		  stat[i] = stat[i-1];
	       }else{
		  break;
	       }
	    }
	    list[i] = name;
	    stat[i] = SLACK_BASIC;
	 }
      }
      cut->rhs = can->rhs[keep];
      if ((cut->sense = can->sense[keep]) == 'R')
	 cut->range = can->range[keep];
      cut->branch = CUT_BRANCHED_ON | can->branch[keep];
      constrain_row_set(lp_data, 1, &branch_row);
      lp_data->rows[branch_row].free = FALSE;
      break;
   }

   /* Since this is a child we dived into, we know that TM stores the stati of
      extra vars/rows wrt the parent */
   p->desc->basis.extravars.type = WRT_PARENT;
   p->desc->basis.extrarows.type = WRT_PARENT;

   free_candidate_completely(&can);
   
   /* the new p->bc_index is received in send_branching_info() */
   p->bc_level++;
   /*p->iter_num = 0;*/

   return(NEW_NODE);
}
   
/*===========================================================================*/

int col_gen_before_branch(lp_prob *p, int *new_vars)
{
   our_col_set *new_cols;
   int dual_feas;

   check_ub(p);
   if (! p->has_ub ||
       (p->colgen_strategy & BEFORE_BRANCH__DO_NOT_GENERATE_COLS) ||
       (p->lp_data->nf_status & NF_CHECK_NOTHING))
      return(DO_BRANCH);

   PRINT(p->par.verbosity, 2, ("Generating cols before branching.\n"));
   p->comp_times.strong_branching += used_time(&p->tt);
   new_cols = price_all_vars(p);
   p->comp_times.pricing += used_time(&p->tt);
   /*price_all_vars sorts by user_ind. We need things sorted by colind */
   colind_sort_extra(p);
   *new_vars = new_cols->num_vars + new_cols->rel_ub + new_cols->rel_lb;
   dual_feas = new_cols->dual_feas;
   free_col_set(&new_cols);
   check_ub(p);
   if (dual_feas == NOT_TDF){
      return(DO_NOT_BRANCH);
   }else{
      if (p->ub - p->par.granularity < p->lp_data->objval ||
	  p->lp_data->termcode == LP_D_OBJLIM ||
	  p->lp_data->termcode == LP_OPT_FEASIBLE){
	 /* If total dual feas and high cost or feasibility ==> fathomable */
	 PRINT(p->par.verbosity, 1, ("Managed to fathom the node.\n"));
	 send_node_desc(p, p->lp_data->termcode == LP_OPT_FEASIBLE ?
			FEASIBLE_PRUNED : OVER_UB_PRUNED);
	 p->comp_times.communication += used_time(&p->tt);
	 return(DO_NOT_BRANCH__FATHOMED);
      }else{
	 return(DO_BRANCH); /* if we got here, then DO_BRANCH */
      }
   }
   return(DO_BRANCH); /* fake return */
}

/*===========================================================================*/

/*****************************************************************************/
/* This is a generic function                                                */
/*****************************************************************************/

void branch_close_to_half(lp_prob *p, int max_cand_num, int *cand_num,
			  branch_obj ***candidates)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol, lpetol1 = 1 - lpetol;
   int *xind = lp_data->tmp.i1; /* n */
   double fracx, *xval = lp_data->tmp.d; /* n */
   branch_obj *cand;
   int i, j, cnt = 0;
   double lim[7] = {.1, .15, .20, .233333, .266667, .3, 1};
   var_desc **vars = lp_data->vars;

   /* first get the fractional values */
   for (i = lp_data->n-1; i >= 0; i--){
      /*FIXME: This is a quick-fix to allow variables without upper bounds*/
      /* Not sure what I meant with this */
      /* if (lp_data->vars[i]->ub < 2.0){ */
      if (vars[i]->is_int){
	 if (x[i] > vars[i]->lb && x[i] < vars[i]->ub){
	    fracx = x[i] - floor(x[i]);
	    if (fracx > lpetol && fracx < lpetol1){
	       xind[cnt] = i;
	       xval[cnt++] = fabs(fracx - .5);
	    }
	 }
      }
      /*}*/
   }
   qsort_di(xval, xind, cnt);

   if (p->bc_level>p->par.strong_br_all_candidates_level || 
         p->par.user_set_strong_branching_cand_num) {
      for (j = 0, i = 0; i < cnt;){
         if (xval[i] > lim[j]){
            if (i == 0){
               j++; continue;
            }else{
               break;
            }
         }else{
            i++;
         }
      }
      cnt = i;
      *cand_num = MIN(max_cand_num, cnt);
   } else {
      *cand_num = cnt;
   }

   if (!*candidates)
      *candidates = (branch_obj **) malloc(*cand_num * sizeof(branch_obj *));
   for (i=*cand_num-1; i>=0; i--){
      cand = (*candidates)[i] = (branch_obj *) calloc(1, sizeof(branch_obj) );
      cand->type = CANDIDATE_VARIABLE;
      cand->child_num = 2;
      cand->position = xind[i];
      cand->sense[0] = 'L';
      cand->sense[1] = 'G';
      cand->rhs[0] = floor(x[xind[i]]);
      cand->rhs[1] = cand->rhs[0] + 1;
      cand->range[0] = cand->range[1] = 0;
   }
}

/*===========================================================================*/

/*****************************************************************************/
/* This is a generic function                                                */
/*****************************************************************************/

void branch_close_to_half_and_expensive(lp_prob *p, int max_cand_num,
					int *cand_num, branch_obj ***candidates)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol, lpetol1 = 1 - lpetol;
   int *xind = lp_data->tmp.i1; /* n */
   double fracx, *xval = lp_data->tmp.d; /* n */
   branch_obj *cand;
   int i, j, cnt = 0;
   double lim[7] = {.1, .15, .20, .233333, .266667, .3, 1};

   /* first get the fractional values */
   for (i = lp_data->n-1; i >= 0; i--){
      fracx = x[i] - floor(x[i]);
      if (fracx > lpetol && fracx < lpetol1){
	 xind[cnt] = i;
	 xval[cnt++] = fabs(fracx - .5);
       }
   }
   qsort_di(xval, xind, cnt);

   for (j=0, i=0; i<cnt; i++){
      if (xval[i] > lim[j]){
	 if (i == 0){
	    j++; continue;
	 }else{
	    break;
	 }
      }
   }
   cnt = i;

   if (max_cand_num >= cnt){
      *cand_num = cnt;
   }else{
      for (i=cnt-1; i>=0; i--){
	 get_objcoef(p->lp_data, xind[i], xval+i);
	 xval[i] *= -1;
      }
      qsort_di(xval, xind, cnt);
      *cand_num = max_cand_num;
   }

   if (!*candidates)
      *candidates = (branch_obj **) malloc(*cand_num * sizeof(branch_obj *));
   for (i=*cand_num-1; i>=0; i--){
      cand = (*candidates)[i] = (branch_obj *) calloc(1, sizeof(branch_obj) );
      cand->type = CANDIDATE_VARIABLE;
      cand->child_num = 2;
      cand->position = xind[i];
      cand->sense[0] = 'L';
      cand->sense[1] = 'G';
      cand->rhs[0] = floor(x[xind[i]]);
      cand->rhs[1] = cand->rhs[0] + 1;
      cand->range[0] = cand->range[1] = 0;
   }
}

/*===========================================================================*/

/*****************************************************************************/
/* This works only for 0/1 problems!!!                                       */
/*****************************************************************************/

void branch_close_to_one_and_cheap(lp_prob *p, int max_cand_num, int *cand_num,
				   branch_obj ***candidates)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol, lpetol1 = 1 - lpetol;
   int *xind = lp_data->tmp.i1; /* n */
   double *xval = lp_data->tmp.d; /* n */
   branch_obj *cand;
   int i, j, cnt = 0;
   double lim[8] = {.1, .2, .25, .3, .333333, .366667, .4, 1};

   /* first get the fractional values */
   for (i = lp_data->n-1; i >= 0; i--){
      if (x[i] > lpetol && x[i] < lpetol1){
	 xind[cnt] = i;
	 xval[cnt++] = 1 - x[i];
      }
   }
   qsort_di(xval, xind, cnt);

   for (j=0, i=0; i<cnt; i++){
      if (xval[i] > lim[j]){
	 if (i == 0){
	    j++; continue;
	 }else{
	    break;
	 }
      }
   }
   cnt = i;

   if (max_cand_num >= cnt){
      *cand_num = cnt;
   }else{
      for (i=cnt-1; i>=0; i--){
	 get_objcoef(p->lp_data, xind[i], xval+i);
      }
      qsort_di(xval, xind, cnt);
      *cand_num = max_cand_num;
   }

   if (!*candidates)
      *candidates = (branch_obj **) malloc(*cand_num * sizeof(branch_obj *));
   for (i=*cand_num-1; i>=0; i--){
      cand = (*candidates)[i] = (branch_obj *) calloc(1, sizeof(branch_obj) );
      cand->type = CANDIDATE_VARIABLE;
      cand->child_num = 2;
      cand->position = xind[i];
      cand->sense[0] = 'L';
      cand->sense[1] = 'G';
      cand->rhs[0] = floor(x[xind[i]]);
      cand->rhs[1] = cand->rhs[0] + 1;
      cand->range[0] = cand->range[1] = 0;
   }
}

/*===========================================================================*/
int should_continue_strong_branching(lp_prob *p, int i, int cand_num,
                                     double st_time, int total_iters, 
                                     int *should_continue)
{
   double allowed_time = 0;
   *should_continue = TRUE;
   int min_cands;
   int verbosity = p->par.verbosity;
   //verbosity = 20;

   if (p->bc_level<p->par.strong_br_all_candidates_level) {
      allowed_time = p->comp_times.lp/pow(2.0,p->bc_level+1)/10;
      min_cands = MIN(cand_num,p->par.strong_branching_cand_num_max);
      if (allowed_time < 2) {
         allowed_time = 2;
      }
   } else {
      allowed_time = p->comp_times.lp/2 - p->comp_times.strong_branching;
      min_cands = MIN(cand_num,p->par.strong_branching_cand_num_min);
   }
   PRINT(verbosity,10,("allowed_time = %f\n",allowed_time));

   if (st_time/(i+1)*cand_num < allowed_time) {
      /* all cands can be evaluated in given time */
      *should_continue = TRUE;
   } else if (i >= min_cands-1 && st_time>allowed_time) {
      /* time is up and min required candidates have been evaluated */
      *should_continue = FALSE;
   } else if (p->par.user_set_max_presolve_iter == TRUE) {
      /* user specified a limit and we wont change it */
      *should_continue = TRUE;
   } else {
      /* we will not be able to evaluate all candidates in given time. we
       * reduce the number of iterations */
      double min_iters = 
         (allowed_time-st_time)*total_iters/st_time/(cand_num-i+1);
      if (min_iters<10) {
         /*
          * cant evaluate all candidates in given time with just 10 iters
          * as well. we have a choice: increase iters and do min_cands
          * or do ten iters and try to evaluate max possible num of cands.
          * we like the second option more.
          */
         min_iters = 10;
      }
      if (p->par.use_hot_starts && !p->par.branch_on_cuts) {
         set_itlim_hotstart(p->lp_data, (int) min_iters);
      } else {
         set_itlim(p->lp_data, (int) min_iters);
      }
      PRINT(verbosity,6, ("iteration limit set to %d\n", (int )min_iters));
      *should_continue = TRUE;
   }
   PRINT(verbosity,29, ("strong branching i = %d\n",i));
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
