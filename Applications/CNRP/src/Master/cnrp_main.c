/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* Capacitated Network Routing Problems.                                     */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER
#define USER_MAIN

/* system include files */
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* SYMPHONY include files */
#include "proccomm.h"
#include "timemeas.h"
#include "messages.h"
#include "BB_types.h"
#include "BB_macros.h"
#include "pack_cut.h"
#include "pack_array.h"
#include "master.h"

/* CNRP include files */
#include "cnrp_types.h"

/*===========================================================================*/

typedef struct SOLUTION_DATA{
   double fixed_cost;
   double variable_cost;
   double gamma;
   double tau;
   int *tree;
}solution_data;

typedef struct SOLUTION_PAIRS{
   int solution1;
   int solution2;
#ifdef BINARY_SEARCH
   double gamma1;
   double gamma2;
#endif
}solution_pairs;

/*===========================================================================*/

void cnrp_solve(problem *p, cut_pool **cp, base_desc *base, node_desc *root);

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the CNRP application.
\*===========================================================================*/

int main(int argc, char **argv)
{
   int i;
   problem *p;
   cut_pool **cp = NULL;
   double gamma, gamma0, gamma1, tau, slope;
   double start_time, t = 0;

   solution_data utopia1;
   solution_data utopia2;
   solution_data solutions[100];
   int numsolutions = 0;
   solution_pairs pairs[100];
   int numpairs = 0;
   int *tree;
   int solution1, solution2;
   double utopia_fixed, utopia_variable, ub = 0.0;
   cnrp_problem *cnrp;
   node_desc *root= NULL;
   base_desc *base = NULL;
   double tolerance = .2;
   
   start_time = wall_clock(NULL);

   setvbuf(stdout, (char *)NULL, _IOLBF, 0);

   printf("\n");
   printf("*******************************************************\n");
   printf("*   This is SYMPHONY Version 4.0                      *\n");
   printf("*   Copyright 2000-2003 Ted Ralphs                    *\n");
   printf("*   All Rights Reserved.                              *\n");
   printf("*   Distributed under the Common Public License 1.0   *\n");
   printf("*******************************************************\n");
   printf("\n");

   /* Initialize */

   (void) used_time(&t);

   p = get_problem_ptr(TRUE);

   initialize_u(p);

   cnrp = (cnrp_problem *)(p->user);
   
   /* Set the parameters */
   readparams_u(p, argc, argv);
   
   /* Get the problem data */
   io_u(p);
   
   /* Start up the graphics window*/
   init_draw_graph_u(p);
   
   p->comp_times.readtime += used_time(&t);
   
   /* Finds the upper and lower bounds for the problem */
   start_heurs_u(p);

   /*---------------------------------------------------------------------*\
    * Generate the base and root description
   \*---------------------------------------------------------------------*/
   
   base = (base_desc *) calloc(1, sizeof(base_desc));
   root = (node_desc *) calloc(1, sizeof(node_desc));
   
   initialize_root_node_u(p, base, root);

#ifdef SAVE_CUT_POOL
   if (p->par.tm_par.max_cp_num){
      cp = (cut_pool **) malloc(p->par.tm_par.max_cp_num*sizeof(cut_pool *));
      for (i = 0; i < p->par.tm_par.max_cp_num; i++){
	 cp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
	 cp[i]->par = p->par.cp_par;
	 CALL_USER_FUNCTION( user_send_cp_data(p->user, &cp[i]->user) );
      }
      get_cp_ptr(cp, 0);
   }
#endif
   
   /* First, calculate the utopia point */
   
   cnrp->lp_par.gamma = 1.0;
   cnrp->cg_par.tau = cnrp->lp_par.tau = 0.0;
      
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = 1.0 tau = 0.0 \n", gamma, tau);  
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Solve */
   cnrp_solve(p, cp, base, root);

   /* Store the solution */
   tree = solutions[numsolutions].tree = (int *) calloc(cnrp->vertnum-1,ISIZE);
   memcpy((char *)tree, cnrp->cur_sol_tree, cnrp->vertnum-1);
   solutions[numsolutions].gamma = 1.0;
   solutions[numsolutions].tau = 0.0;
   solutions[numsolutions].fixed_cost = cnrp->fixed_cost;
   solutions[numsolutions++].variable_cost = cnrp->variable_cost;
   utopia_fixed = cnrp->fixed_cost;
      
   cnrp->lp_par.gamma = 0.0;
   cnrp->cg_par.tau = cnrp->lp_par.tau = 1.0;
      
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = 0.0 tau = 1.0 \n", gamma, tau);  
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Solve */
   cnrp_solve(p, cp, base, root);
      
   /* Store the solution */
   tree = solutions[numsolutions].tree = (int *) calloc(cnrp->vertnum-1,ISIZE);
   memcpy((char *)tree, cnrp->cur_sol_tree, cnrp->vertnum-1);
   solutions[numsolutions].gamma = 0.0;
   solutions[numsolutions].tau = 1.0;
   solutions[numsolutions].fixed_cost = cnrp->fixed_cost;
   solutions[numsolutions++].variable_cost = cnrp->variable_cost;
   utopia_variable = cnrp->variable_cost;
   
   cnrp->utopia_variable = utopia_variable;
   cnrp->utopia_fixed = utopia_fixed;
   
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Utopia point has fixed cost %.3f and variable cost %.3f \n",
	  utopia_fixed, utopia_variable);
   printf("***************************************************\n");
   printf("***************************************************\n\n");
   
   /* Add the first pair to the list */
#ifdef BINARY_SEARCH
   pairs[numpairs].gamma1 = 1.0;
   pairs[numpairs].gamma2 = 0.0;
#endif
   pairs[numpairs].solution1 = 0;
   pairs[numpairs++].solution2 = 1;

   while (TRUE){

      if (!numpairs) break;
      
      solution1 = pairs[--numpairs].solution1;
      solution2 = pairs[numpairs].solution2;

#ifdef BINARY_SEARCH
      gamma = (pairs[numpairs].gamma1 + pairs[numpairs].gamma2)/2;
      printf("DEBUG: Processing interval [%.5f, %.5f] length = %.5f\n",
	     pairs[numpairs].gamma1, pairs[numpairs].gamma2,
	     fabs(pairs[numpairs].gamma1 - pairs[numpairs].gamma2));
#elif defined(FIND_NONDOMINATED_SOLUTIONS)
      gamma = (utopia_variable - solutions[solution1].variable_cost)/
	 (utopia_fixed - solutions[solution2].fixed_cost +
	  utopia_variable - solutions[solution1].variable_cost);
#else
      slope = (solutions[solution1].variable_cost -
	       solutions[solution2].variable_cost)/
	      (solutions[solution2].fixed_cost -
	       solutions[solution1].fixed_cost);
      gamma = slope/(1+slope);
#endif
      tau = 1 - gamma;
      
      cnrp->lp_par.gamma = gamma;
      cnrp->cg_par.tau = cnrp->lp_par.tau = tau;

      /* Find upper bound */

      p->has_ub = FALSE;
      p->ub = MAXDOUBLE;
#ifndef BINARY_SEARCH
      for (i = 0; i < numsolutions; i++){
#ifdef FIND_NONDOMINATED_SOLUTIONS
	 ub = MAX(gamma*(solutions[i].fixed_cost - utopia_fixed),
		  tau*(solutions[i].variable_cost - utopia_variable));
#else
	 ub = gamma*solutions[i].fixed_cost + tau*solutions[i].variable_cost;
#endif 
	 if (ub < p->ub){
	    p->has_ub = TRUE;
	    p->ub = ub - .001;
	 }
      }
#endif
      
      printf("***************************************************\n");
      printf("***************************************************\n");
      printf("Now solving with gamma = %.6f tau = %.6f \n", gamma, tau);  
      printf("***************************************************\n");
      printf("***************************************************\n\n");
      
      cnrp->fixed_cost = cnrp->variable_cost = 0.0;
      
      cnrp_solve(p, cp, base, root);
      
#ifdef BINARY_SEARCH
      if (cnrp->fixed_cost - solutions[solution1].fixed_cost < .001 &&
	  solutions[solution1].variable_cost - cnrp->variable_cost < .001){
	 if (pairs[numpairs].gamma1 - gamma > tolerance){
	    pairs[numpairs].solution1 = solution1;
	    pairs[numpairs].solution2 = solution2;
	    pairs[numpairs++].gamma1 = gamma;
	 }
	 continue;
      }
      if (solutions[solution2].fixed_cost - cnrp->fixed_cost < .001 &&
	  cnrp->variable_cost - solutions[solution2].variable_cost < .001){
	 if (gamma - pairs[numpairs].gamma2 > tolerance){
	    pairs[numpairs].solution1 = solution1;
	    pairs[numpairs].solution2 = solution2;
	    pairs[numpairs++].gamma2 = gamma;
	 }
	 continue;
      }
#else
      if (cnrp->fixed_cost == 0.0 && cnrp->variable_cost == 0.0){
	 continue;
      }
#endif
      
      if (numsolutions == 100){
	 printf("Maximum number of solutions exceeded\n\n");
	 exit(0);
      }
      /* Insert new solution */
      if ((fabs(solutions[0].fixed_cost - cnrp->fixed_cost) < .001) &&
	  (fabs(solutions[0].variable_cost - cnrp->variable_cost) >
	   .001)){
	 tree = solutions[0].tree;
	 memcpy((char *)tree, cnrp->cur_sol_tree, cnrp->vertnum-1);	 
	 solutions[0].fixed_cost = cnrp->fixed_cost;
	 solutions[0].variable_cost = cnrp->variable_cost;
	 if (numpairs + 2 > 100){
	    printf("Maximum number of solution pairs exceeded\n\n");
	    exit(0);
	 }
#ifdef BINARY_SEARCH
	 pairs[numpairs].gamma1 = gamma;
#endif
	 pairs[numpairs].solution1 = 0;
	 pairs[numpairs++].solution2 = 1;
      }else if ((fabs(solutions[numsolutions-1].variable_cost -
		      cnrp->variable_cost) < .001) &&
		(fabs(solutions[numsolutions-1].fixed_cost -
		      cnrp->fixed_cost) > .001)){
	 tree = solutions[numsolutions-1].tree;
	 memcpy((char *)tree, cnrp->cur_sol_tree, cnrp->vertnum-1);	 
	 solutions[numsolutions-1].fixed_cost = cnrp->fixed_cost;
	 solutions[numsolutions-1].variable_cost = cnrp->variable_cost;
	 if (numpairs + 2 > 100){
	    printf("Maximum number of solution pairs exceeded\n\n");
	    exit(0);
	 }
#ifdef BINARY_SEARCH
	 pairs[numpairs].gamma2 = gamma;
#endif
	 pairs[numpairs].solution1 = numsolutions-2;
	 pairs[numpairs++].solution2 = numsolutions-1;
      }else{
	 printf("DEBUG: Found new solution. \n");
#ifdef BINARY_SEARCH
	 pairs[numpairs+1].gamma2 = pairs[numpairs].gamma2;
	 pairs[numpairs+1].gamma1 = gamma;
	 pairs[numpairs].gamma2 = gamma;
#endif
	 pairs[numpairs].solution1 = solution1;
	 pairs[numpairs++].solution2 = solution2;
	 pairs[numpairs].solution1 = solution2;
	 pairs[numpairs++].solution2 = solution2+1;
	 for (i = numsolutions; i > solution2; i--){
	    solutions[i] = solutions[i-1];
	 }
	 numsolutions++;
	 tree = solutions[solution2].tree =
	    (int *) calloc(cnrp->vertnum-1, ISIZE);
	 memcpy((char *)tree, cnrp->cur_sol_tree, cnrp->vertnum-1);
	 solutions[solution2].gamma = gamma;
	 solutions[solution2].tau = tau;
	 solutions[solution2].fixed_cost = cnrp->fixed_cost;
	 solutions[solution2].variable_cost = cnrp->variable_cost;
	 if (numpairs + 2 > 100){
	    printf("Maximum number of solution pairs exceeded\n\n");
	    exit(0);
	 }
      }
   }
   
   printf("\n********************************************************\n");
#ifdef FIND_NONDOMINATED_SOLUTIONS
   printf(  "* Found complete set of non-dominated solutions!!!!!!! *\n");
#else
   printf(  "* Found complete set of supported solutions!!!!!!!     *\n");
#endif
   printf(  "* Now displaying stats...                              *\n");
   printf(  "********************************************************\n\n");

#ifdef SAVE_CUT_POOL
   for (i = 0; i < p->par.tm_par.max_cp_num; i++){
      p->comp_times.bc_time.cut_pool += cp[i]->cut_pool_time;
      p->stat.cuts_in_pool += cp[i]->cut_num;
      cp[i]->msgtag = YOU_CAN_DIE;
      cp_close(cp[i]);
   }
#endif
   
   print_statistics(&(p->comp_times.bc_time), &(p->stat), 0.0, 0.0, 0,
		    start_time);

   printf("***************************************************\n");
   printf("***************************************************\n");
#ifdef FIND_NONDOMINATED_SOLUTIONS
   printf("Displaying non-dominated solution values and breakpoints\n");  
#else
   printf("Displaying supported solution values and breakpoints\n");  
#endif
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   gamma0 = 1.0;
   for (i = 0; i < numsolutions - 1; i++){
#ifdef FIND_NONDOMINATED_SOLUTIONS
      gamma1 = (utopia_variable - solutions[i].variable_cost)/
	 (utopia_fixed - solutions[i+1].fixed_cost +
	  utopia_variable - solutions[i].variable_cost);
#else
      slope = (solutions[i].variable_cost -
	       solutions[i+1].variable_cost)/
	      (solutions[i+1].fixed_cost -
	       solutions[i].fixed_cost);
      gamma1 = slope/(1+slope);
#endif
      printf("Fixed Cost: %.3f Variable Cost: %.3f ",
	     solutions[i].fixed_cost, solutions[i].variable_cost);
      printf("Optimal Range: %.2f - %.2f\n", gamma1, gamma0);
      gamma0 = gamma1;
   }
   printf("Fixed Cost: %.3f Variable Cost: %.3f ",
	  solutions[i].fixed_cost, solutions[i].variable_cost);
   printf("Optimal Range: %.2f - %.2f\n", 0.0, gamma0);
   
   FREE(root->desc);
   FREE(root->uind.list);
   FREE(root->not_fixed.list);
   FREE(root->cutind.list);
   FREE(root);
   FREE(base->userind);
   FREE(base);
   for (i = 0 ; i < numsolutions; i++){
      FREE(solutions[i].tree);
   }
   FREE(cp);

   return(0);
}

void cnrp_solve(problem *p, cut_pool **cp, base_desc *base, node_desc *root)
{
   int termcode;
   double t = 0, start_time;

   tm_prob *tm;

   /*---------------------------------------------------------------------*\
    * Initialize
   \*---------------------------------------------------------------------*/

   start_time = wall_clock(NULL);
   
   send_lp_data_u(p, 0, base);
   send_cg_data_u(p, 0);
   tm = get_tm_ptr(FALSE);
#ifdef SAVE_CUT_POOL
   tm->cpp = cp;
#else
   send_cp_data_u(p, 0);
#endif
   
   tm_initialize(base, root, 0);

   /*---------------------------------------------------------------------*\
    * Solve the problem and receive solutions                         
   \*---------------------------------------------------------------------*/
   
   tm->start_time = start_time;

   termcode = solve(tm);

#ifdef SAVE_CUT_POOL
   /* Save the cut pool from being wiped out */
   tm->cpp = NULL;
#endif
   
   tm_close(tm, termcode);
   
   /*---------------------------------------------------------------------*\
    * Display the the results and solution data                               
   \*---------------------------------------------------------------------*/
   
   if (p->par.verbosity > 0){
      if (termcode == TM_FINISHED){
	 printf("\n****************************************************\n");
	 printf(  "* Branch and Cut Finished!!!!!!!                   *\n");
	 printf(  "* Now displaying stats and optimal solution...     *\n");
	 printf(  "****************************************************\n\n");
      }else{
	 printf("\n****************************************************\n");
	 printf(  "* Time Limit Exceeded :(                           *\n");
	 printf(  "* Now displaying stats and best solution...        *\n");
	 printf(  "****************************************************\n\n");
      }
      
      print_statistics(&(tm->comp_times), &(tm->stat), tm->ub, tm->lb, 0,
		       start_time);
   }

   display_solution_u(p, 0);

   /* Keep cumulative statistics */
   if (tm->stat.max_depth > p->stat.max_depth){
      p->stat.max_depth = tm->stat.max_depth;
   }
   p->stat.chains += tm->stat.chains;
   p->stat.diving_halts += tm->stat.diving_halts;
   p->stat.tree_size += tm->stat.tree_size;
   p->stat.created += tm->stat.created;
   p->stat.analyzed += tm->stat.analyzed;
   p->stat.leaves_before_trimming += tm->stat.leaves_before_trimming;
   p->stat.leaves_after_trimming += tm->stat.leaves_after_trimming;
   p->stat.vars_not_priced += tm->stat.vars_not_priced;

   p->comp_times.bc_time.communication += tm->comp_times.communication;
   p->comp_times.bc_time.lp += tm->comp_times.lp;
   p->comp_times.bc_time.separation += tm->comp_times.separation;
   p->comp_times.bc_time.fixing += tm->comp_times.fixing;
   p->comp_times.bc_time.pricing += tm->comp_times.pricing;
   p->comp_times.bc_time.strong_branching += tm->comp_times.strong_branching;
   p->comp_times.bc_time.wall_clock_lp += tm->comp_times.wall_clock_lp;
   p->comp_times.bc_time.ramp_up_tm += tm->comp_times.ramp_up_tm;
   p->comp_times.bc_time.ramp_up_lp += tm->comp_times.ramp_up_lp;
   p->comp_times.bc_time.ramp_down_time += tm->comp_times.ramp_down_time;
   p->comp_times.bc_time.idle_diving += tm->comp_times.idle_diving;
   p->comp_times.bc_time.idle_node += tm->comp_times.idle_node;
   p->comp_times.bc_time.idle_names += tm->comp_times.idle_names;
   p->comp_times.bc_time.idle_cuts += tm->comp_times.idle_cuts;
   p->comp_times.bc_time.start_node += tm->comp_times.start_node;

#if 0
   if (p->par.do_draw_graph){
      s_bufid = init_send(DataInPlace);
      send_msg(p->dg_tid, CTOI_YOU_CAN_DIE);
      freebuf(s_bufid);
   }
#endif

   free_tm(tm);
}

