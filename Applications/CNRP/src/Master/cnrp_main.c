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

/* system include files */
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

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
   int fixed_cost;
   int variable_cost;
   double gamma;
   double tau;
   int *tree;
}solution_data;

typedef struct SOLUTION_PAIRS{
   int solution1;
   int solution2;
}solution_pairs;

/*===========================================================================*/

void cnrp_solve(problem *p, int argc, char **argv, double gamma, double tau);

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the CNRP application.
\*===========================================================================*/

int main(int argc, char **argv)
{
   int i;
   problem *p;
   double gamma, tau;

   solution_data utopia1;
   solution_data utopia2;
   solution_data solutions[100];
   int numsolutions = 0;
   solution_pairs pairs[100];
   int numpairs = 0;
   int *tree;
   int solution1, solution2;
   double utopia_fixed, utopia_variable;
   vrp_problem *vrp;
   
   setvbuf(stdout, (char *)NULL, _IOLBF, 0);

   printf("\n");
   printf("*******************************************************\n");
   printf("*   This is SYMPHONY Version 4.0                      *\n");
   printf("*   Copyright 2000-2003 Ted Ralphs                    *\n");
   printf("*   All Rights Reserved.                                     *\n");
   printf("*   Distributed under the Common Public License 1.0   *\n");
   printf("*******************************************************\n");
   printf("\n");

#if 0
   for (i = 0; i < 100; i++){
      solutions[i] = (solution_data *) calloc(1, sizeof(solution_data));
   }
#endif
   
   /* First, calculate the utopia point */
   
   p = get_problem_ptr(TRUE);

   gamma = 1.0;
   tau = 1 - gamma;
      
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = %.2f tau = %.2f \n", gamma, tau);  
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Solve */
   cnrp_solve(p, argc, argv,  gamma, tau);

   /* Store the solution */
   vrp = (vrp_problem *)(p->user);
   tree = solutions[numsolutions].tree = (int *) calloc(vrp->vertnum-1, ISIZE);
   memcpy((char *)tree, vrp->cur_sol_tree, vrp->vertnum-1);
   solutions[numsolutions].gamma = gamma;
   solutions[numsolutions].tau = tau;
   solutions[numsolutions].fixed_cost = (int) vrp->fixed_cost;
   solutions[numsolutions++].variable_cost = (int) vrp->variable_cost;
   utopia_fixed = vrp->fixed_cost;
      
   free_master_u(p);      
  
   p = get_problem_ptr(TRUE);

   gamma = 0.0;
   tau = 1 - gamma;
      
   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = %.2f tau = %.2f \n", gamma, tau);  
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Solve */
   cnrp_solve(p, argc, argv,  gamma, tau);
      
   /* Store the solution */
   vrp = (vrp_problem *)(p->user);
   tree = solutions[numsolutions].tree = (int *) calloc(vrp->vertnum-1, ISIZE);
   memcpy((char *)tree, vrp->cur_sol_tree, vrp->vertnum-1);
   solutions[numsolutions].gamma = gamma;
   solutions[numsolutions].tau = tau;
   solutions[numsolutions].fixed_cost = (int) vrp->fixed_cost;
   solutions[numsolutions++].variable_cost = (int) vrp->variable_cost;
   utopia_variable = vrp->variable_cost;

   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Utopia point has fixed cost %i and variable cost %i \n",
	  (int) utopia_fixed, (int) utopia_variable);
   printf("***************************************************\n");
   printf("***************************************************\n\n");
   
   free_master_u(p);

   /* Add the first pair to the list */
   pairs[numpairs].solution1 = 0;
   pairs[numpairs++].solution2 = 1;
   
   while (TRUE){

      if (!numpairs) break;
      
      p = get_problem_ptr(TRUE);

      solution1 = pairs[--numpairs].solution1;
      solution2 = pairs[numpairs].solution2;
	 
      gamma = (utopia_variable - solutions[solution1].variable_cost)/
	 (utopia_fixed - solutions[solution2].fixed_cost +
	  utopia_variable - solutions[solution1].variable_cost);
      
      tau = 1 - gamma;
      
      printf("***************************************************\n");
      printf("***************************************************\n");
      printf("Now solving with gamma = %.2f tau = %.2f \n", gamma, tau);  
      printf("***************************************************\n");
      printf("***************************************************\n\n");

      cnrp_solve(p, argc, argv,  gamma, tau);

      vrp = (vrp_problem *)(p->user);
      
      if ((int) vrp->fixed_cost > solutions[solution1].fixed_cost &&
	  (int) vrp->variable_cost < solutions[solution1].variable_cost &&
	  (int) vrp->fixed_cost < solutions[solution2].fixed_cost &&
	  (int) vrp->variable_cost > solutions[solution2].variable_cost ){
	 if (numsolutions == 100){
	    printf("Maximum number of solutions exceeded\n\n");
	    exit(0);
	 }
	 /* Insert new solution */
	 for (i = numsolutions; i > solution2; i--){
	    solutions[i] = solutions[i-1];
	 }
	 numsolutions++;
	 tree = solutions[solution2].tree =
	    (int *) calloc(vrp->vertnum-1, ISIZE);
	 memcpy((char *)tree, vrp->cur_sol_tree, vrp->vertnum-1);
	 solutions[solution2].gamma = gamma;
	 solutions[solution2].tau = tau;
	 solutions[solution2].fixed_cost = (int) vrp->fixed_cost;
	 solutions[solution2].variable_cost = (int) vrp->variable_cost;
	 if (numpairs + 2 > 100){
	    printf("Maximum number of solution pairs exceeded\n\n");
	    exit(0);
	 }
	 pairs[numpairs].solution1 = solution1;
	 pairs[numpairs++].solution2 = solution2;
	 pairs[numpairs].solution1 = solution2;
	 pairs[numpairs++].solution2 = solution2+1;
      }
	 
      free_master_u(p);      
   }

   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now printing non-dominated solution values\n");  
   printf("***************************************************\n");
   printf("***************************************************\n\n");
   
   for (i = 0; i < numsolutions; i++){
      printf("Fixed Cost: %5i Variable Cost: %5i\n", solutions[i].fixed_cost,
	     solutions[i].variable_cost);
   }
   
   return(0);
}

void cnrp_solve(problem *p, int argc, char **argv, double gamma, double tau)
{
   int termcode;
   double t = 0, total_time=0;
   double start_time;

   node_desc *root= NULL;
   base_desc *base = NULL;

   tm_prob *tm;

   /*---------------------------------------------------------------------*\
    *                         program starts                              
   \*---------------------------------------------------------------------*/

   (void) used_time(&t);
   
   start_time = wall_clock(NULL);
   
   initialize_u(p);
   
   /* Set the parameters */
   readparams_u(p, argc, argv);
   
   /* Get the problem data */
   io_u(p);
   
   /* Start up the graphics window*/
   init_draw_graph_u(p);
   
   p->comp_times.readtime = used_time(&t);
   
   /* Finds the upper and lower bounds for the problem */
   start_heurs_u(p);
   
   (void) used_time(&t);
   
   /*---------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*---------------------------------------------------------------------*/
   
   base = (base_desc *) calloc(1, sizeof(base_desc));
   root = (node_desc *) calloc(1, sizeof(node_desc));
   
   initialize_root_node_u(p, base, root);
   
   /*---------------------------------------------------------------------*\
    * Send out problem data if needed
   \*---------------------------------------------------------------------*/
   
   vrp_problem *vrp = (vrp_problem *)(p->user);
   vrp->lp_par.gamma = gamma;
   vrp->cg_par.tau = vrp->lp_par.tau = tau;
   
   send_lp_data_u(p, 0, base);
   send_cg_data_u(p, 0);
   send_cp_data_u(p, 0);
   
   tm = tm_initialize(base, root, 0);
   
   /*---------------------------------------------------------------------*\
    * Solve the problem and receive solutions                         
   \*---------------------------------------------------------------------*/
   
   tm->start_time += start_time;
   termcode = solve(tm);
   
   tm_close(tm, termcode);
   
   /*---------------------------------------------------------------------*\
    * Display the the results and solution data                               
   \*---------------------------------------------------------------------*/
   
#if 0
   if (termcode == TM_FINISHED){
      printf("\n****************************************************\n");
      printf(  "* Branch and Cut Finished!!!!!!!                   *\n");
      printf(  "* Now displaying stats and optimal solution...     *\n");
      printf(  "****************************************************\n\n");
   }else if (termcode == TIME_LIMIT_EXCEEDED){
      printf("\n****************************************************\n");
      printf(  "* Time Limit Exceeded :(                           *\n");
      printf(  "* Now displaying stats and best solution...        *\n");
      printf(  "****************************************************\n\n");
   }else{
      printf(
	     "***********Something has died -- halting the machine\n\n");
      printf(
	     "***********Printing out partial data\n\n");
   }
   
   total_time  = p->comp_times.readtime;
   total_time += p->comp_times.ub_overhead + p->comp_times.ub_heurtime;
   total_time += p->comp_times.lb_overhead + p->comp_times.lb_heurtime;
   
   printf( "====================== Misc Timing ========================\n");
   printf( "  Problem IO        %.3f\n", p->comp_times.readtime);
   printf( "  UB overhead:      %.3f\n", p->comp_times.ub_overhead);
   printf( "  UB runtime:       %.3f\n", p->comp_times.ub_heurtime);
   printf( "  LB overhead:      %.3f\n", p->comp_times.lb_overhead);
   printf( "  LB runtime:       %.3f\n", p->comp_times.lb_heurtime);
   
   if (tm->lb > p->lb) p->lb = tm->lb;

   print_statistics(&(tm->comp_times), &(tm->stat), tm->ub, p->lb,
		    total_time, start_time);

#endif
   
   display_solution_u(p, 0);

#if 0
   if (p->par.do_draw_graph){
      s_bufid = init_send(DataInPlace);
      send_msg(p->dg_tid, CTOI_YOU_CAN_DIE);
      freebuf(s_bufid);
   }
#endif
   
   FREE(root);
   FREE(base);
   free_tm(tm);
}

