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

/* system include files */
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* SYMPHONY include files */
#include "OsiSymSolverInterface.hpp"
#include "BB_macros.h"
#include "master.h"

/* CNRP include files */
#include "cnrp_const.h"
#include "cnrp_types.h"
#include "cnrp_io.h"
#include "small_graph.h"

void cnrp_load_problem(sym_environment *env, cnrp_problem *cnrp);

/*===========================================================================*\
 * This is main() for the CNRP
\*===========================================================================*/

int main(int argc, char **argv)
{
   cnrp_problem *cnrp = (cnrp_problem *) calloc(1, sizeof(cnrp_problem));

   sym_environment *env = sym_open_environment();

   sym_parse_command_line(env, argc, argv);

   sym_set_user_data(env, cnrp);

   /* FIXME: Get rid of env->par.param_file argument */
   cnrp_readparams(cnrp, env->par.param_file, argc, argv);

   cnrp_io(cnrp, cnrp->par.infile);

   if (cnrp->par.use_small_graph == LOAD_SMALL_GRAPH){
      read_small_graph(cnrp);
   }

   if (!cnrp->numroutes && cnrp->par.prob_type == VRP){
      printf("\nError: Number of trucks not specified or computed "
	     "for VRP\n\n");
      exit(1);
   }
   
   if (cnrp->numroutes > 1){
      printf("NUMBER OF TRUCKS: \t%i\n", cnrp->numroutes);
      printf("TIGHTNESS: \t\t%.2f\n",
	     cnrp->demand[0]/(cnrp->capacity*(double)cnrp->numroutes));
   }
   
   /* Selects the cheapest edges adjacent to each node for the base set */

   if (cnrp->par.use_small_graph == SAVE_SMALL_GRAPH){
      if (!cnrp->g) make_small_graph(cnrp, 0);
      save_small_graph(cnrp);
   }else if (!cnrp->g){
      make_small_graph(cnrp, 0);
   }

   /* FIXME: Eliminate this use of the environment */
   if (env->ub > 0){
      cnrp->cur_tour->cost = (int) env->ub;
   }else{
      cnrp->cur_tour->cost = MAXINT;
   }
   cnrp->cur_tour->numroutes = cnrp->numroutes;
   
   if (cnrp->par.use_small_graph == LOAD_SMALL_GRAPH){
      if (env->ub <= 0 && cnrp->cur_tour->cost > 0)
	 env->ub = (int)(cnrp->cur_tour->cost);
      cnrp->numroutes = cnrp->cur_tour->numroutes;
   }

#if 0
   if(cnrp->par.prob_tpye == BPP)
      env->ub = 1;
#endif
   
   if (env->ub > 0 && !(cnrp->par.prob_type == BPP)){
      printf("INITIAL UPPER BOUND: \t%i\n\n", (int)(env->ub));
   }else if (!(cnrp->par.prob_type == BPP)){
      printf("INITIAL UPPER BOUND: \tNone\n\n");
   }else{
      printf("\n\n");
   }

   cnrp_load_problem(env, cnrp);
   
   sym_solve(env);

   sym_close_environment(env);

   return(0);
}

/*===========================================================================*\
 * This is the function that creates the formulation
\*===========================================================================*/

void cnrp_load_problem(sym_environment* env, cnrp_problem *cnrp)
{
   int n, m, nz;
   char *is_int, *sense;
   int *matbeg, *matind;
   double *matval, *obj, *obj2, *rhs, *rngval;
   double *lb, *ub; 
   int *costs = cnrp->dist.cost;
   int *edges = cnrp->edges;
   int i, j, maxvars = 0;
   char resize = FALSE;
   int vertnum = cnrp->vertnum;
   int total_edgenum = vertnum*(vertnum-1)/2;
   char prob_type = cnrp->par.prob_type, od_const = FALSE, d_x_vars = FALSE;
   int v0, v1;
   double flow_capacity;
#ifdef ADD_CAP_CUTS 
   int basecutnum = (2 + od_const)*vertnum - 1 + 2*total_edgenum;
#elif defined(ADD_FLOW_VARS)
   int basecutnum = (2 + od_const)*vertnum - 1;
#else
   int basecutnum = (1 + od_const)*vertnum;
#endif
#ifdef ADD_X_CUTS
   basecutnum += total_edgenum;
#endif
#if defined(DIRECTED_X_VARS) && !defined(ADD_FLOW_VARS)
   int edgenum = total_edgenum/2;
#elif defined(ADD_FLOW_VARS)
#ifdef DIRECTED_X_VARS
   int edgenum = total_edgenum/4;

   flow_capacity = cnrp->capacity;
#else
   int edgenum = total_edgenum/3;

   if (cnrp->par.prob_type == CSTP || cnrp->par.prob_type == CTP)
      flow_capacity = cnrp->capacity;
   else
      flow_capacity = cnrp->capacity/2;
#endif
#else
   int edgenum = total_edgenum;
#endif
   
   /* set up the inital LP data */

   n = total_edgenum;
   
   /*Estimate the number of nonzeros*/
#ifdef ADD_CAP_CUTS
   nz = 12*edgenum;
#elif defined(ADD_FLOW_VARS)
   nz = 8*edgenum;
#else
   nz = 3*edgenum;
#endif
#ifdef ADD_X_CUTS
   nz += 2*edgenum;
#endif

   /* Allocate the arrays. These are owned by SYMPHONY after returning. */
   matbeg  = (int *) malloc((n + 1) * ISIZE);
   matind  = (int *) malloc((nz) * ISIZE);
   matval  = (double *) malloc((nz) * DSIZE);
   obj     = (double *) malloc(n * DSIZE);
   obj2    = (double *) malloc(n * DSIZE);
   ub      = (double *) malloc(n * DSIZE);
   lb      = (double *) calloc(n, DSIZE); /* zero lower bounds */
   rhs     = (double *) malloc(m * DSIZE);
   sense   = (char *) malloc(m * CSIZE);
   rngval  = (double *) calloc(m, DSIZE);
   is_int  = (char *) calloc(n, CSIZE);

#ifdef DIRECTED_X_VARS
   /*whether or not we will have out-degree constraints*/
   od_const = (prob_type == VRP || prob_type == TSP || prob_type == BPP);
   d_x_vars = TRUE;
#endif
   
   for (i = 0, j = 0; i < total_edgenum; i++){
      is_int[i]    = TRUE;
      ub[i]        = 1.0;
      matbeg[i]    = j;
      obj[i]       = (double) costs[i];
      if (prob_type == CSTP || prob_type == CTP){
	 /*cardinality constraint*/
	 matind[j] = 0;
	 matval[j++] = 1.0;
      }
      /*in-degree constraint*/
      matval[j]    = 1.0;
      matind[j++]  = edges[2*i+1];
#ifdef DIRECTED_X_VARS
      /*out-degree constraint*/
      if (od_const){
	 matval[j]   = 1.0;
	 matind[j++] = vertnum + edges[2*i];
      }
#else
      if (prob_type == VRP || prob_type == TSP ||
	  prob_type == BPP || edges[2*i]){
	 matval[j]   = 1.0;
	 matind[j++] = edges[2*i];
      }
#endif	 
#ifdef ADD_CAP_CUTS
      v0 = edges[2*i];
      matval[j]    = -flow_capacity + (v0 ? cnrp->demand[v0] : 0);
      matind[j++]  = (2 + od_const)*vertnum - 1 + i;
#ifndef DIRECTED_X_VARS
      matval[j]    = -flow_capacity + cnrp->demand[edges[2*i + 1]];
      matind[j++]  = 2*cnrp->vertnum - 1 + total_edgenum + i;
#endif
#endif
#ifdef ADD_X_CUTS
      matval[j]    = 1.0;
      matind[j++]  = (2 + od_const)*vertnum-1 + 2*total_edgenum + i;
#endif
   }
#ifdef DIRECTED_X_VARS
   for (; i < 2*total_edgenum; i++){
      is_int[i]    = TRUE;
      ub[i]        = 1.0;
      matbeg[i]    = j;
      obj[i]       = (double) costs[i - total_edgenum];
      if (prob_type == CSTP || prob_type == CTP){
	 /*cardinality constraint*/
	 matind[j] = 0;
	 matval[j++] = 1.0;
      }
      /*in-degree constraint*/
      if (od_const || edges[2*(i - total_edgenum)]){
	 matval[j]   = 1.0;
	 matind[j++] = edges[2*(i - total_edgenum)];
      }
      /*out-degree constraint*/
      if (od_const){
	 matval[j]    = 1.0;
	 matind[j++]  = vertnum + edges[2*(i - total_edgenum)+1];
      }
#ifdef ADD_CAP_CUTS
      matval[j]    = -flow_capacity +
	 cnrp->demand[edges[2*(i - total_edgenum) + 1]];
      matind[j++]  = (2 + od_const)*vertnum - 1 + i;
#endif
#ifdef ADD_X_CUTS
      matval[j]    = 1.0;
      matind[j++]  = (2 + od_const)* vertnum - 1 + 2 * total_edgenum +
	 i - total_edgenum;
#endif
   }
#endif
   for (; i < (2+d_x_vars)*total_edgenum; i++){
      is_int[i] = FALSE;
      v0 = edges[2*(i-(1+d_x_vars)*total_edgenum)];
      ub[i] = flow_capacity - (v0 ? cnrp->demand[v0] : 0);
      matbeg[i]    = j;
      obj2[i]      = (double) costs[i- (1+d_x_vars)*total_edgenum];
#ifdef ADD_CAP_CUTS
      matval[j]    = 1.0;
      matval[j+1]  = 1.0;
      if (edges[2*(i-(1+d_x_vars)*total_edgenum)])
	 matval[j+2] = -1.0;
      matind[j++]  = (2 + od_const)*vertnum - 1 + i -
	 (1+d_x_vars)*total_edgenum;
      matind[j++]  = (1+od_const)*vertnum + edges[2*(i -
			  (1+d_x_vars)*total_edgenum) + 1] - 1;
      if (edges[2*(i - (1+d_x_vars)*total_edgenum)])
	 matind[j++] = (1+od_const)*vertnum + edges[2*(i -
			    (1+d_x_vars)*total_edgenum)] - 1;
#else
      matval[j]  = 1.0;
      if (edges[2*(i-(1+d_x_vars)*total_edgenum)])
	 matval[j+1] = -1.0;
      matind[j++]  = (1+od_const)*vertnum + edges[2*(i -
			  (1+d_x_vars)*total_edgenum) + 1] - 1;
      if (edges[2*(i - (1+d_x_vars)*total_edgenum)])
	 matind[j++] = (1+od_const)*vertnum + edges[2*(i -
			    (1+d_x_vars)*total_edgenum)] - 1;
#endif	 
   }
   for (; i < (3+d_x_vars)*total_edgenum; i++){
      is_int[i] = FALSE;
      v1 = edges[2*(i-(2+d_x_vars)*total_edgenum) + 1];
      ub[i] = flow_capacity - cnrp->demand[v1];
      matbeg[i]    = j;
      obj2[i]      = (double) costs[i- (2+d_x_vars)*total_edgenum];
#ifdef ADD_CAP_CUTS
      matval[j]    = 1.0;
      matval[j+1]  = -1.0;
      if (edges[2*(i - (2+d_x_vars)*total_edgenum)])
	 matval[j+2] = 1.0;
      matind[j++]  = (2+od_const)*vertnum - 1 + i -
	 (1+d_x_vars)*total_edgenum;
      matind[j++]  = (1+od_const)*vertnum + edges[2*(i -
			  (2+d_x_vars)*total_edgenum)+1] - 1;
      if (edges[2*(i - (2+d_x_vars)*total_edgenum)])
	 matind[j++] = (1+od_const)*vertnum + edges[2*(i -
			    (2+d_x_vars)*total_edgenum)] - 1;
#else
      matval[j]  = -1.0;
      if (edges[2*(i - (2+d_x_vars)*total_edgenum)])
	 matval[j+1] = 1.0;
      matind[j++]  = (1+od_const)*vertnum + edges[2*(i -
			  (2+d_x_vars)*total_edgenum)+1] - 1;
      if (edges[2*(i - (2+d_x_vars)*total_edgenum)])
	 matind[j++] = (1+od_const)*vertnum + edges[2*(i -
			    (2+d_x_vars)*total_edgenum)] - 1;
#endif
   }
   matbeg[i] = j;
   
   /* set the initial right hand side */
   if (od_const){
      /*degree constraints for the depot*/
#if 0
      rhs[0] = cnrp->numroutes;
      sense[0] = 'E';
      rhs[vertnum] = cnrp->numroutes;
      sense[vertnum] = 'E';
#else
      rhs[0] = 1.0;
      sense[0] = 'G';
      rhs[vertnum] = 1.0;
      sense[vertnum] = 'G';
#endif      
   }else if (prob_type == VRP || prob_type == TSP || prob_type == BPP){
      (rhs[0]) = 2*cnrp->numroutes;
      sense[0] = 'E';
   }else{
      /*cardinality constraint*/
      rhs[0] = vertnum - 1;
      sense[0] = 'E';
   }
   for (i = vertnum - 1; i > 0; i--){
      switch (prob_type){
       case VRP:
       case TSP:
       case BPP:
	 if (od_const){
	    rhs[i] = 1.0;
	    sense[i] = 'E';
	    rhs[i+vertnum] = 1.0;
	    sense[i+vertnum] = 'E';
	 }else{
	    rhs[i] = 2.0;
	    sense[i] = 'E';
	 }
	 break;
       case CSTP:
       case CTP:
	 rhs[i] = 1.0;
#ifdef DIRECTED_X_VARS
	 sense[i] = 'E';
#else
	 sense[i] = 'G';
#endif
	 break;
      }
#ifdef ADD_FLOW_VARS
      rhs[(1+od_const)*vertnum + i - 1] = cnrp->demand[i];
      sense[(1+od_const)*vertnum + i - 1] = 'E';
#endif
   }
#ifdef ADD_CAP_CUTS
   for (i = (2+od_const)*vertnum - 1;
	i < (2+od_const)*vertnum - 1 + 2*total_edgenum; i++){
      rhs[i] = 0.0;
      sense[i] = 'L';
   }
#endif
#ifdef ADD_X_CUTS
   for (i = (2+od_const)*vertnum-1+2*total_edgenum;
	i < (2+od_const)*vertnum-1+3*total_edgenum; i++){
      rhs[i] = 1;
      sense[i] = 'L';
   }
#endif
   
   sym_explicit_load_problem(env, n, m, matbeg, matind, matval, lb, ub, obj,
			     obj2, sense, rhs, NULL, FALSE);
   
}      
