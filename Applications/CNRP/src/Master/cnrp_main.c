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
#include "cnrp_types.h"

/*===========================================================================*\
 * This is main() for the CNRP
\*===========================================================================*/

int main(int argc, char **argv)
{
   cnrp_problem *cnrp = (cnrp_problem *) calloc(1, sizeof(cnrp_problem));

   sym_environment *env = sym_open_environment();

   sym_parse_command_line(env, argc, argv);

   sym_set_user_data(p, (void *) cnrp);

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
   
   if (*ub > 0){
      cnrp->cur_tour->cost = (int) (*ub);
   }else{
      cnrp->cur_tour->cost = MAXINT;
   }
   cnrp->cur_tour->numroutes = cnrp->numroutes;
   
   if (cnrp->par.use_small_graph == LOAD_SMALL_GRAPH){
      if (*ub <= 0 && cnrp->cur_tour->cost > 0)
	 *ub = (int)(cnrp->cur_tour->cost);
      cnrp->numroutes = cnrp->cur_tour->numroutes;
   }

#if 0
   if(cnrp->par.prob_tpye == BPP)
      *ub = 1;
#endif
   
   if (*ub > 0 && !(cnrp->par.prob_type == BPP)){
      printf("INITIAL UPPER BOUND: \t%i\n\n", (int)(*ub));
   }else if (!(cnrp->par.prob_type == BPP)){
      printf("INITIAL UPPER BOUND: \tNone\n\n");
   }else{
      printf("\n\n");
   }

   cnrp_load_problem(cnrp)
   
   sym_solve(env);

   sym_close_environment(env);

   return(0);
}

/*===========================================================================*\
 * This is the function that creates the formulation
\*===========================================================================*/

int cnrp_load_problem(cnrp_problem *cnrp)
{
   int *costs = cnrp->costs;
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
#ifdef FIND_NONDOMINATED_SOLUTIONS
   int edgenum = (mip->n - 1)/2;
#else
   int edgenum = (mip->n)/2;
#endif
#elif defined(ADD_FLOW_VARS)
#ifdef DIRECTED_X_VARS
#ifdef FIND_NONDOMINATED_SOLUTIONS
   int edgenum = (mip->n - 1)/4;
#else
   int edgenum = (mip->n)/4;
#endif

   flow_capacity = cnrp->capacity;
#else
#ifdef FIND_NONDOMINATED_SOLUTIONS
   int edgenum = (mip->n-1)/3;
#else
   int edgenum = (mip->n)/3;
#endif

   if (cnrp->par.prob_type == CSTP || cnrp->par.prob_type == CTP)
      flow_capacity = cnrp->capacity;
   else
      flow_capacity = cnrp->capacity/2;
#endif
#else
#ifdef FIND_NONDOMINATED_SOLUTIONS
   int edgenum = mip->n - 1;
#else
   int edgenum = mip->n;
#endif
#endif
   
   /* set up the inital LP data */

   /*Estimate the number of nonzeros*/
#ifdef ADD_CAP_CUTS
   mip->nz = 12*edgenum;
#elif defined(ADD_FLOW_VARS)
   mip->nz = 8*edgenum;
#else
   mip->nz = 3*edgenum;
#endif
#ifdef ADD_X_CUTS
   mip->nz += 2*edgenum;
#endif
#ifdef FIND_NONDOMINATED_SOLUTIONS
   mip->nz += mip->n + 1;
#endif
   *maxm = MAX(100, 3 * mip->m);
#ifdef ADD_FLOW_VARS
   *maxn = 3*total_edgenum;
#else
   *maxn = total_edgenum;
#endif
#ifdef DIRECTED_X_VARS
   *maxn += total_edgenum;
#endif
   *maxnz = mip->nz + ((*maxm) * (*maxn) / 10);

   /* Allocate the arrays. These are owned by SYMPHONY after returning. */
   mip->matbeg  = (int *) malloc((mip->n + 1) * ISIZE);
   mip->matind  = (int *) malloc((mip->nz) * ISIZE);
   mip->matval  = (double *) malloc((mip->nz) * DSIZE);
   mip->obj     = (double *) malloc(mip->n * DSIZE);
   mip->ub      = (double *) malloc(mip->n * DSIZE);
   mip->lb      = (double *) calloc(mip->n, DSIZE); /* zero lower bounds */
   mip->rhs     = (double *) malloc(mip->m * DSIZE);
   mip->sense   = (char *) malloc(mip->m * CSIZE);
   mip->rngval  = (double *) calloc(mip->m, DSIZE);
   mip->is_int  = (char *) calloc(mip->n, CSIZE);

#ifdef DIRECTED_X_VARS
   /*whether or not we will have out-degree constraints*/
   od_const = (prob_type == VRP || prob_type == TSP || prob_type == BPP);
   d_x_vars = TRUE;
#endif
   
   for (i = 0, j = 0; i < mip->n; i++){
#ifdef FIND_NONDOMINATED_SOLUTIONS
      if (indices[i] == mip->n - 1){
	 mip->is_int[i]    = FALSE;
	 mip->ub[i]        = MAXINT;
	 mip->matbeg[i]    = j;
	 mip->obj[i]       = 1.0;
	 mip->matval[j]    = -1.0;
	 mip->matind[j++]  = basecutnum;
	 mip->matval[j]    = -1.0;
	 mip->matind[j++]  = basecutnum + 1;
	 continue;
      }
#endif
      if (indices[i] < total_edgenum){
	 mip->is_int[i]    = TRUE;
	 mip->ub[i]        = 1.0;
	 mip->matbeg[i]    = j;
#ifdef FIND_NONDOMINATED_SOLUTIONS
	 mip->obj[i]       = cnrp->par.rho*((double) costs[indices[i]]);
	 mip->matval[j]    = cnrp->par.gamma*((double) costs[indices[i]]);
	 mip->matind[j++]  = basecutnum;
#else
	 mip->obj[i]       = cnrp->par.gamma*((double) costs[indices[i]]);
#endif
	 if (prob_type == CSTP || prob_type == CTP){
	    /*cardinality constraint*/
	    mip->matind[j] = 0;
	    mip->matval[j++] = 1.0;
	 }
	 /*in-degree constraint*/
	 mip->matval[j]    = 1.0;
	 mip->matind[j++]  = edges[2*indices[i]+1];
#ifdef DIRECTED_X_VARS
	 /*out-degree constraint*/
	 if (od_const){
	    mip->matval[j]   = 1.0;
	    mip->matind[j++] = vertnum + edges[2*indices[i]];
	 }
#else
	 if (prob_type == VRP || prob_type == TSP ||
	     prob_type == BPP || edges[2*indices[i]]){
	    mip->matval[j]   = 1.0;
	    mip->matind[j++] = edges[2*indices[i]];
	 }
#endif	 
#ifdef ADD_CAP_CUTS
	 v0 = edges[2*indices[i]];
	 mip->matval[j]    = -flow_capacity + (v0 ? cnrp->demand[v0] : 0);
	 mip->matind[j++]  = (2 + od_const)*vertnum - 1 + indices[i];
#ifndef DIRECTED_X_VARS
	 mip->matval[j]    = -flow_capacity +
	    cnrp->demand[edges[2*indices[i] + 1]];
	 mip->matind[j++]  = 2*cnrp->vertnum - 1 + total_edgenum +
	    indices[i];
#endif
#endif
#ifdef ADD_X_CUTS
	 mip->matval[j]    = 1.0;
	 mip->matind[j++]  = (2 + od_const)*vertnum-1 + 2*total_edgenum +
	    indices[i];
#endif
#ifdef DIRECTED_X_VARS
      }else if (indices[i] < 2*total_edgenum){
	 mip->is_int[i]    = TRUE;
	 mip->ub[i]        = 1.0;
	 mip->matbeg[i]    = j;
#ifdef FIND_NONDOMINATED_SOLUTIONS
	 mip->obj[i]       = cnrp->par.rho*((double) costs[indices[i] -
							  total_edgenum]);
	 mip->matval[j]    = cnrp->par.gamma*((double) costs[indices[i] -
							    total_edgenum]);
	 mip->matind[j++]  = basecutnum;
#else
	 mip->obj[i]       = cnrp->par.gamma*((double)costs[indices[i] -
							  total_edgenum]);
#endif
	 if (prob_type == CSTP || prob_type == CTP){
	    /*cardinality constraint*/
	    mip->matind[j] = 0;
	    mip->matval[j++] = 1.0;
	 }
	 /*in-degree constraint*/
	 if (od_const || edges[2*(indices[i] - total_edgenum)]){
	    mip->matval[j]   = 1.0;
	    mip->matind[j++] = edges[2*(indices[i] - total_edgenum)];
	 }
	 /*out-degree constraint*/
	 if (od_const){
	    mip->matval[j]    = 1.0;
	    mip->matind[j++]  = vertnum + edges[2*(indices[i] -
						   total_edgenum)+1];
	 }
#ifdef ADD_CAP_CUTS
	 mip->matval[j]    = -flow_capacity +
	    cnrp->demand[edges[2*(indices[i] - total_edgenum) + 1]];
	 mip->matind[j++]  = (2 + od_const)*vertnum - 1 + indices[i];
#endif
#ifdef ADD_X_CUTS
	 mip->matval[j]    = 1.0;
	 mip->matind[j++]  = (2 + od_const)*vertnum-1 + 2*total_edgenum +
	    indices[i] - total_edgenum;
#endif
#endif
      }else if (indices[i] < (2+d_x_vars)*total_edgenum){
	 mip->is_int[i] = FALSE;
	 v0 = edges[2*(indices[i]-(1+d_x_vars)*total_edgenum)];
	 mip->ub[i] = flow_capacity - (v0 ? cnrp->demand[v0] : 0);
	 mip->matbeg[i]    = j;
#ifdef FIND_NONDOMINATED_SOLUTIONS
	 mip->obj[i]       = cnrp->par.rho*((double) costs[indices[i] -
					(1+d_x_vars)*total_edgenum]);
	 mip->matval[j]    = cnrp->par.tau*((double) costs[indices[i]-
					(1+d_x_vars)*total_edgenum]);
	 mip->matind[j++]  = basecutnum + 1;
#else
	 mip->obj[i]       =
	    cnrp->par.tau*((double) costs[indices[i]-
					(1+d_x_vars)*total_edgenum]);
#endif
#ifdef ADD_CAP_CUTS
	 mip->matval[j]    = 1.0;
	 mip->matval[j+1]  = 1.0;
	 if (edges[2*(indices[i]-(1+d_x_vars)*total_edgenum)])
	    mip->matval[j+2] = -1.0;
	 mip->matind[j++]  = (2 + od_const)*vertnum - 1 + indices[i] -
	    (1+d_x_vars)*total_edgenum;
	 mip->matind[j++]  = (1+od_const)*vertnum + edges[2*(indices[i] -
				(1+d_x_vars)*total_edgenum) + 1] - 1;
	 if (edges[2*(indices[i] - (1+d_x_vars)*total_edgenum)])
	    mip->matind[j++] = (1+od_const)*vertnum + edges[2*(indices[i] -
				(1+d_x_vars)*total_edgenum)] - 1;
#else
	 mip->matval[j]  = 1.0;
	 if (edges[2*(indices[i]-(1+d_x_vars)*total_edgenum)])
	    mip->matval[j+1] = -1.0;
	 mip->matind[j++]  = (1+od_const)*vertnum + edges[2*(indices[i] -
				(1+d_x_vars)*total_edgenum) + 1] - 1;
	 if (edges[2*(indices[i] - (1+d_x_vars)*total_edgenum)])
	    mip->matind[j++] = (1+od_const)*vertnum + edges[2*(indices[i] -
				(1+d_x_vars)*total_edgenum)] - 1;
#endif	 
      }else{
	 mip->is_int[i] = FALSE;
	 v1 = edges[2*(indices[i]-(2+d_x_vars)*total_edgenum) + 1];
	 mip->ub[i] = flow_capacity - cnrp->demand[v1];
	 mip->matbeg[i]    = j;
#ifdef FIND_NONDOMINATED_SOLUTIONS
	 mip->obj[i]       = cnrp->par.rho*((double) costs[indices[i] -
					(2+d_x_vars)*total_edgenum]);
	 mip->matval[j]    = cnrp->par.tau*((double) costs[indices[i]-
					(2+d_x_vars)*total_edgenum]);
	 mip->matind[j++]  = basecutnum + 1;
#else
	 mip->obj[i]       =
	    cnrp->par.tau*((double) costs[indices[i]-
					(2+d_x_vars)*total_edgenum]);
#endif
#ifdef ADD_CAP_CUTS
	 mip->matval[j]    = 1.0;
	 mip->matval[j+1]  = -1.0;
	 if (edges[2*(indices[i] - (2+d_x_vars)*total_edgenum)])
	    mip->matval[j+2] = 1.0;
	 mip->matind[j++]  = (2+od_const)*vertnum - 1 + indices[i] -
	    (1+d_x_vars)*total_edgenum;
	 mip->matind[j++]  = (1+od_const)*vertnum + edges[2*(indices[i] -
				(2+d_x_vars)*total_edgenum)+1] - 1;
	 if (edges[2*(indices[i] - (2+d_x_vars)*total_edgenum)])
	    mip->matind[j++] = (1+od_const)*vertnum + edges[2*(indices[i] -
				(2+d_x_vars)*total_edgenum)] - 1;
#else
	 mip->matval[j]  = -1.0;
	 if (edges[2*(indices[i] - (2+d_x_vars)*total_edgenum)])
	    mip->matval[j+1] = 1.0;
	 mip->matind[j++]  = (1+od_const)*vertnum + edges[2*(indices[i] -
				(2+d_x_vars)*total_edgenum)+1] - 1;
	 if (edges[2*(indices[i] - (2+d_x_vars)*total_edgenum)])
	    mip->matind[j++] = (1+od_const)*vertnum + edges[2*(indices[i] -
				(2+d_x_vars)*total_edgenum)] - 1;
#endif
      }
   }
   mip->matbeg[i] = j;
   
   /* set the initial right hand side */
   if (od_const){
      /*degree constraints for the depot*/
#if 0
      mip->rhs[0] = cnrp->numroutes;
      mip->sense[0] = 'E';
      mip->rhs[vertnum] = cnrp->numroutes;
      mip->sense[vertnum] = 'E';
#else
      mip->rhs[0] = 1.0;
      mip->sense[0] = 'G';
      mip->rhs[vertnum] = 1.0;
      mip->sense[vertnum] = 'G';
#endif      
   }else if (prob_type == VRP || prob_type == TSP || prob_type == BPP){
      (mip->rhs[0]) = 2*cnrp->numroutes;
      mip->sense[0] = 'E';
   }else{
      /*cardinality constraint*/
      mip->rhs[0] = vertnum - 1;
      mip->sense[0] = 'E';
   }
   for (i = vertnum - 1; i > 0; i--){
      switch (prob_type){
       case VRP:
       case TSP:
       case BPP:
	 if (od_const){
	    mip->rhs[i] = 1.0;
	    mip->sense[i] = 'E';
	    mip->rhs[i+vertnum] = 1.0;
	    mip->sense[i+vertnum] = 'E';
	 }else{
	    mip->rhs[i] = 2.0;
	    mip->sense[i] = 'E';
	 }
	 break;
       case CSTP:
       case CTP:
	 mip->rhs[i] = 1.0;
#ifdef DIRECTED_X_VARS
	 mip->sense[i] = 'E';
#else
	 mip->sense[i] = 'G';
#endif
	 break;
      }
#ifdef ADD_FLOW_VARS
      mip->rhs[(1+od_const)*vertnum + i - 1] = cnrp->demand[i];
      mip->sense[(1+od_const)*vertnum + i - 1] = 'E';
#endif
   }
#ifdef ADD_CAP_CUTS
   for (i = (2+od_const)*vertnum - 1;
	i < (2+od_const)*vertnum - 1 + 2*total_edgenum; i++){
      mip->rhs[i] = 0.0;
      mip->sense[i] = 'L';
   }
#endif
#ifdef ADD_X_CUTS
   for (i = (2+od_const)*vertnum-1+2*total_edgenum;
	i < (2+od_const)*vertnum-1+3*total_edgenum; i++){
      mip->rhs[i] = 1;
      mip->sense[i] = 'L';
   }
#endif
#ifdef FIND_NONDOMINATED_SOLUTIONS
   mip->rhs[basecutnum] = cnrp->par.gamma*cnrp->utopia_fixed;   
   mip->sense[basecutnum] = 'L';
   mip->rhs[basecutnum+1] = cnrp->par.tau*cnrp->utopia_variable;   
   mip->sense[basecutnum+1] = 'L';
#endif
   return(USER_SUCCESS);
}      
