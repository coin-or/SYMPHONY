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
#ifndef WIN32
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <sys/types.h>

/* SYMPHONY include files */
#include "BB_constants.h"
#include "proccomm.h"
#include "qsortucb.h"

/* CNRP include files */
#include "cnrp_const.h"
#include "cnrp_master_functions.h"
#include "cnrp_macros.h"

/*===========================================================================*/

int is_same_edge(const void *ed0, const void *ed1)
{
   return((INDEX(((edge_data *)ed0)->v0, ((edge_data *)ed0)->v1)) -
	  (INDEX(((edge_data *)ed1)->v0, ((edge_data *)ed1)->v1)));
}

/*===========================================================================*/

/*__BEGIN_EXPERIMENTAL_SECTION__*/
#if 0
/* This comparison function can be used to sort lexicographically */
/*===========================================================================*/

int is_same_edge(const void *ed0, const void *ed1)
{
   return(((edge_data *)ed0)->v0 - ((edge_data *)ed1)->v0 ?
	  ((edge_data *)ed0)->v0 - ((edge_data *)ed1)->v0 :
	  ((edge_data *)ed0)->v1 - ((edge_data *)ed1)->v1);
}

/*===========================================================================*/
#endif
/*__END_EXPERIMENTAL_SECTION__*/
void delete_dup_edges(small_graph *g)
{
   edge_data *ed0, *ed1;
   int pos;
   
   qsort((char *)g->edges, g->edgenum, sizeof(edge_data), is_same_edge);
   for (pos=0, ed0=ed1=g->edges ; pos < g->edgenum; pos++, ed1++){
      if ( memcmp((char *)ed0, (char *)ed1, 2*sizeof(int)) ){
	 ed0++;
	 (void)memcpy((char *)ed0, (char *)ed1, sizeof(edge_data));
      }
   }
   pos = ((int)ed0 - (int)g->edges)/sizeof(edge_data) + 1;
   g->allocated_edgenum -= g->edgenum - pos;
   g->edges = (edge_data *) realloc
      ((char *)(g->edges), g->allocated_edgenum * sizeof(edge_data));
   g->edgenum = pos;
}

/*===========================================================================*/

int *create_edge_list(cnrp_problem *cnrp, int *varnum, char which_edges)
{
   int i, j, k;
   int zero_varnum, edgenum, new_ind;
   int *zero_vars, *uind = NULL;
   int total_edgenum = cnrp->vertnum*(cnrp->vertnum-1)/2;
#ifdef DIRECTED_X_VARS
   char d_x_vars = TRUE;
#else
   char d_x_vars = FALSE;
#endif

   /*DIFF: This routine has to be modified to include the flow variables*/

   switch(which_edges){
    case CHEAP_EDGES:

      cnrp->zero_vars = zero_vars = (int *) calloc(total_edgenum, sizeof(int));
      
      /*first determine which variables can be fixed to zero permanently*/
      for (zero_varnum=0, i=2; i<cnrp->vertnum; i++){
	 for (j=1; j<i; j++){
	    if (cnrp->demand[i] + cnrp->demand[j] > cnrp->capacity){
	       zero_vars[zero_varnum++] = INDEX(i,j);
	    }
	 }
      }
      
      edgenum = cnrp->par.add_all_edges ?
	 cnrp->vertnum*(cnrp->vertnum-1)/2 : cnrp->g->edgenum;
      
      /*First, we construct the index lists*/
#ifdef ADD_FLOW_VARS
      uind = (int *) malloc((3+d_x_vars) * edgenum * ISIZE);
#else
      uind = (int *) malloc((1+d_x_vars) * edgenum * ISIZE);
#endif
      
      *varnum = 0;
      switch(cnrp->par.add_all_edges){
       case FALSE:
	 for (i = 0, j = 0; i<edgenum && j<zero_varnum; i++){
	    new_ind = INDEX(cnrp->g->edges[i].v0, cnrp->g->edges[i].v1);
	    if (new_ind < zero_vars[j]){
	       uind[(*varnum)++] = new_ind;                 /*edge var*/
#ifdef DIRECTED_X_VARS
	       uind[(*varnum)++] = total_edgenum + new_ind; /*edge var*/
#endif
#ifdef ADD_FLOW_VARS
	       /*flow var*/
	       uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + new_ind;
	       /*flow var*/
	       uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + new_ind;
#endif
	    }else{
	       while (j < zero_varnum && new_ind > zero_vars[j])
		  j++;
	       if (j == zero_varnum){
		  uind[(*varnum)++] = new_ind;                   /*edge var*/
#ifdef DIRECTED_X_VARS
		  uind[(*varnum)++] = total_edgenum + new_ind;   /*edge var*/
#endif
#ifdef ADD_FLOW_VARS
		  /*flow var*/
		  uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + new_ind;
		  /*flow var*/
		  uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + new_ind;
#endif
	       }else if (new_ind < zero_vars[j]){
		  uind[(*varnum)++] = new_ind;                   /*edge var*/
#ifdef DIRECTED_X_VARS
		  uind[(*varnum)++] = total_edgenum + new_ind;   /*edge var*/
#endif
#ifdef ADD_FLOW_VARS
		  /*flow var*/
		  uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + new_ind;
		  /*flow var*/
		  uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + new_ind;
#endif
	       }else
		  j++;
	    }
	 }
	 /*Now we have exhausted all the zero edges*/
	 for (; i<edgenum; i++){
	    uind[(*varnum)++] =
	       INDEX(cnrp->g->edges[i].v0, cnrp->g->edges[i].v1);
#ifdef DIRECTED_X_VARS
	    uind[(*varnum)++] =
	       total_edgenum +
	       INDEX(cnrp->g->edges[i].v0, cnrp->g->edges[i].v1);
#endif
#ifdef ADD_FLOW_VARS
	    uind[(*varnum)++] =
	       (1+d_x_vars)*total_edgenum+INDEX(cnrp->g->edges[i].v0,
						cnrp->g->edges[i].v1);
	    uind[(*varnum)++] =
	       (2+d_x_vars)*total_edgenum+INDEX(cnrp->g->edges[i].v0,
						cnrp->g->edges[i].v1);
#endif
	 }
	 break;
       case TRUE:
	 for (i = 0, j = 0; j<zero_varnum; i++){
	    if (zero_vars[j] == i){
	       j++;
	       continue;
	    }
	    uind[(*varnum)++] = i;                    /*edge variable*/
#ifdef DIRECTED_X_VARS
	    uind[(*varnum)++] = total_edgenum + i;    /*edge variable*/
#endif
#ifdef ADD_FLOW_VARS
	    /*flow var*/
	    uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + i;
	    /*flow var*/
	    uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + i;
#endif
	 }/*Now, we have exhausted all the zero edges*/
	 for (; i < edgenum; i++){
	    uind[(*varnum)++] = i;                    /*edge variable*/
#ifdef DIRECTED_X_VARS
	    uind[(*varnum)++] = total_edgenum + i;    /*edge variable*/
#endif
#ifdef ADD_FLOW_VARS
	    /*flow var*/
	    uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + i;
	    /*flow var*/
	    uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + i;
#endif
	 }
	 break;
      }

      if (cnrp->par.verbosity > 0)
	 printf("Fixed %i edges in root creation\n\n", zero_varnum);
      
      cnrp->zero_varnum = zero_varnum;

      break;
      
    case REMAINING_EDGES:

      /*In this case, we are adding all variables at the root, but the small
	graph edges are base and the rest are extra*/

      zero_varnum = cnrp->zero_varnum;
      zero_vars = cnrp->zero_vars;
      edgenum = cnrp->g->edgenum;

#ifdef ADD_FLOW_VARS
      uind = (int *) malloc((3 + d_x_vars) *
			    (total_edgenum-edgenum+cnrp->vertnum -1)* ISIZE);
#else
      uind = (int *) malloc((1 + d_x_vars) *
			    (total_edgenum-edgenum+cnrp->vertnum-1) * ISIZE);
#endif
      *varnum = 0;
      for (i = 0, j = 0, k = 0; i < edgenum; i++, k++){
	 /*In this loop, we check each edge to se if it is in the small
	   graph and whether it is a zero edge*/
	 new_ind = INDEX(cnrp->g->edges[i].v0, cnrp->g->edges[i].v1);
	 for (; k < new_ind; k++){
	    if ((j < zero_varnum && k < zero_vars[j]) || j >= zero_varnum){
	       uind[(*varnum)++] = k;                    /*edge variable*/
#ifdef DIRECTED_X_VARS
	    uind[(*varnum)++] = total_edgenum + k;       /*edge variable*/
#endif
#ifdef ADD_FLOW_VARS
	    /*flow var*/
	    uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + k;
	    /*flow var*/
	    uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + k;
#endif
	    }else{ /*curent edge is a zero edge so don't add it*/
	       j++;
	    }
	 }
	 /*k == new_ind here so we don't want to add that edge */
      }
      /*Now, we have exhausted the small graph so just add all non-zero
	edges*/
      for (; k < total_edgenum && j < zero_varnum; k++)
	 if (k < zero_vars[j]){
	    uind[(*varnum)++] = k;                    /*edge variable*/
#ifdef DIRECTED_X_VARS
	    uind[(*varnum)++] = total_edgenum + k;    /*edge variable*/
#endif
#ifdef ADD_FLOW_VARS
	    /*flow var*/
	    uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + k;
	    /*flow var*/
	    uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + k;
#endif
	 }else{
	    j++;
	 }
      /* Now, there are no more non-zero edges either so add the rest*/
      for (; k < total_edgenum; k++){
	 uind[(*varnum)++] = k;                    /*edge variable*/
#ifdef DIRECTED_X_VARS
	 uind[(*varnum)++] = total_edgenum + k;    /*edge variable*/
#endif
#ifdef ADD_FLOW_VARS
	 /*flow var*/
	 uind[(*varnum)++] = (1+d_x_vars)*total_edgenum + k;
	 /*flow var*/
	 uind[(*varnum)++] = (2+d_x_vars)*total_edgenum + k;
#endif
      }
      break;
   }
   qsortucb_i(uind, *varnum);

   return(uind);
}
