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

#ifndef _CUT_GEN_USER_H
#define _CUT_GEN_USER_H

/* system include files */
#include <stdio.h>

/* SYMPHONY include files */
#include "BB_types.h"
#include "proto.h"

/* CNRP include files */
#include "network.h"
#include "cnrp_cg_params.h"

typedef struct CG_VRP_SPEC{
   cg_user_params par;
   int           dg_id;   /*contains the tid of the graphics window*/
   int           vertnum;  /*the number of nodes in the problem,
			      including the depot                */
   int          *demand;   /* alist of the customer demands*/
   int           capacity; /*the capacity of the trucks*/
   int           numroutes;/*contains the number of routes that the problem
			      is to be solved with. can be prespecified.  */
   int          *edges;    /*contains a list of the edges in the current
			      subproblem*/
   network      *n;
   int           orig_edgenum;
   int          *cost;
   int          *ref;      /* the last five  are for the shrinking routines; */
   char         *in_set;   /* They are here to optimize/speed up things */
   int          *new_demand;
   double       *cut_val;
   char         *cut_list;

/*__BEGIN_EXPERIMENTAL_SECTION__*/
   int          *dec_data;
   int           last_decomp_index;
   double        last_objval;
   FILE         *decomp_res; 
   /* the next four arrays pertain to storing no-columns cuts - kind of an
      auxiliary  cutpool*/ 
   int         **data;
   char        **indicators;
   int          *ones;
   int          *size;
   int           num_nocolscuts;
/*___END_EXPERIMENTAL_SECTION___*/

#ifdef CHECK_CUT_VALIDITY
   int           feas_sol_size;
   int          *feas_sol;
#endif
}cg_vrp_spec;

/*===========================================================================*/
/*========================= Other user subroutines =========================*/
/*===========================================================================*/

int check_connectivity PROTO((network *n, double etol, int capacity,
			      int numroutes, char mult));

int check_flow_connectivity PROTO((network *n, double etol, int capacity,
				   int numroutes, char mult));

/*===========================================================================*/
/*=============================== shrink.c ==================================*/
/*===========================================================================*/

int reduce_graph PROTO((network *n, double etol, int *demand, int capacity,
			 int mult, cut_data *new_cut));
int greedy_shrinking1 PROTO((network *n, double truck_cap, double etol,
			     int max_num_cuts, cut_data *new_cut,
			     int *compnodes, int *compmembers, int compnum,
			     char *in_set, double *cut_val,int *ref,
			     char *cut_list, int *demand, int mult));
int greedy_shrinking6 PROTO((network *n, double truck_cap,
			     double etol, cut_data *new_cut,
			     int *compnodes,
			     int *compmembers, int compnum, char *in_set,
			     double *cut_val,int *ref, char *cut_list,
			     int max_num_cuts, int *demand, int trial_num,
			     double prob, int mult));
int greedy_shrinking1_one PROTO((network *n, double truck_cap,
				 double etol, int max_num_cuts,
				 cut_data *new_cut, char *in_set,
				 double *cut_val, char *cut_list,
				 int num_routes, int *demand, int mult));
int greedy_shrinking6_one PROTO((network *n, double truck_cap,
				 double etol, cut_data *new_cut,
				 char *in_set, double *cut_val, int num_routes,
				 char *cut_list, int max_num_cuts,
				 int *demand,int trial_num, double prob,
				 int mult));
int greedy_shrinking2_one PROTO((network *n, double truck_cap,
				 double etol, cut_data *new_cut,
				 char *in_set, double *cut_val, int num_routes,
				 int *demand, int mult));

/*===========================================================================*/
/*============================ biconnected.c ================================*/
/*===========================================================================*/

void depth_first_search PROTO((vertex *v, int *count1, int *count2));
int biconnected PROTO((network *n, int *compnodes, int
			   *compdemands, double *compcuts));
void compute_comp_nums PROTO((vertex *v, int parent_comp, int *num_comps,
		       char parent_is_art_point));
int tsp_cuts PROTO((network *n, int verbosity, char tsp_prob, int which_cuts));

#endif
