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

#ifndef _MASTER_H
#define _MASTER_H

#define COMPILING_FOR_MASTER

#include "BB_types.h"
#include "master_params.h"
#include "lp_solver.h"
#ifdef COMPILE_IN_TM
#include "tm.h"
#endif

/*===========================================================================*\
 * The problem data structure contains the data for a problem instance, as
 * well as some of the tours that have been generated.
\*===========================================================================*/

typedef struct PROBLEM{
   void            *user;
   int              tm_tid;
   int              dg_tid;
   params           par;         /* problem parameters */
   prob_times       comp_times;  /* keeps track of the computation times for
				    the problem */
   char             has_ub;
   double           ub;
   lp_sol           best_sol;
   char             has_mc_ub;
   double           mc_ub;
   double           obj[2];
   double           utopia[2];
   char             has_ub_estimate;
   double           ub_estimate;
   double           lb;

   MIPdesc         *mip; /* for holding the description when read in from MPS */

   char             probname[81];

   base_desc       *base;
   node_desc       *rootdesc;

   int              termcode;

#ifdef COMPILE_IN_TM
   tm_prob         *tm;
   warm_start_desc *warm_start;
#ifdef COMPILE_IN_CP
   cut_pool       **cp;
#endif
#endif
}problem;

/*===========================================================================*/
/*=================== Master I/O functions (readparams.c) ===================*/
/*===========================================================================*/

int parse_command_line PROTO((problem *p, int argc, char **argv));
void read_string PROTO((char *target, char *line, int maxlen));
void print_statistics PROTO((node_times *tim, problem_stat *stat, double ub,
			     double lb, double initial_time,
			     double start_time, double obj_offset,
			     char obj_sense, char has_ub));

/*===========================================================================*/
/*=============== Master wrapper functions (master_wrapper.c) ===============*/
/*===========================================================================*/

int initialize_u PROTO((problem *p));
int free_master_u PROTO((problem *p));
int readparams_u PROTO((problem *p, int argc, char **argv));
int io_u PROTO((problem *p));
int init_draw_graph_u PROTO((problem *p));
int start_heurs_u PROTO((problem *p));
int display_solution_u PROTO((problem *p, int thread_num));
int initialize_root_node_u PROTO((problem *p));
int receive_feasible_solution_u PROTO((problem *p, int msgtag));
int send_lp_data_u PROTO((problem *p, int sender));
int send_cg_data_u PROTO((problem *p, int sender));
int send_cp_data_u PROTO((problem *p, int sender));
int send_sp_data_u PROTO((problem *p, int sender));
int process_own_messages_u PROTO((problem *p, int msgtag));

/*===========================================================================*/
/*=================== Master helper functions (master_func.c) ===============*/
/*===========================================================================*/

int resolve_node PROTO((problem *p, bc_node * node));
void update_tree_bound PROTO((problem *p, bc_node *root, int change_type));

int copy_node PROTO((bc_node * n_to, bc_node *n_from));
int copy_tree PROTO((bc_node *root_to, bc_node *root_from));
int read_node PROTO((bc_node * node, FILE * f));
int read_tree PROTO((bc_node * root, FILE *f));
int write_node PROTO((bc_node *node, FILE*f));
int write_tree PROTO((bc_node *root, FILE *f));

int set_param PROTO((problem *p,  char *line));

warm_start_desc *create_copy_warm_start PROTO((warm_start_desc * ws));
MIPdesc *create_copy_mip_desc PROTO((MIPdesc *mip));
problem *create_copy_problem PROTO((problem *p));

double get_lb_for_new_rhs PROTO((bc_node *root, MIPdesc *mip, int cnt, 
				 int *ind, double *val));
int check_feasibility_new_rhs PROTO((bc_node * node, MIPdesc * mip, 
					int cnt, int *ind, double *val));
#endif
