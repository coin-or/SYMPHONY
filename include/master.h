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

#ifndef _MASTER_H
#define _MASTER_H

#define COMPILING_FOR_MASTER

#include "proto.h"
#include "BB_types.h"
#include "tm_params.h"
#include "cp_params.h"
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
#include "sp_params.h"
#endif
/*___END_EXPERIMENTAL_SECTION___*/
#include "cg_params.h"
#include "lp_params.h"
#include "dg_params.h"
#include "master_u.h"
#ifdef COMPILE_IN_TM
#include "tm.h"
#ifdef COMPILE_IN_LP
#include "lp.h"
#endif
#endif

/*===========================================================================*\
 * Here we keep track of the computation time for each of the various
 * parts of the computation
\*===========================================================================*/

typedef struct PROB_TIMES{
   double     readtime;    /* time spent reading in the problem*/
   node_times bc_time;
   double     ub_overhead; /* overhead time used doing the upper bounding */
   double     ub_heurtime; /* actual comp time doing the upper bounding */
   double     lb_overhead; /* overhead time doing the lower bounding */
   double     lb_heurtime; /* actual comp time doing the lower bounding */
}prob_times;

/*===========================================================================*\
 * The params structure contains all of the user-specified parameters
 * to be read in from the parameter file. See the README file for an
 * explanation of the parameters
\*===========================================================================*/

typedef struct PARAMS{
   int        warm_start;
   int        verbosity;
   char       param_file[MAX_FILE_NAME_LENGTH +1];
   int        random_seed;
   cp_params  cp_par;
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   sp_params  sp_par;
#endif
/*___END_EXPERIMENTAL_SECTION___*/
   cg_params  cg_par;
   lp_params  lp_par;
   tm_params  tm_par;
   dg_params  dg_par;
   char       tm_exe[MAX_FILE_NAME_LENGTH +1];
   char       dg_exe[MAX_FILE_NAME_LENGTH +1];
   int        tm_debug;
   int        dg_debug;
   int        tm_machine_set;
   char       tm_machine[MACH_NAME_LENGTH +1];
   int        dg_machine_set;
   char       dg_machine[MACH_NAME_LENGTH +1];
   int        pvm_trace;
   int        do_branch_and_cut;
   int        do_draw_graph;
}params;

/*===========================================================================*\
 * The problem data structure contains the data for a problem instance, as
 * well as some of the tours that have been generated.
\*===========================================================================*/

typedef struct PROBLEM{
   void      *user;
#ifdef COMPILE_IN_TM
   tm_prob   *tm;
#endif
   int        tm_tid;
   int        dg_tid;
   params     par;         /* problem parameters */
   prob_times comp_times;  /* keeps track of the computation times for the
			      problem */
   char       has_ub;
   double     ub;
   lp_sol     best_sol;
   char       has_ub_estimate;
   double     ub_estimate;
   double     lb;

   tm_stat    stat;
}problem;

/*===========================================================================*/
/*==================== Master basic functions (master.c) ====================*/
/*===========================================================================*/

problem *get_problem_ptr PROTO((void));
void print_statistics PROTO((node_times *tim, tm_stat *stat, double ub,
			     double lb, double initial_time,
			     double start_time));

/*===========================================================================*/
/*=================== Master I/O functions (readparams.c) ===================*/
/*===========================================================================*/

void bc_readparams PROTO((problem *p, int argc, char **argv));
void read_string PROTO((char *target, char *line, int maxlen));

/*===========================================================================*/
/*=============== Master wrapper functions (master_wrapper.c) ===============*/
/*===========================================================================*/

void initialize_u PROTO((problem *p));
void free_master_u PROTO((problem *p));
void readparams_u PROTO((problem *p, int argc, char **argv));
void io_u PROTO((problem *p));
void init_draw_graph_u PROTO((problem *p));
void start_heurs_u PROTO((problem *p));
void display_solution_u PROTO((problem *p, int thread_num));
node_desc *create_root_u PROTO((problem *p));
base_desc *set_base_u PROTO((problem *p));
void receive_feasible_solution_u PROTO((problem *p, int msgtag));
void send_lp_data_u PROTO((problem *p, int sender, base_desc *base));
void send_cg_data_u PROTO((problem *p, int sender));
void send_cp_data_u PROTO((problem *p, int sender));
void send_sp_data_u PROTO((problem *p, int sender));
void process_own_messages_u PROTO((problem *p, int msgtag));

#endif
