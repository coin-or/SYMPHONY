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

#ifndef _API_H
#define _API_H

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
#ifdef COMPILE_IN_CP
#include "cp.h"
#endif
#endif
#include "lp_solver.h" /* For MIPdesc */

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
   char       infile[MAX_FILE_NAME_LENGTH +1]; /* For MPS file name
						  or GNUMP modelfile */
   char       datafile[MAX_FILE_NAME_LENGTH +1]; /* GNUMP datafile */
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

   MIPdesc   *mip;     /* For holding the description when read in from MPS */
   char       probname[81];

   base_desc *base;
   node_desc *root;

   tm_stat    stat;
}problem;

/*===========================================================================*/
/*===================== Interface functions (master.c) ======================*/
/*===========================================================================*/

problem *sym_open_environment PROTO((void));
void sym_load_problem PROTO((problem *p, int argc, char **argv));
void sym_solve PROTO((problem *p));
void sym_close_environment PROTO((problem *p));

#endif
