/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*
 * TODO:
 * change ifdef _FEASI.. to ifdef _HEURISTICS
 */

#ifndef _PRIMAL_HEURISTICS_H
#define _PRIMAL_HEURISTICS_H
#include "sym_lp_solver.h"
#include "sym_lp.h"
#include "sym_types.h"

/* feasibility pump */
typedef struct FP_VARS {
   char          isBin;
   char          isInt;
   int           xplus;
   int           xminus;
}FPvars;

typedef struct FP_DATA {
   FPvars        **fp_vars;	/* an array of fp_vars */
   int           n0;		/* no. of vars in orignial lp */
   int           m0;
   int           numNonBinInts;
   int           numInts;
}FPdata;

/*  solution pool */
int sp_add_solution PROTO((lp_prob *p, int cnt, int *indices, double *values, double obj_value, int bc_index));
int sp_delete_solution PROTO((sp_desc *sp, int position));
int sp_is_solution_in_sp PROTO((lp_prob *p, int cnt, int *indices, double *values, double obj_value));
int sp_initialize(tm_prob *tm);
int sp_free_sp(sp_desc *sp);

/* feasibility pump */
int feasibility_pump (lp_prob *p, char *found_better_solution, double &solution_value, double *betterSolution);
int fp_round (double *x_lp, double *x_ip, FPvars **vars, const int n);
int fp_is_feasible (LPdata *lp_data, double *x, const CoinPackedMatrix *matrix,  const double *r_low, const double *r_up, FPdata *fp_data );
int fp_initialize_lp_solver(LPdata *lp_data, LPdata *new_lp_data, FPdata *fp_data);
int fp_solve_lp(double *x_lp, double *x_ip, int flip_rand, double T, LPdata *lp_data, int* indexList, FPdata *fp_data) ;
int fp_get_mip_desc(LPdata *lp_data, LPdata *newdata);
int fp_should_call_fp(lp_prob *p, int branching);
int fp_add_obj_row(LPdata *new_lp_data, int n, const double *obj, double rhs);
#endif
