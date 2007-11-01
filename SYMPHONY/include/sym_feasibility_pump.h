/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/* by asm4                                                                   */
/*===========================================================================*/

#ifdef PRIMAL_HEURISTICS
#ifndef _FEASIBILITY_PUMP_H
#define _FEASIBILITY_PUMP_H
/* feasibility pump */
typedef struct FP_VARS {
   char          isBin;
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

int feasibility_pump (lp_prob *p, double &solution_value, double *betterSolution);
int fp_round (double *x_lp, double *x_ip, var_desc **vars, const int n);
int fp_is_feasible (LPdata *lp_data, double *x, int *Rmatbeg,
		    int *Rmatind, double *Rmatval, double *Rlower,
		    double *Rupper, FPdata *fp_data );
int fp_initialize_lp_solver(LPdata *lp_data, LPdata *new_lp_data, FPdata *fp_data);
int fp_solve_lp(double *x_lp, double *x_ip, int flip_rand, double T, LPdata *lp_data, int* indexList, FPdata *fp_data) ;
int fp_get_mip_desc(LPdata *lp_data, LPdata *newdata);
int fp_should_call_fp(lp_prob *p, int branching);
#endif
#endif
