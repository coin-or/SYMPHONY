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
/*                                                                           */
/*===========================================================================*/

/* $ID$ */

#ifndef _ROUNDING_H_
#define _ROUNDING_H_
/* rounding related */
#define  SYM_RND_FAIL  10

typedef struct ROUNDING_PROB {
   /* problem description */
   int 		n;
   int 		m;
   int          nz;
   int 		num_ints;
   int 		*int_vars;
   int          *is_int;
   double 	*xval;
   double       *obj;
   double       objval;
   double 	*rowActivity;
   double	*rowLower;
   double	*rowUpper;
   double 	*rowMin;
   double 	*rowMax;
   double       *rowAct;
   int		*rowStart;
   int    	*rowIndex;
   int		*rowLength;
   int		*rowIndices;
   double       *rowVal;
   int    	*colStart;
   int    	*colIndex;
   int    	*colLength;
   double	*colVal;
   double       *colAct;
   double       *lb;
   double       *ub;
   double 	lpetol;
   int          verbosity;
}rounding_problem;

int rnd_test(lp_prob *p);
int rnd_simple_backtrack(lp_prob *p);
int rnd_initialize_data(rounding_problem *rp, lp_prob *p);
int rnd_find_row_bounds(rounding_problem *rp);
int rnd_initialize_rp(rounding_problem *rp, lp_prob *p);
int rnd_delete_rp(rounding_problem *rp);
int rnd_is_constr_feas(rounding_problem *rp);
int rnd_grp_is_integral(int n, int *intvars, double *xval, double lpetol, int verbosity);
int rnd_lexicographic(rounding_problem *rp, lp_prob *p);
int rnd_find_feas_rounding(rounding_problem *rp, int grp_cnt, int *var_search_grp, int varnum, int recur_level);
int rnd_update_rp(rounding_problem *rp, int xind, double newval, double newlb, double newub);
int rnd_if_feas_impossible(int m, double *rowActivity, const double *rowLower, const double *rowUpper, double *rowMax, double *rowMin, double lpetol, int verbosity);
int sym_rnd_order(rounding_problem *rp, int **var_list, int *rn);
int rnd_free_rounding_problem(rounding_problem *rp);
int rnd_is_fully_integral(rounding_problem *rp);
int rnd_calculate_obj(rounding_problem *rp);
#endif
