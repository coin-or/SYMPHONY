/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
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

#ifndef SYM_PREP__H
#define SYM_PREP__H

#include "symphony.h"
#include "sym_proccomm.h" 
#include "sym_master.h" 
#include "sym_constants.h"
#include "sym_prep_params.h"
#include <list>

/* return codes of the basic preprocessor */

#define PREP_UNMODIFIED 0	/* preprocessor exited without modifying the
				   MIP in any way*/
#define PREP_MODIFIED 1		/* preprocessor modified the MIP in some way */
#define PREP_INFEAS 2		/* preprocessor found the MIP infeasible */
#define PREP_SOLVED 3		/* preprocessor found the MIP unbounded */
#define PREP_UNBOUNDED 4	/* preprocessor found the MIP unbounded */

typedef struct ROWPACKEDARRAY
{
   int    *matbeg;
   int    *matind;
   double *matval;
   char   *isdeleted;
   int     rows_deleted;        /* Number of rows that have been marked
				   as deleted */
} rowpackedarray;

typedef struct LHSPARAMS
{
   /* all arrays of size m */
   /* LHS is a structure which contains the bounds on the constraints */
   int *ub_is_infinite;         /* is 1 if there is no finite upperbound on a
				   constraint's value */
   int *lb_is_infinite;         /* is 1 if there is no finite lowerbound on a
				   constraint's value */
   int *ub_inf_var;             /* is -1 if more than one variables cause
				   infinite upperbound, otherwise has the index
				   of the variable which causes the same. */
   int *lb_inf_var;             /* is -1 if more than one variables cause
				   infinite lowerbound, otherwise has the index
				   of the variable which causes the same. */
   double *ubound;              /*  contains lower bounds of the value of a
				    constraint.*/
   double *lbound;              /*  contains lower bounds of the value of a
				    constraint.*/
}lhsparams;

typedef struct PREP_STATS
{
   int rows_deleted;
   int vars_fixed;
   int coeffs_changed;
   int bounds_tightened;
}prep_stats;

class implication
{
 public:
  implication();
  int fixed_col;		/* variable which is fixed */
  double fixed_val;		/* value at which the variable is fixed */
  implication(int col_num, double value):
    fixed_col(col_num),
    fixed_val(value)
    {
    };
  ~implication()
    {
    };
  /* stores one implication of fixing certain var to a certain value */
};

typedef std::list<implication> impList;

int prep_basic(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
	       prep_stats *stats, prep_params *prep_par);
int prep_light(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
	       prep_stats *stats, prep_params *prep_par);
int prep_advanced(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		  impList *implistL, impList *implistU, prep_stats *stats,
		  prep_params *prep_par);
int prep_integerize_bounds(MIPdesc *P, int & bounds_integerized,
			   int verbosity);
int prep_create_row_rep(MIPdesc *P, rowpackedarray *row_P);
int prep_find_lhs_params(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs);
int prep_tighten_bounds(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			const int max_iterations, int & bounds_tightened,
			int & vars_fixed, int verbosity);
int prep_fix_imp_vars(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		      int &tighten_bounds_out);
int prep_find_imp_sibling(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			  const int col_num, int verbosity);
int prep_check_feas(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs);
int prep_check_redundancy(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			  int verbosity);
int prep_fix_variables(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		       int & variables_fixed, int verbosity);
int prep_BinImproveCoeffs(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			  int & coeffs_changed, int verbosity);
int prep_delete_row(int row_num, MIPdesc *P, rowpackedarray *row_P,
		    lhsparams *lhs);
int prep_purge_del_rows(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			int & rows_purged);
int prep_purge_del_rows2(sym_environment *env, MIPdesc *P, rowpackedarray *row_P, int & rows_purged);
int prep_display_row_mat(MIPdesc *P , rowpackedarray *row_P);
int prep_declare_solved(MIPdesc *P, int verbosity);
int prep_update_const_bounds(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			     char bnd_sense, int var_num, double new_bound);
int prep_display_mip(MIPdesc *P);
int prep_changeCoeff(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		     int col_num, int row_num, const double newCoeff,
		     double newRhs);
int prep_isBinary(MIPdesc *P, int col_num);
void prep_coeffChangeBounds(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			    int col_num, int row_num, double oldCoeff,
			    double newCoeff);
int prep_free_row_structs(lhsparams *lhs, rowpackedarray *row_P);
int prep_find_logical_imp(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			  impList *implistL, impList *implistU,
			  prep_stats *stats, prep_params *prep_par);
int prep_reset_mip_data (MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			 MIPdesc *fP, rowpackedarray *row_fP,
			 lhsparams *f_lhs);
int prep_delete_structs (MIPdesc *fP, rowpackedarray *row_fP,
			 lhsparams *f_lhs, prep_stats *stats);
int prep_disp_stats (prep_stats *stats);
int prep_add_imp(MIPdesc *P, MIPdesc *fP, const int col_num, impList *il,
		 int verbosity);
int prep_list_imp(MIPdesc *P, impList *il);
int prep_apply_imp(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		   int col_num, double col_val, impList *il,
		   prep_stats *stats);
bool prep_isImplied(const int col_num, const int imp_col, impList *il,
		    double &imp_val);
#endif
