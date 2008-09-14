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
#define PREP_NUMERIC_ERROR 5

#define SR_NO_UPDATES 0
#define SR_BOUNDS_UPDATED 1
#define SR_INFEAS 2

#define POS_VAL 0
#define ZERO_VAL 1
#define NEG_VAL 2

#define SR_MIN 0
#define SR_MAX 1

#define RND_FLOOR 0
#define RND_CEIL 1

#define FIX_NO_BOUND 0
#define FIX_BINARY 1
#define FIX_OTHER  2
#define FIX_FIXABLE 3
#define FIX_UB 4
#define FIX_LB 5

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
   int * isdeleted;             /* deleted constraints */
   int num_deleted; 
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
   int bounds_integerized;
}prep_stats;


typedef struct SRDESC{

   int prob_type;
   char sense;
   double rhs;
   
   int max_n;
   double *obj_max;
   double *matval_max;
   double *ratio_max;
   int *matind_max;
   //  int *ratio_type_max;
   double ub_offset;
   double rhs_max;
   double sum_c_max;
   double sum_a_max;

   char ub_updated; 
   double ub;

   int min_n;
   double *obj_min;
   double *matval_min;
   double *ratio_min;
   int *matind_min;
   //   int *ratio_type_min;
   double lb_offset;
   double rhs_min;
   double sum_c_min;
   double sum_a_min;

   char lb_updated; 
   double lb;

   /* for sorting purposes */
   int * fixed_ind; 
   int * tmp_ind;


   /* for variable fixing, bound tightening purpose*/
   double * lb_var; /* for binary, variable is fixed to obtain this value 
		       for others, obj coef is fixed */
   
   double * ub_var;  /* just useful for binary case, 
			variable is fixed for this case */   
}SRdesc;



#if 0
typedef struct SRDESC{
    
   int n;
   int rhs;
   int sense;
   
   double *obj;
   double *matval;
   double *ratios;
   int *matind;

   char no_upper; 
   char no_lower;

   double ub_offset;
   double lb_offset;
   
   double rhs_lb_offset;
   double rhs_ub_offset;  

   double ub;
   double lb;

}SRdesc;
#endif

#if 0
typedef struct SRRELAX
{
   char exists;
   int row_ind; 
   int constr_row_ind; /* -1 for aggregated rows*/

   int sub_n;
   double *obj;
   double *matval;
   int *matind; /* indices of colums appear in this problem */
   double *ratios;  

   /* filled if it is an aggregated row */
   double rhs;
   char sense;
     
   double lb;
   double ub;
   
   int lb_sol_size;
   int *lb_sol; 

   int ub_sol_size;
   int *ub_sol; 
}SRrlx;
#endif 

typedef struct PREPDesc
{
   MIPdesc * mip; 
   prep_stats stats; 
   prep_params params;

   /* trying single/aggr row relaxations to improve bounds*/
   int max_sr_cnt; 
   int max_aggr_cnt; 
   SRdesc *sr; /* for 'L', 'G' constraints */
   SRdesc *d_sr; /* additionally, for 'E' constraints */    
   /* for subproblems checking purposes */
   char *rows_checked;    

}PREPdesc;


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
    }
  ~implication()
    {
    }
  /* stores one implication of fixing certain var to a certain value */
};

typedef std::list<implication> impList;

int preprocess_mip(MIPdesc *mip, prep_params prep_par, char imply_changes, 
		   char keeptrack);

int prep_initialize_mipinfo(MIPdesc *mip,  prep_params prep_par);

int prep_fill_row_ordered(MIPdesc *mip);

int prep_integerize_bounds(PREPdesc *P, char keeptrack);

int prep_basic(PREPdesc *P, char keeptrack);

double prep_rnd_integral(double val, double etol, char rnd_type);

int prep_fix_variable(MIPdesc *mip, int col_ind, int row_ind, int a_loc, 
		      double etol);

int prep_check_redundancy(MIPdesc *mip, int row_ind,  double etol, 
			  char use_sr_bounds);

void prep_fixed_col_update_info(MIPdesc *mip, int col_ind, double fixed_bound,  
				      int fix_type);
void prep_deleted_row_update_info(MIPdesc *mip, int row_ind);


int prep_solve_sr_rlx(PREPdesc *P, int row_cnt, int *row_indices); 
void sr_initialize(SRdesc **sr, int n);
void sr_allocate(SRdesc **sr, int n);
int sr_solve_bounded_prob(SRdesc *sr, SRdesc *d_sr, int obj_ind, int row_ind, 
			  int *r_matbeg, int *r_matind, double *r_matval, 
			  COLinfo *cols, double *ub, double *lb, double etol);

int sr_add_new_col(SRdesc *sr, SRdesc *d_sr, double c_val, double a_val, 
		   int col_ind, char col_var_type, double col_ub, 
		   double col_lb, char sense, int col_type, 
		   int col_bound_type);

int add_new_bounded_col(SRdesc *sr, double c_val, double a_val, int col_ind, 
			double rhs_ub_offset, double rhs_lb_offset, 
			double obj_ub_offset, double obj_lb_offset,
			double col_ub, double col_lb, int obj_sense);

int sr_find_opt_bounded(SRdesc *sr, double *ub, double *lb);

int sr_solve_open_prob(SRdesc *sr, int obj_ind, int row_ind, int *r_matbeg, 
		       int *r_matind, double *r_matval, COLinfo *cols, 
		       double *ub, double *lb, double etol);



int prep_light(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
	       prep_stats *stats, prep_params *prep_par);
int prep_advanced(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		  impList *implistL, impList *implistU, prep_stats *stats,
		  prep_params *prep_par);

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
//int prep_check_redundancy(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
//			  int verbosity);
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
