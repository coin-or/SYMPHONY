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

#ifndef _LPSOLVER_H
#define _LPSOLVER_H

#include "proto.h"
#include "BB_types.h"

#ifdef __CPLEX__

/*****************************************************************************/
/*******              here are the definitions for CPLEX               *******/
/*****************************************************************************/

#include <cplex.h>

/* The second comment indicates where the arrays are resized: BB or LP,
 * BB is BlackBox, LP is the lp solver specific part */

typedef struct TEMPORARY{
   char      *c;           /* max(2m,n) */
   int       *i1;          /* 3m+2n */
   int       *i2;          /* m */
   double    *d;           /* max(2m,2n) */
   void     **p1;          /* m */
   void     **p2;          /* m */

   char      *cv;          /* variable */
   int        cv_size;
   int       *iv;          /* variable (>= */
   int        iv_size;
   double    *dv;          /* variable */
   int        dv_size;
}temporary;

typedef struct LPdata{
   CPXENVptr  cpxenv;
   double     lpetol;
   CPXLPptr   lp;
   char       lp_is_modified;
   char       col_set_changed;
   double     objval;
   int        termcode;
   int        alloc_m;
   int        alloc_mplusn;
   int        alloc_mplusnz;
   int        n;           /* number of columns without slacks */
   int        maxn;
   int        m;           /* number of rows */
   int        maxm;
   int        nz;          /* number of nonzeros */
   int        maxnz;       /* space is allocated for this many nonzeros */
   int       *matbeg;      /* maxn + maxm + 1 */
   int       *matcnt;      /* maxn + maxm */
   int       *matind;      /* maxnz + maxm */
   double    *matval;      /* maxnz + maxm*/
   double    *obj;         /* maxn + maxm */
   double    *rhs;         /* maxm */
   double    *rngval;      /* maxm */
   char      *sense;       /* maxm */
   double    *lb;          /* maxn + maxm */
   double    *ub;          /* maxn + maxm */

   char       ordering;    /* COLIND_AND_USERIND_ORDERED, COLIND_ORDERED or
			      USERIND_ORDERED */
   var_desc **vars;        /* maxn */ /* BB */

   int        not_fixed_num;
   int       *not_fixed;
   int        nf_status;

   char      *status;      /* maxn */ /* BB */
   double    *x;           /* maxn */ /* BB */
   double    *dj;          /* maxn */ /* BB */
   double    *dualsol;     /* maxm */ /* BB */
   double    *slacks;      /* maxm */
   int       *bhead;       /* maxm */ /* BB */
   int        bhead_is_valid;
   double    *xbzero;      /* maxm */ /* BB */

   constraint  *rows;      /* maxm */

   temporary   tmp;
#ifdef PSEUDO_COSTS
   double     *pseudo_costs_one;
   double     *pseudo_costs_zero;
#endif
}LPdata;

void CPX_check_error PROTO((const char *erring_func));

#elif defined(__OSL__)
/*****************************************************************************/
/*******              here are the definitions for OSL                 *******/
/*****************************************************************************/
#include <ekk_c_api.h>

void OSL_check_error PROTO((const char *erring_func));

/* The second comment indicates where the arrays are resized: BB or LP,
 * BB is BlackBox, LP is the lp solver specific part */

typedef struct TEMPORARY{
   char      *c;           /* max(2m,n) */
   int       *i1;          /* 3m+2n */
   int       *i2;          /* m */
   double    *d;           /* max(2m,2n) */
   void     **p1;          /* m */
   void     **p2;          /* m */

   char      *cv;          /* variable */
   int        cv_size;
   int       *iv;          /* variable (>= */
   int        iv_size;
   double    *dv;          /* variable */
   int        dv_size;
}temporary;

typedef struct LPdata{
   EKKContext *env;
   double     lpetol;
   EKKModel   *lp;
   char       lp_is_modified;
   char       col_set_changed;
   double     objval;
   int        termcode;
   int        alloc_m;
   int        alloc_mplusn;
   int        alloc_mplusnz;
   int        n;           /* number of columns without slacks */
   int        maxn;
   int        m;           /* number of rows */
   int        maxm;
   int        nz;          /* number of nonzeros */
   int        maxnz;       /* space is allocated for this many nonzeros */
   int       *matbeg;      /* maxn + maxm + 1 */
   int       *matcnt;      /* maxn + maxm */
   int       *matind;      /* maxnz + maxm */
   double    *matval;      /* maxnz + maxm*/
   double    *obj;         /* maxn + maxm */
   double    *rhs;         /* maxm */
   double    *rngval;      /* maxm */
   char      *sense;       /* maxm */
   double    *lb;          /* maxn + maxm */
   double    *ub;          /* maxn + maxm */

   char       ordering;    /* COLIND_AND_USERIND_ORDERED, COLIND_ORDERED or
			      USERIND_ORDERED */
   var_desc **vars;        /* maxn */ /* BB */

   int        not_fixed_num;
   int       *not_fixed;
   int        nf_status;

   char      *status;      /* maxn */ /* BB */
   double    *x;           /* maxn */ /* BB */
   double    *dj;          /* maxn */ /* BB */
   double    *dualsol;     /* maxm */ /* BB */
   double    *slacks;      /* maxm */

   int       *bhead;       /* maxm */ /* BB */
   int        bhead_is_valid;
   double    *xbzero;      /* maxm */ /* BB */

   constraint  *rows;      /* maxm */

   temporary   tmp;
#ifdef PSEUDO_COSTS
   double     *pseudo_costs_one;
   double     *pseudo_costs_zero;
#endif
}LPdata;

#else

#error Unknown LP solver

#endif 

/*****************************************************************************/
/*******                  end LP solver definitions                    *******/
/*****************************************************************************/

/*****************************************************************************/
/*******                    common definitions                         *******/
/*****************************************************************************/

double dot_product PROTO((double *val, int *ind, int collen, double *col));
void open_lp_solver PROTO((LPdata *lp_data));
void close_lp_solver PROTO((LPdata *lp_data));
#ifdef COMPILE_CHECK_LP
void check_lp PROTO((LPdata *lp_data));
#endif
void load_lp_prob PROTO((LPdata *lp_data, int scaling, int fastmip));
void unload_lp_prob PROTO((LPdata *lp_data));
void load_basis PROTO((LPdata *lp_data, int *cstat, int *rstat));
void refactorize PROTO((LPdata *lp_data));
void add_rows PROTO((LPdata *lp_data, int rcnt, int nzcnt, double *rhs,
		     char *sense, int *rmatbeg, int *rmatind,double *rmatval));
void add_cols PROTO((LPdata *lp_data, int ccnt, int nzcnt, double *obj,
		     int *cmatbeg, int *cmatind, double *cmatval,
		     double *lb, double *ub, char *where_to_move));
void change_row PROTO((LPdata *lp_data, int row_ind,
		       char sense, double rhs, double range));
void change_col PROTO((LPdata *lp_data, int col_ind,
		       char sense, double lb, double ub));
int dual_simplex PROTO((LPdata *lp_data, int *iterd));
void btran PROTO((LPdata *lp_data, double *col));
void get_binvcol PROTO((LPdata *lp_data, int j, double *col));
void get_binvrow PROTO((LPdata *lp_data, int i, double *row));
void get_basis_header PROTO((LPdata *lp_data));
void get_basis PROTO((LPdata *lp_data, int *cstat, int *rstat));
void set_obj_upper_lim PROTO((LPdata *lp_data, double lim));
void set_itlim PROTO((LPdata *lp_data, int itlim));
void match_lp_solver_arrays_to_user PROTO((LPdata *lp_data, int allocm,
					   int allocn, int allocnz));
void resize_lp_solver_arrays PROTO((LPdata *lp_data));
void get_column PROTO((LPdata *lp_data, int j,
		       double *colval, int *colind, int *collen, double *cj));
void get_row PROTO((LPdata *lp_data, int i,
		    double *rowval, int *rowind, int *rowlen));
int get_proof_of_infeas PROTO((LPdata *lp_data, int *infind));
void get_x PROTO((LPdata *lp_data));
void get_dj_pi PROTO((LPdata *lp_data));
#if 0
void get_dualsol PROTO((LPdata *lp_data));
void get_reduced_costs PROTO((LPdata *lp_data));
#endif
void get_slacks PROTO((LPdata *lp_data));
void change_range PROTO((LPdata *lp_data, int rowind, double value));
void change_rhs PROTO((LPdata *lp_data,
		       int rownum, int *rhsind, double *rhsval));
void change_sense PROTO((LPdata *lp_data, int cnt, int *index, char *sense));
void change_bounds PROTO((LPdata *lp_data,
			  int cnt, int *index, char *lu, double *bd));
void change_lbub PROTO((LPdata *lp_data, int j, double lb, double ub));
void change_ub PROTO((LPdata *lp_data, int j, double ub));
void change_lb PROTO((LPdata *lp_data, int j, double lb));
void get_ub PROTO((LPdata *lp_data, int j, double *ub));
void get_lb PROTO((LPdata *lp_data, int j, double *lb));
void get_objcoef PROTO((LPdata *lp_data, int j, double *objcoef));
void delete_rows PROTO((LPdata *lp_data, int deletable, int *free_rows));
int delete_cols PROTO((LPdata *lp_data, int delnum, int *delstat));
void release_var PROTO((LPdata *lp_data, int j, int where_to_move));
void free_row_set PROTO((LPdata *lp_data, int length, int *index));
void constrain_row_set PROTO((LPdata *lp_data, int length, int *index));
void free_lp_solver_data PROTO((LPdata *lp_data, char arrays_too));
void write_mps PROTO((LPdata *lp_data, char *fname));
void write_sav PROTO((LPdata *lp_data, char *fname));

#endif
