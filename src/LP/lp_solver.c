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

#include <stdlib.h>              /* free() is here on AIX ... */
#include <malloc.h>
#include <math.h>

#include "lp_solver.h"
#include "BB_constants.h"
#include "BB_macros.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the interface with the LP Solver.
\*===========================================================================*/

double dot_product(double *val, int *ind, int collen, double *col)
{
   const int* lastind = ind + collen;
   double prod = 0;
   while (ind != lastind)
      prod += (*val++) * col[*ind++];
   return(prod);
}

#ifdef __OSL__

/*****************************************************************************/
/*****************************************************************************/
/*******                                                               *******/
/*******                  routines when OSL is used                    *******/
/*******                                                               *******/
/*******       WARNING! Not well tested. Please, report bugs.          *******/
/*****************************************************************************/
/*****************************************************************************/

/*============================================================================
 * - bhead a xbzero of LPdata structure are not filled
 * - no fastmip is used
 * - no scaling is used - cannot test it
 * - possible problems with getting reduced costs
 * - LPdata->tmp field mostly not used for safe. malloc and free for temporary
 *   fields are used instead => slow down
 *============================================================================
*/

/* Default value of OSL has set */
#define LP_MAX_ITER 9999999
/*
   lp_data->bhead and lp_data->xbzero are never used
   in the other code, so I don't load and maintain them
   anywhere in the lp_solver.c
*/

static int osllib_status;

#include <memory.h>

void OSL_check_error(const char *erring_func)
{
  if (osllib_status){
    printf("!!! OSL status is nonzero !!! [%s, %i]\n",
	   erring_func, osllib_status);
  }
}

/*===========================================================================*/
void open_lp_solver(LPdata *lp_data)
{
  EKKModel *baseModel;

  lp_data->env = ekk_initializeContext();
  osllib_status = (lp_data->env == NULL);
  OSL_check_error("open_lp_solver - ekk_initializeContext");
  baseModel = ekk_baseModel(lp_data->env);
  osllib_status = (baseModel == NULL);
  OSL_check_error("open_lp_solver - ekk_baseModel");
  ekk_setDebug(baseModel, -1, 0);
  ekk_setIloglevel(baseModel, 2);
/*  1    - 2999 informational messsages
    3000 - 5999 warn
    6000 - 6999 error, but keep running
    7000 - 8999 error and stop running */
/*    osllib_status = ekk_messagesPrintOn(baseModel, 1, 8999); */
/*    OSL_check_error("open_lp_solver - ekk_messagePrintOn"); */
  osllib_status = ekk_messagesPrintOff(baseModel, 1, 5999);
  OSL_check_error("open_lp_solver - ekk_messagePrintOff");

  /* default is to minimize */
/*    osllib_status = ekk_setMinimize(baseModel); */
/*    OSL_check_error("open_lp_solver - ekk_setMinimize"); */
  /* This should be infeasibility tolerance.*/
  lp_data->lpetol = ekk_getRtoldinf(baseModel);

  /* Speed up for large sparse problems. Test it, if it's faster or not. */
  osllib_status = ekk_setIuseRowCopy(baseModel, 1);
  OSL_check_error("open_lp_solver - ekk_setIuseRowCopy");
}

/*===========================================================================*/
void close_lp_solver(LPdata *lp_data)
{
  ekk_endContext(lp_data->env);

  /* Free up the lp_data parts */
  FREE(lp_data->matbeg);
  FREE(lp_data->matcnt);
  FREE(lp_data->matind);
  FREE(lp_data->matval);
  FREE(lp_data->obj);
  FREE(lp_data->rhs);
  FREE(lp_data->rngval);
  FREE(lp_data->sense);
  FREE(lp_data->lb);
  FREE(lp_data->ub);

  FREE(lp_data->not_fixed);
  FREE(lp_data->status);
  FREE(lp_data->x);
  FREE(lp_data->dj);
  FREE(lp_data->dualsol);
  FREE(lp_data->slacks);
  FREE(lp_data->bhead);
  FREE(lp_data->xbzero);
  FREE(lp_data->tmp.c);
  FREE(lp_data->tmp.i1);
  FREE(lp_data->tmp.i2);
  FREE(lp_data->tmp.d);
  FREE(lp_data->tmp.p1);
  FREE(lp_data->tmp.p2);
  FREE(lp_data->tmp.cv);
  FREE(lp_data->tmp.iv);
  FREE(lp_data->tmp.dv);
}

/*===========================================================================*/

#ifdef COMPILE_CHECK_LP
/* This implementation does not checks elements of matrix yet. */
void check_lp(LPdata *lp_data)
{
   int i;
   const double *d1, *d2, *ub, *lb;
   
   /* Chcecking the structure, not content */
   osllib_status =  ekk_validateModel(lp_data->lp);
   if (osllib_status != 0) {
      OSL_check_error("check_lp ekk_validateModel");
      exit(-1);
   }
   /* Check coeficients of objective function */
   osllib_status = memcmp(ekk_objective(lp_data->lp), lp_data->obj,
			  lp_data->n * DSIZE);
   OSL_check_error("check_lp objective");
   /* Check upper bounds of (structural) variables */
   osllib_status = memcmp(ekk_colupper(lp_data->lp), lp_data->ub,
			  lp_data->n * DSIZE);
   OSL_check_error("check_lp column upper bounds");
   /* Check upper bounds of (structural) variables */
   osllib_status = memcmp(ekk_collower(lp_data->lp), lp_data->lb,
			  lp_data->n * DSIZE);
   OSL_check_error("check_lp column lower bounds");
   
   osllib_status = -1;
   
   /* check the sense, upper bound and lower bound of rows*/
   d1 = lp_data->rhs;
   d2 = lp_data->rngval;
   ub = ekk_rowupper(lp_data->lp);
   lb = ekk_rowlower(lp_data->lp);
   for (i = lp_data->m - 1; i >= 0; i--) {
      switch (lp_data->sense[i]) {
       case 'E':
	 if ( (ub[i] != d1[i]) && (lb[i] != d1[i])) {
	    OSL_check_error("check_lp not an 'E' row");
	    return;
	 }
	 break;
       case 'L':
	 if ( (ub[i] != d1[i]) && (lb[i] != - OSL_INFINITY)) {
	    OSL_check_error("check_lp not an 'L' row");
	    return;
	 }
	 break;
       case 'G':
	 if ( (lb[i] != d1[i]) && (ub[i] != OSL_INFINITY)) {
	    OSL_check_error("check_lp not an 'G' row");
	    return;
	 }
	 break;
       case 'R':
	 if ( (ub[i] != d1[i]) && (ub[i] - lb[i] !=  d2[i]) ) {
	    OSL_check_error("check_lp not an 'R' row");
	    return;
	 }
	 break;
       default:
	 OSL_check_error("check_lp unknown sense");
	 return;
      }
   }
   osllib_status = 0;
}

#endif

/*===========================================================================*\
 * This function loads the data of an lp into the lp solver. 
\*===========================================================================*/
void load_lp_prob(LPdata *lp_data, int scaling, int fastmip)
{
   int i, *matcnt, *matbeg;
   double *lr = lp_data->tmp.d, *ur = lp_data->tmp.d + lp_data->n;
   
   lp_data->lp = ekk_newModel(lp_data->env, NULL);
   osllib_status = (lp_data->env == NULL);
   OSL_check_error("open_lp_solver - ekk_newModel");
   
   /* realloc_lp_arrays(lp_data); */
   
   matcnt = lp_data->matcnt;
   matbeg = lp_data->matbeg;
   
   /* OSL doesn't need matcnt as CPLEX.
      SYMPHONY need it only when writing MPS file. */
   for (i=lp_data->n-1; i>=0; i--)
      matcnt[i] = matbeg[i+1] - matbeg[i];
   
   for (i = 0; i < lp_data->m; i++) {
      switch (lp_data->sense[i]) {
       case 'E': lr[i] = ur[i] = lp_data->rhs[i]; break;
       case 'L': lr[i] = - OSL_INFINITY; ur[i] = lp_data->rhs[i]; break;
       case 'G': lr[i] = lp_data->rhs[i]; ur[i] = OSL_INFINITY; break;
       case 'R':
	 if (lp_data->rngval[i] >= 0) {
	    ur[i] = lp_data->rhs[i]; lr[i] = ur[i] - lp_data->rngval[i];
	 } else {
	    ur[i] = lp_data->rhs[i]; lr[i] = ur[i] + lp_data->rngval[i];
	 }
	 break;
       default: /* This should never happen ... */
	 osllib_status = -1;
	 OSL_check_error("load_lp - unknown sense");
      }
   }
   osllib_status =
      ekk_loadRimModel(lp_data->lp, lp_data->m, lr, ur,
		       lp_data->n, lp_data->obj, lp_data->lb, lp_data->ub);
   OSL_check_error("load_lp - ekk_loadRimModel");
   osllib_status =
      ekk_addColumnElementBlock(lp_data->lp, lp_data->n, lp_data->matind,
				lp_data->matbeg, lp_data->matval);
   OSL_check_error("load_lp - ekk_addColumnElementBlock");
   /* Not sure we need this since there's only one block */
   osllib_status = ekk_mergeBlocks(lp_data->lp, 1);
   OSL_check_error("load_lp - ekk_mergeBlocks");
   
   /* lp_data->scaling = scaling; */
}

/*===========================================================================*/

void unload_lp_prob(LPdata *lp_data)
{
   osllib_status = ekk_deleteModel(lp_data->lp);
   OSL_check_error("unload_lp - ekk_deleteModel");
   lp_data->lp = NULL;
}

/*===========================================================================*/

void load_basis(LPdata *lp_data, int *cstat, int *rstat)
{
   int *stat, i;
   
   if (cstat != NULL) {
      stat = ekk_getColstat(lp_data->lp);
      for (i = lp_data->n - 1; i >= 0; i--) {
	 stat[i] &= 0x1fffffff;
	 switch (cstat[i]) {
	  case VAR_BASIC: stat[i] |= 0x80000000; break;
	  case VAR_FREE: stat[i] |= 0x60000000; break;
	  case VAR_AT_UB: stat[i] |= 0x40000000; break;
	  case VAR_AT_LB: stat[i] |= 0x20000000; break;
	  case VAR_FIXED: stat[i] |= 0x00000000; break;
	  default: break; /* should never happen */
	 }
      }
      osllib_status = ekk_setColstat(lp_data->lp, stat);
      OSL_check_error("load_basis - ekk_setColstat");
      ekk_free(stat);
   }
   if (rstat != NULL) {
      stat = ekk_getRowstat(lp_data->lp);
      for (i = lp_data->m - 1; i >= 0; i--) {
	 stat[i] &= 0x1fffffff;
	 switch (rstat[i]) {
	  case SLACK_BASIC: stat[i] |= 0x80000000; break;
	  case SLACK_FREE: stat[i] |= 0x60000000; break;
	  case SLACK_AT_UB: stat[i] |= 0x40000000; break;
	  case SLACK_AT_LB: stat[i] |= 0x20000000; break;
	  case SLACK_FIXED: stat[i] |= 0x00000000; break;
	 }
      }
      osllib_status = ekk_setRowstat(lp_data->lp, stat);
      OSL_check_error("load_basis - ekk_setRowstat");
      ekk_free(stat);
   }
   lp_data->lp_is_modified = LP_HAS_NOT_BEEN_MODIFIED;
}

/*===========================================================================*/

void refactorize(LPdata *lp_data)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/

void add_rows(LPdata *lp_data, int rcnt, int nzcnt, double *rhs,
	      char *sense, int *rmatbeg, int *rmatind, double *rmatval)
{
   int i;
   double *lr, *ur;
   /* double *lr = lp_data->tmp.d, *ur = lp_data->tmp.d + lp_data->n; */

   lr = malloc(rcnt * DSIZE);
   ur = malloc(rcnt * DSIZE);
   for (i = rcnt - 1; i >= 0; i--) {
      switch (sense[i]) {
       case 'E': lr[i] = ur[i] = rhs[i]; break;
       case 'L': lr[i] = - OSL_INFINITY; ur[i] = rhs[i]; break;
       case 'G': lr[i] = rhs[i]; ur[i] = OSL_INFINITY; break;
       case 'R': lr[i] = ur[i] = lp_data->rhs[i]; break;
	 /* Range will be added later in change_range */
       default: /*This should never happen ... */
	 osllib_status = -1;
	 OSL_check_error("add_rows - unknown sense");
      }
   }
   osllib_status = ekk_addRows(lp_data->lp, rcnt, lr, ur, rmatbeg, rmatind,
			       rmatval);
   OSL_check_error("add_rows - ekk_addRows");
   
   /* Merge block can make comutation faster */
   osllib_status = ekk_mergeBlocks(lp_data->lp, 1);
   OSL_check_error("add_rows - ekk_mergeBlocks");
   
   FREE(lr);
   FREE(ur);
   
   lp_data->m += rcnt;
   lp_data->nz += nzcnt;
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void add_cols(LPdata *lp_data, int ccnt, int nzcnt, double *obj,
	      int *cmatbeg, int *cmatind, double *cmatval,
	      double *lb, double *ub, char *where_to_move)
{
   osllib_status = ekk_addColumns(lp_data->lp, ccnt, obj, lb, ub,
				  cmatbeg, cmatind, cmatval);
   OSL_check_error("add_cols - ekk_addColumns");
   osllib_status = ekk_mergeBlocks(lp_data->lp, 1);
   OSL_check_error("add_cols - ekk_mergeBlocks");
   lp_data->n += ccnt;
   lp_data->nz += nzcnt;
}

/*===========================================================================*/

void change_row(LPdata *lp_data, int row_ind,
		char sense, double rhs, double range)
{
   /*can be sped up using ekk_rowlower - direct acces to internal data*/
   double lr, ur;
   switch (sense) {
    case 'E': lr = ur = rhs; break;
    case 'L': lr = - OSL_INFINITY; ur = rhs; break;
    case 'G': lr = rhs; ur = OSL_INFINITY; break;
    case 'R':
      if (range >= 0) {
	 lr = rhs; ur = lr + range;
      } else {
	 ur = rhs; lr = ur + range;
      }
      break;
    default: /*This should never happen ... */
      osllib_status = -1;
      OSL_check_error("change_row - default");
   }
   osllib_status = ekk_copyRowlower(lp_data->lp, &lr, row_ind, row_ind + 1);
   OSL_check_error("change_row - ekk_copyRowlower");
   osllib_status = ekk_copyRowupper(lp_data->lp, &ur, row_ind, row_ind + 1);
   OSL_check_error("change_row - ekk_copyRowupper");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_col(LPdata *lp_data, int col_ind,
		char sense, double lb, double ub)
{
   switch (sense){
    case 'E': change_lbub(lp_data, col_ind, lb, ub); break;
    case 'R': change_lbub(lp_data, col_ind, lb, ub); break;
    case 'G': change_lb(lp_data, col_ind, lb); break;
    case 'L': change_ub(lp_data, col_ind, ub); break;
    default: /*This should never happen ... */
      osllib_status = -1;
      OSL_check_error("change_col - default");
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Solve the lp specified in lp_data->lp with dual simplex. The number of
 * iterations is returned in 'iterd'. The return value of the function is
 * the termination code of the dual simplex method.
\*===========================================================================*/
/* Basis head in the end of this function not finished yet */
int dual_simplex(LPdata *lp_data, int *iterd)
{
   int term;

   /*PreSolve seems to cause some problems -- not sure exactly why, but we
     leave it turned off for now. */
#if 0
   if ((osllib_status = ekk_preSolve(lp_data->lp, 3, NULL)) == 1){
      /* This means infeasibility was detected during preSolve */
      term = lp_data->termcode = D_UNBOUNDED;
      *iterd = 0;
      lp_data->lp_is_modified = LP_HAS_NOT_BEEN_MODIFIED;
      return(term);
   }
#endif
   
   if (lp_data->lp_is_modified == LP_HAS_BEEN_ABANDONED) {
      /* osllib_status = ekk_crash(lp_data->lp, 2); */
      /* OSL_check_error("dual_simplex - ekk_crash"); */
      osllib_status = ekk_allSlackBasis(lp_data->lp);
      OSL_check_error("dual_simplex - ekk_allSlackBasis");
   }
   
   term = ekk_dualSimplex(lp_data->lp);

   /*Turn postSolve back on if we figure out preSolve problem */
#if 0
   osllib_status = ekk_postSolve(lp_data->lp, NULL);*/
#endif
   
   /* We don't need this if we are not using preSolve. Not sure we need it
      anyway... */
#if 0
   /* Once more without preSolve. Dual simplex is run again
      only if solution is not optimal */
   if ((term == 1) || (term == 2)) {
      term = ekk_dualSimplex(lp_data->lp); 
      term = ekk_primalSimplex(lp_data->lp, 3);
   }
#endif
   
#if 0
   If (term == 2) {
      /* Dual infeas. This is impossible, so we must have had iteration
       * limit AND bound shifting AND dual feasibility not restored within
       * the given iteration limit. */
      maxiter = ekk_getImaxiter(lp_data->lp);
      osllib_status = ekk_setImaxiter(lp_data->lp, LP_MAX_ITER);
      OSL_check_error("dual_simplex - ekk_setImaxiter");
      term = ekk_dualSimplex(lp_data->lp);
      osllib_status = ekk_setImaxiter(lp_data->lp, maxiter);
      OSL_check_error("dual_simplex - ekk_setImaxiter");
   }
#endif

   switch (term) {
    case 0: term = OPTIMAL; break;
    case 1: term = D_UNBOUNDED;
      ekk_infeasibilities(lp_data->lp, 1, 1, NULL, NULL);
      break;
    case 2: term = D_INFEASIBLE; break;
    case 3: term = D_ITLIM; break;
    case 4:
      osllib_status = -1;
      OSL_check_error("osllib_status-ekk_dualSimplex found no solution!");
      exit(-1);
    case 5: D_OBJLIM; break;
    case 6:
      osllib_status = -1;
      OSL_check_error("osllib_status-ekk_dualSimplex lack of dstorage file"
			 "space!");
      exit(-1);
    default: term = ABANDONED;break;
   }
   
   lp_data->termcode = term;
   
   if (term != ABANDONED){
      *iterd = ekk_getIiternum(lp_data->lp);
      lp_data->objval = ekk_getRobjvalue(lp_data->lp);
      lp_data->lp_is_modified = LP_HAS_NOT_BEEN_MODIFIED;
   }else{
      lp_data->lp_is_modified = LP_HAS_BEEN_ABANDONED;
   }
   return(term);
}

/*===========================================================================*/

void btran(LPdata *lp_data, double *col)
{
   osllib_status = ekk_formBInverseTransposeb(lp_data->lp, col);
   OSL_check_error("btran - ekk_formBInverseTransposeb");
}

/*===========================================================================*/
/* This function is not used currently ...                                   */

void get_binvcol(LPdata *lp_data, int j, double *col)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/
/* This function is used only together with get_proof_of_infeasibility...    */

void get_binvrow(LPdata *lp_data, int i, double *row)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/
/* This function is never called either...                                   */

void get_basis_header(LPdata *lp_data)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/

void get_basis(LPdata *lp_data, int *cstat, int *rstat)
{
   int i, temp_stat;
   const int *stat;
   
   if (cstat != NULL) {
      stat = ekk_colstat(lp_data->lp);
      for (i = lp_data->n - 1; i >= 0; i--) {
	 if ((stat[i] & 0x80000000) != 0) {
	    cstat[i] = VAR_BASIC;
	 } else {
	    temp_stat = stat[i] & 0x60000000;
	    switch (temp_stat) {
	     case 0x60000000: cstat[i] = VAR_FREE; break;
	     case 0x40000000: cstat[i] = VAR_AT_UB; break;
	     case 0x20000000: cstat[i] = VAR_AT_LB; break;
	     case 0x00000000: cstat[i] = VAR_FIXED; break;
	    }
	 }
      }
   }
   if (rstat != NULL) {
      stat = ekk_rowstat(lp_data->lp);
      for (i = lp_data->m - 1; i >= 0; i--) {
	 if ((stat[i] & 0x80000000) != 0) {
	    rstat[i] = SLACK_BASIC;
	 } else {
	    temp_stat = stat[i] & 0x60000000;
	    switch (temp_stat) {
	     case 0x60000000: rstat[i] = SLACK_FREE; break;
	     case 0x40000000: rstat[i] = SLACK_AT_UB; break;
	     case 0x20000000: rstat[i] = SLACK_AT_LB; break;
	     case 0x00000000: rstat[i] = SLACK_FIXED; break;
	    }
	 }
      }
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Set an upper limit on the objective function value.
\*===========================================================================*/void set_obj_upper_lim(LPdata *lp_data, double lim)
{
   /* Not sure how to do this in OSL */
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/

/*===========================================================================*\
 * Set an upper limit on the number of iterations. If itlim < 0 then set
 * it to the maximum.
\*===========================================================================*/

void set_itlim(LPdata *lp_data, int itlim)
{
   if (itlim < 0) itlim = LP_MAX_ITER;
   osllib_status = ekk_setImaxiter(lp_data->lp, itlim);
   OSL_check_error("set_itlim - ekk_setImaxiter");
}

/*===========================================================================*/
void match_lp_solver_arrays_to_user(LPdata *lp_data,
				    int allocm, int allocn, int allocnz)
{
   if (allocm  != lp_data->alloc_m){
      FREE(lp_data->slacks);
      lp_data->slacks = (double *) malloc(allocm * DSIZE);
      lp_data->alloc_m = allocm;
   }
   
   if (allocn != lp_data->alloc_mplusn){
      FREE(lp_data->matcnt);
      FREE(lp_data->lb);
      FREE(lp_data->ub);
      lp_data->matcnt = (int *) malloc(allocn * ISIZE);
      lp_data->lb = (double *) malloc(allocn * DSIZE);
      lp_data->ub = (double *) malloc(allocn * DSIZE);
      lp_data->alloc_mplusn = allocn;
   }
   
   if (allocnz != lp_data->alloc_mplusnz)
      lp_data->alloc_mplusnz = allocnz;
}

/*===========================================================================*/

void resize_lp_solver_arrays(LPdata *lp_data)
{
   int resize = FALSE;
   
   if (lp_data->maxm > lp_data->alloc_m){
      lp_data->alloc_m = lp_data->maxm;
      lp_data->rhs = (double *) realloc((void *)lp_data->rhs,
					lp_data->alloc_m * DSIZE);
      lp_data->rngval = (double *) realloc((void *)lp_data->rngval,
					   lp_data->alloc_m * DSIZE);
      lp_data->sense = (char *) realloc((void *)lp_data->sense,
					lp_data->alloc_m * CSIZE);
      /*slacks must be allocated for cplex (for XMP it only points into x)*/
      lp_data->slacks = (double *) realloc((void *)lp_data->slacks,
					   lp_data->alloc_m * DSIZE);
      if (lp_data->lp != NULL) resize = TRUE;
   }
   if (lp_data->maxm + lp_data->maxn + 1 > lp_data->alloc_mplusn){
      lp_data->alloc_mplusn = lp_data->maxm + lp_data->maxn + 1;
      lp_data->matbeg = (int *) realloc((void *)lp_data->matbeg,
					lp_data->alloc_mplusn * ISIZE);
      lp_data->matcnt = (int *) realloc((void *)lp_data->matcnt,
					lp_data->alloc_mplusn * ISIZE);
      lp_data->obj = (double *) realloc((void *)lp_data->obj,
					lp_data->alloc_mplusn * DSIZE);
      lp_data->lb = (double *) realloc((void *)lp_data->lb,
				       lp_data->alloc_mplusn * DSIZE);
      lp_data->ub = (double *) realloc((void *)lp_data->ub,
				       lp_data->alloc_mplusn * DSIZE);
      if (lp_data->lp != NULL) resize = TRUE;
   }
   if (lp_data->maxnz + lp_data->maxm > lp_data->alloc_mplusnz){
      lp_data->alloc_mplusnz = lp_data->maxm + lp_data->maxnz;
      lp_data->matind = (int *) realloc((void *)lp_data->matind,
					lp_data->alloc_mplusnz * ISIZE);
      lp_data->matval = (double *) realloc((void *)lp_data->matval,
					   lp_data->alloc_mplusnz * DSIZE);
      if (lp_data->lp != NULL) resize = TRUE;
   }
   
   if (resize) {
      osllib_status = ekk_resizeModelExact(lp_data->lp, lp_data->m, lp_data->n,
					   lp_data->maxm, lp_data->maxn,
					   FALSE);
      OSL_check_error("resize_lp_solver_arrays ekk_resizeModelExact");
   }
}

/*===========================================================================*/

void get_column(LPdata *lp_data, int j,
		double *colval, int *colind, int *collen, double *cj)
{
   EKKVector vec;
   vec = ekk_getColumn(lp_data->lp, j);
   *collen = vec.numNonZero;
   memmove(colind, vec.index, *collen * ISIZE);
   memmove(colval, vec.element, *collen * DSIZE);
   ekk_freeVector(&vec);
   get_objcoef(lp_data, j, cj);
}

/*===========================================================================*/
void get_row(LPdata *lp_data, int i,
	     double *rowval, int *rowind, int *rowlen)
{
   EKKVector vec;
   vec = ekk_getRow(lp_data->lp, i);
   *rowlen = vec.numNonZero;
   memmove(rowind, vec.index, *rowlen * ISIZE);
   memmove(rowval, vec.element, *rowlen * DSIZE);
   ekk_freeVector(&vec);
}

/*===========================================================================*/
/* This routine returns the index of a row which proves the lp to be primal
 * infeasible. It is only needed when column generation is used.             */
/*===========================================================================*/int get_proof_of_infeas(LPdata *lp_data, int *infind)
{
  fprintf(stderr, "Function not implemented yet.");
  exit(-1);
}

/*===========================================================================*\
 * Get the solution (values of the structural variables in an optimal
 * solution) to the lp (specified by lp_data->lp) into the vector
 * lp_data->x.
\*===========================================================================*/
void get_x(LPdata *lp_data)
{
   memmove(lp_data->x, ekk_colsol(lp_data->lp), lp_data->n * DSIZE);
}

/*===========================================================================*/
void get_dj_pi(LPdata *lp_data)
{
   /*If scaling, fast integer or compress is used, maybe some changes will be
     needed */
   /* OSL returns changed sign - is it good or not? */
   memmove(lp_data->dualsol, ekk_rowduals(lp_data->lp), lp_data->m * DSIZE);

# if 0
   /* changing the sign */
   for (i = lp_data->m - 1; i >= 0; i --) {
      lp_data->dualsol[i] = - lp_data->dualsol[i];
   }
#endif
   
   memmove(lp_data->dj, ekk_colrcosts(lp_data->lp), lp_data->n * DSIZE);

#if 0
   for (i = lp_data->n - 1; i >= 0; i --) {
      lp_data->dj[i] = - lp_data->dj[i];
   }
#endif
}

/*===========================================================================*/
/* Possible improper implementetion. */

void get_slacks(LPdata *lp_data)
{
   constraint *rows = lp_data->rows;
   double *slacks = lp_data->slacks;
   const double *racts;
   int i, m = lp_data->m;
   
   racts = ekk_rowacts(lp_data->lp);
   
   for (i = m - 1; i >= 0; i--) {
      if ((rows[i].cut->sense == 'R') && (rows[i].cut->range < 0) ) {
	 slacks[i] = - rows[i].cut->rhs + racts[i];
      } else {
	 slacks[i] = rows[i].cut->rhs - racts[i];
      }
   }
}

/*===========================================================================*/

void change_range(LPdata *lp_data, int rowind, double value)
{
   const double *lrow, *urow;
   double lr, ur;
   lrow = ekk_rowlower(lp_data->lp);
   urow = ekk_rowupper(lp_data->lp);
   if (value >= 0) {
      lr = urow[rowind] - value;
   } else {
      lr = lrow[rowind] + value;
   }
   osllib_status = ekk_copyRowlower(lp_data->lp, &ur, rowind, rowind + 1);
   OSL_check_error("change_range - ekk_copyRowupper");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/
/* This function is never called ... */

void change_rhs(LPdata *lp_data, int rownum, int *rhsind, double *rhsval)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
   /* lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED; */
}

/*===========================================================================*/
/* This function is never called ...*/

void change_sense(LPdata *lp_data, int cnt, int *index, char *sense)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
   /* lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED; */
}

/*===========================================================================*/

void change_bounds(LPdata *lp_data, int cnt, int *index, char *lu, double *bd)
{
   double *lb, *ub;
   int i, j;
   lb = ekk_getCollower(lp_data->lp);
   ub = ekk_getColupper(lp_data->lp);
   for (i = cnt - 1; i >= 0; i--) {
      j = index[i];
      switch (lu[i]) {
      case 'L': lb[j] = bd[i];break;
      case 'U': ub[j] = bd[i];break;
      case 'B': lb[j] = ub[j] = bd[i];break;
      default: /*This should never happen ... */
	 osllib_status = -1;
	 OSL_check_error("change_bounds - default");
      }
   }
   osllib_status = ekk_setCollower(lp_data->lp, lb);
   OSL_check_error("change_bounds - ekk_setCollower");
   ekk_free(lb);
   osllib_status = ekk_setColupper(lp_data->lp, ub);
   OSL_check_error("change_bounds - ekk_setColupper");
   ekk_free(ub);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_lbub(LPdata *lp_data, int j, double lb, double ub)
{
   osllib_status = ekk_copyColupper(lp_data->lp, &ub, j, j + 1);
   OSL_check_error("change_lbub - ekk_copyColupper");
   osllib_status = ekk_copyCollower(lp_data->lp, &lb, j, j + 1);
   OSL_check_error("change_lbub - ekk_copyCollower");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_ub(LPdata *lp_data, int j, double ub)
{
   osllib_status = ekk_copyColupper(lp_data->lp, &ub, j, j + 1);
   OSL_check_error("change_ub - ekk_copyColupper");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_lb(LPdata *lp_data, int j, double lb)
{
   osllib_status = ekk_copyCollower(lp_data->lp, &lb, j, j + 1);
   OSL_check_error("change_lb - ekk_copyCollower");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void get_ub(LPdata *lp_data, int j, double *ub)
{
   /* Maybe some range checking could be added ...*/
   const double *uc = ekk_colupper(lp_data->lp);
   *ub = uc[j];
}

/*===========================================================================*/

void get_lb(LPdata *lp_data, int j, double *lb)
{
   /* Maybe some range checking could be added ...*/
   const double *lc = ekk_collower(lp_data->lp);
   *lb = lc[j];
}

/*===========================================================================*/

void get_objcoef(LPdata *lp_data, int j, double *objcoef)
{
   /* Maybe some range checking could be added ...*/
   const double *oc = ekk_objective(lp_data->lp);
   *objcoef = oc[j];
}

/*===========================================================================*/

void delete_rows(LPdata *lp_data, int delnum, int *free_rows)
{
   int i, j, m = lp_data->m, *which = lp_data->tmp.i1 + lp_data->m;

   /* which = calloc(delnum, ISIZE); */
   for (i = m - 1, j = 0; i >= 0; i--) if (free_rows[i]) which[j++] = i;
   osllib_status = ekk_deleteRows(lp_data->lp, j, which);
   OSL_check_error("delete_rows - ekk_deleteRows");
   /* FREE(which); */

   /* Make result as CPLEX does*/
   for (i = 0, j = 0; i < m; i++){
      if (free_rows[i])
	 free_rows[i] = -1;
      else
	 free_rows[i] = j++;
   }
   
   lp_data->m = j;
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

int delete_cols(LPdata *lp_data, int delnum, int *delstat)
{
   int i, j, n = lp_data->n, *which;
   
   which = calloc(delnum, ISIZE);
   for (i = n - 1, j = 0; i >= 0; i--) {
      if (delstat[i]) {
	 which[j++] = i;
      }
   }
   osllib_status = ekk_deleteColumns(lp_data->lp, j, which);
   OSL_check_error("delete_cols - ekk_deleteCols");
   FREE(which);
   
   lp_data->nz = ekk_getInumels(lp_data->lp);
   OSL_check_error("delete_cols - ekk_getInumels");
   
   /* make result as CPLEX does */
   for (i = 0, j = 0; i < lp_data->n; i++){
      if (delstat[i])
	 delstat[i] = -1;
      else
	 delstat[i] = j++;
   }
   
   lp_data->n = j;
   return j;
}

/*===========================================================================*/
/* Original (CPLEX) implementation is nothing :-)                            */
/*===========================================================================*/

void release_var(LPdata *lp_data, int j, int where_to_move)
{
#if 0
   switch (where_to_move){
   case MOVE_TO_UB:
      lp_data->lpbas.cstat[j] = 2; break; /* non-basic at its upper bound */
   case MOVE_TO_LB:
      lp_data->lpbas.cstat[j] = 0; break; /* non-basic at its lower bound */
   }
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
#endif
}

/*===========================================================================*/
/* There were some side effects setting "temp" fields of lp_data. */

void free_row_set(LPdata *lp_data, int length, int *index)
{
   int i, j;
   double *lb = (double *) ekk_getRowlower(lp_data->lp);
   double *ub = (double *) ekk_getRowupper(lp_data->lp);
   
   for (i = length - 1; i >= 0; i--) {
      j = index[i];
      lb[j] = - OSL_INFINITY;
      ub[j] = OSL_INFINITY;
   }
   osllib_status = ekk_setRowlower(lp_data->lp, lb);
   OSL_check_error("free_row_set ekk_setRowLower");
   ekk_free(lb);
   osllib_status = ekk_setRowupper(lp_data->lp, ub);
   OSL_check_error("free_row_set ekk_setRowUpper");
   ekk_free(ub);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/
/* There were some side effects setting "temp" fileds of lp_data. */

void constrain_row_set(LPdata *lp_data, int length, int *index)
{
   int i, j = 0;
   double *lb = ekk_getRowlower(lp_data->lp);
   double *ub = ekk_getRowupper(lp_data->lp);
   constraint *rows = lp_data->rows;
   cut_data *cut;
   
   for (i = length - 1; i >= 0; i--) {
      j = index[i];
      cut = rows[j].cut;
      switch(cut->sense) {
      case 'E': lb[j] = ub[j] = cut->rhs; break;
      case 'L': lb[j] = - OSL_INFINITY; ub[j] = cut->rhs; break;
      case 'G': lb[j] = cut->rhs; ub[j] = OSL_INFINITY; break;
      case 'R':
	 if (lp_data->rngval[j] >= 0) {
	    ub[j] = cut->rhs; lb[j] = ub[j] - lp_data->rngval[j];
	 } else {
	    ub[j] = cut->rhs; lb[j] = ub[j] + lp_data->rngval[j];
	 }
	 break;
      default: /*This should never happen ... */
	 osllib_status = -1;
	 OSL_check_error("load_lp - unknown type of constraint");
      }
   }
   
   j = 0;
   if (j)
      ekk_free(lb);
   
   osllib_status = ekk_setRowlower(lp_data->lp, lb);
   OSL_check_error("constrain_row_set ekk_setRowLower");
   ekk_free(lb);
   osllib_status = ekk_setRowupper(lp_data->lp, ub);
   OSL_check_error("constrain_row_set ekk_setRowUpper");
   ekk_free(ub);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/
void free_lp_solver_data(LPdata *lp_data, char arrays_too)
{
   if (lp_data){
      if (lp_data->lp != NULL) {
	 osllib_status = ekk_deleteModel(lp_data->lp);
	 OSL_check_error("free_lp_solver_data - ekk_deleteModel");
	 lp_data->lp = NULL;
      }
      if (arrays_too){
	 FREE(lp_data->matbeg);
	 FREE(lp_data->matcnt);
	 FREE(lp_data->matind);
	 FREE(lp_data->matval);
	 FREE(lp_data->obj);
	 FREE(lp_data->rhs);
	 FREE(lp_data->rngval);
	 FREE(lp_data->sense);
	 FREE(lp_data->lb);
	 FREE(lp_data->ub);
	 FREE(lp_data->slacks);
      }
   }
}

/*===========================================================================*/
void write_mps(LPdata *lp_data, char *fname)
{
   osllib_status = ekk_exportModel(lp_data->lp, fname, 1, 2);
   OSL_check_error("write_mps");
}

void write_sav(LPdata *lp_data, char *fname)
{
   osllib_status = ekk_saveModel(lp_data->lp, fname);
   OSL_check_error("write_sav");
}

#endif /* __OSL__ */


#ifdef __CPLEX__

/*****************************************************************************/
/*****************************************************************************/
/*******                                                               *******/
/*******                  routines when CPLEX is used                  *******/
/*******                                                               *******/
/*****************************************************************************/
/*****************************************************************************/

static int cpx_status;

#include <memory.h>

/*===========================================================================*/

void CPX_check_error(const char *erring_func)
{
   if (cpx_status){
      printf("!!! Cplex status is nonzero !!! [%s, %i]\n",
	     (char *)erring_func, cpx_status);
   }
}

/*===========================================================================*/

void open_lp_solver(LPdata *lp_data)
{
   int i;

   i = CPX_OFF;
   lp_data->cpxenv = CPXopenCPLEX(&cpx_status);
   CPX_check_error("open_lp_solver - error opening environment");
   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_SCRIND, i);
   CPX_check_error("open_lp_solver - CPXsetintparam, SCRIND");
   cpx_status = CPXgetdblparam(lp_data->cpxenv, CPX_PARAM_EPRHS,
			       &lp_data->lpetol);
   CPX_check_error("open_lp_solver - CPXgetdblparam");
}

/*===========================================================================*/

void close_lp_solver(LPdata *lp_data)
{
   cpx_status = CPXcloseCPLEX(&lp_data->cpxenv);
   CPX_check_error("close_lp_solver");

   /* Free up the lp_data parts */
   FREE(lp_data->matbeg);
   FREE(lp_data->matcnt);
   FREE(lp_data->matind);
   FREE(lp_data->matval);
   FREE(lp_data->obj);
   FREE(lp_data->rhs);
   FREE(lp_data->rngval);
   FREE(lp_data->sense);
   FREE(lp_data->lb);
   FREE(lp_data->ub);

   FREE(lp_data->not_fixed);
   FREE(lp_data->status);
   FREE(lp_data->x);
   FREE(lp_data->dj);
   FREE(lp_data->dualsol);
   FREE(lp_data->slacks);
   FREE(lp_data->bhead);
   FREE(lp_data->xbzero);
   FREE(lp_data->tmp.c);
   FREE(lp_data->tmp.i1);
   FREE(lp_data->tmp.i2);
   FREE(lp_data->tmp.d);
   FREE(lp_data->tmp.p1);
   FREE(lp_data->tmp.p2);
   FREE(lp_data->tmp.cv);
   FREE(lp_data->tmp.iv);
   FREE(lp_data->tmp.dv);
}

/*===========================================================================*/

#ifdef COMPILE_CHECK_LP

void check_lp(LPdata *lp_data)
{
#if CPX_VERSION <= 600 
   cpx_status = CPXchecklp(lp_data->cpxenv,
		(char *) "BB_prob", lp_data->n, lp_data->m, 1, lp_data->obj,
		lp_data->rhs, lp_data->sense, lp_data->matbeg, lp_data->matcnt,
		lp_data->matind, lp_data->matval, lp_data->lb, lp_data->ub,
		lp_data->rngval, lp_data->alloc_mplusn, lp_data->alloc_m,
		lp_data->alloc_mplusnz);
#else
   cpx_status = CPXcheckcopylp(lp_data->cpxenv, lp_data->lp,
		lp_data->n, lp_data->m, 1, lp_data->obj,
		lp_data->rhs, lp_data->sense, lp_data->matbeg, lp_data->matcnt,
		lp_data->matind, lp_data->matval, lp_data->lb, lp_data->ub,
		lp_data->rngval);
#endif
   CPX_check_error("load_lp");
}

#endif

/*===========================================================================*\
 * This function loads the data of an lp into the lp solver. This involves
 * transforming the data into CPLEX format and calling the CPLEX function
 * 'loadlp'.
\*===========================================================================*/

void load_lp_prob(LPdata *lp_data, int scaling, int fastmip)
{
   int i, *matcnt, *matbeg;

   /* realloc_lp_arrays(lp_data); */

   matcnt = lp_data->matcnt;
   matbeg = lp_data->matbeg;
   for (i=lp_data->n-1; i>=0; i--)
      matcnt[i] = matbeg[i+1] - matbeg[i];

   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_SCAIND, -1);
   CPX_check_error("load_lp - CPXsetintparam - SCAIND");

   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_FASTMIP, fastmip);
   CPX_check_error("load_lp - CPXsetintparam - FASTMIP");

   /* essentially disable basis snapshots */
   cpx_status =
      CPXsetintparam(lp_data->cpxenv, CPX_PARAM_BASINTERVAL, 2100000000);
   CPX_check_error("load_lp - CPXsetintparam - BASINTERVAL");

/* This is for the old memory model (user manages memory) */
#if CPX_VERSION <= 600 
   lp_data->lp =
      CPXloadlp(lp_data->cpxenv,
		(char *) "BB_prob", lp_data->n, lp_data->m, 1, lp_data->obj,
		lp_data->rhs, lp_data->sense, lp_data->matbeg, lp_data->matcnt,
		lp_data->matind, lp_data->matval, lp_data->lb, lp_data->ub,
		lp_data->rngval, lp_data->alloc_mplusn, lp_data->alloc_m,
		lp_data->alloc_mplusnz);
   CPX_check_error("load_lp - CPXloadlp");
#else /* This is for the new memory model (CPLEX manages memory) */
   lp_data->lp = CPXcreateprob(lp_data->cpxenv,&cpx_status,(char *) "BB_prob");
   CPX_check_error("load_lp - CPXcreateprob");
   cpx_status = CPXcopylp(lp_data->cpxenv, lp_data->lp,
		lp_data->n, lp_data->m, 1, lp_data->obj,
		lp_data->rhs, lp_data->sense, lp_data->matbeg, lp_data->matcnt,
		lp_data->matind, lp_data->matval, lp_data->lb, lp_data->ub,
		lp_data->rngval);
   CPX_check_error("load_lp - CPXcopylp");
#endif
}

/*===========================================================================*/

void unload_lp_prob(LPdata *lp_data)
{
   cpx_status = CPXfreeprob(lp_data->cpxenv, &lp_data->lp);
   CPX_check_error("unload_lp - CPXfreeprob");
   lp_data->lp = NULL;
}

/*===========================================================================*/

void load_basis(LPdata *lp_data, int *cstat, int *rstat)
{
#if CPX_VERSION <= 600 
   cpx_status = CPXloadbase(lp_data->cpxenv, lp_data->lp, cstat, rstat);
#else
   cpx_status = CPXcopybase(lp_data->cpxenv, lp_data->lp, cstat, rstat);
#endif
   CPX_check_error("load_basis - CPXloadbase");

   lp_data->lp_is_modified = LP_HAS_NOT_BEEN_MODIFIED;
}

/*===========================================================================*/

/* There should be something nicer... */
void refactorize(LPdata *lp_data)
{
   int itlim;

   cpx_status = CPXgetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, &itlim);
   CPX_check_error("refactorize - CPXgetintparam");
   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, 0);
   CPX_check_error("refactorize - CPXsetintparam");
#if CPX_VERSION <= 600 
   cpx_status = CPXoptimize(lp_data->cpxenv, lp_data->lp);
#else
   cpx_status = CPXprimopt(lp_data->cpxenv, lp_data->lp);
#endif
   CPX_check_error("refactorize - CPXoptimize");
   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, itlim);
   CPX_check_error("refactorize - CPXsetintparam");
}

/*===========================================================================*/

void add_rows(LPdata *lp_data, int rcnt, int nzcnt, double *rhs,
	      char *sense, int *rmatbeg, int *rmatind, double *rmatval)
{
   int i, m = lp_data->m, j, indicator = FALSE;
   int *bhead = lp_data->bhead + lp_data->m;

   if (indicator)
      for (i = 0; i < rcnt; i++){
	 printf("\n");
	 printf("%c %1f\n", sense[i], rhs[i]);
	 for (j = rmatbeg[i]; j < rmatbeg[i+1]; j++){
	    printf("%i ", rmatind[j]);
	 }
	 printf("\n");
	 for (j = rmatbeg[i]; j < rmatbeg[i+1]; j++){
	    printf("%1f ", rmatval[j]);
	 }
      }

   cpx_status = CPXaddrows(lp_data->cpxenv, lp_data->lp, 0, rcnt, nzcnt,
			   rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
   CPX_check_error("add_rows");
   lp_data->m += rcnt;
   lp_data->nz += nzcnt;
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;

   /* modify bhead */
   for (i = lp_data->m-1; i >= m; i--)
      bhead[i] = -i-1;
}

/*===========================================================================*/

void add_cols(LPdata *lp_data, int ccnt, int nzcnt, double *obj,
	      int *cmatbeg, int *cmatind, double *cmatval,
	      double *lb, double *ub, char *where_to_move)
{
   cpx_status = CPXaddcols(lp_data->cpxenv, lp_data->lp, ccnt, nzcnt,
	      obj, cmatbeg, cmatind, cmatval, lb, ub, NULL);
   CPX_check_error("add_cols");
   lp_data->n += ccnt;
   lp_data->nz += nzcnt;
}

/*===========================================================================*/

void change_row(LPdata *lp_data, int row_ind,
		char sense, double rhs, double range)
{
   cpx_status = CPXchgsense(lp_data->cpxenv, lp_data->lp, 1, &row_ind, &sense);
   CPX_check_error("change_row - CPXchgsense");
   cpx_status = CPXchgcoef(lp_data->cpxenv, lp_data->lp, row_ind, -1, rhs);
   CPX_check_error("change_row - CPXchgcoef");
   if (sense == 'R'){
      cpx_status = CPXchgcoef(lp_data->cpxenv, lp_data->lp, row_ind, -2,range);
      CPX_check_error("change_row - CPXchgcoef");
   }
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_col(LPdata *lp_data, int col_ind,
		char sense, double lb, double ub)
{
   switch (sense){
    case 'E': change_lbub(lp_data, col_ind, lb, ub); break;
    case 'R': change_lbub(lp_data, col_ind, lb, ub); break;
    case 'G': change_lb(lp_data, col_ind, lb); break;
    case 'L': change_ub(lp_data, col_ind, ub); break;
   }
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

/*===========================================================================*\
 * Solve the lp specified in lp_data->lp with dual simplex. The number of
 * iterations is returned in 'iterd'. The return value of the function is
 * the termination code of the dual simplex method.
\*===========================================================================*/

int dual_simplex(LPdata *lp_data, int *iterd)
{
   int real_term, term, itlim, defit, minit, maxit;
   double objulim, objllim, defobj;

   if (lp_data->lp_is_modified == LP_HAS_BEEN_ABANDONED){
      cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ADVIND, CPX_OFF);
      CPX_check_error("dual_simplex - CPXsetintparam, ADVIND");
   }

   term = CPXdualopt(lp_data->cpxenv, lp_data->lp);
   if (term == CPXERR_PRESLV_INForUNBD){
      cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_PREIND, CPX_OFF);
      CPX_check_error("dual_simplex - CPXsetintparam");
      term = CPXdualopt(lp_data->cpxenv, lp_data->lp);
      CPX_check_error("dual_simplex - CPXdualopt");
      cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_PREIND, CPX_ON);
      CPX_check_error("dual_simplex - CPXsetintparam");
   }

   term = CPXgetstat(lp_data->cpxenv,lp_data->lp);
#if CPX_VERSION >= 800
   if (term == CPX_STAT_INFEASIBLE){
#else
   if (term == CPX_INFEASIBLE){
#endif
      /* Dual infeas. This is impossible, so we must have had iteration
       * limit AND bound shifting AND dual feasibility not restored within
       * the given iteration limit. */
      cpx_status = CPXgetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, &itlim);
      CPX_check_error("dual_simplex - CPXgetintparam, ITLIM");
      cpx_status = CPXinfointparam(lp_data->cpxenv, CPX_PARAM_ITLIM,
				   &defit, &minit, &maxit);
      CPX_check_error("dual_simplex - CPXinfointparam, ITLIM");
      cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, defit);
      CPX_check_error("dual_simplex - CPXsetintparam, ITLIM");
      cpx_status = CPXgetdblparam(lp_data->cpxenv, CPX_PARAM_OBJULIM,&objulim);
      CPX_check_error("dual_simplex - CPXgetdblparam, OBJULIM");
      cpx_status = CPXgetdblparam(lp_data->cpxenv, CPX_PARAM_OBJULIM,&objllim);
      CPX_check_error("dual_simplex - CPXgetdblparam, OBJULIM");
      defobj = 1e75;
      cpx_status = CPXsetdblparam(lp_data->cpxenv, CPX_PARAM_OBJULIM, defobj);
      CPX_check_error("dual_simplex - CPXsetdblparam, OBJULIM");
      defobj = -1e75;
      cpx_status = CPXsetdblparam(lp_data->cpxenv, CPX_PARAM_OBJLLIM, defobj);
      CPX_check_error("dual_simplex - CPXsetdblparam, OBJLLIM");
      term = CPXdualopt(lp_data->cpxenv, lp_data->lp);
      cpx_status = CPXsetdblparam(lp_data->cpxenv, CPX_PARAM_OBJULIM, objulim);
      CPX_check_error("dual_simplex - CPXsetdblparam, OBJULIM");
      cpx_status = CPXsetdblparam(lp_data->cpxenv, CPX_PARAM_OBJLLIM, objllim);
      CPX_check_error("dual_simplex - CPXsetdblparam, OBJLLIM");
      cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, itlim);
      CPX_check_error("dual_simplex - CPXsetintparam, ITLIM");
   }

#if CPX_VERSION >= 800
   switch (real_term = CPXgetstat(lp_data->cpxenv,lp_data->lp)){
    case CPX_STAT_OPTIMAL:                        term = OPTIMAL; break;
    case CPX_STAT_INFEASIBLE:                     term = D_UNBOUNDED; break;
    case CPX_STAT_UNBOUNDED:                      term = D_INFEASIBLE; break;
    case CPX_STAT_ABORT_OBJ_LIM:                  term = D_OBJLIM; break;
    case CPX_STAT_ABORT_IT_LIM:                   term = D_ITLIM; break;
    default:                                      term = ABANDONED; break;
   }
#else
   switch (real_term = CPXgetstat(lp_data->cpxenv,lp_data->lp)){
    case CPX_OPTIMAL:                             term = OPTIMAL; break;
    case CPX_INFEASIBLE:                          term = D_INFEASIBLE; break;
    case CPX_UNBOUNDED:                           term = D_UNBOUNDED; break;
    case CPX_OBJ_LIM:                             term = D_OBJLIM; break;
    case CPX_IT_LIM_FEAS: case CPX_IT_LIM_INFEAS: term = D_ITLIM; break;
    default:                                      term = ABANDONED; break;
   }
#endif

   lp_data->termcode = term;

   if (term != ABANDONED){
      cpx_status = CPXgetbhead(lp_data->cpxenv, lp_data->lp,
			       lp_data->bhead, lp_data->xbzero);
      lp_data->bhead_is_valid = TRUE;
      CPX_check_error("dual_simplex - CPXgetbhead");
#if CPX_VERSION <= 600 
      *iterd = CPXgetitc(lp_data->cpxenv, lp_data->lp);
#else
      *iterd = CPXgetitcnt(lp_data->cpxenv, lp_data->lp);
#endif
      cpx_status = CPXgetobjval(lp_data->cpxenv,lp_data->lp, &lp_data->objval);
      CPX_check_error("dual_simplex - CPXgetobjval");
      cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ADVIND, CPX_ON);
      CPX_check_error("dual_simplex - CPXsetintparam, ADVIND");
      lp_data->lp_is_modified = LP_HAS_NOT_BEEN_MODIFIED;
   }else{
      lp_data->lp_is_modified = LP_HAS_BEEN_ABANDONED;
      printf("CPLEX Abandoned calculation: Code %i \n\n", real_term);
   }
   return(term);
}

/*===========================================================================*/

void btran(LPdata *lp_data, double *col)
{
   cpx_status = CPXbtran(lp_data->cpxenv, lp_data->lp, col);
   CPX_check_error("btran");
}

/*===========================================================================*/

void get_binvcol(LPdata *lp_data, int j, double *col)
{
   cpx_status = CPXbinvcol(lp_data->cpxenv, lp_data->lp, j, col);
   CPX_check_error("get_binvcol");
}

/*===========================================================================*/

void get_binvrow(LPdata *lp_data, int i, double *row)
{
   cpx_status = CPXbinvrow(lp_data->cpxenv, lp_data->lp, i, row);
   CPX_check_error("get_binvrow");
}

/*===========================================================================*/

void get_basis_header(LPdata *lp_data)
{
   cpx_status = CPXgetbhead(lp_data->cpxenv, lp_data->lp,
			    lp_data->bhead, lp_data->xbzero);
   CPX_check_error("get_basis_header - CPXgetbhead");
}

/*===========================================================================*/

void get_basis(LPdata *lp_data, int *cstat, int *rstat)
{
   cpx_status = CPXgetbase(lp_data->cpxenv, lp_data->lp, cstat, rstat);
   CPX_check_error("get_basis");
}

/*===========================================================================*/

/*===========================================================================*\
 * Set an upper limit on the objective function value. Call the 'setobjulim'
 * CPLEX function.
\*===========================================================================*/

void set_obj_upper_lim(LPdata *lp_data, double lim)
{
   cpx_status = CPXsetdblparam(lp_data->cpxenv, CPX_PARAM_OBJULIM, lim);
   CPX_check_error("set_obj_upper_lim");
}

/*===========================================================================*/

/*===========================================================================*\
 * Set an upper limit on the number of iterations. If itlim < 0 then set
 * it to the maximum.
\*===========================================================================*/

void set_itlim(LPdata *lp_data, int itlim)
{
   if (itlim < 0)
      cpx_status = CPXinfointparam(lp_data->cpxenv,
				   CPX_PARAM_ITLIM, &itlim, NULL, NULL);
   CPX_check_error("set_itlim - CPXinfointparam");
   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, itlim);
   CPX_check_error("set_itlim - CPXsetintparam");
}

/*===========================================================================*/

void match_lp_solver_arrays_to_user(LPdata *lp_data,
				    int allocm, int allocn, int allocnz)
{
   if (allocm  != lp_data->alloc_m){
      FREE(lp_data->slacks);
      lp_data->slacks = (double *) malloc(allocm * DSIZE);
      lp_data->alloc_m = allocm;
   }

   if (allocn != lp_data->alloc_mplusn){
      FREE(lp_data->matcnt);
      FREE(lp_data->lb);
      FREE(lp_data->ub);
      lp_data->matcnt = (int *) malloc(allocn * ISIZE);
      lp_data->lb = (double *) malloc(allocn * DSIZE);
      lp_data->ub = (double *) malloc(allocn * DSIZE);
      lp_data->alloc_mplusn = allocn;
   }

   if (allocnz != lp_data->alloc_mplusnz)
      lp_data->alloc_mplusnz = allocnz;
}

/*===========================================================================*/

void resize_lp_solver_arrays(LPdata *lp_data)
{
   char resize = FALSE;

   if (! lp_data->lp){

      if (lp_data->maxm > lp_data->alloc_m){
	 lp_data->alloc_m = lp_data->maxm;
	 lp_data->rhs = (double *) realloc((void *)lp_data->rhs,
					   lp_data->alloc_m * DSIZE);
	 lp_data->rngval = (double *) realloc((void *)lp_data->rngval,
					      lp_data->alloc_m * DSIZE);
	 lp_data->sense = (char *) realloc((void *)lp_data->sense,
					   lp_data->alloc_m * CSIZE);
	 /*slacks must be allocated for cplex (for XMP it only points into x)*/
	 lp_data->slacks = (double *) realloc((void *)lp_data->slacks,
					      lp_data->alloc_m * DSIZE);
      }
      if (lp_data->maxm + lp_data->maxn + 1 > lp_data->alloc_mplusn){
	 lp_data->alloc_mplusn = lp_data->maxm + lp_data->maxn + 1;
	 lp_data->matbeg = (int *) realloc((void *)lp_data->matbeg,
					   lp_data->alloc_mplusn * ISIZE);
	 lp_data->matcnt = (int *) realloc((void *)lp_data->matcnt,
					   lp_data->alloc_mplusn * ISIZE);
	 lp_data->obj = (double *) realloc((void *)lp_data->obj,
					   lp_data->alloc_mplusn * DSIZE);
	 lp_data->lb = (double *) realloc((void *)lp_data->lb,
					  lp_data->alloc_mplusn * DSIZE);
	 lp_data->ub = (double *) realloc((void *)lp_data->ub,
					  lp_data->alloc_mplusn * DSIZE);
      }
      if (lp_data->maxnz + lp_data->maxm > lp_data->alloc_mplusnz){
	 lp_data->alloc_mplusnz = lp_data->maxm + lp_data->maxnz;
	 lp_data->matind = (int *) realloc((void *)lp_data->matind,
					   lp_data->alloc_mplusnz * ISIZE);
	 lp_data->matval = (double *) realloc((void *)lp_data->matval,
					      lp_data->alloc_mplusnz * DSIZE);
      }

   }else{

      if (lp_data->maxm > lp_data->alloc_m){
	 lp_data->alloc_m = lp_data->maxm;
	 lp_data->slacks =
	    (double *) realloc(lp_data->slacks, lp_data->alloc_m * DSIZE);
	 resize = TRUE;
      }
      if (lp_data->maxm + lp_data->maxn + 1 > lp_data->alloc_mplusn){
	 lp_data->alloc_mplusn = lp_data->maxm + lp_data->maxn + 1;
	 resize = TRUE;
      }
      if (lp_data->maxnz + lp_data->maxm > lp_data->alloc_mplusnz){
	 lp_data->alloc_mplusnz = lp_data->maxm + lp_data->maxnz;
	 resize = TRUE;
      }
      if (resize){
#if CPX_VERSION <= 600 
	 /* This is only needed in the old memory model */
	 cpx_status =
	    CPXreallocprob(lp_data->cpxenv, lp_data->lp, &lp_data->obj,
			   &lp_data->rhs, &lp_data->sense, &lp_data->matbeg,
			   &lp_data->matcnt, &lp_data->matind,&lp_data->matval,
			   &lp_data->lb, &lp_data->ub, &lp_data->rngval,
			   NULL, NULL, NULL, NULL, NULL,
			   lp_data->alloc_mplusn, lp_data->alloc_m,
			   lp_data->alloc_mplusnz, 0, 0);
	 CPX_check_error("resize_lp");
#endif
      }
   }
}

/*===========================================================================*/

void get_column(LPdata *lp_data, int j,
		double *colval, int *colind, int *collen, double *cj)
{
   int matbeg, surplus;
   /* If there was no scaling, then we could probably copy the data out
    * directly. Try sometime... */
   cpx_status = CPXgetcols(lp_data->cpxenv, lp_data->lp, collen, &matbeg,
			   colind, colval, lp_data->m, &surplus, j, j);
   CPX_check_error("get_column - CPXgetcols");
   cpx_status = CPXgetobj(lp_data->cpxenv, lp_data->lp, cj, j, j);
   CPX_check_error("get_column - CPXgetobj");
}

/*===========================================================================*/

void get_row(LPdata *lp_data, int i,
	     double *rowval, int *rowind, int *rowlen)
{
   int rmatbeg, surplus;
   /* If there was no scaling, then we could probably copy the data out
    * directly. Try sometime... */
   cpx_status = CPXgetrows(lp_data->cpxenv, lp_data->lp, rowlen, &rmatbeg,
			   rowind, rowval, lp_data->n, &surplus, i, i);
   CPX_check_error("get_row - CPXgetrows");
}

/*===========================================================================*/

/* This routine returns the index of a row which proves the lp to be primal
 * infeasible. There must be one, or this function wouldn't be called. */
/* There MUST be something better than this...
 * A function call perhaps... Ask CPLEX... */
int get_proof_of_infeas(LPdata *lp_data, int *infind)
{
   int idiv, jdiv;
   double bd;

#if 0
   /*something like this should work...*/
   CPXdualfarkas(lp_data->cpxenv, lp_data->lp, ...);
   CPX_check_error("get_proof_of_infeas - CPXdualfarkas");
#endif

   CPXgetijdiv(lp_data->cpxenv, lp_data->lp, &idiv, &jdiv);
   CPX_check_error("get_proof_of_infeas - CPXgetijdiv");
   cpx_status = CPXgetijrow(lp_data->cpxenv, lp_data->lp, idiv, jdiv, infind);
   CPX_check_error("get_proof_of_infeas - CPXgetijrow");
   if (cpx_status)
      return(0);
   if (jdiv < 0){ /* the diverging variable is a slack/range */
      if (lp_data->slacks)
	 return(lp_data->slacks[idiv] < 0 ? LOWER_THAN_LB : HIGHER_THAN_UB);
   }else{ /* the diverging variable is structural */
      cpx_status = CPXgetlb(lp_data->cpxenv, lp_data->lp, &bd, jdiv, jdiv);
      CPX_check_error("get_proof_of_infeas - CPXgetlb");
      if(lp_data->x)
	 return(bd < lp_data->x[jdiv] ? LOWER_THAN_LB : HIGHER_THAN_UB);
   }
   return(0); /* fake return */
}

#if 0
/*****************************************************************************/
/*                 The old version below is commented out !!!!               */
/*****************************************************************************/
int get_proof_of_infeas(LPdata *lp_data, int *infind)
{
   int i, j, k, n = lp_data->n;
   double lb, ub, etol, *xbzero = lp_data->xbzero;
   int *bhead = lp_data->bhead;
   double *z = lp_data->dtmpn;
   /* Let's hope these are correctly updated */
   int *cstat = lp_data->lpbas.cstat;

   cpx_status = CPXgetdblparam(lp_data->cpxenv, CPX_PARAM_EPRHS, &etol);
   CPX_check_error("get_proof_of_infeas - CPXgetdblparam");

   refactorize(lp_data, NULL);

   cpx_status = CPXgetbhead(lp_data->cpxenv, lp_data->lp, bhead, xbzero);
   CPX_check_error("get_proof_of_infeas - CPXgetbhead");
   for (j = 0, i = lp_data->m -1; i >= 0; i--){
      if (bhead[i] >= 0){ /* a structural variable */
	 cpx_status = CPXgetlb(lp_data->cpxenv, lp_data->lp, &lb, bhead[i],
			       bhead[i]);
	 CPX_check_error("get_proof_of_infeas - CPXgetlb");
	 cpx_status = CPXgetub(lp_data->cpxenv, lp_data->lp, &ub, bhead[i],
			       bhead[i]);
	 CPX_check_error("get_proof_of_infeas - CPXgetub");
      }else{ /* slack variable */
	 lb = 0;
	 if (lp_data->sense[i] != 'R'){
	    ub = INFINITY;
	 }else{
	    /* the range is the ub on a slack variable */
	    cpx_status = CPXgetcoef(lp_data->cpxenv, lp_data->lp, i, -2, &ub);
	    CPX_check_error("get_proof_of_infeas - CPXgetcoef");
	 }
      }

      if (xbzero[i] <= lb - etol){ /* LOWER_THAN_LB */
	 cpx_status = CPXbinvarow(lp_data->cpxenv, lp_data->lp, i, z);
	 CPX_check_error("get_proof_of_infeas - CPXbinvarow");
	 for (k=n-1; k>=0; k--){
	    /* We can pivot if an out of basis variable is
	     * either at its lower bound and the matrix element is negative
	     * or at is upper bound and the matrix element is positive */
	    if ((cstat[k] == 0 && z[k] < -etol) ||
		(cstat[k] == 2 && z[k] > etol))
	       break;
	 }
	 if (k < 0){ /* We have proof of infeasibility */
	    *infind = i;
	    return(LOWER_THAN_LB);
	 }
 	 continue;
      }

      if (xbzero[i] >= ub + etol){ /* HIGHER_THAN_UB */
	 cpx_status = CPXbinvarow(lp_data->cpxenv, lp_data->lp, i, z);
	 CPX_check_error("get_proof_of_infeas - CPXbinvarow");
	 for (k=n-1; k>=0; k--){
	    /* We can pivot if an out of basis variable is
	     * either at its lower bound and the matrix element is positive
	     * or at is upper bound and the matrix element is negative */
	    if ((cstat[k] == 0 && z[k] > etol) ||
		(cstat[k] == 2 && z[k] < -etol))
	       break;
	 }
	 if (k < 0){ /* We have proof of infeasibility */
	    *infind = i;
	    return(HIGHER_THAN_UB);
	 }
      }
   }
   return(0); /* fake return */
}
#endif

/*===========================================================================*/

/*===========================================================================*\
 * Get the solution (values of the structural variables in an optimal
 * solution) to the lp (specified by lp_data->lp) into the vector
 * lp_data->x. This can be done by calling the 'getx' CPLEX function.
\*===========================================================================*/

void get_x(LPdata *lp_data)
{
   cpx_status = CPXgetx(lp_data->cpxenv, lp_data->lp, lp_data->x, 0,
			lp_data->n-1);
   CPX_check_error("get_x");
}

/*===========================================================================*/

void get_dj_pi(LPdata *lp_data)
{
#if CPX_VERSION <= 600 
   /* This mess is only needed for old versions of CPLEX */
   int i;
   cpx_status = CPXgetintparam(lp_data->cpxenv, CPX_PARAM_FASTMIP, &i);
   CPX_check_error("get_dj_pi - CPXgetintparam, FASTMIP");
   if (i == CPX_OFF) {
      cpx_status = CPXgetpi(lp_data->cpxenv, lp_data->lp, lp_data->dualsol, 0,
			    lp_data->m-1);
      CPX_check_error("get_dj_pi - CPXgetpi");
      cpx_status = CPXgetdj(lp_data->cpxenv, lp_data->lp, lp_data->dj, 0,
			    lp_data->n-1);
      CPX_check_error("get_dj_pi - CPXgetdj");
   } else {
      double *pi = lp_data->dualsol;
      double *dj = lp_data->dj;
      int *bhead = lp_data->bhead;
      int scaling;
      cpx_status = CPXgetintparam(lp_data->cpxenv, CPX_PARAM_SCAIND, &scaling);
      CPX_check_error("get_dj_pi - CPXgetintparam, SCAIND");
      if (scaling == -1) {
	 /* No scaling. We can pick the data directly from obj and the matrix*/
	 double *obj = lp_data->obj;
	 int *matbeg = lp_data->matbeg;
	 int *matind = lp_data->matind;
	 int *matcnt = lp_data->matcnt;
	 double *matval = lp_data->matval;
	 /* Get the duals first */
	 for (i = lp_data->m - 1; i >= 0; --i)
	    pi[i] = bhead[i] >= 0 ? obj[bhead[i]] : 0;
	 cpx_status = CPXbtran(lp_data->cpxenv, lp_data->lp, pi);
	 CPX_check_error("get_dj_pi - CPXbtran");
	 /* Now get the reduced costs */
	 for (i = lp_data->n - 1; i >= 0; --i) {
	    dj[i] = obj[i] - dot_product(matval+matbeg[i], matind+matbeg[i],
					 matcnt[i], pi);
	 }
      } else {
	 /* Scaling was done. We must get the data through cplex */
	 /* Get the duals first */
	 double *obj = (double *) malloc(lp_data->n * DSIZE);
	 int collen;
	 int matbeg, surplus; /* "dummy" variables */
	 int *colind = (int *) malloc(lp_data->m * ISIZE);
	 double *colval = (double *) malloc(lp_data->m * DSIZE);
	 /* Get the duals first */
	 cpx_status = CPXgetobj(lp_data->cpxenv, lp_data->lp, obj, 0,
				lp_data->n - 1);
	 CPX_check_error("get_dj_pi - CPXgetobj");
	 for (i = lp_data->m - 1; i >= 0; --i)
	    pi[i] = bhead[i] >= 0 ? obj[bhead[i]] : 0;
	 cpx_status = CPXbtran(lp_data->cpxenv, lp_data->lp, pi);
	 CPX_check_error("get_dj_pi - CPXbtran");
	 /* Now get the reduced costs */
	 for (i = lp_data->n - 1; i >= 0; --i) {
	    cpx_status = CPXgetcols(lp_data->cpxenv, lp_data->lp, &collen,
				    &matbeg, colind, colval, lp_data->m,
				    &surplus, i, i);
	    CPX_check_error("get_dj_pi - CPXgetcols");
	    dj[i] = obj[i] - dot_product(colval, colind, collen, pi);
	 }
      }
   }
#else /* We have a newer version of CPLEX */
   cpx_status = CPXgetpi(lp_data->cpxenv, lp_data->lp, lp_data->dualsol, 0,
			 lp_data->m-1);
   CPX_check_error("get_dj_pi - CPXgetpi");
   cpx_status = CPXgetdj(lp_data->cpxenv, lp_data->lp, lp_data->dj, 0,
			 lp_data->n-1);
   CPX_check_error("get_dj_pi - CPXgetdj");
#endif
}

/*===========================================================================*/
#if 0
void get_dualsol(LPdata *lp_data)
{
   cpx_status = CPXgetpi(lp_data->cpxenv, lp_data->lp, lp_data->dualsol, 0,
			 lp_data->m-1);
   CPX_check_error("get_dualsol");
}

/*===========================================================================*/

void get_reduced_costs(LPdata *lp_data)
{
   cpx_status = CPXgetdj(lp_data->cpxenv, lp_data->lp, lp_data->dj, 0,
			 lp_data->n-1);
   CPX_check_error("get_reduced_costs");
}
#endif
/*===========================================================================*/

void get_slacks(LPdata *lp_data)
{
   constraint *rows = lp_data->rows;
   double *slacks = lp_data->slacks;
   int i, m = lp_data->m;
   cpx_status = CPXgetslack(lp_data->cpxenv, lp_data->lp, lp_data->slacks, 0,
			    lp_data->m-1);
   CPX_check_error("get_slacks");
   /* Compute the real slacks for the free rows */
   for (i=m-1; i>=0; i--){
      if (rows[i].free){
	 switch (rows[i].cut->sense){
	  case 'E': slacks[i] +=  rows[i].cut->rhs - INFINITY; break;
	  case 'L': slacks[i] +=  rows[i].cut->rhs - INFINITY; break;
	  case 'G': slacks[i] +=  rows[i].cut->rhs + INFINITY; break;
	  case 'R': slacks[i] += -rows[i].cut->rhs - INFINITY; break;
	 }
      }
   }

}

/*===========================================================================*/

void change_range(LPdata *lp_data, int rowind, double value)
{
   cpx_status = CPXchgcoef(lp_data->cpxenv, lp_data->lp, rowind, -2, value);
   CPX_check_error("change_range");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_rhs(LPdata *lp_data, int rownum, int *rhsind, double *rhsval)
{
   cpx_status = CPXchgrhs(lp_data->cpxenv, lp_data->lp, rownum, rhsind,rhsval);
   CPX_check_error("change_rhs");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_sense(LPdata *lp_data, int cnt, int *index, char *sense)
{
   cpx_status = CPXchgsense(lp_data->cpxenv, lp_data->lp, cnt, index, sense);
   CPX_check_error("change_sense");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_bounds(LPdata *lp_data, int cnt, int *index, char *lu, double *bd)
{
   cpx_status = CPXchgbds(lp_data->cpxenv, lp_data->lp, cnt, index, lu, bd);
   CPX_check_error("change_bounds");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_lbub(LPdata *lp_data, int j, double lb, double ub)
{
   int ind[2];
   double bd[2];
   ind[0] = ind[1] = j;
   bd[0] = lb; bd[1] = ub;
   cpx_status =
      CPXchgbds(lp_data->cpxenv, lp_data->lp, 2, ind, (char *)"LU", bd);
   CPX_check_error("change_lbub");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_ub(LPdata *lp_data, int j, double ub)
{
   cpx_status = CPXchgbds(lp_data->cpxenv, lp_data->lp, 1, &j, (char *)"U",
			  &ub);
   CPX_check_error("change_ub");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_lb(LPdata *lp_data, int j, double lb)
{
   cpx_status = CPXchgbds(lp_data->cpxenv, lp_data->lp, 1, &j, (char *)"L",
			  &lb);
   CPX_check_error("change_lb");
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void get_ub(LPdata *lp_data, int j, double *ub)
{
   cpx_status = CPXgetub(lp_data->cpxenv, lp_data->lp, ub, j, j);
   CPX_check_error("get_ub");
}

/*===========================================================================*/

void get_lb(LPdata *lp_data, int j, double *lb)
{
   cpx_status = CPXgetlb(lp_data->cpxenv, lp_data->lp, lb, j, j);
   CPX_check_error("get_lb");
}

/*===========================================================================*/

void get_objcoef(LPdata *lp_data, int j, double *objcoef)
{
   cpx_status = CPXgetobj(lp_data->cpxenv, lp_data->lp, objcoef, j, j);
   CPX_check_error("get_objcoef");
}

/*===========================================================================*/

void delete_rows(LPdata *lp_data, int deletable, int *free_rows)
{
   int i, m = lp_data->m, *bhead = lp_data->bhead;
   for (i = 0; i < m; ){
      if (bhead[i] < 0 && free_rows[-bhead[i]-1]){
	 bhead[i] = bhead[--m];
      }else{
	 i++;
      }
   }
   if (m + deletable != lp_data->m)
      lp_data->bhead_is_valid = FALSE;
   cpx_status = CPXdelsetrows(lp_data->cpxenv, lp_data->lp, free_rows);
   CPX_check_error("delete_rows");
   lp_data->m = m;
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;

}

/*===========================================================================*/

int delete_cols(LPdata *lp_data, int delnum, int *delstat)
{
   int i, *bhead = lp_data->bhead;
#if CPX_VERSION <= 600
   int j;
#endif

   /* this test can be made similar to the one in delete_cols if we pivot
      in the slacks corresponding to the colums to be deleted */
   for (i = lp_data->m-1; i >= 0; i--){
      if (bhead[i] >= 0)
	 if (delstat[bhead[i]]){
	    lp_data->bhead_is_valid = FALSE;
	    break;
	 }
   }
   cpx_status = CPXdelsetcols(lp_data->cpxenv, lp_data->lp, delstat);
   CPX_check_error("delete_cols - CPXdelsetcols");
   lp_data->nz = CPXgetnumnz(lp_data->cpxenv, lp_data->lp);
   CPX_check_error("delete_cols - CPXgetnumnz");
#if CPX_VERSION <= 600
   for (i = 0, j = 0; i < lp_data->n; i++){
      if (delstat[i])
	 delstat[i] = -1;
      else
	 delstat[i] = j++;
   }
#endif
   lp_data->n -= delnum;

   return(delnum);
}

/*===========================================================================*/

void release_var(LPdata *lp_data, int j, int where_to_move)
{
#if 0
   switch (where_to_move){
   case MOVE_TO_UB:
      lp_data->lpbas.cstat[j] = 2; break; /* non-basic at its upper bound */
   case MOVE_TO_LB:
      lp_data->lpbas.cstat[j] = 0; break; /* non-basic at its lower bound */
   }
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
#endif
}

/*===========================================================================*/

void free_row_set(LPdata *lp_data, int length, int *index)
{
   int i, j;
   constraint *rows = lp_data->rows;
   double *rhsval = lp_data->tmp.d; /* m */
   int *ind_e = lp_data->tmp.i1 + 2 * lp_data->m; /* m (now) */
   /* See comment in check_row_effectiveness why the shift! */
   char *sen_e = lp_data->tmp.c; /* m (now) */

   for (j=0, i=length-1; i>=0; i--){
      switch (rows[index[i]].cut->sense){
       case 'E': rhsval[i] = INFINITY; ind_e[j++] = index[i]; break;
       case 'L': rhsval[i] = INFINITY; break;
       case 'R':
       cpx_status = CPXchgcoef(lp_data->cpxenv, lp_data->lp, index[i], -2,
                               2*INFINITY);
       CPX_check_error("free_row_set - CPXchgcoef");
       case 'G': rhsval[i] = -INFINITY; break;
      }
   }
   cpx_status = CPXchgrhs(lp_data->cpxenv, lp_data->lp, length, index, rhsval);
   CPX_check_error("free_row_set - CPXchgrhs");
   if (j > 0){
      memset(sen_e, 'L', j);
      cpx_status = CPXchgsense(lp_data->cpxenv, lp_data->lp, j, ind_e, sen_e);
      CPX_check_error("free_row_set - CPXchgsense");
   }
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void constrain_row_set(LPdata *lp_data, int length, int *index)
{
   int i;
   constraint *rows = lp_data->rows;
   cut_data *cut;
   double *rhsval = lp_data->tmp.d; /* m (now) */
   char *sense = lp_data->tmp.c + lp_data->m; /* m (now) */
   char range_constraint = FALSE;

   for (i = length-1; i >= 0; i--){
      cut = rows[index[i]].cut;
      rhsval[i] = cut->rhs;
      if ((sense[i] = cut->sense) == 'R'){
	 range_constraint = TRUE;
      }
   }
   cpx_status = CPXchgrhs(lp_data->cpxenv, lp_data->lp, length, index, rhsval);
   CPX_check_error("constrain_row_set - CPXchgrhs");
   cpx_status=CPXchgsense(lp_data->cpxenv, lp_data->lp, length, index, sense);
   CPX_check_error("constrain_row_set - CPXchgsense");
   if (range_constraint){
      for (i = length-1; i >= 0; i--){
	 if (sense[i] == 'R'){
	    cpx_status = CPXchgcoef(lp_data->cpxenv,lp_data->lp, index[i], -2,
				    rows[index[i]].cut->range);
	    CPX_check_error("constrain_row_set - CPXchgcoef");
	 }
      }
   }
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void free_lp_solver_data(LPdata *lp_data, char arrays_too)
{
   if (lp_data){
      CPXfreeprob(lp_data->cpxenv, &(lp_data->lp));
      lp_data->lp = NULL;
      if (arrays_too){
	 FREE(lp_data->matbeg);
	 FREE(lp_data->matcnt);
	 FREE(lp_data->matind);
	 FREE(lp_data->matval);
	 FREE(lp_data->obj);
	 FREE(lp_data->rhs);
	 FREE(lp_data->rngval);
	 FREE(lp_data->sense);
	 FREE(lp_data->lb);
	 FREE(lp_data->ub);
	 FREE(lp_data->slacks);
      }
   }
}

/*===========================================================================*/

void write_mps(LPdata *lp_data, char *fname)
{
   cpx_status = CPXmpswrite(lp_data->cpxenv, lp_data->lp, fname);
   CPX_check_error("write_mps");
}

void write_sav(LPdata *lp_data, char *fname)
{
   cpx_status = CPXsavwrite(lp_data->cpxenv, lp_data->lp, fname);
   CPX_check_error("write_sav");
}

#endif /* CPLEX_xx */
