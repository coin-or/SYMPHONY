/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "proccomm.h"
#include "qsortucb.h"
#include "lp.h"
#include "messages.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "BB_types.h"
#include "pack_cut.h"
#include "rounding.h"
#include "symphony_api.h"


/*===========================================================================*/
/* Test function
*/
int rnd_test()
{
   printf("Rounding: Successfully compiled.\n");
   return 0;
}

/*===========================================================================*/
/* 
 * take the lp_problem from symphony and try to construct a rounding problem.
 * A rounding problem can be constructed in several ways.
 * 1. Just take lexicographic ordering of all int variables.
 * 2. Take lexicographic ordering of all int variables with fractional values. 
*/

int rnd_create_rnd_problem(lp_prob *p)
{



}


/*===========================================================================*/
/* take a rounding problem and try to round the variables in var_search_grp
 * feasibiliy. we start rounding variables varnum, varnum+1 ... till we
 * completly scan the whole var_search_grp.
*/
int rnd_find_feas_rounding(rounding_problem *rp, int grp_cnt, int *var_search_grp, int varnum)
{
   int n = rp->n;
   int m = rp->m;
   int num_ints = rp->num_ints;
   int *int_vars = rp->int_vars;
   double *rowActivity = rp->rowActivity;
   const double *rowLower = rp->rowLower;
   const double *rowUpper = rp->rowUpper;
   double *xval = rp->xval;
   double lpetol = rp->lpetol;
   const int *rowStart = rp->rowStart;
   const int *rowLength = rp->rowLength;
   const double *matrixByCol = rp->matrixByCol;

   double *rowMin = rp->rowMin;
   double *rowMax = rp->rowMax;
   double *lb = rp->lb;
   double *ub = rp->ub;

   int xind = var_search_grp[varnum];
   
   /* find if already feasible or if beyond any repairs*/
   if (rnd_check_integrality(grp_cnt, var_search_grp, xval, lpetol) == TRUE && 
	 rnd_check_constr_feas(n, m, rowActivity, rowLower, rowUpper, lpetol) == 
	 TRUE) {
      return TRUE;
   } else if (rnd_check_if_feas_impossible(m, rowActivity, rowLower, rowUpper, rowMax, rowMin, lpetol)){
      return FALSE;
   } else if (varnum>=grp_cnt) {
      return FALSE;
   }

   double tmplb = lb[xind];
   double tmpub = ub[xind];
   double tmpx = xval[xind];
   int tmpint = FALSE;
   double fixedval;
   
   /* if integer, fix at this value. then round down and then round up */
   if (rp->is_int[xind]) {
      tmpint = TRUE;
      fixedval = floor(tmpx+0.5);
      if (fixedval<tmpub+lpetol && fixedval>tmplb-lpetol) {
	 rnd_update_rp(rp, xind, fixedval, fixedval, fixedval); 
	 if (rnd_find_feas_rounding(rp,grp_cnt,var_search_grp,varnum+1)==TRUE) {
	    return TRUE;
	 }
      }
   }
   
   /*round down */
   if (tmpint) {
      fixedval=floor(tmpx+0.5)-1;
   } else {
      fixedval = floor(tmpx);
   }
   if (fixedval<tmpub+lpetol && fixedval>tmplb-lpetol) {
      rnd_update_rp(rp, xind, fixedval, fixedval, fixedval); 
      if (rnd_find_feas_rounding(rp,grp_cnt,var_search_grp,varnum+1)==TRUE) {
	 return TRUE;
      }
   }

   /* round up */
   if (tmpint) {
      fixedval=floor(tmpx+0.5)+1;
   } else {
      fixedval = ceil(tmpx);
   }
   if (fixedval<tmpub+lpetol && fixedval>tmplb-lpetol) {
      rnd_update_rp(rp, xind, fixedval, fixedval, fixedval); 
      if (rnd_find_feas_rounding(rp,grp_cnt,var_search_grp,varnum+1)==TRUE) {
	 return TRUE;
      }
   }

   /* restore */
   rnd_update_rp(rp,xind, tmpx, tmplb, tmpub);
   return FALSE;

}

/*===========================================================================*/
/* update the data structure rp if one of the variable bounds are modified
*/
int rnd_update_rp(rounding_problem *rp, int xind, double newval, double newlb, double newub)
{
   int n = rp->n;
   int m = rp->m;
   int num_ints = rp->num_ints;
   int *int_vars = rp->int_vars;
   int *is_int = rp->is_int;
   double *rowActivity = rp->rowActivity;
   const double *rowLower = rp->rowLower;
   const double *rowUpper = rp->rowUpper;
   double *xval = rp->xval;
   double lpetol = rp->lpetol;
   const int *rowStart = rp->rowStart;
   const int *rowLength = rp->rowLength;
   const int *colLength = rp->colLength;
   const int *colStart = rp->colStart;
   const int *colIndex = rp->colIndex;
   const double *matrixByCol = rp->matrixByCol;
   double *lb = rp->lb;
   double *ub = rp->ub;

   double *rowMin = rp->rowMin;
   double *rowMax = rp->rowMax;

   double xtmp = xval[xind]; /* old value */
   double lbtmp = lb[xind];
   double ubtmp = ub[xind];
   
   for (int i=colStart[xind];i<colStart[xind]+colLength[xind];i++) {
      int row = colIndex[i];
      double row_delta = (newval-xval[xind])*matrixByCol[i]; 
      /* update rowActivity */
      rowActivity[row] = rowActivity[row]+row_delta;
      /* update rowMin and rowMax */
      if (matrixByCol[i]<0) {
	 rowMin[row] = rowMin[row]+matrixByCol[i]*(newub-ub[xind]);
	 rowMax[row] = rowMax[row]+matrixByCol[i]*(newlb-lb[xind]);
      } else {
	 rowMin[row] = rowMin[row]+matrixByCol[i]*(newlb-lb[xind]);
	 rowMax[row] = rowMax[row]+matrixByCol[i]*(newub-ub[xind]);
      }
   }

   xval[xind] = newval;
   lb[xind] = newlb;
   ub[xind] = newub;
}

/*===========================================================================*/
int rnd_check_if_feas_impossible(int m, double *rowActivity, const double *rowLower, const double *rowUpper, double *rowMax, double *rowMin, double lpetol)
{
   for (int i=0; i<m; i++) {
      if (rowMax[i]<rowLower[i]-lpetol || rowMin[i]>rowUpper[i]+lpetol) {
	 return TRUE;
      }
   }
   return FALSE;
}


/*===========================================================================*/
int rnd_find_row_bounds(rounding_problem *rp)
{
   int n = rp->n;
   int m = rp->m;
   const int *colLength = rp->colLength;
   const int *colStart = rp->colStart;
   const int *colIndex = rp->colIndex;
   const double *matrixByCol = rp->matrixByCol;
   double *lb = rp->lb;
   double *ub = rp->ub;

   double *rowMin = rp->rowMin;
   double *rowMax = rp->rowMax;

   double tmp_infinity = 10e20;

   for (int i=0; i<n; i++) {
      if (lb[i]<-tmp_infinity || ub[i]>tmp_infinity) {
	 printf("in trouble: lb[i] = %g\t ub[i] = %g\n",lb[i],ub[i]);
      }
      for (int index=colStart[i]; index<colStart[i]+colLength[i]; index++) {
	 int row = colIndex[index];
	 if (matrixByCol[index]>0) {
	    rowMax[row] = rowMax[row]+matrixByCol[index]*ub[i];
	    rowMin[row] = rowMin[row]+matrixByCol[index]*lb[i];
	 } else {
	    rowMax[row] = rowMax[row]+matrixByCol[index]*lb[i];
	    rowMin[row] = rowMin[row]+matrixByCol[index]*ub[i];
	 }
      }
   }
}

/*===========================================================================*/
int rnd_check_constr_feas(int n, int m, double *rowActivity, const double 
      *rowLower, const double *rowUpper, double lpetol)
{
   return 0;
}

/*===========================================================================*/
int rnd_check_integrality(int n, int *intvars, double *xval, double lpetol) 
{

}

