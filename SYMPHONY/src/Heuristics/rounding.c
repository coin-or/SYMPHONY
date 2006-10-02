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
 * TODO: rename this function
*/

int rnd_create_rnd_problem(lp_prob *p)
{

   rounding_problem *rp = (rounding_problem*) malloc(sizeof(rounding_problem));
   rnd_initialize_data(rp, p);

   /* create a search order of int variables */
   /* for now its just the lex. order of ALL int variables */
   int **var_list = (int **)malloc(sizeof(int*));
   int group_size = 10;
   int rn;
   int *var_search_grp = (int *)malloc(group_size*ISIZE);
   rnd_order(rp, var_list, rn);
   int grp_first_ind = 0;
   int grp_last_ind;
   int grp_cnt;
   int num_fixed_vars = 0;
   while(grp_first_ind<rn) {
      grp_cnt = group_size;
      grp_last_ind = grp_first_ind+group_size-1;
      if (num_fixed_vars+group_size>rn) {
	grp_last_ind = rn-1;
	grp_cnt = grp_last_ind-grp_first_ind+1;
      }
      memcpy(var_search_grp,var_list[grp_first_ind],grp_cnt*ISIZE);
      if (rnd_find_feas_rounding(rp, grp_cnt, var_search_grp, 0) == 
	    SYM_RND_FAIL) {
	 break;
      }
      grp_first_ind = grp_first_ind+grp_cnt;
      num_fixed_vars = num_fixed_vars+grp_cnt;
   }
   return 0;
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
   const double *colVal = rp->colVal;

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
   const double *colVal = rp->colVal;
   double *lb = rp->lb;
   double *ub = rp->ub;

   double *rowMin = rp->rowMin;
   double *rowMax = rp->rowMax;

   double xtmp = xval[xind]; /* old value */
   double lbtmp = lb[xind];
   double ubtmp = ub[xind];
   
   for (int i=colStart[xind];i<colStart[xind]+colLength[xind];i++) {
      int row = colIndex[i];
      double row_delta = (newval-xval[xind])*colVal[i]; 
      /* update rowActivity */
      rowActivity[row] = rowActivity[row]+row_delta;
      /* update rowMin and rowMax */
      if (colVal[i]<0) {
	 rowMin[row] = rowMin[row]+colVal[i]*(newub-ub[xind]);
	 rowMax[row] = rowMax[row]+colVal[i]*(newlb-lb[xind]);
      } else {
	 rowMin[row] = rowMin[row]+colVal[i]*(newlb-lb[xind]);
	 rowMax[row] = rowMax[row]+colVal[i]*(newub-ub[xind]);
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
   const double *colVal = rp->colVal;
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
	 if (colVal[index]>0) {
	    rowMax[row] = rowMax[row]+colVal[index]*ub[i];
	    rowMin[row] = rowMin[row]+colVal[index]*lb[i];
	 } else {
	    rowMax[row] = rowMax[row]+colVal[index]*lb[i];
	    rowMin[row] = rowMin[row]+colVal[index]*ub[i];
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


/*===========================================================================*/
int rnd_initialize_data(rounding_problem *rp, lp_prob *p)
{
   /* initialize the rounding problem data structure */
   LPdata *lp_data = p->lp_data;
   var_desc **vars = lp_data->vars;
   rp->n = lp_data->n;
   rp->m = lp_data->m;
   rp->num_ints = 0;
   rp->nz = lp_data->nz;

   rp->is_int = (int *)calloc(rp->n,ISIZE);
   rp->xval = (double *)malloc(rp->n*DSIZE);
   rp->obj = (double *)malloc(rp->n*DSIZE);
   rp->rowActivity = (double *)malloc(rp->m*DSIZE);
   rp->rowLower = (double *)malloc(rp->m*DSIZE);
   rp->rowUpper = (double *)malloc(rp->m*DSIZE);
   rp->rowMin = (double *)malloc(rp->m*DSIZE);
   rp->rowMax = (double *)malloc(rp->m*DSIZE);
   rp->rowStart = (int *)malloc((rp->m+1)*DSIZE);
   rp->rowIndex = (int *)malloc(rp->nz*ISIZE);
   rp->rowLength = (int *)malloc(rp->m*ISIZE);
   rp->rowVal = (double *)malloc(rp->nz*DSIZE);
   rp->colStart = (int *)malloc((rp->n+1)*ISIZE);
   rp->colIndex = (int *)malloc(rp->n*ISIZE);
   rp->colLength = (int *)malloc(rp->n*ISIZE);
   rp->colVal = (double *)malloc(rp->nz*DSIZE);
   rp->lb = (double *)malloc(rp->n*DSIZE);
   rp->ub = (double *)malloc(rp->n*DSIZE);
   rp->lpetol = lp_data->lpetol;

   /* query lp solver to retrieve bounds, rowActivities, etc. */
   get_bounds(lp_data);
   memcpy(rp->lb,lp_data->lb,rp->n*DSIZE);
   memcpy(rp->ub,lp_data->ub,rp->n*DSIZE);

   /* TODO: implement a better function in lp_solver.c */
   for (int j=0; j<rp->n; j++){
      get_column(lp_data, j, &rp->colVal[rp->colStart[j]], 
	    &rp->colIndex[rp->colStart[j]], &rp->colLength[j], &rp->obj[j]);
      rp->colStart[j+1] = rp->colStart[j] + rp->colLength[j];
   }

   /* TODO: implement better function in lp_solver.c */
   for (int i=0; i<rp->m; i++){
      get_row(lp_data, i, &rp->rowVal[rp->rowStart[i]],
	    &rp->rowIndex[rp->rowStart[i]], &rp->rowLength[i], &rp->rowUpper[i],
	    &rp->rowLower[i]);
      rp->rowStart[i+1] = rp->rowStart[i] + rp->rowLength[i];
   }

   for (int j=0; j<rp->n; j++) {
      /* search for integer variables */
      if (vars[j]->is_int) {
	 rp->is_int[j] = TRUE;
	 rp->is_int++;
      }
   }

   PRINT(p->par.verbosity,-1,("Rounding: number of integers = %d\n",rp->num_ints));

   return 0;
   /* Initialization complete */
}


/*===========================================================================*/
int rnd_order(rounding_problem *rp, int **var_list_ptr, int &rn)
{

   /* for now we include all ints in lex. order */
   int *var_list = *var_list_ptr;
   int order_type = 0;
   switch(order_type) {
    case 0:
      int count = 0;
      rn = rp->num_ints;
      var_list = (int *)malloc(rp->n*ISIZE);
      for (int j=0; j<rp->n; j++) {
	 if (rp->is_int[j]) {
	    var_list[count] = j;
	    count++;
	 }
      }
      if (count != rp->num_ints) {
	 printf("Rounding: Error in lexicographic ordering. Aborting.\n");
      }
      break;
   }
   return 0;
}

