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

#include "sym_proccomm.h"
#include "qsortucb.h"
#include "sym_lp.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_pack_cut.h"
#include "sym_rounding.h"
#include "symphony.h"


/*===========================================================================*/
/* Test function
*/
int rnd_test(lp_prob *p)
{
   printf("Rounding: Successfully compiled.\n");
   rnd_create_rnd_problem(p);
   exit(0);
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
   rnd_find_row_bounds(rp);
   /* create a search order of int variables */
   /* for now its just the lex. order of ALL int variables */
   int **var_list = (int **)malloc(sizeof(int*));
   int group_size = 10;
   int rn;
   int *var_search_grp = (int *)malloc(group_size*ISIZE);
   printf("Entering rnd_order\n");
   sym_rnd_order(rp, var_list, &rn);
   printf("Exiting rnd_orderaslfjasdflksdjflaksdjflkajk\n");
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
      memcpy(var_search_grp,&(*var_list)[grp_first_ind],grp_cnt*ISIZE);
      printf("FOOBARSDLI:FJ\n");
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
   int verbosity = rp->verbosity;

   int xind = var_search_grp[varnum];
   printf("rounding: inside find_feas_rounding\n"); 
   /* find if already feasible or if beyond any repairs*/
   if (rnd_is_integral(grp_cnt, var_search_grp, xval, lpetol, verbosity) == TRUE
	 && rnd_is_constr_feas(n, m, rowActivity, rowLower, rowUpper, lpetol, 
	    verbosity) == TRUE) {
      return TRUE;
   } else if (rnd_if_feas_impossible(m, rowActivity, rowLower, rowUpper, rowMax, rowMin, lpetol, verbosity)){
      return FALSE;
   } else if (varnum>=grp_cnt) {
      return FALSE;
   }

   double tmplb = lb[xind];
   double tmpub = ub[xind];
   double tmpx = xval[xind];
   int tmpint = FALSE;
   double fixedval;
   
   /* if already integer, fix at this value. 
    * and then round down and then round up */
   if (fabs(tmpx-floor(tmpx+0.5))<lpetol) {
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
int rnd_if_feas_impossible(int m, double *rowActivity, const double *rowLower, const double *rowUpper, double *rowMax, double *rowMin, double lpetol, int verbosity)
{
   for (int i=0; i<m; i++) {
      printf("%f\t%f\t%f\t%f\n",rowLower[i], rowUpper[i], rowMax[i], rowMin[i]);
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
   const double *rowVal = rp->rowVal;
   const int *rowIndex = rp->rowIndex;
   const int *rowStart = rp->rowStart;
   const int *rowLength = rp->rowLength;
   double *lb = rp->lb;
   double *ub = rp->ub;

   double *rowMin = rp->rowMin;
   double *rowMax = rp->rowMax;
   double lpetol = rp->lpetol;

   double tmp_infinity = 10e20;

   for (int row=0; row<m; row++) {
      rowMax[row] = rowMin[row] = 0;
      for (int index=rowStart[row]; index<rowStart[row]+rowLength[row]; 
	    index++) {
	 int col = rowIndex[index];
	 if (rowVal[index]>lpetol) {
	    if (rowMax[row]<tmp_infinity && ub[col]<tmp_infinity) {
	       rowMax[row] = rowMax[row]+rowVal[index]*ub[col];
	    }
	    if (rowMin[row]>-tmp_infinity && lb[col]>-tmp_infinity) {
	       rowMin[row] = rowMin[row]+colVal[index]*lb[col];
	    }
	 } else if (colVal[index]<-lpetol) {
	    if (rowMax[row]<tmp_infinity && lb[col]>-tmp_infinity) {
	       rowMax[row] = rowMax[row]+colVal[index]*lb[col];
	    }
	    if (rowMin[row]>-tmp_infinity && ub[col]<tmp_infinity) {
	       rowMin[row] = rowMin[row]+colVal[index]*ub[col];
	    }
	 }
      }
   }
}

/*===========================================================================*/
int rnd_is_constr_feas(int n, int m, double *rowActivity, const double 
      *rowLower, const double *rowUpper, double lpetol, int verbosity)
{
   for (int row=0; row<m; row++) {
      if (rowActivity[row]>rowUpper[row]+lpetol || 
	    rowActivity[row]<rowLower[row]-lpetol) {
	 PRINT(verbosity,-1,("Rounding: constraint %d has lb %f, ub %f, activity %f\n.",row, rowLower[row], rowUpper[row], rowActivity[row]));
	 return FALSE;
      }
   }
   return TRUE;
}

/*===========================================================================*/
int rnd_is_integral(int n, int *intvars, double *xval, double lpetol, int verbosity) 
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
   rp->colIndex = (int *)malloc(rp->nz*ISIZE);
   rp->colLength = (int *)malloc(rp->n*ISIZE);
   rp->colVal = (double *)malloc(rp->nz*DSIZE);
   rp->lb = (double *)malloc(rp->n*DSIZE);
   rp->ub = (double *)malloc(rp->n*DSIZE);
   rp->lpetol = lp_data->lpetol;
   rp->verbosity = p->par.verbosity;

   printf("Hellow orld\n");
   /* query lp solver to retrieve bounds, rowActivities, etc. */
   get_bounds(lp_data);
   memcpy(rp->lb,lp_data->lb,rp->n*DSIZE);
   memcpy(rp->ub,lp_data->ub,rp->n*DSIZE);
	
   rp->colStart[0] = 0;
   /* TODO: implement a better function in lp_solver.c */
   for (int j=0; j<rp->n; j++){
      printf ("j = %d, colStart = %d\n",j,rp->colStart[j]);
      get_column(lp_data, j, &rp->colVal[rp->colStart[j]], 
	    &rp->colIndex[rp->colStart[j]], &rp->colLength[j], &rp->obj[j]);
      rp->colStart[j+1] = rp->colStart[j] + rp->colLength[j];
   }

   /* TODO: implement better function in lp_solver.c */
   rp->rowStart[0] = 0;
   for (int i=0; i<rp->m; i++){
      printf ("i = %d\n",i);
      get_row(lp_data, i, &rp->rowVal[rp->rowStart[i]],
	    &rp->rowIndex[rp->rowStart[i]], &rp->rowLength[i], &rp->rowUpper[i],
	    &rp->rowLower[i]);
      rp->rowStart[i+1] = rp->rowStart[i] + rp->rowLength[i];
   }

   for (int j=0; j<rp->n; j++) {
      /* search for integer variables */
      if (vars[j]->is_int) {
	 rp->is_int[j] = TRUE;
	 rp->num_ints++;
      }
   }

   PRINT(p->par.verbosity,-1,("Rounding: number of integer variables = %d\n",rp->num_ints));
   printf("foobar\n");
   return 0;
   /* Initialization complete */
}


/*===========================================================================*/
int sym_rnd_order(rounding_problem *rp, int **var_list, int *rn)
{

   /* for now we include all ints in lex. order */
   printf("rounding: inside order\nfoasdfasdf\n");
   int order_type = 0;
   switch(order_type) {
    case 0:
      int count = 0;
      *rn = rp->num_ints;
      int *vlist = (int *)malloc(rp->n*ISIZE);
      printf ("SDLJF\n");
      for (int j=0; j<rp->n; j++) {
	 if (rp->is_int[j]) {
	    vlist[count] = j;
	    count++;
	    PRINT(rp->verbosity,-1,("Rounding: variable %d added to list\n",j));
	 }
      }

      *var_list = vlist;
      if (count != rp->num_ints) {
	 printf("Rounding: Error in lexicographic ordering. Aborting.\n");
      }
      break;
   }
   return 0;
}

