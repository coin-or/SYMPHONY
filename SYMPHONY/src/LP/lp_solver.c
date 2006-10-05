/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* The OSI interface in this file was written by Menal Guzelsoy.             */
/* The OSL interface was written by Ondrej Medek.                            */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>              /* free() is here on AIX ... */
#include <math.h>
#include <string.h>

#include "lp_solver.h"
#include "BB_constants.h"
#include "BB_macros.h"


/*===========================================================================*/

/*===========================================================================*\
 * This file contains the interface with the LP Solver.
 * The first few routines are independent of what LP solver is being used.
\*===========================================================================*/

double dot_product(double *val, int *ind, int collen, double *col)
{
   const int* lastind = ind + collen;
   double prod = 0;
   while (ind != lastind)
      prod += (*val++) * col[*ind++];
   return(prod);
}

/*===========================================================================*/

void free_lp_arrays(LPdata *lp_data)
{
   FREE(lp_data->not_fixed);
   FREE(lp_data->status);
   FREE(lp_data->x);
   FREE(lp_data->dj);
   FREE(lp_data->dualsol);
   FREE(lp_data->slacks);
#ifdef __CPLEX__
   FREE(lp_data->lb);
   FREE(lp_data->ub);
#endif
   FREE(lp_data->vars);
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

void free_mip_desc(MIPdesc *mip)
{
   int j;
   
   FREE(mip->matbeg);
   FREE(mip->matind);
   FREE(mip->matval);
#if 0
   FREE(mip->col_lengths);
   FREE(mip->row_matbeg);
   FREE(mip->row_matind);
   FREE(mip->row_matval);
   FREE(mip->row_lengths);
#endif
   FREE(mip->obj);
   FREE(mip->obj1);
   FREE(mip->obj2);
   FREE(mip->rhs);
   FREE(mip->rngval);
   FREE(mip->sense);
   FREE(mip->lb);
   FREE(mip->ub);
   FREE(mip->is_int);
   if (mip->colname){
      for (j = 0; j < mip->n; j++){
	 FREE(mip->colname[j]);
      }
      FREE(mip->colname);
   }
}

/*===========================================================================*/

void size_lp_arrays(LPdata *lp_data, char do_realloc, char set_max,
		    int row_num, int col_num, int nzcnt)
{
   char resize_m = FALSE;
   char resize_n = FALSE;
   int maxm, maxn, maxnz, maxmax;

   if (set_max){
      maxm = row_num;
      maxn = col_num;
      maxnz = nzcnt;
   }else{
      maxm = lp_data->m + row_num;
      maxn = lp_data->n + col_num;
      maxnz = lp_data->nz + nzcnt;
   }

   if (maxm > lp_data->maxm){
      resize_m = TRUE;
      lp_data->maxm = maxm + (set_max ? 0 : BB_BUNCH);
      if (! do_realloc){
         FREE(lp_data->dualsol);
         lp_data->dualsol = (double *) malloc(lp_data->maxm * DSIZE);
	 FREE(lp_data->slacks);
	 lp_data->slacks  = (double *) malloc(lp_data->maxm * DSIZE);
     }else{
         lp_data->dualsol = (double *) realloc((char *)lp_data->dualsol,
                                               lp_data->maxm * DSIZE);
	 lp_data->slacks  = (double *) realloc((void *)lp_data->slacks,
					       lp_data->maxm * DSIZE);
      }
      /* rows is realloc'd in either case just to keep the base constr */
      lp_data->rows = (row_data *) realloc((char *)lp_data->rows,
                                             lp_data->maxm*sizeof(row_data));
   }
   if (maxn > lp_data->maxn){
      int oldmaxn = MAX(lp_data->maxn, lp_data->n);
      resize_n = TRUE;
      lp_data->maxn = maxn + (set_max ? 0 : 5 * BB_BUNCH);
      if (! do_realloc){
         FREE(lp_data->x);
         lp_data->x = (double *) malloc(lp_data->maxn * DSIZE);
         FREE(lp_data->dj);
         lp_data->dj = (double *) malloc(lp_data->maxn * DSIZE);
         FREE(lp_data->status);
         lp_data->status = (char *) malloc(lp_data->maxn * CSIZE);
#ifdef __CPLEX__
	 FREE(lp_data->lb);
	 lp_data->lb = (double *) malloc(lp_data->maxn * DSIZE);
	 FREE(lp_data->ub);
	 lp_data->ub = (double *) malloc(lp_data->maxn * DSIZE);
#endif
      }else{
         lp_data->x = (double *) realloc((char *)lp_data->x,
                                         lp_data->maxn * DSIZE);
         lp_data->dj = (double *) realloc((char *)lp_data->dj,
                                          lp_data->maxn * DSIZE);
         lp_data->status = (char *) realloc((char *)lp_data->status,
                                            lp_data->maxn * CSIZE);
#ifdef __CPLEX__
	 lp_data->lb = (double *) realloc((char *)lp_data->lb,
					  lp_data->maxn * DSIZE);
	 lp_data->ub = (double *) realloc((char *)lp_data->ub,
					  lp_data->maxn * DSIZE);
#endif
      }
   }
   if (maxnz > lp_data->maxnz){
      lp_data->maxnz = maxnz + (set_max ? 0 : 20 * BB_BUNCH);
   }

   /* re(m)alloc the tmp arrays */
   if (resize_m || resize_n){
      temporary *tmp = &lp_data->tmp;
      maxm = lp_data->maxm;
      maxn = lp_data->maxn;
      maxmax = MAX(maxm, maxn);
      /* anything with maxm and maxn in it has to be resized */
      FREE(tmp->c);
      FREE(tmp->i1);
      FREE(tmp->d);
      tmp->c = (char *) malloc(CSIZE * maxmax);
      tmp->i1 = (int *) malloc(ISIZE * MAX(3*maxm, 2*maxn + 1));
      tmp->d = (double *) malloc(DSIZE * 2 * maxmax);
      /* These have to be resized only if maxm changes */
      if (resize_m){
	 FREE(tmp->i2);
	 FREE(tmp->p1);
	 FREE(tmp->p2);
	 tmp->i2 = (int *) malloc(maxm * ISIZE);
	 tmp->p1 = (void **) malloc(maxm * sizeof(void *));
	 tmp->p2 = (void **) malloc(maxm * sizeof(void *));
      }
   }
}

#ifdef USE_GLPMPL

/*This function reads in the GNU MathProg model file and returns either 1 if 
  it is succeded or 0 otherwise.*/

int read_gmpl(MIPdesc *mip, char *modelfile, char *datafile, char *probname)
{
   MPL * mpl;         
   int errors;
   int i, j, k, obj_index, length, type, count, nonzeros;
   double *matval;
   int *matind;
   int * matbeg;

   double *row_lb;
   double *row_ub;

   int *indices;
   double *values;   
   double inf = SYM_INFINITY;

   mpl = mpl_initialize();  /* initialize the translator */
    
   //   errors = mpl_read_model(mpl, modelfile);   /* read the model file */

   //version > 4.1
   errors = mpl_read_model(mpl, modelfile, FALSE); 
 
   /*if the data is not in the model file and will be given seperately, 
     then errors=1! */
   if (errors == 1){  
      errors = mpl_read_data(mpl, datafile);
   }

   /*if there are some errors in reading the file(s): then errors!=2*/ 
   if (errors != 2){
      printf("\nError in reading the model (or data) file!"); 
      mpl_terminate(mpl);   /* free all the mpl related stuff */
      return(0);
   }

   /*Generate the model variables, constraints, objective by storing in the 
     translator database.It is possible to capture the messages in a file by 
     passing the filename instead of NULL.*/

   errors = mpl_generate(mpl, NULL);  
   if (errors != 3){           
      printf("\nError in generating the model!");  
      mpl_terminate(mpl);
      return(0);
   }
   
   strncpy(probname, mpl_get_prob_name(mpl), 80); /* name the problem */

   /* get num of rows and cols */
   mip->m  = mpl_get_num_rows(mpl)-1; /* subtract the objective row */
   mip->n  = mpl_get_num_cols(mpl);
   mip->nz = 0; /* for now... */

   /*Indices and values of nonzeros will return beginning with indices[1] and
     values[1]. Also note that row and column indices begin with 1 in glpmpl*/

   /*get mip->nz and mip->obj*/
   mip->obj    = (double *) calloc(DSIZE, mip->n);
   mip->obj1   = (double *) calloc(DSIZE, mip->n);
   mip->obj2   = (double *) calloc(DSIZE, mip->n);
   
   indices = (int *) malloc(ISIZE * (mip->n + 1));
   values = (double *) malloc(DSIZE * (mip->n + 1));

   count = 0;

   for(i = 0; i < mip->m + 1; i++){

      type = mpl_get_row_kind(mpl, i+1);
      if (type == MPL_ST){  /* constraints */
	 /* mpl_get_mat_row returns the # of nonzeros in the row i+1. */
	 mip->obj_sense = SYM_MINIMIZE; 
	 mip->nz += mpl_get_mat_row(mpl, i+1, NULL, NULL); 
      }else{
	 obj_index = i;
	 mip->obj_offset = mpl_get_row_c0(mpl, i+1);
	 length = mpl_get_mat_row(mpl, i+1, indices, values);
	 if (type == MPL_MAX){
	    mip->obj_sense = SYM_MAXIMIZE; 

	    /* this is done in create_subproblem_u function! */
	    //	    for (j = 1; j <= length; j++){  
	    //   mip->obj[indices[j]-1] = -values[j];
	    //	    }
	 }else{   /* type == MPL_MIN */
	    for( j = 1; j <= length; j++ ) 
	       mip->obj[indices[j]-1] = values[j]; /* assign the obj coeff. */
	 }
	 count++;
	 if (count > 1){
	    printf("\nError due to having more than one objective function!"); 
	    printf("\nSYMPHONY can solve problems with just one "); 
	    printf("objective function\n");
	    return(0);
	 }  
      }      
   }

   if (count < 1){
      printf("\nError in reading the objective function!"); 
      printf("\nSYMPHONY requires an objective function!");
      return(0);
   }

   /* Define a row ordered dummy constraint matrix since glpmpl returns the 
      constraint definitions as row ordered, we will change its order later. */

   /* fill the dummy matbeg, matind, matval, row_lb and row_ub arrays */
   matbeg = (int *) malloc(ISIZE * (mip->m + 1));
   matind = (int *) malloc(ISIZE * mip->nz);
   matval = (double *) malloc(DSIZE * mip->nz);

   row_ub = (double *) malloc(DSIZE * mip->m);
   row_lb = (double *) malloc(DSIZE * mip->m);

   matbeg[0] = 0;
   nonzeros = 0;
   for(i = 0, k = 0; i < mip->m+1; i++){
      if(i != obj_index){
	 /* read the nonzeros in row i+1 */     
	 length = mpl_get_mat_row(mpl, i+1, indices, values); 
	 /* get the row bounds. we use k instead of i since we have the obj 
	    row somewhere. */
	 type = mpl_get_row_bnds(mpl, i+1, &row_lb[k],&row_ub[k]); 
	 switch(type)                                           
	    {
	    case MPL_FR:  /* free */
	       row_lb[k] = -inf;
	       row_ub[k] =  inf;
	       break;
	    case MPL_LO:  /* has lower bound */
	       row_ub[k] =  inf;
	       break;	 
	    case MPL_UP:  /* has upper bound */
	       row_lb[k] = -inf;
	       break;
	    default: /* is bounded from both sides or is an equality */
	       break;
	    }
	 for (j = 0; j < length; j++){
	    matind[matbeg[k]+j] = indices[j+1] - 1;
	    matval[matbeg[k]+j] = values[j+1];
	 }
	 nonzeros += length;
	 k++;
	 matbeg[k] = nonzeros;
      }
   }

   /* fill the column related definitions: ub, lb, is_int and colname arrays */

   mip->ub      = (double *) malloc(DSIZE * mip->n);
   mip->lb      = (double *) malloc(DSIZE * mip->n);
   mip->is_int  = (char *)   calloc(CSIZE, mip->n);
   mip->colname = (char **)  malloc(sizeof(char *) * mip->n);   

   for (j = 0; j < mip->n; j++){
      type = mpl_get_col_bnds(mpl, j+1, &mip->lb[j], &mip->ub[j]);
      switch(type){
	case  MPL_FR: /* free */
	   mip->lb[j] = -inf;
	   mip->ub[j] =  inf;
	   break;
	case MPL_LO:  /* has lower bound */
	    mip->ub[j] =  inf;
	    break;
	case MPL_UP:  /* has upper bound */
	    mip->lb[j] = -inf;
	    break;
      default:  /* has both lower and upper bound or is a fixed variable */
	    break;
	 }

      type = mpl_get_col_kind(mpl, j+1);
      if(type == MPL_INT || type == MPL_BIN){
	 mip->is_int[j] = TRUE;
      }  
      /* bounds for binary variables were probably not assigned.So assign them! */
      if(type == MPL_BIN){ 
	 mip->ub[j] = 1.0;	 
	 mip->lb[j] = 0.0;

      }

      mip->colname[j] = (char *) malloc(CSIZE * 255); 
      strncpy(mip->colname[j], mpl_get_col_name(mpl,j+1),255);
      mip->colname[j][254] = 0;  /* ??? */
   }

   /*load the definitions to a CoinPackedMatrix as row ordered and get the 
     column ordered matrix after reversing its order in order to fill the 
     matrix definitons as column ordered*/

   mip->matbeg = (int *)    calloc(ISIZE, (mip->n + 1));
   mip->matval = (double *) malloc(DSIZE * mip->nz);
   mip->matind = (int *)    malloc(ISIZE * mip->nz);

#if 0
   //make CoinPackedMatrix help us for now!!!
   CoinPackedMatrix matrixByCol (false, mip->n, 
			   mip->m, mip->nz, matval, matind, matbeg, 0);
   matrixByCol.reverseOrdering();


   memcpy(mip->matbeg, const_cast<int *>(matrixByCol.getVectorStarts()),
	  ISIZE * (mip->n + 1));   
   memcpy(mip->matval, const_cast<double *> (matrixByCol.getElements()),
	  DSIZE * mip->nz);  
   memcpy(mip->matind, const_cast<int *> (matrixByCol.getIndices()), 
	  ISIZE * mip->nz);  
#endif

   /* what if the user doesn't have COIN, is that possible?:) */
   nonzeros = 0;
   for(j = 0; j < mip->n; j++){
      for(i = 0; i < mip->m; i++){
	 for(k = matbeg[i]; k < matbeg[i+1]; k++){
	    if(matind[k] == j){	   
	       mip->matind[nonzeros] = i;
	       mip->matval[nonzeros] = matval[k];
	       nonzeros++;	  
	       break;
	    }
	 }
      } 
      mip->matbeg[j+1] = nonzeros;
   }

   /*get the other definitions: rhs, sense and rngval from row_lb and row_ub*/

   mip->rhs    = (double *) malloc(DSIZE * mip->m);
   mip->sense  = (char *)   malloc(CSIZE * mip->m);
   mip->rngval = (double *) malloc(DSIZE * mip->m);

   /* convertBoundToSense: stolen from COIN :) */
   for(i = 0; i < mip->m; i++) {      
      mip->rngval[i] = 0.0;
      if (row_lb[i] > -inf) {
	 if (row_ub[i] < inf) {
	    mip->rhs[i] = row_ub[i];
	    if (row_lb[i] == row_ub[i]) {
	       mip->sense[i] = 'E';
	    } else {
	       mip->sense[i] = 'R';
	       mip->rngval[i] = row_ub[i] - row_lb[i];
	    }
	 }else{
	    mip->sense[i] = 'G';
	    mip->rhs[i] = row_lb[i];
	 }
      }else{
	 if (row_ub[i] < inf) {
	    mip->sense[i] = 'L';
	    mip->rhs[i] = row_ub[i];
	 }else{
	    mip->sense[i] = 'N';
	    mip->rhs[i] = 0.0;
	 }
      }
   }

   FREE(matind);
   FREE(matval);
   FREE(matbeg);
   FREE(row_lb);
   FREE(row_ub);
   
   /* if you could reach here by chance, then you are safe anymore:) */
   return(1);
}
#endif /* USE_GLPMPL */

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
 * - no fastmip is used
 * - no scaling is used - cannot test it
 * - possible problems with getting reduced costs
 * - LPdata->tmp field mostly not used for safe. malloc and free for temporary
 *   fields are used instead => slow down
 *============================================================================
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
  if (lp_data->lp != NULL) {
    osllib_status = ekk_deleteModel(lp_data->lp);
    OSL_check_error("close_lp_solver - ekk_deleteModel");
    lp_data->lp = NULL;
  }
  ekk_endContext(lp_data->env);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function loads the data of an lp into the lp solver. 
\*===========================================================================*/
void load_lp_prob(LPdata *lp_data, int scaling, int fastmip)
{
   int i;
   double *lr = lp_data->tmp.d, *ur = lp_data->tmp.d + lp_data->n;
   
   lp_data->lp = ekk_newModel(lp_data->env, NULL);
   osllib_status = (lp_data->env == NULL);
   OSL_check_error("load_lp_prob - ekk_newModel");
   
   for (i = 0; i < lp_data->m; i++) {
      switch (lp_data->mip->sense[i]) {
       case 'E': lr[i] = ur[i] = lp_data->mip->rhs[i]; break;
       case 'L': lr[i] = - OSL_INFINITY; ur[i] = lp_data->mip->rhs[i]; break;
       case 'G': lr[i] = lp_data->mip->rhs[i]; ur[i] = OSL_INFINITY; break;
       case 'R':
	 if (lp_data->mip->rngval[i] >= 0) {
	    ur[i] = lp_data->mip->rhs[i]; 
	    lr[i] = ur[i] - lp_data->mip->rngval[i];
	 } else {
	    ur[i] = lp_data->mip->rhs[i]; 
	    lr[i] = ur[i] + lp_data->mip->rngval[i];
	 }
	 break;
       default: /* This should never happen ... */
	 osllib_status = -1;
	 OSL_check_error("load_lp - unknown sense");
      }
   }
   osllib_status =
      ekk_loadRimModel(lp_data->lp, lp_data->m, lr, ur, lp_data->n, 
		       lp_data->mip->obj, lp_data->mip->lb,
		       lp_data->mip->ub);
   OSL_check_error("load_lp - ekk_loadRimModel");
   osllib_status =
      ekk_addColumnElementBlock(lp_data->lp, lp_data->n, lp_data->mip->matind,
				lp_data->mip->matbeg, lp_data->mip->matval);
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

   lp_data->m = lp_data->n = lp_data->nz = 0;
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

   lr = (double *) malloc(rcnt * DSIZE);
   ur = (double *) malloc(rcnt * DSIZE);
   for (i = rcnt - 1; i >= 0; i--) {
      switch (sense[i]) {
       case 'E': lr[i] = ur[i] = rhs[i]; break;
       case 'L': lr[i] = - OSL_INFINITY; ur[i] = rhs[i]; break;
       case 'G': lr[i] = rhs[i]; ur[i] = OSL_INFINITY; break;
       case 'R': lr[i] = ur[i] = lp_data->mip->rhs[i]; break;
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

   ekk_mergeBlocks(lp_data->lp, 1);
   ekk_setIiternum(lp_data->lp, 0);

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
   
#if 1
   ekk_dualSimplex(lp_data->lp); 
#else
   ekk_simplex(lp_data->lp, 256 + 32); // no presolve and no scaling
#endif

   term = ekk_getIprobstat(lp_data->lp);
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
    case 0:
      term = LP_OPTIMAL;
      break;
    case 1:
      term = LP_D_UNBOUNDED;
      ekk_infeasibilities(lp_data->lp, 1, 1, NULL, NULL);
      break;
    case 2:
      term = LP_D_INFEASIBLE;
      break;
    case 3:
      term = LP_D_ITLIM;
      break;
    case 4:
      osllib_status = -1;
      OSL_check_error("osllib_status-ekk_dualSimplex found no solution!");
      term = LP_ABANDONED;
      break;
    case 5:
      LP_D_OBJLIM;
      break;
    case 6:
      osllib_status = -1;
      OSL_check_error("osllib_status-ekk_dualSimplex lack of dstorage file"
			 "space!");
      term = LP_ABANDONED;
      break;
    default: term = LP_ABANDONED;
      break;
   }
   
   lp_data->termcode = term;
   
   if (term != LP_ABANDONED){
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
\*===========================================================================*/
void set_obj_upper_lim(LPdata *lp_data, double lim)
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

void get_column(LPdata *lp_data, int j,
		double *colval, int *colind, int *collen, double *cj)
{
   EKKVector vec;
   vec = ekk_getColumn(lp_data->lp, j);
   *collen = vec.numNonZero;
   memcpy(colind, vec.index, *collen * ISIZE);
   memcpy(colval, vec.element, *collen * DSIZE);
   ekk_freeVector(&vec);
   get_objcoef(lp_data, j, cj);
}

/*===========================================================================*/
void get_row(LPdata *lp_data, int i,
	     double *rowval, int *rowind, int *rowlen,
	     double *rowub, double *rowlb)
{
   EKKVector vec;
   vec = ekk_getRow(lp_data->lp, i);
   *rowlen = vec.numNonZero;
   memcpy(rowind, vec.index, *rowlen * ISIZE);
   memcpy(rowval, vec.element, *rowlen * DSIZE);
   ekk_freeVector(&vec);
}

/*===========================================================================*/
/* This routine returns the index of a row which proves the lp to be primal
 * infeasible. It is only needed when column generation is used.             */
/*===========================================================================*/

int get_proof_of_infeas(LPdata *lp_data, int *infind)
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
   memcpy(lp_data->x, ekk_colsol(lp_data->lp), lp_data->n * DSIZE);
}

/*===========================================================================*/
void get_dj_pi(LPdata *lp_data)
{
   /*If scaling, fast integer or compress is used, maybe some changes will be
     needed */
   /* OSL returns changed sign - is it good or not? */
   memcpy(lp_data->dualsol, ekk_rowduals(lp_data->lp), lp_data->m * DSIZE);

# if 0
   /* changing the sign */
   for (i = lp_data->m - 1; i >= 0; i --) {
      lp_data->dualsol[i] = - lp_data->dualsol[i];
   }
#endif
   
   memcpy(lp_data->dj, ekk_colrcosts(lp_data->lp), lp_data->n * DSIZE);

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
   row_data *rows = lp_data->rows;
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
}

/*===========================================================================*/
/* This function is never called ...*/

void change_sense(LPdata *lp_data, int cnt, int *index, char *sense)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
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

void get_bounds(LPdata *lp_data)
{
   lp_data->ub = const_cast<double *>(ekk_colupper(lp_data->lp));
   lp_data->lb = const_cast<double *>(ekk_collower(lp_data->lp));
}

/*===========================================================================*/

void get_objcoef(LPdata *lp_data, int j, double *objcoef)
{
   /* Maybe some range checking could be added ...*/
   const double *oc = ekk_objective(lp_data->lp);
   *objcoef = oc[j];
}

/*===========================================================================*/

void delete_rows(LPdata *lp_data, int deletable, int *free_rows)
{
   int i, m = lp_data->m;
   int *which = lp_data->tmp.i1 + lp_data->m;
   int delnum = 0;

   /* which = calloc(delnum, ISIZE); */
   for (i = m - 1, delnum = 0; i >= 0; i--){
      if (free_rows[i]){
	 which[delnum++] = i;
      }
   }
   osllib_status = ekk_deleteRows(lp_data->lp, delnum, which);
   OSL_check_error("delete_rows - ekk_deleteRows");

#if 0
   /* Make result as CPLEX does*/
   for (i = 0, delnum = 0; i < m; i++){
      if (free_rows[i])
	 free_rows[i] = -1;
      else
	 free_rows[i] = delnum++;
   }
#endif

   lp_data->m -= delnum;
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

int delete_cols(LPdata *lp_data, int delnum, int *delstat)
{
   int i, n = lp_data->n;
   int *which = (int *) calloc(delnum, ISIZE);
   int num_to_delete = 0, num_to_keep = 0;
   double *dj = lp_data->dj;
   double *x = lp_data->x;
   char *status = lp_data->status;

   for (i = n - 1, num_to_delete = 0; i >= 0; i--) {
      if (delstat[i]) {
	 which[num_to_delete++] = i;
      }
   }

   if (!num_to_delete) return(0);

   osllib_status = ekk_deleteColumns(lp_data->lp, num_to_delete, which);
   OSL_check_error("delete_cols - ekk_deleteCols");
   FREE(which);
   
   lp_data->nz = ekk_getInumels(lp_data->lp);
   OSL_check_error("delete_cols - ekk_getInumels");
   
   for (i = 0, num_to_keep = 0; i < lp_data->n; i++){
      if (delstat[i]){
	 delstat[i] = -1;
      }else{
	 delstat[i] = num_to_keep++;
	 dj[delstat[i]] = dj[i];
	 x[delstat[i]] = x[i];
	 status[delstat[i]] = status[i];
      }
   }
   
   lp_data->n = num_to_keep;
   return (num_to_delete);
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
   row_data *rows = lp_data->rows;
   cut_data *cut;
   
   for (i = length - 1; i >= 0; i--) {
      j = index[i];
      cut = rows[j].cut;
      switch (cut->sense){
       case 'E': lb[j] = ub[j] = cut->rhs; break;
       case 'L': lb[j] = - OSL_INFINITY; ub[j] = cut->rhs; break;
       case 'G': lb[j] = cut->rhs; ub[j] = OSL_INFINITY; break;
       case 'R':
	 if (lp_data->mip->rngval[j] >= 0) {
	    ub[j] = cut->rhs; lb[j] = ub[j] - lp_data->mip->rngval[j];
	 } else {
	    ub[j] = cut->rhs; lb[j] = ub[j] + lp_data->mip->rngval[j];
	 }
	 break;
       default: /*This should never happen ... */
	 osllib_status = -1;
	 OSL_check_error("load_lp - unknown type of constraint");
      }
   }
   
   j = 0;
   if (j){
      ekk_free(lb);
   }
   
   osllib_status = ekk_setRowlower(lp_data->lp, lb);
   OSL_check_error("constrain_row_set ekk_setRowLower");
   ekk_free(lb);
   osllib_status = ekk_setRowupper(lp_data->lp, ub);
   OSL_check_error("constrain_row_set ekk_setRowUpper");
   ekk_free(ub);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

int read_mps(MIPdesc *desc, char *infile, char *probname)
{
   printf("\nMps-format file can be read only through OSI interface.\n");

   return(1);
}

/*===========================================================================*/

void write_mps(LPdata *lp_data, char *fname)
{
   osllib_status = ekk_exportModel(lp_data->lp, fname, 1, 2);
   OSL_check_error("write_mps");
}

/*===========================================================================*/

void write_sav(LPdata *lp_data, char *fname)
{
   osllib_status = ekk_saveModel(lp_data->lp, fname);
   OSL_check_error("write_sav");
}

/*===========================================================================*/

#ifdef USE_CGL_CUTS
void generate_cgl_cuts(LPdata *lp_data, int *num_cuts, cut_data ***cuts,
		       char send_to_pool)
{
   return;
}
#endif

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
#if 0
   lp_data->lpetol = 1e-09;
   cpx_status = CPXsetdblparam(lp_data->cpxenv, CPX_PARAM_EPRHS,
			       lp_data->lpetol);
   CPX_check_error("open_lp_solver - CPXsetdblparam");
#else
   cpx_status = CPXgetdblparam(lp_data->cpxenv, CPX_PARAM_EPRHS,
			       &lp_data->lpetol);
   CPX_check_error("open_lp_solver - CPXgetdblparam");
#endif
}

/*===========================================================================*/

void close_lp_solver(LPdata *lp_data)
{
   if (lp_data->lp){
      cpx_status = CPXfreeprob(lp_data->cpxenv, &(lp_data->lp));
      CPX_check_error("close_lp_solver");
      lp_data->lp = NULL;
   }
   cpx_status = CPXcloseCPLEX(&lp_data->cpxenv);
   CPX_check_error("close_lp_solver");
}

/*===========================================================================*/

/*===========================================================================*\
 * This function loads the data of an lp into the lp solver. This involves
 * transforming the data into CPLEX format and calling the CPLEX function
 * 'loadlp'.
\*===========================================================================*/

void load_lp_prob(LPdata *lp_data, int scaling, int fastmip)
{
   int i, *matcnt, *matbeg;

   /* realloc_lp_arrays(lp_data); */

   matcnt = (int *) malloc (lp_data->n*ISIZE);
   matbeg = lp_data->mip->matbeg;
   for (i = lp_data->n - 1; i >= 0; i--)
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
   printf("\nSorry, CPLEX versions 6.0 and earlier are no longer supported");
   printf("due to incompatibilities in memory management.");
   printf("Please use SYMPHONY 3.0.1 or earlier.\n\n");
   FREE(matcnt);
   exit(-1);
   /* legacy code left for posterity */
#else /* This is for the new memory model (CPLEX manages memory) */
   lp_data->lp = CPXcreateprob(lp_data->cpxenv,&cpx_status,(char *) "BB_prob");
   CPX_check_error("load_lp - CPXcreateprob");
   cpx_status = CPXcopylp(lp_data->cpxenv, lp_data->lp,
		lp_data->n, lp_data->m, 1, lp_data->mip->obj,
		lp_data->mip->rhs, lp_data->mip->sense,lp_data->mip->matbeg,
                matcnt, lp_data->mip->matind, lp_data->mip->matval,
		lp_data->mip->lb, lp_data->mip->ub, lp_data->mip->rngval);
   CPX_check_error("load_lp - CPXcopylp");
   FREE(matcnt);
#endif
}

/*===========================================================================*/

void unload_lp_prob(LPdata *lp_data)
{
   cpx_status = CPXfreeprob(lp_data->cpxenv, &lp_data->lp);
   CPX_check_error("unload_lp - CPXfreeprob");
   lp_data->lp = NULL;

   lp_data->m = lp_data->n = lp_data->nz = 0;
}

/*===========================================================================*/

void load_basis(LPdata *lp_data, int *cstat, int *rstat)
{
   cpx_status = CPXcopybase(lp_data->cpxenv, lp_data->lp, cstat, rstat);
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
   cpx_status = CPXprimopt(lp_data->cpxenv, lp_data->lp);
   CPX_check_error("refactorize - CPXoptimize");
   cpx_status = CPXsetintparam(lp_data->cpxenv, CPX_PARAM_ITLIM, itlim);
   CPX_check_error("refactorize - CPXsetintparam");
}

/*===========================================================================*/

void add_rows(LPdata *lp_data, int rcnt, int nzcnt, double *rhs,
	      char *sense, int *rmatbeg, int *rmatind, double *rmatval)
{
   int i, m = lp_data->m, j, indicator = FALSE;

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
   if (term == CPX_STAT_UNBOUNDED){
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
    case CPX_STAT_OPTIMAL:                        term = LP_OPTIMAL; break;
    case CPX_STAT_INFEASIBLE:                     term = LP_D_UNBOUNDED; break;
    case CPX_STAT_UNBOUNDED:                      term = LP_D_INFEASIBLE; break;
    case CPX_STAT_ABORT_OBJ_LIM:                  term = LP_D_OBJLIM; break;
    case CPX_STAT_ABORT_IT_LIM:                   term = LP_D_ITLIM; break;
    default:                                      term = LP_ABANDONED; break;
   }
#else
   switch (real_term = CPXgetstat(lp_data->cpxenv,lp_data->lp)){
    case CPX_OPTIMAL:                             term = LP_OPTIMAL; break;
    case CPX_INFEASIBLE:                          term = LP_D_INFEASIBLE; break;
    case CPX_UNBOUNDED:                           term = LP_D_UNBOUNDED; break;
    case CPX_OBJ_LIM:                             term = LP_D_OBJLIM; break;
    case CPX_IT_LIM_FEAS: case CPX_IT_LIM_INFEAS: term = LP_D_ITLIM; break;
    default:                                      term = LP_ABANDONED; break;
   }
#endif

   lp_data->termcode = term;

   if (term != LP_ABANDONED){
      *iterd = CPXgetitcnt(lp_data->cpxenv, lp_data->lp);
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
	     double *rowval, int *rowind, int *rowlen,
	     double *rowub, double *rowlb)
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
      if (lp_data->x){
	 return(bd < lp_data->x[jdiv] ? LOWER_THAN_LB : HIGHER_THAN_UB);
      }
   }
   return(0); /* fake return */
}

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
   cpx_status = CPXgetpi(lp_data->cpxenv, lp_data->lp, lp_data->dualsol, 0,
			 lp_data->m-1);
   CPX_check_error("get_dj_pi - CPXgetpi");
   cpx_status = CPXgetdj(lp_data->cpxenv, lp_data->lp, lp_data->dj, 0,
			 lp_data->n-1);
   CPX_check_error("get_dj_pi - CPXgetdj");
}

/*===========================================================================*/

void get_slacks(LPdata *lp_data)
{
   row_data *rows = lp_data->rows;
   double *slacks = lp_data->slacks;
   int i, m = lp_data->m;
   cpx_status = CPXgetslack(lp_data->cpxenv, lp_data->lp, lp_data->slacks, 0,
			    lp_data->m-1);
   CPX_check_error("get_slacks");
   /* Compute the real slacks for the free rows */
   for (i=m-1; i>=0; i--){
      if (rows[i].free){
	 switch (rows[i].cut->sense){
	  case 'E': slacks[i] +=  rows[i].cut->rhs - SYM_INFINITY; break;
	  case 'L': slacks[i] +=  rows[i].cut->rhs - SYM_INFINITY; break;
	  case 'G': slacks[i] +=  rows[i].cut->rhs + SYM_INFINITY; break;
	  case 'R': slacks[i] += -rows[i].cut->rhs - SYM_INFINITY; break;
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

void get_bounds(LPdata *lp_data)
{
   if (!lp_data->lb){
   }
   cpx_status = CPXgetlb(lp_data->cpxenv, lp_data->lp, lp_data->lb, 0,
			 lp_data->n-1);
   CPX_check_error("get_lb");
   cpx_status = CPXgetub(lp_data->cpxenv, lp_data->lp, lp_data->ub, 0,
			 lp_data->n-1);
   CPX_check_error("get_ub");
   
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
   cpx_status = CPXdelsetrows(lp_data->cpxenv, lp_data->lp, free_rows);
   CPX_check_error("delete_rows");
   lp_data->m -= deletable;
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

int delete_cols(LPdata *lp_data, int delnum, int *delstat)
{
   double *dj = lp_data->dj;
   double *x = lp_data->x;
   char *status = lp_data->status;
   int i, num_to_keep;

   cpx_status = CPXdelsetcols(lp_data->cpxenv, lp_data->lp, delstat);
   CPX_check_error("delete_cols - CPXdelsetcols");
   lp_data->nz = CPXgetnumnz(lp_data->cpxenv, lp_data->lp);
   CPX_check_error("delete_cols - CPXgetnumnz");

   for (i = 0, num_to_keep = 0; i < lp_data->n; i++){
      if (delstat[i] != -1){
	 dj[delstat[i]] = dj[i];
	 x[delstat[i]] = x[i];
	 status[delstat[i]] = status[i];
	 num_to_keep++;
      }
   }

   lp_data->n = num_to_keep;

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
   row_data *rows = lp_data->rows;
   double *rhsval = lp_data->tmp.d; /* m */
   int *ind_e = lp_data->tmp.i1 + 2 * lp_data->m; /* m (now) */
   /* See comment in check_row_effectiveness why the shift! */
   char *sen_e = lp_data->tmp.c; /* m (now) */

   for (j=0, i=length-1; i>=0; i--){
      switch (rows[index[i]].cut->sense){
       case 'E': rhsval[i] = SYM_INFINITY; ind_e[j++] = index[i]; break;
       case 'L': rhsval[i] = SYM_INFINITY; break;
       case 'R':
       cpx_status = CPXchgcoef(lp_data->cpxenv, lp_data->lp, index[i], -2,
                               2*SYM_INFINITY);
       CPX_check_error("free_row_set - CPXchgcoef");
       case 'G': rhsval[i] = -SYM_INFINITY; break;
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
   row_data *rows = lp_data->rows;
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

int read_mps(MIPdesc *desc, char *infile, char *probname)
{
   printf("\nMps-format file can be read only through OSI interface.\n");

   return(1);
}

/*===========================================================================*/

void write_mps(LPdata *lp_data, char *fname)
{
   cpx_status = CPXmpswrite(lp_data->cpxenv, lp_data->lp, fname);
   CPX_check_error("write_mps");
}

/*===========================================================================*/

void write_sav(LPdata *lp_data, char *fname)
{
   cpx_status = CPXsavwrite(lp_data->cpxenv, lp_data->lp, fname);
   CPX_check_error("write_sav");
}

/*===========================================================================*/

#ifdef USE_CGL_CUTS
void generate_cgl_cuts(LPdata *lp_data, int *num_cuts, cut_data ***cuts,
		       char send_to_pool)
{
   return;
}
#endif

#endif /* __CPLEX__ */

#if defined(__OSI_CPLEX__) || defined(__OSI_OSL__) || defined(__OSI_CLP__) \
|| defined(__OSI_XPRESS__) || defined(__OSI_SOPLEX__) || defined(__OSI_VOL__) \
|| defined(__OSI_DYLP__) || defined (__OSI_GLPK__)

static bool retval = false;

void open_lp_solver(LPdata *lp_data)
{
   lp_data->si = new OsiXSolverInterface();

   /* Turn off the OSL messages (There are LOTS of them) */
   lp_data->si->setHintParam(OsiDoReducePrint);
   lp_data->si->messageHandler()->setLogLevel(0);
#ifdef __OSI_GLPK__
   lp_data->lpetol = 1e-07; /* glpk doesn't return the value of this param */ 
#else   
   lp_data->si->getDblParam(OsiPrimalTolerance, lp_data->lpetol);
#endif
}

/*===========================================================================*/

void close_lp_solver(LPdata *lp_data)
{
   delete lp_data->si;
}

/*===========================================================================*/

/*===========================================================================*\
 * This function loads the data of an lp into the lp solver. This involves
 * transforming the data into CPLEX format and calling the CPLEX function
 * 'loadlp'.
\*===========================================================================*/

void load_lp_prob(LPdata *lp_data, int scaling, int fastmip)
{

   /* Turn off scaling for CLP */
   //lp_data->si->setHintParam(OsiDoScale,false,OsiHintDo);

   lp_data->si->loadProblem(lp_data->n, lp_data->m,
			    lp_data->mip->matbeg, lp_data->mip->matind,
			    lp_data->mip->matval, lp_data->mip->lb,
			    lp_data->mip->ub, lp_data->mip->obj,
			    lp_data->mip->sense, lp_data->mip->rhs,
			    lp_data->mip->rngval);
}

/*===========================================================================*/

void unload_lp_prob(LPdata *lp_data)
{

   lp_data->si->reset();
  
   /* Set parameters as in open_lp_solver() (do these persist?) */
   lp_data->si->setHintParam(OsiDoReducePrint);
   lp_data->si->messageHandler()->setLogLevel(0);

   lp_data->m = lp_data->n = lp_data->nz = 0;
}

/*===========================================================================*/

void load_basis(LPdata *lp_data, int *cstat, int *rstat)
{
   
   CoinWarmStartBasis *warmstart = new CoinWarmStartBasis;

   int numcols = lp_data->n;
   int numrows = lp_data->m;
   int i;
   
   warmstart->setSize(numcols, numrows);
   
   for (i = 0; i < numrows; i++){
      switch (rstat[i]){
       case SLACK_AT_LB:
	 warmstart->setArtifStatus(i,CoinWarmStartBasis::atLowerBound);
	 break;
       case SLACK_BASIC:
	 warmstart->setArtifStatus(i,CoinWarmStartBasis::basic);
	 break;
       case SLACK_AT_UB:
	 warmstart->setArtifStatus(i,CoinWarmStartBasis::atUpperBound);
	 break;
       case SLACK_FREE:
	 warmstart->setArtifStatus(i,CoinWarmStartBasis::isFree);
	 break;
       default:
	 break;
      }
   }
   
   for (i = 0; i < numcols; i++){
      switch (cstat[i]){
       case VAR_AT_LB:
	 warmstart->setStructStatus(i,CoinWarmStartBasis::atLowerBound);
	 break;
       case VAR_BASIC:
	 warmstart->setStructStatus(i,CoinWarmStartBasis::basic);
	 break;
       case VAR_AT_UB:
	 warmstart->setStructStatus(i,CoinWarmStartBasis::atUpperBound);
	 break;
       case VAR_FREE:
	 warmstart->setStructStatus(i,CoinWarmStartBasis::isFree);
	 break;
       default:
	 break;
      }
   }

   retval = lp_data->si->setWarmStart(warmstart);

   delete warmstart;
}

/*===========================================================================*/

void add_rows(LPdata *lp_data, int rcnt, int nzcnt, double *rhs,
	      char *sense, int *rmatbeg, int *rmatind, double *rmatval)
{
   int i, j, m = lp_data->m;
   
   for (i = 0; i < rcnt; i++){
      CoinPackedVector new_row;
      for (j = rmatbeg[i]; j < rmatbeg[i+1]; j++){
	 new_row.insert(rmatind[j], rmatval[j]);
      }
      lp_data->si->addRow(new_row, sense[i], rhs[i], 0);
   }

   //  lp_data->si->addRows(rcnt,rows,sense,rhs,rowrng);

   lp_data->m += rcnt;
   lp_data->nz += nzcnt;
   lp_data->lp_is_modified=LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void add_cols(LPdata *lp_data, int ccnt, int nzcnt, double *obj,
	      int *cmatbeg, int *cmatind, double *cmatval,
	      double *lb, double *ub, char *where_to_move)
{
   int i, j;
   for (i = 0; i < ccnt; i++){
      CoinPackedVector col;
      for (j = cmatbeg[i]; j < cmatbeg[i+1]; j++)
	 col.insert(cmatind[j], cmatval[j]);
      lp_data->si->addCol(col, lb[i], ub[i], obj[i]);
   }
   
   lp_data->n += ccnt;
   lp_data->nz += nzcnt;
}

/*===========================================================================*/

void change_row(LPdata *lp_data, int row_ind,
		char sense, double rhs, double range)
{
   lp_data->si->setRowType(row_ind, sense, rhs, range);
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
}

/*===========================================================================*/

/*===========================================================================*\
 * Solve the lp specified in lp_data->lp with dual simplex. The number of
 * iterations is returned in 'iterd'. The return value of the function is
 * the termination code of the dual simplex method.
\*===========================================================================*/

int dual_simplex(LPdata *lp_data, int *iterd)
{
   
   //int term = LP_ABANDONED;
   int term;
    
   lp_data->si->resolve();
   
   if (lp_data->si->isProvenDualInfeasible())
      term = LP_D_INFEASIBLE;
   else if (lp_data->si->isProvenPrimalInfeasible())
      term = LP_D_UNBOUNDED;
   else if (lp_data->si->isProvenOptimal())
      term = LP_OPTIMAL;
   else if (lp_data->si->isDualObjectiveLimitReached())
      term = LP_D_OBJLIM;
   else if (lp_data->si->isIterationLimitReached())
      term = LP_D_ITLIM;
   else if (lp_data->si->isAbandoned())
      term = LP_ABANDONED;
   
   /* if(term == D_UNBOUNDED){
      retval=lp_data->si->getIntParam(OsiMaxNumIteration, itlim); 
      CAN NOT GET DEFAULT, MIN VALUES in OSI of CPXinfointparam() */
   
   lp_data->termcode = term;
   
   if (term != LP_ABANDONED){
      
      *iterd = lp_data->si->getIterationCount();
      
      lp_data->objval = lp_data->si->getObjValue();
      
      lp_data->lp_is_modified = LP_HAS_NOT_BEEN_MODIFIED;
   }   
   else{
      lp_data->lp_is_modified = LP_HAS_BEEN_ABANDONED;
      printf("OSI Abandoned calculation: Code %i \n\n", term);
   }
   
   return(term);
}


/*===========================================================================*/

/*===========================================================================*/
/* This function is used only together with get_proof_of_infeasibility...    */

void get_binvrow(LPdata *lp_data, int i, double *row)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/

void get_basis(LPdata *lp_data, int *cstat, int *rstat)
{
   CoinWarmStart * warmstart = lp_data->si->getWarmStart();
   
   CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>(warmstart);
   
   int numcols = ws->getNumStructural();   /* has to be <= lp_data->n */
   int numrows = ws->getNumArtificial();   /* has to be <= lp_data->m */
   int i;                                  /* hence an assert? */  
   
   if (rstat){
      for (i = 0; i < numrows; i++){
	 switch (ws->getArtifStatus(i)){
	  case CoinWarmStartBasis::basic:
	    rstat[i] = SLACK_BASIC;
	    break;
	  case CoinWarmStartBasis::atLowerBound:
	    rstat[i] = SLACK_AT_LB;
	    break;
	  case CoinWarmStartBasis::atUpperBound:
	    rstat[i] = SLACK_AT_UB;
	    break;
	  case CoinWarmStartBasis::isFree:     //can it happen?
	    rstat[i] = SLACK_FREE;
	    break;
	  default:
	    break;                            //can it happen?
	 }
      }
   }
   
   if (cstat){
      for (i = 0; i < numcols; i++){
	 switch (ws->getStructStatus(i)){
	  case CoinWarmStartBasis::basic:
	    cstat[i] = VAR_BASIC;
	    break;
	  case CoinWarmStartBasis::atLowerBound:
	    cstat[i] = VAR_AT_LB;
	    break;
	  case CoinWarmStartBasis::atUpperBound:
	    cstat[i] = VAR_AT_UB;
	    break;
	  case CoinWarmStartBasis::isFree:  
	    cstat[i] = VAR_FREE;
	    break;
	  default:
	    break;                            //can it happen?
	 }
      }
   }

   delete ws;
}

/*===========================================================================*/

/*===========================================================================*\
 * Set an upper limit on the objective function value. Call the 'setobjulim'
 * CPLEX function.
\*===========================================================================*/

void set_obj_upper_lim(LPdata *lp_data, double lim)
{
   OsiDblParam key = OsiDualObjectiveLimit;
   
   retval = lp_data->si->setDblParam(key, lim);
}

/*===========================================================================*/

/*===========================================================================*\
 * Set an upper limit on the number of iterations. If itlim < 0 then set
 * it to the maximum.
\*===========================================================================*/

void set_itlim(LPdata *lp_data, int itlim)
{
   if (itlim < 0) itlim = LP_MAX_ITER;

   OsiIntParam key = OsiMaxNumIteration;
   
   retval = lp_data->si->setIntParam(key, itlim);
}

/*===========================================================================*/

void get_column(LPdata *lp_data, int j,
		double *colval, int *colind, int *collen, double *cj)
{
   const CoinPackedMatrix *matrixByCol = lp_data->si->getMatrixByCol();
   
   int nc = matrixByCol->getNumCols(), i;
   
   const double *matval = matrixByCol->getElements();
   const int *matind = matrixByCol->getIndices(); 
   const int *matbeg = matrixByCol->getVectorStarts();

   *collen = matrixByCol->getVectorSize(j);
   
   for (i = 0; i < (*collen); i++){
      colval[i] = matval[matbeg[j] + i];
      colind[i] = matind[matbeg[j] + i];
   }
   
   const double * objval = lp_data->si->getObjCoefficients();
   
   *cj = objval[j];
}

/*===========================================================================*/

void get_row(LPdata *lp_data, int i,
	     double *rowval, int *rowind, int *rowlen, 
	     double *rowub, double *rowlb)
{
   const CoinPackedMatrix * matrixByRow = lp_data->si->getMatrixByRow();
  
   int nr = matrixByRow->getNumRows(), j = 0;
   
   const double *matval = matrixByRow->getElements();  
   const int *matind = matrixByRow->getIndices(); 
   const int *matbeg = matrixByRow->getVectorStarts();

   *rowlen = matrixByRow->getVectorSize(i);
   *rowub = lp_data->si->getRowUpper()[i];
   *rowlb = lp_data->si->getRowLower()[i];      

   for (j = 0; j < (*rowlen); j++){
      rowval[j] = matval[matbeg[i] + j];
      rowind[j] = matind[matbeg[i] + j];
   }
}

/*===========================================================================*/

/* This routine returns the index of a row which proves the lp to be primal
 * infeasible. There must be one, or this function wouldn't be called. */

int get_proof_of_infeas(LPdata *lp_data, int *infind)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/

/*===========================================================================*\
 * Get the solution (values of the structural variables in an optimal
 * solution) to the lp (specified by lp_data->lp) into the vector
 * lp_data->x. This can be done by calling the 'getx' CPLEX function.
\*===========================================================================*/

void get_x(LPdata *lp_data)
{
   memcpy(lp_data->x, const_cast<double *>(lp_data->si->getColSolution()),
	  lp_data->n * DSIZE);
}

/*===========================================================================*/

void get_dj_pi(LPdata *lp_data)
{
   memcpy(lp_data->dualsol, const_cast<double *>(lp_data->si->getRowPrice()),
	  lp_data->m * DSIZE);
   memcpy(lp_data->dj, const_cast<double *>(lp_data->si->getReducedCost()),
	  lp_data->n * DSIZE);
}

/*===========================================================================*/

void get_slacks(LPdata *lp_data)
{
   int m = lp_data->m, i = 0;
   double * slacks = lp_data->slacks;
   row_data *rows = lp_data->rows;
   
#ifndef __OSI_CPLEX__
   
   const double * rowActivity = lp_data->si->getRowActivity();
   
   for (i = m - 1; i >= 0; i--) {
      if ((rows[i].cut->sense == 'R') && (rows[i].cut->range < 0) ) {
	 slacks[i] = - rows[i].cut->rhs + rowActivity[i];
      } else {
	 slacks[i] = rows[i].cut->rhs - rowActivity[i];
      }
   }

#else
   
   CPXgetslack(lp_data->si->getEnvironmentPtr(), lp_data->si->getLpPtr(),
	       lp_data->slacks, 0, lp_data->m-1);
   /* Compute the real slacks for the free rows */
   for (i =m - 1; i >= 0; i--){
      if (rows[i].free){
	 switch (rows[i].cut->sense){
	  case 'E': slacks[i] +=  rows[i].cut->rhs - SYM_INFINITY; break;
	  case 'L': slacks[i] +=  rows[i].cut->rhs - SYM_INFINITY; break;
	  case 'G': slacks[i] +=  rows[i].cut->rhs + SYM_INFINITY; break;
	  case 'R': slacks[i] += -rows[i].cut->rhs - SYM_INFINITY; break;
	 }
      }
   }

#endif
}

/*===========================================================================*/

void change_range(LPdata *lp_data, int rowind, double value)
{

   double rhs = lp_data->si->getRightHandSide()[rowind]; 

   lp_data->si->setRowType(rowind,'R', rhs, value);
}

/*===========================================================================*/

void change_rhs(LPdata *lp_data, int rownum, int *rhsind, double *rhsval)
{
   char *sense = lp_data->tmp.c; 
   double *range = lp_data->tmp.d; 
   int i;

   for (i = 0; i < rownum; i++){
      sense[i] = lp_data->si->getRowSense()[rhsind[i]];
      if (sense[i] == 'R'){
	 range[i] = lp_data->si->getRowRange()[rhsind[i]];
      }
   }
   
   lp_data->si->setRowSetTypes(rhsind, rhsind + rownum, sense, rhsval, range);
}

/*===========================================================================*/

void change_sense(LPdata *lp_data, int cnt, int *index, char *sense)
{
  double *rhs = lp_data->tmp.d; 
  double *range = (double *) calloc(cnt, DSIZE);
  char range_used = FALSE;
  int i; 

  for (i = 0; i < cnt; i++){
     rhs[i] = lp_data->si->getRightHandSide()[index[i]];
     if (sense[i] == 'R')
	range[i] = lp_data->si->getRowRange()[index[i]];
  }

  if (!range_used){
     FREE(range);
  }
  
  lp_data->si->setRowSetTypes(index, index + cnt, sense, rhs, range);
}

/*===========================================================================*/

void change_bounds(LPdata *lp_data, int cnt, int *index, char *lu, double *bd)
{
   int i;
 
   for (i = 0; i < cnt; i++){
      switch (lu[i]){
       case 'L':
	 lp_data->si->setColLower(index[i], bd[i]); 
	 break;
       case 'U':
	 lp_data->si->setColUpper(index[i], bd[i]); 
	 break;
       default:
	 /* default: can't happen */
	 break;
      }
   }
   
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_lbub(LPdata *lp_data, int j, double lb, double ub)
{
   lp_data->si->setColBounds(j,lb,ub);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_ub(LPdata *lp_data, int j, double ub)
{
   lp_data->si->setColUpper(j,ub);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void change_lb(LPdata *lp_data, int j, double lb)
{
   lp_data->si->setColLower(j,lb);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
}

/*===========================================================================*/

void get_ub(LPdata *lp_data, int j, double *ub)
{
   *ub=lp_data->si->getColUpper()[j];
}

/*===========================================================================*/

void get_lb(LPdata *lp_data, int j, double *lb)
{
   *lb=lp_data->si->getColLower()[j];
}

/*===========================================================================*/

void get_bounds(LPdata *lp_data)
{
   lp_data->ub = const_cast<double *>(lp_data->si->getColUpper());
   lp_data->lb = const_cast<double *>(lp_data->si->getColLower());
}

/*===========================================================================*/

void get_objcoef(LPdata *lp_data, int j, double *objcoef)
{
   *objcoef=lp_data->si->getObjCoefficients()[j];
}

/*===========================================================================*/

void delete_rows(LPdata *lp_data, int deletable, int *free_rows)
{
   
   int i, m = lp_data->m;
   int *which = lp_data->tmp.i1 + lp_data->m;
   int delnum = 0;
   
   CoinFillN(which, deletable, 0);
   
   for (i = 0; i < m; i++){
      if (free_rows[i]){
	 which[delnum++] = i;
      }
   }
   
   lp_data->si->deleteRows(delnum, which);
   lp_data->nz = lp_data->si->getNumElements();
   lp_data->m -= delnum;
}

/*===========================================================================*/

int delete_cols(LPdata *lp_data, int delnum, int *delstat)
{
   int i, m = lp_data->m, n = lp_data->n;
   int *which = (int *) calloc(delnum, ISIZE);
   int num_to_delete = 0, num_to_keep = 0;
   double *dj = lp_data->dj;
   double *x = lp_data->x;
   char *status = lp_data->status;

   for (i = n - 1; i >= 0; i--){
      if (delstat[i]){
	 which[num_to_delete++]=i;
      }
   }

   lp_data->si->deleteCols(num_to_delete, which);
   lp_data->nz = lp_data->si->getNumElements();
   FREE(which);

   /* make result as CPLEX does */
   for (i = 0, num_to_keep = 0; i < lp_data->n; i++){
      if (delstat[i]){
	 delstat[i] = -1;
      }else{
	 delstat[i] = num_to_keep++;
	 dj[delstat[i]] = dj[i];
	 x[delstat[i]] = x[i];
	 status[delstat[i]] = status[i];
      }
   }

   lp_data->n = num_to_keep;

   return(num_to_delete);
}

/*===========================================================================*/

void release_var(LPdata *lp_data, int j, int where_to_move)
{
  fprintf(stderr, "Function not implemented yet.");
  exit(-1);
}

/*===========================================================================*/

void free_row_set(LPdata *lp_data, int length, int *index)
{
   char *sense = lp_data->tmp.c; /* m (now) */
   double *rhs = lp_data->tmp.d; /* m */
   double *range = (double *) calloc(length, DSIZE);
   char range_used = FALSE;
   int i; 
   
   for (i = 0; i < length; i++){
      rhs[i] = lp_data->si->getRightHandSide()[index[i]];
      sense[i] = lp_data->si->getRowSense()[index[i]];
      if (sense[i] =='R'){
	 range[i] = lp_data->si->getRowRange()[index[i]];
      }
   }
   
   for (i = 0; i < length; i++) {
     //     range[i]=0;
     switch (sense[i]){
     case 'E':
       rhs[i] = lp_data->si->getInfinity();
       sense[i] = 'L';
       break;
     case 'L':
       rhs[i] = lp_data->si->getInfinity(); 
       break;
     case 'R':
       range[i] = 2*lp_data->si->getInfinity();
       range_used = TRUE;
       break;
     case 'G':
       rhs[i] = -lp_data->si->getInfinity();
      }
   }

   lp_data->si->setRowSetTypes(index, index + length, sense, rhs, range);

     FREE(range);
}

/*===========================================================================*/

void constrain_row_set(LPdata *lp_data, int length, int *index)
{
   char *sense = lp_data->tmp.c; /* m (now) */
   double *rhs = lp_data->tmp.d; /* m */
   double *range = (double *) calloc(length, DSIZE);
   row_data *rows = lp_data->rows;
   cut_data *cut;
   char range_used = FALSE;
   int i;
   
   for (i = length - 1; i >= 0; i--){
      cut = rows[index[i]].cut;
      rhs[i] = cut->rhs;  
      if ((sense[i] = cut->sense) == 'R'){
	 range[i] = cut->range;
	 range_used = TRUE;
      }
   }

   lp_data->si->setRowSetTypes(index, index + length, sense, rhs, range);

   FREE(range);

}

/*===========================================================================*/

int read_mps(MIPdesc *mip, char *infile, char *probname)
{
   int i, j, k, last_dot = 0, last_dir = 0, nonzeros = 0;
   char fname[80] = "", ext[10] = "", slash;

   CoinMpsIO mps;
   int errors;
    
#ifdef WIN32
	slash = '\\';
#else
	slash = '/';
#endif
	
   for (j = 0;; j++){
      if (infile[j] == '\0')
	 break;
      if (infile[j] == '.') {
	    last_dot = j;
	  }
	  if(infile[j] == slash){
		last_dir = j;
	  }
   }
   
   if(last_dir < last_dot){
	   memcpy(fname, infile, CSIZE*last_dot);
	   memcpy(ext, infile + last_dot + 1, CSIZE*(j - last_dot - 1)); 
   }
   else{
	   memcpy(fname, infile, CSIZE*j);
   }


   mps.setInfinity(mps.getInfinity());
   
   if (errors = mps.readMps(fname,ext)){
      return(errors);
   }
   
   strncpy(probname, const_cast<char *>(mps.getProblemName()), 80);
   
   mip->m  = mps.getNumRows();
   mip->n  = mps.getNumCols();
   mip->nz = mps.getNumElements();
   
   mip->obj    = (double *) malloc(DSIZE * mip->n);
   mip->obj1   = (double *) calloc(mip->n, DSIZE);
   mip->obj2   = (double *) calloc(mip->n, DSIZE);
   mip->rhs    = (double *) malloc(DSIZE * mip->m);
   mip->sense  = (char *)   malloc(CSIZE * mip->m);
   mip->rngval = (double *) malloc(DSIZE * mip->m);
   mip->ub     = (double *) malloc(DSIZE * mip->n);
   mip->lb     = (double *) malloc(DSIZE * mip->n);
   mip->is_int = (char *)   calloc(CSIZE, mip->n);
   
   memcpy(mip->obj, const_cast <double *> (mps.getObjCoefficients()),
	  DSIZE * mip->n); 
   memcpy(mip->rhs, const_cast <double *> (mps.getRightHandSide()),
	  DSIZE * mip->m); 
   memcpy(mip->sense, const_cast <char *> (mps.getRowSense()),
	  CSIZE * mip->m); 
   memcpy(mip->rngval, const_cast <double *> (mps.getRowRange()),
	  DSIZE * mip->m); 
   memcpy(mip->ub, const_cast <double *> (mps.getColUpper()),
	  DSIZE * mip->n); 
   memcpy(mip->lb, const_cast <double *> (mps.getColLower()),
	  DSIZE * mip->n); 
   
   //user defined matind, matval, matbeg--fill as column ordered
   
   const CoinPackedMatrix * matrixByCol= mps.getMatrixByCol();
   
   mip->matbeg = (int *) malloc(ISIZE * (mip->n + 1));
   memcpy(mip->matbeg, const_cast<int *>(matrixByCol->getVectorStarts()),
	  ISIZE * (mip->n + 1));
   
   mip->matval = (double *) malloc(DSIZE*mip->matbeg[mip->n]);
   mip->matind = (int *)    malloc(ISIZE*mip->matbeg[mip->n]);
   
   memcpy(mip->matval, const_cast<double *> (matrixByCol->getElements()),
	  DSIZE * mip->matbeg[mip->n]);  
   memcpy(mip->matind, const_cast<int *> (matrixByCol->getIndices()), 
	  ISIZE * mip->matbeg[mip->n]);  
   
   for (j = 0; j < mip->n; j++){
      mip->is_int[j] = mps.isInteger(j);
   }

   mip->obj_offset = -mps.objectiveOffset();

   mip->colname = (char **) malloc(sizeof(char *) * mip->n);   
   
   for (j = 0; j < mip->n; j++){
      mip->colname[j] = (char *) malloc(CSIZE * 9);
      strncpy(mip->colname[j], const_cast<char*>(mps.columnName(j)), 9);
      mip->colname[j][8] = 0;
   }
   
   return(errors);
}

/*===========================================================================*/

void write_mps(LPdata *lp_data, char *fname)
{
   char * extension = "MPS";
   double ObjSense = lp_data->si->getObjSense();
   int i;
   
   for (i = 0; i < lp_data->n; i++) {
      lp_data->si->setContinuous(i);
   }

   lp_data->si->writeMps(fname, extension, ObjSense);
}

/*===========================================================================*/

void write_sav(LPdata *lp_data, char *fname)
{
   fprintf(stderr, "Function not implemented yet.");
   exit(-1);
}

/*===========================================================================*/

#ifdef USE_CGL_CUTS

#include "qsortucb.h"

void generate_cgl_cuts(LPdata *lp_data, int *num_cuts, cut_data ***cuts,
		       char send_to_pool, int is_rootnode)
{
   OsiCuts cutlist;
   OsiRowCut cut;
   int i = 0, j = 0, k = 0; 
   int *matind;
   double *matval;
   cgl_params *par = &(lp_data->cgl);
   int cut_num = 0;

#ifndef COMPILE_IN_LP
   par->probing_generated_in_root =
   par->gomory_generated_in_root =
   par->oddhole_generated_in_root =
   par->mir_generated_in_root = 
   par->clique_generated_in_root = 
   par->flow_and_cover_generated_in_root =
   par->rounding_generated_in_root =
   par->lift_and_project_generated_in_root = TRUE;
#endif
      

   /* Set proper variables to be integer */
   for (i = 0; i < lp_data->n; i++) {
      if (lp_data->vars[i]->is_int) { // integer or binary
	 lp_data->si->setInteger(i);
      }
   }  

   /* create CGL probing cuts */
   if(par->generate_cgl_probing_cuts){
     if(is_rootnode || par->probing_generated_in_root){
       CglProbing *probe = new CglProbing;
       //#if 0
       probe->setUsingObjective(true);
       probe->setMaxPass(3);
       probe->setMaxLook(50);
       probe->setRowCuts(3);
       //#endif
       probe->generateCuts(*(lp_data->si), cutlist);
       if (!par->probing_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->probing_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();
       delete probe;
     }
   }

   /* create CGL gomory cuts */
   if(par->generate_cgl_gomory_cuts){
     if(is_rootnode || par->gomory_generated_in_root){
       CglGomory *gomory = new CglGomory;
       gomory->generateCuts(*(lp_data->si), cutlist);
       if (!par->gomory_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->gomory_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete gomory;
     }
   }

   /* create CGL knapsack cuts */
   if(par->generate_cgl_knapsack_cuts){
     if(is_rootnode || par->knapsack_generated_in_root){
       CglKnapsackCover *knapsack = new CglKnapsackCover;
       knapsack->generateCuts(*(lp_data->si), cutlist);
       if (!par->knapsack_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->knapsack_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete knapsack;
     }
   }

   /* create CGL odd hole cuts */
   if(par->generate_cgl_oddhole_cuts){
     if(is_rootnode || par->oddhole_generated_in_root){
       CglOddHole *oddhole = new CglOddHole;
       //#if 0
       oddhole->setMinimumViolation(0.005);
       oddhole->setMinimumViolationPer(0.00002);
       oddhole->setMaximumEntries(200);
       //#endif
       oddhole->generateCuts(*(lp_data->si), cutlist);
       if (!par->oddhole_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->oddhole_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete oddhole;
     }
   }

   /* create CGL mir cuts */
   if(par->generate_cgl_mir_cuts){
     if(is_rootnode || par->mir_generated_in_root){
       CglMixedIntegerRounding *mir = new CglMixedIntegerRounding;
       mir->generateCuts(*(lp_data->si), cutlist);
       if (!par->mir_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->mir_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete mir;
     }
   }

   /* create CGL clique cuts */
   if(par->generate_cgl_clique_cuts){
     if(is_rootnode || par->clique_generated_in_root){
       CglClique *clique = new CglClique;
       clique->generateCuts(*(lp_data->si), cutlist);
       if (!par->clique_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->clique_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete clique;
     }
   }
   
   /* create CGL flow cover cuts */
   if(par->generate_cgl_flow_and_cover_cuts){
     if(is_rootnode || par->flow_and_cover_generated_in_root){
       CglFlowCover *flow = new CglFlowCover;
       flow->generateCuts(*(lp_data->si), cutlist);
       if (!par->flow_and_cover_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->flow_and_cover_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete flow;
      //printf("%i\n", cutlist.sizeRowCuts());
     }
   }

   /* create CGL simple rounding cuts */
   if(par->generate_cgl_rounding_cuts){
     if(is_rootnode || par->rounding_generated_in_root){
       CglSimpleRounding * rounding = new CglSimpleRounding;
       rounding->generateCuts(*(lp_data->si), cutlist);
       if (!par->rounding_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->rounding_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete rounding;
       //printf("%i\n", cutlist.sizeRowCuts());
     }
   }
   
   /* create CGL liftandproject cuts (currently buggy) */     
   if(par->generate_cgl_lift_and_project_cuts){
     if(is_rootnode || par->lift_and_project_generated_in_root){
       CglLiftAndProject *liftandproject = new CglLiftAndProject;
       liftandproject->generateCuts(*(lp_data->si), cutlist);
       if (!par->lift_and_project_generated_in_root){
	 if ((cutlist.sizeRowCuts() - cut_num) > 0) {
	   par->lift_and_project_generated_in_root = TRUE;
	 }
       }
       cut_num = cutlist.sizeRowCuts();       
       delete liftandproject;
     }
   }

   if (cutlist.sizeRowCuts() > 0){
      if (*cuts){
	 *cuts = (cut_data **)realloc(*cuts, (*num_cuts+cutlist.sizeRowCuts())
				      * sizeof(cut_data *));
      }else{
	 *cuts = (cut_data **)malloc(cutlist.sizeRowCuts()*sizeof(cut_data *));
      }
      for (i = 0, j = *num_cuts; i < cutlist.sizeRowCuts(); i++){
	 int num_elements;
	 int *indices;
	 double *elements;
	 cut = cutlist.rowCut(i);
	 (*cuts)[j] =  (cut_data *) calloc(1, sizeof(cut_data));
	 num_elements = cut.row().getNumElements();
	 indices = const_cast<int *> (cut.row().getIndices());
	 elements = const_cast<double *> (cut.row().getElements());
	 (*cuts)[j]->type = EXPLICIT_ROW;
	 if (((*cuts)[j]->sense = cut.sense()) == 'R'){
	    FREE((*cuts)[j]);
	    continue; /* This must be a bug. */
	 }
	 (*cuts)[j]->rhs = cut.rhs();
	 (*cuts)[j]->range = cut.range();
	 (*cuts)[j]->size = ISIZE + num_elements * (ISIZE + DSIZE);
	 (*cuts)[j]->coef = (char *) malloc ((*cuts)[j]->size);
	 ((int *) ((*cuts)[j]->coef))[0] = num_elements;
	 matind = (int *) ((*cuts)[j]->coef + ISIZE);
	 for (k = 0; k < num_elements; k++){
	    matind[k] = lp_data->vars[indices[k]]->userind;
	 }
	 matval = (double *) ((*cuts)[j]->coef + (num_elements + 1) * ISIZE);
	 memcpy((char *)matval, (char *)elements, num_elements * DSIZE);
	 qsortucb_id(matind, matval, num_elements);
	 (*cuts)[j]->branch = DO_NOT_BRANCH_ON_THIS_ROW;
	 (*cuts)[j]->deletable = TRUE;
	 if (send_to_pool){
	    (*cuts)[j++]->name = CUT__SEND_TO_CP;
	 }else{
	    (*cuts)[j++]->name = CUT__DO_NOT_SEND_TO_CP;
	 }	    
      }
      *num_cuts = j;
   }
   
   return;
}
#endif

#endif /* __OSI_xxx__ */
