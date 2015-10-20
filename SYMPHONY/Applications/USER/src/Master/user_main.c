/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/*===========================================================================*/

#define CALL_FUNCTION(f) \
if ((termcode = f) < 0){                                                    \
   printf("Error detected: termcode = %i\n", termcode);                     \
   printf("Exiting...\n\n");                                                \
   exit(termcode);                                                          \
}

/*===========================================================================*\
   This file contains the main() for the master process.

   Note that, if you want to use the OSI SYMPHONY interface, you should set the
   USE_OSI_INTERFACE flag and define the COINROOT path in the SYMPHONY 
   Makefile. Otherwise, the C callable library functions will be used by 
   default. See below for the usage.
\*===========================================================================*/

#if defined(USE_OSI_INTERFACE)

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{
   OsiSymSolverInterface si;

   /* Parse the command line */
   si.parseCommandLine(argc, argv);
   
   /* Read in the problem */
   si.loadProblem();

   /* Find a priori problem bounds */
   si.findInitialBounds();

   /* Solve the problem */
   si.branchAndBound();
   
   return(0);
}

#else

#include "symphony.h"
#include "sym_master.h"
#include <stdlib.h>

#include "user.h"

#define USER_FUNC_ERROR -1

int user_test(sym_environment *env);

int main(int argc, char **argv)
{

   int termcode;
   char * infile;

   /* Create a SYMPHONY environment */
   sym_environment *env = sym_open_environment();

   /* Print version info */
   sym_version();
   
   if (!env){
      printf("Error initializing environement\n");
      exit(0);
   }

   /* Create the data structure for storing the problem instance.*/
   user_problem *prob = (user_problem *)calloc(1, sizeof(user_problem));
   prob->mip = (MIPdesc *)calloc(1, sizeof(MIPdesc));
   
   CALL_FUNCTION( sym_set_user_data(env, (void *)prob) );

   CALL_FUNCTION( sym_parse_command_line(env, argc, argv) );

   CALL_FUNCTION( user_read_data(prob, prob->par.infile) );
   
   CALL_FUNCTION( user_load_problem(env, prob) );

   CALL_FUNCTION( sym_solve(env) );
     
   // Added by Suresh for debugging
   printf("\n\n");
   double obj_val = 0.0;
   int n = env->best_sol.xlength;
   int i = 0;
   for (i = 0; i < n; i++) {
      if (env->best_sol.xind[i] < prob->origvar_num) {
         obj_val += env->best_sol.xval[i]*prob->origobj_coeffs[env->best_sol.xind[i]];
      } else {
         continue;
      }
   }
   printf("\n\nORIGINAL OBJECTIVE FUNCTION VALUE = %06f\n\n", obj_val);


   CALL_FUNCTION( sym_close_environment(env) );
   
   return(0);
   
}

/*===========================================================================*\
\*===========================================================================*/

int user_read_data(user_problem *prob, char *infile)
{
   int j, counter;
   CoinMpsIO mps;

   mps.messageHandler()->setLogLevel(0);
   
   mps.setInfinity(mps.getInfinity()); // TODO: What exactly is this doing here?

   if (mps.readMps(infile,"")){
      return(USER_FUNC_ERROR);
   }
   
   prob->mip->m  = mps.getNumRows();
   prob->mip->n  = mps.getNumCols();
   prob->mip->obj_sense = 1.0; // TODO: Remove this 'minimize' assumption.
   prob->infty = mps.getInfinity();
   
   prob->mip->obj    = (double *) malloc(DSIZE * prob->mip->n);
   prob->mip->rhs    = (double *) malloc(DSIZE * prob->mip->m);
   prob->mip->sense  = (char *)   malloc(CSIZE * prob->mip->m);
   prob->mip->rngval = (double *) malloc(DSIZE * prob->mip->m);
   prob->mip->ub     = (double *) malloc(DSIZE * prob->mip->n);
   prob->mip->lb     = (double *) malloc(DSIZE * prob->mip->n);
   
   memcpy(prob->mip->obj, const_cast <double *> (mps.getObjCoefficients()),
	  DSIZE * prob->mip->n);
   memcpy(prob->mip->rhs, const_cast <double *> (mps.getRightHandSide()),
	  DSIZE * prob->mip->m);
   memcpy(prob->mip->sense, const_cast <char *> (mps.getRowSense()),
	  CSIZE * prob->mip->m);
   memcpy(prob->mip->rngval, const_cast <double *> (mps.getRowRange()),
	  DSIZE * prob->mip->m);
   memcpy(prob->mip->ub, const_cast <double *> (mps.getColUpper()),
	  DSIZE * prob->mip->n);
   memcpy(prob->mip->lb, const_cast <double *> (mps.getColLower()),
	  DSIZE * prob->mip->n);

   // Added by Suresh for debugging
   prob->origvar_num = prob->mip->n;
   prob->origobj_coeffs = (double *) malloc(DSIZE * prob->mip->n);
   memcpy(prob->origobj_coeffs, const_cast <double *> (mps.getObjCoefficients()),
	  DSIZE * prob->mip->n);

   //user defined matind, matval, matbeg--fill as column ordered
   const CoinPackedMatrix * matrixByCol= mps.getMatrixByCol();
   
   prob->mip->matbeg = (int *) malloc(ISIZE * (prob->mip->n + 1));
   memcpy(prob->mip->matbeg, const_cast<int *>(matrixByCol->getVectorStarts()),
	  ISIZE * (prob->mip->n + 1));
   
   prob->mip->matval = (double *) malloc(DSIZE*prob->mip->matbeg[prob->mip->n]);
   prob->mip->matind = (int *)    malloc(ISIZE*prob->mip->matbeg[prob->mip->n]);
   
   memcpy(prob->mip->matval, const_cast<double *> (matrixByCol->getElements()),
	  DSIZE * prob->mip->matbeg[prob->mip->n]);
   memcpy(prob->mip->matind, const_cast<int *> (matrixByCol->getIndices()), 
	  ISIZE * prob->mip->matbeg[prob->mip->n]);

   //user defined matind_row, matval_row, matbeg_row--fill as row ordered
   const CoinPackedMatrix * matrixByRow= mps.getMatrixByRow();
   
   prob->matbeg_row = (int *) malloc(ISIZE * (prob->mip->m + 1));
   memcpy(prob->matbeg_row, const_cast<int *>(matrixByRow->getVectorStarts()),
	  ISIZE * (prob->mip->m + 1));
   
   prob->matval_row = (double *) malloc(DSIZE*prob->matbeg_row[prob->mip->m]);
   prob->matind_row = (int *)    malloc(ISIZE*prob->matbeg_row[prob->mip->m]);
   
   memcpy(prob->matval_row, const_cast<double *> (matrixByRow->getElements()),
	  DSIZE * prob->matbeg_row[prob->mip->m]);
   memcpy(prob->matind_row, const_cast<int *> (matrixByRow->getIndices()), 
	  ISIZE * prob->matbeg_row[prob->mip->m]);
  
   //count number of different type constraints.
   prob->con_sense_e = 0;
   prob->con_sense_l = 0;
   prob->con_sense_g = 0;
   prob->ubinfty = 0;
   prob->lbinfty = 0;
   prob->infubind     = (int *) calloc(prob->mip->n, ISIZE);
   prob->inflbind     = (int *) calloc(prob->mip->n, ISIZE);
   prob->infubsofar   = (int *) malloc(ISIZE * prob->mip->n);
   prob->inflbsofar   = (int *) malloc(ISIZE * prob->mip->n);
   for (j = 0; j < prob->mip->m; j++) {
      if (prob->mip->sense[j] == 'E') {
         prob->con_sense_e++;
      } else if (prob->mip->sense[j] == 'L') {
         prob->con_sense_l++;
      } else if (prob->mip->sense[j] == 'G') {
         prob->con_sense_g++;
      } else {
         printf("\nNOOOOOOO!! ERROR!! ERROR!!\n");
      }
   }

   // count number of infinity UBs and -infinity LBs, and also their indicators.
   for (j = 0; j < prob->mip->n; j++) {
      prob->infubsofar[j] = prob->ubinfty;
      if (prob->mip->ub[j] >= prob->infty) {
         prob->ubinfty++;
         prob->infubind[j] = 1;
      }
      prob->inflbsofar[j] = prob->lbinfty;
      if (prob->mip->lb[j] <= -prob->infty) {
         prob->lbinfty++;
         prob->inflbind[j] = 1;
      }
   }

   prob->tempub     = (double *) malloc(DSIZE * (prob->mip->n - prob->ubinfty));
   prob->templb     = (double *) malloc(DSIZE * (prob->mip->n - prob->lbinfty));

   // Fill temporary UB arrays
   counter = 0;
   for (j = 0; j < prob->mip->n; j++) {
      if (prob->mip->ub[j] >= prob->infty) {
         continue;
      } else {
         prob->tempub[counter] = prob->mip->ub[j];
         counter++;
      }
   }

   // Fill temporary LB arrays
   counter = 0;
   for (j = 0; j < prob->mip->n; j++) {
      if (prob->mip->lb[j] <= -1.0 * prob->infty) {
         continue;
      } else {
         prob->templb[counter] = prob->mip->lb[j];
         counter++;
      }
   }

   return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*\
\*===========================================================================*/

int user_load_problem(sym_environment *env, user_problem *prob){
   
   int i, j, index, index1, n, m, nz, nz_index = 0, *column_starts, *matrix_indices;
   double *matrix_values, *lb, *ub, *obj, *rhs, *rngval;
   char *sense, *is_int, obj_sense = 1.0;
   
   /* set up the inital LP data */
   n = 3 * prob->mip->n + prob->mip->m - prob->ubinfty - prob->lbinfty;
   m = prob->mip->m + 3 * prob->mip->n - prob->ubinfty - prob->lbinfty;
   nz = 4 * prob->mip->n - 2 * prob->ubinfty - 2 * prob->lbinfty + 2 * prob->mip->matbeg[prob->mip->n];
   prob->colnum = n;
   prob->rownum = m;

   /* Allocate the arrays */
   column_starts  = (int *) malloc((n + 1) * ISIZE);
   matrix_indices = (int *) malloc((nz) * ISIZE);
   matrix_values  = (double *) malloc((nz) * DSIZE);
   obj            = (double *) malloc(n * DSIZE);
   lb             = (double *) malloc(n * DSIZE);
   ub             = (double *) malloc(n * DSIZE);
   rhs            = (double *) malloc(m * DSIZE);
   sense          = (char *) malloc(m * CSIZE);
   rngval         = (double *) calloc(m, DSIZE); /* TODO:Correct the assumption that this is zero always */
   is_int         = (char *) malloc(n * CSIZE);
   
   /* Fill out the appropriate data structures */
   if (prob->mip->obj_sense > 0.0) {
      for (i = 0; i < n; i++) {
         /*
         if (i < prob->mip->n) {
            obj[i] = prob->mip->obj[i];
         } else if (i < prob->mip->n + prob->mip->m) {
            obj[i] = -1.0 * prob->mip->rhs[i - prob->mip->n];
         } else if (i < 2 * prob->mip->n + prob->mip->m - prob->ubinfty) {
            obj[i] = -1.0 * prob->tempub[i - prob->mip->n - prob->mip->m];
         } else {
            obj[i] = -1.0 * prob->templb[i - 2 * prob->mip->n - prob->mip->m + prob->ubinfty];
         }
         */
         obj[i] = 0.0;
      }
   } else {
      for (i = 0; i < n; i++) {
         /*
         if (i < prob->mip->n) {
            obj[i] = -1.0 * prob->mip->obj[i];
         } else if (i < prob->mip->n + prob->mip->m) {
            obj[i] = prob->mip->rhs[i - prob->mip->n];
         } else if (i < 2 * prob->mip->n + prob->mip->m - prob->ubinfty) {
            obj[i] = prob->tempub[i - prob->mip->n - prob->mip->m];
         } else {
            obj[i] = prob->templb[i - 2 * prob->mip->n - prob->mip->m + prob->ubinfty];
         }
         */
         obj[i] = 0.0;
      }
   }
   for (i = 0; i < n; i++) {
      is_int[i] = FALSE;
      if (i < prob->mip->n) {
         ub[i] = prob->infty;
         lb[i] = -1.0 * prob->infty;
      } else if (i < prob->mip->n + prob->mip->m) {
         if (prob->mip->sense[i - prob->mip->n] == 'L') {
            ub[i] = 0;
            lb[i] = -1.0 * prob->infty;
         } else if (prob->mip->sense[i - prob->mip->n] == 'G') {
            ub[i] = prob->infty;
            lb[i] = 0;
         } else {
            ub[i] = prob->infty;
            lb[i] = -1.0 * prob->infty;
         }
      } else if (i < 2 * prob->mip->n + prob->mip->m - prob->ubinfty) {
         ub[i] = 0;
         lb[i] = -1.0 * prob->infty;
      } else {
         ub[i] = prob->infty;
         lb[i] = 0;
      }
   }
   for (i = 0; i < m; i++) {
      if (i < prob->mip->m) {
         sense[i] = prob->mip->sense[i];
         rhs[i] = prob->mip->rhs[i];
      } else if (i < (prob->mip->m + prob->mip->n - prob->ubinfty)) {
         sense[i] = 'L';
         rhs[i] = prob->tempub[i - prob->mip->m];
      } else if (i < (prob->mip->m + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty)) {
         sense[i] = 'G';
         rhs[i] = prob->templb[i - (prob->mip->m + prob->mip->n - prob->ubinfty)];
      } else {
         sense[i] = 'E';
         rhs[i] = prob->mip->obj[i - (prob->mip->m + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty)];
      }
   }
   column_starts[0] = 0;
   for (i = 0; i < prob->mip->n; i++) {
//      column_starts[i] = prob->mip->matbeg[i] + (2 - prob->infubind[i] - prob->inflbind[i]) * i;
      column_starts[i+1] = column_starts[i] + prob->mip->matbeg[i+1] - prob->mip->matbeg[i] + 2 - prob->infubind[i] - prob->inflbind[i];
      if (prob->mip->matbeg[i + 1] - prob->mip->matbeg[i] > 0) {
         for (j = (prob->mip->matbeg[i]); j < (prob->mip->matbeg[i + 1]); j++) {
            matrix_values[nz_index] = prob->mip->matval[j];
            matrix_indices[nz_index] = prob->mip->matind[j];
            nz_index++;
         }
      }
      if (!prob->infubind[i]) {
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + (i - prob->infubsofar[i]);
         nz_index++;
      }
      if (!prob->inflbind[i]) {
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + prob->mip->n - prob->ubinfty + (i - prob->inflbsofar[i]);
         nz_index++;
      }
   }
//   column_starts[prob->mip->n] = prob->mip->matbeg[prob->mip->n] + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty;
   index = 0;
   for (i = 0; i < prob->mip->m; i++) {
      column_starts[prob->mip->n + 1 + index] = column_starts[prob->mip->n + index] + prob->matbeg_row[i + 1] - prob->matbeg_row[i];
      index++;
      for (j = prob->matbeg_row[i]; j < prob->matbeg_row[i+1]; j++) {
         matrix_values[nz_index] = prob->matval_row[j];
         matrix_indices[nz_index] = prob->mip->m + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty + prob->matind_row[j];
         nz_index++;
      }
   }
   for (i = 0; i < prob->mip->n; i++) {
      if (!prob->infubind[i]) {
         column_starts[prob->mip->n + 1 + index] = column_starts[prob->mip->n + index] + 1;
         index++;
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty + (i - prob->infubsofar[i]);
         nz_index++;
      }
   }
   for (i = 0; i < prob->mip->n; i++) {
      if (!prob->inflbind[i]) {
         column_starts[prob->mip->n + 1 + index] = column_starts[prob->mip->n + index] + 1;
         index++;
         matrix_values[nz_index] = 1.0;
         matrix_indices[nz_index] = prob->mip->m + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty + (i - prob->inflbsofar[i]);
         nz_index++;
      }
   }

   /* Assign memory and fill complementarity constraint indices for variables */
   prob->ccind          =   (int *) malloc(n * ISIZE);
   index = 0;
   index1 = 0;
   for (i = 0; i < n;) {
      if (i < prob->mip->n) {
         prob->ccind[i] = prob->mip->m + 2 * prob->mip->n - prob->ubinfty - prob->lbinfty + i;
         i++;
      } else {
         prob->ccind[i] = index1;
         i++;
         index1++;
      }
   }

   /* Load the problem to SYMPHONY */   
   sym_explicit_load_problem(env, n, m, column_starts, matrix_indices,
			     matrix_values, lb, ub, is_int, obj, 0, sense, rhs,
			     rngval, true);

   /* Change prob->mip values to final problem values */
   prob->mip->n = n;
   prob->mip->m = m;
   prob->mip->nz = nz;
   prob->mip->obj_sense = obj_sense;

   prob->mip->obj    = (double *) realloc(prob->mip->obj, DSIZE * prob->mip->n);
   prob->mip->rhs    = (double *) realloc(prob->mip->rhs, DSIZE * prob->mip->m);
   prob->mip->sense  = (char *)   realloc(prob->mip->sense, CSIZE * prob->mip->m);
   prob->mip->rngval = (double *) realloc(prob->mip->rngval, DSIZE * prob->mip->m);
   prob->mip->ub     = (double *) realloc(prob->mip->ub, DSIZE * prob->mip->n);
   prob->mip->lb     = (double *) realloc(prob->mip->lb, DSIZE * prob->mip->n);
   prob->mip->is_int = (char *)   malloc(CSIZE * prob->mip->n);
   /* Default values for vvind, vvnum and feasible */
   prob->feasible    = USER__DO_NOT_BRANCH;
   prob->vvind       = (int *)    calloc(prob->mip->n, ISIZE);
   prob->vvnum       = 0;
   
   memcpy(prob->mip->obj, obj, DSIZE * prob->mip->n);
   memcpy(prob->mip->rhs, rhs, DSIZE * prob->mip->m);
   memcpy(prob->mip->sense, sense, CSIZE * prob->mip->m);
   memset(prob->mip->rngval, 0, DSIZE * n);                     // TODO: Fix this assumption.
   memcpy(prob->mip->ub, ub, DSIZE * prob->mip->n);
   memcpy(prob->mip->lb, lb, DSIZE * prob->mip->n);
   memcpy(prob->mip->is_int, is_int, CSIZE * prob->mip->n);

   prob->mip->matbeg = (int *) realloc(prob->mip->matbeg, ISIZE * (prob->mip->n + 1));
   memcpy(prob->mip->matbeg, column_starts, ISIZE * (prob->mip->n + 1));
   
   prob->mip->matval = (double *) realloc(prob->mip->matval, DSIZE*prob->mip->matbeg[prob->mip->n]);
   prob->mip->matind = (int *)    realloc(prob->mip->matind, ISIZE*prob->mip->matbeg[prob->mip->n]);
   
   memcpy(prob->mip->matval, matrix_values, DSIZE * prob->mip->matbeg[prob->mip->n]);
   memcpy(prob->mip->matind, matrix_indices, ISIZE * prob->mip->matbeg[prob->mip->n]);

   /* TODO: Delete mat*_row vectors here? */
   FREE(column_starts);
   FREE(matrix_indices);
   FREE(matrix_values);
   FREE(lb);
   FREE(ub);
   FREE(obj);
   FREE(sense);
   FREE(rhs);
   FREE(rngval);
   FREE(is_int);
   /* TODO: Is it good to free tempub, templb, infubind, inflbind, ifubsofar, inflbsofar here? */
   FREE(prob->tempub);
   FREE(prob->templb);
   FREE(prob->inflbsofar);
   FREE(prob->infubsofar);
   FREE(prob->inflbind);
   FREE(prob->infubind);
   FREE(prob->matval_row);
   FREE(prob->matind_row);
   FREE(prob->matbeg_row);

   return (FUNCTION_TERMINATED_NORMALLY);
}

#endif
