/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER

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

#if defined(USE_OSI_INTERFACE) && !defined(USER_MAIN) 

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

#elif !defined(USER_MAIN)

#include "symphony_api.h"
#include "user.h"
#include <stdlib.h>

int main(int argc, char **argv)
{
   int termcode;
   user_problem *prob = (user_problem *)calloc(1, sizeof(user_problem));

   /* FIXME! sym_open_environment calls user_initialize! */
   sym_environment *env = sym_open_environment();

   CALL_FUNCTION( sym_parse_command_line(env, argc, argv) );

   CALL_FUNCTION( match_read_data((void *) prob, env->par.infile));

   CALL_FUNCTION( sym_set_user_data(env, (void *)prob));
   
   CALL_FUNCTION( match_load_problem(env, (void *) prob ));

   CALL_FUNCTION( sym_find_initial_bounds(env) );

   CALL_FUNCTION( sym_solve(env) );

   CALL_FUNCTION( sym_close_environment(env) );

   return(0);
}

int match_read_data(void *user, char *infile)
{
   int i, j;
   FILE *f = NULL;
   user_problem *prob = (user_problem *) user;

   if ((f = fopen(infile, "r")) == NULL){
      printf("main(): user file %s can't be opened\n", infile);
      return(ERROR__USER); 
   }

   /* Read in the costs */
   fscanf(f,"%d",&(prob->nnodes));
   for (i = 0; i < prob->nnodes; i++)
      for (j = 0; j < prob->nnodes; j++)
	 fscanf(f, "%d", &(prob->cost[i][j]));
   
   prob->colnum = (prob->nnodes)*(prob->nnodes-1)/2;
   prob->rownum = prob->nnodes;

   return (FUNCTION_TERMINATED_NORMALLY);
}


int match_load_problem(sym_environment *env, void *user){
   
   int i, j, index, n, m, nz, *matbeg, *matind;
   double *matval, *lb, *ub, *obj, *rhs, *rngval;
   char *sense, *is_int;
   user_problem *prob = (user_problem *) user;

   /* set up the inital LP data */
   n = prob->colnum;
   m = prob->rownum;
   nz = 2 * n;

   /*  Estimate the maximum number of nonzeros 
   *maxm = 2 * mip->m;
   *maxn = mip->n;
   *maxnz = mip->nz + ((*maxm) * (*maxn) / 10);
   */   

   /* Allocate the arrays */
   matbeg  = (int *) malloc((n + 1) * ISIZE);
   matind  = (int *) malloc((nz) * ISIZE);
   matval  = (double *) malloc((nz) * DSIZE);
   obj     = (double *) malloc(n * DSIZE);
   lb      = (double *) calloc(n, DSIZE);
   ub      = (double *) malloc(n * DSIZE);
   rhs     = (double *) malloc(m * DSIZE);
   sense   = (char *) malloc(m * CSIZE);
   rngval  = (double *) calloc(m, DSIZE);
   is_int  = (char *) malloc(n * CSIZE);
   
   /* Fill out the appropriate data structures -- each column has
      exactly two entries */
   index = 0;
   for (i = 0; i < prob->nnodes; i++) {
      for (j = i+1; j < prob->nnodes; j++) {
	 prob->node1[index] = i; /* The first node of assignment 'index' */
	 prob->node2[index] = j; /* The second node of assignment 'index' */
	 obj[index] = prob->cost[i][j]; /* Cost of assignment (i, j) */
	 is_int[index] = TRUE;
	 matbeg[index] = 2*index;
	 matval[2*index] = 1;
	 matval[2*index+1] = 1;
	 matind[2*index] = i;
	 matind[2*index+1] = j;
	 ub[index] = 1.0;
	 index++;
      }
   }
   matbeg[n] = 2 * n;
   
   /* set the initial right hand side */
   for (i = 0; i < prob->nnodes; i++) {
      rhs[i] = 1;
      sense[i] = 'E';
   }
   
   sym_explicit_load_problem(env, n, m, matbeg, matind, matval, lb, ub, 
			     is_int, obj, 0, sense, rhs, rngval, true);
			     
#if 0
   for(i = 0; i<n; i++){
      if (is_int[i]){
	 sym_set_integer(env, i);
	 }
   }
#endif

   FREE(matbeg);
   FREE(matind);
   FREE(matval);
   FREE(lb);
   FREE(ub);
   FREE(obj);
   FREE(sense);
   FREE(rhs);
   FREE(rngval);

	return(0);

}

#endif
