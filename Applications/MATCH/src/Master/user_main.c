/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2004 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
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

#include "symphony_api.h"
#include "user.h"
#include <stdlib.h>

int main(int argc, char **argv)
{

   int termcode;
   char * infile;

   /* Create a SYMPHONY environment */
   sym_environment *env = sym_open_environment();

   /* Create a user problem structure to read in the data and then pass it to  
      SYMPHONY. 
   */
   user_problem *prob = (user_problem *)calloc(1, sizeof(user_problem));

   CALL_FUNCTION( sym_parse_command_line(env, argc, argv) );

   CALL_FUNCTION( sym_get_str_param(env, "infile_name", &infile));

   CALL_FUNCTION( match_read_data(env, (void *) prob, infile));
   
   CALL_FUNCTION( match_load_problem(env, (void *) prob ));

   CALL_FUNCTION( sym_solve(env) );

   CALL_FUNCTION( sym_close_environment(env) );

   return(0);

}

int match_read_data(sym_environment *env, void *user, char *infile)
{
   int i, j;
   FILE *f = NULL;
   /* This gives you access to the user data structure. */
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

   /* This will pass the user data in to SYMPHONY*/
   sym_set_user_data(env, (void *)prob);

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
   
   /* Load the problem to SYMPHONY */   
   sym_explicit_load_problem(env, n, m, matbeg, matind, matval, lb, ub, 
			     is_int, obj, 0, sense, rhs, rngval, true);
			     
   FREE(matbeg);
   FREE(matind);
   FREE(matval);
   FREE(lb);
   FREE(ub);
   FREE(obj);
   FREE(sense);
   FREE(rhs);
   FREE(rngval);

   return (FUNCTION_TERMINATED_NORMALLY);

}

