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

#include <stdio.h>

#include "BB_constants.h"
#include "BB_macros.h"
#include "lp_u.h"
#include "user.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions for the LP process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_lp_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_LP or COMPILE_IN_TM is FALSE. For 
 * sequential computation, nothing is needed here.
\*===========================================================================*/

int user_receive_lp_data(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must create the initial LP relaxation for
 * each search node. Basically, this involves constructing the base matrix in 
 * column ordered format. See the documentation for an explanation of how to 
 * fill out this function.
\*===========================================================================*/

int user_create_lp(void *user, int varnum, var_desc **vars, int rownum,
		   int cutnum, cut_data **cuts, int *nz, int **matbeg,
		   int **matind, double **matval, double **obj, double **rhs,
		   char **sense, double **rngval, int *maxn, int *maxm,
		   int *maxnz, int *allocn, int *allocm, int *allocnz)
{
   user_problem *prob = (user_problem *) user;
   int i, j, index;
   int resize;

   /* set up the inital LP data */

   *nz = 2 * varnum;

   /* We have to check to make sure there is enough space allocated
      for the matrix we are going to build */
   if (2 * rownum > *maxm){
      *maxm = 2 * rownum;
      resize = TRUE;
   }
   
   /* Allocate space for all edges up front since we have small problems */
   if (prob->colnum != *maxn){
      *maxn = prob->colnum;
      resize = TRUE;
   }
   
   if (*nz + ((*maxm) * (*maxn) / 10) > *maxnz){
      *maxnz = *nz + ((*maxm) * (*maxn) / 10);
      resize = TRUE;
   }
   
   /* If there was not enough space, the allocate more */
   if (resize){
      /*re-malloc all the arrays*/
      FREE(*matbeg);
      FREE(*matind);
      FREE(*matval);
      FREE(*obj);
      FREE(*rhs);
      FREE(*sense);
      FREE(*rngval);
      *allocm  = *maxm;
      *allocn  = *maxm + *maxn + 1;
      *allocnz = *maxnz + *maxm;
      *matbeg  = (int *) malloc(*allocn * ISIZE);
      *matind  = (int *) malloc(*allocnz * ISIZE);
      *matval  = (double *) malloc(*allocnz * DSIZE);
      *obj     = (double *) malloc(*allocn * DSIZE);
      *rhs     = (double *) malloc(*allocm * DSIZE);
      *sense   = (char *) malloc(*allocm * CSIZE);
      *rngval  = (double *) calloc(*allocm, DSIZE);
   }
   
   /* Fill out the appropriate data structures -- each column has
      exactly two entried*/
   index = 0;
   for (i = 0; i < prob->nnodes; i++) {
      for (j = i+1; j < prob->nnodes; j++) {
	 prob->node1[index] = i; /* The first node of assignment 'index' */
	 prob->node2[index] = j; /* The second node of assignment 'index' */
	 (*obj)[index] = prob->cost[i][j]; /* Cost of assignment (i, j) */
	 (*matbeg)[index] = 2*index;
	 (*matval)[2*index] = 1;
	 (*matval)[2*index+1] = 1;
	 (*matind)[2*index] = i;
	 (*matind)[2*index+1] = j;
	 index++;
      }
   }
   (*matbeg)[varnum] = 2*varnum;
   
   /* set the initial right hand side */
   for (i = 0; i < prob->nnodes; i++) {
      (*rhs)[i] = 1;
      (*sense)[i] = 'E';
   }

   return(USER_NO_PP);
}      


/*===========================================================================*/

/*===========================================================================*\
 * This function takes an LP solution and checks it for feasibility. By 
 * default, SYMPHONY checks for integrality. If any integral solution for your 
 * problem is feasible, then nothing needs to be done here.
\*===========================================================================*/

int user_is_feasible(void *user, double lpetol, int varnum, int *indices,
		     double *values, int *feasible, double *objval)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, the user can specify a special routine for sending back the feasible
 * solution. This need not be used unless there is a special format the user
 * wants the solution in. For sequential computation, you can use this routine
 * to interpret and store the feasible solution whenever one is found.
\*===========================================================================*/

int user_send_feasible_solution(void *user, double lpetol, int varnum,
				int *indices, double *values)
{
   return(DEFAULT);
}


/*===========================================================================*/

/*===========================================================================*\
 * This function graphically displays the current fractional solution
 * This is done using the Interactive Graph Drawing program, if it is used.
\*===========================================================================*/

int user_display_lp_solution(void *user, int which_sol, int varnum,
			     int *indices, double *values)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You can add whatever information you want about a node to help you
 * recreate it. I don't have a use for it, but maybe you will.
\*===========================================================================*/

int user_add_to_desc(void *user, int *desc_size, char **desc)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Compare cuts to see if they are the same. We use the default, which
 * is just comparing byte by byte.
\*===========================================================================*/

int user_same_cuts(void *user, cut_data *cut1, cut_data *cut2, int *same_cuts)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function receives a cut, unpacks it, and adds it to the set of
 * rows to be added to the LP. Only used if cutting planes are generated.
\*===========================================================================*/

int user_unpack_cuts(void *user, int from, int type, int varnum,
		     var_desc **vars, int cutnum, cut_data **cuts,
		     int *new_row_num, waiting_row ***new_rows)
{
   user_problem *prob = (user_problem *) user;
   
   int i, j, nzcnt;
   int *cutval;
   waiting_row **row_list;
   
   *new_row_num = cutnum;
   if (cutnum > 0)
      *new_rows =
	 row_list = (waiting_row **) calloc (cutnum*sizeof(waiting_row *));
   
   for (j = 0; j < cutnum; j++){
      row_list[j] = (waiting_row *) malloc(sizeof(waiting_row));
      switch (cuts[j]->type){
	 
      case TRIANGLE:
	 cutval = (int *) (cuts[j]->coef);
	 row_list[j]->cut = cuts[j];
	 row_list[j]->matind = (int *) malloc(varnum * ISIZE);
	 row_list[j]->matval = (double *) malloc(varnum * DSIZE);
	 row_list[j]->nzcnt = 0;
	 for (nzcnt = 0, i = 0; i < varnum; i++){
	    if (cutval[prob->node1[vars[i]->userind]] &&
		cutval[prob->node2[vars[i]->userind]]){
	       row_list[j]->matval[nzcnt] = 1.0;
	       row_list[j]->matind[nzcnt++] = vars[i]->userind;
	    }
	 }
	 row_list[j]->nzcnt = nzcnt;
	 break;

       default:
	 printf("Unrecognized cut type!\n");
      }
   }
   
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * If the user wants to fill in a customized routine for sending and receiving
 * the LP solution, it can be done here. For most cases, the default routines
 * are fine.
\*===========================================================================*/

int user_send_lp_solution(void *user, int varnum, var_desc **vars, double *x,
			  int where)
{
   return(SEND_NONZEROS);
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine does logical fixing of variables
\*===========================================================================*/

int user_logical_fixing(void *user, int varnum, var_desc **vars, double *x,
			char *status, int *num_fixed)
{
   *num_fixed = 0;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function generates the 'next' column. Only used for column generation.
\*===========================================================================*/

int user_generate_column(void *user, int generate_what, int cutnum,
			 cut_data **cuts, int prevind, int nextind,
			 int *real_nextind, double *colval, int *colind,
			 int *collen, double *obj)
{
   switch (generate_what){
    case GENERATE_NEXTIND:
      /* Here we just have to generate the specified column. */
      break;
    case GENERATE_REAL_NEXTIND:
      /* In this case, we have to determine what the "real" next edge is*/
      break;
   }

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to print some statistics on the types and quantities
 * of cuts or something like that.
\*===========================================================================*/

int user_print_stat_on_cuts_added(void *user, int rownum, waiting_row **rows)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to eliminate rows from the local pool based on
 * knowledge of problem structure.
\*===========================================================================*/

int user_purge_waiting_rows(void *user, int rownum, waiting_row **rows,
			    char *delete)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user has to generate the ubber bounds for the specified
 * variables. Lower bounds are always assumed (w.l.o.g.) to be zero.
\*===========================================================================*/

int user_get_upper_bounds(void *user, int varnum, int *indices, double *bd)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user might want to generate cuts in the LP using information
 * about the current tableau, etc. This is for advanced users only.
\*===========================================================================*/

int user_generate_cuts_in_lp(void *user, int varnum, var_desc **vars, double *x,
			     int *new_row_num, waiting_row ***new_rows)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Free all the user data structures
\*===========================================================================*/

int user_free_lp(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

