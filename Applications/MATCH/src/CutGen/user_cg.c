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

#include "BB_constants.h"
#include "BB_macros.h"
#include "cg_u.h"
#include "user.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains user-written functions used by the cut generator
 * process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_cg_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_CG, COMPILE_IN_LP, or COMPILE_IN_TM is
 * FALSE. For sequential computation, nothing is needed here.
\*===========================================================================*/

/*__BEGIN_EXPERIMENTAL_SECTION__*/
int user_receive_cg_data(void **user, int dg_id, int *varnum)
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
int user_receive_cg_data(void **user, int dg_id)
#endif
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * If the user wants to fill in a customized routine for sending and receiving
 * the LP solution, it can be done here. For most cases, the default routines
 * are fine.
\*===========================================================================*/

int user_receive_lp_solution_cg(void *user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Find cuts violated by a particular LP solution. This can be a fairly
 * involved function but the bottom line is that an LP solution comes in
 * and cuts go out. Remember, use the function cg_send_cut() to send cuts out
 * when they are found.
\*===========================================================================*/

/*__BEGIN_EXPERIMENTAL_SECTION__*/
int user_find_cuts(void *user, int varnum, int iter_num, int level,
		   int index, double objval, int *indices, double *values,
		   double ub, double etol, int *cutnum, char *status)
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
int user_find_cuts(void *user, int varnum, int iter_num, int level,
		   int index, double objval, int *indices, double *values,
		   double ub, double etol, int *cutnum)
#endif
{
   user_problem *prob = (user_problem *) user;
   double edge_val[200][200];
   int i, j, k, index;
   int *cuts;
   cut_data cut;
   
   *cutnum = 0;

   cuts = malloc(prob->nnodes*ISIZE);
   for (i = 0; i < prob->nnodes; i++) {
      for (j = 0; j < prob->nnodes; j++){
	 edge_val[i][j] = 0.0;
      }
   }
   
   for (index = 0; index < varnum; index++) {
      edge_val[prob->node1[indices[index]]][prob->node2[indices[index]]] 
	 = values[index];
   }
   for (i = 0; i < prob->nnodes; i++){
      for (j = i+1; j < prob->nnodes; j++){
	 for (k = j+1; k < prob->nnodes; k++) {
	    if (edge_val[i][j]+edge_val[j][k]+edge_val[i][k] > 1.0 + etol) {
	       for (index = 0; index < prob->nnodes; index++) {
		  cuts[index] = 0;
	       }
	       cuts[i] = 1; 
	       cuts[j] = 1;
	       cuts[k] = 1;
	       cut.size = (prob->nnodes)*ISIZE;
	       cut.coef = (char *) cuts;
	       cut.rhs = 1.0;
	       cut.range = 0.0;
	       cut.type = 1;
	       cut.sense = 'G';
	       cut.deletable = TRUE;
	       cut.branch = ALLOWED_TO_BRANCH_ON;
	       cg_send_cut(&cut);
	       (*cutnum)++;
	       
	    }
	 }
      }
   }
   
   return(USER_NO_PP);

}

/*===========================================================================*/

/*===========================================================================*\
 * Free the user data structure. If the default setup is used with sequential
 * computation, nothing needs to be filled in here.
\*===========================================================================*/

int user_free_cg(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is an undocumented (for now) debugging feature which can allow the user
 * to identify the cut which cuts off a particular known feasible solution.
\*===========================================================================*/

#ifdef CHECK_CUT_VALIDITY
int user_check_validity_of_cut(void *user, cut_data *new_cut)
{
  return(USER_NO_PP);
}
#endif
