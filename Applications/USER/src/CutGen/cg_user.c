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

/*===========================================================================*/

/*===========================================================================*\
 * This file contains user-written functions used by the cut generator
 * process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_cg_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_CG, COMPILE_IN_LP, or COMPILE_IN_TM is
 * FALSE.
\*===========================================================================*/

/*__BEGIN_EXPERIMENTAL_SECTION__*/
int user_receive_cg_data(void **user, int dg_id, int *varnum)
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
int user_receive_cg_data(void **user, int dg_id)
#endif
{
   *user = NULL;
   
   return(USER_NO_PP);
}

/*===========================================================================*/

int user_receive_lp_solution_cg(void *user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Free the user data structure
\*===========================================================================*/

int user_free_cg(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Find cuts violated by a particular LP solution. This can be a fairly
 * involved function but the bottom line is that an LP solution comes in
 * and cuts go out. Remember, use the function cg_send_cut() to send cuts out
 * when they are found
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
   *cutnum = 0;

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
