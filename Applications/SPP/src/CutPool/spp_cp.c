/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* Set Partitioning Problems.                                                */
/*                                                                           */
/* (c) Copyright 2003 Marta Eso and Ted Ralphs. All Rights Reserved.         */
/*                                                                           */
/* This application was originally developed by Marta Eso and was modified   */
/* Ted Ralphs (tkralphs@lehigh.edu)                                          */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/* system include files */
#include <stdio.h>

/* SYMPHONY include files */
#include "BB_constants.h"
#include "cp_u.h"

/* SPP include files */
#include "spp.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions of the cut pool process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_cp_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_CP, COMPILE_IN_LP, or COMPILE_IN_TM is
 * FALSE. For sequential computation, nothing is needed here.
\*===========================================================================*/

int user_receive_cp_data(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * If the user wants to fill in a customized routine for sending and receiving
 * the LP solution, it can be done here. For most cases, the default routines
 * are fine.
\*===========================================================================*/

int user_receive_lp_solution_cp(void *user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * When a new solution arrives to the cut pool, this function is invoked
 * so that the user can prepare for checking many cuts (probably set up
 * some data structures that make ckecking more efficient). 
\*===========================================================================*/

int user_prepare_to_check_cuts(void *user, int varnum, int *indices,
			       double *values)
{
   return(USER_NO_PP);
}


/*===========================================================================*/

/*===========================================================================*\
 * Check to see whether a particular cut is violated by the current LP sol.
\*===========================================================================*/
      
int user_check_cut(void *user, double etol, int number, int *indices,
		   double *values, cut_data *cut, int *is_violated,
		   double *quality)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function is invoked when all cuts that needed to be checked for
 * the current solution have been checked already. (Disassemble the
 * data structures built up in 'user_prepare_to_check_cuts'.
\*===========================================================================*/

int user_finished_checking_cuts(void *user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we free up the data structures. If the default setup is used with 
 * sequential computation, nothing needs to be filled in here.
\*===========================================================================*/

int user_free_cp(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

