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
#include "master_u.h"
#ifdef COMPILE_IN_TM
#ifdef COMPILE_IN_LP
/* fill these in for sequentail compilation */
#ifdef COMPILE_IN_CG
/* fill these in for sequentail compilation */
#endif
#ifdef COMPILE_IN_CP
/* fill these in for sequentail compilation */
#endif
#endif
#endif

/*===========================================================================*/

void user_usage(void){
         printf("master [ -H ] \n\n\t-H: help");
}

/*===========================================================================*\
 * This file contains the user-written functions for the master process.
\*===========================================================================*/

/*===========================================================================*\
 * Initialize user-defined data structures. In this case, I store all
 * problem-specific data such as the location of the customers, edge costs,
 * etc. in this data-structure.
\*===========================================================================*/

int user_initialize(void **user)
{
   *user = NULL;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Read in parameters from the parameter file given on the command line.
\*===========================================================================*/

int user_readparams(void *user, char *filename, int argc, char **argv)
{
   FILE *f;
   char line[50], key[50], value[50];
   int c;
   
   if ((f = fopen(filename, "r")) == NULL){
      printf("SYMPHONY: file %s can't be opened\n", filename);
      exit(1); /*error check for existence of parameter file*/
   }
   
   while(NULL != fgets(line, 50, f)){  /*read in parameter settings*/
      strcpy(key, "");
      sscanf(line, "%s%s", key, value);
   }      

   while ((c = getopt(argc, argv, "H")) != -1){
      switch (c) {
       case 'H':
	 user_usage();
	 exit(0);
	 break;
      };
   }

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Read in the data file, whose name was given in the parameter file.
 * This file contains instance data.
\*===========================================================================*/

int user_io(void *user)
{
   return(USER_NO_PP);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * Here is where the heuristics are performed and an upper bound is calculated.
 * An upper bound can also be specified in the parameter file. 
\*===========================================================================*/

int user_start_heurs(void *user, double *ub, double *ub_estimate)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * If graph drawing will be use, the user must initialize the drawing
 * window here.
\*===========================================================================*/

int user_init_draw_graph(void *user, int dg_id)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * In this routine, I build the initial edge set for the root. There are
 * several things going on here. First, there is a user-defined parameter
 * defining whether or not to just go ahead and add all variables to the
 * problem up front (vrp->par.add_all_edges). Currently, this seems to be the
 * best option since the problems are small anyway. Further, I am doing some
 * preprocessing here by eliminating edges for which the sum of the demands of
 * their endpoints is greater than the capacity since these edges cannot
 * be in any feasible solution.
 *
 * Notice that there are several options programmed for which set
 * of edges should be in the base set. The
 * base constraints are just the degree constraints from the IP
 * formulation. These do not have to be specified explicitly, just the
 * number of them given.
\*===========================================================================*/

int user_set_base(void *user, int *basevarnum, int **basevars, double **lb,
		  double **ub, int *basecutnum, int *colgen_strat)
{
   *basevarnum = 0;
   *basevars = NULL;
   *basecutnum = 0;

   return(USER_NO_PP);
}

/*===========================================================================*\
 * This is the second step in the process where the user specifies
 * which variables should be active in the root in addition to the base
 * set specified above
\*===========================================================================*/

/*===========================================================================*/

int user_create_root(void *user, int *extravarnum, int **extravars)
{
   *extravarnum = 0;
   *extravars  = NULL;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Receive the feasible solution
\*===========================================================================*/

int user_receive_feasible_solution(void *user, int msgtag, double cost,
				   int numvars, int *indices, double *values)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the LP process. Notice that
 * there are two cases to deal with. If the LP or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. 
\*===========================================================================*/

int user_send_lp_data(void *user, void **user_lp)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* This is is the case when we are copying data directly because
      the LP is not running separately. This code should be virtually
      identical to that of user_receive_lp_data() in the LP process.*/
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_lp_data() in the LP process */
#endif
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the CG process. Notice that
 * there are two cases to deal with. If the CG, LP, or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. 
\*===========================================================================*/

int user_send_cg_data(void *user, void **user_cg)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CG)
   /* This is is the case when we are copying data directly because
      the CG is not running separately. This code should be virtually
      identical to that of user_receive_cg_data() in the CG process.*/
#ifdef CHECK_CUT_VALIDITY
   /* Send the feasible solution here */
#endif
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cg_data() in the CG process */
#ifdef CHECK_CUT_VALIDITY
   /* Send the feasible solution here */
#endif
#endif
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the CP process. Notice that
 * there are two cases to deal with. If the CP, LP, or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. 
\*===========================================================================*/

int user_send_cp_data(void *user, void **user_cp)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CP)
   /* This is is the case when we are copying data directly because
      the LP is not running separately. This code should be virtually
      identical to that of user_receive_cp_data() in the CP process.*/
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cp_data() in the CP process */
#endif
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Generally, this function is not needed but you might find some use
 * for it. Someone did :).
\*===========================================================================*/

int user_process_own_messages(void *user, int msgtag)
{
   switch (msgtag){
    default:
      fprintf(stderr, "\nMaster: unknown message type %i!!!\n\n", msgtag);
      exit(1);
   }

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the user's chance to display the solution in whatever
 * manner desired. 
\*===========================================================================*/

int user_display_solution(void *user)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* In this case, the LP user data structure is passed in */
#else
   /* In this case, it is the master user data structure */
#endif
   return(USER_NO_PP);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * This is a debugging feature which might
 * allow you to find out why a known feasible solution is being cut off.
\*===========================================================================*/

int user_send_feas_sol(void *user, int *feas_sol_size, int **feas_sol)
{
#ifdef TRACE_PATH

#endif
   return(USER_NO_PP);
}   

/*===========================================================================*/

/*===========================================================================*\
 * This function frees everything.
\*===========================================================================*/

int user_free_master(void **user)
{
   return(USER_NO_PP);
}






