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
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "BB_constants.h"
#include "BB_macros.h"
#include "master_u.h"
#include "user.h"
#ifdef COMPILE_IN_TM
#ifdef COMPILE_IN_LP
/* fill these in for sequential compilation if needed. */
#ifdef COMPILE_IN_CG
/* fill these in for sequential compilation if needed. */
#endif
#ifdef COMPILE_IN_CP
/* fill these in for sequential compilation if needed. */
#endif
#endif
#endif

/*===========================================================================*\
 * This file contains stubs for the user-written functions for the master 
 * process. The primary function that has to be filled in here is user_io(),
 * where the data for the instance is read in and the user data structure
 * that stores the instance data filled out (this data structure is defined 
 * in user.h). Other than that, the default routines should work fine.
\*===========================================================================*/

/*===========================================================================*/

/*===========================================================================*\
 * This function gives help on command-line switches defined by the user.
 * All user switches have capital letters by convention.
\*===========================================================================*/

void user_usage(void){
         printf("master [ -H ] [ -F file ] \n\t%s\n\t%s\n",
		"-H: help (user switches)",
		"-F file: problem instance data is in 'file'");
}

/*===========================================================================*/

/*===========================================================================*\
 * Initialize user-defined data structures. This basically consists of 
 * allocating the memory. If you are using the default data structure,
 * nothing needs to be changed here.
\*===========================================================================*/

int user_initialize(void **user)
{
   user_problem *prob = (user_problem *) calloc(1, sizeof(user_problem));

   *user = prob;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Parse the user options Read in parameters from the parameter file given on the command line
\*===========================================================================*/

int user_readparams(void *user, char *filename, int argc, char **argv)
{
   FILE *f;
   char line[50], key[50], value[50], c, tmp;
   int i;
   /* This gives you access to the user data structure*/
   user_problem *prob = (user_problem *) user;
   user_parameters *par = &(prob->par);
   
   if ((f = fopen(filename, "r")) == NULL){
      printf("SYMPHONY: file %s can't be opened\n", filename);
      exit(1); /*error check for existence of parameter file*/
   }
   
   /* Here you can read in the parameter settings from the file. See the 
      function bc_readparams() for an example of how this is done. */
   while(NULL != fgets(line, MAX_LINE_LENGTH, f)){  /*read in parameters*/
      strcpy(key, "");
      sscanf(line, "%s%s", key, value);

      if (strcmp(key, "input_file") == 0){
	 par->infile[MAX_FILE_NAME_LENGTH] = 0;
	 strncpy(par->infile, value, MAX_FILE_NAME_LENGTH);
      }
   }      

   fclose(f);
   
   /* Here you can parse the command line for options. By convention, the
      users options should be capital letters */

   for (i = 1; i < argc; i++){
      sscanf(argv[i], "%c %c", &tmp, &c);
      if (tmp != '-')
	 continue;
      switch (c) {
       case 'H':
	 user_usage();
	 exit(0);
	 break;
       case 'F':
	 strncpy(par->infile, argv[++i], MAX_FILE_NAME_LENGTH);
	 break;
      };
   }

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Read in the data file, whose name was given in the parameter file.
 * This file contains instance data. Right now, this function is set up to 
 * read in just the number of columns and number of rows from the file.
 * Add more data as needed to describe the instance and set up the LP
 * relaxation.
\*===========================================================================*/

int user_io(void *user)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;
   user_parameters *par = &(prob->par);
   char *infile = par->infile;
   FILE *f = NULL;
   char line[MAX_LINE_LENGTH], key[50], value[50];

   if ((f = fopen(infile, "r")) == NULL){
      printf("Readparams: file %s can't be opened\n", infile);
      exit(1); /*error check for existence of parameter file*/
   }

   /* Here you can read in the data for the problem instance. For the default
      setup, the user should set the colnum and rownum here. */
   while(NULL != fgets( line, MAX_LINE_LENGTH, f)){  /*read in problem data*/
      strcpy(key, "");
      sscanf(line, "%s%s", key, value);
      if (strcmp(key, "colnum") == 0){ /* Read in the number of rows */
	 READ_INT_PAR(prob->colnum);
      }
      else if (strcmp(key, "rownum") == 0){ /* Read in the number of columns */
	 READ_INT_PAR(prob->rownum);
      }
   }

   return(USER_NO_PP);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * Here is where the heuristics are performed and an upper bound is calculated.
 * An upper bound can also be specified in the parameter file. This function
 * need not be filled in if no upper bounding is done.
\*===========================================================================*/

int user_start_heurs(void *user, double *ub, double *ub_estimate)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   *ub = MAXINT;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * If graph drawing will be used, the user must initialize the drawing
 * window here. This function need not be filled in.
\*===========================================================================*/

int user_init_draw_graph(void *user, int dg_id)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the subroutine where the user specifies what variables are to be in
 * the base set. To begin with, a good bet is just to put all the variables in
 * the base set. In this case, this function need not be modified.
\*===========================================================================*/

int user_set_base(void *user, int *basevarnum, int **basevars, double **lb,
		  double **ub, int *basecutnum, int *colgen_strat)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;
   int i;
   int *vars, varnum;

   /* Set the number of variables*/
   varnum = *basevarnum = prob->colnum;
 
   /* Allocate memory for the uper and lower bounds. */
   /* Lower bounds are (probably) all zero so calloc those. */
   *lb = (double *) calloc (varnum, DSIZE);
   *ub = (double *) malloc (varnum * DSIZE);

   /* This puts all the variable in the base set and fills out the 
      upper bounds */
   vars = *basevars = (int *) malloc(varnum * ISIZE);
   for (i = 0; i < varnum; i++){
     vars[i] = i;
     (*ub)[i] = 1; /* If the upper bounds are not 1, change this line. */
     /* (*lb)[i] = 0; /* If the lower bounds are not 0, uncomment this line. */
   }

   /* Set the number of rows in the base */
   *basecutnum = prob->rownum;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the second step in the process, where the user specifies
 * which variables should be active in the root in addition to the base
 * set specified above. The set of extra variable would be empty if all
 * variables are in the base, as above.
\*===========================================================================*/

int user_create_root(void *user, int *extravarnum, int **extravars)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   *extravarnum = 0;
   *extravars  = NULL;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Receive the feasible solution. Doesn't need to be filled in.
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
 * configurations. If running sequentially and using the default data
 * structure, nothing needs to be modified in here.
\*===========================================================================*/

int user_send_lp_data(void *user, void **user_lp)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* This is is the case when we are copying data directly because the LP is
      not running separately. The easiest thing to do here is just to use the
      same user data structure in both the master and the LP. Then this
      subroutine would simply consist of the line
      
      *user_lp = user;

      Otherwise, this code should be virtually
      identical to that of user_receive_lp_data() in the LP process.*/

   *user_lp = user;
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
 * configurations. If running sequentially and using the default data
 * structure, nothing needs to be modified in here.
\*===========================================================================*/

int user_send_cg_data(void *user, void **user_cg)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CG)
   /* This is is the case when we are copying data directly because
      the CG is not running separately. The easiest thing to do here is just
      to use the same user data structure in both the master and the cut
      generator. Then this subroutine would simply consist of 
      
      *user_cg = user;

      Otherwise, this code should be virtually
      identical to that of user_receive_cg_data() in the CG process.*/

   *user_cg = user;
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
 * configurations. If running sequentially and using the default data
 * structure, nothing needs to be modified in here.
\*===========================================================================*/

int user_send_cp_data(void *user, void **user_cp)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CP)
   /* This is is the case when we are copying data directly because
      the CP is not running separately. The easiest thing to do here is just
      to use the same user data structure in both the master and the cut
      pool. Then this subroutine would simply consist of 
      
      *user_cp = user;

      Otherwise, this code should be virtually
      identical to that of user_receive_cp_data() in the CP process.*/

   *user_cp = user;
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
 * manner desired. Change the return value to USER_NO_PP if you want to
 * display the solution yourself. A return value of USER_AND_PP will cause the
 * default solution display routine to be executed, even if the user displays
 * the solution as well.
\*===========================================================================*/

int user_display_solution(void *user, double lpetol, int varnum, int *indices,
			  double *values, double objval)
{
   return(DEFAULT);
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
   user_problem *prob = (user_problem *) calloc(1, sizeof(user_problem));

   FREE(prob);

   return(USER_NO_PP);
}






