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

/* system include files */
#include <stdio.h>
#include <math.h>

/* SYMPHONY include files */
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_lp_u.h"

/* User include files */
#include "user.h"

#define TOL 1e-5
double cone_feasibility(int type, int size, int * members, double * sol);
int integer_feasibility(void *user, double lpetol, int varnum, int *indices,
			double *values, int *feasible, double *objval,
			char branching, double *heur_solution);

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions for the LP process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_lp_data() and set up data structures. Note that this function is
 * only called if one of SYM_COMPILE_IN_LP or SYM_COMPILE_IN_TM is FALSE. For
 * sequential computation, nothing is needed here.
\*===========================================================================*/

int user_receive_lp_data(void **user)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must create the initial LP relaxation for
 * each search node. Basically, this involves constructing the base matrix in
 * column ordered format. See the documentation for an explanation of how to
 * fill out this function.
\*===========================================================================*/

int user_create_subproblem(void *user, int *indices, MIPdesc *mip,
			   int *maxn, int *maxm, int *maxnz)
{
   user_problem *prob = (user_problem *) user;

   return(USER_DEFAULT);
}


/*===========================================================================*/

/*===========================================================================*\
 * This function takes an LP solution and checks it for feasibility. By
 * default, SYMPHONY checks for integrality. If any integral solution for your
 * problem is feasible, then nothing needs to be done here.
\*===========================================================================*/

int user_is_feasible(void *user, double lpetol, int varnum, int *indices,
		     double *values, int *feasible, double *objval,
		     char branching, double *heur_solution) {
  //todo(aykut) fill this
  user_problem *prob = (user_problem *) user;
  int num_cones = prob->num_cones;
  int * type = prob->cone_type;
  int * size = prob->cone_size;
  int ** members = prob->cone_members;
  // check feasibility of conic constraints
  int i;
  double feas;
  *feasible = IP_FEASIBLE;
  // allocate solution
  free(prob->curr_solution);
  int colnum = prob->colnum;
  prob->curr_solution = (double *) calloc(colnum, sizeof(double));
  for (i=0; i<varnum; ++i) {
    prob->curr_solution[indices[i]] = values[i];
  }
  for (i=0; i<num_cones; ++i) {
    feas = cone_feasibility(type[i], size[i], members[i], prob->curr_solution);
    if (feas<-TOL) {
      *feasible = IP_INFEASIBLE;
      return USER_SUCCESS;
      //break;
    }
  }
  // check integrality
  int integrality;
  // returns true if feasible
  integrality = integer_feasibility(user, lpetol, varnum, indices, values,
				    feasible, objval, branching,
				    heur_solution);
  if (!integrality)
    *feasible = IP_INFEASIBLE;
  return USER_SUCCESS;
}

// returns true if feasible
int integer_feasibility(void *user, double lpetol, int varnum, int *indices,
		     double *values, int *feasible, double *objval,
		     char branching, double *heur_solution) {
  int i;
  user_problem *prob = (user_problem *) user;
  for (i=0; i<varnum; ++i) {
    if (prob->is_int[indices[i]]) {
      double floor_i = floor(values[i]);
      double ceil_i = floor_i + 1;
      if (((values[i]-floor_i)>TOL) && ((ceil_i-values[i])>TOL))
	return 0;
    }
  }
  return 1;
}

double cone_feasibility(int type, int size, int * members, double * sol) {
  int i;
  double feas = 0.0;
  if (type==0) {
    double term1 = sol[members[0]];
    double term2 = 0.0;
    for (i=1; i<size; ++i) {
      term2 += sol[members[i]]*sol[members[i]];
    }
    feas = term1 - sqrt(term2);
  }
  else {
    double term1 = 2.0*sol[members[0]]*sol[members[1]];
    double term2 = 0.0;
    for (i=2; i<size; ++i) {
      term2 += sol[members[i]]*sol[members[i]];
    }
    feas = term1 - term2;
  }
  return feas;
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, the user can specify a special routine for sending back the feasible
 * solution. This need not be used unless there is a special format the user
 * wants the solution in. This function is only called in sequential mode.
\*===========================================================================*/

int user_send_feasible_solution(void *user, double lpetol, int varnum,
				int *indices, double *values)
{
   return(USER_DEFAULT);
}


/*===========================================================================*/

/*===========================================================================*\
 * This function graphically displays the current fractional solution
 * This is done using the Interactive Graph Drawing program, if it is used.
\*===========================================================================*/

int user_display_lp_solution(void *user, int which_sol, int varnum,
			     int *indices, double *values)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You can add whatever information you want about a node to help you
 * recreate it. I don't have a use for it, but maybe you will.
\*===========================================================================*/

int user_add_to_desc(void *user, int *desc_size, char **desc)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Compare cuts to see if they are the same. We use the default, which
 * is just comparing byte by byte.
\*===========================================================================*/

int user_same_cuts(void *user, cut_data *cut1, cut_data *cut2, int *same_cuts)
{
   return(USER_DEFAULT);
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
   /* This code is just here as a template for customization. Uncomment to use.*/
#if 0
   int j;

   user_problem *prob = (user_problem *) user;

   *new_row_num = 0;
   for (j = 0; j < cutnum; j++){
      switch (cuts[j]->type){

       default:
	 printf("Unrecognized cut type!\n");
      }
   }
#endif

   return(USER_DEFAULT);
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
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine does logical fixing of variables
\*===========================================================================*/

int user_logical_fixing(void *user, int varnum, var_desc **vars, double *x,
			char *status, int *num_fixed)
{
   *num_fixed = 0;

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function generates the 'next' column. Only used for column generation.
\*===========================================================================*/

int user_generate_column(void *user, int generate_what, int cutnum,
			 cut_data **cuts, int prevind, int nextind,
			 int *real_nextind, double *colval, int *colind,
			 int *collen, double *obj, double *lb, double *ub)
{
   switch (generate_what){
    case GENERATE_NEXTIND:
      /* Here we just have to generate the specified column. */
      break;
    case GENERATE_REAL_NEXTIND:
      /* In this case, we have to determine what the "real" next edge is*/
      break;
   }

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to print some statistics on the types and quantities
 * of cuts or something like that.
\*===========================================================================*/

int user_print_stat_on_cuts_added(void *user, int rownum, waiting_row **rows)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to eliminate rows from the local pool based on
 * knowledge of problem structure.
\*===========================================================================*/

int user_purge_waiting_rows(void *user, int rownum, waiting_row **rows,
			    char *delete_rows)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user might want to generate cuts in the LP using information
 * about the current tableau, etc. This is for advanced users only.
\*===========================================================================*/

int user_generate_cuts_in_lp(void *user, LPdata *lp_data, int varnum,
			     var_desc **vars, double *x,
			     int *new_row_num, cut_data ***cuts)
{
   return(GENERATE_CGL_CUTS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Free all the user data structures
\*===========================================================================*/

int user_free_lp(void **user)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

