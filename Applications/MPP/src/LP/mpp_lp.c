/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Mixed Postman Problem.                                                */
/*                                                                           */
/* (c) Copyright 2003 Lehigh University. All Rights Reserved.                */
/*                                                                           */
/* This application was originally developed by Andrew Hofmann and was       */
/* modified by  Ted Ralphs (tkralphs@lehigh.edu)                             */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/* system include files */
#include <stdio.h>
#include <malloc.h>

/* SYMPHONY include files */
#include "BB_constants.h"
#include "BB_macros.h"
#include "lp_u.h"

/* MPP include files */
#include "mpp.h"

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

int user_create_lp(void *user, LPdesc *desc, int *indices, 
		   int *maxn, int *maxm, int *maxnz)
{
   mpp_problem *mpp = (mpp_problem *) user;

   int i, j, ind;
   char first;

   /* set up the inital LP data */

   desc->nz = (6 * mpp->numedges)+ (2 * mpp->numarcs);
 
   /* Estimate the maximum number of nonzeros */
   *maxm = 2 * desc->m;
   *maxn = desc->n;
   *maxnz = desc->nz + ((*maxm) * (*maxn) / 10);

   /* Allocate the arrays. These are owned by SYMPHONY after returning. */
   desc->matbeg  = (int *) malloc((desc->n + 1) * ISIZE);
   desc->matind  = (int *) malloc((desc->nz) * ISIZE);
   desc->matval  = (double *) malloc((desc->nz) * DSIZE);
   desc->obj     = (double *) malloc(desc->n * DSIZE);
   desc->lb      = (double *) calloc(desc->n, DSIZE);
   desc->ub      = (double *) malloc(desc->n * DSIZE);
   desc->rhs     = (double *) malloc(desc->m * DSIZE);
   desc->sense   = (char *) malloc(desc->m * CSIZE);
   desc->rngval  = (double *) calloc(desc->m, DSIZE);

   for (i = 0, ind = 0; i < desc->n; i++){
      desc->matbeg[i] = ind;
      /* indegree equals outdegree constraint */
      for (j = 0; j <= mpp->numnodes - 1; j++){
	 /* checks to see if node i is the start node of every edge arc */
	 if (mpp->head[i] == j){
	    desc->matind[ind] = j;
	    desc->matval[ind++] = 1;
	 }else if (mpp->tail[i] == j){
	    desc->matind[ind] = j;
	    desc->matval[ind++] = -1;
	 }
      }
      
      /* Now the constraint that each edge must be traversed at least once */
      if (i >= mpp->numarcs){ /* Check to see if it is an edge */
	 if (i < mpp->numarcs + mpp->numedges){
	    desc->matind[ind] = mpp->numnodes + i - mpp->numarcs;
	 }else{
	    desc->matind[ind] = mpp->numnodes + i - (mpp->numarcs +
						     mpp->numedges);
	 }
	 desc->matval[ind++] = 1;
	 /* desc->lb[i] = 0; */ /* Already set to zero from calloc */
	 desc->ub[i] = (double) (mpp->numarcs + mpp->numedges);
      }else{
	 desc->lb[i] = 1.0;
	 desc->ub[i] = (double) (mpp->numarcs + mpp->numedges);
      }
      desc->obj[i] = (double) (mpp->cost[i]);
   }
   desc->matbeg[i] = ind;
   
   /* set the initial right hand side */
   for (i = 0; i <= mpp->numnodes-1 ; i++){
      desc->rhs[i]   = 0;
      desc->sense[i] = 'E';
   }
   for (i = mpp->numnodes; i <= mpp->numnodes+mpp->numedges-1 ; i++){
      desc->rhs[i]   = 1;
      desc->sense[i] = 'G';
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
   return(DISP_NZ_INT);
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
   mpp_problem *mpp = (mpp_problem *)user;
   int * node_checker;
   int * cut_holder=NULL;
   int edge_direction=0;
   int i, j, nzcnt = 0;
   waiting_row **row_list = NULL;
   int *matind = NULL;
   cut_data *cut;
   int rhs_count;
   char *coef;
   double *matval = NULL;
   *new_row_num = cutnum;
   node_checker= (int *) calloc(mpp->numnodes, sizeof(int));
  
   if (cutnum > 0)
      *new_rows = row_list = (waiting_row **) calloc (cutnum,
						      sizeof(waiting_row *));

   for (j = 0; j < cutnum; j++){
      coef = (cut = cuts[j])->coef;
      cut_holder = (int *) cut->coef;

      cuts[j] = NULL;
      (row_list[j] = (waiting_row *) malloc(sizeof(waiting_row)))->cut = cut;
      switch (cut->type){
       case ODD_CUT:
	 matind = (int *) malloc(varnum * ISIZE);
	 nzcnt = 0;
	 rhs_count = 0;
	 for (i = 0; i <= mpp->numnodes - 1; i++){
	    node_checker[i]=0;
	 }
	 /*make array for 1 if node is in cut, 0 if not*/
	 for (i = 0; i < (cut->size/4); i++){
	    node_checker[cut_holder[i]]=1;
	 }
	 for (i = 0; i <= varnum - 1; i++){
	    if (vars[i]->userind >= mpp->numarcs+mpp->numedges){
	       if (node_checker[mpp->tail[(vars[i]->userind)-mpp->numedges]]==1&&
		   node_checker[mpp->head[(vars[i]->userind)-mpp->numedges]]==0){
		  matind[nzcnt] = i;
		  nzcnt++;
		  rhs_count++;
	       }
	    }else if ((vars[i]->userind >= mpp->numarcs)){
	       if ((node_checker[mpp->tail[vars[i]->userind]] == 1 &&
		    node_checker[mpp->head[vars[i]->userind]] == 0)){
		  matind[nzcnt] = i;
		  nzcnt++;
		  rhs_count++;
	       }
	    }else if (node_checker[mpp->tail[vars[i]->userind]] == 1 &&
		      node_checker[mpp->head[vars[i]->userind]] == 0){
	       matind[nzcnt] = i;
	       nzcnt++;
	       rhs_count++;
	    } else if (node_checker[mpp->tail[vars[i]->userind]] == 0 &&
		       node_checker[mpp->head[vars[i]->userind]] == 1){
	       rhs_count++;
	    }
	 }
	 
	 row_list[j]->matind = matind =
	    (int *) realloc((char *)matind, nzcnt*ISIZE);
	 cut->rhs=(rhs_count+1)/2;
	 row_list[j]->nzcnt = nzcnt;
	 row_list[j]->matval = matval = (double *) malloc(nzcnt * DSIZE);
	 for (i = nzcnt-1; i >= 0; i--)
	    matval[i] = 1;
	 cut->branch = ALLOWED_TO_BRANCH_ON;
	 break;	
	 
   	 
       default:
	 printf("Unrecognized cut type!\n");
      }
   }
   FREE(node_checker);
   
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
			    char *delete_rows)
{
   return(DEFAULT);
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

