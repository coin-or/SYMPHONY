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
#include <malloc.h>
#include <memory.h>

/* SYMPHONY include files */
#include "BB_constants.h"
#include "BB_macros.h"
#include "proccomm.h"
#include "lp_u.h"

/* SPP include files */
#include "spp.h"
#include "spp_common.h"

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
   spp_problem *spp;
   col_ordered *m;
   int colnum, info;

   spp = (spp_problem *) calloc(1, sizeof(spp_problem));
   *user = spp;

   spp->par = (spp_parameters *) calloc(1, sizeof(spp_parameters));
   
   receive_char_array((char *)spp->par, sizeof(spp_parameters));
   m = spp->cmatrix = (col_ordered *) calloc(1, sizeof(col_ordered));
   receive_int_array(&colnum, 1);
   colnum = m->active_colnum = m->colnum;
   receive_int_array(&m->rownum, 1);
   receive_int_array(&m->nzcnt, 1);
   m->colnames = (int *) malloc(colnum * ISIZE);
   m->col_deleted = (char *) calloc(colnum/BITSPERBYTE + 1, CSIZE); /*calloc!*/
   m->obj = (double *) malloc(colnum * DSIZE);
   m->matbeg = (int *) malloc((colnum + 1) * ISIZE);
   m->matind = (row_ind_type *) malloc(m->nzcnt * sizeof(row_ind_type));
   receive_int_array(m->colnames, colnum);
   receive_dbl_array(m->obj, colnum);
   receive_int_array(m->matbeg, (colnum + 1));
   receive_char_array((char *)m->matind, m->nzcnt * sizeof(row_ind_type));
   
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
   spp_problem *spp = (spp_problem *) user;
   col_ordered *cm = spp->cmatrix;
   int i;

   desc->nz = cm->nzcnt;
   *maxn = desc->n;   /* note that the number of columns cannot increase */
   *maxm = 2 * desc->m;
   *maxnz = desc->nz + ((*maxm) * (*maxn) / 100);

   desc->matbeg = (int *) malloc((desc->n + 1) * ISIZE);
   desc->matind = (int *) malloc(desc->nz * ISIZE);
   desc->matval = (double *) malloc(desc->nz * DSIZE);
   desc->obj    = (double *) malloc(desc->n * DSIZE);
   desc->lb     = (double *) calloc(desc->n, DSIZE);
   desc->ub     = (double *) malloc(desc->n * DSIZE);
   desc->rhs    = (double *) malloc(desc->m * DSIZE);
   desc->sense  = (char *) malloc(desc->m * CSIZE);
   desc->rngval = (double *) malloc(desc->m * DSIZE);

   memcpy((char *) desc->matbeg, (char *) cm->matbeg, (cm->colnum+1) * ISIZE);   
   memcpy((char *) desc->obj, (char *) cm->obj, cm->colnum * DSIZE);      

   for (i = cm->nzcnt - 1; i >= 0; i--) {
      desc->matind[i] = cm->matind[i];   /* cannot memcpy b/c int vs. short */
      desc->matval[i] = 1.0;
   }

   for (i = desc->n - 1; i >= 0; --i){
      desc->ub[i] = 1.0;
      /* desc->lb[i] = 0.0; */ /* Set by calloc */
   }

   for (i = desc->m - 1; i >= 0; i--) {
      desc->rhs[i] = 1.0;
      desc->sense[i] = 'E';
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
   return(TEST_ZERO_ONE);
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
   /* note that names of variables in the current solution will not be
      displayed by the default option */
   return(which_sol == DISP_RELAXED_SOLUTION ? DISP_NOTHING : DISP_NZ_INT);
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
   /* No cutting planes */

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
   /* Needed for cutting planes */
   return(where == LP_SOL_TO_CP ? SEND_NONZEROS : SEND_FRACTIONS);
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
   spp_problem *spp = (spp_problem *)(*user);

   FREE(spp->par);
   spp_free_cmatrix(spp->cmatrix);
   FREE(*user);

   return(USER_NO_PP);

}

/*===========================================================================*/

