/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/*===========================================================================*/
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sym_preprocessor.h"
#include "sym_prep_params.h"
#include "sym_master_params.h"
#include "sym_master.h" 
#include "sym_constants.h" 
#include "symphony_api.h"
#include "symphony.h"

/*===========================================================================*/
/*===========================================================================*/
int sym_preprocess (sym_environment *env)
{
/*
   This function is the master of the preprocessing part. It calls and 
   controls other functions to perform preprocessing jobs.
*/
   int termcode;		/* return status of this function, 0 normal, 1
				   error */
   int termstatus;		/* return status of functions called herein */
   prep_params *prep_par = &env->par.prep_par;
   int verbosity = prep_par->prep_verbosity;

   if (prep_par->do_prep==0) {
      printf ("Skipping Preprocessor\n");
      termcode = 0;
      return(termcode);
   }
   MIPdesc *P;
   P = sym_create_copy_mip_desc (env);
   rowpackedarray *row_P;	/* rowmajor representation of constr. matrix */
   lhsparams *lhs;		/* bounds on constraints based on var bounds */
   impList *implistL, *implistU;
   prep_stats *stats;
   int rows_purged = 0;
   int bounds_integerized = 0;
   
   implistL = new impList[P->n];
   implistU = new impList[P->n]; /* arrays of implication-lists */
   stats    = (prep_stats *)malloc(sizeof(prep_stats));
   stats->rows_deleted = 0;
   stats->vars_fixed = 0;
   stats->coeffs_changed = 0;
   stats->bounds_tightened = 0;

   prep_display_mip(env->mip);
   /* Integerize variable bounds for Int Variables */
   termstatus = prep_integerize_bounds(P, bounds_integerized, verbosity);
   //stats->bounds_tightened = bounds_integerized;

   /* create a row-major representation of the MIP */
   row_P = (rowpackedarray *)malloc(sizeof(rowpackedarray));
   termstatus = prep_create_row_rep(P, row_P);

   /* find lhs parameters */
   lhs = (lhsparams *)malloc(sizeof(lhsparams));
   termstatus = prep_find_lhs_params(P, row_P, lhs);


   /* Start with Basic Preprocessing */
   if (prep_par->do_prep) {
      if (verbosity>=1) {
	 printf("Starting Basic Preprocessing.\n");
      }
      termstatus = prep_basic(P, row_P, lhs, stats, prep_par);
      if (verbosity>=1) {
	 printf("End Basic Preprocessing.\n");
      }
      if (termstatus == PREP_SOLVED) {
	 prep_free_row_structs(lhs, row_P);
	 free(row_P);
	 free(lhs);
	 free(P);
	 if (verbosity>=1) {
	    printf("Basic Preprocessing solved the problem.\n");
	 }
	 return PREP_SOLVED;
      }
   /* Do advanced Preprocessing */
      if (prep_par->do_probe) {
	 if (verbosity>=1) {
	    printf("Starting Advanced Preprocessing ...\n");
	 }
	 //termstatus = prep_advanced(P, row_P, lhs, implistL, implistU, stats,
	//			    prep_par);
	 if (verbosity>=1) {
	    printf("End Advanced Preprocessing.\n");
	 }
	 if (termstatus == PREP_SOLVED) {
	    prep_free_row_structs(lhs, row_P);
	    free(row_P);
	    free(lhs);
	    free(P);
	    if (verbosity>=1) {
	       printf("Advanced Preprocessing solved the problem.\n");
	    }
	    return PREP_SOLVED;
	 }
      }/* advanced preprocessing */
      
      /* new environment */
      sym_environment *env2 = sym_open_environment();
      sym_explicit_load_problem(env2, P->n, P->m, P->matbeg, P->matind, P->matval, P->lb, P->ub, P->is_int, P->obj, P->obj2, P->sense, P->rhs, P->rngval, TRUE);
      prep_display_mip(env2->mip);
      
      /* Comment out these lines to see magic */
      int *indices = (int *)malloc(ISIZE);
      indices[0] = 2;
      termcode = sym_delete_rows(env2, 1, indices);
      printf("tercode from delete_rows = %d\n", termcode);
      prep_display_mip(env2->mip);
      sym_set_int_param(env2,"verbosity",0);
      sym_solve(env2);
      sym_close_environment(env2);
      exit(0);
      /* end new environment */
      
      /* old environment */
      /*
      free_mip_desc(env->mip);
      free(env->base);
      free(env->rootdesc->uind.list);
      free(env->rootdesc);
      sym_explicit_load_problem(env2, P->n, P->m, P->matbeg, P->matind, P->matval, P->lb, P->ub, P->is_int, P->obj, P->obj2, P->sense, P->rhs, P->rngval, TRUE);
      prep_purge_del_rows2(env2, P, row_P, rows_purged);
      */
      /* end old environment */
      
      stats->rows_deleted = stats->rows_deleted + rows_purged;
   } /* preprocessing ends */

   prep_free_row_structs(lhs, row_P);
   free(row_P);
   free(lhs);
   for (int i=0; i<P->n; i++) {
      implistL[i].clear();
      implistU[i].clear();
   }
   delete[] implistL;
   delete[] implistU;
   if (verbosity>=1||prep_par->display_stats==1) {
      prep_disp_stats(stats);
   }
   free_mip_desc(P);
   free(P);

   free(stats);

   printf("Leaving Preprocessor\n");
   return 0; 
}


/*===========================================================================*/
/*===========================================================================*/
