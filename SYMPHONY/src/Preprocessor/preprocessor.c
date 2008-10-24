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
#include "sym_macros.h"
#include "symphony.h"

/*===========================================================================*/
/*===========================================================================*/

/* this will take control of the mip desc, it has to be copied before 
   not to lose the original mip*/
int preprocess_mip (sym_environment *env)
{
/*
   This function is the master of the preprocessing part. It calls and 
   controls other functions to perform preprocessing jobs.
*/
   int termcode;		/* return status of this function, 0 normal, 1
				   error */
   int termstatus;		/* return status of functions called herein */
   
   /* FIXME get a copy of mip, 
      in case of an error, ignore presolve and continue with that...*/

   MIPdesc *mip = env->mip;
   prep_params params = env->par.prep_par; 
   
   int verbosity = params.verbosity;
   int p_level = params.level;
   
   if(env->par.verbosity < 0){
      verbosity = -1;
   }
   
   if (p_level <= 0) {
      if(verbosity >= 0){
	 printf ("Skipping Preprocessor\n");
      }
      termcode = 0;
      return(termcode);
   }

   /* need to fill in the row ordered vars of mip */

   double start_time = wall_clock(NULL);
   //double mark_time; 

   /* Start with Basic Preprocessing */

   PRINT(verbosity, 0, ("Starting Preprocessing...\n")); 

   
   //mark_time = wall_clock(NULL);


   /* initialize prep desc */
   /* have a func for these */
   PREPdesc * P = (PREPdesc *)calloc(1, sizeof(PREPdesc)); 

   /* get a copy of orig mip */

   if(env->prep_mip){
      free_mip_desc(env->prep_mip);
      FREE(env->prep_mip);
   }
   
   P->orig_mip = env->orig_mip = create_copy_mip_desc(mip);
   P->mip = env->prep_mip = mip;
   P->params = params;
   
   P->stats.nz_coeff_changed = (char *)calloc(CSIZE ,mip->nz);
   
   /* these will be needed for both bassic and advanced prep functions
      so we call them here */
   termcode = prep_fill_row_ordered(P);
   
   termcode = prep_initialize_mipinfo(P);//mip, params, &(P->stats));   

   /* no changes so far on column based mip*/

   if(!prep_quit(termcode)){
      termcode = prep_basic(P);
   }
   //PRINT(verbosity, 1, ("Basic Prep time: %f...\n", 
   //		wall_clock(NULL) - mark_time));

   if(env->par.verbosity >= 0){
      prep_report(P, termcode);
   }
   
   if(params.reduce_mip){
      prep_restore_rootdesc(env);
   }
   

   /* debug */
   /*----------*/
   if(params.write_mps || params.write_lp){
      char file_name[80] = "";
      sprintf(file_name, "%s_prep", env->probname);
      
      if(params.write_mps){
	 sym_write_mps(env, file_name);
      }
      if(params.write_lp){
	 sym_write_lp(env, file_name);
      }
   }
   
   /*----------*/

   /* since we use the original mip desc */
   P->mip = 0;
   P->orig_mip = 0;
   
   prep_close_desc(P);   

   PRINT(verbosity, 0, ("Total Presolve Time: %f...\n\n", 
			wall_clock(NULL) - start_time));   

   //exit(0);
  /* Do advanced Preprocessing */

#if 0
   impList *implistL, *implistU;
   implistL = new impList[mip->n];
   implistU = new impList[mip->n]; /* arrays of implication-lists */

   if (prep_par.do_probe) {
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
   /*
     sym_environment *env2 = sym_open_environment();
     sym_explicit_load_problem(env2, P->n, P->m, P->matbeg, P->matind, P->matval, P->lb, P->ub, P->is_int, P->obj, P->obj2, P->sense, P->rhs, P->rngval, TRUE);
   */
   /* Comment out these lines to see sometimes erratic behaviour */
   /*
     int *indices = (int *)malloc(ISIZE);
     indices[0] = 2;
     termcode = sym_delete_rows(env2, 1, indices);
     printf("termcode from delete_rows = %d\n", termcode);
     prep_display_mip(env2->mip);
     sym_set_int_param(env2,"verbosity",0);
     sym_solve(env2);
     sym_close_environment(env2);
     exit(0);
   */
   /* end erratic behaviour */
   
   /* old environment */
   free_mip_desc(env->mip);
   free(env->base);
   free(env->rootdesc->uind.list);
   free(env->rootdesc);
   sym_explicit_load_problem(env, P->n, P->m, P->matbeg, P->matind, P->matval, P->lb, P->ub, P->is_int, P->obj, P->obj2, P->sense, P->rhs, P->rngval, TRUE);
   sym_set_col_names(env,P->colname); 
   prep_purge_del_rows2(env, P, row_P, rows_purged);
   /* end old environment */
   
   stats->rows_deleted = stats->rows_deleted + rows_purged;
   /* preprocessing ends */
   
   prep_free_row_structs(lhs, row_P);
   free(row_P);
   free(lhs);
   for (int i=0; i<P->n; i++) {
      implistL[i].clear();
      implistU[i].clear();
   }
   delete[] implistL;
   delete[] implistU;
   if (verbosity>=1||prep_par.display_stats==1) {
      prep_disp_stats(stats);
   }
   free_mip_desc(P);
   free(P);

   free(stats);
   if (verbosity>0) {
      printf("Leaving Preprocessor\n");
   }
#endif 

   return termcode; 
}


/*===========================================================================*/
/*===========================================================================*/
