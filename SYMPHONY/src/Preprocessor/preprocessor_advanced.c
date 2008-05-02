/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*===========================================================================*/
/* This file is a part of the symphony-processor                             */
/* This file was written by Ashutosh (asm4@lehigh.edu)                       */
/*===========================================================================*/

#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <list>

#include "preprocessor.h"
#include "preprocessor_constants.h"
#include "symphony_api.h"
#include "master.h" 
#include "BB_constants.h"

/*===========================================================================*/
/*===========================================================================*/
int prep_advanced(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		  impList *implistL, impList *implistU, prep_stats *stats,
		  prep_params *prep_par)
{
   int termcode, termstatus, termstatus2, iteration_num;
   int verbosity = prep_par->probe_verbosity;
   
   //   prep_display_mip(P);
   termcode = PREP_UNMODIFIED;
   termstatus = PREP_MODIFIED;
   iteration_num = 0;
   /* find logical implications of the form x_i = 1 => x_j = 0 etc */
   while (termstatus == PREP_MODIFIED) {
      if (verbosity>=1) {
	 printf ("Advanced prep Iteration number = %d\n", ++iteration_num);
      }
      termstatus = prep_find_logical_imp(P, row_P, lhs, implistL, implistU,
					 stats, prep_par);
      if (termstatus==PREP_MODIFIED) {
	 termstatus2 = prep_basic(P, row_P, lhs, stats, prep_par);
	 if (termstatus2==PREP_INFEAS) {
	    for (int i=0; i<P->n; i++) {
	       implistL[i].clear();
	       implistU[i].clear();
	    }
	    printf ("Advanced Preprocessing declared problem infeasible\n");
	    return PREP_INFEAS;
	    /* more to do here */
	 }
      }
      else if (termstatus==PREP_INFEAS) {
	 if (termstatus2==PREP_INFEAS) {
	    for (int i=0; i<P->n; i++) {
	       implistL[i].clear();
	       implistU[i].clear();
	    }
	    printf ("Advanced Preprocessing declared problem infeasible\n");
	    return PREP_INFEAS;
	    /* more to do here */
	 }
      }
   }
   //   prep_display_mip(P);
   if (verbosity >= 2) {
      printf ("implications of fixing at 0: \n");
      termstatus = prep_list_imp(P, implistL);
      printf ("implications of fixing at 1: \n");
      termstatus = prep_list_imp(P, implistU);
   }
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_find_logical_imp (MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			   impList *implistL, impList *implistU,
			   prep_stats *stats, prep_params *prep_par)
{
   int termcode = PREP_UNMODIFIED;
   int termstatus;
   int col_num;
   int verbosity = prep_par->probe_verbosity;
   
   /* copies of our data structures in which some variables are fixed */
   MIPdesc *fP;
   rowpackedarray *row_fP;
   lhsparams *flhs;
   prep_stats *fstats;
   bool reset_fP = TRUE;

   if (verbosity>=1) {
      printf ("Finding Logical implications\n");
   }
   /* initialize them */
   fP              = create_copy_mip_desc(P);

   row_fP          = (rowpackedarray *) malloc(sizeof(rowpackedarray));
   row_fP->matbeg  = (int *)    malloc(ISIZE * (P->m + 1));
   row_fP->matind  = (int *)    malloc(ISIZE*P->nz);
   row_fP->matval  = (double *) malloc(DSIZE*P->nz);
   row_fP->isdeleted  = (char *) malloc((CSIZE*P->m+1));
   
   flhs                  = (lhsparams *) malloc(sizeof(lhsparams));
   flhs->ub_is_infinite  = (int *)calloc(P->m, ISIZE);
   flhs->lb_is_infinite  = (int *)calloc(P->m, ISIZE);
   flhs->ub_inf_var      = (int *)calloc(P->m, ISIZE);
   flhs->lb_inf_var      = (int *)calloc(P->m, ISIZE);
   flhs->ubound          = (double *) malloc(P->m * DSIZE);
   flhs->lbound          = (double *) malloc(P->m * DSIZE);

   fstats                = (prep_stats *)malloc(sizeof(prep_stats));
   
   /* for each binary variable, create of copy of MIP */
   for (col_num=0; col_num < P->n; col_num++){
      if (prep_isBinary(P, col_num)) {
	 if (verbosity >= 3) {
	    printf ("Analyzing implications of fixing variable %d\n", col_num);
	 }
	 
	 /* first fix col_num at zero */
	 if (reset_fP) {
	    termstatus = prep_reset_mip_data(P, row_P, lhs, fP, row_fP, flhs);
	 }
	 reset_fP = TRUE;
	 termstatus = prep_apply_imp(fP, row_fP, flhs, col_num, 0, implistL,
				     fstats);
	 prep_update_const_bounds(fP, row_fP, flhs, 'U', col_num, 0);
	 fP->ub[col_num] = 0;
	 if (prep_par->probe_level==2) {
	    termstatus = prep_basic(fP, row_fP, flhs, fstats, prep_par);
	 }
	 else if (prep_par->probe_level==1) {
	    termstatus = prep_find_imp_sibling(fP, row_fP, flhs, col_num,
					       verbosity);
	    if (termstatus==PREP_MODIFIED) {
	       termstatus = prep_check_feas(fP, row_fP, flhs);
	    }
	 }
	 if (termstatus==PREP_INFEAS) {
	    /*
	      variable col_num may be fixed at 1. first update all the
	      lhs-parameters based on the new bound. Then, update the bound.
	    */
	    if (verbosity >= 2) {
	       printf("fixing variable %d to 1\n", col_num);
	    }
	    termstatus = prep_apply_imp(P, row_P, lhs, col_num, 1, implistU,
					stats);
	    prep_update_const_bounds(P, row_P, lhs, 'L', col_num, 1);
	    P->lb[col_num] = 1;
	    termcode = PREP_MODIFIED;
	    reset_fP = TRUE;
	    continue;
	 }
	 else if (termstatus == PREP_MODIFIED) {
	    termstatus = prep_add_imp(P, fP, col_num, implistL, verbosity);
	    reset_fP = TRUE;
	 }
	 else if (termstatus == PREP_UNMODIFIED && implistL[col_num].empty()) {
	    /*
	      undo changes in fP so that its ready for the next binary
	      variable, and we do not need to remove and create it again
	    */
	    prep_update_const_bounds(fP, row_fP, flhs, 'U', col_num, 1);
	    fP->ub[col_num] = 1;
	    reset_fP = FALSE;
	 }
	 else {
	    reset_fP = TRUE;
	 }

	 
	 /* now fix col_num at one */
	 if (reset_fP) {
	    termstatus = prep_reset_mip_data(P, row_P, lhs, fP, row_fP, flhs);
	 }
	 reset_fP = TRUE;
	 termstatus = prep_apply_imp(fP, row_fP, flhs, col_num, 1, implistU,
				     fstats);
	 prep_update_const_bounds(fP, row_fP, flhs, 'L', col_num, 1);
	 fP->lb[col_num] = 1;
	 if (prep_par->probe_level==2) {
	    termstatus = prep_basic(fP, row_fP, flhs, fstats, prep_par);
	 }
	 else if (prep_par->probe_level==1) {
	    termstatus = prep_find_imp_sibling(fP, row_fP, flhs, col_num,
					       verbosity);
	    if (termstatus==PREP_MODIFIED) {
	       termstatus = prep_check_feas(fP, row_fP, flhs);
	       printf ("breaking here*************************\n");
	    }
	 }
	 if (termstatus==PREP_INFEAS) {
	    /*
	      variable col_num may be fixed at 0. first update all the
	      lhs-parameters based on the new bound. Then, update the bound.
	    */
	    if (verbosity >= 2) {
	       printf("fixing variable %d to 0\n", col_num);
	    }
	    termstatus = prep_apply_imp(P, row_P, lhs, col_num, 0, implistL,
					stats);
	    prep_update_const_bounds(P, row_P, lhs, 'U', col_num, 0);
	    P->ub[col_num] = 0;
	    termcode = PREP_MODIFIED;
	    reset_fP = TRUE;
	    continue;
	 }
	 else if (termstatus == PREP_MODIFIED) {
	    termstatus = prep_add_imp(P, fP, col_num, implistU, verbosity);
	    reset_fP = TRUE;
	 }
	 else if (termstatus == PREP_UNMODIFIED && implistU[col_num].empty()) {
	    /*
	      undo changes in fP so that its ready for the next binary
	      variable, and we do not need to remove and create it again
	    */
	    prep_update_const_bounds(fP, row_fP, flhs, 'L', col_num, 0);
	    fP->lb[col_num] = 0;
	    reset_fP = FALSE;
	 }
	 else {
	    reset_fP = TRUE;
	 }
      }
   }

   if (fP!=0) {
      termstatus = prep_delete_structs(fP, row_fP, flhs, fstats);
      fP = 0;
      row_fP = 0;
      flhs = 0;
   }
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_add_imp(MIPdesc *P, MIPdesc *fP, const int col_num, impList *il,
		 int verbosity)
{
   int is_imp = 0;
   double imp_val;		/* if an implication already exists,
				   this is the value of the variable */
   
   for (int i=0; i<P->n; i++) {
      if (fP->lb[i]==fP->ub[i] && P->lb[i]<P->ub[i] && (i!=col_num)) {
	 is_imp = prep_isImplied(col_num, i, il, imp_val);
	 if (is_imp==0) {
	    /* safely add the new implication */
	    if (verbosity >= 3) {
	       printf ("fixing var %d at %f fixed var %d at %f.\n", col_num,
		       fP->ub[col_num], i, fP->ub[i]);
	    }
	    implication temp_imp(i, fP->ub[i]);
	    il[col_num].push_back(temp_imp);
	 }
	 else if (imp_val!=fP->ub[i]){
	    printf ("Error detected while preprocessing. ");
	    printf ("You should not believe in the output now\n.");
	 }
      }
   }
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
bool prep_isImplied(const int col_num, const int imp_col, impList *il,
		    double &imp_val)
{
   int termcode = 0;
   /*
     function to check if fixing col_num in the list il already has an
     implication for the variable imp_col. if yes, imp_val will then have the
     implied value stored in the list.

     we assume that there are no duplicate entries in a list.
     decent assumption, since multiple entries with different fixed_val implies
     that the col_num cannot be fixed as done in this list.
   */
   std::list<implication>::const_iterator iter;
   for (iter=il[col_num].begin(); iter != il[col_num].end(); iter++) {
      if (iter->fixed_col == imp_col ) {
	 imp_val = iter->fixed_val;
	 termcode = 1;
	 break;
      }
   }
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_reset_mip_data (MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			 MIPdesc *fP, rowpackedarray *row_fP, lhsparams *flhs)
{
   int termcode = 0;
   /*
     ASSUMPTION: P and fP have the same size.
     only some parts of mip are reset.
   */

   memcpy(fP->ub,     P->ub,     DSIZE * P->n); 
   memcpy(fP->lb,     P->lb,     DSIZE * P->n);    
   memcpy(fP->matbeg, P->matbeg, ISIZE * (P->n + 1));
   memcpy(fP->rhs, P->rhs,       DSIZE * P->m); 
   memcpy(fP->rngval, P->rngval, DSIZE * P->m);
   memcpy(fP->matval, P->matval, DSIZE * P->nz);  
   memcpy(fP->matind, P->matind, ISIZE * P->nz);
   
   memcpy(row_fP->matbeg, row_P->matbeg, ISIZE*(P->m + 1));
   memcpy(row_fP->matind, row_P->matind, ISIZE*P->nz);
   memcpy(row_fP->matval, row_P->matval, DSIZE*P->nz);
   memcpy(row_fP->isdeleted, row_P->isdeleted, CSIZE*(P->m+1));
   row_fP->rows_deleted = row_P->rows_deleted;

   memcpy(flhs->ub_is_infinite, lhs->ub_is_infinite, ISIZE*P->m);
   memcpy(flhs->lb_is_infinite, lhs->lb_is_infinite, ISIZE*P->m);
   memcpy(flhs->ub_inf_var, lhs->ub_inf_var, ISIZE*P->m);
   memcpy(flhs->lb_inf_var, lhs->lb_inf_var, ISIZE*P->m);
   memcpy(flhs->ubound, lhs->ubound, DSIZE*P->m);
   memcpy(flhs->lbound, lhs->lbound, DSIZE*P->m);
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_apply_imp (MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		    int col_num, double col_val, impList *il,
		    prep_stats *stats)
{
   /*
     search for any available logical implications of such a fix
     apply those implications as well.     
     this function does DOES NOT do the following:
     fix col_num at col_val
   */
   
   std::list<implication>::const_iterator iter;

   /* this is applicable only to binary variables */
   if (!(prep_isBinary(P, col_num))) {
      return 0;
   }

   for (iter=il[col_num].begin(); iter != il[col_num].end(); iter++) {
      if (iter->fixed_val==P->lb[iter->fixed_col]) {
	 /* the variable is fixed at lb */
	 prep_update_const_bounds(P, row_P, lhs, 'U', iter->fixed_col,
				  iter->fixed_val);
	 P->ub[iter->fixed_col] = iter->fixed_val;
	 stats->vars_fixed++;
      }
      else if (iter->fixed_val==P->ub[iter->fixed_col]) {
	 /* the variable is fixed at ub */
	 prep_update_const_bounds(P, row_P, lhs, 'L', iter->fixed_col,
				  iter->fixed_val);
	 P->lb[iter->fixed_col] = iter->fixed_val;
	 stats->vars_fixed++;
      }
      else {
	 /* the variable is fixed at some other value */
	 prep_update_const_bounds(P, row_P, lhs, 'L', iter->fixed_col,
				  iter->fixed_val);
	 P->lb[iter->fixed_col] = iter->fixed_val;

	 prep_update_const_bounds(P, row_P, lhs, 'U', iter->fixed_col,
				  iter->fixed_val);
	 P->ub[iter->fixed_col] = iter->fixed_val;
	 stats->vars_fixed++;
      }
   }
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_list_imp(MIPdesc *P, impList *il)
{
   int termcode = 0;
   /* function to list the logical implications found so far */
   std::list<implication>::const_iterator iter;
   for (int i=0; i<P->n; i++) {
      for (iter=il[i].begin(); iter != il[i].end(); iter++) {
	 printf (" variable %d => variable %d = %f\n", i, iter->fixed_col,
		 iter->fixed_val);
      }
   }
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
