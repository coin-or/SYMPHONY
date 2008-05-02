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

#define ERR_TOL 1e-6		/* error tolerance, should really be a
				   parameter set based on the processor and
				   input method; unused so far */
/*===========================================================================*/
/*===========================================================================*/
#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sym_preprocessor.h"
#include "sym_prep_params.h"
#include "sym_master.h" 
#include "sym_constants.h" 

/*
#include "preprocessor_constants.h"
#include "symphony_api.h"
#include "proccomm.h" 
#include "master.h" 
#include "BB_constants.h" 
*/

/*===========================================================================*/
/*===========================================================================*/
int prep_basic(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
	       prep_stats *stats, prep_params *prep_par)
{
   /*
     This function is the master of the preprocessing part. It calls and 
     controls other functions to perform preprocessing jobs.
   */
   int termstatus;		/* return status of each function called
				   herein */
   int iterate_more;		/* =1 if we need more iterations of prep'ing */

   /* variables which collect stats*/
   int bounds_tightened=0;
   int vars_fixed=0;
   int bounds_integerized=0;
   int coeffs_changed=0;
   int verbosity = prep_par->prep_verbosity;
   
   /* let the show begin */
   
   /* Display the whole mip - for debugging */
#ifdef DISPLAY_MIP
   termstatus = prep_display_mip(P);
#endif   

   int iteration_number = 0;
   iterate_more = 1;
   /* main preprocessing loop */

   while (iterate_more) {
      iterate_more = 0;
      iteration_number++;
      if (verbosity>=1) {
	 printf("Basic iteration number: %d\n", iteration_number);
      }

      if (verbosity>=1) {
	 printf("Tightening bounds\n");
      }
      vars_fixed = 0;
      //termstatus = 0;
      termstatus = prep_tighten_bounds(P, row_P, lhs, -1, bounds_tightened, 
          vars_fixed, verbosity);
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 stats->bounds_tightened = stats->bounds_tightened+bounds_tightened;
	 stats->vars_fixed = stats->vars_fixed+vars_fixed;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (verbosity>=1) {
	 printf("Checking redundancy\n");
      }
      //termstatus = 0;
      termstatus = prep_check_redundancy(P, row_P, lhs, verbosity); 
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (verbosity>=1) {
	 printf("Fixing variables\n");
      }      
      vars_fixed = 0;
      termstatus = 0;
      termstatus = prep_fix_variables(P, row_P, lhs, vars_fixed, verbosity);
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 stats->vars_fixed = stats->vars_fixed+vars_fixed;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (verbosity>=1) {
	 printf("Improving Binary Coefficients\n");
      } 
      coeffs_changed = 0;
      termstatus = 0;
      termstatus = prep_BinImproveCoeffs(P, row_P, lhs, coeffs_changed,
				 verbosity);
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 stats->coeffs_changed = stats->coeffs_changed+coeffs_changed;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (P->m == 0) {
	 termstatus = prep_declare_solved(P, verbosity);
      }
   }

#ifdef DISPLAY_MIP
   termstatus = prep_display_mip(current_mip);
#endif
   
   if (iteration_number>1) {
      return PREP_MODIFIED;
   }
   else {
      return PREP_UNMODIFIED;
   }
   /* exit basic preprocessor */
}


/*===========================================================================*/
/*===========================================================================*/
