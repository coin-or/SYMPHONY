/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "sym_proccomm.h"
#include "qsortucb.h"
#include "sym_lp.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_pack_cut.h"
#include "sym_rounding.h"
#include "symphony.h"
#include "sym_warm_search.h"
#include "sym_sp.h"
#include "sym_lp_solver.h"
/*===========================================================================*/
int warm_search (lp_prob *p, int * indices, double *values, int cnt, double lpetol, double *heur_solution, double &new_obj_val, int &is_feasible)
{
   
   var_desc **vars;

   /* 
    * if there is a known feasible solution this will be used to guide the
    * warm search
    */
   sp_solution *reference_solution = NULL;
   int have_reference = FALSE;
   
   /* data from current lp solution */
   LPdata *current_lp = p->lp_data;
   char *may_be_fixed;
   double *fixed_val;
   int var_index, var_index2;
   int num_fixed = 0;
   int num_fixable = 0;

   int termcode = IP_INFEASIBLE, termstatus = 0;
   int num_ints=0;
   
   /*
    * double time_left = (wall_clock(NULL)-p->tm->par.time_limit)/
    * p->tm->par.time_limit; 
    */
   tm_prob *tm = p->tm;
   double gap = (tm->has_ub)?fabs((tm->ub-tm->lb)/tm->ub):10;
   sp_desc *sp = p->tm->sp;


   double start_time = wall_clock(NULL);
   if (tm->stat.analyzed%p->tm->par.warm_search_frequency != 1 || is_feasible) {
      return IP_INFEASIBLE;
   }

   PRINT(p->par.verbosity,15,("warm_search: entering warm_search\n"));
   double current_ub=p->tm->ub;
   int n;
   MIPdesc *mip2;
   sym_environment *env2;
   
   if (p->tm->par.warm_search_enabled==2) {
      env2 = sym_open_environment();
      mip2 = create_copy_mip_desc (p->mip);
      sym_explicit_load_problem(env2, mip2->n, mip2->m, mip2->matbeg, 
	    mip2->matind, mip2->matval, mip2->lb, mip2->ub, mip2->is_int, 
	    mip2->obj, mip2->obj2, mip2->sense, mip2->rhs, mip2->rngval, 
	    TRUE);
   } else {
      env2 = p->tm->warm_search_env;
      mip2 = env2->mip;
      sym_set_int_param(env2,"keep_warm_start", TRUE);
      sym_set_int_param(env2,"generate_cgl_cuts", FALSE);
      sym_set_int_param(env2,"warm_search_enabled", FALSE);
   }
   n = mip2->n;
 
   /* we have a solution to lp relaxation which is not IP feasible. We explore
      its nbhd. */
   vars = current_lp->vars;

   /* if we want to use the original MIP while searching */
   num_ints=0;
   for (int i=0;i<n;i++) {
      if (mip2->is_int[i]) {
	 num_ints++;
      }
   }

   /* see if we can use a reference solution from the solution pool */
   if (sp->num_solutions>0) {
      reference_solution = sp->solutions[sp->num_solutions-1];
      have_reference = TRUE;
      PRINT(p->par.verbosity,-1,("warm_search: reference solution value = %f\n", reference_solution->objval));
      PRINT(p->par.verbosity,-1,("warm_search: reference solution length = %i\n", reference_solution->xlength));
      
   }

   if (have_reference) {
      /* we have a reference solution and an LP solution. Fix those variables
       * which have the same value in both
       */
      int c1, c2;
      c1=c2=0;
      fixed_val = (double *)malloc(n*DSIZE);
      may_be_fixed = (char *)calloc(n,CSIZE);
      /*    print indices */
      /*    for (int i=0;i<cnt;i++) { */
      /*       printf("%d\t%d\n",i,indices[i]); */
      /*    } */
      /*    printf("\n\nuserindices:\n"); */
      /*    for (int i=0;i<local_mip->n;i++) { */
      /*       printf("%d\t%d\n",i,vars[i]->userind); */
      /*    } */

      /* see which variables can be fixed */
      for (int i=0;i<p->lp_data->n;i++) {
	 if (vars[i]->is_int) {
	    /* fix the variable if it has the same value in the reference
	     * solution
	     */
	    if (var_is_non_zero(vars[i]->userind, indices, cnt, var_index)) {
	       if (fabs(values[var_index]-floor(values[var_index]+0.5))<lpetol) {
		  fixed_val[vars[i]->userind] = floor(values[var_index]+0.5);
		  if (var_is_non_zero(vars[i]->userind,
			   reference_solution->xind,
			   reference_solution->xlength, var_index2) &&
			fixed_val[vars[i]->userind] == reference_solution->xval[var_index2]) {
		     may_be_fixed[vars[i]->userind] = TRUE; 
		     num_fixable++;
		  }
	       }
	    } else {
	       /* fix zero variables also */
	       fixed_val[vars[i]->userind]=0;
	       if (!is_feasible && !var_is_non_zero(vars[i]->userind, reference_solution->xind, reference_solution->xlength, var_index2)) {
		  may_be_fixed[vars[i]->userind] = TRUE;
		  num_fixable++;
	       }
	    }
	 }
      }
   } else {
      PRINT(p->par.verbosity,10,("Dont know what to do, no reference solution exists\n"));
      return(0);
   }

   PRINT(p->par.verbosity,-1,("warm_search: warm_search_fix_fraction = %f\n",p->tm->par.warm_search_fix_fraction));

   if (num_fixable >= p->tm->par.warm_search_fix_fraction*num_ints) {
      int to_be_fixed = (int) floor(p->tm->par.warm_search_fix_fraction*num_ints);
      num_fixed = 0;

      for (int i=0;i<n;i++) {
	 /* first reset all variable bounds */
	 sym_set_col_lower(env2,i,p->mip->lb[i]);
	 sym_set_col_upper(env2,i,p->mip->ub[i]);
      }
      /* fix only to_be_fixed no. of variables */
      if (p->tm->par.warm_search_enabled==1) {
	 for (int i=0;i<n;i++) {
	 /* first reset all variable bounds */
	 sym_set_col_lower(env2,i,p->mip->lb[i]);
	 sym_set_col_upper(env2,i,p->mip->ub[i]);
	 }
      }
      
      /* fix first to_be_fixed variables which can be fixed */
      for (int i=0;i<n;i++) {
	 int ui = vars[i]->userind;
	 if (vars[i]->is_int && may_be_fixed[ui]) {
	    sym_set_col_lower(env2, ui, fixed_val[ui]);
	    sym_set_col_upper(env2, ui, fixed_val[ui]);
	    num_fixed++;
	    if (num_fixed>=to_be_fixed) {
	       break;
	    }
	 }
      }
      /* TODO: randomize or select best candidates rather than first
       * to_be_fixed variables for fixing
       */

      PRINT(p->par.verbosity,-1,("warm_search: fixed %d out of %d fixable and %d total ints.\n",num_fixed,num_fixable,num_ints));
      
      /* update usage statistic */
      /* TODO */
      /* solve the small mip */
      tm->stat.warm_search_calls++;
      PRINT(p->par.verbosity,-1,("warm search: current upper bound = %f\n",current_ub));
      //sym_set_primal_bound(env2, current_ub);
      sym_set_int_param(env2, "verbosity", -1);
      sym_set_dbl_param(env2, "time_limit", p->tm->par.warm_search_time_limit);
      double solve_start_time = wall_clock(NULL);
      sym_warm_solve(env2);
      double solve_stop_time = wall_clock(NULL);
      PRINT(p->par.verbosity,-1,("warm_search: Time taken = %f\n",solve_stop_time-solve_start_time));
      termstatus = sym_get_status(env2);
      if ((sym_get_obj_val(env2, &new_obj_val) == FUNCTION_TERMINATED_NORMALLY) 
	    && (new_obj_val < current_ub - p->par.granularity))  {
	 /* heur_solution is sent back to lp_wrapper.c in symphony */
	 sym_get_col_solution(env2, heur_solution);
	 termcode = IP_HEUR_FEASIBLE;
	 PRINT(p->par.verbosity,-1,("warm_search: found better soln: %12.8f, %f\n", new_obj_val,p->par.granularity));
	 int *indices2 = p->lp_data->tmp.i1; /* n */
	 double *values2 = p->lp_data->tmp.d; /* n */
	 int cnt = collect_nonzeros(p, heur_solution, indices2, values2);
	 sp_add_solution(p,cnt,indices2,values2,new_obj_val,p->bc_index);
	 p->tm->stat.warm_search_successes++;
      } else if (termstatus==TM_TIME_LIMIT_EXCEEDED) {
	 /*
	   this means that the time ran out before a solution could be found.
	 */
	 PRINT(p->par.verbosity,-1,("warm_search: Time limit reached.\n"));
	 p->tm->stat.warm_search_tl_reached++;
	 p->tm->par.warm_search_fix_fraction *= (1+p->tm->par.warm_search_fix_frac_incr);
      } else if (termstatus==TM_NO_SOLUTION) {
	 PRINT(p->par.verbosity,-1,("warm_search: no soln found"));
	 if (solve_stop_time-solve_start_time<0.5*p->tm->par.warm_search_time_limit) {
	    p->tm->par.warm_search_fix_fraction = p->tm->par.warm_search_fix_fraction*(1.0-p->tm->par.warm_search_fix_frac_decr);
	 } 
      } else {
	 PRINT(p->par.verbosity,-1,("warm_search: bad soln %f\n", new_obj_val));
      }
   } else {
      PRINT(p->par.verbosity,-1,("warm_search: insufficient fixation (%d,%f) in warm_search. Leaving without doing anything.\n",num_fixable,p->tm->par.warm_search_fix_fraction*num_ints));
      p->tm->par.warm_search_fix_fraction = p->tm->par.warm_search_fix_fraction*(1.0-p->tm->par.warm_search_fix_frac_decr/10);
   }

   PRINT(p->par.verbosity,-1,("warm_search: new warm_search_fix_fraction = %f\n",p->tm->par.warm_search_fix_fraction));
   
   /* free the data structures */
   FREE(fixed_val);
   FREE(may_be_fixed);
   if (p->tm->par.warm_search_enabled==2) {
      sym_close_environment(env2);
      free_mip_desc(mip2);
   }
   
   p->tm->stat.warm_search_time = p->tm->stat.warm_search_time+wall_clock(NULL)-start_time;
   PRINT(p->par.verbosity,15,("Leaving warm_search\n"));
   return termcode;
}

/*===========================================================================*/
bool var_is_non_zero(int i, int *indices, int varnum, int &var_index)
{
  /* in the given array: indices, find if we have i somewhere */
  /* DO NOT assume indices is sorted */
  /* Return true if i is in indices else false */
  bool termcode = FALSE;
  var_index = -1;

  for (int j=0;j<varnum;j++) {
    if (i==indices[j]) {
      var_index = j;
      termcode = TRUE;
      break;
    }
  }
  return termcode;
}



