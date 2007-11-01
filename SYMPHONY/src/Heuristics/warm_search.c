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
int warm_search (lp_prob *p, int * indices, double *values, int cnt, 
      double lpetol, double *heur_solution, double &new_obj_val, 
      int &is_feasible)
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
   tm_prob *tm = p->tm;
   tm_params tm_par = tm->par;
   double gap = (tm->has_ub)?fabs((tm->ub-tm->lb)/(fabs(tm->ub)+lpetol)):10;
   double elapsed_time, solve_stop_time, real_obj_value;
   sp_desc *sp = tm->sp;
   double start_time = wall_clock(NULL);
   double current_ub=tm->ub;
   double solve_start_time = wall_clock(NULL);
   int n, cnt2, ui, to_be_fixed, i;
   MIPdesc *mip2;
   sym_environment *env2;
   int *indices2;
   double *values2;
   const int verbosity=p->par.verbosity;

   if (p->bc_index%tm_par.warm_search_frequency != 1 || 
         p->node_iter_num != 1 || is_feasible) {
      /* use this heuristic only when bc_index is a multiple of
       * warm_search_frequency and we are doing first iteration on that node
       * and the solution is not IP feasible.
       */
      return IP_INFEASIBLE;
   }
  //printf("feasible\n");

   if (gap<tm_par.warm_search_min_gap) {
      /* Use this heuristic only when gap is large */
      return IP_INFEASIBLE;
   }
   //printf("min gap\n");

   elapsed_time = wall_clock(NULL) - tm->start_time;
   if (tm->stat.warm_search_time/elapsed_time>
         tm_par.warm_search_max_time_frac) {
      /* Use this heuristic only when a small amount of time has been spent on
       * this heuristic so far.
       */
      //printf("min time");
      return IP_INFEASIBLE;
   }

   PRINT(verbosity,1,("warm_search: entering warm_search\n"));

   /* see if we can use a reference solution from the solution pool */
   if (sp->num_solutions>0) {
      reference_solution = sp->solutions[sp->num_solutions-1];
      have_reference = TRUE;
      PRINT(verbosity,1,("warm_search: reference solution value = %f\n", 
               reference_solution->objval));
      PRINT(verbosity,1,("warm_search: reference solution length = %i\n", 
               reference_solution->xlength));
   }

   if (have_reference) {
      /* 
       * we have a reference solution and an LP solution. Fix those variables
       * which have the same value in both
       */
      if (tm_par.warm_search_enabled==2) {
         env2 = sym_open_environment();
         mip2 = create_copy_mip_desc (p->mip);
         sym_explicit_load_problem(env2, mip2->n, mip2->m, mip2->matbeg, 
               mip2->matind, mip2->matval, mip2->lb, mip2->ub, mip2->is_int, 
               mip2->obj, mip2->obj2, mip2->sense, mip2->rhs, mip2->rngval, 
               TRUE);
      } else {
         env2 = tm->warm_search_env;
         mip2 = env2->mip;
         sym_set_int_param(env2,"keep_warm_start", TRUE);
         sym_set_int_param(env2,"generate_cgl_cuts", FALSE);
         sym_set_int_param(env2,"use_reduced_cost_fixing", FALSE);
         sym_set_int_param(env2,"warm_search_enabled", FALSE);
         sym_set_int_param(env2,"warm_start_node_level_ratio", 
               tm_par.warm_search_node_level_ratio);
      }
      n = mip2->n;

      /* we have a solution to lp relaxation which is not IP feasible. We
       * explore its nbhd. */
      vars = current_lp->vars;

      /* if we want to use the original MIP while searching */
      num_ints=0;
      for (i=0;i<n;i++) {
         if (mip2->is_int[i]) {
            num_ints++;
         }
      }

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
      for (i=0;i<p->lp_data->n;i++) {
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
                        fixed_val[vars[i]->userind] == 
                        reference_solution->xval[var_index2]) {
                     may_be_fixed[vars[i]->userind] = TRUE; 
                     num_fixable++;
                  }
               }
            } else {
               /* fix zero variables also */
               fixed_val[vars[i]->userind]=0;
               if (!var_is_non_zero(vars[i]->userind, reference_solution->xind,
                        reference_solution->xlength, var_index2)) {
                  may_be_fixed[vars[i]->userind] = TRUE;
                  num_fixable++;
               }
            }
         }
      }

      PRINT(verbosity,1,("warm_search: warm_search_fix_fraction = %f\n",
               tm_par.warm_search_fix_fraction));

      if (num_fixable >= tm_par.warm_search_fix_fraction*num_ints) {
         to_be_fixed = (int) floor(tm_par.warm_search_fix_fraction*num_ints);
         num_fixed = 0;

         for (i=0;i<n;i++) {
            /* first reset all variable bounds */
            sym_set_col_lower(env2,i,p->mip->lb[i]);
            sym_set_col_upper(env2,i,p->mip->ub[i]);
         }
         /* fix only to_be_fixed no. of variables */
         if (tm_par.warm_search_enabled==1) {
            for (i=0;i<n;i++) {
               /* first reset all variable bounds */
               sym_set_col_lower(env2,i,p->mip->lb[i]);
               sym_set_col_upper(env2,i,p->mip->ub[i]);
            }
         }

         /* fix first to_be_fixed variables which can be fixed */
         for (i=0;i<n;i++) {
            ui = vars[i]->userind;
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

         PRINT(verbosity,1,("warm_search: fixed %d out of %d fixable and %d "
                  "total ints.\n",num_fixed,num_fixable,num_ints));

         /* update usage statistic */
         /* TODO */
         /* solve the small mip */
         tm->stat.warm_search_calls++;
         PRINT(verbosity,1,("warm search: current upper bound = %f\n",
                  current_ub));
         if (tm_par.warm_search_enabled==2) {
            sym_set_primal_bound(env2, current_ub);
         }
         sym_set_int_param(env2, "verbosity", -5);
         sym_set_dbl_param(env2, "time_limit", tm_par.warm_search_time_limit);
         solve_start_time = wall_clock(NULL);
         sym_warm_solve(env2);
         solve_stop_time = wall_clock(NULL);
         PRINT(verbosity,1,("warm_search: Time taken = %f\n",
                  solve_stop_time-solve_start_time));
         termstatus = sym_get_status(env2);
         if ((sym_get_obj_val(env2, &new_obj_val) ==
                  FUNCTION_TERMINATED_NORMALLY) && (new_obj_val < current_ub -
                     p->par.granularity))  {
            /* heur_solution is sent back to lp_wrapper.c in symphony */
            sym_get_col_solution(env2, heur_solution);
            termcode = IP_HEUR_FEASIBLE;
            PRINT(verbosity,1,("warm_search: found better soln: %12.8f, %f\n", 
                     new_obj_val,p->par.granularity));
            indices2 = p->lp_data->tmp.i1; /* n */
            values2 = p->lp_data->tmp.d; /* n */
            cnt2 = collect_nonzeros(p, heur_solution, indices2, values2);
            sp_add_solution(p,cnt2,indices2,values2,new_obj_val,p->bc_index);
            tm->stat.warm_search_successes++;

            if (p->mip->obj_sense == SYM_MAXIMIZE){
               real_obj_value=-new_obj_val+p->mip->obj_offset;
            } else {
               real_obj_value=new_obj_val+p->mip->obj_offset;
            }
            elapsed_time = wall_clock(NULL)-p->tm->start_time;
            PRINT(verbosity,-1,("warm search: found solution = %10.2f time = "
                     "%10.2f\n", real_obj_value,elapsed_time));

         } else if (termstatus==TM_TIME_LIMIT_EXCEEDED) {
            /*
             * this means that the time ran out before a solution could be
             * found.
            */
            PRINT(verbosity,1,("warm_search: Time limit reached.\n"));
            tm->stat.warm_search_tl_reached++;
            tm_par.warm_search_fix_fraction *= (1+
                  tm_par.warm_search_fix_frac_incr);
         } else if (termstatus==TM_NO_SOLUTION) {
            PRINT(verbosity,1,("warm_search: no soln found\n"));
            if (solve_stop_time-solve_start_time<0.5*
                  tm_par.warm_search_time_limit) {
               tm_par.warm_search_fix_fraction = 
                  tm_par.warm_search_fix_fraction*
                  (1.0-tm_par.warm_search_fix_frac_decr);
            } 
         } else {
            PRINT(verbosity,1,("warm_search: bad soln %f\n", new_obj_val));
         }
      } else {
         PRINT(verbosity,1,("warm_search: insufficient fixation (%d,%f) in "
                  "warm_search. Leaving without doing anything.\n",num_fixable,
                  tm_par.warm_search_fix_fraction*num_ints));
         tm_par.warm_search_fix_fraction = tm_par.warm_search_fix_fraction*
            (1.0-tm_par.warm_search_fix_frac_decr/10);
      }

      PRINT(verbosity,1,("warm_search: new warm_search_fix_fraction = %f\n",
               tm_par.warm_search_fix_fraction));

      /* free the data structures */
      FREE(fixed_val);
      FREE(may_be_fixed);
      if (tm_par.warm_search_enabled==2) {
         sym_close_environment(env2);
         free_mip_desc(mip2);
         FREE(mip2);
      }
   } else {
      PRINT(verbosity,1,("Dont know what to do, no reference solution "
               "exists\n"));
      termcode = IP_INFEASIBLE;
   }

   tm->stat.warm_search_time = tm->stat.warm_search_time+
      wall_clock(NULL)-start_time;
   PRINT(verbosity,1,("Leaving warm_search\n"));
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
   int j;

   for (j=0;j<varnum;j++) {
      if (i==indices[j]) {
         var_index = j;
         termcode = TRUE;
         break;
      }
   }
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
