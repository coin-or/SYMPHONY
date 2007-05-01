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

#include "proccomm.h"
#include "qsortucb.h"
#include "lp.h"
#include "messages.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "BB_types.h"
#include "pack_cut.h"
#include "solution_pool.h"
#include "warm_search.h"

#include "symphony_api.h"


/*===========================================================================*/
int warm_search(lp_prob *p, int * indices, double *values, int cnt, double lpetol, double *heur_solution, double &new_obj_val, int is_feasible)
{
   /*
     p: lp-problem data structure
     indices: indices of vars which are non-zero in the current lp-solution.
     values:  values of these variables.
     cnt:     no. of non-zero variables.
     lpetol:  tolerance.
     heur_sol: array of size n containing a heuristic solution if any. this is
     initialized in is_feasible_u()
     new_obj_val: objective value for the heuristic solution.
     is_feasible: 1 if the current lp_relaxation is ip feasible. in this case,
     we fix some (random) of the variables. we dont match them against the
     incumbent. 
   */

   sym_environment *env2;	/* pointer to symphony's new environment for
				   the sub-mip */
   MIPdesc *local_mip;		/* new mip which we will create */
   var_desc **vars;
   int count;
   sp_sol *reference_sol = NULL; /* solution from the pool which is used to fix
				    the variables*/
   int reference_sol_num = -1;	/* the index of ref_sol in the pool */
   LPdata *current_lp = p->lp_data;
   char *may_be_fixed;
   double *fixed_val;
   int var_index, var_index2;
   int num_fixed = 0;
   int num_fixable = 0;
   int do_local;

   int termcode = IP_INFEASIBLE, termstatus = 0;
   int num_nzfixed = 0, num_zfixed = 0;
   int num_ints=0;
   double time_limit;
   double time_left = (wall_clock(NULL)-p->tm->par.time_limit)/
                       p->tm->par.time_limit; 
   tm_prob *tm = p->tm;
   double gap = (tm->has_ub)?fabs((tm->ub-tm->lb)/tm->ub):10;
   int local_n = 0;		/* number of vars in the current lp relaxation
				 */
   sol_pool_desc *solpool = p->tm->solpool;

   double current_ub=sp_get_best_obj_value(p);

   double start_time = wall_clock(NULL);
   if (tm->stat.analyzed%p->par.rins_freq != 1 && !is_feasible) {
      return IP_INFEASIBLE;
   }

  /* we have a solution to lp relaxation which is not IP feasible. We explore
      its nbhd. */

   vars = current_lp->vars;
   OsiSolverInterface * model = current_lp->si;

   if (!is_feasible && p->par.rins_do_local) {
      do_local=1;
   }
   else {
      do_local=0;
   }

   /* either use the lp for the constraints and bounds */
   if (do_local) {
      /* allocate space for new mip */
      local_mip = (MIPdesc *) calloc(1, sizeof (MIPdesc));


      /* query OSI to get the info about the lp solved */
      
      const CoinPackedMatrix * matrix = model->getMatrixByCol();

      local_mip->n      = SCOOPS_getNumCols(model);
      local_mip->m      = SCOOPS_getNumRows(model);
      local_mip->nz     = SCOOPS_getNumElements(model);
      const int *vectorLengths   = SCOOPS_getColMatLen(matrix);
      const int *vectorStarts    = SCOOPS_getColMatBeg(matrix);
      const int *vectorIndices   = SCOOPS_getColMatInd(matrix);
      const double *vectorValues = SCOOPS_getColMatVal(matrix);

      /* allocate space for new MIP */
      local_mip->is_int = (char *)  calloc (local_mip->n, CSIZE); 
      local_mip->matbeg = (int *)   malloc ((local_mip->n+1)*ISIZE);
      local_mip->matind = (int *)   malloc (local_mip->nz*ISIZE);
      local_mip->matval = (double *)malloc (local_mip->nz*DSIZE);
      local_mip->obj    = (double *)malloc (local_mip->n*DSIZE);
      local_mip->rhs    = (double *)malloc (local_mip->m*DSIZE);
      local_mip->rngval = (double *)malloc (local_mip->m*DSIZE);
      local_mip->sense  = (char *)  malloc ((local_mip->m)*CSIZE);
      local_mip->lb     = (double *)malloc (local_mip->n*DSIZE);
      local_mip->ub     = (double *)malloc (local_mip->n*DSIZE);

      memcpy(local_mip->obj, SCOOPS_getObjCoefficients(model), 
             local_mip->n*DSIZE);
      memcpy(local_mip->rhs, model->getRightHandSide(), local_mip->m*DSIZE);
      memcpy(local_mip->rngval, model->getRowRange(), local_mip->m*DSIZE);
      memcpy(local_mip->sense, model->getRowSense(), (local_mip->m)*CSIZE);
      memcpy(local_mip->lb, SCOOPS_getColLB(model), local_mip->n*DSIZE);
      memcpy(local_mip->ub, SCOOPS_getColUB(model), local_mip->n*DSIZE);

      /* copy values into the new MIP */
      count = 0;
      for (int i=0;i<local_mip->n;i++) {
	 local_mip->matbeg[i] = count;
	 memcpy(&(local_mip->matind[local_mip->matbeg[i]]),
		&(vectorIndices[vectorStarts[i]]), vectorLengths[i]*ISIZE);
	 memcpy(&(local_mip->matval[local_mip->matbeg[i]]),
		&(vectorValues[vectorStarts[i]]), vectorLengths[i]*DSIZE);
	 count = count+vectorLengths[i];
      }
   
      local_mip->matbeg[local_mip->n] = local_mip->nz;

      /* set which variables are integers */
      for (int i=0;i<local_mip->n;i++) {
	 if (vars[i]->is_int) {
	    local_mip->is_int[i] = TRUE;
	    num_ints++;
	 }
	 else {
	    local_mip->is_int[i] = FALSE;
	 }
      }
      local_n = local_mip->n;
   }
   else {
      /* if we want to use the original MIP while searching */
      local_mip = p->mip;
      local_n   = SCOOPS_getNumCols(model);
      num_ints=0;
      for (int i=0;i<local_mip->n;i++) {
	 if (local_mip->is_int[i]) {
	    num_ints++;
	 }
      }
   }


   if (do_local) {
      env2 = sym_open_environment();
   } else {
      env2 = p->tm->rins_sym_env;
   }
      
   if (p->tm->par.rins_initialized==FALSE || do_local) {

      /* make sure the new mip-solver doesnt change our mip */
      char make_copy=TRUE;

      /* load problem and initiate a new instance of SYMPHONY*/
      termstatus = sym_explicit_load_problem(env2, local_mip->n, local_mip->m, local_mip->matbeg, local_mip->matind, local_mip->matval, local_mip->lb, local_mip->ub, local_mip->is_int, local_mip->obj, local_mip->obj2, local_mip->sense, local_mip->rhs, local_mip->rngval, make_copy);

      if (do_local==FALSE) {
	 p->tm->par.rins_initialized = TRUE;
      
	 cut_pool *cp = p->tm->cpp[0];
	 /*       printf("Total cuts = %d\n",cp->cut_num); */
	 /*       printf("level = %d\tindex = %d\n",p->bc_level,p->bc_index);
		  */ 
	 for (int cut_num=0;cut_num<cp->cut_num;cut_num++) {
	    cp_cut_data *cpcutdata = cp->cuts[cut_num];
	    if (cpcutdata->level > 0) {
	       continue;
	    }
	 
	    cut_data cut = cpcutdata->cut;

	    int numelems  = ((int *)cut.coef)[0]; /* no. of elements in the cut
						     */ 
	    char rowsen   = cut.sense;    /* sense of the row: L,G,R  */
	    double rowrhs = cut.rhs;      /* rhs value */
	    double rowrng = cut.range;    /* range */
	 
	    int *matind = (int*) malloc(numelems*ISIZE); /* indices of vars
								*/ 
	    double *matval = (double*) malloc(numelems*DSIZE); /* values */
	    memcpy (matind, cut.coef+ISIZE, numelems*ISIZE);
	    memcpy (matval, cut.coef+(numelems+1)*ISIZE, numelems*DSIZE);
	 
	    /* 	 printf("Cut Number = %d\n",cut_num); */
	    /* 	 int t1; */
	    /* 	 for (t1 = 0;t1<numelems;t1++) { */
	    /* 	    printf("%d\t%d\t%g\n",t1,matind[t1], matval[t1]); */
	    /* 	 } */

	    termstatus = sym_add_row(env2, numelems, matind, matval, rowsen,
				     rowrhs, rowrng);
	    if (termstatus != FUNCTION_TERMINATED_NORMALLY) {
	       printf("RINS: Error adding cuts to new environment. Chickening out\n");
	       exit(0);
	    }
	    FREE(matind);
	    FREE(matval);
	 }
      }
   } else {
      /* do_global and env2 is already initialized */
      /* in this case we simply reset the bounds */
      for (int i=0;i<p->mip->n;i++) {
	 sym_set_col_lower(env2, i, p->mip->lb[i]);
	 sym_set_col_upper(env2, i, p->mip->ub[i]);
      }
   }

   if (p->par.rins_use_sol_pool && !is_feasible) {

      /* first choose the best solution randomly from the solpool */
      double rand01 = CoinDrand48();
      sp_assign_probability(solpool, *local_mip->lb);
      if (p->par.verbosity>0) {
	 sp_display_sol_pool(solpool);
      } 
      for (int i=0;i<solpool->active_sols;i++) {
	 if (rand01 <= solpool->sols[i]->selection_cum_prob) {
	    reference_sol = solpool->sols[i];
	    reference_sol_num = i;
	    break;
	 }
      }
      /* in the event that above assignment didnt work */
      if (!reference_sol) {
	 reference_sol = solpool->sols[0];
	 reference_sol_num = 0;
      }
   } else if (is_feasible) {
	 reference_sol_num = -1;
	 reference_sol = NULL;
   } else {
      reference_sol_num = 0;
      reference_sol = solpool->sols[0];
   }
   
   PRINT(p->par.verbosity,1,("RINS: using solution %d of %d in the pool.\n",reference_sol_num+1, solpool->active_sols)); 
   /* fix variables */
   int c1, c2;
   c1=c2=0;
   fixed_val = (double *)malloc(local_n*DSIZE);
   may_be_fixed = (char *)calloc(local_n,CSIZE);
/*    print indices */
/*    for (int i=0;i<cnt;i++) { */
/*       printf("%d\t%d\n",i,indices[i]); */
/*    } */
/*    printf("\n\nuserindices:\n"); */
/*    for (int i=0;i<local_mip->n;i++) { */
/*       printf("%d\t%d\n",i,vars[i]->userind); */
/*    } */

   /* see which variables can be fixed */
   for (int i=0;i<local_n;i++) {
      int j = (do_local)?i:vars[i]->userind;
      if (vars[i]->is_int) {
	 /* fix non-zero var if its value is same as that in the incumbent */
	 if (rins_is_non_zero(vars[i]->userind, indices, cnt, var_index)) {
	    if (fabs(values[var_index]-floor(values[var_index]+0.5))<lpetol) {
	       fixed_val[i] = floor(values[var_index]+0.5);
	       if (!is_feasible && rins_is_non_zero(vars[i]->userind,
		   reference_sol->xind, reference_sol->xlength, var_index2)
		   && fixed_val[i] == reference_sol->xval[var_index2]) {
		  may_be_fixed[i] = TRUE;
		  num_fixable++;
	       }
	       else if (is_feasible) {
		  may_be_fixed[i] = TRUE;
		  num_fixable++;
	       }
	    }
	 }
	 /* fix zero variables also */
	 else {
	    fixed_val[i]=0;
	    if (is_feasible) {
	       may_be_fixed[i] = TRUE;
	       num_fixable++;
	    }
	    else if (!is_feasible && !rins_is_non_zero(vars[i]->userind, reference_sol->xind, reference_sol->xlength, var_index2)) {
	       may_be_fixed[i] = TRUE;
	       num_fixable++;
	    }
	 }
      }
   }

/*    printf("Num fixed = %d nonzero = %d zero = %d total ints = %d\n", */
/* 	  num_fixed, num_nzfixed, num_zfixed, num_ints); */
/*    printf ("fixation: (%d,%d) in RINS.\n",num_fixed,num_ints); */
   PRINT(p->par.verbosity,1,("Rins: rins_fix_fraction = %f\n",p->tm->par.rins_fix_fraction));
   if (num_fixable >= p->tm->par.rins_fix_fraction*num_ints) {
      int to_be_fixed = (int) floor(p->tm->par.rins_fix_fraction*num_ints);
      double rand01;
      num_fixed = 0;
      /* fix only to_be_fixed no. of variables */
      while (num_fixed<to_be_fixed) {
	 for (int i=0;i<local_n;i++) {
	    int j = (do_local)?i:vars[i]->userind;
	    if (vars[i]->is_int && may_be_fixed[i] && (double)num_fixed/to_be_fixed <= CoinDrand48()) {
	       sym_set_col_lower(env2, j, fixed_val[i]);
	       sym_set_col_upper(env2, j, fixed_val[i]);
	       num_fixed++;
	       if (num_fixed>=to_be_fixed) {
		  break;
	       }
	    }
	 }
      }

      PRINT(p->par.verbosity,-1,("Rins: fixed %d out of %d fixable and %d total ints.\n",num_fixed,num_fixable,num_ints));
      
      /* update usage statistic */
      sp_update_usage(solpool,reference_sol_num); 
      /* solve the small mip */
      tm->stat.rins_instances++;
      PRINT(p->par.verbosity,-1,("Rins: current upper bound = %f\n",current_ub));
      sym_set_primal_bound(env2, current_ub);
      sym_set_int_param(env2, "do_naive_rounding", 1);
      sym_set_int_param(env2, "rounding_freq", 10);
      sym_set_int_param(env2, "verbosity", -1);
      sym_set_dbl_param(env2, "time_limit", p->par.rins_time_limit);
      sym_set_dbl_param(env2, "rc_fix_factor",p->par.rins_rc_fix_factor);
      double solve_start_time = wall_clock(NULL);
      sym_solve(env2);
      double solve_stop_time = wall_clock(NULL);
      PRINT(p->par.verbosity,-1,("Rins: Time taken = %f\n",solve_stop_time-solve_start_time));
      termstatus = sym_get_status(env2);
      if ((sym_get_obj_val(env2, &new_obj_val) ==
	   FUNCTION_TERMINATED_NORMALLY) &&
	  (new_obj_val < current_ub - p->par.granularity))  {
	 /* heur_solution is sent back to lp_wrapper.c in symphony */
	 sym_get_col_solution(env2, heur_solution);
	 termcode = IP_HEUR_FEASIBLE;
	 PRINT(p->par.verbosity,-1,("Rins: found better soln: %12.8f, %f", new_obj_val,p->par.granularity));
	 int *indices2 = p->lp_data->tmp.i1; /* n */
	 double *values2 = p->lp_data->tmp.d; /* n */
	 int cnt = collect_nonzeros(p, heur_solution, indices2, values2);
	 sp_add_to_solpool(p,cnt,indices2,values2,new_obj_val,p->bc_index);
	 p->tm->stat.rins_successes++;
	 
	 /* not used */
	 /* if (p->par.heur_auto_level) { */
	 /*   p->par.rins_time_limit=p->par.rins_time_limit/1.05; */
	 /* } */
      } 
      else if (termstatus==TM_TIME_LIMIT_EXCEEDED) {
	 /*
	   this means that the time ran out before a solution could be found.
	 */
	 
	 /* not used */
	 /* if (p->par.heur_auto_level) { */
	 /*    p->par.rins_time_limit=p->par.rins_time_limit*1.05; */
	 /* } */

	 PRINT(p->par.verbosity,-1,("Rins: Time limit reached."));
	 p->tm->stat.rins_tl_reached++;
	 p->tm->par.rins_fix_fraction *= (1+p->par.rins_fix_frac_incr);
      }
      else if (termstatus==TM_NO_SOLUTION) {
	 /* not used */
	 /* if (p->par.heur_auto_level) { */
	 /*    p->par.rins_time_limit=p->par.rins_time_limit/1.05; */
	 /* } */
	 PRINT(p->par.verbosity,-1,("Rins: no soln found"));
	 if (solve_stop_time-solve_start_time<0.5*p->par.rins_time_limit) {
	    p->tm->par.rins_fix_fraction *= (1.0-p->par.rins_fix_frac_decr);
	 } 
      }
      else {
	 PRINT(p->par.verbosity,-1,("Rins: bad soln %f", new_obj_val));
      }
   }
   else {
      PRINT(p->par.verbosity,1,("Rins: insufficient fixation (%d,%f) in RINS. Leaving without doing anything.\n",num_fixable,p->tm->par.rins_fix_fraction*num_ints));
   }

   PRINT(p->par.verbosity,1,("Rins: new rins_fix_fraction = %f\n",p->tm->par.rins_fix_fraction));
   FREE(fixed_val);
   FREE(may_be_fixed);
   
   if (do_local) {
      free_mip_desc(local_mip);
      FREE(local_mip);
      sym_close_environment(env2);
   }

   PRINT(p->par.verbosity,1,("Leaving RINS\n"));

   p->tm->stat.rins_time = p->tm->stat.rins_time+wall_clock(NULL)-start_time;
   return termcode;
}

/*===========================================================================*/
bool rins_is_non_zero(int i, int *indices, int varnum, int &var_index)
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



