/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2005-2019 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <cstdio>

#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"
#include <iostream>
int main(int argc, char **argv)
{

  OsiSymSolverInterface si;
  si.parseCommandLine(argc, argv);
  si.loadProblem();

  si.setSymParam(OsiSymSensitivityAnalysis, true);

  si.initialSolve();

  int ind[2];
  double val[2];
  ind[0] = 4; val[0] = 7000;
  ind[1] = 7; val[1] = 6000;
  
  double lb = si.getLbForNewRhs(2, ind, val);
  double ub =  si.getUbForNewRhs(2, ind, val);

  printf("\nBounds for the new rhs:\n lb: %f\n ub: %f \n\n", lb, ub);

  return(0);

}

#else

#include "symphony.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
  
int main(int argc, char **argv)
{    
     
   sym_environment *env = sym_open_environment();   
   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);

   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "sensitivity_bounds", TRUE);
   sym_set_int_param(env, "sensitivity_rhs", FALSE);
   sym_set_int_param(env, "do_primal_heuristic", FALSE);
   sym_set_int_param(env, "should_use_rel_br", FALSE);
   sym_set_int_param(env, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env, "use_hot_starts", FALSE);
   sym_set_int_param(env, "should_warmstart_node", FALSE);
   sym_set_int_param(env, "prep_level", 0);
   sym_set_int_param(env, "tighten_root_bounds", FALSE);
   sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
   sym_set_int_param(env, "generate_cgl_knapsack_cuts", GENERATE_ONLY_IN_ROOT);
   sym_set_int_param(env, "generate_cgl_probing_cuts", GENERATE_ONLY_IN_ROOT);
   sym_set_int_param(env, "generate_cgl_clique_cuts", GENERATE_ONLY_IN_ROOT);
   sym_set_int_param(env, "generate_cgl_twomir_cuts", GENERATE_ONLY_IN_ROOT);
   sym_set_int_param(env, "generate_cgl_flowcover_cuts", GENERATE_ONLY_IN_ROOT);
   
   sym_solve(env);

   int ind[2];
   double val[2];
   ind[0] = 4; val[0] = 0;
   ind[1] = 7; val[1] = 0;
   
   int num_cols = 0;
   sym_get_num_cols(env, &num_cols);
   double *lb_lbs = (double *) malloc(DSIZE*num_cols);
   double *ub_lbs = (double *) malloc(DSIZE*num_cols);
   double *objval = (double *) calloc(num_cols, DSIZE);
   double *objval2 = (double *) calloc(num_cols, DSIZE);
   double *objval3 = (double *) calloc(num_cols, DSIZE);
   double *objval4 = (double *) calloc(num_cols, DSIZE);
   double *new_lb = (double *) malloc(DSIZE*num_cols);
   double *new_ub = (double *) malloc(DSIZE*num_cols);
   int *new_lb_inds = (int *) malloc(ISIZE*num_cols);
   int *new_ub_inds = (int *) malloc(ISIZE*num_cols);
   int num_fixed = 3;
   for (int i = 0; i < num_cols; i++){
      new_ub[i] = 0;
      new_ub_inds[i] = i;
   }
   for (int i = 0; i < num_cols; i++){
      new_lb[i] = 0;
      new_lb_inds[i] = i;
   }
   for (int i = num_fixed; i < num_cols; i++){
      sym_get_lb_for_new_rhs(env, 0, NULL, NULL,
			     num_fixed, new_lb_inds+i-num_fixed,
			     new_lb+i-num_fixed,
			     0, NULL, NULL,
			     lb_lbs+i);
      sym_get_lb_for_new_rhs(env, 0, NULL, NULL,
			     0, NULL, NULL,
			     num_fixed, new_ub_inds+i-num_fixed,
			     new_ub+i-num_fixed,
			     ub_lbs+i);
   }

   //sym_set_int_param(env, "verbosity", -2);

   for (int i = 0; i < num_fixed - 1; i++){
      sym_set_col_upper(env, i, 0);
   }
   for (int i = 3; i < num_cols; i++){
      sym_set_col_upper(env, i, 0);

      sym_warm_solve(env);
      sym_get_obj_val(env, objval+i);
      //sym_solve(env);
      sym_get_obj_val(env, objval2+i);
      assert(fabs(objval2[i] - objval[i] <= .1)); 

      sym_set_col_upper(env, i-3, 1);
   }
   for (int i = num_fixed; i < num_cols; i++){
      printf("LB: %.3f UB: %.3f %.3f\n", i, lb_lbs[i], objval[i],
	     objval2[i]);
   }

   for (int i = 0; i < num_fixed; i++){
      sym_set_col_upper(env, num_cols-i-1, 1);
   }

   printf("\n\n");

   for (int i = 0; i < num_fixed - 1; i++){
      sym_set_col_lower(env, i, 1);
   }
   for (int i = 3; i < num_cols; i++){
      sym_set_col_lower(env, i, 1);

      sym_warm_solve(env);
      sym_get_obj_val(env, objval+i);
      //sym_solve(env);
      sym_get_obj_val(env, objval2+i);
      assert(fabs(objval2[i] - objval[i] <= .1)); 

      sym_set_col_lower(env, i-3, 0);
   }
   for (int i = num_fixed; i < num_cols; i++){
      printf("LB: %.3f UB: %.3f %.3f\n", i, lb_lbs[i], objval[i],
	     objval2[i]);
   }

   sym_close_environment(env);
  
   return(0);
}  

#endif

