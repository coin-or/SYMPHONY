<<<<<<< .working
/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _LP_PARAMS_H
#define _LP_PARAMS_H

#include "sym_constants.h"
#include "sym_timemeas.h"

/*---------------------------------------------------------------------------*\
| The list of parameters associated with processing a node in the branch and  |
| cut tree. See the README file for an explanation of the parameters          |
\*---------------------------------------------------------------------------*/

typedef struct CUT_TIME_OUT{
   double            first_cut_time_out;
   double            all_cuts_time_out;
}cut_time_out;

typedef struct CGL_PARAMS{
   /* Cut generation in LP */
   int               generate_cgl_cuts;
   int               generate_cgl_gomory_cuts;
   int               generate_cgl_knapsack_cuts;
   int               generate_cgl_oddhole_cuts;
   int               generate_cgl_probing_cuts;
   int               generate_cgl_mir_cuts;
   int               generate_cgl_clique_cuts;
   int               generate_cgl_flow_and_cover_cuts;
   int               generate_cgl_rounding_cuts;
   int               generate_cgl_lift_and_project_cuts;

   int               gomory_generated_in_root;
   int               knapsack_generated_in_root;
   int               oddhole_generated_in_root;
   int               probing_generated_in_root;
   int               mir_generated_in_root;
   int               clique_generated_in_root;
   int               flow_and_cover_generated_in_root;
   int               rounding_generated_in_root;
   int               lift_and_project_generated_in_root;
}cgl_params;

typedef struct LP_PARAMS{
   int               verbosity;
   double            granularity;
   int               use_cg;
   int               set_obj_upper_lim;
   int               do_primal_heuristic;
   double            time_limit;

   /* these two are passed directly to the lp solver */
   int               scaling;
   int               fastmip;

   int               try_to_recover_from_error;
   /* ZERO_ONE_PROBLEM / INTEGER_PROBLEM / MIXED_INTEGER_PROBLEM */
   int               problem_type;
   int               keep_description_of_pruned;

   int               not_fixed_storage_size;

   int               cut_pool_check_freq;

   int               load_balance_level;
   int               load_balance_iterations;
   int               load_balance_compare_candidates;

   double            fractional_diving_ratio;
   int               fractional_diving_num;

   /* parameters constraining the growth of the matrix */
   double            max_non_dual_feas_to_add_frac;
   int               max_non_dual_feas_to_add_min;
   int               max_non_dual_feas_to_add_max;
   double            max_not_fixable_to_add_frac;
   int               max_not_fixable_to_add_min;
   int               max_not_fixable_to_add_max;

   int               mat_col_compress_num;
   double            mat_col_compress_ratio;
   int               mat_row_compress_num;
   double            mat_row_compress_ratio;

   /* parameters governing tailing off checking */
   int               tailoff_gap_backsteps;
   double            tailoff_gap_frac;
   int               tailoff_obj_backsteps;
   double            tailoff_obj_frac;
   double            tailoff_absolute;

   int               ineff_cnt_to_delete;
   int               eff_cnt_before_cutpool;
   int               ineffective_constraints;
   int               base_constraints_always_effective;

   int               branch_on_cuts;   /* TRUE / FALSE */
   int               discard_slack_cuts;

   cut_time_out      first_lp;
   cut_time_out      later_lp;

   int               max_cut_num_per_iter;

   /* Reduced cost and logical fixing parameters */
   int               do_reduced_cost_fixing;
   double            gap_as_ub_frac;
   double            gap_as_last_gap_frac;
   int               do_logical_fixing;
   int               fixed_to_ub_before_logical_fixing; /* OK */
   double            fixed_to_ub_frac_before_logical_fixing; /* OK */

   /* CGL parameters */
   cgl_params        cgl;
   
   /* Parameters affecting branching */
   int               max_presolve_iter;

   /*Defaults for the user supplied routines*/
   int               is_feasible_default;
   int               send_feasible_solution_default;
   int               display_solution_default;
   int               shall_we_branch_default;
   int               select_candidates_default;
   int               strong_branching_cand_num_min;
   int               strong_branching_cand_num_max;
   double            strong_branching_red_ratio;
   int               compare_candidates_default;
   int               select_child_default;
   int               pack_lp_solution_default;

   /* Multi-criteria parameters */
   int               multi_criteria;
   int               mc_find_supported_solutions;
   int               mc_add_optimality_cuts;
   double            mc_rho;   /* For augmented Chebyshev norm */
   double            mc_gamma; /* Weight on first objective */
   double            mc_tau;   /* Weight on second objective */

   int               sensitivity_analysis;

}lp_params;

#endif
=======
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

#ifndef _LP_PARAMS_H
#define _LP_PARAMS_H

#include "sym_constants.h"
#include "sym_timemeas.h"

/*---------------------------------------------------------------------------*\
| The list of parameters associated with processing a node in the branch and  |
| cut tree. See the README file for an explanation of the parameters          |
\*---------------------------------------------------------------------------*/

typedef struct CUT_TIME_OUT{
   double            first_cut_time_out;
   double            all_cuts_time_out;
}cut_time_out;

typedef struct CGL_PARAMS{
   /* Cut generation in LP */
   int               generate_cgl_cuts;
   int               generate_cgl_gomory_cuts;
   int               generate_cgl_redsplit_cuts;
   int               generate_cgl_knapsack_cuts;
   int               generate_cgl_oddhole_cuts;
   int               generate_cgl_probing_cuts;
   int               generate_cgl_mir_cuts;
   int               generate_cgl_twomir_cuts;
   int               generate_cgl_clique_cuts;
   int               generate_cgl_flow_and_cover_cuts;
   int               generate_cgl_rounding_cuts;
   int               generate_cgl_lift_and_project_cuts;
   int               generate_cgl_landp_cuts;

   int               generate_cgl_gomory_cuts_freq;
   int               generate_cgl_redsplit_cuts_freq;
   int               generate_cgl_knapsack_cuts_freq;
   int               generate_cgl_oddhole_cuts_freq;
   int               generate_cgl_probing_cuts_freq;
   int               generate_cgl_mir_cuts_freq;
   int               generate_cgl_twomir_cuts_freq;
   int               generate_cgl_clique_cuts_freq;
   int               generate_cgl_flow_and_cover_cuts_freq;
   int               generate_cgl_rounding_cuts_freq;
   int               generate_cgl_lift_and_project_cuts_freq;
   int               generate_cgl_landp_cuts_freq;

   int               gomory_generated_in_root;
   int               redsplit_generated_in_root;
   int               knapsack_generated_in_root;
   int               oddhole_generated_in_root;
   int               probing_generated_in_root;
   int               mir_generated_in_root;
   int               twomir_generated_in_root;
   int               clique_generated_in_root;
   int               flow_and_cover_generated_in_root;
   int               rounding_generated_in_root;
   int               lift_and_project_generated_in_root;
   int               landp_generated_in_root;
}cgl_params;

typedef struct LP_PARAMS{
   int               verbosity;
   double            granularity;
   int               use_cg;
   int               set_obj_upper_lim;
   int               do_primal_heuristic;
   double            time_limit;

   /* these two are passed directly to the lp solver */
   int               scaling;
   int               fastmip;

   int               try_to_recover_from_error;
   /* ZERO_ONE_PROBLEM / INTEGER_PROBLEM / MIXED_INTEGER_PROBLEM */
   int               problem_type;
   int               keep_description_of_pruned;

   int               not_fixed_storage_size;

   int               cut_pool_check_freq;

   int               load_balance_level;
   int               load_balance_iterations;
   int               load_balance_compare_candidates;

   double            fractional_diving_ratio;
   int               fractional_diving_num;

   /* parameters constraining the growth of the matrix */
   double            max_non_dual_feas_to_add_frac;
   int               max_non_dual_feas_to_add_min;
   int               max_non_dual_feas_to_add_max;
   double            max_not_fixable_to_add_frac;
   int               max_not_fixable_to_add_min;
   int               max_not_fixable_to_add_max;

   int               mat_col_compress_num;
   double            mat_col_compress_ratio;
   int               mat_row_compress_num;
   double            mat_row_compress_ratio;

   /* parameters governing tailing off checking */
   int               tailoff_gap_backsteps;
   double            tailoff_gap_frac;
   int               tailoff_obj_backsteps;
   double            tailoff_obj_frac;
   double            tailoff_absolute;

   int               ineff_cnt_to_delete;
   int               eff_cnt_before_cutpool;
   int               ineffective_constraints;
   int               base_constraints_always_effective;

   int               branch_on_cuts;   /* TRUE / FALSE */
   int               discard_slack_cuts;

   cut_time_out      first_lp;
   cut_time_out      later_lp;

   int               max_cut_num_per_iter;

   /* Reduced cost and logical fixing parameters */
   int               do_reduced_cost_fixing;
   double            gap_as_ub_frac;
   double            gap_as_last_gap_frac;
   int               do_logical_fixing;
   int               fixed_to_ub_before_logical_fixing; /* OK */
   double            fixed_to_ub_frac_before_logical_fixing; /* OK */

   /* CGL parameters */
   cgl_params        cgl;
   
   /* Parameters affecting branching */
   int               max_presolve_iter;

   /*Defaults for the user supplied routines*/
   int               is_feasible_default;
   int               send_feasible_solution_default;
   int               display_solution_default;
   int               shall_we_branch_default;
   int               select_candidates_default;
   int               strong_branching_cand_num_min;
   int               strong_branching_cand_num_max;
   double            strong_branching_red_ratio;
   int               compare_candidates_default;
   int               select_child_default;
   int               pack_lp_solution_default;

   /* Multi-criteria parameters */
   int               multi_criteria;
   int               mc_find_supported_solutions;
   int               mc_add_optimality_cuts;
   double            mc_rho;   /* For augmented Chebyshev norm */
   double            mc_gamma; /* Weight on first objective */
   double            mc_tau;   /* Weight on second objective */

   int               sensitivity_analysis;

}lp_params;

#endif
>>>>>>> .merge-right.r1068
