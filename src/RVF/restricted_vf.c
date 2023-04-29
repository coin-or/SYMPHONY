#include <cstddef>
#include <cstdio>
#include <cstdlib>

// feb223
// should we use the Osi Interface?
// #include "OsiSymSolverInterface.hpp"

#include "symphony.h"
  
int main(int argc, char **argv)
{    
   double inf = 1e30;
   int i = 0, termcode;
   double lowerRhs = -4;
   double upperRhs = 4;
   double step = .5;
   double rhs = lowerRhs;

   int n = 4;

   int length = ((upperRhs - lowerRhs) / step) + 1;
   double *vf = (double *) malloc(length * sizeof(double));
   double *currObjVal = (double *) malloc(sizeof(double));
   double *currSol = (double *) malloc(n * sizeof(double));
   double *collb = (double *) malloc(n * sizeof(double));

   sym_environment *env = sym_open_environment();   
   
   // Maybe we need to store the warm_start_desc
   // warm_start_desc * ws; 

   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);

   sym_set_int_param(env, "verbosity", -2);
   // feb223
   // Those are the parameters to be set in order to 
   // keep the branch-and-bound tree valid for RHS changes
   // (E.g. Cuts and reduced cost fixing are not RHS-invariant)
   // We should check on a real example that each of this setting
   // is actually needed
   sym_set_int_param(env, "prep_level", -1);
   sym_set_int_param(env, "keep_warm_start", TRUE);
   sym_set_int_param(env, "should_use_rel_br", FALSE);
   sym_set_int_param(env, "use_hot_starts", FALSE);
   sym_set_int_param(env, "should_warmstart_node", TRUE);
   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "sensitivity_bounds", TRUE);
   sym_set_int_param(env, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env, "do_primal_heuristic", FALSE);
   sym_set_int_param(env, "tighten_root_bounds", FALSE);
   sym_set_int_param(env, "max_sp_size", 100);
   sym_set_int_param(env, "generate_cgl_cuts", FALSE);
   sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
   // sym_set_int_param(env, "max_active_nodes", 1);

   // sym_write_lp(env, "LP_initial");

   // Which RHS to start from? (0?)

   // First solve
   if ((termcode = sym_solve(env)) < 0){
      printf("PROBLEM INFEASIBLE!\n");
      return(0);
   }

   // -----------------------------------------------------
   // Main loop
   // -----------------------------------------------------
   while (rhs <= upperRhs){

      // Set the new RHS
      sym_set_row_upper(env, 0, rhs);
      sym_set_row_lower(env, 0, rhs);
      
      printf("SOLVING RHS = %.2f! \n", rhs);
      sym_write_lp(env, "LP");

      // Warm solve
      if ((termcode = sym_warm_solve(env)) < 0){
         printf("PROBLEM INFEASIBLE!\n");
         vf[i] = inf;
      } else {
         // You can get here every information about the solution
         // e.g. the integer part and the value function
         sym_get_col_solution(env, currSol);
         sym_get_col_lower(env, collb);

         sym_get_obj_val(env, currObjVal);
         vf[i] = *currObjVal;
      }
      
      // Branching decision on the RHS
      // right now just increase RHS by stepsize
      rhs += step;
      i++;
   }

   sym_close_environment(env);

   // Here we will have the value function
   // print something
   rhs = lowerRhs;
   for (i = 0; i < length; i++ ){
      printf("rhs = %.2f | vf(rhs) = %.2f\n", rhs, vf[i]);
      rhs += step;
   }

   return 0;
}  