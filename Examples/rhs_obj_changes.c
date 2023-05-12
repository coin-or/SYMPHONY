#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <time.h>

#include "symphony.h"


// Helper function
double* generate_rand_array(int n, int l, int u) {
  double* res = (double *)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    res[i] = (rand() % (u - l + 1)) + l;
  }
  return res;
}


int main(int argc, char **argv)
{   
      
   int termcode;
   int n = 4; // number of variables
   int m = 2; // number of constraints 

   double *rhs = (double *) malloc(m * sizeof(double));
   double *coeff = (double *) malloc(n * sizeof(double));

   double *coldObjVal = (double *) malloc(sizeof(double));
   double *warmObjVal = (double *) malloc(sizeof(double));
   
   // Create SYMPHONY environment and load the problem 
   sym_environment *env = sym_open_environment(); 
   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);

   // Set the parameters for the warm starting
   sym_set_int_param(env, "keep_warm_start", TRUE);
   sym_set_int_param(env, "should_use_rel_br", FALSE);
   sym_set_int_param(env, "use_hot_starts", FALSE);
   sym_set_int_param(env, "should_warmstart_node", TRUE);
   sym_set_int_param(env, "sensitivity_analysis", TRUE);
   sym_set_int_param(env, "set_obj_upper_lim", FALSE);
   sym_set_int_param(env, "do_primal_heuristic", FALSE);
   sym_set_int_param(env, "prep_level", -1);
   sym_set_int_param(env, "tighten_root_bounds", FALSE);
   sym_set_int_param(env, "max_sp_size", 100);
   sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env, "generate_cgl_cuts", FALSE);
   sym_set_int_param(env, "max_active_nodes", 1);

   // First solve the original problem
   // This will create a warm start for the next solve
   if ((termcode = sym_solve(env)) < 0){
      printf("PROBLEM INFEASIBLE!\n");
      return (1);
   }

   // Set the new RHS
   rhs = generate_rand_array(m, -1000, -1);
   for (int i = 0; i < m; i++){
      sym_set_row_upper(env, i, rhs[i]);
   }

   // Set the new Obj Coeff
   coeff = generate_rand_array(n, 1, 1000);
   for (int i = 0; i < n; i++){
      sym_set_obj_coeff(env, i, coeff[i]);
   }

   // Warm solve
   if ((termcode = sym_warm_solve(env)) < 0){
      printf("PROBLEM INFEASIBLE!\n");
      return(1);
   } 
   // Get the objective value
   sym_get_obj_val(env, warmObjVal);

   // Cold solve to check if the values coincide
   if ((termcode = sym_solve(env)) < 0){
      printf("PROBLEM INFEASIBLE!\n");
      return(1);
   } 

   sym_get_obj_val(env, coldObjVal);

   // These must be equal
   assert((*warmObjVal) == (*coldObjVal));

   sym_close_environment(env);
   return 0;
}  