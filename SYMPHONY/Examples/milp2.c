/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2005-2006 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*-------------------------------------------------------------------------*/
/*
  This is an example of using SYMPHONY to construct and solve 
  a simple MILP.

  optimal solution: x* = (1,1)
  
  minimize -1 x0 - 1 x1
  s.t       1 x0 + 2 x1 <= 3
            2 x0 + 1 x1 <= 3
              x0        >= 0 integer
              x1        >= 0 integer
*/

/*-------------------------------------------------------------------------*/

#include "symphony.h"
#include <iostream>

int main(int argc, char* argv[]){

   /* Create a SYMPHONY environment */
   
   sym_environment *env = sym_open_environment();

   int n_cols = 2;
   double * objective    = new double[n_cols];//the objective coefficients
   double * col_lb       = new double[n_cols];//the column lower bounds
   double * col_ub       = new double[n_cols];//the column upper bounds
    
   //Define the objective coefficients.
   //minimize -1 x0 - 1 x1
   objective[0] = -1.0;
   objective[1] = -1.0;
   
   //Define the variable lower/upper bounds.
   // x0 >= 0   =>  0 <= x0 <= infinity
   // x1 >= 0   =>  0 <= x1 <= infinity
   col_lb[0] = 0.0;
   col_lb[1] = 0.0;
   col_ub[0] = sym_get_infinity();
   col_ub[1] = sym_get_infinity();
   
   int n_rows = 2;
   char * row_sense = new char[n_rows]; //the row senses
   double * row_rhs = new double[n_rows]; //the row right-hand-sides
   double * row_range = NULL; //the row ranges   
   row_sense[0] = 'L';
   row_rhs[0] = 3;
   row_sense[1] = 'L';
   row_rhs[1] = 3;

   /* Constraint matrix definitions */
   int non_zeros = 4;
   int * start = new int[n_cols + 1]; 
   int * index = new int[non_zeros];
   double * value = new double[non_zeros];

   start[0] = 0; 
   start[1] = 2;
   start[2] = 4;

   index[0] = 0;
   index[1] = 1;
   index[2] = 0;
   index[3] = 1;

   value[0] = 1;
   value[1] = 2;
   value[2] = 2;
   value[3] = 1;

   //define the integer variables

   char * int_vars = new char[n_cols];

   int_vars[0] = TRUE;
   int_vars[1] = TRUE;

   //load the problem to environment
   sym_explicit_load_problem(env, n_cols, n_rows, start, index, value, col_lb, 
			     col_ub, int_vars, objective, NULL, row_sense, 
			     row_rhs, row_range, TRUE);
 
   //solve the integer program
   sym_solve(env);
   
   //get, print the solution
   double * solution = new double[n_cols];
   double objective_value = 0.0;

   sym_get_col_solution(env, solution);
   sym_get_obj_val(env, &objective_value);

   std::cout<<"The optimal solution is"
	<< " x0 = " << solution[0]
	<< " x1 = " << solution[1]
	<< " with objective value = " << objective_value << std::endl;
   
   //free the memory
   sym_close_environment(env);

   if(objective != 0)   { delete [] objective; objective = 0; }
   if(col_lb != 0)      { delete [] col_lb; col_lb = 0; }
   if(col_ub != 0)      { delete [] col_ub; col_ub = 0; }
   if(row_rhs != 0)      { delete [] row_rhs; row_rhs = 0; }
   if(row_sense != 0)      { delete [] row_sense; row_sense = 0; }
   if(row_range != 0)      { delete [] row_range; row_range = 0; }
   if(index != 0)      { delete [] index; index = 0; }
   if(start != 0)      { delete [] start; start = 0; }
   if(value != 0)      { delete [] value; value = 0; }
   if(int_vars != 0)      { delete [] int_vars; int_vars = 0; }

   return 0;

};

