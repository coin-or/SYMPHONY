/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2004 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"

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
  ind[1] = 7; val[0] = 6000;
  
  double lb = si.getLbForNewRhs(2, ind, val);
  double ub =  si.getUbForNewRhs(2, ind, val);

  printf("\nBounds for the new rhs:\n lb: %f\n ub: %f \n\n", lb, ub);
  
  return(0);

}

#else

#include "symphony_api.h"
  
int main(int argc, char **argv)
{    
     
   sym_environment *env = sym_open_environment();   
   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);

   sym_set_int_param(env, "sensitivity_analysis", TRUE);

   sym_solve(env);

   int ind[2];
   double val[2];
   ind[0] = 4; val[0] = 7000;
   ind[1] = 7; val[0] = 6000;
   
   double lb, ub; 
   sym_get_lb_for_new_rhs(env, 2, ind, val, &lb);
   sym_get_ub_for_new_rhs(env, 2, ind, val, &ub);

   printf("\nBounds for the new rhs:\n lb: %f\n ub: %f \n\n", lb, ub);

   sym_close_environment(env);
  
   return(0);
}  

#endif

