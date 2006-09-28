/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
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

#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{
   int termcode;
   OsiSymSolverInterface si;

   si.parseCommandLine(argc, argv);
   si.loadProblem();

   si.setSymParam(OsiSymTimeLimit, 10);
   si.setSymParam(OsiKeepWarmStart, 1);
   termcode = si.initialSolve();

   while (termcode != TM_OPTIMAL_SOLUTION_FOUND){
      printf("Starting problem again from warm start...\n");
      termcode = si.resolve();
   }

   return(0);
}

#else

#include "symphony_api.h"
  
int main(int argc, char **argv)
{
   int termcode;
   
   sym_environment *env = sym_open_environment();
   sym_parse_command_line(env, argc, argv);
   sym_load_problem(env);

   sym_set_int_param(env, "time_limit", 10);
   sym_set_int_param(env, "keep_warm_start", 1);
   termcode = sym_solve(env);

   while (termcode != TM_OPTIMAL_SOLUTION_FOUND){
      printf("Starting problem again from warm start...\n");
      termcode = sym_warm_solve(env);
   }

   sym_close_environment(env);
}

#endif

