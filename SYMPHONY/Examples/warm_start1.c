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

  OsiSymSolverInterface si;
  si.parseCommandLine(argc, argv);
  si.loadProblem();
  si.setSymParam(OsiSymKeepWarmStart, true);
  si.setSymParam(OsiSymFindFirstFeasible, true);
  si.setSymParam(OsiSymSearchStrategy, DEPTH_FIRST_SEARCH);

  si.initialSolve();

  si.setSymParam(OsiSymFindFirstFeasible, false);
  si.setSymParam(OsiSymSearchStrategy, BEST_FIRST_SEARCH);

  si.resolve();
  
  return(0);

}

#else

#include "symphony_api.h"
  
int main(int argc, char **argv)
{    
     
   sym_environment *env = sym_open_environment();   
   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);
   
   sym_set_int_param(env, "keep_warm_start", TRUE);

   sym_set_int_param(env, "find_first_feasible", TRUE);
   sym_set_int_param(env, "node_selection_rule", DEPTH_FIRST_SEARCH);

   sym_solve(env);

   sym_set_int_param(env, "find_first_feasible", FALSE);
   sym_set_int_param(env, "node_selection_rule", BEST_FIRST_SEARCH);


   sym_warm_solve(env);
   sym_close_environment(env);
  
   return(0);

}  

#endif

