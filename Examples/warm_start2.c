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
  CoinWarmStart ws;

  si.parseCommandLine(argc, argv);
  si.loadProblem();

  si.setSymParam(OsiSymNodeLimit, 100);

  si.initialSolve();
  ws = si.getWarmStart();
  si.setSymParam(OsiSymNodeLimit, -1);

  si.resolve();

  si.setObjCoeff(0, 1);
  si.setObjCoeff(200, 150);
  si.setWarmStart(ws);

  si.resolve();

  return(0);

}

#else

#include "symphony_api.h"
  
int main(int argc, char **argv)
{    
     
   sym_environment *env = sym_open_environment();   
   warm_start_desc * ws; 

   sym_parse_command_line(env, argc, argv);   
   sym_load_problem(env);

   sym_set_int_param(env, "node_limit", 100);

   sym_solve(env);
   ws = sym_get_warm_start(env, true);

   sym_set_int_param(env, "node_limit", -1);

   sym_warm_solve(env);

   sym_set_obj_coeff(env, 0, 1);
   sym_set_obj_coeff(env, 200, 150);

   sym_set_warm_start(env, ws);
   sym_warm_solve(env);

   sym_close_environment(env);
  
   return(0);

}  

#endif

