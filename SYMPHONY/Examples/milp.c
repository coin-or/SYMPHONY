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

   /* Create an OsiSym object */
   OsiSymSolverInterface si;

   /* Parse the command line */
   si.parseCommandLine(argc, argv);
   
   /* Read in the problem */
   si.loadProblem();

   /* Solve the problem */
   si.branchAndBound();

   return(0);
}

#else

#include "symphony_api.h"
  
int main(int argc, char **argv)
{    
     
   sym_environment *env = sym_open_environment();
   
   sym_parse_command_line(env, argc, argv);
   
   sym_load_problem(env);
   
   sym_solve(env);
   
   sym_close_environment(env);
  
   return(0);
}  

#endif
