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

/*===========================================================================*\
 * This file contains the main() for the SYMPHONY generic MIP solver.
 * Note that, if you want to use the OSI SYMPHONY interface, you should set the
 * USE_OSI_INTERFACE flag and define the COINROOT path in the SYMPHONY 
 * Makefile. Otherwise, the C callable library functions will be used by 
 * default. See below for the usage.
\*===========================================================================*/

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

   /* Find a priori problem bounds */
   si.findInitialBounds();

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

   sym_find_initial_bounds(env);

   sym_solve(env);

   sym_close_environment(env);

   return(0);
}

#endif
