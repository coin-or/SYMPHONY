/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER
//#define TEST_MULTI_CRITERIA
//#define TEST_RESOLVE
//#define TEST_SENS_ANALYSIS

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


#ifdef TEST_MULTI_CRITERIA

   /* Test for mc.mps */
   si.setObj2Coeff(1, 1);

   si.setSymParam(OsiSymMultiCriteriaFindNondominatedSolutions, FALSE);
   
   /* Solve the multi-criteria problem */
   si.multiCriteriaBranchAndBound();
#endif
   
#if defined TEST_RESOLVE || defined TEST_WARM_START
   si.setSymParam(OsiSymKeepWarmStart, TRUE); 
   si.setSymParam(OsiSymDoReducedCostFixing, FALSE);
#endif

#if defined TEST_WARM_START
   /* testing for a generic problem */
   si.setSymParam(OsiSymFindFirstFeasible, true);
   si.setSymParam(OsiSymSearchStrategy, DEPTH_FIRST_SEARCH);
#endif

#if defined TEST_SENS_ANALYSIS
   si.setSymParam(OsiSymSensitivityAnalysis, TRUE);
#endif

#ifndef TEST_MULTI_CRITERIA
   si.branchAndBound();

   if(si.isIterationLimitReached()){
      cout<<"Node Limit Reached!"<<endl;
   }
   if(si.isProvenOptimal()){
      cout<<"Optimal Solution Found!"<<endl;
   }
#endif

#ifdef TEST_RESOLVE
   /* test for MIPLIB's p0201 */
   si.setObjCoeff(0, 100);
   si.setObjCoeff(200, 150);

   printf("RESOLVING...\n");
   si.resolve();
#endif

#ifdef TEST_WARM_START
   /* test for a generic problem */
   CoinWarmStart * sWS = si.getWarmStart();
   si.setSymParam(OsiSymFindFirstFeasible, false);
   printf("SOLVING THE ORIGINAL PROBLEM TO OPTIMALITY\n"); 
   si.resolve();
   printf("WARMSTARTING...\n");   
   si.setSymParam(OsiSymSearchStrategy, BEST_FIRST_SEARCH);
   si.setWarmStart(sWS);
   si.resolve();
#endif

#ifdef TEST_SENS_ANALYSIS

   /* test for MIPLIB's flugpl */
   int cnt = 2;
   int * ind = (int*) malloc(ISIZE*cnt);
   double * val = (double*) malloc (DSIZE*cnt); 
   double lb;
   double ub;

   ind[0] = 4;  val[0] = 7000;
   ind[1] = 7;  val[1] = 6000;

   //   lb = si.getLbForNewRhs(cnt, ind, val);
   ub = si.getUbForNewRhs(cnt, ind, val);

   //   printf("LB obtained for new rhs problem: %f \n\n\n",lb);
   printf("UB obtained for new rhs problem: %f \n\n\n",ub);


   lb = si.getLbForNewObj(cnt, ind, val);
   ub = si.getUbForNewObj(cnt, ind, val);
      
   printf("LB obtained for new obj problem: %f \n\n\n",lb);
   printf("UB obtained for new obj problem: %f \n\n\n",ub);
	
#endif

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
