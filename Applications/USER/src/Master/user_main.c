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
#define TEST_MULTI_CRITERIA

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the SYMPHONY generic MIP solver.
\*===========================================================================*/

#include "OsiSymSolverInterface.hpp"

#if 0
int main(int argc, char **argv)
{
   OsiSymSolverInterface si;

   si.parseCommandLine(argc, argv);

   si.loadProblem();

   si.setSymParam(OsiSymKeepDescOfPruned, KEEP_IN_MEMORY); 

   si.setSymParam(OsiSymNodeLimit, 100);

   si.branchAndBound();

   CoinWarmStart * sWS = si.getWarmStart();

   si.setSymParam(OsiSymNodeLimit, 0);
   
   si.branchAndBound();

   //si.setWarmStart(sWS);

   /*test for p0201 */
   si.setObjCoeff(0, 100);
   si.setObjCoeff(200, 150);

   //si.resolve();

   si.branchAndBound();

   return(0);
}
#endif

int main(int argc, char **argv)
{


  OsiSymSolverInterface si;
   //si.setSymParam(OsiSymVerbosity, -1);

   /* Parse the command line */
   si.parseCommandLine(argc, argv);
   
   /* Read in the problem */
   si.loadProblem();

   /* Find a priori problem bounds */
   si.findInitialBounds();

#ifdef TEST_MULTI_CRITERIA


   /* Test for p0033   
   int i;
   for(i = 0; i<10; i++)
      si.setObj2Coeff(i, 50);

   for(i = 10; i<33; i++)
      si.setObj2Coeff(i, -50);
   */

   /* Test for p0033   
      si.setObj2Coeff(0, 200); 
      si.setObj2Coeff(6, 200);  
      si.setObj2Coeff(20, 300);   
      si.setObj2Coeff(25, 300); 
   */

   /* test for lseu 
   si.setObj2Coeff(0, 200);
   si.setObj2Coeff(10, -200);
   si.setObj2Coeff(50, 200);
   si.setObj2Coeff(88, -10);
   */

   /* Test for flugpl 
   si.setObj2Coeff(0, 2000);
   si.setObj2Coeff(3, 2000);
   si.setObj2Coeff(6, 2000);
   si.setObj2Coeff(17, 40);
   */

   /* Test for dc_multi */
   //si.setObjCoeff(1, 8);
   si.setObj2Coeff(1, -1); 
   //   si.setObj2Coeff(30, 100);

   si.branchAndBound();

   si.setSymParam(OsiSymMultiCriteriaFindNondominatedSolutions, FALSE);
   
   /* Solve the multi-criteria problem */
   //si.multiCriteriaBranchAndBound();
#endif
   
#if defined TEST_RESOLVE || defined TEST_SENS_ANALYSIS || \
   defined TEST_WARM_START
   si.setSymParam(OsiSymKeepDescOfPruned, KEEP_IN_MEMORY); 
#endif

#if defined TEST_WARM_START
   si.setSymParam(OsiSymNodeLimit, 5);
#endif

   //si.branchAndBound();

#ifdef TEST_RESOLVE
   si.setSymParam(OsiSymWarmStart, TRUE);    



  /*test for dsbmip 
   si.setObjCoeff(0, -0.001);
   si.setObjCoeff(1693, -0.001);
  */ 

   /*test for p0201 
   si.setObjCoeff(0, 100);
   si.setObjCoeff(200, 150);
   */

   /* test for flugpl 
   si.setObjCoeff(0, 2000);
   si.setObjCoeff(17, 10);
   */


   /* test for flugpl    
   si.setObjCoeff(0, 2000);
   si.setObjCoeff(1, 1400);
   si.setObjCoeff(3, 2000);
   si.setObjCoeff(6, 2000);
   si.setObjCoeff(13, 1000);
   si.setObjCoeff(17, 40);
   */

   printf("RESOLVING...\n");
   //   si.resolve();
#endif

#ifdef TEST_WARM_START
   CoinWarmStart * sWS = si.getWarmStart();
   si.setWarmStart(sWS);
   si.setSymParam(OsiSymNodeLimit, 100);
   si.branchAndBound();
#endif

#ifdef TEST_SENS_ANALYSIS

   int cnt = 2;
   int * ind = (int*) malloc(ISIZE*cnt);
   double * val = (double*) malloc (DSIZE*cnt); 
   double lb;
   
   //   ind[0] = 1;  val[0] = (atoi)(argv[1]);
   //   ind[1] = 4;  val[1] = (atoi)(argv[2]);

   ind[0] = 1;  val[0] = (atoi)(argv[1]);
   ind[1] = 4;  val[1] = (atoi)(argv[2]);

   lb = si.getLbForNewRhs(cnt, ind, val);
	
   printf("LB obtained for new rhs problem: %f \n\n\n",lb);
		      
#endif

   return(0);
}

#if 0
#include "master.h"

int main(int argc, char **argv)
{
   problem *p = sym_open_environment();

   sym_parse_command_line(p, argc, argv);

   sym_load_problem(p);

   sym_find_initial_bounds(p);

   sym_solve(p);

   sym_close_environment(p);

   return(0);
}

#endif
