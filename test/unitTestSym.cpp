/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2006-2022 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "SymConfig.h"

#include <iostream>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "OsiUnitTests.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSymSolverInterface.hpp"

using namespace OsiUnitTest;

#include "symphony.h"

/* The LP interface tests call internal SYMPHONY functions, which are not
   exported from a Windows DLL, and need the LP compiled into the library */
#if defined(SYM_COMPILE_IN_LP) && !defined(_WIN32)
#define SYM_TEST_LP_INTERFACE
#include <cmath>
#include <cstring>
#include "sym_lp.h"
#endif

#ifdef SYM_TEST_LP_INTERFACE
/*===========================================================================*/
/* Check the bounds of row 'row' in the LP solver. Use lb = -COIN_DBL_MAX for
   a row without a lower bound */

static bool rowBoundsAre(LPdata *lp_data, int row, double lb, double ub)
{
   OsiSolverInterface *si = lp_data->si;
   if (row >= si->getNumRows()){
      return false;
   }
   double rlb = si->getRowLower()[row], rub = si->getRowUpper()[row];
   bool lb_ok = (lb <= -COIN_DBL_MAX ? rlb <= -si->getInfinity() :
		 fabs(rlb - lb) < 1e-9);
   if (!lb_ok || fabs(rub - ub) >= 1e-9){
      std::cerr << "row " << row << " has bounds [" << rlb << ", " << rub
		<< "], expected [" << lb << ", " << ub << "]" << std::endl;
      return false;
   }
   return true;
}

/*===========================================================================*/
/* SYMPHONY stores a ranged row ('R') of a cut or branching object as
   rhs <= ax <= rhs + range if range >= 0 and as rhs + range <= ax <= rhs
   otherwise. Osi expects rhs to be the upper bound and range to be
   upper - lower. Check that the LP interface translates between the two
   when it changes rows, restores them from their cut data and adds them. */

static void testLpInterfaceRangedRows()
{
   /* One row, 5 <= x0 + x1 <= 8, loaded in Osi's convention */
   int matbeg[3] = {0, 1, 2};
   int matind[2] = {0, 0};
   double matval[2] = {1.0, 1.0};
   double obj[2] = {0.0, 0.0}, lb[2] = {0.0, 0.0}, ub[2] = {10.0, 10.0};
   double rhs[1] = {8.0}, rngval[1] = {3.0};
   char sense[1] = {'R'};
   MIPdesc mip;
   memset(&mip, 0, sizeof(MIPdesc));
   mip.n = 2;
   mip.m = 1;
   mip.nz = 2;
   mip.matbeg = matbeg;
   mip.matind = matind;
   mip.matval = matval;
   mip.obj = obj;
   mip.lb = lb;
   mip.ub = ub;
   mip.rhs = rhs;
   mip.rngval = rngval;
   mip.sense = sense;

   LPdata *lp_data = (LPdata *) calloc(1, sizeof(LPdata));
   open_lp_solver(lp_data);
   lp_data->mip = &mip;
   lp_data->n = mip.n;
   lp_data->m = mip.m;
   lp_data->nz = mip.nz;
   size_lp_arrays(lp_data, FALSE, TRUE, mip.m, mip.n, mip.nz);
   load_lp_prob(lp_data, FALSE, FALSE);
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 0, 5.0, 8.0), {},
			    "symphony", "LP interface: load ranged row");

   /* change_row() with both forms of a ranged row */
   change_row(lp_data, 0, 'R', 8.0, -3.0);
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 0, 5.0, 8.0), {},
			    "symphony", "LP interface: change_row, range < 0");
   change_row(lp_data, 0, 'R', 2.0, 4.0);
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 0, 2.0, 6.0), {},
			    "symphony", "LP interface: change_row, range > 0");

   /* constrain_row_set() restores a row from its cut data */
   cut_data cut;
   memset(&cut, 0, sizeof(cut_data));
   cut.sense = 'R';
   cut.rhs = 8.0;
   cut.range = -3.0;
   lp_data->rows[0].cut = &cut;
   int row = 0;
   constrain_row_set(lp_data, 1, &row);
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 0, 5.0, 8.0), {},
			    "symphony",
			    "LP interface: constrain_row_set, range < 0");
   cut.rhs = 5.0;
   cut.range = 3.0;
   constrain_row_set(lp_data, 1, &row);
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 0, 5.0, 8.0), {},
			    "symphony",
			    "LP interface: constrain_row_set, range > 0");

   /* add_waiting_rows() adds the rows and then sets the ranges of the
      ranged ones */
   lp_prob *p = (lp_prob *) calloc(1, sizeof(lp_prob));
   p->lp_data = lp_data;
   int ind_x0[1] = {0}, ind_x1[1] = {1}, ind_both[2] = {0, 1};
   double val_one[1] = {1.0}, val_both[2] = {1.0, 1.0};
   cut_data cuts[3];
   waiting_row wrows[3], *wrow_ptrs[3];
   memset(cuts, 0, sizeof(cuts));
   memset(wrows, 0, sizeof(wrows));
   cuts[0].sense = 'R';   /* 2 <= x0 <= 6 */
   cuts[0].rhs = 2.0;
   cuts[0].range = 4.0;
   cuts[1].sense = 'L';   /* x1 <= 9 */
   cuts[1].rhs = 9.0;
   cuts[2].sense = 'R';   /* 5 <= x0 + x1 <= 7 */
   cuts[2].rhs = 7.0;
   cuts[2].range = -2.0;
   wrows[0].matind = ind_x0;
   wrows[0].matval = val_one;
   wrows[0].nzcnt = 1;
   wrows[1].matind = ind_x1;
   wrows[1].matval = val_one;
   wrows[1].nzcnt = 1;
   wrows[2].matind = ind_both;
   wrows[2].matval = val_both;
   wrows[2].nzcnt = 2;
   for (int i = 0; i < 3; i++){
      wrows[i].cut = cuts + i;
      wrow_ptrs[i] = wrows + i;
   }
   add_waiting_rows(p, wrow_ptrs, 3);
   OSIUNITTEST_ASSERT_ERROR(lp_data->m == 4 &&
			    lp_data->si->getNumRows() == 4, {}, "symphony",
			    "LP interface: add_waiting_rows, row count");
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 0, 5.0, 8.0), {},
			    "symphony",
			    "LP interface: add_waiting_rows, existing row");
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 1, 2.0, 6.0), {},
			    "symphony",
			    "LP interface: add_waiting_rows, range > 0");
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 2, -COIN_DBL_MAX, 9.0), {},
			    "symphony",
			    "LP interface: add_waiting_rows, 'L' row");
   OSIUNITTEST_ASSERT_ERROR(rowBoundsAre(lp_data, 3, 5.0, 7.0), {},
			    "symphony",
			    "LP interface: add_waiting_rows, range < 0");

   free(p);
   lp_data->mip = NULL;
   free(lp_data->rows);
   lp_data->rows = NULL;
   free_lp_arrays(lp_data);
   close_lp_solver(lp_data);
   free(lp_data);
}
#endif

/*===========================================================================*/

int main (int argc, const char *argv[])
{
#ifdef _OPENMP
   omp_set_dynamic(FALSE);
   omp_set_num_threads(1);
#endif

  std::string miplib3Dir;
  /*
    Start off with various bits of initialisation that don't really belong
    anywhere else.

    Synchronise C++ stream i/o with C stdio. This makes debugging
    output a bit more comprehensible. It still suffers from interleave of cout
    (stdout) and cerr (stderr), but -nobuf deals with that.
  */
  std::ios::sync_with_stdio() ;
  /*
    Suppress an popup window that Windows shows in response to a crash. See
    note at head of file.
  */
  WindowsErrorPopupBlocker();

  /*
    Process command line parameters.
  */
  std::map<std::string,std::string> parms;
  if (processParameters(argc,argv,parms) == false)
  { return 1; }

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;
  miplib3Dir = parms["-miplib3Dir"];

#ifdef SYM_TEST_LP_INTERFACE
  testingMessage( "Now testing ranged rows in the LP interface\n\n" );
  OSIUNITTEST_CATCH_ERROR(testLpInterfaceRangedRows(), {}, "symphony",
			  "LP interface ranged rows");
#endif

  /*
    Test Osi{Row,Col}Cut routines.
   */
  {
    OsiSymSolverInterface symSi;
    symSi.setSymParam(OsiSymVerbosity, -1);
    testingMessage( "Now testing the OsiRowCut class with OsiSymSolverInterface\n\n");
    OSIUNITTEST_CATCH_ERROR(OsiRowCutUnitTest(&symSi,mpsDir), {}, "symphony", "rowcut unittest");
  }
  {
    OsiSymSolverInterface symSi;
    symSi.setSymParam(OsiSymVerbosity, -1);
    testingMessage( "Now testing the OsiColCut class with OsiSymSolverInterface\n\n" );
    OSIUNITTEST_CATCH_ERROR(OsiColCutUnitTest(&symSi,mpsDir), {}, "symphony", "colcut unittest");
  }
  {
    OsiSymSolverInterface symSi;
    symSi.setSymParam(OsiSymVerbosity, -1);
    testingMessage( "Now testing the OsiRowCutDebugger class with OsiSymSolverInterface\n\n" );
    OSIUNITTEST_CATCH_ERROR(OsiRowCutDebuggerUnitTest(&symSi,mpsDir), {}, "symphony", "rowcut debugger unittest");
  }

  /*
    Run the OsiSym class test. This will also call OsiSolverInterfaceCommonUnitTest.
   */
  testingMessage( "Now testing OsiSymSolverInterface\n\n" );
  OSIUNITTEST_CATCH_ERROR(OsiSymSolverInterfaceUnitTest(mpsDir,netlibDir), {}, "symphony", "Osi unittest");

  /*
    We have run the specialised unit test.
    Check now to see if we need to run through the Netlib problems.
   */
  if (parms.find("-testOsiSolverInterface") != parms.end())
  {
    // Create vector of solver interfaces
    OsiSymSolverInterface* symSi = new OsiSymSolverInterface();
    symSi->setSymParam(OsiSymVerbosity, -1);
    std::vector<OsiSolverInterface*> vecSi(1, symSi);

    testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
    OSIUNITTEST_CATCH_ERROR(OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir), {}, "symphony", "Netlib unittest");

    delete vecSi[0];
  }
  else
  {
    testingMessage( "***Skipped Testing of OsiSymSolverInterface on Netlib problems, use -testOsiSolverInterface to run them.***\n" );
  }

  if (miplib3Dir.length() > 0) {
    int test_status;
    int symargc;
    char* symargv[7];
    testingMessage( "Testing MIPLIB files\n" );

    sym_environment *env = sym_open_environment();
    /* assemble arguments for symphony: -T miplibdir, and -p 2 if we run the punittest */
    symargc = 5;
    symargv[0] = CoinStrdup(argv[0]);
    symargv[1] = "-T";
    symargv[2] = CoinStrdup(miplib3Dir.c_str());
    if( argv[0][0] == 'p' || argv[0][0] == 'P' ) {
       symargv[3] = "-p";
       symargv[4] = "2";
    }else{
       symargv[3] = "-p";
       symargv[4] = "0";
    }
    //sym_parse_command_line(env, symargc, const_cast<char**>(symargv));
    sym_set_int_param(env, "verbosity", -10);
    sym_test(env, symargc, symargv, &test_status);

    OSIUNITTEST_ASSERT_WARNING(test_status == 0, {}, "symphony", "testing MIPLIB");
  }

  /*
    We're done. Report on the results.
  */
  std::cout.flush();
  outcomes.print();

  int nerrors;
  int nerrors_expected;
  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);

  if (nerrors > nerrors_expected)
    std::cerr << "Tests completed with " << nerrors - nerrors_expected << " unexpected errors." << std::endl ;
  else
    std::cerr << "All tests completed successfully\n";

  return nerrors - nerrors_expected;
}
