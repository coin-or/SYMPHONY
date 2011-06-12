/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2006-2011 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "CoinPragma.hpp"
#include "SymConfig.h"

#include <iostream>
#include <cstdio>

#ifdef COIN_HAS_OSITESTS
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiUnitTests.hpp"
#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSymSolverInterface.hpp"

using namespace OsiUnitTest;

#else
void testingMessage( const char * const msg ) {
  std::cout.flush() ;
  std::cerr <<msg;
}

#endif

#include "symphony.h"

int main (int argc, const char *argv[])
{
  std::string miplib3Dir;
  bool exception = false;
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

#ifdef COIN_HAS_OSITESTS
  outcomes.clear();

  /*
    Process command line parameters.
  */
  std::map<std::string,std::string> parms;
  std::map<std::string,int> ignorekeywords;
  ignorekeywords["-p"] = 1;
  if (processParameters(argc,argv,parms,ignorekeywords) == false)
  { return (1) ; }

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;
  miplib3Dir = parms["-miplib3Dir"];

  try {
    /*
      Test Osi{Row,Col}Cut routines.
    */
    {
      OsiSymSolverInterface symSi;
      symSi.setSymParam(OsiSymVerbosity, -1);
      testingMessage( "Now testing the OsiRowCut class with OsiSymSolverInterface\n\n");
      OsiRowCutUnitTest(&symSi,mpsDir);
    }
    {
      OsiSymSolverInterface symSi;
      symSi.setSymParam(OsiSymVerbosity, -1);
      testingMessage( "Now testing the OsiColCut class with OsiSymSolverInterface\n\n" );
      OsiColCutUnitTest(&symSi,mpsDir);
    }
    {
      OsiSymSolverInterface symSi;
      symSi.setSymParam(OsiSymVerbosity, -1);
      testingMessage( "Now testing the OsiRowCutDebugger class with OsiSymSolverInterface\n\n" );
      OsiRowCutDebuggerUnitTest(&symSi,mpsDir);
    }

    /*
      Run the OsiSym class test. This will also call OsiSolverInterfaceCommonUnitTest.
    */
    testingMessage( "Now testing OsiSymSolverInterface\n\n" );
    OsiSymSolverInterfaceUnitTest(mpsDir,netlibDir);

    /*
      We have run the specialised unit test.
      Check now to see if we need to run through the Netlib problems.
    */
    if (parms.find("-testOsiSolverInterface") != parms.end())
    {
      // Create vector of solver interfaces
      std::vector<OsiSolverInterface*> vecSi(1, new OsiSymSolverInterface);

      testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
      OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir);

      delete vecSi[0];
    }
    else {
      testingMessage( "***Skipped Testing of OsiSymSolverInterface on Netlib problems***\n" );
      testingMessage( "***use -testOsiSolverInterface to run them.***\n" );
    }
  } catch (CoinError& error) {
    std::cout.flush();
    std::cerr << "Caught CoinError exception: ";
    error.print(true);
    exception = true;
  }
#else
  /* a very light version of "parameter processing": check if user call with -miplib3Dir=<dir> */
  if( argc >= 2 && strncmp(argv[1], "-miplib3Dir", 11) == 0 )
    miplib3Dir = argv[1]+11;
#endif

  if (miplib3Dir.length() > 0) {
    int test_status;
    testingMessage( "Testing MIPLIB files\n" );

    sym_environment *env = sym_open_environment();
    sym_parse_command_line(env, argc, const_cast<char**>(argv));
    sym_test(env, &test_status);

#ifdef COIN_HAS_OSITESTS
    OSIUNITTEST_ASSERT_WARNING(test_status == 0, {}, "symphony", "testing MIPLIB");
#else
    if (test_status > 0)
      testingMessage( "Warning: some instances may not have returned a correct solution\n");
#endif
  }

  /*
    We're done. Report on the results.
  */
#ifdef COIN_HAS_OSITESTS
  std::cout.flush();
  outcomes.print();

  int nerrors;
  int nerrors_expected;
  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);

  if (nerrors > nerrors_expected)
    std::cerr << "Tests completed with " << nerrors - nerrors_expected << " unexpected errors." << std::endl ;
  else if( exception )
    std::cerr << "Tests completed with exception\n";
  else
    std::cerr << "All tests completed successfully\n";

  return (nerrors - nerrors_expected) + (exception ? 1 : 0);
#else

  testingMessage( "All tests completed successfully\n" );

  return 0;
#endif
}
