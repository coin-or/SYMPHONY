/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2006-2008 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifdef HAVE_WINDOWS_H
#include <windows.h>
#endif

#include <cassert>
#include <iostream>
#include <cstdio>

#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSymSolverInterface.hpp"
#include "symphony.h"

void testingMessage( const char * const msg );

int main (int argc, const char *argv[])
{
   int i;
   const char dirsep =  CoinFindDirSeparator();
   int test_status = 0;
   
#ifdef HAVE_WINDOWS_H
   SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
#endif

   // define valid parameter keywords
   std::set<std::string> definedKeyWords;
   definedKeyWords.insert("-mpsDir");
   definedKeyWords.insert("-netlibDir");
   definedKeyWords.insert("-testOsiSolverInterface");

   std::map<std::string,std::string> parms;
   for ( i=1; i<argc; i++ ) {
      std::string parm(argv[i]);
      std::string key,value;
      unsigned int  eqPos = parm.find('=');
      
      // Does parm contain an '='
      if ( eqPos==std::string::npos ) {
	 //Parm does not contain '='
	 key = parm;
      }
      else {
	 key=parm.substr(0,eqPos);
	 value=parm.substr(eqPos+1);
      }
      parms[key]=value;
   }
   
   std::string mpsDir;
   std::string netlibDir;
   
   if (parms.find("-mpsDir") != parms.end()){
      mpsDir=parms["-mpsDir"] + dirsep;
   }else{
      //#ifdef _MSC_VER
      //      mpsDir = "..\\..\\Data\\Sample\\";
      //#else
      mpsDir = dirsep =='/' ? "../../Data/Sample/" : "..\\..\\Data\\Sample\\";
      //#endif  
   } 
   if (parms.find("-netlibDir") != parms.end()){
      netlibDir=parms["-netlibDir"] + dirsep;
   }else{ 
      //#ifdef _MSC_VER
      //     netlibDir = "..\\..\\Data\\Netlib\\";
      //#else
      netlibDir = dirsep == '/' ? "../../Data/Netlib/" : 
	 "..\\..\\Data\\Netlib\\";
   //#endif
   }

   {
      OsiSymSolverInterface symSi;
      symSi.setSymParam(OsiSymVerbosity, -1);
      testingMessage( "Now testing the OsiRowCut class with " );
      testingMessage( "OsiSymSolverInterface\n\n" );
      OsiRowCutUnitTest(&symSi,mpsDir);
   }
   {
      OsiSymSolverInterface symSi;
      symSi.setSymParam(OsiSymVerbosity, -1);
      testingMessage( "Now testing the OsiColCut class with " );
      testingMessage( "OsiSymSolverInterface\n\n" );
      OsiColCutUnitTest(&symSi,mpsDir);
   }

   {
      OsiSymSolverInterface symSi;
      symSi.setSymParam(OsiSymVerbosity, -1);
      testingMessage( "Now testing the OsiRowCutDebugger class with " );
      testingMessage( "OsiSymSolverInterface\n\n" );
      OsiRowCutDebuggerUnitTest(&symSi,mpsDir);
   }
   
   testingMessage( "Now testing OsiSymSolverInterface\n\n" );
   OsiSymSolverInterfaceUnitTest(mpsDir,netlibDir);
   
   
   if (parms.find("-testOsiSolverInterface") != parms.end()) {
      
      // Create vector of solver interfaces
      std::vector<OsiSolverInterface*> vecSi;
      
      OsiSolverInterface * symSi = new OsiSymSolverInterface;
      vecSi.push_back(symSi);
      
      testingMessage( "Testing OsiSolverInterface\n" );
      OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir);

      for (i=0; i<vecSi.size(); i++){
	 delete vecSi[i];
      }
            
   }     

   if (parms.find("-T") != parms.end()){
      testingMessage( "Testing MIPLIB files\n" );

      sym_environment *env = sym_open_environment();
      sym_parse_command_line(env, argc, const_cast<char**>(argv));
      sym_test(env, &test_status);
      if (test_status>0) {
         testingMessage( "warning: some instances may not have returned a ");
         testingMessage( "correct solution\n" );
      }
   }

   testingMessage( "All tests completed successfully\n" );
  
   return 0;
}

void testingMessage( const char * const msg )
{
  std::cout.flush() ;
  std::cerr <<msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}
