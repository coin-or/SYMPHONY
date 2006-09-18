// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
 
//#include "OsiConfig.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include "CoinTime.hpp"
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdio>

/*
  A utility definition which allows for easy suppression of unused variable
  warnings from GCC. Handy in this environment, where we're constantly def'ing
  things in and out.
*/
#ifndef UNUSED
# if defined(__GNUC__)
#   define UNUSED __attribute__((unused))
# else
#   define UNUSED
# endif
#endif

#include "OsiSolverInterface.hpp"
#ifdef COIN_HAS_VOL
#include "OsiVolSolverInterface.hpp"
#endif
#ifdef COIN_HAS_DYLP
#include "OsiDylpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_GLPK
#include "OsiGlpkSolverInterface.hpp"
#endif
#ifdef COIN_HAS_XPR
#include "OsiXprSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_SPX
#include "OsiSpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_OSL
#include "OsiOslSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_SYM
#include "OsiSymSolverInterface.hpp"
#endif
#ifdef COIN_HAS_FMP
#include "OsiFmpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CBC
#include "OsiCbcSolverInterface.hpp"
#endif
#include "CoinFloatEqual.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"
#include "OsiPresolve.hpp"
#ifdef COIN_OPBDP
#include "OsiOpbdpSolve.hpp"
#endif
//--------------------------------------------------------------------------
// A helper function to write out a message about a test failure
static void
failureMessage(
               const std::string & solverName,
               const std::string & message )
{
  std::string messageText;
  messageText = "*** ";
  messageText += solverName + "SolverInterface testing issue: ";
  messageText += message;
  std::cerr <<messageText.c_str() <<std::endl;
}

static void
failureMessage(
               const OsiSolverInterface & si,
               const std::string & message )
{
  std::string solverName;
  si.getStrParam(OsiSolverName,solverName);
  failureMessage(solverName,message);
}

//--------------------------------------------------------------------------
// A helper function to compare the equivalence of two vectors 
static bool
equivalentVectors(const OsiSolverInterface * si1,
		  const OsiSolverInterface * si2,
		  double tol,
		  const double * v1,
		  const double * v2,
		  int size)
{
  bool retVal = true;
  CoinRelFltEq eq(tol);
  int i;
  for ( i=0; i<size; i++ ) {
    
    // If both are equal to infinity then iterate
    if ( fabs(v1[i])==si1->getInfinity() && fabs(v2[i])==si2->getInfinity() )
       continue;
    
    // Test to see if equal
    if ( !eq(v1[i],v2[i]) ) {
      std::cerr <<"eq " <<i <<" " <<v1[i] <<" " <<v2[i] <<std::endl;
      //const OsiOslSolverInterface * oslSi = dynamic_cast<const OsiOslSolverInterface*>(si1);
      //assert(oslSi != NULL);
      //EKKModel * m = ((OsiOslSolverInterface*)oslSi)->modelPtr();
      //ekk_printSolution(m);
      retVal = false;
      break;
    }

  }
  return retVal;
}

//--------------------------------------------------------------------------
bool test1VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedMatrix m;

	m.transpose();

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	m.appendRow(r0);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	m.appendRow(r1);

	int numcol = 2;

	double *obj = new double[numcol];
	obj[0] = 3;
	obj[1] = 1;

	double *collb = new double[numcol];
	collb[0] = 0;
	collb[1] = 0;

	double *colub = new double[numcol];
	colub[0] = inf;
	colub[1] = inf;

	int numrow = 2;

	double *rowlb = new double[numrow];
	rowlb[0] = -inf;
	rowlb[1] = -inf;

	double *rowub = new double[numrow];
	rowub[0] = 10;
	rowub[1] = 15;

	s->loadProblem(m, collb, colub, obj, rowlb, rowub);

	delete [] obj;
	delete [] collb;
	delete [] colub;

	delete [] rowlb;
	delete [] rowub;

	s->setObjSense(-1);

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001, s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 5};
	ret = ret && equivalentVectors(s,s,0.0001, s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001, s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15};
	ret = ret && equivalentVectors(s,s,0.0001, s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------
bool test2VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedMatrix m;

	m.transpose();

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	m.appendRow(r0);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	m.appendRow(r1);

	CoinPackedVector r2;
	r2.insert(0, 1);
	r2.insert(1, 1);
	m.appendRow(r2);

	int numcol = 2;

	double *obj = new double[numcol];
	obj[0] = 3;
	obj[1] = 1;

	double *collb = new double[numcol];
	collb[0] = 0;
	collb[1] = 0;

	double *colub = new double[numcol];
	colub[0] = inf;
	colub[1] = inf;

	int numrow = 3;

	double *rowlb = new double[numrow];
	rowlb[0] = -inf;
	rowlb[1] = -inf;
	rowlb[2] = 1;

	double *rowub = new double[numrow];
	rowub[0] = 10;
	rowub[1] = 15;
	rowub[2] = inf;

	s->loadProblem(m, collb, colub, obj, rowlb, rowub);

	delete [] obj;
	delete [] collb;
	delete [] colub;

	delete [] rowlb;
	delete [] rowub;

	s->setObjSense(-1);

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001, s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 5, 5};
	ret = ret && equivalentVectors(s,s,0.0001, s->getRowActivity(), activity1, 3);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001, s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15, 7};
	ret = ret && equivalentVectors(s,s,0.0001, s->getRowActivity(), activity2, 3);

  return ret;
}


//--------------------------------------------------------------------------
bool test3VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	//double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, 0, 10, 3);
	s->addCol(empty, 0, 10, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, 0, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, 0, 15);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------
bool test4VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, 0, inf, 3);
	s->addCol(empty, 0, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, -inf, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, -inf, 15);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------

/*
  Constructs the system

     max    3x1 +  x2

    -inf <= 2x1 +  x2 <= 10
    -inf <=  x1 + 3x2 <= 15

  The optimal solution is unbounded. Objective is then changed to

     max     x1 +  x2

  which has a bounded optimum at x1 = 3, x2 = 4.
*/

bool test5VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, -inf, inf, 3);
	s->addCol(empty, -inf, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, -inf, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, -inf, 15);

	s->setObjSense(-1);

	s->writeMps("test");
        s->initialSolve();

	ret = ret && !s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && s->isProvenDualInfeasible();

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------


bool test6VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, 0, inf, 3);
	s->addCol(empty, 0, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, 0, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, 0, 15);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------


bool test7VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, 4, inf, 3);
	s->addCol(empty, 3, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, 0, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, 0, 15);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && !s->isProvenOptimal();
	ret = ret && s->isProvenPrimalInfeasible();

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && !s->isProvenOptimal();
	ret = ret && s->isProvenPrimalInfeasible();

	return ret;
}

//--------------------------------------------------------------------------

bool test8VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, -inf, inf, 3);
	s->addCol(empty, -inf, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, 0, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, 0, 15);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {6, -2};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------

bool test9VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, -inf, inf, 3);
	s->addCol(empty, -inf, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 2);
	r0.insert(1, 1);
	s->addRow(r0, 0, 10);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, 0, 15);

	CoinPackedVector r2;
	r2.insert(0, 1);
	r2.insert(1, 4);
	s->addRow(r2, 12, inf);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {4, 2};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {10, 10, 12};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 3);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {10, 15, 19};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 3);

	return ret;
}

//--------------------------------------------------------------------------

bool test10VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	int numcols = 2;
	int numrows = 2;
	const CoinBigIndex start[] = {0, 2, 4};
	const int index[] = {0, 1, 0, 1};
	const double value[] = {4, 1, 2, 3};
	const double collb[] = {0, 0};
	const double colub[] = {inf, inf};
	double obj[] = {3, 1};
	char rowsen[] = {'R', 'R'};
	double rowrhs[] = {20, 15};
	double rowrng[] = {20, 15};

	s->loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowsen, rowrhs, rowrng);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {20, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {20, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------
bool test11VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	int numcols = 2;
	int numrows = 2;
	const CoinBigIndex start[] = {0, 2, 4};
	const int index[] = {0, 1, 0, 1};
	const double value[] = {4, 1, 2, 3};
	const double collb[] = {0, 0};
	const double colub[] = {inf, inf};
	double obj[] = {3, 1};
	double rowlb[] = {0, 0};
	double rowub[] = {20, 15};

	s->loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {20, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {20, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------

bool test12VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedMatrix m;

	m.transpose();

	CoinPackedVector r0;
	r0.insert(0, 4);
	r0.insert(1, 2);
	m.appendRow(r0);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	m.appendRow(r1);

	int numcol = 2;

	double *obj = new double[numcol];
	obj[0] = 3;
	obj[1] = 1;

	double *collb = new double[numcol];
	collb[0] = 0;
	collb[1] = 0;

	double *colub = new double[numcol];
	colub[0] = inf;
	colub[1] = inf;

	int numrow = 2;

	double *rowlb = new double[numrow];
	rowlb[0] = 0;
	rowlb[1] = 0;

	double *rowub = new double[numrow];
	rowub[0] = 20;
	rowub[1] = 15;

	s->loadProblem(m, collb, colub, obj, rowlb, rowub);

	delete [] obj;
	delete [] collb;
	delete [] colub;

	delete [] rowlb;
	delete [] rowub;

	s->setObjSense(-1);

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {20, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {20, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------

bool test13VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedMatrix m;

	CoinPackedVector c0;
	c0.insert(0, 4);
	c0.insert(1, 1);
	m.appendCol(c0);

	CoinPackedVector c1;
	c1.insert(0, 2);
	c1.insert(1, 3);
	m.appendCol(c1);

	int numcol = 2;

	double *obj = new double[numcol];
	obj[0] = 3;
	obj[1] = 1;

	double *collb = new double[numcol];
	collb[0] = 0;
	collb[1] = 0;

	double *colub = new double[numcol];
	colub[0] = inf;
	colub[1] = inf;

	int numrow = 2;

	double *rowlb = new double[numrow];
	rowlb[0] = 0;
	rowlb[1] = 0;

	double *rowub = new double[numrow];
	rowub[0] = 20;
	rowub[1] = 15;

	s->loadProblem(m, collb, colub, obj, rowlb, rowub);

	delete [] obj;
	delete [] collb;
	delete [] colub;

	delete [] rowlb;
	delete [] rowub;

	s->setObjSense(-1);

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {20, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {20, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------

bool test14VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addCol(empty, 0, inf, 3);
	s->addCol(empty, 0, inf, 1);

	CoinPackedVector r0;
	r0.insert(0, 4);
	r0.insert(1, 2);
	s->addRow(r0, 0, 20);

	CoinPackedVector r1;
	r1.insert(0, 1);
	r1.insert(1, 3);
	s->addRow(r1, 0, 15);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {20, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {20, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------

bool test15VivianDeSmedt(OsiSolverInterface *s)
{
	bool ret = true;

	double inf = s->getInfinity();

	CoinPackedVector empty;

	s->addRow(empty, 0, 20);
	s->addRow(empty, 0, 15);

	CoinPackedVector c0;
	c0.insert(0, 4);
	c0.insert(1, 1);
	s->addCol(c0, 0, inf, 3);

	CoinPackedVector c1;
	c1.insert(0, 2);
	c1.insert(1, 3);
	s->addCol(c1, 0, inf, 1);

	s->setObjSense(-1);

	s->writeMps("test");

	s->initialSolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution1[] = {5, 0};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution1, 2);

	const double activity1[] = {20, 5};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity1, 2);

	s->setObjCoeff(0, 1);
	s->setObjCoeff(1, 1);

	s->resolve();

	ret = ret && s->isProvenOptimal();
	ret = ret && !s->isProvenPrimalInfeasible();
	ret = ret && !s->isProvenDualInfeasible();

	const double solution2[] = {3, 4};
	ret = ret && equivalentVectors(s,s,0.0001,s->getColSolution(), solution2, 2);

	const double activity2[] = {20, 15};
	ret = ret && equivalentVectors(s,s,0.0001,s->getRowActivity(), activity2, 2);

	return ret;
}

//--------------------------------------------------------------------------


/*! \brief Run solvers on NetLib problems.

  The routine creates a vector of NetLib problems (problem name, objective,
  various other characteristics), and a vector of solvers to be tested.
  
  Each solver is run on each problem. The run is deemed successful if the
  solver reports the correct problem size after loading and returns the
  correct objective value after optimization.

  If multiple solvers are available, the results are compared pairwise against
  the results reported by adjacent solvers in the solver vector. Due to
  limitations of the volume solver, it must be the last solver in vecEmptySiP.
*/

void OsiSolverInterfaceMpsUnitTest
  (const std::vector<OsiSolverInterface*> & vecEmptySiP,
   const std::string & mpsDir)

{ int i ;
  unsigned int m ;

/*
  Vectors to hold test problem names and characteristics. The objective value
  after optimization (objValue) must agree to the specified tolerance
  (objValueTol).
*/
  std::vector<std::string> mpsName ;
  std::vector<bool> min ;
  std::vector<int> nRows ;
  std::vector<int> nCols ;
  std::vector<double> objValue ;
  std::vector<double> objValueTol ;
/*
  And a macro to make the vector creation marginally readable.
*/
#define PUSH_MPS(zz_mpsName_zz,zz_min_zz,\
		 zz_nRows_zz,zz_nCols_zz,zz_objValue_zz,zz_objValueTol_zz) \
  mpsName.push_back(zz_mpsName_zz) ; \
  min.push_back(zz_min_zz) ; \
  nRows.push_back(zz_nRows_zz) ; \
  nCols.push_back(zz_nCols_zz) ; \
  objValueTol.push_back(zz_objValueTol_zz) ; \
  objValue.push_back(zz_objValue_zz) ;

/*
  Load up the problem vector. Note that the row counts here include the
  objective function.
*/
  PUSH_MPS("25fv47",true,822,1571,5.5018458883E+03,1.0e-10)
  PUSH_MPS("80bau3b",true,2263,9799,9.8722419241E+05,1.e-10)
  PUSH_MPS("adlittle",true,57,97,2.2549496316e+05,1.e-10)
  PUSH_MPS("afiro",true,28,32,-4.6475314286e+02,1.e-10)
  PUSH_MPS("agg",true,489,163,-3.5991767287e+07,1.e-10)
  PUSH_MPS("agg2",true,517,302,-2.0239252356e+07,1.e-10)
  PUSH_MPS("agg3",true,517,302,1.0312115935e+07,1.e-10)
  PUSH_MPS("bandm",true,306,472,-1.5862801845e+02,1.e-10)
  PUSH_MPS("beaconfd",true,174,262,3.3592485807e+04,1.e-10)
  PUSH_MPS("blend",true,75,83,-3.0812149846e+01,1.e-10)
  PUSH_MPS("bnl1",true,644,1175,1.9776295615E+03,1.e-10)
  PUSH_MPS("bnl2",true,2325,3489,1.8112365404e+03,1.e-10)
  PUSH_MPS("boeing1",true,/*351*/352,384,-3.3521356751e+02,1.e-10)
  PUSH_MPS("boeing2",true,167,143,-3.1501872802e+02,1.e-10)
  PUSH_MPS("bore3d",true,234,315,1.3730803942e+03,1.e-10)
  PUSH_MPS("brandy",true,221,249,1.5185098965e+03,1.e-10)
  PUSH_MPS("capri",true,272,353,2.6900129138e+03,1.e-10)
  PUSH_MPS("cycle",true,1904,2857,-5.2263930249e+00,1.e-9)
  PUSH_MPS("czprob",true,930,3523,2.1851966989e+06,1.e-10)
  PUSH_MPS("d2q06c",true,2172,5167,122784.21557456,1.e-7)
  PUSH_MPS("d6cube",true,416,6184,3.1549166667e+02,1.e-8)
  PUSH_MPS("degen2",true,445,534,-1.4351780000e+03,1.e-10)
  PUSH_MPS("degen3",true,1504,1818,-9.8729400000e+02,1.e-10)
  PUSH_MPS("dfl001",true,6072,12230,1.1266396047E+07,1.e-5)
  PUSH_MPS("e226",true,224,282,(-18.751929066+7.113),1.e-10) // NOTE: Objective function has constant of 7.113
  PUSH_MPS("etamacro",true,401,688,-7.5571521774e+02 ,1.e-6)
  PUSH_MPS("fffff800",true,525,854,5.5567961165e+05,1.e-6)
  PUSH_MPS("finnis",true,498,614,1.7279096547e+05,1.e-6)
  PUSH_MPS("fit1d",true,25,1026,-9.1463780924e+03,1.e-10)
  PUSH_MPS("fit1p",true,628,1677,9.1463780924e+03,1.e-10)
  PUSH_MPS("fit2d",true,26,10500,-6.8464293294e+04,1.e-10)
  PUSH_MPS("fit2p",true,3001,13525,6.8464293232e+04,1.e-9)
  PUSH_MPS("forplan",true,162,421,-6.6421873953e+02,1.e-6)
  PUSH_MPS("ganges",true,1310,1681,-1.0958636356e+05,1.e-5)
  PUSH_MPS("gfrd-pnc",true,617,1092,6.9022359995e+06,1.e-10)
  PUSH_MPS("greenbea",true,2393,5405,/*-7.2462405908e+07*/-72555248.129846,1.e-10)
  PUSH_MPS("greenbeb",true,2393,5405,/*-4.3021476065e+06*/-4302260.2612066,1.e-10)
  PUSH_MPS("grow15",true,301,645,-1.0687094129e+08,1.e-10)
  PUSH_MPS("grow22",true,441,946,-1.6083433648e+08,1.e-10)
  PUSH_MPS("grow7",true,141,301,-4.7787811815e+07,1.e-10)
  PUSH_MPS("israel",true,175,142,-8.9664482186e+05,1.e-10)
  PUSH_MPS("kb2",true,44,41,-1.7499001299e+03,1.e-10)
  PUSH_MPS("lotfi",true,154,308,-2.5264706062e+01,1.e-10)
  PUSH_MPS("maros",true,847,1443,-5.8063743701e+04,1.e-10)
  PUSH_MPS("maros-r7",true,3137,9408,1.4971851665e+06,1.e-10)
  PUSH_MPS("modszk1",true,688,1620,3.2061972906e+02,1.e-10)
  PUSH_MPS("nesm",true,663,2923,1.4076073035e+07,1.e-5)
  PUSH_MPS("perold",true,626,1376,-9.3807580773e+03,1.e-6)
  PUSH_MPS("pilot",true,1442,3652,/*-5.5740430007e+02*/-557.48972927292,5.e-5)
  PUSH_MPS("pilot4",true,411,1000,-2.5811392641e+03,1.e-6)
  PUSH_MPS("pilot87",true,2031,4883,3.0171072827e+02,1.e-4)
  PUSH_MPS("pilotnov",true,976,2172,-4.4972761882e+03,1.e-10)
  // ?? PUSH_MPS("qap8",true,913,1632,2.0350000000e+02,1.e-10)
  // ?? PUSH_MPS("qap12",true,3193,8856,5.2289435056e+02,1.e-10)
  // ?? PUSH_MPS("qap15",true,6331,22275,1.0409940410e+03,1.e-10)
  PUSH_MPS("recipe",true,92,180,-2.6661600000e+02,1.e-10)
  PUSH_MPS("sc105",true,106,103,-5.2202061212e+01,1.e-10)
  PUSH_MPS("sc205",true,206,203,-5.2202061212e+01,1.e-10)
  PUSH_MPS("sc50a",true,51,48,-6.4575077059e+01,1.e-10)
  PUSH_MPS("sc50b",true,51,48,-7.0000000000e+01,1.e-10)
  PUSH_MPS("scagr25",true,472,500,-1.4753433061e+07,1.e-10)
  PUSH_MPS("scagr7",true,130,140,-2.3313892548e+06,1.e-6)
  PUSH_MPS("scfxm1",true,331,457,1.8416759028e+04,1.e-10)
  PUSH_MPS("scfxm2",true,661,914,3.6660261565e+04,1.e-10)
  PUSH_MPS("scfxm3",true,991,1371,5.4901254550e+04,1.e-10)
  PUSH_MPS("scorpion",true,389,358,1.8781248227e+03,1.e-10)
  PUSH_MPS("scrs8",true,491,1169,9.0429998619e+02,1.e-5)
  PUSH_MPS("scsd1",true,78,760,8.6666666743e+00,1.e-10)
  PUSH_MPS("scsd6",true,148,1350,5.0500000078e+01,1.e-10)
  PUSH_MPS("scsd8",true,398,2750,9.0499999993e+02,1.e-8)
  PUSH_MPS("sctap1",true,301,480,1.4122500000e+03,1.e-10)
  PUSH_MPS("sctap2",true,1091,1880,1.7248071429e+03,1.e-10)
  PUSH_MPS("sctap3",true,1481,2480,1.4240000000e+03,1.e-10)
  PUSH_MPS("seba",true,516,1028,1.5711600000e+04,1.e-10)
  PUSH_MPS("share1b",true,118,225,-7.6589318579e+04,1.e-10)
  PUSH_MPS("share2b",true,97,79,-4.1573224074e+02,1.e-10)
  PUSH_MPS("shell",true,537,1775,1.2088253460e+09,1.e-10)
  PUSH_MPS("ship04l",true,403,2118,1.7933245380e+06,1.e-10)
  PUSH_MPS("ship04s",true,403,1458,1.7987147004e+06,1.e-10)
  PUSH_MPS("ship08l",true,779,4283,1.9090552114e+06,1.e-10)
  PUSH_MPS("ship08s",true,779,2387,1.9200982105e+06,1.e-10)
  PUSH_MPS("ship12l",true,1152,5427,1.4701879193e+06,1.e-10)
  PUSH_MPS("ship12s",true,1152,2763,1.4892361344e+06,1.e-10)
  PUSH_MPS("sierra",true,1228,2036,1.5394362184e+07,1.e-10)
  PUSH_MPS("stair",true,357,467,-2.5126695119e+02,1.e-10)
  PUSH_MPS("standata",true,360,1075,1.2576995000e+03,1.e-10)
  // GUB PUSH_MPS("standgub",true,362,1184,1257.6995,1.e-10) 
  PUSH_MPS("standmps",true,468,1075,1.4060175000E+03,1.e-10) 
  PUSH_MPS("stocfor1",true,118,111,-4.1131976219E+04,1.e-10)
  PUSH_MPS("stocfor2",true,2158,2031,-3.9024408538e+04,1.e-10)
  // ?? PUSH_MPS("stocfor3",true,16676,15695,-3.9976661576e+04,1.e-10)
  // ?? PUSH_MPS("truss",true,1001,8806,4.5881584719e+05,1.e-10)
  PUSH_MPS("tuff",true,334,587,2.9214776509e-01,1.e-10)
  PUSH_MPS("vtpbase",true,199,203,1.2983146246e+05,1.e-10)
  PUSH_MPS("wood1p",true,245,2594,1.4429024116e+00,5.e-5)
  PUSH_MPS("woodw",true,1099,8405,1.3044763331E+00,1.e-10)

#undef PUSH_MPS

/*
  Create a vector of solver interfaces that we can use to run the test
  problems. The strategy is to create a fresh clone of the `empty' solvers
  from vecEmptySiP for each problem, then proceed in stages: read the MPS
  file, solve the problem, check the solution. If there are multiple
  solvers in vecSiP, the results of each solver are compared with its
  neighbors in the vector.
*/
  std::vector<OsiSolverInterface*> vecSiP(vecEmptySiP.size()) ;

  // Create vector to store a name for each solver interface
  // and a count on the number of problems the solver intface solved.
  std::vector<std::string> siName;
  std::vector<int> numProbSolved;
  std::vector<double> timeTaken;
  const int vecsize = vecSiP.size();
  for ( i=0; i<vecsize; i++ ) {
    siName.push_back("unknown");
    numProbSolved.push_back(0);
    timeTaken.push_back(0.0);
  }
  
  
  //Open the main loop to step through the MPS problems.
  for (m = 0 ; m < mpsName.size() ; m++) { 
    std::cerr << "  processing mps file: " << mpsName[m] 
      << " (" << m+1 << " out of " << mpsName.size() << ")" << std::endl ;
    bool allSolversReadMpsFile = true;
    
    
    //Stage 1: Read the MPS file into each solver interface.
    //Fill vecSiP with fresh clones of the solvers and read in the MPS file. As
    //a basic check, make sure the size of the constraint matrix is correct.
    for (i = vecSiP.size()-1 ; i >= 0 ; --i) { 
      vecSiP[i] = vecEmptySiP[i]->clone() ;
      
      vecSiP[i]->getStrParam(OsiSolverName,siName[i]);
      
      std::string fn = mpsDir+mpsName[m] ;
      vecSiP[i]->readMps(fn.c_str(),"mps") ;
      
      if (min[m])
        vecSiP[i]->setObjSense(1.0) ;
      else
        vecSiP[i]->setObjSense(-1.0) ;
      
      int nr = vecSiP[i]->getNumRows() ;
      int nc = vecSiP[i]->getNumCols() ;
      assert(nr == nRows[m]-1) ;
      assert(nc == nCols[m]) ; 
    } 
    
    //If we have multiple solvers, compare the representations.
    if ( allSolversReadMpsFile )
      for (i = vecSiP.size()-1 ; i > 0 ; --i) { 
        CoinPackedVector vim1,vi ;
        
        // Compare col lowerbounds
        assert(
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getColLower(),vecSiP[i  ]->getColLower(),
          vecSiP[i  ]->getNumCols() )
          ) ;
        
        // Compare col upperbounds
        assert(
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getColUpper(),vecSiP[i  ]->getColUpper(),
          vecSiP[i  ]->getNumCols() )
          ) ;
        
        // Compare row lowerbounds
        assert(
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getRowLower(),vecSiP[i  ]->getRowLower(),
          vecSiP[i  ]->getNumRows() )
          ) ;
        
        // Compare row upperbounds
        assert( 
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getRowUpper(),vecSiP[i  ]->getRowUpper(),
          vecSiP[i  ]->getNumRows() )
          ) ;
        
        // Compare row sense
        {
          const char * rsm1 = vecSiP[i-1]->getRowSense() ;
          const char * rs   = vecSiP[i  ]->getRowSense() ;
          int nr = vecSiP[i]->getNumRows() ;
          int r ;
          for (r = 0 ; r < nr ; r++) assert (rsm1[r] == rs[r]) ;
        }
        
        // Compare rhs
        assert( 
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getRightHandSide(),vecSiP[i  ]->getRightHandSide(),
          vecSiP[i  ]->getNumRows() )
          ) ;
        
        // Compare range
        assert( 
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getRowRange(),vecSiP[i  ]->getRowRange(),
          vecSiP[i  ]->getNumRows() )
          ) ;
        
        // Compare objective sense
        assert( vecSiP[i-1]->getObjSense() == vecSiP[i  ]->getObjSense() ) ;
        
        // Compare objective coefficients
        assert( 
          equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
          vecSiP[i-1]->getObjCoefficients(),vecSiP[i  ]->getObjCoefficients(),
          vecSiP[i  ]->getNumCols() )
          ) ;
        
        // Compare number of elements
        assert( vecSiP[i-1]->getNumElements() == vecSiP[i]->getNumElements() ) ;
        
        // Compare constraint matrix
        { 
          const CoinPackedMatrix * rmm1=vecSiP[i-1]->getMatrixByRow() ;
          const CoinPackedMatrix * rm  =vecSiP[i  ]->getMatrixByRow() ;
          assert( rmm1->isEquivalent(*rm) ) ;
          
          const CoinPackedMatrix * cmm1=vecSiP[i-1]->getMatrixByCol() ;
          const CoinPackedMatrix * cm  =vecSiP[i  ]->getMatrixByCol() ;
          assert( cmm1->isEquivalent(*cm) ) ; 
        } 
      }
      
      //If we have multiple solvers, compare the variable type information      
      if ( allSolversReadMpsFile )
        for (i = vecSiP.size()-1 ; i > 0 ; --i){ 
          CoinPackedVector vim1,vi ;
          int c ;
          
          { 
            OsiVectorInt sm1 = vecSiP[i-1]->getFractionalIndices() ;
            OsiVectorInt s   = vecSiP[i  ]->getFractionalIndices() ;
            assert( sm1.size() == s.size() ) ;
            for (c = s.size()-1 ; c >= 0 ; --c) assert( sm1[c] == s[c] ) ; 
          }
          
          { 
            int nc = vecSiP[i]->getNumCols() ;
            for (c = 0 ; c < nc ; c++){ 
              assert(
                vecSiP[i-1]->isContinuous(c) == vecSiP[i]->isContinuous(c)
                ) ;
              assert(
                vecSiP[i-1]->isBinary(c) == vecSiP[i]->isBinary(c)
                ) ;
              assert(
                vecSiP[i-1]->isIntegerNonBinary(c) ==
                vecSiP[i  ]->isIntegerNonBinary(c)
                ) ;
              assert(
                vecSiP[i-1]->isFreeBinary(c) == vecSiP[i]->isFreeBinary(c)
                ) ;
              assert(
                vecSiP[i-1]->isInteger(c) == vecSiP[i]->isInteger(c)
                ) ; 
            } 
          } 
        }
        
      //Stage 2: Call each solver to solve the problem.
      //
      // We call each solver, then check the return code and objective.
      //  
      //    Note that the volume solver can't handle the Netlib cases. The strategy is
      //    to require that it be the last solver in vecSiP and then break out of the
      //    loop. This ensures that all previous solvers are run and compared to one
      //    another.      
      for (i = 0 ; i < static_cast<int>(vecSiP.size()) ; ++i) {
        double startTime = CoinCpuTime();
        
#     ifdef COIN_HAS_VOL
        { 
          OsiVolSolverInterface * si =
            dynamic_cast<OsiVolSolverInterface *>(vecSiP[i]) ;
          if (si != NULL )  { 
            // VOL does not solve netlib cases so don't bother trying to solve
            break ; 
          }
        }
#     endif

	try
	{ vecSiP[i]->initialSolve() ; }
	catch (CoinError &thrownErr)
	{ std::cerr << thrownErr.className() << "::" << thrownErr.methodName()
		    << ": " << thrownErr.message() << std::endl ; }
        double timeOfSolution = CoinCpuTime()-startTime;
        if (vecSiP[i]->isProvenOptimal()) { 
          double soln = vecSiP[i]->getObjValue();       
          CoinRelFltEq eq(objValueTol[m]) ;
          if (eq(soln,objValue[m])) { 
            std::cerr 
              <<siName[i]<<"SolverInterface "
              << soln << " = " << objValue[m] <<", "
	      << vecSiP[i]->getIterationCount() << " iters"
	      << "; okay";
            numProbSolved[i]++;
          } else  { 
            std::cerr <<siName[i] <<" " <<soln << " != " <<objValue[m] << "; error=" ;
            std::cerr <<fabs(objValue[m] - soln); 
          }
        } else {
	   if (vecSiP[i]->isProvenPrimalInfeasible()) 
	      std::cerr << "error; primal infeasible" ;
#ifndef COIN_HAS_SYM 
	   else if (vecSiP[i]->isProvenDualInfeasible())
	      std::cerr << "error; dual infeasible" ;
#endif
	   else if (vecSiP[i]->isIterationLimitReached()) 
	      std::cerr << "error; iteration limit" ;
	   else if (vecSiP[i]->isAbandoned()) 
	      std::cerr << "error; abandoned" ;
	   else  
	      std::cerr << "error; unknown" ;
        }
        std::cerr<<" - took " <<timeOfSolution<<" seconds."<<std::endl; 
        timeTaken[i] += timeOfSolution;
      }
      /*
      Delete the used solver interfaces so we can reload fresh clones for the
      next problem.
      */
      for (i = vecSiP.size()-1 ; i >= 0 ; --i) delete vecSiP[i] ;
  }

  const int siName_size = siName.size();
  for ( i=0; i<siName_size; i++ ) {
    std::cerr 
      <<siName[i] 
      <<" solved " 
      <<numProbSolved[i]
      <<" out of "
      <<objValue.size()
      <<" and took "
      <<timeTaken[i]
      <<" seconds."
      <<std::endl;
  } 
}


//#############################################################################
//#############################################################################

static bool testIntParam(OsiSolverInterface * si, int k, int val)
{
  int i = 123456789, orig = 123456789;
  bool ret;
  OsiIntParam key = static_cast<OsiIntParam>(k);
  si->getIntParam(key, orig);
  if (si->setIntParam(key, val)) {
    ret = (si->getIntParam(key, i) == true) && (i == val);
  } else {
    ret = (si->getIntParam(key, i) == true) && (i == orig);
  }
  return ret;
}
static bool testDblParam(OsiSolverInterface * si, int k, double val)
{
  double d = 123456789.0, orig = 123456789.0;
  bool ret;
  OsiDblParam key = static_cast<OsiDblParam>(k);
  si->getDblParam(key, orig);
  if (si->setDblParam(key, val)) {
    ret = (si->getDblParam(key, d) == true) && (d == val);
  } else {
    ret = (si->getDblParam(key, d) == true) && (d == orig);
  }
  return ret;
}

/*
  Original model

static bool testHintParam(OsiSolverInterface * si, int k, bool val,
			  OsiHintStrength strength)
{
  bool i = true, orig = true;
  OsiHintStrength i2 = OsiHintDo, orig2 = OsiHintDo;
  bool ret;
  OsiHintParam key = static_cast<OsiHintParam>(k);
  si->getHintParam(key, orig, orig2);
  if (si->setHintParam(key, val, strength)) {
    ret = (si->getHintParam(key, i, i2) == true) && (i2 == strength);
  } else {
    ret = (si->getHintParam(key, i, i2) == true) && (i2 == orig2);
  }
  return ret;
}
*/

static bool testHintParam(OsiSolverInterface * si, int k, bool sense,
			  OsiHintStrength strength, int *throws)
/*
  Tests for proper behaviour of [set,get]HintParam methods. The initial get
  tests the return value to see if the hint is implemented; the values
  returned for sense and strength are not checked.
  
  If the hint is implemented, a pair of set/get calls is performed at the
  strength specified by the parameter. The set can return true or, at
  strength OsiForceDo, throw an exception if the solver cannot comply. The
  rationale would be that only OsiForceDo must be obeyed, so anything else
  should return true regardless of whether the solver followed the hint.

  The test checks that the value and strength returned by getHintParam matches
  the previous call to setHintParam. This is arguably wrong --- one can argue
  that it should reflect the solver's ability to comply with the hint. But
  that's how the OSI interface standard has evolved up to now.

  If the hint is not implemented, attempting to set the hint should return
  false, or throw an exception at strength OsiForceDo.

  The testing code which calls testHintParam is set up so that a successful
  return is defined as true if the hint is implemented, false if it is not.
  Information printing is suppressed; uncomment and recompile if you want it.
*/
{ bool post_sense ;
  OsiHintStrength post_strength ;
  bool ret ;
  OsiHintParam key = static_cast<OsiHintParam>(k) ;

  if (si->getHintParam(key,post_sense,post_strength))
  { ret = false ;
    try
    { if (si->setHintParam(key,sense,strength))
      { ret = (si->getHintParam(key,post_sense,post_strength) == true) &&
	      (post_strength == strength) && (post_sense == sense) ; } }
    catch (CoinError &thrownErr)
    { // std::ostringstream msg ;
      // msg << "setHintParam throw for hint " << key << " sense " << sense <<
      //      " strength " << strength ;
      // failureMessage(*si,msg.str()) ;
      // std::cerr << thrownErr.className() << "::" << thrownErr.methodName() <<
      //	": " << thrownErr.message() << std::endl ;
      (*throws)++ ;
      ret = (strength == OsiForceDo) ; } }
  else
  { ret = true ;
    try
    { ret = si->setHintParam(key,sense,strength) ; }
    catch (CoinError &thrownErr)
    { // std::ostringstream msg ;
      // msg << "setHintParam throw for hint " << key << " sense " << sense <<
      //      " strength " << strength ;
      // failureMessage(*si,msg.str()) ;
      // std::cerr << thrownErr.className() << "::" << thrownErr.methodName() <<
      //	": " << thrownErr.message() << std::endl ;
      (*throws)++ ;
      ret = !(strength == OsiForceDo) ; } }
  
  return ret ; }

//#############################################################################
//#############################################################################

CoinPackedMatrix &BuildExmip1Mtx ()
/*
  Simple function to build a packed matrix for the exmip1 example used in
  tests. The function exists solely to hide the intermediate variables.
  Probably could be written as an initialised declaration.
  See COIN/Mps/Sample/exmip1.mps for a human-readable presentation.

  Ordered triples seem easiest. They're listed in row-major order.
*/

{ int rowndxs[] = { 0, 0, 0, 0, 0,
		    1, 1,
		    2, 2,
		    3, 3,
		    4, 4, 4 } ;
  int colndxs[] = { 0, 1, 3, 4, 7,
		    1, 2,
		    2, 5,
		    3, 6,
		    0, 4, 7 } ;
  double coeffs[] = { 3.0, 1.0, -2.0, -1.0, -1.0,
		      2.0, 1.1,
		      1.0, 1.0,
		      2.8, -1.2,
		      5.6, 1.0, 1.9 } ;

  static CoinPackedMatrix exmip1mtx =
    CoinPackedMatrix(true,&rowndxs[0],&colndxs[0],&coeffs[0],14) ;

  return (exmip1mtx) ; }


void
OsiSolverInterfaceCommonUnitTest(const OsiSolverInterface* emptySi,
				 const std::string & mpsDir, 
         const std::string & netlibDir)
{
  
  int i;
  CoinRelFltEq eq;

  std::string fn = mpsDir+"exmip1";
  OsiSolverInterface * exmip1Si = emptySi->clone(); 
  exmip1Si->readMps(fn.c_str(),"mps");

  // Test that solverInterface knows its name.
  // The name is used for displaying messages when testing
  std::string solverName;
  {
    OsiSolverInterface * si = emptySi->clone();
    bool supportsSolverName = si->getStrParam(OsiSolverName,solverName);
    assert( supportsSolverName );
    assert( solverName != "Unknown Solver" );
    delete si;
  }

  // Determine if this emptySi is an OsiVolSolverInterface
  bool volSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_VOL
    const OsiVolSolverInterface * si =
      dynamic_cast<const OsiVolSolverInterface *>(emptySi);
    if ( si != NULL ) volSolverInterface = true;
#endif
  }

  // Determine if this emptySi is an OsiOslSolverInterface
  bool oslSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_OSL
    const OsiOslSolverInterface * si =
      dynamic_cast<const OsiOslSolverInterface *>(emptySi);
    if ( si != NULL ) oslSolverInterface = true;
#endif
  }

  // Determine if this emptySi is an OsiDylpSolverInterface
  bool dylpSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_DYLP
    const OsiDylpSolverInterface * si =
      dynamic_cast<const OsiDylpSolverInterface *>(emptySi);
    if ( si != NULL ) dylpSolverInterface = true;
#endif
  }

  // Determine if this emptySi is an OsiGlpkSolverInterface
  bool glpkSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_GLPK
    const OsiGlpkSolverInterface * si =
      dynamic_cast<const OsiGlpkSolverInterface *>(emptySi);
    if ( si != NULL ) glpkSolverInterface = true;
#endif
  }

  // Determine if this emptySi is an OsiFmpSolverInterface
  bool fmpSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_FMP
    const OsiFmpSolverInterface * si =
      dynamic_cast<const OsiFmpSolverInterface *>(emptySi);
    if ( si != NULL ) fmpSolverInterface = true;
#endif
  }

  // Determine if this emptySi is an OsiXprSolverInterface
  bool xprSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_XPR
    const OsiXprSolverInterface * si =
      dynamic_cast<const OsiXprSolverInterface *>(emptySi);
    if ( si != NULL ) xprSolverInterface = true;
#endif
  }

  // Determine if this is the emptySi is an OsiSymSolverInterface
  bool symSolverInterface UNUSED = false;
  {
#ifdef COIN_HAS_SYM
     const OsiSymSolverInterface * si =
	dynamic_cast<const OsiSymSolverInterface *>(emptySi);
     if ( si != NULL ) symSolverInterface = true;
#endif
  }

  // Test that solverInterface knows about constants
  // in objective function.
  // Do not perform test if Vol solver, because it
  // requires problems of a special form and can not
  // solve netlib e226.
  if ( !volSolverInterface ) {
    OsiSolverInterface * si = emptySi->clone();
    std::string fn = netlibDir+"e226";
    int mpsRc = si->readMps(fn.c_str(),"mps");
    assert(mpsRc==0);
    si->initialSolve();
    double objValue = si->getObjValue(); 
    if( !eq(objValue,-18.751929066+7.113) )
      failureMessage(solverName,"getObjValue with constant in objective function");
    delete si;
  }

  // Test that values returned from an empty solverInterface
  {
    OsiSolverInterface * si = emptySi->clone();
    if( si->getNumRows()!=0 )
      failureMessage(solverName,"getNumRows with empty solverInterface");
    if( si->getNumCols()!=0 )
      failureMessage(solverName,"getNumCols with empty solverInterface");
    if( si->getNumElements()!=0 )
      failureMessage(solverName,"getNumElements with empty solverInterface");
    if( si->getColLower()!=NULL )
      failureMessage(solverName,"getColLower with empty solverInterface");
    if( si->getColUpper()!=NULL )
      failureMessage(solverName,"getColUpper with empty solverInterface");
    if( si->getColSolution()!=NULL )
      failureMessage(solverName,"getColSolution with empty solverInterface");
    if( si->getObjCoefficients()!=NULL )
      failureMessage(solverName,"getObjCoefficients with empty solverInterface");
    if( si->getRowRange()!=NULL )
      failureMessage(solverName,"getRowRange with empty solverInterface");
    if( si->getRightHandSide()!=NULL )
      failureMessage(solverName,"getRightHandSide with empty solverInterface");
    if( si->getRowSense()!=NULL )
      failureMessage(solverName,"getRowSense with empty solverInterface");
    if( si->getRowLower()!=NULL )
      failureMessage(solverName,"getRowLower with empty solverInterface");
    if( si->getRowUpper()!=NULL )
      failureMessage(solverName,"getRowUpper with empty solverInterface");
    delete si;
  }
  
  
  // Test that problem was loaded correctly

  { const char   * exmip1Sirs  = exmip1Si->getRowSense();

    assert( exmip1Sirs[0]=='G' );
    assert( exmip1Sirs[1]=='L' );
    assert( exmip1Sirs[2]=='E' );
    assert( exmip1Sirs[3]=='R' );
    assert( exmip1Sirs[4]=='R' );
    
    const double * exmip1Sirhs = exmip1Si->getRightHandSide();
    assert( eq(exmip1Sirhs[0],2.5) );
    assert( eq(exmip1Sirhs[1],2.1) );
    assert( eq(exmip1Sirhs[2],4.0) );
    assert( eq(exmip1Sirhs[3],5.0) );
    assert( eq(exmip1Sirhs[4],15.) ); 
    
    const double * exmip1Sirr  = exmip1Si->getRowRange();
    assert( eq(exmip1Sirr[0],0.0) );
    assert( eq(exmip1Sirr[1],0.0) );
    assert( eq(exmip1Sirr[2],0.0) );
    assert( eq(exmip1Sirr[3],5.0-1.8) );
    assert( eq(exmip1Sirr[4],15.0-3.0) );

    CoinPackedMatrix goldmtx ;
    goldmtx.reverseOrderedCopyOf(BuildExmip1Mtx()) ;
    CoinPackedMatrix pm;
    pm.setExtraGap(0.0);
    pm.setExtraMajor(0.0);
    pm = *exmip1Si->getMatrixByRow();
    pm.removeGaps();
    assert(goldmtx.isEquivalent(pm)) ;
    
    int nc = exmip1Si->getNumCols();
    int nr = exmip1Si->getNumRows();
    const double * cl = exmip1Si->getColLower();
    const double * cu = exmip1Si->getColUpper();
    const double * rl = exmip1Si->getRowLower();
    const double * ru = exmip1Si->getRowUpper();
    assert( nc == 8 );
    assert( nr == 5 );
    assert( eq(cl[0],2.5) );
    assert( eq(cl[1],0.0) );
    assert( eq(cl[2],0.0) );
    assert( eq(cl[3],0.0) );
    assert( eq(cl[4],0.5) );
    assert( eq(cl[5],0.0) );
    assert( eq(cl[6],0.0) );
    assert( eq(cl[7],0.0) );
    //    assert( eq(cu[0],exmip1Si->getInfinity()) );
    assert( eq(cu[1],4.1) );
    assert( eq(cu[2],1.0) );
    assert( eq(cu[3],1.0) );
    assert( eq(cu[4],4.0) );
    //   assert( eq(cu[5],exmip1Si->getInfinity()) );
    // assert( eq(cu[6],exmip1Si->getInfinity()) );
    assert( eq(cu[7],4.3) );

    assert( eq(rl[0],2.5) );
    assert( eq(rl[1],-exmip1Si->getInfinity()) );
    assert( eq(rl[2],4.0) );
    assert( eq(rl[3],1.8) );
    assert( eq(rl[4],3.0) );
    assert( eq(ru[0],exmip1Si->getInfinity()) );
    assert( eq(ru[1],2.1) );
    assert( eq(ru[2],4.0) );
    assert( eq(ru[3],5.0) );
    assert( eq(ru[4],15.0) );
    
    // make sure col solution is something reasonable,
    // that is between upper and lower bounds
    const double * cs = exmip1Si->getColSolution();
    int c;
    bool okColSol=true;
    //double inf = exmip1Si->getInfinity();
    for ( c=0; c<nc; c++ ) {
      // if colSol is not between column bounds then 
      // colSol is unreasonable.
      if( !(cl[c]<=cs[c] && cs[c]<=cu[c]) ) okColSol=false;
      // if at least one column bound is not infinite,
      // then it is unreasonable to have colSol as infinite
      // FIXME: temporarily commented out pending some group thought on the
      //	semantics of this test. -- lh, 03.04.29 --
      // if ( (cl[c]<inf || cu[c]<inf) && cs[c]>=inf ) okColSol=false;
    }
    if( !okColSol )
      failureMessage(solverName,"getColSolution before solve");
    
    // Test value of objective function coefficients
    const double * objCoef = exmip1Si->getObjCoefficients();
    assert( eq( objCoef[0],  1.0) );
    assert( eq( objCoef[1],  0.0) );
    assert( eq( objCoef[2],  0.0) );
    assert( eq( objCoef[3],  0.0) );
    assert( eq( objCoef[4],  2.0) );
    assert( eq( objCoef[5],  0.0) );
    assert( eq( objCoef[6],  0.0) );
    assert( eq( objCoef[7], -1.0) );

    // Test that objective value is correct
    double correctObjValue = CoinPackedVector(nc,objCoef).dotProduct(cs);
    double siObjValue = exmip1Si->getObjValue();
    if( !eq(correctObjValue,siObjValue) ) {
       // FIXME: the test checks the primal value. vol fails this, because vol
       // considers the dual value to be the objective value
       failureMessage(solverName,"getObjValue before solve (OK for vol)");
    }


  }
  
  
  // Test matrixByCol method
  {
    CoinPackedMatrix &goldmtx = BuildExmip1Mtx() ;
    OsiSolverInterface & si = *exmip1Si->clone();
    CoinPackedMatrix sm = *si.getMatrixByCol();
    sm.removeGaps();
    bool getByColOK = goldmtx.isEquivalent(sm) ;

    if (!getByColOK)
      failureMessage(solverName,"getMatrixByCol()") ;
    
    // Test getting and setting of objective offset
    double objOffset;
    bool supportOsiObjOffset = si.getDblParam(OsiObjOffset,objOffset);
    assert( supportOsiObjOffset );
    assert( eq( objOffset, 0.0 ) );
    supportOsiObjOffset = si.setDblParam(OsiObjOffset, 3.21);
    assert( supportOsiObjOffset );
    si.getDblParam(OsiObjOffset,objOffset);
    assert( eq( objOffset, 3.21 ) );
    
    delete &si;
  }
   
  // Test clone
  {
    OsiSolverInterface * si2;  
    int ad = 13579;
    {
      OsiSolverInterface * si1 = exmip1Si->clone(); 
      int ad = 13579;
      si1->setApplicationData(&ad);
      assert( *((int *)(si1->getApplicationData())) == ad );
      si2 = si1->clone();
      delete si1;
    }
    
    if( *((int *)(si2->getApplicationData())) != ad )
      failureMessage(solverName,"getApplicationData on cloned solverInterface");

    const char   * exmip1Sirs  = si2->getRowSense();
    assert( exmip1Sirs[0]=='G' );
    assert( exmip1Sirs[1]=='L' );
    assert( exmip1Sirs[2]=='E' );
    assert( exmip1Sirs[3]=='R' );
    assert( exmip1Sirs[4]=='R' );
    
    const double * exmip1Sirhs = si2->getRightHandSide();
    assert( eq(exmip1Sirhs[0],2.5) );
    assert( eq(exmip1Sirhs[1],2.1) );
    assert( eq(exmip1Sirhs[2],4.0) );
    assert( eq(exmip1Sirhs[3],5.0) );
    assert( eq(exmip1Sirhs[4],15.) ); 
    
    const double * exmip1Sirr  = si2->getRowRange();
    assert( eq(exmip1Sirr[0],0.0) );
    assert( eq(exmip1Sirr[1],0.0) );
    assert( eq(exmip1Sirr[2],0.0) );
    assert( eq(exmip1Sirr[3],5.0-1.8) );
    assert( eq(exmip1Sirr[4],15.0-3.0) );
    
    CoinPackedMatrix goldmtx ;
    goldmtx.reverseOrderedCopyOf(BuildExmip1Mtx()) ;
    CoinPackedMatrix pm;
    pm.setExtraGap(0.0);
    pm.setExtraMajor(0.0);
    pm = *si2->getMatrixByRow();
    assert(goldmtx.isEquivalent(pm)) ;
    
    int nc = si2->getNumCols();
    int nr = si2->getNumRows();
    const double * cl = si2->getColLower();
    const double * cu = si2->getColUpper();
    const double * rl = si2->getRowLower();
    const double * ru = si2->getRowUpper();
    assert( nc == 8 );
    assert( nr == 5 );
    assert( eq(cl[0],2.5) );
    assert( eq(cl[1],0.0) );
    assert( eq(cl[2],0.0) );
    assert( eq(cl[3],0.0) );
    assert( eq(cl[4],0.5) );
    assert( eq(cl[5],0.0) );
    assert( eq(cl[6],0.0) );
    assert( eq(cl[7],0.0) );
    //    assert( eq(cu[0],si2->getInfinity()) );
    assert( eq(cu[1],4.1) );
    assert( eq(cu[2],1.0) );
    assert( eq(cu[3],1.0) );
    assert( eq(cu[4],4.0) );
    //  assert( eq(cu[5],si2->getInfinity()) );
    // assert( eq(cu[6],si2->getInfinity()) );
    assert( eq(cu[7],4.3) );

    assert( eq(rl[0],2.5) );
    assert( eq(rl[1],-si2->getInfinity()) );
    assert( eq(rl[2],4.0) );
    assert( eq(rl[3],1.8) );
    assert( eq(rl[4],3.0) );
    assert( eq(ru[0],si2->getInfinity()) );
    assert( eq(ru[1],2.1) );
    assert( eq(ru[2],4.0) );
    assert( eq(ru[3],5.0) );
    assert( eq(ru[4],15.0) );
        
    // make sure col solution is something reasonable,
    // that is between upper and lower bounds
    const double * cs = exmip1Si->getColSolution();
    int c;
    bool okColSol=true;
    //double inf = exmip1Si->getInfinity();
    for ( c=0; c<nc; c++ ) {
      // if colSol is not between column bounds then 
      // colSol is unreasonable.
      if( !(cl[c]<=cs[c] && cs[c]<=cu[c]) ) okColSol=false;
      // if at least one column bound is not infinite,
      // then it is unreasonable to have colSol as infinite
      // FIXME: temporarily commented out pending some group thought on the
      //	semantics of this test. -- lh, 03.04.29 --
      // if ( (cl[c]<inf || cu[c]<inf) && cs[c]>=inf ) okColSol=false;
    }
    if( !okColSol )
      failureMessage(solverName,"getColSolution before solve on cloned solverInterface");
    
    assert( eq( si2->getObjCoefficients()[0],  1.0) );
    assert( eq( si2->getObjCoefficients()[1],  0.0) );
    assert( eq( si2->getObjCoefficients()[2],  0.0) );
    assert( eq( si2->getObjCoefficients()[3],  0.0) );
    assert( eq( si2->getObjCoefficients()[4],  2.0) );
    assert( eq( si2->getObjCoefficients()[5],  0.0) );
    assert( eq( si2->getObjCoefficients()[6],  0.0) );
    assert( eq( si2->getObjCoefficients()[7], -1.0) );
    
    // Test getting and setting of objective offset
    double objOffset;
    bool supported = si2->getDblParam(OsiObjOffset,objOffset);
    assert( supported );
    if( !eq( objOffset, 0.0 ) )
      failureMessage(solverName,"getDblParam OsiObjOffset on cloned solverInterface");
    delete si2;
  }
  // end of clone testing
  
  // Test apply cuts method
  {      
    OsiSolverInterface & im = *(exmip1Si->clone()); 
    OsiCuts cuts;
    
    // Generate some cuts 
    {
      // Get number of rows and columns in model
      int nr=im.getNumRows();
      int nc=im.getNumCols();
      assert( nr == 5 );
      assert( nc == 8 );
      
      // Generate a valid row cut from thin air
      int c;
      {
        int *inx = new int[nc];
        for (c=0;c<nc;c++) inx[c]=c;
        double *el = new double[nc];
        for (c=0;c<nc;c++) el[c]=((double)c)*((double)c);
        
        OsiRowCut rc;
        rc.setRow(nc,inx,el);
        rc.setLb(-100.);
        rc.setUb(100.);
        rc.setEffectiveness(22);
        
        cuts.insert(rc);
        delete[]el;
        delete[]inx;
      }
      
      // Generate valid col cut from thin air
      {
        const double * oslColLB = im.getColLower();
        const double * oslColUB = im.getColUpper();
        int *inx = new int[nc];
        for (c=0;c<nc;c++) inx[c]=c;
        double *lb = new double[nc];
        double *ub = new double[nc];
        for (c=0;c<nc;c++) lb[c]=oslColLB[c]+0.001;
        for (c=0;c<nc;c++) ub[c]=oslColUB[c]-0.001;
        
        OsiColCut cc;
        cc.setLbs(nc,inx,lb);
        cc.setUbs(nc,inx,ub);
        
        cuts.insert(cc);
        delete [] ub;
        delete [] lb;
        delete [] inx;
      }
      
      {
        // Generate a row and column cut which are ineffective
        OsiRowCut * rcP= new OsiRowCut;
        rcP->setEffectiveness(-1.);
        cuts.insert(rcP);
        assert(rcP==NULL);
        
        OsiColCut * ccP= new OsiColCut;
        ccP->setEffectiveness(-12.);
        cuts.insert(ccP);
        assert(ccP==NULL);
      }
      {
        //Generate inconsistent Row cut
        OsiRowCut rc;
        const int ne=1;
        int inx[ne]={-10};
        double el[ne]={2.5};
        rc.setRow(ne,inx,el);
        rc.setLb(3.);
        rc.setUb(4.);
        assert(!rc.consistent());
        cuts.insert(rc);
      }
      {
        //Generate inconsistent col cut
        OsiColCut cc;
        const int ne=1;
        int inx[ne]={-10};
        double el[ne]={2.5};
        cc.setUbs(ne,inx,el);
        assert(!cc.consistent());
        cuts.insert(cc);
      }
      {
        // Generate row cut which is inconsistent for model m
        OsiRowCut rc;
        const int ne=1;
        int inx[ne]={10};
        double el[ne]={2.5};
        rc.setRow(ne,inx,el);
        assert(rc.consistent());
        assert(!rc.consistent(im));
        cuts.insert(rc);
      }
      {
        // Generate col cut which is inconsistent for model m
        OsiColCut cc;
        const int ne=1;
        int inx[ne]={30};
        double el[ne]={2.0};
        cc.setLbs(ne,inx,el);
        assert(cc.consistent());
        assert(!cc.consistent(im));
        cuts.insert(cc);
      }
      {
        // Generate col cut which is infeasible
        OsiColCut cc;
        const int ne=1;
        int inx[ne]={0};
        double el[ne]={2.0};
        cc.setUbs(ne,inx,el);
        cc.setEffectiveness(1000.);
        assert(cc.consistent());
        assert(cc.consistent(im));
        assert(cc.infeasible(im));
        cuts.insert(cc);
      }
    }
    assert(cuts.sizeRowCuts()==4);
    assert(cuts.sizeColCuts()==5);

   {  
      OsiSolverInterface::ApplyCutsReturnCode rc = im.applyCuts(cuts);
      assert( rc.getNumIneffective() == 2 );
      assert( rc.getNumApplied() == 2 );
      assert( rc.getNumInfeasible() == 1 );
      assert( rc.getNumInconsistentWrtIntegerModel() == 2 );
      assert( rc.getNumInconsistent() == 2 );
      assert( cuts.sizeCuts() == rc.getNumIneffective() +
        rc.getNumApplied() +
        rc.getNumInfeasible() +
        rc.getNumInconsistentWrtIntegerModel() +
        rc.getNumInconsistent() );
    }

    delete &im;
  }
  // end of apply cut method testing
    

  // Test setting solution
#if 0
  {
    OsiSolverInterface & m1 = *(exmip1Si->clone());
    int i;
    
    double * cs = new double[m1.getNumCols()];
    for ( i = 0;  i < m1.getNumCols();  i++ ) 
      cs[i] = i + .5;
    m1.setColSolution(cs);
    for ( i = 0;  i < m1.getNumCols();  i++ ) 
      assert(m1.getColSolution()[i] == i + .5);
    

    double * rs = new double[m1.getNumRows()];
    for ( i = 0;  i < m1.getNumRows();  i++ ) 
      rs[i] = i - .5;
  if ( !symSolverInterface ) {
    m1.setRowPrice(rs);
    for ( i = 0;  i < m1.getNumRows();  i++ ) 
      assert(m1.getRowPrice()[i] == i - .5);
  }
    delete [] cs;
    delete [] rs;
    delete &m1;
  }
#endif
  // Test column type methods

  if ( volSolverInterface ) {
     // Test for vol since it does not support this function
     failureMessage(solverName,
		    "column type methods all report continuous (OK for vol)");
  }
  else {
    OsiSolverInterface & fim = *(emptySi->clone());
    std::string fn = mpsDir+"exmip1";
    fim.readMps(fn.c_str(),"mps");

    // exmip1.mps has 2 integer variables with index 2 & 3
    assert(  fim.getNumIntegers() == 2 ) ;

    assert(  fim.isContinuous(0) );
    assert(  fim.isContinuous(1) );
    assert( !fim.isContinuous(2) );
    assert( !fim.isContinuous(3) );
    assert(  fim.isContinuous(4) );
    
    assert( !fim.isInteger(0) );
    assert( !fim.isInteger(1) );
    assert(  fim.isInteger(2) );
    assert(  fim.isInteger(3) );
    assert( !fim.isInteger(4) );
    
    assert( !fim.isBinary(0) );
    assert( !fim.isBinary(1) );
    assert(  fim.isBinary(2) );
    assert(  fim.isBinary(3) );
    assert( !fim.isBinary(4) );
    
    assert( !fim.isIntegerNonBinary(0) );
    assert( !fim.isIntegerNonBinary(1) );
    assert( !fim.isIntegerNonBinary(2) );
    assert( !fim.isIntegerNonBinary(3) );
    assert( !fim.isIntegerNonBinary(4) );
    
    // Test fractionalIndices

#if 0
    {
      double sol[]={1.0, 2.0, 2.9, 3.0, 4.0,0.0,0.0,0.0};
      fim.setColSolution(sol);
      OsiVectorInt fi = fim.getFractionalIndices(1e-5);
      assert( fi.size() == 1 );
      assert( fi[0]==2 );
      
      // Set integer variables very close to integer values
      sol[2]=5 + .00001/2.;
      sol[3]=8 - .00001/2.;
      fim.setColSolution(sol);
      fi = fim.getFractionalIndices(1e-5);
      assert( fi.size() == 0 );
      
      // Set integer variables close, but beyond tolerances
      sol[2]=5 + .00001*2.;
      sol[3]=8 - .00001*2.;
      fim.setColSolution(sol);
      fi = fim.getFractionalIndices(1e-5);
      assert( fi.size() == 2 );
      assert( fi[0]==2 );
      assert( fi[1]==3 );
    }
#endif    
    // Change data so column 2 & 3 are integerNonBinary
    fim.setColUpper(2,5.0);
    assert( eq(fim.getColUpper()[2],5.0) );
    fim.setColUpper(3,6.0);
    assert( eq(fim.getColUpper()[3],6.0) );
    assert( !fim.isBinary(0) );
    assert( !fim.isBinary(1) );
    if( fim.isBinary(2) )
      failureMessage(solverName,"isBinary or setColUpper");
    if( fim.isBinary(3) )
      failureMessage(solverName,"isBinary or setColUpper");
    assert( !fim.isBinary(4) );
    
    if (fim.getNumIntegers() != 2)
      failureMessage(solverName,"getNumIntegers");

    assert( !fim.isIntegerNonBinary(0) );
    assert( !fim.isIntegerNonBinary(1) );
    if( !fim.isIntegerNonBinary(2) )
      failureMessage(solverName,"isIntegerNonBinary or setColUpper");
    if( !fim.isIntegerNonBinary(3) )
      failureMessage(solverName,"isIntegerNonBinary or setColUpper");
    assert( !fim.isIntegerNonBinary(4) );
    
    delete &fim;
  }

    
  // Test load and assign problem
  {
    {   
      OsiSolverInterface * base = exmip1Si->clone();
      OsiSolverInterface *  si1 = emptySi->clone(); 
      OsiSolverInterface *  si2 = emptySi->clone(); 
      OsiSolverInterface *  si3 = emptySi->clone(); 
      OsiSolverInterface *  si4 = emptySi->clone();  
      OsiSolverInterface *  si5 = emptySi->clone();  
      OsiSolverInterface *  si6 = emptySi->clone(); 
      OsiSolverInterface *  si7 = emptySi->clone();  
      OsiSolverInterface *  si8 = emptySi->clone(); 
        
      si1->loadProblem(*base->getMatrixByCol(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowSense(),base->getRightHandSide(),
		       base->getRowRange());
      si2->loadProblem(*base->getMatrixByRow(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowSense(),base->getRightHandSide(),
		       base->getRowRange());
      si3->loadProblem(*base->getMatrixByCol(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowLower(),base->getRowUpper() );
      si4->loadProblem(*base->getMatrixByCol(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowLower(),base->getRowUpper() );
      {
        double objOffset;
        base->getDblParam(OsiObjOffset,objOffset);
        si1->setDblParam(OsiObjOffset,objOffset);
        si2->setDblParam(OsiObjOffset,objOffset);
        si3->setDblParam(OsiObjOffset,objOffset);
        si4->setDblParam(OsiObjOffset,objOffset);
        si5->setDblParam(OsiObjOffset,objOffset);
        si6->setDblParam(OsiObjOffset,objOffset);
        si7->setDblParam(OsiObjOffset,objOffset);
        si8->setDblParam(OsiObjOffset,objOffset);
      }
      CoinPackedMatrix * pm = new CoinPackedMatrix(*base->getMatrixByCol());
      double * clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      double * cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      double * objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      double * rlb = new double[base->getNumRows()];
      std::copy(base->getRowLower(),
		base->getRowLower()+base->getNumRows(),rlb);
      double * rub = new double[base->getNumRows()];
      std::copy(base->getRowUpper(),
		base->getRowUpper()+base->getNumRows(),rub);
      si5->assignProblem(pm,clb,cub,objc,rlb,rub);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rlb==NULL);
      assert(rub==NULL);
        
      pm = new CoinPackedMatrix(*base->getMatrixByRow());
      clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      rlb = new double[base->getNumRows()];
      std::copy(base->getRowLower(),
		base->getRowLower()+base->getNumRows(),rlb);
      rub = new double[base->getNumRows()];
      std::copy(base->getRowUpper(),
		base->getRowUpper()+base->getNumRows(),rub);
      si6->assignProblem(pm,clb,cub,objc,rlb,rub);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rlb==NULL);
      assert(rub==NULL);      
        
      pm = new CoinPackedMatrix(*base->getMatrixByCol());
      clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      char * rsen = new char[base->getNumRows()];
      std::copy(base->getRowSense(),
		base->getRowSense()+base->getNumRows(),rsen);
      double * rhs = new double[base->getNumRows()];
      std::copy(base->getRightHandSide(),
		base->getRightHandSide()+base->getNumRows(),rhs);
      double * rng = new double[base->getNumRows()];
      std::copy(base->getRowRange(),
		base->getRowRange()+base->getNumRows(),rng);
      si7->assignProblem(pm,clb,cub,objc,rsen,rhs,rng);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rsen==NULL);
      assert(rhs==NULL);
      assert(rng==NULL);
        
      pm = new CoinPackedMatrix(*base->getMatrixByCol());
      clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      rsen = new char[base->getNumRows()];
      std::copy(base->getRowSense(),
		base->getRowSense()+base->getNumRows(),rsen);
      rhs = new double[base->getNumRows()];
      std::copy(base->getRightHandSide(),
		base->getRightHandSide()+base->getNumRows(),rhs);
      rng = new double[base->getNumRows()];
      std::copy(base->getRowRange(),
		base->getRowRange()+base->getNumRows(),rng);
      si8->assignProblem(pm,clb,cub,objc,rsen,rhs,rng);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rsen==NULL);
      assert(rhs==NULL);
      assert(rng==NULL);
     
        
      // Create an indices vector        
      CoinPackedVector basePv,pv;
      assert(base->getNumCols()<10);
      assert(base->getNumRows()<10);
      int indices[10];
      int i;
      for (i=0; i<10; i++) indices[i]=i;
        
      // Test solve methods.
      try {
	base->initialSolve();
	si1->initialSolve();
	si2->initialSolve();
	si3->initialSolve();
	si4->initialSolve();
	si5->initialSolve();
	si6->initialSolve();
	si7->initialSolve();
	si8->initialSolve();
      }        
      catch (CoinError e) {
#ifdef COIN_HAS_VOL
	// Vol solver interface is expected to throw
	// an error if the data has a ranged row.
          
	// Check that using Vol SI
	OsiVolSolverInterface * vsi =
	  dynamic_cast<OsiVolSolverInterface *>(base);
	assert( vsi != NULL );        
          
	// Test that there is non-zero range
	basePv.setFull(base->getNumRows(),base->getRowRange());
	pv.setConstant( base->getNumRows(), indices, 0.0 );
	assert(!basePv.isEquivalent(pv)); 
#else
	assert(0==1);
#endif
      }
      
      // Test WriteMps
      
      {

	OsiSolverInterface *  si1 = emptySi->clone(); 
	OsiSolverInterface *  si2 = emptySi->clone(); 
	si1->readMps(fn.c_str(),"mps");
	si1->writeMpsNative("test.out",NULL,NULL);
	si1->writeMps("test2","out");
	si2->readMps("test.out","");
	bool solved = true;
	try {
	   si1->initialSolve();
	}
	catch (CoinError e) {
	   if (e.className() != "OsiVolSolverInterface") {
	      printf("Couldn't solve initial LP in testing WriteMps\n");
	      abort();
	   }
	   solved = false;
	}
	if (solved) {
	   si2->initialSolve();
	   double soln = si1->getObjValue();       
	   CoinRelFltEq eq(1.0e-8) ;
	   assert( eq(soln,si2->getObjValue()));       
	}
	delete si1;
	delete si2;
      }
        
      // Test WriteLp
      
      {
#if 0
	OsiSolverInterface *  si1 = emptySi->clone(); 
	OsiSolverInterface *  si2 = emptySi->clone(); 
	si1->readMps(fn.c_str(),"mps");
	si1->writeLpNative("test.lp",NULL,NULL,1.0e-9,10,8);
	si1->writeLp("test2");
	si2->readLp("test.lp");
	bool solved = true;
	try {
	   si1->initialSolve();
	}
	catch (CoinError e) {
	   if (e.className() != "OsiVolSolverInterface") {
	      printf("Couldn't solve initial LP in testing WriteMps\n");
	      abort();
	   }
	   solved = false;
	}
	if (solved) {
	   si2->initialSolve();
	   double soln = si1->getObjValue();       
	   CoinRelFltEq eq(1.0e-8) ;
	   assert( eq(soln,si2->getObjValue()));       
	}
	delete si1;
	delete si2;
#endif
      }
        
      // Test collower
      basePv.setVector(base->getNumCols(),indices,base->getColLower());
      pv.setVector( si1->getNumCols(),indices, si1->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumCols(),indices, si2->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumCols(),indices, si3->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumCols(),indices, si4->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumCols(),indices, si5->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumCols(),indices, si6->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumCols(),indices, si7->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumCols(),indices, si8->getColLower());
      assert(basePv.isEquivalent(pv));
        
      // Test colupper
      basePv.setVector(base->getNumCols(),indices,base->getColUpper());
      pv.setVector( si1->getNumCols(),indices, si1->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumCols(),indices, si2->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumCols(),indices, si3->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumCols(),indices, si4->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumCols(),indices, si5->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumCols(),indices, si6->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumCols(),indices, si7->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumCols(),indices, si8->getColUpper());
      assert(basePv.isEquivalent(pv));
        
      // Test getObjCoefficients
      basePv.setVector(base->getNumCols(),indices,base->getObjCoefficients());
      pv.setVector( si1->getNumCols(),indices, si1->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumCols(),indices, si2->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumCols(),indices, si3->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumCols(),indices, si4->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumCols(),indices, si5->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumCols(),indices, si6->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumCols(),indices, si7->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumCols(),indices, si8->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
                
      // Test rowrhs
      basePv.setFull(base->getNumRows(),base->getRightHandSide());
      pv.setFull( si1->getNumRows(), si1->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si2->getNumRows(), si2->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si3->getNumRows(), si3->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si4->getNumRows(), si4->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si5->getNumRows(), si5->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si6->getNumRows(), si6->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si7->getNumRows(), si7->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si8->getNumRows(), si8->getRightHandSide());
      assert(basePv.isEquivalent(pv));
        
      // Test rowrange
      basePv.setFull(base->getNumRows(),base->getRowRange());
      pv.setFull( si1->getNumRows(), si1->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si2->getNumRows(), si2->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si3->getNumRows(), si3->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si4->getNumRows(), si4->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si5->getNumRows(), si5->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si6->getNumRows(), si6->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si7->getNumRows(), si7->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si8->getNumRows(), si8->getRowRange());
      assert(basePv.isEquivalent(pv));
        
      // Test row sense
      {
	const char * cb = base->getRowSense();
	const char * c1 = si1->getRowSense();
	const char * c2 = si2->getRowSense();
	const char * c3 = si3->getRowSense();
	const char * c4 = si4->getRowSense();
	const char * c5 = si5->getRowSense();
	const char * c6 = si6->getRowSense();
	const char * c7 = si7->getRowSense();
	const char * c8 = si8->getRowSense();
	int nr = base->getNumRows();
	for ( i=0; i<nr; i++ ) {
	  assert( cb[i]==c1[i] );
	  assert( cb[i]==c2[i] );
	  assert( cb[i]==c3[i] );
	  assert( cb[i]==c4[i] );
	  assert( cb[i]==c5[i] );
	  assert( cb[i]==c6[i] );
	  assert( cb[i]==c7[i] );
	  assert( cb[i]==c8[i] );
	}
      }
        
      // Test rowlower
      basePv.setVector(base->getNumRows(),indices,base->getRowLower());
      pv.setVector( si1->getNumRows(),indices, si1->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumRows(),indices, si2->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumRows(),indices, si3->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumRows(),indices, si4->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumRows(),indices, si5->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumRows(),indices, si6->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumRows(),indices, si7->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumRows(),indices, si8->getRowLower());
      assert(basePv.isEquivalent(pv));
        
      // Test rowupper
      basePv.setVector(base->getNumRows(),indices,base->getRowUpper());
      pv.setVector( si1->getNumRows(),indices, si1->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumRows(),indices, si2->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumRows(),indices, si3->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumRows(),indices, si4->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumRows(),indices, si5->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumRows(),indices, si6->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumRows(),indices, si7->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumRows(),indices, si8->getRowUpper());
      assert(basePv.isEquivalent(pv)); 
        
      // Test Constraint Matrix
      assert( base->getMatrixByCol()->isEquivalent(*si1->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si1->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si2->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si2->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si3->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si3->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si4->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si4->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si5->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si5->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si6->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si6->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si7->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si7->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si8->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si8->getMatrixByRow()) );
        
      // Test Objective Value
      assert( eq(base->getObjValue(),si1->getObjValue()) );
      assert( eq(base->getObjValue(),si2->getObjValue()) );
      assert( eq(base->getObjValue(),si3->getObjValue()) );
      assert( eq(base->getObjValue(),si4->getObjValue()) );
      assert( eq(base->getObjValue(),si5->getObjValue()) );
      assert( eq(base->getObjValue(),si6->getObjValue()) );
      assert( eq(base->getObjValue(),si7->getObjValue()) );
      assert( eq(base->getObjValue(),si8->getObjValue()) );
        
      // Clean-up
      delete si8;
      delete si7;
      delete si6;
      delete si5;
      delete si4;
      delete si3;
      delete si2;
      delete si1;
      delete base;
    }     
    // Test load/assign with null parms
    {
      //Load problem with row bounds and all rims at defaults
      {
	OsiSolverInterface *  si = emptySi->clone(); 
          
	si->loadProblem(*exmip1Si->getMatrixByCol(),NULL,NULL,NULL,NULL,NULL);
          
	// Test column settings
	assert(si->getNumCols()==exmip1Si->getNumCols() );
	for ( i=0; i<si->getNumCols(); i++ ) {
	  assert( eq(si->getColLower()[i],0.0) );
	  //	  assert( eq(si->getColUpper()[i],si->getInfinity()) );
	  assert( eq(si->getObjCoefficients()[i],0.0) );
	}
	// Test row settings
	assert(si->getNumRows()==exmip1Si->getNumRows() );
	const double * rh = si->getRightHandSide();
	const double * rr = si->getRowRange();
	const char * rs = si->getRowSense();
	const double * rl = si->getRowLower();
	const double * ru = si->getRowUpper();
	for ( i=0; i<si->getNumRows(); i++ ) {
	  assert( eq(rh[i],0.0) );
	  assert( eq(rr[i],0.0) );
	  assert( 'N'==rs[i] );
	  assert( eq(rl[i],-si->getInfinity()) );
	  assert( eq(ru[i], si->getInfinity()) );
	}
          
	delete si;
      }
      //Load problem with row rhs and all rims at defaults
      {
	OsiSolverInterface *  si = emptySi->clone(); 
          
	si->loadProblem(*exmip1Si->getMatrixByRow(),
			NULL,NULL,NULL,
			exmip1Si->getRowSense(),
			exmip1Si->getRightHandSide(),
			exmip1Si->getRowRange());
	// Test column settings
	assert(si->getNumCols()==exmip1Si->getNumCols() );
	for ( i=0; i<si->getNumCols(); i++ ) {
	  assert( eq(si->getColLower()[i],0.0) );
	  //	  assert( eq(si->getColUpper()[i],si->getInfinity()) );
	  assert( eq(si->getObjCoefficients()[i],0.0) );
	}
	// Test row settings
	assert(si->getNumRows()==exmip1Si->getNumRows() );
	for ( i=0; i<si->getNumRows(); i++ ) {
	  char s = si->getRowSense()[i];
	  assert( eq(si->getRightHandSide()[i],
		     exmip1Si->getRightHandSide()[i]) );
	  assert( eq(si->getRowRange()[i],
		     exmip1Si->getRowRange()[i]) );
	  assert( s==exmip1Si->getRowSense()[i] );
            
	  if ( s=='G' ) {
	    assert( eq(si->getRowLower()[i],
		       exmip1Si->getRightHandSide()[i]) );
	    assert( eq(si->getRowUpper()[i],
		       si->getInfinity()) );
	  }
	  else if ( s=='L' ) {
	    assert( eq(si->getRowLower()[i],
		       -si->getInfinity()) );
	    assert( eq(si->getRowUpper()[i],
		       exmip1Si->getRightHandSide()[i]) );
	  }
	  else if ( s=='E' ) {
	    assert( eq(si->getRowLower()[i],
		       si->getRowUpper()[i]) );
	    assert( eq(si->getRowUpper()[i],
		       exmip1Si->getRightHandSide()[i]) );
	  }
	  else if ( s=='N' ) {
	    assert( eq(si->getRowLower()[i], -si->getInfinity()) );
	    assert( eq(si->getRowUpper()[i],  si->getInfinity()) );
	  }
	  else if ( s=='R' ) {
	    assert( eq(si->getRowLower()[i],
		       exmip1Si->getRightHandSide()[i] -
		       exmip1Si->getRowRange()[i]) );
	    assert( eq(si->getRowUpper()[i],
		       exmip1Si->getRightHandSide()[i]) );
	  }
	}
          
	delete si;
      }
      // Test adding rows to NULL
      {
	OsiSolverInterface *  si = emptySi->clone();
	int i;

	//Matrix
	int column[]={0,1,2};
	double row1E[]={4.0,7.0,5.0};
	double row2E[]={7.0,4.0,5.0};
	CoinPackedVector row1(3,column,row1E);
	CoinPackedVector row2(3,column,row2E);

	double objective[]={5.0,6.0,5.5};

  {
	  // Add empty columns
	  for (i=0;i<3;i++) 
	    si->addCol(CoinPackedVector(),0.0,10.0,objective[i]);

	  // Add rows
	  si->addRow(row1,2.0,100.0);
	  si->addRow(row2,2.0,100.0);

    // Vol can not solve problem of this form
    if ( !volSolverInterface ) {
	    // solve
	    si->initialSolve();

      CoinRelFltEq eq(1.0e-7) ;
      double objValue = si->getObjValue();
	    if ( !eq(objValue,2.0) )
        failureMessage(solverName,"getObjValue after adding empty cols and then rows.");;
    }
  }

  delete si;
      }
      // Test adding rows to NULL - alternative format
      {
	OsiSolverInterface *  si = emptySi->clone();
	int i;

	//Matrix
	int column[]={0,1,2,0,1,2};
	double row1E[]={4.0,7.0,5.0};
	double row2E[]={7.0,4.0,5.0};
	double row12E[]={4.0,7.0,5.0,7.0,4.0,5.0};
	int starts[]={0,3,6};
	double ub[]={100.0,100.0};
	//CoinPackedVector row1(3,column,row1E);
	//CoinPackedVector row2(3,column,row2E);

	double objective[]={5.0,6.0,5.5};

	{
	  // Add empty columns
	  for (i=0;i<3;i++) 
	    si->addCol(CoinPackedVector(),0.0,10.0,objective[i]);
	  
	  // Add rows
	  si->addRow(3,column,row1E,2.0,100.0);
	  si->addRow(3,column,row2E,2.0,100.0);
	  // and again
	  si->addRows(2,starts,column,row12E,NULL,ub);
	  
	  // Vol can not solve problem of this form
	  if ( !volSolverInterface ) {
	    // solve
	    si->initialSolve();
	    
	    CoinRelFltEq eq(1.0e-7) ;
	    double objValue = si->getObjValue();
	    if ( !eq(objValue,2.0) )
	      failureMessage(solverName,"getObjValue after adding empty cols and then rows.");;
	  }
	}
	
	delete si;
      }
      // Test adding columns to NULL
      {
	OsiSolverInterface *  si = emptySi->clone();
	int i;

	//Matrix
	int row[]={0,1};
	double col1E[]={4.0,7.0};
	double col2E[]={7.0,4.0};
	double col3E[]={5.0,5.0};
	CoinPackedVector col1(2,row,col1E);
	CoinPackedVector col2(2,row,col2E);
	CoinPackedVector col3(2,row,col3E);

	double objective[]={5.0,6.0,5.5};
  {
     // Add empty rows
     for (i=0;i<2;i++) 
	si->addRow(CoinPackedVector(),2.0,100.0);

     // Add columns
     if ( volSolverInterface ) {
	// FIXME: this test could be done w/ the volume, but the rows must not
	// be ranged.
	failureMessage(solverName,"addCol add columns to null");
     }
     else {
	si->addCol(col1,0.0,10.0,objective[0]);
	si->addCol(col2,0.0,10.0,objective[1]);
	si->addCol(col3,0.0,10.0,objective[2]);

	// solve
	si->initialSolve();

	CoinRelFltEq eq(1.0e-7) ;      
	double objValue = si->getObjValue();
	if ( !eq(objValue,2.0) )
	   failureMessage(solverName,"getObjValue after adding empty rows and then cols.");

     }
  }
	delete si;
      }
      // Test Simplex interface
      {
        OsiSolverInterface * si = emptySi->clone();
        std::string solverName;
        si->getStrParam(OsiSolverName,solverName);
        
        if (si->canDoSimplexInterface()==2) {
          // solve an lp by hand
          
          std::string fn = mpsDir+"p0033";
          si->readMps(fn.c_str(),"mps");
          si->setObjSense(-1.0);
          si->initialSolve();
          si->setObjSense(1.0);
          // enable special mode
          si->enableSimplexInterface(true);
          // we happen to know that variables are 0-1 and rows are L
          int numberIterations=0;
          int numberColumns = si->getNumCols();
          int numberRows = si->getNumRows();
          double * fakeCost = new double[numberColumns];
          double * duals = new double [numberRows];
          double * djs = new double [numberColumns];
          const double * solution = si->getColSolution();
          memcpy(fakeCost,si->getObjCoefficients(),numberColumns*sizeof(double));
          while (1) {
            const double * dj;
            const double * dual;
            if ((numberIterations&1)==0) {
              // use given ones
              dj = si->getReducedCost();
              dual = si->getRowPrice();
            } else {
              // create
              dj = djs;
              dual = duals;
              si->getReducedGradient(djs,duals,fakeCost);
            }
            int i;
            int colIn=9999;
            int direction=1;
            double best=1.0e-6;
            // find most negative reduced cost
            // Should check basic - but should be okay on this problem
            for (i=0;i<numberRows;i++) {
              double value=dual[i];
              if (value>best) {
                direction=-1;
                best=value;
                colIn=-i-1;
              }
            }
            for (i=0;i<numberColumns;i++) {
              double value=dj[i];
              if (value<-best&&solution[i]<1.0e-6) {
                direction=1;
                best=-value;
                colIn=i;
              } else if (value>best&&solution[i]>1.0-1.0e-6) {
                direction=-1;
                best=value;
                colIn=i;
              }
            }
            if (colIn==9999)
              break; // should be optimal
            int colOut;
            int outStatus;
            double theta;
            assert(!si->primalPivotResult(colIn,direction,colOut,outStatus,theta,NULL));
            printf("out %d, direction %d theta %g\n",
                   colOut,outStatus,theta);
            numberIterations++;
          }
          delete [] fakeCost;
          delete [] duals;
          delete [] djs;
          // exit special mode
          si->disableSimplexInterface();
          si->resolve();
          assert (!si->getIterationCount());
          si->setObjSense(-1.0);
          si->initialSolve();
          std::cout<<solverName<<" passed OsiSimplexInterface test"<<std::endl;
        } else {
          std::cout<<solverName<<" has no OsiSimplexInterface"<<std::endl;
        }
        delete si;
      }

      // Test adding columns to NULL - alternative format
      {
	OsiSolverInterface *  si = emptySi->clone();
	int i;

	//Matrix
	int row[]={0,1};
	double col1E[]={4.0,7.0};
	double col23E[]={7.0,4.0,5.0,5.0};
	int row23E[]={0,1,0,1};
	int start23E[]={0,2,4};
	double ub23E[]={10.0,10.0};
	//CoinPackedVector col1(2,row,col1E);
	//CoinPackedVector col2(2,row,col2E);
	//CoinPackedVector col3(2,row,col3E);

	double objective[]={5.0,6.0,5.5};
	{
	  // Add empty rows
	  for (i=0;i<2;i++) 
	    si->addRow(CoinPackedVector(),2.0,100.0);
	  
	  // Add columns
	  if ( volSolverInterface ) {
	    // FIXME: this test could be done w/ the volume, but the rows must not
	    // be ranged.
	    failureMessage(solverName,"addCol add columns to null");
	  }
	  else {
	    si->addCol(2,row,col1E,0.0,10.0,objective[0]);
	    si->addCols(2,start23E,row23E,col23E,NULL,ub23E,objective+1);
	    
	    // solve
	    si->initialSolve();
	    
	    CoinRelFltEq eq(1.0e-7) ;      
	    double objValue = si->getObjValue();
	    if ( !eq(objValue,2.0) )
	      failureMessage(solverName,"getObjValue after adding empty rows and then cols.");
	    
	  }
	}
	delete si;
      }
    }
  }
#ifdef COIN_OPBDP
  // test Opbdp interface
  {
    OsiSolverInterface * si = emptySi->clone();
    std::string solverName;
    si->getStrParam(OsiSolverName,solverName);
    std::string fn = mpsDir+"p0033";
    si->readMps(fn.c_str(),"mps");
    int numberFound = solveOpbdp(si);
    assert (numberFound==1);
    printf("Objective value %g\n",si->getObjValue());
    unsigned int ** array = solveOpbdp(si,numberFound);
    printf("%d solutions found\n",numberFound);
    int numberColumns = si->getNumCols();
    for (int i=0;i<numberFound;i++) {
      unsigned int * thisArray = array[i];
      for (int j=0;j<numberColumns;j++) 
        if (atOne(j,thisArray))
          printf(" %d",j);
      printf("\n");
      delete [] thisArray;
    }
    delete [] array;
    delete si;
  }
#endif

  // Add a Laci suggested test case
  // Load in a problem as column ordered matrix, 
  // extract the row ordered copy, 
  // add a row, 
  // extract the row ordered copy again and test whether it's ok. 
  // (the same can be done with reversing the role
  //  of row and column ordered.)
  {
    OsiSolverInterface *  si = emptySi->clone(); 
      
    si->loadProblem(
		    *(exmip1Si->getMatrixByCol()),
		    exmip1Si->getColLower(),
		    exmip1Si->getColUpper(),
		    exmip1Si->getObjCoefficients(),
		    exmip1Si->getRowSense(),
		    exmip1Si->getRightHandSide(),
		    exmip1Si->getRowRange() );

    CoinPackedMatrix pm1 = *(si->getMatrixByRow());

    // Get a row of the matrix to make a cut
    CoinPackedVector pv =exmip1Si->getMatrixByRow()->getVector(1);
    pv.setElement(0,3.14*pv.getElements()[0]);

    OsiRowCut rc;
    rc.setRow( pv );
    rc.setLb( exmip1Si->getRowLower()[1]-0.5 );
    rc.setUb( exmip1Si->getRowUpper()[1]-0.5 );

    OsiCuts cuts;
    cuts.insert(rc);

    si->applyCuts(cuts);
      
    CoinPackedMatrix pm2 = *(si->getMatrixByRow());

    assert(pm1.getNumRows()==pm2.getNumRows()-1);
    int i;
    for( i=0; i<pm1.getNumRows(); ++i ) {
      assert( pm1.getVector(i) == pm2.getVector(i) );
    }
    // Test that last row of pm2 is same as added cut
    assert( pm2.getVector(pm2.getNumRows()-1).isEquivalent(pv) );

    delete si;
  }
  {
    OsiSolverInterface *  si = emptySi->clone(); 
      
    si->loadProblem(
		    *(exmip1Si->getMatrixByRow()),
		    exmip1Si->getColLower(),
		    exmip1Si->getColUpper(),
		    exmip1Si->getObjCoefficients(),
		    exmip1Si->getRowLower(),
		    exmip1Si->getRowUpper() );

    CoinPackedMatrix pm1 = *(si->getMatrixByCol());

    // Get a row of the matrix to make a cut
    CoinPackedVector pv =exmip1Si->getMatrixByRow()->getVector(1);
    pv.setElement(0,3.14*pv.getElements()[0]);

    OsiRowCut rc;
    rc.setRow( pv );
    rc.setLb( exmip1Si->getRowLower()[1]-0.5 );
    rc.setUb( exmip1Si->getRowUpper()[1]-0.5 );

    OsiCuts cuts;
    cuts.insert(rc);

    si->applyCuts(cuts);
      
    CoinPackedMatrix pm2 = *(si->getMatrixByCol());

    assert( pm1.isColOrdered() );
    assert( pm2.isColOrdered() );
    assert( pm1.getNumRows()==pm2.getNumRows()-1 );

    CoinPackedMatrix pm1ByRow;
    pm1ByRow.reverseOrderedCopyOf(pm1);
    CoinPackedMatrix pm2ByRow;
    pm2ByRow.reverseOrderedCopyOf(pm2);

    assert( !pm1ByRow.isColOrdered() );
    assert( !pm2ByRow.isColOrdered() );      
    assert( pm1ByRow.getNumRows()==pm2ByRow.getNumRows()-1 );
    assert( pm1.getNumRows() == pm1ByRow.getNumRows() );
    assert( pm2.getNumRows() == pm2ByRow.getNumRows() );

    int i;
    for( i=0; i<pm1ByRow.getNumRows(); ++i ) {
      assert( pm1ByRow.getVector(i) == pm2ByRow.getVector(i) );
    }
    // Test that last row of pm2 is same as added cut
    assert( pm2ByRow.getVector(pm2ByRow.getNumRows()-1).isEquivalent(pv) );

    delete si;
  }
    
  delete exmip1Si;

  {
    // Testing parameter settings
    OsiSolverInterface *  si = emptySi->clone();
    int i;
    int ival;
    double dval;
    bool hint;
    OsiHintStrength hintStrength;
    assert(si->getIntParam(OsiLastIntParam, ival) == false);
    assert(si->getDblParam(OsiLastDblParam, dval) == false);
    assert(si->getHintParam(OsiLastHintParam, hint) == false);
    assert(si->setIntParam(OsiLastIntParam, 0) == false);
    assert(si->setDblParam(OsiLastDblParam, 0) == false);
    assert(si->setHintParam(OsiLastHintParam, false) == false);
    
    for (i = 0; i < OsiLastIntParam; ++i) {
      const bool exists = si->getIntParam(static_cast<OsiIntParam>(i), ival);
      // existence and test should result in the same
      assert(!exists ^ testIntParam(si, i, -1));
      assert(!exists ^ testIntParam(si, i, 0));
      assert(!exists ^ testIntParam(si, i, 1));
      assert(!exists ^ testIntParam(si, i, 9999999));
      assert(!exists ^ testIntParam(si, i, INT_MAX));
      if (exists)
        assert(si->getIntParam(static_cast<OsiIntParam>(i), ival));
    }
    // Test for glpk since it dies on this test
    if ( glpkSolverInterface ) {
      failureMessage(solverName,"[g,s]etDblParam");
    }
    else {
      for (i = 0; i < OsiLastDblParam; ++i) {
        const bool exists = si->getDblParam(static_cast<OsiDblParam>(i), dval);
        // existence and test should result in the same
        assert(!exists ^ testDblParam(si, i, -1e50));
        assert(!exists ^ testDblParam(si, i, -1e10));
        assert(!exists ^ testDblParam(si, i, -1));
        assert(!exists ^ testDblParam(si, i, -1e-4));
        assert(!exists ^ testDblParam(si, i, -1e-15));
        assert(!exists ^ testDblParam(si, i, 1e50));
        assert(!exists ^ testDblParam(si, i, 1e10));
        assert(!exists ^ testDblParam(si, i, 1));
        assert(!exists ^ testDblParam(si, i, 1e-4));
        assert(!exists ^ testDblParam(si, i, 1e-15));
        if (exists)
          assert(si->setDblParam(static_cast<OsiDblParam>(i), dval));
      }
    }

    // test hints --- see testHintParam for detailed explanation.

    { int throws = 0 ;

      for (i = 0 ; i < OsiLastHintParam ; ++i)
      { const bool exists =
	  si->getHintParam(static_cast<OsiHintParam>(i),hint,hintStrength) ;

	assert(!exists ^ testHintParam(si,i,true,OsiHintIgnore,&throws)) ;
	assert(!exists ^ testHintParam(si,i,true,OsiHintTry,&throws)) ;
	assert(!exists ^ testHintParam(si,i,false,OsiHintTry,&throws)) ;
	assert(!exists ^ testHintParam(si,i,true,OsiHintDo,&throws)) ;
	assert(!exists ^ testHintParam(si,i,false,OsiHintDo,&throws)) ;
	assert(!exists ^ testHintParam(si,i,true,OsiForceDo,&throws)) ;
	assert(!exists ^ testHintParam(si,i,false,OsiForceDo,&throws)) ; }
      
      std::cerr << "Checked " << OsiLastHintParam <<
		   " hints x (true, false) at strength OsiForceDo; " <<
		   throws << " throws." << std::endl ;
    }

    delete si;
  }
  
  // Test case submitted by Vivian De Smedt (slightly modifed to work with
  // Vol Algorithm).

  {
    OsiSolverInterface *s = emptySi->clone();
    double dEmpty = 0;
    int iEmpty = 0;
    CoinBigIndex iEmpty2 = 0;
    //char cEmpty = '?';
    
    s->loadProblem(0, 0, &iEmpty2, &iEmpty, &dEmpty, &dEmpty, &dEmpty, &dEmpty, &dEmpty, &dEmpty);
    double inf = s->getInfinity();
    CoinPackedVector c;
    
    s->addCol(c, 0, 10, 3);
    s->addCol(c, 0, 10, 1);
    
    CoinPackedVector r1;
    r1.insert(0, 2);
    r1.insert(1, 1);
    s->addRow(r1, -inf, 10);
    
    CoinPackedVector r2;
    r2.insert(0, 1);
    r2.insert(1, 3);
    s->addRow(r2, -inf, 15);
    
    s->setObjSense(-1);

    s->initialSolve() ;
    const double * colSol = s->getColSolution();
    // Don't test for exact answer, because Vol algorithm
    // only returns an appoximate solution
    if ( colSol[0]<4.5 )
      failureMessage(*s,"colsol[0] bad value");
    if ( colSol[1]>0.5 )
      failureMessage(*s,"colsol[1] bad value");
    
    s->setObjCoeff(0, 1);
    s->setObjCoeff(1, 1);
    
    s->resolve();    
    colSol = s->getColSolution();
    // Don't test for exact answer, because Vol algorithm
    // only returns an appoximate solution
    if( colSol[0]<2.3 || colSol[0]>3.7 )
      failureMessage(*s,"colsol[0] bad value");
    if( colSol[1]<3.5 || colSol[1]>4.5 )
      failureMessage(*s,"colsol[1] bad value");
    delete s;
  }

  // Test presolve
  if ( !volSolverInterface /*&& !fmpSolverInterface*/ && !symSolverInterface ) {
    OsiSolverInterface * si = emptySi->clone();
    std::string fn = netlibDir+"25fv47";
    si->readMps(fn.c_str(),"mps");
    OsiSolverInterface * presolvedModel;
    OsiPresolve pinfo;
    int numberPasses=5; // can change this
    /* Use a tolerance of 1.0e-8 for feasibility, treat problem as 
       not being integer, do "numberpasses" passes */
    presolvedModel = pinfo.presolvedModel(*si,1.0e-8,false,numberPasses);
    assert(presolvedModel);
    // switch off presolve
    presolvedModel->setHintParam(OsiDoPresolveInInitial,false);
    presolvedModel->initialSolve();
    double objValue = presolvedModel->getObjValue(); 
    if( !eq(objValue,5.5018458883e+03) )
      failureMessage(solverName,"OsiPresolved model has wrong objective");
    pinfo.postsolve(true);


    delete presolvedModel;
    si->setHintParam(OsiDoPresolveInResolve,false);
    si->setHintParam(OsiDoDualInResolve,false);
    si->resolve();
    objValue = si->getObjValue(); 
    if( !eq(objValue,5.5018458883e+03) )
      failureMessage(solverName,"OsiPresolve - final objective wrong");
    if (si->getIterationCount())
      failureMessage(solverName,"OsiPresolve - minor error, needs iterations");
    delete si;
  }

  // Perform tests that are embodied in functions
  if ( !volSolverInterface && !symSolverInterface)
  {
    
    typedef bool (*TestFunction)(OsiSolverInterface*);
    std::vector<std::pair<TestFunction, const char*> > test_functions;
    test_functions.push_back(std::pair<TestFunction, const char*>(&test1VivianDeSmedt, "test1VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test2VivianDeSmedt, "test2VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test3VivianDeSmedt, "test3VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test4VivianDeSmedt, "test4VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test5VivianDeSmedt, "test5VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test6VivianDeSmedt, "test6VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test7VivianDeSmedt, "test7VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test8VivianDeSmedt, "test8VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test9VivianDeSmedt, "test9VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test10VivianDeSmedt,"test10VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test11VivianDeSmedt,"test11VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test12VivianDeSmedt,"test12VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test13VivianDeSmedt,"test13VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test14VivianDeSmedt,"test14VivianDeSmedt"));
    test_functions.push_back(std::pair<TestFunction, const char*>(&test15VivianDeSmedt,"test15VivianDeSmedt"));
    
    unsigned int i;
    for (i = 0; i < test_functions.size(); ++i) {
      OsiSolverInterface *s = emptySi->clone();
      const char * testName = test_functions[i].second;
      {
        bool test = test_functions[i].first(s);
        if (!test)
          failureMessage(*s, testName);
      }
      delete s;
    }
  }
  /*
    Orphan comment? If anyone happens to poke at the code that this belongs
    to, move it. My guess is it should go somewhere in the deSmedt tests.

    With this matrix we have a primal/dual infeas problem. Leaving the first
    row makes it primal feas, leaving the first col makes it dual feas.
    All vars are >= 0

    obj: -1  2 -3  4 -5 (min)

          0 -1  0  0 -2  >=  1
          1  0 -3  0  4  >= -2
          0  3  0 -5  0  >=  3
          0  0  5  0 -6  >= -4
          2 -4  0  6  0  >=  5
  */
}

