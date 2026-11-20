/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2026 Lehigh University. All Rights Reserved.                */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*===========================================================================*\
 * A C-style entry point to the feasibility jump reference implementation in
 * feasibilityjump.hh, which is kept unchanged. Everything C++ about
 * feasibility jump (its classes, std::vector, the callback) stays in these
 * two files; SYMPHONY's own code calls fj_solve() with plain arrays.
\*===========================================================================*/

#ifndef FEASIBILITYJUMP_ADAPTER_HH
#define FEASIBILITYJUMP_ADAPTER_HH

#include <chrono>

#include "CoinHelperFunctions.hpp" /* for CoinBigIndex */
#include "feasibilityjump.hh"

/*===========================================================================*\
 * Run feasibility jump on
 *    min obj x  s.t.  the m rows, lb <= x <= ub, x_j integer if is_int[j],
 * started from init. The rows are given row-wise (rowbeg, rowlen, rowind,
 * rowval) with sense 'E', 'L', 'G' or 'R' and right-hand side rhs; an 'R'
 * row is rhs - range <= ax <= rhs. Rows with any other sense are skipped.
 *
 * It stops at the first feasible solution, when the total effort reaches
 * max_effort, when the effort since the last improvement exceeds stall_limit
 * (if stall_limit > 0), or after timeleft seconds (if timeleft >= 0).
 *
 * Returns 1 and the solution in sol (n values) if it found one, 0 if not.
 * The effort used is returned in *effort.
\*===========================================================================*/

static int fj_solve(int n, const char *is_int, const double *lb,
		    const double *ub, const double *obj, int m,
		    const CoinBigIndex *rowbeg, const int *rowlen,
		    const int *rowind, const double *rowval, const char *sense,
		    const double *rhs, const double *range, const double *init,
		    double tol, long long max_effort, long long stall_limit,
		    double timeleft, double *sol, long long *effort)
{
   FeasibilityJumpSolver fj(0, 0, 1.0, tol, tol);
   int found = 0;

   for (int j = 0; j < n; j++){
      fj.addVar(is_int[j] ? VarType::Integer : VarType::Continuous,
		lb[j], ub[j], obj[j]);
   }
   for (int i = 0; i < m; i++){
      int len = rowlen[i];
      int *ind = const_cast<int *>(rowind + rowbeg[i]);
      double *val = const_cast<double *>(rowval + rowbeg[i]);
      switch (sense[i]){
       case 'E':
	 fj.addConstraint(RowType::Equal, rhs[i], len, ind, val, 0);
	 break;
       case 'L':
	 fj.addConstraint(RowType::Lte, rhs[i], len, ind, val, 0);
	 break;
       case 'G':
	 fj.addConstraint(RowType::Gte, rhs[i], len, ind, val, 0);
	 break;
       case 'R':
	 fj.addConstraint(RowType::Lte, rhs[i], len, ind, val, 0);
	 fj.addConstraint(RowType::Gte, rhs[i] - range[i], len, ind, val, 0);
	 break;
       default:
	 break;
      }
   }

   const auto start = std::chrono::steady_clock::now();
   *effort = 0;
   fj.solve(const_cast<double *>(init), [&](FJStatus status)
	    -> CallbackControlFlow {
      *effort = status.totalEffort;
      if (status.solution != nullptr){
	 for (int j = 0; j < n && j < status.numVars; j++){
	    sol[j] = status.solution[j];
	 }
	 found = 1;
	 return CallbackControlFlow::Terminate;
      }
      if (timeleft >= 0.0 &&
	  std::chrono::duration<double>(std::chrono::steady_clock::now() -
					start).count() >= timeleft){
	 return CallbackControlFlow::Terminate;
      }
      if (status.totalEffort >= max_effort){
	 return CallbackControlFlow::Terminate;
      }
      if (stall_limit > 0 && status.effortSinceLastImprovement > stall_limit){
	 return CallbackControlFlow::Terminate;
      }
      return CallbackControlFlow::Continue;
   });
   return found;
}

#endif
