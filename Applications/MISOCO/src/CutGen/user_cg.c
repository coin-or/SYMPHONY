/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/* SYMPHONY include files */
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_cg_u.h"

/* User include files */
#include "user.h"

/* C STL */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#define CONIC_APPROX_CUT 1001


#define TOL 1e-5

int separate_lorentz_cone(int size, double * sol, double * coef);
int separate_rotated_cone(int size, double * sol, double * coef);
void simple_separation(int size, double * p, double * coef);
double lorentz_cone_feasibility(int size, double const * point);
double rotated_cone_feasibility(int size, double const * point);

/*===========================================================================*/

/*===========================================================================*\
 * This file contains user-written functions used by the cut generator
 * process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_cg_data() and set up data structures. Note that this function is
 * only called if one of SYM_COMPILE_IN_CG, SYM_COMPILE_IN_LP, or SYM_COMPILE_IN_TM is
 * FALSE. For sequential computation, nothing is needed here.
\*===========================================================================*/

int user_receive_cg_data(void **user, int dg_id)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * If the user wants to fill in a customized routine for sending and receiving
 * the LP solution, it can be done here. For most cases, the default routines
 * are fine.
\*===========================================================================*/

int user_receive_lp_solution_cg(void *user)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Find cuts violated by a particular LP solution. This can be a fairly
 * involved function but the bottom line is that an LP solution comes in
 * and cuts go out. Remember, use the function cg_send_cut() to send cuts out
 * when they are found.
\*===========================================================================*/

int user_find_cuts(void *user, int varnum, int iter_num, int level,
		   int index, double objval, int *indices, double *values,
		   double ub, double etol, int *num_cuts, int *alloc_cuts,
		   cut_data ***cuts) {
   user_problem *prob = (user_problem *) user;
   int i;
   int j;
   int num_cones = prob->num_cones;
   int * cone_type = prob->cone_type;
   int * cone_size = prob->cone_size;
   int ** cone_members = prob->cone_members;
   int n = prob->colnum;
   int m = prob->rownum;
   /* cut_data * new_cut = NULL; */
   double * coef;
   int res;
   // cg_send_cut(new_cut, num_cuts, alloc_cuts, cuts);
   // fill sol using values and indices
   double * par_sol;
   for (i=0; i<num_cones; ++i) {
     /* new_cut = (cut_data *) calloc(1, sizeof(cut_data)); */
     /* new_cut->size = cone_size[i]; */
     /* new_cut->rhs = 0.0; */
     /* new_cut->range = ; */
     /* new_cut->type = CONIC_APPROX_CUT; */
     /* new_cut->sense = 'G'; */
     /* new_cut->deletable = ; */
     /* new_cut->branch = ; */
     /* new_cut->name = -2; */
     // fill par_sol
     coef = (double *) calloc(cone_size[i], sizeof(double));
     par_sol = (double *) calloc(cone_size[i], sizeof(double));
     for (j=0; j<cone_size[i]; j++) {
       par_sol[j] = prob->curr_solution[cone_members[i][j]];
     }
     if (cone_type[i]==0) {
       res = separate_lorentz_cone(cone_size[i], par_sol, coef);
     }
     else if (cone_type[i]==1) {
       res = separate_rotated_cone(cone_size[i], par_sol, coef);
     }
     else {
       fprintf(stderr, "Unknown cone type!\n");
       free(par_sol);
       return (USER_ERROR);
     }
     if (res==0) {
       // add cut
       /* cg_send_cut(new_cut, num_cuts, alloc_cuts, cuts); */
       // fprintf(stdout, "Cone %d is infeasible. Generating cut.\n", i);
       cg_add_explicit_cut(cone_size[i], cone_members[i], coef, 0.0, 0, 'G',
			   TRUE, num_cuts, alloc_cuts, cuts);
     }
     // free par_sol
     free(par_sol);
     /* free(new_cut); */
     free(coef);
   }
   return(USER_SUCCESS);
}

// returns 0 if point is not epsilon feasible, nonzero otherwise
// fills cut->coef
int separate_lorentz_cone(int size, double * sol, double * coef) {
  double feas = lorentz_cone_feasibility(size, sol);
  if (feas>-TOL)
    return 1;
  double * p = (double *) calloc(size, sizeof(double));
  for (int i=0; i<size; ++i) {
    p[i] = sol[i];
  }
  // 2. compute point on cone and coefficient from the point
  simple_separation(size, p, coef);
  // closest_point_separation(p, coef);
  // check if we actually cut the point
  double term1 = 0.0;
  int i;
  for (i=0; i<size; ++i) {
    term1 += coef[i]*sol[i];
  }
  if (term1> TOL) {
    fprintf(stderr, "Generated plane does not cut point.\n");
    assert(0);
  }
  // point is not feasible, add cut to cuts and return false
  // rhs is allways 0.0
  free(p);
  return 0;
}

void simple_separation(int size, double * p, double * coef) {
  // update p[0]
  int i;
  double term1 = 0.0;
  for (i=1; i<size; ++i) {
    term1 += p[i]*p[i];
  }
  p[0] = sqrt(term1);
  // coef array, [2p1, -2p2, -2p3, ... -2pn]
  coef[0] = 2.0*p[0];
  for (i=1; i<size; ++i) {
    coef[i] = -2.0*p[i];
  }
}

// fills cut->coef
int separate_rotated_cone(int size, double * sol, double * coef) {
/*   double feas = feasibility(point); */
/*   if (feas>-options()->get_dbl_option(TOL)) */
/*     return 1; */
/*   // coef array, [2x1, -2x2, -2x3, ... -2xn] */
/*   double * coef = new double[size_]; */
/*   double * p = new double[size_]; */
/*   for (int i=0; i<size_; ++i) { */
/*     p[i] = point[members_[i]]; */
/*   } */
/*   // 2. compute point on cone and coefficient from the point */
/*   simple_separation(p, coef); */
/*   // closest_point_separation(p, coef); */
/*   // check if we actually cut the point */
/*   double term1 = std::inner_product(coef, coef+size_, p, 0.0); */
/*   if (term1< -options()->get_dbl_option(TOL)) { */
/*     std::cerr << "Generated plane does not cut point." << std::endl; */
/*     throw std::exception(); */
/*   } */
/*   // point is not feasible, add cut to cuts_ and return false */
/*   // index is cone_members */
/*   // rhs is allways 0.0 */
/*   cut = new CoinPackedVector(size_, members_, coef); */
/*   delete[] coef; */
/*   delete[] p; */
  return 0;
}

// return the feasibility of point
// for Lorentz cones; x1-|x_2:n| or 2x1x2-|x_3:n|^2
double lorentz_cone_feasibility(int size, double const * point) {
  int i;
  double const * p = point;
  double feas;
  double term1 = p[0];
  double term2 = 0.0;
  for (i=1; i<size; ++i) {
    term2+=p[i]*p[i];
  }
  term2 = sqrt(term2);
  feas = term1-term2;
  return feas;
}

double rotated_cone_feasibility(int size, double const * point) {
  int i;
  double const * p = point;
  double feas;
  double term1 = 2.0*p[0]*p[1];
  double term2 = 0.0;
  for (i=2; i<size; ++i) {
    term2+=p[i]*p[i];
  }
  feas = term1-term2;
  return feas;
}
/*===========================================================================*/

/*===========================================================================*\
 * Free the user data structure. If the default setup is used with sequential
 * computation, nothing needs to be filled in here.
\*===========================================================================*/

int user_free_cg(void **user)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is an undocumented (for now) debugging feature which can allow the user
 * to identify the cut which cuts off a particular known feasible solution.
\*===========================================================================*/

#ifdef CHECK_CUT_VALIDITY
int user_check_validity_of_cut(void *user, cut_data *new_cut)
{
  return(USER_DEFAULT);
}
#endif
