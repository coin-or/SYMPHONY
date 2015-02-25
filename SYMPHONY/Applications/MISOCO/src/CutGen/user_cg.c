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

/*===========================================================================*/

/*===========================================================================*\
 * This file contains user-written functions used by the cut generator
 * process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_cg_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_CG, COMPILE_IN_LP, or COMPILE_IN_TM is
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
   /* Fill in cut generation here if desired */
   int i;
   //int n = prob->si->getNumCols();
   int num_cones = prob->num_cones;
   int * cone_type = prob->cone_type;
   int * cone_size = prob->cone_size;
   int ** cone_members = prob->cone_members;
   //double * sol = calloc(n, sizeof(double));
   /* // fill sol */
   /* for (i=0; i<varnum; ++i) */
   /* double * par_sol; */
   /* for (i=0; i<num_cones; ++i) { */
   /*   // fill par_sol */
   /*   if (cone_type[i]==0) { */
   /*     separate_lorentz_cone(); */
   /*   } */
   /*   else if (cone_type[i]==1) { */
   /*     separate_rotated_cone(); */
   /*   } */
   /*   else { */
   /*     fprintf(stderr, "Unknown cone type!\n"); */
   /*     return (USER__ERROR); */
   /*   } */
   /*   // free par_sol */
   /*   free(par_sol); */
   /* } */
   // add cut to the model
   //free(sol);
   return(USER_SUCCESS);
}

// returns 0 if point is not epsilon feasible, nonzero otherwise
/* int separate_lorentz_cone() { */
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
/*   return 0; */
/* } */

/* int separate_rotated_cone() { */
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
/*   return 0; */
/* } */


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
