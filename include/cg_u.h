/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2002 Ted Ralphs. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _CG_U_H
#define _CG_U_H

#include "BB_types.h"

int cg_send_cut PROTO((cut_data *new_cut));

/*===========================================================================*/
/*======================= User supplied functions ===========================*/
/*===========================================================================*/

/*__BEGIN_EXPERIMENTAL_SECTION__*/
int user_receive_cg_data PROTO((void **user, int dg_id, int *varnum));
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT THIS FOR PRODUCTION CODE*/
#if 0
int user_receive_cg_data PROTO((void **user, int dg_id));
#endif
int user_free_cg PROTO((void **user));
/*__BEGIN_EXPERIMENTAL_SECTION__*/
int user_find_cuts PROTO((void *user, int varnum, int iter_num, int level,
			  int index, double objval, int *indices,
			  double *values, double ub, double lpetol,
			  int *cutnum, char *status));
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
int user_find_cuts PROTO((void *user, int varnum, int iter_num, int level,
			  int index, double objval, int *indices,
			  double *values, double ub, double lpetol,
			  int *cutnum));
#endif
int user_receive_lp_solution_cg PROTO((void *user));
#ifdef CHECK_CUT_VALIDITY
int user_check_validity_of_cut PROTO((void *user, cut_data *new_cut));
#endif

#endif
