/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*                                                                           */
/*===========================================================================*/

/* $ID$ */

#ifndef _WARMSEARCH_H_
#define _WARMSEARCH_H_
/* rounding related */
#define  SYM_RND_FAIL  10
int warm_search(lp_prob *p, int * indices, double *values, int cnt, double lpetol, double *heur_solution, double &new_obj_val, int is_feasible);
bool var_is_non_zero(int i, int *indices, int varnum, int &var_index);
#endif
