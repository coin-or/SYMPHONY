/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2005 Ted Ralphs. All Rights Reserved.                       */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifdef PRIMAL_HEURISTICS
#ifndef _SOLUTIONS_H_
#define _SOLUTIONS_H_
#include "sym_proto.h"
#include "sym_types.h"
#include "sym_lp.h"

/* feasible solutions managed by symphony */
int sp_add_solution PROTO((lp_prob *p, int cnt, int *indices, double *values, double obj_value, int bc_index));
int sp_delete_solution PROTO((sp_desc *sp, int position));
int sp_is_solution_in_sp PROTO((lp_prob *p, int cnt, int *indices, double *values, double obj_value));
int sp_initialize(tm_prob *tm);
#endif
#endif
