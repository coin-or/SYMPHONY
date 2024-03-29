/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2005 Ted Ralphs. All Rights Reserved.                       */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _API_H
#define _API_H

#define COMPILING_FOR_MASTER

#include "sym_proto.h"
#include "sym_master.h"
#include "sym_messages.h"

/*===========================================================================*/
/*===================== Interface functions (master.c) ======================*/
/*===========================================================================*/

sym_environment *sym_open_environment PROTO((void));
int sym_set_defaults PROTO((sym_environment *env));
int sym_parse_command_line PROTO((sym_environment *env, int argc, 
				  char **argv));
int sym_set_user_data PROTO((sym_environment *env, void *user));
int sym_get_user_data PROTO((sym_environment *env, void **user));
int sym_read_mps PROTO((sym_environment *env, char *infile));
int sym_read_gmpl PROTO((sym_environment *env, char *modelfile, 
			 char *datafile));
int sym_load_problem PROTO((sym_environment *env));
int sym_find_initial_bounds PROTO((sym_environment *env));

int sym_solve PROTO((sym_environment *env));
int sym_warm_solve PROTO((sym_environment *env));
int sym_mc_solve PROTO((sym_environment *env));

int sym_create_permanent_cut_pools PROTO((sym_environment *env, int *cp_num));
int sym_close_environment PROTO((sym_environment *env));
int sym_explicit_load_problem PROTO((sym_environment *env, int numcols, 
				     int numrows, int *start, int *index, 
				     double *value, double *collb,
				     double *colub, char *is_int, double *obj,
				     double *obj2, char *rowsen,
				     double *rowrhs, double *rowrng,
				     char make_copy));   

int sym_is_abandoned PROTO((sym_environment *env));
int sym_is_proven_optimal PROTO((sym_environment *env));
int sym_is_proven_primal_infeasible PROTO((sym_environment *env));	 
int sym_is_iteration_limit_reached PROTO((sym_environment *env)); 
int sym_is_time_limit_reached PROTO((sym_environment *env));
int sym_is_target_gap_achieved PROTO((sym_environment *env));

int sym_get_status PROTO((sym_environment *env));
int sym_get_num_cols PROTO((sym_environment *env, int *numcols));
int sym_get_num_rows PROTO((sym_environment *env, int *numrows));
int sym_get_num_elements PROTO((sym_environment *env, int *numelems));
int sym_get_col_lower PROTO((sym_environment *env, double *collb));
int sym_get_col_upper PROTO((sym_environment *env, double *colub));
int sym_get_row_sense PROTO((sym_environment *env, char *rowsen));
int sym_get_rhs PROTO((sym_environment *env, double *rowrhs));
int sym_get_matrix PROTO((sym_environment *env, int *nz, int *matbeg, 
			  int *matind, double *matval));
int sym_get_row_range PROTO((sym_environment *env, double *rowrng));
int sym_get_row_lower PROTO((sym_environment *env, double *rowlb));
int sym_get_row_upper PROTO((sym_environment *env, double *rowub));
int sym_get_obj_coeff PROTO((sym_environment *env, double *obj));
int sym_get_obj2_coeff PROTO((sym_environment *env, double *obj2));
int sym_get_obj_sense PROTO((sym_environment *env, int *sense));

int sym_is_continuous PROTO((sym_environment *env, int index, int *value));
int sym_is_binary PROTO((sym_environment *env, int index, int *value));
int sym_is_integer PROTO((sym_environment *env, int index, int *value));

double sym_get_infinity PROTO(());

int sym_get_col_solution PROTO((sym_environment *env, double *colsol));
int sym_get_row_activity PROTO((sym_environment *env, double *rowact));
int sym_get_obj_val PROTO((sym_environment *env, double *objval));
int sym_get_primal_bound PROTO((sym_environment *env, double *ub));
int sym_get_iteration_count PROTO((sym_environment *env, int *numnodes));

int sym_set_obj_coeff PROTO((sym_environment *env, int index, double value));
int sym_set_obj2_coeff PROTO((sym_environment *env, int index, double value));
int sym_set_col_lower PROTO((sym_environment *env, int index, double value));
int sym_set_col_upper PROTO((sym_environment *env, int index, double value));
int sym_set_row_lower PROTO((sym_environment *env, int index, double value));
int sym_set_row_upper PROTO((sym_environment *env, int index, double value));
int sym_set_row_type PROTO((sym_environment *env, int index, char rowsense, 
			     double rowrhs, double rowrng));
int sym_set_obj_sense PROTO((sym_environment *env, int sense));
int sym_set_col_solution PROTO((sym_environment *env, double * colsol));
int sym_set_primal_bound PROTO((sym_environment *env, double bound));
int sym_set_continuous PROTO((sym_environment *env, int index));
int sym_set_integer PROTO((sym_environment *env, int index));
int sym_set_col_names PROTO((sym_environment *env, char **colname));
int sym_add_col PROTO((sym_environment *env, int numelems, int *indices, 
		       double *elements, double collb, double colub,
		       double obj, char is_int, char *name));
int sym_add_row PROTO((sym_environment *env, int numelems, int *indices, 
		       double *elements, char rowsen, double rowrhs,
		       double rowrng));
int sym_delete_cols PROTO((sym_environment *env, int num, int * indices));
int sym_delete_rows PROTO((sym_environment *env, int num, int * indices));

int sym_write_warm_start_desc PROTO((warm_start_desc *ws, char *file));
warm_start_desc *sym_read_warm_start PROTO((char *file));

void sym_delete_warm_start PROTO((warm_start_desc *ws));
warm_start_desc *sym_get_warm_start PROTO((sym_environment *env, 
					   int copy_warm_start));

int sym_set_warm_start PROTO((sym_environment *env, warm_start_desc *ws));

int sym_set_int_param PROTO((sym_environment *env, const char *key, int value));
int sym_set_dbl_param PROTO((sym_environment *env, const char *key, double value));
int sym_set_str_param PROTO((sym_environment *env, const char *key, const char *value));

int sym_get_int_param PROTO((sym_environment *env, const char *key, int *value));
int sym_get_dbl_param PROTO((sym_environment *env, const char *key, double *value));
int sym_get_str_param PROTO((sym_environment *env, const char *key, char **value));
//Anahita
int sym_get_dual_pruned PRROTO((sym_environment *env,
				double** dual_pieces, int* num_pieces,
				int MAX_ALLOWABLE_NUM_PIECES)));

int sym_get_lb_for_new_rhs PROTO((sym_environment *env, int cnt, 
				  int *new_rhs_ind, double *new_rhs_val,
				  double *lb_for_new_rhs));
int sym_get_ub_for_new_rhs PROTO((sym_environment *env, int cnt, 
				  int *new_rhs_ind, double *new_rhs_val,
				  double *ub_for_new_rhs));
#if 0
int sym_get_lb_for_new_obj PROTO((sym_environment *env, int cnt, 
				  int *new_obj_ind, double *new_obj_val,
				  double *lb_for_new_obj));
#endif
int sym_get_ub_for_new_obj PROTO((sym_environment *env, int cnt, 
				  int *new_obj_ind, double *new_obj_val,
				  double *ub_for_new_obj));
				     
warm_start_desc *sym_create_copy_warm_start PROTO((warm_start_desc * ws));
MIPdesc *sym_create_copy_mip_desc PROTO((sym_environment *env));
sym_environment * sym_create_copy_environment PROTO((sym_environment *env));

int sym_test PROTO((sym_environment *env));

#endif
