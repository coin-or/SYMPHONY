/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _API_H
#define _API_H

#define COMPILING_FOR_MASTER

#include "proto.h"
#include "master.h"
#include "messages.h"

/*===========================================================================*/
/*===================== Interface functions (master.c) ======================*/
/*===========================================================================*/

sym_environment *sym_open_environment PROTO((void));
int sym_set_defaults PROTO((sym_environment *env));
int sym_parse_command_line PROTO((sym_environment *env, int argc, 
				  char **argv));
int sym_set_user_data PROTO((sym_environment *env, void *user));
int sym_load_problem PROTO((sym_environment *env));
int sym_find_initial_bounds PROTO((sym_environment *env));
int sym_solve PROTO((sym_environment *env));
int sym_resolve PROTO((sym_environment *env));
int sym_initial_solve PROTO((sym_environment *env));
int sym_mc_solve PROTO((sym_environment *env));
int sym_create_permanent_cut_pools PROTO((sym_environment *env));
int sym_close_environment PROTO((sym_environment *env));
int sym_load_problem_user PROTO((sym_environment *env, int numcols, 
				 int numrows, int *start, int *index, 
				 double *value, double *collb, double *colub, 
				 double *obj, char *rowsen, double *rowrhs, 
				 double *rowrng));   
int sym_is_abandoned PROTO((sym_environment *env));
int sym_is_proven_optimal PROTO((sym_environment *env));
int sym_is_proven_primal_infeasible PROTO((sym_environment *env));	 
int sym_is_primal_objective_limit_reached PROTO((sym_environment *env));
int sym_is_iteration_limit_reached PROTO((sym_environment *env)); 
int sym_is_time_limit_reached PROTO((sym_environment *env));
int sym_is_target_gap_achieved PROTO((sym_environment *env));
int sym_get_num_cols PROTO((sym_environment *env));
int sym_get_num_rows PROTO((sym_environment *env));
int sym_get_num_elements PROTO((sym_environment *env));
double *sym_get_col_lower PROTO((sym_environment *env));
double *sym_get_col_upper PROTO((sym_environment *env));
char *sym_get_row_sense PROTO((sym_environment *env));
double *sym_get_rhs PROTO((sym_environment *env));
double *sym_get_row_range PROTO((sym_environment *env));
double *sym_get_row_lower PROTO((sym_environment *env));
double *sym_get_row_upper PROTO((sym_environment *env));
double *sym_get_obj_coeff PROTO((sym_environment *env));
int sym_get_obj_sense PROTO((sym_environment *env));
int sym_is_continuous PROTO((sym_environment *env, int index));
int sym_is_binary PROTO((sym_environment *env, int index));
int sym_is_integer PROTO((sym_environment *env, int index));
double sym_get_infinity PROTO(());
double *sym_get_col_solution PROTO((sym_environment *env));
double *sym_get_row_activity PROTO((sym_environment *env));
double sym_get_obj_val PROTO((sym_environment *env));
int sym_get_iteration_count PROTO((sym_environment *env));
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
int sym_set_continuous PROTO((sym_environment *env, int index));
int sym_set_integer PROTO((sym_environment *env, int index));
int sym_set_col_names PROTO((sym_environment *env, char **colname));
int sym_add_col PROTO((sym_environment *env, int num_elements, int *indices, 
			double *elements, double collb, double colub,
			double obj, char *name));
int sym_add_row PROTO((sym_environment *env, int num_elements, int *indices, 
			double *elements, char rowsen, double rowrhs,
			double rowrng));
int sym_delete_cols PROTO((sym_environment *env, int num, int * indices));
int sym_delete_rows PROTO((sym_environment *env, int num, int * indices));
int sym_write_warm_start_desc PROTO((warm_start_desc *ws, char *file));
warm_start_desc *sym_read_warm_start PROTO((char *file));

void sym_trim_tree PROTO((bc_node *node));
void sym_delete_warm_start PROTO((warm_start_desc *ws));

warm_start_desc *sym_get_warm_start PROTO((sym_environment *env, 
						int copy_warm_start));
int sym_set_warm_start PROTO((sym_environment *env, warm_start_desc *ws));
void sym_set_int_param PROTO((sym_environment *env,  char *key, int value));
void sym_set_dbl_param PROTO((sym_environment *env,  char *key, double value));
void sym_set_str_param PROTO((sym_environment *env,  char *key, char *value));
int  sym_get_int_param PROTO((sym_environment *env,  char *key));
double sym_get_dbl_param PROTO((sym_environment *env,  char *key));
char *sym_get_str_param PROTO((sym_environment *env,  char *key));

warm_start_desc *sym_create_copy_warm_start PROTO((warm_start_desc * ws));
MIPdesc *sym_create_copy_mip_desc PROTO((sym_environment *env));
sym_environment * sym_create_copy_environment PROTO((sym_environment *env));

double sym_get_lb_for_new_rhs PROTO((sym_environment *env, int cnt, 
				     int *new_rhs_ind, double *new_rhs_val));
				     
#endif
