/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2005-2019 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _SYM_API_H
#define _SYM_API_H

#define COMPILING_FOR_MASTER

#ifdef PROTO
#undef PROTO
#endif
#define PROTO(x) x

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  Return Values                       **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*----------------------- Global return codes -------------------------------*/
#define FUNCTION_TERMINATED_NORMALLY      0
#define FUNCTION_TERMINATED_ABNORMALLY   -1
#define ERROR__USER                      -100

/*-------------- Return codes for sym_parse_comand_line() -------------------*/
#define ERROR__OPENING_PARAM_FILE        -110
#define ERROR__PARSING_PARAM_FILE        -111

/*----------------- Return codes for sym_load_problem() ---------------------*/
#define ERROR__READING_GMPL_FILE         -120
#define ERROR__READING_WARM_START_FILE   -121
#define ERROR__READING_MPS_FILE          -122
#define ERROR__READING_LP_FILE           -123

/*-------------------- Return codes for sym_solve() -------------------------*/
#define TM_NO_PROBLEM                     225
#define TM_NO_SOLUTION                    226
#define TM_OPTIMAL_SOLUTION_FOUND         227
#define TM_TIME_LIMIT_EXCEEDED            228
#define TM_NODE_LIMIT_EXCEEDED            229
#define TM_ITERATION_LIMIT_EXCEEDED       230
#define TM_TARGET_GAP_ACHIEVED            231
#define TM_FOUND_FIRST_FEASIBLE           232
#define TM_FINISHED                       233
#define TM_UNFINISHED                     234
#define TM_FEASIBLE_SOLUTION_FOUND        235
#define TM_SIGNAL_CAUGHT                  236
#define TM_UNBOUNDED                      237
#define PREP_OPTIMAL_SOLUTION_FOUND       238
#define PREP_NO_SOLUTION                  239
#define TM_ERROR__NO_BRANCHING_CANDIDATE -250
#define TM_ERROR__ILLEGAL_RETURN_CODE    -251
#define TM_ERROR__NUMERICAL_INSTABILITY  -252
#define TM_ERROR__COMM_ERROR             -253
#define TM_ERROR__USER                   -275
#define PREP_ERROR                       -276

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  General Constants                   **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef ANYONE
#define ANYONE   -1
#endif
#ifndef ANYTHING
#define ANYTHING -1
#endif

#define DSIZE sizeof(double)
#define ISIZE sizeof(int)
#define CSIZE sizeof(char)

#ifndef BITSPERBYTE
#define	BITSPERBYTE 8
#endif
#ifndef BITS
#define	BITS(type) (BITSPERBYTE * (int)sizeof (type))
#endif

#ifdef  HIBITI
#undef  HIBITI
#endif
#define	HIBITI (1U << (BITS(int) - 1))
#ifdef  MAXINT
#undef  MAXINT
#endif
#define	MAXINT ((int)(~(HIBITI)))
#ifdef MAXDOUBLE
#undef MAXDOUBLE
#endif
#define	MAXDOUBLE 1.79769313486231570e+308

#define SYM_INFINITY                 1e20

#define BIG_DBL                      1e40

#define SYM_MINIMIZE                 0
#define SYM_MAXIMIZE                 1 

#define MAX_NAME_SIZE                255

/*--------------------- return values for user-written functions ------------*/
#define USER_ERROR              -5
#define USER_SUCCESS            -4
#define USER_NO_PP              -3
#define USER_AND_PP             -2
#define USER_DEFAULT            -1

/*------------ search order options for multi-criteria problems -------------*/
#define MC_FIFO                  0
#define MC_LIFO                  1

/*------------  warm_starting options for multi-criteria problems -------------*/
#define MC_WS_UTOPIA_FIRST               0
#define MC_WS_UTOPIA_BOTH_FIXED          1
#define MC_WS_UTOPIA_BOTH                2
#define MC_WS_BEST_CLOSE                 3

/*------------------------ compare_candidates -------------------------------*/
#define BIGGEST_DIFFERENCE_OBJ   0
#define LOWEST_LOW_OBJ           1
#define HIGHEST_LOW_OBJ          2
#define LOWEST_HIGH_OBJ          3
#define HIGHEST_HIGH_OBJ         4
#define HIGH_LOW_COMBINATION     9

/*--------------------------- select_child ----------------------------------*/
#define PREFER_LOWER_OBJ_VALUE   0
#define PREFER_HIGHER_OBJ_VALUE  1

/*-------------------- generate_cuts_in_lp defaults -------------------------*/
#define GENERATE_CGL_CUTS                  20
#define DO_NOT_GENERATE_CGL_CUTS           21

/*-------------------- xxx_cuts_generation_levels ---------------------------*/
#define DO_NOT_GENERATE        -1
#define GENERATE_DEFAULT        0
#define GENERATE_IF_IN_ROOT     1    
#define GENERATE_ONLY_IN_ROOT   2
#define GENERATE_ALWAYS         3 
#define GENERATE_PERIODICALLY   4

/*------------------------- node selection rules ----------------------------*/
#define LOWEST_LP_FIRST       0
#define HIGHEST_LP_FIRST      1
#define BREADTH_FIRST_SEARCH  2
#define DEPTH_FIRST_SEARCH    3
#define BEST_FIRST_SEARCH     4
#define DEPTH_FIRST_THEN_BEST_FIRST 5

/*-------------------------- diving_strategy --------------------------------*/
#define BEST_ESTIMATE         0
#define COMP_BEST_K           1
#define COMP_BEST_K_GAP       2

/*--------------- parameter values for feasibility pump heuristic -----------*/
#define SYM_FEAS_PUMP_DEFAULT    1       /* use fp using the default rules   */
#define SYM_FEAS_PUMP_REPEATED   2       /* use fp till the end of solve     */
#define SYM_FEAS_PUMP_TILL_SOL   3       /* use fp till a solution is found  */
#define SYM_FEAS_PUMP_DISABLE   -1       /* dont use fp */

typedef struct MIPDESC MIPdesc;
typedef struct WARM_START_DESC warm_start_desc;
typedef struct SYM_ENVIRONMENT sym_environment;
typedef struct CUT_POOL cut_pool;

#include "SymConfig.h"

#ifdef __cplusplus
extern "C" // Export C names when in C
{
#endif
   
/*===========================================================================*/
/*===================== Interface functions (master.c) ======================*/
/*===========================================================================*/

SYMPHONYLIB_EXPORT void sym_version(void);
SYMPHONYLIB_EXPORT sym_environment * sym_open_environment(void);
SYMPHONYLIB_EXPORT int sym_close_environment(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_reset_environment(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_set_defaults(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_parse_command_line(sym_environment *env,
                                                  int argc, char **argv);
SYMPHONYLIB_EXPORT int sym_set_user_data(sym_environment *env, void *user);
SYMPHONYLIB_EXPORT int sym_get_user_data(sym_environment *env, void **user);
SYMPHONYLIB_EXPORT int sym_read_mps(sym_environment *env, char *infile);
SYMPHONYLIB_EXPORT int sym_read_lp(sym_environment *env, char *infile);
SYMPHONYLIB_EXPORT int sym_read_gmpl(sym_environment *env, char *modelfile, 
                                         char *datafile);
SYMPHONYLIB_EXPORT int sym_write_mps(sym_environment *env, char *infile);
SYMPHONYLIB_EXPORT int sym_write_lp(sym_environment *env, char *infile);
SYMPHONYLIB_EXPORT int sym_load_problem(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_find_initial_bounds(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_solve(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_warm_solve(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_mc_solve(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_create_permanent_cut_pools(sym_environment *env,
                                                          int *cp_num); 
SYMPHONYLIB_EXPORT cut_pool **sym_get_permanent_cut_pools(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_explicit_load_problem(sym_environment
                                     *env, int numcols, 
				     int numrows, int *start, int *index, 
				     double *value, double *collb,
				     double *colub, char *is_int, double *obj,
				     double *obj2, char *rowsen,
				     double *rowrhs, double *rowrng,
				     char make_copy);   

SYMPHONYLIB_EXPORT int sym_is_abandoned(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_is_proven_optimal(sym_environment *env);
SYMPHONYLIB_EXPORT int
sym_is_proven_primal_infeasible(sym_environment *env);	  
SYMPHONYLIB_EXPORT int 
sym_is_iteration_limit_reached(sym_environment *env); 
SYMPHONYLIB_EXPORT int sym_is_time_limit_reached(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_is_target_gap_achieved(sym_environment *env);

SYMPHONYLIB_EXPORT int sym_get_status(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_get_num_cols(sym_environment *env, int *numcols);
SYMPHONYLIB_EXPORT int sym_get_num_rows(sym_environment *env, int *numrows);
SYMPHONYLIB_EXPORT int sym_get_num_elements(sym_environment *env,
                                                int *numelems);
SYMPHONYLIB_EXPORT int sym_get_col_lower(sym_environment *env,
                                             double *collb);
SYMPHONYLIB_EXPORT int sym_get_col_upper(sym_environment *env,
                                             double *colub);
SYMPHONYLIB_EXPORT int sym_get_row_sense(sym_environment *env,
                                             char *rowsen);
SYMPHONYLIB_EXPORT int sym_get_rhs(sym_environment *env,
                                       double *rowrhs);
SYMPHONYLIB_EXPORT int sym_get_matrix(sym_environment *env, int *nz,
                                          int *matbeg, 
                                          int *matind, double *matval);
SYMPHONYLIB_EXPORT int sym_get_row_range(sym_environment *env,
                                             double *rowrng);
SYMPHONYLIB_EXPORT int sym_get_row_lower(sym_environment *env,
                                             double *rowlb);
SYMPHONYLIB_EXPORT int sym_get_row_upper(sym_environment *env,
                                             double *rowub);
SYMPHONYLIB_EXPORT int sym_get_obj_coeff(sym_environment *env,
                                             double *obj);
SYMPHONYLIB_EXPORT int sym_get_obj2_coeff(sym_environment *env,
                                              double *obj2);
SYMPHONYLIB_EXPORT int sym_get_obj_sense(sym_environment *env, int *sense);
SYMPHONYLIB_EXPORT int sym_is_continuous(sym_environment *env, int index,
                                             int *value);
SYMPHONYLIB_EXPORT int sym_is_binary(sym_environment *env, int index,
                                         int *value);
SYMPHONYLIB_EXPORT int sym_is_integer(sym_environment *env, int index,
                                          char *value);
SYMPHONYLIB_EXPORT double sym_get_infinity();

SYMPHONYLIB_EXPORT int sym_get_col_solution(sym_environment *env,
                                                double *colsol);
SYMPHONYLIB_EXPORT int sym_get_sp_size(sym_environment *env, int *size);
SYMPHONYLIB_EXPORT int sym_get_sp_solution(sym_environment *env, int index,
                                               double *colsol, double *objval);
SYMPHONYLIB_EXPORT int 
sym_get_row_activity(sym_environment *env, double *rowact);
SYMPHONYLIB_EXPORT int sym_get_obj_val(sym_environment *env,
                                           double *objval);
SYMPHONYLIB_EXPORT int sym_get_primal_bound(sym_environment *env,
                                                double *ub);
SYMPHONYLIB_EXPORT int sym_get_iteration_count(sym_environment *env,
                                                   int *numnodes);
SYMPHONYLIB_EXPORT int sym_set_obj_coeff(sym_environment *env,
                                             int index, double value);
SYMPHONYLIB_EXPORT int sym_set_obj2_coeff(sym_environment *env,
                                              int index, double value);
SYMPHONYLIB_EXPORT int sym_set_col_lower(sym_environment *env,
                                             int index, double value);
SYMPHONYLIB_EXPORT int sym_set_col_upper(sym_environment *env,
                                             int index, double value);
SYMPHONYLIB_EXPORT int sym_set_row_lower(sym_environment *env,
                                             int index, double value);
SYMPHONYLIB_EXPORT int sym_set_row_upper(sym_environment *env,
                                             int index, double value);
SYMPHONYLIB_EXPORT int sym_set_row_type(sym_environment *env,
                                            int index, char rowsense, 
                                            double rowrhs, double rowrng);
SYMPHONYLIB_EXPORT int sym_set_obj_sense(sym_environment *env, int sense);
SYMPHONYLIB_EXPORT int sym_set_col_solution(sym_environment *env,
                                                double * colsol);
SYMPHONYLIB_EXPORT int sym_set_primal_bound(sym_environment *env,
                                                double bound);
SYMPHONYLIB_EXPORT int sym_set_continuous(sym_environment *env, int index);
SYMPHONYLIB_EXPORT int sym_set_integer(sym_environment *env, int index);
SYMPHONYLIB_EXPORT int sym_set_col_names(sym_environment *env,
                                             char **colname);
SYMPHONYLIB_EXPORT int sym_add_col(sym_environment *env, int numelems,
                                       int *indices, double *elements,
                                       double collb, double colub,
                                       double obj, char is_int, char *name);
SYMPHONYLIB_EXPORT int sym_add_row(sym_environment *env, int numelems,
                                       int *indices, double *elements,
                                       char rowsen, double rowrhs,
                                       double rowrng);
SYMPHONYLIB_EXPORT int sym_delete_cols(sym_environment *env, int num,
                                           int * indices);
SYMPHONYLIB_EXPORT int sym_delete_rows(sym_environment *env, int num,
                                           int * indices);
SYMPHONYLIB_EXPORT int sym_write_warm_start_desc(warm_start_desc *ws,
                                                     char *file);
SYMPHONYLIB_EXPORT warm_start_desc * sym_read_warm_start(char *file);
SYMPHONYLIB_EXPORT void sym_delete_warm_start(warm_start_desc *ws);
SYMPHONYLIB_EXPORT warm_start_desc *
sym_get_warm_start(sym_environment *env, int copy_warm_start);
SYMPHONYLIB_EXPORT int sym_set_warm_start(sym_environment *env,
                                              warm_start_desc *ws);
SYMPHONYLIB_EXPORT int sym_set_int_param(sym_environment *env,
                                             const char *key, int value);
SYMPHONYLIB_EXPORT int sym_set_dbl_param(sym_environment *env,
                                             const char *key, double value);
SYMPHONYLIB_EXPORT int sym_set_str_param(sym_environment *env,
                                             const char *key,
                                             const char *value);
SYMPHONYLIB_EXPORT int sym_get_int_param(sym_environment *env,
                                             const char *key, int *value);
SYMPHONYLIB_EXPORT int sym_get_dbl_param(sym_environment *env,
                                             const char *key, double *value);
SYMPHONYLIB_EXPORT int sym_get_str_param(sym_environment *env,
                                             const char *key, char **value);
SYMPHONYLIB_EXPORT int sym_get_lb_for_new_rhs(sym_environment *env,
                                              int rhs_cnt, int *new_rhs_ind,
                                              double *new_rhs_val,
                                              int lb_cnt, int *new_lb_ind,
                                              double *new_lb_val,
                                              int ub_cnt, int *new_ub_ind,
                                              double *new_ub_val,
                                              double *lb_for_new_rhs);
SYMPHONYLIB_EXPORT int sym_get_dual_pruned(sym_environment *env,
                                           double ** dual_pieces,
                                           int* num_pieces,
                                           int MAX_ALLOWABLE_NUM_PIECES);
SYMPHONYLIB_EXPORT int sym_get_ub_for_new_rhs(sym_environment *env,
                                                  int cnt, int *new_rhs_ind,
                                                  double *new_rhs_val,
                                                  double *ub_for_new_rhs);
#if 0
SYMPHONYLIB_EXPORT int sym_get_lb_for_new_obj(sym_environment *env,
                                                  int cnt, int *new_obj_ind,
                                                  double *new_obj_val,
                                                  double *lb_for_new_obj);
#endif
SYMPHONYLIB_EXPORT int sym_get_ub_for_new_obj(sym_environment *env,
                                                  int cnt, int *new_obj_ind,
                                                  double *new_obj_val,
                                                  double *ub_for_new_obj);
SYMPHONYLIB_EXPORT warm_start_desc *
sym_create_copy_warm_start(warm_start_desc * ws); 
SYMPHONYLIB_EXPORT MIPdesc * 
sym_create_copy_mip_desc(sym_environment *env);
SYMPHONYLIB_EXPORT MIPdesc * 
sym_get_presolved_mip_desc(sym_environment *env); 
SYMPHONYLIB_EXPORT sym_environment *
sym_create_copy_environment(sym_environment *env);
SYMPHONYLIB_EXPORT int sym_test(sym_environment *env, int argc,
                                    char **argv, int *test_status);
SYMPHONYLIB_EXPORT void sym_print_statistics(sym_environment *env,
                                             double start_time,
                                             double finish_time);
SYMPHONYLIB_EXPORT double sym_wall_clock(double *T);
SYMPHONYLIB_EXPORT int sym_set_param(sym_environment *env, char *line);
SYMPHONYLIB_EXPORT int sym_free_env(sym_environment *env);

#ifdef __cplusplus
} // end extern "C"
#endif

#endif
