/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2008 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _BB_CONSTANTS_H
#define _BB_CONSTANTS_H

#include "symphony.h"

#define BB_BUNCH 127 * sizeof(double)

/*----------------- Error codes for process_chain() ------------------------*/
#define ERROR__NO_BRANCHING_CANDIDATE  -101
#define ERROR__ILLEGAL_RETURN_CODE     -102
#define ERROR__NUMERICAL_INSTABILITY   -103
#define ERROR__ILLEGAL_BRANCHING       -104
#define ERROR__COMM_ERROR              -105

/*---------------------------- type of the problem --------------------------*/
#define ZERO_ONE_PROBLEM         0
#define INTEGER_PROBLEM          1
#define MIXED_INTEGER_PROBLEM    2

/*---------------------------- input format ---------------------------------*/
#define MPS_FORMAT               0 
#define LP_FORMAT                1 
#define GMPL_FORMAT              2 

/*--------------------------- modes of giving a list ------------------------*/
#define WRT_PARENT               0
#define EXPLICIT_LIST            1
#define NO_DATA_STORED           2

/*----------------- possible stati of a node in the search tree -------------*/
#define NODE_STATUS__CANDIDATE    0
#define NODE_STATUS__BRANCHED_ON  1
#define NODE_STATUS__HELD         2
#define NODE_STATUS__ROOT         3
#define NODE_STATUS__PRUNED       4
#define NODE_STATUS__INTERRUPTED  5
#define NODE_STATUS__WARM_STARTED 6
#define NODE_STATUS__WSPRUNED     7
/*------------------------------ not_fixed stati ----------------------------*/
#define NF_CHECK_ALL             0x00
#define NF_CHECK_AFTER_LAST      0x01
#define NF_CHECK_UNTIL_LAST      0x02
#define NF_CHECK_NOTHING         0x04

/*------------------- options on whether to dive or not ---------------------*/
#define DO_NOT_DIVE              0
#define DO_DIVE                  1
#define CHECK_BEFORE_DIVE        2

/*-------------- possible node types when sending a node to TM --------------*/
#define ROOT_NODE                       0
#define NODE_BRANCHED_ON                1
#define INFEASIBLE_HOLD_FOR_NEXT_PHASE  2
#define OVER_UB_HOLD_FOR_NEXT_PHASE     3
#define INFEASIBLE_PRUNED               4
#define FEASIBLE_PRUNED                 5
#define OVER_UB_PRUNED                  6
#define DISCARDED_NODE                  7
#define INTERRUPTED_NODE                8
#define REPRICED_NODE                   9
#define MC_FEASIBLE_PRUNED              10
/*to be used when warm_started*/
#define PRUNED_HAS_CAN_SOLUTION        11
#define NOT_PRUNED_HAS_CAN_SOLUTION    12

/*------------------- possible node types for VBC Tool ----------------------*/
#define VBC_INTERIOR_NODE       1 /*Dark Red*/ 
#define VBC_PRUNED              2 /*Green*/
#define VBC_ACTIVE_NODE         3 /*White*/
#define VBC_CAND_NODE           4 /*Light Red*/
#define VBC_FEAS_SOL_FOUND      5 /*Blue*/
#define VBC_PRUNED_INFEASIBLE   6 /*color ??, used only when vbc_emulation = 3*/
#define VBC_PRUNED_FATHOMED     7 /*do*/
#define VBC_IGNORE              8 /*do*/


/*------------------ what to do with pruned nodes in the TM -----------------*/
#define DISCARD                  0
#define KEEP_ON_DISK_FULL        1
#define KEEP_ON_DISK_VBC_TOOL    2
#define KEEP_IN_MEMORY           3

/*---------------------- logging options in the TM --------------------------*/
#define NO_LOGGING               0
#define FULL_LOGGING             1
#define VBC_TOOL                 2

/*-------------------- emulation options for vbctool ------------------------*/
#define NO_VBC_EMULATION         0
#define VBC_EMULATION_FILE       1
#define VBC_EMULATION_LIVE       2
#define VBC_EMULATION_FILE_NEW   3

/*------------ search order options for multi-criteria problems -------------*/
#define MC_FIFO                  0
#define MC_LIFO                  1

/*------------  warm_starting options for multi-criteria problems -------------*/
#define MC_WS_UTOPIA_FIRST               0
#define MC_WS_UTOPIA_BOTH_FIXED          1
#define MC_WS_UTOPIA_BOTH                2
#define MC_WS_BEST_CLOSE                 3

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************       Constants related to the LP process            **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*****************************************************************************/
/*******         Constants related to solving of an lp                 *******/
/*****************************************************************************/

/*--------------------- Constants related to the basis ----------------------*/
#  define INVALID_BASIS_STATUS 255
#  define VAR_AT_LB    0
#  define VAR_BASIC    1
#  define VAR_AT_UB    2
#  define VAR_FREE     3
#  define VAR_FIXED    4   /* OSLLIB uses this */
#  define SLACK_AT_LB  0
#  define SLACK_BASIC  1
#  define SLACK_AT_UB  2
#  define SLACK_FREE   3
#  define SLACK_FIXED  4   /* OSLLIB uses this for 'E' rows */

/*----------------- LP Solver return codes (dual simplex) -------------------*/
#define LP_OPTIMAL                   0
#define LP_D_INFEASIBLE              1
#define LP_D_UNBOUNDED               2
#define LP_D_ITLIM                   3
#define LP_D_OBJLIM                  4
#define LP_OPT_FEASIBLE              5
#define LP_OPT_FEASIBLE_BUT_CONTINUE 6
#define LP_ABANDONED                 7

#define MOVE_TO_LB               0
#define MOVE_TO_UB               1

#define LOWER_THAN_LB            0 
#define HIGHER_THAN_UB           1

/*----------------------------- LP status -----------------------------------*/
#define LP_HAS_BEEN_ABANDONED    0
#define LP_HAS_NOT_BEEN_MODIFIED 1
#define LP_HAS_BEEN_MODIFIED     2

/*--------------------- how the variables are ordered -----------------------*/
#define COLIND_ORDERED             0
#define USERIND_ORDERED            1
#define COLIND_AND_USERIND_ORDERED 2

/*---------------- The possible total dual feasibility stati ----------------*/
#define NOT_TDF     0
#define TDF_NOT_ALL 1
#define TDF_HAS_ALL 2

/*----- for marking an lp_cut that it was copied to slack_cuts already ------*/
#define DO_NOT_BRANCH_ON_THIS_ROW  0x01
#define ALLOWED_TO_BRANCH_ON       0x02
#define CANDIDATE_FOR_BRANCH       0x04
#define SWITCH_CANDIDATE_ALLOWED   0x06
#define CUT_BRANCHED_ON            0x08

/*------ which constraints considered ineffective (ineffective_constraints) -*/
#define NO_CONSTRAINT_IS_INEFFECTIVE      0
#define NONZERO_SLACKS_ARE_INEFFECTIVE    1
#define BASIC_SLACKS_ARE_INEFFECTIVE      2
#define ZERO_DUAL_VALUES_ARE_INEFFECTIVE  3

/*------ whether a cut should be sent to CP if effective for long enough ----*/
#define CUT__DO_NOT_SEND_TO_CP     -1
#define CUT__SEND_TO_CP            -2
#define BASE_CONSTRAINT            -3

/*---------------- source of cut --------------------------------------------*/
#define LEFTOVER                    0
#define INTERNAL_CUT_POOL          -1
#define EXTERNAL_CUT_POOL          -2
#define INTERNAL_CUT_GEN           -3
#define EXTERNAL_CUT_GEN           -4

/*---------------- what is the status of a particular constraint ------------*/
#define VIOLATED_ROW                0
#define TIGHT_ROW                   1
#define SLACK_ROW                   2

/*--------------------- flags in the lp_data->status ------------------------*/
/************* variables NOT_REMOVABLE always remain in the LP ***************/

#define NOT_FIXED                   0x01

#define TEMP_FIXED_TO_LB            0x02
#define PERM_FIXED_TO_LB            0x04
#define TEMP_PERM_LB__SWITCH        0x06
#define NOT_FIXED__TEMP_LB__SWITCH  0x03
#define NOT_FIXED__PERM_LB__SWITCH  0x05 

#define TEMP_FIXED_TO_UB            0x08
#define PERM_FIXED_TO_UB            0x10
#define TEMP_PERM_UB__SWITCH        0x18
#define NOT_FIXED__TEMP_UB__SWITCH  0x09
#define NOT_FIXED__PERM_UB__SWITCH  0x11

#define BASE_VARIABLE               0x20
#define VARIABLE_BRANCHED_ON        0x40

#define NOT_REMOVABLE               0x60

/*----------------------- col_gen -------------------------------------------*
 * last 2 digits refer to what to do when we would fathom:
 *     DO_NOT_GENERATE__COLS_DISCARD
 *     DO_NOT_GENERATE_COLS_SEND (for next phase)
 *     GENERATE_COLS__RESOLVE
 * next 2 digits refer to what to do right before branching
 ----------------------------------------------------------------------------*/
#define FATHOM__DO_NOT_GENERATE_COLS__DISCARD    0x00
#define FATHOM__DO_NOT_GENERATE_COLS__SEND       0x01
#define FATHOM__GENERATE_COLS__RESOLVE           0x02
#define COLGEN__FATHOM                           0x03

#define BEFORE_BRANCH__DO_NOT_GENERATE_COLS      0x04
#define BEFORE_BRANCH__GENERATE_COLS__RESOLVE    0x08

#define COLGEN_REPRICING                         0x10

/*-------------------- argument for generate_column_u -----------------------*/
#define GENERATE_NEXTIND            0
#define GENERATE_REAL_NEXTIND       1

/*---------------- Where the unpack_cuts routine is called from -------------*/
#define CUT_FROM_CG                 0
#define CUT_FROM_CP                 1
#define CUT_FROM_TM                 2
#define CUT_LEFTOVER                3
#define CUT_NOT_IN_MATRIX_SLACK     5
#if 0
#define CUT_VIOLATED_SLACK          4
#endif

/*--------------------- how should a cut be unpacked ------------------------*/
#define UNPACK_CUTS_MULTIPLE        0
#define UNPACK_CUTS_SINGLE          1

/*------------------------- built-in cut types ------------------------------*/
#define EXPLICIT_ROW                100
#define OPTIMALITY_CUT_FIRST        101
#define OPTIMALITY_CUT_SECOND       102
#define ORIGINAL_CONSTRAINT         103

/*----------------- possible types of candidate objects ---------------------*/
#define CANDIDATE_VARIABLE          0
#define CANDIDATE_CUT_IN_MATRIX     1
#define CANDIDATE_CUT_NOT_IN_MATRIX 2
#define VIOLATED_SLACK              3
#define SLACK_TO_BE_DISCARDED       4

/*----------------- possible types of branching objects ---------------------*/
#define BRANCHING_VARIABLE          0
#define BRANCHING_CUT               1

/*------------ possible return values of select_candidates_u() --------------*/
#define DO_BRANCH                   0
#define DO_NOT_BRANCH               1
#define DO_NOT_BRANCH__FATHOMED     2

/*---------------- possible return values of branch() -----------------------*/
#define NEW_NODE                     -1
#define FATHOMED_NODE                -2

/*------------- normal return value of various functions --------------------*/

/*---------------------------------------------------------------------------*\
 * anything non-negative means continue the node and the value is the
 * number of cuts added. if there were vars added then that fact has already
 * been printed out
\*---------------------------------------------------------------------------*/

/*----------------------- discard_slack_cuts --------------------------------*/
#define DISCARD_SLACKS_BEFORE_NEW_ITERATION   0
#define DISCARD_SLACKS_WHEN_STARTING_NEW_NODE 1

/*---------------- What is the solution to be displayed ---------------------*/
#define DISP_FEAS_SOLUTION          0
#define DISP_RELAXED_SOLUTION       1
#define DISP_FINAL_RELAXED_SOLUTION 2

/*****************************************************************************/
/************ Default options and results of user defined functions **********/
/*****************************************************************************/

/*---------------------------- colgen_strat ---------------------------------*/
#define COLGEN_STR_SIZE 5
#define COLGEN_STR_ARRAY {				\
      { "FATHOM__DO_NOT_GENERATE_COLS__DISCARD",	\
	   FATHOM__DO_NOT_GENERATE_COLS__DISCARD },	\
      { "FATHOM__DO_NOT_GENERATE_COLS__SEND",		\
	   FATHOM__DO_NOT_GENERATE_COLS__SEND    },	\
      { "FATHOM__GENERATE_COLS__RESOLVE",		\
	   FATHOM__GENERATE_COLS__RESOLVE        },	\
      { "BEFORE_BRANCH__DO_NOT_GENERATE_COLS",		\
	   BEFORE_BRANCH__DO_NOT_GENERATE_COLS   },	\
      { "BEFORE_BRANCH__GENERATE_COLS__RESOLVE",	\
	   BEFORE_BRANCH__GENERATE_COLS__RESOLVE }	\
}

/*------------------------- candidate selection -----------------------------*/
#ifdef COMPILE_FRAC_BRANCHING
#define COMPARE_CAN_STR_SIZE 9
#define COMPARE_CAN_STR_ARRAY {				\
   { "LOWEST_LOW_FRAC", LOWEST_LOW_FRAC },			\
   { "HIGHEST_LOW_FRAC", HIGHEST_LOW_FRAC },			\
   { "LOWEST_HIGH_FRAC", LOWEST_HIGH_FRAC },			\
   { "HIGHEST_HIGH_FRAC", HIGHEST_HIGH_FRAC },			\
   { "BIGGEST_DIFFERENCE_OBJ", BIGGEST_DIFFERENCE_OBJ },	\
   { "LOWEST_LOW_OBJ", LOWEST_LOW_OBJ },			\
   { "HIGHEST_LOW_OBJ", HIGHEST_LOW_OBJ },			\
   { "LOWEST_HIGH_OBJ", LOWEST_HIGH_OBJ },			\
   { "HIGHEST_HIGH_OBJ", HIGHEST_HIGH_OBJ }			\
}
#else
#define COMPARE_CAN_STR_SIZE 5
#define COMPARE_CAN_STR_ARRAY {				\
   { "BIGGEST_DIFFERENCE_OBJ", BIGGEST_DIFFERENCE_OBJ },	\
   { "LOWEST_LOW_OBJ", LOWEST_LOW_OBJ },			\
   { "HIGHEST_LOW_OBJ", HIGHEST_LOW_OBJ },			\
   { "LOWEST_HIGH_OBJ", LOWEST_HIGH_OBJ },			\
   { "HIGHEST_HIGH_OBJ", HIGHEST_HIGH_OBJ }			\
}
#endif

/*---------------------------- same_cuts ------------------------------------*/
/* only default is to compare byte by byte */

#define DIFFERENT_CUTS           1
#define SAME_CUTS                2
#define FIRST_CUT_BETTER         3
#define SECOND_CUT_BETTER        4

/*--------------------------- is_feasible -----------------------------------*/
#define TEST_ZERO_ONE            0
#define TEST_INTEGRALITY         1

#define IP_INFEASIBLE               0
#define IP_FEASIBLE                 1
#define IP_FEASIBLE_BUT_CONTINUE    2
#define IP_ALMOST_FEASIBLE          3
#define IP_FEASIBILITY_NOT_KNOWN    4
#define IP_HEUR_FEASIBLE            5

/*------------------------- display_solution --------------------------------*/
#define DISP_NOTHING             0
#define DISP_NZ_INT              1 
#define DISP_NZ_HEXA             2
#define DISP_FRAC_INT            3 
#define DISP_FRAC_HEXA           4
/* no result */

/*--------------------- pack_feasible_solution ------------------------------*/
/* #define SEND_NONZEROS  0*/
/* no result */

/*------------------------- pack_lp_solution --------------------------------*/
#define SEND_NONZEROS            0
#define SEND_FRACTIONS           1

#define LP_SOL_TO_CG             0
#define LP_SOL_TO_CP             1
#define LP_SOL_WITHIN_LP         2
/* no result */

/*------------------------ compare_candidates -------------------------------*/
#define BIGGEST_DIFFERENCE_OBJ   0
#define LOWEST_LOW_OBJ           1
#define HIGHEST_LOW_OBJ          2
#define LOWEST_HIGH_OBJ          3
#define HIGHEST_HIGH_OBJ         4
#ifdef COMPILE_FRAC_BRANCHING
#define HIGHEST_LOW_FRAC         5
#define LOWEST_LOW_FRAC          6
#define HIGHEST_HIGH_FRAC        7
#define LOWEST_HIGH_FRAC         8
#endif

#define FIRST_CANDIDATE_BETTER                    0
#define FIRST_CANDIDATE_BETTER_AND_BRANCH_ON_IT   1
#define SECOND_CANDIDATE_BETTER                   2
#define SECOND_CANDIDATE_BETTER_AND_BRANCH_ON_IT  3
#define BRANCH_ON_IT                              1

/*--------------------------- select_child ----------------------------------*/
#define PREFER_LOWER_OBJ_VALUE   0
#define PREFER_HIGHER_OBJ_VALUE  1
#ifdef COMPILE_FRAC_BRANCHING
#define PREFER_MORE_FRACTIONAL   2
#define PREFER_LESS_FRACTIONAL   3
#endif

#define PRUNE_THIS_CHILD            0
#define RETURN_THIS_CHILD           1
#define KEEP_THIS_CHILD             2
#define PRUNE_THIS_CHILD_FATHOMABLE 3

/*to be used to differentiate the fathomed nodes */
#define PRUNE_THIS_CHILD_INFEASIBLE 4

/*--------------------- shall_we_branch defaults ----------------------------*/
#define USER__DO_NOT_BRANCH      0
#define USER__DO_BRANCH          1
#define USER__BRANCH_IF_MUST     2 /*default*/
#define USER__BRANCH_IF_TAILOFF  3

/*--------------------- select_candidates defaults --------------------------*/
#define USER__CLOSE_TO_HALF                10
#define USER__CLOSE_TO_HALF_AND_EXPENSIVE  11
#define USER__CLOSE_TO_ONE_AND_CHEAP       12

/*-------------------- generate_cuts_in_lp defaults -------------------------*/
#define GENERATE_CGL_CUTS                  20
#define DO_NOT_GENERATE_CGL_CUTS           21

/*-------------------- xxx_cuts_generation_levels ---------------------------*/
#define DO_NOT_GENERATE        -1
#define GENERATE_DEFAULT        0
#define GENERATE_IF_IN_ROOT     1    /*default*/
#define GENERATE_ONLY_IN_ROOT   2
#define GENERATE_ALWAYS         3 
#define GENERATE_PERIODICALLY   4

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************       Constants related to the Tree Manager          **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*--------------------------- tree manager ----------------------------------*/
#define NEW_NODE__NONE          -1
#define NEW_NODE__STARTED       -2
#define NEW_NODE__ERROR         -3

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

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************        Constants related to the Cut Pool             **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*----------------------- cut pool warm start -------------------------------*/
#define NO_WARM_START            0
#define READ_CP_LIST             1 
#define READ_TM_LIST             2

/*--------------- parameter values for "check_which_cuts" -------------------*/
#define CHECK_ALL_CUTS           0
#define CHECK_LEVEL              1
#define CHECK_TOUCHES            2
#define CHECK_LEVEL_AND_TOUCHES  3

/*--------------- parameter values for "delete_which_cuts" ------------------*/
#define DELETE_BY_QUALITY                 1
#define DELETE_BY_TOUCHES                 2
/*These are for backward compatibility -- they default to the old style*/
#define DELETE_DUPLICATES                 2
#define DELETE_DUPLICATE_AND_INEFFECTIVE  2



/*--------------- parameter values for restart/sens analysis ----------------*/

#define NOTHING_CHANGED                   0
#define RHS_CHANGED                       1
#define OBJ_COEFF_CHANGED                 2
#define CONSTRAINT_MATRIX_CHANGED         3 
#define COL_BOUNDS_CHANGED                4  
#define OBJ_SENSE_CHANGED                 5
#define RHS_SENSE_CHANGED                 6 
#define COLS_ADDED                        7

/*--------------- parameter values for restart/sens analysis ----------------*/
#define DO_NOT_TRIM          0
#define TRIM_LEVEL           1
#define TRIM_INDEX           2
#define ON_CRU_VARS          3

/*--------------- parameter values for primal heuristics --------------------*/
#define SYM_FEAS_PUMP_DEFAULT    1       /* use fp using the default rules   */
#define SYM_FEAS_PUMP_REPEATED   2       /* use fp till the end of solve     */
#define SYM_FEAS_PUMP_TILL_SOL   3       /* use fp till a solution is found  */
#define SYM_FEAS_PUMP_DISABLE   -1       /* dont use fp */

#endif
