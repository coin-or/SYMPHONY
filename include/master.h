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

#ifndef _MASTER_H
#define _MASTER_H

#define COMPILING_FOR_MASTER

#include "symphony_api.h"

/*===========================================================================*/
/*=================== Master I/O functions (readparams.c) ===================*/
/*===========================================================================*/

int parse_command_line PROTO((problem *p, int argc, char **argv));
void read_string PROTO((char *target, char *line, int maxlen));
void print_statistics PROTO((node_times *tim, tm_stat *stat, double ub,
			     double lb, double initial_time,
			     double start_time));

/*===========================================================================*/
/*=============== Master wrapper functions (master_wrapper.c) ===============*/
/*===========================================================================*/

int initialize_u PROTO((problem *p));
int free_master_u PROTO((problem *p));
int readparams_u PROTO((problem *p, int argc, char **argv));
int io_u PROTO((problem *p));
int init_draw_graph_u PROTO((problem *p));
int start_heurs_u PROTO((problem *p));
int display_solution_u PROTO((problem *p, int thread_num));
int initialize_root_node_u PROTO((problem *p));
int receive_feasible_solution_u PROTO((problem *p, int msgtag));
int send_lp_data_u PROTO((problem *p, int sender));
int send_cg_data_u PROTO((problem *p, int sender));
int send_cp_data_u PROTO((problem *p, int sender));
int send_sp_data_u PROTO((problem *p, int sender));
int process_own_messages_u PROTO((problem *p, int msgtag));

#endif
