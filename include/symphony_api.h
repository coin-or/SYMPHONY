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
#include "lp_solver.h"
#include "messages.h"

/*===========================================================================*/
/*===================== Interface functions (master.c) ======================*/
/*===========================================================================*/

problem *sym_open_environment PROTO((void));
int sym_set_defaults PROTO((problem *p));
int sym_parse_command_line PROTO((problem *p, int argc, char **argv));
int sym_load_problem PROTO((problem *p));
int sym_find_initial_bounds PROTO((problem *p));
int sym_solve PROTO((problem *p));
int sym_mc_solve PROTO((problem *p));
int sym_create_permanent_cut_pools PROTO((problem *p));
int sym_close_environment PROTO((problem *p));

#endif
