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

#ifndef _MASTER_PARAMS_H
#define _MASTER_PARAMS_H

#include "tm_params.h"
#include "cp_params.h"
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
#include "sp_params.h"
#endif
/*___END_EXPERIMENTAL_SECTION___*/
#include "cg_params.h"
#include "lp_params.h"
#include "dg_params.h"

/*===========================================================================*\
 * The params structure contains all of the user-specified parameters
 * to be read in from the parameter file. See the README file for an
 * explanation of the parameters
\*===========================================================================*/

typedef struct PARAMS{
   int        warm_start;
   int        verbosity;
   char       param_file[MAX_FILE_NAME_LENGTH +1];
   int        random_seed;
   cp_params  cp_par;
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   sp_params  sp_par;
#endif
/*___END_EXPERIMENTAL_SECTION___*/
   cg_params  cg_par;
   lp_params  lp_par;
   tm_params  tm_par;
   dg_params  dg_par;
   char       tm_exe[MAX_FILE_NAME_LENGTH +1];
   char       dg_exe[MAX_FILE_NAME_LENGTH +1];
   int        tm_debug;
   int        dg_debug;
   int        tm_machine_set;
   char       tm_machine[MACH_NAME_LENGTH +1];
   int        dg_machine_set;
   char       dg_machine[MACH_NAME_LENGTH +1];
   int        pvm_trace;
   int        do_branch_and_cut;
   int        do_draw_graph;
   char       infile[MAX_FILE_NAME_LENGTH +1]; /* For MPS file name
						  or GNUMP modelfile */
   char       datafile[MAX_FILE_NAME_LENGTH +1]; /* GNUMP datafile */
}params;

#endif
