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

#ifndef _CUT_GEN_H
#define _CUT_GEN_H

#include "BB_types.h"
#include "cg_params.h"
#include "cg_u.h"
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
#include "lp_solver.h"
#include "decomp_types.h"
#endif

/*===========================================================================*/

#ifdef COMPILE_DECOMP
typedef struct DECOMP_DATA{
   LPdata        *lp_data;
   int            iter_num;
   double        *unbdd_row;
   int            dunbr;
   int            timeout;
}decomp_data;
#endif
/*___END_EXPERIMENTAL_SECTION___*/

/*===========================================================================*/

/*stores the data needed by the cut_generator*/
typedef struct CG_PROB{
   int            proc_index;
   void          *user;
   int            msgtag;
   int            master;
   int            draw_graph;    /* the tid of DrawGraph */
   int            tree_manager;  /* the tid of the tree manager */
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int            sol_pool;      /* the tid of the solution pool */
   /*___END_EXPERIMENTAL_SECTION___*/
   cg_params      par;           /* the parameters for the cut generator */
   char           has_ub;        /* is there an upper bound */
   double         ub;            /* the current best upper bound if there
				    is one */
   double         tt;
   lp_sol         cur_sol;
#ifdef COMPILE_IN_CG
   int           cuts_to_add_num;
   cut_data    **cuts_to_add;
   int           cuts_to_add_size;
#endif
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   decomp_data    dcmp_data;
#endif
/*___END_EXPERIMENTAL_SECTION___*/
}cg_prob;

/*===========================================================================*/
/*==================== CG basic functions (cg_func.c) =======================*/
/*===========================================================================*/

cg_prob *get_cg_ptr PROTO((cg_prob **cg_list));
void cg_initialize PROTO((cg_prob *p, int master_tid));
void cg_close PROTO((cg_prob * p));

/*===========================================================================*/
/*=============== CG communication functions (cg_proccomm.c) ================*/
/*===========================================================================*/

int cg_process_message PROTO((cg_prob *p, int r_bufid));
int cg_send_cut PROTO((cut_data *new_cut));

#endif
