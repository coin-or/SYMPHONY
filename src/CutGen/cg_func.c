/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>

#include "pack_cut.h"
#include "messages.h"
#include "proccomm.h"
#include "timemeas.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "cg.h"

/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
#   include "decomp.h"
#endif

/*___END_EXPERIMENTAL_SECTION___*/
/*===========================================================================*/

/*===========================================================================*\
 * This file contains general functions used by the cut generator process.
\*===========================================================================*/

/*===========================================================================*/

/*===========================================================================*\
 * This function receives the problem data (if we are running in parallel)   
 * and intitializes the data structures.                                     
\*===========================================================================*/


void cg_initialize(cg_prob *p, int master_tid)
{
#ifndef COMPILE_IN_CG
   int bytes, msgtag;
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   int s_bufid, r_bufid, info;
#endif
#endif
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP)
   int r_bufid, s_bufid = 0;
#endif

   /*------------------------------------------------------------------------
    * Receive the problem data 
    *------------------------------------------------------------------------*/

#ifdef COMPILE_IN_CG

   p->master = master_tid;

#else

   /* We only need to do this part if the CG is running as a separate process*/
   /* Otherwise, all of this setup is done in the master in the function     */
   /* pack_cg_data_u()*/
   
   /* set stdout to be line buffered */
   setvbuf(stdout, (char *)NULL, _IOLBF, 0);
   
   register_process();

   r_bufid = receive_msg(ANYONE, MASTER_TID_INFO);
   bufinfo(r_bufid, &bytes, &msgtag, &p->tree_manager);
   receive_int_array(&p->master, 1);
   receive_int_array(&p->proc_index, 1);
   freebuf(r_bufid);
   
#endif
   
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) || \
   !defined(COMPILE_IN_CG)
      
   /* This part, we only need to do if we are not running in full serial mode*/

   s_bufid = init_send(DataInPlace);
   send_msg(p->master, REQUEST_FOR_CG_DATA);
   freebuf(s_bufid);

   receive_cg_data_u(p);
   
#endif
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   if (p->par.do_decomp)
      /* This doesn't work anymore because varnum is not defined */
      open_decomp_lp(p, varnum);
#endif
/*___END_EXPERIMENTAL_SECTION___*/
   
   (void) used_time(&p->tt);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function is provided for the user to send cuts. This function is
 * retained for backwards compatibility, but is deprecated. See
 * cg_add_user_cut() below.                       
\*===========================================================================*/

int cg_send_cut(cut_data *new_cut, int *num_cuts, int *alloc_cuts,
		cut_data ***cuts)
{
#ifdef COMPILE_IN_CG

   int i;
   cut_data *tmp_cut;

   for (i = 0; i < *num_cuts; i++){
      if (new_cut->type != (*cuts)[i]->type ||
	  new_cut->size != (*cuts)[i]->size ||
	  new_cut->rhs != (*cuts)[i]->rhs){
	 continue;
      }
      if (!new_cut->coef){
	 return(0);
      }
      if (memcmp(new_cut->coef, (*cuts)[i]->coef,
		 new_cut->size) == 0){
	 return(0);
      }
   }
   if (new_cut->name != CUT__DO_NOT_SEND_TO_CP)
      new_cut->name = CUT__SEND_TO_CP;
   tmp_cut = (cut_data *) malloc (sizeof(cut_data));
   memcpy((char *)tmp_cut, (char *)new_cut, sizeof(cut_data));
   if (new_cut->size >0){
      tmp_cut->coef = (char *) malloc (new_cut->size * sizeof(char));
      memcpy((char *)tmp_cut->coef, (char *)new_cut->coef,
	     new_cut->size * sizeof(char));
   }
   REALLOC((*cuts), cut_data *, (*alloc_cuts), (*num_cuts + 1), BB_BUNCH);
   (*cuts)[(*num_cuts)++] = tmp_cut;
   
#else

   int s_bufid;
   
   if (new_cut->name != CUT__DO_NOT_SEND_TO_CP)
      new_cut->name = CUT__SEND_TO_CP;
   s_bufid = init_send(DataInPlace);
   pack_cut(new_cut);
   send_msg(p->cur_sol.lp, PACKED_CUT);
   freebuf(s_bufid);
   
#endif

#ifdef CHECK_CUT_VALIDITY
   check_validity_of_cut_u(p->user, new_cut);
#endif
 
   return(1);
}

/*===========================================================================*/

cut_data *create_explicit_cut(int nzcnt, int *indices, double *values,
			      double rhs, double range, char sense,
			      char send_to_cp)
{
   cut_data *cut = (cut_data *) calloc(1, sizeof(cut_data));

   cut->type = EXPLICIT_ROW;
   cut->sense = sense;
   cut->rhs = rhs;
   cut->range = range;
   cut->size = ISIZE + nzcnt * (ISIZE + DSIZE);
   cut->coef = (char *) malloc (cut->size);
   ((int *) cut->coef)[0] = nzcnt;
   memcpy(cut->coef + ISIZE, (char *)indices, nzcnt*ISIZE);
   memcpy(cut->coef + (nzcnt + 1) * ISIZE, (char *)values, nzcnt * DSIZE);
   cut->branch = DO_NOT_BRANCH_ON_THIS_ROW;
   cut->deletable = TRUE;
   cut->name = send_to_cp ? CUT__SEND_TO_CP : CUT__DO_NOT_SEND_TO_CP;

   return(cut);
}

/*===========================================================================*/

int cg_add_explicit_cut(int nzcnt, int *indices, double *values,
			double rhs, double range, char sense,
			char send_to_cp, int *num_cuts, int *alloc_cuts,
			cut_data ***cuts)
{
   cut_data *cut = (cut_data *) calloc(1, sizeof(cut_data));

   cut->type = EXPLICIT_ROW;
   cut->sense = sense;
   cut->rhs = rhs;
   cut->range = range;
   cut->size = ISIZE + nzcnt * (ISIZE + DSIZE);
   cut->coef = (char *) malloc (cut->size);
   ((int *) cut->coef)[0] = nzcnt;
   memcpy(cut->coef + ISIZE, (char *)indices, nzcnt*ISIZE);
   memcpy(cut->coef + (nzcnt + 1) * ISIZE, (char *)values, nzcnt * DSIZE);
   cut->branch = DO_NOT_BRANCH_ON_THIS_ROW;
   cut->deletable = TRUE;
   cut->name = send_to_cp ? CUT__SEND_TO_CP : CUT__DO_NOT_SEND_TO_CP;

   return(cg_add_user_cut(cut, num_cuts, alloc_cuts, cuts));
}

/*===========================================================================*/

int cg_add_user_cut(cut_data *new_cut, int *num_cuts, int *alloc_cuts,
		    cut_data ***cuts)
{
#ifdef COMPILE_IN_CG

   int i;
   cut_data *tmp_cut;

   for (i = 0; i < *num_cuts; i++){
      if (new_cut->size != (*cuts)[i]->size){
	 continue;
      }
      if (memcmp(new_cut->coef, (*cuts)[i]->coef, new_cut->size) == 0){
	 return(0);
      }
   }
   if (new_cut->name != CUT__DO_NOT_SEND_TO_CP)
      new_cut->name = CUT__SEND_TO_CP;
   tmp_cut = (cut_data *) malloc (sizeof(cut_data));
   memcpy((char *)tmp_cut, (char *)new_cut, sizeof(cut_data));
   if (new_cut->size >0){
      tmp_cut->coef = (char *) malloc (new_cut->size * sizeof(char));
      memcpy((char *)tmp_cut->coef, (char *)new_cut->coef,
	     new_cut->size * sizeof(char));
   }
   REALLOC((*cuts), cut_data *, (*alloc_cuts), (*num_cuts + 1), BB_BUNCH);
   (*cuts)[(*num_cuts)++] = tmp_cut;

#else

   int s_bufid;
   
   if (new_cut->name != CUT__DO_NOT_SEND_TO_CP)
      new_cut->name = CUT__SEND_TO_CP;
   s_bufid = init_send(DataInPlace);
   pack_cut(new_cut);
   send_msg(p->cur_sol.lp, PACKED_CUT);
   freebuf(s_bufid);
   
#endif

#ifdef CHECK_CUT_VALIDITY
   check_validity_of_cut_u(p->user, new_cut);
#endif
 
   return(1);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function frees data structures
\*===========================================================================*/

void cg_close(cg_prob *p)
{
   free_cg_u(p);
}
