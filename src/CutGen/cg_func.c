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

/*===========================================================================*\
 * This small function simply returns a pointer to the current cut generator 
 * problem structure so that it can be accessed from anywhere within the CG  
 * without explicitly passing a pointer. This is essentially only uselful    
 * to the more advanced user sho wishes to have access to the internal data  
 * structures                                                                
\*===========================================================================*/

cg_prob *get_cg_ptr(cg_prob **cg_list)
{
#ifdef _OPENMP
   int thread_num = omp_get_thread_num();
#else
   int thread_num = 0;
#endif
   static cg_prob **cg;

   if (cg_list){
      cg = cg_list;
   }

   return(cg[thread_num]);
}

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
/*__BEGIN_EXPERIMENTAL_SECTION__*/
   int varnum = 0;
/*___END_EXPERIMENTAL_SECTION___*/

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

   r_bufid = receive_msg(p->master, CG_DATA);
   receive_char_array((char *)&p->par, sizeof(cg_params));
   receive_int_array(&p->draw_graph, 1);
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   CALL_USER_FUNCTION( user_receive_cg_data(&p->user, p->draw_graph, &varnum));
   /*___END_EXPERIMENTAL_SECTION___*/
   /*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
   CALL_USER_FUNCTION( user_receive_cg_data(&p->user, p->draw_graph) );
#endif
   freebuf(r_bufid);
#endif
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   if (p->par.do_decomp)
      open_decomp_lp(p, varnum);
#endif
/*___END_EXPERIMENTAL_SECTION___*/
   
   (void) used_time(&p->tt);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function is provided for the user to send cuts                       
\*===========================================================================*/

int cg_send_cut(cut_data *new_cut)
{
   cg_prob *p = get_cg_ptr(NULL);

#ifdef COMPILE_IN_CG

   cut_data *tmp_cut;
     
   if (new_cut->name != CUT__DO_NOT_SEND_TO_CP)
      new_cut->name = CUT__SEND_TO_CP;
   tmp_cut = (cut_data *) malloc (sizeof(cut_data));
   memcpy((char *)tmp_cut, (char *)new_cut, sizeof(cut_data));
   if (new_cut->size >0){
      tmp_cut->coef = (char *) malloc (new_cut->size * sizeof(char));
      memcpy((char *)tmp_cut->coef, (char *)new_cut->coef,
	     new_cut->size * sizeof(char));
   }
   REALLOC(p->cuts_to_add, cut_data *, p->cuts_to_add_size,
	   p->cuts_to_add_num + 1, BB_BUNCH);
   p->cuts_to_add[p->cuts_to_add_num++] = tmp_cut;

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
   user_check_validity_of_cut(p->user, new_cut);
#endif
 
   return(1);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function frees data structures
\*===========================================================================*/

void cg_close(cg_prob *p)
{
#ifdef COMPILE_IN_CG
   FREE(p->cuts_to_add);
#else
   FREE(p->cur_sol.xind);
   FREE(p->cur_sol.xval);
#endif   
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   if (p->par.do_decomp)
      close_decomp_lp(p);
#endif
/*___END_EXPERIMENTAL_SECTION___*/
   user_free_cg(&p->user);
   FREE(p);
}
