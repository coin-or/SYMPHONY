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

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "BB_types.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "proccomm.h"
#include "messages.h"
#include "cg.h"
#include "cg_u.h"
#ifdef USE_CGL_CUTS
#include "lp_solver.h" /* For CGL cut generation */
#include "lp.h" /* for free_cuts() */
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains CG wrapper functions that interface with the user.
\*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_receive_cg_data that
 * receives the initial data from the Master process. Returns TRUE if
 * succeeded, FALSE otherwise.
\*===========================================================================*/

int receive_cg_data_u(cg_prob *p)
{
   int r_bufid;
   
   r_bufid = receive_msg(p->master, CG_DATA);
   receive_char_array((char *)&p->par, sizeof(cg_params));
   receive_int_array(&p->draw_graph, 1);
   switch( user_receive_cg_data(&p->user, p->draw_graph) ){
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      freebuf(r_bufid);
      break;
    case ERROR:
    default:
      freebuf(r_bufid);
      /* Unexpected return value. Do something!! */
      return(FALSE);
   }
   return(TRUE);
}

/*===========================================================================*/

int receive_lp_solution_cg_u(cg_prob *p)
{
   return(user_receive_lp_solution_cg(&p->user));
}

/*===========================================================================*/

void free_cg_u(cg_prob *p)
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
   CALL_USER_FUNCTION( user_free_cg(&p->user) );
   FREE(p);
}

/*===========================================================================*/

void find_cuts_u(cg_prob *p, LPdata *lp_data, int *num_cuts)
{
   cut_data **cuts;
   int i, explicit_num_cuts;
   
/*__BEGIN_EXPERIMENTAL_SECTION__*/
   CALL_USER_FUNCTION( user_find_cuts(p->user, p->cur_sol.xlength,
				      p->cur_sol.xiter_num, p->cur_sol.xlevel,
				      p->cur_sol.xindex, p->cur_sol.objval,
				      p->cur_sol.xind, p->cur_sol.xval, p->ub,
				      p->cur_sol.lpetol, num_cuts, NULL) );
/*___END_EXPERIMENTAL_SECTION___*/
/*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
   CALL_USER_FUNCTION( user_find_cuts(p->user, p->cur_sol.xlength,
				      p->cur_sol.xiter_num, p->cur_sol.xlevel,
				      p->cur_sol.xindex, p->cur_sol.objval,
				      p->cur_sol.xind, p->cur_sol.xval,
				      p->ub, p->cur_sol.lpetol, num_cuts) );
#endif

#ifdef USE_CGL_CUTS
   
   if ((explicit_num_cuts = generate_cgl_cuts(lp_data, &cuts)) > 0){
      for (i = 0; i < explicit_num_cuts; i++){
	 cg_send_cut(cuts[i]);
      }
      
      free_cuts(cuts, explicit_num_cuts);
   }

#else
   
   explicit_num_cuts = 0;

#endif

   *num_cuts += explicit_num_cuts;
   return;
}

/*===========================================================================*/

#ifdef CHECK_VALIDITY
void check_validity_of_cut_u(cg_prob *p, cut_data *new_cut)
{
   switch(new_cut->type){

    case EXPLICIT_ROW:

      /* Not implemented yet */
	 
       return;

    default:

      CALL_USER_FUNCTION( user_check_validity_of_cut(p->user, new_cut) );
      return;
   }
}
#endif
