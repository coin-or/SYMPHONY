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

#ifdef USE_SYM_APPLICATION 
  switch( user_receive_cg_data(&p->user, p->draw_graph) ){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
    case USER_DEFAULT:
      freebuf(r_bufid);
      return(TRUE);
    case USER_ERROR:
    default:
      freebuf(r_bufid);
      /* Unexpected return value. Do something!! */
      return(FALSE);
   }
#else
  freebuf(r_bufid);
  return(TRUE);
#endif
}

/*===========================================================================*/

int receive_lp_solution_cg_u(cg_prob *p)
{
#ifdef USE_SYM_APPLICATION
   return(user_receive_lp_solution_cg(&p->user));
#else
   return(USER_DEFAULT);
#endif
}

/*===========================================================================*/

int find_cuts_u(cg_prob *p, LPdata *lp_data, int *num_cuts)
{
   cut_data **cuts;
   int i;
   int termcode;
   int tmp = p->cuts_to_add_num; 
   
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_find_cuts(p->user, p->cur_sol.xlength,
				      p->cur_sol.xiter_num, p->cur_sol.xlevel,
				      p->cur_sol.xindex, p->cur_sol.objval,
				      p->cur_sol.xind, p->cur_sol.xval,
				      p->ub, p->cur_sol.lpetol,
				      &p->cuts_to_add_num, &p->cuts_to_add_size,
				      &p->cuts_to_add) );
#endif

   *num_cuts += p->cuts_to_add_num - tmp;
   
   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int free_cg_u(cg_prob *p)
{
   int termcode;

#ifdef COMPILE_IN_CG
   FREE(p->cuts_to_add);
#else
   FREE(p->cur_sol.xind);
   FREE(p->cur_sol.xval);
#endif   

#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_free_cg(&p->user) );
#endif
   FREE(p);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

#ifdef CHECK_VALIDITY
void check_validity_of_cut_u(cg_prob *p, cut_data *new_cut)
{
   int termcode;

   switch(new_cut->type){

    case EXPLICIT_ROW:

      /* Not implemented yet */
	 
       return (FUNCTION_TERMINATED_NORMALLY);

    default:
#ifdef USE_APPLICATION
      CALL_USER_FUNCTION( user_check_validity_of_cut(p->user, new_cut) );
#endif
   }
}
#endif
