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

#ifndef COMPILE_IN_LP 

#include <stdlib.h> 
#include <math.h>
#include <memory.h>
#include <malloc.h>

#include "lp.h"
#include "proccomm.h"
#include "messages.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "BB_types.h"
#include "lp_solver.h"

/*===========================================================================*/

/*===========================================================================*\
 * This is the main() that is used if the LP is running as a separate        
 * process. This file is only used in that case.                             
\*===========================================================================*/

int main(void)
{
   lp_prob *p;
   int r_bufid;
   double timeout = 10, time, diff;
   double start_time;
   first_node_rec = FALSE;
   
   start_time = wall_clock(NULL);
   
   p = (lp_prob *) calloc(1, sizeof(lp_prob));

   get_lp_ptr(&p);
   
   lp_initialize(p, 0);
   
   /*------------------------------------------------------------------------*\
    * Continue receiving node data and fathoming branches until this
    * process is killed
   \*------------------------------------------------------------------------*/

   p->phase = 0;
   while (TRUE){
      p->lp_data->col_set_changed = TRUE;
      /*---------------------------------------------------------------------*\
       * waits for an active node message but if there's anything left after
       * receiving that, those messages are processed, before going to
       * process_chain().
      \*---------------------------------------------------------------------*/
      time = wall_clock(NULL);
      do{
	 r_bufid = treceive_msg(ANYONE, ANYTHING, &timeout);
      }while (! process_message(p, r_bufid, NULL, NULL) );
      diff = wall_clock(NULL) - time;
      if (first_node_rec){
	 p->comp_times.idle_node += diff;
      }else{
	 first_node_rec = TRUE;
	 p->comp_times.ramp_up_lp += diff;
      }	 
      do{
	 r_bufid = nreceive_msg(ANYONE, ANYTHING);
	 if (r_bufid)
	    process_message(p, r_bufid, NULL, NULL);
      }while (r_bufid);

      p->comp_times.communication += used_time(&p->tt);

      if (!process_chain(p)){
	 printf("\nUser had problems creating LP! Exiting now.\n\n");
	 /* User had problems creating initial LP. Abandon node. */
	 exit(-1);
      }
   }

   p->comp_time.wall_clock_lp = wall_clock(NULL) - start_time;
   
   lp_exit(p);

   return(0);
}

#endif
