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

#define COMPILING_FOR_MASTER

#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __PVM__
#include <pvmtev.h>
#endif

#include "proccomm.h"
#include "timemeas.h"
#include "messages.h"
#include "BB_types.h"
#include "BB_macros.h"
#include "pack_cut.h"
#include "pack_array.h"
#include "master.h"

#ifndef TEV_INIT_MASK
/* We must have pvm3.4 where it is called TEV_MASK_INIT */
#  define TEV_INIT_MASK(m)  TEV_MASK_INIT(m)
#  define TEV_SET_MASK(m,k)  TEV_MASK_SET(m,k)
#  define TEV_MCAST0  TEV_MCAST
#  define TEV_RECV0   TEV_RECV
#  define TEV_SEND0   TEV_SEND
#  define TEV_NRECV0  TEV_NRECV
#endif

/*===========================================================================*/

problem *get_problem_ptr()
{
   static problem* p;

   if (!p) p = (problem *) calloc(1, sizeof(problem));

   return(p);
}

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the master process.
\*===========================================================================*/

int main(int argc, char **argv)
{
   problem *p;
   int s_bufid, r_bufid, bytes, msgtag = 0, sender, termcode;
   
   node_desc *root= NULL;
   base_desc *base = NULL;

   double t = 0, total_time=0;
   double start_time, lb;
   struct timeval timeout = {10, 0};

   char lp_data_sent = FALSE, cg_data_sent = FALSE, cp_data_sent = FALSE;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   char sp_data_sent = TRUE; /*for now, we are not using this one*/
   /*___END_EXPERIMENTAL_SECTION___*/

#if (!defined(COMPILE_IN_LP) || !defined(COMPILE_IN_CG) || \
   !defined(COMPILE_IN_CP)) && defined(__PVM__)
   int xpvm_tid;
   Pvmtmask trace_mask;
#endif
#ifndef COMPILE_IN_TM
   char repricing, node_type;
#else
   tm_prob *tm;
#endif
   
   /*------------------------------------------------------------------------*\
    *                         program starts                                 
   \*------------------------------------------------------------------------*/
   (void) used_time(&t);
   start_time = wall_clock(NULL);

   setvbuf(stdout, (char *)NULL, _IOLBF, 0);

   printf("\n");
   printf("*******************************************************\n");
   printf("*   This is SYMPHONY Version 4.0                      *\n");
   printf("*   Copyright 2000-2003 Ted Ralphs                    *\n");
   printf("*   All Rights Reserved                               *\n");
   printf("*   Distributed under the Common Public License 1.0   *\n");
   printf("*******************************************************\n");
   printf("\n");
   
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)
       
   register_process();   /* Enroll this process */

#ifdef __PVM__
   pvm_catchout(stdout); /* Tells PVM to treat all output from the children of
			    this process as output from this process itself*/
#endif
#endif
   
   p = get_problem_ptr();

   initialize_u(p);

   /* Set the parameters */
   readparams_u(p, argc, argv);

   /* This next set of commands has to be executed if we want to create a PVM
      trace file for viewing in xpvm (this is a very slow process) */

#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
   if (p->par.pvm_trace){
      if ((xpvm_tid = pvm_gettid((char *)"xpvm", 0)) > 0){
	 pvm_setopt(PvmSelfTraceTid, xpvm_tid);
	 pvm_setopt(PvmSelfTraceCode, 666);
	 pvm_setopt(PvmSelfOutputTid, xpvm_tid);
	 pvm_setopt(PvmSelfOutputCode, 667);
	 pvm_setopt(PvmTraceTid, xpvm_tid);
	 pvm_setopt(PvmTraceCode, 666);
	 pvm_setopt(PvmOutputTid, xpvm_tid);
	 pvm_setopt(PvmOutputCode, 667);
	 TEV_INIT_MASK(trace_mask);
	 TEV_SET_MASK(trace_mask, TEV_MCAST0);
	 TEV_SET_MASK(trace_mask, TEV_RECV0);
	 TEV_SET_MASK(trace_mask, TEV_SEND0);
	 TEV_SET_MASK(trace_mask, TEV_NRECV0);
	 pvm_settmask(PvmTaskSelf, trace_mask);
	 pvm_settmask(PvmTaskChild, trace_mask);
      }else{
	 PVM_ERROR(xpvm_tid);
      }
   }
#endif
  
   /* Get the problem data */
   io_u(p);

   /* Start up the graphics window*/
#ifndef WIN32
   init_draw_graph_u(p);
#endif

   p->comp_times.readtime = used_time(&t);

   /* Finds the upper and lower bounds for the problem */
   start_heurs_u(p);

   if (!p->par.do_branch_and_cut){
      printf("\n****************************************************\n");
      printf(  "* Heuristics Finished!!!!!!!                       *\n");
      printf(  "* Now displaying stats and best solution....       *\n");
      printf(  "****************************************************\n\n");
      total_time += p->comp_times.ub_overhead + p->comp_times.ub_heurtime;
      total_time += p->comp_times.lb_overhead + p->comp_times.lb_heurtime;
#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
      printf( "  Problem IO     %.3f\n", p->comp_times.readtime);
      printf( "  Overhead: UB   %.3f\n", p->comp_times.ub_overhead);
      printf( "            LB   %.3f\n", p->comp_times.lb_overhead);
      printf( "  Runtime:  UB   %.3f\n", p->comp_times.ub_heurtime);
      printf( "            LB   %.3f\n", p->comp_times.lb_heurtime);
      printf( "  Total User Time    %.3f\n", total_time);
#endif
      printf( "  Total Real Time    %.3f\n\n", wall_clock(NULL) - start_time);
      printf( "Upper Bound: %.3f\n", p->ub);
      display_solution_u(p, 0);
      free_master_u(p); /* free the problem data structures */
      if (p->par.tm_par.lp_machs)
	 FREE(p->par.tm_par.lp_machs[0]);
      FREE(p->par.tm_par.lp_machs);
      FREE(p);
#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
      pvm_catchout(0);
      comm_exit();
#endif
      exit(1);
   }

   (void) used_time(&t);
   
   /*------------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*------------------------------------------------------------------------*/

   base = (base_desc *) calloc(1, sizeof(base_desc));
   root = (node_desc *) calloc(1, sizeof(node_desc));
   
   initialize_root_node_u(p, base, root);

#ifndef COMPILE_IN_TM
   /*------------------------------------------------------------------------*\
    * Start the tree manager and send the parameters
   \*------------------------------------------------------------------------*/

   if (p->par.tm_machine_set){
      spawn(p->par.tm_exe, (char **)NULL, p->par.tm_debug | TaskHost,
	    p->par.tm_machine, 1, &p->tm_tid);
   }else{
      spawn(p->par.tm_exe, (char **)NULL, p->par.tm_debug, (char *)NULL, 1,
	    &p->tm_tid);
   }
   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&p->par.tm_par), sizeof(tm_params));
   send_char_array(&p->has_ub, 1);
   if (p->has_ub)
      send_dbl_array(&p->ub, 1);
   send_char_array(&p->has_ub_estimate, 1);
   if (p->has_ub_estimate)
      send_dbl_array(&p->ub_estimate, 1);
   if (p->par.tm_par.lp_mach_num)
      send_char_array(p->par.tm_par.lp_machs[0],
		      p->par.tm_par.lp_mach_num*MACH_NAME_LENGTH);
   if (p->par.tm_par.cg_mach_num)
      send_char_array(p->par.tm_par.cg_machs[0],
		      p->par.tm_par.cg_mach_num*MACH_NAME_LENGTH);
   if (p->par.tm_par.cp_mach_num)
      send_char_array(p->par.tm_par.cp_machs[0],
		      p->par.tm_par.cp_mach_num*MACH_NAME_LENGTH);
   send_int_array(&base->varnum, 1);
   send_int_array(&base->cutnum, 1);
#ifdef TRACE_PATH
   {
      int feas_sol;
      int *feas_sol_size;

      if (user_send_feas_sol(p->user, &feas_sol_size, &feas_sol)==USER_NO_PP){
	 send_int_array(&feas_sol_size, 1);
	 if (feas_sol_size){
	    send_int_array(feas_sol, feas_sol_size);
	 }
      }
   }
#endif   
   send_msg(p->tm_tid, TM_DATA);
      
   /*------------------------------------------------------------------------*\
    * Send out the root node
   \*------------------------------------------------------------------------*/

   if (!p->par.warm_start){
      repricing = FALSE;
      node_type = ROOT_NODE;
      
      s_bufid = init_send(DataInPlace);
      send_char_array(&repricing, 1);
      send_char_array(&node_type, 1);
      send_dbl_array(&p->lb, 1);
      send_int_array(&root->nf_status, 1);
      pack_array_desc(&root->uind);
      if (root->nf_status == NF_CHECK_AFTER_LAST ||
	  root->nf_status == NF_CHECK_UNTIL_LAST)
	 pack_array_desc(&root->not_fixed);
      pack_array_desc(&root->cutind);
      pack_basis(&root->basis, TRUE);
      send_int_array(&root->desc_size, 1);
      if (root->desc_size)
	 send_char_array(root->desc, root->desc_size);
      if (root->cutind.size > 0){ /* Hey, we have cuts! Pack them, too. */
	 /* Pack their number again, so we can call unpack_cut_set in TM */
	 int i;
	 send_int_array(&root->cutind.size, 1);
	 for (i = 0; i < root->cutind.size; i++)
	    pack_cut(root->cuts[i]);
      }
      send_msg(p->tm_tid, TM_ROOT_DESCRIPTION);
      freebuf(s_bufid);
      
      FREE(root->desc);
      FREE(root->uind.list);
      FREE(root->not_fixed.list);
      FREE(root->cutind.list);
      FREE(root);
   }
#else
   
   /*------------------------------------------------------------------------*\
    * Send out problem data if needed
   \*------------------------------------------------------------------------*/

#ifdef COMPILE_IN_LP
   send_lp_data_u(p, 0, base);
   lp_data_sent = TRUE;
#ifdef COMPILE_IN_CG
   send_cg_data_u(p, 0);
   cg_data_sent = TRUE;
#endif
#ifdef COMPILE_IN_CP
   send_cp_data_u(p, 0);
   cp_data_sent = TRUE;
#endif
#endif
   
   tm = tm_initialize(base, root, 0);

#ifdef TRACE_PATH
   {
      int feas_sol_size;
      int *feas_sol;
      
      if (user_send_feas_sol(p->user, &feas_sol_size, &feas_sol)==USER_NO_PP){
	 tm->feas_sol_size = feas_sol_size;
	 tm->feas_sol = (int *) calloc (tm->feas_sol_size, sizeof(int));
	 memcpy((char *)tm->feas_sol, (char *)feas_sol, feas_sol_size * ISIZE);
      }
   }
#endif
#endif   
   
   /*------------------------------------------------------------------------*\
    * Wait for messages
   \*------------------------------------------------------------------------*/
   
#ifdef COMPILE_IN_TM
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   while (!lp_data_sent || !cg_data_sent || !cp_data_sent || !sp_data_sent){
   /*___END_EXPERIMENTAL_SECTION___*/
   /*UNCOMMENT FOR PRODUCTION CODE*/
#if 0
   while (!lp_data_sent || !cg_data_sent || !cp_data_sent){
#endif
#else
   do{
#endif
      r_bufid = treceive_msg(ANYONE, ANYTHING, &timeout);
      if (r_bufid == 0){
#ifndef COMPILE_IN_TM
	 if (pstat(p->tm_tid) != PROCESS_OK){
	    printf("\nThe treemanager has died :-(\n\n");
#else
	 if (!processes_alive(p->tm)){
#endif
	    termcode = msgtag = SOMETHING_DIED;
	    break;
	 }else{
	    continue;
	 }
      }
      bufinfo(r_bufid, &bytes, &msgtag, &sender);

      switch (msgtag){
       case FEASIBLE_SOLUTION_NONZEROS:
       case FEASIBLE_SOLUTION_USER:
	 receive_feasible_solution_u(p, msgtag);
	 if (p->par.verbosity > 0){
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
	    display_solution_u(p, p->tm->opt_thread_num);
#else
	    display_solution_u(p, 0);
#endif
	 }
	 break;

       case REQUEST_FOR_LP_DATA:
	 /* An LP process has been started and asks for all necessary data */
	 send_lp_data_u(p, sender, base);
	 lp_data_sent = TRUE;
	 break;

       case REQUEST_FOR_CG_DATA:
	 /* A CG process has been started and asks for all necessary data */
	 send_cg_data_u(p, sender);
	 cg_data_sent = TRUE;
	 break;

       case REQUEST_FOR_CP_DATA:
	 /* A CP process has been started and asks for all necessary data */
	 send_cp_data_u(p, sender);
	 cp_data_sent = TRUE;
	 break;

       /*__BEGIN_EXPERIMENTAL_SECTION__*/
       case REQUEST_FOR_SP_DATA:
	 /* An SP process has been started and asks for all necessary data */
	 send_sp_data_u(p, sender);
	 sp_data_sent = TRUE;
	 break;

       /*___END_EXPERIMENTAL_SECTION___*/
       case TM_FIRST_PHASE_FINISHED:
	 receive_char_array((char *)(&p->comp_times.bc_time),
			     sizeof(node_times));
	 receive_dbl_array(&lb, 1);
	 if (lb > p->lb) p->lb = lb;
	 receive_char_array((char *)&p->stat, sizeof(tm_stat));
	 printf( "\n");
	 printf( "****************************************************\n");
	 printf( "* Branch and Cut First Phase Finished!!!!          *\n");
	 printf( "* Now displaying stats and best solution...        *\n");
	 printf( "****************************************************\n\n");

	 print_statistics(&(p->comp_times.bc_time), &(p->stat), p->ub, p->lb,
			  0, start_time);
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
	 display_solution_u(p, p->tm->opt_thread_num);
#else
	 display_solution_u(p, 0);
#endif
	 break;

       case SOMETHING_DIED:
       case TIME_LIMIT_EXCEEDED:
       case TM_FINISHED:
	 receive_char_array((char *)(&p->comp_times.bc_time),
			    sizeof(node_times));
	 receive_dbl_array(&lb, 1);
	 if (lb > p->lb) p->lb = lb;
	 receive_char_array((char *)&p->stat, sizeof(tm_stat));
	 break;

       default:
	 process_own_messages_u(p, msgtag);
	 break;
      }
      freebuf(r_bufid);

#ifndef COMPILE_IN_TM
   }while (msgtag != TM_FINISHED && msgtag != SOMETHING_DIED &&
	   msgtag != TIME_LIMIT_EXCEEDED);

   termcode = msgtag;
#else
   }
   
   /*------------------------------------------------------------------------*\
    * Solve the problem and receive solutions                         
   \*------------------------------------------------------------------------*/

   tm->start_time += start_time;
   termcode = solve(tm);

   tm_close(tm, termcode);

#ifndef COMPILE_IN_LP
   if (termcode == TM_FINISHED || termcode == TIME_LIMIT_EXCEEDED){
      do{
	 r_bufid = receive_msg(ANYONE, ANYTHING);
	 if (r_bufid == 0){
	    printf("\nError receiving solution ...\n");
	    break;
	 }
	 bufinfo(r_bufid, &bytes, &msgtag, &sender);
	 if (msgtag == FEASIBLE_SOLUTION_NONZEROS ||
	     msgtag == FEASIBLE_SOLUTION_USER){
	    receive_feasible_solution_u(p, msgtag);
	 }
      }while (msgtag != FEASIBLE_SOLUTION_NONZEROS &&
	      msgtag != FEASIBLE_SOLUTION_USER);
   }
#endif
#endif

   /*------------------------------------------------------------------------*\
    * Display the the results and solution data                               
   \*------------------------------------------------------------------------*/

   if (termcode == TM_FINISHED){
      printf("\n****************************************************\n");
      printf(  "* Branch and Cut Finished!!!!!!!                   *\n");
      printf(  "* Now displaying stats and optimal solution...     *\n");
      printf(  "****************************************************\n\n");
   }else if (termcode == TIME_LIMIT_EXCEEDED){
      printf("\n****************************************************\n");
      printf(  "* Time Limit Exceeded :(                           *\n");
      printf(  "* Now displaying stats and best solution...        *\n");
      printf(  "****************************************************\n\n");
   }else{
      printf(
	      "***********Something has died -- halting the machine\n\n");
      printf(
	      "***********Printing out partial data\n\n");
   }

   total_time  = p->comp_times.readtime;
   total_time += p->comp_times.ub_overhead + p->comp_times.ub_heurtime;
   total_time += p->comp_times.lb_overhead + p->comp_times.lb_heurtime;

#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
   printf( "====================== Misc Timing =========================\n");
   printf( "  Problem IO        %.3f\n", p->comp_times.readtime);
   printf( "  UB overhead:      %.3f\n", p->comp_times.ub_overhead);
   printf( "  UB runtime:       %.3f\n", p->comp_times.ub_heurtime);
   printf( "  LB overhead:      %.3f\n", p->comp_times.lb_overhead);
   printf( "  LB runtime:       %.3f\n", p->comp_times.lb_heurtime);
#endif
   
#ifdef COMPILE_IN_TM
   if (tm->lb > p->lb) p->lb = tm->lb;
   print_statistics(&(tm->comp_times), &(tm->stat), tm->ub, p->lb, total_time,
		    start_time);
#ifdef COMPILE_IN_LP
   display_solution_u(p, p->tm->opt_thread_num);
#else
   display_solution_u(p, 0);
#endif
#else
   print_statistics(&(p->comp_times.bc_time), &(p->stat), p->ub, p->lb, 0,
		    start_time);
   display_solution_u(p, 0);
#endif

   if (p->par.do_draw_graph){
      s_bufid = init_send(DataInPlace);
      send_msg(p->dg_tid, CTOI_YOU_CAN_DIE);
      freebuf(s_bufid);
   }

#ifdef COMPILE_IN_TM
   FREE(root);
   FREE(base);
   free_tm(tm);
#else
   FREE(base->lb);
   FREE(base->ub);
   FREE(base->userind);
   FREE(base);
#endif   

   free_master_u(p);

   if (p->par.tm_par.lp_machs)
      FREE(p->par.tm_par.lp_machs[0]);
   FREE(p->par.tm_par.lp_machs);
   if (p->par.tm_par.cg_machs)
      FREE(p->par.tm_par.cg_machs[0]);
   FREE(p->par.tm_par.cg_machs);
   if (p->par.tm_par.cp_machs)
      FREE(p->par.tm_par.cp_machs[0]);
   FREE(p->par.tm_par.cp_machs);
   FREE(p);

#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
   pvm_catchout(0);
   comm_exit();
#endif
   
   return(0);
}

/*===========================================================================*/

void print_statistics(node_times *tim, tm_stat *stat, double ub, double lb,
		      double initial_time, double start_time)
{
   static str_int nfstatus[4] = {
      {"NF_CHECK_ALL"           , NF_CHECK_ALL }
      , {"NF_CHECK_AFTER_LAST"    , NF_CHECK_AFTER_LAST }
      , {"NF_CHECK_UNTIL_LAST"    , NF_CHECK_UNTIL_LAST }
      , {"NF_CHECK_NOTHING"       , NF_CHECK_NOTHING }
   };

   initial_time += tim->communication;
   initial_time += tim->lp;
   initial_time += tim->separation;
   initial_time += tim->fixing;
   initial_time += tim->pricing;
   initial_time += tim->strong_branching;
   initial_time += tim->cut_pool;
#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
   printf("====================== CP Timing =========================\n");
   printf("  Cut Pool                  %.3f\n", tim->cut_pool);
#endif
   printf("====================== LP/CG Timing =========================\n");
#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
   printf("  LP: Solution Time         %.3f\n", tim->lp);
   printf("      Variable Fixing       %.3f\n", tim->fixing);
   printf("      Pricing               %.3f\n", tim->pricing);
   printf("      Strong Branching      %.3f\n", tim->strong_branching);
   printf("      Communication         %.3f\n", tim->communication);
#ifndef COMPILE_IN_LP
   printf("      Ramp Up Time (TM)     %.3f\n", tim->ramp_up_tm);
   printf("      Ramp Up Time (LP)     %.3f\n", tim->ramp_up_lp);
   printf("      Ramp Down Time        %.3f\n", tim->ramp_down_time);
#endif
   printf("      Idle Time (Node Pack) %.3f\n", tim->start_node);
   printf("      Idle Time (Nodes)     %.3f\n", tim->idle_node);
   printf("      Idle Time (Names)     %.3f\n", tim->idle_names);
   printf("      Idle Time (Diving)    %.3f\n", tim->idle_diving);
   printf("      Idle Time (Cuts)      %.3f\n", tim->idle_cuts);
   printf("  Separation                %.3f\n", tim->separation); 
   printf("  Total User Time              %.3f\n", initial_time);
#endif
   printf("  Total Real Time              %.3f\n\n", wall_clock(NULL)-
	  start_time);
   printf("====================== Statistics =========================\n");
   printf("Number of created nodes :       %i\n", stat->created);
   printf("Number of analyzed nodes:       %i\n", stat->analyzed);
   printf("Depth of tree:                  %i\n", stat->max_depth);
   printf("Size of the tree:               %i\n", stat->tree_size);
   printf("Leaves before trimming:         %i\n", stat->leaves_before_trimming);
   printf("Leaves after trimming:          %i\n", stat->leaves_after_trimming);
   printf("Repriced root's nf_status:      %s\n",
	  nfstatus[(int)stat->nf_status].str);
   printf("Not fixed variable num:         %i\n", stat->vars_not_priced);
   printf("Number of Chains:               %i\n", stat->chains);
   printf("Number of Diving Halts:         %i\n", stat->diving_halts);
   printf("Number of cuts in cut pool:     %i\n", stat->cuts_in_pool);
   printf("Lower Bound in Root:            %.3f\n", stat->root_lb);
   if (lb > 0){
      printf("\nCurrent Upper Bound:         %.3f", ub);
      printf("\nCurrent Lower Bound:         %.3f", lb);
      printf("\nGap Percentage:              %.2f\n", 100*(ub-lb)/ub);
   }else{
      printf("\nUpper Bound:        %.3f\n", ub);
   }
}
