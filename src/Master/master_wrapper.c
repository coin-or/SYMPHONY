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
#include <stdlib.h>          /* malloc() is defined here in AIX ... */
#include <stdio.h>

#include "qsortucb.h"
#include "messages.h"
#include "proccomm.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "master.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the wrapper functions for the master process.
\*===========================================================================*/

void initialize_u(problem *p)
{
   CALL_USER_FUNCTION( user_initialize(&p->user) );
}

/*===========================================================================*/

void free_master_u(problem *p)
{
   CALL_USER_FUNCTION( user_free_master(&p->user) );
}

/*===========================================================================*/

void readparams_u(problem *p, int argc, char **argv)
{
   bc_readparams(p, argc, argv);
   CALL_USER_FUNCTION(user_readparams(p->user, p->par.param_file, argc, argv));
}

/*===========================================================================*/

void io_u(problem *p)
{
   CALL_USER_FUNCTION( user_io(p->user) );
}

/*===========================================================================*/

void init_draw_graph_u(problem *p)
{
   if (p->par.do_draw_graph){ /*start up the graphics window*/
      int s_bufid;
      if (p->par.dg_machine_set){
	 spawn(p->par.dg_exe, (char **)NULL, p->par.dg_debug | TaskHost,
	       p->par.dg_machine, 1, &p->dg_tid);
      }else{
	 spawn(p->par.dg_exe, (char **)NULL, p->par.dg_debug, (char *)NULL, 1,
	       &p->dg_tid);
      }
      s_bufid = init_send(DataInPlace);
      send_char_array((char *)&p->par.dg_par, sizeof(dg_params));
      send_msg(p->dg_tid, DG_DATA);
      freebuf(s_bufid);

      if (p->dg_tid)
	 CALL_USER_FUNCTION( user_init_draw_graph(p->user, p->dg_tid) );
   }
}

/*===========================================================================*/

void start_heurs_u(problem *p)
{
   double ub = p->has_ub ? p->ub : -MAXDOUBLE;
   double ub_estimate = p->has_ub_estimate ? p->ub_estimate : -MAXDOUBLE;

   CALL_USER_FUNCTION( user_start_heurs(p->user, &ub, &ub_estimate) );

   if (!p->has_ub){
      if (ub > -MAXDOUBLE){
	 p->has_ub = TRUE;
	 p->ub = ub;
      }
   }else if (ub < p->ub){
      p->ub = ub;
   }
   if (!p->has_ub_estimate){
      if (ub_estimate > -MAXDOUBLE){
	 p->has_ub_estimate = TRUE;
	 p->ub_estimate = ub_estimate;
      }
   }else if (ub_estimate < p->ub_estimate){
      p->ub_estimate = ub_estimate;
   }
   if (p->par.tm_par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f = NULL;
      if (!(f = fopen(p->par.tm_par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 fprintf(f, "00:00:00.00 U %.2f \n", p->ub);
	 fclose(f); 
      }
   }else if (p->par.tm_par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$U %.2f\n", p->ub);
   }
}

/*===========================================================================*/

base_desc *set_base_u(problem *p)
{
   base_desc *base = (base_desc *) calloc(1, sizeof(base_desc));

   switch (user_set_base(p->user, &base->varnum, &base->userind,
			 &base->lb, &base->ub, &base->cutnum,
			 p->par.tm_par.colgen_strat)){
    case ERROR:
      printf("\n\n*********User error detected -- aborting***********\n\n");
      exit(1000);
    case USER_NO_PP:
      if (base->varnum)
	 qsortucb_i(base->userind, base->varnum);
    case USER_AND_PP:
    default:
      break;
   }
      
   return(base);
}

/*===========================================================================*/

node_desc *create_root_u(problem *p)
{
   node_desc *root = (node_desc *) calloc(1, sizeof(node_desc));
   
   /* set some defaults for root first */
   root->uind.type = EXPLICIT_LIST;
   root->cutind.type = EXPLICIT_LIST;
   root->not_fixed.type = EXPLICIT_LIST;
   root->basis.basis_exists = FALSE;
   root->nf_status = NF_CHECK_NOTHING;

   switch (user_create_root(p->user, &root->uind.size, &root->uind.list)){
    case ERROR:
      printf("\n\n*********User error detected -- aborting***********\n\n");
      exit(1000);
    case USER_NO_PP:
      if (root->uind.size)
	 qsortucb_i(root->uind.list, root->uind.size);
    case USER_AND_PP:
    default:
      break;
   }
   root->nf_status = (p->par.tm_par.colgen_strat[0] & COLGEN__FATHOM) ?
                      NF_CHECK_ALL : NF_CHECK_NOTHING;
      
   return(root);
}

/*===========================================================================*/

void receive_feasible_solution_u(problem *p, int msgtag)
{
   receive_int_array(&(p->best_sol.xlevel), 1);
   receive_int_array(&(p->best_sol.xindex), 1);
   receive_int_array(&(p->best_sol.xiter_num), 1);
   receive_dbl_array(&(p->best_sol.lpetol), 1);
   receive_dbl_array(&(p->best_sol.objval), 1);
   receive_int_array(&(p->best_sol.xlength), 1);
   if (p->best_sol.xlength > 0){
      FREE(p->best_sol.xind);
      FREE(p->best_sol.xval);
      p->best_sol.xind = (int *) malloc(p->best_sol.xlength * ISIZE);
      p->best_sol.xval = (double *) malloc(p->best_sol.xlength * DSIZE);
      receive_int_array(p->best_sol.xind, p->best_sol.xlength);
      receive_dbl_array(p->best_sol.xval, p->best_sol.xlength);
   }
   if (!p->has_ub || p->best_sol.objval < p->ub){
      p->has_ub = TRUE;
      p->ub = p->best_sol.objval;
   }
   
   switch (msgtag){
    case FEASIBLE_SOLUTION_NONZEROS:
      break;

    case FEASIBLE_SOLUTION_USER:
      /* A feasible solution has been found in the LP process, and
       * it was packed by the user */
      CALL_USER_FUNCTION( user_receive_feasible_solution(p->user, msgtag,
							 p->best_sol.objval,
							 p->best_sol.xlength,
							 p->best_sol.xind,
							 p->best_sol.xval) );
      break;
   }
}

/*===========================================================================*/

void send_lp_data_u(problem *p, int sender, base_desc *base)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   int i;
   tm_prob *tm = p->tm = get_tm_ptr();

   tm->par.max_active_nodes = p->par.tm_par.max_active_nodes;

#ifdef _OPENMP
   omp_set_dynamic(FALSE);
   omp_set_num_threads(tm->par.max_active_nodes);
#else
   tm->par.max_active_nodes = 1;
#endif

   tm->lpp = (lp_prob **) malloc(tm->par.max_active_nodes * sizeof(lp_prob *));

#pragma omp parallel for
   for (i = 0; i < tm->par.max_active_nodes; i ++){
      tm->lpp[i] = (lp_prob *) calloc(1, sizeof(lp_prob));
      tm->lpp[i]->proc_index = i;
      tm->lpp[i]->par = p->par.lp_par;

      if ((tm->lpp[i]->has_ub = p->has_ub))
	 tm->lpp[i]->ub = p->ub;
      else
	 p->ub = - (MAXDOUBLE / 2);
      
      tm->lpp[i]->draw_graph = p->dg_tid;
      tm->lpp[i]->base = *base;
      
      CALL_USER_FUNCTION( user_send_lp_data(p->user, &(tm->lpp[i]->user)) );
   }
   get_lp_ptr(tm->lpp);
#else   
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&p->par.lp_par), sizeof(lp_params));
   send_char_array(&p->has_ub, 1);
   if (p->has_ub)
      send_dbl_array(&p->ub, 1);
   send_int_array(&p->dg_tid, 1);
   send_int_array(&base->varnum, 1);
   if (base->varnum){
      send_int_array(base->userind, base->varnum);
      send_dbl_array(base->lb, base->varnum);
      send_dbl_array(base->ub, base->varnum);
   }
   send_int_array(&base->cutnum, 1);
   CALL_USER_FUNCTION( user_send_lp_data(p->user, NULL) );
   send_msg(sender, LP_DATA);
   freebuf(s_bufid);
#endif
}

/*===========================================================================*/

void send_cg_data_u(problem *p, int sender)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined(COMPILE_IN_CG)
   int i;
   tm_prob *tm = p->tm;
   cg_prob **cg_list = (cg_prob **)
                       malloc(tm->par.max_active_nodes*sizeof(cg_prob *));
#pragma omp parallel for
   for (i = 0; i < tm->par.max_active_nodes; i++){
      tm->lpp[i]->cgp = cg_list[i] = (cg_prob *) calloc(1, sizeof(cg_prob));
      
      cg_list[i]->par = p->par.cg_par;
      
      cg_list[i]->draw_graph = p->dg_tid;
      
      CALL_USER_FUNCTION( user_send_cg_data(p->user,
					    &(tm->lpp[i]->cgp->user)) );
   }
   get_cg_ptr(cg_list);
#else
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&p->par.cg_par), sizeof(cg_params));
   send_int_array(&p->dg_tid, 1);
   CALL_USER_FUNCTION( user_send_cg_data(p->user, NULL) );
   send_msg(sender, CG_DATA);
   freebuf(s_bufid);
#endif
}

/*===========================================================================*/

void send_cp_data_u(problem *p, int sender)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP)
   int i;
   tm_prob *tm = p->tm;

   tm->cpp = (cut_pool **) malloc(p->par.tm_par.max_cp_num*sizeof(cut_pool *));
   for (i = 0; i < p->par.tm_par.max_cp_num; i++){
      tm->cpp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
      tm->cpp[i]->par = p->par.cp_par;
      CALL_USER_FUNCTION( user_send_cp_data(p->user, &p->tm->cpp[i]->user) );
   }
   if (p->par.tm_par.max_cp_num)
      get_cp_ptr(tm->cpp, 0);
#else
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&p->par.cp_par), sizeof(cp_params));
   CALL_USER_FUNCTION( user_send_cp_data(p->user, NULL) );
   send_msg(sender, CP_DATA);
   freebuf(s_bufid);
#endif
}

/*__BEGIN_EXPERIMENTAL_SECTION__*/
/*===========================================================================*/

void send_sp_data_u(problem *p, int sender)
{
#ifdef COMPILE_DECOMP
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&p->par.sp_par), sizeof(sp_params));
   CALL_USER_FUNCTION( user_send_sp_data(p->user) );
   send_msg(sender, SP_DATA);
   freebuf(s_bufid);
#endif
}

/*___END_EXPERIMENTAL_SECTION___*/
/*===========================================================================*/

void display_solution_u(problem *p, int thread_num)
{
   int user_res, i;
   lp_sol sol;

   sol.xlength = 0;
   
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   if (p->tm && p->tm->lpp[thread_num]){
      sol = p->tm->lpp[thread_num]->best_sol;
   }
#else
   sol = p->best_sol;
#endif
   
   if (!sol.xlength){
      printf("\nNo Solution Found\n\n");
      return;
   }

   printf("\nSolution Found: Node %i, Level %i\n", sol.xindex, sol.xlevel);
   printf("Solution Cost: %.3f\n", p->tm->ub);
   qsortucb_id(sol.xind, sol.xval, sol.xlength);
   
   user_res = user_display_solution(p->user, sol.lpetol, sol.xlength, sol.xind,
				    sol.xval, sol.objval);
   
   switch(user_res){
    case USER_NO_PP:
      return;
    case USER_AND_PP:
    case DEFAULT:
      if (sol.xlength){
	 printf("\nUser indices and values of nonzeros in the solution:\n");
	 printf("\nINDEX     VALUE\n");
	 printf("=====     =====\n");
	 for (i = 0; i < sol.xlength; i++)
	    printf("%5d %10.3f\n", sol.xind[i], sol.xval[i]);
	 return;
      }
   }
}

/*===========================================================================*/

void process_own_messages_u(problem *p, int msgtag)
{
   CALL_USER_FUNCTION( user_process_own_messages(p->user, msgtag) );
}

