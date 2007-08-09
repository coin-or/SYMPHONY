/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_TM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sym_tm.h"
#include "sym_constants.h"
#include "sym_types.h"
#include "sym_macros.h"
#include "sym_messages.h"
#include "sym_proccomm.h"
#include "sym_timemeas.h"
#include "sym_pack_cut.h"
#include "sym_pack_array.h"
#ifdef COMPILE_IN_LP
#include "sym_lp.h"
#endif
#ifdef COMPILE_IN_TM
#include "sym_master.h"
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains function related primarily to process communication
 * to/from the treemanager.
\*===========================================================================*/

/*===========================================================================*\
 * This function starts up a number of the specified processes
\*===========================================================================*/

process_set start_processes(tm_prob *tm,
			    int procnum, char *procname, int procdebug,
			    int machnum, char **mach)
{
   int s_bufid, i;
   process_set pset;

   pset.procnum = procnum;
   pset.procs = (int *) malloc(procnum * ISIZE);
   pset.free_num = procnum;
   pset.free_ind = (int *) malloc(procnum * ISIZE);
   for (i = procnum - 1; i >= 0; i--)
      pset.free_ind[i] = i;

   if (machnum){
      for (i = 0; i < procnum; i++){
	 spawn(procname, (char **)NULL, procdebug + TaskHost,
	       mach[i % machnum], 1, pset.procs + i);
      }
   }else{
      spawn(procname, (char **)NULL, procdebug, (char *)NULL,
	    procnum, pset.procs);
   }
   /*------------------------------------------------------------------------*\
    * Send the master tid info.
   \*------------------------------------------------------------------------*/
   s_bufid = init_send(DataInPlace);
   send_int_array(&tm->master, 1);
   send_int_array(&i, 1);
   msend_msg(pset.procs, procnum, MASTER_TID_INFO);

   return(pset);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function stops a process set.
\*===========================================================================*/

void stop_processes(process_set *pset)
{
   if (pset->procnum){
      int s_bufid;
      s_bufid = init_send(DataInPlace);
      msend_msg(pset->procs, pset->procnum, YOU_CAN_DIE);
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Checks to make sure that all the processes are alive.
\*===========================================================================*/

char processes_alive(tm_prob *tm)
{
   int i;

#ifndef COMPILE_IN_LP
   for (i = tm->lp.procnum-1; i>=0; i--){
      if (pstat(tm->lp.procs[i]) != PROCESS_OK){
	 printf("\nLP process has died -- halting machine\n\n");
	 return(FALSE);
      }
   }
#endif
   for (i = tm->cg.procnum-1; i>=0; i--){
      if (pstat(tm->cg.procs[i]) != PROCESS_OK){
	 printf("\nCG process has died -- halting machine\n\n");
	 return(FALSE);
      }
   }
#ifndef COMPILE_IN_CP
   for (i = tm->cp.procnum-1; i>=0; i--){
      if (pstat(tm->cp.procs[i]) != PROCESS_OK){
	 printf("\nCP process has died -- halting machine\n\n");
	 return(FALSE);
      }
   }
#endif
   return(TRUE);
}

/*===========================================================================*/

/*===========================================================================*\
 * Send the active node to the LP
\*===========================================================================*/

void send_active_node(tm_prob *tm, bc_node *node, char colgen_strat,
		      int thread_num)
{
#ifdef COMPILE_IN_LP
   lp_prob **lp = tm->lpp;
   node_desc *new_desc;
#else
   int length;
   char dive;
   int s_bufid;
#endif
   int level, i, j;
   bc_node **path, *n;
   branch_desc *branch_path, *bpath;
   branch_obj *bobj;
   node_desc *desc = &node->desc;
   int varexp_ind = 0, cutexp_ind = 0, nfexp_ind = 0;
   int bv_ind = 0, br_ind = 0, ev_ind = 0, er_ind = 0;
   array_desc extravar = { EXPLICIT_LIST, 0, 0, NULL };
   array_desc extrarow = { EXPLICIT_LIST, 0, 0, NULL };
   array_desc not_fixed = { EXPLICIT_LIST, 0, 0, NULL };
   basis_desc basis;

   int *list, *stat;

   char deal_with_nf = (desc->nf_status == NF_CHECK_AFTER_LAST ||
			desc->nf_status == NF_CHECK_UNTIL_LAST);

   if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f;
#pragma omp critical(write_vbc_emulation_file)
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 PRINT_TIME(tm, f);
	 fprintf(f, "P %i %i\n",node->bc_index+1,VBC_ACTIVE_NODE);
	 fclose(f); 
      }
   }else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$P %i %i\n", node->bc_index+1, VBC_ACTIVE_NODE);
   }

   memset((char *)(&basis), 0, sizeof(basis_desc));

   /*------------------------------------------------------------------------*\
    * First go up in the search tree to the root and record for every field
    * the node where the first explicit description occurs.
   \*------------------------------------------------------------------------*/
   level = node->bc_level;
   REMALLOC(tm->rpath[thread_num], bc_node *, tm->rpath_size[thread_num],
	    2*(level+1), BB_BUNCH);
   path = tm->rpath[thread_num];
   REMALLOC(tm->bpath[thread_num], branch_desc, tm->bpath_size[thread_num],
	    2*(level+1), BB_BUNCH);
   branch_path = tm->bpath[thread_num];

   if (desc->uind.type == NO_DATA_STORED){
      varexp_ind = -1;
   }else{
      for (i = level, n = node; !varexp_ind && i > 0; n = n->parent, i--)
	 if (n->desc.uind.type == EXPLICIT_LIST)
	    varexp_ind = i;
   }
   if (desc->cutind.type == NO_DATA_STORED){
      cutexp_ind = -1;
   }else{
      for (i = level, n = node; !cutexp_ind && i > 0; n = n->parent, i--)
	 if (n->desc.cutind.type == EXPLICIT_LIST)
	    cutexp_ind = i;
   }
   if (deal_with_nf){
      for (i = level, n = node; !nfexp_ind && i > 0; n = n->parent, i--)
	 if (n->desc.not_fixed.type == EXPLICIT_LIST)
	    nfexp_ind = i;
   }
   if ((basis.basis_exists = desc->basis.basis_exists) == TRUE){
      for (i = level, n = node; !bv_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.basevars.type == EXPLICIT_LIST)
	    bv_ind = i;
      for (i = level, n = node; !br_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.baserows.type == EXPLICIT_LIST)
	    br_ind = i;
      for (i = level, n = node; !ev_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.extravars.type == EXPLICIT_LIST)
	    ev_ind = i;
      for (i = level, n = node; !er_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.extrarows.type == EXPLICIT_LIST)
	    er_ind = i;
   }else{
      ev_ind = er_ind = level;
   }

   for (i = level, n = node; i >= 0; n = n->parent, i--)
      path[i] = n;

   /* An upper estimate on the total length of arrays */
   if (varexp_ind >= 0){
      extravar.size = (n = path[varexp_ind]) -> desc.uind.size;
      for (i = varexp_ind + 1; i <= level; i++)
	 extravar.size += path[i]->desc.uind.added;
   }
   if (cutexp_ind >= 0){
      extrarow.size = (n = path[cutexp_ind]) -> desc.cutind.size;
      for (i = cutexp_ind + 1; i <= level; i++)
	 extrarow.size += path[i]->desc.cutind.added;
   }
   if (deal_with_nf && nfexp_ind >= 0){
      not_fixed.size = (n = path[nfexp_ind]) -> desc.not_fixed.size;
      for (i = nfexp_ind + 1; i <= level; i++)
	 not_fixed.size += path[i]->desc.not_fixed.added;
   }else{
      not_fixed.size = 0;
   }
#ifdef COMPILE_IN_LP
   /* If the LP function is compiled into the tree manager as a single
      executable, then we allocate new memory for these arrays since
      these arrays will be used directly instead of being passed
      through PVM. */
   if (extravar.size){
      extravar.list = (int *) malloc(extravar.size*ISIZE);
      if (basis.basis_exists)
	 basis.extravars.stat = (int *) malloc(extravar.size*ISIZE);
   }
   if (extrarow.size){
      extrarow.list = (int *) malloc(extrarow.size*ISIZE);
      if (basis.basis_exists)
	 basis.extrarows.stat = (int *) malloc(extrarow.size*ISIZE);
   }
   if (not_fixed.size)
      not_fixed.list = (int *) malloc(not_fixed.size*ISIZE);
   if (tm->bvarnum && basis.basis_exists)
      basis.basevars.stat = (int *) malloc(tm->bvarnum*ISIZE);
   if (tm->bcutnum && basis.basis_exists)
      basis.baserows.stat = (int *) malloc(tm->bcutnum*ISIZE);
#else
   /* Here, we just use temporary arrays */
   length = 2 * (extravar.size + extrarow.size + not_fixed.size) +
      tm->bvarnum + tm->bcutnum;
   REMALLOC(tm->tmp.i, int, tm->tmp.i_size, length, BB_BUNCH);
   extravar.list = tm->tmp.i;
   extrarow.list = extravar.list + extravar.size;
   not_fixed.list = extrarow.list + extrarow.size;
   basis.basevars.stat = not_fixed.list + not_fixed.size;
   basis.extravars.stat = basis.basevars.stat + tm->bvarnum;
   basis.baserows.stat = basis.extravars.stat + extravar.size;
   basis.extrarows.stat = basis.baserows.stat + tm->bcutnum;
#endif
   
   /* The extra variables (uind) and the corresponding basis part */
   if (varexp_ind >= 0){
      extravar.size = (n = path[varexp_ind]) -> desc.uind.size;
      if (extravar.size > 0)
	 memcpy(extravar.list, n->desc.uind.list, ISIZE * extravar.size);
      for (i = varexp_ind + 1; i <= ev_ind; i++)
	 modify_list(&extravar, &path[i]->desc.uind);
      if (basis.basis_exists){
	 /* at this point i == ev_ind */
#ifdef DO_TESTS
	 if (extravar.size != path[ev_ind]->desc.basis.extravars.size){
	    printf("problem with extravar.size!!! \n\n");
	    exit(-6);
	 }
#endif
	 if (path[ev_ind]->desc.basis.extravars.size > 0)
	    memcpy(basis.extravars.stat,
		   path[ev_ind]->desc.basis.extravars.stat,
		   path[ev_ind]->desc.basis.extravars.size * ISIZE);
	 for (i = ev_ind + 1; i <= level; i++){
	    modify_list_and_stat(&extravar, basis.extravars.stat,
				 &path[i]->desc.uind,
				 &path[i]->desc.basis.extravars);
	 }
	 /* Although we send an explicit list, the type is sent over to show
	    to the LP process how extravars are stored in TM */
	 basis.extravars.type = node->desc.basis.extravars.type;
	 basis.extravars.size = extravar.size;
	 basis.extravars.list = NULL;
	 /* Now extravar.list/extravar.size and basis.extravars are OK */
	 /* Fix basis.basevars */
	 basis.basevars.type = EXPLICIT_LIST;
	 basis.basevars.size = path[bv_ind]->desc.basis.basevars.size;
	 basis.basevars.list = NULL;
	 if (basis.basevars.size > 0){
	    memcpy(basis.basevars.stat,
		   path[bv_ind]->desc.basis.basevars.stat,
		   basis.basevars.size * ISIZE);
	    for (i = bv_ind + 1; i <= level; i++){
	       list = path[i]->desc.basis.basevars.list;
	       stat = path[i]->desc.basis.basevars.stat;
	       for (j = path[i]->desc.basis.basevars.size - 1; j >= 0; j--)
		  basis.basevars.stat[list[j]] = stat[j];
	    }
	 }
      }
   }


   /* Now take care of cutind and the corresponding basis part */
   if (cutexp_ind >= 0){
      extrarow.size = (n = path[cutexp_ind]) -> desc.cutind.size;
      if (extrarow.size > 0)
	 memcpy(extrarow.list, n->desc.cutind.list, ISIZE * extrarow.size);
      for (i = cutexp_ind + 1; i <= er_ind; i++)
	 modify_list(&extrarow, &path[i]->desc.cutind);
      if (basis.basis_exists){
	 /* at this point i == er_ind */
#ifdef DO_TESTS
	 if (extrarow.size != path[er_ind]->desc.basis.extrarows.size){
	    printf("problem with extrarow.size!!! \n\n");
	    exit(-6);
	 }
#endif
	 if (path[er_ind]->desc.basis.extrarows.size > 0)
	    memcpy(basis.extrarows.stat,
		   path[er_ind]->desc.basis.extrarows.stat,
		   path[er_ind]->desc.basis.extrarows.size * ISIZE);
	 for (i = er_ind + 1; i <= level; i++){
	    modify_list_and_stat(&extrarow, basis.extrarows.stat,
				 &path[i]->desc.cutind,
				 &path[i]->desc.basis.extrarows);
	 }
	 /* Same trick as above */
	 basis.extrarows.type = node->desc.basis.extrarows.type;
	 basis.extrarows.size = extrarow.size;
	 basis.extrarows.list = NULL;
	 /* Now extrarow.list/extrarow.size and basis.extrarows are OK */
	 /* Fix basis.baserows */
	 basis.baserows.type = EXPLICIT_LIST;
	 basis.baserows.size = path[br_ind]->desc.basis.baserows.size;
	 basis.baserows.list = NULL;
	 if (basis.baserows.size > 0){
	    memcpy(basis.baserows.stat,
		   path[br_ind]->desc.basis.baserows.stat,
		   basis.baserows.size * ISIZE);
	    for (i = br_ind + 1; i <= level; i++){
	       list = path[i]->desc.basis.baserows.list;
	       stat = path[i]->desc.basis.baserows.stat;
	       for (j = path[i]->desc.basis.baserows.size - 1; j >= 0; j--)
		  basis.baserows.stat[list[j]] = stat[j];
	    }
	 }
      }
   }

   /* Finally the not fixed ones */
   if (deal_with_nf){
      not_fixed.size = (n = path[nfexp_ind]) -> desc.not_fixed.size;
      if (not_fixed.size > 0)
	 memcpy(not_fixed.list, n->desc.not_fixed.list, ISIZE*not_fixed.size);
      for (i = nfexp_ind + 1; i <= level; i++)
	 modify_list(&not_fixed, &path[i]->desc.not_fixed);
   }

   for (bpath = branch_path, i = 0; i < level; i++, bpath++){
      for (j = path[i]->bobj.child_num - 1; j >= 0; j--)
	 if (path[i]->children[j] == path[i+1])
	    break;
      bobj = &path[i]->bobj;
      bpath->type = bobj->type;
      bpath->name = bobj->name;
      bpath->sense = bobj->sense[j];
      bpath->rhs = bobj->rhs[j];
      bpath->range = bobj->range[j];
      bpath->branch = bobj->branch[j];
   }
   
#ifdef COMPILE_IN_LP

#if 0
   if(!tm->par.sensitivity_analysis){   
      /* Again, here, we need to do some things directly if the LP
	 function is being performed within the tree manager. Otherwise,
	 we just send out the data below */
      if (! (colgen_strat & COLGEN_REPRICING) &&
	  lp[thread_num]->has_ub && node->lower_bound >
	  lp[thread_num]->ub - lp[thread_num]->par.granularity){
	 if (desc->nf_status == NF_CHECK_NOTHING ||
	     (colgen_strat & FATHOM__DO_NOT_GENERATE_COLS__DISCARD)){
	    tm->active_nodes[thread_num]->node_status = NODE_STATUS__PRUNED;
	    if (lp[thread_num]->par.verbosity > 0){
	       printf("***************************************************\n");
	       printf("* Immediately pruning NODE %i LEVEL %i\n",
		   node->bc_index, node->bc_level);
	       printf("***************************************************\n");
	    }
	    return;
	 }
	 if (colgen_strat & FATHOM__DO_NOT_GENERATE_COLS__SEND){
	    tm->active_nodes[thread_num]->node_status = NODE_STATUS__HELD;
	    REALLOC(tm->nextphase_cand, bc_node *,
		    tm->nextphase_cand_size, tm->nextphase_candnum+1,
		    BB_BUNCH);
	    tm->nextphase_cand[tm->nextphase_candnum++] =
	       tm->active_nodes[thread_num];
	    if (lp[thread_num]->par.verbosity > 0){
	       printf("***************************************************\n");
	       printf("* Sending back NODE %i LEVEL %i\n",
		      node->bc_index, node->bc_level);
	       printf("***************************************************\n");
	    }
	    return;
	 }
      }
   }
#endif
   new_desc = lp[thread_num]->desc = (node_desc *) calloc(1,sizeof(node_desc));

   lp[thread_num]->cut_pool = node->cp;
   lp[thread_num]->bc_index = node->bc_index;
   lp[thread_num]->bc_level = node->bc_level;
   lp[thread_num]->lp_data->objval = node->lower_bound;
   lp[thread_num]->colgen_strategy = colgen_strat;
   
   new_desc->nf_status = desc->nf_status;
   new_desc->basis = basis;
   if (deal_with_nf)
      new_desc->not_fixed = not_fixed;
   new_desc->uind = extravar;
   new_desc->cutind = extrarow;

   /* The cuts themselves */
   if (extrarow.size > 0){
      new_desc->cuts = (cut_data **)
	 malloc(extrarow.size*sizeof(cut_data *));
      for (i = 0; i < extrarow.size; i++){
	 new_desc->cuts[i] = tm->cuts[extrarow.list[i]];
      }
   }

   if (level > 0)
      lp[thread_num]->bdesc = branch_path;

   /* diving instructions */
   //lp[thread_num]->dive = shall_we_dive(tm, node->lower_bound);

   /* User defined description */
   new_desc->desc_size = desc->desc_size;
   if (new_desc->desc_size > 0)
      memcpy((char *)new_desc->desc, (char *)desc->desc, new_desc->desc_size);
   
#else
   
   /*------------------------------------------------------------------------*\
    * Now put together the message
   \*------------------------------------------------------------------------*/

   s_bufid = init_send(DataInPlace);
   send_int_array(&node->cp, 1);
   send_int_array(&node->bc_index, 1);
   send_int_array(&node->bc_level, 1);
   send_dbl_array(&node->lower_bound, 1);
   send_char_array(&colgen_strat, 1);
   send_int_array(&desc->nf_status, 1);

   pack_basis(&basis, TRUE);
   if (deal_with_nf)
      pack_array_desc(&not_fixed);
   pack_array_desc(&extravar);
   pack_array_desc(&extrarow);

   /* The cuts themselves */
   for (i = 0; i < extrarow.size; i++)
      pack_cut(tm->cuts[extrarow.list[i]]);

   if (level > 0)
      send_char_array((char *)branch_path, sizeof(branch_desc)*level);

   /* diving instructions */
   dive = shall_we_dive(tm, node->lower_bound);
   send_char_array(&dive, 1);

   /* User defined description */
   send_int_array(&node->desc.desc_size, 1);
   if (node->desc.desc_size){
      send_char_array(node->desc.desc, node->desc.desc_size);
   }

   /* send it out... */
   send_msg(node->lp, LP__ACTIVE_NODE_DATA);
   freebuf(s_bufid);
#endif
}

/*===========================================================================*/

/*===========================================================================*\
 * Receive the description of the node back from the LP process.
\*===========================================================================*/
   
void receive_node_desc(tm_prob *tm, bc_node *n)
{
   char node_type, repricing;
   node_desc *desc = &n->desc;
   node_desc *newdesc;
#ifdef DO_TESTS
   double old_lower_bound  = n->lower_bound;
#endif

#ifdef SENSITIVITY_ANALYSIS
   if (tm->par.sensitivity_analysis){
      if (n->sol){
	 FREE(n->sol);
	 FREE(n->duals);
      }
      receive_int_array(&n->sol_size, 1);
      n->sol = (double *) malloc (DSIZE * n->sol_size);
      receive_dbl_array(n->sol, tm->rootnode->desc.uind.size);
      n->duals = (double *) malloc (DSIZE * tm->bcutnum);
      send_dbl_array(n->duals, tm->bcutnum);
   }
#endif
   receive_char_array(&repricing, 1);
   receive_char_array(&node_type, 1);
   if (node_type == INFEASIBLE_PRUNED || node_type == OVER_UB_PRUNED ||
       node_type == DISCARDED_NODE || node_type == FEASIBLE_PRUNED){
      n->node_status = NODE_STATUS__PRUNED;
      if (node_type == FEASIBLE_PRUNED) {
	 if (!tm->par.sensitivity_analysis){ 
	    receive_int_array(&n->sol_size, 1);
	    n->sol = (double *) malloc (DSIZE * n->sol_size);
	    receive_dbl_array(n->sol, n->sol_size);
	 }
      }
#ifdef TRACE_PATH
      if (n->optimal_path){
	 printf("\n\nAttempting to prune the optimal path!!!!!!!!!\n\n");
	 sleep(600);
	 if (tm->par.logging){
	    write_tm_info(tm, tm->par.tree_log_file_name, NULL, FALSE);
	    write_subtree(tm->rootnode, tm->par.tree_log_file_name, NULL,
			  TRUE, tm->par.logging);
	    write_tm_cut_list(tm, tm->par.cut_log_file_name, FALSE);
	 }
	 exit(1);
      }
#endif
      if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL)
	 write_pruned_nodes(tm, n);
      if (tm->par.keep_description_of_pruned == DISCARD ||
	  tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
	 purge_pruned_nodes(tm, n, node_type == FEASIBLE_PRUNED ?
			    VBC_FEAS_SOL_FOUND : VBC_PRUNED);
      }
      return;
   }

   /* This function is called either when a a node was finished and we really
      do have to unpack the differences OR when we have repriced the root.
      In the later case the LP sends explicit description (and this function
      is called with an empty 'n') so we can still use this function. */
   receive_dbl_array(&n->lower_bound, 1);
#ifdef DO_TESTS
   if (n->lower_bound < old_lower_bound - 10){
      printf("#####Error: lower bound descrease in node from %.3f to %.3f\n",
	     old_lower_bound, n->lower_bound);
   }
#endif

   if (node_type == INTERRUPTED_NODE){
      n->node_status = NODE_STATUS__INTERRUPTED;
      insert_new_node(tm, n);
      return;
   }
   
   newdesc = (node_desc *) calloc(1, sizeof(node_desc));
   /* Unpack the new description */
   receive_int_array(&newdesc->nf_status, 1);
   unpack_array_desc(&newdesc->uind);
   if (newdesc->nf_status == NF_CHECK_AFTER_LAST ||
       newdesc->nf_status == NF_CHECK_UNTIL_LAST)
      unpack_array_desc(&newdesc->not_fixed);
   unpack_array_desc(&newdesc->cutind);
   unpack_basis(&newdesc->basis, FALSE);
   receive_int_array(&desc->desc_size, 1);
   FREE(desc->desc);
   if (desc->desc_size){
      desc->desc = (char *) malloc(desc->desc_size);
      receive_char_array(desc->desc, desc->desc_size);
   }

   /* Now merge the old and the new together... */
   merge_descriptions(desc, newdesc);

   FREE(newdesc);
   
   if (tm->par.verbosity > 10){
      printf("TM: node %4i: ", n->bc_index);
      if (desc->uind.type == WRT_PARENT){
	 printf("uind:WRT(%i,%i) ", desc->uind.size, desc->uind.added);
      }else{
	 printf("uind:EXP(%i) ", desc->uind.size);
      }
      printf("nf:%s ",
	     ((desc->nf_status == NF_CHECK_AFTER_LAST ||
	       desc->nf_status == NF_CHECK_UNTIL_LAST) ?
	      (desc->not_fixed.type == EXPLICIT_LIST ? "EXP":"WRT") : "N/A"));
      if (desc->cutind.type == WRT_PARENT){
	 printf("cind:WRT(%i,%i)\n", desc->cutind.size, desc->cutind.added);
      }else{
	 printf("cind:EXP(%i)\n", desc->cutind.size);
      }
      printf("               bvar:%s evar:%s brow:%s erow:%s\n",
	     desc->basis.basevars.type == EXPLICIT_LIST ? "EXP" : "WRT",
	     desc->basis.extravars.type == EXPLICIT_LIST ? "EXP" : "WRT",
	     desc->basis.baserows.type == EXPLICIT_LIST ? "EXP" : "WRT",
	     desc->basis.extrarows.type == EXPLICIT_LIST ? "EXP" : "WRT");
   }
   if (! repricing){
      /* If it's not repricing then we have to insert the node into the
       * appropriate heap */
      switch (node_type){
       case INFEASIBLE_HOLD_FOR_NEXT_PHASE:
       case OVER_UB_HOLD_FOR_NEXT_PHASE:
	 n->node_status = NODE_STATUS__HELD;
	 REALLOC(tm->nextphase_cand, bc_node *,
		 tm->nextphase_cand_size, tm->nextphase_candnum+1, BB_BUNCH);
	 tm->nextphase_cand[tm->nextphase_candnum++] = n;
	 /* update the nodes_per_... stuff */
	 /* the active_nodes_per_... will be updated when the LP__IS_FREE
	    message comes */
	 if (n->cp)
#ifdef COMPILE_IN_CP
	    tm->nodes_per_cp[n->cp]++;
#else
	    tm->nodes_per_cp[find_process_index(&tm->cp, n->cp)]++;
#endif
	 break;
       case NODE_BRANCHED_ON:
	 n->node_status = NODE_STATUS__BRANCHED_ON;
	 if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
	    FILE *f;
#pragma omp critical(write_vbc_emulation_file)
	    if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	       printf("\nError opening vbc emulation file\n\n");
	    }else{
	       PRINT_TIME(tm, f);
	       fprintf(f, "P %i %i\n", n->bc_index + 1,
		       VBC_INTERIOR_NODE);
	       fclose(f); 
	    }
	 }else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
	    printf("$P %i %i\n", n->bc_index + 1, VBC_INTERIOR_NODE);
	 }
	 break;
       case ROOT_NODE:
	 tm->rootnode = n;
	 n->bc_index = tm->stat.created++;
	 tm->stat.tree_size++;
	 /* these are set by calloc:
	    n->bc_level = 0;
	    n->lp = n->cg = n->cp = n->sp = 0;
	    n->parent = NULL;
	    */
	 n->node_status = NODE_STATUS__ROOT;
	 insert_new_node(tm, n);
	 break;
      }
   }
   if (n->node_status == NODE_STATUS__PRUNED){
#ifdef TRACE_PATH
      if (n->optimal_path){
	 printf("\n\nAttempting to prune the optimal path!!!!!!!!!\n\n");
	 sleep(600);
	 if (tm->par.logging){
	    write_tm_info(tm, tm->par.tree_log_file_name, NULL, FALSE);
	    write_subtree(tm->rootnode, tm->par.tree_log_file_name, NULL,
			  TRUE, tm->par.logging);
	    write_tm_cut_list(tm, tm->par.cut_log_file_name, FALSE);
	 }
	 exit(1);
      }
#endif
      if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_FULL ||
	  tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
	 write_pruned_nodes(tm, n);
	 purge_pruned_nodes(tm, n, node_type == FEASIBLE_PRUNED ?
			    VBC_FEAS_SOL_FOUND : VBC_PRUNED);
      }
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Receives the description of the branching object and processes this
 * information appropriately.
\*===========================================================================*/

void process_branching_info(tm_prob *tm, bc_node *node)
{
   int s_bufid;
   int old_cut_name = 0;
   branch_obj *bobj = &node->bobj;
   char *action;
   int *feasible;
   double *objval;
   int oldkeep, keep;
   char olddive, dive;
   int new_branching_cut = FALSE, lp, i;

   receive_char_array(&bobj->type, 1);
   receive_int_array(&bobj->name, 1);
   if (bobj->type == CANDIDATE_CUT_IN_MATRIX){
      receive_int_array(&new_branching_cut, 1);
      if ( (old_cut_name = bobj->name) == -tm->bcutnum-1){
	 bobj->name = add_cut_to_list(tm, unpack_cut(NULL));
      }
   }
   receive_int_array(&bobj->child_num, 1);
   REMALLOC(tm->tmp.c, char, tm->tmp.c_size, bobj->child_num, BB_BUNCH);
   REMALLOC(tm->tmp.i, int, tm->tmp.i_size, bobj->child_num, BB_BUNCH);
   REMALLOC(tm->tmp.d, double, tm->tmp.d_size, bobj->child_num, BB_BUNCH);
   action = tm->tmp.c;
   feasible = tm->tmp.i;
   objval = tm->tmp.d;

#ifndef MAX_CHILDREN_NUM
   bobj->sense = malloc(bobj->child_num * CSIZE);
   bobj->rhs = (double *) malloc(bobj->child_num * DSIZE);
   bobj->range = (double *) malloc(bobj->child_num * DSIZE);
   bobj->branch = (int *) malloc(bobj->child_num * ISIZE);
#endif
   receive_char_array(bobj->sense, bobj->child_num);
   receive_dbl_array(bobj->rhs, bobj->child_num);
   receive_dbl_array(bobj->range, bobj->child_num);
   receive_int_array(bobj->branch, bobj->child_num);

   receive_dbl_array(objval, bobj->child_num);
   receive_int_array(feasible, bobj->child_num);
   bobj->solutions = (double **) calloc(bobj->child_num, sizeof(double *));
   for (i = 0; i < bobj->child_num; i++){
      if (feasible[i]){
#if 0
	 bobj->solutions[i] = (double *)
	    malloc(DSIZE*tm->rootnode->desc.uind.size); 
	 receive_dbl_array(bobj->solutions[i], tm->rootnode->desc.uind.size);
#endif
      }
   }
   receive_char_array(action, bobj->child_num);

   receive_char_array(&olddive, 1);
   receive_int_array(&keep, 1);
   oldkeep = keep;
   lp = node->lp;
   
   dive = generate_children(tm, node, bobj, objval, feasible, action, olddive,
			    &keep, new_branching_cut); 

   if (oldkeep >= 0 && (olddive == CHECK_BEFORE_DIVE || olddive == DO_DIVE)){
      /* We have to reply */
      s_bufid = init_send(DataInPlace);
      send_char_array(&dive, 1);
      if (dive == DO_DIVE || dive == CHECK_BEFORE_DIVE){
	 /* Give the index of the node kept and also the index of the
	  * branching cut if necessary */
	 send_int_array(&node->children[keep]->bc_index, 1);
	 if (bobj->type == CANDIDATE_CUT_IN_MATRIX &&
	     old_cut_name == -tm->bcutnum-1)
	    send_int_array(&bobj->name, 1);
	 node->children[keep]->lp = node->lp;
	 node->children[keep]->cg = node->cg;
	 /* update the info which node is processed by that lp process */
	 tm->active_nodes[find_process_index(&tm->lp, node->lp)] =
	    node->children[keep];
	 tm->stat.analyzed++;
      }
      send_msg(lp, LP__DIVING_INFO);
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Processes the messages received by the tree manager
\*===========================================================================*/

char process_messages(tm_prob *tm, int r_bufid)
{
   int bytes, msgtag, sender;
   int lp, cp;
   bc_node *node;

   do{
      bufinfo(r_bufid, &bytes, &msgtag, &sender);

      switch(msgtag){

#ifdef COMPILE_IN_TM
       case FEASIBLE_SOLUTION_NONZEROS:
       case FEASIBLE_SOLUTION_USER:
	 receive_int_array(&(tm->best_sol.xlevel), 1);
	 receive_int_array(&(tm->best_sol.xindex), 1);
	 receive_int_array(&(tm->best_sol.xiter_num), 1);
	 receive_dbl_array(&(tm->best_sol.lpetol), 1);
	 receive_dbl_array(&(tm->best_sol.objval), 1);
	 receive_int_array(&(tm->best_sol.xlength), 1);
	 if (tm->best_sol.xlength > 0){
	    FREE(tm->best_sol.xind);
	    FREE(tm->best_sol.xval);
	    tm->best_sol.xind = (int *) malloc(tm->best_sol.xlength*ISIZE);
	    tm->best_sol.xval = (double *) malloc(tm->best_sol.xlength*DSIZE);
	    receive_int_array(tm->best_sol.xind, tm->best_sol.xlength);
	    receive_dbl_array(tm->best_sol.xval, tm->best_sol.xlength);
	 }
	 if (!tm->has_ub || tm->best_sol.objval < tm->ub){
	    tm->has_ub = TRUE;
	    tm->ub = tm->best_sol.objval;
	 }
	 tm->best_sol.has_sol = TRUE;
	 break;
#endif
	 
       case UPPER_BOUND:
	 process_ub_message(tm);
	 break;

       case LP__IS_FREE:
	 receive_int_array(&cp, 1);
	 tm->stat.chains++;
	 mark_lp_process_free(tm, find_process_index(&tm->lp, sender), cp);
	 break;


       case LP__NODE_DESCRIPTION:
	 node = tm->active_nodes[find_process_index(&tm->lp, sender)];
	 receive_node_desc(tm, node);
	 break;


       case LP__BRANCHING_INFO:
	 node = tm->active_nodes[find_process_index(&tm->lp, sender)];
	 process_branching_info(tm, node);
	 break;


       case LP__CUT_NAMES_REQUESTED:
	 unpack_cut_set(tm, sender, 0, NULL);
	 break;


       case LP__NODE_RESHELVED: /* implies LP__IS_FREE !! */
	 lp = find_process_index(&tm->lp, sender);
	 tm->active_nodes[lp]->node_status = NODE_STATUS__HELD;
	 REALLOC(tm->nextphase_cand, bc_node *,
		 tm->nextphase_cand_size, tm->nextphase_candnum+1, BB_BUNCH);
	 tm->nextphase_cand[tm->nextphase_candnum++] = tm->active_nodes[lp];
	 mark_lp_process_free(tm, lp, tm->active_nodes[lp]->cp);
	 break;


       case LP__NODE_DISCARDED: /* implies LP__IS_FREE !! */
	 lp = find_process_index(&tm->lp, sender);
	 tm->active_nodes[lp]->node_status = NODE_STATUS__PRUNED;
	 mark_lp_process_free(tm, lp, tm->active_nodes[lp]->cp);
	 break;


       case SOMETHING_DIED:
	 printf("Something has died... Halting the machine.\n\n");
	 return(FALSE);


       default:
	 printf("Unknown message type %i\n\n", msgtag);
	 return(FALSE);
      }
      freebuf(r_bufid);
   }while((r_bufid = nreceive_msg(ANYONE, ANYTHING)));
   return(TRUE);
}

/*===========================================================================*/

void process_ub_message(tm_prob *tm)
{
   int s_bufid;
   double new_ub;
   /* A new best solution has been found. The solution is sent
    * to the master, but the bound comes here, too.*/
   receive_dbl_array(&new_ub, 1);
   if ((!tm->has_ub) || (tm->has_ub && new_ub < tm->ub)){
      tm->has_ub = TRUE;
      tm->ub = new_ub;
      s_bufid = init_send(DataInPlace);
      send_dbl_array(&tm->ub, 1);
      msend_msg(tm->lp.procs, tm->lp.procnum, UPPER_BOUND);
      freebuf(s_bufid);
   }
   if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f;
#pragma omp critical(write_vbc_emulation_file)
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 PRINT_TIME(tm, f);
	 fprintf(f, "U %.2f \n", tm->ub);
	 fclose(f); 
      }
   }else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$U %.2f\n", tm->ub);
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Receives and adds a set of cuts
\*===========================================================================*/

void unpack_cut_set(tm_prob *tm, int sender, int cutnum, row_data *rows)
{
   int old_cutnum = tm->cut_num, new_cutnum = cutnum, *itmp, i;
   cut_data **cuts;
#ifndef COMPILE_IN_LP   
   int s_bufid;

   /* If the LP solver exists as a separate process, we have to
      receive the cuts through PVM. Otherwise, we can access them
      directly. */
   receive_int_array(&new_cutnum, 1);
#endif
   REALLOC(tm->cuts, cut_data *, tm->allocated_cut_num, old_cutnum +
	   new_cutnum, (old_cutnum / tm->stat.created + 5) * BB_BUNCH);
   cuts = tm->cuts;
   REMALLOC(tm->tmp.i, int, tm->tmp.i_size, new_cutnum, BB_BUNCH);
   itmp = tm->tmp.i;
   tm->cut_num += new_cutnum;
   for (i = 0; i < new_cutnum; i++){
#ifdef COMPILE_IN_LP
      cuts[(itmp[i] = old_cutnum + i)] = rows[i].cut;
      cuts[itmp[i]]->name = itmp[i];
#else   
      cuts[(itmp[i] = old_cutnum + i)] = unpack_cut(NULL);
      cuts[itmp[i]]->name = itmp[i];
#endif
   }
#ifndef COMPILE_IN_LP
   if (sender){ /* Do we have to return the names? */
      s_bufid = init_send(DataInPlace);
      send_int_array(itmp, new_cutnum);
      send_msg(sender, LP__CUT_NAMES_SERVED);
   }
#endif
}

/*===========================================================================*/

int receive_lp_timing(tm_prob *tm)
{
   char something_died = FALSE;
#ifndef COMPILE_IN_LP
   int i, r_bufid = 0, msgtag, bytes, sender;
   node_times tim;
   struct timeval timeout = {5, 0};
   double ramp_up_tm = tm->comp_times.ramp_up_tm;
   double ramp_down_time = tm->comp_times.ramp_down_time;
   double start_node = tm->comp_times.start_node;
   int lp, cp;
   bc_node *node;
   
   memset(&tm->comp_times, 0, sizeof(node_times));
   tm->comp_times.ramp_up_tm = ramp_up_tm;
   tm->comp_times.ramp_down_time = ramp_down_time;
   tm->comp_times.start_node = start_node;

   for (i = 0; i < tm->lp.procnum; i++){
      msgtag = 0;
      while (msgtag != LP__TIMING){
	 r_bufid = treceive_msg(tm->lp.procs[i], ANYTHING, &timeout);
	 if (r_bufid > 0){
	    bufinfo(r_bufid, &bytes, &msgtag, &sender);
	    switch(msgtag){
	       
#ifdef COMPILE_IN_TM
	     case FEASIBLE_SOLUTION_NONZEROS:
	     case FEASIBLE_SOLUTION_USER:
	       receive_int_array(&(tm->best_sol.xlevel), 1);
	       receive_int_array(&(tm->best_sol.xindex), 1);
	       receive_int_array(&(tm->best_sol.xiter_num), 1);
	       receive_dbl_array(&(tm->best_sol.lpetol), 1);
	       receive_dbl_array(&(tm->best_sol.objval), 1);
	       receive_int_array(&(tm->best_sol.xlength), 1);
	       if (tm->best_sol.xlength > 0){
		  FREE(tm->best_sol.xind);
		  FREE(tm->best_sol.xval);
		  tm->best_sol.xind =
		     (int *) malloc(tm->best_sol.xlength*ISIZE);
		  tm->best_sol.xval =
		     (double *) malloc(tm->best_sol.xlength*DSIZE);
		  receive_int_array(tm->best_sol.xind, tm->best_sol.xlength);
		  receive_dbl_array(tm->best_sol.xval, tm->best_sol.xlength);
	       }
	       if (!tm->has_ub || tm->best_sol.objval < tm->ub){
		  tm->has_ub = TRUE;
		  tm->ub = tm->best_sol.objval;
	       }
	       tm->best_sol.has_sol = TRUE;
	       break;
#endif
	       
	     case UPPER_BOUND:
	       process_ub_message(tm);
	       break;
	       
	     case LP__IS_FREE:
	       receive_int_array(&cp, 1);
	       tm->stat.chains++;
	       mark_lp_process_free(tm,find_process_index(&tm->lp, sender),cp);
	       break;
	       
	       
	     case LP__NODE_DESCRIPTION:
	       node = tm->active_nodes[find_process_index(&tm->lp, sender)];
	       receive_node_desc(tm, node);
	       break;
	       
	       
	     case LP__BRANCHING_INFO:
	       node = tm->active_nodes[find_process_index(&tm->lp, sender)];
	       process_branching_info(tm, node);
	       break;
	       
	       
	     case LP__CUT_NAMES_REQUESTED:
	       unpack_cut_set(tm, sender, 0, NULL);
	       break;
	       
	       
	     case LP__NODE_RESHELVED: /* implies LP__IS_FREE !! */
	       lp = find_process_index(&tm->lp, sender);
	       tm->active_nodes[lp]->node_status = NODE_STATUS__HELD;
	       REALLOC(tm->nextphase_cand, bc_node *, tm->nextphase_cand_size,
		       tm->nextphase_candnum+1, BB_BUNCH);
	       tm->nextphase_cand[tm->nextphase_candnum++] =
		  tm->active_nodes[lp];
	       mark_lp_process_free(tm, lp, tm->active_nodes[lp]->cp);
	       break;
	       
	       
	     case LP__NODE_DISCARDED: /* implies LP__IS_FREE !! */
	       lp = find_process_index(&tm->lp, sender);
	       tm->active_nodes[lp]->node_status = NODE_STATUS__PRUNED;
	       mark_lp_process_free(tm, lp, tm->active_nodes[lp]->cp);
	       break;
	       
	       
	     case SOMETHING_DIED:
	       printf("Something has died... Halting the machine.\n\n");
	       return(FALSE);

	     case LP__TIMING:
	       receive_char_array((char *)&tim, sizeof(node_times));
	       tm->comp_times.communication    += tim.communication;
	       tm->comp_times.lp               += tim.lp;
	       tm->comp_times.separation       += tim.separation;
	       tm->comp_times.fixing           += tim.fixing;
	       tm->comp_times.pricing          += tim.pricing;
	       tm->comp_times.strong_branching += tim.strong_branching;
	       tm->comp_times.wall_clock_lp    += tim.wall_clock_lp;
	       tm->comp_times.ramp_up_lp       += tim.ramp_up_lp;
	       tm->comp_times.idle_diving      += tim.idle_diving;
	       tm->comp_times.idle_node        += tim.idle_node;
	       tm->comp_times.idle_names       += tim.idle_names;
	       tm->comp_times.idle_cuts        += tim.idle_cuts;
	       tm->comp_times.cut_pool         += tim.cut_pool;
	       break;

	     default:
	       printf("Unknown message type %i\n\n", msgtag);
	       break;
	    }
	    freebuf(r_bufid);
	 }else{
	    if (pstat(tm->lp.procs[i]) != PROCESS_OK){
#if 0
	       /* Probably don't need this */
	       stop_processes(&tm->lp);
	       stop_processes(&tm->cg);
	       stop_processes(&tm->cp);
#endif
	       something_died = TRUE;
	       break;
	    }
	 }
      }
   }
#endif
   
   return(something_died ? FUNCTION_TERMINATED_ABNORMALLY :
	  FUNCTION_TERMINATED_NORMALLY);
}


