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

#include "lp.h"
#include "BB_types.h"
#include "BB_macros.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP functions related to freeing data structures
\*===========================================================================*/

void free_cut(cut_data **cut)
{
   if (*cut){
      if ((*cut)->coef){
	 FREE((*cut)->coef);
      }
      FREE(*cut);
   }
}

/*===========================================================================*/

void free_cuts(cut_data **cuts, int cut_num)
{
   int i;
   if (cuts)
      for (i=cut_num-1; i>=0; i--)
	 if (cuts[i])
#ifdef COMPILE_IN_LP
	    if (cuts[i]->name < 0 || cuts[i]->branch & CUT_BRANCHED_ON)
#endif
	       free_cut(cuts+i);
}

/*===========================================================================*/

void free_col_set(our_col_set **colset)
{
   if (*colset){
      our_col_set *cols = *colset;
      FREE(cols->rel_lb_ind);
      FREE(cols->rel_ub_ind);
      FREE(cols->userind);
      FREE(cols->objx);
      FREE(cols->lb);
      FREE(cols->ub);
      FREE(cols->matbeg);
      FREE(cols->matind);
      FREE(cols->matval);
      FREE(*colset);
   }
}

/*===========================================================================*/

void free_candidate(branch_obj **cand)
{
   if (*cand){
      branch_obj *can = *cand;
#ifdef COMPILE_FRAC_BRANCHING
      int i;
      for (i = can->child_num-1; i >= 0; i--){
	 if (can->frac_num[i]){
	    FREE(can->frac_ind[i]);
	    FREE(can->frac_val[i]);
	 }
      }
#endif
      free_waiting_row(&(can->row));
#ifndef MAX_CHILDREN_NUM
      FREE(can->sense); FREE(can->rhs); FREE(can->range); FREE(can->branch);
#endif
      FREE(*cand);
   }
}

/*===========================================================================*/

void free_candidate_completely(branch_obj **cand)
{
   if (*cand){
#ifndef MAX_CHILDREN_NUM
      branch_obj *can = *cand;
#endif
#ifndef MAX_CHILDREN_NUM
      FREE(can->objval);
      FREE(can->termcode);
      FREE(can->feasible);
      FREE(can->iterd);
#  ifdef COMPILE_FRAC_BRANCHING
      FREE(can->frac_num); FREE(can->frac_ind); FREE(can->frac_val);
#  endif
#endif
      free_candidate(cand);
   }
}

/*===========================================================================*/

void free_waiting_row(waiting_row **wrow)
{
   waiting_row *wr = *wrow;
   if (wr){
      FREE(wr->matval);
      FREE(wr->matind);
      free_cut(&wr->cut);
      free(wr);
      *wrow = NULL;
   }
}
   
/*===========================================================================*/

void free_waiting_rows(waiting_row **rows, int row_num)
{
   int i;
   if (rows)
      for (i=row_num-1; i>=0; i--)
	 free_waiting_row(rows+i);
}
   
/*===========================================================================*/

void free_waiting_row_array(waiting_row ***rows, int row_num)
{
   free_waiting_rows(*rows, row_num);
   FREE(*rows);
}

/*===========================================================================*/

void free_node_desc(node_desc **desc)
{
   if (*desc){
      node_desc *n = *desc;
      FREE(n->cutind.list);
      FREE(n->uind.list);
      if (n->nf_status == NF_CHECK_AFTER_LAST ||
	  n->nf_status == NF_CHECK_UNTIL_LAST)
	 FREE(n->not_fixed.list);
      if (n->basis.basis_exists){
	 FREE(n->basis.basevars.list);
	 FREE(n->basis.basevars.stat);
	 FREE(n->basis.extravars.list);
	 FREE(n->basis.extravars.stat);
	 FREE(n->basis.baserows.list);
	 FREE(n->basis.baserows.stat);
	 FREE(n->basis.extrarows.list);
	 FREE(n->basis.extrarows.stat);
      }
      if (n->desc_size > 0)
	 FREE(n->desc);
      FREE(*desc);
   }
}

/*===========================================================================*/

void free_node_dependent(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   int i;

   free_node_desc(&p->desc);
   for (i = p->base.cutnum; i < lp_data->m; i++){
#ifdef COMPILE_IN_LP
      if (lp_data->rows[i].cut->name < 0 ||
	  lp_data->rows[i].cut->branch & CUT_BRANCHED_ON)
#endif
	 free_cut(&lp_data->rows[i].cut);
#ifdef COMPILE_IN_LP
      else
	 lp_data->rows[i].cut = NULL;
#endif
   }
   if (p->par.branch_on_cuts && p->slack_cut_num > 0){
      free_cuts(p->slack_cuts, p->slack_cut_num);
      p->slack_cut_num = 0;
   }
   unload_lp_prob(lp_data);
}

/*===========================================================================*/

void free_lp(lp_prob *p)
{
   int i;

   free_prob_dependent_u(p);
   free_waiting_row_array(&p->waiting_rows, p->waiting_row_num);
   for (i = p->lp_data->maxn - 1; i >= 0; i--)
      FREE(p->lp_data->vars[i]);
   FREE(p->lp_data->vars);
#ifdef COMPILE_IN_LP
   for (i = p->base.cutnum - 1; i >= 0; i--)
      free_cut(&(p->lp_data->rows[i].cut));
   free_node_desc(&p->desc);
#else
   for (i = p->lp_data->m - 1; i >= 0; i--)
      free_cut(&(p->lp_data->rows[i].cut));
   FREE(p->bdesc);
#endif
   FREE(p->lp_data->rows);
   close_lp_solver(p->lp_data);
   free_lp_arrays(p->lp_data);
   free_mip_desc(p->lp_data->mip);
   FREE(p->lp_data->mip);
   FREE(p->lp_data);
#if !(defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP))
   free_mip_desc(p->mip);
   FREE(p->mip);
#endif
   FREE(p->base.userind);
   FREE(p->best_sol.xind);
   FREE(p->best_sol.xval);
   if (p->par.branch_on_cuts){
      FREE(p->slack_cuts);
   }
   FREE(p->obj_history);
   FREE(p);
}

/*===========================================================================*/

void free_mip_desc(MIPdesc *mip)
{
   int j;
   
   FREE(mip->matbeg);
   FREE(mip->matind);
   FREE(mip->matval);
   FREE(mip->obj);
   FREE(mip->rhs);
   FREE(mip->rngval);
   FREE(mip->sense);
   FREE(mip->lb);
   FREE(mip->ub);
   FREE(mip->is_int);
   if (mip->colname){
      for (j = 0; j < mip->n; j++){
	 FREE(mip->colname[j]);
      }
      FREE(mip->colname);
   }
}

/*===========================================================================*/

void free_lp_arrays(LPdata *lp_data)
{
   int i;
   
   FREE(lp_data->not_fixed);
   FREE(lp_data->status);
   FREE(lp_data->x);
   FREE(lp_data->dj);
   FREE(lp_data->dualsol);
   FREE(lp_data->slacks);
#ifdef __CPLEX__
   FREE(lp_data->lb);
   FREE(lp_data->ub);
#endif
   FREE(lp_data->vars);
   FREE(lp_data->tmp.c);
   FREE(lp_data->tmp.i1);
   FREE(lp_data->tmp.i2);
   FREE(lp_data->tmp.d);
   FREE(lp_data->tmp.p1);
   FREE(lp_data->tmp.p2);
   FREE(lp_data->tmp.cv);
   FREE(lp_data->tmp.iv);
   FREE(lp_data->tmp.dv);
}

