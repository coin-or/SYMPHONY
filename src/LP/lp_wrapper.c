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
#include <memory.h>
#include <math.h>
#include <malloc.h>

#include "lp.h"
#include "proccomm.h"
#include "qsortucb.h"
#include "messages.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "BB_types.h"
#include "lp_solver.h"
#ifdef USE_CGL_CUTS
#include "cg.h"
#endif
#if defined (COMPILE_IN_LP) && defined (COMPILE_IN_TM)
#include "master_u.h"
#endif
#ifdef COMPILE_IN_CP
#include "cp.h"
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP wrapper functions that interface with the user.
\*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_receive_lp_data that
 * receives the initial data from the Master process. Returns TRUE if
 * succeeded, FALSE otherwise.
\*===========================================================================*/

int receive_lp_data_u(lp_prob *p)
{
   int r_bufid;
   char has_desc;
   char has_colnames;
   int i;

   r_bufid = receive_msg(p->master, LP_DATA);
   receive_char_array((char *)(&p->par), sizeof(lp_params));
   receive_char_array(&p->has_ub, 1);
   if (p->has_ub){
      receive_dbl_array(&p->ub, 1);
   }else{
      p->ub = - (MAXDOUBLE / 2);
   }
   receive_int_array(&p->draw_graph, 1);
   receive_int_array(&p->base.varnum, 1);
   if (p->base.varnum > 0){
      p->base.userind = (int *) malloc(p->base.varnum * ISIZE);
      receive_int_array(p->base.userind, p->base.varnum);
   }
   receive_int_array(&p->base.cutnum, 1);
   receive_char_array(&has_desc, 1);
   if (has_desc){
      MIPdesc *mip = p->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
      receive_int_array(&(mip->m), 1);
      receive_int_array(&(mip->n), 1);
      receive_int_array(&(mip->nz), 1);
      receive_char_array(&(mip->obj_sense), 1);
      receive_dbl_array(&(mip->obj_offset), 1);

      /* Allocate memory */
      mip->matbeg = (int *) malloc(ISIZE * (mip->n + 1));
      mip->matval = (double *) malloc(DSIZE * mip->matbeg[mip->n]);
      mip->matind = (int *)    malloc(ISIZE * mip->matbeg[mip->n]);
      mip->obj    = (double *) malloc(DSIZE * mip->n);
      mip->rhs    = (double *) malloc(DSIZE * mip->m);
      mip->sense  = (char *)   malloc(CSIZE * mip->m);
      mip->rngval = (double *) malloc(DSIZE * mip->m);
      mip->ub     = (double *) malloc(DSIZE * mip->n);
      mip->lb     = (double *) malloc(DSIZE * mip->n);
      mip->is_int = (char *)   calloc(CSIZE, mip->n);

      /* Receive the problem description */
      receive_int_array(mip->matbeg, mip->n);
      receive_int_array(mip->matind, mip->nz);
      receive_dbl_array(mip->matval, mip->nz);
      receive_dbl_array(mip->obj, mip->n);
      receive_dbl_array(mip->rhs, mip->m);
      receive_char_array(mip->sense, mip->m);
      receive_dbl_array(mip->rngval, mip->m);
      receive_dbl_array(mip->ub, mip->n);
      receive_dbl_array(mip->lb, mip->n);
      receive_char_array(mip->is_int, mip->n);
      receive_char_array(&has_colnames, 1);
      if (has_colnames){
	 mip->colname = (char **) malloc(sizeof(char *) * mip->n);   
	 for (i = 0; i < mip->n; i++){
	    mip->colname[i] = (char *) malloc(CSIZE * 9);
	    receive_char_array(mip->colname[i], 8);
	    mip->colname[i][8] = 0;
	 }
      }
   }
   
   switch( user_receive_lp_data(&p->user)){
    case USER_ERROR:
      freebuf(r_bufid);
      return(FALSE);
    case USER_SUCCESS:
    case USER_AND_PP: 
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      break;
    default:
      freebuf(r_bufid);
      /* Unexpected return value. Do something!! */
      return(FALSE);
   }
   return(TRUE);
}

/*===========================================================================*/

int comp_cut_name(const void *c0, const void *c1)
{
   return((*((cut_data **)c0))->name - (*((cut_data **)c1))->name);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_create_subproblem that
 * creates the problem matrix.
\*===========================================================================*/

int create_subproblem_u(lp_prob *p)
{
   node_desc *desc = p->desc;

   LPdata *lp_data = p->lp_data;
   int i, j, maxm, maxn, maxnz;
   row_data *row, *rows;

   int bvarnum = p->base.varnum;
   int bcutnum = p->base.cutnum;
   var_desc **vars;

   int *d_uind = NULL, *d_cind = NULL; /* just to keep gcc quiet */

   double *rhs, *rngval;
   char *sense, *status;
   cut_data *cut;
   branch_desc *bobj;

   int new_row_num;
   waiting_row **new_rows;

   int user_res;
   int *userind;

   double *lb, *ub;
   char *is_int;

   lp_data->n = bvarnum + desc->uind.size;
   lp_data->m = p->base.cutnum + desc->cutind.size;

   maxm = lp_data->maxm;
   maxn = lp_data->maxn;
   maxnz = lp_data->maxnz;

   lp_data->nf_status = desc->nf_status;
   if (desc->nf_status == NF_CHECK_AFTER_LAST ||
       desc->nf_status == NF_CHECK_UNTIL_LAST){
      lp_data->not_fixed_num = desc->not_fixed.size;
      memcpy(lp_data->not_fixed, desc->not_fixed.list,
	     lp_data->not_fixed_num * ISIZE);
   }

   if (desc->uind.size > 0){ /* fill up the rest of lp_data->vars */
      if (MAX(maxn, bvarnum) < lp_data->n){
	 lp_data->vars = (var_desc **)
	    realloc(lp_data->vars, lp_data->n * sizeof(var_desc *));
	 vars = lp_data->vars;
	 for (i = lp_data->n - 1; i >= MAX(maxn, bvarnum); i--)
	    vars[i] = (var_desc *) malloc( sizeof(var_desc) );
      }
      vars = lp_data->vars;
      d_uind = desc->uind.list;
      for (i = desc->uind.size - 1; i >= 0; i--){
	 vars[i + bvarnum]->userind = d_uind[i];
	 vars[i + bvarnum]->colind = bvarnum + i;
      }
   }
   lp_data->ordering = COLIND_AND_USERIND_ORDERED;

   lp_data->mip->n = lp_data->n;
   lp_data->mip->m = lp_data->m;


   /* Create the list of indices to pass to the user */
   userind = (int *) malloc(lp_data->n * ISIZE);
   vars = lp_data->vars;
   for (i = lp_data->n - 1; i >= 0; --i){
      userind[i] = vars[i]->userind;
   }
   
   user_res = user_create_subproblem(p->user,
       /* list of base and extra variables */
       userind,
       /* description of the LP relaxation to be filled out by the user */
       lp_data->mip, 
       /* max sizes (estimated by the user) */
       &maxn, &maxm, &maxnz);
   
   switch (user_res){
      
    case USER_DEFAULT:
       
      if (!p->mip->matbeg){
	 printf("Illegal return code.\n");
	 printf("Trying to use default user_create_subproblem without");
	 printf("reading an MPS or AMPL file. Exiting...\n\n");
	 exit(2000);
      }

      lp_data->mip->nz = p->mip->nz;      
      /* Allocate the arrays.*/
      lp_data->mip->matbeg  = (int *) malloc((lp_data->mip->n + 1) * ISIZE);
      lp_data->mip->matind  = (int *) malloc((lp_data->mip->nz) * ISIZE);
      lp_data->mip->matval  = (double *) malloc((lp_data->mip->nz) * DSIZE);
      lp_data->mip->obj     = (double *) malloc(lp_data->mip->n * DSIZE);
      lp_data->mip->ub      = (double *) malloc(lp_data->mip->n * DSIZE);
      lp_data->mip->lb      = (double *) calloc(lp_data->mip->n, DSIZE); 
      lp_data->mip->rhs     = (double *) malloc(lp_data->mip->m * DSIZE);
      lp_data->mip->sense   = (char *)   malloc(lp_data->mip->m * CSIZE);
      lp_data->mip->rngval  = (double *) calloc(lp_data->mip->m, DSIZE);
      lp_data->mip->is_int  = (char *)   calloc(lp_data->mip->n, CSIZE);

      /* Fill out the appropriate data structures*/
      lp_data->mip->matbeg[0] = 0;
      for (i = 0; i < lp_data->mip->n; i++){
	 lp_data->mip->obj[i]        = p->mip->obj[userind[i]];
	 lp_data->mip->ub[i]         = p->mip->ub[userind[i]];
	 lp_data->mip->lb[i]         = p->mip->lb[userind[i]];
	 lp_data->mip->is_int[i]     = p->mip->is_int[userind[i]];
	 lp_data->mip->matbeg[i+1]   = lp_data->mip->matbeg[i] + 
	                                (p->mip->matbeg[userind[i]+1] -
	                                 p->mip->matbeg[userind[i]]);
	 for (j = 0; j < (lp_data->mip->matbeg[i+1]-lp_data->mip->matbeg[i]);
	      j++){
	    lp_data->mip->matind[lp_data->mip->matbeg[i]+j]=
	       p->mip->matind[p->mip->matbeg[userind[i]]+j];
	    lp_data->mip->matval[lp_data->mip->matbeg[i]+j]=
	       p->mip->matval[p->mip->matbeg[userind[i]]+j];
	 }
      }
      for (i = 0; i < p->mip->m; i++){
	 lp_data->mip->rhs[i] = p->mip->rhs[i];
	 lp_data->mip->sense[i] = p->mip->sense[i];
	 lp_data->mip->rngval[i] = p->mip->rngval[i];
      }
      
      maxm = p->mip->m;
      maxn = p->mip->n;
      maxnz = p->mip->nz;

      /* Fall through to next case */

    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
       
      /* User function terminated without problems. In the post-processing
       * the extra cuts are added. HOWEVER, this might not be done until the
       * problem is loaded into the lp solver (for cplex it is not possible).
       * So for now just reset lp_data->m, do everything to load in the
       * stuff into the lp solver then come back to adding the cuts. */
      lp_data->m = p->base.cutnum;
      break;

    case USER_ERROR:
       
      /* Error. The search tree node will not be processed. */
      FREE(userind);
      return(FALSE);
      
    default:
       
      /* Unexpected return value. Do something!! */
      FREE(userind);
      return(FALSE);
   }
   
   FREE(userind); /* No longer needed */

   /*------------------------------------------------------------------------*\
    * Let's see about reallocing...
   \*----------------------------------------------------------------------- */

   if (maxm  < lp_data->m)  maxm  = lp_data->m;
   if (maxn  < lp_data->n)  maxn  = lp_data->n;
   if (maxnz < lp_data->nz) maxnz = lp_data->nz;

   size_lp_arrays(lp_data, FALSE, TRUE, maxm, maxn, maxnz);

   /* Default status of every variable is NOT_FIXED */
   if (bvarnum > 0)
      memset(lp_data->status, NOT_FIXED | BASE_VARIABLE, bvarnum);
   if (bvarnum < lp_data->n)
      memset(lp_data->status + bvarnum, NOT_FIXED, lp_data->n - bvarnum);

   /*------------------------------------------------------------------------*\
    * Set the necessary fields in rows
   \*----------------------------------------------------------------------- */

   rows = lp_data->rows;
   rhs = lp_data->mip->rhs;
   rngval = lp_data->mip->rngval;
   sense = lp_data->mip->sense;
   for (i = bcutnum - 1; i >= 0; i--){
      row = rows + i;
      cut = row->cut;
      cut->rhs = rhs[i];
      cut->range = rngval[i];
      cut->branch = (((cut->sense = sense[i]) != 'E') ?
		     ALLOWED_TO_BRANCH_ON : DO_NOT_BRANCH_ON_THIS_ROW);
      cut->size = 0;
      row->eff_cnt = 1;
      row->free = FALSE;
      cut->name = BASE_CONSTRAINT;
   }

   /*------------------------------------------------------------------------*\
    * Set the upper and lower bounds and integer status
   \*----------------------------------------------------------------------- */

   lb = lp_data->mip->lb;
   ub = lp_data->mip->ub;
   is_int = lp_data->mip->is_int;

   vars = lp_data->vars;
   for (i = lp_data->n - 1; i >= 0; i--){
      vars[i]->lb = lb[i];
      vars[i]->ub = ub[i];
      vars[i]->is_int = is_int[i];
   }
   
   /*------------------------------------------------------------------------*\
    * Load the lp problem (load_lp is an lp solver dependent routine).
   \*----------------------------------------------------------------------- */

   load_lp_prob(lp_data, p->par.scaling, p->par.fastmip);

   /* Free the user's description */
   free_mip_desc(lp_data->mip);

   if (desc->cutind.size > 0){
      unpack_cuts_u(p, CUT_FROM_TM, UNPACK_CUTS_SINGLE,
		    desc->cutind.size, desc->cuts, &new_row_num, &new_rows);
      add_row_set(p, new_rows, new_row_num);
      FREE(new_rows);
   }

   /* We don't need the cuts anymore. Free them. */
   if (desc->cutind.size > 0){
#ifndef COMPILE_IN_LP /*If we are using shared memory, we don't need to free*/
      free_cuts(desc->cuts, desc->cutind.size);
#endif 
      FREE(desc->cuts);
   }else{
      desc->cuts = NULL;
   }

   /*------------------------------------------------------------------------*\
    * Now go through the branching stuff
   \*----------------------------------------------------------------------- */

   d_cind = desc->cutind.list;
   vars = lp_data->vars;
   rows = lp_data->rows;
   if (p->bc_level){
      status = lp_data->status;
      for (i = 0; i < p->bc_level; i++){
	 bobj = p->bdesc + i;
	 if (bobj->type == BRANCHING_VARIABLE){
	    j = bobj->name < 0 ? /* base variable : extra variable */
	       -bobj->name-1 :
	       bfind(bobj->name, d_uind, desc->uind.size) + bvarnum;
	    switch (bobj->sense){
	     case 'E':
	       change_lbub(lp_data, j, bobj->rhs, bobj->rhs);
	       vars[j]->lb = vars[j]->ub = bobj->rhs;
	       break;
	     case 'L':
	       change_ub(lp_data, j, bobj->rhs);
	       vars[j]->ub = bobj->rhs;
	       break;
	     case 'G':
	       change_lb(lp_data, j, bobj->rhs);
	       vars[j]->lb = bobj->rhs;
	       break;
	     case 'R':
	       change_lbub(lp_data, j, bobj->rhs, bobj->rhs + bobj->range);
	       vars[j]->lb = bobj->rhs;
	       vars[j]->ub = bobj->rhs + bobj->range;
	       break;
	    }
	    status[j] |= VARIABLE_BRANCHED_ON;
	 }else{ /* BRANCHING_CUT */
	    j = bobj->name < 0 ? /* base constraint : extra constraint */
	       -bobj->name-1 :
	       bfind(bobj->name, d_cind, desc->cutind.size) + bcutnum;
	    change_row(lp_data, j, bobj->sense, bobj->rhs, bobj->range);
#ifdef COMPILE_IN_LP
	    /* Because these cuts are shared with the treemanager, we have to
	       make a copy before changing them if the LP is compiled in */
	    cut = (cut_data *) malloc(sizeof(cut_data));
	    memcpy((char *)cut, (char *)rows[j].cut, sizeof(cut_data));
	    if (cut->size){
	       cut->coef = (char *) malloc(cut->size);
	       memcpy((char *)cut->coef, (char *)rows[j].cut->coef,
		      cut->size);
	    }
	    rows[j].cut = cut;
#else      
	    cut = rows[j].cut;
#endif
	    cut->rhs = bobj->rhs;
	    cut->range = bobj->range;
	    cut->sense = bobj->sense;
	    cut->branch |= CUT_BRANCHED_ON;
	 }
      }
   }

   /*------------------------------------------------------------------------*\
    * The final step: load in the basis.
    * This is cplex style. sorry about it... Still, it
    * might be ok if {VAR,SLACK}_{B,LB,UB} are properly defined
   \*----------------------------------------------------------------------- */

   if (desc->basis.basis_exists == TRUE){
      int *rstat, *cstat;
      if (desc->basis.extravars.size == 0){
	 cstat = desc->basis.basevars.stat;
      }else if (desc->basis.basevars.size == 0){
	 cstat = desc->basis.extravars.stat;
      }else{ /* neither is zero */
	 cstat = lp_data->tmp.i1; /* n */
	 memcpy(cstat,
		desc->basis.basevars.stat, desc->basis.basevars.size *ISIZE);
	 memcpy(cstat + desc->basis.basevars.size,
		desc->basis.extravars.stat, desc->basis.extravars.size *ISIZE);
      }
      if (desc->basis.extrarows.size == 0){
	 rstat = desc->basis.baserows.stat;
      }else if (desc->basis.baserows.size == 0){
	 rstat = desc->basis.extrarows.stat;
      }else{ /* neither is zero */
	 rstat = lp_data->tmp.i2; /* m */
	 memcpy(rstat,
		desc->basis.baserows.stat, desc->basis.baserows.size *ISIZE);
	 memcpy(rstat + desc->basis.baserows.size,
		desc->basis.extrarows.stat, desc->basis.extrarows.size *ISIZE);
      }
      load_basis(lp_data, cstat, rstat);
   }

   return(TRUE);
}

/*===========================================================================*/

int is_feasible_u(lp_prob *p)
{
#ifndef COMPILE_IN_LP
   int s_bufid;
#endif
   int user_res;
   int feasible;
   double new_ub, true_objval = 0;
   LPdata *lp_data = p->lp_data;
   double lpetol = lp_data->lpetol, lpetol1 = 1 - lpetol;
   int *indices;
   double *values, valuesi;
   int cnt, i;

   get_x(lp_data); /* maybe just fractional -- parameter ??? */

   indices = lp_data->tmp.i1; /* n */
   values = lp_data->tmp.d; /* n */

   cnt = collect_nonzeros(p, lp_data->x, indices, values);

   user_res = user_is_feasible(p->user, lpetol, cnt, indices, values,
			       &feasible, &true_objval);
   switch (user_res){
    case USER_ERROR: /* Error. Consider as feasibility not recognized. */
      return(FALSE);
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    case USER_DEFAULT: /* set the default */
      user_res = p->par.is_feasible_default;
      break;
    default:
      break;
   }

   switch (user_res){
    case TEST_ZERO_ONE: /* User wants us to test 0/1 -ness. */
      for (i=cnt-1; i>=0; i--)
	 if (values[i] < lpetol1) break;
      feasible = i < 0 ? IP_FEASIBLE : IP_INFEASIBLE;
      break;
    case TEST_INTEGRALITY:
      for (i = lp_data->n - 1; i >= 0; i--){
	 if (!lp_data->vars[i]->is_int)
	    continue; /* Not an integer variable */
	 valuesi = lp_data->x[i];
	 if (valuesi-floor(valuesi) > lpetol && ceil(valuesi)-valuesi > lpetol)
	    break;
      }
      feasible = i < 0 ? IP_FEASIBLE : IP_INFEASIBLE;
      break;
    default:
      break;
   }

   if (feasible == IP_FEASIBLE){
      /* Send the solution value to the treemanager */
      new_ub = true_objval > 0 ? true_objval : lp_data->objval;
      if (!p->has_ub || new_ub < p->ub){
	 p->has_ub = TRUE;
	 p->ub = new_ub;
	 if (p->par.set_obj_upper_lim)
	    set_obj_upper_lim(p->lp_data, p->ub - p->par.granularity);
	 p->best_sol.xlevel = p->bc_level;
	 p->best_sol.xindex = p->bc_index;
	 p->best_sol.xiter_num = p->iter_num;
	 p->best_sol.xlength = cnt;
	 p->best_sol.lpetol = lpetol;
	 p->best_sol.objval = new_ub;
	 FREE(p->best_sol.xind);
	 FREE(p->best_sol.xval);
	 p->best_sol.xind = (int *) malloc(cnt*ISIZE);
	 p->best_sol.xval = (double *) malloc(cnt*DSIZE);
	 memcpy((char *)p->best_sol.xind, (char *)indices, cnt*ISIZE);
	 memcpy((char *)p->best_sol.xval, (char *)values, cnt*DSIZE);
#ifdef COMPILE_IN_LP
	 p->tm->has_ub = TRUE;
	 p->tm->ub = p->ub;
	 p->tm->opt_thread_num = p->proc_index;
	 if (p->tm->par.vbc_emulation == VBC_EMULATION_FILE){
	    FILE *f;
#pragma omp critical(write_vbc_emulation_file)
	    if (!(f = fopen(p->tm->par.vbc_emulation_file_name, "a"))){
	       printf("\nError opening vbc emulation file\n\n");
	    }else{
	       PRINT_TIME(p->tm, f);
	       fprintf(f, "U %.2f\n", p->ub);
	       fclose(f); 
	    }
	 }else if (p->tm->par.vbc_emulation == VBC_EMULATION_LIVE){
	    printf("$U %.2f\n", p->ub);
	 }
#else
	 s_bufid = init_send(DataInPlace);
	 send_dbl_array(&new_ub, 1);
	 send_msg(p->tree_manager, UPPER_BOUND);
	 freebuf(s_bufid);
#endif
	 PRINT(p->par.verbosity, -1,
	       ("\n****** Found Better Feasible Solution !\n"));
	 if (p->mip->obj_sense == MAXIMIZE){
	    PRINT(p->par.verbosity, -1, ("****** Cost: %f\n\n", -new_ub 
					 + p->mip->obj_offset));
	 }else{
	    PRINT(p->par.verbosity, -1, ("****** Cost: %f\n\n", new_ub
					 + p->mip->obj_offset));
	 }
      }else{
	 PRINT(p->par.verbosity, -1,
	       ("\n* Found Another Feasible Solution.\n"));
	 if (p->mip->obj_sense == MAXIMIZE){
	    PRINT(p->par.verbosity, -1, ("* Cost: %f\n\n", -new_ub
					 + p->mip->obj_offset));
	 }else{
	    PRINT(p->par.verbosity, -1, ("****** Cost: %f\n\n", new_ub
					 + p->mip->obj_offset));
	 }
      }
#ifndef COMPILE_IN_TM
      send_feasible_solution_u(p, p->bc_level, p->bc_index, p->iter_num,
			       lpetol, new_ub, cnt, indices, values);
#endif
      display_lp_solution_u(p, DISP_FEAS_SOLUTION);
      lp_data->termcode = LP_OPT_FEASIBLE;
   }

   return(feasible);
}

/*===========================================================================*/

void send_feasible_solution_u(lp_prob *p, int xlevel, int xindex,
			      int xiter_num, double lpetol, double new_ub,
			      int cnt, int *xind, double *xval)
{
   int s_bufid, msgtag, user_res;

   /* Send to solution to the master */
   s_bufid = init_send(DataInPlace);
   send_int_array(&xlevel, 1);
   send_int_array(&xindex, 1);
   send_int_array(&xiter_num, 1);
   send_dbl_array(&lpetol, 1);
   send_dbl_array(&new_ub, 1);
   send_int_array(&cnt, 1);
   if (cnt > 0){
      send_int_array(xind, cnt);
      send_dbl_array(xval, cnt);
   }
   user_res = user_send_feasible_solution(p->user, lpetol, cnt, xind, xval);
   switch (user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    case USER_ERROR: /* Error. Do the default */
    case USER_DEFAULT: /* set the default */
	 user_res = p->par.send_feasible_solution_default;
      break;
   }
   switch (user_res){
    case SEND_NONZEROS:
      msgtag = FEASIBLE_SOLUTION_NONZEROS;
      break;
    default: /* Otherwise the user packed it */
      msgtag = FEASIBLE_SOLUTION_USER;
      break;
   }
   send_msg(p->master, msgtag);
   freebuf(s_bufid);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_display_solution
 * that (graphically) displays the current solution. 
\*===========================================================================*/

void display_lp_solution_u(lp_prob *p, int which_sol)
{
   int user_res;
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol;

   int number = 0;
   int i, *xind = lp_data->tmp.i1; /* n */
   double tmpd, *xval = lp_data->tmp.d; /* n */

   number = collect_nonzeros(p, x, xind, xval);

   /* Invoke user written function. */
   user_res = user_display_lp_solution(p->user, which_sol, number, xind, xval);
   
   switch(user_res){
    case USER_ERROR:
      /* SYMPHONY ignores error message. */
      return;
    case USER_AND_PP:
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      return;
    case USER_DEFAULT:
      user_res = p->par.display_solution_default;
      break;
    default:
      break;
   }

   switch(user_res){
    case DISP_NOTHING:
      break;
    case DISP_NZ_INT:
      if (p->mip->colname){ 
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" Column names and values of nonzeros in the solution\n");
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    printf("%8s %10.7f\n", p->mip->colname[xind[i]], xval[i]);
	 }
	 printf("\n");
      }else{
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" User indices and values of nonzeros in the solution\n");
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    printf("%7d %10.7f\n", xind[i], xval[i]);
	 }
	 printf("\n");
      }
      break;
    case DISP_NZ_HEXA:
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf(" User indices (hexa) and values of nonzeros in the solution\n");
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      for (i = 0; i < number; ){
	 printf("%7x %10.7f ", xind[i], xval[i]);
	 if (!(++i & 3)) printf("\n"); /* new line after every four pair*/
      }
      printf("\n");
      break;
    case DISP_FRAC_INT:
      if (p->mip->colname){ 
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" Column names and values of fractional vars in the solution\n");
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    tmpd = xval[i];
	    if ((tmpd > floor(tmpd)+lpetol) && (tmpd < ceil(tmpd)-lpetol)){
	       printf("%8s %10.7f\n", p->mip->colname[xind[i]], tmpd);
	    }
	 }
	 printf("\n");
      }else{
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" User indices and values of fractional vars in the solution\n");
	 printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; ){
	    tmpd = xval[i];
	    if ((tmpd > floor(tmpd)+lpetol) && (tmpd < ceil(tmpd)-lpetol)){
	       printf("%7d %10.7f ", xind[i], tmpd);
	       if (!(++i & 3)) printf("\n"); /* new line after every four pair*/
	    }
	 }
      }
      printf("\n");
      break;
    case DISP_FRAC_HEXA:
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf(" User indices (hexa) and values of frac vars in the solution\n");
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      for (i = 0; i < number; ){
	 tmpd = xval[i];
	 if ((tmpd > floor(tmpd)+lpetol) && (tmpd < ceil(tmpd)-lpetol)){
	    printf("%7x %10.7f ", xind[i], tmpd);
	    if (!(++i & 3)) printf("\n"); /* new line after every four pair*/
	 }
      }
      printf("\n");
      break;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_branch that selects
 * candidates to branch on. It receives a number of arguments:
 *  sim_num : slacks in matrix number (the 
\*===========================================================================*/

int select_candidates_u(lp_prob *p, int *cuts, int *new_vars,
			int *cand_num, branch_obj ***candidates)
{
   int user_res, action = USER__BRANCH_IF_MUST;
   LPdata *lp_data = p->lp_data;
   row_data *rows = lp_data->rows;
   double lpetol = lp_data->lpetol;
   int i, j = 0, m = lp_data->m;
   int *candidate_rows;
   branch_obj *can;
   cut_data **slacks_in_matrix = NULL; /* just to keep gcc quiet */

   /* If the user might need to generate rows, we better have the
    * columns COLIND_ORDERED */
   colind_sort_extra(p);

   candidate_rows = lp_data->tmp.i2; /* m */
   if (p->par.branch_on_cuts){
      slacks_in_matrix = (cut_data **)lp_data->tmp.p2; /* m */
      /* get a list of row that are candidates for branching */
      for (i=0; i<m; i++){ /* can't branch on original rows */
	 if ((rows[i].cut->branch & CANDIDATE_FOR_BRANCH)){
	    slacks_in_matrix[j] = rows[i].cut;
	    candidate_rows[j++] = i;
	 }
      }
   }

   /* First decide if we are going to branch or not */
   user_res = user_shall_we_branch(p->user, lpetol, *cuts, j, slacks_in_matrix,
				   p->slack_cut_num, p->slack_cuts, lp_data->n,
				   lp_data->vars, lp_data->x, lp_data->status, 
				   cand_num, candidates, &action);
   switch (user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    case USER_ERROR:   /* In case of error, default is used. */
    case USER_DEFAULT:
      action = p->par.shall_we_branch_default;
      break;
   }

   if (p->bc_level <= p->par.load_balance_level &&
       p->node_iter_num >= p->par.load_balance_iterations)
      action = USER__DO_BRANCH;

   if ((action == USER__DO_NOT_BRANCH) ||
       (action == USER__BRANCH_IF_TAILOFF && *cuts > 0 && !check_tailoff(p)) ||
       (action == USER__BRANCH_IF_MUST && *cuts > 0))
      return(DO_NOT_BRANCH);

   action = col_gen_before_branch(p, new_vars);
   /* vars might have been added, so tmp arrays might be freed/malloc'd,
      but only those where maxn plays any role in the size. Therefore tmp.i2
      and tmp.p2 does NOT change. Phew... */

   if (action == DO_NOT_BRANCH__FATHOMED)
      return(DO_NOT_BRANCH__FATHOMED);

   /* In the other two cases we may have to re-generate the rows
      corresponding to slacks not in the matrix (whether violated
      slacks or branching candidates), depending on new_vars */
   if (*new_vars > 0 && *cand_num > 0){
      cut_data **regen_cuts = (cut_data **) malloc(*cand_num*sizeof(cut_data));
      for (j = 0, i = 0; i < *cand_num; i++){
	 can = (*candidates)[i];
	 if (can->type == VIOLATED_SLACK ||
	     can->type == CANDIDATE_CUT_NOT_IN_MATRIX){
	    regen_cuts[j++] = can->row->cut;
	 }
      }
      if (j > 0){
	 int new_row_num;
	 waiting_row **new_rows;
	 unpack_cuts_u(p, CUT_FROM_TM, UNPACK_CUTS_SINGLE,
		       j, regen_cuts, &new_row_num, &new_rows);
	 for (j = 0, i = 0; i < *cand_num; i++){
	    can = (*candidates)[i];
	    if (can->type == VIOLATED_SLACK ||
		can->type == CANDIDATE_CUT_NOT_IN_MATRIX){
	       free_waiting_row(&can->row);
	       can->row = new_rows[j++];
	    }
	 }
	 FREE(new_rows);
      }
      FREE(regen_cuts);
   }

   if (action == DO_NOT_BRANCH)
      return(DO_NOT_BRANCH);

   /* So the action from col_gen_before_branch is DO_BRANCH */

   action = USER__DO_BRANCH;

   /* OK, so we got to branch */
   user_res = user_select_candidates(p->user, lpetol, *cuts, j,
				     slacks_in_matrix, p->slack_cut_num,
				     p->slack_cuts, lp_data->n, lp_data->vars,
				     lp_data->x, lp_data->status, cand_num,
				     candidates, &action, p->bc_level);
   /* Get rid of any contsraint from slack_cuts which is listed in candidates
    * and rewrite the position of the CANDIDATE_CUT_IN_MATRIX ones */
   if (p->par.branch_on_cuts){
      for (i = 0; i < *cand_num; ){
	 can = (*candidates)[i];
	 switch (can->type){
	  case CANDIDATE_VARIABLE:
	    i++;
	    break;
	  case CANDIDATE_CUT_IN_MATRIX:
	    can->position = candidate_rows[can->position];
	    i++;
	    break;
	  case VIOLATED_SLACK:
	  case CANDIDATE_CUT_NOT_IN_MATRIX:
	    free_cut(p->slack_cuts + can->position);
	    i++;
	    break;
	  case SLACK_TO_BE_DISCARDED:
	    free_cut(p->slack_cuts + can->position);
	    free_candidate(*candidates + i);
	    (*candidates)[i] = (*candidates)[--(*cand_num)];
	    break;
	 }
      }
      compress_slack_cuts(p);
   }

   if (action == USER__DO_NOT_BRANCH)
      return(DO_NOT_BRANCH);

   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      if (! *cand_num){
	 printf("Error! User didn't select branching candidates!!!\n");
	 exit(-1);
      }
      return(DO_BRANCH);
    case USER_ERROR:    /* In case of error, default is used. */
    case USER_DEFAULT:
      user_res = p->par.select_candidates_default;
      break;
    default:
      break;
   }

   i = (int) (p->par.strong_branching_cand_num_max -
      p->par.strong_branching_red_ratio * p->bc_level);
   i = MAX(i, p->par.strong_branching_cand_num_min);

   switch(user_res){
    case USER__CLOSE_TO_HALF:
      branch_close_to_half(i, cand_num, candidates);
      break;
    case USER__CLOSE_TO_HALF_AND_EXPENSIVE:
      branch_close_to_half_and_expensive(i, cand_num, candidates);
      break;
    case USER__CLOSE_TO_ONE_AND_CHEAP:
      branch_close_to_one_and_cheap(i, cand_num, candidates);
      break;

    default:
      /* Unexpected return value. Do something!! */
      break;
   }

   if (! *cand_num){
      PRINT(p->par.verbosity, 2,
	    ("No branching candidates found using default rule...\n"));
      return(DO_NOT_BRANCH);
   }
   return(DO_BRANCH);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_compare_candidates
 * that compares to branching candidates. 
\*===========================================================================*/

int compare_candidates_u(lp_prob *p, double oldobjval,
			 branch_obj *best, branch_obj *can)
{
   int user_res;
   int i;
   double low0, low1, high0, high1;
#ifdef COMPILE_FRAC_BRANCHING
   int frl0, frl1, frh0, frh1;
#endif
   for (i = can->child_num-1; i >= 0; i--){
      switch (can->termcode[i]){
       case LP_OPTIMAL:
	 break;
       case LP_OPT_FEASIBLE:
       case LP_D_UNBOUNDED:
       case LP_D_OBJLIM:
	 can->objval[i] = MAXDOUBLE / 2;
	 break;
       case LP_D_ITLIM:
	 can->objval[i] = MAX(can->objval[i], oldobjval);
	 break;
       case LP_D_INFEASIBLE:
       case LP_ABANDONED:
	 can->objval[i] = oldobjval;
	 break;
      }
   }

   /*------------------------------------------------------------------------*\
    * If ALL descendants in cand terminated with primal infeasibility
    * or high cost, that proves that the current node can be fathomed,
    * so we select cand and force branching on it.
    *
    * MAYBE THIS SHOULD BE LEFT TO THE USER ?????????????????
   \*------------------------------------------------------------------------*/

   for (i = can->child_num-1; i >= 0; i--)
      if (! (can->termcode[i] == LP_D_UNBOUNDED ||
	     can->termcode[i] == LP_D_OBJLIM ||
	     can->termcode[i] == LP_OPT_FEASIBLE ||
	     (can->termcode[i] == LP_OPTIMAL && p->has_ub &&
	      can->objval[i] > p->ub - p->par.granularity))) break;

   if (i < 0){
      /* i.e., we did not break, i.e., we'll select this cand */
      return(SECOND_CANDIDATE_BETTER_AND_BRANCH_ON_IT);
   }

   /* If this is the first, keep it */
   if (best == NULL){
      return(SECOND_CANDIDATE_BETTER);
   }

   /* Otherwise, first give the choice to the user */

   user_res = user_compare_candidates(p->user, best, can, p->ub,
				      p->par.granularity, &i);

   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
       /* User function terminated without problems. No post-processing. */
      return(i);
    case USER_ERROR:
      /* In case of error, default is used. */
    case USER_DEFAULT:
      user_res = p->par.compare_candidates_default;
      break;
    default:
      break;
   }

   /* Well, the user let us make the choice.
    *
    * If something had gone wrong with at least one descendant in
    * can, then prefer to choose something else. */
   for (i = can->child_num-1; i >= 0; i--)
      if (can->termcode[i] == LP_ABANDONED)
	 return(FIRST_CANDIDATE_BETTER);

   /* OK, so all descendants in can finished fine. Just do whatever
    * built-in was asked */
#ifdef COMPILE_FRAC_BRANCHING
   for (frl0 = frh0 = best->frac_num[0], i = best->child_num-1; i; i--){
      frl0 = MIN(frl0, best->frac_num[i]);
      frh0 = MAX(frh0, best->frac_num[i]);
   }
   for (frl1 = frh1 = can->frac_num[0], i = can->child_num-1; i; i--){
      frl1 = MIN(frl1, can->frac_num[i]);
      frh1 = MAX(frh1, can->frac_num[i]);
   }
#endif
   for (low0 = high0 = best->objval[0], i = best->child_num-1; i; i--){
      low0 = MIN(low0, best->objval[i]);
      high0 = MAX(high0, best->objval[i]);
   }
   for (low1 = high1 = can->objval[0], i = can->child_num-1; i; i--){
      low1 = MIN(low1, can->objval[i]);
      high1 = MAX(high1, can->objval[i]);
   }

   switch(user_res){
    case BIGGEST_DIFFERENCE_OBJ:
      i = (high0 - low0 >= high1 - low1) ? 0 : 1;
      break;
    case LOWEST_LOW_OBJ:
      i = (low0 == low1) ? (high0 <= high1 ? 0 : 1) : (low0 < low1 ? 0 : 1);
      break;
    case HIGHEST_LOW_OBJ:
      i = (low0 == low1) ? (high0 >= high1 ? 0 : 1) : (low0 > low1 ? 0 : 1);
      break;
    case LOWEST_HIGH_OBJ:
      i = (high0 == high1) ? (low0 <= low1 ? 0 : 1) : (high0 < high1 ? 0 : 1);
      break;
    case HIGHEST_HIGH_OBJ:
      i = (high0 == high1) ? (low0 >= low1 ? 0 : 1) : (high0 > high1 ? 0 : 1);
      break;
#ifdef COMPILE_FRAC_BRANCHING
    case HIGHEST_LOW_FRAC:
      i = (frl0 == frl1) ? (frh0 >= frh1 ? 0 : 1) : (frl0 > frl1 ? 0 : 1);
      break;
    case LOWEST_LOW_FRAC:
      i = (frl0 == frl1) ? (frh0 <= frh1 ? 0 : 1) : (frl0 < frl1 ? 0 : 1);
      break;
    case HIGHEST_HIGH_FRAC:
      i = (frh0 == frh1) ? (frl0 >= frl1 ? 0 : 1) : (frh0 > frh1 ? 0 : 1);
      break;
    case LOWEST_HIGH_FRAC:
      i = (frh0 == frh1) ? (frl0 <= frl1 ? 0 : 1) : (frh0 < frh1 ? 0 : 1);
      break;
#endif
    default: /* Unexpected return value. Do something!! */
      break;
   }
   return(i == 0 ? FIRST_CANDIDATE_BETTER : SECOND_CANDIDATE_BETTER);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_select_child that
 * selects one of the candidates after branching for further processing.
\*===========================================================================*/

void select_child_u(lp_prob *p, branch_obj *can, char *action)
{
   int user_res;
   int ind, i;

#ifdef DO_TESTS
   char sense;
   for (i = can->child_num-1; i >= 0; i--){
      sense = can->sense[i];
      if (sense != 'E' && sense != 'L' && sense != 'G' && sense != 'R'){
	 printf("The sense of a child doesn't make sense! (nonexistent)\n\n");
	 exit(-212);
      }
   }
#endif

   user_res = user_select_child(p->user, p->ub, can, action);

   switch(user_res){
    case USER_NO_PP:
    case USER_AND_PP:
      /* User function terminated without problems. Skip post-processing. */
      break;
    case USER_ERROR:
      /* In case of error, default is used. */
    case USER_DEFAULT:
      user_res = p->par.select_child_default;
      break;
    default:
      break;
   }

   switch(user_res){
    case PREFER_LOWER_OBJ_VALUE:
      for (i = can->child_num-1; i >= 0; i--)
	 action[i] = RETURN_THIS_CHILD;
      for (ind = 0, i = can->child_num-1; i; i--){
	 if (can->objval[i] < can->objval[ind])
	    ind = i;
      }
      if (!p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that if the lowest objval child is fathomed then everything is */
      break;

    case PREFER_HIGHER_OBJ_VALUE:
      for (i = can->child_num-1; i >= 0; i--)
	 action[i] = RETURN_THIS_CHILD;
      for (ind = 0, i = can->child_num-1; i; i--){
	 if ((can->objval[i] > can->objval[ind]) &&
	     (! p->has_ub ||
	      (p->has_ub && can->objval[i] < p->ub - p->par.granularity)))
	    ind = i;
      }
      if (! p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that this selects the highest objval child NOT FATHOMED, thus
       * if the highest objval child is fathomed then so is everything */
      break;
      
#ifdef COMPILE_FRAC_BRANCHING
    case PREFER_MORE_FRACTIONAL:
      for (i = can->child_num-1; i >= 0; i--)
	 action[i] = RETURN_THIS_CHILD;
      for (ind = 0, i = can->child_num-1; i; i--){
	 if ((can->frac_num[i] > can->frac_num[ind]) &&
	     (! p->has_ub ||
	      (p->has_ub && can->objval[i] < p->ub - p->par.granularity)))
	    ind = i;
      }
      if (! p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that this selects the most fractional child NOT FATHOMED, thus
       * if that child is fathomed then so is everything */
      break;

    case PREFER_LESS_FRACTIONAL:
      for (i = can->child_num-1; i >= 0; i--)
	 action[i] = RETURN_THIS_CHILD;
      for (ind = 0, i = can->child_num-1; i; i--){
	 if ((can->frac_num[i] < can->frac_num[ind]) &&
	     (! p->has_ub ||
	      (p->has_ub && can->objval[i] < p->ub - p->par.granularity)))
	    ind = i;
      }
      if (! p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that this selects the least fractional child NOT FATHOMED, thus
       * if that child is fathomed then so is everything */
      break;
#endif

    case USER_NO_PP:
    case USER_AND_PP:
      break;

    default:
      /* Unexpected return value. Do something!! */
      break;
   }

   /* Throw out the fathomable ones. */
   if (p->lp_data->nf_status == NF_CHECK_NOTHING && p->has_ub){
      for (ind = 0, i = can->child_num-1; i >= 0; i--)
	 if (can->objval[i] > p->ub - p->par.granularity)
	    action[i] = PRUNE_THIS_CHILD_FATHOMABLE;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This function prints whatever statistics we want on branching
\*===========================================================================*/

void print_branch_stat_u(lp_prob *p, branch_obj *can, char *action)
{
   int i;
   
   if (can->type == CANDIDATE_VARIABLE){
      if (p->mip){
	 if (p->mip->colname){
	    printf("Branching on variable %s \n   children: ",
		   p->mip->colname[p->lp_data->vars[can->position]->userind]);
	 }
      }else{
	 printf("Branching on variable %i ( %i )\n   children: ",
		can->position, p->lp_data->vars[can->position]->userind);
      }
   }else{ /* must be CANDIDATE_CUT_IN_MATRIX */
      printf("Branching on a cut %i\n   children: ",
	     p->lp_data->rows[can->position].cut->name);
   }
   for (i=0; i<can->child_num; i++){
      if (can->objval[i] != MAXDOUBLE / 2){
	 printf("[%.3f, %i,%i]  ",
		can->objval[i], can->termcode[i], can->iterd[i]);
      }else{
	 printf("[-1, %i,%i]  ", can->termcode[i], can->iterd[i]);
      }
   }
   printf("\n");

   if (can->type == CANDIDATE_VARIABLE){
      CALL_USER_FUNCTION( user_print_branch_stat(p->user, can, NULL,
			     p->lp_data->n, p->lp_data->vars, action) );
   }else{
      CALL_USER_FUNCTION( user_print_branch_stat(p->user, can,
			     p->lp_data->rows[can->position].cut,
			     p->lp_data->n, p->lp_data->vars, action) );
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Append additional information to the description of an active node 
 * before it is sent back to the tree manager. 
\*===========================================================================*/

void add_to_desc_u(lp_prob *p, node_desc *desc)
{
   desc->desc_size = 0;
   desc->desc = NULL;

   CALL_USER_FUNCTION( user_add_to_desc(p->user, &desc->desc_size,
					&desc->desc) );
}

/*===========================================================================*/

int same_cuts_u(lp_prob *p, waiting_row *wrow1, waiting_row *wrow2)
{
   int user_res;
   int same_cuts = DIFFERENT_CUTS;
   cut_data *rcut1 = NULL, *rcut2 = NULL;

   user_res = user_same_cuts(p->user, wrow1->cut, wrow2->cut, &same_cuts);
   switch (user_res){
    case USER_SUCCESS:
    case USER_NO_PP:
    case USER_AND_PP:
      break;
    case USER_ERROR: /* Error. Use the default */
    case USER_DEFAULT: /* the only default is to compare byte by byte */
      rcut1 = wrow1->cut;
      rcut2 = wrow2->cut;
      if (rcut1->type != rcut2->type || rcut1->sense != rcut2->sense ||
	  rcut1->size != rcut2->size ||
	  memcmp(rcut1->coef, rcut2->coef, rcut1->size))
	 break; /* if LHS is different, then just break out. */

      /* Otherwise the two cuts have the same left hand side. Test which
       * one is stronger */
      /********* something should be done about ranged constraints ***********/
      /* FIXME! */
      if (rcut1->sense == 'L'){
	 same_cuts = rcut1->rhs > rcut2->rhs - p->lp_data->lpetol ?
	    SECOND_CUT_BETTER : FIRST_CUT_BETTER;
	 break;
      }else if (rcut1->sense == 'G'){
	 same_cuts = rcut1->rhs < rcut2->rhs + p->lp_data->lpetol ?
	    SECOND_CUT_BETTER : FIRST_CUT_BETTER;
	 break;
      }
      same_cuts = wrow1->source_pid < wrow2->source_pid ?
	 SECOND_CUT_BETTER : FIRST_CUT_BETTER;
      break;
   }

   switch(same_cuts){
    case SECOND_CUT_BETTER: /* effective replace the old with the new, then..*/
      same_cuts = SAME_CUTS;
      wrow1->violation += fabs(rcut1->rhs - rcut2->rhs);
      rcut1->rhs = rcut2->rhs;
      rcut1->name = rcut2->name;
    case SAME_CUTS:
    case FIRST_CUT_BETTER:  /* delete the new */
      FREE(rcut2->coef);
      break;

    case DIFFERENT_CUTS:
      break;
   }      

   return(same_cuts);
}

/*===========================================================================*/

void unpack_cuts_u(lp_prob *p, int from, int type,
		   int cut_num, cut_data **cuts,
		   int *new_row_num, waiting_row ***new_rows)
{
   LPdata *lp_data = p->lp_data;
   int user_res;
   int i, j, k, l = 0, nzcnt, real_nzcnt, explicit_row_num = 0;
   int *matbeg, *matind;
   double *matval;
   waiting_row **row_list;
   
   colind_sort_extra(p);

   if (cut_num > 0)
      row_list = (waiting_row **) calloc (cut_num, sizeof(waiting_row *));

   /* First SYMPHONY looks for cut types that it knows */
   for (i = 0; i < cut_num; i++){

      switch (cuts[i]->type){
	 
      case EXPLICIT_ROW:
	 row_list[explicit_row_num] =
	    (waiting_row *) malloc(sizeof(waiting_row));
	 row_list[explicit_row_num]->cut = cuts[i];
	 nzcnt = ((int *) (cuts[i]->coef))[0];
	 matind = (int *) (cuts[i]->coef + ISIZE);
	 matval = (double *) (cuts[i]->coef + (1 + nzcnt) * ISIZE);
	 row_list[explicit_row_num]->matind = (int *) malloc(nzcnt * ISIZE);
	 row_list[explicit_row_num]->matval = (double *) malloc(nzcnt * DSIZE);
	 real_nzcnt = 0;
	 for (j = 0; j < lp_data->n; j++){
	    for (k = 0; k < nzcnt; k++){
	       if (matind[k] == lp_data->vars[j]->userind){
		  row_list[explicit_row_num]->matind[real_nzcnt] = j;
		  row_list[explicit_row_num]->matval[real_nzcnt++] = matval[k];
	       }
	    }
	 }
	 row_list[explicit_row_num++]->nzcnt = real_nzcnt;
	 cuts[i] = NULL;
	 break;

      default: /* A user cut type */
	 if (l != i){
	    cuts[l++] = cuts[i];
	    cuts[i] = NULL;
	 }else{
	    l++;
	 }
	 break;
      }
   }

   *new_row_num = 0;
   user_res = user_unpack_cuts(p->user, from, type,
			       lp_data->n, lp_data->vars,
			       l, cuts, new_row_num, new_rows);

   for (i = 0; i < l; i++){
      if (cuts[i]){
	 (*new_rows)[i]->cut = cuts[i];
	 cuts[i] = NULL;
      }
   }
       
   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
    case USER_DEFAULT:

      /* Combine the user's rows with SYMPHONY's rows */ 
      if (*new_row_num == 0 && explicit_row_num == 0){
	 FREE(row_list);
	 *new_row_num = 0;
	 *new_rows = NULL;
      }else if (*new_row_num == 0 && explicit_row_num > 0){
	 *new_row_num = explicit_row_num;
	 *new_rows = row_list;
      }else if (*new_row_num > 0 && explicit_row_num > 0){
	 if (*new_row_num + explicit_row_num > cut_num){
	    row_list = (waiting_row **) realloc(row_list, *new_row_num +
						explicit_row_num);
	 }
	 for (i = explicit_row_num; i < *new_row_num + explicit_row_num; i++){
	    row_list[i] = (*new_rows)[i - explicit_row_num];
	 }
	 
	 FREE(*new_rows);
	 *new_rows = row_list;
      }else{
	 FREE(row_list);
      }
      
      break;

    case USER_ERROR: /* Error. ??? what will happen ??? */
      *new_row_num = 0;
      FREE(*new_rows);

      break;
       
    default: /* No builtin possibility. Counts as ERROR. */
      break;
   }

   free_cuts(cuts, cut_num);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user packs together and sends a message to the cut generator or
 * cut pool process to obtain violated cuts.
 * Default options: SEND_NONZEROS, SEND_FRACTIONS.
 * The function return 1 or 0, depending on whether the sending of the
 * lp solution was successful or not.
\*===========================================================================*/

int send_lp_solution_u(lp_prob *p, int tid)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   int user_res, nzcnt, s_bufid, msgtag = ANYTHING;
   int *xind = lp_data->tmp.i1; /* n */
   double *xval = lp_data->tmp.d; /* n */

   s_bufid = init_send(DataInPlace);
   send_int_array(&p->bc_level, 1);
   send_int_array(&p->bc_index, 1);
   send_int_array(&p->iter_num, 1);
   send_dbl_array(&lp_data->lpetol, 1);
   if (tid == p->cut_gen){
      send_dbl_array(&lp_data->objval, 1);
      send_char_array(&p->has_ub, 1);
      if (p->has_ub)
	 send_dbl_array(&p->ub, 1);
   }
   colind_sort_extra(p);
   user_res = user_send_lp_solution(p->user, lp_data->n, lp_data->vars, x,
				    tid == p->cut_gen ?
				    LP_SOL_TO_CG : LP_SOL_TO_CP);
   switch (user_res){
    case USER_ERROR: /* Error. Consider as couldn't send to cut_gen, i.e.,
		   equivalent to NO_MORE_CUTS_FOUND */
      freebuf(s_bufid);
      return(0);
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      msgtag = LP_SOLUTION_USER;
      break;
    case USER_DEFAULT: /* set the default */
      user_res = p->par.pack_lp_solution_default; /* SEND_NONZEROS */
      break;
   }

   if (msgtag == LP_SOLUTION_USER){
      send_msg(tid, LP_SOLUTION_USER);
      freebuf(s_bufid);
      return(1);
   }

   switch(user_res){
    case SEND_NONZEROS:
      nzcnt = collect_nonzeros(p, x, xind, xval);
      msgtag = LP_SOLUTION_NONZEROS;
      break;
    case SEND_FRACTIONS:
      nzcnt = collect_fractions(p, x, xind, xval);
      msgtag = LP_SOLUTION_FRACTIONS;
      break;
   }
   /* send the data */
   send_int_array(&nzcnt, 1);
   send_int_array(xind, nzcnt);
   send_dbl_array(xval, nzcnt);
   send_msg(tid, msgtag);
   freebuf(s_bufid);

   return(1);
}

/*===========================================================================*/

void logical_fixing_u(lp_prob *p)
{
   char *status = p->lp_data->tmp.c; /* n */
   char *lpstatus = p->lp_data->status;
   char *laststat = status + p->lp_data->n;
   int fixed_num = 0, user_res;

   colind_sort_extra(p);
   memcpy(status, lpstatus, p->lp_data->n);

   user_res = user_logical_fixing(p->user, p->lp_data->n, p->lp_data->vars,
				  p->lp_data->x, status, &fixed_num);
   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      if (fixed_num > 0){
	 while (status != laststat) {
	    *lpstatus &= NOT_REMOVABLE;
	    *lpstatus++ |= (*status++ & (NOT_FIXED |
					 TEMP_FIXED_TO_LB | TEMP_FIXED_TO_UB |
					 PERM_FIXED_TO_LB | PERM_FIXED_TO_UB));
	 }
      }
    case USER_DEFAULT:
      break;
   }
}

/*===========================================================================*/

int generate_column_u(lp_prob *p, int lpcutnum, cut_data **cuts,
		      int prevind, int nextind, int generate_what,
		      double *colval, int *colind, int *collen, double *obj,
		      double *lb, double *ub)
{
   int real_nextind = nextind;
   CALL_USER_FUNCTION( user_generate_column(p->user, generate_what,
					    p->lp_data->m - p->base.cutnum,
					    cuts, prevind, nextind,
					    &real_nextind, colval, colind,
					    collen, obj, lb, ub) );
   return(real_nextind);
}

/*===========================================================================*/

void generate_cuts_in_lp_u(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   int user_res, new_row_num = 0;
   waiting_row **new_rows = NULL;
   char deleted_cut;
   cut_data **cuts;
   
   int i, j;
   waiting_row **wrows = p->waiting_rows;
   
   colind_sort_extra(p);
   
#if defined(COMPILE_IN_CG) || defined(COMPILE_IN_CP) 
   {
#ifdef COMPILE_IN_CP
      int cp_new_row_num = 0;
      waiting_row **cp_new_rows = NULL;
#endif
#ifdef COMPILE_IN_CG
      int cg_new_row_num = 0;
      waiting_row **cg_new_rows = NULL;
#endif
      int user_res2, xlength = 0, *xind = NULL;
      lp_sol *cur_sol = &(p->cgp->cur_sol);
      double *xval = NULL, lpetol = 0;
      
      user_res2 = user_send_lp_solution(p->user,
					lp_data->n, lp_data->vars, x,
					LP_SOL_WITHIN_LP);
      
      if (user_res2 == USER_DEFAULT)
	 user_res2 = p->par.pack_lp_solution_default;
      
      switch (user_res2){
       case USER_ERROR: 
	 return;
       case USER_SUCCESS:
       case USER_AND_PP:
       case USER_NO_PP:
	 break;
       case SEND_NONZEROS:
       case SEND_FRACTIONS:
	 cur_sol->xind = xind = lp_data->tmp.i1; /* n */
	 cur_sol->xval = xval = lp_data->tmp.d; /* n */
	 cur_sol->lpetol = lpetol = lp_data->lpetol;
	 cur_sol->xlevel = p->bc_level;
	 cur_sol->xindex = p->bc_index;
	 cur_sol->xiter_num = p->iter_num;
	 cur_sol->objval = lp_data->objval;
	 if (p->has_ub)
	    p->cgp->ub = p->ub;
	 cur_sol->xlength = xlength = user_res2 == SEND_NONZEROS ?
	                             collect_nonzeros(p, x, xind, xval) :
		                     collect_fractions(p, x, xind, xval);
	 break;
      }
#ifdef COMPILE_IN_CG      
      if (p->cgp->par.do_findcuts && !new_row_num)
	 find_cuts_u(p->cgp, p->lp_data, &cg_new_row_num);
#endif

      if (cg_new_row_num){
	 unpack_cuts_u(p, CUT_FROM_CG, UNPACK_CUTS_MULTIPLE,
		       p->cgp->cuts_to_add_num, p->cgp->cuts_to_add,
		       &cg_new_row_num, &cg_new_rows);
	 p->cgp->cuts_to_add_num = 0;
	 if (cg_new_row_num){
	    for (i = 0; i < cg_new_row_num; i++){
	       if (cg_new_rows[i]->cut->name != CUT__SEND_TO_CP)
		  cg_new_rows[i]->cut->name = CUT__DO_NOT_SEND_TO_CP;
	       cg_new_rows[i]->source_pid = INTERNAL_CUT_GEN;
	       for (j = p->waiting_row_num - 1; j >= 0; j--){
		  if (same_cuts_u(p, p->waiting_rows[j],
				  cg_new_rows[i]) !=
		      DIFFERENT_CUTS){
		     free_waiting_row(cg_new_rows+i);
		     break;
		  }
	       }
	       if (j < 0){
		  add_new_rows_to_waiting_rows(p, cg_new_rows+i, 1);
	       }
	    }
	    FREE(cg_new_rows);
	 }
      }
#if defined(COMPILE_IN_CP) && defined(COMPILE_IN_LP)
      
      if ((p->iter_num == 1 && (p->bc_level > 0 || p->phase==1)) ||
	  (p->iter_num % p->par.cut_pool_check_freq == 0) ||
	  (!cg_new_row_num)){
	 cut_pool *cp = p->tm->cpp[p->cut_pool];
	 p->comp_times.separation += used_time(&p->tt);
	 cur_sol->lp = 0;
#pragma omp critical(cut_pool)
	 if (cp){
	    cp_new_row_num = check_cuts_u(cp, cur_sol);
	    if (++cp->reorder_count % 10 == 0){
	       delete_duplicate_cuts(cp);
	       order_cuts_by_quality(cp);
	       cp->reorder_count = 0;
	    }
	    if (cp_new_row_num){
	       unpack_cuts_u(p, CUT_FROM_CG, UNPACK_CUTS_MULTIPLE,
			     cp->cuts_to_add_num, cp->cuts_to_add,
			     &cp_new_row_num, &cp_new_rows);
	       cp->cuts_to_add_num = 0;
	    }
	 }
	 if (cp_new_row_num){
	    for (i = 0; i < cp_new_row_num; i++){
	       if (cp_new_rows[i]->cut->name != CUT__SEND_TO_CP)
		  cp_new_rows[i]->cut->name = CUT__DO_NOT_SEND_TO_CP;
	       cp_new_rows[i]->source_pid = INTERNAL_CUT_POOL;
	       for (j = p->waiting_row_num - 1; j >= 0; j--){
		  if (same_cuts_u(p, p->waiting_rows[j],
				  cp_new_rows[i]) !=
		      DIFFERENT_CUTS){
		     free_waiting_row(cp_new_rows+i);
		     break;
		  }
	       }
	       if (j < 0){
		  add_new_rows_to_waiting_rows(p, cp_new_rows+i, 1);
	       }
	    }
	    FREE(cp_new_rows);
	 }
	 p->comp_times.cut_pool += used_time(&p->tt);
      }
#endif
   }
#endif      
   
   user_res = user_generate_cuts_in_lp(p->user, lp_data, lp_data->n,
				       lp_data->vars, x,
				       &new_row_num, &cuts);
   
   switch(user_res){
    case USER_ERROR:
      return;
    case GENERATE_CGL_CUTS:
    case USER_DEFAULT:
      /* Add to the user's list of cuts */
#ifdef USE_CGL_CUTS
      if (p->par.generate_cgl_cuts){
	 generate_cgl_cuts(lp_data, &new_row_num, &cuts);
      }
#endif
      /* Fall through to next case */
      
    case DO_NOT_GENERATE_CGL_CUTS:
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      /* Process the generated cuts */
      if (new_row_num){
	 unpack_cuts_u(p, CUT_FROM_CG, UNPACK_CUTS_MULTIPLE,
		       new_row_num, cuts, &new_row_num, &new_rows);
	 for (i = 0; i < new_row_num; i++){
	    if (new_rows[i]->cut->name != CUT__SEND_TO_CP)
	       new_rows[i]->cut->name = CUT__DO_NOT_SEND_TO_CP;
	    new_rows[i]->source_pid = INTERNAL_CUT_GEN;
	 }
      }
      /* Test whether any of the new cuts are identical to any of
         the old ones. */
      if (p->waiting_row_num && new_row_num){
	 for (i = 0, deleted_cut = FALSE; i < new_row_num - 1;
	      deleted_cut = FALSE){
	    for (j = p->waiting_row_num - 1; j >= 0; j--){
	       if (same_cuts_u(p, wrows[j], new_rows[i]) !=
		   DIFFERENT_CUTS){
		  free_waiting_row(new_rows+i);
		  new_rows[i] = new_rows[--new_row_num];
		  deleted_cut = TRUE;
		  break;
	       }
	    }
	    if (!deleted_cut) i++;
	 }
      }
      if (new_row_num){
	 add_new_rows_to_waiting_rows(p, new_rows, new_row_num);
	 FREE(new_rows);
      }
      return;
    default:
      /* Unexpected return value. Do something!! */
      return;
   }      
}

/*===========================================================================*/

void print_stat_on_cuts_added_u(lp_prob *p, int added_rows)
{
   int user_res;
   
   user_res = user_print_stat_on_cuts_added(p->user, added_rows,
					    p->waiting_rows);
   switch(user_res){
    case USER_ERROR:
    case USER_DEFAULT:
      /* print out how many cuts have been added */
      PRINT(p->par.verbosity, 5,
	    ("Number of cuts added to the problem: %i\n", added_rows));
      break;
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }      
}

/*===========================================================================*/

void purge_waiting_rows_u(lp_prob *p)
{
   int user_res, i, j;
   waiting_row **wrows = p->waiting_rows;
   int wrow_num = p->waiting_row_num;
   char *delete_rows;

   REMALLOC(p->lp_data->tmp.cv, char, p->lp_data->tmp.cv_size, wrow_num,
	    BB_BUNCH);
   delete_rows = p->lp_data->tmp.cv; /* wrow_num */

   memset(delete_rows, 0, wrow_num);
   
   user_res = user_purge_waiting_rows(p->user, wrow_num, wrows, delete_rows);
   switch (user_res){
    case USER_ERROR: /* purge all */
      free_waiting_rows(wrows, wrow_num);
      p->waiting_row_num = 0;
      break;
    case USER_DEFAULT: /* the only default is to keep enough for one
			  iteration */
      if (wrow_num - p->par.max_cut_num_per_iter > 0){
	 free_waiting_rows(wrows + p->par.max_cut_num_per_iter,
			   wrow_num-p->par.max_cut_num_per_iter);
	 p->waiting_row_num = p->par.max_cut_num_per_iter;
      }
      break;
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      for (i = j = 0; i < wrow_num; i++){
	 if (delete_rows[i]){
	    free_waiting_row(wrows + i);
	 }else{
	    wrows[j++] = wrows[i];
	 }
      }
      p->waiting_row_num = j;
      break;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
}

/*===========================================================================*/
/*===========================================================================*\
 * This function invokes the user written function user_free_prob_dependent
 * that deallocates the user defined part of the data structure. Returns TRUE
 * if succeeded, FALSE otherwise.
\*===========================================================================*/

void free_prob_dependent_u(lp_prob *p)
{
   switch (user_free_lp(&p->user)){
    case USER_ERROR:
      /* BlackBox ignores error message */
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      return;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
}

/*===========================================================================*/

