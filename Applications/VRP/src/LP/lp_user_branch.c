/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2001 Ted Ralphs. All Rights Reserved.           */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <math.h>
#include <malloc.h>
#include <memory.h>
#include <stdio.h>

#include "BB_constants.h"
#include "BB_macros.h"
#include "qsortucb.h"
#include "lp_user.h"
#include "lp_u.h"
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#include "lp.h"
/*___END_EXPERIMENTAL_SECTION___*/
#include "vrp_macros.h"
#include "vrp_const.h"
#include "timemeas.h"
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#include "lp.h"
/*___END_EXPERIMENTAL_SECTION___*/

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions for the LP process related
 * to branching.
\*===========================================================================*/

/*===========================================================================*\
 * This function determines whether to branch. You can eseentially
 * leave this up to SYMPHONY unless there is some compelling reason not to.
\*===========================================================================*/

int user_shall_we_branch(void *user, double lpetol, int cutnum,
			 int slacks_in_matrix_num, cut_data **slacks_im_matrix,
			 int slack_cut_num, cut_data **slack_cuts, int varnum,
			 var_desc **vars, double *x, char *status,
			 int *cand_num, branch_obj ***candidates,
			 int *action)
{
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#if 0
   lp_prob *p = get_lp_ptr();
   double *obj_hist = p->obj_history;
   int i;
   int backsteps = p->par.tailoff_obj_backsteps;
#endif
/*___END_EXPERIMENTAL_SECTION___*/
   vrp_spec *vrp = (vrp_spec *) user;

   if (!vrp->par.detect_tailoff){
      *action = USER__BRANCH_IF_MUST;
   }else{
      *action = USER__BRANCH_IF_TAILOFF;
   }

   return(USER_NO_PP);

/*__BEGIN_EXPERIMENTAL_SECTION__*/
#if 0
   for (i = MIN(p->iter_num-1, backsteps) - 1; i >= 0; i--)
      obj_hist[i+1] = obj_hist[i];
   obj_hist[0] = p->lp_data->objval;

#ifdef DO_TESTS
   /*The solution cannot be integral at this point*/
   for (i = 0; i < varnum; i++){
      xvali = x[i];
      if (xvali - floor(xvali) > lpetol && ceil(xvali) - xvali > lpetol)
	 break;
   }
   if (i >= varnum){
      /* solution is integral , yet infeasible*/
      *action = USER__DO_NOT_BRANCH;
      return(USER_NO_PP);
   }
#endif
   
   /* if the objective function differences are below tailoff_obj_frac
      tailoff_obj_backsteps  times in a row, branch */
   for (i = 1; i <= backsteps; i++)
     if ((obj_hist[i-1]-obj_hist[i]) > p->par.tailoff_obj_frac) break;  

   if  (i < backsteps){
      /* no tailoff */
      *action = USER__BRANCH_IF_MUST;
      return(USER_NO_PP);
   }else{
      *action = USER__DO_BRANCH;
      return(USER_NO_PP);
   }
   return(DEFAULT);
#endif
/*___END_EXPERIMENTAL_SECTION___*/
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we select the branching candidates. This can essentially be
 * left to SYMPHONY too using one of the built-in functions, but here, I
 * demonstrate how to branch on cuts, which must be done by the user.
\*===========================================================================*/

int user_select_candidates(void *user, double lpetol, int cutnum,
			   int slacks_in_matrix_num,
			   cut_data **slacks_in_matrix, int slack_cut_num,
			   cut_data **slack_cuts, int varnum, var_desc **vars,
			   double *x, char *status, int *cand_num,
			   branch_obj ***candidates, int *action,
			   int bc_level)

{
   vrp_spec *vrp = (vrp_spec *)user;
   cut_data *cut;  
   branch_obj **cand_list, *can;
   int i, candnum, found_violated = FALSE;
   p_w_l *pwl;
   double left_hand_side, lpetol1 = 1 - lpetol, lpetol5=.95, slack, fracx; 
   waiting_row **new_rows;
   int new_row_num;
   int *userind;
   int sim_cand_num = 0, sc_cand_num = 0;
   /* numbers of slacks_in_matrix, slack_cuts and variables, that are used
      as branching cands */

   if (!vrp->par.branch_on_cuts && vrp->par.branching_rule == 2)
      /* use the built-in rule */
      return(DEFAULT);

   *cand_num = 0;
   candnum = 0;
   /* allocate also memory for the basic vars */
   *candidates = cand_list = (branch_obj **)
      malloc((varnum + (vrp->par.branch_on_cuts ?
		       (slacks_in_matrix_num + slack_cut_num) : 0)) *
      sizeof(branch_obj *));
   switch (vrp->par.branch_on_cuts){
    case TRUE:
      pwl = (p_w_l *) malloc((slacks_in_matrix_num + slack_cut_num)*
			     sizeof(p_w_l));
      userind = (int *) malloc(varnum*ISIZE);
      
      for (i = varnum - 1; i >= 0; i--)
	 userind[i] = vars[i]->userind;
      
      /* First go through the slack cuts and enlist the violated ones */
      for (i = 0; i < slack_cut_num; i++){
	 left_hand_side = compute_lhs(varnum, userind, x, cut = slack_cuts[i],
				      vrp->vertnum);
	 switch (cut->type){
	  case SUBTOUR_ELIM_SIDE:
	    slack = cut->rhs - left_hand_side;
	    if (slack < -lpetol){/*---------------- This cut became violated */
	       found_violated = TRUE;
	       can = cand_list[candnum++] =
		  (branch_obj *) calloc(1, sizeof(branch_obj));
	       can->type = VIOLATED_SLACK;
	       can->position = i;
	    }else if ((!found_violated) && (slack>lpetol) &&
		      !(cut->sense == 'R' || cut->sense == 'E')){
	       pwl[sc_cand_num].lhs = left_hand_side;
	       pwl[sc_cand_num++].position = i;
	    }
	    break;
	    
	  case SUBTOUR_ELIM_ACROSS:
	    slack = left_hand_side  - cut->rhs;
	    if (slack < -lpetol){/*---------------- This cut became violated */
	       found_violated = TRUE;
	       can = cand_list[candnum++] =
		  (branch_obj *) calloc(1, sizeof(branch_obj));
	       can->type = VIOLATED_SLACK;
	       can->position = i;
	    }else if ((!found_violated) && (slack>lpetol) &&
		      !(cut->sense=='R' || cut->sense=='E')){
	       pwl[sc_cand_num].lhs = left_hand_side;
	       pwl[sc_cand_num++].position = i;
	    }
	    break;
	    
	 default:
	    break;
	 }
      }
      if (found_violated){
	 *cand_num = candnum;
	 FREE(pwl);
	 *action = USER__DO_NOT_BRANCH;
	 return(USER_NO_PP);
      }
      
      /* now go through the slack rows in the matrix and add the potential
	 candidates to the end of the pos/weight lists */
      for (i = 0; i < slacks_in_matrix_num; i++)
	 if ((cut = slacks_in_matrix[i]) &&
	     (cut->type==SUBTOUR_ELIM_SIDE || cut->type== SUBTOUR_ELIM_ACROSS)
	     && !(cut->sense=='R' || cut->sense=='E')){
	    left_hand_side = compute_lhs(varnum, userind, x, cut,
					 vrp->vertnum);
	    switch (cut->type){
	     case SUBTOUR_ELIM_SIDE:
	       slack = cut->rhs - left_hand_side;
	       /* if (slack > lpetol && slack < lpetol1){*/
	       if (slack > lpetol ){
		  pwl[sc_cand_num+sim_cand_num].lhs = left_hand_side;
		  pwl[sc_cand_num+sim_cand_num].position = i;
		  sim_cand_num++;    
	       }
	       break;
	       
	     case SUBTOUR_ELIM_ACROSS:
	       slack = left_hand_side - cut->rhs;
	       /* if (slack > lpetol && slack < 2-lpetol){*/
	       if (slack > lpetol ){
		  pwl[sc_cand_num+sim_cand_num].lhs = left_hand_side;
		  pwl[sc_cand_num+sim_cand_num].position = i;
		  sim_cand_num++;
	       }
	       break;
	       
	     default:
	       break;
	    }
	 }
      
      /* set the children's rhs etc */
      for (i = 0 ; i < sc_cand_num + sim_cand_num; i++){
	 can = cand_list[candnum++] =
	    (branch_obj *) calloc(1, sizeof(branch_obj));
	 if (i < sc_cand_num ){
	    can->type = CANDIDATE_CUT_NOT_IN_MATRIX;
	    cut = slack_cuts[can->position = pwl[i].position ];
	    user_unpack_cuts(user, CUT_NOT_IN_MATRIX_SLACK,
			     UNPACK_CUTS_SINGLE, varnum, vars, 1,
			     slack_cuts+can->position, &new_row_num,
			     &new_rows);
	    can->row = *new_rows;
	    cut = can->row->cut;
	    FREE(new_rows);
	 }else{
	    can->type = CANDIDATE_CUT_IN_MATRIX;
	    cut = slacks_in_matrix[can->position = pwl[i].position  ];
	 }
	 can->child_num = 2;
	 can->lhs = pwl[i].lhs;
	 /* no need to allocate these. they are of fixed length
	    can->sense = (char *) malloc(2 * CSIZE);
	    can->rhs = (double *) malloc(2 * DSIZE);
	    can->range = (double *) malloc(2 * DSIZE);
	    */
	 switch (cut->type){
	  case SUBTOUR_ELIM_SIDE:
	    if ((slack = cut->rhs - can->lhs) < lpetol1){ 
	       can->sense[0] = 'E';
	       can->rhs[0] = cut->rhs;
	       can->range[0] = 0;
	       can->branch[0] = DO_NOT_BRANCH_ON_THIS_ROW;
	       can->sense[1] = 'L';
	       can->rhs[1] = cut->rhs - 1;
	       can->range[1] = 0;
	       can->branch[1] = ALLOWED_TO_BRANCH_ON;
	    }else{
	       can->sense[0] = 'R';
	       can->rhs[0] = cut->rhs;
	       can->range[0] = -floor(slack) ;
	       can->branch[0] = DO_NOT_BRANCH_ON_THIS_ROW;
	       can->sense[1] = 'L';
	       can->rhs[1] = cut->rhs - ceil(slack);
	       can->range[1] = 0;
	       can->branch[1] = ALLOWED_TO_BRANCH_ON;
	    }
	    break;
	  case SUBTOUR_ELIM_ACROSS:
	    if ((slack = can->lhs - cut->rhs) < 2-lpetol){
	       can->sense[0] = 'E';
	       can->rhs[0] = cut->rhs;
	       can->range[0] = 0;
	       can->branch[0] = DO_NOT_BRANCH_ON_THIS_ROW;
	       can->sense[1] = 'G';
	       can->rhs[1] = cut->rhs + 2;
	       can->range[1] = 0;
	       can->branch[1] = ALLOWED_TO_BRANCH_ON;
	    }else{
	       can->sense[0] = 'R';
	       can->rhs[0] = cut->rhs;
	       can->range[0] = 2*floor(slack/2);
	       can->branch[0] = DO_NOT_BRANCH_ON_THIS_ROW;
	       can->sense[1] = 'G';
	       can->rhs[1] = cut->rhs + 2*floor(slack/2) +2;
	       can->range[1] = 0;
	       can->branch[1] = ALLOWED_TO_BRANCH_ON;
	    }
	    break;
	 /*__BEGIN_EXPERIMENTAL_SECTION__*/
	 case FARKAS:
	    break;
	 case NO_COLUMNS:
	    break;
	 /*___END_EXPERIMENTAL_SECTION___*/
	 }
      }
      FREE(pwl);
      FREE(userind);
      *cand_num = candnum;
      cand_list = *candidates + *cand_num;

    case FALSE:

      switch (((vrp_spec *)user)->par.branching_rule){
       case 0:
	 {
	    int *xind = (int *) malloc(varnum*ISIZE);
	    double *xval = (double *) malloc(varnum*DSIZE);
	    int cnt = 0;
	    
	    for (i = varnum-1; i >= 0; i--){
	       fracx = x[i] - floor(x[i]);
	       if (fracx > lpetol && fracx < 1-lpetol){
		  xind[cnt] = i;
		  xval[cnt++] = fabs(fracx - .5);
	       }
	    }
	    qsortucb_di(xval, xind, cnt);

	    candnum = vrp->par.strong_branching_cand_num_max -
	       vrp->par.strong_branching_red_ratio * bc_level;
	    candnum = MAX(candnum, vrp->par.strong_branching_cand_num_min);
	    candnum = MIN(candnum, cnt);
	    
	    for (i = candnum-1; i >= 0; i--){
	       can=cand_list[i]=(branch_obj *) calloc(1, sizeof(branch_obj));
	       can->type = CANDIDATE_VARIABLE;
	       can->child_num = 2;
	       can->position = xind[i];
	       can->sense[0] = 'L';
	       can->sense[1] = 'G';
	       can->rhs[0] = floor(x[xind[i]]);
	       can->rhs[1] = can->rhs[0] + 1;
	       can->range[0] = can->range[1] = 0;
	    }
	    FREE(xind);
	    FREE(xval);
	 }
       break;
       
       case 1:
	 candnum = 0;
	 for (i = varnum-1; i >= 0; i--){
	    fracx = x[i] - floor(x[i]);
	    if (fracx > lpetol && fracx < lpetol5){
	       can = cand_list[candnum++] =
		  (branch_obj *) calloc(1, sizeof(branch_obj));
	       can->type = CANDIDATE_VARIABLE;
	       can->child_num = 2;
	       can->position = i;
	       can->sense[0] = 'L';
	       can->sense[1] = 'G';
	       can->rhs[0] = floor(x[i]);
	       can->rhs[1] = can->rhs[0] + 1;
	       can->range[0] = can->range[1] = 0;
	    }
	 }
	 break;
	 
      case 2:
	candnum = vrp->par.strong_branching_cand_num_max -
	   vrp->par.strong_branching_red_ratio * bc_level;
	candnum = MAX(candnum, vrp->par.strong_branching_cand_num_min);

	branch_close_to_half(candnum, &candnum, &cand_list);
	break;
      }
      *cand_num += candnum;
      *action = USER__DO_BRANCH;
   }
   return(USER_NO_PP);
}

/*===========================================================================*/
      
double compute_lhs(int number, int *indices, double *values, cut_data *cut,
		   int vertnum)
{
   char *coef;
   int v0, v1;
   double lhs = 0;
   int i;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   int num_arcs, edge_index ; 
   char *cpt; 
   int *arcs ;
   char *indicators;
   double bigM, *weights ;
   int jj;
   /*___END_EXPERIMENTAL_SECTION___*/
 
   switch (cut->type){
    
    case SUBTOUR_ELIM_SIDE:
      coef = (char *)(cut->coef);
      for (i = 0, lhs = 0; i<number; i++){
	 BOTH_ENDS(indices[i], &v1, &v0);
	 if (coef[v0 >> DELETE_POWER] & (1 << (v0 & DELETE_AND)) &&
	     (coef[v1 >> DELETE_POWER]) & (1 << (v1 & DELETE_AND)))
	    lhs += values[i];
	   
      }
      return(lhs);

    case SUBTOUR_ELIM_ACROSS:
      coef = (char *)(cut->coef);
      for (lhs = 0, i = 0; i<number; i++){
	 BOTH_ENDS(indices[i], &v1, &v0);
	 if ((coef[v0 >> DELETE_POWER] >> (v0 & DELETE_AND) & 1) ^
	     (coef[v1 >> DELETE_POWER] >> (v1 & DELETE_AND) & 1))
	    lhs += values[i];
      }
      
      return(lhs);
    /*__BEGIN_EXPERIMENTAL_SECTION__*/

    case FARKAS:
      coef = (char *)(cut->coef);
      arcs = (int *) calloc( ((vertnum*vertnum)/2), sizeof(int));
      weights = (double *) calloc( ((vertnum*vertnum)/2), sizeof(double)) ;
      
      cpt=coef+ ((vertnum >> DELETE_POWER) + 1); 
      memcpy((char *)&num_arcs, cpt, ISIZE);
      memcpy((char *)arcs, cpt, ISIZE*(num_arcs +1));
      arcs++;
      cpt+=(ISIZE*(num_arcs +1))/CSIZE;
      memcpy((char *)weights, cpt, DSIZE*(num_arcs +1));
      bigM=(*(double *)weights);
      weights++;
      
      for (i = 0, lhs=0 ; i<number; i++){
	 BOTH_ENDS(indices[i], &v1, &v0);
	 edge_index=indices[i];
	    
	 if (isset(coef, v1) || isset(coef,v0)){
	    for (jj=0; jj<num_arcs; jj++){
	       if ( arcs[jj]==edge_index){
		       lhs+=weights[jj]*values[i];
		       break;
	       }
	    }
	    if (jj==num_arcs) lhs+=bigM*values[i];
	    
	 }
      }
      arcs--;
      weights--;
      FREE(arcs);
      FREE(weights);
      return(lhs); 
     
      
    case NO_COLUMNS:
      coef = (char *)(cut->coef);
      arcs = (int *) malloc( ((vertnum*vertnum)/2)* sizeof(int));
      indicators = (char *) malloc( ((vertnum*vertnum)/2)* sizeof(char));
	
      cpt=coef+ ((vertnum >> DELETE_POWER) + 1); 
      memcpy((char *)&num_arcs, cpt, ISIZE);
      cpt+= ISIZE /CSIZE;
      memcpy((char *)arcs, cpt, ISIZE*(num_arcs));
      cpt+=(ISIZE*(num_arcs ))/CSIZE;
      memcpy((char *)indicators, cpt, CSIZE*(num_arcs));
	
      for (i = 0, lhs=0 ; i<number; i++){
	 BOTH_ENDS(indices[i], &v1, &v0);
	 edge_index=indices[i];
	 	 
	 if (isset(coef, v1) || isset(coef,v0)){
	    for (jj=0; jj<num_arcs; jj++){
	       if ( arcs[jj]==edge_index){
		  lhs+=values[i]* ( indicators[jj] ? 1.0:0.0);
		  break;
	       }
	    }
	    if (jj==num_arcs) lhs-=values[i];
	 }
      }
      FREE(arcs);
      FREE(indicators);
      return(lhs); 
    /*___END_EXPERIMENTAL_SECTION___*/
      
    default:
      printf("Cut type's not recognized! \n\n");
      return(0);
   }
}

/*__BEGIN_EXPERIMENTAL_SECTION__*/
#if 0
/*===========================================================================*/

int user_select_candidates(lp_prob *p, void *user, double lpetol, int cutnum,
			   int slacks_in_matrix_num,
			   cut_data **slacks_im_matrix, int slack_cut_num,
			   cut_data **slack_cuts, int varnum, var_desc **vars,
			   double *x, char *status, int *cand_num,
			   branch_obj ***candidates, int *action)

{
   return(USER__CLOSE_TO_HALF);
}

/*===========================================================================*/

int user_select_candidates(lp_prob *p, void *user, double lpetol, int cutnum,
			   int slacks_in_matrix_num,
			   cut_data **slacks_im_matrix, int slack_cut_num,
			   cut_data **slack_cuts, int varnum, var_desc **vars,
			   double *x, char *status, int *cand_num,
			   branch_obj ***candidates, int *action)

{
   vrp_spec *vrp = (vrp_spec *)user;
   int i;
   int *sorted_cand_list, j = 0, k = 0;
   double *sorted_diffs;
   int cand_num = vrp->par.max_branch_cand;
   int branch_rule = vrp->par.branching_rule;
   branch_obj **candidates, *cand;

   candidates =
      *candidates = (branch_obj **) malloc(cand_num * sizeof(branch_obj *));

   sorted_cand_list = (int *) calloc (varnum, sizeof(int));
   sorted_diffs     = (double *) calloc (varnum, sizeof(double));
   for (i = 0; i < varnum; i++){
      sorted_cand_list[i] = i;
      sorted_diffs[i] = fabs(x[i] - .5);
   }

   cand_sort(sorted_cand_list, sorted_diffs, varnum);

   switch(vrp->par.branching_rule){ /* which branching rule to use */
      
   case DEPOTS_AT_HALF_BRANCH_RIGHT:
   case DEPOTS_AT_HALF:
   case VARS_AT_HALF_PREFER_DEPOT_BRANCH_RIGHT:
   case VARS_AT_HALF_PREFER_DEPOT:
      
      for (i=0, j=0; sorted_diffs[i]<lpetol && i<varnum && j<cand_num; i++){
	 /* Here we find the depot edges that are at a half*/
	 if (vrp->edges[(vars[sorted_cand_list[i]]->userind) << 1] == 0){
	    /* This userind serves the purpose of checking whether the
	       variable corresponds to an edge connected to the depot */
	    cand = candidates[i] = (branch_obj *)calloc(1, sizeof(branch_obj));
	    cand->position = sorted_cand_list[i];
	    cand->type = CANDIDATE_VARIABLE;
	    cand->child_num = 2;
	    cand->sense[0] = 'L';
	    cand->sense[1] = 'G';
	    cand->rhs[0] = 0;
	    cand->rhs[1] = 1;
	    cand->range[0] = cand->range[1] = 0;
	 }
      }
      if (j && (branch_rule == DEPOTS_AT_HALF_BRANCH_RIGHT ||
		branch_rule == DEPOTS_AT_HALF)){
	 *cand_num = j;
	 break;
      }
      for (k = 0; k < i && j < cand_num; k++){
	 if (vrp->edges[(vars[sorted_cand_list[k]]->userind) << 1] != 0)
	    cand = candidates[j++] = (branch_obj *)
	                               calloc(1, sizeof(branch_obj));
	    cand->position = sorted_cand_list[i];
	    cand->type = CANDIDATE_VARIABLE;
	    cand->child_num = 2;
	    cand->sense[0] = 'L';
	    cand->sense[1] = 'G';
	    cand->rhs[0] = 0;
	    cand->rhs[1] = 1;
	    cand->range[0] = cand->range[1] = 0;
      }
      if (j){
	 *cand_num = j;
	 break;
      }
      /*otherwise, fall through*/

    case DEPOTS_CLOSEST_TO_HALF:
    case DEPOTS_CLOSEST_TO_HALF_BRANCH_RIGHT:
    case VARS_CLOSEST_TO_HALF_PREFER_DEPOT:
    case VARS_CLOSEST_TO_HALF_PREFER_DEPOT_BRANCH_RIGHT:

      for (i = 0, j = 0;
	   i < varnum && sorted_diffs[i]<.5-lpetol && j < cand_num; i++){
	 if (vrp->edges[(vars[sorted_cand_list[i]]->userind) << 1] == 0){
	    cand = candidates[j++] = (branch_obj *)
	       calloc(1, sizeof(branch_obj));
	    cand->position = sorted_cand_list[i];
	    cand->type = CANDIDATE_VARIABLE;
	    cand->child_num = 2;
	    cand->sense[0] = 'L';
	    cand->sense[1] = 'G';
	    cand->rhs[0] = 0;
	    cand->rhs[1] = 1;
	    cand->range[0] = cand->range[1] = 0;
	 }
      }
      if (j && (branch_rule == DEPOTS_CLOSEST_TO_HALF ||
		branch_rule == DEPOTS_CLOSEST_TO_HALF_BRANCH_RIGHT)){
	 *cand_num = j;
	 break;
      }
      for (k = 0; k < i && j < cand_num; k++){
	 cand = candidates[j++] = (branch_obj *)calloc(1, sizeof(branch_obj));
	 cand->position = sorted_cand_list[i];
	 cand->type = CANDIDATE_VARIABLE;
	 cand->child_num = 2;
	 cand->sense[0] = 'L';
	 cand->sense[1] = 'G';
	 cand->rhs[0] = 0;
	 cand->rhs[1] = 1;
	 cand->range[0] = cand->range[1] = 0;
      }

      if (j){
	 *cand_num = j;
	 break;
      }
      /*otherwise, fall through*/
      
    case VARS_CLOSEST_TO_HALF:
      
      for (i = 0, j = 0;
	   sorted_diffs[i]<.5-lpetol && i < varnum && j < cand_num; i++){
	 cand = candidates[j++] = (branch_obj *)calloc(1, sizeof(branch_obj));
	 cand->position = sorted_cand_list[i];
	 cand->type = CANDIDATE_VARIABLE;
	 cand->child_num = 2;
	 cand->sense[0] = 'L';
	 cand->sense[1] = 'G';
	 cand->rhs[0] = 0;
	 cand->rhs[1] = 1;
	 cand->range[0] = cand->range[1] = 0;
      }
      if (j){
	 *cand_num = j;
	 break;
      }
      /*otherwise, fall through*/

    case BEST_K:
#if 0
      for (i = 0, j = 0; i < varnum && j < cand_num; i++){
	 if (p->lp_data->fixed[sorted_cand_list[i]] == NOT_FIXED){
	    cand = candidates[j++] = (branch_obj *)
	                                 calloc(1, sizeof(branch_obj));
	    cand->position = sorted_cand_list[i];
	    cand->type = CANDIDATE_VARIABLE;
	    cand->child_num = 2;
	    cand->sense[0] = 'L';
	    cand->sense[1] = 'G';
	    cand->rhs[0] = 0;
	    cand->rhs[1] = 1;
	    cand->range[0] = cand->range[1] = 0;
	 }
      }

      *cand_num = j;
      break;
#endif
      *cand_num = 0;
      break;

    case VARS_AT_HALF:

      for (j = 0; sorted_diffs[j] < lpetol && j < cand_num; j++){
	 cand = candidates[j++] = (branch_obj *)calloc(1, sizeof(branch_obj));
	 cand->position = sorted_cand_list[i];
	 cand->type = CANDIDATE_VARIABLE;
	 cand->child_num = 2;
	 cand->sense[0] = 'L';
	 cand->sense[1] = 'G';
	 cand->rhs[0] = 0;
	 cand->rhs[1] = 1;
	 cand->range[0] = cand->range[1] = 0;
      }

      *cand_num = j;
      break;

    default:
      FREE(sorted_cand_list);
      FREE(sorted_diffs);
      return(DEFAULT);
   }
   FREE(sorted_cand_list);
   FREE(sorted_diffs);

   return(USER_NO_PP);
}
#endif

/*___END_EXPERIMENTAL_SECTION___*/
/*===========================================================================*/

/*===========================================================================*\
 * I wrote my own function to compare candidates. Maybe this should go
 * into SYMPHONY.
\*===========================================================================*/

int user_compare_candidates(void *user, branch_obj *can1, branch_obj *can2,
			    double ub, double granularity,
			    int *which_is_better)
{
   int i, j;
   double low1, low2;
   
   for (i = 1; i >= 0; i--)
      if (can1->termcode[i] == OPT_FEASIBLE || can1->termcode[i] == D_OBJLIM ||
	  can1->termcode[i] == D_UNBOUNDED ||
	  (can1->termcode[i] == OPTIMAL && can1->objval[i] > ub - granularity))
	 break;

   for (j = 1; j >= 0; j--)
      if (can2->termcode[j] == OPT_FEASIBLE || can2->termcode[j] == D_OBJLIM ||
	  can2->termcode[j] == D_UNBOUNDED ||
	  (can2->termcode[j] == OPTIMAL && can2->objval[j] > ub - granularity))
	 break;

   if (i < 0 && j < 0)
      return(DEFAULT);
      
   if (i < 0 && j > 0){
      *which_is_better = SECOND_CANDIDATE_BETTER;
      return(USER_NO_PP);
   }

   if (i > 0 && j < 0){
      *which_is_better = FIRST_CANDIDATE_BETTER;
      return(USER_NO_PP);
   }

   low1 = i ? can1->objval[0] : can1->objval[1];
   low2 = j ? can2->objval[0] : can2->objval[1];

   *which_is_better = low1 > low2 ? FIRST_CANDIDATE_BETTER :
                                    SECOND_CANDIDATE_BETTER;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * You can let SYMPHONY choose which child to retain. The default is
 * to keep the one with the lower objective function value.
\*===========================================================================*/

int user_select_child(void *user, double ub, branch_obj *can, char *action)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, you can print out a more identifiable description of the
 * branching object than just "variable 51". I print out the end points
 * of the edge if an edge is branched on.
\*===========================================================================*/

int user_print_branch_stat(void *user, branch_obj *can, cut_data *cut,
			   int n, var_desc **vars, char *action)
{
   int v0, v1, i;
   char *coef;
   
   if (cut){
      switch(cut->type){
       case SUBTOUR_ELIM_SIDE:
	 coef = (char *)(cut->coef);
	 for (i = 0; i < n; i++){
	    BOTH_ENDS(vars[i]->userind, &v1, &v0);
	    if (coef[v0 >> DELETE_POWER] & (1 << (v0 & DELETE_AND)) &&
		(coef[v1 >> DELETE_POWER]) & (1 << (v1 & DELETE_AND)))
	       printf("Edge (%i, %i)\n", v0, v1);
	 }

       case SUBTOUR_ELIM_ACROSS:
         coef = (char *)(cut->coef);
	 for (i = 0; i < n; i++){
	    BOTH_ENDS(vars[i]->userind, &v1, &v0);
	    if ((coef[v0 >> DELETE_POWER] >> (v0 & DELETE_AND) & 1) ^
		(coef[v1 >> DELETE_POWER] >> (v1 & DELETE_AND) & 1))
	       printf("Edge (%i, %i)\n", v0, v1);
	 } 
      }
   }else{
      BOTH_ENDS(vars[can->position]->userind, &v1, &v0);
      printf("Edge (%i, %i)\n", v0, v1);
   }

   return(USER_NO_PP);
}
