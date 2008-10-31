/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*===========================================================================*/
/* This file is a part of the symphony-processor                             */
/*===========================================================================*/

#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sym_preprocessor.h"
#include "sym_prep_params.h"
#include "sym_master.h" 
#include "sym_constants.h" 
#include "sym_qsort.h"

/*===========================================================================*/
/*===========================================================================*/
int prep_basic(PREPdesc *P)
{

   /*
     This function is the master of the preprocessing part. It calls and 
     controls other functions to perform preprocessing jobs.
   */

   int termcode;     /* return status of each function called
		       herein */
   int iter_cnt = 0, iter_cnt_limit, p_level;    
   int verbosity;
   double a_val, etol, new_bound;// min_ub, max_lb; 
   int do_sr_rlx;
   char do_aggr_row_rlx;// fix_var; 
   int i, j, m, n, nz, *r_matbeg, *r_matind, *matbeg, *matind; 
   //int i, max_size, min_size, *max_ind, *min_ind; 
   double *obj, *rhs, *r_matval, *matval, *ub, *lb;
   char *sense; 
   //int * updated_cols_ind, updated_cols_cnt;
   //int * updated_rows_ind, updated_rows_cnt;
   int col_ind, row_ind, fix_type;
   //int * modified_cols_ind, * modified_rows_ind;
   //   int modified_cols_cnt, modified_rows_cnt;
   //char *is_col_updated, *is_row_updated;
   char can_impl = FALSE, bin_type = FALSE;
   int dive_level, impl_dive_level, impl_limit;
   int old_changes_cnt, new_changes_cnt, init_changes_cnt = 0;
   /* first initialize P and mip etc */
   int old_others_cnt, new_others_cnt, mark_others_cnt = 0;
   //double start_time = wall_clock(NULL);
   double mark_time, impl_time = 0.0;
   prep_stats *stats = &(P->stats);
   prep_params params = P->params;
   
   stats->row_infeas_ind = stats->col_infeas_ind = -1;
   
   verbosity = params.verbosity; /* debug */
   p_level = params.level;
   etol = params.etol;
   iter_cnt_limit = params.iteration_limit;
   do_sr_rlx = params.do_single_row_rlx;      
   do_aggr_row_rlx = params.do_aggregate_row_rlx;
   dive_level = params.dive_level;
   impl_dive_level = params.impl_dive_level;
   impl_limit = params.impl_limit;

   MIPdesc *mip = P->mip;
   MIPinfo *mip_inf = mip->mip_inf;
   COLinfo * cols = mip_inf->cols;
   ROWinfo *rows = mip_inf->rows;

   if(mip_inf->prob_type == CONTINUOUS_TYPE){
      /* no need for prep, just quit */
      return PREP_UNMODIFIED;
   }
   
   /* first integerize the bounds */
   /* can be embedded somewhere in basic prep down*/
   /* for now let it be */
   /* also we integerize the integerizable vars if p_level > 2 */
   termcode = prep_integerize_bounds(P);
   
   if(prep_quit(termcode)){
      return termcode;
   }
   
   /* Display the whole mip - for debugging */
#ifdef DISPLAY_MIP
   termstatus = prep_display_mip(P->mip);
#endif   
   
   m = mip->m;
   n = mip->n;
   nz = mip->nz;

   r_matbeg = mip->row_matbeg;
   r_matind = mip->row_matind;
   r_matval = mip->row_matval;

   matbeg = mip->matbeg;
   matind = mip->matind;
   matval = mip->matval;

   ub = mip->ub;
   lb = mip->lb;

   obj = mip->obj;
   sense = mip->sense;
   rhs = mip->rhs;


   char need_reset, *impl_vars = NULL;//, *impl_vars_checked = NULL;
   //int var_ind, impl_cnt, impl_cnt_limit;
   //int *ind_list = NULL, *impl_vars_weight = NULL;
   
#if 0
   updated_cols_ind = (int *) malloc(ISIZE*n);
   updated_rows_ind = (int *) malloc(ISIZE*m);
   modified_cols_ind = (int *) malloc(ISIZE*n);
   modified_rows_ind = (int *) malloc(ISIZE*m);

   is_col_updated = (char *) calloc(CSIZE,n);
   is_row_updated = (char *) calloc(CSIZE,m);
   
   
   updated_rows_cnt = m; 
   updated_cols_cnt = n;
   

   if(m <= n){
      max_size = n;
      min_size = m;
      max_ind = updated_cols_ind;
      min_ind = updated_rows_ind;
   }else{
      max_size = m;
      min_size = n;
      max_ind = updated_rows_ind;
      min_ind = updated_cols_ind;
   }

   for(i = 0; i < min_size; i++){
      max_ind[i] = i;
      min_ind[i] = i;
   }
   for(; i < max_size; i++){
      max_ind[i] = i;
   }
#endif 

   if(mip_inf->prob_type == BINARY_TYPE ||
      mip_inf->prob_type == BIN_CONT_TYPE ||
      mip_inf->prob_type == BIN_INT_TYPE ||
      mip_inf->prob_type == ALL_MIXED_TYPE){
      
      bin_type = TRUE;
   } 
   
   /* first check duplicate rows, cols */
   
   termcode = prep_delete_duplicate_rows_cols(P, TRUE, TRUE);
      
   if(prep_quit(termcode)){
      return termcode;
   } 

   init_changes_cnt = stats->vars_fixed + stats->rows_deleted;
   
   if(p_level >= 5 && can_impl){ /* disabled now */ 
      if(bin_type){	 
	 /* for now, just between binary variables */
	 P->impl_rows = (ROWinfo *)malloc(sizeof(ROWinfo)*m); 
	 P->impl_cols = (COLinfo *)malloc(sizeof(COLinfo)*n); 
	 P->impl_ub = (double *) malloc(DSIZE*n);
	 P->impl_lb = (double *) malloc(DSIZE*n);
	 P->impl_var_ind = (int *)malloc(ISIZE * n);
	 P->impl_var_stat = (char *)malloc(CSIZE * n);
	 
	 P->ulist_checked = (char *)malloc(CSIZE * n);
	 P->llist_checked = (char *)malloc(CSIZE * n);
	 
	 P->impl_limit = impl_limit;
	 can_impl = TRUE;

	 /* get the list of columns to apply impl on */

	 
	 /* not effective */
	 impl_vars = (char *) calloc(CSIZE,n);
#if 0
	 impl_vars_weight = (int *) calloc(ISIZE,n);
	 impl_vars_checked = (char *) malloc(CSIZE*n);
	 ind_list = (int *) malloc(ISIZE*n);

	 for(col_ind = 0; col_ind < n; col_ind++){
	    ind_list[col_ind] = col_ind;
	    if(cols[col_ind].var_type != 'B'){
	       impl_vars_weight[col_ind] = 0;
	       continue;
	    }
	    memset(impl_vars_checked, FALSE, CSIZE*n);
	    for(i = matbeg[col_ind]; i < matbeg[col_ind + 1]; i++){
	       row_ind = matind[i];
	       for(j = r_matbeg[row_ind]; j < r_matbeg[row_ind + 1]; j++){
		  var_ind = r_matind[j];
		  if(var_ind != col_ind && 
		     cols[var_ind].var_type == 'B' && 
		     !impl_vars_checked[var_ind]){
		     impl_vars_checked[var_ind] = TRUE;
		     impl_vars_weight[col_ind]++;
		  }
	       }
	    }
	 }
	  
	 qsort_ii(impl_vars_weight, ind_list, n);
	 /* fixme this ugly thing here -very important, choose the 
	    vars varefully */
	 if(p_level >= 6 || mip_inf->binary_var_num < 2000){
	    impl_cnt_limit = n;
	 }else{
	    impl_cnt_limit = 1000 + (int)(mip_inf->binary_var_num/2);
	    if(impl_cnt_limit > mip_inf->binary_var_num){
	       impl_cnt_limit = n;
	    }
	 }
	 impl_cnt_limit = n;
	 impl_cnt = 0;
	 //for(i = n - 1; i >= 0 && impl_cnt < impl_cnt_limit; i--){
	 //  if(impl_vars_weight[i] == 0){
#if 0
	 for(i = 0; i < n && impl_cnt < impl_cnt_limit; i++){
	    if(impl_vars_weight[i] == nz + 1){
	       break;
	    }
	    impl_vars[ind_list[i]] = TRUE;
	    impl_cnt++;
	 }
#endif
	 for(i = n - 1; i >= 0 && impl_cnt < impl_cnt_limit; i--){
	    if(impl_vars_weight[i] == 0){
	       break;
	    }
	    impl_vars[ind_list[i]] = TRUE;
	    impl_cnt++;
	 }
#endif

	 for(i = n - 1; i >= 0; i--){
	    impl_vars[i] = TRUE;

	 }
#if 0
	 FREE(impl_vars_weight);
	 FREE(ind_list);
	 FREE(impl_vars_checked);
#endif
      }
   }
   
   /* main preprocessing loop */

   old_changes_cnt = new_changes_cnt = 0;
   old_others_cnt = new_others_cnt = 0;
   while(iter_cnt < iter_cnt_limit){// && updated_rows_cnt > 0){
      
      iter_cnt++;
      //modified_rows_ind = 0;
      //modified_cols_ind = 0;

      PRINT(verbosity, 1, ("Basic iteration number: %d\n", iter_cnt));

      /* check the updated bounds and cols to iterate on*/
      /* for each updated row and column do bla bla...*/
      /* while iterating update the new cols and rows indices*/
      
      /* do this only once at the end of first iteration*/
      //if(prep_level > 2 && (do_sr_rlx || do_aggr_row_rlx) && 
      // iter_cnt == 1){
      /*=====================================================================*/
      /*=====================================================================*/
      
      for(col_ind = 0; col_ind < n; col_ind++){
	 //printf("working on var %i\n", col_ind);	 
	 if(cols[col_ind].var_type == 'F'){
	    continue;
	 }
	 /* can we fix it? first check implications */	 
	 if(can_impl && iter_cnt < 2 && impl_time < params.time_limit){
	    //printf("current impl time %f\n", impl_time); 
	    if(cols[col_ind].var_type == 'B' && 
	       impl_vars[col_ind]){// && col_ind < n/10){
	       //printf("impl_process var %i\n", col_ind);
	       
	       /* fist copy initial info 
		  fixme! - work on this! 
	       */
	       
	       /* do once for each variable */
	       mark_time = wall_clock(NULL);
	       fix_type = FIX_NO_BOUND;

	       need_reset = FALSE;
	       mark_time = wall_clock(NULL);
	       memcpy(P->impl_rows, rows, sizeof(ROWinfo)*m); 
	       memcpy(P->impl_cols, cols, sizeof(COLinfo)*n); 
	       memcpy(P->impl_ub, ub, DSIZE*n);
	       memcpy(P->impl_lb, lb, DSIZE*n);
	       P->impl_stats = P->stats;
	       P->alloc_time += wall_clock(NULL) - mark_time;

	       if(cols[col_ind].sign_type != ALL_NEG_VEC){
		  need_reset = TRUE;
		  if(!cols[col_ind].ulist){
		     cols[col_ind].ulist = 
			(IMPlist *)calloc(sizeof(IMPlist),1);
		  }	       
		  
		  P->list = cols[col_ind].ulist;	      
		  /* fix it to 1.0 and see if that causes any infeasibility 
		     otherwise get the impllist and continue*/
		  /* get the implication list */
		  
		  //	       termcode = prep_get_impl_list(P, col_ind);
		  /* fix this column, update row bounds of this column
		     check for redundancy, */	       	       
		  P->impl_col_ind = col_ind;
		  //if(col_ind == 6989){
		  //  P->params.verbosity = 3;
		  //  impl_dive_level = 100;
		  //}
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind, -1,
							    impl_dive_level,
							    1.0,
							    FIX_BINARY, TRUE, 
							    TRUE);
		  //prep_delete_imp_list(&(cols[col_ind].ulist));
		  if(termcode == PREP_INFEAS){
		     //printf("infeasibility detected! %i %f \n", col_ind, 1.0);
		     prep_delete_imp_list(&(cols[col_ind].ulist));
		     /*then this column is fixable to its lower bound! */
		     new_bound = 0.0;
		     fix_type = FIX_BINARY;
		  }
	       }
	       
	       if(fix_type != FIX_BINARY && 
		  cols[col_ind].sign_type != ALL_POS_VEC){
		  mark_time = wall_clock(NULL);
		  /* reset what we had */
		  /* need to reset? */
		  if(need_reset){
		     memcpy(rows, P->impl_rows,sizeof(ROWinfo)*mip->m); 
		     memcpy(cols, P->impl_cols, sizeof(COLinfo)*mip->n); 
		     
		     memcpy(ub, P->impl_ub, DSIZE*mip->n);
		     memcpy(lb, P->impl_lb, DSIZE*mip->n);
		     P->stats = P->impl_stats;
		     P->alloc_time += wall_clock(NULL) - mark_time;
		  }
		  if(!cols[col_ind].llist){
		     cols[col_ind].llist = 
			(IMPlist *)calloc(sizeof(IMPlist),1);
		  }
		  
		  P->list = cols[col_ind].llist;	      
		  P->impl_col_ind = col_ind;

		  termcode = prep_modified_cols_update_info(P, 1, &col_ind, -1,
							    impl_dive_level,
							    0.0, FIX_BINARY, 
							    TRUE, TRUE);
		  if(termcode == PREP_INFEAS){
		     //printf("infeasibility detected! %i %f\n", col_ind, 0.0);
		     /*then this column is fixable to its lower bound! */
		     prep_delete_imp_list(&(cols[col_ind].llist));
		     new_bound = 1.0;
		     fix_type = FIX_BINARY;
		  }
	       }
	       
	       /* now get back */
	       mark_time = wall_clock(NULL);
	       memcpy(rows, P->impl_rows,sizeof(ROWinfo)*mip->m); 
	       memcpy(cols, P->impl_cols, sizeof(COLinfo)*mip->n); 	    
	       memcpy(ub, P->impl_ub, DSIZE*mip->n);
	       memcpy(lb, P->impl_lb, DSIZE*mip->n);
	       P->stats = P->impl_stats;
	       P->alloc_time += wall_clock(NULL) - mark_time;    
	       /* and now check if we can fix anything */
	       //prep_delete_imp_list(&(cols[col_ind].llist));	       
	       
	       impl_time += wall_clock(NULL) - mark_time;
	       if(fix_type == FIX_BINARY){
		  termcode = prep_modified_cols_update_info(P, 1, 
							    &col_ind, -1, 
							    dive_level,
							    new_bound, 
							    fix_type, TRUE, 
							    FALSE);
		  if(prep_quit(termcode)){
		     return termcode;
		  }
		  continue;
	       }
	       //printf("\tdone %i\n", col_ind);
	    }
	 }

	 /* couldnt fix it, continue */

	 for(j = matbeg[col_ind]; j < matbeg[col_ind+1]; j++){
	    row_ind = matind[j];
	    
	    if(rows[row_ind].is_redundant){
	       continue;
	    }else{
	       if(rows[row_ind].ub >= INF && rows[row_ind].lb <= -INF){
		  continue;
	       }

	       termcode = prep_check_redundancy(P, row_ind, FALSE, 0.0, 0.0,
						FALSE, 0);
	       if(prep_quit(termcode)){
		  return termcode;
	       }
	       if(rows[row_ind].is_redundant){
		  continue;
	       }
	    }
	    
	    a_val = matval[j];

	    if(rows[row_ind].ub_inf_var_num <= 1 ||
	       rows[row_ind].lb_inf_var_num <= 1){
	       if(ub[col_ind] >= INF || lb[col_ind] <= -INF){
		  if((a_val > etol && ub[col_ind] >= INF) ||
		     (a_val < -etol && lb[col_ind] <= -INF) ||
		     sense[row_ind] == 'E'){
		     termcode = 
			prep_force_row_bounds(P, row_ind, col_ind, j);
		     if(prep_quit(termcode)){
			return termcode;
		     }

		     if(rows[row_ind].is_redundant){
			continue;
		     }
		  }
	       }
	    }

	    if((a_val > etol || a_val < -etol) && 
	       cols[col_ind].var_type != 'F'){
	       
	       termcode = prep_improve_variable(P, col_ind, 
						row_ind, j, 
						dive_level, TRUE, 
						FALSE, FALSE, 0.0, 
						0.0, COL_ORDERED);
	       if(prep_quit(termcode)){
		  return termcode;
	       }
	    }
	 }
      }	       
     
      new_changes_cnt = stats->rows_deleted + 
	 stats->vars_fixed - stats->vars_aggregated; 
      new_others_cnt = stats->coeffs_changed + 
	 stats->bounds_tightened;
      /* and check duplicacy */
      /* fix this, we dont want to check duplicacy for unmodified 
	 cols and rows...
      */
      if(new_changes_cnt > old_changes_cnt){
	 old_changes_cnt = new_changes_cnt;
	 old_others_cnt = new_others_cnt;
      }else{
	 if(new_others_cnt > old_others_cnt){
	    old_others_cnt = new_others_cnt;
	    mark_others_cnt++;
	 }else{
	    break;
	 }
	 if(mark_others_cnt > 3){
	    break;
	 }
      }
      //if(old_changes_cnt == new_changes_cnt){
      //	 break;
      //}else{
      //	 old_changes_cnt = new_changes_cnt;
      // }
   }
   
   //prep_report(P, termcode);
   //exit(0);
   
   if(do_sr_rlx){ 
      for(row_ind = 0; row_ind < m; row_ind++){
	 if(!rows[row_ind].is_redundant){  
	    termcode = prep_solve_sr_rlx(P, 1, &row_ind);	
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 }
      }
   }   

  /* fixme work on these */
   /* dont need to do for each col, keep an array to track them down */
   if(stats->rows_deleted > 0){
      for(col_ind = 0; col_ind < n; col_ind++){
	 if(cols[col_ind].var_type != 'F' && cols[col_ind].col_size == 0){
	    termcode = prep_improve_variable(P, col_ind, -1, 0, 
					     dive_level, TRUE, FALSE, 
					     FALSE,
					     0.0,0.0, COL_ORDERED);
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 }
      }
   }	
   
   /* similary for 0 sized rows, eliminate them */
   /* fixme, same issue as above */
   if(stats->vars_fixed > 0){
      for(row_ind = 0; row_ind < m; row_ind++){
	 if(!rows[row_ind].is_redundant &&  
	    rows[row_ind].fixed_var_num  >= rows[row_ind].size - 1){
	    termcode = prep_check_redundancy(P, row_ind, FALSE, 0.0, 0.0,
					     FALSE, dive_level);
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 } 
      }
   }
   
   if(new_changes_cnt > init_changes_cnt){
      termcode = prep_delete_duplicate_rows_cols(P, TRUE, TRUE);
   }
   if(verbosity >= 2){
      printf("total alloc time: %f\n", P->alloc_time);
      printf("total alloc time2: %f\n", P->alloc2_time);
      printf("total impl time2: %f\n", impl_time);
      printf("total impl_cols_time: %f\n", P->impl_cols_time);
      printf("total impl_rows_time: %f\n", P->impl_rows_time);
   }

   FREE(impl_vars);

   if(new_changes_cnt + stats->coeffs_changed + mip_inf->fixed_var_num > 0){
      termcode = prep_cleanup_desc(P);
   }
   //prep_report(P, termcode);
   //exit(0);
   //termcode = prep_update_mip(mip, *stats, params, etol);
   if(prep_quit(termcode)){
      return termcode;
   }   
   
   if(stats->rows_deleted + 
      stats->vars_fixed + 
      stats->bounds_integerized + 
      stats->coeffs_changed + 
      stats->bounds_tightened > 0){
      return PREP_MODIFIED;
   }   

   return PREP_UNMODIFIED;

   /* exit basic preprocessor */

}

/*===========================================================================*/
/*===========================================================================*/

int prep_delete_duplicate_rows_cols(PREPdesc *P, char check_rows, 
				    char check_cols){

   int termcode = PREP_UNMODIFIED;

   if(!check_cols && !check_rows){
      return termcode;
   }
 
   int i, j, k, l, delete_ind, l_ind, r_ind, cr_ind, cl_ind;
   int obj_ind, col_ind, row_ind, end, obj_size, row_size, delete_row_ind;
   char can_iterate, in_conflict, delete_row; 
   int fix_type, diff_cnt, diff_ind;
   double new_bound, diff_val, diff_obj_val, diff_row_val, rhs_obj, rhs_row;
   double delete_val;
   //double mark_time, start_time = wall_clock(NULL); 

   MIPdesc *mip = P->mip;
   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   double etol = P->params.etol;
   int dive_level = 0; //P->params.dive_level;
   int verbosity = P->params.verbosity;
   prep_stats *stats = &(P->stats);
   
   int m = mip->m;
   int n = mip->n;

   int *matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;

   double *ub = mip->ub, new_row_ub;
   double *lb = mip->lb, new_row_lb;
   
   char *sense = mip->sense;
   double *rhs = mip->rhs;
   double *obj = mip->obj;

   double *col_sum = NULL, *col_factor = NULL;
   double *row_sum = NULL, *row_factor = NULL;
   int last_lloc, last_rloc, *r_loc = NULL, *c_loc = NULL;
 
   int * col_del_ind = NULL, col_del_cnt = 0;
   int * col_fix_type = NULL, dup_type;   
   double *col_fix_val = NULL;
   char *col_orig_type = NULL, type_l, type_r, bin_type; 
   double obj_l, obj_r, obj_diff; 

   if(check_rows){
      col_factor = (double *)malloc(n*DSIZE);
      row_sum = (double *)calloc(m,DSIZE);
      for(i = 0; i < n; i++){
	 col_factor[i] = 1 + (double(rand()) / 
			      RAND_MAX * (10 - 1));;
	 if((double(rand()) / RAND_MAX) < 0.5) col_factor[i] 
						  = -col_factor[i];
      }
      
      r_loc = (int *)malloc(m*ISIZE);
      memcpy(r_loc, P->user_row_ind, ISIZE*m);
  }

   if(check_cols){   
      col_del_ind = (int *)malloc(n*ISIZE);
      col_fix_type = (int *)malloc(n*ISIZE);
      col_fix_val = (double *)malloc(n*DSIZE);
      col_orig_type = (char *)malloc(n*CSIZE);
      row_factor = (double *)malloc(m*DSIZE);
      col_sum = (double *)calloc(n,DSIZE);

      for(i = 0; i < m; i++){
	 row_factor[i] = 1 + (double(rand()) / 
			      RAND_MAX * (10 - 1));;
	 if((double(rand()) / RAND_MAX) < 0.5) row_factor[i] = 
						  -row_factor[i];
      }

      c_loc = (int *)malloc(n*ISIZE);
      memcpy(c_loc, P->user_col_ind, ISIZE*n);
   }

   if(check_rows && check_cols){
      for(col_ind = 0; col_ind < n; col_ind++){
	 end = matbeg[col_ind + 1];
	 for(j = matbeg[col_ind]; j < end; j++){
	    row_ind = matind[j];
	    row_sum[row_ind] += matval[j]*col_factor[col_ind];
	    col_sum[col_ind] += matval[j]*row_factor[row_ind];
	 }
      }
   }else if(check_rows){
      for(col_ind = 0; col_ind < n; col_ind++){
	 end = matbeg[col_ind + 1];
	 for(j = matbeg[col_ind]; j < end; j++){
	    row_ind = matind[j];
	    row_sum[row_ind] += matval[j]*col_factor[col_ind];
	 }
      }
   }else{
      for(col_ind = 0; col_ind < n; col_ind++){
	 end = matbeg[col_ind + 1];
	 for(j = matbeg[col_ind]; j < end; j++){
	    row_ind = matind[j];
	    col_sum[col_ind] += matval[j]*row_factor[row_ind];
	 }
      }
   }

   /* first columns */
   /* for same cols, 
      -if objs are same, then aggregate them (though for binary, 
            check if they conflict first)
      -otherwise, 			 
      -for binary cols, check if they conflict
      -for others, we can delete the duplicate col only under some 
      specific requirements 
    */

   if(check_cols){
      qsort_di(col_sum, c_loc, n);      
      last_lloc = last_rloc = 0;
      while(TRUE){
	 if(last_lloc == n - 1){
	    break;
	 }
	 /* search for same cols */
	 for(i = last_lloc; i < n - 1; i++){
	    if(cols[c_loc[i]].var_type == 'F' ||
	       cols[c_loc[i]].var_type == 'U' ||
	       cols[c_loc[i]].var_type == 'L'){
	       continue;
	    }
	    if(prep_is_equal(col_sum[i], col_sum[i+1], etol)){
	       last_rloc = i+1;
	       if( i < n - 2 ){
		  for(j = i+2; j < n; j++){
		     if(!prep_is_equal(col_sum[i], col_sum[j], etol)){
			last_rloc = j;
			break;
		     }
		  }
	       }
	       break;
	    }	    
	 }
	 
	 if(i == n - 1){
	    break;
	 }

	 /* now we got 2 candidate cols*/
	 l_ind = i;
	 r_ind = l_ind + 1;
	 last_lloc = last_rloc;

	 //printf("starting while loop - cols section \n");  	 
	 while(l_ind < last_rloc - 1){
	    
	    cl_ind = c_loc[l_ind];
	    cr_ind = c_loc[r_ind];

	    //printf("processing cl_ind, %i cr_ind %i\n", cl_ind, cr_ind);
	    
	    if(r_ind == last_rloc || cols[cl_ind].var_type == 'F' ||
	       cols[cl_ind].var_type == 'U' ||
	       cols[cl_ind].var_type == 'L'){
	       //cols[cl_ind].var_type != 'B'){
	       l_ind++;
	       r_ind = l_ind + 1;
	       continue;
	    }
	    
	    //if(cols[cr_ind].var_type != 'B'){
	    if(cols[cr_ind].var_type == 'F' ||
	       cols[cl_ind].var_type == 'U' ||
	       cols[cl_ind].var_type == 'L'){
	       r_ind++;
	       continue;
	    }

	    if(cols[cl_ind].col_size == 0){
	       /* fix this here */
	       col_orig_type[col_del_cnt] = cols[cl_ind].var_type;
	       cols[cl_ind].var_type = 'F';
	       col_del_ind[col_del_cnt] = cl_ind;	       
	       col_fix_type[col_del_cnt] = FIX_OTHER;
	       if(obj[cl_ind] >= 0){
		  col_fix_val[col_del_cnt] = lb[cl_ind];
	       }else{
		  col_fix_val[col_del_cnt] = ub[cl_ind];
	       }
	       col_del_cnt++;
	       l_ind++;
	       r_ind = l_ind + 1;
	       continue;
	    }
	    
	    if(cols[cr_ind].col_size == 0){
	       /* fix this here */
	       col_orig_type[col_del_cnt] = cols[cr_ind].var_type;
	       cols[cr_ind].var_type = 'F';
	       col_del_ind[col_del_cnt] = cr_ind;	       
	       col_fix_type[col_del_cnt] = FIX_OTHER;
	       if(obj[cr_ind] >= 0){
		  col_fix_val[col_del_cnt] = lb[cr_ind];
	       }else{
		  col_fix_val[col_del_cnt] = ub[cr_ind];
	       }
	       col_del_cnt++;
	       r_ind++;
	       continue;
	    }

	    /* if sizes are different, then they are definitely not 
	       same cols 
	       -remember, col_size = real_size - redundant_rows	       
	    */
	    /* also I dont want to mess with diff type of columns now, 
	       so skip if one is int and the other is cont column */

	    type_l = cols[cl_ind].var_type;
	    type_r = cols[cr_ind].var_type;

	    if((cols[cl_ind].col_size != cols[cr_ind].col_size) || 
	       (type_l == 'C' && type_r != 'C') ||
	       (type_l != 'C' && type_r == 'C')){
	       r_ind++;
	       continue;
	    }

	    /* now we have two cols with same size, but this is not enough, 
	       we can iterate only if we have two binary cols or one of 
	       the following conditions are satisfied */

	    obj_l = obj[cl_ind];
	    obj_r = obj[cr_ind];
	    obj_diff = obj_l - obj_r;
	    
	    if(type_l == 'B' && type_r == 'B'){
	       bin_type = TRUE;
	    }else{
	       bin_type = FALSE;
	    }

	    dup_type = -1;
	    
	
    /* 
	       dup_type: -1 -> others, 
	                  0 -> obj_l = obj_r, obj_diff = 0,

	                  1 -> obj_l, obj_r >= 0, obj_diff > 0, l->lb
                               lb[cl_ind] > -INF, ub[cr_ind] = INF
                           
	                  2 -> obj_l, obj_r <= 0, obj_diff > 0, r->lb
                               lb[cl_ind] = -INF, ub[cr_ind] < INF

	                  3 -> obj_1, obj_r >= 0, obj_diff < 0, r->ub
                               ub[cl_ind] = INF, lb[cr_ind] > -INF

	                  4 -> obj_l, obj_r <= 0, obj_diff < 0, l->ub
                               ub[cl_ind] < INF, lb[cr_ind] = -INF

			      --- unbounded cases -- 

			  5 -> ub[cl_ind] = INF, lb[cr_ind] = -INF    
                               obj_l <= 0, obj_r >= 0

		          6 -> ub[cl_ind] = INF, lb[cr_ind] = -INF
			       obj_l >= 0, obj_r >= 0, obj_diff < 0

		          7 -> ub[cl_ind] = INF, lb[cr_ind] = -INF
			       obj_l <= 0, obj_r <= 0, obj_diff < 0

                         
		          8 -> lb[cl_ind] = -INF, ub[cr_ind] = INF
                               obj_l >= 0, obj_r <= 0 

		          9 -> lb[cl_ind] = -INF, ub[cr_ind] = INF
                               obj_l >= 0, obj_r >= 0, obj_diff > 0 

		         10 -> lb[cl_ind] = -INF, ub[cr_ind] = INF
                               obj_l <= 0, obj_r <= 0, obj_diff > 0 

	    */
	    /* get dup type */

	    if(obj_l == obj_r){ //fixme ??
	       /* fixme, make this more efficient */
	       //if(prep_is_equal(obj_r, obj_l, etol)){
	       dup_type = 0;
	    }else if(!bin_type){
	       /* at least one inf */
	       /* fix this ugly thing here */
	       if(ub[cl_ind] >= INF || lb[cl_ind] <= -INF || 
		  ub[cr_ind] >= INF || lb[cl_ind] <= -INF){
		  if(obj_diff > 0.0){
		     if(lb[cl_ind] > -INF && ub[cr_ind] >= INF){
			if(obj_l >= 0 && obj_r >= 0) dup_type = 1;
		     }else if(lb[cl_ind] <= -INF && ub[cr_ind] < INF){
			if(obj_l <= 0.0 && obj_r <= 0.0) dup_type = 2;
		     }else if(lb[cl_ind] <= -INF && ub[cr_ind] >= INF){
			if(obj_l >= 0.0 && obj_r <= 0.0) dup_type = 8;
			else if(obj_l >= 0.0 && obj_r >= 0.0) dup_type = 9;
			else if(obj_l <= 0.0 && obj_r <= 0.0) dup_type = 10;
		     }
		  }else{
		     if(ub[cl_ind] >= INF && lb[cr_ind] > -INF){
			if(obj_l >= 0 && obj_r >= 0) dup_type = 3;
		     }else if(ub[cl_ind] < INF && lb[cr_ind] <= -INF){
			if(obj_l <= 0.0 && obj_r <= 0.0) dup_type = 4;
		     }else if(ub[cl_ind] >= INF && lb[cr_ind] <= -INF){
			if(obj_l <= 0.0 && obj_r >= 0.0) dup_type = 5;
			else if(obj_l >= 0.0 && obj_r >= 0.0) dup_type = 6;
			else if(obj_l <= 0.0 && obj_r <= 0.0) dup_type = 7;
		     }
		  }
	       }
	    }
	    
	    /* fixme - check other cases! 
	       cant iterate in this case for now */
	    if(dup_type < 0 && !bin_type){
	       r_ind++;
	       continue;
	    }

	    /* so now, we have a case*/
	    /* first check if we have same cols here */
	    /* also check conflict if it is binary */

	    in_conflict = FALSE;
	    can_iterate = TRUE;
	    
	    for( k = matbeg[cl_ind], l = matbeg[cr_ind];;){
	       if(k < matbeg[cl_ind + 1] && 
		  (matind[k] < matind[l] ||
		   l >= matbeg[cr_ind + 1])){
		  if(!rows[matind[k]].is_redundant){
		     can_iterate = FALSE;
		     break;
		  }
		  k++;
	       }else if(l < matbeg[cr_ind + 1] && 
			(matind[k] > matind[l] ||
			 k >= matbeg[cl_ind+1])){ 
		  if(!rows[matind[l]].is_redundant){
		     can_iterate = FALSE;
		     break;
		  }	 
		  l++;
	       }else{
		  if(!rows[matind[l]].is_redundant){
		     if(!prep_is_equal(matval[l], matval[k], etol)){
			can_iterate = FALSE;
			break;
		     }	       
		     if(bin_type){
			if(!in_conflict){
			   new_row_lb = rows[matind[l]].lb;
			   new_row_ub = rows[matind[l]].ub;
			   if(matval[l] > 0.0){
			      new_row_lb += matval[l];
			   }else{
			      new_row_ub += matval[l];
			   }
			   if(matval[k] > 0.0){
			      new_row_lb += matval[k];
			   }else{
			      new_row_ub += matval[k];
			   }
			   
			   switch(sense[matind[l]]){
			    case 'E':
			       if(new_row_lb > rhs[matind[l]] + etol ||
				  new_row_ub < rhs[matind[l]] - etol){
				  in_conflict = TRUE;
			       }
			       break;
			    case 'L':
			       if(new_row_lb > rhs[matind[l]] + etol){
				  in_conflict = TRUE;
			       }
			       break;
			   }
			}
		     }		     
		  } 
		  k++;
		  l++;
	       }
	       if((k == matbeg[cl_ind + 1] && l == matbeg[cr_ind + 1])){
		  break;
	       }
	    }
	    
	    if(!can_iterate){
	       /* so columns are not same, go back */
	       r_ind++;
	       continue;
	    }
	    
	    /* so we have same columns, check if we can aggregate/ 
	       delete one of them */

	    /* first check if we have conflict in binary case */
	    delete_ind = -1;
	    if(bin_type && in_conflict){
	       if(obj_diff > 0.0){
		  delete_ind = cl_ind; 
		  delete_val = 0.0;
	       }else{
		  delete_ind = cr_ind;
		  delete_val = 0.0;
	       }
	       fix_type = FIX_BINARY;
	    }else{
	       /* so not binary or are not in conflict */
	       /*check if we can aggregate first*/
	       if(dup_type == 0){
		  /* same obj, same col */
		  /* just merge it to the one on the left and 
		     make the one on the right invisible
		     row bounds wont change in this case */
		  
		  /* first merge */
		  /* just the bounds */
		  if(lb[cl_ind] > -INF){
		     if(lb[cr_ind] <= -INF){
			lb[cl_ind] = -INF;
		     }else{
			lb[cl_ind] += lb[cr_ind];
		     }
		  }
		  if(ub[cl_ind] < INF){
		     if(ub[cr_ind] >= INF){
			ub[cl_ind] = INF;
		     }else{
			ub[cl_ind] += ub[cr_ind];
		     }
		  }

		  if(type_l == 'B'){
		     cols[cl_ind].var_type = 'I';
		  }

		  stats->vars_aggregated++;
		  /* now cleanup*/
		  /* well just fix it to 0.0 so that we 
		     will take care of other attr with this row*/
		  
		  delete_ind = cr_ind;
		  delete_val = 0.0;
		  fix_type = FIX_AGGREGATE;
	       }else{
		  if(bin_type){
		     /* cant do anything with this pair*/
		     r_ind++;
		     continue;
		  }
		  
		  /* now check dup_type and see what we can do */
		  fix_type = FIX_OTHER;
		  switch(dup_type){
		   case 1:
		      delete_ind = cl_ind;
		      delete_val = lb[cl_ind];
		      break;
		   case 2:
		      delete_ind = cr_ind;
		      delete_val = ub[cr_ind];
		      break;
		   case 3:
		      delete_ind = cr_ind;
		      delete_val = lb[cr_ind];
		      break;
		   case 4:
		      delete_ind = cl_ind;
		      delete_val = ub[cl_ind];
		      break;
		   case 5:
		   case 6:
		   case 10:
		      stats->col_unbound_ind = cr_ind;
		      return PREP_UNBOUNDED;
		   case 7:
		   case 8:
		   case 9:
		      stats->col_unbound_ind = cl_ind;
		      return PREP_UNBOUNDED;
		   default:
		      printf("error in prep_delete_duplicate_row_cols()\n");
		      return PREP_OTHER_ERROR;
		  }
	       }
	    }
		  
	    if(delete_ind >= 0){
	       if(delete_ind == cl_ind){
		  l_ind++;
		  r_ind = l_ind + 1;
	       }else{
		  r_ind++;
	       }
	       col_orig_type[col_del_cnt] = cols[delete_ind].var_type;
	       cols[delete_ind].var_type = 'F';
	       col_del_ind[col_del_cnt] = delete_ind;
	       col_fix_val[col_del_cnt] = delete_val;
	       col_fix_type[col_del_cnt++] = fix_type;	       
	    }
	 }
      }

      /* ok, now fix each of these duplicate cols */
      for(i = 0; i < col_del_cnt; i++){
	 cols[col_del_ind[i]].var_type = col_orig_type[i];
	 termcode = prep_modified_cols_update_info(P, 1, &col_del_ind[i], 
						   -1, dive_level, 
						   col_fix_val[i], 
						   col_fix_type[i], 
						   TRUE, FALSE);
	 if(prep_quit(termcode)){
	    return termcode;
	 }
      }
   }
   
   //  printf("Aggre cnt: %i\n", stats->vars_aggregated);

#if 0
   mark_time = wall_clock(NULL);
   printf("Total duplicate cols Time: %f...\n\n", 
	  mark_time - start_time);    
#endif   

   /* now same for rows */

   if(check_rows){
      qsort_di(row_sum, r_loc, m);
      last_lloc = last_rloc = 0;   
      while(TRUE && check_rows){
	 
	 if(last_lloc == m - 1){
	    break;
	 }
	 
	 for(i = last_lloc; i < m - 1; i++){
	    if(prep_is_equal(row_sum[i], row_sum[i+1], etol)){
	       last_rloc = i+1;
	       if( i < m - 2 ){
		  for(j = i+2; j < m; j++){
		     if(!prep_is_equal(row_sum[i], row_sum[j], etol)){
			last_rloc = j;
			break;
		     }
		  }
	       }
	       break;
	    }	    
	 }
	 
	 if(i == m - 1){
	    break;
	 }
	 
	 l_ind = i;
	 r_ind = l_ind + 1;
	 last_lloc = last_rloc;
	 while(l_ind < last_rloc - 1){
	    
	    obj_ind = r_loc[l_ind];
	    row_ind = r_loc[r_ind];
	    
	    if(r_ind == last_rloc || 
	       rows[obj_ind].is_redundant){
	       l_ind++;
	       r_ind = l_ind + 1;
	       continue;
	    }
	    
	    if(rows[row_ind].is_redundant){
	       r_ind++;
	       continue;
	    }
	    
	    obj_size = rows[obj_ind].size - rows[obj_ind].fixed_var_num;
	    row_size = rows[row_ind].size - rows[row_ind].fixed_var_num;
	    
	    if(obj_size - row_size > 2 ||
	       obj_size - row_size < -2){
	       r_ind++;
	       continue;
	    }
	    
	    delete_row = FALSE;
	    
	    if(obj_size == 0){
	       delete_row = TRUE;
	       delete_row_ind = obj_ind;
	       l_ind++;
	       r_ind = l_ind + 1;
	    }
	    
	    if(!delete_row){
	       if(row_size == 0){
		  delete_row = TRUE;
		  delete_row_ind = row_ind;
		  r_ind++;
	       }
	    }
	    
	    if(!delete_row){
	       
	       /* now check if rows are same */
	       diff_cnt = 0;
	       diff_ind = 0;
	       diff_obj_val = 0;
	       diff_row_val = 0; 
	       
	       for( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
		  if(k < r_matbeg[obj_ind + 1] && 
		     (r_matind[k] < r_matind[l] ||
		      l >= r_matbeg[row_ind + 1])){
		     
		     if(cols[r_matind[k]].var_type != 'F'){
			diff_ind = r_matind[k];
			diff_obj_val = r_matval[k];
			diff_cnt++;
		     }
		     k++;
		  }else if(l < r_matbeg[row_ind + 1] && 
			   (r_matind[k] > r_matind[l] ||
			    k >= r_matbeg[obj_ind+1])){ 
		     if(cols[r_matind[l]].var_type != 'F'){
			diff_ind = r_matind[l];
			diff_row_val = r_matval[l];
			diff_cnt++;
		     }	 
		     l++;
		  }else{
		     if(cols[r_matind[l]].var_type != 'F'){
			if(!prep_is_equal(r_matval[l], r_matval[k], etol)){
			   diff_ind = r_matind[k];
			   diff_obj_val = r_matval[k];
			   diff_row_val = r_matval[l];
			   diff_cnt++;
			}	 	
		     } 
		     k++;
		     l++;
		  }
		  if(diff_cnt > 1 || 
		     (k == r_matbeg[obj_ind + 1] && 
		      l == r_matbeg[row_ind + 1])){
		     break;
		  }
	       }
	       
	       if(diff_cnt < 2 &&(diff_cnt == 0 || 
				  prep_is_equal(diff_obj_val - 
						diff_row_val, 
						0.0, 
						etol))){
		  rhs_obj = rhs[obj_ind] - rows[obj_ind].fixed_lhs_offset;
		  rhs_row = rhs[row_ind] + rows[row_ind].fixed_lhs_offset;
		  delete_row = TRUE;
		  if(sense[obj_ind] == 'E'){
		     if(sense[row_ind] == 'E'){
			if(!prep_is_equal(rhs_obj, rhs_row, etol)){
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
		     }else{
			if(rhs_row < rhs_obj - etol){
			   stats->row_infeas_ind = obj_ind;
			   return PREP_INFEAS;
			}
		     }
		     delete_row_ind = row_ind;
		     r_ind++;		  
		  }else{
		     if(sense[row_ind] == 'E'){
			if(rhs_row > rhs_obj + etol){
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
			delete_row_ind = obj_ind;
			l_ind++;		  
			r_ind = l_ind + 1;
		     }else{
			if(rhs_row < rhs_obj - etol){
			   delete_row_ind = obj_ind;
			   l_ind++;
			   r_ind = l_ind + 1;
			}else{
			   delete_row_ind = row_ind;
			   r_ind++;
			}
		     }
		  }
	       }else if(diff_cnt == 1){
		  rhs_obj = rhs[obj_ind] - rows[obj_ind].fixed_lhs_offset;
		  rhs_row = rhs[row_ind] + rows[row_ind].fixed_lhs_offset;
		  diff_val = diff_obj_val - diff_row_val;
		  new_bound = (rhs_obj - rhs_row)/diff_val;
		  
		  if(sense[obj_ind] == 'E'){
		     if(sense[row_ind] == 'E'){
			if(obj_size > row_size){
			   delete_row_ind = row_ind;
			   r_ind++;
			}else{
			   delete_row_ind = obj_ind;
			   l_ind++;
			   r_ind = l_ind + 1;
			}	
			//printf("obj-row:%i %i\n", obj_ind, row_ind); 
			termcode = prep_modified_cols_update_info(P, 1, 
								  &diff_ind, 
								  -1, 
								  dive_level, 
								  new_bound, 
								  FIX_OTHER, 
								  TRUE, FALSE);
			if(prep_quit(termcode)){
			   return termcode;
			}	       
			delete_row = TRUE;
		     }else{
			fix_type = FIX_NO_BOUND;
			if(obj_size > row_size){
			   if(diff_val > 0.0){
			      if(lb[diff_ind] < new_bound - etol){
				 fix_type = IMPROVE_LB;
			      }
			   }else{
			      if(ub[diff_ind] > new_bound + etol){
				 fix_type = IMPROVE_UB;
			      }
			   }
			}else{
			   if(diff_val > 0.0){
			      fix_type = IMPROVE_UB;
			   }else{
			      fix_type = IMPROVE_LB;
			   }
			}
			if(fix_type == IMPROVE_UB){
			   if(ub[diff_ind] < new_bound - etol){
			      fix_type = FIX_NO_BOUND;
			   }
			}else if(fix_type == IMPROVE_LB){
			   if(lb[diff_ind] > new_bound + etol){
			      fix_type = FIX_NO_BOUND;
			   }
			}
			if(fix_type != FIX_NO_BOUND){
			   //printf("obj-row:%i %i\n", obj_ind, row_ind); 
			   termcode = 
			      prep_modified_cols_update_info(P, 1, 
							     &diff_ind,
							     -1, 
							     dive_level, 
							     new_bound, 
							     fix_type, 
							     TRUE, FALSE);
			   if(prep_quit(termcode)){
			      return termcode;
			   }
			}
			r_ind++;
		     }
		  }else{
		     if(sense[obj_ind] != 'E'){
			fix_type = FIX_NO_BOUND;
			if(obj_size > row_size){
			   if(diff_val > 0.0){
			      if(ub[diff_ind] > new_bound + etol){
				 fix_type = IMPROVE_UB;
			      }
			   }else{
			      if(lb[diff_ind] < new_bound - etol){
				 fix_type = IMPROVE_LB;
			      }
			   }
			}else{
			   if(diff_val > 0.0){
			      fix_type = IMPROVE_LB;
			   }else{
			      fix_type = IMPROVE_UB;
			   }
			}
			
			if(fix_type == IMPROVE_UB){
			   if(ub[diff_ind] < new_bound - etol){
			      fix_type = FIX_NO_BOUND;
			   }
			}else if(fix_type == IMPROVE_LB){
			   if(lb[diff_ind] > new_bound + etol){
			      fix_type = FIX_NO_BOUND;
			   }
			}
			if(fix_type != FIX_NO_BOUND){
			   //printf("obj-row:%i %i\n", obj_ind, row_ind); 
			   termcode = 
			      prep_modified_cols_update_info(P, 1, 
							     &diff_ind, 
							     -1, 
							     dive_level, 
							     new_bound, 
							     fix_type, 
							     TRUE, FALSE);
			   if(prep_quit(termcode)){
			      return termcode;
			   }
			}
		     }
		     r_ind++;
		  }
	       }else{
		  r_ind++;
	       }
	       
	    }
	    
	    if(delete_row){
	       stats->rows_deleted++;
	       if(verbosity >= 2){
		  prep_declare_redundant_row(rows[delete_row_ind], 
					     delete_row_ind, 
					     sense[delete_row_ind], 
					     rhs[delete_row_ind]);
	       }
	       termcode = prep_deleted_row_update_info(mip, delete_row_ind);
	       if(prep_quit(termcode)){
		  return termcode;
	       }
	    }
	 }
      }
   }

#if 0
   printf("Total duplicate rows Time: %f...\n\n", 
	  wall_clock(NULL) - mark_time);    
#endif
 
   if(check_cols){
      FREE(col_orig_type);
      FREE(col_del_ind);
      FREE(col_fix_type);
      FREE(col_fix_val);

      FREE(col_sum);
      FREE(col_factor);
      FREE(c_loc);
   }
   if(check_rows){
      FREE(row_sum);
      FREE(row_factor);
      FREE(r_loc);
   }
   
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
#if 0
int prep_delete_duplicate_rows(PREPdesc *P){

   int termcode = PREP_MODIFIED;
   char delete_row; 
   int k, l, diff_ind, diff_cnt, obj_size, row_size;
   int fix_type, delete_row_ind; 
   double new_bound, rhs_obj, rhs_row, diff_val, diff_obj_val, diff_row_val;
   int obj_ind, row_ind;
   
   MIPdesc *mip = P->mip;
   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   prep_stats *stats = &(P->stats);
   double etol = P->params.etol;
   int dive_level = P->params.dive_level;
   int verbosity = P->params.verbosity;
   int m = mip->m;
   
   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;
   
   double *ub = mip->ub;
   double *lb = mip->lb;
   
   char *sense = mip->sense;
   double *rhs = mip->rhs;


   /* delete duplicate rows if there are any */
   for(obj_ind = 0, row_ind = 1; obj_ind < m - 1;){
      //printf("obj-row:%i %i\n", obj_ind, row_ind); 

      if(row_ind == m){
	 obj_ind++;
	 row_ind = obj_ind + 1;
	 continue;
      }

      if(obj_ind == row_ind){
	 row_ind ++;
	 continue;
      }

      if(rows[obj_ind].is_redundant){
	 obj_ind++;
	 row_ind = obj_ind + 1;
	 continue;
      }
      if(rows[row_ind].is_redundant){
	 row_ind++;
	 continue;
      }


      diff_cnt = 0;
      diff_ind = 0;
      diff_obj_val = 0;
      diff_row_val = 0;


      obj_size = rows[obj_ind].size - rows[obj_ind].fixed_var_num;
      row_size = rows[row_ind].size - rows[row_ind].fixed_var_num;
      
      if(obj_size == 0 || row_size == 0){
	 printf("error in duplicacy check: prep_basic()\n");
	 return PREP_OTHER_ERROR;
      }

      if(obj_size - row_size > 2 ||
	 obj_size - row_size < -2){
	 row_ind++;
	 continue;
      }

      for( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
	 if(k < r_matbeg[obj_ind + 1] && 
	    (r_matind[k] < r_matind[l] ||
	     l >= r_matbeg[row_ind + 1])){

	    if(cols[r_matind[k]].var_type != 'F'){
	       diff_ind = r_matind[k];
	       diff_obj_val = r_matval[k];
	       diff_cnt++;
	    }
	 k++;
      }else if(l < r_matbeg[row_ind + 1] && 
	       (r_matind[k] > r_matind[l] ||
		k >= r_matbeg[obj_ind+1])){ 
	 if(cols[r_matind[l]].var_type != 'F'){
	    diff_ind = r_matind[l];
	    diff_row_val = r_matval[l];
	    diff_cnt++;
	 }	 
	 l++;
      }else{
	 if(cols[r_matind[l]].var_type != 'F'){
	    if(!prep_is_equal(r_matval[l], r_matval[k], etol)){
	       diff_ind = r_matind[k];
	       diff_obj_val = r_matval[k];
	       diff_row_val = r_matval[l];
	       diff_cnt++;
	    }	 	
	 } 
	 k++;
	 l++;
      }
	 if(diff_cnt > 1 || 
	    (k == r_matbeg[obj_ind + 1] && l == r_matbeg[row_ind + 1])){
	    break;
	 }
      }

      delete_row = FALSE;

      if(diff_cnt < 2 &&(diff_cnt == 0 || prep_is_equal(diff_obj_val - 
							diff_row_val, 0.0, 
							etol))){
	 rhs_obj = rhs[obj_ind] - rows[obj_ind].fixed_lhs_offset;
	 rhs_row = rhs[row_ind] + rows[row_ind].fixed_lhs_offset;
	 delete_row = TRUE;
	 if(sense[obj_ind] == 'E'){
	    if(sense[row_ind] == 'E'){
	       if(!prep_is_equal(rhs_obj, rhs_row, etol)){
		  stats->row_infeas_ind = row_ind;
		  return PREP_INFEAS;
	       }
	    }else{
	       if(rhs_row < rhs_obj - etol){
		  stats->row_infeas_ind = obj_ind;
		  return PREP_INFEAS;
	       }
	    }
	    delete_row_ind = row_ind;
	    row_ind++;		  
	 }else{
	    if(sense[row_ind] == 'E'){
	       if(rhs_row > rhs_obj + etol){
		  stats->row_infeas_ind = row_ind;
		  return PREP_INFEAS;
	       }
	       delete_row_ind = obj_ind;
	       obj_ind++;		  
	       row_ind = obj_ind + 1;
	    }else{
	       if(rhs_row < rhs_obj - etol){
		  delete_row_ind = obj_ind;
		  obj_ind++;
		  row_ind = obj_ind + 1;
	       }else{
		  delete_row_ind = row_ind;
		  row_ind++;
	       }
	    }
	 }
      }else if(diff_cnt == 1){
	 rhs_obj = rhs[obj_ind] - rows[obj_ind].fixed_lhs_offset;
	 rhs_row = rhs[row_ind] + rows[row_ind].fixed_lhs_offset;
	 diff_val = diff_obj_val - diff_row_val;
	 new_bound = (rhs_obj - rhs_row)/diff_val;
	    
	 if(sense[obj_ind] == 'E'){
	    if(sense[row_ind] == 'E'){
	       if(obj_size > row_size){
		  delete_row_ind = row_ind;
		  row_ind++;
	       }else{
		  delete_row_ind = obj_ind;
		  obj_ind++;
		  row_ind = obj_ind + 1;
	       }	
	       //printf("obj-row:%i %i\n", obj_ind, row_ind); 
	       termcode = prep_modified_cols_update_info(P, 1, &diff_ind, 
						    -1, dive_level, 
						    new_bound, 
						    FIX_OTHER, 
						    TRUE, FALSE);
	       if(prep_quit(termcode)){
		  return termcode;
	       }	       
	       delete_row = TRUE;
	    }else{
	       fix_type = FIX_NO_BOUND;
	       if(obj_size > row_size){
		  if(diff_val > 0.0){
		     if(lb[diff_ind] < new_bound - etol){
			fix_type = IMPROVE_LB;
		     }
		  }else{
		     if(ub[diff_ind] > new_bound + etol){
			fix_type = IMPROVE_UB;
		     }
		  }
	       }else{
		  if(diff_val > 0.0){
		     fix_type = IMPROVE_UB;
		  }else{
		     fix_type = IMPROVE_LB;
		  }
	       }
	       if(fix_type == IMPROVE_UB){
		  if(ub[diff_ind] < new_bound - etol){
		     fix_type = FIX_NO_BOUND;
		  }
	       }else if(fix_type == IMPROVE_LB){
		  if(lb[diff_ind] > new_bound + etol){
		     fix_type = FIX_NO_BOUND;
		  }
	       }
	       if(fix_type != FIX_NO_BOUND){
		  //printf("obj-row:%i %i\n", obj_ind, row_ind); 
		  termcode = prep_modified_cols_update_info(P, 1, &diff_ind, 
						       -1, dive_level, 
						       new_bound, 
						       fix_type, 
						       TRUE, FALSE);
		  if(prep_quit(termcode)){
		     return termcode;
		  }
	       }
	       row_ind++;
	    }
	 }else{
	    if(sense[obj_ind] != 'E'){
	       fix_type = FIX_NO_BOUND;
	       if(obj_size > row_size){
		  if(diff_val > 0.0){
		     if(ub[diff_ind] > new_bound + etol){
			fix_type = IMPROVE_UB;
		     }
		  }else{
		     if(lb[diff_ind] < new_bound - etol){
			fix_type = IMPROVE_LB;
		     }
		  }
	       }else{
		  if(diff_val > 0.0){
		     fix_type = IMPROVE_LB;
		  }else{
		     fix_type = IMPROVE_UB;
		  }
	       }

	       if(fix_type == IMPROVE_UB){
		  if(ub[diff_ind] < new_bound - etol){
		     fix_type = FIX_NO_BOUND;
		  }
	       }else if(fix_type == IMPROVE_LB){
		  if(lb[diff_ind] > new_bound + etol){
		     fix_type = FIX_NO_BOUND;
		  }
	       }
	       if(fix_type != FIX_NO_BOUND){
		  //printf("obj-row:%i %i\n", obj_ind, row_ind); 
		  termcode = prep_modified_cols_update_info(P, 1, &diff_ind, 
						       -1, dive_level, 
						       new_bound, 
						       fix_type, 
						       TRUE, FALSE);
		  if(prep_quit(termcode)){
		     return termcode;
		  }
	       }
	    }
	    row_ind++;
	 }
      }else{
	 row_ind++;
      }

      if(delete_row){
	 stats->rows_deleted++;
	 if(verbosity >= 2){
	    prep_declare_redundant_row(rows[delete_row_ind], 
				       delete_row_ind, sense[delete_row_ind], 
				       rhs[delete_row_ind]);
	 }
	 termcode = prep_deleted_row_update_info(mip, delete_row_ind);
	 if(prep_quit(termcode)){
	    return termcode;
	 }
      }
   }
   
   return termcode;
}
#endif
/*===========================================================================*/
/*===========================================================================*/

int prep_delete_imp_list(IMPlist **list)
{

   IMPvar *imp_var;
   IMPvar *tmp_var;

   if(*list){
      for(imp_var = (*list)->head; imp_var != 0;){
	 tmp_var = imp_var->right;
	 FREE(imp_var);
	 imp_var = tmp_var;
      }
      FREE(*list);
      *list = 0;
   }   
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
/* only used when this row has only 1 lb_inf_var_num or 1 ub_inf_var_num */
int prep_force_row_bounds(PREPdesc *P, int row_ind, int col_ind, int a_loc) 
{
   //   COLinfo *cols = mip->mip_inf->cols;

   int termcode;
   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;
   
   //int i;// *r_matbeg = mip->row_matbeg;
   //   int *r_matind = mip->row_matind;
   //  double *r_matval = mip->row_matval; 

   //int *matbeg = mip->matbeg;
   //int *matind = mip->matind;

   double *rhs = mip->rhs;
   double *ub = mip->ub;
   double *lb = mip->lb;
   char sense = mip->sense[row_ind];
   
   double new_bound;// old_ub = ub[col_ind];
   //double old_lb = lb[col_ind];
   
   char row_lb_improved = FALSE; 
   char row_ub_improved = FALSE;
   int fix_type = FIX_NO_BOUND;

   //int verbosity = P->params.verbosity;
   //double etol = P->params.etol;
   //prep_stats *stats = &(P->stats);
   double etol = P->params.etol;
   int dive_level = 0;
   if(rows[row_ind].lb <= -INF && rows[row_ind].ub >= INF){
      return PREP_UNMODIFIED;
   }

   double a_val = mip->matval[a_loc];

   /* debug check */

   if(sense != 'E'){
      if(!((a_val > 0.0 && ub[col_ind] >= INF) || 
	   (a_val < 0.0 && lb[col_ind] <= -INF))){
	 printf("error in prep_force_row_bounds()\n");
	 return PREP_OTHER_ERROR;
      }
   }else{
      if(!((a_val > 0.0 && ub[col_ind] >= INF) || 
	   (a_val < 0.0 && lb[col_ind] <= -INF)|| 
	   (a_val < 0.0 && ub[col_ind] >= INF) || 
	   (a_val > 0.0 && lb[col_ind] <= -INF))){
	 printf("error -1 in prep_force_row_bounds()\n");
	 return PREP_OTHER_ERROR;
      }
   }

   if(rows[row_ind].ub_inf_var_num <= 1){
      if(a_val > etol && ub[col_ind] >= INF){
	 if(rows[row_ind].lb > -INF){
	    //ub[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].ub >= INF){
	    // printf("error -2 in prep_force_row_bounds()\n");
	    //  return PREP_OTHER_ERROR;
	    //}	 
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].lb + 
				  lb[col_ind]*a_val)/a_val);
	    
	    fix_type = IMPROVE_UB;
	    row_ub_improved = TRUE;
	 }
      }else if(a_val < -etol && lb[col_ind] <= -INF){
	 if(rows[row_ind].lb > -INF){
	    // lb[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].ub >= INF){
	    //   printf("error -2 in prep_force_row_bounds()\n");
	    //  return PREP_OTHER_ERROR;
	    //}
	    
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].lb + 
				  ub[col_ind]*a_val)/a_val);
	    
	    fix_type = IMPROVE_LB;
	    row_ub_improved = TRUE;
	 }
      }      
   }else if(sense == 'E'){
      if(a_val > etol && lb[col_ind] <= -INF){
	 if(rows[row_ind].ub < INF){
	    
	    //lb[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].lb <= -INF){
	    //   printf("error -2 in prep_force_row_bounds()\n");
	    //   return PREP_OTHER_ERROR;
	    // }
	    
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].ub + 
				  ub[col_ind]*a_val)/a_val);
	    
	    fix_type = IMPROVE_LB;
	    row_lb_improved = TRUE;	    
	 }      
      }else if(a_val < -etol && ub[col_ind] >= INF){
	 if(rows[row_ind].ub < INF){
	    //ub[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].lb <= -INF){
	    //  printf("error -2 in prep_force_row_bounds()\n");
	    //  return PREP_OTHER_ERROR;
	    //}
	    
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].ub + 
				  lb[col_ind]*a_val)/a_val);

	    fix_type = IMPROVE_UB;
	    row_lb_improved = TRUE;
	 }
      }
   } 
   
   if(row_lb_improved || row_ub_improved){
      //ub[col_ind] = old_ub;
      //lb[col_ind] = old_lb;
      
      /* we have already obtained this rows bound, however we need to  
	 apply this bound to other rows */
	 
      //rows[row_ind].is_redundant = TRUE;
      termcode = prep_modified_cols_update_info(P, 1, &col_ind, row_ind, 
					       dive_level, 
					       new_bound, 
					       fix_type, TRUE, FALSE);
				    
      if(prep_quit(termcode)){
	 return termcode;
      }
      
      //rows[row_ind].is_redundant = FALSE;

      /* now back to current row*/
  /*     if(fix_type == IMPROVE_LB){ */
/* 	 lb[col_ind] = new_bound; */
/* 	 if(row_lb_improved){ */
/* 	    rows[row_ind].lb += lb[col_ind]*a_val; */
/* 	 }else{ */
/* 	    rows[row_ind].ub += lb[col_ind]*a_val; */
/* 	 } */
/*       }else{ */
/*       	 ub[col_ind] = new_bound; */
/*       	 if(row_lb_improved){ */
/*       	    rows[row_ind].lb += ub[col_ind]*a_val; */
/*       	 }else{ */
/* 	    rows[row_ind].ub += ub[col_ind]*a_val; */
/*       	 } */
/*       } */
      
/*       if(row_lb_improved){ */
/*       	 rows[row_ind].lb_inf_var_num--; */
/*       }else{ */
/*       	 rows[row_ind].ub_inf_var_num--; */
/*       } */

      return PREP_MODIFIED;
   }
   
   return PREP_UNMODIFIED;

}

/*===========================================================================*/
/*===========================================================================*/
#if 0
int prep_initialize_impl_lists(PREPdesc *P){
   
   int i, j, k;
   MIPdesc *mip = P->mip;
   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   int m = mip->m;
   
   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval; 
   
   double *ub = mip->ub;
   double *lb = mip->lb;
   double *obj = mip->obj;

   char is_int = mip->is_int[col_ind];
   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double a_val = r_matval[a_loc];

   for(i = 0; i < m; i++){
      if(rows[i].type == BIN_CONT_TYPE ||
	 rows[i].type == BIN_INT_TYPE ||
	 rows[i].type == BINARY_TYPE){	 	 
      for(j = r_matbeg[i]; j < r_matbeg[i+1]; j++){
	 col_ind = r_matind[j];
	 a_val = r_matvaj[j];
	 /* first assume it is fixed to its upper bound */
	 if(rows[i].type == BINARY_TYPE){
	    if(rows[i].sign_type == ALL_POS_VEC){
	       if(rows[i].coef_type != FRACTIONAL_VEC){
		  if(rhs[i] - a_val < 1 - etol){
		     /* all others can be fixed to their lower bounds */
		     if(sense[i] == 'L'){
			prep_add_impl_list(cols[col_ind].ulist, IMP_ROW, 
					   i, 'L', -1);
		     }
		     
		  }
	       }
	    }
	 }else{
	    
	    

	 }
      }else if(rows[i].sign_type == ALL_NEG_VEC){
	 
      }
	 
	 
      }
      
      
   }
}
#endif
/*===========================================================================*/
/*===========================================================================*/
int prep_improve_variable(PREPdesc *P, int col_ind, int row_ind, int a_loc, 
			  int dive_level, char check_improve, char impl_mode, 
			  char use_sr_bounds, double sr_ub, double sr_lb, 
			  int use_mip) 
			  
{

   int i, fix_type = FIX_NO_BOUND, termcode = PREP_UNMODIFIED;

   MIPdesc *mip = P->mip;

   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   //ROWinfo row = rows[row_ind];

   //   int *r_matbeg = mip->row_matbeg;
   // int *r_matind = mip->row_matind;

   double *maj_matval;

   if(use_mip == COL_ORDERED){
      maj_matval = mip->matval; 
   }else{
      maj_matval = mip->row_matval;
   }

   double *ub = mip->ub;
   double *lb = mip->lb;
   double *obj = mip->obj;

   char is_int = mip->is_int[col_ind];
   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double a_val = maj_matval[a_loc];
   
   // int *matbeg = mip->matbeg;
   // int *matind = mip->matind;
   // double *matval = mip->matval;
   
   char fix_to_lb, fix_to_ub, improve_coef;
   char col_lb_unbounded, col_ub_unbounded;

   double new_bound, new_ub, new_lb;// improve_offset; 

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   prep_stats *stats = &(P->stats);
   

   if(cols[col_ind].var_type != 'U' &&
      cols[col_ind].var_type != 'L'){
      if(cols[col_ind].col_size <= 1){
	 //if(cols[col_ind].col_size == 1){
	 //a_val = r_matval[a_loc];
	 if(obj[col_ind] >= 0.0 && ( (sense == 'G' && a_val < -etol) || 
				     (sense == 'L' && a_val > etol) ||
				     (cols[col_ind].col_size == 0) ) ){ 
	    if(lb[col_ind] <= -INF){
	       stats->col_unbound_ind = col_ind;
	       return PREP_UNBOUNDED;
	    }else{
	       new_bound = lb[col_ind];
	       if(cols[col_ind].var_type == 'B'){
		  fix_type = FIX_BINARY;
	       }else{
		  fix_type = FIX_OTHER;	    
	       }
	    }
	 }else if(obj[col_ind] <= 0.0 && ( (sense == 'G' && a_val > etol) || 
					   (sense == 'L' && a_val < -etol) || 
					   (cols[col_ind].col_size == 0)) ){
	    if(ub[col_ind] >= INF){
	       stats->col_unbound_ind = col_ind;
	       return PREP_UNBOUNDED;
	    }else{
	       new_bound = ub[col_ind];
	       if(cols[col_ind].var_type == 'B'){
		  fix_type = FIX_BINARY;
	       }else{
		  fix_type = FIX_OTHER;	    
	       }
	    }
	 }
	 if(fix_type != FIX_NO_BOUND){
	    termcode = PREP_MODIFIED;
	 }      
      }

   }

   if(termcode == PREP_UNMODIFIED){
      if(cols[col_ind].var_type == 'B'){
	 
	 /* see if we can fix this variable */
	 /* or improve the coefficients */
	 
	 fix_to_lb = FALSE;
	 fix_to_ub = FALSE;
	 improve_coef = FALSE;
	 //a_val = r_matval[a_loc];
	 if(a_val > etol){	
	    switch(sense){
	     case 'L':
		/* try to fix */
		if(rows[row_ind].lb > -INF){		   
		   if(use_sr_bounds){
		      new_lb = sr_lb;
		   }else{
		      new_lb = rows[row_ind].lb + a_val;
		   }
		   if(new_lb > rhs + etol)
		      fix_to_lb = TRUE;
		}			       	       
		/* try to improve */
		if(!fix_to_lb && check_improve && !impl_mode){
		   if(rows[row_ind].ub < INF){
		      /* this row might have been redundant all the 
			 way up here */
		      
		      //if(rows[row_ind].ub < rhs - etol){
		   // termcode = prep_check_redundancy(P, row_ind, FALSE);
		   // if(prep_quit(termcode)){
		   ///    return termcode;
		   //	 }
		   //	 termcode = PREP_UNMODIFIED;
			 /* do nothing below */
		      //      }else{
		      if(use_sr_bounds){
			 new_ub = sr_ub;
			 if(sr_ub < rhs - etol){
			    new_ub = sr_ub - rhs;
			    maj_matval[a_loc] -= new_ub;
			    mip->rhs[row_ind] -= new_ub;

			    if(prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
#if 0
			       /*fixme! */
			       (rows[row_ind].fixed_var_num)++;
			       (rows[row_ind].bin_var_num)--;
			       if(!prep_is_integral(a_val, etol)){
				  (rows[row_ind].frac_coef_num)--;
			       }
			       (cols[col_ind].col_size)--;

			       if(rows[row_ind].fixed_var_num >=
				  rows[row_ind].size -1){
				  printf("assigned to 0\n");
				  printf("row size 1 %i %i %i\n", row_ind,
					 rows[row_ind].fixed_var_num,
					 rows[row_ind].size);
			       }
#endif
			    }
			    rows[row_ind].ub += 
			       (maj_matval[a_loc] - a_val) * 
			       ub[col_ind];
			    
			    improve_coef = TRUE; 
			 }			 
		      }else{
			 new_ub = rows[row_ind].ub - a_val;
			 if(new_ub < rhs - etol){
			    
			    /* update coefs */			 
			    maj_matval[a_loc] = rows[row_ind].ub - rhs;
			    mip->rhs[row_ind] = new_ub;
			    
			    /* debug */
			    if(maj_matval[a_loc] < -etol){
			       printf("error -0 in prep_improve_variable()\n");
			       return PREP_OTHER_ERROR;
			    } 
			    
			    if(prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
#if 0
			       (rows[row_ind].fixed_var_num)++;
			       (rows[row_ind].bin_var_num)--;
			       if(!prep_is_integral(a_val, etol)){
				  (rows[row_ind].frac_coef_num)--;
			       }
			       (cols[col_ind].col_size)--;

			       if(rows[row_ind].fixed_var_num >=
				  rows[row_ind].size -1){
				  printf("assigned to 0\n");
				  printf("row size 1 %i %i %i\n", row_ind,
					 rows[row_ind].fixed_var_num,
					 rows[row_ind].size);
			       }
#endif
			    }
			    
			    /* update bounds */
			    rows[row_ind].ub += 
			       (maj_matval[a_loc] - a_val) * 
			       ub[col_ind];

			    improve_coef = TRUE; 
			 }
		      }
		   }
		}
		break;
	     case 'G':
		/* debug */
		/* cant have 'G' row */
		printf("error -2 in prep_improve_variable()\n");
		return PREP_OTHER_ERROR;		
#if 0
		if(rows[row_ind].ub < INF){
		   new_ub = rows[row_ind].ub - a_val;
		   if(new_ub < rhs){
		      fix_to_ub = TRUE;
		   }
		}
		if(!fix_to_ub){
		   if(rows[row_ind].lb > -INF){
		      new_lb = rows[row_ind].lb + a_val;
		      if(new_lb > rhs){
			 improve_offset = rhs - new_lb; 
			 improve_coef = TRUE;
		      }
		   }
		}
		break;
#endif
	     case 'E':
		if(rows[row_ind].lb > -INF){
		   if(use_sr_bounds){
		      new_lb = sr_lb;
		   }else{
		      new_lb = rows[row_ind].lb + a_val;		   
		   }
		   if(new_lb > rhs){
		      fix_to_lb = TRUE;
		   }
		}
		if(rows[row_ind].ub < INF){
		   if(use_sr_bounds){
		      new_ub = sr_ub;
		   }else{
		      new_ub = rows[row_ind].ub - a_val;
		   }
		   if(new_ub < rhs){
		      fix_to_ub = TRUE;
		   }
		}
		if(fix_to_lb && fix_to_ub){
		   stats->col_infeas_ind = col_ind;
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
		break;
	    }
	 }else if(a_val < -etol){
	    switch(sense){
	     case 'L':
		if(rows[row_ind].lb > -INF){
		   if(use_sr_bounds){
		      new_lb = sr_lb;
		   }else{
		      new_lb = rows[row_ind].lb - a_val; 	 	     
		   }
		   if(new_lb > rhs)
		      fix_to_ub = TRUE;
		}
		if(!fix_to_ub && check_improve && !impl_mode){
		   if(rows[row_ind].ub < INF){
		      if(use_sr_bounds){
			 new_ub = sr_ub;
			 if(sr_ub < rhs - etol){
			    new_ub = sr_ub - rhs;
			    maj_matval[a_loc] -= new_ub;

			    if(prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
			       // printf("assigned to 0\n");
			    }
			    rows[row_ind].lb += 
			       (maj_matval[a_loc] - a_val) * 
			       ub[col_ind];
			    
			    improve_coef = TRUE; 
			 }			 
		      }else{
			 new_ub = rows[row_ind].ub + a_val;
			 if(new_ub < rhs - etol){
			    //improve_offset = rhs - new_ub; 
			    
			    /* update coef*/			 
			    maj_matval[a_loc] -= new_ub - rhs;
			    
			    /* debug */
			    if(maj_matval[a_loc] > etol){
			       printf("error -3 in prep_improve_variable()\n");
			       return PREP_OTHER_ERROR;
			    } 
			    
			    if(prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
			       //	    printf("assigned to 0\n");
			    }
			    
			    /* update bounds */
			    if(rows[row_ind].lb > -INF){
			       rows[row_ind].lb += 
				  (maj_matval[a_loc] - a_val) * 
				  ub[col_ind];
			    }
			    
			    improve_coef = TRUE;
			 }
		      }
		   }
		}		
		break;
	     case 'G':
		/* debug */
		/* cant have 'G' row */
		printf("error -5 in prep_improve_variable()\n");
		return PREP_OTHER_ERROR;
#if 0
		if(rows[row_ind].ub < INF){
		   new_ub = rows[row_ind].ub + a_val;	 
		   if(new_ub < rhs){
		      fix_to_lb = TRUE;
		   }
		}
		if(!fix_to_lb){
		   if(rows[row_ind].lb > -INF){
		      new_lb = rows[row_ind].lb - a_val;
		      if(new_lb > rhs){
			 improve_offset = rhs - new_lb; 
			 improve_coef = TRUE;
		      }
		   }
		}
#endif
		break;
	     case 'E':
		if(rows[row_ind].lb > -INF){
		   if(use_sr_bounds){
		      new_lb = sr_lb;
		   }else{
		      new_lb = rows[row_ind].lb - a_val; 	 	     
		   }
		   if(new_lb > rhs){
		      fix_to_ub = TRUE;
		   }
		}
		if(rows[row_ind].ub < INF){
		   if(use_sr_bounds){
		      new_ub = sr_ub;
		   }else{
		      new_ub = rows[row_ind].ub + a_val;	 
		   }
		   if(new_ub < rhs){
		      fix_to_lb = TRUE;
		   }
		}
		if(fix_to_lb && fix_to_ub){
		   stats->col_infeas_ind = col_ind;
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
		break;
	    }
	 }
	 
	 if(fix_to_lb || fix_to_ub){
	    if(fix_to_lb){
	       new_bound = 0.0;
	       //	    ub[col_ind] = lb[col_ind] = 0.0;
	    }else{
	       new_bound = 1.0;
	       // ub[col_ind] = lb[col_ind] = 1.0;
	       // rows[row_ind].fixed_lhs_offset += r_matval[a_loc];
	    }	
	    //cols[col_ind].var_type = 'F';
	    //cols[col_ind].fix_row_ind = row_ind;

	    //rows[row_ind].fixed_var_num++;
	    //rows[row_ind].bin_var_num--;
	    
	    //if(rows[row_ind].bin_var_num < 0 || 
	    //   rows[row_ind].fixed_var_num + rows[row_ind].fixable_var_num > 
	    //   rows[row_ind].size){
	    /* debug */
	    // printf("error in fixing vars 1, prep_fix_variable()\n");
	    // return PREP_OTHER_ERROR;
	    //}	 
	    fix_type = FIX_BINARY;
	    termcode = PREP_MODIFIED;
	 }else if(improve_coef){
	    /* so we can improve a_val and rhs */

	    //printf("coef improved:\n");
	    //printf("var %i - %s - old coef: %f row %i old rhs %f\n",
	    //col_ind,  
	    //   mip->colname[col_ind], a_val, row_ind, rhs); 
	 
	    /* need to update row bounds again here */
	    /* debug -fixme */
	    /* i really dont like this brute forcing here, try to fix it*/

	    if(use_mip == COL_ORDERED){	    
	       for(i = mip->row_matbeg[row_ind]; i < 
		      mip->row_matbeg[row_ind + 1]; i++){
		  if(mip->row_matind[i] == col_ind){
		     mip->row_matval[i] = maj_matval[a_loc];
		     break;			       
		  }
	       }
	       /* debug */
	       if(i == mip->row_matbeg[row_ind + 1]){
		  printf("error -1 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;			    
	       }
	    }else{
	       for(i = mip->matbeg[col_ind]; i < 
		      mip->matbeg[col_ind + 1]; i++){
		  if(mip->matind[i] == row_ind){
		     mip->matval[i] = maj_matval[a_loc];
		     break;			       
		  }
	       }
	       /* debug */
	       if(i == mip->matbeg[col_ind + 1]){
		  printf("error -6 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;			    
	       }
	    }

	    if(verbosity >=3){
	       if(mip->colname){
		  prep_declare_coef_change(row_ind, col_ind, 
					   mip->colname[col_ind], 
					   maj_matval[a_loc], 
					   mip->rhs[row_ind]);
	       }else{
		  prep_declare_coef_change(row_ind, col_ind, 
					   NULL,
					   maj_matval[a_loc], 
					   mip->rhs[row_ind]);
	       }
	    }	       
	    if(!(stats->nz_coeff_changed[a_loc])){
	       stats->nz_coeff_changed[a_loc] = TRUE;
	       stats->coeffs_changed++;
	    }

	    /* since row bound(s) has been updated, check for redundancy */
	    /* usually this shouldnt be here, but we have updated row bounds
	       here */

	    termcode = prep_check_redundancy(P, row_ind, FALSE, 0.0, 0.0,
					     impl_mode, dive_level);
	    
	    if(prep_quit(termcode)){
	       return termcode;
	    }else if(cols[col_ind].var_type == 'F'){
	       return PREP_MODIFIED;
	    }

	    fix_type = IMPROVE_COEF; 
	    termcode = PREP_MODIFIED;
	 }else if(FALSE && !impl_mode &&
		  ((a_val > etol && !P->ulist_checked[col_ind]) ||
		   a_val < -etol && !P->llist_checked[col_ind])){
	    //printf("implication of col: %i\n", col_ind);
	    /* for now - just among the binary variables*/
	    /* so cant fix it, cant improve it and binary 
	       try logical fixing - with dive-level = 1 */
	    
	    /* fist copy initial info 
	       fixme! - work on this! 
	    */

	    /* do once for each variable */
	    int impl_dive_level = 2;	    
	    memcpy(P->impl_rows, rows, sizeof(ROWinfo)*mip->m); 
	    memcpy(P->impl_cols, cols, sizeof(COLinfo)*mip->n); 
	    memcpy(P->impl_ub, ub, DSIZE*mip->n);
	    memcpy(P->impl_lb, lb, DSIZE*mip->n);
	    P->impl_stats = P->stats;
	    if(a_val > etol){

	       if(!cols[col_ind].ulist){
		  cols[col_ind].ulist = (IMPlist *)calloc(sizeof(IMPlist),1);
	       }	       

	       P->list = cols[col_ind].ulist;	      
	       P->ulist_checked[col_ind] = TRUE;
	       /* fix it to 1.0 and see if that causes any infeasibility 
		  otherwise get the impllist and continue*/
	       /* get the implication list */

	       //	       termcode = prep_get_impl_list(P, col_ind);
	       /* fix this column, update row bounds of this column
		  check for redundancy, */
	       


	       termcode = prep_modified_cols_update_info(P, 1, &col_ind, 
							row_ind, impl_dive_level,
							1.0,
							IMPROVE_LB, TRUE, TRUE);
	       if(termcode == PREP_INFEAS){
		  printf("infeasibility detected!\n");
		  /*then this column is fixable to its lower bound! */
		  new_bound = 0.0;
		  fix_type = FIX_BINARY;
		  termcode = PREP_MODIFIED;
	       }else{

		  termcode = PREP_UNMODIFIED;
		  /* here check other problems and 
		     check if we have redundancies and modify coefficients */

		  /*FIXME */
		  
	       }
	    }else if (a_val < etol){

	       if(!cols[col_ind].llist){
		  cols[col_ind].llist = (IMPlist *)calloc(sizeof(IMPlist),1);
	       }

	       P->list = cols[col_ind].llist;	      
	       P->llist_checked[col_ind] = TRUE;
	       /* fix it to 1.0 and see if that causes any infeasibility 
		  otherwise get the impllist and continue*/
	       
	       termcode = prep_modified_cols_update_info(P, 1, &col_ind, row_ind,
							impl_dive_level,
							0.0, IMPROVE_UB, TRUE, TRUE);
	       
	       if(termcode == PREP_INFEAS){
		  printf("infeasibility detected!\n");
		  /*then this column is fixable to its lower bound! */
		  new_bound = 1.0;
		  fix_type = FIX_BINARY;
		  termcode = PREP_MODIFIED;
	       }else{
		  termcode = PREP_UNMODIFIED;
		  /* here check if we have redundancies and modify coefficients */
		  
		  /*FIXME */
	       }
	    }

	    /* now get back */
	    memcpy(rows, P->impl_rows,sizeof(ROWinfo)*mip->m); 
	    memcpy(cols, P->impl_cols, sizeof(COLinfo)*mip->n); 
	    memcpy(ub, P->impl_ub, DSIZE*mip->n);
	    memcpy(lb, P->impl_lb, DSIZE*mip->n);
	    P->stats = P->impl_stats;
	 }
      }else if(cols[col_ind].var_type == 'U' ||
	       cols[col_ind].var_type == 'L'){
	 if(cols[col_ind].var_type == 'U'){
	    new_bound = ub[col_ind];
	    if(is_int){
	       new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
	    }
	 }else{
	    new_bound = lb[col_ind];
	    if(is_int){
	       new_bound = prep_rnd_integral(new_bound, etol, RND_CEIL);
	    }
	 }
	 
	 //      cols[col_ind].var_type = 'F';
	 //rows[row_ind].fixed_var_num++;
	 //rows[row_ind].fixable_var_num--;
	 //rows[row_ind].fixed_lhs_offset += r_matval[a_loc] * new_bound;
	 
	 /* debug */
	 //if(rows[row_ind].fixed_var_num + rows[row_ind].fixable_var_num > 
	 // rows[row_ind].size){
	 // printf("error in fixing vars 2, prep_fix_variable()\n");
	 // return PREP_OTHER_ERROR;
	 //}	 
	 fix_type = FIX_OTHER; 
	 termcode = PREP_MODIFIED;

      }else{
	 
	 /* not binary, not fixable etc. */
	 /* now try bounds improvement */

	 col_lb_unbounded = FALSE;
	 col_ub_unbounded = FALSE;
	 //a_val = r_matval[a_loc];
	 
	 if(a_val > etol){
	    if(lb[col_ind] <= -INF){
	       
	       /* debug */
	       if(rows[row_ind].lb > -INF){
		  printf("error -7 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }
	       col_lb_unbounded = TRUE;
	       /* can we fix it? */
	       /* is fixable if sense = 'E' and ub < INF*/
	       if(sense == 'E' && rows[row_ind].ub < INF){
		  new_bound = (double)((rhs - rows[row_ind].ub + 
					a_val*ub[col_ind])/a_val);
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_CEIL);
		  }
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind, 
							   row_ind, dive_level,
							   new_bound,
							   IMPROVE_LB, TRUE,
							    impl_mode);
		  if(prep_quit(termcode)){
		     return termcode;
		  }else if(rows[row_ind].is_redundant){
		     return PREP_MODIFIED;
		  }
		  termcode = PREP_UNMODIFIED;
		  col_lb_unbounded = FALSE;		  
	       }
	    }
	    if(!col_lb_unbounded){
	       if(rows[row_ind].lb > -INF){
		  new_bound = (double)((rhs - rows[row_ind].lb + 
					a_val*lb[col_ind])/a_val);
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
		  }	       
		  
		  if(new_bound < ub[col_ind]){
		     termcode = PREP_MODIFIED;
		     fix_type = IMPROVE_UB;
		     //prep_modified_col_update_info(mip, col_ind, new_bound,
		     //			   0.0, IMPROVE_UB);
		  }
	       }
	    }
	 }else if(a_val < -etol){
	    if(ub[col_ind] >= INF){

	       /* debug */
	       if(rows[row_ind].lb > -INF){
		  printf("error -2 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }
	       
	       col_ub_unbounded = TRUE;
	       /* can we fix it? */
	       /* is fixable if sense = 'E' and ub < INF*/
	       if(sense == 'E' && rows[row_ind].ub < INF){
		  new_bound = (double)((rhs - rows[row_ind].ub + 
					a_val*lb[col_ind])/a_val);
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
		  }
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind, 
							   row_ind, dive_level,
							   new_bound,
							   IMPROVE_UB, TRUE, impl_mode);
		  if(prep_quit(termcode)){
		     return termcode;
		  }else{
		     if(rows[row_ind].is_redundant){
			return PREP_MODIFIED;
		     }
		  }
		  termcode = PREP_UNMODIFIED;
		  col_ub_unbounded = FALSE;		  
	       }
	    }
	    if(!col_ub_unbounded){
	       if(rows[row_ind].lb > -INF){
		  new_bound = (double)((rhs - rows[row_ind].lb + 
					a_val*ub[col_ind])/a_val);
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_CEIL);
		  }	       
		  
		  if(new_bound > lb[col_ind]){
		     termcode = PREP_MODIFIED;
		     fix_type = IMPROVE_LB;
			//prep_modified_col_update_info(mip, col_ind, new_bound,
			//		   0.0, IMPROVE_LB);
		  }
	       }
	    }
	 }
      }
   }
   /* now check if we need to update row bounds */
   if(termcode == PREP_MODIFIED && fix_type != IMPROVE_COEF){
      /* set col type to 'T', set it to F after you have visited all 
	 other rows? */
      /* have col.fix_row_ind to mark on which row you have fixed it? */      
      /* isnt worth it, update row bounds here */
      //cols[col_ind].fix_row_ind = row_ind;
      /* this col might have been improved above */
      if(cols[col_ind].var_type != 'F'){
	 termcode = prep_modified_cols_update_info(P, 1, &col_ind, row_ind, 
						  dive_level,
						  new_bound, 
						  fix_type, TRUE, impl_mode);
	 if(prep_quit(termcode)){
	    return termcode;
	 }else{
	    return PREP_MODIFIED;
	 }	    
      }
   }
   
   
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
   int prep_add_to_impl_list(IMPlist *list, int ind, int fix_type, 
			  double val){
   
   if(!list){
      printf("error in prep_add_to_impl_list\n");
      exit(0);      
   }

   IMPvar * var = (IMPvar *)calloc(sizeof(IMPvar),1);
 
   var->ind = ind;
   var->fix_type = fix_type;
   var->val = val;

   if(!list->head){
      list->head = list->tail = var;      
   }else{
      list->tail->right = var;
      list->tail = var; 
   }

   list->size++;
   return 0;

}
/*===========================================================================*/
/*===========================================================================*/
int prep_modified_cols_update_info(PREPdesc *P, int col_cnt, int *col_start,
				   int row_ind, int dive_level, 
				   double fixed_bound,  int intl_fix_type,
				  char check_redundancy, char impl_mode)
{   
   /* fix_type 
      0 FIX_NO_BOUND
      1 FIX_BINARY
      2 FIX_FIXABLE 
      3 FIX_OTHER
      4 IMPROVE_UB
      5 IMPROVE_LB
      6 IMPROVE_COEF
      7 FIX_ALL_LB
      8 FIX_ALL_UB
   */
   /* debug */
   //printf("working on row %i\n", row_ind);
   
   int termcode = PREP_UNMODIFIED, i, j, k, r_ind, end;
   int col_ind, a_loc_ref = 0;
   MIPdesc *mip = P->mip;
   int *matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;
   char *is_int = mip->is_int;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;
   //double * obj = mip->obj;

   double old_ub;
   double old_lb;

   double *ub = mip->ub;
   double *lb = mip->lb;

   double a_val; 
   char get_row_ubounds;
   char get_row_lbounds;

   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   //double r_etol = P->params.etol;
   prep_stats *stats = &(P->stats);

   //char *is_row_updated;
   int fix_type = 0, row_cnt = 0;
   //int * row_updated_ind;
   char row_updated;
   char *is_row_updated = NULL;
   int *row_updated_ind = NULL; 
   IMPlist *imp_list;
   IMPvar *imp_var;

    /*first, if in impl mode, add these variables to current impl_list */
   /* do this above */
   int can_iterate = FALSE;

#if 0   
   if(fix_type == FIX_BINARY){
      /* if in impl_mode, add this to the current list */
      if(impl_mode){
	 prep_add_to_impl_list(P->list, col_ind, 'F', 
			       fixed_bound);
      }

      /* now imply the logical fixings of this variable if it is binary and 
	 there have been logical fixing so far */      
      
      IMPlist * list = NULL;
      IMPvar * imp;
      if(fixed_bound > 1 - etol){
	 list = cols[col_ind].ulist;
      }else{
	 list = cols[col_ind].llist;
      }

      if(list){
	 if(list->size > 0){
	    for(imp = list->head; imp != NULL; imp=imp->right){
	       if(cols[imp->ind].var_type != 'F'){ /* already fixed... */
		  termcode =
		     prep_modified_cols_update_info(P, imp->ind, mip->m+1,
						    imp->val, 0, FIX_BINARY,
						    impl_mode);   
	       }
	    }
	 }
      }
   }
   
   if(check_redundancy){
      is_row_updated = P->check_rows_updated;
      memset(is_row_updated, FALSE, CSIZE*mip->m);
   }
#endif 
 
   double mark_time = wall_clock(NULL);
   //if(intl_fix_type != FIX_AGGREGATE){
   is_row_updated = (char *)calloc(CSIZE,mip->m);
   row_updated_ind = (int *)malloc(ISIZE*mip->m);
   //}

   P->alloc2_time += wall_clock(NULL) - mark_time;  
  
   if(intl_fix_type == FIX_ROW_LB ||
      intl_fix_type == FIX_ROW_UB){
      a_loc_ref = r_matbeg[row_ind];
      //row_updated_ind = (int*)malloc(ISIZE*mip->m);
   }

   mark_time = wall_clock(NULL);

   for(j = 0; j < col_cnt; j++){
      
      col_ind = col_start[j];      
      
      if(cols[col_ind].var_type == 'F'){
	 if(intl_fix_type != FIX_ROW_LB &&
	    intl_fix_type != FIX_ROW_UB){
	    if(!prep_is_equal(ub[col_ind], fixed_bound, etol)){
	       if(!impl_mode){
		  stats->col_infeas_ind = col_ind;
	       }
	       termcode = PREP_INFEAS;
	       can_iterate = FALSE;
	       break;
	    }
	 }
	 continue;
      }

      fix_type = intl_fix_type;

      old_ub = ub[col_ind];
      old_lb = lb[col_ind];

      /* debug */
#if 0      
      if(a_loc_ref > r_matbeg[row_ind + 1]){
	 printf("error in prep_modified_cols_update_info()\n");
	 exit(0);
      }
#endif

      if(fix_type == FIX_ROW_LB || fix_type == FIX_ROW_UB){
	 a_val = r_matval[a_loc_ref + j];	 
	 if(fix_type == FIX_ROW_LB){
	    if(a_val > etol){
	       fixed_bound = ub[col_ind] = lb[col_ind];
	    }else if (a_val < -etol){
	       fixed_bound = lb[col_ind] = ub[col_ind];
	    }else{
	       continue;
	    }
	    if(cols[col_ind].var_type == 'B'){
	       fix_type = FIX_BINARY;
	    }else{
	       fix_type = FIX_OTHER;
	    }
	 }else{
	    if(a_val > etol){
	       fixed_bound = lb[col_ind] = ub[col_ind];
	    }else if (a_val < -etol){
	       fixed_bound = ub[col_ind] = lb[col_ind];
	    }else{
	       continue;
	    }
	    if(cols[col_ind].var_type == 'B'){
	       fix_type = FIX_BINARY;
	    }else{
	       fix_type = FIX_OTHER;
	    }
	 }
      }else{
	 if(fix_type != IMPROVE_LB && fix_type != FIX_AGGREGATE){
	    if(ub[col_ind] < fixed_bound - etol){
	       if(!impl_mode){
		  stats->col_infeas_ind = col_ind;
		  stats->row_infeas_ind = row_ind;
	       }
	       termcode = PREP_INFEAS;
	       can_iterate = FALSE;
	       break;
	    }  
	    ub[col_ind] = fixed_bound;     
	    //(stats->bounds_tightened)++;
	 }
	 if(fix_type != IMPROVE_UB && fix_type != FIX_AGGREGATE){
	    
	    if(lb[col_ind] > fixed_bound + etol){
	       if(!impl_mode){
		  stats->col_infeas_ind = col_ind;
		  stats->row_infeas_ind = row_ind;
	       }
	       termcode = PREP_INFEAS;
	       can_iterate = FALSE;
	       break;
	    }      
	    lb[col_ind] = fixed_bound;
	    // if(fix_type == IMPROVE_LB){
	    //	 (stats->bounds_tightened)++;
	    //}
	 }
      }
      
      if(verbosity >= 3){
	 if(fix_type == FIX_AGGREGATE){
	    if(mip->colname){
	       printf("var %s [%i] is aggregated: \n", mip->colname[col_ind],
		      col_ind);
	    }else{
	       printf("var [%i] is aggregated: \n", col_ind);	       
	    }
	 }else{
	    if(mip->colname){
	       printf("var %s [%i] bounds are improved: ",
		      mip->colname[col_ind], col_ind); 
	    }else{
	       printf("var [%i] bounds are improved: ", col_ind); 
	    }
	    if(lb[col_ind] > -INF){
	       printf("\t lb:%f", lb[col_ind]);
	    }
	    if(ub[col_ind] < INF){
	       printf("\t ub:%f ", ub[col_ind]);
	    }
	    printf("\n");
	 }
      }      
      
      if(fix_type != FIX_BINARY){	 
	 if(fix_type != FIX_AGGREGATE){
	    if(prep_is_equal(ub[col_ind], lb[col_ind], etol)){
	       if(cols[col_ind].var_type == 'B'){
		  fix_type = FIX_BINARY;	      	    
	       }else{
		  fix_type = FIX_OTHER;
	       }
	       cols[col_ind].var_type = 'F';      
	    } 
	 }else{
	    cols[col_ind].var_type = 'F';      
	 }
      }else{
	 cols[col_ind].var_type = 'F';      
      }      
      
      /* now add to impl list if in impl_mode */

      if(fix_type != FIX_AGGREGATE){
	 if(impl_mode && P->impl_col_ind != col_ind){
	    if(P->list->size < P->impl_limit){
	       prep_add_to_impl_list(P->list, col_ind, fix_type, 
				     fixed_bound);
	    }	 
	 }
      }

      if(cols[col_ind].var_type == 'F'){
	 if(fix_type != FIX_AGGREGATE){
	 
	 /* first see if you can fix any other variables from the 
	    impl list of this variable */
	    if(fix_type == FIX_BINARY){
	       if(lb[col_ind] >= 1.0 - etol){
		  imp_list = cols[col_ind].ulist;
	       }else{
		  imp_list = cols[col_ind].llist;
	       }
	       if(imp_list){
		  if(imp_list->size > 0){
		     for(imp_var = imp_list->head; imp_var != 0; 
			 imp_var = imp_var->right){
			termcode = prep_modified_cols_update_info(P, 1, 
								  &imp_var->ind,
								  -1, 0, 
								  imp_var->val,
								  FIX_BINARY, 
								  FALSE, FALSE);
			if(prep_quit(termcode)){
			   can_iterate = FALSE;
			   break;
			}
		     }
		  }
	       }
	    }

	    if(verbosity >= 2){
	       if(mip->colname){
		  prep_declare_fixed_var(col_ind, mip->colname[col_ind], 
					 ub[col_ind]);
	       }else{
		  prep_declare_fixed_var(col_ind, NULL, 
					 ub[col_ind]);		  
	       }
	    }

	    if(!impl_mode){
	       mip->mip_inf->sum_obj_offset += mip->obj[col_ind]*ub[col_ind];
	    }
	 }
	
	 (stats->vars_fixed)++;
	 
      }else{
	 (stats->bounds_tightened)++;
      }
      
      if(cols[col_ind].col_size == 0 ){
	 continue;
      }else if(cols[col_ind].col_size < 0){
	 printf("error -00 in prep_fixed_col_update_info()\n");
	 termcode = PREP_OTHER_ERROR;
	 can_iterate = FALSE;
	 break;
      }
      
      end = matbeg[col_ind + 1];
      
      for(i = matbeg[col_ind]; i < end; i++){
	 if(!(rows[matind[i]].is_redundant)){	 
	    a_val = matval[i];
#if 0
	    if(a_val ==  0.0){ /* we have set it to 0.0 */
	       continue;
	    }
#endif
	    get_row_ubounds = FALSE;
	    get_row_lbounds = FALSE;
	    r_ind = matind[i];
	    row_updated = FALSE;
	    if(fix_type != IMPROVE_UB && fix_type != IMPROVE_LB){	       
	       rows[r_ind].fixed_var_num++; 
	       if(!is_int[col_ind]){
		  (rows[r_ind].cont_var_num)--;
	       }
	       
	       if(!prep_is_integral(a_val, etol)){
		  (rows[r_ind].frac_coef_num)--;
	       }
	       
	       if(fix_type == FIX_BINARY){
		  rows[r_ind].bin_var_num--;
		  //}else if(fix_type == FIX_FIXABLE){
		  //rows[r_ind].fixable_var_num--;
	       }
	       /* debug */
	       if(rows[r_ind].bin_var_num < 0 || 
		  rows[r_ind].fixable_var_num < 0){
		  printf("error -0 in prep_fixed_col_update_info()\n");
		  termcode = PREP_OTHER_ERROR;
		  break;
	       }
	       if(fix_type != FIX_AGGREGATE){
		  rows[r_ind].fixed_lhs_offset += a_val * fixed_bound;
	       }
	    }
	    
	    if(old_ub >= INF && fix_type != IMPROVE_LB){
	       if(a_val > etol){
		  rows[r_ind].ub_inf_var_num--;
		  if(rows[r_ind].ub_inf_var_num == 0){
		     get_row_ubounds = TRUE;
		  } 
	       }else if(a_val < -etol){
		  rows[r_ind].lb_inf_var_num--;
		  if(rows[r_ind].lb_inf_var_num == 0){
		     get_row_ubounds = TRUE;
		  }
	       }
	    }
	    
	    if(old_lb <= -INF && fix_type != IMPROVE_UB){
	       if(a_val > etol){
		  rows[r_ind].lb_inf_var_num--;
		  if(rows[r_ind].lb_inf_var_num == 0){
		     get_row_lbounds = TRUE;
		  }
	       }else if(a_val < -etol){
		  rows[r_ind].ub_inf_var_num--;
		  if(rows[r_ind].lb_inf_var_num == 0){
		     get_row_lbounds = TRUE;
		  }
	       }
	    }
	    
	    if(fix_type != FIX_AGGREGATE){
	       if(a_val > etol){
		  if(fix_type != IMPROVE_LB){
		     if(rows[r_ind].ub < INF){
			/* debug */
			if(old_ub >= INF){
			   printf("error -1 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if(fixed_bound != old_ub){
			   rows[r_ind].ub += a_val*(fixed_bound - old_ub);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }
		  if(fix_type != IMPROVE_UB && fix_type != FIX_AGGREGATE){
		     if(rows[r_ind].lb > -INF){
			/* debug */
			if(old_lb <= -INF){
			   printf("error -2 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if(fixed_bound != old_lb){
			   rows[r_ind].lb += a_val*(fixed_bound - old_lb);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }
	       }else if(a_val < -etol){
		  if(fix_type != IMPROVE_UB){	    
		     if(rows[r_ind].ub < INF){
			/* debug */
			if(old_lb <= -INF){
			   printf("error -3 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if(fixed_bound != old_lb){
			   rows[r_ind].ub += a_val*(fixed_bound - old_lb);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }

		  if(fix_type != IMPROVE_LB){
		     if(rows[r_ind].lb > -INF){
			/* debug */
			if(old_ub >= INF){
			   printf("error -4 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if(fixed_bound != old_ub){
			   rows[r_ind].lb += a_val*(fixed_bound - old_ub);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }
	       }
	    }
	    
	    /* debug */
	    if(rows[r_ind].fixed_var_num + 
	       rows[r_ind].fixable_var_num > rows[r_ind].size){
	       printf("error in fixing vars 2, prep_fix_variable()\n");
	       termcode = PREP_OTHER_ERROR;
	       break;
	    }	 
	    
	    if(get_row_lbounds || get_row_ubounds){
	       rows[r_ind].is_updated = TRUE;
	       row_updated = TRUE;
	       prep_get_row_bounds(mip, r_ind, etol);
	    }      
	    if(row_updated){
	       if(!is_row_updated[r_ind]){
		  is_row_updated[r_ind] = TRUE;
		  row_updated_ind[row_cnt++] = r_ind;	    
	       }
	    }      
	 }	 
      }      
   }
   if(impl_mode){
      P->impl_cols_time += wall_clock(NULL) - mark_time;
   }
   /* if row_updated_cnt > 0 also just rows updated? this is inefficient*/
   mark_time = wall_clock(NULL);

   if(!prep_quit(termcode) && check_redundancy && 
      fix_type != IMPROVE_LB &&
      fix_type != IMPROVE_UB && fix_type){

      for(i = 0; i < row_cnt; i++){
	 r_ind = row_updated_ind[i];
	 if(rows[r_ind].is_updated && !rows[r_ind].is_redundant){
	    if(impl_mode){
	       printf("processing row %i\n", r_ind);
	    }
	    termcode = prep_check_redundancy(P, r_ind, FALSE,
					     0.0, 0.0, impl_mode, dive_level);	    
	    if(prep_quit(termcode)){
	       break;
	    }
	    
	    rows[r_ind].is_updated = FALSE;
	    
	    /* now do we want to dive on variables of the rows those share
	       a comman variable 
	       with these fixed column(s)? */
	    if(dive_level > 0){
	       for(k = r_matbeg[r_ind]; k < r_matbeg[r_ind + 1]; k++){
		  if(rows[r_ind].is_redundant){
		     break;
		  }
		  col_ind = r_matind[k];
		  //if(rows[r_ind].vars_checked){
		  // break;
		  //}
		  
		  if(cols[col_ind].var_type != 'F'){
		     termcode = prep_improve_variable(P, col_ind, 
						      r_ind, k, 
						      (dive_level - 1), 
						      TRUE, impl_mode, FALSE,
						      0.0,0.0, ROW_ORDERED); 
		     if(prep_quit(termcode)){
			break;
		     }
		  }
	       }
	    }
	 }
	 if(prep_quit(termcode)){
	    break;	 
	 }
      }
   }

   if(impl_mode){
      P->impl_rows_time += wall_clock(NULL) - mark_time;
   }

   FREE(row_updated_ind);
   FREE(is_row_updated);


   if(prep_quit(termcode)){
      return termcode;
   }

   //#endif
#if 0
   if(rows[r_ind].vars_checked){
      rows[r_ind].vars_checked = FALSE;
   }else
      
      rows[r_ind].vars_checked = TRUE;
#endif 
   
   // mip->ub[col_ind] = mip->lb[col_ind] = fixed_bound;      
   
   /*debug */
   return PREP_MODIFIED;
}


/*===========================================================================*/
/*===========================================================================*/

int prep_get_row_bounds(MIPdesc *mip, int r_ind, double etol)
{    

   //COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;

   int j, c_ind, *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval; 
   //double * obj = mip->obj;

   double *ub = mip->ub;
   double *lb = mip->lb;

   double a_val;
  
   rows[r_ind].ub = rows[r_ind].lb = 0.0;
   //   rows[r_ind].fixed_obj_offset =  rows[r_ind].fixed_lhs_offset = 0.0;
   for(j = r_matbeg[r_ind]; j < r_matbeg[r_ind + 1]; j++){
      a_val = r_matval[j];
      c_ind = r_matind[j];
      if(a_val > etol){ 
	 if(rows[r_ind].ub < INF){
	    if(ub[c_ind] >= INF){
	       rows[r_ind].ub = INF;
	    }else{
	       rows[r_ind].ub += a_val * ub[c_ind];
	    }
	 }
	 if(rows[r_ind].lb > -INF){
	    if(lb[c_ind] <= -INF){
	       rows[r_ind].lb = -INF;
	    }else{
	       rows[r_ind].lb += a_val * lb[c_ind];
	    }
	 }
      }else if(a_val < -etol){
	 if(rows[r_ind].ub < INF){
	    if(lb[c_ind] <= -INF){
	       rows[r_ind].ub = INF;
	    }else{
	       rows[r_ind].ub += a_val * lb[c_ind];
	    }
	 }
	 if(rows[r_ind].lb > -INF){
	    if(ub[c_ind] >= INF){
	       rows[r_ind].lb = -INF;
	    }else{
	       rows[r_ind].lb += a_val * ub[c_ind];
	    }
	 }
      }
      
      // if(cols[c_ind].var_type == 'F'){
      //	 rows[r_ind].fixed_obj_offset = obj[c_ind]*ub[c_ind];
      //	 rows[r_ind].fixed_lhs_offset = a_val * ub[c_ind];
      // }
   }

   return 0;

}

/*===========================================================================*/
/*===========================================================================*/

double prep_rnd_integral(double val, double etol, char rnd_type)
{

   double new_bound = 0.0;
   
   if(rnd_type == RND_FLOOR){
      new_bound = ceil(val);
      if(val < new_bound - etol){
	 new_bound = floor(val);
      }
   }else{
      new_bound = floor(val);
      if(val > new_bound + etol){
	 new_bound = ceil(val);
      }	    
   }
   
   return new_bound;
}
/*===========================================================================*/
/*===========================================================================*/

int prep_solve_sr_rlx(PREPdesc *P, int row_cnt, int *row_indices)
{

   int i, j, k, l; 
   int termcode = SR_NO_UPDATES; 
   
   MIPdesc * mip = P->mip;
   prep_params params = P->params;
   MIPinfo *mip_inf = mip->mip_inf;

   COLinfo *cols = mip_inf->cols;
   ROWinfo *rows = mip_inf->rows;

   int n = mip->n, m = mip->m;
   int *c_matbeg = mip->matbeg;
   int *c_matind = mip->matind;   

   int * r_matbeg = mip->row_matbeg;
   int * r_matind = mip->row_matind;
   double * r_matval = mip->row_matval;
   double *rhs = mip->rhs;
   char *sense = mip->sense;

   double *ub = mip->ub;
   double *lb = mip->lb;

   //  SRrlx ** srows = (P->srows = (SRrlx **)malloc(m* sizeof(SRrlx*)));
   
   int max_sr_cnt, max_aggr_cnt, verbosity; //max_aggr_row_num, verbosity;
   char p_level, do_sr_rlx, do_aggr_row_rlx;  
   int obj_ind, tot_sub_pr;
   char can_iterate = TRUE;
   double etol;

   do_sr_rlx = params.do_single_row_rlx;      
   do_aggr_row_rlx = params.do_aggregate_row_rlx;
   p_level = params.level;
   etol = params.etol;
   verbosity = params.verbosity;

   /* get the max iteration limits */
   max_sr_cnt = params.max_sr_cnt;
   max_aggr_cnt = 0;
#if 0
   if(p_level > 2){
      if(do_sr_rlx){
	 max_sr_cnt = (int)(m*params.single_row_rlx_ratio);
	 if(max_sr_cnt > mip_inf->max_col_size - 1){
	    max_sr_cnt = mip_inf->max_col_size - 1;
	 }else if (max_sr_cnt < 1){
	    max_sr_cnt = 1;
	 }	
      }
      if(do_aggr_row_rlx){
	 /* how many rows to be aggregated */
	 max_aggr_row_num = (int)(m*params.max_aggr_row_ratio);
	 if(max_aggr_row_num > mip_inf->max_col_size - 1){
	    max_aggr_row_num = mip_inf->max_col_size - 1;
	 }else if (max_aggr_row_num < 2){
	    max_aggr_row_num = 2;
	 }
	 /* how many different aggr attempts for a single row */
	 max_aggr_cnt = max_sr_cnt; /* for now assume same 
				       with max_sr_cnt*/
      }
      
     }
#endif
   /* initialize arrays to be used for each subproblem*/

   
   SRdesc * sr, *d_sr;
   
   if(!(P->rows_checked)){
      P->rows_checked = (char *)malloc(m* CSIZE);
   }

   char *rows_checked = P->rows_checked; 
   double old_bound;  
   //char no_upper, no_lower;
   int row_ind;// const_row_ind; 
   //int updated_lb_cnt = 0, updated_ub_cnt = 0;
   int last_col_loc, last_row_loc; 

   /*max_sr_cnt should be up to some ratio!!!! may be too many to handle*/
   tot_sub_pr = max_sr_cnt + max_aggr_cnt; 

   /* do this only for all unbounded or all bounded rows for now...*/
   /* fixme... extend this if results happen to be good...*/

   for(i = 0; i < row_cnt; i++){

      obj_ind = row_indices[i];

      if(rows[obj_ind].bound_type == MIXED_BOUNDED_ROW || 
	 rows[obj_ind].is_redundant){
	 continue;
      }

      rows[obj_ind].orig_ub = rows[obj_ind].sr_ub = rows[obj_ind].ub;
      rows[obj_ind].orig_lb = rows[obj_ind].sr_lb = rows[obj_ind].lb;      
      
      if(verbosity >=4){
	 printf("init bounds: row: %i", i);
	 printf("\told_lb:");
	 if(rows[obj_ind].sr_lb > -INF){
	 printf("%f", rows[obj_ind].sr_lb);
	 }else{
	    printf("-inf");
	 }
	 printf("\told_ub:");
	 if(rows[obj_ind].sr_ub < INF){
	    printf("%f", rows[obj_ind].sr_ub);
	 }else{
	    printf("inf");
	 }
	 printf("\n");
      }
      

      //     srows[i] = (SRrlx *)calloc(tot_sub_pr, sizeof(SRrlx));
      memset(rows_checked, FALSE, CSIZE*m); 
      last_col_loc = r_matbeg[obj_ind];
      last_row_loc = c_matbeg[r_matind[last_col_loc]];

      for(j = 0; j < tot_sub_pr; j++){
	 /* first do single row stuff if can*/
	 //	 is_open_prob = 1;
	 row_ind = -1;
	 can_iterate = TRUE;	 
	 if(j < max_sr_cnt){
	    //int max_shared_row_ind = -1; 
	    //int max_shared_size = 0; 
	    //int row_search_iter_cnt = 0;	  
	    /* get something smarter here */
	    /*find a row that has the most common shared vars with this
	      one */	    
#if 0
	    for(k = last_col_loc; k < r_matbeg[obj_ind+1]; k++){
	       for(l = last_row_loc; l < c_matbeg[r_matind[k]+1]; 
		   l++){
		  if(!rows[c_matind[l]].is_redundant && 
		     !rows_checked[c_matind[l]]){
		     rows_checked[c_matind[l]] = TRUE;
		     if(rows[obj_ind].bound_type == rows[c_matind[l]].bound_type
			&& c_matind[l] != obj_ind) {	
			if(rows[c_matind[l]].size < max_shared_size){
			   row_ind = c_matind[l];			
			   max_shared_row_ind = row_ind;
			   row_search_iter_cnt++;
			}
			if(row_search_iter_cnt > 100){
			   break;
			}
		     }		     
		  }
	       }
	    }
#endif	       
	    
	    /*find a row to be used as a constraint*/	    
	    for(k = last_col_loc; k < r_matbeg[obj_ind+1]; k++){
	       for(l = last_row_loc; l < c_matbeg[r_matind[k]+1]; 
		   l++){
		  if(!rows[c_matind[l]].is_redundant && 
		     !rows_checked[c_matind[l]]){
		     rows_checked[c_matind[l]] = TRUE;
		     if(rows[obj_ind].bound_type == 
			rows[c_matind[l]].bound_type
			&& c_matind[l] != obj_ind) {	
			row_ind = c_matind[l];
			break;
		     }		     
		  }
	       }
	       
	       if(row_ind >= 0){
		  last_col_loc = k;
		  last_row_loc = l;
		  break;
	       }
	    }
	    
	    if(row_ind >= 0){
	       
	       sr_initialize(&(P->sr), n);	 
	       sr = P->sr;
#if 0
	       if(!sr){
		  sr = (SRdesc *)calloc(1, sizeof(SRdesc));
	       }	  

	       /* so the obj is row i, const is row row_ind */
	       /* gather the data for this problem */

	       sr->max_n = sr->min_n = 0;
	       sr->ub = sr->lb = 0.0;
	       sr->ub_offset = sr->lb_offset = 0.0;
	       sr->sum_a_max = sr->sum_a_min = 0.0;
	       sr->sum_c_max = sr->sum_c_min = 0.0;
	       sr->ub_updated = sr->lb_updated = FALSE;
#endif 	       
	       sr->prob_type = rows[obj_ind].bound_type;
	       sr->rhs = rhs[row_ind];
	       sr->sense = sense[row_ind];
	       
	       /* convert the problem to <= constraint*/
	       /* or solve it if it is unbounded */

	       switch(rows[obj_ind].bound_type){
		case OPEN_ROW: /* easiest case */
		   
		   sr->rhs_max = sr->rhs_min = sr->rhs;
		   
                   sr_solve_open_prob(P, sr, obj_ind, row_ind, r_matbeg, 
				      r_matind, r_matval, cols, ub, lb, etol);
		   
		   break;
		case ALL_BOUNDED_ROW:
		   if(rows[obj_ind].ub_inf_var_num + 
		      rows[obj_ind].lb_inf_var_num +
		      rows[obj_ind].free_var_num > 0 ||
		      rows[row_ind].ub_inf_var_num + 
		      rows[row_ind].lb_inf_var_num +
		      rows[row_ind].free_var_num > 0){
		      /* debug */
		      /* fixme, get rid of this bug*/
		      printf("something is wrong -case all_bounded_row-"
			     "prep_solve_sr_rlx(), exiting...\n");
		      return PREP_OTHER_ERROR;
		   }

		   /* always convert the problem into ax <= b ineq for max 
		      ax >= b for min */
		   
		   /* so,
 
		      _max arrays for solving [max cx st. ax <= b] 
		      _min arrays for solving [min cx st. ax >= b] 
		      
		      if sense = E;
		      
		      d_max arrays for solving [max cx st. -ax <= -b]
		      d_min arrays for solving [min cx st. -ax >= -b]

		   */

		   if(!sr->obj_max && rows[obj_ind].bound_type != OPEN_ROW){
		      sr_allocate(&sr, n);
		   }
#if 0
		   sr->obj_max = (double *)malloc(DSIZE*n);
		   sr->matval_max = (double *)malloc(DSIZE*n);
		   sr->matind_max = (int *)malloc(ISIZE*n);
		   sr->ratio_max = (double *)malloc(DSIZE*n);
		   
		   sr->obj_min = (double *)malloc(DSIZE*n);
		   sr->matval_min = (double *)malloc(DSIZE*n);
		   sr->matind_min = (int *)malloc(ISIZE*n);
		   sr->ratio_min = (double *)malloc(DSIZE*n);
		   

		   /* debug, get something smart instead of these */
		   sr->tmp_ind = (int *)malloc(ISIZE*n);
		   sr->fixed_ind = (int *)malloc(ISIZE*n);
		   
		   for(k = 0; k < n; k++){
		      sr->fixed_ind[k] = k;
		   }
#endif		   
		   switch(sr->sense){
		    case 'G': 
		       sr->rhs_max = -sr->rhs;
		       sr->rhs_min = sr->rhs;
		       break;
		    case 'L':
		       sr->rhs_max = sr->rhs;
		       sr->rhs_min = -sr->rhs;
		       break;
		    case 'E':
		       sr->rhs_max = sr->rhs;
		       sr->rhs_min = -sr->rhs;
		       
		       sr_initialize(&(P->d_sr), n);
		       d_sr = P->d_sr;
#if 0
		       if(!d_sr){
			  d_sr = (SRdesc *)calloc(1, sizeof(SRdesc));
		       }	  
		       
		       d_sr->max_n = d_sr->min_n = 0;
		       d_sr->ub = d_sr->lb = 0.0;
		       d_sr->ub_offset = d_sr->lb_offset = 0.0;
		       d_sr->sum_a_max = d_sr->sum_a_min = 0.0;
		       d_sr->sum_c_max = d_sr->sum_c_min = 0.0;
		       d_sr->ub_updated = d_sr->lb_updated = FALSE;
#endif
		       d_sr->prob_type = rows[obj_ind].bound_type;	       
		       d_sr->rhs = rhs[row_ind];
		       d_sr->sense = sense[row_ind];
		       
		       d_sr->rhs_max = -d_sr->rhs;
		       d_sr->rhs_min = d_sr->rhs;
		       
		       if(!d_sr->obj_max){
			  sr_allocate(&d_sr, n);
		       }
		       
#if 0
		        d_sr->obj_max = (double *)malloc(DSIZE*n);
			d_sr->matval_max = (double *)malloc(DSIZE*n);
			d_sr->matind_max = (int *)malloc(ISIZE*n);
			d_sr->ratio_max = (double *)malloc(DSIZE*n);	  
			
			d_sr->obj_min = (double *)malloc(DSIZE*n);
			d_sr->matval_min = (double *)malloc(DSIZE*n);
			d_sr->matind_min = (int *)malloc(ISIZE*n);
			d_sr->ratio_min = (double *)malloc(DSIZE*n);
			
			/* debug, get something smart instead of these */
			d_sr->tmp_ind = (int *)malloc(ISIZE*n);
			d_sr->fixed_ind = (int *)malloc(ISIZE*n);
			
			for(k = 0; k < n; k++){
			   d_sr->fixed_ind[k] = k;
			}
#endif
			break;
		   }
		   
	           sr_solve_bounded_prob(P, sr, d_sr, obj_ind, row_ind, 
					 r_matbeg, r_matind, r_matval, 
					 cols, ub, lb, etol);
		   if(!rows[obj_ind].is_redundant){
		      if(sr->sense == 'E'){
			 if(sr->ub > d_sr->ub){
			    sr->ub = d_sr->ub;
			 }
			 if(sr->lb < d_sr->lb){
			    sr->lb = d_sr->lb;
			 }		      
		      }
		      
		      sr->lb_updated = sr->ub_updated = TRUE;
		   }else{
		      break;
		   }
	       }
	       
	       /* check for any progress */
	       if(sr->lb_updated){
		  if(rows[obj_ind].sr_lb < sr->lb){
		     old_bound = rows[obj_ind].sr_lb;
		     rows[obj_ind].sr_lb = sr->lb;
		     /* debug */
		     if(termcode != SR_BOUNDS_UPDATED){
			termcode = SR_BOUNDS_UPDATED;
		     }
		     
		     if(verbosity >=5){
			printf("lb improved, " 
			       "row: %i \told_lb:%f \tnew_lb:%f\n", 
			       obj_ind, old_bound <= -INF ? 1 : old_bound, sr->lb);
		     }
		  }else if (rows[obj_ind].orig_lb > sr->lb + etol){
		     /* debug */
		     printf("error-lb, row: %i \told_lb:%f \tnew_lb:%f\n", 
			    obj_ind, rows[obj_ind].orig_lb, sr->lb);
		  }
	       }
	       if(sr->ub_updated){
		  if(rows[obj_ind].sr_ub > sr->ub){
		     old_bound = rows[obj_ind].sr_ub;
		     rows[obj_ind].sr_ub = sr->ub;
		     /* debug */
		     if(termcode != SR_BOUNDS_UPDATED){
			termcode = SR_BOUNDS_UPDATED;
		     }
		     if(verbosity >=5){
			printf("ub improved, " 
			       "row: %i \told_ub:%f \tnew_ub:%f\n", 
			       obj_ind, old_bound >= INF ? -1 : old_bound, sr->ub);
		     }
		  }else if(rows[obj_ind].orig_ub < sr->ub - etol){
		     /*debug*/
		     //		     if(verbosity >=5){
		     printf("error-ub, row: %i \told_ub:%f \tnew_ub:%f\n", 
			    obj_ind, rows[obj_ind].orig_ub, sr->ub);
			//		     }
		  }
		  if(sr->lb_updated){
		     if(sr->ub < sr->lb - etol){
			/* debug */
			printf("bounds err : " 
			       "row: %i \tnew_ub:%f \tnew_lb:%f\n", 
			       obj_ind, sr->ub, sr->lb);
			termcode = SR_INFEAS;
			break;
		     }
		  }
	       }
	       
	    }
	 }
      }

      /* debug */
      if(termcode == SR_INFEAS){
	 break;
      }

      if(verbosity >=4){
	 printf("finl bounds: row: %i", i);
	 printf("\tnew_lb:");
	 if(rows[obj_ind].sr_lb > -INF){
	    printf("%f", rows[obj_ind].sr_lb);
	 }else{
	    printf("-inf");
      }
	 printf("\tnew_ub:");
	 if(rows[obj_ind].sr_ub < INF){
	 printf("%f", rows[obj_ind].sr_ub);
	 }else{
	    printf("inf");
	 }
	 printf("\n\n");

      }
      //      printf("finl bounds: " 
      //	     "row: %i \tnew_lb:%f \tnew_ub:%f\n\n", 
      //	     i, rows[obj_ind].lb, rows[obj_ind].ub);       
   }
   
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
void sr_initialize(SRdesc **sr, int n){

   int do_clean = FALSE;

   if(!(*sr)){
      *sr = (SRdesc *)calloc(1, sizeof(SRdesc));
      do_clean = TRUE;
   }
   
   if(!do_clean){
      (*sr)->prob_type = 0;
      (*sr)->max_n = (*sr)->min_n = 0;
      (*sr)->ub = (*sr)->lb = 0.0;
      (*sr)->ub_offset = (*sr)->lb_offset = 0.0;
      (*sr)->sum_a_max = (*sr)->sum_a_min = 0.0;
      (*sr)->sum_c_max = (*sr)->sum_c_min = 0.0;
      (*sr)->ub_updated = (*sr)->lb_updated = FALSE;
      (*sr)->rhs = (*sr)->rhs_max = (*sr)->rhs_min = 0.0;
      (*sr)->sense = ' ';
      if((*sr)->obj_max){
	 memset((*sr)->reversed_max, FALSE, CSIZE*n);
	 memset((*sr)->reversed_min, FALSE, CSIZE*n);
	 memset((*sr)->var_stat_max, SR_VAR_IN, ISIZE*n);
	 memset((*sr)->var_stat_min, SR_VAR_IN, ISIZE*n);
      }
   }
}

/*===========================================================================*/
/*===========================================================================*/

void sr_allocate(SRdesc **sr, int n){
   
   int k;
   (*sr)->obj_max = (double *)malloc(DSIZE*n);
   (*sr)->matval_max = (double *)malloc(DSIZE*n);
   (*sr)->matind_max = (int *)malloc(ISIZE*n);
   (*sr)->ratio_max = (double *)malloc(DSIZE*n);
   (*sr)->reversed_max = (char *)malloc(CSIZE*n);
   
   (*sr)->obj_min = (double *)malloc(DSIZE*n);
   (*sr)->matval_min = (double *)malloc(DSIZE*n);
   (*sr)->matind_min = (int *)malloc(ISIZE*n);
   (*sr)->ratio_min = (double *)malloc(DSIZE*n);
   (*sr)->reversed_min = (char *)malloc(CSIZE*n);   
   
   /* for variable fixing, tightening etc... */

   (*sr)->var_max_opt = (double *)malloc(n* DSIZE);
   (*sr)->var_min_opt = (double *)malloc(n* DSIZE);
   (*sr)->var_stat_max = (int *)malloc(ISIZE*n);
   (*sr)->var_stat_min = (int *)malloc(n* ISIZE);
   (*sr)->var_obj_max = (double *)malloc(n* DSIZE);
   (*sr)->var_obj_min = (double *)malloc(n* DSIZE);
   (*sr)->var_matval_max = (double *)malloc(n* DSIZE);
   (*sr)->var_matval_min = (double *)malloc(n* DSIZE);   

   /* debug, get something smart instead of these */
   (*sr)->tmp_ind = (int *)malloc(ISIZE*n);
   (*sr)->fixed_ind = (int *)malloc(ISIZE*n);
   
   for(k = 0; k < n; k++){
      (*sr)->fixed_ind[k] = k;
   }
}
/*===========================================================================*/
/*===========================================================================*/

int sr_solve_bounded_prob(PREPdesc *P, SRdesc *sr, SRdesc *d_sr, 
			  int obj_ind, int row_ind, 
			  int *r_matbeg, int *r_matind, double *r_matval, 
			  COLinfo *cols, double *ub, double *lb, double etol)
{

   int k, l, col_ind;
   double c_val, a_val;

   for( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
      if(k < r_matbeg[obj_ind + 1] && 
	 (r_matind[k] < r_matind[l] ||
	  l >= r_matbeg[row_ind + 1])){
	 c_val = r_matval[k];
	 col_ind = r_matind[k];
	 sr_add_new_col(sr, d_sr, c_val, 0.0, col_ind,
			cols[col_ind].var_type, ub[col_ind], lb[col_ind], 
			sr->sense, 1, 1);
#if 0			 
	 if(c_val > 0.0){
	    ub_offset += ub[col_ind] * c_val;
	    lb_offset += lb[col_ind] * c_val;
	 }else if(r_matval[k] < 0.0){
	    lb_offset += ub[col_ind] * c_val;
	    ub_offset += lb[col_ind] * c_val;
	 }
#endif 
	 k++;
      }else if(l < r_matbeg[row_ind + 1] && 
	       (r_matind[k] > r_matind[l] ||
		k >= r_matbeg[obj_ind+1])){ 
	 a_val = r_matval[l];
	 col_ind = r_matind[l];
	 
	 sr_add_new_col(sr, d_sr, 0.0, a_val, col_ind,
			cols[col_ind].var_type, ub[col_ind], lb[col_ind],
			sr->sense, 0, 1);
	 l++;
      }else{
	 /* now the indices are equal, fill in the arrays */
	 c_val = r_matval[k];
	 a_val = r_matval[l];
	 col_ind = r_matind[k];
	 
	 if(c_val == 0.0 || a_val == 0.0){
	    printf("not nonzero???" 
		   "numerical issues -case bounded row-"
		   "sr_solve_bounded_prob(), exiting...\n");
	    return PREP_OTHER_ERROR; 
	 }
	 
	 sr_add_new_col(sr, d_sr, c_val, a_val, col_ind,
			cols[col_ind].var_type, ub[col_ind], lb[col_ind], 
			sr->sense, 2, 1);
	 k++;
	 l++;
      }
      if(k == r_matbeg[obj_ind + 1] && l == r_matbeg[row_ind + 1]){
	 break;
      }
   }  

   /* now solve the problem */
   if(!P->mip->mip_inf->rows[obj_ind].is_redundant){
      sr_find_opt_bounded(P, sr, obj_ind, ub, lb);      
   }

   if(!P->mip->mip_inf->rows[obj_ind].is_redundant){
      if(sr->sense == 'E'){
	 sr_find_opt_bounded(P, d_sr, obj_ind, ub, lb);
      }
   }
   

   int termcode = 0; 
   ROWinfo *rows = P->mip->mip_inf->rows;
   double min_ub = sr->ub;
   double max_lb = sr->lb;

   if(sr->sense == 'E'){
      if(!P->mip->mip_inf->rows[obj_ind].is_redundant){
	 if(min_ub > d_sr->ub){
	    min_ub = d_sr->ub;
	 }
	 
	 if(max_lb < d_sr->lb){
	    max_lb = d_sr->lb;
	 }
      }
   }
   if(rows[obj_ind].ub > min_ub || rows[obj_ind].lb < max_lb){
      termcode = prep_check_redundancy(P, obj_ind, TRUE, min_ub, max_lb,
				       FALSE, 0);
   }   
   
   if(prep_quit(termcode)){
      return termcode;
   }

#if 0
   if(termcode != PREP_MODIFIED){
      /* now see if we can improve the bounds */
      for(j= r_matbeg[obj_ind]; j < r_matbeg[obj_ind + 1]; j++){
	 col_ind = r_matind[j];
	 if(cols[col_ind].var_type == 'B'){
	    termcode = prep_improve_variable(P, col_ind, obj_ind, j, TRUE, 
					     impl_mode, TRUE, 
					     sr->var_max_opt[col_ind], 
					     sr->var_min_opt[col_ind]);
	    
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 }
      }

   }
#endif   
   
   return(0);

}

/*===========================================================================*/
/*===========================================================================*/

int sr_find_opt_bounded(PREPdesc *P, SRdesc *sr, int obj_ind, 
			 double *ub, double *lb)

{
   int i, last_ind, col_loc, col_ind, *var_stat; //,j, var_ind;
   //int var_fixed_ind, start_ind;
   char max_solved = FALSE, min_solved = FALSE;
   double lhs, ax, var_frac_val;
   /* get opt for each column (col fixed ub in min solved and 
      lb in max solved - check also a_vals)*/
   double bound; //*var_opt, var_offset, var_rhs, var_lhs;
   //char solve_for_fixed_var = TRUE; 


   //MIPdesc *mip = P->mip;

   //int *r_matbeg = mip->row_matbeg;
   //int *r_matind = mip->row_matind;
   //double *r_matval = mip->row_matval;


   //COLinfo *cols = mip->mip_inf->cols;
   //ROWinfo *rows = mip->mip_inf->rows;

   int * tmp_ind = sr->tmp_ind;
   double etol = P->params.etol;

   if(sr->sum_a_max < sr->rhs_max +etol || sr->max_n <= 0){
      sr->ub += sr->sum_c_max + sr->ub_offset;
      max_solved = TRUE;
   }
   
   if(sr->sum_a_min > sr->rhs_min - etol|| sr->min_n <= 0){
      sr->lb += sr->sum_c_min + sr->lb_offset;
      min_solved = TRUE;
   }

   if(max_solved && min_solved){
      /* check redundancy */
      //termcode = prep_check_redundancy(P, obj_ind, TRUE);
      
      /* redundant and useless row*/
      return PREP_UNMODIFIED;
   }

   if(!max_solved){ /* otherwise, the row is redundant and useless */

      var_stat = sr->var_stat_max;
      memcpy(tmp_ind, sr->fixed_ind, ISIZE*sr->max_n);
      qsort_di(sr->ratio_max, tmp_ind, sr->max_n);
      /* now fill in knapsack */
      lhs = 0;      
      for(i = sr->max_n - 1; i >=0; i--){
	 col_loc = tmp_ind[i];
	 col_ind = sr->matind_max[col_loc];
	 bound = ub[col_ind] - lb[col_ind];

	 if(lhs > sr->rhs_max - etol){
	    break;
	 }

	 ax = sr->matval_max[col_loc] * bound;

	 if(lhs + ax < sr->rhs_max - etol){
	    sr->ub += bound * 
	       sr->obj_max[col_loc];	    
	    lhs += ax;
	    var_stat[col_ind] = SR_VAR_IN_FIXED_UB;
	 }else{
	    var_frac_val = sr->obj_max[col_loc] *
	       (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    sr->ub += var_frac_val; //sr->obj_max[col_loc] * var_frac_val;
	    var_stat[col_ind] = SR_VAR_IN_FRAC;	    
	    last_ind = i;
	    break;
	 }
      }
      sr->ub += sr->ub_offset;

#if 0
      var_opt = sr->var_max_opt;
      //   var_stat = sr->var_stat_max;
      for(j= r_matbeg[obj_ind]; j < r_matbeg[obj_ind + 1]; j++){
	 
	 var_ind = r_matind[j];
	 /* for now just try binary variables */
	 if(cols[var_ind].var_type != 'B'){
	    continue;
	 }
	 
	 var_opt[var_ind] = sr->ub - sr->ub_offset;
	 var_offset = sr->ub_offset;
	 var_rhs = sr->rhs_max;
	 var_lhs = start_ind = 0;
	 solve_for_fixed_var = TRUE;	 	 
	 if(cols[var_ind].var_type == 'B'){
	    if(r_matval[j] > etol){
	       /* we want this to be fixed to its lower bound */
	       if(var_stat[var_ind] == SR_VAR_IN ||
		  var_stat[var_ind] == SR_VAR_FIXED_LB ||
		  var_stat[var_ind] == SR_VAR_IN_FIXED_LB){
		  /* same optimal for this var */
		  solve_for_fixed_var = FALSE;
	       }else if(var_stat[var_ind] == SR_VAR_IN_FRAC){
		  /* continue from where we left */
		  var_opt[var_ind] -= var_frac_val; 
		  start_ind = last_ind;
		  var_lhs = lhs;
		  solve_for_fixed_var = TRUE;
	       }else{
		  /* start from beginning */
		  if(var_stat[var_ind] == SR_VAR_FIXED_UB){
		     var_offset -= sr->var_obj_max[col_ind];
		     var_rhs += sr->var_matval_max[col_ind];
		  }
		  var_opt[var_ind] = 0.0;
		  solve_for_fixed_var = TRUE;
	       }
	    }else{
	       /* we want this to be fixed to its upper bound */
	       if(var_stat[var_ind] == SR_VAR_FIXED_UB ||
		  var_stat[var_ind] == SR_VAR_IN_FIXED_UB){
		  /* same optimal for this var */
		  solve_for_fixed_var = FALSE;
	       }else {
		  /* start from beginning */
		  var_offset += sr->var_obj_max[col_ind];
		  var_rhs -= sr->var_matval_max[col_ind];
		  var_opt[var_ind] = 0.0;
		  solve_for_fixed_var = TRUE;
	       }
	    }
	 }
	 
	 if(solve_for_fixed_var){
	    for(i =  start_ind; i >= 0; i--){
	       col_loc = tmp_ind[i];
	       col_ind = sr->matind_max[col_loc];
	       
	       if(col_ind != var_ind){
		  
		  if(var_lhs > var_rhs - etol){
		     break;
		  }
		  
		  bound = ub[col_ind] - lb[col_ind];		  
		  ax = sr->matval_max[col_loc] * bound;
		  
		  if(var_lhs + ax < var_rhs - etol){
		     var_opt[var_ind] += bound * 
			sr->obj_max[col_loc];	    
		     var_lhs += ax;
		  }else{
		     var_opt[var_ind] += sr->obj_max[col_loc] *
			(var_rhs - var_lhs)/sr->matval_max[col_loc];
		     break;
		  }
	       }
	    }
	 }
	 
	 var_opt[var_ind] += var_offset;	 	 
      }
#endif
   }

   if(!min_solved){ /* otherwise this row is redundant and useless */
      memcpy(tmp_ind, sr->fixed_ind, ISIZE*sr->min_n);
      qsort_di(sr->ratio_min, tmp_ind, sr->min_n);
      /* now fill in knapsack */
      lhs = 0;
      var_stat = sr->var_stat_min;
      for(i = 0; i < sr->min_n; i++){
	 col_loc = tmp_ind[i];
	 col_ind = sr->matind_min[col_loc];
	 bound = ub[col_ind] - lb[col_ind];
	 ax = sr->matval_min[col_loc] * bound;

	 if(lhs > sr->rhs_min - etol){
	    break;
	 }

	 if(lhs + ax < sr->rhs_min - etol){
	    sr->lb += bound * 
	       sr->obj_min[col_loc];	    
	    lhs += ax;
	    var_stat[col_ind] = SR_VAR_IN_FIXED_UB;
	 }else{
	    //	    ax = (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    sr->lb += sr->obj_min[col_loc] *
	       (sr->rhs_min - lhs)/sr->matval_min[col_loc];
	    var_stat[col_ind] = SR_VAR_IN_FIXED_UB;
	    last_ind = i;
	    break;
	 }
      }
      sr->lb += sr->lb_offset;
#if 0
      var_opt = sr->var_min_opt;
      //   var_stat = sr->var_stat_min; 
      
      for(j= r_matbeg[obj_ind]; j < r_matbeg[obj_ind + 1]; j++){
	 
	 var_ind = r_matind[j];
	 /* for now just try binary variables */
	 if(cols[var_ind].var_type != 'B'){
	    continue;
	 }
	 
	 var_opt[var_ind] = sr->lb - sr->lb_offset;
	 var_offset = sr->lb_offset;
	 var_rhs = sr->rhs_min;
	 var_lhs = start_ind = 0;
	 solve_for_fixed_var = TRUE;	 	 
	 if(cols[var_ind].var_type == 'B'){
	    if(r_matval[j] < etol){
	       /* we want this to be fixed to its lower bound */
	       if(var_stat[var_ind] == SR_VAR_IN ||
		  var_stat[var_ind] == SR_VAR_FIXED_LB ||
		  var_stat[var_ind] == SR_VAR_IN_FIXED_LB){
		  /* same optimal for this var */
		  solve_for_fixed_var = FALSE;
	       }else if(var_stat[var_ind] == SR_VAR_IN_FRAC){
		  /* continue from where we left */
		  var_opt[var_ind] -= var_frac_val; 
		  start_ind = last_ind;
		  var_lhs = lhs;
		  solve_for_fixed_var = TRUE;
	       }else{
		  /* start from beginning */
		  if(var_stat[var_ind] == SR_VAR_FIXED_UB){
		     var_offset -= sr->var_obj_min[col_ind];
		     var_rhs += sr->var_matval_min[col_ind];
		  }
		  var_opt[var_ind] = 0.0;
		  solve_for_fixed_var = TRUE;
	       }
	    }else{
	       /* we want this to be fixed to its upper bound */
	       if(var_stat[var_ind] == SR_VAR_FIXED_UB ||
		  var_stat[var_ind] == SR_VAR_IN_FIXED_UB){
		  /* same optimal for this var */
		  solve_for_fixed_var = FALSE;
	       }else {
		  /* start from beginning */
		  var_offset += sr->var_obj_min[col_ind];
		  var_rhs -= sr->var_matval_min[col_ind];
		  var_opt[var_ind] = 0.0;
		  solve_for_fixed_var = TRUE;
	       }
	    }
	 }
	 if(solve_for_fixed_var){
	    for(i = start_ind; i < sr->min_n; i++){
	       col_loc = tmp_ind[i];
	       col_ind = sr->matind_min[col_loc];
	       if(col_ind != var_ind){
		  
		  if(var_lhs > var_rhs - etol){
		     break;
		  }
		  
		  bound = ub[col_ind] - lb[col_ind];
		  ax = sr->matval_min[col_loc] * bound;
		  
		  if(var_lhs + ax < var_rhs - etol){
		     var_opt[var_ind] += bound * 
			sr->obj_min[col_loc];	    
		     var_lhs += ax;
		  }else{
		     //	    ax = (sr->rhs_max - lhs)/sr->matval_max[col_loc];
		     var_opt[var_ind] += sr->obj_min[col_loc] *
			(var_rhs - var_lhs)/sr->matval_min[col_loc];
		     break;
		  }
	       }
	    }
	    var_opt[var_ind] += var_offset;
	 }    
      }
#endif
   }


   /* first check redundancy */
#if 0
   int termcode; 


   termcode = prep_check_redundancy(P, obj_ind, TRUE, sr->ub, sr->lb, FALSE);

   if(prep_quit(termcode)){
      return termcode;
   }

   if(termcode != PREP_MODIFIED){
      /* now see if we can improve the bounds */
      for(j= r_matbeg[obj_ind]; j < r_matbeg[obj_ind + 1]; j++){
	 col_ind = r_matind[j];
	 if(cols[col_ind].var_type == 'B'){
	    termcode = prep_improve_variable(P, col_ind, obj_ind, j, TRUE, 
					     impl_mode, TRUE, 
					     sr->var_max_opt[col_ind], 
					     sr->var_min_opt[col_ind]);	    
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 }
      }

   }
#endif

   return 0;

}

/*===========================================================================*/
/*===========================================================================*/

/* will add the column to problem if necessary */

int sr_add_new_col(SRdesc *sr, SRdesc *d_sr, double c_val, double a_val, 
		   int col_ind, char var_type, double col_ub, 
		   double col_lb, char sense, 
		   int col_type, int col_bound_type)
{
   /* col_type = 
      0 => c_val = 0, a_val != 0
      1 => c_val != 0, a_val = 0
      2 => c_val != 0, a_val != 0
   */

   /* col_bound_type = 
      0 => open row
      1 => all bounded row
      2 => mixed bounded row
   */
   
   double rhs_ub_offset = a_val * col_ub;
   double rhs_lb_offset = a_val * col_lb;
   
   double obj_ub_offset = c_val * col_ub;
   double obj_lb_offset = c_val * col_lb;	
 
   if(col_bound_type == 1){
      if(col_type >= 0){
	 if(var_type != 'F'){
	    switch(sense){
	     case 'L':
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX, var_type);
		add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,	
				    col_ub, col_lb, 
				    SR_MIN, var_type);
		break;
	     case 'G':
		add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX, col_type);
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MIN, var_type);
		break;
	     case 'E':
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX, var_type);
		add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MIN, var_type);
		add_new_bounded_col(d_sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX, var_type);
		add_new_bounded_col(d_sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MIN, var_type);
		break;
	    }
	 }else{
	    
	    sr->ub_offset += obj_ub_offset;
	    sr->lb_offset += obj_ub_offset;
	    sr->rhs_max -= rhs_ub_offset;
	    sr->rhs_min -= rhs_ub_offset;

	    if(sense == 'E'){
	       d_sr->ub_offset += obj_ub_offset;
	       d_sr->lb_offset += obj_ub_offset;
	       d_sr->rhs_max -= rhs_ub_offset;
	       d_sr->rhs_min -= rhs_ub_offset;
	    }
	 }
      }
   }
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/

/* will add the column to problem if necessary, otherwise will update the 
   offset values. 
   will assume the sense is L for max and G for min, so 
   the a_val has to be sent after being updated. 
   For E, this function will be called twice each for max and min*/

int add_new_bounded_col(SRdesc *sr, double c_val, double a_val, int col_ind, 
			double rhs_ub_offset, double rhs_lb_offset, 
			double obj_ub_offset, double obj_lb_offset,
			double col_ub, double col_lb, int obj_sense, 
			char var_type)
{ 
   /* 
      ratio_type = 
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */

   //  int n;// = sr->max_n, min_n = sr->min_n; 		
 
   /* we will convert the vars so that u-l >= x >= 0 */


   int ratio_type = 0; 

   if(c_val > 0.0){
      if(a_val <= 0.0){
	 ratio_type = 1;
      }
   }else if(c_val < 0.0){
      if(a_val >= 0.0){
	 ratio_type = 2;
      }else{
	 ratio_type = 3;
      }
   }else{
      if(a_val <= 0.0){
	 ratio_type = 1;
      }else{
	 ratio_type = 2;
      }
   }
   
   int *n, *matind, *var_stat; 
   double *obj, *matval, *rhs, *obj_offset, *sum, *obj_sum, *ratios; 
   double *var_matval, *var_obj;
   char *is_reversed;
   if(obj_sense == SR_MAX){
      n = &(sr->max_n);
      obj_offset = &(sr->ub_offset);
      sum = &(sr->sum_a_max);
      obj_sum = &(sr->sum_c_max);
      rhs = &(sr->rhs_max);
      obj = sr->obj_max;
      matind = sr->matind_max;
      matval = sr->matval_max;
      ratios = sr->ratio_max;
      is_reversed = sr->reversed_max;
      var_stat = sr->var_stat_max;
      var_matval = sr->var_matval_max;
      var_obj = sr->var_obj_max;
      //var_lhs_offset = sr->opt_ub_var_offset;
   }else{
      n = &(sr->min_n);
      obj_offset = &(sr->lb_offset);
      sum = &(sr->sum_a_min);
      obj_sum = &(sr->sum_c_min);
      rhs = &(sr->rhs_min);
      obj = sr->obj_min;
      matind = sr->matind_min;
      matval = sr->matval_min;
      ratios = sr->ratio_min;
      is_reversed = sr->reversed_min;
      var_stat = sr->var_stat_min;
      var_matval = sr->var_matval_min;
      var_obj = sr->var_obj_min;
   }

#if 0   
   sr->var_ub_lhs_offset[col_ind] = -rhs_ub_offset;
   sr->var_lb_lhs_offset[col_ind] = -rhs_lb_offset;
   sr->var_lb_obj_offset[col_ind] = obj_lb_offset;
   sr->var_ub_obj_offset[col_ind] = obj_ub_offset;
#endif

   if(ratio_type == 0){
      obj[*n] = c_val;
      matval[*n] = a_val;
      matind[*n] = col_ind;
      ratios[*n] = c_val/a_val;
      if(obj_sense == SR_MAX){
	 *sum += (rhs_ub_offset - rhs_lb_offset);
	 *obj_sum += (obj_ub_offset - obj_ub_offset);
      }else{
	 /* since bounds are converted to be 
	    u - l >= x >= 0 */
	 *sum += 0.0;//rhs_lb_offset;   
	 *obj_sum += 0.0;//obj_lb_offset;
      }
      (*n)++;
      /* to solve bounds problem */      
      *rhs += -(rhs_lb_offset); /* conversion by x = y + l */
      *obj_offset += obj_lb_offset;

   }else if((ratio_type == 1 && obj_sense == SR_MAX) ||
	    (ratio_type == 2 && obj_sense == SR_MIN)){
      *rhs += -rhs_ub_offset;
      *obj_offset += obj_ub_offset;
      var_stat[col_ind] = SR_VAR_FIXED_UB;
      var_matval[col_ind] = a_val;
      var_obj[col_ind] = c_val;
   }else if((ratio_type == 1 && obj_sense == SR_MIN) ||
	    (ratio_type == 2 && obj_sense == SR_MAX)){
      *rhs += -rhs_lb_offset;
      *obj_offset += obj_lb_offset;
      var_stat[col_ind] = SR_VAR_FIXED_LB;
      var_matval[col_ind] = a_val;
      var_obj[col_ind] = c_val;
   }else{
      obj[*n] = -c_val;
      matval[*n] = -a_val;
      matind[*n] = col_ind;
      ratios[*n] = c_val/a_val;
      is_reversed[*n] = TRUE;
      if(obj_sense == SR_MAX){
	 *sum += -rhs_ub_offset + +rhs_lb_offset;
	 *obj_sum += -obj_ub_offset + rhs_lb_offset;
      }else{
	 *sum += 0.0; //-rhs_lb_offset;
	 *obj_sum += 0.0; //-obj_lb_offset;
      }
      (*n)++;
      /* to solve bounds problem */      
      *rhs += -(rhs_ub_offset); /* conversion by x = -y + u */
      *obj_offset += obj_ub_offset;
   }

#if 0
   if(obj_sense == SR_MAX){
      n = sr->max_n;
      /* boundlar da degisiyor!!!! shoot!*/
      /* always u-l??? */

      if(ratio_type == 0){
	 sr->obj_max[n] = obj_val;
	 sr->matval_max[n] = a_val;
	 sr->matind_max[n] = col_ind;
	 n++;
	 	 
      }else if(ratio_type == 1){

	 sr->rhs_max += -rhs_ub_offset;
	 sr->ub_offset += obj_ub_offset;
	 
      }else if(ratio_type == 2){
	 sr->rhs_max += -rhs_lb_offset;
	 sr->ub_offset += obj_lb_offset;

      }else{

	 sr->obj_max[n] = -obj_val;
	 sr->matval_max[n] = -a_val;
	 sr->matind_max[n] = col_ind;
	 //	 sr->bounds_changed[n] = TRUE;
	 n++;
	 
	 sr->rhs_max += -(rhs_ub_offset + rhs_lb_offset);
	 sr->ub_offset += obj_ub_offset + obj_lb_offset;	 
      }
      sr->max_n = n;

   }else if(sense == SR_MIN){
      n = sr->min_n;

      if(ratio_type == 0){

	 sr->obj_min[n] = obj_val;
	 sr->matval_min[n] = a_val;
	 sr->matind_min[n] = col_ind;
	 n++;
	 
      }else if(ratio_type == 1){

	 sr->rhs_min += -rhs_lb_offset;
	 sr->lb_offset += obj_lb_offset;

      }else if(ratio_type == 2){
	 
	 sr->rhs_min += -rhs_ub_offset;
	 sr->lb_offset += obj_ub_offset;
	 
      }else{
		
	 sr->obj_min[n] = -obj_val;
	 sr-.matval_min[n] = -a_val;
	 sr->matind_min[n] = col_ind;
	 //	 sr->bounds_changed_min[n] = TRUE;
	 n++;
	 
	 sr->rhs_min += -(rhs_ub_offset + rhs_lb_offset) ;
	 sr->lb_offset += obj_ub_offset + obj_lb_offset; /* this should solve 
							    the bounds 
							    problem */
      }
      sr->min_n = n;
   }
#endif
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/

/* will modify the constraint to E and solve it */

int sr_solve_open_prob(PREPdesc *P, SRdesc *sr, int obj_ind, 
		       int row_ind, int *r_matbeg, 
		       int *r_matind, double *r_matval, COLinfo *cols, 
		       double *ub, double *lb, double etol)
{

   int l, k, col_ind;

   double max_dual_ub = INF, min_dual_ub = INF;
   double max_dual_lb = -INF, min_dual_lb = -INF;
   double d_ratio, obj_val, a_val; 

   char no_upper = FALSE, no_lower = FALSE, is_free_column;
   char is_fixed_column = FALSE;
   char can_iterate = TRUE, prob_infeasible = FALSE, is_null_obj;

   double *ub_offset = &(sr->ub_offset);
   double *lb_offset = &(sr->lb_offset);
   double rhs = sr->rhs;
   char sense = sr->sense;   

   double obj_ub_offset;
   double obj_lb_offset;


   //  sr->prob_type = OPEN_PROB;

   for( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
      if(k < r_matbeg[obj_ind + 1] && 
	 (r_matind[k] < r_matind[l] ||
	  l >= r_matbeg[row_ind + 1])){
	 if(r_matval[k] > 0.0){	    
	    if(!no_upper){
	       if(ub[r_matind[k]] >= INF){
		  no_upper = TRUE;
	       }else{
		  *ub_offset += ub[r_matind[k]] * r_matval[k];
	       }
	    }
	    if(!no_lower){
	       if(lb[r_matind[k]] <= -INF){
		  no_lower = TRUE;
	       }else{
		  *lb_offset += lb[r_matind[k]] * r_matval[k];
	       }
	    }
	 }else if (r_matval[k] < 0.0){
	    if(!no_lower){
	       if(ub[r_matind[k]] >= INF){
		  no_lower = TRUE;
	       }else{
		  *lb_offset += ub[r_matind[k]] * r_matval[k];
	       }
	    }
	    if(!no_upper){
	       if(lb[r_matind[k]] <= -INF){
		  no_upper = TRUE;
	       }else{
		  *ub_offset += lb[r_matind[k]] * r_matval[k];
	       }
	    }
	 }
	 k++;
/*       }else if(l < r_matbeg[row_ind + 1] &&  */
/* 	       (r_matind[k] > r_matind[l] || */
/* 		k >= r_matbeg[obj_ind+1])){  */
	 
/* 	 if(ub[col_ind] < INF && lb[col_ind] > -INF){ */
/* 	    /\* debug - get vars.type here *\/ */
/* 	    /\* we check this very same thing down too, fix this! *\/ */
/* 	    if(ub[col_ind] > lb[col_ind] + etol){ */
/* 	       /\* debug *\/ */
/* 	       printf("bounded column-first case -case all open row-" */
/* 		      "sr_solve_open_prob(), exiting...\n"); */
/* 	       return PREP_OTHER_ERROR;  */
/* 	    }else{ */
/* 	       a_val = r_matval[l]; */
/* 	       col_ind = r_matind[l]; */
/* 	       /\* fix column, take care of it here *\/ */
/* 	       rhs += -(a_val * lb[col_ind]); */
/* 	    } */
/* 	 }else{ */
/* 	    if(r_matval[l] > 0.0){ */
/* 	       if(min_dual_ub > 0.0){ */
/* 		  min_dual_ub = 0.0; */
/* 	       } */
/* 	       if(max_dual_ub > 0.0){ */
/* 		  max_dual_ub = 0.0; */
/* 	       } */
/* 	    }else if(r_matval[l] < 0.0){ */
/* 	       if(min_dual_lb < 0.0){ */
/* 		  min_dual_lb = 0.0; */
/* 	       } */
/* 	       if(max_dual_lb < 0.0){ */
/* 		  max_dual_lb = 0.0; */
/* 	       } */
/* 	    } */
/* 	 } */
/* 	 l++; */
      }else{
	 if(l < r_matbeg[row_ind + 1] && 
	    (r_matind[k] > r_matind[l] ||
	     k >= r_matbeg[obj_ind+1])) {
	    is_null_obj = TRUE;
	    obj_val = 0.0;
	 }else{
	    is_null_obj = FALSE;
	    obj_val = r_matval[k];
	 }

	 /* first convert the column into 0 <= x <= inf */

	 a_val = r_matval[l];
	 col_ind = r_matind[l];
	 is_free_column = FALSE;
	 is_fixed_column = FALSE;
	 if(ub[col_ind] < INF && lb[col_ind] > -INF){
	    /* debug - get vars.type here */
	    if(ub[col_ind] > lb[col_ind] + etol){
		  /* debug */
	       printf("bounded column -case all open row-"
		      "sr_solve_open_prob(), exiting...\n");
	       return PREP_OTHER_ERROR; 
	    }else{
	       /* fix column, take care of it here */
	       if(!is_null_obj){
		  obj_lb_offset = obj_val * lb[col_ind]; 
		  if(!no_upper){
		     *ub_offset += obj_lb_offset;		  
		  }
		  if(!no_lower){
		     *lb_offset += obj_lb_offset;
		  }
	       }
	       rhs += -(a_val *lb[col_ind]);
	       is_fixed_column = TRUE;
	    }
	 }else if(ub[col_ind] >= INF){
	    if(lb[col_ind] > -INF){
	       if(!is_null_obj){
		  obj_lb_offset = obj_val * lb[col_ind];     
		  if(!no_upper){
		     *ub_offset += obj_lb_offset;
		  }
		  if(!no_lower){
		     *lb_offset += obj_lb_offset;
		  }
	       }
	       rhs += -(a_val * lb[col_ind]);
	    }else{
	       is_free_column = TRUE;
	    }
	 }else{
	    if(!is_null_obj){
	       obj_ub_offset = obj_val * ub[col_ind];
	       if(!no_upper){
		  *ub_offset += obj_ub_offset;
	       }
	       if(!no_lower){
		  *lb_offset += obj_ub_offset;
	       }
	    }
	    rhs += -(a_val * ub[col_ind]);
	    obj_val = -obj_val;
	    a_val = -a_val; 
	 }
	 
	 /* now get dual bounds */
	 if(!is_fixed_column){
	    if(a_val > 0.0 || a_val < 0.0){
	       d_ratio = obj_val/a_val;
	       if(a_val > 0.0){
		  if(d_ratio < min_dual_ub){
		     min_dual_ub = d_ratio;
		  }
		  if(-d_ratio < max_dual_ub){
		     max_dual_ub = -d_ratio;
		  }
		  if(is_free_column){
		     if(d_ratio > min_dual_lb){
			min_dual_lb = d_ratio;
		     }
		     if(-d_ratio > max_dual_lb){
			max_dual_lb = -d_ratio;
		     }
		  }
	       }else{
		  if(d_ratio > min_dual_lb){
		     min_dual_lb = d_ratio;
		  }
		  if(-d_ratio > max_dual_lb){
		     max_dual_lb = -d_ratio;
		  }
		  
		  if(is_free_column){
		     if(d_ratio < min_dual_ub){
			min_dual_ub = d_ratio;
		     }
		     if(-d_ratio < max_dual_ub){
			max_dual_ub = -d_ratio;
		     }
		  }
	       }
	       
	       if(min_dual_lb > min_dual_ub){
		  no_lower = TRUE;
	       }
	       if(max_dual_lb > max_dual_ub){
		  no_upper = TRUE;
		  /* debug */
		  //printf("unbounded or infeasible problem?" 
		  // "-case all open row-"
		  // "prep_solve_sr_rlx(), exiting...\n");
		  //return PREP_OTHER_ERROR; 
	       }
	    }else{
	       /* debug */
	       printf("not nonzero???" 
		      "numerical issues -case all open row-"
		      "prep_solve_sr_rlx(), exiting...\n");
	       return PREP_OTHER_ERROR; 
	    }
	 }
	 l++;
	 if(!is_null_obj){
	    k++;
	 }
      }
      
      if((no_upper && no_lower)){
	 can_iterate = FALSE;
	 break;
      }
      if(k == r_matbeg[obj_ind + 1] && l == r_matbeg[row_ind + 1]){
	 break;
      }
   }
   
   if(can_iterate){
      /* update the bounds for this row */
      
      switch(sense){
       case 'L': 
	  if(max_dual_ub > 0.0){
	     max_dual_ub = 0.0;
	  }
	  if(min_dual_ub > 0.0){
	     min_dual_ub = 0.0;
	  }
	  break;
       case 'G':
	  if(max_dual_lb < 0.0){
	     max_dual_lb = 0.0;
	  }
	  if(min_dual_lb < 0.0){
	     min_dual_lb = 0.0;
	  }
	  break;
      }   
      
      /* check again */
      //   if(min_dual_lb > min_dual_ub ||/
      //	 max_dual_lb > max_dual_ub){/
      //	 printf("unbounded or infeasible problem?" 
      //		"-case all open row-"
      //		"prep_solve_sr_rlx(), exiting...\n");
      //	 return PREP_OTHER_ERROR; 
      //  }
      
      if(!no_lower){
	 if(rhs >= 0){
	    if(min_dual_ub < INF){
	       sr->lb = min_dual_ub * rhs;
	    }else{
	       prob_infeasible = TRUE;
	    }
	 }else{
	    if(min_dual_lb > -INF){
	       sr->lb = min_dual_lb *rhs;
	    }else{
	       prob_infeasible = TRUE;
	    }
	 }
	 if(!prob_infeasible){
	    sr->lb += *lb_offset; 
	    sr->lb_updated = TRUE;
	 }
#if 0
	 if(rows[obj_ind].lb < sub_lb){
	    rows[obj_ind].lb = sub_lb;
	    updated_lb_cnt++;
	 }
#endif
      }

      if(!prob_infeasible){
	 if(!no_upper){
	    if(rhs >= 0){
	       if(max_dual_ub < INF){
		  sr->ub = -(max_dual_ub * rhs);
	       }else{
		  prob_infeasible = TRUE;
	       }
	    }else{
	       if(max_dual_lb > -INF){
		  sr->ub = -(max_dual_lb *rhs);
	       }else{
		  prob_infeasible = TRUE;
	       }
	    }
	    if(!prob_infeasible){
	       sr->ub += *ub_offset;
	       sr->ub_updated = TRUE;
	    }
#if 0
	    if(rows[obj_ind].ub > sub_ub){
	       rows[obj_ind].ub = sub_ub;
	       updated_ub_cnt++;
	    }
#endif
	 }
      }
   }
   
   //   return(sr->lb_updated || sr->ub_updated);
   return(prob_infeasible);
}

/*===========================================================================*/
/*===========================================================================*/

//int prep_check_redundancy(PREPdesc *P, int row_cnt, int *row_start, 
//			  char use_sr_bounds, 
//			  double sr_ub, double sr_lb, char impl_mode,
//			  int dive_level)
int prep_check_redundancy(PREPdesc *P, int row_ind, 
			  char use_sr_bounds, 
			  double sr_ub, double sr_lb, char impl_mode,
			  int dive_level)
{
   int i, termcode = PREP_UNMODIFIED, fix_type = FIX_NO_BOUND;
   int fixed_row = FALSE, fix_all_lb = FALSE, fix_all_ub = FALSE, col_ind;
   int debug_cnt = 0; 
   double a_val, ub, lb, new_bound = 0.0, rnd_floor, rnd_ceil;

   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols; 

   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double *c_ub = mip->ub;
   double *c_lb = mip->lb;
   double *obj = mip->obj;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind; 
   double *r_matval = mip->row_matval;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   prep_stats *stats = &(P->stats);
   //int dive_level = 1;

  /* 
      ratio_type = 
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */
   //   if(row.fixed_var_num + row.fixable_var_num >= row.size){

   
   //   for(i = 0; i < row_cnt; i++){
   //      row_ind = row_start[i];

      if(!use_sr_bounds && rows[row_ind].fixed_var_num  >= rows[row_ind].size){
	 if((sense == 'L' && rows[row_ind].fixed_lhs_offset > rhs + etol) ||
	    (sense == 'E' && 
	     !prep_is_equal(rows[row_ind].fixed_lhs_offset, rhs, etol))){
	    stats->row_infeas_ind = row_ind;
	    return PREP_INFEAS;
	 }									       
	 //rows[row_ind].is_redundant = TRUE;
	 termcode = PREP_MODIFIED;
      }else if(sense != 'R' && !use_sr_bounds && 
	       rows[row_ind].fixed_var_num >= rows[row_ind].size - 1){
	 for(i = r_matbeg[row_ind]; i < r_matbeg[row_ind + 1]; i++){
	    col_ind = r_matind[i];
	    a_val = r_matval[i];
	    if(cols[col_ind].var_type != 'F'){	    
	       debug_cnt++;
	       if(a_val > etol || a_val < -etol){
		  new_bound = (double)
		     ((rhs - rows[row_ind].fixed_lhs_offset)/a_val);	 
		  
		  if(sense == 'E'){
		     if(new_bound > c_ub[col_ind] + etol || 
			new_bound < c_lb[col_ind] - etol){
			stats->col_infeas_ind = col_ind;
			stats->row_infeas_ind = row_ind;
			return PREP_INFEAS;
		     }
		     if(cols[col_ind].var_type != 'C'){
			rnd_floor = floor(new_bound);
			rnd_ceil = ceil(new_bound);
			if(new_bound >= rnd_floor + etol &&
			   new_bound <= rnd_ceil - etol){
			   stats->col_infeas_ind = col_ind;
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}else{
			   if(new_bound < rnd_floor + etol){
			      new_bound = rnd_floor;
			   }else{
			      new_bound = rnd_ceil;
			   }
			}				  
		     }
		     fix_type = FIX_OTHER;
		  }else{
		     if(cols[col_ind].col_size > 1){
			if((sense == 'G' && a_val > etol) ||
			   (sense == 'L' && a_val < -etol)){
			   if(new_bound > c_ub[col_ind] + etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if(new_bound > c_lb[col_ind] + etol ){
			      if(cols[col_ind].var_type != 'C'){
				 new_bound = prep_rnd_integral(new_bound, etol, 
							       RND_CEIL);
			      }
			      fix_type = IMPROVE_LB;
			   }
			}else{
			   if(new_bound < c_lb[col_ind] - etol){
			      stats->col_infeas_ind = col_ind; 
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if(new_bound < c_ub[col_ind] - etol){
			      if(cols[col_ind].var_type != 'C'){
				 new_bound= prep_rnd_integral(new_bound, etol, 
							      RND_FLOOR);
			      }
			      fix_type = IMPROVE_UB;
			   }
			}
		     }else{
			/*so we can fix this column and delete this row*/
			/* and sense == 'L' */
			if(a_val > etol){
			   if(new_bound < c_lb[col_ind] - etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if(obj[col_ind] >= 0.0){
			      new_bound = c_lb[col_ind];
			   }else{
			      if(new_bound > c_ub[col_ind] + etol){
				 new_bound = c_ub[col_ind];
			      }else{
				 if(cols[col_ind].var_type != 'C'){
				    new_bound= prep_rnd_integral(new_bound, etol, 
								 RND_FLOOR);
				 }
			      }
			   }
			}else if(a_val < -etol){
			   if(new_bound > c_ub[col_ind] + etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if(obj[col_ind] <= 0.0){
			      new_bound = c_ub[col_ind];
			   }else{
			      if(new_bound <  c_lb[col_ind] - etol){
				 new_bound = c_lb[col_ind];
			      }else{
				 if(cols[col_ind].var_type != 'C'){
				    new_bound= prep_rnd_integral(new_bound, etol, 
								 RND_CEIL);
				 }
			      }
			   }
			}
			
			fix_type = FIX_OTHER;
		     }
		  }
	       }else{
		  /* so a_val has been fixed to 0.0 */
		  /* this row is redundant */
		  /* check if we can fix this column*/
		  if(cols[col_ind].col_size == 1){
		     if(sense == 'E'){
			if(!prep_is_equal(rows[row_ind].fixed_lhs_offset, rhs, 
					  etol)){
			   stats->col_infeas_ind = col_ind;
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
		     }else{ /* sense is 'L' */
			if(rows[row_ind].fixed_lhs_offset > rhs + etol){
			   stats->col_infeas_ind = col_ind;
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
			
			if(obj[col_ind] >= 0.0){
			   new_bound = c_lb[col_ind];
			}else{
			   new_bound = c_ub[col_ind];
			}
		     }
		     fix_type = FIX_OTHER;
		  }else{
		     /* just declare this row to be redundant 
			here*/
		     //cols[col_ind].col_size--;
		     termcode = PREP_MODIFIED;
		  }
	       }
	       //c_ub[col_ind] = c_lb[col_ind] = new_bound;
	       //cols[col_ind].var_type = 'F';	       
	       
	       //cols[col_ind].fix_row_ind = row_ind;
	       //rows[row_ind].fixed_var_num++;
	       //if(cols[col_ind].var_type == 'B'){
	       //   rows[row_ind].bin_var_num--;
	       // }
	       
	       //rows[row_ind].is_redundant = TRUE;
	       if(fix_type != FIX_NO_BOUND){
		  if(cols[col_ind].var_type == 'B'){
		     fix_type = FIX_BINARY;
		  }
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind,
							    row_ind, 
							    dive_level,
							    new_bound, 
							    fix_type, TRUE,
							    impl_mode);
		  if(prep_quit(termcode)){
		     return termcode;		  
		  }else if(rows[row_ind].is_redundant){
		     return PREP_MODIFIED;
		  }
	       }
	       termcode = PREP_MODIFIED;
	       /* debug - need to break here */
	       break;
	    }
	 }
	 /* debug */
	 if(debug_cnt > 1){
	    printf("error in prep_check_redundancy()\n");
	    return PREP_OTHER_ERROR;
	 }else if(debug_cnt == 0){
	    /* end point of recursive procedure */
	    /*
	      probably many variables have been fixed during recursive 
	      procedure, so just declare redundancy here 
	      no need to do anything */	 
	    termcode = PREP_MODIFIED;
	 }
      }
      
      if(termcode == PREP_UNMODIFIED){
	 if(use_sr_bounds){
	    ub = sr_ub;
	    lb = sr_lb;
	 }else{
	    ub = rows[row_ind].ub;
	    lb = rows[row_ind].lb;
	 }
	 
	 /* 
	    if(rows[row_ind].type != CONTINUOUS_TYPE && 
	    rows[row_ind].type != BIN_CONT_TYPE && 
	    rows[row_ind].type != INT_CONT_TYPE && 
	    rows[row_ind].type != ALL_MIXED_TYPE && 	       
	    rows[row_ind].coef_type != FRACTIONAL_VEC){
	    lb = prep_rnd_integral(lb, etol, RND_CEIL);
	    ub = prep_rnd_integral(ub, etol, RND_FLOOR);
	    }
	 */
	 
	 /*check redundancy and infeasibility*/
	 if(lb > ub + etol){
	    stats->row_infeas_ind = row_ind;
	    /* debug, can this happen? */
	    return PREP_INFEAS;
	 }else if (lb > ub - etol){
	    fixed_row = TRUE;
	    fix_all_ub = TRUE;
	    /* check infeasibility here */
	    if(lb > rhs + etol){
	       stats->row_infeas_ind = row_ind;
	       return PREP_INFEAS;
	    }
	    if(sense == 'E'){
	       if(ub < rhs - etol){
		  stats->row_infeas_ind = row_ind;
		  return PREP_INFEAS;
	       }
	    }
	 }
	 
	 if(!fixed_row){
	    switch(sense){
	     case 'L': 
		if(lb > rhs + etol){
		   stats->row_infeas_ind = row_ind;
		   /* prob infeasible */
		   return PREP_INFEAS;
		}
		if(ub < rhs - etol){ 
		   //rows[row_ind].is_redundant = TRUE;
		   termcode = PREP_MODIFIED;
		}
		if(lb > rhs - etol){
		   fix_all_lb = TRUE;
		}
		break;
#if 0
	     case 'G':
		if(ub < rhs - etol){
		   stats->row_infeas_ind = row_ind;
		   /* prob infeasible */
		   return PREP_INFEAS;
		}
		if(lb > rhs + etol || fixed_row){
		   //rows[row_ind].is_redundant = TRUE;
		   termcode =  PREP_MODIFIED;
		}		
		break;
#endif
	     case 'E':
		if(lb > rhs + etol || 
		   ub < rhs - etol){
		   /* prob infeasible */
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
		if(prep_is_equal(lb, rhs, etol)){
		   fix_all_lb = TRUE;
		}
		if(prep_is_equal(ub, rhs, etol)){
		   fix_all_ub = TRUE;
		}
		break;
	     default:
		/* debug */		
		break;
	    }
	 }
      }
      
      if(fix_all_lb || fix_all_ub){
	 if(!use_sr_bounds){
	    
	    /* first temporarily fix the columns 
	       this is to detect infeasibilities during
	       recursive procedure */
	    
	    /* keep bounds for recursive infeasibility detection */
	    /* not very efficient */
#if 0
	    new_bounds = (double *)malloc(DSIZE * rows[row_ind].size)
	       
	       if(fix_all_lb){
		  for(k = 0, i = r_matbeg[row_ind]; i < r_matbeg[row_ind +1]; i++){
		     col_ind = r_matind[i];
		     a_val = r_matval[i];
		     if(cols[col_ind].var_type != 'F' && 
			(a_val > etol || a_val <-etol)){
			if(a_val > etol){
			   new_bound[k] = c_lb[col_ind];
			}else{
			   new_bound[k] = c_ub[col_ind];
			}
			k++;
		     }
		  }
	       }else{
		  for(k = 0, i = r_matbeg[row_ind]; i < r_matbeg[row_ind +1];i++){
		     col_ind = r_matind[i];
		     a_val = r_matval[i];
		     if(cols[col_ind].var_type != 'F' && 
			(a_val > etol || a_val <-etol)){
			if(a_val > etol){
			   new_bound [k]= c_ub[col_ind];
			}else{
			   new_bound[k] = c_lb[col_ind];
			}
			k++;
		     }
		  }
	       }
	    
	    /* now go ahead and imply them */
	    
	    for(i = r_matbeg[row_ind]; i < r_matbeg[row_ind +1]; i++){
	       col_ind = r_matind[i];
	       a_val = r_matval[i];
	       if(cols[col_ind].var_type != 'F' && 
		  (a_val > etol || a_val <-etol)){
	       }
	    }
	    
#endif 	       

	    rows[row_ind].is_redundant = TRUE;
	    
	    if(fix_all_lb){
	       fix_type = FIX_ROW_LB;
	    }else{
	       fix_type = FIX_ROW_UB;
	    }
	    
	    //fix_type = FIX_OTHER;
	    
	    /* now go ahead and imply them */
#if 0	 
	    for(i = r_matbeg[row_ind]; i < r_matbeg[row_ind +1]; i++){
	       col_ind = r_matind[i];
	       a_val = r_matval[i];
	       if(cols[col_ind].var_type != 'F' && 
		  (a_val > etol || a_val <-etol)){
		  if(a_val > etol){
		     if(fix_all_lb){
			new_bound = c_lb[col_ind];
		     }else{
			new_bound = c_ub[col_ind];
		     }
		  }else{ 
		     if(fix_all_lb){
			new_bound = c_ub[col_ind];
		     }else{
			new_bound = c_lb[col_ind];
		     }
		  }
	       }
	    }
#endif
	    termcode = prep_modified_cols_update_info(P, r_matbeg[row_ind + 1]  - 
						      r_matbeg[row_ind], 
						      &r_matind[r_matbeg[row_ind]],
						      row_ind, 
						      dive_level, 0.0, 
						      fix_type, TRUE, impl_mode);
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	    
	    termcode = PREP_MODIFIED;
	    
	 }else{
	    if(fix_all_lb && fix_all_ub){
	       printf("sr bounds are equal to rhs - row redundant!\n");
	       termcode =  PREP_MODIFIED;
	    }	 
	 }
      }
      
      if(termcode == PREP_MODIFIED){
	 stats->rows_deleted++;
	 if(verbosity >= 2){
	    prep_declare_redundant_row(rows[row_ind], row_ind, sense, 
				       rhs);
	 }
	 termcode = prep_deleted_row_update_info(mip, row_ind);
	 if(prep_quit(termcode)){
	    return termcode;
	 }else{
	    return PREP_MODIFIED;
	 }
	 
      }
      //   }

   return termcode;
}
 
/*===========================================================================*/
/*===========================================================================*/
/* make these more robust */
int prep_declare_redundant_row(ROWinfo row, int row_ind, char sense, 
			       double rhs)
{
   printf("row [%i] is redundant: ", row_ind); 
   printf("ub: ");
   if(row.ub < INF){
      printf("%f",row.ub);
   }else{
      printf("INF");
   }
   printf("\t lb: ");
   if(row.lb > -INF){
      printf("%f",row.lb);
   }else{
      printf("-INF");
   }
   printf("\t sense: %c \t rhs: %f\n", sense, rhs);
   
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_declare_fixed_var(int col_ind, char *name, double fixed_bound){

   if(name){
      printf("var %s [%i] is fixed to %f\n",
	     name, col_ind, fixed_bound);
   }else{
      printf("var [%i] is fixed to %f\n",
	     col_ind, fixed_bound);
   }
   return 0;

}
/*===========================================================================*/
/*===========================================================================*/
int prep_declare_coef_change(int row_ind, int col_ind, 
			     char *name, double a_val, 
			     double rhs){   
   if(name){
      printf("row [%i] with rhs %f: col %s [%i]: coeff improved to %f\n",
	     row_ind, rhs, name, col_ind, a_val);     
   }else{
      printf("row [%i] with rhs %f: col [%i]: coeff improved to %f\n",
	     row_ind, rhs, col_ind, a_val);     
   }
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_deleted_row_update_info(MIPdesc *mip, int row_ind)
{
   /* nothing fancy to do here */
   mip->mip_inf->rows[row_ind].is_redundant = TRUE;
   
   /* and update col sizes */
   COLinfo *cols = mip->mip_inf->cols; 

   int i, end, *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   //   double *r_matval = mip->row_matval; 

   end = r_matbeg[row_ind + 1];
   for(i = r_matbeg[row_ind]; i < end; i++){
      if(cols[r_matind[i]].var_type != 'F'){
	 (cols[r_matind[i]].col_size)--;
	 /* debug */
	 if(cols[r_matind[i]].col_size < 0){
	    printf("error in prep_deleted_row_update_info()\n");
	    return PREP_OTHER_ERROR;
	 }
      }
   }
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_initialize_mipinfo(PREPdesc *P)

{
   int i, j;
   double coef_val, fixed_obj_offset;  
   int row_ind, cont_var_cnt = 0, bin_var_cnt = 0, fixed_var_cnt = 0;   
   int  row_unbounded_cnt, max_row_size, max_col_size;
   int * row_coef_bin_cnt = NULL, *row_sign_pos_cnt = NULL;   
   char is_binary, is_bounded, unbounded_below, unbounded_above;
   int gen_type; /* 0 fractional, 1 integer, 2 binary */
   int col_size, col_coef_bin_cnt, col_coef_frac_cnt, col_sign_pos_cnt; 
   int *rows_integerized_var_ind = NULL, integerizable_var_num;
   char is_opt_val_integral = TRUE;
   int is_col_all_neg; /* if we convert all constraints to 'L' */
   int is_col_all_pos; /* if we convert all constraints to 'L' */
   int obj_size; /* number of nonzeros in objective function */

   MIPdesc *mip = P->mip;
   prep_stats *stats = &(P->stats);
   prep_params params = P->params;

   double etol = params.etol;
   double coeff_etol = 1e-15;
   int verbosity = params.verbosity;
   //   int p_level = prep_par.prep_level;
   //   int termcode; 

   /* fixme! objsense min max issue!!! will always assume that it is a 
      min problem here!!!! 
   */

   if(!mip){
      if(verbosity >= 3){
	 printf("prep_initialize_mipinfocollect_mipinfo(): Empty mip description...\n");
      }
      return(PREP_OTHER_ERROR);
   }
   
   int n = mip->n;
   int m = mip->m;
   int * matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;
   double * ub = mip->ub;
   double *lb = mip->lb;
   char *is_int = mip->is_int;   
   double *obj = mip->obj;
   char *sense = mip->sense;
   double *rhs = mip->rhs;
   //   char obj_sense = env->mip->obj_sense;

   MIPinfo *mip_inf = (MIPinfo *)calloc (1, sizeof(MIPinfo));
   COLinfo *cols = NULL;
   ROWinfo *rows = NULL;

   if(m > 0){
      rows = (ROWinfo *)calloc(m, sizeof(ROWinfo));   
      row_coef_bin_cnt = (int *)calloc(ISIZE,m);
      row_sign_pos_cnt = (int *)calloc(ISIZE,m);
      /* will be done for rows with only one cont var */
      /* this might be further improved in advanced option of 
	 preprocessor */
      integerizable_var_num = 0;
      rows_integerized_var_ind = (int *)malloc(ISIZE*m);
   }
   if(n > 0){
      cols = (COLinfo *)calloc(n, sizeof(COLinfo));           
   }

   max_col_size = 0;
   fixed_obj_offset = 0;
   obj_size = 0;
   
   for(i = 0; i < n; i++){
      is_binary = FALSE;
      is_bounded = FALSE;
      unbounded_below = FALSE;
      unbounded_above = FALSE;

      if (fabs(obj[i]) > etol) {
         obj_size++;
      }
      
      cols[i].var_type = 'I';
      if(lb[i] >= ub[i] + etol){
	 stats->col_infeas_ind = i;
	 return(PREP_INFEAS);
      }else if(lb[i] > (ub[i] - etol)){	 
	 cols[i].var_type = 'F';
	 fixed_obj_offset += obj[i]*ub[i];
	 fixed_var_cnt++;	 
	 if(ub[i] >= INF || ub[i] <= -INF){
	    stats->col_numeric_ind = i;
	    return(PREP_NUMERIC_ERROR);
	 }else{
	    is_bounded = TRUE;
	 }
      } else if (is_int[i]){
	 if(lb[i] > (-1.0 + etol) && 
	    ub[i] < (2.0 - etol)){
	    is_binary = is_bounded = TRUE; //is_lb_zero = TRUE;	    
	    cols[i].var_type = 'B';
	    bin_var_cnt++;		    
	 }else if(lb[i] > (-2.0 + etol) && 
		  ub[i] < (1.0 - etol)){
	    is_binary = is_bounded = TRUE; //is_lb_zero = TRUE;	    
	    cols[i].var_type = 'R';
	    bin_var_cnt++;		    
	 }
      }else {
	 cols[i].var_type = 'C';
	 cont_var_cnt++;
      }

      if(!is_bounded){
	 if(lb[i] <= -INF){
	    unbounded_below = TRUE;
	 }
	 if(ub[i] >= INF){
	    unbounded_above = TRUE;
	 }
	 if(!unbounded_below && !unbounded_above){
	    is_bounded = TRUE;
	 }
      }
      
      col_coef_bin_cnt = 0;
      col_coef_frac_cnt = 0;
      col_sign_pos_cnt = 0;
      is_col_all_neg = TRUE;
      is_col_all_pos = TRUE;

      for(j = matbeg[i]; j < matbeg[i+1]; j++){
	 row_ind = matind[j];
	 coef_val = matval[j];
	 rows[row_ind].size++;
	 /* for row types */	 
	 if(cols[i].var_type != 'F'){
	    if(is_int[i]){
	       if(is_binary){
		  rows[row_ind].bin_var_num++;
	       }
	    }else{
	       rows[row_ind].cont_var_num++;
	       if(rows[row_ind].cont_var_num < 2){
		  rows_integerized_var_ind[row_ind] = i;
	       }
	    }

	    //if(cols[i].var_type == 'U' ||
	    //  cols[i].var_type == 'L'){
	    //   rows[row_ind].fixable_var_num++;
	    // }
	 }else{
	    rows[row_ind].fixed_var_num++;
	 }
	 
	 /* for bound types */
	 if(!is_bounded){
	    if(unbounded_above){
	       if(coef_val > 0.0){
		  rows[row_ind].ub_inf_var_num++;
	       }else{
		  rows[row_ind].lb_inf_var_num++;
	       }
	    }
	    if(unbounded_below){
	       if(coef_val > 0.0){
		  rows[row_ind].lb_inf_var_num++;
	       }else{
		  rows[row_ind].ub_inf_var_num++;
	       }
	    }
	 }

	 /* for coef types */
	 if (!(coef_val-floor(coef_val) > coeff_etol &&
	     ceil(coef_val)-coef_val > coeff_etol)){
	    if((coef_val > 1.0 - coeff_etol && coef_val < 1.0 + coeff_etol) ||
	       (coef_val > -1.0 - coeff_etol &&
		coef_val < -1.0 + coeff_etol)) {	       
	       row_coef_bin_cnt[row_ind]++;
	       col_coef_bin_cnt++;
	    }
	 }else{
	    rows[row_ind].frac_coef_num++;
	    col_coef_frac_cnt++;
	 }

	 /* for sign types  and update bounds */
	 
	 if(coef_val > 0.0){
	    row_sign_pos_cnt[row_ind]++;
	    col_sign_pos_cnt++;
	    if(rows[row_ind].ub < INF){
	       if(ub[i] >= INF){
		  rows[row_ind].ub = INF;
	       }else{
		  rows[row_ind].ub += coef_val * ub[i];
	       }
	    }
	    if(rows[row_ind].lb > -INF){
	       if(lb[i] <= -INF){
		  rows[row_ind].lb = -INF;
	       }else{
		  rows[row_ind].lb += coef_val * lb[i];
	       }
	    }

	    if(is_col_all_neg){
	       if(sense[row_ind] != 'G'){
		  is_col_all_neg = FALSE;
	       }
	    }
	    if(is_col_all_pos){
	       if(sense[row_ind] != 'L'){
		  is_col_all_pos = FALSE;
	       }
	    }

	 }else if(coef_val < 0.0){
	    if(rows[row_ind].ub < INF){
	       if(lb[i] <= -INF){
		  rows[row_ind].ub = INF;
	       }else{
		  rows[row_ind].ub += coef_val * lb[i];
	       }
	    }
	    if(rows[row_ind].lb > -INF){
	       if(ub[i] >= INF){
		  rows[row_ind].lb = -INF;
	       }else{
		  rows[row_ind].lb += coef_val * ub[i];
	       }
	    }

	    if(is_col_all_neg){
	       if(sense[row_ind] != 'L'){ 
		  is_col_all_neg = FALSE;
	       }else{
		  

	       }
	    }
	    if(is_col_all_pos){
	       if(sense[row_ind] != 'G'){
		  is_col_all_pos = FALSE;
	       }
	    }
	 }

	 if(cols[i].var_type == 'F'){
	    rows[row_ind].fixed_obj_offset += obj[i]*ub[i];
	    rows[row_ind].fixed_lhs_offset += coef_val * ub[i];
	 }
      }


      /* check if this column is fixable */
      /* if so set type, update the obj_offset etc */

      col_size = cols[i].col_size = matbeg[i+1] - matbeg[i];

      if(is_col_all_neg || col_size <= 0){
	 //if((obj[i] > 0.0 && obj_sense == SYM_MAXIMIZE) ||
	 //(obj[i] < 0.0 && obj_sense == SYM_MINIMIZE)){
	 if(obj[i] < 0.0){
	    if(ub[i] >= INF){
	       stats->col_unbound_ind = i;
	       return(PREP_UNBOUNDED); /* fixme: unbounded return code */
	    }else{
	       /* total obj offset here is for fixed variables
		  later(if prep.is used) will include 'U' and 'L' 
		  variables */
	       //total_obj_offset += obj[i] * ub[i];
	       if(verbosity >= 2){
		  if(mip->colname){
		     printf("var %s [%i] is fixable to its upper bound: %f\n", 
			    mip->colname[i], i, ub[i]);
		  }else{
		     printf("var [%i] is fixable to its upper bound: %f\n", 
			    i, ub[i]);
		  }
		  cols[i].var_type = 'U';
	       }
	    }
	 }
      }
      if(is_col_all_pos || col_size <= 0){
	 if(obj[i] > 0.0){
	    if(lb[i] <= -INF){
	       stats->col_unbound_ind = i;
	       return(PREP_UNBOUNDED);
	    }else{
	       //total_obj_offset += obj[i] * lb[i];
	       if(verbosity >= 2){
		  if(mip->colname){
		     printf("var %s [%i] is fixable to its lower bound: %f\n", 
			    mip->colname[i], i, lb[i]);
		  }else{
		     printf("var [%i] is fixable to its lower bound: %f\n",
			    i, lb[i]);
		  }		     
		  cols[i].var_type = 'L';
	       }
	    }  
	 }
      }
      
      /* fill in col info - if not a fixed variable */
      if(col_size && cols[i].var_type != 'F'){
	 
	 if(col_size > max_col_size){
	    max_col_size = col_size;
	 }
	 
	 gen_type = 0;
	 if(col_coef_frac_cnt > 0){
	    gen_type = FRACTIONAL_VEC;
	 }else{
	    if(col_coef_bin_cnt < col_size){
	       gen_type = ALL_INTEGER_VEC;
	    }else {
	       gen_type = ALL_BINARY_VEC;
	    }
	 }	 
	 cols[i].coef_type = gen_type;

	 gen_type = 0;
	 if(col_sign_pos_cnt > 0){
	    if(col_sign_pos_cnt < col_size){
	       gen_type = MIXED_TYPE_VEC;
	    } else{
	       gen_type = ALL_POS_VEC;
	    }
	 } else{
	    gen_type = ALL_NEG_VEC;
	 }
	 cols[i].sign_type = gen_type;      
      }
   }
 
   max_row_size = 0;
   /* fill in row info */
   for(j = 0; j < m; j++){
      
      if(rows[j].size > max_row_size){
	 max_row_size = rows[j].size;
      }

      gen_type = 0;
      if(rows[j].cont_var_num > 0 ){
	 if(rows[j].bin_var_num > 0){
	    if( rows[j].cont_var_num + rows[j].bin_var_num + 
		rows[j].fixed_var_num 
		>= rows[j].size){
	       gen_type = BIN_CONT_TYPE;
	    } else{
	       gen_type = ALL_MIXED_TYPE;
	    }
	 } else {
	    if (rows[j].cont_var_num + rows[j].fixed_var_num < rows[j].size){
	       gen_type = INT_CONT_TYPE;
	    } else{
	       gen_type = CONTINUOUS_TYPE;
	    }
	 }
      }else{
	 if(rows[j].bin_var_num > 0){
	    if(rows[j].bin_var_num + rows[j].fixed_var_num < rows[j].size){
	       gen_type = BIN_INT_TYPE;
	    } else {
	       gen_type = BINARY_TYPE;
	    }
	 }else{
	    gen_type = INTEGER_TYPE;
	 }
      }
      
      rows[j].type = gen_type;
      gen_type = 0;
      
      row_unbounded_cnt = rows[j].lb_inf_var_num +
	 rows[j].ub_inf_var_num;

      if(row_unbounded_cnt == 0){
	 gen_type = ALL_BOUNDED_ROW;
      }else if(row_unbounded_cnt + rows[j].fixed_var_num < rows[j].size){
	 gen_type = MIXED_BOUNDED_ROW;
      }else{
	 gen_type = OPEN_ROW;
      }

      rows[j].bound_type = gen_type;
      gen_type = 0;

      if(rows[j].frac_coef_num > 0){
	 gen_type = FRACTIONAL_VEC;
      }else {
	 if(row_coef_bin_cnt[j] +  rows[j].fixed_var_num < rows[j].size){
	    gen_type = ALL_INTEGER_VEC;
	 }else {
	    gen_type = ALL_BINARY_VEC;
	 }
      }

      rows[j].coef_type = gen_type;
      gen_type = 0;
      
      if(row_sign_pos_cnt[j] > 0){
	 if(row_sign_pos_cnt[j] + rows[j].fixed_var_num < rows[j].size){
	    gen_type = MIXED_TYPE_VEC;
	 }else{
	    gen_type = ALL_POS_VEC;
	 }
      }else {
	 gen_type = ALL_NEG_VEC;
      }

      rows[j].sign_type = gen_type;      

      /* work on integerizable vars */
      
  
      if(sense[j] == 'E' && rows[j].cont_var_num == 1 && 
	 rows[j].coef_type != FRACTIONAL_VEC && 
	 prep_is_integral(rhs[j], coeff_etol) && 
	 prep_is_integral(rows[j].fixed_lhs_offset, coeff_etol)){
	 if(cols[rows_integerized_var_ind[j]].var_type != 'Z'){
	    cols[rows_integerized_var_ind[j]].var_type = 'Z';
	    integerizable_var_num++;
	 }
      }
      
      rows[j].sr_ub = rows[j].ub;
      rows[j].sr_lb = rows[j].lb;
   }

  /* work on obj */

   if(!(cont_var_cnt - integerizable_var_num)){
      for(i = 0; i < n; i++){
	 coef_val = obj[i];
	 if (!prep_is_integral(coef_val, coeff_etol)){
	    if(cols[i].var_type == 'F'){
	       if(ub[i] < etol && ub[i] > -etol){
		  continue;
	       }
	    }
	    is_opt_val_integral = FALSE;
	    break;
	 }
      }
   }else{
      is_opt_val_integral = FALSE;
   }

   //   printf("is_opt_val_integral = %i \n", is_opt_val_integral);


   /* finally prob type */
   gen_type = 0;
   if(cont_var_cnt > 0 ){
      if(bin_var_cnt > 0){
	 if(cont_var_cnt + bin_var_cnt + fixed_var_cnt < n){
	    gen_type = ALL_MIXED_TYPE;
	 }else{
	    gen_type = BIN_CONT_TYPE;
	 }
      }else{
	 if(cont_var_cnt + fixed_var_cnt < n){
	    gen_type = INT_CONT_TYPE;
	 } else{
	    gen_type = CONTINUOUS_TYPE;
	 }
      }
   }else{ 
      if(bin_var_cnt > 0){
	 if(bin_var_cnt + fixed_var_cnt < n){
	    gen_type  = BIN_INT_TYPE;
	 } else{
	    gen_type = BINARY_TYPE;
	 }
      }else{
	 gen_type = INTEGER_TYPE;
      }
   }

   mip_inf->prob_type = gen_type;
   mip_inf->cont_var_num = cont_var_cnt;
   mip_inf->binary_var_num = bin_var_cnt;
   mip_inf->fixed_var_num = fixed_var_cnt;
   mip_inf->max_row_size = max_row_size;
   mip_inf->max_col_size = max_col_size;   
   mip_inf->obj_size = obj_size;
   mip_inf->mat_density = mip->nz/(n*m);
   mip_inf->integerizable_var_num = integerizable_var_num;
   mip_inf->is_opt_val_integral = is_opt_val_integral;
   mip_inf->sum_obj_offset = fixed_obj_offset;

   if(mip->mip_inf){
      //      free_mipinfo(env->mip->mip_inf);
   }

   mip_inf->rows = rows;
   mip_inf->cols = cols;
   mip->mip_inf = mip_inf;  

   FREE(row_coef_bin_cnt);
   FREE(row_sign_pos_cnt);
   FREE(rows_integerized_var_ind);

   return(PREP_MODIFIED); 
}

/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
int prep_integerize_bounds(PREPdesc *P)
{
   /* Change the bounds of integer variables to floor/ceiling appropriately */
   int termcode = 0;
   MIPdesc *mip = P->mip;
   MIPinfo *mip_inf = mip->mip_inf;
   COLinfo *cols = mip_inf->cols;
   int i, b_cnt = 0, n = mip->n;
   double *ub = mip->ub;
   double *lb = mip->lb;
   double temp_fl, temp_cl;
   double diff_ub, diff_lb;
   double etol = P->params.etol;
   int verbosity = P->params.verbosity;

   //   int * bounds_updated = (int *)calloc(ISIZE,n);

   if(P->params.level > 2 && mip_inf->integerizable_var_num){
      for (i = 0; i < n; i++) {
	 if (cols[i].var_type == 'Z'){
	    termcode = prep_integerize_var(P, i);
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 }
      }
   }

   for (i = 0; i < n; i++) {
      if (cols[i].var_type != 'F' && 
	  cols[i].var_type != 'C') {
	 diff_ub = diff_lb = 0.0;
	 if (ub[i] < INF) {
	    temp_fl = floor(ub[i]);
	    temp_cl = ceil(ub[i]);
	    if(temp_cl - ub[i] < etol ){
	       ub[i] = temp_cl;
	    }else{
	       diff_ub = ub[i] - temp_fl; 
	       ub[i] = temp_fl;	       
	    }
	 }
	 if(lb[i] > -INF){
	    temp_fl = floor(lb[i]);
	    temp_cl = ceil(lb[i]);
	    if(lb[i] - temp_fl < etol){
	       lb[i] = temp_fl;
	    }else {
	       diff_lb = temp_cl - lb[i];
	       lb[i] = temp_cl;
	    }
	 }
	 if(diff_ub >= etol || diff_lb >= etol ){
	    if(ub[i] > lb[i] - etol && ub[i] < lb[i] + etol){
	       if(cols[i].var_type =='B'){
		  mip_inf->binary_var_num--;
	       }
	       mip_inf->fixed_var_num--;
	       cols[i].var_type = 'F';
	    }
	    b_cnt++;
	    if (verbosity>=3) {
	       if(mip->colname){
		  printf("integerized bounds [lb-ub] of variable %s:"
			 "%f - %f\n",
			 mip->colname[i],lb[i],ub[i]);
	       }else{
		  printf("integerized bounds [lb-ub] of variable: "
			 "%f - %f\n",
			 lb[i],ub[i]);
	       }
	    }
	 }
      }
   }

#if 0   
   if(keeptrack){
      P->mip_diff->bounds_integerized_num = b_cnt;
      P->mip_diff->bounds_integerized_ind = bounds_updated;
   }else{
      FREE(bounds_updated);
   }
#endif
   P->stats.bounds_integerized = b_cnt;
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

int prep_integerize_var(PREPdesc *P, int col_ind) {
   
   int j, k, row_ind, termcode = PREP_MODIFIED;
   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;
   double etol = P->params.etol;
   double coeff_etol = 1e-15;

   if(P->params.verbosity >= 2 ){
      printf("col %i is integerized\n", col_ind);
   }

   (P->stats.vars_integerized)++;
   mip->is_int[col_ind] = TRUE;
   cols[col_ind].var_type = 'I';
   if(mip->lb[col_ind] > (-1.0 + etol) && 
      mip->ub[col_ind] < (2.0 - etol)){
      cols[col_ind].var_type = 'B';
   }
   for(j = mip->matbeg[col_ind];
       j < mip->matbeg[col_ind+1]; j++){
      row_ind = mip->matind[j];
      if(cols[col_ind].var_type == 'B'){
	 rows[row_ind].bin_var_num++;
      }
      rows[row_ind].cont_var_num--;
      if(rows[row_ind].cont_var_num < 0){
	 printf("error: prep_integerize_var()\n");
	 return PREP_OTHER_ERROR;
      }else if(rows[row_ind].cont_var_num < 1){
	 if(rows[row_ind].bin_var_num){
	    if(rows[row_ind].bin_var_num + 
	       rows[row_ind].fixed_var_num 
	       >= rows[row_ind].size){
	       rows[row_ind].type = BINARY_TYPE;
	    }else{
	       rows[row_ind].type = BIN_INT_TYPE;
	    }
	 }else{
	    rows[row_ind].type = INTEGER_TYPE;
	 }
      }else if(rows[row_ind].cont_var_num == 1){
	 if(mip->sense[row_ind] == 'E' && 
	    rows[row_ind].coef_type != FRACTIONAL_VEC && 
	    prep_is_integral(mip->rhs[row_ind], coeff_etol) && 
	    prep_is_integral(rows[row_ind].fixed_lhs_offset, coeff_etol)){
	    for(k = mip->row_matbeg[row_ind];
		k < mip->row_matbeg[row_ind + 1]; k++){
	       if(cols[mip->row_matind[k]].var_type == 'C'){
		  termcode = prep_integerize_var(P, mip->row_matind[k]);
		  break;
	       }
	    }
	 }
      }
      if(prep_quit(termcode)){
	 break;
      }
   }
   return termcode;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_fill_row_ordered(PREPdesc *P)
{ 
   /*
     recreates 'A' matrix using three matrices just like the standard
     notation. However, matrices contain row representations rather than
     column
   */ 

   int i, j, *o_ind;
   int row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg, *r_lengths; 
   double * matval, *r_matval, *rhs;
   MIPdesc *mip = P->mip;
   int n = mip->n;
   int m = mip->m;
   int nz = mip->nz;
   char *sense, *o_sense;
   int *u_col_ind, *u_row_ind;
   //   ROWinfo * rows = mip->mip_inf->rows;

   //  if(!rows){
      /* debug */
   //  printf("error in prep_fill_row_ordered - 1");
   //  return PREP_OTHER_ERROR;
   // }

   matval = mip->matval;
   matind = mip->matind;
   matbeg = mip->matbeg;
   rhs = mip->rhs;
   sense = mip->sense;

   /* allocate space for different arrays */

   r_matval = (mip->row_matval = (double *)malloc(nz*DSIZE)); 
   r_matind = (mip->row_matind = (int *)malloc(nz*ISIZE)); 
   r_matbeg = (mip->row_matbeg = (int *)malloc((m+1)*ISIZE));
   r_lengths = (mip->row_lengths = (int *)calloc(m,ISIZE));
   o_sense = (mip->orig_sense = (char *)malloc(m *CSIZE));
   o_ind = (mip->orig_ind = (int *)malloc(n*ISIZE));
   u_col_ind = (P->user_col_ind) = (int *)malloc(n*ISIZE);
   u_row_ind = (P->user_row_ind) = (int *)malloc(m*ISIZE);
   /* these are initialized here, we have to visit this function anyway */
   
   srand ( time(NULL) ); 

   /* first get row legths */   
   for(i = 0; i < n; i++){
      /* get orig indices here */
      o_ind[i] = u_col_ind[i] = i;
      for(j = matbeg[i]; j < matbeg[i+1]; j++){
	 r_lengths[matind[j]]++;
      }
   }

   r_matbeg[0] = 0;

   /* fill in matbegs */
   for(i = 0; i < m; i++){
      u_row_ind[i] = i;
      r_matbeg[i + 1] = r_matbeg[i] + r_lengths[i];
   }

   /* get matrix, change 'G' rows to 'L'*/
   for(i = 0; i < n; i++){
      for(j = matbeg[i]; j < matbeg[i+1]; j++){
	 row_ind = matind[j];
	 elem_ind = r_matbeg[row_ind];
	 r_matind[elem_ind] = i;
	 if(sense[row_ind] == 'G'){
	    matval[j] = -matval[j];
	 }
	 r_matval[elem_ind] = matval[j];
	 r_matbeg[row_ind] = elem_ind + 1;
      }
   }
   /* and update matbegs, rhs, and rows with 'G' to 'L'*/
   memcpy(o_sense, sense, CSIZE*m);   
   
   for(i = 0; i < m; i++){
      r_matbeg[i] -= r_lengths[i];
      if(sense[i] == 'G'){
	 sense[i] = 'L';
	 rhs[i] = -rhs[i];
      }
   }   

   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_cleanup_desc(PREPdesc *P)
{ 
   
   int i, j, col_nz, col_num, row_num, fixed_nz, *fixed_ind, *o_ind;
   int row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg, *r_lengths; 
   double *ub, *lb, *matval, *r_matval, *obj, *rhs, *rngval, *fixed_val;
   double obj_offset, debug_offset;
   int new_del_cnt;

   MIPdesc *mip = P->mip;
   int n = mip->n;
   int r_ind, m = mip->m;
   int nz = mip->nz;
   char *is_int, *sense, *o_sense, **colnames;
   ROWinfo * rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;
   //  if(!rows){
      /* debug */
   //  printf("error in prep_fill_row_ordered - 1");
   //  return PREP_OTHER_ERROR;
   // }

   prep_stats * stats = &(P->stats);
   prep_params params = P->params;

   double etol = params.etol;
   double coeff_etol = 1e-15;
   int deleted_row_cnt = stats->rows_deleted;
   //int vars_fixed = stats->vars_fixed;
   //int keep_row_ordered = params.keep_row_ordered;
   //int reduce_mip = params.reduce_mip;

   int old_start, *row_new_inds = NULL;
   int fixed_num = stats->vars_fixed + mip->mip_inf->fixed_var_num;

   obj_offset = 0.0;
   debug_offset = 0.0;
   sense = mip->sense;
   rngval = mip->rngval;
   o_sense = mip->orig_sense;
   rhs = mip->rhs;
   obj = mip->obj;
   ub = mip->ub;
   lb = mip->lb;
   is_int = mip->is_int;
   
   fixed_nz = 0;
   fixed_ind = mip->fixed_ind = (int *)malloc(n*ISIZE);
   fixed_val = mip->fixed_val = (double *)malloc(n*DSIZE);
   o_ind = mip->orig_ind;
   
   if(!params.reduce_mip || fixed_num == n || stats->rows_deleted == m){
      if(fixed_num == n || stats->rows_deleted == m){
	 /* get fixed nz vals */	 
	 for(i = 0; i < n; i++){
	    if(cols[i].var_type == 'F'){	 
	       if(!prep_is_equal(mip->ub[i], 0.0, etol)){
		  fixed_ind[fixed_nz] = i;
		  fixed_val[fixed_nz] = mip->ub[i];
	       }
	    }else{
	       if(obj[i] > 0.0){
		  if(lb[i] <= -INF){
		     return PREP_UNBOUNDED;
		  }else{
		     if(!prep_is_equal(mip->lb[i], 0.0, etol)){
			fixed_ind[fixed_nz] = i;
			fixed_val[fixed_nz] = mip->lb[i];
			obj_offset += obj[i]*mip->lb[i];
		     }
		  }
	       }else if(obj[i] < 0.0){
		  if(ub[i] >= INF){
		     return PREP_UNBOUNDED;
		  }else{
		     if(!prep_is_equal(mip->ub[i], 0.0, etol)){
			fixed_ind[fixed_nz] = i;
			fixed_val[fixed_nz] = mip->ub[i];
			obj_offset += obj[i]*mip->ub[i];			
		     }
		  }
	       }
	    }
	    fixed_nz++;
	 }
	 mip->fixed_n = fixed_nz;
	 return PREP_SOLVED;
      }
      return PREP_UNMODIFIED;
   }

   if(!fixed_num && !stats->rows_deleted){
      return PREP_UNMODIFIED;
   }
   
   row_new_inds = (int *)calloc(m,ISIZE);
   
   mip->alloc_n = n;
   mip->alloc_m = m;
   mip->alloc_nz = nz;
   
   matval = mip->matval;
   matind = mip->matind;
   matbeg = mip->matbeg;

   colnames = mip->colname;
   
   /* first get new row indices */
   
   row_num = 0;
   
   for(i = 0; i < m; i++){
      if(!(rows[i].is_redundant)){
	 row_new_inds[i] = row_num;
	 row_num++;
      }
      rows[i].size = 0; /* might have zero matval vals that we want 
			   to exclude, so we will re-evaluate */
   }
   
   /* debug */
   /* ------ */
   if(row_num != m - deleted_row_cnt){
      printf("error: missing rows \n");
      return PREP_OTHER_ERROR;
   }
   /* ------ */
   
   /* first fill in col ordered */
   
   col_nz = col_num = 0;
   old_start = 0;
   for(i = 0; i < n; i++){
      if(cols[i].var_type != 'F'){
	 for(j = old_start; j < matbeg[i+1]; j++){
	    r_ind = matind[j];
	    if(!(rows[r_ind].is_redundant)){
	       if(!prep_is_equal(matval[j], 0.0, coeff_etol)){
		  matind[col_nz] = row_new_inds[r_ind];
		  matval[col_nz] = matval[j];
		  rows[r_ind].size++;
		  col_nz++;
	       }
	    }
	 }
	 
	 /* check if this column has any nonzeros, otherwise 
	    fix it*/
	 
	 if(col_nz == matbeg[col_num]){
	    cols[i].var_type = 'F';
	    if(obj[i] >= 0.0){
	       obj_offset += obj[i]*lb[i];
	       fixed_val[fixed_nz] = lb[i];
	    }else{
	       obj_offset += obj[i]*ub[i];
	       fixed_val[fixed_nz] = ub[i];
	    }
	    if(!prep_is_equal(fixed_val[fixed_nz], 0.0, etol)){
	       fixed_ind[fixed_nz++] = i;
	    }
	    stats->vars_fixed++;
	 }else{
	    o_ind[col_num] = i;
	    obj[col_num] = obj[i];
	    ub[col_num] = ub[i];
	    lb[col_num] = lb[i];
	    is_int[col_num] = is_int[i];
	    if(col_num != i){
	       cols[col_num] = cols[i];
	       if(colnames){
		  strcpy(colnames[col_num], colnames[i]);
	       }
	    }
	    cols[col_num].col_size = col_nz - matbeg[col_num];
	    old_start = matbeg[i+1];
	    matbeg[(++col_num)] = col_nz;
	    /* debug */
	    /* ----------- */
	    if(cols[col_num - 1].col_size <= 0){
	       printf("error: empty size column \n");
	       return PREP_OTHER_ERROR;
	    }
	    /* ---------- */
	 }
      }else{
	 
	 /* debug */
	 /*---- -*/
	 if(stats->vars_aggregated <= 0){
	    if(!prep_is_equal(ub[i], lb[i], etol)){
	       printf("error: not fixed column? \n");
	       return PREP_OTHER_ERROR;
	    }
	 }
	 /* ----- */
	 debug_offset += obj[i]*ub[i];
	 old_start = matbeg[i+1];
	 if(!prep_is_equal(ub[i], 0.0, etol)){
	    fixed_ind[fixed_nz] = i;
	    fixed_val[fixed_nz++] = ub[i];
	 }	    
      }
   }
      
   /* debug */
   /* -------------- */
   if(col_num != n - (stats->vars_fixed + mip->mip_inf->fixed_var_num)){
      printf("error: missing cols \n");
      return PREP_OTHER_ERROR;
   }
   /* ---------- */
   
   /* now update row_info */
   row_num = 0; 
   new_del_cnt = 0;
   for(i = 0; i < m; i++){     
      if(!(rows[i].is_redundant) && rows[i].size >= 0){
	 row_new_inds[row_num + new_del_cnt] = row_num;
	 if(rows[i].size == 0){	 
	    rows[i].is_redundant = TRUE;
	    new_del_cnt++;
	 }else{      
	    if(i != row_num){
	       rows[row_num] = rows[i];
	       sense[row_num] = sense[i];
	       if(sense[row_num] == 'R'){
		  rngval[row_num] = rngval[i];
	       }
	    }
	    rhs[row_num] = rhs[i] - rows[i].fixed_lhs_offset;
	    row_num++;
	 }
      }
   }
   
   stats->rows_deleted += new_del_cnt;

   /* now convert it row ordered if asked*/
   /* get row_lengths and fill in r_matbeg, sense, rhs etc */
   
   r_matbeg = mip->row_matbeg;
   r_matind = mip->row_matind;
   r_matval = mip->row_matval;
   r_lengths = mip->row_lengths;      
   
   for(i = 0; i < row_num; i++){
      r_lengths[i] = rows[i].size;
      r_matbeg[i+1] = r_matbeg[i] + 
	 rows[i].size;
   }
   
   /* debug */
   if(r_matbeg[row_num ] != col_nz){
      printf("error; missing nonzeros\n");
      return PREP_OTHER_ERROR;
   }
   
   /* */

   for(i = 0; i < col_num; i++){
      for(j = matbeg[i]; j < matbeg[i+1]; j++){
	 matind[j] = row_new_inds[matind[j]];
	 row_ind = matind[j];
	 elem_ind = r_matbeg[row_ind];
	 r_matind[elem_ind] = i;
	 r_matval[elem_ind] = matval[j];
	 r_matbeg[row_ind] = elem_ind + 1;
      }
   }

   for(i = 0; i < row_num; i++){
      r_matbeg[i] -= r_lengths[i];
   }

   /*}else{
     FREE(mip->row_matbeg);
     FREE(mip->row_matind);
     FREE(mip->row_matval);
     FREE(mip->row_lengths);
     FREE(mip->orig_sense);
     }*/
      
   mip->n = col_num;
   mip->m = row_num;
   mip->nz = col_nz;
   mip->obj_offset = mip->mip_inf->sum_obj_offset + obj_offset;     
   mip->fixed_n = fixed_nz;
      
   FREE(row_new_inds);

   if(mip->n <= 0 || mip->m <= 0){
      return PREP_SOLVED;
   }
   
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
char prep_is_equal(double lval, double rval, double etol)
{

   double diff = lval - rval; 
   
   if(diff < etol &&
      diff > -etol){
      return TRUE;
   }else{
      return FALSE;
   }
}
/*===========================================================================*/
/*===========================================================================*/

char prep_is_integral(double val, double etol)
{

   if(val - floor(val) < etol ||
      ceil(val) - val < etol){
      return TRUE;
   }else{
      return FALSE;
   }
   
}

/*===========================================================================*/
/*===========================================================================*/

char prep_quit(int termcode){
   
   if(termcode != PREP_UNMODIFIED &&
      termcode != PREP_MODIFIED) {
      return TRUE;
   }else {
      return FALSE;
   }
}

/*===========================================================================*/
/*===========================================================================*/

int prep_report(PREPdesc *P, int termcode)
{

   MIPdesc *mip = P->mip;
   int i;
   prep_stats stats = P->stats;
   
   switch(termcode){
    case PREP_INFEAS:
      printf("Preprocessing detected infeasibility...");
      if(stats.col_infeas_ind >= 0 ||
	 stats.row_infeas_ind >= 0){
	 printf("while improving bounds of \n\t");	  
	 if(stats.col_infeas_ind >= 0){
	    printf("variable ");
	    if(mip->colname){
	       printf("%s ", mip->colname[stats.col_infeas_ind]); 
	    }
	    printf("[%i]", stats.col_infeas_ind);
	    if(stats.row_infeas_ind >= 0){
	       printf(" on the ");
	    }
	 }
	 if(stats.row_infeas_ind >= 0){
	    printf("row [%i]", stats.row_infeas_ind);
	 }
	 printf("\n");	     
      }
      break;
    case PREP_UNBOUNDED:
      printf("Preprocessing detected unbounded problem...");	  
      if(stats.col_unbound_ind >= 0){
	 printf("while improving bounds on \n");	  
	 if(mip->colname){
	    printf("variable %s [%i]\n", 
		   mip->colname[stats.col_unbound_ind], 
		   stats.col_unbound_ind);
	 }else{
	    printf("variable [%i]\n", 
		   stats.col_unbound_ind);		
	 }
      }
      break;
    case PREP_NUMERIC_ERROR:
      printf("Preprocessing detected numerical problems ");	  
      if(stats.col_numeric_ind >= 0){
	 printf("while improving bounds on \n");	  
	 if(mip->colname){
	    printf("variable %s [%i]\n", 
		   mip->colname[stats.col_numeric_ind], 
		   stats.col_numeric_ind);
	 }else{
	    printf("variable [%i]\n", 
		   stats.col_numeric_ind);		
	 }
      }
      break;
    case PREP_OTHER_ERROR:
      printf("Preprocessing - unknown error.. ignoring presolve...\n");
      break;
    case PREP_SOLVED:
      printf("Preprocessing found the optimum:\n");	  	  
      printf("Solution Cost: %f\n:", mip->obj_offset);
      if (mip->colname){ 
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf("Column names and values of nonzeros in the solution\n");
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < mip->fixed_n; i++){
	    printf("%8s %10.3f\n", P->orig_mip->colname[mip->fixed_ind[i]],
		   mip->fixed_val[i]);
	 }
	 printf("\n");
      }else{
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf("User indices and values of nonzeros in the solution\n");
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < mip->fixed_n; i++){
	    printf("%7d %10.3f\n", mip->fixed_ind[i], mip->fixed_val[i]);
	 }
	 printf("\n");
      }
      break;	  
    default:
      printf("Preprocessing finished...\n ");	  
      
      if(stats.coeffs_changed + 
	 stats.bounds_tightened + 
	 stats.rows_deleted + 
	 stats.vars_fixed + 
	 stats.vars_aggregated +
	 stats.vars_integerized > 0){
	 if(stats.coeffs_changed > 0){
	    printf("\t coefficients modified: %i\n",
		   stats.coeffs_changed);		       
	 }
	 if(stats.bounds_tightened > 0){
	    printf("\t bounds improved: %i\n", 
		   stats.bounds_tightened);
	 }	     
	 if(stats.rows_deleted + 
	    stats.vars_fixed > 0){
	    if(stats.rows_deleted > 0){
	       printf("\t constraints removed: %i\n", 
		      stats.rows_deleted);
	       //printf("\t %i remained\n", mip->m);
	    }
	    if(stats.vars_fixed > 0){
	       printf("\t variables fixed: %i\n", stats.vars_fixed);
	       //printf("\t %i remained\n", mip->n);
	    }	     
	 }
	 if(stats.vars_aggregated > 0){
	    printf("\t variables aggregated: %i\n", stats.vars_aggregated);
	 }
	 if(stats.vars_integerized > 0){
	    printf("\t variables integerized: %i\n", stats.vars_integerized);
	 }
	 
      }else{
	 printf("\t with no modifications...\n");
      }
      printf("Problem has \n"
	     "\t %i constraints \n"
	     "\t %i variables \n"
	     "\t %i nonzero coefficients\n", 
	     mip->m, mip->n, mip->nz);	     	  
      break;
   }
   printf("\n");
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_restore_rootdesc(sym_environment *env)

{
      //int bvarnum = env->base->varnum, bind = 0;
   int i, user_size = env->rootdesc->uind.size;// uind = 0;
   //int *bvar_ind = env->base->userind; 
   int *user_ind = env->rootdesc->uind.list;
   
   env->base->cutnum = env->mip->m;

   if(user_size == env->mip->n){
      return PREP_UNMODIFIED;
   }else{
      for(i = 0; i < env->mip->n; i++){
	 user_ind[i] = i;
      }
   }
   
   env->rootdesc->uind.size = env->mip->n;
   
   return PREP_MODIFIED;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_close_desc(PREPdesc *P){

   
   if(P){
      if(P->sr){
	 prep_free_sr_desc(P->sr);
      }
      if(P->d_sr){
	 prep_free_sr_desc(P->d_sr);
      }

      if(P->mip){
	 free_mip_desc(P->mip);
      }
      /* fixme - add impl stuff here - disabled now*/
      
      FREE(P->user_col_ind);
      FREE(P->user_row_ind);
      FREE(P->stats.nz_coeff_changed);
      FREE(P);
   }
   
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_free_sr_desc(SRdesc *sr){

   if(sr){
      FREE(sr->obj_max);
      FREE(sr->matval_max);
      FREE(sr->matind_max);
      FREE(sr->ratio_max);

      FREE(sr->obj_min);
      FREE(sr->matval_min);
      FREE(sr->matind_min);
      FREE(sr->ratio_min);

      FREE(sr->fixed_ind);
      FREE(sr->tmp_ind);
      
      FREE(sr);
   }

   return 0;
}
/*===========================================================================*/
/*===========================================================================*/

#if 0
char prep_get_rnd_row_lb(ROWinfo row, double var_offset, 
			 int var_status, int is_int, double a_val, 
			 double *lb)
{
   double row_lb = row.lb;
   double total_offset;

   if(row.frac_coef_num > 1){
      return row_lb;
   }else if(row.frac_coef_num == 1 && var_status == VAR_IN_BOUND && 
	    prep_is_integral(a_val)){
      return row_lb;
   }
   
   if(row.cont_var_num > 1 ){
      return return_lb;
   }else if{row.cont_var_num == 1 && var_status == VAR_IN_BOUND &&
	       is_int){
      return return_lb;
   }
   
   if(var_status != VAR_IN_BOUND){
      var_offset = 0.0;
   }
   
   total_offset = row->fixed_lhs_offset + var_offset;
   row_lb -= total_offset;
   row_lb = prep_rnd_integral(row_lb, etol, RND_CEIL);
   
   return (row_lb + total_offset); 
}

/*===========================================================================*/
/*===========================================================================*/

char prep_get_rnd_row_ub(ROWinfo row, double var_offset, 
			 int var_status, int is_int, double a_val, 
			 double *ub)
{
   double row_ub = row.ub;
   double total_offset;

   if(row.frac_coef_num > 1){
      return 0;
   }else if(row.frac_coef_num == 1 && var_status == VAR_IN_BOUND && 
	    prep_is_integral(a_val)){
      return 0;
   }
   
   if(row.cont_var_num > 1 ){
      return 0;
   }else if{row.cont_var_num == 1 && var_status == VAR_IN_BOUND &&
	       is_int){
      return 0;
   }
   
   if(var_status != VAR_IN_BOUND){
      var_offset = 0.0;
   }
   
   total_offset = row->fixed_lhs_offset + var_offset;
   row_lb -= total_offset;
   row_lb = prep_rnd_integral(row_lb, etol, RND_FLOOR);
   
   *lb = row_lb + total_offset;
   
   return 0;  
}
#endif
/*===========================================================================*/
/*===========================================================================*/

#if 0
int find_sign(double val, double etol){

   if(val > etol){
      return POS_VAL;
   }else if(val < -etol){
      return NEG_VAL;
   }else{
      return ZERO_VAL;
   }
   
}
#endif
