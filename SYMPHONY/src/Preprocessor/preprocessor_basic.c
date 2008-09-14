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
   int can_iterate = TRUE, iter_cnt = 0, iter_cnt_limit, p_level; 
   
   int verbosity;
   double a_val, etol;// min_ub, max_lb; 
   char do_sr_rlx, do_aggr_row_rlx;// fix_var; 
   int j, m, n, nz, *r_matbeg, *r_matind; 
   //int i, max_size, min_size, *max_ind, *min_ind; 
   double *obj, *rhs, *r_matval, *ub, *lb;// old_val, old_lb, old_ub;  
   char *sense; 
   //int * updated_cols_ind, updated_cols_cnt;
   //int * updated_rows_ind, updated_rows_cnt;
   int col_ind, row_ind;
   //int * modified_cols_ind, * modified_rows_ind;
   //   int modified_cols_cnt, modified_rows_cnt;
   //char *is_col_updated, *is_row_updated;
   char var_type;// sr_termcode; 


   /* first initialize P and mip etc */

   prep_stats *stats = &(P->stats);
   prep_params params = P->params;
   
   verbosity = params.verbosity; /* debug */
   p_level = params.level;
   etol = params.etol;
   iter_cnt_limit = params.iteration_limit;
   do_sr_rlx = params.do_single_row_rlx;      
   do_aggr_row_rlx = params.do_aggregate_row_rlx;
   
   MIPdesc *mip = P->mip;
   MIPinfo *mip_inf = mip->mip_inf;
   COLinfo * cols = mip_inf->cols;
   ROWinfo row, *rows = mip_inf->rows;
   
   /* first integerize the bounds */
   /* can be embedded somewhere in basic prep down*/
   /* for now let it be */
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
   ub = mip->ub;
   lb = mip->lb;

   obj = mip->obj;
   sense = mip->sense;
   rhs = mip->rhs;
   
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

   /* main preprocessing loop */

   while(iter_cnt < iter_cnt_limit){// && updated_rows_cnt > 0){
      
      iter_cnt++;
      //modified_rows_ind = 0;
      //modified_cols_ind = 0;

      PRINT(verbosity, 2, ("Basic iteration number: %d\n", iter_cnt));

      /* check the updated bounds and cols to iterate on*/
      /* for each updated row and column do bla bla...*/
      /* while iterating update the new cols and rows indices*/
      
      /* do this only once at the end of first iteration*/
      //if(prep_level > 2 && (do_sr_rlx || do_aggr_row_rlx) && 
      // iter_cnt == 1){
      /*=====================================================================*/
      /*=====================================================================*/
      
      for(row_ind = 0; row_ind < m; row_ind++){

	 row = rows[row_ind];

	 if(row.is_redundant){
	    continue;
	 }

	 can_iterate = FALSE;
	 if(row.ub >= INF && row.lb <= -INF){
	    //    if(prep_solve_sr_rlx(P, 1, &row_ind) != SR_BOUNDS_UPDATED){
	    continue;
	 }else{
	    can_iterate = TRUE;
	 }
	 
	 if(can_iterate){
	    
	    termcode = prep_check_redundancy(P, row_ind, FALSE);
					     
	    if(prep_quit(termcode)){
	       break;
	    }
	    if(row.is_redundant){
	       continue;
	    }
	 }
	 if(can_iterate){

	    for(j = r_matbeg[row_ind]; j < r_matbeg[row_ind+1]; j++){
	       
	       /* bounds improvement, variable fixing, coeff improvement 
		  all here */
	       /* while checking bounds improvement, update bounds */	       
	       /* when an unbounded column is bounded, keep the info that 
		  which rows' bounds should be updated? */
	       col_ind = r_matind[j];
	       a_val = r_matval[j];
	       var_type = cols[col_ind].var_type;
	       if(var_type != 'F'){
		  if(row.ub_inf_var_num <= 1 ||
		     row.lb_inf_var_num <= 1){
		     if(ub[col_ind] >= INF || lb[col_ind] <= -INF){
			if((a_val > 0.0 && ub[col_ind] >= INF) ||
			   (a_val < 0.0 && lb[col_ind] <= -INF) ||
			   sense[row_ind] == 'E'){
			   termcode = 
			      prep_force_row_bounds(P, row_ind, col_ind, j);
			   if(prep_quit(termcode)){
			      return termcode;
			   }
			}
		     }
		  }
		  termcode = prep_improve_variable(P, col_ind, row_ind, j);
		  if(prep_quit(termcode)){
		     return termcode;
		  }
	       }    
	    }
	    
#if 0	 
	    if(do_sr_rlx){
	       if(iter_cnt <= 2){
		  sr_termcode = prep_solve_sr_rlx(P, 1, &row_ind);
	       }
	    }
#endif
	    //if(iter_cnt == iter_cnt_limit - 1){
	    /* do gcd thing and see for infeasiblity and further bound 
	       tightening:  ax <= b
	       gcd(a) = p, make rhs = kp <= b
	    */
	    //}
	 }
      

	 /* one last time to detect more redundancy */
#if 0
	 if(do_sr_rlx){
	    for(row_ind = 0; row_ind < m; row_ind++){
	       
	       if(row.is_redundant){
		  continue;
	       }
	       sr_termcode = prep_solve_sr_rlx(P, 1, &row_ind);
	       if(sr_termcode == SR_BOUNDS_UPDATED){
		  prep_check_redundancy(rows, row_ind, sense[row_ind], 
					rhs[row_ind], 
					cols, r_matind, r_matbeg, r_matval, 
					ub, lb, 
					etol, TRUE);
		  if(row.is_redundant){
		     redundant_row_cnt++;
		     if(verbosity >= 1){
			printf("row %i is redundant: ub: %f \t lb: %f \t" 
			       "sense: %c \t rhs: %f \n", 
			       row_ind, row.ub >= INF ? -1 : row.ub, 
			       row.lb <= -INF ? 1 : row.lb, sense[row_ind], 
			       rhs[row_ind]);
		     }
		  }
	       }
	    }
	 }
#endif 
	 
      }

      /* now check if we have any col with col.size = 0 
	 because we cant catch them above */
      for(col_ind = 0; col_ind < n; col_ind++){
	 if(cols[col_ind].var_type != 'F' && cols[col_ind].col_size == 0){
	    termcode = prep_improve_variable(P, col_ind, row_ind, j);
	    if(prep_quit(termcode)){
	       return termcode;
	    }
	 }
      }	       

   }
      
   /* check and imply if we need to reduce mip by 
      removing fixed rows and cols */      
   termcode = prep_update_mip(mip, *stats, params, etol);
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

   if(rows[row_ind].lb <= -INF && rows[row_ind].ub >= INF){
      return PREP_UNMODIFIED;
   }

   double a_val = mip->row_matval[a_loc];

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
      if(a_val > 0.0 && ub[col_ind] >= INF){
	 if(rows[row_ind].lb > -INF){
	    //ub[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].ub >= INF){
	    // printf("error -2 in prep_force_row_bounds()\n");
	    //  return PREP_OTHER_ERROR;
	    //}	 
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].lb)/a_val);
	    
	    fix_type = IMPROVE_UB;
	    row_ub_improved = TRUE;
	 }
      }else if(a_val < 0.0 && lb[col_ind] <= -INF){
	 if(rows[row_ind].lb > -INF){
	    // lb[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].ub >= INF){
	    //   printf("error -2 in prep_force_row_bounds()\n");
	    //  return PREP_OTHER_ERROR;
	    //}
	    
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].lb)/a_val);
	    
	    fix_type = IMPROVE_LB;
	    row_ub_improved = TRUE;
	 }
      }      
   }else if(sense == 'E'){
      if(a_val > 0.0 && lb[col_ind] <= -INF){
	 if(rows[row_ind].ub < INF){
	    
	    //lb[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].lb <= -INF){
	    //   printf("error -2 in prep_force_row_bounds()\n");
	    //   return PREP_OTHER_ERROR;
	    // }
	    
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].ub)/a_val);
	    
	    fix_type = IMPROVE_LB;
	    row_lb_improved = TRUE;	    
	 }      
      }else if(a_val < 0.0 && ub[col_ind] >= INF){
	 if(rows[row_ind].ub < INF){
	    //ub[col_ind] = 0.0;
	    //prep_get_row_bounds(mip, row_ind);
	    
	    /* debug check */
	    //if(rows[row_ind].lb <= -INF){
	    //  printf("error -2 in prep_force_row_bounds()\n");
	    //  return PREP_OTHER_ERROR;
	    //}
	    
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].ub)/a_val);

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
      termcode = prep_modified_col_update_info(P, col_ind, new_bound, 
					       fix_type);
				    
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
int prep_improve_variable(PREPdesc *P, int col_ind, int row_ind, int a_loc) 
			  
{

   int i, fix_type = FIX_NO_BOUND, termcode = PREP_UNMODIFIED;

   MIPdesc *mip = P->mip;

   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   //ROWinfo row = rows[row_ind];

   //   int *r_matbeg = mip->row_matbeg;
   // int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval; 
   
   double *ub = mip->ub;
   double *lb = mip->lb;
   double *obj = mip->obj;

   char is_int = mip->is_int[col_ind];
   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double a_val = r_matval[a_loc];
   
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
	 if(obj[col_ind] >= 0.0 && ( (sense == 'G' && a_val < 0.0) || 
				     (sense == 'L' && a_val > 0.0) ||
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
	 }else if(obj[col_ind] <= 0.0 && ( (sense == 'G' && a_val > 0.0) || 
					   (sense == 'L' && a_val < 0.0) || 
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
	 if(a_val > 0){	
	    switch(sense){
	     case 'L':
		/* try to fix */
		if(rows[row_ind].lb > -INF){
		   new_lb = rows[row_ind].lb + a_val;
		   if(new_lb > rhs)
		      fix_to_lb = TRUE;
		}
		/* try to improve */
		if(!fix_to_lb){
		   if(rows[row_ind].ub < INF){
		      new_ub = rows[row_ind].ub - a_val;
		      if(new_ub < rhs){
			 
			 /* update coefs */			 
			 r_matval[a_loc] = rows[row_ind].ub - rhs;
			 mip->rhs[row_ind] = new_ub;

			 /* debug */
			 if(r_matval[a_loc] < -etol){
			    printf("error in prep_improve_variable()\n");
			    return PREP_OTHER_ERROR;
			 } 
			 
			 /* update bounds */
			 rows[row_ind].ub += 
			    (r_matval[a_loc] - a_val) * 
			    ub[col_ind];
			 
			 improve_coef = TRUE; 
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
		   new_lb = rows[row_ind].lb + a_val;		   
		   if(new_lb > rhs){
		      fix_to_lb = TRUE;
		   }
		}
		if(rows[row_ind].ub < INF){
		   new_ub = rows[row_ind].ub - a_val;
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
	 }else{
	    switch(sense){
	     case 'L':
		if(rows[row_ind].lb > -INF){
		   new_lb = rows[row_ind].lb - a_val; 	 	     
		   if(new_lb > rhs)
		      fix_to_ub = TRUE;
		}
		if(!fix_to_ub){
		   if(rows[row_ind].ub < INF){
		      new_ub = rows[row_ind].ub + a_val;
		      if(new_ub < rhs){
			 //improve_offset = rhs - new_ub; 

			 /* update coef*/			 
			 r_matval[a_loc] -= new_ub - rhs;

			 /* debug */
			 if(r_matval[a_loc] > etol){
			    printf("error -3 in prep_improve_variable()\n");
			    return PREP_OTHER_ERROR;
			 } 
			 
			 /* update bounds */
			 rows[row_ind].lb += 
			    (r_matval[a_loc] - a_val) * 
			    lb[col_ind];
			 
			 improve_coef = TRUE;
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
		   new_lb = rows[row_ind].lb - a_val; 	 	     
		   if(new_lb > rhs){
		      fix_to_ub = TRUE;
		   }
		}
		if(rows[row_ind].ub < INF){
		   new_ub = rows[row_ind].ub + a_val;	 
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
	    
	    for(i = mip->matbeg[col_ind]; i < 
		   mip->matbeg[col_ind + 1]; i++){
	       if(mip->matind[i] == row_ind){
		  mip->matval[i] = r_matval[a_loc];
		  break;			       
	       }
	    }

	    /* debug */
	    if(i == mip->matbeg[col_ind + 1]){
	       printf("error -1 in prep_improve_variable()\n");
	       return PREP_OTHER_ERROR;			    
	    }
	    
	    if(verbosity >=3){
	       prep_declare_coef_change(row_ind, col_ind, 
					mip->colname[col_ind], 
					r_matval[a_loc], 
					mip->rhs[row_ind]);
	    }




	    
	    if(!(stats->nz_coeff_changed[a_loc])){
	       stats->nz_coeff_changed[a_loc] = TRUE;
	       stats->coeffs_changed++;
	    }

	    fix_type = IMPROVE_COEF; 
	    termcode = PREP_MODIFIED;
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
	 
	 if(a_val > 0.0){
	    if(lb[col_ind] <= -INF){

	       /* debug */
	       if(rows[row_ind].lb > -INF){
		  printf("error in prep_improve_variable()\n");
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
		  termcode = prep_modified_col_update_info(P, col_ind, 
							   new_bound,
							   IMPROVE_LB);
		  if(prep_quit(termcode)){
		     return termcode;
		  }

		  col_lb_unbounded = FALSE;		  
	       }
	    }
	    if(!col_lb_unbounded){
	       if(rows[row_ind].lb > -INF){
		  new_bound = (double)((rhs - rows[row_ind].lb + a_val*lb[col_ind])/a_val);
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
	 }else{
	    if(ub[col_ind] >= INF){

	       /* debug */
	       if(rows[row_ind].lb > -INF){
		  printf("error in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }
	       
	       col_ub_unbounded = TRUE;
	       /* can we fix it? */
	       /* is fixable if sense = 'E' and ub < INF*/
	       if(sense == 'E' && rows[row_ind].ub < INF){
		  new_bound = (double)((rhs - rows[row_ind].ub + a_val*lb[col_ind])/a_val);
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
		  }
		  termcode = prep_modified_col_update_info(P, col_ind, 
							   new_bound,
							   IMPROVE_UB);
		  if(prep_quit(termcode)){
		     return termcode;
		  }
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
      termcode = prep_modified_col_update_info(P, col_ind, new_bound, 
					       fix_type);
		
      if(prep_quit(termcode)){
	 return termcode;
      }else{
	 return PREP_MODIFIED;
      }	    
   }
   
   
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
#if 0
int prep_fix_variable(MIPdesc *mip, int col_ind, int row_ind, int a_loc, double etol){

   int fix_type, termcode = PREP_UNMODIFIED;

   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   //ROWinfo row = rows[row_ind];

   //   int *r_matbeg = mip->row_matbeg;
   // int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval; 
   
   double *ub = mip->ub;
   double *lb = mip->lb;
   char is_int = mip->is_int[col_ind];
   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];

   // int *matbeg = mip->matbeg;
   // int *matind = mip->matind;
   // double *matval = mip->matval;
   
   char fix_low_var, fix_up_var;
   double new_bound, new_ub, new_lb; 

   if(cols[col_ind].col_size <= 1){
      /* it appears on this row, so take care of it here */
      

   }
   
   if(cols[col_ind].var_type == 'B'){
      fix_low_var = FALSE;
      fix_up_var = FALSE;
      if(r_matval[a_loc] > 0){	
	 switch(sense){
	  case 'L':
	     if(rows[row_ind].lb > -INF){
		new_lb = rows[row_ind].lb + r_matval[a_loc];
		if(new_lb > rhs)
		   fix_low_var = TRUE;
	     }
	     break;
	  case 'G':
	     if(rows[row_ind].ub < INF){
		new_ub = rows[row_ind].ub - r_matval[a_loc];
		if(new_ub < rhs){
		   fix_up_var = TRUE;
		}
	     }
	     break;
	  case 'E':
	     if(rows[row_ind].lb > -INF && rows[row_ind].ub < INF){
		new_lb = rows[row_ind].lb + r_matval[a_loc];
		new_ub = rows[row_ind].ub - r_matval[a_loc];
		if(new_lb > rhs){
		   fix_low_var = TRUE;
		}
		if(new_ub < rhs){
		   fix_up_var = TRUE;
		}
		if(fix_low_var && fix_up_var){
		   stats->col_infeas_ind = col_ind;
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
	     }
	     break;
	 }
      }else{
	 
	 switch(sense){
	  case 'L':
	     if(rows[row_ind].lb > -INF){
		new_lb = rows[row_ind].lb - r_matval[a_loc]; 	 	     
		if(new_lb > rhs)
		   fix_up_var = TRUE;
	     }
	     break;
	  case 'G':
	     if(rows[row_ind].ub < INF){
		new_ub = rows[row_ind].ub + r_matval[a_loc];	 
		if(new_ub < rhs){
		   fix_low_var = TRUE;
		}
	     }
	     break;
	  case 'E':
	     if(rows[row_ind].lb > -INF && rows[row_ind].ub < INF){
		new_lb = rows[row_ind].lb - r_matval[a_loc]; 	 	     
		new_ub = rows[row_ind].ub + r_matval[a_loc];	 
		if(new_lb > rhs){
		   fix_up_var = TRUE;
		}
		if(new_ub < rhs){
		   fix_low_var = TRUE;
		}
		if(fix_low_var && fix_up_var){
		   stats->col_infeas_ind = col_ind;
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
	     }
	     break;
	 }
      }
      
      if(fix_low_var || fix_up_var){
	 if(fix_low_var){
	    new_bound = 0.0;
	    //	    ub[col_ind] = lb[col_ind] = 0.0;
	 }else{
	    new_bound = 1.0;
	    // ub[col_ind] = lb[col_ind] = 1.0;
	    // rows[row_ind].fixed_lhs_offset += r_matval[a_loc];
	 }	
	 //cols[col_ind].var_type = 'F';
	 //cols[col_ind].fix_row_ind = row_ind;
	 fix_type = FIX_BINARY;
	 //rows[row_ind].fixed_var_num++;
	 //rows[row_ind].bin_var_num--;
	 	 
	 //if(rows[row_ind].bin_var_num < 0 || rows[row_ind].fixed_var_num + rows[row_ind].fixable_var_num > rows[row_ind].size){
	    /* debug */
	 // printf("error in fixing vars 1, prep_fix_variable()\n");
	 // return PREP_OTHER_ERROR;
	 //}	 
	 termcode = PREP_MODIFIED;
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
      fix_type = FIX_OTHER; 
      //rows[row_ind].fixed_var_num++;
      //rows[row_ind].fixable_var_num--;
      //rows[row_ind].fixed_lhs_offset += r_matval[a_loc] * new_bound;

      /* debug */
      //if(rows[row_ind].fixed_var_num + rows[row_ind].fixable_var_num > rows[row_ind].size){
      // printf("error in fixing vars 2, prep_fix_variable()\n");
      // return PREP_OTHER_ERROR;
      //}	 
      termcode = PREP_MODIFIED;
   }

   /* now check if we need to update row bounds */
   if(termcode == PREP_MODIFIED){
      /* set col type to 'T', set it to F after you have visited all other rows? */
      /* have col.fix_row_ind to mark on which row you have fixed it? */      
      /* isnt worth it, update row bounds here */
      cols[col_ind].fix_row_ind = row_ind;
      termcode = prep_modified_col_update_info(P, col_ind, new_bound, 
					       fix_type);
      if(prep_quit(termcode)){
	 return termcode;
      }else
	 return PREP_MODIFIED;				 
   }
   
   return termcode;
}
#endif
/*===========================================================================*/
/*===========================================================================*/

int prep_modified_col_update_info(PREPdesc *P, int col_ind,  
				  double fixed_bound,  int fix_type)
{
   
   /* fix_type 
      0 FIX_NO_BOUND
      1 FIX_BINARY
      2 FIX_FIXABLE 
      3 FIX_OTHER
      4 IMPROVE_UB
      5 IMPROVE_LB
      6 IMPROVE_COEF
   */

   int i, r_ind, end;

   MIPdesc *mip = P->mip;
   int *matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;
   char *is_int = mip->is_int;

   //int *r_matbeg = mip->row_matbeg;
   //int *r_matind = mip->row_matind;
   //double *r_matval = mip->row_matval;
   //double * obj = mip->obj;

   double *ub = mip->ub;
   double *lb = mip->lb;

   double old_ub = ub[col_ind];
   double old_lb = lb[col_ind];
   double a_val; 
   char get_row_ubounds;
   char get_row_lbounds;

   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   prep_stats *stats = &(P->stats);


   if(fix_type != IMPROVE_LB){
      ub[col_ind] = fixed_bound;
      (stats->bounds_tightened)++;
   }
   if(fix_type != IMPROVE_UB){
      lb[col_ind] = fixed_bound;
      if(fix_type == IMPROVE_LB){
	 (stats->bounds_tightened)++;
      }
   }
   
   if(verbosity >= 3){
      printf("var %i - %s - bounds are improved:",
	     col_ind, mip->colname[col_ind]); 
      if(lb[col_ind] > -INF){
	 printf("lb:%f", fixed_bound);
      }
      if(ub[col_ind] < INF){
	 printf("ub:%f ", fixed_bound);
      }
      printf("\n");
   }

   
   if(fix_type != FIX_BINARY){
      if(prep_is_equal(ub[col_ind], lb[col_ind], etol)){
	 if(fix_type != FIX_BINARY){
	    if(cols[col_ind].var_type == 'B'){
	       fix_type = FIX_BINARY;
	    }else{
	       fix_type = FIX_OTHER;
	    }
	 }	 
	 cols[col_ind].var_type = 'F';      
      } 
   }else{
      cols[col_ind].var_type = 'F';      

   }
   
   if(cols[col_ind].var_type == 'F'){
      (stats->vars_fixed)++;
      if(verbosity >= 2){
	 prep_declare_fixed_var(col_ind, mip->colname[col_ind], ub[col_ind]);
      }
      
      mip->mip_inf->sum_obj_offset += mip->obj[col_ind]*ub[col_ind];
   }
   
   if(cols[col_ind].col_size ==0 ){
      return PREP_MODIFIED;
   }else if(cols[col_ind].col_size < 0){
      printf("error -00 in prep_fixed_col_update_info()\n");
      return PREP_OTHER_ERROR;

   }


   end = matbeg[col_ind + 1];

   for(i = matbeg[col_ind]; i < end; i++){
      get_row_ubounds = FALSE;
      get_row_lbounds = FALSE;
      if(!(rows[matind[i]].is_redundant)){
	 r_ind = matind[i];
	 a_val = matval[i];
	 if(fix_type != IMPROVE_UB && fix_type != IMPROVE_LB){
	    rows[r_ind].fixed_var_num++; 
	    if(!is_int[col_ind]){
	       rows[r_ind].cont_var_num--;
	    }
	    
	    if(!prep_is_integral(a_val, etol)){
	       rows[r_ind].frac_coef_num--;
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
	       return PREP_OTHER_ERROR;
	    }

	    rows[r_ind].fixed_lhs_offset += a_val * fixed_bound;
	 }

	 if(old_ub >= INF && fix_type != IMPROVE_LB){
	    if(a_val > 0.0){
	       rows[r_ind].ub_inf_var_num--;
	       if(rows[r_ind].ub_inf_var_num == 0){
		  get_row_ubounds = TRUE;
	       } 
	    }else{
	       rows[r_ind].lb_inf_var_num--;
	       if(rows[r_ind].lb_inf_var_num == 0){
		  get_row_ubounds = TRUE;
	       }
	    }
	 }

	 if(old_lb <= -INF && fix_type != IMPROVE_UB){
	    if(a_val > 0.0){
	       rows[r_ind].lb_inf_var_num--;
	       if(rows[r_ind].lb_inf_var_num == 0){
		  get_row_lbounds = TRUE;
	       }
	    }else{
	       rows[r_ind].ub_inf_var_num--;
	       if(rows[r_ind].lb_inf_var_num == 0){
		  get_row_lbounds = TRUE;
	       }
	    }
	 }

	 if(a_val > 0){
	    if(fix_type != IMPROVE_LB){
	       if(rows[r_ind].ub < INF){
		  /* debug */
		  if(old_ub >= INF){
		     printf("error -1 in prep_fixed_col_update_info()\n");
		     return PREP_OTHER_ERROR;
		  }
		  if(fixed_bound != old_ub){
		     rows[r_ind].ub += a_val*(fixed_bound - old_ub);
		  }
	       }
	    }
	    if(fix_type != IMPROVE_UB){
	       if(rows[r_ind].lb > -INF){
		  /* debug */
		  if(old_lb <= -INF){
		     printf("error -2 in prep_fixed_col_update_info()\n");
		     return PREP_OTHER_ERROR;
		  }
		  if(fixed_bound != old_lb){
		     rows[r_ind].lb += a_val*(fixed_bound - old_lb);
		  }
	       }
	    }
	 }else{
	    if(fix_type != IMPROVE_UB){	    
	       if(rows[r_ind].ub < INF){
		  /* debug */
		  if(old_lb <= -INF){
		     printf("error -3 in prep_fixed_col_update_info()\n");
		     return PREP_OTHER_ERROR;
		  }
		  if(fixed_bound != old_lb){
		     rows[r_ind].ub += a_val*(fixed_bound - old_lb);
		  }
	       }
	    }
	    if(fix_type != IMPROVE_LB){
	       if(rows[r_ind].lb > -INF){
		  /* debug */
		  if(old_ub >= INF){
		     printf("error -4 in prep_fixed_col_update_info()\n");
		     return PREP_OTHER_ERROR;
		  }
		  if(fixed_bound != old_ub){
		     rows[r_ind].lb += a_val*(fixed_bound - old_ub);
		  }
	       }
	    }
	 }
	 
	 /* debug */
	 if(rows[r_ind].fixed_var_num + 
	    rows[r_ind].fixable_var_num > rows[r_ind].size){
	    printf("error in fixing vars 2, prep_fix_variable()\n");
	    return PREP_OTHER_ERROR;
	 }	 
      }
      
      if(get_row_lbounds || get_row_ubounds){
	 prep_get_row_bounds(mip, r_ind);
      }
   }
   
   // mip->ub[col_ind] = mip->lb[col_ind] = fixed_bound;      
   
   /*debug */
   return PREP_MODIFIED;
}


/*===========================================================================*/
/*===========================================================================*/

int prep_get_row_bounds(MIPdesc *mip, int r_ind)
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
      if(a_val > 0.0){ 
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
      }else{
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

   double new_bound;
   
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

      rows[obj_ind].orig_ub = rows[obj_ind].sr_ub;
      rows[obj_ind].orig_lb = rows[obj_ind].sr_lb;      
      
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

	    /*find a row to be used as a constraint*/
	    for(k = last_col_loc; k < r_matbeg[obj_ind+1]; k++){
	       for(l = last_row_loc; l < c_matbeg[r_matind[k]+1]; 
		   l++){
		  if(!rows[c_matind[l]].is_redundant && !rows_checked[c_matind[l]]){
		     rows_checked[c_matind[l]] = TRUE;
		     if(rows[obj_ind].bound_type == rows[c_matind[l]].bound_type
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
		   
                   sr_solve_open_prob(sr, obj_ind, row_ind, r_matbeg, 
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
		       sr->rhs_min = sr->rhs;
		       
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
		       d_sr->rhs_min = -d_sr->rhs;
		       
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
		   
	           sr_solve_bounded_prob(sr, d_sr, obj_ind, row_ind, 
					 r_matbeg, r_matind, r_matval, 
					 cols, ub, lb, etol);
		   if(sr->sense == 'E'){
		      if(sr->ub > d_sr->ub){
			 sr->ub = d_sr->ub;
		      }
		      if(sr->lb < d_sr->lb){
			 sr->lb = d_sr->lb;
		      }		      
		   }
		   
		   sr->lb_updated = sr->ub_updated = TRUE;
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
      (*sr)->lb_var = (double *)malloc(n* DSIZE);
      (*sr)->ub_var = (double *)malloc(n* DSIZE);
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
   
   (*sr)->obj_min = (double *)malloc(DSIZE*n);
   (*sr)->matval_min = (double *)malloc(DSIZE*n);
   (*sr)->matind_min = (int *)malloc(ISIZE*n);
   (*sr)->ratio_min = (double *)malloc(DSIZE*n);
   
   
   /* debug, get something smart instead of these */
   (*sr)->tmp_ind = (int *)malloc(ISIZE*n);
   (*sr)->fixed_ind = (int *)malloc(ISIZE*n);
   
   for(k = 0; k < n; k++){
      (*sr)->fixed_ind[k] = k;
   }
}
/*===========================================================================*/
/*===========================================================================*/

int sr_solve_bounded_prob(SRdesc *sr, SRdesc *d_sr, int obj_ind, int row_ind, 
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
   
   sr_find_opt_bounded(sr, ub, lb);
   if(sr->sense == 'E'){
      sr_find_opt_bounded(d_sr, ub, lb);
   }

   return(0);

}

/*===========================================================================*/
/*===========================================================================*/

int sr_find_opt_bounded(SRdesc *sr, double *ub, double *lb)

{
   int i, col_loc;
   char max_solved = FALSE, min_solved = FALSE;
   double lhs, ax;

   int * tmp_ind = sr->tmp_ind;

   if(sr->sum_a_max <= sr->rhs_max || sr->max_n <= 0){
      sr->ub += sr->sum_c_max + sr->ub_offset;
      max_solved = TRUE;
   }
   
   if(sr->sum_a_min >= sr->rhs_min || sr->min_n <= 0){
      sr->lb += sr->sum_c_min + sr->lb_offset;
      min_solved = TRUE;
   }

   if(max_solved && min_solved){
      return 0;      
   }

   /* debug */
   /* if tmp_sum[i] = rhs, and ratios[i] > ratios[i+1], then ignore i+1 and 
      pass to i+2 */
   //hmmmm

   if(!max_solved){

      memcpy(tmp_ind, sr->fixed_ind, ISIZE*sr->max_n);
      qsort_di(sr->ratio_max, tmp_ind, sr->max_n);
      /* now fill in knapsack */
      lhs = 0;
      for(i = sr->max_n - 1; i >=0; i--){
	 col_loc = tmp_ind[i];
	 ax = sr->matval_max[col_loc] * ub[sr->matind_max[col_loc]];
	 if(lhs + ax < sr->rhs_max){
	    sr->ub += ub[sr->matind_max[col_loc]] * 
	       sr->obj_max[col_loc];	    
	    lhs += ax;
	 }else{
	    //	    ax = (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    sr->ub += sr->obj_max[col_loc] *
	       (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    break;
	 }
      }
      sr->ub += sr->ub_offset;
   }

   if(!min_solved){
      memcpy(tmp_ind, sr->fixed_ind, ISIZE*sr->min_n);

      qsort_di(sr->ratio_min, tmp_ind, sr->min_n);
      /* now fill in knapsack */
      lhs = 0;
      for(i = 0; i < sr->min_n; i++){
	 col_loc = tmp_ind[i];
	 ax = sr->matval_min[col_loc] * ub[sr->matind_min[col_loc]];
	 if(lhs + ax <= sr->rhs_min){
	    sr->lb += ub[sr->matind_min[col_loc]] * 
	       sr->obj_min[col_loc];	    
	    lhs += ax;
	 }else{
	    //	    ax = (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    sr->lb += sr->obj_min[col_loc] *
	       (sr->rhs_min - lhs)/sr->matval_min[col_loc];
	    break;
	 }
      }
      sr->lb += sr->lb_offset;
   }   
   
   return 0;

}

/*===========================================================================*/
/*===========================================================================*/

/* will add the column to problem if necessary */

int sr_add_new_col(SRdesc *sr, SRdesc *d_sr, double c_val, double a_val, 
		   int col_ind, char col_var_type, double col_ub, 
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
	 if(col_var_type != 'F'){
	    switch(sense){
	     case 'L':
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX);
		add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,	
				    col_ub, col_lb, 
				    SR_MIN);
		break;
	     case 'G':
		add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX);
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MIN);
		break;
	     case 'E':
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX);
		add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MIN);
		add_new_bounded_col(d_sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MAX);
		add_new_bounded_col(d_sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb, 
				    SR_MIN);
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
			double col_ub, double col_lb, int obj_sense)
{ 
   /* 
      ratio_type = 
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */

   //  int n;// = sr->max_n, min_n = sr->min_n; 		
 
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
   
   int *n, *matind; 
   double *obj, *matval, *rhs, *obj_offset, *sum, *obj_sum, *ratios; 

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
   }

   if(ratio_type == 0){
      obj[*n] = c_val;
      matval[*n] = a_val;
      matind[*n] = col_ind;
      ratios[*n] = c_val/a_val;
      if(obj_sense == SR_MAX){
	 *sum += rhs_ub_offset;
	 *obj_sum += obj_ub_offset;
      }else{
	 *sum += rhs_lb_offset;
	 *obj_sum += obj_lb_offset;
      }
      (*n)++;
   }else if((ratio_type == 1 && obj_sense == SR_MAX) ||
	    (ratio_type == 2 && obj_sense == SR_MIN)){
      *rhs += -rhs_ub_offset;
      *obj_offset += obj_ub_offset;
   }else if((ratio_type == 1 && obj_sense == SR_MIN) ||
	    (ratio_type == 2 && obj_sense == SR_MAX)){
      *rhs += -rhs_lb_offset;
      *obj_offset += obj_lb_offset;
   }else{
      obj[*n] = -c_val;
      matval[*n] = -a_val;
      matind[*n] = col_ind;
      ratios[*n] = c_val/a_val;
      if(obj_sense == SR_MAX){
	 *sum += -rhs_ub_offset;
	 *obj_sum += -obj_ub_offset;
      }else{
	 *sum += -rhs_lb_offset;
	 *obj_sum += -obj_lb_offset;
      }
      (*n)++;
      /* to solve bounds problem */      
      *rhs += -(rhs_ub_offset + rhs_lb_offset);
      *obj_offset += obj_ub_offset + obj_lb_offset;
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

int sr_solve_open_prob(SRdesc *sr, int obj_ind, int row_ind, int *r_matbeg, 
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

int prep_check_redundancy(PREPdesc *P, int row_ind, char use_sr_bounds)
{
   int i, termcode = PREP_UNMODIFIED, fix_type = FIX_NO_BOUND;
   int fixed_row = FALSE, col_ind;
   int debug_cnt = 0; 
   double a_val, ub, lb, new_bound, rnd_floor, rnd_ceil;

   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols; 

   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double *c_ub = mip->ub;
   double *c_lb = mip->lb;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind; 
   double *r_matval = mip->row_matval;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   prep_stats *stats = &(P->stats);

  /* 
      ratio_type = 
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */
   //   if(row.fixed_var_num + row.fixable_var_num >= row.size){
   if(rows[row_ind].fixed_var_num  >= rows[row_ind].size){
      if((sense == 'L' && rows[row_ind].fixed_lhs_offset > rhs + etol) ||
	 (sense == 'E' && 
	  !prep_is_equal(rows[row_ind].fixed_lhs_offset, rhs, etol))){
	 stats->row_infeas_ind = row_ind;
	 return PREP_INFEAS;
      }									       
      //rows[row_ind].is_redundant = TRUE;
      termcode = PREP_MODIFIED;
   }else if(rows[row_ind].fixed_var_num >= rows[row_ind].size - 1){
      for(i = r_matbeg[row_ind]; i < r_matbeg[row_ind + 1]; i++){
	 col_ind = r_matind[i];
	 a_val = r_matval[i];
	 if(cols[col_ind].var_type != 'F'){	    
	    debug_cnt++;
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
	    }else if((sense == 'G' && a_val > 0.0) ||
		     (sense == 'L' && a_val < 0.0)){
	       if(new_bound > c_ub[col_ind]){
		  stats->col_infeas_ind = col_ind;
		  stats->row_infeas_ind = row_ind;
		  return PREP_INFEAS;
	       }
	       if(new_bound > c_lb[col_ind]){
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, 
						   RND_CEIL);
		     fix_type = IMPROVE_LB;
		  }
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
		     fix_type = IMPROVE_UB;
		  }
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
	       
	       termcode = prep_modified_col_update_info(P, col_ind, new_bound, 
							fix_type);

	       if(prep_quit(termcode)){
		  return termcode;
	       }
	    }
	    
	    termcode =  PREP_MODIFIED;
	    /* debug - need to break here 
	       break; */
	 }
      }
      /* debug */
      if(debug_cnt != 1){
	 printf("error in prep_check_redundancy()\n");
	 return PREP_OTHER_ERROR;
      }
   }
   
   if(termcode == PREP_UNMODIFIED){
      if(use_sr_bounds){
	 ub = rows[row_ind].sr_ub;
	 lb = rows[row_ind].sr_lb;
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
      if(lb >= ub + etol){
	 stats->row_infeas_ind = row_ind;
	 /* debug, can this happen? */
	 return PREP_INFEAS;
      }else if (lb > ub - etol){
	 fixed_row = TRUE;
      }
      
      switch(sense){
       case 'L': 
	  if(lb > rhs + etol){
	     stats->row_infeas_ind = row_ind;
	     /* prob infeasible */
	     return PREP_INFEAS;
	  }
	  if(ub < rhs - etol || fixed_row){
	     //rows[row_ind].is_redundant = TRUE;
	     termcode = PREP_MODIFIED;
	  }
	  break;
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
       case 'E':
	  if(lb > rhs + etol || 
	     ub < rhs - etol){
	     /* prob infeasible */
	     stats->row_infeas_ind = row_ind;
	     return PREP_INFEAS;
	  }
	  if(fixed_row){
	     //rows[row_ind].is_redundant = TRUE;
	     termcode =  PREP_MODIFIED;
	  }
	  break;
       default:
	  /* debug */		
	  break;
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

   printf("col %s [%i] is fixed to %f\n",
	  name, col_ind, fixed_bound);
   return 0;

}
/*===========================================================================*/
/*===========================================================================*/
int prep_declare_coef_change(int row_ind, int col_ind, 
			     char *name, double a_val, 
			     double rhs){   
   printf("row [%i] with rhs %f: col %s [%i]: coeff improved to %f\n",
	  row_ind, rhs, name, col_ind, a_val);     
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
int prep_initialize_mipinfo(MIPdesc *mip,  prep_params params, 
			    prep_stats *stats) {
   int i, j;
   double coef_val, fixed_obj_offset;  
   int row_ind, cont_var_cnt = 0, bin_var_cnt = 0, fixed_var_cnt = 0;   
   int  row_unbounded_cnt, max_row_size, max_col_size;
   int * row_coef_bin_cnt, *row_sign_pos_cnt;   
   char is_binary, is_bounded, unbounded_below, unbounded_above;
   int gen_type; /* 0 fractional, 1 integer, 2 binary */
   int col_size, col_coef_bin_cnt, col_coef_frac_cnt, col_sign_pos_cnt; 
   int *rows_integerized_var_ind, integerizable_var_num;
   char is_opt_val_integral = TRUE;
   int is_col_all_neg; /* if we convert all constraints to 'L' */
   int is_col_all_pos; /* if we convert all constraints to 'L' */
   double etol = params.etol; 
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
   COLinfo *cols;
   ROWinfo *rows;

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
   
   for(i = 0; i < n; i++){
      is_binary = FALSE;
      is_bounded = FALSE;
      unbounded_below = FALSE;
      unbounded_above = FALSE;
      
      cols[i].var_type = 'I';
      if(lb[i] >= ub[i] + etol){
	 stats->col_infeas_ind = i;
	 return(PREP_INFEAS);
      }else if(lb[i] > (ub[i] - etol)){	 
	 cols[i].var_type = 'F';
	 fixed_obj_offset += ub[i];
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
	 if (!(coef_val-floor(coef_val) > etol &&
	     ceil(coef_val)-coef_val > etol)){
	    if((coef_val > 1.0 - etol && coef_val < 1.0 + etol) ||
	       (coef_val > -1.0 - etol && coef_val < -1.0 + etol)) {	       
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
		  printf("var %i is fixable to its upper bound %f\n", i, ub[i]);
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
		  printf("var %i is fixable to its lower bound %f\n", i, lb[i]);
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
	 prep_is_integral(rhs[j], etol) && 
	 prep_is_integral(rows[j].fixed_lhs_offset, etol)){
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
	 if (!prep_is_integral(coef_val, etol)){
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
	       printf("integerized bounds [lb-ub] of variable %s: %f - %f\n",
		      mip->colname[i],lb[i],ub[i]);
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
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_fill_row_ordered(MIPdesc *mip)
{ 
   /*
     recreates 'A' matrix using three matrices just like the standard
     notation. However, matrices contain row representations rather than
     column
   */ 

   int i, j;
   int row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg, *r_lengths; 
   double * matval, *r_matval, *rhs;
   int n = mip->n;
   int m = mip->m;
   int nz = mip->nz;
   char *sense, *o_sense;
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

   /* first get row legths */   
   for(i = 0; i < n; i++){
      for(j = matbeg[i]; j < matbeg[i+1]; j++){
	 r_lengths[matind[j]]++;
      }
   }

   r_matbeg[0] = 0;

   /* fill in matbegs */
   for(i = 0; i < m; i++){
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
int prep_update_mip(MIPdesc *mip, prep_stats stats, prep_params params, 
		    double etol)
{ 
   
   int i, j, col_nz, col_num, row_num;
   int row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg, *r_lengths; 
   double *ub, *lb, *matval, *r_matval, *obj, *rhs, obj_offset;
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

   int deleted_row_cnt = stats.rows_deleted;
   int vars_fixed = stats.vars_fixed;
   int keep_row_ordered = params.keep_row_ordered;
   int reduce_mip = params.reduce_mip;

   int *row_new_inds = NULL;
      
   /* debug */
   /* ---- */
   double *fixed_lhs_offset = NULL; 
   int *row_lengths = NULL;
   
   obj_offset = mip->mip_inf->sum_obj_offset;

   if(keep_row_ordered && (!reduce_mip || 
			   (deleted_row_cnt <= 0 && vars_fixed))){
      return PREP_UNMODIFIED;      
   }


   sense = mip->sense;
   o_sense = mip->orig_sense;
   rhs = mip->rhs;
   
   if(reduce_mip){
      
      mip->alloc_n = n;
      mip->alloc_m = m;
      mip->alloc_nz = nz;
      
      matval = mip->matval;
      matind = mip->matind;
      matbeg = mip->matbeg;
                 
      colnames = mip->colname;
      obj = mip->obj;
      ub = mip->ub;
      lb = mip->lb;
      is_int = mip->is_int;
      
      int *row_new_inds = (int *)malloc(m*ISIZE);
      
      /* debug */
      /* ---- */
      double *fixed_lhs_offset = (double *)malloc(m*DSIZE);
      int *row_lengths = (int *)malloc(m*ISIZE);
      for(i= 0; i < m; i++){
	 fixed_lhs_offset[i] = rows[i].fixed_lhs_offset;
	 row_lengths[i] = rows[i].size - rows[i].fixed_var_num;
	 rows[i].fixed_lhs_offset = 0.0;
	 rows[i].size = 0;
      }
      /* ---- */
      
      /* first get new row indices */
      
      row_num = 0;
      
      if(deleted_row_cnt > 0){
	 for(i = 0; i < m; i++){
	    if(!(rows[i].is_redundant)){
	       row_new_inds[i] = row_num;
	       row_num++;
	    }
	 }
	 
	 /* debug */
	 if(row_num != m - deleted_row_cnt){
	    printf("error: missing rows \n");
	    return PREP_OTHER_ERROR;
	 }
      }else{
	 for(i = 0; i < m; i++){
	    row_new_inds[i] = i;
	    /* debug */
	    if(!(rows[i].is_redundant)){
	       printf("error: missing redundants \n");
	       return PREP_OTHER_ERROR;
	    }
	 }
	 row_num = m;
      }
      

      /* first fill in col ordered */
      /* debug  - fixme: original senses, need to convert nonzeros too*/
      col_nz = col_num = 0;
      obj_offset = 0;
      for(i = 0; i < n; i++){
	 if(cols[i].var_type != 'F'){
	    for(j = matbeg[i]; j < matbeg[i+1]; j++){
	       r_ind = matind[j];
	       if(!(rows[r_ind].is_redundant)){
		  if(!prep_is_equal(matval[j], 0.0, etol)){
		     matind[col_nz] = row_new_inds[r_ind];
		     matval[col_nz] = o_sense[r_ind] == 'G' ? -matval[j] : 
			matval[j];
		     col_nz++;
		     (rows[r_ind].size)++;
		  }
	       }
	    }
	    
	    obj[col_num] = obj[i];
	    ub[col_num] = ub[i];
	    lb[col_num] = lb[i];
	    is_int[col_num] = is_int[i];
	    cols[col_num] = cols[i];
	    cols[col_num].col_size = col_nz - matbeg[col_num];
	    strcpy(colnames[col_num], colnames[i]);
	    matbeg[(++col_num)] = col_nz;
	    /* debug */
	    if(cols[col_num - 1].col_size <= 0){
	       printf("error: empty size column \n");
	       return PREP_OTHER_ERROR;
	    }
	 }else{
	    /* debug */
	    /*---- -*/
	    if(!prep_is_equal(ub[i], lb[i], etol)){
	       printf("error: not fixed column? \n");
	       return PREP_OTHER_ERROR;
	    }
	    /* ----- */
	    
	    obj_offset += ub[i] * obj[i];
	    /* debug */
	    /* ---- */
	    for(j = matbeg[i]; j < matbeg[i+1]; j++){
	       r_ind = matind[j];
	       if(!(rows[r_ind].is_redundant)){
		  rows[r_ind].fixed_lhs_offset += ub[i]* matval[j];	       
	       }
	    }
	    /* ---- */
	 }
      }
      
      /* debug */
      /* ----- */
      if(col_num != n - vars_fixed){
	 printf("error: missing cols \n");
	 return PREP_OTHER_ERROR;
      }
      
      for(i= 0; i < m; i++){
	 if(!(rows[i].is_redundant)){
	    if(!prep_is_equal(fixed_lhs_offset[i], rows[i].fixed_lhs_offset, 
			      etol)){
	       printf("error: missing offsets\n");
	       return PREP_OTHER_ERROR;
	    }
	    //if(coeffs_nulled <= 0){
	    //  if(rows[i].size != row_lengths[i]){
	    //     printf("error: missing row size\n");
	    //     return PREP_OTHER_ERROR;
	    //  }
	    //}
	 }
      }   
      /*----*/
 
      /* now convert it row ordered if asked*/
      /* get row_lengths and fill in r_matbeg, sense, rhs etc */

      if(keep_row_ordered){

	 r_matbeg = mip->row_matbeg;
	 r_matind = mip->row_matind;
	 r_matval = mip->row_matval;
	 r_lengths = mip->row_lengths;      
	 
	 row_num = 0;
	 for(i = 0; i < m; i++){
	    if(!(rows[i].is_redundant)){
	       rows[row_num] = rows[i];
	       r_lengths[row_num] = rows[row_num].size;
	       r_matbeg[row_num + 1] = r_matbeg[row_num] + 
		  r_lengths[row_num];	 
	       sense[row_num] = o_sense[row_num] = o_sense[i];
	       rhs[row_num] -= rows[row_num].fixed_lhs_offset;
	       rhs[row_num] = sense[row_num] == 'G'? -rhs[i] : rhs[i];
	       row_num++;
	    }	 
	 }
	 
      /* debug */
	 if(r_matbeg[row_num ] != col_nz){
	    printf("error; missing nonzeros\n");
	    return PREP_OTHER_ERROR;
	 }
	 
	 /* */
	 
	 for(i = 0; i < n; i++){
	    for(j = matbeg[i]; j < matbeg[i+1]; j++){
	       row_ind = matind[j];
	       elem_ind = r_matbeg[row_ind];
	       r_matind[elem_ind] = i;
	       r_matval[elem_ind] = matval[j];
	       r_matbeg[row_ind] = elem_ind + 1;
	    }
	 }
	 
	 for(i = 0; i < m; i++){
	    r_matbeg[i] -= r_lengths[i];
	 }
      }
   }else{
      /* just update rows, sense and rhs and delete row ordered */

      /* and get fixed offset value */
      for(i = 0; i < n; i++){
	 if(cols[i].var_type == 'F'){
	    obj_offset += ub[i] * obj[i];
	    //for(j = matbeg[i]; j < matbeg[i+1]; j++){
	    //  r_ind = matind[j];
	    //  if(!(rows[r_ind].is_redundant)){
	    //	  rows[r_ind].fixed_lhs_offset;	       
	    //  }
	 }
      }
      
      for(i = 0; i < m; i++){
	 if(!(rows[i].is_redundant)){
	    rows[row_num] = rows[i];
	    sense[row_num] = o_sense[i];
	    rhs[row_num] -= rows[row_num].fixed_lhs_offset;
	    rhs[row_num] = sense[row_num] == 'G'? -rhs[i] : rhs[i];
	    row_num++;
	 }	 
      }      

      if(!keep_row_ordered){
	 FREE(mip->row_matbeg);
	 FREE(mip->row_matind);
	 FREE(mip->row_matval);
	 FREE(mip->row_lengths);
	 FREE(mip->orig_sense);
      }
   }
 
   mip->obj_offset = obj_offset;
     
   if(reduce_mip){
      mip->n = col_num;
      mip->m = row_num;
      mip->nz = col_nz;

      /* debug */
      if(!prep_is_equal(mip->obj_offset, obj_offset, etol)){
	 printf("error; deviance on obj fixed offset\n");
	 return PREP_OTHER_ERROR;
      }
       
      FREE(row_new_inds);
      FREE(fixed_lhs_offset);
      FREE(row_lengths);
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
   int verbosity = P->params.verbosity;
   prep_stats stats = P->stats;

   if(verbosity >= 1){
      switch(termcode){
       case PREP_INFEAS:
	  printf("Preprocessing detected infeasibility...");
	  if(stats.col_infeas_ind >= 0 ||
	     stats.row_infeas_ind >= 0){
	     printf("while improving bounds of \n\t");	  
	     if(stats.col_infeas_ind >= 0){
		printf("variable %s [%i] ", 
		       mip->colname[stats.col_infeas_ind], 
		       stats.col_infeas_ind);
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
	     printf("variable %s [%i]\n", 
		    mip->colname[stats.col_unbound_ind], 
		    stats.col_unbound_ind);
	  }
	  break;
       case PREP_NUMERIC_ERROR:
	  printf("Preprocessing detected numerical problems ");	  
	  if(stats.col_numeric_ind >= 0){
	     printf("while improving bounds on \n");	  
	     printf("variable %s [%i]\n", 
		    mip->colname[stats.col_numeric_ind], 
		    stats.col_numeric_ind);
	  }
	  break;
       case PREP_OTHER_ERROR:
	  printf("Preprocessing - unknown error.. ignoring presolve...\n");
	  break;
       case PREP_SOLVED:
	  printf("Preprocessing found the optimum:");	  	  
	  /* fixme declare solution */
	  break;	  
       default:
	  printf("Preprocessing finished...\n ");	  

	  if(stats.coeffs_changed + 
	     stats.bounds_tightened + 
	     stats.rows_deleted + 
	     stats.vars_fixed > 0){
	     if(stats.coeffs_changed > 0){
		printf("\t modified %i coefficients \n",
		       stats.coeffs_changed);		       
	     }
	     if(stats.bounds_tightened > 0){
		printf("\t improved %i bounds \n", 
		       stats.bounds_tightened);
	     }	     
	     if(stats.rows_deleted + 
		stats.vars_fixed > 0){
		if(stats.rows_deleted > 0){
		   printf("\t removed %i constraints\n", 
			  stats.rows_deleted);
		   //printf("\t %i remained\n", mip->m);
		}
		if(stats.vars_fixed > 0){
		   printf("\t fixed %i variables\n", stats.vars_fixed);
		   //printf("\t %i remained\n", mip->n);
		}	     
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
   }
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
      FREE(P->stats.nz_coeff_changed);
      if(P->sr){
	 prep_free_sr_desc(P->sr);
      }
      if(P->d_sr){
	 prep_free_sr_desc(P->d_sr);
      }
      
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
      FREE(sr->lb_var);
      FREE(sr->ub_var);

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
