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
/* This file was written by Ashutosh (asm4@lehigh.edu)                       */
/*===========================================================================*/

#define ERR_TOL 1e-6		/* error tolerance, should really be a
				   parameter set based on the processor and
				   input method; unused so far */
/*===========================================================================*/
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
/*
#include "preprocessor_constants.h"
#include "symphony_api.h"
#include "proccomm.h" 
#include "master.h" 
#include "BB_constants.h" 
*/

/*===========================================================================*/
/*===========================================================================*/
int prep_basic(PREPdesc *P, char keep_track)
{

   /*
     This function is the master of the preprocessing part. It calls and 
     controls other functions to perform preprocessing jobs.
   */

   int termstatus;		/* return status of each function called
				   herein */
   int can_iterate = TRUE, iter_cnt = 0, iter_cnt_limit, p_level; 
   int redundant_row_cnt = 0, fixed_col_cnt = 0;
   
   /* variables which collect stats*/
   //int bounds_tightened = 0, vars_fixed=0, coeffs_changed=0, 
   int verbosity;
   double etol;// min_ub, max_lb; 
   char do_sr_rlx, do_aggr_row_rlx;// fix_var; 
   int i, j, m, n, nz, *r_matbeg, *r_matind, max_size, min_size, *max_ind, *min_ind; 
   double *obj, *rhs, *r_matval, *ub, *lb;// old_val, old_lb, old_ub;  
   char *sense; 
   int * updated_cols_ind, updated_cols_cnt, * updated_rows_ind, updated_rows_cnt;
   int col_ind, row_ind;
   int * modified_cols_ind, * modified_rows_ind;
   //   int modified_cols_cnt, modified_rows_cnt;
   char *is_col_updated, *is_row_updated;
   char var_type;// sr_termcode; 
   MIPdesc *mip = P->mip;
   MIPinfo *mip_inf = mip->mip_inf;
   COLinfo * cols = mip_inf->cols;
   ROWinfo *row, *rows = mip_inf->rows;
   prep_stats stats = P->stats;
   prep_params params = P->params;
   
   params.prep_verbosity=2;

   verbosity = params.prep_verbosity; /* debug */
   p_level = params.prep_level;
   etol = params.etol;
   iter_cnt_limit = params.iteration_limit;
   do_sr_rlx = params.do_single_row_rlx;      
   do_aggr_row_rlx = params.do_aggregate_row_rlx;
   
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
 

   /* main preprocessing loop */

   while(iter_cnt < iter_cnt_limit && updated_rows_cnt > 0){
      
      iter_cnt++;
      modified_rows_ind = 0;
      modified_cols_ind = 0;

      if (verbosity>=2) {
	 printf("Basic iteration number: %d\n", iter_cnt);
      }

      /* check the updated bounds and cols to iterate on*/
      /* for each updated row and column do bla bla...*/
      /* while iterating update the new cols and rows indices*/
      
      /* do this only once at the end of first iteration*/
      //if(prep_level > 2 && (do_sr_rlx || do_aggr_row_rlx) && 
      // iter_cnt == 1){
      /*===========================================================================*/
      /*===========================================================================*/
      
      for(row_ind = 0; row_ind < m; row_ind++){

	 row = &(rows[row_ind]);

	 if(row->is_redundant){
	    continue;
	 }

	 can_iterate = FALSE;
	 if(row->ub >= DBL_MAX && row->lb <= -DBL_MAX){
	    //    if(prep_solve_sr_rlx(P, 1, &row_ind) != SR_BOUNDS_UPDATED){
	       continue;
	       // }else{
	       // can_iterate = TRUE;
	       //  }
	    /* call sr to see if this can be changed */
	 }else{
	    can_iterate = TRUE;
	 }
	 
	 if(can_iterate){
	    
	    if(row->type != CONTINUOUS_TYPE &&
	       row->type != ALL_MIXED_TYPE && 
	       row->coef_type != FRACTIONAL_VEC){
	       row->lb = prep_rnd_integral(row->lb, etol, RND_CEIL);
	       row->ub = prep_rnd_integral(row->ub, etol, RND_FLOOR);
	    }
	    
	    if(prep_check_redundancy(mip, row_ind, etol, FALSE)
	       == PREP_INFEAS){	       
				     
	       if(verbosity >= 1){
		  printf("infeasibility detected\n"); 
		  printf("row %i is redundant: ub: %f \t lb: %f \t sense:" 
			 "%c \t rhs: %f \n", 
			 row_ind, row->ub >= DBL_MAX ? -1 : row->ub, 
			 row->lb <= -DBL_MAX ? 1 : row->lb, sense[row_ind], 
			 rhs[row_ind]);
	       }
	       /* debug*/
	       exit(0);
	       return PREP_INFEAS;
	    }
	    if(row->is_redundant){
	       redundant_row_cnt++;
	       if(verbosity >= 1){
		  printf("row %i is redundant: ub: %f \t lb: %f \t sense: %c \t rhs: %f \n", 
			 row_ind, row->ub >= DBL_MAX ? -1 : row->ub, 
			 row->lb <= -DBL_MAX ? 1 : row->lb, sense[row_ind], 
			 rhs[row_ind]);
	       }
	       continue;
	    }
	 }
	 if(can_iterate){/* || (!(can_iterate) && 
			    (row.ub_inf_var_num < 1 || 
			    row.lb_inf_var_num < 1))){*/

	    //row_ub = row.sr_ub;
	    //row_lb = row.sr_lb;
	    //new_ub = 0;
	    //new_lb = 0;
	    // prep_solve_sr_rlx(P, 1, &row_ind);

	    for(j = r_matbeg[row_ind]; j < r_matbeg[row_ind+1]; j++){
	       
	       /* bounds improvement, variable fixing, coeff improvement all here */
	       /* while checking bounds improvement, update bounds */
	       
	       /* when an unbounded column is bounded, keep the info that 
		  which rows' bounds should be updated */
	       col_ind = r_matind[j];
	       var_type = cols[col_ind].var_type;
	       if(var_type != 'F'){
		  if(prep_fix_variable(mip, col_ind, row_ind, j, etol) == 
		     //if(prep_improve_variable(mip, col_ind, row_ind, j, etol) == 
		     PREP_MODIFIED){
		     fixed_col_cnt++;
		     printf("var %i is fixed\n", col_ind);
		  }
	       }	    
	    }
	    


	 /*===========================================================================*/
	 /*===========================================================================*/
	    if (verbosity>=10) {
	       printf("Checking redundancy\n");
	    }
	    
	 
	 
	  /*===========================================================================*/
	  /*===========================================================================*/
	    if (verbosity>=10) {
	       printf("Fixing variables\n");
	    } 
	    
	 




	 /*===========================================================================*/
	 /*===========================================================================*/
	    if (verbosity>=10) {
	       printf("Improving Binary Coefficients\n");
	    } 

#if 0	 
	    if(do_sr_rlx){
	       if(iter_cnt <= 2){
		  sr_termcode = prep_solve_sr_rlx(P, 1, &row_ind);
	       }
	    }
#endif
	 /* update bounds and bla bla */	 
	 //      }

	    if(iter_cnt == iter_cnt_limit - 1){
	    /* do gcd thing and see for infeasiblity and further bound 
	       tightening:  ax <= b
	       gcd(a) = p, make rhs = kp <= b
	    */
	    }
	 
	 }
      
	 /* one last time to detect more redundancy */
#if 0
	 if(do_sr_rlx){
	    for(row_ind = 0; row_ind < m; row_ind++){
	       
	       if(row->is_redundant){
		  continue;
	       }
	       sr_termcode = prep_solve_sr_rlx(P, 1, &row_ind);
	       if(sr_termcode == SR_BOUNDS_UPDATED){
		  prep_check_redundancy(rows, row_ind, sense[row_ind], rhs[row_ind], 
					cols, r_matind, r_matbeg, r_matval, ub, lb, 
					etol, TRUE);
		  if(row->is_redundant){
		     redundant_row_cnt++;
		     if(verbosity >= 1){
			printf("row %i is redundant: ub: %f \t lb: %f \t sense: %c \t rhs: %f \n", 
			       row_ind, row->ub >= DBL_MAX ? -1 : row->ub, 
			       row->lb <= -DBL_MAX ? 1 : row->lb, sense[row_ind], 
			       rhs[row_ind]);
		     }
		  }
	       }
	    }
	 }
#endif 
      }
   }
      
   if (verbosity>=1) {
      printf("Redundant row cnt: %i\n", redundant_row_cnt);
      printf("Fixed col cnt: %i\n", fixed_col_cnt);
   }

   
#if 0
   while (iterate_more) {
      iterate_more = 0;
      iteration_number++;
      if (verbosity>=1) {
	 printf("Basic iteration number: %d\n", iteration_number);
      }

      if (verbosity>=1) {
	 printf("Tightening bounds\n");
      }
      vars_fixed = 0;
      //termstatus = 0;
      termstatus = prep_tighten_bounds(P, row_P, lhs, -1, bounds_tightened, 
          vars_fixed, verbosity);
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 stats->bounds_tightened = stats->bounds_tightened+bounds_tightened;
	 stats->vars_fixed = stats->vars_fixed+vars_fixed;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (verbosity>=1) {
	 printf("Checking redundancy\n");
      }
      //termstatus = 0;
      termstatus = prep_check_redundancy(P, row_P, lhs, verbosity); 
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (verbosity>=1) {
	 printf("Fixing variables\n");
      }      
      vars_fixed = 0;
      termstatus = 0;
      termstatus = prep_fix_variables(P, row_P, lhs, vars_fixed, verbosity);
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 stats->vars_fixed = stats->vars_fixed+vars_fixed;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (verbosity>=1) {
	 printf("Improving Binary Coefficients\n");
      } 
      coeffs_changed = 0;
      termstatus = 0;
      termstatus = prep_BinImproveCoeffs(P, row_P, lhs, coeffs_changed,
				 verbosity);
      switch (termstatus) {
       case PREP_MODIFIED:
	 iterate_more = 1;
	 stats->coeffs_changed = stats->coeffs_changed+coeffs_changed;
	 break;
       case PREP_INFEAS:
	 return PREP_INFEAS;
	 break;
       case PREP_UNBOUNDED:
	 return PREP_UNBOUNDED;
	 break;
       case PREP_SOLVED:
	 termstatus = prep_declare_solved(P, verbosity);
	 return PREP_SOLVED;
	 break;
      }

      if (P->m == 0) {
	 termstatus = prep_declare_solved(P, verbosity);
      }
   }

#ifdef DISPLAY_MIP
   termstatus = prep_display_mip(current_mip);
#endif
   
   if (iteration_number>1) {
      return PREP_MODIFIED;
   }
   else {
      return PREP_UNMODIFIED;
   }

#endif
   return termstatus;

   /* exit basic preprocessor */
}


/*===========================================================================*/
/*===========================================================================*/
#if 0
int prep_improve_variable(MIPdesc *mip, int col_ind, int row_ind, int a_loc, 
			  double etol){

   int fix_type, termcode = PREP_UNMODIFIED;

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
   
   char fix_low_var, fix_up_var;
   double new_bound, new_ub, new_lb; 

   if(cols[col_ind].col_size <= 1){

      ratio_type = prep_get_ratio_type(obj[col_ind], a_val);

      /* 
	 ratio_type = 
	 0 c_val >0, a_val>0
	 1 c_val >= 0, a_val <= 0
	 2 c_val <= 0, a_val >= 0
	 3 c_val < 0, a_val < 0
      */
      
      switch(sense){
       case 'L': 
	  if(ratio_type == 0){
	     

	  }
      



      if(obj[col_ind] > 0.0){
	 if(a_val > 0.0){
	    
	 }


      }




      switch(sense){
       case 'L':
	  if(a_val > 0.0){
	     
	     


	  }

      }



      new_bound = (double)((rhs - rows[row_ind].fixed_lhs_offset)/(r_matval[i]));
      if(new_bound > c_ub[col_ind] || new_bound < c_lb[col_ind]){
	 return PREP_INFEAS;
      }
	    if(cols[col_ind].var_type != 'C'){
	       if(new_bound > floor(new_bound) + etol &&
		  new_bound < ceil(new_bound) - etol){
		  return PREP_INFEAS;
	       }
	    }
	    
	    if(cols[col_ind].var_type == 'B'){
	       fix_type = FIX_BINARY;
	    }else{
	       fix_type = FIX_OTHER;
	    }

	    //c_ub[col_ind] = c_lb[col_ind] = new_bound;
	    //cols[col_ind].var_type = 'F';	       

	    //cols[col_ind].fix_row_ind = row_ind;
	    //rows[row_ind].fixed_var_num++;
	    //if(cols[col_ind].var_type == 'B'){
	    //   rows[row_ind].bin_var_num--;
	    // }
	    
	    rows[row_ind].is_redundant = TRUE;
	    prep_fixed_col_update_info(mip, col_ind, new_bound, fix_type);

	    termcode =  PREP_MODIFIED;
	    debug_cnt++;
      }
      



      








      
      
      
      
      if(a_val > 0.0){
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
   }

   if(cols[col_ind].var_type == 'B'){
      fix_low_var = FALSE;
      fix_up_var = FALSE;
      if(r_matval[a_loc] > 0){	
	 switch(sense){
	  case 'L':
	     if(rows[row_ind].lb > -DBL_MAX){
		new_lb = rows[row_ind].lb + r_matval[a_loc];
		if(new_lb > rhs)
		   fix_low_var = TRUE;
	     }
	     break;
	  case 'G':
	     if(rows[row_ind].ub < DBL_MAX){
		new_ub = rows[row_ind].ub - r_matval[a_loc];
		if(new_ub < rhs){
		   fix_up_var = TRUE;
		}
	     }
	     break;
	  case 'E':
	     if(rows[row_ind].lb > -DBL_MAX && rows[row_ind].ub < DBL_MAX){
		new_lb = rows[row_ind].lb + r_matval[a_loc];
		new_ub = rows[row_ind].ub - r_matval[a_loc];
		if(new_lb > rhs){
		   fix_low_var = TRUE;
		}
		if(new_ub < rhs){
		   fix_up_var = TRUE;
		}
		if(fix_low_var && fix_up_var){
		   return PREP_INFEAS;
		}
	     }
	     break;
	 }
      }else{

	 switch(sense){
	  case 'L':
	     if(rows[row_ind].lb > -DBL_MAX){
		new_lb = rows[row_ind].lb - r_matval[a_loc]; 	 	     
		if(new_lb > rhs)
		   fix_up_var = TRUE;
	     }
	     break;
	  case 'G':
	     if(rows[row_ind].ub < DBL_MAX){
		new_ub = rows[row_ind].ub + r_matval[a_loc];	 
		if(new_ub < rhs){
		   fix_low_var = TRUE;
		}
	     }
	     break;
	  case 'E':
	     if(rows[row_ind].lb > -DBL_MAX && rows[row_ind].ub < DBL_MAX){
		new_lb = rows[row_ind].lb - r_matval[a_loc]; 	 	     
		new_ub = rows[row_ind].ub + r_matval[a_loc];	 
		if(new_lb > rhs){
		   fix_up_var = TRUE;
		}
		if(new_ub < rhs){
		   fix_low_var = TRUE;
		}
		if(fix_low_var && fix_up_var){
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
	 // exit(0);
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
      fix_type = FIX_FIXABLE; 
      //rows[row_ind].fixed_var_num++;
      //rows[row_ind].fixable_var_num--;
      //rows[row_ind].fixed_lhs_offset += r_matval[a_loc] * new_bound;

      /* debug */
      //if(rows[row_ind].fixed_var_num + rows[row_ind].fixable_var_num > rows[row_ind].size){
      // printf("error in fixing vars 2, prep_fix_variable()\n");
      // exit(0);
      //}	 
      termcode = PREP_MODIFIED;
   }

   /* now check if we need to update row bounds */
   if(termcode == PREP_MODIFIED){
      /* set col type to 'T', set it to F after you have visited all other rows? */
      /* have col.fix_row_ind to mark on which row you have fixed it? */      
      /* isnt worth it, update row bounds here */
      cols[col_ind].fix_row_ind = row_ind;
      prep_fixed_col_update_info(mip, col_ind, new_bound, 
				 fix_type);
   }


   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
#endif

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
	     if(rows[row_ind].lb > -DBL_MAX){
		new_lb = rows[row_ind].lb + r_matval[a_loc];
		if(new_lb > rhs)
		   fix_low_var = TRUE;
	     }
	     break;
	  case 'G':
	     if(rows[row_ind].ub < DBL_MAX){
		new_ub = rows[row_ind].ub - r_matval[a_loc];
		if(new_ub < rhs){
		   fix_up_var = TRUE;
		}
	     }
	     break;
	  case 'E':
	     if(rows[row_ind].lb > -DBL_MAX && rows[row_ind].ub < DBL_MAX){
		new_lb = rows[row_ind].lb + r_matval[a_loc];
		new_ub = rows[row_ind].ub - r_matval[a_loc];
		if(new_lb > rhs){
		   fix_low_var = TRUE;
		}
		if(new_ub < rhs){
		   fix_up_var = TRUE;
		}
		if(fix_low_var && fix_up_var){
		   return PREP_INFEAS;
		}
	     }
	     break;
	 }
      }else{

	 switch(sense){
	  case 'L':
	     if(rows[row_ind].lb > -DBL_MAX){
		new_lb = rows[row_ind].lb - r_matval[a_loc]; 	 	     
		if(new_lb > rhs)
		   fix_up_var = TRUE;
	     }
	     break;
	  case 'G':
	     if(rows[row_ind].ub < DBL_MAX){
		new_ub = rows[row_ind].ub + r_matval[a_loc];	 
		if(new_ub < rhs){
		   fix_low_var = TRUE;
		}
	     }
	     break;
	  case 'E':
	     if(rows[row_ind].lb > -DBL_MAX && rows[row_ind].ub < DBL_MAX){
		new_lb = rows[row_ind].lb - r_matval[a_loc]; 	 	     
		new_ub = rows[row_ind].ub + r_matval[a_loc];	 
		if(new_lb > rhs){
		   fix_up_var = TRUE;
		}
		if(new_ub < rhs){
		   fix_low_var = TRUE;
		}
		if(fix_low_var && fix_up_var){
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
	 // exit(0);
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
      fix_type = FIX_FIXABLE; 
      //rows[row_ind].fixed_var_num++;
      //rows[row_ind].fixable_var_num--;
      //rows[row_ind].fixed_lhs_offset += r_matval[a_loc] * new_bound;

      /* debug */
      //if(rows[row_ind].fixed_var_num + rows[row_ind].fixable_var_num > rows[row_ind].size){
      // printf("error in fixing vars 2, prep_fix_variable()\n");
      // exit(0);
      //}	 
      termcode = PREP_MODIFIED;
   }

   /* now check if we need to update row bounds */
   if(termcode == PREP_MODIFIED){
      /* set col type to 'T', set it to F after you have visited all other rows? */
      /* have col.fix_row_ind to mark on which row you have fixed it? */      
      /* isnt worth it, update row bounds here */
      cols[col_ind].fix_row_ind = row_ind;
      prep_fixed_col_update_info(mip, col_ind, new_bound, 
				 fix_type);
   }


   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

void prep_fixed_col_update_info(MIPdesc *mip, int col_ind,  
				double fixed_bound, int fix_type)
{

   /* fix_type 
      0 FIX_NO_BOUND
      1 FIX_BINARY
      2 FIX_FIXABLE 
      3 FIX_OTHER
      4 FIX_UB
      5 FIX_LB
   */

   int i, j, c_ind, r_ind, end;
   int *matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;
   double * obj = mip->obj;

   double *ub = mip->ub;
   double *lb = mip->lb;

   double old_ub = ub[col_ind];
   double old_lb = lb[col_ind];
   double a_val; 
   char get_row_ubounds;
   char get_row_lbounds;

   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;

   if(fix_type != FIX_LB){
      ub[col_ind] = fixed_bound;
   }
   if(fix_type != FIX_UB){
      lb[col_ind] = fixed_bound;
   }
   
   if(fix_type != FIX_LB && fix_type != FIX_UB){
      cols[col_ind].var_type = 'F';
   }

   end = matbeg[col_ind + 1];

   for(i = matbeg[col_ind]; i < end; i++){
      get_row_ubounds = FALSE;
      get_row_lbounds = FALSE;
      if(!(rows[matind[i]].is_redundant)){
	 r_ind = matind[i];
	 a_val = matval[i];
	 if(fix_type != FIX_UB && fix_type != FIX_LB){
	    rows[r_ind].fixed_var_num++; 

	    if(fix_type == FIX_BINARY){
	       rows[r_ind].bin_var_num--;
	    }else if(fix_type == FIX_FIXABLE){
	       rows[r_ind].fixable_var_num--;
	    }
	    /* debug */
	    if(rows[r_ind].bin_var_num < 0 || 
	       rows[r_ind].fixable_var_num < 0){
	       printf("error -0 in prep_fixed_col_update_info()\n");
	       exit(0);
	    }

	    rows[r_ind].fixed_lhs_offset += a_val * fixed_bound;
	 }

	 if(old_ub >= DBL_MAX && fix_type != FIX_LB){
	    rows[r_ind].ub_inf_var_num--;
	    if(old_lb <= -DBL_MAX){
	       rows[r_ind].free_var_num--;
	    }
	    if(rows[r_ind].ub_inf_var_num == 0 && rows[r_ind].ub >= DBL_MAX ){
	       get_row_ubounds = TRUE;
	    }
	 }

	 if(old_lb <= -DBL_MAX && fix_type != FIX_UB){
	    rows[r_ind].lb_inf_var_num--;
	    if(old_ub >= DBL_MAX){
	       rows[r_ind].free_var_num--;
	    }
	    if(rows[r_ind].lb_inf_var_num == 0 && rows[r_ind].lb <= -DBL_MAX ){
	       get_row_lbounds = TRUE;
	    }
	 }

	 if(a_val > 0){
	    if(fix_type != FIX_LB){
	       if(rows[r_ind].ub < DBL_MAX){
		  /* debug */
		  if(old_ub >= DBL_MAX){
		     printf("error -1 in prep_fixed_col_update_info()\n");
		     exit(0);
		  }
		  if(fixed_bound != old_ub){
		     rows[r_ind].ub += a_val*(fixed_bound - old_ub);
		  }
	       }
	    }
	    if(fix_type != FIX_UB){
	       if(rows[r_ind].lb > -DBL_MAX){
		  /* debug */
		  if(old_lb <= -DBL_MAX){
		     printf("error -2 in prep_fixed_col_update_info()\n");
		     exit(0);
		  }
		  if(fixed_bound != old_lb){
		     rows[r_ind].lb += a_val*(fixed_bound - old_lb);
		  }
	       }
	    }
	 }else{
	    if(fix_type != FIX_UB){	    
	       if(rows[r_ind].ub < DBL_MAX){
		  /* debug */
		  if(old_lb <= -DBL_MAX){
		     printf("error -3 in prep_fixed_col_update_info()\n");
		     exit(0);
		  }
		  if(fixed_bound != old_lb){
		     rows[r_ind].ub += a_val*(fixed_bound - old_lb);
		  }
	       }
	    }
	    if(fix_type != FIX_LB){
	       if(rows[r_ind].lb > -DBL_MAX){
		  /* debug */
		  if(old_ub >= DBL_MAX){
		     printf("error -4 in prep_fixed_col_update_info()\n");
		     exit(0);
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
	    exit(0);
	 }	 
      }

      if(get_row_lbounds || get_row_ubounds){
	 rows[r_ind].ub = rows[r_ind].lb = 0.0;
	 for(j = r_matbeg[r_ind]; j < r_matbeg[r_ind + 1]; j++){
	    a_val = r_matval[j];
	    c_ind = r_matind[j];
	    if(a_val > 0.0){ 
	       if(rows[r_ind].ub < DBL_MAX){
		  if(ub[c_ind] >= DBL_MAX){
		     rows[r_ind].ub = DBL_MAX;
		  }else{
		     rows[r_ind].ub += a_val * ub[c_ind];
		  }
	       }
	       if(rows[r_ind].lb > -DBL_MAX){
		  if(lb[c_ind] <= -DBL_MAX){
		     rows[r_ind].lb = -DBL_MAX;
		  }else{
		     rows[r_ind].lb += a_val * lb[c_ind];
		  }
	       }
	    }else{
	       if(rows[r_ind].ub < DBL_MAX){
		  if(lb[c_ind] <= -DBL_MAX){
		     rows[r_ind].ub = DBL_MAX;
		  }else{
		     rows[r_ind].ub += a_val * lb[c_ind];
		  }
	       }
	       if(rows[r_ind].lb > -DBL_MAX){
		  if(ub[c_ind] >= DBL_MAX){
		     rows[r_ind].lb = -DBL_MAX;
		  }else{
		     rows[r_ind].lb += a_val * ub[c_ind];
		  }
	       }
	    }

	    if(cols[c_ind].var_type == 'F'){
	       rows[r_ind].fixed_obj_offset = obj[c_ind]*ub[c_ind];
	       rows[r_ind].fixed_lhs_offset = a_val * ub[c_ind];
	    }
	 }
      }
   }
  
   // mip->ub[col_ind] = mip->lb[col_ind] = fixed_bound;      

   /*debug */

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
   p_level = params.prep_level;
   etol = params.etol;
   verbosity = params.prep_verbosity;

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
	 if(rows[obj_ind].sr_lb > -DBL_MAX){
	 printf("%f", rows[obj_ind].sr_lb);
	 }else{
	    printf("-inf");
	 }
	 printf("\told_ub:");
	 if(rows[obj_ind].sr_ub < DBL_MAX){
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
		      exit(0);
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
			       obj_ind, old_bound <= -DBL_MAX ? 1 : old_bound, sr->lb);
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
			       obj_ind, old_bound >= DBL_MAX ? -1 : old_bound, sr->ub);
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
	 if(rows[obj_ind].sr_lb > -DBL_MAX){
	    printf("%f", rows[obj_ind].sr_lb);
	 }else{
	    printf("-inf");
      }
	 printf("\tnew_ub:");
	 if(rows[obj_ind].sr_ub < DBL_MAX){
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
	    exit(0); 
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

   double max_dual_ub = DBL_MAX, min_dual_ub = DBL_MAX;
   double max_dual_lb = -DBL_MAX, min_dual_lb = -DBL_MAX;
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
	       if(ub[r_matind[k]] >= DBL_MAX){
		  no_upper = TRUE;
	       }else{
		  *ub_offset += ub[r_matind[k]] * r_matval[k];
	       }
	    }
	    if(!no_lower){
	       if(lb[r_matind[k]] <= -DBL_MAX){
		  no_lower = TRUE;
	       }else{
		  *lb_offset += lb[r_matind[k]] * r_matval[k];
	       }
	    }
	 }else if (r_matval[k] < 0.0){
	    if(!no_lower){
	       if(ub[r_matind[k]] >= DBL_MAX){
		  no_lower = TRUE;
	       }else{
		  *lb_offset += ub[r_matind[k]] * r_matval[k];
	       }
	    }
	    if(!no_upper){
	       if(lb[r_matind[k]] <= -DBL_MAX){
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
	 
/* 	 if(ub[col_ind] < DBL_MAX && lb[col_ind] > -DBL_MAX){ */
/* 	    /\* debug - get vars.type here *\/ */
/* 	    /\* we check this very same thing down too, fix this! *\/ */
/* 	    if(ub[col_ind] > lb[col_ind] + etol){ */
/* 	       /\* debug *\/ */
/* 	       printf("bounded column-first case -case all open row-" */
/* 		      "sr_solve_open_prob(), exiting...\n"); */
/* 	       exit(0);  */
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
	 if(ub[col_ind] < DBL_MAX && lb[col_ind] > -DBL_MAX){
	    /* debug - get vars.type here */
	    if(ub[col_ind] > lb[col_ind] + etol){
		  /* debug */
	       printf("bounded column -case all open row-"
		      "sr_solve_open_prob(), exiting...\n");
	       exit(0); 
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
	 }else if(ub[col_ind] >= DBL_MAX){
	    if(lb[col_ind] > -DBL_MAX){
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
		  //exit(0); 
	       }
	    }else{
	       /* debug */
	       printf("not nonzero???" 
		      "numerical issues -case all open row-"
		      "prep_solve_sr_rlx(), exiting...\n");
	       exit(0); 
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
      //	 exit(0); 
      //  }
      
      if(!no_lower){
	 if(rhs >= 0){
	    if(min_dual_ub < DBL_MAX){
	       sr->lb = min_dual_ub * rhs;
	    }else{
	       prob_infeasible = TRUE;
	    }
	 }else{
	    if(min_dual_lb > -DBL_MAX){
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
	       if(max_dual_ub < DBL_MAX){
		  sr->ub = -(max_dual_ub * rhs);
	       }else{
		  prob_infeasible = TRUE;
	       }
	    }else{
	       if(max_dual_lb > -DBL_MAX){
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

int prep_check_redundancy(MIPdesc *mip, int row_ind,  double etol, 
			  char use_sr_bounds)  
{
   
   int i, termcode = PREP_UNMODIFIED, fix_type = FIX_NO_BOUND;
   int fixed_row = FALSE, col_ind;
   int debug_cnt = 0; 
   double a_val, ub, lb, new_bound, rnd_floor, rnd_ceil;
   
   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols; 

   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double *c_ub = mip->ub;
   double *c_lb = mip->lb;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind; 
   double *r_matval = mip->row_matval;

  /* 
      ratio_type = 
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */
   //   if(row.fixed_var_num + row.fixable_var_num >= row.size){
   if(rows[row_ind].fixed_var_num  >= rows[row_ind].size){
      rows[row_ind].is_redundant = TRUE;
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
	       if(new_bound > c_ub[col_ind] || new_bound < c_lb[col_ind]){
		  return PREP_INFEAS;
	       }
	       if(cols[col_ind].var_type != 'C'){
		  rnd_floor = floor(new_bound);
		  rnd_ceil = ceil(new_bound);
		  if(new_bound >= rnd_floor + etol &&
		     new_bound <= rnd_ceil - etol){
		     return PREP_INFEAS;
		  }else{
		      if(new_bound < rnd_floor + etol){
			 new_bound = rnd_floor;
		      }else{
			 new_bound = rnd_ceil;
		      }
		  }

		  fix_type = FIX_OTHER;

	       }
	    }else if((sense == 'G' && a_val > 0.0) ||
		     (sense == 'L' && a_val < 0.0)){
	       if(new_bound > c_ub[col_ind]){
		  return PREP_INFEAS;
	       }
	       if(new_bound > c_lb[col_ind]){
		  if(cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, 
						   RND_CEIL);
		     fix_type = FIX_LB;
		  }
	       }
	    }else{
	       if(new_bound < c_lb[col_ind]){
		  return PREP_INFEAS;
	       }
	       if(new_bound < c_ub[col_ind]){
		  if(cols[col_ind].var_type != 'C'){
		     new_bound= prep_rnd_integral(new_bound, etol, 
						  RND_FLOOR);
		     fix_type = FIX_UB;
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
	    
	    rows[row_ind].is_redundant = TRUE;
	    if(fix_type != FIX_NO_BOUND){
	       if(cols[col_ind].var_type == 'B'){
		  fix_type = FIX_BINARY;
	       }
	       prep_fixed_col_update_info(mip, col_ind, new_bound, fix_type);
	    }
	    
	    termcode =  PREP_MODIFIED;
	 }
      }
      /* debug */
      if(debug_cnt != 1){
	 printf("error in prep_check_redundancy()\n");
	 exit(0);
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
      
      /*check redundancy and infeasibility*/
      if(lb >= ub + etol){
	 /* debug, can this happen? */
	 return PREP_INFEAS;
      }else if (lb > ub - etol){
	 fixed_row = TRUE;
      }
      
      switch(sense){
       case 'L': 
	  if(lb > rhs + etol){
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
      prep_deleted_row_update_info(mip, row_ind);
   }
  
   return termcode;
}
 
/*===========================================================================*/
/*===========================================================================*/
void prep_deleted_row_update_info(MIPdesc *mip, int row_ind)
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
	    exit(0);
	 }
      }
   }
}
/*===========================================================================*/
/*===========================================================================*/
int prep_initialize_mipinfo(MIPdesc *mip,  prep_params prep_par) {
   int i, j;
   double coef_val, fixed_obj_offset;  
   int row_ind, cont_var_cnt = 0, bin_var_cnt = 0, fixed_var_cnt = 0;   
   int  row_unbounded_cnt, max_row_size, max_col_size;
   int * row_coef_bin_cnt, *row_coef_frac_cnt, *row_sign_pos_cnt;   
   char is_binary, is_bounded, unbounded_below, unbounded_above;
   int gen_type; /* 0 fractional, 1 integer, 2 binary */
   int col_size, col_coef_bin_cnt, col_coef_frac_cnt, col_sign_pos_cnt; 
   int *rows_integerized_var_ind, integerizable_var_num;
   char is_opt_val_integral = TRUE;
   int is_col_all_neg; /* if we convert all constraints to 'L' */
   int is_col_all_pos; /* if we convert all constraints to 'L' */
   double etol = prep_par.etol; 
   int verbosity = prep_par.prep_verbosity;
   //   int p_level = prep_par.prep_level;
   //   int termcode; 

   /* fixme! objsense min max issue!!! will always assume that it is a 
      min problem here!!!! 
   */

   if(!mip){
      if(verbosity >= 2){
	 printf("prep_initialize_mipinfocollect_mipinfo(): Empty mip description...\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
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
      row_coef_frac_cnt = (int *)calloc(ISIZE,m);      
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
	 return(PREP_INFEAS);
      }else if(lb[i] > (ub[i] - etol)){	 
	 cols[i].var_type = 'F';
	 fixed_obj_offset += ub[i];
	 fixed_var_cnt++;	 
	 if(ub[i] >= DBL_MAX || ub[i] <= -DBL_MAX){
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
	 if(lb[i] <= -DBL_MAX){
	    unbounded_below = TRUE;
	 }
	 if(ub[i] >= DBL_MAX){
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

	    if(cols[i].var_type == 'U' ||
	       cols[i].var_type == 'L'){
	       rows[row_ind].fixable_var_num++;
	    }
	 }else{
	    rows[row_ind].fixed_var_num++;
	 }
	 
	 /* for bound types */
	 if(!is_bounded){
	    if(unbounded_above){
	       rows[row_ind].ub_inf_var_num++;
	       if(unbounded_below){
		  rows[row_ind].lb_inf_var_num++;
		  rows[row_ind].free_var_num++;
	       }
	    }else{
	       rows[row_ind].lb_inf_var_num++;
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
	    row_coef_frac_cnt[row_ind]++;
	    col_coef_frac_cnt++;
	 }

	 /* for sign types  and update bounds */
	 
	 if(coef_val > 0.0){
	    row_sign_pos_cnt[row_ind]++;
	    col_sign_pos_cnt++;
	    if(rows[row_ind].ub < DBL_MAX){
	       if(ub[i] >= DBL_MAX){
		  rows[row_ind].ub = DBL_MAX;
	       }else{
		  rows[row_ind].ub += coef_val * ub[i];
	       }
	    }
	    if(rows[row_ind].lb > -DBL_MAX){
	       if(lb[i] <= -DBL_MAX){
		  rows[row_ind].lb = -DBL_MAX;
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
	    if(rows[row_ind].ub < DBL_MAX){
	       if(lb[i] <= -DBL_MAX){
		  rows[row_ind].ub = DBL_MAX;
	       }else{
		  rows[row_ind].ub += coef_val * lb[i];
	       }
	    }
	    if(rows[row_ind].lb > -DBL_MAX){
	       if(ub[i] >= DBL_MAX){
		  rows[row_ind].lb = -DBL_MAX;
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
	    rows[row_ind].fixed_obj_offset = obj[i]*ub[i];
	    rows[row_ind].fixed_lhs_offset = coef_val * ub[i];
	 }
      }


      /* check if this column is fixable */
      /* if so set type, update the obj_offset etc */

      col_size = cols[i].col_size = matbeg[i+1] - matbeg[i];

      if(is_col_all_neg || col_size <= 0){
	 //if((obj[i] > 0.0 && obj_sense == SYM_MAXIMIZE) ||
	 //(obj[i] < 0.0 && obj_sense == SYM_MINIMIZE)){
	 if(obj[i] < 0.0){
	    if(ub[i] >= DBL_MAX){
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
	    if(lb[i] <= -DBL_MAX){
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

      if(row_coef_frac_cnt[j] > 0){
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
      
  
      if(sense[j] == 'E' && rows[j].cont_var_num == 1){
	 if(!(rhs[j]-floor(rhs[j]) > etol &&
	      ceil(rhs[j])-rhs[j] > etol)){
	    if(rows[j].coef_type != FRACTIONAL_VEC){
	       if(cols[rows_integerized_var_ind[j]].var_type != 'Z'){
		  cols[rows_integerized_var_ind[j]].var_type = 'Z';
		  integerizable_var_num++;
	       }
	    }
	 }
      }
      
      rows[j].sr_ub = rows[j].ub;
      rows[j].sr_lb = rows[j].lb;
   }

  /* work on obj */

   if(!(cont_var_cnt - integerizable_var_num)){
      for(i = 0; i < n; i++){
	 coef_val = obj[i];
	 if ((coef_val-floor(coef_val) > etol &&
	      ceil(coef_val)-coef_val > etol)){  
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
   FREE(row_coef_frac_cnt);
   FREE(row_sign_pos_cnt);
   FREE(rows_integerized_var_ind);

   return(PREP_MODIFIED); 
}

/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
int prep_integerize_bounds(PREPdesc *P, char keeptrack)
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
   int verbosity = P->params.prep_verbosity;

   //   int * bounds_updated = (int *)calloc(ISIZE,n);

   for (i = 0; i < n; i++) {
      if (cols[i].var_type != 'F' && 
	  cols[i].var_type != 'C') {
	 diff_ub = diff_lb = 0.0;
	 if (ub[i] < DBL_MAX) {
	    temp_fl = floor(ub[i]);
	    temp_cl = ceil(ub[i]);
	    if(temp_cl - ub[i] < etol ){
	       ub[i] = temp_cl;
	    }else{
	       diff_ub = ub[i] - temp_fl; 
	       ub[i] = temp_fl;	       
	    }
	 }
	 if(lb[i] > -DBL_MAX){
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
	    if (verbosity>=2) {
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

   int i, j, row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg; 
   double * matval, *r_matval;
   int n = mip->n;
   int m = mip->m;
   int nz = mip->nz;
   ROWinfo * rows = mip->mip_inf->rows;

   if(!rows){
      /* debug */
      printf("error in prep_fill_row_ordered - 1");
      exit(0);
   }

   matval = mip->matval;
   matind = mip->matind;
   matbeg = mip->matbeg;
   
   /* allocate space for different arrays */

   r_matval = (mip->row_matval = (double *)malloc(nz*DSIZE)); 
   r_matind = (mip->row_matind = (int *)malloc(nz*ISIZE)); 
   r_matbeg = (mip->row_matbeg = (int *)malloc((m+1)*ISIZE));

   r_matbeg[0] = 0;

   /* got from coin */
   /* debug, get rid of rowinfo dependency here!... */

   for(i = 0; i < m; i++){
      r_matbeg[i + 1] = r_matbeg[i] + rows[i].size;
   }

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
      r_matbeg[i] -= rows[i].size;
   }
   

   /* debug check */
#if 0
   int r_nz = 0;

   for(i = 0; i < m; i++){
      for(j = 0; j < n; j++){
	 for(k = matbeg[j]; k < matbeg[j+1]; k++){
	    if(matind[k] == i){
	       if(r_matind[r_nz] != j || 
		  r_matval[r_nz] != matval[k]){
		  printf("error in prep_fill_row_ordered - 2");
		  exit(0);
	       }
	       r_nz++;
	       break;
	    }
	 }
      }
      if(r_matbeg[i+1] != r_nz){
	 printf("error in prep_fill_row_ordered - 3");
	 exit(0);
      }
   }
   //#endif
#endif 

   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/

#if 0

			 if(r_matval[l] > etol){
			    if(sub_rhs_ub_offset < DBL_MAX){
			       if(lb[r_matval[l]] <= -DBL_MAX){
				  sub_rhs_ub_offset = DBL_MAX;
			       }else{
				  sub_rhs_ub_offset += -(lb[r_matval[l]] * 
							 r_matval[l]);
			       }
			    }
			    if(sub_rhs_lb_offset > -DBL_MAX){
			       if(ub[r_matval[l]] >= DBL_MAX){
				  sub_rhs_lb_offset = -DBL_MAX;
			       }else{
				  sub_rhs_lb_offset += -(ub[r_matval[l]] * 
							 r_matval[l]);
			       }
			    }
			 }else if(r_matval[l] < -etol){
			    if(sub_rhs_lb_offset > -DBL_MAX){
			       if(lb[r_matval[l]] <= -DBL_MAX){
				  sub_rhs_lb_offset = DBL_MAX;
			       }else{
				  sub_rhs_lb_offset += -(lb[r_matval[l]] * 
							 r_matval[l]);
			       }
			    }
			    if(sub_rhs_ub_offset < DBL_MAX){
			       if(ub[r_matval[l]] >= DBL_MAX){
				  sub_rhs_ub_offset = DBL_MAX;
			       }else{
				  sub_rhs_ub_offset += -(ub[r_matval[l]] * 
							 r_matval[l]);
			       }
			    }
			 }		     
			 //		     sub_obj[sub_n] = 0.0;
			 //  sub_matval[sub_n] = r_matval[l];
			 //   sub_matind[sub_n] = r_matind[l];
			 //   sub_n++;
			 l++;
		      }else{
		     
		     if(r_matval[k] < -etol && r_matval[l] > etol){
			

		     }else if(r_matval[k] > etol && r_matval[l] < -etol){



		     }else if(r_matval[k] < -etol && r_matval[l] < -etol){

		     }
		     
		     

		     sub_obj[sub_n] = r_matval[k];
		     sub_matval[sub_n] = r_matval[l];
		     sub_matind[sub_n] = r_matind[l];
		     sub_n++;
		     l++;
		     k++;
		  }
		  
		  if((no_upper && no_lower) ||
		     (sub_rhs_lb_offset <= -DBL_MAX 
		      && sub_rhs_ub_offset >= DBL_MAX)){
		     can_iterate = FALSE;
		     break;
		  }
		  if(k == r_matbeg[obj_ind + 1] && l == r_matbeg[row_ind + 1]){
		     break;
		  }
	       }
	       
	       if(!can_iterate){
		  j--;
		  continue;
	       }else{
		  /* now process this problem */
		  termcode = process_sr_sub_prob(sub_n, sub_obj, sub_matval, 
						 sub_matind, sub_ratios, 
						 sub_rhs, 
						 sub_sense, ub_offset, 
						 lb_offset, sub_rhs_lb_offset, 
						 sub_rhs_ub_offset, ub, lb, 
						 &sub_lb, &sub_ub);
	       }
	    }else{
	       /* couldnt find any other rows to process */
	       break;
	       
	    }
		     
	 }else{
	    /* then do aggregation */
	 }


      }
      
      //    srows[i][j].exists = TRUE;
      //      srows[i][j].row_ind = i; 
      
   }
   return termcode; 
}

int  process_sr_new_col(sr, r_matval[k], r_matval[l], 
			ub[r_matind[k]], lb[r_matind[k]], ub[r_matind[l]],
			lb[r_matind[l]]);

int  process_sr_new_col(SRdec *sr, double c_val, double a_val, int c_sign, 
			int a_sign, double c_ub, 
			double c_lb, double a_ub, double a_lb, double etol){
   
   int c_sign = find_sign(c_val, etol);
   int a_sign = find_sign(a_val, etol);

   if(c_sign == ZERO_VAL && a_sign == ZERO_VAL){
      /* dont add, do nothing */
      return 0;
   }


   if((c_sign == POS_VAL || c_sign == NEG_VAL) && 
      (a_sign == ZERO_VAL)){
  
   }else if((c_sign == ZERO_VAL) && 
	    (a_sign == NEG_VAL || a_sign == POS_VAL)){
      
   }else if(c_sign == POS_VAL && a_sign == NEG_VAL){

   }else if(c_sign == NEG_VAL && a_sign == POS_VAL){
      

   }else if(c_sign == NEG_VAL && a_sign == NEG_VAL){
      

   }else{ /* both pos */
      
   }




}

void prep_compare_val_ub(double * val, double new_ub){

   if(*val > new_ub){
      *val = new_ub;
   }
}
			    




   





	       for( k = r_matbeg[o_ind], l = r_matbeg[c_ind];;){
		  		  
		  if(k < r_matbeg[o_ind + 1] && 
		      (r_matind[k] < r_matind[l] ||
		       l >= r_matbeg[c_ind + 1])){
		     if(r_matval[k] > etol){
			if(!no_upper){
			   if(ub[r_matind[k]] >= DBL_MAX){
			      no_upper = TRUE;
			   }else{
			      ub_offset += ub[r_matind[k]] * r_matval[k];
			   }
			}
			if(!no_lower){
			   if(lb[r_matind[k]] <= -DBL_MAX){
			      no_lower = TRUE;
			   }else{
			      lb_offset += lb[r_matind[k]] * r_matval[k];
			   }
			}
		     }else if(r_matval[k] < -etol){
			if(!no_lower){
			   if(ub[r_matind[k]] >= DBL_MAX){
			      no_lower = TRUE;
			   }else{
			      lb_offset += ub[r_matind[k]] * r_matval[k];
			   }
			}
			if(!no_upper){
			   if(lb[r_matind[k]] <= -DBL_MAX){
			      no_upper = TRUE;
			   }else{
			      ub_offset += lb[r_matind[k]] * r_matval[k];
			   }
			}
		     }
		     k++;
		  }else if(l < r_matbeg[c_ind + 1] && 
			    (r_matind[k] > r_matind[l] ||
			     k >= r_matbeg[o_ind+1])){ 
		     if(r_matval[l] > etol){
			if(sub_rhs_ub_offset < DBL_MAX){
			   if(lb[r_matval[l]] <= -DBL_MAX){
			      sub_rhs_ub_offset = DBL_MAX;
			   }else{
			      sub_rhs_ub_offset += -(lb[r_matval[l]] * 
						     r_matval[l]);
			   }
			}
			if(sub_rhs_lb_offset > -DBL_MAX){
			   if(ub[r_matval[l]] >= DBL_MAX){
			      sub_rhs_lb_offset = -DBL_MAX;
			   }else{
			      sub_rhs_lb_offset += -(ub[r_matval[l]] * 
						     r_matval[l]);
			   }
			}
		     }else if(r_matval[l] < -etol){
			if(sub_rhs_lb_offset > -DBL_MAX){
			   if(lb[r_matval[l]] <= -DBL_MAX){
			      sub_rhs_lb_offset = DBL_MAX;
			   }else{
			      sub_rhs_lb_offset += -(lb[r_matval[l]] * 
						     r_matval[l]);
			   }
			}
			if(sub_rhs_ub_offset < DBL_MAX){
			   if(ub[r_matval[l]] >= DBL_MAX){
			      sub_rhs_ub_offset = DBL_MAX;
			   }else{
			      sub_rhs_ub_offset += -(ub[r_matval[l]] * 
						     r_matval[l]);
			   }
			}
		     }		     
		     //		     sub_obj[sub_n] = 0.0;
		     //  sub_matval[sub_n] = r_matval[l];
		     //   sub_matind[sub_n] = r_matind[l];
		     //   sub_n++;
		     l++;
		  }else{
		     
		     if(r_matval[k] < -etol && r_matval[l] > etol){
			

		     }else if(r_matval[k] > etol && r_matval[l] < -etol){



		     }else if(r_matval[k] < -etol && r_matval[l] < -etol){

		     }
		     
		     

		     sub_obj[sub_n] = r_matval[k];
		     sub_matval[sub_n] = r_matval[l];
		     sub_matind[sub_n] = r_matind[l];
		     sub_n++;
		     l++;
		     k++;
		  }





}


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

