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
#include "sym_master_params.h"
#include "sym_master.h" 
#include "sym_constants.h" 

/*
#include "preprocessor.h"
#include "preprocessor_constants.h"
#include "symphony_api.h"
#include "proccomm.h" 
#include "master.h" 
#include "BB_constants.h" 
*/

/*===========================================================================*/
/*===========================================================================*/
int prep_integerize_bounds(MIPdesc *P, int & bounds_integerized, int verbosity)
{
   /* Change the bounds of integer variables to floor/ceiling appropriately */
   int var_num = 0;
   double ceil_bound = 0;
   double floor_bound = 0;
   for (var_num=0;var_num<P->n;var_num++) {
      if (P->is_int[var_num]==TRUE) {
	 if (P->ub[var_num]<DBL_MAX) {
	    P->ub[var_num] = floor(P->ub[var_num]);
	    if (verbosity>=2) {
	       printf("new upperbound of variable %s: %f\n",
		      P->colname[var_num],P->ub[var_num]);
	    }
	    bounds_integerized++;
	 }
	 if (P->lb[var_num]> -1*DBL_MAX) {
	    P->lb[var_num] = ceil(P->lb[var_num]);
	    if (verbosity>=2) {
	       printf("new lowerbound of variable %s: %f\n",
		      P->colname[var_num],P->lb[var_num]);
	    }
	    bounds_integerized++;
	 }
      }
   }
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_create_row_rep(MIPdesc *P, rowpackedarray *row_P) 
{ 
   /*
     recreates 'A' matrix using three matrices just like the standard
     notation. However, matrices contain row representations rather than
     column
   */ 
   int row_num = 0;		/* the number of the constraint which is being
				   added to the matrices */ 
   int row_start_index = 0; 	/* start index, of the current row, in the
				   array row_matval and row_matind*/ 
   int num_var_in_row = 0; 	/* number of variables with nonzero coeffs in
				   the current row.*/ 
   int row_matind_index = 0;	/* index of the array row_matind */
   int i = 0;			/* counter for loops */
   
   /* allocate space for different arrays */
   row_P->matval = (double *)malloc(P->nz*DSIZE); 
   row_P->matind = (int *)malloc(P->nz*ISIZE); 
   row_P->matbeg = (int *)malloc((P->m+1)*ISIZE);
   row_P->isdeleted = (char *)calloc(P->m+1, CSIZE);
   row_P->matbeg[0] = 0;
   row_P->rows_deleted = 0;

   row_P->isdeleted[P->m]='\0';
   while (row_num<P->m) { 
      /*
	in row number 'row_num', search for variables by checking entries of
	each column
      */  
      int col_num = 0; 
      int j = 0; 
      num_var_in_row = 0; 
      while (col_num<P->n) {	/* columns 0 to (n-1) */
 	 int j = 0;
	 int col_start_index, col_stop_index;
 	 col_start_index = P->matbeg[col_num]; 
 	 col_stop_index = P->matbeg[col_num+1]; 
 	 for (j=col_start_index;j<col_stop_index;j++) { 
 	    if (P->matind[j]==row_num) { 
 	       num_var_in_row++; 
 	       row_P->matind[row_matind_index] = col_num; 
 	       row_P->matval[row_matind_index] = P->matval[j]; 
 	       row_matind_index++; 
 	       break; 
 	    } 
 	 } 
 	 col_num++; 
      } 
      /* column n */ 
      row_num++; 
      row_P->matbeg[row_num]=row_P->matbeg[row_num-1]+num_var_in_row; 
   }
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_find_lhs_params(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs)
{
   /*
     allocate fresh space for lhs structures
   */ 
   int row_num;			/* counter */
   int i;			/* counter */
   /* initialize lhs */
   lhs->ub_is_infinite = (int *)calloc(P->m, ISIZE);
   lhs->lb_is_infinite = (int *)calloc(P->m, ISIZE);
   lhs->ub_inf_var =     (int *)malloc(P->m*ISIZE);
   lhs->lb_inf_var =     (int *)malloc(P->m*ISIZE);
   lhs->ubound =         (double *)calloc(P->m, DSIZE);
   lhs->lbound =         (double *)calloc(P->m, DSIZE);

   /* for each row, find bounds */
   for (row_num=0;row_num<P->m;row_num++) {
      int row_matindex = 0;
      for (row_matindex=row_P->matbeg[row_num];
	   row_matindex<row_P->matbeg[row_num+1]; row_matindex++) {
	 
	 if (row_P->matval[row_matindex]>0) {
	    /* work with upperbound */
	    if (P->ub[row_P->matind[row_matindex]]==DBL_MAX) {
	       /* upperbound is infinite */
	       if(lhs->ub_is_infinite[row_num]==1){
		  /* we have more than one contributors to infinite bounds */
		  lhs->ub_inf_var[row_num] = -1;
	       }
	       else {
		  lhs->ub_inf_var[row_num] = row_P->matind[row_matindex];
	       }
	       lhs->ub_is_infinite[row_num] = 1;
	    }
	    else {
	       /* update upperbound */
	       lhs->ubound[row_num] = lhs->ubound[row_num] +
		  row_P->matval[row_matindex]*
		  P->ub[row_P->matind[row_matindex]];
	    }

	    /* work with lowerbound */
	    if (P->lb[row_P->matind[row_matindex]]==-DBL_MAX) {
	       /* lowerbound is infinite */
	       if(lhs->lb_is_infinite[row_num]==1) {
		  /* we have more than one contributors to infinite bounds */
		  lhs->lb_inf_var[row_num] = -1;
	       }
	       else {
		  lhs->lb_inf_var[row_num] = row_P->matind[row_matindex];
	       }
	       lhs->lb_is_infinite[row_num] = 1;
	    }
	    else {
	       /* update lowerbound */
	       lhs->lbound[row_num] = lhs->lbound[row_num] +
		  row_P->matval[row_matindex]*
		  P->lb[row_P->matind[row_matindex]];
	    }
	 }
	 
	 else if (row_P->matval[row_matindex]<0) {
	    /* case when coefficient is negative */
	    /* first work with upperbound        */
	    if (P->lb[row_P->matind[row_matindex]]==-DBL_MAX) {
	       /* upperbound is infinite */
	       if(lhs->ub_is_infinite[row_num]==1) {
		  /* we have more than one contributors to infinite bounds */
		  lhs->ub_inf_var[row_num] = -1;
	       }
	       else {
		  lhs->ub_inf_var[row_num] = row_P->matind[row_matindex];
	       }
	       lhs->ub_is_infinite[row_num] = 1;
	    }
	    else {
	       /* update upperbound */
	       lhs->ubound[row_num] = lhs->ubound[row_num] +
		  row_P->matval[row_matindex]*
		  P->lb[row_P->matind[row_matindex]];
	    }
	    
	    /* work with lowerbound */
	    if (P->ub[row_P->matind[row_matindex]]==DBL_MAX) {
	       /* lowerbound is infinite */
	       if(lhs->lb_is_infinite[row_num]==1) {
		  /* we have more than one contributors to infinite bounds */
		  lhs->lb_inf_var[row_num] = -1;
	       }
	       else{
		  lhs->lb_inf_var[row_num] = row_P->matind[row_matindex];
	       }
	       lhs->lb_is_infinite[row_num] = 1;
	    }
	    
	    else {
	       /* update lowerbound */
	       lhs->lbound[row_num] = lhs->lbound[row_num] +
		  row_P->matval[row_matindex]*
		  P->ub[row_P->matind[row_matindex]];
	    }
	 }
      }	/* end for row_matindex */
   } /* end for row_num=0... */

   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_tighten_bounds(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			const int max_iterations, int & bounds_tightened, 
			int & vars_fixed, int verbosity)
{
   /* tightens the lower and upper bound of each variable by substituting
      lb/ub of other variables in each constraint and repeating the process
      over and over again untill no further improvements are seen In the
      process, also checks if the problem is infeasible
   */      
   int has_tightened = 1;	/* is 1 if bounds tighten, else 0. steps are
				   repeated until this is 0 */  
   int num_iterations = 0;      /* number of rounds of tightening. in one round
				   bounds of each variable may be updated */  
   int termcode = PREP_UNMODIFIED;
   
   while (has_tightened==1 && num_iterations!=max_iterations) {
      int row_num = 0;	        /* row number whose variable are to be
				   tightened */  
      has_tightened = 0; 

      /* for each row carry out the following */ 
      for (row_num=0; row_num<P->m; row_num++) {
	 if (row_P->isdeleted[row_num]=='1') {
	    continue;
	 }
	 /* type of constraint */
	 if (P->sense[row_num]=='L'||P->sense[row_num]=='E'||
	     P->sense[row_num]=='R') {
	    /* check if tightening of bounds using this row is possible or
	       not */ 
	    if (lhs->lb_is_infinite[row_num] == 1) {
	       if (lhs->lb_inf_var[row_num] != -1) {
		  /* only one variable(lhs->lb_inf_var) bounds may be tightened
		     using this row. find index of the coefficient, of that
		     variable, in rowpacked matrix */ 
		  int i = 0;	/* counter */
		  int position_coeff; /* index of the required coefficient in
					 row_num */ 
		  for (i=row_P->matbeg[row_num];
		       i<row_P->matbeg[row_num+1]; i++) {
		     if (row_P->matind[i]==lhs->lb_inf_var[row_num]) {
			break;
		     }
		  }

		  /*
		    i is now the index of the required coefficient in row_num
		  */ 
		  position_coeff = i;
		  if (row_P->matval[position_coeff]>0) {
		     /* calculate the upper bound */
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (P->rhs[row_num]-lhs->lbound[row_num])/
			row_P->matval[position_coeff];
		     if (new_bound<P->ub[lhs->lb_inf_var[row_num]] - ERR_TOL) {
			if (new_bound<P->lb[lhs->lb_inf_var[row_num]]){
			   if (verbosity>=2) {
			      printf("infeasibility detected... ");
			      printf("probably because of constraint %d\n",
				     row_num);
			   }
			   return PREP_INFEAS;
			}
			
			/* update upperbound */
			if (P->is_int[lhs->lb_inf_var[row_num]]==TRUE) {
			   /* if its an integer, take the floor of the bound */
			   new_bound = floor(new_bound);
			}
			prep_update_const_bounds(P, row_P, lhs, 'U',
						 lhs->lb_inf_var[row_num],
						 new_bound);
			if (verbosity>=2) {
			   printf("new upperbound of variable %s: ",P->colname[lhs->lb_inf_var[row_num]]);
			   printf(" %f, old bound was: ", new_bound);
			   if (P->ub[lhs->lb_inf_var[row_num]]==DBL_MAX) {
			      printf ("INFINITY");
			   }
			   else {
			      printf("%f", P->ub[lhs->lb_inf_var[row_num]]);
			   }
			   printf ("\n\n");
			}
			P->ub[lhs->lb_inf_var[row_num]] = new_bound;
			has_tightened = 1;
			if (P->lb[lhs->lb_inf_var[row_num]] == P->ub[lhs->lb_inf_var[row_num]]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
		  
		  else if (row_P->matval[position_coeff]<0) {
		     /* calculate the lower bound */
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (P->rhs[row_num]-lhs->lbound[row_num])/
			row_P->matval[position_coeff];
		     if (new_bound>P->lb[lhs->lb_inf_var[row_num]] + ERR_TOL) {
			if (new_bound>P->ub[lhs->lb_inf_var[row_num]]){
			   if (verbosity>=2) {
			      printf("infeasibility detected...");
			      printf("probably because of constraint number ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			/* update lowerbound */
			if (P->is_int[lhs->lb_inf_var[row_num]]==TRUE) {
			   /* if its integer, take the ceiling of the bound */
			   new_bound = ceil(new_bound);
			}
			prep_update_const_bounds(P, row_P, lhs, 'L',
						 lhs->lb_inf_var[row_num],
						 new_bound);
			if (verbosity>=2) {
			   printf("new lowerbound of variable %s:",
				  P->colname[lhs->lb_inf_var[row_num]]);
			   printf(" %f, old bound was: ", new_bound);
			   if (P->lb[lhs->lb_inf_var[row_num]]==-DBL_MAX) {
			      printf ("-INFINITY");
			   }
			   else {
			      printf("%f", P->lb[lhs->lb_inf_var[row_num]]);
			   }
			   printf ("\n\n");
			}
			P->lb[lhs->lb_inf_var[row_num]] = new_bound;
			if (P->lb[lhs->lb_inf_var[row_num]] == P->ub[lhs->lb_inf_var[row_num]]) {
			   vars_fixed++;
			   printf("fixed variables = %d\n",vars_fixed);
			}
			else {
			   bounds_tightened++;
			}
			has_tightened = 1;
			termcode = PREP_MODIFIED;
		     }
		  }
	       }
	       /* else no updates for this row */
	    }
	    
	    else {
	       /* if lhs->lb_is_infinite == 0 */
	       int row_index = 0; /* index for row_P->matind */
	       int var_num = 0;	/* variable whose bounds we wish to update */
	       
	       /* new bounds can be checked for each variable present in the
		  row */ 
	       for (row_index=row_P->matbeg[row_num];
		    row_index<row_P->matbeg[row_num+1]; row_index++) {
		  var_num = row_P->matind[row_index];
		  if (row_P->matval[row_index]>0) {
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (P->rhs[row_num]-lhs->lbound[row_num])/
			row_P->matval[row_index]+P->lb[var_num];
		     if (new_bound<P->ub[var_num] - ERR_TOL) {
			if (new_bound<P->lb[var_num]) {
			   if (verbosity>=2) {
			      printf("infeasibility detected ... ");
			      printf("probably because of constraint number ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			/* update the upperbound */
			if (P->is_int[var_num]==TRUE) {
			   /* if its an integer, take the floor of the bound */
			   new_bound = floor(new_bound);
			   printf("Taking floor\n");
			}
			
			prep_update_const_bounds(P, row_P, lhs, 'U', var_num,
						 new_bound);
			if (verbosity>=2) {
			   printf("new upperbound of variable %s:",
				  P->colname[var_num]);
			   printf(" %f, old bound was: ", new_bound);
			   if (P->ub[var_num]==DBL_MAX) {
			      printf ("INFINITY");
			   }
			   else {
			      printf("%f", P->ub[var_num]);
			   }
			   printf ("\n\n");
			}
			P->ub[var_num] = new_bound;
			has_tightened = 1;
			if (P->lb[var_num] == P->ub[var_num]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
		  else if (row_P->matval[row_index]<0) {
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (P->rhs[row_num]-lhs->lbound[row_num])/
			row_P->matval[row_index]+P->ub[var_num];
		     if (new_bound>P->lb[var_num] + ERR_TOL) {
			if (new_bound>P->ub[var_num]) {
			   if (verbosity>=2) {
			      printf("infeasibility detected... ");
			      printf("probably because of constraint number ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			
			/* update the lowerbound */
			if (P->is_int[var_num]==TRUE) {
			   /* if its an integer, take the ceil of the bound */
			   new_bound = ceil(new_bound);
			}
			if (verbosity>=2) {
			   printf("new lowerbound of variable %s: ",
				  P->colname[var_num]);
			   printf(" %f, old bound was: ",
				  new_bound);
			   if (P->lb[var_num]==-DBL_MAX) {
			      printf ("-INFINITY");
			   }
			   else {
			      printf("%f", P->lb[var_num]);
			   }
			   printf ("\n\n");
			}
			prep_update_const_bounds(P, row_P, lhs, 'L', var_num,
						 new_bound);
			P->lb[var_num] = new_bound;
			has_tightened = 1;
			if (P->lb[var_num] == P->ub[var_num]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
	       }
	    }
	 }
	 if (P->sense[row_num]=='G'||P->sense[row_num]=='E'||
	     P->sense[row_num]=='R') {
	    double effective_rhs = 0;
	    effective_rhs = (P->sense[row_num]=='R') ?
	       P->rhs[row_num]-P->rngval[row_num] : P->rhs[row_num];
	    /*
	      check if tightening of bounds using this row is possible or not
	    */ 
	    if (lhs->ub_is_infinite[row_num] == 1) {
	       if (lhs->ub_inf_var[row_num] != -1) {
		  /*
		    only one variable(lhs->ub_inf_var) bounds may be tightened
		    using this row.
		    find index of the coefficient, of that variable, in
		    rowpacked matrix
		  */ 
		  int i = 0;	      /* counter */
		  int position_coeff; /* index of the required coefficient in
					 row_num */ 
		  for (i=row_P->matbeg[row_num];i<row_P->matbeg[row_num+1];
		       i++) {
		     if (row_P->matind[i]==lhs->ub_inf_var[row_num]) {
			break;
		     }
		  }
		  
		  /*
		    i is now the index of the required coefficient in row_num
		  */ 
		  position_coeff = i;
		  if (row_P->matval[position_coeff]>0) {
		     /* calculate the lower bound */
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (effective_rhs-lhs->ubound[row_num])/
			row_P->matval[position_coeff];
		     if (new_bound>P->lb[lhs->ub_inf_var[row_num]] + ERR_TOL) {
			if (new_bound>P->ub[lhs->ub_inf_var[row_num]]){
			   if (verbosity>=2) {
			      printf("infeasibility detected... ");
			      printf("probably because of constraint number ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			
			/* update lowerbound */
			if (P->is_int[lhs->ub_inf_var[row_num]]==TRUE) {
			   /*
			     if its an integer, take the ceiling of the bound
			   */ 
			   new_bound = ceil(new_bound);
			}
			prep_update_const_bounds(P, row_P, lhs, 'L',
						 lhs->ub_inf_var[row_num],
						 new_bound);
			if (verbosity>=2) {
			   printf("new lowerbound of variable %s:",
				  P->colname[lhs->ub_inf_var[row_num]]);
			   printf(" %f, old bound was: ",
				  new_bound);
			   if (P->lb[lhs->ub_inf_var[row_num]]==DBL_MAX) {
			      printf ("-INFINITY");
			   }
			   else {
			      printf("%f", P->lb[lhs->ub_inf_var[row_num]]);
			   }
			   printf ("\n\n");
			}
			P->lb[lhs->ub_inf_var[row_num]] = new_bound;
			has_tightened = 1;
			if (P->lb[lhs->ub_inf_var[row_num]] == P->ub[lhs->ub_inf_var[row_num]]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
		  
		  else if (row_P->matval[position_coeff]<0) {
		     /* calculate the upper bound */
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (effective_rhs-lhs->ubound[row_num])/
			row_P->matval[position_coeff];
		     if (new_bound<P->ub[lhs->ub_inf_var[row_num]] - ERR_TOL) {
			/* update upperbound */
			if (new_bound<P->lb[lhs->ub_inf_var[row_num]]){
			   if (verbosity>=2) {
			      printf("infeasibility detected... ");
			      printf("probably because of constraint number ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			
			if (P->is_int[lhs->ub_inf_var[row_num]]==TRUE) {
			   /* if its an integer, take the floor of the bound */
			   new_bound = floor(new_bound);
			}

			prep_update_const_bounds(P, row_P, lhs, 'U',
						 lhs->ub_inf_var[row_num],
						 new_bound);
			if (verbosity>=2) {
			   printf("new upperbound of variable %s:",
				  P->colname[lhs->ub_inf_var[row_num]]);
			   printf(" %f, old bound was: ", new_bound);
			   if (P->ub[lhs->ub_inf_var[row_num]]==DBL_MAX) {
			      printf ("INFINITY");
			   }
			   else {
			      printf("%f", P->ub[lhs->ub_inf_var[row_num]]);
			   }
			   printf ("\n\n");
			}
			P->ub[lhs->ub_inf_var[row_num]] = new_bound;
			has_tightened = 1;
			if (P->lb[lhs->ub_inf_var[row_num]] == P->ub[lhs->ub_inf_var[row_num]]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
	       }
	       /* else no updates for this row */
	    }
	    
	    else {
	       /* if lhs->ub_is_infinite == 0 */
	       int row_index = 0; /* index for row_P->matind */
	       int var_num = 0;	/* variable whose bounds we wish to update */ 

	       /*
		 new bounds can be checked for each variable present in the
		 row
	       */ 
	       for (row_index=row_P->matbeg[row_num];
		    row_index<row_P->matbeg[row_num+1];row_index++) {
		  var_num = row_P->matind[row_index];
		  if (row_P->matval[row_index]>0) {
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (effective_rhs-lhs->ubound[row_num])/
			row_P->matval[row_index]+P->ub[var_num];
		     if (new_bound>P->lb[var_num] + ERR_TOL) {
			if (new_bound>P->ub[var_num]) {
			   if (verbosity>=2) {
			      printf("infeasibility detected... ");
			      printf("probably because of constraint number ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			
			/* update the lowerbound */
			if (P->is_int[var_num]==TRUE) {
			   /* if its an integer, take the ceil of the bound */
			   new_bound = ceil(new_bound);
			}
			prep_update_const_bounds(P, row_P, lhs, 'L',
						 var_num, new_bound);
			if (verbosity>=2) {
			   printf("new lowerbound of variable %s:",
				  P->colname[var_num]);
			   printf(" %f, old bound was: ", new_bound);
			   if (P->lb[var_num]==-DBL_MAX) {
			      printf ("-INFINITY");
			   }
			   else {
			      printf("%f", P->lb[var_num]);
			   }
			   printf ("\n\n");
			}
			P->lb[var_num] = new_bound;
			has_tightened = 1;
			if (P->lb[var_num] == P->ub[var_num]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
		  else if (row_P->matval[row_index]<0) {
		     double new_bound = 0; /* new bound on the variable */
		     new_bound = (effective_rhs-lhs->ubound[row_num])/
			row_P->matval[row_index]+P->lb[var_num];
		     if (new_bound<P->ub[var_num] - ERR_TOL) {
			if (new_bound<P->lb[var_num]) {
			   if (verbosity>=2) {
			      printf("infeasibility detected... ");
			      printf("probably because of constraint no. ");
			      printf("%d\n", row_num);
			   }
			   return PREP_INFEAS;
			}
			/* update the upperbound */
			if (P->is_int[var_num]==TRUE) {
			   /* if its an integer, take the floor of the
			      bound */
			   new_bound = floor(new_bound);
			}
			prep_update_const_bounds(P, row_P, lhs, 'U',
						 var_num, new_bound);
			if (verbosity>=2) {
			   printf("new upperbound of variable %s:",
				  P->colname[var_num]);
			   printf(" %f, old bound was:",
				  new_bound);
			   if (P->ub[var_num]==DBL_MAX) {
			      printf ("INFINITY");
			   }
			   else {
			      printf("%f", P->ub[var_num]);
			   }
			   printf ("\n\n");
			}
			P->ub[var_num] = new_bound;
			has_tightened = 1;
			if (P->lb[var_num] == P->ub[var_num]) {
			   vars_fixed++;
			}
			else {
			   bounds_tightened++;
			}
			termcode = PREP_MODIFIED;
		     }
		  }
	       }
	    }
	 }
      }	/* for row_num=0... */
      num_iterations++;
   } /* while has_tightened==1 */
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_find_imp_sibling(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			  const int col_num, int verbosity)
{
   int col_start, col_stop, col_index, col_num2;
   int row_start, row_stop, row_index, row_num;
   int termstatus;
   int termcode = PREP_UNMODIFIED;
   double rhs_val;

   return PREP_UNMODIFIED;
   
   col_start = P->matbeg[col_num];
   col_stop = P->matbeg[col_num+1];
   for (col_index=col_start; col_index<col_stop; col_index++) {
      row_num = P->matind[col_index];
      if (row_P->isdeleted[row_num]=='1') {
	 continue;
      }
      /* for each variable in row_num, try to tighten bound */
      row_start = row_P->matbeg[row_num];
      row_stop = row_P->matbeg[row_num+1];
      if (P->sense[row_num] == 'L' || P->sense[row_num] == 'E' ||
	  P->sense[row_num] == 'R') {
	 if (lhs->lb_is_infinite[row_num]) {
	    continue;
	 }
	 for (row_index=row_start; row_index<row_stop; row_index++) {
	    col_num2 = row_P->matind[row_index];
	    if ( !prep_isBinary(P,col_num2) || col_num2==col_num) {
	       continue;
	    }
	    if (row_P->matval[row_index] > 0) {
	       /* see effect of setting col2 at 1*/
	       if (lhs->lbound[row_num]+row_P->matval[row_index] >
		   P->rhs[row_num]) {
		  /* col_num2 has to be fixed at 0 */
		  P->ub[col_num2] = 0;
		  termcode = PREP_MODIFIED;
		  if (verbosity >= 2) {
		     printf ("variable %d fixed at %f after fixing variable");
		     printf (" %d\n", col_num2, P->ub[col_num2], col_num);
		  }
	       }
	    }
	    else if (row_P->matval[row_index] < 0) {
	       /* see effect of setting col2 at 0*/
	       if (lhs->lbound[row_num]+row_P->matval[row_index] > P->rhs[row_num]) {
		  /* col_num2 has to be fixed at 1 */
		  P->lb[col_num2] = 1;
		  termcode = PREP_MODIFIED;
		  if (verbosity >= 2) {
		     printf ("variable %d fixed at %f after fixing variable");
		     printf (" %d\n", col_num2, P->ub[col_num2], col_num);
		  }
	       }
	    }
	 }
      }
      if (P->sense[row_num] == 'G' || P->sense[row_num] == 'E' ||
	  P->sense[row_num] == 'R') {
	 if (lhs->ub_is_infinite[row_num]) {
	    continue;
	 }
	 if (P->sense[row_num] == 'R') {
	    rhs_val = P->rhs[row_num]-P->rngval[row_num];
	 }
	 else {
	    rhs_val = P->rhs[row_num];
	 }
	 for (row_index=row_start; row_index<row_stop; row_index++) {
	    col_num2 = row_P->matind[row_index];
	    if ( !prep_isBinary(P,col_num2) || col_num2==col_num) {
	       continue;
	    }
	    if (row_P->matval[row_index] > 0) {
	       /* see effect of setting col2 at 0*/
	       if (lhs->ubound[row_num]-row_P->matval[row_index] < rhs_val) {
		  /* col_num2 has to be fixed at 1 */
		  P->lb[col_num2] = 1;
		  termcode = PREP_MODIFIED;
		  if (verbosity >= 2) {
		     printf ("variable %d fixed at %f after fixing variable");
		     printf ("%d\n", col_num2, P->ub[col_num2], col_num);
		  }
	       }
	    }
	    else if (row_P->matval[row_index] < 0) {
	       /* see effect of setting col2 at 1*/
	       if (lhs->ubound[row_num]-row_P->matval[row_index] < rhs_val) {
		  /* col_num2 has to be fixed at 0 */
		  P->ub[col_num2] = 0;
		  termcode = PREP_MODIFIED;
		  if (verbosity >= 2) {
		     printf ("variable %d fixed at %f after fixing variable");
		     printf (" %d\n", col_num2, P->ub[col_num2], col_num);
		  }
	       }
	    }
	 }
      }
   }
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_check_feas(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs)
{
   int row_num = 0;
   int termcode = PREP_UNMODIFIED;
   int termstatus = 0;

   /* for each row, check if the lhs bounds are in harmony with the rhs val */
   for (row_num=0; row_num<P->m; row_num++) {
      switch (P->sense[row_num]) {
       case 'L':
	 if (!(lhs->lb_is_infinite[row_num]) &&
	     lhs->lbound[row_num]>P->rhs[row_num]) {
	    return PREP_INFEAS;
	 }
	 break;
       case 'G':
	 if (!(lhs->ub_is_infinite[row_num]) &&
	     lhs->ubound[row_num]<P->rhs[row_num]) {
	    return PREP_INFEAS;
	 }
	 break;
       case 'E':
	 if ( (!(lhs->lb_is_infinite[row_num]) &&
	       lhs->lbound[row_num]>P->rhs[row_num]) ||
	      (!(lhs->ub_is_infinite[row_num]) &&
	       lhs->ubound[row_num]<P->rhs[row_num]) ) {
	    return PREP_INFEAS;
	 }
	 break;
       case 'R':
	 if ( (!(lhs->lb_is_infinite[row_num]) &&
	       lhs->lbound[row_num]>P->rhs[row_num]) ||
	      (!(lhs->ub_is_infinite[row_num]) &&
	       lhs->ubound[row_num] < P->rhs[row_num]-P->rngval[row_num]) ) {
	    return PREP_INFEAS;
	 }
	 break;
      }
   }
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_check_redundancy(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			  int verbosity)
     /*
       checks if any of the constraints are redundant by substituting ub/lb of
       variables into each constraint.
     */
{
   int row_num = 0;
   int termcode = 0;
   int termstatus = 0;
   
   /* for each row... */   
   for (row_num=0;row_num<P->m;row_num++) {
      if (row_P->isdeleted[row_num]=='1') {
	 continue;
      }
      if (P->sense[row_num]=='L') {
	 if ( !(lhs->ub_is_infinite[row_num]) && lhs->ubound[row_num] <=
	      P->rhs[row_num]) {
	    if (verbosity>=2) {
	       printf("Constraint number %d is redundant... Deleting row\n",
		      row_num);
	    }
	    termcode = PREP_MODIFIED;
	    termstatus = prep_delete_row(row_num, P, row_P, lhs);
	    if (termstatus == PREP_SOLVED) {
	       return PREP_SOLVED;
	    }
	 }
      }
      else if (P->sense[row_num]=='G') {
	 if ( !(lhs->lb_is_infinite[row_num]) && lhs->lbound[row_num] >=
	      P->rhs[row_num]) {
	    if (verbosity>=2) {
	       printf("Constraint number %d is redundant... Deleting row\n",
		      row_num);
	    }
	    termcode = PREP_MODIFIED;
	    termstatus = prep_delete_row(row_num, P, row_P, lhs);
	    if (termstatus == PREP_SOLVED) {
	       return PREP_SOLVED;
	    }
	 }
      }
      else if (P->sense[row_num]=='R') {
	 if ( !(lhs->ub_is_infinite[row_num]) &&
	      !(lhs->lb_is_infinite[row_num]) &&
	      lhs->ubound[row_num] <= P->rhs[row_num] && lhs->lbound[row_num]
	      >= P->rhs[row_num] - P->rngval[row_num]) {
	    if (verbosity>=2) {
	       printf("Constraint number %d is redundant... Deleting row\n",
		      row_num);
	    }
	    termcode = PREP_MODIFIED;
	    termstatus = prep_delete_row(row_num, P, row_P, lhs);
	    if (termstatus == PREP_SOLVED) {
	       return PREP_SOLVED;
	    }
	 }
      }
   }
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_fix_variables(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		       int & variables_fixed, int verbosity)
     /*
       Fix variables if possible to their lower or upper bounds depending on
       (a) function is MAX/MIN (b) sign of cost of that variable in the obj
       function (c) sign of coefficient of that variable in each constraint
     */

     /*
       assume that Symphony only deals with minimization. note that there is a
       obj_sense char in the mipdesc structure, but seems empty.
     */

{
   int i = 0;
   int termcode = 0;

   for (i=0;i<P->n;i++) {
      /* i is the variable we are trying to fix */
      int fixing_possible = 1;	/*  is 1 iff the variable i can be fixed*/
      int matind_startindex = P->matbeg[i]; /* the information regarding ith
					       var starts at matind_startindex
					       in array matind and matval */
      int matind_stopindex = P->matbeg[i+1]-1; /* ... stops at
						  matind_stopindex*/ 
      int matind_index = 0;	/* counter */

      
      if (P->obj[i]>=0) {
	 /* try to fix the variable at its lower bound */
	 if (P->lb[i]<=P->ub[i] && P->lb[i]>=P->ub[i]) {
	    fixing_possible = 0; /* variable is already fixed */
	    continue;
	 }
	 for (matind_index=matind_startindex; matind_index<=matind_stopindex;
	      matind_index++) {
	    int row_num = P->matind[matind_index]; /* the row where variable i
						      exists  */
	    if (row_P->isdeleted[row_num]=='1') {
	       continue;
	    }
	    if (P->rngval[row_num]>0 || P->rngval[row_num]<0) {
	       fixing_possible = 0;
	       break;
	    }
	    
	    else if (P->sense[row_num]=='L' &&
		     P->matval[matind_index]<0) {
	       fixing_possible = 0;
	       break;
	    }
	    
	    else if (P->sense[row_num]=='G' &&
		     P->matval[matind_index]>0) {
	       fixing_possible = 0;
	       break;
	    }
	    
	    else if (P->sense[row_num]=='E') {
	       fixing_possible = 0;
	       break;
	    }
	 }

	 
	 if (fixing_possible==1) {
	    if (P->lb[i] == -DBL_MAX) {
	       if (verbosity >= 2) {
		  printf("Problem unbounded, probably ");
		  printf("because of variable number: %d\n", i);
	       }
	       return PREP_UNBOUNDED;
	    }
	    else {
	       if (P->lb[i] != P->ub[i]) {
		  prep_update_const_bounds(P, row_P, lhs, 'U', i, P->lb[i]);
		  P->ub[i]=P->lb[i];
		  variables_fixed++;
		  termcode = PREP_MODIFIED;
		  if (verbosity >= 2) {
		     printf("Variable[%d] fixed at: %f, fixed vars = %d\n", i, P->ub[i], variables_fixed);
		  }
	       }
	    }
	 }
      }
      else if (P->obj[i]<=0) {
	 /* try to fix the variable at its upper bound */
	 if (P->lb[i]<=P->ub[i] && P->lb[i]>=P->ub[i]) {
	    fixing_possible = 0; /* variable is already fixed */
	    continue;
	 }
	 for (matind_index=matind_startindex;matind_index<=matind_stopindex;
	      matind_index++) {
	    int row_num = P->matind[matind_index]; /* the row where variable i
						      exists */
	    if (row_P->isdeleted[row_num]=='1') {
	       continue;
	    }
	    if (P->rngval[row_num]<0 || P->rngval[row_num]>0) {
	       fixing_possible = 0;
	       break;
	    }
	    
	    if (P->sense[row_num]=='L' && P->matval[matind_index]>0) {
	       fixing_possible = 0;
	       break;
	    }
	    
	    else if (P->sense[row_num]=='G'
		     && P->matval[matind_index]<0) {
	       fixing_possible = 0;
	       break;
	    }
	    
	    else if (P->sense[row_num]=='E') {
	       fixing_possible = 0;
	       break;
	    }
	 }

	 if (fixing_possible==1) {
	    if (P->ub[i]==DBL_MAX) {
	       if (verbosity >= 2) {
		  printf("Problem unbounded, probably ");
		  printf("because of variable number: %d\n", i);
	       }
	       return PREP_UNBOUNDED;
	    }
	    else {
	       if (P->lb[i] != P->ub[i]) {
		  prep_update_const_bounds(P, row_P, lhs, 'L', i, P->ub[i]);
		  P->lb[i]=P->ub[i];
		  variables_fixed++;
		  termcode = PREP_MODIFIED;
		  if (verbosity>=2) {
		     printf("Variable[%d] fixed at: %f, fixed vars = %d\n", i, P->lb[i], variables_fixed);
		  }
	       }
	    }
	 }
      }
   }
   return termcode;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_delete_row(int row_num, MIPdesc *P, rowpackedarray *row_P,
		    lhsparams *lhs)
{
   row_P->isdeleted[row_num]='1';
   row_P->rows_deleted = row_P->rows_deleted+1;
   if (row_P->rows_deleted==P->m) {
      return PREP_SOLVED;
   }
   else {
      return PREP_MODIFIED;
   }
}


/*===========================================================================*/
/*===========================================================================*/
int prep_purge_del_rows(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			int & rows_purged)
{
   /* deletes all the rows marked deleted from the MIP structures */
   int vals_deleted = 0;	/* number of values deleted from P->matval */
   
   int    *new_row_matbeg;	/* to replace row_matbeg */
   int    *new_row_matind;	/* to replace row_matind */
   double *new_row_matval;	/* to replace row_matval */
   int    *new_matbeg;		/* to replace P->matval */
   int    *new_matind;		/* to replace P->matind */
   double *new_matval;		/* to replace P->matval */
   double *new_rhs;		/* to replace P->rhs */
   char   *new_sense;		/* to replace P->sense */
   double *new_rngval;		/* to replace P->rngval */
   double *new_lhs_lbound;	/* to replace lhs->lbound */
   double *new_lhs_ubound;	/* to replace lhs->ubound */
   int    *new_lhs_lb_infinite;	/* to replace lhs->lb_is_infinite */
   int    *new_lhs_ub_infinite;	/* to replace lhs->ub_is_infinite */
   int    *new_lhs_lb_inf_var;	/* to replace lhs->lb_inf_var */
   int    *new_lhs_ub_inf_var;	/* to replace lhs->ub_inf_var */
   int    *deleted_rows;	/* array of row numbers to be deleted */
   int length_row = 0;		/* number of coefficients in row which is
				   being deleted */
   int col_num = 0;		/* the column number*/ 
   int row_num = 0;		/* the row number*/ 
   int col_index;		/* indices of col_num in P->matind */
   int row_index;		/* indices of row_num in row_P->matind */
   int elems_to_delete = 0;	/* number of elements to delete */
   int row_count = 0;
   int new_m = P->m-row_P->rows_deleted; /* the new no. of rows */
   int new_nz = 0;		/* the new no. of nonzeros */
   int i,j,k;

   rows_purged = rows_purged+row_P->rows_deleted;
   if (row_P->rows_deleted==0) {
      return PREP_UNMODIFIED;
   }

   if (row_P->rows_deleted==P->m) {
      free (P->matind);
      free (P->matval);
      free (P->rhs);
      free (P->rngval);
      free (P->sense);
      P->nz = 0;
      for (i=0;i<=P->n;i++) {
	 P->matbeg[i]=0;
      }
      P->m = 0;

      free (row_P->matbeg);
      row_P->matbeg = (int *)calloc(1, ISIZE);
      free (row_P->matind);
      free (row_P->matval);
      free (row_P->isdeleted);
      row_P->rows_deleted = 0;

      free (lhs->ub_is_infinite);
      free (lhs->lb_is_infinite);
      free (lhs->ub_inf_var);
      free (lhs->lb_inf_var);
      free (lhs->ubound);
      free (lhs->lbound);
      return PREP_SOLVED;
   }

   for (row_num=0; row_num<P->m; row_num++) {
      if (row_P->isdeleted[row_num]=='1') {
	 row_count++;
	 elems_to_delete = elems_to_delete+(row_P->matbeg[row_num+1]-
					    row_P->matbeg[row_num]);
	 if (row_count==row_P->rows_deleted) {
	    break;
	 }
      }
   }

   new_nz = P->nz - elems_to_delete;
   
   /* allocate space for the replacement arrays (row packed arrays)  */
   new_row_matbeg = (int *)malloc((new_m+1)*ISIZE);
   new_row_matind = (int *)malloc(new_nz*ISIZE);
   new_row_matval = (double *)malloc(new_nz*DSIZE);
   
   /* allocate space for the replacement arrays (column packed arrays)  */   
   new_matbeg = (int *)malloc((P->n+1)*ISIZE);
   new_matind = (int *)malloc(new_nz*ISIZE);
   new_matval = (double *)malloc(new_nz*DSIZE);
   new_rhs =    (double *)malloc(new_m*DSIZE);
   new_sense =  (char *)malloc((new_m+1)*CSIZE); /* includes the last '\0'
						    character as well */ 
   new_rngval = (double *)malloc(new_m*DSIZE);

   /* allocate space for the replacement arrays of lhsparams (constraint
      bounds) */ 
   new_lhs_lbound =      (double *)malloc(new_m*DSIZE);
   new_lhs_ubound =      (double *)malloc(new_m*DSIZE);
   new_lhs_lb_infinite = (int *)malloc(new_m*ISIZE);
   new_lhs_ub_infinite = (int *)malloc(new_m*ISIZE);
   new_lhs_lb_inf_var =  (int *)malloc(new_m*ISIZE);
   new_lhs_ub_inf_var =  (int *)malloc(new_m*ISIZE);

   
   /* first update the row packed arrays */
   row_count = 0;
   i = 0;			/* count of number of elements deleted */
   j = 0;			/* count for new rows */
   k = 0;
   for (row_num=0;row_num<P->m;row_num++) {
      if (row_P->isdeleted[row_num]=='1') {
	 i = i+row_P->matbeg[row_num+1]-row_P->matbeg[row_num];
      }
      else {
	 new_row_matbeg[j] = row_P->matbeg[row_num]-i;
	 for (int h=row_P->matbeg[row_num];h<row_P->matbeg[row_num+1]; h++) {
	    new_row_matind[k] = row_P->matind[h];
	    new_row_matval[k] = row_P->matval[h];
	    k++;
	 }
	 j++;
      }
   }
   new_row_matbeg[j]=row_P->matbeg[row_num]-i;


   /* now update column packed matrix */

   /* recreate column major representation from row major representation from
      scratch. sounds easier than changing the existing representation */

   col_num = 0;
   k = 0;
   new_matbeg[0]=0;
   while (col_num<P->n) {
      row_num = 0;
      j = 0;
      i = 0;
      while (row_num<new_m) {
	 j = 0;
	 for (j=new_row_matbeg[row_num]; j<new_row_matbeg[row_num+1]; j++) {
	    if (new_row_matind[j]==col_num) {
	       i++;
	       new_matind[k] = row_num;
	       new_matval[k] = new_row_matval[j];
	       k++;
	    }
	 }
	 row_num++;
      }
      col_num++;
      new_matbeg[col_num]=new_matbeg[col_num-1]+i;
   }


   /* update P */
   free(P->matind);
   free(P->matval);
   P->nz = new_nz;		/* update nz */
   P->matind = (int *)malloc(P->nz*ISIZE);
   P->matval = (double *)malloc(P->nz*DSIZE);
   for (i=0;i<P->nz;i++) {
      P->matind[i] = new_matind[i];
      P->matval[i] = new_matval[i];
   }
   for (i=0;i<=P->n;i++) {
      P->matbeg[i] = new_matbeg[i];
   }
   
   /* update rhs, sense and rangevalues. Also update the lhs arrays */
   i = 0;
   j = 0;
   for (row_num=0; row_num<P->m; row_num++) {
      if (row_P->isdeleted[row_num]=='1') {
	 continue;
      }
      new_rhs[i]=P->rhs[row_num];
      new_sense[i]=P->sense[row_num];
      new_rngval[i]=P->rngval[row_num];

      new_lhs_lbound[i] = lhs->lbound[row_num];
      new_lhs_ubound[i] = lhs->ubound[row_num];
      new_lhs_lb_infinite[i] = lhs->lb_is_infinite[row_num];
      new_lhs_ub_infinite[i] = lhs->ub_is_infinite[row_num];
      new_lhs_lb_inf_var[i] = lhs->lb_inf_var[row_num];
      new_lhs_ub_inf_var[i] = lhs->ub_inf_var[row_num];
      i++;
   }
   
   /* update m */
   P->m = P->m - row_P->rows_deleted;

   free(P->rhs);
   free(P->sense);
   free(P->rngval);
   free(lhs->lbound);
   free(lhs->ubound);
   free(lhs->lb_is_infinite);
   free(lhs->ub_is_infinite);
   free(lhs->lb_inf_var);
   free(lhs->ub_inf_var);
   
   P->rhs =              (double *)malloc(P->m*DSIZE);
   P->rngval =           (double *)malloc(P->m*DSIZE);
   P->sense =            (char *)malloc((P->m+1)*CSIZE);
   lhs->lbound =         (double *)malloc(P->m*DSIZE);
   lhs->ubound =         (double *)malloc(P->m*DSIZE);
   lhs->lb_is_infinite = (int *)malloc(P->m*ISIZE);
   lhs->ub_is_infinite = (int *)malloc(P->m*ISIZE);
   lhs->lb_inf_var =     (int *)malloc(P->m*ISIZE);
   lhs->ub_inf_var =     (int *)malloc(P->m*ISIZE);

   
   for (i=0;i<P->m;i++) {
      P->rhs[i] = new_rhs[i];
      P->sense[i] = new_sense[i];
      P->rngval[i] = new_rngval[i];
      lhs->lbound[i] = new_lhs_lbound[i];
      lhs->ubound[i] = new_lhs_ubound[i];
      lhs->lb_is_infinite[i] = new_lhs_lb_infinite[i];
      lhs->ub_is_infinite[i] = new_lhs_ub_infinite[i];
      lhs->lb_inf_var[i] = new_lhs_lb_inf_var[i];
      lhs->ub_inf_var[i] = new_lhs_ub_inf_var[i];
   }
   P->sense[P->m] = '\0';

   /* update row-packed vectors */
   free(row_P->matind);
   free(row_P->matbeg);
   free(row_P->matval);
   row_P->matbeg = (int *)malloc((P->m+1)*ISIZE);
   row_P->matind = (int *)malloc(P->nz*ISIZE);
   row_P->matval = (double *)malloc(P->nz*DSIZE);

   for (i=0;i<=P->m;i++) {
      row_P->matbeg[i]=new_row_matbeg[i];  
      row_P->isdeleted[i]='0';
   }
   row_P->isdeleted[P->m]='\0';
   
   for (i=0;i<P->nz;i++) {
      row_P->matind[i]=new_row_matind[i];
      row_P->matval[i]=new_row_matval[i];
   }
   row_P->rows_deleted = 0;
   
   free(new_matval);
   free(new_matind);
   free(new_matbeg);
   free(new_rhs);
   free(new_sense);
   free(new_rngval);
   free(new_row_matval);
   free(new_row_matind);
   free(new_row_matbeg);
   free(new_lhs_lbound);
   free(new_lhs_ubound);
   free(new_lhs_lb_infinite);
   free(new_lhs_ub_infinite);
   free(new_lhs_lb_inf_var);
   free(new_lhs_ub_inf_var);
   return PREP_MODIFIED;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_purge_del_rows2(sym_environment *env, MIPdesc *P, rowpackedarray *row_P, int & rows_purged)
{
   int r_count = 0;
   int *r_indices = (int *) malloc(ISIZE*row_P->rows_deleted);
   for (int r=0; r<P->m; r++) {
      if (row_P->isdeleted[r]=='1') {
	 /* add this row to the array: r_indices */
         r_indices[r_count] = r;
	 r_count++;
      }
   }
   sym_delete_rows(env, r_count, r_indices);
   rows_purged = r_count;
   printf("Purged %d rows\n.",r_count);
   free(r_indices);
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_declare_solved(MIPdesc *P, int verbosity)
{
   /*
     declares that the preprocessor solved the problem fully.
   */

   int i = 0;			/* counter */
   if (verbosity>=1) {
      printf("Problem solved by preprocessor alone. Displaying results:\n\n");
      printf("column name\tvalue\n");
      for (i=0;i<P->n;i++) {
	 printf("%s\t\t%f\n",P->colname[i],P->lb[i]);
      }
   }
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_update_const_bounds(MIPdesc *P, rowpackedarray *row_P,
			     lhsparams *lhs, char bnd_sense, int var_num,
			     double new_bound)
{
   /*
     update the constraint bounds as and when bounds of a variable are
     tightened.

     really big function in terms of number of lines, but seems efficient. size
     mainly on account of lot of cases/subcases
   */
   
   int i; 			/* counter */
   int row_num = 0;
   for (i=P->matbeg[var_num];i<P->matbeg[var_num+1];i++) {
      int j = 0;		/* counter */
      int row_index = 0;	/* index of coefficient of the var_num in
				   rowpackedarray */ 
      row_num = P->matind[i];	/* we have to update bounds for this row */

      if (row_P->isdeleted[row_num]=='1') {
	 continue;
      }
      
      /* find index of coefficient of var_num in row_num in the row_P */
      for(j=row_P->matbeg[row_num];j<row_P->matbeg[row_num+1];j++) {
	 if (row_P->matind[j]==var_num) {
	    break;
	 }
      }
      
      /* j is now the required index */
      row_index = j;
      if (row_P->matval[row_index]>0) {
	 switch (bnd_sense) {
	  case 'L':
	    /* update constraint bounds using the new lowerbndlower */
	    if (lhs->lb_is_infinite[row_num]==1) {
	       if (lhs->lb_inf_var[row_num]==-1) {
		  /* in case of infinite bounds, recalculate lb entirely */
		  /* we will store original P->lb in a temporary variable, make
		     P->lb 
		     the new_bound and then restore the old value to P->lb  */
		  /* untested */
		  double tmp = P->lb[var_num];
		  P->lb[var_num] = new_bound;
		  /* end untested */
		  /* could be erroneous */
		  int row_matindex;
		  lhs->lb_is_infinite[row_num] = 0;
		  lhs->lbound[row_num] = 0;
		  for (row_matindex=row_P->matbeg[row_num];
		       row_matindex<row_P->matbeg[row_num+1]; row_matindex++){
		     if (row_P->matval[row_matindex] > 0) {
			if (P->lb[row_P->matind[row_matindex]]==-DBL_MAX) {
			   /* lowerbound is infinite */
			   if(lhs->lb_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to infinite
				bounds 
			      */  
			      lhs->lb_inf_var[row_num] = -1;
			      break;
			   }   
			   else {
			      lhs->lb_inf_var[row_num] =
				 row_P->matind[row_matindex];
			   }
			   
			   lhs->lb_is_infinite[row_num] = 1;
			}
			else {
			   /* update lowerbound */
			   lhs->lbound[row_num] = lhs->lbound[row_num] +
			      row_P->matval[row_matindex]*
			      P->lb[row_P->matind[row_matindex]];
			}
		     }
		     else if (row_P->matval[row_matindex]<0) {
			if (P->ub[row_P->matind[row_matindex]]==DBL_MAX) {
			   /* lowerbound is infinite */
			   if(lhs->lb_is_infinite[row_num]==1) {
			      /* we have more than one contributors to
				 infinite bounds */ 
			      lhs->lb_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->lb_inf_var[row_num] =
				 row_P->matind[row_matindex]; 
			   }
			   lhs->lb_is_infinite[row_num] = 1;
			}
			else {
			   /* update lowerbound */
			   lhs->lbound[row_num] = lhs->lbound[row_num]+
			      row_P->matval[row_matindex]*
			      P->ub[row_P->matind[row_matindex]];
			}
		     }
		  }
		  /* restore the value of P->lb to its original value */
		  /* untested */
		  P->lb[var_num] = tmp;
		  /* end untested */
		  /* could be erroneous */
		  
	       }
	       else if (lhs->lb_inf_var[row_num]==var_num) {
		  /* infinite bound is replaced by a finite one */
		  lhs->lb_is_infinite[row_num] = 0;
		  lhs->lbound[row_num] = lhs->lbound[row_num] +
		     new_bound*row_P->matval[row_index];
	       }
	       else {
		  /*
		    bound is still infinite, but the remaining-sum is updated
		  */ 
		  lhs->lbound[row_num] = lhs->lbound[row_num] +
		     (new_bound-P->lb[var_num])*row_P->matval[row_index];
	       }
	    }
	    else {
	       /* update finite lbound */
	       lhs->lbound[row_num] = lhs->lbound[row_num] +
		  (new_bound-P->lb[var_num])*row_P->matval[row_index];
	    }
	    break;

	  case 'U':
	    /* update upper bound of the constraint */
	    if (lhs->ub_is_infinite[row_num]==1) {
	       if (lhs->ub_inf_var[row_num]==-1) {

		  /* in case of infinite bounds, recalculate ub entirely */
		  /* we will store original P->ub in a temporary variable, make
		     P->ub the new_bound and then restore the old value to
		     P->ub*/
		  
		  double tmp = P->ub[var_num];
		  P->ub[var_num] = new_bound;

		  int row_matindex;
		  lhs->ub_is_infinite[row_num] = 0;
		  lhs->ubound[row_num] = 0;
		  for (row_matindex=row_P->matbeg[row_num];
		       row_matindex<row_P->matbeg[row_num+1];
		       row_matindex++) {
		     if (row_P->matval[row_matindex]>0) {
			if (P->ub[row_P->matind[row_matindex]]==DBL_MAX) {
			   /* upperbound is infinite */
			   if(lhs->ub_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to infinite
				bounds
			      */  
			      lhs->ub_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->ub_inf_var[row_num] =
				 row_P->matind[row_matindex]; 
			   }
			   lhs->ub_is_infinite[row_num] = 1;
			}
			else {
			   /* update upperbound */
			   lhs->ubound[row_num] = lhs->ubound[row_num] +
			      row_P->matval[row_matindex]*
			      P->ub[row_P->matind[row_matindex]];
			}
		     }
		     
		     else if (row_P->matval[row_matindex]<0) {
			if (P->lb[row_P->matind[row_matindex]]==-DBL_MAX) {
			   /* upperbound is infinite */
			   if(lhs->ub_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to infinite
				bounds
			      */
			      lhs->ub_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->ub_inf_var[row_num] =
				 row_P->matind[row_matindex]; 
			   }
			   lhs->ub_is_infinite[row_num] = 1;
			}
			else {
			   /* update upperbound */
			   lhs->ubound[row_num] = lhs->ubound[row_num] +
			      row_P->matval[row_matindex]*
			      P->lb[row_P->matind[row_matindex]];
			}
		     }
		  }
		  /* restore original P->ub to its old value  */
		  /* untested */
		  P->ub[var_num] = tmp;
		  /* end untested */
		  /* could be erroneous */
	       }
	       
	       else if (lhs->ub_inf_var[row_num]==var_num) {
		  /* infinite bound is replaced by a finite one */
		  lhs->ub_is_infinite[row_num] = 0;
		  lhs->ubound[row_num] = lhs->ubound[row_num] +
		     new_bound*row_P->matval[row_index];
	       }
	       else {
		  /*
		    bound is still infinite, but the remaining-sum is updated
		  */
		  lhs->ubound[row_num] = lhs->ubound[row_num] +
		     (new_bound-P->ub[var_num])*row_P->matval[row_index];
	       }
	    }
	    
	    else {
	       /* update finite ubound */
	       lhs->ubound[row_num] = lhs->ubound[row_num] +
		  (new_bound-P->ub[var_num])*row_P->matval[row_index];
	    }
	    break;
	 } /* end switch bnd_sense */
      }	/* end if (coeff>0) */
      
      else {			/* if the coefficient is -tive */
	 switch (bnd_sense) {
	  case 'L':
	    /* update upper bound of the constraint */
	    /* remember that the coeff is negative */
	    if (lhs->ub_is_infinite[row_num]==1) {
	       if (lhs->ub_inf_var[row_num]==-1) {
		  /* in case of infinite bounds, recalculate lb entirely */
		  /* we will store original P->lb in a temporary variable, make
		     P->lb the new_bound and then restore the old value to
		     P->lb
		  */ 
		  double tmp = P->lb[var_num];
		  P->lb[var_num] = new_bound;

		  int row_matindex;
		  lhs->ub_is_infinite[row_num] = 0;
		  lhs->ubound[row_num] = 0;
		  for (row_matindex=row_P->matbeg[row_num];
		       row_matindex<row_P->matbeg[row_num+1];
		       row_matindex++) {
		     if (row_P->matval[row_matindex]>0) {
			if (P->ub[row_P->matind[row_matindex]]==DBL_MAX) {
			   /* upperbound is infinite */
			   if(lhs->ub_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to
				infinite bounds
			      */ 
			      lhs->ub_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->ub_inf_var[row_num] =
				 row_P->matind[row_matindex];
			   }
			   lhs->ub_is_infinite[row_num] = 1;
			}
			
			else {
			   /* update upperbound */
			   lhs->ubound[row_num] = lhs->ubound[row_num] +
			      row_P->matval[row_matindex]*
			      P->ub[row_P->matind[row_matindex]];
			}
		     }
		     
		     else if (row_P->matval[row_matindex]<0) {
			if (P->lb[row_P->matind[row_matindex]]==-DBL_MAX) {
			   /* upperbound is infinite */
			   if(lhs->ub_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to
				infinite bounds
			      */
			      lhs->ub_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->ub_inf_var[row_num] =
				 row_P->matind[row_matindex];
			   }
			   lhs->ub_is_infinite[row_num] = 1;
			}
			else {
			   /* update upperbound */
			   lhs->ubound[row_num] = lhs->ubound[row_num] +
			      row_P->matval[row_matindex]*
			      P->lb[row_P->matind[row_matindex]];
			}
		     }
		  }
		  /*restore the old value to P->lb  */
		  /* untested */
		  P->lb[var_num] = tmp;
		  /* end untested */
		  /* could be erroneous */
	       }
	       
	       else if (lhs->ub_inf_var[row_num]==var_num) {
		  /* infinite bound is replaced by a finite one */
		  lhs->ub_is_infinite[row_num] = 0;
		  lhs->ubound[row_num] = lhs->ubound[row_num] +
		     new_bound*row_P->matval[row_index];
	       }
	       
	       else {
		  /*
		    bound is still infinite, but the remaining-sum is updated
		  */ 
		  lhs->ubound[row_num] = lhs->ubound[row_num] +
		     (new_bound-P->lb[var_num])*row_P->matval[row_index];
	       }
	    }
	    else {
	       /* update finite ubound */
	       lhs->ubound[row_num] = lhs->ubound[row_num] +
		  (new_bound-P->lb[var_num])*row_P->matval[row_index];
	    }
	    break;

	  case 'U':
	    /* update lower bound of the constraint */
	    if (lhs->lb_is_infinite[row_num]==1) {
	       if (lhs->lb_inf_var[row_num]==-1) {
		  /* in case of infinite bounds, recalculate lb entirely */

		  /* we will store original P->lb in a temporary variable, make
		     P->lb the new_bound and then restore the old value to
		     P->lb
		  */ 
		  double tmp = P->ub[var_num];
		  int row_matindex;
		  P->ub[var_num] = new_bound;
		  lhs->lb_is_infinite[row_num] = 0;
		  lhs->lbound[row_num] = 0;
		  for (row_matindex=row_P->matbeg[row_num];
		       row_matindex < row_P->matbeg[row_num+1];
		       row_matindex++) {
		     if (row_P->matval[row_matindex]>0) {
			if (P->lb[row_P->matind[row_matindex]]==-DBL_MAX) {
			   /* lowerbound is infinite */
			   if(lhs->lb_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to infinite
				bounds
			      */
			      lhs->lb_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->lb_inf_var[row_num] =
				 row_P->matind[row_matindex];
			      lhs->lb_is_infinite[row_num] = 1;
			   }
			}
			else {
			   /* update lowerbound */
			   lhs->lbound[row_num] = lhs->lbound[row_num] +
			      row_P->matval[row_matindex]*
			      P->lb[row_P->matind[row_matindex]];
			}
		     }
		     
		     else if (row_P->matval[row_matindex]<0) {
			if (P->ub[row_P->matind[row_matindex]]==DBL_MAX) {
			   /* lowerbound is infinite */
			   if(lhs->lb_is_infinite[row_num]==1) {
			      /*
				we have more than one contributors to infinite
				bounds
			      */
			      lhs->lb_inf_var[row_num] = -1;
			      break;
			   }
			   else {
			      lhs->lb_inf_var[row_num] =
				 row_P->matind[row_matindex];
			   }
			   lhs->lb_is_infinite[row_num] = 1;
			}
			else {
			   /* update lowerbound */
			   lhs->lbound[row_num] = lhs->lbound[row_num] +
			      row_P->matval[row_matindex]*
			      P->ub[row_P->matind[row_matindex]];
			}
		     }
		  }
		  /* restore the value of P->ub  */
		  /* untested */
		  P->ub[var_num] = tmp;
		  /* end untested */
		  /* could be erroneous */
		  
	       }
	       
	       else if (lhs->lb_inf_var[row_num]==var_num) {
		  /* infinite bound is replaced by a finite one */
		  lhs->lb_is_infinite[row_num] = 0;
		  lhs->lbound[row_num] = lhs->lbound[row_num] +
		     new_bound*row_P->matval[row_index];
	       }
	       else {
		  /*
		    bound is still infinite, but the remaining-sum is updated
		  */
		  lhs->lbound[row_num] = lhs->lbound[row_num] +
		     (new_bound-P->ub[var_num])*
		     row_P->matval[row_index];
	       }
	    }
	    else {
	       /* update finite lbound */
	       lhs->lbound[row_num] = lhs->lbound[row_num] +
		  (new_bound-P->ub[var_num])*
		  row_P->matval[row_index];
	    }
	    break;
	 } /* end switch bnd_sense */
      }	/* end if (coeff<0) */
   }
}


/*===========================================================================*/
/*===========================================================================*/
int prep_BinImproveCoeffs(MIPdesc *P, rowpackedarray *row_P, lhsparams
			  *lhs, int & coeffs_changed, int verbosity)
{ 
   /*
     For each binary variable, change the coefficient in the A matrix. if
     possible
   */
	 
   int row_num;		        /* row number */ 
   int col_num;			/* column numbeer */ 
   int last_row_index;		/* last column index row->matbeg */ 
   int row_index;		/* column index found in row->matbeg */ 
   double newCoeff;		/* coefficient of binary variable*/ 
   double coeff;		/* coefficient of binary variable*/ 
   double delta;		/* change in constraint bound on fixing a bin 
 				   variable */
   double zk;			/* upper/lower bound of a constraint without
				   the kth var */
   double newRhs;		/* new rhs value of a constraint */
   int termcode;		/* termcode returned by this function */

   termcode = 0;

   /* for each constraint, do the following  */ 
   for (row_num=0; row_num < P->m; row_num++) {
      if (row_P->isdeleted[row_num]=='1') {
	 continue;
      }
      if (P->sense[row_num]=='L') { 
	 if (!lhs->ub_is_infinite[row_num]) { 
	    row_index = row_P->matbeg[row_num]; 
	    last_row_index = row_P->matbeg[row_num+1]; 
	    /* for each binary variable in this row, do the following */ 
	    for(; row_index < last_row_index; row_index++) { 
	       col_num = row_P->matind[row_index]; 
	       if (!prep_isBinary(P, col_num)) {
		  /* skip this variable */
		  continue; 
	       }
	       /* this var. is binary. try changing its coefficient and the
		  corresponding rhs */
	       coeff = row_P->matval[row_index];
	       zk = lhs->ubound[row_num] - coeff;

	       if ( zk < P->rhs[row_num]) { 
		  if (coeff>0) {
		     /* x[col_num]=0 makes this constraint redundant */ 
		     delta = P->rhs[row_num] - zk; 
		     newCoeff = coeff - delta;
		     newRhs = zk; /* rhs has changed */
		     coeffs_changed++;
		     termcode = PREP_MODIFIED;
		     if (verbosity >= 2) {
			printf("coeff for col %d, row %d", col_num, row_num);
			printf("changed from %g to %g\n", coeff, newCoeff);
			printf("\nRow %d has %d elements",row_num, last_row_index-row_P->matbeg[row_num]+1);
		     }
		     /* Change coefficent, rhs and the related bounds for that
			constraint */
		     prep_changeCoeff(P, row_P, lhs, col_num, row_num,
				      newCoeff, newRhs);
		  } 
		   
		  else if (coeff<0) {
		     /* x[col_num]=1 makes this constraint redundant */
		     delta = P->rhs[row_num] - zk;
		     newCoeff = coeff + delta;
		     newRhs = P->rhs[row_num]; /* rhs has not changed */
		     coeffs_changed++;
		     termcode = PREP_MODIFIED;
		     if (verbosity >= 2) {
			printf("coeff for col %d, row %d", col_num, row_num);
			printf("changed from %g to %g\n", coeff, newCoeff);
			printf("\nRow %d has %d elements", row_num, last_row_index-row_P->matbeg[row_num]+1);
		     }
		     /* Change coefficent, rhs and the related bounds for that
			constraint */
		     prep_changeCoeff(P, row_P, lhs, col_num, row_num,
				      newCoeff, newRhs);
		  }
	       }
	    }
	 }
      }
      else if (P->sense[row_num]=='G') { 
	 if (!lhs->lb_is_infinite[row_num]) { 
	    row_index = row_P->matbeg[row_num]; 
	    last_row_index = row_P->matbeg[row_num+1]; 
	    /* for each binary variable in this row, do the following */ 
	    for(; row_index < last_row_index; row_index++) { 
	       col_num = row_P->matind[row_index]; 
	       if (!prep_isBinary(P, col_num)) {
		  /* skip this variable */
		  continue; 
	       }
	       /* this var. is binary. try changing its coefficient and the
		  corresponding rhs */
	       coeff = row_P->matval[row_index];
	       zk = lhs->lbound[row_num] + coeff;
	       if ( zk > P->rhs[row_num]) {
		  if (coeff>0) {
		     /* x[col_num]=1 makes this constraint redundant */ 
		     delta = zk - P->rhs[row_num]; 
		     newCoeff = coeff - delta;
		     newRhs = P->rhs[row_num]; /* rhs has not changed */
		     coeffs_changed++;
		     termcode = PREP_MODIFIED;
		     if (verbosity >= 2) {
			printf("coeff for col %d, row %d", col_num, row_num);
			printf("changed from %g to %g\n", coeff, newCoeff);
			printf("\nRow %d has %d elements",row_num, last_row_index-row_P->matbeg[row_num]+1);
		     }
		     /* Change coefficent, rhs and the related bounds for that
			constraint */
		     prep_changeCoeff(P, row_P, lhs, col_num, row_num,
				      newCoeff, newRhs);
		  }
		  else if (coeff < 0) {
		     zk = lhs->lbound[row_num];
		     /* x[col_num]=0 makes this constraint redundant */
		     delta = zk - P->rhs[row_num];
		     newCoeff = coeff + delta;
		     newRhs = zk; /* new and old are same */
		     coeffs_changed++;
		     termcode = PREP_MODIFIED;
		     if (verbosity >= 2) {
			printf("coeff for col %d, row %d", col_num, row_num);
			printf("changed from %g to %g\n", coeff, newCoeff);
			printf("\nRow %d has %d elements",row_num, last_row_index-row_P->matbeg[row_num]+1);
		     }
		     /* Change coefficent, rhs and the related bounds for that
			constraint */
		     prep_changeCoeff(P, row_P, lhs, col_num, row_num,
				      newCoeff, newRhs);
		  }
	       }
	    } 
	 }
      } 
      /* Nothing to be done for the case when sense = 'R' */
   } 
   return termcode; 
} 


/*===========================================================================*/
/*===========================================================================*/
int prep_isBinary(MIPdesc *P, int col_num) 
{ 
   if (P->is_int[col_num]==TRUE && P->ub[col_num] == 1 && P->lb[col_num] == 0) { 
      return 1; 
   }
   else { 
      return 0; 
   } 
} 


/*===========================================================================*/
/*===========================================================================*/
int prep_changeCoeff(MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
		     int col_num, int row_num, const double newCoeff,
		     const double newRhs) 
{ 
   /* 
      change coefficient to newCoeff in P and row_P. 
      update the bounds of lhs_params. 
   */ 
   double oldCoeff; 		/* original coefficient in P and row_P */
   double oldRhs; 		/* original value of rhs */
   int col_index = P->matbeg[col_num]; 
   int col_last_index = P->matbeg[col_num+1]; 
   int row_index = row_P->matbeg[row_num]; 
   int row_last_index = row_P->matbeg[row_num+1]; 

   /* first change the old coefficient in P */ 
   for (; col_index < col_last_index; col_index++) { 
      if (P->matind[col_index]==row_num) { 
 	 oldCoeff = P->matval[col_index]; 
	 /* just update the coefficient */
	 P->matval[col_index] = newCoeff; 
	 break; 
      } 
   } 


   /* now change the old coefficient in row_P */
   for (; row_index < row_last_index; row_index++) { 
      if (row_P->matind[row_index]==col_num) { 
	 /* just update the coefficient */
	 row_P->matval[row_index] = newCoeff; 
	 break; 
      } 
   }


   /* now change the rhs */
   oldRhs = P->rhs[row_num];
   P->rhs[row_num] = newRhs;


   /* since lhs has changed, change the lhs_params  */
   prep_coeffChangeBounds(P, row_P, lhs, col_num, row_num, oldCoeff, newCoeff);
   return 0; 
} 


/*===========================================================================*/
/*===========================================================================*/

void prep_coeffChangeBounds (MIPdesc *P, rowpackedarray *row_P, lhsparams *lhs,
			     int col_num, int row_num, double oldCoeff,
			     double newCoeff) 
{ 
   /* 
      here we are changing lhsparams because somehow the coefficients in lhs
      were changed. It is assumed that this function is called only when a
      coefficient of a binary variable is changed. So we are not checking for
      unboundedness. Could be fatal if this function is called for something
      else.
   */ 

   int row_index;			/* row index */ 
   int col_index;			/* column index */

   /* first find the position of the element in P  */
   for (col_index=P->matbeg[col_num]; col_index<P->matbeg[col_num];
	col_index++){ 
      if (P->matind[col_index]==row_num) { 
 	 break; 
      } 
   } 
      
   /* eleminate the contribution of old coefficient  */
   if(oldCoeff>0) { 
      if (!(lhs->lb_is_infinite[row_num])) {
	 lhs->lbound[row_num] - oldCoeff*P->lb[col_num];
      }
      if (!(lhs->ub_is_infinite[row_num])) {
	 lhs->ubound[row_num] - oldCoeff*P->ub[col_num];
      }
   } 
   else if (oldCoeff < 0){ 
      if (!(lhs->lb_is_infinite[row_num])) {
	 lhs->lbound[row_num] - oldCoeff*P->ub[col_num];
      }
      if (!(lhs->ub_is_infinite[row_num])) {
	 lhs->ubound[row_num] - oldCoeff*P->lb[col_num];
      }
   } 


   /* add in the contribution of new coefficient  */
   if(newCoeff>0) { 
      if (!(lhs->lb_is_infinite[row_num])) {
	 lhs->lbound[row_num] + newCoeff*P->lb[col_num];
      }
      if (!(lhs->ub_is_infinite[row_num])) {
	 lhs->ubound[row_num] + newCoeff*P->ub[col_num];
      }
   } 
   else if (newCoeff < 0){ 
      if (!(lhs->lb_is_infinite[row_num])) {
	 lhs->lbound[row_num] + newCoeff*P->ub[col_num];
      }
      if (!(lhs->ub_is_infinite[row_num])) {
	 lhs->ubound[row_num] + newCoeff*P->lb[col_num];
      }
   } 
}


/*===========================================================================*/
/*===========================================================================*/
int prep_delete_structs (MIPdesc *fP, rowpackedarray *row_fP, lhsparams *flhs,
			 prep_stats *stats)
{ 
   /* 
      free memory from the copies of MIP. 
   */ 

   free_mip_desc(fP);
   
   prep_free_row_structs(flhs, row_fP);

   free(fP);
   free(row_fP); 
   free(flhs);

   free(stats);
} 


/*===========================================================================*/
/*===========================================================================*/
int prep_free_row_structs(lhsparams *lhs, rowpackedarray *row_P)
{
   /* remove the row packed matrices */
   free(lhs->ub_is_infinite);
   free(lhs->lb_is_infinite);
   free(lhs->ub_inf_var);
   free(lhs->lb_inf_var);
   free(lhs->ubound);
   free(lhs->lbound);     

   free(row_P->matbeg);
   free(row_P->matind);
   free(row_P->matval);
   free(row_P->isdeleted);
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_disp_stats (prep_stats *stats)
{
   printf ("======== Statistics from preprocessor ======== \n");
   printf ("Rows deleted: %d\n", stats->rows_deleted);
   printf ("Variables fixed: %d\n", stats->vars_fixed);
   printf ("Coefficients changed: %d\n", stats->coeffs_changed);
   printf ("Bounds Tightened: %d\n", stats->bounds_tightened);
   printf ("============================================== \n");
}


/*===========================================================================*/
/*===========================================================================*/
int prep_display_mip(MIPdesc *current_mip)
{
   /* Display _everything_about the current problem, for debugging only*/ 
   int i = 0;			/* counter for different for loops */
   
   /* print size of the problem */
   printf("Objective Sense: ");
   putchar(current_mip->obj_sense);
   printf("\nnumber of variables: %d\n", current_mip->n);
   printf("number of constraints: %d\n", current_mip->m);
   printf("number of nonzeroes: %d\n", current_mip->nz);

   /* print which variables are integers */
   printf("is_int array: {");
   for(i=0;i<current_mip->n;i++) {
      printf("%d ", current_mip->is_int[i]);
   }
   printf("}\n");

   /* print lowerbounds */
   printf("\nlowerbounds={");
   for (i=0;i<current_mip->n;i++) {
      if (current_mip->lb[i]==-DBL_MAX) {
	 printf("(-)infinity ");
      }
      else {
	 printf("%f ",current_mip->lb[i]);
      }
   }
   printf("}\n");

   /* display upperbounds */
   printf("\nupperbounds={");
   for (i=0;i<current_mip->n;i++)
      {
	 if (current_mip->ub[i]==DBL_MAX)
	    printf("infinity ");
	 else
	    printf("%f ",current_mip->ub[i]);
      }
   printf("}\n");

   /* display RHS */
   printf("\nRHS={");
   for (i=0;i<current_mip->m;i++) {
      printf("%f ",current_mip->rhs[i]);
   }
   printf("}\n");

   /* display range values of the constraints */
   printf("\nRange values={");
   for (i=0;i<current_mip->m;i++) {
      printf("%f ",current_mip->rngval[i]);
   }
   printf("}\n");

   /* Objective function coefficients */
   printf("\nObjective fn={");
   for (i=0;i<current_mip->n;i++) {
      printf("%f ",current_mip->obj[i]);
   }
   printf("}\n");

   /* Array matbeg */
   printf("\nmatbeg={");
   for (i=0;i<current_mip->n+1;i++){
      printf("%d ",current_mip->matbeg[i]);
   }
   printf("}\n");

   /* Array matind */
   printf("\nmatind={");
   for (i=0;i<current_mip->nz;i++){
      printf("%d ",current_mip->matind[i]);
   }
   printf("}\n");

   /* Array matval */
   printf("\nmatval={");
   for (i=0;i<current_mip->nz;i++) {
      printf("%f ",current_mip->matval[i]);
   }
   printf("}\n");

   /* Array column names */
   printf ("\nColumn names: ");
   for (i=0;i<current_mip->n;i++) {
      printf(current_mip->colname[i]);
      printf(" ");
   }

   /* sense of constraints */
   printf("\nSense of constraints:");
   printf(current_mip->sense);
   printf("\n");

   /* exit */
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
int prep_display_row_mat(MIPdesc *P , rowpackedarray *row_P)
{
   /*
     display the problem using row packed matrices - only for debugging
   */

   int i; 			/* counter */

   printf("row_matbeg = {");  
   for (i=0;i<=P->m;i++) {
      printf("%d ", row_P->matbeg[i]);
   }
   printf("\b}\n");

   printf("row_matind = {");  
   for (i=0;i<P->nz;i++) {
      printf("%d ", row_P->matind[i]);  
   }

   printf("\b}\nrow_matval = {"); 
   for (i=0;i<P->nz;i++) { 
      printf("%f ", row_P->matval[i]); 
   }
   printf("\b}\n");
   
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/
