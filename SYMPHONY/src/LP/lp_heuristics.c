/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Ashutosh Mahajan                               */
/*                                                                           */
/* (c) Copyright 2006-2010 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sym_qsort.h"
#include "sym_lp.h"
#include "sym_constants.h"
#include "sym_lp_solver.h"
#include "sym_primal_heuristics.h"
#include "sym_prep.h"
#include "sym_macros.h"
#include "sym_master.h"
/*===========================================================================*/
/*===========================================================================*\
 * This file contains heuristics to find an integral solution after an LP
 * is solved.
 \*===========================================================================*/
/*===========================================================================*/
/*
 * TODO: 
 * make independent of solver
 */

int feasibility_pump (lp_prob *p, char *found_better_solution, double &solution_value, 
		      double *colSolution, double *betterSolution)
{
   int                      termcode    = FUNCTION_TERMINATED_NORMALLY;
   LPdata                  *lp_data     = p->lp_data;
   /* no. of max pumping cycles */
   int                      max_iter    = p->par.fp_max_cycles;
   int                      n           = lp_data->n;
   /* use OSI to get lp data */
   OsiSolverInterface      *model       = p->lp_data->si;
   const CoinPackedMatrix  *matrix      = model->getMatrixByRow();
   const double            *lp_r_low    = model->getRowLower();
   const double            *lp_r_up     = model->getRowUpper();
   int                      i, r, iter, cnt, verbosity;
   int                     *indices;
   double                  *values;
   double                   fp_time, last_fp_time, real_obj_value, target_ub;
   FPvars                 **vars;
   double                   gap           = model->getInfinity();
   double                   obj_lb        = lp_data->objval;
   double                   total_time    = 0;
   const double            *mip_obj       = model->getObjCoefficients();
   char                     is_feasible   = FALSE;
   double                  *x_ip, *x_lp, new_solution_value;
   const double             fp_display_interval = p->par.fp_display_interval;
   /* number of solutions with obj value more than the best */
   int                      num_poor_sols = 0;
   int                      num_better_sols = 0;
   const double             lpetol = p->lp_data->lpetol;
   int                      fp_poor_sol_lim = p->par.fp_poor_sol_lim_fac;
   int                      total_iter_cnt = 0;
   fp_time                                = used_time(&total_time);
   char                     fp_abandoned = FALSE;
   *found_better_solution = FALSE;

   if (p->lp_stat.fp_calls < 1) {
      CoinSeedRandom(17000);
   }

   int lp_iter_limit = 1e5;

   if(p->bc_index > 0){
     int reg_limit = MAX(1e4, (int)(((int)(1.0*p->tm->stat.analyzed/100) + 1)*1e8/lp_data->nz));
     lp_iter_limit = reg_limit - p->lp_stat.fp_num_iter;
   }

   if(lp_iter_limit < 0){
     if(p->has_ub || solution_value < SYM_INFINITY/2 || p->tm->stat.analyzed > 100) 
       return termcode;     
     lp_iter_limit = p->lp_stat.fp_num_iter/p->lp_stat.fp_calls;
   }

   //printf("fp : ind - iter limit %i %i\n", p->bc_index, lp_iter_limit);

   FPdata                  *fp_data     = (FPdata*) malloc(sizeof(FPdata));
   LPdata                  *new_lp_data = (LPdata *)calloc(1,sizeof(LPdata));	

   fp_data->total_iter_limit = lp_iter_limit; 
   fp_data->single_iter_limit = 1000;
   if(p->bc_index < 1){
     fp_data->single_iter_limit = 5000;
   }

   /* total_time and fp_time both now have total time used by symphony's lp
    * process */
   fp_time                                = used_time(&total_time);
   last_fp_time                           = fp_time;
   /* fp_time should now be zero and total_time be still the same */

   verbosity = fp_data->verbosity         = p->par.verbosity;
   if (p->bc_index<1) {
      PRINT(verbosity, 0, ("starting feasibility pump\n"));
   }

   fp_data->mip_obj       = (double *)malloc(n*DSIZE);
   fp_data->flip_fraction = p->par.fp_flip_fraction;
   fp_data->sos_row_filled = 0;
   fp_data->sos_var_fixed_zero = 0;
   fp_data->can_check_sos = FALSE;
   
   if(p->mip->matbeg && p->mip->mip_inf && 
      p->mip->mip_inf->binary_sos_row_num > 0){
      fp_data->can_check_sos = TRUE;
      fp_data->sos_row_filled = (char *)malloc(p->mip->m*CSIZE);
      //fp_data->sos_var_fixed_zero = (char *)malloc(p->mip->n*CSIZE);      
   }

   memcpy(fp_data->mip_obj,mip_obj,n*DSIZE);

   /* initialize the lp solver. load the current basis */
   fp_initialize_lp_solver(p, new_lp_data, fp_data, 
			   (solution_value < SYM_INFINITY/2 ? colSolution : NULL));
   x_ip = fp_data->x_ip;
   x_lp = fp_data->x_lp;

   if (p->has_ub) {
      solution_value = p->ub-p->par.granularity;
   }
   else {
      solution_value = model->getInfinity();
   }

   if (p->has_ub && p->mip->mip_inf && 
         (p->mip->mip_inf->obj_size <= p->mip->mip_inf->max_row_size || 
          p->mip->mip_inf->obj_size < n/10)) {
      fp_add_obj_row(new_lp_data, n, mip_obj, p->ub-p->par.granularity);
   } 
   /* round the x_lp and store as x_ip, it will usually become infeasible */
   vars = fp_data->fp_vars;

   /* do the following max_iter times */
   fp_time += used_time(&total_time);
   int fp_override_cnt = 0;
   /*
   if(p->lp_stat.fp_calls == 1){
      p->par.fp_time_limit += 10;
      p->par.fp_max_initial_time += 10;
   }else if (p->lp_stat.fp_calls == 2){
      p->par.fp_time_limit -= 10;    
      p->par.fp_max_initial_time -= 10;
   }
   */   
   if(p->lp_stat.fp_calls < 1){
      p->par.fp_time_limit += 20;
   }else if(p->lp_stat.fp_calls < 2){
      p->par.fp_time_limit -= 20;
      p->par.fp_max_initial_time += 20;
   }else if(p->lp_stat.fp_calls < 3){
      p->par.fp_max_initial_time -= 20;
   }
      
   for (iter=0; (iter<max_iter && fp_time<p->par.fp_time_limit &&
		 fp_time + p->comp_times.fp < p->par.fp_max_initial_time) ||
	   fp_override_cnt > 0; iter++) {
      if (fp_time - last_fp_time > fp_display_interval || verbosity > 5) {
         PRINT(verbosity, 0, 
               ("feasibility pump: starting iteration %d, time used = %.2f\n",
                iter, fp_time));
         last_fp_time = fp_time;
      }

      is_feasible = FALSE;
      /* solve an lp */
       fp_round(p, fp_data, new_lp_data);
      if (fp_data->x_bar_len[fp_data->iter] == -1) {
         /*
          * the cost and reference point are same as some other iteration. we
          * should stop here because we are cycling
          */
         PRINT(verbosity,5,("fp: leaving because of cycling\n"));
         fp_data->iter++;
         break;
      }
      fp_is_feasible (lp_data,matrix,lp_r_low,lp_r_up,fp_data,&is_feasible);

      if (is_feasible == TRUE) {
         new_solution_value = 0;
         for (i=0;i<n;i++) {
            new_solution_value += x_ip[i]*mip_obj[i];
         }
         if (new_solution_value<solution_value-p->par.granularity-lpetol) {
	    /* we found what we wanted */
	    memcpy(betterSolution, x_ip, n*DSIZE);

	    /* we found what we wanted */
	    memcpy(betterSolution, x_ip, n*DSIZE);

	    solution_value = new_solution_value;
            indices = p->lp_data->tmp.i1;          /* n */
            values  = p->lp_data->tmp.d;           /* n */
            cnt     = collect_nonzeros(p, betterSolution, indices, values);
            gap     = (solution_value -
                      obj_lb)/(fabs(solution_value)+0.001)*100;
            p->lp_stat.fp_num_sols++;
            num_better_sols++;
            PRINT(verbosity,5,("fp: found solution with value = %f\n",
                     solution_value));
            PRINT(verbosity,5,("fp: gap = %f\n", gap));
            //sp_add_solution(p,cnt,indices,values,
	    //    solution_value+p->mip->obj_offset,p->bc_index);
            if (gap <= p->par.fp_min_gap) {
               *found_better_solution = TRUE;
               fp_data->iter++;
               break;
            }
            target_ub = (obj_lb + solution_value)/2;
            if (p->mip->mip_inf && (p->mip->mip_inf->obj_size <= 
                     p->mip->mip_inf->max_row_size
                  || p->mip->mip_inf->obj_size < n/10)) {
               if (*found_better_solution != TRUE && p->has_ub==FALSE) {
                  // add another objective function constraint to lower the
                  // objective value.
                  fp_add_obj_row(new_lp_data, n, mip_obj, target_ub);
               } else {
                  r = new_lp_data->m-1;
                  change_rhs(new_lp_data, 1, &r, &target_ub);
               }
            }
            *found_better_solution = TRUE;
            fp_poor_sol_lim = p->par.fp_poor_sol_lim_fac *
                              num_better_sols;
	    /* menal ---*/ 
	    if(p->bc_level > 0) {
              fp_data->iter++;
              break;	    
            }
	    /* --- */
         } else {
            num_poor_sols++;
            /*
            PRINT(verbosity,5,("fp: rejecting poor solution with value = %f\n",
                     solution_value));
            PRINT(verbosity,5,("fp: number of poor sols = %d, better sols = %d, limit=%d\n",
                     num_poor_sols, num_better_sols, fp_poor_sol_lim));
            */
            if (num_poor_sols > fp_poor_sol_lim) {
            /*
               PRINT(verbosity,5,("fp: breaking because of too many (%d) poor"
                       " solutions\n", num_poor_sols));
            */
               fp_data->iter++;
               break;
            }
         }
      } 

      PRINT(verbosity,5,("fp: solve lp %d\n",iter));
      p->lp_stat.lp_calls++;
      p->lp_stat.fp_lp_calls++;
     
      if (fp_solve_lp(new_lp_data, fp_data, &is_feasible) != 
            FUNCTION_TERMINATED_NORMALLY) {
	fp_abandoned = TRUE;
      }

      fp_data->iter++;
      fp_time += used_time(&total_time);
      total_iter_cnt += fp_data->iterd;
      //if(fp_abandoned) break;
      if(fp_abandoned || total_iter_cnt > fp_data->total_iter_limit) break;
   }

   p->lp_stat.fp_poor_sols = num_poor_sols;
   p->lp_stat.fp_lp_total_iter_num += total_iter_cnt;
   close_lp_solver(new_lp_data);
   /* free all the allocated memory */
   FREE(new_lp_data->x);
   FREE(new_lp_data->lb);
   FREE(new_lp_data->ub);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->dualsol);
   FREE(new_lp_data->dj);
   FREE(new_lp_data->tmp.c);
   FREE(new_lp_data->tmp.d);
   FREE(new_lp_data->tmp.i1);
   FREE(new_lp_data);
   
   for (i=0;i<n;i++) {
      FREE(fp_data->fp_vars[i]);
   }
  
   for (i=0;i<fp_data->iter;i++) {
      FREE(fp_data->x_bar_val[i]);
      FREE(fp_data->x_bar_ind[i]);
   }
   FREE(fp_data->x_bar_val);
   FREE(fp_data->x_bar_ind);
   FREE(fp_data->x_bar_len);
   FREE(fp_data->fp_vars);
   FREE(fp_data->sos_row_filled);
   FREE(fp_data->sos_var_fixed_zero);
   FREE(fp_data->obj);
   FREE(fp_data->mip_obj);
   FREE(fp_data->x_lp);
   FREE(fp_data->x_ip);
   FREE(fp_data->index_list);
   FREE(fp_data->alpha_p);
   FREE(fp_data);

   /* update stats */
   fp_time                        += used_time(&total_time);
   p->comp_times.fp               += fp_time;
   p->lp_stat.fp_calls++;
   p->lp_stat.fp_last_call_ind = p->bc_index; 
   p->lp_stat.fp_num_iter += total_iter_cnt;
   if (*found_better_solution==TRUE) {
      if (p->mip->obj_sense == SYM_MAXIMIZE){
         real_obj_value=-solution_value+p->mip->obj_offset;
      } else {
         real_obj_value=solution_value+p->mip->obj_offset;
      }
      PRINT(verbosity,5,("fp: found solution = %10.2f time = %10.2f\n",
               real_obj_value,total_time));
   }

   if (p->bc_index<1 || verbosity > 5) {
      PRINT(verbosity, 0, ("leaving feasibility pump.\n"));
   }
 
   return termcode;
}


/*===========================================================================*/
int fp_is_feasible (LPdata *lp_data, const CoinPackedMatrix *matrix,
		    const double *r_low, const double *r_up, FPdata *fp_data,
		    char *is_feasible )		    
{
   /* check if x is a integer feasible solution to problem in p */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double lpetol = lp_data->lpetol;
   //int n = fp_data->n0;
   int m = fp_data->m0;
   //FPvars **vars = fp_data->fp_vars;
   int i,c,j;
   double Ractivity;
   const int *r_matbeg = matrix->getVectorStarts();
   const int *r_matlen = matrix->getVectorLengths();
   const int *r_matind = matrix->getIndices();
   const double *r_matval = matrix->getElements();
   double *x = fp_data->x_ip;

   *is_feasible = TRUE;
   /* some int variable is non-integral */
   /* is not possible, since this function is called after rounding */

   /* check feasibility of constraints */
   for (i=0;i<m;i++) {
      Ractivity = 0;
      c=0;			/* column */
      for (j=r_matbeg[i];j<r_matbeg[i]+r_matlen[i];j++) {
         c=r_matind[j];
         Ractivity += x[c]*r_matval[j];
      }
      //      printf("Ractivity[%d] = \t%f\n",i,Ractivity);
      if (Ractivity>r_up[i]+lpetol || Ractivity<r_low[i]-lpetol) {
         /* constraint infeasibility is possible since we call this func. after
            rounding */
         *is_feasible = FALSE;
         //printf("constraint %d activity = %f, down = %g, up = %g\n",
	 //i, Ractivity, r_low[i], r_up[i]);
         break;
      }
   }

   return termcode;
}

/*===========================================================================*/
int fp_initialize_lp_solver(lp_prob *p, LPdata *new_lp_data, FPdata *fp_data, 
			    double *colSolution)
{
   /*
      create a copy of lp_data into new_lp_data
      for general mixed int programs, we will have to add 2 new vars for each
      non-binary integer var. (x_j+ and x_j-)
   */

   /* first create an exact copy of lp_data */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   LPdata *lp_data  = p->lp_data;
   new_lp_data->lpetol = lp_data->lpetol;
   int n = lp_data->n;
   int m = lp_data->m;
   int i, k, *outrhsind;
   //int *rstat,*cstat;

   double one=1.0;
   char sense='G';
   char where_to_move='E';	/* redundant */
   int col_number = n;
   int *rmatbeg = (int *) malloc(2*ISIZE);
   int *cmatbeg = (int *) malloc(2*ISIZE);
   int *rmatind = (int *) malloc(3*ISIZE);
   double *rmatval = (double *) malloc(3*DSIZE);
   int *cmatind = NULL;
   double *cmatval = NULL;
   double rhs;
   double lb, ub;
   double lpetol = lp_data->lpetol;
   double *lp_lb, *lp_ub, *fp_obj;
   double norm_c = 0;
   double *mip_obj = fp_data->mip_obj;
   int verbosity = fp_data->verbosity;
   int *index_list;
   int fp_max_length_cuts = 1; 
   row_data *rows = lp_data->rows;
  
   /* used because we can not call si directly */
   copy_lp_data(lp_data,new_lp_data);
   new_lp_data->si->setupForRepeatedUse(3,0); 

#ifdef COMPILE_IN_LP
   if(p->mip->matbeg){
     double mat_den = (1.0)*p->mip->nz/(p->mip->m * p->mip->n + 1);
     if(p->mip->nz > 1e5 && mat_den > 0.01){
       new_lp_data->si->setupForRepeatedUse(0,0); 
     }
   }
#endif

   if(p->par.fp_fix_ratio > 0.0 && p->mip->mip_inf){
     double *x = lp_data->x;
     double *x_rank = lp_data->tmp.d; 
     double *x_rank2 = lp_data->tmp.d + n; 
     int ind, *x_ind = lp_data->tmp.i1;
     int *x_ind2 = lp_data->tmp.i1;
     int vars_eff_cnt, int_cnt = 0, int_cnt2 = 0;
     double bd, x_obj, min_obj = DBL_MAX, big_number = 1e20; 
     double * obj = const_cast <double *> (lp_data->si->getObjCoefficients());  
     for(i = 0; i < n; i++){
       if(obj[i] < min_obj) min_obj = obj[i];
     }

     min_obj = fabs(min_obj);

     for(i = 0; i < n; i++){
       if(lp_data->vars[i]->is_int && lp_data->ub[i] > lp_data->lb[i] + lpetol){
	 if(obj[i] >= 0.0 && x[i] < lp_data->lb[i] + lpetol){
	   x_obj = obj[i] + min_obj + 1e-4; 
	   vars_eff_cnt = MAX(p->mip->mip_inf->cols[i].sos_num,  p->mip->mip_inf->cols[i].col_size) + 1;
	   x_rank[int_cnt] = big_number*x_obj/vars_eff_cnt; 
	   if(colSolution && x[i] < colSolution[i] + lpetol && x[i] > colSolution[i] - lpetol){
	     x_rank[int_cnt] /= 2;
	   }
	   x_ind[int_cnt] = i;
	   int_cnt++;
	 }else if(obj[i] <= 0.0 && x[i] > lp_data->ub[i] - lpetol){
	   x_obj = obj[i] + min_obj + 1e-4; 
	   vars_eff_cnt = MAX(p->mip->mip_inf->cols[i].sos_num,  p->mip->mip_inf->cols[i].col_size) + 1;
	   x_rank2[int_cnt2] = big_number*x_obj/vars_eff_cnt; 
	   if(colSolution && x[i] < colSolution[i] + lpetol && x[i] > colSolution[i] - lpetol){
	     x_rank2[int_cnt2] /= 2;
	   }
	   x_ind2[int_cnt2] = i;
	   int_cnt2++;	   
	 }
       }
     }
	
     qsort_di(x_rank, x_ind, int_cnt);
     qsort_di(x_rank2, x_ind2, int_cnt2);
    
     int fix_cnt = MIN((int)(0.5*int_cnt), (int)(p->par.fp_fix_ratio*int_cnt));
     int fix_cnt2 = MIN((int)(0.5*int_cnt2), (int)(p->par.fp_fix_ratio*int_cnt2));
     
     //printf("F-cnt : %i\n", fix_cnt);
     for(i = 0; i < fix_cnt; i++) {
       ind = x_ind[i];
       bd = floor(x[ind] + lpetol);
       change_lbub(new_lp_data, ind, bd ,bd);
     }     
     for(i = 0; i < fix_cnt2; i++) {
       ind = x_ind2[i];
       bd = floor(x[ind] + lpetol);
       change_lbub(new_lp_data, ind, bd ,bd);
     }     
   }
   
   lp_lb = new_lp_data->lb;
   lp_ub = new_lp_data->ub;

   /* delete cuts that are long as they slow down the lp */
   outrhsind = (int *)calloc(m, ISIZE);
   k = 0;

#ifdef COMPILE_IN_LP   
   if(p->bc_level < 1 && p->mip->mip_inf && p->mip->mip_inf->cont_var_num <= 0){
      fp_max_length_cuts = 100;   
   }
#endif
   
   for (i = p->base.cutnum; i < m; i++){
      if (((int *)rows[i].cut->coef)[0] > fp_max_length_cuts) {
         outrhsind[k] = i;
         k++;
      }
   }
   PRINT(verbosity, 5, ("feasibility pump: cuts discarded = %d\n", k));
   delete_rows_with_ind(new_lp_data, k, outrhsind);
   m -= k;
   //   printf("m: %i \n",m);
   /* set up fp_data */
   fp_data->alpha           = 0.8;
   fp_data->alpha_decr      = 0.7;
   fp_data->n0 = fp_data->n = n;
   
   fp_data->m0              = m;
   fp_data->iter            = 0;

   /* count how many binary variables */
   fp_data->fp_vars         = (FPvars **) malloc(sizeof(FPvars *)*n);
   fp_data->x_ip            = (double *) calloc(n,DSIZE);
   fp_data->x_lp            = (double *) calloc(n,DSIZE);
   fp_data->index_list      = (int *)    calloc(n,DSIZE);
   fp_data->x_bar_ind       = (int **)   calloc(p->par.fp_max_cycles,
                                                sizeof(int*));
   fp_data->x_bar_val       = (double **)calloc(p->par.fp_max_cycles,
                                                sizeof(double*));
   fp_data->x_bar_len       = (int *)    calloc(p->par.fp_max_cycles,ISIZE);
   fp_data->alpha_p         = (double *) malloc(p->par.fp_max_cycles*DSIZE);
   FPvars **fp_vars         = fp_data->fp_vars;
   fp_data->numNonBinInts   = 0;
   fp_data->numInts         = 0;

   index_list = fp_data->index_list;
   for (i=0;i<n;i++) {
      index_list[i]=i;
      fp_vars[i] = (FPvars *)malloc(sizeof(FPvars));
      if (lp_data->vars[i]->is_int) {
         fp_data->numInts++;
         fp_vars[i]->is_int = TRUE;
         if (lp_data->vars[i]->lb<-lpetol||lp_data->vars[i]->ub>1+lpetol) {
            fp_vars[i]->is_bin = FALSE;
            fp_data->numNonBinInts++;
         }
         else {
            fp_vars[i]->is_bin = TRUE;
         }
      } else {
         fp_vars[i]->is_int = fp_vars[i]->is_bin = FALSE;
      }
      /* calculate ||C|| */
      norm_c += mip_obj[i]*mip_obj[i];
   }
   
   norm_c = sqrt(norm_c);
   PRINT(verbosity, 20, ("fp: norm_c = %f\n",norm_c));

   fp_data->n       = n+fp_data->numNonBinInts;
   fp_data->m       = m+2*fp_data->numNonBinInts;
   fp_data->obj     = (double *)malloc(fp_data->n*DSIZE);
   new_lp_data->x   = (double *)calloc(fp_data->n,DSIZE);
   memcpy(fp_data->x_lp,p->lp_data->x,DSIZE*n);

   if (norm_c>lpetol) {
      for (i=0;i<n;i++) {
         mip_obj[i] = mip_obj[i]/norm_c;
      }
   }
   
   /* load basis */
   //rstat = (int *) malloc(m * ISIZE);
   //cstat = (int *) malloc(n * ISIZE);

   //get_basis(lp_data,cstat,rstat);
   //load_basis (new_lp_data,cstat,rstat);

   //FREE(rstat);
   //FREE(cstat);

   /* add 1 columns and 2 rows for each nonBinary Integer */
   /* 
    * min d_i
    * s.t.
    * d_i - x_i >= -x_i^0
    * d_i + x_i >=  x_i^0
    */
   rmatbeg[0] =  0;
   rmatbeg[1] =  2;
   cmatbeg[0] =  0;
   cmatbeg[1] =  0;
   rmatval[0] =  1.0;
   lb         = -SYM_INFINITY;
   ub         =  SYM_INFINITY;
   fp_obj     =  fp_data->obj;

   /* used by change_rhs*/
   new_lp_data->tmp.c = (char *)malloc(2*CSIZE);
   new_lp_data->tmp.d = (double *)malloc(DSIZE*(n+2)); /* +2 for add_rows */
   new_lp_data->tmp.i1 = (int *)malloc(ISIZE*n);

   for (i=0;i<n;i++) {
      if (fp_vars[i]->is_int && !fp_vars[i]->is_bin) {
         /* add d_i */
	 add_cols(new_lp_data, 1, 0, &one, cmatbeg, cmatind, cmatval, &lb, &ub, 
               &where_to_move);
         fp_vars[i]->xplus = col_number;

         /* now add two rows */
         /* d_i - x_i >= -x_i^0 */
         rhs        = -1*lp_data->x[i];
         rmatind[0] =  col_number;
         rmatind[1] =  i;
         rmatval[1] = -1.0;
         add_rows(new_lp_data, 1, 2, &rhs, &sense, rmatbeg, rmatind, rmatval);

         /* d_i - x_i >= -x_i^0 */
         rhs = lp_data->x[i];
         rmatval[1] = 1.0;
         add_rows(new_lp_data, 1, 2, &rhs, &sense, rmatbeg, rmatind, rmatval);
         
         fp_obj[col_number] = 1.0;
         col_number++;
      }
   }

   FREE(rmatval);
   FREE(rmatind);
   FREE(cmatbeg);
   FREE(rmatbeg);
   FREE(outrhsind);

   return termcode;
}

/*===========================================================================*/
int fp_solve_lp(LPdata *lp_data, FPdata *fp_data, char* is_feasible) 
{
   /* construct an lp based on x_ip. solve it. store the result in x_lp */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double *objcoeff= fp_data->obj;
   int n = fp_data->n;
   //int iterd;
   int termstatus;
   int i;
   double delta_x;
   double norm = 0;
   FPvars **fp_vars = fp_data->fp_vars;
   double *mip_obj  = fp_data->mip_obj;
   int verbosity = fp_data->verbosity;
   int  *index_list = fp_data->index_list;
   double *x_ip = fp_data->x_ip;
   double *x_lp = fp_data->x_lp;
   double alpha = fp_data->alpha;
   double one_minus_alpha = 1-fp_data->alpha;
   int n0 = fp_data->n0;
   double *lp_data_x = lp_data->x;
   double etol = lp_data->lpetol;
      
   is_feasible = FALSE;
   memset ((char *)(objcoeff),0,DSIZE*n);
   for (i=0;i<n0;i++) {
      if (fp_vars[i]->is_int) {
         if (fp_vars[i]->is_bin) {
            if (x_ip[i] <= 0.0 + etol && x_ip[i] >= 0.0 - etol) {
               objcoeff[i] = 10.0;
	    } else if (x_ip[i] >= 1.0 - etol && x_ip[i] <= 1.0 + etol ) {
	       objcoeff[i] = -10.0;
            }
         } else {
            objcoeff[i] = 0.0;
            objcoeff[fp_vars[i]->xplus] = 1.0;
         }
      } else {
         objcoeff[i] = 0.0;
      }
      /* calculate ||coeff||, norm is not zero because otherwise x_ip is
       * feasible */
   }

   if (fp_data->iter < 1) {
      norm = 0;
      for (i=0; i < n0; i++) {
         norm += objcoeff[i]*objcoeff[i]; /* stays the same every iteration */
      }
      norm = sqrt(norm);
      fp_data->norm = norm;
   } else {
      norm = fp_data->norm;
   }

   //norm = 0;
   PRINT(verbosity, 15, ("fp: norm = %f\n",norm));
   for (i=0;i<n0;i++) {
      objcoeff[i] = 
      one_minus_alpha*objcoeff[i]+alpha*mip_obj[i]*norm;
   }
  /*
   for (i=fp_data->n0;i<fp_data->n;i++) {
      objcoeff[i] = (1-alpha)*objcoeff[i];
   }
   alpha = alpha*fp_data->alpha_decr;
   for (i=0;i<n0;i++) {
      if (fp_vars[i]->is_int) {
         lp_data->si->setInteger(i);
      }
   }
   */

   change_objcoeff(lp_data, index_list, &index_list[n-1], objcoeff);
   if (fp_data->iter > 0) { 
     set_itlim(lp_data, fp_data->single_iter_limit);   
     termstatus = dual_simplex(lp_data, &fp_data->iterd);
   } else {
     set_itlim(lp_data, 5*fp_data->single_iter_limit);   
     termstatus = initial_lp_solve(lp_data, &fp_data->iterd);
   }
   //printf("iter - %i\n", fp_data->iterd); 

   if (termstatus != LP_OPTIMAL) {
     //PRINT(verbosity,0,("Feasibility Pump: Unable to solve LP. Pump malfunction.\n"));
      return FUNCTION_TERMINATED_ABNORMALLY;
   }

   get_x(lp_data);

   delta_x = 0;
   memcpy(x_lp,lp_data_x,DSIZE*n0);

   /*
   for (i=0;i<n0;i++) {
      if (fp_vars[i]->is_int) {
         delta_x += fabs(x_lp[i]-x_ip[i]);
      }
   }
   PRINT(verbosity, 15, ("fp: delta_x = %f\n",delta_x));
   */

   return termcode;
}

/*===========================================================================*/
int fp_add_obj_row(LPdata *new_lp_data, int n, const double *obj, double rhs)
{
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   char sense = 'L';
   int *rmatbeg, *rmatind;
   double *rmatval;
   int i, count, nz;
   double lpetol = new_lp_data->lpetol;

   // count non zeros
   // we dont trust p->mip->mip_inf->obj_size because it is the size before
   // preprocessing.
   nz = 0;
   for (i=0;i<n;i++) {
      if (fabs(obj[i])>lpetol) {
         nz++;
      }
   }

   rmatbeg = (int *) malloc(2*ISIZE);
   rmatind = (int *) malloc(nz*ISIZE);
   rmatval = (double *) malloc(nz*DSIZE);

   count = 0;
   for (i=0;i<n;i++) {
      if (fabs(obj[i])>lpetol) {
         rmatval[count] = obj[i];
         rmatind[count] = i;
         count++;
      }
   }
   rmatbeg[0] = 0;
   rmatbeg[1] = nz;
   add_rows(new_lp_data, 1, nz, &rhs, &sense, rmatbeg, rmatind, rmatval);
   FREE(rmatbeg);
   FREE(rmatind);
   FREE(rmatval);
   return termcode;
}

/*===========================================================================*/
int fp_round(lp_prob *p, FPdata *fp_data, LPdata *lp_data)
{
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double *x_ip = fp_data->x_ip;
   double *x_lp = fp_data->x_lp;
   int i,j, has_changed;
   int n = fp_data->n0;
   double lpetol = lp_data->lpetol;
   int *tind = lp_data->tmp.i1; /* n */
   double *tx = lp_data->tmp.d; /* n */
   int cnt = 0;
   int *index = fp_data->index_list;
   double **x_bar_val_p = fp_data->x_bar_val;
   double *x_bar_val;
   int **x_bar_ind_p = fp_data->x_bar_ind;
   int *x_bar_ind;
   int *x_bar_len = fp_data->x_bar_len;
   double flip_fraction = fp_data->flip_fraction;
   FPvars **vars = fp_data->fp_vars;
   int fp_iter = fp_data->iter;
   double *alpha_p = fp_data->alpha_p;
   int sos_row_filled_cnt = 0;
 
   if(fp_data->can_check_sos){
      memset(fp_data->sos_row_filled, 0, CSIZE*p->mip->m);
      //memset(fp_data->sos_var_fixed_zero, 0, CSIZE*p->mip->n); 
   }
   
   for (i=0;i<n;i++) {
      if (vars[i]->is_int) {
         /* round x_lp[i] and put into x_ip[i] */
         x_ip[i]=floor(x_lp[i]+0.5);
	 /*
	 if(vars[i]->is_bin && fp_data->can_check_sos && x_ip[i] == 1.0 && 
	    p->mip->mip_inf->cols[i].sos_num){
	    if(fp_data->sos_var_fixed_zero[i]) x_ip[i] = 0;
	    else fp_fix_sos_var(p, fp_data, i);
	 }
	 */
	 if(vars[i]->is_bin && fp_data->can_check_sos && x_ip[i] == 1.0 && 
	    p->mip->mip_inf->cols[i].sos_num){
	    if(!(fp_can_sos_var_fix(p, fp_data, i, &sos_row_filled_cnt))){
	       x_ip[i] = 0.0;
	    }
	 }	 
      }
      else {
         x_ip[i]=x_lp[i];
      }
   }

   // TODO: make it work for '0'
   //       remove randomness
   while (1) {
      cnt = 0;
      for (i = 0; i < n; i++){
         if (vars[i]->is_int && (x_ip[i] > lpetol || x_ip[i] < -lpetol)){
            tind[cnt] = index[i];
            tx[cnt++] = x_ip[i];
         }
      }
      /* order indices and values according to indices */
      qsort_id(tind, tx, cnt);

      /* go through all 'iter' points and check if x_ip already exists */
      has_changed = TRUE;
      for (i=0; i<fp_iter; i++) {
         //printf("alpha = %f, len = %d\n", alpha_p[i], x_bar_len[i]);
         if (x_bar_len[i] == cnt && alpha_p[i] < 0.08) {
            x_bar_val = x_bar_val_p[i];
            x_bar_ind = x_bar_ind_p[i];
            for (j=0; j<cnt; j++) {
               if (tind[j]!=x_bar_ind[j] || fabs(tx[j]-x_bar_val[j])>lpetol) {
                  break;
               }
            }
            if (j==cnt) {
               PRINT(fp_data->verbosity,5,("fp: same as %d\n",i));
               break; //its same
            }
         }
      }

      if (i<fp_iter) {
         /* flip some vars in x_ip */	 
	 //if(fp_data->can_check_sos){
	 //  memset(fp_data->sos_row_filled, 0, CSIZE*p->mip->m); 
	 //  sos_row_filled_cnt = 0;	 
	 //}
	 
         int num_flipped = 0;
	 
         has_changed = FALSE;
         PRINT(fp_data->verbosity,5,("fp: flipping\n"));

	 for (j=0; j<n; j++) {
	    if (vars[j]->is_bin) {
	       
	       if (CoinDrand48()<flip_fraction) {
		  x_ip[j] = 1-x_ip[j];
		  num_flipped++;
	       }
	       // if(fp_data->can_check_sos && x_ip[j] == 1.0 && 
	       //  p->mip->mip_inf->cols[j].sos_num){
		  //if(!(fp_can_sos_var_fix(p, fp_data, j, &sos_row_filled_cnt))){
		  //   x_ip[j] = 0.0;
		  // }
	       // }	       
	    } else if (vars[j]->is_int) {
	       if (CoinDrand48()<flip_fraction) {
		  x_ip[j] = floor(x_lp[j]) + 
		     floor(ceil(x_lp[j]) - x_lp[j] + 0.5); /*round and flip*/
	       }
	    }	    
	 }

	 PRINT(fp_data->verbosity,5,("fp: flipping %d\n", num_flipped));
         if (num_flipped==0) {
            // TODO: dont know what to do
            break;
         }
      } else {
         break;
      }
   }

   /*
   int k;
   if(fp_data->can_check_sos && p->mip->mip_inf->binary_sos_row_num > sos_row_filled_cnt){
      int fix_col = 0;
      int row_ind = 0;
      for(k = 0; k < p->mip->m; k++){
	 if(p->mip->mip_inf->rows[k].is_sos_row &&
	    !(fp_data->sos_row_filled[k])){
	    fix_col = p->mip->row_matind[p->mip->row_matbeg[k]];
	    for(j = p->mip->matbeg[fix_col]; j < p->mip->matbeg[fix_col + 1];
		j++){
	       row_ind = p->mip->matind[j];
	       if(p->mip->mip_inf->rows[row_ind].is_sos_row){		     
		  fp_data->sos_row_filled[row_ind] = TRUE;
		  sos_row_filled_cnt++;
	       }
	    }
	    x_ip[fix_col] = 1.0;
	    if(sos_row_filled_cnt >= p->mip->mip_inf->binary_sos_row_num){
	       break;
	    }
	 }
      }
   }
   */

   if (has_changed==TRUE || fp_data->alpha>0) {
      fp_data->x_bar_ind[fp_iter] = (int *)malloc(ISIZE*cnt);
      fp_data->x_bar_val[fp_iter] = (double *)malloc(DSIZE*cnt);
      x_bar_len[fp_iter] = cnt;
      memcpy(fp_data->x_bar_ind[fp_iter],tind,ISIZE*cnt);
      memcpy(fp_data->x_bar_val[fp_iter],tx,DSIZE*cnt);
      fp_data->alpha = fp_data->alpha*fp_data->alpha_decr;
      if (fp_data->alpha<0.08) {
         fp_data->alpha = 0;
      }
      fp_data->alpha_p[fp_iter] = fp_data->alpha;
   } else {
      x_bar_len[fp_iter] = -1;
   }
   return termcode;
}
/*===========================================================================*/

int fp_fix_sos_var(lp_prob *p, FPdata *fp_data, int ind)
{

   int k, j, row_ind, col_ind;
   for(k = p->mip->matbeg[ind]; k < p->mip->matbeg[ind+1]; k++){
      row_ind = p->mip->matind[k];
      for(j = p->mip->row_matbeg[row_ind + 1] - 1; j >= p->mip->row_matbeg[row_ind] ; j--){
	 col_ind = p->mip->row_matind[j];
	 if(col_ind <= ind) break;
	 else fp_data->sos_var_fixed_zero[col_ind] = TRUE;
      }
   }

   return 0;
}

/*===========================================================================*/

int fp_can_sos_var_fix(lp_prob *p, FPdata *fp_data, int ind, int *filled_row_cnt)
{
   int k, row_ind;
   
   for(k = p->mip->matbeg[ind]; k < p->mip->matbeg[ind+1]; k++){
      row_ind = p->mip->matind[k];
      if(p->mip->mip_inf->rows[row_ind].is_sos_row){
	 if(fp_data->sos_row_filled[row_ind]){
	    return FALSE;
	 }			   
      }
   }
   for(k = p->mip->matbeg[ind]; k < p->mip->matbeg[ind+1]; k++){
      row_ind = p->mip->matind[k];
      if(p->mip->mip_inf->rows[row_ind].is_sos_row){
	 fp_data->sos_row_filled[row_ind] = TRUE; 
	 (*filled_row_cnt)++;
      }
   }

   return TRUE; 
}
/*===========================================================================*/
int fp_should_call_fp(lp_prob *p, int branching, int *should_call, 
		      char is_last_iter, double t_lb)
{
   int        termcode = FUNCTION_TERMINATED_NORMALLY;
   
   *should_call = FALSE;
   if (is_last_iter==FALSE || (p->has_ub && p->lp_stat.fp_calls > 100)){
			       //(p->ub-p->lp_data->objval)/(fabs(p->ub)+0.0001)*100 >
			       //2*p->par.fp_min_gap)){
      return termcode;
   }

   int fp_freq_base = p->bc_level;
#ifdef COMPILE_IN_LP
   //   fp_freq_base = p->tm->stat.analyzed - 1;
#endif

   int orig_fp_freq = p->par.fp_frequency;
   if(!p->has_ub && p->lp_stat.fp_calls < 3 &&
      p->lp_stat.lp_total_iter_num/(p->lp_stat.lp_calls -
				    p->lp_stat.str_br_lp_calls -
				    p->lp_stat.fp_lp_calls + 1) > 1000){
      p->par.fp_frequency = 5;
   }

   if (p->par.fp_enabled>0 && !branching) {
      if (p->par.fp_enabled == SYM_FEAS_PUMP_REPEATED && 
	  (fp_freq_base)%p->par.fp_frequency==0) {
         *should_call = TRUE;
      } else if (p->has_ub==FALSE && p->par.fp_enabled==SYM_FEAS_PUMP_TILL_SOL
            && p->bc_level%p->par.fp_frequency==0) {
         *should_call = TRUE;
      } else if (  (p->has_ub==FALSE|| 
		    (p->ub-t_lb)/(fabs(p->ub)+0.0001)*100>
		    p->par.fp_min_gap) &&
		   (p->comp_times.fp < p->par.fp_max_initial_time) &&// ||
		   //p->lp_stat.fp_lp_total_iter_num <
		   //0.025*(p->lp_stat.lp_total_iter_num +
		   //   p->lp_stat.fp_lp_total_iter_num)) &&
		   //p->comp_times.fp < 0.025*p->tt) &&
		   //p->comp_times.fp < 0.5*p->tt && //menal, also index->level
		   fp_freq_base%p->par.fp_frequency == 0){// && 
	           //(int)(0.1*1e9/p->lp_data->nz) - p->lp_stat.fp_num_iter > 0 ) {
         *should_call = TRUE;
      }
   }
   
   if(p->bc_level < 1 && p->lp_stat.fp_calls > 0 &&
      p->comp_times.fp >= 0.5*p->par.fp_time_limit){
      *should_call = FALSE;
   }else if (!should_call){
      if(p->bc_level > 0 && !p->has_ub && 
	 (p->lp_stat.fp_calls <= 3 || p->tm->stat.analyzed >= 100)){
	if(p->lp_stat.fp_calls <= 3){
	  *should_call = TRUE;
	}else{
	  if(p->tm->stat.analyzed%p->par.fp_frequency == 0 && p->tm->stat.analyzed <= 1000){
	    *should_call = TRUE;
	    if((p->tm->stat.analyzed - 
		(int)(1.0*p->tm->stat.analyzed/100)*100)/p->par.fp_frequency <= 1){
	      p->par.fp_max_initial_time += 20; 
	    }
	  }
	}
      }
   }
      
   p->par.fp_frequency = orig_fp_freq;
   
   if (*should_call == TRUE) {
      p->lp_stat.num_fp_calls_in_path++;
   }
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

int diving_search(lp_prob *p, double *solutionValue, double *colSolution,
		  double *betterSolution, char is_last_iter, double t_lb)
{
  
  int i, iter_cnt = 0, lp_iter_limit = 0;
  LPdata *lp_data = p->lp_data;
  LPdata *diving_lp = (LPdata *)calloc(1,sizeof(LPdata)); 
  double *x = NULL;//lp_data->d;
  //int fixed_n = 0, *fixed_cols = NULL;
  double incr_ratio = p->par.ds_incr_ratio; 
  double ip_solve_col_ratio = p->par.ds_solve_ip_col_ratio; 
  double ip_solve_gap = p->par.ds_solve_ip_min_gap;
  int fix_incr_cnt = 0;
  //double etol = lp_data->lpetol; 
  //double etol100 = etol*100;
  int n = lp_data->n;
  int has_ub = FALSE;
  double ub, lb, obj_ub = 0;
  int d_cnt = 0; 
  int diving_type[DIVING_HEURS_CNT];
  int d_type;

  if(*solutionValue < SYM_INFINITY/2){
    has_ub = TRUE;
    obj_ub = *solutionValue;
  }else if((has_ub = p->has_ub)) obj_ub = p->ub;

  double dual_gap = 100;

  if(has_ub){
    dual_gap = d_gap(obj_ub, t_lb, 
		     p->mip->obj_offset, p->mip->obj_sense);
  }

  int reg_limit;
  if(p->tm->stat.analyzed < 500){
    reg_limit = MAX(1e4, (int)(((int)(1.0*p->tm->stat.analyzed/100) + 1)*1e8/lp_data->nz));
  }else{
    reg_limit = MAX(1e4, (int)(0.5e9/lp_data->nz));
  }

  lp_iter_limit = reg_limit - p->lp_stat.ds_num_iter;
  char iter_extended = FALSE; 
  if(lp_iter_limit < 0){
    if((has_ub && dual_gap < 10.0) || p->lp_stat.ds_calls > 20 || p->par.rs_mode_enabled)
      return FALSE;
    lp_iter_limit = p->lp_stat.ds_num_iter/p->lp_stat.ds_calls;
    iter_extended = TRUE;
  }
  
  //lp_iter_limit = 10*(p->lp_stat.lp_total_iter_num/(p->lp_stat.lp_calls -
  //					    p->lp_stat.str_br_lp_calls -
  //					    p->lp_stat.fp_lp_calls + 1) + 1);  
  
  /* fix-me \ get rid of some of the cuts? */

  char is_ip_feasible = FALSE, is_return_feasible = FALSE; 
  int tot_lp_iter, dive_depth, dive_depth_limit; 
  //int * is_fixed = (int*)malloc(ISIZE*n);
  int ip_vars_cnt = 0, frac_ip_cnt = 0, int_ip_cnt = 0, init_int_ip_cnt;
  int no_impr_cnt, no_prog_cnt, init_frac_ip_cnt,
     prev_frac_ip_cnt, min_frac_ip_cnt; 
  int no_better_cnt;
  //double force_lp_ratio = 0.05; //ds_solve_lp_col_ratio 
  char can_iterate;
  //double valuesi;
  sym_environment * env = NULL;

  //double *init_fix_weights = lp_data->tmp.d;
  //char *direction = lp_data->tmp.c;  

  //ub = (double*)malloc(DSIZE*n);
  //lb = (double*)malloc(DSIZE*n);

  int * init_frac_ind = lp_data->tmp.i1;
  int * frac_ind = lp_data->tmp.i1 + n;

  ds_get_frac_vars(lp_data, lp_data->x, init_frac_ind, &init_frac_ip_cnt,
		   &init_int_ip_cnt);  
  if(has_ub){
     if(p->bc_level > 0 && init_frac_ip_cnt > 1000) return FALSE; /* not worth it */
     if(!is_last_iter && init_frac_ip_cnt > 100) return FALSE; /* not worth it */
  }else{
     //if(p->bc_level > 0 && init_frac_ip_cnt > 2500) return FALSE; /* not worth it */
     if(!is_last_iter && init_frac_ip_cnt > 500) return FALSE; /* not worth it */
  }

  int d_factor = 1;
  char check_init = TRUE;
  char check_fix = TRUE;
  if(!iter_extended)
    d_factor = (int)(dual_gap/5.0) + 1;

  //double fixable_ratio = 1.0*init_int_ip_cnt/(init_int_ip_cnt +init_frac_ip_cnt);
  if((has_ub && dual_gap > 10.0) || !is_last_iter || (p->bc_index < 1 && init_frac_ip_cnt > 100) || 
     (p->bc_index >= 1 && init_frac_ip_cnt > 50) || p->par.rs_mode_enabled)
     check_init = FALSE; /*probably isn't worth*/
  
  //if(init_frac_ip_cnt > 1000){
  // if((p->bc_index < 1 && fixable_ratio < 0.5) || 
  //   (p->bc_index >= 1 && fixable_ratio < 0.8)) check_fix = FALSE; 
  //}
    
  for(d_type = 0; d_type < DIVING_HEURS_CNT; d_type++){
    switch(d_type){
    case VLENGTH_FIX_DIVING:
    case GUIDED_FIX_DIVING:
    case CROSSOVER_FIX_DIVING:
      if(check_fix && 
	 (p->lp_stat.ds_type_num_sols[d_type] > 0 ||
	  p->lp_stat.ds_type_calls[d_type] < 5*d_factor))
        diving_type[d_cnt++] = d_type;
      break;
    case EUC_FIX_DIVING:
    case RANK_FIX_DIVING:
      if(!check_fix || p->tm->stat.analyzed < 10 || p->par.rs_mode_enabled) break;
    case FRAC_FIX_DIVING:
      if(check_fix && 
	 (p->lp_stat.ds_type_num_sols[d_type] > 0 ||
	  p->lp_stat.ds_type_calls[d_type] < 5*d_factor) && !p->par.rs_mode_enabled)
	 diving_type[d_cnt++] = d_type;
      break;
    case VLENGTH_DIVING:
    case GUIDED_DIVING:
    case CROSSOVER_DIVING:
      if(check_init && 
	 (p->lp_stat.ds_type_num_sols[d_type] > 0 ||
	  p->lp_stat.ds_type_calls[d_type] < 2*d_factor))
        diving_type[d_cnt++] = d_type;
      break;
    case EUC_DIVING:
    case RANK_DIVING:
      if(!check_init || p->tm->stat.analyzed < 5) break;
    case FRAC_DIVING:
      if(check_init && 
	 (p->lp_stat.ds_type_num_sols[d_type] > 0 ||
	  p->lp_stat.ds_type_calls[d_type] < d_factor))
        diving_type[d_cnt++] = d_type;
      break;
    default:
      break;
    }
  }

  if(d_cnt < 1) return FALSE; 

  double start_time = wall_clock(NULL);
  double mark_time = 0;
  
  lp_iter_limit = (int)(1.0*lp_iter_limit/d_cnt) + 1;

  //if(p->bc_level > 0){// &&
     //(init_frac_ip_cnt > 200 || p->mip->mip_inf->prob_type == BINARY_TYPE)){
  //  if(!((has_ub && p->tm->stat.analyzed < 501) ||  
  //  (!has_ub && p->tm->stat.analyzed < 1001)))
  //return FALSE;
  //}
  for (i = 0; i < n; i++){
    if (lp_data->vars[i]->is_int) {
      ip_vars_cnt++;
    }
    ///get_lb(lp_data, i, &lb[i]);
    //get_ub(lp_data, i, &ub[i]);
  }

  int *cstat = lp_data->tmp.i1 + 2*n;//(int*)calloc(ISIZE*n);
  int *rstat = lp_data->tmp.i2;//(int*)calloc(ISIZE*m);
  //int *cstat = (int*)malloc(ISIZE*n);
  //int *rstat = (int*)malloc(ISIZE*lp_data->m);
  
  get_basis(lp_data, cstat, rstat);  
  copy_lp_data(lp_data, diving_lp);
  load_basis(diving_lp, cstat, rstat);
  //diving_lp->x = (double*)malloc(DSIZE*n);
  diving_lp->x = lp_data->tmp.d + 3*n; //betterSolution;
  
  //incr_ratio = MIN((p->bc_level % p->par.ds_frequency) * 0.025, 0.3);

  fix_incr_cnt = (int)(init_frac_ip_cnt*incr_ratio);// + p->bc_level % p->par.ds_frequency;
  
  //double frac_ratio = 1.0*init_frac_ip_cnt/ip_vars_cnt;

  //if(frac_ratio > 0.05){
  //}

  dive_depth_limit = 200;

  //int min_depth_limit = 100;
  //int max_depth_limit = 400;
  int no_impr_cnt_limit = 10; //ds_no_impr_cnt_limit
  int no_prog_cnt_limit = 10; //ds_no_impr_cnt_limit
  int no_better_cnt_limit = 10; //ds_no_impr_cnt_limit
  int no_impr_cnt_limit2 = 10;  

  if(fix_incr_cnt < 1) fix_incr_cnt = 1;
  else if(fix_incr_cnt > init_frac_ip_cnt) fix_incr_cnt = init_frac_ip_cnt;
  
  //dive_depth_limit = (int)(2.0*init_frac_ip_cnt/fix_incr_cnt) + 1;  
  
  //if(init_frac_ip_cnt > 200){
  // dive_depth_limit = MIN(500, (int)(2.0*init_frac_ip_cnt/


  //dive_depth_limit = 3*((int)(sqrt(init_frac_ip_cnt/fix_incr_cnt + 4)) + 1);

  frac_ip_cnt = init_frac_ip_cnt; 
  int_ip_cnt = init_int_ip_cnt;
  memcpy(frac_ind, init_frac_ind, ISIZE*frac_ip_cnt);  

  if(has_ub) set_obj_upper_lim(diving_lp, obj_ub - p->par.granularity + 
			       lp_data->lpetol);
  
  //int disc_code = 0;

  int min_iter = 1000000000, max_iter = -min_iter, tot_iter_cnt = 0, 
    avg_iter_cnt = 0; 
  int allowed_iter_num;
  int d_fixed_limit; //dual itlim fixed limit

  //int remaining_iter_limit;
  //lp_iter_limit = MAX(10*frac_ip_cnt, 10*(int)(1.0*1000000000/lp_data->nz));
  if(p->bc_index < 1 && is_last_iter) {
    lp_iter_limit *= 5; //100000;
    allowed_iter_num = 2000;
    //max_allowed_iter_num = 1000;
    d_fixed_limit = 20*fix_incr_cnt;
  }else{
    //lp_iter_limit = 10000;
    allowed_iter_num = 1000;
    d_fixed_limit = 5*fix_incr_cnt;
    dive_depth_limit = 100;
    no_impr_cnt_limit = 5; 
    no_prog_cnt_limit = 5; 
    no_better_cnt_limit = 5;
    no_impr_cnt_limit2 = 5; 
  }

  if(p->par.rs_mode_enabled){     
     allowed_iter_num = 200;     
     dive_depth_limit = 20;
     no_impr_cnt_limit = 3; 
     no_prog_cnt_limit = 3; 
     no_better_cnt_limit = 3;
     no_impr_cnt_limit2 = 3;
  }
  
  //lp_iter_limit = MAX(0, ((int)(2.5*1e9/lp_data->nz) - p->lp_stat.ds_num_iter));  

  //if(lp_iter_limit > 0) lp_iter_limit = (int)(1.0*lp_iter_limit/d_cnt) + 1; 
  //else return FALSE;

  double * obj = const_cast <double *> (lp_data->si->getObjCoefficients());
  //MAX(10*frac_ip_cnt, 10*(int)(1.0*1000000000/lp_data->nz));
  int single_iter_limit = -1; //(int)(1.0*lp_iter_limit/frac_ip_cnt); 

  double *x_rank = lp_data->tmp.d;
  char min_dir, * direction = lp_data->tmp.c;
  int min_ind, d_fixed_cnt;
  int fixed_cnt; 
  int better_expected_cnt = MAX(1, (int)(init_frac_ip_cnt*3.0/dive_depth_limit));	    
  char abandon_lp;
  char other_tried;
  //printf("n %i frac_ip_cnt %i nz %i\n", lp_data->n, frac_ip_cnt, lp_data->nz);
  //x = lp_data->x;
  int ignore_type = -1; 
  for(int k = 0; k < d_cnt; k++){

     if(k == ignore_type) continue; 

     d_type = diving_type[k];
     //printf("DIVING - %i\n", d_type);
     mark_time = wall_clock(NULL);
     tot_lp_iter = 0;
     dive_depth = 0;
     no_impr_cnt = 0;
     int no_impr_cnt2 = 0;
     no_prog_cnt = 0;
     no_better_cnt = 0;
     d_fixed_cnt = 0;
     single_iter_limit = allowed_iter_num;
     x = lp_data->x;
     can_iterate = TRUE;
     //disc_code = 0;
     min_frac_ip_cnt = frac_ip_cnt;
     iter_cnt = 0;
     
     if(has_ub) {
	double adj_obj_ub = obj_ub - p->par.granularity + lp_data->lpetol; 
	if(adj_obj_ub > obj_ub - 100*lp_data->lpetol) adj_obj_ub = obj_ub - 100*lp_data->lpetol; 
	//printf("obj_ub - adj_ub %f %f \n", obj_ub, adj_obj_ub);
	set_obj_upper_lim(diving_lp, adj_obj_ub);
     }
     while(true){
	
	//fix_incr_cnt += MIN(4, (int)(1.0*frac_ip_cnt/100));
	fix_incr_cnt = 1;
	
	if(frac_ip_cnt > 1000) fix_incr_cnt += 100;
	
	if(dive_depth > dive_depth_limit/2){	  
	   fix_incr_cnt =(int)(1.0*(frac_ip_cnt - d_fixed_cnt)/
			       (dive_depth_limit > (dive_depth ? dive_depth_limit - dive_depth : 1)));
	   if(fix_incr_cnt < 1) fix_incr_cnt = 1;
	}//else{
	// if(frac_ip_cnt > 200) fix_incr_cnt += (int)(frac_ip_cnt/200.0); 
	//}
	
	if(dive_depth > dive_depth_limit && lp_iter_limit - tot_lp_iter < 2*iter_cnt){
	   fix_incr_cnt = (int)((frac_ip_cnt - d_fixed_cnt)/(dive_depth/4.0));
	}
	
	//printf("%i %i %i %i %i %i %i\n", dive_depth, d_type, frac_ip_cnt, fix_incr_cnt, 
	//     single_iter_limit, d_fixed_cnt, tot_lp_iter);
	
	/*
	  switch(d_type){
	  case VLENGTH_FIX_DIVING:
	  case GUIDED_FIX_DIVING:
	  case CROSSOVER_FIX_DIVING:
	  case EUC_FIX_DIVING:
	  case RANK_FIX_DIVING:
	  case FRAC_FIX_DIVING:
	  if(p->bc_level < 1 && frac_ip_cnt > 1000) fix_incr_cnt = frac_ip_cnt - 1000; 
	  else if (p->bc_level > 0 && frac_ip_cnt > 500) fix_incr_cnt = frac_ip_cnt - 500; 
	  break;
	  case VLENGTH_DIVING:
	  case GUIDED_DIVING:
	  case CROSSOVER_DIVING:
	  case EUC_DIVING:
	  case RANK_DIVING:
	  case FRAC_DIVING:
	  if(p->bc_level < 1 && frac_ip_cnt > 500) fix_incr_cnt = frac_ip_cnt - 500; 
	  else if (p->bc_level > 0 && frac_ip_cnt > 100) fix_incr_cnt = frac_ip_cnt - 100; 
	  default:
	  break;
	  }
	*/
	if((fixed_cnt =
	    ds_fix_vars(p, diving_lp, x, frac_ind, frac_ip_cnt, 
			d_fixed_cnt, fix_incr_cnt, d_type, obj, 
			has_ub ? colSolution : NULL, x_rank, direction, &min_ind,
			&min_dir, TRUE)) < 0){
	   //disc_code = 1;
	   break;
	}
	
	if(frac_ip_cnt <= d_fixed_cnt + fixed_cnt){
	   single_iter_limit = MAX(5*allowed_iter_num, 
				   lp_iter_limit - tot_lp_iter);
	}
	
	if(frac_ip_cnt > 1){
	   single_iter_limit = 2*allowed_iter_num;
	}
	//if(d_fixed_cnt > 0 && d_fixed_cnt % 5 == 0){
	//  single_iter_limit = 1000;
	//}
	//if(d_fixed_cnt > 0) 
	// single_iter_limit += d_fixed_cnt * allowed_iter_num;
	//}else{
	//single_iter_limit = MAX(5*max_allowed_iter_num, 
	//			 lp_iter_limit - tot_lp_iter);
	// }
	
	set_itlim(diving_lp, single_iter_limit);      
	int termcode;
	abandon_lp = FALSE;
	other_tried = FALSE;
	while(true){
	   if((termcode = dual_simplex(diving_lp, &iter_cnt)) != LP_OPTIMAL){
	      //printf("termcode %i\n", termcode);
	      //disc_code  = 2;
	      if(termcode == LP_D_ITLIM){
		 d_fixed_cnt += fixed_cnt;
		 if(frac_ip_cnt <= d_fixed_cnt){
		    abandon_lp = TRUE;
		 }
		 single_iter_limit += allowed_iter_num;
		 //continue;
		 abandon_lp = TRUE;
	      }else{
		 if(fixed_cnt < 2 && !other_tried){
		    double bd = (min_dir == 'L' ? ceil(x[min_ind]) : floor(x[min_ind]));
		    diving_lp->si->setColLower(min_ind, bd);
		    diving_lp->si->setColUpper(min_ind, bd);
		    other_tried = TRUE;
		    tot_lp_iter += iter_cnt; 
		    continue;
		 }else{
		    abandon_lp = TRUE;
		 }
	      }
	      //	continue_iterate = FALSE;
	   }else{
	      single_iter_limit = iter_cnt + allowed_iter_num;
	      d_fixed_cnt = 0;
	   }
	   break;
	}    
	
	if(dive_depth <= 1){
	   //lp_iter_limit = ((int)(iter_cnt/5.0) + 5)*dive_depth_limit;
	   if(min_iter > iter_cnt) min_iter = iter_cnt;
	   if(max_iter < iter_cnt) max_iter = iter_cnt;
	   tot_iter_cnt += iter_cnt;
	   avg_iter_cnt++;
	   //break;
	}
	dive_depth++;
	
	if(abandon_lp) break;
	
	if(d_fixed_cnt < 1){
	   get_x(diving_lp);
	   x = diving_lp->x;
	   
	   prev_frac_ip_cnt = frac_ip_cnt; 
	   ds_get_frac_vars(lp_data, x, frac_ind, &frac_ip_cnt, &int_ip_cnt);
	   
	   if(frac_ip_cnt >=
	      init_frac_ip_cnt -
	      (int)(0.7*fix_incr_cnt*dive_depth)) no_impr_cnt++;
	   if(frac_ip_cnt < min_frac_ip_cnt){
	      min_frac_ip_cnt = frac_ip_cnt;
	      no_impr_cnt2 = 0;
	   }else{
	      no_impr_cnt2++;
	   }
	   
	   if(frac_ip_cnt > prev_frac_ip_cnt - fix_incr_cnt) no_prog_cnt++;	 
	   if(dive_depth % 3 == 0){
	      int d_times = (int)(1.0*dive_depth/3);
	      if(init_frac_ip_cnt - frac_ip_cnt < d_times * better_expected_cnt)
		 no_better_cnt++;
	   }
	   
	   //printf("min: %i no_impr: %i no_impr2:%i no_prog: %i no_better: %i\n\n",
	   //min_frac_ip_cnt, no_impr_cnt, no_impr_cnt2, no_prog_cnt, no_better_cnt);
	   
	   //else no_prog_cnt--;
	   
	   if(frac_ip_cnt <= 0){
	      //printf("found_feas nz cnt %i %i\n", lp_data->nz, p->lp_stat.ds_type_calls[d_type]);
	      is_ip_feasible = TRUE;
	      *solutionValue = diving_lp->objval; 
	      printf("DS-FEAS: %i - %i ----- %f %f\n", p->bc_index, d_type, has_ub? obj_ub : 0.0, *solutionValue);
	      //printf("sol: %f\n", *solutionValue);
	      memcpy(betterSolution, diving_lp->x, DSIZE*n);
	   }
	}
	
	//      if(is_ip_feasible || d_fixed_cnt > d_fixed_limit ||  
	//	 (((tot_lp_iter += iter_cnt) >= lp_iter_limit &&
	//	   dive_depth >= dive_depth_limit &&
	//	   no_impr_cnt >= no_impr_cnt_limit &&
	//	   no_prog_cnt >= no_prog_cnt_limit) &&
	//	  1.0*frac_ip_cnt/init_frac_ip_cnt > force_lp_ratio)){
	
	if(p->bc_level < 1){
	   if((//no_impr_cnt >= no_impr_cnt_limit &&
	       no_prog_cnt >= no_prog_cnt_limit &&
	       no_impr_cnt2 >= no_impr_cnt_limit2)&&  
	      no_better_cnt >= no_better_cnt_limit){// && 
	      //tot_lp_iter > lp_iter_limit/20){
	      can_iterate = FALSE;
	   }
	}else{
	   if(((//no_impr_cnt >= no_impr_cnt_limit &&
		//no_prog_cnt >= no_prog_cnt_limit &&
		no_impr_cnt2 >= no_impr_cnt_limit2) && // ||   
	       no_better_cnt >= no_better_cnt_limit)){// && tot_lp_iter > lp_iter_limit/10) {
	      can_iterate = FALSE;
	   }
	}
	
	//if(is_ip_feasible || d_fixed_cnt > d_fixed_limit ||  
	// (tot_lp_iter += iter_cnt) >= lp_iter_limit ||
	// (dive_depth >= dive_depth_limit &&
	//  no_impr_cnt >= no_impr_cnt_limit &&
	//  no_prog_cnt >= no_prog_cnt_limit)){ //&&
	//1.0*frac_ip_cnt/init_frac_ip_cnt > force_lp_ratio)){
	//if(is_ip_feasible){printf("sol: %f\n", *solutionValue);}
	double frac_ip_ratio = 1.0*frac_ip_cnt/n; 
	tot_lp_iter += iter_cnt; 
	if(is_ip_feasible || d_fixed_cnt > d_fixed_limit ||  !can_iterate ||
	   //((frac_ip_ratio < 0.2 || p->bc_index < 1) &&
	   (tot_lp_iter >= lp_iter_limit && 
	    dive_depth >= dive_depth_limit) || tot_lp_iter >= 2*lp_iter_limit){
	   //|| 
	   //((frac_ip_ratio >= 0.2 || p->bc_index >= 1) && 
	   //(tot_lp_iter > lp_iter_limit || dive_depth >= dive_depth_limit))){
#if 0
	  if(is_ip_feasible){printf("sol: %f\n", *solutionValue); disc_code = 3;}
	  if((tot_lp_iter += iter_cnt) >= lp_iter_limit && 
	     1.0*frac_ip_cnt/init_frac_ip_cnt > force_lp_ratio) disc_code = 4;
	  if(dive_depth >= dive_depth_limit && 
	     1.0*frac_ip_cnt/init_frac_ip_cnt > force_lp_ratio) disc_code = 5;
	  if(no_impr_cnt >= no_impr_cnt_limit &&
	     1.0*frac_ip_cnt/init_frac_ip_cnt > force_lp_ratio) disc_code = 6;
#endif
	  break; 
	}
	
	if(p->par.ds_solve_ip && 
	   frac_ip_cnt  <= ip_solve_col_ratio * ip_vars_cnt + 1 && 
	   (!has_ub || 
	    (has_ub && 100*(obj_ub - diving_lp->objval)/(fabs(obj_ub) + 0.0001) > 
	     ip_solve_gap))){
	   /* solve ip here */
	   if(!env){
	      env = lp_to_sym(p, diving_lp, FALSE, 0, NULL, NULL, NULL, NULL, NULL);
	      for(i = 0; i < n; i++){
		 if(lp_data->vars[i]->is_int){
		    sym_set_integer(env, i);
		 }
	      }
	      sym_set_int_param(env, "verbosity", -2);
	      sym_set_int_param(env, "ds_solve_ip", FALSE);
	      sym_set_int_param(env, "ds_enabled", FALSE);
	      sym_set_int_param(env, "fr_enabled", FALSE);
	      //sym_set_int_param(env, "find_first_feasible", TRUE);
	   }else{
	      double * lb = const_cast <double *> (diving_lp->si->getColLower());
	      double * ub = const_cast <double *> (diving_lp->si->getColUpper());
	      for(i = 0; i < n; i++){
		 sym_set_col_lower(env, i, lb[i]);
		 sym_set_col_upper(env, i, ub[i]);
	      }
	   }	  
	   if(has_ub) sym_set_dbl_param(env, "upper_bound", obj_ub - 
					p->par.granularity + lp_data->lpetol);	
	   sym_solve(env);
	   
	   if(sym_is_iteration_limit_reached(env) || 
	      sym_is_proven_optimal(env)){
	      //if(env->warm_start->stat.created > 1){
	      sym_get_col_solution(env, betterSolution);
	      sym_get_obj_val(env, solutionValue);  
	      is_ip_feasible = TRUE;
	   }
	   break;
	}
     }
     
     //if(d_type == 4){
     //  printf("iters lp:%i avg:%i max:%i min:%i set:%i n:%i ip_cnt:%i frac_cnt:%i frac_ratio:%.2f nz:%i\n",
     //      p->lp_stat.lp_total_iter_num, 
     //      tot_iter_cnt/avg_iter_cnt, max_iter, min_iter, lp_iter_limit/avg_iter_cnt,
     //      n, ip_vars_cnt, frac_ip_cnt, 1.0*frac_ip_cnt/ip_vars_cnt, lp_data->nz); 
     //  abort();
     //}
     //printf("DIVING DONE - %i\n", d_type);
     
     p->lp_stat.ds_type_num_iter[d_type] += tot_lp_iter;
     p->lp_stat.ds_num_iter += tot_lp_iter;
     if(dive_depth > 0){
	(p->lp_stat.ds_type_calls[d_type])++;
     }
     p->comp_times.ds_type[d_type] += wall_clock(NULL) - mark_time; 
     
     //if(disc_code != 1){
     // printf("d_code %i\n", disc_code);
     //}
     if(is_ip_feasible) {
	//if(!has_ub) has_ub = TRUE;
	//obj_ub = *solutionValue;
	(p->lp_stat.ds_type_num_sols[d_type])++;
	p->lp_stat.ds_num_sols++;
	//is_ip_feasible = FALSE;
	//is_feas = TRUE;
	is_return_feasible = TRUE; 
	if(p->bc_index <= 0){
	   if(d_gap(*solutionValue, t_lb, p->mip->obj_offset,
		    p->mip->obj_sense) > p->par.ds_min_gap){
	      ignore_type = k; 
	      k = -1;
	      memcpy(colSolution, betterSolution, DSIZE*n);
	      has_ub = TRUE;
	      obj_ub = *solutionValue;	     
	      is_ip_feasible = FALSE;
	   }else{
	      break;
	   }
	}else{
	   break;
	}
     }
     if(dive_depth > 0){
	frac_ip_cnt = init_frac_ip_cnt;
	int_ip_cnt = init_int_ip_cnt;
	memcpy(frac_ind, init_frac_ind, ISIZE*frac_ip_cnt);
	load_basis(diving_lp, cstat, rstat);
	
	for(i = 0; i < n; i++){
	   get_ub(lp_data, i, &ub);
	   get_lb(lp_data, i, &lb);
	   //change_lbub(diving_lp, frac_ind[i], lb[frac_ind[i]], ub[frac_ind[i]]);
	   change_lbub(diving_lp, i, lb, ub);
	}
	//for(i = 0; i < frac_ip_cnt; i++){
	//  get_ub(lp_data, frac_ind[i], &ub);
	//  get_lb(lp_data, frac_ind[i], &lb);
	  //change_lbub(diving_lp, frac_ind[i], lb[frac_ind[i]], ub[frac_ind[i]]);
	// change_lbub(diving_lp, frac_ind[i], lb, ub);
	//}
     }
  }
  
  diving_lp->x = 0;
  //FREE(is_fixed);
  //  FREE(ub);
  //  FREE(lb);
  //  FREE(frac_ind);
  // FREE(init_frac_ind);  
  //FREE(cstat);
  //FREE(rstat);
  
  close_lp_solver(diving_lp);
  FREE(diving_lp->lb);
  FREE(diving_lp->ub);
  FREE(diving_lp);

  if(env) sym_close_environment(env);

  p->lp_stat.ds_calls++;
  p->lp_stat.ds_last_call_ind = p->bc_index; 
  //printf("ds_call cnt %i\n", p->lp_stat.ds_calls);

  p->comp_times.ds += wall_clock(NULL) - start_time; 

  //abort(); 
  return is_return_feasible;
  //return is_ip_feasible;
  //return is_feas;
}

/*===========================================================================*/  
  /* fractional diving */
  /* vectorlength  diving */
  /* euc-search diving*/
  /* guided diving */
  /* rootsolution diving */
  /* coefficient diving */
  /* pseudocost diving */
  /* ranking diving */

/*===========================================================================*/  

int ds_fix_vars(lp_prob *p, LPdata *diving_lp, double *x, 
		int *frac_ind, int frac_cnt, int d_fixed_cnt, 
		int fix_incr_cnt, int d_type, double *obj, 
		double *ip_sol, double *x_rank, char *direction,
		int *min_ind, char *min_dir, char should_fix)
{   
   int i, j, n = diving_lp->n;
   int ind;
   //double *x_rank = p->lp_data->tmp.d;
   //char * direction = p->lp_data->tmp.c;
   int vars_eff_cnt = 0;
   double *root_lp = 0;
   double base_rank;
   double min_x_rank = SYM_INFINITY;
   char min_x_dir = 'L';
   int min_x_ind = 0;
   int fixed_cnt = fix_incr_cnt;
   double x_obj, val, etol = p->lp_data->lpetol;
   // double x_rank_1257;
   // double x_rank_1006;
   sp_desc *sp = p->tm->sp;
   sp_solution *sol;

   int *x_dir_cnt = p->lp_data->tmp.i2 + n;
   double *x_diff = p->lp_data->tmp.d + n;   
   double *feas_sol = p->lp_data->tmp.d + 2*n;

   double min_obj = DBL_MAX;

   if(frac_cnt < 1 || fix_incr_cnt < 1 || 
      d_fixed_cnt >= frac_cnt) {
      printf("problem ds_fix_vars\n");
      return -1;
   }

   if(fix_incr_cnt + d_fixed_cnt > frac_cnt) 
     fixed_cnt = frac_cnt - d_fixed_cnt;

   if(d_fixed_cnt < 1){
      switch(d_type){
      case FRAC_FIX_DIVING: 	 
	if(!p->par.ds_fractional_fix_enabled)// || !ip_sol) 
	  return -1;
	ds_fix_common_vars(diving_lp, p->lp_data->vars, ip_sol, x);
       case FRAC_DIVING:
	 if(d_type != FRAC_FIX_DIVING && !p->par.ds_fractional_enabled) return -1;	 
	 for(i = 0; i < frac_cnt; i++){
	    ind = frac_ind[i];
	    x_rank[i] = x[ind] - floor(x[ind]);
	    direction[ind] = 'L';
	    if(x_rank[i] > 0.5) {
	       x_rank[i] = 1 - x_rank[i];
	       direction[ind] = 'U';
	    }
	    //x_rank[i] = -x_rank[i];
	    if(x_rank[i] < min_x_rank){
	       min_x_rank = x_rank[i];
	       min_x_ind = ind;
	       min_x_dir = direction[ind];
	    }
	 }
	 break;
       case VLENGTH_FIX_DIVING: 	 
	 if(!p->par.ds_vlength_fix_enabled || !p->mip->mip_inf)// || !ip_sol) 
	    return -1;
	 ds_fix_common_vars(diving_lp, p->lp_data->vars, ip_sol, x);
       case VLENGTH_DIVING:
	 if(d_type != VLENGTH_FIX_DIVING && 
	    (!p->par.ds_vlength_enabled || !p->mip->mip_inf)) return -1;
	 vars_eff_cnt = 0;
	 //obj = const_cast <double *> (diving_lp->si->getObjCoefficients());
	 for(i = 0; i < frac_cnt; i++){
	    ind = frac_ind[i];
	    //vars_eff_cnt = MAX(p->mip->mip_inf->cols[ind].sos_num,
	    //	       p->mip->mip_inf->cols[ind].col_size);
	    vars_eff_cnt = p->mip->mip_inf->cols[ind].nz;
	    //vars_eff_cnt = p->mip->mip_inf->cols[ind].col_size;
	    if(obj[ind] < 0){
	       x_rank[i] = (obj[ind]*(ceil(x[ind]) - x[ind]))/vars_eff_cnt;
	       direction[ind] = 'U';
	    }
	    else if(obj[ind] > 0){
	       x_rank[i] = (-obj[ind]*(x[ind] - floor(x[ind])))/vars_eff_cnt;
	       direction[ind] = 'L';
	    }else{
	       x_rank[i] = ceil(x[ind]) - x[ind];
	       direction[ind] = 'U';
	       if(x_rank[i] > 0.5) {
		  x_rank[i] = 1 - x_rank[i];
		  direction[ind] = 'L';
	       }
	       x_rank[i] /= vars_eff_cnt;
	    }
#if 0			   
	    if(obj[ind] > 0){
	       x_rank[i] = (obj[ind]*(ceil(x[ind]) - x[ind]))/vars_eff_cnt;
	       direction[ind] = 'U';
	    }
	    else if(obj[ind] < 0){
	       x_rank[i] = (-obj[ind]*(x[ind] - floor(x[ind])))/vars_eff_cnt;
	       direction[ind] = 'L';
	    }else{
	       x_rank[i] = ceil(x[ind]) - x[ind];
	       direction[ind] = 'U';
	       if(x_rank[i] > 0.5) {
		  x_rank[i] = 1 - x_rank[i];
		  direction[ind] = 'L';
	       }
	       x_rank[i] /= vars_eff_cnt;
	    }
#endif
	    if(x_rank[i] < min_x_rank){
	       min_x_rank = x_rank[i];
	       min_x_ind = ind;
	       min_x_dir = direction[ind];
	    }
	    // if(ind == 1257) x_rank_1257=x_rank[i];
	    // if(ind == 1006) x_rank_1006=x_rank[i];
	 }
	 break;
       case EUC_FIX_DIVING: 	 
	 if(!p->par.ds_euc_fix_enabled || !p->root_lp || //!ip_sol || 
	    (p->bc_index < 0 && p->lp_stat.lp_calls < 2)) return -1;
	 ds_fix_common_vars(diving_lp, p->lp_data->vars, ip_sol, x);
       case EUC_DIVING:
	 if(d_type != EUC_FIX_DIVING && 
	    (!p->par.ds_euc_enabled || !p->root_lp ||
	     (p->bc_index < 0 && p->lp_stat.lp_calls < 2))) return -1;
	 root_lp = p->root_lp;
	 for(i = 0; i < frac_cnt; i++){
	    ind = frac_ind[i];
	    base_rank = pow(x[ind] - root_lp[ind], 2);
	    x_rank[i] = base_rank +
	       pow(ceil(x[ind]) - root_lp[ind], 2);
	    base_rank +=
	       pow(floor(x[ind]) - root_lp[ind], 2);
	    direction[ind] = 'U';
	    if(base_rank < x_rank[i]){
	       x_rank[i] = base_rank;
	       direction[ind] = 'L';
	    }
	    if(x_rank[i] < min_x_rank){
	       min_x_rank = x_rank[i];
	       min_x_ind = ind;
	       min_x_dir = direction[ind];
	    }
	 }
	 break;
       case GUIDED_FIX_DIVING:
	 if(!p->par.ds_guided_fix_enabled || !ip_sol) return -1;
	 ds_fix_common_vars(diving_lp, p->lp_data->vars, ip_sol, x);
       case GUIDED_DIVING:
	 if(d_type != GUIDED_FIX_DIVING && 
	    (!p->par.ds_guided_enabled || !ip_sol)) return -1;
	 for(i = 0; i < frac_cnt; i++){
	    ind = frac_ind[i];
	    x_rank[i] = fabs(x[ind] - ip_sol[ind]);
	    //x_fix[ind] = ip_sol[ind]; -FIXME
	    direction[ind] = 'L';
	    if(ip_sol[ind] > x[ind] + etol) direction[ind] = 'U';
	    if(x_rank[i] < min_x_rank){
	       min_x_rank = x_rank[i];
	       min_x_ind = ind;
	       min_x_dir = direction[ind];
	    }
	 }
	 break;
      case CROSSOVER_FIX_DIVING: 	 
	if(!p->par.ds_crossover_fix_enabled || !sp || sp->num_solutions < 2 || 
	   !p->mip->mip_inf) return -1;
	ds_fix_common_vars(diving_lp, p->lp_data->vars, ip_sol, x);
       case CROSSOVER_DIVING: 
	 if(d_type != CROSSOVER_FIX_DIVING && 
	    (!p->par.ds_crossover_enabled || !sp || sp->num_solutions < 2 || 
	     !p->mip->mip_inf)) return -1;
	 for(i = 0; i < n; i++){
	   x_dir_cnt[i] = 0;
	   x_diff[i] = DBL_MAX;
	   if(obj[i] < min_obj) min_obj = obj[i];
	 }
	 min_obj = fabs(min_obj);
	 
	 for(i = 0; i < sp->num_solutions; i++){
	   sol = sp->solutions[i];
	   memset(feas_sol, 0, DSIZE*n);
	   for(j = 0; j < sol->xlength; j++){
	     feas_sol[sol->xind[j]] = sol->xval[j];
	   }
	   
	   for(j = 0; j < frac_cnt; j++){
	     ind = frac_ind[j];
	     val = feas_sol[ind];
	     if(val >= x[ind] + etol) x_dir_cnt[ind]++;
	     else if (val <= x[ind] - etol) x_dir_cnt[ind]--;
	     if(fabs(val - x[ind]) < x_diff[ind]){
	       x_diff[ind] = val - x[ind];
	     }
	   }
	 }
	 
	 for(i = 0; i < frac_cnt; i++){
	   ind = frac_ind[i];
	   vars_eff_cnt = p->mip->mip_inf->cols[ind].nz;
	   x_obj = obj[ind] + min_obj + 1e-4;
	   direction[ind] = 'U';
	   if(x_dir_cnt[ind] < 0 ||
	      (x_dir_cnt[ind] == 0 && x_diff[ind] < 0.0)){
	       direction[ind] = 'L';
	   }
	   
	   if(direction[ind] == 'U'){
	     x_rank[i] = x_obj*
	       (ceil(x[ind]) - x[ind])/vars_eff_cnt;
	   }else{
	     x_rank[i] = x_obj*
	       (x[ind] - floor(x[ind]))/vars_eff_cnt;
	   }
	   if(x_rank[i] < min_x_rank){
	     min_x_rank = x_rank[i];
	     min_x_ind = ind;
	     min_x_dir = direction[ind];
	   }
	 }
	 break;
       case ROOT_DIVING:
	 return -1;
	 if(!p->par.ds_root_enabled)
	    break;
       case COEFF_DIVING:
	 return -1;
	 if(!p->par.ds_coeff_enabled) return -1;
	 break;
       case PC_DIVING:
	 return -1;
	 if(!p->par.ds_pc_enabled || !p->root_lp) return -1;
#if 0
	 double *pcost_down = p->pcost_down;
	 double *pcost_up = p->pcost_up;
	 for(i = 0; i < frac_cnt; i++){
	    ind = frac_ind[i];
	    x_rank[i] = pcost_down[ind] * (x[ind] - floor(x[ind]));
	    direction[ind] = 'L';
	    up_chg = pcost_up[ind] * (ceil(x[ind]) - x[ind]);
	    if(up_chg < x_rank[i]) {
	       x_rank[i] = up_chg;
	       direction[ind] = 'U';
	    }
	 }
#endif
	 break;
       case RANK_FIX_DIVING: 	 
	 if(!p->par.ds_rank_fix_enabled || !ip_sol || 
	    (p->bc_index < 0 && p->lp_stat.lp_calls < 2) || p->var_rank_cnt < 1) return -1;
	 ds_fix_common_vars(diving_lp, p->lp_data->vars, ip_sol, x);
       case RANK_DIVING:
	 if(d_type != RANK_FIX_DIVING && 
	    (!p->par.ds_rank_enabled ||
	     (p->bc_index < 0 && p->lp_stat.lp_calls < 2) || p->var_rank_cnt < 1)) return -1;
	 for(i = 0; i < frac_cnt; i++){
	    ind = frac_ind[i];
	    x_rank[i] = p->var_rank[ind]/p->var_rank_cnt;//p->lp_stat.lp_node_calls;
	    //(x[ind] - floor(x[ind])))/(p->lp_stat.lp_node_calls + 1);
	    direction[ind] = 'L';
	    if(x_rank[i] > 0.5) direction[ind] = 'U';
	    if(x_rank[i] < min_x_rank){
	       min_x_rank = x_rank[i];
	       min_x_ind = ind;
	       min_x_dir = direction[ind];
	    }
	 }
	 break;
       default:
	 return -1;
      }
   }

   if(should_fix){
      *min_dir = min_x_dir;
      *min_ind = min_x_ind;
      if(d_fixed_cnt < 1 && fix_incr_cnt < 2){
	 //frac_ind[0] = min_x_ind;
	 //direction[min_x_ind] = min_x_dir;
	double bd = (min_x_dir == 'L' ? 
		     floor(x[min_x_ind]) : ceil(x[min_x_ind]));
	 diving_lp->si->setColLower(min_x_ind, bd);
	 diving_lp->si->setColUpper(min_x_ind, bd);
	 //    if(d_type == VLENGTH_DIVING && p->lp_data->nz == 10032){
	 // printf("\nmin_x_ind %i %c %f %f %f\n",
	 //	min_x_ind, min_x_dir, bd, x[min_x_ind], diving_lp->objval);
	 // printf("x[1257]:%f x[1006]:%f x_rank[1257]:%f x_rank[1006]:%f\n",
	 //	x[1257], x[1006], x_rank_1257, x_rank_1006);
	 //}
	 return 1;
      }else if(d_fixed_cnt <= fix_incr_cnt){
	 qsort_di(x_rank, frac_ind, frac_cnt);
      } //else x_rank should already be sorted when d_fixed_cnt <= fix_incr_cnt
      
      double bd;
      
      for(i = d_fixed_cnt; i < d_fixed_cnt + fixed_cnt; i++){
	 ind = frac_ind[i];
	 bd = (direction[ind] == 'L' ? floor(x[ind]) : ceil(x[ind]));
	 diving_lp->si->setColLower(ind, bd);
	 diving_lp->si->setColUpper(ind, bd);
      }
      return fixed_cnt;
   }
   
   return 0;

}
		     
/*===========================================================================*/

int ds_get_frac_vars(LPdata *lp_data, double *x, int *indices, 
		     int *frac_ip_cnt, int *int_ip_cnt)
{

  int i, n = lp_data->n;
  double floorx, ceilx; 
  double etol = lp_data->lpetol;
  *frac_ip_cnt = *int_ip_cnt = 0;

  for(i = 0; i < n; i++){
    floorx = floor(x[i]); 
    ceilx = ceil(x[i]);
    if (lp_data->vars[i]->is_int){ 
       if(x[i] > floorx + etol && x[i] < ceilx - etol){      
	  indices[*frac_ip_cnt] = i;
	  (*frac_ip_cnt)++;
       }else (*int_ip_cnt)++;
    }
  }
  
  return 0;//(*frac_ip_cnt <= 0 ? TRUE : FALSE);

}



/*===========================================================================*/
// menal - adapted from cbc
// See if rounding will give solution



int round_solution(lp_prob *p, double *solutionValue, double *betterSolution, double t_lb)
{

   int numberColumns = p->mip->n;
   int numberRows = p->mip->m; 
   int nz = p->mip->nz;
   int returnCode = 0, numberIntegers = 0;
   double primalTolerance = p->lp_data->lpetol,
      integerTolerance = primalTolerance;
   double *lower, *upper, *solution, *objective;
   double direction = p->mip->obj_sense == SYM_MINIMIZE ? 1: -1 ;
   double newSolutionValue = direction*p->lp_data->objval;
   double *element, *elementByRow;
   int * integerVariable, row_ind, elem_ind;
   int *row, *column, *columnStart, *rowStart, *columnLength, *rowLength;
   int i, j;
   double total_time = 0;
   total_time = used_time(&total_time);

   if(!(p->mip->matbeg)){
     return returnCode;
   }
   
   get_bounds(p->lp_data);
   get_x(p->lp_data);
   
   lower = p->lp_data->lb;
   upper = p->lp_data->ub;
   solution = p->lp_data->x;

   element = p->mip->matval;
   row = p->mip->matind;
   columnStart = p->mip->matbeg;
   objective = p->mip->obj;

   columnLength = p->mip->col_lengths;

   if(!columnLength){
      columnLength=(p->mip->col_lengths = (int *)calloc(numberColumns,ISIZE));
      elementByRow = (p->mip->row_matval = (double *)malloc(nz*DSIZE)); 
      column = (p->mip->row_matind = (int *)malloc(nz*ISIZE)); 
      rowStart = (p->mip->row_matbeg = (int *)malloc((numberRows+1)*ISIZE));
      rowLength = (p->mip->row_lengths = (int *)calloc(numberRows,ISIZE));

      /* first get row legths */   
      for(i = 0; i < numberColumns; i++){
	 /* get orig indices here */
	 for(j = columnStart[i]; j < columnStart[i+1]; j++){
	    rowLength[row[j]]++;
	 }
	 columnLength[i] = columnStart[i+1] - columnStart[i];
      }
      
      rowStart[0] = 0;
      
      /* fill in matbegs */
      for(i = 0; i < numberRows; i++){
	 rowStart[i + 1] = rowStart[i] + rowLength[i];
      }

      /* get matrix, change 'G' rows to 'L'*/
      for(i = 0; i < numberColumns; i++){
	 for(j = columnStart[i]; j < columnStart[i+1]; j++){
	    row_ind = row[j];
	    elem_ind = rowStart[row_ind];
	    column[elem_ind] = i;

	    elementByRow[elem_ind] = element[j];
	    rowStart[row_ind] = elem_ind + 1;
	 }
      }

      for(i = 0; i < numberRows; i++){
	 rowStart[i] -= rowLength[i];
      }      
   }else{   
   
      elementByRow = p->mip->row_matval;
      column = p->mip->row_matind;
      rowStart = p->mip->row_matbeg;
      rowLength = p->mip->row_lengths;
   }


   const double * rowUpper = p->lp_data->si->getRowUpper();
   const double * rowLower = p->lp_data->si->getRowLower();

   integerVariable = p->lp_data->tmp.i1; //new int[numberColumns];
   
   for (i = 0; i<numberColumns; i++){
      if (p->mip->is_int[i]){
	 integerVariable[numberIntegers++] = i;
      }
   }
   
   // Get solution array for heuristic solution
   
   double * newSolution = p->lp_data->tmp.d;//new double [numberColumns];
   memcpy(newSolution,solution,numberColumns*sizeof(double));
   
   double * rowActivity = p->lp_data->tmp.d + numberColumns;//new double[numberRows];
   memset(rowActivity,0,numberRows*sizeof(double));
   for (i=0;i<numberColumns;i++) {
      int j;
      double value = newSolution[i];
      if (value) {
	 for (j=columnStart[i];
	      j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    //	printf("rowind %i: %i \n", j, iRow);
	    //	if(j < 5){
	    //	printf("element %i: %f \n", j, element[j]);
	    //	}
	    rowActivity[iRow] += value*element[j];
	 }
      }
   }
   // check was feasible - if not adjust (cleaning may move)
   for (i=0;i<numberRows;i++) {
      if(rowActivity[i]<rowLower[i]) {
	 //assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
	 rowActivity[i]=rowLower[i];
      } else if(rowActivity[i]>rowUpper[i]) {
	 //assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
	 rowActivity[i]=rowUpper[i];
      }
   }
   for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      double value=newSolution[iColumn];
      if (fabs(floor(value+0.5)-value)>integerTolerance) {
	 double below = floor(value);
	 double newValue=newSolution[iColumn];
	 double cost = direction * objective[iColumn];
	 double move;
	 if (cost>0.0) {
	    // try up
	    move = 1.0 -(value-below);
	 } else if (cost<0.0) {
	    // try down
	    move = below-value;
	 } else {
	    // won't be able to move unless we can grab another variable
	    // just for now go down
	    move = below-value;
	 }
	 newValue += move;
	 newSolution[iColumn] = newValue;
	 newSolutionValue += move*cost;
	 int j;
	 for (j=columnStart[iColumn];
	      j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int iRow = row[j];
	    rowActivity[iRow] += move*element[j];
	 }
      }
   }
   
   double penalty=0.0;
   
   // see if feasible
   for (i=0;i<numberRows;i++) {
      double value = rowActivity[i];
      double thisInfeasibility=0.0;
      if (value<rowLower[i]-primalTolerance)
	 thisInfeasibility = value-rowLower[i];
      else if (value>rowUpper[i]+primalTolerance)
	 thisInfeasibility = value-rowUpper[i];
      if (thisInfeasibility) {
	 // See if there are any slacks I can use to fix up
	 // maybe put in coding for multiple slacks?
	 double bestCost = 1.0e50;
	 int k;
	 int iBest=-1;
	 double addCost=0.0;
	 double newValue=0.0;
	 double changeRowActivity=0.0;
	 double absInfeasibility = fabs(thisInfeasibility);
	 for (k=rowStart[i];k<rowStart[i]+rowLength[i];k++) {
	    int iColumn = column[k];
	    if (columnLength[iColumn]==1) {
	       double currentValue = newSolution[iColumn];
	       double elementValue = elementByRow[k];
	       double lowerValue = lower[iColumn];
	       double upperValue = upper[iColumn];
	       double gap = rowUpper[i]-rowLower[i];
	       double absElement=fabs(elementValue);
	       if (thisInfeasibility*elementValue>0.0) {
		  // we want to reduce
		  if ((currentValue-lowerValue)*absElement>=absInfeasibility) {
		     // possible - check if integer
		     double distance = absInfeasibility/absElement;
		     double thisCost = -direction*objective[iColumn]*distance;
		     if (p->mip->is_int[iColumn]) {
			distance = ceil(distance-primalTolerance);
			if (currentValue-distance>=lowerValue-primalTolerance) {
			   if (absInfeasibility-distance*absElement< -gap-primalTolerance)
			      thisCost=1.0e100; // no good
			   else
			      thisCost = -direction*objective[iColumn]*distance;
			} else {
			   thisCost=1.0e100; // no good
			}
		     }
		     if (thisCost<bestCost) {
			bestCost=thisCost;
			iBest=iColumn;
			addCost = thisCost;
			newValue = currentValue-distance;
			changeRowActivity = -distance*elementValue;
		     }
		  }
	       } else {
		  // we want to increase
		  if ((upperValue-currentValue)*absElement>=absInfeasibility) {
		     // possible - check if integer
		     double distance = absInfeasibility/absElement;
		     double thisCost = direction*objective[iColumn]*distance;
		     if (p->mip->is_int[iColumn]) {
			distance = ceil(distance-primalTolerance);
		//assert (currentValue-distance<=upperValue+primalTolerance);
			if (absInfeasibility-distance*absElement< -gap-primalTolerance)
			   thisCost=1.0e100; // no good
			else
			   thisCost = direction*objective[iColumn]*distance;
		     }
		     if (thisCost<bestCost) {
			bestCost=thisCost;
			iBest=iColumn;
			addCost = thisCost;
			newValue = currentValue+distance;
			changeRowActivity = distance*elementValue;
		     }
		  }
	       }
	    }
	 }
	 if (iBest>=0) {
	    /*printf("Infeasibility of %g on row %d cost %g\n",
	      thisInfeasibility,i,addCost);*/
	    newSolution[iBest]=newValue;
	    thisInfeasibility=0.0;
	    newSolutionValue += addCost;
	    rowActivity[i] += changeRowActivity;
	 }
	 penalty += fabs(thisInfeasibility);
      }
   }
   
   // Could also set SOS (using random) and repeat
   if (!penalty) {
      // See if we can do better
      //seed_++;
      //CoinSeedRandom(seed_);
      // Random number between 0 and 1.
      double randomNumber = CoinDrand48();
      int iPass;
      int start[2];
      int end[2];
      int iRandom = (int) (randomNumber*((double) numberIntegers));
      start[0]=iRandom;
      end[0]=numberIntegers;
      start[1]=0;
      end[1]=iRandom;
      for (iPass=0;iPass<2;iPass++) {
	 int i;
	 for (i=start[iPass];i<end[iPass];i++) {
	    int iColumn = integerVariable[i];
	    //double value=newSolution[iColumn];
	    //assert (fabs(floor(value+0.5)-value)<integerTolerance);
	    double cost = direction * objective[iColumn];
	    double move=0.0;
	    if (cost>0.0)
	       move = -1.0;
	    else if (cost<0.0)
	       move=1.0;
	    while (move) {
	       bool good=true;
	       double newValue=newSolution[iColumn]+move;
	       if (newValue<lower[iColumn]-primalTolerance||
		   newValue>upper[iColumn]+primalTolerance) {
		  move=0.0;
	       } else {
		  // see if we can move
		  int j;
		  for (j=columnStart[iColumn];
		       j<columnStart[iColumn]+columnLength[iColumn];j++) {
		     int iRow = row[j];
		     double newActivity = rowActivity[iRow] + move*element[j];
		     if (newActivity<rowLower[iRow]-primalTolerance||
			 newActivity>rowUpper[iRow]+primalTolerance) {
			good=false;
			break;
		     }
		  }
		  if (good) {
		     newSolution[iColumn] = newValue;
		     newSolutionValue += move*cost;
		     int j;
		     for (j=columnStart[iColumn];
			  j<columnStart[iColumn]+columnLength[iColumn];j++) {
			int iRow = row[j];
			rowActivity[iRow] += move*element[j];
		     }
		  } else {
		     move=0.0;
		  }
	       }
	    }
	 }
      }
      if (newSolutionValue < *solutionValue) {
	 // paranoid check
	 memset(rowActivity,0,numberRows*sizeof(double));
	 for (i=0;i<numberColumns;i++) {
	    int j;
	    double value = newSolution[i];
	    if (value) {
	       for (j=columnStart[i];
		    j<columnStart[i]+columnLength[i];j++) {
		  int iRow=row[j];
		  rowActivity[iRow] += value*element[j];
	       }
	    }
	 }
	 // check was approximately feasible
	 bool feasible=true;
	 for (i=0;i<numberRows;i++) {
	    if(rowActivity[i]<rowLower[i]) {
	       if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
		  feasible = false;
	    } else if(rowActivity[i]>rowUpper[i]) {
	       if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
		  feasible = false;
	    }
	 }
	 if (feasible) {
	    // new solution
	    memcpy(betterSolution, newSolution, numberColumns*DSIZE);
	    *solutionValue = newSolutionValue;
	    //printf("** Solution of %g found by rounding\n",newSolutionValue);
	    returnCode=1;
	    p->lp_stat.rh_num_sols++;
	 } else {
	    // Can easily happen
	    //printf("Debug CbcRounding giving bad solution\n");
	 }
      }
   }
   //delete [] integerVariable;
 
   //delete [] newSolution;
   //delete [] rowActivity;

  p->comp_times.rh += used_time(&total_time);
  p->lp_stat.rh_calls++;
  p->lp_stat.rh_last_call_ind = p->bc_index; 
  return returnCode;
}

/*===========================================================================*/
/* --menal
  -adapted from cbc
  -added variable fixing
  -fixed errors and bugs
*/
/*===========================================================================*/
int local_search(lp_prob *p, double *solutionValue, double *colSolution,
		 double *betterSolution)
{
 
   LPdata *lp_data = p->lp_data;
   int numberColumns = p->mip->n;
   int numberRows = p->mip->m, nz = p->mip->nz;
   int returnCode = 0, numberIntegers = 0;
   double primalTolerance = lp_data->lpetol;
   double *solution = colSolution, *objective;
   double direction = p->mip->obj_sense == SYM_MINIMIZE ? 1: -1 ;
   double newSolutionValue = direction*(*solutionValue);
   double *element, *elementByRow;
   int * integerVariable, *rowStart;
   int *row, *columnStart, *columnLength, *column, *rowLength;
   int i, j, row_ind, elem_ind, n = numberColumns;

   double total_time = 0;
   total_time = used_time(&total_time);

  element = p->mip->matval;
  row = p->mip->matind;
  columnStart = p->mip->matbeg;
  objective = p->mip->obj;

  columnLength = p->mip->col_lengths;
  
  if(!columnLength){
     columnLength=(p->mip->col_lengths = (int *)calloc(numberColumns,ISIZE));
     elementByRow = (p->mip->row_matval = (double *)malloc(nz*DSIZE)); 
     column = (p->mip->row_matind = (int *)malloc(nz*ISIZE)); 
     rowStart = (p->mip->row_matbeg = (int *)malloc((numberRows+1)*ISIZE));
     rowLength = (p->mip->row_lengths = (int *)calloc(numberRows,ISIZE));
     
     /* first get row legths */   
     for(i = 0; i < numberColumns; i++){
	/* get orig indices here */
	for(j = columnStart[i]; j < columnStart[i+1]; j++){
	   rowLength[row[j]]++;
	}
	columnLength[i] = columnStart[i+1] - columnStart[i];
     }
     
     rowStart[0] = 0;
     
     /* fill in matbegs */
     for(i = 0; i < numberRows; i++){
	rowStart[i + 1] = rowStart[i] + rowLength[i];
     }
     
     /* get matrix, change 'G' rows to 'L'*/
     for(i = 0; i < numberColumns; i++){
	for(j = columnStart[i]; j < columnStart[i+1]; j++){
	   row_ind = row[j];
	   elem_ind = rowStart[row_ind];
	   column[elem_ind] = i;
	   
	   elementByRow[elem_ind] = element[j];
	   rowStart[row_ind] = elem_ind + 1;
	}
     }
     
     for(i = 0; i < numberRows; i++){
	rowStart[i] -= rowLength[i];
     }      
  }else{   
     
     elementByRow = p->mip->row_matval;
     column = p->mip->row_matind;
     rowStart = p->mip->row_matbeg;
     rowLength = p->mip->row_lengths;
  }
  

  const double * rowUpper = p->lp_data->si->getRowUpper();
  const double * rowLower = p->lp_data->si->getRowLower();

  double *lb, *ub; 
  double lpetol = lp_data->lpetol;

  /* fix some vars here */
  
  if(p->par.ls_fix_ratio > 0.0 && p->mip->mip_inf){

    ub = lp_data->tmp.d; 
    lb = lp_data->tmp.d + n;
    //memcpy(lb, const_cast <double *> (lp_data->si->getColLower()), DSIZE*n);
    //memcpy(ub, const_cast <double *> (lp_data->si->getColUpper()), DSIZE*n);

    memcpy(lb, p->mip->lb, DSIZE*n);
    memcpy(ub, p->mip->ub, DSIZE*n);

    double *x = lp_data->x;
    double *x_rank = lp_data->tmp.d + 2*n;
    double *x_rank2 = lp_data->tmp.d + 3*n;
    int ind, *x_ind = lp_data->tmp.i1;
    int *x_ind2 = lp_data->tmp.i1;
    int vars_eff_cnt, int_cnt = 0, int_cnt2 = 0;
    double bd, x_obj, min_obj = DBL_MAX, big_number = 1e20;
    double * obj = const_cast <double *> (lp_data->si->getObjCoefficients());
    for(i = 0; i < n; i++){
      if(obj[i] < min_obj) min_obj = obj[i];
    }

    min_obj = fabs(min_obj);

    for(i = 0; i < n; i++){
      if(lp_data->vars[i]->is_int && lp_data->ub[i] > lp_data->lb[i] + lpetol){
	if(obj[i] >= 0.0 && x[i] < lp_data->lb[i] + lpetol){
	  x_obj = obj[i] + min_obj + 1e-4;
	  vars_eff_cnt = MAX(p->mip->mip_inf->cols[i].sos_num,  p->mip->mip_inf->cols[i].col_size) + 1;
	  x_rank[int_cnt] = big_number*x_obj/vars_eff_cnt;
	  if(colSolution && x[i] < colSolution[i] + lpetol && x[i] > colSolution[i] - lpetol){
	    x_rank[int_cnt] /= 2;
	  }
          x_ind[int_cnt] = i;
	  int_cnt++;
	}else if(obj[i] <= 0.0 && x[i] > lp_data->ub[i] - lpetol){
	  x_obj = obj[i] + min_obj + 1e-4;
	  vars_eff_cnt = MAX(p->mip->mip_inf->cols[i].sos_num,  p->mip->mip_inf->cols[i].col_size) + 1;
	  x_rank2[int_cnt2] = big_number*x_obj/vars_eff_cnt;
	  if(colSolution && x[i] < colSolution[i] + lpetol && x[i] > colSolution[i] - lpetol){
	    x_rank2[int_cnt2] /= 2;
	  }
	  x_ind2[int_cnt2] = i;
	  int_cnt2++;
	}
      }
    }

    qsort_di(x_rank, x_ind, int_cnt);
    qsort_di(x_rank2, x_ind2, int_cnt2);

    int fix_cnt = MIN((int)(0.5*int_cnt), (int)(p->par.ls_fix_ratio*int_cnt));
    int fix_cnt2 = MIN((int)(0.5*int_cnt2), (int)(p->par.ls_fix_ratio*int_cnt2));

    //printf("F-cnt : %i\n", fix_cnt);
    for(i = 0; i < fix_cnt; i++) {
      ind = x_ind[i];
      bd = floor(x[ind] + lpetol);
      lb[ind] = ub[ind] = bd;
      //change_lbub(new_lp_data, ind, bd ,bd);
    }
    for(i = 0; i < fix_cnt2; i++) {
      ind = x_ind2[i];
      bd = floor(x[ind] + lpetol);
      lb[ind] = ub[ind] = bd;
      //change_lbub(new_lp_data, ind, bd ,bd);
    }
  }else{

    //lb =const_cast <double *> (lp_data->si->getColLower());
    //ub = const_cast <double *> (lp_data->si->getColUpper());
    ub = p->mip->ub;
    lb = p->mip->lb;
  }

  integerVariable = lp_data->tmp.i1;//new int[numberColumns];
  char *way = lp_data->tmp.c;
  char *mark = lp_data->tmp.c + numberColumns;

  for (i = 0; i<numberColumns; i++){
    if (lp_data->vars[i]->is_int && ub[i] > lb[i] + 100*lpetol){//p->mip->is_int[i]){
      integerVariable[numberIntegers++] = i;
    }
  }
   
  // Column copy
  /* 
  const double * element = matrix.getElements();
  const int * row = matrix.getIndices();
  const CoinBigIndex * columnStart = matrix.getVectorStarts();
  const int * columnLength = matrix.getVectorLengths();
  */

  // Get solution array for heuristic solution
  double * newSolution = lp_data->tmp.d; //new double [numberColumns];
  memcpy(newSolution,solution,numberColumns*sizeof(double));

  // way is 1 if down possible, 2 if up possible, 3 if both possible
  //char * way = new char[numberIntegers];
  // corrected costs
  double * cost = lp_data->tmp.d + numberColumns; //new double[numberIntegers];
  // for array to mark infeasible rows after iColumn branch
  //char * mark = new char[numberRows];
  memset(mark,0,numberRows);
  // space to save values so we don't introduce rounding errors
  double * save = lp_data->tmp.d + (numberColumns + numberIntegers);//new double[numberRows];

  // clean solution
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    
    // get original bounds
    //    double originalLower = lp_data->vars[iColumn]->lb; //p->mip->lb[iColumn];
    // double originalUpper = lp_data->vars[iColumn]->ub; //p->mip->ub[iColumn];

    double originalLower = lb[iColumn];//p->mip->lb[iColumn];
    double originalUpper = ub[iColumn];//p->mip->ub[iColumn];

    double value=newSolution[iColumn];

    if (value<originalLower) {
       value=originalLower;
       newSolution[iColumn]=value;
    } else if (value>originalUpper) {
       value=originalUpper;
       newSolution[iColumn]=value;
    }

    double nearest=floor(value+0.5);
    //assert(fabs(value-nearest)<10.0*primalTolerance);
    value=nearest;
    newSolution[iColumn]=nearest;
    // if away from lower bound mark that fact
    if (nearest>originalLower) {
      //      used_[iColumn]=1;
    }
    cost[i] = direction*objective[iColumn];
    int iway=0;
    
    if (value>originalLower+0.5) 
      iway = 1;
    if (value<originalUpper-0.5)
       iway |= 2;
    way[i]=(char)iway;
  }
  // get row activities
  double * rowActivity = lp_data->tmp.d + (numberColumns + numberIntegers + numberRows);//new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));

  for (i=0;i<numberColumns;i++) {
    int j;
    double value = newSolution[i];
    if (value) {
      for (j=columnStart[i];
	   j<columnStart[i]+columnLength[i];j++) {
	int iRow=row[j];
	rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  // if very infeasible then give up
  bool tryHeuristic=true;
  for (i=0;i<numberRows;i++) {
     if(rowActivity[i]<rowLower[i]) {
      if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowLower[i];
    } else if(rowActivity[i]>rowUpper[i]) {
      if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowUpper[i];
    }
  }
  //printf("In LS\n");
  if (tryHeuristic) {
     //printf("Trying LS %i %f\n", p->bc_index, *solutionValue);    
    // best change in objective
    double bestAllChange = 0.0;
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      double bestChange=0.0;      
      double objectiveCoefficient = cost[i];
      int k;
      int j;
      int goodK=-1;
      int wayK=-1,wayI=-1;
      if ((way[i]&1)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] -= element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try down
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (-objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=-1;
		bestChange = -objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (-objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=-1;
		bestChange = -objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if ((way[i]&2)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] += element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try up
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=1;
		bestChange = objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=1;
		bestChange = objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if (goodK>=0) {
	// we found something - update solution
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayI * element[j];
	}
	newSolution[iColumn] += wayI;
	int kColumn = integerVariable[goodK];
	for (j=columnStart[kColumn];
	     j<columnStart[kColumn]+columnLength[kColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayK * element[j];
	}
	newSolution[kColumn] += wayK;
	// See if k can go further ?
	// get original bounds
	double originalLower = lb[kColumn];//p->mip->lb[kColumn];
	double originalUpper = ub[kColumn];//p->mip->ub[kColumn];

	double value=newSolution[kColumn];
	int iway=0;
	if (value>originalLower+0.5) 
	  iway = 1;
	if (value<originalUpper-0.5) 
	  iway |= 2;
	way[goodK]=(char)iway;

	bestAllChange += bestChange; 
      }
    }
    if (bestAllChange+newSolutionValue<*solutionValue) {
       // paranoid check
      memset(rowActivity,0,numberRows*sizeof(double));
      double new_opt = 0.0;
      for (i=0;i<numberColumns;i++) {
	int j;
	double value = newSolution[i];
	new_opt += value*objective[i];
	if (value) {
	  for (j=columnStart[i];
	       j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    rowActivity[iRow] += value*element[j];
	  }
	}
      }
      int numberBad=0;
      double sumBad=0.0;
      
      // check was approximately feasible
      for (i=0;i<numberRows;i++) {
	 if(rowActivity[i]<rowLower[i]) {
	    sumBad += rowLower[i]-rowActivity[i];
	    if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	       numberBad++;
	 } else if(rowActivity[i]>rowUpper[i]) {
	    sumBad += rowUpper[i]-rowActivity[i];
	    if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
	       numberBad++;
	 }
      }
      if (!numberBad) {
	/*
	 for (i=0;i<numberIntegers;i++) {
	    int iColumn = integerVariable[i];
	    // get original bounds
	    double originalLower = lb[iColumn];//p->mip->lb[iColumn];
	    //double originalUpper = ub[iColumn];//p->mip->ub[iColumn];
	    
	    double value=newSolution[iColumn];
	    // if away from lower bound mark that fact
	    if (value>originalLower) {
	        //used_[iColumn]=1;
		}
	  }
	*/
	 // new solution
	 memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
	 p->lp_stat.ls_num_sols++;
	 returnCode=1;
	 printf("LS-FEAS: %i --- %f ", p->bc_index, *solutionValue);
	 *solutionValue = newSolutionValue + bestAllChange;
	 printf("%f\n", *solutionValue);
	 //printf("LS-new_opt: %f\n", new_opt);
      } else {
	 // bad solution - should not happen so debug if see message
	 printf("Local search got bad solution with %d infeasibilities"
		"summing to %g\n",
		numberBad,sumBad);
      }
    }
  }
  
  //delete [] integerVariable;

  //delete [] newSolution;
  //delete [] rowActivity;
  //delete [] way;
  //delete [] cost;
  //delete [] save;
  //delete [] mark;

  p->comp_times.ls += used_time(&total_time);
  p->lp_stat.ls_calls++;
  p->lp_stat.ls_last_call_ind = p->bc_index; 
  return returnCode;
}

/*===========================================================================*/
/*===========================================================================*/

int apply_local_search(lp_prob *p, double *solutionValue, double *colSolution,
		       double *betterSolution, double *dual_gap, double t_lb)
{

   int is_ip_feasible = FALSE; 
   
   while(true){

      char new_sol_found = FALSE;
      
      if(*dual_gap > p->par.ls_min_gap && p->par.ls_enabled &&
	 local_search(p, solutionValue, colSolution, betterSolution)){
	 memcpy(colSolution, betterSolution, DSIZE*p->lp_data->n);
	 if(*solutionValue > t_lb + 100*p->lp_data->lpetol){
	    *dual_gap = d_gap(*solutionValue, t_lb, p->mip->obj_offset,
			      p->mip->obj_sense);
	 }else{
	    *dual_gap = MIN(1e-4, 0.1*p->par.ls_min_gap);
	 }
	 new_sol_found = TRUE; 
	 is_ip_feasible = TRUE; 
      }
      
      if(!new_sol_found) break; 
   }
      
   return is_ip_feasible;
}

/*===========================================================================*/
/*===========================================================================*/

int lbranching_search(lp_prob *p, double *solutionValue, double *colSolution,
		      double *betterSolution, double t_lb)
{
   int i, j, row_ind, is_ip_feasible = FALSE;
   double moved_bd, coeff, etol = p->lp_data->lpetol;  


   double total_time = 0;
   total_time = used_time(&total_time);
   double obj_ub = *solutionValue; 
      
   MIPdesc *p_mip = p->mip;       
   int n = p_mip->n, m = p_mip->m;// + 2; // + 3?
   int nz = p_mip->nz;
   
   //double *p_obj = p_mip->obj; 
   char *p_is_int = p_mip->is_int;

   double *obj     = (double *) malloc(n * DSIZE);
   double *ub      = (double *) malloc(n * DSIZE);
   double *lb      = (double *) malloc(n* DSIZE);
   double *rhs     = (double *) malloc((m + 1) * DSIZE);
   char *is_int  = (char *)   malloc(n* CSIZE);
   double *rngval  = (double *) malloc((m + 1)* DSIZE);
   char *sense   = (char *)   malloc((m + 1) * CSIZE);
   
   int *matbeg  = (int *) malloc((n + 1) * ISIZE);
   int *matind  = (int *) malloc((nz + n) * ISIZE); 
   double *matval  = (double *) malloc((nz + n) * DSIZE);

   double *col_moved_bd = p->lp_data->tmp.d;
   int *ncol_ind = p->lp_data->tmp.i1; 
   
   for(i = 0; i < m; i++){
      rhs[i] = p_mip->rhs[i];
      rngval[i] = p_mip->rngval[i];
      sense[i] = p_mip->sense[i];
   }

   int col_nz, tot_nz = 0; 
   int srow_size = 0;
   double srow_offset = 0.0, obj_offset = 0.0; 
   double col_sol; 
   char add_r1;//, add_r2;  // first row is for search space, second is for obj
   char add_srow_offset;
   int int_cnt = 0;
   
   matbeg[0] = 0; 
   for(i = 0; i < n; i++){

      is_int[i] = p_is_int[i]; 
      tot_nz = matbeg[i]; 
      col_nz = p_mip->matbeg[i+1] - p_mip->matbeg[i]; 
      add_r1 = FALSE; //add_r2 = FALSE;
      add_srow_offset = FALSE;
      
      obj[i] = p_mip->obj[i]; 
      
      //if(obj[i] > lpetol || obj[i] < lpetol) add_r2 = TRUE;      
      
      moved_bd = 0.0;
      if(is_int[i] && (p_mip->ub[i] > p_mip->lb[i] + etol)){
	 add_r1 = TRUE; 
	 int_cnt++;
	 /* get moved_bd */
	 col_sol = floor(colSolution[i] + 100*etol); 
	 
	 //if(p_mip->ub[i] - p_mip->lb[i] > 1.0 + 0.001){
	 if(colSolution[i] > p_mip->ub[i] - etol){
	    moved_bd  = col_sol - 1.0;
	    add_srow_offset = TRUE; 
	 }else if(colSolution[i] < p_mip->lb[i] + etol){
	    moved_bd = col_sol;	       
	 }else{
	    if(obj[i] > 0.0){
	       moved_bd = col_sol - 1.0;		  	       
	       add_srow_offset = TRUE;  
	    }else{
	       moved_bd = col_sol;
	    }
	 }
	 lb[i] = 0.0;
	 ub[i] = 1.0;	 
      }else{
	 lb[i] = p_mip->lb[i];
	 ub[i] = p_mip->ub[i];
      }      

      for(j = p_mip->matbeg[i]; j < p_mip->matbeg[i+1]; j++){
	 row_ind = p_mip->matind[j]; 
	 coeff = p_mip->matval[j]; 
	 matind[tot_nz] = row_ind; 
	 matval[tot_nz] = coeff; 
	 tot_nz++;

	 //rhs[row_ind] = p_mip->rhs[row_ind];
	 //rngval[row_ind] = p_mip->rngval[row_ind]; 
	 
	 if(moved_bd > etol || moved_bd < -etol){
	    double moved_value = moved_bd*coeff;
	    rhs[row_ind] -= moved_value; 
	    //if(rngval[row_ind] > lpetol || rngval[row_ind] < -lpetol){
	    // rngval[row_ind] -= moved_value; 
	    //}	    
	 }
      }

      //if(add_r2){
      // matind[tot_nz] = m;
      // matval[tot_nz] = obj[i];	    
      // tot_nz++; 
      // col_nz++;
      //}

      if(add_r1){
	 ncol_ind[i] = tot_nz; 
	 matind[tot_nz] = m;
	 matval[tot_nz] = 1.0;
	 if(add_srow_offset){
	    matval[tot_nz] = -1.0; 
	    srow_offset -= 1.0; 
	 }
	 
	 srow_size++; 
	 tot_nz++;
	 col_nz++;	 
      }else{
	 ncol_ind[i] = -1; 
      }
      
      col_moved_bd[i] = moved_bd;
      obj_offset += obj[i]*moved_bd; 
      matbeg[i+1] = matbeg[i] + col_nz;

   }

   rngval[m] = 0.0; 
   rhs[m] = p->par.lb_search_k + srow_offset;
   sense[m] ='L'; 
   m++;

   sym_environment * env = sym_open_environment();
   //printf("HERE1\n");
   sym_explicit_load_problem(env, n, m, matbeg, matind, matval, 
			     lb, ub, is_int, obj, NULL, sense, 
			     rhs, rngval, FALSE);
   
   int node_limit = 100; 
   double gap_limit = 1.0; 
   
   sym_set_int_param(env, "node_limit", node_limit);
   sym_set_dbl_param(env, "gap_limit", gap_limit); 
   
   sym_set_dbl_param(env, "upper_bound", *solutionValue - obj_offset - 
		     p->par.granularity + etol);	

   if(p->par.lb_first_feas_enabled){
      sym_set_int_param(env, "find_first_feasible", TRUE);
   }
   
   sym_set_int_param(env, "verbosity", -2);
   //sym_set_int_param(env, "fr_enabled", 0);
   sym_set_int_param(env, "lb_enabled", 0);
   //sym_set_int_param(env, "prep_level", 0);

   //printf("\nint_cnt: %i\n", int_cnt);
   //int_cnt = 0;
   //for(i = 245; i <= 462; i++){
   // printf("x%i %f\n", i, colSolution[i]);
   // int_cnt++;
   //}
   //printf("\nint_cnt: %i\n", int_cnt);
   //sym_write_lp(env, "lb_test");
   sym_solve(env);
   
   int termcode = sym_get_status(env);
   
   if(termcode == TM_OPTIMAL_SOLUTION_FOUND ||
      termcode == PREP_OPTIMAL_SOLUTION_FOUND ||
      termcode == TM_FOUND_FIRST_FEASIBLE ||
      env->best_sol.has_sol){
      //if(env->warm_start->stat.created > 1){
      double *new_sol = p->lp_data->tmp.d + n;
      sym_get_col_solution(env, new_sol);
      for(i = 0; i < n; i++){
	 betterSolution[i] = new_sol[i] + col_moved_bd[i];
      }
      sym_get_obj_val(env, solutionValue);
      *solutionValue += obj_offset;
      is_ip_feasible = TRUE;
      
      //if(!sol_check_feasible(p->mip, betterSolution, 10*p->lp_data->lpetol)){
      //	 printf("LB - feasibility error... exiting\n");
      //	 exit(0);
      //}
   }
   
   //for(i = 245; i <= 462; i++){
   // printf("x%i %f\n", i, betterSolution[i]);
   //}
   
   p->lp_stat.lb_calls++;
   p->lp_stat.lb_last_call_ind = p->bc_index;
   p->comp_times.lb += used_time(&total_time);
   if(is_ip_feasible){
      p->lp_stat.lb_num_sols++;
      p->lp_stat.lb_last_sol_call = p->lp_stat.lb_calls;
      printf("LB-FEAS: %i --- %f %f\n", p->bc_index, obj_ub, *solutionValue);
      //p->lp_stat.lb_analyzed_nodes += analyzed_nodes;
   }

   if(env) sym_close_environment(env);
   
   return is_ip_feasible;
   
}
/*===========================================================================*/
/*===========================================================================*/   
int restricted_search(lp_prob *p, double *solutionValue, double *colSolution,
		      double *betterSolution, int fr_mode, double t_lb)
{

   // if(!(p->mip->mip_inf->binary_var_num)){
   // return FALSE;
   // }    

  LPdata *lp_data  = p->lp_data;
  double etol = lp_data->lpetol;
  int n = lp_data->n;
  int i, j, ind; 

  double total_time = 0; 
  double max_int_fixed_ratio = p->par.fr_max_int_fixed_ratio;
  double min_int_fixed_ratio = p->par.fr_min_int_fixed_ratio;
  double max_c_fixed_ratio = p->par.fr_max_c_fixed_ratio;
  double min_c_fixed_ratio = p->par.fr_min_c_fixed_ratio;
  double incr_ratio = p->par.fr_incr_ratio; 
  int first_feas_enabled = p->par.fr_first_feas_enabled;
  //int fix_option = 1;//p->par.fr_fix_option; 
  //  char *is_int = (char *)   calloc(CSIZE, n);

  total_time = used_time(&total_time);

  if(fr_mode == RINS_SEARCH){
    max_int_fixed_ratio = p->par.rs_min_int_fixed_ratio;
    max_c_fixed_ratio = p->par.rs_min_c_fixed_ratio;
  }

  int has_ub = FALSE;
  double obj_ub = 0;  

  double *x = lp_data->x;
  char * direction = lp_data->tmp.c;
  
  //int *int_ind = lp_data->tmp.i1;//(int*)malloc(sizeof(int)*n);
  //int *v_type = lp_data->tmp.i2;//(char*)malloc(sizeof(char)*n);
  //double *f_obj_coeff  = (double*)malloc(sizeof(double)*n);
  //char *is_int = env->mip->is_int;
  double * lb = const_cast <double *> (lp_data->si->getColLower());
  double * ub = const_cast <double *> (lp_data->si->getColUpper());
  double * obj = const_cast <double *> (lp_data->si->getObjCoefficients());

  double * feas_sol = lp_data->tmp.d;
  double * x_rank = lp_data->tmp.d;
  //double * x_lp_rank = lp_data->tmp.d + n;
  double * x_diff = lp_data->tmp.d + 2*n;
  
  //double floorx, ceilx;
  //int d_type, int_num = 0, c_num = 0; //frac_ip_cnt = 0, int_num = 0;
  //int * frac_ind = lp_data->tmp.i1;

  int *x_fix_cnt = lp_data->tmp.i1;
  int *x_dir_cnt = lp_data->tmp.i1 + n;
  int *x_ind = lp_data->tmp.i1 + 2*n;

  int exp_int_fix = 0, max_int_fix = 1; 
  int exp_c_fix = 0, max_c_fix = 0;
  int relax_int_num = 0, relax_c_num = 0;

  int max_sol_length = 0;
  int max_int_sol_length = 0;
  int cross_cnt = 0, feas_sol_cnt = 0;

  if(fr_mode != FR_SEARCH &&  
     fr_mode != RINS_SEARCH) return FALSE; 

  double dual_gap = 100.0;
  if(p->has_ub || *solutionValue < SYM_INFINITY/2){
    has_ub = TRUE;
    obj_ub = *solutionValue;
    dual_gap = d_gap(obj_ub, t_lb, p->mip->obj_offset,
		     p->mip->obj_sense);
  }else if(fr_mode == RINS_SEARCH) return false; 

  if(p->bc_index > 1){
    if((fr_mode == FR_SEARCH && (p->lp_stat.fr_analyzed_nodes > p->par.fr_max_nodes)) ||
       (fr_mode == RINS_SEARCH && (p->lp_stat.rs_analyzed_nodes > p->par.rs_max_nodes))){
      if((fr_mode == FR_SEARCH && p->lp_stat.fr_calls > 10) || 
	 (fr_mode == RINS_SEARCH && p->lp_stat.rs_calls > 10) || 
	 //p->tm->stat.analyzed > 100 || 
	 (has_ub && dual_gap < 5.0)){
	 //printf("FR/RS -return 1\n");
	return FALSE;       
      }
    }
  }

#if 0
  if((fr_mode == FR_SEARCH && p->lp_stat.fr_calls >= 10) || 
     (fr_mode == RINS_SEARCH && p->lp_stat.rs_calls >= 10)){
     double ratio_incr = 0.05;
     int no_impt_cnt = 0;
     if(fr_mode == FR_SEARCH){
	no_impr_cnt = p->lp_stat.fr_calls - p->lp_stat.fr_last_sol_call;
	ratio_incr = 0.05*(int)(p->lp_stat.fr_calls/10) + 0.05*(int)(no_impr_cnt/10);
     }else{
	no_impr_cnt = p->lp_stat.ds_calls - p->lp_stat.ds_last_sol_call;
	ratio_incr = 0.05*(int)(p->lp_stat.rs_calls/10) + 0.05*(int)(no_impr_cnt/10);
     }
     
     max_int_fixed_ratio = MIN(0.9, max_int_fixed_ratio + ratio_incr);
     max_c_fixed_ratio = MIN(0.4, max_c_fixed_ratio + ratio_incr);      
  }
#endif

  int int_num = 0, c_num = 0;
  int init_int_cnt = 0, init_c_cnt = 0;
  int init_fixed_int_cnt = 0, init_fixed_c_cnt = 0;
  int c_nz_cnt = 0, int_nz_cnt = 0;
  
  if(fr_mode == FR_SEARCH){
  
    for(i = 0; i < n; i++){
      x_fix_cnt[i] = 0;
      x_dir_cnt[i] = 0;
      x_diff[i] = DBL_MAX;
    }  
  
    sp_desc *sp = p->tm->sp;  
    if(sp && sp->num_solutions > 1){	    
      feas_sol_cnt += sp->num_solutions;
      sp_solution *sol;
      double val;  
      int int_sol_length;
      for(i = 0; i < sp->num_solutions; i++){
	sol = sp->solutions[i];
	if(sol->xlength > max_sol_length) max_sol_length = sol->xlength;
	int_sol_length = 0;
	memset(feas_sol, 0, DSIZE*n);
	for(j = 0; j < sol->xlength; j++){
	  //ind = sol->xind[j];
	  //val = sol->xval[j];
	  feas_sol[sol->xind[j]] = sol->xval[j];
	  if(lp_data->vars[sol->xind[j]]->is_int) int_sol_length++;
	}
	
	if(int_sol_length > max_int_sol_length) max_int_sol_length =
	  int_sol_length;
	for(ind = 0; ind < n; ind++){
	  val = feas_sol[ind];
	  if(x[ind] < val + etol && x[ind] > val - etol) {
	    x_fix_cnt[ind]--; 
	    cross_cnt++;
	  }else{
	    
	    if(val >= x[ind] + etol) x_dir_cnt[ind]++;
	    else x_dir_cnt[ind]--;
	    
	    if(fabs(val - x[ind]) < x_diff[ind]){
	      x_diff[ind] = val - x[ind]; 
	    }
	  }
	}
      }
    }
    //  if(feas_sol_cnt){
    // if(cross_cnt){
    double big_number = 1e20;
    int rank1_cnt = 0, rank2_cnt = 0, vars_eff_cnt;
    //double up_cor, down_cor, fix_cor; 
    double min_obj = DBL_MAX;
    for(ind = 0; ind < n; ind++){
      x_ind[ind] = ind;
      if(obj[ind] < min_obj) min_obj = obj[ind];
    }
    min_obj = fabs(min_obj);    
    
    for(ind = 0; ind < n; ind++){
      
      if(lp_data->vars[ind]->is_int){
	int_num++;	
	if(ub[ind] <= lb[ind] + etol){
	  init_fixed_int_cnt++;
	}
      }else{
	c_num++;
	//c_nz_cnt += p->mip->mip_inf->cols[ind].col_size; 
	if(ub[ind] <= lb[ind] + etol){	
	  init_fixed_c_cnt++;
	}
      }
      if(p->mip->mip_inf){
	vars_eff_cnt = MAX(MAX(p->mip->mip_inf->cols[ind].sos_num,
			       p->mip->mip_inf->cols[ind].col_size), 
			   p->mip->mip_inf->cols[ind].nz);
      }else{
	vars_eff_cnt = lp_data->si->getMatrixByCol()->getVectorSize(ind);
      }
      int is_int = lp_data->vars[ind]->is_int;
      double x_obj = obj[ind] + min_obj + 1e-4;
      if(x_fix_cnt[ind] < 0 || !is_int || 
	 (is_int && fabs(x[ind] - floor(x[ind] + etol)) < etol)){
        //(!is_int)){// && 
        //(x[ind] < lb[ind] + etol || x[ind] > ub[ind]
        // - etol))){
        direction[ind] = 'F';
	x_fix_cnt[ind]--;	
	if(!is_int){	  
	  if(x[ind] > lb[ind] + etol && x[ind] < ub[ind] - etol){
	    x_fix_cnt[ind]++;
	    direction[ind] ='C';
	  }else{
	    init_c_cnt++;
	    //  direction[ind] = 'C';
	  }
	}else{
	  init_int_cnt++;
	}
	
	if(x_fix_cnt[ind] < 0){
	  x_rank[ind] = x_fix_cnt[ind]*big_number*//x[ind]*
	     (x_obj)/vars_eff_cnt;	
	  if(x_fix_cnt[ind] < -1){
	    rank1_cnt++;
	  }else 
	    rank2_cnt++;
	}else{
	  x_rank[ind] = //x[ind]*
	    (x_obj)/vars_eff_cnt;
	}
	//}
	//else if((lp_data->vars[ind]->is_int &&		 
	//fabs(x[ind] - floorx(x[ind] + etol)) < etol) ||
	//		   (!lp_data->vars[ind]->is_int && 
	//	    (x[ind] < lb[ind] + etol || x[ind] > ub[ind]
	//	     - etol))){
	//  x_fix_cnt[ind]--;
	//  rank2_cnt++;
	//  x_rank[ind] = x_fix_cnt[ind]*big_number * 
	//		 (obj[ind] + min_obj + 1e-4)/vars_eff_cnt;
	//  direction[ind] = 'F';
      }else{
	direction[ind] = 'U';
	if(feas_sol_cnt > 0){
	  if(x_dir_cnt[ind] < 0 ||
	     (x_dir_cnt[ind] == 0 && x_diff[ind] < 0.0)){
	    direction[ind] = 'L';
	  }
	}else{
	  if(obj[ind] > 0.0 || (obj[ind] == 0.0 &&
				ceil(x[ind]) - x[ind] > 0.5)){
	    direction[ind] = 'L';
	  }
	}
	
	if(direction[ind] == 'U'){
	  x_rank[ind] = x_obj*
	    (ceil(x[ind]) - x[ind])/vars_eff_cnt;
	}else{
	  x_rank[ind] = x_obj*
	    (x[ind] - floor(x[ind]))/vars_eff_cnt;
	}
      }
    }
  
    double int_gap = 1.0, c_gap = 1.0; 
    
    if(p->bc_level >= 5){
      if(p->bc_level < 30){
	int_gap = 0.1*(int)(1.0*p->bc_level/5);
	c_gap = 0.05*(int)(1.0*p->bc_level/5);
      }else{
	int_gap = 0.6;
	c_gap = 0.7;
      }

      if((int_num > 10 && (init_int_cnt + init_fixed_int_cnt)/(1.0*int_num) < int_gap) ||  
	 (c_num > 100 && (init_c_cnt + init_fixed_c_cnt)/(1.0*(c_num)) < c_gap)){
	//printf("FR-EXIT: %i -- %f %f -- %f %f \n", p->bc_index, 
	//      (init_int_cnt + init_fixed_int_cnt)/(1.0*int_num), int_gap, 
	//      (init_c_cnt + init_fixed_c_cnt)/(1.0*(c_num+1)), c_gap);
	 //printf("FR -return 1\n");
	 return FALSE;
      }
      if(has_ub){
	 if(int_num - init_fixed_int_cnt - init_int_cnt > 100){ 
	    //printf("FR -return 2\n");
	    return FALSE; 
	 }
      }else{
	 if(int_num - init_fixed_int_cnt - init_int_cnt > 500){ 
	    //printf("FR -return 3\n");
	    return FALSE; 
	 }
      }
    }
  }else{
     for(ind = 0; ind < n; ind++){     
	x_rank[ind] = 0;
	x_ind[ind] = ind;
	direction[ind] = 'L';
	int add_fixed = 0;      
	if (ub[ind] < lb[ind] + etol) add_fixed = 1;
	else if(x[ind] > colSolution[ind] - etol && x[ind] < colSolution[ind] + etol){
	   x_rank[ind] = -1.0;
	   direction[ind] = 'F';
	}
	
	if(lp_data->vars[ind]->is_int) {
	   int_num++;
	   init_int_cnt += -(int)(x_rank[ind]); 
	   init_fixed_int_cnt += add_fixed; 
	}else{
	   c_num++;
	   //c_nz_cnt += p->mip->mip_inf->cols[ind].col_size; 
	   init_c_cnt += -(int)(x_rank[ind]);
	   init_fixed_c_cnt += add_fixed;
	   if(x_rank[ind] < 0) direction[ind] = 'S';
	}
     }

     int dec_fac = MIN(3, (int)(dual_gap/10.0));     
     max_int_fixed_ratio -= dec_fac*0.05; 
     
     if((int_num > 10 && (init_int_cnt + init_fixed_int_cnt)/(1.0*int_num) < max_int_fixed_ratio) || 
	(c_num > 100 && (init_c_cnt + init_fixed_c_cnt)/(1.0*(c_num+1)) < max_c_fixed_ratio)){
	//printf("RS-EXIT: %i -- %f %f -- %f %f \n", p->bc_index, 
	//     (init_int_cnt + init_fixed_int_cnt)/(1.0*int_num), p->par.rs_min_int_fixed_ratio, 
	//     (init_c_cnt + init_fixed_c_cnt)/(1.0*(c_num+1)), p->par.rs_min_c_fixed_ratio);
	//printf("RS -return 1 %f %f \n", (init_int_cnt + init_fixed_int_cnt)/(1.0*int_num), max_int_fixed_ratio);
	return FALSE;
     }    
     
     max_int_fix = init_int_cnt; 
     max_c_fix = init_c_cnt;

     if(has_ub){
	if(int_num - init_fixed_int_cnt - max_int_fix > 1000){
	   //printf("RS -return 2\n");
	   return FALSE; 
	}
     }else{
	if(int_num - init_fixed_int_cnt - max_int_fix > 2500){
	   //printf("RS -return 2\n");
	   return FALSE; 
	}
     }
  }
  
  qsort_di(x_rank, x_ind, n);

  int_num = 0, c_num = 0;
  int fixed_int_cnt = 0, fixed_c_cnt = 0;
  int *int_ind = lp_data->tmp.i1;
  int *c_ind = lp_data->tmp.i1 + n;

  char *sym_fixed_type = lp_data->tmp.c + n; /* - N don't touch
					       - F fixed to val and don't keep in model 
					       - T fixed to val and keep in model
					       - U update lower bound to val
					       - L update upper bound to val 
					    */
  double *sym_fixed_val = lp_data->tmp.d; 
  int sym_fixed_int_cnt = 0, sym_fixed_c_cnt = 0, sym_upd_int_cnt = 0;// sym_upd_c_cnt = 0;
  double sym_fixed_offset = 0.0;


  for(i = 0; i < n; i++){
     ind = x_ind[i];
     sym_fixed_type[ind] = 'N';
     sym_fixed_val[ind] = 0.0;     
     if(ub[ind] > lb[ind] + etol){ 
	if(lp_data->vars[ind]->is_int){
	   //sym_set_integer(env, ind);
	   int_ind[int_num] = ind;
	   int_num++;
	   int_nz_cnt += p->mip->mip_inf->cols[ind].col_size; 
	}else if(direction[ind] == 'F' || direction[ind] == 'S'){
	  c_ind[c_num] = ind;	  
	  c_num++;
	  c_nz_cnt += p->mip->mip_inf->cols[ind].col_size; 
	}
     }else{
	//sym_fixed_offset += obj[ind]*lb[ind]; 
	sym_fixed_type[ind] = 'B'; //coming from branching decisions 
	sym_fixed_val[ind] = lb[ind];
	if(lp_data->vars[ind]->is_int){
	   //sym_set_integer(env, ind);
	   fixed_int_cnt++;
	}else{
	   fixed_c_cnt++;
	}
     }
  }
  sym_fixed_int_cnt = fixed_int_cnt;
  sym_fixed_c_cnt = fixed_c_cnt;

  if(fr_mode == FR_SEARCH){

     max_int_fix = (int)(max_int_fixed_ratio*int_num) + 1;
     if(c_num)
	max_c_fix = (int)(max_c_fixed_ratio*c_num) + 1;

  }

  if(max_int_fix > int_num) max_int_fix = int_num;
  
  if(c_num){
     if(max_c_fix > c_num) max_c_fix = c_num;
  }
  
  relax_int_num = 0;
  relax_c_num = 0;
  exp_int_fix = max_int_fix;
  exp_c_fix = max_c_fix;

  if(fr_mode == FR_SEARCH){
    relax_int_num = max_int_fix - (int)(min_int_fixed_ratio*max_int_fix) + 1;   
    if(relax_int_num < 0) relax_int_num = 0;
    if(relax_int_num > max_int_fix) relax_int_num = max_int_fix;
    if(c_num){
      relax_c_num = max_c_fix - (int)(min_c_fixed_ratio*max_c_fix) + 1;   
      if(relax_c_num < 0) relax_c_num = 0;
      if(relax_c_num > max_c_fix) relax_c_num = max_c_fix;
    }

    if(feas_sol_cnt){
      exp_int_fix = (int)(1.1*max_int_sol_length) + 1;
      exp_c_fix = (int)(0.1*(max_sol_length - max_int_sol_length))+1;
    }
  }
  
  int nz_int_fix = 0, nz_c_fix = 0;
  double new_lb, new_ub; 
  char fix_all_z = FALSE; 

  if(fr_mode == FR_SEARCH){
     int int_updated = FALSE;
  
     if(int_nz_cnt > 10000){
	int fix_int_nz = 0;
	int new_max_int_fix = 0;
	for(i = 0; i < int_num; i++){
	   ind = int_ind[i]; 
	   fix_int_nz += p->mip->mip_inf->cols[ind].col_size;
	   new_max_int_fix = i;
	   if(int_nz_cnt - fix_int_nz < 10000){
	      break;
	   }
	}
	if(new_max_int_fix > max_int_fix){
	   max_int_fix = new_max_int_fix;
	   int_updated = TRUE;
	}
     }  
     
     if(c_nz_cnt > 50000 && false){
	int_updated = TRUE;
     }

     if(int_updated){
	if(int_num > 50 && max_int_fix >= int_num - 20) max_int_fix = int_num - 20;
	relax_int_num = max_int_fix - (int)(min_int_fixed_ratio*max_int_fix) + 1;   
	if(relax_int_num < 0) relax_int_num = 0;
	if(relax_int_num > max_int_fix) relax_int_num = max_int_fix;
     }
  }
  
  for(i = 0; i < max_int_fix; i ++){
     ind = int_ind[i];
     double bd;
     if(!fix_all_z || lb[ind] > etol || ub[ind] < -etol){
	if(direction[ind] == 'U'){        
	   bd = ceil(x[ind]);
	   if(ub[ind] < bd + etol){
	      sym_fixed_type[ind] = 'F';
	      sym_fixed_int_cnt++;
	   }else{
	      sym_fixed_type[ind] = 'U';
	      sym_upd_int_cnt++;
	   }
	   sym_fixed_val[ind] = bd;	   
	   //sym_set_col_lower(env, ind, bd);
	   new_lb = bd;
	   new_ub = ub[ind]; 
	}else if(direction[ind] == 'L'){
	   bd = floor(x[ind]);
	   if(lb[ind] > bd - etol){
	      sym_fixed_type[ind] = 'F';
	      sym_fixed_int_cnt++;
	      //sym_fixed_offset += obj[ind]*bd;
	   }else{
	      sym_fixed_type[ind] = 'L';
	      sym_upd_int_cnt++;
	   }
	   sym_fixed_val[ind] = bd;	   
	   //sym_set_col_upper(env, ind, bd);
	   new_lb = lb[ind];
	   new_ub = bd;
	}else{
	   bd = floor(x[ind] + 100*etol);
	   //sym_fixed_offset += obj[ind]*bd; 
	   sym_fixed_type[ind] = 'F';
	   sym_fixed_val[ind] = bd;
	   sym_fixed_int_cnt++;
	   //sym_set_col_lower(env, ind, bd);
	   //sym_set_col_upper(env, ind, bd);
	   new_lb = new_ub = bd; 
	}	
     }else{
	bd = 0.0;
	if(bd < lb[ind] - etol) bd = lb[ind];
	else if(bd > ub[ind] + etol) bd = ub[ind];
	//sym_fixed_offset += obj[ind]*bd; 
	sym_fixed_type[ind] = 'F';
	sym_fixed_val[ind] = bd;	   	
	sym_fixed_int_cnt++;
	//sym_set_col_lower(env, ind, bd);
	//sym_set_col_upper(env, ind, bd);      
	new_lb = new_ub = bd;
     }
     /* sanity check */
     if(bd < lb[ind] - etol || bd > ub[ind] + etol){
	printf("ERROR in FR-I routine...\n");
     }
     
     if(!fix_all_z){
	//sym_get_col_lower(env, ind, &new_lb);
	//sym_get_col_upper(env, ind, &new_ub);
	//new_lb = env->mip->lb[ind];
	//new_ub = env->mip->ub[ind];
	
	if(new_lb > etol || new_ub < -etol){
	   if(nz_int_fix++ >= exp_int_fix) fix_all_z = TRUE;
	}
     }
  }
  /*-----------------------*/

  //  max_c_fix = 0; 
  //c_num = 0; 
  //relax_c_num = 0;

  /* check if we need to fix cont cols */
  if(fr_mode == FR_SEARCH || fr_mode == RINS_SEARCH){

     int c_updated = FALSE;  
     if(c_nz_cnt > 50000){// && fr_mode == RINS_SEARCH){
	int fix_c_nz = 0;
	int new_max_c_fix = 0;
	for(i = 0; i < c_num; i++){
	   ind = c_ind[i]; 
	   fix_c_nz += p->mip->mip_inf->cols[ind].col_size;
	   new_max_c_fix = i;
	   if(c_nz_cnt - fix_c_nz < 40000){
	      break;
	   }
	}
	if(new_max_c_fix > 0 && fix_c_nz > 5000){
	   c_updated = TRUE;
	   if(fr_mode == RINS_SEARCH && new_max_c_fix > max_c_fix)
	      printf("error in c_fixing for rins...\n");
	   max_c_fix = new_max_c_fix;	
	   relax_c_num = max_c_fix - (int)(min_c_fixed_ratio*max_c_fix) + 1;   
	   if(relax_c_num < 0) relax_c_num = 0;
	   if(relax_c_num > max_c_fix) relax_c_num = max_c_fix;
	}
     }
  
     if(!c_updated){
	max_c_fix = 0;
	c_num = 0;
	relax_c_num = 0;
     }
     //}else{
     //if(c_nz_cnt < 50000){
     //	max_c_fix = 0;
     //c_num = 0;
     //	relax_c_num = 0;
     // }
  }

  //printf("nz_int - nz_c --- fix_int_nz - fix_c_nz: %i %i %i %i\n", int_nz_cnt, c_nz_cnt, fix_int_nz, fix_c_nz);
  /*-----------------------*/
  
  fix_all_z = FALSE;
  for(i = 0; i < max_c_fix; i++){
     ind = c_ind[i];
     if(ind >= lp_data->n || ind < 0) 
       printf("ERROR- %i %i %i\n",p->bc_index, ind, lp_data->n); 
     double bd; 
     //if(direction[ind] == 'C'){
     //   bd = x[ind];
     //	sym_set_col_lower(env, ind, bd);
     //	sym_set_col_upper(env, ind, bd);
     //}else{     
     if(direction[ind] == 'F'){
       if(!fix_all_z || lb[ind] > etol || ub[ind] < -etol){
	 bd = floor(x[ind] + 100*etol);
	 if(x[ind] < lb[ind] + etol) bd = lb[ind]; 
	 else bd = ub[ind];
	 //sym_set_col_lower(env, ind, bd);
	 //sym_set_col_upper(env, ind, bd);
       }else{
	 bd = 0.0;
	 if(bd < lb[ind] - etol) bd = lb[ind];
	 else if(bd > ub[ind] + etol) bd = ub[ind];
	 //sym_set_col_lower(env, ind, bd);
	 //sym_set_col_upper(env, ind, bd);
       }
       //sym_fixed_offset += obj[ind]*bd;
       sym_fixed_type[ind] = 'F';
       sym_fixed_val[ind] = bd;
       sym_fixed_c_cnt++;
       new_lb = new_ub = bd; 
       /* sanity check */
       if(bd < lb[ind] - etol || bd > ub[ind] + etol){
	 printf("ERROR in FR-C routine...\n");
       }
       
       if(!fix_all_z){
	 //sym_get_col_lower(env, ind, &new_lb);
	 //sym_get_col_upper(env, ind, &new_ub);
	  //new_lb = env->mip->lb[ind];
	  //new_ub = env->mip->ub[ind];	 
	  if(new_lb > etol || new_ub < -etol){
	     if(nz_c_fix++ >= exp_c_fix) fix_all_z = TRUE;
	  }
       }       
     }else if(direction[ind] == 'S'){
       bd = x[ind];
       if(bd < lb[ind]) bd = lb[ind];
       if(bd > ub[ind]) bd = ub[ind];
       //sym_fixed_offset += obj[ind]*bd;
       sym_fixed_type[ind] = 'F';
       sym_fixed_val[ind] = bd;
       sym_fixed_c_cnt++;
       //sym_set_col_lower(env, ind, bd);
       //sym_set_col_upper(env, ind, bd);
     }     
  }
  
  //sym_set_int_param(env, "node_limit", -1);   
  //sym_set_dbl_param(env, "time_limit", 100);

  int rest_ns_cnt = 1; 
  int imp_sol_found_cnt = 0;
  if(incr_ratio > 1.0) incr_ratio = 1.0;
  
  int unfix_int_inc_cnt = (int)(incr_ratio * max_int_fix) + 1; 
  if(unfix_int_inc_cnt < 0) unfix_int_inc_cnt = 1;
  if(unfix_int_inc_cnt > max_int_fix) unfix_int_inc_cnt = (int)(0.5*max_int_fix) + 1;

  //int orig_unfix_int_inc_cnt = unfix_int_inc_cnt; 
  //int wide_unfix_int_inc_cnt = MIN(0.1*max_int_fix, 5*unfix_int_inc_cnt);

  int unfix_c_inc_cnt = 0;
  if(c_num){
     unfix_c_inc_cnt = (int)(2*incr_ratio * max_c_fix) + 1; 
     if(unfix_c_inc_cnt < 0) unfix_c_inc_cnt = 1;
     if(unfix_c_inc_cnt > max_c_fix) unfix_c_inc_cnt = (int)(0.5*max_c_fix) + 1;
     //if(unfix_c_inc_cnt > 200) unfix_c_inc_cnt = 200;
  }

  int iter_cnt_limit = 20;
  int extended_iter_cnt_limit = 0;

  /*decide which vars we should keep in model */
  if(fr_mode == FR_SEARCH){
     if(p->bc_level < 1){
	if(p->lp_stat.fr_calls < 3){
	   iter_cnt_limit = 100;
	   //extended_iter_cnt_limit = 2;
	}else{
	   iter_cnt_limit = 50;
	   //extended_iter_cnt_limit = 1;
	}
     }

     /*first force the model be feasible - only by variable bounds- by unfixing some variables */
     int old_max_int_fix = max_int_fix;
     int old_max_c_fix = max_c_fix;

     fr_force_feasible(p, FALSE, &sym_fixed_int_cnt, &sym_fixed_c_cnt,
		       sym_fixed_type, sym_fixed_val, &max_int_fix, &max_c_fix);
     if(max_int_fix < 1) return FALSE; 

     if(max_int_fix < old_max_int_fix){
	for(i = 0, j = 0;i < int_num; i++){
	   ind = int_ind[i];
	   if(sym_fixed_type[ind] != 'N' && sym_fixed_type[ind] != 'B'){
	      int_ind[j] = int_ind[i];	      
	      j++;
	   }
	}

	int_num = j;
	
	relax_int_num = max_int_fix - (int)(min_int_fixed_ratio*max_int_fix) + 1;   
	if(relax_int_num < 0) relax_int_num = 0;
	if(relax_int_num > max_int_fix) relax_int_num = max_int_fix;
     }
     
     if(max_c_fix < old_max_c_fix){
	for(i = 0, j = 0;i < c_num; i++){
	   ind = c_ind[i];
	   if(sym_fixed_type[ind] != 'N' && sym_fixed_type[ind] != 'B'){
	      c_ind[j] = c_ind[i];
	      j++;
	   }
	}

	c_num = j;

	relax_c_num = max_c_fix - (int)(min_c_fixed_ratio*max_c_fix) + 1;
	if(relax_c_num < 0) relax_c_num = 0;
	if(relax_c_num > max_c_fix) relax_c_num = max_c_fix;
     }
     
     //int int_keep_in = MIN(sym_fixed_int_cnt + sym_upd_int_cnt, (iter_cnt_limit + 1)* unfix_int_inc_cnt);
     int int_keep_in = MIN(max_int_fix, (iter_cnt_limit + 1)* unfix_int_inc_cnt);
     relax_int_num = int_keep_in;
     //int c_keep_in = MIN(max_c_fix, (iter_cnt_limit + 1)* unfix_c_inc_cnt);
     //relax_c_num = c_keep_in;
     int c_keep_in = 0;
     relax_c_num = 0;

     for(i = 0, j = int_num - 1; i < int_keep_in && j >= 0; j--){
	ind = int_ind[j];
	if(sym_fixed_type[ind] != 'N' && sym_fixed_type[ind] != 'B'){
	   if(sym_fixed_type[ind] == 'F'){
	      sym_fixed_type[ind] = 'T';
	      sym_fixed_int_cnt--;
	   }
	   i++;
	}
     }
     for(i = 0, j = c_num - 1; i < c_keep_in && j >= 0; j--){
	ind = c_ind[j];
	if(sym_fixed_type[ind] != 'N' && sym_fixed_type[ind] != 'B'){
	   if(sym_fixed_type[ind] == 'F'){
	      sym_fixed_type[ind] = 'T';
	      sym_fixed_c_cnt--;
	   }
	   i++;
	}
     }
  }

  int *new_ind = lp_data->tmp.i1 + 2*n;
  int unfix_nz = 0;
  sym_environment * env = lp_to_sym(p, lp_data, FALSE, sym_fixed_int_cnt + sym_fixed_c_cnt, sym_fixed_type, 
				    sym_fixed_val, &sym_fixed_offset, &unfix_nz, new_ind);  
				    
  
  //int unfix_inc_cnt = (int)(0.01 * mip->n) + 1; 
  int unfix_int_cnt = 0; //b_num - max_fix;
  int unfix_c_cnt = 0;
  int n_int_unfix_cnt = 0;
  int n_c_unfix_cnt = 0;
     
  int is_ip_feasible = FALSE;

  sym_set_int_param(env, "verbosity", -2);
  //sym_set_int_param(env, "out_mode", 1);
  sym_set_int_param(env, "fr_enabled", FALSE);
  sym_set_int_param(env, "fp_enabled", -1);
  sym_set_int_param(env, "rs_enabled", FALSE);
  sym_set_dbl_param(env, "ds_min_gap", 5.0);
  sym_set_int_param(env, "ds_frequency", 10000);
  sym_set_dbl_param(env, "ls_min_gap", 0.0001);

  sym_set_int_param(env, "ds_guided_enabled", FALSE);
  sym_set_int_param(env, "ds_vlength_enabled", FALSE);
  sym_set_int_param(env, "ds_crossover_enabled", FALSE);
  sym_set_int_param(env, "ds_rank_enabled", FALSE);
  sym_set_int_param(env, "ds_fractional_enabled", FALSE);
  sym_set_int_param(env, "ds_euc_enabled", FALSE);
  sym_set_int_param(env, "fp_max_cycles", 10);
  sym_set_dbl_param(env, "fp_fix_ratio", 0.5); 
  sym_set_int_param(env, "use_branching_prep", 1);
  
  sym_set_int_param(env, "probing_max_depth", 5);
  sym_set_int_param(env, "gomory_max_depth", 5);
  sym_set_int_param(env, "generate_cgl_flowcover_cuts", 2);
  sym_set_int_param(env, "clique_max_depth", 3);
  sym_set_int_param(env, "knapsack_max_depth", 3);
  sym_set_int_param(env, "rel_br_cand_threshold", 2);
  sym_set_int_param(env, "rel_br_threshold", 2);
  //sym_set_int_param(env, "generate_cgl_probing_cuts", 2);
  //sym_set_int_param(env, "generate_cgl_gomory_cuts", 2);
  //sym_set_int_param(env, "generate_cgl_twomir_cuts", 2);
  //sym_set_int_param(env, "generate_cgl_flowcover_cuts", 2);
  //sym_set_int_param(env, "generate_cgl_clique_cuts", 2);
  //sym_set_int_param(env, "generate_cgl_knapsack_cuts", 2);
  //sym_set_int_param(env, "generate_cgl_cuts", 0);
  // sym_set_int_param(env, "prep_level", 0);
  //sym_set_int_param(env, "reduce_mip", 0);
  sym_set_int_param(env, "min_root_cut_rounds", 6); 
  sym_set_int_param(env, "max_cut_num_per_iter_root", 50); 
  sym_set_int_param(env, "use_branching_prep", 1); 

  //if(fr_mode == FR_SEARCH){
  int node_limit = 0;
  double gap_limit = dual_gap/2; 
  
  if(p->bc_level > 5){
    node_limit = MAX(100, MIN(200, 2*(int_num - max_int_fix - unfix_int_inc_cnt)));
    gap_limit = MAX(1.0, MIN(dual_gap/4, 10.0));
    if(c_num) gap_limit = MAX(2.0, MIN(dual_gap/2, 10.0));
    //sym_set_int_param(env, "node_limit", 
    //	      MAX(200, MIN(500, 2*(int_num - max_int_fix - unfix_int_inc_cnt))));
    //sym_set_dbl_param(env, "gap_limit",  2.0);
  }else{
    node_limit = MAX(200, MIN(500, 2*(int_num - max_int_fix - unfix_int_inc_cnt)));
    gap_limit = MAX(2.5, MIN(dual_gap/4, 10.0));
    if(c_num) gap_limit = MAX(4.0, MIN(dual_gap/2, 10.0));
    //sym_set_int_param(env, "node_limit", 
    //	      MAX(500, MIN(2000, 2*(int_num - max_int_fix - unfix_int_inc_cnt))));
    //sym_set_dbl_param(env, "gap_limit", 5.0);
  }

  //if(env->mip->nz > 1e5){
  //  node_limit = MIN(50, node_limit);
  // }
  
  sym_set_int_param(env, "node_limit", node_limit);
  sym_set_dbl_param(env, "gap_limit", gap_limit); 
  sym_set_int_param(env, "rs_mode_enabled", TRUE);
  if(p->bc_level < 1){
     sym_set_int_param(env, "rs_lp_iter_limit", (int)(1e9/unfix_nz));
  }else{
     sym_set_int_param(env, "rs_lp_iter_limit", MIN(2000, (int)(2.0*1e8/unfix_nz)));
  }

  if(fr_mode == FR_SEARCH && p->lp_stat.fr_calls < 1) 
    p->par.fr_max_nodes = 100*node_limit; 
  if(fr_mode == RINS_SEARCH && p->lp_stat.rs_calls < 1) 
    p->par.rs_max_nodes = 100*node_limit; 

    //}else{
    //if(p->bc_level > 10){
    //  sym_set_int_param(env, "node_limit", 1000);			
    //  sym_set_dbl_param(env, "gap_limit",  0.5);
    //}else{
    // sym_set_int_param(env, "node_limit", 2000);			
    // sym_set_dbl_param(env, "gap_limit", 1.0);
    //}   
    //}

  if(first_feas_enabled){
     sym_set_int_param(env, "find_first_feasible", TRUE);
  }
  
  //printf("START -FR1: %i\n", p->bc_index);
  
  int iter_cnt = 0, extended_iter_cnt = 0, analyzed_nodes = 0;  
  char search_extended = FALSE; 
  int orig_ind; 
  c_num = 0;

  while(int_num && imp_sol_found_cnt < rest_ns_cnt && unfix_nz < 1e5){
     //printf("%i - %i %i\n", p->bc_index, unfix_int_cnt, unfix_c_cnt);
     if(fr_mode == FR_SEARCH && iter_cnt > 0){
	n_int_unfix_cnt = unfix_int_cnt + unfix_int_inc_cnt;     
	if(n_int_unfix_cnt > relax_int_num){
	   if(n_int_unfix_cnt - unfix_int_inc_cnt <= relax_int_num ){
	      n_int_unfix_cnt = relax_int_num;
	   }else{
	      break;
	   }
	}
      
	for(i = unfix_int_cnt; i < n_int_unfix_cnt; i++){
	   orig_ind = int_ind[max_int_fix - i - 1]; 
	   ind = new_ind[orig_ind];	   
	   
	   if(ind < 0){
	      printf("error - in new indices restricted search...%i \n", max_int_fix - i - 1);
	      continue;
	   }
	   unfix_nz += env->mip->matbeg[ind + 1] - env->mip->matbeg[ind]; 
	   sym_set_col_upper(env, ind, ub[orig_ind]);
	   sym_set_col_lower(env, ind, lb[orig_ind]);
	   //env->mip->ub[b_ind[max_fix - i]] = 1.0;
	   //env->mip->lb[b_ind[max_fix - i]] = 0.0;
	}
	
	unfix_int_cnt = n_int_unfix_cnt;
	
	if(c_num && ((iter_cnt%5) == 0 || n_int_unfix_cnt <= relax_int_num)){
	   n_c_unfix_cnt = unfix_c_cnt + unfix_c_inc_cnt;     
	   
	   if(n_c_unfix_cnt > relax_c_num){
	      if(n_c_unfix_cnt - unfix_c_inc_cnt <= relax_c_num ){
		 n_c_unfix_cnt = relax_c_num;
	      }else{
		 unfix_c_cnt = n_c_unfix_cnt;
	      }	      
	      //  break;
	      //	}
	   }
	   
	   for(i = unfix_c_cnt; i < n_c_unfix_cnt; i++){
	      orig_ind = c_ind[max_c_fix - i - 1];
	      ind = new_ind[orig_ind];
	      unfix_nz += env->mip->matbeg[ind + 1] - env->mip->matbeg[ind]; 
	      sym_set_col_upper(env, ind, ub[orig_ind]);
	      sym_set_col_lower(env, ind, lb[orig_ind]);
	      //env->mip->ub[b_ind[max_fix - i]] = 1.0;
	      //env->mip->lb[b_ind[max_fix - i]] = 0.0;
	   }
	   
	   unfix_c_cnt = n_c_unfix_cnt;
	}
     }
     //printf("%i %i : %i %i -- %i %i\n", p->bc_index, iter_cnt, int_num, 
     //   max_int_fix - unfix_int_cnt, c_num, max_c_fix - unfix_c_cnt);
     
     if(has_ub) sym_set_dbl_param(env, "upper_bound", obj_ub - sym_fixed_offset - 
				  p->par.granularity + lp_data->lpetol);	

     iter_cnt++;
     sym_solve(env);
     
     int termcode = sym_get_status(env);

     if(env->warm_start){
	analyzed_nodes += env->warm_start->stat.created;
	//printf("%i node_cnt - node_limit - termcode - unfixed - fixed: %i %i %i %i %i\n", iter_cnt, 
	//     env->warm_start->stat.analyzed, p->par.fr_max_nodes, termcode, unfix_int_cnt, max_int_fix);
     }

     if(termcode == TM_OPTIMAL_SOLUTION_FOUND ||
	termcode == PREP_OPTIMAL_SOLUTION_FOUND ||
	termcode == TM_FOUND_FIRST_FEASIBLE){
	//if(env->warm_start->stat.created > 1){
	double *new_sol = lp_data->tmp.d + n; 
	sym_get_col_solution(env, new_sol);
	for(i = 0; i < n; i++){
	   if(new_ind[i] < 0){
	      betterSolution[i] = sym_fixed_val[i];
	   }else{
	      betterSolution[i] = new_sol[new_ind[i]];
	   }
	}
	sym_get_obj_val(env, solutionValue);  
	*solutionValue += sym_fixed_offset;
	imp_sol_found_cnt++;      
	//p->lp_stat.fr_num_sols++;
	is_ip_feasible = TRUE;
	//if(!sol_check_feasible(p->mip, betterSolution, 10*p->lp_data->lpetol)){
	// if(fr_mode == FR_SEARCH){
	//    printf("FR - feasibility error... exiting\n");
	// }else{
	//    printf("RS - feasibility error... exiting\n");
	// }
	// exit(0);
	//}
	//if(env->warm_start){
	//  printf("%i -- %i \n", p->bc_level, env->warm_start->stat.analyzed);
	//}
	//printf("%i : %i %i -- %i %i\n", p->bc_index, int_num, 
	//      max_int_fix - unfix_int_cnt, c_num, max_c_fix - unfix_c_cnt);
	
	break;
	//}
     }else {// if(termcode == TM_NODE_LIMIT_EXCEEDED){
	//printf("unfix_cnt - max_cnt: %i %i\n", unfix_int_cnt, max_int_fix);	
	if(unfix_nz > 1e5) break;
	if(fr_mode == FR_SEARCH && termcode != TM_NODE_LIMIT_EXCEEDED && termcode != TM_TARGET_GAP_ACHIEVED && 
	   !search_extended){
	   if(!env->warm_start || env->warm_start->stat.analyzed < 2){
	      if(!env->warm_start){
		 if(iter_cnt > iter_cnt_limit) break;
	      }else{
		 unfix_int_inc_cnt = MIN(100, (int)(0.1*max_int_fix) + 1);
		 search_extended = TRUE;
	      }
	   }else{
	      extended_iter_cnt++;
	      if(extended_iter_cnt > extended_iter_cnt_limit ||
		 1.0*unfix_int_cnt > 0.8*max_int_fix) break;	      
	   }
	}else{
	   if(search_extended){
	      extended_iter_cnt++;
	      if(extended_iter_cnt > extended_iter_cnt_limit ||
		 1.0*unfix_int_cnt > 0.5*max_int_fix) break;
	   }else{
	      break;
	   }
	}
     }
     if(n_int_unfix_cnt == relax_int_num || relax_int_num == 0) break;
  }
  
  //printf("HEUR-n-size-finished: %i \n", n);
  //if(imp_sol_found_cnt < rest_ns_cnt){
  // printf("probably all solved in root, exiting\n");
  //return FALSE;
  //}   

  sym_close_environment(env);
  //FREE(b_ind);
  //FREE(is_bin);
  //FREE(f_obj_coeff);
  //FREE(x);
  //FREE(fix_direction);
  //printf("END -FR: %i\n", p->bc_index);
  
  if(fr_mode == FR_SEARCH){
     p->lp_stat.fr_calls++;
     p->lp_stat.fr_last_call_ind = p->bc_index; 
     p->comp_times.fr += used_time(&total_time);
     if(is_ip_feasible){
	p->lp_stat.fr_num_sols++;
	p->lp_stat.fr_last_sol_call = p->lp_stat.fr_calls;
	printf("FR-FEAS: %i --- %f %f\n", p->bc_index, has_ub? obj_ub : 0.0, *solutionValue);
     }
     p->lp_stat.fr_analyzed_nodes += analyzed_nodes; 
  }else{
     p->lp_stat.rs_calls++;
     p->lp_stat.rs_last_call_ind = p->bc_index; 
     p->comp_times.rs += used_time(&total_time);
     if(is_ip_feasible){
	p->lp_stat.rs_num_sols++;
	p->lp_stat.rs_last_sol_call = p->lp_stat.rs_calls;
	printf("RS-FEAS: %i --- %f %f\n", p->bc_index, has_ub? obj_ub : 0.0, *solutionValue);
     }
     p->lp_stat.rs_analyzed_nodes += analyzed_nodes; 
  }
  
  return is_ip_feasible;
  
}
/*===========================================================================*/
/*===========================================================================*/

int fr_force_feasible(lp_prob *p, char use_base, int *sym_fixed_int_cnt, int *sym_fixed_c_cnt, 
		      char *sym_fixed_type, double *sym_fixed_val, int *max_int_fix, int *max_c_fix)
{

   /* now update max_int_fix and max_c_fix so that we get something feasible at least by var bounds */
   LPdata *lp_data = p->lp_data; 
   int n = lp_data->n;
   int m; 
   if(!use_base){
      m  = lp_data->m;      
   }else{
      m = p->mip->m; //p->base.cutnum;
   }
   
   double *row_lb = lp_data->tmp.d + n;
   double *row_ub = lp_data->tmp.d + n + m;     
   double new_lb, new_ub, coeff, etol = lp_data->lpetol;   
   int i, j, r_ind, c_ind, v_end; 

   const double *ub = lp_data->si->getColUpper();
   const double *lb = lp_data->si->getColLower();
   
   const CoinPackedMatrix * matrix = lp_data->si->getMatrixByCol();
   const double *matval = matrix->getElements();
   const int *matind = matrix->getIndices();
   const int *matbeg = matrix->getVectorStarts();
   const int *len = matrix->getVectorLengths();
   
   const double *si_row_ub = lp_data->si->getRowUpper();
   const double *si_row_lb = lp_data->si->getRowLower();
   const double inf = lp_data->si->getInfinity();
   for(r_ind = 0; r_ind < m; r_ind++){
      row_ub[r_ind] = row_lb[r_ind] = 0.0;
   }
   
   for(c_ind = 0; c_ind < n; c_ind++){
      v_end = matbeg[c_ind] + len[c_ind];
      new_lb = lb[c_ind];
      new_ub = ub[c_ind];
      if(sym_fixed_type[c_ind] == 'L'){
	 new_ub = sym_fixed_val[c_ind];
      }else if(sym_fixed_type[c_ind] == 'U'){
	 new_lb = sym_fixed_val[c_ind];
      }else if(sym_fixed_type[c_ind] == 'F'){
	 new_lb = new_ub = sym_fixed_val[c_ind];
      }
      
      for(i = matbeg[c_ind]; i < v_end; i++){
	 r_ind = matind[i];
	 coeff = matval[i];
	 if(row_ub[r_ind] < inf){
	    if(coeff >= 0.0 && new_ub < inf){
	       row_ub[r_ind] += new_ub*coeff;
	    }else if(coeff < 0.0 && new_lb > -inf){
	       row_ub[r_ind] += new_lb*coeff;
	    }else{
	       row_ub[r_ind] = inf;
	    }
	 }
	 if(row_lb[r_ind] > -inf){
	    if(coeff >= 0.0 && new_lb > -inf){
	       row_lb[r_ind] += new_lb*coeff;
	    }else if(coeff < 0.0 && new_ub < inf){
	       row_lb[r_ind] += new_ub*coeff;
	    }else{
	       row_lb[r_ind] = -inf;
	    }
	 }
      }
   }

   /*now we have row bounds, check if they satisfy the constraints */
   char * is_row_violated = lp_data->tmp.c + 2*n;   
   int * violated_row = lp_data->tmp.i1 + 2*n;
   
   for(i = 0, r_ind = 0; r_ind < m; r_ind++){
      if(row_lb[r_ind] > si_row_ub[r_ind] + etol ||
	 row_ub[r_ind] < si_row_lb[r_ind] - etol){
	 is_row_violated[r_ind] = TRUE;
	 violated_row[i] = r_ind;
	 i++;
      }else{
	 is_row_violated[r_ind] = FALSE;
      }
   }

   const double *r_matval;
   const int *r_matind; 
   const int *r_matbeg;
   const int *r_len;
   
   if(i > 0){
      matrix = lp_data->si->getMatrixByRow();
      r_matval = matrix->getElements();
      r_matind = matrix->getIndices();
      r_matbeg = matrix->getVectorStarts();
      r_len = matrix->getVectorLengths();
   }else{
      return 0;
   }

   int violated_cnt = i; 
   int unfix_var_cnt; 
   int * unfix_var = lp_data->tmp.i1 + 2*n + m; 
   char * is_row_added = lp_data->tmp.c + 2*n + m;
   
   while(violated_cnt > 0 && *max_int_fix + *max_c_fix > 0){
      unfix_var_cnt = 0 ;
      for(j = 0; j < violated_cnt; j++){	 
	 r_ind = violated_row[j];
	 if(!is_row_violated[r_ind]) continue;
	 v_end = r_matbeg[r_ind] + r_len[r_ind];
	 char violated_found = FALSE; 
	 for(i = r_matbeg[r_ind]; i < v_end; i++){
	    c_ind = r_matind[i];
	    if(sym_fixed_type[c_ind] != 'N' && sym_fixed_type[c_ind] != 'B'){
	       violated_found = TRUE;
	       break;
	    }
	 }
	 if(!violated_found){
	    printf("error in restricted_search - fr-feas row_violation evaluation...\n");
	    break;
	 }
	 
	 unfix_var[unfix_var_cnt] = c_ind;
	 unfix_var_cnt++;
      }
      double old_lb, old_ub;
      int new_violated_cnt = 0;
      int new_r_ind;
      memset(is_row_added, 0, CSIZE*m);
      for(i = 0; i < unfix_var_cnt; i++){
	 c_ind = unfix_var[i];
	 if(sym_fixed_type[c_ind] == 'N' || sym_fixed_type[c_ind] == 'B') continue; 
	 
	 old_lb = lb[c_ind];
	 old_ub = ub[c_ind];
	 
	 if(sym_fixed_type[c_ind] == 'F'){
	    old_lb = old_ub = sym_fixed_val[c_ind];
	    if(lp_data->vars[c_ind]->is_int) (*sym_fixed_int_cnt)--;
	    else (*sym_fixed_c_cnt)--;
	 }else if(sym_fixed_type[c_ind] == 'U'){
	    old_lb = sym_fixed_val[c_ind];
	 }else{
	    old_ub = sym_fixed_val[c_ind];
	 }
	 
	 /*unfix this variable */
	 
	 sym_fixed_type[c_ind] = 'N';
	 if(lp_data->vars[c_ind]->is_int) (*max_int_fix)--;
	 else (*max_c_fix)--;
	 
	 if(ub[c_ind] < inf && ub[c_ind] < old_ub - etol)
	    printf("error with ub in restricted_search violation evaluation..\n");
	 if(lb[c_ind] > -inf && lb[c_ind] > old_lb + etol)
	    printf("error with lb in restricted_search violation evaluation..\n");
	 
	 v_end = matbeg[c_ind] + len[c_ind];
	 for(j = matbeg[c_ind]; j < v_end; j++){
	    new_r_ind = matind[j]; 
	    coeff = matval[j]; 
	    if(is_row_violated[new_r_ind]){
	       if(row_ub[new_r_ind] < inf){
		  if(coeff >= 0.0 && ub[c_ind] < inf){
		     row_ub[new_r_ind] += coeff*(ub[c_ind] - old_ub);
		  }else if(coeff < 0.0 && lb[c_ind] > -inf){
		     row_ub[new_r_ind] += coeff*(lb[c_ind] - old_lb);
		  }else{
		     row_ub[new_r_ind] = inf;
		  }
	       }
	       if(row_lb[new_r_ind] > -inf){
		  if(coeff >= 0.0 && lb[c_ind] > -inf){
		     row_lb[new_r_ind] += coeff*(lb[c_ind] - old_lb);
		  }else if(coeff < 0.0 && ub[i] < inf){
		     row_lb[new_r_ind] += coeff*(ub[c_ind] - old_ub);
		  }else{
		     row_lb[new_r_ind] = -inf;
		  }
	       }
	       if((row_lb[new_r_ind] > -inf && si_row_ub[new_r_ind] < inf &&
		   row_lb[new_r_ind] > si_row_ub[new_r_ind] + etol) ||
		  (row_ub[new_r_ind] < inf && si_row_lb[new_r_ind] > -inf &&
		   row_ub[new_r_ind] < si_row_lb[new_r_ind] - etol)){
		  if(!is_row_added[new_r_ind]){
		     violated_row[new_violated_cnt] = new_r_ind;
		     new_violated_cnt++;
		     is_row_added[new_r_ind] = TRUE; 
		  }
	       }else{
		  is_row_violated[new_r_ind] = FALSE;
	       }
	    }
	 }
      }
      
      if(new_violated_cnt > violated_cnt){
	 printf("error in new_violated_row_cnt...\n");	      
      }else{
	 violated_cnt = new_violated_cnt; 
      }
   }

   return 0;
}

/*===========================================================================*/
/*===========================================================================*/

sym_environment * lp_to_sym(lp_prob *p, LPdata *lp_data, char use_base, int sym_fixed_cnt,
			    char *sym_fixed_type, double *sym_fixed_val, 
			    double *sym_fixed_offset, int *unfix_nz, int *new_ind)
{   

   int n = lp_data->n - (sym_fixed_cnt = MAX(sym_fixed_cnt, 0));
   int m, nz; 
   if(!use_base){
      m  = lp_data->m;
      nz = lp_data->nz;
   }else{
      m = p->mip->m; //p->base.cutnum;
      nz = p->mip->nz;
   }
   
   if(lp_data->n != lp_data->si->getNumCols()){
      printf("col-error %i %i\n", n, lp_data->si->getNumCols());
   }
   if(!use_base){
      if(m != lp_data->si->getNumRows()){
	 printf("row-error %i %i\n", m, lp_data->si->getNumRows());
      }
      if(nz != lp_data->si->getNumElements()){
	 printf("nz-error %i %i\n", nz, lp_data->si->getNumElements());
      }
   }

   *unfix_nz = 0; 
   char * is_row_used = lp_data->tmp.c + 2*(lp_data->n);   

   sym_environment * env = sym_open_environment();
   
   double *obj    = (double *) malloc(DSIZE * n);
   double *rhs    = (double *) malloc(DSIZE * m);
   char *sense  = (char *)   malloc(CSIZE * m);
   double *rngval = (double *) malloc(DSIZE * m);
   double *ub     = (double *) malloc(DSIZE * n);
   double *lb     = (double *) malloc(DSIZE * n);
   char *is_int  = (char *)   malloc(CSIZE * n);
   
   int *matbeg = (int *) malloc(ISIZE * (n + 1));
   double *matval = (double *) malloc(DSIZE*nz);
   int *matind = (int *)    malloc(ISIZE*nz);

   const double *si_obj = lp_data->si->getObjCoefficients();
   const double *si_ub = lp_data->si->getColUpper();
   const double *si_lb = lp_data->si->getColLower();
   const double *si_rhs = lp_data->si->getRightHandSide();
   const char *si_sense = lp_data->si->getRowSense();
   const double *si_rngval = lp_data->si->getRowRange();

   const CoinPackedMatrix * matrixByCol= lp_data->si->getMatrixByCol();
   const int *si_matbeg = matrixByCol->getVectorStarts();
   const double *si_matval = matrixByCol->getElements();
   const int *si_matind = matrixByCol->getIndices();
   const int *si_length = matrixByCol->getVectorLengths();   
   int i, j, n_nz, si_end; 
   
   for(i = 0; i < m; i++){
      is_row_used[i] = FALSE;
      rhs[i] = si_rhs[i];
      sense[i] = si_sense[i];
      rngval[i] = si_rngval[i];
   }
   
   //int col_length; 
   matbeg[0] = 0;
   n_nz = 0;
   /*
     if(sym_fixed_cnt < 1){
     
     for(i = 0; i < n; i++){
     obj[i] = si_obj[i];
     ub[i] = si_ub[i];
     lb[i] = si_lb[i];
     
     si_end = si_matbeg[i] + si_length[i];
     col_length = n_nz;
     for(j = si_matbeg[i]; j < si_end; j++){
     if(si_matind[j] < m){
     matval[n_nz] = si_matval[j];
     matind[n_nz] = si_matind[j];
     n_nz++;	       
     }
     }
     is_int[i] = lp_data->vars[i]->is_int; 
     matbeg[i+1] = n_nz;	 
     }
     }else{
   */

   int new_n = 0, used_row_cnt = 0;
   char add_nz; 
   matbeg[0] = 0;
   for(i = 0; i < lp_data->n; i++){
      add_nz = TRUE; 
      if(sym_fixed_type[i] != 'F' && sym_fixed_type[i] != 'B'){
	 obj[new_n] = si_obj[i];
	 if(sym_fixed_type[i] == 'N'){
	    lb[new_n] = si_lb[i];
	    ub[new_n] = si_ub[i];
	 }else if(sym_fixed_type[i] == 'T'){
	    lb[new_n] = ub[new_n] = sym_fixed_val[i];
	    add_nz = FALSE;
	 }else if(sym_fixed_type[i] == 'L'){
	    lb[new_n] = si_lb[i];
	    ub[new_n] = sym_fixed_val[i];	    
	 }else if(sym_fixed_type[i] == 'U'){
	    lb[new_n] = sym_fixed_val[i];
	    ub[new_n] = si_ub[i];
	 }else{
	    printf("error in lp_to_sym... %c \n", sym_fixed_type[i]);
	 }
	 
	 si_end = si_matbeg[i] + si_length[i];
	 for(j = si_matbeg[i]; j < si_end; j++){
	    if(si_matind[j] < m){
	       matval[n_nz] = si_matval[j];
	       matind[n_nz] = si_matind[j];
	       n_nz++;
	       if(!is_row_used[si_matind[j]]){
		  used_row_cnt++;
		  is_row_used[si_matind[j]] = TRUE;
	       }
	    }
	 }
	 if(add_nz) (*unfix_nz) += n_nz - matbeg[new_n];
	 is_int[new_n] = lp_data->vars[i]->is_int;	    
	 matbeg[new_n + 1] = n_nz;
	 new_ind[i] = new_n; 
	 new_n++;
      }else{
	 *sym_fixed_offset += si_obj[i] * sym_fixed_val[i];
	 new_ind[i] = -1; 
	 si_end = si_matbeg[i] + si_length[i];
	 for(j = si_matbeg[i]; j < si_end; j++){
	    if(si_matind[j] < m){
	       rhs[si_matind[j]] -= si_matval[j]*sym_fixed_val[i];
	    }
	 }
      }
   }
      
   if(new_n != n){
      printf("error in num cols in lp_to_sym...\n");
   }
   // }

   if(used_row_cnt < m){
      int *new_row_ind = lp_data->tmp.i1 + 3*(lp_data->n);
      double etol = lp_data->lpetol;
      
      for(i = 0, j = 0; j < m; j++){
	 new_row_ind[j] = -1;
	 if(is_row_used[j]){
	    rhs[i] = rhs[j];
	    sense[i] = sense[j];
	    rngval[i] = rngval[j];
	    new_row_ind[j] = i;
	    i++;
	 }else{
	    int is_feas = TRUE;
	    switch(sense[j]){
	     case 'E':
	       if(rhs[j] > etol || rhs[j] < -etol) is_feas = FALSE; 
	       break;
	     case 'L':
	       if(rhs[j] < -etol) is_feas = FALSE;
	       break;
	     case 'G':
	       if(rhs[j] > etol) is_feas = FALSE;
	       break;
	     case 'R':	       
	       if(si_rhs[j] - rhs[j] < lp_data->si->getRowLower()[j] - etol ||
		  si_rhs[j] - rhs[j] > lp_data->si->getRowUpper()[j] + etol) is_feas = FALSE;
	       break;
	     default:
	       break;
	    }

	    if(!is_feas)
	       printf("error in lp_to_sym feasibility %i %i\n", p->bc_index, j);
	 }
      }

      if(i != used_row_cnt) printf("error in lp_to_sym feasibility row num: %i %i \n", i, used_row_cnt);
      for(i = 0; i < n_nz; i++){
	 if(new_row_ind[matind[i]] < 0){
	    printf("error in lp_to_sym feasibility row new index: %i \n", i);
	    continue;
	 }
	 matind[i] = new_row_ind[matind[i]];
      }      
      m = used_row_cnt; 
   }
   
   sym_explicit_load_problem(env, n, m, matbeg, matind, matval, 
			     lb, ub, is_int, obj, NULL, sense, 
			     rhs, rngval, FALSE);
   
   return env;   
}

#if 0
sym_environment * lp_to_sym(lp_prob *p, LPdata *lp_data, char use_base)
{
  
  int n = lp_data->n;
  int m, nz; 
  if(!use_base){
    m  = lp_data->m;
    nz = lp_data->nz;
  }else{
    m = p->mip->m; //p->base.cutnum;
    nz = p->mip->nz;
  }

  if(n != lp_data->si->getNumCols()){
    printf("col-error %i %i\n", n, lp_data->si->getNumCols());
  }
  if(!use_base){
    if(m != lp_data->si->getNumRows()){
      printf("row-error %i %i\n", m, lp_data->si->getNumRows());
    }
    if(nz != lp_data->si->getNumElements()){
      printf("nz-error %i %i\n", nz, lp_data->si->getNumElements());
    }
  }

  sym_environment * env = sym_open_environment();
  
  double *obj    = (double *) malloc(DSIZE * n);
  double *rhs    = (double *) malloc(DSIZE * m);
  char *sense  = (char *)   malloc(CSIZE * m);
  double *rngval = (double *) malloc(DSIZE * m);
  double *ub     = (double *) malloc(DSIZE * n);
  double *lb     = (double *) malloc(DSIZE * n);

  int *matbeg = (int *) malloc(ISIZE * (n + 1));
  double *matval = (double *) malloc(DSIZE*nz);
  int *matind = (int *)    malloc(ISIZE*nz);


  memcpy(obj, const_cast <double *> (lp_data->si->getObjCoefficients()),
	 DSIZE * n); 
  memcpy(ub, const_cast <double *> (lp_data->si->getColUpper()),
	 DSIZE * n); 
  memcpy(lb, const_cast <double *> (lp_data->si->getColLower()),
	 DSIZE * n); 
  memcpy(rhs, const_cast <double *> (lp_data->si->getRightHandSide()),
	 DSIZE * m); 
  memcpy(sense, const_cast <char *> (lp_data->si->getRowSense()),
	 CSIZE * m); 
  memcpy(rngval, const_cast <double *> (lp_data->si->getRowRange()),
	 DSIZE * m); 
  
  if(!use_base){

    int i, j, n_nz, si_end; 
  
    const CoinPackedMatrix * matrixByCol= lp_data->si->getMatrixByCol();
    const int *si_matbeg = matrixByCol->getVectorStarts();
    const double *si_matval = matrixByCol->getElements();
    const int *si_matind = matrixByCol->getIndices();
    const int *si_length = matrixByCol->getVectorLengths();

    matbeg[0] = 0;
    n_nz = 0;
    for(i = 0; i < n; i++){
      matbeg[i+1] = matbeg[i] + si_length[i]; 
      si_end = si_matbeg[i] + si_length[i];
      for(j = si_matbeg[i]; j < si_end; j++){
	matval[n_nz] = si_matval[j];
	matind[n_nz] = si_matind[j];
	n_nz++;
      }
    }
  }else{

    //memcpy(rhs, p->mip->rhs, DSIZE * m); 
    //memcpy(sense, p->mip->sense, CSIZE * m); 
    //memcpy(rngval, p->mip->rngval, DSIZE * m); 
    memcpy(matbeg, p->mip->matbeg, ISIZE * (n + 1));
    memcpy(matval, p->mip->matval, DSIZE * nz); 	   
    memcpy(matind, p->mip->matind, ISIZE * nz);  	   
  }

   //   char ** colname = (char **) malloc(sizeof(char *) * n);  

   //for (j = 0; j < n; j++){
   //  if (lp_data->vars[j]->is_int) {
   //	is_int[j] = TRUE;
   //   }
	//colname[j] = (char *) malloc(CSIZE * 9);
	//strncpy(colname[j], p->mip->colname[j], 9);
	//colname[j][8] = 0;
   // }

   sym_explicit_load_problem(env, n, m, matbeg, matind, matval, 
			     lb, ub, NULL, obj, NULL, sense, 
			     rhs, rngval, FALSE);

   return env; 
}
#endif
/*===========================================================================*/
/*===========================================================================*/
int ds_fix_common_vars(LPdata * lp_data, var_desc **vars, double *ip_sol, double *x)
{
   int i, n = lp_data->n;
   double ub, lb, etol = lp_data->lpetol;
   
   for(i = 0; i < n; i++){
      get_ub(lp_data, i, &ub);
      get_lb(lp_data, i, &lb);
      if(ub > lb + etol){
	 if(ip_sol){
	    if(x[i] < ip_sol[i] + etol && x[i] > ip_sol[i] - etol){
	       change_lbub(lp_data, i, ip_sol[i], ip_sol[i]);
	    }
	 }else{
	    if(vars[i]->is_int){
	       double bd = floor(x[i] + etol); 
	       if(fabs(x[i] - bd) < etol){
		  change_lbub(lp_data, i, bd, bd);
	       }
	    }
	 }
      }
   }
   

   return 0;
}
/*===========================================================================*/
/*===========================================================================*/

#if 0
double get_dgap(double obj_ub, double obj_lb, double obj_offset, char obj_sense){

  double t_ub = obj_ub + obj_offset, t_lb = obj_lb + obj_offset;

  if(obj_sense == SYM_MAXIMIZE){
    t_lb -= (obj_ub + obj_lb);
    t_ub -= (obj_lb + obj_ub);
  }
  
  return (t_ub ? (t_ub - t_lb)/fabs(t_ub)*100 : 100.0);

}
#endif







