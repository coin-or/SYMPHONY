/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2008 Ted Ralphs. All Rights Reserved.                  */
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
#include "sym_macros.h"
#include "sym_return_values.h"

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

int feasibility_pump (lp_prob *p, char *found_better_solution,
      double &solution_value, double *betterSolution)
{
   int                      termcode    = FUNCTION_TERMINATED_NORMALLY;
   FPdata                  *fp_data     = (FPdata*) malloc(sizeof(FPdata));
   LPdata                  *lp_data     = p->lp_data;
   LPdata                  *new_lp_data = (LPdata *)calloc(1,sizeof(LPdata));	
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
   double                   fp_time, real_obj_value, target_ub;
   FPvars                 **vars;
   double                   gap           = model->getInfinity();
   double                   obj_lb        = lp_data->objval;
   double                   total_time    = 0;
   const double            *mip_obj       = model->getObjCoefficients();
   char                     is_feasible   = FALSE;
   double                  *x_ip, *x_lp, new_solution_value;

   fp_time                                = used_time(&total_time);
   /* total_time and fp_time both now have total time used by symphony's lp
    * process */
   fp_time                                = used_time(&total_time);
   /* fp_time should now be zero and total_time be still the same */

   *found_better_solution = FALSE;
   verbosity = fp_data->verbosity     = p->par.verbosity;
   //verbosity = 10;
   fp_data->mip_obj       = (double *)malloc(n*DSIZE);
   fp_data->flip_fraction = p->par.fp_flip_fraction;
   memcpy(fp_data->mip_obj,mip_obj,n*DSIZE);

   /* initialize the lp solver. load the current basis */
   fp_initialize_lp_solver(p, new_lp_data, fp_data);
   x_ip = fp_data->x_ip;
   x_lp = fp_data->x_lp;
   if (p->has_ub) {
      solution_value = p->ub-p->mip->obj_offset;
      fp_add_obj_row(new_lp_data, n, mip_obj, p->ub-p->par.granularity);
   } else {
      solution_value = model->getInfinity();
   }

   /* round the x_lp and store as x_ip, it will usually become infeasible */
   vars = fp_data->fp_vars;

   /* do the following max_iter times */
   fp_time += used_time(&total_time);
   for (iter=0; iter<max_iter && fp_time<p->par.fp_time_limit; iter++) {
      PRINT(verbosity,5,("fp: iteration %d\n",iter));
      is_feasible = FALSE;
      /* solve an lp */
      fp_round(fp_data, new_lp_data);
      if (fp_data->x_bar_len[fp_data->iter] == -1) {
         /*
          * the cost and reference point are same as some other iteration. we
          * should stop here because we are cycling
          */
         PRINT(verbosity,5,("fp: leaving because of cycling\n"));
         break;
      }
      fp_is_feasible (lp_data,matrix,lp_r_low,lp_r_up,fp_data,&is_feasible);

      if (is_feasible == TRUE) {
         /* we found what we wanted */
         memcpy(betterSolution, x_ip, n*DSIZE);

         new_solution_value = 0;
         for (i=0;i<n;i++) {
            new_solution_value += betterSolution[i]*mip_obj[i];
         }
         if (new_solution_value<solution_value-p->par.granularity) {
            solution_value = new_solution_value;
            indices = p->lp_data->tmp.i1;          /* n */
            values  = p->lp_data->tmp.d;           /* n */
            cnt     = collect_nonzeros(p, betterSolution, indices, values);
            gap     = (solution_value -
                      obj_lb)/(fabs(solution_value)+0.001)*100;
            p->lp_stat.fp_num_sols++;
            PRINT(verbosity,5,("fp: found solution with value = %f\n",
                     solution_value));
            PRINT(verbosity,5,("fp: gap = %f\n", gap));
            sp_add_solution(p,cnt,indices,values,
                  solution_value+p->mip->obj_offset,p->bc_index);
            if (gap <= p->par.fp_min_gap) {
               *found_better_solution = TRUE;
               break;
            }
            target_ub = (obj_lb + solution_value)/2;
            if (*found_better_solution != TRUE && p->has_ub==FALSE) {
               // add another objective function constraint to lower the
               // objective value.
               fp_add_obj_row(new_lp_data, n, mip_obj, target_ub);
            } else {
               r = new_lp_data->m-1;
               change_rhs(new_lp_data, 1, &r, &target_ub);
            }
            *found_better_solution = TRUE;
         }
      } 

      PRINT(verbosity,10,("fp: solve lp %d\n",iter));
      p->lp_stat.lp_calls++;
      if (fp_solve_lp(new_lp_data, fp_data, &is_feasible) != 
            FUNCTION_TERMINATED_NORMALLY) {
         break;
      }

      fp_data->iter++;
      fp_time += used_time(&total_time);
   }
   close_lp_solver(new_lp_data);
   /* free all the allocated memory */
   FREE(new_lp_data->x);
   FREE(new_lp_data->lb);
   FREE(new_lp_data->ub);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->dualsol);
   FREE(new_lp_data->dj);
   FREE(new_lp_data->slacks);
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
   FREE(fp_data->obj);
   FREE(fp_data->mip_obj);
   FREE(fp_data->x_lp);
   FREE(fp_data->x_ip);
   FREE(fp_data->index_list);
   FREE(fp_data->x_bar_len);
   FREE(fp_data->x_bar_val);
   FREE(fp_data->x_bar_ind);
   FREE(fp_data->alpha_p);
   FREE(fp_data);

   /* update stats */
   fp_time                        += used_time(&total_time);
   p->comp_times.fp               += fp_time;
   p->lp_stat.fp_calls++;
   if (*found_better_solution==TRUE) {
      if (p->mip->obj_sense == SYM_MAXIMIZE){
         real_obj_value=-solution_value+p->mip->obj_offset;
      } else {
         real_obj_value=solution_value+p->mip->obj_offset;
      }
      PRINT(verbosity,5,("fp: found solution = %10.2f time = %10.2f\n",
               real_obj_value,total_time));
   }

   PRINT(verbosity,5,("Leaving Feasibility Pump.\n"));
    //exit(0);
   return termcode;
}


/*===========================================================================*/
int fp_is_feasible (LPdata *lp_data, const CoinPackedMatrix *matrix, const double *r_low, const double *r_up, FPdata *fp_data, char *is_feasible )
{
   /* check if x is a integer feasible solution to problem in p */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double lpetol = lp_data->lpetol;
   int n = fp_data->n0;
   int m = fp_data->m0;
   FPvars **vars = fp_data->fp_vars;
   int i,c,j;
   double Ractivity;
   const int *r_matbeg = matrix->getVectorStarts();
   const int *r_matlen = matrix->getVectorLengths();
   const int *r_matind = matrix->getIndices();
   const double *r_matval = matrix->getElements();
   double *x = fp_data->x_ip;

   *is_feasible = TRUE;
   for (i=0;i<n;i++) {
      if (vars[i]->is_int) {
         if (x[i]-floor(x[i])>lpetol && ceil(x[i])-x[i]>lpetol) {
            /* some int variable is non-integral */
            /* is not possible, since this function is called after rounding */
            printf("Bok!\n");
            return FUNCTION_TERMINATED_ABNORMALLY;
         }
      }
   }

   /* check feasibility of constraints */
   for (i=0;i<m;i++) {
      Ractivity = 0;
      c=0;			/* column */
      for (j=r_matbeg[i];j<r_matbeg[i]+r_matlen[i];j++) {
         c=r_matind[j];
         Ractivity = Ractivity+x[c]*r_matval[j];
      }
      //      printf("Ractivity[%d] = \t%f\n",i,Ractivity);
      if (Ractivity>r_up[i]+lpetol || Ractivity<r_low[i]-lpetol) {
         /* constraint infeasibility is possible since we call this func. after
            rounding */
         *is_feasible = FALSE;
         break;
      }
   }

   return termcode;
}

/*===========================================================================*/
int fp_initialize_lp_solver(lp_prob *p, LPdata *new_lp_data, FPdata *fp_data)
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
   int i;
   int *rstat,*cstat;

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

   /* used because we can not call si directly */
   copy_lp_data(lp_data,new_lp_data);
   lp_lb = new_lp_data->lb;
   lp_ub = new_lp_data->ub;

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
   rstat = (int *) malloc(m * ISIZE);
   cstat = (int *) malloc(n * ISIZE);

   get_basis(lp_data,cstat,rstat);
   load_basis (new_lp_data,cstat,rstat);

   FREE(rstat);
   FREE(cstat);

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

   /* used by change_rhs */
   new_lp_data->tmp.c = (char *)malloc(2*CSIZE);
   new_lp_data->tmp.d = (double *)malloc(DSIZE*n);
   new_lp_data->tmp.i1 = (int *)malloc(ISIZE*n);

   FREE(rmatval);
   FREE(rmatind);
   FREE(cmatbeg);
   FREE(rmatbeg);

   return termcode;
}

/*===========================================================================*/
int fp_solve_lp(LPdata *lp_data, FPdata *fp_data, char* is_feasible) 
{
   /* construct an lp based on x_ip. solve it. store the result in x_lp */
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double *objcoeff= fp_data->obj;
   int n = fp_data->n;
   int iterd;
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

   is_feasible = FALSE;
   memset ((char *)(objcoeff),0,DSIZE*n);
   for (i=0;i<fp_data->n0;i++) {
      if (fp_vars[i]->is_int) {
         if (fp_vars[i]->is_bin) {
            if (x_ip[i]==0) {
               objcoeff[i] = 1.0;
            } else if (x_ip[i]==1) {
               objcoeff[i] = -1.0;
            }
         } else {
            objcoeff[i] = 0;
            objcoeff[fp_vars[i]->xplus] = 1;
            norm += 1.0; /* stays the same every iteration */
         }
      } else {
         objcoeff[i] = 0;
      }
      /* calculate ||coeff||, norm is not zero because otherwise x_ip is
       * feasible */
      norm += objcoeff[i]*objcoeff[i]; /* stays the same every iteration */
   }

   norm = sqrt(norm);
   //norm = 0;
   PRINT(verbosity, 15, ("fp: norm = %f\n",norm));
   for (i=0;i<fp_data->n0;i++) {
      objcoeff[i] = (1-alpha)*objcoeff[i]+alpha*mip_obj[i]*norm;
   }
   for (i=fp_data->n0;i<fp_data->n;i++) {
      objcoeff[i] = (1-alpha)*objcoeff[i];
   }
   alpha = alpha*fp_data->alpha_decr;

   change_objcoeff(lp_data, index_list, &index_list[n-1], objcoeff);
   //lp_data->si->writeLp("fp.lp");
   termstatus = dual_simplex(lp_data, &iterd);

   if (termstatus != LP_OPTIMAL) {
      PRINT(verbosity,0,("Feasibility Pump: Unable to solve LP. Pump malfunction.\n"));
      return FUNCTION_TERMINATED_ABNORMALLY;
   }

   get_x(lp_data);

   delta_x = 0;
   for (i=0;i<fp_data->n0;i++) {
      x_lp[i]=lp_data->x[i];
      if (fp_data->fp_vars[i]->is_int) {
         delta_x = delta_x+fabs(x_lp[i]-x_ip[i]);
      }
   }
   PRINT(verbosity, 15, ("fp: delta_x = %f\n",delta_x));

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
int fp_round(FPdata *fp_data, LPdata *lp_data)
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
   double **x_bar_val = fp_data->x_bar_val;
   int **x_bar_ind = fp_data->x_bar_ind;
   int *x_bar_len = fp_data->x_bar_len;
   double flip_fraction = fp_data->flip_fraction;
   FPvars **vars = fp_data->fp_vars;

   for (i=0;i<n;i++) {
      if (vars[i]->is_int) {
         /* round x_lp[i] and put into x_ip[i] */
         x_ip[i]=floor(x_lp[i]+0.5);
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
         if (x_ip[i] > lpetol || x_ip[i] < -lpetol){
            tind[cnt] = index[i];
            tx[cnt++] = x_ip[i];
         }
      }
      /* order indices and values according to indices */
      qsort_id(tind, tx, cnt);

      /* go through all 'iter' points and check if x_ip already exists */
      has_changed = TRUE;
      for (i=0; i<fp_data->iter; i++) {
         if (x_bar_len[i] == cnt && fp_data->alpha_p[i] < 0.08) {
            for (j=0; j<cnt; j++) {
               if (tind[j]!=x_bar_ind[i][j] || 
                     fabs(tx[j]-x_bar_val[i][j])>lpetol) {
                  break;
               }
            }
            if (j==cnt) {
               PRINT(fp_data->verbosity,5,("fp: same as %d\n",i));
               break; //its same
            }
         }
      }
      if (i<fp_data->iter) {
         /* flip some vars in x_ip */
         int num_flipped = 0;
         has_changed = FALSE;
         PRINT(fp_data->verbosity,5,("fp: flipping\n"));
         for (j=0; j<n; j++) {
            if (CoinDrand48()<flip_fraction) {
               if (vars[j]->is_bin) {
                  x_ip[j] = 1-x_ip[j];
                  num_flipped++;
               } else if (vars[j]->is_int) {
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

   if (has_changed==TRUE || fp_data->alpha>0) {
      fp_data->x_bar_ind[fp_data->iter] = (int *)malloc(ISIZE*cnt);
      fp_data->x_bar_val[fp_data->iter] = (double *)malloc(DSIZE*cnt);
      x_bar_len[fp_data->iter] = cnt;
      memcpy(fp_data->x_bar_ind[fp_data->iter],tind,ISIZE*cnt);
      memcpy(fp_data->x_bar_val[fp_data->iter],tx,DSIZE*cnt);
      fp_data->alpha = fp_data->alpha*fp_data->alpha_decr;
      if (fp_data->alpha<0.08) {
         fp_data->alpha = 0;
      }
      fp_data->alpha_p[fp_data->iter] = fp_data->alpha;
   } else {
      x_bar_len[fp_data->iter] = -1;
   }
   return termcode;
}

/*===========================================================================*/
int fp_should_call_fp(lp_prob *p, int branching, int *should_call, 
      char is_last_iter)
{
   int        termcode = FUNCTION_TERMINATED_NORMALLY;

   *should_call = FALSE;
   if (is_last_iter==FALSE) {
      return termcode;
   }
   if (p->par.fp_enabled>0 && !branching) {
      if (p->par.fp_enabled == SYM_FEAS_PUMP_REPEATED && 
            p->bc_index%p->par.fp_frequency==0) {
         *should_call = TRUE;
      } else if (p->has_ub==FALSE && p->par.fp_enabled==SYM_FEAS_PUMP_TILL_SOL
            && p->bc_index%p->par.fp_frequency==0) {
         *should_call = TRUE;
      } else if ( (p->has_ub==FALSE||
                   (p->ub-p->lp_data->objval)/(fabs(p->ub)+0.0001)*100>
                   p->par.fp_min_gap) &&
                 p->comp_times.fp < p->par.fp_max_total_time &&
                 p->comp_times.fp < 0.5*p->tt &&
                 p->bc_index%p->par.fp_frequency == 0) {
         *should_call = TRUE;
      }
   }
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

