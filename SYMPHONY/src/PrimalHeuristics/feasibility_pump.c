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

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sym_proccomm.h"
#include "sym_qsort.h"
#include "sym_lp.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_pack_cut.h"
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
 * TODO: termcode
 * make independent of solver
 * remove newlpdata
 * Total solutions found. add a statistic to sp
 * change rhs, lb and ub for gen int
 * change isBin to is_bin, isInt to is_int
 */

int feasibility_pump (lp_prob *p, char *found_better_solution,
      double &solution_value, double *betterSolution)
{
   int                      termcode    = FUNCTION_TERMINATED_NORMALLY;
   FPdata                  *fp_data     = (FPdata*) malloc(sizeof(FPdata));
   LPdata                  *lp_data     = p->lp_data;
   LPdata                  *new_lp_data = (LPdata *)calloc(1,sizeof(LPdata));	
   /* no. of max pumping cycles */
   int                      MaxIter     = p->par.fp_max_cycles;
   int                      n           = lp_data->n;
   /* use OSI to get lp data */
   OsiSolverInterface      *model       = p->lp_data->si;
   const CoinPackedMatrix  *matrix      = model->getMatrixByRow();
   const double            *lp_r_low    = model->getRowLower();
   const double            *lp_r_up     = model->getRowUpper();
   int                      count, i, r, iter, cnt, verbosity;
   int                     *indices;
   double                  *values;
   int                      flip_rand = FALSE;
   double                   fp_time, real_obj_value, target_ub;
   FPvars                 **vars;
   double                   gap           = model->getInfinity();
   double                   obj_lb        = lp_data->objval;
   double                   total_time    = 0;
   const double            *mip_obj       = model->getObjCoefficients();
   char                     is_feasible   = FALSE;
   double                  *x_ip,*x_lp;

   fp_time                                = used_time(&total_time);
   /* total_time and fp_time both now have total time used by symphony's lp
    * process */
   fp_time                                = used_time(&total_time);
   /* fp_time should now be zero and total_time be still the same */

   *found_better_solution = FALSE;
   verbosity = fp_data->verbosity     = 111; //p->par.verbosity;
   fp_data->mip_obj       = (double *)malloc(n*DSIZE);
   fp_data->flip_fraction = p->par.fp_flip_fraction;
   memcpy(fp_data->mip_obj,mip_obj,n*DSIZE);

   /* initialize the lp solver. load the current basis */
   fp_initialize_lp_solver(p, new_lp_data, fp_data);
   x_ip = fp_data->x_ip;
   x_lp = fp_data->x_lp;

   /* round the x_lp and store as x_ip, it will usually become infeasible */
   vars = fp_data->fp_vars;

   /* do the following MaxIter times */
   fp_time += used_time(&total_time);
   for (iter=0; iter<MaxIter && fp_time<p->par.fp_time_limit; iter++) {
      PRINT(verbosity,5,("fp: iteration %d\n",iter));
      is_feasible = FALSE;
      /* solve an lp */
      fp_round(fp_data, new_lp_data);
      fp_is_feasible (lp_data, matrix, lp_r_low, lp_r_up, fp_data, &is_feasible);

      if (is_feasible == TRUE) {
         /* we found what we wanted */
         memcpy(betterSolution, x_ip, n*DSIZE);

         solution_value = 0;
         for (i=0;i<n;i++) {
            solution_value = solution_value + betterSolution[i]*mip_obj[i];
         }
         indices = p->lp_data->tmp.i1;          /* n */
         values  = p->lp_data->tmp.d;           /* n */
         cnt     = collect_nonzeros(p, betterSolution, indices, values);
         gap     = (solution_value - obj_lb)/(fabs(solution_value)+0.001)*100;
         p->lp_stat.fp_num_sols++;
         PRINT(verbosity,5,("fp: found solution with value = %f\n",
                  solution_value));
         PRINT(verbosity,5,("fp: gap = %f\n", gap));
         sp_add_solution(p,cnt,indices,values,solution_value,p->bc_index);
         if (gap <= p->par.fp_min_gap) {
            *found_better_solution = TRUE;
            break;
         }
         target_ub = (obj_lb + solution_value)/2;
         if (*found_better_solution != TRUE) {
            // add another objective function constraint to lower the
            // objective value.
            fp_add_obj_row(new_lp_data, n, mip_obj, target_ub);
            *found_better_solution = TRUE;
         } else {
            r = new_lp_data->m-1;
            change_rhs(new_lp_data, 1, &r, &target_ub);
         }
      } 

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
   p->comp_times.primal_heur      += fp_time;
   p->lp_stat.fp_calls++;
   if (*found_better_solution==TRUE) {
      if (p->mip->obj_sense == SYM_MAXIMIZE){
         real_obj_value=-solution_value+p->mip->obj_offset;
      } else {
         real_obj_value=solution_value+p->mip->obj_offset;
      }
      PRINT(verbosity,-1,("fp: found solution = %10.2f time = %10.2f\n",
               real_obj_value,total_time));
   }

   PRINT(verbosity,5,("Leaving Feasibility Pump.\n"));
    //exit(0);
   return termcode;
}


/*===========================================================================*/
int fp_round (FPdata *fp_data, LPdata *lp_data)
{
   int      i;
   int      n    = fp_data->n;
   FPvars **vars = fp_data->fp_vars;
   double  *x_ip = fp_data->x_ip;
   double  *x_lp = fp_data->x_lp;

   /* rounds x_lp and returns value as x_ip */
   /* add x_ip to list of solutions */
   fp_add_rounded_point(fp_data, lp_data);
   return 0;
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
      if (vars[i]->isInt) {
         if (x[i]-floor(x[i])>lpetol && ceil(x[i])-x[i]>lpetol) {
            /* some int variable is non-integral */
            /* is not possible, since this function is called after rounding */
            printf("Bok!\n");
            return termcode;
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
int fp_initialize_lp_solver(lp_prob *p, LPdata *new_data, FPdata *fp_data)
{
   /*
      create a copy of lp_data into new_data
      for general mixed int programs, we will have to add 2 new vars for each
      non-binary integer var. (x_j+ and x_j-)
      */

   /* first create an exact copy of lp_data */
   LPdata *lp_data  = p->lp_data;
   new_data->lpetol = lp_data->lpetol;
   int n = lp_data->n;
   int m = lp_data->m;
   int i,newn, newm;
   int *rstat,*cstat;

   double one=1.0;
   char sense='E';
   char where_to_move='E';	/* redundant */
   int colNumber = n;
   int *rmatbeg = (int *) malloc(2*ISIZE);
   int *cmatbeg = (int *) malloc(2*ISIZE);
   int *rmatind = (int *) malloc(3*ISIZE);
   double *rmatval = (double *) malloc(3*DSIZE);
   int *cmatind;
   double *cmatval;
   double rhs;
   double lb, ub;
   double lpetol = lp_data->lpetol;
   double *x = lp_data->x;
   double *lp_lb;
   double *lp_ub;
   double norm_c = 0;
   double *mip_obj = fp_data->mip_obj;
   int verbosity = fp_data->verbosity;
   int *index_list;

   /* used because we can not call si directly */
   copy_lp_data(lp_data,new_data);
   lp_lb = new_data->lb;
   lp_ub = new_data->ub;

   /* set up fp_data */
   fp_data->alpha         = 0.8;
   fp_data->alpha_decr    = 0.7;
   fp_data->n0            = n;
   fp_data->m0            = m;
   fp_data->iter          = 0;

   /* count how many binary variables */
   fp_data->fp_vars       = (FPvars **) malloc(sizeof(FPvars *)*n);
   fp_data->x_ip          = (double *) calloc(n,DSIZE);
   fp_data->x_lp          = (double *) calloc(n,DSIZE);
   fp_data->index_list    = (int *)    calloc(n,DSIZE);
   fp_data->x_bar_ind     = (int **)   calloc(p->par.fp_max_cycles,
                                              sizeof(int*));
   fp_data->x_bar_val     = (double **)calloc(p->par.fp_max_cycles,
                                              sizeof(double*));
   fp_data->x_bar_len     = (int *)    calloc(p->par.fp_max_cycles,ISIZE);
   fp_data->alpha_p       = (double *) malloc(p->par.fp_max_cycles*DSIZE);
   FPvars **fp_vars       = fp_data->fp_vars;
   fp_data->numNonBinInts = 0;
   fp_data->numInts       = 0;

   // TODO: make this work when new cols are added
   memcpy(fp_data->x_lp,p->lp_data->x,DSIZE*n);
   
   index_list = fp_data->index_list;
   for (i=0;i<n;i++) {
      index_list[i]=i;
      fp_vars[i] = (FPvars *)malloc(sizeof(FPvars));
      if (lp_data->vars[i]->is_int) {
         fp_data->numInts++;
         fp_vars[i]->isInt = TRUE;
         if (lp_data->vars[i]->lb<-lpetol||lp_data->vars[i]->ub>1+lpetol) {
            fp_vars[i]->isBin = FALSE;
            fp_data->numNonBinInts++;
         }
         else {
            fp_vars[i]->isBin = TRUE;
         }
      } else {
         fp_vars[i]->isInt = fp_vars[i]->isBin = FALSE;
      }
      /* calculate ||C|| */
      norm_c += mip_obj[i]*mip_obj[i];
   }
   
   norm_c = sqrt(norm_c);
   PRINT(verbosity, 20, ("fp: norm_c = %f\n",norm_c));

   fp_data->n       = n+2*fp_data->numNonBinInts;
   fp_data->m       = m;
   fp_data->obj     = (double *)malloc(n*DSIZE);

   if (norm_c>lpetol) {
      for (i=0;i<n;i++) {
         mip_obj[i] = mip_obj[i]/norm_c;
      }
   }
   
   /* load basis */
   rstat = (int *) malloc(m * ISIZE);
   cstat = (int *) malloc(n * ISIZE);

   get_basis(lp_data,cstat,rstat);
   load_basis (new_data,cstat,rstat);

   FREE(rstat);
   FREE(cstat);

#if 0

   /* add 2 columns for each nonBinary Integer */
   rmatbeg[0] = 0;
   rmatbeg[1] = 3;
   cmatbeg[0] = 0;
   cmatbeg[1] = 0;
   rmatval[0] = 1;
   rmatval[1] = -1;
   rmatval[2] = 1;

   /* copy the ub and lb of the original vars */
   /*
   get_bounds(lp_data);   
   memcpy(new_data->lb,lp_data->lb,lp_data->n*DSIZE);
   memcpy(new_data->ub,lp_data->ub,lp_data->n*DSIZE);
   */

   //get_objcoeffs(new_data); /* get_objcoeffs is available in symphony */
   //TODO: this is not entirely correct
   for (i=0;i<n;i++) {
      if (fp_vars[i]->isInt && !fp_vars[i]->isBin) {
         /*
            obj coeff. is 1, lb is zero
            constr is:
            x_i - x_i+ + x_i- = x_i(bar).
            first add 2 empty cols.
         */
         lb = 0;
         /* add x_j+ */
         ub = (fabs(x[i]-lp_ub[i])<lpetol) ? 0 : lp_ub[i];
         add_cols(new_data, 1, 0, &one, cmatbeg, cmatind, cmatval, &lb, &ub, 
               &where_to_move);
         fp_vars[i]->xplus = colNumber;
         colNumber++;

         /* add x_j- */
         ub = (fabs(x[i]-lp_lb[i])<lpetol) ? 0 : lp_ub[i];
         add_cols(new_data, 1, 0, &one, cmatbeg, cmatind, cmatval, &lb, &ub, 
               &where_to_move);
         fp_vars[i]->xminus = colNumber;
         colNumber++;

         /* now add a row */
         rhs = lp_data->x[i];
         rmatind[0] = i;
         rmatind[1] = fp_vars[i]->xplus;
         rmatind[2] = fp_vars[i]->xminus;
         add_rows(new_data, 1, 3, &rhs, &sense, rmatbeg, rmatind, rmatval);
      }
   }
#endif
   /* used by change_rhs */
   new_data->tmp.c = (char *)malloc(2*CSIZE);
   new_data->tmp.d = (double *)malloc(DSIZE*n);
   new_data->tmp.i1 = (int *)malloc(ISIZE*n);

   FREE(rmatval);
   FREE(rmatind);
   FREE(cmatbeg);
   FREE(rmatbeg);

   return 0;
}

/*===========================================================================*/
int fp_solve_lp(LPdata *lp_data, FPdata *fp_data, char* is_feasible) 
{
   /* construct an lp based on x_ip. solve it. store the result in x_lp */
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
   double flip_fraction = fp_data->flip_fraction;
   char flip_rand = FALSE;
   int  *index_list = fp_data->index_list;
   double *x_ip = fp_data->x_ip;
   double *x_lp = fp_data->x_lp;
   double alpha = fp_data->alpha;

   is_feasible = FALSE;
   memset ((char *)(objcoeff),0,DSIZE*n);
   for (i=0;i<fp_data->n0;i++) {
      if (fp_vars[i]->isInt) {
         if (fp_vars[i]->isBin) {
            if (flip_rand && CoinDrand48()<flip_fraction) {
               x_ip[i]=1-x_ip[i];
            }
            if (x_ip[i]==0) {
               objcoeff[i] = 1.0;
            }
            else if (x_ip[i]==1) {
               objcoeff[i] = -1.0;
            }
         }
         else {
            if (flip_rand && CoinDrand48()<flip_fraction && lp_data->ub[i]>=lp_data->lb[i]+1) {
               if (x_ip[i]==lp_data->ub[i]) {
                  x_ip[i]=x_ip[i]--;
               }
               else if (x_ip[i]==lp_data->lb[i]) {
                  x_ip[i]=x_ip[i]++;
               }
               else if (CoinDrand48()<0.5) {
                  x_ip[i]=x_ip[i]--;
               }
               else {
                  x_ip[i]=x_ip[i]++;
               }
            }
            if (x_ip[i]!=lp_data->ub[i]) {
               objcoeff[fp_vars[i]->xplus] = 1;
            }
            if (x_ip[i]!=lp_data->lb[i]) {
               objcoeff[fp_vars[i]->xminus] = 1;
            }
         }
      }
      /* calculate ||coeff||, norm is not zero because otherwise x_ip is
       * feasible */
      norm += objcoeff[i]*objcoeff[i];
   }

   norm = sqrt(norm);
   //norm = 0;
   PRINT(verbosity, 15, ("fp: norm = %f\n",norm));
   for (i=0;i<fp_data->n0;i++) {
      objcoeff[i] = (1-alpha)*objcoeff[i]+alpha*mip_obj[i]*norm;
   }
   alpha = alpha*fp_data->alpha_decr;

   change_objcoeff(lp_data, index_list, &index_list[n-1], objcoeff);
   termstatus = dual_simplex(lp_data, &iterd);
   if (termstatus != LP_OPTIMAL) {
      PRINT(verbosity,0,("Feasibility Pump: Unable to solve LP. Pump malfunction.\n"));
      return FUNCTION_TERMINATED_ABNORMALLY;
   }

   get_x(lp_data);

   delta_x = 0;
   for (i=0;i<fp_data->n0;i++) {
      x_lp[i]=lp_data->x[i];
      if (fp_data->fp_vars[i]->isInt) {
         delta_x = delta_x+fabs(x_lp[i]-x_ip[i]);
      }
   }
   PRINT(verbosity, 15, ("fp: delta_x = %f\n",delta_x));

   return 0;
}

/*===========================================================================*/
int fp_get_mip_desc(LPdata *lp_data, LPdata *new_data)
{
   /* get mip desc for lp_data using get_*() functions and build mip for
      new_data */
   MIPdesc *newmip = new_data->mip;
   MIPdesc *oldmip = lp_data->mip;

   newmip->n = new_data->n;
   newmip->m = new_data->m;
   newmip->nz = new_data->nz;

   newmip->is_int = (char *) malloc(newmip->n*CSIZE);
   newmip->matbeg = (int *)  malloc((newmip->n+1)*ISIZE);
   newmip->matind = (int *)  malloc(newmip->nz*ISIZE);
   newmip->matval = (double *) malloc(newmip->nz*DSIZE);
   newmip->obj    = (double *) malloc(newmip->n*DSIZE);
   newmip->rhs    = (double *) malloc(newmip->m*DSIZE);
   newmip->rngval = (double *) calloc(newmip->m,DSIZE);
   newmip->sense  = (char *) malloc(newmip->m*CSIZE);
   newmip->lb     = (double *) malloc(newmip->n*DSIZE);
   newmip->ub     = (double *) malloc(newmip->n*DSIZE);
   //TODO: move these up
   int *colLength = (int *) malloc(newmip->n*ISIZE);
   int i;

   /* get data by querying the lp-solver */
   get_objcoeffs(lp_data);
   memcpy(newmip->obj, oldmip->obj, newmip->n*DSIZE);

   newmip->matbeg[0] = 0;
   for (i=0; i<newmip->n; i++) {
      get_column(lp_data, i, &newmip->matval[newmip->matbeg[i]], &newmip->matind[newmip->matbeg[i]], &colLength[i], &newmip->obj[i]);     
      newmip->matbeg[i+1] = newmip->matbeg[i] + colLength[i];
   }

   lp_data->mip->rhs = (double *) malloc(lp_data->m*DSIZE);
   lp_data->mip->rngval = (double *) malloc(lp_data->m*DSIZE);
   lp_data->mip->sense = (char *) malloc(lp_data->m*CSIZE);
   get_rhs_rng_sense(lp_data);
   memcpy(newmip->rhs, oldmip->rhs, lp_data->m*DSIZE);
   memcpy(newmip->rngval, oldmip->rngval, lp_data->m*DSIZE);
   memcpy(newmip->sense, oldmip->sense, lp_data->m*CSIZE);

   get_bounds(lp_data);   
   memcpy(newmip->lb,lp_data->lb,newmip->n*DSIZE);
   memcpy(newmip->ub,lp_data->ub,newmip->n*DSIZE);

   FREE(colLength);
   FREE(lp_data->mip->rhs);
   FREE(lp_data->mip->rngval);
   FREE(lp_data->mip->sense);   

   return 0;
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
int fp_should_call_fp(lp_prob *p, int branching)
{
   var_desc **vars = p->lp_data->vars;
   int i;

   if (p->par.fp_enabled>0 && 
       !branching &&
       p->bc_index%p->par.fp_frequency == 0 &&
       p->has_ub==FALSE  &&
       p->comp_times.fp < p->par.fp_max_total_time) {
      /*
       * check if it has ints
       * TODO: remove this check. make fp work for general ints
       */
      for (i=0;i<p->lp_data->n;i++) {
         if (vars[i]->is_int && (vars[i]->lb<0 || vars[i]->ub>1)) {
            p->par.fp_enabled = 0;
            return FALSE;
         }
      }
      return TRUE;
   } else {
      return FALSE;
   }

   return FALSE;
}

/*===========================================================================*/
int fp_add_rounded_point(FPdata *fp_data, LPdata *lp_data)
{
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   double *x_ip = fp_data->x_ip;
   double *x_lp = fp_data->x_lp;
   int i,j, has_changed;
   int n = fp_data->n;
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
      if (vars[i]->isInt) {
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
      for (i=0; i<fp_data->iter; i++) {
         if (fp_data->x_bar_len[i] == cnt && fp_data->alpha_p[i] < 0.08) {
            for (j=0; j<cnt; j++) {
               if (tind[j]!=x_bar_ind[i][j] || fabs(tx[j]-x_bar_val[i][j])>lpetol) {
                  break;
               }
            }
            if (j==cnt) {
               PRINT(fp_data->verbosity,15,("fp: same as %d\n",i));
               break; //its same
            }
         }
      }
      if (i<fp_data->iter) {
         /* flip some vars in x_ip */
         int num_flipped = 0;
         has_changed = FALSE;
         PRINT(fp_data->verbosity,15,("fp: flipping\n"));
         for (j=0; j<n; j++) {
            if (CoinDrand48()<flip_fraction && vars[j]->isBin) {
               x_ip[j] = 1-x_ip[j];
               num_flipped++;
            }
         }
         PRINT(fp_data->verbosity,15,("fp: flipping %d\n", num_flipped));
      } else {
         break;
      }
   }

   fp_data->x_bar_ind[fp_data->iter] = (int *)malloc(ISIZE*cnt);
   fp_data->x_bar_val[fp_data->iter] = (double *)malloc(DSIZE*cnt);
   fp_data->x_bar_len[fp_data->iter] = cnt;
   memcpy(fp_data->x_bar_ind[fp_data->iter],tind,ISIZE*cnt);
   memcpy(fp_data->x_bar_val[fp_data->iter],tx,DSIZE*cnt);
   fp_data->alpha = fp_data->alpha*fp_data->alpha_decr;
   if (fp_data->alpha<0.08) {
      fp_data->alpha = 0;
   }
   fp_data->alpha_p[fp_data->iter] = fp_data->alpha;
   return (termcode);
}

/*===========================================================================*/
int fp_should_change_int_point(FPdata *fp_data)
{
   int termcode = FUNCTION_TERMINATED_NORMALLY;
   int i;


   return termcode;
}
/*===========================================================================*/
/*===========================================================================*/

