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
#include <malloc.h>
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
   const double            *mip_obj     = model->getObjCoefficients();
   int                      count,i,r,iter,cnt;
   int                      verbosity = p->par.verbosity;
   int                     *indexList, *indices;
   double                  *values;
   int                      flip_rand = FALSE;
   /* x_* has only n original values. */
   /* solution of an lp */
   double                  *x_lp = (double *) calloc(n,DSIZE);
   /* rounding of an lp sol */
   double                  *x_ip = (double *) malloc(n*DSIZE);
   double                  *x_temp = (double *) malloc(n*DSIZE);
   double                   fp_time, real_obj_value, target_ub;
   FPvars                 **vars;
   int                      min_verbosity = 5;
   double                   gap           = model->getInfinity();
   double                   obj_lb        = lp_data->objval;
   double                   total_time    = 0;

   fp_time                                = used_time(&total_time);
   /* total_time and fp_time both now have total time used by symphony's lp
    * process */

   fp_time                                = used_time(&total_time);
   /* fp_time should now be zero and total_time be still the same */


   *found_better_solution = FALSE;

   /* initialize the lp solver. load the current basis */
   fp_initialize_lp_solver(lp_data, new_lp_data, fp_data);
   indexList = (int *) malloc(new_lp_data->n*ISIZE);
   for (i=0;i<n;i++) {
      x_lp[i]=lp_data->x[i];
      /* initialize the indexList */
      indexList[i]=i;
   }
   for (i=n;i<new_lp_data->n;i++) {
      indexList[i]=i;
   }

   /* round the x_lp and store as x_ip, it will usually become infeasible */
   vars = fp_data->fp_vars;
   fp_round(x_lp,x_ip,vars,n);

   /* do the following MaxIter times */
   fp_time += used_time(&total_time);
   for (iter=0; iter<MaxIter && fp_time<p->par.fp_time_limit; iter++) {
      PRINT(verbosity,min_verbosity,("fp: iteration %d\n",iter));
      if (fp_is_feasible(new_lp_data, x_ip, matrix, lp_r_low, lp_r_up,
               fp_data)) { 
         /* we found what we wanted */
         memcpy(betterSolution, x_ip, n*DSIZE);

         solution_value = 0;
         for (i=0;i<n;i++) {
            solution_value = solution_value+
               betterSolution[i]*mip_obj[i];
         }
         indices = p->lp_data->tmp.i1;          /* n */
         values  = p->lp_data->tmp.d;           /* n */
         cnt     = collect_nonzeros(p, betterSolution, indices, values);
         gap     = (solution_value - obj_lb)/(fabs(solution_value)+0.001)*100;
         p->lp_stat.fp_num_sols++;
         PRINT(verbosity,min_verbosity,("fp: found solution with value = %f\n",
                  solution_value));
         PRINT(verbosity,min_verbosity,("fp: gap = %f\n", gap));
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

      /* solve an lp */
      if (fp_solve_lp(x_lp,x_ip,flip_rand,p->par.fp_flip_fraction,
               new_lp_data,indexList, fp_data) != FUNCTION_TERMINATED_NORMALLY
         ) {
         break;
      }
      
      memcpy(x_temp,x_ip,n*DSIZE);
      fp_round(x_lp,x_ip,vars,n);
      flip_rand = TRUE;
      //      printf("DIFFS: \n\n");
      for (i=0;i<n;i++) {
         //printf ("%d\t%f\t%f\t",i,x_temp[i],x_ip[i]);
         if (vars[i]->isInt && x_temp[i]!=x_ip[i]) {
            flip_rand = FALSE;
            //printf("different\n");
         } else if (vars[i]->isInt) {
            //printf("same\n");
         } else {
            //printf("\n");
         }
      }
      fp_time += used_time(&total_time);
   }
   close_lp_solver(new_lp_data);
   /* free all the allocated memory */
   FREE(x_lp);
   FREE(x_ip);
   FREE(x_temp);
   FREE(indexList);
   FREE(x_temp);
   FREE(new_lp_data->x);
   FREE(new_lp_data->lb);
   FREE(new_lp_data->ub);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->dualsol);
   FREE(new_lp_data->dj);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->tmp.c);
   FREE(new_lp_data);
   for (i=0;i<n;i++) {
      FREE(fp_data->fp_vars[i]);
   }
   FREE(fp_data->fp_vars);
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

   PRINT(verbosity,min_verbosity,("Leaving Feasibility Pump.\n"));
    //exit(0);
   return termcode;
}


/*===========================================================================*/
int fp_round (double *x_lp, double *x_ip, FPvars **vars, const int n)
{
   int i;
   /* rounds x_lp and returns value as x_ip */
   for (i=0;i<n;i++) {
      if (vars[i]->isInt) {
         /* round x_lp[i] and put into x_ip[i] */
         x_ip[i]=rint(x_lp[i]);
      }
      else {
         x_ip[i]=x_lp[i];
      }
   }
   return 0;
}

/*===========================================================================*/
int fp_is_feasible (LPdata *lp_data, double *x, const CoinPackedMatrix *matrix, const double *r_low, const double *r_up, FPdata *fp_data )
{
   /* check if x is a integer feasible solution to problem in p */
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

   for (i=0;i<n;i++) {
      if (vars[i]->isInt) {
         if (x[i]-floor(x[i])>lpetol && ceil(x[i])-x[i]>lpetol) {
            /* some int variable is non-integral */
            /* is not possible, since this function is called after rounding */
            printf("Bok!\n");
            return FALSE;
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
         return FALSE;
      }
   }

   return TRUE;
}

/*===========================================================================*/
int fp_initialize_lp_solver(LPdata *lp_data, LPdata *new_data, FPdata *fp_data)
{
   /*
      create a copy of lp_data into new_data
      for general mixed int programs, we will have to add 2 new vars for each
      non-binary integer var. (x_j+ and x_j-)
      */

   /* first create an exact copy of lp_data */
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
   
   /* used because we can not call si directly */
   copy_lp_data(lp_data,new_data);
   lp_lb = new_data->lb;
   lp_ub = new_data->ub;

   /* set up fp_data */
   fp_data->n0 = n;
   fp_data->m0 = m;

   /* count how many binary variables */
   fp_data->fp_vars = (FPvars **) malloc(sizeof(FPvars *)*n);
   FPvars **fp_vars = fp_data->fp_vars;
   fp_data->numNonBinInts = 0;
   fp_data->numInts = 0;
   for (i=0;i<n;i++) {
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
   }

   newn = n+2*fp_data->numNonBinInts;
   newm = m+fp_data->numNonBinInts;

   rstat = (int *) malloc(m * ISIZE);
   cstat = (int *) malloc(n * ISIZE);

   get_basis(lp_data,cstat,rstat);
   load_basis (new_data,cstat,rstat);

   FREE(rstat);
   FREE(cstat);

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

   new_data->tmp.c = (char *)malloc(2*CSIZE);

   FREE(rmatval);
   FREE(rmatind);
   FREE(cmatbeg);
   FREE(rmatbeg);

   return 0;
}

/*===========================================================================*/
int fp_solve_lp(double *x_lp, double *x_ip, int flip_rand, double flip_fraction, LPdata *lp_data, int* indexList, FPdata *fp_data) 
{
   /* construct an lp based on x_ip. solve it. store the result in x_lp */

   /* construct new objcoeff */
   //TODO: remove this
   double *coeffList= (double *) calloc(lp_data->n,DSIZE);
   int n = lp_data->n;
   int iterd;
   int termstatus;
   int i;
   double delta_x;

   /* TODO: make this more efficient */

   /*
   if (flip_rand) {
      PRINT(5,0,("flipping\n"));
   }
   */

   for (i=0;i<fp_data->n0;i++) {
      if (fp_data->fp_vars[i]->isInt) {
         if (fp_data->fp_vars[i]->isBin) {
            if (flip_rand && CoinDrand48()<flip_fraction) {
               x_ip[i]=1-x_ip[i];
            }
            if (x_ip[i]==0) {
               coeffList[i] = 1.0;
            }
            else if (x_ip[i]==1) {
               coeffList[i] = -1.0;
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
               coeffList[fp_data->fp_vars[i]->xplus] = 1;
            }
            if (x_ip[i]!=lp_data->lb[i]) {
               coeffList[fp_data->fp_vars[i]->xminus] = 1;
            }
         }
      }
   }

   change_objcoeff(lp_data, indexList, &indexList[n-1], coeffList);
   termstatus = dual_simplex(lp_data, &iterd);
   if (termstatus != LP_OPTIMAL) {
      printf("Feasibility Pump: Unable to solve LP. Pump malfunction.\n");
      FREE(coeffList);
      return FUNCTION_TERMINATED_ABNORMALLY;
   }

   //TODO: remove this
   FREE(lp_data->x);
   lp_data->x = (double *) malloc(lp_data->n*DSIZE);
   get_x(lp_data);

   delta_x = 0;
   for (i=0;i<fp_data->n0;i++) {
      x_lp[i]=lp_data->x[i];
      if (fp_data->fp_vars[i]->isInt) {
         delta_x = delta_x+fabs(x_lp[i]-x_ip[i]);
      }
   }
   //PRINT(5,1,("delta_x = %f\n",delta_x));

   FREE(coeffList);

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
/*===========================================================================*/
