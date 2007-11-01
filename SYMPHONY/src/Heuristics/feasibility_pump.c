/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "sym_proccomm.h"
#include "qsortucb.h"
#include "sym_lp.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_pack_cut.h"
#include "sym_rounding.h"
#include "symphony.h"
#include "sym_sp.h"
#include "sym_lp_solver.h"
#include "sym_feasibility_pump.h"

/*===========================================================================*/
/*===========================================================================*\
 * This file contains heuristics to find an integral solution after an LP
 * is solved.
 \*===========================================================================*/
/*===========================================================================*/

int feasibility_pump (lp_prob *p, double &solution_value, 
      double *betterSolution)
{
   int termcode = IP_INFEASIBLE;
   FPdata *fp_data = (FPdata*) malloc(sizeof(FPdata));
   LPdata *lp_data = p->lp_data;
   LPdata *new_lp_data = (LPdata *)malloc(sizeof(LPdata));
   int MaxIter = p->par.fp_max_cycles;		/* no. of max pumping cycles */
   int n = lp_data->n;
   int m = lp_data->m;
   int nz = lp_data->nz;
   var_desc **vars = lp_data->vars;
   double pump_start_time = wall_clock(NULL);
   /* row ordered matrices */
   double *Rmatval = (double *) calloc(nz,DSIZE);
   int *Rmatind = (int *) calloc(nz,ISIZE);
   int *Rmatbeg = (int *) calloc(m+1,ISIZE);
   int *Rlength = (int *) calloc(m,ISIZE);
   double *Rupper  = (double *)calloc(m,DSIZE);
   double *Rlower  = (double *)calloc(m,DSIZE);
   /* use OSI to get lp data */
   OsiSolverInterface * model = p->lp_data->si;
   const CoinPackedMatrix * matrix = model->getMatrixByRow();
   const int    *lpRowLengths    = matrix->getVectorLengths();
   const int    *lpRowStarts     = matrix->getVectorStarts();
   const int    *lpRowIndices    = matrix->getIndices();
   const double *lpRowValues     = matrix->getElements();
   const double *lpRowLower      = model->getRowLower();
   const double *lpRowUpper      = model->getRowUpper();
   int count,i,iter,cnt;
   int verbosity = p->par.verbosity;
   int *indexList;
   int *indices;
   double *values;
   int flip_rand=FALSE;
   /* x_* has only n original values. */
   double * x_lp = (double *) calloc(n,DSIZE); /* solution of an lp */
   double * x_ip = (double *) malloc(n*DSIZE); /* rounding of an lp sol */
   double * x_temp = (double *) malloc(n*DSIZE);
   double elapsed_time, real_obj_value;


   PRINT(verbosity,5,("Entering Feasibility Pump.\n"));
   memcpy (Rlength, lpRowLengths, ISIZE*m);
   memcpy (Rupper, lpRowUpper, DSIZE*m);
   memcpy (Rlower, lpRowLower, DSIZE*m);

   count = 0;
   for (i=0;i<m;i++) {
      Rmatbeg[i] = count;
      memcpy (&(Rmatind[Rmatbeg[i]]), &(lpRowIndices[lpRowStarts[i]]), lpRowLengths[i]*ISIZE);
      memcpy (&(Rmatval[Rmatbeg[i]]), &(lpRowValues[lpRowStarts[i]]), lpRowLengths[i]*DSIZE);
      count = count+lpRowLengths[i];
   }
   Rmatbeg[m] = count;		/* always easy to forget this line ;-( */

   lp_data->mip->obj = (double *) malloc(lp_data->n*DSIZE);

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
   fp_round(x_lp,x_ip,vars,n);

   /* retrieve the coeffs in objective */
   get_objcoeffs(lp_data); /* get_objcoeffs is available in symphony */


   /* do the following MaxIter times */
   for (iter=0; iter<MaxIter &&
         wall_clock(NULL)-pump_start_time<p->par.fp_time_limit; iter++) {
      PRINT(verbosity,2,("FP: iteration %d\n",iter));
      if (fp_is_feasible(new_lp_data,x_ip,Rmatbeg,Rmatind,Rmatval,Rlower,
               Rupper,fp_data)) { 
         /* we found what we wanted */
         termcode = IP_HEUR_FEASIBLE;
         memcpy(betterSolution, x_ip, n*DSIZE);

         solution_value = 0;
         for (i=0;i<n;i++) {
            solution_value = solution_value+betterSolution[i]*lp_data->mip->obj[i];
         }
         indices = p->lp_data->tmp.i1; /* n */
         values = p->lp_data->tmp.d; /* n */
         cnt = collect_nonzeros(p, betterSolution, indices, values);
         sp_add_solution(p,cnt,indices,values,solution_value,p->bc_index);
         break;
      }
      else {
         /* solve an lp */
         if (fp_solve_lp(x_lp,x_ip,flip_rand,p->par.fp_flip_fraction,new_lp_data,indexList, fp_data) != 0) {
            break;
         }
      }
      memcpy(x_temp,x_ip,n*DSIZE);
      fp_round(x_lp,x_ip,vars,n);
      flip_rand = TRUE;
      //      printf("DIFFS: \n\n");
      for (i=0;i<n;i++) {
         //printf ("%d\t%f\t%f\n",i,x_lp[i],x_ip[i]);
         if (vars[i]->is_int && x_temp[i]!=x_ip[i]) {
            flip_rand = FALSE;
         }
      }
   }
   close_lp_solver(new_lp_data);
   /* free all the allocated memory */
   FREE(x_lp);
   FREE(x_ip);
   FREE(x_temp);
   FREE(indexList);
   FREE(x_temp);
   FREE(Rmatbeg);
   FREE(Rmatval);
   FREE(Rmatind);
   FREE(Rlower);
   FREE(Rupper);
   FREE(Rlength);
   FREE(lp_data->mip->obj);
   FREE(new_lp_data->x);
   FREE(new_lp_data->lb);
   FREE(new_lp_data->ub);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->dualsol);
   FREE(new_lp_data->dj);
   FREE(new_lp_data->slacks);
   FREE(new_lp_data->mip->is_int);
   FREE(new_lp_data->mip->matbeg);
   FREE(new_lp_data->mip->matind);
   FREE(new_lp_data->mip->matval);
   FREE(new_lp_data->mip->obj);
   FREE(new_lp_data->mip->rhs);
   FREE(new_lp_data->mip->rngval);
   FREE(new_lp_data->mip->sense);
   FREE(new_lp_data->mip->lb);
   FREE(new_lp_data->mip->ub);
   FREE(new_lp_data->mip);
   FREE(new_lp_data);
   for (i=0;i<n;i++) {
      FREE(fp_data->fp_vars[i]);
   }
   FREE(fp_data->fp_vars);
   FREE(fp_data);

   //   printf("Feasibility Pump: exiting\n");

   /* update stats */
   p->tm->stat.fp_time += wall_clock(NULL)-pump_start_time;
   p->tm->stat.fp_calls++;
   if (termcode==IP_HEUR_FEASIBLE) {
      p->tm->stat.fp_successes++;
      if (p->mip->obj_sense == SYM_MAXIMIZE){
         real_obj_value=-solution_value+p->mip->obj_offset;
      } else {
         real_obj_value=solution_value+p->mip->obj_offset;
      }
      elapsed_time = wall_clock(NULL)-p->tm->start_time;
      PRINT(verbosity,-1,("fp: found solution = %10.2f time = %10.2f\n",
               real_obj_value,elapsed_time));
   }

   PRINT(verbosity,5,("Leaving Feasibility Pump.\n"));
   //   exit(0);
   return termcode;
}


/*===========================================================================*/
int fp_round (double *x_lp, double *x_ip, var_desc **vars, const int n)
{
   int i;
   /* rounds x_lp and returns value as x_ip */
   for (i=0;i<n;i++) {
      if (vars[i]->is_int) {
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
int fp_is_feasible (LPdata *lp_data, double *x, int *Rmatbeg,
      int *Rmatind, double *Rmatval, double *Rlower,
      double *Rupper,FPdata *fp_data )
{
   /* check if x is a integer feasible solution to problem in p */
   double lpetol = lp_data->lpetol;
   int n = fp_data->n0;
   int m = fp_data->m0;
   var_desc **vars = lp_data->vars;
   int i,c,j;
   double Ractivity;

   for (i=0;i<n;i++) {
      if (vars[i]->is_int) {
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
      for (j=Rmatbeg[i];j<Rmatbeg[i+1];j++) {
         c=Rmatind[j];
         Ractivity = Ractivity+x[c]*Rmatval[j];
      }
      //      printf("Ractivity[%d] = \t%f\n",i,Ractivity);
      if (Ractivity>Rupper[i]+lpetol || Ractivity<Rlower[i]-lpetol) {
         /* constraint infeasibility is possible since we call this func. after
            rounding */
         return FALSE;
      }
   }

   return TRUE;
}

/*===========================================================================*/
int fp_initialize_lp_solver(LPdata *lp_data, LPdata *newdata, FPdata *fp_data)
{
   /*
      create a copy of lp_data into newdata
      for general mixed int programs, we will have to add 2 new vars for each
      non-binary integer var. (x_j+ and x_j-)
      */


   /* first create an exact copy of lp_data */
   newdata->lpetol = lp_data->lpetol;
   int n = lp_data->n;
   int m = lp_data->m;
   newdata->n = n;
   newdata->m = m;
   newdata->nz = lp_data->nz;
   newdata->maxn = lp_data->maxn;
   newdata->maxm = lp_data->maxm;
   newdata->maxnz = lp_data->maxnz;
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
   //   newdata->mip = lp_data->mip;	/* mip currently empty, we also need
   //			   fp_get_mip_desc */ 

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
         if (lp_data->vars[i]->lb!=0||lp_data->vars[i]->ub!=1) {
            fp_vars[i]->isBin = FALSE;
            fp_data->numNonBinInts++;
         }
         else {
            fp_vars[i]->isBin = TRUE;
         }
      }
   }

   newn = n+2*fp_data->numNonBinInts;
   newm = m+fp_data->numNonBinInts;
   newdata->vars = lp_data->vars;
   newdata->rows = lp_data->rows;
   newdata->dualsol = (double*)malloc(newm*DSIZE);
   newdata->dj = (double*)malloc(newn*DSIZE);
   newdata->slacks = (double*)malloc(newm*DSIZE);
   newdata->x = (double *) malloc(newn*DSIZE);
   newdata->mip = (MIPdesc *) calloc(sizeof(MIPdesc),1);
   newdata->lb = (double*) malloc(newn*DSIZE);
   newdata->ub = (double*) malloc(newn*DSIZE);

   /* copy mipdesc of lp_data into that of newdata */
   fp_get_mip_desc(lp_data, newdata);

   open_lp_solver(newdata);

   /* the second parameter is redundant, third only used for cplex: FASTMIP,
      whatever that is */
   load_lp_prob(newdata, 666, 0);

   rstat = (int *) malloc(m * ISIZE);
   cstat = (int *) malloc(n * ISIZE);

   get_basis(lp_data,cstat,rstat);
   load_basis (newdata,cstat,rstat);

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
   memcpy(newdata->lb,lp_data->lb,lp_data->n*DSIZE);
   memcpy(newdata->ub,lp_data->ub,lp_data->n*DSIZE);

   for (i=0;i<n;i++) {
      if (lp_data->vars[i]->is_int && !fp_vars[i]->isBin) {
         /*
            obj coeff. is 1, lb is zero
            constr is:
            x_i - x_i+ + x_i- = x_i(bar).
            first add 2 empty cols.
            */
         lb = 0;
         /* add x_j+ */
         ub = (fabs(lp_data->x[i]-lp_data->ub[i])<newdata->lpetol)?0:lp_data->vars[i]->ub;
         add_cols(newdata, 1, 0, &one, cmatbeg, cmatind, cmatval, &lb, &ub, &where_to_move);
         fp_vars[i]->xplus = colNumber;
         colNumber++;

         /* add x_j- */
         ub = (fabs(lp_data->x[i]-lp_data->lb[i])<newdata->lpetol)?0:lp_data->vars[i]->ub;
         add_cols(newdata, 1, 0, &one, cmatbeg, cmatind, cmatval, &lb, &ub, &where_to_move);
         fp_vars[i]->xminus = colNumber;
         colNumber++;

         /* now add a row */
         rhs = lp_data->x[i];
         rmatind[0] = i;
         rmatind[1] = fp_vars[i]->xplus;
         rmatind[2] = fp_vars[i]->xminus;
         add_rows(newdata, 1, 3, &rhs, &sense, rmatbeg, rmatind, rmatval);
      }
   }
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
   double *coeffList= (double *) calloc(lp_data->n,DSIZE);
   int n = lp_data->n;
   int iterd;
   int termstatus;
   int i;
   double delta_x;

   for (i=0;i<fp_data->n0;i++) {
      if (lp_data->vars[i]->is_int) {
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
               if (x_ip[i]==lp_data->vars[i]->ub) {
                  x_ip[i]=x_ip[i]--;
               }
               else if (x_ip[i]==lp_data->vars[i]->lb) {
                  x_ip[i]=x_ip[i]++;
               }
               else if (CoinDrand48()<0.5) {
                  x_ip[i]=x_ip[i]--;
               }
               else {
                  x_ip[i]=x_ip[i]++;
               }
            }
            if (x_ip[i]!=lp_data->vars[i]->ub) {
               coeffList[fp_data->fp_vars[i]->xplus] = 1;
            }
            if (x_ip[i]!=lp_data->vars[i]->lb) {
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
      return 1;
   }

   FREE(lp_data->x);
   lp_data->x = (double *) malloc(lp_data->n*DSIZE);
   get_x(lp_data);

   delta_x = 0;
   for (i=0;i<fp_data->n0;i++) {
      x_lp[i]=lp_data->x[i];
      if (lp_data->vars[i]->is_int) {
         delta_x = delta_x+fabs(x_lp[i]-x_ip[i]);
      }
   }
   //   printf("delta_x = %f\n",delta_x);

   FREE(coeffList);

   return 0;
}

/*===========================================================================*/
int fp_get_mip_desc(LPdata *lp_data, LPdata *newdata)
{
   /* get mip desc for lp_data using get_*() functions and build mip for
      newdata */
   MIPdesc *newmip = newdata->mip;
   MIPdesc *oldmip = lp_data->mip;

   newmip->n = newdata->n;
   newmip->m = newdata->m;
   newmip->nz = newdata->nz;

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
int fp_should_call_fp(lp_prob *p, int branching)
{
   if (p->par.fp_enabled>0 && 
       !branching &&
       p->bc_index%p->par.fp_frequency == 1 &&
       p->has_ub==FALSE  &&
       p->tm->stat.fp_time < p->par.fp_max_total_time) {
      return TRUE;
   } else {
      return FALSE;
   }

   return FALSE;
}
/*===========================================================================*/
/*===========================================================================*/
