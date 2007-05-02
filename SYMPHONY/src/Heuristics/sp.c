/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include "sym_constants.h"
#include "sym_sp.h"
#include "sym_macros.h"

/* Functions related to solution pool */

/*===========================================================================*/
/*===========================================================================*/
int sp_add_solution (lp_prob *p, int cnt, int *indices, double *values, double obj_value, int bc_index)
{
   sp_desc *sp = p->tm->sp;
   PRINT(p->par.verbosity,-1,("sp: solution pool size = %d. \n", sp->num_solutions));
   
   if (sp->num_solutions == sp->max_solutions) {
      /* delete first solution and move everything up by 1 */
      sp_delete_solution(sp,0);
      for (int i=0;i<(sp->max_solutions-1);i++) {
	 sp->solutions[i] = sp->solutions[i+1];
      }
   }
   sp_solution *sol = sp->solutions[sp->num_solutions];
   sol->objval = obj_value;
   sol->xlength = cnt;
   sol->xind = (int *) malloc(ISIZE*cnt);
   memcpy(sol->xind,indices,ISIZE*cnt);
   sol->xval = (double *) malloc(DSIZE*cnt);
   memcpy(sol->xval,values,DSIZE*cnt);
   sol->node_index = bc_index;
   sp->num_solutions++;
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_delete_solution (sp_desc *sp, int position)
{
   /* not implemented yet */
   FREE(sp->solutions[position]->xind);
   FREE(sp->solutions[position]->xval);
   for (int i=position; i<sp->max_solutions-1; i++) {
      sp->solutions[i]=sp->solutions[i+1];
   }
   sp->num_solutions--;
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_is_solution_in_sp (lp_prob *p, int cnt, int *indices, double *values, double obj_value)
{
   /* not implemented yet */
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int sp_initialize(tm_prob *tm)
{
   tm->sp = (sp_desc*)malloc(sizeof(sp_desc));
   sp_desc *sp = tm->sp;
   sp->max_solutions = 10;
   sp->num_solutions = 0;
   sp->solutions = (sp_solution **) malloc(sp->max_solutions*sizeof(sp_solution*));
   for (int i=0;i<sp->max_solutions;i++) {
      sp->solutions[i] = (sp_solution *) malloc(sizeof(sp_solution));
   }

}
