/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2009 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/* last modified: Feb 09, menal*/

#include <stdio.h>
#include <stdlib.h>

#include "sym_master.h"
#include "sym_macros.h"
#include "sym_prep.h"

/*===========================================================================*/
/* Accessing preprocessor through sym_environment */
/*===========================================================================*/

int sym_presolve(sym_environment *env)
{
   int termcode = 0;

   PREPdesc * P = (PREPdesc *)calloc(1, sizeof(PREPdesc));
   
   if(env->prep_mip){
      free_mip_desc(env->prep_mip);
      FREE(env->prep_mip);
   }

   P->orig_mip = env->orig_mip = create_copy_mip_desc(env->mip);
   P->mip = env->prep_mip = env->mip;
   P->params = env->par.prep_par;
   
   if(P->mip){
      termcode = prep_solve_desc(P);
   }
   
   if(termcode > -1 && P->params.reduce_mip){
      sym_restore_rootdesc(env);
   }

   /* debug */
   /*----------*/
   if(P->params.write_mps || P->params.write_lp){
      char file_name[80] = "";
      sprintf(file_name, "%s_prep", env->probname);
      
      if(P->params.write_mps){
	 sym_write_mps(env, file_name);
      }
      if(P->params.write_lp){
	 sym_write_lp(env, file_name);
      }
   }

   /* since we use the original mip desc */
   P->mip = 0;
   P->orig_mip = 0;
   
   free_prep_desc(P);

   return termcode;
}
/*===========================================================================*/
/*===========================================================================*/
int sym_restore_rootdesc(sym_environment *env)

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
/* open and initialize an environment */
/*===========================================================================*/

prep_environment * prep_open_environment()
{

   prep_environment *prep =
      (prep_environment *)calloc(1, sizeof(prep_environment));   
   prep->P = (PREPdesc *) calloc(1, sizeof(PREPdesc));
   prep->P->mip = (MIPdesc *)calloc(1, sizeof(MIPdesc));   

   /*set defaults here */

   prep_params *prep_par = &prep->params;
   prep_par->level = 5;
   prep_par->dive_level = 5;
   prep_par->impl_dive_level = 0;
   prep_par->impl_limit = 50;
   prep_par->do_probe = 1;
   prep_par->verbosity = 1;
   prep_par->reduce_mip = 1;
   prep_par->probe_verbosity = 0;
   prep_par->probe_level = 1;
   prep_par->display_stats = 0;
   prep_par->iteration_limit = 10;
   prep_par->etol = 1e-07;
   prep_par->do_single_row_rlx = 0;
   prep_par->single_row_rlx_ratio = 0.1;
   prep_par->max_sr_cnt = 5;
   prep_par->do_aggregate_row_rlx = 0;
   prep_par->max_aggr_row_cnt = 0;
   prep_par->max_aggr_row_ratio = 0.1;
   prep_par->keep_row_ordered = 1;
   prep_par->keep_track = 0;
   prep_par->time_limit = 100;
   prep_par->write_mps = 0;
   prep_par->write_lp = 0;
   
   return prep;
}

/*===========================================================================*/

int prep_solve(prep_environment *prep){

   int termcode = 0;
   PREPdesc * P = prep->P;
   
   if(P->mip){
      P->orig_mip = create_copy_mip_desc(P->mip);
      P->params = prep->params;
      termcode = prep_solve_desc(P);
   }
   
   return termcode;
}

/*===========================================================================*/
/*
   This function is the master of the preprocessing part. It calls and 
   controls other functions to perform preprocessing jobs.
*/
/*===========================================================================*/

int prep_solve_desc (PREPdesc * P) 
{

   int termcode = 0;		/* return status of this function, 0 normal, -1
				   error */
   MIPdesc *mip = P->mip;
   prep_params params = P->params;
   
   int verbosity = params.verbosity;
   int p_level = params.level;

   if (p_level <= 0) {
      if(verbosity >= 0){
	 printf ("Skipping Preprocessor\n");
      }
      return(termcode);
   }

   double start_time = wall_clock(NULL);

   /* Start with Basic Preprocessing */

   PRINT(verbosity, -2, ("Starting Preprocessing...\n"));
   
   P->stats.nz_coeff_changed = (char *)calloc(CSIZE ,mip->nz);

   /* need to fill in the row ordered vars of mip */

   /* these will be needed for both basic and advanced prep functions
      so we call them here */
   termcode = prep_fill_row_ordered(P);
   termcode = prep_initialize_mipinfo(P);//mip, params, &(P->stats));   

   /* no changes so far on column based mip*/
   /* call the main sub function of presolver */
   if(!PREP_QUIT(termcode)){
      termcode = prep_basic(P);
   }

   /* report what we have done */
   if(verbosity > -2){
      prep_report(P, termcode);
   }

   PRINT(verbosity, 0, ("Total Presolve Time: %f...\n\n", 
			wall_clock(NULL) - start_time));   

   return termcode; 
}
 
/*===========================================================================*/
/* Load a MIP model with arrays */
/*===========================================================================*/
int prep_load_problem(prep_environment *prep, int numcols, int numrows,
		      int *start, int *index, double *value,         
		      double *collb, double *colub, char *is_int,    
		      double *obj, double obj_offset, char *rowsen,       
		      double *rowrhs, double *rowrng, char make_copy)
{
   int termcode = 0;   
   double inf = INF;
   int i = 0;
   MIPdesc *mip;
   
   if ((!numcols && !numrows) || numcols < 0 || numrows <0){
      printf("prep_load_problem():The given problem description is"
	     "empty or incorrect ");
      return(PREP_FUNC_ERROR);
   }

   mip = prep->P->mip;
   
   mip->m  = numrows;
   mip->n  = numcols;

   if (make_copy){      
      
      if(numcols){
	 mip->obj    = (double *) calloc(numcols, DSIZE);
	 mip->ub     = (double *) calloc(numcols, DSIZE);
	 mip->lb     = (double *) calloc(numcols, DSIZE);
	 mip->is_int = (char *)   calloc(CSIZE, numcols);

	 if (obj){
	    memcpy(mip->obj,  obj,  DSIZE * numcols);
	 }

	 if (colub){
	    memcpy(mip->ub, colub, DSIZE * numcols); 
	 }else{
	    for(i = 0; i<mip->n; i++){
	       mip->ub[i] = inf;
	    }
	 }
	 
	 if(collb){
	    memcpy(mip->lb, collb, DSIZE * numcols);
	 }
	 
	 if (is_int){
	    memcpy(mip->is_int, is_int, CSIZE * numcols);
	 }
      }

      if(numrows){

	 mip->rhs    = (double *) calloc(numrows, DSIZE);
	 mip->sense  = (char *)   malloc(CSIZE * numrows);
	 mip->rngval = (double *) calloc(numrows, DSIZE);

	 if (rowsen){
	    memcpy(mip->sense, rowsen, CSIZE * numrows); 
	 }else{
	    memset(mip->sense, 'N', CSIZE *numrows);
	 }
	 
	 if(rowrhs){
	    memcpy(mip->rhs, rowrhs, DSIZE * numrows);
	 }
	 
	 if (rowrng){
	    memcpy(mip->rngval, rowrng, DSIZE * numrows);
	 }
      }
      
      //user defined matind, matval, matbeg--fill as column ordered
      
      if(start){      

	 mip->nz = start[numcols];
	 mip->matbeg = (int *) calloc(ISIZE, (numcols + 1));
	 mip->matval = (double *) calloc(DSIZE,start[numcols]);
	 mip->matind = (int *)    calloc(ISIZE,start[numcols]);
	 
	 memcpy(mip->matbeg, start, ISIZE *(numcols + 1));
	 memcpy(mip->matval, value, DSIZE *start[numcols]);  
	 memcpy(mip->matind, index, ISIZE *start[numcols]);  
      }
      
   }else{
      
      if (obj){
	 mip->obj = obj;
      }else{
	 mip->obj    = (double *) calloc(numcols, DSIZE);	 
      }

      if (rowsen){
	 mip->sense = rowsen;
      }else{
	 mip->sense  = (char *) malloc(CSIZE * numrows);
	 memset(mip->sense, 'N', CSIZE *numrows);
      }

      if(rowrhs){
	 mip->rhs = rowrhs;
      }else{
	 mip->rhs = (double *) calloc(numrows, DSIZE);	 
      }

      if (rowrng){
	 mip->rngval = rowrng;
      }else{
	 mip->rngval = (double *) calloc(numrows, DSIZE);
      }

      if (colub){
	 mip->ub = colub;
      }else{
	 mip->ub = (double *) calloc(numcols, DSIZE);
	 for(i = 0; i<mip->n; i++){
	    mip->ub[i] = inf;
	 }
      }

      if (collb){
	 mip->lb = collb;
      }else{
	 mip->lb = (double *) calloc(numcols, DSIZE);	 
      }

      if (is_int){
	 mip->is_int = is_int;
      }else{
	 mip->is_int = (char *)   calloc(CSIZE, numcols);
      }

      if(start){
	 mip->nz = start[numcols];
	 mip->matbeg = start;
	 mip->matval = value;
	 mip->matind = index;
      }
   }

   mip->obj_offset = -obj_offset;

   return termcode;
}

/*************************************************************************
 ***                     preprocessing - parameters                    ***
 *************************************************************************/ 

int prep_set_param(prep_environment *prep, char *key, int value)
{

   prep_params *prep_par = &prep->params;
   
   //if (strcmp(key, "prep_do_preprocessing") == 0){
   //  prep_par->do_prep = value;
   //  return(0);
   //}
   if (strcmp(key, "prep_level") == 0){
      prep_par->level = value;
      return(0);
   }
   else if (strcmp(key, "prep_dive_level") == 0){
      prep_par->dive_level = value;
      return(0);
   }
   else if (strcmp(key, "prep_impl_dive_level") == 0){
      prep_par->impl_dive_level = value;
      return(0);
   }
   else if (strcmp(key, "prep_impl_limit") == 0){
      prep_par->impl_limit = value;
      return(0);
   }
   else if (strcmp(key, "prep_iter_limit") == 0){
      prep_par->iteration_limit = value;
      return(0);
   }
   else if (strcmp(key, "prep_do_probing") == 0){
      prep_par->do_probe = value;
      return(0);
   }
   else if (strcmp(key, "prep_do_sr") == 0){
      prep_par->do_single_row_rlx = value;
      return(0);
   }
   else if (strcmp(key, "prep_verbosity") == 0){
      prep_par->verbosity = value;
      return(0);
   }
   else if (strcmp(key, "prep_reduce_mip") == 0){
      prep_par->reduce_mip = value;
      return(0);
   }
   else if (strcmp(key, "prep_probing_verbosity") == 0){
      prep_par->probe_verbosity = value;
      return(0);
   }
   else if (strcmp(key, "prep_probing_level") == 0){
      prep_par->probe_level = value;
      return(0);
   }
   else if (strcmp(key, "prep_display_stats") == 0){
      prep_par->display_stats = value;
      return(0);
   }
   else if (strcmp(key, "max_sr_cnt") == 0){
      prep_par->max_sr_cnt = value;
      return(0);
   }
   else if (strcmp(key, "max_aggr_row_cnt") == 0){
      prep_par->max_aggr_row_cnt = value;
      return(0);
   }
   else if (strcmp(key, "keep_row_ordered") == 0){
      prep_par->keep_row_ordered = value;
      return(0);
   }
   else if (strcmp(key, "write_mps") == 0){
      prep_par->write_mps = value;
      return(0);
   }
   else if (strcmp(key, "write_lp") == 0){
      prep_par->write_lp = value;
      return(0);
   }
   else if (strcmp(key, "prep_time_limit") == 0){
      prep_par->time_limit = value;
      return(0);
   }

   return(PREP_FUNC_ERROR);
}

/*===========================================================================*/
/*===========================================================================*/
/* -not used here */
#if 0
int prep_read_mps(prep_environment *prep, char *infile)
{
   int j, errors;
   MIPdesc *mip = prep->P->mip;
   CoinMpsIO mps;
   
   mps.messageHandler()->setLogLevel(0);
   
#if 0

   int j, last_dot = 0, last_dir = 0;
   char fname[80] = "", ext[10] = "";

   size_t size = 1000;
   char* buf = 0;

   while (true) {
      buf = (char*)malloc(CSIZE*size);
      if (getcwd(buf, size))
	 break;
      FREE(buf);
      buf = 0;
      size = 2*size;
   }
   char slash = buf[0] == '/' ? '/' : '\\';
   FREE(buf);
   
   for (j = 0;; j++){
      if (infile[j] == '\0')
	 break;
      if (infile[j] == '.') {
	    last_dot = j;
	  }
	  if(infile[j] == slash){
		last_dir = j;
	  }
   }
   
   if(last_dir < last_dot){
	   memcpy(fname, infile, CSIZE*last_dot);
	   memcpy(ext, infile + last_dot + 1, CSIZE*(j - last_dot - 1)); 
   }
   else{
	   memcpy(fname, infile, CSIZE*j);
   }
#endif

   mps.setInfinity(mps.getInfinity());

   if (mps.readMps(infile,"")){
      return(PREP_FUNC_ERROR);
   }
   
   //strncpy(probname, const_cast<char *>(mps.getProblemName()), 80);
   
   mip->m  = mps.getNumRows();
   mip->n  = mps.getNumCols();
   mip->nz = mps.getNumElements();
   
   mip->obj    = (double *) malloc(DSIZE * mip->n);
   mip->rhs    = (double *) malloc(DSIZE * mip->m);
   mip->sense  = (char *)   malloc(CSIZE * mip->m);
   mip->rngval = (double *) malloc(DSIZE * mip->m);
   mip->ub     = (double *) malloc(DSIZE * mip->n);
   mip->lb     = (double *) malloc(DSIZE * mip->n);
   mip->is_int = (char *)   calloc(CSIZE, mip->n);
   
   memcpy(mip->obj, const_cast <double *> (mps.getObjCoefficients()),
	  DSIZE * mip->n); 
   memcpy(mip->rhs, const_cast <double *> (mps.getRightHandSide()),
	  DSIZE * mip->m); 
   memcpy(mip->sense, const_cast <char *> (mps.getRowSense()),
	  CSIZE * mip->m); 
   memcpy(mip->rngval, const_cast <double *> (mps.getRowRange()),
	  DSIZE * mip->m); 
   memcpy(mip->ub, const_cast <double *> (mps.getColUpper()),
	  DSIZE * mip->n); 
   memcpy(mip->lb, const_cast <double *> (mps.getColLower()),
	  DSIZE * mip->n); 
   
   //user defined matind, matval, matbeg--fill as column ordered
   
   const CoinPackedMatrix * matrixByCol= mps.getMatrixByCol();
   
   mip->matbeg = (int *) malloc(ISIZE * (mip->n + 1));
   memcpy(mip->matbeg, const_cast<int *>(matrixByCol->getVectorStarts()),
	  ISIZE * (mip->n + 1));
   
   mip->matval = (double *) malloc(DSIZE*mip->matbeg[mip->n]);
   mip->matind = (int *)    malloc(ISIZE*mip->matbeg[mip->n]);
   
   memcpy(mip->matval, const_cast<double *> (matrixByCol->getElements()),
	  DSIZE * mip->matbeg[mip->n]);  
   memcpy(mip->matind, const_cast<int *> (matrixByCol->getIndices()), 
	  ISIZE * mip->matbeg[mip->n]);  

   mip->colname = (char **) malloc(sizeof(char *) * mip->n);  

   for (j = 0; j < mip->n; j++){
      mip->is_int[j] = mps.isInteger(j);
      mip->colname[j] = (char *) malloc(CSIZE * 9);
      strncpy(mip->colname[j], const_cast<char*>(mps.columnName(j)), 9);
      mip->colname[j][8] = 0;
   }

#if 0
   if (mip->obj_sense == SYM_MAXIMIZE){
      for (j = 0; j < mip->n; j++){
	 mip->obj[j] *= -1.0;
      }
   }
#endif
   
   mip->obj_offset = -mps.objectiveOffset();

   return(0);
}

/*===========================================================================*/
/*===========================================================================*/

int prep_read_lp(prep_environment *prep, char *infile)
{

   int j;
   CoinLpIO lp;
   MIPdesc *mip = prep->P->mip;
   
   lp.readLp(infile);

   /*
   if (lp.readLp(infile)){
      return(PREP_FUNC_ERROR);
   }
   */
   //strncpy(probname, const_cast<char *>(lp.getProblemName()), 80);
   
   mip->m  = lp.getNumRows();
   mip->n  = lp.getNumCols();
   mip->nz = lp.getNumElements();
   
   mip->obj    = (double *) malloc(DSIZE * mip->n);
   mip->rhs    = (double *) malloc(DSIZE * mip->m);
   mip->sense  = (char *)   malloc(CSIZE * mip->m);
   mip->rngval = (double *) malloc(DSIZE * mip->m);
   mip->ub     = (double *) malloc(DSIZE * mip->n);
   mip->lb     = (double *) malloc(DSIZE * mip->n);
   mip->is_int = (char *)   calloc(CSIZE, mip->n);
   
   memcpy(mip->obj, const_cast <double *> (lp.getObjCoefficients()),
	  DSIZE * mip->n); 
   memcpy(mip->rhs, const_cast <double *> (lp.getRightHandSide()),
	  DSIZE * mip->m); 
   memcpy(mip->sense, const_cast <char *> (lp.getRowSense()),
	  CSIZE * mip->m); 
   memcpy(mip->rngval, const_cast <double *> (lp.getRowRange()),
	  DSIZE * mip->m); 
   memcpy(mip->ub, const_cast <double *> (lp.getColUpper()),
	  DSIZE * mip->n); 
   memcpy(mip->lb, const_cast <double *> (lp.getColLower()),
	  DSIZE * mip->n); 
   
   //user defined matind, matval, matbeg--fill as column ordered
   
   const CoinPackedMatrix * matrixByCol= lp.getMatrixByCol();
   
   mip->matbeg = (int *) malloc(ISIZE * (mip->n + 1));
   memcpy(mip->matbeg, const_cast<int *>(matrixByCol->getVectorStarts()),
	  ISIZE * (mip->n + 1));
   
   mip->matval = (double *) malloc(DSIZE*mip->matbeg[mip->n]);
   mip->matind = (int *)    malloc(ISIZE*mip->matbeg[mip->n]);
   
   memcpy(mip->matval, const_cast<double *> (matrixByCol->getElements()),
	  DSIZE * mip->matbeg[mip->n]);  
   memcpy(mip->matind, const_cast<int *> (matrixByCol->getIndices()), 
	  ISIZE * mip->matbeg[mip->n]);  

   mip->colname = (char **) malloc(sizeof(char *) * mip->n); 

   for (j = 0; j < mip->n; j++){
      mip->is_int[j] = lp.isInteger(j);
      mip->colname[j] = (char *) malloc(CSIZE * 9);
      strncpy(mip->colname[j], const_cast<char*>(lp.columnName(j)), 9);
      mip->colname[j][8] = 0;
   }

#if 0
   if (mip->obj_sense == SYM_MAXIMIZE){
      for (j = 0; j < mip->n; j++){
	 mip->obj[j] *= -1.0;
      }
   }
#endif
   
   mip->obj_offset = -lp.objectiveOffset();

   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
void prep_write_mps(prep_environment *prep, char *outfile)
{
   int i;
   char filename[80] = "";
   CoinMpsIO mps;
   MIPdesc *mip = prep->P->mip;
   
   CoinPackedMatrix mip_matrix(true, mip->m, mip->n, mip->nz, mip->matval, 
			       mip->matind, mip->matbeg, 0);

   mps.setMpsData(mip_matrix, mps.getInfinity(), mip->lb, mip->ub, mip->obj, 
		  mip->is_int, mip->sense, mip->rhs, mip->rngval, 
		  mip->colname, NULL);
   mps.setObjectiveOffset(mip->obj_offset); 

   sprintf(filename, "%s%s%s", outfile, ".","MPS");
   mps.writeMps(filename);
}

/*===========================================================================*/
/*===========================================================================*/
void prep_write_lp(prep_environment *prep, char *outfile)
{
   int i;
   double * rlb, * rub, infinity; 
   char filename[80] = "";
   CoinLpIO lp;
   MIPdesc *mip = prep->P->mip;
   
   CoinPackedMatrix mip_matrix(true, mip->m, mip->n, mip->nz, mip->matval, 
			       mip->matind, mip->matbeg, 0);

   rlb = (double *) malloc(DSIZE*mip->m);
   rub = (double *) malloc(DSIZE*mip->m);
   infinity = lp.getInfinity();

   /* convert sense to bound */
   for(i = 0; i < mip->m; i++){
      switch (mip->sense[i]){
       case 'E':
	  rlb[i] = rub[i] = mip->rhs[i];
	  break;
       case 'L':
	  rlb[i] = -infinity;
	  rub[i] = mip->rhs[i];
	  break;
       case 'G':
	  rlb[i] = mip->rhs[i];
	  rub[i] = infinity;
	  break;
       case 'R':
	  rlb[i] = mip->rhs[i] - mip->rngval[i];
	  rub[i] = mip->rhs[i];
	  break;
       case 'N':
	  rlb[i] = -infinity;
	  rub[i] = infinity;
	  break;
      }
   }
   
   lp.setLpDataWithoutRowAndColNames(mip_matrix, mip->lb, mip->ub, mip->obj, 
				     mip->is_int, rlb, rub);
   lp.setObjectiveOffset(mip->obj_offset); 
   lp.setLpDataRowAndColNames(NULL, mip->colname);
   sprintf(filename, "%s%s%s", outfile, ".","LPT");
   lp.writeLp(filename);

   FREE(rlb);
   FREE(rub);
}
#endif
/*===========================================================================*/
/*===========================================================================*/
void prep_close_environment(prep_environment *prep)
{
   free_prep_desc(prep->P);
   FREE(prep);
}

/*===========================================================================*/
