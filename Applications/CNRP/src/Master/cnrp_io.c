/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* Capacitated Network Routing Problems.                                     */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/* system include files */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

/* SYMPHONY include files */
#include "BB_macros.h"
#include "BB_types.h"
#include "qsortucb.h"
#include "master_u.h"
#include "lp_params.h"

/* CNRP include files */
#include "cnrp_io.h"
#include "compute_cost.h"
#include "cnrp_types.h"
#include "cnrp_const.h"
#include "cnrp_macros.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user I/O functions for the master process.
\*===========================================================================*/

/*===========================================================================*\
 * This first function reads in the data instance.
\*===========================================================================*/

void vrp_io(vrp_problem *vrp, char *infile)
{
  static char keywords[KEY_NUM][22] = {
    "NAME", 
    "NAME:",                 /*This section lists the names of the */
    "TYPE",                  /*possible fields in the data file    */
    "TYPE:",
    "COMMENT",
    "COMMENT:",
    "DIMENSION",
    "DIMENSION:",
    "CAPACITY",
    "CAPACITY:",
    "EDGE_WEIGHT_TYPE",
    "EDGE_WEIGHT_TYPE:",
    "EDGE_WEIGHT_FORMAT", 
    "EDGE_WEIGHT_FORMAT:", 
    "DISPLAY_DATA_TYPE",
    "DISPLAY_DATA_TYPE:",
    "EDGE_WEIGHT_SECTION", 
    "EDGE_WEIGHT_SECTION:", 
    "DISPLAY_DATA_SECTION", 
    "DISPLAY_DATA_SECTION:",
    "NODE_COORD_SECTION",
    "NODE_COORD_SECTION:",
    "NODE_COORD_TYPE",
    "NODE_COORD_TYPE:",
    "DEPOT_SECTION",
    "DEPOT_SECTION:",
    "CAPACITY_VOL",
    "CAPACITY_VOL:",
    "DEMAND_SECTION",
    "DEMAND_SECTION:",
    "TIME_WINDOW_SECTION",
    "TIME_WINDOW_SECTION:",
    "STANDTIME_SECTION",
    "STANDTIME_SECTION:",
    "PICKUP_SECTION",
    "PICKUP_SECTION:",
    "EOF",
    "EOF.",
    "",
    "",
    "NO_MORE_TYPE"
  };

#define NCTYPE_NUM 3

  static char nctypes[NCTYPE_NUM][14] = {
    "TWOD_COORDS",
    "THREED_COORDS",     /*This section lists the possible node*/
    "NO_COORDS"          /*coordinate data types               */
  };

#define WTYPE_NUM 8

  static char wtypes[WTYPE_NUM][9] = {
    "EXPLICIT",
    "EUC_2D",            /*This is a list of the possible data types for */
    "EUC_3D",            /*edge weights                                  */
    "MAX_2D",
    "MAX_3D",
    "MAN_2D",
    "MAN_3D",
    "GEO"
  };

#define WFORMAT_NUM 9

  static char wformats[WFORMAT_NUM][20] = {
    "UPPER_ROW",
    "LOWER_ROW",          /*This is a list of the possible formats that*/
    "UPPER_DIAG_ROW",     /*the edge weight matrix could be given in   */
    "LOWER_DIAG_ROW",     /*if it is given explicitly                  */
    "UPPER_COL",
    "LOWER_COL",
    "UPPER_DIAG_COL",
    "LOWER_DIAG_COL",
    "FULL_MATRIX"
  };

#define DTYPE_NUM 3

  static char dtypes[DTYPE_NUM][14] = {
    "COORD_DISPLAY",
    "TWOD_DISPLAY",     /*This is a list of the various display data*/
    "NO_DISPLAY"        /*types                                     */
  };

  char line[LENGTH], line1[LENGTH], key[30], tmp[80];
  int wformat=-1, dtype=-1, nctype=-1;
  double fdummy;
  int i, j = 0;
  int l, m, *coef2;
  FILE *f;
  int node;
  double deg, min, coord_x, coord_y, coord_z;
  double x, y;
  int capacity_vol = FALSE;
  int k;
  register int vertnum = 0;
  distances *dist = &vrp->dist;

  if (!strcmp(infile, "")){
     printf("\nVrp I/O: No problem data file specified\n\n");
     exit(1);
  }
  
  if ((f = fopen(infile, "r")) == NULL){
     fprintf(stderr, "Vrp I/O: file '%s' can't be opened\n", infile);
     exit(1);
  }
  
  /*This loop reads in the next line of the data file and compares it
    to the list of possible keywords to determine what data will follow.
    It then reads the data into the appropriate field and iterates */
  
  while(NULL != fgets( line1, LENGTH, f)){
     strcpy(key,"");
     sscanf(line1,"%s",key); /*read in next keyword*/
     
     for (k = 0; k < KEY_NUM; k++) /*which data field comes next?*/
	if (strcmp(keywords[k], key) == 0) break;
     
     if (k == KEY_NUM){
	continue;
	fprintf(stderr, "Unknown keyword! bye.\n");
	exit(1); /*error check for acceptable data field*/
     }
     
     k >>= 1; /* This is a bit shift operation that divides k by 2    */
              /* since in the list of keywords, there are two possible*/
              /* formats for the keyword                              */
    
     if (strchr(line1,':')){
	strcpy(line, strchr(line1, ':')+1);
     }
     
     switch (k){
	
      case 0: /* NAME */
	if (!sscanf(line, "%s", vrp->name))
	   fprintf(stderr, "\nVrp I/O: error reading NAME\n\n");
	printf("PROBLEM NAME: \t\t%s\n", vrp->name);
	break;

      case 1 : /*TYPE*/
	sscanf(line, "%s", tmp);
	if (!vrp->par.prob_type){
	   if (strcmp("CVRP", tmp) == 0){
	      vrp->par.prob_type = vrp->cg_par.prob_type =
		 vrp->lp_par.prob_type = VRP;
	   }else if (strcmp("TSP", tmp) == 0){
	      vrp->par.prob_type = vrp->cg_par.prob_type = 
		 vrp->lp_par.prob_type = TSP;
	   }else if (strcmp("BPP", tmp) == 0){
	      vrp->par.prob_type = vrp->cg_par.prob_type = 
		 vrp->lp_par.prob_type = BPP;
	   }else if (strcmp("CSTP", tmp) == 0){
	      vrp->par.prob_type = vrp->cg_par.prob_type = 
		 vrp->lp_par.prob_type = CSTP;
	   }else if (strcmp("CTP", tmp) == 0){
	      vrp->par.prob_type = vrp->cg_par.prob_type = 
		 vrp->lp_par.prob_type = CTP;
	   }else{
	      fprintf(stderr, "This is not a recognized problem type!\n");
	      exit(1);
	   }
	}
	switch(vrp->par.prob_type){
	 case VRP:
	   printf("TYPE: \t\t\tCVRP\n");
	   break;
	   
	 case TSP:
	   printf("TYPE: \t\t\tTSP\n");
	   break;

	 case BPP:
	   printf("TYPE: \t\t\tBPP\n");
	   break;

	 case CSTP:
	   printf("TYPE: \t\t\tCSTP\n");
	   break;

	 case CTP:
	   printf("TYPE: \t\t\tCTP\n");
	   break;
	}
	break;

      case 2 : /*COMMENT*/
#if 0
	if (!strncpy(tmp, line, 80))
	   fprintf(stderr, "\nVrp I/O: error reading COMMENT\n\n");
	printf("DESCRIPTION: \t\t%s\n", tmp);
#endif
	break;
      case 3 : /* DIMENSION */
	if (!sscanf(line, "%i", &k)){
	   fprintf(stderr, "Vrp I/O: error reading DIMENSION\n\n");
	   exit(1);
	}
	vertnum = vrp->vertnum = (int) k;
	vrp->edgenum = (int) vertnum * (vertnum - 1)/2;
	printf("DIMENSION: \t\t%i\n", k);
	break;
      case 4 : /*CAPACITY*/
	if (!sscanf(line, "%i", &k)){
	   fprintf(stderr, "Vrp I/O: error reading CAPACITY\n\n");
	   exit(1);
	}
	vrp->capacity = (int) k;
	break;
      case 5 : /* EDGE_WEIGHT_TYPE */
	sscanf(line, "%s", tmp);
	for (dist->wtype = 0; dist->wtype < WTYPE_NUM; (dist->wtype)++)
	   if (strcmp(wtypes[dist->wtype], tmp) == 0) break;
	if (dist->wtype == WTYPE_NUM) {
	   fprintf(stderr, "Unknown weight type : %s !!!\n", tmp);
	   exit(1);
	}
	break;
      case 6 : /* EDGE_WEIGHT_FORMAT */
	sscanf(line, "%s", tmp);
	for (wformat = 0; wformat < WFORMAT_NUM; wformat++)
	   if (strcmp(wformats[wformat], tmp) == 0) break;
	if (wformat == WFORMAT_NUM) {
	   fprintf(stderr, "Unknown weight type : %s !!!\n", tmp);
	   exit(1);
	}
	break;
      case 7 : /* DISPLAY_DATA_TYPE */
	sscanf(line, "%s", tmp);
	for (dtype = 0; dtype < DTYPE_NUM; dtype++)
	   if (strcmp(dtypes[dtype], tmp) == 0) break;
	if (dtype == DTYPE_NUM) {
	   fprintf(stderr, "Unknown display type : %s !!!\n", tmp);
	   exit(1);
	}
	break;
      case 8: /* EDGE_WEIGHT_SECTION */
	/*------------------------break if not EXPLICIT -*/
	if (dist->wtype != _EXPLICIT) break; 
	dist->cost = (int *) malloc (vrp->edgenum*sizeof(int));
	switch (wformat){
	 case 1 : /* LOWER_ROW */
	 case 4 : /* UPPER_COL */
	 case 3 : /* LOWER_DIAG_ROW */
	 case 6 : /* UPPER_DIAG_COL */
	   for (i=0, coef2=dist->cost; i<vertnum; i++){
	      for (j=0; j<i; j++, coef2++){
		 if (!fscanf(f,"%lf", &fdummy)){
		    fprintf(stderr, "Not enough data -- DIMENSION or "
			    "EDGE_WEIGHT_TYPE declared wrong\n");
		    exit(1);
		 }
		 else *coef2 = (int) fdummy;
	      }
	      if ((wformat==3 || wformat==6) && 
		  !fscanf(f,"%lf", &fdummy)){
		 fprintf(stderr, "Not enough data -- DIMENSION or "
			 "EDGE_WEIGHT_TYPE declared wrong\n");
		 exit(1);
	      }
	   }
	   if (fscanf(f,"%lf", &fdummy)){
	      fprintf(stderr, "Too much data -- DIMENSION or "
		      "EDGE_WEIGHT_TYPE declared wrong\n");
	      exit(1);
	   }
	   break;
	 case 0 : /* UPPER_ROW */
	 case 5 : /* LOWER_COL */
	 case 2 : /* UPPER_DIAG_ROW */
	 case 7 : /* LOWER_DIAG_COL */
	   for (i=0, coef2=dist->cost; i<vertnum; i++){
	      if (wformat==2 || wformat==7) 
		 if (!fscanf(f,"%lf", &fdummy)){
		    fprintf(stderr, "Not enough data -- DIMENSION or "
			    "EDGE_WEIGHT_TYPE declared wrong");
		    exit(1);
		 }
	      for (j=i+1; j<vertnum; j++){
		 if (!fscanf(f,"%lf", &fdummy)){
		    fprintf(stderr, "Not enough data -- DIMENSION or "
			    "EDGE_WEIGHT_TYPE declared wrong");
		    exit(1);
		 }
		 else coef2[j*(j-1)/2+i] = (int) fdummy;
	      }
	   }
	   if (fscanf(f,"%lf", &fdummy)){
	      fprintf(stderr, "Too much data -- DIMENSION or "
		      "EDGE_WEIGHT_TYPE declared wrong\n");
	      exit(1);
	   }
	   break;
	 case 8 : /* FULL_MATRIX */
	   for (i=0, coef2=dist->cost; i<vertnum; i++){
	      for (j=0; j<=i; j++)
		 if(!fscanf(f,"%lf", &fdummy)){
		    fprintf(stderr, "Not enough data -- DIMENSION or "
			    "EDGE_WEIGHT_TYPE declared wrong");
		    exit(1);
		 }
	      for (j=i+1; j<vertnum; j++){
		 if(!fscanf(f,"%lf", &fdummy)){
		    fprintf(stderr, "Not enough data -- DIMENSION or "
			    "EDGE_WEIGHT_TYPE declared wrong");
		    exit(1);
		 }
		 coef2[j*(j-1)/2+i] = (int) fdummy;
	      }
	   }
	   if (fscanf(f,"%lf", &fdummy)){
	      fprintf(stderr, "Too much data -- DIMENSION or "
		      "EDGE_WEIGHT_TYPE declared wrong\n");
	      exit(1);
	   }
	   break;
	}
	break;
      case 9 : /* DISPLAY_DATA_SECTION */
	/*--------------------- break if NO_DISPLAY -*/
	if (dtype != 1){
	   fprintf(stderr, "DISPLAY_DATA_SECTION exists"
		   "but not TWOD_DISPLAY!\n");
	   exit(1);
	}
	/* posx, posy -*/
	vrp->posx = (int *) malloc (vertnum*sizeof(int));
	vrp->posy = (int *) malloc (vertnum*sizeof(int));
	for (i=0; i<vertnum; i++){
	   if ((k = fscanf(f,"%i%lf%lf", &node, &x, &y)) != 3){
	      fprintf(stderr, "\nVrp I/O: error reading DISPLAY_DATA\n");
	      break;
	   }
	   vrp->posx[node-1] = (int)(x + 0.5);
	   vrp->posy[node-1] = (int)(y + 0.5);
	}
	if (fscanf(f,"%lf", &fdummy)){
	   fprintf(stderr, "\nVrp I/O: too much display data\n");
	   break;
	}
	break;
      case 10 : /* NODE_COORD_SECTION */
	if (nctype == -1) nctype = 0;  /*if not given: TWOD_COORDS*/
	if (dtype == -1 && ((dist->wtype == _EUC_2D) || /*display type*/
			    (dist->wtype == _MAX_2D) ||  /*not defd yet*/
			    (dist->wtype == _MAN_2D)   ))/*&& can disp.*/
	   dtype = 0;                               /* COORD_DISPLAY */
	if (dtype == 0){
	   vrp->posx = (int *) malloc (vertnum*sizeof(int));
	   vrp->posy = (int *) malloc (vertnum*sizeof(int));
	}
	dist->coordx = (double *) malloc (vertnum*sizeof(double));
	dist->coordy = (double *) malloc (vertnum*sizeof(double));
	if (nctype == 1)
	   dist->coordz = (double *) malloc (vertnum*sizeof(double));
	for (i=0; i<vertnum; i++){
	   if (nctype == 0)          /* TWOD_COORDS */
	      if (fscanf(f,"%i%lf%lf", &node, &coord_x, &coord_y) != 3){
		 fprintf(stderr, "\nVrp I/O: error reading NODE_COORD\n\n");
		 exit(1);
	      }
	   if (nctype == 1)          /* THREED_COORDS */
	      if (fscanf(f,"%i%lf%lf%lf", &node, &coord_x, &coord_y,
			 &coord_z) != 4){
		 fprintf(stderr, "\nVrp I/O: error reading NODE_COORD\n\n");
		 exit(1);
	      }
	   dist->coordx[node-1] = coord_x;
	   dist->coordy[node-1] = coord_y;
	   /*since position is an integer and coord is a double, I must
	     round off here if dtype is EXPLICIT*/
	   if (dtype == 0){
	      vrp->posx[node-1] = (int)coord_x;
	      vrp->posy[node-1] = (int)coord_y;
	   }
	   if (nctype == 1) dist->coordz[node-1] = coord_z;
	   if (dist->wtype == 7){ /* GEO */
	      /*--- latitude & longitude for i ------------*/
	      deg = floor(dist->coordx[i]);
	      min = dist->coordx[i] - deg;
	      dist->coordx[i] = MY_PI * (deg + 5.0*min/3.0 ) / 180.0;
	      deg = floor(dist->coordy[i]);
	      min = dist->coordy[i] - deg;
	      dist->coordy[i] = MY_PI * (deg + 5.0*min/3.0 ) / 180.0;
	      /*--- latitude & longitude for j ------------*/
	      deg = floor(dist->coordx[j]);
	      min = dist->coordx[j] - deg;
	      dist->coordx[j] = MY_PI * (deg + 5.0*min/3.0 ) / 180.0;
	      deg = floor(dist->coordy[j]);
	      min = dist->coordy[j] - deg;
	      dist->coordy[j] = MY_PI * (deg + 5.0*min/3.0 ) / 180.0;
	   }
	}
	if (fscanf(f,"%i%lf%lf%lf", &node, &coord_x, &coord_y, &coord_z)){
	   fprintf(stderr, "\nVrp I/O: too much data in NODE_COORD\n\n");
	   exit(1);
	}
	break;
      case 11: /* NODE_COORD_TYPE */
	sscanf(line, "%s", tmp);
	for (nctype = 0; nctype < NCTYPE_NUM; nctype++)
	   if (strcmp(nctypes[nctype], tmp) == 0) break;
	if (nctype == NCTYPE_NUM) {
	   fprintf(stderr, "Unknown node_coord_type : %s !!!\n", tmp);
	   exit(1);
	}
	break;
      case 12: /*DEPOT_SECTION*/
	fscanf(f, "%i", &k);
	if (k != 1){
	   fprintf(stderr, "Error in data: depot must be node 1");
	   exit(1);
	}
	vrp->depot = k - 1;
	while (-1 != k) fscanf(f, "%i", &k);
	break;
      case 13: /*CAPACITY_VOL*/
	sscanf(line, "%i", &k);
	capacity_vol = TRUE;
	break;
      case 14: /*DEMAND_SECTION*/
	vrp->demand = (int *) malloc(vertnum*sizeof(int));
	for (i = 0; i < vertnum; i++){
	   if (capacity_vol){
	      if (fscanf(f, "%i%i%i", &k, &l, &m) != 3){
		 fprintf(stderr,"\nVrp I/O: error reading DEMAND_SECTION\n\n");
		 exit(1);
	      }
	   }
	   else if (fscanf(f, "%i%i", &k, &l) != 2){
	      fprintf(stderr, "\nVrp I/O: error reading DEMAND_SECTION\n\n");
	      exit(1);
	   }
	   vrp->demand[k-1] = l;
	   vrp->demand[0] += l;
	}
	if (fscanf(f, "%i%i", &k, &l)){
	   fprintf(stderr, "\nVrp I/O: too much data in DEMAND_SECTION\n\n");
	   exit(1);
	}
	break;
      case 15: /*TIME_WINDOW_SECTION*/  /*These sections are not used*/
	while (fscanf(f, "%d %*d:%*d %*d:%*d", &k));
	break;
      case 16: /*STANDTIME_SECTION*/
	while (fscanf(f, "%d%*d", &k));
	break;
      case 17: /*PICKUP_SECTION*/       
	while (fscanf(f, "%d%*d%*d", &k));
	break;
      case 18: /*  EOF  */
      default: break;
     }
  }
  
  if (f != stdin)
     fclose(f);
  
  vrp->cur_tour = (best_tours *) calloc(1, sizeof(best_tours));
  vrp->cur_tour->tour = (_node *) calloc(vertnum, sizeof(_node));
  
  /*calculate all the distances explcitly and then use distance type EXPLICIT*/

  if (vrp->par.prob_type != VRP && vrp->par.prob_type != BPP &&
      vrp->par.prob_type != TSP && vrp->par.prob_type != CSTP &&
      vrp->par.prob_type != CTP){
     printf("Unknown problem type! Exiting...\n");
     exit(1);
  }
  
  if (vrp->par.prob_type == BPP){
     dist->cost = (int *) calloc (vrp->edgenum, sizeof(int));
     for (i = 1, k = 0; i < vertnum; i++){
	for (j = 0; j < i; j++){
	   dist->cost[k++] = vrp->demand[i]+vrp->demand[j];
	}
     }
  }else if (dist->wtype != _EXPLICIT){
     dist->cost = (int *) calloc (vrp->edgenum, sizeof(int));
     for (i = 1, k = 0; i < vertnum; i++){
	for (j = 0; j < i; j++){
	   dist->cost[k++] = ICOST(dist, i, j);
	}
     }
  }
  dist->wtype = _EXPLICIT;
  
  if (vrp->par.k_closest < 0){
     vrp->par.k_closest = (int) ceil(0.1 * vrp->vertnum);
     if (vrp->par.k_closest < vrp->par.min_closest ) 
	vrp->par.k_closest = vrp->par.min_closest;
     if (vrp->par.k_closest > vrp->par.max_closest) 
	vrp->par.k_closest = vrp->par.max_closest;
     if (vrp->par.k_closest > vertnum-1) 
	vrp->par.k_closest = vertnum-1;
  }
  if (vrp->par.prob_type == TSP || vrp->par.prob_type == CTP){
     vrp->numroutes = 1;
     if (!vrp->demand)
	vrp->demand = (int *) malloc (vertnum * ISIZE);
     for (i = vertnum - 1; i > 0; i--)
	vrp->demand[i] = 1;
     if (!vrp->cg_par.which_tsp_cuts)
	vrp->cg_par.which_tsp_cuts = ALL_TSP_CUTS;
  }
  
  if (vrp->par.prob_type == TSP)
     vrp->capacity = vrp->demand[0] = vertnum;
  else if (vrp->par.prob_type == CTP)
     vrp->capacity = vrp->demand[0] = vertnum-1;
}

/*===========================================================================*/

/*===========================================================================*\
 * This second function reads in the parameters from the parameter file.
\*===========================================================================*/

void vrp_readparams(vrp_problem *vrp, char *filename, int argc, char **argv)
{
   int i, j;
   char line[LENGTH], key[50], value[50], c, tmp;
   FILE *f = NULL;
   str_int colgen_str[COLGEN_STR_SIZE] = COLGEN_STR_ARRAY;
   
   vrp_params *par = &vrp->par;
   lp_user_params *lp_par = &vrp->lp_par;
   cg_user_params *cg_par = &vrp->cg_par;

   vrp->numroutes = 0;
#if defined(CHECK_CUT_VALIDITY) || defined(TRACE_PATH)
   vrp->feas_sol_size = 0;
   vrp->feas_sol = NULL;
#endif
   par->prob_type = NONE;
   par->k_closest = -1;
   par->min_closest = 4;
   par->max_closest = 10;
   par->add_all_edges = TRUE;
   /*For now, I think this is necessary to make it work*/
   par->base_variable_selection = SOME_ARE_BASE; 
   par->use_small_graph = FALSE;
   par->colgen_strat[0] = 0;
   par->colgen_strat[1] = 0;
   par->verbosity = 0;

   lp_par->verbosity = 0;
   lp_par->branching_rule = 2;
   lp_par->branch_on_cuts = FALSE;
   lp_par->strong_branching_cand_num_max = 7;
   lp_par->strong_branching_cand_num_min = 7;
   lp_par->strong_branching_red_ratio = 0;
   lp_par->detect_tailoff  = 0;
   lp_par->child_compar_obj_tol = .01;
   lp_par->gamma = 1; /*Determines the fixed cost*/
   lp_par->tau = 0;   /*Determines the variable cost*/

   cg_par->verbosity = 0;
   cg_par->do_greedy = TRUE;
   cg_par->greedy_num_trials = 5;
   cg_par->do_extra_in_root = FALSE;
   cg_par->which_tsp_cuts = NO_TSP_CUTS;
   cg_par->which_connected_routine = BOTH;
   cg_par->max_num_cuts_in_shrink = 200;
#if defined(DIRECTED_X_VARS) && !defined(ADD_X_CUTS)
   cg_par->generate_x_cuts = TRUE;
#else
   cg_par->generate_x_cuts = FALSE;
#endif
#if defined(ADD_FLOW_VARS) && !defined(ADD_CAP_CUTS)
   cg_par->generate_cap_cuts = TRUE;
   cg_par->generate_tight_cap_cuts = TRUE;
#else
   cg_par->generate_cap_cuts = FALSE;
   cg_par->generate_tight_cap_cuts = FALSE;
#endif
   /*For minimum cut*/
   cg_par->do_mincut = TRUE;
   cg_par->update_contr_above = 0;
   cg_par->shrink_one_edges = FALSE;
   cg_par->do_extra_checking = FALSE;
   
   if (!strcmp(filename, ""))
      goto EXIT;
   
   if ((f = fopen(filename, "r")) == NULL){
      printf("VRP Readparams: file %s can't be opened\n", filename);
      exit(1); /*error check for existence of parameter file*/
   }

   while(NULL != fgets( line, LENGTH, f)){  /*read in parameter settings*/
      strcpy(key, "");
      sscanf(line, "%s%s", key, value);

      if (strcmp(key, "input_file") == 0){
	 par->infile[MAX_FILE_NAME_LENGTH] = 0;
	 strncpy(par->infile, value, MAX_FILE_NAME_LENGTH);
      }
      else if (strcmp(key, "prob_type") == 0){
	 if (strcmp("VRP", value) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = VRP;
	 }else if (strcmp("TSP", value) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = TSP;
	 }else if (strcmp("BPP", value) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = BPP;
	 }else if (strcmp("CSTP", value) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = CSTP;
	 }else if (strcmp("CTP", value) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = CTP;
	 }else{
	    fprintf(stderr, "Unknown problem type!\n");
	    exit(1);
	 }
      }
      else if (strcmp(key, "k_closest") == 0){
	 READ_INT_PAR(par->k_closest);
      }
      else if (strcmp(key, "k_closest_minimum") == 0){
	 READ_INT_PAR(par->min_closest);
      }
      else if (strcmp(key, "k_closest_maximum") == 0){
	 READ_INT_PAR(par->max_closest);
      }
      else if (strcmp(key, "add_all_edges") == 0){
	 READ_INT_PAR(par->add_all_edges);
      }
      else if (strcmp(key, "base_variable_selection") == 0){
	 READ_INT_PAR(par->base_variable_selection);
      }
      else if (strcmp(key, "use_small_graph") == 0){
	 READ_INT_PAR(par->use_small_graph);
	 if (par->use_small_graph){
	    if (fgets( line, LENGTH, f) == NULL){
	       printf("No small graph file!/n/n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "small_graph_file_name") != 0){
	       printf("Need small_graph_file_name next!!!/n/n");
	       exit(1);
	    }
	    strcpy(par->small_graph_file, value);
	 }
      }
      else if (strcmp(key, "colgen_in_first_phase") == 0 ||
	       strcmp(key, "TM_colgen_in_first_phase") == 0){
	 READ_INT_PAR(par->colgen_strat[0]);
      }
      else if (strcmp(key, "colgen_in_second_phase") == 0 ||
	       strcmp(key, "TM_colgen_in_second_phase") == 0){
	 READ_INT_PAR(par->colgen_strat[1]);
      }
      else if (strcmp(key, "colgen_in_first_phase_str") == 0 ||
	       strcmp(key, "TM_colgen_in_first_phase_str") == 0){
	 READ_STRINT_PAR(par->colgen_strat[0],
			 colgen_str, COLGEN_STR_SIZE, value);
      }
      else if (strcmp(key, "colgen_in_second_phase_str") == 0 ||
	       strcmp(key, "TM_colgen_in_second_phase_str") == 0){
	 READ_STRINT_PAR(par->colgen_strat[1],
			 colgen_str, COLGEN_STR_SIZE, value);
      }
      else if (strcmp(key, "numroutes") == 0){
	 READ_INT_PAR(j);
	 vrp->numroutes = j;
      }

      /************************ lp parameters *******************************/
      else if (strcmp(key, "branching_rule") == 0){
	 READ_INT_PAR(lp_par->branching_rule);
      }
      else if (strcmp(key, "branch_on_cuts") == 0){
	 READ_INT_PAR(lp_par->branch_on_cuts);
      }
      else if (strcmp(key, "strong_branching_cand_num") == 0){
	 READ_INT_PAR(lp_par->strong_branching_cand_num_max);
	 lp_par->strong_branching_cand_num_min =
	    lp_par->strong_branching_cand_num_max;
	 lp_par->strong_branching_red_ratio = 0;
      }
      else if (strcmp(key, "strong_branching_cand_num_min") == 0){
	 READ_INT_PAR(lp_par->strong_branching_cand_num_min);
      }
      else if (strcmp(key, "strong_branching_cand_num_max") == 0){
	 READ_INT_PAR(lp_par->strong_branching_cand_num_max);
      }
      else if (strcmp(key, "strong_branching_red_ratio") == 0){
	 READ_INT_PAR(lp_par->strong_branching_red_ratio);
      }
      else if (strcmp(key, "child_compar_obj_tol") == 0){
	 READ_FLOAT_PAR(lp_par->child_compar_obj_tol);
      }
      else if (strcmp(key, "detect_tailoff") == 0){
	 READ_INT_PAR(lp_par->detect_tailoff);
      }
      else if (strcmp(key, "gamma") == 0){
	 READ_DBL_PAR(lp_par->gamma);
      }
      else if (strcmp(key, "tau") == 0){
	 READ_DBL_PAR(lp_par->tau);
	 cg_par->tau = lp_par->tau;
      }

      /************************* cutgen parameters ***************************/

      else if (strcmp(key, "verbosity") == 0){
	 READ_INT_PAR(par->verbosity);
	 lp_par->verbosity = cg_par->verbosity = par->verbosity;
      }
      else if (strcmp(key, "do_greedy") == 0){
	 READ_INT_PAR(cg_par->do_greedy);
      }
      else if (strcmp(key, "greedy_num_trials") == 0){
	 READ_INT_PAR(cg_par->greedy_num_trials);
      }
      else if (strcmp(key, "do_extra_in_root") == 0){
	 READ_INT_PAR(cg_par->do_extra_in_root);
      }
      else if (strcmp(key, "which_tsp_cuts") == 0){
	 READ_INT_PAR(cg_par->which_tsp_cuts);
      }
#if defined(CHECK_CUT_VALIDITY) || defined(TRACE_PATH)
      else if (strcmp(key, "feasible_solution_edges") == 0){
	 READ_INT_PAR(vrp->feas_sol_size);
	 if (vrp->feas_sol_size){
	    int cur_node, prev_node = 0;
	    char value1[10], value2[10];
	    
	    vrp->feas_sol = (int *)calloc(vrp->feas_sol_size, sizeof(int));
	    for (i = 0; i < vrp->feas_sol_size; i++){
	       if (!fgets( line, LENGTH, f)){
		  fprintf(stderr,
			  "\nVrp I/O: error reading in feasible solution\n\n");
		  exit(1);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s%s", key, value1, value2);
	       if (strcmp(key, "edge")){
		  fprintf(stderr,
			  "\nVrp I/O: error reading in feasible solution\n\n");
		  exit(1);
	       }
	       if (sscanf(value1, "%i", &prev_node) != 1){
		  fprintf(stderr,
			  "\nVrp I/O: error reading in feasible solution %s\n\n",
			  key);
		  exit(1);
	       }
	       if (sscanf(value2, "%i", &cur_node) != 1){
		  fprintf(stderr,
			  "\nVrp I/O: error reading in feasible solution %s\n\n",
			  key);
		  exit(1);
	       }
	       vrp->feas_sol[i] = INDEX(prev_node, cur_node);
	    }
	 }
      }
      else if (strcmp(key, "feasible_solution_nodes") == 0){
	 READ_INT_PAR(vrp->feas_sol_size);
	 if (vrp->feas_sol_size){
	    int cur_node, prev_node = 0;
	    
	    vrp->feas_sol = (int *)calloc(vrp->feas_sol_size, sizeof(int));
	    for (i=0; i<vrp->feas_sol_size; i++){
	       if (!fgets( line, LENGTH, f)){
		  fprintf(stderr,
			  "\nVrp I/O: error reading in feasible solution\n\n");
		  exit(1);
	       }
	       sscanf(line, "%s", value);
	       if (sscanf(value, "%i", &cur_node) != 1){
		  fprintf(stderr,
			"\nVrp I/O: error reading in feasible solution %s\n\n",
			  key);
		  exit(1);
	       }else{
		  vrp->feas_sol[i] = INDEX(prev_node, cur_node);
		  prev_node = cur_node;
	       }
	    }
	 }
      }
#endif
      else if (strcmp(key, "which_connected_routine") == 0){
	 READ_INT_PAR(cg_par->which_connected_routine);
      }
      else if (strcmp(key, "max_num_cuts_in_shrink") == 0){
	 READ_INT_PAR(cg_par->max_num_cuts_in_shrink);
      }
      else if (strcmp(key, "generate_x_cuts") == 0){
	 READ_INT_PAR(cg_par->generate_x_cuts);
      }
      else if (strcmp(key, "generate_cap_cuts") == 0){
	 READ_INT_PAR(cg_par->generate_cap_cuts);
      }
      else if (strcmp(key, "generate_tight_cap_cuts") == 0){
	 READ_INT_PAR(cg_par->generate_tight_cap_cuts);
      }
      /*For minimum cut*/
      else if (strcmp(key, "do_mincut") == 0){
	 READ_INT_PAR(cg_par->do_mincut);
      }
      else if (strcmp(key, "update_contr_above") == 0){
	 READ_INT_PAR(cg_par->update_contr_above);
      }
      else if (strcmp(key, "shrink_one_edges") == 0){
	 READ_INT_PAR(cg_par->shrink_one_edges);
      }
      else if (strcmp(key, "do_extra_checking") == 0){
	 READ_INT_PAR(cg_par->do_extra_checking);
      }
   }

EXIT:
   
   for (i = 1; i < argc; i++){
      sscanf(argv[i], "%c %c", &tmp, &c);
      if (tmp != '-')
	 continue;
      switch (c) {
       case 'D':
	 sscanf(argv[++i], "%i", &lp_par->verbosity);
	 break;
       case 'H':
	 user_usage();
	 exit(0);
       case 'E':
	 par->add_all_edges = FALSE;
	 break;
       case 'S':
	 strncpy(par->small_graph_file, argv[++i], MAX_FILE_NAME_LENGTH);
	 par->use_small_graph = LOAD_SMALL_GRAPH;
	 break;
       case 'F':
	 strncpy(par->infile, argv[++i], MAX_FILE_NAME_LENGTH);
	 break;
       case 'B':
	 sscanf(argv[++i], "%i", &lp_par->branching_rule);
	 break;
       case 'V':
	 sscanf(argv[++i], "%i", &par->base_variable_selection);
	 break;
       case 'K':
	 sscanf(argv[++i], "%i", &par->k_closest);
	 break;
       case 'N':
	 sscanf(argv[++i], "%i", &vrp->numroutes);
	 break;
       case 'M':
	 cg_par->do_mincut = FALSE;
	 break;
       case 'X':
	 sscanf(argv[++i], "%i", &cg_par->generate_x_cuts);
	 break;
       case 'Y':
	 sscanf(argv[++i], "%i", &cg_par->generate_cap_cuts);
	 break;
       case 'Z':
	 sscanf(argv[++i], "%i", &cg_par->generate_tight_cap_cuts);
	 break;
       case 'G':
	 sscanf(argv[++i], "%lf", &lp_par->tau);
	 cg_par->tau = lp_par->tau;
	 break;
       case 'T':
	 i++;
	 if (strcmp("VRP", argv[i]) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = VRP;
	 }else if (strcmp("TSP", argv[i]) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = TSP;
	 }else if (strcmp("BPP", argv[i]) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = BPP;
	 }else if (strcmp("CSTP", argv[i]) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = CSTP;
	 }else if (strcmp("CTP", argv[i]) == 0){
	    par->prob_type = cg_par->prob_type = 
		 vrp->lp_par.prob_type = CTP;
	 }else{
	    fprintf(stderr, "Unknown problem type!\n");
	    exit(1);
	 }
	 break;

       case 'C':
	 sscanf(argv[++i], "%i", &vrp->capacity);
	 break;
      };
   }

   if (f)
      fclose(f);
}
