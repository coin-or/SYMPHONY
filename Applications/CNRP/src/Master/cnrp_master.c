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
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

/* SYMPHONY include files */
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#include "master.h"
/*___END_EXPERIMENTAL_SECTION___*/
#include "BB_macros.h"
#include "BB_constants.h"
#include "proccomm.h"
#include "qsortucb.h"
#include "dg_params.h"
#include "master_u.h"

/* CNRP include files */
#include "cnrp_const.h"
#include "cnrp_types.h"
#include "cnrp_io.h"
#include "compute_cost.h"
#include "cnrp_master_functions.h"
#include "cnrp_dg_functions.h"
#include "cnrp_macros.h"
#include "small_graph.h"
#ifdef COMPILE_IN_TM
#ifdef COMPILE_IN_LP
#include "cnrp_lp.h"
#ifdef COMPILE_IN_CG
#include "cnrp_cg.h"
#endif
#ifdef COMPILE_IN_CP
#include "cnrp_cp.h"
#endif
#endif
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions for the master process.
\*===========================================================================*/

void user_usage(void){
         printf("master [ -HEP ] [ -S file ] [ -F file ] [ -B rule ]\n"
		"[ -V sel ] [ -K closest ] [ -N routes ] [ -C capacity ]\n"
		"[ -D level ] [ -M ] [ -X toggle ] [ -Y toggle ] \n"
		"[ -Z toggle ] [-G tau] \n"
		"\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n"
		"\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\n",
		"-H: help",
		"-E: use sparse edge set",
		"-D level: verbosity level for displaying LP solutions",
		"-P type: specify problem type",
		"-S file: load sparse graph from 'file'",
		"-F file: problem data is in 'file'",
		"-B i: which candidates to check in strong branching",
		"-V i: how to construct the base set of variables",
		"-K k: use 'k' closest edges to build sparse graph",
		"-N n: use 'n' routes",
		"-M  : use min cut subroutine",
		"-C c: use capacity 'c'",
		"-X t: toggles generation of X cuts",
		"-Y t: toggles generation of capacity cuts",
		"-Z t: toggles generation of tight capacity cuts",
		"-G t: set tau to 't'");
}

/*===========================================================================*\
 * Initialize user-defined data structures. In this case, I store all
 * problem-specific data such as the location of the customers, edge costs,
 * etc. in this data-structure.
\*===========================================================================*/

int user_initialize(void **user)
{
   vrp_problem *vrp = (vrp_problem *) calloc(1, sizeof(vrp_problem));

   *user = vrp;

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * In this function, I set up the user parameters. The first step is to cast
 * the void pointer in order to access my data. In the readparams() function,
 * I read in parameters from the parameter file given on the command line.
\*===========================================================================*/

int user_readparams(void *user, char *filename, int argc, char **argv)
{
   vrp_problem *vrp = (vrp_problem *)user;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   problem *p = get_problem_ptr(FALSE);

   p->par.tm_par.granularity = p->par.lp_par.granularity = .00001;
   p->par.lp_par.problem_type = INTEGER_PROBLEM;
   strcpy(p->par.dg_par.source_path, "/home/tkr/BlackBox/DrawGraph/IGD_1.0/");
   /*___END_EXPERIMENTAL_SECTION___*/

   vrp_readparams(vrp, filename, argc, argv);

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * After I've read in the parameters, I can now read in the data file, whose
 * name was given in the parameter file. This file contains instance data.
\*===========================================================================*/

int user_io(void *user)
{
   vrp_problem *vrp = (vrp_problem *)user;

   vrp_io(vrp, vrp->par.infile);

   return(USER_SUCCESS);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * Here is where the heuristics are performed and an upper bound is calculated.
 * An upper bound can also be specified in the parameter file. The
 * other thing I do in this routine is build up a graph of the
 * cheapest k edges adjacent to the each node plus any edges chosen
 * during the heuristics to comprise my base set later.
\*===========================================================================*/

int user_start_heurs(void *user, double *ub, double *ub_estimate)
{
   vrp_problem *vrp = (vrp_problem *)user;

   if (*ub > 0){
      vrp->cur_tour->cost = (int) (*ub);
   }else{
      vrp->cur_tour->cost = MAXINT;
   }

   vrp->cur_tour->numroutes = vrp->numroutes;
   
   if (vrp->par.use_small_graph == LOAD_SMALL_GRAPH){
      read_small_graph(vrp);
      if (*ub <= 0 && vrp->cur_tour->cost > 0)
	 *ub = (int)(vrp->cur_tour->cost);
      vrp->numroutes = vrp->cur_tour->numroutes;
   }

#if 0
   if(vrp->par.prob_tpye == BPP)
      *ub = 1;
#endif
   
   if (!vrp->numroutes && vrp->par.prob_type == VRP){
      printf("\nError: Number of trucks not specified or computed "
	     "for VRP\n\n");
      exit(1);
   }
   
   if (vrp->numroutes > 1){
      printf("NUMBER OF TRUCKS: \t%i\n", vrp->numroutes);
      printf("TIGHTNESS: \t\t%.2f\n",
     (double)vrp->demand[0]/((double)vrp->capacity*(double)vrp->numroutes));
   }
   
   if (*ub > 0 && !(vrp->par.prob_type == BPP))
      printf("INITIAL UPPER BOUND: \t%i\n\n", (int)(*ub));
   else if (!(vrp->par.prob_type == BPP))
      printf("INITIAL UPPER BOUND: \tNone\n\n");
   else
      printf("\n\n");
   
   /* Selects the cheapest edges adjacent to each node for the base set */

   if (vrp->par.use_small_graph == SAVE_SMALL_GRAPH){
      if (!vrp->g) make_small_graph(vrp, 0);
      save_small_graph(vrp);
   }else if (!vrp->g){
      make_small_graph(vrp, 0);
   }

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * If graph drawing will be use, the user must initialize the drawing
 * window here.
\*===========================================================================*/

int user_init_draw_graph(void *user, int dg_id)
{
#ifndef WIN32   /* FIXME : None of this works in Windows */
   vrp_problem *vrp = (vrp_problem *)user;
   int s_bufid;
      
   if (!(vrp->posx && vrp->posy)) return(USER_SUCCESS);
   if ( (vrp->dg_id = dg_id) ){
      int i, zero = 0, eight = 0x08;
      char node_place[MAX_NAME_LENGTH] = {"node_placement"};
      char weight[5];
      int *posx = vrp->posx, *posy = vrp->posy;
      int minx=MAXINT, miny=MAXINT, maxx=-MAXINT, maxy=-MAXINT, xx, yy;
      int width = 1000, height = 700;
#if 0
      int width=p->par.dg_par.canvas_width, height=p->par.dg_par.canvas_height;
#endif
      double mult;

      for (i = vrp->vertnum - 1; i >= 0; i--){
	 if (posx[i] < minx) minx = posx[i];
	 if (posx[i] > maxx) maxx = posx[i];
	 if (posy[i] < miny) miny = posy[i];
	 if (posy[i] > maxy) maxy = posy[i];
      }
      xx = maxx - minx;
      yy = maxy - miny;
      mult = (int) MIN((width - 20.0)/xx, (height-20.0)/yy);
      width = (int) (xx * mult + 30);
      height = (int) (yy * mult + 30);
      for (i = vrp->vertnum-1; i >= 0; i--){
	 posx[i] = (int) ((posx[i] - minx) * mult + 10);
	 posy[i] = (int) ((maxy - posy[i]) * mult + 10);
      }

      init_window(dg_id, node_place, width, height);
      /* Now pack the placement of the nodes of the graph */
      s_bufid = init_send(DataInPlace);
      send_str(node_place);
      send_int_array(&vrp->vertnum, 1);
      for (i = 0; i < vrp->vertnum; i++){
	 send_int_array(&i, 1);
	 send_int_array(posx + i, 1);
	 send_int_array(posy + i, 1);
	 send_int_array(&eight, 1);
	 sprintf(weight, "%i", vrp->demand[i]);
	 send_str(weight);
      }
      /* No edges are passed to the default graph */
      send_int_array(&zero, 1);
      send_msg(dg_id, CTOI_SET_GRAPH);
      freebuf(s_bufid);
      
      display_graph(dg_id, node_place);
   }
#endif

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * In this routine, I build the initial edge set for the root. There are
 * several things going on here. First, there is a user-defined parameter
 * defining whether or not to just go ahead and add all variables to the
 * problem up front (vrp->par.add_all_edges). Currently, this seems to be the
 * best option since the problems are small anyway. Further, I am doing some
 * preprocessing here by eliminating edges for which the sum of the demands of
 * their endpoints is greater than the capacity since these edges cannot
 * be in any feasible solution.
 *
 * Notice that there are several options programmed for which set
 * of edges should be in the base set. The
 * base constraints are just the degree constraints from the IP
 * formulation. These do not have to be specified explicitly, just the
 * number of them given.
\*===========================================================================*/

int user_initialize_root_node(void *user, int *basevarnum, int **basevars,
			      int *basecutnum, int *extravarnum,
			      int **extravars, char *obj_sense,
			      double *obj_offset, char ***colnames,
			      int *colgen_strat)
{
   vrp_problem *vrp = (vrp_problem *)user;
   int base_varnum = 0, i, j, k, l;
   int zero_varnum, *zero_vars;
   int *edges;
   int vertnum = vrp->vertnum;
   
#ifdef DIRECTED_X_VARS
   /*whether or not we will have the out-degree constraints*/
   char od_const = (vrp->par.prob_type == TSP || vrp->par.prob_type == VRP ||
		    vrp->par.prob_type == BPP);
   char d_x_vars = TRUE;
#else
   char od_const = FALSE;
   char d_x_vars = FALSE;
#endif
#if defined(ADD_CAP_CUTS) || defined(ADD_X_CUTS) 
   int total_edgenum = vertnum*(vertnum - 1)/2;
#endif
#ifdef ADD_FLOW_VARS
   int v0, v1;
   double flow_capacity;
#ifdef DIRECTED_X_VARS
   flow_capacity = (double) vrp->capacity;
#else
   if (vrp->par.prob_type == CSTP || vrp->par.prob_type == CTP)
      flow_capacity = (double) vrp->capacity;
   else
      flow_capacity = ((double)vrp->capacity)/2;
#endif
#endif
   
#ifdef ADD_CAP_CUTS 
   *basecutnum = (2 + od_const)*vertnum - 1 + 2*total_edgenum;
#elif defined(ADD_FLOW_VARS)
   *basecutnum = (2 + od_const)*vertnum - 1;
#else
   *basecutnum = (1 + od_const)*vertnum;
#endif
#ifdef ADD_X_CUTS
   *basecutnum += total_edgenum;
#endif

   switch(vrp->par.base_variable_selection){
    case SOME_ARE_BASE:
      if (vrp->par.add_all_edges == FALSE)
	 /*If we are not adding all the edges, then really EVERYTHING_IS_BASE*/
	 vrp->par.base_variable_selection = EVERYTHING_IS_BASE;
      else /*Otherwise, all we need to do is set this and then fall through --
	     the remaining edges get added in user_create_root()*/
	 vrp->par.add_all_edges = FALSE;


    case EVERYTHING_IS_BASE:
      *basevars = create_edge_list(vrp, &base_varnum, CHEAP_EDGES);
      
      *basevars = (int *) realloc((char *)(*basevars), base_varnum * ISIZE);
      *basevarnum = base_varnum;
      
      /*First, set base_varnum equal to the number of edges in the base set*/
#ifdef ADD_FLOW_VARS
      base_varnum = (*basevarnum)/(3+d_x_vars);
#else
      base_varnum = (*basevarnum)/(1+d_x_vars);
#endif
      break;

    case EVERYTHING_IS_EXTRA:
      *basevarnum = 0;
      break;
   }

   if (!vrp->par.colgen_strat[0]){
      if (vrp->par.add_all_edges ||
	  vrp->par.base_variable_selection == SOME_ARE_BASE){
	 colgen_strat[0]=(FATHOM__DO_NOT_GENERATE_COLS__DISCARD |
			  BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
      }else{
	 colgen_strat[0] = (FATHOM__DO_NOT_GENERATE_COLS__SEND  |
			    BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
      }
   }else{
      colgen_strat[0] = vrp->par.colgen_strat[0];
   }
   if (!vrp->par.colgen_strat[1]){
      if (vrp->par.add_all_edges ||
	  vrp->par.base_variable_selection == SOME_ARE_BASE){
	 colgen_strat[1]=(FATHOM__DO_NOT_GENERATE_COLS__DISCARD |
			  BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
      }else{
	 colgen_strat[1] = (FATHOM__GENERATE_COLS__RESOLVE  |
			    BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
      }
   }else{
      colgen_strat[1] = vrp->par.colgen_strat[1];
   }
   
#if 0
   if (vrp->par.prob_tpye == BPP){
      for (i = 0; i < *basevarnum; i++){
	 vrp->dist.cost[(*basevars)[i]] = 10;
      }
   }
#endif
       
   /*create the edge list (we assume a complete graph) The edge is set to
     (0,0) in the edge list if it was eliminated in preprocessing*/
   edges = vrp->edges = (int *) calloc (vertnum*(vertnum-1), sizeof(int));
   zero_varnum = vrp->zero_varnum;
   zero_vars = vrp->zero_vars;
   for (i = 1, k = 0, l = 0; i < vertnum; i++){
      for (j = 0; j < i; j++){
	 if (l < zero_varnum && k == zero_vars[l]){
	    /*This is one of the zero edges*/
	    edges[2*k] = edges[2*k+1] = 0;
	    l++;
	    k++;
	    continue;
	 }
	 edges[2*k] = j;
	 edges[2*k+1] = i;
	 k++;
      }
   }

   switch(vrp->par.base_variable_selection){
    case EVERYTHING_IS_EXTRA:

      *extravars  = create_edge_list(vrp, extravarnum, CHEAP_EDGES);
      
      break;

    case SOME_ARE_BASE:
      
      vrp->par.add_all_edges = TRUE; /*We turned this off in user_set_base()
				       -- now we need to turn it back on*/

      *extravars  = create_edge_list(vrp, extravarnum, REMAINING_EDGES);

      break;

    case EVERYTHING_IS_BASE:

      break;
   }

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * In my case, a feasible solution is specified most compactly by
 * essentially a permutation of the customers along with routes numbers,
 * specifying the order of the customers on their routes. This is just
 * sent as a character array which then gets cast to an array of
 * structures, one for each customers specifying the route number and
 * the next customer on the route.
\*===========================================================================*/

int user_receive_feasible_solution(void *user, int msgtag, double cost,
				   int numvars, int *indices, double *values)
{
   vrp_problem *vrp = (vrp_problem *)user;

   if (vrp->par.prob_type == TSP || vrp->par.prob_type == VRP ||
       vrp->par.prob_type == BPP)
      receive_char_array((char *)vrp->cur_tour->tour,
			 vrp->vertnum*sizeof(_node));
   else
      receive_int_array(vrp->cur_sol_tree, vrp->vertnum);

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the LP process. Notice that
 * there are two cases to deal with. If the LP or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. 
\*===========================================================================*/

int user_send_lp_data(void *user, void **user_lp)
{
   vrp_problem *vrp = (vrp_problem *)user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* This is the case when we are copying data directly because
      the LP is not running separately. This code should be virtually
      identical to that of user_receive_lp_data() in the LP process.*/
   
   vrp_spec *vrp_lp = (vrp_spec *) calloc(1, sizeof(vrp_spec));
   int zero_varnum = vrp->zero_varnum;
   int *zero_vars = vrp->zero_vars;
   int vertnum, i, j, k, l;

   *user_lp = (void *)vrp_lp;
   
   vrp_lp->par = vrp->lp_par;
   vrp_lp->window = vrp->dg_id;
   vrp_lp->numroutes = vrp->numroutes;
   vertnum = vrp_lp->vertnum = vrp->vertnum;
   vrp_lp->edges = vrp->edges;
   vrp_lp->demand = vrp->demand;
   vrp_lp->capacity = vrp->capacity;
   vrp_lp->costs = vrp->dist.cost;

   if (vrp->par.prob_type == VRP || vrp->par.prob_type == TSP ||
       vrp->par.prob_type == BPP){
      vrp_lp->cur_sol = (_node *) calloc (vrp->vertnum, sizeof(_node));
   }else{
      vrp_lp->cur_sol_tree = (int *) calloc (vrp->vertnum - 1, ISIZE);
   }
/*__BEGIN_EXPERIMENTAL_SECTION__*/
   if (vrp_lp->window){
      copy_node_set(vrp_lp->window, TRUE, (char *)"Weighted solution");
      copy_node_set(vrp_lp->window, TRUE, (char *)"Flow solution");
   }
/*___END_EXPERIMENTAL_SECTION___*/
   
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_lp_data() in the LP process */
   
   send_char_array((char *)(&vrp->lp_par), sizeof(lp_user_params));
   send_int_array(&vrp->dg_id, 1);
   send_int_array(&vrp->numroutes, 1);
   send_int_array(&vrp->vertnum, 1);
   send_int_array(vrp->demand, vrp->vertnum);
   send_int_array(&vrp->capacity, 1);
   send_int_array(vrp->dist.cost, vrp->edgenum);
   send_int_array(&vrp->zero_varnum, 1);
   if (vrp->zero_varnum){
      send_int_array(vrp->zero_vars, vrp->zero_varnum);
   }
#endif

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the CG process. Notice that
 * there are two cases to deal with. If the CG, LP, or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. 
\*===========================================================================*/

int user_send_cg_data(void *user, void **user_cg)
{
   vrp_problem *vrp = (vrp_problem *)user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CG)
   /* This is is the case when we are copying data directly because
      the CG is not running separately. This code should be virtually
      identical to that of user_receive_cg_data() in the CG process.*/
   
   cg_vrp_spec *vrp_cg = (cg_vrp_spec *) malloc (sizeof(cg_vrp_spec));
   int edgenum, vertnum, i, j, k;
   
   *user_cg = (void *)vrp_cg;

   vrp_cg->par = vrp->cg_par;
   vrp_cg->numroutes = vrp->numroutes;
   vertnum = vrp_cg->vertnum = vrp->vertnum;
   vrp_cg->demand = vrp->demand;
   vrp_cg->capacity = vrp->capacity;
   vrp_cg->dg_id = vrp->dg_id;
   
   edgenum = vrp->vertnum*(vrp->vertnum-1)/2;
      
   vrp_cg->in_set = (char *) calloc(vrp->vertnum, sizeof(char));
   vrp_cg->ref = (int *) malloc(vrp->vertnum*sizeof(int));
   vrp_cg->new_demand = (int *) malloc(vrp->vertnum*sizeof(int));
   vrp_cg->cut_val = (double *) calloc(vrp->vertnum, sizeof(double));
   vrp_cg->cut_list = (char *) malloc(((vrp->vertnum >> DELETE_POWER)+1)*
				   (vrp->cg_par.max_num_cuts_in_shrink + 1)*
				   sizeof(char));

   vrp_cg->edges = (int *) calloc (2*edgenum, sizeof(int));
   
   /*create the edge list (we assume a complete graph)*/
   for (i = 1, k = 0; i < vertnum; i++){
      for (j = 0; j < i; j++){
	 vrp_cg->edges[2*k] = j;
	 vrp_cg->edges[2*k+1] = i;
	 k++;
      }
   }

#ifdef CHECK_CUT_VALIDITY
   if ((vrp_cg->feas_sol_size = vrp->feas_sol_size)){
      vrp_cg->feas_sol = vrp->feas_sol;
   }
#endif
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cg_data() in the CG process */
   
   send_char_array((char *)&vrp->cg_par, sizeof(cg_user_params));
   send_int_array(&vrp->dg_id, 1);
   send_int_array(&vrp->numroutes, 1);
   send_int_array(&vrp->vertnum, 1);
   send_int_array(vrp->demand, vrp->vertnum);
   send_int_array(&vrp->capacity, 1);
#ifdef CHECK_CUT_VALIDITY
   send_int_array(&vrp->feas_sol_size, 1);
   if (vrp->feas_sol_size){
      send_int_array(vrp->feas_sol, vrp->feas_sol_size);
   }
#endif
#endif

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the CP process. Notice that
 * there are two cases to deal with. If the CP, LP, or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. 
\*===========================================================================*/

int user_send_cp_data(void *user, void **user_cp)
{
   vrp_problem *vrp = (vrp_problem *)user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CP)
   /* This is is the case when we are copying data directly because
      the LP is not running separately. This code should be virtually
      identical to that of user_receive_cp_data() in the CP process.*/
   
   vrp_spec_cp *vrp_cp = (vrp_spec_cp *) malloc (sizeof(vrp_spec_cp));
   int i, j, k;

   vrp_cp->vertnum = vrp->vertnum;

   *user_cp = (void *)vrp_cp;

   vrp_cp->edgenum = vrp_cp->vertnum*(vrp_cp->vertnum-1)/2 + vrp_cp->vertnum-1;
   vrp_cp->edges = (int *) calloc ((int)2*vrp_cp->edgenum, sizeof(int));
     
   /* create the edge list (we assume a complete graph) */
   for (i = 1, k = 0; i < vrp_cp->vertnum; i++){
      for (j = 0; j < i; j++){
	 vrp_cp->edges[2*k] = j;
	 vrp_cp->edges[2*k+1] = i;
	 k++;
      }
   }

   /* now add the duplicate copies of the depot edges to allow for
      routes with one customer */
   for (i = 1; i < vrp_cp->vertnum; i++){
      vrp_cp->edges[2*k] = 0;
      vrp_cp->edges[2*k+1] = i;
      k++;
   }
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cp_data() in the CP process */
   
   send_int_array(&vrp->vertnum, 1);
#endif

   return(USER_SUCCESS);
}

/*__BEGIN_EXPERIMENTAL_SECTION__*/
/*===========================================================================*/

int user_send_sp_data(void *user)
{
   return(USER_SUCCESS);
}

/*___END_EXPERIMENTAL_SECTION___*/
/*===========================================================================*/

/*===========================================================================*\
 * Generally, this function is not needed but you might find some use
 * for it. Someone did :).
\*===========================================================================*/

int user_process_own_messages(void *user, int msgtag)
{
   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the user's chance to display the solution in whatever
 * manner desired. In my case, I can display a text version, which is
 * just a list of the customers on each route, or use the interactive
 * Graph Drawing application to graphically display the solution.
\*===========================================================================*/

int user_display_solution(void *user, double lpetol, int varnum, int *indices,
			  double *values, double objval)
{
   vrp_problem *vrp = (vrp_problem *)user;
   _node *tour = vrp->cur_tour->tour;
   int cur_vert = 0, prev_vert = 0, cur_route, i, count;
   elist *cur_route_start = NULL;
   edge *edge_data;
   vertex *verts;
   double fixed_cost = 0.0, variable_cost = 0.0;
   int window = vrp->dg_id;
   int vertnum = vrp->vertnum, v0, v1;
   int total_edgenum =  vertnum*(vertnum-1)/2;
   network *n;
   
#if 0
   if (tour && vrp->cur_tour->cost > (int) objval){
      node = tour[0].next;
      
      printf("\nSolution Found:\n");
      if (tour[0].route == 1)
	 printf("\n0 ");
      while (node != 0){
	 if (tour[prev_node].route != tour[node].route){
	    printf("\nRoute #%i: ", tour[node].route);
	    count = 0;
	 }
	 printf("%i ", node);
	 count++;
	 if (count > 15){
	    printf("\n");
	    count = 0;
	 }
	 prev_node = node;
	 node = tour[node].next;
      }
      printf("\n\n");
      
      if (window){
	 char name[MAX_NAME_LENGTH] = {"feas_solution"};
	 disp_vrp_tour(window, TRUE, name, tour, vrp->vertnum, vrp->numroutes,
		       CTOI_WAIT_FOR_CLICK_AND_REPORT);
      }
   }

   if (tree){
      printf("\nSolution Found:\n");
      printf("Edge List:\n");
      for (i = 0; i < vertnum - 1; i++){
	 BOTH_ENDS(tree[i], &v0, &v1);
	 printf("%i %i\n", v0, v1);
      }
      printf("\n\n");

      if (window){
	 char name[MAX_NAME_LENGTH] = {"feas_solution"};
	 copy_node_set(window, TRUE, name);
	 draw_edge_set_from_userind(window, name, vertnum - 1, tree);
	 display_graph(window, name);
	 wait_for_click(window, name, CTOI_WAIT_FOR_CLICK_NO_REPORT);
      }
   }
#endif

   /*Otherwise, construct the solution from scratch*/

#ifdef ADD_FLOW_VARS
#ifndef ADD_CAP_CUTS
      n = create_flow_net(indices, values, varnum, lpetol, vrp->edges,
			  vrp->demand, vertnum);
#else
#ifdef DIRECTED_X_VARS
   for (i = 0; i < varnum && indices[i] < 2*total_edgenum; i++);
#else
   for (i = 0; i < varnum && indices[i] < total_edgenum; i++);
#endif   
   varnum = i;
   
   n = create_net(indices, values, varnum, lpetol, vrp->edges, vrp->demand,
		  vertnum);
#endif   
#else
   n = create_net(indices, values, varnum, lpetol, vrp->edges, vrp->demand,
		  vertnum);
#endif

   for (i = 0; i < n->edgenum; i++){
      fixed_cost += vrp->dist.cost[INDEX(n->edges[i].v0, n->edges[i].v1)];
#ifdef ADD_FLOW_VARS
      variable_cost += (n->edges[i].flow1+n->edges[i].flow2)*
	 vrp->dist.cost[INDEX(n->edges[i].v0, n->edges[i].v1)];
#endif
   }
   vrp->fixed_cost = fixed_cost;
   vrp->variable_cost = variable_cost;
   
   printf("\nSolution Found:\n");
#ifdef ADD_FLOW_VARS
   printf("Solution Fixed Cost: %.1f\n", fixed_cost);
   printf("Solution Variable Cost: %.1f\n", variable_cost);
#else
   printf("Solution Cost: %.0f\n", fixed_cost);
#endif
   
   if (vrp->par.prob_type == TSP || vrp->par.prob_type == VRP ||
       vrp->par.prob_type == BPP){ 

      verts = n->verts;
   
     /*construct the tour corresponding to this solution vector*/
      for (cur_route_start = verts[0].first, cur_route = 1,
	      edge_data = cur_route_start->data; cur_route <= vrp->numroutes;
	   cur_route++){
	 edge_data = cur_route_start->data;
	 edge_data->scanned = TRUE;
	 cur_vert = edge_data->v1;
	 tour[prev_vert].next = cur_vert;
	 tour[cur_vert].route = cur_route;
	 prev_vert = 0;
	 while (cur_vert){
	    if (verts[cur_vert].first->other_end != prev_vert){
	       prev_vert = cur_vert;
	       edge_data = verts[cur_vert].first->data;
	       cur_vert = verts[cur_vert].first->other_end;
	    }
	    else{
	       prev_vert = cur_vert;
	       edge_data = verts[cur_vert].last->data; /*This statement
							 could possibly
							 be taken out to speed
							 things up a bit*/
	       cur_vert = verts[cur_vert].last->other_end;
	    }
	    tour[prev_vert].next = cur_vert;
	    tour[cur_vert].route = cur_route;
	 }
	 edge_data->scanned = TRUE;
	 
	 while (cur_route_start->data->scanned){
	    if (!(cur_route_start = cur_route_start->next_edge)) break;
	 }
      }
      
      /* Display the solution */
      
      cur_vert = tour[0].next;
      
      if (tour[0].route == 1)
	 printf("\n0 ");
      while (cur_vert != 0){
	 if (tour[prev_vert].route != tour[cur_vert].route){
	    printf("\nRoute #%i: ", tour[cur_vert].route);
	    count = 0;
	 }
	 printf("%i ", cur_vert);
	 count++;
	 if (count > 15){
	    printf("\n");
	    count = 0;
	 }
	 prev_vert = cur_vert;
	 cur_vert = tour[cur_vert].next;
      }
      printf("\n");
   }else{
      
      for (i = 0; i < n->edgenum; i++){
	 vrp->cur_sol_tree[i] = INDEX(n->edges[i].v0, n->edges[i].v1);
      }
      
      /* Display the solution */
      
      for (i = 0; i < n->edgenum; i++){
	 printf("%i %i\n", n->edges[i].v0, n->edges[i].v1);
      }
      printf("\n");
   }

   return(USER_SUCCESS);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * This is a debugging feature which might
 * allow you to find out why a known feasible solution is being cut off.
\*===========================================================================*/

int user_send_feas_sol(void *user, int *feas_sol_size, int **feas_sol)
{
#ifdef TRACE_PATH
   vrp_problem *vrp = (vrp_problem *)user;

   *feas_sol_size = vrp->feas_sol_size;
   *feas_sol = vrp->feas_sol;
#endif
   return(USER_SUCCESS);
}   

/*===========================================================================*/

/*===========================================================================*\
 * This function frees everything.
\*===========================================================================*/

int user_free_master(void **user)
{
   vrp_problem *vrp = (vrp_problem *)(*user);

   if (vrp->cur_tour){
      FREE(vrp->cur_tour->tour);
      FREE(vrp->cur_tour->route_info);
      FREE(vrp->cur_tour);
   }
   FREE(vrp->cur_sol_tree);
   FREE(vrp->posy);
   FREE(vrp->posx);
   FREE(vrp->dist.coordx);
   FREE(vrp->dist.coordy);
   FREE(vrp->dist.coordz);
#if !(defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP))
   FREE(vrp->dist.cost);
   FREE(vrp->edges);
   FREE(vrp->demand);
#endif
   if (vrp->g){
      FREE(vrp->g->edges);
      FREE(vrp->g);
   }
#ifdef CHECK_CUT_VALIDITY
   FREE(vrp->feas_sol);
#endif
   FREE(vrp->zero_vars);
   FREE(vrp);

   return(USER_SUCCESS);
}







