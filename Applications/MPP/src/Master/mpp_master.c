/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2002 Ted Ralphs. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "BB_constants.h"
#include "BB_macros.h"
#include "master_u.h"
#include "mpp.h"
#ifdef COMPILE_IN_TM
#ifdef COMPILE_IN_LP
/* fill these in for sequentail compilation */
#ifdef COMPILE_IN_CG
/* fill these in for sequentail compilation */
#endif
#ifdef COMPILE_IN_CP
/* fill these in for sequentail compilation */
#endif
#endif
#endif

/*===========================================================================*/

void user_usage(void){
  printf("master [ -H ] [ -F file ] \n\t%s\n\t%s\n",
	 "-H: help (user switches)",
	 "-F file: problem instance data is in 'file'");
}

/*===========================================================================*\
 * This file contains the user-written functions for the master process.
\*===========================================================================*/

/*===========================================================================*\
 * Initialize user-defined data structures. In this case, I store all
 * problem-specific data such as the location of the customers, edge costs,
 * etc. in this data-structure.
\*===========================================================================*/

int user_initialize(void **user)
{
   mpp_problem *mpp = (mpp_problem *) calloc(1, sizeof(mpp_problem));
  
   *user = mpp;
   
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Parse the user options and read in parameters from the parameter file 
 * given on the command line
\*===========================================================================*/

int user_readparams(void *user, char *filename, int argc, char **argv)
{
   FILE *f;
   mpp_problem *mpp=(mpp_problem *)user;
  
   char line[50], key[50], value[50];
   int c;
  
   strcpy(mpp->infile,argv[1]);
 /*  if ((f = fopen(filename, "r")) == NULL){
      printf("SYMPHONY: file %s can't be opened\n", filename);
      exit(1); /*error check for existence of parameter file
} 
   
   while(NULL != fgets(line, 50, f)){  /*read in parameter settings
      strcpy(key, "");
      sscanf(line, "%s%s", key, value);
   } */     

  /* while ((c = getopt(argc, argv, "H")) != -1){
      switch (c) {
       case 'H':
	 user_usage();
	 exit(0);
	 break;
      }; 
   }*/

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Read in the data file, whose name was given in the parameter file.
 * This file contains instance data.
\*===========================================================================*/



void mpp_io(mpp_problem *mpp, char *infile)
{
	int i;
    char line1[80]; 
	int holder;

  FILE *f;
	/* input will be of the following format
	//  number of nodes
	//  number of arcs
	//  number of edges
	//  a "start node" "end node" "cost"
	//    --this until all edges are declared
	//  e "start node" "end node" "cost"
	//    --this until all arcs are declared */


//sets aside memory for the number for variables
//laid out in the following way:
// ****changed the order to arcs before edges for simplicity
//vars[0] to vars [arcs-1] is the weights from start to end on arcs
//vars[arcs] to vars [arcs+edges-1] is the weights from start to end on edges
//vars[arcs+edges] to vars [2*edges+ arcs-1] is the weights from end to start on edges
//  */


/*makes sure the file exists and can be opened*/
  if (!strcmp(infile, "")){
     printf("\nMpp I/O: No problem data file specified\n\n");
     exit(1);
  }
  
  if ((f = fopen(infile, "r")) == NULL){
     fprintf(stderr, "Mpp I/O: file '%s' can't be opened\n", infile);
     exit(1);
  }

  /* input first 3 ints  as number of nodes,arcs and edges respectively */

  fgets( line1, 80, f);
  sscanf(line1,"%d",&(mpp->nodes)); /*read in number of nodes*/


  fgets( line1, 80, f);
  sscanf(line1,"%d",&(mpp->arcs)); /*read in number of arcs*/
  

  fgets( line1, 80, f);
  sscanf(line1,"%d",&(mpp->edges)); /*read in number of edges*/
  holder=(mpp->edges+mpp->arcs);

  mpp->types= (char *) malloc (holder*sizeof(char));

  mpp->cost=(int *) malloc (holder*sizeof(int));

  mpp->start_node=(int *) malloc (holder*sizeof(int));

  mpp->end_node=(int *) malloc (holder*sizeof(int));

  mpp->odd_checker=0;

for(i = 0; i< mpp->edges + mpp->arcs ;i++)
	{ 
	 fgets(line1, 80, f);
	 sscanf(line1,"%c %d %d %d",&(mpp->types[i]),&(mpp->start_node[i]),&(mpp->end_node[i]),&(mpp->cost[i])); /*read in next type*/

	/* sscanf(line1,"%d", &(mpp->start_node[i])); /*read in next startnode*/
   holder=mpp->start_node[i];
 /* 	 sscanf(line1,"%d",&(mpp->end_node[i])); /*read in next endnode*/
    holder=mpp->end_node[i];
/* 	 sscanf(line1,"%d",&(mpp->cost[i])); /*read in next cost*/
 holder=mpp->cost[i]; 
	}
 
}
   
/*===========================================================================*/

/*===========================================================================*\
 * Here is where the heuristics are performed and an upper bound is calculated.
 * An upper bound can also be specified in the parameter file. 
\*===========================================================================*/

int user_start_heurs(void *user, double *ub, double *ub_estimate)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * If graph drawing will be use, the user must initialize the drawing
 * window here.
\*===========================================================================*/
int user_io(void *user)
{
   mpp_problem *mpp = (mpp_problem *)user;
   
   mpp_io(mpp, mpp->infile);

   return(USER_NO_PP);
}
   

int user_init_draw_graph(void *user, int dg_id)
{
   return(USER_NO_PP);
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

int user_set_base(void *user, int *basevarnum, int **basevars, double **lb,
		  double **ub, int *basecutnum, int *colgen_strat)
	{
	mpp_problem *mpp=(mpp_problem *)user;
int i;
int base_varnum=0;
	/*
	RALPHS stuff
	*basevarnum = 0;
	*basevars = NULL;
	*basecutnum = 0;

	return(USER_NO_PP);
	*/
    base_varnum=2*mpp->edges+mpp->arcs;
	*lb = (double *) malloc (base_varnum*sizeof(double));
    *ub = (double *) malloc (base_varnum*sizeof(double));
   
	/* i want to set the number of variables in the base*/ 
    
	*basevarnum=base_varnum;
	*basevars = (int *) realloc((char *)(*basevars), base_varnum * ISIZE);
	/* the cuts in the base go as follows: 
	//    (nodes) for each node indegree=outdegree
	//    (edges) for each edge the sum of weights in each direction >= 1
	*/
	*basecutnum=(mpp->nodes)+(mpp->edges);
/*  nonzeros
4 for each edge
2 for each arc
2 for each edge
  */
    for (i = 0; i <= *basevarnum-1; i++)
		{
		(*basevars)[i]=i;
		}
	for (i = 0; i <= mpp->arcs-1; i++)
		{
		(*lb)[i]=1;
		(*ub)[i]=mpp->arcs+mpp->edges;
		}

	for (i = mpp->arcs; i <= mpp->arcs+mpp->edges-1; i++)
		{
		(*lb)[i]=0;
		(*ub)[i]=mpp->arcs+mpp->edges;
		(*lb)[i+mpp->edges]=0;
		(*ub)[i+mpp->edges]=mpp->arcs+mpp->edges;
		}
    return(USER_NO_PP);
	}

			
            


/*===========================================================================*\
 * This is the second step in the process where the user specifies
 * which variables should be active in the root in addition to the base
 * set specified above
\*===========================================================================*/

/*===========================================================================*/

int user_create_root(void *user, int *extravarnum, int **extravars)
{
   mpp_problem *mpp=(mpp_problem *)user;
   *extravarnum = 0;
   *extravars  = NULL;
 
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Receive the feasible solution
\*===========================================================================*/

int user_receive_feasible_solution(void *user, int msgtag, double cost,
				   int numvars, int *indices, double *values)
{
   return(USER_NO_PP);
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
/*#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* This is is the case when we are copying data directly because
      the LP is not running separately. This code should be virtually
      identical to that of user_receive_lp_data() in the LP process.
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_lp_data() in the LP process 
#endif*/
	mpp_problem *mpp=(mpp_problem *)user;
	*user_lp = (void *)mpp;
   return(USER_NO_PP);
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
mpp_problem *mpp=(mpp_problem *)user;
	*user_cg = (void *)mpp;
   return(USER_NO_PP);
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
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CP)
   /* This is is the case when we are copying data directly because
      the LP is not running separately. This code should be virtually
      identical to that of user_receive_cp_data() in the CP process.*/
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cp_data() in the CP process */
#endif
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Generally, this function is not needed but you might find some use
 * for it. Someone did :).
\*===========================================================================*/

int user_process_own_messages(void *user, int msgtag)
{
   switch (msgtag){
    default:
      fprintf(stderr, "\nMaster: unknown message type %i!!!\n\n", msgtag);
      exit(1);
   }

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the user's chance to display the solution in whatever
 * manner desired. 
\*===========================================================================*/

int user_display_solution(void *user, double lpetol, int varnum, int *indices,
			  double *values, double objval)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* In this case, the LP user data structure is passed in */
#else
   /* In this case, it is the master user data structure */
#endif
	return(USER_NO_PP);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * This is a debugging feature which might
 * allow you to find out why a known feasible solution is being cut off.
\*===========================================================================*/

int user_send_feas_sol(void *user, int *feas_sol_size, int **feas_sol)
{
#ifdef TRACE_PATH

#endif
   return(USER_NO_PP);
}   

/*===========================================================================*/

/*===========================================================================*\
 * This function frees everything.
\*===========================================================================*/

int user_free_master(void **user)
{
   mpp_problem *mpp=(mpp_problem *)user;
   FREE(mpp->cost);   /* an array containing the cost for each edge*/
   FREE(mpp->start_node);    /* an array containing the start for each edge*/
   FREE(mpp->end_node);    /* an array containing the end for each edge*/
   FREE(mpp->types);  /* an array containing the types for each edge */
   
   return(USER_NO_PP);
}






