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

#include <stdio.h>
#include <malloc.h>
#include "BB_constants.h"
#include "BB_macros.h"
#include "lp_u.h"
#include "mpp.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions for the LP process.
\*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_lp_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_LP or COMPILE_IN_TM is FALSE.
\*===========================================================================*/

int user_receive_lp_data(void **user)
{
   mpp_problem *mpp = (mpp_problem *) calloc(1, sizeof(mpp_problem));

   *user = mpp;
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Free all the user data structures
\*===========================================================================*/

int user_free_lp(void **user)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must create the initial LP relaxation for
 * each search node. See the comments below.
\*===========================================================================*/

int user_create_lp(void *user, int varnum, var_desc **vars, int rownum,
		   int cutnum, cut_data **cuts, int *nz, int **matbeg,
		   int **matind, double **matval, double **obj, double **rhs,
		   char **sense, double **rngval, int *maxn, int *maxm,
		   int *maxnz, int *allocn, int *allocm, int *allocnz)
{
   mpp_problem *mpp = (mpp_problem *)user;
   int i;
   char resize = FALSE;
   int x;
	int y;
	int first;
	double switcher;
	int adjuster;
	int j;
    double holder;
	int total_edgenum;
	total_edgenum=2*mpp->edges+mpp->arcs;
   /* set up the inital LP data */

   /* Our base constraints are that the sum of the weights of the
      edges adjacent to each node (except the depot) is 2. This means
      that each column will have exactly two nonzeros, one in each of
      the rows corresponding to its end points. Hence the total
      nonzeros is 2*n (n is the number of active variables). */
   *nz = (6 * mpp->edges)+ (2 * mpp->arcs);

   /* We have to check to make sure there is enough space allocated
      for the matrix we are going to build */
   if (2 * rownum > *maxm){
      *maxm = 2 * rownum;
      resize = TRUE;
   }
     /* Allocate space for all edges up front since we have small problems */
   if (total_edgenum != *maxn){
      *maxn = total_edgenum;
      resize = TRUE;
  
	 } 
   
   if (*nz + ((*maxm) * (*maxn) / 10) > *maxnz){
      *maxnz = *nz + ((*maxm) * (*maxn) / 10);
      resize = TRUE;
   }

   /* If there was not enough space, the allocate more */
   if (resize){
      /*re-malloc all the arrays*/
      FREE(*matbeg);
      FREE(*matind);
      FREE(*matval);
      FREE(*obj);
      FREE(*rhs);
      FREE(*sense);
      FREE(*rngval);
      *allocm  = *maxm;
      *allocn  = (*maxm + *maxn + 1);
      *allocnz = *maxnz + *maxm;
      *matbeg  = (int *) malloc(*allocn * ISIZE);
      *matind  = (int *) malloc(*allocnz * ISIZE);
      *matval  = (double *) malloc(*allocnz * DSIZE);
      *obj     = (double *) malloc(*allocn * DSIZE);
      *rhs     = (double *) malloc(*allocm * DSIZE);
      *sense   = (char *) malloc(*allocm * CSIZE);
      *rngval  = (double *) calloc(*allocm, DSIZE);
   }

   /* Fill out the appropriate data structures -- each column has
      exactly two entried*/
   

	x=0;
	y=0;
	first=0;
	switcher=1; /* used in a shortcut to make variables for both directions of edges*/
	adjuster=0;
    /*checks each edge and arc to make the in degree equals out degree constraint*/
	for(i=0;i<=2*(mpp->edges)+(mpp->arcs)-1;i++)
		{
		first=0;     
		if((i>=mpp->edges+mpp->arcs))
			{
			switcher=1;
		    adjuster=mpp->edges;
			}
		else
			{
			switcher=-1;
			adjuster=0;
			}
		
	    for(j=0;j<=mpp->nodes-1;j++)
			{
		/*/checks to see if node i is the start node of every edge arc*/
			if(mpp->start_node[i-adjuster]==j) 
				{
				(*matind)[x]=j;
				(*matval)[x]=switcher;
				holder=(*matval)[x];
				x++;
				if(first==0)
					{
					first=1;
					(*matbeg)[i]=x-1;
					}
					
				}
			else if(mpp->end_node[i-adjuster]==j)
				{
			    (*matind)[x]=j;
				(*matval)[x] = switcher*-1;
			
				holder=(*matval)[x];
				x++;
				if(first==0)
					{
					first=1;
					(*matbeg)[i]=x-1;
					}
				}
			}

	   
		/* now i need to make the constraint for the sum of the edge weights in each direction
		// checks to see if it is an edge */
		if(i>=(mpp->arcs))
			{
			(*matind)[x]=mpp->nodes+i-mpp->arcs-adjuster;
			(*matval)[x] = 1;
			holder=(*matval)[x];
			x++;
			}
		}
 (*matbeg)[2*(mpp->edges)+(mpp->arcs)]=*nz;
   for (i = 0; i <= mpp->arcs-1; i++)
		{
		(*obj)[i] = mpp->cost[i];
		}

	for (i = mpp->arcs; i <= mpp->arcs+mpp->edges-1; i++)
		{
		(*obj)[i] = mpp->cost[i];
		(*obj)[i+mpp->edges] = mpp->cost[i];
		}
   
   /*Someday, I should add in all the cuts in here too, but for now, they get
     added in during post-processing in create_lp_u(), which is less efficient
     but easier for me :)*/
   
   /* set the initial right hand side */
  
	for (i = 0; i <= mpp->nodes-1 ; i++)
		{
		(*rhs)[i]   = 0;
		(*sense)[i] = 'E';
		}
	for (i = mpp->nodes; i <= mpp->nodes+mpp->edges-1 ; i++)
		{
		(*rhs)[i]   = 1;
		(*sense)[i] = 'G';
		}

return(USER_NO_PP);
   }

  


/*===========================================================================*/

/*===========================================================================*\
 * This function takes an LP solution and checks it for feasibility.
\*===========================================================================*/

int user_is_feasible(void *user, double lpetol, int varnum, int *indices,
		     double *values, int *feasible, double *objval)
{
   /* i think here should just be something that checks all the variables and
      makes sure they are integral */
	
	return(DEFAULT);
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

int user_send_feasible_solution(void *user, double lpetol, int varnum,
				int *indices, double *values)
{
   return(DEFAULT);
}


/*===========================================================================*/

/*===========================================================================*\
 * This function graphically displays the current fractional solution
 * This is done using the Interactie Graph Drawing program.
\*===========================================================================*/

int user_display_lp_solution(void *user, int which_sol, int varnum,
			     int *indices, double *values)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You can add whatever information you want about a node to help you
 * recreate it. I don't have a use for it, but maybe you will.
\*===========================================================================*/

int user_add_to_desc(void *user, int *desc_size, char **desc)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Compare cuts to see if they are the same. We use the default, which
 * is just comparing byte by byte.
\*===========================================================================*/

int user_same_cuts(void *user, cut_data *cut1, cut_data *cut2, int *same_cuts)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function receives a cut, unpacks it, and adds it to the set of
 * rows to be added to the LP.
\*===========================================================================*/

int user_unpack_cuts(void *user, int from, int type, int varnum,
		     var_desc **vars, int cutnum, cut_data **cuts,
		     int *new_row_num, waiting_row ***new_rows)
{
   mpp_problem *mpp = (mpp_problem *)user;
   int * node_checker;
   int * cut_holder=NULL;
   int edge_direction=0;
   int i, j, nzcnt = 0;
   waiting_row **row_list = NULL;
   int *matind = NULL;
   cut_data *cut;
   int rhs_count;
   char *coef;
   double *matval = NULL;
   *new_row_num = cutnum;
	node_checker= (int *) calloc(mpp->nodes, sizeof(int));
  

   if (cutnum > 0)
      *new_rows = row_list = (waiting_row **) calloc (cutnum,
						     sizeof(waiting_row *));

  for (j = 0; j < cutnum; j++)
	 {
     coef = (cut = cuts[j])->coef;
     cut_holder=(int*)cut->coef;

	 cuts[j] = NULL;
     (row_list[j] = (waiting_row *) malloc(sizeof(waiting_row)))->cut = cut;
     switch (cut->type)
		{
		case 1:
	    matind = (int *) malloc(varnum * ISIZE);
		nzcnt=0;
		rhs_count=0;
		for(i=0;i<=mpp->nodes-1;i++)
			{
			node_checker[i]=0;
			}
        /*make array for 1 if node is in cut, 0 if not*/
		for(i=0;i<(cut->size/4);i++)
			{
			node_checker[cut_holder[i]]=1;
			}
        for(i=0;i<=varnum-1;i++)
			{
			if(vars[i]->userind >= mpp->arcs+mpp->edges)
				{
				if(node_checker[mpp->end_node[(vars[i]->userind)-mpp->edges]]==1 && 
				   node_checker[mpp->start_node[(vars[i]->userind)-mpp->edges]]==0)
					{
					matind[nzcnt]=i;
					nzcnt++;
					rhs_count++;
					}
				}
			else if((vars[i]->userind >= mpp->arcs))
				{
				if((node_checker[mpp->start_node[vars[i]->userind]]==1 && node_checker[mpp->end_node[vars[i]->userind]]==0))
					{
					matind[nzcnt]=i;
					nzcnt++;
					rhs_count++;
					}
				}
			else if
				(node_checker[mpp->start_node[vars[i]->userind]]==1 && node_checker[mpp->end_node[vars[i]->userind]]==0)
				{
				matind[nzcnt]=i;
				nzcnt++;
				rhs_count++;
				}
			else if(node_checker[mpp->start_node[vars[i]->userind]]==0 && node_checker[mpp->end_node[vars[i]->userind]]==1)
				{
				rhs_count++;
				}
			
			}
       
		row_list[j]->matind = matind =
	    (int *) realloc((char *)matind, nzcnt*ISIZE);
        cut->rhs=(rhs_count+1)/2;
		row_list[j]->nzcnt = nzcnt;
		row_list[j]->matval = matval = (double *) malloc(nzcnt * DSIZE);
		for (i = nzcnt-1; i >= 0; i--)
			matval[i] = 1;
		cut->branch = ALLOWED_TO_BRANCH_ON;
		break;	

   	 
         default:
		 printf("Unrecognized cut type!\n");
		}
   }
  FREE(node_checker);
   return(USER_NO_PP);
}

/*===========================================================================*/

int user_send_lp_solution(void *user, int varnum, var_desc **vars, double *x,
			  int where)
{
   return(SEND_NONZEROS);
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine does logical fixing of variables
\*===========================================================================*/

int user_logical_fixing(void *user, int varnum, var_desc **vars, double *x,
			char *status, int *num_fixed)
{
   *num_fixed = 0;

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function generates the 'next' column
\*===========================================================================*/

int user_generate_column(void *user, int generate_what, int cutnum,
			 cut_data **cuts, int prevind, int nextind,
			 int *real_nextind, double *colval, int *colind,
			 int *collen, double *obj)
{
   switch (generate_what){
    case GENERATE_NEXTIND:
      /* Here we just have to generate the specified column. */
      break;
    case GENERATE_REAL_NEXTIND:
      /* In this case, we have to determine what the "real" next edge is*/
      break;
   }

   return(USER_NO_PP);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to print some statistics on the types and quantities
 * of cuts or something like that.
\*===========================================================================*/

int user_print_stat_on_cuts_added(void *user, int rownum, waiting_row **rows)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to eliminate rows from the local pool based on
 * knowledge of problem structure.
\*===========================================================================*/

int user_purge_waiting_rows(void *user, int rownum, waiting_row **rows,
			    char *delete)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user has to generate the ubber bounds for the specified
 * variables. Lower bounds are always assumed (w.l.o.g.) to be zero.
\*===========================================================================*/

int user_get_upper_bounds(void *user, int varnum, int *indices, double *bd)
{
   return(DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user might want to generate cuts in the LP using information
 * about the current tableau, etc. This is for advanced users only.
\*===========================================================================*/

int user_generate_cuts_in_lp(void *user, int varnum, var_desc **vars, double *x,
			     int *new_row_num, waiting_row ***new_rows)
{
   return(USER_NO_PP);
}

/*===========================================================================*/

