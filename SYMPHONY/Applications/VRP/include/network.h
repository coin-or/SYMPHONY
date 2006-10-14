/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _NETWORK
#define _NETWORK

/* SYMPHONY include files */
#include "proto.h"
#include "qsortucb.h"

/* VRP include files */
#include "vrp_const.h"

#define NOT_INTEGRAL -1
#define OTHER_END(cur_edge, v) \
        (cur_edge->data->v0 == v) ? cur_edge->data->v1 : cur_edge->data->v0

/*-----------------------------------------------------------------------*\
| These are data tructures used in constructing the solution graph used   |
| by the cut generator to locate cuts. The graph is stored using          |
| adjacency lists                                                         |
\*-----------------------------------------------------------------------*/

typedef struct EDGE{
   int v0;      
   int v1;
   int cost;
   double weight;  
   char scanned;
   char tree_edge;
   char deleted;
}edge;

typedef struct ELIST{
   struct ELIST  *next_edge; /* next edge in the edgelist */
   struct EDGE   *data;      /* the data of the edge */
   int            other_end; /* the other end of the edge */
   struct VERTEX *other;
}elist;

typedef struct VERTEX{
  int           enodenum;   /* the node number in the contracted graph */
  int           orignodenum;/* the node number in the original graph */
  struct ELIST *first; /* points to the first edge in the adjacency list */
  struct ELIST *last;  /* points to the last edge in the adjacency list */
  int           comp;  /* contains the component number if the graph is
			  disconnected */
  char          scanned;
  int           demand; /* contains the demand for this node */
  int           degree; /* contains the degree of the node in the graph */
  int           orig_node_list_size;
  int          *orig_node_list; /* contains a list of the nodes that have been
				   contracted into this node to make a
				   "super node" */
  int           dfnumber;
  int           low;
  char          is_art_point;
  char          deleted; 
}vertex;

typedef struct NETWORK{
  int             vertnum; /* the number of vertices in the graph */
  char            is_integral; /* indicates whether the graph is integral or
				  not */
  int             edgenum; /* the number of edges in the graph */
  float           mincut;  /* the value of the current mincut */
  struct ELIST   *adjlist; /* the array containing the adajacency lists for
			      each node */
  struct EDGE    *edges; /* the list of edges in the graph */
  struct VERTEX  *verts; /* the list of vertices */
  int            *compnodes;  /* number of nodes in each component */
  double         *new_demand; /* the amounts of demand for each node that gets
				 added to it when the network is contracted */
}network;

network *createnet PROTO((int *xind, double *xval, int edgenum, double etol,
			  int *edges, int *demands, int vertnum));
int connected PROTO((network *n, int *compnodes, int *compdemands,
		   int *compmembers, double *compcuts, double *compdensity));
void free_net PROTO((network *n));

#endif
