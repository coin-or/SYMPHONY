#ifndef _S_PATH_H
#define _S_PATH_H

#include "proto.h"
#include "heur_types.h"

typedef struct ADJ_LIST{
  int custnum;
  int cost;
  struct ADJ_LIST *next;
}adj_list;

        int *sp PROTO((adj_list **adj, int numnodes, int origin, int dest));
        void make_routes PROTO((heur_prob *p, _node *tsp_tour, int start, 
				 best_tours *new_tour));

#endif
