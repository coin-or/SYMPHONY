#ifndef _SAVINGS_H
#define _SAVINGS_H

#include "proto.h"
#include "heur_types.h"
#include "vrp_macros.h"

#define SAV(d, a, b, c) (p->par.savings_par.lamda) * ICOST(d, 0, c) - \
                       (ICOST(d,a,c) + ICOST(d,b,c) -  \
			(p->par.savings_par.mu) * ICOST(d,a,b))

void find_max PROTO((
	      int *ins_cust, int *savings, int *node1,
	      int *node2, _node *tour, int *intour,
	      int prev_route_end, heur_prob *p));
void print_routes PROTO((_node *tour));
void insert_cust PROTO((
	      int cust_num, _node *tour, int node1,
	      int node2, int cur_route, int prev_route_end));
int new_start PROTO((int *intour, heur_prob *p, 
	      int start, int num_cust));

#endif
