#ifndef _SAVINGS2_H
#define _SAVINGS2_H

#include "proto.h"
#include "binomial.h"
#include "vrp_common_types.h"
#include "heur_types.h"

void print_routes PROTO((_node *tour));
void insert_cust PROTO((int cust_num, _node *tour, int node1,
			int node2, int cur_route,
			int prev_route_end));
tree_node *start_new_route PROTO((heur_prob *p, tree_node *head,
				  int starter));
tree_node *update_savings PROTO(( heur_prob *p, tree_node *head,
				 tree_node *mav_ptr, _node *tour, 
				 int prev_route_end));
int new_start PROTO((int *intour, heur_prob *p,
			 int start, int num_cust));
void update PROTO((tree_node *cur_node, int savings, int node1,
		   int node2));
int new_savings PROTO((heur_prob *p, tree_node *max_ptr, tree_node *head,
		       _node *tour, int prev_route_end, int *node1,
		       int *node2));

#endif
