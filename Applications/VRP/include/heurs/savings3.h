#ifndef _SAVINGS_3_H
#define _SAVINGS_3_H

#include "binomial.h"
#include "proto.h"
#include "vrp_common_types.h"
#include "heur_types.h"

void insert_cust PROTO((int cust_num, _node *tour, int node1,
			int node2, int cur_route, int *demand,
			route_data *route_info));
int new_savings PROTO((heur_prob *p, tree_node *max_ptr, tree_node *head,
		       _node *tour, int *node1, int *node2,
		       route_data *route_info));
tree_node *update_savings PROTO((heur_prob *p, tree_node *head,
				 tree_node *max_ptr, _node *tour,
				 route_data *route_info));
int find_new_ins_route PROTO((heur_prob *p, int ins_node, _node *tour,
			      int *node1, int *node2,
			      route_data *route_info));

#endif



