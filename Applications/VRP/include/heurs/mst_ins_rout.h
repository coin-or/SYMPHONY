#ifndef _MST_INS_ROUT_H
#define _MST_INS_ROUT_H

#include "proto.h"
#include "lb_types.h"
#include "heur_types.h"

int make_k_tree PROTO((lb_prob *p, int *tree, int *lamda, int k));
int closest PROTO((neighbor *nbtree, int *intree, int *last, int *host));
void ni_insert_edges PROTO((lb_prob *p, int new_node, neighbor *nbtree,
	int *intree, int *last, int *lamda, int mu));
int new_lamda PROTO((lb_prob *p, int upper_bound, int cur_bound, int *lamda,
	int numroutes, int *tree, edge_data *cur_edges, int alpha));

#endif
