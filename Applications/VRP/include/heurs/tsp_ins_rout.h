#ifndef _TSP_INS_ROUT_H
#define _TSP_INS_ROUT_H

#include "heur_types.h"
#include "proto.h"

int farthest_ins_from_to PROTO((
			 heur_prob *p, _node *tour, int cost,
			 int from_size, int to_size,
			 int starter, neighbor *nbtree,
			 int *intour, int *last));
int nearest_ins_from_to PROTO((
			heur_prob *p, _node *tour, int cost,
			int from_size, int to_size,
			int starter, neighbor *nbtree, 
			int *intour, int *last));
int closest PROTO((
		neighbor *nbtree, int *intour, int *last));
void ni_insert_edges PROTO((
		     heur_prob *p, int new_node, neighbor *nbtree,
		     int *intour, int *last));
int farthest PROTO((
		 neighbor *nbtree, int *intour, int *last));
void fi_insert_edges PROTO((
		     heur_prob *p, int new_node, neighbor *nbtree,
		     int *intour, int *last));
int insert_into_tour PROTO((
		     heur_prob *p, _node *tour, int starter,
		     int size, int new_node));

#endif
