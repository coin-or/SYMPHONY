#ifndef _INS_ROUTINES2_H
#define _INS_ROUTINES2_H

#include "proto.h"
#include "heur_types.h"
#include "vrp_common_types.h"

void nearest_ins PROTO((
          heur_prob *p, _node *tour, route_data *route_info, int from_size, 
          int to_size, neighbor *nbtree, int *intour, int *last));
int closest PROTO(( neighbor *nbtree, int *intour, int *last, int *host));
void ni_insert_edges PROTO((
	  heur_prob *p, int new_node, neighbor *nbtree, int *intour, int *last,
	  _node *tour, route_data *route_info));
int insert_into_tour PROTO((
          heur_prob *p, _node *tour, int new_node, route_data *route_info));
void new_host PROTO((
	  heur_prob *p, int node, neighbor *nbtree, int *intour, int *last,
	  _node *tour, route_data *route_info));
void seeds PROTO((
	  heur_prob *p, int *numroutes, int *intour, neighbor *nbtree));
void farthest_ins_from_to PROTO((
	  heur_prob *p, _node *tour, int from_size, 
	  int to_size, neighbor *nbtree, 
	  int *intour, int *last));
int farthest PROTO((
	  neighbor *nbtree, int *intour, int *last));
void fi_insert_edges PROTO((
	  heur_prob *p, int new_node, neighbor *nbtree, int *intour, int *last));

#endif
