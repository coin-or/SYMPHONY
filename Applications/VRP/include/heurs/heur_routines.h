#ifndef _HEUR_ROUTINES_H
#define _HEUR_ROUTINES_H

#include "proto.h"
#include "heur_types.h"
#include "lb_types.h"

int receive PROTO((heur_prob *p));
void send_tour PROTO((_node *tour, int cost, int numroutes, int algorithm, 
		      double cpu_time, int parent, int vertnum,
		      int routes, route_data *route_info));
void free_heur_prob PROTO((heur_prob *p));
void free_lb_prob PROTO((lb_prob *p));

#endif
