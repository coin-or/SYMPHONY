#ifndef _RECEIVE_ROUT_H
#define _RECEIVE_ROUT_H

#include "proto.h"
#include "vrp_types.h"

double receive_tours PROTO((vrp_problem *vrp, heurs *hh, int *last, char print,
			    char routes, char add_edges, char win));
double receive_lbs PROTO((vrp_problem *vrp, heurs *hh,
			  char win, int numroutes));

#endif
