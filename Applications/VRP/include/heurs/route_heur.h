#ifndef _ROUTE_HEUR_H
#define _ROUTE_HEUR_H

#include "proto.h"
#include "vrp_types.h"

void route_heur PROTO((vrp_problem *vrp, heur_params *heur_par,
		       heurs *rh, int trials));

#endif
