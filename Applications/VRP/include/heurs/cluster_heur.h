#ifndef _CLUSTER_HEUR_H
#define _CLUSTER_HEUR_H

#include "proto.h"
#include "vrp_types.h"

void cluster_heur PROTO((vrp_problem *vrp, heur_params *heur_par,
			 heurs *ch, int trials));
void generate_starter PROTO((int vertnum, int *starter, 
			     int *startpos, int num));

#endif
