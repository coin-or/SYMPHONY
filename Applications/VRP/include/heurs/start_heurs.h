#ifndef _START_HEURS_H
#define _START_HEURS_H

#include "proto.h"
#include "vrp_types.h"
#include "heur_types.h"
#include "lb_types.h"

void start_heurs PROTO((vrp_problem *vrp, heur_params *heur_par,
			lb_params *lb_par, double *ub, char no_windows));

#endif
