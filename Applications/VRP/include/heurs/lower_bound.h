#ifndef _LOWER_BOUND_H
#define _LOWER_BOUND_H

#include "proto.h"
#include "vrp_types.h"
#include "heur_types.h"
#include "lb_params.h"

void lower_bound PROTO((vrp_problem *vrp, lb_params *lb_par,
			heurs *lh, int ub));

#endif
