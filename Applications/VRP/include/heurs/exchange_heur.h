#ifndef _EXCHANGE_HEUR_H
#define _EXCHANGE_HEUR_H

#include "proto.h"
#include "vrp_types.h"

void exchange_heur PROTO((vrp_problem *vrp, heurs *eh, int trials, int which));

#endif
