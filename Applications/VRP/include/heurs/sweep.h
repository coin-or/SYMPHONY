#ifndef _SWEEP_H
#define _SWEEP_H

#include "heur_types.h"
#include "proto.h"

typedef struct SWEEP_DATA{
   float angle;
   int cust;
}sweep_data;

void make_tour PROTO((heur_prob *p, sweep_data *data, best_tours *final_tour));

#endif
