#include <malloc.h>
#include <pvm3.h>

#include "sweep.h"
#include "messages.h"
#include "qsort.h"
#include "timemeas.h"
#include "compute_cost.h"
#include "pvm_error.h"
#include "vrp_const.h"
#include "heur_routines.h"

/*-----------------------------------------------------------------------*\
| Make_tour receives the nodes in sorted order and performs the algorithm |
| described below, starting with a number of equally spaced nodes as      |
| determined by the parameter p->par.sweep_trials.                        |
\*-----------------------------------------------------------------------*/

void make_tour(heur_prob *p, sweep_data *data, best_tours *final_tour)
{
  int i, k, weight = 0, j=1, l;
  int interval, start=0;
  int cost;
  _node *tour;
  int *demand = p->demand;
  int vertnum = p->vertnum;
  int capacity = p->capacity;

  if (p->par.sweep_trials>vertnum - 1) p->par.sweep_trials = vertnum-1;

  final_tour->cost = MAXINT;

  tour = (_node *) calloc (vertnum, sizeof(_node));

  interval = vertnum/p->par.sweep_trials;

  for (l=0; l < p->par.sweep_trials; l++){

    tour[0].next = data[start].cust;
    tour[0].route = 0;

    for (k=1; k < vertnum-1; k++){

      i = (k-1) + start;
      if (i>vertnum-2) i-=(vertnum-1);
      
      if (weight + demand[data[i].cust] <= capacity){
	weight += demand[data[i].cust];
	if (i != vertnum - 2)
	  tour[data[i].cust].next = data[i+1].cust;
	else
	  tour[data[i].cust].next = data[0].cust;
	tour[data[i].cust].route=j;
      }
      else{
	weight = demand[data[i].cust];
	j++;
	if (i != vertnum - 2)
	  tour[data[i].cust].next = data[i+1].cust;
	else
	  tour[data[i].cust].next = data[0].cust;
	tour[data[i].cust].route=j;
      }
    }
    i = (k-1) + start;
    if (i>vertnum-2) i-=(vertnum-1);
    if (weight + demand[data[i].cust] <= capacity){
      tour[data[i].cust].next = 0;
      tour[data[i].cust].route = j;
    }
    else{
      j++;
      tour[data[i].cust].next = 0;
      tour[data[i].cust].route = j;
    }
    
    cost = compute_tour_cost (&p->dist, tour);
    if (cost < final_tour->cost){
      memcpy(final_tour->tour, tour, vertnum*sizeof(_node));
      final_tour->cost = cost;
      final_tour->numroutes=j;
    }
    
    j=1;
    weight=0;
    start += interval;

  }
  free((char *)tour);
}

/*===========================================================================*/

/*-----------------------------------------------------------------*\
| The sweep algorithm is a very simple heuristic for clustering.    |
| We simply order the nodes radially about the depot and then add   |
| to the current route in this order until capacity is exceeded and |
| then we start a new route. Depending on where we start in the     |
| cyclic ordering, we get different solutions.                      |
\*-----------------------------------------------------------------*/

void main(void)
{
  heur_prob *p;
  int mytid, info, r_bufid, parent;
  int i;
  int vertnum;
  sweep_data *data;
  float depotx, depoty;
  float tempx, tempy;
  double t;

  (void) used_time(&t);

  mytid = pvm_mytid();

  p = (heur_prob *) calloc(1, sizeof(heur_prob));

  /*-----------------------------------------------------------------------*\
  |                     Receive the VRP data                                |
  \*-----------------------------------------------------------------------*/

  parent = receive(p);

  PVM_FUNC(r_bufid, pvm_recv(-1, SWEEP_TRIALS));
  PVM_FUNC(info, pvm_upkint(&(p->par.sweep_trials), 1, 1));
  PVM_FUNC(r_bufid, pvm_recv(-1, COORD_DATA));
  p->dist.coordx = (float *) calloc(p->vertnum, sizeof(float));
  p->dist.coordy = (float *) calloc(p->vertnum, sizeof(float));
  PVM_FUNC(info, pvm_upkfloat(p->dist.coordx, (int)p->vertnum, 1));
  PVM_FUNC(info, pvm_upkfloat(p->dist.coordy, (int)p->vertnum, 1));

  /*-----------------------------------------------------------------------*/

  vertnum = p->vertnum;
  p->cur_tour = (best_tours *)calloc(1, sizeof(best_tours));
  p->cur_tour->tour = (_node *)calloc(vertnum, sizeof(_node));

  data = (sweep_data *)calloc(vertnum-1, sizeof(sweep_data));
  if (p->dist.coordx && p->dist.coordy){
     depotx = p->dist.coordx[0];
     depoty = p->dist.coordy[0];
  }else{
     exit(1);
  }

  /*calculate angles for sorting*/
  for (i=0; i<vertnum-1; i++){
    tempx = p->dist.coordx[i+1] - depotx;
    tempy = p->dist.coordy[i+1] - depoty;
    data[i].angle = (float) atan2(tempy, tempx);
    if (data[i].angle < 0) data[i].angle += 2*M_PI;
    data[i].cust=i+1;
  }

  quicksort(data, vertnum-1);

  make_tour(p, data, p->cur_tour);

  /*-----------------------------------------------------------------------*\
  |               Transmit the tour back to the parent                      |
  \*-----------------------------------------------------------------------*/

  send_tour(p->cur_tour->tour, p->cur_tour->cost, p->cur_tour->numroutes,
	    SWEEP, used_time(&t), parent, vertnum, 0, NULL);

  if (data) free((char *) data);

  free_heur_prob(p);

  PVM_FUNC(r_bufid, pvm_recv(parent, YOU_CAN_DIE));
  PVM_FUNC(info, pvm_freebuf(r_bufid));
  PVM_FUNC(info, pvm_exit());
}

  

    

