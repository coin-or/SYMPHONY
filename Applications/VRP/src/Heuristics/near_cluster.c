#include <pvm3.h>
#include <malloc.h>

#include "ins_routines2.h"
#include "timemeas.h"
#include "messages.h"
#include "vrp_const.h"
#include "heur_routines.h"
#include "pvm_error.h"
#include "compute_cost.h"

void main(void)
{
  heur_prob *p;
  int vertnum;
  int mytid, info, r_bufid, parent;
  int cost;
  int *intour, last = 0;
  int numroutes, i;
  neighbor *nbtree;
  _node *tour;
  route_data *route_info;
  int cur_route;
  best_tours *tours;
  double t;

  (void) used_time(&t);

  mytid = pvm_mytid();
	
  p = (heur_prob *) calloc ((int)1, sizeof(heur_prob));
  tours = p->cur_tour = (best_tours *) calloc (1, sizeof(best_tours));
	
  /*-----------------------------------------------------------------------*\
  |                    Receive the VRP data                                 |
  \*-----------------------------------------------------------------------*/

  parent = receive(p);

  PVM_FUNC(r_bufid, pvm_recv(-1, NUMROUTES));
  PVM_FUNC(info, pvm_upkint(&p->numroutes, 1, 1));
  if (!p->numroutes){
     PVM_FUNC(info, pvm_recv(-1, COORD_DATA));
     p->dist.coordx = (float *) calloc (p->vertnum, sizeof(float));
     p->dist.coordy = (float *) calloc (p->vertnum, sizeof(float));
     PVM_FUNC(info, pvm_upkfloat(p->dist.coordx, p->vertnum, 1));
     PVM_FUNC(info, pvm_upkfloat(p->dist.coordy, p->vertnum, 1));
  }

  vertnum = p->vertnum;
  numroutes = p->numroutes;

  /*-----------------------------------------------------------------------*\
  |                    Allocate arrays                                      |
  \*-----------------------------------------------------------------------*/

  nbtree   = (neighbor *) malloc (vertnum * sizeof(neighbor));
  intour   = (int *)      calloc (vertnum, sizeof(int));
  tour = p->cur_tour->tour
       = (_node *) calloc (vertnum, sizeof(_node));

  /*-----------------------------------------------------------------------*\
  | First we generate seed customers for all the routes                     |
  \*-----------------------------------------------------------------------*/

  seeds(p, &numroutes, intour, nbtree);

  route_info = p->cur_tour->route_info;

  p->cur_tour->numroutes = numroutes;

  /*-----------------------------------------------------------------------*\
  | This algorithm builds all routes simultaneously. For each customer not  |
  | already in the solution, it keeps track of the route that it is closest |
  | to (i.e. the route # of the customer already in the solution that it is |
  | closest to). At each step, the customer in closest proximity to its host|
  | route is added to that route and all distances are updated. If it is    |
  | found to be infeasible to add the customer to the route, then a new host|
  | is found and we start the process again. Note that as we go through the |
  | algorithm, the last node of each route always points towards the depot  |
  \*-----------------------------------------------------------------------*/

  if (nbtree) free ((char *) nbtree);

  nbtree = (neighbor *) calloc (vertnum, sizeof(neighbor));

  for( i = 0; i<vertnum; i++){
    if (intour[i] != IN_TOUR)
      intour[i] = 0;
    tour[i].next = 0;
  }
  
  for (cur_route = 1; cur_route<=numroutes; cur_route++)
    ni_insert_edges(p, route_info[cur_route].first, nbtree, intour, &last,
		    tour, route_info);
  
  /* Form the routes by nearest insertion as described above */
    
  nearest_ins(p, tour, route_info, numroutes+1, vertnum,
		      nbtree, intour, &last);

  tour[0].next = route_info[1].first;

  
  /*-------------------------------------------------------------------------*\
  | This loop points the last node of each route to the first node of the next|
  | route. At the end of this procedure, the last node of each route is       |
  | pointing at the depot, which is not what we want.                         |
  \*-------------------------------------------------------------------------*/

  for (cur_route = 1; cur_route < numroutes; cur_route++)
    tour[route_info[cur_route].last].next = route_info[cur_route+1].first;

  cost = compute_tour_cost(&p->dist, tour);
	
  /*-----------------------------------------------------------------------*\
  |               Transmit the tour back to the parent                      |
  \*-----------------------------------------------------------------------*/
  
  send_tour(tour, cost, numroutes, NEAR_CLUSTER, used_time(&t), parent,
	    vertnum, 1, 
	    route_info);
	
  if ( nbtree ) free ((char *) nbtree);
  if ( intour ) free ((char *) intour);

  free_heur_prob(p);
	
  PVM_FUNC(r_bufid, pvm_recv(parent, YOU_CAN_DIE));
  PVM_FUNC(info, pvm_freebuf(r_bufid));
  PVM_FUNC(info, pvm_exit());
	
}
