#include <pvm3.h>
#include <malloc.h>
#include <stdlib.h>

#include "savings.h"
#include "messages.h"
#include "timemeas.h"
#include "pvm_error.h"
#include "heur_routines.h"
#include "vrp_const.h"
#include "compute_cost.h"

void main(void)
{
  heur_prob *p;
  _node *tour;
  int mytid, info, r_bufid, parent;
  int i, capacity;
  int vertnum;
  int cur_route=1;
  int weight, *demand;
  int savings, start, *intour;
  int node1, node2, ins_cust, starter;
  int cur_route_end = 0, prev_route_end = 0;
  double t;

  (void) used_time(&t);

  mytid = pvm_mytid();

  p = (heur_prob *) calloc(1, sizeof(heur_prob));

  /*-----------------------------------------------------------------------*\
  |                     Receive the VRP data                                |
  \*-----------------------------------------------------------------------*/

  parent = receive(p);

  /*-----------------------------------------------------------------------*\
  |                     Receive the parameters                              |
  \*-----------------------------------------------------------------------*/
  PVM_FUNC(r_bufid, pvm_recv(-1, SAVINGS_DATA));
  PVM_FUNC(info, pvm_upkfloat(&p->par.savings_par.mu, 1, 1));
  PVM_FUNC(info, pvm_upkfloat(&p->par.savings_par.lamda, 1, 1));
  PVM_FUNC(info, pvm_upkint(&start, 1, 1));
  capacity = p->capacity;
  demand = p->demand;
  vertnum = p->vertnum;

  if (start != FAR_INS) srand(start);

  p->cur_tour = (best_tours *) calloc (1, sizeof(best_tours));
  tour = p->cur_tour->tour = (_node *) calloc (vertnum, sizeof(_node));
  intour = (int *) calloc (vertnum, sizeof(int));

  /*------------------------------------------------------------------*\  
  |  Initiate the first route by placing the initial starter (either   |
  |  the farthest node from the depot of a random node) and the depot  |
  |  on the first route                                                |
  \*------------------------------------------------------------------*/
  starter = new_start(intour, p, start, vertnum - 1); 
  tour[0].next = starter;
  tour[starter].route = cur_route;
  intour[starter] = IN_TOUR;
  intour[0] = IN_TOUR;
  cur_route_end = starter;
  weight = demand[starter];

  /*------------------------------------------------------------------*\  
  |  The following loop first finds the customer with the maximum      |
  |  savings among all customers not yet put on routes and returns the |
  |  the position in which to insert the customer in order to achieve  |
  |  that savings. It then checks the feasibility of inserting that    |
  |  that customer on the current route and if it is infeasible, starts|
  |  a new route.                                                      |
  \*------------------------------------------------------------------*/
  for (i=0; i<vertnum - 2; i++){
    find_max(&ins_cust, &savings, &node1, &node2,
	     tour, intour, prev_route_end, p);
    if (weight + demand[ins_cust] <= capacity){
      intour[ins_cust] = IN_TOUR;
      weight += demand[ins_cust];
      insert_cust(ins_cust, tour, node1, node2, 
		  cur_route, prev_route_end);
      if (node2 == 0)
	cur_route_end = ins_cust;
    }
    else{
      cur_route++;
      tour[cur_route_end].next = new_start(intour, p, start, vertnum - 2 - i);
      prev_route_end = cur_route_end;
      cur_route_end = tour[prev_route_end].next;
      tour[cur_route_end].route = cur_route;
      intour[cur_route_end] = IN_TOUR;
      weight = demand[cur_route_end];
    }
  }

  p->cur_tour->cost = compute_tour_cost(&p->dist, tour);
  p->cur_tour->numroutes = cur_route;

  /*-----------------------------------------------------------------------*\
  |               Transmit the tour back to the parent                      |
  \*-----------------------------------------------------------------------*/

  send_tour(tour, p->cur_tour->cost, p->cur_tour->numroutes, SAVINGS,
	    used_time(&t),
	    parent, vertnum, 0, NULL);

  if (intour) free ((char *) intour);

  free_heur_prob(p);

  PVM_FUNC(r_bufid, pvm_recv(parent, YOU_CAN_DIE));
  PVM_FUNC(info, pvm_freebuf(r_bufid));
  PVM_FUNC(info, pvm_exit());
}
