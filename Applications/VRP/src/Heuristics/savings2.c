#include <malloc.h>
#include <stdlib.h>

#include "savings2.h"
#include "timemeas.h"
#include "messages.h"
#include "BB_constants.h"
#include "heur_routines.h"
#include "compute_cost.h"
#include "vrp_const.h"
#include "proccomm.h"

void main(void)
{
  heur_prob *p;
  _node *tour;
  int mytid, info, r_bufid, parent;
  int i, capacity;
  int vertnum;
  int cur_route=1;
  int weight = 0, *demand;
  int starter;
  tree_node *head, *max_ptr;
  int savings, start, *intour;
  int node1, node2, cust_num, num_cust = 0;
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
  vertnum = p->vertnum;
  demand = p->demand;

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
  tour[starter].next = 0;
  tour[0].route = 0;
  tour[starter].route = cur_route;
  intour[0] = IN_TOUR;
  intour[starter] = IN_TOUR;
  cur_route_end = starter;
  weight += demand[starter];

  /*------------------------------------------------------------------*\
  | Initialize the heap with node 1 or node 2, whichever is not the    |
  | starter.                                                           |
  \*------------------------------------------------------------------*/

  if (starter != 1){
    savings = (int) (SAV(&p->dist, 0, starter, 1));
    node1 = 0;
    node2 = starter;
    cust_num = 1;
  }
  else{
    savings = (int) (SAV(&p->dist, 0, 1, 2));
    node1 = 0;
    node2 = 1;
    cust_num = 2;
  }

  head = make_heap(cust_num, savings, node1, node2);

  /*----------------------------------------------------------*\
  |          Make the initial heap                             |
  \*----------------------------------------------------------*/

  if (cust_num == 1){
    for (i=2; i<vertnum; i++){
      if (i != starter){
	savings = (int) (SAV(&p->dist, 0, starter, i));
	cust_num = i;
	head = heap_insert(head, cust_num, savings, node1, node2);
      }
    }
  }
  else{
    for (i=3; i<vertnum; i++){
      savings = (int) (SAV(&p->dist, 0, 1, i));
      cust_num = i;
      head = heap_insert(head, cust_num, savings, node1, node2);
    }
  }

  /*------------------------------------------------------------------*\  
  |  The following loop first finds the customer with the maximum      |
  |  savings among all customers not yet put on routes and returns the |
  |  the position in which to insert the customer in order to achieve  |
  |  that savings. It then checks the feasibility of inserting that    |
  |  that customer on the current route and if it is infeasible, starts|
  |  a new route.                                                      |
  \*------------------------------------------------------------------*/
    
  while (head != NULL){
    max_ptr = find_max(head);
    if (weight + demand[max_ptr->cust_num] <= capacity){
      head = extract_max(head, max_ptr);
      intour[max_ptr->cust_num] = IN_TOUR;
      weight += demand[max_ptr->cust_num];
      insert_cust(max_ptr->cust_num, tour, max_ptr->node1,
		  max_ptr->node2, cur_route, prev_route_end);
      if (max_ptr->node2 == 0)
	cur_route_end = max_ptr->cust_num;
      head = update_savings(p, head, max_ptr, tour, prev_route_end);
      free(max_ptr);
      num_cust++;
    }
    else{
      cur_route++;
      tour[cur_route_end].next = new_start(intour, p, 
					   start, vertnum-2-num_cust);
      prev_route_end = cur_route_end;
      cur_route_end = tour[prev_route_end].next;
      tour[cur_route_end].route = cur_route;
      intour[cur_route_end] = IN_TOUR;
      head = start_new_route(p, head, cur_route_end);
      weight = demand[cur_route_end];
      num_cust++;
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
      
