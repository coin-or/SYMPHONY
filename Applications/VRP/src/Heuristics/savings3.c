#include <malloc.h>
#include <pvm3.h>

#include "savings3.h"
#include "timemeas.h"
#include "messages.h"
#include "pvm3.h"
#include "BB_constants.h"
#include "heur_routines.h"
#include "compute_cost.h"
#include "ins_routines2.h"
#include "pvm_error.h"
#include "vrp_const.h"

void main(void)
{
   heur_prob *p;
   _node *tour;
   int mytid, info, r_bufid, parent;
   int i, capacity;
   int vertnum;
   int cur_route=1;
   int *demand;
   tree_node *head, *max_ptr;
   int savings, *intour, max_savings;
   int numroutes, v0, v1, ins_route, cust_num;
   double t;
   neighbor *nbtree;
   route_data *route_info;
   
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
   PVM_FUNC(r_bufid, pvm_recv(-1, NUMROUTES));
   PVM_FUNC(info, pvm_upkint(&p->numroutes, 1, 1));
   if (!p->numroutes){
      PVM_FUNC(info, pvm_recv(-1, COORD_DATA));
      p->dist.coordx = (float *) calloc (p->vertnum, sizeof(float));
      p->dist.coordy = (float *) calloc (p->vertnum, sizeof(float));
      PVM_FUNC(info, pvm_upkfloat(p->dist.coordx, p->vertnum, 1));
      PVM_FUNC(info, pvm_upkfloat(p->dist.coordy, p->vertnum, 1));
   }
   PVM_FUNC(r_bufid, pvm_recv(-1, SAVINGS_DATA));
   PVM_FUNC(info, pvm_upkfloat(&p->par.savings_par.mu, 1, 1));
   PVM_FUNC(info, pvm_upkfloat(&p->par.savings_par.lamda, 1, 1));
   capacity = p->capacity;
   vertnum = p->vertnum;
   demand = p->demand;
   numroutes = p->numroutes;
   
   p->cur_tour = (best_tours *) calloc (1, sizeof(best_tours));
   tour = p->cur_tour->tour = (_node *) calloc (vertnum, sizeof(_node));
   intour   = (int *)      calloc (vertnum, sizeof(int));
   nbtree   = (neighbor *) malloc (vertnum * sizeof(neighbor));
   
   seeds(p, &numroutes, intour, nbtree);

   p->numroutes = numroutes;
   
   route_info = p->cur_tour->route_info;
   p->cur_tour->numroutes = numroutes;
   
   for( i = 0; i<vertnum; i++){
      if (intour[i] != IN_TOUR)
	 intour[i] = 0;
      tour[i].next = 0;
   }
   
   for (i = 0; intour[i] == IN_TOUR; i++);
   for (max_savings = MAXINT, cur_route = 1; cur_route <= numroutes;
	cur_route++){
      if ((savings = SAV(&p->dist, 0, route_info[cur_route].first, i))
	  < max_savings){
	 max_savings = savings;
	 v1 = route_info[cur_route].first;
      }
   }
   head = make_heap(i, max_savings, 0, v1);
   
   /*----------------------------------------------------------*\
   |          Make the initial heap                             |
   \*----------------------------------------------------------*/
   
   for (i++; i<vertnum; i++){
      if (intour[i] == IN_TOUR)
	 continue;
      for (max_savings = MAXINT, cur_route = 1; cur_route <= numroutes;
	   cur_route++){
	 if (((savings = SAV(&p->dist, 0, route_info[cur_route].first, i))
	      < max_savings) && (route_info[cur_route].weight +
				 demand[i] <= capacity)){
	    max_savings = savings;
	    v1 = route_info[cur_route].first;
	 }
      }
      head = heap_insert(head, i, max_savings, 0, v1);
   }
      
   /*------------------------------------------------------------------*\  
   |  The following loop first finds the customer with the maximum      |
   |  savings among all customers not yet put on routes and returns the |
   |  the position in which to insert the customer in order to achieve  |
   |  that savings. It then checks the feasibility of inserting that    |
   |  that customer on the current route and if it is infeasible, starts|
   |  a new route.                                                      |
   \*------------------------------------------------------------------*/

   cust_num = numroutes+1;
   while (head != NULL){
      max_ptr = find_max(head);
      head = extract_max(head, max_ptr);
      ins_route = max_ptr->node1 ? tour[max_ptr->node1].route :
	 tour[max_ptr->node2].route;
      if (route_info[ins_route].weight+demand[max_ptr->cust_num] <= capacity){
	 insert_cust(max_ptr->cust_num, tour, max_ptr->node1,
		     max_ptr->node2, ins_route, demand, route_info);
	 head = update_savings(p, head, max_ptr, tour, route_info);
	 cust_num++;
      }else{
	 savings = find_new_ins_route(p, max_ptr->cust_num, tour, &v0, &v1,
				      route_info);
	 if (savings == -MAXINT){
	    head = NULL;
	    continue;
	 }
	 head = heap_insert(head, max_ptr->cust_num, savings, v0, v1);
      }
      FREE(max_ptr);
   }

   tour[0].next = route_info[1].first;
   
   /*------------------------------------------------------------------------*\
   | This loop points the last node of each route to the first node of the    |
   | next route. At the end of this procedure, the last node of each route is |
   | pointing at the depot, which is not what we want.                        |
   \*------------------------------------------------------------------------*/

   for (cur_route = 1; cur_route < numroutes; cur_route++)
      tour[route_info[cur_route].last].next = route_info[cur_route+1].first;
      
   if (cust_num != vertnum){
      printf(
       "\n\nError: customer cannot be inserted on any route .... aborting!\n");
      p->cur_tour->cost = 0;
   }else{
      p->cur_tour->cost = compute_tour_cost(&p->dist, tour);
   }
   
   /*-----------------------------------------------------------------------*\
   |               Transmit the tour back to the parent                      |
   \*-----------------------------------------------------------------------*/
   
   send_tour(tour, p->cur_tour->cost, p->cur_tour->numroutes, SAVINGS3,
	     used_time(&t), parent, vertnum, 0, NULL);

   FREE(intour);
   FREE(nbtree);
   
   free_heur_prob(p);
   
   PVM_FUNC(r_bufid, pvm_recv(parent, YOU_CAN_DIE));
   PVM_FUNC(info, pvm_freebuf(r_bufid));
   PVM_FUNC(info, pvm_exit());
}
      
