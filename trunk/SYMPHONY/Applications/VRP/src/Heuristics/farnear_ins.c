/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/* This file was modified by Ali Pilatin January, 2005 (alp8@lehigh.edu)     */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "farnear_ins.h"
#include <stdio.h>
#ifndef COMPAR
#define COMPAR
static int compar(const void *elem1, const void *elem2){
   return(((neighbor *)elem1)->cost - ((neighbor *)elem2)->cost);
}
#endif
void farnear_ins(int parent, heur_prob *p)
{
  printf("\nIn farnear_ins....\n\n");
   int mytid, info, r_bufid;
   int farnode, *starter;
   int *intour;
   int i, last, cost, numroutes;
   neighbor *nbtree;
   _node *tour;
   route_data *route_info;
   int cur_route, start;
   best_tours *tours;
   double t=0;

   mytid = pvm_mytid();   
   (void) used_time(&t);
   
   tours = p->cur_tour = (best_tours *) calloc (1, sizeof(best_tours));

   
   /*-----------------------------------------------------------------------*\
   |                    Receive the VRP data                                 |
   \*-----------------------------------------------------------------------*/	

   PVM_FUNC(r_bufid, pvm_recv(-1, ROUTE_FNINS_VRP_DATA));
   PVM_FUNC(info, pvm_upkbyte((char *)tours, sizeof(best_tours), 1));
   tour = p->cur_tour->tour = (_node *) calloc (p->vertnum, sizeof(_node));
   PVM_FUNC(info, pvm_upkbyte((char *)tour, (p->vertnum)*sizeof(_node), 1));
   numroutes = p->cur_tour->numroutes;
   starter = (int *) calloc (numroutes, sizeof(int));
   route_info = tours->route_info
      = (route_data *) calloc (numroutes+1, sizeof(route_data));
   
   PVM_FUNC(r_bufid, pvm_recv(-1, ROUTE_FNINS_START_RULE));
   PVM_FUNC(info, pvm_upkint(&start, 1, 1));/*receive the start
							    rule*/

   PVM_FUNC(r_bufid, pvm_recv(parent, FINI_RATIO));
   PVM_FUNC(info, pvm_upkfloat(&p->par.fini_ratio, 1, 1));
   
   if (start != FAR_INS) srand(start); /*if the start rule is random, then*\
				       \*initialize the random number gen.*/
   starters(p, starter, route_info, start);/*generate the route starters for*\
					   \*all the clusters.              */
   /*-----------------------------------------------------------------------*\
   |                     Allocate arrays                                     |
   \*-----------------------------------------------------------------------*/

   nbtree   = (neighbor *) malloc (p->vertnum * sizeof(neighbor));
   intour   = (int *)      calloc (p->vertnum, sizeof(int));
   
   for(cur_route = 1; cur_route<=numroutes; cur_route++){
      /*---------------------------------------------------------------------*\
      |   Find the first 'farside node with farthest insertion from 'starter' |
      \*---------------------------------------------------------------------*/
      if (route_info[cur_route].numcust <= 1) continue;
      cost = 0;
      last = 0;
      intour[0] = 0;
      intour[starter[cur_route-1]] = IN_TOUR;
      fi_insert_edges(p, starter[cur_route-1], nbtree, intour, &last, tour,
		      cur_route);
      farnode = farthest(nbtree, intour, &last);
      intour[farnode] = IN_TOUR;
      fi_insert_edges(p, farnode, nbtree, intour, &last, tour, cur_route);
      tour[starter[cur_route-1]].next = farnode;
      tour[farnode].next = starter[cur_route-1];
      if (starter[cur_route - 1] == 0)
	 route_info[cur_route].first = route_info[cur_route].last = farnode;
      if (farnode == 0)
	 route_info[cur_route].first = route_info[cur_route].last
	    = starter[cur_route-1];
      
      cost = 2 * ICOST(&p->dist, starter[cur_route-1], farnode);
      cost = farthest_ins_from_to(p, tour, cost, 
				  2, (int)ceil(route_info[cur_route].
					       numcust * p->par.fini_ratio)+1,
				  starter[cur_route-1], nbtree, intour,
				  &last, route_info, cur_route);

      /*---------------------------------------------------------------------*\
      | Order the elements in nbtree (and fix the intour values) so after that|
      | nbtree is suitable to continue with nearest insertion.                |
      \*---------------------------------------------------------------------*/
      qsort((char *)(nbtree+1), last, sizeof(neighbor), compar);
      for (i=1; i<=last; i++)
	 intour[nbtree[i].nbor] = i;

      /*---------------------------------------------------------------------*\
      |              Continue with nearest insertion                          |
      \*---------------------------------------------------------------------*/
      cost = nearest_ins_from_to(p, tour, cost, 
				 MAX(((int)ceil(route_info[cur_route].numcust*
					   p->par.fini_ratio)+1), 2),
				 route_info[cur_route].numcust+1,
				 starter[cur_route-1],
				 nbtree, intour, &last, route_info, cur_route);
    
      route_info[cur_route].cost = cost;
   }
   
   tour[0].next = route_info[1].first;
   
   /*------------------------------------------------------------------------*\
   | This loop points the last node of each route to the first node of the    |
   | next route. At the end of this procedure, the last node of each route is |
   | pointing at the depot, which is not what we want.                        |
   \*------------------------------------------------------------------------*/
   for (cur_route = 1; cur_route< numroutes; cur_route++)
      tour[route_info[cur_route].last].next = route_info[cur_route+1].first;
   
   cost = compute_tour_cost(&p->dist, tour);
   
   /*-----------------------------------------------------------------------*\
   |             Transmit the tour back to the parent                        |
   \*-----------------------------------------------------------------------*/
  
   send_tour(tour, cost, numroutes, tours->algorithm, used_time(&t), parent,
	     p->vertnum, 1, route_info);
   
   if ( nbtree ) free ((char *) nbtree);
   if ( intour ) free ((char *) intour);
   if ( starter ) free ((char *) starter);
   
   free_heur_prob(p);
   
} 
