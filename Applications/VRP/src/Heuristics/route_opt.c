#include <math.h>
#include <malloc.h>

#include "BB_constants.h"
#include "BB_macros.h"
#include "messages.h"
#include "vrp_const.h"
#include "vrp_common_types.h"
#include "vrp_master_functions.h"
#include "heur_routines.h"
#include "compute_cost.h"
#include "receive_rout.h"
#include "timemeas.h"
#include "proccomm.h"
#include "vrp_macros.h"

/*----------------------------------------------------------------------*\
| This program takes clusterings of nodes (determined by the initial     |
| solution which is given to the program) and tries to improve their     |
| routings using simple TSP heuristics such as far insert and near insert|
| For each cluster of customers, it takes the cheapest routing found by  |
| any of the heuristics and adds it to the tour to be sent back to the   |
| parent.                                                                |
\*----------------------------------------------------------------------*/

void main(void)
{
  vrp_problem *p;
  _node *tour;
  int i, mintour = 0, tournum;
  int *tids, *fi_tid, *ni_tid, *fini_tid;
  int jobs, fi_jobs, ni_jobs, fini_jobs;
  int trials;
  char *finished;
  int mytid, info, s_bufid, r_bufid, bytes, msgtag, parent;
  heurs *ro;
  int cur_node, cur_route, numroutes;
  int *tourorder, pos, cost, last = 0;
  best_tours *tours;
  int mincost, temp_cost, temp;
  route_data *route_info;
  double routing_time = 0, t;
  heur_params heur_par;

  (void) used_time(&t);
  
  mytid = pvm_mytid();

  p = (vrp_problem *) calloc(1, sizeof(vrp_problem));
  tours = p->cur_tour = (best_tours *) calloc (1, sizeof(best_tours));
  ro = (heurs *) calloc (1, sizeof(heurs));

  /*-----------------------------------------------------------------------*\
  |                    Receive the VRP data                                 |
  \*-----------------------------------------------------------------------*/

  PVM_FUNC(r_bufid, pvm_recv(-1, VRP_DATA));
  PVM_FUNC(info, pvm_bufinfo(r_bufid, &bytes, &msgtag, &parent));
  PVM_FUNC(info, pvm_upkint(&(p->dist.wtype), 1, 1));
  PVM_FUNC(info, pvm_upkint(&(p->vertnum), 1, 1));
  PVM_FUNC(info, pvm_upkint(&(p->depot), 1, 1));
  PVM_FUNC(info, pvm_upkint(&p->capacity, 1, 1));
  p->demand = (int *) calloc (p->vertnum, sizeof(int));
  PVM_FUNC(info, pvm_upkint(p->demand, p->vertnum, 1));
  p->edgenum = p->vertnum*(p->vertnum-1)/2;
  if (p->dist.wtype){ /* not EXPLICIT */
    p->dist.coordx = (double *) calloc(p->vertnum, sizeof(double));
    p->dist.coordy = (double *) calloc(p->vertnum, sizeof(double));
    PVM_FUNC(info, pvm_upkdouble(p->dist.coordx, (int)p->vertnum, 1));
    PVM_FUNC(info, pvm_upkdouble(p->dist.coordy, (int)p->vertnum, 1));
    if ((p->dist.wtype == _EUC_3D) || (p->dist.wtype == _MAX_3D) || 
		    (p->dist.wtype == _MAN_3D)){
      p->dist.coordz = (double *) calloc(p->vertnum, sizeof(double));
      PVM_FUNC(info, pvm_upkdouble(p->dist.coordz, (int)p->vertnum, 1));
    }
  }
  else{ /* EXPLICIT */
    p->dist.cost = (int *) malloc ((int)p->edgenum*sizeof(int));
    PVM_FUNC(info, pvm_upkint(p->dist.cost, (int)p->edgenum, 1));
  }

  PVM_FUNC(r_bufid, pvm_recv(-1, VRP_DATA));
  PVM_FUNC(info, pvm_upkfloat(&heur_par.fini_ratio, 1, 1));
  PVM_FUNC(info, pvm_upkint(&heur_par.ni_trials, 1, 1));
  PVM_FUNC(info, pvm_upkint(&heur_par.fi_trials, 1, 1));
  PVM_FUNC(info, pvm_upkint(&heur_par.fini_trials, 1, 1));
  PVM_FUNC(info, pvm_upkbyte((char *)&p->par.debug, sizeof(debugging), 1));
  PVM_FUNC(info, pvm_upkint(&p->par.tours_to_keep, 1, 1));
  PVM_FUNC(info, pvm_upkint(&p->par.time_out.ub, 1, 1));
  p->par.rand_seed  = (int *) calloc (NUM_RANDS, sizeof(int));
  PVM_FUNC(info, pvm_upkint(p->par.rand_seed, NUM_RANDS, 1));
  PVM_FUNC(info, pvm_upkstr(p->par.executables.nearest_ins));
  PVM_FUNC(info, pvm_upkstr(p->par.executables.farthest_ins));
  PVM_FUNC(info, pvm_upkstr(p->par.executables.farnear_ins));

  /*-----------------------------------------------------------------------*\
  |                     Receive the starting tour                           |
  \*-----------------------------------------------------------------------*/

  PVM_FUNC(r_bufid, pvm_recv(-1, HEUR_TOUR));
  PVM_FUNC(info, pvm_upkbyte((char *)tours, sizeof(best_tours), 1));
  tour = p->cur_tour->tour = (_node *) calloc (p->vertnum, sizeof(_node));
  PVM_FUNC(info, pvm_upkbyte((char *)tour, (p->vertnum)*sizeof(_node), 1));
  numroutes = tours->numroutes;
  tours->route_info = (route_data *) calloc (numroutes+1, sizeof(route_data));
  route_calc(&p->dist, tour, numroutes, tours->route_info, p->demand);

  trials = heur_par.fi_trials + heur_par.ni_trials +
    heur_par.fini_trials;
  tids = ro->tids = (int *) calloc (trials, sizeof(int));
  fi_jobs = ni_jobs = fini_jobs = 0;

  fi_tid = tids;
  if (heur_par.fi_trials)
    fi_jobs = pvm_spawn(p->par.executables.farthest_ins, (char **)NULL,
			p->par.debug.farthest_ins, 
			(char *)NULL, heur_par.fi_trials, fi_tid);
  
  ni_tid = fi_tid + fi_jobs;
  if ( heur_par.ni_trials )
    ni_jobs = pvm_spawn(p->par.executables.nearest_ins, (char **)NULL,
			p->par.debug.nearest_ins, 
			(char *)NULL, heur_par.ni_trials, ni_tid);

  fini_tid = ni_tid + ni_jobs;
  if ( heur_par.fini_trials )
    fini_jobs = pvm_spawn(p->par.executables.farnear_ins, (char **)NULL,
			  p->par.debug.farnear_ins, 
			  (char *)NULL, heur_par.fini_trials, fini_tid);

  ro->jobs = jobs = fi_jobs + ni_jobs + fini_jobs;

  if (jobs == 0){
    printf("No jobs started... :-(\n");
  }
  finished = ro->finished = (char *) calloc ((int)jobs, sizeof(char));

  /*-----------------------------------------------------------------------*\
  |                  Broadcast data  to the processes                       |
  \*-----------------------------------------------------------------------*/

  broadcast(p, tids, jobs);

  PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
  PVM_FUNC(info, pvm_pkbyte((char *)tours, sizeof(best_tours), 1));
  PVM_FUNC(info, pvm_pkbyte((char *)tour, (p->vertnum)*sizeof(_node), 1));
  PVM_FUNC(info, pvm_mcast(tids, jobs, VRP_DATA));

  /*--------------------------------------------------------------------*\
  | Here we broadcast the starting point rules to the various heuristics |
  | If there is only one trial of that particular heuristic, then we just|
  | use far insert as the starting point rule. For all the remaining     |
  | trials, we use random starting points.                               |
  \*--------------------------------------------------------------------*/

  temp = FAR_INS;
  if (fi_jobs){
    PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
    PVM_FUNC(info, pvm_pkint(&temp, 1, 1));
    PVM_FUNC(info, pvm_send(fi_tid[0], VRP_DATA));
  }
  if (ni_jobs){
    PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
    PVM_FUNC(info, pvm_pkint(&temp, 1, 1));
    PVM_FUNC(info, pvm_send(ni_tid[0], VRP_DATA));
  }
  if (fini_jobs){
    PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
    PVM_FUNC(info, pvm_pkint(&temp, 1, 1));
    PVM_FUNC(info, pvm_send(fini_tid[0], VRP_DATA));
  }
  temp = MAX(fi_jobs, ni_jobs);
  temp = MAX(temp, fini_jobs);
  for (i=1; i<temp; i++){
    if (fi_jobs>i){
      PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
      PVM_FUNC(info, pvm_pkint(&p->par.rand_seed[(i-1)%NUM_RANDS], 1, 1));
      PVM_FUNC(info, pvm_send(fi_tid[i], VRP_DATA));
    }
    if (ni_jobs>i){
      PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
      PVM_FUNC(info, pvm_pkint(&p->par.rand_seed[(i-1)%NUM_RANDS], 1,
                                1));
      PVM_FUNC(info, pvm_send(ni_tid[i], VRP_DATA));
    }
    if (fini_jobs>i){
      PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
      PVM_FUNC(info, pvm_pkint(&p->par.rand_seed[(i-1)%NUM_RANDS], 1,
                                1));
      PVM_FUNC(info, pvm_send(fini_tid[i], VRP_DATA));
    }
  }
  /*-----------------------------------------------------------------------*\
  |         Broadcast the far/near ratio to the 'fini' processes            |
  \*-----------------------------------------------------------------------*/

  if ( fini_jobs ){
    PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
    PVM_FUNC(info, pvm_pkfloat(&heur_par.fini_ratio, 1, 1));
    PVM_FUNC(info, pvm_mcast(fini_tid, fini_jobs, VRP_DATA));
  }

  /*-----------------------------------------------------------------------*\
  |         Receive the tours back from the processes                       |
  \*-----------------------------------------------------------------------*/

  tours = p->tours = (best_tours *) calloc (p->par.tours_to_keep, sizeof(best_tours));
  p->tours[0].cost = MAXINT;
  p->tours[0].tour = (_node *) malloc ((int)p->par.tours_to_keep * p->vertnum
				    * sizeof(_node));
  for (pos=1; pos<p->par.tours_to_keep; pos++){
    tours[pos].tour = tours[pos-1].tour + p->vertnum;
    tours[pos].cost = MAXINT;
  }
  tourorder = p->tourorder = (int *) calloc((int)p->par.tours_to_keep+1,
					    sizeof(int));

  routing_time = receive_tours(p, ro, &last, FALSE, TRUE, FALSE, FALSE);

  /*-----------------------------------------------------------------------*\
  |               Form the final tour                                       |
  \*-----------------------------------------------------------------------*/

  /*--------------------------------------------------------------------*\
  | This section of the program compares all the routings for a          |
  | particular cluster and chooses the best one to be in the final       |
  | solution.                                                            |
  \*--------------------------------------------------------------------*/
  tournum = p->tournum;
  route_info = p->cur_tour->route_info;

  cost = 0;
  cur_node = 0;
  for (cur_route=1; cur_route<=numroutes; cur_route++){
    mincost = MAXINT;
    for (i = 0; i<=tournum; i++)
      if ((temp_cost = tours[i].route_info[cur_route].cost) < mincost){
	mincost = temp_cost;
	mintour = i;
      }
    if (mincost < route_info[cur_route].cost){
      tour[cur_node].next = tours[mintour].route_info[cur_route].first;
      do{
	cur_node = tour[cur_node].next;
	tour[cur_node].next = tours[mintour].tour[cur_node].next;
      }while(cur_node != tours[mintour].route_info[cur_route].last);
      cost += tours[mintour].route_info[cur_route].cost;
    }
    else{
      tour[cur_node].next = route_info[cur_route].first;
      cost += route_info[cur_route].cost;
      cur_node = route_info[cur_route].last;
    }
  }
  
  /*-----------------------------------------------------------------------*\
  |              Transmit the final tour back to the parent                 |
  \*-----------------------------------------------------------------------*/

  send_tour (tour, cost, numroutes, p->cur_tour->algorithm,
	     used_time(&t)+routing_time,parent, p->vertnum, 0, NULL);

  if (ro->tids) free ((char *) ro->tids);
  if (ro->finished) free ((char *) ro->finished);
  if (ro) free((char *) ro);
			  
/*   free_vrp_prob(p); */

  PVM_FUNC(r_bufid, pvm_recv(parent, YOU_CAN_DIE));
  PVM_FUNC(info, pvm_freebuf(r_bufid));
  PVM_FUNC(info, pvm_exit());
}

  





