#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

#include "proccomm.h"
#include "route_heur.h"
#include "vrp_types.h"
#include "vrp_const.h"
#include "vrp_routines.h"

/*-----------------------------------------------------------------*\
| Here we take the tours stored in p->tours and try to improve the  |
| routing of the customers using various simple tsp heuristics. The |
| process route_opt coordinates this                                |
\*-----------------------------------------------------------------*/

void route_heur(vrp_problem *vrp, heur_params *heur_par,
		heurs *rh, int trials)
{
  int *tids;
  int jobs = 0;
  int i;
  int s_bufid;
  _node *tour;
  best_tours tours;

  if (vrp->par.verbosity>1)
     printf("\nNow beginning routing heuristics ....\n\n");

  tids = rh->tids = (int *) calloc (trials, sizeof(int));

  if (trials)
     jobs = spawn(vrp->par.executables.route_opt, (char **)NULL,
		  vrp->par.debug.route_opt,
		  (char *)NULL, trials, tids);

  rh->jobs = jobs;

  if (jobs == 0){
     fprintf(stderr, "\nNo jobs started... \n\n");
     return;
  }
  else if (vrp->par.verbosity >2)
     printf("\n%i jobs started ....\n\n", jobs);
  
  rh->finished = (char *) calloc (jobs, sizeof(char));
  
  /*-----------------------------------------------------------------------*\
  |                  Broadcast data to the processes                        |
  \*-----------------------------------------------------------------------*/

  broadcast(vrp, tids, jobs);

  s_bufid = init_send(DataInPlace);
  send_float_array(&heur_par->fini_ratio, 1);
  send_int_array(&heur_par->ni_trials, 1);
  send_int_array(&heur_par->fi_trials, 1);
  send_int_array(&heur_par->fini_trials, 1);
  send_char_array((char *)&vrp->par.debug, sizeof(debugging));
  send_int_array(&vrp->par.tours_to_keep, 1);
  send_int_array(&vrp->par.time_out.ub, 1);
  send_int_array(vrp->par.rand_seed, NUM_RANDS);
  send_str(vrp->par.executables.nearest_ins);
  send_str(vrp->par.executables.farthest_ins);
  send_str(vrp->par.executables.farnear_ins);
  msend_msg(tids, jobs, VRP_DATA);

  /*-----------------------------------------------------------------------*\
  |                  Broadcast best tours to the processes                  |
  \*-----------------------------------------------------------------------*/

  for (i=0; i<jobs; i++){
    s_bufid = init_send(DataInPlace);
    tours = vrp->tours[vrp->tourorder[i%(vrp->tournum+1)]];
    tour = tours.tour;
    send_char_array((char *)&tours, sizeof(best_tours));
    send_char_array((char *)tour, (vrp->vertnum)*sizeof(_node));
    send_msg(tids[i], HEUR_TOUR);
  }
}
    
    
