#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

#include "exchange_heur.h"
#include "vrp_routines.h"
#include "vrp_const.h"
#include "proccomm.h"

void exchange_heur(vrp_problem *vrp, heurs *eh, int trials, int which)
{
  int *tids;
  int jobs = 0;
  int i;
  int s_bufid;
  _node *tour;
  best_tours tours;

  if (vrp->par.verbosity >1){
     if (which == FIRST_SET)
	printf("\nNow beginning first set of exchange heuristics ....\n\n");
     else if (which == SECOND_SET)
	printf("\nNow beginning second set of exchange heuristics ....\n\n");
  }

  tids = eh->tids = (int *) calloc (trials, sizeof(int));

  if (trials)
    jobs = spawn(which ? vrp->par.executables.exchange:
		 vrp->par.executables.exchange2, (char **)NULL,
		 which ? vrp->par.debug.exchange:
		 vrp->par.debug.exchange2, (char *)NULL, trials,
		 tids);

  eh->jobs = jobs;

  if (jobs == 0){
    fprintf(stderr, "\nNo jobs started... \n\n");
    return;
  }
  else if (vrp->par.verbosity >2)
     printf("\n%i jobs started ....\n\n", jobs);

  eh->finished = (char *) calloc (jobs, sizeof(char));

  /*-----------------------------------------------------------------------*\
  |                  Broadcast data to the processes                        |
  \*-----------------------------------------------------------------------*/

  broadcast(vrp, tids, jobs);

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

