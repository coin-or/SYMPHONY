#include <malloc.h>
#include <stdio.h>

#include "BB_constants.h"
#include "BB_macros.h"
#include "start_heurs.h"
#include "cluster_heur.h"
#include "route_heur.h"
#include "exchange_heur.h"
#include "small_graph.h"
#include "receive_rout.h"
#include "lower_bound.h"
#include "timemeas.h"
#include "vrp_const.h"
#include "dg_params.h"
#include "vrp_dg.h"

void start_heurs(vrp_problem *vrp, heur_params *heur_par, lb_params *lb_par,
		 double *ub, char windows)
{
  int pos, last;
  int trials = 0, i;
  best_tours *tours = NULL;
  int *tourorder = NULL, ub_pos;
  heurs *heur_out = NULL;
  double heurtime = 0, lbtime = 0, t = 0;

  (void)used_time(&t);

  /*-----------------------------------------------------------------------*\
  |                   start the cluster heuristics                          |
  \*-----------------------------------------------------------------------*/

  heur_out = (heurs *) calloc (1, sizeof(heurs));
     
  if (vrp->par.do_heuristics){
    if (heur_par->sweep_trials && vrp->dist.coordx && vrp->dist.coordy){
      trials += 1;
    }
    if (heur_par->savings3_par.savings_trials && (vrp->numroutes ||
	(!vrp->numroutes && vrp->dist.coordx && vrp->dist.coordy))){
      trials += heur_par->savings3_par.grid_size
	* heur_par->savings3_par.grid_size;
    }
    if (heur_par->near_cluster_trials && (vrp->numroutes ||
        (!vrp->numroutes && vrp->dist.coordx && vrp->dist.coordy))){
      trials += 1;
    }
     
    trials += ((heur_par->savings_par.savings_trials +
		heur_par->savings_par.savings2_trials)
	       * (heur_par->savings_par.grid_size *
		  heur_par->savings_par.grid_size) +
	       heur_par->tsp.ni_trials + heur_par->tsp.fi_trials +
	       heur_par->tsp.fini_trials);
  }else{
     trials = 0;
  }
     
  if (trials)
    cluster_heur(vrp, heur_par, heur_out, trials);

  /*-----------------------------------------------------------------------*\
  |                     receive the heuristical tours                       |
  \*-----------------------------------------------------------------------*/

  if (trials){
    tours = vrp->tours;
    vrp->tours[0].cost = MAXINT;
    vrp->tours[0].tour = (_node *) calloc
      (vrp->par.tours_to_keep, vrp->vertnum * sizeof(_node));
    for (pos=1; pos<vrp->par.tours_to_keep; pos++){
      tours[pos].tour = tours[pos-1].tour + vrp->vertnum;
      tours[pos].cost = MAXINT;
    }
    tourorder = vrp->tourorder = (int *)
      calloc(vrp->par.tours_to_keep+1, sizeof(int));
    vrp->tournum = -1;

    heurtime += receive_tours(vrp, heur_out, &last, TRUE, FALSE, FALSE,
			      windows);
                                             /*receives the tours output by *\
				             | the clustering processes and  |
				             | orders them , keeping the     |
				             \*vrp->par.tours_to_keep best  */
  }

  /*-----------------------------------------------------------------------*\
  |       Start the route heuristics                                        |
  \*-----------------------------------------------------------------------*/
  
  trials = trials ? MIN(vrp->tournum+1, heur_par->route_opt1) : 0;
  if (trials){
    route_heur(vrp, heur_par, heur_out, trials);
    heurtime += receive_tours(vrp, heur_out, &last, TRUE, FALSE, FALSE,
			      windows);
  }

  /*-----------------------------------------------------------------------*\
  |       Start the exchange heuristics                                     |
  \*-----------------------------------------------------------------------*/

  trials = trials ? MIN(vrp->tournum+1, heur_par->exchange) : 0;
  if (trials){
    exchange_heur(vrp, heur_out, trials, FIRST_SET);
    heurtime += receive_tours(vrp, heur_out, &last, TRUE, FALSE, FALSE,
			      windows);
  }

  /*-----------------------------------------------------------------------*\
  |       Run the route heuristics again                                    |
  \*-----------------------------------------------------------------------*/

  trials = trials ? MIN(vrp->tournum+1, heur_par->route_opt2) : 0;
  if (trials){
    route_heur(vrp, heur_par, heur_out, trials);
    heurtime += receive_tours(vrp, heur_out, &last, TRUE, FALSE, FALSE,
			      windows);
  }

  /*-----------------------------------------------------------------------*\
  |       Start the second set of exchange heuristics                       |
  \*-----------------------------------------------------------------------*/
  
  trials = trials ? MIN(vrp->tournum+1, heur_par->exchange2) : 0;
  if (trials){
    exchange_heur(vrp, heur_out, trials, SECOND_SET);
    heurtime += receive_tours(vrp, heur_out, &last, TRUE, FALSE, FALSE,
			      windows);
  }
  /*-----------------------------------------------------------------------*\
  |       Run route heuristics once more                                    |
  \*-----------------------------------------------------------------------*/

  trials = trials ? MIN(vrp->tournum+1, heur_par->route_opt3) : 0;
  if (trials){
    route_heur(vrp, heur_par, heur_out, trials);
    make_small_graph(vrp, 2*heur_out->jobs*vrp->vertnum); 
    heurtime += receive_tours(vrp, heur_out, &last, TRUE, FALSE, TRUE,
			      windows);
  }

  /*-----------------------------------------------------------------------*\
  |      Compute the upper and lower bounds                                 |
  \*-----------------------------------------------------------------------*/
  vrp->lb->lower_bound = (double) 0;
  ub_pos = -1;
  
  if (vrp->tournum>=0 && !*ub){
     if (!vrp->numroutes){
	vrp->numroutes = tours[tourorder[0]].numroutes;
	*ub = (double) tours[tourorder[0]].cost;
	ub_pos = 0;
     }
     else{
	for (i=0; i<= vrp->tournum; i++)
	   if (tours[tourorder[i]].numroutes == vrp->numroutes){
	      *ub = (double) tours[tourorder[i]].cost;
	      ub_pos = i;
	      break;
	   }
     }
     if (ub_pos >= 0){
	/*vrp->cur_tour = &tours[tourorder[ub_pos]];*/
	vrp->cur_tour->numroutes = tours[tourorder[ub_pos]].numroutes;
	vrp->cur_tour->cost = tours[tourorder[ub_pos]].cost;
	memcpy((char *)vrp->cur_tour->tour,
	       (char *)tours[tourorder[ub_pos]].tour,
	       vrp->vertnum*sizeof(_node));
	vrp->cur_tour->algorithm = 0;
	vrp->cur_tour->solve_time = 0;
	FREE(vrp->cur_tour->route_info);
     }else{
	printf("No tours found with %i routes -- no upper bound found\n\n",
	       vrp->numroutes);
     }
  }
  
  vrp->bd_time.ub_overhead = used_time(&t);
  vrp->bd_time.ub_heurtime = heurtime;

  if (vrp->numroutes && lb_par->lower_bound){
    lower_bound(vrp, lb_par, heur_out, (int)(*ub));
    lbtime += receive_lbs(vrp, heur_out, 1, vrp->numroutes);
  }

  vrp->bd_time.lb_overhead += used_time(&t);
  vrp->bd_time.lb_heurtime = lbtime;

  if (heur_out) free((char *)heur_out);

  if (ub_pos>=0 && vrp->par.verbosity>0){
    printf("\nThe best upper bound found with %i trucks is %i\n\n",
	    vrp->numroutes, (int)(*ub));
    if (vrp->par.verbosity >1){
       if (windows){
	  /*if (!vrp->window)
	     vrp->window = init_win(vrp->par.executables.winprog,
				    vrp->par.debug.winprog,
				    vrp->posx, vrp->posy, vrp->vertnum);
	  vrp->window->windata.draw = TRUE;
	  disp_tour(vrp->window, tours[tourorder[ub_pos]].tour, TRUE);*/
       }
    }
  }
  else if (*ub && vrp->par.verbosity > 2){
    printf("\nUsing upper bound of %i given in parameter file\n\n",
	   (int)(*ub));
  }

  if (vrp->lb->lower_bound && vrp->par.verbosity>0){
    printf("\nThe best lower bound found with %i trucks is %.1f\n\n",
	   vrp->numroutes, vrp->lb->lower_bound);
#if 0
    if (vrp->par.verbosity >1 && vrp->dg_id){
      char name[20] = {"Lower Bound"};
      disp_lb(vrp->dg_id, TRUE, name, vrp->lb->tree, vrp->lb->best_edges,
	      vrp->vertnum, vrp->numroutes, CTOI_WAIT_FOR_CLICK_AND_REPORT);
    }
#endif
  }

  /*
     p->comp_times.bc_overhead = used_time(&t);
   */
}

