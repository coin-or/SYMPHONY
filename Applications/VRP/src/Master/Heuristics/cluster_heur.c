#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cluster_heur.h"
#include "vrp_master_functions.h"
#include "vrp_const.h"
#include "proccomm.h"
#include "BB_macros.h"


/*===========================================================================*/

void generate_starter(int vertnum, int *starter, int *startpos,
		      int num)
{
   int i, rans = *startpos,  pos = *startpos;
   int ran;
   
   for (i=0; i<num; i++){
      do{
	 ran = (RANDOM() % vertnum)+1;
	 for (pos=*startpos; pos < rans; pos++)
	    if (starter[pos] == ran) break;
      }while (pos < rans);
      starter[rans++] = ran;
   }
   *startpos = rans;
}

/*===========================================================================*/

/*----------------------------------------------------------------------*\
| This function spawns the various cluster heuristics which are the first|
| step in the process of finding reasonable upper bounds                 |
\*----------------------------------------------------------------------*/

void cluster_heur(vrp_problem *vrp, heur_params *heur_par,
		  heurs *ch, int trials)
{
  vrp_params *par = &vrp->par;
  int *tids;
  int *sweep_tid, *savings_tid, *savings2_tid, *near_cluster_tid, *temp_tid;
  int *tsp_ni_tid, *tsp_fi_tid, *tsp_fini_tid, *savings3_tid;
  int jobs = 0, sweep_job=0, savings_jobs=0;
  int savings2_jobs = 0, savings3_jobs = 0, near_cluster_job=0;
  int tsp_ni_jobs = 0, tsp_fi_jobs = 0, tsp_fini_jobs = 0;
  int i;
  int s_bufid;
  int vertnum = vrp->vertnum;
  float lamda, mu;
  float temp_mu, temp_lamda;
  float lamda_int, mu_int;
  int grid_size, count = 0, j, l;
  int start;
  int farside = (int) (vertnum * heur_par->fini_ratio);
  int *starter;
  int rand_nums;

  tids = ch->tids = (int *) calloc (trials, sizeof(int));

  if (vrp->par.verbosity > 1)
     printf("\nNow beginning cluster heuristics ....\n\n");

  /*-----------------------------------------------------------------------*\
  |       Enlisting this process and starting the heuristics processes      |
  \*-----------------------------------------------------------------------*/

  sweep_tid = tids;

  if (heur_par->sweep_trials && vrp->dist.coordx && vrp->dist.coordy)
     sweep_job = spawn(par->executables.sweep, (char **)NULL,
		       par->debug.sweep, (char *)NULL, 1, sweep_tid);
  
  savings_tid = sweep_tid + sweep_job;
  if (heur_par->savings_par.grid_size && heur_par->savings_par.savings_trials)
     savings_jobs = spawn(par->executables.savings, (char **)NULL,
			  par->debug.savings, (char *)NULL,
			  heur_par->savings_par.savings_trials *
			  (heur_par->savings_par.grid_size *
			   heur_par->savings_par.grid_size), savings_tid);
  
  savings2_tid = savings_tid + savings_jobs;
  if (heur_par->savings_par.grid_size && heur_par->savings_par.savings2_trials)
     savings2_jobs = spawn(par->executables.savings2,
			   (char **)NULL, par->debug.savings, (char *)NULL,
			   (heur_par->savings_par.savings2_trials *
			    (heur_par->savings_par.grid_size *
			     heur_par->savings_par.grid_size)),
			   savings2_tid);
  
  savings3_tid = savings2_tid + savings2_jobs;
  if (heur_par->savings3_par.grid_size &&
      heur_par->savings3_par.savings_trials &&
      (vrp->numroutes ||
       (!vrp->numroutes && vrp->dist.coordx && vrp->dist.coordy)))
     savings3_jobs = spawn(par->executables.savings3,
			   (char **)NULL, par->debug.savings3, (char *)NULL,
			   heur_par->savings3_par.grid_size *
			   heur_par->savings3_par.grid_size,
			   savings3_tid);
  
  near_cluster_tid = savings3_tid + savings3_jobs;
  if (heur_par->near_cluster_trials &&
      (vrp->numroutes ||
       (!vrp->numroutes && vrp->dist.coordx && vrp->dist.coordy)))
     near_cluster_job = spawn(par->executables.near_cluster,(char **)NULL,
			      par->debug.near_cluster, (char *)NULL,
			      1, near_cluster_tid);
  
  tsp_ni_tid = near_cluster_tid + near_cluster_job;
  if (heur_par->tsp.ni_trials)
     tsp_ni_jobs = spawn(par->executables.tsp_ni, (char **)NULL,
			 par->debug.tsp_ni, (char *)NULL,
			 heur_par->tsp.ni_trials, tsp_ni_tid);
  
  tsp_fi_tid = tsp_ni_tid + tsp_ni_jobs;
  if (heur_par->tsp.fi_trials)
     tsp_fi_jobs = spawn(par->executables.tsp_fi, (char **)NULL,
			 par->debug.tsp_fi, (char *)NULL,
			 heur_par->tsp.fi_trials, tsp_fi_tid);
  
  tsp_fini_tid = tsp_fi_tid + tsp_fi_jobs;
  if (heur_par->tsp.fini_trials)
     tsp_fini_jobs = spawn(par->executables.tsp_fini, (char **)NULL,
			   par->debug.tsp_fini, (char *)NULL,
			   heur_par->tsp.fini_trials, tsp_fini_tid);
  
  ch->jobs = jobs = sweep_job + savings_jobs + savings2_jobs + savings3_jobs +
     near_cluster_job + tsp_ni_jobs + tsp_fi_jobs + tsp_fini_jobs;
  
  if (jobs == 0){
     fprintf(stderr, "\nNo jobs started... \n\n");
     return;
  }
  else if (vrp->par.verbosity > 2)
     printf("\n%i jobs started ....\n\n", jobs);
  
  ch->finished = (char *) calloc (jobs, sizeof(char));
  
  /*-----------------------------------------------------------------------*\
  |                  Broadcast data  to the processes                       |
  \*-----------------------------------------------------------------------*/

  broadcast(vrp, tids, jobs);

  /*-----------------------------------------------------------------------*\
  |         Broadcast number of trials to sweep and tsp routines            |
  \*-----------------------------------------------------------------------*/
  if (sweep_job){
     s_bufid = init_send(DataInPlace);
     send_int_array(&heur_par->sweep_trials, 1);
     send_msg(*sweep_tid, SWEEP_TRIALS);
  }
  if (tsp_ni_jobs+tsp_fi_jobs+tsp_fini_jobs){
     s_bufid = init_send(DataInPlace);
     send_int_array(&heur_par->tsp.num_starts, 1);
     msend_msg(tsp_ni_tid, tsp_ni_jobs+tsp_fi_jobs+tsp_fini_jobs, TSP_TRIALS);
  }
  if (savings3_jobs || near_cluster_job){
     s_bufid = init_send(DataInPlace);
     send_int_array(&vrp->numroutes, 1);
     if (savings3_jobs)
	msend_msg(savings3_tid, savings3_jobs, NUMROUTES);
     if (near_cluster_job)
	send_msg(*near_cluster_tid, NUMROUTES);
  }
  
  if (sweep_job || (!vrp->numroutes && (savings3_jobs || near_cluster_job))){
     s_bufid = init_send(DataInPlace);
     send_dbl_array(vrp->dist.coordx, vrp->vertnum);
     send_dbl_array(vrp->dist.coordy, vrp->vertnum);
     if (sweep_job)
	send_msg(*sweep_tid, COORD_DATA);
     if (!vrp->numroutes && savings3_jobs)
	msend_msg(savings3_tid, savings3_jobs, COORD_DATA);
     if (!vrp->numroutes && near_cluster_job)
	send_msg(*near_cluster_tid, COORD_DATA);
  }
    
  /*------------------------------------------------------------------------*\
  |        Broadcast the far/near ratio to the 'fini' processes              |
  \*------------------------------------------------------------------------*/
  if (tsp_fini_jobs){
     s_bufid = init_send(DataInPlace);
     send_int_array(&farside, 1);
     msend_msg(tsp_fini_tid, tsp_fini_jobs, VRP_DATA);
  }

  /*------------------------------------------------------------------------*\
  |           Generate the random starting points and broadcast              |
  \*------------------------------------------------------------------------*/
  SRANDOM(vrp->par.rand_seed[0]);
  starter = ch->starter = 
     (int *) calloc (tsp_fi_jobs+tsp_ni_jobs+tsp_fini_jobs,sizeof(int));
  rand_nums = 0;
  generate_starter(vertnum-1, starter, &rand_nums, 
		   tsp_fi_jobs+tsp_ni_jobs+tsp_fini_jobs);
  if (tsp_fi_jobs)
     starter[0] = vertnum;
  if (tsp_ni_jobs)
     starter[tsp_fi_jobs] = vertnum;
  if (tsp_fini_jobs)
     starter[tsp_fi_jobs+tsp_ni_jobs] = vertnum;
  
  for (i=0; i<tsp_fi_jobs+tsp_ni_jobs+tsp_fini_jobs; i++){
     s_bufid = init_send(DataInPlace);
     send_int_array(starter+i, 1);
     send_msg(tsp_ni_tid[i], HEUR_START_POINT);
  }
  
  /*-----------------------------------------------------------------------*\
  |                Broadcast parameters to savings                          |
  \*-----------------------------------------------------------------------*/
  
  /*--------------------------------------------------------------------*\
  | Here is where the parameter settings are given to the various        |
  | savings proceses. The parameters mu and lamda from the parameter     |
  | file indicate maximum values for these particular parameters. Then   |
  | depending on the grid-size that is specified, I partition the        |
  | square (0, mu) X (0, lamda) into a grid and send various             |
  | combinations of parameters to the various processes. Also, I give    |
  | the processes different starting point rules. If I send FAR_INS to   |
  | a process, then the frathest customer from the depot is always       |
  | used to intialize new routes. Otherwise, I send a random seed that   |
  | can be used to randonly select the next customer to be added to a    |
  | route. The first "block" of trials receives the FAR_INS signal.      |
  | The remaining blocks, launched if savings_trials >1, receive         |
  | random seeds                                                         |
  \*--------------------------------------------------------------------*/

  if ( savings_jobs ){
     temp_tid = savings_tid;
     grid_size = heur_par->savings_par.grid_size;
     lamda = heur_par->savings_par.lamda;
     mu = heur_par->savings_par.mu;
     if (grid_size == 1){
	lamda_int = 0;
	mu_int = 0;
     }
     else{
	lamda_int = ((float)lamda)/((float)(grid_size - 1));
	mu_int = ((float)mu)/((float)(grid_size -1));
     }
     for (l=0; l<heur_par->savings_par.savings_trials; l++){
	switch(l){
	 case 0: start = FAR_INS;
	   break;  /*send FAR_INS rule to the first "block"*/
	 default: start = vrp->par.rand_seed[(l-1)%6];
	}   /*send random seeds to the remaining blocks*/
	for (i=0; i < grid_size; i++){
	   temp_lamda = i * lamda_int;
	   if (i == grid_size-1) temp_lamda = (float) lamda; 
	   for (j=0; j < grid_size; j++){
	      temp_mu = j * mu_int;
	      if (j == grid_size-1) temp_mu = (float) mu; 
	      s_bufid = init_send(DataInPlace);
	      send_float_array(&temp_mu, 1);
	      send_float_array(&temp_lamda, 1);
	      send_int_array(&start, 1);
	      msend_msg(savings_tid+count, 1, SAVINGS_DATA);
	      count++;
	      if (count > savings_jobs) break;
	   }
	   if (count > savings_jobs) break;
	}
     }
  }
  count = 0; /*This is the same as above except for savings2 instead *\
	     \*of savings.                                           */
  if ( savings2_jobs ){
     temp_tid = savings2_tid;
     grid_size = heur_par->savings_par.grid_size;
     lamda = heur_par->savings_par.lamda;
     mu = heur_par->savings_par.mu;
     if (grid_size == 1){
	lamda_int = 0;
	mu_int = 0;
     }
     else{
	lamda_int = ((float)lamda)/((float)(grid_size - 1));
	mu_int = ((float)mu)/((float)(grid_size -1));
     }
     for(l = 0; l < heur_par->savings_par.savings2_trials; l++){
	switch(l){
	 case 0: start = FAR_INS;
	   break;
	 default: start = vrp->par.rand_seed[(l-1)%6];
	}
	for (i = 0; i < grid_size; i++){
	   temp_lamda = i * lamda_int;
	   if (i == grid_size-1) temp_lamda = (float) lamda; 
	   for (j = 0; j < grid_size; j++){
	      temp_mu = j * mu_int;
	      if (j == grid_size-1) temp_mu = (float) mu; 
	      s_bufid = init_send(DataInPlace);
	      send_float_array(&temp_mu, 1);
	      send_float_array(&temp_lamda, 1);
	      send_int_array(&start, 1);
	      msend_msg(savings2_tid+count, 1, SAVINGS_DATA);
	      count++;
	      if (count > savings2_jobs) break;
	   }
	   if (count > savings2_jobs) break;
	}
     }
  }
  if ( savings3_jobs ){
     count = 0;
     temp_tid = savings3_tid;
     grid_size = heur_par->savings3_par.grid_size;
     lamda = heur_par->savings3_par.lamda;
     mu = heur_par->savings3_par.mu;
     if (grid_size == 1){
	lamda_int = 0;
	mu_int = 0;
     }
     else{
	lamda_int = ((float)lamda)/((float)(grid_size - 1));
	mu_int = ((float)mu)/((float)(grid_size -1));
     }
     for (i=0; i < grid_size; i++){
	temp_lamda = i * lamda_int;
	if (i == grid_size-1) temp_lamda = (float) lamda; 
	for (j=0; j < grid_size; j++){
	   temp_mu = j * mu_int;
	   if (j == grid_size-1) temp_mu = (float) mu; 
	   s_bufid = init_send(DataInPlace);
	   send_float_array(&temp_mu, 1);
	   send_float_array(&temp_lamda, 1);
	   send_int_array(&start, 1);
	   msend_msg(savings3_tid+count, 1, SAVINGS_DATA);
	   count++;
	   if (count > savings3_jobs) break;
	}
	if (count > savings3_jobs) break;
     }
  }
}

