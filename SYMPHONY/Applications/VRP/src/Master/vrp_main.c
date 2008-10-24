/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the master process.
\*===========================================================================*/

#include "symphony.h"
#include "sym_master.h"
#include "vrp_types.h"
#include "vrp_const.h"

#define MAX_FN_LENGTH  100

int vrp_test(sym_environment *env, int argc, char **argv);
int vrp_load_problem(sym_environment *env, vrp_problem *vrp);

int main(int argc, char **argv)
{

   //if (argc != 4) {
   // printf("Usage: ./vrp -F file_name -N numroutes vrp_param_file");
   // abort();
   // }

   /* set initial params */

   int scenario_cnt = 10; /* scneratio count */
   double prob_show_up = 0.9; /*probability of a customer to show up */
   double prob_demand_change = 0.3; /* probability that a customer
				       change their demand */
   double demand_deviation = 0.2; /* deviation ratio in customer demand */
   double ws_node_level = 0.2; /* node level to start for ws */
   double time_limit = 10000; /* total time limit */
   
   /* read param file if given */

   char paramfile[MAX_FN_LENGTH+1];
   FILE *f;
   char line[100];
   char key[MAX_LINE_LENGTH +1], value[MAX_LINE_LENGTH +1];
   int tmp_id;
   double tmp_val; 
   
   strcpy(paramfile, "");
   if(argc > 5){
      strcpy(paramfile, argv[5]);
   }

   if(strcmp(paramfile, "") != 0){
      if ((f = fopen(paramfile, "r")) == NULL){
	 printf("param file can't be opened...exiting...\n");
	 exit(0);
      }
      while (NULL != fgets(line, MAX_LINE_LENGTH, f)){
	 strcpy(key,"");
	 sscanf(line,"%s%s", key, value);
	 if(strcmp(key, "scenario_cnt") == 0){
	    sscanf(value, "%i", &tmp_id);
	    scenario_cnt = tmp_id;
	 } else if(strcmp(key, "prob_show_up") == 0){
	    sscanf(value, "%lf", &tmp_val);
	    prob_show_up = tmp_val;
	 } else if(strcmp(key, "prob_demand_change") == 0){
	    sscanf(value, "%lf", &tmp_val);
	    prob_demand_change = tmp_val;
	 } else if(strcmp(key, "demand_deviation") == 0){
	    sscanf(value, "%lf", &tmp_val);
	    demand_deviation = tmp_val;
	 } else if(strcmp(key, "ws_node_level") == 0){
	    sscanf(value, "%lf", &tmp_val);
	    ws_node_level = tmp_val;
	 } else if(strcmp(key, "time_limit") == 0){
	    sscanf(value, "%lf", &tmp_val);
	    time_limit = tmp_val;
	 }
      }
      fclose(f);
   }
   /* timing stuff */
   double start_time = wall_clock(NULL), time_remained = 0.0;
   double ws_time, solve_time; 
   double init_solve_time, total_ws_time = 0.0, total_solve_time = 0.0; 
   double tol = 1e-06;
   
   /* initialize symphony*/
   sym_environment *env = sym_open_environment();
   version();
   sym_parse_command_line(env, argc, argv);
   sym_load_problem(env);
   sym_find_initial_bounds(env);
   
   //  vrp_load_problem(env, (vrp_problem *)env->user);
   
   sym_set_str_param(env, "lp_executable_name", "vrp_lp_cg");
   sym_set_str_param(env, "cp_executable_name", "vrp_cp");
   sym_set_int_param(env, "generate_cgl_cuts", FALSE);
   sym_set_int_param(env, "fp_enabled", -1);
   //sym_set_int_param(env, "strong_branching_cand_num_min", 20);
   //sym_set_int_param(env, "strong_branching_cand_num_max", 20);
   sym_set_int_param(env, "strong_branching_cand_num", 20);
   sym_set_int_param(env, "max_presolve_iter", 40);
   sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
   sym_set_int_param(env, "keep_warm_start", TRUE);      
   sym_set_dbl_param(env, "warm_start_node_level_ratio", ws_node_level);
      
   /* initial solve */
   double mark_time = wall_clock(NULL);
   sym_set_dbl_param(env, "time_limit", time_limit);
   sym_solve(env);
   init_solve_time = wall_clock(NULL) - mark_time;
   time_remained = time_limit - init_solve_time;

   if(time_remained < 5){
      printf("probably couldnt solve initial problem in given time...\n");
      exit(0);
   }

   
   /*get ws_info */
   warm_start_desc *ws = sym_get_warm_start(env, FALSE);
   sym_set_int_param(env, "keep_warm_start", FALSE);      
   
   /* get vrp data */   
   vrp_problem *vrp = (vrp_problem *)env->user;
   int i, vertnum = vrp->vertnum;
   int * act_demand = (int *)malloc(ISIZE*vertnum);
   int * new_demand = (int *)malloc(ISIZE*vertnum);
   memcpy(act_demand, vrp->demand, ISIZE*vertnum);
   int act_numroutes = vrp->numroutes;
   int capacity = vrp->capacity;
   int deviation;

   //int * cust_del_ind = (int *)malloc(ISIZE*vertnum);
   int del_num; 
   
   double prob, ws_objval, solve_objval;
   int random_int, batch = 1, demand_chg_cnt = 0, done_scn_cnt = 0, scn_num;  
   char ws_infeas, solve_infeas;
   for(scn_num = 0; scn_num < scenario_cnt; scn_num++){
      del_num = 0;
      new_demand[0] = 0;   
      srand ( time(NULL) );      
      for(i = 1; i < vertnum; i++){
	 prob = (double) rand()/(RAND_MAX+1.0);
	 deviation = (int)(vrp->demand[i] * demand_deviation);
	 if(deviation < 1){
	    deviation = 1;
	 }
	 if(prob < prob_show_up){
	    prob = (double) rand()/(RAND_MAX+1.0);
	    if(prob < prob_demand_change){
	       demand_chg_cnt++;
	       random_int = int((2*deviation+1)*(rand()/(RAND_MAX + 1.0)));
	       new_demand[i] = act_demand[i] + (random_int - deviation) *
		  batch;
	       if(new_demand[i] < 1)
		  new_demand[i] = batch;
	       else if(new_demand[i] > capacity)
		  new_demand[i] = capacity;
	    }else{
	       new_demand[i] = act_demand[i];
	    }
	 }else{
	    new_demand[i] = capacity;	    
	    //cust_del_ind[del_num] = i;
	    del_num++;
	 }
	 new_demand[0] += new_demand[i];
      }
      
      memcpy(vrp->demand, new_demand, ISIZE*vertnum);
      vrp->numroutes += del_num;

      printf("scenario %i\n", scn_num);
      printf("total # of cancelled customers %i\n", del_num);	
      printf("total # of updated demands %i\n", demand_chg_cnt);
      printf("total demand: %i\n", vrp->demand[0]);
      printf("available demand: %i\n", vrp->numroutes * capacity);      
      
      env->mip->change_num = 1;
      env->mip->change_type[0] = RHS_CHANGED;     
      
      //  sym_set_int_param(env, "verbosity", 3);
      mark_time = wall_clock(NULL);     
      sym_set_dbl_param(env, "time_limit", time_remained);      
      ws_infeas = FALSE;
      sym_set_warm_start(env, ws);
      sym_warm_solve(env);

      if(sym_is_proven_optimal(env)){
	 sym_get_obj_val(env, &ws_objval);
      }else{
	 ws_infeas = TRUE;
      }

      ws_time = wall_clock(NULL) - mark_time;
      time_remained -= ws_time;
      if(time_remained < 5){
	 break;
      }
      
      mark_time = wall_clock(NULL);
      //sym_set_int_param(env, "verbosity", 0);
      solve_infeas = FALSE;
      sym_set_dbl_param(env, "time_limit", time_remained);      
      sym_solve(env);

      if(sym_is_proven_optimal(env)){
	 sym_get_obj_val(env, &solve_objval);
      }else{
	 solve_infeas = TRUE;
      }      

      solve_time = wall_clock(NULL) - mark_time;
      time_remained -= solve_time;
      if(time_remained < 5){
	 break;
      }      
      /* check if correct optimal solution */
      if(ws_infeas != solve_infeas){
	 printf("feasibility error in ws... exiting...\n");
	 printf("scn_cnt: %i \t time remained: %f", done_scn_cnt,
		time_remained);
	 exit(0);	 
      }else{
	 if(!ws_infeas){
	    if(ws_objval > solve_objval + tol ||
	       ws_objval < solve_objval - tol){
	       printf("diff objvals...exiting...\n");
	       printf("ws_objval: %f \t solve_objval: %f",
		      ws_objval, solve_objval);
	       printf("scn_cnt: %i \t time remained: %f", done_scn_cnt,
		      time_remained);		      
	       exit(0);	 
	    }
	 }
      }

      /* so successful run on this scenario */
      /* get back */
      memcpy(vrp->demand, act_demand, ISIZE*vertnum); 
      vrp->numroutes = act_numroutes;	 
      done_scn_cnt++;
      total_ws_time += ws_time;
      total_solve_time += solve_time;
   }
   
   if(done_scn_cnt > 0){
      printf("done scenario count %i\n", done_scn_cnt);
      printf("total time: %f\n", wall_clock(NULL) - start_time);
      printf("initial time: %f\n", init_solve_time);
      printf("total solve time: %f\n", total_solve_time);
      printf("total warm_start time: %f\n", total_ws_time);
   }else{
      printf("probably didnt have enough time to solve the problems...\n");
      printf("initial time: %f\n", init_solve_time);      
   }

   FREE(act_demand);
   FREE(new_demand);
   
   sym_close_environment(env);
   
   return(0);
}

/*===========================================================================*/
/*===========================================================================*/

int vrp_load_problem(sym_environment *env, vrp_problem *vrp){

   //vrp_lp_problem *vrp = (vrp_lp_problem *)user;
   int *costs = vrp->dist.cost;
   int *edges = vrp->edges;
   int i;
   //int total_edgenum = vrp->vertnum*(vrp->vertnum-1)/2;
   int cap_check = vrp->capacity*(vrp->numroutes-1);
   int total_demand = vrp->demand[0];
   MIPdesc *mip = env->mip;

   node_desc *desc = env->rootdesc;

   int bvarnum = env->base->varnum;
   int bcutnum = env->base->cutnum;

   //int *d_uind = NULL, *d_cind = NULL; /* just to keep gcc quiet */

   //int user_res;
   int *indices;

   mip->n = bvarnum + desc->uind.size;
   mip->m = bcutnum + desc->cutind.size;

   indices = (int *) malloc(mip->n * ISIZE);
   if (bvarnum > 0){
      for (i = bvarnum - 1; i >= 0; i--){
	 indices[i] = env->base->userind[i];
      }
   }
   
   if (desc->uind.size > 0){ /* fill up the rest of lp_data->vars */
      for (i = desc->uind.size - 1; i >= 0; i--){
	 indices[i+bvarnum] = desc->uind.list[i];
      }
   }
      
   //int n = mip->n;
   //int m = mip->m;
   
   /* set up the inital LP data */

   /* Our base constraints are that the sum of the weights of the
      edges adjacent to each node (except the depot) is 2. This means
      that each column will have exactly two nonzeros, one in each of
      the rows corresponding to its end points. Hence the total
      nonzeros is 2*n (n is the number of active variables). */
   mip->nz = 2 * mip->n;

   /* Estimate the maximum number of nonzeros (not needed, but helpful for
      efficiency*/
   //*maxm = 2 * mip->m;
   //*maxn = total_edgenum;
   //*maxnz = mip->nz + ((*maxm) * (*maxn) / 10);

   /* Allocate the arrays. These are owned by SYMPHONY after returning. */
   mip->matbeg  = (int *) malloc((mip->n + 1) * ISIZE);
   mip->matind  = (int *) malloc((mip->nz) * ISIZE);
   mip->matval  = (double *) malloc((mip->nz) * DSIZE);
   mip->obj     = (double *) malloc(mip->n * DSIZE);
   mip->ub      = (double *) malloc(mip->n * DSIZE);
   mip->lb      = (double *) calloc(mip->n, DSIZE); /* zero lower bounds */
   mip->rhs     = (double *) malloc(mip->m * DSIZE);
   mip->sense   = (char *) malloc(mip->m * CSIZE);
   mip->rngval  = (double *) calloc(mip->m, DSIZE);
   mip->is_int  = (char *) calloc(mip->n, CSIZE);

   /* Fill out the appropriate data structures -- each column has
      exactly two entries */
   for (i = 0; i < mip->n; i++){
      mip->obj[i]        = (double) costs[indices[i]];
      mip->matbeg[i]     = 2*i;
      mip->matval[2*i]   = 1.0;
      mip->matval[2*i+1] = 1.0;
      mip->matind[2*i]   = edges[2*indices[i]];
      mip->matind[2*i+1] = edges[2*indices[i] + 1];
      mip->is_int[i]     = TRUE;
      
      /* Set the upper and lower bounds */
      if (edges[indices[i] << 1] == 0){
	 /*This is a depot edge and we have to check what its ub should be*/
	 /*Recall vrp->demand[0] == total_demand*/
	 if (total_demand - vrp->demand[edges[(indices[i] << 1) + 1]] >
	     cap_check)
	    /*in this case, the node cannot be on its own route*/
	    mip->ub[i] = 1;
	 else
	    mip->ub[i] = 2;
      }else{
	 mip->ub[i] = 1;
      }
   }
   mip->matbeg[mip->n] = 2*mip->n;
   
   /* set the initial right hand side */
   mip->rhs[0] = 2*vrp->numroutes;
   mip->sense[0] = 'E';
   for (i = vrp->vertnum - 1; i > 0; i--){
      if(i == 50){
	 mip->rhs[i] = (double) 2; //0.0;
      }else{
	 mip->rhs[i]   = (double) 2;
      }
      mip->sense[i] = 'E';
   }
   return 0;
}

int vrp_test(sym_environment *env, int argc, char **argv)
{

   int termcode, i, file_num = 34;
   char input_files[34][MAX_FILE_NAME_LENGTH +1] = {"A/A-n34-k5",
						   "A/A-n32-k5",
						   "A/A-n33-k5",
						   "E/E-n13-k4",
						   "E/E-n22-k4",
						   "E/E-n23-k3",
						   "E/E-n30-k3",
						   "E/E-n33-k4",
						   "V/att-n48-k4",
						   "E/E-n51-k5",
						   "A/A-n33-k6",
						   "A/A-n36-k5",
						   "A/A-n37-k5",
						   "A/A-n38-k5",
						   "A/A-n39-k5",
						   "A/A-n39-k6",
						   "A/A-n45-k6",
						   "A/A-n46-k7",
						   "B/B-n31-k5",
						   "B/B-n34-k5",
						   "B/B-n35-k5",
						   "B/B-n38-k6",
						   "B/B-n39-k5",
						   "B/B-n41-k6",
						   "B/B-n43-k6",
						   "B/B-n44-k7",
						   "B/B-n45-k5",
						   "B/B-n50-k7",
						   "B/B-n51-k7",
						   "B/B-n52-k7",
						   "B/B-n56-k7",
						   "B/B-n64-k9",
						   "A/A-n48-k7",
						   "A/A-n53-k7"};   

   double sol[34] = {778, 784, 661, 247, 375, 569, 534, 835, 40002, 521, 742,
		     799, 669, 730, 822, 831, 944, 914, 672, 788, 955, 805,
		     549, 829, 742, 909, 751, 741, 1032, 747, 707, 861,
		     1073, 1010};
   
   char *input_dir = (char*)malloc(CSIZE*(MAX_FILE_NAME_LENGTH+1));
   char *infile = (char*)malloc(CSIZE*(MAX_FILE_NAME_LENGTH+1));
   char *sgfile = (char*)malloc(CSIZE*(MAX_FILE_NAME_LENGTH+1));
   double *obj_val = (double *)calloc(DSIZE,file_num);
   double tol = 1e-06, ub;

   vrp_problem *vrp = (vrp_problem *) env->user;

   if (strcmp(vrp->par.test_dir, "") == 0){ 
     strcpy(input_dir, "../../../VRPLIB");
   } else{
     strcpy(input_dir, vrp->par.test_dir);
   }

   for(i = 0; i<file_num; i++){
     
     strcpy(infile, "");
     strcpy(sgfile, "");
     sprintf(infile, "%s%s%s%s", input_dir, "/", input_files[i], ".vrp");
     sprintf(sgfile, "%s%s%s", "./small_graph/", input_files[i], ".sg");


     strcpy(vrp->par.infile, infile);
     strcpy(vrp->par.small_graph_file, sgfile);
     vrp->par.use_small_graph = LOAD_SMALL_GRAPH;

     sym_load_problem(env);
     
     sym_find_initial_bounds(env);

     printf("Solving %s...\n", input_files[i]); 
        
     sym_solve(env);
     
     sym_get_obj_val(env, &obj_val[i]);
     
     if((obj_val[i] < sol[i] + tol) && 
	(obj_val[i] > sol[i] - tol)){
       printf("Success!\n");
     } else {
       printf("Failure!(%f, %f) \n", obj_val[i], sol[i]);
     }

     if(i+1 < file_num){
       sym_close_environment(env);

       env = sym_open_environment();       
       sym_parse_command_line(env, argc, argv);
     }
     
   }
   return (0);
}

