/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000-2010 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (ted@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/* This file is a modified version of SYMPHONY's vrp_main.c file that 
   adds the ability to determine an initial upper bound by running
   heuristics from the VRPH library. This capability can be turned off by
   setting USE_VRPH 0 below.

   Modifications by Chris Groer (cgroer@gmail.com) */

#define COMPILING_FOR_MASTER

#define USE_VRPH 0

/* We will have VRPH generate NUM_VRPH_SOLUTIONS */
#define NUM_VRPH_SOLUTIONS 10

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the master process.
\*===========================================================================*/

#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{
   OsiSymSolverInterface si;

   //si.setSymParam(OsiSymVerbosity, 3);
   si.setSymParam(OsiSymGranularity, 0.9999); 
   si.setSymParam("generate_cgl_cuts", FALSE);
   si.setSymParam("lp_executable_name", "vrp_lp_cg");
   si.setSymParam("cp_executable_name", "vrp_cp");
   
   /* Parse the command line */
   si.parseCommandLine(argc, argv);
   
   /* Read in the problem */
   si.loadProblem();

   /* Find a priori problem bounds */
   si.findInitialBounds();
   
   /* Solve the problem */
   si.branchAndBound();

   return(0);
}

#else

#include "symphony.h"
#include "sym_master.h"
#include "vrp_types.h"
#include "vrp_const.h"

#if USE_VRPH
#include "VRPH.h"
#endif

int vrp_test(sym_environment *env, int argc, char **argv);

int main(int argc, char **argv)
{

#if USE_VRPH
   int i;
   double best_sol=VRP_INFINITY;
   int best_sol_buff[500];
   clock_t start, vrph_stop, final_stop;
#endif

   vrp_problem *vrp;

   sym_environment *env = sym_open_environment();

   version();

   sym_parse_command_line(env, argc, argv);

   sym_get_user_data(env, (void**)&vrp);
   
#if USE_VRPH
   
   start=clock();

   /* Get the size of the problem in the input file */
   int n=VRPGetDimension(vrp->par.infile);

   /* Declare a VRP object of size n */
   VRP V(n);

   /* Declare a ClarkeWright object of size n */
   ClarkeWright CW(n);
 
   /* Populate the VRP object with the input file */
   V.read_TSPLIB_file(vrp->par.infile);

   /* Now create NUM_VRPH_SOLUTIONS solutions using VRPH and set the
      upper bound to the best solution discovered */

   for (i = 0; i < NUM_VRPH_SOLUTIONS; i++){

      /* Create default routes - each customer on its own route */
      V.create_default_routes();

      /* Create a new random feasible solution with Clarke Wright */

      CW.Construct(&V, .5+ lcgrand(1),false);

      /* Improve it with the RTR heuristic */
      V.RTR_solve(ONE_POINT_MOVE+TWO_POINT_MOVE+TWO_OPT+THREE_OPT,
		  30,5,1,.01,25,VRPH_LI_PERTURB,VRPH_BEST_ACCEPT,false);

      if (V.get_total_route_length()-V.get_total_service_time()<best_sol){
	 best_sol=V.get_total_route_length()-V.get_total_service_time();
	 V.export_canonical_solution_buff(best_sol_buff);
      }
      
      /* Reset VRPH's internal data structures */
      V.reset();
   }
   vrph_stop=clock();

   /* Import the best solution and display it - if SYMPHONY claims an
      infeasibility because the VRPH solution is optimal, we wouldn't see it
      otherwise! */
   printf("VRPH set SYMPHONY upper bound to %f based on solution:\n",best_sol);
   V.import_solution_buff(best_sol_buff);
   V.summary();

   /* Set the upper bound using VRPH solution by accessing SYMPHONY's
      internal data structures */
   env->has_ub=1;
   env->ub=best_sol;

#if 0
   /* Note that this might be incorrect if the VRPH solution is not optimal
      So the # of trucks still needs to be passed in on the command line! */
   vrp->numroutes=V.get_total_number_of_routes();
#endif
#endif

   
   if (vrp->par.test){

#if USE_VRPH
      printf("Can't run in test mode with VRPH. Exiting...");
#else
      vrp_test(env, argc, argv);
#endif
   }else{

     sym_load_problem(env);
     
     sym_find_initial_bounds(env);
     
     sym_set_str_param(env, "lp_executable_name", "vrp_lp_cg");
     sym_set_str_param(env, "cp_executable_name", "vrp_cp");
     sym_set_int_param(env, "generate_cgl_cuts", FALSE);

     sym_solve(env);

   }

#if USE_VRPH
   final_stop=clock();

   /* Print timings to stderr: file vrph_time, sym_time */
   fprintf(stderr,"%s %5.3f %5.3f\n", vrp->par.infile,
	   (double)(vrph_stop-start)/CLOCKS_PER_SEC,
	   (double)(final_stop-vrph_stop)/CLOCKS_PER_SEC);
#endif

   sym_close_environment(env);
     
   return(0);
}

/*===========================================================================*/
/*===========================================================================*/

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

#endif
