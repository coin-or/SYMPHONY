/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*===========================================================================*\
 * This file contains the main() for the SYMPHONY generic MIP solver.
 * Note that, if you want to use the OSI SYMPHONY interface, you should set the
 * USE_OSI_INTERFACE flag and define the COINROOT path in the SYMPHONY 
 * Makefile. Otherwise, the C callable library functions will be used by 
 * default. See below for the usage.
\*===========================================================================*/


#ifdef USE_OSI_INTERFACE

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{

   /* Create an OsiSym object */
   OsiSymSolverInterface si;

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

#include "symphony_api.h"
  
int sym_help(char *key);
   
int main(int argc, char **argv)
{    
     
   sym_environment *env = sym_open_environment();
   
   if (argc > 1){
   
     sym_parse_command_line(env, argc, argv);

     if (env->par.test){

       sym_test(env);

     } else {
     
       sym_load_problem(env);
       
       sym_find_initial_bounds(env);
     
       sym_solve(env);
     }
   
   } else{
     FILE *f = NULL;
     char key[MAX_LINE_LENGTH +1], value[MAX_LINE_LENGTH +1];
     char infile1[MAX_LINE_LENGTH +1], infile2[MAX_LINE_LENGTH +1];     
     char ext[4], *line = NULL;     
     int last_dot = 0, j, terminate = FALSE, termcode = 0, int_value = 0;
     char * is_int = NULL;
     double *colsol = NULL, objval = 0.0, initial_time = 0.0, start_time = 0.0;
     double finish_time = 0.0, dbl_value = 0;

     strcpy(key, "");
     strcpy(value, "");
     strcpy(infile2, "");
     
     printf("***** WELCOME TO SYMPHONY INTERACTIVE MIP SOLVER ******\n\n"
	    "Please type 'help'/'?' to see the main commands!\n\n");

     env->par.verbosity = -1;

     while(true){
       printf("SYMPHONY: ");
       scanf("%s", &key);
       
       if (strcmp(key, "help") == 0 || strcmp(key, "?") == 0) {
	 sym_help("main_help");
       } else if (strcmp(key, "load") == 0){ 
	 if(env->mip->n){
	   free_master_u(env);
	   strcpy(env->par.infile, "");
	   strcpy(env->par.datafile, "");
	   env->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
	 }
	 printf("Name of the file: ");
	 strcpy(infile1, "");
	 scanf("%s", &infile1);       
	 
	 if (fopen(infile1, "r") == NULL){
	    printf("Input file '%s' can't be opened\n",
		  infile1);
	   continue;
	 }

	 /* check to see if SYMPHONY knows the input type! */

	 last_dot = 0;
	 for (j = 0;; j++){
	   if (infile1[j] == '\0')
	     break;
	   if (infile1[j] == '.') {
	     last_dot = j;
	   }
	 } 
	 strcpy(ext, infile1 + last_dot + 1);

	 if(!(strcmp(ext, "mod") == 0 || strcmp(ext, "mps") == 0)){
	   while(true){
	     printf("Type of the file ('mps'/'ampl'/'gmpl'): ");
	     scanf("%s", &ext);
	     if(!(strcmp(ext, "mps") == 0 || strcmp(ext, "ampl") == 0 ||
		  strcmp(ext, "gmpl") == 0)){
	       printf("Unknown type!\n");
	       continue; 
	     } else {
	       break;
	     }
	   }
	 }
	 
	 if (strcmp(ext, "mps") == 0){
	   if(sym_read_mps(env, infile1)){
	     continue;
	   }
	 } else {
	   printf("Name of the data file: ");
	   strcpy(infile2, "");
	   scanf("%s", &infile2); 
	   
	   if(strcmp(infile2, "") != 0 ) {
	     if(fopen(infile2, "r") == NULL){
	       printf("Data file '%s' can't be opened\n",
		      infile2);
	       continue;
	     }
	   }
	   if(sym_read_gmpl(env, infile1, infile2)){
	     continue;
	   }
	 }
       } else if(strcmp(key, "solve") == 0 || strcmp(key, "lpsolve") == 0){
	 if(!env->mip->n){
	   printf("No loaded problem. Use 'load' to read in a problem!\n");
	   continue;
	 } 
	 if(strcmp(key, "solve") == 0){	  
	   start_time = wall_clock(NULL);
	   termcode = sym_solve(env);
	   finish_time = wall_clock(NULL);
	 } else {
	   is_int = env->mip->is_int;
	   env->mip->is_int  = (char *)   calloc(CSIZE, env->mip->n);
	   start_time = wall_clock(NULL);
	   termcode = sym_solve(env);
	   finish_time = wall_clock(NULL);
	   env->mip->is_int = is_int;
	   is_int = 0;
	 }
       } else if (strcmp(key, "display") == 0){
	 printf("Please type 'help'/'?' to see the display options!\n");
	 while (true){
	   printf("SYMPHONY\\Display: ");	 
	   scanf("%s", &key);
	   if (strcmp(key, "help") == 0 || strcmp(key, "?") == 0) {
	     sym_help("display_help");
	   } else if (strcmp(key, "solution") == 0 || strcmp(key, "obj") == 0 
		      || strcmp(key, "stats") == 0){
	     if(!env->mip->n){
	       printf("No loaded problem! "
		      "Use 'load' in the main menu to read in a problem!\n");
	       continue;
	     } 
	     if(strcmp(key, "solution") == 0){
	       if(colsol) FREE(colsol);
	       colsol = (double *) malloc(DSIZE * env->mip->n);
	       if(sym_get_col_solution(env, colsol)){
		 printf("Error in displaying solution! The problem is either "
			"infeasible or has not been solved yet!\n");
		 continue;
	       } else {
		 if (env->mip->colname){ 
		   printf("+++++++++++++++++++++++++++++++++++++++++++++++\n");
		   printf("Column names and values in the solution\n");
		   printf("+++++++++++++++++++++++++++++++++++++++++++++++\n");
		   for(j = 0; j<env->mip->n; j++){
		     printf("%8s %10.3f\n", env->mip->colname[j],
			    colsol[j]);
		   }
		   printf("\n");
		 }else{
		   printf("+++++++++++++++++++++++++++++++++++++++++++++++\n");
		   printf("User indices and values in the solution\n");
		   printf("+++++++++++++++++++++++++++++++++++++++++++++++\n");
		   for(j = 0; j<env->mip->n; j++){
		     printf("%7d %10.3f\n", j, colsol[j]);
		   }
		   printf("\n");
		 }
	       }
	     } else if (strcmp(key, "obj") == 0){
	       if(sym_get_obj_val(env, &objval)){
		 printf("Error in displaying objective value!" 
			"The problem is either infeasible" 
			"or has not been solved yet!\n");
		 continue;
	       } else { 
		 printf("Objective Value: %f\n", objval);
	       }
	     } else {
	       initial_time  = env->comp_times.readtime;
	       initial_time += env->comp_times.ub_overhead + 
		 env->comp_times.ub_heurtime;
	       initial_time += env->comp_times.lb_overhead + 
		 env->comp_times.lb_heurtime;
	       
	       print_statistics(&(env->warm_start->comp_times), 
				&(env->warm_start->stat),
				env->ub, env->lb, initial_time, 
				start_time, finish_time,
				env->mip->obj_offset, env->mip->obj_sense,
				env->has_ub);
	       printf("\n");
	     }
	   } else if (strcmp(key, "parameter") == 0){
	     printf("Please type 'help'/'?' " 
		    "to see the list of available parameters!\n");
	     while(true){
	       printf("SYMPHONY\\Display\\Parameter: ");
	       scanf("%s", &key);
	       if (strcmp(key, "help") == 0 || strcmp(key, "?") == 0) {
		 sym_help("display_param_help");
	       } else if (strcmp(key, "back") == 0){
		 break;
	       } else if (strcmp(key, "quit") == 0){
		 terminate = TRUE;
		 break;
	       } else {
		 if (sym_get_int_param(env, key, &int_value) == 0){
		   printf("The value of %s: %i\n", key, int_value);
		 } else if ( sym_get_dbl_param(env, key, &dbl_value) == 0){
		   printf("The value of %s: %f\n", key, dbl_value);
		 }else {
		   printf("Unknown parameter/command!\n");
		   continue;
		 }
	       }
	     }
	   } else if (strcmp(key, "back") == 0){
	     break;
	   } else if (strcmp(key, "quit") == 0){
	     terminate = TRUE;
	     break;
	   } else {
	     printf("Unknown command!\n");
	     continue;	
	   }     
	   if(terminate) break;
	 }
       } else if (strcmp(key, "set") == 0){
	 printf("Please type 'help'/'?' to see the list of parameters!\n");
	 line = (char*) malloc(CSIZE*(MAX_LINE_LENGTH+1));	     
	 while (true){
	   printf("SYMPHONY\\Set: ");	 
	   scanf("%s", &key);

	   if (strcmp(key, "help") == 0 || strcmp(key, "?") == 0) {
	     sym_help("set_help");
	   } else if (strcmp(key, "back") == 0){
	     break;
	   } else if (strcmp(key, "quit") == 0){
	     terminate = TRUE;
	     break;
	   } else if (strcmp(key, "param_file") == 0){
	     printf("Name of the parameter file: ");
	     scanf("%s", infile1);
	     if ((f = fopen(infile1, "r")) == NULL){
	       printf("Parameter file '%s' can't be opened\n",
		      infile1);
	       continue;
	     }
	     while(NULL != fgets(line, MAX_LINE_LENGTH, f)){  /*read in parameters*/
	       strcpy(key, "");
	       strcpy(value, "");
	       sscanf(line,"%s%s", key, value);
	       if(set_param(env, line) == 0){
		 printf("Setting %s to: %s\n", key, value); 
	       } else {
		 printf("Unknown parameter %s: !\n", key);
		 continue;
	       }	     
	     }
	     fclose(f);
	   } else {
	     printf("Value of the parameter: ");
	     scanf("%s", value);
	     sprintf(line, "%s %s", key, value);  
	     if(set_param(env, line) == 0){
	       printf("Setting %s to: %s\n", key, value); 
	     } else {
	       printf("Unknown parameter/command!\n");
	       continue;
	     }
	   }
         }
         if(line) FREE(line);       
       } else if (strcmp(key, "reset") == 0){
	 printf("Resetting...\n");
	 sym_close_environment(env);
	 env = sym_open_environment();
	 env->par.verbosity = -1;
       } else if (strcmp(key, "quit") == 0){
	 break;
       } else {
	 printf("Unknown command!\n");
	 continue;
       }

       if(terminate) break;

     }
   } 
   
   sym_close_environment(env);
  
   return(0);
}    
     
int sym_help(char *key)
{    
  if(strcmp(key, "main_help") == 0){

    printf("\nList of main commands: \n\n");
    printf("load      : read a problem in mps or ampl format\n"
	   "solve     : solve the problem\n"
	   "lpsolve   : solve the lp relaxation of the problem\n"
	   "set       : set a parameter\n"
	   "display   : display optimization results and stats\n"
	   "reset     : restart the optimizer\n"
	   "help      : show the available commands/params/options\n\n"

	   "quit      : leave the optimizer\n\n");
    
  } else if (strcmp(key, "set_help") == 0 || strcmp(key, "display_param_help") == 0){

    printf("\n\nList of parameters: \n\n"); 
    printf("verbosity                          : set verbosity (default: 1)\n"
	   "upper_bound                        : use an initial upper bound\n"
	   "find_first_feasible                : whether to find the first feasible solution or\n"
	   "                                     to solve the optimality (default: 0) \n"
	   "generate_cgl_cuts                  : whether or not to use cgl cuts (default: 1)\n"
	   "generate_cgl_gomory_cuts           : whether or not to use cgl gomory cuts (default: 1)\n"
	   "generate_cgl_knapsack_cuts         : whether or not to use cgl knapsack cuts (default: 1)\n"
	   "generate_cgl_oddhole_cuts          : whether or not to use cgl oddhole cuts (default: 1)\n"
	   "generate_cgl_probing_cuts          : whether or not to use cgl probing cuts (default: 1)\n"
	   "generate_cgl_clique_cuts           : whether or not to use cgl clique cuts (default: 1)\n"	   
	   "generate_cgl_mir_cuts              : whether or not to use cgl mixed integer rounding cuts\n" 
           "                                     (default: 0)\n"
 	   "generate_cgl_flow_and_cover_cuts   : whether or not to use cgl flow and cover cuts (default: 1)\n"
	   "generate_cgl_rounding_cuts         : whether or not to use cgl rounding cuts (default: 0)\n"
	   "generate_cgl_lift_and_project_cuts : whether or not to use cgl lift and project cuts (default: 0)\n"
	   "node_selection_rule                : set the node selection rule/search strategy (default: 5)\n"
	   "strong_branching_candidate_num     : set the stong branching candidates number (default: var)\n"
	   "compare_candidadates_dafult        : set the rule to compare the candidates (defualt: 2)\n"
	   "select_child_default               : set the rule to select the children (default: 0)\n"
	   "diving_threshold                   : set diving threshold (default: 0)\n"
	   "diving_strategy                    : set diving strategy (default: 0)\n"
	   "do_reduced_cost_fixing             : whether ot not to use reduced cost fixing (default: 1)\n"
	   "time_limit                         : set the time limit\n"
	   "node_limit                         : set the node limit\n"
	   "gap_limit                          : set the target gap between the lower and upper bound\n"
           "param_file                         : read parameters from a parameter file\n\n"

	   "back                               : leave this menu\n"
	   "quit                               : leave the optimizer\n\n");
					    
  } else if (strcmp(key, "display_help") == 0){

    printf("\nList of display options: \n\n");
    printf("solution     : display the column values\n"
	   "obj          : display the objective value\n"
	   "stats        : display the statistics\n"
	   "parameter    : display the value of a parameter\n\n"

	   "back         : leave this menu\n"
	   "quit         : leave the optimizer\n\n");
  }

  return(0);
}

#endif

