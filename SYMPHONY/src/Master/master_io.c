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

#include <malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "timemeas.h"
#include "BB_constants.h"
#include "BB_macros.h"
#include "BB_types.h"
#include "lp_params.h"
#include "master.h"
#include "master_u.h"

void usage(void);

/*===========================================================================*/

/*===========================================================================*\
 * This file contains I/O functions for the master process.
\*===========================================================================*/

/*===========================================================================*/

void usage(void)
{
   printf("Generic switches:\n\n");
   printf("master [ -hagrtbd ] [ -u ub ] [ -p procs ] [ -n rule ]\n\t"
	  "[ -v level ] [ -s cands ] [ -c rule ] [ -k rule ] \n\t"
	  "[ -m max ] [ -l pools ] [ -i iters ] "
	  "[ -f parameter_file_name ] [-j 0/1]"
	  "\n\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n"
	  "\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n"
	  "\t%s\n\t%s\n\t%s\n\t%s\n\n",
	  "-h: help",
	  "-a: no cut timeout",
	  "-d: enable graph drawing",
	  "-g: use cut generator",
	  "-r: do repricing in root",
	  "-t: trim the tree",
	  "-b: don't perform branch and cut",
	  "-u ub: use upper bound 'ub'",
	  "-p procs: allow 'procs' active nodes",
	  "-n i: use node selection rule 'i'",
	  "-v i: set verbosity to level 'i'",
	  "-s cands: use 'cands' candidates for strong branching",
	  "-c i: use rule 'i' to compare candidates",
	  "-k i: use rule 'i' to select child",
	  "-m n: allow a max of 'n' cuts to enter per iteration",
	  "-e n: allow a max of 'n' cut pools",
	  "-l n k: load balance level 'n' and iterations 'k'",
	  "-i n: allow a max of 'n' iterations in presolve",
	  "-f file: read parameters from parameter file 'file'",
	  "-j 0/1: whether or not to generate cgl cuts",
	  "-z n: set diving threshold to 'n'");
   printf("Solver-specific switches:\n\n");
#ifdef USE_SYM_APPLICATION
   user_usage();
#else
   printf("master [ -H ] [ -F file ] \n\n\t%s\n\t%s\n\t%s\n\t%s\n\n",
	  "-H: help (solver-specific switches)",
	  "-F model: model should be read in from file 'model'",
	  "          (MPS format is assumed unless -D is also present)",
	  "-D data: model is in AMPL format and data is in file 'data'");
#endif   
}

/*===========================================================================*/
/*===========================================================================*/

int parse_command_line(sym_environment *env, int argc, char **argv)
{
   int i, tmpi;
   double tmpd;
   char line[MAX_LINE_LENGTH +1], tmp, c;
   char key[MAX_LINE_LENGTH +1], value[MAX_LINE_LENGTH +1];
   FILE *f = NULL, *f1 = NULL;
   //   str_int colgen_str[COLGEN_STR_SIZE] = COLGEN_STR_ARRAY;
   tm_params *tm_par = &env->par.tm_par;
   lp_params *lp_par = &env->par.lp_par;
   cg_params *cg_par = &env->par.cg_par;
   cp_params *cp_par = &env->par.cp_par;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   sp_params *sp_par = &env->par.sp_par;
#endif
   /*___END_EXPERIMENTAL_SECTION___*/
   dg_params *dg_par = &env->par.dg_par;

   if (argc < 2){
      usage();
      exit(0);
   }
   
   printf("SYMPHONY was called with the following arguments:\n");
   printf("%s ", argv[0]);
   for (i = 1; i < argc; i++){
      sscanf(argv[i], "%c", &tmp);
      if (tmp == '-')
	 printf("\n");
      printf("%s ", argv[i]);
   }
   printf("\n\n");
   
   for (i = 0; i < argc; i++){
      if (!strcmp(argv[i], "-f"))
	 break;
   }
   
   if (i == argc){
      goto EXIT;
   }else{
      strncpy(env->par.param_file, argv[i+1], MAX_FILE_NAME_LENGTH);
   }
   
   if ((f = fopen(env->par.param_file, "r")) == NULL){
      (void) fprintf(stderr, "Readparams: file '%s' can't be opened\n\n",
		     env->par.param_file);
      return(ERROR__OPENING_PARAM_FILE);
   }

   printf("============= Other Parameter Settings =============\n\n");

   while (NULL != fgets(line, MAX_LINE_LENGTH, f)){  /* read in parameters */

      set_param(env, line);

      printf("%s", line);
      strcpy(key,"");
      sscanf(line,"%s%s", key, value);

      if (strcmp(key, "lp_mach_num") == 0 ||
	  strcmp(key, "TM_lp_mach_num") == 0){
	 if (tm_par->lp_mach_num){
	    char *lp_machs = (char *) malloc
	       (tm_par->lp_mach_num * (MACH_NAME_LENGTH + 1));
	    tm_par->lp_machs =
	       (char **) malloc(tm_par->lp_mach_num * sizeof(char *));
	    for (i=0; i<tm_par->lp_mach_num; i++)
	       tm_par->lp_machs[i] = lp_machs + i * (MACH_NAME_LENGTH+1);
	    for (i=0; i<tm_par->lp_mach_num; i++){
	       if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
		  fprintf(stderr, "\nio: error reading lp_machine list\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "TM_lp_machine") != 0){
		  fprintf(stderr, "\nio: error reading lp_machine list\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       read_string(tm_par->lp_machs[i], line, MACH_NAME_LENGTH);
	       printf("%s", line);
	    }
	 }
      }
      else if (strcmp(key, "cg_mach_num") == 0 ||
	       strcmp(key, "TM_cg_mach_num") == 0){
	 if (tm_par->cg_mach_num){
	    char *cg_machs = (char *) malloc
	       (tm_par->cg_mach_num * (MACH_NAME_LENGTH + 1));
	    tm_par->cg_machs =
	       (char **) malloc(tm_par->cg_mach_num * sizeof(char *));
	    for (i=0; i<tm_par->cg_mach_num; i++)
	       tm_par->cg_machs[i] = cg_machs + i * (MACH_NAME_LENGTH+1);
	    for (i=0; i<tm_par->cg_mach_num; i++){
	       if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
		  fprintf(stderr, "\nio: error reading cg_machine list\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "TM_cg_machine") != 0){
		  fprintf(stderr, "\nio: error reading cg_machine list\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       read_string(tm_par->cg_machs[i], line, MACH_NAME_LENGTH);
	       printf("%s", line);
	    }
	 }
      }
      else if (strcmp(key, "cp_mach_num") == 0 ||
	       strcmp(key, "TM_cp_mach_num") == 0){
	 if (tm_par->cp_mach_num){
	    char *cp_machs = (char *) malloc
	       (tm_par->cp_mach_num * (MACH_NAME_LENGTH + 1));
	    tm_par->cp_machs =
	       (char **) malloc(tm_par->cp_mach_num * sizeof(char *));
	    for (i=0; i<tm_par->cp_mach_num; i++)
	       tm_par->cp_machs[i] = cp_machs + i * (MACH_NAME_LENGTH+1);
	    for (i=0; i<tm_par->cp_mach_num; i++){
	       if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
		  fprintf(stderr, "\nio: error reading cp_machine list\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "TM_cp_machine") != 0){
		  fprintf(stderr, "\nio: error reading cp_machine list\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       read_string(tm_par->cp_machs[i], line, MACH_NAME_LENGTH);
	       printf("%s", line);
	    }
	 }
      }
      else if (strcmp(key, "keep_description_of_pruned") == 0 ||
	       strcmp(key, "TM_keep_description_of_pruned") == 0){
	 if (tm_par->keep_description_of_pruned == KEEP_ON_DISK_FULL ||
	     tm_par->keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No pruned node file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "pruned_node_file_name") != 0){
	       printf("Need pruned_node_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(tm_par->pruned_node_file_name, value);
	    if (!(f1 = fopen(tm_par->pruned_node_file_name, "w"))){
	       printf("\nError opening pruned node file\n\n");
	    }else{
	       if (tm_par->keep_description_of_pruned == KEEP_ON_DISK_FULL){
		  fprintf(f1, "******* Pruned Node Log File *******\n\n");
	       }else{
		  fprintf(f1, "#TYPE: COMPLETE TREE\n");
		  fprintf(f1, "#TIME: NOT\n");
		  fprintf(f1, "#BOUNDS: NONE\n");
		  fprintf(f1, "#INFORMATION: EXCEPTION\n");
		  fprintf(f1, "#NODE_NUMBER: NONE\n");
	       }
	       fclose(f1);
	    }
	 }
      }
      else if (strcmp(key, "warm_start") == 0 ||
	       strcmp(key, "TM_warm_start") == 0){
	 if ((env->par.warm_start = tm_par->warm_start)){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No warm start tree file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "warm_start_tree_file_name") != 0){
	       printf("Need warm_start_tree_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(tm_par->warm_start_tree_file_name, value);
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No warm start cut file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "warm_start_cut_file_name") != 0){
	       printf("Need warm_start_cut_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(tm_par->warm_start_cut_file_name, value);
	 }
      }
      else if (strcmp(key, "vbc_emulation") == 0 ||
	       strcmp(key, "TM_vbc_emulation") == 0){
	 if (tm_par->vbc_emulation == VBC_EMULATION_FILE){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No vbc emulation file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "vbc_emulation_file_name") != 0){
	       printf("Need vbc_emulation_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(tm_par->vbc_emulation_file_name, value);
	    if (!(f1 = fopen(tm_par->vbc_emulation_file_name, "w"))){
	       printf("\nError opening vbc emulation file\n\n");
	    }else{
	       fprintf(f1, "#TYPE: COMPLETE TREE\n");
	       fprintf(f1, "#TIME: SET\n");
	       fprintf(f1, "#BOUNDS: NONE\n");
	       fprintf(f1, "#INFORMATION: STANDARD\n");
	       fprintf(f1, "#NODE_NUMBER: NONE\n");
	       fprintf(f1, "00:00:00.00 N 0 1 %i\n", VBC_CAND_NODE);
	       fclose(f1);
	    }
	 }else if (tm_par->vbc_emulation == VBC_EMULATION_LIVE){
	    printf("$#TYPE: COMPLETE TREE\n");
	    printf("$#TIME: SET\n");
	    printf("$#BOUNDS: NONE\n");
	    printf("$#INFORMATION: STANDARD\n");
	    printf("$#NODE_NUMBER: NONE\n");
	    printf("$N 0 1 %i\n", VBC_CAND_NODE);
	 }
      }
      else if (strcmp(key, "logging") == 0 ||
	       strcmp(key, "TM_logging") == 0){
	 if (tm_par->logging){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No tree log file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "tree_log_file_name") != 0){
	       printf("tree_log_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(tm_par->tree_log_file_name, value);
	    if (tm_par->logging != VBC_TOOL){
	       if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
		  printf("No cut log file!\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "cut_log_file_name") != 0){
		  printf("Need cut_log_file_name next!!!\n\n");
		  return(ERROR__PARSING_PARAM_FILE);
	       }
	       strcpy(tm_par->cut_log_file_name, value);
	    }
	 }
      }
      else if (strcmp(key, "cp_warm_start") == 0 ||
	       strcmp(key, "CP_warm_start") == 0){
	 if (cp_par->warm_start){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No cut pool warm start file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "cp_warm_start_file_name") != 0){
	       printf("Need cp_warm_start_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(cp_par->warm_start_file_name, value);
	 }
      }
      else if (strcmp(key, "cp_logging") == 0 ||
	       strcmp(key, "CP_logging") == 0){
	 if ((tm_par->cp_logging = cp_par->logging)){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No cut pool log file!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "cp_log_file_name") != 0){
	       printf("Need cp_log_file_name next!!!\n\n");
	       return(ERROR__PARSING_PARAM_FILE);
	    }
	    strcpy(cp_par->log_file_name, value);
	 }
      }
   }

   printf("\n====================================================\n\n");
   
 EXIT:
   
   for (i = 1; i < argc; i++){
      sscanf(argv[i], "%c %c", &tmp, &c);
      if (tmp != '-')
	 continue;
      switch (c) {
       case 'h':
	 usage();
	 exit(0);
       case 'H':
#ifdef USE_SYM_APPLICATION
	  user_usage();
#else
	  printf("master [ -H ] [ -F file ] \n\n\t%s\n\t%s\n\t%s\n\t%s\n\n",
		 "-H: help (solver-specific switches)",
		 "-F model: model should be read in from file 'model'",
		 "          (MPS format is assumed unless -D is also present)",
		 "-D data: model is in AMPL format and data is in file 'data'");
#endif 
	  exit(0);
       case 'a':
	 lp_par->first_lp.first_cut_time_out = 0;
	 lp_par->first_lp.all_cuts_time_out = 0;
	 lp_par->later_lp.first_cut_time_out = 0;
	 lp_par->later_lp.all_cuts_time_out = 0;
	 /*__BEGIN_EXPERIMENTAL_SECTION__*/
	 cg_par->decomp_dynamic_timeout = 6000;
	 /*___END_EXPERIMENTAL_SECTION___*/
	 break;
       case 'd':
	 env->par.do_draw_graph = TRUE;
	 break;
       case 'g':
	 lp_par->use_cg = tm_par->use_cg = TRUE;
	 break;
       case 'r':
	 tm_par->price_in_root = TRUE;
	 break;
       case 't':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%lf", &tmpd)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->time_limit = tm_par->time_limit = tmpd;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'b':
	 env->par.do_branch_and_cut = FALSE;
	 break;
       case 'u':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%lf", &tmpd)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       env->ub = tmpd;
	       env->has_ub = TRUE;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'p':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       tm_par->max_active_nodes = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'n':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       tm_par->node_selection_rule = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'v':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       tm_par->verbosity = lp_par->verbosity = cg_par->verbosity =
		  /*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
		  sp_par->verbosity =
#endif 
		  /*___END_EXPERIMENTAL_SECTION___*/
		  cp_par->verbosity = env->par.verbosity = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
 	 break;
       case 's':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->strong_branching_cand_num_min =
	       lp_par->strong_branching_cand_num_max = tmpi;
	       lp_par->strong_branching_red_ratio = 0;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'c':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->compare_candidates_default = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'k':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->select_child_default = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'm':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->max_cut_num_per_iter = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'e':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       tm_par->max_cp_num = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'l':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->load_balance_level = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->load_balance_iterations = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'i':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->max_presolve_iter = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'f':
	 if (i < argc - 1){
	    sscanf(argv[i+1], "%c", &tmp);
	    if (tmp == '-'){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       strncpy(env->par.param_file, argv[i+1], MAX_FILE_NAME_LENGTH);
	       i++;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'j':
       if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%i", &tmpi)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       lp_par->cgl.generate_cgl_cuts = tmpi;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       case 'z':
	 if (i < argc - 1){
	    if (!sscanf(argv[i+1], "%lf", &tmpd)){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	    }else{
	       i++;
	       tm_par->diving_threshold = tmpd;
	    }
	 }else{
	    printf("Warning: Missing argument to command-line switch -%c\n",c);
	 }
	 break;
       default:
	 if (c < 'A'){
	    printf("Warning: Ignoring unrecognized command-line switch -%c\n",
		   c);
	 }
	 break;
      };
   }

   /*Sanity checks*/

   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   if (cg_par->decomp_max_col_num_per_iter >
       cg_par->decomp_col_block_size){
      printf("io: decomp_max_col_num_per_iter is greater than\n");
      printf("    decomp_col_block_size -- adjusting\n");
      cg_par->decomp_max_col_num_per_iter = cg_par->decomp_col_block_size;
   }
	 
   /*___END_EXPERIMENTAL_SECTION___*/
   if (cp_par->block_size >cp_par->max_number_of_cuts){
      printf("io: Cut pool block size is too big -- adjusting\n");
      cp_par->block_size = cp_par->max_number_of_cuts;
   }

   if (cp_par->min_to_delete > cp_par->max_number_of_cuts -
                               cp_par->cuts_to_check){
      printf("io: Cut pool min to delete is too big -- adjusting\n");
      cp_par->min_to_delete = cp_par->max_number_of_cuts -
	                      cp_par->cuts_to_check;
   }

   /*if (tm_par->price_in_root &&
       tm_par->colgen_strat[0] != (FATHOM__DO_NOT_GENERATE_COLS__SEND |
				   BEFORE_BRANCH__DO_NOT_GENERATE_COLS)){
      printf("io: pricing in root is asked for but colums are to be\n");
      printf("    generated in the 1st phase -- adjusting colgen_strat[0]\n");
      tm_par->colgen_strat[0] = (FATHOM__DO_NOT_GENERATE_COLS__SEND |
				 BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   }*/

   if (f)
      fclose(f);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

void read_string(char *target, char *line, int maxlen)
{
   char key[MAX_LINE_LENGTH +1], value[MAX_LINE_LENGTH +1], *quote1, *quote2;
   int len;

   if (sscanf(line, "%s%s", key, value) != 2)
      READPAR_ERROR(key);

   if (value[0] != '"'){ /* the string is not quoted */
      quote1 = value;
      len = strlen(quote1);
   }else{ /* the string is quoted */
      quote1 = strchr(line, '"');
      quote2 = strrchr(line,'"');
      if (quote1 == quote2)
	 READPAR_ERROR(key);
      quote1++;
      len = quote2 - quote1;
   }
   
   if (len > maxlen)
      READPAR_ERROR(key);
   if (len > 0)
      strncpy(target, quote1, len);
   target[len] = 0;
   if (strchr(target, '{') || strchr(target, '}'))
      READPAR_ERROR(key);
}

/*===========================================================================*/
/*===========================================================================*/

void print_statistics(node_times *tim, problem_stat *stat, double ub,
		      double lb, double initial_time, double start_time,
		      double finish_time, double obj_offset, char obj_sense, 
		      char has_ub)
{
   double gap = 0.0;

   static str_int nfstatus[4] = {
      {"NF_CHECK_ALL"           , NF_CHECK_ALL }
      , {"NF_CHECK_AFTER_LAST"    , NF_CHECK_AFTER_LAST }
      , {"NF_CHECK_UNTIL_LAST"    , NF_CHECK_UNTIL_LAST }
      , {"NF_CHECK_NOTHING"       , NF_CHECK_NOTHING }
   };

   initial_time += tim->communication;
   initial_time += tim->lp;
   initial_time += tim->separation;
   initial_time += tim->fixing;
   initial_time += tim->pricing;
   initial_time += tim->strong_branching;
   initial_time += tim->cut_pool;
#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
   printf("======================= CP Timing ===========================\n");
   printf("  Cut Pool                  %.3f\n", tim->cut_pool);
#endif
   printf("====================== LP/CG Timing =========================\n");
#ifndef WIN32  /* FIXME: CPU timing doesn't work in Windows */
   printf("  LP Solution Time          %.3f\n", tim->lp);
   printf("  Variable Fixing           %.3f\n", tim->fixing);
   printf("  Pricing                   %.3f\n", tim->pricing);
   printf("  Strong Branching          %.3f\n", tim->strong_branching);
   printf("  Separation                %.3f\n", tim->separation); 
#ifndef COMPILE_IN_LP
   printf("=================== Parallel Overhead ======================\n");
   printf("  Communication         %.3f\n", tim->communication);
   printf("  Ramp Up Time (TM)     %.3f\n", tim->ramp_up_tm);
   printf("  Ramp Up Time (LP)     %.3f\n", tim->ramp_up_lp);
   printf("  Ramp Down Time        %.3f\n", tim->ramp_down_time);
   printf("  Idle Time (Node Pack) %.3f\n", tim->start_node);
   printf("  Idle Time (Nodes)     %.3f\n", tim->idle_node);
   printf("  Idle Time (Names)     %.3f\n", tim->idle_names);
   printf("  Idle Time (Diving)    %.3f\n", tim->idle_diving);
   printf("  Idle Time (Cuts)      %.3f\n", tim->idle_cuts);
#endif
   printf("  Total User Time              %.3f\n", initial_time);
#endif
   printf("  Total Wallclock Time         %.3f\n\n", finish_time -
	  start_time);
   printf("====================== Statistics =========================\n");
   printf("Number of created nodes :       %i\n", stat->created);
   printf("Number of analyzed nodes:       %i\n", stat->analyzed);
   printf("Depth of tree:                  %i\n", stat->max_depth);
   printf("Size of the tree:               %i\n", stat->tree_size);
#if 0
   printf("Leaves before trimming:         %i\n",
	  stat->leaves_before_trimming);
   printf("Leaves after trimming:          %i\n", stat->leaves_after_trimming);
   printf("Repriced root's nf_status:      %s\n",
	  nfstatus[(int)stat->nf_status].str);
   printf("Not fixed variable num:         %i\n", stat->vars_not_priced);
#endif
   printf("Number of Chains:               %i\n", stat->chains);
   printf("Number of Diving Halts:         %i\n", stat->diving_halts);
   printf("Number of cuts in cut pool:     %i\n", stat->cuts_in_pool);
   if (stat->root_lb > -MAXDOUBLE){
      if (obj_sense == SYM_MAXIMIZE){
	 printf("Upper Bound in Root:            %.3f\n",
		-stat->root_lb + obj_offset);
      }else{
	 printf("Lower Bound in Root:            %.3f\n",
		stat->root_lb + obj_offset);
      }
   }

   if (has_ub){
     gap = abs(100*(ub-lb)/ub);
   }

   if (obj_sense == SYM_MAXIMIZE){
     if (gap > -1e-07 && gap < 0){
       printf("\nCurrent Lower Bound:         %.3f", -ub + obj_offset);
       printf("\nCurrent Upper Bound:         %.3f", -lb + obj_offset);
       printf("\nGap Percentage:              %.2f\n", -gap);
     } else if (!has_ub) {
       printf("\nCurrent Upper Bound:         %.3f\n", -lb + obj_offset);
     }
   }else{
     if (gap > 1e-07){
       printf("\nCurrent Upper Bound:         %.3f", ub + obj_offset);
       printf("\nCurrent Lower Bound:         %.3f", lb + obj_offset);
       printf("\nGap Percentage:              %.2f\n", gap);
     } else if (!has_ub){
       printf("\nCurrent Lower Bound:         %.3f\n", lb + obj_offset);
     }
   }
}
