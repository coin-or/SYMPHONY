/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2002 Ted Ralphs. All Rights Reserved.           */
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
#include "lp_params.h"
#include "master.h"

void usage(void);

/*===========================================================================*/

/*===========================================================================*\
 * This file contains I/O functions for the master process.
\*===========================================================================*/

/*===========================================================================*/

void usage(void)
{
         printf("master [ -hagrtbd ] [ -u ub ] [ -p procs ] [ -n rule ]\n\t"
		"[ -v level ] [ -s cands ] [ -c rule ] [ -k rule ] \n\t"
		"[ -m max ] [ -l pools ] [ -i iters ] "
		"[ -f parameter_file_name ]"
		"\n\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n"
		"\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n"
		"\t%s\n\t%s\n\t%s\n\n",
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
   		"-z n: set diving threshold to 'n'",
	 "-f file: read parameters from parameter file 'file'");
	 printf("Type 'master -H' to get help for user options\n\n");
}

/*===========================================================================*/

void bc_readparams(problem *p, int argc, char **argv)
{
   int i;
   char line[MAX_LINE_LENGTH +1], tmp, c;
   char key[MAX_LINE_LENGTH +1], value[MAX_LINE_LENGTH +1];
   FILE *f = NULL, *f1 = NULL;
   double timeout;
   str_int colgen_str[COLGEN_STR_SIZE] = COLGEN_STR_ARRAY;
   str_int compare_can_str[COMPARE_CAN_STR_SIZE] = COMPARE_CAN_STR_ARRAY;
   
   tm_params *tm_par = &p->par.tm_par;
   lp_params *lp_par = &p->par.lp_par;
   cg_params *cg_par = &p->par.cg_par;
   cp_params *cp_par = &p->par.cp_par;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
   sp_params *sp_par = &p->par.sp_par;
#endif
   /*___END_EXPERIMENTAL_SECTION___*/
   dg_params *dg_par = &p->par.dg_par;

   /************************* Global defaults ********************************/
   p->ub = 0;
   p->lb = 0;
   p->par.verbosity = 0;
   p->par.random_seed = 17;
   p->par.tm_machine_set = FALSE;
   p->par.dg_machine_set = FALSE;
   strcpy(p->par.tm_exe, "tm");
#ifdef COMPILE_IN_LP
   strcat(p->par.tm_exe, "_lp");
#ifdef COMPILE_IN_CG
   strcat(p->par.tm_exe, "_cg");
#endif
#endif
#ifdef COMPILE_IN_CP
   strcat(p->par.tm_exe, "_cp");
#endif   
   strcpy(p->par.dg_exe, "dg");
   p->par.tm_debug = 0;
   p->par.dg_debug = 0;
   p->par.pvm_trace = 0;
   p->par.do_branch_and_cut = 1;
   p->par.do_draw_graph = FALSE;

   /************************** treemanager defaults **************************/
   tm_par->verbosity = 0;
   tm_par->granularity = 0.000001;
   strcpy(tm_par->lp_exe, "lp");
#ifdef COMPILE_IN_CG
   strcat(tm_par->lp_exe, "_cg");
#endif
   strcpy(tm_par->cg_exe, "cg");
   strcpy(tm_par->cp_exe, "cp");
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   strcpy(tm_par->sp_exe, "sp");
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->lp_debug = 0;
   tm_par->cg_debug = 0;
   tm_par->cp_debug = 0;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   tm_par->sp_debug = 0;
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->max_active_nodes = 1;
   tm_par->max_cp_num = 1;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   tm_par->max_sp_num = 0;
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->lp_mach_num = 0;
   tm_par->lp_machs = NULL;
   tm_par->cg_mach_num = 0;
   tm_par->cg_machs = NULL;
   tm_par->cp_mach_num = 0;
   tm_par->cp_machs = NULL;

   tm_par->use_cg = FALSE;
   tm_par->random_seed = 17;
   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   tm_par->do_decomp = FALSE;
   /*___END_EXPERIMENTAL_SECTION___*/
   tm_par->unconditional_dive_frac = .1;
   tm_par->diving_strategy = BEST_ESTIMATE;
   tm_par->diving_k = 1;
   tm_par->diving_threshold = 0;
   tm_par->node_selection_rule = LOWEST_LP_FIRST;

   tm_par->keep_description_of_pruned = DISCARD;
   tm_par->warm_start = FALSE;
   tm_par->logging = NO_LOGGING;
   tm_par->logging_interval = 1800;
   tm_par->vbc_emulation = NO_VBC_EMULATION;
   tm_par->price_in_root = FALSE;
   tm_par->trim_search_tree = FALSE;
   tm_par->colgen_strat[0] = (FATHOM__DO_NOT_GENERATE_COLS__DISCARD  |
			      BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   tm_par->colgen_strat[1] = (FATHOM__DO_NOT_GENERATE_COLS__DISCARD  |
			      BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   tm_par->not_fixed_storage_size = 2048;
   tm_par->time_limit = 0;

   /************************** lp defaults ***********************************/
   lp_par->verbosity = 0;
   lp_par->granularity = tm_par->granularity;
   lp_par->use_cg = tm_par->use_cg;
   lp_par->set_obj_upper_lim = FALSE;
   lp_par->scaling = -1; /* CPLEX'ism ... don't scale */
   lp_par->fastmip = 1; /* CPLEX'ism ... set it to 1 */
   lp_par->try_to_recover_from_error = TRUE;
   lp_par->problem_type = ZERO_ONE_PROBLEM;
   lp_par->keep_description_of_pruned = tm_par->keep_description_of_pruned;
   lp_par->not_fixed_storage_size = tm_par->not_fixed_storage_size;
   lp_par->cut_pool_check_freq = 10;
   lp_par->load_balance_level = -1;
   lp_par->load_balance_iterations = -1;
   lp_par->load_balance_compare_candidates = HIGHEST_LOW_OBJ;
   lp_par->fractional_diving_ratio = 0.02;
   lp_par->fractional_diving_num = 0;
   lp_par->max_non_dual_feas_to_add_frac = 0.05;
   lp_par->max_non_dual_feas_to_add_min = 20;
   lp_par->max_non_dual_feas_to_add_max = 200;
   lp_par->max_not_fixable_to_add_frac = 0.1;
   lp_par->max_not_fixable_to_add_min = 100;
   lp_par->max_not_fixable_to_add_max = 500;
   lp_par->mat_col_compress_num = 50;
   lp_par->mat_col_compress_ratio = .05;
   lp_par->mat_row_compress_num = 20;
   lp_par->mat_row_compress_ratio = .05;
   lp_par->tailoff_gap_backsteps = 2;
   lp_par->tailoff_gap_frac = .99;
   lp_par->tailoff_obj_backsteps = 3;
   lp_par->tailoff_obj_frac = .75;
   lp_par->ineff_cnt_to_delete = 0;
   lp_par->eff_cnt_before_cutpool = 3;
   lp_par->ineffective_constraints = BASIC_SLACKS_ARE_INEFFECTIVE;
   lp_par->base_constraints_always_effective = TRUE;
   lp_par->branch_on_cuts = FALSE;
   lp_par->discard_slack_cuts = DISCARD_SLACKS_BEFORE_NEW_ITERATION;
   lp_par->first_lp.first_cut_time_out = 0;
   lp_par->first_lp.all_cuts_time_out = 0;
   lp_par->later_lp.first_cut_time_out = 5;
   lp_par->later_lp.first_cut_time_out = 0;
   lp_par->later_lp.all_cuts_time_out = 1;
   lp_par->later_lp.all_cuts_time_out = 0;
   lp_par->max_cut_num_per_iter = 20;
   lp_par->do_reduced_cost_fixing = TRUE;
   lp_par->gap_as_ub_frac = .1;
   lp_par->gap_as_last_gap_frac = .7;
   lp_par->do_logical_fixing = 1;
   lp_par->fixed_to_ub_before_logical_fixing = 1;
   lp_par->fixed_to_ub_frac_before_logical_fixing = .01;

   lp_par->max_presolve_iter = 50;

   lp_par->is_feasible_default = TEST_INTEGRALITY;
   lp_par->send_feasible_solution_default = SEND_NONZEROS;
   lp_par->display_solution_default = DISP_NOTHING;
   lp_par->shall_we_branch_default = USER__BRANCH_IF_MUST;
   lp_par->select_candidates_default = USER__CLOSE_TO_HALF;
   lp_par->strong_branching_cand_num_max = 25;
   lp_par->strong_branching_cand_num_min = 5;
   lp_par->strong_branching_red_ratio = 1;
   lp_par->compare_candidates_default = HIGHEST_LOW_OBJ;
   lp_par->select_child_default = PREFER_LOWER_OBJ_VALUE;
   lp_par->pack_lp_solution_default = SEND_NONZEROS;

   /************************** cut_gen defaults *****************************/
   cg_par->verbosity = 0;
   cg_par->do_findcuts = TRUE;

   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   cg_par->do_decomp = FALSE;
   cg_par->decomp_sol_pool_check_freq = 10;
   cg_par->decomp_wait_for_cols = TRUE;
   cg_par->decomp_max_col_num_per_iter = 1000;
   cg_par->decomp_col_block_size = 1000;
   cg_par->decomp_mat_block_size = 100000;
   cg_par->decomp_initial_timeout = 5;
   cg_par->decomp_dynamic_timeout = 5;
   cg_par->decomp_complete_enum = TRUE;

   /*___END_EXPERIMENTAL_SECTION___*/
   /************************** cutpool defaults ******************************/
   cp_par->verbosity = 0;
   cp_par->warm_start = FALSE;
   cp_par->logging = FALSE;
   cp_par->block_size = 5000;
   cp_par->max_size = 2000000;
   cp_par->max_number_of_cuts = 10000;
   cp_par->cuts_to_check = 1000;
   cp_par->delete_which = DELETE_BY_QUALITY;
   cp_par->touches_until_deletion = 10;
   cp_par->min_to_delete = 1000;
   cp_par->check_which = CHECK_ALL_CUTS;

   /*__BEGIN_EXPERIMENTAL_SECTION__*/
   /************************** solpool defaults ******************************/
#ifdef COMPILE_DECOMP
   sp_par->verbosity = 0;
   sp_par->etol = 0.000001;
   sp_par->block_size = 1000;
   sp_par->max_size = 1000000;
   sp_par->max_number_of_sols = 10000;
   sp_par->delete_which = DELETE_DUPLICATE_COLS;
   sp_par->touches_until_deletion = 10;
   sp_par->min_to_delete = 100;
   sp_par->compress_num = 10;
   sp_par->compress_ratio = .01;
   sp_par->check_which = CHECK_COL_LEVEL_AND_TOUCHES;
#endif
   /*___END_EXPERIMENTAL_SECTION___*/
   /********************** draw_graph defaults  ******************************/
   strcpy(dg_par->source_path, ".");
   dg_par->echo_commands = FALSE;
   dg_par->canvas_width = 1000;
   dg_par->canvas_height = 700;
   dg_par->viewable_width = 600;
   dg_par->viewable_height = 400;
   dg_par->disp_nodelabels = 1;
   dg_par->disp_nodeweights = 1;
   dg_par->disp_edgeweights = 1;
   dg_par->node_dash[0] = 0;
   dg_par->edge_dash[0] = 0;
   dg_par->node_radius = 8;
   dg_par->interactive_mode = 1;
   dg_par->mouse_tracking = 1;
   dg_par->scale_factor = 1;
   strcpy(dg_par->nodelabel_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");
   strcpy(dg_par->nodeweight_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");
   strcpy(dg_par->edgeweight_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");

   if (argc < 2){
      usage();
      exit(1);
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
      strncpy(p->par.param_file, argv[i+1], MAX_FILE_NAME_LENGTH);
   }
   
   if ((f = fopen(p->par.param_file, "r")) == NULL){
      (void) fprintf(stderr, "Readparams: file '%s' can't be opened\n\n",
		     p->par.param_file);
      exit(1);
   }

   printf("============= Other Parameter Settings =============\n\n");

   while (NULL != fgets(line, MAX_LINE_LENGTH, f)){  /* read in parameters */
      printf("%s", line);
      strcpy(key,"");
      sscanf(line,"%s%s", key, value);

      /***********************************************************************
       ***                    Global parameters                            ***
       ***********************************************************************/
      if (strcmp(key, "verbosity") == 0){
	 READ_INT_PAR(p->par.verbosity);
	 tm_par->verbosity = lp_par->verbosity = cg_par->verbosity =
	 /*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
	    sp_par->verbosity =
#endif 
	 /*___END_EXPERIMENTAL_SECTION___*/
	    cp_par->verbosity = p->par.verbosity;
      }
      else if (strcmp(key, "random_seed") == 0){
	 READ_INT_PAR(p->par.random_seed);
	 tm_par->random_seed = p->par.random_seed;
      }
      else if (strcmp(key, "granularity") == 0){
	 READ_DBL_PAR(tm_par->granularity);
	 lp_par->granularity = tm_par->granularity;
      }
      /*__BEGIN_EXPERIMENTAL_SECTION__*/
      else if (strcmp(key, "do_decomp") == 0 ||
	       strcmp(key, "CG_do_decomp") == 0 ||
	       strcmp(key, "TM_do_decomp") == 0){
	 READ_INT_PAR(tm_par->do_decomp);
	 cg_par->do_decomp = tm_par->do_decomp;
      }
      /*___END_EXPERIMENTAL_SECTION___*/

      /***********************************************************************
       ***                    Master parameters                            ***
       ***********************************************************************/
      else if (strcmp(key, "upper_bound") == 0 ||
	       strcmp(key, "M_upper_bound") == 0){
	 READ_DBL_PAR(p->ub);
	 p->has_ub = TRUE;
      }
      else if (strcmp(key, "upper_bound_estimate") == 0 ||
	       strcmp(key, "M_upper_bound_estimate") == 0){
	 READ_DBL_PAR(p->ub_estimate);
	 p->has_ub_estimate = TRUE;
      }
      else if (strcmp(key, "lower_bound") == 0 ||
	       strcmp(key, "M_lower_bound") == 0){
	 READ_DBL_PAR(p->lb);
      }

      else if (strcmp(key, "M_verbosity") == 0){
	 READ_INT_PAR(p->par.verbosity);
      }
      else if (strcmp(key, "M_random_seed") == 0){
	 READ_INT_PAR(p->par.random_seed);
      }

      else if (strcmp(key, "tm_executable_name") == 0 ||
	       strcmp(key, "tm_exe") == 0 ||
	       strcmp(key, "M_tm_exe") == 0 ||
	       strcmp(key, "M_tm_executable_name") == 0){
	 read_string(p->par.tm_exe, line, MAX_FILE_NAME_LENGTH);
      }
      else if (strcmp(key, "dg_executable_name") == 0 ||
	       strcmp(key, "dg_exe") == 0 ||
	       strcmp(key, "M_dg_exe") == 0 ||
	       strcmp(key, "M_dg_executable_name") == 0){
	 read_string(p->par.dg_exe, line, MAX_FILE_NAME_LENGTH);
      }
      else if (strcmp(key, "tm_debug") == 0 ||
	       strcmp(key, "M_tm_debug") == 0){
	 READ_INT_PAR(p->par.tm_debug);
	 if (p->par.tm_debug) p->par.tm_debug = 4;
      }
      else if (strcmp(key, "dg_debug") == 0 ||
	       strcmp(key, "M_dg_debug") == 0){
	 READ_INT_PAR(p->par.dg_debug);
	 if (p->par.dg_debug) p->par.dg_debug = 4;
      }
      else if (strcmp(key, "tm_machine") == 0 ||
	       strcmp(key, "M_tm_machine") == 0){
	 read_string(p->par.tm_machine, line, MACH_NAME_LENGTH);
	 p->par.tm_machine_set = TRUE;
      }
      else if (strcmp(key, "dg_machine") == 0 ||
	       strcmp(key, "M_dg_machine") == 0){
	 read_string(p->par.dg_machine, line, MACH_NAME_LENGTH);
	 p->par.dg_machine_set = TRUE;
      }

      else if (strcmp(key, "pvm_trace") == 0 ||
	       strcmp(key, "M_pvm_trace") == 0){
	 READ_INT_PAR(p->par.pvm_trace);
      }
      else if (strcmp(key, "do_branch_and_cut") == 0 ||
	       strcmp(key, "M_do_branch_and_cut") == 0){
	 READ_INT_PAR(p->par.do_branch_and_cut);
      }
      else if (strcmp(key, "do_draw_graph") == 0 ||
	       strcmp(key, "M_do_draw_graph") == 0){
	 READ_INT_PAR(p->par.do_draw_graph);
      }

      /***********************************************************************
       ***                 DrawGraph parameters                            ***
       ***********************************************************************/

      else if (strcmp(key, "source_path") == 0 ||
	       strcmp(key, "DG_source_path") == 0){
	 read_string(dg_par->source_path, line, MAX_FILE_NAME_LENGTH);
      }
      else if (strcmp(key, "echo_commands") == 0 ||
	       strcmp(key, "DG_echo_commands") == 0){
	 READ_INT_PAR(dg_par->echo_commands);
      }
      else if (strcmp(key, "canvas_width") == 0 ||
	       strcmp(key, "DG_canvas_width") == 0){
	 READ_INT_PAR(dg_par->canvas_width);
      }
      else if (strcmp(key, "canvas_height") == 0 ||
	       strcmp(key, "DG_canvas_height") == 0){
	 READ_INT_PAR(dg_par->canvas_height);
      }
      else if (strcmp(key, "viewable_width") == 0 ||
	       strcmp(key, "DG_viewable_width") == 0){
	 READ_INT_PAR(dg_par->viewable_width);
      }
      else if (strcmp(key, "viewable_height") == 0 ||
	       strcmp(key, "DG_viewable_height") == 0){
	 READ_INT_PAR(dg_par->viewable_width);
      }
      else if (strcmp(key, "disp_nodelabels") == 0 ||
	       strcmp(key, "DG_disp_nodelabels") == 0){
	 READ_INT_PAR(dg_par->disp_nodelabels);
      }
      else if (strcmp(key, "disp_nodeweights") == 0 ||
	       strcmp(key, "DG_disp_nodeweights") == 0){
	 READ_INT_PAR(dg_par->disp_nodeweights);
      }
      else if (strcmp(key, "disp_edgeweights") == 0 ||
	       strcmp(key, "DG_disp_edgeweights") == 0){
	 READ_INT_PAR(dg_par->disp_edgeweights);
      }
      else if (strcmp(key, "node_dash") == 0 ||
	       strcmp(key, "DG_node_dash") == 0){
	 read_string(dg_par->node_dash, line, MAX_DASH_PATTERN_LENGTH);
      }
      else if (strcmp(key, "edge_dash") == 0 ||
	       strcmp(key, "DG_edge_dash") == 0){
	 read_string(dg_par->edge_dash, line, MAX_DASH_PATTERN_LENGTH);
      }
      else if (strcmp(key, "node_radius") == 0 ||
	       strcmp(key, "DG_node_radius") == 0){
	 READ_INT_PAR(dg_par->node_radius);
      }
      else if (strcmp(key, "interactive_mode") == 0 ||
	       strcmp(key, "DG_interactive_mode") == 0){
	 READ_INT_PAR(dg_par->interactive_mode);
      }
      else if (strcmp(key, "mouse_tracking") == 0 ||
	       strcmp(key, "DG_mouse_tracking") == 0){
	 READ_INT_PAR(dg_par->mouse_tracking);
      }
      else if (strcmp(key, "scale_factor") == 0 ||
	       strcmp(key, "DG_scale_factor") == 0){
	 READ_DBL_PAR(dg_par->scale_factor);
      }
      else if (strcmp(key, "nodelabel_font") == 0 ||
	       strcmp(key, "DG_nodelabel_font") == 0){
	 read_string(dg_par->nodelabel_font, line, MAX_FONT_LENGTH);
      }
      else if (strcmp(key, "nodeweight_font") == 0 ||
	       strcmp(key, "DG_nodeweight_font") == 0){
	 read_string(dg_par->nodeweight_font, line, MAX_FONT_LENGTH);
      }
      else if (strcmp(key, "edgeweight_font") == 0 ||
	       strcmp(key, "DG_edgeweight_font") == 0){
	 read_string(dg_par->edgeweight_font, line, MAX_FONT_LENGTH);
      }

      /***********************************************************************
       ***                  Treemanager parameters                         ***
       ***********************************************************************/
      else if (strcmp(key, "TM_verbosity") == 0){
	 READ_INT_PAR(tm_par->verbosity);
      }
      else if (strcmp(key, "TM_granularity") == 0){
	 READ_DBL_PAR(tm_par->granularity);
	 lp_par->granularity = tm_par->granularity;
      }
      else if (strcmp(key, "lp_executable_name") == 0 ||
	       strcmp(key, "lp_exe") == 0 ||
	       strcmp(key, "TM_lp_exe") == 0 ||
	       strcmp(key, "TM_lp_executable_name") == 0){
	 read_string(tm_par->lp_exe, line, MAX_FILE_NAME_LENGTH);
      }
      else if (strcmp(key, "cg_executable_name") == 0 ||
	       strcmp(key, "cg_exe") == 0 ||
	       strcmp(key, "TM_cg_exe") == 0 ||
	       strcmp(key, "TM_cg_executable_name") == 0){
	 read_string(tm_par->cg_exe, line, MAX_FILE_NAME_LENGTH);
      }
      else if (strcmp(key, "cp_executable_name") == 0 ||
	       strcmp(key, "cp_exe") == 0 ||
	       strcmp(key, "TM_cp_exe") == 0 ||
	       strcmp(key, "TM_cp_executable_name") == 0){
	 read_string(tm_par->cp_exe, line, MAX_FILE_NAME_LENGTH);
      }
      /*__BEGIN_EXPERIMENTAL_SECTION__*/
      else if (strcmp(key, "sp_executable_name") == 0 ||
	       strcmp(key, "sp_exe") == 0 ||
	       strcmp(key, "TM_sp_exe") == 0 ||
	       strcmp(key, "TM_sp_executable_name") == 0){
	 read_string(tm_par->sp_exe, line, MAX_FILE_NAME_LENGTH);
      }
      /*___END_EXPERIMENTAL_SECTION___*/
      else if (strcmp(key, "lp_debug") == 0 ||
	       strcmp(key, "TM_lp_debug") == 0){
	 READ_INT_PAR(tm_par->lp_debug);
	 if (tm_par->lp_debug) tm_par->lp_debug = 4;
      }
      else if (strcmp(key, "cg_debug") == 0 ||
	       strcmp(key, "TM_cg_debug") == 0){
	 READ_INT_PAR(tm_par->cg_debug);
	 if (tm_par->cg_debug) tm_par->cg_debug = 4;
      }
      else if (strcmp(key, "cp_debug") == 0 ||
	       strcmp(key, "TM_cp_debug") == 0){
	 READ_INT_PAR(tm_par->cp_debug);
	 if (tm_par->cp_debug) tm_par->cp_debug = 4;
      }
      /*__BEGIN_EXPERIMENTAL_SECTION__*/
      else if (strcmp(key, "sp_debug") == 0 ||
	       strcmp(key, "TM_sp_debug") == 0){
	 READ_INT_PAR(tm_par->sp_debug);
	 if (tm_par->sp_debug) tm_par->sp_debug = 4;
      }
      /*___END_EXPERIMENTAL_SECTION___*/
      else if (strcmp(key, "max_active_nodes") == 0 ||
	       strcmp(key, "TM_max_active_nodes") == 0){
	 READ_INT_PAR(tm_par->max_active_nodes);
      }
      else if (strcmp(key, "max_cp_num") == 0 ||
	       strcmp(key, "TM_max_cp_num") == 0){
	 READ_INT_PAR(tm_par->max_cp_num);
      }
      /*__BEGIN_EXPERIMENTAL_SECTION__*/
      else if (strcmp(key, "max_sp_num") == 0 ||
	       strcmp(key, "TM_max_sp_num") == 0){
	 READ_INT_PAR(tm_par->max_sp_num);
      }
      /*___END_EXPERIMENTAL_SECTION___*/
      else if (strcmp(key, "lp_mach_num") == 0 ||
	       strcmp(key, "TM_lp_mach_num") == 0){
	 READ_INT_PAR(tm_par->lp_mach_num);
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
		  exit(1);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "TM_lp_machine") != 0){
		  fprintf(stderr, "\nio: error reading lp_machine list\n\n");
		  exit(1);
	       }
	       read_string(tm_par->lp_machs[i], line, MACH_NAME_LENGTH);
	       printf("%s", line);
	    }
	 }
      }
      else if (strcmp(key, "cg_mach_num") == 0 ||
	       strcmp(key, "TM_cg_mach_num") == 0){
	 READ_INT_PAR(tm_par->cg_mach_num);
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
		  exit(1);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "TM_cg_machine") != 0){
		  fprintf(stderr, "\nio: error reading cg_machine list\n\n");
		  exit(1);
	       }
	       read_string(tm_par->cg_machs[i], line, MACH_NAME_LENGTH);
	       printf("%s", line);
	    }
	 }
      }
      else if (strcmp(key, "cp_mach_num") == 0 ||
	       strcmp(key, "TM_cp_mach_num") == 0){
	 READ_INT_PAR(tm_par->cp_mach_num);
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
		  exit(1);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "TM_cp_machine") != 0){
		  fprintf(stderr, "\nio: error reading cp_machine list\n\n");
		  exit(1);
	       }
	       read_string(tm_par->cp_machs[i], line, MACH_NAME_LENGTH);
	       printf("%s", line);
	    }
	 }
      }
#ifndef COMPILE_IN_CG
      else if (strcmp(key, "use_cg") == 0 ||
	       strcmp(key, "TM_use_cg") == 0 ||
	       strcmp(key, "LP_use_cg") == 0){
	 READ_INT_PAR(tm_par->use_cg);
	 lp_par->use_cg = tm_par->use_cg;
      }
#endif
      else if (strcmp(key, "TM_random_seed") == 0){
	 READ_INT_PAR(tm_par->random_seed);
      }
      else if (strcmp(key, "unconditional_dive_frac") == 0 ||
	       strcmp(key, "TM_unconditional_dive_frac") == 0){
	 READ_DBL_PAR(tm_par->unconditional_dive_frac);
      }
      else if (strcmp(key, "diving_strategy") == 0 ||
	       strcmp(key, "TM_diving_strategy") == 0){
	 READ_INT_PAR(tm_par->diving_strategy);
      }
      else if (strcmp(key, "diving_k") == 0 ||
	       strcmp(key, "TM_diving_k") == 0){
	 READ_INT_PAR(tm_par->diving_k);
      }
      else if (strcmp(key, "diving_threshold") == 0 ||
	       strcmp(key, "TM_diving_threshold") == 0){
	 READ_DBL_PAR(tm_par->diving_threshold);
      }
      else if (strcmp(key, "node_selection_rule") == 0 ||
	       strcmp(key, "TM_node_selection_rule") == 0){
	 READ_INT_PAR(tm_par->node_selection_rule);
      }
      else if (strcmp(key, "keep_description_of_pruned") == 0 ||
	       strcmp(key, "TM_keep_description_of_pruned") == 0){
	 READ_INT_PAR(tm_par->keep_description_of_pruned);
	 if (tm_par->keep_description_of_pruned == KEEP_ON_DISK_FULL ||
	     tm_par->keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No pruned node file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "pruned_node_file_name") != 0){
	       printf("Need pruned_node_file_name next!!!\n\n");
	       exit(1);
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
	 READ_INT_PAR(tm_par->warm_start);
	 if ((p->par.warm_start = tm_par->warm_start)){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No warm start tree file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "warm_start_tree_file_name") != 0){
	       printf("Need warm_start_tree_file_name next!!!\n\n");
	       exit(1);
	    }
	    strcpy(tm_par->warm_start_tree_file_name, value);
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No warm start cut file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "warm_start_cut_file_name") != 0){
	       printf("Need warm_start_cut_file_name next!!!\n\n");
	       exit(1);
	    }
	    strcpy(tm_par->warm_start_cut_file_name, value);
	 }
      }
      else if (strcmp(key, "vbc_emulation") == 0 ||
	       strcmp(key, "TM_vbc_emulation") == 0){
	 READ_INT_PAR(tm_par->vbc_emulation);
	 if (tm_par->vbc_emulation == VBC_EMULATION_FILE){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No vbc emulation file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "vbc_emulation_file_name") != 0){
	       printf("Need vbc_emulation_file_name next!!!\n\n");
	       exit(1);
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
      else if (strcmp(key, "logging_interval") == 0 ||
	       strcmp(key, "TM_logging_interval") == 0){
	 READ_INT_PAR(tm_par->logging_interval);
      }
      else if (strcmp(key, "logging") == 0 ||
	       strcmp(key, "TM_logging") == 0){
	 READ_INT_PAR(tm_par->logging);
	 if (tm_par->logging){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No tree log file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "tree_log_file_name") != 0){
	       printf("tree_log_file_name next!!!\n\n");
	       exit(1);
	    }
	    strcpy(tm_par->tree_log_file_name, value);
	    if (tm_par->logging != VBC_TOOL){
	       if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
		  printf("No cut log file!\n\n");
		  exit(1);
	       }
	       strcpy(key, "");
	       sscanf(line, "%s%s", key, value);
	       if (strcmp(key, "cut_log_file_name") != 0){
		  printf("Need cut_log_file_name next!!!\n\n");
		  exit(1);
	       }
	       strcpy(tm_par->cut_log_file_name, value);
	    }
	 }
      }
      else if (strcmp(key, "price_in_root") == 0 ||
	       strcmp(key, "TM_price_in_root") == 0){
	 READ_INT_PAR(tm_par->price_in_root);
      }
      else if (strcmp(key, "trim_search_tree") == 0 ||
	       strcmp(key, "TM_trim_search_tree") == 0){
	 READ_INT_PAR(tm_par->trim_search_tree);
      }
      else if (strcmp(key, "colgen_in_first_phase") == 0 ||
	       strcmp(key, "TM_colgen_in_first_phase") == 0){
	 READ_INT_PAR(tm_par->colgen_strat[0]);
      }
      else if (strcmp(key, "colgen_in_second_phase") == 0 ||
	       strcmp(key, "TM_colgen_in_second_phase") == 0){
	 READ_INT_PAR(tm_par->colgen_strat[1]);
      }
      else if (strcmp(key, "colgen_in_first_phase_str") == 0 ||
	       strcmp(key, "TM_colgen_in_first_phase_str") == 0){
	 READ_STRINT_PAR(tm_par->colgen_strat[0],
			 colgen_str, COLGEN_STR_SIZE, value);
      }
      else if (strcmp(key, "colgen_in_second_phase_str") == 0 ||
	       strcmp(key, "TM_colgen_in_second_phase_str") == 0){
	 READ_STRINT_PAR(tm_par->colgen_strat[1],
			 colgen_str, COLGEN_STR_SIZE, value);
      }
      else if (strcmp(key, "time_limit") == 0 ||
	       strcmp(key, "TM_time_limit") == 0){
	 READ_DBL_PAR(tm_par->time_limit);
      }

      /***********************************************************************
       ***                      LP parameters                              ***
       ***********************************************************************/
      else if (strcmp(key, "LP_verbosity") == 0){
	 READ_INT_PAR(lp_par->verbosity);
      }
      else if (strcmp(key, "LP_granularity") == 0){
	 READ_DBL_PAR(lp_par->granularity);
	 tm_par->granularity = lp_par->granularity;
      }
      else if (strcmp(key, "set_obj_upper_lim") == 0 ||
	       strcmp(key, "LP_set_obj_upper_lim") == 0){
	 READ_INT_PAR(lp_par->set_obj_upper_lim);
      }

      else if (strcmp(key, "scaling") == 0 ||
	       strcmp(key, "LP_scaling") == 0){
	 READ_INT_PAR(lp_par->scaling);
      }
      else if (strcmp(key, "fastmip") == 0 ||
	       strcmp(key, "LP_fastmip") == 0){
	 READ_INT_PAR(lp_par->fastmip);
      }
      else if (strcmp(key, "try_to_recover_from_error") == 0 ||
	       strcmp(key, "LP_try_to_recover_from_error") == 0){
	 READ_INT_PAR(lp_par->try_to_recover_from_error);
      }
      else if (strcmp(key, "problem_type") == 0 ||
	       strcmp(key, "LP_problem_type") == 0){
	 READ_INT_PAR(lp_par->problem_type);
      }
      else if (strcmp(key, "not_fixed_storage_size") == 0 ||
	       strcmp(key, "LP_not_fixed_storage_size") == 0 ||
	       strcmp(key, "TM_not_fixed_storage_size") == 0 ){
	 READ_INT_PAR(lp_par->not_fixed_storage_size);
	 tm_par->not_fixed_storage_size = lp_par->not_fixed_storage_size;
      }
      else if (strcmp(key, "cut_pool_check_frequency") == 0 ||
	       strcmp(key, "LP_cut_pool_check_frequency") == 0){
	 READ_INT_PAR(lp_par->cut_pool_check_freq);
      }
      else if (strcmp(key, "load_balance_level") == 0 ||
	       strcmp(key, "LP_load_balance_level") == 0){
	 READ_INT_PAR(lp_par->load_balance_level);
      }
      else if (strcmp(key, "load_balance_iterations") == 0 ||
	       strcmp(key, "LP_load_balance_iterations") == 0){
	 READ_INT_PAR(lp_par->load_balance_iterations);
      }
      else if (strcmp(key, "load_balance_compare_candidates") == 0 ||
	       strcmp(key, "LP_load_balance_compare_candidates") == 0){
	 READ_INT_PAR(lp_par->load_balance_compare_candidates);
      }
      else if (strcmp(key, "fractional_diving_ratio") == 0 ||
	       strcmp(key, "LP_fractional_diving_ratio") == 0){
	 READ_DBL_PAR(lp_par->fractional_diving_ratio);
      }
      else if (strcmp(key, "fractional_diving_num") == 0 ||
	       strcmp(key, "LP_fractional_diving_num") == 0){
	 READ_INT_PAR(lp_par->fractional_diving_num);
      }
      else if (strcmp(key, "max_non_dual_feas_to_add_frac") == 0 ||
	       strcmp(key, "LP_max_non_dual_feas_to_add_frac") == 0){
	 READ_DBL_PAR(lp_par->max_non_dual_feas_to_add_frac);
      }
      else if (strcmp(key, "max_cols_to_add_min") == 0 ||
	       strcmp(key, "LP_max_non_dual_feas_to_add_min") == 0){
	 READ_INT_PAR(lp_par->max_non_dual_feas_to_add_min);
      }
      else if (strcmp(key, "max_non_dual_feas_to_add_max") == 0 ||
	       strcmp(key, "LP_max_non_dual_feas_to_add_max") == 0){
	 READ_INT_PAR(lp_par->max_non_dual_feas_to_add_max);
      }
      else if (strcmp(key, "max_not_fixable_to_add_frac") == 0 ||
	       strcmp(key, "LP_max_not_fixable_to_add_frac") == 0){
	 READ_DBL_PAR(lp_par->max_not_fixable_to_add_frac);
      }
      else if (strcmp(key, "max_not_fixable_to_add_min") == 0 ||
	       strcmp(key, "LP_max_not_fixable_to_add_min") == 0){
	 READ_INT_PAR(lp_par->max_not_fixable_to_add_min);
      }
      else if (strcmp(key, "max_not_fixable_to_add_max") == 0 ||
	       strcmp(key, "LP_max_not_fixable_to_add_max") == 0){
	 READ_INT_PAR(lp_par->max_not_fixable_to_add_max);
      }

      else if (strcmp(key, "mat_col_compress_num") == 0 ||
	       strcmp(key, "LP_mat_col_compress_num") == 0){
	 READ_INT_PAR(lp_par->mat_col_compress_num);
      }
      else if (strcmp(key, "mat_col_compress_ratio") == 0 ||
	       strcmp(key, "LP_mat_col_compress_ratio") == 0){
	 READ_DBL_PAR(lp_par->mat_col_compress_ratio);
      }
      else if (strcmp(key, "mat_row_compress_num") == 0 ||
	       strcmp(key, "LP_mat_row_compress_num") == 0){
	 READ_INT_PAR(lp_par->mat_row_compress_num);
      }
      else if (strcmp(key, "mat_row_compress_ratio") == 0 ||
	       strcmp(key, "LP_mat_row_compress_ratio") == 0){
	 READ_DBL_PAR(lp_par->mat_row_compress_ratio);
      }

      else if (strcmp(key, "tailoff_gap_backsteps") == 0 ||
	       strcmp(key, "LP_tailoff_gap_backsteps") == 0){
	 READ_INT_PAR(lp_par->tailoff_gap_backsteps);
      }
      else if (strcmp(key, "tailoff_obj_backsteps") == 0 ||
	       strcmp(key, "LP_tailoff_obj_backsteps") == 0){
	 READ_INT_PAR(lp_par->tailoff_obj_backsteps);
      }
      else if (strcmp(key, "tailoff_gap_frac") == 0 ||
	       strcmp(key, "LP_tailoff_gap_frac") == 0){
	 READ_DBL_PAR(lp_par->tailoff_gap_frac);
      }
      else if (strcmp(key, "tailoff_obj_frac") == 0 ||
	       strcmp(key, "LP_tailoff_obj_frac") == 0){
	 READ_DBL_PAR(lp_par->tailoff_obj_frac);
      }

     else if (strcmp(key, "ineff_cnt_to_delete") == 0 ||
	      strcmp(key, "LP_ineff_cnt_to_delete") == 0){
	 READ_INT_PAR(lp_par->ineff_cnt_to_delete);
      }
      else if (strcmp(key, "eff_cnt_before_cutpool") == 0 ||
	       strcmp(key, "LP_eff_cnt_before_cutpool") == 0){
	 READ_INT_PAR(lp_par->eff_cnt_before_cutpool);
      }
      else if (strcmp(key, "ineffective_constraints") == 0 ||
	       strcmp(key, "LP_ineffective_constraints") == 0){
	 READ_INT_PAR(lp_par->ineffective_constraints);
      }
      else if (strcmp(key, "base_constraints_always_effective") == 0 ||
	       strcmp(key, "LP_base_constraints_always_effective") == 0){
	 READ_INT_PAR(lp_par->base_constraints_always_effective);
      }

      else if (strcmp(key, "branch_on_cuts") == 0 ||
	       strcmp(key, "LP_branch_on_cuts") == 0){
	 READ_INT_PAR(lp_par->branch_on_cuts);
      }
      else if (strcmp(key, "discard_slack_cuts") == 0 ||
	       strcmp(key, "LP_discard_slack_cuts") == 0){
	 READ_INT_PAR(lp_par->discard_slack_cuts);
      }

      /* timeouts on receiving cuts */
      else if (strcmp(key, "first_lp_first_cut_time_out") == 0 ||
	       strcmp(key, "LP_first_lp_first_cut_time_out") == 0){
	 READ_DBL_PAR(timeout);
	 if (timeout == -1){
	    lp_par->first_lp.first_cut_time_out = 0;
	 }else{
	    lp_par->first_lp.first_cut_time_out = timeout;
	 }
      }
      else if (strcmp(key, "first_lp_all_cuts_time_out") == 0 ||
	       strcmp(key, "LP_first_lp_all_cuts_time_out") == 0){
	 READ_DBL_PAR(timeout);
	 if (timeout == -1){
	    lp_par->first_lp.all_cuts_time_out = 0;
	 }else{
	    lp_par->first_lp.all_cuts_time_out = timeout;
	 }
      }
      else if (strcmp(key, "later_lp_first_cut_time_out") == 0 ||
	       strcmp(key, "LP_later_lp_first_cut_time_out") == 0){
	 READ_DBL_PAR(timeout);
	 if (timeout == -1){
	    lp_par->later_lp.first_cut_time_out = 0;
	 }else{
	   lp_par->later_lp.first_cut_time_out = timeout;
	 }
      }
      else if (strcmp(key, "later_lp_all_cuts_time_out") == 0 ||
	       strcmp(key, "LP_later_lp_all_cuts_time_out") == 0){
	 READ_DBL_PAR(timeout);
	 if (timeout == -1){
	    lp_par->later_lp.all_cuts_time_out = 0;
	 }else{
	    lp_par->later_lp.all_cuts_time_out = timeout;
	 }
      }

      else if (strcmp(key, "no_cut_timeout") == 0 ||
	       strcmp(key, "LP_no_cut_timeout") == 0){
	 lp_par->first_lp.first_cut_time_out = 0;
	 lp_par->first_lp.all_cuts_time_out = 0;
	 lp_par->later_lp.first_cut_time_out = 0;
	 lp_par->later_lp.all_cuts_time_out = 0;
	 /*__BEGIN_EXPERIMENTAL_SECTION__*/
	 cg_par->decomp_dynamic_timeout = 6000;
	 /*___END_EXPERIMENTAL_SECTION___*/
      }
      else if (strcmp(key, "all_cut_timeout") == 0 ||
	       strcmp(key, "LP_all_cut_timeout") == 0){
	 READ_DBL_PAR(timeout);
	 lp_par->first_lp.first_cut_time_out = timeout;
	 lp_par->first_lp.all_cuts_time_out = timeout;
	 lp_par->later_lp.first_cut_time_out= timeout;
	 lp_par->later_lp.all_cuts_time_out = timeout;
	 /*__BEGIN_EXPERIMENTAL_SECTION__*/
	 cg_par->decomp_dynamic_timeout = timeout;
	 /*___END_EXPERIMENTAL_SECTION___*/
      }

      else if (strcmp(key, "max_cut_num_per_iter") == 0 ||
	       strcmp(key, "LP_max_cut_num_per_iter") == 0){
	 READ_INT_PAR(lp_par->max_cut_num_per_iter);
      }

      /* variable fixing parameters */
      else if (strcmp(key, "do_reduced_cost_fixing") == 0 ||
	       strcmp(key, "LP_do_reduced_cost_fixing") == 0){
	 READ_INT_PAR(lp_par->do_reduced_cost_fixing);
      }
      else if (strcmp(key, "gap_as_ub_frac") == 0 ||
	       strcmp(key, "LP_gap_as_ub_frac") == 0){
	 READ_DBL_PAR(lp_par->gap_as_ub_frac);
      }
      else if (strcmp(key, "gap_as_last_gap_frac") == 0 ||
	       strcmp(key, "LP_gap_as_last_gap_frac") == 0){
	 READ_DBL_PAR(lp_par->gap_as_last_gap_frac);
      }
      else if (strcmp(key, "do_logical_fixing") == 0 ||
	       strcmp(key, "LP_do_logical_fixing") == 0){
	 READ_INT_PAR(lp_par->do_logical_fixing);
      }
      else if (strcmp(key, "fixed_to_ub_before_logical_fixing") == 0 ||
	       strcmp(key, "LP_fixed_to_ub_before_logical_fixing") == 0){
	 READ_INT_PAR(lp_par->fixed_to_ub_before_logical_fixing);
      }
      else if (strcmp(key, "fixed_to_ub_frac_before_logical_fixing")==0 ||
	       strcmp(key, "LP_fixed_to_ub_frac_before_logical_fixing")==0){
	 READ_DBL_PAR(lp_par->fixed_to_ub_frac_before_logical_fixing);
      }

      else if (strcmp(key, "max_presolve_iter") == 0 ||
	       strcmp(key, "LP_max_presolve_iter") == 0){
	 READ_INT_PAR(lp_par->max_presolve_iter);
      }

      /* user-defined function defaults */
      else if (strcmp(key, "is_feasible_default") == 0 ||
	       strcmp(key, "LP_is_feasible_default") == 0){
	 READ_INT_PAR(lp_par->is_feasible_default);
      }
      else if (strcmp(key, "send_feasible_solution_default") == 0 ||
	       strcmp(key, "LP_send_feasible_solution_default") == 0){
	 READ_INT_PAR(lp_par->send_feasible_solution_default);
      }
      else if (strcmp(key, "display_solution_default") == 0 ||
	       strcmp(key, "LP_display_solution_default") == 0){
	 READ_INT_PAR(lp_par->display_solution_default);
      }
      else if (strcmp(key, "shall_we_branch_default") == 0 ||
	       strcmp(key, "LP_shall_we_branch_default") == 0){
	 READ_INT_PAR(lp_par->shall_we_branch_default);
      }
      else if (strcmp(key, "select_candidates_default") == 0 ||
	       strcmp(key, "LP_select_candidates_default") == 0){
	 READ_INT_PAR(lp_par->select_candidates_default);
      }
      else if (strcmp(key, "strong_branching_cand_num") == 0){
	 READ_INT_PAR(lp_par->strong_branching_cand_num_max);
	 lp_par->strong_branching_cand_num_min =
	    lp_par->strong_branching_cand_num_max;
	 lp_par->strong_branching_red_ratio = 0;
      }
      else if (strcmp(key, "strong_branching_cand_num_max") == 0 ||
	       strcmp(key, "LP_strong_branching_cand_num_max") == 0){
	 READ_INT_PAR(lp_par->strong_branching_cand_num_max);
      }
      else if (strcmp(key, "strong_branching_cand_num_min") == 0 ||
	       strcmp(key, "LP_strong_branching_cand_num_min") == 0){
	 READ_INT_PAR(lp_par->strong_branching_cand_num_min);
      }
      else if (strcmp(key,"strong_branching_red_ratio") == 0 ||
	       strcmp(key,"LP_strong_branching_red_ratio") == 0){
	 READ_DBL_PAR(lp_par->strong_branching_red_ratio);
      }
      else if (strcmp(key, "compare_candidates_default") == 0 ||
	       strcmp(key, "LP_compare_candidates_default") == 0){
	 READ_INT_PAR(lp_par->compare_candidates_default);
      }
      else if (strcmp(key, "compare_candidates_default_str") == 0 ||
	       strcmp(key, "LP_compare_candidates_default_str") == 0){
	 READ_STRINT_PAR(lp_par->compare_candidates_default,
			 compare_can_str, COMPARE_CAN_STR_SIZE, value);
      }
      else if (strcmp(key, "select_child_default") == 0 ||
	       strcmp(key, "LP_select_child_default") == 0){
	 READ_INT_PAR(lp_par->select_child_default);
      }
      else if (strcmp(key, "pack_lp_solution_default") == 0 ||
	       strcmp(key, "LP_pack_lp_solution_default") == 0){
	 READ_INT_PAR(lp_par->pack_lp_solution_default);
      }

      /***********************************************************************
       ***                     cut_gen parameters                          ***
       ***********************************************************************/
      else if (strcmp(key, "CG_verbosity") == 0){
	 READ_INT_PAR(cg_par->verbosity);
      }
      else if (strcmp(key, "do_findcuts") == 0 ||
	       strcmp(key, "CG_do_findcuts") == 0){
	 READ_INT_PAR(cg_par->do_findcuts);
      }
      /*__BEGIN_EXPERIMENTAL_SECTION__*/
      else if (strcmp(key, "decomp_sol_pool_check_freq") == 0 ||
	       strcmp(key, "CG_decomp_sol_pool_check_freq") == 0){
	 READ_INT_PAR(cg_par->decomp_sol_pool_check_freq);
      }
      else if (strcmp(key, "decomp_wait_for_cols") == 0 ||
	       strcmp(key, "CG_decomp_wait_for_cols") == 0){
	 READ_INT_PAR(cg_par->decomp_wait_for_cols);
      }
      else if (strcmp(key, "decomp_max_col_num_per_iter") == 0 ||
	       strcmp(key, "CG_decomp_max_col_num_per_iter") == 0){
	 READ_INT_PAR(cg_par->decomp_max_col_num_per_iter);
      }
      else if (strcmp(key, "decomp_col_block_size") == 0 ||
	       strcmp(key, "CG_decomp_col_block_size") == 0){
	 READ_INT_PAR(cg_par->decomp_col_block_size);
      }
      else if (strcmp(key, "decomp_mat_block_size") == 0 ||
	       strcmp(key, "CG_decomp_mat_block_size") == 0){
	 READ_INT_PAR(cg_par->decomp_mat_block_size);
      }
      else if (strcmp(key, "decomp_initial_timeout") == 0 ||
	       strcmp(key, "CG_decomp_initial_timeout") == 0){
	 READ_DBL_PAR(cg_par->decomp_initial_timeout);
      }
      else if (strcmp(key, "decomp_dynamic_timeout") == 0 ||
	       strcmp(key, "CG_decomp_dynamic_timeout") == 0){
	 READ_DBL_PAR(cg_par->decomp_dynamic_timeout);
      }
      else if (strcmp(key, "decomp_complete_enum") == 0 ||
	       strcmp(key, "CG_decomp_complete_enum") == 0){
	 READ_INT_PAR(cg_par->decomp_complete_enum);
      }
      /*___END_EXPERIMENTAL_SECTION___*/

      /***********************************************************************
       ***                      cutpool parameters                         ***
       ***********************************************************************/
      else if (strcmp(key, "CP_verbosity") == 0){
	 READ_INT_PAR(cp_par->verbosity);
      }
      else if (strcmp(key, "cp_warm_start") == 0 ||
	       strcmp(key, "CP_warm_start") == 0){
	 READ_INT_PAR(cp_par->warm_start);
	 if (cp_par->warm_start){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No cut pool warm start file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "cp_warm_start_file_name") != 0){
	       printf("Need cp_warm_start_file_name next!!!\n\n");
	       exit(1);
	    }
	    strcpy(cp_par->warm_start_file_name, value);
	 }
      }
      else if (strcmp(key, "cp_logging") == 0 ||
	       strcmp(key, "CP_logging") == 0){
	 READ_INT_PAR(cp_par->logging);
	 if ((tm_par->cp_logging = cp_par->logging)){
	    if (fgets(line, MAX_LINE_LENGTH, f) == NULL){
	       printf("No cut pool log file!\n\n");
	       exit(1);
	    }
	    strcpy(key, "");
	    sscanf(line, "%s%s", key, value);
	    if (strcmp(key, "cp_log_file_name") != 0){
	       printf("Need cp_log_file_name next!!!\n\n");
	       exit(1);
	    }
	    strcpy(cp_par->log_file_name, value);
	 }
      }
      else if (strcmp(key, "block_size") == 0 ||
	       strcmp(key, "CP_block_size") == 0){
	 READ_INT_PAR(cp_par->block_size);
      }
      else if (strcmp(key, "max_size") == 0 ||
	       strcmp(key, "CP_max_size") == 0){
	 READ_INT_PAR(cp_par->max_size);
      }
      else if (strcmp(key, "max_number_of_cuts") == 0 ||
	       strcmp(key, "CP_max_number_of_cuts") == 0){
	 READ_INT_PAR(cp_par->max_number_of_cuts);
      }
      else if (strcmp(key, "cuts_to_check") == 0 ||
	       strcmp(key, "cuts_to_check") == 0){
	 READ_INT_PAR(cp_par->cuts_to_check);
      }
      else if (strcmp(key, "delete_which") == 0 ||
	       strcmp(key, "CP_delete_which") == 0){
	 READ_INT_PAR(cp_par->delete_which);
      }
      else if (strcmp(key, "touches_until_deletion") == 0 ||
	       strcmp(key, "CP_touches_until_deletion") == 0){
	 READ_INT_PAR(cp_par->touches_until_deletion);
      }
      else if (strcmp(key, "min_to_delete") == 0 ||
	       strcmp(key, "CP_min_to_delete") == 0){
	 READ_INT_PAR(cp_par->min_to_delete);
      }
      else if (strcmp(key, "check_which") == 0 ||
	       strcmp(key, "CP_check_which") == 0){
	 READ_INT_PAR(cp_par->check_which);
      }
      /*__BEGIN_EXPERIMENTAL_SECTION__*/

      /***********************************************************************
       ***                     solpool parameters                          ***
       ***********************************************************************/
#ifdef COMPILE_DECOMP
      else if (strcmp(key, "SP_verbosity") == 0){
	 READ_INT_PAR(sp_par->verbosity);
      }
      else if (strcmp(key, "SP_etol") == 0){
	 READ_DBL_PAR(sp_par->etol);
      }

      else if (strcmp(key, "SP_block_size") == 0){
	 READ_INT_PAR(sp_par->block_size);
      }
      else if (strcmp(key, "SP_max_size") == 0){
	 READ_INT_PAR(sp_par->max_size);
      }
      else if (strcmp(key, "max_number_of_sols") == 0 ||
	       strcmp(key, "SP_max_number_of_sols") == 0){
	 READ_INT_PAR(sp_par->max_number_of_sols);
      }
      else if (strcmp(key, "SP_delete_which") == 0){
	 READ_INT_PAR(sp_par->delete_which);
      }
      else if (strcmp(key, "SP_touches_until_deletion") == 0){
	 READ_INT_PAR(sp_par->touches_until_deletion);
      }
      else if (strcmp(key, "SP_min_to_delete") == 0){
	 READ_INT_PAR(sp_par->min_to_delete);
      }
      else if (strcmp(key, "SP_compress_num") == 0){
	 READ_INT_PAR(sp_par->compress_num);
      }
      else if (strcmp(key, "SP_compress_ratio") == 0){
	 READ_DBL_PAR(sp_par->compress_ratio);
      }
      else if (strcmp(key, "SP_check_which") == 0){
	 READ_INT_PAR(sp_par->check_which);
      }
#endif
      /*___END_EXPERIMENTAL_SECTION___*/
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
	 user_usage();
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
	 p->par.do_draw_graph = TRUE;
	 break;
       case 'g':
	 lp_par->use_cg = tm_par->use_cg = TRUE;
	 break;
       case 'r':
	 tm_par->price_in_root = TRUE;
	 break;
       case 't':
	 tm_par->trim_search_tree = TRUE;
	 break;
       case 'b':
	 p->par.do_branch_and_cut = FALSE;
	 break;
       case 'u':
	 sscanf(argv[++i], "%lf", &p->ub);
	 p->has_ub = TRUE;
	 break;
       case 'p':
	 sscanf(argv[++i], "%i", &tm_par->max_active_nodes);
	 break;
       case 'n':
	 sscanf(argv[++i], "%i", &tm_par->node_selection_rule);
	 break;
       case 'v':
	 sscanf(argv[++i], "%i", &p->par.verbosity);
	 tm_par->verbosity = lp_par->verbosity = cg_par->verbosity =
	 /*__BEGIN_EXPERIMENTAL_SECTION__*/
#ifdef COMPILE_DECOMP
	    sp_par->verbosity =
#endif 
	 /*___END_EXPERIMENTAL_SECTION___*/
	    cp_par->verbosity = p->par.verbosity;
 	 break;
       case 's':
	 sscanf(argv[++i], "%i",
		&lp_par->strong_branching_cand_num_max);
	 lp_par->strong_branching_cand_num_min =
	    lp_par->strong_branching_cand_num_max;
	 lp_par->strong_branching_red_ratio = 0;
	 break;
       case 'c':
	 sscanf(argv[++i], "%i", &lp_par->compare_candidates_default);
	 break;
       case 'k':
	 sscanf(argv[++i], "%i", &lp_par->select_child_default);
	 break;
       case 'm':
	 sscanf(argv[++i], "%i", &lp_par->max_cut_num_per_iter);
	 break;
       case 'e':
	 sscanf(argv[++i], "%i", &tm_par->max_cp_num);
	 break;
       case 'l':
	 sscanf(argv[++i], "%i", &lp_par->load_balance_level);
	 sscanf(argv[++i], "%i", &lp_par->load_balance_iterations);
	 break;
       case 'i':
	 sscanf(argv[++i], "%i", &lp_par->max_presolve_iter);
	 break;
       case 'f':
	 strncpy(p->par.param_file, argv[++i], MAX_FILE_NAME_LENGTH);
	 break;
       case 'z':
	 sscanf(argv[++i], "%lf", &tm_par->diving_threshold);
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
}

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
