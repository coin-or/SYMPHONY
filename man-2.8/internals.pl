# LaTeX2HTML 98.1p1 release (March 2nd, 1998)
# Associate internals original text with physical files.


$key = q/user_check_validity_of_cut/;
$ref_files{$key} = "$dir".q|node75.html|; 
$noresave{$key} = "$nosave";

$key = q//;
$ref_files{$key} = "$dir".q|node91.html|; 
$noresave{$key} = "$nosave";

$key = q/user_select_child/;
$ref_files{$key} = "$dir".q|node60.html|; 
$noresave{$key} = "$nosave";

$key = q/diving/;
$ref_files{$key} = "$dir".q|node93.html|; 
$noresave{$key} = "$nosave";

$key = q/user_receive_cg_data/;
$ref_files{$key} = "$dir".q|node72.html|; 
$noresave{$key} = "$nosave";

$key = q/fractional_diving/;
$ref_files{$key} = "$dir".q|node93.html|; 
$noresave{$key} = "$nosave";

$key = q/tm_params/;
$ref_files{$key} = "$dir".q|node93.html|; 
$noresave{$key} = "$nosave";

$key = q/cut_pool_params/;
$ref_files{$key} = "$dir".q|node96.html|; 
$noresave{$key} = "$nosave";

$key = q/waiting_row/;
$ref_files{$key} = "$dir".q|node48.html|; 
$noresave{$key} = "$nosave";

$key = q/debugging/;
$ref_files{$key} = "$dir".q|node22.html|; 
$noresave{$key} = "$nosave";

$key = q/debugging-PVM/;
$ref_files{$key} = "$dir".q|node20.html|; 
$noresave{$key} = "$nosave";

$key = q/user_display_solution/;
$ref_files{$key} = "$dir".q|node42.html|; 
$noresave{$key} = "$nosave";

$key = q/user_send_feas_sol/;
$ref_files{$key} = "$dir".q|node43.html|; 
$noresave{$key} = "$nosave";

$key = q/user_usage/;
$ref_files{$key} = "$dir".q|node29.html|; 
$noresave{$key} = "$nosave";

$key = q/user_send_cp_data/;
$ref_files{$key} = "$dir".q|node41.html|; 
$noresave{$key} = "$nosave";

$key = q/user_compare_candidates/;
$ref_files{$key} = "$dir".q|node59.html|; 
$noresave{$key} = "$nosave";

$key = q/exe_names/;
$ref_files{$key} = "$dir".q|node17.html|; 
$noresave{$key} = "$nosave";

$key = q/user_receive_feasible_solution/;
$ref_files{$key} = "$dir".q|node38.html|; 
$noresave{$key} = "$nosave";

$key = q/shared/;
$ref_files{$key} = "$dir".q|node13.html|; 
$noresave{$key} = "$nosave";

$key = q/configuration/;
$ref_files{$key} = "$dir".q|node16.html|; 
$noresave{$key} = "$nosave";

$key = q/user_shall_we_branch/;
$ref_files{$key} = "$dir".q|node57.html|; 
$noresave{$key} = "$nosave";

$key = q/user_check_cuts/;
$ref_files{$key} = "$dir".q|node82.html|; 
$noresave{$key} = "$nosave";

$key = q/user_send_lp_data/;
$ref_files{$key} = "$dir".q|node39.html|; 
$noresave{$key} = "$nosave";

$key = q/upper_bound_estimate/;
$ref_files{$key} = "$dir".q|node91.html|; 
$noresave{$key} = "$nosave";

$key = q/cuts_to_check/;
$ref_files{$key} = "$dir".q|node96.html|; 
$noresave{$key} = "$nosave";

$key = q/diving_strategy_diving_strategy_--_integer__BEST_ESTIMATE_0__._The_strategy_employed_when_deciding_whether_to_dive_or_not.____The_BEST_ESTIMATE__0___strategy_continues_to_dive_until_the_lower_bound_in_the_child_to_be_dived_into_exceeds_the_parameter__htmlrefupper_bound_estimate_labelparameter_file_The_parameter_file_name_is_passed_to__SYMPHONY__as_the_only_command_line_argument_to_the_master_process_which_is_started_by_the_user._Each_line_of_the_parameter_file_contains_either_a_comment_or_two_words_--_a_keyword_and_a_value__separated_by_white_space._If_the_first_word__sequence_of_non-white-space_characters__on_a_line_is_not_a_keyword__then_the_line_is_considered_a_comment_line._Otherwise_the_parameter_corresponding_to_the_keyword_is_set_to_the_listed_value._Usually_the_keyword_is_the_same_as_the_parameter_name_in_the_source_code._Here_we_list_the_keywords__the_type_of_value_that_should_be_given_with_the_keywords_and_the_default_value._A_parameter_corresponding_to_keyword___K___in_process___P___can_also_be_set_by_using_the_keyword___P__K__.___To_make_this_list_shorter__occasionally_a_comma_separated_list_of_parameters_is_given_if_the_meanings_of_those_parameters_are_strongly_connected._For_clarity__the_constant_name_is_sometimes_given_instead_of_the_numerical_value_for_default_settings_and_options._The_corresponding_value_is_given_in_curly_braces_for_convenience.______item_TM_verbosity_--_integer__0_.__The_verbosity_of_the_TM_process.____item_lp_exe__cg_exe__cp_exe_--_strings____lp______cg______cp___.__The_name_of_the_LP__CG__and_CP_process_binaries._Note:_when_running_in_parallel_using_PVM__these_executables__or_links_to_them__must_reside_in_the_PVM_ROOT_bin_PVM_ARCH__directory._Also__be_sure_to_note_that_the_executable_names_may_have_extensions_that_depend_on_the_configuration_of_the_modules__but_the_defaults_will_always_be_set_to_the_name_that_the_make_file_produce.____item_lp_debug__cg_debug__cp_debug_--_boolean__all_FALSE_.__Whether_the_processes_should_be_started_under_a_debugger_or_not.____item_max_active_nodes_--_integer__1_.__The_maximum_number_of_active_search_tree_nodes---equal_to_the_number_of_LP_and_CG_tandems_to_be_started_up.____item_max_cp_num_--_integer__0_.__The_maximum_number_of_cut_pools_to_be_used.____item_lp_mach_num__cg_mach_num__cp_mach_num_--_integers__all_0_.____The_number_of_processors_in_the_virtual_machine_to_run_LP__CG__CP__processes._If_this_value_is_0_then_the_processes_will_be_assigned_to_processors_in_round-robin_order._Otherwise_the_next_xx_mach_num_lines_describe_the_processors_where_the_LP__CG__CP__processes_must_run._The_keyword_--_value_pairs_on_these_lines_must_be_TM_xx_machine_and_the_name_or_IP_address_of_a_processor__the_processor_names_need_not_be_distinct_._In_this_case_the_actual_processes_are_assigned_in_a_round_robin_fashion_to_the_processors_on_this_list.________This_feature_is_useful_if_a_specific_software_package_is_needed_for_some_process__but_that_software_is_not_licensed_for_every_node_of_the_virtual_machine_or_if_a_certain_process_must_run_on_a_certain_type_of_machine_due_to_resource_requirements.____item_use_cg_--_boolean__FALSE_.__Whether_to_use_a_cut_generator_or_not._____item_TM_random_seed_--_integer__17_.__The_random_seed_used_in_the_TM.____item_unconditional_dive_frac_--_double__0.1_.__The_fraction_of_the_nodes_on_which__SYMPHONY__randomly_dives_unconditionally_into_one_of_the_children.____labeldiving_strategy/;
$ref_files{$key} = "$dir".q|node91.html|; 
$noresave{$key} = "$nosave";

$key = q/parameter_file/;
$ref_files{$key} = "$dir".q|node89.html|; 
$noresave{$key} = "$nosave";

$key = q/user-written-lp/;
$ref_files{$key} = "$dir".q|node45.html|; 
$noresave{$key} = "$nosave";

$key = q/delete_which/;
$ref_files{$key} = "$dir".q|node96.html|; 
$noresave{$key} = "$nosave";

$key = q/cut_data/;
$ref_files{$key} = "$dir".q|node47.html|; 
$noresave{$key} = "$nosave";

$key = q/strong_branching/;
$ref_files{$key} = "$dir".q|node94.html|; 
$noresave{$key} = "$nosave";

$key = q/getting_started/;
$ref_files{$key} = "$dir".q|node5.html|; 
$noresave{$key} = "$nosave";

$key = q/user_send_cg_data/;
$ref_files{$key} = "$dir".q|node40.html|; 
$noresave{$key} = "$nosave";

$key = q/user_prepare_to_check_cuts/;
$ref_files{$key} = "$dir".q|node81.html|; 
$noresave{$key} = "$nosave";

$key = q/user_receive_cp_data/;
$ref_files{$key} = "$dir".q|node78.html|; 
$noresave{$key} = "$nosave";

$key = q/user_select_candidates/;
$ref_files{$key} = "$dir".q|node58.html|; 
$noresave{$key} = "$nosave";

$key = q/output/;
$ref_files{$key} = "$dir".q|node25.html|; 
$noresave{$key} = "$nosave";

$key = q/user_start_heurs/;
$ref_files{$key} = "$dir".q|node35.html|; 
$noresave{$key} = "$nosave";

$key = q/PVM/;
$ref_files{$key} = "$dir".q|node12.html|; 
$noresave{$key} = "$nosave";

$key = q/diving_strategy/;
$ref_files{$key} = "$dir".q|node93.html|; 
$noresave{$key} = "$nosave";

$key = q/IGD/;
$ref_files{$key} = "$dir".q|node23.html|; 
$noresave{$key} = "$nosave";

$key = q/user_receive_lp_data/;
$ref_files{$key} = "$dir".q|node50.html|; 
$noresave{$key} = "$nosave";

$key = q/communication/;
$ref_files{$key} = "$dir".q|node11.html|; 
$noresave{$key} = "$nosave";

$key = q/resources/;
$ref_files{$key} = "$dir".q|node26.html|; 
$noresave{$key} = "$nosave";

$key = q/user_unpack_cuts/;
$ref_files{$key} = "$dir".q|node64.html|; 
$noresave{$key} = "$nosave";

$key = q/user_send_lp_solution/;
$ref_files{$key} = "$dir".q|node65.html|; 
$noresave{$key} = "$nosave";

1;

