#ifndef _SPP_CG_CLIQUE_H
#define _SPP_CG_CLIQUE_H

int find_violated_star_cliques PROTO((spp_cg_problem *spp, double etol));
void spp_delete_node PROTO((spp_cg_problem *spp, int del_ind,
			    int *pcurrent_nodenum, int *current_indices,
			    int *current_degrees, double *current_values));
int choose_next_node PROTO((spp_cg_problem *spp, int current_nodenum,
			    int *current_indices, int *current_degrees,
			    double *current_values));
int find_violated_row_cliques PROTO((spp_cg_problem *spp, double etol));
int enumerate_maximal_cliques PROTO((spp_cg_problem *spp, int pos, double etol));
int greedy_maximal_clique PROTO((spp_cg_problem *spp, cut_data *new_cut,
				 int length, int *indices, int pos, double etol));
#endif
