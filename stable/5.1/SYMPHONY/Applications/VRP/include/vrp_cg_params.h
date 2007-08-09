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

#ifndef _VRP_CG_PARAMS_H
#define _VRP_CG_PARAMS_H

/* which_connected_routine choices */
#define CONNECTED    0
#define BICONNECTED  1
#define BOTH         2

typedef struct VRP_CG_PARAMS{
   int            verbosity;
   char           tsp_prob;
   int            do_greedy;
   int            greedy_num_trials;
   int            do_extra_in_root;
   int            which_tsp_cuts;
   int            which_connected_routine;
   int            max_num_cuts_in_shrink;
}vrp_cg_params;

#endif
