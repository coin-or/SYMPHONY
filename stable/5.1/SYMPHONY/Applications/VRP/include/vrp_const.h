/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _VRP_CONST_H
#define _VRP_CONST_H

#define LENGTH 255
#define KEY_NUM 43
#define DEAD 2
#define NEAR_INS -1
#define FAR_INS -2
#define DEPOT_PENALTY 20
#define RRR 6378.388
#define MY_PI 3.141592
#define LINE_LEN 80

/*---------------- distance types -------------------------------------------*/
#define _EXPLICIT 0
#define _EUC_2D   1
#define _EUC_3D   2
#define _MAX_2D   3
#define _MAX_3D   4
#define _MAN_2D   5
#define _MAN_3D   6
#define _CEIL_2D  7
#define _GEO      8
#define _ATT      9

/*---------------- message types --------------------------------------------*/
#define VRP_LB_DATA                1
#define VRP_LB_DATA2               2
#define VRP_BROADCAST_DATA         3
#define EXCHANGE_HEUR_TOUR         4
#define ROUTE_FINS_START_RULE      5
#define ROUTE_NINS_START_RULE      6
#define ROUTE_FNINS_START_RULE     7
#define FINI_RATIO                 8
#define TSP_FINI_RATIO             9
#define ROUTE_FINS_VRP_DATA        10
#define ROUTE_NINS_VRP_DATA        11
#define ROUTE_FNINS_VRP_DATA       12
#define SWEEP_TRIALS               13
#define TSP_NI_TRIALS              14
#define TSP_FI_TRIALS              15
#define TSP_FINI_TRIALS            16
#define S3_NUMROUTES               17
#define NC_NUMROUTES               18
#define TSP_START_POINT            19
#define SAVINGS_DATA               20
#define SAVINGS2_DATA              21
#define SAVINGS3_DATA              22
#define DISPLAY_DATA               23
#define STOP                       24


/*---------------- cut types ------------------------------------------------*/
#define SUBTOUR_ELIM_SIDE    0
#define SUBTOUR_ELIM_ACROSS  1
#define SUBTOUR_ELIM         2
#define CLIQUE               3

/*---------------- tsp cut routines -----------------------------------------*/

#define NO_TSP_CUTS    0
#define SUBTOUR        1
#define BLOSSOM        2
#define COMB           4
#define ALL_TSP_CUTS   7

#define NUM_RANDS 6

#define ACTIVE_NODE_LIST_BLOCK_SIZE 100
#define DELETE_POWER 3
#define DELETE_AND 0x07

/*-------------- base variable selection rules ------------------------------*/
#define EVERYTHING_IS_EXTRA 0
#define SOME_ARE_BASE       1
#define EVERYTHING_IS_BASE  2

/*--------- constants used in creating the edges lists for the root ---------*/
#define CHEAP_EDGES      0
#define REMAINING_EDGES  1

/*--------- constants for saving the small graph ----------------------------*/
#define SAVE_SMALL_GRAPH 1
#define LOAD_SMALL_GRAPH 2

/*--------- constants for defining which set of exchange heuristics to do --*/     
#define FIRST_SET        1
#define SECOND_SET       2

#endif
