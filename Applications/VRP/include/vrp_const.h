/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Vehicle Routing Problem and the Traveling Salesman Problem.           */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2001 Ted Ralphs. All Rights Reserved.           */
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
#define KEY_NUM 41
#define DEAD 2
#define NEAR_INS -1
#define FAR_INS -2
#define DEPOT_PENALTY 20
#define RRR 6378.388
#define MY_PI 3.141592
#define LINE_LEN 80
#ifndef M_PI
#define M_PI 3.141592
#endif

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
#define VRP_DATA         1
#define DISPLAY_DATA     2
#define NUMROUTES        3
#define COORD_DATA       4
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#define HEUR_START_POINT 5
#define HEUR_TOUR        6
#define SWEEP_TRIALS     7
#define SAVINGS_DATA     8
#define TSP_TRIALS       9
#define LOWER_BOUND      10

/*--------------- algorithms ------------------------------------------------*/
#define SWEEP         0
#define SAVINGS       1
#define SAVINGS3      2
#define NEAR_CLUSTER  3
#define TSP_NI        4
#define TSP_FI        5
#define TSP_FINI      6

#define IN_TOUR      -1
#define IN_TREE      -1
#define NOT_NEIGHBOR  0
/*___END_EXPERIMENTAL_SECTION___*/

/*---------------- cut types ------------------------------------------------*/
#define SUBTOUR_ELIM_SIDE    0
#define SUBTOUR_ELIM_ACROSS  1
#define SUBTOUR_ELIM         2
#define CLIQUE               3
/*__BEGIN_EXPERIMENTAL_SECTION__*/
#define FARKAS               4
#define NO_COLUMNS           5
#define GENERAL_NONZEROS     6
/*___END_EXPERIMENTAL_SECTION___*/

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
