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

#ifndef _MPP_TYPES_H
#define _MPP_TYPES_H

#include "proto.h"

/*---------------------------------------------------------------------------*\
 * The problem data structure contains the data for a problem instance, as
 * well as some of the tours that have been generated.
\*---------------------------------------------------------------------------*/

typedef struct MPP_PROBLEM{
   char            name[100];  /* the name of the problem instance */
   char infile[12];
 //  int             dg_id;     /* drawgraph process id */
   
   int nodes;   /* the number of nodes in the problem*/
   int edges;   /* number of edges in the problem */
   int arcs;    /* number of arcs in the problem */
   int *cost;   /* an array containing the cost for each edge*/
   int *start_node;    /* an array containing the start for each edge*/
   int *end_node;    /* an array containing the end for each edge*/
   char *types;  /* an array containing the types for each edge */
   int *is_odd;  /* arrat containing for each node a 1 if odd 0 if not*/
   int odd_checker; /*indicates if odds have been checked */

}mpp_problem; 

#endif
