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

#ifndef _COMPUTE_COST_H
#define _COMPUTE_COST_H

#include "proto.h"
#include "vrp_common_types.h"

int compute_icost PROTO((distance *dist, int v0, int v1));
void canonical_tour PROTO((distance *dist, best_tours *cur_tour,
			   int vertnum, int capacity, int *demand));
int route_calc PROTO((distance *dist, _node *tour, int numroutes, 
		      route_data *route_info, int *demand));
int compute_tour_cost PROTO((distance *dist, _node *tour));
double ECOST PROTO((double *cost, int v0, int v1, int vertnum));
int ICOST PROTO((distance *dist, int v0, int v1));
int MCOST PROTO((distance *dist, int v0, int v1, int *lamda));
int TCOST PROTO((distance *dist, int v0, int v1, int *lamda, int mu));

#endif
