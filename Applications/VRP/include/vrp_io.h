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

#ifndef _VRP_IO_H
#define _VRP_IO_H

/* SYMPHONY include files */
#include "proto.h"

/* VRP include files */
#include "vrp_types.h"

void vrp_readparams PROTO((vrp_problem *vrp, char *filename, int argc,
		       char **argv));
void vrp_io PROTO((vrp_problem *vrp, char *infile));

#endif
