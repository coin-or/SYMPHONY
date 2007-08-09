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

#ifndef _VRP_MACROS_H
#define _VRP_MACROS_H

/* SYMPHONY include files */
#include "sym_proto.h"

#define CHECK_DEBUG_PAR(x, y) \
if (x != 0 && x != 4) {                                               \
   (void) fprintf(stderr, "\nio: illegal debug parameter %s\n\n", y); \
   exit(1);                                                           \
}

#define READ_FLOAT_PAR(par)						\
if (sscanf(value, "%f", &(par)) != 1){					\
   (void) fprintf(stderr, "\nio: error reading parameter %s\n\n", key);	\
   exit(1);								\
}

void BOTH_ENDS PROTO((int index, int *vh, int *vl));
int NEAREST_INT PROTO((double num));
int INDEX PROTO((int v0, int v1));
int BINS PROTO((int weight, int capacity));
int RHS PROTO((int cust_num, int weight, int capacity));

#endif
