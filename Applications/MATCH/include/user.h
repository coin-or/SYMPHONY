/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Matching Problem.                                                     */
/*                                                                           */
/* (c) Copyright 2004 Michael Trick and Ted Ralphs. All Rights Reserved.     */
/*                                                                           */
/* This application was originally written by Michael Trick and was modified */
/* by Ted Ralphs (tkralphs@lehigh.edu).                                      */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _USER_H
#define _USER_H

#include "master.h"
#include "BB_macros.h"

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the instance data after it is read in.
\*---------------------------------------------------------------------------*/

typedef struct USER_PROBLEM{
   int              colnum;         /* Number of rows in base matrix */
   int              rownum;         /* Number of columns in base matrix */
   int		    nnodes;         /* Number of nodes */
   int		    cost[200][200]; /* Cost of assigning i to j */ 
   int		    node1[20000];   /* node1[i] is the first component of
				       the assignment with index 'i' */
   int		    node2[20000];   /* node2[i] is the second component of
				       the assignment with index 'i' */
   int              index[200][200];/* index[j][k] is the index of the variable
				       associated with assigning 'j' to 'k'*/
}user_problem;


int match_read_data PROTO((sym_environment *env, void *user, char *infile));
int match_load_problem PROTO((sym_environment *env, void *user));

#endif
