/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* the Matching Problem.                                                     */
/*                                                                           */
/* (c) Copyright 2003 Michael Trick and Ted Ralphs. All Rights Reserved.     */
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
/* Cut types */

#define TRIANGLE 1

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the value of any run-time parameters.
\*---------------------------------------------------------------------------*/

typedef struct USER_PARAMETERS{
   /* Name of file containing the instance data */
   char             infile[MAX_FILE_NAME_LENGTH + 1];
}user_parameters;

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the instance data after it is read in.
\*---------------------------------------------------------------------------*/

typedef struct USER_PROBLEM{
   int              colnum;         /* Number of rows in base matrix */
   int              rownum;         /* Number of columns in base matrix */
   user_parameters  par;            /* Parameters */
   int		    nnodes;         /* Number of nodes */
   int		    cost[200][200]; /* Cost of assigning i to j */ 
   int		    node1[20000];   /* First index of each variable */
   int		    node2[20000];   /* Second index of each variable */
}user_problem;


int match_read_data PROTO((void *user, char *infile));
int match_load_problem PROTO((sym_environment *env, void *user));

#endif
