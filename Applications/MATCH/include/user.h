/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000, 2001, 2002 Ted Ralphs. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _USER_H
#define _USER_H

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

#endif
