/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2013 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _USER_H
#define _USER_H

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the value of any run-time parameters.
\*---------------------------------------------------------------------------*/

typedef struct USER_PARAMETERS {
  /* Name of file containing the instance data */
  char             infile[MAX_FILE_NAME_LENGTH + 1];
} user_parameters;

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the instance data after it is read in.
\*---------------------------------------------------------------------------*/

typedef struct USER_PROBLEM {
  int              colnum; /* Number of columns in base matrix */
  int              rownum; /* Number of rows in base matrix */
  user_parameters  par;    /* Parameters */
  // number of cones
  int num_cones;
  // type of cone 0 for LORENTZ cone and 1 for ROTATED LORENTZ cone
  int * cone_type;
  int * cone_size;
  int ** cone_members;
  char * is_int;
  double * curr_solution;
} user_problem;

#endif
