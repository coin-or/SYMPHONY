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

/*---------------------------------------------------------------------------*\
 * 
\*---------------------------------------------------------------------------*/

typedef struct USER_PARAMETERS{
   char             infile[MAX_FILE_NAME_LENGTH + 1];
}user_parameters;

/*---------------------------------------------------------------------------*\
 * 
\*---------------------------------------------------------------------------*/

typedef struct USER_PROBLEM{
   int              colnum;
   int              rownum;
   user_parameters  par;
}user_problem;

#endif
