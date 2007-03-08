/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _PREP_PARAMS_H
#define _PREP_PARAMS_H

/*---------------------------------------------------------------------------*\
| The list of parameters associated with pre-processing                       |
|                                                                             |
\*---------------------------------------------------------------------------*/
typedef struct PREP_PARAMS{
   int               do_prep;
   int               do_probe;
   int               prep_verbosity;
   int               probe_verbosity;
   int               probe_level;
   int               display_stats;
}prep_params;

#endif
