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
   int               prep_level;
   int               do_prep;
   int               do_probe;
   int               prep_verbosity;
   int               probe_verbosity;
   int               probe_level;
   int               display_stats;
   double            etol; 
   int               keep_row_ordered; 
   char              do_single_row_rlx; 
   double            single_row_rlx_ratio;
   int               max_sr_cnt;
   char              do_aggregate_row_rlx; 
   double            max_aggr_row_ratio;   
   int               max_aggr_row_cnt;
   int               iteration_limit; 
}prep_params;

#endif
