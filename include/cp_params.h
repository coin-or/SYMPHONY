/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _CUT_POOL_PARAMS_H
#define _CUT_POOL_PARAMS_H

/*===========================================================================*\
 * Contains the parameters necessary for the functioning of the cut pool
\*===========================================================================*/

typedef struct CP_PARAMS{
   int     verbosity;
   int     warm_start;
   char    warm_start_file_name[MAX_FILE_NAME_LENGTH +1];
   int     logging;
   char    log_file_name[MAX_FILE_NAME_LENGTH +1];
   int     block_size;
   int     max_size;
   int     max_number_of_cuts;
   int     cuts_to_check;
   int     delete_which;
   int     touches_until_deletion;
   int     min_to_delete;
   int     check_which;
}cp_params;

#endif
