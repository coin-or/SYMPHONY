/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2005 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef DG_U_H
#define DG_U_H

#include <stdio.h>

#include "proto.h"

struct WINDOW;

/*===========================================================================*/
/*====================== User supplied functions ============================*/
/*===========================================================================*/

int user_dg_process_message PROTO((void *user, struct WINDOW *win,
				   FILE *write_to));
int user_dg_init_window PROTO((void **user, struct WINDOW *win));
int user_dg_free_window PROTO((void **user, struct WINDOW *win));
int user_initialize_dg PROTO((void **user));
int user_free_dg PROTO((void **user));
int user_interpret_text PROTO((void *user, int text_length, char *text,
			       int owner_tid));

#endif


