/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* Capacitated Network Routing Problems.                                     */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This application was developed by Ted Ralphs (tkralphs@lehigh.edu)        */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _CNRP_IO_H
#define _CNRP_IO_H

/* SYMPHONY include files */
#include "proto.h"

/* CNRP include files */
#include "cnrp_types.h"

void vrp_readparams PROTO((vrp_problem *vrp, char *filename, int argc,
		       char **argv));
void vrp_io PROTO((vrp_problem *vrp, char *infile));

#endif
