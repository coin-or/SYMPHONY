/*===========================================================================*/
/*                                                                           */
/* This file is part of a demonstration application for use with the         */
/* SYMPHONY Branch, Cut, and Price Library. This application is a solver for */
/* Set Partitioning Problems.                                                */
/*                                                                           */
/* (c) Copyright 2003 Marta Eso and Ted Ralphs. All Rights Reserved.         */
/*                                                                           */
/* This application was originally developed by Marta Eso and was modified   */
/* Ted Ralphs (tkralphs@lehigh.edu)                                          */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _SPP_COMMON_H_
#define _SPP_COMMON_H_

#include <stdio.h>

#include "proto.h"
#include "spp_types.h"

int *spp_lexsort PROTO((col_ordered *matrix));
int spp_lexcompare PROTO((const void *i, const void *j));
void spp_free_cmatrix PROTO((col_ordered *m));
void spp_free_rmatrix PROTO((row_ordered *m));
void spp_column_to_row PROTO((col_ordered *cm, row_ordered *rm, int *i_tmpm,
			      int **istar_tmpm));
void spp_row_to_column PROTO((row_ordered *rm, col_ordered *cm, int *i_tmpn,
			      row_ind_type **rowindstar_tmpn));
FILE *get_filehandler PROTO((int counter, char *filename, char *fname));

#endif
