/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2004 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _PACK_ARRAY_H
#define _PACK_ARRAY_H

#include "proto.h"
#include "BB_types.h"

void pack_array_desc PROTO((array_desc *adesc));
array_desc *unpack_array_desc PROTO((array_desc *padesc));
void pack_double_array_desc PROTO((double_array_desc *dad,
				   char explicit_packing));
void unpack_double_array_desc PROTO((double_array_desc *dad,
				     char explicit_packing));
void pack_basis PROTO((basis_desc *basis, char explicit_packing));
basis_desc *unpack_basis PROTO((basis_desc *pbasis, char explicit_packing));

#endif
