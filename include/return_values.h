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

#ifndef _RETURN_VALUES_H
#define _RETURN_VALUES_H

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  Return Values                       **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*----------------------- Global return codes -------------------------------*/
#define FUNCTION_TERMINATED_NORMALLY    0
#define ERROR__USER                    -100

/*------------------- Return codes for sym_solve() --------------------------*/
#define ERROR__NO_BRANCHING_CANDIDATE  -101
#define ERROR__ILLEGAL_RETURN_CODE     -102
#define ERROR__NUMERICAL_INSTABILITY   -103
#define ERROR__ILLEGAL_BRANCHING       -104
#define ERROR__COMM_ERROR              -105

/*-------------- Return codes for sym_parse_comand_line() -------------------*/
#define ERROR__OPENING_PARAM_FILE      -110
#define ERROR__PARSING_PARAM_FILE      -111

/*----------------- Return codes for sym_load_problem() ---------------------*/
#define ERROR__READING_GMPL_FILE       -120
#define ERROR__READING_WARM_START_FILE -121

#endif
