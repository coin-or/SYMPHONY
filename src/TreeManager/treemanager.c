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

#ifndef COMPILE_IN_TM

#define COMPILING_FOR_TM

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tm.h"
#include "proccomm.h"
#include "timemeas.h"

/*===========================================================================*/

/*===========================================================================*\
 * This is the main() that is used if the TM is running as a separate        
 * process. This file is only used in that case.                             
\*===========================================================================*/

int main(void)
{
   tm_prob *tm;
   int termcode;
   
   /* set stdout to be line buffered */
   setvbuf(stdout, (char *)NULL, _IOLBF, 0);
   
   register_process();  /*Enroll this program in PVM*/

   tm = tm_initialize(NULL, NULL, 0);
   
   tm->start_time = wall_clock(NULL);
   
   termcode = tm_close(tm, solve(tm));

   comm_exit();

   return(termcode);
}

#endif
