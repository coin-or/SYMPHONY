/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2013 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _USER_H
#define _USER_H

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the value of any run-time parameters.
\*---------------------------------------------------------------------------*/

typedef struct USER_PARAMETERS{
   /* Name of file containingthe instance data */
   char             infile[MAX_FILE_NAME_LENGTH + 1];
   int              test;
   char             test_dir[MAX_FILE_NAME_LENGTH +1]; /* Test files directory */
}user_parameters;

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the instance data after it is read in.
\*---------------------------------------------------------------------------*/

typedef struct USER_PROBLEM{
   user_parameters  par;    /* Parameters */
   int              colnum; /* Number of columns in base matrix */
   int              rownum; /* Number of rows in base matrix */

   /* Various parameters to store data read from mps file.
    * Note that data read from mps file is not the final problem that
    * is loaded into SYMPHONY.
    * These parameters are used in creating final problem data */
   int             *matbeg_row; /* */
   int             *matind_row; /* */
   double          *matval_row; /* */
   int              con_sense_e;/* */
   int              con_sense_g;/* */
   int              con_sense_l;/* */
   int              ubinfty;    /* number of infinity UBs */
   int              lbinfty;    /* number of -infinity LBs */
   double          *tempub;     /* condensed finite UBs */
   double          *templb;     /* condensed finite LBs */
   int             *infubsofar; /* number of infinite UBs so far */
   int             *inflbsofar; /* number of infinite LBs so far */
   int             *infubind;   /* indicator of a infinite UB */
   int             *inflbind;   /* indicator of a -infinite LB */
   double           infty;      /* */
   int             *ccind;      /* Index of complementarity constraint */
   MIPdesc         *mip;        /* MIP data structure to store problem data */

   int              feasible;   /* Feasibility status of current LP relaxation */
   int             *vvind;      /* Var indices of violated complementarity constraints */
   int              vvnum;      /* Number of violated complementarity cons. Size of vvind */

   double          *rowact;     /* row activities of constraints of current LP relaxation */

}user_problem;

int user_read_data PROTO((user_problem *prob, char *infile));
int user_load_problem PROTO((sym_environment *env, user_problem *prob));

#endif
