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
   /* Name of file containing the instance data */
   char             infile[MAX_FILE_NAME_LENGTH + 1];
   /* Name of the file containing auxiliary data related to instance */
   char             auxfile[MAX_FILE_NAME_LENGTH + 1];
   /* Indicator representing if infile contains a bilevel instance (1) or a single
      level instance (0) */
   // TODO: Suresh - this is kind of redundant with auxfile
   bool             bilevel;
   int              test;
   char             test_dir[MAX_FILE_NAME_LENGTH +1]; /* Test files directory */
}user_parameters;

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the auxiliary data after it is read in.
\*---------------------------------------------------------------------------*/

typedef struct AUX_DATA{
   int         num_lower_vars;       /* Number of lower level variables */
   int         num_lower_cons;       /* Number of lower level constraints */
   int        *index_lower_vars;     /* Indices of lower level variables */
   int        *index_lower_cons;     /* Indices of lower level constraints */
   double     *coeff_lower_obj;      /* Coefficients of lower level obj fn. */
   int         sense_lower_obj;      /* Lower level objective sense: +1 = min, -1 = max */
   double     *coeff_upper_cons;     /* Coefficients of upper level constraint
                                        for interdiction problems */
   double      rhs_upper_cons;       /* RHS of upper level constraint for
                                        interdiction problems */
}aux_data;

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
   int             origvar_num; /* original # of vars */
   double      *origobj_coeffs; /* original obj fn coeffs */
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
   int             ccnum;       /* number of complementarity conditions */
   MIPdesc         *mip;        /* MIP data structure to store problem data */

   int              feasible;   /* Feasibility status of current LP relaxation */
   int             *vvind;      /* Var indices of violated complementarity constraints */
   int              vvnum;      /* Number of violated complementarity cons. Size of vvind */

   double          *rowact;     /* row activities of constraints of current LP relaxation */

#if defined(CHECK_CUT_VALIDITY) || defined(TRACE_PATH)
   int             feas_sol_size;
   double            *feas_sol;
#endif
   aux_data        aux;         /* Auxiliary data for bilevel instance */

}user_problem;

int user_read_data PROTO((user_problem *prob, char *infile));
int user_load_problem PROTO((sym_environment *env, user_problem *prob));
int user_read_aux_data PROTO((user_problem *prob, char *infile));
int user_load_bilevel_problem PROTO((sym_environment *env, user_problem *prob));

#endif
