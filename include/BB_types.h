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

#ifndef _BB_TYPES_H
#define _BB_TYPES_H

#define MAX_CHILDREN_NUM 4

/*===========================================================================*\
 * This file contains those type definitions which are used in more than one
 * process of the black box model.
\*===========================================================================*/

#ifndef WIN32
#include <unistd.h>            /* this defines sleep() */
#endif

#include "proto.h"

typedef struct LP_SOL{
   int            lp;          /* the tid of the lp process asssociated with
				  the current solution */
   int            xlength;     /* the number of nonzeros in the lp solution
				  currently being processed*/
   int            xlevel;      /* the level at which the current solution was
				  generated*/
   int            xindex;
   int            xiter_num;
   int            max_sol_length;
   int           *xind;        /* the indices of the nonzeros in the current
				  solution*/
   double        *xval;        /* the values of the nonzeros in the current
				  solution*/
   double         objval;      /* the objective function value of the current
				  relaxation*/
   double         lpetol;
}lp_sol;

typedef struct BASE_DESC{
   int            varnum;
   int           *userind;
   double        *lb;          /* even if there are global lb and ub, we */
   double        *ub;          /* fill these arrays out */
   int            cutnum;
}base_desc;

/*===========================================================================*\
 * Stores the data corresponding to a particular cut
\*===========================================================================*/

typedef struct CUT_DATA{
   int            size;        /* the size of the coef array */
   char          *coef;        /* an array which contains the data necessary to
				  construct the cut -- it is stored in a
				  packed form. The types of the cut tells how
				  to "unpack" it */
   double         rhs;         /* the right hand side for the constraint*/
   double         range;   
   char           type;        /* the type of cut */
   char           sense;       /* sense of the cut constraint */
   char           deletable;   /* whether or not this cut should be removed
				  from the LP after being added */
   char           branch;      /* shows whether we can branch on its cut if
				  this row becomes slack */
   int            name;        /* internal to the BB. The identifier of the
				  cut. >=0 if exists, -1 if does not exist yet,
				  but the cuT is sent to the cutpool, -2 if
				  no name & no pool */
}cut_data;

typedef struct CONSTRAINT{
   cut_data      *cut;
   int            ineff_cnt;
   int            eff_cnt;
   char           free;
   char           deletable;
}constraint;

typedef struct WAITING_ROW{
   int            source_pid;
   cut_data      *cut;
   int           *matind;
   double        *matval;
   int            nzcnt;
   double         violation;
}waiting_row;

/*===========================================================================*\
 * The following three definitions are used to describe the search tree
 * nodes.
\*===========================================================================*/

typedef struct VAR_DESC{
   int            userind;
   int            colind;
   double         lb;
   double         ub;
}var_desc;

typedef struct BRANCH_DESC{
   int            name;        /* the userind/cut name depending on the type */
   char           type;        /* description of the child. All of them are
				  natural for branching cuts. For branching
				  variables they should be interpreted as if
				  we were adding a cut with a single variable
				  on the left hand side */
   char           sense;
   double         rhs;
   double         range;
   int            branch;
}branch_desc;

typedef struct ARRAY_DESC{
   char           type;        /* NO_DATA_STORED, EXPLICIT_LIST, WRT_PARENT */
   int            size;
   int            added; 
   int           *list;
}array_desc;

typedef struct DOUBLE_ARRAY_DESC{
   char           type;        /* NO_DATA_STORED, EXPLICIT_LIST, WRT_PARENT */
   int            size;        /* the size of list, stat */
   int           *list;
   int           *stat;
}double_array_desc;

typedef struct BASIS_DESC{
   char                basis_exists;
   
   /*========================================================================*\    * Notes:
    *  1) for base...:
    *     if list is non-NULL then it refers to col/row inds, not userinds
    *	    or cut names.
    *  2) the stat field of extra... ponts into the stat field of base...
    *  3) if extra... is EXPLICIT_LIST then the node_desc structure's
    *        cutind/uind fields should be used.
    *  4) EXPLICIT_LIST in uind implies that extravars is explicit,
    *     EXPLICIT_LIST in cutind implies that extrarows is explicit.
   \*========================================================================*/
   double_array_desc   basevars;
   double_array_desc   extravars;
   double_array_desc   baserows;
   double_array_desc   extrarows;
}basis_desc;

typedef struct NODE_DESC{
   /*========================================================================*\
    * The userindices of variables in this node (but not for the base
    * variables); The basis header for this node; The not-yet-permanently-fixed
    * variables (again, no base variable is listed here); and the cuts at
    * this node
   \*========================================================================*/
   array_desc     uind;
   basis_desc     basis;
   array_desc     not_fixed;
   int            nf_status;   /* NF_CHECK_ALL, NF_CHECK_AFTER_LAST,
				  NF_CHECK_UNTIL_LAST, NF_CHECK_NOTHING */
   array_desc     cutind;
#if defined(COMPILING_FOR_LP) || defined(COMPILING_FOR_MASTER) || defined(COMPILE_IN_LP)
   cut_data     **cuts;        /* this is not used in TM anyway. */
#endif

   /* Any additional info the user might want to pass */
   int           desc_size;
   char         *desc;
}node_desc;

typedef struct BRANCH_OBJ{
   char          type;         /* Type of the candidate */
#if defined(COMPILING_FOR_LP) || defined(COMPILE_IN_LP) 
   int           position;     /* The position of the candidate */
   waiting_row  *row;          /* Description of the left hand side; makes
				  sense only for branching cuts */
#endif
   int           child_num;    /* Number of kids */
#if defined(COMPILING_FOR_TM) || defined(COMPILE_IN_LP) 
   int           name;         /* userind for VAR, the index for CUT */
#endif

   /*========================================================================*\
    * Description of the children.
    * All of them are natural for branching cuts.
    * For branching variables they should be interpreted as if we were adding
    * a cut with a single variable on the left hand side
   \*========================================================================*/

#ifdef MAX_CHILDREN_NUM
   char          sense[MAX_CHILDREN_NUM];
   double        rhs[MAX_CHILDREN_NUM];
   double        range[MAX_CHILDREN_NUM];
   int           branch[MAX_CHILDREN_NUM];
#ifdef COMPILE_FRAC_BRANCHING
   int           frac_num[MAX_CHILDREN_NUM];
   int          *frac_ind[MAX_CHILDREN_NUM];
   double       *frac_val[MAX_CHILDREN_NUM];
#endif
#else
   char         *sense;
   double       *rhs;
   double       *range;
   int          *branch;
#ifdef COMPILE_FRAC_BRANCHING
   int          *frac_num;
   int         **frac_ind;
   double      **frac_val;
#endif
#endif

#if defined(COMPILING_FOR_LP) || defined(COMPILE_IN_LP) 
   double        lhs;          /* purely for the user */

#ifdef MAX_CHILDREN_NUM
   double        objval[MAX_CHILDREN_NUM];   /* arrays of size 'number' */
   int           termcode[MAX_CHILDREN_NUM];
   int           iterd[MAX_CHILDREN_NUM];
   int           feasible[MAX_CHILDREN_NUM];
#else
   double       *objval;   /* arrays of size 'number' */
   int          *termcode;
   int          *iterd;
   int          *feasible;
#endif

#endif
}branch_obj;

/*===========================================================================*/

typedef struct STR_INT{
#ifdef _OPENMP
   char      *str;
#else
   char       str[MAX_LINE_LENGTH +1];
#endif
   int        code;
}str_int;

/*===========================================================================*\
 * This is the time measurement structure for an LP node
\*===========================================================================*/

typedef struct NODE_TIMES{
   double        communication;
   double        lp;
   double        separation;
   double        fixing;
   double        pricing;
   double        strong_branching;
   double        wall_clock_lp;
   double        ramp_up_tm;
   double        ramp_up_lp;
   double        ramp_down_time;
   double        idle_diving;
   double        idle_node;
   double        idle_names;
   double        idle_cuts;
   double        start_node;
   double        cut_pool;
}node_times;

#endif
