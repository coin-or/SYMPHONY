#ifndef _SPP_LP_H_
#define _SPP_LP_H_

#include "proto.h"
#include "BB_types.h"
#include "dg_params.h"
#include "lp.h"

#include "spp_types.h"
#include "spp_lp_params.h"

typedef struct SPP_LP_TMP {
   char               *ctmp_2nD;       /* length: 2*n*DSIZE */
   double             *dtmp_m;         /* length rownum */
   double             *dtmp_n;
   int                *itmp_m;
   int                *itmp_2n;
}spp_lp_tmp;


typedef struct SPP_LP_PROBLEM {
   spp_lp_params      *par;
   spp_lp_tmp         *tmp;
   col_ordered        *cmatrix;
   char                wname[MAX_NAME_LENGTH +1];
                          /* name of window in which frac solns are dispd */
}spp_lp_problem;

#endif
