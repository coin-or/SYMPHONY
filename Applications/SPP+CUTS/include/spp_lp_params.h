#ifndef _SPP_LP_PARAMS_H
#define _SPP_LP_PARAMS_H

typedef struct SPP_LP_PARAMS{
   double           granularity;      /* inherited from air_params */
   int              do_lift_in_lp;
   double           lp_dj_threshold_frac;
   double           lp_dj_threshold_abs;
   int              lanmax;
   int              which_atilde;
}spp_lp_params;

#endif
