/* glprng.h */

/*----------------------------------------------------------------------
-- Copyright (C) 2000, 2001, 2002, 2003, 2004 Andrew Makhorin,
-- Department for Applied Informatics, Moscow Aviation Institute,
-- Moscow, Russia. All rights reserved. E-mail: <mao@mai2.rcnet.ru>.
--
-- This file is part of GLPK (GNU Linear Programming Kit).
--
-- GLPK is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2, or (at your option)
-- any later version.
--
-- GLPK is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
-- or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
-- License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with GLPK; see the file COPYING. If not, write to the Free
-- Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
-- 02111-1307, USA.
----------------------------------------------------------------------*/

#ifndef _GLPRNG_H
#define _GLPRNG_H

#define rng_create_rand       glp_rng_create_rand
#define rng_init_rand         glp_rng_init_rand
#define rng_next_rand         glp_rng_next_rand
#define rng_unif_rand         glp_rng_unif_rand
#define rng_delete_rand       glp_rng_delete_rand

typedef struct RNG RNG;

struct RNG
{     /* Knuth's portable pseudo-random number generator */
      int A[56];
      /* pseudo-random values */
      int *fptr;
      /* the next A value to be exported */
};

RNG *rng_create_rand(void);
/* create pseudo-random number generator */

void rng_init_rand(RNG *rand, int seed);
/* initialize pseudo-random number generator */

int rng_next_rand(RNG *rand);
/* obtain pseudo-random integer in [0, 2^31-1] */

int rng_unif_rand(RNG *rand, int m);
/* obtain pseudo-random integer in [0, m-1] */

void rng_delete_rand(RNG *rand);
/* delete pseudo-random number generator */

#endif

/* eof */
