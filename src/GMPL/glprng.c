/* glprng.c */

/*----------------------------------------------------------------------
-- This file is part of GLPK (GNU Linear Programming Kit).
--
-- This module is a modified version of the module GB_FLIP, a portable
-- pseudo-random number generator. The original version of GB_FLIP is a
-- part of The Stanford GraphBase developed by Donald E. Knuth (see
-- http://www-cs-staff.stanford.edu/~knuth/sgb.html).
--
-- Note that all changes concern only external names, so this modified
-- version provides exactly the same results as the original version.
--
-- Changes were made by Andrew Makhorin <mao@mai2.rcnet.ru>.
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

#include "glplib.h"
#include "glprng.h"

#if 0
int A[56] = { -1 };
#else
#define A (rand->A)
#endif
/* pseudo-random values */

#if 0
int *fptr = A;
#else
#define fptr (rand->fptr)
#endif
/* the next A value to be exported */

#define mod_diff(x, y) (((x) - (y)) & 0x7FFFFFFF)
/* difference modulo 2^31 */

static int flip_cycle(RNG *rand)
{     /* this is an auxiliary routine to do 55 more steps of the basic
         recurrence, at high speed, and to reset fptr */
      int *ii, *jj;
      for (ii = &A[1], jj = &A[32]; jj <= &A[55]; ii++, jj++)
         *ii = mod_diff(*ii, *jj);
      for (jj = &A[1]; ii <= &A[55]; ii++, jj++)
         *ii = mod_diff(*ii, *jj);
      fptr = &A[54];
      return A[55];
}

/*----------------------------------------------------------------------
-- rng_create_rand - create pseudo-random number generator.
--
-- *Synopsis*
--
-- #include "glprng.h"
-- RNG *rng_create_rand(void);
--
-- *Description*
--
-- The routine rng_create_rand creates and initializes a pseudo-random
-- number generator.
--
-- *Returns*
--
-- The routine returns a pointer to the generator created. */

RNG *rng_create_rand(void)
{     RNG *rand;
      int i;
      rand = umalloc(sizeof(RNG));
      A[0] = -1;
      for (i = 1; i <= 55; i++) A[i] = 0;
      fptr = A;
      rng_init_rand(rand, 0);
      return rand;
}

/*----------------------------------------------------------------------
-- rng_init_rand - initialize pseudo-random number generator.
--
-- *Synopsis*
--
-- #include "glprng.h"
-- void rng_init_rand(RNG *rand, int seed);
--
-- *Description*
--
-- The routine rng_init_rand initializes the pseudo-random number
-- generator. The parameter seed may be any integer number. Note that
-- on creating the generator this routine is called with the parameter
-- seed equal to zero. */

void rng_init_rand(RNG *rand, int seed)
{     int i;
      int prev = seed, next = 1;
      seed = prev = mod_diff(prev, 0);
      A[55] = prev;
      for (i = 21; i; i = (i + 21) % 55)
      {  A[i] = next;
         next = mod_diff(prev, next);
         if (seed & 1)
            seed = 0x40000000 + (seed >> 1);
         else
            seed >>= 1;
         next = mod_diff(next, seed);
         prev = A[i];
      }
      flip_cycle(rand);
      flip_cycle(rand);
      flip_cycle(rand);
      flip_cycle(rand);
      flip_cycle(rand);
      return;
}

/*----------------------------------------------------------------------
-- rng_next_rand - obtain pseudo-random integer in [0, 2^31-1].
--
-- *Synopsis*
--
-- #include "glprng.h"
-- int rng_next_rand(RNG *rand);
--
-- *Returns*
--
-- The routine rng_next_rand returns a next pseudo-random integer which
-- is uniformly distributed between 0 and 2^31-1, inclusive. The period
-- length of the generated numbers is 2^85 - 2^30. The low order bits of
-- the generated numbers are just as random as the high-order bits. */

int rng_next_rand(RNG *rand)
{     return
         *fptr >= 0 ? *fptr-- : flip_cycle(rand);
}

/*----------------------------------------------------------------------
-- rng_unif_rand - obtain pseudo-random integer in [0, m-1].
--
-- *Synopsis*
--
-- #include "glprng.h"
-- int rng_unif_rand(RNG *rand, int m);
--
-- *Returns*
--
-- The routine rng_unif_rand returns a next pseudo-random integer which
-- is uniformly distributed between 0 and m-1, inclusive, where m is any
-- positive integer less than 2^31. */

#define two_to_the_31 ((unsigned int)0x80000000)

int rng_unif_rand(RNG *rand, int m)
{     unsigned int t = two_to_the_31 - (two_to_the_31 % m);
      int r;
      do { r = rng_next_rand(rand); } while (t <= (unsigned int)r);
      return r % m;
}

/*----------------------------------------------------------------------
-- rng_delete_rand - delete pseudo-random number generator.
--
-- *Synopsis*
--
-- #include "glprng.h"
-- void rng_delete_rand(RNG *rand);
--
-- *Description*
--
-- The routine rng_delete_rand frees all the memory allocated to the
-- pseudo-random number generator. */

void rng_delete_rand(RNG *rand)
{     ufree(rand);
      return;
}

/*--------------------------------------------------------------------*/

#if 0
int main(void)
{     /* to be sure that this version provides the same results as the
         original version, run this validation program */
      RNG *rand;
      int j;
      rand = rng_create_rand();
      rng_init_rand(rand, -314159);
      if (rng_next_rand(rand) != 119318998)
      {  print("Failure on the first try!");
         return -1;
      }
      for (j = 1; j <= 133; j++) rng_next_rand(rand);
      if (rng_unif_rand(rand, 0x55555555) != 748103812)
      {  print("Failure on the second try!");
         return -2;
      }
      rng_delete_rand(rand);
      print("OK, the random-number generator routines seem to work!");
      return 0;
}
#endif

/* eof */
