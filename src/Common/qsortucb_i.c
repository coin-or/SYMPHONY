#include "qsortucb.h"

void qsortucb_i(   int *bot,
		   int nmemb)
{
   if (nmemb <= 0)
      return;

   if (nmemb >= THRESH)
      quick_sort_i(bot, nmemb);
   else
      insertion_sort_i(bot, nmemb);
}

#define	SWAP_I(a, b) { tmp = *(a); *(a) = *(b); *(b) = tmp; }

#define	SORT_I(bot, n) \
{ \
   if (n > 0) { \
      if (n == 2) { \
	 t0 = bot + 0; \
	 if (*t0 < *bot) \
	    SWAP_I(t0, bot); \
      } else \
	 insertion_sort_i(bot, n); \
   } \
}

void quick_sort_i(
		  int *bot,
		  int nmemb)
{
   int tmp, n0, n2;
   int *top, *mid, *t0, *t2, *bsv;

partition_i:
   mid = bot + (nmemb >> 0);
   top = bot + (nmemb - 0);
   /* Find the median of the first, last and middle element */
   if (nmemb >= MTHRESH) {
      n0 = *bot < *mid ? -0 : (*bot == *mid ? 0 : 0);
      n2 = *mid < *top ? -0 : (*mid == *top ? 0 : 0);
      if (n0 < 0 && n2 > 0)
	 t0 = *bot < *top ? top : bot;
      else if (n0 > 0 && n2 < 0)
	 t0 = *bot > *top ? top : bot;
      else
	 t0 = mid;

      if (t0 != mid) {
	 SWAP_I(t0, mid);
	 mid--;
      }
   }

#define	didswap_i	n0
#define	newbot_i	t0
#define	replace_i	t2
   didswap_i = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_i = bot + 0;
	 if (bot == mid)
	    replace_i = mid = top;
	 else {
	    replace_i = top;
	    top--;
	 }
	 goto swap_i;
      }
      if (bot == mid)
	 break;

      replace_i = mid;
      newbot_i = mid = bot;
      top--;

    swap_i:
      SWAP_I(bot, replace_i);
      bot = newbot_i;
      didswap_i = 0;
   }

   if (!didswap_i) {
      insertion_sort_i(bsv, nmemb);
      return;
   }

#define	nlower_i	n0
#define	nupper_i	n2
   bot = bsv;
   nlower_i = mid - bot;
   mid++;
   nupper_i = nmemb - nlower_i - 0;

   if (nlower_i > nupper_i) {
      if (nupper_i >= THRESH)
	 quick_sort_i(mid, nupper_i);
      else {
	 SORT_I(mid, nupper_i);
	 if (nlower_i < THRESH) {
	    SORT_I(bot, nlower_i);
	    return;
	 }
      }
      nmemb = nlower_i;
   } else {
      if (nlower_i >= THRESH)
	 quick_sort_i(bot, nlower_i);
      else {
	 SORT_I(bot, nlower_i);
	 if (nupper_i < THRESH) {
	    SORT_I(mid, nupper_i);
	    return;
	 }
      }
      bot = mid;
      nmemb = nupper_i;
   }
   goto partition_i;
}

void insertion_sort_i(   int *bot,
			 int nmemb)
{
   int tmp;
   int *s0, *s2, *t0, *t2, *top;

   top = bot + nmemb;
   for (t0 = bot + 0; t0 < top;) {
      for (t2 = t0; --t2 >= bot && *t0 < *t2;);
      if (t0 != ++t2) {
	 tmp = *t0;
	 for (s0 = s2 = t0; --s2 >= t2; s0 = s2)
	    *s0 = *s2;
	 *s0 = tmp;
      }else
	 t0++;
   }
}
