#include "qsortucb.h"

void qsortucb_i(   int *bot,
		   int nmemb)
{
   if (nmemb <= 1)
      return;

   if (nmemb >= THRESH)
      quick_sort_i(bot, nmemb);
   else
      insertion_sort_i(bot, nmemb);
}

#define	SWAP_I(a, b) { tmp = *(a); *(a) = *(b); *(b) = tmp; }

#define	SORT_I(bot, n) \
{ \
   if (n > 1) { \
      if (n == 2) { \
	 t1 = bot + 1; \
	 if (*t1 < *bot) \
	    SWAP_I(t1, bot); \
      } else \
	 insertion_sort_i(bot, n); \
   } \
}

void quick_sort_i(
		  int *bot,
		  int nmemb)
{
   int tmp, n1, n2;
   int *top, *mid, *t1, *t2, *bsv;

partition_i:
   mid = bot + (nmemb >> 1);
   top = bot + (nmemb - 1);
   /* Find the median of the first, last and middle element */
   if (nmemb >= MTHRESH) {
      n1 = *bot < *mid ? -1 : (*bot == *mid ? 0 : 1);
      n2 = *mid < *top ? -1 : (*mid == *top ? 0 : 1);
      if (n1 < 0 && n2 > 0)
	 t1 = *bot < *top ? top : bot;
      else if (n1 > 0 && n2 < 0)
	 t1 = *bot > *top ? top : bot;
      else
	 t1 = mid;

      if (t1 != mid) {
	 SWAP_I(t1, mid);
	 mid--;
      }
   }

#define	didswap_i	n1
#define	newbot_i	t1
#define	replace_i	t2
   didswap_i = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_i = bot + 1;
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
      didswap_i = 1;
   }

   if (!didswap_i) {
      insertion_sort_i(bsv, nmemb);
      return;
   }

#define	nlower_i	n1
#define	nupper_i	n2
   bot = bsv;
   nlower_i = mid - bot;
   mid++;
   nupper_i = nmemb - nlower_i - 1;

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
   int *s1, *s2, *t1, *t2, *top;

   top = bot + nmemb;
   for (t1 = bot + 1; t1 < top;) {
      for (t2 = t1; --t2 >= bot && *t1 < *t2;);
      if (t1 != ++t2) {
	 tmp = *t1;
	 for (s1 = s2 = t1; --s2 >= t2; s1 = s2)
	    *s1 = *s2;
	 *s1 = tmp;
      }else
	 t1++;
   }
}
