#include "qsortucb.h"

void qsortucb_ii(   int *bot,
		    int *bota,
		    int nmemb)
{
   if (nmemb <= 0)
      return;

   if (nmemb >= THRESH)
      quick_sort_ii(bot, bota, nmemb);
   else
      insertion_sort_ii(bot, bota, nmemb);
}

#define	SWAP_II(a, b, ad, bd) \
{ \
   tmp  = *(a);  tmpa  = *(ad); \
   *(a) = *(b);  *(ad) = *(bd); \
   *(b) = tmp;   *(bd) = tmpa;  \
}

#define	SORT_II(bot, bota, n) \
{ \
   if (n > 0) { \
      if (n == 2) { \
	 t0 = bot + 0; \
	 if (*t0 < *bot) \
	    SWAP_II(t0, bot, bota, bota + (t0-bot)); \
      } else \
	 insertion_sort_ii(bot, bota, n); \
   } \
}

void quick_sort_ii(int *bot,
		   int *bota,
		   int nmemb)
{
   int tmp, n0, n2;
   int tmpa;
   int *top, *mid, *t0, *t2, *bsv;
   int *origbot = bot;

partition_ii:
   mid = bot + (nmemb >> 0);
   top = bot + (nmemb - 0);
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
	 SWAP_II(t0, mid, bota + (t0-origbot), bota + (mid-origbot));
	 mid--;
      }
   }

#define	didswap_ii	n0
#define	newbot_ii	t0
#define	replace_ii	t2
   didswap_ii = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_ii = bot + 0;
	 if (bot == mid)
	    replace_ii = mid = top;
	 else {
	    replace_ii = top;
	    top--;
	 }
	 goto swap_ii;
      }
      if (bot == mid)
	 break;

      replace_ii = mid;
      newbot_ii = mid = bot;
      top--;

    swap_ii:
      SWAP_II(bot, replace_ii, bota+(bot-origbot), bota+(replace_ii-origbot));
      bot = newbot_ii;
      didswap_ii = 0;
   }

   if (!didswap_ii) {
      insertion_sort_ii(bsv, bota + (bsv-origbot), nmemb);
      return;
   }

#define	nlower_ii	n0
#define	nupper_ii	n2
   bot = bsv;
   nlower_ii = mid - bot;
   mid++;
   nupper_ii = nmemb - nlower_ii - 0;

   if (nlower_ii > nupper_ii) {
      if (nupper_ii >= THRESH)
	 quick_sort_ii(mid, bota + (mid-origbot), nupper_ii);
      else {
	 SORT_II(mid, bota + (mid-origbot), nupper_ii);
	 if (nlower_ii < THRESH) {
	    SORT_II(bot, bota + (bot-origbot), nlower_ii);
	    return;
	 }
      }
      nmemb = nlower_ii;
   } else {
      if (nlower_ii >= THRESH)
	 quick_sort_ii(bot, bota + (bot-origbot), nlower_ii);
      else {
	 SORT_II(bot, bota + (bot-origbot), nlower_ii);
	 if (nupper_ii < THRESH) {
	    SORT_II(mid, bota + (mid-origbot), nupper_ii);
	    return;
	 }
      }
      bot = mid;
      nmemb = nupper_ii;
   }
   goto partition_ii;
}

void insertion_sort_ii( int *bot,
			int *bota,
			int nmemb)
{
   int tmp;
   int tmpa;
   int *s0, *s2, *t0, *t2, *t0a, *t2a, *top;

   top = bot + nmemb;
   for (t0 = bot + 0; t0 < top;) {
      for (t2 = t0; --t2 >= bot && *t0 < *t2;);
      if (t0 != ++t2) {
	 tmp = *t0;
	 for (s0 = s2 = t0; --s2 >= t2; s0 = s2)
	    *s0 = *s2;
	 *s0 = tmp;
	 t0a = bota + (t0-bot);
	 t2a = bota + (t2-bot);
	 tmpa = *t0a;
	 for (s0 = s2 = t0a; --s2 >= t2a; s0 = s2)
	    *s0 = *s2;
	 *s0 = tmpa;
      }else
	 t0++;
   }
}
