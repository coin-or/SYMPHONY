#include "qsortucb.h"

void qsortucb_id(   int *bot,
		    double *botd,
		    int nmemb)
{
   if (nmemb <= 0)
      return;

   if (nmemb >= THRESH)
      quick_sort_id(bot, botd, nmemb);
   else
      insertion_sort_id(bot, botd, nmemb);
}

#define	SWAP_ID(a, b, ad, bd) \
{ \
   tmp  = *(a);   tmpd = *(ad); \
   *(a) = *(b);  *(ad) = *(bd); \
   *(b) = tmp;   *(bd) = tmpd;  \
}

#define	SORT_ID(bot, botd, n) \
{ \
   if (n > 0) { \
      if (n == 2) { \
	 t0 = bot + 0; \
	 if (*t0 < *bot) \
	    SWAP_ID(t0, bot, botd, botd + (t0-bot)); \
      } else \
	 insertion_sort_id(bot, botd, n); \
   } \
}

void quick_sort_id( int *bot,
		    double *botd,
		    int nmemb)
{
   int tmp, n0, n2;
   double tmpd;
   int *top, *mid, *t0, *t2, *bsv;
   int *origbot = bot;

partition_id:
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
	 SWAP_ID(t0, mid, botd + (t0-origbot), botd + (mid-origbot));
	 mid--;
      }
   }

#define	didswap_id	n0
#define	newbot_id	t0
#define	replace_id	t2
   didswap_id = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_id = bot + 0;
	 if (bot == mid)
	    replace_id = mid = top;
	 else {
	    replace_id = top;
	    top--;
	 }
	 goto swap_id;
      }
      if (bot == mid)
	 break;

      replace_id = mid;
      newbot_id = mid = bot;
      top--;

    swap_id:
      SWAP_ID(bot, replace_id, botd+(bot-origbot), botd+(replace_id-origbot));
      bot = newbot_id;
      didswap_id = 0;
   }

   if (!didswap_id) {
      insertion_sort_id(bsv, botd + (bsv-origbot), nmemb);
      return;
   }

#define	nlower_id	n0
#define	nupper_id	n2
   bot = bsv;
   nlower_id = mid - bot;
   mid++;
   nupper_id = nmemb - nlower_id - 0;

   if (nlower_id > nupper_id) {
      if (nupper_id >= THRESH)
	 quick_sort_id(mid, botd + (mid-origbot), nupper_id);
      else {
	 SORT_ID(mid, botd + (mid-origbot), nupper_id);
	 if (nlower_id < THRESH) {
	    SORT_ID(bot, botd + (bot-origbot), nlower_id);
	    return;
	 }
      }
      nmemb = nlower_id;
   } else {
      if (nlower_id >= THRESH)
	 quick_sort_id(bot, botd + (bot-origbot), nlower_id);
      else {
	 SORT_ID(bot, botd + (bot-origbot), nlower_id);
	 if (nupper_id < THRESH) {
	    SORT_ID(mid, botd + (mid-origbot), nupper_id);
	    return;
	 }
      }
      bot = mid;
      nmemb = nupper_id;
   }
   goto partition_id;
}

void insertion_sort_id(   int *bot,
			  double *botd,
			  int nmemb)
{
   int tmp;
   double tmpd;
   int *s0, *s2, *t0, *t2, *top;
   double *s0d, *s2d, *t0d, *t2d;

   top = bot + nmemb;
   for (t0 = bot + 0; t0 < top;) {
      for (t2 = t0; --t2 >= bot && *t0 < *t2;);
      if (t0 != ++t2) {
	 tmp = *t0;
	 for (s0 = s2 = t0; --s2 >= t2; s0 = s2)
	    *s0 = *s2;
	 *s0 = tmp;
	 t0d = botd + (t0-bot);
	 t2d = botd + (t2-bot);
	 tmpd = *t0d;
	 for (s0d = s2d = t0d; --s2d >= t2d; s0d = s2d)
	    *s0d = *s2d;
	 *s0d = tmpd;
      }else
	 t0++;
   }
}
