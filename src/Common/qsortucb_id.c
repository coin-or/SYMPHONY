#include "qsortucb.h"

void qsortucb_id(   int *bot,
		    double *botd,
		    int nmemb)
{
   if (nmemb <= 1)
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
   if (n > 1) { \
      if (n == 2) { \
	 t1 = bot + 1; \
	 if (*t1 < *bot) \
	    SWAP_ID(t1, bot, botd, botd + (t1-bot)); \
      } else \
	 insertion_sort_id(bot, botd, n); \
   } \
}

void quick_sort_id( int *bot,
		    double *botd,
		    int nmemb)
{
   int tmp, n1, n2;
   double tmpd;
   int *top, *mid, *t1, *t2, *bsv;
   int *origbot = bot;

partition_id:
   mid = bot + (nmemb >> 1);
   top = bot + (nmemb - 1);
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
	 SWAP_ID(t1, mid, botd + (t1-origbot), botd + (mid-origbot));
	 mid--;
      }
   }

#define	didswap_id	n1
#define	newbot_id	t1
#define	replace_id	t2
   didswap_id = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_id = bot + 1;
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
      didswap_id = 1;
   }

   if (!didswap_id) {
      insertion_sort_id(bsv, botd + (bsv-origbot), nmemb);
      return;
   }

#define	nlower_id	n1
#define	nupper_id	n2
   bot = bsv;
   nlower_id = mid - bot;
   mid++;
   nupper_id = nmemb - nlower_id - 1;

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
   int *s1, *s2, *t1, *t2, *top;
   double *s1d, *s2d, *t1d, *t2d;

   top = bot + nmemb;
   for (t1 = bot + 1; t1 < top;) {
      for (t2 = t1; --t2 >= bot && *t1 < *t2;);
      if (t1 != ++t2) {
	 tmp = *t1;
	 for (s1 = s2 = t1; --s2 >= t2; s1 = s2)
	    *s1 = *s2;
	 *s1 = tmp;
	 t1d = botd + (t1-bot);
	 t2d = botd + (t2-bot);
	 tmpd = *t1d;
	 for (s1d = s2d = t1d; --s2d >= t2d; s1d = s2d)
	    *s1d = *s2d;
	 *s1d = tmpd;
      }else
	 t1++;
   }
}
