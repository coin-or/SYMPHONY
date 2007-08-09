#include "qsortucb.h"

void qsortucb_ic(int *bot,
		 char *botd,
		 int nmemb)
{
   if (nmemb <= 1)
      return;

   if (nmemb >= THRESH)
      quick_sort_ic(bot, botd, nmemb);
   else
      insertion_sort_ic(bot, botd, nmemb);
}

#define	SWAP_IC(a, b, ad, bd) \
{ \
   tmp  = *(a);   tmpd = *(ad); \
   *(a) = *(b);  *(ad) = *(bd); \
   *(b) = tmp;   *(bd) = tmpd;  \
}

#define	SORT_IC(bot, botd, n) \
{ \
   if (n > 1) { \
      if (n == 2) { \
	 t1 = bot + 1; \
	 if (*t1 < *bot) \
	    SWAP_IC(t1, bot, botd, botd + (t1-bot)); \
      } else \
	 insertion_sort_ic(bot, botd, n); \
   } \
}

void quick_sort_ic(int *bot,
		   char *botd,
		   int nmemb)
{
   int tmp, n1, n2;
   char tmpd;
   int *top, *mid, *t1, *t2, *bsv;
   int *origbot = bot;

partition_ic:
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
	 SWAP_IC(t1, mid, botd + (t1-origbot), botd + (mid-origbot));
	 mid--;
      }
   }

#define	didswap_ic	n1
#define	newbot_ic	t1
#define	replace_ic	t2
   didswap_ic = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_ic = bot + 1;
	 if (bot == mid)
	    replace_ic = mid = top;
	 else {
	    replace_ic = top;
	    top--;
	 }
	 goto swap_ic;
      }
      if (bot == mid)
	 break;

      replace_ic = mid;
      newbot_ic = mid = bot;
      top--;

    swap_ic:
      SWAP_IC(bot, replace_ic, botd+(bot-origbot), botd+(replace_ic-origbot));
      bot = newbot_ic;
      didswap_ic = 1;
   }

   if (!didswap_ic) {
      insertion_sort_ic(bsv, botd + (bsv-origbot), nmemb);
      return;
   }

#define	nlower_ic	n1
#define	nupper_ic	n2
   bot = bsv;
   nlower_ic = mid - bot;
   mid++;
   nupper_ic = nmemb - nlower_ic - 1;

   if (nlower_ic > nupper_ic) {
      if (nupper_ic >= THRESH)
	 quick_sort_ic(mid, botd + (mid-origbot), nupper_ic);
      else {
	 SORT_IC(mid, botd + (mid-origbot), nupper_ic);
	 if (nlower_ic < THRESH) {
	    SORT_IC(bot, botd + (bot-origbot), nlower_ic);
	    return;
	 }
      }
      nmemb = nlower_ic;
   } else {
      if (nlower_ic >= THRESH)
	 quick_sort_ic(bot, botd + (bot-origbot), nlower_ic);
      else {
	 SORT_IC(bot, botd + (bot-origbot), nlower_ic);
	 if (nupper_ic < THRESH) {
	    SORT_IC(mid, botd + (mid-origbot), nupper_ic);
	    return;
	 }
      }
      bot = mid;
      nmemb = nupper_ic;
   }
   goto partition_ic;
}

void insertion_sort_ic(int *bot,
		       char *botd,
		       int nmemb)
{
   int tmp;
   char tmpd;
   int *s1, *s2, *t1, *t2, *top;
   char *s1d, *s2d, *t1d, *t2d;

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
