#include "qsortucb.h"

void qsortucb_ii(bot, bota, nmemb)
   int *bot;
   int *bota;
   int nmemb;
{
   if (nmemb <= 1)
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
   if (n > 1) { \
      if (n == 2) { \
	 t1 = bot + 1; \
	 if (*t1 < *bot) \
	    SWAP_II(t1, bot, bota, bota + (t1-bot)); \
      } else \
	 insertion_sort_ii(bot, bota, n); \
   } \
}

void quick_sort_ii(bot, bota, nmemb)
   int *bot;
   int *bota;
   int nmemb;
{
   int tmp, n1, n2;
   int tmpa;
   int *top, *mid, *t1, *t2, *bsv;
   int *origbot = bot;

partition_ii:
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
	 SWAP_II(t1, mid, bota + (t1-origbot), bota + (mid-origbot));
	 mid--;
      }
   }

#define	didswap_ii	n1
#define	newbot_ii	t1
#define	replace_ii	t2
   didswap_ii = 0;
   for (bsv = bot;;) {
      for (; bot < mid && *bot <= *mid; bot++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_ii = bot + 1;
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
      didswap_ii = 1;
   }

   if (!didswap_ii) {
      insertion_sort_ii(bsv, bota + (bsv-origbot), nmemb);
      return;
   }

#define	nlower_ii	n1
#define	nupper_ii	n2
   bot = bsv;
   nlower_ii = mid - bot;
   mid++;
   nupper_ii = nmemb - nlower_ii - 1;

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

void insertion_sort_ii(bot, bota, nmemb)
   int *bot;
   int *bota;
   int nmemb;
{
   int tmp;
   int tmpa;
   int *s1, *s2, *t1, *t2, *t1a, *t2a, *top;

   top = bot + nmemb;
   for (t1 = bot + 1; t1 < top;) {
      for (t2 = t1; --t2 >= bot && *t1 < *t2;);
      if (t1 != ++t2) {
	 tmp = *t1;
	 for (s1 = s2 = t1; --s2 >= t2; s1 = s2)
	    *s1 = *s2;
	 *s1 = tmp;
	 t1a = bota + (t1-bot);
	 t2a = bota + (t2-bot);
	 tmpa = *t1a;
	 for (s1 = s2 = t1a; --s2 >= t2a; s1 = s2)
	    *s1 = *s2;
	 *s1 = tmpa;
      }else
	 t1++;
   }
}
