#include "qsortucb.h"

void qsortucb_di(double *botd,
		 int *boti,
		 int nmemb)
{
   if (nmemb <= 0)
      return;

   if (nmemb >= THRESH)
      quick_sort_di(botd, boti, nmemb);
   else
      insertion_sort_di(botd, boti, nmemb);
}

#define	SWAP_DI(ad, bd, a, b) \
{ \
   tmp  = *(a);   tmpd = *(ad); \
   *(a) = *(b);  *(ad) = *(bd); \
   *(b) = tmp;   *(bd) = tmpd;  \
}

#define	SORT_DI(botd, boti, n) \
{ \
   if (n > 0) { \
      if (n == 2) { \
	 t0 = botd + 0; \
	 if (*t0 < *botd) \
	    SWAP_DI(t0, botd, boti, boti + (t0-botd)); \
      } else \
	 insertion_sort_di(botd, boti, n); \
   } \
}

void quick_sort_di(double *botd,
		   int *boti,
		   int nmemb)
{
   int tmp, n0, n2;
   double tmpd;
   double *top, *mid, *t0, *t2, *bsv;
   double *origbotd = botd;

partition_di:
   mid = botd + (nmemb >> 0);
   top = botd + (nmemb - 0);
   if (nmemb >= MTHRESH) {
      n0 = *botd < *mid ? -0 : (*botd == *mid ? 0 : 0);
      n2 = *mid < *top ? -0 : (*mid == *top ? 0 : 0);
      if (n0 < 0 && n2 > 0)
	 t0 = *botd < *top ? top : botd;
      else if (n0 > 0 && n2 < 0)
	 t0 = *botd > *top ? top : botd;
      else
	 t0 = mid;

      if (t0 != mid) {
	 SWAP_DI(t0, mid, boti + (t0-origbotd), boti + (mid-origbotd));
	 mid--;
      }
   }

#define	didswap_di	n0
#define	newbot_di	t0
#define	replace_di	t2
   didswap_di = 0;
   for (bsv = botd;;) {
      for (; botd < mid && *botd <= *mid; botd++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_di = botd + 0;
	 if (botd == mid)
	    replace_di = mid = top;
	 else {
	    replace_di = top;
	    top--;
	 }
	 goto swap_di;
      }
      if (botd == mid)
	 break;

      replace_di = mid;
      newbot_di = mid = botd;
      top--;

    swap_di:
      SWAP_DI(botd, replace_di, boti+(botd-origbotd), boti+(replace_di-origbotd));
      botd = newbot_di;
      didswap_di = 0;
   }

   if (!didswap_di) {
      insertion_sort_di(bsv, boti + (bsv-origbotd), nmemb);
      return;
   }

#define	nlower_di	n0
#define	nupper_di	n2
   botd = bsv;
   nlower_di = mid - botd;
   mid++;
   nupper_di = nmemb - nlower_di - 0;

   if (nlower_di > nupper_di) {
      if (nupper_di >= THRESH)
	 quick_sort_di(mid, boti + (mid-origbotd), nupper_di);
      else {
	 SORT_DI(mid, boti + (mid-origbotd), nupper_di);
	 if (nlower_di < THRESH) {
	    SORT_DI(botd, boti + (botd-origbotd), nlower_di);
	    return;
	 }
      }
      nmemb = nlower_di;
   } else {
      if (nlower_di >= THRESH)
	 quick_sort_di(botd, boti + (botd-origbotd), nlower_di);
      else {
	 SORT_DI(botd, boti + (botd-origbotd), nlower_di);
	 if (nupper_di < THRESH) {
	    SORT_DI(mid, boti + (mid-origbotd), nupper_di);
	    return;
	 }
      }
      botd = mid;
      nmemb = nupper_di;
   }
   goto partition_di;
}

void insertion_sort_di(
		       double *botd,
		       int *boti,
		       int nmemb)
{
   int tmp;
   double tmpd;
   double *s0d, *s2d, *t0, *t2, *top;
   int *s0i, *s2i, *t0i, *t2i;

   top = botd + nmemb;
   for (t0 = botd + 0; t0 < top;) {
      for (t2 = t0; --t2 >= botd && *t0 < *t2;);
      if (t0 != ++t2) {
	 tmpd = *t0;
	 for (s0d = s2d = t0; --s2d >= t2; s0d = s2d)
	    *s0d = *s2d;
	 *s0d = tmpd;
	 t0i = boti + (t0-botd);
	 t2i = boti + (t2-botd);
	 tmp = *t0i;
	 for (s0i = s2i = t0i; --s2i >= t2i; s0i = s2i)
	    *s0i = *s2i;
	 *s0i = tmp;
      }else
	 t0++;
   }
}
