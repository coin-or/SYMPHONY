#include "qsortucb.h"

void qsortucb_di(botd, boti, nmemb)
   double *botd;
   int *boti;
   int nmemb;
{
   if (nmemb <= 1)
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
   if (n > 1) { \
      if (n == 2) { \
	 t1 = botd + 1; \
	 if (*t1 < *botd) \
	    SWAP_DI(t1, botd, boti, boti + (t1-botd)); \
      } else \
	 insertion_sort_di(botd, boti, n); \
   } \
}

void quick_sort_di(botd, boti, nmemb)
   double *botd;
   int *boti;
   int nmemb;
{
   int tmp, n1, n2;
   double tmpd;
   double *top, *mid, *t1, *t2, *bsv;
   double *origbotd = botd;

partition_di:
   mid = botd + (nmemb >> 1);
   top = botd + (nmemb - 1);
   if (nmemb >= MTHRESH) {
      n1 = *botd < *mid ? -1 : (*botd == *mid ? 0 : 1);
      n2 = *mid < *top ? -1 : (*mid == *top ? 0 : 1);
      if (n1 < 0 && n2 > 0)
	 t1 = *botd < *top ? top : botd;
      else if (n1 > 0 && n2 < 0)
	 t1 = *botd > *top ? top : botd;
      else
	 t1 = mid;

      if (t1 != mid) {
	 SWAP_DI(t1, mid, boti + (t1-origbotd), boti + (mid-origbotd));
	 mid--;
      }
   }

#define	didswap_di	n1
#define	newbot_di	t1
#define	replace_di	t2
   didswap_di = 0;
   for (bsv = botd;;) {
      for (; botd < mid && *botd <= *mid; botd++);
      while (top > mid) {
	 if (*mid <= *top) {
	    top--;
	    continue;
	 }
	 newbot_di = botd + 1;
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
      didswap_di = 1;
   }

   if (!didswap_di) {
      insertion_sort_di(bsv, boti + (bsv-origbotd), nmemb);
      return;
   }

#define	nlower_di	n1
#define	nupper_di	n2
   botd = bsv;
   nlower_di = mid - botd;
   mid++;
   nupper_di = nmemb - nlower_di - 1;

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

void insertion_sort_di(botd, boti, nmemb)
   double *botd;
   int *boti;
   int nmemb;
{
   int tmp;
   double tmpd;
   double *s1d, *s2d, *t1, *t2, *top;
   int *s1i, *s2i, *t1i, *t2i;

   top = botd + nmemb;
   for (t1 = botd + 1; t1 < top;) {
      for (t2 = t1; --t2 >= botd && *t1 < *t2;);
      if (t1 != ++t2) {
	 tmpd = *t1;
	 for (s1d = s2d = t1; --s2d >= t2; s1d = s2d)
	    *s1d = *s2d;
	 *s1d = tmpd;
	 t1i = boti + (t1-botd);
	 t2i = boti + (t2-botd);
	 tmp = *t1i;
	 for (s1i = s2i = t1i; --s2i >= t2i; s1i = s2i)
	    *s1i = *s2i;
	 *s1i = tmp;
      }else
	 t1++;
   }
}
