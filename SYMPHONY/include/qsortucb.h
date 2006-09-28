/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2006 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _QSORTUCB_H
#define _QSORTUCB_H

/*
 * MTHRESH is the smallest partition for which we compare for a median
 * value instead of using the middle value.
 */
#define	MTHRESH 6

/*
 * THRESH is the minimum number of entries in a partition for continued
 * partitioning.
 */
#define	THRESH  4

void qsortucb(char *bot, unsigned int nmemb, int size,
	      int (*compar)(const void *, const void *));

void insertion_sort(char *bot, unsigned int nmemb, int size,
		    int (*compar)(const void *, const void *));

void quick_sort(char *bot, unsigned int nmemb, int size,
		int (*compar)(const void *, const void *));

void qsortucb_i(int *bot, int nmemb);
void quick_sort_i(int *bot, int nmemb);
void insertion_sort_i(int *bot, int nmemb);

void qsortucb_id(int *bot, double *botd, int nmemb);
void quick_sort_id(int *bot, double *botd, int nmemb);
void insertion_sort_id(int *bot, double *botd, int nmemb);

void qsortucb_ic(int *bot, char *botc, int nmemb);
void quick_sort_ic(int *bot, char *botc, int nmemb);
void insertion_sort_ic(int *bot, char *botc, int nmemb);

void qsortucb_ii(int *bot, int *bota, int nmemb);
void quick_sort_ii(int *bot, int *bota, int nmemb);
void insertion_sort_ii(int *bot, int *bota, int nmemb);

void qsortucb_di(double *botd, int *boti, int nmemb);
void quick_sort_di(double *botd, int *boti, int nmemb);
void insertion_sort_di(double *botd, int *boti, int nmemb);

#endif
