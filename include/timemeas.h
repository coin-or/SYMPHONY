/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2004 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef __TIMEMEAS_H
#define __TIMEMEAS_H

#ifdef WIN32
#include "win32_time.h"
#else
#include <sys/time.h>
#endif

#include "proto.h"

#define PRINT_TIME(tm, f) { /* Print the elapsed time in vbctool format*/    \
   double elapsed = wall_clock(NULL) - tm->start_time;                       \
   int hours, minutes, seconds, msec;                                        \
   hours = (int)(elapsed/3600.0);                                            \
   elapsed -= hours*3600.0;                                                  \
   minutes = (int)(elapsed/60.0);                                            \
   elapsed -= minutes*60.0;                                                  \
   seconds = (int)elapsed;                                                   \
   elapsed -= (double)seconds;                                               \
   msec = (int)(elapsed*100.0);                                              \
   fprintf(f, "%.2d:%.2d:%.2d:%.2d ", hours, minutes, seconds, msec);         \
}

#define	TVCLEAR(tvp)	(tvp.tv_sec = tvp.tv_usec = 0)
#define	PTVCLEAR(tvp)	((tvp)->tv_sec = (tvp)->tv_usec = 0)

#define	TVISSET(tvp)	(tvp.tv_sec || tvp.tv_usec)
#define	PTVISSET(tvp)	((tvp)->tv_sec || (tvp)->tv_usec)

#define	TVXLTY(xtv, ytv) 						\
   ( (xtv.tv_sec < ytv.tv_sec) ||					\
     (xtv.tv_sec == ytv.tv_sec && xtv.tv_usec < ytv.tv_usec))
#define	PTVXLTY(xtv, ytv) 						  \
   ( ((xtv)->tv_sec < (ytv)->tv_sec) ||					  \
     ((xtv)->tv_sec == (ytv)->tv_sec && (xtv)->tv_usec < (ytv)->tv_usec))

#define	TVXADDY(ztv, xtv, ytv)						\
     if ((ztv.tv_usec = xtv.tv_usec + ytv.tv_usec) < 1000000) {		\
	ztv.tv_sec = xtv.tv_sec + ytv.tv_sec;				\
     } else {								\
	ztv.tv_usec -= 1000000;						\
	ztv.tv_sec = xtv.tv_sec + ytv.tv_sec + 1;			\
     }
#define	PTVXADDY(ztv, xtv, ytv)						 \
     if (((ztv)->tv_usec = (xtv)->tv_usec + (ytv)->tv_usec) < 1000000) { \
	(ztv)->tv_sec = (xtv)->tv_sec + (ytv)->tv_sec;			 \
     } else {								 \
	(ztv)->tv_usec -= 1000000;					 \
	(ztv)->tv_sec = (xtv)->tv_sec + (ytv)->tv_sec + 1;		 \
     }

#define	TVXSUBY(ztv, xtv, ytv)						\
     if (xtv.tv_usec >= ytv.tv_usec) {					\
	ztv.tv_sec = xtv.tv_sec - ytv.tv_sec;				\
	ztv.tv_usec = xtv.tv_usec - ytv.tv_usec;			\
     } else {								\
	ztv.tv_sec = xtv.tv_sec - ytv.tv_sec - 1;			\
	ztv.tv_usec = xtv.tv_usec + 1000000 - ytv.tv_usec;		\
     }
#define	PTVXSUBY(ztv, xtv, ytv)						 \
     if ((xtv)->tv_usec >= (ytv)->tv_usec) {				 \
	(ztv)->tv_sec = (xtv)->tv_sec - (ytv)->tv_sec;			 \
	(ztv)->tv_usec = (xtv)->tv_usec - (ytv)->tv_usec;		 \
     } else {								 \
	(ztv)->tv_sec = (xtv)->tv_sec - (ytv)->tv_sec - 1;		 \
	(ztv)->tv_usec = (xtv)->tv_usec + 1000000 - (ytv)->tv_usec;	 \
     }

#define TVTODBL(tvp)  ((double)tvp.tv_sec + ((double)tvp.tv_usec)/1000000 )
#define TVPTODBL(tvp) ((double)(tvp)->tv_sec+((double)(tvp)->tv_usec)/1000000)

#define DBLTOTV(x, tv)						\
     tv.tv_sec = (int) floor(x);				\
     tv.tv_usec = (int) floor(1000000 * (x - tv.tv_sec));
#define DBLTOPTV(x, tvp)						\
     (tvp)->tv_sec = (int) floor(x);					\
     (tvp)->tv_usec = (int) floor(1000000 * (x - (tvp)->tv_sec));

void start_time PROTO((void));
double used_time PROTO((double *T));
double wall_clock PROTO((double *T));

#endif
