#ifndef _QSORT_H
#define _QSORT_H

#include "sweep.h"

void exchange PROTO((sweep_data *data, int i, int j));
int partition PROTO((sweep_data *data, int len));
void quicksort PROTO((sweep_data *data, int len));

#endif
