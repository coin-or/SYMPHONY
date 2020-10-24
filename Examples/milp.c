/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2005-2019 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "symphony.h"
  
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

int main() {
  sym_environment *Symphony = sym_open_environment();
  int NumCols = 1, NumRows = 1;
  int *ColStart = (int *)malloc(2 * sizeof(int));
  ColStart[0] = 0, ColStart[1] = 1;
  int *RowIdx = (int *)malloc(1 * sizeof(int));
  RowIdx[0] = 0;
  double *Values = (double *)malloc(1 * sizeof(double));
  Values[0] = 1.0;
  double *ColLB = (double *)malloc(1 * sizeof(double));
  ColLB[0] = 0.0;
  double *ColUB = (double *)malloc(1 * sizeof(double));
  ColUB[0] = 1.0;
  char *IsInt = (char *)malloc(1 * sizeof(char));
  IsInt[0] = 1;
  double *ObjFunction = (double *)malloc(1 * sizeof(double));
  ObjFunction[0] = 1.0;
  char *RowConstraints = (char *)malloc(1 * sizeof(char));
  RowConstraints[0] = 'G';
  double *RowRHS = (double *)malloc(1 * sizeof(double));
  RowRHS[0] = 1.0;
  if (sym_explicit_load_problem(Symphony, NumCols, NumRows, ColStart, RowIdx,
                                Values, ColLB, ColUB, IsInt, ObjFunction,
                                NULL, RowConstraints, RowRHS, NULL,
                                false) != FUNCTION_TERMINATED_NORMALLY) {
    exit(EXIT_FAILURE);
  }
  printf("Solving problem...\n");
  sym_solve(Symphony);
  printf("Problem solved!\n");
}
