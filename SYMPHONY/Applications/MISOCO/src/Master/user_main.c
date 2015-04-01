/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/*===========================================================================*/

#define CALL_FUNCTION(f) \
if ((termcode = f) < 0){                                                    \
   printf("Error detected: termcode = %i\n", termcode);                     \
   printf("Exiting...\n\n");                                                \
   exit(termcode);                                                          \
}

// symphony headers
#include "symphony.h"
#include "sym_master.h"
// application header
#include "user.h"
// C library headers
#include <stdlib.h>

int user_read_data(user_problem * prob, char * infile);
int user_load_problem(sym_environment *env, user_problem *prob);

int main(int argc, char **argv) {
  int termcode;
  sym_environment *env = sym_open_environment();
  if (!env){
    printf("Error initializing environement\n");
    exit(0);
  }
  sym_version();
  user_problem *prob = (user_problem *) calloc(1, sizeof(user_problem));
  CALL_FUNCTION( sym_set_user_data(env, (void *)prob) );
  CALL_FUNCTION( sym_parse_command_line(env, argc, argv) );




  CALL_FUNCTION( user_read_data(prob, prob->par.infile) );
  CALL_FUNCTION( user_load_problem(env, prob) );
  /* CALL_FUNCTION( sym_find_initial_bounds(env) ); */
  CALL_FUNCTION( sym_solve(env) );
  CALL_FUNCTION( sym_close_environment(env) );
  return(0);
}

// here you are supposed to fill user_prob instance.
int user_read_data(user_problem * prob, char * infile) {
  // fill just the conic related part of the problem, this is what we will
  // store in (user_prob *) prob
  int i;
  int j;
  int errors;
  CoinMpsIO mps;
  // read cones from the input file
  int nOfCones = 0;
  int * coneStart = NULL;
  int * coneIdx = NULL;
  int * coneType = NULL;
  // store infile in prob->par->infile
  strncpy(prob->par.infile, infile, MAX_FILE_NAME_LENGTH);
  prob->par.infile[MAX_FILE_NAME_LENGTH-1] = 0;
  mps.messageHandler()->setLogLevel(0);
  mps.setInfinity(mps.getInfinity());
  errors = mps.readMps(infile,"");
  if (!(errors==0 || errors==-3)) {
     return(errors);
  }
  errors = mps.readConicMps(NULL, coneStart, coneIdx, coneType, nOfCones);
  // when there is no conic section errors is -3.
  if (errors==-3) {
    fprintf(stdout, "OsiConic: No conic section is mps file.\n");
  }
  else {
    fprintf(stderr, "OsiConic: readConicMps returned code %d.\n", errors);
    assert (!errors);
  }
  // set number of cones
  prob->num_cones = nOfCones;
  // allocate memory for type, size, members
  prob->cone_type = (int *) malloc(nOfCones*sizeof(int));
  prob->cone_size = (int *) malloc(nOfCones*sizeof(int));
  prob->cone_members = (int **) malloc(nOfCones*sizeof(int*));
  for (i=0; i<nOfCones; ++i) {
    if (coneType[i]!=1 and coneType[i]!=2) {
      fprintf(stderr, "Invalid cone type!\n");
      assert(0);
    }
    int num_members = coneStart[i+1]-coneStart[i];
    if (coneType[i]==2 && num_members<3) {
      fprintf(stderr, "Rotated cones should have at least 3 members!\n");
      assert(0);
    }
    // set size of cone i
    prob->cone_size[i] = num_members;
    // set type for cone i
    if (coneType[i]==1)
      prob->cone_type[i] = 0;
    else
      prob->cone_type[i] = 1;
    // set members of cone i
    prob->cone_members[i] = (int *) malloc(num_members*sizeof(int));
    int k=0;
    for (j=coneStart[i]; j<coneStart[i+1]; ++j) {
      prob->cone_members[i][k] = coneIdx[j];
      k++;
    }
  }
  // check log level and print accordingly
  if (nOfCones) {
    printf("Conic section has %d cones\n",nOfCones);
    for (int iCone=0;iCone<nOfCones;iCone++) {
      printf("Cone %d has %d entries (type %d) ",iCone,coneStart[iCone+1]-coneStart[iCone],
  	      coneType[iCone]);
      for (int j=coneStart[iCone];j<coneStart[iCone+1];j++)
  	 printf("%d ",coneIdx[j]);
      printf("\n");
    }
  }
  free(coneStart);
  free(coneIdx);
  free(coneType);
  if (errors==0 || errors==-3)
    return (FUNCTION_TERMINATED_NORMALLY);
  else
    return (ERROR__USER);
}

// read problem using prob->infile and load it to env
int user_load_problem(sym_environment *env, user_problem *prob) {
  int j;
  int errors;
  int n;
  int m;
  int nz;
  double * obj = 0;
  double * lb = 0;
  double * ub = 0;
  char * is_int = 0;
  int * matbeg = 0;
  char ** colname = 0;
  double * rhs = 0;
  char * sense = 0;
  double * rngval = 0;
  double * matval = 0;
  int * matind = 0;
  double obj_offset = 0.0;
  CoinMpsIO mps;
  mps.messageHandler()->setLogLevel(0);
  mps.setInfinity(mps.getInfinity());
  errors = mps.readMps(prob->par.infile,"");
  // when there is no conic section errors is -3.
  if (errors==-3) {
    fprintf(stdout, "OsiConic: No conic section is mps file.\n");
  }
  else {
    fprintf(stderr, "OsiConic: readConicMps returned code %d.\n", errors);
    assert (!errors);
  }
  //strncpy(probname, mps.getProblemName(), 80);
  m  = mps.getNumRows();
  n  = mps.getNumCols();
  nz = mps.getNumElements();
  prob->rownum = m;
  prob->colnum = n;
  const CoinPackedMatrix * matrixByCol = mps.getMatrixByCol();
  if (n) {
     obj    = (double *) malloc(DSIZE * n);
     ub     = (double *) malloc(DSIZE * n);
     lb     = (double *) malloc(DSIZE * n);
     is_int = (char *) calloc(CSIZE, n);
     memcpy(obj, mps.getObjCoefficients(), DSIZE * n);
     memcpy(ub, mps.getColUpper(), DSIZE * n);
     memcpy(lb, mps.getColLower(), DSIZE * n);
     matbeg = (int *) malloc(ISIZE * (n + 1));
     memcpy(matbeg, matrixByCol->getVectorStarts(), ISIZE * (n + 1));
     colname = (char **) malloc(sizeof(char *) * n);
  }
  if (m) {
     rhs    = (double *) malloc(DSIZE * m);
     sense  = (char *)   malloc(CSIZE * m);
     rngval = (double *) malloc(DSIZE * m);
     memcpy(rhs, mps.getRightHandSide(), DSIZE * m);
     memcpy(sense, mps.getRowSense(), CSIZE * m);
     memcpy(rngval, mps.getRowRange(), DSIZE * m);
  }
  //user defined matind, matval, matbeg--fill as column ordered
  if (nz) {
     matval = (double *) malloc(DSIZE*matbeg[n]);
     matind = (int *)    malloc(ISIZE*matbeg[n]);
     memcpy(matval, matrixByCol->getElements(), DSIZE * matbeg[n]);
     memcpy(matind, matrixByCol->getIndices(), ISIZE * matbeg[n]);
  }
  for (j=0; j<n; j++){
     is_int[j] = mps.isInteger(j);
     colname[j] = (char *) malloc(CSIZE * MAX_NAME_SIZE);
     strncpy(colname[j], mps.columnName(j), MAX_NAME_SIZE);
     colname[j][MAX_NAME_SIZE-1] = 0;
  }
  prob->is_int = is_int;
  // the following is not relevant here, since all mps problems are
  // in minimization
  /* if (obj_sense == SYM_MAXIMIZE){ */
  /*    for (j = 0; j < mip->n; j++){ */
  /* 	 mip->obj[j] *= -1.0; */
  /*    } */
  /* } */
  obj_offset = -mps.objectiveOffset();
  /* Load the problem to SYMPHONY */
  sym_explicit_load_problem(env, n, m, matbeg, matind,
			    matval, lb, ub, is_int, obj, 0, sense,
			    rhs, rngval, true);
  return (FUNCTION_TERMINATED_NORMALLY);
}
