/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Callable         */
/* Library.                                                                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "OsiSymSolverInterface.hpp"
#include "symphony_api.h"

/*===========================================================================*/
/*===========================================================================*/

OsiSymSolverInterface::OsiSymSolverInterface()
{

   env_ = sym_open_environment();

}   
   
/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem()
{

   sym_load_problem(env_);

   setApplicationData((void *) (env_->user));

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::branchAndBound()
{

   sym_solve(env_);

}

/*===========================================================================*/
/*===========================================================================*/

OsiSymSolverInterface::~OsiSymSolverInterface()
{

   sym_close_environment(env_);

   env_ = 0;
   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::reset()
{
   
   sym_close_environment(env_);

   env_ = sym_open_environment();

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setIntParam(OsiIntParam key, int value)
{

   switch(key) {

    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
       env_->par.tm_par.node_limit = value;
       break;

    case OsiLastIntParam:
       return false;

    default:
       return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(OsiSymIntParam key, int value)
{
   switch(key){

    case OsiSymVerbosity:
      env_->par.verbosity =
	 env_->par.tm_par.verbosity =
	 env_->par.lp_par.verbosity =
	 env_->par.cg_par.verbosity = value;
      return true;

    case OsiSymWarmStart:
      env_->par.tm_par.warm_start = value;
      return true;

    case OsiSymNodeLimit:
      env_->par.tm_par.node_limit = value;
      return true;

    case OsiSymFindFirstFeasible:
      env_->par.tm_par.find_first_feasible = value;
      return true;

    case OsiSymUsePermanentCutPools:
      env_->par.use_permanent_cut_pools = value;
      return true;
      
    default: 
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setDblParam(OsiDblParam key, double value)
{
   switch(key){
      
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    case OsiDualTolerance:
    case OsiPrimalTolerance:
    case OsiLastDblParam:
       return false;

    case OsiObjOffset:
      env_->mip->obj_offset = value;
      return true;

    default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(OsiSymDblParam key, double value)
{
   switch(key){

    case OsiSymGranularity:
      env_->par.lp_par.granularity = env_->par.tm_par.granularity = value;
      return true;

    case OsiSymTimeLimit:
      env_->par.tm_par.time_limit = value;
      return true;

    case OsiSymGapLimit:
      env_->par.tm_par.gap_limit = value;
      return true;

   default: 
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setStrParam(OsiStrParam key, 
					const std::string & value)
{

   switch(key) {

    case OsiProbName:
       memcpy(env_->probname, value.c_str(), CSIZE*value.length());
       return true;

    case OsiSolverName:
    case OsiLastStrParam:
       return false;

    default:
       return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(OsiSymStrParam key, 
					   const std::string & value)
{
   switch(key){

   default: 
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getIntParam(OsiIntParam key, int& value) const
{

   switch(key) {
      
    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
       value = env_->par.tm_par.node_limit;
       break;

    case OsiLastIntParam:
       return false;
    default:
       return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(OsiSymIntParam key, int& value)
     const
{
   switch(key){

    case OsiSymVerbosity:
      value = env_->par.verbosity;
      return true;

    case OsiSymWarmStart:
      value = env_->par.tm_par.warm_start;
      return true;

    case OsiSymNodeLimit:
      value = env_->par.tm_par.node_limit;
      return true;

    case OsiSymFindFirstFeasible:
      value = env_->par.tm_par.find_first_feasible;
      return true;
      
    case OsiSymUsePermanentCutPools:
      env_->par.use_permanent_cut_pools = value;
      return true;
      
   default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getDblParam(OsiDblParam key, double& value) const
{

   switch(key){
      
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    case OsiDualTolerance:
    case OsiPrimalTolerance:
    case OsiLastDblParam:
       return false;
       
    case OsiObjOffset:
       value = env_->mip->obj_offset;
      return true;

    default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(OsiSymDblParam key, 
					double& value) const
{
   switch(key){

    case OsiSymGranularity:
      value = env_->par.lp_par.granularity;
      return true;

    case OsiSymTimeLimit:
      value = env_->par.tm_par.time_limit;
      return true;

    case OsiSymGapLimit:
      value = env_->par.tm_par.gap_limit;
      return true;

   default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getStrParam(OsiStrParam key, 
					std::string& value) const
{

   switch(key) {

    case OsiProbName:
       value = env_->probname;
       return true;

    case OsiSolverName:
    case OsiLastStrParam:
       return false;

    default:
       return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(OsiSymStrParam key, 
					   std::string& value) const
{
   switch(key){

   default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setInitialData()
{

   sym_set_defaults(env_);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::parseCommandLine(int argc, char **argv)
{

   sym_parse_command_line(env_, argc, argv);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::findInitialBounds()
{

   sym_find_initial_bounds(env_);

}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::createPermanentCutPools()
{

   return(sym_create_permanent_cut_pools(env_));
   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::initialSolve() //FIX_ME
{

   assert(env_->mip != 0);
   int j;

   char * copyVarTypes = new char[env_->mip->n];

   memcpy(copyVarTypes, env_->mip->is_int, CSIZE * env_->mip->n); 
   
   for (j = 0; j < env_->mip->n; j++){
      env_->mip->is_int[j] = FALSE;
   }
   
   sym_solve(env_); 
   
   memcpy(env_->mip->is_int, copyVarTypes, CSIZE * env_->mip->n); 
   
   delete [] copyVarTypes;
}   

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
					const double* collb, 
					const double* colub, const double* obj,
					const double* rowlb, 
					const double* rowub)
{      
   const double inf = INFINITY;
   
   int nrows = matrix.getNumRows();
   char   * rowSense = new char  [nrows];
   double * rowRhs   = new double[nrows];
   double * rowRange = new double[nrows];
   
   int i;
   for ( i = nrows - 1; i >= 0; --i )
      {
	 const double lower = rowlb ? rowlb[i] : -inf;
	 const double upper = rowub ? rowub[i] : inf;
	 convertBoundToSense( lower, upper, rowSense[i], rowRhs[i], 
			      rowRange[i] );
      }
   
   loadProblem( matrix, collb, colub, obj, rowSense, rowRhs, rowRange ); 
   delete [] rowSense;
   delete [] rowRhs;
   delete [] rowRange;
}  

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
					  double*& collb, double*& colub, 
					  double*& obj, double*& rowlb, 
					  double*& rowub)
					  
{
   loadProblem( *matrix, collb, colub, obj, rowlb, rowub );
   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowlb;  rowlb = 0;
   delete[] rowub;  rowub = 0;   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
					const double* collb, 
					const double* colub,const double* obj,
					const char* rowsen, 
					const double* rowrhs,   
					const double* rowrng)
{

   CoinPackedMatrix * symMatrix;
   bool isColOrdered = true;
   
   if ( !matrix.isColOrdered() ) 
      {
	 symMatrix = new CoinPackedMatrix();
	 symMatrix->reverseOrderedCopyOf(matrix);
	 isColOrdered = false;
      }
   else
      symMatrix = const_cast<CoinPackedMatrix *>(&matrix);

   int nCols = symMatrix->getNumCols();
   int nRows = symMatrix->getNumRows();
   
   if(nCols == 0 || nRows == 0){
      cout<<"loadProblem():The given matrix is empty!"<<endl;
      exit(1);
   }

   const int * matbeg = symMatrix->getVectorStarts();
   const int * matind = symMatrix->getIndices();
   const double * matval = symMatrix->getElements();
   
   loadProblem(nCols,nRows, matbeg, matind, matval, collb, colub, obj, rowsen,
	       rowrhs, rowrng);

   if(!isColOrdered)
      delete symMatrix;
   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
					  double*& collb, double*& colub, 
					  double*& obj, char*& rowsen, 
					  double*& rowrhs, double*& rowrng)
{
   loadProblem( *matrix, collb, colub, obj, rowsen, rowrhs, rowrng );

   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowsen; rowsen = 0;
   delete[] rowrhs; rowrhs = 0;
   delete[] rowrng; rowrng = 0;   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const int numcols, const int numrows,
					const CoinBigIndex * start, 
					const int* index, const double* value,
					const double* collb, 
					const double* colub, const double* obj,
					const double* rowlb, 
					const double* rowub)
{

   if(numcols == 0 || numrows == 0){
      cout<<"loadProblem():The given problem is empty!"<<endl;
      exit(1);
   }

   const double inf = INFINITY;
   
   char   * rowSense = new char  [numrows];
   double * rowRhs   = new double[numrows];
   double * rowRange = new double[numrows];
   
   int i;
   for ( i = numrows - 1; i >= 0; --i )
      {
	 const double lower = rowlb ? rowlb[i] : -inf;
	 const double upper = rowub ? rowub[i] : inf;
	 convertBoundToSense( lower, upper, rowSense[i], rowRhs[i], 
			      rowRange[i] );
      }
   
   loadProblem(numcols,numrows, start, index, value, collb, colub, obj, 
	       rowSense, rowRhs, rowRange);
	          
   delete [] rowSense;
   delete [] rowRhs;
   delete [] rowRange;
   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const int numcols, const int numrows,
					const CoinBigIndex * start, 
					const int* index, const double* value,
					const double* collb,
					const double* colub, const double* obj,
					const char* rowsen, 
					const double* rowrhs, 
					const double* rowrng)   
{

   if(numcols == 0 || numrows == 0){
      cout<<"loadProblem():The given problem is empty!"<<endl;
      exit(1);
   }

   //Assuming all the pointers are always given NOT null, except rowrng!

   double * rowRange = new double[numrows];
   double t =0;
   int j;
   (void)used_time(&t);

   if(rowrng)
      rowRange = const_cast<double *> (rowrng);
   else
      CoinFillN(rowRange, numrows, 0.0);
   
   env_->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
 
   env_->mip->m  = numrows;
   env_->mip->n  = numcols;
   env_->mip->nz = start[numcols];
   
   env_->mip->obj    = (double *) malloc(DSIZE * numcols);
   env_->mip->rhs    = (double *) malloc(DSIZE * numrows);
   env_->mip->sense  = (char *)   malloc(CSIZE * numrows);
   env_->mip->rngval = (double *) malloc(DSIZE * numrows);
   env_->mip->ub     = (double *) malloc(DSIZE * numcols);
   env_->mip->lb     = (double *) malloc(DSIZE * numcols);
   env_->mip->is_int = (char *)   calloc(CSIZE, numcols);
   
   memcpy(env_->mip->obj, const_cast <double *> (obj), DSIZE * numcols); 
   memcpy(env_->mip->rhs, const_cast <double *> (rowrhs), DSIZE * numrows); 
   memcpy(env_->mip->sense, const_cast <char *> (rowsen), CSIZE * numrows); 
   memcpy(env_->mip->rngval, rowRange, DSIZE * numrows); 	  
   memcpy(env_->mip->ub, const_cast <double *> (colub), DSIZE * numcols); 
   memcpy(env_->mip->lb, const_cast <double *> (collb), DSIZE * numcols); 
	     
   //user defined matind, matval, matbeg--fill as column ordered
   
   env_->mip->matbeg = (int *) malloc(ISIZE * numcols + 1);
   env_->mip->matval = (double *) malloc(DSIZE*start[numcols]);
   env_->mip->matind = (int *)    malloc(ISIZE*start[numcols]);
   
   memcpy(env_->mip->matbeg, const_cast<int *>(start), ISIZE * (numcols + 1));
   memcpy(env_->mip->matval, const_cast<double *> (value), DSIZE * 
	  start[numcols]);  
   memcpy(env_->mip->matind, const_cast<int *> (index), ISIZE * 
	  start[numcols]);  

#if 0   //FIX_ME may be having a funct. called setColNames()???
   mip->colname = (char **) malloc(sizeof(char *) * mip->n);   
   
   for (j = 0; j < mip->n; j++){
      mip->colname[j] = (char *) malloc(CSIZE * 9);
      strncpy(mip->colname[j], const_cast<char*>(mps.columnName(j)), 9);
      mip->colname[j][8] = 0;
   }
#endif
  
   /* Start up the graphics window*/
#ifndef WIN32
   init_draw_graph_u(env_);   
#endif

   /*------------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*------------------------------------------------------------------------*/

   initialize_root_node_u(env_); 

   env_->comp_times.readtime = used_time(&t);

 
   delete [] rowRange;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isAbandoned() const 
{
   switch(env_->termcode){

    case SOMETHING_DIED:
    case TM_ERROR__NUMERICAL_INSTABILITY:    
       return true;

    default:
       break;
   }

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isProvenOptimal() const
{
 
   switch(env_->termcode){

    case TM_OPTIMAL_SOLUTION_FOUND:
       return true;

    default:
       break;
   }
   return false;
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isProvenPrimalInfeasible() const
{
   switch(env_->termcode){

    case TM_OPTIMAL_SOLUTION_FOUND:
    case TM_FOUND_FIRST_FEASIBLE:
       return true;

    default:
       break;
   }
   return false;
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isPrimalObjectiveLimitReached() const
{
   switch(env_->termcode){

    case TM_OPTIMAL_SOLUTION_FOUND:
    case TM_FOUND_FIRST_FEASIBLE:
       return true;

    default:
       break;
   }
   return false;
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isIterationLimitReached() const
{
   switch(env_->termcode){

    case TM_NODE_LIMIT_EXCEEDED:
       return true;

    default:
       break;
   }
   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isTimeLimitReached() const
{
   switch(env_->termcode){

    case TM_TIME_LIMIT_EXCEEDED:
       return true;

    default:
       break;
   }
   return false;
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isTargetGapReached() const
{
   switch(env_->termcode){

    case TM_TARGET_GAP_ACHIEVED:
       return true;

    default:
       break;
   }
   return false;
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface:: getNumCols() const 
{

   assert (env_->mip);

   return env_->mip->n;

}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getNumRows() const
{
   assert (env_->mip);
   return env_->mip->m;
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getNumElements() const
{
   assert (env_->mip);
   return env_->mip->nz;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColLower() const
{
   assert (env_->mip);
   return env_->mip->lb;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColUpper() const
{
   assert (env_->mip);
   return env_->mip->ub;
}

/*===========================================================================*/
/*===========================================================================*/

const char * OsiSymSolverInterface::getRowSense() const
{
   assert (env_->mip);
   return env_->mip->sense;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRightHandSide() const
{
   assert (env_->mip);
   return env_->mip->rhs;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowRange() const
{
   assert (env_->mip);
   return env_->mip->rngval;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowLower() const
{
   assert (env_->mip);
   
   double * lower = new double[env_->mip->m];
   double upper;
   int i;

   for ( i = env_->mip->m - 1; i >= 0; --i )
      {
	 convertBoundToSense( lower[i], upper, env_->mip->sense[i], 
			      env_->mip->rhs[i], env_->mip->rngval[i] ); 
      }

   return lower;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowUpper() const
{
   assert (env_->mip);
   
   double * upper = new double[env_->mip->m];
   double lower;
   int i;

   for ( i = env_->mip->m - 1; i >= 0; --i )
      {
	 convertBoundToSense( lower, upper[i], env_->mip->sense[i], 
			      env_->mip->rhs[i], env_->mip->rngval[i] ); 
      }

   return upper;
}

/*===========================================================================*/
/*===========================================================================*/
  
const double * OsiSymSolverInterface::getObjCoefficients() const
{
   assert (env_->mip);
   return env_->mip->obj;

}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getObjSense() const
{
   assert (env_->mip);
   
   if(env_->mip->obj_sense == MINIMIZE)
      return 1;
   else if(env_->mip->obj_sense == MAXIMIZE)
      return -1;
   else
      return 1; //if it was not set, then, behave as having a min problem.
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isContinuous(int colIndex) const
{

   assert (env_->mip);
      
   if(!env_->mip->is_int[colIndex])
      return true;
   else
      return false;
}      

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isBinary(int colIndex) const
{
   assert (env_->mip);
   
   if(env_->mip->is_int[colIndex] && env_->mip->lb[colIndex] == 0.0 &&
      env_->mip->ub[colIndex] == 1.0)
      return true;
   else
      return false;
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isInteger(int colIndex) const
{
   assert (env_->mip);
   
   if(env_->mip->is_int[colIndex])
      return true;
   else
      return false;
}      

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isIntegerNonBinary(int colIndex) const
{

   if(!isBinary(colIndex) && isInteger(colIndex))
      return true;
   else
      return false;
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isFreeBinary(int colIndex) const //FIX_ME
{
   assert (env_->mip);
   
   if(env_->mip->is_int[colIndex] && env_->mip->lb[colIndex] == 0.0 &&
      env_->mip->ub[colIndex] == 1.0)
      return true;
   else
      return false;

}

/*===========================================================================*/
/*===========================================================================*/

const CoinPackedMatrix * OsiSymSolverInterface::getMatrixByRow() const
{
   assert (env_->mip);

   CoinPackedMatrix * rowMatrix = 
      const_cast<CoinPackedMatrix *>(getMatrixByCol()); 
   
   rowMatrix->reverseOrdering();

   return rowMatrix;

}
/*===========================================================================*/
/*===========================================================================*/

const CoinPackedMatrix * OsiSymSolverInterface::getMatrixByCol() const
{
   assert (env_->mip);
   
   CoinPackedMatrix * colMatrix = 
      new CoinPackedMatrix(true, env_->mip->m, env_->mip->n, env_->mip->nz,
			   env_->mip->matval, env_->mip->matind, 
			   env_->mip->matbeg, 0);
   return colMatrix;
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getInfinity() const
{
   /* FIXME; Make the use of INFINITY consistent */
   return INFINITY;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColSolution() const
{

   int i;
   double * colSol;
   lp_sol sol = env_->best_sol;

   if(!sol.xlength){
      printf("\nNo Solution Found\n\n");
      return 0;
   }
   else{
      colSol = new double[env_->mip->n];
      CoinFillN(colSol, env_->mip->n, 0.0);
      for( i = 0; i<sol.xlength; i++){
	 colSol[sol.xind[i]] = sol.xval[i];
      }
      return (colSol);
   }
}

/*===========================================================================*/
/*===========================================================================*/

 /* FIX_ME, will return the rowAct of the root_ rows all of which may not be 
    binding
 */					 
const double * OsiSymSolverInterface::getRowActivity() const
{

   assert(env_->mip);

   double * rowAct = new double [env_->mip->m];
   CoinFillN(rowAct, env_->mip->m, 0.0);
   const double * colSol = getColSolution();  
   int i, j;

   if(colSol){

      const CoinPackedMatrix * rowMatrix = getMatrixByRow();
      const int * matBeg = rowMatrix->getVectorStarts();
      const double * matVal = rowMatrix->getElements();
      const int * matInd = rowMatrix->getIndices();
      for(i = 0; i<env_->mip->m; i++){
	 for(j = matBeg[i]; j<matBeg[i+1]; j++){
	    rowAct[i] += matVal[j] * colSol[matInd[j]];
	 }
      }
      return rowAct;
   }
   else
      return 0;
         
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getObjValue() const
{

   if(env_->best_sol.xlength)
      return (env_->best_sol.objval);
   else{
      printf("\nNo Solution Found\n\n");
      return INFINITY;
   }
}

/*===========================================================================*/
/*===========================================================================*/

/* returns the number of analyzed nodes */
int OsiSymSolverInterface::getIterationCount() const
{
   return env_->warm_start->stat.analyzed;  
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObjCoeff( int elementIndex, 
					 double elementValue )
{

   assert(env_->mip);

   env_->mip->obj[elementIndex] = elementValue;
   //   env_->desc_modification = OBJ_COEF_CHANGED;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColLower( int elementIndex, 
					 double elementValue )
{
   assert(env_->mip);

   env_->mip->lb[elementIndex] = elementValue;
   //   env_->desc_modification = COL_BOUNDS_CHANGED;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColUpper( int elementIndex, 
					 double elementValue )
{
   assert(env_->mip);

   env_->mip->ub[elementIndex] = elementValue;
   //   env_->desc_modification = COL_BOUNDS_CHANGED;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowLower( int elementIndex, 
					 double elementValue )
{

   double rhs   = getRightHandSide()[elementIndex];
   double range = getRowRange()[elementIndex];
   char   sense = getRowSense()[elementIndex];
   double lower, upper;
   
   convertSenseToBound( sense, rhs, range, lower, upper );
   if( lower != elementValue ) {
      convertBoundToSense( elementValue, upper, sense, rhs, range );
      setRowType( elementIndex, sense, rhs, range );
   }
   //   env_->desc_modification = ROW_TYPE_CHANGED;
}   
   
/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowUpper( int elementIndex, 
					 double elementValue )
{

   double rhs   = getRightHandSide()[elementIndex];
   double range = getRowRange()[elementIndex];
   char   sense = getRowSense()[elementIndex];
   double lower, upper;
   
   convertSenseToBound( sense, rhs, range, lower, upper );
   if( upper != elementValue ) {
      convertBoundToSense( lower, elementValue, sense, rhs, range );
      setRowType(elementIndex, sense, rhs, range );
   }
}   

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowType(int index, char sense, 
				       double rightHandSide, double range)
{

   /* assuming row_lb = rightHandSide - range */
   assert(env_->mip);

   env_->mip->sense[index] = sense;   
   env_->mip->rhs[index] = rightHandSide;
   env_->mip->rngval[index] = range;

   //   env_->desc_modification = ROW_TYPE_CHANGED;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObjSense(double s)
{
   assert(env_->mip);

   if(s == 1)
      env_->mip->obj_sense = MINIMIZE;
   else if (s==-1)
      env_->mip->obj_sense = MAXIMIZE;
   else
      /* assume it to be min problem */
      env_->mip->obj_sense = MINIMIZE;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColSolution(const double *colsol)
{

   assert(env_->mip);   
   assert(colsol);

   int i, j, nz =0;
   double value; 
   bool feasible;

   /* Feasibility Check*/

   /* step 1. check for bounds and integrality */   
   for (i = env_->mip->n - 1; i >= 0; i--){
      if(colsol[i] < env_->mip->lb[i] || colsol[i] > env_->mip->ub[i])
	 break;
      if(colsol[i] !=0.0)
	 nz++;
      if (!env_->mip->is_int[i])
	 continue; /* Not an integer variable */
      value = colsol[i];
      if (colsol[i] > env_->mip->lb[i] && colsol[i] < env_->mip->ub[i]
	  && colsol[i]-floor(colsol[i]) > env_->par.lp_par.granularity &&
	  ceil(colsol[i])-colsol[i] > env_->par.lp_par.granularity){
	 break;   //FIX_ME, can we use granularity here?
      }
   }

   feasible = i < 0 ? true : false;
   
   /* step 2. check for the constraint matrix */
   
   if(feasible){      
      double * rowAct = new double [env_->mip->m];
      CoinFillN(rowAct, env_->mip->m, 0.0);
      int i, j;
      const CoinPackedMatrix * rowMatrix = getMatrixByRow();
      const int * matBeg = rowMatrix->getVectorStarts();
      const double * matVal = rowMatrix->getElements();
      const int * matInd = rowMatrix->getIndices();
      for(i = 0; i<env_->mip->m; i++){
	 for(j = matBeg[i]; j<matBeg[i+1]; j++){
	    rowAct[i] += matVal[j] * colsol[matInd[j]];
	 }
	 
	 switch(env_->mip->sense[i]){
	  case 'L': 
	     if(rowAct[i] > env_->mip->rhs[i])
		feasible = false;
	     break;
	  case 'G':
	     if(rowAct[i] < env_->mip->rhs[i])
		feasible = false;
	     break;
	  case 'E':
	     if(rowAct[i] != env_->mip->rhs[i])
		feasible = false;
	     break;
	  case 'R':
	     if(rowAct[i] > env_->mip->rhs[i] || 
		rowAct[i] < env_->mip->rhs[i] - env_->mip->rngval[i])
		feasible = false;
	     break;
	  case 'N':
	  default:
	     break;
	 }

	 if(!feasible) 
	    break;
      }
   }

   if(feasible){
      /* now, it is feasible, set the best_sol to colsol */
      //FIX_ME
      
      lp_sol * sol = new lp_sol;      
      sol->xlength = nz;
      sol->xind = new int[nz];
      sol->xval = new double[nz];
      
      for(i = 0, j =0 ; i<env_->mip->n; i++){
	 if(colsol[i] != 0.0){
	    sol->xind[j] = i;
	    sol->xval[j] = colsol[i];
	    sol->objval += colsol[i] * env_->mip->obj[i];	   
	    j++;
	 }      
      }

      if(env_->has_ub_estimate){
	 if(env_->ub_estimate > sol->objval)
	    env_->ub_estimate = sol->objval; //no need for this, I guess.
      }
      else{
	 env_->has_ub_estimate = TRUE;
	 env_->ub_estimate = sol->objval; //no need for this, I guess.
      }

      if(env_->has_ub){
	 if(env_->ub > sol->objval)
	    env_->ub = sol->objval;
      }
      else{
	 env_->has_ub = TRUE;
	 env_->ub = sol->objval;
      }
      
      if(env_->best_sol.xlength){
	 if(env_->best_sol.objval > sol->objval)
	    env_->best_sol = *sol;
      }
      else
	 env_->best_sol = *sol;     
   }
   else{
      cout<<"The given colSol is not feasible..."<<endl;
   }  
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setContinuous(int index)
{
   env_->mip->is_int[index] = FALSE;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setInteger(int index)
{
   env_->mip->is_int[index] = TRUE;
}
/*===========================================================================*/
/*===========================================================================*/
void OsiSymSolverInterface::addCol(const CoinPackedVectorBase& vec,
				   const double collb, const double colub,   
				   const double obj)
{

   CoinPackedMatrix * colMat = 
      const_cast<CoinPackedMatrix*> (getMatrixByCol());
   colMat->appendCol(vec);

   int n = colMat->getNumCols(), m = colMat->getNumRows();
   int nz=colMat->getNumElements();
   
   assert(env_->mip->m == m);

   int * matBeg = new int[n + 1];
   int * matInd = new int[nz];   
   double * matVal = new double[nz];
   double * colLow = new double[n];
   double * colUp = new double[n];
   double * objN = new double[n];
   char * isInt = new char[n];


   memcpy(matBeg, const_cast<int *>(colMat->getVectorStarts()), 
	  ISIZE * (n + 1));
   memcpy(matVal, const_cast<double *> (colMat->getElements()), 
	  DSIZE * nz);
   memcpy(matInd, const_cast<int *> (colMat->getIndices()), 
	  ISIZE * nz);

   memcpy(colLow, env_->mip->lb, DSIZE * n-1);
   memcpy(colUp, env_->mip->ub, DSIZE * n-1);
   memcpy(objN, env_->mip->obj, DSIZE * n-1);
   memcpy(isInt, env_->mip->is_int, CSIZE * n-1);

   colLow[n-1] = collb;
   colUp[n-1] = colub;
   objN[n-1] = obj;

   delete [] env_->mip->matbeg;
   delete [] env_->mip->matval;
   delete [] env_->mip->matind;

   delete [] env_->mip->lb;
   delete [] env_->mip->ub;
   delete [] env_->mip->obj;
   delete [] env_->mip->is_int;

   env_->mip->n = n;
   env_->mip->nz = nz;

   env_->mip->matbeg = matBeg;
   env_->mip->matind = matInd;
   env_->mip->matval = matVal;

   env_->mip->lb = colLow;
   env_->mip->ub = colUp;
   env_->mip->obj = objN;
   env_->mip->is_int = isInt;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::addRow(const CoinPackedVectorBase& vec,
				   const double rowlb, const double rowub)
{
   char    rowSen;
   double  rowRhs;
   double  rowRange;
   
   convertBoundToSense( rowlb, rowub, rowSen, rowRhs, rowRange);
   addRow(vec, rowSen, rowRhs, rowRange);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::addRow(const CoinPackedVectorBase& vec,
				   const char rowsen, const double rowrhs,   
				   const double rowrng)
{

   CoinPackedMatrix * colMat = 
      const_cast<CoinPackedMatrix*> (getMatrixByCol());   
   colMat->appendRow(vec);

   int n = colMat->getNumCols(), m = colMat->getNumRows();
   int nz=colMat->getNumElements();

   assert(env_->mip->n == n);
   
   int * matBeg = new int[n + 1];
   int * matInd = new int[nz];   
   double * matVal = new double[nz];

   char * rowSense = new char[m];
   double * rowRhs = new double[m];
   double * rowRange = new double[m];

   memcpy(matBeg, const_cast<int *>(colMat->getVectorStarts()), 
	  ISIZE * (n + 1));
   memcpy(matVal, const_cast<double *> (colMat->getElements()), 
	  DSIZE * nz);
   memcpy(matInd, const_cast<int *> (colMat->getIndices()), 
	  ISIZE * nz);

   memcpy(rowSense, env_->mip->sense, CSIZE * m-1);
   memcpy(rowRhs, env_->mip->rhs, DSIZE * m-1);
   memcpy(rowRange, env_->mip->rngval, DSIZE * m-1);

   rowSense[m-1] = rowsen;
   rowRhs[m-1] = rowrhs;
   rowRange[m-1] = rowrng;

   delete [] env_->mip->matbeg;
   delete [] env_->mip->matval;
   delete [] env_->mip->matind;

   delete [] env_->mip->sense;
   delete [] env_->mip->rhs;
   delete [] env_->mip->rngval;

   env_->mip->m = m;
   env_->mip->nz = nz;

   env_->mip->matbeg = matBeg;
   env_->mip->matind = matInd;
   env_->mip->matval = matVal;

   env_->mip->sense = rowSense;
   env_->mip->rhs = rowRhs;
   env_->mip->rngval = rowRange;
   
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::deleteCols(const int num, const int * colIndices)
{

   CoinPackedMatrix * colMat = 
      const_cast<CoinPackedMatrix*> (getMatrixByCol());   
   colMat->deleteCols(num, colIndices);
   
   int n = colMat->getNumCols(), m = colMat->getNumRows(); 
   int nz=colMat->getNumElements();
   int i, j, k;

   assert(env_->mip->m == m);
   assert(env_->mip->n == n + num);

   int * matBeg = new int[n + 1];
   int * matInd = new int[nz];   
   double * matVal = new double[nz];
 
   double * colLow = new double[n];
   double * colUp = new double[n];
   double * objN = new double[n];
   char * isInt = new char[n];

   memcpy(matBeg, const_cast<int *>(colMat->getVectorStarts()), 
	  ISIZE * (n + 1));
   memcpy(matVal, const_cast<double *> (colMat->getElements()), 
	  DSIZE * nz);
   memcpy(matInd, const_cast<int *> (colMat->getIndices()), 
	  ISIZE * nz);


   for(i = 0, j = 0, k = 0; i < env_->mip->n; i++){
      if(i != colIndices[j]){
	 colLow[k] = env_->mip->lb[i];
	 colUp[k] = env_->mip->ub[i];
	 objN[k] = env_->mip->obj[i];
	 isInt[k] = env_->mip->is_int[i];
	 k++;
      }
      else{
	 j++;
      }
   }
      
   if(k != n){
      cout<<"Unknown problem in deleteCols()!"<<endl;
      exit(1);
   }      

   delete [] env_->mip->matbeg;
   delete [] env_->mip->matval;
   delete [] env_->mip->matind;

   delete [] env_->mip->lb;
   delete [] env_->mip->ub;
   delete [] env_->mip->obj;
   delete [] env_->mip->is_int;

   env_->mip->n = n;
   env_->mip->nz = nz;

   env_->mip->matbeg = matBeg;
   env_->mip->matind = matInd;
   env_->mip->matval = matVal;

   env_->mip->lb = colLow;
   env_->mip->ub = colUp;
   env_->mip->obj = objN;  
   env_->mip->is_int = isInt;

}   

/*===========================================================================*/
/*===========================================================================*/
void OsiSymSolverInterface::deleteRows(const int num, const int * rowIndices)
{

   CoinPackedMatrix * colMat = const_cast<CoinPackedMatrix*> (getMatrixByCol());   
   colMat->deleteRows(num, rowIndices);
   
   int n = colMat->getNumCols(), m = colMat->getNumRows();
   int nz=colMat->getNumElements();
   int i, j, k;

   assert(env_->mip->n == n);
   assert(env_->mip->m == m + num);

   int * matBeg = new int[n + 1];
   int * matInd = new int[nz];   
   double * matVal = new double[nz];
 
   char * rowSense = new char[m];
   double * rowRhs = new double[m];
   double * rowRange = new double[m];

   memcpy(matBeg, const_cast<int *>(colMat->getVectorStarts()), 
	  ISIZE * (n + 1));
   memcpy(matVal, const_cast<double *> (colMat->getElements()), 
	  DSIZE * nz);
   memcpy(matInd, const_cast<int *> (colMat->getIndices()), 
	  ISIZE * nz);

   for(i = 0, j = 0, k = 0; i < env_->mip->m; i++){
      if(i != rowIndices[j]){
	 rowSense[k] = env_->mip->sense[i];
	 rowRhs[k] = env_->mip->rhs[i];
	 rowRange[k] = env_->mip->rngval[i];
	 k++;
      }
      else{
	 j++;
      }
   }

   if(k != m){
      cout<<"Unknown problem in deleteRows()!"<<endl;
      exit(1);
   }      

   delete [] env_->mip->matbeg;
   delete [] env_->mip->matval;
   delete [] env_->mip->matind;

   delete [] env_->mip->sense;
   delete [] env_->mip->rhs;
   delete [] env_->mip->rngval;

   env_->mip->m = m;
   env_->mip->nz = nz;

   env_->mip->matbeg = matBeg;
   env_->mip->matind = matInd;
   env_->mip->matval = matVal;

   env_->mip->sense = rowSense;
   env_->mip->rhs = rowRhs;
   env_->mip->rngval = rowRange;
   
}   

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::writeMps(const char *filename,
				     const char *extension,
				     double objSense) const
{

   const CoinPackedMatrix * colMat = getMatrixByCol();
      
   CoinMpsIO mps;
   mps.setMpsData(*colMat, INFINITY, env_->mip->lb, env_->mip->ub, 
		 env_->mip->obj, env_->mip->is_int, 
		 env_->mip->sense, env_->mip->rhs,
		 env_->mip->rngval, NULL, NULL);

   string f(filename);
   string e(extension);
   string fullname = f + "." + e;

   mps.writeMps(fullname.c_str());

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::applyRowCut( const OsiRowCut & rc) 
{
   double lb, ub;
   CoinPackedVector rowVec;

   rowVec = rc.row();
   lb = rc.lb();
   ub = rc.ub();

   addRow(rowVec, lb, ub);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::applyColCut( const OsiColCut & cc) 
{

   /* assuming the given bounds are feasible */ //FIX_ME

   int i;
   const CoinPackedVector & lbs = cc.lbs();
   const CoinPackedVector & ubs = cc.ubs();

   const int * colInd = lbs.getIndices();
   const double * colVal = lbs.getElements();
   
   for(i = 0; i<lbs.getNumElements(); i++){
      env_->mip->lb[colInd[i]] = colVal[i];
   }		    

   colInd = ubs.getIndices();
   colVal = ubs.getElements();

   for(i = 0; i<ubs.getNumElements(); i++){
      env_->mip->ub[colInd[i]] = colVal[i];
   }		    
}

/*===========================================================================*/
/*===========================================================================*/
CoinWarmStart * OsiSymSolverInterface::getWarmStart() const
{

#ifdef COMPILE_IN_TM

   /* FIX_ME! Add an SymIntParam to determine whether keep the tree
      in case a getWarmStart() is called or not!
   */

   SymWarmStart * symWS = new SymWarmStart(env_->warm_start, false);
   return symWS;
   
#else
   return 0;
#endif
}   

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{

   CoinWarmStart * wsC = const_cast<CoinWarmStart*> (warmstart);
   SymWarmStart *symWS = dynamic_cast<SymWarmStart*>(wsC);
   
   warm_start_desc * ws  = 
      const_cast<warm_start_desc*>(symWS->getWarmStartDesc());
   if(!ws){
      cout<<"setWarmStart(): An empty warmstart was given!"<<endl;
      return false;
   }
   
#if 0
   SymWarmStart *copySymWS = symWS->clone();
   env_->warm_start = 
      const_cast<warm_start_desc*>(copySymWS->getWarmStartDesc());
#endif

   env_->warm_start = ws;
   env_->par.tm_par.warm_start = TRUE;   
   return true;   
}

/*===========================================================================*/
/*===========================================================================*/
