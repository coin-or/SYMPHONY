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

void OsiSymSolverInterface::resolve()
{

   sym_resolve(env_);

}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getLbForNewRhs(int cnt, int *index, 
						  double * value)
{
   return sym_get_lb_for_new_rhs(env_, cnt, index, value);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::multiCriteriaBranchAndBound()
{

   sym_mc_solve(env_);

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
       sym_set_int_param(env_, "node limit", value);
       return true;
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
       sym_set_int_param(env_, "verbosity", value);
       return true;

    case OsiSymWarmStart:
       sym_set_int_param(env_, "warm_start", value);
       return true;
       
    case OsiSymNodeLimit:
       sym_set_int_param(env_, "node_limit", value);
       return true;
       
    case OsiSymFindFirstFeasible:
       sym_set_int_param(env_, "find_first_feasible", value);
       return true;
       
   case OsiSymSearchStrategy:
      sym_set_int_param(env_, "node_selection_strategy", value);
      return true;

    case OsiSymUsePermanentCutPools:
       sym_set_int_param(env_, "use_permanent_cut_pools", value);
       return true;
 
    case OsiSymKeepDescOfPruned:
       sym_set_int_param(env_, "keep_description_of_pruned", value);
       return true;

    case OsiSymMultiCriteriaFindNondominatedSolutions:
       sym_set_int_param(env_, "mc_find_nondominated_solutions", value);
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
       sym_set_dbl_param(env_, "granularity", value);
       return true;
       
    case OsiSymTimeLimit:
       sym_set_dbl_param(env_, "time_limit", value);
       return true;
       
    case OsiSymGapLimit:
       sym_set_dbl_param(env_, "gap_limit", value);
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
       sym_set_str_param(env_, "problem_name", 
			     const_cast<char*>(value.c_str()));
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
       value = sym_get_int_param(env_, "node_limit");
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
      value = sym_get_int_param(env_, "verbosity");
      return true;

    case OsiSymWarmStart:
      value = sym_get_int_param(env_, "warm_start");
      return true;

    case OsiSymNodeLimit:
      value = sym_get_int_param(env_, "node_limit");
      return true;

    case OsiSymFindFirstFeasible:
      value = sym_get_int_param(env_, "find_first_feasible");
      return true;
      
    case OsiSymUsePermanentCutPools:
       value = sym_get_int_param(env_, "use_permanent_cut_pools");
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
      value = sym_get_dbl_param(env_, "granularity");
      return true;

    case OsiSymTimeLimit:
      value = sym_get_dbl_param(env_, "time_limit");
      return true;

    case OsiSymGapLimit:
      value = sym_get_dbl_param(env_, "gap_limit");
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
       value = sym_get_str_param(env_, "problem_name");
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
	 const double lower = rowlb[i] ? rowlb[i] : -inf;
	 const double upper = rowub[i] ? rowub[i] : inf;
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

   /* Assuming no extra gap is put in given matrix */
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
   
   char   * sense = new char  [numrows];
   double * rhs   = new double[numrows];
   double * range = new double[numrows];
   
   int i;
   for ( i = numrows - 1; i >= 0; --i ){
      const double lower = rowlb ? rowlb[i] : -inf;
      const double upper = rowub ? rowub[i] : inf;
      
      /* convertBountToSense */
      range[i] = 0.0;
      if (lower > -inf) {
	 if (upper < inf) {
	    rhs[i] = upper;
	    if (upper==lower) {
	       sense[i] = 'E';
	    } else {
	       sense[i] = 'R';
	       range[i] = upper - lower;
	    }
	 } else {
	    sense[i] = 'G';
	    rhs[i] = lower;
	 }
      } else {
	 if (upper < inf) {
	    sense[i] = 'L';
	    rhs[i] = upper;
	 } else {
	    sense[i] = 'N';
	    rhs[i] = 0.0;
	 }
      }    
      /*	 convertBoundToSense( lower, upper, rowSense[i], rowRhs[i], 
		 rowRange[i] ); */
   }

   loadProblem(numcols,numrows, start, index, value, collb, colub, obj, 
	       sense, rhs, range);
	          
   delete [] sense;
   delete [] rhs;
   delete [] range;
   
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

   sym_load_problem_user(env_, numcols, numrows, const_cast<int*>(start), 
			 const_cast<int*>(index), const_cast<double*>(value), 
			 const_cast<double*>(collb), 
			 const_cast<double*>(colub), const_cast<double*>(obj), 
			 const_cast<char*>(rowsen), 
			 const_cast<double*>(rowrhs), 
			 const_cast<double*>(rowrng)); 
			 
   setApplicationData((void *) (env_->user));
    
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isAbandoned() const 
{

   if(sym_is_abandoned(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isProvenOptimal() const
{
   if(sym_is_proven_optimal(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isProvenPrimalInfeasible() const
{
   if(sym_is_proven_primal_infeasible(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isPrimalObjectiveLimitReached() const
{
   if(sym_is_primal_objective_limit_reached(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isIterationLimitReached() const
{
   if(sym_is_iteration_limit_reached(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isTimeLimitReached() const
{
   if(sym_is_time_limit_reached(env_)){
      return true;
   }
   else{
      return false;
   }
}
/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isTargetGapReached() const
{
   if(sym_is_target_gap_achieved(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface:: getNumCols() const 
{
   return sym_get_num_cols(env_);
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getNumRows() const
{
   return sym_get_num_rows(env_);
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getNumElements() const
{
   return sym_get_num_elements(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColLower() const
{
   return sym_get_col_lower(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColUpper() const
{
   return sym_get_col_upper(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const char * OsiSymSolverInterface::getRowSense() const
{
   return sym_get_row_sense(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRightHandSide() const
{

   return sym_get_rhs(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowRange() const
{

   return sym_get_row_range(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowLower() const
{
   return sym_get_row_lower(env_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowUpper() const
{
   return sym_get_row_upper(env_);
}

/*===========================================================================*/
/*===========================================================================*/
  
const double * OsiSymSolverInterface::getObjCoefficients() const
{
   return sym_get_obj_coeff(env_);
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getObjSense() const
{
   return (double)sym_get_obj_sense(env_);
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isContinuous(int colIndex) const
{
   if(sym_is_continuous(env_, colIndex)){
      return true;
   }
   else{
      return false;
   }
}      

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isBinary(int colIndex) const
{
   if(sym_is_binary(env_, colIndex)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isInteger(int colIndex) const
{
   if(sym_is_integer(env_, colIndex)){
      return true;
   }
   else{
      return false;
   }
}      

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isIntegerNonBinary(int colIndex) const
{

   if(!isBinary(colIndex) && isInteger(colIndex)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isFreeBinary(int colIndex) const //FIX_ME
{

   if(isBinary(colIndex)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

const CoinPackedMatrix * OsiSymSolverInterface::getMatrixByRow() const
{

   CoinPackedMatrix * rowMatrix = 
      const_cast<CoinPackedMatrix *>(getMatrixByCol()); 
   
   rowMatrix->reverseOrdering();

   return rowMatrix;

}
/*===========================================================================*/
/*===========================================================================*/

const CoinPackedMatrix * OsiSymSolverInterface::getMatrixByCol() const
{

   assert(env_->mip);
   
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
   return sym_get_infinity();
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColSolution() const
{
   return (sym_get_col_solution(env_));
}

/*===========================================================================*/
/*===========================================================================*/

 /* FIX_ME, will return the rowAct of the root_ rows all of which may not be 
    binding
 */					 
const double * OsiSymSolverInterface::getRowActivity() const
{
   return (sym_get_row_activity(env_));         
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getObjValue() const
{
   return(sym_get_obj_val(env_));
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getIterationCount() const
{
   return sym_get_iteration_count(env_);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObjCoeff( int elementIndex, 
					 double elementValue )
{

   sym_set_obj_coeff(env_, elementIndex, elementValue);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObj2Coeff( int elementIndex, 
					 double elementValue )
{

   sym_set_obj2_coeff(env_, elementIndex, elementValue);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColLower( int elementIndex, 
					 double elementValue )
{
   sym_set_col_lower(env_, elementIndex, elementValue);
   //   env_->desc_modification = COL_BOUNDS_CHANGED;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColUpper( int elementIndex, 
					 double elementValue )
{
   sym_set_col_upper(env_, elementIndex, elementValue);
   //   env_->desc_modification = COL_BOUNDS_CHANGED;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowLower( int elementIndex, 
					 double elementValue )
{
   sym_set_row_lower(env_, elementIndex, elementValue);
}   
   
/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowUpper( int elementIndex, 
					 double elementValue )
{
   sym_set_row_upper(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowType(int index, char sense, 
				       double rightHandSide, double range)
{
   sym_set_row_type(env_, index, sense, rightHandSide, range);

   //   env_->desc_modification = ROW_TYPE_CHANGED;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObjSense(double s)
{
   sym_set_obj_sense(env_, (int)s);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColSolution(const double *colsol)
{
   sym_set_col_solution(env_, const_cast<double*>(colsol));
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setContinuous(int index)
{
   sym_set_continuous(env_, index);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setInteger(int index)
{
   sym_set_integer(env_, index);
}
/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColName(char **colname)
{
   sym_set_col_names (env_, colname);
}

/*=====================x======================================================*/
/*===========================================================================*/
void OsiSymSolverInterface::addCol(const CoinPackedVectorBase& vec,
				   const double collb, const double colub,   
				   const double obj)
{
   int numElements = vec.getNumElements();
   int * indices = const_cast<int*>(vec.getIndices());
   double * elements = const_cast<double*>(vec.getElements());
   
   sym_add_col(env_, numElements, indices, elements, collb, colub, obj, NULL);
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
   int i, j, m, n, nz;
   int numElements = vec.getNumElements();
   int * indices = const_cast<int*>(vec.getIndices());
   double * elements = const_cast<double*>(vec.getElements());

   sym_add_row(env_, numElements, indices, elements, rowsen, rowrhs, rowrng);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::deleteCols(const int num, const int * colIndices)
{
   sym_delete_cols(env_, num, const_cast<int*>(colIndices));
}   

/*===========================================================================*/
/*===========================================================================*/
void OsiSymSolverInterface::deleteRows(const int num, const int * rowIndices)
{
   sym_delete_rows(env_, num, const_cast<int*>(rowIndices));   
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
#ifdef USE_CGL_CUTS
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
#endif

/*===========================================================================*/
/*===========================================================================*/

CoinWarmStart * OsiSymSolverInterface::getWarmStart() const
{

#ifdef COMPILE_IN_TM

   /* FIX_ME! Add a SymIntParam to determine whether keep the tree
      in case a getWarmStart() is called or not!
   */

   SymWarmStart * symWS = 
      new SymWarmStart(sym_get_warm_start(env_, false));
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
   
   sym_set_warm_start(env_, ws);
   return true;   
}

/*===========================================================================*/
/* copy constructor, clone and assignment operator
/*===========================================================================*/


OsiSymSolverInterface::OsiSymSolverInterface(const OsiSolverInterface & si)
{

   OsiSolverInterface * si_copy = const_cast<OsiSolverInterface *>(&si);
   OsiSymSolverInterface * sym = dynamic_cast<OsiSymSolverInterface*>(si_copy);

   env_= sym_create_copy_problem(sym->getSymphonyEnvironment());
}

/*===========================================================================*/
/*===========================================================================*/

OsiSolverInterface * OsiSymSolverInterface::clone(bool copyData) const
{
   return new OsiSymSolverInterface(*this);
}

/*===========================================================================*/
/*===========================================================================*/

OsiSymSolverInterface & OsiSymSolverInterface::operator=(const OsiSolverInterface& rhs)
{
   OsiSolverInterface * si_copy = const_cast<OsiSolverInterface *>(&rhs);
   OsiSymSolverInterface * sym = dynamic_cast<OsiSymSolverInterface*>(si_copy);

   if(this != &rhs){
      env_= sym_create_copy_problem(sym->getSymphonyEnvironment());
   }

   return *this;
}

/*===========================================================================*/
/*===========================================================================*/
