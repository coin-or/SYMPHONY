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

#ifndef OsiSymSolverInterface_hpp
#define OsiSymSolverInterface_hpp

#include "OsiSolverInterface.hpp"
#include "OsiSymSolverParameters.hpp"
#include "symphony_api.h"

//#############################################################################

/** OSI Solver Interface for SYMPHONY

  Many OsiSolverInterface query methods return a const pointer to the
  requested read-only data. If the model data is changed or the solver
  is called, these pointers may no longer be valid and should be 
  refreshed by invoking the member function to obtain an updated copy
  of the pointer.
  For example:
  \code
      OsiSolverInterface solverInterfacePtr ;
      const double * ruBnds = solverInterfacePtr->getRowUpper();
      solverInterfacePtr->applyCuts(someSetOfCuts);
      // ruBnds is no longer a valid pointer and must be refreshed
      ruBnds = solverInterfacePtr->getRowUpper();
  \endcode

  Querying a problem that has no data associated with it will result in
  zeros for the number of rows and columns, and NULL pointers from
  the methods that return vectors.
*/

class OsiSymSolverInterface : public OsiSolverInterface {

public:
  ///@name Solve methods 
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve(){
       throw CoinError("Error: SYMPHONY is only meant to solve MIPs, not LPs",
		       "initialSolve", "OsiSymSolverInterface");
    }; 

    /// Resolve an LP relaxation after problem modification
    virtual void resolve(){
       throw CoinError("Error: SYMPHONY is only meant to solve MIPs, not LPs",
		       "resolve", "OsiSymSolverInterface");
    };

    /// Invoke solver's built-in enumeration algorithm
    virtual void branchAndBound();
  //@}

  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. There can be various reasons for failure: the given
     parameter is not applicable for the solver (e.g., refactorization
     frequency for the volume algorithm), the parameter is not yet implemented
     for the solver or simply the value of the parameter is out of the range
     the solver accepts. If a parameter setting call returns false check the
     details of your solver.

     The get methods return true if the given parameter is applicable for the
     solver and is implemented. In this case the value of the parameter is
     returned in the second argument. Otherwise they return false.

     \note
     There is a default implementation of the set/get
     methods, namely to store/retrieve the given value using an array in the
     base class. A specific solver implementation can use this feature, for
     example, to store parameters that should be used later on. Implementors
     of a solver interface should overload these functions to provide the
     proper interface to and accurately reflect the capabilities of a
     specific solver.

     The format for hints is slightly different in that the value is 
     boolean and there is an enum to show strength of hint.
     There is also an optional void pointer to allow for any eventuality.
     Hints should be initialised when a solver is instantiated.
     (See OsiSolverParameters.hpp for defined hint parameters and strength.)
     A value of true means to work with the hint, false to work against it.
     For example,
     <ul>
       <li> \code setHintParam(OsiDoScale,true,OsiHintTry) \endcode
	    is a mild suggestion to the solver to scale the constraint
	    system.
       <li> \code setHintParam(OsiDoScale,false,OsiForceDo) \endcode
	    tells the solver to disable scaling, or throw an exception if
	    it cannot comply.
     </ul>
     As another example, a solver interface could use the value and strength
     of the \c OsiDoReducePrint hint to adjust the amount of information
     printed by the interface and/or solver.  The extent to which a solver
     obeys hints is left to the solver.  The value and strength returned by
     \c getHintParam will match the most recent call to \c setHintParam,
     and will not necessarily reflect the solver's ability to comply with the
     hint.  If the hint strength is \c OsiForceDo, the solver is required to
     throw an exception if it cannot perform the specified action.

     \note
     As with the other set/get methods, there is a default implementation
     which maintains arrays in the base class for hint value and strength.
     The default implementation does not store the void pointer, and always
     throws an exception for strength \c OsiForceDo. Implementors of a solver
     interface should overload these functions to provide the proper interface
     to and accurately reflect the capabilities of a specific solver.
  */
  //@{
    // Set an integer parameter
    virtual bool setIntParam(OsiIntParam key, int value);

    // Set SYMPHONY int parameter
    bool setSymParam(OsiSymIntParam key, int value);

    // Set an double parameter
    virtual bool setDblParam(OsiDblParam key, double value);

    // Set SYMPHONY double parameter
    bool setSymParam(OsiSymDblParam key, double value);

    // Set an string parameter
    virtual bool setStrParam(OsiStrParam key, const std::string & value);

    // Set SYMPHONY string parameter
    bool setSymParam(OsiSymStrParam key, const std::string & value);

    // Get an integer parameter
    virtual bool getIntParam(OsiIntParam key, int& value) const;

    // Get SYMPHONY int parameter
    bool getSymParam(OsiSymIntParam key, int& value) const;

    // Get an double parameter
    virtual bool getDblParam(OsiDblParam key, double& value) const;

    // Get SYMPHONY double parameter
    bool getSymParam(OsiSymDblParam key, double& value) const;

    // Get a string parameter
    virtual bool getStrParam(OsiStrParam key, std::string& value) const;

    // Get SYMPHONY string parameter
    bool getSymParam(OsiSymStrParam key, std::string & value) const;

  //@}

  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
    /// Are there numerical difficulties?
    virtual bool isAbandoned() const{
       throw CoinError("Error: Function not implemented",
		       "isAbandoned", "OsiSymSolverInterface");
    };
    /// Is optimality proven?
    virtual bool isProvenOptimal() const{
       throw CoinError("Error: Function not implemented",
		       "isProvenOptimal", "OsiSymSolverInterface");
    };
    /// Is primal infeasiblity proven?
    virtual bool isProvenPrimalInfeasible() const{
       throw CoinError("Error: Function not implemented",
		       "isProvenPrimalInfeasible", "OsiSymSolverInterface");
    };
    /// Is dual infeasiblity proven?
    virtual bool isProvenDualInfeasible() const{
       throw CoinError("Error: Function not implemented",
		       "", "OsiSymSolverInterface");
    };
    /// Is the given primal objective limit reached?
    virtual bool isPrimalObjectiveLimitReached() const{
       throw CoinError("Error: Function not implemented",
		       "isProvenDualInfeasible", "OsiSymSolverInterface");
    };
    /// Is the given dual objective limit reached?
    virtual bool isDualObjectiveLimitReached() const{
       throw CoinError("Error: Function not implemented",
		       "isDualObjectiveLimitReached", "OsiSymSolverInterface");
    };
    /// Iteration limit reached?
    virtual bool isIterationLimitReached() const{
       throw CoinError("Error: Function not implemented",
		       "isIterationLimitReached", "OsiSymSolverInterface");
    };
  //@}

  //---------------------------------------------------------------------------
  /**@name Warm start methods */
  //@{
    /*! \brief Get an empty warm start object
      
      This routine returns an empty warm start object. Its purpose is
      to provide a way to give a client a warm start object of the
      appropriate type, which can resized and modified as desired.
    */

    virtual CoinWarmStart *getEmptyWarmStart () const{
       throw CoinError("Error: Function not implemented",
		       "getEmptyWarmStart", "OsiSymSolverInterface");
    };

    /** Get warm start information.

      If there is no valid solution, an empty warm start object (0 rows, 0
      columns) wil be returned.
    */
    virtual CoinWarmStart* getWarmStart() const{
       throw CoinError("Error: Function not implemented",
		       "getWarmStart", "OsiSymSolverInterface");
    };

    /** Set warm start information.
    
      Return true/false depending on whether the warm start information was
      accepted or not. */
    virtual bool setWarmStart(const CoinWarmStart* warmstart){
       throw CoinError("Error: Function not implemented",
		       "setWarmStart", "OsiSymSolverInterface");
    };
  //@}

  //---------------------------------------------------------------------------
    /**@name Problem query methods

     Querying a problem that has no data associated with it will result in
     zeros for the number of rows and columns, and NULL pointers from
     the methods that return vectors.
     
     Const pointers returned from any data-query method are valid as
     long as the data is unchanged and the solver is not called.
    */
    //@{
      /// Get pointer to SYMPHONY environment (eventually we won't need this)
      problem* getSymphonyEnvironment() const {return env_;}

      /// Get number of columns
      virtual int getNumCols() const{
       throw CoinError("Error: Function not implemented",
		       "getNumCols", "OsiSymSolverInterface");
    };
  
      /// Get number of rows
      virtual int getNumRows() const{
       throw CoinError("Error: Function not implemented",
		       "getNumRows", "OsiSymSolverInterface");
    };
  
      /// Get number of nonzero elements
      virtual int getNumElements() const{
       throw CoinError("Error: Function not implemented",
		       "getNumElements", "OsiSymSolverInterface");
    };
  
      /// Get pointer to array[getNumCols()] of column lower bounds
      virtual const double * getColLower() const{
       throw CoinError("Error: Function not implemented",
		       "getColLower", "OsiSymSolverInterface");
    };
  
      /// Get pointer to array[getNumCols()] of column upper bounds
      virtual const double * getColUpper() const{
       throw CoinError("Error: Function not implemented",
		       "getColUpper", "OsiSymSolverInterface");
    };
  
      /** Get pointer to array[getNumRows()] of row constraint senses.
  	<ul>
  	<li>'L': <= constraint
  	<li>'E': =  constraint
  	<li>'G': >= constraint
  	<li>'R': ranged constraint
  	<li>'N': free constraint
  	</ul>
      */
      virtual const char * getRowSense() const{
       throw CoinError("Error: Function not implemented",
		       "getRowSense", "OsiSymSolverInterface");
    };
  
      /** Get pointer to array[getNumRows()] of row right-hand sides
  	<ul>
  	  <li> if getRowSense()[i] == 'L' then
	       getRightHandSide()[i] == getRowUpper()[i]
  	  <li> if getRowSense()[i] == 'G' then
	       getRightHandSide()[i] == getRowLower()[i]
  	  <li> if getRowSense()[i] == 'R' then
	       getRightHandSide()[i] == getRowUpper()[i]
  	  <li> if getRowSense()[i] == 'N' then
	       getRightHandSide()[i] == 0.0
  	</ul>
      */
      virtual const double * getRightHandSide() const{
       throw CoinError("Error: Function not implemented",
		       "getRightHandSide", "OsiSymSolverInterface");
    };
  
      /** Get pointer to array[getNumRows()] of row ranges.
  	<ul>
            <li> if getRowSense()[i] == 'R' then
                    getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
            <li> if getRowSense()[i] != 'R' then
                    getRowRange()[i] is 0.0
          </ul>
      */
      virtual const double * getRowRange() const{
       throw CoinError("Error: Function not implemented",
		       "getRowRange", "OsiSymSolverInterface");
    };
  
      /// Get pointer to array[getNumRows()] of row lower bounds
      virtual const double * getRowLower() const{
       throw CoinError("Error: Function not implemented",
		       "getRowLower", "OsiSymSolverInterface");
    };
  
      /// Get pointer to array[getNumRows()] of row upper bounds
      virtual const double * getRowUpper() const{
       throw CoinError("Error: Function not implemented",
		       "getRowUpper", "OsiSymSolverInterface");
    };
  
      /// Get pointer to array[getNumCols()] of objective function coefficients
      virtual const double * getObjCoefficients() const{
       throw CoinError("Error: Function not implemented",
		       "getObjCoefficients", "OsiSymSolverInterface");
    };
  
      /// Get objective function sense (1 for min (default), -1 for max)
   virtual double getObjSense() const{
       throw CoinError("Error: Function not implemented",
		       "getObjSense", "OsiSymSolverInterface");
    };
  
      /// Return true if variable is continuous
      virtual bool isContinuous(int colIndex) const{
       throw CoinError("Error: Function not implemented",
		       "isContinuous", "OsiSymSolverInterface");
    };
  
      /// Return true if variable is binary
      virtual bool isBinary(int colIndex) const{
       throw CoinError("Error: Function not implemented",
		       "isBinary", "OsiSymSolverInterface");
    };
  
      /** Return true if column is integer.
          Note: This function returns true if the the column
          is binary or a general integer.
      */
      virtual bool isInteger(int colIndex) const{
       throw CoinError("Error: Function not implemented",
		       "isInteger", "OsiSymSolverInterface");
    };
  
      /// Return true if variable is general integer
      virtual bool isIntegerNonBinary(int colIndex) const{
       throw CoinError("Error: Function not implemented",
		       "isIntegerNonBinary", "OsiSymSolverInterface");
    };
  
      /// Return true if variable is binary and not fixed at either bound
      virtual bool isFreeBinary(int colIndex) const{
       throw CoinError("Error: Function not implemented",
		       "isFreeBinary", "OsiSymSolverInterface");
    }; 
    
      /// Get pointer to row-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByRow() const{
       throw CoinError("Error: Function not implemented",
		       "getMatrixByRow", "OsiSymSolverInterface");
    };
  
      /// Get pointer to column-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByCol() const{
       throw CoinError("Error: Function not implemented",
		       "getMatrixByCol", "OsiSymSolverInterface");
    };
  
      /// Get solver's value for infinity
      virtual double getInfinity() const{
       throw CoinError("Error: Function not implemented",
		       "getInfinity", "OsiSymSolverInterface");
    };
    //@}
    
    /**@name Solution query methods */
    //@{
      /// Get pointer to array[getNumCols()] of primal variable values
      virtual const double * getColSolution() const{
       throw CoinError("Error: Function not implemented",
		       "getColSolution", "OsiSymSolverInterface");
    };
  
      /// Get pointer to array[getNumRows()] of dual variable values
      virtual const double * getRowPrice() const{
       throw CoinError("Error: Function not implemented",
		       "getRowPrice", "OsiSymSolverInterface");
    };
  
      /// Get a pointer to array[getNumCols()] of reduced costs
      virtual const double * getReducedCost() const{
       throw CoinError("Error: Function not implemented",
		       "getReducedCost", "OsiSymSolverInterface");
    };
  
      /** Get pointer to array[getNumRows()] of row activity levels (constraint
  	matrix times the solution vector). */
      virtual const double * getRowActivity() const{
       throw CoinError("Error: Function not implemented",
		       "getRowActivity", "OsiSymSolverInterface");
    };
  
      /// Get objective function value
      virtual double getObjValue() const{
       throw CoinError("Error: Function not implemented",
		       "getObjValue", "OsiSymSolverInterface");
    };

      /** Get the number of iterations it took to solve the problem (whatever
	  ``iteration'' means to the solver). */
      virtual int getIterationCount() const{
       throw CoinError("Error: Function not implemented",
		       "getIterationCount", "OsiSymSolverInterface");
    };
  
      /** Get as many dual rays as the solver can provide. In case of proven
          primal infeasibility there should be at least one.
     
          \note
	  Implementors of solver interfaces note that
          the double pointers in the vector should point to arrays of length
          getNumRows() and they should be allocated via new[].
     
          \note
	  Clients of solver interfaces note that
          it is the client's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getDualRays(int maxNumRays) const{
       throw CoinError("Error: Function not implemented",
		       "getDualRays", "OsiSymSolverInterface");
    };
      /** Get as many primal rays as the solver can provide. (In case of proven
          dual infeasibility there should be at least one.)
     
          <strong>NOTE for implementers of solver interfaces:</strong> <br>
          The double pointers in the vector should point to arrays of length
          getNumCols() and they should be allocated via new[]. <br>
     
          <strong>NOTE for users of solver interfaces:</strong> <br>
          It is the user's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getPrimalRays(int maxNumRays) const{
       throw CoinError("Error: Function not implemented",
		       "getPrimalRays", "OsiSymSolverInterface");
    };
  
    //@}

    //-------------------------------------------------------------------------
    /**@name Methods to modify the objective, bounds, and solution

       For functions which take a set of indices as parameters
       (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
       \c setRowSetTypes()), the parameters follow the C++ STL iterator
       convention: \c indexFirst points to the first index in the
       set, and \c indexLast points to a position one past the last index
       in the set.
    
    */
    //@{
      /** Set an objective function coefficient */
      virtual void setObjCoeff( int elementIndex, double elementValue ){
       throw CoinError("Error: Function not implemented",
		       "setObjCoeff", "OsiSymSolverInterface");
    };

      /** Set a single column lower bound.
    	  Use -getInfinity() for -infinity. */
      virtual void setColLower( int elementIndex, double elementValue ){
       throw CoinError("Error: Function not implemented",
		       "setColLower", "OsiSymSolverInterface");
    };
      
      /** Set a single column upper bound.
    	  Use getInfinity() for infinity. */
      virtual void setColUpper( int elementIndex, double elementValue ){
       throw CoinError("Error: Function not implemented",
		       "setColUpper", "OsiSymSolverInterface");
    };
      
       /** Set a single row lower bound.
    	  Use -getInfinity() for -infinity. */
      virtual void setRowLower( int elementIndex, double elementValue ){
       throw CoinError("Error: Function not implemented",
		       "setRowLower", "OsiSymSolverInterface");
    };
      
      /** Set a single row upper bound.
    	  Use getInfinity() for infinity. */
      virtual void setRowUpper( int elementIndex, double elementValue ){
       throw CoinError("Error: Function not implemented",
		       "setRowUpper", "OsiSymSolverInterface");
    };
    
      /** Set the type of a single row */
      virtual void setRowType(int index, char sense, double rightHandSide,
    			      double range){
       throw CoinError("Error: Function not implemented",
		       "setRowType", "OsiSymSolverInterface");
    };
    
    /// Set the objective function sense.
    /// (1 for min (default), -1 for max)
    virtual void setObjSense(double s){
       throw CoinError("Error: Function not implemented",
		       "setObjSense", "OsiSymSolverInterface");
    };
    
    /** Set the primal solution variable values
    
	colsol[getNumCols()] is an array of values for the primal variables.
	These values are copied to memory owned by the solver interface object
	or the solver.  They will be returned as the result of getColSolution()
	until changed by another call to setColSolution() or by a call to any
	solver routine.  Whether the solver makes use of the solution in any
	way is solver-dependent.
    */
    virtual void setColSolution(const double *colsol){
       throw CoinError("Error: Function not implemented",
		       "setColSolution", "OsiSymSolverInterface");
    };

    /** Set dual solution variable values

	rowprice[getNumRows()] is an array of values for the dual
	variables. These values are copied to memory owned by the solver
	interface object or the solver.  They will be returned as the result of
	getRowPrice() until changed by another call to setRowPrice() or by a
	call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */

   virtual void setRowPrice(const double * rowprice){
       throw CoinError("Error: Function not implemented",
		       "setRowPrice", "OsiSymSolverInterface");
    };

    //@}

    //-------------------------------------------------------------------------
    /**@name Methods to set variable type */
    //@{
      /** Set the index-th variable to be a continuous variable */
      virtual void setContinuous(int index){
       throw CoinError("Error: Function not implemented",
		       "setContinuous", "OsiSymSolverInterface");
    };
      /** Set the index-th variable to be an integer variable */
      virtual void setInteger(int index){
       throw CoinError("Error: Function not implemented",
		       "setInteger", "OsiSymSolverInterface");
    };
    //@}
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    /**@name Methods to expand a problem.

       Note that new columns are added as continuous variables.

    */
    //@{
      /** Add a column (primal variable) to the problem. */
      virtual void addCol(const CoinPackedVectorBase& vec,
			  const double collb, const double colub,   
			  const double obj){
       throw CoinError("Error: Function not implemented",
		       "addCol", "OsiSymSolverInterface");
    };

      /** Remove a set of columns (primal variables) from the problem.  */
      virtual void deleteCols(const int num, const int * colIndices){
       throw CoinError("Error: Function not implemented",
		       "deleteCols", "OsiSymSolverInterface");
    };
    
      /** Add a row (constraint) to the problem. */
      virtual void addRow(const CoinPackedVectorBase& vec,
    			  const double rowlb, const double rowub){
       throw CoinError("Error: Function not implemented",
		       "addRow", "OsiSymSolverInterface");
    };
      /** */
      virtual void addRow(const CoinPackedVectorBase& vec,
    			  const char rowsen, const double rowrhs,   
    			  const double rowrng){
       throw CoinError("Error: Function not implemented",
		       "addRow", "OsiSymSolverInterface");
    };

      /** Delete a set of rows (constraints) from the problem. */
      virtual void deleteRows(const int num, const int * rowIndices){
       throw CoinError("Error: Function not implemented",
		       "deleteRows", "OsiSymSolverInterface");
    };
    
    //@}

  //---------------------------------------------------------------------------

  /**@name Methods to input a problem */
  //@{

    virtual void loadProblem();
   
    /** Load in an problem by copying the arguments (the constraints on the
        rows are given by lower and upper bounds). If a pointer is 0 then the
        following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>
    */
    virtual void loadProblem(const CoinPackedMatrix& matrix,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub){
       throw CoinError("Error: Function not implemented",
		       "loadProblem", "OsiSymSolverInterface");
    };
			    
    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by lower and upper bounds).
	For default values see the previous method.

	\warning
	The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       double*& rowlb, double*& rowub){
       throw CoinError("Error: Function not implemented",
		       "assignProblem", "OsiSymSolverInterface");
    };

    /** Load in an problem by copying the arguments (the constraints on the
	rows are given by sense/rhs/range triplets). If a pointer is 0 then the
	following values are the default:
	<ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
	  <li> <code>obj</code>: all variables have 0 objective coefficient
          <li> <code>rowsen</code>: all rows are >=
          <li> <code>rowrhs</code>: all right hand sides are 0
          <li> <code>rowrng</code>: 0 for the ranged rows
        </ul>
    */
    virtual void loadProblem(const CoinPackedMatrix& matrix,
			     const double* collb, const double* colub,
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng){
       throw CoinError("Error: Function not implemented",
		       "loadProblem", "OsiSymSolverInterface");
    };

    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by sense/rhs/range triplets). For
        default values see the previous method.

	\warning
	The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       char*& rowsen, double*& rowrhs,
			       double*& rowrng){
       throw CoinError("Error: Function not implemented",
		       "assignProblem", "OsiSymSolverInterface");
    };

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const CoinBigIndex * start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub){
       throw CoinError("Error: Function not implemented",
		       "loadProblem", "OsiSymSolverInterface");
    };

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const CoinBigIndex * start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng){
       throw CoinError("Error: Function not implemented",
		       "loadProblem", "OsiSymSolverInterface");
    };

    /** Write the problem in MPS format to the specified file.

      If objSense is non-zero, a value of -1.0 causes the problem to be
      written with a maximization objective; +1.0 forces a minimization
      objective. If objSense is zero, the choice is left to implementation.
    */
    virtual void writeMps(const char *filename,
			  const char *extension = "mps",
			  double objSense=0.0) const{
       throw CoinError("Error: Function not implemented",
		       "writeMps", "OsiSymSolverInterface");
    };

   void parseCommandLine(int argc, char **argv);

   void findInitialBounds();

   int createPermanentCutPools();

  //@}

  //---------------------------------------------------------------------------

  ///@name Constructors and destructors
  //@{
    /// Default Constructor
    OsiSymSolverInterface(); 
    
    /** Clone

      The result of calling clone(false) is defined to be equivalent to
      calling the default constructor OsiSolverInterface().
    */
    virtual OsiSolverInterface * clone(bool copyData = true) const{
       throw CoinError("Error: Function not implemented",
		       "clone", "OsiSymSolverInterface");
    };
  
    /// Copy constructor 
    OsiSymSolverInterface(const OsiSolverInterface &);
  
    /// Assignment operator 
    OsiSymSolverInterface & operator=(const OsiSolverInterface& rhs);
  
    /// Destructor 
    virtual ~OsiSymSolverInterface ();

    /** Reset the solver interface.

    A call to reset() returns the solver interface to the same state as
    it would have if it had just been constructed by calling the default
    constructor OsiSolverInterface().
    */
    virtual void reset();
  //@}

  //---------------------------------------------------------------------------

protected:
  ///@name Protected methods
  //@{
    /** Apply a row cut (append to the constraint matrix). */
    virtual void applyRowCut( const OsiRowCut & rc ){
       throw CoinError("Error: Function not implemented",
		       "applyRowCut", "OsiSymSolverInterface");
    };

    /** Apply a column cut (adjust the bounds of one or more variables). */
    virtual void applyColCut( const OsiColCut & cc ){
       throw CoinError("Error: Function not implemented",
		       "applyColCut", "OsiSymSolverInterface");
    };

    /** Set OsiSolverInterface object state for default constructor

      This routine establishes the initial values of data fields in the
      OsiSolverInterface object when the object is created using the
      default constructor.
    */
    void setInitialData();
  //@}

private:

  /**@name Private member data */
  //@{
   /// The pointer to the SYMPHONY problem environment
   problem* env_;
  //@}
};

#endif
