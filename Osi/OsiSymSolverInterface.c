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

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymIntParam(OsiSymIntParam key, int value)
{
   switch(key){

    case OsiSymVerbosity:
      env_->par.verbosity = value;
      return true;

    case OsiSymWarmStart:
      env_->par.tm_par.warm_start = value;
      return true;

    case OsiSymNodeLimit:
      env_->par.tm_par.node_limit = value;
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
      
    case OsiObjOffset:
      env_->mip->obj_offset = value;
      return true;

    default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymDblParam(OsiSymDblParam key, double value)
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

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymStrParam(OsiSymStrParam key, 
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

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymIntParam(OsiSymIntParam key, int& value)
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

   default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
   switch (key){
      
    case OsiObjOffset:
      value = env_->mip->obj_offset;
      return true;

   default:
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymDblParam(OsiSymDblParam key, 
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

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymStrParam(OsiSymStrParam key, 
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
