
#include "SymWarmStart.hpp"

//#############################################################################

SymWarmStart::SymWarmStart(warm_start_desc * ws)
{
   warmStart_ = sym_create_copy_warm_start(ws);
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::SymWarmStart(char * fileName)
{
   
   warmStart_ = sym_read_warm_start(fileName);
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::SymWarmStart(const SymWarmStart & symWS)
{
   warm_start_desc * wsCopy;
   SymWarmStart * sWS = const_cast<SymWarmStart *>(&symWS);
   wsCopy = const_cast<warm_start_desc *>(sWS->getCopyOfWarmStartDesc());

   warmStart_ = wsCopy;
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::~SymWarmStart()
{
   sym_delete_warm_start(warmStart_);
}

/*===========================================================================*/
/*===========================================================================*/

CoinWarmStart * SymWarmStart::clone () const
{
   return new SymWarmStart(*this);
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc * SymWarmStart::getCopyOfWarmStartDesc()
{

   if(warmStart_){
      return(sym_create_copy_warm_start(warmStart_));
   }
   else{
      cout<<"getWarmStart(): No loaded warm start desc. to return!"<<endl;
      return 0;
   }
}

/*===========================================================================*/
/*===========================================================================*/

int SymWarmStart::writeToFile(char * fileName)
{
   return(sym_write_warm_start_desc(warmStart_, fileName));
}

/*===========================================================================*/
/*===========================================================================*/
