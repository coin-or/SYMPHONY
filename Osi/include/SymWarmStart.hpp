
#ifndef SymWarmStart_H
#define SymWarmStart_H

#include "CoinHelperFunctions.hpp"
#include "CoinWarmStart.hpp"
#include <iostream>
#include "symphony_api.h"

using namespace std;

//#############################################################################

class SymWarmStart : public CoinWarmStart 
{

public:

   /* Default constructor. Will do nothing! */
   SymWarmStart(){}
   
   /* Initialize the warmStart_ using the given warm start. If dominate
      WarmStart is set, then, SymWarmStart will take the control of the 
      given description, otherwise, will copy everything.
   */
   SymWarmStart(warm_start_desc * ws);
   
   /*Get the warmStart info from a file*/
   SymWarmStart(char *f);
   
   /* Copy constructor */
   SymWarmStart(const SymWarmStart & symWS);

   /* Destructor */
   virtual ~SymWarmStart();

   /* Clone the warmstart */
   virtual CoinWarmStart * clone() const; 

   /* Get the pointer to the loaded warmStart_ */
   virtual const warm_start_desc * getWarmStartDesc();

   /* Trim the sub tree rooted at node */
   virtual void trimTree(bc_node * node);

   /* Move the pointer to the rootnode of the warmStart to another
      node which will change the underlying tree 
   */
   virtual void setRoot(bc_node *root) {} //FIX_ME! Ask Prof. Ralphs.

   /* Write the current warm start info to a file */
   virtual bool writeToFile(char * f);

private:

   /* Private warm start desc. to keep everything */
   warm_start_desc *warmStart_;

};

#endif
