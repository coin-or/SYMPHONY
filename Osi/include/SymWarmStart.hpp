
#ifndef CoinSymWarmStartBasis_H
#define CoinSymWarmStartBasis_H

#include "CoinHelperFunctions.hpp"
#include "CoinWarmStart.hpp"
#include "BB_types.h"
#include "BB_constants.h"
#include "tm.h"
#include <iostream>
#include <fstream>

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
   SymWarmStart(warm_start_desc *& ws, bool dominateWarmStart);
   
   /*Get the warmStart info from a file*/
   SymWarmStart(char *f);
   
   /* Copy constructor */
   SymWarmStart(const SymWarmStart & symWS);

   /* Destructor */
   virtual ~SymWarmStart();

   /* Clone the warmstart */
   virtual CoinWarmStart * clone() const; 

   /* Get the pointer to the loaded warmStart_ */
   virtual const warm_start_desc * getWarmStartDesc() const;

   /* Trim the sub tree rooted at node */
   virtual void trimTree(bc_node *& node);

   /* Move the pointer to the rootnode of the warmStart to another
      node which will change the underlying tree 
   */
   virtual void setRoot(bc_node *root) {} //FIX_ME! Ask Prof. Ralphs.

   /* Write the current warm start info to a file */
   virtual bool writeToFile(char * f);

private:

   /* Private warm start desc. to keep everything */
   warm_start_desc *warmStart_;

    /* Write the underlying tree to the given file */
   void writeTreeToFile_(bc_node * root, ofstream & os);

   /* Write the node to the given file */
   void writeNode_(bc_node *node, ofstream & os);
  
   /* Get the warmStart info from the given file */
   void readWarmStartFromFile_(ifstream & is);

   /* Get the tree desc. rooted at 'root' from the given file */
   void readTreeFromFile_(bc_node * root, ifstream & is);

   /* Read in the next description in the given file to the 'node'*/
   void readNode_(bc_node *node, ifstream & is);

   /* Copy the tree rooted at 'rootFrom' to memory without interrupting
      the given tree
   */
   void copyTree_(bc_node *rootTo, bc_node * rootFrom);

   /* Copy the node 'nFrom' to 'nTo' */
   void copyNode_(bc_node * nTo, bc_node * nFrom);

   /* Delete the given node */
   void deleteNode_(bc_node *n);
};

#endif
