#ifndef BACH_H_
#define BACH_H_

#include "sym_lp.h"

#define BACH_BRANCH 531
#define BACH_CONTINUE 387
#define BACH_ALPHA 0.001
#define BACH_BETA 0.1
#define BACH_GAMMA 0.1


 int bach_init_node(int nodeID, int verbosity, int method, bach_node** bachnode );
 bach_node* bach_record_value(int ID, bach_node* bachnode, double objval, double conditionnumber);
 double bach_condition_number(LPdata* lpdata);
 int bach_should_we_branch(bach_node* bachnode);
 int bach_print_all_values(bach_node* bachnode);
 int bach_free_node(bach_node* bachnode );
 //int bach_print(char* message, int verbosity);
 
#endif

