#ifndef BACH_TYPES_H_
#define BACH_TYPES_H_

// Structures
typedef struct BACH_NODE{
    double* objvalues;
    double* conditions;
    int id;
    int counter;
    int size;
    int verbosity;
    int method;
 }bach_node;
 
#endif

