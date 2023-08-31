
/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/*
    cont.h: included in contml & contrast
*/


typedef struct cont_node 
{
  node node;
  phenotype3 view;
  long totalleles;

} cont_node;


#ifndef OLDC
/*function prototypes*/
node* cont_node_new(node_type, long);
void alloctree(pointarray *, long);
void freetree(pointarray *, long);
void setuptree(tree *, long);
void allocview(tree *, long, long);
void freeview(tree *, long);
void standev2(long, long, long, long, double, double *, double **, longer);
void cont_node_copy(node* src, node* dst);
void cont_node_init(cont_node* n, boolean tip, long index);
/*function prototypes*/
#endif
