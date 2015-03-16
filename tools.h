#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED
#include "setup.h"

void Error_no_free(char *);// Exit program printing error message but don't try to free any memory

void Error(char *); // Exit program printing error message and trying to free any allocated memory
// resize 2D arrays
int **resize_2D_int(int **, int , int , int , int );

int *resize_1D_int(int *, int , int ) ;

double *resize_1D_double(double *, int , int ) ;

void links() ;  // sorts all the particles into cells, result given by head-of-chain and linked list arrays


void Dyn_add(int *, int f, int , int *, int *, int* **, int* **, int , int *, int n, int , int , int* **, int *);
void Dyn_add_6A(int , int *, int f, int , int *, int *, int* **, int* **, int , int *, int n, int , int , int* **, int *) ;
void Dyn_add_8A(int *, int f, int , int *, int *, int* **, int* **, int , int *, int n, int , int , int* **, int *);
int quickSort(int *, int );
#endif
