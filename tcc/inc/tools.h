#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED
#include "setup.h"
#include "stdio.h"

FILE* open_file(char* file_name, char* mode);

int make_directory(const char* name);

void Error_no_free(char *);// Exit program printing error message but don't try to free any memory

void Error(char *); // Exit program printing error message and trying to free any allocated memory

int **resize_2D_int(int **, int , int , int , int );

int *resize_1D_int(int *, int , int ) ;

void links() ;  // sorts all the particles into cells, result given by head-of-chain and linked list arrays

int quickSort(int *, int );
#endif
