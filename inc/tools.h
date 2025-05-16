#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED
#include "setup.h"
#include "stdio.h"
#include "globals.h"

long get_long_from_string(const char *buff, int *validLong);

double get_double_from_string(const char *buff, int *validDouble);

int try_read_line_from_file(FILE *file_name);

FILE* open_file(char* file_name, char* mode);

int make_directory(const char* name);

void Error_no_free(char *);// Exit program printing error message but don't try to free any memory

void Error(char *); // Exit program printing error message and trying to free any allocated memory

int **resize_2D_int(int **, int , int , int , int );

int *resize_1D_int(int *, int , int ) ;

int sort_list_of_ints(const void *lhs, const void *rhs);

int sort_list_of_lists_of_ints(const void *lhs, const void *rhs);

int quickSort(int *, int );
#endif
