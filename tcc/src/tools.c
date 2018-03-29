#include "tools.h"
#include "globals.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#if defined(__linux__) || defined(__APPLE__)
    #include <sys/stat.h>
#endif

#if defined(_WIN32)
    #include "direct.h"
#endif

long get_max_particle_number(struct xyz_info input_xyz_info) {
    int i;
    long max_particle_number = 0;

    for(i=0; i < input_xyz_info.total_frames; i++) {
        if(input_xyz_info.num_particles[i] > max_particle_number) {
            max_particle_number = input_xyz_info.num_particles[i];
        }
    }
    return max_particle_number;
}

long get_long_from_string(const char *buff, int *validLong) {
    char *end;
    errno = 0;
    long converted_number = 0;
    *validLong = 1;

    const long sl = strtol(buff, &end, 10);

    if (end == buff) {
        fprintf(stderr, "%s: not a decimal number\n", buff);
        *validLong = 0;
    } else if ((sl == LONG_MIN || sl == LONG_MAX) && errno == ERANGE) {
        fprintf(stderr, "%s out of range of type long\n", buff);
        *validLong = 0;
    } else {
        converted_number = sl;
    }
    return converted_number;
}

double get_double_from_string(const char *buff, int *validDouble) {
    char *end;
    errno = 0;
    double converted_number = 0;
    *validDouble = 1;

    const double sl = strtod(buff, &end);

    if (end == buff) {
        fprintf(stderr, "%s: not a valid double number\n", buff);
        *validDouble = 0;
    } else if ((sl == HUGE_VAL || sl == -HUGE_VAL) && errno == ERANGE) {
        fprintf(stderr, "%s out of range of type float\n", buff);
        *validDouble = 0;
    } else {
        converted_number = sl;
    }
    return converted_number;
}

int try_read_line_from_file(FILE *file_name) {
    // Try to read a line from a file. If line successfuly read return 1 else return 0.
    char line[1000];
    size_t line_length;

    if (fgets(line, 1000, file_name) == NULL) {
        printf("EOF reached.");
    }
    else {
        line_length = strlen(line);
        if(line_length != 0 && line[line_length-1] == '\n') {
            return 1;
        }
        else {
            return 0;
        }
    }
}

FILE* open_file(char* file_name, char* mode) {
    FILE *file_pointer;

    file_pointer = fopen(file_name, mode);

    return file_pointer;
}

int make_directory(const char* name) {
    int error_number;
    char errMsg[100];

    #if defined(_WIN32)
        if(_mkdir(name) != 0) {
            error_number = errno;
        }
        else return 0;
    #elif defined(__linux__) || defined(__APPLE__)
        if(mkdir(name, 0744) != 0) {
            error_number = errno;
        }
        else return 0;
    #else
        Error("Compiler Not Recognised");
    #endif

    if(error_number == 17) {
        return 0;
    }
    else {
        sprintf(errMsg,"Error creating directory %s.\n",name);
        Error(errMsg);
    }
}

void Error_no_free(char *msg) { // Exit program printing error message but don't try to free any memory
    printf("\n%s\n", msg);
    exit(1);
}

void Error(char *msg) { // Exit program printing error message and trying to free any allocated memory
    printf("\n%s\n",msg);

    Free_All_Variables();
    exit(1);
}

int **resize_2D_int(int **the_array, int old_row_size, int new_row_size, int new_col_size, int value) {
    int i, j;
    char errMsg[1000];
    
    the_array=realloc(the_array,new_row_size*sizeof(int *));
    if (the_array == NULL) { sprintf(errMsg,"resize_2D_int(): the_array[] out of memory old_row_size %d new_row_size %d new_col_size %d\n",old_row_size,new_row_size,new_col_size); Error_no_free(errMsg); }
    for (i=old_row_size; i<new_row_size; i++) {
        the_array[i] = malloc(new_col_size*sizeof(int));
        if (the_array[i] == NULL) { sprintf(errMsg,"resize_2D_int(): the_array[][] out of memory\n"); Error_no_free(errMsg); }
    }
    for (i=old_row_size; i<new_row_size; i++) {
        for (j=0; j<new_col_size; j++) {
            the_array[i][j]=value;
        }
    }
    
    return the_array;
}

int *resize_1D_int(int *the_array, int old_col_size, int new_col_size) {
    int i;
    char errMsg[1000];
    
    the_array=realloc(the_array,new_col_size*sizeof(int));
    if (the_array == NULL) { sprintf(errMsg,"resize_1D_int(): the_array[] out of memory old_col_size %d new_col_size %d\n",old_col_size,new_col_size); Error_no_free(errMsg); }
    for (i=old_col_size; i<new_col_size; i++) the_array[i]=-1;
    
    return the_array;
}

int sort_list_of_lists_of_ints(const void *lhs, const void *rhs) {
    int* n1 = *(int **) lhs;
    int* n2 = *(int **) rhs;
    int i;

    for(i=0; i < num_sort_columns; i++) {
        if (n1[i] < n2[i]) return -1;
        else if (n2[i] < n1[i]) return 1;
    }
    return 0;
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BEGIN QUICKSORT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * Returns YES if sort was successful, or NO if the nested
//    pivots went too deep, in which case your array will have
//    been re-ordered, but probably not sorted correctly.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7

int quickSort(int *arr, int elements) {
    int piv, beg[10], end[10], i, L, R;

    i=0;
    beg[0]=0; 
    end[0]=elements;
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R) {
            piv=arr[L]; 
            if (i==10-1) return 0;
            while (L<R) {
                while (arr[R]>=piv && L<R) R--; 
                if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++; 
                if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv; 
            beg[i+1]=L+1; 
            end[i+1]=end[i]; 
            end[i++]=L; 
        }
        else i--;
    }
    return 1; 
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END QUICKSORT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


