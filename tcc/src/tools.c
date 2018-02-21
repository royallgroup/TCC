#include "tools.h"
#include "globals.h"
#ifdef _WIN32
    #include "direct.h"
#endif

FILE* open_file(char* file_name, char* mode) {
    FILE *file_pointer;

    file_pointer = fopen(file_name, mode);

    return file_pointer;
}

int make_directory(const char* name) {
    int error_number;
    char errMsg[100];

    #ifdef __linux__
        if(mkdir(name, 0744) != 0) {
            error_number = errno;
        }
        else return 0;
    #else
        if(_mkdir(name) != 0) {
            error_number = errno;
        }
        else return 0;
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
    
    Setup_FreeStaticVars();
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

void links() {  // sorts all the particles into cells, result given by head-of-chain and linked list arrays
    int i, ic;
    for (ic=1;ic<=ncells;ic++) head[ic]=0;
    for (i=1;i<=N;i++) {
        ic = 1 + (int)((x[i-1]+ halfSide)*invcellSide) + M*((int)((y[i-1]+halfSide)*invcellSide)) + M*M*((int)((z[i-1]+halfSide)*invcellSide));
        if (ic > ncells || ic <= 0) {
            printf("i %d r_x %lg r_y %lg r_z %lg side %lg halfSide %lg ic %d ncells %d\n",i-1,x[i-1],y[i-1],z[i-1],side,halfSide,ic,ncells);
            Error("links(): ic > ncells, i.e. particle coord no longer in simulation box!!\n");
        }
        llist[i]=head[ic];
        head[ic]=i;
    }
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


