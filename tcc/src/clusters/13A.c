#include <globals.h>
#include "clusters.h"
#include "13A.h"

void Clust_Write_13A() {
    int i;
    s13A[hc13A[n13A][0]] = 'S';
    if(s13A[hc13A[n13A][1]] != 'S') s13A[hc13A[n13A][1]] = 'O';
    if(s13A[hc13A[n13A][2]] != 'S') s13A[hc13A[n13A][2]] = 'O';
    for(i=3; i<13; i++){
        if (s13A[hc13A[n13A][i]] == 'C') s13A[hc13A[n13A][i]] = 'B';
    }

    ++n13A;
}