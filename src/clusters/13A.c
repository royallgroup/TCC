#include <globals.h>
#include "13A.h"

//!  A 13A clusters is the intersection of a 12B and a 7A cluster. Topologically it is a regular icosahedron.
/*!
*  Find 13A clusters
*  A 13A is constructed from a 12B and a 7A where:
*      - The 7A cluster has one spindle given by sc of the 12B cluster, and one spindle that is distinct from the 12B particles.
*      - The sp5 ring particles of the 7A cluster are distinct from the sp5 ring particles of the central 7A cluster in 12B.
*
*  Cluster output: SOOBBBBBBBBBB
*  Storage order: unknown
*
*/
void Clust_Write_13A() {
    s13A[hc13A[n13A][0]] = 'S';
    if(s13A[hc13A[n13A][1]] != 'S') s13A[hc13A[n13A][1]] = 'O';
    if(s13A[hc13A[n13A][2]] != 'S') s13A[hc13A[n13A][2]] = 'O';
    for(int i = 3; i < 13; i++){
        if (s13A[hc13A[n13A][i]] == 'C') s13A[hc13A[n13A][i]] = 'B';
    }

    ++n13A;
}