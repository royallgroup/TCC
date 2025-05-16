#include <globals.h>
#include <tools.h>
#include <bonds.h>
#include <clusters/simple_cluster_methods.h>
#include "10A.h"

//!  An 10A cluster is the two bonded sp4b clusters with no common particles.
/*!
*  Find 10A clusters
*  An 10A is 2 sp4b clusters where:
*      - There are no common particles.
*      - The spindles are not bonded.
*      - Each particle in each sp4 ring is bonded to exactly two particles in the other sp4 ring
*
*  Cluster output: BBBBBBBOO
*  Storage order: first_sp4b_ring x 4, second_sp4b_ring x 4, firrst_sp4b_spindle, second_sp4b_spindle
*/
void Clusters_Get10A() {
    int neighbour_pointer, k, l, m;
    int *used_sp4b;

    used_sp4b = malloc(nsp4b * sizeof(int));
    if (used_sp4b == NULL) {
        Error("Clusters_Get10A(): used_sp4b[] malloc out of memory\n");
    }

    for (int first_sp4b_id = 0; first_sp4b_id < nsp4b; ++first_sp4b_id) {
        int *first_sp4b_cluster = hcsp4b[first_sp4b_id];
        zero_used_array(used_sp4b);
        used_sp4b[first_sp4b_id] = 1;
        for (neighbour_pointer = 0; neighbour_pointer < num_bonds[first_sp4b_cluster[0]]; ++neighbour_pointer) {
            for (int neighbour_id = 0; neighbour_id < nmem_sp4b[bond_list[first_sp4b_cluster[0]][neighbour_pointer]]; ++neighbour_id) {
                int second_sp4b_id = mem_sp4b[bond_list[first_sp4b_cluster[0]][neighbour_pointer]][neighbour_id];
                if (second_sp4b_id > first_sp4b_id) {
                    int *second_sp4b_cluster = hcsp4b[second_sp4b_id];
                    if (used_sp4b[second_sp4b_id] == 0) {
                        used_sp4b[second_sp4b_id] = 1;

                        if (Bonds_BondCheck(first_sp4b_cluster[4], second_sp4b_cluster[4]) == 0) {

                            if (are_clusters_distinct(first_sp4b_cluster, second_sp4b_cluster, 5, 5) == 1) {

                                //check each ring particle is bonded to two others
                                for (k = 0; k < 4; ++k) {
                                    m = 0;
                                    for (l = 0; l < 4; ++l)
                                        if (Bonds_BondCheck(first_sp4b_cluster[k], second_sp4b_cluster[l]))++m;
                                    if (m != 2) break;
                                }
                                // ERROR: Need to check converse, i. e. each SP4 ring particle from sp4b_j bonded to to exactly two particles from sp4b_i
                                if (k == 4) {
                                    Cluster_Write_10A(first_sp4b_cluster, second_sp4b_cluster);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    free(used_sp4b);
}

void zero_used_array(int *used_sp4b) {
    for (int i=0; i < nsp4b; ++i) {
        used_sp4b[i] = 0;
    }
}

void Cluster_Write_10A(int *first_sp4b_cluster, int *second_sp4b_cluster) {
    int clusSize = 10;

    if(n10A == m10A) {
        hc10A= resize_2D_int(hc10A, m10A, m10A + incrStatic, clusSize, -1);
        m10A= m10A + incrStatic;
    }

    hc10A[n10A][0]= first_sp4b_cluster[0];
    hc10A[n10A][1]= first_sp4b_cluster[1];
    hc10A[n10A][2]= first_sp4b_cluster[2];
    hc10A[n10A][3]= first_sp4b_cluster[3];
    hc10A[n10A][4]= second_sp4b_cluster[0];
    hc10A[n10A][5]= second_sp4b_cluster[1];
    hc10A[n10A][6]= second_sp4b_cluster[2];
    hc10A[n10A][7]= second_sp4b_cluster[3];
    hc10A[n10A][8]= first_sp4b_cluster[4];
    hc10A[n10A][9]= second_sp4b_cluster[4];

    for (int i = 0; i < 8; i++) {
        if (s10A[hc10A[n10A][i]] == 'C') s10A[hc10A[n10A][i]] = 'B';
    }
    s10A[hc10A[n10A][8]] = 'O';
    s10A[hc10A[n10A][9]] = 'O';

    ++n10A;
}