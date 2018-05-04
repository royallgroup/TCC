#include <globals.h>
#include <tools.h>
#include <bonds.h>
#include "10A.h"

void Clusters_Get10A() {

    //!  An 10A cluster is the two bonded sp4b clusters with no common particles.
    /*!
   *  Find 10A clusters
   *  An 10A is 2 sp4b clusters where:
   *      - There are no common particles.
   *      - The spindles are not bonded.
   *      - Each particle in each sp4 ring is bonded to exactly two particles in the other sp4 ring
   *
   *  Cluster output: BBBBBBBOO
   *  Storage order: first_cluster_ring x 4, second_cluster_ring x 4, first_cluster_spindle, second_cluster_spindle
   */

    int k, l, m;

    for (int first_sp4b_id = 0; first_sp4b_id < nsp4b; ++first_sp4b_id) {
        int *first_sp4b_cluster = hcsp4b[first_sp4b_id];
        for (int neighbour_pointer = 0; neighbour_pointer < num_bonds[first_sp4b_cluster[0]]; ++neighbour_pointer) {
            for (int second_sp4b_pointer = 0; second_sp4b_pointer < nmem_sp4b[bNums[first_sp4b_cluster[0]][neighbour_pointer]]; ++second_sp4b_pointer) {
                int second_sp4b_id = mem_sp4b[bNums[first_sp4b_cluster[0]][neighbour_pointer]][second_sp4b_pointer];
                int *second_sp4b_cluster = hcsp4b[second_sp4b_id];
                if (second_sp4b_id > first_sp4b_id) {

                    if (first_sp4b_cluster[4] != second_sp4b_cluster[4]) {
                        if (Bonds_BondCheck(first_sp4b_cluster[4], second_sp4b_cluster[4]) == 0) {

                            // Check that there are no shared particles between two clusters
                            for (k = 0; k < 5; ++k) {
                                for (l = 0; l < 5; ++l) {
                                    if (first_sp4b_cluster[k] == second_sp4b_cluster[l]) break;
                                }
                                if (l < 5) break;
                            }
                            if (k < 5) continue;

                            // Check there are two bonds to each ring particle
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

void Cluster_Write_10A(const int *first_sp4b_cluster, const int *second_sp4b_cluster) {
    int clusSize = 10;

    if (n10A == m10A) {
        hc10A = resize_2D_int(hc10A, m10A, m10A + incrStatic, clusSize, -1);
        m10A = m10A + incrStatic;
    }

    hc10A[n10A][0] = first_sp4b_cluster[0];
    hc10A[n10A][1] = first_sp4b_cluster[1];
    hc10A[n10A][2] = first_sp4b_cluster[2];
    hc10A[n10A][3] = first_sp4b_cluster[3];
    hc10A[n10A][4] = second_sp4b_cluster[0];
    hc10A[n10A][5] = second_sp4b_cluster[1];
    hc10A[n10A][6] = second_sp4b_cluster[2];
    hc10A[n10A][7] = second_sp4b_cluster[3];
    hc10A[n10A][8] = first_sp4b_cluster[4];
    hc10A[n10A][9] = second_sp4b_cluster[4];


    for(int i = 0; i < 8; i++) {
        if (s10A[hc10A[n10A][i]] == 'C') s10A[hc10A[n10A][i]] = 'B';
    }
    s10A[hc10A[n10A][8]] = 'O';
    s10A[hc10A[n10A][9]] = 'O';

    ++n10A;
}