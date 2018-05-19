#include "simple_cluster_methods.h"
#include "11A.h"
#include "globals.h"
#include "tools.h"
#include "bonds.h"

void Clusters_Get11A() {

    //!  An 11A cluster is the intersection two 6A clusters with a common spindle
    /*!
   *  Find 11A clusters
   *  An 11A is two 6A clusters where:
   *      - There is one common spindle particle
   *      - All other particles are distinct.
   *      - Each particle in each ring is bonded to exactly 2 particles in the other ring
   *
   *  Cluster output: BBBBBBBBOOS
   *  Storage order: sp4_i x 4, sp4_j x 4, uncommon_spindles x 2, common_spindle
   */
    int common_spindle_id[2], uncommon_spindle_ids[2];
    int common_ring_particles[2];

    for (int first_6A_id = 0; first_6A_id < nsp4c; first_6A_id++) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int first_6A_spindle_pointer = 0; first_6A_spindle_pointer < 2; first_6A_spindle_pointer++) {
            int first_6A_spindle_ID = first_6A_cluster[first_6A_spindle_pointer + 4];
            for (int mem_pointer = 0; mem_pointer < nmem_sp4c[first_6A_spindle_ID]; mem_pointer++) {
                int second_6A_id = mem_sp4c[first_6A_spindle_ID][mem_pointer];
                if (second_6A_id > first_6A_id) {
                    int *second_6A_cluster = hcsp4c[second_6A_id];
                    if (count_common_spindle_particles(first_6A_cluster, second_6A_cluster, 6, 6, common_spindle_id) == 1) {
                        uncommon_spindle_ids[0] = get_uncommon_spindle(first_6A_cluster, 6, common_spindle_id[0]);
                        uncommon_spindle_ids[1] = get_uncommon_spindle(second_6A_cluster, 6, common_spindle_id[0]);
                        if (count_common_ring_particles(first_6A_cluster, second_6A_cluster, 4, 4, common_ring_particles) == 0) {
                            if (Check_6A_rings_bonded(first_6A_cluster, second_6A_cluster) == 1) {
                                Cluster_Write_11A(first_6A_cluster, second_6A_cluster, uncommon_spindle_ids, common_spindle_id[0]);
                            }
                        }
                    }
                }
            }
        }
    }
}

int Check_6A_rings_bonded(const int *first_6A_cluster, const int *second_6A_cluster) {
    int first_ring_pointer, second_ring_pointer, num_bonds;
    // Check if there are two bonds between each particle in ring 1 and particles in ring 2
    // Returns 1 if all ring 1 particles have 2 bonds to ring 2 particles, return 0 if not

    for(first_ring_pointer = 0; first_ring_pointer < 4; first_ring_pointer++) {
        num_bonds = 0;
        for(second_ring_pointer = 0; second_ring_pointer < 4; second_ring_pointer++) {
            if(Bonds_BondCheck(first_6A_cluster[first_ring_pointer], second_6A_cluster[second_ring_pointer])) {
                num_bonds++;
            }
        }
        if(num_bonds != 2) {
            return 0;
        }
    }
    return 1;
}

void Cluster_Write_11A(const int *first_6A, const int *second_6A, const int *sother, const int scom) {
    int clusSize=11;

    if (n11A == m11A) {
        hc11A = resize_2D_int(hc11A, m11A, m11A + incrStatic, clusSize, -1);
        m11A = m11A + incrStatic;
    }

    hc11A[n11A][0] = first_6A[0];
    hc11A[n11A][1] = first_6A[1];
    hc11A[n11A][2] = first_6A[2];
    hc11A[n11A][3] = first_6A[3];
    hc11A[n11A][4] = second_6A[0];
    hc11A[n11A][5] = second_6A[1];
    hc11A[n11A][6] = second_6A[2];
    hc11A[n11A][7] = second_6A[3];
    hc11A[n11A][8] = sother[0];
    hc11A[n11A][9] = sother[1];
    hc11A[n11A][10] = scom;

    for(int i = 0; i < 8; i++) {
        if (s11A[hc11A[n11A][i]] == 'C') s11A[hc11A[n11A][i]] = 'B';
    }
    if (s11A[hc11A[n11A][8]] != 'S') s11A[hc11A[n11A][8]] = 'O';
    if (s11A[hc11A[n11A][9]] != 'S') s11A[hc11A[n11A][9]] = 'O';
    s11A[hc11A[n11A][10]] = 'S';
    ++n11A;
}