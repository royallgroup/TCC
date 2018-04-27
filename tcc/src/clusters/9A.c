#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "simple_cluster_methods.h"
#include "9A.h"

void Clusters_Get9A() {
    //!  A 9A is made of three overlapping sp4b clusters which each share 2 ring particles with each of the others.
    /*!
    *  Find 9A clusters
    *  9A is made from three overlapping sp4b clusters where:
    *  - All of the spindles are uncommon
    *  - Two particles are common between the rings of sp4b_i and sp4b_j
    *  - The ring of sp4b_j is made from the uncommon ring particles of sp4b_i and sp4b_j
    *
    *  Cluster output: OOOBBBBBB
    *  Storage order: common spindle, other spindle x 2, common ring x 2, other ring x 2)
    */

    int i_j_common_ring_particles[8], i_j_uncommon_ring_particles[8];

    for (int first_sp4b_id = 0; first_sp4b_id < nsp4b; ++first_sp4b_id) {
        int *first_sp4b_cluster = hcsp4b[first_sp4b_id];
        for (int second_sp4b_pointer = 0; second_sp4b_pointer < nmem_sp4b[first_sp4b_cluster[0]]; ++second_sp4b_pointer) {
            int second_sp4b_id = mem_sp4b[first_sp4b_cluster[0]][second_sp4b_pointer];
            int *second_sp4b_cluster = hcsp4b[second_sp4b_id];
            if (second_sp4b_id > first_sp4b_id) {
                // Spindles should not be common
                if (first_sp4b_cluster[4] != second_sp4b_cluster[4]) {
                    // Spindles should not be bonded
                    if (Bonds_BondCheck(first_sp4b_cluster[4], second_sp4b_cluster[4]) == 0) {

                        // Check there are two common particles between the rings of the two clusters
                        if (count_common_ring_particles(first_sp4b_cluster, second_sp4b_cluster, 4, 4, i_j_common_ring_particles) == 2) {

                            if (count_bonded_ring_particles(first_sp4b_cluster, second_sp4b_cluster, i_j_common_ring_particles, i_j_uncommon_ring_particles) == 4) {

                                // POSSIBLE IMPROVEMENT!! could make detection faster here by looping over sp4b clusters bonded to uncommon particles to sp4b_i
                                for (int third_sp4b_id = second_sp4b_id + 1; third_sp4b_id < nsp4b; third_sp4b_id++) {
                                    int *third_sp4b_cluster = hcsp4b[third_sp4b_id];
                                    // ERROR!! need to check spindle of sp4b_k distinct from spindles of sp4b_i and sp4b_j and no bonds between these three particles

                                    int tmp[4];
                                    if (count_common_ring_particles(third_sp4b_cluster, i_j_uncommon_ring_particles, 4, 4, tmp) == 4) {

                                        Cluster_Write_9A(first_sp4b_cluster, second_sp4b_cluster, third_sp4b_cluster, i_j_common_ring_particles, i_j_uncommon_ring_particles);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int count_bonded_ring_particles(const int *first_sp4b_cluster, const int *second_sp4b_cluster, const int *db, int *ob) {

    int num_bonded_ring_particles = 0;
    for (int first_cluster_pointer = 0; first_cluster_pointer < 4; ++first_cluster_pointer) {
        // find particles in SP4 ring of sp4b_i not common to SP4 ring of sp4b_j
        if (first_sp4b_cluster[first_cluster_pointer] != db[0] && first_sp4b_cluster[first_cluster_pointer] != db[1]) {
            for (int second_cluster_pointer = 0; second_cluster_pointer < 4; ++second_cluster_pointer) {
                // find particles in SP4 ring of sp4b_j not common to SP4 ring of sp4b_i
                if (second_sp4b_cluster[second_cluster_pointer] != db[0] && second_sp4b_cluster[second_cluster_pointer] != db[1]) {
                    // check non-common SP4 ring particles from sp4b_i and sp4b_j are bonded
                    if (Bonds_BondCheck(first_sp4b_cluster[first_cluster_pointer], second_sp4b_cluster[second_cluster_pointer])) {
                        ob[num_bonded_ring_particles] = first_sp4b_cluster[first_cluster_pointer];
                        num_bonded_ring_particles++;
                        ob[num_bonded_ring_particles] = second_sp4b_cluster[second_cluster_pointer];
                        num_bonded_ring_particles++;
                    }
                }
            }
        }
    }
    return num_bonded_ring_particles;
}

void Cluster_Write_9A(int *first_sp4b_cluster, int *second_sp4b_cluster, int *third_sp4b_cluster,
                      int *i_j_common_ring_particles, int *i_j_uncommon_ring_particles) {

    int clusSize = 9;

    if (n9A == m9A) {
        hc9A = resize_2D_int(hc9A, m9A, m9A + incrStatic, clusSize, -1);
        m9A = m9A + incrStatic;
    }

    hc9A[n9A][0] = first_sp4b_cluster[4];
    hc9A[n9A][1] = second_sp4b_cluster[4];
    hc9A[n9A][2] = third_sp4b_cluster[4];
    hc9A[n9A][3] = i_j_common_ring_particles[0];
    hc9A[n9A][4] = i_j_common_ring_particles[1];
    hc9A[n9A][5] = i_j_uncommon_ring_particles[0];
    hc9A[n9A][6] = i_j_uncommon_ring_particles[1];
    hc9A[n9A][7] = i_j_uncommon_ring_particles[2];
    hc9A[n9A][8] = i_j_uncommon_ring_particles[3];

    s9A[hc9A[n9A][0]] = 'O';
    s9A[hc9A[n9A][1]] = 'O';
    s9A[hc9A[n9A][2]] = 'O';
    for (int i = 3; i < 9; ++i) {
        if (s9A[hc9A[n9A][i]] == 'C') s9A[hc9A[n9A][i]] = 'B';
    }

    ++n9A;
}