#include "7K.h"
#include "globals.h"
#include "tools.h"
#include "simple_cluster_methods.h"

void Clusters_Get7K() {

    //!  A 7K is made of two overlapping 5A clusters which have one common spindle.
    /*!
    *  Find 7K clusters
    *  7K is made from two overlapping 5A particles where:
    *  - 5Ai and 5Aj have one common spindle.
    *  - The other spindle of 5Ai is distinct from all the particles in 5Aj.
    *  - The other spindle of 5Aj is distinct from all the particles in 5Ai.
    *  - There are two common particles between the sp3 rings of 5Ai and 5Aj.
    *
    *  Cluster output: OOOOBBB
    *  Storage order: common spindle, other spindle x 2, common ring x 2, other ring x 2)
    */

    int common_spindle_id[2], other_spindle_ids[2], common_ring_ids[5], uncommon_ring_particles[2];

    for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
        int *first_5A_cluster = hcsp3c[first_5A_id];
        for (int first_5A_spindle_pointer = 3; first_5A_spindle_pointer < 5; ++first_5A_spindle_pointer) {
            for (int second_5A_pointer = 0; second_5A_pointer < nmem_sp3c[first_5A_cluster[first_5A_spindle_pointer]]; ++second_5A_pointer) {  // loop over all 5A_j common with spindle of 5A_i
                int second_5A_id = mem_sp3c[first_5A_cluster[first_5A_spindle_pointer]][second_5A_pointer];
                int *second_5A_cluster = hcsp3c[second_5A_id];
                if (second_5A_id > first_5A_id) {
                    // exactly one common spindle between 5A_i and 5A_j
                    if (count_common_spindle_particles(first_5A_cluster, second_5A_cluster, 5, 5, common_spindle_id) == 1) {

                        other_spindle_ids[0] = get_uncommon_spindle(first_5A_cluster, 5, common_spindle_id[0]);
                        other_spindle_ids[1] = get_uncommon_spindle(second_5A_cluster, 5, common_spindle_id[0]);

                        if (is_particle_in_cluster(second_5A_cluster, 5, other_spindle_ids[0]) == 0) {
                            if (is_particle_in_cluster(first_5A_cluster, 5, other_spindle_ids[1]) == 0) {

                                if (count_common_ring_particles(first_5A_cluster, second_5A_cluster, 3, 3, common_ring_ids) == 2) {

                                    count_uncommon_ring_particles(first_5A_cluster, second_5A_cluster, 3, 3, uncommon_ring_particles);

                                    Cluster_Write_7K(common_spindle_id[0], other_spindle_ids, common_ring_ids, uncommon_ring_particles);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Cluster_Write_7K(int common_spindle_id, int *other_spindle_ids, int *common_ring_particles, int *uncommon_ring_particles) {
    int clusSize = 7;

    if (n7K == m7K) {
        hc7K = resize_2D_int(hc7K, m7K, m7K + incrStatic, clusSize, -1);
        m7K = m7K + incrStatic;
    }

    hc7K[n7K][0] = common_spindle_id;
    hc7K[n7K][1] = other_spindle_ids[0];
    hc7K[n7K][2] = other_spindle_ids[1];
    hc7K[n7K][3] = common_ring_particles[0];
    hc7K[n7K][4] = common_ring_particles[1];
    hc7K[n7K][5] = uncommon_ring_particles[0];
    hc7K[n7K][6] = uncommon_ring_particles[1];

    s7K[hc7K[n7K][0]] = 'O';
    s7K[hc7K[n7K][1]] = 'O';
    s7K[hc7K[n7K][2]] = 'O';
    if (s7K[hc7K[n7K][3]] == 'C') s7K[hc7K[n7K][3]] = 'B';
    if (s7K[hc7K[n7K][4]] == 'C') s7K[hc7K[n7K][4]] = 'B';
    if (s7K[hc7K[n7K][5]] == 'C') s7K[hc7K[n7K][5]] = 'B';
    if (s7K[hc7K[n7K][6]] == 'C') s7K[hc7K[n7K][6]] = 'B';

    ++n7K;
}