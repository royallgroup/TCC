#include "simple_cluster_methods.h"
#include "8A.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "string.h"

//!  An 8A cluster is one of 3 possible topological combinations of sp5b/c clusters.
/*!
*  Find 8A clusters
*  There are 3 methods used for 8A detection
*  - A pair of sp5b where the spindles are distinct and share 4 ring particles
*  - A pair of 7A clusters where
*      - Both 7Ai spindle particles are common with the 7Aj spindles.
*      - There are four common particles between sp5 rings of 7Ai and 7Aj .
*  - A sp5b cluster and a 7A cluster where:
*      - One 7A spindle is common with the sp5b spindle.
*      - The other 7A spindle is distinct from all the sp5b particles.
*      - There are four common particles between sp5 rings of sp5b and 7A.
*
*  Cluster output: BBBBBBBB
*  Storage order: particles ordered by id x 8)
*/
void Clusters_Get8A() {
    get_first_8A_type();
    get_second_8A_type();
    get_third_8A_type();
}

void get_first_8A_type() {

    int uncommon_ring_particle_ids[10];
    int common_ring_particle_ids[5];
    int trial[8];

    for (int first_sp5b_id = 0; first_sp5b_id < nsp5b; ++first_sp5b_id) {
        int *first_sp5b_cluster = hcsp5b[first_sp5b_id];

        for (int first_sp5b_ring_pointer = 0; first_sp5b_ring_pointer < 5; ++first_sp5b_ring_pointer) {
            for (int second_sp5b_pointer = 0; second_sp5b_pointer < nmem_sp5b[first_sp5b_cluster[first_sp5b_ring_pointer]]; ++second_sp5b_pointer) {
                int second_sp5b_id = mem_sp5b[first_sp5b_cluster[first_sp5b_ring_pointer]][second_sp5b_pointer];
                int *second_sp5b_cluster = hcsp5b[second_sp5b_id];

                if (second_sp5b_id > first_sp5b_id) {

                    // Check rings share 4 particles
                    if (count_common_ring_particles(first_sp5b_cluster, second_sp5b_cluster, 5, 5, common_ring_particle_ids) == 4) {

                        // Check for distinct spindles
                        if (first_sp5b_cluster[5] != second_sp5b_cluster[5]) {

                            count_uncommon_ring_particles(first_sp5b_cluster, second_sp5b_cluster, 5, 5, uncommon_ring_particle_ids);

                            // build up trial cluster
                            trial[0] = first_sp5b_cluster[5];
                            trial[1] = second_sp5b_cluster[5];
                            trial[2] = uncommon_ring_particle_ids[0];
                            trial[3] = uncommon_ring_particle_ids[1];
                            trial[4] = common_ring_particle_ids[0];
                            trial[5] = common_ring_particle_ids[1];
                            trial[6] = common_ring_particle_ids[2];
                            trial[7] = common_ring_particle_ids[3];

                            quickSort(&trial[0], 8);

                            if (check_unique_cluster(trial, 8, hc8A, n8A) == 0) {
                                Cluster_Write_8A(trial);
                            }
                        }
                    }
                }
            }
        }
    }
}

void get_second_8A_type() {

    int uncommon_ring_particle_ids[10];
    int common_particle_ids[5];
    int spindle_ids[2];
    int trial[8];

    for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (int first_7A_spindle_pointer = 5; first_7A_spindle_pointer < 6; ++first_7A_spindle_pointer) {
            for (int second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[first_7A_cluster[first_7A_spindle_pointer]]; ++second_7A_pointer) {
                int second_7A_id = mem_sp5c[first_7A_cluster[first_7A_spindle_pointer]][second_7A_pointer];
                int *second_7A_cluster = hcsp5c[second_7A_id];

                if (second_7A_id <= first_7A_id) {

                    // exactly four members of the SP5 rings of 7A_i and 7A_j in common
                    if (count_common_ring_particles(first_7A_cluster, second_7A_cluster, 5, 5, common_particle_ids) == 4) {

                        if (count_common_spindle_particles(first_7A_cluster, second_7A_cluster, 7, 7, spindle_ids) == 2) {

                            count_uncommon_ring_particles(first_7A_cluster, second_7A_cluster, 5, 5, uncommon_ring_particle_ids);
                            // build up trial cluster
                            trial[0] = first_7A_cluster[5];
                            trial[1] = first_7A_cluster[6];
                            trial[2] = uncommon_ring_particle_ids[0];
                            trial[3] = uncommon_ring_particle_ids[1];
                            trial[4] = common_particle_ids[0];
                            trial[5] = common_particle_ids[1];
                            trial[6] = common_particle_ids[2];
                            trial[7] = common_particle_ids[3];

                            quickSort(&trial[0], 8);

                            if (check_unique_cluster(trial, 8, hc8A, n8A) == 0) {
                                Cluster_Write_8A(trial);
                            }
                        }
                    }
                }
            }
        }
    }
}

void get_third_8A_type() {

    int uncommon_ring_particle_ids[10];
    int common_ring_particle_ids[5];
    int trial[8];

    for (int first_sp5b_id = 0; first_sp5b_id < nsp5b; ++first_sp5b_id) {
        int *first_sp5b_cluster = hcsp5b[first_sp5b_id];
        for (int first_sp5b_spindle_pointer = 5; first_sp5b_spindle_pointer < 6; ++first_sp5b_spindle_pointer) {
            for (int first_7A_id = 0; first_7A_id < nmem_sp5c[first_sp5b_cluster[first_sp5b_spindle_pointer]]; ++first_7A_id) {
                int *first_7A_cluster = hcsp5c[mem_sp5c[first_sp5b_cluster[first_sp5b_spindle_pointer]][first_7A_id]];

                if (count_common_ring_particles(first_sp5b_cluster, first_7A_cluster, 5, 5, common_ring_particle_ids) == 4) {

                    // sp5b_i spindle common with one of 7A_j spindles
                    if (first_sp5b_cluster[5] == first_7A_cluster[5] || first_sp5b_cluster[5] == first_7A_cluster[6]) {

                        count_uncommon_ring_particles(first_sp5b_cluster, first_7A_cluster, 5, 5, uncommon_ring_particle_ids);

                        // build up trial cluster
                        trial[0] = first_7A_cluster[5];
                        trial[1] = first_7A_cluster[6];
                        trial[2] = uncommon_ring_particle_ids[0];
                        trial[3] = uncommon_ring_particle_ids[1];
                        trial[4] = common_ring_particle_ids[0];
                        trial[5] = common_ring_particle_ids[1];
                        trial[6] = common_ring_particle_ids[2];
                        trial[7] = common_ring_particle_ids[3];

                        quickSort(&trial[0], 8);

                        if (check_unique_cluster(trial, 8, hc8A, n8A) == 0) {
                            Cluster_Write_8A(trial);
                        }
                    }
                }
            }
        }
    }
}


void Cluster_Write_8A(const int *trial) {// hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
    int clusSize = 8;

    if (n8A == m8A) {
        hc8A = resize_2D_int(hc8A, m8A, m8A + incrStatic, clusSize, -1);
        m8A = m8A + incrStatic;
    }

    for (int i = 0; i < 8; ++i) {
        hc8A[n8A][i] = trial[i];
    }

    for (int i = 0; i < 8; ++i) {
        s8A[hc8A[n8A][i]] = 'B';
    }

    ++n8A;
}
