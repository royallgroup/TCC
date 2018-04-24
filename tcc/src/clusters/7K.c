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

    int first_5A_id, first_5A_spindle_pointer;
    int second_5A_id, second_5A_pointer;
    int *first_5A_cluster, *second_5A_cluster;
    int common_spindle_id, other_spindle_ids[2], common_ring_ids[2], uncommon_ring_particles[2];

    for (first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {  // loop over all 5A_i
        first_5A_cluster = hcsp3c[first_5A_id];
        for (first_5A_spindle_pointer = 3;
             first_5A_spindle_pointer < 5; ++first_5A_spindle_pointer) {    // loop over both spindles of 5A_i
            for (second_5A_pointer = 0; second_5A_pointer < nmem_sp3c[first_5A_cluster[first_5A_spindle_pointer]]; ++second_5A_pointer) {  // loop over all 5A_j common with spindle of 5A_i
                second_5A_id = mem_sp3c[first_5A_cluster[first_5A_spindle_pointer]][second_5A_pointer];
                second_5A_cluster = hcsp3c[second_5A_id];
                if (second_5A_id > first_5A_id) {
                    // exactly one common spindle between 5A_i and 5A_j
                    if (count_common_spindles_between_5As(first_5A_cluster, second_5A_cluster, &common_spindle_id) == 1) {

                        get_other_spindle_ids(first_5A_cluster, second_5A_cluster, common_spindle_id, other_spindle_ids);

                        if (is_particle_in_5A(second_5A_cluster, other_spindle_ids[0]) == 0) {
                            if (is_particle_in_5A(first_5A_cluster, other_spindle_ids[1]) == 0) {

                                if (count_common_ring_particles(first_5A_cluster, second_5A_cluster, common_ring_ids) == 2) {

                                    uncommon_ring_particles[0] = get_uncommon_ring_particle(first_5A_cluster, common_ring_ids);
                                    uncommon_ring_particles[1] = get_uncommon_ring_particle(second_5A_cluster, common_ring_ids);

                                    Cluster_Write_7K(common_spindle_id, other_spindle_ids, common_ring_ids, uncommon_ring_particles);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void get_other_spindle_ids(const int *first_5A_cluster, const int *second_5A_cluster, int common_spindle_id, int *other_spindle_ids) {
    if (first_5A_cluster[3] == common_spindle_id) {
        other_spindle_ids[0] = first_5A_cluster[4];
    } else {
        other_spindle_ids[0] = first_5A_cluster[3];
    }
    if (second_5A_cluster[3] == common_spindle_id) {
        other_spindle_ids[1] = second_5A_cluster[4];
    } else {
        other_spindle_ids[1] = second_5A_cluster[3];
    }
}

int is_particle_in_5A(const int *five_A_cluster, int particle_id) {
    for (int five_A_pointer = 0; five_A_pointer < 5; five_A_pointer++) {
        if (particle_id == five_A_cluster[five_A_pointer]) {
            return 1;
        }
    }
    return 0;
}

int count_common_ring_particles(const int *first_5A_cluster, const int *second_5A_cluster, int *common_ring_ids) {
    int num_common_ring_particles = 0;
    int common_particle_ids[3];

    for (int first_5A_pointer = 0; first_5A_pointer < 3; first_5A_pointer++) {
        for (int second_5A_pointer = 0; second_5A_pointer < 3; second_5A_pointer++) {
            if (first_5A_cluster[first_5A_pointer] == second_5A_cluster[second_5A_pointer]) {
                common_particle_ids[num_common_ring_particles] = first_5A_cluster[first_5A_pointer];
                num_common_ring_particles++;
            }
        }
    }
    if (num_common_ring_particles == 2) {
        common_ring_ids[0] = common_particle_ids[0];
        common_ring_ids[1] = common_particle_ids[1];
    }
    return num_common_ring_particles;
}

int get_uncommon_ring_particle(const int *first_5A_cluster, const int *common_ring_ids) {
    for (int first_5A_pointer = 0; first_5A_pointer < 3; first_5A_pointer++) {
        if (first_5A_cluster[first_5A_pointer] != common_ring_ids[0] &&
            first_5A_cluster[first_5A_pointer] != common_ring_ids[1]) {
            return (first_5A_cluster[first_5A_pointer]);
        }
    }
    Error("uncommon ring particle not found.");
    return 0;
}

void Cluster_Write_7K(int scom, int *sother, int *sp3_com, int *uncommon_ring_particles) {
    int clusSize = 7;

    if (n7K == m7K) {
        hc7K = resize_2D_int(hc7K, m7K, m7K + incrStatic, clusSize, -1);
        m7K = m7K + incrStatic;
    }
    // Now we have found the 7K cluster

    hc7K[n7K][0] = scom;
    hc7K[n7K][1] = sother[0];
    hc7K[n7K][2] = sother[1];
    hc7K[n7K][3] = sp3_com[0];
    hc7K[n7K][4] = sp3_com[1];
    hc7K[n7K][5] = uncommon_ring_particles[0];
    hc7K[n7K][6] = uncommon_ring_particles[1];

    quickSort(&hc7K[n7K][1], 2);
    quickSort(&hc7K[n7K][3], 2);
    quickSort(&hc7K[n7K][5], 2);


    // hc7K key: (scom, sother, ring_com, ring_other)
    s7K[hc7K[n7K][0]] = 'O';
    s7K[hc7K[n7K][1]] = 'O';
    s7K[hc7K[n7K][2]] = 'O';
    if (s7K[hc7K[n7K][3]] == 'C') s7K[hc7K[n7K][3]] = 'B';
    if (s7K[hc7K[n7K][4]] == 'C') s7K[hc7K[n7K][4]] = 'B';
    if (s7K[hc7K[n7K][5]] == 'C') s7K[hc7K[n7K][5]] = 'B';
    if (s7K[hc7K[n7K][6]] == 'C') s7K[hc7K[n7K][6]] = 'B';

    ++n7K;
}