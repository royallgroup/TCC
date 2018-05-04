#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include <clusters/10B.h>
#include <clusters/11B.h>
#include <clusters/11E.h>
#include "simple_cluster_methods.h"
#include "9B.h"

void Clusters_Get9B_10B_11B_11E_12D() {

    //!  An 9B cluster is the intersection of 7A clusters shraing a spindle and two ring particles.
    /*!
    *  Find 9B clusters
    *  An 9B is 2 7A clusters where:
    *      - There is one common spindle.
    *      - The uncommon spindles are bonded.
    *      - Each distinct spindle is common with a ring particle of the other 7A
    *      - There are two common particles between the two 7A rings
    *
    *  Cluster output: BBBBBBOOS
    *  Storage order: uncommon_ring_particles x 4, common_ring_particles x 2, uncommon_spindles x 2, common_spindle
    */

    int common_spindle_ids[2];
    int uncommon_spindles[2];
    int common_ring_particles[2];
    int uncommon_ring_particles[10];

    for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            int spindle_particle_id = first_7A_cluster[spindle_pointer + 5];
            for (int second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[spindle_particle_id]; ++second_7A_pointer) {
                int second_7A_id = mem_sp5c[spindle_particle_id][second_7A_pointer];
                int *second_7A_cluster = hcsp5c[second_7A_id];

                if (count_common_spindle_particles(first_7A_cluster, second_7A_cluster, 7, 7, common_spindle_ids) == 1) {
                    uncommon_spindles[0] = get_uncommon_spindle(first_7A_cluster, 7, common_spindle_ids[0]);
                    uncommon_spindles[1] = get_uncommon_spindle(second_7A_cluster, 7, common_spindle_ids[0]);


                    if (is_particle_in_cluster(second_7A_cluster, 7, uncommon_spindles[0]) == 1) {
                        if (is_particle_in_cluster(first_7A_cluster, 7, uncommon_spindles[1]) == 1) {

                            if (count_common_ring_particles(first_7A_cluster, second_7A_cluster, 7, 7, common_ring_particles) == 2) {

                                count_uncommon_ring_particles(first_7A_cluster, second_7A_cluster, 5, 5, uncommon_ring_particles);

                                Cluster_Write_9B(uncommon_ring_particles, common_ring_particles, uncommon_spindles, common_spindle_ids[0]);

                                if (do10B == 1) Clusters_Get10B(first_7A_id, second_7A_id);
                                if (do11B == 1) {
                                    if (Clusters_Get11B()) {
                                        s11B[hc9B[n9B][8]] = 'S';
                                        ++n11B;
                                    }
                                }
                                if (do11E == 1)
                                    Clusters_Get11E_12D(first_7A_id, second_7A_id, common_spindle_ids[0], uncommon_spindles[0], uncommon_spindles[1]);

                                ++n9B;
                            }
                        }
                    }
                }
            }
        }
    }
}

void Cluster_Write_9B(const int* uncommon_ring_particles, const int *common_ring_particles,
                      const int *uncommon_spindles, int common_spindle_id) {
    int clusSize = 9;

    // Now we have found the 9B C2v cluster
    if (n9B == m9B) {
        hc9B = resize_2D_int(hc9B, m9B, m9B + incrStatic, clusSize, -1);
        m9B = m9B + incrStatic;
    }
    hc9B[n9B][0] = uncommon_ring_particles[1];
    hc9B[n9B][1] = uncommon_ring_particles[2];
    hc9B[n9B][2] = uncommon_ring_particles[3];
    hc9B[n9B][3] = uncommon_ring_particles[4];
    hc9B[n9B][4] = common_ring_particles[0];
    hc9B[n9B][5] = common_ring_particles[1];
    hc9B[n9B][6] = uncommon_spindles[0];
    hc9B[n9B][7] = uncommon_spindles[1];
    hc9B[n9B][8] = common_spindle_id;

    for(int i = 0; i < 6; i++) {
        if (s9B[hc9B[n9B][i]] == 'C') s9B[hc9B[n9B][i]] = 'B';
    }
    if (s9B[hc9B[n9B][6]] != 'S') s9B[hc9B[n9B][6]] = 'O';
    if (s9B[hc9B[n9B][7]] != 'S') s9B[hc9B[n9B][7]] = 'O';
    s9B[hc9B[n9B][8]] = 'S';
}