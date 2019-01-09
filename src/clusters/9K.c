#include <globals.h>
#include <tools.h>
#include "9K.h"
#include "simple_cluster_methods.h"

//!  An 9K cluster is the intersection of two 6A clusters sharing a spindle and two ring particles.
/*!
*  Find 9K clusters
*  An 9K is 2 6A clusters where:
*      - There is one common spindle.
*      - The uncommon spindles are not in the opposite cluster.
*      - There are two common particles between the two 6A rings
*      - The uncommon sp4 ring particle are bonded
*
*  Cluster output: BBBBBBOOO
*  Storage order: common_ring_particles x 2, uncommon_ring_particles x 4, uncommon_spindles x 2, common_spindle
*/
void Clusters_Get9K() {
    int common_ring_particle_ids[5], common_spindle_ids[2], uncommon_spindle_particles[2], uncommon_ring_particles[4];

    for(int first_6A_id = 0; first_6A_id < nsp4c; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int spindle_pointer = 4; spindle_pointer < 6; spindle_pointer++) {
            for (int j = 0; j < nmem_sp4c[first_6A_cluster[spindle_pointer]]; ++j) {
                int second_6A_id = mem_sp4c[first_6A_cluster[spindle_pointer]][j];
                int *second_6A_cluster = hcsp4c[second_6A_id];
                if (second_6A_id > first_6A_id) {

                    if (count_common_spindle_particles(first_6A_cluster, second_6A_cluster, 6, 6, common_spindle_ids) == 1) {

                        uncommon_spindle_particles[0] = get_uncommon_spindle(first_6A_cluster, 6, common_spindle_ids[0]);
                        uncommon_spindle_particles[1] = get_uncommon_spindle(second_6A_cluster, 6, common_spindle_ids[0]);

                        if (is_particle_in_cluster(second_6A_cluster, 6, uncommon_spindle_particles[0]) == 0) {
                            if (is_particle_in_cluster(first_6A_cluster, 6, uncommon_spindle_particles[1]) == 0) {

                                if (count_common_particles(first_6A_cluster, second_6A_cluster, 4, 4,
                                                           common_ring_particle_ids) == 2) {

                                    count_uncommon_particles(first_6A_cluster, second_6A_cluster, 4, 4, uncommon_ring_particles);

                                    Cluster_Write_9k(common_ring_particle_ids, uncommon_ring_particles,
                                                     uncommon_spindle_particles, common_spindle_ids[0]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Cluster_Write_9k(const int *common_ring_particle_ids, const int *uncommon_ring_particles,
                      const int *uncommon_spindle_particles, int common_spindle_id) {
    int clusSize = 9;

    if(n9K == m9K) {
        hc9K= resize_2D_int(hc9K, m9K, m9K + incrStatic, clusSize, -1);
        m9K= m9K + incrStatic;
    }

    hc9K[n9K][0] = common_ring_particle_ids[0];
    hc9K[n9K][1] = common_ring_particle_ids[1];
    hc9K[n9K][2] = uncommon_ring_particles[0];
    hc9K[n9K][3] = uncommon_ring_particles[1];
    hc9K[n9K][4] = uncommon_ring_particles[2];
    hc9K[n9K][5] = uncommon_ring_particles[3];
    hc9K[n9K][6] = uncommon_spindle_particles[0];
    hc9K[n9K][7] = uncommon_spindle_particles[1];
    hc9K[n9K][8] = common_spindle_id;


    for(int i = 0; i < 6; i++) {
        if (s9K[hc9K[n9K][i]] == 'C') s9K[hc9K[n9K][i]] = 'B';
    }
    s9K[hc9K[n9K][6]] = 'O';
    s9K[hc9K[n9K][7]] = 'O';
    s9K[hc9K[n9K][8]] = 'O';

    ++n9K;
}