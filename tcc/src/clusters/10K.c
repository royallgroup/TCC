#include "simple_cluster_methods.h"
#include "10K.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get10K() {

    //!  An 10K cluster is 9K with a SINGLE particle bonded to the common spindle of 9K.
    /*!
   *  Find 10K clusters
   *  An 10K is a 9K with one extra particle where:
   *      -  The extra particle is bonded to common spindle of 9K.
   *      -  The 9K common spindle has no other extra neighbours.
   *
   *  Cluster output: BBBBBOOOO
   *  Storage order: as_for_9K x 9, extra_particle
   */

    int extra_particle_id;

    for (int first_9K_id = 0; first_9K_id < n9K; first_9K_id++) {
        int *first_9K_cluster = hc9K[first_9K_id];
        int first_9K_common_id = first_9K_cluster[8];

        if (num_bonds[first_9K_common_id] < 10) {
            if (count_extra_particles(first_9K_cluster, first_9K_common_id, &extra_particle_id) == 1) {
                Cluster_Write_10K(first_9K_id, extra_particle_id);
            }
        }
    }
}

int count_extra_particles(const int *first_9K_cluster, int first_9K_common_id, int *extra_particle_id) {
    int num_extra_particles = 0;

    for (int bonded_particle_pointer = 0; bonded_particle_pointer < num_bonds[first_9K_common_id]; bonded_particle_pointer++) {
        int bonded_particle_id = bNums[first_9K_common_id][bonded_particle_pointer];
        if (is_particle_in_cluster(first_9K_cluster, 9, bonded_particle_id) == 0) {
            num_extra_particles++;
            (*extra_particle_id) = bonded_particle_id;
        }
    }
    return num_extra_particles;
}

void Cluster_Write_10K(int id_9k, int extra_particle) {

    int clusSize=10;

    if(n10K == m10K) {
        hc10K= resize_2D_int(hc10K, m10K, m10K + incrStatic, clusSize, -1);
        m10K= m10K + incrStatic;
    }

    for(int i = 0; i < 9; i++) {
        hc10K[n10K][i] = hc9K[id_9k][i];
    }
    hc10K[n10K][9] = extra_particle;

    for(int i = 0; i < 6; i++) {
        if (s10K[hc10K[n10K][i]] == 'C') s10K[hc10K[n10K][i]] = 'B';
    }
    s10K[hc10K[n10K][6]] = 'O';
    s10K[hc10K[n10K][7]] = 'O';
    s10K[hc10K[n10K][8]] = 'O';
    s10K[hc10K[n10K][9]] = 'O';

    n10K++;
}