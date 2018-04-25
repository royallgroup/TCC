#include "simple_cluster_methods.h"
#include "10K.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get10K() { // Detect 10K clusters
    // A 10K is a 9K with a SINGLE particle bonded to the common spindle of 9K.
    int bonded_particle_pointer, num_extra_particles, extra_particle = 0;
    int parent_9K_id, parent_9K_common_id;

    for (parent_9K_id = 0; parent_9K_id < n9K; parent_9K_id++) {
        int *parent_9K_cluster = hc9K[parent_9K_id];
        parent_9K_common_id = parent_9K_cluster[8];
        if (num_bonds[parent_9K_common_id] < 10) {
            num_extra_particles = 0;
            for (bonded_particle_pointer = 0; bonded_particle_pointer < num_bonds[parent_9K_common_id]; bonded_particle_pointer++) {
                int bonded_particle_id = bNums[parent_9K_common_id][bonded_particle_pointer];
                if (is_particle_in_cluster(parent_9K_cluster, 9, bonded_particle_id) == 0) {
                    num_extra_particles++;
                    extra_particle = bonded_particle_id;
                }
            }
            if (num_extra_particles == 1) {
                Cluster_Write_10K(parent_9K_id, extra_particle);
            }
        }
    }
}

void Cluster_Write_10K(int id_9k, int extra_particle) {
    // hc10K key: (common_SP4_1, common_SP4_2, other_SP4*4, other_spindle_1, other_spindle_2, scom, ep)

    int i;
    int clusSize=10;

    if(n10K == m10K) {
        hc10K= resize_2D_int(hc10K, m10K, m10K + incrStatic, clusSize, -1);
        m10K= m10K + incrStatic;
    }

    for(i=0; i<9; i++) {
        hc10K[n10K][i]= hc9K[id_9k][i];
    }
    hc10K[n10K][9]=extra_particle;

    for(i=0; i<6; i++) {
        if (s10K[hc10K[n10K][i]] == 'C') s10K[hc10K[n10K][i]] = 'B';
    }
    s10K[hc10K[n10K][6]] = 'O';
    s10K[hc10K[n10K][7]] = 'O';
    s10K[hc10K[n10K][8]] = 'O';
    s10K[hc10K[n10K][9]] = 'O';

    n10K++;
}