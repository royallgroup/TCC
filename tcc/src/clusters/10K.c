#include "10K.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get10K() { // Detect 10K clusters
    // A 10K is a 9K with a SINGLE particle bonded to the common spindle of 9K.
    int bonded_to_spindle_id, num_extra_particles, extra_particle = 0;
    int parent_9K_id, id_9K_common;

    for (parent_9K_id=0; parent_9K_id < n9K; parent_9K_id++) {
        id_9K_common = hc9K[parent_9K_id][8];
        if (num_bonds[id_9K_common] < 10) {
            num_extra_particles = 0;
            for (bonded_to_spindle_id = 0; bonded_to_spindle_id < num_bonds[id_9K_common]; bonded_to_spindle_id++) {
                if (is_particle_in_9K(parent_9K_id, bNums[id_9K_common][bonded_to_spindle_id])) {
                    num_extra_particles++;
                    extra_particle = bNums[id_9K_common][bonded_to_spindle_id];
                }
            }
            if (num_extra_particles == 1) {
                Cluster_Write_10K(parent_9K_id, extra_particle);
            }
        }
    }
}

int is_particle_in_9K(int id_9K, int id_particle){
    // Returns 0 if particle is not in the specified 9K, else returns 1
    int i;

    for (i=0; i<9; i++) {
        if (id_particle == hc9K[id_9K][i]){
            return 0;
        }
    }
    return 1;
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