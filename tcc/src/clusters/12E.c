#include <globals.h>
#include <tools.h>
#include "12E.h"

int Clusters_Get12E(int j) {
    // Return 1 if 11F is also 12E
    // An 11F with an extra 5A bonded
    int new_5A_id, new_5A_ring_pointer, m;
    int num_common_particles, common_particle_ids[4];
    int uncommon_sp3_ring_particle;
    int* new_5A_cluster;


    for(new_5A_id = j + 1; new_5A_id < nsp3c; new_5A_id++) { // loop through 5A clusters with ID > j
        new_5A_cluster = hcsp3c[new_5A_id];

        // Check spindles of new 5A are common with uncommon spindles of 6As in 11F
        if (new_5A_cluster[3] == hc11F[n11F][1] || new_5A_cluster[3] == hc11F[n11F][2]) {
            if (new_5A_cluster[4] == hc11F[n11F][1] || new_5A_cluster[4] == hc11F[n11F][2]) {
                // Check two new 5A ring particles are common with 11F ring particles
                num_common_particles = 0;
                for (new_5A_ring_pointer = 0; new_5A_ring_pointer < 3; new_5A_ring_pointer++) { // loop through sp3 ring particles
                    for (m = 7; m < 11; m++) { // loop through ring particles of 11F
                        if(new_5A_cluster[new_5A_ring_pointer] == hc11F[n11F][m]) {
                            common_particle_ids[num_common_particles] = new_5A_cluster[new_5A_ring_pointer];
                            num_common_particles++;
                        }
                    }
                }
                if (num_common_particles == 2) {
                    // Find which of the 3 ring particles in new_5A is the uncommon one
                    uncommon_sp3_ring_particle = get_uncommon_5A_ring_particle(common_particle_ids, new_5A_cluster);
                    Raw_Write_12E(uncommon_sp3_ring_particle);
                    Cluster_Write_12E();

                    return 1;
                }
            }
        }
    }
    return 0;
}

int get_uncommon_5A_ring_particle(const int *common, const int *new_5A_cluster) {
    int i;

    for (i = 0; i < 3; i++) {
        if (new_5A_cluster[i] != common[0] && new_5A_cluster[i] != common[1]) {
            return(new_5A_cluster[i]);
        }
    }
}

void Raw_Write_12E(int uncommon_sp3_ring_particle) {// now we have found the 12E
    int i;
    int clusSize = 12;

    if (n12E == m12E) {
        hc12E = resize_2D_int(hc12E, m12E, m12E + incrStatic, clusSize, -1);
        m12E = m12E + incrStatic;
    }

    hc12E[n12E][11] = uncommon_sp3_ring_particle;
    for (i = 0; i < 10; ++i) {
        hc12E[n12E][i] = hc11F[n11F][i + 1];
    }
    hc12E[n12E][10] = hc11F[n11F][0];
    quickSort(&hc12E[n12E][0], 6);
    quickSort(&hc12E[n12E][6], 6);


}

void Cluster_Write_12E() {
    int i;
    for(i=0; i<6; i++){
        s12E[hc12E[n12E][i]] = 'O';
    }
    for(i=6; i<12; i++){
        if(s12E[hc12E[n12E][i]] == 'C') s12E[hc12E[n12E][i]] = 'B';
    }
}