#include <globals.h>
#include <tools.h>
#include "12E.h"

void Clusters_Get12E() {
    //!  An 12E cluster is the intersection of an 11F and 5A cluster.
    /*!
   *  Find 12E clusters
   *  A 12E cluster is an 11F and 5A cluster where:
   *      - The spindle atoms of the 5A cluster are common with the uncommon spindle atoms of the 6A clusters constituting the 11F cluster.
   *      - Of the SP3 ring particles in the 5A cluster, two are common with rd1 and rd2 from the 11F cluster, and one is new.
   *
   *  Cluster output: BOOOOOOBBBBB
   *  Storage order: as for 11F x 11, new_particle)
   */
    int new_5A_ring_pointer, m;
    int num_common_particles, common_particle_ids[4];
    int uncommon_sp3_ring_particle;
    int common_particle_overlap;

    for (int first_11F_id = 0; first_11F_id < n11F; first_11F_id++) {
        int *first_11F_cluster = hc11F[first_11F_id];
        for (int first_5A_id = 0; first_5A_id < nsp3c; first_5A_id++) { // loop through 5A clusters with ID > j
            int *first_5A_cluster = hcsp3c[first_5A_id];

            // Check spindles of new 5A are common with uncommon spindles of 6As in 11F
            if (first_5A_cluster[3] == first_11F_cluster[1] || first_5A_cluster[3] == first_11F_cluster[2]) {
                if (first_5A_cluster[4] == first_11F_cluster[1] || first_5A_cluster[4] == first_11F_cluster[2]) {
                    // Check two new 5A ring particles are common with 11F ring particles
                    num_common_particles = 0;
                    for (new_5A_ring_pointer = 0; new_5A_ring_pointer < 3; new_5A_ring_pointer++) { // loop through sp3 ring particles
                        for (m = 7; m < 11; m++) { // loop through ring particles of 11F
                            if (first_5A_cluster[new_5A_ring_pointer] == first_11F_cluster[m]) {
                                common_particle_ids[num_common_particles] = first_5A_cluster[new_5A_ring_pointer];
                                num_common_particles++;
                            }
                        }
                    }
                    if (num_common_particles == 2) {
                        // Find which of the 3 ring particles in new_5A is the uncommon one
                        uncommon_sp3_ring_particle = get_uncommon_5A_ring_particle(common_particle_ids, first_5A_cluster);
                        // check that new particle is not alread in 11F
                        common_particle_overlap = 0;
                        for (int i = 0; i < 11; i++) {
                            if (uncommon_sp3_ring_particle == first_11F_cluster[i]) {
                                common_particle_overlap = 1;
                                break;
                            }
                        }
                        if (common_particle_overlap == 0) {
                            Raw_Write_12E(first_11F_cluster, uncommon_sp3_ring_particle);
                            Cluster_Write_12E();
                        }
                    }
                }
            }
        }
    }
}

int get_uncommon_5A_ring_particle(const int *common_particle_ids, const int *new_5A_cluster) {

    for (int i = 0; i < 3; i++) {
        if (new_5A_cluster[i] != common_particle_ids[0] && new_5A_cluster[i] != common_particle_ids[1]) {
            return(new_5A_cluster[i]);
        }
    }
    Error("No uncommon ring particle found.");
    return 0;
}

void Raw_Write_12E(const int *parent11F, int uncommon_sp3_ring_particle) {// now we have found the 12E
    int i;
    int clusSize = 12;

    if (n12E == m12E) {
        hc12E = resize_2D_int(hc12E, m12E, m12E + incrStatic, clusSize, -1);
        m12E = m12E + incrStatic;
    }


    for (i = 0; i < 11; i++) {
        hc12E[n12E][i] = parent11F[i];
    }
    hc12E[n12E][11] = uncommon_sp3_ring_particle;
}

void Cluster_Write_12E() {

    if(s12E[hc12E[n12E][0]] == 'C') s12E[hc12E[n12E][0]] = 'B';
    s12E[hc12E[n12E][1]] = 'O';
    s12E[hc12E[n12E][2]] = 'O';
    s12E[hc12E[n12E][3]] = 'O';
    s12E[hc12E[n12E][4]] = 'O';
    s12E[hc12E[n12E][5]] = 'O';
    s12E[hc12E[n12E][6]] = 'O';
    if(s12E[hc12E[n12E][7]] == 'C') s12E[hc12E[n12E][7]] = 'B';
    if(s12E[hc12E[n12E][8]] == 'C') s12E[hc12E[n12E][8]] = 'B';
    if(s12E[hc12E[n12E][9]] == 'C') s12E[hc12E[n12E][9]] = 'B';
    if(s12E[hc12E[n12E][10]] == 'C') s12E[hc12E[n12E][10]] = 'B';
    if(s12E[hc12E[n12E][11]] == 'C') s12E[hc12E[n12E][11]] = 'B';

    n12E++;
}