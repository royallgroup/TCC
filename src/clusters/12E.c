#include <globals.h>
#include <tools.h>
#include "simple_cluster_methods.h"
#include "12E.h"

//!  An 12E cluster is the intersection of an 11F and 5A cluster.
/*!
*  Find 12E clusters
*  A 12E cluster is an 11F and 5A cluster where:
*      - The spindle atoms of the 5A cluster are common with the uncommon spindle atoms
*        of the 6A clusters constituting the 11F cluster.
*      - Of the SP3 ring particles in the 5A cluster, two are common with rd1 and rd2
*        from the 11F cluster, and one is new.
*
*  Cluster output: BBBBBBBBBBBB
*  Storage order: particles ordered by id x 12
*/
void Clusters_Get12E() {
    int common_particle_ids[4], trial[12];

    for (int id_11F = 0; id_11F < n11F; id_11F++) {
        int *cluster_11F = hc11F[id_11F];

        for (int pointer_5A = 0; pointer_5A < nmem_sp3c[cluster_11F[1]]; pointer_5A++) {
            int id_5A = mem_sp3c[cluster_11F[1]][pointer_5A];
            int *cluster_5A = hcsp3c[id_5A];

            // Check spindles of new 5A are common with uncommon spindles of 6As in 11F
            if (are_5A_spindles_common(cluster_5A, cluster_11F) == 0) continue;

            //if (are_5A_ring_particles_common(first_5A_cluster, first_11F_cluster) == 0) continue;
            // Check two new 5A ring particles are common with 11F ring particles
            int num_common_particles = 0;
            for (int ring_pointer = 0;
                // loop through 5A ring particles
                ring_pointer < 3; ring_pointer++) {
                for (int m = 7; m < 11; m++) { // loop through ring particles of 11F
                    if (cluster_5A[ring_pointer] == cluster_11F[m]) {
                        common_particle_ids[num_common_particles] = cluster_5A[ring_pointer];
                        num_common_particles++;
                    }
                }
            }
            if (num_common_particles != 2) continue;

            // Find which of the 3 ring particles in new_5A is the uncommon one
            int uncommon_sp3_ring_particle = get_uncommon_5A_ring_particle(common_particle_ids, cluster_5A);
            // check that new particle is not already in 11F
            if (is_particle_in_cluster(cluster_11F, 11, uncommon_sp3_ring_particle) == 0) {

                // Now need to search the current cluster list to make sure the detected cluster is uniques
                for (int i = 0; i < 11; ++i) {
                    trial[i] = cluster_11F[i];
                }
                trial[11] = uncommon_sp3_ring_particle;

                quickSort(&trial[0], 12);

                if (check_unique_cluster(trial, 12, hc12E, n12E) == 0) {
                    Write_12E(trial);
                }
            }
        }
    }
}

int are_5A_ring_particles_common(const int *first_5A_cluster, const int *first_11F_cluster) {
    return 0;
}

int are_5A_spindles_common(const int *first_5A_cluster, const int *first_11F_cluster) {
    // Check whether the spindles of the 5A are common with the uncommon spindle atoms
    // of the 6A clusters constituting in the 11F cluster
    if (first_5A_cluster[3] == first_11F_cluster[1] || first_5A_cluster[3] == first_11F_cluster[2]) {
        if (first_5A_cluster[4] == first_11F_cluster[1] || first_5A_cluster[4] == first_11F_cluster[2]) {
            return 1;
        }
    }
    return 0;
}

int get_uncommon_5A_ring_particle(const int *common_particle_ids, const int *first_5A_cluster) {

    for (int i = 0; i < 3; i++) {
        if (first_5A_cluster[i] != common_particle_ids[0] && first_5A_cluster[i] != common_particle_ids[1]) {
            return(first_5A_cluster[i]);
        }
    }
    Error("No uncommon ring particle found.");
    return 0;
}

void Write_12E(const int *trial) {
    int clusSize = 12;

    if (n12E == m12E) {
        hc12E = resize_2D_int(hc12E, m12E, m12E + incrStatic, clusSize, -1);
        m12E = m12E + incrStatic;
    }

    for (int i = 0; i < 12; i++) {
        hc12E[n12E][i] = trial[i];
        s12E[hc12E[n12E][i]] = 'B';
    }

    n12E++;
}