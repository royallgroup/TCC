#include <clusters/simple_cluster_methods.h>
#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "11F.h"
#include "13K.h"

//!  An 11F cluster is the intersection of two 5A and two 6A clusters
/*!
*  Find 11F clusters
*  An 11F is constructed from two 5A and two 6A clusters where:
*      - Each spindle of 5Ai is bonded to a spindle of 5Aj.
*      - There is one common ring particle between the 5A clusters.
*      - There is one bonded pair of ring particle between the 5A clusters.
*      - 6A_i has a distinct spindle and a spindle which is the common ring particle of the 5As
 *     - 6A_j has a distinct spindle and a spindle which is the common ring particle of the 5As
*      - Each 6A has two bonded 5A spindles and two bonded 5A ring particles as its ring.
*
*  Cluster output: BOOOOOOBBBB
*  Storage order: 5A_common particle, 6A_uncommon_spindle x 2, 5A_spindles x 4, 5A_ring_particles x 4
*
*/
void Clusters_Get11F_13K() {

    int bonded_6A_id;              // The id of the first 6A cluster to pass to the 13A method
    int trial_cluster[11];         // The ids of the particles as we build the cluster
    // The trial cluster is built up with particles to form the 11F
    // Trial[0] is the common ring particle of the 5As, this is also the common spindle of the 6As
    // Trial[1 and 2] are the uncommon 6A spindles. These are the only particles not in the 5As
    // Trial[3 and 4] are the spindles of the first 5A cluster.
    // Trial[5 and 6] are the spindles of the second 5A cluster.
    // Trial[3] is bonded to trial[5] and trial[4] is bonded to trial[6].
    // Trial[7 and 8] are the distinct but bonded ring particles of the 2 5As. These are the common ring particles of the two 6As
    // Trial[9 and 10] are the uncommon unbonded ring particles of the 5As. These are not in the 6As.

    for(int first_5A_id = 0; first_5A_id < nsp3c; first_5A_id++) {
        int *first_5A_cluster = hcsp3c[first_5A_id];
        // loop over only the rings of the 5A clusters
        for (int first_5A_ring_pointer = 0; first_5A_ring_pointer < 3; first_5A_ring_pointer++) {
            int first_5A_ring_particle = first_5A_cluster[first_5A_ring_pointer];
            for (int second_5A_pointer = 0; second_5A_pointer < nmem_sp3c[first_5A_ring_particle]; second_5A_pointer++) {
                int second_5A_id = mem_sp3c[first_5A_ring_particle][second_5A_pointer];
                int *second_5A_cluster = hcsp3c[second_5A_id];
                // Don't detect clusters twice
                if (second_5A_id <= first_5A_id) continue;

                // Check that the 5A spindles are bonded
                if (are_spindles_bonded(first_5A_id, second_5A_id, trial_cluster) != 1) continue;

                // There must be only 1 particle common between the two 5As
                int return_array[6];
                if (count_common_particles(first_5A_cluster, second_5A_cluster, 3, 3, return_array) != 1) continue;
                trial_cluster[0] = return_array[0];

                // There must be exactly one pair of bonded ring particles between the 5As
                if (count_bonded_ring_particles_11F(first_5A_cluster, second_5A_cluster, trial_cluster) != 1) continue;

                if (get_bonded_6As(&bonded_6A_id, trial_cluster)) {
                    get_unbonded_5A_particles(trial_cluster, first_5A_cluster, second_5A_cluster);

                    write_11F(trial_cluster);

                    if (do13K == 1) {
                        if (Clusters_Get13K(first_5A_id, second_5A_id, bonded_6A_id)) {
                            ++n13K;
                        }
                    }
                    ++n11F;
                }
            }
        }
    }
}

void get_unbonded_5A_particles(int *trial_cluster, const int* first_5A_cluster, const int* second_5A_cluster) {
    // We know all particles except for the unbonded uncommon 5A ring particles.
    // Finding these particles completes the cluster.

    int known_ring_particles[3];
    known_ring_particles[0] = trial_cluster[0];
    known_ring_particles[1] = trial_cluster[7];
    known_ring_particles[2] = trial_cluster[8];

    for (int i = 0; i < 3; ++i) {
        if(is_particle_in_cluster(known_ring_particles, 3, first_5A_cluster[i]) == 0)
            trial_cluster[9] = first_5A_cluster[i];
        if(is_particle_in_cluster(known_ring_particles, 3, second_5A_cluster[i]) == 0)
            trial_cluster[10] = second_5A_cluster[i];
    }
}

int are_spindles_bonded(int first_5A_id, int second_5A_id, int *trial_cluster) {

    int *first_5A_cluster = hcsp3c[first_5A_id];
    int *second_5A_cluster = hcsp3c[second_5A_id];

    // First check that the spindles are unique
    if (first_5A_cluster[3] == second_5A_cluster[3]) return 0;
    if (first_5A_cluster[3] == second_5A_cluster[4]) return 0;
    if (first_5A_cluster[4] == second_5A_cluster[4]) return 0;

    // Now check if both pairs of spindles are bonded
    if (Bonds_BondCheck(first_5A_cluster[3], second_5A_cluster[3]) == 1 &&
        Bonds_BondCheck(first_5A_cluster[4], second_5A_cluster[4]) == 1) {
        trial_cluster[3] = first_5A_cluster[3];
        trial_cluster[4] = first_5A_cluster[4];
        trial_cluster[5] = second_5A_cluster[3];
        trial_cluster[6] = second_5A_cluster[4];
        return 1;
    }
    else if (Bonds_BondCheck(first_5A_cluster[3], second_5A_cluster[4]) == 1 &&
               Bonds_BondCheck(first_5A_cluster[4], second_5A_cluster[3]) == 1) {
        trial_cluster[3] = first_5A_cluster[3];
        trial_cluster[4] = first_5A_cluster[4];
        trial_cluster[5] = second_5A_cluster[4];
        trial_cluster[6] = second_5A_cluster[3];
        return 1;
    }
    return 0;
}

int count_bonded_ring_particles_11F(const int *first_5A, const int *second_5A, int *trial_cluster) {
    int num_bonded_pairs = 0;
    int common_particle = trial_cluster[0]; // The particle common to all clusters

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (first_5A[i] != common_particle && second_5A[j] != common_particle) {
                if (Bonds_BondCheck(first_5A[i], second_5A[j]) == 1) {
                    trial_cluster[7] = first_5A[i];
                    trial_cluster[8] = second_5A[j];
                    num_bonded_pairs++;
                }
            }
        }
    }
    return num_bonded_pairs;
}

int get_bonded_6As(int *bonded_6A_id, int *trial_cluster) {
    // Returns 1 if both 6A clusters are successfully found, 0 if not.

    int first_6A_detected = 0;
    int second_6A_detected = 0;
    int common_particle = trial_cluster[0];     // Particle which is common to all clusters

    // We know the identity of the ring particles already since they are in the 5As
    // so to find the 6As we just have to search for a 6A containing these ring particles
    int first_6A_ring[4];
    int second_6A_ring[4];
    setup_6A_rings(trial_cluster, first_6A_ring, second_6A_ring);

    // Loop through all 6A clusters
    // We know that the 6As must have the 5A common particle as a spindle so we can use the mem arrays
    for (int mem_pointer = 0; mem_pointer < nmem_sp4c[common_particle]; ++mem_pointer) {
        int potential_6A_id = mem_sp4c[common_particle][mem_pointer];
        int *potential_6A_cluster = hcsp4c[potential_6A_id];

        if (first_6A_detected == 0) {
            if (check_6A(trial_cluster, potential_6A_cluster, first_6A_ring, 1) == 1) {
                first_6A_detected = 1;
                *bonded_6A_id = potential_6A_id;
            }
        }
        if (second_6A_detected == 0) {
            if (check_6A(trial_cluster, potential_6A_cluster, second_6A_ring, 2) == 1) {
                second_6A_detected = 1;
            }
        }
        if (first_6A_detected == 1 && second_6A_detected == 1) {
            return 1;
        }
    }
    return 0;
}

int check_6A(int *trial_cluster, const int *potential_6A_cluster, const int *ring_particles, int which_6A) {
    int common_particle = trial_cluster[0];
    int common_list[4];

    if (count_common_particles(ring_particles, potential_6A_cluster, 4, 4, common_list) == 4) {
        if (potential_6A_cluster[4] == common_particle) {
            trial_cluster[which_6A] = potential_6A_cluster[5];
        }
        else {
            trial_cluster[which_6A] = potential_6A_cluster[4];
        }
        return 1;
    }
    return 0;
}

void setup_6A_rings(const int *trial_cluster, int *first_6A_ring, int *second_6A_ring) {
    first_6A_ring[0] = trial_cluster[3];    // bonded spindles of 5As
    first_6A_ring[1] = trial_cluster[5];    // bonded spindles of 5As
    first_6A_ring[2] = trial_cluster[7];    // bonded ring particles of 5A
    first_6A_ring[3] = trial_cluster[8];    // bonded ring particles of 5A

    second_6A_ring[0] = trial_cluster[4];    // bonded spindles of 5As
    second_6A_ring[1] = trial_cluster[6];    // bonded spindles of 5As
    second_6A_ring[2] = trial_cluster[7];    // bonded ring particles of 5A
    second_6A_ring[3] = trial_cluster[8];    // bonded ring particles of 5A
}

void write_11F(const int *trial_cluster) {
    int clusSize = 11;

    if (n11F == m11F) {
        hc11F = resize_2D_int(hc11F, m11F, m11F + incrStatic, clusSize, -1);
        m11F = m11F + incrStatic;
    }

    for (int i = 0; i < 11; ++i) {
        hc11F[n11F][i] = trial_cluster[i];
    }

    if (s11F[hc11F[n11F][0]] == 'C') s11F[hc11F[n11F][0]] = 'B';
    for(int i = 1; i < 7; i++) {
        s11F[hc11F[n11F][i]] = 'O';
    }
    for(int i = 7; i < 11; i++) {
        if (s11F[hc11F[n11F][i]] == 'C') s11F[hc11F[n11F][i]] = 'B';
    }
}