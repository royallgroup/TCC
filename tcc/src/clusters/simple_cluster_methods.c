#include "simple_cluster_methods.h"

int count_common_spindles_between_5As(const int *first_5A_cluster, const int *second_5A_cluster, int *scom) {
    int num_common_spindles = 0;
    for (int ring_1_pointer = 3; ring_1_pointer < 5; ring_1_pointer++) {
        for (int ring_2_pointer = 3; ring_2_pointer < 5; ring_2_pointer++) {
            if (first_5A_cluster[ring_1_pointer] == second_5A_cluster[ring_2_pointer]) {
                *scom = first_5A_cluster[ring_1_pointer];
                num_common_spindles++;
            }
        }
    }
    return num_common_spindles;
}

int count_common_ring_particles(const int *cluster_1, const int *cluster_2, int num_particles_in_ring, int* common_particle_ids) {
    //!  Function to count number of common ring particles between two clusters
    /*!
   *  \param cluster_1 - a pointer to a cluster stored in an hc memory array
   *  \param cluster_2 - a pointer to a cluster stored in an hc memory array
   *  \param num_particles_in_ring - the number of particles in the ring of the cluster, this is 3, 4 or 5 for the basic clusters
   *  \param common_particle_ids - a pointer to an array of length num_particles_in_ring, ids of common particles will be written to this array
   *  \return an integer giving the number of common ring particles between the clusters
   */

    int num_common_particles = 0;
    for (int first_ring_pointer = 0; first_ring_pointer < num_particles_in_ring; ++first_ring_pointer) {
        for (int second_ring_pointer = 0; second_ring_pointer < num_particles_in_ring; ++second_ring_pointer) {
            if (cluster_1[first_ring_pointer] == cluster_2[second_ring_pointer]) {
                common_particle_ids[num_common_particles] = cluster_1[first_ring_pointer];
                num_common_particles++;
                break;
            }
        }
    }
    return num_common_particles;

}