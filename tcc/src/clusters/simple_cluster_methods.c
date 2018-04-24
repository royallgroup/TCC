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
    //!  Function to count number of common ring particles between two clusters and get thier ids
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

int count_uncommon_ring_particles(const int *cluster_1, const int *cluster_2, int num_in_ring_1, int num_in_ring_2,
                                  int *common_particle_ids) {
    //!  Function to count number of uncommon particles between two clusters and get their ids
    /*!
   *  \param cluster_1 - a pointer to a cluster stored in an hc memory array
   *  \param cluster_2 - a pointer to a cluster stored in an hc memory array
   *  \param num_in_ring_1 - the number of particles in the ring of cluster_1, this is 3, 4 or 5 for the basic clusters
   *  \param num_in_ring_2 - the number of particles in the ring of cluster_1, this is 3, 4 or 5 for the basic clusters
   *  \param uncommon_particle_ids - a pointer to an array of length num_in_ring_1 + num_in_ring_2, ids of uncommon particles will be written to this array
   *  \return an integer giving the number of uncommon ring particles between the clusters
   */
    int num_uncommon_particles = 0;

    for (int first_ring_pointer = 0; first_ring_pointer < num_in_ring_1; ++first_ring_pointer) {
        int is_uncommon_particle = 0;
        for (int second_ring_pointer = 0; second_ring_pointer < num_in_ring_2; ++second_ring_pointer) {
            if (cluster_1[first_ring_pointer] == cluster_2[second_ring_pointer]) {
                is_uncommon_particle = 1;
                break;
            }
        }
        if(is_uncommon_particle == 0) {
            common_particle_ids[num_uncommon_particles] = cluster_1[first_ring_pointer];
            num_uncommon_particles++;
        }
    }

    for (int first_ring_pointer = 0; first_ring_pointer < num_in_ring_1; ++first_ring_pointer) {
        int is_uncommon_particle = 0;
        for (int second_ring_pointer = 0; second_ring_pointer < num_in_ring_2; ++second_ring_pointer) {
            if (cluster_2[first_ring_pointer] == cluster_1[second_ring_pointer]) {
                is_uncommon_particle = 1;
                break;
            }
        }
        if(is_uncommon_particle == 0) {
            common_particle_ids[num_uncommon_particles] = cluster_2[first_ring_pointer];
            num_uncommon_particles++;
        }
    }

    return num_uncommon_particles;

}