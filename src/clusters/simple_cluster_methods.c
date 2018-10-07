#include "simple_cluster_methods.h"

/*!  Count number of common ring particles between two clusters and get thier ids
*
*  @param cluster_1 - a pointer to a cluster stored in an hc memory array
*  @param cluster_2 - a pointer to a cluster stored in an hc memory array
*  @param cluster_1_ring_particles - the number of ring particles in cluster 1
*  @param cluster_2_ring_particles - the number of ring particles in cluster 2
*  @param common_particle_ids - a pointer to an array of length num_particles_in_ring, ids of common particles will be written to this array
*  @return an integer giving the number of common ring particles between the clusters
*
*  This function assumes that the n ring particles are the first n particles of the cluster. This is true for the
*  basic clusters but may not be true for other clusters.
*/
int count_common_ring_particles(const int *cluster_1, const int *cluster_2, int cluster_1_ring_particles, int cluster_2_ring_particles,
                                int *common_particle_ids) {

    int num_common_particles = 0;
    for (int first_ring_pointer = 0; first_ring_pointer < cluster_1_ring_particles; ++first_ring_pointer) {
        for (int second_ring_pointer = 0; second_ring_pointer < cluster_2_ring_particles; ++second_ring_pointer) {
            if (cluster_1[first_ring_pointer] == cluster_2[second_ring_pointer]) {
                common_particle_ids[num_common_particles] = cluster_1[first_ring_pointer];
                num_common_particles++;
                break;
            }
        }
    }
    return num_common_particles;

}

//!  Count number of uncommon particles between two clusters and get their ids
/*!
*  @param cluster_1 - a pointer to a cluster stored in an hc memory array
*  @param cluster_2 - a pointer to a cluster stored in an hc memory array
*  @param num_in_ring_1 - the number of particles in the ring of cluster_1, this is 3, 4 or 5 for the basic clusters
*  @param num_in_ring_2 - the number of particles in the ring of cluster_1, this is 3, 4 or 5 for the basic clusters
*  @param uncommon_particle_ids - a pointer to an array of length num_in_ring_1 + num_in_ring_2, ids of uncommon particles will be written to this array
*  @return an integer giving the number of uncommon ring particles between the clusters
*/
int count_uncommon_ring_particles(const int *cluster_1, const int *cluster_2, int num_in_ring_1, int num_in_ring_2,
                                  int *uncommon_particle_ids) {

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
            uncommon_particle_ids[num_uncommon_particles] = cluster_1[first_ring_pointer];
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
            uncommon_particle_ids[num_uncommon_particles] = cluster_2[first_ring_pointer];
            num_uncommon_particles++;
        }
    }

    return num_uncommon_particles;

}

//!  Count number of common spindle particles between two clusters and get thier ids
/*!
*  @param cluster_1 - a pointer to a cluster stored in an hc memory array
*  @param cluster_2 - a pointer to a cluster stored in an hc memory array
*  @param cluster_1_size - the number of particles in cluster_1
*  @param cluster_2_size - the number of particles in cluster_2
*  @param common_spindle_ids - a pointer to an array of length 2 where common spindle ids will be written
*  @return an integer giving the number of common ring particles between the clusters
*
*  This function assumes that the spindle particles are the final two particles in the cluster, this is true for
*  all of the basic clusters but may not be true for larger clusters.
*/
int count_common_spindle_particles(const int *cluster_1, const int *cluster_2, int cluster_1_size, int cluster_2_size, int *common_spindle_ids) {

    int num_common_spindles = 0;

    for (int i = cluster_1_size - 2; i < cluster_1_size; i++) {
        for (int j = cluster_2_size - 2; j < cluster_2_size; j++) {
            if (cluster_1[i] == cluster_2[j]) {
                common_spindle_ids[num_common_spindles] = cluster_1[i];
                num_common_spindles++;
            }
        }
    }
    return num_common_spindles;
}

//!  Get the id of the uncommon spindle in an sp3c, sp4c or sp5c cluster
/*!
*  @param cluster - a pointer to a cluster stored in an hc memory array
*  @param cluster_size - the number of particles in cluster
*  @param common_spindle_id - the id of the common particle
*  @return the id_of the uncommon spindle
*/
int get_uncommon_spindle(const int *cluster, int cluster_size, int common_spindle_id) {


    if (common_spindle_id == cluster[cluster_size - 2]) {
        return cluster[cluster_size - 1];
    }
    else {
        return cluster[cluster_size - 2];
    }
}

//!  Determine if a particle id exists within a cluster
/*!
*  @param cluster - a pointer to a cluster stored in an hc memory array
*  @param cluster_size - the number of particles in cluster
*  @param particle_id - the id of the particle to check
*  @return 1 if particle_id is in cluster, 0 if particle_id is not in cluster
*/
int is_particle_in_cluster(const int *cluster, int cluster_size, int particle_id) {


    for (int i = 0; i < cluster_size; i++) {
        if (cluster[i] == particle_id) {
            return 1;
        }
    }
    return 0;
}

//!  Determine if there are any shared particles between two clusters
/*!
*  @param cluster_1 - a pointer to a cluster stored in an hc memory array
*  @param cluster_2 - a pointer to a cluster stored in an hc memory array
*  @param cluster_1_size - the number of particles in cluster_1
*  @param cluster_2_size - the number of particles in cluster_2
*  @return 1 if there are no common particles, 0 if there is at least 1 common particle
*/

int are_clusters_distinct(const int *cluster_1, const int *cluster_2, int cluster_1_size, int cluster_2_size) {

    for (int i = 0; i < cluster_1_size; i++) {
        if (is_particle_in_cluster(cluster_2, cluster_2_size, cluster_1[i])) {
            return 0;
        }
    }
    return 1;
}

//!  Determine if a trial cluster already exists in list of clusters
/*!
*  @param trial - a pointer to a cluster, particle ids sorted in ascending order
*  @param cluster_size - the number of particles in the trial cluster
*  @param cluster_list - a pointer to a list of clusters with the particle ids making up each cluster sorted in ascending order
*  @param num_clusters - the number of clusters in cluster list
*  @return 0 if the cluster is unique, 1 if the cluster already exists in the list
*/
int check_unique_cluster(const int *trial, const int cluster_size, int **cluster_list, const int num_clusters) {


    int i;
    for (int existing_cluster_pointer = 0; existing_cluster_pointer < num_clusters; ++existing_cluster_pointer) {
        for (i = 0; i < cluster_size; ++i) {
            if (trial[i] != cluster_list[existing_cluster_pointer][i]) {
                break;
            }
        }
        if (i == cluster_size) {
            return 1;
        }
    }
    return 0;
}
