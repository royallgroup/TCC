#ifndef TCC_SIMPLE_CLUSTER_METHODS_H
#define TCC_SIMPLE_CLUSTER_METHODS_H

int count_common_spindles_between_5As(const int *first_5A_cluster, const int *second_5A_cluster, int *scom);

int count_common_ring_particles(const int *cluster_1, const int *cluster_2, int num_particles_in_ring, int* common_particle_ids);

int count_uncommon_ring_particles(const int *cluster_1, const int *cluster_2, int num_in_ring_1, int num_in_ring_2, int *common_particle_ids);

#endif
