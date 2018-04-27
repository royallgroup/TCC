#ifndef TCC_11C_H
#define TCC_11C_H

void Clusters_Get11C();

int count_particles_bonded_to_common(const int *cluster, const int *common_particles, int *bonded_particles);

int count_bonded_ring_particles_11C(const int *common_ring_particles, const int *first_7A_cluster, const int *second_7A_cluster);

void Cluster_Write_11C(int *trial);

#endif