#ifndef TCC_12E_H
#define TCC_12E_H

void Clusters_Get12E();

int are_5A_spindles_common(const int *first_5A_cluster, const int *first_11F_cluster);

int are_5A_ring_particles_common(const int *cluster_5A, const int *cluster_11F);

int get_uncommon_5A_ring_particle(const int *cluster_11F, const int *cluster_5A);

int is_12E_unique(const int *cluster_11F, int uncommon_sp3_ring_particle, int *trial);

void Write_12E(const int *trial);

#endif