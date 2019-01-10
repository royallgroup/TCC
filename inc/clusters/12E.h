#ifndef TCC_12E_H
#define TCC_12E_H

void Clusters_Get12E();

int is_11F_unique(int *trial, const int *first_11F_cluster, int uncommon_sp3_ring_particle);

int are_5A_ring_particles_common(const int *first_5A_cluster, const int *first_11F_cluster);

int are_5A_spindles_common(const int *first_5A_cluster, const int *first_11F_cluster);

int get_uncommon_5A_ring_particle(const int *common_particle_ids, const int *first_5A_cluster);

void Write_12E(const int *trial);

#endif