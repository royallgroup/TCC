#ifndef TCC_7K_H
#define TCC_7K_H

void Clusters_Get7K();

void get_other_spindle_ids(const int *first_5A_cluster, const int *second_5A_cluster, int common_spindle_id, int *other_spindle_ids);

int is_particle_in_5A(const int *five_A_cluster, int particle_id);

int count_common_ring_particles_7K(const int *first_5A_cluster, const int *second_5A_cluster, int *common_ring_ids);

int get_uncommon_ring_particle(const int *first_5A_cluster, const int *common_ring_ids);

void Cluster_Write_7K(int scom, int *sother, int *sp3_com, int *uncommon_ring_particles);

#endif
