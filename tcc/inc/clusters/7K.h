#ifndef TCC_7K_H
#define TCC_7K_H

void Clusters_Get7K();

int count_common_spindles_between_5As(const int *first_5A_cluster, const int *second_5A_cluster, int *scom);

void get_other_spindle_ids(const int *first_5A_cluster, const int *second_5A_cluster, int scom, int *sother);

int is_particle_in_5A(const int *five_A_cluster, int particle_id);

int count_common_ring_particles(const int *first_5A_cluster, const int *second_5A_cluster, int *sp3_com);

void Cluster_Write_7K();

#endif
