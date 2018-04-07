#ifndef TCC_6Z_H
#define TCC_6Z_H

void Clusters_Get6Z();

int check_for_common_spindle_particles(const int *first_5A_cluster, const int *second_5A_cluster);

int count_spindles_in_ring(const int *first_5A_cluster, const int *second_5A_cluster, int *spindles);

int get_bonds_between_rings(const int *first_5A_cluster, const int *second_5A_cluster, int *common_ring_particles);

void Cluster_Write_6Z(const int *first_spindles, const int *second_spindles, const int *common_ring_particles);

#endif
