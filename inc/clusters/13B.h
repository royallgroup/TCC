#ifndef TCC_13B_H
#define TCC_13B_H

void Clusters_Get13B();

int check_rings_are_bonded(const int *first_7A_cluster, const int *second_7A_cluster);

int count_bonds_to_ring(int particle_id, const int *first_7A_cluster);

int check_rings_are_uncommon(const int *first_7A_cluster, const int *second_7A_cluster);

void Cluster_Write_13B(const int *first_7A_cluster, const int *second_7A_cluster, int common_spindle_id,
                       const int *uncommon_spindle_ids);

#endif
