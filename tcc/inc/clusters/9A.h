#ifndef TCC_9A_H
#define TCC_9A_H

void Clusters_Get9A();

int check_spindles_are_uncommon_and_unbonded(int *cluster_1, int *cluster_2);

int count_bonded_ring_particles(const int *first_sp4b_cluster, const int *second_sp4b_cluster, const int *db, int *ob);

void Cluster_Write_9A(int *first_sp4b_cluster, int *second_sp4b_cluster, int *third_sp4b_cluster,
                      int *i_j_common_ring_particles, int *i_j_uncommon_ring_particles);

#endif
