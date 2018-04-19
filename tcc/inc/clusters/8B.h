#ifndef TCC_8B_H
#define TCC_8B_H

void Clusters_Get8B();

int is_particle_in_7A(int first_7A_id, int new_particle_id);

int count_bonds_to_7A_ring(int first_7A_id, int new_particle_id);

void Cluster_Write_8B(int *first_7A_cluster, int new_particle_id);

#endif
