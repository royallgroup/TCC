#ifndef TCC_12K_H
#define TCC_12K_H

void Clusters_Get12K();

void get_12K_ring_bonds(int ptr_11A, int (*sp3_rings)[3]);

void find_12K_cluster(int ptr_11A, const int *sp3_ring);

int is_particle_in_11A(int ptr_11A, int particle_id);

void Cluster_Write_12K(int ep, int id_11A);

#endif