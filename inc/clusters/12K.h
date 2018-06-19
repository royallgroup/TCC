#ifndef TCC_12K_H
#define TCC_12K_H

void Clusters_Get12K();

void get_12K_ring_bonds(int ptr_11A, int (*sp3_rings)[3]);

void find_12K_cluster(int *parent_11A_cluster, const int *sp3_ring);

void Cluster_Write_12K(int ep, const int *parent_11A_cluster);

#endif