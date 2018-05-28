#ifndef TCC_11A_H
#define TCC_11A_H

void Clusters_Get11A();

int Check_6A_rings_bonded(const int *first_6A_cluster, const int *second_6A_cluster);

void Cluster_Write_11A(const int *first_6A, const int *second_6A, const int *sother, int scom);

#endif
