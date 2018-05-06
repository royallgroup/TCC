#ifndef TCC_HCP_H
#define TCC_HCP_H

void Clusters_GetHCP();

int are_spindles_bonded(const int *first_cluster, const int *second_cluster);

void Cluster_Write_HCP(int i, int j, int j2, int k);

#endif
