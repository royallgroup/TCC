#ifndef TCC_12E_H
#define TCC_12E_H

void Clusters_Get12E();

int get_uncommon_5A_ring_particle(const int *common, const int *new_5A_cluster);

void Raw_Write_12E(const int* parent11F, int uncommon_sp3_ring_particle);

void Cluster_Write_12E();

#endif