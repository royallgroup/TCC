#ifndef TCC_12E_H
#define TCC_12E_H

void Clusters_Get12E();

int get_uncommon_5A_ring_particle(const int *common_particle_ids, const int *first_5A_cluster);

void Write_12E(const int *first_11F, int uncommon_sp3_ring_particle);

#endif