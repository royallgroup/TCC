#ifndef TCC_7T_H
#define TCC_7T_H

void Clusters_Get7T();

int check_ring_bonds(const int *new_5A_cluster);

void check_7T_type(int bond_counter, const int *old_6Z, int new_particle);

void add_7T_a(const int *old_6Z, int new_particle);

void add_7T_s(const int *old_6Z, int new_particle);

#endif
