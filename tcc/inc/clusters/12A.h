#ifndef TCC_12A_H
#define TCC_12A_H

void Clusters_Get12A();

int get_12A_extra_particle(int *parent_11C_cluster);

int bond_check_12A_extra_particle(int *first_11C_cluster, int ep);

void Write_12A(const int *first_11C_cluster, int ep);

#endif
