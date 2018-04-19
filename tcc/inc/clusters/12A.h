#ifndef TCC_12A_H
#define TCC_12A_H

void Clusters_Get12A();

int get_12A_extra_particle(int id_11C);

int is_particle_in_11C(int particle_id, int id_11C);

int bond_check_12A_extra_particle(int id_11C, int ep);

void resize_hc12A();

void populate_hc12A(int id_11C, int ep);

void populate_s12A();

#endif
