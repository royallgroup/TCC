#ifndef TCC_11W_H
#define TCC_11W_H

void Clusters_Get11W();

int get_11W_extra_particle(int *parent_10B_cluster, int spindle_10B);

int is_particle_bonded_to_7As(int id_10B, int extra_particle);

void Write_11W(int id_10B, int extra_particle);

#endif