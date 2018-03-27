#ifndef TCC_11W_H
#define TCC_11W_H

void Clusters_Get11W();

int get_11W_extra_particle(int id_10B, int spindle_10B);

int is_particle_in_10B(int particle_id, int id_10B);

int is_particle_bonded_to_7As(int id_10B, int extra_particle);

void resize_hc11W();

void populate_hc11W(int id_10B, int extra_particle);

void populate_s11W();

#endif //TCC_11W_H
