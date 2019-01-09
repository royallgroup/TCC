#ifndef TCC_11F_H
#define TCC_11F_H

void Clusters_Get11F_13K();

int are_spindles_bonded(int first_5A_id, int second_5A_id, int *bonded_pairs);

int count_bonded_ring_particles_11F(int cp, const int *first_5A, const int *second_5A, int *bonded_particle_ids);

int get_bonded_6As(int common_particle, int *bonded_particle_ids, int *extra_particles, int *bonded_6A_id, const int *bonded_pairs);

void write_11F(int common_particle, const int *extra_particles, const int *first_5A, const int *second_5A);

#endif
