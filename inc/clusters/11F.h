#ifndef TCC_11F_H
#define TCC_11F_H

void Clusters_Get11F_12E_13K();

int do_5As_have_distinct_spindles(const int *first_5A, const int *second_5A);

void check_common_particle(int first_5A_cluster_id, int second_5A_pointer, const int *bonded_pairs);

int get_common_particle_between_5As(const int *first_5A, const int *second_5A);

int get_bonded_ring_particles(int cp, const int *first_5A, const int *second_5A, int *bpi, int *bpj);

int get_bonded_6As(int common_particle, int bpi, int bpj, int *ep1, int *ep2, int *bonded_6A_id, const int *bonded_pairs);

int get_first_6A(const int *first_6A, int bpi, int bpj, int common_particle, int *ep1, int *bonded_6A_id, int potential_6A_pointer, const int *bonded_pairs);

int get_second_6A(const int *first_6A, int bpi, int bpj, int common_particle, int *ep2, const int *bonded_pairs);

void write_11F(int common_particle, int ep1, int ep2, const int *first_5A, const int *second_5A);

#endif
