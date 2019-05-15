#ifndef TCC_11F_H
#define TCC_11F_H

void Clusters_Get11F_13K();

void get_unbonded_5A_particles(int *trial_cluster, const int *first_5A_cluster, const int *second_5A_cluster);

int are_spindles_bonded(int first_5A_id, int second_5A_id, int *trial_cluster);

int count_bonded_ring_particles_11F(const int *first_5A, const int *second_5A, int *trial_cluster);

int get_bonded_6As(int *bonded_6A_id, int *trial_cluster);

void setup_6A_rings(const int *trial_cluster, int *first_6A_ring, int *second_6A_ring);

int check_6A(int *trial_cluster, const int *potential_6A_cluster, const int *ring_particles, int which_6A);

void write_11F(const int *trial_cluster);

#endif
