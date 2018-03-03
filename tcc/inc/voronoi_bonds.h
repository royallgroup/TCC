#ifndef TCC_VORONOI_BONDS_H
#define TCC_VORONOI_BONDS_H

void Get_Bonds_With_Voronoi();

double is_particle_bonded(int p1, int p2, int p3);

void check_bond_cut_offs(int particle_1, int num_particle_1_neighbours, const int *sorted_particle_1_neighbours,
                         const double *sorted_particle_1_bond_lengths, int *Sb);

void add_new_voronoi_bond(int particle_1, int num_particle_1_neighbours, const int *sorted_particle_1_neighbours,
                          const double *store_dr2, const int *particle_1_bonds);

void Remove_Unbonded_Neighbours(int particle_1, int num_particle_1_neighbours, const int *sorted_particle_1_neighbours,
                                int *Sb);

int get_particle_1_neighbours(int particle_1, int max_allowed_bonds, int *particle_1_bonds,
                              double *particle_1_bond_lengths, double *store_dr2);

void Insertion_Sort_Bond_Lengths(int num_particle_1_neighbours, const int *particle_1_neighbours, int *sorted_particle_1_neighbours,
                                 const double *particle_1_bond_lengths, double *sorted_particle_1_bond_lengths);

#endif
