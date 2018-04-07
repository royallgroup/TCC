#ifndef TCC_CELL_LIST_H
#define TCC_CELL_LIST_H

int cell_list_get_particle_1_neighbours(int i, int *particle_1_neighbours, int *particle_1_bonds,
                                        double *particle_1_bond_lengths, double *store_dr2);

void get_all_particle_neighbours();

void loop_over_neighbouring_cells(int cell_x, int cell_y, int cell_z, int current_cell_index);

void loop_over_particles_in_cell(int current_cell_index, int neighbour_cell_index);

void fill_cell_list();  // sorts all the particles into cells, result given by head-of-chain and linked list arrays

int get_scalar_cell_index(int tix, int tiy, int tiz);

void set_up_cell_list();

void free_cell_list();

#endif

