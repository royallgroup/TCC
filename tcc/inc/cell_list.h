#ifndef TCC_CELL_LIST_H
#define TCC_CELL_LIST_H

void Get_Bonds_With_Voronoi_And_Cell_List();

int cell_list_get_particle_1_neighbours(int i, int *particle_1_neighbours, int *particle_1_bonds,
                                        double *particle_1_bond_lengths, double *store_dr2);

void cell_list_get_neigbours(int max_allowed_bonds, int *temp_cnb, int **temp_bNums);

void fill_cell_list();  // sorts all the particles into cells, result given by head-of-chain and linked list arrays

int get_scalar_cell_index(int tix, int tiy, int tiz);

void Setup_Cell_List();

#endif

