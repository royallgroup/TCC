#ifndef TCC_CELL_LIST_H
#define TCC_CELL_LIST_H

void Get_Bonds_With_Voronoi_And_Cell_List();

int cell_list_get_particle_1_neighbours(int i, int num_particle_1_neighbours, int *particle_1_neighbours,
                                        int *particle_1_bonds, double *particle_1_bond_lengths, double *store_dr2,
                                        const int *temp_cnb, int *const *temp_bNums);

void cell_list_get_neigbours(int max_allowed_bonds, int *temp_cnb, int **temp_bNums);

void links();  // sorts all the particles into cells, result given by head-of-chain and linked list arrays

int icell(int tix, int tiy, int tiz);

void Setup_Cell_List();

#endif

