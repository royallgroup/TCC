#include "cell_list.h"
#include "voronoi_bonds.h"
#include "globals.h"
#include "tools.h"
#include "bonds.h"
#include <math.h>


int cell_list_get_particle_1_neighbours(int i, int *particle_1_neighbours, int *particle_1_bonds,
                                        double *particle_1_bond_lengths, double *store_dr2) {
    int neighbour_pointer, neighbour_ID;
    double squared_distance;

    for (neighbour_pointer = 0; neighbour_pointer < num_bonds[i]; ++neighbour_pointer) {
        neighbour_ID = bNums[i][neighbour_pointer];
        squared_distance = Get_Interparticle_Distance(i, neighbour_ID);
        particle_1_neighbours[neighbour_pointer] = neighbour_ID;
        particle_1_bonds[neighbour_pointer] = 1;
        particle_1_bond_lengths[neighbour_pointer] = squared_distance;
        store_dr2[neighbour_pointer] = squared_distance;
    }
    return num_bonds[i];
}

void fill_cell_list() {
    int particle;
    int x_index, y_index, z_index;
    int scalar_index;

    for(particle=0; particle<N; particle++) {
        x_index = (int) floor(x[particle] / cell_len_x);
        y_index = (int) floor(y[particle] / cell_len_y);
        z_index = (int) floor(z[particle] / cell_len_z);

        scalar_index = get_scalar_cell_index(x_index, y_index, z_index);
        linked_list[particle] = head[scalar_index];
        head[scalar_index] = particle;
    }
}

int get_scalar_cell_index(int x_index, int y_index, int z_index) {
    return ((x_index + n_cells_x) % n_cells_x) * n_cells_y * n_cells_z +
           ((y_index + n_cells_y) % n_cells_y) * n_cells_z +
           ((z_index + n_cells_z) % n_cells_z);
}

void Setup_Cell_List() {
    int i;

    n_cells_x = (int)(sidex/rcutAA);
    n_cells_y = (int)(sidey/rcutAA);
    n_cells_z = (int)(sidez/rcutAA);
    if (n_cells_x < 3 || n_cells_y < 3 || n_cells_z < 3) Error_no_free("main(): M<3, too few cells");
    n_cells_total = n_cells_x*n_cells_y*n_cells_z;

    printf("x_cells %d y_cells %d z_cells %d total_cells %d\n", n_cells_x, n_cells_y, n_cells_z, n_cells_total);

    head=malloc((n_cells_total)*sizeof(int));
    linked_list=malloc((current_frame_particle_number)*sizeof(int));

    for (i=0; i<n_cells_total; i++) {
        head[i]=-1;
    }
    for (i=0; i<current_frame_particle_number; i++) {
        linked_list[i]=0;
    }


}
