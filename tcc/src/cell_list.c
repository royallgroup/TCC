#include "cell_list.h"
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
        squared_distance = squared_bondlengths[i][neighbour_pointer];
        particle_1_neighbours[neighbour_pointer] = neighbour_ID;
        particle_1_bonds[neighbour_pointer] = 1;
        particle_1_bond_lengths[neighbour_pointer] = squared_distance;
        store_dr2[neighbour_ID] = squared_distance;
    }
    return num_bonds[i];
}

void get_all_particle_neighbours() {
    int cell_x, cell_y, cell_z;
    int current_cell_index;

    for(cell_x = 0; cell_x<n_cells_x; cell_x++) {
        for (cell_y = 0; cell_y<n_cells_y; cell_y++) {
            for (cell_z = 0; cell_z < n_cells_z; cell_z++) {
                current_cell_index = get_scalar_cell_index(cell_x, cell_y, cell_z);
                loop_over_neighbouring_cells(cell_x, cell_y, cell_z, current_cell_index);
            }
        }
    }
}

void loop_over_neighbouring_cells(int cell_x, int cell_y, int cell_z, int current_cell_index) {

    int neighbour_cell_index;
    int neighbour_x, neighbour_y, neighbour_z;

    for(neighbour_x = cell_x - 1; neighbour_x < cell_x + 2; neighbour_x++) {
        for (neighbour_y = cell_y - 1; neighbour_y < cell_y + 2; neighbour_y++) {
            for (neighbour_z = cell_z - 1; neighbour_z < cell_z + 2; neighbour_z++) {
                neighbour_cell_index = get_scalar_cell_index(neighbour_x, neighbour_y, neighbour_z);
                loop_over_particles_in_cell(current_cell_index, neighbour_cell_index);
            }
        }
    }
}

void loop_over_particles_in_cell(int current_cell_index, int neighbour_cell_index) {
    int particle_1, particle_2;
    double sq_dist;
    
    particle_1 = head[current_cell_index];
    while (particle_1 != -1) {
        particle_2 = head[neighbour_cell_index];
        while (particle_2 != -1) {
            if (particle_1 < particle_2) {
                sq_dist = Get_Interparticle_Distance(particle_1, particle_2);
                Check_For_Valid_Bond(particle_1, particle_2, sq_dist);
            }
            particle_2 = linked_list[particle_2];
        }
        particle_1 = linked_list[particle_1];
    }
}

void fill_cell_list() {
    int particle;
    int x_index, y_index, z_index;
    int scalar_index;

    for(particle=0; particle<particles_in_current_frame; particle++) {
        x_index = (int) floor(x[particle] / cell_len_x);
        y_index = (int) floor(y[particle] / cell_len_y);
        z_index = (int) floor(z[particle] / cell_len_z);

        scalar_index = get_scalar_cell_index(x_index, y_index, z_index);
        if(scalar_index < 0 || scalar_index >= n_cells_total) {
            Error("Scalar index out of range");
        }
        linked_list[particle] = head[scalar_index];
        head[scalar_index] = particle;
    }
}

int get_scalar_cell_index(int x_index, int y_index, int z_index) {
    return ((x_index + n_cells_x) % n_cells_x) * n_cells_y * n_cells_z +
           ((y_index + n_cells_y) % n_cells_y) * n_cells_z +
           ((z_index + n_cells_z) % n_cells_z);
}

void set_up_cell_list() {
    int i;

    n_cells_x = (int)(sidex/rcutAA);
    n_cells_y = (int)(sidey/rcutAA);
    n_cells_z = (int)(sidez/rcutAA);
    if (n_cells_x < 3 || n_cells_y < 3 || n_cells_z < 3) Error_no_free("main(): M<3, too few cells");
    n_cells_total = n_cells_x*n_cells_y*n_cells_z;
    cell_len_x = sidex/n_cells_x;
    cell_len_y = sidey/n_cells_y;
    cell_len_z = sidez/n_cells_z;

    printf("x_cells %d y_cells %d z_cells %d total_cells %d\n", n_cells_x, n_cells_y, n_cells_z, n_cells_total);

    head=malloc((n_cells_total)*sizeof(int));
    linked_list=malloc((particles_in_current_frame)*sizeof(int));

    for (i=0; i<n_cells_total; i++) {
        head[i]=-1;
    }
    for (i=0; i<particles_in_current_frame; i++) {
        linked_list[i]=0;
    }
}
