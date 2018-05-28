#include "globals.h"
#include "tools.h"
#include "bonds.h"
#include "voronoi_bonds.h"
#include "cell_list.h"

void Get_Bonds_With_Voronoi() {
    int i, particle_1;
    const int max_allowed_bonds = 4 * max_num_bonds;

    int num_particle_1_neighbours;
    int *particle_1_neighbours;
    double *particle_1_bond_lengths;

    double *store_dr2;
    int *particle_1_bonds;

    particle_1_neighbours = malloc(max_allowed_bonds*sizeof(int));
    particle_1_bonds = malloc(max_allowed_bonds*sizeof(int));
    particle_1_bond_lengths = malloc(max_allowed_bonds*sizeof(double));
    store_dr2 = malloc(particles_in_current_frame*sizeof(double));

    printf("Calculating Voronoi bond network for %ld particles.\n",particles_in_current_frame);

    for (i=0; i<particles_in_current_frame; i++) {
        num_bonds[i] = 0;
    }

    if (use_cell_list == 1) {
        fill_cell_list();
        get_all_particle_neighbours();
    }

    for (particle_1 = 0; particle_1 < particles_in_current_frame; particle_1++) {
        if(use_cell_list == 1) {

            num_particle_1_neighbours = cell_list_get_particle_1_neighbours(particle_1, particle_1_neighbours,
                                                                            particle_1_bonds, particle_1_bond_lengths,
                                                                            store_dr2);
        }
        else{
            num_particle_1_neighbours = get_particle_1_neighbours(particle_1, max_allowed_bonds, particle_1_neighbours,
                                                                  particle_1_bond_lengths, store_dr2);
        }
        Insertion_Sort_Bond_Lengths(num_particle_1_neighbours, particle_1_neighbours, particle_1_bond_lengths);

        for (i=0; i < num_particle_1_neighbours; i++) {
            particle_1_bonds[i] = 1;
        }

        Remove_Unbonded_Neighbours(particle_1, num_particle_1_neighbours, particle_1_neighbours, particle_1_bonds);

        check_bond_cut_offs(particle_1, num_particle_1_neighbours, particle_1_neighbours, particle_1_bond_lengths, particle_1_bonds);

        add_new_voronoi_bond(particle_1, num_particle_1_neighbours, particle_1_neighbours, store_dr2, particle_1_bonds);
    }

    free(particle_1_neighbours);
    free(particle_1_bonds);
    free(particle_1_bond_lengths);
    free(store_dr2);
}

void add_new_voronoi_bond(int particle_1, int num_particle_1_neighbours, const int *sorted_particle_1_neighbours,
                          const double *store_dr2, const int *particle_1_bonds) {
    int particle_2_pointer, particle_2;
    // Overwrite the simple distance bonds
    num_bonds[particle_1] = 0;

    for (particle_2_pointer = 0; particle_2_pointer < num_particle_1_neighbours; ++particle_2_pointer) {
        if (particle_1_bonds[particle_2_pointer] == 1) {
            particle_2 = sorted_particle_1_neighbours[particle_2_pointer];
            if (num_bonds[particle_1] < max_num_bonds) {
                Add_New_Bond(particle_1, particle_2, store_dr2[particle_2]);
            }
            else {
                too_many_bonds(particle_1, particle_2, __func__);
            }
        }
    }
}

void check_bond_cut_offs(int particle_1, int num_particle_1_neighbours, const int *sorted_particle_1_neighbours,
                         const double *sorted_particle_1_bond_lengths, int *particle_1_bonds) {
    int i, particle_2;
    for (i = 0; i < num_particle_1_neighbours; ++i) {
        particle_2 = sorted_particle_1_neighbours[i];
        if (particle_type[particle_1] == 2 && particle_type[particle_2] == 2) {
            if (sorted_particle_1_bond_lengths[i] > rcutBB2) {
                particle_1_bonds[i] = 0;
            }
        } else if (particle_type[particle_1] == 2 || particle_type[particle_2] == 2) {
            if (sorted_particle_1_bond_lengths[i] > rcutAB2) {
                particle_1_bonds[i] = 0;
            }
        }
    }
}

void Remove_Unbonded_Neighbours(int particle_1, const int num_particle_1_neighbours, const int *sorted_particle_1_neighbours, int *particle_1_bonds) {
    int particle_2, particle_3, p3_pointer, p2_pointer;

    for (p3_pointer = 0; p3_pointer < num_particle_1_neighbours - 1; ++p3_pointer) {
        particle_3 = sorted_particle_1_neighbours[p3_pointer];
        for (p2_pointer = p3_pointer + 1; p2_pointer < num_particle_1_neighbours; ++p2_pointer) {
            particle_2 = sorted_particle_1_neighbours[p2_pointer];
            if(is_particle_bonded(particle_1, particle_2, particle_3) == 0) {
                particle_1_bonds[p2_pointer] = 0;
            }
        }
    }
}

double is_particle_bonded(int p1, int p2, int p3) {
    double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
    double x1, x2;

    rijx = x[p1] - x[p2];
    rijy = y[p1] - y[p2];
    rijz = z[p1] - z[p2];
    rikx = x[p1] - x[p3];
    riky = y[p1] - y[p3];
    rikz = z[p1] - z[p3];
    rjkx = x[p2] - x[p3];
    rjky = y[p2] - y[p3];
    rjkz = z[p2] - z[p3];

    if (PBCs == 1) {
        enforce_PBCs(&rijx, &rijy, &rijz);
        enforce_PBCs(&rikx, &riky, &rikz);
        enforce_PBCs(&rjkx, &rjky, &rjkz);
    }

    x1 = rijx * rikx + rijy * riky + rijz * rikz;
    x1 -= rijx * rjkx + rijy * rjky + rijz * rjkz;
    x2 = rikx * rikx + riky * riky + rikz * rikz;
    x2 += rjkx * rjkx + rjky * rjky + rjkz * rjkz;
    x1 = x1 / x2;
    if (x1 - fc > EPS) {
        return 0;
    }
    else {
        return 1;
    }
}

int get_particle_1_neighbours(int particle_1, const int max_allowed_bonds, int *particle_1_bonds,
                                double *particle_1_bond_lengths, double *store_dr2) {
    int particle_2, num_particle_1_neighbours;
    double squared_distance;

    num_particle_1_neighbours = 0;
    for (particle_2 = 0; particle_2 < particles_in_current_frame; ++particle_2) {
        if (particle_1 != particle_2) {
            squared_distance = Get_Interparticle_Distance(particle_1, particle_2);
            if (squared_distance < rcutAA2) {
                if (num_particle_1_neighbours < max_allowed_bonds) {
                    particle_1_bonds[num_particle_1_neighbours] = particle_2;
                    particle_1_bond_lengths[num_particle_1_neighbours] = squared_distance;
                    store_dr2[particle_2] = squared_distance;
                    num_particle_1_neighbours++;
                }
                else {
                    too_many_bonds(particle_1, particle_2, __func__);
                }
            }
        }
    }
    return num_particle_1_neighbours;
}

void Insertion_Sort_Bond_Lengths(int num_particle_1_neighbours, int *particle_1_neighbours,
                                 double *particle_1_bond_lengths) {
    char error_message[200];
    int i, k, l;
    int cnbs2;

    int *sorted_particle_1_neighbours;
    double  *sorted_particle_1_bond_lengths;
    sorted_particle_1_neighbours = malloc(num_particle_1_neighbours*sizeof(int));
    sorted_particle_1_bond_lengths = malloc(num_particle_1_neighbours*sizeof(double));

    cnbs2 = 0;
    for (i = 0; i < num_particle_1_neighbours; ++i) {
        for (k = 0; k < cnbs2; ++k) { // find spot to insert particle_1_neighbours[i]
            if (particle_1_bond_lengths[i] < sorted_particle_1_bond_lengths[k]) {
                for (l = cnbs2; l > k; --l) {
                    sorted_particle_1_neighbours[l] = sorted_particle_1_neighbours[l - 1];
                    sorted_particle_1_bond_lengths[l] = sorted_particle_1_bond_lengths[l - 1];
                }
                sorted_particle_1_neighbours[k] = particle_1_neighbours[i];
                sorted_particle_1_bond_lengths[k] = particle_1_bond_lengths[i];
                break;
            }
        }
        if (k == cnbs2) {
            sorted_particle_1_neighbours[cnbs2] = particle_1_neighbours[i];
            sorted_particle_1_bond_lengths[cnbs2] = particle_1_bond_lengths[i];
        }
        ++cnbs2;
    }
    if (num_particle_1_neighbours != cnbs2) {
        sprintf(error_message, "%s: part particle_1_bond_num %d does not equal cnbs2 %d \n", __func__, num_particle_1_neighbours, cnbs2);
        Error(error_message);
    }

    for(i=0; i<num_particle_1_neighbours; i++) {
        particle_1_neighbours[i] = sorted_particle_1_neighbours[i];
        particle_1_bond_lengths[i] = sorted_particle_1_bond_lengths[i];
    }

    free(sorted_particle_1_neighbours);
    free(sorted_particle_1_bond_lengths);
}