#include "simple_cluster_methods.h"
#include "13B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get13B() {
    int common_spindle_id[2];
    int uncommon_spindle_ids[2];

    for(int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id){
        int *first_7A_cluster = hcsp5c[first_7A_id];

        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            int spindle_id = first_7A_cluster[5 + spindle_pointer];
            for (int second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[spindle_id]; second_7A_pointer++) {
                int second_7A_id = mem_sp5c[spindle_id][second_7A_pointer];
                if (first_7A_id > second_7A_id) {
                    int *second_7A_cluster = hcsp5c[second_7A_id];

                    if (count_common_spindle_particles(first_7A_cluster, second_7A_cluster, 7, 7, common_spindle_id) == 1) {
                        get_uncommon_spindle_particles(first_7A_cluster, second_7A_cluster, common_spindle_id[0], uncommon_spindle_ids);

                        if (Bonds_BondCheck(uncommon_spindle_ids[0], uncommon_spindle_ids[1]) == 0) {

                            if (check_rings_are_uncommon(first_7A_cluster, second_7A_cluster) == 1) {

                                if (check_rings_are_bonded(first_7A_cluster, second_7A_cluster) == 1) {

                                    if (check_rings_are_bonded(second_7A_cluster, first_7A_cluster) == 1) {

                                        Cluster_Write_13B(first_7A_cluster, second_7A_cluster, common_spindle_id[0], uncommon_spindle_ids);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int check_rings_are_bonded(const int *first_7A_cluster, const int *second_7A_cluster) {
    // If each particle in 7Ai is bonded to exactly 1 particle in 7Aj returns 1 else returns 0
    int first_ring_pointer;

    for (first_ring_pointer = 0; first_ring_pointer < 5; ++first_ring_pointer) {
        int first_ring_particle = first_7A_cluster[first_ring_pointer];

        if (count_bonds_to_ring(first_ring_particle, second_7A_cluster) != 1) {
            break;
        }
    }
    if(first_ring_pointer == 5) {
        return 1;
    } else {
        return 0;
    }
}

int count_bonds_to_ring(int particle_id, const int *first_7A_cluster) {
    // Counts the number of bonds between particle and the ring of a 7A cluster
    int num_bonds = 0;
    for (int second_ring_pointer = 0; second_ring_pointer < 5; ++second_ring_pointer) {
        if (Bonds_BondCheck(particle_id, first_7A_cluster[second_ring_pointer]) == 1) {
            ++num_bonds;
        }
    }
    return num_bonds;
}

int check_rings_are_uncommon(const int *first_7A_cluster, const int *second_7A_cluster) {
    // Return 1 if all ring particles in first 7A are distinct from all ring particles in second 7A

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (first_7A_cluster[i] == second_7A_cluster[j]) {
            return 0;
            }
        }
    }
    return 1;
}

void get_uncommon_spindle_particles(const int *first_7A_cluster, const int *second_7A_cluster, int common_spindle_id, int *uncommon_spindle_ids) {
    if (first_7A_cluster[5] == common_spindle_id) {
        uncommon_spindle_ids[0] = first_7A_cluster[6];
    } else {
        uncommon_spindle_ids[0] = first_7A_cluster[5];
    }
    if (second_7A_cluster[5] == common_spindle_id) {
        uncommon_spindle_ids[1] = second_7A_cluster[6];
    } else {
        uncommon_spindle_ids[1] = second_7A_cluster[5];
    }
}

void Cluster_Write_13B(const int *first_7A_cluster, const int *second_7A_cluster, int common_spindle_id, const int *uncommon_spindle_ids) {
    int clusSize = 13;

    if (n13B == m13B) {
        hc13B = resize_2D_int(hc13B, m13B, m13B + incrStatic, clusSize, -1);
        m13B = m13B + incrStatic;
    }

    hc13B[n13B][0] = common_spindle_id;

    hc13B[n13B][1] = uncommon_spindle_ids[0];
    hc13B[n13B][2] = uncommon_spindle_ids[1];

    for (int i = 0; i < 5; ++i) {
        hc13B[n13B][i + 3] = first_7A_cluster[i];
    }
    for (int i = 0; i < 5; ++i) {
        hc13B[n13B][i + 8] = second_7A_cluster[i];
    }

    s13B[hc13B[n13B][0]] = 'S';
    if(s13B[hc13B[n13B][1]] != 'S') s13B[hc13B[n13B][1]] = 'O';
    if(s13B[hc13B[n13B][2]] != 'S') s13B[hc13B[n13B][2]] = 'O';
    for(int i = 3; i < 13; ++i) {
        if (s13B[hc13B[n13B][i]] == 'C') s13B[hc13B[n13B][i]] = 'B';
    }

    ++n13B;
}