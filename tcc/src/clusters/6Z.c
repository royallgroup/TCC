#include "6Z.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get6Z() {    // Detect 6Z clusters from 2 5A clusters

    int j, first_5A_spindle_id;

    int first_spindles[2], second_spindles[2], common_ring_particles[2];


    int *first_5A_cluster, *second_5A_cluster;
    int first_5A_id, second_5A_id;

    for (first_5A_id = 0; first_5A_id < nsp3c - 1; first_5A_id++) {
        first_5A_cluster = hcsp3c[first_5A_id];
        for (first_5A_spindle_id = 0; first_5A_spindle_id < 1; first_5A_spindle_id++) {
            for (j = 0; j < nmem_sp3c[first_5A_cluster[first_5A_spindle_id]]; j++) {
                second_5A_id = mem_sp3c[first_5A_cluster[first_5A_spindle_id]][j];
                second_5A_cluster = hcsp3c[second_5A_id];
                if (second_5A_id > first_5A_id) {
                    if (check_for_common_spindle_particles(first_5A_cluster, second_5A_cluster) == 0) {
                        if (count_spindles_in_ring(first_5A_cluster, second_5A_cluster, first_spindles) == 1) {
                            if (count_spindles_in_ring(second_5A_cluster, first_5A_cluster, second_spindles) == 1) {
                                if (Bonds_BondCheck(first_spindles[0], second_spindles[0]) == 1) {
                                    if (get_bonds_between_rings(first_5A_cluster, second_5A_cluster, common_ring_particles) == 2) {

                                        Cluster_Write_6Z(first_spindles, second_spindles, common_ring_particles);
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

int check_for_common_spindle_particles(const int *first_5A_cluster, const int *second_5A_cluster) {
    if (first_5A_cluster[3] == second_5A_cluster[3]) return 1;
    else if (first_5A_cluster[3] == second_5A_cluster[4]) return 1;
    else if (first_5A_cluster[4] == second_5A_cluster[3]) return 1;
    else if (first_5A_cluster[4] == second_5A_cluster[4]) return 1;
    else return 0;
}

int count_spindles_in_ring(const int *first_5A_cluster, const int *second_5A_cluster, int *spindles) {

    int spindles_in_ring = 0;
    spindles_in_ring += first_5A_cluster[3] == second_5A_cluster[0];
    spindles_in_ring += first_5A_cluster[3] == second_5A_cluster[1];
    spindles_in_ring += first_5A_cluster[3] == second_5A_cluster[2];

    if (spindles_in_ring == 1) {
        spindles[0] = first_5A_cluster[3];
        spindles[1] = first_5A_cluster[4];
        return 1;
    }

    spindles_in_ring += first_5A_cluster[4] == second_5A_cluster[0];
    spindles_in_ring += first_5A_cluster[4] == second_5A_cluster[1];
    spindles_in_ring += first_5A_cluster[4] == second_5A_cluster[2];

    if (spindles_in_ring == 1) {
        spindles[0] = first_5A_cluster[4];
        spindles[1] = first_5A_cluster[3];
        return 1;
    }
    else {
        return 0;
    }
}

int get_bonds_between_rings(const int *first_5A_cluster, const int *second_5A_cluster, int *common_ring_particles) {
    int bonds_between_rings = 0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (first_5A_cluster[i] == second_5A_cluster[j]) {
                common_ring_particles[bonds_between_rings] = first_5A_cluster[i];
                bonds_between_rings++;
                break;
            }
        }
    }
    return bonds_between_rings;
}

void Cluster_Write_6Z(const int *first_spindles, const int *second_spindles, const int *common_ring_particles) {
    int clusSize = 6;

    if (n6Z == m6Z) {
        hc6Z = resize_2D_int(hc6Z, m6Z, m6Z + incrStatic, clusSize, -1);
        m6Z = m6Z + incrStatic;
    }

    hc6Z[n6Z][0] = first_spindles[0];
    hc6Z[n6Z][1] = first_spindles[1];
    hc6Z[n6Z][2] = second_spindles[0];
    hc6Z[n6Z][3] = second_spindles[1];
    hc6Z[n6Z][4] = common_ring_particles[0];
    hc6Z[n6Z][5] = common_ring_particles[1];

    // hc6Z key: (5A spindles * 4, 5A ring particles * 2)
    s6Z[hc6Z[n6Z][0]] = 'O';
    s6Z[hc6Z[n6Z][1]] = 'O';
    s6Z[hc6Z[n6Z][2]] = 'O';
    s6Z[hc6Z[n6Z][3]] = 'O';
    if (s6Z[hc6Z[n6Z][4]] == 'C') s6Z[hc6Z[n6Z][4]] = 'B';
    if (s6Z[hc6Z[n6Z][5]] == 'C') s6Z[hc6Z[n6Z][5]] = 'B';

    ++n6Z;
}