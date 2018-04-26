#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include <clusters/simple_cluster_methods.h>
#include "9A.h"

int count_bonded_ring_particles(const int *first_sp4b_cluster, const int *second_sp4b_cluster, const int *db, int *ob);

void Clusters_Get9A() {    // Detect 9A D3h clusters
    int first_sp4b_id, second_sp4b_pointer, first_sp4b_spindle_pointer, third_sp4b_id, l, m, n;
    int i_j_common_ring_particles[8], i_j_uncommon_ring_particles[8];
    int flg;
    int clusSize = 9;

    for (first_sp4b_id = 0; first_sp4b_id < nsp4b - 2; ++first_sp4b_id) {  // loop over all sp4b_i
        int *first_sp4b_cluster = hcsp4b[first_sp4b_id];
        for (first_sp4b_spindle_pointer = 0; first_sp4b_spindle_pointer < 1; first_sp4b_spindle_pointer++) {
            for (second_sp4b_pointer = 0; second_sp4b_pointer <
                                          nmem_sp4b[first_sp4b_cluster[first_sp4b_spindle_pointer]]; ++second_sp4b_pointer) { // loop over all sp4b_j
                int second_sp4b_id = mem_sp4b[first_sp4b_cluster[first_sp4b_spindle_pointer]][second_sp4b_pointer];
                int *second_sp4b_cluster = hcsp4b[second_sp4b_id];
                if (second_sp4b_id > first_sp4b_id) {
                    // Spindles should not be common
                    if (first_sp4b_cluster[4] != second_sp4b_cluster[4]) {
                        // Spindles should not be bonded
                        if (Bonds_BondCheck(first_sp4b_cluster[4], second_sp4b_cluster[4]) == 0) {

                            // Check there are two common particles between the rings of the two clusters
                            if (count_common_ring_particles(first_sp4b_cluster, second_sp4b_cluster, 4, 4, i_j_common_ring_particles) == 2) {

                                if (count_bonded_ring_particles(first_sp4b_cluster, second_sp4b_cluster, i_j_common_ring_particles, i_j_uncommon_ring_particles) != 4) continue;
                                // POSSIBLE IMPROVEMENT!! could make detection faster here by looping over sp4b clusters bonded to uncommon particles to sp4b_i
                                for (third_sp4b_id = second_sp4b_id + 1; third_sp4b_id < nsp4b; third_sp4b_id++) {
                                    int *third_sp4b_cluster = hcsp4b[third_sp4b_id];
                                    // ERROR!! need to check spindle of sp4b_k distinct from spindles of sp4b_i and sp4b_j and no bonds between these three particles
                                    int tmp[4];

                                    if (count_common_ring_particles(third_sp4b_cluster, i_j_uncommon_ring_particles, 4, 4, tmp) == 4) {

                                        // Now we have found the 9A D3h cluster
                                        if (n9A == m9A) {
                                            hc9A = resize_2D_int(hc9A, m9A, m9A + incrStatic, clusSize, -1);
                                            m9A = m9A + incrStatic;
                                        }
                                        if (first_sp4b_cluster[4] < second_sp4b_cluster[4] &&
                                            first_sp4b_cluster[4] < third_sp4b_cluster[4]) {
                                            hc9A[n9A][6] = first_sp4b_cluster[4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                            hc9A[n9A][0] = first_sp4b_cluster[0];
                                            hc9A[n9A][1] = first_sp4b_cluster[1];
                                            hc9A[n9A][2] = first_sp4b_cluster[2];
                                            hc9A[n9A][3] = first_sp4b_cluster[3];

                                            if (Bonds_BondCheck(second_sp4b_cluster[4], hc9A[n9A][0])) {
                                                hc9A[n9A][7] = second_sp4b_cluster[4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                                hc9A[n9A][8] = third_sp4b_cluster[4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                            } else {
                                                hc9A[n9A][7] = third_sp4b_cluster[4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                                hc9A[n9A][8] = second_sp4b_cluster[4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                            }

                                            for (l = 0; l < 4; l++) {
                                                flg = 1;
                                                for (m = 0; m < 4; m++) {
                                                    if (i_j_uncommon_ring_particles[l] == hc9A[n9A][m]) {
                                                        flg = 0;
                                                        break;
                                                    }
                                                }
                                                if (flg == 1 &&
                                                    Bonds_BondCheck(i_j_uncommon_ring_particles[l], hc9A[n9A][0])) {
                                                    hc9A[n9A][4] = i_j_uncommon_ring_particles[l];
                                                } else if (flg == 1) {
                                                    hc9A[n9A][5] = i_j_uncommon_ring_particles[l];
                                                }
                                            }
                                        } else if (second_sp4b_cluster[4] < first_sp4b_cluster[4] &&
                                                   second_sp4b_cluster[4] <
                                                   third_sp4b_cluster[4]) {
                                            hc9A[n9A][6] = second_sp4b_cluster[4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                            hc9A[n9A][0] = second_sp4b_cluster[0];
                                            hc9A[n9A][1] = second_sp4b_cluster[1];
                                            hc9A[n9A][2] = second_sp4b_cluster[2];
                                            hc9A[n9A][3] = second_sp4b_cluster[3];

                                            if (Bonds_BondCheck(first_sp4b_cluster[4], hc9A[n9A][0])) {
                                                hc9A[n9A][7] = first_sp4b_cluster[4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                                hc9A[n9A][8] = third_sp4b_cluster[4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                            } else {
                                                hc9A[n9A][7] = third_sp4b_cluster[4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                                hc9A[n9A][8] = first_sp4b_cluster[4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                            }

                                            for (l = 0; l < 4; l++) {
                                                flg = 1;
                                                for (m = 0; m < 4; m++) {
                                                    if (i_j_uncommon_ring_particles[l] == hc9A[n9A][m]) {
                                                        flg = 0;
                                                        break;
                                                    }
                                                }
                                                if (flg == 1 &&
                                                    Bonds_BondCheck(i_j_uncommon_ring_particles[l], hc9A[n9A][0])) {
                                                    hc9A[n9A][4] = i_j_uncommon_ring_particles[l];
                                                } else if (flg == 1) {
                                                    hc9A[n9A][5] = i_j_uncommon_ring_particles[l];
                                                }
                                            }
                                        } else {
                                            hc9A[n9A][6] = third_sp4b_cluster[4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                            hc9A[n9A][0] = third_sp4b_cluster[0];
                                            hc9A[n9A][1] = third_sp4b_cluster[1];
                                            hc9A[n9A][2] = third_sp4b_cluster[2];
                                            hc9A[n9A][3] = third_sp4b_cluster[3];

                                            if (Bonds_BondCheck(first_sp4b_cluster[4], hc9A[n9A][0])) {
                                                hc9A[n9A][7] = first_sp4b_cluster[4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                                hc9A[n9A][8] = second_sp4b_cluster[4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                            } else {
                                                hc9A[n9A][7] = second_sp4b_cluster[4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                                hc9A[n9A][8] = first_sp4b_cluster[4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                            }

                                            if (Bonds_BondCheck(i_j_common_ring_particles[0], hc9A[n9A][0])) {
                                                hc9A[n9A][4] = i_j_common_ring_particles[0];
                                                hc9A[n9A][5] = i_j_common_ring_particles[1];
                                            } else {
                                                hc9A[n9A][4] = i_j_common_ring_particles[1];
                                                hc9A[n9A][5] = i_j_common_ring_particles[0];
                                            }
                                        }
                                        Cluster_Write_9A();
                                        break;

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

int count_bonded_ring_particles(const int *first_sp4b_cluster, const int *second_sp4b_cluster, const int *db, int *ob) {

    int num_bonded_ring_particles = 0;
    for (int first_cluster_pointer = 0; first_cluster_pointer < 4; ++first_cluster_pointer) {
        // find particles in SP4 ring of sp4b_i not common to SP4 ring of sp4b_j
        if (first_sp4b_cluster[first_cluster_pointer] != db[0] && first_sp4b_cluster[first_cluster_pointer] != db[1]) {
            for (int second_cluster_pointer = 0; second_cluster_pointer < 4; ++second_cluster_pointer) {
                // find particles in SP4 ring of sp4b_j not common to SP4 ring of sp4b_i
                if (second_sp4b_cluster[second_cluster_pointer] != db[0] && second_sp4b_cluster[second_cluster_pointer] != db[1]) {
                    // check non-common SP4 ring particles from sp4b_i and sp4b_j are bonded
                    if (Bonds_BondCheck(first_sp4b_cluster[first_cluster_pointer], second_sp4b_cluster[second_cluster_pointer])) {
                        ob[num_bonded_ring_particles] = first_sp4b_cluster[first_cluster_pointer];
                        num_bonded_ring_particles++;
                        ob[num_bonded_ring_particles] = second_sp4b_cluster[second_cluster_pointer];
                        num_bonded_ring_particles++;
                    }
                }
            }
        }
    }
    return num_bonded_ring_particles;
}

void Cluster_Write_9A() {
    // hc9A key: (SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_to_0_in_SP4_lowest_s, SP4_to_1_in_SP4_lowest_s, s_lowest, s_to_0_in_SP4_lowest_s,s_to_2_in_SP4_lowest_s)
    if (s9A[hc9A[n9A][0]] == 'C') s9A[hc9A[n9A][0]] = 'B';
    if (s9A[hc9A[n9A][1]] == 'C') s9A[hc9A[n9A][1]] = 'B';
    if (s9A[hc9A[n9A][3]] == 'C') s9A[hc9A[n9A][3]] = 'B';
    if (s9A[hc9A[n9A][4]] == 'C') s9A[hc9A[n9A][4]] = 'B';
    if (s9A[hc9A[n9A][6]] == 'C') s9A[hc9A[n9A][6]] = 'B';
    if (s9A[hc9A[n9A][7]] == 'C') s9A[hc9A[n9A][7]] = 'B';
    s9A[hc9A[n9A][2]] = 'O';
    s9A[hc9A[n9A][5]] = 'O';
    s9A[hc9A[n9A][8]] = 'O';
    ++n9A;
}