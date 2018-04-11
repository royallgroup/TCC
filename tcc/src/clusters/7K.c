#include "7K.h"
#include "globals.h"
#include "tools.h"

int count_common_spindles_between_5As(const int *first_5A_cluster, const int *second_5A_cluster, int *scom);

void get_other_spindle_ids(const int *first_5A_cluster, const int *second_5A_cluster, int scom, int *sother);

void Clusters_Get7K() {    // Detect 7K clusters from 2 5A clusters
    int first_5A_id, second_5A_pointer, first_5A_spindle_pointer, k, l, m;
    int second_5A_id;
    int *first_5A_cluster, *second_5A_cluster;
    int scom, sother[2], sp3_com[2], sp3c_i_other, sp3c_j_other;
    int clusSize=7;
    int common_spindles;

    scom=sother[0]=sother[1]=sp3_com[0]=sp3_com[1]=sp3c_i_other=sp3c_j_other=-1;

    for (first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {  // loop over all 5A_i
        first_5A_cluster = hcsp3c[first_5A_id];
        for (first_5A_spindle_pointer = 3; first_5A_spindle_pointer < 5; ++first_5A_spindle_pointer) {    // loop over both spindles of 5A_i
            for (second_5A_pointer = 0; second_5A_pointer < nmem_sp3c[first_5A_cluster[first_5A_spindle_pointer]]; ++second_5A_pointer) {  // loop over all 5A_j common with spindle of 5A_i
                second_5A_id = mem_sp3c[first_5A_cluster[first_5A_spindle_pointer]][second_5A_pointer];
                second_5A_cluster = hcsp3c[second_5A_id];
                if (second_5A_id > first_5A_id) {
                    // exactly one common spindle between 5A_i and 5A_j
                    if (count_common_spindles_between_5As(first_5A_cluster, second_5A_cluster, &scom) == 1) {

                        get_other_spindle_ids(first_5A_cluster, second_5A_cluster, scom, sother);

                        m = 0;
                        for (k = 0; k < 5; k++) {
                            if (sother[0] == second_5A_cluster[k]) {
                                m++;
                                break;
                            }
                        }
                        if (m != 0 && k != 5) continue; // other spindle of 5A_i distinct from whole 5A_j

                        m = 0;
                        for (k = 0; k < 5; k++) {
                            if (sother[1] == first_5A_cluster[k]) {
                                m++;
                                break;
                            }
                        }
                        if (m != 0 && k != 5) continue; // other spindle of 5A_j distinct from whole 5A_i

                        m = 0;
                        for (k = 0; k < 3; k++) {
                            for (l = 0; l < 3; l++) {
                                if (first_5A_cluster[k] == second_5A_cluster[l]) {
                                    if (m >= 2) {
                                        m++;
                                        break;
                                    }
                                    sp3_com[m] = first_5A_cluster[k];
                                    m++;
                                }
                            }
                            if (m >= 3) break;
                        }
                        if (m != 2) continue; // exactly two common particles in SP3 rings of 5A_i and 5A_j

                        m = 0;
                        for (k = 0; k < 3; k++) {
                            for (l = 0; l < 2; l++) {
                                if (first_5A_cluster[k] == sp3_com[l]) break;
                            }
                            if (l == 2) {
                                if (m >= 1) {
                                    m++;
                                    break;
                                }
                                sp3c_i_other = first_5A_cluster[k];
                                m++;
                            }
                        }
                        if (m != 1) continue; // found other uncommon particle from SP3 ring of 5A_i

                        m = 0;
                        for (k = 0; k < 3; k++) {
                            for (l = 0; l < 2; l++) {
                                if (second_5A_cluster[k] == sp3_com[l]) break;
                            }
                            if (l == 2) {
                                if (m >= 1) {
                                    m++;
                                    break;
                                }
                                sp3c_j_other = second_5A_cluster[k];
                                m++;
                            }
                        }
                        if (m != 1) continue; // found other uncommon particle from SP3 ring of 5A_i

                        m = 0;
                        for (k = 0; k < 5; k++) {
                            if (sp3c_i_other == second_5A_cluster[k]) {
                                m++;
                                break;
                            }
                        }
                        if (m != 0 && k != 5) continue; // other ring of 5A_i distinct from whole 5A_j

                        m = 0;
                        for (k = 0; k < 5; k++) {
                            if (sp3c_j_other == first_5A_cluster[k]) {
                                m++;
                                break;
                            }
                        }
                        if (m != 0 && k != 5) continue; // other ring of 5A_j distinct from whole 5A_i

                        if (n7K == m7K) {
                            hc7K = resize_2D_int(hc7K, m7K, m7K + incrStatic, clusSize, -1);
                            m7K = m7K + incrStatic;
                        }
                        // Now we have found the 7K cluster

                        hc7K[n7K][0] = scom;
                        hc7K[n7K][1] = sother[0];
                        hc7K[n7K][2] = sother[1];
                        hc7K[n7K][3] = sp3_com[0];
                        hc7K[n7K][4] = sp3_com[1];
                        hc7K[n7K][5] = sp3c_i_other;
                        hc7K[n7K][6] = sp3c_j_other;

                        quickSort(&hc7K[n7K][1], 2);
                        quickSort(&hc7K[n7K][3], 2);
                        quickSort(&hc7K[n7K][5], 2);

                        Cluster_Write_7K();
                    }
                }
            }
        }
    }
}

void get_other_spindle_ids(const int *first_5A_cluster, const int *second_5A_cluster, int scom, int *sother) {
    if (first_5A_cluster[3] == scom) {
        sother[0] = first_5A_cluster[4];
    }
    else {
        sother[0] = first_5A_cluster[3];
    }
    if (second_5A_cluster[3] == scom) {
        sother[1] = second_5A_cluster[4];
    }
    else {
        sother[1] = second_5A_cluster[3];
    }
}

int count_common_spindles_between_5As(const int *first_5A_cluster, const int *second_5A_cluster, int *scom) {
    int num_common_spindles = 0;
    for (int ring_1_pointer = 3; ring_1_pointer < 5; ring_1_pointer++) {
        for (int ring_2_pointer = 3; ring_2_pointer < 5; ring_2_pointer++) {
            if (first_5A_cluster[ring_1_pointer] == second_5A_cluster[ring_2_pointer]) {
                *scom = first_5A_cluster[ring_1_pointer];
                num_common_spindles++;
            }
        }
    }
    return num_common_spindles;
}

void Cluster_Write_7K() {
    // hc7K key: (scom, sother, ring_com, ring_other)
    s7K[hc7K[n7K][0]] = 'O';
    s7K[hc7K[n7K][1]] = 'O';
    s7K[hc7K[n7K][2]] = 'O';
    if (s7K[hc7K[n7K][3]] == 'C') s7K[hc7K[n7K][3]] = 'B';
    if (s7K[hc7K[n7K][4]] == 'C') s7K[hc7K[n7K][4]] = 'B';
    if (s7K[hc7K[n7K][5]] == 'C') s7K[hc7K[n7K][5]] = 'B';
    if (s7K[hc7K[n7K][6]] == 'C') s7K[hc7K[n7K][6]] = 'B';

    ++n7K;
}