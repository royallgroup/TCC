#include "6Z.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

int check_for_common_spindle_particles(const int *first_5A_cluster, const int *second_5A_cluster);

void Clusters_Get6Z() {    // Detect 6Z clusters from 2 5A clusters
    int flg;
    int j, first_5A_spindle_id, k, l;
    int cnt;
    int s1a, s2a, s1b, s2b;
    int clusSize=6;

    s1a=s2a=s1b=s2b=-1;

    int *first_5A_cluster, *second_5A_cluster;
    int first_5A_id, second_5A_id;

    for (first_5A_id = 0; first_5A_id < nsp3c - 1; first_5A_id++) {
        first_5A_cluster = hcsp3c[first_5A_id];
        for (first_5A_spindle_id = 0; first_5A_spindle_id < 1; first_5A_spindle_id++) {
            for (j=0; j < nmem_sp3c[first_5A_cluster[first_5A_spindle_id]]; j++) {
                second_5A_id = mem_sp3c[first_5A_cluster[first_5A_spindle_id]][j];
                second_5A_cluster = hcsp3c[second_5A_id];

                if (second_5A_id > first_5A_id) {

                    if (check_for_common_spindle_particles(first_5A_cluster, second_5A_cluster) == 0) {

                        // for each cluster one spindle particle is in sp3 ring of other 5A and other spindle isn't
                        cnt = 0;    // check one 5A_i spindle are in sp3 ring of 5A_j
                        flg = first_5A_cluster[3] == second_5A_cluster[0] ||
                              first_5A_cluster[3] == second_5A_cluster[1] ||
                              first_5A_cluster[3] == second_5A_cluster[2];
                        if (flg == 1) {
                            s1a = first_5A_cluster[3];
                            s2a = first_5A_cluster[4];
                            ++cnt;
                        }
                        flg = first_5A_cluster[4] == second_5A_cluster[0] ||
                              first_5A_cluster[4] == second_5A_cluster[1] ||
                              first_5A_cluster[4] == second_5A_cluster[2];
                        if (flg == 1) {
                            s1a = first_5A_cluster[4];
                            s2a = first_5A_cluster[3];
                            ++cnt;
                        }
                        if (cnt != 1) continue;
                        cnt = 0;    // check one 5A_j spindle are in sp3 ring of 5A_i
                        flg = second_5A_cluster[3] == first_5A_cluster[0] ||
                              second_5A_cluster[3] == first_5A_cluster[1] ||
                              second_5A_cluster[3] == first_5A_cluster[2];
                        if (flg == 1) {
                            s1b = second_5A_cluster[3];
                            s2b = second_5A_cluster[4];
                            ++cnt;
                        }
                        flg = second_5A_cluster[4] == first_5A_cluster[0] ||
                              second_5A_cluster[4] == first_5A_cluster[1] ||
                              second_5A_cluster[4] == first_5A_cluster[2];
                        if (flg == 1) {
                            s1b = second_5A_cluster[4];
                            s2b = second_5A_cluster[3];
                            ++cnt;
                        }
                        if (cnt != 1) continue;
                        if (!Bonds_BondCheck(s1a, s1b))
                            continue;   // check spindles of i and j in sp3 ring of j and i respectively are bonded

                        cnt = 0;    // check 2 particles in the sp3 rings of 5A_i and 5A_j are common
                        for (k = 0; k < 3; ++k) {
                            for (l = 0; l < 3; ++l) {
                                if (first_5A_cluster[k] == second_5A_cluster[l]) {
                                    ++cnt;
                                    break;
                                }
                            }
                        }
                        if (cnt != 2) continue;

                        if (n6Z == m6Z) {
                            hc6Z = resize_2D_int(hc6Z, m6Z, m6Z + incrStatic, clusSize, -1);
                            m6Z = m6Z + incrStatic;
                        }
                        // Now we have found the 6Z cluster
                        if (s1a < s1b) {
                            hc6Z[n6Z][0] = s1a;    // insert cluster
                            hc6Z[n6Z][1] = s1b;
                            hc6Z[n6Z][2] = s2a;
                            hc6Z[n6Z][3] = s2b;
                        } else {
                            hc6Z[n6Z][0] = s1b;    // insert cluster
                            hc6Z[n6Z][1] = s1a;
                            hc6Z[n6Z][2] = s2b;
                            hc6Z[n6Z][3] = s2a;
                        }
                        cnt = 4;
                        for (k = 0; k < 3; ++k) {
                            flg = 1;
                            for (l = 0; l < 4; ++l) {
                                if (first_5A_cluster[k] == hc6Z[n6Z][l]) {
                                    flg = 0;
                                    break;
                                }
                            }
                            if (flg == 1) {
                                hc6Z[n6Z][cnt] = first_5A_cluster[k];
                                cnt++;
                            }
                        }
                        if (hc6Z[n6Z][5] < hc6Z[n6Z][4]) {
                            k = hc6Z[n6Z][5];
                            hc6Z[n6Z][5] = hc6Z[n6Z][4];
                            hc6Z[n6Z][4] = k;
                        }
                        Cluster_Write_6Z();
                    }
                }
            }
        }
    }
}

int check_for_common_spindle_particles(const int *first_5A_cluster, const int *second_5A_cluster) {
    if (first_5A_cluster[3] == second_5A_cluster[3]) return 1;
    else if(first_5A_cluster[3] == second_5A_cluster[4]) return 1;
    else if(first_5A_cluster[4] == second_5A_cluster[3]) return 1;
    else if (first_5A_cluster[4] == second_5A_cluster[4]) return 1;
    else return 0;
}

void Cluster_Write_6Z() {
    // hc6Z key: (5A_i_s_in_SP3_j, 5A_j_s_in_SP3_i, 5A_i_s_oth, 5A_j_s_oth, SP3_i_j_com_1, SP3_i_j_com_2)
    s6Z[hc6Z[n6Z][0]] = 'O';
    s6Z[hc6Z[n6Z][1]] = 'O';
    s6Z[hc6Z[n6Z][2]] = 'O';
    s6Z[hc6Z[n6Z][3]] = 'O';
    if (s6Z[hc6Z[n6Z][4]] == 'C') s6Z[hc6Z[n6Z][4]] = 'B';
    if (s6Z[hc6Z[n6Z][5]] == 'C') s6Z[hc6Z[n6Z][5]] = 'B';

    ++n6Z;
}