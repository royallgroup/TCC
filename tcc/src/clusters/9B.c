#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include <clusters/10B.h>
#include <clusters/11B.h>
#include <clusters/11E.h>
#include "9B.h"

void Clusters_Get9B_10B_11B_11E_12D() {

    //!  An 9B cluster is the intersection of 7A clusters shraing a spindle and two ring particles.
    /*!
   *  Find 9B clusters
   *  An 9B is 2 7A clusters where:
   *      - There is one common spindle.
   *      - The uncommon spindles are bonded.
   *      - Each distinct spindle is common with a ring particle of the other 7A
   *      - There are two common particles between the two 7A rings
   *
   *  Cluster output: BBBBBBOOS
   *  Storage order: SP5_lowerd_to_4, SP5_lowerd_to_5, SP5_higherd_to_4, SP5_higherd_to_5, SP5_i_j_com_lower, SP5_i_j_com_higher, sp5c_d1_lower, sp5c_d2_higher, s_com
   */

    int sp1, sp2i, sp2j;
    int sp5com[2];
    int k, l, m;
    int flg, fb1, fb2;
    int clusSize=9;

    sp1=sp2i=sp2j=-1;

    for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            int spindle_particle_id = first_7A_cluster[spindle_pointer + 5];
            for (int second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[spindle_particle_id]; ++second_7A_pointer) {
                int second_7A_id = mem_sp5c[spindle_particle_id][second_7A_pointer];
                int *second_7A_cluster = hcsp5c[second_7A_id];

                flg = 0;
                if (first_7A_cluster[5] == second_7A_cluster[5] && first_7A_cluster[6] != second_7A_cluster[6]) {
                    if (Bonds_BondCheck(first_7A_cluster[6], second_7A_cluster[6])) { // spindle particles arranged
                        flg = 1;
                        sp1 = first_7A_cluster[5];   // s_com common spindle
                        sp2i = first_7A_cluster[6];  // 2nd spindle particle of cluster 7A_i
                        sp2j = second_7A_cluster[6];  // 2nd spindle particle of cluster 7A_j
                    }
                }
                if (first_7A_cluster[6] == second_7A_cluster[6] && first_7A_cluster[5] != second_7A_cluster[5]) {
                    if (Bonds_BondCheck(first_7A_cluster[5], second_7A_cluster[5])) {
                        flg = 1;
                        sp1 = first_7A_cluster[6];   // s_com common spindle
                        sp2i = first_7A_cluster[5];  // 2nd spindle particle of cluster 7A_i
                        sp2j = second_7A_cluster[5];  // 2nd spindle particle of cluster 7A_j
                    }
                }
                if (first_7A_cluster[5] == second_7A_cluster[6] && first_7A_cluster[6] != second_7A_cluster[5]) {
                    if (Bonds_BondCheck(first_7A_cluster[6], second_7A_cluster[5])) {
                        flg = 1;
                        sp1 = first_7A_cluster[5];   // s_com common spindle
                        sp2i = first_7A_cluster[6];  // 2nd spindle particle of cluster 7A_i
                        sp2j = second_7A_cluster[5];  // 2nd spindle particle of cluster 7A_j
                    }
                }
                if (first_7A_cluster[6] == second_7A_cluster[5] && first_7A_cluster[5] != second_7A_cluster[6]) {
                    if (Bonds_BondCheck(first_7A_cluster[5], second_7A_cluster[6])) {
                        flg = 1;
                        sp1 = first_7A_cluster[6];   // s_com common spindle
                        sp2i = first_7A_cluster[5];  // 2nd spindle particle of cluster 7A_i
                        sp2j = second_7A_cluster[6];  // 2nd spindle particle of cluster 7A_j
                    }
                }
                if (flg == 0) continue;

                fb1 = fb2 = 1;  // ensure the two distinct spindle particles are part of the other SP5 ring
                for (k = 0; k < 5; ++k) {
                    if (sp2i == second_7A_cluster[k]) fb1 = 0;
                    if (sp2j == first_7A_cluster[k]) fb2 = 0;
                }
                if (fb1 || fb2) continue;

                m = 0;  // check for two common SP5 particles
                for (k = 0; k < 5; ++k) {
                    for (l = 0; l < 5; ++l) {
                        if (first_7A_cluster[k] == second_7A_cluster[l]) {
                            if (m == 2) {
                                m++;
                                break;
                            }
                            sp5com[m] = first_7A_cluster[k];
                            ++m;
                        }
                    }
                }
                if (m != 2) continue;

                // Now we have found the 9B C2v cluster
                if (n9B == m9B) {
                    hc9B = resize_2D_int(hc9B, m9B, m9B + incrStatic, clusSize, -1);
                    m9B = m9B + incrStatic;
                }
                if (sp5com[0] < sp5com[1]) {
                    hc9B[n9B][4] = sp5com[0];
                    hc9B[n9B][5] = sp5com[1];
                } else {
                    hc9B[n9B][4] = sp5com[1];
                    hc9B[n9B][5] = sp5com[0];
                }

                if (sp2i < sp2j) {
                    hc9B[n9B][6] = sp2i;
                    hc9B[n9B][7] = sp2j;

                    for (k = 0; k < 5; ++k) {
                        if (Bonds_BondCheck(first_7A_cluster[k], hc9B[n9B][4]) && first_7A_cluster[k] != hc9B[n9B][7] &&
                            first_7A_cluster[k] !=
                            hc9B[n9B][4]) {
                            hc9B[n9B][0] = first_7A_cluster[k];
                        }
                        if (Bonds_BondCheck(first_7A_cluster[k], hc9B[n9B][5]) && first_7A_cluster[k] != hc9B[n9B][7] &&
                            first_7A_cluster[k] !=
                            hc9B[n9B][5]) {
                            hc9B[n9B][1] = first_7A_cluster[k];
                        }
                        if (Bonds_BondCheck(second_7A_cluster[k], hc9B[n9B][4]) &&
                            second_7A_cluster[k] != hc9B[n9B][6] && second_7A_cluster[k] !=
                                                                    hc9B[n9B][4]) {
                            hc9B[n9B][2] = second_7A_cluster[k];
                        }
                        if (Bonds_BondCheck(second_7A_cluster[k], hc9B[n9B][5]) &&
                            second_7A_cluster[k] != hc9B[n9B][6] && second_7A_cluster[k] !=
                                                                    hc9B[n9B][5]) {
                            hc9B[n9B][3] = second_7A_cluster[k];
                        }
                    }
                } else {
                    hc9B[n9B][6] = sp2j;
                    hc9B[n9B][7] = sp2i;

                    for (k = 0; k < 5; ++k) {
                        if (Bonds_BondCheck(second_7A_cluster[k], hc9B[n9B][4]) &&
                            second_7A_cluster[k] != hc9B[n9B][7] && second_7A_cluster[k] !=
                                                                    hc9B[n9B][4]) {
                            hc9B[n9B][0] = second_7A_cluster[k];
                        }
                        if (Bonds_BondCheck(second_7A_cluster[k], hc9B[n9B][5]) &&
                            second_7A_cluster[k] != hc9B[n9B][7] && second_7A_cluster[k] !=
                                                                    hc9B[n9B][5]) {
                            hc9B[n9B][1] = second_7A_cluster[k];
                        }
                        if (Bonds_BondCheck(first_7A_cluster[k], hc9B[n9B][4]) && first_7A_cluster[k] != hc9B[n9B][6] &&
                            first_7A_cluster[k] !=
                            hc9B[n9B][4]) {
                            hc9B[n9B][2] = first_7A_cluster[k];
                        }
                        if (Bonds_BondCheck(first_7A_cluster[k], hc9B[n9B][5]) && first_7A_cluster[k] != hc9B[n9B][6] &&
                            first_7A_cluster[k] !=
                            hc9B[n9B][5]) {
                            hc9B[n9B][3] = first_7A_cluster[k];
                        }
                    }
                }
                hc9B[n9B][8] = sp1;
                Cluster_Write_9B();

                if (do10B == 1) Clusters_Get10B(first_7A_id, second_7A_id);
                if (do11B == 1) {
                    if (Clusters_Get11B()) {
                        s11B[hc9B[n9B][8]] = 'S';
                        ++n11B;
                    }
                }
                if (do11E == 1) Clusters_Get11E_12D(first_7A_id, second_7A_id, sp1, sp2i, sp2j);

                ++n9B;
            }
        }
    }
}

void Cluster_Write_9B() {
    int i;
    for(i=0; i<6; i++) {
        if (s9B[hc9B[n9B][i]] == 'C') s9B[hc9B[n9B][i]] = 'B';
    }
    if (s9B[hc9B[n9B][6]] != 'S') s9B[hc9B[n9B][6]] = 'O';
    if (s9B[hc9B[n9B][7]] != 'S') s9B[hc9B[n9B][7]] = 'O';
    s9B[hc9B[n9B][8]] = 'S';
}