#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include <clusters/10B.h>
#include <clusters/11B.h>
#include <clusters/11E.h>
#include "9B.h"

void Clusters_Get9B_10B_11B_11E_12D() {    // Detect 9B, 10B, 11B, 11E & 12D
    int sp1, sp2i, sp2j;
    int sp5com[2];
    int i, j, k, l, m;
    int flg, fb1, fb2;
    int clusSize=9;

    sp1=sp2i=sp2j=-1;

    for (i=0; i < nsp5c - 1; ++i) {  // loop over all 7A_i
        // POSSIBLE IMPROVEMENT!! - 2 loops: over all 7A clusters which each spindle is in
        for (j=i+1; j < nsp5c; ++j) {  // loop over all 7A_j
            flg = 0;
            if (hcsp5c[i][5] == hcsp5c[j][5] && hcsp5c[i][6] != hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][6])) { // spindle particles arranged
                    flg = 1;
                    sp1 = hcsp5c[i][5];   // s_com common spindle
                    sp2i = hcsp5c[i][6];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][6];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (hcsp5c[i][6] == hcsp5c[j][6] && hcsp5c[i][5] != hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][5], hcsp5c[j][5])) {
                    flg = 1;
                    sp1 = hcsp5c[i][6];   // s_com common spindle
                    sp2i = hcsp5c[i][5];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][5];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (hcsp5c[i][5] == hcsp5c[j][6] && hcsp5c[i][6] != hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][5])) {
                    flg = 1;
                    sp1 = hcsp5c[i][5];   // s_com common spindle
                    sp2i = hcsp5c[i][6];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][5];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (hcsp5c[i][6] == hcsp5c[j][5] && hcsp5c[i][5] != hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][5], hcsp5c[j][6])) {
                    flg = 1;
                    sp1 = hcsp5c[i][6];   // s_com common spindle
                    sp2i = hcsp5c[i][5];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][6];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (flg==0) continue;

            fb1 = fb2 = 1;  // ensure the two distinct spindle particles are part of the other SP5 ring
            for (k=0; k<5; ++k) {
                if (sp2i == hcsp5c[j][k]) fb1 = 0;
                if (sp2j == hcsp5c[i][k]) fb2 = 0;
            }
            if (fb1 || fb2) continue;

            m = 0;  // check for two common SP5 particles
            for (k=0; k<5; ++k) {
                for (l=0; l<5; ++l) {
                    if (hcsp5c[i][k] == hcsp5c[j][l]) {
                        if (m==2) {m++; break; }
                        sp5com[m]= hcsp5c[i][k];
                        ++m;
                    }
                }
            }
            if (m!=2) continue;

            // Now we have found the 9B C2v cluster
            if (n9B == m9B) {
                hc9B= resize_2D_int(hc9B, m9B, m9B + incrStatic, clusSize, -1);
                m9B= m9B + incrStatic;
            }
            if (sp5com[0]<sp5com[1]) {
                hc9B[n9B][4]=sp5com[0];
                hc9B[n9B][5]=sp5com[1];
            }
            else {
                hc9B[n9B][4]=sp5com[1];
                hc9B[n9B][5]=sp5com[0];
            }

            if (sp2i<sp2j) {
                hc9B[n9B][6]=sp2i;
                hc9B[n9B][7]=sp2j;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[i][k], hc9B[n9B][4]) && hcsp5c[i][k] != hc9B[n9B][7] && hcsp5c[i][k] !=
                                                                                                       hc9B[n9B][4]) {
                        hc9B[n9B][0]= hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k], hc9B[n9B][5]) && hcsp5c[i][k] != hc9B[n9B][7] && hcsp5c[i][k] !=
                                                                                                       hc9B[n9B][5]) {
                        hc9B[n9B][1]= hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[j][k], hc9B[n9B][4]) && hcsp5c[j][k] != hc9B[n9B][6] && hcsp5c[j][k] !=
                                                                                                       hc9B[n9B][4]) {
                        hc9B[n9B][2]= hcsp5c[j][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[j][k], hc9B[n9B][5]) && hcsp5c[j][k] != hc9B[n9B][6] && hcsp5c[j][k] !=
                                                                                                       hc9B[n9B][5]) {
                        hc9B[n9B][3]= hcsp5c[j][k];
                    }
                }
            }
            else {
                hc9B[n9B][6]=sp2j;
                hc9B[n9B][7]=sp2i;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[j][k], hc9B[n9B][4]) && hcsp5c[j][k] != hc9B[n9B][7] && hcsp5c[j][k] !=
                                                                                                       hc9B[n9B][4]) {
                        hc9B[n9B][0]= hcsp5c[j][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[j][k], hc9B[n9B][5]) && hcsp5c[j][k] != hc9B[n9B][7] && hcsp5c[j][k] !=
                                                                                                       hc9B[n9B][5]) {
                        hc9B[n9B][1]= hcsp5c[j][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k], hc9B[n9B][4]) && hcsp5c[i][k] != hc9B[n9B][6] && hcsp5c[i][k] !=
                                                                                                       hc9B[n9B][4]) {
                        hc9B[n9B][2]= hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k], hc9B[n9B][5]) && hcsp5c[i][k] != hc9B[n9B][6] && hcsp5c[i][k] !=
                                                                                                       hc9B[n9B][5]) {
                        hc9B[n9B][3]= hcsp5c[i][k];
                    }
                }
            }
            hc9B[n9B][8]=sp1;
            Cluster_Write_9B();

            if (do10B == 1) Clusters_Get10B(i, j);
            if (do11B == 1) {
                if (Clusters_Get11B()) {
                    s11B[hc9B[n9B][8]] = 'S';
                    ++n11B;
                }
            }
            if (do11E == 1) Clusters_Get11E_12D(i, j, sp1, sp2i, sp2j);

            ++n9B;
        }
    }
}

void Cluster_Write_9B() {
    // hc9B key: (SP5_lowerd_to_4, SP5_lowerd_to_5, SP5_higherd_to_4, SP5_higherd_to_5, SP5_i_j_com_lower, SP5_i_j_com_higher, sp5c_d1_lower, sp5c_d2_higher, s_com)
    int i;
    for(i=0; i<6; i++) {
        if (s9B[hc9B[n9B][i]] == 'C') s9B[hc9B[n9B][i]] = 'B';
    }
    if (s9B[hc9B[n9B][6]] != 'S') s9B[hc9B[n9B][6]] = 'O';
    if (s9B[hc9B[n9B][7]] != 'S') s9B[hc9B[n9B][7]] = 'O';
    s9B[hc9B[n9B][8]] = 'S';
}