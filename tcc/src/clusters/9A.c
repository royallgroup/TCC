#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include "9A.h"

void Clusters_Get9A() {    // Detect 9A D3h clusters
    int i, j, j2, k, l, m, n;
    int db[2], ob[4];
    int flg;
    int clusSize=9;

    for (i=0; i < nsp4b - 2; ++i) {  // loop over all sp4b_i
        for (j2=0; j2<1; j2++) {
            for (j=0; j < nmem_sp4b[hcsp4b[i][j2]]; ++j) { // loop over all sp4b_j
                if (mem_sp4b[hcsp4b[i][j2]][j] <= i) continue;
                if (hcsp4b[i][4] == hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4]) continue;  // if spindles common continue
                if (Bonds_BondCheck(hcsp4b[i][4], hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4])) continue;   // if spindles bonded continue
                m = 0;
                for(k=0; k<4; ++k) {
                    for(l=0; l<4; ++l) {
                        if(hcsp4b[i][k] == hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l]){
                            if(m<2) db[m] = hcsp4b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=2) continue;         // two common particles between SP4 rings of sp4b_i and sp4b_j

                m = 0;
                for (k=0; k<4; ++k) {
                    if(hcsp4b[i][k] == db[0] || hcsp4b[i][k] == db[1]) continue;  // find particles in SP4 ring of sp4b_i not common to SP4 ring of sp4b_j
                    for(l=0; l<4; ++l){
                        if(hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l] == db[0] || hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l] == db[1]) continue;    // find particles in SP4 ring of sp4b_j not common to SP4 ring of sp4b_i
                        if(Bonds_BondCheck(hcsp4b[i][k], hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l])) {    // check non-common SP4 ring particles from sp4b_i and sp4b_j are bonded
                            if(m<4) ob[m] = hcsp4b[i][k];
                            ++m;
                            if(m<4) ob[m] = hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue;
                // POSSIBLE IMPROVEMENT!! could make detection faster here by looping over sp4b clusters bonded to uncommon particles to sp4b_i
                for(k= mem_sp4b[hcsp4b[i][j2]][j] + 1; k < nsp4b; k++) {    // loop over all sp4b_k
                    // ERROR!! need to check spindle of sp4b_k distinct from spindles of sp4b_i and sp4b_j and no bonds between these three particles
                    n = 0;
                    for(l=0; l<4; ++l){
                        for(m=0; m<4; ++m){
                            if(hcsp4b[k][l] == ob[m]){
                                ++n;
                                break;
                            }
                        }
                    }
                    if (n != 4) continue;

                    // Now we have found the 9A D3h cluster
                    if (n9A == m9A) {
                        hc9A= resize_2D_int(hc9A, m9A, m9A + incrStatic, clusSize, -1);
                            m9A= m9A + incrStatic;
                        }
                        if (hcsp4b[i][4] < hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4] && hcsp4b[i][4] < hcsp4b[k][4]) {
                            hc9A[n9A][6]= hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                            hc9A[n9A][0]= hcsp4b[i][0];
                            hc9A[n9A][1]= hcsp4b[i][1];
                            hc9A[n9A][2]= hcsp4b[i][2];
                            hc9A[n9A][3]= hcsp4b[i][3];

                            if (Bonds_BondCheck(hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4], hc9A[n9A][0])) {
                                hc9A[n9A][7]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                hc9A[n9A][8]= hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                            }
                            else {
                                hc9A[n9A][7]= hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                hc9A[n9A][8]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                            }

                            for (l=0;l<4;l++) {
                                flg=1;
                                for (m=0;m<4;m++) {
                                    if (ob[l] == hc9A[n9A][m]) {
                                        flg=0;
                                        break;
                                    }
                                }
                                if (flg==1 && Bonds_BondCheck(ob[l], hc9A[n9A][0])) {
                                    hc9A[n9A][4]=ob[l];
                                }
                                else if (flg==1) {
                                    hc9A[n9A][5]=ob[l];
                                }
                            }
                        }

                        else if (hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4] < hcsp4b[i][4] && hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4] <
                                                                                         hcsp4b[k][4]) {
                            hc9A[n9A][6]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                            hc9A[n9A][0]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][0];
                            hc9A[n9A][1]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][1];
                            hc9A[n9A][2]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][2];
                            hc9A[n9A][3]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][3];

                            if (Bonds_BondCheck(hcsp4b[i][4], hc9A[n9A][0])) {
                                hc9A[n9A][7]= hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                hc9A[n9A][8]= hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                            }
                            else {
                                hc9A[n9A][7]= hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                hc9A[n9A][8]= hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                            }

                            for (l=0;l<4;l++) {
                                flg=1;
                                for (m=0;m<4;m++) {
                                    if (ob[l] == hc9A[n9A][m]) {
                                        flg=0;
                                        break;
                                    }
                                }
                                if (flg==1 && Bonds_BondCheck(ob[l], hc9A[n9A][0])) {
                                    hc9A[n9A][4]=ob[l];
                                }
                                else if (flg==1) {
                                    hc9A[n9A][5]=ob[l];
                                }
                            }
                        }

                        else {
                            hc9A[n9A][6]= hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                            hc9A[n9A][0]= hcsp4b[k][0];
                            hc9A[n9A][1]= hcsp4b[k][1];
                            hc9A[n9A][2]= hcsp4b[k][2];
                            hc9A[n9A][3]= hcsp4b[k][3];

                            if (Bonds_BondCheck(hcsp4b[i][4], hc9A[n9A][0])) {
                                hc9A[n9A][7]= hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                hc9A[n9A][8]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                            }
                            else {
                                hc9A[n9A][7]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                hc9A[n9A][8]= hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                            }

                            if (Bonds_BondCheck(db[0], hc9A[n9A][0])) {
                                hc9A[n9A][4]=db[0];
                                hc9A[n9A][5]=db[1];
                            }
                            else {
                                hc9A[n9A][4]=db[1];
                                hc9A[n9A][5]=db[0];
                            }
                        }
                    Cluster_Write_9A();
                    break;
                }
            }
        }
    }
}

void Cluster_Write_9A() {
    // hc9A key: (SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_to_0_in_SP4_lowest_s, SP4_to_1_in_SP4_lowest_s, s_lowest, s_to_0_in_SP4_lowest_s,s_to_2_in_SP4_lowest_s)
    if(s9A[hc9A[n9A][0]] == 'C') s9A[hc9A[n9A][0]] = 'B';
    if(s9A[hc9A[n9A][1]] == 'C') s9A[hc9A[n9A][1]] = 'B';
    if(s9A[hc9A[n9A][3]] == 'C') s9A[hc9A[n9A][3]] = 'B';
    if(s9A[hc9A[n9A][4]] == 'C') s9A[hc9A[n9A][4]] = 'B';
    if(s9A[hc9A[n9A][6]] == 'C') s9A[hc9A[n9A][6]] = 'B';
    if(s9A[hc9A[n9A][7]] == 'C') s9A[hc9A[n9A][7]] = 'B';
    s9A[hc9A[n9A][2]] = 'O';
    s9A[hc9A[n9A][5]] = 'O';
    s9A[hc9A[n9A][8]] = 'O';
    ++n9A;
}