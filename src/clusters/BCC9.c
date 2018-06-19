#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include "BCC9.h"

void Clusters_GetBCC_9() {

    //!  A BCC_9 cluster is the intersection of sp4b and sp4c clusters
    /*!
   *  Find BCC_9 clusters
   *  A BCC_9 is constructed from some combination of sp4b and sp4c clusters
   *      -
   *      -
   *      -
   *      -
   *
   *  Cluster output: SBBBBBBBB
   *  Storage order: unknown
   *
   */

    int i, j, j2, k, l, m;
    int flg;
    int s_com=-1;
    int trial[9];
    int clusSize=9;

    for (i=0; i < nsp4b - 1; i++) {
        for (j2=4; j2<5; j2++) {
            for (j=0; j < nmem_sp4b[hcsp4b[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4b[hcsp4b[i][j2]][j] <= i) continue;
                if (hcsp4b[i][4] != hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4]) continue;
                s_com= hcsp4b[i][4];

                flg=0;
                for (k=0; k<4; k++) {
                    for (l=0; l<4; l++) {
                        if (hcsp4b[i][k] == hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4b[i][k], hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][k], hcsp4b[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]= hcsp4b[i][k];
                    trial[k+5]= hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k < nBCC_9; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l] != hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9 == mBCC_9) {
                    hcBCC_9= resize_2D_int(hcBCC_9, mBCC_9, mBCC_9 + incrStatic, clusSize, -1);
                    mBCC_9= mBCC_9 + incrStatic;
                }

                for (k=0; k<9; ++k) hcBCC_9[nBCC_9][k]=trial[k];

                Cluster_Write_BCC9();
            }
        }
    }

    for (i=0; i < nsp4c - 1; i++) {
        for (j2=4; j2<6; j2++) {
            for (j=0; j < nmem_sp4c[hcsp4c[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4c[hcsp4c[i][j2]][j] <= i) continue;
                m=0;
                for (k=4; k<6; k++) {
                    for (l=4; l<6; l++) {
                        if (hcsp4c[i][k] == hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]) {
                            s_com= hcsp4c[i][k];
                            m++;
                        }
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<6; k++) {
                    if (hcsp4c[i][k] == s_com) continue;
                    for (l=0; l<6; l++) {
                        if (hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l] == s_com) continue;
                        if (hcsp4c[i][k] == hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[i][k], hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][k], hcsp4c[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]= hcsp4c[i][k];
                    trial[k+5]= hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k < nBCC_9; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l] != hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9 == mBCC_9) {
                    hcBCC_9= resize_2D_int(hcBCC_9, mBCC_9, mBCC_9 + incrStatic, clusSize, -1);
                    mBCC_9= mBCC_9 + incrStatic;
                }
                for (k=0; k<9; ++k) hcBCC_9[nBCC_9][k]=trial[k];

                Cluster_Write_BCC9();
            }
        }
    }

    for (i=0; i < nsp4b; i++) {
        for (j2=4; j2<5; j2++) {
            for (j=0; j < nmem_sp4c[hcsp4b[i][j2]]; ++j) { // loop over all sp3c_j
                m=0;
                for (k=4; k<6; k++) {
                    if (hcsp4b[i][4] == hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][k]) {
                        s_com= hcsp4b[i][4];
                        m++;
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    for (l=0; l<6; l++) {
                        if (hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][l] == s_com) continue;
                        if (hcsp4b[i][k] == hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4b[i][k], hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][k], hcsp4b[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]= hcsp4b[i][k];
                    trial[k+5]= hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k < nBCC_9; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l] != hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9 == mBCC_9) {
                    hcBCC_9= resize_2D_int(hcBCC_9, mBCC_9, mBCC_9 + incrStatic, clusSize, -1);
                    mBCC_9= mBCC_9 + incrStatic;
                }
                for (k=0; k<9; ++k) hcBCC_9[nBCC_9][k]=trial[k];

                Cluster_Write_BCC9();

            }
        }
    }
}

void Cluster_Write_BCC9() {
    int i;
    sBCC_9[hcBCC_9[nBCC_9][0]] = 'S';
    for (i = 1; i< 9; i++){
        if (sBCC_9[hcBCC_9[nBCC_9][i]] == 'C') sBCC_9[hcBCC_9[nBCC_9][i]] = 'B';
    }

    ++nBCC_9;
}