#include "clusters.h"
#include "8K.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get8K() {    // Detect 8K clusters
    int i, j, j2, k, l, m, n;
    int cp[2], unc[3], scom, sother[2];
    int clusSize=8;

    cp[0]=cp[1]=unc[0]=unc[1]=unc[2]=scom=sother[0]=sother[1]=-1;

    for (i=0; i < nsp3c - 2; ++i) {  // loop over all sp3c_i
        for (j2=0; j2<3; j2++) {    // loop over all particles in SP3 ring of sp3c_i
            for (j=0; j < nmem_sp3c[hcsp3c[i][j2]] - 1; ++j) { // loop over all sp3c_j which sp3c[i][j2] is a member of
                if (mem_sp3c[hcsp3c[i][j2]][j] <= i) continue;  // check not used mem_sp3c[sp3c[i][j2]][j] before

                m = 0;  // check j2 from SP3 ring of sp3c_i is in SP3 ring of mem_sp3c[sp3c[i][j2]][j]
                for(k=0; k<3; ++k) {
                    if (hcsp3c[i][j2] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k]) {
                        cp[0]= hcsp3c[i][j2];
                        m++;
                    }
                }
                if (m!=1) continue;

                m = 0;  // find extra 1 particle from SP3 ring of sp3c_i which is also in SP3 ring of mem_sp3c[sp3c[i][j2]][j]
                for(k=0; k<3; ++k) {
                    if (k==j2) continue;    // don't check j2 again
                    if (hcsp3c[i][k] < hcsp3c[i][j2]) continue; // will have found before or do not find after when using different j2 particle
                    for(l=0; l<3; ++l) {
                        if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            cp[1]= hcsp3c[i][k];
                            m++;
                        }
                    }
                }
                if (m!=1) continue;

                m = 0;  // check exactly one common spindle between sp3c_i and mem_sp3c[sp3c[i][j2]][j]
                for(k=3; k<5; ++k) {
                    for(l=3; l<5; ++l) {
                        if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            scom= hcsp3c[i][k];
                            m++;
                        }
                    }
                }
                if (m!=1) continue;

                if (hcsp3c[i][3] == scom) sother[0]= hcsp3c[i][4];
                else sother[0]= hcsp3c[i][3];
                if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == scom) sother[1]= hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                else sother[1]= hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];

                k=0;    // find particle from SP3 ring of sp3c_i which is common with 5A mem_sp3c[sp3c[i][j2]][j]
                for(l=0; l<3; ++l) {
                    if (hcsp3c[i][l] == cp[0] || hcsp3c[i][l] == cp[1]) continue;
                    for(m=0; m<5; ++m) {
                        if (hcsp3c[i][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][m]) break;
                    }
                    if (m==5) {
                        if (k==1) {
                            k++;
                            break;
                        }
                        unc[0]= hcsp3c[i][l];
                        k++;
                    }
                }
                if (k!=1) continue;

                k=0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                for(l=0; l<3; ++l) {
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == cp[0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == cp[1]) continue;
                    for(m=0; m<5; ++m) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == hcsp3c[i][m]) break;
                    }
                    if (m==5) {
                        if (k==1) {
                            k++;
                            break;
                        }
                        unc[1]= hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                        k++;
                    }
                }
                if (k!=1) continue;

                for (k=j+1; k < nmem_sp3c[hcsp3c[i][j2]]; ++k) {
                    if (mem_sp3c[hcsp3c[i][j2]][k] <= i) continue;  // higher index of sp3c cluster than i
                    if (mem_sp3c[hcsp3c[i][j2]][k] <= mem_sp3c[hcsp3c[i][j2]][j]) continue;   // higher index of sp3c cluster than mem_sp3c[sp3c[i][j2]][j]

                    n = 0;  // check common SP3 ring particles are exactly cp
                    for(l=0; l<3; ++l) {
                        for(m=0; m<2; ++m) {
                            if (cp[m] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]) {
                                n++;
                            }
                        }
                    }
                    if (n!=2) continue;

                    n = 0;  // check spindles are exactly sother
                    for(l=3; l<5; ++l) {
                        for(m=0; m<2; ++m) {
                            if (sother[m] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]) {
                                n++;
                            }
                        }
                    }
                    if (n!=2) continue;

                    n=0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                    for(l=0; l<3; ++l) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] == cp[0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] == cp[1]) continue;
                        for(m=0; m<5; ++m) {
                            if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] == hcsp3c[i][m] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] ==
                                                                                         hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) break;
                        }
                        if (m==5) {
                            if (n==1) {
                                n++;
                                break;
                            }
                            unc[2]= hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                            n++;
                        }
                    }
                    if (n!=1) continue;

                    // Now we have found the 8K cluster
                    if (n8K == m8K) {
                        hc8K= resize_2D_int(hc8K, m8K, m8K + incrStatic, clusSize, -1);
                        m8K= m8K + incrStatic;
                    }
                    // hc8K key: (SP3_common_1, SP3_common_2, spindle_1, spindle_2, spindle_3, other_SP3_1, other_SP3_2, other_SP3_3)
                    hc8K[n8K][0]=cp[0];
                    hc8K[n8K][1]=cp[1];
                    hc8K[n8K][2]=scom;
                    hc8K[n8K][3]=sother[0];
                    hc8K[n8K][4]=sother[1];
                    hc8K[n8K][5]=unc[0];
                    hc8K[n8K][6]=unc[1];
                    hc8K[n8K][7]=unc[2];

                    quickSort(&hc8K[n8K][0], 2);
                    quickSort(&hc8K[n8K][2], 3);
                    quickSort(&hc8K[n8K][5], 3);

                    Cluster_Write_8K();
                }
            }
        }
    }
}

void Cluster_Write_8K() {
    if(s8K[hc8K[n8K][2]] == 'C') s8K[hc8K[n8K][2]] = 'B';
    if(s8K[hc8K[n8K][3]] == 'C') s8K[hc8K[n8K][3]] = 'B';
    if(s8K[hc8K[n8K][4]] == 'C') s8K[hc8K[n8K][4]] = 'B';
    if(s8K[hc8K[n8K][5]] == 'C') s8K[hc8K[n8K][5]] = 'B';
    if(s8K[hc8K[n8K][6]] == 'C') s8K[hc8K[n8K][6]] = 'B';
    if(s8K[hc8K[n8K][7]] == 'C') s8K[hc8K[n8K][7]] = 'B';
    s8K[hc8K[n8K][0]] = 'O';
    s8K[hc8K[n8K][1]] = 'O';
    ++n8K;
}