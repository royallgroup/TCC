#include "7K.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get7K() {    // Detect 7K clusters from 2 5A clusters
    int i, j, j2, k, l, m;
    int scom, sother[2], sp3_com[2], sp3c_i_other, sp3c_j_other;
    int clusSize=7;

    scom=sother[0]=sother[1]=sp3_com[0]=sp3_com[1]=sp3c_i_other=sp3c_j_other=-1;

    for (i=0; i<nsp3c-1; ++i) {  // loop over all 5A_i
        for (j2=3; j2<5; ++j2) {    // loop over both spindles of 5A_i
            for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]; ++j) {  // loop over all 5A_j common with spindle of 5A_i
                if (mem_sp3c[hcsp3c[i][j2]][j]<=i) continue;  // don't find 7K twice

                m=0;
                for (k=3; k<5; k++) {
                    for (l=3; l<5; l++) {
                        if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            if (m>=1) {
                                m++;
                                break;
                            }
                            scom=hcsp3c[i][k];
                            m++;
                        }
                    }
                    if (m>=2) break;
                }
                if (m!=1) continue; // exactly one common spindle between 5A_i and 5A_j

                if (hcsp3c[i][3] == scom) sother[0]=hcsp3c[i][4];
                else sother[0]=hcsp3c[i][3];
                if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == scom) sother[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                else sother[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];

                m=0;
                for (k=0; k<5; k++) {
                    if (sother[0]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other spindle of 5A_i distinct from whole 5A_j

                m=0;
                for (k=0; k<5; k++) {
                    if (sother[1]==hcsp3c[i][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other spindle of 5A_j distinct from whole 5A_i

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<3; l++) {
                        if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            if (m>=2) {
                                m++;
                                break;
                            }
                            sp3_com[m]=hcsp3c[i][k];
                            m++;
                        }
                    }
                    if (m>=3) break;
                }
                if (m!=2) continue; // exactly two common particles in SP3 rings of 5A_i and 5A_j

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<2; l++) {
                        if (hcsp3c[i][k] == sp3_com[l]) break;
                    }
                    if (l==2) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        sp3c_i_other=hcsp3c[i][k];
                        m++;
                    }
                }
                if (m!=1) continue; // found other uncommon particle from SP3 ring of 5A_i

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<2; l++) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k] == sp3_com[l]) break;
                    }
                    if (l==2) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        sp3c_j_other=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k];
                        m++;
                    }
                }
                if (m!=1) continue; // found other uncommon particle from SP3 ring of 5A_i

                m=0;
                for (k=0; k<5; k++) {
                    if (sp3c_i_other==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other ring of 5A_i distinct from whole 5A_j

                m=0;
                for (k=0; k<5; k++) {
                    if (sp3c_j_other==hcsp3c[i][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other ring of 5A_j distinct from whole 5A_i

                if (n7K==m7K) {
                    hc7K=resize_2D_int(hc7K,m7K,m7K+incrStatic,clusSize,-1);
                    m7K=m7K+incrStatic;
                }
                // Now we have found the 7K cluster

                hc7K[n7K][0]=scom;
                hc7K[n7K][1]=sother[0];
                hc7K[n7K][2]=sother[1];
                hc7K[n7K][3]=sp3_com[0];
                hc7K[n7K][4]=sp3_com[1];
                hc7K[n7K][5]=sp3c_i_other;
                hc7K[n7K][6]=sp3c_j_other;

                quickSort(&hc7K[n7K][1],2);
                quickSort(&hc7K[n7K][3],2);
                quickSort(&hc7K[n7K][5],2);

                Cluster_Write_7K();
            }
        }
    }
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