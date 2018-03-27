#include "6Z.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get6Z() {    // Detect 6Z clusters from 2 5A clusters
    int flg;
    int i, j, j2, k, l;
    int cnt;
    int s1a, s2a, s1b, s2b;
    int clusSize=6;

    s1a=s2a=s1b=s2b=-1;


    for (i=0; i<nsp3c-1; ++i) {  // loop over all 5A_i
        for (j2=0; j2<1; ++j2) {
            for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]; ++j) {  // loop over all 5A_j
                if (mem_sp3c[hcsp3c[i][j2]][j]<=i) continue;
                if (hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) continue;   // no spindles of 5A_i and 5A_j are common
                if (hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) continue;

                // for each cluster one spindle particle is in sp3 ring of other 5A and other spindle isn't
                cnt = 0;    // check one 5A_i spindle are in sp3 ring of 5A_j
                flg = hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][0] || hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][1] || hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][2];
                if (flg==1){
                    s1a = hcsp3c[i][3];
                    s2a = hcsp3c[i][4];
                    ++cnt;
                }
                flg = hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][0] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][1] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][2];
                if (flg==1){
                    s1a = hcsp3c[i][4];
                    s2a = hcsp3c[i][3];
                    ++cnt;
                }
                if (cnt != 1) continue;
                cnt = 0;    // check one 5A_j spindle are in sp3 ring of 5A_i
                flg = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[i][0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[i][1] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[i][2];
                if (flg==1){
                    s1b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                    s2b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                    ++cnt;
                }
                flg = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[i][0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[i][1] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[i][2];
                if (flg==1){
                    s1b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                    s2b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                    ++cnt;
                }
                if (cnt != 1) continue;
                if (!Bonds_BondCheck(s1a, s1b)) continue;   // check spindles of i and j in sp3 ring of j and i respectively are bonded

                cnt = 0;    // check 2 particles in the sp3 rings of 5A_i and 5A_j are common
                for (k=0; k<3; ++k){
                    for (l=0; l<3; ++l){
                        if(hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]){
                            ++cnt;
                            break;
                        }
                    }
                }
                if (cnt != 2) continue;

                if (n6Z==m6Z) {
                    hc6Z=resize_2D_int(hc6Z,m6Z,m6Z+incrStatic,clusSize,-1);
                    m6Z=m6Z+incrStatic;
                }
                // Now we have found the 6Z cluster
                if (s1a<s1b) {
                    hc6Z[n6Z][0]=s1a;    // insert cluster
                    hc6Z[n6Z][1]=s1b;
                    hc6Z[n6Z][2]=s2a;
                    hc6Z[n6Z][3]=s2b;
                }
                else {
                    hc6Z[n6Z][0]=s1b;    // insert cluster
                    hc6Z[n6Z][1]=s1a;
                    hc6Z[n6Z][2]=s2b;
                    hc6Z[n6Z][3]=s2a;
                }
                cnt=4;
                for (k=0; k<3; ++k) {
                    flg=1;
                    for (l=0; l<4; ++l){
                        if (hcsp3c[i][k]==hc6Z[n6Z][l]) {
                            flg=0;
                            break;
                        }
                    }
                    if (flg==1) { hc6Z[n6Z][cnt]=hcsp3c[i][k]; cnt++; }
                }
                if (hc6Z[n6Z][5]<hc6Z[n6Z][4]) {
                    k=hc6Z[n6Z][5];
                    hc6Z[n6Z][5]=hc6Z[n6Z][4];
                    hc6Z[n6Z][4]=k;
                }
                Cluster_Write_6Z();
            }
        }
    }
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