#include <globals.h>
#include <tools.h>
#include <bonds.h>
#include "10A.h"

void Clusters_Get10A() { // Detect 10A D4d clusters
    int i, j, j2, k, l, m;
    char errMsg[1000];
    int clusSize=10;
    int *used_sp4b;

    used_sp4b= malloc(nsp4b * sizeof(int)); if (used_sp4b == NULL) { sprintf(errMsg, "Clusters_Get10A(): used_sp4b[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i < nsp4b; ++i) used_sp4b[i] = 0;

    for (i=0; i < nsp4b - 1; ++i) {  // loop over all sp4b_i
        for (j2=0; j2 < nsp4b; ++j2) used_sp4b[j2] = 0;
        used_sp4b[i]=1;
        for (j2=0; j2 < num_bonds[hcsp4b[i][0]]; ++j2) {
            for (j=0; j < nmem_sp4b[bNums[hcsp4b[i][0]][j2]]; ++j) {    // loop over sp4b_j
                if (mem_sp4b[bNums[hcsp4b[i][0]][j2]][j] <= i) continue;
                if (used_sp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]] == 1) continue;
                used_sp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]]=1;

                if(hcsp4b[i][4] == hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][4]) continue;
                if(Bonds_BondCheck(hcsp4b[i][4], hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][4])) continue;
                for(k=0; k<5; ++k){
                    for(l=0; l<5; ++l){
                        if(hcsp4b[i][k] == hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][l]) break;
                    }
                    if(l<5) break;
                }
                if(k<5) continue;
                for(k=0; k<4; ++k){
                    m = 0;
                    for(l=0; l<4; ++l) if(Bonds_BondCheck(hcsp4b[i][k], hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][l])) ++m;
                    if(m!=2) break;
                }
                // ERROR: Need to check converse, i. e. each SP4 ring particle from sp4b_j bonded to to exactly two particles from sp4b_i
                if(k==4) {
                    if(n10A == m10A) {
                        hc10A= resize_2D_int(hc10A, m10A, m10A + incrStatic, clusSize, -1);
                        m10A= m10A + incrStatic;
                    }

                    // hc10A key: (SP4s going up, spindles going up)

                    hc10A[n10A][0]= hcsp4b[i][0];
                    hc10A[n10A][1]= hcsp4b[i][1];
                    hc10A[n10A][2]= hcsp4b[i][2];
                    hc10A[n10A][3]= hcsp4b[i][3];
                    hc10A[n10A][4]= hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][0];
                    hc10A[n10A][5]= hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][1];
                    hc10A[n10A][6]= hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][2];
                    hc10A[n10A][7]= hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][3];
                    quickSort(&hc10A[n10A][0], 8);
                    hc10A[n10A][8]= hcsp4b[i][4];
                    hc10A[n10A][9]= hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][4];
                    quickSort(&hc10A[n10A][8], 2);

                    Cluster_Write_10A();
                }
            }
        }
    }
    free(used_sp4b);
}

void Cluster_Write_10A() {
    int i;

    for(i=0; i<8; i++) {
        if (s10A[hc10A[n10A][i]] == 'C') s10A[hc10A[n10A][i]] = 'B';
    }
    s10A[hc10A[n10A][8]] = 'O';
    s10A[hc10A[n10A][9]] = 'O';

    ++n10A;
}