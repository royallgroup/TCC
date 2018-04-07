#include "13B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get13B() {   // Detect 13B D5h clusters, i.e. twisted icosahedra
    int cp;
    int i, j, k, l, m;
    int flg;
    int clusSize=13;

    cp=-1;

    for(i=0; i < nsp5c - 1; ++i){
        for(j=i+1; j < nsp5c; ++j){
            flg = 0;
            if (hcsp5c[i][5] == hcsp5c[j][5] && hcsp5c[i][6] != hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][6]) == 1) continue;
                cp = hcsp5c[i][5];
                flg = 1;
            }
            if (hcsp5c[i][6] == hcsp5c[j][6] && hcsp5c[i][5] != hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][5], hcsp5c[j][5]) == 1) continue;
                cp = hcsp5c[i][6];
                flg = 1;
            }
            if (hcsp5c[i][5] == hcsp5c[j][6] && hcsp5c[i][6] != hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][5]) == 1) continue;
                cp = hcsp5c[i][5];
                flg = 1;
            }
            if (hcsp5c[i][6] == hcsp5c[j][5] && hcsp5c[i][5] != hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][5], hcsp5c[j][6]) == 1) continue;
                cp = hcsp5c[i][6];
                flg = 1;
            }
            if (flg==1) {
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if(hcsp5c[i][k] == hcsp5c[j][l]) break;
                    }
                    if(l<5) break;
                }
                if (k<5) continue;
                for (k=0; k<5; ++k) {
                    m = 0;
                    for (l=0; l<5; ++l) {
                        if (Bonds_BondCheck(hcsp5c[i][k], hcsp5c[j][l]) == 1) ++m;
                        if (m == 2) break;
                    }
                    if(m != 1) break;
                }
                if (k<5) continue;
                // ERROR: need to check converse, i.e. every SP5 from 7A_j bonded to one from SP5 from 7A_i

                if (n13B == m13B) {
                    hc13B= resize_2D_int(hc13B, m13B, m13B + incrStatic, clusSize, -1);
                    m13B= m13B + incrStatic;
                }

                hc13B[n13B][0] = cp;
                k = 1;
                if(hcsp5c[i][5] != cp) hc13B[n13B][k++] = hcsp5c[i][5];
                if(hcsp5c[i][6] != cp) hc13B[n13B][k++] = hcsp5c[i][6];
                if(hcsp5c[j][5] != cp) hc13B[n13B][k++] = hcsp5c[j][5];
                if(hcsp5c[j][6] != cp) hc13B[n13B][k++] = hcsp5c[j][6];
                for(l=0; l<5; ++l) hc13B[n13B][k++] = hcsp5c[i][l];
                for(l=0; l<5; ++l) hc13B[n13B][k++] = hcsp5c[j][l];
                quickSort(&hc13B[n13B][1], 2);
                quickSort(&hc13B[n13B][3], 10);

                Cluster_Write_13B();
            }
        }
    }
}

void Cluster_Write_13B() {
    int i;
    s13B[hc13B[n13B][0]] = 'S';
    if(s13B[hc13B[n13B][1]] != 'S') s13B[hc13B[n13B][1]] = 'O';
    if(s13B[hc13B[n13B][2]] != 'S') s13B[hc13B[n13B][2]] = 'O';
    for(i=3; i<13; i++) {
        if (s13B[hc13B[n13B][i]] == 'C') s13B[hc13B[n13B][i]] = 'B';
    }

    ++n13B;
}