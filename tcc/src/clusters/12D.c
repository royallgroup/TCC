#include <globals.h>
#include <tools.h>
#include "12D.h"

int Clusters_Get12D(int j, int k, int sp1, int sp2) {  // Return 1 if 12B is also 11E
    int l, m, n, o, p, q;
    int flg1, flg2;
    int break_out;
    int trial[12];
    int clusSize=12;

    if(k > j) m = k;
    else m = j;
    for (l=m+1; l < nsp5c; l++) {
        if ((hcsp5c[l][5] == sp1 && hcsp5c[l][6] == sp2) || (hcsp5c[l][6] == sp1 && hcsp5c[l][5] == sp2)) {
            flg1=flg2=0;
            p=11;
            q=0;
            for (n=0; n<5; n++) {
                if (hcsp5c[l][n] == hc11E[n11E][0]) {
                    flg1=1;
                    continue;
                }
                if (hcsp5c[l][n] == hc11E[n11E][1]) {
                    flg2=1;
                    continue;
                }
                break_out=0;
                for (o=4; o<11; o++) {
                    if (hcsp5c[l][n] == hc11E[n11E][o]) {
                        q++;
                        break_out=1;
                        break;
                    }
                }
                if (break_out==1) continue;

                if (p==12) {
                    p++;
                    break;
                }
                trial[p]= hcsp5c[l][n];
                p++;
            }
            if (flg1==0 || flg2==0 || p!=12 || q!=2) continue;

            for (n=0; n<11; n++) trial[n]= hc11E[n11E][n];
            // hc12D key: (d1_unc, d2_unc, d3_unc, d4_unc, d12_com, d13_com, d24_com, d34_com, s_d1, s_d2, s_d3, s_com)

            if(n12D == m12D) {
                hc12D= resize_2D_int(hc12D, m12D, m12D + incrStatic, clusSize, -1);
                m12D= m12D + incrStatic;
            }
            quickSort(&trial[0],4);
            quickSort(&trial[4],8);

            for(m=0; m<12; ++m) hc12D[n12D][m] = trial[m];

            Cluster_Write_12D();
            return 1;
        }
    }
    return 0;
}

void Cluster_Write_12D() {
    int i;

    for (i = 0; i < 4; i++) {
        s12D[hc12D[n12D][i]] = 'O';
    }
    for (i = 4; i < 12; i++) {
        if (s12D[hc12D[n12D][i]] == 'C') s12D[hc12D[n12D][i]] = 'B';
    }
}