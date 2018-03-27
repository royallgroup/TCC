#include "8B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get8B() { // Detect 8B Cs clusters
    int i;
    int clusSize=8;

    for (i=0; i<nsp5c; ++i) {    // loop over all 7A_i
        Clusters_8B_loop(i, clusSize, hcsp5c[i][5], hcsp5c[i][6]);
        Clusters_8B_loop(i, clusSize, hcsp5c[i][6], hcsp5c[i][5]);
    }
}

void Clusters_8B_loop(int i, int clusSize, int primary_spindle, int secondary_spindle) {

    int j, k, l, m;
    int n1, nbs, unc[3];
    int break_out;

    for (j=0; j < num_bonds[primary_spindle]; ++j) { // loop over all j particles bonded to first spindle of 7A_i
        n1 = bNums[primary_spindle][j];
        for (k=0; k<5; ++k) if (n1 == hcsp5c[i][k]) break;
        if (k<5) continue;
        if (n1 == secondary_spindle) continue; // now is n1 bonded to sp5
        nbs = 0; // number of bonds
        for (k=0; k<5; ++k) if (Bonds_BondCheck(n1, hcsp5c[i][k])) ++nbs;
        if (nbs != 2) continue;

        // Now we have found the 8B Cs cluster
        if (n8B==m8B) {
            hc8B=resize_2D_int(hc8B,m8B,m8B+incrStatic,clusSize,-1);
            m8B=m8B+incrStatic;
        }
        l=0;
        m=3;
        break_out=0;
        for (k=0; k<5; ++k) {
            if (Bonds_BondCheck(n1, hcsp5c[i][k])) {
                if (m==5) {
                    break_out=1;
                    break;
                }
                hc8B[n8B][m]=hcsp5c[i][k];
                m++;
            }
            else {
                if (l==3) {
                    break_out=1;
                    break;
                }
                unc[l]=hcsp5c[i][k];
                l++;
            }
        }
        if (break_out==1 || m<5 || l<3) continue;

        quickSort(&hc8B[n8B][3],2);
        for (k=0;k<3;k++) hc8B[n8B][k]=unc[k];
        quickSort(&hc8B[n8B][0],3);

        hc8B[n8B][5]=secondary_spindle;
        hc8B[n8B][6]=primary_spindle;
        hc8B[n8B][7]=n1;


        Cluster_Write_8B();
    }
}

void Cluster_Write_8B() {
    // hc8B key: (SP5_to_4, SP5_to_0/2, SP5_to_3, SP5_to_n1(lower), SP5_to_n1(greater), s, s_to_n1, n1)
    if (s8B[hc8B[n8B][7]] == 'C') s8B[hc8B[n8B][7]] = 'B';
    if (s8B[hc8B[n8B][0]] == 'C') s8B[hc8B[n8B][0]] = 'B';
    if (s8B[hc8B[n8B][1]] == 'C') s8B[hc8B[n8B][1]] = 'B';
    if (s8B[hc8B[n8B][2]] == 'C') s8B[hc8B[n8B][2]] = 'B';
    if (s8B[hc8B[n8B][3]] == 'C') s8B[hc8B[n8B][3]] = 'B';
    if (s8B[hc8B[n8B][4]] == 'C') s8B[hc8B[n8B][4]] = 'B';
    s8B[hc8B[n8B][5]] = 'O';
    s8B[hc8B[n8B][6]] = 'O';
    ++n8B;
}