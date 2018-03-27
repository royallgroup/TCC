#include "BCC15.h"
#include "globals.h"
#include "tools.h"

void Clusters_GetBCC_15() {    // Detect 15 particle BCC clusters
    int i,j,k,l,m;
    int no_sp4cs,noSP4s;
    int sj[5];
    int clusSize=15;

    for (i=0; i < nsp4c; i++) {
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15 == mBCC_15) {
            hcBCC_15= resize_2D_int(hcBCC_15, mBCC_15, mBCC_15 + incrStatic, clusSize, -1);
            mBCC_15= mBCC_15 + incrStatic;
        }
        for (j=0; j < nBCC_15; j++) if (hcsp4c[i][4] == hcBCC_15[j][0]) break;
        if (j == nBCC_15) {
            hcBCC_15[nBCC_15][0]= hcsp4c[i][4];
            hcBCC_15[nBCC_15][1]= hcsp4c[i][5];
            for (j=0; j<4; j++) hcBCC_15[nBCC_15][j + 7]= hcsp4c[i][j];

            no_sp4cs=0;
            noSP4s=4;
            sj[0]=sj[1]=sj[2]=sj[3]=sj[4]=-1;
            for (j=i+1; j < nsp4c; j++) {
                if (hcsp4c[j][4] == hcBCC_15[nBCC_15][0] || hcsp4c[j][5] == hcBCC_15[nBCC_15][0]) {
                    m=0;
                    for (l=1;l<2+no_sp4cs;l++) {
                        if (hcsp4c[j][4] == hcBCC_15[nBCC_15][l] || hcsp4c[j][5] == hcBCC_15[nBCC_15][l]) m++;
                    }
                    if (m==0) {
                        for (l=7;l<7+noSP4s;l++) {
                            if (hcsp4c[j][4] == hcBCC_15[nBCC_15][l] || hcsp4c[j][5] == hcBCC_15[nBCC_15][l]) m++;
                        }
                        if (m==0) {
                            for (k=0; k<4; k++) {
                                for (l=0;l<2+no_sp4cs;l++) {
                                    if (hcsp4c[j][k] == hcBCC_15[nBCC_15][l]) m++;
                                }
                            }
                            if (m==0) {
                                if (hcsp4c[j][4] == hcBCC_15[nBCC_15][0]) hcBCC_15[nBCC_15][2 + no_sp4cs]= hcsp4c[j][5];
                                else hcBCC_15[nBCC_15][2 + no_sp4cs]= hcsp4c[j][4];
                                for (k=0; k<4; k++) {
                                    m=0;
                                    for (l=7;l<7+noSP4s;l++) {
                                        if (hcsp4c[j][k] == hcBCC_15[nBCC_15][l]) m++;
                                    }
                                    if (m==0) {
                                        if (noSP4s>=8) {
                                            noSP4s++;
                                            break;
                                        }
                                        hcBCC_15[nBCC_15][7 + noSP4s]= hcsp4c[j][k];
                                        noSP4s++;
                                    }
                                }
                                sj[no_sp4cs]=j;
                                no_sp4cs++;
                            }
                        }
                    }
                }
            }
            if (no_sp4cs==5 && noSP4s==8) {
                // We've now found an BCC_15 cluster
                Cluster_Write_BCC_15(clusSize);
            }
        }
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15 == mBCC_15) {
            hcBCC_15= resize_2D_int(hcBCC_15, mBCC_15, mBCC_15 + incrStatic, clusSize, -1);
            mBCC_15= mBCC_15 + incrStatic;
        }
        for (j=0; j < nBCC_15; j++) if (hcsp4c[i][5] == hcBCC_15[j][0]) break;
        if (j < nBCC_15) continue;

        hcBCC_15[nBCC_15][0]= hcsp4c[i][5];
        hcBCC_15[nBCC_15][1]= hcsp4c[i][4];
        for (j=0; j<4; j++) hcBCC_15[nBCC_15][j + 7]= hcsp4c[i][j];

        no_sp4cs=0;
        noSP4s=4;
        sj[0]=sj[1]=sj[2]=sj[3]=sj[4]=-1;
        for (j=i+1; j < nsp4c; j++) {
            if (hcsp4c[j][4] == hcBCC_15[nBCC_15][0] || hcsp4c[j][5] == hcBCC_15[nBCC_15][0]) {
                m=0;
                for (l=1;l<2+no_sp4cs;l++) {
                    if (hcsp4c[j][4] == hcBCC_15[nBCC_15][l] || hcsp4c[j][5] == hcBCC_15[nBCC_15][l]) m++;
                }
                if (m==0) {
                    for (l=7;l<7+noSP4s;l++) {
                        if (hcsp4c[j][4] == hcBCC_15[nBCC_15][l] || hcsp4c[j][5] == hcBCC_15[nBCC_15][l]) m++;
                    }
                    if (m==0) {
                        for (k=0; k<4; k++) {
                            for (l=0;l<2+no_sp4cs;l++) {
                                if (hcsp4c[j][k] == hcBCC_15[nBCC_15][l]) m++;
                            }
                        }
                        if (m==0) {
                            if (hcsp4c[j][4] == hcBCC_15[nBCC_15][0]) hcBCC_15[nBCC_15][2 + no_sp4cs]= hcsp4c[j][5];
                            else hcBCC_15[nBCC_15][2 + no_sp4cs]= hcsp4c[j][4];
                            for (k=0; k<4; k++) {
                                m=0;
                                for (l=7;l<7+noSP4s;l++) {
                                    if (hcsp4c[j][k] == hcBCC_15[nBCC_15][l]) m++;
                                }
                                if (m==0) {
                                    if (noSP4s>=8) {
                                        noSP4s++;
                                        break;
                                    }
                                    hcBCC_15[nBCC_15][7 + noSP4s]= hcsp4c[j][k];
                                    noSP4s++;
                                }
                            }
                            sj[no_sp4cs]=j;
                            no_sp4cs++;
                        }
                    }
                }
            }
        }

        if (no_sp4cs==5 && noSP4s==8) {
            // We've now found an BCC_15 cluster
            Cluster_Write_BCC_15(clusSize);
        }

    }
}

void Cluster_Write_BCC_15(int clusSize) {
    int i;
    if (nBCC_15 == mBCC_15) {
        hcBCC_15 = resize_2D_int(hcBCC_15, mBCC_15, mBCC_15 + incrStatic, clusSize, -1);
        mBCC_15 = mBCC_15 + incrStatic;
    }

    quickSort(&hcBCC_15[nBCC_15][1], 6);
    quickSort(&hcBCC_15[nBCC_15][7], 8);

    sBCC_15[hcBCC_15[nBCC_15][0]] = 'S';

    for (i = 1; i < 7; i++) {
        if(sBCC_15[hcBCC_15[nBCC_15][i]] != 'S') sBCC_15[hcBCC_15[nBCC_15][i]] = 'O';
    }

    for (i = 7; i < 15; i++) {
        if (sBCC_15[hcBCC_15[nBCC_15][i]] == 'C') sBCC_15[hcBCC_15[nBCC_15][i]] = 'B';
    }
    ++nBCC_15;
}