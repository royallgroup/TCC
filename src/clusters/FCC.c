#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include "FCC.h"

//!  An FCC cluster is a 13 particle cluster of the FCC lattice, it is made from either four sp3b clusters or three sp3b clusters and a 5A cluster
/*!
*  Find FCC clusters
*  Method 1: An FCC cluster is constructed from four sp3b clusters
*  Method 2: An FCC cluster is constructed from three sp3b clusters and a 5A cluster
*
*  Cluster output: SBBOOBBOBBOOO
*  Storage order: unknown
*
*/
void Clusters_GetFCC() {
    int i, j, j2, k, l, m, n;
    int i1, i2, i3;
    int cp, bpi, bpj, nbpi, nbpj;
    int flg1, flg2, flg3;
    int clusSize=13;

    cp=bpi=bpj=nbpi=nbpj=i3=-1;


    for (i=0; i < nsp3b - 2; ++i) { // loop over all sp3b_i
        for (j2=0; j2<3; j2++) {
            for (j=0; j < nmem_sp3b[hcsp3b[i][j2]] - 1; ++j) { // loop over all sp3b_j
                if (mem_sp3b[hcsp3b[i][j2]][j] <= i) continue;
                if (Bonds_BondCheck(hcsp3b[i][3], hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3]) == 0) continue; // spindle spindle bond
                m = n = 0;
                for (k=0; k<4; ++k) {
                    for (l=0; l<4; ++l) {
                        if (hcsp3b[i][k] == hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l]) {
                            if (k == 3 || l == 3) m +=2;
                            cp = hcsp3b[i][k];
                            ++m;
                        }
                    }
                }
                if(m != 1) continue; // one common particle between sp3b_i and sp3b_j

                for (k=0; k < nFCC; ++k) {
                    if (hcFCC[k][0] == cp) break;   // check for other degenerate FCC clusters which have cp here
                }
                if (k < nFCC) continue;  // found this fcc cluster before

                for (k=0; k<3; ++k) {
                    if (hcsp3b[i][k] == cp) continue;
                    for (l=0; l<3; ++l) {
                        if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l] == cp) continue;
                        if (Bonds_BondCheck(hcsp3b[i][k], hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l]) == 1) {
                            bpi = hcsp3b[i][k];
                            bpj = hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l];
                            ++n;
                        }
                    }
                }
                if(n != 1) continue;    // one bond between one pair of atoms from SP3_i and SP3_j

                // we may have an FCC cluster, build it into hcFCC then overwrite it later if it aint
                if (nFCC == mFCC) {
                    hcFCC= resize_2D_int(hcFCC, mFCC, mFCC + incrStatic, clusSize, -1);
                    mFCC= mFCC + incrStatic;
                }

                for (k=0; k<3; ++k) {
                    if (hcsp3b[i][k] != cp && hcsp3b[i][k] != bpi) nbpi= hcsp3b[i][k];
                    if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][k] != cp && hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][k] != bpj) nbpj= hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][k];
                }
                hcFCC[nFCC][0] = cp;
                hcFCC[nFCC][1] = nbpi;
                hcFCC[nFCC][2] = bpi;
                hcFCC[nFCC][3] = hcsp3b[i][3];
                hcFCC[nFCC][4] = hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3];
                hcFCC[nFCC][5] = nbpj;
                hcFCC[nFCC][6] = bpj;   // store first two clusters

                for (k=j+1; k < nmem_sp3b[hcsp3b[i][j2]]; ++k) { // loop over all sp3b_k
                    if (mem_sp3b[hcsp3b[i][j2]][k] <= i) continue;
                    for (l=0; l < nFCC; ++l) {
                        if (hcFCC[l][0] == cp) break;   // check for other degenerate FCC clusters which have cp here
                    }
                    if (l < nFCC) break; // found this fcc cluster before

                    if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3] == cp) continue;    // check sp3b_k spindle isnt common particle
                    if (!(Bonds_BondCheck(hcsp3b[i][3], hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3]) == 1 && Bonds_BondCheck(
                            hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3], hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3]) == 1)) continue;  // check sp3b_k spindle is bonded to sp3b_i sp3b_j spindles
                    i1=i2=-1;
                    for (l=0; l<3; ++l) {
                        if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][l] == cp) {
                            i1 = l; // identity of common particle in SP3_k
                        }
                        else {
                            i2 = hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][l]; // identity of a loose particle in SP3_k
                            flg1 = i2 == nbpi || i2 == bpi || i2 == nbpj || i2 == bpj;
                            flg1 = flg1 || i2 == hcsp3b[i][3] || i2 == hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3];
                            if(flg1==1) break;
                        }
                    }
                    if (l<3 || i1==-1 || i2==-1) continue;  // SP3_k must have 2 new particles

                    m = 0;
                    flg1 = flg2 = 0;
                    for (l=0; l<3; ++l) {
                        if (l == i1) continue;
                        n = hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][l];
                        if (Bonds_BondCheck(n, nbpi) == 1 && Bonds_BondCheck(n, nbpj) == 0) {
                            i2 = n;
                            ++m;
                            flg1 =1;
                        }
                        if (Bonds_BondCheck(n, nbpj) == 1 && Bonds_BondCheck(n, nbpi) == 0) {
                            i3 = n;
                            ++m;
                            flg2 = 1;
                        }
                    }
                    if (m != 2) continue;  // SP3_k must have 2 particles bonded to the non-bonded non-common pair from SP3_i and SP3_j

                    if (!(flg1==1 && flg2==1)) continue; // now we have the 6 membered ring

                    for (l=0; l < nsp3b; ++l) { // loop over all sp3b_l
                        if (l==i || l == mem_sp3b[hcsp3b[i][j2]][j] || l == mem_sp3b[hcsp3b[i][j2]][k]) continue;   // must be new sp3b
                        if (hcsp3b[l][3] != cp) continue; // must have spindle as common particle
                        for (m=0; m<3; ++m) {
                            n = hcsp3b[l][m];
                            flg1 = n == nbpi || n == bpi || n == bpj;
                            flg1 = flg1 || n == nbpj || n == i3 || n == i2;
                            flg1 = flg1 || n == hcsp3b[i][3] || n == hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3];
                            flg1 = flg1 || n == hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3];
                            if (flg1==1) break;
                        }
                        if (m<3) continue; // the SP3_l ring particles must be new

                        flg1 = flg2 = flg3 = 0;
                        for (m=0; m<3; ++m) {
                            n = hcsp3b[l][m];
                            if (Bonds_BondCheck(n, bpi) == 1 && Bonds_BondCheck(n, bpj) == 1) {
                                if (flg1==1) break;
                                flg1 = 1;
                                hcFCC[nFCC][10] = n;
                            }
                            if (Bonds_BondCheck(n, nbpj) == 1 && Bonds_BondCheck(n, i3) == 1) {
                                if(flg2==1) break;
                                flg2 = 1;
                                hcFCC[nFCC][11] = n;
                            }
                            if (Bonds_BondCheck(n, i2) == 1 && Bonds_BondCheck(n, nbpi) == 1) {
                                if(flg3==1) break;
                                flg3 = 1;
                                hcFCC[nFCC][12] = n;
                            }
                        }
                        if(m<3) continue;  // SP3_l particles must be bonded correctly to six-membered ring
                        if(flg1==1 && flg2==1 && flg3==1) break;   // SP3_l particles must be bonded correctly to six-membered ring
                    } // end l loop
                    if (l == nsp3b) { // required SP3b cluster not found
                        for (l=0; l < nsp3c; ++l) { // loop over all sp3c_l
                            if (hcsp3c[l][3] != cp && hcsp3c[l][4] != cp) continue;  // common particle must be a spindle
                            for (m=0; m<3; ++m) {
                                n = hcsp3c[l][m];
                                flg1 = n == nbpi || n == bpi || n == bpj;
                                flg1 = flg1 || n == nbpj || n == i3 || n == i2;
                                flg1 = flg1 || n == hcsp3b[i][3] || n == hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3];
                                flg1 = flg1 || n == hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3];
                                if (flg1==1) break;
                            }
                            if (m<3) continue; // the SP3_l ring particles must be new

                            flg1 = flg2 = flg3 = 0;
                            for (m=0; m<3; ++m) {
                                n = hcsp3c[l][m];
                                if (Bonds_BondCheck(n, bpi) == 1 && Bonds_BondCheck(n, bpj) == 1) {
                                    if (flg1==1) break;
                                    flg1 = 1;
                                    hcFCC[nFCC][10] = n;
                                }
                                if (Bonds_BondCheck(n, nbpj) == 1 && Bonds_BondCheck(n, i3) == 1) {
                                    if (flg2==1) break;
                                    flg2 = 1;
                                    hcFCC[nFCC][11] = n;
                                }
                                if (Bonds_BondCheck(n, i2) == 1 && Bonds_BondCheck(n, nbpi) == 1) {
                                    if(flg3==1) break;
                                    flg3 = 1;
                                    hcFCC[nFCC][12] = n;
                                }
                            }
                            if (m<3) continue;  // SP3_l particles must be bonded correctly to six-membered ring
                            if (flg1==1 && flg2==1 && flg3==1) break;   // SP3_l particles must be bonded correctly to six-membered ring
                        }
                        if (l == nsp3c) continue;
                    } // End if statement for SP3c search loop

                    // We've now found an FCC cluster
                    if (nFCC == mFCC) { hcFCC= resize_2D_int(hcFCC, mFCC, mFCC + incrStatic, clusSize, -1);
                        mFCC= mFCC + incrStatic;
                    }

                    hcFCC[nFCC][7] = hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3];
                    hcFCC[nFCC][8] = i3;
                    hcFCC[nFCC][9] = i2;
                    quickSort(&hcFCC[nFCC][1], 12);

                    Cluster_Write_FCC();
                }
            }
        }
    }
}

void Cluster_Write_FCC() {
    sFCC[hcFCC[nFCC][0]] = 'S';
    if(sFCC[hcFCC[nFCC][1]] == 'C') sFCC[hcFCC[nFCC][1]] = 'B';
    if(sFCC[hcFCC[nFCC][2]] == 'C') sFCC[hcFCC[nFCC][2]] = 'B';
    if(sFCC[hcFCC[nFCC][3]] != 'S') sFCC[hcFCC[nFCC][3]] = 'O';
    if(sFCC[hcFCC[nFCC][4]] != 'S') sFCC[hcFCC[nFCC][4]] = 'O';
    if(sFCC[hcFCC[nFCC][5]] == 'C') sFCC[hcFCC[nFCC][5]] = 'B';
    if(sFCC[hcFCC[nFCC][6]] == 'C') sFCC[hcFCC[nFCC][6]] = 'B';
    if(sFCC[hcFCC[nFCC][7]] != 'S') sFCC[hcFCC[nFCC][7]] = 'O';
    if(sFCC[hcFCC[nFCC][8]] == 'C') sFCC[hcFCC[nFCC][8]] = 'B';
    if(sFCC[hcFCC[nFCC][9]] == 'C') sFCC[hcFCC[nFCC][9]] = 'B';
    if(sFCC[hcFCC[nFCC][10]] != 'S') sFCC[hcFCC[nFCC][10]] = 'O';
    if(sFCC[hcFCC[nFCC][11]] != 'S') sFCC[hcFCC[nFCC][11]] = 'O';
    if(sFCC[hcFCC[nFCC][12]] != 'S') sFCC[hcFCC[nFCC][12]] = 'O';

    ++nFCC;
}