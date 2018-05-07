#include "HCP.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_GetHCP() {

    //!  An HCP cluster is a 13 particle cluster of the HCP lattice, it is made from three 5A clusters
    /*!
   *  Find HCP clusters
   *  An HCP cluster is made from three 5A clusters where:
   *      - There is one common ring particle the three 5A clusters.
   *      - The spindles are all distinct and form two sp3 rings above and below the plane created by the ring particles.
   *      - Within the three 5A clusters, the spindle atoms are only bonded to the particles from the cluster’s own sp3 ring.
   *      - The uncommon ring particles form a six-membered ring around the common ring particle
   *
   *  Cluster output: unknown
   *  Storage order: unknown
   *
   */

    int i, j, j2, k, l, m, n;
    int ia[2], ja[2], ka[2];
    int cp, x;
    int h1i, h1j, h2i, h2j;
    int flg1, flg2, flg3;
    int clusSize=13;

    cp=-1;

    for (i=0; i < nsp3c - 2; ++i) { // loop over all sp3c_i
        for (j2=0; j2<3; j2++) {
        for (j=0; j < nmem_sp3c[hcsp3c[i][j2]] - 1; ++j) { // loop over all sp3c_j
            if (mem_sp3c[hcsp3c[i][j2]][j] <= i) continue;
            flg1 = Bonds_BondCheck(hcsp3c[i][3], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) && Bonds_BondCheck(hcsp3c[i][4],
                                                                                                           hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]);
            flg2 = Bonds_BondCheck(hcsp3c[i][3], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) && Bonds_BondCheck(hcsp3c[i][4],
                                                                                                           hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]);
            if (!(flg1==1 || flg2==1)) continue; // sp3c_i has both spindles bonded to sp3c_j spindles
            if (flg1==1) {
                h1i = hcsp3c[i][3]; // h1i bonded to h1j, h2i bonded to h2j
                h1j = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                h2i = hcsp3c[i][4];
                h2j = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
            }
            else {
                h1i = hcsp3c[i][3]; // h1i bonded to h1j, h2i bonded to h2j
                h1j = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                h2i = hcsp3c[i][4];
                h2j = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
            }

            flg1 = hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
            flg2 = hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
            flg1 = flg1 || flg2;
            if (flg1==1) continue; // spindles particles must be distinct, no common spindles

            m = 0;
            for (k=0; k<3; ++k) {
                for(l=0; l<3; ++l) {
                    if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                        cp = hcsp3c[i][k];
                        ++m;
                    }
                }
            }
            if (m != 1) continue; // 1 common central particle between SP3_i and SP3_j
            if (cp != hcsp3c[i][j2]) continue;

            for (k=j+1; k < nmem_sp3c[hcsp3c[i][j2]]; ++k) { // loop over all sp3c_k
                if (mem_sp3c[hcsp3c[i][j2]][k] <= i) continue;
                flg1 = Bonds_BondCheck(h1i, hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) && Bonds_BondCheck(h1j,
                                                                                                      hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]);
                flg2 = Bonds_BondCheck(h2i, hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]) && Bonds_BondCheck(h2j,
                                                                                                      hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                flg1 = flg1 && flg2;
                flg2 = Bonds_BondCheck(h1i, hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]) && Bonds_BondCheck(h1j,
                                                                                                      hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                flg3 = Bonds_BondCheck(h2i, hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) && Bonds_BondCheck(h2j,
                                                                                                      hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]);
                flg2 = flg2 && flg3;
                flg1 = flg1 || flg2;
                if (flg1==0) continue; // spindle particles not bonded correctly
                flg1 = hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4];
                flg2 = hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3];
                flg1 = flg1 || flg2;
                flg2 = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4];
                flg3 = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3];
                flg2 = flg2 || flg3;
                flg1 = flg1 || flg2;
                if (flg1==1) continue; // common spindle particles

                n = 0;
                for (l=0; l<3; ++l) {
                    for (m=0; m<3; ++m) {
                        if (hcsp3c[i][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][m]) {
                            if (cp != hcsp3c[i][l]) break;
                            ++n;
                        }
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][m]) {
                            if(cp != hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) break;
                            ++n;
                        }
                    }
                    if (m<3) break;
                }
                if (l<3) continue;
                if (n != 2) continue; // only one common particle, cp

                for (l=0; l<3; ++l) { // sp3 particles can't be spindle particles
                    x = hcsp3c[i][l];
                    flg1 = x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4];
                    x = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                    flg2 = x == hcsp3c[i][3] || x == hcsp3c[i][4] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4];
                    x = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                    flg3 = x == hcsp3c[i][3] || x == hcsp3c[i][4] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || x == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                    if (flg1==1 || flg2==1 || flg3==1) break;
                }
                if (l<3) continue;

                // Check for bonds between sp3 particles and other spindles
                for (l=0; l<3; ++l) {
                    if (hcsp3c[i][l] != cp) {
                        flg1 = Bonds_BondCheck(hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) || Bonds_BondCheck(
                                hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]);
                        flg1 = flg1 || Bonds_BondCheck(hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) || Bonds_BondCheck(
                                hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                        if (flg1==1) break;
                    }
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] != cp) {
                        flg1 = Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[i][3]) || Bonds_BondCheck(
                                hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[i][4]);
                        flg1 = flg1 || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) || Bonds_BondCheck(
                                hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                        if (flg1==1) break;
                    }
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] != cp) {
                        flg1 = Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[i][3]) || Bonds_BondCheck(
                                hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[i][4]);
                        flg1 = flg1 || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) || Bonds_BondCheck(
                                hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]);
                        if (flg1==1) break;
                    }
                }
                if (l<3) continue;

                m = 0;  // find uncommon particles from 5A_i
                for (l=0; l<3; ++l) {
                    if (hcsp3c[i][l] != cp) {
                        ia[m] = hcsp3c[i][l];
                        m++;
                    }
                }
                flg1 = flg2 = 0;
                ja[0]=ja[1]=-1; // find uncommon particles from 5A_j
                for (l=0; l<3; ++l) {
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == cp) continue;
                    if (Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], ia[0]) == 1) {
                        flg1 = 1;
                        ja[0] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                        if(Bonds_BondCheck(ja[0], ia[1]) == 1) break;
                    }
                    else ja[1] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                }
                if (l<3) continue;
                if (!flg1) {
                    m = ia[0];
                    ia[0] = ia[1];
                    ia[1] = m;
                    for (l=0; l<3; ++l) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == cp) continue;
                        if (Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], ia[0]) == 1) {
                            flg2 = 1;
                            ja[0] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                            if (Bonds_BondCheck(ja[0], ia[1]) == 1) break;
                        }
                        else ja[1] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                    }
                    if (l<3 || flg2==0) continue;
                }
                if (!(flg1==1 || flg2==1)) continue; // found four uncommon particles from 5A_i and 5A_j and found the two that are bonded between these
                if (ja[0]==-1 || ja[1]==-1) continue;
                if (Bonds_BondCheck(ja[1], ia[0]) == 1 || Bonds_BondCheck(ja[1], ia[1]) == 1) continue;
                flg1 = 0;
                ka[0]=ka[1]=-1;
                for (l=0; l<3; ++l) {
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] == cp) continue;
                    if (Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], ia[1])) {
                        flg1 = 1;
                        ka[0] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                        if (Bonds_BondCheck(ka[0], ia[0]) == 1 || Bonds_BondCheck(ka[0], ja[0]) == 1 || Bonds_BondCheck(ka[0], ja[1]) == 1) break;
                    }
                    else ka[1] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                }
                if (l<3 || ka[0]==-1 || ka[1]==-1) continue;
                if (Bonds_BondCheck(ka[1], ja[1]) == 0) continue;
                if (Bonds_BondCheck(ka[1], ia[1]) == 1 || Bonds_BondCheck(ka[1], ia[0]) == 1 || Bonds_BondCheck(ka[1], ja[0]) == 1) continue;

                // We've now found an HCP cluster
                if (nHCP == mHCP) {
                    hcHCP= resize_2D_int(hcHCP, mHCP, mHCP + incrStatic, clusSize, -1);
                    mHCP= mHCP + incrStatic;
                }

                hcHCP[nHCP][0] = cp;
                l = 1;
                for (m=0; m<3; ++m) {
                    if(hcsp3c[i][m] != cp) hcHCP[nHCP][l++] = hcsp3c[i][m];
                    if(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][m] != cp) hcHCP[nHCP][l++] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][m];
                    if(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][m] != cp) hcHCP[nHCP][l++] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][m];
                }
                hcHCP[nHCP][l++] = hcsp3c[i][3];
                hcHCP[nHCP][l++] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                hcHCP[nHCP][l++] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3];
                hcHCP[nHCP][l++] = hcsp3c[i][4];
                hcHCP[nHCP][l++] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                hcHCP[nHCP][l] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4];
                quickSort(&hcHCP[nHCP][1], 6);
                quickSort(&hcHCP[nHCP][7], 6);

                Cluster_Write_HCP(i, j, j2, k);
            }
        }
        }
    }
}

void Cluster_Write_HCP(int i, int j, int j2, int k) {

    for (int counter = 0; counter < 3; counter++){
        if (sHCP[hcsp3c[i][counter]] == 'C') sHCP[hcsp3c[i][counter]] = 'B';
        if (sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][counter]] == 'C') sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][counter]] = 'B';
        if (sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][counter]] == 'C') sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][counter]] = 'B';
    }
    for (int counter = 3; counter < 5; counter++) {
        sHCP[hcsp3c[i][counter]] = 'O';
        sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][counter]] = 'F';
        sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][counter]] = 'H';
    }

    ++nHCP;
}