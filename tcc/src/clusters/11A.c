#include "11A.h"
#include "globals.h"
#include "tools.h"
#include "bonds.h"

void Clusters_Get11A() {
    // Detect 11A clusters

    int first_6A_id, second_6A_id, k, l, m;
    int scom, sother[2];

    scom=sother[0]=sother[1]=-1;

    for(first_6A_id=0; first_6A_id<nsp4c-1; ++first_6A_id){
        for(second_6A_id=first_6A_id+1; second_6A_id<nsp4c; ++second_6A_id) {
            m=0;
            for (k=4; k<6; k++) {
                for (l=4; l<6; l++) {
                    if (hcsp4c[first_6A_id][k] == hcsp4c[second_6A_id][l]) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        scom=hcsp4c[first_6A_id][k];
                        m++;
                    }
                }
                if (m>=2) break;
            }
            if(m==1) {  // one common spindle
                // Identify the non-common spindles
                if (scom == hcsp4c[first_6A_id][4]) {
                    sother[0] = hcsp4c[first_6A_id][5];
                } else {
                    sother[0] = hcsp4c[first_6A_id][4];
                }
                if (scom == hcsp4c[second_6A_id][4]) {
                    sother[1] = hcsp4c[second_6A_id][5];
                } else {
                    sother[1] = hcsp4c[second_6A_id][4];
                }

                if (Check_unique_6A_rings(first_6A_id, second_6A_id) == 0) {
                    if (Check_6A_rings_bonded(first_6A_id, second_6A_id) == 1) {
                        Cluster_Write_11A(first_6A_id, second_6A_id, sother, scom);
                    }
                }
            }
        }
    }
}

void Cluster_Write_11A(int first_6A_id, int second_6A_id, const int sother[], int scom) {
    int clusSize=11;
    int i;

    if (n11A == m11A) {
        hc11A = resize_2D_int(hc11A, m11A, m11A + incrStatic, clusSize, -1);
        m11A = m11A + incrStatic;
    }

    hc11A[n11A][0] = hcsp4c[first_6A_id][0];
    hc11A[n11A][1] = hcsp4c[first_6A_id][1];
    hc11A[n11A][2] = hcsp4c[first_6A_id][2];
    hc11A[n11A][3] = hcsp4c[first_6A_id][3];
    hc11A[n11A][4] = hcsp4c[second_6A_id][0];
    hc11A[n11A][5] = hcsp4c[second_6A_id][1];
    hc11A[n11A][6] = hcsp4c[second_6A_id][2];
    hc11A[n11A][7] = hcsp4c[second_6A_id][3];
    hc11A[n11A][8] = sother[0];
    hc11A[n11A][9] = sother[1];
    hc11A[n11A][10] = scom;

    for(i=0; i<8; i++) {
        if (s11A[hc11A[n11A][i]] == 'C') s11A[hc11A[n11A][i]] = 'B';
    }
    for(i=8; i<10; i++) {
        if (s11A[hc11A[n11A][i]] != 'S') s11A[hc11A[n11A][i]] = 'O';
    }
    s11A[hc11A[n11A][10]] = 'S';
    ++n11A;
}

int Check_6A_rings_bonded(int first_6A_id, int second_6A_id) {
    int i, j, num_bonds;
    // Check if there are two bonds between each particle in ring 1 and particles in ring 2
    // Returns 1 if all ring 1 particles have 2 bonds to ring 2 particles, return 0 if not
    // i loops over the particles in the first ring, j loops over the particles in the second ring
    for(i=0; i<4; ++i) { // loop through all first ring particles
        num_bonds = 0;
        for(j=0; j<4; ++j) { // loop through second ring particles
            if(Bonds_BondCheck(hcsp4c[first_6A_id][i], hcsp4c[second_6A_id][j])) {
                num_bonds++;
            }
        }
        if(num_bonds!=2) {
            return 0;
        }
    }
    return 1;
}

int Check_unique_6A_rings(int first_6A_id, int second_6A_id) {
    // Check that there are no common ring particles between two 6As.
    // Return 1 if the two 6A rings share a particle, else return 0
    int first_ring, second_ring;

    for(first_ring=0; first_ring<4; ++first_ring) {
        for(second_ring=0; second_ring<4; ++second_ring) {
            if(hcsp4c[first_6A_id][first_ring] == hcsp4c[second_6A_id][second_ring]) {
                return 1;
            }
        }
    }
    return 0;
}