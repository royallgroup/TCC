#include "13K.h"
#include "globals.h"
#include "tools.h"

int Clusters_Get13K(int sp3c_i, int sp3c_j, int the6A_i) {
    /* Function Clusters_Get13K - Take an 11F particle and determine if it meets the criteria for the presence of a 13K
     *
     * f: Frame number currently being analysed
     * sp3c_i: The id of a relevant 5A cluster
     * sp3c_i: The id of a different relevant 5A cluster
     * the6A_i: The id of a relevant 6A ring
     *
     * Returns 1 if a 13K is successfully detected
     * Returns 0 if no 13K is detected
     * 13K arrays are edited in place to add new 13K
     */
    int i, j, k, l;
    int sp3c_i_unc, sp3c_j_unc, ep[2], eclus5A[2], tmp;
    int clusSize=13;

    sp3c_i_unc=sp3c_j_unc=ep[0]=ep[1]=eclus5A[0]=eclus5A[1]=-1;

    // Clusters 5A_i and 5A_j are new to the 13K, they have:
    // two spindle particles in the 11F
    // one sp3 particle is the central particle (rc) from 11F
    // one sp3 particle is from a 5A in the 11F (sp3c_i/j_unc)
    // one sp3 particle is distinct from the 11F

    // Identification of sp3c_i_unc
    k=0;
    for (i=0; i<3; i++) {
        // Make sure the sp3 ring of 5A_i does not contain the central particle in the 11F
        if (hcsp3c[sp3c_i][i] != hc11F[n11F][0]) { ;
            // Make sure none of the sp3 ring of 5A_i is in sp4 ring of the specified 6A
            for (j = 0; j < 4; j++) {
                if (hcsp3c[sp3c_i][i] == hcsp4c[the6A_i][j]) {
                    break;
                }
            }
            // if none of the sp3 ring is in sp4
            if (j == 4) {
                if (k >= 1) {
                    k++;
                    break;
                }
                sp3c_i_unc = hcsp3c[sp3c_i][i];
                k++;
            }
        }
    }
    if (k!=1) return 0;

    // Identification of sp3c_j_unc
    k=0;
    for (i=0; i<3; i++) {
        // Make sure the sp3 ring of 5A_i does not contain the central particle in the 11F
        if (hcsp3c[sp3c_j][i] != hc11F[n11F][0]) {
            // Make sure none of the sp3 ring of 5A_j is in sp4 ring of the specified 6A
            for (j = 0; j < 4; j++) {
                if (hcsp3c[sp3c_j][i] == hcsp4c[the6A_i][j]) {
                    break;
                }
            }
            if (j == 4) {
                if (k >= 1) {
                    k++;
                    break;
                }
                sp3c_j_unc = hcsp3c[sp3c_j][i];
                k++;
            }
        }
    }
    if (k!=1) return 0;

    // Try to identify the new particle in 5A_j which is not in 11F
    k=0;
    // loop over all 5A clusters of which hc13K[n13K][0] is a member
    for (i=0; i < nmem_sp3c[hc11F[n11F][0]]; ++i) {
        if (mem_sp3c[hc11F[n11F][0]][i] != sp3c_i && mem_sp3c[hc11F[n11F][0]][i] != sp3c_j) {
            if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][3] == hcsp3c[sp3c_i][3]) {
                if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][4] == hcsp3c[sp3c_i][4]) {
                    for (j = 0; j < 3; j++) {
                        if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][j] != hc11F[n11F][0]) {
                            if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][j] != sp3c_i_unc) {
                                if (k == 1) {
                                    return 0;
                                }
                                else {
                                    tmp = hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][j];
                                    // check that tmp is not already in 11F
                                    for (l = 0; l < 11; l++) {
                                        if (tmp == hc11F[n11F][l]) break;
                                    }
                                    if (l == 11) {
                                        ep[0] = tmp;
                                        eclus5A[0] = mem_sp3c[hc11F[n11F][0]][i];
                                        k++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (k!=1) return 0;
    // have found particle uncommon to 11F which forms a 5A cluster sharing spindles of sp3c_i and ring particles hc13K[n13K][0] and sp3c_i_unc

    // Try to identify the new particle in 5A_j which is not in 11F
    k=0;
    for (i=0; i < nmem_sp3c[hc11F[n11F][0]]; ++i) { // loop over all sp3c which hc13K[n13K][0] is a member
        if (mem_sp3c[hc11F[n11F][0]][i] != sp3c_i && mem_sp3c[hc11F[n11F][0]][i] != sp3c_j) {
            if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][3] == hcsp3c[sp3c_j][3]) {
                if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][4] == hcsp3c[sp3c_j][4]) {
                    for (j = 0; j < 3; j++) {
                        if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][j] != hc11F[n11F][0]) {
                            if (hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][j] != sp3c_j_unc) {
                                if (k == 1) {
                                    return 0;
                                }
                                else {
                                    tmp = hcsp3c[mem_sp3c[hc11F[n11F][0]][i]][j];

                                    for (l = 0; l < 11; l++) {  // check temp not already in 11F
                                        if (tmp == hc11F[n11F][l]) break;
                                    }
                                    if (l == 11) {
                                        ep[1] = tmp;
                                        eclus5A[1] = mem_sp3c[hc11F[n11F][0]][i];
                                        k++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (k!=1) return 0;
    // have found particle uncommon to 11F which forms a 5A cluster sharing spindles of sp3c_j and ring particles hc13K[n13K][0] and sp3c_j_unc

    if(n13K == m13K) {
        hc13K= resize_2D_int(hc13K, m13K, m13K + incrStatic, clusSize, -1);
        m13K= m13K + incrStatic;
    }
    // hc13K key: (11F, extra SP3 ring particle to make 5A #1, extra SP3 ring particle to make 5A #2)

    for (i=0; i<11; i++) hc13K[n13K][i] = hc11F[n11F][i];
    hc13K[n13K][11]=ep[0];
    hc13K[n13K][12]=ep[1];

    quickSort(&hc13K[n13K][11], 2);
    quickSort(&eclus5A[0],2);

    Cluster_Write_13K();

    return 1;
}

void Cluster_Write_13K() {
    int i;
    s13K[hc13K[n13K][0]] = 'S';
    for(i=1; i<11; i++) {
        if (s13K[hc13K[n13K][i]] == 'C') s13K[hc13K[n13K][i]] = 'B';
    }
    if(s13K[hc13K[n13K][11]] != 'S') s13K[hc13K[n13K][11]] = 'O';
    if(s13K[hc13K[n13K][12]] != 'S') s13K[hc13K[n13K][12]] = 'O';
}