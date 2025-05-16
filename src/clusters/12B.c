#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include "13A.h"
#include "12B.h"

//!  A 12B is the intersection of 6 7A clusters.
/*!
*  Find 12B clusters
*  A 12B is 6 7A clusters where:
*      - There is one central 7A with spindles A and B
*      - Every other 7A has a spindle common with the 7A spindle A and a spindle common with a ring particle of the first 7A
*
*  Cluster output: SBOOOOOBBBBB
*  Storage order: unknown
*
*/
void Clusters_Get12B_13A() {
    int j, k, l, m;
    int spindle_1, spindle_2;
    int sj1[5], sj2[5];
    int nSB1, nSB2;
    int flg;
    int break_out;
    int clusSize=12;

    for (int first_7A = 0; first_7A < nsp5c; ++first_7A) {
        int *first_7A_cluster = hcsp5c[first_7A];
        spindle_1 = first_7A_cluster[5];
        spindle_2 = first_7A_cluster[6];
        nSB2 = 0; // count up spindle bonds

        // Counting number of 7A spindle-spindle bonds
        nSB1 = count_7A_spindle_bonds(sj1, first_7A);

        if (nSB1 == 5 && do13A==1) {     // possibly found 13A, definately found 12B, now establish status
            for (j=first_7A+1; j < nsp5c; ++j) {
                if (spindle_1 == hcsp5c[j][5] || spindle_1 == hcsp5c[j][6]) {
                    for (k=0; k<5; ++k) {
                        for (l=0; l<5; ++l) {
                            if (first_7A_cluster[k] == hcsp5c[j][l]) break;
                        }
                        if (l<5) break;
                    }
                    if(k==5) { // got 13A, Check all sp5c[j][] - spindle_1 spindles are less than first_7A
                        for (k=0; k<first_7A; ++k) {
                            for(l=0; l<5; ++l) {
                                if (hcsp5c[j][l] == hcsp5c[k][5] && spindle_1 == hcsp5c[k][6]) break;
                                if (hcsp5c[j][l] == hcsp5c[k][6] && spindle_1 == hcsp5c[k][5]) break;
                            }
                            if(l<5) break; // index k < first_7A present
                        }
                        if(k==first_7A) break; // no index k < first_7A present
                    }
                }
            }
            if (j < nsp5c) { // 13A found
                if (n13A == m13A) {
                    hc13A= resize_2D_int(hc13A, m13A, m13A + incrStatic, clusSize + 1, -1);
                    m13A= m13A + incrStatic;
                }

                hc13A[n13A][0] = spindle_1;
                k = 1;
                if(first_7A_cluster[5] != spindle_1) hc13A[n13A][k++] = first_7A_cluster[5];
                if(first_7A_cluster[6] != spindle_1) hc13A[n13A][k++] = first_7A_cluster[6];
                if(hcsp5c[j][5] != spindle_1) hc13A[n13A][k++] = hcsp5c[j][5];
                if(hcsp5c[j][6] != spindle_1) hc13A[n13A][k] = hcsp5c[j][6];
                hc13A[n13A][3] = first_7A_cluster[0];
                hc13A[n13A][4] = first_7A_cluster[1];
                hc13A[n13A][5] = first_7A_cluster[2];
                hc13A[n13A][6] = first_7A_cluster[3];
                hc13A[n13A][7] = first_7A_cluster[4];
                hc13A[n13A][8] = hcsp5c[j][0];
                hc13A[n13A][9] = hcsp5c[j][1];
                hc13A[n13A][10] = hcsp5c[j][2];
                hc13A[n13A][11] = hcsp5c[j][3];
                hc13A[n13A][12] = hcsp5c[j][4];
                quickSort(&hc13A[n13A][1], 12);

                Clust_Write_13A();
            }
        }

        for (j=first_7A+1; j < nsp5c; ++j) {
            flg = spindle_2 == hcsp5c[j][5] && Bonds_BondCheck(spindle_1, hcsp5c[j][6]);
            flg = flg || (spindle_2 == hcsp5c[j][6] && Bonds_BondCheck(spindle_1, hcsp5c[j][5]));
            if (flg==1) {
                if (nSB2>=5) {
                    nSB2++;
                    break;
                }
                sj2[nSB2++] = j;
            }
        }

        if(nSB2 == 5 && do13A==1) { // possibly found 13A, definately found 12B, now establish status
            for (j=first_7A+1; j < nsp5c; ++j) {
                if (spindle_2 == hcsp5c[j][5] || spindle_2 == hcsp5c[j][6]) {
                    for (k=0; k<5; ++k) {
                        for (l=0; l<5; ++l) {
                            if (first_7A_cluster[k] == hcsp5c[j][l]) break;
                        }
                        if(l<5) break;
                    }
                    if (k==5) { // Check all sp5c[j][] - spindle_2 spindles are less than first_7A
                        for (k=0; k<first_7A; ++k) {
                            for (l=0; l<5; ++l) {
                                if (hcsp5c[j][l] == hcsp5c[k][5] && spindle_2 == hcsp5c[k][6]) break;
                                if (hcsp5c[j][l] == hcsp5c[k][6] && spindle_2 == hcsp5c[k][5]) break;
                            }
                            if (l<5) break;
                        }
                        if (k==first_7A) break;
                    }
                }
            }
            if (j < nsp5c){ // 13A found, is it unique
                if (n13A == m13A) {
                    hc13A= resize_2D_int(hc13A, m13A, m13A + incrStatic, clusSize + 1, -1);
                    m13A= m13A + incrStatic;
                }

                hc13A[n13A][0] = spindle_2;
                k = 1;
                if(first_7A_cluster[5] != spindle_2) hc13A[n13A][k++] = first_7A_cluster[5];
                if(first_7A_cluster[6] != spindle_2) hc13A[n13A][k++] = first_7A_cluster[6];
                if(hcsp5c[j][5] != spindle_2) hc13A[n13A][k++] = hcsp5c[j][5];
                if(hcsp5c[j][6] != spindle_2) hc13A[n13A][k] = hcsp5c[j][6];
                hc13A[n13A][3] = first_7A_cluster[0];
                hc13A[n13A][4] = first_7A_cluster[1];
                hc13A[n13A][5] = first_7A_cluster[2];
                hc13A[n13A][6] = first_7A_cluster[3];
                hc13A[n13A][7] = first_7A_cluster[4];
                hc13A[n13A][8] = hcsp5c[j][0];
                hc13A[n13A][9] = hcsp5c[j][1];
                hc13A[n13A][10] = hcsp5c[j][2];
                hc13A[n13A][11] = hcsp5c[j][3];
                hc13A[n13A][12] = hcsp5c[j][4];
                quickSort(&hc13A[n13A][1], 12);

                Clust_Write_13A();
            }
        }

        if ((nSB1 > 5) && (nSB2 > 5)) continue;

        for (j=0; j<first_7A; ++j) { // keep looking for 12B
            flg = spindle_1 == hcsp5c[j][5] && Bonds_BondCheck(spindle_2, hcsp5c[j][6]);
            flg = flg || (spindle_1 == hcsp5c[j][6] && Bonds_BondCheck(spindle_2, hcsp5c[j][5]));
            if (flg==1) {
                if (nSB1 >= 5) {
                    nSB1++;
                    break;
                }
                sj1[nSB1++] = j;
            }
        }

        if (nSB1 == 5) {
            if (n12B == m12B) {
                hc12B= resize_2D_int(hc12B, m12B, m12B + incrStatic, clusSize, -1);
                m12B= m12B + incrStatic;
            }
            hc12B[n12B][0] = spindle_1;
            hc12B[n12B][1] = spindle_2;
            hc12B[n12B][2] = first_7A_cluster[0];
            hc12B[n12B][3] = first_7A_cluster[1];
            hc12B[n12B][4] = first_7A_cluster[2];
            hc12B[n12B][5] = first_7A_cluster[3];
            hc12B[n12B][6] = first_7A_cluster[4];

            m = 7;
            break_out=0;
            for (j=0; j<5; ++j) {
                for (k=0; k<5; ++k) {
                    for (l=0; l<m; ++l) {
                        if(hc12B[n12B][l] == hcsp5c[sj1[j]][k]) break;
                    }
                    if (l==m) {
                        if (m==12) {
                            break_out=1;
                            break;
                        }
                        hc12B[n12B][m] = hcsp5c[sj1[j]][k];
                        m++;
                    }
                }
                if (break_out==1) break;
            }
            if (break_out==1 || m<12) continue;

            quickSort(&hc12B[n12B][2], 5);
            quickSort(&hc12B[n12B][7], 5);

            Clust_Write_12B();
        }

        for (j=0; j<first_7A; ++j) {
            flg = spindle_2 == hcsp5c[j][5] && Bonds_BondCheck(spindle_1, hcsp5c[j][6]);
            flg = flg || (spindle_2 == hcsp5c[j][6] && Bonds_BondCheck(spindle_1, hcsp5c[j][5]));
            if (flg==1) {
                if (nSB2 >= 5) {
                    nSB2++;
                    break;
                }
                sj2[nSB2++] = j;
            }
        }

        if(nSB2 == 5) {
            if (n12B == m12B) {
                hc12B= resize_2D_int(hc12B, m12B, m12B + incrStatic, clusSize, -1);
                m12B= m12B + incrStatic;
            }

            hc12B[n12B][0] = spindle_2;
            hc12B[n12B][1] = spindle_1;
            hc12B[n12B][2] = first_7A_cluster[0];
            hc12B[n12B][3] = first_7A_cluster[1];
            hc12B[n12B][4] = first_7A_cluster[2];
            hc12B[n12B][5] = first_7A_cluster[3];
            hc12B[n12B][6] = first_7A_cluster[4];

            m = 7;
            break_out=0;
            for(j=0; j<5; ++j){
                for(k=0; k<5; ++k){
                    for(l=0; l<m; ++l) if(hc12B[n12B][l] == hcsp5c[sj2[j]][k]) break;
                    if(l==m) {
                        if (m==12) {
                            break_out=1;
                            break;
                        }
                        hc12B[n12B][m] = hcsp5c[sj2[j]][k];
                        m++;
                    }
                }
                if (break_out==1) break;
            }
            if (break_out==1 || m<12) continue;

            quickSort(&hc12B[n12B][2], 5);
            quickSort(&hc12B[n12B][7], 5);

            Clust_Write_12B();

        }
    }
}

int count_7A_spindle_bonds(int *sj1, const int first_7A) {
    int flg;
    int *first_7A_cluster = hcsp5c[first_7A];
    int num_spindle_bonds = 0;

    for (int second_7A = first_7A + 1; second_7A < nsp5c; ++second_7A) {
        int *second_7A_cluster = hcsp5c[second_7A];
        flg = first_7A_cluster[5] == second_7A_cluster[5] && Bonds_BondCheck(first_7A_cluster[6], second_7A_cluster[6]);
        flg = flg || (first_7A_cluster[5] == second_7A_cluster[6] && Bonds_BondCheck(first_7A_cluster[6], second_7A_cluster[5]));
        if (flg == 1) {
            if (num_spindle_bonds >= 5) {
                num_spindle_bonds++;
                break;
            }
            sj1[num_spindle_bonds++] = second_7A;
        }
    }
    return num_spindle_bonds;
}

void Clust_Write_12B() {

    s12B[hc12B[n12B][0]] = 'S';
    if(s12B[hc12B[n12B][1]] == 'C') s12B[hc12B[n12B][1]] = 'B';
    if(s12B[hc12B[n12B][2]] != 'S') s12B[hc12B[n12B][2]] = 'O';
    if(s12B[hc12B[n12B][3]] != 'S') s12B[hc12B[n12B][3]] = 'O';
    if(s12B[hc12B[n12B][4]] != 'S') s12B[hc12B[n12B][4]] = 'O';
    if(s12B[hc12B[n12B][5]] != 'S') s12B[hc12B[n12B][5]] = 'O';
    if(s12B[hc12B[n12B][6]] != 'S') s12B[hc12B[n12B][6]] = 'O';
    if(s12B[hc12B[n12B][7]] == 'C') s12B[hc12B[n12B][7]] = 'B';
    if(s12B[hc12B[n12B][8]] == 'C') s12B[hc12B[n12B][8]] = 'B';
    if(s12B[hc12B[n12B][9]] == 'C') s12B[hc12B[n12B][9]] = 'B';
    if(s12B[hc12B[n12B][10]] == 'C') s12B[hc12B[n12B][10]] = 'B';
    if(s12B[hc12B[n12B][11]] == 'C') s12B[hc12B[n12B][11]] = 'B';

    ++n12B;
}