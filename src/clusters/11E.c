#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include <clusters/12D.h>
#include "11E.h"

void Clusters_Get11E_12D(int i, int j, int sp1, int sp2i, int sp2j) {

    //!  An 11E cluster is the intersection of a 9B and a 7A
    /*!
   *  Find 11E clusters
   *  An 11E is a 9B and a 7A where:
   *      - One 7A spindle particle is common to one of the uncommon spindle particles, sd1,
   *      in the 7A clusters constituting the 9B cluster.
   *      - The other spindle particle of the additional 7A is labeled sd3 and is bonded to the
   *      other uncommon spindle particle sd2 in 9B and the common spindle particle sc of 9B.
   *      Of the 7A cluster sp5 ring particles, one is common to the common with sc, one is common
   *      with sd2, and one is common to one of the uncommon sp5 ring particles  of the 9B cluster.
   *      The final two sp5 ring particles are distinct from the 9B cluster.
   *
   *  Cluster output: OOOOBBBBBBB
   *  Storage order: unknown
   *
   */

    int k, l, m, n;
    int trial[11];
    int break_out,break_out2;
    int flg1, flg2, flg3;
    int clusSize=11;

    for(k=i+1; k < nsp5c; k++) { // loop over all 7A_k
        if(k == j) continue;
        if(hcsp5c[k][5] == sp2j && hcsp5c[k][6] != sp2i) { // one 7A_k spindle is sp2j, one is not sp2i
            if(Bonds_BondCheck(hcsp5c[k][6], sp1) && Bonds_BondCheck(hcsp5c[k][6], sp2i)) { // non sp2j 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2j;
                trial[2] = sp2i;
                trial[3] = hcsp5c[k][6];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (hcsp5c[k][l] == sp1) {  // one SP5 ring particle of 7A_i common to common spindle of 9B
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l] == sp2i) { // one SP5 ring particle of 7A_i common to the other uncommon spindle of 9B
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {   // one SP5 ring particle of 7A_i common to the uncommon particles of the SP5 rings in the 7A constituting 9B
                        if (hcsp5c[k][l] == hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]= hcsp5c[k][l];
                            break_out2=1;
                            break;
                        }
                    }
                    if (break_out2==1) continue;
                    if (n==6) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l] == trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l] == trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue; // fetched final 4 particles from 9B not in 7A_i

                if(n11E == m11E) {
                    hc11E= resize_2D_int(hc11E, m11E, m11E + incrStatic, clusSize, -1);
                    m11E= m11E + incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D == 1) n12D += Clusters_Get12D(j, k, sp2i, hcsp5c[k][6]);

                ++n11E;
            }
        }
        if(hcsp5c[k][6] == sp2j && hcsp5c[k][5] != sp2i) { // one 7A_k spindle is sp2j, one is not sp2i
            if(Bonds_BondCheck(hcsp5c[k][5], sp1) && Bonds_BondCheck(hcsp5c[k][5], sp2i)) { // non sp2j 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2j;
                trial[2] = sp2i;
                trial[3] = hcsp5c[k][5];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (hcsp5c[k][l] == sp1) {
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l] == sp2i) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (hcsp5c[k][l] == hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]= hcsp5c[k][l];
                            break_out2=1;
                            break;
                        }
                    }
                    if (break_out2==1) continue;
                    if (n==6) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l] == trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l] == trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E == m11E) {
                    hc11E= resize_2D_int(hc11E, m11E, m11E + incrStatic, clusSize, -1);
                    m11E= m11E + incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D == 1) n12D += Clusters_Get12D(j, k, sp2i, hcsp5c[k][5]);
                ++n11E;
            }
        }
    }
    for(k=j+1; k < nsp5c; k++) {
        if(hcsp5c[k][5] == sp2i && hcsp5c[k][6] != sp2j) {  // one 7A_k spindle is sp2i, one is not sp2j
            if(Bonds_BondCheck(hcsp5c[k][6], sp1) && Bonds_BondCheck(hcsp5c[k][6], sp2j)) { // non sp2i 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2i;
                trial[2] = sp2j;
                trial[3] = hcsp5c[k][6];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (hcsp5c[k][l] == sp1) {
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l] == sp2j) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (hcsp5c[k][l] == hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]= hcsp5c[k][l];
                            break_out2=1;
                            break;
                        }
                    }
                    if (break_out2==1) continue;
                    if (n==6) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l] == trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l] == trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E == m11E) {
                    hc11E= resize_2D_int(hc11E, m11E, m11E + incrStatic, clusSize, -1);
                    m11E= m11E + incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D == 1) n12D += Clusters_Get12D(j, k, sp2j, hcsp5c[k][6]);
                ++n11E;
            }
        }
        if(hcsp5c[k][6] == sp2i && hcsp5c[k][5] != sp2j) {  // one 7A_k spindle is sp2i, one is not sp2j
            if(Bonds_BondCheck(hcsp5c[k][5], sp1) && Bonds_BondCheck(hcsp5c[k][5], sp2j)) { // non sp2i 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2i;
                trial[2] = sp2j;
                trial[3] = hcsp5c[k][5];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (hcsp5c[k][l] == sp1) {
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l] == sp2j) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (hcsp5c[k][l] == hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]= hcsp5c[k][l];
                            break_out2=1;
                            break;
                        }
                    }
                    if (break_out2==1) continue;
                    if (n==6) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l] == trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l] == trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]= hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E == m11E) {
                    hc11E= resize_2D_int(hc11E, m11E, m11E + incrStatic, clusSize, -1);
                    m11E= m11E + incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D == 1) n12D += Clusters_Get12D(j, k, sp2j, hcsp5c[k][5]);
                ++n11E;
            }
        }
    }
}

void Clust_Write_11E() {
    int i;

    for(i=0; i<4; i++) {
        s11E[hc11E[n11E][i]] = 'O';
    }
    for(i=4; i<11; i++) {
        if (s11E[hc11E[n11E][i]] == 'C') s11E[hc11E[n11E][i]] = 'B';
    }
}