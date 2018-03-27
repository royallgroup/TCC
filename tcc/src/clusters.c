#include "globals.h"
#include "clusters.h"
#include "bonds.h"
#include "tools.h"
#include "11B.h"

void Clusters_Get9B_10B_11B_11E_12D() {    // Detect 9B, 10B, 11B, 11E & 12D
    int sp1, sp2i, sp2j;
    int sp5com[2];
    int i, j, k, l, m;
    int flg, fb1, fb2;
    int clusSize=9;

    sp1=sp2i=sp2j=-1;

    for (i=0; i<nsp5c-1; ++i) {  // loop over all 7A_i
        // POSSIBLE IMPROVEMENT!! - 2 loops: over all 7A clusters which each spindle is in
        for (j=i+1; j<nsp5c; ++j) {  // loop over all 7A_j
            flg = 0;
            if (hcsp5c[i][5] == hcsp5c[j][5] && hcsp5c[i][6] != hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][6])) { // spindle particles arranged
                    flg = 1;
                    sp1 = hcsp5c[i][5];   // s_com common spindle
                    sp2i = hcsp5c[i][6];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][6];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (hcsp5c[i][6] == hcsp5c[j][6] && hcsp5c[i][5] != hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][5], hcsp5c[j][5])) {
                    flg = 1;
                    sp1 = hcsp5c[i][6];   // s_com common spindle
                    sp2i = hcsp5c[i][5];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][5];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (hcsp5c[i][5] == hcsp5c[j][6] && hcsp5c[i][6] != hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][5])) {
                    flg = 1;
                    sp1 = hcsp5c[i][5];   // s_com common spindle
                    sp2i = hcsp5c[i][6];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][5];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (hcsp5c[i][6] == hcsp5c[j][5] && hcsp5c[i][5] != hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][5], hcsp5c[j][6])) {
                    flg = 1;
                    sp1 = hcsp5c[i][6];   // s_com common spindle
                    sp2i = hcsp5c[i][5];  // 2nd spindle particle of cluster 7A_i
                    sp2j = hcsp5c[j][6];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (flg==0) continue;

            fb1 = fb2 = 1;  // ensure the two distinct spindle particles are part of the other SP5 ring
            for (k=0; k<5; ++k) {
                if (sp2i == hcsp5c[j][k]) fb1 = 0;
                if (sp2j == hcsp5c[i][k]) fb2 = 0;
            }
            if (fb1 || fb2) continue;

            m = 0;  // check for two common SP5 particles
            for (k=0; k<5; ++k) {
                for (l=0; l<5; ++l) {
                    if (hcsp5c[i][k] == hcsp5c[j][l]) {
                        if (m==2) {m++; break; }
                        sp5com[m]=hcsp5c[i][k];
                        ++m;
                    }
                }
            }
            if (m!=2) continue;

            // Now we have found the 9B C2v cluster
            if (n9B==m9B) {
                hc9B=resize_2D_int(hc9B,m9B,m9B+incrStatic,clusSize,-1);
                m9B=m9B+incrStatic;
            }
            if (sp5com[0]<sp5com[1]) {
                hc9B[n9B][4]=sp5com[0];
                hc9B[n9B][5]=sp5com[1];
            }
            else {
                hc9B[n9B][4]=sp5com[1];
                hc9B[n9B][5]=sp5com[0];
            }

            if (sp2i<sp2j) {
                hc9B[n9B][6]=sp2i;
                hc9B[n9B][7]=sp2j;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[i][k],hc9B[n9B][4]) && hcsp5c[i][k]!=hc9B[n9B][7] && hcsp5c[i][k]!=hc9B[n9B][4]) {
                        hc9B[n9B][0]=hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k],hc9B[n9B][5]) && hcsp5c[i][k]!=hc9B[n9B][7] && hcsp5c[i][k]!=hc9B[n9B][5]) {
                        hc9B[n9B][1]=hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[j][k],hc9B[n9B][4]) && hcsp5c[j][k]!=hc9B[n9B][6] && hcsp5c[j][k]!=hc9B[n9B][4]) {
                        hc9B[n9B][2]=hcsp5c[j][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[j][k],hc9B[n9B][5]) && hcsp5c[j][k]!=hc9B[n9B][6] && hcsp5c[j][k]!=hc9B[n9B][5]) {
                        hc9B[n9B][3]=hcsp5c[j][k];
                    }
                }
            }
            else {
                hc9B[n9B][6]=sp2j;
                hc9B[n9B][7]=sp2i;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[j][k],hc9B[n9B][4]) && hcsp5c[j][k]!=hc9B[n9B][7] && hcsp5c[j][k]!=hc9B[n9B][4]) {
                        hc9B[n9B][0]=hcsp5c[j][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[j][k],hc9B[n9B][5]) && hcsp5c[j][k]!=hc9B[n9B][7] && hcsp5c[j][k]!=hc9B[n9B][5]) {
                        hc9B[n9B][1]=hcsp5c[j][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k],hc9B[n9B][4]) && hcsp5c[i][k]!=hc9B[n9B][6] && hcsp5c[i][k]!=hc9B[n9B][4]) {
                        hc9B[n9B][2]=hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k],hc9B[n9B][5]) && hcsp5c[i][k]!=hc9B[n9B][6] && hcsp5c[i][k]!=hc9B[n9B][5]) {
                        hc9B[n9B][3]=hcsp5c[i][k];
                    }
                }
            }
            hc9B[n9B][8]=sp1;
            Cluster_Write_9B();

            if (do10B==1) Clusters_Get10B(i, j);
            if (do11B==1) {
                if (Clusters_Get11B()) {
                    s11B[hc9B[n9B][8]] = 'S';
                    ++n11B;
                }
            }
            if (do11E==1) Clusters_Get11E_12D(i, j, sp1, sp2i, sp2j);

            ++n9B;
        }
    }
}

void Cluster_Write_9B() {
    // hc9B key: (SP5_lowerd_to_4, SP5_lowerd_to_5, SP5_higherd_to_4, SP5_higherd_to_5, SP5_i_j_com_lower, SP5_i_j_com_higher, sp5c_d1_lower, sp5c_d2_higher, s_com)
    int i;
    for(i=0; i<6; i++) {
        if (s9B[hc9B[n9B][i]] == 'C') s9B[hc9B[n9B][i]] = 'B';
    }
    if (s9B[hc9B[n9B][6]] != 'S') s9B[hc9B[n9B][6]] = 'O';
    if (s9B[hc9B[n9B][7]] != 'S') s9B[hc9B[n9B][7]] = 'O';
    s9B[hc9B[n9B][8]] = 'S';
}

void Clusters_Get10B(int i, int j) {        // Return 1 if 9B is also 10B cluster
    int k,l,m;
    int flg1, flg2;
    int trial[10];
    int break_out;
    int clusSize=10;

    for (k=j+1; k<nsp5c; ++k) {  // loop over all 7A_k
        if (hcsp5c[k][5] == hc9B[n9B][8]) {    // check one spindle of 7A_k is the common spindle of 9B (hc9B[id9B][.] at this point)
            if (Bonds_BondCheck(hcsp5c[k][6], hc9B[n9B][6])==0) continue;  // check other spindle of 7A_k is bonded to spindle d1 of 9B
            if (Bonds_BondCheck(hcsp5c[k][6], hc9B[n9B][7])==0) continue;  // check other spindle of 7A_k is bonded to spindle d2 of 9B

            flg1=0;
            flg2=0;
            for (l=0;l<5;l++) {
                if (hcsp5c[k][l]==hc9B[n9B][6]) {
                    flg1=1;
                    continue;
                }
                if (hcsp5c[k][l]==hc9B[n9B][7]) {
                    flg2=1;
                    continue;
                }
            }
            if (flg1==0 || flg2==0) continue;
            trial[6]=hc9B[n9B][6];
            trial[7]=hc9B[n9B][7];
            trial[8]=hcsp5c[k][6];
            trial[9]=hc9B[n9B][8];

            m=0;
            break_out=0;
            for (l=0;l<6;l++) {
                if (hc9B[n9B][l]==hcsp5c[k][6]) continue;
                if (m==5) {
                    m++;
                    break_out=1;
                    break;
                }
                trial[m]=hc9B[n9B][l];
                m++;
            }
            if (break_out==1 || m!=5) continue;

            break_out=0;
            for (l=0;l<5;l++) {
                if (hcsp5c[k][l]==hc9B[n9B][6]) continue;
                if (hcsp5c[k][l]==hc9B[n9B][7]) continue;
                for (m=0;m<5;m++) {
                    if (hcsp5c[k][l]==trial[m]) break;
                }
                if (m==5) {
                    trial[5]=hcsp5c[k][l];
                    break_out++;
                }
            }
            if (break_out!=1) continue;

            if (n10B==m10B) { hc10B=resize_2D_int(hc10B,m10B,m10B+incrStatic,clusSize,-1);
                m10B=m10B+incrStatic;
            }
            // Now we have found the 10B C3v cluster
            // ###### NOTE #####
            // we have sterically assumed that
            // 1) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_i
            // (this member of 7A_i was uncommon to the SP5 ring of 7A_j)
            // 2) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_j
            // (this member of 7A_j was uncommon to the SP5 ring of 7A_i)

            quickSort(&trial[0],6);
            quickSort(&trial[6],3);
            for (l=0;l<10;l++) hc10B[n10B][l]=trial[l];

            Cluster_Write_10B();
        }

        if (hcsp5c[k][6] == hc9B[n9B][8]) {    // check one spindle of 7A_k is the common spindle of 9B (hc9B[id9B][.] at this point)
            if (Bonds_BondCheck(hcsp5c[k][5], hc9B[n9B][6])==0) continue;  // check other spindle of 7A_k is bonded to spindle d1 of 9B
            if (Bonds_BondCheck(hcsp5c[k][5], hc9B[n9B][7])==0) continue;  // check other spindle of 7A_k is bonded to spindle d2 of 9B

            flg1=0;
            flg2=0;
            for (l=0;l<5;l++) {
                if (hcsp5c[k][l]==hc9B[n9B][6]) {
                    flg1=1;
                    continue;
                }
                if (hcsp5c[k][l]==hc9B[n9B][7]) {
                    flg2=1;
                    continue;
                }
            }
            if (flg1==0 || flg2==0) continue;
            trial[6]=hc9B[n9B][6];
            trial[7]=hc9B[n9B][7];
            trial[8]=hcsp5c[k][5];
            trial[9]=hc9B[n9B][8];

            m=0;
            break_out=0;
            for (l=0;l<6;l++) {
                if (hc9B[n9B][l]==hcsp5c[k][5]) continue;
                if (m==5) {
                    m++;
                    break_out=1;
                    break;
                }
                trial[m]=hc9B[n9B][l];
                m++;
            }
            if (break_out==1 || m!=5) continue;

            break_out=0;
            for (l=0;l<5;l++) {
                if (hcsp5c[k][l]==hc9B[n9B][6]) continue;
                if (hcsp5c[k][l]==hc9B[n9B][7]) continue;
                for (m=0;m<5;m++) {
                    if (hcsp5c[k][l]==trial[m]) break;
                }
                if (m==5) {
                    trial[5]=hcsp5c[k][l];
                    break_out++;
                }
            }
            if (break_out!=1) continue;

            if (n10B==m10B) { hc10B=resize_2D_int(hc10B,m10B,m10B+incrStatic,clusSize,-1);
                m10B=m10B+incrStatic;
            }
            // Now we have found the 10B C3v cluster
            // ###### NOTE #####
            // we have sterically assumed that
            // 1) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_i
            // (this member of 7A_i was uncommon to the SP5 ring of 7A_j)
            // 2) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_j
            // (this member of 7A_j was uncommon to the SP5 ring of 7A_i)

            quickSort(&trial[0],6);
            quickSort(&trial[6],3);
            for (l=0;l<10;l++) hc10B[n10B][l]=trial[l];

            // hc10B key: (ordered shell particles, s1, s2, s3 (ordered), s_com)
            Cluster_Write_10B();
        }
    }
}

void Cluster_Write_10B() {
    // hc10B key: (ordered shell particles, s1, s2, s3 (ordered), s_com)
    int i;

    for(i=0; i<6; i++) {
        if (s10B[hc10B[n10B][i]] == 'C') s10B[hc10B[n10B][i]] = 'B';
    }
    for(i=6; i<9; i++) {
        if (s10B[hc10B[n10B][i]] != 'S') s10B[hc10B[n10B][i]] = 'O';
    }
    s10B[hc10B[n10B][9]] = 'S';
    ++n10B;
}

void Clusters_Get11E_12D(int i, int j, int sp1, int sp2i, int sp2j) {    // Returns number of 11Es for a single 9B
    //  ###### NOTE #####
    //  for 11E C2 we sterically assume that given that two members of the SP5 ring of 7A_k are new, the other three are
    // 1) sp1
    // 2) sp2i/j
    // 3) is common with one of the SP5_j/i_unc

    int k, l, m, n;
    int trial[11];
    int break_out,break_out2;
    int flg1, flg2, flg3;
    int clusSize=11;

    for(k=i+1; k<nsp5c; k++) { // loop over all 7A_k
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
                    if (hcsp5c[k][l]==sp1) {  // one SP5 ring particle of 7A_i common to common spindle of 9B
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l]==sp2i) { // one SP5 ring particle of 7A_i common to the other uncommon spindle of 9B
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {   // one SP5 ring particle of 7A_i common to the uncommon particles of the SP5 rings in the 7A constituting 9B
                        if (hcsp5c[k][l]==hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]=hcsp5c[k][l];
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
                    trial[n]=hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue; // fetched final 4 particles from 9B not in 7A_i

                if(n11E == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D += Clusters_Get12D(j, k, sp2i, hcsp5c[k][6]);

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
                    if (hcsp5c[k][l]==sp1) {
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l]==sp2i) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (hcsp5c[k][l]==hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]=hcsp5c[k][l];
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
                    trial[n]=hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D += Clusters_Get12D(j, k, sp2i, hcsp5c[k][5]);
                ++n11E;
            }
        }
    }
    for(k=j+1; k<nsp5c; k++) {
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
                    if (hcsp5c[k][l]==sp1) {
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l]==sp2j) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (hcsp5c[k][l]==hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]=hcsp5c[k][l];
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
                    trial[n]=hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D += Clusters_Get12D(j, k, sp2j, hcsp5c[k][6]);
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
                    if (hcsp5c[k][l]==sp1) {
                        flg1=1;
                        continue;
                    }
                    if (hcsp5c[k][l]==sp2j) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (hcsp5c[k][l]==hc9B[n9B][m]) {
                            flg3=1;
                            trial[6]=hcsp5c[k][l];
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
                    trial[n]=hcsp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E][i]=trial[i];

                Clust_Write_11E();

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D += Clusters_Get12D(j, k, sp2j, hcsp5c[k][5]);
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

int Clusters_Get12D(int j, int k, int sp1, int sp2) {  // Return 1 if 12B is also 11E
    int l, m, n, o, p, q;
    int flg1, flg2;
    int break_out;
    int trial[12];
    int clusSize=12;

    if(k > j) m = k;
    else m = j;
    for (l=m+1; l<nsp5c; l++) {
        if ((hcsp5c[l][5] == sp1 && hcsp5c[l][6] == sp2) || (hcsp5c[l][6] == sp1 && hcsp5c[l][5] == sp2)) {
            flg1=flg2=0;
            p=11;
            q=0;
            for (n=0; n<5; n++) {
                if (hcsp5c[l][n]==hc11E[n11E][0]) {
                    flg1=1;
                    continue;
                }
                if (hcsp5c[l][n]==hc11E[n11E][1]) {
                    flg2=1;
                    continue;
                }
                break_out=0;
                for (o=4; o<11; o++) {
                    if (hcsp5c[l][n]==hc11E[n11E][o]) {
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
                trial[p]=hcsp5c[l][n];
                p++;
            }
            if (flg1==0 || flg2==0 || p!=12 || q!=2) continue;

            for (n=0; n<11; n++) trial[n]=hc11E[n11E][n];
            // hc12D key: (d1_unc, d2_unc, d3_unc, d4_unc, d12_com, d13_com, d24_com, d34_com, s_d1, s_d2, s_d3, s_com)

            if(n12D == m12D) {
                hc12D=resize_2D_int(hc12D,m12D,m12D+incrStatic,clusSize,-1);
                m12D=m12D+incrStatic;
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

void Clusters_Get9K() {
    // 9K are made from 2 6A clusters with a common spindle and two common SP4 ring particles
    int j2, j, k, l, m;
    int cp[2], scom, sother[2];
    int trial[9];
    int id_first_6A, id_second_6A;

    cp[0]=cp[1]=scom=sother[0]=sother[1]=-1;


    for(id_first_6A=0; id_first_6A<nsp4c-1; ++id_first_6A) {   // loop over all sp4c_i
        for (j2=4; j2<6; j2++) {    // loop over all spindles of sp4c_i
            for (j=0; j<nmem_sp4c[hcsp4c[id_first_6A][j2]]; ++j) {
                id_second_6A = mem_sp4c[hcsp4c[id_first_6A][j2]][j];
                if (id_second_6A<=id_first_6A) continue; // don't find again 9K twice

                m=0;        // sp4c_i and sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly one common spindle
                for(k=4; k<6; ++k) {
                    for(l=4; l<6; ++l) {
                        if(hcsp4c[id_first_6A][k] == hcsp4c[id_second_6A][l]) {
                            m++;
                            scom=hcsp4c[id_first_6A][k];
                        }
                    }
                }
                if(m!=1) continue;

                if (hcsp4c[id_first_6A][4]==scom) sother[0]=hcsp4c[id_first_6A][5];
                else sother[0]=hcsp4c[id_first_6A][4];
                if (hcsp4c[id_second_6A][4]==scom) sother[1]=hcsp4c[id_second_6A][5];
                else sother[1]=hcsp4c[id_second_6A][4];

                m=0;        // check sother[0] is not in cluster sp4c_mem_sp4c[sp4c[i][j2]][j]
                for(k=0; k<6; ++k) {
                    if(sother[0] == hcsp4c[id_second_6A][k] || sother[1] == hcsp4c[id_first_6A][k]) {
                        m++;
                    }
                }
                if(m!=0) continue;

                m=0;        // SP4 ring from sp4c_i and SP4 ring from sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly two common particles
                for(k=0; k<4; ++k) {
                    for(l=0; l<4; ++l) {
                        if(hcsp4c[id_first_6A][k] == hcsp4c[id_second_6A][l]) {
                            if (m>=2) {
                                m=3;
                                break;
                            }
                            cp[m]=hcsp4c[id_first_6A][k];
                            m++;
                        }
                    }
                    if (m>2) break;
                }
                if(m!=2) continue;

                // hc9K key: (common_SP4_1, common_SP4_2, other_SP4*4, other_spindle_1, other_spindle_2, scom)

                trial[0]=cp[0];
                trial[1]=cp[1];
                trial[6]=sother[0];
                trial[7]=sother[1];
                trial[8]=scom;
                quickSort(&trial[0],2);
                quickSort(&trial[6],2);

                m=2;        // find uncommon SP4 ring particles in sp4c[i][k]
                for(k=0; k<4; ++k) {
                    for (l=0; l<2; l++) {
                        if (hcsp4c[id_first_6A][k]==cp[l]) break;
                    }
                    if (l==2) {
                        if (m>=4) {
                            m++;
                            break;
                        }
                        trial[m]=hcsp4c[id_first_6A][k];
                        m++;
                    }
                }
                if(m!=4) continue;

                for(k=0; k<4; ++k) { // find uncommon SP4 ring particles in sp4c[mem_sp4c[sp4c[i][j2]][j]][k]
                    for (l=0; l<2; l++) {
                        if (hcsp4c[id_second_6A][k]==cp[l]) break;
                    }
                    if (l==2) {
                        if (m>=6) {
                            m++;
                            break;
                        }
                        trial[m]=hcsp4c[id_second_6A][k];
                        m++;
                    }
                }
                if(m!=6) continue;

                quickSort(&trial[2],4);

                Cluster_Write_9k(trial);
            }
        }
    }
}

void Cluster_Write_9k(const int trial[]) {
    int i;
    int clusSize = 9;

    if(n9K == m9K) {
        hc9K=resize_2D_int(hc9K,m9K,m9K+incrStatic,clusSize,-1);
        m9K=m9K+incrStatic;
    }
    for (i=0; i<9; i++) hc9K[n9K][i]=trial[i];

    for(i=0; i<6; i++) {
        if (s9K[hc9K[n9K][i]] == 'C') s9K[hc9K[n9K][i]] = 'B';
    }
    s9K[hc9K[n9K][6]] = 'O';
    s9K[hc9K[n9K][7]] = 'O';
    s9K[hc9K[n9K][8]] = 'O';

    ++n9K;
}

void Clusters_Get10K() { // Detect 10K clusters
    // A 10K is a 9K with a SINGLE particle bonded to the common spindle of 9K.
    int bonded_to_spindle_id, num_extra_particles, extra_particle = 0;
    int id_9K, id_9K_common;

    for (id_9K=0; id_9K<n9K; id_9K++) {
        id_9K_common = hc9K[id_9K][8];
        if (num_bonds[id_9K_common] < 10) {
            num_extra_particles = 0;
            for (bonded_to_spindle_id = 0; bonded_to_spindle_id < num_bonds[id_9K_common]; bonded_to_spindle_id++) {
                if (is_particle_in_9K(id_9K, bNums[id_9K_common][bonded_to_spindle_id])) {
                    num_extra_particles++;
                    extra_particle = bNums[id_9K_common][bonded_to_spindle_id];
                }
            }
            if (num_extra_particles == 1) {
                Cluster_Write_10K(id_9K, extra_particle);
            }
        }
    }
}

int is_particle_in_9K(int id_9K, int id_particle){
    // Returns 0 if particle is not in the specified 9K, else returns 1
    int i;

    for (i=0; i<9; i++) {
        if (id_particle==hc9K[id_9K][i]){
            return 0;
        }
    }
    return 1;
}

void Cluster_Write_10K(int id_9k, int extra_particle) {
    // hc10K key: (common_SP4_1, common_SP4_2, other_SP4*4, other_spindle_1, other_spindle_2, scom, ep)

    int i;
    int clusSize=10;

    if(n10K == m10K) {
        hc10K=resize_2D_int(hc10K,m10K,m10K+incrStatic,clusSize,-1);
        m10K=m10K+incrStatic;
    }

    for(i=0; i<9; i++) {
        hc10K[n10K][i]=hc9K[id_9k][i];
    }
    hc10K[n10K][9]=extra_particle;

    for(i=0; i<6; i++) {
        if (s10K[hc10K[n10K][i]] == 'C') s10K[hc10K[n10K][i]] = 'B';
    }
    s10K[hc10K[n10K][6]] = 'O';
    s10K[hc10K[n10K][7]] = 'O';
    s10K[hc10K[n10K][8]] = 'O';
    s10K[hc10K[n10K][9]] = 'O';

    n10K++;
}

void Clusters_Get10A() { // Detect 10A D4d clusters
    int i, j, j2, k, l, m;
    char errMsg[1000];
    int clusSize=10;
    int *used_sp4b;

    used_sp4b=malloc(nsp4b*sizeof(int)); if (used_sp4b==NULL) { sprintf(errMsg,"Clusters_Get10A(): used_sp4b[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i<nsp4b; ++i) used_sp4b[i] = 0;

    for (i=0; i<nsp4b-1; ++i) {  // loop over all sp4b_i
        for (j2=0; j2<nsp4b; ++j2) used_sp4b[j2] = 0;
        used_sp4b[i]=1;
        for (j2=0; j2<num_bonds[hcsp4b[i][0]]; ++j2) {
            for (j=0; j<nmem_sp4b[bNums[hcsp4b[i][0]][j2]]; ++j) {    // loop over sp4b_j
                if (mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]<=i) continue;
                if (used_sp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]]==1) continue;
                used_sp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]]=1;

                if(hcsp4b[i][4] == hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][4]) continue;
                if(Bonds_BondCheck(hcsp4b[i][4], hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][4])) continue;
                for(k=0; k<5; ++k){
                    for(l=0; l<5; ++l){
                        if(hcsp4b[i][k] == hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][l]) break;
                    }
                    if(l<5) break;
                }
                if(k<5) continue;
                for(k=0; k<4; ++k){
                    m = 0;
                    for(l=0; l<4; ++l) if(Bonds_BondCheck(hcsp4b[i][k], hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][l])) ++m;
                    if(m!=2) break;
                }
                // ERROR: Need to check converse, i. e. each SP4 ring particle from sp4b_j bonded to to exactly two particles from sp4b_i
                if(k==4) {
                    if(n10A == m10A) {
                        hc10A=resize_2D_int(hc10A,m10A,m10A+incrStatic,clusSize,-1);
                        m10A=m10A+incrStatic;
                    }

                    // hc10A key: (SP4s going up, spindles going up)

                    hc10A[n10A][0]=hcsp4b[i][0];
                    hc10A[n10A][1]=hcsp4b[i][1];
                    hc10A[n10A][2]=hcsp4b[i][2];
                    hc10A[n10A][3]=hcsp4b[i][3];
                    hc10A[n10A][4]=hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][0];
                    hc10A[n10A][5]=hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][1];
                    hc10A[n10A][6]=hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][2];
                    hc10A[n10A][7]=hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][3];
                    quickSort(&hc10A[n10A][0],8);
                    hc10A[n10A][8]=hcsp4b[i][4];
                    hc10A[n10A][9]=hcsp4b[mem_sp4b[bNums[hcsp4b[i][0]][j2]][j]][4];
                    quickSort(&hc10A[n10A][8],2);

                    Cluster_Write_10A();
                }
            }
        }
    }
    free(used_sp4b);
}

void Cluster_Write_10A() {
    int i;

    for(i=0; i<8; i++) {
        if (s10A[hc10A[n10A][i]] == 'C') s10A[hc10A[n10A][i]] = 'B';
    }
    s10A[hc10A[n10A][8]] = 'O';
    s10A[hc10A[n10A][9]] = 'O';

    ++n10A;
}

