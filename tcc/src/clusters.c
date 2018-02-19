#include "globals.h"
#include "clusters.h"
#include "bonds.h"
#include "tools.h"
#include "string.h"

void Clusters_Get6Z_C2v() {    // Detect 6Z clusters from 2 5A clusters
    int flg;
    int i, j, j2, k, l;
    int cnt;
    int s1a, s2a, s1b, s2b;
    int clusSize=6;

    s1a=s2a=s1b=s2b=-1;


    for (i=0; i<nsp3c-1; ++i) {  // loop over all 5A_i
        for (j2=0; j2<1; ++j2) {
            for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]; ++j) {  // loop over all 5A_j
                if (mem_sp3c[hcsp3c[i][j2]][j]<=i) continue;
                if (hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) continue;   // no spindles of 5A_i and 5A_j are common
                if (hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) continue;

                // for each cluster one spindle particle is in sp3 ring of other 5A and other spindle isn't
                cnt = 0;    // check one 5A_i spindle are in sp3 ring of 5A_j
                flg = hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][0] || hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][1] || hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][2];
                if (flg==1){
                    s1a = hcsp3c[i][3];
                    s2a = hcsp3c[i][4];
                    ++cnt;
                }
                flg = hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][0] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][1] || hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][2];
                if (flg==1){
                    s1a = hcsp3c[i][4];
                    s2a = hcsp3c[i][3];
                    ++cnt;
                }
                if (cnt != 1) continue;
                cnt = 0;    // check one 5A_j spindle are in sp3 ring of 5A_i
                flg = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[i][0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[i][1] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == hcsp3c[i][2];
                if (flg==1){
                    s1b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                    s2b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                    ++cnt;
                }
                flg = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[i][0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[i][1] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] == hcsp3c[i][2];
                if (flg==1){
                    s1b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                    s2b = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                    ++cnt;
                }
                if (cnt != 1) continue;
                if (!Bonds_BondCheck(s1a, s1b)) continue;   // check spindles of i and j in sp3 ring of j and i respectively are bonded

                cnt = 0;    // check 2 particles in the sp3 rings of 5A_i and 5A_j are common
                for (k=0; k<3; ++k){
                    for (l=0; l<3; ++l){
                        if(hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]){
                            ++cnt;
                            break;
                        }
                    }
                }
                if (cnt != 2) continue;

                if (n6Z==m6Z) {
                    hc6Z=resize_2D_int(hc6Z,m6Z,m6Z+incrStatic,clusSize,-1);
                    m6Z=m6Z+incrStatic;
                }
                // Now we have found the 6Z cluster
                if (s1a<s1b) {
                    hc6Z[n6Z][0]=s1a;    // insert cluster
                    hc6Z[n6Z][1]=s1b;
                    hc6Z[n6Z][2]=s2a;
                    hc6Z[n6Z][3]=s2b;
                }
                else {
                    hc6Z[n6Z][0]=s1b;    // insert cluster
                    hc6Z[n6Z][1]=s1a;
                    hc6Z[n6Z][2]=s2b;
                    hc6Z[n6Z][3]=s2a;
                }
                cnt=4;
                for (k=0; k<3; ++k) {
                    flg=1;
                    for (l=0; l<4; ++l){
                        if (hcsp3c[i][k]==hc6Z[n6Z][l]) {
                            flg=0;
                            break;
                        }
                    }
                    if (flg==1) { hc6Z[n6Z][cnt]=hcsp3c[i][k]; cnt++; }
                }
                if (hc6Z[n6Z][5]<hc6Z[n6Z][4]) {
                    k=hc6Z[n6Z][5];
                    hc6Z[n6Z][5]=hc6Z[n6Z][4];
                    hc6Z[n6Z][4]=k;
                }
                Cluster_Write_6Z();
            }
        }
    }
}

void Cluster_Write_6Z() {
    // hc6Z key: (5A_i_s_in_SP3_j, 5A_j_s_in_SP3_i, 5A_i_s_oth, 5A_j_s_oth, SP3_i_j_com_1, SP3_i_j_com_2)
    s6Z[hc6Z[n6Z][0]] = 'O';
    s6Z[hc6Z[n6Z][1]] = 'O';
    s6Z[hc6Z[n6Z][2]] = 'O';
    s6Z[hc6Z[n6Z][3]] = 'O';
    if (s6Z[hc6Z[n6Z][4]] == 'C') s6Z[hc6Z[n6Z][4]] = 'B';
    if (s6Z[hc6Z[n6Z][5]] == 'C') s6Z[hc6Z[n6Z][5]] = 'B';

    ++n6Z;
}

void Clusters_Get7K() {    // Detect 7K clusters from 2 5A clusters
    int i, j, j2, k, l, m;
    int scom, sother[2], sp3_com[2], sp3c_i_other, sp3c_j_other;
    int clusSize=7;

    scom=sother[0]=sother[1]=sp3_com[0]=sp3_com[1]=sp3c_i_other=sp3c_j_other=-1;

    for (i=0; i<nsp3c-1; ++i) {  // loop over all 5A_i
        for (j2=3; j2<5; ++j2) {    // loop over both spindles of 5A_i
            for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]; ++j) {  // loop over all 5A_j common with spindle of 5A_i
                if (mem_sp3c[hcsp3c[i][j2]][j]<=i) continue;  // don't find 7K twice

                m=0;
                for (k=3; k<5; k++) {
                    for (l=3; l<5; l++) {
                        if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            if (m>=1) {
                                m++;
                                break;
                            }
                            scom=hcsp3c[i][k];
                            m++;
                        }
                    }
                    if (m>=2) break;
                }
                if (m!=1) continue; // exactly one common spindle between 5A_i and 5A_j

                if (hcsp3c[i][3] == scom) sother[0]=hcsp3c[i][4];
                else sother[0]=hcsp3c[i][3];
                if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] == scom) sother[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                else sother[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];

                m=0;
                for (k=0; k<5; k++) {
                    if (sother[0]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other spindle of 5A_i distinct from whole 5A_j

                m=0;
                for (k=0; k<5; k++) {
                    if (sother[1]==hcsp3c[i][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other spindle of 5A_j distinct from whole 5A_i

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<3; l++) {
                        if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            if (m>=2) {
                                m++;
                                break;
                            }
                            sp3_com[m]=hcsp3c[i][k];
                            m++;
                        }
                    }
                    if (m>=3) break;
                }
                if (m!=2) continue; // exactly two common particles in SP3 rings of 5A_i and 5A_j

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<2; l++) {
                        if (hcsp3c[i][k] == sp3_com[l]) break;
                    }
                    if (l==2) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        sp3c_i_other=hcsp3c[i][k];
                        m++;
                    }
                }
                if (m!=1) continue; // found other uncommon particle from SP3 ring of 5A_i

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<2; l++) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k] == sp3_com[l]) break;
                    }
                    if (l==2) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        sp3c_j_other=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k];
                        m++;
                    }
                }
                if (m!=1) continue; // found other uncommon particle from SP3 ring of 5A_i

                m=0;
                for (k=0; k<5; k++) {
                    if (sp3c_i_other==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other ring of 5A_i distinct from whole 5A_j

                m=0;
                for (k=0; k<5; k++) {
                    if (sp3c_j_other==hcsp3c[i][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other ring of 5A_j distinct from whole 5A_i

                if (n7K==m7K) {
                    hc7K=resize_2D_int(hc7K,m7K,m7K+incrStatic,clusSize,-1);
                    m7K=m7K+incrStatic;
                }
                // Now we have found the 7K cluster

                hc7K[n7K][0]=scom;
                hc7K[n7K][1]=sother[0];
                hc7K[n7K][2]=sother[1];
                hc7K[n7K][3]=sp3_com[0];
                hc7K[n7K][4]=sp3_com[1];
                hc7K[n7K][5]=sp3c_i_other;
                hc7K[n7K][6]=sp3c_j_other;

                quickSort(&hc7K[n7K][1],2);
                quickSort(&hc7K[n7K][3],2);
                quickSort(&hc7K[n7K][5],2);

                Cluster_Write_7K();
            }
        }
    }
}

void Cluster_Write_7K() {
    // hc7K key: (scom, sother, ring_com, ring_other)
    s7K[hc7K[n7K][0]] = 'O';
    s7K[hc7K[n7K][1]] = 'O';
    s7K[hc7K[n7K][2]] = 'O';
    if (s7K[hc7K[n7K][3]] == 'C') s7K[hc7K[n7K][3]] = 'B';
    if (s7K[hc7K[n7K][4]] == 'C') s7K[hc7K[n7K][4]] = 'B';
    if (s7K[hc7K[n7K][5]] == 'C') s7K[hc7K[n7K][5]] = 'B';
    if (s7K[hc7K[n7K][6]] == 'C') s7K[hc7K[n7K][6]] = 'B';

    ++n7K;
}

void Clusters_Get8A_D2d() { // Detect 8A D2d clusters
    int unc[2];
    int com[4];
    int i, j, j2, k, l, m;
    int cnt;
    int flg;
    int break_out;
    int trial[8];
    char errMsg[1000];
    int clusSize=8;
    int *used_sp5b;


    used_sp5b=malloc(nsp5b*sizeof(int)); if (used_sp5b==NULL) { sprintf(errMsg,"Clusters_Get8A_D2d(): used_sp5b[] malloc out of memory\n");  Error(errMsg); }

    for (i=0; i<nsp5b-1; ++i) {  // loop over all sp5b_i
        memset(used_sp5b, 0, nsp5b*sizeof(*used_sp5b));
        used_sp5b[i]=1;
        for (j2=0; j2<5; ++j2) {
            for (j=0; j<nmem_sp5b[hcsp5b[i][j2]]; ++j) {  // loop over all sp5b_j
                if (mem_sp5b[hcsp5b[i][j2]][j]<=i) continue;
                if (used_sp5b[mem_sp5b[hcsp5b[i][j2]][j]]==1) continue;
                used_sp5b[mem_sp5b[hcsp5b[i][j2]][j]]=1;
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if(hcsp5b[i][k] == hcsp5b[mem_sp5b[hcsp5b[i][j2]][j]][l]) {
                            if (m<5) com[m]=hcsp5b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue; // exactly four members of the SP5 rings of sp5b_i and sp5b_j in common

                if (hcsp5b[i][5] == hcsp5b[mem_sp5b[hcsp5b[i][j2]][j]][5]) continue;  // distinct spindles

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5b[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=hcsp5b[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5b[mem_sp5b[hcsp5b[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=hcsp5b[mem_sp5b[hcsp5b[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=hcsp5b[i][5];    // build up trial cluster
                trial[1]=hcsp5b[mem_sp5b[hcsp5b[i][j2]][j]][5];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[4])==1 && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=hcsp5b[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[2])==1 && hcsp5b[i][k]!=trial[2] && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        trial[6]=hcsp5b[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[3])==1 && hcsp5b[i][k]!=trial[3] && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        trial[7]=hcsp5b[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A][k]=trial[k];

                    Cluster_Write_8A();
                }
            }
        }
    }
    for (i=0; i<nsp5c - 1; ++i) {    // loop over all 7A_i
        for (j2=5; j2<6; ++j2) {
            for (j=0; j<nmem_sp5c[hcsp5c[i][j2]]; ++j) {  // loop over all 7A_j
                if (mem_sp5c[hcsp5c[i][j2]][j]<=i) continue;
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if (hcsp5c[i][k] == hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][l]) {
                            if (m<5) com[m]=hcsp5c[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue;     // exactly four members of the SP5 rings of 7A_i and 7A_j in common

                flg = hcsp5c[i][5] == hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][5] && hcsp5c[i][6] == hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][6];
                flg = flg || (hcsp5c[i][5] == hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][6] && hcsp5c[i][6] == hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][5]);
                if (flg!=1) continue; // both spindles common

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5c[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=hcsp5c[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=hcsp5c[mem_sp5c[hcsp5c[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=hcsp5c[i][5];    // build up trial cluster
                trial[1]=hcsp5c[i][6];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[i][k],trial[4])==1 && hcsp5c[i][k]!=trial[4] && hcsp5c[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=hcsp5c[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[i][k],trial[2])==1 && hcsp5c[i][k]!=trial[2] && hcsp5c[i][k]!=trial[4] && hcsp5c[i][k]!=trial[5]) {
                        trial[6]=hcsp5c[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[i][k],trial[3])==1 && hcsp5c[i][k]!=trial[3] && hcsp5c[i][k]!=trial[4] && hcsp5c[i][k]!=trial[5]) {
                        trial[7]=hcsp5c[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A][k]=trial[k];

                    Cluster_Write_8A();
                }
            }
        }
    }
    for (i=0; i<nsp5b; ++i) {    // loop over all sp5b_i
        for (j2=5; j2<6; ++j2) {
            for (j=0; j<nmem_sp5c[hcsp5b[i][j2]]; ++j) {  // loop over all 7A_j
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if (hcsp5b[i][k] == hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][l]) {
                            if (m<5) com[m]=hcsp5b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue; // exactly four members of the SP5 rings of sp5b_i and 7A_j in common

                flg = hcsp5b[i][5] == hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][5] || hcsp5b[i][5] == hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][6];
                if (flg!=1) continue;   // sp5b_i spindle common with one of 7A_j spindles

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5b[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=hcsp5b[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][5]; // build up trial cluster
                trial[1]=hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][6];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[4])==1 && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=hcsp5b[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[2])==1 && hcsp5b[i][k]!=trial[2] && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        trial[6]=hcsp5b[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[3])==1 && hcsp5b[i][k]!=trial[3] && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        trial[7]=hcsp5b[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A][k]=trial[k];
                    Cluster_Write_8A();
                }
            }
        }
    }

    free(used_sp5b);
}

void Cluster_Write_8A() {// hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
    if (s8A[hc8A[n8A][0]] == 'C') s8A[hc8A[n8A][0]] = 'B';
    if (s8A[hc8A[n8A][1]] == 'C') s8A[hc8A[n8A][1]] = 'B';
    if (s8A[hc8A[n8A][2]] == 'C') s8A[hc8A[n8A][2]] = 'B';
    if (s8A[hc8A[n8A][3]] == 'C') s8A[hc8A[n8A][3]] = 'B';
    if (s8A[hc8A[n8A][4]] == 'C') s8A[hc8A[n8A][4]] = 'B';
    if (s8A[hc8A[n8A][5]] == 'C') s8A[hc8A[n8A][5]] = 'B';
    s8A[hc8A[n8A][6]] = 'O';
    s8A[hc8A[n8A][7]] = 'O';

    ++n8A;
}

void Clusters_Get8B_Cs() { // Detect 8B Cs clusters
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

    for (j=0; j < cnb[primary_spindle]; ++j) { // loop over all j particles bonded to first spindle of 7A_i
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

void Clusters_Get8K() {    // Detect 8K clusters
    int i, j, j2, k, l, m, n;
    int cp[2], unc[3], scom, sother[2];
    int clusSize=8;

    cp[0]=cp[1]=unc[0]=unc[1]=unc[2]=scom=sother[0]=sother[1]=-1;

    for (i=0; i<nsp3c-2; ++i) {  // loop over all sp3c_i
        for (j2=0; j2<3; j2++) {    // loop over all particles in SP3 ring of sp3c_i
            for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]-1; ++j) { // loop over all sp3c_j which sp3c[i][j2] is a member of
                if (mem_sp3c[hcsp3c[i][j2]][j]<=i) continue;  // check not used mem_sp3c[sp3c[i][j2]][j] before

                m = 0;  // check j2 from SP3 ring of sp3c_i is in SP3 ring of mem_sp3c[sp3c[i][j2]][j]
                for(k=0; k<3; ++k) {
                    if (hcsp3c[i][j2]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][k]) {
                        cp[0]=hcsp3c[i][j2];
                        m++;
                    }
                }
                if (m!=1) continue;

                m = 0;  // find extra 1 particle from SP3 ring of sp3c_i which is also in SP3 ring of mem_sp3c[sp3c[i][j2]][j]
                for(k=0; k<3; ++k) {
                    if (k==j2) continue;    // don't check j2 again
                    if (hcsp3c[i][k]<hcsp3c[i][j2]) continue; // will have found before or do not find after when using different j2 particle
                    for(l=0; l<3; ++l) {
                        if (hcsp3c[i][k]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            cp[1]=hcsp3c[i][k];
                            m++;
                        }
                    }
                }
                if (m!=1) continue;

                m = 0;  // check exactly one common spindle between sp3c_i and mem_sp3c[sp3c[i][j2]][j]
                for(k=3; k<5; ++k) {
                    for(l=3; l<5; ++l) {
                        if (hcsp3c[i][k]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                            scom=hcsp3c[i][k];
                            m++;
                        }
                    }
                }
                if (m!=1) continue;

                if (hcsp3c[i][3]==scom) sother[0]=hcsp3c[i][4];
                else sother[0]=hcsp3c[i][3];
                if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]==scom) sother[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                else sother[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];

                k=0;    // find particle from SP3 ring of sp3c_i which is common with 5A mem_sp3c[sp3c[i][j2]][j]
                for(l=0; l<3; ++l) {
                    if (hcsp3c[i][l]==cp[0] || hcsp3c[i][l]==cp[1]) continue;
                    for(m=0; m<5; ++m) {
                        if (hcsp3c[i][l]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][m]) break;
                    }
                    if (m==5) {
                        if (k==1) {
                            k++;
                            break;
                        }
                        unc[0]=hcsp3c[i][l];
                        k++;
                    }
                }
                if (k!=1) continue;

                k=0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                for(l=0; l<3; ++l) {
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]==cp[0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]==cp[1]) continue;
                    for(m=0; m<5; ++m) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]==hcsp3c[i][m]) break;
                    }
                    if (m==5) {
                        if (k==1) {
                            k++;
                            break;
                        }
                        unc[1]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                        k++;
                    }
                }
                if (k!=1) continue;

                for (k=j+1; k<nmem_sp3c[hcsp3c[i][j2]]; ++k) {
                    if (mem_sp3c[hcsp3c[i][j2]][k]<=i) continue;  // higher index of sp3c cluster than i
                    if (mem_sp3c[hcsp3c[i][j2]][k]<=mem_sp3c[hcsp3c[i][j2]][j]) continue;   // higher index of sp3c cluster than mem_sp3c[sp3c[i][j2]][j]

                    n = 0;  // check common SP3 ring particles are exactly cp
                    for(l=0; l<3; ++l) {
                        for(m=0; m<2; ++m) {
                            if (cp[m]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]) {
                                n++;
                            }
                        }
                    }
                    if (n!=2) continue;

                    n = 0;  // check spindles are exactly sother
                    for(l=3; l<5; ++l) {
                        for(m=0; m<2; ++m) {
                            if (sother[m]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]) {
                                n++;
                            }
                        }
                    }
                    if (n!=2) continue;

                    n=0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                    for(l=0; l<3; ++l) {
                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]==cp[0] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]==cp[1]) continue;
                        for(m=0; m<5; ++m) {
                            if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]==hcsp3c[i][m] || hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l]==hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) break;
                        }
                        if (m==5) {
                            if (n==1) {
                                n++;
                                break;
                            }
                            unc[2]=hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                            n++;
                        }
                    }
                    if (n!=1) continue;

                    // Now we have found the 8K cluster
                    if (n8K==m8K) {
                        hc8K=resize_2D_int(hc8K,m8K,m8K+incrStatic,clusSize,-1);
                        m8K=m8K+incrStatic;
                    }
                    // hc8K key: (SP3_common_1, SP3_common_2, spindle_1, spindle_2, spindle_3, other_SP3_1, other_SP3_2, other_SP3_3)
                    hc8K[n8K][0]=cp[0];
                    hc8K[n8K][1]=cp[1];
                    hc8K[n8K][2]=scom;
                    hc8K[n8K][3]=sother[0];
                    hc8K[n8K][4]=sother[1];
                    hc8K[n8K][5]=unc[0];
                    hc8K[n8K][6]=unc[1];
                    hc8K[n8K][7]=unc[2];

                    quickSort(&hc8K[n8K][0],2);
                    quickSort(&hc8K[n8K][2],3);
                    quickSort(&hc8K[n8K][5],3);

                    Cluster_Write_8K();
                }
            }
        }
    }
}

void Cluster_Write_8K() {
    if(s8K[hc8K[n8K][2]] == 'C') s8K[hc8K[n8K][2]] = 'B';
    if(s8K[hc8K[n8K][3]] == 'C') s8K[hc8K[n8K][3]] = 'B';
    if(s8K[hc8K[n8K][4]] == 'C') s8K[hc8K[n8K][4]] = 'B';
    if(s8K[hc8K[n8K][5]] == 'C') s8K[hc8K[n8K][5]] = 'B';
    if(s8K[hc8K[n8K][6]] == 'C') s8K[hc8K[n8K][6]] = 'B';
    if(s8K[hc8K[n8K][7]] == 'C') s8K[hc8K[n8K][7]] = 'B';
    s8K[hc8K[n8K][0]] = 'O';
    s8K[hc8K[n8K][1]] = 'O';
    ++n8K;
}

void Clusters_Get9A_D3h() {    // Detect 9A D3h clusters
    int i, j, j2, k, l, m, n;
    int db[2], ob[4];
    int flg;
    int clusSize=9;

    for (i=0; i<nsp4b-2; ++i) {  // loop over all sp4b_i
        for (j2=0; j2<1; j2++) {
            for (j=0; j<nmem_sp4b[hcsp4b[i][j2]]; ++j) { // loop over all sp4b_j
                if (mem_sp4b[hcsp4b[i][j2]][j]<=i) continue;
                if (hcsp4b[i][4] == hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4]) continue;  // if spindles common continue
                if (Bonds_BondCheck(hcsp4b[i][4], hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4])) continue;   // if spindles bonded continue
                m = 0;
                for(k=0; k<4; ++k) {
                    for(l=0; l<4; ++l) {
                        if(hcsp4b[i][k] == hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l]){
                            if(m<2) db[m] = hcsp4b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=2) continue;         // two common particles between SP4 rings of sp4b_i and sp4b_j

                m = 0;
                for (k=0; k<4; ++k) {
                    if(hcsp4b[i][k] == db[0] || hcsp4b[i][k] == db[1]) continue;  // find particles in SP4 ring of sp4b_i not common to SP4 ring of sp4b_j
                    for(l=0; l<4; ++l){
                        if(hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l] == db[0] || hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l] == db[1]) continue;    // find particles in SP4 ring of sp4b_j not common to SP4 ring of sp4b_i
                        if(Bonds_BondCheck(hcsp4b[i][k], hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l])) {    // check non-common SP4 ring particles from sp4b_i and sp4b_j are bonded
                            if(m<4) ob[m] = hcsp4b[i][k];
                            ++m;
                            if(m<4) ob[m] = hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue;
                // POSSIBLE IMPROVEMENT!! could make detection faster here by looping over sp4b clusters bonded to uncommon particles to sp4b_i
                for(k=mem_sp4b[hcsp4b[i][j2]][j]+1; k<nsp4b; k++) {    // loop over all sp4b_k
                    // ERROR!! need to check spindle of sp4b_k distinct from spindles of sp4b_i and sp4b_j and no bonds between these three particles
                    n = 0;
                    for(l=0; l<4; ++l){
                        for(m=0; m<4; ++m){
                            if(hcsp4b[k][l] == ob[m]){
                                ++n;
                                break;
                            }
                        }
                    }
                    if (n != 4) continue;

                    // Now we have found the 9A D3h cluster
                    if (n9A==m9A) {
                        hc9A=resize_2D_int(hc9A,m9A,m9A+incrStatic,clusSize,-1);
                            m9A=m9A+incrStatic;
                        }
                        if (hcsp4b[i][4]<hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4] && hcsp4b[i][4]<hcsp4b[k][4]) {
                            hc9A[n9A][6]=hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                            hc9A[n9A][0]=hcsp4b[i][0];
                            hc9A[n9A][1]=hcsp4b[i][1];
                            hc9A[n9A][2]=hcsp4b[i][2];
                            hc9A[n9A][3]=hcsp4b[i][3];

                            if (Bonds_BondCheck(hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4],hc9A[n9A][0])) {
                                hc9A[n9A][7]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                hc9A[n9A][8]=hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                            }
                            else {
                                hc9A[n9A][7]=hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                hc9A[n9A][8]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                            }

                            for (l=0;l<4;l++) {
                                flg=1;
                                for (m=0;m<4;m++) {
                                    if (ob[l]==hc9A[n9A][m]) {
                                        flg=0;
                                        break;
                                    }
                                }
                                if (flg==1 && Bonds_BondCheck(ob[l],hc9A[n9A][0])) {
                                    hc9A[n9A][4]=ob[l];
                                }
                                else if (flg==1) {
                                    hc9A[n9A][5]=ob[l];
                                }
                            }
                        }

                        else if (hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4]<hcsp4b[i][4] && hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4]<hcsp4b[k][4]) {
                            hc9A[n9A][6]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                            hc9A[n9A][0]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][0];
                            hc9A[n9A][1]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][1];
                            hc9A[n9A][2]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][2];
                            hc9A[n9A][3]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][3];

                            if (Bonds_BondCheck(hcsp4b[i][4],hc9A[n9A][0])) {
                                hc9A[n9A][7]=hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                hc9A[n9A][8]=hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                            }
                            else {
                                hc9A[n9A][7]=hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                                hc9A[n9A][8]=hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                            }

                            for (l=0;l<4;l++) {
                                flg=1;
                                for (m=0;m<4;m++) {
                                    if (ob[l]==hc9A[n9A][m]) {
                                        flg=0;
                                        break;
                                    }
                                }
                                if (flg==1 && Bonds_BondCheck(ob[l],hc9A[n9A][0])) {
                                    hc9A[n9A][4]=ob[l];
                                }
                                else if (flg==1) {
                                    hc9A[n9A][5]=ob[l];
                                }
                            }
                        }

                        else {
                            hc9A[n9A][6]=hcsp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[k][.] cluster
                            hc9A[n9A][0]=hcsp4b[k][0];
                            hc9A[n9A][1]=hcsp4b[k][1];
                            hc9A[n9A][2]=hcsp4b[k][2];
                            hc9A[n9A][3]=hcsp4b[k][3];

                            if (Bonds_BondCheck(hcsp4b[i][4],hc9A[n9A][0])) {
                                hc9A[n9A][7]=hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                                hc9A[n9A][8]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                            }
                            else {
                                hc9A[n9A][7]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[j][.] cluster
                                hc9A[n9A][8]=hcsp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[i][.] cluster
                            }

                            if (Bonds_BondCheck(db[0],hc9A[n9A][0])) {
                                hc9A[n9A][4]=db[0];
                                hc9A[n9A][5]=db[1];
                            }
                            else {
                                hc9A[n9A][4]=db[1];
                                hc9A[n9A][5]=db[0];
                            }
                        }
                    Cluster_Write_9A();
                    break;
                }
            }
        }
    }
}

void Cluster_Write_9A() {
    // hc9A key: (SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_to_0_in_SP4_lowest_s, SP4_to_1_in_SP4_lowest_s, s_lowest, s_to_0_in_SP4_lowest_s,s_to_2_in_SP4_lowest_s)
    if(s9A[hc9A[n9A][0]] == 'C') s9A[hc9A[n9A][0]] = 'B';
    if(s9A[hc9A[n9A][1]] == 'C') s9A[hc9A[n9A][1]] = 'B';
    if(s9A[hc9A[n9A][3]] == 'C') s9A[hc9A[n9A][3]] = 'B';
    if(s9A[hc9A[n9A][4]] == 'C') s9A[hc9A[n9A][4]] = 'B';
    if(s9A[hc9A[n9A][6]] == 'C') s9A[hc9A[n9A][6]] = 'B';
    if(s9A[hc9A[n9A][7]] == 'C') s9A[hc9A[n9A][7]] = 'B';
    s9A[hc9A[n9A][2]] = 'O';
    s9A[hc9A[n9A][5]] = 'O';
    s9A[hc9A[n9A][8]] = 'O';
    ++n9A;
}

void Clusters_Get9B_10B_11B_11E_12D() {    // Detect 9B, 10B, 11A, 11E & 12D
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

            if (do10B==1) Clusters_Get10B_C3v(i, j);
            if (do11B==1) {
                if (Clusters_Get11B_C2v()) {
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

void Clusters_Get10B_C3v(int i, int j) {        // Return 1 if 9B is also 10B cluster
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

int Clusters_Get11B_C2v() { // Detect 11B C2v clusters
    //  11B is very similar to 11C & 11D
    // Call from 9B, the central particle must have 11 particles bonded to it.
    // The 10th & 11th particles are bonded to each other. They are both bonded
    // to 2 other shell particles which aren't bonded to each other. These
    // bonded shell particles are distict i.e. there are 4 of
    // them. The extra 2 particles form 4 sp4 rings in the shell.
    int k, l, m;
    int b1[2], b2[2], nb1, nb2;
    int ep[2]; // The two extra particles
    int flg11, flg12, flg21, flg22;
    int break_out;
    int clusSize=11;

    if(cnb[hc9B[n9B][8]]!= 10) return 0; // s_com has 10 bonds in total (all forming the shell)

    m = 0;
    break_out=0;
    for(k=0; k<10; ++k) {
        for(l=0; l<8; ++l) {
            if(bNums[hc9B[n9B][8]][k] == hc9B[n9B][l]) break;
        }
        if(l==8){
            if(m==2) {
                break_out=1;
                break;
            }
            ep[m++] = bNums[hc9B[n9B][8]][k];    // two extra particles
        }
    }
    if(break_out==1 || m<2) return 0;

    if(Bonds_BondCheck(ep[0], ep[1])==0) return 0;  // extra particles must be bonded
    nb1 = nb2 = 0;
    for(k=0; k<8; ++k){
        if(Bonds_BondCheck(hc9B[n9B][k], ep[0])){
            if(nb1 == 2) return 0;  // extra particle 1 bonded to 2 members of 9B shell
            b1[nb1]=hc9B[n9B][k];
            nb1++;
        }
        if(Bonds_BondCheck(hc9B[n9B][k], ep[1])){
            if(nb2 == 2) return 0;  // extra particle 2 bonded to 2 members of 9B shell
            b2[nb2]=hc9B[n9B][k];
            nb2++;
        }
    }
    if(nb1 != 2 || nb2 != 2) return 0;  // extra particles 1 and 2 bonded to 2 members of 9B shell

    flg11 = b1[0] == b2[0] || b1[0] == b2[1]; // Particles bonded to extra 2 particles b[]
    flg11 = flg11 || b1[1] == b2[0] || b1[1] == b2[1]; // must be distinct.
    if(flg11) return 0;
    flg11 = Bonds_BondCheck(b1[0], b1[1]); // paritcles b1[] mustn't be bonded
    flg22 = Bonds_BondCheck(b2[0], b2[1]);
    if(flg11 || flg22) return 0;
    flg11 = Bonds_BondCheck(b1[0], b2[0]);
    flg12 = Bonds_BondCheck(b1[0], b2[1]);
    flg21 = Bonds_BondCheck(b1[1], b2[0]);
    flg22 = Bonds_BondCheck(b1[1], b2[1]);
    if(!((flg11 && !flg12) || (!flg11 && flg12))) return 0;
    if(!((flg21 && !flg22) || (!flg21 && flg22))) return 0;

    if(n11B == m11B) {
        hc11B=resize_2D_int(hc11B,m11B,m11B+incrStatic,clusSize,-1);
        m11B=m11B+incrStatic;
    }
    hc11B[n11B][0] = hc9B[n9B][0];
    hc11B[n11B][1] = hc9B[n9B][1];
    hc11B[n11B][2] = hc9B[n9B][2];
    hc11B[n11B][3] = hc9B[n9B][3];
    hc11B[n11B][4] = hc9B[n9B][4];
    hc11B[n11B][5] = hc9B[n9B][5];
    hc11B[n11B][6] = hc9B[n9B][6];
    hc11B[n11B][7] = hc9B[n9B][7];
    hc11B[n11B][8] = hc9B[n9B][8];
    if (b1[0]==1) {
        hc11B[n11B][9] = ep[0];
        hc11B[n11B][10] = ep[1];
    }
    else {
        hc11B[n11B][9] = ep[1];
        hc11B[n11B][10] = ep[0];
    }
    Cluster_Write_11B();

    return 1;
}

void Cluster_Write_11B() {
    // hc11B key: (as 9B, ep_1_to_9B_0, ep_2_to_9B_1)
    int i;
    for(i=0; i<6; i++) {
        if (s11B[hc11B[n11B][i]] == 'C') s11B[hc11B[n11B][i]] = 'B';
    }
    if(s11B[hc11B[n11B][6]] != 'S') s11B[hc11B[n11B][6]] = 'O';
    if(s11B[hc11B[n11B][7]] != 'S') s11B[hc11B[n11B][7]] = 'O';
    s11B[hc11B[n11B][8]] = 'S';
    if(s11B[hc11B[n11B][9]] == 'C') s11B[hc11B[n11B][9]] = 'B';
    if(s11B[hc11B[n11B][10]] == 'C') s11B[hc11B[n11B][10]] = 'B';
}

void Clusters_Get11W() {
    // 10B with 1 extra particle. The central particle of 10B must have exactly 10 particles bonded to it.
    // The extra particle is bonded to 10B central particle but not bonded to three shell spindles of 7A in 10B

    int id_10B, spindle_10B;
    int extra_particle;

    for(id_10B=0;id_10B<n10B; id_10B++) {
        spindle_10B = hc10B[id_10B][9];
        if (cnb[spindle_10B] == 10) {   // s_com has 10 bonds in total (all forming the shell)

            extra_particle = get_11W_extra_particle(id_10B, spindle_10B);

            // extra particle must not be bonded to three 7A spindles in shell of 10B
            if (is_particle_bonded_to_7As(id_10B, extra_particle)) continue;

            resize_hc11W();
            populate_hc11W(id_10B, extra_particle);
            populate_s11W();
            n11W++;
        }
    }
}

void populate_hc11W(int id_10B, int extra_particle) {
    int i;
    for (i = 0; i < 10; i++){
                hc11W[n11W][i] = hc10B[id_10B][i];
            }
    hc11W[n11W][10] = extra_particle;
}

void resize_hc11W() {
    int clusSize=11;

    if (n11W == m11W) {
        hc11W = resize_2D_int(hc11W, m11W, m11W + incrStatic, clusSize, -1);
        m11W = m11W + incrStatic;
    }
}

int is_particle_bonded_to_7As(int id_10B, int extra_particle) {
    int i;

    for(i=6; i < 9; i++) {
        if (Bonds_BondCheck(extra_particle, hc10B[id_10B][i]) == 1){
            return 1;
        }
    }
    return 0;
}

int get_11W_extra_particle(int id_10B, int spindle_10B) {
    int i;
    for (i = 0; i < 10; ++i) {
        if (is_particle_in_10B(bNums[spindle_10B][i], id_10B) == 0) {
            return bNums[spindle_10B][i];
        }
    }
}

int is_particle_in_10B(int particle_id, int id_10B) {
    // Returns 1 if particle_id is in id_10B else returns 0
    int i;

    for (i=0; i<9; i++) {
        if (particle_id == hc10B[id_10B][i]) {
            return 1;
        }
    }
    return 0;
}

void populate_s11W() {
    int i;

    for(i=0; i<10; i++) {
        if (s11W[hc11W[n11W][i]] == 'C') s11W[hc11W[n11W][i]] = 'B';
    }
    s11W[hc11W[n11W][10]] = 'O';
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
                if (do12D==1) n12D += Clusters_Get12D_D2d(j, k, sp2i, hcsp5c[k][6]);

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
                if (do12D==1) n12D += Clusters_Get12D_D2d(j, k, sp2i, hcsp5c[k][5]);
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
                if (do12D==1) n12D += Clusters_Get12D_D2d(j, k, sp2j, hcsp5c[k][6]);
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
                if (do12D==1) n12D += Clusters_Get12D_D2d(j, k, sp2j, hcsp5c[k][5]);
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

int Clusters_Get12D_D2d(int j, int k, int sp1, int sp2) {  // Return 1 if 12B is also 11E
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
            break_out=0;
            for (n=0; n<5; n++) {
                if (hcsp5c[l][n]==hc11E[n11E][0]) {
                    flg1=1;
                    continue;
                }
                if (hcsp5c[l][n]==hc11E[n11E][1]) {
                    flg2=1;
                    continue;
                }
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
        if (cnb[id_9K_common] < 10) {
            num_extra_particles = 0;
            for (bonded_to_spindle_id = 0; bonded_to_spindle_id < cnb[id_9K_common]; bonded_to_spindle_id++) {
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

void Clusters_Get10A_C3v() { // Detect 10A D4d clusters
    int i, j, j2, k, l, m;
    char errMsg[1000];
    int clusSize=10;
    int *used_sp4b;

    used_sp4b=malloc(nsp4b*sizeof(int)); if (used_sp4b==NULL) { sprintf(errMsg,"Clusters_Get10A_C3v(): used_sp4b[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i<nsp4b; ++i) used_sp4b[i] = 0;

    for (i=0; i<nsp4b-1; ++i) {  // loop over all sp4b_i
        for (j2=0; j2<nsp4b; ++j2) used_sp4b[j2] = 0;
        used_sp4b[i]=1;
        for (j2=0; j2<cnb[hcsp4b[i][0]]; ++j2) {
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

void Clusters_Get10W() { // Detect 10W clusters
    int i, j, k, l, m;
    int sp5b_clusts[5], shell_parts[9];
    int clusSize=10;

    sp5b_clusts[0]=sp5b_clusts[1]=sp5b_clusts[2]=sp5b_clusts[3]=sp5b_clusts[4]=-1;


    for (i=0; i<nsp5b; ++i) { // loop over all sp5b
        if (cnb[hcsp5b[i][5]]!=9) continue;   // central particle must have coordination number 9

        k=0;    // find 5 other sp5b's with spindle in common with sp5b_i
        for (j=0; j<nmem_sp5b[hcsp5b[i][5]]; ++j) { // note check that spindle of sp5b_i and sp5b_j must be common by later check
            if (mem_sp5b[hcsp5b[i][5]][j]<=i) continue;   // i for sp5b must be lowest of all sp5b indices
            // ERROR !! need to check that spindle of sp5b_j is spindle of sp5b_i
            if (k>=5) {
                k++;
                break;
            }
            sp5b_clusts[k]=mem_sp5b[hcsp5b[i][5]][j];
            k++;
        }
        if (k!=5) continue; // not correct number of sp5b clusters
        // now found exactly 5 sp5b clusters common to spindle of sp5b_i
        for (j=0; j<5; j++) {
            shell_parts[j]=hcsp5b[i][j];
        }

        m=5;
        for (j=0; j<5; j++) {
            for (k=0; k<5; k++) {
                for (l=0; l<m; l++) {
                    if (shell_parts[l]==hcsp5b[mem_sp5b[hcsp5b[i][5]][j]][k]) break;
                }
                if (l==m) {
                    if (m>=9) {
                        m++;
                        break;
                    }
                    shell_parts[m]=hcsp5b[mem_sp5b[hcsp5b[i][5]][j]][k];
                    m++;
                }
            }
            if (m>=10) break;
        }
        if (m!=9) continue; // not all coordination shell particles of sp5b[i][5] are in the SP5 rings of the 5xsp5b clusters we found

        if (n10W == m10W) {
            hc10W=resize_2D_int(hc10W,m10W,m10W+incrStatic,clusSize,-1);
            m10W=m10W+incrStatic;
        }
        // hc10W key: (sp5bs_common_central_spindle_particle, sp5bs_SP5_ring_shell_particles)
        hc10W[n10W][0] = hcsp5b[i][5];
        for (j=0; j<9; j++) hc10W[n10W][j+1]=shell_parts[j];
        quickSort(&hc10W[n10W][1],9);
        Cluster_Write_10W();
    }
}

void Cluster_Write_10W() {
    int i;

    for(i=1; i<10; i++) {
        if (s10W[hc10W[n10W][i]] == 'C') s10W[hc10W[n10W][i]] = 'B';
    }
    s10W[hc10W[n10W][0]] = 'S';

    ++n10W;
}

void Clusters_Get11A() {
    // Detect 11A clusters

    int first_6A_id, second_6A_id, k, l, m;
    int scom, sother[2];

    scom=sother[0]=sother[1]=-1;

    for(first_6A_id=0; first_6A_id<nsp4c-1; ++first_6A_id){
        // POSSIBLE IMPROVEMENT: loop over all sp4c clusters for spindles of sp4c_i
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
        s11A[hc11A[n11A][i]] = 'O';
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

void Clusters_Get12K() {
    // 12K clusters are an 11A with an extra particle bonded to three of the ring particles of the 11A
    // Since there are 8 possible places to attach an extra particle, there may be more than 1 12K for each 11A.

    int ptr_11A, ring_number;
    int sp3_rings[8][3] = {-1};
    int *sp3_ring;

    // Loop over all 11A clusters
    for(ptr_11A=0; ptr_11A<n11A; ptr_11A++) {

        // Find the particle IDs of the 8 sp3 rings made from the rings of 11A
        get_12K_ring_bonds(ptr_11A, sp3_rings);

        // Loop through the 8 sp3 rings made by the sp4s of 11A and see if a particle is attached to any of them
        for(ring_number=0; ring_number<8; ring_number++) {
            sp3_ring = sp3_rings[ring_number];
            find_12K_cluster(ptr_11A, sp3_ring);
        }
    }
}

void find_12K_cluster(int ptr_11A, const int *sp3_ring) {
    // Check if there is a particle attached to the sp3 ring of an 11A, if there is write a 12K
    // It is possible that there are two particles attached to the ring particles in which case we ignore this cluster
    int i, num_attached_particles;
    int particle_id, ep = {-1};

    num_attached_particles = 0;
    // loop through all particles bonded to the first sp3 ring particle
    for (i = 0; i < cnb[sp3_ring[0]]; i++) {
        particle_id = bNums[sp3_ring[0]][i];
        if (Bonds_BondCheck(sp3_ring[1], particle_id) == 1) {
            if (Bonds_BondCheck(sp3_ring[2], particle_id) == 1) {
                if (is_particle_in_11A(ptr_11A, particle_id) == 0) {
                    num_attached_particles++;
                    ep = particle_id;
                }
            }
        }
    }
    if (num_attached_particles == 1) {
        Cluster_Write_12K(ep, ptr_11A);
    }
}

int is_particle_in_11A(int id_11A, int id_particle) {
    // Returns 0 if particle is not in specified 11A
    int i;
    for (i = 0; i < 11; i++) {
        if (id_particle == hc11A[id_11A][i]) {
            return 1;
        }
    }
    return 0;
}

void get_12K_ring_bonds(int id_11A, int (*sp3_rings)[3]) {

    int m;
    int particle_1, particle_2;

    for(particle_1=0; particle_1<4; particle_1++) {
        m = 0;
        sp3_rings[particle_1][m] = hc11A[id_11A][particle_1];
        for(particle_2=4; particle_2<8; particle_2++) {
            if (Bonds_BondCheck(hc11A[id_11A][particle_1], hc11A[id_11A][particle_2])) {
                m++;
                sp3_rings[particle_1][m] = hc11A[id_11A][particle_2];
            }
        }
    }
    for(particle_1=4; particle_1<8; particle_1++) {
        m = 0;
        sp3_rings[particle_1][m] = hc11A[id_11A][particle_1];
        for(particle_2=0; particle_2<4; particle_2++) {
            if (Bonds_BondCheck(hc11A[id_11A][particle_1], hc11A[id_11A][particle_2])) {
                m++;
                sp3_rings[particle_1][m] = hc11A[id_11A][particle_2];
            }
        }
    }
}

void Cluster_Write_12K(int ep, int id_11A) {

    int i;
    int clusSize=12;

    if(n12K == m12K) {
        hc12K=resize_2D_int(hc12K,m12K,m12K+incrStatic,clusSize,-1);
        m12K=m12K+incrStatic;
    }
    // hc12K key: (SP4 going up, sd going up, scom, ep)

    for (i=0; i<11; i++) hc12K[n12K][i] = hc11A[id_11A][i];
    hc12K[n12K][11]=ep;

    for (i=0; i<8; i++) {
        if (s12K[hc12K[n12K][i]] == 'C') s12K[hc12K[n12K][i]] = 'B';
    }
    if(s12K[hc12K[n12K][8]] != 'S') s12K[hc12K[n12K][8]] = 'O';
    if(s12K[hc12K[n12K][9]] != 'S') s12K[hc12K[n12K][9]] = 'O';
    s12K[hc12K[n12K][10]] = 'S';
    if(s12K[hc12K[n12K][11]] == 'C') s12K[hc12K[n12K][11]] = 'B';
    n12K++;
}

void Clusters_Get11C() {
    int ar[2],uncommon_spindle[2];
    int id_first_7A, id_second7A, k, l, m, ncom, common_spindle;
    int break_out;

    int first_spindle_id, first_spindle_pointer, second_7A_pointer;

    ar[0]=ar[1]=uncommon_spindle[0]=uncommon_spindle[1]=common_spindle=-1;

    for (id_first_7A=0; id_first_7A<nsp5c-1; ++id_first_7A) {
        for(first_spindle_pointer=5; first_spindle_pointer<7; first_spindle_pointer++) {
            first_spindle_id = hcsp5c[id_first_7A][first_spindle_pointer];
            for (second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[first_spindle_id]; second_7A_pointer++) {
                id_second7A = mem_sp5c[first_spindle_id][second_7A_pointer];
                if(id_second7A<=id_first_7A) continue; // Dont detect the same cluster twice!

                if (get_11C_spindle_particles(uncommon_spindle, id_first_7A, id_second7A, &common_spindle) == 0) continue;

                ncom = 0;
                // need two common particles from SP5 rings

                ncom = get_bonded_7A_ring_particles(ar, id_first_7A, id_second7A, ncom);

                if (ncom != 2) continue;
                if (Bonds_BondCheck(ar[0], ar[1]) != 1) continue;

                // two common SP5 ring particles are bonded

                ncom = 0;
                for (k = 0; k < 5; ++k) {
                    if (hcsp5c[id_first_7A][k] == ar[0] || hcsp5c[id_first_7A][k] == ar[1]) continue;
                    for (l = 0; l < 5; ++l) {
                        if (hcsp5c[id_second7A][l] == ar[0] || hcsp5c[id_second7A][l] == ar[1]) continue;
                        if (Bonds_BondCheck(hcsp5c[id_first_7A][k], hcsp5c[id_second7A][l])) ++ncom;
                    }
                }
                if (ncom != 2) continue;

                // two bonds between non-common SP5 ring particles
                resize_hc11C();

                // hc11C key: (s_com, s_i, s_j, r_ca, r_cb, d_i, d_i, d_j, d_j, unc_i, unc_j)
                hc11C[n11C][0] = common_spindle;
                hc11C[n11C][1] = uncommon_spindle[0];
                hc11C[n11C][2] = uncommon_spindle[1];
                hc11C[n11C][3] = ar[0];
                hc11C[n11C][4] = ar[1];

                l = 5;
                m = 7;
                break_out = 0;
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_first_7A][k], ar[0]) && hcsp5c[id_first_7A][k] != ar[1]) {
                        if (l == 7) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][l] = hcsp5c[id_first_7A][k];
                        l++;
                    }
                }
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_first_7A][k], ar[1]) && hcsp5c[id_first_7A][k] != ar[0]) {
                        if (l == 7) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][l] = hcsp5c[id_first_7A][k];
                        l++;
                    }
                }
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_second7A][k], ar[0]) && hcsp5c[id_second7A][k] != ar[1]) {
                        if (m == 9) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][m] = hcsp5c[id_second7A][k];
                        m++;
                    }
                }
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_second7A][k], ar[1]) && hcsp5c[id_second7A][k] != ar[0]) {
                        if (m == 9) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][m] = hcsp5c[id_second7A][k];
                        m++;
                    }
                }
                if (break_out == 1 || l < 7 || m < 9) continue;

                // Check that the bonded non-common particles are bonded
                if (Bonds_BondCheck(hc11C[n11C][5], hc11C[n11C][7]) == 0) continue;
                if (Bonds_BondCheck(hc11C[n11C][6], hc11C[n11C][8]) == 0) continue;

                // Get the ID's of the non-common particles
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_first_7A][k], hc11C[n11C][5]) &&
                        Bonds_BondCheck(hcsp5c[id_first_7A][k], hc11C[n11C][6])) {
                        hc11C[n11C][9] = hcsp5c[id_first_7A][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[id_second7A][k], hc11C[n11C][7]) &&
                        Bonds_BondCheck(hcsp5c[id_second7A][k], hc11C[n11C][8])) {
                        hc11C[n11C][10] = hcsp5c[id_second7A][k];
                    }
                }
                quickSort(&hc11C[n11C][1], 2);
                quickSort(&hc11C[n11C][3], 2);
                quickSort(&hc11C[n11C][5], 4);
                quickSort(&hc11C[n11C][9], 2);

                Cluster_Write_11C();

                ++n11C;
            }
        }
    }
}

int get_11C_spindle_particles(int *uncommon_spindle, int id_first_7A, int id_second7A, int *common_spindle) {
    int num_common_spindles = 0;

    if (hcsp5c[id_first_7A][5] == hcsp5c[id_second7A][5]) {
        (*common_spindle) = hcsp5c[id_first_7A][5];
        uncommon_spindle[0] = hcsp5c[id_first_7A][6];
        uncommon_spindle[1] = hcsp5c[id_second7A][6];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][6] == hcsp5c[id_second7A][6]) {
        (*common_spindle) = hcsp5c[id_first_7A][6];
        uncommon_spindle[0] = hcsp5c[id_first_7A][5];
        uncommon_spindle[1] = hcsp5c[id_second7A][5];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][5] == hcsp5c[id_second7A][6]) {
        (*common_spindle) = hcsp5c[id_first_7A][5];
        uncommon_spindle[0] = hcsp5c[id_first_7A][6];
        uncommon_spindle[1] = hcsp5c[id_second7A][5];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][6] == hcsp5c[id_second7A][5]) {
        (*common_spindle) = hcsp5c[id_first_7A][6];
        uncommon_spindle[0] = hcsp5c[id_first_7A][5];
        uncommon_spindle[1] = hcsp5c[id_second7A][6];
        ++num_common_spindles;
    }

    if (num_common_spindles == 1) return 1;
    else return 0;
}

int get_bonded_7A_ring_particles(int *ar, int id_first_7A, int id_second7A, int ncom) {
    int first_ring_pointer, second_ring_pointer;

    for (first_ring_pointer = 0; first_ring_pointer < 5; ++first_ring_pointer) {
        for (second_ring_pointer = 0; second_ring_pointer < 5; ++second_ring_pointer) {
            if (hcsp5c[id_first_7A][first_ring_pointer] == hcsp5c[id_second7A][second_ring_pointer]) {
                if (ncom == 2) {
                    ++ncom;
                    break;
                }
                ar[ncom++] = hcsp5c[id_first_7A][first_ring_pointer];
                break;
            }
        }
        if (ncom > 2) break;
    }
    return ncom;
}

void resize_hc11C() {
    int clusSize=11;

    if(n11C == m11C) {
        hc11C=resize_2D_int(hc11C,m11C,m11C+incrStatic,clusSize,-1);
        m11C=m11C+incrStatic;
    }
}

void Cluster_Write_11C() {
    int i;

    s11C[hc11C[n11C][0]] = 'S';
    if(s11C[hc11C[n11C][1]] != 'S') s11C[hc11C[n11C][1]] = 'O';
    if(s11C[hc11C[n11C][2]] != 'S') s11C[hc11C[n11C][2]] = 'O';
    for(i=3; i< 11; i++) {
        if (s11C[hc11C[n11C][i]] == 'C') s11C[hc11C[n11C][i]] = 'B';
    }
}

int Clusters_Get12A() {
    // A 12A is an 11C with an extra particle bonded to only 2 other specific outer shell particles in the 11C.

    int id_11C;
    int ep;

    for(id_11C=0; id_11C<n11C; id_11C++) {
        if (cnb[hc11C[id_11C][0]] == 11) {

            ep = get_12A_extra_particle(id_11C);

            // Extra particle must be bonded to a specific two particles in the 11C
            if (Bonds_BondCheck(ep, hc11C[id_11C][9]) == 0 || Bonds_BondCheck(ep, hc11C[id_11C][10]) == 0) continue;

            // The extra particle should not be bonded to particles 2-8 of the 11C
            if (bond_check_12A_extra_particle(id_11C, ep) == 1) continue;

            resize_hc12A();
            populate_hc12A(id_11C, ep);
            populate_s12A();
            ++n12A;
        }
    }
}

int get_12A_extra_particle(int id_11C) {
    int i;
    // Returns id of extra particle
    // The extra particle is the one bonded to the 11C center that is not in the 11C,
    for (i = 0; i < 11; ++i) {
        if (is_particle_in_11C(bNums[hc11C[id_11C][0]][i], id_11C) == 0) {
            return bNums[hc11C[id_11C][0]][i]; // The extra particle
        }
    }
}

int bond_check_12A_extra_particle(int id_11C, int extra_particle) {
    // Return 1 if particle is bonded to particles 1-8 of the 11C, else return 0
    int i;

    for (i = 1; i < 9; ++i) {
        if (Bonds_BondCheck(extra_particle, hc11C[id_11C][i])){
            return 1;
        }
    }
    return 0;
}

int is_particle_in_11C(int particle_id, int id_11C) {
    // Return 1 if particle is in 11C, else returns 0
    int i;

    for (i=1; i<11; i++) {
        if (particle_id == hc11C[id_11C][i]) {
            return 1;
        }
    }
    return 0;
}

void populate_hc12A(int id_11C, int ep) {
    int i;

    for (i = 0; i<11; i++) {
        hc12A[n12A][i] = hc11C[id_11C][i];
    }
    hc12A[n12A][11] = ep;
}

void resize_hc12A() {
    int clusSize=12;
    if (n12A == m12A) {
        hc12A = resize_2D_int(hc12A, m12A, m12A + incrStatic, clusSize, -1);
        m12A = m12A + incrStatic;
    }
}

void populate_s12A() {
    int i;
    // hc12A key: (as 11C, extra_s)
    s12A[hc12A[n12A][0]] = 'S';
    if(s12A[hc12A[n12A][1]] != 'S') s12A[hc12A[n12A][1]] = 'O';
    if(s12A[hc12A[n12A][2]] != 'S') s12A[hc12A[n12A][2]] = 'O';
    for(i=3; i<12; i++) {
        if (s12A[hc12A[n12A][i]] == 'C') s12A[hc12A[n12A][i]] = 'B';
    }
}

void Clusters_Get11F_12E_13K() {   // Detect 11F C2v & 12E 3h
    int cp, bpi, bpj, ep1, ep2, the6A_i, the6A_j;
    int i, j, j2, k, l, m;
    int flg, flg1, flg2;
    int break_out;
    int clusSize=11;

    cp=bpi=bpj=ep1=ep2=the6A_i=the6A_j=-1;

    for(i=0; i<nsp3c-1; i++) {   // loop over all sp3c clusters
        for (j2=0; j2<3; j2++) {    // loop over only the rings of the sp3c clusters
            for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]; ++j) {
                if (mem_sp3c[hcsp3c[i][j2]][j]>i) {
                    flg = hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] ||
                          hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3] ||
                          hcsp3c[i][3] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4] ||
                          hcsp3c[i][4] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                    if (flg == 1) continue;
                    if (Bonds_BondCheck(hcsp3c[i][3], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) == 1 &&
                        Bonds_BondCheck(hcsp3c[i][4], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) == 1) {
                        m = 0;
                        for (k = 0; k < 3; ++k) {
                            for (l = 0; l < 3; ++l) {
                                if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                                    cp = hcsp3c[i][k];
                                    ++m;
                                }
                            }
                        }

                        if (m == 1) { // only 1 particle common between the pair
                            m = 0;
                            for (k = 0; k < 3; ++k) {
                                for (l = 0; l < 3; ++l) {
                                    if (hcsp3c[i][k] == cp || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == cp) continue;
                                    if (Bonds_BondCheck(hcsp3c[i][k], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) == 1) {
                                        if (m++) break;
                                        bpi = hcsp3c[i][k];
                                        bpj = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                                    }
                                }
                            }

                            if (m == 1) { // There must be exactly one paticle of sp3_i bonded to exactly one of sp3_j
                                flg1 = flg2 = 0;
                                for (k = 0; k < nsp4c; ++k) {
                                    if (hcsp4c[k][4] == cp || hcsp4c[k][5] == cp) {
                                        if (flg1 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = hcsp4c[k][l] == hcsp3c[i][3] ||
                                                      hcsp4c[k][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                                                flg = flg || hcsp4c[k][l] == bpi || hcsp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg1 = 1;
                                                if (hcsp4c[k][4] == cp) ep1 = hcsp4c[k][5];
                                                else ep1 = hcsp4c[k][4];
                                                the6A_i = k;
                                            }
                                        }
                                        if (flg2 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = hcsp4c[k][l] == hcsp3c[i][4] ||
                                                      hcsp4c[k][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                                                flg = flg || hcsp4c[k][l] == bpi || hcsp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg2 = 1;
                                                if (hcsp4c[k][4] == cp) ep2 = hcsp4c[k][5];
                                                else ep2 = hcsp4c[k][4];
                                                the6A_j = k;
                                            }
                                        }
                                    }
                                    if (flg1 == 1 && flg2 == 1) break;
                                }

                                if (k < nsp4c) { // 11F found
                                    if (n11F == m11F) {
                                        hc11F = resize_2D_int(hc11F, m11F, m11F + incrStatic, clusSize, -1);
                                        m11F = m11F + incrStatic;
                                    }

                                    // hc11F key: (com, ep1(6A extra spindle), ep2(6A extra spindle), 5A_i_s1, 5A_i_s2, 5A_j_s1, 5A_j_s2, 4 ring particles)

                                    hc11F[n11F][0] = cp;
                                    hc11F[n11F][1] = ep1;
                                    hc11F[n11F][2] = ep2;
                                    hc11F[n11F][3] = hcsp3c[i][3];
                                    hc11F[n11F][4] = hcsp3c[i][4];
                                    hc11F[n11F][5] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                                    hc11F[n11F][6] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];

                                    m = 7;
                                    break_out = 0;
                                    for (l = 0; l < 3; ++l) {
                                        if (hcsp3c[i][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F][m] = hcsp3c[i][l];
                                            m++;
                                        }
                                    }
                                    for (l = 0; l < 3; ++l) {
                                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F][m] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                                            m++;
                                        }
                                    }
                                    if (break_out == 1 || m < 11) continue;

                                    quickSort(&hc11F[n11F][1], 2);
                                    quickSort(&hc11F[n11F][3], 4);
                                    quickSort(&hc11F[n11F][7], 4);

                                    Cluster_Write_11F();

                                    if (do12E == 1) {
                                        if (Clusters_Get12E_D3h(mem_sp3c[hcsp3c[i][j2]][j])) ++n12E;
                                    }

                                    if (do13K == 1) {
                                        if (Clusters_Get13K(i, mem_sp3c[hcsp3c[i][j2]][j], the6A_i))
                                            ++n13K;
                                    }
                                    ++n11F;
                                }
                            }
                        }
                    }

                    if (Bonds_BondCheck(hcsp3c[i][3], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) == 1 &&
                        Bonds_BondCheck(hcsp3c[i][4], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) == 1) {
                        m = 0;
                        for (k = 0; k < 3; ++k) {
                            for (l = 0; l < 3; ++l) {
                                if (hcsp3c[i][k] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) {
                                    cp = hcsp3c[i][k];
                                    ++m;
                                }
                            }
                        }

                        if (m == 1) { // only 1 particle common between the pair
                            m = 0;
                            for (k = 0; k < 3; ++k) {
                                for (l = 0; l < 3; ++l) {
                                    if (hcsp3c[i][k] == cp || hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] == cp) continue;
                                    if (Bonds_BondCheck(hcsp3c[i][k], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l]) == 1) {
                                        if (m++) break;
                                        bpi = hcsp3c[i][k];
                                        bpj = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                                    }
                                }
                            }

                            if (m == 1) { // The single bonded particle found
                                flg1 = flg2 = 0;
                                for (k = 0; k < nsp4c; ++k) {
                                    if (hcsp4c[k][4] == cp || hcsp4c[k][5] == cp) {
                                        if (flg1 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = hcsp4c[k][l] == hcsp3c[i][3] ||
                                                      hcsp4c[k][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];
                                                flg = flg || hcsp4c[k][l] == bpi || hcsp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg1 = 1;
                                                if (hcsp4c[k][4] == cp) ep1 = hcsp4c[k][5];
                                                else ep1 = hcsp4c[k][4];
                                                the6A_i = k;
                                            }
                                        }
                                        if (flg2 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = hcsp4c[k][l] == hcsp3c[i][4] ||
                                                      hcsp4c[k][l] == hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                                                flg = flg || hcsp4c[k][l] == bpi || hcsp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg2 = 1;
                                                if (hcsp4c[k][4] == cp) ep2 = hcsp4c[k][5];
                                                else ep2 = hcsp4c[k][4];
                                                the6A_j = k;
                                            }
                                        }
                                    }
                                    if (flg1 == 1 && flg2 == 1) break;
                                }
                                if (k < nsp4c) { // 11F found
                                    if (n11F == m11F) {
                                        hc11F = resize_2D_int(hc11F, m11F, m11F + incrStatic, clusSize, -1);
                                        m11F = m11F + incrStatic;
                                    }

                                    // hc11F key: (com, ep1(6A extra spindle), ep2(6A extra spindle), 5A_i_s1, 5A_i_s2, 5A_j_s1, 5A_j_s2, 4 ring particles)

                                    hc11F[n11F][0] = cp;
                                    hc11F[n11F][1] = ep1;
                                    hc11F[n11F][2] = ep2;
                                    hc11F[n11F][3] = hcsp3c[i][3];
                                    hc11F[n11F][4] = hcsp3c[i][4];
                                    hc11F[n11F][5] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3];
                                    hc11F[n11F][6] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4];

                                    m = 7;
                                    break_out = 0;
                                    for (l = 0; l < 3; ++l) {
                                        if (hcsp3c[i][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F][m] = hcsp3c[i][l];
                                            m++;
                                        }
                                    }
                                    for (l = 0; l < 3; ++l) {
                                        if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F][m] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                                            m++;
                                        }
                                    }
                                    if (break_out == 1 || m < 11) continue;

                                    quickSort(&hc11F[n11F][1], 2);
                                    quickSort(&hc11F[n11F][3], 4);
                                    quickSort(&hc11F[n11F][7], 4);

                                    Cluster_Write_11F();

                                    if (do12K == 1) {
                                        if (Clusters_Get12E_D3h(mem_sp3c[hcsp3c[i][j2]][j])) ++n12E;
                                    }

                                    if (do13K == 1) {
                                        if (Clusters_Get13K(i, mem_sp3c[hcsp3c[i][j2]][j], the6A_i))
                                            ++n13K;
                                    }
                                    ++n11F;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Cluster_Write_11F() {
    int i;
    if (s11F[hc11F[n11F][0]] == 'C') s11F[hc11F[n11F][0]] = 'B';
    for(i=1; i<7; i++) {
        s11F[hc11F[n11F][i]] = 'O';
    }
    for(i=7; i<11; i++) {
        if (s11F[hc11F[n11F][i]] == 'C') s11F[hc11F[n11F][i]] = 'B';
    }
}

int Clusters_Get12E_D3h(int j) {  // Return 1 is 11F is also 12E
    //  Made from three sp3c or 5A clusters
    int k, l, m, ncom, common[2], uncom;
    int flg;
    int clusSize=12;

    for(k=j+1; k<nsp3c; ++k) {
        flg = (hcsp3c[k][3] == hc11F[n11F][1] && hcsp3c[k][4] == hc11F[n11F][2])  || (hcsp3c[k][3] == hc11F[n11F][2] && hcsp3c[k][4] == hc11F[n11F][1]);
        if (flg==1) { // spindles bonded correctly :)
            ncom=common[0]=common[1]=0;
            for (l=0; l<3; ++l) {
                for (m=0; m<5; ++m) {
                    if (m==0) uncom=0;
                    else uncom=m+6;
                    if(hcsp3c[k][l] == hc11F[n11F][uncom]) {
                        if (ncom==2) {
                            ncom++;
                            break;
                        }
                        common[ncom]=hcsp3c[k][l];
                        ncom++;
                    }
                }
                if (ncom>2) break;
            }
            if (ncom!=2) continue;

            uncom=-1;
            for (l=0; l<3; ++l) {
                ncom=0;
                for (m=0; m<2; ++m) {
                    if(hcsp3c[k][l] == common[m]) ncom++;
                }
                if (ncom==0) {
                    uncom=hcsp3c[k][l];
                    break;
                }
            }

            // now we have found the 12E
            if(n12E == m12E) {
                hc12E=resize_2D_int(hc12E,m12E,m12E+incrStatic,clusSize,-1);
                m12E=m12E+incrStatic;
            }

            hc12E[n12E][11] = uncom;
            for(l=0; l<10; ++l) hc12E[n12E][l] = hc11F[n11F][l+1];
            hc12E[n12E][10] = hc11F[n11F][0];
            quickSort(&hc12E[n12E][0],6);
            quickSort(&hc12E[n12E][6],6);

            Cluster_Write_12E();
            return 1;
        }
    }
    return 0;
}

void Cluster_Write_12E() {
    int i;
    for(i=0; i<6; i++){
        s12E[hc12E[n12E][i]] = 'O';
    }
    for(i=6; i<12; i++){
        if(s12E[hc12E[n12E][i]] == 'C') s12E[hc12E[n12E][i]] = 'B';
    }
}

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
        if (hcsp3c[sp3c_j][i]!=hc11F[n11F][0]) {
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
    for (i=0; i<nmem_sp3c[hc11F[n11F][0]]; ++i) {
        if (mem_sp3c[hc11F[n11F][0]][i]!=sp3c_i && mem_sp3c[hc11F[n11F][0]][i]!=sp3c_j) {
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
    for (i=0; i<nmem_sp3c[hc11F[n11F][0]]; ++i) { // loop over all sp3c which hc13K[n13K][0] is a member
        if (mem_sp3c[hc11F[n11F][0]][i]!=sp3c_i && mem_sp3c[hc11F[n11F][0]][i]!=sp3c_j) {
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
        hc13K=resize_2D_int(hc13K,m13K,m13K+incrStatic,clusSize,-1);
        m13K=m13K+incrStatic;
    }
    // hc13K key: (11F, extra SP3 ring particle to make 5A #1, extra SP3 ring particle to make 5A #2)

    for (i=0; i<11; i++) hc13K[n13K][i] = hc11F[n11F][i];
    hc13K[n13K][11]=ep[0];
    hc13K[n13K][12]=ep[1];

    quickSort(&hc13K[n13K][11],2);
    quickSort(&eclus5A[0],2);

    Cluster_Write_13K();

    return 1;
}

void Cluster_Write_13K() {
    int i;
    for(i=1; i<11; i++) {
        if (s13K[hc13K[n13K][i]] == 'C') s13K[hc13K[n13K][i]] = 'B';
    }
    s13K[hc13K[n13K][0]] = 'S';
    if(s13K[hc13K[n13K][11]] != 'S') s13K[hc13K[n13K][11]] = 'O';
    if(s13K[hc13K[n13K][12]] != 'S') s13K[hc13K[n13K][12]] = 'O';
}

void Clusters_Get12B_13A() { // Detect 12B & 13A D5h clusters together
    int i, j, k, l, m;
    int sp1, sp2;
    int sj1[5], sj2[5];
    int nSB1, nSB2;
    int flg;
    int break_out;
    int clusSize=12;

    for (i=0; i<nsp5c; ++i) { //first 7A
        sp1 = hcsp5c[i][5];
        sp2 = hcsp5c[i][6];
        nSB1 = nSB2 = 0; // count up spindle bonds

        for (j=i+1; j<nsp5c; ++j) { // second 7A
            flg = sp1 == hcsp5c[j][5] && Bonds_BondCheck(sp2,hcsp5c[j][6]);
            flg = flg || (sp1 == hcsp5c[j][6] && Bonds_BondCheck(sp2,hcsp5c[j][5]));
            if (flg==1) {
                if (nSB1>=5) {
                    nSB1++;
                    break;
                }
                sj1[nSB1++] = j;
            }
        }

        if (nSB1 == 5 && do13A==1) {     // possibly found 13A, definately found 12B, now establish status
            for (j=i+1; j<nsp5c; ++j) {
                if (sp1 == hcsp5c[j][5] || sp1 == hcsp5c[j][6]) {
                    for (k=0; k<5; ++k) {
                        for (l=0; l<5; ++l) {
                            if (hcsp5c[i][k] == hcsp5c[j][l]) break;
                        }
                        if (l<5) break;
                    }
                    if(k==5) { // got 13A, Check all sp5c[j][] - sp1 spindles are less than i
                        for (k=0; k<i; ++k) {
                            for(l=0; l<5; ++l) {
                                if (hcsp5c[j][l] == hcsp5c[k][5] && sp1 == hcsp5c[k][6]) break;
                                if (hcsp5c[j][l] == hcsp5c[k][6] && sp1 == hcsp5c[k][5]) break;
                            }
                            if(l<5) break; // index k < i present
                        }
                        if(k==i) break; // no index k < i present
                    }
                }
            }
            if (j<nsp5c) { // 13A found
                if (n13A == m13A) {
                    hc13A=resize_2D_int(hc13A,m13A,m13A+incrStatic,clusSize+1,-1);
                    m13A=m13A+incrStatic;
                }

                hc13A[n13A][0] = sp1;
                k = 1;
                if(hcsp5c[i][5] != sp1) hc13A[n13A][k++] = hcsp5c[i][5];
                if(hcsp5c[i][6] != sp1) hc13A[n13A][k++] = hcsp5c[i][6];
                if(hcsp5c[j][5] != sp1) hc13A[n13A][k++] = hcsp5c[j][5];
                if(hcsp5c[j][6] != sp1) hc13A[n13A][k] = hcsp5c[j][6];
                hc13A[n13A][3] = hcsp5c[i][0];
                hc13A[n13A][4] = hcsp5c[i][1];
                hc13A[n13A][5] = hcsp5c[i][2];
                hc13A[n13A][6] = hcsp5c[i][3];
                hc13A[n13A][7] = hcsp5c[i][4];
                hc13A[n13A][8] = hcsp5c[j][0];
                hc13A[n13A][9] = hcsp5c[j][1];
                hc13A[n13A][10] = hcsp5c[j][2];
                hc13A[n13A][11] = hcsp5c[j][3];
                hc13A[n13A][12] = hcsp5c[j][4];
                quickSort(&hc13A[n13A][1],12);

                Clust_Write_13A();
            }
        }

        for (j=i+1; j<nsp5c; ++j) {
            flg = sp2 == hcsp5c[j][5] && Bonds_BondCheck(sp1,hcsp5c[j][6]);
            flg = flg || (sp2 == hcsp5c[j][6] && Bonds_BondCheck(sp1,hcsp5c[j][5]));
            if (flg==1) {
                if (nSB2>=5) {
                    nSB2++;
                    break;
                }
                sj2[nSB2++] = j;
            }
        }

        if(nSB2 == 5 && do13A==1) { // possibly found 13A, definately found 12B, now establish status
            for (j=i+1; j<nsp5c; ++j) {
                if (sp2 == hcsp5c[j][5] || sp2 == hcsp5c[j][6]) {
                    for (k=0; k<5; ++k) {
                        for (l=0; l<5; ++l) {
                            if (hcsp5c[i][k] == hcsp5c[j][l]) break;
                        }
                        if(l<5) break;
                    }
                    if (k==5) { // Check all sp5c[j][] - sp2 spindles are less than i
                        for (k=0; k<i; ++k) {
                            for (l=0; l<5; ++l) {
                                if (hcsp5c[j][l] == hcsp5c[k][5] && sp2 == hcsp5c[k][6]) break;
                                if (hcsp5c[j][l] == hcsp5c[k][6] && sp2 == hcsp5c[k][5]) break;
                            }
                            if (l<5) break;
                        }
                        if (k==i) break;
                    }
                }
            }
            if (j<nsp5c){ // 13A found, is it unique
                if (n13A == m13A) {
                    hc13A=resize_2D_int(hc13A,m13A,m13A+incrStatic,clusSize+1,-1);
                    m13A=m13A+incrStatic;
                }

                hc13A[n13A][0] = sp2;
                k = 1;
                if(hcsp5c[i][5] != sp2) hc13A[n13A][k++] = hcsp5c[i][5];
                if(hcsp5c[i][6] != sp2) hc13A[n13A][k++] = hcsp5c[i][6];
                if(hcsp5c[j][5] != sp2) hc13A[n13A][k++] = hcsp5c[j][5];
                if(hcsp5c[j][6] != sp2) hc13A[n13A][k] = hcsp5c[j][6];
                hc13A[n13A][3] = hcsp5c[i][0];
                hc13A[n13A][4] = hcsp5c[i][1];
                hc13A[n13A][5] = hcsp5c[i][2];
                hc13A[n13A][6] = hcsp5c[i][3];
                hc13A[n13A][7] = hcsp5c[i][4];
                hc13A[n13A][8] = hcsp5c[j][0];
                hc13A[n13A][9] = hcsp5c[j][1];
                hc13A[n13A][10] = hcsp5c[j][2];
                hc13A[n13A][11] = hcsp5c[j][3];
                hc13A[n13A][12] = hcsp5c[j][4];
                quickSort(&hc13A[n13A][1],12);

                Clust_Write_13A();
            }
        }

        if ((nSB1 > 5) && (nSB2 > 5)) continue;

        for (j=0; j<i; ++j) { // keep looking for 12B
            flg = sp1 == hcsp5c[j][5] && Bonds_BondCheck(sp2,hcsp5c[j][6]);
            flg = flg || (sp1 == hcsp5c[j][6] && Bonds_BondCheck(sp2,hcsp5c[j][5]));
            if (flg==1) {
                if (nSB1 >= 5) {
                    nSB1++;
                    break;
                }
                sj1[nSB1++] = j;
            }
        }

        if (nSB1 == 5) {
            if (n12B==m12B) {
                hc12B=resize_2D_int(hc12B,m12B,m12B+incrStatic,clusSize,-1);
                m12B=m12B+incrStatic;
            }
            hc12B[n12B][0] = sp1;
            hc12B[n12B][1] = sp2;
            hc12B[n12B][2] = hcsp5c[i][0];
            hc12B[n12B][3] = hcsp5c[i][1];
            hc12B[n12B][4] = hcsp5c[i][2];
            hc12B[n12B][5] = hcsp5c[i][3];
            hc12B[n12B][6] = hcsp5c[i][4];

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
            quickSort(&hc12B[n12B][7],5);

            Clust_Write_12B();
        }

        for (j=0; j<i; ++j) {
            flg = sp2 == hcsp5c[j][5] && Bonds_BondCheck(sp1,hcsp5c[j][6]);
            flg = flg || (sp2 == hcsp5c[j][6] && Bonds_BondCheck(sp1,hcsp5c[j][5]));
            if (flg==1) {
                if (nSB2 >= 5) {
                    nSB2++;
                    break;
                }
                sj2[nSB2++] = j;
            }
        }

        if(nSB2 == 5) {
            if (n12B==m12B) {
                hc12B=resize_2D_int(hc12B,m12B,m12B+incrStatic,clusSize,-1);
                m12B=m12B+incrStatic;
            }

            hc12B[n12B][0] = sp2;
            hc12B[n12B][1] = sp1;
            hc12B[n12B][2] = hcsp5c[i][0];
            hc12B[n12B][3] = hcsp5c[i][1];
            hc12B[n12B][4] = hcsp5c[i][2];
            hc12B[n12B][5] = hcsp5c[i][3];
            hc12B[n12B][6] = hcsp5c[i][4];

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
            quickSort(&hc12B[n12B][7],5);

            Clust_Write_12B();

        }
    }
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

void Clust_Write_13A() {
    int i;
    s13A[hc13A[n13A][0]] = 'S';
    if(s13A[hc13A[n13A][1]] != 'S') s13A[hc13A[n13A][1]] = 'O';
    if(s13A[hc13A[n13A][2]] != 'S') s13A[hc13A[n13A][2]] = 'O';
    for(i=3; i<13; i++){
        if (s13A[hc13A[n13A][i]] == 'C') s13A[hc13A[n13A][i]] = 'B';
    }

    ++n13A;
}

void Clusters_Get13B_D5h() {   // Detect 13B D5h clusters, i.e. twisted icosahedra
    int cp;
    int i, j, k, l, m;
    int flg;
    int clusSize=13;

    cp=-1;

    for(i=0; i<nsp5c-1; ++i){
        // POSSIBLE IMPROVEMENT: use 7A clusters of spindles of 7A_i
        for(j=i+1; j<nsp5c; ++j){
            flg = 0;
            if (hcsp5c[i][5]==hcsp5c[j][5] && hcsp5c[i][6]!=hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][6], hcsp5c[j][6])==1) continue;
                cp = hcsp5c[i][5];
                flg = 1;
            }
            if (hcsp5c[i][6]==hcsp5c[j][6] && hcsp5c[i][5]!=hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][5],hcsp5c[j][5])==1) continue;
                cp = hcsp5c[i][6];
                flg = 1;
            }
            if (hcsp5c[i][5]==hcsp5c[j][6] && hcsp5c[i][6]!=hcsp5c[j][5]) {
                if (Bonds_BondCheck(hcsp5c[i][6],hcsp5c[j][5])==1) continue;
                cp = hcsp5c[i][5];
                flg = 1;
            }
            if (hcsp5c[i][6]==hcsp5c[j][5] && hcsp5c[i][5]!=hcsp5c[j][6]) {
                if (Bonds_BondCheck(hcsp5c[i][5],hcsp5c[j][6])==1) continue;
                cp = hcsp5c[i][6];
                flg = 1;
            }
            if (flg==1) {
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if(hcsp5c[i][k] == hcsp5c[j][l]) break;
                    }
                    if(l<5) break;
                }
                if (k<5) continue;
                for (k=0; k<5; ++k) {
                    m = 0;
                    for (l=0; l<5; ++l) {
                        if (Bonds_BondCheck(hcsp5c[i][k],hcsp5c[j][l])==1) ++m;
                        if (m == 2) break;
                    }
                    if(m != 1) break;
                }
                if (k<5) continue;
                // ERROR: need to check converse, i.e. every SP5 from 7A_j bonded to one from SP5 from 7A_i

                if (n13B == m13B) {
                    hc13B=resize_2D_int(hc13B,m13B,m13B+incrStatic,clusSize,-1);
                    m13B=m13B+incrStatic;
                }

                hc13B[n13B][0] = cp;
                k = 1;
                if(hcsp5c[i][5] != cp) hc13B[n13B][k++] = hcsp5c[i][5];
                if(hcsp5c[i][6] != cp) hc13B[n13B][k++] = hcsp5c[i][6];
                if(hcsp5c[j][5] != cp) hc13B[n13B][k++] = hcsp5c[j][5];
                if(hcsp5c[j][6] != cp) hc13B[n13B][k++] = hcsp5c[j][6];
                for(l=0; l<5; ++l) hc13B[n13B][k++] = hcsp5c[i][l];
                for(l=0; l<5; ++l) hc13B[n13B][k++] = hcsp5c[j][l];
                quickSort(&hc13B[n13B][1],2);
                quickSort(&hc13B[n13B][3],10);

                Cluster_Write_13B();
            }
        }
    }
}

void Cluster_Write_13B() {
    int i;
    s13B[hc13B[n13B][0]] = 'S';
    if(s13B[hc13B[n13B][1]] != 'S') s13B[hc13B[n13B][1]] = 'O';
    if(s13B[hc13B[n13B][2]] != 'S') s13B[hc13B[n13B][2]] = 'O';
    for(i=3; i<13; i++) {
        if (s13B[hc13B[n13B][i]] == 'C') s13B[hc13B[n13B][i]] = 'B';
    }

    ++n13B;
}

void Clusters_GetFCC() {   // Detect 13 particle FCC clusters
    int i, j, j2, k, l, m, n;
    int i1, i2, i3;
    int cp, bpi, bpj, nbpi, nbpj;
    int flg1, flg2, flg3;
    int clusSize=13;

    cp=bpi=bpj=nbpi=nbpj=i3=-1;


    for (i=0; i<nsp3b-2; ++i) { // loop over all sp3b_i
        for (j2=0; j2<3; j2++) {
            for (j=0; j<nmem_sp3b[hcsp3b[i][j2]]-1; ++j) { // loop over all sp3b_j
                if (mem_sp3b[hcsp3b[i][j2]][j]<=i) continue;
                if (Bonds_BondCheck(hcsp3b[i][3], hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3])==0) continue; // spindle spindle bond
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

                for (k=0; k<nFCC; ++k) {
                    if (hcFCC[k][0] == cp) break;   // check for other degenerate FCC clusters which have cp here
                }
                if (k < nFCC) continue;  // found this fcc cluster before

                for (k=0; k<3; ++k) {
                    if (hcsp3b[i][k] == cp) continue;
                    for (l=0; l<3; ++l) {
                        if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l] == cp) continue;
                        if (Bonds_BondCheck(hcsp3b[i][k],hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l])==1) {
                            bpi = hcsp3b[i][k];
                            bpj = hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][l];
                            ++n;
                        }
                    }
                }
                if(n != 1) continue;    // one bond between one pair of atoms from SP3_i and SP3_j

                // we may have an FCC cluster, build it into hcFCC then overwrite it later if it aint
                if (nFCC == mFCC) {
                    hcFCC=resize_2D_int(hcFCC,mFCC,mFCC+incrStatic,clusSize,-1);
                    mFCC=mFCC+incrStatic;
                }

                for (k=0; k<3; ++k) {
                    if (hcsp3b[i][k]!=cp && hcsp3b[i][k]!=bpi) nbpi=hcsp3b[i][k];
                    if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][k]!=cp && hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][k]!=bpj) nbpj=hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][k];
                }
                hcFCC[nFCC][0] = cp;
                hcFCC[nFCC][1] = nbpi;
                hcFCC[nFCC][2] = bpi;
                hcFCC[nFCC][3] = hcsp3b[i][3];
                hcFCC[nFCC][4] = hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3];
                hcFCC[nFCC][5] = nbpj;
                hcFCC[nFCC][6] = bpj;   // store first two clusters

                for (k=j+1; k<nmem_sp3b[hcsp3b[i][j2]]; ++k) { // loop over all sp3b_k
                    if (mem_sp3b[hcsp3b[i][j2]][k]<=i) continue;
                    for (l=0; l<nFCC; ++l) {
                        if (hcFCC[l][0] == cp) break;   // check for other degenerate FCC clusters which have cp here
                    }
                    if (l < nFCC) break; // found this fcc cluster before

                    if (hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3]==cp) continue;    // check sp3b_k spindle isnt common particle
                    if (!(Bonds_BondCheck(hcsp3b[i][3],hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3])==1 && Bonds_BondCheck(hcsp3b[mem_sp3b[hcsp3b[i][j2]][j]][3], hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3])==1)) continue;  // check sp3b_k spindle is bonded to sp3b_i sp3b_j spindles
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
                        if (Bonds_BondCheck(n, nbpi)==1 && Bonds_BondCheck(n,nbpj)==0) {
                            i2 = n;
                            ++m;
                            flg1 =1;
                        }
                        if (Bonds_BondCheck(n, nbpj)==1 && Bonds_BondCheck(n, nbpi)==0) {
                            i3 = n;
                            ++m;
                            flg2 = 1;
                        }
                    }
                    if (m != 2) continue;  // SP3_k must have 2 particles bonded to the non-bonded non-common pair from SP3_i and SP3_j

                    if (!(flg1==1 && flg2==1)) continue; // now we have the 6 membered ring

                    for (l=0; l<nsp3b; ++l) { // loop over all sp3b_l
                        if (l==i || l==mem_sp3b[hcsp3b[i][j2]][j] || l==mem_sp3b[hcsp3b[i][j2]][k]) continue;   // must be new sp3b
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
                            if (Bonds_BondCheck(n,bpi)==1 && Bonds_BondCheck(n,bpj)==1) {
                                if (flg1==1) break;
                                flg1 = 1;
                                hcFCC[nFCC][10] = n;
                            }
                            if (Bonds_BondCheck(n,nbpj)==1 && Bonds_BondCheck(n,i3)==1) {
                                if(flg2==1) break;
                                flg2 = 1;
                                hcFCC[nFCC][11] = n;
                            }
                            if (Bonds_BondCheck(n,i2)==1 && Bonds_BondCheck(n,nbpi)==1) {
                                if(flg3==1) break;
                                flg3 = 1;
                                hcFCC[nFCC][12] = n;
                            }
                        }
                        if(m<3) continue;  // SP3_l particles must be bonded correctly to six-membered ring
                        if(flg1==1 && flg2==1 && flg3==1) break;   // SP3_l particles must be bonded correctly to six-membered ring
                    } // end l loop
                    if (l==nsp3b) { // required SP3b cluster not found
                        for (l=0; l<nsp3c; ++l) { // loop over all sp3c_l
                            if (hcsp3c[l][3]!=cp && hcsp3c[l][4]!=cp) continue;  // common particle must be a spindle
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
                                if (Bonds_BondCheck(n,bpi)==1 && Bonds_BondCheck(n,bpj)==1) {
                                    if (flg1==1) break;
                                    flg1 = 1;
                                    hcFCC[nFCC][10] = n;
                                }
                                if (Bonds_BondCheck(n,nbpj)==1 && Bonds_BondCheck(n,i3)==1) {
                                    if (flg2==1) break;
                                    flg2 = 1;
                                    hcFCC[nFCC][11] = n;
                                }
                                if (Bonds_BondCheck(n,i2)==1 && Bonds_BondCheck(n,nbpi)==1) {
                                    if(flg3==1) break;
                                    flg3 = 1;
                                    hcFCC[nFCC][12] = n;
                                }
                            }
                            if (m<3) continue;  // SP3_l particles must be bonded correctly to six-membered ring
                            if (flg1==1 && flg2==1 && flg3==1) break;   // SP3_l particles must be bonded correctly to six-membered ring
                        }
                        if (l==nsp3c) continue;
                    } // End if statement for SP3c search loop

                    // We've now found an FCC cluster
                    if (nFCC == mFCC) { hcFCC=resize_2D_int(hcFCC,mFCC,mFCC+incrStatic,clusSize,-1);
                        mFCC=mFCC+incrStatic;
                    }

                    hcFCC[nFCC][7] = hcsp3b[mem_sp3b[hcsp3b[i][j2]][k]][3];
                    hcFCC[nFCC][8] = i3;
                    hcFCC[nFCC][9] = i2;
                    quickSort(&hcFCC[nFCC][1],12);

                    Cluster_Write_FCC();
                }
            }
        }
    }
}

void Cluster_Write_FCC() {
    if(sFCC[hcFCC[nFCC][1]] == 'C') sFCC[hcFCC[nFCC][1]] = 'B';
    if(sFCC[hcFCC[nFCC][2]] == 'C') sFCC[hcFCC[nFCC][2]] = 'B';
    if(sFCC[hcFCC[nFCC][5]] == 'C') sFCC[hcFCC[nFCC][5]] = 'B';
    if(sFCC[hcFCC[nFCC][6]] == 'C') sFCC[hcFCC[nFCC][6]] = 'B';
    if(sFCC[hcFCC[nFCC][8]] == 'C') sFCC[hcFCC[nFCC][8]] = 'B';
    if(sFCC[hcFCC[nFCC][9]] == 'C') sFCC[hcFCC[nFCC][9]] = 'B';
    if(sFCC[hcFCC[nFCC][3]]  != 'S') sFCC[hcFCC[nFCC][3]] = 'O';
    if(sFCC[hcFCC[nFCC][4]]  != 'S') sFCC[hcFCC[nFCC][4]] = 'O';
    if(sFCC[hcFCC[nFCC][7]]  != 'S') sFCC[hcFCC[nFCC][7]] = 'O';
    if(sFCC[hcFCC[nFCC][10]] != 'S') sFCC[hcFCC[nFCC][10]] = 'O';
    if(sFCC[hcFCC[nFCC][11]] != 'S') sFCC[hcFCC[nFCC][11]] = 'O';
    if(sFCC[hcFCC[nFCC][12]] != 'S') sFCC[hcFCC[nFCC][12]] = 'O';
    sFCC[hcFCC[nFCC][0]] = 'S';

    ++nFCC;
}

void Clusters_GetHCP() {   // Detect 13 particle HCP clusters
    int i, j, j2, k, l, m, n;
    int ia[2], ja[2], ka[2];
    int cp, x;
    int h1i, h1j, h2i, h2j;
    int flg1, flg2, flg3;
    int clusSize=13;

    cp=-1;

    for (i=0; i<nsp3c-2; ++i) { // loop over all sp3c_i
        for (j2=0; j2<3; j2++) {
        for (j=0; j<nmem_sp3c[hcsp3c[i][j2]]-1; ++j) { // loop over all sp3c_j
            if (mem_sp3c[hcsp3c[i][j2]][j]<=i) continue;
            flg1 = Bonds_BondCheck(hcsp3c[i][3],hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) && Bonds_BondCheck(hcsp3c[i][4],hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]);
            flg2 = Bonds_BondCheck(hcsp3c[i][3],hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]) && Bonds_BondCheck(hcsp3c[i][4],hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]);
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
            if (cp!=hcsp3c[i][j2]) continue;

            for (k=j+1; k<nmem_sp3c[hcsp3c[i][j2]]; ++k) { // loop over all sp3c_k
                if (mem_sp3c[hcsp3c[i][j2]][k]<=i) continue;
                flg1 = Bonds_BondCheck(h1i,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) && Bonds_BondCheck(h1j,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]);
                flg2 = Bonds_BondCheck(h2i,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]) && Bonds_BondCheck(h2j,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                flg1 = flg1 && flg2;
                flg2 = Bonds_BondCheck(h1i,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]) && Bonds_BondCheck(h1j,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                flg3 = Bonds_BondCheck(h2i,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) && Bonds_BondCheck(h2j,hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]);
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
                        flg1 = Bonds_BondCheck(hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) || Bonds_BondCheck(hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]);
                        flg1 = flg1 || Bonds_BondCheck(hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) || Bonds_BondCheck(hcsp3c[i][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                        if (flg1==1) break;
                    }
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l] != cp) {
                        flg1 = Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[i][3]) || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[i][4]);
                        flg1 = flg1 || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][3]) || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][4]);
                        if (flg1==1) break;
                    }
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] != cp) {
                        flg1 = Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[i][3]) || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[i][4]);
                        flg1 = flg1 || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][3]) || Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l], hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][4]);
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
                    if (Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l], ia[0])==1) {
                        flg1 = 1;
                        ja[0] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                        if(Bonds_BondCheck(ja[0], ia[1])==1) break;
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
                        if (Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l],ia[0])==1) {
                            flg2 = 1;
                            ja[0] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                            if (Bonds_BondCheck(ja[0],ia[1])==1) break;
                        }
                        else ja[1] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][l];
                    }
                    if (l<3 || flg2==0) continue;
                }
                if (!(flg1==1 || flg2==1)) continue; // found four uncommon particles from 5A_i and 5A_j and found the two that are bonded between these
                if (ja[0]==-1 || ja[1]==-1) continue;
                if (Bonds_BondCheck(ja[1],ia[0])==1 || Bonds_BondCheck(ja[1],ia[1])==1) continue;
                flg1 = 0;
                ka[0]=ka[1]=-1;
                for (l=0; l<3; ++l) {
                    if (hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l] == cp) continue;
                    if (Bonds_BondCheck(hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l],ia[1])) {
                        flg1 = 1;
                        ka[0] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                        if (Bonds_BondCheck(ka[0],ia[0])==1 || Bonds_BondCheck(ka[0],ja[0])==1 || Bonds_BondCheck(ka[0],ja[1])==1) break;
                    }
                    else ka[1] = hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][l];
                }
                if (l<3 || ka[0]==-1 || ka[1]==-1) continue;
                if (Bonds_BondCheck(ka[1],ja[1])==0) continue;
                if (Bonds_BondCheck(ka[1],ia[1])==1 || Bonds_BondCheck(ka[1],ia[0])==1 || Bonds_BondCheck(ka[1],ja[0])==1) continue;

                // We've now found an HCP cluster
                if (nHCP == mHCP) {
                    hcHCP=resize_2D_int(hcHCP,mHCP,mHCP+incrStatic,clusSize,-1);
                    mHCP=mHCP+incrStatic;
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
                quickSort(&hcHCP[nHCP][1],6);
                quickSort(&hcHCP[nHCP][7],6);

                Cluster_Write_HCP(i, j, j2, k);
            }
        }
        }
    }
}

void Cluster_Write_HCP(int i, int j, int j2, int k) {
    int counter;

    for (counter=0; counter<3; counter++){
        if (sHCP[hcsp3c[i][counter]] == 'C') sHCP[hcsp3c[i][counter]] = 'B';
        if (sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][counter]] == 'C') sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][counter]] = 'B';
        if (sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][counter]] == 'C') sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][counter]] = 'B';
    }
    for (counter=3; counter<5; counter++) {
        sHCP[hcsp3c[i][counter]] = 'O';
        sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][j]][counter]] = 'F';
        sHCP[hcsp3c[mem_sp3c[hcsp3c[i][j2]][k]][counter]] = 'H';
    }

    ++nHCP;
}

void Clusters_GetBCC_9() {
    int i, j, j2, k, l, m;
    int flg;
    int s_com=-1;
    int trial[9];
    int clusSize=9;

    for (i=0; i<nsp4b-1; i++) {
        for (j2=4; j2<5; j2++) {
            for (j=0; j<nmem_sp4b[hcsp4b[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4b[hcsp4b[i][j2]][j]<=i) continue;
                if (hcsp4b[i][4]!=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][4]) continue;
                s_com=hcsp4b[i][4];

                flg=0;
                for (k=0; k<4; k++) {
                    for (l=0; l<4; l++) {
                        if (hcsp4b[i][k]==hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4b[i][k],hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][k],hcsp4b[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]=hcsp4b[i][k];
                    trial[k+5]=hcsp4b[mem_sp4b[hcsp4b[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k<nBCC_9; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l]!=hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9 == mBCC_9) {
                    hcBCC_9=resize_2D_int(hcBCC_9,mBCC_9,mBCC_9+incrStatic,clusSize,-1);
                    mBCC_9=mBCC_9+incrStatic;
                }

                for (k=0; k<9; ++k) hcBCC_9[nBCC_9][k]=trial[k];

                Cluster_Write_BCC9();
            }
        }
    }

    for (i=0; i<nsp4c-1; i++) {
        for (j2=4; j2<6; j2++) {
            for (j=0; j<nmem_sp4c[hcsp4c[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4c[hcsp4c[i][j2]][j]<=i) continue;
                m=0;
                for (k=4; k<6; k++) {
                    for (l=4; l<6; l++) {
                        if (hcsp4c[i][k]==hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]) {
                            s_com=hcsp4c[i][k];
                            m++;
                        }
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<6; k++) {
                    if (hcsp4c[i][k]==s_com) continue;
                    for (l=0; l<6; l++) {
                        if (hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]==s_com) continue;
                        if (hcsp4c[i][k]==hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[i][k],hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][k],hcsp4c[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]=hcsp4c[i][k];
                    trial[k+5]=hcsp4c[mem_sp4c[hcsp4c[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k<nBCC_9; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l]!=hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9 == mBCC_9) {
                    hcBCC_9=resize_2D_int(hcBCC_9,mBCC_9,mBCC_9+incrStatic,clusSize,-1);
                    mBCC_9=mBCC_9+incrStatic;
                }
                for (k=0; k<9; ++k) hcBCC_9[nBCC_9][k]=trial[k];

                Cluster_Write_BCC9();
            }
        }
    }

    for (i=0; i<nsp4b; i++) {
        for (j2=4; j2<5; j2++) {
            for (j=0; j<nmem_sp4c[hcsp4b[i][j2]]; ++j) { // loop over all sp3c_j
                m=0;
                for (k=4; k<6; k++) {
                    if (hcsp4b[i][4]==hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][k]) {
                        s_com=hcsp4b[i][4];
                        m++;
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    for (l=0; l<6; l++) {
                        if (hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][l]==s_com) continue;
                        if (hcsp4b[i][k]==hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][l]) {
                            flg=1;
                            break;
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4b[i][k],hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    m=0;
                    for (l=0; l<4; l++) {
                        if (Bonds_BondCheck(hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][k],hcsp4b[i][l])) {
                            m++;
                            if (m==2) {
                                flg=1;
                                break;
                            }
                        }
                    }
                    if (flg==1) break;
                }
                if (flg==1) continue;

                trial[0]=s_com;
                for (k=0; k<4; k++) {
                    trial[k+1]=hcsp4b[i][k];
                    trial[k+5]=hcsp4c[mem_sp4c[hcsp4b[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k<nBCC_9; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l]!=hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9 == mBCC_9) {
                    hcBCC_9=resize_2D_int(hcBCC_9,mBCC_9,mBCC_9+incrStatic,clusSize,-1);
                    mBCC_9=mBCC_9+incrStatic;
                }
                for (k=0; k<9; ++k) hcBCC_9[nBCC_9][k]=trial[k];

                Cluster_Write_BCC9();

            }
        }
    }
}

void Cluster_Write_BCC9() {
    int i;
    for (i = 1; i< 9; i++){
        if (sBCC_9[hcBCC_9[nBCC_9][i]] == 'C') sBCC_9[hcBCC_9[nBCC_9][i]] = 'B';
    }
    sBCC_9[hcBCC_9[nBCC_9][0]] = 'S';

    ++nBCC_9;
}

void Clusters_GetBCC_15() {    // Detect 15 particle BCC clusters
    int i,j,k,l,m;
    int no_sp4cs,noSP4s;
    int sj[5];
    int clusSize=15;

    for (i=0; i<nsp4c; i++) {
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15 == mBCC_15) {
            hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
            mBCC_15=mBCC_15+incrStatic;
        }
        for (j=0; j<nBCC_15; j++) if (hcsp4c[i][4]==hcBCC_15[j][0]) break;
        if (j==nBCC_15) {
            hcBCC_15[nBCC_15][0]=hcsp4c[i][4];
            hcBCC_15[nBCC_15][1]=hcsp4c[i][5];
            for (j=0; j<4; j++) hcBCC_15[nBCC_15][j+7]=hcsp4c[i][j];

            no_sp4cs=0;
            noSP4s=4;
            sj[0]=sj[1]=sj[2]=sj[3]=sj[4]=-1;
            for (j=i+1; j<nsp4c; j++) {
                if (hcsp4c[j][4]==hcBCC_15[nBCC_15][0] || hcsp4c[j][5]==hcBCC_15[nBCC_15][0]) {
                    m=0;
                    for (l=1;l<2+no_sp4cs;l++) {
                        if (hcsp4c[j][4]==hcBCC_15[nBCC_15][l] || hcsp4c[j][5]==hcBCC_15[nBCC_15][l]) m++;
                    }
                    if (m==0) {
                        for (l=7;l<7+noSP4s;l++) {
                            if (hcsp4c[j][4]==hcBCC_15[nBCC_15][l] || hcsp4c[j][5]==hcBCC_15[nBCC_15][l]) m++;
                        }
                        if (m==0) {
                            for (k=0; k<4; k++) {
                                for (l=0;l<2+no_sp4cs;l++) {
                                    if (hcsp4c[j][k]==hcBCC_15[nBCC_15][l]) m++;
                                }
                            }
                            if (m==0) {
                                if (hcsp4c[j][4]==hcBCC_15[nBCC_15][0]) hcBCC_15[nBCC_15][2+no_sp4cs]=hcsp4c[j][5];
                                else hcBCC_15[nBCC_15][2+no_sp4cs]=hcsp4c[j][4];
                                for (k=0; k<4; k++) {
                                    m=0;
                                    for (l=7;l<7+noSP4s;l++) {
                                        if (hcsp4c[j][k]==hcBCC_15[nBCC_15][l]) m++;
                                    }
                                    if (m==0) {
                                        if (noSP4s>=8) {
                                            noSP4s++;
                                            break;
                                        }
                                        hcBCC_15[nBCC_15][7+noSP4s]=hcsp4c[j][k];
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
            hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
            mBCC_15=mBCC_15+incrStatic;
        }
        for (j=0; j<nBCC_15; j++) if (hcsp4c[i][5]==hcBCC_15[j][0]) break;
        if (j<nBCC_15) continue;

        hcBCC_15[nBCC_15][0]=hcsp4c[i][5];
        hcBCC_15[nBCC_15][1]=hcsp4c[i][4];
        for (j=0; j<4; j++) hcBCC_15[nBCC_15][j+7]=hcsp4c[i][j];

        no_sp4cs=0;
        noSP4s=4;
        sj[0]=sj[1]=sj[2]=sj[3]=sj[4]=-1;
        for (j=i+1; j<nsp4c; j++) {
            if (hcsp4c[j][4]==hcBCC_15[nBCC_15][0] || hcsp4c[j][5]==hcBCC_15[nBCC_15][0]) {
                m=0;
                for (l=1;l<2+no_sp4cs;l++) {
                    if (hcsp4c[j][4]==hcBCC_15[nBCC_15][l] || hcsp4c[j][5]==hcBCC_15[nBCC_15][l]) m++;
                }
                if (m==0) {
                    for (l=7;l<7+noSP4s;l++) {
                        if (hcsp4c[j][4]==hcBCC_15[nBCC_15][l] || hcsp4c[j][5]==hcBCC_15[nBCC_15][l]) m++;
                    }
                    if (m==0) {
                        for (k=0; k<4; k++) {
                            for (l=0;l<2+no_sp4cs;l++) {
                                if (hcsp4c[j][k]==hcBCC_15[nBCC_15][l]) m++;
                            }
                        }
                        if (m==0) {
                            if (hcsp4c[j][4]==hcBCC_15[nBCC_15][0]) hcBCC_15[nBCC_15][2+no_sp4cs]=hcsp4c[j][5];
                            else hcBCC_15[nBCC_15][2+no_sp4cs]=hcsp4c[j][4];
                            for (k=0; k<4; k++) {
                                m=0;
                                for (l=7;l<7+noSP4s;l++) {
                                    if (hcsp4c[j][k]==hcBCC_15[nBCC_15][l]) m++;
                                }
                                if (m==0) {
                                    if (noSP4s>=8) {
                                        noSP4s++;
                                        break;
                                    }
                                    hcBCC_15[nBCC_15][7+noSP4s]=hcsp4c[j][k];
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