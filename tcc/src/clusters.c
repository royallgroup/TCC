#include "globals.h"
#include "clusters.h"
#include "bonds.h"

void Clusters_Get6Z_C2v(int f) {    // Detect 6Z clusters from 2 5A clusters
    int flg;
    int i, j, j2, k, l;
    int cnt;
    int s1a, s2a, s1b, s2b;
    char *ach, errMsg[1000];
    int clusSize=6;

    s1a=s2a=s1b=s2b=-1;
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get6Z_C2v(): ach[] malloc out of memory\n");  Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = 'C';

    for (i=0; i<nsp3c[f]-1; ++i) {  // loop over all 5A_i
        for (j2=0; j2<1; ++j2) {
            for (j=0; j<nmem_sp3c[sp3c[i][j2]]; ++j) {  // loop over all 5A_j
                if (mem_sp3c[sp3c[i][j2]][j]<=i) continue;
                if (sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3] || sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4]) continue;   // no spindles of 5A_i and 5A_j are common
                if (sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4]) continue;

                // for each cluster one spindle particle is in sp3 ring of other 5A and other spindle isn't
                cnt = 0;    // check one 5A_i spindle are in sp3 ring of 5A_j
                flg = sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][0] || sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][1] || sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][2];
                if (flg==1){
                    s1a = sp3c[i][3];
                    s2a = sp3c[i][4];
                    ++cnt;
                }
                flg = sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][0] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][1] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][2];
                if (flg==1){
                    s1a = sp3c[i][4];
                    s2a = sp3c[i][3];
                    ++cnt;
                }
                if (cnt != 1) continue;
                cnt = 0;    // check one 5A_j spindle are in sp3 ring of 5A_i
                flg = sp3c[mem_sp3c[sp3c[i][j2]][j]][3] == sp3c[i][0] || sp3c[mem_sp3c[sp3c[i][j2]][j]][3] == sp3c[i][1] || sp3c[mem_sp3c[sp3c[i][j2]][j]][3] == sp3c[i][2];
                if (flg==1){
                    s1b = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                    s2b = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                    ++cnt;
                }
                flg = sp3c[mem_sp3c[sp3c[i][j2]][j]][4] == sp3c[i][0] || sp3c[mem_sp3c[sp3c[i][j2]][j]][4] == sp3c[i][1] || sp3c[mem_sp3c[sp3c[i][j2]][j]][4] == sp3c[i][2];
                if (flg==1){
                    s1b = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                    s2b = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                    ++cnt;
                }
                if (cnt != 1) continue;
                if (!Bonds_BondCheck(s1a, s1b)) continue;   // check spindles of i and j in sp3 ring of j and i respectively are bonded

                cnt = 0;    // check 2 particles in the sp3 rings of 5A_i and 5A_j are common
                for (k=0; k<3; ++k){
                    for (l=0; l<3; ++l){
                        if(sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]){
                            ++cnt;
                            break;
                        }
                    }
                }
                if (cnt != 2) continue;

                if (n6Z[f]==m6Z) {
                    hc6Z=resize_2D_int(hc6Z,m6Z,m6Z+incrStatic,clusSize,-1);
                    m6Z=m6Z+incrStatic;
                }
                // Now we have found the 6Z cluster
                if (s1a<s1b) {
                    hc6Z[n6Z[f]][0]=s1a;    // insert cluster
                    hc6Z[n6Z[f]][1]=s1b;
                    hc6Z[n6Z[f]][2]=s2a;
                    hc6Z[n6Z[f]][3]=s2b;
                }
                else {
                    hc6Z[n6Z[f]][0]=s1b;    // insert cluster
                    hc6Z[n6Z[f]][1]=s1a;
                    hc6Z[n6Z[f]][2]=s2b;
                    hc6Z[n6Z[f]][3]=s2a;
                }
                cnt=4;
                for (k=0; k<3; ++k) {
                    flg=1;
                    for (l=0; l<4; ++l){
                        if (sp3c[i][k]==hc6Z[n6Z[f]][l]) {
                            flg=0;
                            break;
                        }
                    }
                    if (flg==1) { hc6Z[n6Z[f]][cnt]=sp3c[i][k]; cnt++; }
                }
                if (hc6Z[n6Z[f]][5]<hc6Z[n6Z[f]][4]) {
                    k=hc6Z[n6Z[f]][5];
                    hc6Z[n6Z[f]][5]=hc6Z[n6Z[f]][4];
                    hc6Z[n6Z[f]][4]=k;
                }
                Cluster_Write_6Z(f, ach);
            }
        }
    }
    for (i=0; i<N; ++i) s6Z[i]=ach[i];
    free(ach);
}

void Cluster_Write_6Z(int f, char *ach) {// hc6Z key: (5A_i_s_in_SP3_j, 5A_j_s_in_SP3_i, 5A_i_s_oth, 5A_j_s_oth, SP3_i_j_com_1, SP3_i_j_com_2)
    ach[hc6Z[n6Z[f]][0]] = 'O';
    ach[hc6Z[n6Z[f]][1]] = 'O';
    ach[hc6Z[n6Z[f]][2]] = 'O';
    ach[hc6Z[n6Z[f]][3]] = 'O';
    if (ach[hc6Z[n6Z[f]][4]] == 'C') ach[hc6Z[n6Z[f]][4]] = 'B';
    if (ach[hc6Z[n6Z[f]][5]] == 'C') ach[hc6Z[n6Z[f]][5]] = 'B';

    ++n6Z[f];
}

void Clusters_Get7K(int f) {    // Detect 7K clusters from 2 5A clusters
    int i, j, j2, k, l, m;
    int scom, sother[2], sp3_com[2], sp3c_i_other, sp3c_j_other;
    char *ach, errMsg[1000];
    int clusSize=7;

    scom=sother[0]=sother[1]=sp3_com[0]=sp3_com[1]=sp3c_i_other=sp3c_j_other=-1;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get7K(): ach[] malloc out of memory\n");  Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = 'C';

    for (i=0; i<nsp3c[f]-1; ++i) {  // loop over all 5A_i
        for (j2=3; j2<5; ++j2) {    // loop over both spindles of 5A_i
            for (j=0; j<nmem_sp3c[sp3c[i][j2]]; ++j) {  // loop over all 5A_j common with spindle of 5A_i
                if (mem_sp3c[sp3c[i][j2]][j]<=i) continue;  // don't find 7K twice

                m=0;
                for (k=3; k<5; k++) {
                    for (l=3; l<5; l++) {
                        if (sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                            if (m>=1) {
                                m++;
                                break;
                            }
                            scom=sp3c[i][k];
                            m++;
                        }
                    }
                    if (m>=2) break;
                }
                if (m!=1) continue; // exactly one common spindle between 5A_i and 5A_j

                if (sp3c[i][3] == scom) sother[0]=sp3c[i][4];
                else sother[0]=sp3c[i][3];
                if (sp3c[mem_sp3c[sp3c[i][j2]][j]][3] == scom) sother[1]=sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                else sother[1]=sp3c[mem_sp3c[sp3c[i][j2]][j]][3];

                m=0;
                for (k=0; k<5; k++) {
                    if (sother[0]==sp3c[mem_sp3c[sp3c[i][j2]][j]][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other spindle of 5A_i distinct from whole 5A_j

                m=0;
                for (k=0; k<5; k++) {
                    if (sother[1]==sp3c[i][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other spindle of 5A_j distinct from whole 5A_i

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<3; l++) {
                        if (sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                            if (m>=2) {
                                m++;
                                break;
                            }
                            sp3_com[m]=sp3c[i][k];
                            m++;
                        }
                    }
                    if (m>=3) break;
                }
                if (m!=2) continue; // exactly two common particles in SP3 rings of 5A_i and 5A_j

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<2; l++) {
                        if (sp3c[i][k] == sp3_com[l]) break;
                    }
                    if (l==2) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        sp3c_i_other=sp3c[i][k];
                        m++;
                    }
                }
                if (m!=1) continue; // found other uncommon particle from SP3 ring of 5A_i

                m=0;
                for (k=0; k<3; k++) {
                    for (l=0; l<2; l++) {
                        if (sp3c[mem_sp3c[sp3c[i][j2]][j]][k] == sp3_com[l]) break;
                    }
                    if (l==2) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        sp3c_j_other=sp3c[mem_sp3c[sp3c[i][j2]][j]][k];
                        m++;
                    }
                }
                if (m!=1) continue; // found other uncommon particle from SP3 ring of 5A_i

                m=0;
                for (k=0; k<5; k++) {
                    if (sp3c_i_other==sp3c[mem_sp3c[sp3c[i][j2]][j]][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other ring of 5A_i distinct from whole 5A_j

                m=0;
                for (k=0; k<5; k++) {
                    if (sp3c_j_other==sp3c[i][k]) {
                        m++;
                        break;
                    }
                }
                if (m!=0 && k!=5) continue; // other ring of 5A_j distinct from whole 5A_i

                if (n7K[f]==m7K) {
                    hc7K=resize_2D_int(hc7K,m7K,m7K+incrStatic,clusSize,-1);
                    m7K=m7K+incrStatic;
                }
                // Now we have found the 7K cluster

                hc7K[n7K[f]][0]=scom;
                hc7K[n7K[f]][1]=sother[0];
                hc7K[n7K[f]][2]=sother[1];
                hc7K[n7K[f]][3]=sp3_com[0];
                hc7K[n7K[f]][4]=sp3_com[1];
                hc7K[n7K[f]][5]=sp3c_i_other;
                hc7K[n7K[f]][6]=sp3c_j_other;

                quickSort(&hc7K[n7K[f]][1],2);
                quickSort(&hc7K[n7K[f]][3],2);
                quickSort(&hc7K[n7K[f]][5],2);

                Cluster_Write_7K(f, ach);
            }
        }
    }
    for (i=0; i<N; ++i) s7K[i]=ach[i];
    free(ach);
}

void Cluster_Write_7K(int f, char *ach) {// hc7K key: (scom, sother, ring_com, ring_other)
    ach[hc7K[n7K[f]][0]] = 'O';
    ach[hc7K[n7K[f]][1]] = 'O';
    ach[hc7K[n7K[f]][2]] = 'O';
    if (ach[hc7K[n7K[f]][3]] == 'C') ach[hc7K[n7K[f]][3]] = 'B';
    if (ach[hc7K[n7K[f]][4]] == 'C') ach[hc7K[n7K[f]][4]] = 'B';
    if (ach[hc7K[n7K[f]][5]] == 'C') ach[hc7K[n7K[f]][5]] = 'B';
    if (ach[hc7K[n7K[f]][6]] == 'C') ach[hc7K[n7K[f]][6]] = 'B';

    ++n7K[f];
}

void Clusters_Get8A_D2d(int f)  { // Detect 8A D2d clusters
    int unc[2];
    int com[4];
    int i, j, j2, k, l, m;
    int cnt;
    int flg;
    int break_out;
    int trial[8];
    char *ach, errMsg[1000];
    int clusSize=8;
    int *used_sp5b;

    unc[0]=unc[1]=com[0]=com[1]=com[2]=com[3]=-1;
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get8A_D2d(): ach[] malloc out of memory\n");  Error(errMsg); }
    for (i=0; i<N; ++i) ach[i] = 'C';

    used_sp5b=malloc(nsp5b[f]*sizeof(int)); if (used_sp5b==NULL) { sprintf(errMsg,"Clusters_Get8A_D2d(): used_sp5b[] malloc out of memory\n");  Error(errMsg); }

    for (i=0; i<nsp5b[f]-1; ++i) {  // loop over all sp5b_i
        memset(used_sp5b, 0, nsp5b[f]*sizeof(*used_sp5b));
        used_sp5b[i]=1;
        for (j2=0; j2<5; ++j2) {
            for (j=0; j<nmem_sp5b[sp5b[i][j2]]; ++j) {  // loop over all sp5b_j
                if (mem_sp5b[sp5b[i][j2]][j]<=i) continue;
                if (used_sp5b[mem_sp5b[sp5b[i][j2]][j]]==1) continue;
                used_sp5b[mem_sp5b[sp5b[i][j2]][j]]=1;
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if(sp5b[i][k] == sp5b[mem_sp5b[sp5b[i][j2]][j]][l]) {
                            if (m<5) com[m]=sp5b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue; // exactly four members of the SP5 rings of sp5b_i and sp5b_j in common

                if (sp5b[i][5] == sp5b[mem_sp5b[sp5b[i][j2]][j]][5]) continue;  // distinct spindles

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (sp5b[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=sp5b[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (sp5b[mem_sp5b[sp5b[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=sp5b[mem_sp5b[sp5b[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A[f]==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=sp5b[i][5];    // build up trial cluster
                trial[1]=sp5b[mem_sp5b[sp5b[i][j2]][j]][5];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5b[i][k],trial[4])==1 && sp5b[i][k]!=trial[4] && sp5b[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=sp5b[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5b[i][k],trial[2])==1 && sp5b[i][k]!=trial[2] && sp5b[i][k]!=trial[4] && sp5b[i][k]!=trial[5]) {
                        trial[6]=sp5b[i][k];
                    }
                    if (Bonds_BondCheck(sp5b[i][k],trial[3])==1 && sp5b[i][k]!=trial[3] && sp5b[i][k]!=trial[4] && sp5b[i][k]!=trial[5]) {
                        trial[7]=sp5b[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A[f]; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A[f]][k]=trial[k];

                    Cluster_Write_8A(f, ach);
                }
            }
        }
    }
    for (i=0; i<nsp5c[f] - 1; ++i) {    // loop over all 7A_i
        for (j2=5; j2<6; ++j2) {
            for (j=0; j<nmem_sp5c[sp5c[i][j2]]; ++j) {  // loop over all 7A_j
                if (mem_sp5c[sp5c[i][j2]][j]<=i) continue;
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if (sp5c[i][k] == sp5c[mem_sp5c[sp5c[i][j2]][j]][l]) {
                            if (m<5) com[m]=sp5c[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue;     // exactly four members of the SP5 rings of 7A_i and 7A_j in common

                flg = sp5c[i][5] == sp5c[mem_sp5c[sp5c[i][j2]][j]][5] && sp5c[i][6] == sp5c[mem_sp5c[sp5c[i][j2]][j]][6];
                flg = flg || (sp5c[i][5] == sp5c[mem_sp5c[sp5c[i][j2]][j]][6] && sp5c[i][6] == sp5c[mem_sp5c[sp5c[i][j2]][j]][5]);
                if (flg!=1) continue; // both spindles common

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (sp5c[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=sp5c[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (sp5c[mem_sp5c[sp5c[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=sp5c[mem_sp5c[sp5c[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A[f]==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=sp5c[i][5];    // build up trial cluster
                trial[1]=sp5c[i][6];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5c[i][k],trial[4])==1 && sp5c[i][k]!=trial[4] && sp5c[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=sp5c[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5c[i][k],trial[2])==1 && sp5c[i][k]!=trial[2] && sp5c[i][k]!=trial[4] && sp5c[i][k]!=trial[5]) {
                        trial[6]=sp5c[i][k];
                    }
                    if (Bonds_BondCheck(sp5c[i][k],trial[3])==1 && sp5c[i][k]!=trial[3] && sp5c[i][k]!=trial[4] && sp5c[i][k]!=trial[5]) {
                        trial[7]=sp5c[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A[f]; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A[f]][k]=trial[k];

                    Cluster_Write_8A(f, ach);
                }
            }
        }
    }
    for (i=0; i<nsp5b[f]; ++i) {    // loop over all sp5b_i
        for (j2=5; j2<6; ++j2) {
            for (j=0; j<nmem_sp5c[sp5b[i][j2]]; ++j) {  // loop over all 7A_j
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if (sp5b[i][k] == sp5c[mem_sp5c[sp5b[i][j2]][j]][l]) {
                            if (m<5) com[m]=sp5b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue; // exactly four members of the SP5 rings of sp5b_i and 7A_j in common

                flg = sp5b[i][5] == sp5c[mem_sp5c[sp5b[i][j2]][j]][5] || sp5b[i][5] == sp5c[mem_sp5c[sp5b[i][j2]][j]][6];
                if (flg!=1) continue;   // sp5b_i spindle common with one of 7A_j spindles

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (sp5b[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=sp5b[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (sp5c[mem_sp5c[sp5b[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=sp5c[mem_sp5c[sp5b[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A[f]==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=sp5c[mem_sp5c[sp5b[i][j2]][j]][5]; // build up trial cluster
                trial[1]=sp5c[mem_sp5c[sp5b[i][j2]][j]][6];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5b[i][k],trial[4])==1 && sp5b[i][k]!=trial[4] && sp5b[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=sp5b[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5b[i][k],trial[2])==1 && sp5b[i][k]!=trial[2] && sp5b[i][k]!=trial[4] && sp5b[i][k]!=trial[5]) {
                        trial[6]=sp5b[i][k];
                    }
                    if (Bonds_BondCheck(sp5b[i][k],trial[3])==1 && sp5b[i][k]!=trial[3] && sp5b[i][k]!=trial[4] && sp5b[i][k]!=trial[5]) {
                        trial[7]=sp5b[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A[f]; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A[f]][k]=trial[k];
                    Cluster_Write_8A(f, ach);
                }
            }
        }
    }

    for(i=0; i<N; ++i) s8A[i]=ach[i];
    free(ach);
    free(used_sp5b);
}

void Cluster_Write_8A(int f, char *ach) {// hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
    if (ach[hc8A[n8A[f]][0]] == 'C') ach[hc8A[n8A[f]][0]] = 'B';
    if (ach[hc8A[n8A[f]][1]] == 'C') ach[hc8A[n8A[f]][1]] = 'B';
    if (ach[hc8A[n8A[f]][2]] == 'C') ach[hc8A[n8A[f]][2]] = 'B';
    if (ach[hc8A[n8A[f]][3]] == 'C') ach[hc8A[n8A[f]][3]] = 'B';
    if (ach[hc8A[n8A[f]][4]] == 'C') ach[hc8A[n8A[f]][4]] = 'B';
    if (ach[hc8A[n8A[f]][5]] == 'C') ach[hc8A[n8A[f]][5]] = 'B';
    ach[hc8A[n8A[f]][6]] = 'O';
    ach[hc8A[n8A[f]][7]] = 'O';

    ++n8A[f];
}

void Clusters_Get8B_Cs(int f) { // Detect 8B Cs clusters
    int i;
    char *ach, errMsg[1000];
    int clusSize=8;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get8B_Cs(): ach[] malloc out of memory\n");   Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = 'C';

    for (i=0; i<nsp5c[f]; ++i) {    // loop over all 7A_i
        Clusters_8B_loop(f, i, ach, clusSize, sp5c[i][5], sp5c[i][6]);
        Clusters_8B_loop(f, i, ach, clusSize, sp5c[i][6], sp5c[i][5]);
    }
    for(i=0; i<N; ++i) s8B[i]=ach[i];
    free(ach);
}

void Clusters_8B_loop(int f, int i, char *ach, int clusSize, int primary_spindle, int secondary_spindle) {

    int j, k, l, m;
    int n1, nbs, unc[3];
    int break_out;

    for (j=0; j < cnb[primary_spindle]; ++j) { // loop over all j particles bonded to first spindle of 7A_i
        n1 = bNums[primary_spindle][j];
        for (k=0; k<5; ++k) if (n1 == sp5c[i][k]) break;
        if (k<5) continue;
        if (n1 == secondary_spindle) continue; // now is n1 bonded to sp5
        nbs = 0; // number of bonds
        for (k=0; k<5; ++k) if (Bonds_BondCheck(n1, sp5c[i][k])) ++nbs;
        if (nbs != 2) continue;

        // Now we have found the 8B Cs cluster
        if (n8B[f]==m8B) {
            hc8B=resize_2D_int(hc8B,m8B,m8B+incrStatic,clusSize,-1);
            m8B=m8B+incrStatic;
        }
        l=0;
        m=3;
        break_out=0;
        for (k=0; k<5; ++k) {
            if (Bonds_BondCheck(n1, sp5c[i][k])) {
                if (m==5) {
                    break_out=1;
                    break;
                }
                hc8B[n8B[f]][m]=sp5c[i][k];
                m++;
            }
            else {
                if (l==3) {
                    break_out=1;
                    break;
                }
                unc[l]=sp5c[i][k];
                l++;
            }
        }
        if (break_out==1 || m<5 || l<3) continue;

        quickSort(&hc8B[n8B[f]][3],2);
        for (k=0;k<3;k++) hc8B[n8B[f]][k]=unc[k];
        quickSort(&hc8B[n8B[f]][0],3);

        hc8B[n8B[f]][5]=secondary_spindle;
        hc8B[n8B[f]][6]=primary_spindle;
        hc8B[n8B[f]][7]=n1;


        Cluster_Write_8B(f, ach);
    }
}

void Cluster_Write_8B(int f, char *ach) {
    // hc8B key: (SP5_to_4, SP5_to_0/2, SP5_to_3, SP5_to_n1(lower), SP5_to_n1(greater), s, s_to_n1, n1)
    if (ach[hc8B[n8B[f]][7]] == 'C') ach[hc8B[n8B[f]][7]] = 'B';
    if (ach[hc8B[n8B[f]][0]] == 'C') ach[hc8B[n8B[f]][0]] = 'B';
    if (ach[hc8B[n8B[f]][1]] == 'C') ach[hc8B[n8B[f]][1]] = 'B';
    if (ach[hc8B[n8B[f]][2]] == 'C') ach[hc8B[n8B[f]][2]] = 'B';
    if (ach[hc8B[n8B[f]][3]] == 'C') ach[hc8B[n8B[f]][3]] = 'B';
    if (ach[hc8B[n8B[f]][4]] == 'C') ach[hc8B[n8B[f]][4]] = 'B';
    ach[hc8B[n8B[f]][5]] = 'O';
    ach[hc8B[n8B[f]][6]] = 'O';
    ++n8B[f];
}


void Clusters_Get8K(int f) {    // Detect 8K clusters
    int i, j, j2, k, l, m, n;
    int cp[2], unc[3], scom, sother[2];
    int trial[8];
    char *ach, errMsg[1000];
    int clusSize=8;

    cp[0]=cp[1]=unc[0]=unc[1]=unc[2]=scom=sother[0]=sother[1]=-1;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get8K(): ach[] malloc out of memory\n");  Error(errMsg); }
    for (i=0; i<N; ++i) ach[i] = 'C';

    for (i=0; i<nsp3c[f]-2; ++i) {  // loop over all sp3c_i
        for (j2=0; j2<3; j2++) {    // loop over all particles in SP3 ring of sp3c_i
            for (j=0; j<nmem_sp3c[sp3c[i][j2]]-1; ++j) { // loop over all sp3c_j which sp3c[i][j2] is a member of
                if (mem_sp3c[sp3c[i][j2]][j]<=i) continue;  // check not used mem_sp3c[sp3c[i][j2]][j] before

                m = 0;  // check j2 from SP3 ring of sp3c_i is in SP3 ring of mem_sp3c[sp3c[i][j2]][j]
                for(k=0; k<3; ++k) {
                    if (sp3c[i][j2]==sp3c[mem_sp3c[sp3c[i][j2]][j]][k]) {
                        cp[0]=sp3c[i][j2];
                        m++;
                    }
                }
                if (m!=1) continue;

                m = 0;  // find extra 1 particle from SP3 ring of sp3c_i which is also in SP3 ring of mem_sp3c[sp3c[i][j2]][j]
                for(k=0; k<3; ++k) {
                    if (k==j2) continue;    // don't check j2 again
                    if (sp3c[i][k]<sp3c[i][j2]) continue; // will have found before or do not find after when using different j2 particle
                    for(l=0; l<3; ++l) {
                        if (sp3c[i][k]==sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                            cp[1]=sp3c[i][k];
                            m++;
                        }
                    }
                }
                if (m!=1) continue;

                m = 0;  // check exactly one common spindle between sp3c_i and mem_sp3c[sp3c[i][j2]][j]
                for(k=3; k<5; ++k) {
                    for(l=3; l<5; ++l) {
                        if (sp3c[i][k]==sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                            scom=sp3c[i][k];
                            m++;
                        }
                    }
                }
                if (m!=1) continue;

                if (sp3c[i][3]==scom) sother[0]=sp3c[i][4];
                else sother[0]=sp3c[i][3];
                if (sp3c[mem_sp3c[sp3c[i][j2]][j]][3]==scom) sother[1]=sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                else sother[1]=sp3c[mem_sp3c[sp3c[i][j2]][j]][3];

                k=0;    // find particle from SP3 ring of sp3c_i which is common with 5A mem_sp3c[sp3c[i][j2]][j]
                for(l=0; l<3; ++l) {
                    if (sp3c[i][l]==cp[0] || sp3c[i][l]==cp[1]) continue;
                    for(m=0; m<5; ++m) {
                        if (sp3c[i][l]==sp3c[mem_sp3c[sp3c[i][j2]][j]][m]) break;
                    }
                    if (m==5) {
                        if (k==1) {
                            k++;
                            break;
                        }
                        unc[0]=sp3c[i][l];
                        k++;
                    }
                }
                if (k!=1) continue;

                k=0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                for(l=0; l<3; ++l) {
                    if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l]==cp[0] || sp3c[mem_sp3c[sp3c[i][j2]][j]][l]==cp[1]) continue;
                    for(m=0; m<5; ++m) {
                        if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l]==sp3c[i][m]) break;
                    }
                    if (m==5) {
                        if (k==1) {
                            k++;
                            break;
                        }
                        unc[1]=sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                        k++;
                    }
                }
                if (k!=1) continue;

                for (k=j+1; k<nmem_sp3c[sp3c[i][j2]]; ++k) {
                    if (mem_sp3c[sp3c[i][j2]][k]<=i) continue;  // higher index of sp3c cluster than i
                    if (mem_sp3c[sp3c[i][j2]][k]<=mem_sp3c[sp3c[i][j2]][j]) continue;   // higher index of sp3c cluster than mem_sp3c[sp3c[i][j2]][j]

                    n = 0;  // check common SP3 ring particles are exactly cp
                    for(l=0; l<3; ++l) {
                        for(m=0; m<2; ++m) {
                            if (cp[m]==sp3c[mem_sp3c[sp3c[i][j2]][k]][l]) {
                                n++;
                            }
                        }
                    }
                    if (n!=2) continue;

                    n = 0;  // check spindles are exactly sother
                    for(l=3; l<5; ++l) {
                        for(m=0; m<2; ++m) {
                            if (sother[m]==sp3c[mem_sp3c[sp3c[i][j2]][k]][l]) {
                                n++;
                            }
                        }
                    }
                    if (n!=2) continue;

                    n=0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                    for(l=0; l<3; ++l) {
                        if (sp3c[mem_sp3c[sp3c[i][j2]][k]][l]==cp[0] || sp3c[mem_sp3c[sp3c[i][j2]][k]][l]==cp[1]) continue;
                        for(m=0; m<5; ++m) {
                            if (sp3c[mem_sp3c[sp3c[i][j2]][k]][l]==sp3c[i][m] || sp3c[mem_sp3c[sp3c[i][j2]][k]][l]==sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) break;
                        }
                        if (m==5) {
                            if (n==1) {
                                n++;
                                break;
                            }
                            unc[2]=sp3c[mem_sp3c[sp3c[i][j2]][k]][l];
                            n++;
                        }
                    }
                    if (n!=1) continue;

                    // Now we have found the 8K cluster
                    if (n8K[f]==m8K) {
                        hc8K=resize_2D_int(hc8K,m8K,m8K+incrStatic,clusSize,-1);
                        m8K=m8K+incrStatic;
                    }
                    // hc8K key: (SP3_common_1, SP3_common_2, spindle_1, spindle_2, spindle_3, other_SP3_1, other_SP3_2, other_SP3_3)
                    trial[0]=cp[0];
                    trial[1]=cp[1];
                    trial[2]=scom;
                    trial[3]=sother[0];
                    trial[4]=sother[1];
                    trial[5]=unc[0];
                    trial[6]=unc[1];
                    trial[7]=unc[2];

                    quickSort(&trial[0],2);
                    quickSort(&trial[2],3);
                    quickSort(&trial[5],3);
                    for (l=0; l<clusSize; l++) hc8K[n8K[f]][l]=trial[l];

                    Cluster_Write_8K(f, ach);
                }
            }
        }
    }

    for(i=0; i<N; ++i) s8K[i]=ach[i];
    free(ach);
}

void Cluster_Write_8K(int f, char *ach) {
    if(ach[hc8K[n8K[f]][2]] == 'C') ach[hc8K[n8K[f]][2]] = 'B';
    if(ach[hc8K[n8K[f]][3]] == 'C') ach[hc8K[n8K[f]][3]] = 'B';
    if(ach[hc8K[n8K[f]][4]] == 'C') ach[hc8K[n8K[f]][4]] = 'B';
    if(ach[hc8K[n8K[f]][5]] == 'C') ach[hc8K[n8K[f]][5]] = 'B';
    if(ach[hc8K[n8K[f]][6]] == 'C') ach[hc8K[n8K[f]][6]] = 'B';
    if(ach[hc8K[n8K[f]][7]] == 'C') ach[hc8K[n8K[f]][7]] = 'B';
    ach[hc8K[n8K[f]][0]] = 'O';
    ach[hc8K[n8K[f]][1]] = 'O';
    ++n8K[f];
}

void Clusters_Get9A_D3h(int f) {    // Detect 9A D3h clusters
    int i, j, j2, k, l, m, n;
    int db[2], ob[4];
    int flg;
    char *ach, errMsg[1000];
    int clusSize=9;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get9A_D3h(): ach[] malloc out of memory\n");  Error(errMsg); }
    for (i=0; i<N; ++i) ach[i] = 'C';

    for (i=0; i<nsp4b[f]-2; ++i) {  // loop over all sp4b_i
        for (j2=0; j2<1; j2++) {
            for (j=0; j<nmem_sp4b[sp4b[i][j2]]; ++j) { // loop over all sp4b_j
                if (mem_sp4b[sp4b[i][j2]][j]<=i) continue;
                if (sp4b[i][4] == sp4b[mem_sp4b[sp4b[i][j2]][j]][4]) continue;  // if spindles common continue
                if (Bonds_BondCheck(sp4b[i][4], sp4b[mem_sp4b[sp4b[i][j2]][j]][4])) continue;   // if spindles bonded continue
                m = 0;
                for(k=0; k<4; ++k) {
                    for(l=0; l<4; ++l) {
                        if(sp4b[i][k] == sp4b[mem_sp4b[sp4b[i][j2]][j]][l]){
                            if(m<2) db[m] = sp4b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=2) continue;         // two common particles between SP4 rings of sp4b_i and sp4b_j

                m = 0;
                for (k=0; k<4; ++k) {
                    if(sp4b[i][k] == db[0] || sp4b[i][k] == db[1]) continue;  // find particles in SP4 ring of sp4b_i not common to SP4 ring of sp4b_j
                    for(l=0; l<4; ++l){
                        if(sp4b[mem_sp4b[sp4b[i][j2]][j]][l] == db[0] || sp4b[mem_sp4b[sp4b[i][j2]][j]][l] == db[1]) continue;    // find particles in SP4 ring of sp4b_j not common to SP4 ring of sp4b_i
                        if(Bonds_BondCheck(sp4b[i][k], sp4b[mem_sp4b[sp4b[i][j2]][j]][l])) {    // check non-common SP4 ring particles from sp4b_i and sp4b_j are bonded
                            if(m<4) ob[m] = sp4b[i][k];
                            ++m;
                            if(m<4) ob[m] = sp4b[mem_sp4b[sp4b[i][j2]][j]][l];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue;
                // POSSIBLE IMPROVEMENT!! could make detection faster here by looping over sp4b clusters bonded to uncommon particles to sp4b_i
                for(k=mem_sp4b[sp4b[i][j2]][j]+1; k<nsp4b[f]; k++) {    // loop over all sp4b_k
                    // ERROR!! need to check spindle of sp4b_k distinct from spindles of sp4b_i and sp4b_j and no bonds between these three particles
                    n = 0;
                    for(l=0; l<4; ++l){
                        for(m=0; m<4; ++m){
                            if(sp4b[k][l] == ob[m]){
                                ++n;
                                break;
                            }
                        }
                    }
                    if (n != 4) continue;

                    // Now we have found the 9A D3h cluster
                    if (n9A[f]==m9A) {
                        hc9A=resize_2D_int(hc9A,m9A,m9A+incrStatic,clusSize,-1);
                            m9A=m9A+incrStatic;
                        }
                        if (sp4b[i][4]<sp4b[mem_sp4b[sp4b[i][j2]][j]][4] && sp4b[i][4]<sp4b[k][4]) {
                            hc9A[n9A[f]][6]=sp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[f][i][.] cluster
                            hc9A[n9A[f]][0]=sp4b[i][0];
                            hc9A[n9A[f]][1]=sp4b[i][1];
                            hc9A[n9A[f]][2]=sp4b[i][2];
                            hc9A[n9A[f]][3]=sp4b[i][3];

                            if (Bonds_BondCheck(sp4b[mem_sp4b[sp4b[i][j2]][j]][4],hc9A[n9A[f]][0])) {
                                hc9A[n9A[f]][7]=sp4b[mem_sp4b[sp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[f][j][.] cluster
                                hc9A[n9A[f]][8]=sp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[f][k][.] cluster
                            }
                            else {
                                hc9A[n9A[f]][7]=sp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[f][k][.] cluster
                                hc9A[n9A[f]][8]=sp4b[mem_sp4b[sp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[f][j][.] cluster
                            }

                            for (l=0;l<4;l++) {
                                flg=1;
                                for (m=0;m<4;m++) {
                                    if (ob[l]==hc9A[n9A[f]][m]) {
                                        flg=0;
                                        break;
                                    }
                                }
                                if (flg==1 && Bonds_BondCheck(ob[l],hc9A[n9A[f]][0])) {
                                    hc9A[n9A[f]][4]=ob[l];
                                }
                                else if (flg==1) {
                                    hc9A[n9A[f]][5]=ob[l];
                                }
                            }
                        }

                        else if (sp4b[mem_sp4b[sp4b[i][j2]][j]][4]<sp4b[i][4] && sp4b[mem_sp4b[sp4b[i][j2]][j]][4]<sp4b[k][4]) {
                            hc9A[n9A[f]][6]=sp4b[mem_sp4b[sp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[f][j][.] cluster
                            hc9A[n9A[f]][0]=sp4b[mem_sp4b[sp4b[i][j2]][j]][0];
                            hc9A[n9A[f]][1]=sp4b[mem_sp4b[sp4b[i][j2]][j]][1];
                            hc9A[n9A[f]][2]=sp4b[mem_sp4b[sp4b[i][j2]][j]][2];
                            hc9A[n9A[f]][3]=sp4b[mem_sp4b[sp4b[i][j2]][j]][3];

                            if (Bonds_BondCheck(sp4b[i][4],hc9A[n9A[f]][0])) {
                                hc9A[n9A[f]][7]=sp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[f][i][.] cluster
                                hc9A[n9A[f]][8]=sp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[f][k][.] cluster
                            }
                            else {
                                hc9A[n9A[f]][7]=sp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[f][k][.] cluster
                                hc9A[n9A[f]][8]=sp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[f][i][.] cluster
                            }

                            for (l=0;l<4;l++) {
                                flg=1;
                                for (m=0;m<4;m++) {
                                    if (ob[l]==hc9A[n9A[f]][m]) {
                                        flg=0;
                                        break;
                                    }
                                }
                                if (flg==1 && Bonds_BondCheck(ob[l],hc9A[n9A[f]][0])) {
                                    hc9A[n9A[f]][4]=ob[l];
                                }
                                else if (flg==1) {
                                    hc9A[n9A[f]][5]=ob[l];
                                }
                            }
                        }

                        else {
                            hc9A[n9A[f]][6]=sp4b[k][4]; // spindle of sp4 ring ob[0],ob[2],ob[1],ob[3] of sp4b[f][k][.] cluster
                            hc9A[n9A[f]][0]=sp4b[k][0];
                            hc9A[n9A[f]][1]=sp4b[k][1];
                            hc9A[n9A[f]][2]=sp4b[k][2];
                            hc9A[n9A[f]][3]=sp4b[k][3];

                            if (Bonds_BondCheck(sp4b[i][4],hc9A[n9A[f]][0])) {
                                hc9A[n9A[f]][7]=sp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[f][i][.] cluster
                                hc9A[n9A[f]][8]=sp4b[mem_sp4b[sp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[f][j][.] cluster
                            }
                            else {
                                hc9A[n9A[f]][7]=sp4b[mem_sp4b[sp4b[i][j2]][j]][4];  // spindle of sp4 ring db[0],db[1],ob[1],ob[3] of sp4b[f][j][.] cluster
                                hc9A[n9A[f]][8]=sp4b[i][4]; // spindle of sp4 ring db[0],db[1],ob[0],ob[2] of sp4b[f][i][.] cluster
                            }

                            if (Bonds_BondCheck(db[0],hc9A[n9A[f]][0])) {
                                hc9A[n9A[f]][4]=db[0];
                                hc9A[n9A[f]][5]=db[1];
                            }
                            else {
                                hc9A[n9A[f]][4]=db[1];
                                hc9A[n9A[f]][5]=db[0];
                            }
                        }
                    Cluster_Write_9A(f, ach);
                    break;
                }
            }
        }
    }

    for(i=0; i<N; ++i) s9A[i]=ach[i];
    free(ach);
}

void Cluster_Write_9A(int f, char *ach) {// hc9A key: (SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_to_0_in_SP4_lowest_s, SP4_to_1_in_SP4_lowest_s, s_lowest, s_to_0_in_SP4_lowest_s,s_to_2_in_SP4_lowest_s)
    if(ach[hc9A[n9A[f]][0]] == 'C') ach[hc9A[n9A[f]][0]] = 'B';
    if(ach[hc9A[n9A[f]][1]] == 'C') ach[hc9A[n9A[f]][1]] = 'B';
    if(ach[hc9A[n9A[f]][3]] == 'C') ach[hc9A[n9A[f]][3]] = 'B';
    if(ach[hc9A[n9A[f]][4]] == 'C') ach[hc9A[n9A[f]][4]] = 'B';
    if(ach[hc9A[n9A[f]][6]] == 'C') ach[hc9A[n9A[f]][6]] = 'B';
    if(ach[hc9A[n9A[f]][7]] == 'C') ach[hc9A[n9A[f]][7]] = 'B';
    ach[hc9A[n9A[f]][2]] = 'O';
    ach[hc9A[n9A[f]][5]] = 'O';
    ach[hc9A[n9A[f]][8]] = 'O';
    ++n9A[f];
}

void Clusters_Get9B_10B_11B_11E_12D(int f) {    // Detect 9B, 10B, 11A, 11E & 12D
    int sp1, sp2i, sp2j;
    int sp5com[2];
    int i, j, k, l, m;
    int flg, fb1, fb2;
    char *ach1, *ach1_cen, *ach1_shell, *ach2, *ach2_cen, *ach2_shell, *ach3, *ach4, *ach5, *ach5_cen, *ach5_shell, *ach6, *ach6_cen, *ach6_shell, errMsg[1000];
    int clusSize=9;

    sp1=sp2i=sp2j=-1;
    ach1=malloc(N*sizeof(char));    if (ach1==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach1[] malloc out of memory\n");    Error(errMsg); }
    ach1_cen=malloc(N*sizeof(char));    if (ach1_cen==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach1_cen[] malloc out of memory\n");    Error(errMsg); }
    ach1_shell=malloc(N*sizeof(char));  if (ach1_shell==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach1_shell[] malloc out of memory\n");    Error(errMsg); }
    ach2=malloc(N*sizeof(char));    if (ach2==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach2[] malloc out of memory\n");    Error(errMsg); }
    ach2_cen=malloc(N*sizeof(char));    if (ach2_cen==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach2_cen[] malloc out of memory\n");    Error(errMsg); }
    ach2_shell=malloc(N*sizeof(char));  if (ach2_shell==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach2_shell[] malloc out of memory\n");    Error(errMsg); }
    ach3=malloc(N*sizeof(char));    if (ach3==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach3[] malloc out of memory\n");    Error(errMsg); }
    ach4=malloc(N*sizeof(char));    if (ach4==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach4[] malloc out of memory\n");    Error(errMsg); }
    ach5=malloc(N*sizeof(char));    if (ach5==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach5[] malloc out of memory\n");    Error(errMsg); }
    ach5_cen=malloc(N*sizeof(char));    if (ach5_cen==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach5_cen[] malloc out of memory\n");    Error(errMsg); }
    ach5_shell=malloc(N*sizeof(char));  if (ach5_shell==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach5_shell[] malloc out of memory\n");    Error(errMsg); }
    ach6=malloc(N*sizeof(char));    if (ach6==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach6[] malloc out of memory\n");    Error(errMsg); }
    ach6_cen=malloc(N*sizeof(char));    if (ach6_cen==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach6_cen[] malloc out of memory\n");    Error(errMsg); }
    ach6_shell=malloc(N*sizeof(char));  if (ach6_shell==NULL) { sprintf(errMsg,"Clusters_Get9B_10B_11B_11E_12D(): ach6_shell[] malloc out of memory\n");    Error(errMsg); }
    for(i=0; i<N; ++i) ach1[i] = ach1_cen[i] = ach1_shell[i] = ach2[i] = ach2_cen[i] = ach2_shell[i] = ach3[i] = ach4[i] = ach5[i] = ach5_cen[i] = ach5_shell[i] = ach6[i] = ach6_cen[i] = ach6_shell[i] = 'C';

    for (i=0; i<nsp5c[f]-1; ++i) {  // loop over all 7A_i
        // POSSIBLE IMPROVEMENT!! - 2 loops: over all 7A clusters which each spindle is in
        for (j=i+1; j<nsp5c[f]; ++j) {  // loop over all 7A_j
            flg = 0;
            if (sp5c[i][5] == sp5c[j][5] && sp5c[i][6] != sp5c[j][6]) {
                if (Bonds_BondCheck(sp5c[i][6], sp5c[j][6])) { // spindle particles arranged
                    flg = 1;
                    sp1 = sp5c[i][5];   // s_com common spindle
                    sp2i = sp5c[i][6];  // 2nd spindle particle of cluster 7A_i
                    sp2j = sp5c[j][6];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (sp5c[i][6] == sp5c[j][6] && sp5c[i][5] != sp5c[j][5]) {
                if (Bonds_BondCheck(sp5c[i][5], sp5c[j][5])) {
                    flg = 1;
                    sp1 = sp5c[i][6];   // s_com common spindle
                    sp2i = sp5c[i][5];  // 2nd spindle particle of cluster 7A_i
                    sp2j = sp5c[j][5];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (sp5c[i][5] == sp5c[j][6] && sp5c[i][6] != sp5c[j][5]) {
                if (Bonds_BondCheck(sp5c[i][6], sp5c[j][5])) {
                    flg = 1;
                    sp1 = sp5c[i][5];   // s_com common spindle
                    sp2i = sp5c[i][6];  // 2nd spindle particle of cluster 7A_i
                    sp2j = sp5c[j][5];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (sp5c[i][6] == sp5c[j][5] && sp5c[i][5] != sp5c[j][6]) {
                if (Bonds_BondCheck(sp5c[i][5], sp5c[j][6])) {
                    flg = 1;
                    sp1 = sp5c[i][6];   // s_com common spindle
                    sp2i = sp5c[i][5];  // 2nd spindle particle of cluster 7A_i
                    sp2j = sp5c[j][6];  // 2nd spindle particle of cluster 7A_j
                }
            }
            if (flg==0) continue;

            fb1 = fb2 = 1;  // ensure the two distinct spindle particles are part of the other SP5 ring
            for (k=0; k<5; ++k) {
                if (sp2i == sp5c[j][k]) fb1 = 0;
                if (sp2j == sp5c[i][k]) fb2 = 0;
            }
            if (fb1 || fb2) continue;

            m = 0;  // check for two common SP5 particles
            for (k=0; k<5; ++k) {
                for (l=0; l<5; ++l) {
                    if (sp5c[i][k] == sp5c[j][l]) {
                        if (m==2) {m++; break; }
                        sp5com[m]=sp5c[i][k];
                        ++m;
                    }
                }
            }
            if (m!=2) continue;

            // Now we have found the 9B C2v cluster
            if (n9B[f]==m9B) {
                hc9B=resize_2D_int(hc9B,m9B,m9B+incrStatic,clusSize,-1);
                m9B=m9B+incrStatic;
            }
            if (sp5com[0]<sp5com[1]) {
                hc9B[n9B[f]][4]=sp5com[0];
                hc9B[n9B[f]][5]=sp5com[1];
            }
            else {
                hc9B[n9B[f]][4]=sp5com[1];
                hc9B[n9B[f]][5]=sp5com[0];
            }

            if (sp2i<sp2j) {
                hc9B[n9B[f]][6]=sp2i;
                hc9B[n9B[f]][7]=sp2j;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5c[i][k],hc9B[n9B[f]][4]) && sp5c[i][k]!=hc9B[n9B[f]][7] && sp5c[i][k]!=hc9B[n9B[f]][4]) {
                        hc9B[n9B[f]][0]=sp5c[i][k];
                    }
                    if (Bonds_BondCheck(sp5c[i][k],hc9B[n9B[f]][5]) && sp5c[i][k]!=hc9B[n9B[f]][7] && sp5c[i][k]!=hc9B[n9B[f]][5]) {
                        hc9B[n9B[f]][1]=sp5c[i][k];
                    }
                    if (Bonds_BondCheck(sp5c[j][k],hc9B[n9B[f]][4]) && sp5c[j][k]!=hc9B[n9B[f]][6] && sp5c[j][k]!=hc9B[n9B[f]][4]) {
                        hc9B[n9B[f]][2]=sp5c[j][k];
                    }
                    if (Bonds_BondCheck(sp5c[j][k],hc9B[n9B[f]][5]) && sp5c[j][k]!=hc9B[n9B[f]][6] && sp5c[j][k]!=hc9B[n9B[f]][5]) {
                        hc9B[n9B[f]][3]=sp5c[j][k];
                    }
                }
            }
            else {
                hc9B[n9B[f]][6]=sp2j;
                hc9B[n9B[f]][7]=sp2i;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(sp5c[j][k],hc9B[n9B[f]][4]) && sp5c[j][k]!=hc9B[n9B[f]][7] && sp5c[j][k]!=hc9B[n9B[f]][4]) {
                        hc9B[n9B[f]][0]=sp5c[j][k];
                    }
                    if (Bonds_BondCheck(sp5c[j][k],hc9B[n9B[f]][5]) && sp5c[j][k]!=hc9B[n9B[f]][7] && sp5c[j][k]!=hc9B[n9B[f]][5]) {
                        hc9B[n9B[f]][1]=sp5c[j][k];
                    }
                    if (Bonds_BondCheck(sp5c[i][k],hc9B[n9B[f]][4]) && sp5c[i][k]!=hc9B[n9B[f]][6] && sp5c[i][k]!=hc9B[n9B[f]][4]) {
                        hc9B[n9B[f]][2]=sp5c[i][k];
                    }
                    if (Bonds_BondCheck(sp5c[i][k],hc9B[n9B[f]][5]) && sp5c[i][k]!=hc9B[n9B[f]][6] && sp5c[i][k]!=hc9B[n9B[f]][5]) {
                        hc9B[n9B[f]][3]=sp5c[i][k];
                    }
                }
            }
            hc9B[n9B[f]][8]=sp1;

            // hc9B key: (SP5_lowerd_to_4, SP5_lowerd_to_5, SP5_higherd_to_4, SP5_higherd_to_5, SP5_i_j_com_lower, SP5_i_j_com_higher, sp5c_d1_lower, sp5c_d2_higher, s_com)
            if (ach1[hc9B[n9B[f]][0]] == 'C') ach1[hc9B[n9B[f]][0]] = ach1_shell[hc9B[n9B[f]][0]] = 'B';
            if (ach1[hc9B[n9B[f]][1]] == 'C') ach1[hc9B[n9B[f]][1]] = ach1_shell[hc9B[n9B[f]][1]] = 'B';
            if (ach1[hc9B[n9B[f]][2]] == 'C') ach1[hc9B[n9B[f]][2]] = ach1_shell[hc9B[n9B[f]][2]] = 'B';
            if (ach1[hc9B[n9B[f]][3]] == 'C') ach1[hc9B[n9B[f]][3]] = ach1_shell[hc9B[n9B[f]][3]] = 'B';
            if (ach1[hc9B[n9B[f]][4]] == 'C') ach1[hc9B[n9B[f]][4]] = ach1_shell[hc9B[n9B[f]][4]] = 'B';
            if (ach1[hc9B[n9B[f]][5]] == 'C') ach1[hc9B[n9B[f]][5]] = ach1_shell[hc9B[n9B[f]][5]] = 'B';
            ach1[hc9B[n9B[f]][6]] = ach1_shell[hc9B[n9B[f]][6]] = 'O';
            ach1[hc9B[n9B[f]][7]] = ach1_shell[hc9B[n9B[f]][7]] = 'O';
            ach1[hc9B[n9B[f]][8]] = ach1_cen[hc9B[n9B[f]][8]] = 'O';

            if (do10B==1) Clusters_Get10B_C3v(f,i,j,ach2,ach2_cen,ach2_shell,ach6,ach6_cen,ach6_shell);
            if (do11B==1) {
                if (Clusters_Get11B_C2v(f,ach5,ach5_cen,ach5_shell)) {
                    ach5_cen[hc9B[n9B[f]][8]] = 'O';
                    ++n11B[f];
                }
            }
            if (do11E==1) Clusters_Get11E_12D(f,i,j,sp1,sp2i,sp2j,ach3,ach4);

            ++n9B[f];
        }
    }

    for(i=0; i<N; ++i) {
        s9B[i]=ach1[i];
        s9B_cen[i]=ach1_cen[i];
        s9B_shell[i]=ach1_shell[i];
        s10B[i]=ach2[i];
        s10B_cen[i]=ach2_cen[i];
        s10B_shell[i]=ach2_shell[i];
        s11B[i]=ach5[i];
        s11B_cen[i]=ach5_cen[i];
        s11B_shell[i]=ach5_shell[i];
        s11E[i]=ach3[i];
        s11W[i]=ach6[i];
        s11W_cen[i]=ach6_cen[i];
        s11W_shell[i]=ach6_shell[i];
        s12D[i]=ach4[i];
    }
    free(ach1);
    free(ach1_cen);
    free(ach1_shell);
    free(ach2);
    free(ach2_cen);
    free(ach2_shell);
    free(ach3);
    free(ach4);
    free(ach5);
    free(ach5_cen);
    free(ach5_shell);
    free(ach6);
    free(ach6_cen);
    free(ach6_shell);
}

void Clusters_Get10B_C3v(int f, int i, int j, char *ach, char *ach_cen, char *ach_shell, char *ach1, char *ach1_cen, char *ach1_shell) {        // Return 1 if 9B is also 10B cluster
    int k,l,m;
    int flg1, flg2;
    int trial[10];
    int break_out;
    int clusSize=10;

    for (k=j+1; k<nsp5c[f]; ++k) {  // loop over all 7A_k
        if (sp5c[k][5] == hc9B[n9B[f]][8]) {    // check one spindle of 7A_k is the common spindle of 9B (hc9B[id9B][.] at this point)
            if (Bonds_BondCheck(sp5c[k][6], hc9B[n9B[f]][6])==0) continue;  // check other spindle of 7A_k is bonded to spindle d1 of 9B
            if (Bonds_BondCheck(sp5c[k][6], hc9B[n9B[f]][7])==0) continue;  // check other spindle of 7A_k is bonded to spindle d2 of 9B

            flg1=0;
            flg2=0;
            for (l=0;l<5;l++) {
                if (sp5c[k][l]==hc9B[n9B[f]][6]) {
                    flg1=1;
                    continue;
                }
                if (sp5c[k][l]==hc9B[n9B[f]][7]) {
                    flg2=1;
                    continue;
                }
            }
            if (flg1==0 || flg2==0) continue;
            trial[6]=hc9B[n9B[f]][6];
            trial[7]=hc9B[n9B[f]][7];
            trial[8]=sp5c[k][6];
            trial[9]=hc9B[n9B[f]][8];

            m=0;
            break_out=0;
            for (l=0;l<6;l++) {
                if (hc9B[n9B[f]][l]==sp5c[k][6]) continue;
                if (m==5) {
                    m++;
                    break_out=1;
                    break;
                }
                trial[m]=hc9B[n9B[f]][l];
                m++;
            }
            if (break_out==1 || m!=5) continue;

            break_out=0;
            for (l=0;l<5;l++) {
                if (sp5c[k][l]==hc9B[n9B[f]][6]) continue;
                if (sp5c[k][l]==hc9B[n9B[f]][7]) continue;
                for (m=0;m<5;m++) {
                    if (sp5c[k][l]==trial[m]) break;
                }
                if (m==5) {
                    trial[5]=sp5c[k][l];
                    break_out++;
                }
            }
            if (break_out!=1) continue;

            if (n10B[f]==m10B) { hc10B=resize_2D_int(hc10B,m10B,m10B+incrStatic,clusSize,-1);
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
            for (l=0;l<10;l++) hc10B[n10B[f]][l]=trial[l];

            // hc10B key: (ordered shell particles, s1, s2, s3 (ordered), s_com)
            if(ach[hc10B[n10B[f]][0]] == 'C') ach[hc10B[n10B[f]][0]] = ach_shell[hc10B[n10B[f]][0]] = 'B';
            if(ach[hc10B[n10B[f]][1]] == 'C') ach[hc10B[n10B[f]][1]] = ach_shell[hc10B[n10B[f]][1]] = 'B';
            if(ach[hc10B[n10B[f]][2]] == 'C') ach[hc10B[n10B[f]][2]] = ach_shell[hc10B[n10B[f]][2]] = 'B';
            if(ach[hc10B[n10B[f]][3]] == 'C') ach[hc10B[n10B[f]][3]] = ach_shell[hc10B[n10B[f]][3]] = 'B';
            if(ach[hc10B[n10B[f]][4]] == 'C') ach[hc10B[n10B[f]][4]] = ach_shell[hc10B[n10B[f]][4]] = 'B';
            if(ach[hc10B[n10B[f]][5]] == 'C') ach[hc10B[n10B[f]][5]] = ach_shell[hc10B[n10B[f]][5]] = 'B';
            ach[hc10B[n10B[f]][6]] = ach_shell[hc10B[n10B[f]][6]] = 'O';
            ach[hc10B[n10B[f]][7]] = ach_shell[hc10B[n10B[f]][7]] = 'O';
            ach[hc10B[n10B[f]][8]] = ach_shell[hc10B[n10B[f]][8]] = 'O';
            ach[hc10B[n10B[f]][9]] = ach_cen[hc10B[n10B[f]][9]] = 'O';

            if (do11W==1) Clusters_Get11W_Cs(f,ach1,ach1_cen,ach1_shell);
            ++n10B[f];
        }

        if (sp5c[k][6] == hc9B[n9B[f]][8]) {    // check one spindle of 7A_k is the common spindle of 9B (hc9B[id9B][.] at this point)
            if (Bonds_BondCheck(sp5c[k][5], hc9B[n9B[f]][6])==0) continue;  // check other spindle of 7A_k is bonded to spindle d1 of 9B
            if (Bonds_BondCheck(sp5c[k][5], hc9B[n9B[f]][7])==0) continue;  // check other spindle of 7A_k is bonded to spindle d2 of 9B

            flg1=0;
            flg2=0;
            for (l=0;l<5;l++) {
                if (sp5c[k][l]==hc9B[n9B[f]][6]) {
                    flg1=1;
                    continue;
                }
                if (sp5c[k][l]==hc9B[n9B[f]][7]) {
                    flg2=1;
                    continue;
                }
            }
            if (flg1==0 || flg2==0) continue;
            trial[6]=hc9B[n9B[f]][6];
            trial[7]=hc9B[n9B[f]][7];
            trial[8]=sp5c[k][5];
            trial[9]=hc9B[n9B[f]][8];

            m=0;
            break_out=0;
            for (l=0;l<6;l++) {
                if (hc9B[n9B[f]][l]==sp5c[k][5]) continue;
                if (m==5) {
                    m++;
                    break_out=1;
                    break;
                }
                trial[m]=hc9B[n9B[f]][l];
                m++;
            }
            if (break_out==1 || m!=5) continue;

            break_out=0;
            for (l=0;l<5;l++) {
                if (sp5c[k][l]==hc9B[n9B[f]][6]) continue;
                if (sp5c[k][l]==hc9B[n9B[f]][7]) continue;
                for (m=0;m<5;m++) {
                    if (sp5c[k][l]==trial[m]) break;
                }
                if (m==5) {
                    trial[5]=sp5c[k][l];
                    break_out++;
                }
            }
            if (break_out!=1) continue;

            if (n10B[f]==m10B) { hc10B=resize_2D_int(hc10B,m10B,m10B+incrStatic,clusSize,-1);
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
            for (l=0;l<10;l++) hc10B[n10B[f]][l]=trial[l];

            // hc10B key: (ordered shell particles, s1, s2, s3 (ordered), s_com)
            if(ach[hc10B[n10B[f]][0]] == 'C') ach[hc10B[n10B[f]][0]] = ach_shell[hc10B[n10B[f]][0]] = 'B';
            if(ach[hc10B[n10B[f]][1]] == 'C') ach[hc10B[n10B[f]][1]] = ach_shell[hc10B[n10B[f]][1]] = 'B';
            if(ach[hc10B[n10B[f]][2]] == 'C') ach[hc10B[n10B[f]][2]] = ach_shell[hc10B[n10B[f]][2]] = 'B';
            if(ach[hc10B[n10B[f]][3]] == 'C') ach[hc10B[n10B[f]][3]] = ach_shell[hc10B[n10B[f]][3]] = 'B';
            if(ach[hc10B[n10B[f]][4]] == 'C') ach[hc10B[n10B[f]][4]] = ach_shell[hc10B[n10B[f]][4]] = 'B';
            if(ach[hc10B[n10B[f]][5]] == 'C') ach[hc10B[n10B[f]][5]] = ach_shell[hc10B[n10B[f]][5]] = 'B';
            ach[hc10B[n10B[f]][6]] = ach_shell[hc10B[n10B[f]][6]] = 'O';
            ach[hc10B[n10B[f]][7]] = ach_shell[hc10B[n10B[f]][7]] = 'O';
            ach[hc10B[n10B[f]][8]] = ach_shell[hc10B[n10B[f]][8]] = 'O';
            ach[hc10B[n10B[f]][9]] = ach_cen[hc10B[n10B[f]][9]] = 'O';

            if (do11W==1) Clusters_Get11W_Cs(f,ach1,ach1_cen,ach1_shell);
            ++n10B[f];
        }
    }
}

int Clusters_Get11B_C2v(int f, char *ach, char *ach_cen, char *ach_shell) { // Detect 11B C2v clusters
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

    if(cnb[hc9B[n9B[f]][8]]!= 10) return 0; // s_com has 10 bonds in total (all forming the shell)

    m = 0;
    break_out=0;
    for(k=0; k<10; ++k) {
        for(l=0; l<8; ++l) {
            if(bNums[hc9B[n9B[f]][8]][k] == hc9B[n9B[f]][l]) break;
        }
        if(l==8){
            if(m==2) {
                break_out=1;
                break;
            }
            ep[m++] = bNums[hc9B[n9B[f]][8]][k];    // two extra particles
        }
    }
    if(break_out==1 || m<2) return 0;

    if(Bonds_BondCheck(ep[0], ep[1])==0) return 0;  // extra particles must be bonded
    nb1 = nb2 = 0;
    for(k=0; k<8; ++k){
        if(Bonds_BondCheck(hc9B[n9B[f]][k], ep[0])){
            if(nb1 == 2) return 0;  // extra particle 1 bonded to 2 members of 9B shell
            b1[nb1]=hc9B[n9B[f]][k];
            nb1++;
        }
        if(Bonds_BondCheck(hc9B[n9B[f]][k], ep[1])){
            if(nb2 == 2) return 0;  // extra particle 2 bonded to 2 members of 9B shell
            b2[nb2]=hc9B[n9B[f]][k];
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

    if(n11B[f] == m11B) {
        hc11B=resize_2D_int(hc11B,m11B,m11B+incrStatic,clusSize,-1);
        m11B=m11B+incrStatic;
    }
    hc11B[n11B[f]][0] = hc9B[n9B[f]][0];
    hc11B[n11B[f]][1] = hc9B[n9B[f]][1];
    hc11B[n11B[f]][2] = hc9B[n9B[f]][2];
    hc11B[n11B[f]][3] = hc9B[n9B[f]][3];
    hc11B[n11B[f]][4] = hc9B[n9B[f]][4];
    hc11B[n11B[f]][5] = hc9B[n9B[f]][5];
    hc11B[n11B[f]][6] = hc9B[n9B[f]][6];
    hc11B[n11B[f]][7] = hc9B[n9B[f]][7];
    hc11B[n11B[f]][8] = hc9B[n9B[f]][8];
    if (b1[0]==1) {
        hc11B[n11B[f]][9] = ep[0];
        hc11B[n11B[f]][10] = ep[1];
    }
    else {
        hc11B[n11B[f]][9] = ep[1];
        hc11B[n11B[f]][10] = ep[0];
    }
    // hc11B key: (as 9B, ep_1_to_9B_0, ep_2_to_9B_1)
    if(ach[hc11B[n11B[f]][0]] == 'C') ach[hc11B[n11B[f]][0]] = ach_shell[hc11B[n11B[f]][0]] = 'B';
    if(ach[hc11B[n11B[f]][1]] == 'C') ach[hc11B[n11B[f]][1]] = ach_shell[hc11B[n11B[f]][1]] = 'B';
    if(ach[hc11B[n11B[f]][2]] == 'C') ach[hc11B[n11B[f]][2]] = ach_shell[hc11B[n11B[f]][2]] = 'B';
    if(ach[hc11B[n11B[f]][3]] == 'C') ach[hc11B[n11B[f]][3]] = ach_shell[hc11B[n11B[f]][3]] = 'B';
    if(ach[hc11B[n11B[f]][4]] == 'C') ach[hc11B[n11B[f]][4]] = ach_shell[hc11B[n11B[f]][4]] = 'B';
    if(ach[hc11B[n11B[f]][5]] == 'C') ach[hc11B[n11B[f]][5]] = ach_shell[hc11B[n11B[f]][5]] = 'B';
    if(ach[hc11B[n11B[f]][9]] == 'C') ach[hc11B[n11B[f]][9]] = ach_shell[hc11B[n11B[f]][9]] = 'B';
    if(ach[hc11B[n11B[f]][10]] == 'C') ach[hc11B[n11B[f]][10]] = ach_shell[hc11B[n11B[f]][10]] = 'B';
    ach[hc11B[n11B[f]][6]] = ach_shell[hc11B[n11B[f]][6]] = 'O';
    ach[hc11B[n11B[f]][7]] = ach_shell[hc11B[n11B[f]][7]] = 'O';
    ach[hc11B[n11B[f]][8]] = ach_cen[hc11B[n11B[f]][8]] = 'O';

    return 1;
}

int Clusters_Get11W_Cs(int f, char *ach, char *ach_cen, char *ach_shell) {  // Detect 11W C2s clusters
    //  11W is the ground state of the 11 Wahnstrom particles
    // Call from 10B, the central particle must have exactly 10 particles bonded to it.
    // One extra particle bonded to 10B central particle but not bonded to three shell spindles of 7A in 10B

    int k, l, m;
    int ep=-1;
    int break_out;
    int clusSize=11;

    if(cnb[hc10B[n10B[f]][9]]!= 10) return 0;   // s_com has 10 bonds in total (all forming the shell)

    m = 0;  // find extra particle
    break_out=0;
    for(k=0; k<10; ++k) {
        for(l=0; l<9; ++l) {
            if(bNums[hc10B[n10B[f]][9]][k] == hc10B[n10B[f]][l]) break;
        }
        if(l==9){
            if(m==1) {
                break_out=1;
                break;
            }
            ep= bNums[hc10B[n10B[f]][9]][k];    // two extra particles
            m++;
        }
    }

    if(break_out==1 || m<1) return 0;
    if(Bonds_BondCheck(ep, hc10B[n10B[f]][6])==1) return 0; // extra particles must not be bonded to three 7A spindles in shell of 10B
    if(Bonds_BondCheck(ep, hc10B[n10B[f]][7])==1) return 0; // extra particles must not be bonded to three 7A spindles in shell of 10B
    if(Bonds_BondCheck(ep, hc10B[n10B[f]][8])==1) return 0; // extra particles must not be bonded to three 7A spindles in shell of 10B
    // Found 11W add as hc11W key: (as 10B, ep)
    if(n11W[f] == m11W) {
        hc11W=resize_2D_int(hc11W,m11W,m11W+incrStatic,clusSize,-1);
        m11W=m11W+incrStatic;
    }
    hc11W[n11W[f]][0] = hc10B[n10B[f]][0];
    hc11W[n11W[f]][1] = hc10B[n10B[f]][1];
    hc11W[n11W[f]][2] = hc10B[n10B[f]][2];
    hc11W[n11W[f]][3] = hc10B[n10B[f]][3];
    hc11W[n11W[f]][4] = hc10B[n10B[f]][4];
    hc11W[n11W[f]][5] = hc10B[n10B[f]][5];
    hc11W[n11W[f]][6] = hc10B[n10B[f]][6];
    hc11W[n11W[f]][7] = hc10B[n10B[f]][7];
    hc11W[n11W[f]][8] = hc10B[n10B[f]][8];
    hc11W[n11W[f]][9] = hc10B[n10B[f]][9];
    hc11W[n11W[f]][10] = ep;

    if(ach[hc11W[n11W[f]][0]] == 'C') ach[hc11W[n11W[f]][0]] = ach_shell[hc11W[n11W[f]][0]] = 'B';
    if(ach[hc11W[n11W[f]][1]] == 'C') ach[hc11W[n11W[f]][1]] = ach_shell[hc11W[n11W[f]][1]] = 'B';
    if(ach[hc11W[n11W[f]][2]] == 'C') ach[hc11W[n11W[f]][2]] = ach_shell[hc11W[n11W[f]][2]] = 'B';
    if(ach[hc11W[n11W[f]][3]] == 'C') ach[hc11W[n11W[f]][3]] = ach_shell[hc11W[n11W[f]][3]] = 'B';
    if(ach[hc11W[n11W[f]][4]] == 'C') ach[hc11W[n11W[f]][4]] = ach_shell[hc11W[n11W[f]][4]] = 'B';
    if(ach[hc11W[n11W[f]][5]] == 'C') ach[hc11W[n11W[f]][5]] = ach_shell[hc11W[n11W[f]][5]] = 'B';
    if(ach[hc11W[n11W[f]][6]] == 'C') ach[hc11W[n11W[f]][6]] = ach_shell[hc11W[n11W[f]][6]] = 'B';
    if(ach[hc11W[n11W[f]][7]] == 'C') ach[hc11W[n11W[f]][7]] = ach_shell[hc11W[n11W[f]][7]] = 'B';
    if(ach[hc11W[n11W[f]][8]] == 'C') ach[hc11W[n11W[f]][8]] = ach_shell[hc11W[n11W[f]][8]] = 'B';
    if(ach[hc11W[n11W[f]][9]] == 'C') ach[hc11W[n11W[f]][9]] = ach_cen[hc11W[n11W[f]][9]] = 'B';
    ach[hc11W[n11W[f]][10]] = ach_shell[hc11W[n11W[f]][10]] = 'O';

    n11W[f]++;

    return 1;
}

void Clusters_Get11E_12D(int f, int i, int j, int sp1, int sp2i, int sp2j, char *ach1, char *ach2) {    // Returns number of 11Es for a single 9B
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

    for(k=i+1; k<nsp5c[f]; k++) { // loop over all 7A_k
        if(k == j) continue;
        if(sp5c[k][5] == sp2j && sp5c[k][6] != sp2i) { // one 7A_k spindle is sp2j, one is not sp2i
            if(Bonds_BondCheck(sp5c[k][6], sp1) && Bonds_BondCheck(sp5c[k][6], sp2i)) { // non sp2j 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2j;
                trial[2] = sp2i;
                trial[3] = sp5c[k][6];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (sp5c[k][l]==sp1) {  // one SP5 ring particle of 7A_i common to common spindle of 9B
                        flg1=1;
                        continue;
                    }
                    if (sp5c[k][l]==sp2i) { // one SP5 ring particle of 7A_i common to the other uncommon spindle of 9B
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {   // one SP5 ring particle of 7A_i common to the uncommon particles of the SP5 rings in the 7A constituting 9B
                        if (sp5c[k][l]==hc9B[n9B[f]][m]) {
                            flg3=1;
                            trial[6]=sp5c[k][l];
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
                    trial[n]=sp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B[f]][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B[f]][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B[f]][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue; // fetched final 4 particles from 9B not in 7A_i

                if(n11E[f] == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E[f]][i]=trial[i];

                Clust_Write_11E(f, ach1);

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2i, sp5c[k][6], ach2);

                ++n11E[f];
            }
        }
        if(sp5c[k][6] == sp2j && sp5c[k][5] != sp2i) { // one 7A_k spindle is sp2j, one is not sp2i
            if(Bonds_BondCheck(sp5c[k][5], sp1) && Bonds_BondCheck(sp5c[k][5], sp2i)) { // non sp2j 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2j;
                trial[2] = sp2i;
                trial[3] = sp5c[k][5];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (sp5c[k][l]==sp1) {
                        flg1=1;
                        continue;
                    }
                    if (sp5c[k][l]==sp2i) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (sp5c[k][l]==hc9B[n9B[f]][m]) {
                            flg3=1;
                            trial[6]=sp5c[k][l];
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
                    trial[n]=sp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B[f]][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B[f]][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B[f]][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E[f] == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E[f]][i]=trial[i];

                Clust_Write_11E(f, ach1);

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2i, sp5c[k][5], ach2);
                ++n11E[f];
            }
        }
    }
    for(k=j+1; k<nsp5c[f]; k++) {
        if(sp5c[k][5] == sp2i && sp5c[k][6] != sp2j) {  // one 7A_k spindle is sp2i, one is not sp2j
            if(Bonds_BondCheck(sp5c[k][6], sp1) && Bonds_BondCheck(sp5c[k][6], sp2j)) { // non sp2i 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2i;
                trial[2] = sp2j;
                trial[3] = sp5c[k][6];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (sp5c[k][l]==sp1) {
                        flg1=1;
                        continue;
                    }
                    if (sp5c[k][l]==sp2j) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (sp5c[k][l]==hc9B[n9B[f]][m]) {
                            flg3=1;
                            trial[6]=sp5c[k][l];
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
                    trial[n]=sp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B[f]][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B[f]][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B[f]][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E[f] == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E[f]][i]=trial[i];

                Clust_Write_11E(f, ach1);

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2j, sp5c[k][6], ach2);
                ++n11E[f];
            }
        }
        if(sp5c[k][6] == sp2i && sp5c[k][5] != sp2j) {  // one 7A_k spindle is sp2i, one is not sp2j
            if(Bonds_BondCheck(sp5c[k][5], sp1) && Bonds_BondCheck(sp5c[k][5], sp2j)) { // non sp2i 7A_k spindle is bonded to sp1 and sp2i
                trial[0] = sp1;
                trial[1] = sp2i;
                trial[2] = sp2j;
                trial[3] = sp5c[k][5];

                flg1=flg2=flg3=0;
                n=4;
                break_out=0;
                for (l=0; l<5; ++l) {
                    if (sp5c[k][l]==sp1) {
                        flg1=1;
                        continue;
                    }
                    if (sp5c[k][l]==sp2j) {
                        flg2=1;
                        continue;
                    }
                    break_out2=0;
                    for (m=0; m<4; ++m) {
                        if (sp5c[k][l]==hc9B[n9B[f]][m]) {
                            flg3=1;
                            trial[6]=sp5c[k][l];
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
                    trial[n]=sp5c[k][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || flg3==0 || n!=6 || break_out==1) continue;

                n=7;
                flg1=flg2=0;
                break_out=0;
                for (l=0; l<6; ++l) {
                    if (hc9B[n9B[f]][l]==trial[6]) {
                        flg1=1;
                        continue;
                    }
                    if (hc9B[n9B[f]][l]==trial[3]) {
                        flg2=1;
                        continue;
                    }
                    if (n==11) {
                        n++;
                        break_out=1;
                        break;
                    }
                    trial[n]=hc9B[n9B[f]][l];
                    n++;
                }
                if (flg1==0 || flg2==0 || n!=11 || break_out==1) continue;

                if(n11E[f] == m11E) {
                    hc11E=resize_2D_int(hc11E,m11E,m11E+incrStatic,clusSize,-1);
                    m11E=m11E+incrStatic;
                }

                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (i=0;i<11;i++) hc11E[n11E[f]][i]=trial[i];

                Clust_Write_11E(f, ach1);

                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2j, sp5c[k][5], ach2);
                ++n11E[f];
            }
        }
    }
}

void Clust_Write_11E(int f, char *ach1) {
    if(ach1[hc11E[n11E[f]][4]] == 'C') ach1[hc11E[n11E[f]][4]] = 'B';
    if(ach1[hc11E[n11E[f]][5]] == 'C') ach1[hc11E[n11E[f]][5]] = 'B';
    if(ach1[hc11E[n11E[f]][6]] == 'C') ach1[hc11E[n11E[f]][6]] = 'B';
    if(ach1[hc11E[n11E[f]][7]] == 'C') ach1[hc11E[n11E[f]][7]] = 'B';
    if(ach1[hc11E[n11E[f]][8]] == 'C') ach1[hc11E[n11E[f]][8]] = 'B';
    if(ach1[hc11E[n11E[f]][9]] == 'C') ach1[hc11E[n11E[f]][9]] = 'B';
    if(ach1[hc11E[n11E[f]][10]] == 'C') ach1[hc11E[n11E[f]][10]] = 'B';
    ach1[hc11E[n11E[f]][0]] = 'O';
    ach1[hc11E[n11E[f]][1]] = 'O';
    ach1[hc11E[n11E[f]][2]] = 'O';
    ach1[hc11E[n11E[f]][3]] = 'O';
}

int Clusters_Get12D_D2d(int f, int i, int j, int k, int sp1, int sp2, char *ach) {  // Return 1 if 12B is also 11E
    int l, m, n, o, p, q;
    int flg1, flg2;
    int break_out;
    int trial[12];
    int clusSize=12;

    if(k > j) m = k;
    else m = j;
    for (l=m+1; l<nsp5c[f]; l++) {
        if ((sp5c[l][5] == sp1 && sp5c[l][6] == sp2) || (sp5c[l][6] == sp1 && sp5c[l][5] == sp2)) {
            flg1=flg2=0;
            p=11;
            q=0;
            break_out=0;
            for (n=0; n<5; n++) {
                if (sp5c[l][n]==hc11E[n11E[f]][0]) {
                    flg1=1;
                    continue;
                }
                if (sp5c[l][n]==hc11E[n11E[f]][1]) {
                    flg2=1;
                    continue;
                }
                break_out=0;
                for (o=4; o<11; o++) {
                    if (sp5c[l][n]==hc11E[n11E[f]][o]) {
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
                trial[p]=sp5c[l][n];
                p++;
            }
            if (flg1==0 || flg2==0 || p!=12 || q!=2) continue;

            for (n=0; n<11; n++) trial[n]=hc11E[n11E[f]][n];
            // hc12D key: (d1_unc, d2_unc, d3_unc, d4_unc, d12_com, d13_com, d24_com, d34_com, s_d1, s_d2, s_d3, s_com)

            if(n12D[f] == m12D) {
                hc12D=resize_2D_int(hc12D,m12D,m12D+incrStatic,clusSize,-1);
                m12D=m12D+incrStatic;
            }
            quickSort(&trial[0],4);
            quickSort(&trial[4],8);

            for(m=0; m<12; ++m) hc12D[n12D[f]][m] = trial[m];

            if(ach[hc12D[n12D[f]][4]] == 'C') ach[hc12D[n12D[f]][4]] = 'B';
            if(ach[hc12D[n12D[f]][5]] == 'C') ach[hc12D[n12D[f]][5]] = 'B';
            if(ach[hc12D[n12D[f]][6]] == 'C') ach[hc12D[n12D[f]][6]] = 'B';
            if(ach[hc12D[n12D[f]][7]] == 'C') ach[hc12D[n12D[f]][7]] = 'B';
            if(ach[hc12D[n12D[f]][8]] == 'C') ach[hc12D[n12D[f]][8]] = 'B';
            if(ach[hc12D[n12D[f]][9]] == 'C') ach[hc12D[n12D[f]][9]] = 'B';
            if(ach[hc12D[n12D[f]][10]] == 'C') ach[hc12D[n12D[f]][10]] = 'B';
            if(ach[hc12D[n12D[f]][11]] == 'C') ach[hc12D[n12D[f]][11]] = 'B';
            ach[hc12D[n12D[f]][0]] = 'O';
            ach[hc12D[n12D[f]][1]] = 'O';
            ach[hc12D[n12D[f]][2]] = 'O';
            ach[hc12D[n12D[f]][3]] = 'O';
            return 1;
        }
    }
    return 0;
}

void Clusters_Get9K_10K(int f)  { // Detect 9K & 10K clusters
    // Made from 2 sp4c clusters with a common sp4c spindle
    // particl and two common SP4 ring particles
    char *ach_1, *ach_cen_1, *ach_shell_1;
    int i, j2, j, k, l, m;
    int cp[2], scom, sother[2];
    int trial[9];
    char errMsg[1000];
    int clusSize=9;

    cp[0]=cp[1]=scom=sother[0]=sother[1]=-1;

    ach_1=malloc(N*sizeof(char));   if (ach_1==NULL) { sprintf(errMsg,"Clusters_Get9K_10K(): ach_1[] malloc out of memory\n");  Error(errMsg); }
    ach_cen_1=malloc(N*sizeof(char));   if (ach_cen_1==NULL) { sprintf(errMsg,"Clusters_Get9K_10K): ach_cen_1[] malloc out of memory\n");   Error(errMsg); }
    ach_shell_1=malloc(N*sizeof(char)); if (ach_shell_1==NULL) { sprintf(errMsg,"Clusters_Get9K_10K(): ach_shell_1[] malloc out of memory\n");  Error(errMsg); }

    for (i=0; i<N; ++i) {
        ach_1[i] = 'C';
        ach_cen_1[i] = 'C';
        ach_shell_1[i] = 'C';
    }

    for(i=0; i<nsp4c[f]-1; ++i) {   // loop over all sp4c_i
        for (j2=4; j2<6; j2++) {    // loop over all spindles of sp4c_i
            for (j=0; j<nmem_sp4c[sp4c[i][j2]]; ++j) {
                if (mem_sp4c[sp4c[i][j2]][j]<=i) continue; // don't find again 9K twice

                m=0;        // sp4c_i and sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly one common spindle
                for(k=4; k<6; ++k) {
                    for(l=4; l<6; ++l) {
                        if(sp4c[i][k] == sp4c[mem_sp4c[sp4c[i][j2]][j]][l]) {
                            m++;
                            scom=sp4c[i][k];
                        }
                    }
                }
                if(m!=1) continue;

                if (sp4c[i][4]==scom) sother[0]=sp4c[i][5];
                else sother[0]=sp4c[i][4];
                if (sp4c[mem_sp4c[sp4c[i][j2]][j]][4]==scom) sother[1]=sp4c[mem_sp4c[sp4c[i][j2]][j]][5];
                else sother[1]=sp4c[mem_sp4c[sp4c[i][j2]][j]][4];

                m=0;        // check sother[0] is not in cluster sp4c_mem_sp4c[sp4c[i][j2]][j]
                for(k=0; k<6; ++k) {
                    if(sother[0] == sp4c[mem_sp4c[sp4c[i][j2]][j]][k]) {
                        m++;
                    }
                }
                if(m!=0) continue;

                m=0;        // check sother[1] is not in cluster sp4c_mem_sp4c[sp4c[i][j2]][j]
                for(k=0; k<6; ++k) {
                    if(sother[1] == sp4c[i][k]) {
                        m++;
                    }
                }
                if(m!=0) continue;

                m=0;        // SP4 ring from sp4c_i and SP4 ring from sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly two common particles
                for(k=0; k<4; ++k) {
                    for(l=0; l<4; ++l) {
                        if(sp4c[i][k] == sp4c[mem_sp4c[sp4c[i][j2]][j]][l]) {
                            if (m>=2) {
                                m=3;
                                break;
                            }
                            cp[m]=sp4c[i][k];
                            m++;
                        }
                    }
                    if (m>2) break;
                }
                if(m!=2) continue;

                if(n9K[f] == m9K) {
                    hc9K=resize_2D_int(hc9K,m9K,m9K+incrStatic,clusSize,-1);
                    m9K=m9K+incrStatic;
                }
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
                        if (sp4c[i][k]==cp[l]) break;
                    }
                    if (l==2) {
                        if (m>=4) {
                            m++;
                            break;
                        }
                        trial[m]=sp4c[i][k];
                        m++;
                    }
                }
                if(m!=4) { sprintf(errMsg,"Clusters_Get9K_10K(): n9K[%d] %d, too many or too few uncommon SP4 ring particles found in cluster sp4c[i=%d][k]\n",f,n9K[f],i); Error(errMsg); }

                for(k=0; k<4; ++k) { // find uncommon SP4 ring particles in sp4c[mem_sp4c[sp4c[i][j2]][j]][k]
                    for (l=0; l<2; l++) {
                        if (sp4c[mem_sp4c[sp4c[i][j2]][j]][k]==cp[l]) break;
                    }
                    if (l==2) {
                        if (m>=6) {
                            m++;
                            break;
                        }
                        trial[m]=sp4c[mem_sp4c[sp4c[i][j2]][j]][k];
                        m++;
                    }
                }
                if(m!=6) { sprintf(errMsg,"Clusters_Get9K_10K(): n9K[%d] %d, too many or too few uncommon SP4 ring particles found in cluster sp4c[mem_sp4c[sp4c[i=%d][j2=%d]][j=%d]][k]\n",f,n9K[f],i,j2,j);   Error(errMsg); }

                quickSort(&trial[2],4);
                for (k=0; k<9; k++) hc9K[n9K[f]][k]=trial[k];

                if(ach_1[hc9K[n9K[f]][0]]  == 'C') ach_1[hc9K[n9K[f]][0]] = ach_shell_1[hc9K[n9K[f]][0]] = 'B';
                if(ach_1[hc9K[n9K[f]][1]]  == 'C') ach_1[hc9K[n9K[f]][1]] = ach_shell_1[hc9K[n9K[f]][1]] = 'B';
                if(ach_1[hc9K[n9K[f]][2]]  == 'C') ach_1[hc9K[n9K[f]][2]] = ach_shell_1[hc9K[n9K[f]][2]] = 'B';
                if(ach_1[hc9K[n9K[f]][3]]  == 'C') ach_1[hc9K[n9K[f]][3]] = ach_shell_1[hc9K[n9K[f]][3]] = 'B';
                if(ach_1[hc9K[n9K[f]][4]]  == 'C') ach_1[hc9K[n9K[f]][4]] = ach_shell_1[hc9K[n9K[f]][4]] = 'B';
                if(ach_1[hc9K[n9K[f]][5]]  == 'C') ach_1[hc9K[n9K[f]][5]] = ach_shell_1[hc9K[n9K[f]][5]] = 'B';
                ach_1[hc9K[n9K[f]][6]] = ach_shell_1[hc9K[n9K[f]][6]] = 'O';
                ach_1[hc9K[n9K[f]][7]] = ach_shell_1[hc9K[n9K[f]][7]] = 'O';
                ach_1[hc9K[n9K[f]][8]] = ach_cen_1[hc9K[n9K[f]][8]] = 'O';

                ++n9K[f];
            }
        }
    }

    for (i=0; i<N; ++i) {
        s9K[i]=ach_1[i];
        s9K_cen[i]=ach_cen_1[i];
        s9K_shell[i]=ach_shell_1[i];

    }
    free(ach_1);
    free(ach_cen_1);
    free(ach_shell_1);
}

void Clusters_Get10K(int f) { // Detect 10K clusters
    // A 10K is a 9K with a SINGLE particle bonded to the common spindle of 9K.
    int bonded_to_spindle_id, num_extra_particles, extra_particle = 0;
    int id_9K, id_9K_common;

    for (id_9K=0; id_9K<n9K[f]; id_9K++) {
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
                Cluster_Write_10K(f, id_9K, extra_particle);
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

void Cluster_Write_10K(int f, int id_9k, int extra_particle) {
    // hc10K key: (common_SP4_1, common_SP4_2, other_SP4*4, other_spindle_1, other_spindle_2, scom, ep)

    int i;
    int clusSize=10;

    if(n10K[f] == m10K) {
        hc10K=resize_2D_int(hc10K,m10K,m10K+incrStatic,clusSize,-1);
        m10K=m10K+incrStatic;
    }

    for(i=0; i<9; i++) {
        hc10K[n10K[f]][i]=hc9K[id_9k][i];
    }
    hc10K[n10K[f]][9]=extra_particle;

    for(i=0; i<6; i++) {
        if (s10K[hc10K[n10K[f]][i]] == 'C') s10K[hc10K[n10K[f]][i]] = s10K_shell[hc10K[n10K[f]][i]] = 'B';
    }
    s10K[hc10K[n10K[f]][6]] = s10K_shell[hc10K[n10K[f]][6]] = 'O';
    s10K[hc10K[n10K[f]][7]] = s10K_shell[hc10K[n10K[f]][7]] = 'O';
    s10K[hc10K[n10K[f]][8]] = s10K_cen[hc10K[n10K[f]][8]] = 'O';
    s10K[hc10K[n10K[f]][9]] = s10K_shell[hc10K[n10K[f]][9]] = 'O';

    n10K[f]++;
}

void Clusters_Get10A_C3v(int f) { // Detect 10A D4d clusters
    char *ach;
    int i, j, j2, k, l, m;
    char errMsg[1000];
    int clusSize=10;
    int *used_sp4b;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get10A_C3v(): ach[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i<N; ++i) ach[i] = 'C';

    used_sp4b=malloc(nsp4b[f]*sizeof(int)); if (used_sp4b==NULL) { sprintf(errMsg,"Clusters_Get10A_C3v(): used_sp4b[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i<nsp4b[f]; ++i) used_sp4b[i] = 0;

    for (i=0; i<nsp4b[f]-1; ++i) {  // loop over all sp4b_i
        for (j2=0; j2<nsp4b[f]; ++j2) used_sp4b[j2] = 0;
        used_sp4b[i]=1;
        for (j2=0; j2<cnb[sp4b[i][0]]; ++j2) {
            for (j=0; j<nmem_sp4b[bNums[sp4b[i][0]][j2]]; ++j) {    // loop over sp4b_j
                if (mem_sp4b[bNums[sp4b[i][0]][j2]][j]<=i) continue;
                if (used_sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]]==1) continue;
                used_sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]]=1;

                if(sp4b[i][4] == sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][4]) continue;
                if(Bonds_BondCheck(sp4b[i][4], sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][4])) continue;
                for(k=0; k<5; ++k){
                    for(l=0; l<5; ++l){
                        if(sp4b[i][k] == sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][l]) break;
                    }
                    if(l<5) break;
                }
                if(k<5) continue;
                for(k=0; k<4; ++k){
                    m = 0;
                    for(l=0; l<4; ++l) if(Bonds_BondCheck(sp4b[i][k], sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][l])) ++m;
                    if(m!=2) break;
                }
                // ERROR: Need to check converse, i. e. each SP4 ring particle from sp4b_j bonded to to exactly two particles from sp4b_i
                if(k==4) {
                    if(n10A[f] == m10A) {
                        hc10A=resize_2D_int(hc10A,m10A,m10A+incrStatic,clusSize,-1);
                        m10A=m10A+incrStatic;
                    }

                    // hc10A key: (SP4s going up, spindles going up)

                    hc10A[n10A[f]][0]=sp4b[i][0];
                    hc10A[n10A[f]][1]=sp4b[i][1];
                    hc10A[n10A[f]][2]=sp4b[i][2];
                    hc10A[n10A[f]][3]=sp4b[i][3];
                    hc10A[n10A[f]][4]=sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][0];
                    hc10A[n10A[f]][5]=sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][1];
                    hc10A[n10A[f]][6]=sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][2];
                    hc10A[n10A[f]][7]=sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][3];
                    quickSort(&hc10A[n10A[f]][0],8);
                    hc10A[n10A[f]][8]=sp4b[i][4];
                    hc10A[n10A[f]][9]=sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]][4];
                    quickSort(&hc10A[n10A[f]][8],2);

                    Cluster_Write_10A(f, ach);
                }
            }
        }
    }

    for (i=0; i<N; ++i) s10A[i]=ach[i];
    free(ach);
    free(used_sp4b);
}

void Cluster_Write_10A(int f, char *ach) {
    int i;

    for(i=0; i<8; i++) {
        if (ach[hc10A[n10A[f]][i]] == 'C') ach[hc10A[n10A[f]][i]] = 'B';
    }
    ach[hc10A[n10A[f]][8]] = 'O';
    ach[hc10A[n10A[f]][9]] = 'O';

    ++n10A[f];
}

void Clusters_Get10W(int f) { // Detect 10W clusters
    int i, j, k, l, m;
    int sp5b_clusts[5], shell_parts[9];
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int clusSize=10;

    sp5b_clusts[0]=sp5b_clusts[1]=sp5b_clusts[2]=sp5b_clusts[3]=sp5b_clusts[4]=-1;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get10W(): ach[] malloc out of memory\n"); Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_Get10W(): ach_cen[] malloc out of memory\n"); Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_Get10W(): ach_shell[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for (i=0; i<nsp5b[f]; ++i) { // loop over all sp5b
        if (cnb[sp5b[i][5]]!=9) continue;   // central particle must have coordination number 9

        k=0;    // find 5 other sp5b's with spindle in common with sp5b_i
        for (j=0; j<nmem_sp5b[sp5b[i][5]]; ++j) { // note check that spindle of sp5b_i and sp5b_j must be common by later check
            if (mem_sp5b[sp5b[i][5]][j]<=i) continue;   // i for sp5b must be lowest of all sp5b indices
            // ERROR !! need to check that spindle of sp5b_j is spindle of sp5b_i
            if (k>=5) {
                k++;
                break;
            }
            sp5b_clusts[k]=mem_sp5b[sp5b[i][5]][j];
            k++;
        }
        if (k!=5) continue; // not correct number of sp5b clusters
        // now found exactly 5 sp5b clusters common to spindle of sp5b_i
        for (j=0; j<5; j++) {
            shell_parts[j]=sp5b[i][j];
        }

        m=5;
        for (j=0; j<5; j++) {
            for (k=0; k<5; k++) {
                for (l=0; l<m; l++) {
                    if (shell_parts[l]==sp5b[mem_sp5b[sp5b[i][5]][j]][k]) break;
                }
                if (l==m) {
                    if (m>=9) {
                        m++;
                        break;
                    }
                    shell_parts[m]=sp5b[mem_sp5b[sp5b[i][5]][j]][k];
                    m++;
                }
            }
            if (m>=10) break;
        }
        if (m!=9) continue; // not all coordination shell particles of sp5b[i][5] are in the SP5 rings of the 5xsp5b clusters we found

        if (n10W[f] == m10W) {
            hc10W=resize_2D_int(hc10W,m10W,m10W+incrStatic,clusSize,-1);
            m10W=m10W+incrStatic;
        }
        // hc10W key: (sp5bs_common_central_spindle_particle, sp5bs_SP5_ring_shell_particles)
        hc10W[n10W[f]][0] = sp5b[i][5];
        for (j=0; j<9; j++) hc10W[n10W[f]][j+1]=shell_parts[j];
        quickSort(&hc10W[n10W[f]][1],9);
        Cluster_Write_10W(f, ach, ach_cen, ach_shell);
    }

    for(i=0; i<N; ++i) {
        s10W[i]=ach[i];
        s10W_cen[i]=ach_cen[i];
        s10W_shell[i]=ach_shell[i];
    }
    free(ach);
    free(ach_cen);
    free(ach_shell);
}

void Cluster_Write_10W(int f, char *ach, char *ach_cen, char *ach_shell) {
    int i;

    for(i=1; i<10; i++) {
        if (ach[hc10W[n10W[f]][i]] == 'C') ach[hc10W[n10W[f]][i]] = ach_shell[hc10W[n10W[f]][i]] = 'B';
    }
    ach[hc10W[n10W[f]][0]] = ach_cen[hc10W[n10W[f]][0]] = 'O';

    ++n10W[f];
}

void Clusters_Get11A(int f) {
    // Detect 11A clusters

    int first_6A_id, second_6A_id, k, l, m;
    int scom, sother[2];

    scom=sother[0]=sother[1]=-1;

    for(first_6A_id=0; first_6A_id<nsp4c[f]-1; ++first_6A_id){
        // POSSIBLE IMPROVEMENT: loop over all sp4c clusters for spindles of sp4c_i
        for(second_6A_id=first_6A_id+1; second_6A_id<nsp4c[f]; ++second_6A_id) {
            m=0;
            for (k=4; k<6; k++) {
                for (l=4; l<6; l++) {
                    if (sp4c[first_6A_id][k] == sp4c[second_6A_id][l]) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        scom=sp4c[first_6A_id][k];
                        m++;
                    }
                }
                if (m>=2) break;
            }
            if(m==1) {  // one common spindle
                // Identify the non-common spindles
                if (scom == sp4c[first_6A_id][4]) {
                    sother[0] = sp4c[first_6A_id][5];
                } else {
                    sother[0] = sp4c[first_6A_id][4];
                }
                if (scom == sp4c[second_6A_id][4]) {
                    sother[1] = sp4c[second_6A_id][5];
                } else {
                    sother[1] = sp4c[second_6A_id][4];
                }

                if (Check_unique_6A_rings(first_6A_id, second_6A_id) == 0) {
                    if (Check_6A_rings_bonded(first_6A_id, second_6A_id) == 1) {
                        Cluster_Write_11A(f, first_6A_id, second_6A_id, sother, scom);
                    }
                }
            }
        }
    }
}

void Cluster_Write_11A(int f, int first_6A_id, int second_6A_id, const int sother[], int scom) {
    int clusSize=11;
    int i;

    if (n11A[f] == m11A) {
        hc11A = resize_2D_int(hc11A, m11A, m11A + incrStatic, clusSize, -1);
        m11A = m11A + incrStatic;
    }

    hc11A[n11A[f]][0] = sp4c[first_6A_id][0];
    hc11A[n11A[f]][1] = sp4c[first_6A_id][1];
    hc11A[n11A[f]][2] = sp4c[first_6A_id][2];
    hc11A[n11A[f]][3] = sp4c[first_6A_id][3];
    hc11A[n11A[f]][4] = sp4c[second_6A_id][0];
    hc11A[n11A[f]][5] = sp4c[second_6A_id][1];
    hc11A[n11A[f]][6] = sp4c[second_6A_id][2];
    hc11A[n11A[f]][7] = sp4c[second_6A_id][3];
    hc11A[n11A[f]][8] = sother[0];
    hc11A[n11A[f]][9] = sother[1];
    hc11A[n11A[f]][10] = scom;

    for(i=0; i<8; i++) {
        if (s11A[hc11A[n11A[f]][i]] == 'C') s11A[hc11A[n11A[f]][i]] = s11A_shell[hc11A[n11A[f]][i]] = 'B';
    }
    for(i=8; i<11; i++) {
        s11A[hc11A[n11A[f]][i]] = s11A_shell[hc11A[n11A[f]][i]] = 'O';
    }
    ++n11A[f];
}

int Check_6A_rings_bonded(int first_6A_id, int second_6A_id) {
    int i, j, num_bonds;
    // Check if there are two bonds between each particle in ring 1 and particles in ring 2
    // Returns 1 if all ring 1 particles have 2 bonds to ring 2 particles, return 0 if not
    // i loops over the particles in the first ring, j loops over the particles in the second ring
    for(i=0; i<4; ++i) { // loop through all first ring particles
        num_bonds = 0;
        for(j=0; j<4; ++j) { // loop through second ring particles
            if(Bonds_BondCheck(sp4c[first_6A_id][i], sp4c[second_6A_id][j])) {
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
                    if(sp4c[first_6A_id][first_ring] == sp4c[second_6A_id][second_ring]) {
                        return 1;
                    }
                }
            }
    return 0;
}

void Clusters_Get12K(int f) {
    // 12K clusters are an 11A with an extra particle bonded to three of the ring particles of the 11A
    // Since there are 8 possible places to attach an extra particle, there may be more than 1 12K for each 11A.

    int ptr_11A, ring_number;
    int sp3_rings[8][3] = {-1};
    int *sp3_ring;

    // Loop over all 11A clusters
    for(ptr_11A=0; ptr_11A<n11A[f]; ptr_11A++) {

        // Find the particle IDs of the 8 sp3 rings made from the rings of 11A
        get_12K_ring_bonds(ptr_11A, sp3_rings);

        // Loop through the 8 sp3 rings made by the sp4s of 11A and see if a particle is attached to any of them
        for(ring_number=0; ring_number<8; ring_number++) {
            sp3_ring = sp3_rings[ring_number];
            find_12K_cluster(f, ptr_11A, sp3_ring);
        }
    }
}

void find_12K_cluster(int f, int ptr_11A, const int *sp3_ring) {
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
        Cluster_Write_12K(f, ep, ptr_11A);
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

void Cluster_Write_12K(int f, int ep, int id_11A) {

    int i;
    int clusSize=12;

    if(n12K[f] == m12K) {
        hc12K=resize_2D_int(hc12K,m12K,m12K+incrStatic,clusSize,-1);
        m12K=m12K+incrStatic;
    }
    // hc12K key: (SP4 going up, sd going up, scom, ep)

    for (i=0; i<11; i++) hc12K[n12K[f]][i] = hc11A[id_11A][i];
    hc12K[n12K[f]][11]=ep;

    if(s12K[hc12K[n12K[f]][0]]  == 'C') s12K[hc12K[n12K[f]][0]] = s12K_shell[hc12K[n12K[f]][0]] = 'B';
    if(s12K[hc12K[n12K[f]][1]]  == 'C') s12K[hc12K[n12K[f]][1]] = s12K_shell[hc12K[n12K[f]][1]] = 'B';
    if(s12K[hc12K[n12K[f]][2]]  == 'C') s12K[hc12K[n12K[f]][2]] = s12K_shell[hc12K[n12K[f]][2]] = 'B';
    if(s12K[hc12K[n12K[f]][3]]  == 'C') s12K[hc12K[n12K[f]][3]] = s12K_shell[hc12K[n12K[f]][3]] = 'B';
    if(s12K[hc12K[n12K[f]][4]]  == 'C') s12K[hc12K[n12K[f]][4]] = s12K_shell[hc12K[n12K[f]][4]] = 'B';
    if(s12K[hc12K[n12K[f]][5]]  == 'C') s12K[hc12K[n12K[f]][5]] = s12K_shell[hc12K[n12K[f]][5]] = 'B';
    if(s12K[hc12K[n12K[f]][6]]  == 'C') s12K[hc12K[n12K[f]][6]] = s12K_shell[hc12K[n12K[f]][6]] = 'B';
    if(s12K[hc12K[n12K[f]][7]]  == 'C') s12K[hc12K[n12K[f]][7]] = s12K_shell[hc12K[n12K[f]][7]] = 'B';
    s12K[hc12K[n12K[f]][8]] = s12K_shell[hc12K[n12K[f]][8]] = 'O';
    s12K[hc12K[n12K[f]][9]] = s12K_shell[hc12K[n12K[f]][9]] = 'O';
    s12K[hc12K[n12K[f]][10]] = s12K_cen[hc12K[n12K[f]][10]] = 'O';
    s12K[hc12K[n12K[f]][11]] = s12K_shell[hc12K[n12K[f]][11]] = 'O';
    n12K[f]++;
}

void Clusters_Get11C_12A(int f) { // Detect 11C Cs & 12A C2v clusters
    char *ach1, *ach1_cen, *ach1_shell, *ach2, *ach2_cen, *ach2_shell;
    int ar[2],sd[2];
    int i, j, k, l, m, ncom, spc;
    int flg;
    int break_out;
    char errMsg[1000];
    int clusSize=11;

    ar[0]=ar[1]=sd[0]=sd[1]=spc=-1;
    ach1=malloc(N*sizeof(char));    if (ach1==NULL) { sprintf(errMsg,"Clusters_Get11C_12A(): ach1[] malloc out of memory\n");   Error(errMsg); }
    ach1_cen=malloc(N*sizeof(char));    if (ach1_cen==NULL) { sprintf(errMsg,"Clusters_Get11C_12A(): ach1_cen[] malloc out of memory\n");   Error(errMsg); }
    ach1_shell=malloc(N*sizeof(char));  if (ach1_shell==NULL) { sprintf(errMsg,"Clusters_Get11C_12A(): ach1_shell[] malloc out of memory\n");   Error(errMsg); }
    ach2=malloc(N*sizeof(char));    if (ach2==NULL) { sprintf(errMsg,"Clusters_Get11C_12A(): ach2[] malloc out of memory\n");   Error(errMsg); }
    ach2_cen=malloc(N*sizeof(char));    if (ach2_cen==NULL) { sprintf(errMsg,"Clusters_Get11C_12A(): ach2_cen[] malloc out of memory\n");   Error(errMsg); }
    ach2_shell=malloc(N*sizeof(char));  if (ach2_shell==NULL) { sprintf(errMsg,"Clusters_Get11C_12A(): ach2_shell[] malloc out of memory\n");   Error(errMsg); }
    for(i=0; i<N; ++i) ach1[i]=ach1_cen[i]=ach1_shell[i]=ach2[i]=ach2_cen[i]=ach2_shell[i]='C';

    for (i=0; i<nsp5c[f]-1; ++i) {
        // POSSIBLE IMPROVEMENT: loop over all spindles of 7A_i
        for (j=i+1; j<nsp5c[f]; ++j) {
            ncom = 0;
            if (sp5c[i][5] == sp5c[j][5]) {
                spc = sp5c[i][5];
                sd[0] = sp5c[i][6];
                sd[1] = sp5c[j][6];
                ++ncom;
            }
            if (sp5c[i][6] == sp5c[j][6]) {
                spc = sp5c[i][6];
                sd[0] = sp5c[i][5];
                sd[1] = sp5c[j][5];
                ++ncom;
            }
            if (sp5c[i][5] == sp5c[j][6]) {
                spc = sp5c[i][5];
                sd[0] = sp5c[i][6];
                sd[1] = sp5c[j][5];
                ++ncom;
            }
            if (sp5c[i][6] == sp5c[j][5]) {
                spc = sp5c[i][6];
                sd[0] = sp5c[i][5];
                sd[1] = sp5c[j][6];
                ++ncom;
            }
            if (ncom == 1) { // One common spindle particle
                ncom = 0;
                // need two common particles from SP5 rings
                break_out=0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if (sp5c[i][k] == sp5c[j][l]) {
                            if (ncom >= 2) {
                                ++ncom;
                                break_out=1;
                                break;
                            }
                            ar[ncom++] = sp5c[i][k];
                            break;
                        }
                    }
                    if (break_out==1) break;
                }
                flg = ncom == 2 && Bonds_BondCheck(ar[0],ar[1]) && break_out==0; // two common SP5 ring particles are bonded
                if(flg==1) {
                    ncom = 0;
                    for(k=0; k<5; ++k) {
                        if(sp5c[i][k] == ar[0] || sp5c[i][k] == ar[1]) continue;
                        for(l=0; l<5; ++l) {
                            if(sp5c[j][l] == ar[0] || sp5c[j][l] == ar[1]) continue;
                            if(Bonds_BondCheck(sp5c[i][k], sp5c[j][l])) ++ncom;
                        }
                    }
                    if(ncom != 2) flg = 0;
                }
                if(flg==1) { // two bonds between non-common SP5 ring particles
                    if(n11C[f] == m11C) {
                        hc11C=resize_2D_int(hc11C,m11C,m11C+incrStatic,clusSize,-1);
                        m11C=m11C+incrStatic;
                    }

                    // hc11C key: (s_com, s_i, s_j, r_ca, r_cb, d_i, d_i, d_j, d_j, unc_i, unc_j)

                    hc11C[n11C[f]][0]=spc;
                    hc11C[n11C[f]][1]=sd[0];
                    hc11C[n11C[f]][2]=sd[1];
                    hc11C[n11C[f]][3]=ar[0];
                    hc11C[n11C[f]][4]=ar[1];

                    l=5;
                    m=7;
                    break_out=0;
                    for(k=0; k<5; ++k) {
                        if(Bonds_BondCheck(sp5c[i][k], ar[0]) && sp5c[i][k]!=ar[0] && sp5c[i][k]!=ar[1]) {
                            if (l==7) {
                                break_out=1;
                                break;
                            }
                            hc11C[n11C[f]][l] = sp5c[i][k];
                            l++;
                        }
                        if(Bonds_BondCheck(sp5c[i][k], ar[1]) && sp5c[i][k]!=ar[0] && sp5c[i][k]!=ar[1]) {
                            if (l==7) {
                                break_out=1;
                                break;
                            }
                            hc11C[n11C[f]][l] = sp5c[i][k];
                            l++;
                        }
                        if(Bonds_BondCheck(sp5c[j][k], ar[0]) && sp5c[j][k]!=ar[0] && sp5c[j][k]!=ar[1]) {
                            if (m==9) {
                                break_out=1;
                                break;
                            }
                            hc11C[n11C[f]][m] = sp5c[j][k];
                            m++;
                        }
                        if(Bonds_BondCheck(sp5c[j][k], ar[1]) && sp5c[j][k]!=ar[0] && sp5c[j][k]!=ar[1]) {
                            if (m==9) {
                                break_out=1;
                                break;
                            }
                            hc11C[n11C[f]][m] = sp5c[j][k];
                            m++;
                        }
                    }
                    if (break_out==1 || l<7 || m<9) continue;

                    for(k=0; k<5; ++k) {
                        if(Bonds_BondCheck(sp5c[i][k], hc11C[n11C[f]][5]) && Bonds_BondCheck(sp5c[i][k], hc11C[n11C[f]][6])) {
                            hc11C[n11C[f]][9]=sp5c[i][k];
                        }
                        if(Bonds_BondCheck(sp5c[j][k], hc11C[n11C[f]][7]) && Bonds_BondCheck(sp5c[j][k], hc11C[n11C[f]][8])) {
                            hc11C[n11C[f]][10]=sp5c[j][k];
                        }
                    }
                    quickSort(&hc11C[n11C[f]][1],2);
                    quickSort(&hc11C[n11C[f]][3],2);
                    quickSort(&hc11C[n11C[f]][5],4);
                    quickSort(&hc11C[n11C[f]][9],2);

                    if(ach1[hc11C[n11C[f]][3]] == 'C') ach1[hc11C[n11C[f]][3]] = ach1_shell[hc11C[n11C[f]][3]] = 'B';
                    if(ach1[hc11C[n11C[f]][4]] == 'C') ach1[hc11C[n11C[f]][4]] = ach1_shell[hc11C[n11C[f]][4]] = 'B';
                    if(ach1[hc11C[n11C[f]][5]] == 'C') ach1[hc11C[n11C[f]][5]] = ach1_shell[hc11C[n11C[f]][5]] = 'B';
                    if(ach1[hc11C[n11C[f]][6]] == 'C') ach1[hc11C[n11C[f]][6]] = ach1_shell[hc11C[n11C[f]][6]] = 'B';
                    if(ach1[hc11C[n11C[f]][7]] == 'C') ach1[hc11C[n11C[f]][7]] = ach1_shell[hc11C[n11C[f]][7]] = 'B';
                    if(ach1[hc11C[n11C[f]][8]] == 'C') ach1[hc11C[n11C[f]][8]] = ach1_shell[hc11C[n11C[f]][8]] = 'B';
                    if(ach1[hc11C[n11C[f]][9]] == 'C') ach1[hc11C[n11C[f]][9]] = ach1_shell[hc11C[n11C[f]][9]] = 'B';
                    if(ach1[hc11C[n11C[f]][10]] == 'C') ach1[hc11C[n11C[f]][10]] = ach1_shell[hc11C[n11C[f]][10]] = 'B';
                    ach1[hc11C[n11C[f]][0]] = ach1_cen[hc11C[n11C[f]][0]] = 'O';
                    ach1[hc11C[n11C[f]][1]] = ach1_shell[hc11C[n11C[f]][1]] = 'O';
                    ach1[hc11C[n11C[f]][2]] = ach1_shell[hc11C[n11C[f]][2]] = 'O';

                    if (do12A==1) {
                        if(Clusters_Get12A_C2v(f, ach2, ach2_cen, ach2_shell)) {
                            ach2_cen[hc11C[n11C[f]][0]] = 'O';
                            ++n12A[f];
                        }
                    }
                    ++n11C[f];
                }
            }
        }
    }

    for(i=0; i<N; ++i) {
        s11C[i]=ach1[i];
        s11C_cen[i]=ach1_cen[i];
        s11C_shell[i]=ach1_shell[i];
        s12A[i]=ach2[i];
        s12A_cen[i]=ach2_cen[i];
        s12A_shell[i]=ach2_shell[i];
    }
    free(ach1);
    free(ach1_cen);
    free(ach1_shell);
    free(ach2);
    free(ach2_cen);
    free(ach2_shell);
}

int Clusters_Get12A_C2v(int f, char *ach, char *ach_cen, char *ach_shell) { // Return 1 if 11C C2v is also 12A C2v
    //  the central particle must have 11 particles bonded to it. The 11th
    // particle is only bonded to 2 other outer shell particles.
    int k, l;
    int ep;
    int clusSize=12;

    ep=-1;

    if(cnb[hc11C[n11C[f]][0]] != 11) return 0;

    for(k=0; k<11; ++k){
        for(l=1; l<11; ++l) if(bNums[hc11C[n11C[f]][0]][k] == hc11C[n11C[f]][l]) break;
        if(l == 11){
            ep = bNums[hc11C[n11C[f]][0]][k]; // The extra particle
            break;
        }
    }

    if (Bonds_BondCheck(ep,hc11C[n11C[f]][9])==0 || Bonds_BondCheck(ep,hc11C[n11C[f]][10])==0) return 0;

    for(k=2; k<9; ++k){
        if(Bonds_BondCheck(ep,hc11C[n11C[f]][k])) return 0;
    }

    if(n12A[f] == m12A) {
        hc12A=resize_2D_int(hc12A,m12A,m12A+incrStatic,clusSize,-1);
        m12A=m12A+incrStatic;
    }

    // hc12A key: (as 11C, extra_s)

    for(k=0;k<11;k++) hc12A[n12A[f]][k]=hc11C[n11C[f]][k];
    hc12A[n12A[f]][11]=ep;

    if(ach[hc12A[n12A[f]][3]] == 'C') ach[hc12A[n12A[f]][3]] = ach_shell[hc12A[n12A[f]][3]] = 'B';
    if(ach[hc12A[n12A[f]][4]] == 'C') ach[hc12A[n12A[f]][4]] = ach_shell[hc12A[n12A[f]][4]] = 'B';
    if(ach[hc12A[n12A[f]][5]] == 'C') ach[hc12A[n12A[f]][5]] = ach_shell[hc12A[n12A[f]][5]] = 'B';
    if(ach[hc12A[n12A[f]][6]] == 'C') ach[hc12A[n12A[f]][6]] = ach_shell[hc12A[n12A[f]][6]] = 'B';
    if(ach[hc12A[n12A[f]][7]] == 'C') ach[hc12A[n12A[f]][7]] = ach_shell[hc12A[n12A[f]][7]] = 'B';
    if(ach[hc12A[n12A[f]][8]] == 'C') ach[hc12A[n12A[f]][8]] = ach_shell[hc12A[n12A[f]][8]] = 'B';
    if(ach[hc12A[n12A[f]][9]] == 'C') ach[hc12A[n12A[f]][9]] = ach_shell[hc12A[n12A[f]][9]] = 'B';
    if(ach[hc12A[n12A[f]][10]] == 'C') ach[hc12A[n12A[f]][10]] = 'B';
    ach[hc12A[n12A[f]][0]] = ach_cen[hc12A[n12A[f]][0]] = 'O';
    ach[hc12A[n12A[f]][1]] = ach_shell[hc12A[n12A[f]][1]] = 'O';
    ach[hc12A[n12A[f]][2]] = ach_shell[hc12A[n12A[f]][2]] = 'O';
    ach[hc12A[n12A[f]][11]] = ach_shell[hc12A[n12A[f]][11]] = 'O';

    return 1;
}

void Clusters_Get11F_12E_13K(int f) {   // Detect 11F C2v & 12E 3h
    char *ach1, *ach2, *ach3, *ach3_cen, *ach3_shell;
    int cp, bpi, bpj, ep1, ep2, the6A_i, the6A_j;
    int i, j, j2, k, l, m;
    int flg, flg1, flg2;
    int break_out;
    char errMsg[1000];
    int clusSize=11;

    cp=bpi=bpj=ep1=ep2=the6A_i=the6A_j=-1;
    ach1=malloc(N*sizeof(char));    if (ach1==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach1[] malloc out of memory\n");   Error(errMsg); }
    ach2=malloc(N*sizeof(char));    if (ach2==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach2[] malloc out of memory\n");   Error(errMsg); }
    ach3=malloc(N*sizeof(char));    if (ach3==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach3[] malloc out of memory\n");   Error(errMsg); }
    ach3_cen=malloc(N*sizeof(char));    if (ach3_cen==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach3_cen[] malloc out of memory\n");   Error(errMsg); }
    ach3_shell=malloc(N*sizeof(char));  if (ach3_shell==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach3_shell[] malloc out of memory\n");   Error(errMsg); }
    for(i=0; i<N; ++i) ach1[i] = ach2[i] = ach3[i] = ach3_cen[i] = ach3_shell[i] ='C';

    for(i=0; i<nsp3c[f]-1; i++) {   // loop over all sp3c clusters
        for (j2=0; j2<3; j2++) {    // loop over only the rings of the sp3c clusters
            for (j=0; j<nmem_sp3c[sp3c[i][j2]]; ++j) {
                if (mem_sp3c[sp3c[i][j2]][j]>i) {
                    flg = sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4] ||
                          sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3] ||
                          sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4] ||
                          sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                    if (flg == 1) continue;
                    if (Bonds_BondCheck(sp3c[i][3], sp3c[mem_sp3c[sp3c[i][j2]][j]][3]) == 1 &&
                        Bonds_BondCheck(sp3c[i][4], sp3c[mem_sp3c[sp3c[i][j2]][j]][4]) == 1) {
                        m = 0;
                        for (k = 0; k < 3; ++k) {
                            for (l = 0; l < 3; ++l) {
                                if (sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                                    cp = sp3c[i][k];
                                    ++m;
                                }
                            }
                        }

                        if (m == 1) { // only 1 particle common between the pair
                            m = 0;
                            for (k = 0; k < 3; ++k) {
                                for (l = 0; l < 3; ++l) {
                                    if (sp3c[i][k] == cp || sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == cp) continue;
                                    if (Bonds_BondCheck(sp3c[i][k], sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) == 1) {
                                        if (m++) break;
                                        bpi = sp3c[i][k];
                                        bpj = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                                    }
                                }
                            }

                            if (m == 1) { // There must be exactly one paticle of sp3_i bonded to exactly one of sp3_j
                                flg1 = flg2 = 0;
                                for (k = 0; k < nsp4c[f]; ++k) {
                                    if (sp4c[k][4] == cp || sp4c[k][5] == cp) {
                                        if (flg1 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = sp4c[k][l] == sp3c[i][3] ||
                                                      sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                                                flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg1 = 1;
                                                if (sp4c[k][4] == cp) ep1 = sp4c[k][5];
                                                else ep1 = sp4c[k][4];
                                                the6A_i = k;
                                            }
                                        }
                                        if (flg2 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = sp4c[k][l] == sp3c[i][4] ||
                                                      sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                                                flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg2 = 1;
                                                if (sp4c[k][4] == cp) ep2 = sp4c[k][5];
                                                else ep2 = sp4c[k][4];
                                                the6A_j = k;
                                            }
                                        }
                                    }
                                    if (flg1 == 1 && flg2 == 1) break;
                                }

                                if (k < nsp4c[f]) { // 11F found
                                    if (n11F[f] == m11F) {
                                        hc11F = resize_2D_int(hc11F, m11F, m11F + incrStatic, clusSize, -1);
                                        m11F = m11F + incrStatic;
                                    }

                                    // hc11F key: (com, ep1(6A extra spindle), ep2(6A extra spindle), 5A_i_s1, 5A_i_s2, 5A_j_s1, 5A_j_s2, 4 ring particles)

                                    hc11F[n11F[f]][0] = cp;
                                    hc11F[n11F[f]][1] = ep1;
                                    hc11F[n11F[f]][2] = ep2;
                                    hc11F[n11F[f]][3] = sp3c[i][3];
                                    hc11F[n11F[f]][4] = sp3c[i][4];
                                    hc11F[n11F[f]][5] = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                                    hc11F[n11F[f]][6] = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];

                                    m = 7;
                                    break_out = 0;
                                    for (l = 0; l < 3; ++l) {
                                        if (sp3c[i][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F[f]][m] = sp3c[i][l];
                                            m++;
                                        }
                                    }
                                    for (l = 0; l < 3; ++l) {
                                        if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F[f]][m] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                                            m++;
                                        }
                                    }
                                    if (break_out == 1 || m < 11) continue;

                                    quickSort(&hc11F[n11F[f]][1], 2);
                                    quickSort(&hc11F[n11F[f]][3], 4);
                                    quickSort(&hc11F[n11F[f]][7], 4);

                                    if (ach1[hc11F[n11F[f]][0]] == 'C') ach1[hc11F[n11F[f]][0]] = 'B';
                                    if (ach1[hc11F[n11F[f]][7]] == 'C') ach1[hc11F[n11F[f]][7]] = 'B';
                                    if (ach1[hc11F[n11F[f]][8]] == 'C') ach1[hc11F[n11F[f]][8]] = 'B';
                                    if (ach1[hc11F[n11F[f]][9]] == 'C') ach1[hc11F[n11F[f]][9]] = 'B';
                                    if (ach1[hc11F[n11F[f]][10]] == 'C') ach1[hc11F[n11F[f]][10]] = 'B';
                                    ach1[hc11F[n11F[f]][1]] = 'O';
                                    ach1[hc11F[n11F[f]][2]] = 'O';
                                    ach1[hc11F[n11F[f]][3]] = 'O';
                                    ach1[hc11F[n11F[f]][4]] = 'O';
                                    ach1[hc11F[n11F[f]][5]] = 'O';
                                    ach1[hc11F[n11F[f]][6]] = 'O';

                                    if (do12E == 1) {
                                        if (Clusters_Get12E_D3h(f, mem_sp3c[sp3c[i][j2]][j], ach2)) ++n12E[f];
                                    }

                                    if (do13K == 1) {
                                        if (Clusters_Get13K(f, i, mem_sp3c[sp3c[i][j2]][j], the6A_i, ach3, ach3_cen,
                                                            ach3_shell))
                                            ++n13K[f];
                                    }
                                    ++n11F[f];
                                }
                            }
                        }
                    }

                    if (Bonds_BondCheck(sp3c[i][3], sp3c[mem_sp3c[sp3c[i][j2]][j]][4]) == 1 &&
                        Bonds_BondCheck(sp3c[i][4], sp3c[mem_sp3c[sp3c[i][j2]][j]][3]) == 1) {
                        m = 0;
                        for (k = 0; k < 3; ++k) {
                            for (l = 0; l < 3; ++l) {
                                if (sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                                    cp = sp3c[i][k];
                                    ++m;
                                }
                            }
                        }

                        if (m == 1) { // only 1 particle common between the pair
                            m = 0;
                            for (k = 0; k < 3; ++k) {
                                for (l = 0; l < 3; ++l) {
                                    if (sp3c[i][k] == cp || sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == cp) continue;
                                    if (Bonds_BondCheck(sp3c[i][k], sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) == 1) {
                                        if (m++) break;
                                        bpi = sp3c[i][k];
                                        bpj = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                                    }
                                }
                            }

                            if (m == 1) { // The single bonded particle found
                                flg1 = flg2 = 0;
                                for (k = 0; k < nsp4c[f]; ++k) {
                                    if (sp4c[k][4] == cp || sp4c[k][5] == cp) {
                                        if (flg1 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = sp4c[k][l] == sp3c[i][3] ||
                                                      sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                                                flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg1 = 1;
                                                if (sp4c[k][4] == cp) ep1 = sp4c[k][5];
                                                else ep1 = sp4c[k][4];
                                                the6A_i = k;
                                            }
                                        }
                                        if (flg2 == 0) { // check for first sp4c
                                            for (l = 0; l < 4; ++l) {
                                                flg = sp4c[k][l] == sp3c[i][4] ||
                                                      sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                                                flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                                if (flg == 0) break;
                                            }
                                            if (l == 4) {
                                                flg2 = 1;
                                                if (sp4c[k][4] == cp) ep2 = sp4c[k][5];
                                                else ep2 = sp4c[k][4];
                                                the6A_j = k;
                                            }
                                        }
                                    }
                                    if (flg1 == 1 && flg2 == 1) break;
                                }
                                if (k < nsp4c[f]) { // 11F found
                                    if (n11F[f] == m11F) {
                                        hc11F = resize_2D_int(hc11F, m11F, m11F + incrStatic, clusSize, -1);
                                        m11F = m11F + incrStatic;
                                    }

                                    // hc11F key: (com, ep1(6A extra spindle), ep2(6A extra spindle), 5A_i_s1, 5A_i_s2, 5A_j_s1, 5A_j_s2, 4 ring particles)

                                    hc11F[n11F[f]][0] = cp;
                                    hc11F[n11F[f]][1] = ep1;
                                    hc11F[n11F[f]][2] = ep2;
                                    hc11F[n11F[f]][3] = sp3c[i][3];
                                    hc11F[n11F[f]][4] = sp3c[i][4];
                                    hc11F[n11F[f]][5] = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                                    hc11F[n11F[f]][6] = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];

                                    m = 7;
                                    break_out = 0;
                                    for (l = 0; l < 3; ++l) {
                                        if (sp3c[i][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F[f]][m] = sp3c[i][l];
                                            m++;
                                        }
                                    }
                                    for (l = 0; l < 3; ++l) {
                                        if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l] != cp) {
                                            if (m == 11) {
                                                break_out = 1;
                                                break;
                                            }
                                            hc11F[n11F[f]][m] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                                            m++;
                                        }
                                    }
                                    if (break_out == 1 || m < 11) continue;

                                    quickSort(&hc11F[n11F[f]][1], 2);
                                    quickSort(&hc11F[n11F[f]][3], 4);
                                    quickSort(&hc11F[n11F[f]][7], 4);

                                    if (ach1[hc11F[n11F[f]][0]] == 'C') ach1[hc11F[n11F[f]][0]] = 'B';
                                    if (ach1[hc11F[n11F[f]][7]] == 'C') ach1[hc11F[n11F[f]][7]] = 'B';
                                    if (ach1[hc11F[n11F[f]][8]] == 'C') ach1[hc11F[n11F[f]][8]] = 'B';
                                    if (ach1[hc11F[n11F[f]][9]] == 'C') ach1[hc11F[n11F[f]][9]] = 'B';
                                    if (ach1[hc11F[n11F[f]][10]] == 'C') ach1[hc11F[n11F[f]][10]] = 'B';
                                    ach1[hc11F[n11F[f]][1]] = 'O';
                                    ach1[hc11F[n11F[f]][2]] = 'O';
                                    ach1[hc11F[n11F[f]][3]] = 'O';
                                    ach1[hc11F[n11F[f]][4]] = 'O';
                                    ach1[hc11F[n11F[f]][5]] = 'O';
                                    ach1[hc11F[n11F[f]][6]] = 'O';

                                    if (do12K == 1) {
                                        if (Clusters_Get12E_D3h(f, mem_sp3c[sp3c[i][j2]][j], ach2)) ++n12E[f];
                                    }

                                    if (do13K == 1) {
                                        if (Clusters_Get13K(f, i, mem_sp3c[sp3c[i][j2]][j], the6A_i, ach3, ach3_cen,
                                                            ach3_shell))
                                            ++n13K[f];
                                    }
                                    ++n11F[f];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for(i=0; i<N; ++i) s11F[i]=ach1[i];
    for(i=0; i<N; ++i) s12E[i]=ach2[i];
    for(i=0; i<N; ++i) s13K[i]=ach3[i];
    for(i=0; i<N; ++i) s13K_cen[i]=ach3_cen[i];
    for(i=0; i<N; ++i) s13K_shell[i]=ach3_shell[i];
    free(ach1);
    free(ach2);
    free(ach3);
    free(ach3_cen);
    free(ach3_shell);
}

int Clusters_Get12E_D3h(int f, int j, char *ach) {  // Return 1 is 11F is also 12E
    //  Made from three sp3c or 5A clusters
    int k, l, m, ncom, common[2], uncom;
    int flg;
    int clusSize=12;

    for(k=j+1; k<nsp3c[f]; ++k) {
        flg = (sp3c[k][3] == hc11F[n11F[f]][1] && sp3c[k][4] == hc11F[n11F[f]][2])  || (sp3c[k][3] == hc11F[n11F[f]][2] && sp3c[k][4] == hc11F[n11F[f]][1]);
        if (flg==1) { // spindles bonded correctly :)
            ncom=common[0]=common[1]=0;
            for (l=0; l<3; ++l) {
                for (m=0; m<5; ++m) {
                    if (m==0) uncom=0;
                    else uncom=m+6;
                    if(sp3c[k][l] == hc11F[n11F[f]][uncom]) {
                        if (ncom==2) {
                            ncom++;
                            break;
                        }
                        common[ncom]=sp3c[k][l];
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
                    if(sp3c[k][l] == common[m]) ncom++;
                }
                if (ncom==0) {
                    uncom=sp3c[k][l];
                    break;
                }
            }

            // now we have found the 12E
            if(n12E[f] == m12E) {
                hc12E=resize_2D_int(hc12E,m12E,m12E+incrStatic,clusSize,-1);
                m12E=m12E+incrStatic;
            }

            hc12E[n12E[f]][11] = uncom;
            for(l=0; l<10; ++l) hc12E[n12E[f]][l] = hc11F[n11F[f]][l+1];
            hc12E[n12E[f]][10] = hc11F[n11F[f]][0];
            quickSort(&hc12E[n12E[f]][0],6);
            quickSort(&hc12E[n12E[f]][6],6);

            if(ach[hc12E[n12E[f]][6]] == 'C') ach[hc12E[n12E[f]][6]] = 'B';
            if(ach[hc12E[n12E[f]][7]] == 'C') ach[hc12E[n12E[f]][7]] = 'B';
            if(ach[hc12E[n12E[f]][8]] == 'C') ach[hc12E[n12E[f]][8]] = 'B';
            if(ach[hc12E[n12E[f]][9]] == 'C') ach[hc12E[n12E[f]][9]] = 'B';
            if(ach[hc12E[n12E[f]][10]] == 'C') ach[hc12E[n12E[f]][10]] = 'B';
            if(ach[hc12E[n12E[f]][11]] == 'C') ach[hc12E[n12E[f]][11]] = 'B';
            ach[hc12E[n12E[f]][0]] = 'O';
            ach[hc12E[n12E[f]][1]] = 'O';
            ach[hc12E[n12E[f]][2]] = 'O';
            ach[hc12E[n12E[f]][3]] = 'O';
            ach[hc12E[n12E[f]][4]] = 'O';
            ach[hc12E[n12E[f]][5]] = 'O';
            return 1;
        }
    }
    return 0;
}

int Clusters_Get13K(int f, int sp3c_i, int sp3c_j, int the6A_i, char *ach, char *ach_cen, char *ach_shell) {
    /* Function Clusters_Get13K - Take an 11F particle and determine if it meets the criteria for the presence of a 13K
     *
     * f: Frame number currently being analysed
     * sp3c_i: The id of a relevant 5A cluster
     * sp3c_i: The id of a different relevant 5A cluster
     * the6A_i: The id of a relevant 6A ring
     * ach: A list keeping track of which particles are in a 13K
     * ach_cen: A list keeping track of which particles are in the center of a 13K
     * ach_shell: A list keeping track of which particles are in the shell of a 13K
     *
     * Returns 1 if a 13K is successfully detected
     * Returns 0 if no 13K is detected
     * 13K arrays are edited in place to add new 13K
     */
    int i, j, k, l;
    int sp3c_i_unc, sp3c_j_unc, ep[2], eclus5A[2], tmp;
    int clusSize=13;

    sp3c_i_unc=sp3c_j_unc=ep[0]=ep[1]=eclus5A[0]=eclus5A[1]=tmp=-1;

    // Clusters 5A_i and 5A_j are new to the 13K, they have:
    // two spindle particles in the 11F
    // one sp3 particle is the central particle (rc) from 11F
    // one sp3 particle is from a 5A in the 11F (sp3c_i/j_unc)
    // one sp3 particle is distinct from the 11F

    // Identification of sp3c_i_unc
    k=0;
    for (i=0; i<3; i++) {
        // Make sure the sp3 ring of 5A_i does not contain the central particle in the 11F
        if (sp3c[sp3c_i][i] != hc11F[n11F[f]][0]) { ;
            // Make sure none of the sp3 ring of 5A_i is in sp4 ring of the specified 6A
            for (j = 0; j < 4; j++) {
                if (sp3c[sp3c_i][i] == sp4c[the6A_i][j]) {
                    break;
                }
            }
            // if none of the sp3 ring is in sp4
            if (j == 4) {
                if (k >= 1) {
                    k++;
                    break;
                }
                sp3c_i_unc = sp3c[sp3c_i][i];
                k++;
            }
        }
    }
    if (k!=1) return 0;

    // Identification of sp3c_j_unc
    k=0;
    for (i=0; i<3; i++) {
        // Make sure the sp3 ring of 5A_i does not contain the central particle in the 11F
        if (sp3c[sp3c_j][i]!=hc11F[n11F[f]][0]) {
            // Make sure none of the sp3 ring of 5A_j is in sp4 ring of the specified 6A
            for (j = 0; j < 4; j++) {
                if (sp3c[sp3c_j][i] == sp4c[the6A_i][j]) {
                    break;
                }
            }
            if (j == 4) {
                if (k >= 1) {
                    k++;
                    break;
                }
                sp3c_j_unc = sp3c[sp3c_j][i];
                k++;
            }
        }
    }
    if (k!=1) return 0;

    // Try to identify the new particle in 5A_j which is not in 11F
    k=0;
    // loop over all 5A clusters of which hc13K[n13K[f]][0] is a member
    for (i=0; i<nmem_sp3c[hc11F[n11F[f]][0]]; ++i) {
        if (mem_sp3c[hc11F[n11F[f]][0]][i]!=sp3c_i && mem_sp3c[hc11F[n11F[f]][0]][i]!=sp3c_j) {
            if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][3] == sp3c[sp3c_i][3]) {
                if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][4] == sp3c[sp3c_i][4]) {
                    for (j = 0; j < 3; j++) {
                        if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][j] != hc11F[n11F[f]][0]) {
                            if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][j] != sp3c_i_unc) {
                                if (k == 1) {
                                    return 0;
                                }
                                else {
                                    tmp = sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][j];
                                    // check that tmp is not already in 11F
                                    for (l = 0; l < 11; l++) {
                                        if (tmp == hc11F[n11F[f]][l]) break;
                                    }
                                    if (l == 11) {
                                        ep[0] = tmp;
                                        eclus5A[0] = mem_sp3c[hc11F[n11F[f]][0]][i];
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
    // have found particle uncommon to 11F which forms a 5A cluster sharing spindles of sp3c_i and ring particles hc13K[n13K[f]][0] and sp3c_i_unc

    // Try to identify the new particle in 5A_j which is not in 11F
    k=0;
    for (i=0; i<nmem_sp3c[hc11F[n11F[f]][0]]; ++i) { // loop over all sp3c which hc13K[n13K[f]][0] is a member
        if (mem_sp3c[hc11F[n11F[f]][0]][i]!=sp3c_i && mem_sp3c[hc11F[n11F[f]][0]][i]!=sp3c_j) {
            if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][3] == sp3c[sp3c_j][3]) {
                if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][4] == sp3c[sp3c_j][4]) {
                    for (j = 0; j < 3; j++) {
                        if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][j] != hc11F[n11F[f]][0]) {
                            if (sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][j] != sp3c_j_unc) {
                                if (k == 1) {
                                    return 0;
                                }
                                else {
                                    tmp = sp3c[mem_sp3c[hc11F[n11F[f]][0]][i]][j];

                                    for (l = 0; l < 11; l++) {  // check temp not already in 11F
                                        if (tmp == hc11F[n11F[f]][l]) break;
                                    }
                                    if (l == 11) {
                                        ep[1] = tmp;
                                        eclus5A[1] = mem_sp3c[hc11F[n11F[f]][0]][i];
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
    // have found particle uncommon to 11F which forms a 5A cluster sharing spindles of sp3c_j and ring particles hc13K[n13K[f]][0] and sp3c_j_unc

    if(n13K[f] == m13K) {
        hc13K=resize_2D_int(hc13K,m13K,m13K+incrStatic,clusSize,-1);
        m13K=m13K+incrStatic;
    }
    // hc13K key: (11F, extra SP3 ring particle to make 5A #1, extra SP3 ring particle to make 5A #2)

    for (i=0; i<11; i++) hc13K[n13K[f]][i] = hc11F[n11F[f]][i];
    hc13K[n13K[f]][11]=ep[0];
    hc13K[n13K[f]][12]=ep[1];

    quickSort(&hc13K[n13K[f]][11],2);
    quickSort(&eclus5A[0],2);

    if(ach[hc13K[n13K[f]][1]]  == 'C') ach[hc13K[n13K[f]][1]] = ach_shell[hc13K[n13K[f]][1]] = 'B';
    if(ach[hc13K[n13K[f]][2]]  == 'C') ach[hc13K[n13K[f]][2]] = ach_shell[hc13K[n13K[f]][2]] = 'B';
    if(ach[hc13K[n13K[f]][3]]  == 'C') ach[hc13K[n13K[f]][3]] = ach_shell[hc13K[n13K[f]][3]] = 'B';
    if(ach[hc13K[n13K[f]][4]]  == 'C') ach[hc13K[n13K[f]][4]] = ach_shell[hc13K[n13K[f]][4]] = 'B';
    if(ach[hc13K[n13K[f]][5]]  == 'C') ach[hc13K[n13K[f]][5]] = ach_shell[hc13K[n13K[f]][5]] = 'B';
    if(ach[hc13K[n13K[f]][6]]  == 'C') ach[hc13K[n13K[f]][6]] = ach_shell[hc13K[n13K[f]][6]] = 'B';
    if(ach[hc13K[n13K[f]][7]]  == 'C') ach[hc13K[n13K[f]][7]] = ach_shell[hc13K[n13K[f]][7]] = 'B';
    if(ach[hc13K[n13K[f]][8]]  == 'C') ach[hc13K[n13K[f]][8]] = ach_shell[hc13K[n13K[f]][8]] = 'B';
    if(ach[hc13K[n13K[f]][9]]  == 'C') ach[hc13K[n13K[f]][9]] = ach_shell[hc13K[n13K[f]][9]] = 'B';
    if(ach[hc13K[n13K[f]][10]]  == 'C') ach[hc13K[n13K[f]][10]] = ach_shell[hc13K[n13K[f]][10]] = 'B';
    ach[hc13K[n13K[f]][0]] = ach_cen[hc13K[n13K[f]][0]] = 'O';
    ach[hc13K[n13K[f]][11]] = ach_shell[hc13K[n13K[f]][11]] = 'O';
    ach[hc13K[n13K[f]][12]] = ach_shell[hc13K[n13K[f]][12]] = 'O';

    return 1;
}

void Clusters_Get12B_13A(int f) { // Detect 12B & 13A D5h clusters together
    int i, j, k, l, m;
    int sp1, sp2;
    int sj1[5], sj2[5];
    int nSB1, nSB2;
    char *ach1, *ach1_cen, *ach1_shell, *ach2, *ach2_cen, *ach2_shell;
    int flg;
    int break_out;
    char errMsg[1000];
    int clusSize=12;


    ach1=malloc(N*sizeof(char));    if (ach1==NULL) { sprintf(errMsg,"Clusters_Get12B_13A(): ach1[] malloc out of memory\n");   Error(errMsg); }
    ach1_cen=malloc(N*sizeof(char));    if (ach1_cen==NULL) { sprintf(errMsg,"Clusters_Get12B_13A(): ach1_cen[] malloc out of memory\n");   Error(errMsg); }
    ach1_shell=malloc(N*sizeof(char));  if (ach1_shell==NULL) { sprintf(errMsg,"Clusters_Get12B_13A(): ach1_shell[] malloc out of memory\n");   Error(errMsg); }
    ach2=malloc(N*sizeof(char));    if (ach2==NULL) { sprintf(errMsg,"Clusters_Get12B_13A(): ach2[] malloc out of memory\n");   Error(errMsg); }
    ach2_cen=malloc(N*sizeof(char));    if (ach2_cen==NULL) { sprintf(errMsg,"Clusters_Get12B_13A(): ach2_cen[] malloc out of memory\n");   Error(errMsg); }
    ach2_shell=malloc(N*sizeof(char));  if (ach2_shell==NULL) { sprintf(errMsg,"Clusters_Get12B_13A(): ach2_shell[] malloc out of memory\n");   Error(errMsg); }
    for (i=0; i<N; ++i) ach1[i] = ach1_cen[i] = ach1_shell[i] = ach2[i] = ach2_cen[i] = ach2_shell[i] = 'C';

    for (i=0; i<nsp5c[f]; ++i) { //first 7A
        sp1 = sp5c[i][5];
        sp2 = sp5c[i][6];
        nSB1 = nSB2 = 0; // count up spindle bonds

        for (j=i+1; j<nsp5c[f]; ++j) { // second 7A
            flg = sp1 == sp5c[j][5] && Bonds_BondCheck(sp2,sp5c[j][6]);
            flg = flg || (sp1 == sp5c[j][6] && Bonds_BondCheck(sp2,sp5c[j][5]));
            if (flg==1) {
                if (nSB1>=5) {
                    nSB1++;
                    break;
                }
                sj1[nSB1++] = j;
            }
        }

        if (nSB1 == 5 && do13A==1) {     // possibly found 13A, definately found 12B, now establish status
            for (j=i+1; j<nsp5c[f]; ++j) {
                if (sp1 == sp5c[j][5] || sp1 == sp5c[j][6]) {
                    for (k=0; k<5; ++k) {
                        for (l=0; l<5; ++l) {
                            if (sp5c[i][k] == sp5c[j][l]) break;
                        }
                        if (l<5) break;
                    }
                    if(k==5) { // got 13A, Check all sp5c[j][] - sp1 spindles are less than i
                        for (k=0; k<i; ++k) {
                            for(l=0; l<5; ++l) {
                                if (sp5c[j][l] == sp5c[k][5] && sp1 == sp5c[k][6]) break;
                                if (sp5c[j][l] == sp5c[k][6] && sp1 == sp5c[k][5]) break;
                            }
                            if(l<5) break; // index k < i present
                        }
                        if(k==i) break; // no index k < i present
                    }
                }
            }
            if (j<nsp5c[f]) { // 13A found
                if (n13A[f] == m13A) {
                    hc13A=resize_2D_int(hc13A,m13A,m13A+incrStatic,clusSize+1,-1);
                    m13A=m13A+incrStatic;
                }

                hc13A[n13A[f]][0] = sp1;
                k = 1;
                if(sp5c[i][5] != sp1) hc13A[n13A[f]][k++] = sp5c[i][5];
                if(sp5c[i][6] != sp1) hc13A[n13A[f]][k++] = sp5c[i][6];
                if(sp5c[j][5] != sp1) hc13A[n13A[f]][k++] = sp5c[j][5];
                if(sp5c[j][6] != sp1) hc13A[n13A[f]][k] = sp5c[j][6];
                hc13A[n13A[f]][3] = sp5c[i][0];
                hc13A[n13A[f]][4] = sp5c[i][1];
                hc13A[n13A[f]][5] = sp5c[i][2];
                hc13A[n13A[f]][6] = sp5c[i][3];
                hc13A[n13A[f]][7] = sp5c[i][4];
                hc13A[n13A[f]][8] = sp5c[j][0];
                hc13A[n13A[f]][9] = sp5c[j][1];
                hc13A[n13A[f]][10] = sp5c[j][2];
                hc13A[n13A[f]][11] = sp5c[j][3];
                hc13A[n13A[f]][12] = sp5c[j][4];
                quickSort(&hc13A[n13A[f]][1],12);

                Clust_Write_13A(f, ach2, ach2_cen, ach2_shell);
            }
        }

        for (j=i+1; j<nsp5c[f]; ++j) {
            flg = sp2 == sp5c[j][5] && Bonds_BondCheck(sp1,sp5c[j][6]);
            flg = flg || (sp2 == sp5c[j][6] && Bonds_BondCheck(sp1,sp5c[j][5]));
            if (flg==1) {
                if (nSB2>=5) {
                    nSB2++;
                    break;
                }
                sj2[nSB2++] = j;
            }
        }

        if(nSB2 == 5 && do13A==1) { // possibly found 13A, definately found 12B, now establish status
            for (j=i+1; j<nsp5c[f]; ++j) {
                if (sp2 == sp5c[j][5] || sp2 == sp5c[j][6]) {
                    for (k=0; k<5; ++k) {
                        for (l=0; l<5; ++l) {
                            if (sp5c[i][k] == sp5c[j][l]) break;
                        }
                        if(l<5) break;
                    }
                    if (k==5) { // Check all sp5c[j][] - sp2 spindles are less than i
                        for (k=0; k<i; ++k) {
                            for (l=0; l<5; ++l) {
                                if (sp5c[j][l] == sp5c[k][5] && sp2 == sp5c[k][6]) break;
                                if (sp5c[j][l] == sp5c[k][6] && sp2 == sp5c[k][5]) break;
                            }
                            if (l<5) break;
                        }
                        if (k==i) break;
                    }
                }
            }
            if (j<nsp5c[f]){ // 13A found, is it unique
                if (n13A[f] == m13A) {
                    hc13A=resize_2D_int(hc13A,m13A,m13A+incrStatic,clusSize+1,-1);
                    m13A=m13A+incrStatic;
                }

                hc13A[n13A[f]][0] = sp2;
                k = 1;
                if(sp5c[i][5] != sp2) hc13A[n13A[f]][k++] = sp5c[i][5];
                if(sp5c[i][6] != sp2) hc13A[n13A[f]][k++] = sp5c[i][6];
                if(sp5c[j][5] != sp2) hc13A[n13A[f]][k++] = sp5c[j][5];
                if(sp5c[j][6] != sp2) hc13A[n13A[f]][k] = sp5c[j][6];
                hc13A[n13A[f]][3] = sp5c[i][0];
                hc13A[n13A[f]][4] = sp5c[i][1];
                hc13A[n13A[f]][5] = sp5c[i][2];
                hc13A[n13A[f]][6] = sp5c[i][3];
                hc13A[n13A[f]][7] = sp5c[i][4];
                hc13A[n13A[f]][8] = sp5c[j][0];
                hc13A[n13A[f]][9] = sp5c[j][1];
                hc13A[n13A[f]][10] = sp5c[j][2];
                hc13A[n13A[f]][11] = sp5c[j][3];
                hc13A[n13A[f]][12] = sp5c[j][4];
                quickSort(&hc13A[n13A[f]][1],12);

                Clust_Write_13A(f, ach2, ach2_cen, ach2_shell);
            }
        }

        if ((nSB1 > 5) && (nSB2 > 5)) continue;

        for (j=0; j<i; ++j) { // keep looking for 12B
            flg = sp1 == sp5c[j][5] && Bonds_BondCheck(sp2,sp5c[j][6]);
            flg = flg || (sp1 == sp5c[j][6] && Bonds_BondCheck(sp2,sp5c[j][5]));
            if (flg==1) {
                if (nSB1 >= 5) {
                    nSB1++;
                    break;
                }
                sj1[nSB1++] = j;
            }
        }

        if (nSB1 == 5) {
            if (n12B[f]==m12B) {
                hc12B=resize_2D_int(hc12B,m12B,m12B+incrStatic,clusSize,-1);
                m12B=m12B+incrStatic;
            }
            hc12B[n12B[f]][0] = sp1;
            hc12B[n12B[f]][1] = sp2;
            hc12B[n12B[f]][2] = sp5c[i][0];
            hc12B[n12B[f]][3] = sp5c[i][1];
            hc12B[n12B[f]][4] = sp5c[i][2];
            hc12B[n12B[f]][5] = sp5c[i][3];
            hc12B[n12B[f]][6] = sp5c[i][4];

            m = 7;
            break_out=0;
            for (j=0; j<5; ++j) {
                for (k=0; k<5; ++k) {
                    for (l=0; l<m; ++l) {
                        if(hc12B[n12B[f]][l] == sp5c[sj1[j]][k]) break;
                    }
                    if (l==m) {
                        if (m==12) {
                            break_out=1;
                            break;
                        }
                        hc12B[n12B[f]][m] = sp5c[sj1[j]][k];
                        m++;
                    }
                }
                if (break_out==1) break;
            }
            if (break_out==1 || m<12) continue;

            quickSort(&hc12B[n12B[f]][2], 5);
            quickSort(&hc12B[n12B[f]][7],5);

            Clust_Write_12B(f, ach1, ach1_cen, ach1_shell);
        }

        for (j=0; j<i; ++j) {
            flg = sp2 == sp5c[j][5] && Bonds_BondCheck(sp1,sp5c[j][6]);
            flg = flg || (sp2 == sp5c[j][6] && Bonds_BondCheck(sp1,sp5c[j][5]));
            if (flg==1) {
                if (nSB2 >= 5) {
                    nSB2++;
                    break;
                }
                sj2[nSB2++] = j;
            }
        }

        if(nSB2 == 5) {
            if (n12B[f]==m12B) {
                hc12B=resize_2D_int(hc12B,m12B,m12B+incrStatic,clusSize,-1);
                m12B=m12B+incrStatic;
            }

            hc12B[n12B[f]][0] = sp2;
            hc12B[n12B[f]][1] = sp1;
            hc12B[n12B[f]][2] = sp5c[i][0];
            hc12B[n12B[f]][3] = sp5c[i][1];
            hc12B[n12B[f]][4] = sp5c[i][2];
            hc12B[n12B[f]][5] = sp5c[i][3];
            hc12B[n12B[f]][6] = sp5c[i][4];

            m = 7;
            break_out=0;
            for(j=0; j<5; ++j){
                for(k=0; k<5; ++k){
                    for(l=0; l<m; ++l) if(hc12B[n12B[f]][l] == sp5c[sj2[j]][k]) break;
                    if(l==m) {
                        if (m==12) {
                            break_out=1;
                            break;
                        }
                        hc12B[n12B[f]][m] = sp5c[sj2[j]][k];
                        m++;
                    }
                }
                if (break_out==1) break;
            }
            if (break_out==1 || m<12) continue;

            quickSort(&hc12B[n12B[f]][2], 5);
            quickSort(&hc12B[n12B[f]][7],5);

            Clust_Write_12B(f, ach1, ach1_cen, ach1_shell);

        }
    }

    for(i=0; i<N; ++i) {
        s12B[i]=ach1[i];
        s12B_cen[i]=ach1_cen[i];
        s12B_shell[i]=ach1_shell[i];
        s13A[i]=ach2[i];
        s13A_cen[i]=ach2_cen[i];
        s13A_shell[i]=ach2_shell[i];
    }
    free(ach1);
    free(ach1_cen);
    free(ach1_shell);
    free(ach2);
    free(ach2_cen);
    free(ach2_shell);
}

void Clust_Write_12B(int f, char *ach1, char *ach1_cen, char *ach1_shell) {

    ach1[hc12B[n12B[f]][0]] = ach1_cen[hc12B[n12B[f]][0]] = 'O';
    ach1[hc12B[n12B[f]][2]] = ach1_shell[hc12B[n12B[f]][2]] = 'O';
    ach1[hc12B[n12B[f]][3]] = ach1_shell[hc12B[n12B[f]][3]] = 'O';
    ach1[hc12B[n12B[f]][4]] = ach1_shell[hc12B[n12B[f]][4]] = 'O';
    ach1[hc12B[n12B[f]][5]] = ach1_shell[hc12B[n12B[f]][5]] = 'O';
    ach1[hc12B[n12B[f]][6]] = ach1_shell[hc12B[n12B[f]][6]] = 'O';
    if(ach1[hc12B[n12B[f]][1]] == 'C') ach1[hc12B[n12B[f]][1]] = ach1_shell[hc12B[n12B[f]][1]] = 'B';
    if(ach1[hc12B[n12B[f]][7]] == 'C') ach1[hc12B[n12B[f]][7]] = ach1_shell[hc12B[n12B[f]][7]] = 'B';
    if(ach1[hc12B[n12B[f]][8]] == 'C') ach1[hc12B[n12B[f]][8]] = ach1_shell[hc12B[n12B[f]][8]] = 'B';
    if(ach1[hc12B[n12B[f]][9]] == 'C') ach1[hc12B[n12B[f]][9]] = ach1_shell[hc12B[n12B[f]][9]] = 'B';
    if(ach1[hc12B[n12B[f]][10]] == 'C') ach1[hc12B[n12B[f]][10]] = ach1_shell[hc12B[n12B[f]][10]] = 'B';
    if(ach1[hc12B[n12B[f]][11]] == 'C') ach1[hc12B[n12B[f]][11]] = ach1_shell[hc12B[n12B[f]][11]] = 'B';

    ++n12B[f];
}

void Clust_Write_13A(int f, char *ach2, char *ach2_cen, char *ach2_shell) {

    if(ach2[hc13A[n13A[f]][3]] == 'C') ach2[hc13A[n13A[f]][3]] = ach2_shell[hc13A[n13A[f]][3]] = 'B';
    if(ach2[hc13A[n13A[f]][4]] == 'C') ach2[hc13A[n13A[f]][4]] = ach2_shell[hc13A[n13A[f]][4]] = 'B';
    if(ach2[hc13A[n13A[f]][5]] == 'C') ach2[hc13A[n13A[f]][5]] = ach2_shell[hc13A[n13A[f]][5]] = 'B';
    if(ach2[hc13A[n13A[f]][6]] == 'C') ach2[hc13A[n13A[f]][6]] = ach2_shell[hc13A[n13A[f]][6]] = 'B';
    if(ach2[hc13A[n13A[f]][7]] == 'C') ach2[hc13A[n13A[f]][7]] = ach2_shell[hc13A[n13A[f]][7]] = 'B';
    if(ach2[hc13A[n13A[f]][8]] == 'C') ach2[hc13A[n13A[f]][8]] = ach2_shell[hc13A[n13A[f]][8]] = 'B';
    if(ach2[hc13A[n13A[f]][9]] == 'C') ach2[hc13A[n13A[f]][9]] = ach2_shell[hc13A[n13A[f]][9]] = 'B';
    if(ach2[hc13A[n13A[f]][10]] == 'C') ach2[hc13A[n13A[f]][10]] = ach2_shell[hc13A[n13A[f]][10]] = 'B';
    if(ach2[hc13A[n13A[f]][11]] == 'C') ach2[hc13A[n13A[f]][11]] = ach2_shell[hc13A[n13A[f]][11]] = 'B';
    if(ach2[hc13A[n13A[f]][12]] == 'C') ach2[hc13A[n13A[f]][12]] = ach2_shell[hc13A[n13A[f]][12]] = 'B';
    ach2[hc13A[n13A[f]][0]] = ach2_cen[hc13A[n13A[f]][0]] = 'O';
    ach2[hc13A[n13A[f]][1]] = ach2_shell[hc13A[n13A[f]][1]] = 'O';
    ach2[hc13A[n13A[f]][2]] = ach2_shell[hc13A[n13A[f]][2]] = 'O';

    ++n13A[f];
}

void Clusters_Get13B_D5h(int f) {   // Detect 13B D5h clusters, i.e. twisted icosahedra
    int cp;
    int i, j, k, l, m;
    int flg;
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int clusSize=13;

    cp=-1;
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get13B_D5h(): ach[] malloc out of memory\n"); Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_Get13B_D5h(): ach_cen[] malloc out of memory\n"); Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_Get13B_D5h(): ach_shell[] malloc out of memory\n"); Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for(i=0; i<nsp5c[f]-1; ++i){
        // POSSIBLE IMPROVEMENT: use 7A clusters of spindles of 7A_i
        for(j=i+1; j<nsp5c[f]; ++j){
            flg = 0;
            if (sp5c[i][5]==sp5c[j][5] && sp5c[i][6]!=sp5c[j][6]) {
                if (Bonds_BondCheck(sp5c[i][6], sp5c[j][6])==1) continue;
                cp = sp5c[i][5];
                flg = 1;
            }
            if (sp5c[i][6]==sp5c[j][6] && sp5c[i][5]!=sp5c[j][5]) {
                if (Bonds_BondCheck(sp5c[i][5],sp5c[j][5])==1) continue;
                cp = sp5c[i][6];
                flg = 1;
            }
            if (sp5c[i][5]==sp5c[j][6] && sp5c[i][6]!=sp5c[j][5]) {
                if (Bonds_BondCheck(sp5c[i][6],sp5c[j][5])==1) continue;
                cp = sp5c[i][5];
                flg = 1;
            }
            if (sp5c[i][6]==sp5c[j][5] && sp5c[i][5]!=sp5c[j][6]) {
                if (Bonds_BondCheck(sp5c[i][5],sp5c[j][6])==1) continue;
                cp = sp5c[i][6];
                flg = 1;
            }
            if (flg==1) {
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if(sp5c[i][k] == sp5c[j][l]) break;
                    }
                    if(l<5) break;
                }
                if (k<5) continue;
                for (k=0; k<5; ++k) {
                    m = 0;
                    for (l=0; l<5; ++l) {
                        if (Bonds_BondCheck(sp5c[i][k],sp5c[j][l])==1) ++m;
                        if (m == 2) break;
                    }
                    if(m != 1) break;
                }
                if (k<5) continue;
                // ERROR: need to check converse, i.e. every SP5 from 7A_j bonded to one from SP5 from 7A_i

                if (n13B[f] == m13B) {
                    hc13B=resize_2D_int(hc13B,m13B,m13B+incrStatic,clusSize,-1);
                    m13B=m13B+incrStatic;
                }

                hc13B[n13B[f]][0] = cp;
                k = 1;
                if(sp5c[i][5] != cp) hc13B[n13B[f]][k++] = sp5c[i][5];
                if(sp5c[i][6] != cp) hc13B[n13B[f]][k++] = sp5c[i][6];
                if(sp5c[j][5] != cp) hc13B[n13B[f]][k++] = sp5c[j][5];
                if(sp5c[j][6] != cp) hc13B[n13B[f]][k++] = sp5c[j][6];
                for(l=0; l<5; ++l) hc13B[n13B[f]][k++] = sp5c[i][l];
                for(l=0; l<5; ++l) hc13B[n13B[f]][k++] = sp5c[j][l];
                quickSort(&hc13B[n13B[f]][1],2);
                quickSort(&hc13B[n13B[f]][3],10);

                Cluster_Write_13B(f, ach, ach_cen, ach_shell);
            }
        }
    }
    for(i=0; i<N; ++i) {
        s13B[i]=ach[i];
        s13B_cen[i]=ach_cen[i];
        s13B_shell[i]=ach_shell[i];
    }
    free(ach);
    free(ach_cen);
    free(ach_shell);
}

void Cluster_Write_13B(int f, char *ach, char *ach_cen, char *ach_shell) {
    if(ach[hc13B[n13B[f]][3]] == 'C') ach[hc13B[n13B[f]][3]] = ach_shell[hc13B[n13B[f]][3]] = 'B';
    if(ach[hc13B[n13B[f]][4]] == 'C') ach[hc13B[n13B[f]][4]] = ach_shell[hc13B[n13B[f]][4]] = 'B';
    if(ach[hc13B[n13B[f]][5]] == 'C') ach[hc13B[n13B[f]][5]] = ach_shell[hc13B[n13B[f]][5]] = 'B';
    if(ach[hc13B[n13B[f]][6]] == 'C') ach[hc13B[n13B[f]][6]] = ach_shell[hc13B[n13B[f]][6]] = 'B';
    if(ach[hc13B[n13B[f]][7]] == 'C') ach[hc13B[n13B[f]][7]] = ach_shell[hc13B[n13B[f]][7]] = 'B';
    if(ach[hc13B[n13B[f]][8]] == 'C') ach[hc13B[n13B[f]][8]] = ach_shell[hc13B[n13B[f]][8]] = 'B';
    if(ach[hc13B[n13B[f]][9]] == 'C') ach[hc13B[n13B[f]][9]] = ach_shell[hc13B[n13B[f]][9]] = 'B';
    if(ach[hc13B[n13B[f]][10]] == 'C') ach[hc13B[n13B[f]][10]] = ach_shell[hc13B[n13B[f]][10]] = 'B';
    if(ach[hc13B[n13B[f]][11]] == 'C') ach[hc13B[n13B[f]][11]] = ach_shell[hc13B[n13B[f]][11]] = 'B';
    if(ach[hc13B[n13B[f]][12]] == 'C') ach[hc13B[n13B[f]][12]] = ach_shell[hc13B[n13B[f]][12]] = 'B';
    ach[hc13B[n13B[f]][0]] = ach_cen[hc13B[n13B[f]][0]] = 'O';
    ach[hc13B[n13B[f]][1]] = ach_shell[hc13B[n13B[f]][1]] = 'O';
    ach[hc13B[n13B[f]][2]] = ach_shell[hc13B[n13B[f]][2]] = 'O';

    ++n13B[f];
}

void Clusters_GetFCC(int f) {   // Detect 13 particle FCC clusters
    int i, j, j2, k, l, m, n;
    int i1, i2, i3;
    int cp, bpi, bpj, nbpi, nbpj;
    int flg1, flg2, flg3;
    int l_clust_type; // 0 if l-clust is sp3b, 1 if is sp3c
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int clusSize=13;

    cp=bpi=bpj=nbpi=nbpj=i3=-1;
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_GetFCC(): ach[] malloc out of memory\n"); Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_GetFCC(): ach_cen[] malloc out of memory\n"); Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_GetFCC(): ach_shell[] malloc out of memory\n"); Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for (i=0; i<nsp3b[f]-2; ++i) { // loop over all sp3b_i
        for (j2=0; j2<3; j2++) {
            for (j=0; j<nmem_sp3b[sp3b[i][j2]]-1; ++j) { // loop over all sp3b_j
                if (mem_sp3b[sp3b[i][j2]][j]<=i) continue;
                if (Bonds_BondCheck(sp3b[i][3], sp3b[mem_sp3b[sp3b[i][j2]][j]][3])==0) continue; // spindle spindle bond
                m = n = 0;
                for (k=0; k<4; ++k) {
                    for (l=0; l<4; ++l) {
                        if (sp3b[i][k] == sp3b[mem_sp3b[sp3b[i][j2]][j]][l]) {
                            if (k == 3 || l == 3) m +=2;
                            cp = sp3b[i][k];
                            ++m;
                        }
                    }
                }
                if(m != 1) continue; // one common particle between sp3b_i and sp3b_j

                for (k=0; k<nFCC[f]; ++k) {
                    if (hcFCC[k][0] == cp) break;   // check for other degenerate FCC clusters which have cp here
                }
                if (k < nFCC[f]) continue;  // found this fcc cluster before

                for (k=0; k<3; ++k) {
                    if (sp3b[i][k] == cp) continue;
                    for (l=0; l<3; ++l) {
                        if (sp3b[mem_sp3b[sp3b[i][j2]][j]][l] == cp) continue;
                        if (Bonds_BondCheck(sp3b[i][k],sp3b[mem_sp3b[sp3b[i][j2]][j]][l])==1) {
                            bpi = sp3b[i][k];
                            bpj = sp3b[mem_sp3b[sp3b[i][j2]][j]][l];
                            ++n;
                        }
                    }
                }
                if(n != 1) continue;    // one bond between one pair of atoms from SP3_i and SP3_j

                // we may have an FCC cluster, build it into hcFCC then overwrite it later if it aint
                if (nFCC[f] == mFCC) {
                    hcFCC=resize_2D_int(hcFCC,mFCC,mFCC+incrStatic,clusSize,-1);
                    mFCC=mFCC+incrStatic;
                }

                for (k=0; k<3; ++k) {
                    if (sp3b[i][k]!=cp && sp3b[i][k]!=bpi) nbpi=sp3b[i][k];
                    if (sp3b[mem_sp3b[sp3b[i][j2]][j]][k]!=cp && sp3b[mem_sp3b[sp3b[i][j2]][j]][k]!=bpj) nbpj=sp3b[mem_sp3b[sp3b[i][j2]][j]][k];
                }
                hcFCC[nFCC[f]][0] = cp;
                hcFCC[nFCC[f]][1] = nbpi;
                hcFCC[nFCC[f]][2] = bpi;
                hcFCC[nFCC[f]][3] = sp3b[i][3];
                hcFCC[nFCC[f]][4] = sp3b[mem_sp3b[sp3b[i][j2]][j]][3];
                hcFCC[nFCC[f]][5] = nbpj;
                hcFCC[nFCC[f]][6] = bpj;   // store first two clusters

                for (k=j+1; k<nmem_sp3b[sp3b[i][j2]]; ++k) { // loop over all sp3b_k
                    if (mem_sp3b[sp3b[i][j2]][k]<=i) continue;
                    for (l=0; l<nFCC[f]; ++l) {
                        if (hcFCC[l][0] == cp) break;   // check for other degenerate FCC clusters which have cp here
                    }
                    if (l < nFCC[f]) break; // found this fcc cluster before

                    if (sp3b[mem_sp3b[sp3b[i][j2]][k]][3]==cp) continue;    // check sp3b_k spindle isnt common particle
                    if (!(Bonds_BondCheck(sp3b[i][3],sp3b[mem_sp3b[sp3b[i][j2]][k]][3])==1 && Bonds_BondCheck(sp3b[mem_sp3b[sp3b[i][j2]][j]][3], sp3b[mem_sp3b[sp3b[i][j2]][k]][3])==1)) continue;  // check sp3b_k spindle is bonded to sp3b_i sp3b_j spindles
                    i1=i2=-1;
                    for (l=0; l<3; ++l) {
                        if (sp3b[mem_sp3b[sp3b[i][j2]][k]][l] == cp) {
                            i1 = l; // identity of common particle in SP3_k
                        }
                        else {
                            i2 = sp3b[mem_sp3b[sp3b[i][j2]][k]][l]; // identity of a loose particle in SP3_k
                            flg1 = i2 == nbpi || i2 == bpi || i2 == nbpj || i2 == bpj;
                            flg1 = flg1 || i2 == sp3b[i][3] || i2 == sp3b[mem_sp3b[sp3b[i][j2]][j]][3];
                            if(flg1==1) break;
                        }
                    }
                    if (l<3 || i1==-1 || i2==-1) continue;  // SP3_k must have 2 new particles

                    m = 0;
                    flg1 = flg2 = 0;
                    for (l=0; l<3; ++l) {
                        if (l == i1) continue;
                        n = sp3b[mem_sp3b[sp3b[i][j2]][k]][l];
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

                    for (l=0; l<nsp3b[f]; ++l) { // loop over all sp3b_l
                        if (l==i || l==mem_sp3b[sp3b[i][j2]][j] || l==mem_sp3b[sp3b[i][j2]][k]) continue;   // must be new sp3b
                        if (sp3b[l][3] != cp) continue; // must have spindle as common particle
                        for (m=0; m<3; ++m) {
                            n = sp3b[l][m];
                            flg1 = n == nbpi || n == bpi || n == bpj;
                            flg1 = flg1 || n == nbpj || n == i3 || n == i2;
                            flg1 = flg1 || n == sp3b[i][3] || n == sp3b[mem_sp3b[sp3b[i][j2]][j]][3];
                            flg1 = flg1 || n == sp3b[mem_sp3b[sp3b[i][j2]][k]][3];
                            if (flg1==1) break;
                        }
                        if (m<3) continue; // the SP3_l ring particles must be new

                        flg1 = flg2 = flg3 = 0;
                        for (m=0; m<3; ++m) {
                            n = sp3b[l][m];
                            if (Bonds_BondCheck(n,bpi)==1 && Bonds_BondCheck(n,bpj)==1) {
                                if (flg1==1) break;
                                flg1 = 1;
                                hcFCC[nFCC[f]][10] = n;
                            }
                            if (Bonds_BondCheck(n,nbpj)==1 && Bonds_BondCheck(n,i3)==1) {
                                if(flg2==1) break;
                                flg2 = 1;
                                hcFCC[nFCC[f]][11] = n;
                            }
                            if (Bonds_BondCheck(n,i2)==1 && Bonds_BondCheck(n,nbpi)==1) {
                                if(flg3==1) break;
                                flg3 = 1;
                                hcFCC[nFCC[f]][12] = n;
                            }
                        }
                        if(m<3) continue;  // SP3_l particles must be bonded correctly to six-membered ring
                        if(flg1==1 && flg2==1 && flg3==1) break;   // SP3_l particles must be bonded correctly to six-membered ring
                    } // end l loop
                    l_clust_type=0;
                    if (l==nsp3b[f]) { // required SP3b cluster not found
                        for (l=0; l<nsp3c[f]; ++l) { // loop over all sp3c_l
                            if (sp3c[l][3]!=cp && sp3c[l][4]!=cp) continue;  // common particle must be a spindle
                            for (m=0; m<3; ++m) {
                                n = sp3c[l][m];
                                flg1 = n == nbpi || n == bpi || n == bpj;
                                flg1 = flg1 || n == nbpj || n == i3 || n == i2;
                                flg1 = flg1 || n == sp3b[i][3] || n == sp3b[mem_sp3b[sp3b[i][j2]][j]][3];
                                flg1 = flg1 || n == sp3b[mem_sp3b[sp3b[i][j2]][k]][3];
                                if (flg1==1) break;
                            }
                            if (m<3) continue; // the SP3_l ring particles must be new

                            flg1 = flg2 = flg3 = 0;
                            for (m=0; m<3; ++m) {
                                n = sp3c[l][m];
                                if (Bonds_BondCheck(n,bpi)==1 && Bonds_BondCheck(n,bpj)==1) {
                                    if (flg1==1) break;
                                    flg1 = 1;
                                    hcFCC[nFCC[f]][10] = n;
                                }
                                if (Bonds_BondCheck(n,nbpj)==1 && Bonds_BondCheck(n,i3)==1) {
                                    if (flg2==1) break;
                                    flg2 = 1;
                                    hcFCC[nFCC[f]][11] = n;
                                }
                                if (Bonds_BondCheck(n,i2)==1 && Bonds_BondCheck(n,nbpi)==1) {
                                    if(flg3==1) break;
                                    flg3 = 1;
                                    hcFCC[nFCC[f]][12] = n;
                                }
                            }
                            if (m<3) continue;  // SP3_l particles must be bonded correctly to six-membered ring
                            if (flg1==1 && flg2==1 && flg3==1) break;   // SP3_l particles must be bonded correctly to six-membered ring
                        }
                        if (l==nsp3c[f]) continue;
                        l_clust_type=1;
                    } // End if statement for SP3c search loop

                    // We've now found an FCC cluster
                    if (nFCC[f] == mFCC) { hcFCC=resize_2D_int(hcFCC,mFCC,mFCC+incrStatic,clusSize,-1);
                        mFCC=mFCC+incrStatic;
                    }

                    hcFCC[nFCC[f]][7] = sp3b[mem_sp3b[sp3b[i][j2]][k]][3];
                    hcFCC[nFCC[f]][8] = i3;
                    hcFCC[nFCC[f]][9] = i2;
                    quickSort(&hcFCC[nFCC[f]][1],12);

                    Cluster_Write_FCC(f, ach, ach_cen, ach_shell);
                }
            }
        }
    }
    for(i=0; i<N; ++i) {
        sFCC[i]=ach[i];
        sFCC_cen[i]=ach_cen[i];
        sFCC_shell[i]=ach_shell[i];
    }
    free(ach);
    free(ach_cen);
    free(ach_shell);
}

void Cluster_Write_FCC(int f, char *ach, char *ach_cen, char *ach_shell) {
    if(ach[hcFCC[nFCC[f]][1]] == 'C') ach[hcFCC[nFCC[f]][1]] = ach_shell[hcFCC[nFCC[f]][1]] = 'B';
    if(ach[hcFCC[nFCC[f]][2]] == 'C') ach[hcFCC[nFCC[f]][2]] = ach_shell[hcFCC[nFCC[f]][2]] = 'B';
    if(ach[hcFCC[nFCC[f]][5]] == 'C') ach[hcFCC[nFCC[f]][5]] = ach_shell[hcFCC[nFCC[f]][5]] = 'B';
    if(ach[hcFCC[nFCC[f]][6]] == 'C') ach[hcFCC[nFCC[f]][6]] = ach_shell[hcFCC[nFCC[f]][6]] = 'B';
    if(ach[hcFCC[nFCC[f]][8]] == 'C') ach[hcFCC[nFCC[f]][8]] = ach_shell[hcFCC[nFCC[f]][8]] = 'B';
    if(ach[hcFCC[nFCC[f]][9]] == 'C') ach[hcFCC[nFCC[f]][9]] = ach_shell[hcFCC[nFCC[f]][9]] = 'B';
    ach[hcFCC[nFCC[f]][3]] = ach_shell[hcFCC[nFCC[f]][3]] = 'O';
    ach[hcFCC[nFCC[f]][4]] = ach_shell[hcFCC[nFCC[f]][4]] = 'O';
    ach[hcFCC[nFCC[f]][7]] = ach_shell[hcFCC[nFCC[f]][7]] = 'O';
    ach[hcFCC[nFCC[f]][10]] = ach_shell[hcFCC[nFCC[f]][10]] = 'O';
    ach[hcFCC[nFCC[f]][11]] = ach_shell[hcFCC[nFCC[f]][11]] = 'O';
    ach[hcFCC[nFCC[f]][12]] = ach_shell[hcFCC[nFCC[f]][12]] = 'O';
    ach[hcFCC[nFCC[f]][0]] = ach_cen[hcFCC[nFCC[f]][0]] = 'F';

    ++nFCC[f];
}

void Clusters_GetHCP(int f) {   // Detect 13 particle HCP clusters
    int i, j, j2, k, l, m, n;
    int ia[2], ja[2], ka[2];
    int cp, x;
    int h1i, h1j, h2i, h2j;
    int flg1, flg2, flg3;
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int clusSize=13;

    cp=-1;
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_GetHCP(): ach[] malloc out of memory\n"); Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_GetHCP(): ach_cen[] malloc out of memory\n"); Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_GetHCP(): ach_shell[] malloc out of memory\n"); Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for (i=0; i<nsp3c[f]-2; ++i) { // loop over all sp3c_i
        for (j2=0; j2<3; j2++) {
        for (j=0; j<nmem_sp3c[sp3c[i][j2]]-1; ++j) { // loop over all sp3c_j
            if (mem_sp3c[sp3c[i][j2]][j]<=i) continue;
            flg1 = Bonds_BondCheck(sp3c[i][3],sp3c[mem_sp3c[sp3c[i][j2]][j]][3]) && Bonds_BondCheck(sp3c[i][4],sp3c[mem_sp3c[sp3c[i][j2]][j]][4]);
            flg2 = Bonds_BondCheck(sp3c[i][3],sp3c[mem_sp3c[sp3c[i][j2]][j]][4]) && Bonds_BondCheck(sp3c[i][4],sp3c[mem_sp3c[sp3c[i][j2]][j]][3]);
            if (!(flg1==1 || flg2==1)) continue; // sp3c_i has both spindles bonded to sp3c_j spindles
            if (flg1==1) {
                h1i = sp3c[i][3]; // h1i bonded to h1j, h2i bonded to h2j
                h1j = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                h2i = sp3c[i][4];
                h2j = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
            }
            else {
                h1i = sp3c[i][3]; // h1i bonded to h1j, h2i bonded to h2j
                h1j = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                h2i = sp3c[i][4];
                h2j = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
            }

            flg1 = sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
            flg2 = sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
            flg1 = flg1 || flg2;
            if (flg1==1) continue; // spindles particles must be distinct, no common spindles

            m = 0;
            for (k=0; k<3; ++k) {
                for(l=0; l<3; ++l) {
                    if (sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) {
                        cp = sp3c[i][k];
                        ++m;
                    }
                }
            }
            if (m != 1) continue; // 1 common central particle between SP3_i and SP3_j
            if (cp!=sp3c[i][j2]) continue;

            for (k=j+1; k<nmem_sp3c[sp3c[i][j2]]; ++k) { // loop over all sp3c_k
                if (mem_sp3c[sp3c[i][j2]][k]<=i) continue;
                flg1 = Bonds_BondCheck(h1i,sp3c[mem_sp3c[sp3c[i][j2]][k]][3]) && Bonds_BondCheck(h1j,sp3c[mem_sp3c[sp3c[i][j2]][k]][3]);
                flg2 = Bonds_BondCheck(h2i,sp3c[mem_sp3c[sp3c[i][j2]][k]][4]) && Bonds_BondCheck(h2j,sp3c[mem_sp3c[sp3c[i][j2]][k]][4]);
                flg1 = flg1 && flg2;
                flg2 = Bonds_BondCheck(h1i,sp3c[mem_sp3c[sp3c[i][j2]][k]][4]) && Bonds_BondCheck(h1j,sp3c[mem_sp3c[sp3c[i][j2]][k]][4]);
                flg3 = Bonds_BondCheck(h2i,sp3c[mem_sp3c[sp3c[i][j2]][k]][3]) && Bonds_BondCheck(h2j,sp3c[mem_sp3c[sp3c[i][j2]][k]][3]);
                flg2 = flg2 && flg3;
                flg1 = flg1 || flg2;
                if (flg1==0) continue; // spindle particles not bonded correctly
                flg1 = sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][k]][3] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][k]][4];
                flg2 = sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][k]][4] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][k]][3];
                flg1 = flg1 || flg2;
                flg2 = sp3c[mem_sp3c[sp3c[i][j2]][j]][3] == sp3c[mem_sp3c[sp3c[i][j2]][k]][3] || sp3c[mem_sp3c[sp3c[i][j2]][j]][4] == sp3c[mem_sp3c[sp3c[i][j2]][k]][4];
                flg3 = sp3c[mem_sp3c[sp3c[i][j2]][j]][3] == sp3c[mem_sp3c[sp3c[i][j2]][k]][4] || sp3c[mem_sp3c[sp3c[i][j2]][j]][4] == sp3c[mem_sp3c[sp3c[i][j2]][k]][3];
                flg2 = flg2 || flg3;
                flg1 = flg1 || flg2;
                if (flg1==1) continue; // common spindle particles

                n = 0;
                for (l=0; l<3; ++l) {
                    for (m=0; m<3; ++m) {
                        if (sp3c[i][l] == sp3c[mem_sp3c[sp3c[i][j2]][k]][m]) {
                            if (cp != sp3c[i][l]) break;
                            ++n;
                        }
                        if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == sp3c[mem_sp3c[sp3c[i][j2]][k]][m]) {
                            if(cp != sp3c[mem_sp3c[sp3c[i][j2]][j]][l]) break;
                            ++n;
                        }
                    }
                    if (m<3) break;
                }
                if (l<3) continue;
                if (n != 2) continue; // only one common particle, cp

                for (l=0; l<3; ++l) { // sp3 particles can't be spindle particles
                    x = sp3c[i][l];
                    flg1 = x == sp3c[mem_sp3c[sp3c[i][j2]][j]][3] || x == sp3c[mem_sp3c[sp3c[i][j2]][j]][4] || x == sp3c[mem_sp3c[sp3c[i][j2]][k]][3] || x == sp3c[mem_sp3c[sp3c[i][j2]][k]][4];
                    x = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                    flg2 = x == sp3c[i][3] || x == sp3c[i][4] || x == sp3c[mem_sp3c[sp3c[i][j2]][k]][3] || x == sp3c[mem_sp3c[sp3c[i][j2]][k]][4];
                    x = sp3c[mem_sp3c[sp3c[i][j2]][k]][l];
                    flg3 = x == sp3c[i][3] || x == sp3c[i][4] || x == sp3c[mem_sp3c[sp3c[i][j2]][j]][3] || x == sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                    if (flg1==1 || flg2==1 || flg3==1) break;
                }
                if (l<3) continue;

                // Check for bonds between sp3 particles and other spindles
                for (l=0; l<3; ++l) {
                    if (sp3c[i][l] != cp) {
                        flg1 = Bonds_BondCheck(sp3c[i][l], sp3c[mem_sp3c[sp3c[i][j2]][j]][3]) || Bonds_BondCheck(sp3c[i][l], sp3c[mem_sp3c[sp3c[i][j2]][j]][4]);
                        flg1 = flg1 || Bonds_BondCheck(sp3c[i][l], sp3c[mem_sp3c[sp3c[i][j2]][k]][3]) || Bonds_BondCheck(sp3c[i][l], sp3c[mem_sp3c[sp3c[i][j2]][k]][4]);
                        if (flg1==1) break;
                    }
                    if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l] != cp) {
                        flg1 = Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][j]][l], sp3c[i][3]) || Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][j]][l], sp3c[i][4]);
                        flg1 = flg1 || Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][j]][l], sp3c[mem_sp3c[sp3c[i][j2]][k]][3]) || Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][j]][l], sp3c[mem_sp3c[sp3c[i][j2]][k]][4]);
                        if (flg1==1) break;
                    }
                    if (sp3c[mem_sp3c[sp3c[i][j2]][k]][l] != cp) {
                        flg1 = Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][k]][l], sp3c[i][3]) || Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][k]][l], sp3c[i][4]);
                        flg1 = flg1 || Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][k]][l], sp3c[mem_sp3c[sp3c[i][j2]][j]][3]) || Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][k]][l], sp3c[mem_sp3c[sp3c[i][j2]][j]][4]);
                        if (flg1==1) break;
                    }
                }
                if (l<3) continue;

                m = 0;  // find uncommon particles from 5A_i
                for (l=0; l<3; ++l) {
                    if (sp3c[i][l] != cp) {
                        ia[m] = sp3c[i][l];
                        m++;
                    }
                }
                flg1 = flg2 = 0;
                ja[0]=ja[1]=-1; // find uncommon particles from 5A_j
                for (l=0; l<3; ++l) {
                    if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == cp) continue;
                    if (Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][j]][l], ia[0])==1) {
                        flg1 = 1;
                        ja[0] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                        if(Bonds_BondCheck(ja[0], ia[1])==1) break;
                    }
                    else ja[1] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                }
                if (l<3) continue;
                if (!flg1) {
                    m = ia[0];
                    ia[0] = ia[1];
                    ia[1] = m;
                    for (l=0; l<3; ++l) {
                        if (sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == cp) continue;
                        if (Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][j]][l],ia[0])==1) {
                            flg2 = 1;
                            ja[0] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                            if (Bonds_BondCheck(ja[0],ia[1])==1) break;
                        }
                        else ja[1] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                    }
                    if (l<3 || flg2==0) continue;
                }
                if (!(flg1==1 || flg2==1)) continue; // found four uncommon particles from 5A_i and 5A_j and found the two that are bonded between these
                if (ja[0]==-1 || ja[1]==-1) continue;
                if (Bonds_BondCheck(ja[1],ia[0])==1 || Bonds_BondCheck(ja[1],ia[1])==1) continue;
                flg1 = 0;
                ka[0]=ka[1]=-1;
                for (l=0; l<3; ++l) {
                    if (sp3c[mem_sp3c[sp3c[i][j2]][k]][l] == cp) continue;
                    if (Bonds_BondCheck(sp3c[mem_sp3c[sp3c[i][j2]][k]][l],ia[1])) {
                        flg1 = 1;
                        ka[0] = sp3c[mem_sp3c[sp3c[i][j2]][k]][l];
                        if (Bonds_BondCheck(ka[0],ia[0])==1 || Bonds_BondCheck(ka[0],ja[0])==1 || Bonds_BondCheck(ka[0],ja[1])==1) break;
                    }
                    else ka[1] = sp3c[mem_sp3c[sp3c[i][j2]][k]][l];
                }
                if (l<3 || ka[0]==-1 || ka[1]==-1) continue;
                if (Bonds_BondCheck(ka[1],ja[1])==0) continue;
                if (Bonds_BondCheck(ka[1],ia[1])==1 || Bonds_BondCheck(ka[1],ia[0])==1 || Bonds_BondCheck(ka[1],ja[0])==1) continue;

                // We've now found an HCP cluster
                if (nHCP[f] == mHCP) {
                    hcHCP=resize_2D_int(hcHCP,mHCP,mHCP+incrStatic,clusSize,-1);
                    mHCP=mHCP+incrStatic;
                }

                hcHCP[nHCP[f]][0] = cp;
                l = 1;
                for (m=0; m<3; ++m) {
                    if(sp3c[i][m] != cp) hcHCP[nHCP[f]][l++] = sp3c[i][m];
                    if(sp3c[mem_sp3c[sp3c[i][j2]][j]][m] != cp) hcHCP[nHCP[f]][l++] = sp3c[mem_sp3c[sp3c[i][j2]][j]][m];
                    if(sp3c[mem_sp3c[sp3c[i][j2]][k]][m] != cp) hcHCP[nHCP[f]][l++] = sp3c[mem_sp3c[sp3c[i][j2]][k]][m];
                }
                hcHCP[nHCP[f]][l++] = sp3c[i][3];
                hcHCP[nHCP[f]][l++] = sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                hcHCP[nHCP[f]][l++] = sp3c[mem_sp3c[sp3c[i][j2]][k]][3];
                hcHCP[nHCP[f]][l++] = sp3c[i][4];
                hcHCP[nHCP[f]][l++] = sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                hcHCP[nHCP[f]][l] = sp3c[mem_sp3c[sp3c[i][j2]][k]][4];
                quickSort(&hcHCP[nHCP[f]][1],6);
                quickSort(&hcHCP[nHCP[f]][7],6);

                Cluster_Write_HCP(f, i, j, j2, k, ach, ach_cen, ach_shell);
            }
        }
        }
    }
    for(i=0; i<N; ++i) {
        sHCP[i]=ach[i];
        sHCP_cen[i]=ach_cen[i];
        sHCP_shell[i]=ach_shell[i];
    }
    free(ach);
    free(ach_cen);
    free(ach_shell);
}

void Cluster_Write_HCP(int f, int i, int j, int j2, int k, char *ach, char *ach_cen, char *ach_shell) {
    if(ach[sp3c[i][0]] == 'C') ach[sp3c[i][0]] = 'B';
    if(ach[sp3c[i][1]] == 'C') ach[sp3c[i][1]] = 'B';
    if(ach[sp3c[i][2]] == 'C') ach[sp3c[i][2]] = 'B';
    if(ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][0]] == 'C') ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][0]] = 'B';
    if(ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][1]] == 'C') ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][1]] = 'B';
    if(ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][2]] == 'C') ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][2]] = 'B';
    if(ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][0]] == 'C') ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][0]] = 'B';
    if(ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][1]] == 'C') ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][1]] = 'B';
    if(ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][2]] == 'C') ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][2]] = 'B';
    ach[sp3c[i][3]] = 'O';
    ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][3]] = 'F';
    ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][3]] = 'H';
    ach[sp3c[i][4]] = 'O';
    ach[sp3c[mem_sp3c[sp3c[i][j2]][j]][4]] = 'F';
    ach[sp3c[mem_sp3c[sp3c[i][j2]][k]][4]] = 'H';
    ach_cen[hcHCP[nHCP[f]][0]] = 'F';
    ach_shell[hcHCP[nHCP[f]][1]] = 'B';
    ach_shell[hcHCP[nHCP[f]][2]] = 'B';
    ach_shell[hcHCP[nHCP[f]][3]] = 'B';
    ach_shell[hcHCP[nHCP[f]][4]] = 'B';
    ach_shell[hcHCP[nHCP[f]][5]] = 'B';
    ach_shell[hcHCP[nHCP[f]][6]] = 'B';
    ach_shell[hcHCP[nHCP[f]][7]] = 'B';
    ach_shell[hcHCP[nHCP[f]][8]] = 'B';
    ach_shell[hcHCP[nHCP[f]][9]] = 'B';
    ach_shell[hcHCP[nHCP[f]][10]] = 'B';
    ach_shell[hcHCP[nHCP[f]][11]] = 'B';
    ach_shell[hcHCP[nHCP[f]][12]] = 'B';

    ++nHCP[f];
}

void Clusters_GetBCC_9(int f) {
    int i, j, j2, k, l, m;
    int flg;
    int s_com=-1;
    int trial[9];
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int clusSize=9;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_GetBCC_9(): ach[] malloc out of memory\n");   Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_GetBCC_9(): ach_cen[] malloc out of memory\n");   Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_GetBCC_9(): ach_shell[] malloc out of memory\n");   Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for (i=0; i<nsp4b[f]-1; i++) {
        for (j2=4; j2<5; j2++) {
            for (j=0; j<nmem_sp4b[sp4b[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4b[sp4b[i][j2]][j]<=i) continue;
                if (sp4b[i][4]!=sp4b[mem_sp4b[sp4b[i][j2]][j]][4]) continue;
                s_com=sp4b[i][4];

                flg=0;
                for (k=0; k<4; k++) {
                    for (l=0; l<4; l++) {
                        if (sp4b[i][k]==sp4b[mem_sp4b[sp4b[i][j2]][j]][l]) {
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
                        if (Bonds_BondCheck(sp4b[i][k],sp4b[mem_sp4b[sp4b[i][j2]][j]][l])) {
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
                        if (Bonds_BondCheck(sp4b[mem_sp4b[sp4b[i][j2]][j]][k],sp4b[i][l])) {
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
                    trial[k+1]=sp4b[i][k];
                    trial[k+5]=sp4b[mem_sp4b[sp4b[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k<nBCC_9[f]; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l]!=hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9[f] == mBCC_9) {
                    hcBCC_9=resize_2D_int(hcBCC_9,mBCC_9,mBCC_9+incrStatic,clusSize,-1);
                    mBCC_9=mBCC_9+incrStatic;
                }

                for (k=0; k<9; ++k) hcBCC_9[nBCC_9[f]][k]=trial[k];

                Cluster_Write_BCC9(f, ach, ach_cen, ach_shell);
            }
        }
    }

    for (i=0; i<nsp4c[f]-1; i++) {
        for (j2=4; j2<6; j2++) {
            for (j=0; j<nmem_sp4c[sp4c[i][j2]]; ++j) { // loop over all sp3c_j
                if (mem_sp4c[sp4c[i][j2]][j]<=i) continue;
                m=0;
                for (k=4; k<6; k++) {
                    for (l=4; l<6; l++) {
                        if (sp4c[i][k]==sp4c[mem_sp4c[sp4c[i][j2]][j]][l]) {
                            s_com=sp4c[i][k];
                            m++;
                        }
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<6; k++) {
                    if (sp4c[i][k]==s_com) continue;
                    for (l=0; l<6; l++) {
                        if (sp4c[mem_sp4c[sp4c[i][j2]][j]][l]==s_com) continue;
                        if (sp4c[i][k]==sp4c[mem_sp4c[sp4c[i][j2]][j]][l]) {
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
                        if (Bonds_BondCheck(sp4c[i][k],sp4c[mem_sp4c[sp4c[i][j2]][j]][l])) {
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
                        if (Bonds_BondCheck(sp4c[mem_sp4c[sp4c[i][j2]][j]][k],sp4c[i][l])) {
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
                    trial[k+1]=sp4c[i][k];
                    trial[k+5]=sp4c[mem_sp4c[sp4c[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k<nBCC_9[f]; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l]!=hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9[f] == mBCC_9) {
                    hcBCC_9=resize_2D_int(hcBCC_9,mBCC_9,mBCC_9+incrStatic,clusSize,-1);
                    mBCC_9=mBCC_9+incrStatic;
                }
                for (k=0; k<9; ++k) hcBCC_9[nBCC_9[f]][k]=trial[k];

                Cluster_Write_BCC9(f, ach, ach_cen, ach_shell);
            }
        }
    }

    for (i=0; i<nsp4b[f]; i++) {
        for (j2=4; j2<5; j2++) {
            for (j=0; j<nmem_sp4c[sp4b[i][j2]]; ++j) { // loop over all sp3c_j
                m=0;
                for (k=4; k<6; k++) {
                    if (sp4b[i][4]==sp4c[mem_sp4c[sp4b[i][j2]][j]][k]) {
                        s_com=sp4b[i][4];
                        m++;
                    }
                }
                if (m==0 || m>1) continue;

                flg=0;
                for (k=0; k<4; k++) {
                    for (l=0; l<6; l++) {
                        if (sp4c[mem_sp4c[sp4b[i][j2]][j]][l]==s_com) continue;
                        if (sp4b[i][k]==sp4c[mem_sp4c[sp4b[i][j2]][j]][l]) {
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
                        if (Bonds_BondCheck(sp4b[i][k],sp4c[mem_sp4c[sp4b[i][j2]][j]][l])) {
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
                        if (Bonds_BondCheck(sp4c[mem_sp4c[sp4b[i][j2]][j]][k],sp4b[i][l])) {
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
                    trial[k+1]=sp4b[i][k];
                    trial[k+5]=sp4c[mem_sp4c[sp4b[i][j2]][j]][k];
                }
                quickSort(&trial[1],8);

                flg=0;  // check trial cluster not already found
                for (k=0; k<nBCC_9[f]; ++k) {
                    for (l=0; l<9; ++l) {
                        if (trial[l]!=hcBCC_9[k][l]) break;
                    }
                    if (l==9) flg=1;
                }
                if (flg==1) continue;

                if (nBCC_9[f] == mBCC_9) {
                    hcBCC_9=resize_2D_int(hcBCC_9,mBCC_9,mBCC_9+incrStatic,clusSize,-1);
                    mBCC_9=mBCC_9+incrStatic;
                }
                for (k=0; k<9; ++k) hcBCC_9[nBCC_9[f]][k]=trial[k];

                Cluster_Write_BCC9(f, ach, ach_cen, ach_shell);

            }
        }
    }

    for(i=0; i<N; ++i) {
        sBCC_9[i]=ach[i];
        sBCC_9_cen[i]=ach_cen[i];
        sBCC_9_shell[i]=ach_shell[i];
    }
    free(ach);
    free(ach_cen);
    free(ach_shell);
}

void Cluster_Write_BCC9(int f, char *ach, char *ach_cen, char *ach_shell) {
    if(ach[hcBCC_9[nBCC_9[f]][1]] == 'C') ach[hcBCC_9[nBCC_9[f]][1]] = ach_shell[hcBCC_9[nBCC_9[f]][1]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][2]] == 'C') ach[hcBCC_9[nBCC_9[f]][2]] = ach_shell[hcBCC_9[nBCC_9[f]][2]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][3]] == 'C') ach[hcBCC_9[nBCC_9[f]][3]] = ach_shell[hcBCC_9[nBCC_9[f]][3]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][4]] == 'C') ach[hcBCC_9[nBCC_9[f]][4]] = ach_shell[hcBCC_9[nBCC_9[f]][4]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][5]] == 'C') ach[hcBCC_9[nBCC_9[f]][5]] = ach_shell[hcBCC_9[nBCC_9[f]][5]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][6]] == 'C') ach[hcBCC_9[nBCC_9[f]][6]] = ach_shell[hcBCC_9[nBCC_9[f]][6]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][7]] == 'C') ach[hcBCC_9[nBCC_9[f]][7]] = ach_shell[hcBCC_9[nBCC_9[f]][7]] = 'B';
    if(ach[hcBCC_9[nBCC_9[f]][8]] == 'C') ach[hcBCC_9[nBCC_9[f]][8]] = ach_shell[hcBCC_9[nBCC_9[f]][8]] = 'B';
    ach[hcBCC_9[nBCC_9[f]][0]] = ach_cen[hcBCC_9[nBCC_9[f]][0]] = 'F';

    ++nBCC_9[f];
}

void Clusters_GetBCC_15(int f) {    // Detect 15 particle BCC clusters
    int i,j,k,l,m;
    int no_sp4cs,noSP4s;
    int sj[5];
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int clusSize=15;

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_GetBCC_15(): ach[] malloc out of memory\n");  Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_GetBCC_15(): ach_cen[] malloc out of memory\n");  Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_GetBCC_15(): ach_shell[] malloc out of memory\n");  Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for (i=0; i<nsp4c[f]; i++) {
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15[f] == mBCC_15) {
            hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
            mBCC_15=mBCC_15+incrStatic;
        }
        for (j=0; j<nBCC_15[f]; j++) if (sp4c[i][4]==hcBCC_15[j][0]) break;
        if (j==nBCC_15[f]) {
            hcBCC_15[nBCC_15[f]][0]=sp4c[i][4];
            hcBCC_15[nBCC_15[f]][1]=sp4c[i][5];
            for (j=0; j<4; j++) hcBCC_15[nBCC_15[f]][j+7]=sp4c[i][j];

            no_sp4cs=0;
            noSP4s=4;
            sj[0]=sj[1]=sj[2]=sj[3]=sj[4]=-1;
            for (j=i+1; j<nsp4c[f]; j++) {
                if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][0] || sp4c[j][5]==hcBCC_15[nBCC_15[f]][0]) {
                    m=0;
                    for (l=1;l<2+no_sp4cs;l++) {
                        if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][l] || sp4c[j][5]==hcBCC_15[nBCC_15[f]][l]) m++;
                    }
                    if (m==0) {
                        for (l=7;l<7+noSP4s;l++) {
                            if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][l] || sp4c[j][5]==hcBCC_15[nBCC_15[f]][l]) m++;
                        }
                        if (m==0) {
                            for (k=0; k<4; k++) {
                                for (l=0;l<2+no_sp4cs;l++) {
                                    if (sp4c[j][k]==hcBCC_15[nBCC_15[f]][l]) m++;
                                }
                            }
                            if (m==0) {
                                if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][0]) hcBCC_15[nBCC_15[f]][2+no_sp4cs]=sp4c[j][5];
                                else hcBCC_15[nBCC_15[f]][2+no_sp4cs]=sp4c[j][4];
                                for (k=0; k<4; k++) {
                                    m=0;
                                    for (l=7;l<7+noSP4s;l++) {
                                        if (sp4c[j][k]==hcBCC_15[nBCC_15[f]][l]) m++;
                                    }
                                    if (m==0) {
                                        if (noSP4s>=8) {
                                            noSP4s++;
                                            break;
                                        }
                                        hcBCC_15[nBCC_15[f]][7+noSP4s]=sp4c[j][k];
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
                Cluster_Write_BCC_15(f, ach, ach_cen, ach_shell, clusSize);
            }
        }
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15[f] == mBCC_15) {
            hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
            mBCC_15=mBCC_15+incrStatic;
        }
        for (j=0; j<nBCC_15[f]; j++) if (sp4c[i][5]==hcBCC_15[j][0]) break;
        if (j<nBCC_15[f]) continue;

        hcBCC_15[nBCC_15[f]][0]=sp4c[i][5];
        hcBCC_15[nBCC_15[f]][1]=sp4c[i][4];
        for (j=0; j<4; j++) hcBCC_15[nBCC_15[f]][j+7]=sp4c[i][j];

        no_sp4cs=0;
        noSP4s=4;
        sj[0]=sj[1]=sj[2]=sj[3]=sj[4]=-1;
        for (j=i+1; j<nsp4c[f]; j++) {
            if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][0] || sp4c[j][5]==hcBCC_15[nBCC_15[f]][0]) {
                m=0;
                for (l=1;l<2+no_sp4cs;l++) {
                    if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][l] || sp4c[j][5]==hcBCC_15[nBCC_15[f]][l]) m++;
                }
                if (m==0) {
                    for (l=7;l<7+noSP4s;l++) {
                        if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][l] || sp4c[j][5]==hcBCC_15[nBCC_15[f]][l]) m++;
                    }
                    if (m==0) {
                        for (k=0; k<4; k++) {
                            for (l=0;l<2+no_sp4cs;l++) {
                                if (sp4c[j][k]==hcBCC_15[nBCC_15[f]][l]) m++;
                            }
                        }
                        if (m==0) {
                            if (sp4c[j][4]==hcBCC_15[nBCC_15[f]][0]) hcBCC_15[nBCC_15[f]][2+no_sp4cs]=sp4c[j][5];
                            else hcBCC_15[nBCC_15[f]][2+no_sp4cs]=sp4c[j][4];
                            for (k=0; k<4; k++) {
                                m=0;
                                for (l=7;l<7+noSP4s;l++) {
                                    if (sp4c[j][k]==hcBCC_15[nBCC_15[f]][l]) m++;
                                }
                                if (m==0) {
                                    if (noSP4s>=8) {
                                        noSP4s++;
                                        break;
                                    }
                                    hcBCC_15[nBCC_15[f]][7+noSP4s]=sp4c[j][k];
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
            Cluster_Write_BCC_15(f, ach, ach_cen, ach_shell, clusSize);
        }

    }
    for(i=0; i<N; ++i) {
        sBCC_15[i]=ach[i];
        sBCC_15_cen[i]=ach_cen[i];
        sBCC_15_shell[i]=ach_shell[i];
    }
    free(ach);
    free(ach_cen);
    free(ach_shell);
}

void Cluster_Write_BCC_15(int f, char *ach, char *ach_cen, char *ach_shell, int clusSize) {
    if (nBCC_15[f] == mBCC_15) {
                hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
                mBCC_15=mBCC_15+incrStatic;
            }

    quickSort(&hcBCC_15[nBCC_15[f]][1],6);
    quickSort(&hcBCC_15[nBCC_15[f]][7],8);

    if(ach[hcBCC_15[nBCC_15[f]][7]] == 'C') ach[hcBCC_15[nBCC_15[f]][7]] = ach_shell[hcBCC_15[nBCC_15[f]][7]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][8]] == 'C') ach[hcBCC_15[nBCC_15[f]][8]] = ach_shell[hcBCC_15[nBCC_15[f]][8]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][9]] == 'C') ach[hcBCC_15[nBCC_15[f]][9]] = ach_shell[hcBCC_15[nBCC_15[f]][9]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][10]] == 'C') ach[hcBCC_15[nBCC_15[f]][10]] = ach_shell[hcBCC_15[nBCC_15[f]][10]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][11]] == 'C') ach[hcBCC_15[nBCC_15[f]][11]] = ach_shell[hcBCC_15[nBCC_15[f]][11]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][12]] == 'C') ach[hcBCC_15[nBCC_15[f]][12]] = ach_shell[hcBCC_15[nBCC_15[f]][12]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][13]] == 'C') ach[hcBCC_15[nBCC_15[f]][13]] = ach_shell[hcBCC_15[nBCC_15[f]][13]] = 'B';
    if(ach[hcBCC_15[nBCC_15[f]][14]] == 'C') ach[hcBCC_15[nBCC_15[f]][14]] = ach_shell[hcBCC_15[nBCC_15[f]][14]] = 'B';
    ach[hcBCC_15[nBCC_15[f]][1]] = ach_shell[hcBCC_15[nBCC_15[f]][1]] = 'O';
    ach[hcBCC_15[nBCC_15[f]][2]] = ach_shell[hcBCC_15[nBCC_15[f]][2]] = 'O';
    ach[hcBCC_15[nBCC_15[f]][3]] = ach_shell[hcBCC_15[nBCC_15[f]][3]] = 'O';
    ach[hcBCC_15[nBCC_15[f]][4]] = ach_shell[hcBCC_15[nBCC_15[f]][4]] = 'O';
    ach[hcBCC_15[nBCC_15[f]][5]] = ach_shell[hcBCC_15[nBCC_15[f]][5]] = 'O';
    ach[hcBCC_15[nBCC_15[f]][6]] = ach_shell[hcBCC_15[nBCC_15[f]][6]] = 'O';
    ach[hcBCC_15[nBCC_15[f]][0]] = ach_cen[hcBCC_15[nBCC_15[f]][0]] = 'F';

    ++nBCC_15[f];
}