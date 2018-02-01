#include "globals.h"
#include "clusters.h"
#include "bonds.h"

void Clusters_Get6Z_C2v(int f) {    // Detect 6Z clusters from 2 5A clusters
    int flg;
    int i, j, j2, k, l;
    int cnt;
    int s1a, s2a, s1b, s2b;
    char *ach, errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=6;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_6Z=resize_1D_double(bl_mom_6Z,m6Z,m6Z+incrStatic);
                }
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
            // hc6Z key: (5A_i_s_in_SP3_j, 5A_j_s_in_SP3_i, 5A_i_s_oth, 5A_j_s_oth, SP3_i_j_com_1, SP3_i_j_com_2)
            if (doDynamics==1 && dyn_m6Z!=-1) {
                if (doSubClusts==1 && dyn_msp3c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp3c[i];
                    sub[1]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][j]];
                    quickSort(&sub[0],2);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc6Z[n6Z[f]], f, clusSize, &dyn_n6Z, &dyn_m6Z, &dyn_l6Z, &dyn_hc6Z, do_up, dummy_up, n6Z[f], do_sub, n_sub, &dyn_sub_6Z, sub);
            }
            ach[hc6Z[n6Z[f]][0]] = 'O';
            ach[hc6Z[n6Z[f]][1]] = 'O';
            ach[hc6Z[n6Z[f]][2]] = 'O';
            ach[hc6Z[n6Z[f]][3]] = 'O';
            if (ach[hc6Z[n6Z[f]][4]] == 'C') ach[hc6Z[n6Z[f]][4]] = 'B';
            if (ach[hc6Z[n6Z[f]][5]] == 'C') ach[hc6Z[n6Z[f]][5]] = 'B';
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc6Z[n6Z[f]][binAcnt],hc6Z[n6Z[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc6Z[n6Z[f]][binAcnt]][Bonds_cnb_j(hc6Z[n6Z[f]][binAcnt],hc6Z[n6Z[f]][binBcnt])],BLDistro6Z,&BLDistroNoSamples6Z);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc6Z[n6Z[f]][binAcnt]]==1) {
                        nA6Z++;
                        number_of_A++;
                    }
                    else nB6Z++;
                }
                n_distro_6Z[number_of_A]++;
            }
            
            ++n6Z[f]; 
        }
        }
    }
    for (i=0; i<N; ++i) s6Z[i]=ach[i];
    free(ach);
}

void Clusters_Get7K(int f) {    // Detect 7K clusters from 2 5A clusters
    int i, j, j2, k, l, m;
    int scom, sother[2], sp3_com[2], sp3c_i_other, sp3c_j_other;
    char *ach, errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=7;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_7K=resize_1D_double(bl_mom_7K,m7K,m7K+incrStatic);
                }
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
            
            // hc7K key: (scom, sother, ring_com, ring_other)
            
            if (doDynamics==1 && dyn_m7K!=-1) {
                if (doSubClusts==1 && dyn_msp3c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp3c[i];
                    sub[1]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][j]];
                    quickSort(&sub[0],2);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc7K[n7K[f]], f, clusSize, &dyn_n7K, &dyn_m7K, &dyn_l7K, &dyn_hc7K, do_up, dummy_up, n7K[f], do_sub, n_sub, &dyn_sub_7K, sub);
            }
            ach[hc7K[n7K[f]][0]] = 'O';
            ach[hc7K[n7K[f]][1]] = 'O';
            ach[hc7K[n7K[f]][2]] = 'O';
            if (ach[hc7K[n7K[f]][3]] == 'C') ach[hc7K[n7K[f]][3]] = 'B';
            if (ach[hc7K[n7K[f]][4]] == 'C') ach[hc7K[n7K[f]][4]] = 'B';
            if (ach[hc7K[n7K[f]][5]] == 'C') ach[hc7K[n7K[f]][5]] = 'B';
            if (ach[hc7K[n7K[f]][6]] == 'C') ach[hc7K[n7K[f]][6]] = 'B';
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc7K[n7K[f]][binAcnt],hc7K[n7K[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc7K[n7K[f]][binAcnt]][Bonds_cnb_j(hc7K[n7K[f]][binAcnt],hc7K[n7K[f]][binBcnt])],BLDistro7K,&BLDistroNoSamples7K);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc7K[n7K[f]][binAcnt]]==1) {
                        nA7K++;
                        number_of_A++;
                    }
                    else nB7K++;
                }
                n_distro_7K[number_of_A]++;
            }
            
            ++n7K[f]; 
        }
        }
    }
    for (i=0; i<N; ++i) s7K[i]=ach[i];
    free(ach);
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=8;
    int *used_sp5b;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=12;
    int sub[12];
    
    unc[0]=unc[1]=com[0]=com[1]=com[2]=com[3]=-1;
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get8A_D2d(): ach[] malloc out of memory\n");  Error(errMsg); }
    for (i=0; i<N; ++i) ach[i] = 'C';
    
    used_sp5b=malloc(nsp5b[f]*sizeof(int)); if (used_sp5b==NULL) { sprintf(errMsg,"Clusters_Get8A_D2d(): used_sp5b[] malloc out of memory\n");  Error(errMsg); }
    for (i=0; i<nsp5b[f]; ++i) used_sp5b[i] = 0;
    
    for (i=0; i<nsp5b[f]-1; ++i) {  // loop over all sp5b_i
        for (j2=0; j2<nsp5b[f]; ++j2) used_sp5b[j2] = 0;
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
                if (doClusBLDeviation==1) {
                    bl_mom_8A=resize_1D_double(bl_mom_8A,m8A,m8A+incrStatic);
                }
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
                
                if (doDynamics==1 && dyn_m8A!=-1) {
                    if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
                        sub[0]=dyn_up_sp5b[i];
                        sub[1]=dyn_up_sp5b[mem_sp5c[sp5b[i][j2]][j]];
                        sub[2]=-1;
                        sub[3]=-1;
                        sub[4]=-1;
                        sub[5]=-1;
                        sub[6]=-1;
                        sub[7]=-1;
                        sub[8]=-1;
                        sub[9]=-1;
                        sub[10]=-1;
                        sub[11]=-1;
                        quickSort(&sub[0],2);
                        do_sub=1;
                    }
                    else {
                        sub[0]=1;
                        sub[1]=1;
                        sub[2]=-1;
                        sub[3]=-1;
                        sub[4]=-1;
                        sub[5]=-1;
                        sub[6]=-1;
                        sub[7]=-1;
                        sub[8]=-1;
                        sub[9]=-1;
                        sub[10]=-1;
                        sub[11]=-1;
                        do_sub=0;
                    }
                    do_up=0;
                    Dyn_add_8A(trial, f, clusSize, &dyn_n8A, &dyn_m8A, &dyn_l8A, &dyn_hc8A, do_up, dummy_up, n8A[f], do_sub, n_sub, &dyn_sub_8A, sub);
                }
                if (ach[hc8A[n8A[f]][0]] == 'C') ach[hc8A[n8A[f]][0]] = 'B';
                if (ach[hc8A[n8A[f]][1]] == 'C') ach[hc8A[n8A[f]][1]] = 'B';
                if (ach[hc8A[n8A[f]][2]] == 'C') ach[hc8A[n8A[f]][2]] = 'B';
                if (ach[hc8A[n8A[f]][3]] == 'C') ach[hc8A[n8A[f]][3]] = 'B';
                if (ach[hc8A[n8A[f]][4]] == 'C') ach[hc8A[n8A[f]][4]] = 'B';
                if (ach[hc8A[n8A[f]][5]] == 'C') ach[hc8A[n8A[f]][5]] = 'B';
                ach[hc8A[n8A[f]][6]] = 'O';
                ach[hc8A[n8A[f]][7]] = 'O';
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc8A[n8A[f]][binAcnt],hc8A[n8A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc8A[n8A[f]][binAcnt]][Bonds_cnb_j(hc8A[n8A[f]][binAcnt],hc8A[n8A[f]][binBcnt])],BLDistro8A,&BLDistroNoSamples8A);
                            }
                        }
                    }
                }

                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc8A[n8A[f]][binAcnt]]==1) {
                            nA8A++;
                            number_of_A++;
                        }
                        else nB8A++;
                    }
                    n_distro_8A[number_of_A]++;
                }
                
                ++n8A[f];
                
                // hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
            }
            else if (doDynamics==1 && dyn_m8A!=-1) {
                if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    sub[6]=dyn_up_sp5b[i];
                    sub[7]=dyn_up_sp5b[mem_sp5c[sp5b[i][j2]][j]];
                    sub[8]=-1;
                    sub[9]=-1;
                    sub[10]=-1;
                    sub[11]=-1;
                    quickSort(&sub[6],2);
                    do_sub=1;
                }
                else {
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    sub[6]=1;
                    sub[7]=1;
                    sub[8]=-1;
                    sub[9]=-1;
                    sub[10]=-1;
                    sub[11]=-1;
                    do_sub=0;
                }
                do_up=0;
                Dyn_add_8A(trial, f, clusSize, &dyn_n8A, &dyn_m8A, &dyn_l8A, &dyn_hc8A, do_up, dummy_up, n8A[f], do_sub, n_sub, &dyn_sub_8A, sub);
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
                if (doClusBLDeviation==1) {
                    bl_mom_8A=resize_1D_double(bl_mom_8A,m8A,m8A+incrStatic);
                }
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
                
                if (doDynamics==1 && dyn_m8A!=-1) {
                    if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
                        sub[0]=-1;
                        sub[1]=-1;
                        sub[2]=dyn_up_sp5c[i];
                        sub[3]=dyn_up_sp5c[mem_sp5c[sp5b[i][j2]][j]];
                        sub[4]=-1;
                        sub[5]=-1;
                        sub[6]=-1;
                        sub[7]=-1;
                        sub[8]=-1;
                        sub[9]=-1;
                        sub[10]=-1;
                        sub[11]=-1;
                        quickSort(&sub[2],2);
                        do_sub=1;
                    }
                    else {
                        sub[0]=-1;
                        sub[1]=-1;
                        sub[2]=1;
                        sub[3]=1;
                        sub[4]=-1;
                        sub[5]=-1;
                        sub[6]=-1;
                        sub[7]=-1;
                        sub[8]=-1;
                        sub[9]=-1;
                        sub[10]=-1;
                        sub[11]=-1;
                        do_sub=0;
                    }
                    do_up=0;
                    Dyn_add_8A(trial, f, clusSize, &dyn_n8A, &dyn_m8A, &dyn_l8A, &dyn_hc8A, do_up, dummy_up, n8A[f], do_sub, n_sub, &dyn_sub_8A, sub);
                }
                if (ach[hc8A[n8A[f]][0]] == 'C') ach[hc8A[n8A[f]][0]] = 'B';
                if (ach[hc8A[n8A[f]][1]] == 'C') ach[hc8A[n8A[f]][1]] = 'B';
                if (ach[hc8A[n8A[f]][2]] == 'C') ach[hc8A[n8A[f]][2]] = 'B';
                if (ach[hc8A[n8A[f]][3]] == 'C') ach[hc8A[n8A[f]][3]] = 'B';
                if (ach[hc8A[n8A[f]][4]] == 'C') ach[hc8A[n8A[f]][4]] = 'B';
                if (ach[hc8A[n8A[f]][5]] == 'C') ach[hc8A[n8A[f]][5]] = 'B';
                ach[hc8A[n8A[f]][6]] = 'O';
                ach[hc8A[n8A[f]][7]] = 'O';
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc8A[n8A[f]][binAcnt],hc8A[n8A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc8A[n8A[f]][binAcnt]][Bonds_cnb_j(hc8A[n8A[f]][binAcnt],hc8A[n8A[f]][binBcnt])],BLDistro8A,&BLDistroNoSamples8A);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc8A[n8A[f]][binAcnt]]==1) {
                            nA8A++;
                            number_of_A++;
                        }
                        else nB8A++;
                    }
                    n_distro_8A[number_of_A]++;
                }
                
                ++n8A[f];
                
                // hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
            }
            else if (doDynamics==1 && dyn_m8A!=-1) {
                if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    sub[6]=-1;
                    sub[7]=-1;
                    sub[8]=dyn_up_sp5c[i];
                    sub[9]=dyn_up_sp5c[mem_sp5c[sp5b[i][j2]][j]];
                    sub[10]=-1;
                    sub[11]=-1;
                    quickSort(&sub[8],2);
                    do_sub=1;
                }
                else {
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    sub[6]=-1;
                    sub[7]=-1;
                    sub[8]=1;
                    sub[9]=1;
                    sub[10]=-1;
                    sub[11]=-1;
                    do_sub=0;
                }
                do_up=0;
                Dyn_add_8A(trial, f, clusSize, &dyn_n8A, &dyn_m8A, &dyn_l8A, &dyn_hc8A, do_up, dummy_up, n8A[f], do_sub, n_sub, &dyn_sub_8A, sub);
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
                if (doClusBLDeviation==1) {
                    bl_mom_8A=resize_1D_double(bl_mom_8A,m8A,m8A+incrStatic);
                }
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
                
                if (doDynamics==1 && dyn_m8A!=-1) {
                    if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
                        do_sub=1;
                        sub[0]=-1;
                        sub[1]=-1;
                        sub[2]=-1;
                        sub[3]=-1;
                        sub[4]=dyn_up_sp5b[i];
                        sub[5]=dyn_up_sp5c[mem_sp5c[sp5b[i][j2]][j]];
                        sub[6]=-1;
                        sub[7]=-1;
                        sub[8]=-1;
                        sub[9]=-1;
                        sub[10]=-1;
                        sub[11]=-1;
                    }
                    else {
                        do_sub=0;
                        sub[0]=-1;
                        sub[1]=-1;
                        sub[2]=-1;
                        sub[3]=-1;
                        sub[4]=1;
                        sub[5]=1;
                        sub[6]=-1;
                        sub[7]=-1;
                        sub[8]=-1;
                        sub[9]=-1;
                        sub[10]=-1;
                        sub[11]=-1;
                    }
                    do_up=0;
                    Dyn_add_8A(trial, f, clusSize, &dyn_n8A, &dyn_m8A, &dyn_l8A, &dyn_hc8A, do_up, dummy_up, n8A[f], do_sub, n_sub, &dyn_sub_8A, sub);
                }
                if (ach[hc8A[n8A[f]][0]] == 'C') ach[hc8A[n8A[f]][0]] = 'B';
                if (ach[hc8A[n8A[f]][1]] == 'C') ach[hc8A[n8A[f]][1]] = 'B';
                if (ach[hc8A[n8A[f]][2]] == 'C') ach[hc8A[n8A[f]][2]] = 'B';
                if (ach[hc8A[n8A[f]][3]] == 'C') ach[hc8A[n8A[f]][3]] = 'B';
                if (ach[hc8A[n8A[f]][4]] == 'C') ach[hc8A[n8A[f]][4]] = 'B';
                if (ach[hc8A[n8A[f]][5]] == 'C') ach[hc8A[n8A[f]][5]] = 'B';
                ach[hc8A[n8A[f]][6]] = 'O';
                ach[hc8A[n8A[f]][7]] = 'O';
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc8A[n8A[f]][binAcnt],hc8A[n8A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc8A[n8A[f]][binAcnt]][Bonds_cnb_j(hc8A[n8A[f]][binAcnt],hc8A[n8A[f]][binBcnt])],BLDistro8A,&BLDistroNoSamples8A);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc8A[n8A[f]][binAcnt]]==1) {
                            nA8A++;
                            number_of_A++;
                        }
                        else nB8A++;
                    }
                    n_distro_8A[number_of_A]++;
                }
                
                ++n8A[f];
                
                // hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
            }
            else if (doDynamics==1 && dyn_m8A!=-1) {
                if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    sub[6]=-1;
                    sub[7]=-1;
                    sub[8]=-1;
                    sub[9]=-1;
                    sub[10]=dyn_up_sp5b[i];
                    sub[11]=dyn_up_sp5c[mem_sp5c[sp5b[i][j2]][j]];
                }
                else {
                    do_sub=0;
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    sub[6]=-1;
                    sub[7]=-1;
                    sub[8]=-1;
                    sub[9]=-1;
                    sub[10]=1;
                    sub[11]=1;
                }
                do_up=0;
                Dyn_add_8A(trial, f, clusSize, &dyn_n8A, &dyn_m8A, &dyn_l8A, &dyn_hc8A, do_up, dummy_up, n8A[f], do_sub, n_sub, &dyn_sub_8A, sub);
            }
        }
        }
    }

    for(i=0; i<N; ++i) s8A[i]=ach[i];
    free(ach);
    free(used_sp5b);
}

void Clusters_Get8B_Cs(int f) { // Detect 8B Cs clusters
    int i, j, k, l, m;
    int s1, s2, n1, nbs, unc[3];
    int break_out;
    char *ach, errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=8;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=1;
    int sub[1];
    
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_Get8B_Cs(): ach[] malloc out of memory\n");   Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = 'C';
    
    for (i=0; i<nsp5c[f]; ++i) {    // loop over all 7A_i
        s1 = sp5c[i][5];
        s2 = sp5c[i][6];
        for (j=0; j<cnb[s1]; ++j) { // loop over all j particles bonded to first spindle of 7A_i
            n1 = bNums[s1][j];
            for (k=0; k<5; ++k) if (n1 == sp5c[i][k]) break;
            if (k<5) continue;
            if (n1 == s2) continue; // now is n1 bonded to sp5
            nbs = 0; // number of bonds
            for (k=0; k<5; ++k) if (Bonds_BondCheck(n1, sp5c[i][k])) ++nbs;
            if (nbs != 2) continue; 
            
            // Now we have found the 8B Cs cluster
            if (n8B[f]==m8B) { 
                hc8B=resize_2D_int(hc8B,m8B,m8B+incrStatic,clusSize,-1);
                if (doClusBLDeviation==1) {
                    bl_mom_8B=resize_1D_double(bl_mom_8B,m8B,m8B+incrStatic);
                }
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
            
            hc8B[n8B[f]][5]=sp5c[i][6]; // note switch of spindles here so hc8B[.][.][6] spindle is bonded to the extra particle hc8B[.][.][7]
            hc8B[n8B[f]][6]=sp5c[i][5]; // note switch of spindles here so hc8B[.][.][6] spindle is bonded to the extra particle hc8B[.][.][7]
            hc8B[n8B[f]][7]=n1;
            // hc8B key: (SP5_to_4, SP5_to_0/2, SP5_to_3, SP5_to_n1(lower), SP5_to_n1(greater), s, s_to_n1, n1)
            
            if (doDynamics==1 && dyn_m8B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc8B[n8B[f]], f, clusSize, &dyn_n8B, &dyn_m8B, &dyn_l8B, &dyn_hc8B, do_up, dummy_up, n8B[f], do_sub, n_sub, &dyn_sub_8B, sub);
            }
            if (ach[hc8B[n8B[f]][7]] == 'C') ach[hc8B[n8B[f]][7]] = 'B';
            if (ach[hc8B[n8B[f]][0]] == 'C') ach[hc8B[n8B[f]][0]] = 'B';
            if (ach[hc8B[n8B[f]][1]] == 'C') ach[hc8B[n8B[f]][1]] = 'B';
            if (ach[hc8B[n8B[f]][2]] == 'C') ach[hc8B[n8B[f]][2]] = 'B';
            if (ach[hc8B[n8B[f]][3]] == 'C') ach[hc8B[n8B[f]][3]] = 'B';
            if (ach[hc8B[n8B[f]][4]] == 'C') ach[hc8B[n8B[f]][4]] = 'B';
            ach[hc8B[n8B[f]][5]] = 'O';
            ach[hc8B[n8B[f]][6]] = 'O';
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc8B[n8B[f]][binAcnt],hc8B[n8B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc8B[n8B[f]][binAcnt]][Bonds_cnb_j(hc8B[n8B[f]][binAcnt],hc8B[n8B[f]][binBcnt])],BLDistro8B,&BLDistroNoSamples8B);
                        }
                    }
                }
            }

            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc8B[n8B[f]][binAcnt]]==1) {
                        nA8B++;
                        number_of_A++;
                    }
                    else nB8B++;
                }
                n_distro_8B[number_of_A]++;
            }
            
            ++n8B[f];
            
        }
        for (j=0; j<cnb[s2]; ++j) {
            n1 = bNums[s2][j];
            for (k=0; k<5; ++k) if (n1 == sp5c[i][k]) break;
            if (k<5) continue;
            if (n1 == s1) continue; // now is n1 bonded to sp5
            nbs = 0;
            for (k=0; k<5; ++k) if (Bonds_BondCheck(n1,sp5c[i][k])) ++nbs;
            if (nbs != 2) continue;
            
            // Now we have found the 8B Cs cluster
            if (n8B[f]==m8B) { 
                hc8B=resize_2D_int(hc8B,m8B,m8B+incrStatic,clusSize,-1);
                if (doClusBLDeviation==1) {
                    bl_mom_8B=resize_1D_double(bl_mom_8B,m8B,m8B+incrStatic);
                }
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

            hc8B[n8B[f]][5]=sp5c[i][5]; 
            hc8B[n8B[f]][6]=sp5c[i][6]; // no switch of spindles
            hc8B[n8B[f]][7]=n1;
            
            // hc8B key: (SP5_to_4, SP5_to_0/2, SP5_to_3, SP5_to_n1(lower), SP5_to_n1(greater), s, s_to_n1, n1)
            if (doDynamics==1 && dyn_m8B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc8B[n8B[f]], f, clusSize, &dyn_n8B, &dyn_m8B, &dyn_l8B, &dyn_hc8B, do_up, dummy_up, n8B[f], do_sub, n_sub, &dyn_sub_8B, sub);
            }
            if (ach[hc8B[n8B[f]][7]] == 'C') ach[hc8B[n8B[f]][7]] = 'B';
            if (ach[hc8B[n8B[f]][0]] == 'C') ach[hc8B[n8B[f]][0]] = 'B';
            if (ach[hc8B[n8B[f]][1]] == 'C') ach[hc8B[n8B[f]][1]] = 'B';
            if (ach[hc8B[n8B[f]][2]] == 'C') ach[hc8B[n8B[f]][2]] = 'B';
            if (ach[hc8B[n8B[f]][3]] == 'C') ach[hc8B[n8B[f]][3]] = 'B';
            if (ach[hc8B[n8B[f]][4]] == 'C') ach[hc8B[n8B[f]][4]] = 'B';
            ach[hc8B[n8B[f]][5]] = 'O';
            ach[hc8B[n8B[f]][6]] = 'O';
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc8B[n8B[f]][binAcnt],hc8B[n8B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc8B[n8B[f]][binAcnt]][Bonds_cnb_j(hc8B[n8B[f]][binAcnt],hc8B[n8B[f]][binBcnt])],BLDistro8B,&BLDistroNoSamples8B);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc8B[n8B[f]][binAcnt]]==1) {
                        nA8B++;
                        number_of_A++;
                    }
                    else nB8B++;
                }
                n_distro_8B[number_of_A]++;
            }
            
            ++n8B[f];
        }
    } 
    for(i=0; i<N; ++i) s8B[i]=ach[i];
    free(ach);
}

void Clusters_Get8K(int f) {    // Detect 8K clusters
    int i, j, j2, k, l, m, n;
    int cp[2], unc[3], scom, sother[2];
    int trial[8];
    char *ach, errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=8;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=3;
    int sub[3];
    
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
                    if (doClusBLDeviation==1) {
                        bl_mom_8K=resize_1D_double(bl_mom_8K,m8K,m8K+incrStatic);
                    }
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
                
                if (doDynamics==1 && dyn_m8K!=-1) {
                    if (doSubClusts==1 && dyn_msp3c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp3c[i];
                        sub[1]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][j]];
                        sub[2]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][k]];
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc8K[n8K[f]], f, clusSize, &dyn_n8K, &dyn_m8K, &dyn_l8K, &dyn_hc8K, do_up, dummy_up, n8K[f], do_sub, n_sub, &dyn_sub_8K, sub);
                }
                if(ach[hc8K[n8K[f]][2]] == 'C') ach[hc8K[n8K[f]][2]] = 'B';
                if(ach[hc8K[n8K[f]][3]] == 'C') ach[hc8K[n8K[f]][3]] = 'B';
                if(ach[hc8K[n8K[f]][4]] == 'C') ach[hc8K[n8K[f]][4]] = 'B';
                if(ach[hc8K[n8K[f]][5]] == 'C') ach[hc8K[n8K[f]][5]] = 'B';
                if(ach[hc8K[n8K[f]][6]] == 'C') ach[hc8K[n8K[f]][6]] = 'B';
                if(ach[hc8K[n8K[f]][7]] == 'C') ach[hc8K[n8K[f]][7]] = 'B';
                ach[hc8K[n8K[f]][0]] = 'O';
                ach[hc8K[n8K[f]][1]] = 'O';
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc8K[n8K[f]][binAcnt],hc8K[n8K[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc8K[n8K[f]][binAcnt]][Bonds_cnb_j(hc8K[n8K[f]][binAcnt],hc8K[n8K[f]][binBcnt])],BLDistro8K,&BLDistroNoSamples8K);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc8K[n8K[f]][binAcnt]]==1) {
                            nA8K++;
                            number_of_A++;
                        }
                        else nB8K++;
                    }
                    n_distro_8K[number_of_A]++;
                }
                
                ++n8K[f];
            }
        }
        }
    }

    for(i=0; i<N; ++i) s8K[i]=ach[i];
    free(ach);
}

void Clusters_Get9A_D3h(int f) {    // Detect 9A D3h clusters
    int i, j, j2, k, l, m, n;
    int db[2], ob[4];
    int flg;
    char *ach, errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=9;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=3;
    int sub[3];

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
                        if (doClusBLDeviation==1) {
                            bl_mom_9A=resize_1D_double(bl_mom_9A,m9A,m9A+incrStatic);
                        }
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
                
                // hc9A key: (SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_lowest_s, SP4_to_0_in_SP4_lowest_s, SP4_to_1_in_SP4_lowest_s, s_lowest, s_to_0_in_SP4_lowest_s,s_to_2_in_SP4_lowest_s)
                if (doDynamics==1 && dyn_m9A!=-1) {
                    if (doSubClusts==1 && dyn_msp4b!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp4b[i];
                        sub[1]=dyn_up_sp4b[mem_sp4b[sp4b[i][j2]][j]];
                        sub[2]=dyn_up_sp4b[k];
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc9A[n9A[f]], f, clusSize, &dyn_n9A, &dyn_m9A, &dyn_l9A, &dyn_hc9A, do_up, dummy_up, n9A[f], do_sub, n_sub, &dyn_sub_9A, sub);
                }
                if(ach[hc9A[n9A[f]][0]] == 'C') ach[hc9A[n9A[f]][0]] = 'B';
                if(ach[hc9A[n9A[f]][1]] == 'C') ach[hc9A[n9A[f]][1]] = 'B';
                if(ach[hc9A[n9A[f]][3]] == 'C') ach[hc9A[n9A[f]][3]] = 'B';
                if(ach[hc9A[n9A[f]][4]] == 'C') ach[hc9A[n9A[f]][4]] = 'B';
                if(ach[hc9A[n9A[f]][6]] == 'C') ach[hc9A[n9A[f]][6]] = 'B';
                if(ach[hc9A[n9A[f]][7]] == 'C') ach[hc9A[n9A[f]][7]] = 'B';
                ach[hc9A[n9A[f]][2]] = 'O';
                ach[hc9A[n9A[f]][5]] = 'O';
                ach[hc9A[n9A[f]][8]] = 'O';
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc9A[n9A[f]][binAcnt],hc9A[n9A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc9A[n9A[f]][binAcnt]][Bonds_cnb_j(hc9A[n9A[f]][binAcnt],hc9A[n9A[f]][binBcnt])],BLDistro9A,&BLDistroNoSamples9A);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc9A[n9A[f]][binAcnt]]==1) {
                            nA9A++;
                            number_of_A++;
                        }
                        else nB9A++;
                    }
                    n_distro_9A[number_of_A]++;
                }
                
                ++n9A[f];
                
                break;
            }
        }
        }
    }

    for(i=0; i<N; ++i) s9A[i]=ach[i];
    free(ach);
}

void Clusters_Get9B_10B_11B_11E_12D(int f) {    // Detect 9B, 10B, 11A, 11E & 12D
    int sp1, sp2i, sp2j;
    int sp5com[2];
    int i, j, k, l, m;
    int flg, fb1, fb2;
    char *ach1, *ach1_cen, *ach1_shell, *ach2, *ach2_cen, *ach2_shell, *ach3, *ach4, *ach5, *ach5_cen, *ach5_shell, *ach6, *ach6_cen, *ach6_shell, errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=9;
    int do_up=0;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_9B=resize_1D_double(bl_mom_9B,m9B,m9B+incrStatic);
                }
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
            if (doDynamics==1 && dyn_m9B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                    sub[1]=dyn_up_sp5c[j];
                    quickSort(&sub[0],2);
                }
                else do_sub=0;
                if (doSubClusts==1) do_up=1;
                else do_up=0;
                Dyn_add(hc9B[n9B[f]], f, clusSize, &dyn_n9B, &dyn_m9B, &dyn_l9B, &dyn_hc9B, do_up, dyn_up_9B, n9B[f], do_sub, n_sub, &dyn_sub_9B, sub);
            }
            if (ach1[hc9B[n9B[f]][0]] == 'C') ach1[hc9B[n9B[f]][0]] = ach1_shell[hc9B[n9B[f]][0]] = 'B';
            if (ach1[hc9B[n9B[f]][1]] == 'C') ach1[hc9B[n9B[f]][1]] = ach1_shell[hc9B[n9B[f]][1]] = 'B';
            if (ach1[hc9B[n9B[f]][2]] == 'C') ach1[hc9B[n9B[f]][2]] = ach1_shell[hc9B[n9B[f]][2]] = 'B';
            if (ach1[hc9B[n9B[f]][3]] == 'C') ach1[hc9B[n9B[f]][3]] = ach1_shell[hc9B[n9B[f]][3]] = 'B';
            if (ach1[hc9B[n9B[f]][4]] == 'C') ach1[hc9B[n9B[f]][4]] = ach1_shell[hc9B[n9B[f]][4]] = 'B';
            if (ach1[hc9B[n9B[f]][5]] == 'C') ach1[hc9B[n9B[f]][5]] = ach1_shell[hc9B[n9B[f]][5]] = 'B';
            ach1[hc9B[n9B[f]][6]] = ach1_shell[hc9B[n9B[f]][6]] = 'O';
            ach1[hc9B[n9B[f]][7]] = ach1_shell[hc9B[n9B[f]][7]] = 'O';
            ach1[hc9B[n9B[f]][8]] = ach1_cen[hc9B[n9B[f]][8]] = 'O';
            
            if (doBondedCen==1) {
                n_bonded_to_cen_9B+=cnb[hc9B[n9B[f]][8]];
                n_distro_bonded_to_cen_9B[cnb[hc9B[n9B[f]][8]]]++;
            }
            
            if (do10B==1) Clusters_Get10B_C3v(f,i,j,ach2,ach2_cen,ach2_shell,ach6,ach6_cen,ach6_shell); 
            if (do11B==1) {
                if (Clusters_Get11B_C2v(f,ach5,ach5_cen,ach5_shell)) {
                    ach5_cen[hc9B[n9B[f]][8]] = 'O';
                    ++n11B[f];  
                }
            }
            if (do11E==1) Clusters_Get11E_12D(f,i,j,sp1,sp2i,sp2j,ach3,ach4);
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc9B[n9B[f]][binAcnt],hc9B[n9B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc9B[n9B[f]][binAcnt]][Bonds_cnb_j(hc9B[n9B[f]][binAcnt],hc9B[n9B[f]][binBcnt])],BLDistro9B,&BLDistroNoSamples9B);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc9B[n9B[f]][binAcnt]]==1) {
                        nA9B++;
                        number_of_A++;
                    }
                    else nB9B++;
                    
                    if (binAcnt!=8) {
                        if (rtype[hc9B[n9B[f]][binAcnt]]==1) nA_shell_9B++;
                        else nB_shell_9B++;
                    }
                }
                n_distro_9B[number_of_A]++;
                
                if (rtype[hc9B[n9B[f]][8]]==1) {
                    nA_cen_9B++;
                    n_distro_cen_9B[1]++;
                    n_distro_shell_9B[number_of_A-1]++;
                }
                else {
                    nB_cen_9B++;
                    n_distro_cen_9B[0]++;
                    n_distro_shell_9B[number_of_A]++;
                }
            }

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

    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=10;
    int do_up=0;
    int do_sub=0;
    int n_sub=3;
    int sub[3];

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
                if (doClusBLDeviation==1) {
                    bl_mom_10B=resize_1D_double(bl_mom_10B,m10B,m10B+incrStatic);
                }
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
            if (doDynamics==1 && dyn_m10B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                    sub[1]=dyn_up_sp5c[j];
                    sub[2]=dyn_up_sp5c[k];
                    quickSort(&sub[0],3);
                }
                else do_sub=0;
                if (doSubClusts==1) do_up=1;
                else do_up=0;
                Dyn_add(hc10B[n10B[f]], f, clusSize, &dyn_n10B, &dyn_m10B, &dyn_l10B, &dyn_hc10B, do_up, dyn_up_10B, n10B[f], do_sub, n_sub, &dyn_sub_10B, sub);
            }
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
            
            if (doBondedCen==1) {
                n_bonded_to_cen_10B+=cnb[hc10B[n10B[f]][9]];
                n_distro_bonded_to_cen_10B[cnb[hc10B[n10B[f]][9]]]++;
            }
            
            if (do11W==1) Clusters_Get11W_Cs(f,ach1,ach1_cen,ach1_shell);   
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc10B[n10B[f]][binAcnt],hc10B[n10B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc10B[n10B[f]][binAcnt]][Bonds_cnb_j(hc10B[n10B[f]][binAcnt],hc10B[n10B[f]][binBcnt])],BLDistro10B,&BLDistroNoSamples10B);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc10B[n10B[f]][binAcnt]]==1) {
                        nA10B++;
                        number_of_A++;
                    }
                    else nB10B++;
                    
                    if (binAcnt!=9) {
                        if (rtype[hc10B[n10B[f]][binAcnt]]==1) nA_shell_10B++;
                        else nB_shell_10B++;
                    }
                }
                n_distro_10B[number_of_A]++;
                
                if (rtype[hc10B[n10B[f]][9]]==1) {
                    nA_cen_10B++;
                    n_distro_cen_10B[1]++;
                    n_distro_shell_10B[number_of_A-1]++;
                }
                else {
                    nB_cen_10B++;
                    n_distro_cen_10B[0]++;
                    n_distro_shell_10B[number_of_A]++;
                }
            }
            
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
                if (doClusBLDeviation==1) {
                    bl_mom_10B=resize_1D_double(bl_mom_10B,m10B,m10B+incrStatic);
                }
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
            if (doDynamics==1 && dyn_m10B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                    sub[1]=dyn_up_sp5c[j];
                    sub[2]=dyn_up_sp5c[k];
                    quickSort(&sub[0],3);
                }
                else do_sub=0;
                if (doSubClusts==1) do_up=1;
                else do_up=0;
                Dyn_add(hc10B[n10B[f]], f, clusSize, &dyn_n10B, &dyn_m10B, &dyn_l10B, &dyn_hc10B, do_up, dyn_up_10B, n10B[f], do_sub, n_sub, &dyn_sub_10B, sub);
            }
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
            
            if (doBondedCen==1) {
                n_bonded_to_cen_10B+=cnb[hc10B[n10B[f]][9]];
                n_distro_bonded_to_cen_10B[cnb[hc10B[n10B[f]][9]]]++;
            }
            
            if (do11W==1) Clusters_Get11W_Cs(f,ach1,ach1_cen,ach1_shell);   
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc10B[n10B[f]][binAcnt],hc10B[n10B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc10B[n10B[f]][binAcnt]][Bonds_cnb_j(hc10B[n10B[f]][binAcnt],hc10B[n10B[f]][binBcnt])],BLDistro10B,&BLDistroNoSamples10B);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc10B[n10B[f]][binAcnt]]==1) {
                        nA10B++;
                        number_of_A++;
                    }
                    else nB10B++;
                    
                    if (binAcnt!=9) {
                        if (rtype[hc10B[n10B[f]][binAcnt]]==1) nA_shell_10B++;
                        else nB_shell_10B++;
                    }
                }
                n_distro_10B[number_of_A]++;
                
                if (rtype[hc10B[n10B[f]][9]]==1) {
                    nA_cen_10B++;
                    n_distro_cen_10B[1]++;
                    n_distro_shell_10B[number_of_A-1]++;
                }
                else {
                    nB_cen_10B++;
                    n_distro_cen_10B[0]++;
                    n_distro_shell_10B[number_of_A]++;
                }
            }
            
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=11;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=1;
    int sub[1];

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
        if (doClusBLDeviation==1) {
            bl_mom_11B=resize_1D_double(bl_mom_11B,m11B,m11B+incrStatic);
        }
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
    
    if (doDynamics==1 && dyn_m11B!=-1) {
        if (doSubClusts==1 && dyn_m9B!=-1) {
            do_sub=1;
            sub[0]=dyn_up_9B[n9B[f]];
        }
        else do_sub=0;
        do_up=0;
        Dyn_add(hc11B[n11B[f]], f, clusSize, &dyn_n11B, &dyn_m11B, &dyn_l11B, &dyn_hc11B, do_up, dummy_up, n11B[f], do_sub, n_sub, &dyn_sub_11B, sub);
    }
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
    
    if (doBondedCen==1) {
        n_bonded_to_cen_11B+=cnb[hc11B[n11B[f]][8]];
        n_distro_bonded_to_cen_11B[cnb[hc11B[n11B[f]][8]]]++;
    }
    
    if (doClusBLDistros==1) {
        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                if (Bonds_BondCheck(hc11B[n11B[f]][binAcnt],hc11B[n11B[f]][binBcnt])==1) {
                    Bonds_TickBLDistro(bondlengths[hc11B[n11B[f]][binAcnt]][Bonds_cnb_j(hc11B[n11B[f]][binAcnt],hc11B[n11B[f]][binBcnt])],BLDistro11B,&BLDistroNoSamples11B);
                }
            }
        }
    }
    
    if (doClusComp==1) {
        number_of_A=0;
        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
            if (rtype[hc11B[n11B[f]][binAcnt]]==1) {
                nA11B++;
                number_of_A++;
            }
            else nB11B++;
            
            if (binAcnt!=8) {
                if (rtype[hc11B[n11B[f]][binAcnt]]==1) nA_shell_11B++;
                else nB_shell_11B++;
            }
        }
        n_distro_11B[number_of_A]++;
    
        if (rtype[hc11B[n11B[f]][8]]==1) {
            nA_cen_11B++;
            n_distro_cen_11B[1]++;
            n_distro_shell_11B[number_of_A-1]++;
        }
        else {
            nB_cen_11B++;
            n_distro_cen_11B[0]++;
            n_distro_shell_11B[number_of_A]++;
        }
    }
    
    return 1;
}

int Clusters_Get11W_Cs(int f, char *ach, char *ach_cen, char *ach_shell) {  // Detect 11W C2s clusters 
    //  11W is the ground state of the 11 Wahnstrom particles
    // Call from 10B, the central particle must have exactly 10 particles bonded to it. 
    // One extra particle bonded to 10B central particle but not bonded to three shell spindles of 7A in 10B
    
    int k, l, m;
    int ep=-1;
    int break_out;
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=11;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=1;
    int sub[1];

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
        if (doClusBLDeviation==1) {
            bl_mom_11W=resize_1D_double(bl_mom_11W,m11W,m11W+incrStatic);
        }
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
    
    
    if (doDynamics==1 && dyn_m11W!=-1) {
        if (doSubClusts==1 && dyn_m10B!=-1) {
            do_sub=1;
            sub[0]=dyn_up_10B[n10B[f]];
        }
        else do_sub=0;
        do_up=0;
        Dyn_add(hc11W[n11W[f]], f, clusSize, &dyn_n11W, &dyn_m11W, &dyn_l11W, &dyn_hc11W, do_up, dummy_up, n11W[f], do_sub, n_sub, &dyn_sub_11W, sub);
    }
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
    
    if (doBondedCen==1) {
        n_bonded_to_cen_11W+=cnb[hc11W[n11W[f]][9]];
        n_distro_bonded_to_cen_11W[cnb[hc11W[n11W[f]][9]]]++;
    }
    
    if (doClusBLDistros==1) {
        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                if (Bonds_BondCheck(hc11W[n11W[f]][binAcnt],hc11W[n11W[f]][binBcnt])==1) {
                    Bonds_TickBLDistro(bondlengths[hc11W[n11W[f]][binAcnt]][Bonds_cnb_j(hc11W[n11W[f]][binAcnt],hc11W[n11W[f]][binBcnt])],BLDistro11W,&BLDistroNoSamples11W);
                }
            }
        }
    }
    
    if (doClusComp==1) {
        number_of_A=0;
        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
            if (rtype[hc11W[n11W[f]][binAcnt]]==1) {
                nA11W++;
                number_of_A++;
            }
            else nB11W++;
            
            if (binAcnt!=10) {
                if (rtype[hc11W[n11W[f]][binAcnt]]==1) nA_shell_11W++;
                else nB_shell_11W++;
            }
        }
        n_distro_11W[number_of_A]++;
    
        if (rtype[hc11W[n11W[f]][9]]==1) {
            nA_cen_11W++;
            n_distro_cen_11W[1]++;
            n_distro_shell_11W[number_of_A-1]++;
        }
        else {
            nB_cen_11W++;
            n_distro_cen_11W[0]++;
            n_distro_shell_11W[number_of_A]++;
        }
    }
    
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=11;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=3;
    int sub[3];
    
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
                    if (doClusBLDeviation==1) {
                        bl_mom_11E=resize_1D_double(bl_mom_11E,m11E,m11E+incrStatic);
                    }
                    m11E=m11E+incrStatic;
                }
                
                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (l=0;l<11;l++) hc11E[n11E[f]][l]=trial[l];

                if (doDynamics==1 && dyn_m11E!=-1) {
                    if (doSubClusts==1 && dyn_msp5c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp5c[i];
                        sub[1]=dyn_up_sp5c[j];
                        sub[2]=dyn_up_sp5c[k];
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc11E[n11E[f]], f, clusSize, &dyn_n11E, &dyn_m11E, &dyn_l11E, &dyn_hc11E, do_up, dummy_up, n11E[f], do_sub, n_sub, &dyn_sub_11E, sub);
                }
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
                
                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2i, sp5c[k][6], ach2);
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc11E[n11E[f]][binAcnt]][Bonds_cnb_j(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])],BLDistro11E,&BLDistroNoSamples11E);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc11E[n11E[f]][binAcnt]]==1) {
                            nA11E++;
                            number_of_A++;
                        }
                        else nB11E++;
                    }
                    n_distro_11E[number_of_A]++;
                }
                
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
                    if (doClusBLDeviation==1) {
                        bl_mom_11E=resize_1D_double(bl_mom_11E,m11E,m11E+incrStatic);
                    }
                    m11E=m11E+incrStatic;
                }
                
                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (l=0;l<11;l++) hc11E[n11E[f]][l]=trial[l];

                if (doDynamics==1 && dyn_m11E!=-1) {
                    if (doSubClusts==1 && dyn_msp5c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp5c[i];
                        sub[1]=dyn_up_sp5c[j];
                        sub[2]=dyn_up_sp5c[k];
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc11E[n11E[f]], f, clusSize, &dyn_n11E, &dyn_m11E, &dyn_l11E, &dyn_hc11E, do_up, dummy_up, n11E[f], do_sub, n_sub, &dyn_sub_11E, sub);
                }
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
                
                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2i, sp5c[k][5], ach2);
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc11E[n11E[f]][binAcnt]][Bonds_cnb_j(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])],BLDistro11E,&BLDistroNoSamples11E);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc11E[n11E[f]][binAcnt]]==1) {
                            nA11E++;
                            number_of_A++;
                        }
                        else nB11E++;
                    }
                    n_distro_11E[number_of_A]++;
                }
                
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
                    if (doClusBLDeviation==1) {
                        bl_mom_11E=resize_1D_double(bl_mom_11E,m11E,m11E+incrStatic);
                    }
                    m11E=m11E+incrStatic;
                }
                
                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (l=0;l<11;l++) hc11E[n11E[f]][l]=trial[l];

                if (doDynamics==1 && dyn_m11E!=-1) {
                    if (doSubClusts==1 && dyn_msp5c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp5c[i];
                        sub[1]=dyn_up_sp5c[j];
                        sub[2]=dyn_up_sp5c[k];
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc11E[n11E[f]], f, clusSize, &dyn_n11E, &dyn_m11E, &dyn_l11E, &dyn_hc11E, do_up, dummy_up, n11E[f], do_sub, n_sub, &dyn_sub_11E, sub);
                }
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
                
                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2j, sp5c[k][6], ach2);
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc11E[n11E[f]][binAcnt]][Bonds_cnb_j(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])],BLDistro11E,&BLDistroNoSamples11E);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc11E[n11E[f]][binAcnt]]==1) {
                            nA11E++;
                            number_of_A++;
                        }
                        else nB11E++;
                    }
                    n_distro_11E[number_of_A]++;
                }
                
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
                    if (doClusBLDeviation==1) {
                        bl_mom_11E=resize_1D_double(bl_mom_11E,m11E,m11E+incrStatic);
                    }
                    m11E=m11E+incrStatic;
                }
                
                quickSort(&trial[0],2);
                quickSort(&trial[2],2);
                quickSort(&trial[4],7);
                for (l=0;l<11;l++) hc11E[n11E[f]][l]=trial[l];

                if (doDynamics==1 && dyn_m11E!=-1) {
                    if (doSubClusts==1 && dyn_msp5c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp5c[i];
                        sub[1]=dyn_up_sp5c[j];
                        sub[2]=dyn_up_sp5c[k];
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc11E[n11E[f]], f, clusSize, &dyn_n11E, &dyn_m11E, &dyn_l11E, &dyn_hc11E, do_up, dummy_up, n11E[f], do_sub, n_sub, &dyn_sub_11E, sub);
                }
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
                
                // 12D - is there an SP5 spindle, > k & j, with sp5c[k][6] & sp2i
                if (do12D==1) n12D[f] += Clusters_Get12D_D2d(f, i, j, k, sp2j, sp5c[k][5], ach2);
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc11E[n11E[f]][binAcnt]][Bonds_cnb_j(hc11E[n11E[f]][binAcnt],hc11E[n11E[f]][binBcnt])],BLDistro11E,&BLDistroNoSamples11E);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc11E[n11E[f]][binAcnt]]==1) {
                            nA11E++;
                            number_of_A++;
                        }
                        else nB11E++;
                    }
                    n_distro_11E[number_of_A]++;
                }
                
                ++n11E[f];  
            }           
        }
    }
}

int Clusters_Get12D_D2d(int f, int i, int j, int k, int sp1, int sp2, char *ach) {  // Return 1 if 12B is also 11E
    int l, m, n, o, p, q;
    int flg1, flg2;
    int break_out;
    int trial[12];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=12;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=4;
    int sub[4];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_12D=resize_1D_double(bl_mom_12D,m12D,m12D+incrStatic);
                }
                m12D=m12D+incrStatic;
            }
            quickSort(&trial[0],4);
            quickSort(&trial[4],8);
            
            for(m=0; m<12; ++m) hc12D[n12D[f]][m] = trial[m];   
            
            if (doDynamics==1 && dyn_m12D!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                    sub[1]=dyn_up_sp5c[j];
                    sub[2]=dyn_up_sp5c[k];
                    sub[3]=dyn_up_sp5c[l];
                    quickSort(&sub[0],4);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc12D[n12D[f]], f, clusSize, &dyn_n12D, &dyn_m12D, &dyn_l12D, &dyn_hc12D, do_up, dummy_up, n12D[f], do_sub, n_sub, &dyn_sub_12D, sub);
            }
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
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc12D[n12D[f]][binAcnt],hc12D[n12D[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc12D[n12D[f]][binAcnt]][Bonds_cnb_j(hc12D[n12D[f]][binAcnt],hc12D[n12D[f]][binBcnt])],BLDistro12D,&BLDistroNoSamples12D);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc12D[n12D[f]][binAcnt]]==1) {
                        nA12D++;
                        number_of_A++;
                    }
                    else nB12D++;
                }
                n_distro_12D[number_of_A]++;
            }
            
            return 1;
        }
    }
    return 0;
}

void Clusters_Get9K_10K(int f)  { // Detect 9K & 10K clusters
    // Made from 2 sp4c clusters with a common sp4c spindle 
    // particl and two common SP4 ring particles
    char *ach_1, *ach_cen_1, *ach_shell_1, *ach_2, *ach_cen_2, *ach_shell_2;
    int i, j2, j, k, l, m;
    int cp[2], scom, sother[2];
    int trial[9];
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=9;
    int do_up=0;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
    cp[0]=cp[1]=scom=sother[0]=sother[1]=-1;
    
    ach_1=malloc(N*sizeof(char));   if (ach_1==NULL) { sprintf(errMsg,"Clusters_Get9K_10K(): ach_1[] malloc out of memory\n");  Error(errMsg); }
    ach_cen_1=malloc(N*sizeof(char));   if (ach_cen_1==NULL) { sprintf(errMsg,"Clusters_Get9K_10K): ach_cen_1[] malloc out of memory\n");   Error(errMsg); }
    ach_shell_1=malloc(N*sizeof(char)); if (ach_shell_1==NULL) { sprintf(errMsg,"Clusters_Get9K_10K(): ach_shell_1[] malloc out of memory\n");  Error(errMsg); }
    ach_2=malloc(N*sizeof(char));   if (ach_2==NULL) { sprintf(errMsg,"Clusters_Get9K_10K(): ach_2[] malloc out of memory\n");  Error(errMsg); }
    ach_cen_2=malloc(N*sizeof(char));   if (ach_cen_2==NULL) { sprintf(errMsg,"Clusters_Get9K_10K): ach_cen_2[] malloc out of memory\n");   Error(errMsg); }
    ach_shell_2=malloc(N*sizeof(char)); if (ach_shell_2==NULL) { sprintf(errMsg,"Clusters_Get9K_10K(): ach_shell_2[] malloc out of memory\n");  Error(errMsg); }
    for (i=0; i<N; ++i) {
        ach_1[i] = 'C';
        ach_cen_1[i] = 'C';
        ach_shell_1[i] = 'C';
        ach_2[i] = 'C';
        ach_cen_2[i] = 'C';
        ach_shell_2[i] = 'C';
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
                if (doClusBLDeviation==1) {
                    bl_mom_9K=resize_1D_double(bl_mom_9K,m9K,m9K+incrStatic);
                }
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

            if (doDynamics==1 && dyn_m9K!=-1) {
                if (doSubClusts==1 && dyn_m6A!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp4c[i];
                    sub[1]=dyn_up_sp4c[mem_sp4c[sp4c[i][j2]][j]];
                    quickSort(&sub[0],2);
                }
                else do_sub=0;
                if (doSubClusts==1) do_up=1;
                else do_up=0;
                Dyn_add(hc9K[n9K[f]], f, clusSize, &dyn_n9K, &dyn_m9K, &dyn_l9K, &dyn_hc9K, do_up, dyn_up_9K, n9K[f], do_sub, n_sub, &dyn_sub_9K, sub);
            }
            if(ach_1[hc9K[n9K[f]][0]]  == 'C') ach_1[hc9K[n9K[f]][0]] = ach_shell_1[hc9K[n9K[f]][0]] = 'B';
            if(ach_1[hc9K[n9K[f]][1]]  == 'C') ach_1[hc9K[n9K[f]][1]] = ach_shell_1[hc9K[n9K[f]][1]] = 'B';
            if(ach_1[hc9K[n9K[f]][2]]  == 'C') ach_1[hc9K[n9K[f]][2]] = ach_shell_1[hc9K[n9K[f]][2]] = 'B';
            if(ach_1[hc9K[n9K[f]][3]]  == 'C') ach_1[hc9K[n9K[f]][3]] = ach_shell_1[hc9K[n9K[f]][3]] = 'B';
            if(ach_1[hc9K[n9K[f]][4]]  == 'C') ach_1[hc9K[n9K[f]][4]] = ach_shell_1[hc9K[n9K[f]][4]] = 'B';
            if(ach_1[hc9K[n9K[f]][5]]  == 'C') ach_1[hc9K[n9K[f]][5]] = ach_shell_1[hc9K[n9K[f]][5]] = 'B';
            ach_1[hc9K[n9K[f]][6]] = ach_shell_1[hc9K[n9K[f]][6]] = 'O';
            ach_1[hc9K[n9K[f]][7]] = ach_shell_1[hc9K[n9K[f]][7]] = 'O';
            ach_1[hc9K[n9K[f]][8]] = ach_cen_1[hc9K[n9K[f]][8]] = 'O';

            if (do10K==1) n10K[f] += Clusters_Get10K(f, ach_2, ach_cen_2, ach_shell_2);

            if (doBondedCen==1) {
                n_bonded_to_cen_9K+=cnb[hc9K[n9K[f]][8]];
                n_distro_bonded_to_cen_9K[cnb[hc9K[n9K[f]][8]]]++;
            }
                
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc9K[n9K[f]][binAcnt],hc9K[n9K[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc9K[n9K[f]][binAcnt]][Bonds_cnb_j(hc9K[n9K[f]][binAcnt],hc9K[n9K[f]][binBcnt])],BLDistro9K,&BLDistroNoSamples9K);
                        }
                    }
                }
            }
                
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc9K[n9K[f]][binAcnt]]==1) {
                        nA9K++;
                        number_of_A++;
                    }
                    else nB9K++;
                    
                    if (binAcnt!=8) {
                        if (rtype[hc9K[n9K[f]][binAcnt]]==1) nA_shell_9K++;
                        else nB_shell_9K++;
                    }
                }
                n_distro_9K[number_of_A]++;
                
                if (rtype[hc9K[n9K[f]][8]]==1) {
                    nA_cen_9K++;
                    n_distro_cen_9K[1]++;
                    n_distro_shell_9K[number_of_A-1]++;
                }
                else {
                    nB_cen_9K++;
                    n_distro_cen_9K[0]++;
                    n_distro_shell_9K[number_of_A]++;
                }
            }

            ++n9K[f];
        }
        }
    }

    for (i=0; i<N; ++i) {
        s9K[i]=ach_1[i];
        s9K_cen[i]=ach_cen_1[i];
        s9K_shell[i]=ach_shell_1[i];
        s10K[i]=ach_2[i];
        s10K_cen[i]=ach_cen_2[i];
        s10K_shell[i]=ach_shell_2[i];
    }
    free(ach_1);
    free(ach_cen_1);
    free(ach_shell_1);
    free(ach_2);
    free(ach_cen_2);
    free(ach_shell_2);
}

int Clusters_Get10K(int f, char *ach, char *ach_cen, char *ach_shell) { // Detect 10K clusters
    int i, j, k, ep;
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=10;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=1;
    int sub[1];
    
    ep=-1;

    k=0;    // check exactly one extra particle bonded to 9K central particle, extra particle is separate from any of the 9K shell particles
    for (i=0; i<cnb[hc9K[n9K[f]][8]]; i++) {
        for (j=0; j<8; j++) {
            if (bNums[hc9K[n9K[f]][8]][i]==hc9K[n9K[f]][j]) break;
        }
        if (j==8) {
            if (k>=1) {
                k++;
                break;
            }
            ep=bNums[hc9K[n9K[f]][8]][i];
            k++;
        }
        if (k>=2) break;
    }
    if (k!=1) return 0;

    if(n10K[f] == m10K) { 
        hc10K=resize_2D_int(hc10K,m10K,m10K+incrStatic,clusSize,-1);
        if (doClusBLDeviation==1) {
                bl_mom_10K=resize_1D_double(bl_mom_10K,m10K,m10K+incrStatic);
            }
        m10K=m10K+incrStatic;
    }
    // hc10K key: (common_SP4_1, common_SP4_2, other_SP4*4, other_spindle_1, other_spindle_2, scom, ep)
    
    for (i=0; i<9; i++) {
        hc10K[n10K[f]][i]=hc9K[n9K[f]][i];
    }
    hc10K[n10K[f]][9]=ep;

    if (doDynamics==1 && dyn_m10K!=-1) {
        if (doSubClusts==1 && dyn_msp5b!=-1) {
            do_sub=1;
            sub[0]=dyn_up_9K[n9K[f]];
        }
        else do_sub=0;
        do_up=0;
        Dyn_add(hc10K[n10K[f]], f, clusSize, &dyn_n10K, &dyn_m10K, &dyn_l10K, &dyn_hc10K, do_up, dummy_up, n10K[f], do_sub, n_sub, &dyn_sub_10K, sub);
    }
    if(ach[hc10K[n10K[f]][0]]  == 'C') ach[hc10K[n10K[f]][0]] = ach_shell[hc10K[n10K[f]][0]] = 'B';
    if(ach[hc10K[n10K[f]][1]]  == 'C') ach[hc10K[n10K[f]][1]] = ach_shell[hc10K[n10K[f]][1]] = 'B';
    if(ach[hc10K[n10K[f]][2]]  == 'C') ach[hc10K[n10K[f]][2]] = ach_shell[hc10K[n10K[f]][2]] = 'B';
    if(ach[hc10K[n10K[f]][3]]  == 'C') ach[hc10K[n10K[f]][3]] = ach_shell[hc10K[n10K[f]][3]] = 'B';
    if(ach[hc10K[n10K[f]][4]]  == 'C') ach[hc10K[n10K[f]][4]] = ach_shell[hc10K[n10K[f]][4]] = 'B';
    if(ach[hc10K[n10K[f]][5]]  == 'C') ach[hc10K[n10K[f]][5]] = ach_shell[hc10K[n10K[f]][5]] = 'B';
    ach[hc10K[n10K[f]][6]] = ach_shell[hc10K[n10K[f]][6]] = 'O';
    ach[hc10K[n10K[f]][7]] = ach_shell[hc10K[n10K[f]][7]] = 'O';
    ach[hc10K[n10K[f]][8]] = ach_cen[hc10K[n10K[f]][8]] = 'O';
    ach[hc10K[n10K[f]][9]] = ach_shell[hc10K[n10K[f]][9]] = 'O';
    
    if (doBondedCen==1) {
        n_bonded_to_cen_10K+=cnb[hc10K[n10K[f]][8]];
        n_distro_bonded_to_cen_10K[cnb[hc10K[n10K[f]][8]]]++;
    }
                
    if (doClusBLDistros==1) {
        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                if (Bonds_BondCheck(hc10K[n10K[f]][binAcnt],hc10K[n10K[f]][binBcnt])==1) {
                    Bonds_TickBLDistro(bondlengths[hc10K[n10K[f]][binAcnt]][Bonds_cnb_j(hc10K[n10K[f]][binAcnt],hc10K[n10K[f]][binBcnt])],BLDistro10K,&BLDistroNoSamples10K);
                }
            }
        }
    }
                
    if (doClusComp==1) {
        number_of_A=0;
        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
            if (rtype[hc10K[n10K[f]][binAcnt]]==1) {
                nA10K++;
                number_of_A++;
            }
            else nB10K++;
            
            if (binAcnt!=8) {
                if (rtype[hc10K[n10K[f]][binAcnt]]==1) nA_shell_10K++;
                else nB_shell_10K++;
            }
        }
        n_distro_10K[number_of_A]++;
        
        if (rtype[hc10K[n10K[f]][8]]==1) {
            nA_cen_10K++;
            n_distro_cen_10K[1]++;
            n_distro_shell_10K[number_of_A-1]++;
        }
        else {
            nB_cen_10K++;
            n_distro_cen_10K[0]++;
            n_distro_shell_10K[number_of_A]++;
        }
    }

    return 1;
}

void Clusters_Get10A_C3v(int f) { // Detect 10A D4d clusters
    char *ach;
    int i, j, j2, k, l, m;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=10;
    int *used_sp4b;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                    if (doClusBLDeviation==1) {
                            bl_mom_10A=resize_1D_double(bl_mom_10A,m10A,m10A+incrStatic);
                        }
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
                
                
                if (doDynamics==1 && dyn_m10A!=-1) {
                    if (doSubClusts==1 && dyn_msp4b!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp4b[i];
                        sub[1]=dyn_up_sp4b[mem_sp4b[bNums[sp4b[i][0]][j2]][j]];
                        quickSort(&sub[0],2);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc10A[n10A[f]], f, clusSize, &dyn_n10A, &dyn_m10A, &dyn_l10A, &dyn_hc10A, do_up, dummy_up, n10A[f], do_sub, n_sub, &dyn_sub_10A, sub);
                }
                if(ach[hc10A[n10A[f]][0]] == 'C') ach[hc10A[n10A[f]][0]] = 'B';
                if(ach[hc10A[n10A[f]][1]] == 'C') ach[hc10A[n10A[f]][1]] = 'B';
                if(ach[hc10A[n10A[f]][2]] == 'C') ach[hc10A[n10A[f]][2]] = 'B';
                if(ach[hc10A[n10A[f]][3]] == 'C') ach[hc10A[n10A[f]][3]] = 'B';
                if(ach[hc10A[n10A[f]][4]] == 'C') ach[hc10A[n10A[f]][4]] = 'B';
                if(ach[hc10A[n10A[f]][5]] == 'C') ach[hc10A[n10A[f]][5]] = 'B';
                if(ach[hc10A[n10A[f]][6]] == 'C') ach[hc10A[n10A[f]][6]] = 'B';
                if(ach[hc10A[n10A[f]][7]] == 'C') ach[hc10A[n10A[f]][7]] = 'B';
                ach[hc10A[n10A[f]][8]] = 'O';
                ach[hc10A[n10A[f]][9]] = 'O';
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc10A[n10A[f]][binAcnt],hc10A[n10A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc10A[n10A[f]][binAcnt]][Bonds_cnb_j(hc10A[n10A[f]][binAcnt],hc10A[n10A[f]][binBcnt])],BLDistro10A,&BLDistroNoSamples10A);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc10A[n10A[f]][binAcnt]]==1) {
                            nA10A++;
                            number_of_A++;
                        }
                        else nB10A++;
                    }
                    n_distro_10A[number_of_A]++;
                }
                
                ++n10A[f];
            }
        }
        }
    }

    for (i=0; i<N; ++i) s10A[i]=ach[i];
    free(ach);
    free(used_sp4b);
}

void Clusters_Get10W(int f) { // Detect 10W clusters
    int i, j, k, l, m;
    int sp5b_clusts[5], shell_parts[9];
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=10;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=5;
    int sub[5];
    
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
            if (doClusBLDeviation==1) {
                    bl_mom_10W=resize_1D_double(bl_mom_10W,m10W,m10W+incrStatic);
                }
            m10W=m10W+incrStatic;
        }
        // hc10W key: (sp5bs_common_central_spindle_particle, sp5bs_SP5_ring_shell_particles)
        hc10W[n10W[f]][0] = sp5b[i][5]; 
        for (j=0; j<9; j++) hc10W[n10W[f]][j+1]=shell_parts[j];
        quickSort(&hc10W[n10W[f]][1],9);
        
        if (doDynamics==1 && dyn_m10W!=-1) {
            if (doSubClusts==1 && dyn_msp5b!=-1) {
                do_sub=1;
                sub[0]=dyn_up_sp5b[sp5b_clusts[0]];
                sub[1]=dyn_up_sp5b[sp5b_clusts[1]];
                sub[2]=dyn_up_sp5b[sp5b_clusts[2]];
                sub[3]=dyn_up_sp5b[sp5b_clusts[3]];
                sub[4]=dyn_up_sp5b[sp5b_clusts[4]];
                quickSort(&sub[0],5);
            }
            else do_sub=0;
            do_up=0;
            Dyn_add(hc10W[n10W[f]], f, clusSize, &dyn_n10W, &dyn_m10W, &dyn_l10W, &dyn_hc10W, do_up, dummy_up, n10W[f], do_sub, n_sub, &dyn_sub_10W, sub);
        }
        if(ach[hc10W[n10W[f]][1]] == 'C') ach[hc10W[n10W[f]][1]] = ach_shell[hc10W[n10W[f]][1]] = 'B';
        if(ach[hc10W[n10W[f]][2]] == 'C') ach[hc10W[n10W[f]][2]] = ach_shell[hc10W[n10W[f]][2]] = 'B';
        if(ach[hc10W[n10W[f]][3]] == 'C') ach[hc10W[n10W[f]][3]] = ach_shell[hc10W[n10W[f]][3]] = 'B';
        if(ach[hc10W[n10W[f]][4]] == 'C') ach[hc10W[n10W[f]][4]] = ach_shell[hc10W[n10W[f]][4]] = 'B';
        if(ach[hc10W[n10W[f]][5]] == 'C') ach[hc10W[n10W[f]][5]] = ach_shell[hc10W[n10W[f]][5]] = 'B';
        if(ach[hc10W[n10W[f]][6]] == 'C') ach[hc10W[n10W[f]][6]] = ach_shell[hc10W[n10W[f]][6]] = 'B';
        if(ach[hc10W[n10W[f]][7]] == 'C') ach[hc10W[n10W[f]][7]] = ach_shell[hc10W[n10W[f]][7]] = 'B';
        if(ach[hc10W[n10W[f]][8]] == 'C') ach[hc10W[n10W[f]][8]] = ach_shell[hc10W[n10W[f]][8]] = 'B';
        if(ach[hc10W[n10W[f]][9]] == 'C') ach[hc10W[n10W[f]][9]] = ach_shell[hc10W[n10W[f]][9]] = 'B';
        ach[hc10W[n10W[f]][0]] = ach_cen[hc10W[n10W[f]][0]] = 'O';
                
        if (doBondedCen==1) {
            n_bonded_to_cen_10W+=cnb[hc10W[n10W[f]][0]];
            n_distro_bonded_to_cen_10W[cnb[hc10W[n10W[f]][0]]]++;
        }
                
        if (doClusBLDistros==1) {
            for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                    if (Bonds_BondCheck(hc10W[n10W[f]][binAcnt],hc10W[n10W[f]][binBcnt])==1) {
                        Bonds_TickBLDistro(bondlengths[hc10W[n10W[f]][binAcnt]][Bonds_cnb_j(hc10W[n10W[f]][binAcnt],hc10W[n10W[f]][binBcnt])],BLDistro10W,&BLDistroNoSamples10W);
                    }
                }
            }
        }
                
        if (doClusComp==1) {
            number_of_A=0;
            for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                if (rtype[hc10W[n10W[f]][binAcnt]]==1) {
                    nA10W++;
                    number_of_A++;
                }
                else nB10W++;
            
                if (binAcnt!=0) {
                    if (rtype[hc10W[n10W[f]][binAcnt]]==1) nA_shell_10W++;
                    else nB_shell_10W++;
                }
            }
            n_distro_10W[number_of_A]++;
            
            if (rtype[hc10W[n10W[f]][0]]==1) {
                nA_cen_10W++;
                n_distro_cen_10W[1]++;
                n_distro_shell_10W[number_of_A-1]++;
            }
            else {
                nB_cen_10W++;
                n_distro_cen_10W[0]++;
                n_distro_shell_10W[number_of_A]++;
            }
        }
            
        ++n10W[f];
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

void Clusters_Get11A_12K(int f) { // Detect 11A D4d & 12K clusters
    //  Difficult to be desisive about. Made from 2 sp4c clusters with a common sp4 spindle 
    // particle. Big gaps in the two 4 membered rings. Does work if the bond length is large
    // enough.
    char *ach_1, *ach_cen_1, *ach_shell_1, *ach_2, *ach_cen_2, *ach_shell_2;
    int i, j, k, l, m, n;
    int scom, sother[2], shell_SP3[8][3];
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=11;
    int do_up=0;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
    ach_1=malloc(N*sizeof(char));   if (ach_1==NULL) { sprintf(errMsg,"Clusters_Get11A_12K(): ach_1[] malloc out of memory\n"); Error(errMsg); }
    ach_cen_1=malloc(N*sizeof(char));   if (ach_cen_1==NULL) { sprintf(errMsg,"Clusters_Get11A_12K(): ach_cen_1[] malloc out of memory\n"); Error(errMsg); }
    ach_shell_1=malloc(N*sizeof(char)); if (ach_shell_1==NULL) { sprintf(errMsg,"Clusters_Get11A_12K(): ach_shell_1[] malloc out of memory\n"); Error(errMsg); }
    ach_2=malloc(N*sizeof(char));   if (ach_2==NULL) { sprintf(errMsg,"Clusters_Get11A_12K(): ach_2[] malloc out of memory\n"); Error(errMsg); }
    ach_cen_2=malloc(N*sizeof(char));   if (ach_cen_2==NULL) { sprintf(errMsg,"Clusters_Get11A_12K(): ach_cen_2[] malloc out of memory\n"); Error(errMsg); }
    ach_shell_2=malloc(N*sizeof(char)); if (ach_shell_2==NULL) { sprintf(errMsg,"Clusters_Get11A_12K(): ach_shell_2[] malloc out of memory\n"); Error(errMsg); }
    for (i=0; i<N; ++i) {
        ach_1[i] = 'C';
        ach_cen_1[i] = 'C';
        ach_shell_1[i] = 'C';
        ach_2[i] = 'C';
        ach_cen_2[i] = 'C';
        ach_shell_2[i] = 'C';
    }
    
    scom=sother[0]=sother[1]=-1;
    
    for(i=0; i<nsp4c[f]-1; ++i){
        // POSSIBLE IMPROVEMENT: loop over all sp4c clusters for spindles of sp4c_i
        for(j=i+1; j<nsp4c[f]; ++j) {
            m=0;
            for (k=4; k<6; k++) {
                for (l=4; l<6; l++) {
                    if (sp4c[i][k] == sp4c[j][l]) {
                        if (m>=1) {
                            m++;
                            break;
                        }
                        scom=sp4c[i][k];
                        m++;
                    }
                }
                if (m>=2) break;
            }
            if(m!=1) continue;  // one common spindle
            // ERROR !! need to check uncommon spindles are not in other cluster at all
            if (scom==sp4c[i][4]) sother[0]=sp4c[i][5];
            else sother[0]=sp4c[i][4];
            if (scom==sp4c[j][4]) sother[1]=sp4c[j][5];
            else sother[1]=sp4c[j][4];
            
            for(k=0; k<4; ++k) {
                for(l=0; l<4; ++l) {
                    if(sp4c[i][k] == sp4c[j][l]) break;     
                }
                if(l<4) break;
            }
            if(k<4) continue; // no common ring particles
            
            for (k=0; k<8; ++k) {
                for (l=0; l<3; ++l) shell_SP3[k][l]=-1;
            }
            
            n=0;
            for(k=0; k<4; ++k) {
                m = 0;
                shell_SP3[n][m]=sp4c[i][k];
                for(l=0; l<4; ++l) {
                    if(Bonds_BondCheck(sp4c[i][k], sp4c[j][l])) {
                        if (m==2) {
                            m++;
                            break;
                        }
                        m++;
                        shell_SP3[n][m]=sp4c[j][l];
                    }
                }
                if(m!=2) break;
                n++;
            }
            if (k!=4) continue;
            if (n!=4) continue;
            
            for(k=0; k<4; ++k) {
                m = 0;
                shell_SP3[n][m]=sp4c[j][k];
                for(l=0; l<4; ++l) {
                    if(Bonds_BondCheck(sp4c[j][k], sp4c[i][l])) {
                        if (m==2) {
                            m++;
                            break;
                        }
                        m++;
                        shell_SP3[n][m]=sp4c[i][l];
                    }
                }
                if(m!=2) break;
                n++;
            }
            if (k!=4) continue;
            if (n!=8) continue;
            
            if(n11A[f] == m11A) { 
                hc11A=resize_2D_int(hc11A,m11A,m11A+incrStatic,clusSize,-1);
                if (doClusBLDeviation==1) {
                        bl_mom_11A=resize_1D_double(bl_mom_11A,m11A,m11A+incrStatic);
                    }
                m11A=m11A+incrStatic;
            }
        
            // hc11A key: (SP4 going up, sd going up, scom)
            
            hc11A[n11A[f]][0] = sp4c[i][0];
            hc11A[n11A[f]][1] = sp4c[i][1];
            hc11A[n11A[f]][2] = sp4c[i][2];
            hc11A[n11A[f]][3] = sp4c[i][3];
            hc11A[n11A[f]][4] = sp4c[j][0];
            hc11A[n11A[f]][5] = sp4c[j][1];
            hc11A[n11A[f]][6] = sp4c[j][2];
            hc11A[n11A[f]][7] = sp4c[j][3];
            hc11A[n11A[f]][8] = sother[0];
            hc11A[n11A[f]][9] = sother[1];
            hc11A[n11A[f]][10] = scom;                 
        
            quickSort(&hc11A[n11A[f]][0],8);
            quickSort(&hc11A[n11A[f]][8],2);
                
            if (doDynamics==1 && dyn_m11A!=-1) {
                if (doSubClusts==1 && dyn_m6A!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp4c[i];
                    sub[1]=dyn_up_sp4c[j];
                    quickSort(&sub[0],2);
                }
                else do_sub=0;
                if (doSubClusts==1) do_up=1;
                else do_up=0;
                Dyn_add(hc11A[n11A[f]], f, clusSize, &dyn_n11A, &dyn_m11A, &dyn_l11A, &dyn_hc11A, do_up, dyn_up_11A, n11A[f], do_sub, n_sub, &dyn_sub_11A, sub);
            }
            if(ach_1[hc11A[n11A[f]][0]]  == 'C') ach_1[hc11A[n11A[f]][0]] = ach_shell_1[hc11A[n11A[f]][0]] = 'B';
            if(ach_1[hc11A[n11A[f]][1]]  == 'C') ach_1[hc11A[n11A[f]][1]] = ach_shell_1[hc11A[n11A[f]][1]] = 'B';
            if(ach_1[hc11A[n11A[f]][2]]  == 'C') ach_1[hc11A[n11A[f]][2]] = ach_shell_1[hc11A[n11A[f]][2]] = 'B';
            if(ach_1[hc11A[n11A[f]][3]]  == 'C') ach_1[hc11A[n11A[f]][3]] = ach_shell_1[hc11A[n11A[f]][3]] = 'B';
            if(ach_1[hc11A[n11A[f]][4]]  == 'C') ach_1[hc11A[n11A[f]][4]] = ach_shell_1[hc11A[n11A[f]][4]] = 'B';
            if(ach_1[hc11A[n11A[f]][5]]  == 'C') ach_1[hc11A[n11A[f]][5]] = ach_shell_1[hc11A[n11A[f]][5]] = 'B';
            if(ach_1[hc11A[n11A[f]][6]]  == 'C') ach_1[hc11A[n11A[f]][6]] = ach_shell_1[hc11A[n11A[f]][6]] = 'B';
            if(ach_1[hc11A[n11A[f]][7]]  == 'C') ach_1[hc11A[n11A[f]][7]] = ach_shell_1[hc11A[n11A[f]][7]] = 'B';
            ach_1[hc11A[n11A[f]][8]] = ach_shell_1[hc11A[n11A[f]][8]] = 'O';
            ach_1[hc11A[n11A[f]][9]] = ach_shell_1[hc11A[n11A[f]][9]] = 'O';
            ach_1[hc11A[n11A[f]][10]] = ach_cen_1[hc11A[n11A[f]][10]] = 'O';
            
            if (do12K==1) { 
                for (k=0; k<8; k++) {
                    n12K[f]+=Clusters_Get12K(f,shell_SP3[k][0],shell_SP3[k][1],shell_SP3[k][2],ach_2, ach_cen_2, ach_shell_2);
                }
            }
                
            if (doBondedCen==1) {
                n_bonded_to_cen_11A+=cnb[hc11A[n11A[f]][10]];
                n_distro_bonded_to_cen_11A[cnb[hc11A[n11A[f]][10]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc11A[n11A[f]][binAcnt],hc11A[n11A[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc11A[n11A[f]][binAcnt]][Bonds_cnb_j(hc11A[n11A[f]][binAcnt],hc11A[n11A[f]][binBcnt])],BLDistro11A,&BLDistroNoSamples11A);
                        }
                    }
                }
            }
                
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc11A[n11A[f]][binAcnt]]==1) {
                        nA11A++;
                        number_of_A++;
                    }
                    else nB11A++;
                    
                    if (binAcnt!=10) {
                        if (rtype[hc11A[n11A[f]][binAcnt]]==1) nA_shell_11A++;
                        else nB_shell_11A++;
                    }
                }
                n_distro_11A[number_of_A]++;
                
                if (rtype[hc11A[n11A[f]][10]]==1) {
                    nA_cen_11A++;
                    n_distro_cen_11A[1]++;
                    n_distro_shell_11A[number_of_A-1]++;
                }
                else {
                    nB_cen_11A++;
                    n_distro_cen_11A[0]++;
                    n_distro_shell_11A[number_of_A]++;
                }
            }
                
            ++n11A[f];
        }
    }

    for (i=0; i<N; ++i) {
        s11A[i]=ach_1[i];
        s11A_cen[i]=ach_cen_1[i];
        s11A_shell[i]=ach_shell_1[i];
        s12K[i]=ach_2[i];
        s12K_cen[i]=ach_cen_2[i];
        s12K_shell[i]=ach_shell_2[i];
    }
    free(ach_1);
    free(ach_cen_1);
    free(ach_shell_1);
    free(ach_2);
    free(ach_cen_2);
    free(ach_shell_2);
}

int Clusters_Get12K(int f, int SP3_1, int SP3_2, int SP3_3, char *ach, char *ach_cen, char *ach_shell) {    // Detect 12K clusters
    int i, j;
    int ep, nep;
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=12;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=1;
    int sub[1];
    
    ep=-1;
    
    nep=0;
    for (i=0; i<cnb[SP3_1]; i++) {
        if (Bonds_BondCheck(SP3_2,bNums[SP3_1][i])!=1 || Bonds_BondCheck(SP3_3,bNums[SP3_1][i])!=1) continue;
        for (j=0; j<11; j++) {
            if (bNums[SP3_1][i]==hc11A[n11A[f]][j]) break;
        }
        if (j!=11) continue;
        nep++;
        if (nep>=2) break;
        ep=bNums[SP3_1][i];
    }
    if (nep!=1) return 0;
    
    if(n12K[f] == m12K) { 
        hc12K=resize_2D_int(hc12K,m12K,m12K+incrStatic,clusSize,-1);
        if (doClusBLDeviation==1) {
                bl_mom_12K=resize_1D_double(bl_mom_12K,m12K,m12K+incrStatic);
            }
        m12K=m12K+incrStatic;
    }
    // hc12K key: (SP4 going up, sd going up, scom, ep)

    for (i=0; i<11; i++) hc12K[n12K[f]][i] = hc11A[n11A[f]][i];
    hc12K[n12K[f]][11]=ep;

    if (doDynamics==1 && dyn_m11A!=-1) {
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            do_sub=1;
            sub[0]=dyn_up_11A[n11A[f]];
        }
        else do_sub=0;
        do_up=0;
        Dyn_add(hc12K[n12K[f]], f, clusSize, &dyn_n12K, &dyn_m12K, &dyn_l12K, &dyn_hc12K, do_up, dummy_up, n12K[f], do_sub, n_sub, &dyn_sub_12K, sub);
    }
    if(ach[hc12K[n12K[f]][0]]  == 'C') ach[hc12K[n12K[f]][0]] = ach_shell[hc12K[n12K[f]][0]] = 'B';
    if(ach[hc12K[n12K[f]][1]]  == 'C') ach[hc12K[n12K[f]][1]] = ach_shell[hc12K[n12K[f]][1]] = 'B';
    if(ach[hc12K[n12K[f]][2]]  == 'C') ach[hc12K[n12K[f]][2]] = ach_shell[hc12K[n12K[f]][2]] = 'B';
    if(ach[hc12K[n12K[f]][3]]  == 'C') ach[hc12K[n12K[f]][3]] = ach_shell[hc12K[n12K[f]][3]] = 'B';
    if(ach[hc12K[n12K[f]][4]]  == 'C') ach[hc12K[n12K[f]][4]] = ach_shell[hc12K[n12K[f]][4]] = 'B';
    if(ach[hc12K[n12K[f]][5]]  == 'C') ach[hc12K[n12K[f]][5]] = ach_shell[hc12K[n12K[f]][5]] = 'B';
    if(ach[hc12K[n12K[f]][6]]  == 'C') ach[hc12K[n12K[f]][6]] = ach_shell[hc12K[n12K[f]][6]] = 'B';
    if(ach[hc12K[n12K[f]][7]]  == 'C') ach[hc12K[n12K[f]][7]] = ach_shell[hc12K[n12K[f]][7]] = 'B';
    ach[hc12K[n12K[f]][8]] = ach_shell[hc12K[n12K[f]][8]] = 'O';
    ach[hc12K[n12K[f]][9]] = ach_shell[hc12K[n12K[f]][9]] = 'O';
    ach[hc12K[n12K[f]][10]] = ach_cen[hc12K[n12K[f]][10]] = 'O';
    ach[hc12K[n12K[f]][11]] = ach_shell[hc12K[n12K[f]][11]] = 'O';
    
    if (doBondedCen==1) {
        n_bonded_to_cen_12K+=cnb[hc12K[n12K[f]][10]];
        n_distro_bonded_to_cen_12K[cnb[hc12K[n12K[f]][10]]]++;
    }
                
    if (doClusBLDistros==1) {
        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                if (Bonds_BondCheck(hc12K[n12K[f]][binAcnt],hc12K[n12K[f]][binBcnt])==1) {
                    Bonds_TickBLDistro(bondlengths[hc12K[n12K[f]][binAcnt]][Bonds_cnb_j(hc12K[n12K[f]][binAcnt],hc12K[n12K[f]][binBcnt])],BLDistro12K,&BLDistroNoSamples12K);
                }
            }
        }
    }
    
    if (doClusComp==1) {
        number_of_A=0;
        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
            if (rtype[hc12K[n12K[f]][binAcnt]]==1) {
                nA12K++;
                number_of_A++;
            }
            else nB12K++;
            
            if (binAcnt!=10) {
                if (rtype[hc12K[n12K[f]][binAcnt]]==1) nA_shell_12K++;
                else nB_shell_12K++;
            }
        }
        n_distro_12K[number_of_A]++;
        
        if (rtype[hc12K[n12K[f]][10]]==1) {
            nA_cen_12K++;
            n_distro_cen_12K[1]++;
            n_distro_shell_12K[number_of_A-1]++;
        }
        else {
            nB_cen_12K++;
            n_distro_cen_12K[0]++;
            n_distro_shell_12K[number_of_A]++;
        }
    }
    
    return 1;
}

void Clusters_Get11C_12A(int f) { // Detect 11C Cs & 12A C2v clusters
    char *ach1, *ach1_cen, *ach1_shell, *ach2, *ach2_cen, *ach2_shell;
    int ar[2],sd[2];
    int i, j, k, l, m, ncom, spc;
    int flg;
    int break_out;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=11;
    int do_up=0;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                        if (doClusBLDeviation==1) {
                                bl_mom_11C=resize_1D_double(bl_mom_11C,m11C,m11C+incrStatic);
                            }
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
                    
                    if (doDynamics==1 && dyn_m11C!=-1) {
                        if (doSubClusts==1 && dyn_msp5c!=-1) {
                            do_sub=1;
                            sub[0]=dyn_up_sp5c[i];
                            sub[1]=dyn_up_sp5c[j];
                            quickSort(&sub[0],2);
                        }
                        else do_sub=0;
                        if (doSubClusts==1) do_up=1;
                        else do_up=0;
                        Dyn_add(hc11C[n11C[f]], f, clusSize, &dyn_n11C, &dyn_m11C, &dyn_l11C, &dyn_hc11C, do_up, dyn_up_11C, n11C[f], do_sub, n_sub, &dyn_sub_11C, sub);
                    }
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
                    
                    if (doBondedCen==1) {
                        n_bonded_to_cen_11C+=cnb[hc11C[n11C[f]][0]];
                        n_distro_bonded_to_cen_11C[cnb[hc11C[n11C[f]][0]]]++;
                    }
                    
                    if (do12A==1) {
                        if(Clusters_Get12A_C2v(f, ach2, ach2_cen, ach2_shell)) {
                            ach2_cen[hc11C[n11C[f]][0]] = 'O';
                            ++n12A[f];
                        }
                    }
                    
                    if (doClusBLDistros==1) {
                        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                                if (Bonds_BondCheck(hc11C[n11C[f]][binAcnt],hc11C[n11C[f]][binBcnt])==1) {
                                    Bonds_TickBLDistro(bondlengths[hc11C[n11C[f]][binAcnt]][Bonds_cnb_j(hc11C[n11C[f]][binAcnt],hc11C[n11C[f]][binBcnt])],BLDistro11C,&BLDistroNoSamples11C);
                                }
                            }
                        }
                    }

                    if (doClusComp==1) {
                        number_of_A=0;
                        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                            if (rtype[hc11C[n11C[f]][binAcnt]]==1) {
                                nA11C++;
                                number_of_A++;
                            }
                            else nB11C++;
                
                            if (binAcnt!=0) {
                                if (rtype[hc11C[n11C[f]][binAcnt]]==1) nA_shell_11C++;
                                else nB_shell_11C++;
                            }
                        }
                        n_distro_11C[number_of_A]++;
                        
                        if (rtype[hc11C[n11C[f]][0]]==1) {
                            nA_cen_11C++;
                            n_distro_cen_11C[1]++;
                            n_distro_shell_11C[number_of_A-1]++;
                        }
                        else {
                            nB_cen_11C++;
                            n_distro_cen_11C[0]++;
                            n_distro_shell_11C[number_of_A]++;
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=12;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=1;
    int sub[1];
    
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
        if (doClusBLDeviation==1) {
            bl_mom_12A=resize_1D_double(bl_mom_12A,m12A,m12A+incrStatic);
        }
        m12A=m12A+incrStatic;
    }
    
    // hc12A key: (as 11C, extra_s)
    
    for(k=0;k<11;k++) hc12A[n12A[f]][k]=hc11C[n11C[f]][k];
    hc12A[n12A[f]][11]=ep;
    
    if (doDynamics==1 && dyn_m12A!=-1) {
        if (doSubClusts==1 && dyn_m11C!=-1) {
            do_sub=1;
            sub[0]=dyn_up_11C[n11C[f]];
        }
        else do_sub=0;
        do_up=0;
        Dyn_add(hc12A[n12A[f]], f, clusSize, &dyn_n12A, &dyn_m12A, &dyn_l12A, &dyn_hc12A, do_up, dummy_up, n12A[f], do_sub, n_sub, &dyn_sub_12A, sub);
    }
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
    
    if (doBondedCen==1) {
        n_bonded_to_cen_12A+=cnb[hc12A[n12A[f]][0]];
        n_distro_bonded_to_cen_12A[cnb[hc12A[n12A[f]][0]]]++;
    }
    
    if (doClusBLDistros==1) {
        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                if (Bonds_BondCheck(hc12A[n12A[f]][binAcnt],hc12A[n12A[f]][binBcnt])==1) {
                    Bonds_TickBLDistro(bondlengths[hc12A[n12A[f]][binAcnt]][Bonds_cnb_j(hc12A[n12A[f]][binAcnt],hc12A[n12A[f]][binBcnt])],BLDistro12A,&BLDistroNoSamples12A);
                }
            }
        }
    }
    
    if (doClusComp==1) {
        number_of_A=0;
        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
            if (rtype[hc12A[n12A[f]][binAcnt]]==1) {
                nA12A++;
                number_of_A++;
            }
            else nB12A++;
    
            if (binAcnt!=0) {
                if (rtype[hc12A[n12A[f]][binAcnt]]==1) nA_shell_12A++;
                else nB_shell_12A++;
            }
        }
        n_distro_12A[number_of_A]++;
        
        if (rtype[hc12A[n12A[f]][0]]==1) {
            nA_cen_12A++;
            n_distro_cen_12A[1]++;
            n_distro_shell_12A[number_of_A-1]++;
        }
        else {
            nB_cen_12A++;
            n_distro_cen_12A[0]++;
            n_distro_shell_12A[number_of_A]++;
        }
    }
    
    return 1;
}

void Clusters_Get11F_12E_13K(int f) {   // Detect 11F C2v & 12E 3h
    char *ach1, *ach2, *ach3, *ach3_cen, *ach3_shell;
    int cp, bpi, bpj, ep1, ep2, the6A_i, the6A_j;
    int i, j, j2, k, l, m;  
    int flg, flg1, flg2;
    int break_out;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=11;
    int do_up=0;
    int do_sub=0;
    int n_sub=4;
    int sub[4];

    cp=bpi=bpj=ep1=ep2=the6A_i=the6A_j=-1;
    ach1=malloc(N*sizeof(char));    if (ach1==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach1[] malloc out of memory\n");   Error(errMsg); }
    ach2=malloc(N*sizeof(char));    if (ach2==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach2[] malloc out of memory\n");   Error(errMsg); }
    ach3=malloc(N*sizeof(char));    if (ach3==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach3[] malloc out of memory\n");   Error(errMsg); }
    ach3_cen=malloc(N*sizeof(char));    if (ach3_cen==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach3_cen[] malloc out of memory\n");   Error(errMsg); }
    ach3_shell=malloc(N*sizeof(char));  if (ach3_shell==NULL) { sprintf(errMsg,"Clusters_Get11F_12E_13K(): ach3_shell[] malloc out of memory\n");   Error(errMsg); }
    for(i=0; i<N; ++i) ach1[i] = ach2[i] = ach3[i] = ach3_cen[i] = ach3_shell[i] ='C';
    
    for(i=0; i<nsp3c[f]-1; i++) {
        for (j2=0; j2<3; j2++) {
        for (j=0; j<nmem_sp3c[sp3c[i][j2]]; ++j) { // loop over all sp3c_j
            if (mem_sp3c[sp3c[i][j2]][j]<=i) continue;
            flg = sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4] || sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
            flg = flg || sp3c[i][3] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4] || sp3c[i][4] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3]; 
            if(flg==1) continue;
            if(Bonds_BondCheck(sp3c[i][3],sp3c[mem_sp3c[sp3c[i][j2]][j]][3])==1 && Bonds_BondCheck(sp3c[i][4], sp3c[mem_sp3c[sp3c[i][j2]][j]][4])==1) {
                m = 0;
                for(k=0; k<3; ++k) { 
                    for(l=0; l<3; ++l){
                        if(sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]){
                            cp = sp3c[i][k];
                            ++m;
                        }
                    }
                }
                
                if(m==1) { // only 1 particle common between the pair
                    m = 0;
                    for(k=0; k<3; ++k) {
                        for(l=0; l<3; ++l) {
                            if(sp3c[i][k] == cp || sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == cp) continue;
                            if(Bonds_BondCheck(sp3c[i][k],sp3c[mem_sp3c[sp3c[i][j2]][j]][l])==1){ 
                                if(m++) break;
                                bpi = sp3c[i][k];
                                bpj = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                            }                               
                        }
                    }
                    
                    if(m==1) { // There must be exactly one paticle of sp3_i bonded to exactly one of sp3_j
                        flg1 = flg2 = 0;
                        for(k=0; k<nsp4c[f]; ++k) {
                            if(sp4c[k][4] == cp || sp4c[k][5] == cp) {
                                if(flg1==0) { // check for first sp4c
                                    for(l=0; l<4; ++l) {
                                        flg = sp4c[k][l] == sp3c[i][3] || sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                                        flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                        if(flg==0) break;
                                    }   
                                    if(l==4) {
                                        flg1 = 1;
                                        if(sp4c[k][4] == cp) ep1 = sp4c[k][5];
                                        else ep1 = sp4c[k][4];
                                        the6A_i=k;
                                    }
                                }
                                if(flg2==0) { // check for first sp4c
                                    for(l=0; l<4; ++l) {
                                        flg = sp4c[k][l] == sp3c[i][4] || sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                                        flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                        if(flg==0) break;
                                    }   
                                    if(l==4) {
                                        flg2 = 1;
                                        if(sp4c[k][4] == cp) ep2 = sp4c[k][5];
                                        else ep2 = sp4c[k][4];
                                        the6A_j=k;
                                    }
                                }                           
                            }
                            if(flg1==1 && flg2==1) break;                       
                        }
                
                        if(k<nsp4c[f]) { // 11F found 
                            if(n11F[f] == m11F) { 
                                hc11F=resize_2D_int(hc11F,m11F,m11F+incrStatic,clusSize,-1);
                                if (doClusBLDeviation==1) {
                                    bl_mom_11F=resize_1D_double(bl_mom_11F,m11F,m11F+incrStatic);
                                }
                                m11F=m11F+incrStatic;
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
                            break_out=0;
                            for(l=0; l<3; ++l) {
                                if(sp3c[i][l] != cp) {
                                    if (m==11) {
                                        break_out=1;
                                        break;
                                    }
                                    hc11F[n11F[f]][m] = sp3c[i][l];
                                    m++;
                                }
                            }
                            for(l=0; l<3; ++l) {
                                if(sp3c[mem_sp3c[sp3c[i][j2]][j]][l] != cp) {
                                    if (m==11) {
                                        break_out=1;
                                        break;
                                    }
                                    hc11F[n11F[f]][m] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                                    m++;
                                }
                            }
                            if(break_out==1 || m<11) continue;
                            
                            quickSort(&hc11F[n11F[f]][1],2);
                            quickSort(&hc11F[n11F[f]][3],4);
                            quickSort(&hc11F[n11F[f]][7],4);
                            
                            if (doDynamics==1 && dyn_m11F!=-1) {
                                if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m6A!=-1) {
                                    do_sub=1;
                                    sub[0]=dyn_up_sp3c[i];
                                    sub[1]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][j]];
                                    sub[2]=dyn_up_sp4c[the6A_i];
                                    sub[3]=dyn_up_sp4c[the6A_j];
                                    quickSort(&sub[0],2);
                                    quickSort(&sub[2],2);
                                }
                                else do_sub=0;
                                if (doSubClusts==1) do_up=1;
                                else do_up=0;
                                Dyn_add(hc11F[n11F[f]], f, clusSize, &dyn_n11F, &dyn_m11F, &dyn_l11F, &dyn_hc11F, do_up, dyn_up_11F, n11F[f], do_sub, n_sub, &dyn_sub_11F, sub);
                            }
                            if(ach1[hc11F[n11F[f]][0]] == 'C') ach1[hc11F[n11F[f]][0]] = 'B';
                            if(ach1[hc11F[n11F[f]][7]] == 'C') ach1[hc11F[n11F[f]][7]] = 'B';
                            if(ach1[hc11F[n11F[f]][8]] == 'C') ach1[hc11F[n11F[f]][8]] = 'B';
                            if(ach1[hc11F[n11F[f]][9]] == 'C') ach1[hc11F[n11F[f]][9]] = 'B';
                            if(ach1[hc11F[n11F[f]][10]] == 'C') ach1[hc11F[n11F[f]][10]] = 'B';
                            ach1[hc11F[n11F[f]][1]] = 'O';
                            ach1[hc11F[n11F[f]][2]] = 'O';
                            ach1[hc11F[n11F[f]][3]] = 'O';
                            ach1[hc11F[n11F[f]][4]] = 'O';  
                            ach1[hc11F[n11F[f]][5]] = 'O';
                            ach1[hc11F[n11F[f]][6]] = 'O';  
                            
                            if (do12E==1) {
                                if (Clusters_Get12E_D3h(f,mem_sp3c[sp3c[i][j2]][j],ach2)) ++n12E[f];
                            }
                            
                            if (do13K==1) {
                                if (Clusters_Get13K(f,i,mem_sp3c[sp3c[i][j2]][j],the6A_i,ach3,ach3_cen,ach3_shell)) ++n13K[f];
                            }                               

                            if (doClusBLDistros==1) {
                                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                                        if (Bonds_BondCheck(hc11F[n11F[f]][binAcnt],hc11F[n11F[f]][binBcnt])==1) {
                                            Bonds_TickBLDistro(bondlengths[hc11F[n11F[f]][binAcnt]][Bonds_cnb_j(hc11F[n11F[f]][binAcnt],hc11F[n11F[f]][binBcnt])],BLDistro11F,&BLDistroNoSamples11F);
                                        }
                                    }
                                }
                            }
                            
                            if (doClusComp==1) {
                                number_of_A=0;
                                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                                    if (rtype[hc11F[n11F[f]][binAcnt]]==1) {
                                        nA11F++;
                                        number_of_A++;
                                    }
                                    else nB11F++;
                                }
                                n_distro_11F[number_of_A]++;
                            }
                            
                            ++n11F[f];      
                        }
                    }
                }
            }
            if(Bonds_BondCheck(sp3c[i][3], sp3c[mem_sp3c[sp3c[i][j2]][j]][4])==1 && Bonds_BondCheck(sp3c[i][4], sp3c[mem_sp3c[sp3c[i][j2]][j]][3])==1) {
                m = 0;
                for(k=0; k<3; ++k) { 
                    for(l=0; l<3; ++l){
                        if(sp3c[i][k] == sp3c[mem_sp3c[sp3c[i][j2]][j]][l]){
                            cp = sp3c[i][k];
                            ++m;
                        }
                    }
                }
                
                if(m==1) { // only 1 particle common between the pair
                    m = 0;
                    for(k=0; k<3; ++k) {
                        for(l=0; l<3; ++l) {
                            if(sp3c[i][k] == cp || sp3c[mem_sp3c[sp3c[i][j2]][j]][l] == cp) continue;
                            if(Bonds_BondCheck(sp3c[i][k], sp3c[mem_sp3c[sp3c[i][j2]][j]][l])==1) { 
                                if(m++) break;
                                bpi = sp3c[i][k];
                                bpj = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];          
                            }
                        }
                    }
                    
                    if(m==1) { // The single bonded particle found
                        flg1 = flg2 = 0;
                        for(k=0; k<nsp4c[f]; ++k) {
                            if(sp4c[k][4] == cp || sp4c[k][5] == cp) {
                                if(flg1==0) { // check for first sp4c
                                    for(l=0; l<4; ++l) {
                                        flg = sp4c[k][l] == sp3c[i][3] || sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][4];
                                        flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                        if(flg==0) break;
                                    }   
                                    if(l==4) {
                                        flg1 = 1;
                                        if(sp4c[k][4] == cp) ep1 = sp4c[k][5];
                                        else ep1 = sp4c[k][4];
                                        the6A_i=k;
                                    }
                                }
                                if(flg2==0) { // check for first sp4c
                                    for(l=0; l<4; ++l) {
                                        flg = sp4c[k][l] == sp3c[i][4] || sp4c[k][l] == sp3c[mem_sp3c[sp3c[i][j2]][j]][3];
                                        flg = flg || sp4c[k][l] == bpi || sp4c[k][l] == bpj;
                                        if(flg==0) break;
                                    }   
                                    if(l==4) {
                                        flg2 = 1;
                                        if(sp4c[k][4] == cp) ep2 = sp4c[k][5];
                                        else ep2 = sp4c[k][4];
                                        the6A_j=k;
                                    }
                                }                           
                            }
                            if(flg1==1 && flg2==1) break;                       
                        }
                        if(k<nsp4c[f]) { // 11F found 
                            if(n11F[f] == m11F) { 
                                hc11F=resize_2D_int(hc11F,m11F,m11F+incrStatic,clusSize,-1);
                                if (doClusBLDeviation==1) {
                                    bl_mom_11F=resize_1D_double(bl_mom_11F,m11F,m11F+incrStatic);
                                }
                                m11F=m11F+incrStatic;
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
                            break_out=0;
                            for(l=0; l<3; ++l) {
                                if(sp3c[i][l] != cp) {
                                    if (m==11) {
                                        break_out=1;
                                        break;
                                    }
                                    hc11F[n11F[f]][m] = sp3c[i][l];
                                    m++;
                                }
                            }
                            for(l=0; l<3; ++l) {
                                if(sp3c[mem_sp3c[sp3c[i][j2]][j]][l] != cp) {
                                    if (m==11) {
                                        break_out=1;
                                        break;
                                    }
                                    hc11F[n11F[f]][m] = sp3c[mem_sp3c[sp3c[i][j2]][j]][l];
                                    m++;
                                }
                            }
                            if(break_out==1 || m<11) continue;
                            
                            quickSort(&hc11F[n11F[f]][1],2);
                            quickSort(&hc11F[n11F[f]][3],4);
                            quickSort(&hc11F[n11F[f]][7],4);
                            
                            if (doDynamics==1 && dyn_m11F!=-1) {
                                if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m6A!=-1) {
                                    do_sub=1;
                                    sub[0]=dyn_up_sp3c[i];
                                    sub[1]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][j]];
                                    sub[2]=dyn_up_sp4c[the6A_i];
                                    sub[3]=dyn_up_sp4c[the6A_j];
                                    quickSort(&sub[0],2);
                                    quickSort(&sub[2],2);
                                }
                                else do_sub=0;
                                if (doSubClusts==1) do_up=1;
                                else do_up=0;
                                Dyn_add(hc11F[n11F[f]], f, clusSize, &dyn_n11F, &dyn_m11F, &dyn_l11F, &dyn_hc11F, do_up, dyn_up_11F, n11F[f], do_sub, n_sub, &dyn_sub_11F, sub);
                            }
                            if(ach1[hc11F[n11F[f]][0]] == 'C') ach1[hc11F[n11F[f]][0]] = 'B';
                            if(ach1[hc11F[n11F[f]][7]] == 'C') ach1[hc11F[n11F[f]][7]] = 'B';
                            if(ach1[hc11F[n11F[f]][8]] == 'C') ach1[hc11F[n11F[f]][8]] = 'B';
                            if(ach1[hc11F[n11F[f]][9]] == 'C') ach1[hc11F[n11F[f]][9]] = 'B';
                            if(ach1[hc11F[n11F[f]][10]] == 'C') ach1[hc11F[n11F[f]][10]] = 'B';
                            ach1[hc11F[n11F[f]][1]] = 'O';
                            ach1[hc11F[n11F[f]][2]] = 'O';
                            ach1[hc11F[n11F[f]][3]] = 'O';
                            ach1[hc11F[n11F[f]][4]] = 'O';  
                            ach1[hc11F[n11F[f]][5]] = 'O';
                            ach1[hc11F[n11F[f]][6]] = 'O';  
                            
                            if (do12K==1) {
                                if (Clusters_Get12E_D3h(f,mem_sp3c[sp3c[i][j2]][j],ach2)) ++n12E[f];
                            }
                            
                            if (do13K==1) {
                                if (Clusters_Get13K(f,i,mem_sp3c[sp3c[i][j2]][j],the6A_i,ach3,ach3_cen,ach3_shell)) ++n13K[f];
                            }
                            
                            if (doClusBLDistros==1) {
                                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                                        if (Bonds_BondCheck(hc11F[n11F[f]][binAcnt],hc11F[n11F[f]][binBcnt])==1) {
                                            Bonds_TickBLDistro(bondlengths[hc11F[n11F[f]][binAcnt]][Bonds_cnb_j(hc11F[n11F[f]][binAcnt],hc11F[n11F[f]][binBcnt])],BLDistro11F,&BLDistroNoSamples11F);
                                        }
                                    }
                                }
                            }
                        
                            if (doClusComp==1) {
                                number_of_A=0;
                                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                                    if (rtype[hc11F[n11F[f]][binAcnt]]==1) {
                                        nA11F++;
                                        number_of_A++;
                                    }
                                    else nB11F++;
                                }
                                n_distro_11F[number_of_A]++;
                            }
                            
                            ++n11F[f];      
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=12;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_12E=resize_1D_double(bl_mom_12E,m12E,m12E+incrStatic);
                }
                m12E=m12E+incrStatic;
            }
            
            hc12E[n12E[f]][11] = uncom;
            for(l=0; l<10; ++l) hc12E[n12E[f]][l] = hc11F[n11F[f]][l+1];
            hc12E[n12E[f]][10] = hc11F[n11F[f]][0];
            quickSort(&hc12E[n12E[f]][0],6);
            quickSort(&hc12E[n12E[f]][6],6);
            
            if (doDynamics==1 && dyn_m12E!=-1) {
                if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m11F!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_11F[n11F[f]];
                    sub[1]=dyn_up_sp3c[k];
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc12E[n12E[f]], f, clusSize, &dyn_n12E, &dyn_m12E, &dyn_l12E, &dyn_hc12E, do_up, dummy_up, n12E[f], do_sub, n_sub, &dyn_sub_12E, sub);
            }
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
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc12E[n12E[f]][binAcnt],hc12E[n12E[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc12E[n12E[f]][binAcnt]][Bonds_cnb_j(hc12E[n12E[f]][binAcnt],hc12E[n12E[f]][binBcnt])],BLDistro12E,&BLDistroNoSamples12E);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc12E[n12E[f]][binAcnt]]==1) {
                        nA12E++;
                        number_of_A++;
                    }
                    else nB12E++;
                }
                n_distro_12E[number_of_A]++;
            }
            /*if (hc12E[n12E[f]][0]==325) {
                printf("hi i%d j%d k%d n11F[f]%d\n",i,j,k,n11F[f]);
                printf("sp3c[ %d ] %d %d %d %d %d\n",i,sp3c[i][0],sp3c[i][1],sp3c[i][2],sp3c[i][3],sp3c[i][4]);
                printf("sp3c[ %d ] %d %d %d %d %d\n",j,sp3c[j][0],sp3c[j][1],sp3c[j][2],sp3c[j][3],sp3c[j][4]);
                printf("sp3c[ %d ] %d %d %d %d %d\n",k,sp3c[k][0],sp3c[k][1],sp3c[k][2],sp3c[k][3],sp3c[k][4]);
            }*/
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=13;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=3;
    int sub[3];
    
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
        if (doClusBLDeviation==1) {
            bl_mom_13K=resize_1D_double(bl_mom_13K,m13K,m13K+incrStatic);
        }
        m13K=m13K+incrStatic;
    }
    // hc13K key: (11F, extra SP3 ring particle to make 5A #1, extra SP3 ring particle to make 5A #2)

    for (i=0; i<11; i++) hc13K[n13K[f]][i] = hc11F[n11F[f]][i];
    hc13K[n13K[f]][11]=ep[0];
    hc13K[n13K[f]][12]=ep[1];
    
    quickSort(&hc13K[n13K[f]][11],2);
    quickSort(&eclus5A[0],2);

    if (doDynamics==1 && dyn_m13K!=-1) {
        if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m11F!=-1) {
            do_sub=1;
            sub[0]=dyn_up_11F[n11F[f]];
            sub[1]=dyn_up_sp3c[eclus5A[0]];
            sub[2]=dyn_up_sp3c[eclus5A[1]];
            quickSort(&sub[1],2);
        }
        else do_sub=0;
        do_up=0;
        Dyn_add(hc13K[n13K[f]], f, clusSize, &dyn_n13K, &dyn_m13K, &dyn_l13K, &dyn_hc13K, do_up, dummy_up, n13K[f], do_sub, n_sub, &dyn_sub_13K, sub);
    }
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
    
    if (doBondedCen==1) {
        n_bonded_to_cen_13K+=cnb[hc13K[n13K[f]][0]];
        n_distro_bonded_to_cen_13K[cnb[hc13K[n13K[f]][0]]]++;
    }
                
    if (doClusBLDistros==1) {
        for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
            for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                if (Bonds_BondCheck(hc13K[n13K[f]][binAcnt],hc13K[n13K[f]][binBcnt])==1) {
                    Bonds_TickBLDistro(bondlengths[hc13K[n13K[f]][binAcnt]][Bonds_cnb_j(hc13K[n13K[f]][binAcnt],hc13K[n13K[f]][binBcnt])],BLDistro13K,&BLDistroNoSamples13K);
                }
            }
        }
    }
    
    if (doClusComp==1) {
        number_of_A=0;
        for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
            if (rtype[hc13K[n13K[f]][binAcnt]]==1) {
                nA13K++;
                number_of_A++;
            }
            else nB13K++;
            
            if (binAcnt!=0) {
                if (rtype[hc13K[n13K[f]][binAcnt]]==1) nA_shell_13K++;
                else nB_shell_13K++;
            }
        }
        n_distro_13K[number_of_A]++;
        
        if (rtype[hc13K[n13K[f]][0]]==1) {
            nA_cen_13K++;
            n_distro_cen_13K[1]++;
            n_distro_shell_13K[number_of_A-1]++;
        }
        else {
            nB_cen_13K++;
            n_distro_cen_13K[0]++;
            n_distro_shell_13K[number_of_A]++;
        }
    }
    
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
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=12;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=6;
    int sub[6];
    int **dummy_sub=NULL;
    
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
                    if (doClusBLDeviation==1) {
                        bl_mom_13A=resize_1D_double(bl_mom_13A,m13A,m13A+incrStatic);
                    }
                    m13A=m13A+incrStatic;
                }
                
                hc13A[n13A[f]][0] = sp1; 
                k = 1;
                if(sp5c[i][5] != sp1) hc13A[n13A[f]][k++] = sp5c[i][5]; 
                if(sp5c[i][6] != sp1) hc13A[n13A[f]][k++] = sp5c[i][6]; 
                if(sp5c[j][5] != sp1) hc13A[n13A[f]][k++] = sp5c[j][5]; 
                if(sp5c[j][6] != sp1) hc13A[n13A[f]][k++] = sp5c[j][6];  
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
                
                if (doDynamics==1 && dyn_m13A!=-1) {
                    do_sub=0;
                    do_up=0;
                    Dyn_add(hc13A[n13A[f]], f, clusSize+1, &dyn_n13A, &dyn_m13A, &dyn_l13A, &dyn_hc13A, do_up, dummy_up, n13A[f], do_sub, 0, &dummy_sub, sub);
                }
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
                
                if (doBondedCen==1) {
                    n_bonded_to_cen_13A+=cnb[hc13A[n13A[f]][0]];
                    n_distro_bonded_to_cen_13A[cnb[hc13A[n13A[f]][0]]]++;
                }
    
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize+1; binBcnt++) {
                            if (Bonds_BondCheck(hc13A[n13A[f]][binAcnt],hc13A[n13A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc13A[n13A[f]][binAcnt]][Bonds_cnb_j(hc13A[n13A[f]][binAcnt],hc13A[n13A[f]][binBcnt])],BLDistro13A,&BLDistroNoSamples13A);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize+1; binAcnt++) {
                        if (rtype[hc13A[n13A[f]][binAcnt]]==1) {
                            nA13A++;
                            number_of_A++;
                        }
                        else nB13A++;
                    
                        if (binAcnt!=0) {
                            if (rtype[hc13A[n13A[f]][binAcnt]]==1) nA_shell_13A++;
                            else nB_shell_13A++;
                        }
                    }
                    n_distro_13A[number_of_A]++;
                    
                    if (rtype[hc13A[n13A[f]][0]]==1) {
                        nA_cen_13A++;
                        n_distro_cen_13A[1]++;
                        n_distro_shell_13A[number_of_A-1]++;
                    }
                    else {
                        nB_cen_13A++;
                        n_distro_cen_13A[0]++;
                        n_distro_shell_13A[number_of_A]++;
                    }
                }
                
                ++n13A[f];
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
                    if (doClusBLDeviation==1) {
                        bl_mom_13A=resize_1D_double(bl_mom_13A,m13A,m13A+incrStatic);
                    }
                    m13A=m13A+incrStatic;
                }
                
                hc13A[n13A[f]][0] = sp2; 
                k = 1;
                if(sp5c[i][5] != sp2) hc13A[n13A[f]][k++] = sp5c[i][5]; 
                if(sp5c[i][6] != sp2) hc13A[n13A[f]][k++] = sp5c[i][6]; 
                if(sp5c[j][5] != sp2) hc13A[n13A[f]][k++] = sp5c[j][5]; 
                if(sp5c[j][6] != sp2) hc13A[n13A[f]][k++] = sp5c[j][6];  
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
                
                if (doDynamics==1 && dyn_m13A!=-1) {
                    do_sub=0;
                    do_up=0;
                    Dyn_add(hc13A[n13A[f]], f, clusSize+1, &dyn_n13A, &dyn_m13A, &dyn_l13A, &dyn_hc13A, do_up, dummy_up, n13A[f], do_sub, 0, &dummy_sub, sub);
                }
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
                
                if (doBondedCen==1) {
                    n_bonded_to_cen_13A+=cnb[hc13A[n13A[f]][0]];
                    n_distro_bonded_to_cen_13A[cnb[hc13A[n13A[f]][0]]]++;
                }
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize+1; binBcnt++) {
                            if (Bonds_BondCheck(hc13A[n13A[f]][binAcnt],hc13A[n13A[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc13A[n13A[f]][binAcnt]][Bonds_cnb_j(hc13A[n13A[f]][binAcnt],hc13A[n13A[f]][binBcnt])],BLDistro13A,&BLDistroNoSamples13A);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize+1; binAcnt++) {
                        if (rtype[hc13A[n13A[f]][binAcnt]]==1) {
                            nA13A++;
                            number_of_A++;
                        }
                        else nB13A++;
                        
                        if (binAcnt!=0) {
                            if (rtype[hc13A[n13A[f]][binAcnt]]==1) nA_shell_13A++;
                            else nB_shell_13A++;
                        }
                    }
                    n_distro_13A[number_of_A]++;
                    
                    if (rtype[hc13A[n13A[f]][0]]==1) {
                        nA_cen_13A++;
                        n_distro_cen_13A[1]++;
                        n_distro_shell_13A[number_of_A-1]++;
                    }
                    else {
                        nB_cen_13A++;
                        n_distro_cen_13A[0]++;
                        n_distro_shell_13A[number_of_A]++;
                    }
                }
                
                ++n13A[f];
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
                if (doClusBLDeviation==1) {
                    bl_mom_12B=resize_1D_double(bl_mom_12B,m12B,m12B+incrStatic);
                }
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
            
            quickSort(&hc12B[n12B[f]][2],5);
            quickSort(&hc12B[n12B[f]][7],5);
            
            if (doDynamics==1 && dyn_m12B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                    sub[1]=dyn_up_sp5c[sj1[0]];
                    sub[2]=dyn_up_sp5c[sj1[1]];
                    sub[3]=dyn_up_sp5c[sj1[2]];
                    sub[4]=dyn_up_sp5c[sj1[3]];
                    sub[5]=dyn_up_sp5c[sj1[4]];
                    quickSort(&sub[0],6);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc12B[n12B[f]], f, clusSize, &dyn_n12B, &dyn_m12B, &dyn_l12B, &dyn_hc12B, do_up, dummy_up, n12B[f], do_sub, n_sub, &dyn_sub_12B, sub);
            }
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
            
            if (doBondedCen==1) {
                n_bonded_to_cen_12B+=cnb[hc12B[n12B[f]][0]];
                n_distro_bonded_to_cen_12B[cnb[hc12B[n12B[f]][0]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc12B[n12B[f]][binAcnt],hc12B[n12B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc12B[n12B[f]][binAcnt]][Bonds_cnb_j(hc12B[n12B[f]][binAcnt],hc12B[n12B[f]][binBcnt])],BLDistro12B,&BLDistroNoSamples12B);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc12B[n12B[f]][binAcnt]]==1) {
                        nA12B++;
                        number_of_A++;
                    }
                    else nB12B++;
                
                    if (binAcnt!=0) {
                        if (rtype[hc12B[n12B[f]][binAcnt]]==1) nA_shell_12B++;
                        else nB_shell_12B++;
                    }   
                }
                n_distro_12B[number_of_A]++;
                
                if (rtype[hc12B[n12B[f]][0]]==1) {
                    nA_cen_12B++;
                    n_distro_cen_12B[1]++;
                    n_distro_shell_12B[number_of_A-1]++;
                }
                else {
                    nB_cen_12B++;
                    n_distro_cen_12B[0]++;
                    n_distro_shell_12B[number_of_A]++;
                }
            }
                
            ++n12B[f];
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
                if (doClusBLDeviation==1) {
                    bl_mom_12B=resize_1D_double(bl_mom_12B,m12B,m12B+incrStatic);
                }
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
    
            quickSort(&hc12B[n12B[f]][2],5);
            quickSort(&hc12B[n12B[f]][7],5);
            
            if (doDynamics==1 && dyn_m12B!=-1) {
                if (doSubClusts==1 && dyn_msp5c!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp5c[i];
                    sub[1]=dyn_up_sp5c[sj2[0]];
                    sub[2]=dyn_up_sp5c[sj2[1]];
                    sub[3]=dyn_up_sp5c[sj2[2]];
                    sub[4]=dyn_up_sp5c[sj2[3]];
                    sub[5]=dyn_up_sp5c[sj2[4]];
                    quickSort(&sub[0],6);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hc12B[n12B[f]], f, clusSize, &dyn_n12B, &dyn_m12B, &dyn_l12B, &dyn_hc12B, do_up, dummy_up, n12B[f], do_sub, n_sub, &dyn_sub_12B, sub);
            }
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
            
            if (doBondedCen==1) {
                n_bonded_to_cen_12B+=cnb[hc12B[n12B[f]][0]];
                n_distro_bonded_to_cen_12B[cnb[hc12B[n12B[f]][0]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hc12B[n12B[f]][binAcnt],hc12B[n12B[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hc12B[n12B[f]][binAcnt]][Bonds_cnb_j(hc12B[n12B[f]][binAcnt],hc12B[n12B[f]][binBcnt])],BLDistro12B,&BLDistroNoSamples12B);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hc12B[n12B[f]][binAcnt]]==1) {
                        nA12B++;
                        number_of_A++;
                    }
                    else nB12B++;
                
                    if (binAcnt!=0) {
                        if (rtype[hc12B[n12B[f]][binAcnt]]==1) nA_shell_12B++;
                        else nB_shell_12B++;
                    }
                }
                n_distro_12B[number_of_A]++;
            
                if (rtype[hc12B[n12B[f]][0]]==1) {
                    nA_cen_12B++;
                    n_distro_cen_12B[1]++;
                    n_distro_shell_12B[number_of_A-1]++;
                }
                else {
                    nB_cen_12B++;
                    n_distro_cen_12B[0]++;
                    n_distro_shell_12B[number_of_A]++;
                }
            }
            
            ++n12B[f];
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

void Clusters_Get13B_D5h(int f) {   // Detect 13B D5h clusters, i.e. twisted icosahedra
    int cp;
    int i, j, k, l, m; 
    int flg;
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=13;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=2;
    int sub[2];
    
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
                    if (doClusBLDeviation==1) {
                        bl_mom_13B=resize_1D_double(bl_mom_13B,m13B,m13B+incrStatic);
                    }
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
                
                if (doDynamics==1 && dyn_m13B!=-1) {
                    if (doSubClusts==1 && dyn_msp5c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp5c[i];
                        sub[1]=dyn_up_sp5c[j];
                        quickSort(&sub[0],2);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hc13B[n13B[f]], f, clusSize, &dyn_n13B, &dyn_m13B, &dyn_l13B, &dyn_hc13B, do_up, dummy_up, n13B[f], do_sub, n_sub, &dyn_sub_13B, sub);
                }
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
                
                if (doBondedCen==1) {
                    n_bonded_to_cen_13B+=cnb[hc13B[n13B[f]][0]];
                    n_distro_bonded_to_cen_13B[cnb[hc13B[n13B[f]][0]]]++;
                }
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hc13B[n13B[f]][binAcnt],hc13B[n13B[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hc13B[n13B[f]][binAcnt]][Bonds_cnb_j(hc13B[n13B[f]][binAcnt],hc13B[n13B[f]][binBcnt])],BLDistro13B,&BLDistroNoSamples13B);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hc13B[n13B[f]][binAcnt]]==1) {
                            nA13B++;
                            number_of_A++;
                        }
                        else nB13B++;
                    
                        if (binAcnt!=0) {
                            if (rtype[hc13B[n13B[f]][binAcnt]]==1) nA_shell_13B++;
                            else nB_shell_13B++;
                        }
                    }
                    n_distro_13B[number_of_A]++;
                
                    if (rtype[hc13B[n13B[f]][0]]==1) {
                        nA_cen_13B++;
                        n_distro_cen_13B[1]++;
                        n_distro_shell_13B[number_of_A-1]++;
                    }
                    else {
                        nB_cen_13B++;
                        n_distro_cen_13B[0]++;
                        n_distro_shell_13B[number_of_A]++;
                    }
                }
                
                ++n13B[f];
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

void Clusters_GetFCC(int f) {   // Detect 13 particle FCC clusters
    int i, j, j2, k, l, m, n; 
    int i1, i2, i3;
    int cp, bpi, bpj, nbpi, nbpj;
    int flg1, flg2, flg3;
    int l_clust_type; // 0 if l-clust is sp3b, 1 if is sp3c
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=13;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=5;
    int sub[5];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_FCC=resize_1D_double(bl_mom_FCC,mFCC,mFCC+incrStatic);
                }
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
                    if (doClusBLDeviation==1) {
                        bl_mom_FCC=resize_1D_double(bl_mom_FCC,mFCC,mFCC+incrStatic);
                    }
                    mFCC=mFCC+incrStatic;
                }
        
                hcFCC[nFCC[f]][7] = sp3b[mem_sp3b[sp3b[i][j2]][k]][3];
                hcFCC[nFCC[f]][8] = i3;
                hcFCC[nFCC[f]][9] = i2;
                quickSort(&hcFCC[nFCC[f]][1],12);
                
                if (doDynamics==1 && dyn_mFCC!=-1) {
                    if (doSubClusts==1 && dyn_msp3b!=-1 && dyn_msp3c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp3b[i];
                        sub[1]=dyn_up_sp3b[mem_sp3b[sp3b[i][j2]][j]];
                        sub[2]=dyn_up_sp3b[mem_sp3b[sp3b[i][j2]][k]];
                        if (l_clust_type==0) {
                            sub[3]=dyn_up_sp3b[l];
                            sub[4]=-1;
                        }
                        else {
                            sub[3]=-1;
                            sub[4]=dyn_up_sp3c[l];
                        }
                        quickSort(&sub[0],3);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hcFCC[nFCC[f]], f, clusSize, &dyn_nFCC, &dyn_mFCC, &dyn_lFCC, &dyn_hcFCC, do_up, dummy_up, nFCC[f], do_sub, n_sub, &dyn_sub_FCC, sub);
                }
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
                
                if (doBondedCen==1) {
                    n_bonded_to_cen_FCC+=cnb[hcFCC[nFCC[f]][0]];
                    n_distro_bonded_to_cen_FCC[cnb[hcFCC[nFCC[f]][0]]]++;
                }
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hcFCC[nFCC[f]][binAcnt],hcFCC[nFCC[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hcFCC[nFCC[f]][binAcnt]][Bonds_cnb_j(hcFCC[nFCC[f]][binAcnt],hcFCC[nFCC[f]][binBcnt])],BLDistroFCC,&BLDistroNoSamplesFCC);
                            }
                        }
                    }
                }

                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hcFCC[nFCC[f]][binAcnt]]==1) {
                            nAFCC++;
                            number_of_A++;
                        }
                        else nBFCC++;
                    
                        if (binAcnt!=0) {
                            if (rtype[hcFCC[nFCC[f]][binAcnt]]==1) nA_shell_FCC++;
                            else nB_shell_FCC++;
                        }
                    }
                    n_distro_FCC[number_of_A]++;
                    
                    if (rtype[hcFCC[nFCC[f]][0]]==1) {
                        nA_cen_FCC++;
                        n_distro_cen_FCC[1]++;
                        n_distro_shell_FCC[number_of_A-1]++;
                    }
                    else {
                        nB_cen_FCC++;
                        n_distro_cen_FCC[0]++;
                        n_distro_shell_FCC[number_of_A]++;
                    }
                }

                ++nFCC[f];
            }  // End k loop
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

void Clusters_GetHCP(int f) {   // Detect 13 particle HCP clusters
    int i, j, j2, k, l, m, n; 
    int ia[2], ja[2], ka[2];
    int cp, x;
    int h1i, h1j, h2i, h2j;
    int flg1, flg2, flg3;
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=13;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=3;
    int sub[3];
    
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
                    if (doClusBLDeviation==1) {
                        bl_mom_HCP=resize_1D_double(bl_mom_HCP,mHCP,mHCP+incrStatic);
                    }
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
                hcHCP[nHCP[f]][l++] = sp3c[mem_sp3c[sp3c[i][j2]][k]][4];
                quickSort(&hcHCP[nHCP[f]][1],6);
                quickSort(&hcHCP[nHCP[f]][7],6);
    
                if (doDynamics==1 && dyn_mHCP!=-1) {
                    if (doSubClusts==1 && dyn_msp3c!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp3c[i];
                        sub[1]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][j]];
                        sub[2]=dyn_up_sp3c[mem_sp3c[sp3c[i][j2]][k]];
                        quickSort(&sub[0],n_sub);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hcHCP[nHCP[f]], f, clusSize, &dyn_nHCP, &dyn_mHCP, &dyn_lHCP, &dyn_hcHCP, do_up, dummy_up, nHCP[f], do_sub, n_sub, &dyn_sub_HCP, sub);
                }
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
                
                if (doBondedCen==1) {
                    n_bonded_to_cen_HCP+=cnb[hcHCP[nHCP[f]][0]];
                    n_distro_bonded_to_cen_HCP[cnb[hcHCP[nHCP[f]][0]]]++;
                }
                
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hcHCP[nHCP[f]][binAcnt],hcHCP[nHCP[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hcHCP[nHCP[f]][binAcnt]][Bonds_cnb_j(hcHCP[nHCP[f]][binAcnt],hcHCP[nHCP[f]][binBcnt])],BLDistroHCP,&BLDistroNoSamplesHCP);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hcHCP[nHCP[f]][binAcnt]]==1) {
                            nAHCP++;
                            number_of_A++;
                        }
                        else nBHCP++;
                    
                        if (binAcnt!=0) {
                            if (rtype[hcHCP[nHCP[f]][binAcnt]]==1) nA_shell_HCP++;
                            else nB_shell_HCP++;
                        }
                    }
                    n_distro_HCP[number_of_A]++;
                    
                    if (rtype[hcHCP[nHCP[f]][0]]==1) {
                        nA_cen_HCP++;
                        n_distro_cen_HCP[1]++;
                        n_distro_shell_HCP[number_of_A-1]++;
                    }
                    else {
                        nB_cen_HCP++;
                        n_distro_cen_HCP[0]++;
                        n_distro_shell_HCP[number_of_A]++;
                    }
                }
                
                ++nHCP[f];
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

void Clusters_GetBCC_9(int f) {
    int i, j, j2, k, l, m;
    int flg;
    int s_com=-1;
    int trial[9];
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=9;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=6;
    int sub[6];
    
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
                if (doClusBLDeviation==1) {
                    bl_mom_BCC_9=resize_1D_double(bl_mom_BCC_9,mBCC_9,mBCC_9+incrStatic);
                }
                mBCC_9=mBCC_9+incrStatic;
            }
            
            for (k=0; k<9; ++k) hcBCC_9[nBCC_9[f]][k]=trial[k];
            
            if (doDynamics==1 && dyn_mBCC_9!=-1) {
                if (doSubClusts==1 && dyn_msp4b!=-1 && dyn_m6A!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp4b[i];
                    sub[1]=dyn_up_sp4b[mem_sp4b[sp4b[i][j2]][j]];
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=-1;
                    sub[5]=-1;
                    quickSort(&sub[0],2);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hcBCC_9[nBCC_9[f]], f, clusSize, &dyn_nBCC_9, &dyn_mBCC_9, &dyn_lBCC_9, &dyn_hcBCC_9, do_up, dummy_up, nBCC_9[f], do_sub, n_sub, &dyn_sub_BCC_9, sub);
            }
            if(ach[hcBCC_9[nBCC_9[f]][1]] == 'C') ach[hcBCC_9[nBCC_9[f]][1]] = ach_shell[hcBCC_9[nBCC_9[f]][1]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][2]] == 'C') ach[hcBCC_9[nBCC_9[f]][2]] = ach_shell[hcBCC_9[nBCC_9[f]][2]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][3]] == 'C') ach[hcBCC_9[nBCC_9[f]][3]] = ach_shell[hcBCC_9[nBCC_9[f]][3]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][4]] == 'C') ach[hcBCC_9[nBCC_9[f]][4]] = ach_shell[hcBCC_9[nBCC_9[f]][4]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][5]] == 'C') ach[hcBCC_9[nBCC_9[f]][5]] = ach_shell[hcBCC_9[nBCC_9[f]][5]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][6]] == 'C') ach[hcBCC_9[nBCC_9[f]][6]] = ach_shell[hcBCC_9[nBCC_9[f]][6]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][7]] == 'C') ach[hcBCC_9[nBCC_9[f]][7]] = ach_shell[hcBCC_9[nBCC_9[f]][7]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][8]] == 'C') ach[hcBCC_9[nBCC_9[f]][8]] = ach_shell[hcBCC_9[nBCC_9[f]][8]] = 'B';
            ach[hcBCC_9[nBCC_9[f]][0]] = ach_cen[hcBCC_9[nBCC_9[f]][0]] = 'F';
            
            if (doBondedCen==1) {
                n_bonded_to_cen_BCC_9+=cnb[hcBCC_9[nBCC_9[f]][0]];
                n_distro_bonded_to_cen_BCC_9[cnb[hcBCC_9[nBCC_9[f]][0]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hcBCC_9[nBCC_9[f]][binAcnt],hcBCC_9[nBCC_9[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hcBCC_9[nBCC_9[f]][binAcnt]][Bonds_cnb_j(hcBCC_9[nBCC_9[f]][binAcnt],hcBCC_9[nBCC_9[f]][binBcnt])],BLDistroBCC_9,&BLDistroNoSamplesBCC_9);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hcBCC_9[nBCC_9[f]][binAcnt]]==1) {
                        nABCC_9++;
                        number_of_A++;
                    }
                    else nBBCC_9++;
                
                    if (binAcnt!=0) {
                        if (rtype[hcBCC_9[nBCC_9[f]][binAcnt]]==1) nA_shell_BCC_9++;
                        else nB_shell_BCC_9++;
                    }
                }
                n_distro_BCC_9[number_of_A]++;
                
                if (rtype[hcBCC_9[nBCC_9[f]][0]]==1) {
                    nA_cen_BCC_9++;
                    n_distro_cen_BCC_9[1]++;
                    n_distro_shell_BCC_9[number_of_A-1]++;
                }
                else {
                    nB_cen_BCC_9++;
                    n_distro_cen_BCC_9[0]++;
                    n_distro_shell_BCC_9[number_of_A]++;
                }
            }
            
            ++nBCC_9[f];
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
                if (doClusBLDeviation==1) {
                    bl_mom_BCC_9=resize_1D_double(bl_mom_BCC_9,mBCC_9,mBCC_9+incrStatic);
                }
                mBCC_9=mBCC_9+incrStatic;
            }
            for (k=0; k<9; ++k) hcBCC_9[nBCC_9[f]][k]=trial[k];
            
            if (doDynamics==1 && dyn_mBCC_9!=-1) {
                if (doSubClusts==1 && dyn_msp4b!=-1 && dyn_m6A!=-1) {
                    do_sub=1;
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=dyn_up_sp4c[i];
                    sub[3]=dyn_up_sp4c[mem_sp4c[sp4c[i][j2]][j]];
                    sub[4]=-1;
                    sub[5]=-1;
                    quickSort(&sub[2],2);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hcBCC_9[nBCC_9[f]], f, clusSize, &dyn_nBCC_9, &dyn_mBCC_9, &dyn_lBCC_9, &dyn_hcBCC_9, do_up, dummy_up, nBCC_9[f], do_sub, n_sub, &dyn_sub_BCC_9, sub);
            }
            if(ach[hcBCC_9[nBCC_9[f]][1]] == 'C') ach[hcBCC_9[nBCC_9[f]][1]] = ach_shell[hcBCC_9[nBCC_9[f]][1]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][2]] == 'C') ach[hcBCC_9[nBCC_9[f]][2]] = ach_shell[hcBCC_9[nBCC_9[f]][2]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][3]] == 'C') ach[hcBCC_9[nBCC_9[f]][3]] = ach_shell[hcBCC_9[nBCC_9[f]][3]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][4]] == 'C') ach[hcBCC_9[nBCC_9[f]][4]] = ach_shell[hcBCC_9[nBCC_9[f]][4]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][5]] == 'C') ach[hcBCC_9[nBCC_9[f]][5]] = ach_shell[hcBCC_9[nBCC_9[f]][5]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][6]] == 'C') ach[hcBCC_9[nBCC_9[f]][6]] = ach_shell[hcBCC_9[nBCC_9[f]][6]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][7]] == 'C') ach[hcBCC_9[nBCC_9[f]][7]] = ach_shell[hcBCC_9[nBCC_9[f]][7]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][8]] == 'C') ach[hcBCC_9[nBCC_9[f]][8]] = ach_shell[hcBCC_9[nBCC_9[f]][8]] = 'B';
            ach[hcBCC_9[nBCC_9[f]][0]] = ach_cen[hcBCC_9[nBCC_9[f]][0]] = 'F';
            
            if (doBondedCen==1) {
                n_bonded_to_cen_BCC_9+=cnb[hcBCC_9[nBCC_9[f]][0]];
                n_distro_bonded_to_cen_BCC_9[cnb[hcBCC_9[nBCC_9[f]][0]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hcBCC_9[nBCC_9[f]][binAcnt],hcBCC_9[nBCC_9[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hcBCC_9[nBCC_9[f]][binAcnt]][Bonds_cnb_j(hcBCC_9[nBCC_9[f]][binAcnt],hcBCC_9[nBCC_9[f]][binBcnt])],BLDistroBCC_9,&BLDistroNoSamplesBCC_9);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hcBCC_9[nBCC_9[f]][binAcnt]]==1) {
                        nABCC_9++;
                        number_of_A++;
                    }
                    else nBBCC_9++;
                
                    if (binAcnt!=0) {
                        if (rtype[hcBCC_9[nBCC_9[f]][binAcnt]]==1) nA_shell_BCC_9++;
                        else nB_shell_BCC_9++;
                    }
                }
                n_distro_BCC_9[number_of_A]++;
                
                if (rtype[hcBCC_9[nBCC_9[f]][0]]==1) {
                    nA_cen_BCC_9++;
                    n_distro_cen_BCC_9[1]++;
                    n_distro_shell_BCC_9[number_of_A-1]++;
                }
                else {
                    nB_cen_BCC_9++;
                    n_distro_cen_BCC_9[0]++;
                    n_distro_shell_BCC_9[number_of_A]++;
                }
            }
            
            ++nBCC_9[f];
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
                if (doClusBLDeviation==1) {
                    bl_mom_BCC_9=resize_1D_double(bl_mom_BCC_9,mBCC_9,mBCC_9+incrStatic);
                }
                mBCC_9=mBCC_9+incrStatic;
            }
            for (k=0; k<9; ++k) hcBCC_9[nBCC_9[f]][k]=trial[k];
            
            if (doDynamics==1 && dyn_mBCC_9!=-1) {
                if (doSubClusts==1 && dyn_m6A!=-1) {
                    do_sub=1;
                    sub[0]=-1;
                    sub[1]=-1;
                    sub[2]=-1;
                    sub[3]=-1;
                    sub[4]=dyn_up_sp4b[i];
                    sub[5]=dyn_up_sp4c[mem_sp4c[sp4b[i][j2]][j]];
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hcBCC_9[nBCC_9[f]], f, clusSize, &dyn_nBCC_9, &dyn_mBCC_9, &dyn_lBCC_9, &dyn_hcBCC_9, do_up, dummy_up, nBCC_9[f], do_sub, n_sub, &dyn_sub_BCC_9, sub);
            }
            if(ach[hcBCC_9[nBCC_9[f]][1]] == 'C') ach[hcBCC_9[nBCC_9[f]][1]] = ach_shell[hcBCC_9[nBCC_9[f]][1]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][2]] == 'C') ach[hcBCC_9[nBCC_9[f]][2]] = ach_shell[hcBCC_9[nBCC_9[f]][2]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][3]] == 'C') ach[hcBCC_9[nBCC_9[f]][3]] = ach_shell[hcBCC_9[nBCC_9[f]][3]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][4]] == 'C') ach[hcBCC_9[nBCC_9[f]][4]] = ach_shell[hcBCC_9[nBCC_9[f]][4]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][5]] == 'C') ach[hcBCC_9[nBCC_9[f]][5]] = ach_shell[hcBCC_9[nBCC_9[f]][5]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][6]] == 'C') ach[hcBCC_9[nBCC_9[f]][6]] = ach_shell[hcBCC_9[nBCC_9[f]][6]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][7]] == 'C') ach[hcBCC_9[nBCC_9[f]][7]] = ach_shell[hcBCC_9[nBCC_9[f]][7]] = 'B';
            if(ach[hcBCC_9[nBCC_9[f]][8]] == 'C') ach[hcBCC_9[nBCC_9[f]][8]] = ach_shell[hcBCC_9[nBCC_9[f]][8]] = 'B';
            ach[hcBCC_9[nBCC_9[f]][0]] = ach_cen[hcBCC_9[nBCC_9[f]][0]] = 'F';
            
            if (doBondedCen==1) {
                n_bonded_to_cen_BCC_9+=cnb[hcBCC_9[nBCC_9[f]][0]];
                n_distro_bonded_to_cen_BCC_9[cnb[hcBCC_9[nBCC_9[f]][0]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hcBCC_9[nBCC_9[f]][binAcnt],hcBCC_9[nBCC_9[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hcBCC_9[nBCC_9[f]][binAcnt]][Bonds_cnb_j(hcBCC_9[nBCC_9[f]][binAcnt],hcBCC_9[nBCC_9[f]][binBcnt])],BLDistroBCC_9,&BLDistroNoSamplesBCC_9);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hcBCC_9[nBCC_9[f]][binAcnt]]==1) {
                        nABCC_9++;
                        number_of_A++;
                    }
                    else nBBCC_9++;

                    if (binAcnt!=0) {
                        if (rtype[hcBCC_9[nBCC_9[f]][binAcnt]]==1) nA_shell_BCC_9++;
                        else nB_shell_BCC_9++;
                    }
                }
                n_distro_BCC_9[number_of_A]++;
                
                if (rtype[hcBCC_9[nBCC_9[f]][0]]==1) {
                    nA_cen_BCC_9++;
                    n_distro_cen_BCC_9[1]++;
                    n_distro_shell_BCC_9[number_of_A-1]++;
                }
                else {
                    nB_cen_BCC_9++;
                    n_distro_cen_BCC_9[0]++;
                    n_distro_shell_BCC_9[number_of_A]++;
                }
            }
            
            ++nBCC_9[f];
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

void Clusters_GetBCC_15(int f) {    // Detect 15 particle BCC clusters
    int i,j,k,l,m;
    int no_sp4cs,noSP4s;
    int sj[5];
    char *ach, *ach_cen, *ach_shell;
    char errMsg[1000];
    int binAcnt, binBcnt;
    int number_of_A;
    int clusSize=15;
    int do_up=0;
    int *dummy_up=NULL;
    int do_sub=0;
    int n_sub=6;
    int sub[6];
    
    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Clusters_GetBCC_15(): ach[] malloc out of memory\n");  Error(errMsg); }
    ach_cen=malloc(N*sizeof(char)); if (ach_cen==NULL) { sprintf(errMsg,"Clusters_GetBCC_15(): ach_cen[] malloc out of memory\n");  Error(errMsg); }
    ach_shell=malloc(N*sizeof(char));   if (ach_shell==NULL) { sprintf(errMsg,"Clusters_GetBCC_15(): ach_shell[] malloc out of memory\n");  Error(errMsg); }
    for(i=0; i<N; ++i) ach[i] = ach_cen[i] = ach_shell[i] = 'C';

    for (i=0; i<nsp4c[f]; i++) {
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15[f] == mBCC_15) { 
            hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
            if (doClusBLDeviation==1) {
                bl_mom_BCC_15=resize_1D_double(bl_mom_BCC_15,mBCC_15,mBCC_15+incrStatic);
            }
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
                if (nBCC_15[f] == mBCC_15) { 
                    hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
                    if (doClusBLDeviation==1) {
                        bl_mom_BCC_15=resize_1D_double(bl_mom_BCC_15,mBCC_15,mBCC_15+incrStatic);
                    }
                    mBCC_15=mBCC_15+incrStatic;
                }
                
                quickSort(&hcBCC_15[nBCC_15[f]][1],6);
                quickSort(&hcBCC_15[nBCC_15[f]][7],8);  
                
                if (doDynamics==1 && dyn_mBCC_15!=-1) {
                    if (doSubClusts==1 && dyn_m6A!=-1) {
                        do_sub=1;
                        sub[0]=dyn_up_sp4c[i];
                        sub[1]=dyn_up_sp4c[sj[0]];
                        sub[2]=dyn_up_sp4c[sj[1]];
                        sub[3]=dyn_up_sp4c[sj[2]];
                        sub[4]=dyn_up_sp4c[sj[3]];
                        sub[5]=dyn_up_sp4c[sj[4]];
                        quickSort(&sub[0],n_sub);
                    }
                    else do_sub=0;
                    do_up=0;
                    Dyn_add(hcBCC_15[nBCC_15[f]], f, clusSize, &dyn_nBCC_15, &dyn_mBCC_15, &dyn_lBCC_15, &dyn_hcBCC_15, do_up, dummy_up, nBCC_15[f], do_sub, n_sub, &dyn_sub_BCC_15, sub);
                }
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
                
                if (doBondedCen==1) {
                    n_bonded_to_cen_BCC_15+=cnb[hcBCC_15[nBCC_15[f]][0]];
                    n_distro_bonded_to_cen_BCC_15[cnb[hcBCC_15[nBCC_15[f]][0]]]++;
                }
            
                if (doClusBLDistros==1) {
                    for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                        for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                            if (Bonds_BondCheck(hcBCC_15[nBCC_15[f]][binAcnt],hcBCC_15[nBCC_15[f]][binBcnt])==1) {
                                Bonds_TickBLDistro(bondlengths[hcBCC_15[nBCC_15[f]][binAcnt]][Bonds_cnb_j(hcBCC_15[nBCC_15[f]][binAcnt],hcBCC_15[nBCC_15[f]][binBcnt])],BLDistroBCC_15,&BLDistroNoSamplesBCC_15);
                            }
                        }
                    }
                }
                
                if (doClusComp==1) {
                    number_of_A=0;
                    for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                        if (rtype[hcBCC_15[nBCC_15[f]][binAcnt]]==1) {
                            nABCC_15++;
                            number_of_A++;
                        }
                        else nBBCC_15++;
                    
                        if (binAcnt!=0) {
                            if (rtype[hcBCC_15[nBCC_15[f]][binAcnt]]==1) nA_shell_BCC_15++;
                            else nB_shell_BCC_15++;
                        }
                    }
                    n_distro_BCC_15[number_of_A]++;
                    
                    if (rtype[hcBCC_15[nBCC_15[f]][0]]==1) {
                        nA_cen_BCC_15++;
                        n_distro_cen_BCC_15[1]++;
                        n_distro_shell_BCC_15[number_of_A-1]++;
                    }
                    else {
                        nB_cen_BCC_15++;
                        n_distro_cen_BCC_15[0]++;
                        n_distro_shell_BCC_15[number_of_A]++;
                    }
                }

                ++nBCC_15[f];
            }
        }
        // we may have an BCC_15 cluster, build it into hcBCC_15 then overwrite it later if it aint
        if (nBCC_15[f] == mBCC_15) { 
            hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
            if (doClusBLDeviation==1) {
                bl_mom_BCC_15=resize_1D_double(bl_mom_BCC_15,mBCC_15,mBCC_15+incrStatic);
            }
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
            if (nBCC_15[f] == mBCC_15) { 
                hcBCC_15=resize_2D_int(hcBCC_15,mBCC_15,mBCC_15+incrStatic,clusSize,-1);
                if (doClusBLDeviation==1) {
                    bl_mom_BCC_15=resize_1D_double(bl_mom_BCC_15,mBCC_15,mBCC_15+incrStatic);
                }
                mBCC_15=mBCC_15+incrStatic;
            }
            
            quickSort(&hcBCC_15[nBCC_15[f]][1],6);
            quickSort(&hcBCC_15[nBCC_15[f]][7],8);  
            
            if (doDynamics==1 && dyn_mBCC_15!=-1) {
                if (doSubClusts==1 && dyn_m6A!=-1) {
                    do_sub=1;
                    sub[0]=dyn_up_sp4c[i];
                    sub[1]=dyn_up_sp4c[sj[0]];
                    sub[2]=dyn_up_sp4c[sj[1]];
                    sub[3]=dyn_up_sp4c[sj[2]];
                    sub[4]=dyn_up_sp4c[sj[3]];
                    sub[5]=dyn_up_sp4c[sj[4]];
                    quickSort(&sub[0],n_sub);
                }
                else do_sub=0;
                do_up=0;
                Dyn_add(hcBCC_15[nBCC_15[f]], f, clusSize, &dyn_nBCC_15, &dyn_mBCC_15, &dyn_lBCC_15, &dyn_hcBCC_15, do_up, dummy_up, nBCC_15[f], do_sub, n_sub, &dyn_sub_BCC_15, sub);
            }
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
            
            if (doBondedCen==1) {
                n_bonded_to_cen_BCC_15+=cnb[hcBCC_15[nBCC_15[f]][0]];
                n_distro_bonded_to_cen_BCC_15[cnb[hcBCC_15[nBCC_15[f]][0]]]++;
            }
            
            if (doClusBLDistros==1) {
                for (binAcnt=0; binAcnt<clusSize-1; binAcnt++) {
                    for (binBcnt=binAcnt+1; binBcnt<clusSize; binBcnt++) {
                        if (Bonds_BondCheck(hcBCC_15[nBCC_15[f]][binAcnt],hcBCC_15[nBCC_15[f]][binBcnt])==1) {
                            Bonds_TickBLDistro(bondlengths[hcBCC_15[nBCC_15[f]][binAcnt]][Bonds_cnb_j(hcBCC_15[nBCC_15[f]][binAcnt],hcBCC_15[nBCC_15[f]][binBcnt])],BLDistroBCC_15,&BLDistroNoSamplesBCC_15);
                        }
                    }
                }
            }
            
            if (doClusComp==1) {
                number_of_A=0;
                for (binAcnt=0; binAcnt<clusSize; binAcnt++) {
                    if (rtype[hcBCC_15[nBCC_15[f]][binAcnt]]==1) {
                        nABCC_15++;
                        number_of_A++;
                    }
                    else nBBCC_15++;
                    
                    if (binAcnt!=0) {
                        if (rtype[hcBCC_15[nBCC_15[f]][binAcnt]]==1) nA_shell_BCC_15++;
                        else nB_shell_BCC_15++;
                    }
                }
                n_distro_BCC_15[number_of_A]++;
                
                if (rtype[hcBCC_15[nBCC_15[f]][0]]==1) {
                    nA_cen_BCC_15++;
                    n_distro_cen_BCC_15[1]++;
                    n_distro_shell_BCC_15[number_of_A-1]++;
                }
                else {
                    nB_cen_BCC_15++;
                    n_distro_cen_BCC_15[0]++;
                    n_distro_shell_BCC_15[number_of_A]++;
                }
            }

            ++nBCC_15[f];
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