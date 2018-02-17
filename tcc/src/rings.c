#include "string.h"
#include "rings.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Rings_gSP3(int f, int n0) {	// get SP3/4/5 rings including particle n0
    int i,j;
    int n1, n2;

    for (i=0; i<cnb[n0]-1; i++){
        n1=bNums[n0][i];
        if (n1 < n0) continue;
        for (j=i+1; j<cnb[n0]; ++j){
            n2=bNums[n0][j];
            if (n2<n0) continue;
            if (Bonds_BondCheck(n1,n2)) { // is n1 bonded to n2
                if (n1<n2) Rings_aSP3(f,n0,n1,n2); // SP3 found, check type and store
                else Rings_aSP3(f,n0,n2,n1); // SP3 found, check type and store
            }
            else { // not SP3, search for SP4 & SP5
                if (dosp4==1) {
                    if (n1<n2) Rings_gSP4(f,n0,n1,n2);
                    else Rings_gSP4(f,n0,n2,n1);
                }
            }
        }
    }
}

void Rings_gSP4(int f, int n0, int n1, int n2) {    // {n0,n1,n2} is not an SP3 ring, is it an SP4 or SP5 ring?
    int i;
    int n3;
    
    for (i=0; i<cnb[n1]; ++i) {
        n3=bNums[n1][i];
        if (n3 <= n0) continue;
        if (!Bonds_BondCheck(n0,n3)) {  // n1 not bonded to n2 & n0 not bonded to n3
            if (Bonds_BondCheck(n2,n3)) { // 4 membered ring found 
                Rings_aSP4(f,n0,n1,n3,n2); // check SP4 type and store  
            }
            else{ // n1 not bonded to n3
                if (dosp5==1) Rings_gSP5(f,n0,n1,n3,n2);
            }
        }       
    }
}

void Rings_gSP5(int f, int n0, int n1, int n2, int n3) {    // {n0,n1,n2,n3} is not an SP4 ring, is it an SP5 ring?
    int i,j;
    int n4,n5;
    int bond4_1;
    
    for (i=0; i<cnb[n2]; ++i){
        n4=bNums[n2][i];
        if(n4 < n0 || n4 == n3) continue; // Now: is n4 bonded to n1 and not to n2 or n0
        bond4_1 = 0;
        for (j=0; j<cnb[n4]; ++j){
            n5=bNums[n4][j];
            if (n5==n3) bond4_1 = 1;
            if (n5==n1 || n5==n0) break; // Not SP ring
        }
        if (j==cnb[n4] && bond4_1==1) {
            Rings_aSP5(f,n0, n1, n2, n4, n3); // check SP5 type and store 
        }
    }
}

void Rings_aSP3(int f, int n0, int n1, int n2) {    // Take {n0,n1,n2}, check SP3 ring and if so detect SP3a/b/c cluster
    int i, j;
    int type = 0;
    int cp[2];  // common spindles - particles bonded to all members of three membered ring
    int bcheck;

    cp[0]=cp[1]=-1;
    for (i=0; i<cnb[n0]; ++i) {
        j = bNums[n0][i];
        bcheck = j == n1 || j == n2;
        if (bcheck) continue;
        bcheck = Bonds_BondCheck(n1,j)==1 && Bonds_BondCheck(n2,j)==1;
        if (bcheck) {
            if (type<2) {
                cp[type] = j;
                type++;
            }
            else type++;
        }
    }
    
    if (maxto3<type) maxto3=type;
    
    if (type==0 && dosp3a==1) {
        if (nsp3a[f] == msp3a) {
            sp3a=resize_2D_int(sp3a,msp3a,msp3a+incrStatic,3,-1);
            msp3a=msp3a+incrStatic;
        }
        sp3a[nsp3a[f]][0] = n0;
        sp3a[nsp3a[f]][1] = n1;
        sp3a[nsp3a[f]][2] = n2;
        
        ++nsp3a[f];
    }
    else if (type==1 && dosp3b==1) {
        if (nsp3b[f] == msp3b) { 
            sp3b=resize_2D_int(sp3b,msp3b,msp3b+incrStatic,4,-1);
            msp3b=msp3b+incrStatic;
        }
        sp3b[nsp3b[f]][0] = n0;
        sp3b[nsp3b[f]][1] = n1;
        sp3b[nsp3b[f]][2] = n2;
        sp3b[nsp3b[f]][3] = cp[0];

        add_mem_sp3b(n0, f);
        add_mem_sp3b(n1, f);
        add_mem_sp3b(n2, f);
        add_mem_sp3b(cp[0], f);
        
        ++nsp3b[f];
    }
    else if (type==2 && dosp3c==1) {
        if (nsp3c[f] == msp3c) { 
            sp3c=resize_2D_int(sp3c,msp3c,msp3c+incrStatic,5,-1);
            msp3c=msp3c+incrStatic;
        }
        sp3c[nsp3c[f]][0] = n0;
        sp3c[nsp3c[f]][1] = n1;
        sp3c[nsp3c[f]][2] = n2; 
        if (cp[0]<cp[1]) {
            sp3c[nsp3c[f]][3] = cp[0];
            sp3c[nsp3c[f]][4] = cp[1];
        }
        else {
            sp3c[nsp3c[f]][3] = cp[1];
            sp3c[nsp3c[f]][4] = cp[0];
        }

        add_mem_sp3c(n0, f);
        add_mem_sp3c(n1, f);
        add_mem_sp3c(n2, f);
        add_mem_sp3c(cp[0], f);
        add_mem_sp3c(cp[1], f);

        if (Bonds_BondCheck(sp3c[nsp3c[f]][3],sp3c[nsp3c[f]][4])==1) nsp3c_spindlebonds[f]++;
        
        ++nsp3c[f];
    }
    else if (dosp3a==1) {
        if (nsp3a[f] == msp3a) { 
            sp3a=resize_2D_int(sp3a,msp3a,msp3a+incrStatic,3,-1);
            msp3a=msp3a+incrStatic;
        }
        sp3a[nsp3a[f]][0] = n0;
        sp3a[nsp3a[f]][1] = n1;
        sp3a[nsp3a[f]][2] = n2;

        ++nsp3a[f];
        ++nsp3_excess_spindles[f];
    }
    
    ++nsp3[f];
}

void Rings_aSP4(int f, int n0, int n1, int n2, int n3) {    // Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
    int i, j;
    int type = 0;
    int cp[2];  // common spindles - particles bonded to all members of three membered ring
    int bcheck;

    cp[0]=cp[1]=-1;
    for (i=0; i<cnb[n0]; ++i) {
        j = bNums[n0][i];
        bcheck = j == n1 || j == n3;
        if (bcheck) continue;
        bcheck = Bonds_BondCheck(n1,j)==1 && Bonds_BondCheck(n2,j)==1 && Bonds_BondCheck(n3,j)==1;
        if (bcheck) {
            if (type<2) {
                cp[type] = j;
                type++;
            }
            else type++;
        }
    }
    
    if (maxto4<type) maxto4=type;
    
    if (type==0 && dosp4a==1) {
        if (nsp4a[f] == msp4a) { 
            sp4a=resize_2D_int(sp4a,msp4a,msp4a+incrStatic,4,-1);
            msp4a=msp4a+incrStatic;
        }
        sp4a[nsp4a[f]][0] = n0;
        sp4a[nsp4a[f]][1] = n1;
        sp4a[nsp4a[f]][2] = n2;
        sp4a[nsp4a[f]][3] = n3;

        ++nsp4a[f];
    }
    else if (type==1 && dosp4b==1) {
        if (nsp4b[f] == msp4b) { 
            sp4b=resize_2D_int(sp4b,msp4b,msp4b+incrStatic,5,-1);
            msp4b=msp4b+incrStatic;
        }
        sp4b[nsp4b[f]][0] = n0;
        sp4b[nsp4b[f]][1] = n1;
        sp4b[nsp4b[f]][2] = n2;
        sp4b[nsp4b[f]][3] = n3;
        sp4b[nsp4b[f]][4] = cp[0];

        add_mem_sp4b(n0, f);
        add_mem_sp4b(n1, f);
        add_mem_sp4b(n2, f);
        add_mem_sp4b(n3, f);
        add_mem_sp4b(cp[0], f);

        ++nsp4b[f];
    }
    else if (type==2 && dosp4c==1) {
        if (nsp4c[f] == msp4c) { 
            sp4c=resize_2D_int(sp4c,msp4c,msp4c+incrStatic,6,-1);
            msp4c=msp4c+incrStatic;
        }
        sp4c[nsp4c[f]][0] = n0;
        sp4c[nsp4c[f]][1] = n1;
        sp4c[nsp4c[f]][2] = n2;
        sp4c[nsp4c[f]][3] = n3; 
        if (cp[0]<cp[1]) {
            sp4c[nsp4c[f]][4] = cp[0];
            sp4c[nsp4c[f]][5] = cp[1];
        }
        else {
            sp4c[nsp4c[f]][4] = cp[1];
            sp4c[nsp4c[f]][5] = cp[0];
        }

        add_mem_sp4c(n0, f);
        add_mem_sp4c(n1, f);
        add_mem_sp4c(n2, f);
        add_mem_sp4c(n3, f);
        add_mem_sp4c(cp[0], f);
        add_mem_sp4c(cp[1], f);

        if (Bonds_BondCheck(sp4c[nsp4c[f]][4],sp4c[nsp4c[f]][5])==1) nsp4c_spindlebonds[f]++;
            
        // hc6A key: (SP4_1, SP4_2, SP4_3, SP4_4, s1, s2)
        
        ++nsp4c[f];
    }
    else if (dosp4a==1) {
        if (nsp4a[f] == msp4a) { 
            sp4a=resize_2D_int(sp4a,msp4a,msp4a+incrStatic,4,-1);
            msp4a=msp4a+incrStatic;
        }
        sp4a[nsp4a[f]][0] = n0;
        sp4a[nsp4a[f]][1] = n1;
        sp4a[nsp4a[f]][2] = n2;
        sp4a[nsp4a[f]][3] = n3;
        
        ++nsp4a[f];
        ++nsp4_excess_spindles[f];
    }
    
    ++nsp4[f];
}

void Rings_aSP5(int f, int n0, int n1, int n2, int n3, int n4) {    // Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster
    int i, j;
    int type = 0;
    int cp[2];  // common spindles - particles bonded to all members of three membered ring
    int bcheck;

    cp[0]=cp[1]=-1;
    for (i=0; i<cnb[n0]; ++i) {
        j = bNums[n0][i];
        bcheck = j == n1 || j == n4;
        if (bcheck) continue;
        bcheck = Bonds_BondCheck(n1,j)==1 && Bonds_BondCheck(n2,j)==1 && Bonds_BondCheck(n3,j)==1 && Bonds_BondCheck(n4,j)==1;
        if (bcheck) {
            if (type<2) {
                cp[type] = j;
                type++;
            }
            else type++;
        }
    }
    
    if (maxto5<type) maxto5=type;
    
    if (type==0 && dosp5a==1) { // Now store ring
        if (nsp5a[f] == msp5a) { 
            sp5a=resize_2D_int(sp5a,msp5a,msp5a+incrStatic,5,-1);
            msp5a=msp5a+incrStatic;
        }
        sp5a[nsp5a[f]][0] = n0;
        sp5a[nsp5a[f]][1] = n1;
        sp5a[nsp5a[f]][2] = n2;
        sp5a[nsp5a[f]][3] = n3;
        sp5a[nsp5a[f]][4] = n4;
        
        ++nsp5a[f];
    }
    else if (type==1 && dosp5b==1) {
        if (nsp5b[f] == msp5b) { 
            sp5b=resize_2D_int(sp5b,msp5b,msp5b+incrStatic,6,-1);
            msp5b=msp5b+incrStatic;
        }
        sp5b[nsp5b[f]][0] = n0;
        sp5b[nsp5b[f]][1] = n1;
        sp5b[nsp5b[f]][2] = n2;
        sp5b[nsp5b[f]][3] = n3;
        sp5b[nsp5b[f]][4] = n4;
        sp5b[nsp5b[f]][5] = cp[0];

        add_mem_sp5b(n0, f);
        add_mem_sp5b(n1, f);
        add_mem_sp5b(n2, f);
        add_mem_sp5b(n3, f);
        add_mem_sp5b(n4, f);
        add_mem_sp5b(cp[0], f);
        
        ++nsp5b[f];
    }
    else if (type==2 && dosp5c==1) {
        if (nsp5c[f] == msp5c) { 
            sp5c=resize_2D_int(sp5c,msp5c,msp5c+incrStatic,7,-1);
            msp5c=msp5c+incrStatic;
        }
        sp5c[nsp5c[f]][0] = n0;
        sp5c[nsp5c[f]][1] = n1;
        sp5c[nsp5c[f]][2] = n2;
        sp5c[nsp5c[f]][3] = n3; 
        sp5c[nsp5c[f]][4] = n4; 
        if (cp[0]<cp[1]) {
            sp5c[nsp5c[f]][5] = cp[0];
            sp5c[nsp5c[f]][6] = cp[1];
        }
        else {
            sp5c[nsp5c[f]][5] = cp[1];
            sp5c[nsp5c[f]][6] = cp[0];
        }

        add_mem_sp5c(n0, f);
        add_mem_sp5c(n1, f);
        add_mem_sp5c(n2, f);
        add_mem_sp5c(n3, f);
        add_mem_sp5c(n4, f);
        add_mem_sp5c(cp[0], f);
        add_mem_sp5c(cp[1], f);

        if (Bonds_BondCheck(sp5c[nsp5c[f]][5],sp5c[nsp5c[f]][6])==1) nsp5c_spindlebonds[f]++;
        
        ++nsp5c[f];
    }
    else if (dosp5a==1) {   // Now store ring
        if (nsp5a[f] == msp5a) { 
            sp5a=resize_2D_int(sp5a,msp5a,msp5a+incrStatic,5,-1);
            msp5a=msp5a+incrStatic;
        }
        sp5a[nsp5a[f]][0] = n0;
        sp5a[nsp5a[f]][1] = n1;
        sp5a[nsp5a[f]][2] = n2;
        sp5a[nsp5a[f]][3] = n3;
        sp5a[nsp5a[f]][4] = n4;
        
        ++nsp5a[f];
        ++nsp5_excess_spindles[f];
    }
    
    ++nsp5[f];
}

void Rings_setSP3c(int f) { // store cluster 5A D3h from Bonds_aSP3
    int i;
    char *ach, errMsg[1000];

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Rings_setSP3c(): ach[] malloc out of memory\n");   Error(errMsg); }

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp3a[f]; i++) {
        ach[sp3a[i][0]] = 'B';
        ach[sp3a[i][1]] = 'B';
        ach[sp3a[i][2]] = 'B';
    }
    for (i=0; i<N; ++i) ssp3a[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp3b[f]; i++) {
        if (ach[sp3b[i][0]] == 'C') ach[sp3b[i][0]] = 'B';
        if (ach[sp3b[i][1]] == 'C') ach[sp3b[i][1]] = 'B';
        if (ach[sp3b[i][2]] == 'C') ach[sp3b[i][2]] = 'B';
        ach[sp3b[i][3]] = 'O';
    }
    for (i=0; i<N; ++i) ssp3b[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp3c[f]; i++) {
        if (ach[sp3c[i][0]] == 'C') ach[sp3c[i][0]] = 'B';
        if (ach[sp3c[i][1]] == 'C') ach[sp3c[i][1]] = 'B';
        if (ach[sp3c[i][2]] == 'C') ach[sp3c[i][2]] = 'B';
        ach[sp3c[i][3]] = 'O';
        ach[sp3c[i][4]] = 'O';
    }
    for (i=0; i<N; ++i) ssp3c[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp3a[f]; i++) {
        ach[sp3a[i][0]] = 'B';
        ach[sp3a[i][1]] = 'B';
        ach[sp3a[i][2]] = 'B';
    }
    for (i=0; i<nsp3b[f]; i++) {
        ach[sp3b[i][0]] = 'B';
        ach[sp3b[i][1]] = 'B';
        ach[sp3b[i][2]] = 'B';
    }
    for (i=0; i<nsp3c[f]; i++) {
        ach[sp3c[i][0]] = 'B';
        ach[sp3c[i][1]] = 'B';
        ach[sp3c[i][2]] = 'B';
    }
    for (i=0; i<N; ++i) ssp3[i]=ach[i];

    free(ach);
}

void Rings_setSP4c(int f) { // store cluster 6A Oh from Bonds_aSP4()
    int i;
    char *ach, errMsg[1000];

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Rings_setSP4c(): ach[] malloc out of memory\n");   Error(errMsg); }

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp4a[f]; i++) {
        ach[sp4a[i][0]] = 'B';
        ach[sp4a[i][1]] = 'B';
        ach[sp4a[i][2]] = 'B';
        ach[sp4a[i][3]] = 'B';
    }
    for (i=0; i<N; ++i) ssp4a[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp4b[f]; i++) {
        if (ach[sp4b[i][0]] == 'C') ach[sp4b[i][0]] = 'B';
        if (ach[sp4b[i][1]] == 'C') ach[sp4b[i][1]] = 'B';
        if (ach[sp4b[i][2]] == 'C') ach[sp4b[i][2]] = 'B';
        if (ach[sp4b[i][3]] == 'C') ach[sp4b[i][3]] = 'B';
        ach[sp4b[i][4]] = 'O';
    }
    for (i=0; i<N; ++i) ssp4b[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp4c[f]; ++i) {
        if (ach[sp4c[i][0]] == 'C') ach[sp4c[i][0]] = 'B';
        if (ach[sp4c[i][1]] == 'C') ach[sp4c[i][1]] = 'B';
        if (ach[sp4c[i][2]] == 'C') ach[sp4c[i][2]] = 'B';
        if (ach[sp4c[i][3]] == 'C') ach[sp4c[i][3]] = 'B';
        ach[sp4c[i][4]] = 'O';
        ach[sp4c[i][5]] = 'O';
    }
    for (i=0; i<N; ++i) ssp4c[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp4a[f]; i++) {
        ach[sp4a[i][0]] = 'B';
        ach[sp4a[i][1]] = 'B';
        ach[sp4a[i][2]] = 'B';
        ach[sp4a[i][3]] = 'B';
    }
    for (i=0; i<nsp4b[f]; i++) {
        ach[sp4b[i][0]] = 'B';
        ach[sp4b[i][1]] = 'B';
        ach[sp4b[i][2]] = 'B';
        ach[sp4b[i][3]] = 'B';
    }
    for (i=0; i<nsp4c[f]; ++i) {
        ach[sp4c[i][0]] = 'B';
        ach[sp4c[i][1]] = 'B';
        ach[sp4c[i][2]] = 'B';
        ach[sp4c[i][3]] = 'B';
    }
    for (i=0; i<N; ++i) ssp4[i]=ach[i];

    
    free(ach);
}

void Rings_setSP5c(int f) { // store cluster 7A D5h from Bonds_aSP5()
    int i;
    char *ach, errMsg[1000];

    ach=malloc(N*sizeof(char)); if (ach==NULL) { sprintf(errMsg,"Rings_setSP5c(): ach[] malloc out of memory\n");   Error(errMsg); }

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp5a[f]; i++) {
        ach[sp5a[i][0]] = 'B';
        ach[sp5a[i][1]] = 'B';
        ach[sp5a[i][2]] = 'B';
        ach[sp5a[i][3]] = 'B';
        ach[sp5a[i][4]] = 'B';
    }
    for (i=0; i<N; ++i) ssp5a[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp5b[f]; i++) {
        if (ach[sp5b[i][0]] == 'C') ach[sp5b[i][0]] = 'B';
        if (ach[sp5b[i][1]] == 'C') ach[sp5b[i][1]] = 'B';
        if (ach[sp5b[i][2]] == 'C') ach[sp5b[i][2]] = 'B';
        if (ach[sp5b[i][3]] == 'C') ach[sp5b[i][3]] = 'B';
        if (ach[sp5b[i][4]] == 'C') ach[sp5b[i][4]] = 'B';
        ach[sp5b[i][5]] = 'O';
    }
    for (i=0; i<N; ++i) ssp5b[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp5c[f]; ++i) {
        if (ach[sp5c[i][0]] == 'C') ach[sp5c[i][0]] = 'B';
        if (ach[sp5c[i][1]] == 'C') ach[sp5c[i][1]] = 'B';
        if (ach[sp5c[i][2]] == 'C') ach[sp5c[i][2]] = 'B';
        if (ach[sp5c[i][3]] == 'C') ach[sp5c[i][3]] = 'B';
        if (ach[sp5c[i][4]] == 'C') ach[sp5c[i][4]] = 'B';
        ach[sp5c[i][5]] = 'O';
        ach[sp5c[i][6]] = 'O';
    }
    for (i=0; i<N; ++i) ssp5c[i]=ach[i];

    memset(ach, 'C', N*sizeof(*ach));
    for (i=0; i<nsp5a[f]; i++) {
        ach[sp5a[i][0]] = 'B';
        ach[sp5a[i][1]] = 'B';
        ach[sp5a[i][2]] = 'B';
        ach[sp5a[i][3]] = 'B';
        ach[sp5a[i][4]] = 'B';
    }
    for (i=0; i<nsp5b[f]; i++) {
        ach[sp5b[i][0]] = 'B';
        ach[sp5b[i][1]] = 'B';
        ach[sp5b[i][2]] = 'B';
        ach[sp5b[i][3]] = 'B';
        ach[sp5b[i][4]] = 'B';
    }
    for (i=0; i<nsp5c[f]; ++i) {
        ach[sp5c[i][0]] = 'B';
        ach[sp5c[i][1]] = 'B';
        ach[sp5c[i][2]] = 'B';
        ach[sp5c[i][3]] = 'B';
        ach[sp5c[i][4]] = 'B';
    }
    for (i=0; i<N; ++i) ssp5[i]=ach[i];
    
    free(ach);
}

void add_mem_sp3b(int particle_ID, int frame) {
    int binAcnt;

    mem_sp3b[particle_ID][nmem_sp3b[particle_ID]]=nsp3b[frame];
    nmem_sp3b[particle_ID]++;
    if (nmem_sp3b[particle_ID] >= mmem_sp3b) {
        for (binAcnt=0; binAcnt<N; binAcnt++) {
            mem_sp3b[binAcnt]=resize_1D_int(mem_sp3b[binAcnt],mmem_sp3b,mmem_sp3b+incrClustPerPart);
        }
        mmem_sp3b=mmem_sp3b+incrClustPerPart;
    }
}

void add_mem_sp3c(int particle_ID, int frame) {
    int binAcnt;

    mem_sp3c[particle_ID][nmem_sp3c[particle_ID]]=nsp3c[frame];
    nmem_sp3c[particle_ID]++;
    if (nmem_sp3c[particle_ID] >= mmem_sp3c) {
        for (binAcnt=0; binAcnt<N; binAcnt++) {
            mem_sp3c[binAcnt]=resize_1D_int(mem_sp3c[binAcnt],mmem_sp3c,mmem_sp3c+incrClustPerPart);
        }
        mmem_sp3c=mmem_sp3c+incrClustPerPart;
    }
}

void add_mem_sp4b(int particle_ID, int frame) {
    int binAcnt;

    mem_sp4b[particle_ID][nmem_sp4b[particle_ID]] = nsp4b[frame];
    nmem_sp4b[particle_ID]++;
    if (nmem_sp4b[particle_ID] >= mmem_sp4b) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp4b[binAcnt] = resize_1D_int(mem_sp4b[binAcnt], mmem_sp4b, mmem_sp4b + incrClustPerPart);
        }
        mmem_sp4b = mmem_sp4b + incrClustPerPart;
    }
}

void add_mem_sp4c(int particle_ID, int frame) {
    int binAcnt;

    mem_sp4c[particle_ID][nmem_sp4c[particle_ID]] = nsp4c[frame];
    nmem_sp4c[particle_ID]++;
    if (nmem_sp4c[particle_ID] >= mmem_sp4c) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp4c[binAcnt] = resize_1D_int(mem_sp4c[binAcnt], mmem_sp4c, mmem_sp4c + incrClustPerPart);
        }
        mmem_sp4c = mmem_sp4c + incrClustPerPart;
    }
}

void add_mem_sp5b(int particle_ID, int frame) {
    int binAcnt;

    mem_sp5b[particle_ID][nmem_sp5b[particle_ID]] = nsp5b[frame];
    nmem_sp5b[particle_ID]++;
    if (nmem_sp5b[particle_ID] >= mmem_sp5b) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp5b[binAcnt] = resize_1D_int(mem_sp5b[binAcnt], mmem_sp5b, mmem_sp5b + incrClustPerPart);
        }
        mmem_sp5b = mmem_sp5b + incrClustPerPart;
    }
}

void add_mem_sp5c(int particle_ID, int frame) {
    int binAcnt;

    mem_sp5c[particle_ID][nmem_sp5c[particle_ID]] = nsp5c[frame];
    nmem_sp5c[particle_ID]++;
    if (nmem_sp5c[particle_ID] >= mmem_sp5c) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp5c[binAcnt] = resize_1D_int(mem_sp5c[binAcnt], mmem_sp5c, mmem_sp5c + incrClustPerPart);
        }
        mmem_sp5c = mmem_sp5c + incrClustPerPart;
    }
}