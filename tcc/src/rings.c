#include "rings.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Rings_gSP3(int n0) {	// get SP3/4/5 rings including particle n0
    int i,j;
    int n1, n2;

    for (i=0; i<cnb[n0]-1; i++){
        n1=bNums[n0][i];
        if (n1 < n0) continue;
        for (j=i+1; j<cnb[n0]; ++j){
            n2=bNums[n0][j];
            if (n2<n0) continue;
            if (Bonds_BondCheck(n1,n2)) { // is n1 bonded to n2
                if (n1<n2) Rings_aSP3(n0, n1, n2); // SP3 found, check type and store
                else Rings_aSP3(n0, n2, n1); // SP3 found, check type and store
            }
            else { // not SP3, search for SP4 & SP5
                if (dosp4==1) {
                    if (n1<n2) Rings_gSP4(n0, n1, n2);
                    else Rings_gSP4(n0, n2, n1);
                }
            }
        }
    }
}

void Rings_gSP4(int n0, int n1, int n2) {    // {n0,n1,n2} is not an SP3 ring, is it an SP4 or SP5 ring?
    int i;
    int n3;
    
    for (i=0; i<cnb[n1]; ++i) {
        n3=bNums[n1][i];
        if (n3 <= n0) continue;
        if (!Bonds_BondCheck(n0,n3)) {  // n1 not bonded to n2 & n0 not bonded to n3
            if (Bonds_BondCheck(n2,n3)) { // 4 membered ring found 
                Rings_aSP4(n0, n1, n3, n2); // check SP4 type and store
            }
            else{ // n1 not bonded to n3
                if (dosp5==1) Rings_gSP5(n0, n1, n3, n2);
            }
        }       
    }
}

void Rings_gSP5(int n0, int n1, int n2, int n3) {    // {n0,n1,n2,n3} is not an SP4 ring, is it an SP5 ring?
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
            Rings_aSP5(n0, n1, n2, n4, n3); // check SP5 type and store
        }
    }
}

void Rings_aSP3(int n0, int n1, int n2) {    // Take {n0,n1,n2}, check SP3 ring and if so detect SP3a/b/c cluster
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
    
    if (type==0 && dosp3a==1) {
        if (nsp3a == msp3a) {
            hcsp3a=resize_2D_int(hcsp3a,msp3a,msp3a+incrStatic,3,-1);
            msp3a=msp3a+incrStatic;
        }
        hcsp3a[nsp3a][0] = n0;
        hcsp3a[nsp3a][1] = n1;
        hcsp3a[nsp3a][2] = n2;
        
        ++nsp3a;
    }
    else if (type==1 && dosp3b==1) {
        if (nsp3b == msp3b) { 
            hcsp3b=resize_2D_int(hcsp3b,msp3b,msp3b+incrStatic,4,-1);
            msp3b=msp3b+incrStatic;
        }
        hcsp3b[nsp3b][0] = n0;
        hcsp3b[nsp3b][1] = n1;
        hcsp3b[nsp3b][2] = n2;
        hcsp3b[nsp3b][3] = cp[0];

        add_mem_sp3b(n0);
        add_mem_sp3b(n1);
        add_mem_sp3b(n2);
        add_mem_sp3b(cp[0]);
        
        ++nsp3b;
    }
    else if (type==2 && dosp3c==1) {
        if (nsp3c == msp3c) { 
            hcsp3c=resize_2D_int(hcsp3c,msp3c,msp3c+incrStatic,5,-1);
            msp3c=msp3c+incrStatic;
        }
        hcsp3c[nsp3c][0] = n0;
        hcsp3c[nsp3c][1] = n1;
        hcsp3c[nsp3c][2] = n2;
        if (cp[0]<cp[1]) {
            hcsp3c[nsp3c][3] = cp[0];
            hcsp3c[nsp3c][4] = cp[1];
        }
        else {
            hcsp3c[nsp3c][3] = cp[1];
            hcsp3c[nsp3c][4] = cp[0];
        }

        add_mem_sp3c(n0);
        add_mem_sp3c(n1);
        add_mem_sp3c(n2);
        add_mem_sp3c(cp[0]);
        add_mem_sp3c(cp[1]);

        ++nsp3c;
    }
    else if (dosp3a==1) {
        if (nsp3a == msp3a) { 
            hcsp3a=resize_2D_int(hcsp3a,msp3a,msp3a+incrStatic,3,-1);
            msp3a=msp3a+incrStatic;
        }
        hcsp3a[nsp3a][0] = n0;
        hcsp3a[nsp3a][1] = n1;
        hcsp3a[nsp3a][2] = n2;

        ++nsp3a;
    }
}

void Rings_aSP4(int n0, int n1, int n2, int n3) {    // Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
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
    
    if (type==0 && dosp4a==1) {
        if (nsp4a == msp4a) { 
            hcsp4a=resize_2D_int(hcsp4a,msp4a,msp4a+incrStatic,4,-1);
            msp4a=msp4a+incrStatic;
        }
        hcsp4a[nsp4a][0] = n0;
        hcsp4a[nsp4a][1] = n1;
        hcsp4a[nsp4a][2] = n2;
        hcsp4a[nsp4a][3] = n3;

        ++nsp4a;
    }
    else if (type==1 && dosp4b==1) {
        if (nsp4b == msp4b) { 
            hcsp4b=resize_2D_int(hcsp4b,msp4b,msp4b+incrStatic,5,-1);
            msp4b=msp4b+incrStatic;
        }
        hcsp4b[nsp4b][0] = n0;
        hcsp4b[nsp4b][1] = n1;
        hcsp4b[nsp4b][2] = n2;
        hcsp4b[nsp4b][3] = n3;
        hcsp4b[nsp4b][4] = cp[0];

        add_mem_sp4b(n0);
        add_mem_sp4b(n1);
        add_mem_sp4b(n2);
        add_mem_sp4b(n3);
        add_mem_sp4b(cp[0]);

        ++nsp4b;
    }
    else if (type==2 && dosp4c==1) {
        if (nsp4c == msp4c) { 
            hcsp4c=resize_2D_int(hcsp4c,msp4c,msp4c+incrStatic,6,-1);
            msp4c=msp4c+incrStatic;
        }
        hcsp4c[nsp4c][0] = n0;
        hcsp4c[nsp4c][1] = n1;
        hcsp4c[nsp4c][2] = n2;
        hcsp4c[nsp4c][3] = n3;
        if (cp[0]<cp[1]) {
            hcsp4c[nsp4c][4] = cp[0];
            hcsp4c[nsp4c][5] = cp[1];
        }
        else {
            hcsp4c[nsp4c][4] = cp[1];
            hcsp4c[nsp4c][5] = cp[0];
        }

        add_mem_sp4c(n0);
        add_mem_sp4c(n1);
        add_mem_sp4c(n2);
        add_mem_sp4c(n3);
        add_mem_sp4c(cp[0]);
        add_mem_sp4c(cp[1]);

        // hc6A key: (SP4_1, SP4_2, SP4_3, SP4_4, s1, s2)
        
        ++nsp4c;
    }
    else if (dosp4a==1) {
        if (nsp4a == msp4a) { 
            hcsp4a=resize_2D_int(hcsp4a,msp4a,msp4a+incrStatic,4,-1);
            msp4a=msp4a+incrStatic;
        }
        hcsp4a[nsp4a][0] = n0;
        hcsp4a[nsp4a][1] = n1;
        hcsp4a[nsp4a][2] = n2;
        hcsp4a[nsp4a][3] = n3;
        
        ++nsp4a;
    }
}

void Rings_aSP5(int n0, int n1, int n2, int n3, int n4) {    // Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster
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

    if (type==0 && dosp5a==1) { // Now store ring
        if (nsp5a == msp5a) { 
            hcsp5a=resize_2D_int(hcsp5a,msp5a,msp5a+incrStatic,5,-1);
            msp5a=msp5a+incrStatic;
        }
        hcsp5a[nsp5a][0] = n0;
        hcsp5a[nsp5a][1] = n1;
        hcsp5a[nsp5a][2] = n2;
        hcsp5a[nsp5a][3] = n3;
        hcsp5a[nsp5a][4] = n4;
        
        ++nsp5a;
    }
    else if (type==1 && dosp5b==1) {
        if (nsp5b == msp5b) { 
            hcsp5b=resize_2D_int(hcsp5b,msp5b,msp5b+incrStatic,6,-1);
            msp5b=msp5b+incrStatic;
        }
        hcsp5b[nsp5b][0] = n0;
        hcsp5b[nsp5b][1] = n1;
        hcsp5b[nsp5b][2] = n2;
        hcsp5b[nsp5b][3] = n3;
        hcsp5b[nsp5b][4] = n4;
        hcsp5b[nsp5b][5] = cp[0];

        add_mem_sp5b(n0);
        add_mem_sp5b(n1);
        add_mem_sp5b(n2);
        add_mem_sp5b(n3);
        add_mem_sp5b(n4);
        add_mem_sp5b(cp[0]);
        
        ++nsp5b;
    }
    else if (type==2 && dosp5c==1) {
        if (nsp5c == msp5c) { 
            hcsp5c=resize_2D_int(hcsp5c,msp5c,msp5c+incrStatic,7,-1);
            msp5c=msp5c+incrStatic;
        }
        hcsp5c[nsp5c][0] = n0;
        hcsp5c[nsp5c][1] = n1;
        hcsp5c[nsp5c][2] = n2;
        hcsp5c[nsp5c][3] = n3;
        hcsp5c[nsp5c][4] = n4;
        if (cp[0]<cp[1]) {
            hcsp5c[nsp5c][5] = cp[0];
            hcsp5c[nsp5c][6] = cp[1];
        }
        else {
            hcsp5c[nsp5c][5] = cp[1];
            hcsp5c[nsp5c][6] = cp[0];
        }

        add_mem_sp5c(n0);
        add_mem_sp5c(n1);
        add_mem_sp5c(n2);
        add_mem_sp5c(n3);
        add_mem_sp5c(n4);
        add_mem_sp5c(cp[0]);
        add_mem_sp5c(cp[1]);

        ++nsp5c;
    }
    else if (dosp5a==1) {   // Now store ring
        if (nsp5a == msp5a) { 
            hcsp5a=resize_2D_int(hcsp5a,msp5a,msp5a+incrStatic,5,-1);
            msp5a=msp5a+incrStatic;
        }
        hcsp5a[nsp5a][0] = n0;
        hcsp5a[nsp5a][1] = n1;
        hcsp5a[nsp5a][2] = n2;
        hcsp5a[nsp5a][3] = n3;
        hcsp5a[nsp5a][4] = n4;
        
        ++nsp5a;
    }
}

void Rings_setSP3c() { // store cluster 5A D3h from Bonds_aSP3
    int i;

    for (i=0; i<nsp3a; i++) {
        ssp3a[hcsp3a[i][0]] = 'B';
        ssp3a[hcsp3a[i][1]] = 'B';
        ssp3a[hcsp3a[i][2]] = 'B';
    }

    for (i=0; i<nsp3b; i++) {
        if (ssp3b[hcsp3b[i][0]] == 'C') ssp3b[hcsp3b[i][0]] = 'B';
        if (ssp3b[hcsp3b[i][1]] == 'C') ssp3b[hcsp3b[i][1]] = 'B';
        if (ssp3b[hcsp3b[i][2]] == 'C') ssp3b[hcsp3b[i][2]] = 'B';
        ssp3b[hcsp3b[i][3]] = 'O';
    }

    for (i=0; i<nsp3c; i++) {
        if (ssp3c[hcsp3c[i][0]] == 'C') ssp3c[hcsp3c[i][0]] = 'B';
        if (ssp3c[hcsp3c[i][1]] == 'C') ssp3c[hcsp3c[i][1]] = 'B';
        if (ssp3c[hcsp3c[i][2]] == 'C') ssp3c[hcsp3c[i][2]] = 'B';
        ssp3c[hcsp3c[i][3]] = 'O';
        ssp3c[hcsp3c[i][4]] = 'O';
    }
}

void Rings_setSP4c() { // store cluster 6A Oh from Bonds_aSP4()
    int i;

    for (i=0; i<nsp4a; i++) {
        ssp4a[hcsp4a[i][0]] = 'B';
        ssp4a[hcsp4a[i][1]] = 'B';
        ssp4a[hcsp4a[i][2]] = 'B';
        ssp4a[hcsp4a[i][3]] = 'B';
    }

    for (i=0; i<nsp4b; i++) {
        if (ssp4b[hcsp4b[i][0]] == 'C') ssp4b[hcsp4b[i][0]] = 'B';
        if (ssp4b[hcsp4b[i][1]] == 'C') ssp4b[hcsp4b[i][1]] = 'B';
        if (ssp4b[hcsp4b[i][2]] == 'C') ssp4b[hcsp4b[i][2]] = 'B';
        if (ssp4b[hcsp4b[i][3]] == 'C') ssp4b[hcsp4b[i][3]] = 'B';
        ssp4b[hcsp4b[i][4]] = 'O';
    }

    for (i=0; i<nsp4c; ++i) {
        if (ssp4c[hcsp4c[i][0]] == 'C') ssp4c[hcsp4c[i][0]] = 'B';
        if (ssp4c[hcsp4c[i][1]] == 'C') ssp4c[hcsp4c[i][1]] = 'B';
        if (ssp4c[hcsp4c[i][2]] == 'C') ssp4c[hcsp4c[i][2]] = 'B';
        if (ssp4c[hcsp4c[i][3]] == 'C') ssp4c[hcsp4c[i][3]] = 'B';
        ssp4c[hcsp4c[i][4]] = 'O';
        ssp4c[hcsp4c[i][5]] = 'O';
    }
}

void Rings_setSP5c() { // store cluster 7A D5h from Bonds_aSP5()
    int i;

    for (i=0; i<nsp5a; i++) {
        ssp5a[hcsp5a[i][0]] = 'B';
        ssp5a[hcsp5a[i][1]] = 'B';
        ssp5a[hcsp5a[i][2]] = 'B';
        ssp5a[hcsp5a[i][3]] = 'B';
        ssp5a[hcsp5a[i][4]] = 'B';
    }

    for (i=0; i<nsp5b; i++) {
        if (ssp5b[hcsp5b[i][0]] == 'C') ssp5b[hcsp5b[i][0]] = 'B';
        if (ssp5b[hcsp5b[i][1]] == 'C') ssp5b[hcsp5b[i][1]] = 'B';
        if (ssp5b[hcsp5b[i][2]] == 'C') ssp5b[hcsp5b[i][2]] = 'B';
        if (ssp5b[hcsp5b[i][3]] == 'C') ssp5b[hcsp5b[i][3]] = 'B';
        if (ssp5b[hcsp5b[i][4]] == 'C') ssp5b[hcsp5b[i][4]] = 'B';
        ssp5b[hcsp5b[i][5]] = 'O';
    }

    for (i=0; i<nsp5c; ++i) {
        if (ssp5c[hcsp5c[i][0]] == 'C') ssp5c[hcsp5c[i][0]] = 'B';
        if (ssp5c[hcsp5c[i][1]] == 'C') ssp5c[hcsp5c[i][1]] = 'B';
        if (ssp5c[hcsp5c[i][2]] == 'C') ssp5c[hcsp5c[i][2]] = 'B';
        if (ssp5c[hcsp5c[i][3]] == 'C') ssp5c[hcsp5c[i][3]] = 'B';
        if (ssp5c[hcsp5c[i][4]] == 'C') ssp5c[hcsp5c[i][4]] = 'B';
        ssp5c[hcsp5c[i][5]] = 'O';
        ssp5c[hcsp5c[i][6]] = 'O';
    }
}

void add_mem_sp3b(int particle_ID) {
    int binAcnt;

    mem_sp3b[particle_ID][nmem_sp3b[particle_ID]]=nsp3b;
    nmem_sp3b[particle_ID]++;
    if (nmem_sp3b[particle_ID] >= mmem_sp3b) {
        for (binAcnt=0; binAcnt<N; binAcnt++) {
            mem_sp3b[binAcnt]=resize_1D_int(mem_sp3b[binAcnt],mmem_sp3b,mmem_sp3b+incrClustPerPart);
        }
        mmem_sp3b=mmem_sp3b+incrClustPerPart;
    }
}

void add_mem_sp3c(int particle_ID) {
    int binAcnt;

    mem_sp3c[particle_ID][nmem_sp3c[particle_ID]]=nsp3c;
    nmem_sp3c[particle_ID]++;
    if (nmem_sp3c[particle_ID] >= mmem_sp3c) {
        for (binAcnt=0; binAcnt<N; binAcnt++) {
            mem_sp3c[binAcnt]=resize_1D_int(mem_sp3c[binAcnt],mmem_sp3c,mmem_sp3c+incrClustPerPart);
        }
        mmem_sp3c=mmem_sp3c+incrClustPerPart;
    }
}

void add_mem_sp4b(int particle_ID) {
    int binAcnt;

    mem_sp4b[particle_ID][nmem_sp4b[particle_ID]] = nsp4b;
    nmem_sp4b[particle_ID]++;
    if (nmem_sp4b[particle_ID] >= mmem_sp4b) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp4b[binAcnt] = resize_1D_int(mem_sp4b[binAcnt], mmem_sp4b, mmem_sp4b + incrClustPerPart);
        }
        mmem_sp4b = mmem_sp4b + incrClustPerPart;
    }
}

void add_mem_sp4c(int particle_ID) {
    int binAcnt;

    mem_sp4c[particle_ID][nmem_sp4c[particle_ID]] = nsp4c;
    nmem_sp4c[particle_ID]++;
    if (nmem_sp4c[particle_ID] >= mmem_sp4c) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp4c[binAcnt] = resize_1D_int(mem_sp4c[binAcnt], mmem_sp4c, mmem_sp4c + incrClustPerPart);
        }
        mmem_sp4c = mmem_sp4c + incrClustPerPart;
    }
}

void add_mem_sp5b(int particle_ID) {
    int binAcnt;

    mem_sp5b[particle_ID][nmem_sp5b[particle_ID]] = nsp5b;
    nmem_sp5b[particle_ID]++;
    if (nmem_sp5b[particle_ID] >= mmem_sp5b) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp5b[binAcnt] = resize_1D_int(mem_sp5b[binAcnt], mmem_sp5b, mmem_sp5b + incrClustPerPart);
        }
        mmem_sp5b = mmem_sp5b + incrClustPerPart;
    }
}

void add_mem_sp5c(int particle_ID) {
    int binAcnt;

    mem_sp5c[particle_ID][nmem_sp5c[particle_ID]] = nsp5c;
    nmem_sp5c[particle_ID]++;
    if (nmem_sp5c[particle_ID] >= mmem_sp5c) {
        for (binAcnt = 0; binAcnt < N; binAcnt++) {
            mem_sp5c[binAcnt] = resize_1D_int(mem_sp5c[binAcnt], mmem_sp5c, mmem_sp5c + incrClustPerPart);
        }
        mmem_sp5c = mmem_sp5c + incrClustPerPart;
    }
}