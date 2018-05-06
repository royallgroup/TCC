#include "rings.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void get_6A_clusters();

void get_basic_clusters() {	// get SP3/4/5 rings including particle n0
    int n1_pointer, n2_pointer;
    int n0, n1, n2;

    for (n0 = 0; n0 < particles_in_current_frame; n0++) {
        for (n1_pointer = 0; n1_pointer < num_bonds[n0]; n1_pointer++) {
            n1 = bond_list[n0][n1_pointer];
            if (n1 > n0) {
                for (n2_pointer = n1_pointer + 1; n2_pointer < num_bonds[n0]; n2_pointer++) {
                    n2 = bond_list[n0][n2_pointer];
                    if (n2 > n0) {
                        if (Bonds_BondCheck(n1, n2)) {
                            // SP3 found, check type and store
                            get_sp3_clusters(n0, n1, n2);
                        }
                        else { // not SP3, search for SP4 & SP5
                            if (dosp4 == 1) {
                                get_basic_sp4_rings(n0, n1, n2);
                            }
                        }
                    }
                }
            }
        }
    }
}

void get_basic_sp4_rings(int n0, int n1, int n2) {    // {n0,n1,n2} is not an SP3 ring, is it an SP4 or SP5 ring?
    int i;
    int n3;
    int tmp;

    if(n1 > n2) {
        tmp = n2;
        n2 = n1;
        n1 = tmp;
    }

    for (i=0; i<num_bonds[n1]; ++i) {
        n3 = bond_list[n1][i];
        if (n3 > n0){
            if (Bonds_BondCheck(n0, n3) == 0) {  // n1 not bonded to n2 & n0 not bonded to n3
                if (Bonds_BondCheck(n2, n3)) { // 4 membered ring found
                    get_sp4_clusters(n0, n1, n2, n3); // check SP4 type and store
                }
                else { // n1 not bonded to n3
                    if (dosp5 == 1) {
                        get_basic_sp5_rings(n0, n1, n2, n3);
                    }
                }
            }
        }       
    }
}

void get_basic_sp5_rings(int n0, int n1, int n2, int n3) {    // {n0,n1,n3,n2} is not an SP4 ring, is it an SP5 ring?
    int i, j;
    int n4, n5;
    int bond4_1;

    // Loop through all neighbours of n3
    for (i = 0; i < num_bonds[n3]; ++i){
        n4 = bond_list[n3][i];
        if(n4 > n0) {
            if (n4 != n2) {
                bond4_1 = 0;
                // loop through neighbours of n4
                for (j = 0; j < num_bonds[n4]; ++j) {
                    n5 = bond_list[n4][j];
                    // if n4 is bonded to n0 or n1 then it isn't an sp5
                    if (n5 == n1 || n5 == n0) {
                        break;
                    }
                    // if n4 is bonded to n2 then it might be an sp5
                    if (n5 == n2) {
                        bond4_1 = 1;
                    }
                }
                if (j == num_bonds[n4] && bond4_1 == 1) {
                    get_sp5_clusters(n0, n1, n3, n4, n2); // check SP5 type and store
                }
            }
        }
    }
}

void get_sp3_clusters(int n0, int n1, int n2) {    // Take {n0,n1,n2}, check SP3 ring and if so detect SP3a/b/c cluster
    int i, j;
    int type = 0;
    int cp[2];  // common spindles - particles bonded to all members of three membered ring
    int tmp;

    if(n1 > n2) {
        tmp = n2;
        n2 = n1;
        n1 = tmp;
    }

    cp[0]=cp[1]=-1;
    for (i=0; i<num_bonds[n0]; ++i) {
        j = bond_list[n0][i];
        if (j != n1 && j != n2) {
            if (Bonds_BondCheck(n1, j) == 1 && Bonds_BondCheck(n2, j) == 1) {
                if (type < 2) {
                    cp[type] = j;
                    type++;
                }
                else {
                    type++;
                }
            }
        }
    }
    
    if (type==0 && dosp3a==1) {
        Store_sp3a(n0, n1, n2);
    }
    else if (type==1 && dosp3b==1) {
        Store_sp3b(n0, n1, n2, cp);
    }
    else if (type==2 && dosp3c==1) {
        Store_sp3c(n0, n1, n2, cp);
    }
    else if (dosp3a==1) {
        Store_sp3a(n0, n1, n2);
    }
}

void Store_sp3c(int n0, int n1, int n2, const int *cp) {
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

    if (ssp3c[hcsp3c[nsp3c][0]] == 'C') ssp3c[hcsp3c[nsp3c][0]] = 'B';
    if (ssp3c[hcsp3c[nsp3c][1]] == 'C') ssp3c[hcsp3c[nsp3c][1]] = 'B';
    if (ssp3c[hcsp3c[nsp3c][2]] == 'C') ssp3c[hcsp3c[nsp3c][2]] = 'B';
    ssp3c[hcsp3c[nsp3c][3]] = 'O';
    ssp3c[hcsp3c[nsp3c][4]] = 'O';

    add_mem_sp3c(n0);
    add_mem_sp3c(n1);
    add_mem_sp3c(n2);
    add_mem_sp3c(cp[0]);
    add_mem_sp3c(cp[1]);

    ++nsp3c;
}

void Store_sp3b(int n0, int n1, int n2, const int *cp) {
    if (nsp3b == msp3b) {
        hcsp3b=resize_2D_int(hcsp3b,msp3b,msp3b+incrStatic,4,-1);
        msp3b=msp3b+incrStatic;
    }
    hcsp3b[nsp3b][0] = n0;
    hcsp3b[nsp3b][1] = n1;
    hcsp3b[nsp3b][2] = n2;
    hcsp3b[nsp3b][3] = cp[0];

    if (ssp3b[hcsp3b[nsp3b][0]] == 'C') ssp3b[hcsp3b[nsp3b][0]] = 'B';
    if (ssp3b[hcsp3b[nsp3b][1]] == 'C') ssp3b[hcsp3b[nsp3b][1]] = 'B';
    if (ssp3b[hcsp3b[nsp3b][2]] == 'C') ssp3b[hcsp3b[nsp3b][2]] = 'B';
    ssp3b[hcsp3b[nsp3b][3]] = 'O';

    add_mem_sp3b(n0);
    add_mem_sp3b(n1);
    add_mem_sp3b(n2);
    add_mem_sp3b(cp[0]);

    ++nsp3b;
}

void Store_sp3a(int n0, int n1, int n2) {
    if (nsp3a == msp3a) {
        hcsp3a=resize_2D_int(hcsp3a,msp3a,msp3a+incrStatic,3,-1);
        msp3a=msp3a+incrStatic;
    }
    hcsp3a[nsp3a][0] = n0;
    hcsp3a[nsp3a][1] = n1;
    hcsp3a[nsp3a][2] = n2;

    ssp3a[hcsp3a[nsp3a][0]] = 'B';
    ssp3a[hcsp3a[nsp3a][1]] = 'B';
    ssp3a[hcsp3a[nsp3a][2]] = 'B';

    ++nsp3a;
}

void get_sp4_clusters(int n0, int n1, int n2, int n3) {    // Take {n0,n1,n3,n2}, check SP4 ring and if so detect SP4a/b/c cluster
    int i, j;
    int type = 0;
    int cp[2];  // common spindles - particles bonded to all members of three membered ring

    cp[0]=cp[1]=-1;
    for (i=0; i<num_bonds[n0]; ++i) {
        j = bond_list[n0][i];
        if (j != n1 && j != n2) {
            if (Bonds_BondCheck(n1, j) == 1 && Bonds_BondCheck(n3, j) == 1 && Bonds_BondCheck(n2, j) == 1) {
                if (type < 2) {
                    cp[type] = j;
                    type++;
                }
                else {
                    type++;
                }
            }
        }
    }
    
    if (type==0 && dosp4a==1) {
        Store_sp4a(n0, n1, n2, n3);
    }
    else if (type==1 && dosp4b==1) {
        Store_sp4b(n0, n1, n2, n3, cp);
    }
    else if (type==2 && dosp4c==1) {
        Store_sp4c(n0, n1, n2, n3, cp);
    }
    else if (dosp4a==1) {
        Store_sp4a(n0, n1, n3, n2);
    }
}

void Store_sp4c(int n0, int n1, int n2, int n3, const int *cp) {

    if (nsp4c == msp4c) {
        hcsp4c = resize_2D_int(hcsp4c, msp4c, msp4c + incrStatic, 6, -1);
        msp4c = msp4c + incrStatic;
    }
    hcsp4c[nsp4c][0] = n0;
    hcsp4c[nsp4c][1] = n1;
    hcsp4c[nsp4c][2] = n3;
    hcsp4c[nsp4c][3] = n2;
    if (cp[0] < cp[1]) {
        hcsp4c[nsp4c][4] = cp[0];
        hcsp4c[nsp4c][5] = cp[1];
    }
    else {
        hcsp4c[nsp4c][4] = cp[1];
        hcsp4c[nsp4c][5] = cp[0];
    }

    if (ssp4c[hcsp4c[nsp4c][0]] == 'C') ssp4c[hcsp4c[nsp4c][0]] = 'B';
    if (ssp4c[hcsp4c[nsp4c][1]] == 'C') ssp4c[hcsp4c[nsp4c][1]] = 'B';
    if (ssp4c[hcsp4c[nsp4c][2]] == 'C') ssp4c[hcsp4c[nsp4c][2]] = 'B';
    if (ssp4c[hcsp4c[nsp4c][3]] == 'C') ssp4c[hcsp4c[nsp4c][3]] = 'B';
    ssp4c[hcsp4c[nsp4c][4]] = 'O';
    ssp4c[hcsp4c[nsp4c][5]] = 'O';

    add_mem_sp4c(n0);
    add_mem_sp4c(n1);
    add_mem_sp4c(n2);
    add_mem_sp4c(n3);
    add_mem_sp4c(cp[0]);
    add_mem_sp4c(cp[1]);


    get_6A_clusters();

    ++nsp4c;
}

void get_6A_clusters() {
    // The sp4c cluster can be detected multiple times within the same 6 particles by rotation.
    // The 6A cluster is a unique sp4c such that no 6A has identical particles to another 6A.
    // This means that the 6A is a subset of the sp4c

    int trial[6];
    int flg;
    int trial_pointer;

    // Check for sp4c isomers
    for (int i = 0; i < 6; ++i) {
        trial[i] = hcsp4c[nsp4c][i];
    }
    quickSort(&trial[0], 6);
    flg = 0;  // check trial cluster not already found
    for (int existing_6A_pointer = 0; existing_6A_pointer < n6A; ++existing_6A_pointer) {
        int * exisiting_6A_cluster = hc6A[existing_6A_pointer];
        for (trial_pointer = 0; trial_pointer < 6; ++trial_pointer) {
            if (trial[trial_pointer] != exisiting_6A_cluster[trial_pointer]) {
                break;
            }
        }
        if (trial_pointer==6) {
            flg=1;
            break;
        }
    }
    if (flg==0) {
        if (n6A == m6A) {
            hc6A = resize_2D_int(hc6A, m6A, m6A + incrStatic, 6, -1);
            m6A = m6A + incrStatic;
        }
        for (int i = 0; i < 6; ++i) {
            hc6A[n6A][i] = trial[i];
            s6A[trial[i]] = 'B';
        }
        ++n6A;
    }
}

void Store_sp4b(int n0, int n1, int n2, int n3, const int *cp) {
    if (nsp4b == msp4b) {
        hcsp4b=resize_2D_int(hcsp4b,msp4b,msp4b+incrStatic,5,-1);
        msp4b=msp4b+incrStatic;
    }
    hcsp4b[nsp4b][0] = n0;
    hcsp4b[nsp4b][1] = n1;
    hcsp4b[nsp4b][2] = n3;
    hcsp4b[nsp4b][3] = n2;
    hcsp4b[nsp4b][4] = cp[0];

    if (ssp4b[hcsp4b[nsp4b][0]] == 'C') ssp4b[hcsp4b[nsp4b][0]] = 'B';
    if (ssp4b[hcsp4b[nsp4b][1]] == 'C') ssp4b[hcsp4b[nsp4b][1]] = 'B';
    if (ssp4b[hcsp4b[nsp4b][2]] == 'C') ssp4b[hcsp4b[nsp4b][2]] = 'B';
    if (ssp4b[hcsp4b[nsp4b][3]] == 'C') ssp4b[hcsp4b[nsp4b][3]] = 'B';
    ssp4b[hcsp4b[nsp4b][4]] = 'O';

    add_mem_sp4b(n0);
    add_mem_sp4b(n1);
    add_mem_sp4b(n2);
    add_mem_sp4b(n3);
    add_mem_sp4b(cp[0]);

    ++nsp4b;
}

void Store_sp4a(int n0, int n1, int n2, int n3) {
    if (nsp4a == msp4a) {
        hcsp4a=resize_2D_int(hcsp4a,msp4a,msp4a+incrStatic,4,-1);
        msp4a=msp4a+incrStatic;
    }
    hcsp4a[nsp4a][0] = n0;
    hcsp4a[nsp4a][1] = n1;
    hcsp4a[nsp4a][2] = n3;
    hcsp4a[nsp4a][3] = n2;

    ssp4a[hcsp4a[nsp4a][0]] = 'B';
    ssp4a[hcsp4a[nsp4a][1]] = 'B';
    ssp4a[hcsp4a[nsp4a][2]] = 'B';
    ssp4a[hcsp4a[nsp4a][3]] = 'B';

    ++nsp4a;
}

void get_sp5_clusters(int n0, int n1, int n2, int n3, int n4) {    // Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster
    int i, j;
    int type = 0;
    int cp[2];  // common spindles - particles bonded to all members of three membered ring

    // Loop through all neighbours of n0
    for (i=0; i<num_bonds[n0]; ++i) {
        j = bond_list[n0][i];
        // If the n0 neighbour is not n1 or n4
        if (j != n1 || j != n4) {
            // If the neighbour is bonded to all other ring particles
            if (Bonds_BondCheck(n1,j)==1 && Bonds_BondCheck(n2,j)==1 && Bonds_BondCheck(n3,j)==1 && Bonds_BondCheck(n4,j)==1) {
                if (type < 2) {
                    cp[type] = j;
                    type++;
                }
                else {
                    type++;
                }
            }
        }
    }

    if (type==0 && dosp5a==1) {
        Store_sp5a(n0, n1, n2, n3, n4);

    }
    else if (type==1 && dosp5b==1) {
        Store_sp5b(n0, n1, n2, n3, n4, cp);

    }
    else if (type==2 && dosp5c==1) {
        Store_sp5c(n0, n1, n2, n3, n4, cp);

    }
    else if (dosp5a==1) {
        Store_sp5a(n0, n1, n2, n3, n4);
    }
}

void Store_sp5c(int n0, int n1, int n2, int n3, int n4, const int *cp) {
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

    if (ssp5c[hcsp5c[nsp5c][0]] == 'C') ssp5c[hcsp5c[nsp5c][0]] = 'B';
    if (ssp5c[hcsp5c[nsp5c][1]] == 'C') ssp5c[hcsp5c[nsp5c][1]] = 'B';
    if (ssp5c[hcsp5c[nsp5c][2]] == 'C') ssp5c[hcsp5c[nsp5c][2]] = 'B';
    if (ssp5c[hcsp5c[nsp5c][3]] == 'C') ssp5c[hcsp5c[nsp5c][3]] = 'B';
    if (ssp5c[hcsp5c[nsp5c][4]] == 'C') ssp5c[hcsp5c[nsp5c][4]] = 'B';
    ssp5c[hcsp5c[nsp5c][5]] = 'O';
    ssp5c[hcsp5c[nsp5c][6]] = 'O';

    ++nsp5c;
}

void Store_sp5b(int n0, int n1, int n2, int n3, int n4, const int *cp) {
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

    if (ssp5b[hcsp5b[nsp5b][0]] == 'C') ssp5b[hcsp5b[nsp5b][0]] = 'B';
    if (ssp5b[hcsp5b[nsp5b][1]] == 'C') ssp5b[hcsp5b[nsp5b][1]] = 'B';
    if (ssp5b[hcsp5b[nsp5b][2]] == 'C') ssp5b[hcsp5b[nsp5b][2]] = 'B';
    if (ssp5b[hcsp5b[nsp5b][3]] == 'C') ssp5b[hcsp5b[nsp5b][3]] = 'B';
    if (ssp5b[hcsp5b[nsp5b][4]] == 'C') ssp5b[hcsp5b[nsp5b][4]] = 'B';
    ssp5b[hcsp5b[nsp5b][5]] = 'O';

    add_mem_sp5b(n0);
    add_mem_sp5b(n1);
    add_mem_sp5b(n2);
    add_mem_sp5b(n3);
    add_mem_sp5b(n4);
    add_mem_sp5b(cp[0]);

    ++nsp5b;
}

void Store_sp5a(int n0, int n1, int n2, int n3, int n4) {
    if (nsp5a == msp5a) {
        hcsp5a=resize_2D_int(hcsp5a,msp5a,msp5a+incrStatic,5,-1);
        msp5a=msp5a+incrStatic;
    }
    hcsp5a[nsp5a][0] = n0;
    hcsp5a[nsp5a][1] = n1;
    hcsp5a[nsp5a][2] = n2;
    hcsp5a[nsp5a][3] = n3;
    hcsp5a[nsp5a][4] = n4;

    ssp5a[hcsp5a[nsp5a][0]] = 'B';
    ssp5a[hcsp5a[nsp5a][1]] = 'B';
    ssp5a[hcsp5a[nsp5a][2]] = 'B';
    ssp5a[hcsp5a[nsp5a][3]] = 'B';
    ssp5a[hcsp5a[nsp5a][4]] = 'B';

    ++nsp5a;
}

void add_mem_sp3b(int particle_ID) {
    int binAcnt;

    mem_sp3b[particle_ID][nmem_sp3b[particle_ID]]=nsp3b;
    nmem_sp3b[particle_ID]++;
    if (nmem_sp3b[particle_ID] >= mmem_sp3b) {
        for (binAcnt=0; binAcnt < particles_in_current_frame; binAcnt++) {
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
        for (binAcnt=0; binAcnt < particles_in_current_frame; binAcnt++) {
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
        for (binAcnt = 0; binAcnt < particles_in_current_frame; binAcnt++) {
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
        for (binAcnt = 0; binAcnt < particles_in_current_frame; binAcnt++) {
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
        for (binAcnt = 0; binAcnt < particles_in_current_frame; binAcnt++) {
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
        for (binAcnt = 0; binAcnt < particles_in_current_frame; binAcnt++) {
            mem_sp5c[binAcnt] = resize_1D_int(mem_sp5c[binAcnt], mmem_sp5c, mmem_sp5c + incrClustPerPart);
        }
        mmem_sp5c = mmem_sp5c + incrClustPerPart;
    }
}