#include "8A.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "string.h"

void Clusters_Get8A() { // Detect 8A D2d clusters
    int unc[2];
    int com[4];
    int i, j, j2, k, l, m;
    int cnt;
    int flg;
    int break_out;
    int trial[8];
    int clusSize=8;
    int *used_sp5b;


    used_sp5b = malloc(nsp5b * sizeof(int));
    if (used_sp5b==NULL) {
        Error("Clusters_Get8A(): used_sp5b[] malloc out of memory\n");
    }

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
